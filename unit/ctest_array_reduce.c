#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_reduce.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_util.h>
#include <time.h>

void test_dummy()
{
}

void test_reduce()
{
  int ncomp = 3, ncells = 200;
  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, ncomp, ncells);

  for (size_t i=0; i<arr->size; ++i) {
    double *d = gkyl_array_fetch(arr, i);
    for (size_t c=0; c<ncomp; ++c)
      d[c] = 0.5*i + 0.1*c;
  }
  
  double amin[ncomp], amax[ncomp];
  gkyl_array_reduce(&amin[0], arr, GKYL_MIN);
  gkyl_array_reduce(&amax[0], arr, GKYL_MAX);

  for (size_t c=0; c<ncomp; ++c) {
    TEST_CHECK( amin[c] == 0.1*c );
    TEST_CHECK( amax[c] == 0.5*(ncells-1) + 0.1*c );
  }

  gkyl_array_release(arr);
}

void test_reduce_dg(bool use_gpu)
{
  int poly_order = 1;
  int ncomp = 3;
  double lower[] = {-M_PI}, upper[] = {M_PI};
  int cells[] = {20};

  int ndim = sizeof(lower)/sizeof(lower[0]);

  struct gkyl_basis *basis;
  if (use_gpu) {
    basis = gkyl_cu_malloc(sizeof(struct gkyl_basis));
    gkyl_cart_modal_serendip_cu_dev(basis, ndim, poly_order);
  }
  else {
    basis = gkyl_malloc(sizeof(struct gkyl_basis));
    gkyl_cart_modal_serendip(basis, ndim, poly_order);
  }
  struct gkyl_basis basis_ho;
  gkyl_cart_modal_serendip(&basis_ho, ndim, poly_order);

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  int ghost[ndim];
  for (int d=0; d<ndim; d++) ghost[d] = 1;
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  struct gkyl_array *arr = use_gpu? gkyl_array_cu_dev_new(GKYL_DOUBLE, ncomp*basis_ho.num_basis, local_ext.volume) 
                                  : gkyl_array_new(GKYL_DOUBLE, ncomp*basis_ho.num_basis, local_ext.volume);
  struct gkyl_array *arr_ho = use_gpu? gkyl_array_new(GKYL_DOUBLE, arr->ncomp, arr->size) : gkyl_array_acquire(arr);

  // Load 1D Gauss-Legendre nodes.
  int num_quad = poly_order+1;
  double ordinates1[num_quad];
  memcpy(ordinates1, gkyl_gauss_ordinates[num_quad], sizeof(double[num_quad]));

  // Create range to loop over nodes.
  int qshape[GKYL_MAX_DIM];
  for (int i=0; i<ndim; ++i) qshape[i] = num_quad;
  struct gkyl_range qrange;
  gkyl_range_init_from_shape(&qrange, ndim, qshape);

  // Create nodes array.
  int num_nodes = qrange.volume;
  struct gkyl_array *nodes = gkyl_array_new(GKYL_DOUBLE, ndim, num_nodes);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &qrange);
  while (gkyl_range_iter_next(&iter)) {
    long linc = gkyl_range_idx(&qrange, iter.idx);
    double *nod = gkyl_array_fetch(nodes, linc);
    for (int i=0; i<ndim; ++i)
      nod[i] = ordinates1[iter.idx[i]-qrange.lower[i]];
  }

  // Populate arr with a function evaluated at nodes (transformed to modal).
  // Keep track of nodal max, min and sum.
  double arr_max[ncomp], arr_min[ncomp], arr_sum[ncomp];
  for (size_t i=0; i<ncomp; ++i) {
    arr_max[i] = -1e20;
    arr_min[i] = 1e20;
    arr_sum[i] = 0.0;
  }

  for (size_t i=0; i<arr->size; ++i) {

    double *arr_c = gkyl_array_fetch(arr_ho, i);

    int idx[ndim];
    gkyl_range_inv_idx(&local_ext, i, idx);

    double xc[ndim];
    gkyl_rect_grid_cell_center(&grid, idx, xc);

    for (int ci=0; ci<ncomp; ci++) {
      double arr_nodal[num_nodes];
      for (size_t k=0; k<num_nodes; k++) {
        const double *nod = gkyl_array_cfetch(nodes, k);
        double x[ndim];
        for (int d=0; d<ndim; d++) x[d] = xc[d] + 0.5*grid.dx[0]*nod[d];
  
        arr_nodal[k] = (ci+1);
        for (int d=0; d<ndim; d++) arr_nodal[k] *= sin(((d+1)*2.0*M_PI/(upper[d]-lower[d]))*x[d]);
  
        arr_max[ci] = GKYL_MAX2(arr_max[ci], arr_nodal[k]);
        arr_min[ci] = GKYL_MIN2(arr_min[ci], arr_nodal[k]);
        arr_sum[ci] += arr_nodal[k];
      }
  
      for (size_t k=0; k<basis_ho.num_basis; k++) {
        basis_ho.quad_nodal_to_modal(arr_nodal, &arr_c[ci*basis_ho.num_basis], k);
      }
    }
  }
  gkyl_array_copy(arr, arr_ho);

  // Reduce each vector component.
  double *amax, *amin, *asum;
  if (use_gpu) {
    amax = gkyl_cu_malloc(ncomp*sizeof(double));
    amin = gkyl_cu_malloc(ncomp*sizeof(double));
    asum = gkyl_cu_malloc(ncomp*sizeof(double));
  }
  else {
    amax = gkyl_malloc(ncomp*sizeof(double));
    amin = gkyl_malloc(ncomp*sizeof(double));
    asum = gkyl_malloc(ncomp*sizeof(double));
  }

  for (int ci=0; ci<ncomp; ci++) {
    gkyl_array_reducec_dg(&amax[ci], arr, ci, GKYL_MAX, basis);
    gkyl_array_reducec_dg(&amin[ci], arr, ci, GKYL_MIN, basis);
    gkyl_array_reducec_dg(&asum[ci], arr, ci, GKYL_SUM, basis);
  }

  // Check results.
  double amax_ho[ncomp], amin_ho[ncomp], asum_ho[ncomp];
  if (use_gpu) {
    gkyl_cu_memcpy(amax_ho, amax, ncomp*sizeof(double), GKYL_CU_MEMCPY_D2H);
    gkyl_cu_memcpy(amin_ho, amin, ncomp*sizeof(double), GKYL_CU_MEMCPY_D2H);
    gkyl_cu_memcpy(asum_ho, asum, ncomp*sizeof(double), GKYL_CU_MEMCPY_D2H);
  }
  else {
    memcpy(amax_ho, amax, ncomp*sizeof(double));
    memcpy(amin_ho, amin, ncomp*sizeof(double));
    memcpy(asum_ho, asum, ncomp*sizeof(double));
  }

  for (int c=0; c<ncomp; ++c) {
    TEST_CHECK( gkyl_compare(arr_max[c], amax_ho[c], 1e-14) );
    TEST_MSG( "%d MAX: Expected: %g | Got: %g\n", c, arr_max[c], amax_ho[c] );
    TEST_CHECK( gkyl_compare(arr_min[c], amin_ho[c], 1e-14) );
    TEST_MSG( "%d MIN: Expected: %g | Got: %g\n", c, arr_min[c], amin_ho[c] );
    TEST_CHECK( gkyl_compare(arr_sum[c], asum_ho[c], 1e-14) );
    TEST_MSG( "%d SUM: Expected: %g | Got: %g\n", c, arr_sum[c], asum_ho[c] );
  }

  gkyl_array_release(arr_ho);
  gkyl_array_release(arr);
  gkyl_array_release(nodes);
  if (use_gpu) {
    gkyl_cart_modal_basis_release_cu(basis);
    gkyl_cu_free(amax);
    gkyl_cu_free(amin);
    gkyl_cu_free(asum);
  }
  else {
    gkyl_cart_modal_basis_release(basis);
    gkyl_free(amax);
    gkyl_free(amin);
    gkyl_free(asum);
  }
}

void test_reduce_range()
{
  int shape[] = {10, 20};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, 3, range.volume);

  int count = -1000;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);
    double *d = gkyl_array_fetch(arr, loc);
    d[0] = count+0.5;
    d[1] = count+1.5;
    d[2] = count+2.5;

    count += 1;
  }

  double amin[3], amax[3];
  gkyl_array_reduce_range(amin, arr, GKYL_MIN, &range);

  TEST_CHECK( amin[0] == -999.5 );
  TEST_CHECK( amin[1] == -998.5 );
  TEST_CHECK( amin[2] == -997.5 );

  gkyl_array_reduce_range(amax, arr, GKYL_MAX, &range);
  
  TEST_CHECK( amax[0] == -800.5 );
  TEST_CHECK( amax[1] == -799.5 );
  TEST_CHECK( amax[2] == -798.5 );

  gkyl_array_release(arr);
}

void test_sum_reduce_range()
{
  int shape[] = {10, 20};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, 3, range.volume);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);
    double *d = gkyl_array_fetch(arr, loc);
    d[0] = 0.5;
    d[1] = 1.5;
    d[2] = 2.5;
  }

  double asum[3];
  gkyl_array_reduce_range(asum, arr, GKYL_SUM, &range);

  TEST_CHECK( asum[0] == 0.5*range.volume );
  TEST_CHECK( asum[1] == 1.5*range.volume );
  TEST_CHECK( asum[2] == 2.5*range.volume );

  gkyl_array_release(arr);
}

void test_reduce_dg_range(bool use_gpu)
{
  int poly_order = 1;
  int ncomp = 3;
  double lower[] = {-M_PI}, upper[] = {M_PI};
  int cells[] = {20};

  int ndim = sizeof(lower)/sizeof(lower[0]);

  struct gkyl_basis *basis;
  if (use_gpu) {
    basis = gkyl_cu_malloc(sizeof(struct gkyl_basis));
    gkyl_cart_modal_serendip_cu_dev(basis, ndim, poly_order);
  }
  else {
    basis = gkyl_malloc(sizeof(struct gkyl_basis));
    gkyl_cart_modal_serendip(basis, ndim, poly_order);
  }
  struct gkyl_basis basis_ho;
  gkyl_cart_modal_serendip(&basis_ho, ndim, poly_order);

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  int ghost[ndim];
  for (int d=0; d<ndim; d++) ghost[d] = 1;
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  struct gkyl_array *arr = use_gpu? gkyl_array_cu_dev_new(GKYL_DOUBLE, ncomp*basis_ho.num_basis, local_ext.volume) 
                                  : gkyl_array_new(GKYL_DOUBLE, ncomp*basis_ho.num_basis, local_ext.volume);
  struct gkyl_array *arr_ho = use_gpu? gkyl_array_new(GKYL_DOUBLE, arr->ncomp, arr->size) : gkyl_array_acquire(arr);


  // Load 1D Gauss-Legendre nodes.
  int num_quad = poly_order+1;
  double ordinates1[num_quad];
  memcpy(ordinates1, gkyl_gauss_ordinates[num_quad], sizeof(double[num_quad]));

  // Create range to loop over nodes.
  int qshape[GKYL_MAX_DIM];
  for (int i=0; i<ndim; ++i) qshape[i] = num_quad;
  struct gkyl_range qrange;
  gkyl_range_init_from_shape(&qrange, ndim, qshape);

  // Create nodes array.
  int num_nodes = qrange.volume;
  struct gkyl_array *nodes = gkyl_array_new(GKYL_DOUBLE, ndim, num_nodes);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &qrange);
  while (gkyl_range_iter_next(&iter)) {
    long linc = gkyl_range_idx(&qrange, iter.idx);
    double *nod = gkyl_array_fetch(nodes, linc);
    for (int i=0; i<ndim; ++i)
      nod[i] = ordinates1[iter.idx[i]-qrange.lower[i]];
  }

  // Populate arr with a function evaluated at nodes (transformed to modal).
  double arr_max[ncomp], arr_min[ncomp], arr_sum[ncomp];
  for (size_t i=0; i<ncomp; ++i) {
    arr_max[i] = -1e20;
    arr_min[i] = 1e20;
    arr_sum[i] = 0.0;
  }

  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&local, iter.idx);

    double *arr_c = gkyl_array_fetch(arr_ho, linidx);

    double xc[ndim];
    gkyl_rect_grid_cell_center(&grid, iter.idx, xc);

    for (int ci=0; ci<ncomp; ci++) {
      double arr_nodal[num_nodes];
      for (size_t k=0; k<num_nodes; k++) {
        const double *nod = gkyl_array_cfetch(nodes, k);
        double x[ndim];
        for (int d=0; d<ndim; d++) x[d] = xc[d] + 0.5*grid.dx[0]*nod[d];
  
        arr_nodal[k] = (ci+1);
        for (int d=0; d<ndim; d++) arr_nodal[k] *= sin(((d+1)*2.0*M_PI/(upper[d]-lower[d]))*x[d]);
  
        arr_max[ci] = GKYL_MAX2(arr_max[ci], arr_nodal[k]);
        arr_min[ci] = GKYL_MIN2(arr_min[ci], arr_nodal[k]);
        arr_sum[ci] += arr_nodal[k];
      }
  
      for (size_t k=0; k<basis_ho.num_basis; k++) {
        basis_ho.quad_nodal_to_modal(arr_nodal, &arr_c[ci*basis_ho.num_basis], k);
      }
    }
  }
  gkyl_array_copy(arr, arr_ho);
  
  // Reduce each vector component.
  double *amax, *amin, *asum;
  if (use_gpu) {
    amax = gkyl_cu_malloc(ncomp*sizeof(double));
    amin = gkyl_cu_malloc(ncomp*sizeof(double));
    asum = gkyl_cu_malloc(ncomp*sizeof(double));
  }
  else {
    amax = gkyl_malloc(ncomp*sizeof(double));
    amin = gkyl_malloc(ncomp*sizeof(double));
    asum = gkyl_malloc(ncomp*sizeof(double));
  }
  for (int ci=0; ci<ncomp; ci++) {
    gkyl_array_reducec_dg_range(&amax[ci], arr, ci, GKYL_MAX, basis, &local);
    gkyl_array_reducec_dg_range(&amin[ci], arr, ci, GKYL_MIN, basis, &local);
    gkyl_array_reducec_dg_range(&asum[ci], arr, ci, GKYL_SUM, basis, &local);
  }

  // Check results.
  double amax_ho[ncomp], amin_ho[ncomp], asum_ho[ncomp];
  if (use_gpu) {
    gkyl_cu_memcpy(amax_ho, amax, ncomp*sizeof(double), GKYL_CU_MEMCPY_D2H);
    gkyl_cu_memcpy(amin_ho, amin, ncomp*sizeof(double), GKYL_CU_MEMCPY_D2H);
    gkyl_cu_memcpy(asum_ho, asum, ncomp*sizeof(double), GKYL_CU_MEMCPY_D2H);
  }
  else {
    memcpy(amax_ho, amax, ncomp*sizeof(double));
    memcpy(amin_ho, amin, ncomp*sizeof(double));
    memcpy(asum_ho, asum, ncomp*sizeof(double));
  }

  for (int c=0; c<ncomp; ++c) {
    TEST_CHECK( gkyl_compare(arr_max[c], amax_ho[c], 1e-14) );
    TEST_MSG( "%d MAX: Expected: %g | Got: %g\n", c, arr_max[c], amax_ho[c] );
    TEST_CHECK( gkyl_compare(arr_min[c], amin_ho[c], 1e-14) );
    TEST_MSG( "%d MIN: Expected: %g | Got: %g\n", c, arr_min[c], amin_ho[c] );
    TEST_CHECK( gkyl_compare(arr_sum[c], asum_ho[c], 1e-14) );
    TEST_MSG( "%d SUM: Expected: %g | Got: %g\n", c, arr_sum[c], asum_ho[c] );
  }

  gkyl_array_release(arr_ho);
  gkyl_array_release(arr);
  gkyl_array_release(nodes);
  if (use_gpu) {
    gkyl_cart_modal_basis_release_cu(basis);
    gkyl_cu_free(amax);
    gkyl_cu_free(amin);
    gkyl_cu_free(asum);
  }
  else {
    gkyl_cart_modal_basis_release(basis);
    gkyl_free(amax);
    gkyl_free(amin);
    gkyl_free(asum);
  }
}

void test_reduce_dg_ho() {
  test_reduce_dg(false);
}

void test_reduce_dg_range_ho() {
  test_reduce_dg_range(false);
}

// CUDA specific tests
#ifdef GKYL_HAVE_CUDA

void test_cu_array_reduce_max()
{

  unsigned long numComp = 3, numCells = 10;
  // Create host array and device copy
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, numComp, numCells);
  struct gkyl_array *a1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, numComp, numCells);
  // Initialize data
  double *a1_d = a1->data;
  for (unsigned i=0; i<numCells; ++i) {
    for (unsigned k=0; k<numComp; ++k) {
      a1_d[i*numComp+k] = (double)i+(double)k*0.10;
  }}
  gkyl_array_copy(a1_cu, a1);
  // create a device and host arrays to store reduction
  double* a1max = gkyl_malloc(numComp*sizeof(double));
  double* a1max_cu = (double*) gkyl_cu_malloc(numComp*sizeof(double));

  // Component-wise reduce array.
  gkyl_array_reduce(a1max_cu, a1_cu, GKYL_MAX);

  // Copy to host and check values.
  gkyl_cu_memcpy(a1max, a1max_cu, numComp*sizeof(double), GKYL_CU_MEMCPY_D2H);
  for (unsigned k=0; k<numComp; ++k) {
    TEST_CHECK( gkyl_compare(a1max[k], (double)(numCells-1)+(double)k*0.10, 1e-14) );
  }

  gkyl_free(a1max);
  gkyl_cu_free(a1max_cu);
  gkyl_array_release(a1);
  gkyl_array_release(a1_cu);

}

void test_cu_array_reduce_max_big()
{

  unsigned long numComp = 12, numCells = 1000;
  // Create host array and device copy
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, numComp, numCells);
  struct gkyl_array *a1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, numComp, numCells);
  // Initialize data
  double *a1_d = a1->data;
  for (unsigned i=0; i<numCells; ++i) {
    for (unsigned k=0; k<numComp; ++k) {
      a1_d[i*numComp+k] = (double)i+(double)k*0.10;
  }}
  gkyl_array_copy(a1_cu, a1);
  // create a device and host arrays to store reduction
  double* a1max = gkyl_malloc(numComp*sizeof(double));
  double* a1max_cu = (double*) gkyl_cu_malloc(numComp*sizeof(double));

  // Component-wise reduce array.
  gkyl_array_reduce(a1max_cu, a1_cu, GKYL_MAX);

  // Copy to host and check values.
  gkyl_cu_memcpy(a1max, a1max_cu, numComp*sizeof(double), GKYL_CU_MEMCPY_D2H);
  for (unsigned k=0; k<numComp; ++k) {
    // Make print statements to manually check the comparison
    printf("a1max[%d] = %g\n", k, a1max[k]);
    printf("numCells-1+(double)k*0.10 = %g\n", (double)(numCells-1)+(double)k*0.10);
    TEST_CHECK( gkyl_compare(a1max[k], (double)(numCells-1)+(double)k*0.10, 1e-14) );
  }

  gkyl_free(a1max);
  gkyl_cu_free(a1max_cu);
  gkyl_array_release(a1);
  gkyl_array_release(a1_cu);

}

void test_cu_array_reduce_range_1d_max()
{
  unsigned long numComp = 1, numCells = 10;
  // create host array and device copy
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, numComp, numCells);
  struct gkyl_array *a1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, numComp, numCells);
  // initialize data
  double *a1_d = a1->data;
  for (unsigned i=0; i<numCells; ++i) {
    for (unsigned k=0; k<numComp; ++k) {
      a1_d[i*numComp+k] = (double)i+(double)k*0.10;
  }}
  gkyl_array_copy(a1_cu, a1);
  // create a device and host arrays to store reduction
  double* a1max = gkyl_malloc(numComp*sizeof(double));
  double* a1max_cu = (double*) gkyl_cu_malloc(numComp*sizeof(double));

  int lower[] = {1}, upper[] = {numCells};
  int sublower[] = {3}, subupper[] = {7};
  struct gkyl_range range, subrange;
  gkyl_range_init(&range, 1, lower, upper);
  gkyl_sub_range_init(&subrange, &range, sublower, subupper);

  // Component-wise reduce array on range
  gkyl_array_reduce_range(a1max_cu, a1_cu, GKYL_MAX, &range);
  // Copy to host and check values.
  gkyl_cu_memcpy(a1max, a1max_cu, numComp*sizeof(double), GKYL_CU_MEMCPY_D2H);
  for (unsigned k=0; k<numComp; ++k) {
    TEST_CHECK( gkyl_compare(a1max[k], (double)(upper[0]-1)+(double)k*0.10, 1e-14) );
  }

  // Component-wise reduce array on subrange
  gkyl_cu_memset(a1max_cu, 0, sizeof(double)*numComp);
  gkyl_array_reduce_range(a1max_cu, a1_cu, GKYL_MAX, &subrange);

  // Copy to host and check values.
  gkyl_cu_memcpy(a1max, a1max_cu, numComp*sizeof(double), GKYL_CU_MEMCPY_D2H);
  for (unsigned k=0; k<numComp; ++k) {
    TEST_CHECK( gkyl_compare(a1max[k], (double)(subupper[0]-1)+(double)k*0.10, 1e-14) );
  }

  gkyl_free(a1max);
  gkyl_cu_free(a1max_cu);
  gkyl_array_release(a1);
  gkyl_array_release(a1_cu);
}

void test_cu_array_reduce_range_2d_max()
{
  int cells[] = {8, 10};
  int ghost[] = {1, 0};
  double lower[] = {0., -1.};
  double upper[] = {1., 1.};

  struct gkyl_rect_grid grid;
  struct gkyl_range range, range_ext;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);
  gkyl_create_grid_ranges(&grid, ghost, &range_ext, &range);

  unsigned long numComp = 1;
  unsigned long numCells = range_ext.volume;
  // create host array and device copy
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, numComp, numCells);
  struct gkyl_array *a1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, numComp, numCells);
  struct gkyl_array *b1 = gkyl_array_new(GKYL_DOUBLE, 1, numCells);
  struct gkyl_array *b1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, numCells);
  // initialize data
  double *a1_d = a1->data;
  double *b1_d = b1->data;
  for (unsigned i=0; i<numCells; ++i) {
    b1_d[i] = (double) i+1000;
    for (unsigned k=0; k<numComp; ++k) {
      a1_d[i*numComp+k] = (double)numCells-1-i+(double)k*0.10;
  }}
  gkyl_array_copy(a1_cu, a1);
  gkyl_array_copy(b1_cu, b1);
  // create a device and host arrays to store reduction
  double* a1max = gkyl_malloc(numComp*sizeof(double));
  double* a1max_correct = gkyl_malloc(numComp*sizeof(double));
  double* a1max_cu = (double*) gkyl_cu_malloc(numComp*sizeof(double));
  double* b1max = gkyl_malloc(1*sizeof(double));
  double* b1max_cu = (double*) gkyl_cu_malloc(1*sizeof(double));
  gkyl_cu_memset(a1max_cu, 0, sizeof(double)*numComp);
  gkyl_cu_memset(b1max_cu, 0, sizeof(double)*1);

  // Component-wise reduce array on range
  gkyl_array_reduce_range(b1max_cu, b1_cu, GKYL_MAX, &range_ext);
  // Copy to host and check values.
  gkyl_cu_memcpy(b1max, b1max_cu, sizeof(double), GKYL_CU_MEMCPY_D2H);
  TEST_CHECK( gkyl_compare(b1max[0], (double)(1000+numCells-1), 1e-14) );

  // Component-wise reduce array on subrange
  gkyl_cu_memset(a1max_cu, 0, sizeof(double)*numComp);
  gkyl_array_reduce_range(a1max_cu, a1_cu, GKYL_MAX, &range);

  gkyl_array_reduce_range(a1max_correct, a1, GKYL_MAX, &range);

  // Copy to host and check values.
  gkyl_cu_memcpy(a1max, a1max_cu, numComp*sizeof(double), GKYL_CU_MEMCPY_D2H);
  for (unsigned k=0; k<numComp; ++k) {
    TEST_CHECK( gkyl_compare(a1max[k], a1max_correct[k], 1e-14) );
  }

  gkyl_free(a1max);
  gkyl_cu_free(a1max_cu);
  gkyl_array_release(a1);
  gkyl_array_release(a1_cu);

  gkyl_cu_free(b1max);
  gkyl_cu_free(b1max_cu);   
  gkyl_array_release(b1);
  gkyl_array_release(b1_cu);
}

void
test_cu_array_reduce_range_max_timer(int NX, int NY, int VX, int VY)
{
  int cells[] = {NX, NY, VX, VY};
  int ghost[] = {1, 1, 0, 0};
  double lower[] = {0.0, 0.0, 0.0, 0.0};
  double upper[] = {1.0, 1.0, 1.0, 1.0};

  int ndim = 4;
  struct gkyl_rect_grid grid;
  struct gkyl_range range, range_ext;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  gkyl_create_grid_ranges(&grid, ghost, &range_ext, &range);

  unsigned long numComp = 1;
  unsigned long numCells = range_ext.volume;
  // create host array and device copy
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, numComp, numCells);
  struct gkyl_array *a1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, numComp, numCells);
  struct gkyl_array *b1 = gkyl_array_new(GKYL_DOUBLE, 1, numCells);
  struct gkyl_array *b1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, numCells);
  // initialize data
  double *a1_d = a1->data;
  double *b1_d = b1->data;
  for (unsigned i=0; i<numCells; ++i) {
    b1_d[i] = (double) i+1000;
    for (unsigned k=0; k<numComp; ++k) {
      a1_d[i*numComp+k] = (double)numCells-1-i+(double)k*0.10;
  }}
  gkyl_array_copy(a1_cu, a1);
  gkyl_array_copy(b1_cu, b1);
  // create a device and host arrays to store reduction
  double* a1max = gkyl_malloc(numComp*sizeof(double));
  double* a1max_correct = gkyl_malloc(numComp*sizeof(double));
  double* a1max_cu = (double*) gkyl_cu_malloc(numComp*sizeof(double));
    
  double* b1max = gkyl_malloc(1*sizeof(double));
  double* b1max_cu = (double*) gkyl_cu_malloc(1*sizeof(double));
  
  gkyl_cu_memset(a1max_cu, 0, sizeof(double)*numComp);
  gkyl_cu_memset(b1max_cu, 0, sizeof(double)*1);

  // Component-wise reduce array on subrange
  gkyl_cu_memset(a1max_cu, 0, sizeof(double)*numComp);

  struct timespec tm = gkyl_wall_clock();

  for (int i=0; i<100; ++i)
    gkyl_array_reduce_range(a1max_cu, a1_cu, GKYL_MAX, &range);
  double red_tm = gkyl_time_diff_now_sec(tm);

  printf("100 reductions on (%d,%d,%d,%d) took %g sec\n", NX, NY, VX, VY,
    red_tm);

  gkyl_free(a1max);
  gkyl_cu_free(a1max_cu);
  gkyl_array_release(a1);
  gkyl_array_release(a1_cu);

  gkyl_cu_free(b1max);
  gkyl_cu_free(b1max_cu);   
  gkyl_array_release(b1);
  gkyl_array_release(b1_cu);
}

void
test_cu_array_reduce_range_max_timer_32x32x40x40()
{
  test_cu_array_reduce_range_max_timer(32, 32, 40, 40);
}

void
test_cu_array_reduce_range_max_timer_32x32x32x32()
{
  test_cu_array_reduce_range_max_timer(32, 32, 32, 32);
}

void test_reduce_dg_dev() {
  test_reduce_dg(true);
}

void test_reduce_dg_range_dev() {
  test_reduce_dg_range(true);
}

#endif

TEST_LIST = {
  { "dummy", test_dummy },
  { "array_reduce", test_reduce },
  { "array_reduce_dg", test_reduce_dg_ho },
  { "array_reduce_range", test_reduce_range },
  { "array_reduce_sum_range", test_sum_reduce_range },
  { "array_reduce_dg_range", test_reduce_dg_range_ho },
#ifdef GKYL_HAVE_CUDA
  { "cu_array_reduce_max", test_cu_array_reduce_max },
  { "cu_array_reduce_max_big", test_cu_array_reduce_max_big },
  { "cu_array_reduce_range_1d_max", test_cu_array_reduce_range_1d_max },
  { "cu_array_reduce_range_2d_max", test_cu_array_reduce_range_2d_max },
  { "cu_array_reduce_range_max_timer_32x32x40x40", test_cu_array_reduce_range_max_timer_32x32x40x40  },
  { "cu_array_reduce_range_max_timer_32x32x32x32", test_cu_array_reduce_range_max_timer_32x32x32x32  },  
  { "cu_array_reduce_dg", test_reduce_dg_dev },
  { "cu_array_reduce_dg_range", test_reduce_dg_range_dev },
#endif
  { NULL, NULL },
};
