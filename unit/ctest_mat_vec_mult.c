#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_mat.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_basis.h>
#include <gkyl_array.h>
#include <gkyl_array_rio.h>
#include <gkyl_range.h>


void test_mat_mm()
{
  struct gkyl_mat *A = gkyl_mat_new(4, 4, 1);
  struct gkyl_mat *B = gkyl_mat_new(4, 1, 2);
  struct gkyl_mat *C = gkyl_mat_new(4, 1, 0);
  // C = 1.0*A*B + 0.0*C
  gkyl_mat_mm(1.0, 0.0, GKYL_NO_TRANS, A, GKYL_NO_TRANS, B, C);

  for (size_t j=0; j<C->nc; ++j)
    for (size_t i=0; i<C->nr; ++i)
      TEST_CHECK ( 8 == gkyl_mat_get(C, i, j) );
}

void test_mat_mv()
{

  struct gkyl_mat *A = gkyl_mat_new(4, 4, 1);
  struct gkyl_mat *x = gkyl_mat_new(4, 1, 2);
  struct gkyl_mat *y = gkyl_mat_new(4, 1, 1);
  // y = 1.0*A*x + 0.0*C
  gkyl_mat_mv(1.0, 0.0, GKYL_NO_TRANS, A, x, y);
  for (size_t j=0; j<y->nc; ++j)
    for (size_t i=0; i<y->nr; ++i)
      TEST_CHECK ( 8 == gkyl_mat_get(y, i, j) );
}

void test_nmat_mv()
{

  // n_do matrices with shape 4x4
  int n_do = 3;
  struct gkyl_nmat *nmat_A = gkyl_nmat_new(n_do, 4, 4);
  struct gkyl_nmat *nmat_x = gkyl_nmat_new(n_do, 4, 1);
  struct gkyl_nmat *nmat_y = gkyl_nmat_new(n_do, 4, 1);

  //fill each A matrix with its number
  for (size_t n=0; n<nmat_A->num; ++n) {
    struct gkyl_mat m = gkyl_nmat_get(nmat_A, n);
    for (size_t j=0; j<nmat_A->nc; ++j)
      for (size_t i=0; i<nmat_A->nr; ++i)
        gkyl_mat_set(&m, i, j, n);
  }

  //fill each x vecotr with ones
  for (size_t n=0; n<nmat_x->num; ++n) {
    struct gkyl_mat m = gkyl_nmat_get(nmat_x, n);
    for (size_t j=0; j<nmat_x->nc; ++j)
      for (size_t i=0; i<nmat_x->nr; ++i)
        gkyl_mat_set(&m, i, j, 1);
  }

  //fill each y matrix with 0
  for (size_t n=0; n<nmat_y->num; ++n) {
    struct gkyl_mat m = gkyl_nmat_get(nmat_y, n);
    for (size_t j=0; j<nmat_y->nc; ++j)
      for (size_t i=0; i<nmat_y->nr; ++i)
        gkyl_mat_set(&m, i, j, 0);
  }

  enum gkyl_mat_trans transa[n_do];
  double alpha[n_do];
  double beta[n_do];
  for(int i=0; i<n_do;i++){
    transa[i] = GKYL_NO_TRANS;
    alpha[i] = 1.0;
    beta[i] = 0.0;
  }

  gkyl_nmat_mv(alpha, beta, transa, nmat_A, nmat_x, nmat_y);


  //check the expected result
  for (size_t n=0; n<nmat_y->num; ++n) {
    //printf("n = %ld\n", n);
    struct gkyl_mat y = gkyl_nmat_get(nmat_y,n);
    double expected = 4*n;
    //printf("expected = %g\n", expected);
    for (size_t j=0; j<y.nc; ++j)
      for (size_t i=0; i<y.nr; ++i)
        //printf("actual =%g\n",gkyl_mat_get(&y, i, j));
        TEST_CHECK ( expected == gkyl_mat_get(&y, i, j) );
  }


}



//Will complete this test later
//void test_nmat_vec_mult_1d()
//{
//  //create the grid
//  double lower[] = {1.0}, upper[] = {2.5};
//  int cells[] = {10};
//  struct gkyl_rect_grid grid;
//  gkyl_rect_grid_init(&grid, 1, lower, upper, cells);
//  //create the ranges
//  struct gkyl_range local, local_ext;
//  int nghost[GKYL_MAX_CDIM] = {1};
//  gkyl_create_grid_ranges(&grid, nghost, &local_ext, &local);
//  //create the basis
//  struct gkyl_basis basis;
//  int poly_order = 1;
//  gkyl_cart_modal_serendip(&basis, grid.ndim, poly_order);
//
//  // n_do matrices with shape N_bxN_b
//  int n_do = 3;
//  struct gkyl_nmat *nmat = gkyl_nmat_new(n_do, basis.num_basis, basis.num_basis);
//  TEST_CHECK( false == gkyl_nmat_is_cu_dev(nmat) );
//
//  //check num, rows, cols are correct. remove later
//  TEST_CHECK( n_do == nmat->num );
//  TEST_CHECK( basis.num_basis == nmat->nr );
//  TEST_CHECK( basis.num_basis == nmat->nc );
//
//  struct gkyl_mat m = gkyl_nmat_get(nmat, 0);
//  TEST_CHECK( basis.num_basis == m.nr );
//  TEST_CHECK( basis.num_basis == m.nc );
//
//  //fill each matrix with its number
//  for (size_t n=0; n<nmat->num; ++n) {
//    struct gkyl_mat m = gkyl_nmat_get(nmat, n);
//    for (size_t j=0; j<nmat->nc; ++j)
//      for (size_t i=0; i<nmat->nr; ++i)
//        gkyl_mat_set(&m, i, j, n);
//  }
//
//  //now create a target DG field
//}


TEST_LIST = {
  //{ "nmat_vec_mult", test_nmat_vec_mult_1d},
  { "mv", test_mat_mv},
  { "mm", test_mat_mm},
  { "nmat_mv", test_nmat_mv},
  //{ "mm", test_nmat_mm},
  { NULL, NULL },
};
