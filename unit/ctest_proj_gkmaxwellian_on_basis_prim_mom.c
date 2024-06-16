#include "gkyl_array.h"
#include "gkyl_util.h"
#include <acutest.h>

#include <gkyl_array_rio.h>
#include <gkyl_proj_maxwellian_on_basis.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_velocity_map.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_array_ops.h>

// allocate array (filled with zeros)
static struct gkyl_array*
mkarr(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

struct skin_ghost_ranges {
  struct gkyl_range lower_skin[GKYL_MAX_DIM];
  struct gkyl_range lower_ghost[GKYL_MAX_DIM];

  struct gkyl_range upper_skin[GKYL_MAX_DIM];
  struct gkyl_range upper_ghost[GKYL_MAX_DIM];
};

// Create ghost and skin sub-ranges given a parent range
static void
skin_ghost_ranges_init(struct skin_ghost_ranges *sgr,
  const struct gkyl_range *parent, const int *ghost)
{
  int ndim = parent->ndim;

  for (int d=0; d<ndim; ++d) {
    gkyl_skin_ghost_ranges(&sgr->lower_skin[d], &sgr->lower_ghost[d],
      d, GKYL_LOWER_EDGE, parent, ghost);
    gkyl_skin_ghost_ranges(&sgr->upper_skin[d], &sgr->upper_ghost[d],
      d, GKYL_UPPER_EDGE, parent, ghost);
  }
}

void eval_den(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 1.0;
}

void eval_udrift_2v_gk(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 0.5;
}

void eval_vtsq(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double vtsq = 1.0;
  fout[0] = vtsq;
}

void eval_bmag(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 1.0;
}

void eval_jacob_tot(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 1.0;
}

void
test_1x2v_gk(int poly_order, bool use_gpu)
{
  double mass = 1.0;
  double lower[] = {0.1, -6.0, 0.0}, upper[] = {1.0, 6.0, 6.0};
  int cells[] = {2, 16, 16};
  const int vdim = 2;
  const int ndim = sizeof(cells)/sizeof(cells[0]);
  const int cdim = ndim-vdim;

  double confLower[cdim], confUpper[cdim];
  int confCells[cdim];
  for (int d=0; d<cdim; d++) {
    confLower[d] = lower[d];
    confUpper[d] = upper[d];
    confCells[d] = cells[d];
  }
  double vLower[vdim], vUpper[vdim];
  int vCells[vdim];
  for (int d=0; d<vdim; d++) {
    vLower[d] = lower[cdim+d];
    vUpper[d] = upper[cdim+d];
    vCells[d] = cells[cdim+d];
  }

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);
  struct gkyl_rect_grid vGrid;
  gkyl_rect_grid_init(&vGrid, vdim, vLower, vUpper, vCells);

  // basis functions
  struct gkyl_basis basis, confBasis;
  if (poly_order == 1) 
    gkyl_cart_modal_gkhybrid(&basis, cdim, vdim);
  else
    gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  struct gkyl_basis *basis_on_dev, *conf_basis_on_dev;
  if (use_gpu) {
#ifdef GKYL_HAVE_CUDA
    basis_on_dev = gkyl_cu_malloc(sizeof(struct gkyl_basis));
    conf_basis_on_dev = gkyl_cu_malloc(sizeof(struct gkyl_basis));
    if (poly_order == 1) 
      gkyl_cart_modal_gkhybrid_cu_dev(basis_on_dev, cdim, vdim);
    else
      gkyl_cart_modal_serendip_cu_dev(basis_on_dev, ndim, poly_order);
    gkyl_cart_modal_serendip_cu_dev(conf_basis_on_dev, cdim, poly_order);
#endif
  }
  else 
    basis_on_dev = &basis;
    conf_basis_on_dev = &confBasis;

  int confGhost[] = { 1 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);
  struct skin_ghost_ranges confSkin_ghost; // conf-space skin/ghost
  skin_ghost_ranges_init(&confSkin_ghost, &confLocal_ext, confGhost);

  int ghost[] = { confGhost[0], 0, 0 };
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  struct skin_ghost_ranges skin_ghost; // phase-space skin/ghost
  skin_ghost_ranges_init(&skin_ghost, &local_ext, ghost);

  int vGhost[] = {0, 0};
  struct gkyl_range vLocal, vLocal_ext;
  gkyl_create_grid_ranges(&vGrid, vGhost, &vLocal_ext, &vLocal);

  // create primitive moment arrays
  struct gkyl_array *den, *udrift, *vtsq;
  den = mkarr(confBasis.num_basis, confLocal_ext.volume);
  udrift = mkarr(confBasis.num_basis, confLocal_ext.volume);
  vtsq = mkarr(confBasis.num_basis, confLocal_ext.volume);

  gkyl_proj_on_basis *proj_den = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_den, NULL);
  gkyl_proj_on_basis *proj_udrift = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, vdim, eval_udrift_2v_gk, NULL);
  gkyl_proj_on_basis *proj_vtsq = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_vtsq, NULL);

  gkyl_proj_on_basis_advance(proj_den, 0.0, &confLocal, den);
  gkyl_proj_on_basis_advance(proj_udrift, 0.0, &confLocal, udrift);
  gkyl_proj_on_basis_advance(proj_vtsq, 0.0, &confLocal, vtsq);

  // proj_maxwellian expects the primitive moments as a single array.
  struct gkyl_array *prim_moms_ho = mkarr(3*confBasis.num_basis, confLocal_ext.volume);
  gkyl_array_set_offset(prim_moms_ho, 1., den, 0*confBasis.num_basis);
  gkyl_array_set_offset(prim_moms_ho, 1., udrift, 1*confBasis.num_basis);
  gkyl_array_set_offset(prim_moms_ho, 1., vtsq  , 2*confBasis.num_basis);
  struct gkyl_array *prim_moms;
  if (use_gpu) { // copy host array to device
    prim_moms = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*confBasis.num_basis, confLocal_ext.volume);
    gkyl_array_copy(prim_moms, prim_moms_ho);
  } else {
    prim_moms = prim_moms_ho;
  }

  // create bmag and jacob_tot arrays
  struct gkyl_array *bmag, *jacob_tot;
  bmag = mkarr(confBasis.num_basis, confLocal_ext.volume);
  jacob_tot = mkarr(confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *bmag_cu, *jacob_tot_cu;
  if (use_gpu) { // create device copies
    bmag_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
    jacob_tot_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  }
  gkyl_proj_on_basis *proj_bmag = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_bmag, NULL);
  gkyl_proj_on_basis *proj_jac = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_jacob_tot, NULL);
  gkyl_proj_on_basis_advance(proj_bmag, 0.0, &confLocal, bmag);
  gkyl_proj_on_basis_advance(proj_jac, 0.0, &confLocal, jacob_tot);
  if (use_gpu) {  // copy host array to device
    gkyl_array_copy(bmag_cu, bmag);
    gkyl_array_copy(jacob_tot_cu, jacob_tot);
  }

  // create distribution function array
  struct gkyl_array *distf;
  distf = mkarr(basis.num_basis, local_ext.volume);
  struct gkyl_array *distf_cu;
  if (use_gpu)  // create device copy.
    distf_cu  = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);

  // Velocity space mapping.
  struct gkyl_mapc2p_inp c2p_in = { };
  struct gkyl_velocity_map *gvm = gkyl_velocity_map_new(c2p_in, grid, vGrid,
    local, local_ext, vLocal, vLocal_ext, use_gpu);

  // projection updater to compute Maxwellian
  struct gkyl_proj_maxwellian_on_basis_inp inp_proj = {
    .grid = &grid,
    .conf_basis = &confBasis,
    .phase_basis = &basis,
    .phase_basis_on_dev = basis_on_dev, 
    .conf_basis_on_dev = conf_basis_on_dev, 
    .phase_range = &local, 
    .phase_range_ext = &local_ext, 
    .conf_range = &confLocal, 
    .conf_range_ext = &confLocal_ext, 
    .vel_map = gvm,
    .use_gpu = use_gpu,
  };
  gkyl_proj_maxwellian_on_basis *proj_max = gkyl_proj_maxwellian_on_basis_inew(&inp_proj);

  if (use_gpu) {
    printf("flag 0\n");
    gkyl_proj_gkmaxwellian_on_basis_prim_mom(proj_max, &local, &confLocal, prim_moms,
                                             bmag_cu, jacob_tot_cu, mass, distf_cu);
    printf("flag 1\n");
    gkyl_array_copy(distf, distf_cu);
    printf("flag 2\n");
  } else {
    gkyl_proj_gkmaxwellian_on_basis_prim_mom(proj_max, &local, &confLocal, prim_moms,
                                             bmag, jacob_tot, mass, distf);
  }

  // values to compare  at index (1, 9, 9) [remember, lower-left index is (1,1,1)]
  double p1_vals[] = {  
     7.2307139183122714e-03, 0.0000000000000000e+00, 1.9198293226362615e-04, -7.7970439910196674e-04, 0.0000000000000000e+00, 0.0000000000000000e+00,
    -2.0701958137127286e-05, 0.0000000000000000e+00, -1.4953406100022537e-04, 0.0000000000000000e+00, 1.6124599381836546e-05, 0.0000000000000000e+00,
    -8.2719200283232917e-19, 0.0000000000000000e+00, -3.4806248503322844e-20, 0.0000000000000000e+00, };
  double p2_vals[] = { 
    7.2307468609012666e-03, 0.0000000000000000e+00, 1.9198380692343289e-04, -7.8092230706225602e-04, 0.0000000000000000e+00, 0.0000000000000000e+00,
    -2.0734294852987710e-05, 3.6591823321385775e-18, -1.4953474226616330e-04, 3.7739922227981074e-05, 0.0000000000000000e+00, 7.0473141211557788e-19,
    0.0000000000000000e+00, -4.8789097761847700e-19, 1.6149786206441256e-05, 0.0000000000000000e+00, 1.0020339643610290e-06, 5.4210108624275222e-20,
    0.0000000000000000e+00, 0.0000000000000000e+00 };

  const double *fv = gkyl_array_cfetch(distf, gkyl_range_idx(&local_ext, (int[3]) { 1, 9, 9 }));
  if (poly_order == 1) {
    for (int i=0; i<basis.num_basis; ++i)
      TEST_CHECK( gkyl_compare_double(p1_vals[i], fv[i], 1e-12) );
  }

  if (poly_order == 2) {
    for (int i=0; i<basis.num_basis; ++i)
      TEST_CHECK( gkyl_compare_double(p2_vals[i], fv[i], 1e-12) );
  }

  // write distribution function to file
  char fname[1024];
  sprintf(fname, "ctest_proj_maxwellian_on_basis_test_1x2v_p%d.gkyl", poly_order);
  gkyl_grid_sub_array_write(&grid, &local, 0, distf, fname);

  gkyl_array_release(den); 
  gkyl_array_release(udrift); 
  gkyl_array_release(vtsq);
  gkyl_proj_on_basis_release(proj_den);
  gkyl_proj_on_basis_release(proj_udrift);
  gkyl_proj_on_basis_release(proj_vtsq);

  gkyl_array_release(prim_moms_ho);
  gkyl_array_release(bmag); 
  gkyl_array_release(jacob_tot);
  if (use_gpu) {
    gkyl_array_release(prim_moms);
    gkyl_array_release(bmag_cu); 
    gkyl_array_release(jacob_tot_cu);
  }

  gkyl_array_release(distf);
  if (use_gpu)
    gkyl_array_release(distf_cu);
  gkyl_proj_maxwellian_on_basis_release(proj_max);
  gkyl_velocity_map_release(gvm);

#ifdef GKYL_HAVE_CUDA
  if (use_gpu) 
    gkyl_cu_free(basis_on_dev);
    gkyl_cu_free(conf_basis_on_dev);
#endif  
}

void test_1x2v_p1_gk() { test_1x2v_gk(1, false); }
void test_1x2v_p2_gk() { test_1x2v_gk(2, false); }

#ifdef GKYL_HAVE_CUDA
void test_1x2v_p1_gk_gpu() { test_1x2v_gk(1, true); }
void test_1x2v_p2_gk_gpu() { test_1x2v_gk(2, true); }
#endif

TEST_LIST = {
  { "test_1x2v_p1_gk", test_1x2v_p1_gk },
  { "test_1x2v_p2_gk", test_1x2v_p2_gk },

#ifdef GKYL_HAVE_CUDA
  { "test_1x2v_p1_gk_gpu", test_1x2v_p1_gk_gpu },
  { "test_1x2v_p2_gk_gpu", test_1x2v_p2_gk_gpu },
#endif
  { NULL, NULL },
};
