// A test of calculation of moments of a distribution function.
//
#include <acutest.h>
#include <math.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_mom_bcorr.h>
#include <gkyl_lbo_mom_bcorr.h>
#include <gkyl_lbo_mom_bcorr_priv.h>
#include <gkyl_array_rio.h>
#include <gkyl_mom_calc.h>
#include <gkyl_vlasov_mom.h>
#include <gkyl_prim_lbo.h>
#include <gkyl_prim_lbo_calc.h>
#include <gkyl_prim_lbo_vlasov.h>

inline double
maxwellian1D(double n, double vx, double ux, double vth)
{
  double v2 = (vx-ux)*(vx-ux);
  return n/sqrt(2*M_PI*vth*vth)*exp(-v2/(2*vth*vth));
}

inline double
maxwellian2D(double n, double vx, double vy, double ux, double uy, double vth)
{
  double v2 = (vx-ux)*(vx-ux) + (vy-uy)*(vy-uy);
  return n/(2*M_PI*vth*vth)*exp(-v2/(2*vth*vth));
}

void
evalDistFunc1x1v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], vx = xn[1];
  
  fout[0] = maxwellian1D(1.0, vx, 0.0, 1.0);
}

void
evalDistFunc1x2v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], vx = xn[1], vy = xn[2];
  
  fout[0] = maxwellian2D(1.0, vx, vy, 0.0, 0.0, 1.0);
}

// allocate array (filled with zeros)
static struct gkyl_array*
mkarr(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}
// allocate cu_dev array
static struct gkyl_array*
mkarr_cu(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
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

void
test_1x1v_p2()
{
  int poly_order = 2;
  double lower[] = {0.0, -2.0}, upper[] = {1.0, 2.0};
  int cells[] = {4, 24};
  int ndim = sizeof(lower)/sizeof(lower[0]);
  int vdim = 1, cdim = 1;

  double v_bounds[] = {lower[1], upper[1]};
  int v_edge_idx[] = {cells[1]-1};

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  int confCells[] = {cells[0]};

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

  // basis functions
  struct gkyl_basis basis, confBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  int confGhost[] = { 0 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);
  struct skin_ghost_ranges confSkin_ghost; // conf-space skin/ghost
  skin_ghost_ranges_init(&confSkin_ghost, &confLocal_ext, confGhost);

  int ghost[] = { confGhost[0], 0, 0};
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  struct skin_ghost_ranges skin_ghost; // phase-space skin/ghost
  skin_ghost_ranges_init(&skin_ghost, &local_ext, ghost);

  // projection updater for dist-function
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, evalDistFunc1x1v, NULL);

  // create distribution function array
  struct gkyl_array *distf;
  distf = mkarr(basis.num_basis, local_ext.volume);

  // project distribution function on basis
  gkyl_proj_on_basis_advance(projDistf, 0.0, &local, distf);
  
  struct gkyl_mom_type *vmM0_t = gkyl_vlasov_mom_new(&confBasis, &basis, "M0");
  struct gkyl_mom_type *vmM1i_t = gkyl_vlasov_mom_new(&confBasis, &basis, "M1i");
  struct gkyl_mom_type *vmM2_t = gkyl_vlasov_mom_new(&confBasis, &basis, "M2");
  gkyl_mom_calc *m0calc = gkyl_mom_calc_new(&grid, vmM0_t);
  gkyl_mom_calc *m1icalc = gkyl_mom_calc_new(&grid, vmM1i_t);
  gkyl_mom_calc *m2calc = gkyl_mom_calc_new(&grid, vmM2_t);

  // create moment arrays
  struct gkyl_array *m0, *m1i, *m2;
  m0 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  m1i = mkarr(vdim*confBasis.num_basis, confLocal_ext.volume);
  m2 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  
  // compute the moments
  gkyl_mom_calc_advance(m0calc, local, confLocal, distf, m0);
  gkyl_mom_calc_advance(m1icalc, local, confLocal, distf, m1i);
  gkyl_mom_calc_advance(m2calc, local, confLocal, distf, m2);

  struct gkyl_mom_type *F = gkyl_vlasov_lbo_mom_new(&confBasis, &basis, "f", v_bounds);
  struct gkyl_mom_type *VF = gkyl_vlasov_lbo_mom_new(&confBasis, &basis, "vf", v_bounds);

  gkyl_mom_bcorr *fcalc = gkyl_mom_bcorr_new(&grid, F);
  gkyl_mom_bcorr *vFcalc = gkyl_mom_bcorr_new(&grid, VF);
  
  // create moment arrays
  struct gkyl_array *f, *vf;
  f = mkarr(2*confBasis.num_basis, confLocal_ext.volume);
  vf = mkarr(2*confBasis.num_basis, confLocal_ext.volume);

  // compute the moment corrections
  gkyl_mom_bcorr_advance(fcalc, local, confLocal, distf, f);
  gkyl_mom_bcorr_advance(vFcalc, local, confLocal, distf, vf);
  
  struct gkyl_prim_lbo *prim = gkyl_prim_lbo_vlasov_new(&confBasis, &basis);

  TEST_CHECK( prim->cdim == 1 );
  TEST_CHECK( prim->pdim == 2 );
  TEST_CHECK( prim->poly_order == 2 );
  TEST_CHECK( prim->num_config == confBasis.num_basis );
  TEST_CHECK( prim->num_phase == basis.num_basis );

  gkyl_prim_lbo_calc *primcalc = gkyl_prim_lbo_calc_new(&grid, prim);
  
  // create moment arrays
  struct gkyl_array *u, *vth;
  u = mkarr(confBasis.num_basis, confLocal_ext.volume);
  vth = mkarr(confBasis.num_basis, confLocal_ext.volume);

  // compute the moment corrections
  gkyl_prim_lbo_calc_advance(primcalc, confBasis, confLocal, m0, m1i, m2, f, vf, u, vth);

  // Check u
  for (unsigned int i=0; i<cells[0]; ++i) {
    int cidx[] = {i};
    long linc = gkyl_range_idx(&confLocal, cidx);
    double *uptr = gkyl_array_fetch(u, linc);
    for (unsigned int k=0; k<confBasis.num_basis; ++k) {
      TEST_CHECK( gkyl_compare( 0., uptr[k], 1e-12) );
  }}

  // Check vtSq.
  for (unsigned int i=0; i<cells[0]; ++i) {
    int cidx[] = {i};
    long linc = gkyl_range_idx(&confLocal, cidx);
    double *vthptr = gkyl_array_fetch(vth, linc);
    TEST_CHECK( gkyl_compare( 1.4142398195471544, vthptr[0], 1e-12) );
    for (unsigned int k=1; k<confBasis.num_basis; ++k) {
      TEST_CHECK( gkyl_compare( 0., vthptr[k], 1e-12) );
  }}

  // release memory for objects
  gkyl_array_release(m0); gkyl_array_release(m1i); gkyl_array_release(m2);
  gkyl_mom_calc_release(m0calc); gkyl_mom_calc_release(m1icalc); gkyl_mom_calc_release(m2calc);
  gkyl_mom_type_release(vmM0_t); gkyl_mom_type_release(vmM1i_t); gkyl_mom_type_release(vmM2_t);

  gkyl_array_release(f); gkyl_array_release(vf);
  gkyl_mom_bcorr_release(fcalc); gkyl_mom_bcorr_release(vFcalc);
  gkyl_mom_type_release(F); gkyl_mom_type_release(VF);

  gkyl_proj_on_basis_release(projDistf);
  gkyl_array_release(distf);

  gkyl_array_release(u); gkyl_array_release(vth);
  gkyl_prim_lbo_release(prim);
  gkyl_prim_lbo_calc_release(primcalc);
}

#ifdef GKYL_HAVE_CUDA
void
test_1x1v_p2_cu()
{
  int poly_order = 2;
  double lower[] = {0.0, -2.0}, upper[] = {1.0, 2.0};
  int cells[] = {4, 24};
  int ndim = sizeof(lower)/sizeof(lower[0]);
  int vdim = 1, cdim = 1;

  double v_bounds[] = {lower[1], upper[1]};
  int v_edge_idx[] = {cells[1]-1};

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  int confCells[] = {cells[0]};
  
  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

  // basis functions
  struct gkyl_basis basis, confBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  int confGhost[] = { 0 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);
  struct skin_ghost_ranges confSkin_ghost; // conf-space skin/ghost
  skin_ghost_ranges_init(&confSkin_ghost, &confLocal_ext, confGhost);

  int ghost[] = { confGhost[0], 0, 0};
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  struct skin_ghost_ranges skin_ghost; // phase-space skin/ghost
  skin_ghost_ranges_init(&skin_ghost, &local_ext, ghost);

  // projection updater for dist-function
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, evalDistFunc1x1v, NULL);

  // create distribution function array
  struct gkyl_array *distf, *distf_cu;
  distf = mkarr(basis.num_basis, local_ext.volume);
  distf_cu = mkarr_cu(basis.num_basis, local_ext.volume);

  // project distribution function on basis
  gkyl_proj_on_basis_advance(projDistf, 0.0, &local, distf);
  gkyl_array_copy(distf_cu, distf);
  
  struct gkyl_mom_type *vmM0_t = gkyl_vlasov_mom_cu_dev_new(&confBasis, &basis, "M0");
  struct gkyl_mom_type *vmM1i_t = gkyl_vlasov_mom_cu_dev_new(&confBasis, &basis, "M1i");
  struct gkyl_mom_type *vmM2_t = gkyl_vlasov_mom_cu_dev_new(&confBasis, &basis, "M2");
  gkyl_mom_calc *m0calc = gkyl_mom_calc_cu_dev_new(&grid, vmM0_t);
  gkyl_mom_calc *m1icalc = gkyl_mom_calc_cu_dev_new(&grid, vmM1i_t);
  gkyl_mom_calc *m2calc = gkyl_mom_calc_cu_dev_new(&grid, vmM2_t);

  // create moment arrays
  struct gkyl_array *m0_cu, *m1i_cu, *m2_cu;
  m0_cu = mkarr_cu(confBasis.num_basis, confLocal_ext.volume);
  m1i_cu = mkarr_cu(vdim*confBasis.num_basis, confLocal_ext.volume);
  m2_cu = mkarr_cu(confBasis.num_basis, confLocal_ext.volume);
  
  // compute the moment and copy back to host
  gkyl_mom_calc_advance_cu(m0calc, local, confLocal, distf_cu, m0_cu);
  gkyl_mom_calc_advance_cu(m1icalc, local, confLocal, distf_cu, m1i_cu);
  gkyl_mom_calc_advance_cu(m2calc, local, confLocal, distf_cu, m2_cu);

  struct gkyl_mom_type *F = gkyl_vlasov_lbo_mom_cu_dev_new(&confBasis, &basis, "f", v_bounds);
  struct gkyl_mom_type *VF = gkyl_vlasov_lbo_mom_cu_dev_new(&confBasis, &basis, "vf", v_bounds);

  gkyl_mom_bcorr *fcalc = gkyl_mom_bcorr_cu_dev_new(&grid, F);
  gkyl_mom_bcorr *vFcalc = gkyl_mom_bcorr_cu_dev_new(&grid, VF);
  
  // create moment arrays
  struct gkyl_array *f_cu, *vf_cu;
  f_cu = mkarr_cu(2*confBasis.num_basis, confLocal_ext.volume);
  vf_cu = mkarr_cu(2*confBasis.num_basis, confLocal_ext.volume);

  // compute the moment corrections
  gkyl_mom_bcorr_advance_cu(fcalc, local, confLocal, distf_cu, f_cu);
  gkyl_mom_bcorr_advance_cu(vFcalc, local, confLocal, distf_cu, vf_cu);
  
  struct gkyl_prim_lbo *prim = gkyl_prim_lbo_vlasov_cu_dev_new(&confBasis, &basis);

  gkyl_prim_lbo_calc *primcalc = gkyl_prim_lbo_calc_cu_dev_new(&grid, prim);
  
  // create moment arrays
  struct gkyl_array *u, *vth, *u_cu, *vth_cu;
  u = mkarr(confBasis.num_basis, confLocal_ext.volume);
  vth = mkarr(confBasis.num_basis, confLocal_ext.volume);
  u_cu = mkarr_cu(confBasis.num_basis, confLocal_ext.volume);
  vth_cu = mkarr_cu(confBasis.num_basis, confLocal_ext.volume);

  // compute the moment corrections
  gkyl_prim_lbo_calc_advance_cu(primcalc, confBasis, confLocal, m0_cu, m1i_cu, m2_cu, f_cu, vf_cu, u_cu, vth_cu);
  gkyl_array_copy(u, u_cu);
  gkyl_array_copy(vth, vth_cu);

  // Check u
  for (unsigned int i=0; i<cells[0]; ++i) {
    int cidx[] = {i};
    long linc = gkyl_range_idx(&confLocal, cidx);
    double *uptr = gkyl_array_fetch(u, linc);
    for (unsigned int k=0; k<confBasis.num_basis; ++k) {
      TEST_CHECK( gkyl_compare( 0., uptr[k], 1e-12) );
  }}

  // Check vtSq.
  for (unsigned int i=0; i<cells[0]; ++i) {
    int cidx[] = {i};
    long linc = gkyl_range_idx(&confLocal, cidx);
    double *vthptr = gkyl_array_fetch(vth, linc);
    TEST_CHECK( gkyl_compare( 1.4142398195471544, vthptr[0], 1e-12) );
    for (unsigned int k=1; k<confBasis.num_basis; ++k) {
      TEST_CHECK( gkyl_compare( 0., vthptr[k], 1e-12) );
  }}

  // release memory for objects
  gkyl_array_release(m0_cu); gkyl_array_release(m1i_cu); gkyl_array_release(m2_cu);
  gkyl_mom_calc_release(m0calc); gkyl_mom_calc_release(m1icalc); gkyl_mom_calc_release(m2calc);
  gkyl_mom_type_release(vmM0_t); gkyl_mom_type_release(vmM1i_t); gkyl_mom_type_release(vmM2_t);
  
  gkyl_array_release(f_cu); gkyl_array_release(vf_cu);
  gkyl_mom_bcorr_release(fcalc); gkyl_mom_bcorr_release(vFcalc);
  gkyl_mom_type_release(F); gkyl_mom_type_release(VF);

  gkyl_proj_on_basis_release(projDistf);
  gkyl_array_release(distf); gkyl_array_release(distf_cu);

  gkyl_array_release(u); gkyl_array_release(vth);
  gkyl_array_release(u_cu); gkyl_array_release(vth_cu);

  gkyl_prim_lbo_release(prim);
  gkyl_prim_lbo_calc_release(primcalc);
}
#endif

TEST_LIST = {
  { "test_1x1v_p2", test_1x1v_p2 },
#ifdef GKYL_HAVE_CUDA
  { "test_1x1v_p2_cu", test_1x1v_p2_cu },
#endif
  { NULL, NULL },
};
