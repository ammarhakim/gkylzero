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

  // Write the initial distribution array to file.
  const char *fmt = "%s-%s_%d.gkyl";
  char name[] = "bcorr_1x1v";
  char momName[] = "distf";
  int frame = 0;
  int sz = snprintf(0, 0, fmt, name, momName, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, name, momName, frame);
  gkyl_grid_sub_array_write(&grid, &local, distf, fileNm);

  struct gkyl_mom_type *F = gkyl_vlasov_lbo_mom_new(&confBasis, &basis, "f");
  struct gkyl_mom_type *VF = gkyl_vlasov_lbo_mom_new(&confBasis, &basis, "vf");
  gkyl_lbo_mom_set_atLower(F, v_edge_idx);
  gkyl_lbo_mom_set_vBoundary(F, v_bounds);
  gkyl_lbo_mom_set_atLower(VF, v_edge_idx);
  gkyl_lbo_mom_set_vBoundary(VF, v_bounds);
  
  gkyl_mom_bcorr *fcalc = gkyl_mom_bcorr_new(&grid, F);
  gkyl_mom_bcorr *vFcalc = gkyl_mom_bcorr_new(&grid, VF);
  
  // create moment arrays
  struct gkyl_array *f, *vf;
  f = mkarr(2*confBasis.num_basis, confLocal_ext.volume);
  vf = mkarr(2*confBasis.num_basis, confLocal_ext.volume);

  // compute the moments
  gkyl_mom_bcorr_advance(fcalc, local, confLocal, distf, f);
  gkyl_mom_bcorr_advance(vFcalc, local, confLocal, distf, vf);

   // Write the initial distribution array to file.
  char momName1[] = "f";
  sz = snprintf(0, 0, fmt, name, momName1, frame);
  fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, name, momName1, frame);
  gkyl_grid_sub_array_write(&grid, &local, f, fileNm);

  // Write the initial distribution array to file.
  char momName2[] = "vf";
  sz = snprintf(0, 0, fmt, name, momName2, frame);
  fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, name, momName2, frame);
  gkyl_grid_sub_array_write(&grid, &local, vf, fileNm);
 
  // Check M0.
  for (unsigned int i=0; i<cells[0]; ++i) {
    int cidx[] = {i};
    long linc = gkyl_range_idx(&confLocal, cidx);
    double *fptr = gkyl_array_fetch(f, linc);
    TEST_CHECK( gkyl_compare( 0.07635960492981765, fptr[0], 1e-12) );
    for (unsigned int k=1; k<confBasis.num_basis; ++k) {
      TEST_CHECK( gkyl_compare( 0., fptr[k], 1e-12) );
  }}

  // Check M1i.
  for (unsigned int i=0; i<cells[0]; ++i) {
    int cidx[] = {i};
    long linc = gkyl_range_idx(&confLocal, cidx);
    double *vfptr = gkyl_array_fetch(vf, linc);
    TEST_CHECK( gkyl_compare( 0.1527192098596353, vfptr[0], 1e-12) );
    for (unsigned int k=1; k<vf->ncomp; ++k) {
      TEST_CHECK( gkyl_compare( 0., vfptr[k], 1e-12) );
  }}

  // release memory for moment data object
  gkyl_array_release(f); gkyl_array_release(vf);
  gkyl_mom_bcorr_release(fcalc); gkyl_mom_bcorr_release(vFcalc);
  gkyl_mom_type_release(F); gkyl_mom_type_release(VF);

  gkyl_proj_on_basis_release(projDistf);
  gkyl_array_release(distf);
}

void
test_1x2v_p2()
{
  int poly_order = 2;
  double lower[] = {0.0, -2.0, -2.0}, upper[] = {1.0, 2.0, 2.0};
  int cells[] = {4, 24, 24};
  int ndim = sizeof(lower)/sizeof(lower[0]);
  int vdim = 2, cdim = 1;

  double v_bounds[] = {lower[1], lower[2], upper[1], upper[2]};
  int v_edge_idx[] = {cells[1]-1, cells[2]-1};

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

  // Write the initial distribution array to file.
  const char *fmt = "%s-%s_%d.gkyl";
  char name[] = "bcorr_1x2v";
  char momName[] = "distf";
  int frame = 0;
  int sz = snprintf(0, 0, fmt, name, momName, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, name, momName, frame);
  gkyl_grid_sub_array_write(&grid, &local, distf, fileNm);

  struct gkyl_mom_type *F = gkyl_vlasov_lbo_mom_new(&confBasis, &basis, "f");
  struct gkyl_mom_type *VF = gkyl_vlasov_lbo_mom_new(&confBasis, &basis, "vf");
  gkyl_lbo_mom_set_atLower(F, v_edge_idx);
  gkyl_lbo_mom_set_vBoundary(F, v_bounds);
  gkyl_lbo_mom_set_atLower(VF, v_edge_idx);
  gkyl_lbo_mom_set_vBoundary(VF, v_bounds);
  
  gkyl_mom_bcorr *fcalc = gkyl_mom_bcorr_new(&grid, F);
  gkyl_mom_bcorr *vFcalc = gkyl_mom_bcorr_new(&grid, VF);
  
  // create moment arrays
  struct gkyl_array *f, *vf;
  f = mkarr(2*confBasis.num_basis, confLocal_ext.volume);
  vf = mkarr(2*confBasis.num_basis, confLocal_ext.volume);

  // compute the moments
  gkyl_mom_bcorr_advance(fcalc, local, confLocal, distf, f);
  gkyl_mom_bcorr_advance(vFcalc, local, confLocal, distf, vf);

  // Write the initial distribution array to file.
  char momName1[] = "f";
  sz = snprintf(0, 0, fmt, name, momName1, frame);
  fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, name, momName1, frame);
  gkyl_grid_sub_array_write(&grid, &local, f, fileNm);

  // Write the initial distribution array to file.
  char momName2[] = "vf";
  sz = snprintf(0, 0, fmt, name, momName2, frame);
  fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, name, momName2, frame);
  gkyl_grid_sub_array_write(&grid, &local, vf, fileNm);
 
  // Check M0.
  for (unsigned int i=0; i<cells[0]; ++i) {
    int cidx[] = {i};
    long linc = gkyl_range_idx(&confLocal, cidx);
    double *fptr = gkyl_array_fetch(f, linc);
    TEST_CHECK( gkyl_compare( 0.305438419719271, fptr[0], 1e-12) );
    TEST_CHECK( gkyl_compare( 0.0, fptr[1], 1e-12) );
    TEST_CHECK( gkyl_compare( 0.0, fptr[2], 1e-12) );
    TEST_CHECK( gkyl_compare( 1.34986647204995, fptr[3], 1e-12) );
    TEST_CHECK( gkyl_compare( 0.0, fptr[4], 1e-12) );
    TEST_CHECK( gkyl_compare( 0.0, fptr[5], 1e-12) );
  }

  // Check M1i.
  for (unsigned int i=0; i<cells[0]; ++i) {
    int cidx[] = {i};
    long linc = gkyl_range_idx(&confLocal, cidx);
    double *vfptr = gkyl_array_fetch(vf, linc);
    TEST_CHECK( gkyl_compare( 3.31060978353845, vfptr[0], 1e-12) );
    TEST_CHECK( gkyl_compare( 0.0, vfptr[1], 1e-12) );
    TEST_CHECK( gkyl_compare( 0.0, vfptr[2], 1e-12) );
    TEST_CHECK( gkyl_compare( 0.0, vfptr[3], 1e-12) );
    TEST_CHECK( gkyl_compare( 0.0, vfptr[4], 1e-12) );
    TEST_CHECK( gkyl_compare( 0.0, vfptr[5], 1e-12) );
    }

  // release memory for moment data object
  gkyl_array_release(f); gkyl_array_release(vf);
  gkyl_mom_bcorr_release(fcalc); gkyl_mom_bcorr_release(vFcalc);
  gkyl_mom_type_release(F); gkyl_mom_type_release(VF);

  gkyl_proj_on_basis_release(projDistf);
  gkyl_array_release(distf);
}

TEST_LIST = {
  { "test_1x1v_p2", test_1x1v_p2 },
  { "test_1x2v_p2", test_1x2v_p2 },
  { NULL, NULL },
};
