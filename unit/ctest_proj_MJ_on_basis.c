#include "gkyl_array.h"
#include "gkyl_util.h"
#include <acutest.h>

#include <gkyl_array_rio.h>
#include <gkyl_correct_maxwellian.h>
#include <gkyl_correct_MJ.h>
#include <gkyl_MJ_moments.h>
#include <gkyl_proj_MJ_on_basis.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <math.h>

// allocate array (filled with zeros)
static struct gkyl_array*
mkarr(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

// Projection functions for p/(gamma) = v in special relativistic systems
// Simplifies to p/sqrt(1 + p^2) where c = 1
static void
ev_p_over_gamma_1p(double t, const double *xn, double *out, void *ctx)
{
  out[0] = xn[0]/sqrt(1.0 + xn[0]*xn[0]);
}
static void
ev_p_over_gamma_2p(double t, const double *xn, double *out, void *ctx)
{
  out[0] = xn[0]/sqrt(1.0 + xn[0]*xn[0] + xn[1]*xn[1]);
  out[1] = xn[1]/sqrt(1.0 + xn[0]*xn[0] + xn[1]*xn[1]);
}
static void
ev_p_over_gamma_3p(double t, const double *xn, double *out, void *ctx)
{
  out[0] = xn[0]/sqrt(1.0 + xn[0]*xn[0] + xn[1]*xn[1] + xn[2]*xn[2]);
  out[1] = xn[1]/sqrt(1.0 + xn[0]*xn[0] + xn[1]*xn[1] + xn[2]*xn[2]);
  out[2] = xn[2]/sqrt(1.0 + xn[0]*xn[0] + xn[1]*xn[1] + xn[2]*xn[2]);
}

static const evalf_t p_over_gamma_func[3] = {ev_p_over_gamma_1p, ev_p_over_gamma_2p, ev_p_over_gamma_3p};

// Projection functions for gamma = sqrt(1 + p^2) in special relativistic systems
static void
ev_gamma_1p(double t, const double *xn, double *out, void *ctx)
{
  out[0] = sqrt(1.0 + xn[0]*xn[0]);
}
static void
ev_gamma_2p(double t, const double *xn, double *out, void *ctx)
{
  out[0] = sqrt(1.0 + xn[0]*xn[0] + xn[1]*xn[1]);
}
static void
ev_gamma_3p(double t, const double *xn, double *out, void *ctx)
{
  out[0] = sqrt(1.0 + xn[0]*xn[0] + xn[1]*xn[1] + xn[2]*xn[2]);
}

static const evalf_t gamma_func[3] = {ev_gamma_1p, ev_gamma_2p, ev_gamma_3p};

// Projection functions for gamma_inv = 1/sqrt(1 + p^2) in special relativistic systems
static void
ev_gamma_inv_1p(double t, const double *xn, double *out, void *ctx)
{
  out[0] = 1.0/sqrt(1.0 + xn[0]*xn[0]);
}
static void
ev_gamma_inv_2p(double t, const double *xn, double *out, void *ctx)
{
  out[0] = 1.0/sqrt(1.0 + xn[0]*xn[0] + xn[1]*xn[1]);
}
static void
ev_gamma_inv_3p(double t, const double *xn, double *out, void *ctx)
{
  out[0] = 1.0/sqrt(1.0 + xn[0]*xn[0] + xn[1]*xn[1] + xn[2]*xn[2]);
}

static const evalf_t gamma_inv_func[3] = {ev_gamma_inv_1p, ev_gamma_inv_2p, ev_gamma_inv_3p};



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

void eval_M0(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 1.0;
}

void eval_M1i_1v_no_drift(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 0.0;
}

void eval_M2_1v_no_drift(double t, const double *xn, double* restrict fout, void *ctx)
{
  double T = 1.0;
  double x = xn[0];
  fout[0] = T;
}

void eval_M1i_1v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 0.5; //0.5;
}

void eval_M2_1v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double T = 1.0;
  double x = xn[0];
  fout[0] = T;
}

void eval_M1i_2v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 0.5; fout[1] = 0.25;
}

void eval_M2_2v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double T = 1.0;
  double x = xn[0];
  fout[0] = T;
}

void eval_M1i_3v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 0.5; fout[1] = 0.25; fout[2] = -0.5;
}

void eval_M2_3v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double T = 1.0;
  double x = xn[0];
  fout[0] = T;
}

void
test_1x1v_no_drift(int poly_order)
{
  double lower[] = {0.1, -15.0}, upper[] = {1.0, 15.0};
  int cells[] = {2, 32};
  int vdim = 1, cdim = 1;
  int ndim = cdim+vdim;

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  double velLower[] = {lower[1]}, velUpper[] = {upper[1]};
  int confCells[] = {cells[0]};
  int velCells[] = {cells[1]};

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);
  struct gkyl_rect_grid vel_grid;
  gkyl_rect_grid_init(&vel_grid, vdim, velLower, velUpper, velCells);

  // velocity range
  int velGhost[] = { 0 };
  struct gkyl_range velLocal, velLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&vel_grid, velGhost, &velLocal_ext, &velLocal);

  // basis functions
  struct gkyl_basis basis, confBasis, velBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);
  gkyl_cart_modal_serendip(&velBasis, vdim, poly_order);

  int confGhost[] = { 1 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);
  struct skin_ghost_ranges confSkin_ghost; // conf-space skin/ghost
  skin_ghost_ranges_init(&confSkin_ghost, &confLocal_ext, confGhost);

  int ghost[] = { confGhost[0], 0 };
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  struct skin_ghost_ranges skin_ghost; // phase-space skin/ghost
  skin_ghost_ranges_init(&skin_ghost, &local_ext, ghost);

  // create moment arrays
  struct gkyl_array *m0, *m1i, *m2;
  m0 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  m1i = mkarr(vdim*confBasis.num_basis, confLocal_ext.volume);
  m2 = mkarr(confBasis.num_basis, confLocal_ext.volume);

  gkyl_proj_on_basis *proj_m0 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_M0, NULL);
  gkyl_proj_on_basis *proj_m1i = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, vdim, eval_M1i_1v_no_drift, NULL);
  gkyl_proj_on_basis *proj_m2 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_M2_1v_no_drift, NULL);

  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0);
  gkyl_proj_on_basis_advance(proj_m1i, 0.0, &confLocal, m1i);
  gkyl_proj_on_basis_advance(proj_m2, 0.0, &confLocal, m2);

  // create distribution function array
  struct gkyl_array *distf;
  distf = mkarr(basis.num_basis, local_ext.volume);

  // projection updater to compute MJ
  gkyl_proj_MJ_on_basis *proj_MJ = gkyl_proj_MJ_on_basis_new(&grid,
    &confBasis, &basis, poly_order+1);

  gkyl_proj_MJ_on_basis_fluid_stationary_frame_mom(proj_MJ, &local, &confLocal, m0, m1i, m2, distf);

// build p_over_gamma
struct gkyl_array *p_over_gamma;
p_over_gamma = mkarr(vdim*velBasis.num_basis, velLocal.volume);
gkyl_proj_on_basis *p_over_gamma_proj = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
    .grid = &vel_grid,
    .basis = &velBasis,
    .qtype = GKYL_GAUSS_LOBATTO_QUAD,
    .num_quad = 8,
    .num_ret_vals = vdim,
    .eval = p_over_gamma_func[vdim-1],
    .ctx = 0
  });
  gkyl_proj_on_basis_advance(p_over_gamma_proj, 0.0, &velLocal, p_over_gamma);


  // correct the MJ distribution m0 Moment
  gkyl_correct_MJ *corr_MJ = gkyl_correct_MJ_new(&grid,&confBasis,&basis,&confLocal,&velLocal,confLocal.volume,confLocal_ext.volume, false);
  gkyl_correct_MJ_fix(corr_MJ,p_over_gamma,distf,m0,m1i,&local,&confLocal);

  // values to compare  at index (1, 17) [remember, lower-left index is (1,1)]
  double p1_vals[] = {  5.3918752026566863e-01, -1.0910243387206232e-17,  -6.0196985297046972e-02,
    5.0006050167249552e-18 };
  double p2_vals[] = {  5.3922143701031633e-01,  -9.6625223288531320e-18,  -5.7898881215132203e-02,
    7.8842251929957589e-18, 1.9166441863144966e-17, -1.0173903909543560e-02,
    1.7916734900988946e-17, 1.4245174569363429e-18 };

  const double *fv = gkyl_array_cfetch(distf, gkyl_range_idx(&local_ext, (int[2]) { 1, 17 }));

  if (poly_order == 1) {
    for (int i=0; i<basis.num_basis; ++i){
      TEST_CHECK( gkyl_compare_double(p1_vals[i], fv[i], 1e-12) );
    }
  }

  if (poly_order == 2) {
    for (int i=0; i<basis.num_basis; ++i){
      TEST_CHECK( gkyl_compare_double(p2_vals[i], fv[i], 1e-12) );
    }
  }

  // write distribution function to file
  char fname[1024];
  sprintf(fname, "ctest_proj_MJ_on_basis_test_1x1v_p%d_no_drift.gkyl", poly_order);
  gkyl_grid_sub_array_write(&grid, &local, distf, fname);

  // release memory for moment data object
  gkyl_proj_on_basis_release(p_over_gamma_proj);
  gkyl_correct_MJ_release(corr_MJ);
  gkyl_array_release(m0); gkyl_array_release(m1i); gkyl_array_release(m2);
  gkyl_array_release(distf);
  gkyl_proj_MJ_on_basis_release(proj_MJ);
  gkyl_proj_on_basis_release(proj_m0);
  gkyl_proj_on_basis_release(proj_m1i);
  gkyl_proj_on_basis_release(proj_m2);
  gkyl_array_release(p_over_gamma);
}

void test_1x1v_no_drift_p1() { test_1x1v_no_drift(1); }
void test_1x1v_no_drift_p2() { test_1x1v_no_drift(2); }


void
test_1x1v(int poly_order)
{
  double lower[] = {0.1, -15.0}, upper[] = {1.0, 15.0}; // +/- 15 on velocity
  int cells[] = {2, 32};  // default {2, 32}
  int vdim = 1, cdim = 1;
  int ndim = cdim+vdim;

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  double velLower[] = {lower[1]}, velUpper[] = {upper[1]};
  int confCells[] = {cells[0]};
  int velCells[] = {cells[1]};

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);
  struct gkyl_rect_grid vel_grid;
  gkyl_rect_grid_init(&vel_grid, vdim, velLower, velUpper, velCells);

  // velocity range
  int velGhost[] = { 0 };
  struct gkyl_range velLocal, velLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&vel_grid, velGhost, &velLocal_ext, &velLocal);

  // basis functions
  struct gkyl_basis basis, confBasis, velBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);
  gkyl_cart_modal_serendip(&velBasis, vdim, poly_order);

  int confGhost[] = { 1 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);
  struct skin_ghost_ranges confSkin_ghost; // conf-space skin/ghost
  skin_ghost_ranges_init(&confSkin_ghost, &confLocal_ext, confGhost);

  int ghost[] = { confGhost[0], 0 };
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  struct skin_ghost_ranges skin_ghost; // phase-space skin/ghost
  skin_ghost_ranges_init(&skin_ghost, &local_ext, ghost);

  // create moment arrays
  struct gkyl_array *m0, *m1i, *m2;
  m0 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  m1i = mkarr(vdim*confBasis.num_basis, confLocal_ext.volume);
  m2 = mkarr(confBasis.num_basis, confLocal_ext.volume);

  gkyl_proj_on_basis *proj_m0 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_M0, NULL);
  gkyl_proj_on_basis *proj_m1i = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, vdim, eval_M1i_1v, NULL);
  gkyl_proj_on_basis *proj_m2 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_M2_1v, NULL);

  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0);
  gkyl_proj_on_basis_advance(proj_m1i, 0.0, &confLocal, m1i);
  gkyl_proj_on_basis_advance(proj_m2, 0.0, &confLocal, m2);

  // create distribution function array
  struct gkyl_array *distf;
  distf = mkarr(basis.num_basis, local_ext.volume);

  // projection updater to compute MJ
  gkyl_proj_MJ_on_basis *proj_MJ = gkyl_proj_MJ_on_basis_new(&grid,
    &confBasis, &basis, poly_order+1);

  gkyl_proj_MJ_on_basis_fluid_stationary_frame_mom(proj_MJ, &local, &confLocal, m0, m1i, m2, distf);

// build the p_over_gamma
struct gkyl_array *p_over_gamma;
p_over_gamma = mkarr(vdim*velBasis.num_basis, velLocal.volume);
gkyl_proj_on_basis *p_over_gamma_proj = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
    .grid = &vel_grid,
    .basis = &velBasis,
    .qtype = GKYL_GAUSS_LOBATTO_QUAD,
    .num_quad = 8,
    .num_ret_vals = vdim,
    .eval = p_over_gamma_func[vdim-1], //ev_p_over_gamma_1p,
    .ctx = 0
  });
  gkyl_proj_on_basis_advance(p_over_gamma_proj, 0.0, &velLocal, p_over_gamma);

  // build gamma
  struct gkyl_array *gamma;
  gamma = mkarr(velBasis.num_basis, velLocal.volume);
  gkyl_proj_on_basis *gamma_proj = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
      .grid = &vel_grid,
      .basis = &velBasis,
      .qtype = GKYL_GAUSS_LOBATTO_QUAD,
      .num_quad = 8,
      .num_ret_vals = 1,
      .eval = gamma_func[vdim-1],
      .ctx = 0
    });
    gkyl_proj_on_basis_advance(gamma_proj, 0.0, &velLocal, gamma);


  // build gamma_inv
  struct gkyl_array *gamma_inv;
  gamma_inv = mkarr(velBasis.num_basis, velLocal.volume);
  gkyl_proj_on_basis *gamma_inv_proj = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
      .grid = &vel_grid,
      .basis = &velBasis,
      .qtype = GKYL_GAUSS_LOBATTO_QUAD,
      .num_quad = 8,
      .num_ret_vals = 1,
      .eval = gamma_inv_func[vdim-1],
      .ctx = 0
    });
    gkyl_proj_on_basis_advance(gamma_inv_proj, 0.0, &velLocal, gamma_inv);

  // correct the MJ distribution m0 Moment
  gkyl_correct_MJ *corr_MJ = gkyl_correct_MJ_new(&grid,&confBasis,&basis,&confLocal,&velLocal,confLocal.volume,confLocal_ext.volume, false);
  gkyl_correct_MJ_fix(corr_MJ,p_over_gamma,distf,m0,m1i,&local,&confLocal);

  // test accuracy of the projection:
  gkyl_MJ_moments *MJ_moms = gkyl_MJ_moments_new(&grid,&confBasis,&basis,&confLocal,&velLocal,confLocal.volume,confLocal_ext.volume, false);
  gkyl_MJ_moments_advance(MJ_moms,p_over_gamma,gamma,gamma_inv,distf,m0,m1i,m2,&local,&confLocal);


  // values to compare  at index (1, 17) [remember, lower-left index is (1,1)]
  double p2_vals[] = {  5.9020018022791720e-01,  1.8856465819367569e-17,  1.6811060851198739e-02,
    -3.1835607256007792e-18, -2.5559407924922751e-17, -1.6328957375267440e-02,
    3.5362637935871408e-17, -1.1626136441969583e-17 };

  const double *fv = gkyl_array_cfetch(distf, gkyl_range_idx(&local_ext, (int[2]) { 1, 17 }));

  if (poly_order == 2) {
    for (int i=0; i<basis.num_basis; ++i){
      TEST_CHECK( gkyl_compare_double(p2_vals[i], fv[i], 1e-12) );
    }
  }

  // write distribution function to file
  char fname[1024];
  sprintf(fname, "ctest_proj_MJ_on_basis_test_1x1v_p%d.gkyl", poly_order);
  gkyl_grid_sub_array_write(&grid, &local, distf, fname);

  // release memory for moment data object
  gkyl_proj_on_basis_release(p_over_gamma_proj);
  gkyl_proj_on_basis_release(gamma_proj);
  gkyl_proj_on_basis_release(gamma_inv_proj);
  gkyl_correct_MJ_release(corr_MJ);
  gkyl_MJ_moments_release(MJ_moms);
  gkyl_array_release(m0); gkyl_array_release(m1i); gkyl_array_release(m2);
  gkyl_array_release(distf);
  gkyl_proj_MJ_on_basis_release(proj_MJ);
  gkyl_proj_on_basis_release(proj_m0);
  gkyl_proj_on_basis_release(proj_m1i);
  gkyl_proj_on_basis_release(proj_m2);
  gkyl_array_release(p_over_gamma);
}

// special note, the p1 basis does not function
void test_1x1v_p2() { test_1x1v(2); }

void
test_1x2v(int poly_order)
{
  double lower[] = {0.1, -15.0, -15.0}, upper[] = {1.0, 15.0, 15.0};
  int cells[] = {2, 16, 16};
  int vdim = 2, cdim = 1;
  int ndim = cdim+vdim;

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  double velLower[] = {lower[1], lower[2]}, velUpper[] = {upper[1], upper[2]};
  int confCells[] = {cells[0]};
  int velCells[] = {cells[1], cells[2]};

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);
  struct gkyl_rect_grid vel_grid;
  gkyl_rect_grid_init(&vel_grid, vdim, velLower, velUpper, velCells);

  // velocity range
  int velGhost[] = { 0, 0 };
  struct gkyl_range velLocal, velLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&vel_grid, velGhost, &velLocal_ext, &velLocal);

  // basis functions
  struct gkyl_basis basis, confBasis, velBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);
  gkyl_cart_modal_serendip(&velBasis, vdim, poly_order);

  int confGhost[] = { 1 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);
  struct skin_ghost_ranges confSkin_ghost; // conf-space skin/ghost
  skin_ghost_ranges_init(&confSkin_ghost, &confLocal_ext, confGhost);

  int ghost[] = { confGhost[0], 0 , 0};
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  struct skin_ghost_ranges skin_ghost; // phase-space skin/ghost
  skin_ghost_ranges_init(&skin_ghost, &local_ext, ghost);

  // create moment arrays
  struct gkyl_array *m0, *m1i, *m2;
  m0 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  m1i = mkarr(vdim*confBasis.num_basis, confLocal_ext.volume);
  m2 = mkarr(confBasis.num_basis, confLocal_ext.volume);

  gkyl_proj_on_basis *proj_m0 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_M0, NULL);
  gkyl_proj_on_basis *proj_m1i = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, vdim, eval_M1i_2v, NULL);
  gkyl_proj_on_basis *proj_m2 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_M2_2v, NULL);

  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0);
  gkyl_proj_on_basis_advance(proj_m1i, 0.0, &confLocal, m1i);
  gkyl_proj_on_basis_advance(proj_m2, 0.0, &confLocal, m2);

  // create distribution function array
  struct gkyl_array *distf;
  distf = mkarr(basis.num_basis, local_ext.volume);

  // projection updater to compute MJ
  gkyl_proj_MJ_on_basis *proj_MJ = gkyl_proj_MJ_on_basis_new(&grid,
    &confBasis, &basis, poly_order+1);

  gkyl_proj_MJ_on_basis_fluid_stationary_frame_mom(proj_MJ, &local, &confLocal, m0, m1i, m2, distf);

// build the p_over_gamma
struct gkyl_array *p_over_gamma;
p_over_gamma = mkarr(vdim*velBasis.num_basis, velLocal.volume);
gkyl_proj_on_basis *p_over_gamma_proj = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
    .grid = &vel_grid,
    .basis = &velBasis,
    .qtype = GKYL_GAUSS_LOBATTO_QUAD,
    .num_quad = 8,
    .num_ret_vals = vdim,
    .eval = p_over_gamma_func[vdim-1],
    .ctx = 0
  });
  gkyl_proj_on_basis_advance(p_over_gamma_proj, 0.0, &velLocal, p_over_gamma);

  // correct the MJ distribution m0 Moment
  gkyl_correct_MJ *corr_MJ = gkyl_correct_MJ_new(&grid,&confBasis,&basis,&confLocal,&velLocal,confLocal.volume,confLocal_ext.volume, false);
  gkyl_correct_MJ_fix(corr_MJ,p_over_gamma,distf,m0,m1i,&local,&confLocal);

  // values to compare  at index (1, 9, 9) [remember, lower-left index is (1,1,1)]
  double p2_vals[] = { 1.6962401858070214e-01, -7.7409013792616062e-18, -3.9380955867889663e-03,
  -2.8639115438725903e-02, 5.5535324575725082e-19, 1.1489034695079269e-18,
  8.1594568989774036e-03, -1.2454319344697310e-17, -1.0858276718973473e-02,
  -8.3719428179921188e-03, -5.6169974218438426e-19, 1.1852358453425043e-17,
  1.5733767456010547e-18, 1.2600839078583602e-17, 2.0991499099505761e-03,
  9.9534173931022546e-19, -6.4527265157546547e-04, -1.0507606385517198e-17,
  -9.2457560036345989e-19, 7.4145989533791131e-19 };

  const double *fv = gkyl_array_cfetch(distf, gkyl_range_idx(&local_ext, (int[3]) { 1, 9, 9 }));

  if (poly_order == 2) {
    for (int i=0; i<basis.num_basis; ++i){
      TEST_CHECK( gkyl_compare_double(p2_vals[i], fv[i], 1e-12) );
    }
  }

  // write distribution function to file
  char fname[1024];
  sprintf(fname, "ctest_proj_MJ_on_basis_test_1x2v_p%d.gkyl", poly_order);
  gkyl_grid_sub_array_write(&grid, &local, distf, fname);

  // release memory for moment data object
  gkyl_proj_on_basis_release(p_over_gamma_proj);
  gkyl_correct_MJ_release(corr_MJ);
  gkyl_array_release(m0); gkyl_array_release(m1i); gkyl_array_release(m2);
  gkyl_array_release(distf);
  gkyl_proj_MJ_on_basis_release(proj_MJ);
  gkyl_proj_on_basis_release(proj_m0);
  gkyl_proj_on_basis_release(proj_m1i);
  gkyl_proj_on_basis_release(proj_m2);
  gkyl_array_release(p_over_gamma);
}

void test_1x2v_p2() { test_1x2v(2); }



void
test_1x3v(int poly_order)
{
  double lower[] = {0.1, -15.0, -15.0, -15.0}, upper[] = {1.0, 15.0, 15.0, 15.0};
  int cells[] = {2, 16, 16, 16};
  int vdim = 3, cdim = 1;
  int ndim = cdim+vdim;

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  double velLower[] = {lower[1], lower[2], lower[3]}, velUpper[] = {upper[1], upper[2], upper[3]};
  int confCells[] = {cells[0]};
  int velCells[] = {cells[1], cells[2], cells[3]};

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);
  struct gkyl_rect_grid vel_grid;
  gkyl_rect_grid_init(&vel_grid, vdim, velLower, velUpper, velCells);

  // velocity range
  int velGhost[] = { 0, 0, 0 };
  struct gkyl_range velLocal, velLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&vel_grid, velGhost, &velLocal_ext, &velLocal);

  // basis functions
  struct gkyl_basis basis, confBasis, velBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);
  gkyl_cart_modal_serendip(&velBasis, vdim, poly_order);

  int confGhost[] = { 1 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);
  struct skin_ghost_ranges confSkin_ghost; // conf-space skin/ghost
  skin_ghost_ranges_init(&confSkin_ghost, &confLocal_ext, confGhost);

  int ghost[] = { confGhost[0], 0 , 0, 0};
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  struct skin_ghost_ranges skin_ghost; // phase-space skin/ghost
  skin_ghost_ranges_init(&skin_ghost, &local_ext, ghost);

  // create moment arrays
  struct gkyl_array *m0, *m1i, *m2;
  m0 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  m1i = mkarr(vdim*confBasis.num_basis, confLocal_ext.volume);
  m2 = mkarr(confBasis.num_basis, confLocal_ext.volume);

  gkyl_proj_on_basis *proj_m0 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_M0, NULL);
  gkyl_proj_on_basis *proj_m1i = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, vdim, eval_M1i_3v, NULL);
  gkyl_proj_on_basis *proj_m2 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_M2_3v, NULL);

  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0);
  gkyl_proj_on_basis_advance(proj_m1i, 0.0, &confLocal, m1i);
  gkyl_proj_on_basis_advance(proj_m2, 0.0, &confLocal, m2);

  // M0, M1, M2 well defined here

  // create distribution function array
  struct gkyl_array *distf;
  distf = mkarr(basis.num_basis, local_ext.volume);

  // projection updater to compute MJ
  gkyl_proj_MJ_on_basis *proj_MJ = gkyl_proj_MJ_on_basis_new(&grid,
    &confBasis, &basis, poly_order+1);

  gkyl_proj_MJ_on_basis_fluid_stationary_frame_mom(proj_MJ, &local, &confLocal, m0, m1i, m2, distf);

// build the p_over_gamma
struct gkyl_array *p_over_gamma;
p_over_gamma = mkarr(vdim*velBasis.num_basis, velLocal.volume);
gkyl_proj_on_basis *p_over_gamma_proj = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
    .grid = &vel_grid,
    .basis = &velBasis,
    .qtype = GKYL_GAUSS_LOBATTO_QUAD,
    .num_quad = 8,
    .num_ret_vals = vdim,
    .eval = p_over_gamma_func[vdim-1],
    .ctx = 0
  });
  gkyl_proj_on_basis_advance(p_over_gamma_proj, 0.0, &velLocal, p_over_gamma);

  // correct the MJ distribution m0 Moment
  gkyl_correct_MJ *corr_MJ = gkyl_correct_MJ_new(&grid,&confBasis,&basis,&confLocal,&velLocal,confLocal.volume,confLocal_ext.volume, false);
  gkyl_correct_MJ_fix(corr_MJ,p_over_gamma,distf,m0,m1i,&local,&confLocal);

  // values to compare  at index (1, 9, 9, 9) [remember, lower-left index is (1,1,1,1)]
  double p2_vals[] = { 1.0 };

  const double *fv = gkyl_array_cfetch(distf, gkyl_range_idx(&local_ext, (int[4]) { 1, 9, 9, 9 }));

  if (poly_order == 2) {
    for (int i=0; i<basis.num_basis; ++i){
      printf("fv[%d] = %1.16e\n",i,fv[i]);
      //TEST_CHECK( gkyl_compare_double(p2_vals[i], fv[i], 1e-12) );
    }
  }

  // Temporarily print out the n, vb, T: Recall the quantities are *1/sqrt(2)
  // due to them being the coefficients
struct gkyl_range_iter biter;
gkyl_range_iter_init(&biter, &confLocal);
while (gkyl_range_iter_next(&biter)) {
    long midx = gkyl_range_idx(&confLocal, biter.idx);

    // Update the moments
    double *m0_d = gkyl_array_fetch(m0, midx);
    double *m1i_d = gkyl_array_fetch(m1i, midx);
    double *m2_d = gkyl_array_fetch(m2, midx);

    printf("\n----------- Ouptuts Start (Check M0, M1i, M2: ctest_proj_MJ_on_basis.c) ---------\n");
    int i;
    for (i = 0; i < 3; ++i){
      printf("num[%d]: %1.16e\n",i,m0_d[i]);
    }
    for (i = 0; i < 9; ++i){ // 9 for 3d * p2 (3*3=9 numbers)
      printf("vb[%d] : %1.16e\n",i,m1i_d[i]);
    }
    for (i = 0; i < 3; ++i){
      printf("T[%d]: %1.16g\n",i,m2_d[i]);
    }
    printf("\n----------- Ouptuts End (Check M0, M1i, M2: ctest_proj_MJ_on_basis.c) ---------\n");
}


  // write distribution function to file
  char fname[1024];
  sprintf(fname, "ctest_proj_MJ_on_basis_test_1x3v_p%d.gkyl", poly_order);
  gkyl_grid_sub_array_write(&grid, &local, distf, fname);

  // release memory for moment data object
  gkyl_proj_on_basis_release(p_over_gamma_proj);
  gkyl_correct_MJ_release(corr_MJ);
  gkyl_array_release(m0); gkyl_array_release(m1i); gkyl_array_release(m2);
  gkyl_array_release(distf);
  gkyl_proj_MJ_on_basis_release(proj_MJ);
  gkyl_proj_on_basis_release(proj_m0);
  gkyl_proj_on_basis_release(proj_m1i);
  gkyl_proj_on_basis_release(proj_m2);
  gkyl_array_release(p_over_gamma);
}

void test_1x3v_p2() { test_1x3v(2); }


TEST_LIST = {
  { "test_1x1v_no_drift_p1", test_1x1v_no_drift_p1 },
  { "test_1x1v_no_drift_p2", test_1x1v_no_drift_p2 },

  { "test_1x1v_p2", test_1x1v_p2 },

  { "test_1x2v_p2", test_1x2v_p2 },

  { "test_1x3v_p2", test_1x3v_p2 },

  { NULL, NULL },
};
