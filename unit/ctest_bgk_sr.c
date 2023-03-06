#include "gkyl_array.h"
#include "gkyl_util.h"
#include <acutest.h>

#include <gkyl_array_rio.h>
#include <gkyl_correct_maxwellian.h>
#include <gkyl_correct_mj.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_mj_moments.h>
#include <gkyl_proj_mj_on_basis.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <math.h>

#include <gkyl_bgk_collisions.h>
#include <gkyl_bgk_collisions_priv.h>
#include <gkyl_array_ops_priv.h>

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

// waterbag distribution
void evalFunc(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], vx = xn[1];
  fout[0] = 0.0;
  if (vx < 2.0 &&  vx > 0.0){
    fout[0] = 1.0;
  }
}

// Collision distribution all phase space
void evalFunc_nu(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], vx = xn[1];
  // calculate nu, assumed constant
  double dt = 0.01;
  double nu_all = 1.0;
  double nudt = nu_all*dt;
  fout[0] = nudt;
}

// Collision conf space
void evalFunc_nu_conf(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  // calculate nu, assumed constant
  double dt = 0.01;
  double nu_all = 1.0;
  double nudt = nu_all*dt;
  fout[0] = nudt;
}


void
test_1x1v(int poly_order)
{
  double lower[] = {0.1, -32.0}, upper[] = {1.0, 32.0};
  int cells[] = {2, 320};
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

  // create distribution function array
  struct gkyl_array *distf;
  distf = mkarr(basis.num_basis, local_ext.volume);
  struct gkyl_array *distf_mj;
  distf_mj = mkarr(basis.num_basis, local_ext.volume);
  struct gkyl_array *bgk_out;
  bgk_out = mkarr(basis.num_basis, local_ext.volume);
  struct gkyl_array *nudt;
  nudt = mkarr(basis.num_basis, local_ext.volume);
  struct gkyl_array *nudt_conf;
  nudt_conf = mkarr(confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *cflfreq;
  cflfreq = mkarr(basis.num_basis, local_ext.volume);

  // projection updater to compute mj
  gkyl_proj_mj_on_basis *proj_mj = gkyl_proj_mj_on_basis_new(&grid,
    &confBasis, &basis, poly_order+1);

  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc, NULL);

  gkyl_proj_on_basis *projnu = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc_nu, NULL);

  gkyl_proj_on_basis *projnu_conf = gkyl_proj_on_basis_new(&confGrid, &confBasis,
      poly_order+1, 1, evalFunc_nu_conf, NULL);

  // create the bgk operator structure
  gkyl_bgk_collisions *bgk_obj = gkyl_bgk_collisions_new(&confBasis, &basis, false);

  // create the waterbag distribution for distf
  gkyl_proj_on_basis_advance(projDistf, 0.0, &local, distf);

  // create the collision frequency matrix
  gkyl_proj_on_basis_advance(projnu, 0.0, &local, nudt);

  // create the collision frequency matrix in conf spade
  gkyl_proj_on_basis_advance(projnu_conf, 0.0, &confLocal, nudt_conf);

  // timeloop evolving partial_t(f) = -nu(f-f^mj)
  for (int i=0; i<1000; ++i){ //3000

    //printf("\n----------- ************************* ---------\n");
    //printf("----------- Begining timeloop: T = %d ---------\n", i);
    //printf("----------- ************************* ---------\n\n");

    // write distribution function to file
    char fname[1024];
    if ( i == 999 || i == 0){
      sprintf(fname, "ctest_bgk_sr_1x1v_p%d_time_%03d.gkyl", poly_order,i);
    }
    gkyl_grid_sub_array_write(&grid, &local, distf, fname);

    // calculate the moments of the dist (n, vb, T -> m0, m1i, m2)
    gkyl_mj_moments *mj_moms = gkyl_mj_moments_new(&grid,&confBasis,&basis,&confLocal,&velLocal,confLocal.volume,confLocal_ext.volume, false);
    gkyl_mj_moments_advance(mj_moms,p_over_gamma,gamma,gamma_inv,distf,m0,m1i,m2,&local,&confLocal);

    // Update the dist_mj using the moments
    gkyl_proj_mj_on_basis_fluid_stationary_frame_mom(proj_mj, &local, &confLocal, m0, m1i, m2, distf_mj);

    // correct the mj distribution m0 (n) Moment
    gkyl_correct_mj *corr_mj = gkyl_correct_mj_new(&grid,&confBasis,&basis,&confLocal,&velLocal,confLocal.volume,confLocal_ext.volume, false);
    gkyl_correct_mj_fix(corr_mj,p_over_gamma,distf_mj,m0,m1i,&local,&confLocal);

    // calculate nu*f^mj,
    //gkyl_dg_mul_op_range(basis, 0, distf_mj, 0, nudt, 0, distf_mj, &local);
    gkyl_dg_mul_conf_phase_op_range(&confBasis, &basis,
      distf_mj, nudt_conf, distf_mj, &confLocal, &local);


      // TEMPORARY: Verf the mj projections
      //gkyl_mj_moments_advance(mj_moms,p_over_gamma,gamma,gamma_inv,distf_mj,m0,m1i,m2,&local,&confLocal);


    // calculate the BGK contribution
    gkyl_bgk_collisions_advance(bgk_obj,
      &confLocal, &local,
      nudt_conf, distf_mj, distf,
      bgk_out, cflfreq);

    // update f with the BGK operation (Forward Euler update)
    // f^j+1 = f^j + out;  out = -dt*nu*(f - f^mj)
    struct gkyl_range_iter piter;
    gkyl_range_iter_init(&piter, &local);
    while (gkyl_range_iter_next(&piter)) {
      long ploc = gkyl_range_idx(&local, piter.idx);
      double *out_d = gkyl_array_fetch(bgk_out, ploc);
      double *fv = gkyl_array_fetch(distf, ploc);

      // Add out_d to the evolving distribution function
      array_acc1(basis.num_basis, fv, 1., out_d);

    }

    // Release the memory
    gkyl_correct_mj_release(corr_mj);
    gkyl_mj_moments_release(mj_moms);
  }
  // end timeloop



  // values to compare  at index (1, 17) [remember, lower-left index is (1,1)]
  double p2_vals[] = {  0.4908183182853421,-1.779993558391936e-17,
0.01957103464383442,4.791671361253045e-18, -2.283871108067826e-17,
-0.000481110805805566,-9.813947903343648e-19, -1.29480691392842e-17 };

  const double *fv = gkyl_array_cfetch(distf, gkyl_range_idx(&local_ext, (int[2]) { 1, 160 }));

  if (poly_order == 2) {
    for (int i=0; i<basis.num_basis; ++i){
      //printf("fv[%d] = %1.16g\n",i,fv[i]);
      TEST_CHECK( gkyl_compare_double(p2_vals[i], fv[i], 1e-12) );
    }
  }



  // release memory for moment data object
  gkyl_proj_on_basis_release(p_over_gamma_proj);
  gkyl_proj_on_basis_release(gamma_proj);
  gkyl_proj_on_basis_release(gamma_inv_proj);
  gkyl_array_release(m0); gkyl_array_release(m1i); gkyl_array_release(m2);
  gkyl_array_release(distf);
  gkyl_array_release(distf_mj);
  gkyl_array_release(bgk_out);
  gkyl_array_release(nudt);
  gkyl_array_release(nudt_conf);
  gkyl_array_release(cflfreq);
  gkyl_proj_mj_on_basis_release(proj_mj);
  gkyl_proj_on_basis_release(projDistf);
  gkyl_proj_on_basis_release(projnu);
  gkyl_proj_on_basis_release(projnu_conf);
  gkyl_proj_on_basis_release(proj_m0);
  gkyl_proj_on_basis_release(proj_m1i);
  gkyl_proj_on_basis_release(proj_m2);
  gkyl_array_release(p_over_gamma);
  gkyl_bgk_collisions_release(bgk_obj);
}

// special note, the p1 basis does not function
void test_1x1v_p2() { test_1x1v(2); }


TEST_LIST = {
  { "test_1x1v_p2", test_1x1v_p2 },
  { NULL, NULL },
};
