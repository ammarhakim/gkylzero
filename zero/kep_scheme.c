#include <float.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_kep_scheme.h>
#include <gkyl_moment_prim_euler.h>
#include <gkyl_moment_braginskii_priv.h>
#include <gkyl_wv_euler.h>
#include <gkyl_util.h>

static const int dir_shuffle[][3] = {
  {1, 2, 3},
  {2, 3, 1},
  {3, 1, 2}
};

// for indexing non-ideal terms, which are located at cell nodes
static void
create_offsets_centers(const struct gkyl_range *range, long offsets[])
{
  // box spanning stencil
  struct gkyl_range box3;
  gkyl_range_init(&box3, range->ndim, (int[]) { 0, 0, 0 }, (int[]) { 1, 1, 1 });

  struct gkyl_range_iter iter3;
  gkyl_range_iter_init(&iter3, &box3);

  // construct list of offsets
  int count = 0;
  while (gkyl_range_iter_next(&iter3))
    offsets[count++] = gkyl_range_offset(range, iter3.idx);
}

struct gkyl_kep_scheme {
  struct gkyl_rect_grid grid; // grid object
  int ndim; // number of dimensions
  int num_up_dirs; // number of update directions
  int update_dirs[GKYL_MAX_DIM]; // directions to update
  const struct gkyl_wv_eqn *equation; // equation object
  double cfl; // cfl number
};

gkyl_kep_scheme*
gkyl_kep_scheme_new(struct gkyl_kep_scheme_inp inp)
{
  struct gkyl_kep_scheme *up;
  up = gkyl_malloc(sizeof(*up));

  up->grid = *(inp.grid);
  up->ndim = up->grid.ndim;
  
  up->num_up_dirs = inp.num_up_dirs;
  for (int i=0; i<inp.num_up_dirs; ++i)
    up->update_dirs[i] = inp.update_dirs[i];

  up->equation = gkyl_wv_eqn_acquire(inp.equation);

  up->cfl = inp.cfl;
  return up;
}

// compute kinetic energy (no density factor)
static inline double euler_ke(const double v[5])
{
  return 0.5*(v[1]*v[1] + v[2]*v[2] + v[3]*v[3]);
}

// Numerical flux using modified KEP scheme
static void
mkep_flux(int dir, double gas_gamma, const double vm[5], const double vp[5], double flux[5])
{       
#define RHOU d[0]
#define RHOV d[1]
#define RHOW d[2]
#define PR 4

  double vbar[5];
  
  for (int i=0; i<5; ++i) vbar[i] = 0.5*(vm[i]+vp[i]);
  double kebar = 0.5*(euler_ke(vm) + euler_ke(vp));
  
  const int *d = dir_shuffle[dir];
  flux[0] = vbar[RHO]*vbar[RHOU]; // rho*u
  // momentum flux must have this form to ensure KEP property
  flux[RHOU] =  flux[0]*vbar[RHOU] + vbar[PR]; // rho*u*u + pe
  flux[RHOV] =  flux[0]*vbar[RHOV]; // rho*u*v
  flux[RHOW] =  flux[0]*vbar[RHOW]; // rho*u*v
  // following ensure stability of linear perturbations around uniform flow
  flux[ER] = gas_gamma/(gas_gamma-1)*vbar[PR]*vbar[RHOU] + flux[0]*kebar; // (E+p)*u
 
#undef RHOU
#undef RHOV
#undef RHOW
#undef PR  
}

static void
non_ideal_calc_update(const gkyl_kep_scheme *kep,
  const double *brag_d[], double *rhs)
{
  int ndim = kep->ndim;
  double div_pi[3] = {0.0};
  double div_q = {0.0};
  if (ndim == 1) {
    const double dx = kep->grid.dx[0];
    double pi[2][6] = {};
    double q[2][3] = {};

    for (int j = L_1D; j <= U_1D; ++j) {
      pi[j][PIXX] = brag_d[j][PIXX];
      pi[j][PIXY] = brag_d[j][PIXY];
      pi[j][PIXZ] = brag_d[j][PIXZ];
      pi[j][PIYY] = brag_d[j][PIYY];
      pi[j][PIYZ] = brag_d[j][PIYZ];
      pi[j][PIZZ] = brag_d[j][PIZZ];
    }

    div_pi[0] = calc_sym_grad_1D(dx, pi[L_1D][0], pi[U_1D][0]);
    div_pi[1] = calc_sym_grad_1D(dx, pi[L_1D][1], pi[U_1D][1]);
    div_pi[2] = calc_sym_grad_1D(dx, pi[L_1D][2], pi[U_1D][2]);

    rhs[RHO] += 0.0;
    rhs[MX] += -div_pi[0];
    rhs[MY] += -div_pi[1];
    rhs[MZ] += -div_pi[2];

    // Increment heat flux and viscous heating to energy variable
    for (int j = L_1D; j <= U_1D; ++j) {
      q[j][0] = brag_d[j][QX];
      q[j][1] = brag_d[j][QY];
      q[j][2] = brag_d[j][QZ];
    }

    div_q = calc_sym_grad_1D(dx, q[L_1D][0], q[U_1D][0]);
    rhs[ER] += -div_q;
  }
  else if (ndim == 2) {
    const double dx = kep->grid.dx[0];
    const double dy = kep->grid.dx[1];
    double pi[4][6] = {};
    double q[4][3] = {};

    for (int j = LL_2D; j <= UU_2D; ++j) {
      pi[j][PIXX] = brag_d[j][PIXX];
      pi[j][PIXY] = brag_d[j][PIXY];
      pi[j][PIXZ] = brag_d[j][PIXZ];
      pi[j][PIYY] = brag_d[j][PIYY];
      pi[j][PIYZ] = brag_d[j][PIYZ];
      pi[j][PIZZ] = brag_d[j][PIZZ];
    }

    div_pi[0] = calc_sym_gradx_2D(dx, pi[LL_2D][0], pi[LU_2D][0], pi[UL_2D][0], pi[UU_2D][0])
              + calc_sym_grady_2D(dy, pi[LL_2D][1], pi[LU_2D][1], pi[UL_2D][1], pi[UU_2D][1]);
    
    div_pi[1] = calc_sym_gradx_2D(dx, pi[LL_2D][1], pi[LU_2D][1], pi[UL_2D][1], pi[UU_2D][1])
              + calc_sym_grady_2D(dy, pi[LL_2D][3], pi[LU_2D][3], pi[UL_2D][3], pi[UU_2D][3]);
    
    div_pi[2] = calc_sym_gradx_2D(dx, pi[LL_2D][2], pi[LU_2D][2], pi[UL_2D][2], pi[UU_2D][2])
              + calc_sym_grady_2D(dy, pi[LL_2D][4], pi[LU_2D][4], pi[UL_2D][4], pi[UU_2D][4]);

    rhs[RHO] += 0.0;
    rhs[MX] += -div_pi[0];
    rhs[MY] += -div_pi[1];
    rhs[MZ] += -div_pi[2];
    // Increment heat flux and viscous heating to energy variable
    for (int j = LL_2D; j <= UU_2D; ++j){
      q[j][0] = brag_d[j][QX];
      q[j][1] = brag_d[j][QY];
      q[j][2] = brag_d[j][QZ];
    }

    div_q = calc_sym_gradx_2D(dx, q[LL_2D][0], q[LU_2D][0], q[UL_2D][0], q[UU_2D][0])
          + calc_sym_grady_2D(dy, q[LL_2D][1], q[LU_2D][1], q[UL_2D][1], q[UU_2D][1]);
    rhs[ER] += -div_q;
  }
}

struct gkyl_wave_prop_status
gkyl_kep_scheme_advance(const gkyl_kep_scheme *kep, double dt, 
  const struct gkyl_range *update_range, const struct gkyl_range *non_ideal_range,
  const struct gkyl_array *qin, const struct gkyl_array *non_ideal_vars, 
  struct gkyl_array *cflrate, struct gkyl_array *rhs)
{
  int ndim = update_range->ndim;
  long sz[] = { 2, 4, 8 };
  long offsets_centers[sz[ndim-1]];
  create_offsets_centers(update_range, offsets_centers);

  double gas_gamma = gkyl_wv_euler_gas_gamma(kep->equation);

  double cfla = 0.0, cfl = kep->cfl, cflm = 1.1*cfl;

  int idxl[GKYL_MAX_DIM], idxc[GKYL_MAX_DIM], idxr[GKYL_MAX_DIM];

  double vl[5], vc[5], vr[5], fluxl[5], fluxr[5];
  // array of pointers to non-ideal terms
  const double* non_ideal_vars_d[sz[ndim-1]];

  gkyl_array_clear_range(rhs, 0.0, *update_range);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_range);
  while (gkyl_range_iter_next(&iter)) {
    gkyl_copy_int_arr(ndim, iter.idx, idxc);
    long linc = gkyl_range_idx(update_range, idxc);

    const double *qc = gkyl_array_cfetch(qin, linc);
    gkyl_euler_prim_vars(gas_gamma, qc, vc);

    double *cflrate_d = gkyl_array_fetch(cflrate, linc);
    double *out = gkyl_array_fetch(rhs, linc);

    for (int d=0; d<kep->num_up_dirs; ++d) {
      int dir = kep->update_dirs[d];
      double dx = kep->grid.dx[dir];

      gkyl_copy_int_arr(ndim, iter.idx, idxl);
      gkyl_copy_int_arr(ndim, iter.idx, idxr);
      idxl[dir] = idxl[dir]-1; idxr[dir] = idxr[dir]+1;

      long linl = gkyl_range_idx(update_range, idxl); 
      long linr = gkyl_range_idx(update_range, idxr);

      const double *ql = gkyl_array_cfetch(qin, linl);
      const double *qr = gkyl_array_cfetch(qin, linr);

      gkyl_euler_prim_vars(gas_gamma, ql, vl);
      gkyl_euler_prim_vars(gas_gamma, qr, vr);

      mkep_flux(dir, gas_gamma, vl, vc, fluxl);
      mkep_flux(dir, gas_gamma, vc, vr, fluxr);

      out[0] += -(fluxr[0]-fluxl[0])/dx;
      out[1] += -(fluxr[1]-fluxl[1])/dx;
      out[2] += -(fluxr[2]-fluxl[2])/dx;
      out[3] += -(fluxr[3]-fluxl[3])/dx;
      out[4] += -(fluxr[4]-fluxl[4])/dx;

      // TODO: Rotation Is Required!
      cflrate_d[0] += gkyl_euler_max_abs_speed(gas_gamma, qc)/dx;
    }
    // Accumulate non-ideal terms
    long linc_vertex = gkyl_range_idx(non_ideal_range, iter.idx);
    
    for (int i=0; i<sz[ndim-1]; ++i)
      non_ideal_vars_d[i] = gkyl_array_cfetch(non_ideal_vars, linc_vertex + offsets_centers[i]);
    
    non_ideal_calc_update(kep, non_ideal_vars_d, out);
  }
  // compute allowable time-step from this update, but suggest only
  // bigger time-step; (Only way dt can reduce is if the update
  // fails. If the code comes here the update suceeded and so we
  // should not allow dt to reduce).
  gkyl_array_reduce_range(&cfla, cflrate, GKYL_MAX, *update_range);
  double dt_suggested = dt*cfl/fmax(cfla, DBL_MIN);

  return (struct gkyl_wave_prop_status) {
    .success = 1,
    .dt_suggested = dt_suggested > dt ? dt_suggested : dt
  };
}

void
gkyl_kep_scheme_release(gkyl_kep_scheme* up)
{
  gkyl_wv_eqn_release(up->equation);
  gkyl_free(up);
}
