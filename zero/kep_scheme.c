#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_kep_scheme.h>
#include <gkyl_prim_euler.h>
#include <gkyl_wv_euler.h>

static const int dir_shuffle[][3] = {
  {1, 2, 3},
  {2, 3, 1},
  {3, 1, 2}
};

struct gkyl_kep_scheme {
  struct gkyl_rect_grid grid; // grid object
  int ndim; // number of dimensions
  int num_up_dirs; // number of update directions
  int update_dirs[GKYL_MAX_DIM]; // directions to update
  const struct gkyl_wv_eqn *equation; // equation object
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
#define RHO 0     
#define RHOU d[0]
#define RHOV d[1]
#define RHOW d[2]
#define ER 4
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

#undef RHO  
#undef RHOU
#undef RHOV
#undef RHOW
#undef ER
#undef PR  
}

void
gkyl_kep_scheme_advance(const gkyl_kep_scheme *kep, const struct gkyl_range *update_rng,
  const struct gkyl_array *qin, struct gkyl_array *cflrate, struct gkyl_array *rhs)
{
  int ndim = update_rng->ndim;
  double gas_gamma = gkyl_wv_euler_gas_gamma(kep->equation);

  int idxm[GKYL_MAX_DIM], idxp[GKYL_MAX_DIM];

  double vm[5], vp[5], flux[5];

  gkyl_array_clear_range(rhs, 0.0, *update_rng);

  for (int d=0; d<kep->num_up_dirs; ++d) {
    int dir = kep->update_dirs[d];
    double dx = kep->grid.dx[dir];

    int loidx = update_rng->lower[dir];
    int upidx = update_rng->upper[dir]+1; // one more edge than cells

    struct gkyl_range perp_range;
    gkyl_range_shorten(&perp_range, update_rng, dir, 1);
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &perp_range);

    while (gkyl_range_iter_next(&iter)) {

      gkyl_copy_int_arr(ndim, iter.idx, idxm);
      gkyl_copy_int_arr(ndim, iter.idx, idxp);
      
      for (int i=loidx; i<=upidx; ++i) { // note upidx is inclusive
        idxm[dir] = i-1; idxp[dir] = i;

        long linm = gkyl_range_idx(update_rng, idxm);
        long linp = gkyl_range_idx(update_rng, idxp);

        const double *qm = gkyl_array_cfetch(qin, linm);
        const double *qp = gkyl_array_cfetch(qin, linp);

        gkyl_euler_prim_vars(gas_gamma, qm, vm);
        gkyl_euler_prim_vars(gas_gamma, qp, vp);

        mkep_flux(dir, gas_gamma, vm, vp, flux);

        // accumulate contribution of flux to left/right cell
        double *rhsm = gkyl_array_fetch(rhs, linm);
        double *rhsp = gkyl_array_fetch(rhs, linp);
        for (int m=0; m<5; ++m) {
          rhsp[m] += flux[m]/dx;
          rhsm[m] += -flux[m]/dx;
        }

        double *cflrate_d = gkyl_array_fetch(cflrate, linp);
        // TODO: Rotation Is Required!
        cflrate_d[0] += gkyl_euler_max_abs_speed(gas_gamma, qp)/dx;
      }
    }
  }
}

void
gkyl_kep_scheme_release(gkyl_kep_scheme* up)
{
  gkyl_wv_eqn_release(up->equation);
  gkyl_free(up);
}
