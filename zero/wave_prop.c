#include <math.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_util.h>
#include <gkyl_wave_prop.h>

struct gkyl_wave_prop {
    struct gkyl_rect_grid grid; // grid object
    int ndim; // number of dimensions
    int num_up_dirs; // number of update directions
    int update_dirs[GKYL_MAX_DIM]; // directions to update
    enum gkyl_wave_limiter limiter; // limiter to use
    double cfl; // CFL number
    const struct gkyl_wv_eqn *equation; // equation object

    // data for 1D slice update
    struct gkyl_array *waves, *speeds, *flux2;
};

gkyl_wave_prop*
gkyl_wave_prop_new(struct gkyl_wave_prop_inp winp)
{
  gkyl_wave_prop *up = gkyl_malloc(sizeof(gkyl_wave_prop));

  up->grid = *(winp.grid);
  up->ndim = up->grid.ndim;
  
  up->num_up_dirs = winp.num_up_dirs;
  for (int i=0; i<winp.num_up_dirs; ++i)
    up->update_dirs[i] = winp.update_dirs[i];

  up->limiter = winp.limiter;
  up->cfl = winp.cfl;
  up->equation = gkyl_wv_eqn_aquire(winp.equation);

  int nghost[3] = { 2, 2, 2 };
  struct gkyl_range range, ext_range;
  gkyl_create_grid_ranges(&up->grid, nghost, &ext_range, &range);

  int max_1d = 0;
  for (int d=0; d<ext_range.ndim; ++d) {
    int shape = gkyl_range_shape(&ext_range, d);
    max_1d = max_1d > shape ? max_1d : shape;
  }

  // allocate memory to store 1D slices of waves, speeds and
  // second-order correction flux
  int meqn = winp.equation->num_equations, mwaves = winp.equation->num_waves;
  up->waves = gkyl_array_new(GKYL_DOUBLE, meqn*mwaves, max_1d);
  up->speeds = gkyl_array_new(GKYL_DOUBLE, mwaves, max_1d);
  up->flux2 = gkyl_array_new(GKYL_DOUBLE, meqn, max_1d);

  return up;
}

// some helper functions
static inline void
calc_jump(int n, const double *ql, const double *qr, double *restrict jump)
{
  for (int d=0; d<n; ++d) jump[d] = qr[d]-ql[d];
}

static inline void
calc_first_order_update(int meqn, double dtdx, double *restrict ql, double *restrict qr,
  const double *amdq, const double *apdq)
{
  for (int i=0; i<meqn; ++i) {
    qr[i] = qr[i] - dtdx*apdq[i];
    ql[i] = ql[i] - dtdx*amdq[i];
  }
}

static inline double
calc_cfla(int mwaves, double cfla, double dtdx, const double *s)
{
  double c = cfla;
  for (int i=0; i<mwaves; ++i)
    c = fmax(c, dtdx*fabs(s[i]));
  return c;
}

// advance method
struct gkyl_wave_prop_status
gkyl_wave_prop_advance(const gkyl_wave_prop *wv,
  double tm, double dt, const struct gkyl_range *update_range,
  const struct gkyl_array *qin, struct gkyl_array *qout)
{
  int ndim = update_range->ndim;
  int meqn = wv->equation->num_equations, mwaves = wv->equation->num_waves;

  double cfla = 0.0, cfl = wv->cfl, cflm = 1.1*cfl;
  double delta[meqn], waves[meqn*mwaves], s[mwaves];
  double amdq[meqn], apdq[meqn];

  int idxl[GKYL_MAX_DIM], idxr[GKYL_MAX_DIM];
  
  for (int d=0; d<wv->num_up_dirs; ++d) {
    int dir = wv->update_dirs[d];

    double dtdx = dt/wv->grid.dx[dir];

    // upper/lower bounds in direction 'd'. Note these are edge indices
    int loidx = update_range->lower[dir]-1;
    int upidx = update_range->upper[dir]+2;

    struct gkyl_range slice_range;
    gkyl_range_init(&slice_range, 1, (int[]) { loidx }, (int[]) { upidx } );

    struct gkyl_range perp_range;
    gkyl_range_shorten(&perp_range, update_range, dir, 1);
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &perp_range);

    // outer loop is over perpendicular directions, inner loop over 1D
    // slice along that direction
    while (gkyl_range_iter_next(&iter)) {
      
      gkyl_copy_int_arr(ndim, iter.idx, idxl);
      gkyl_copy_int_arr(ndim, iter.idx, idxr);

      for (int i=loidx; i<upidx; ++i) {
        idxl[dir] = i-1; idxr[dir] = i;

        const double *qinl = gkyl_array_cfetch(qin, gkyl_range_idx(update_range, idxl));
        const double *qinr = gkyl_array_cfetch(qin, gkyl_range_idx(update_range, idxr));

        // compute jump across interface
        calc_jump(meqn, qinl, qinr, delta);
        // compute waves and fluctuations
        gkyl_wv_eqn_waves(wv->equation, dir, delta, qinl, qinr, waves, s);
        gkyl_wv_eqn_qfluct(wv->equation, dir, qinl, qinr, waves, s, amdq, apdq);

        double *qoutl = gkyl_array_fetch(qout, gkyl_range_idx(update_range, idxl));
        double *qoutr = gkyl_array_fetch(qout, gkyl_range_idx(update_range, idxr));

        calc_first_order_update(meqn, dtdx, qoutl, qoutr, amdq, apdq);
        cfla = calc_cfla(mwaves, cfla, dtdx, s);
      }

      if (cfla > cflm)
        return (struct gkyl_wave_prop_status) { .success = 0, .dt_suggested = dt*cfl/cfla };
    }
  }

  return (struct gkyl_wave_prop_status) { .success = 1, .dt_suggested = dt };
}

void
gkyl_wave_prop_release(gkyl_wave_prop* up)
{
  gkyl_wv_eqn_release(up->equation);
  gkyl_array_release(up->waves);
  gkyl_array_release(up->speeds);
  gkyl_array_release(up->flux2);
  
  free(up);
}
