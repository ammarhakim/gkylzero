#include <float.h>
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

static inline double
fmax3(double a, double b, double c)
{
  return fmax(fmax(a,b),c);
}

static inline double
fmin3(double a, double b, double c)
{
  return fmin(fmin(a,b),c);
}

// limiter function
static inline double
limiter_function(double r, enum gkyl_wave_limiter limiter)
{
  double theta = 0.0;
  switch (limiter) {
    case GKYL_NO_LIMITER:
        theta = 1.0;
        break;
    
    case GKYL_MIN_MOD:
        theta = fmax(0, fmin(1, r));
        break;

    case GKYL_SUPERBEE:
        theta = fmax3(0.0, fmin(1, 2*r), fmin(2.0, r));
        break;

    case GKYL_VAN_LEER:
        theta = (r+fabs(r))/(1+fabs(r));
        break;

    case GKYL_MONOTONIZED_CENTERED:
        theta = fmax(0.0, fmin3((1.0+r)/2, 2, 2*r));
        break;

    case GKYL_BEAM_WARMING:
        theta = r;
        break;

    case GKYL_ZERO:
        theta = 0;
        break;
  }
  return theta;
}

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
calc_first_order_update(int meqn, double dtdx,
  double *restrict ql, double *restrict qr, const double *amdq, const double *apdq)
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

static inline double
wave_dot_prod(int meqn, const double *restrict wa, const double *restrict wb)
{
  double dot = 0.0;
  for (int i=0; i<meqn; ++i) dot += wa[i]*wb[i];
  return dot;
}

static inline void
wave_rescale(int meqn, double fact, double *w)
{
  for (int i=0; i<meqn; ++i) w[i] *= fact; 
}

static inline void
calc_second_order_flux(int meqn, double dtdx, double s,
  const double *waves, double *restrict flux2)
{
  double sfact = 0.5*fabs(s)*(1-fabs(s)*dtdx);
  for (int i=0; i<meqn; ++i)
    flux2[i] += sfact*waves[i];
}

static inline void
calc_second_order_update(int meqn, double dtdx, double *restrict qout,
  const double *fl, const double *fr)
{
  for (int i=0; i<meqn; ++i)
    qout[i] += -dtdx*(fr[i]-fl[i]);
}

static void
limit_waves(const gkyl_wave_prop *wv, const struct gkyl_range *slice_range,
  int lower, int upper, struct gkyl_array *waves, const struct gkyl_array *speed)
{
  int meqn = wv->equation->num_equations, mwaves = wv->equation->num_waves;

  for (int mw=0; mw<mwaves; ++mw) {
    const double *wl = gkyl_array_cfetch(waves, gkyl_ridx(*slice_range, lower-1));
    const double *wr = gkyl_array_cfetch(waves, gkyl_ridx(*slice_range, lower));

    double dotr = wave_dot_prod(meqn, &wl[mw*meqn], &wr[mw*meqn]);

    for (int i=lower; i<=upper; ++i) {
      double dotl = dotr;
      
      double *restrict wi = gkyl_array_fetch(waves, gkyl_ridx(*slice_range, i));
      const double *restrict wi1 = gkyl_array_cfetch(waves, gkyl_ridx(*slice_range, i+1));
      
      double wnorm2 = wave_dot_prod(meqn, &wi[mw*meqn], &wi[mw*meqn]);
      dotr = wave_dot_prod(meqn, &wi[mw*meqn], &wi1[mw*meqn]);

      if (wnorm2 > 0) {
        const double *s = gkyl_array_cfetch(speed, gkyl_ridx(*slice_range, i));
        double r = s[mw] > 0 ? dotl/wnorm2 : dotr/wnorm2;
        double theta = limiter_function(r, wv->limiter);
        wave_rescale(meqn, theta, &wi[mw*meqn]);
      }
    }
  }
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
  double delta[meqn], amdq[meqn], apdq[meqn];

  int idxl[GKYL_MAX_DIM], idxr[GKYL_MAX_DIM];

  gkyl_array_copy(qout, qin);

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

        long sidx = gkyl_ridx(slice_range, i);
        long lidx = gkyl_range_idx(update_range, idxl);
        long ridx = gkyl_range_idx(update_range, idxr);        

        const double *qinl = gkyl_array_cfetch(qin, lidx);
        const double *qinr = gkyl_array_cfetch(qin, ridx);

        double *waves = gkyl_array_fetch(wv->waves, sidx);
        double *s = gkyl_array_fetch(wv->speeds, sidx);

        calc_jump(meqn, qinl, qinr, delta);

        gkyl_wv_eqn_waves(wv->equation, dir, delta, qinl, qinr, waves, s);
        gkyl_wv_eqn_qfluct(wv->equation, dir, qinl, qinr, waves, s, amdq, apdq);

        double *qoutl = gkyl_array_fetch(qout, lidx);
        double *qoutr = gkyl_array_fetch(qout, ridx);

        calc_first_order_update(meqn, dtdx, qoutl, qoutr, amdq, apdq);
        cfla = calc_cfla(mwaves, cfla, dtdx, s);
      }

      if (cfla > cflm)
        return (struct gkyl_wave_prop_status) { .success = 0, .dt_suggested = dt*cfl/cfla };

      // apply limiters to waves for all edges in update range,
      // including edges that are on the range boundary
      limit_waves(wv, &slice_range,
        update_range->lower[dir], update_range->upper[dir]+1, wv->waves, wv->speeds);

      gkyl_array_clear(wv->flux2, 0.0);
      // compute second-order correction fluxes at each interface:
      // note that there is one extra edge than cell
      for (int i=update_range->lower[dir]; i<=update_range->upper[dir]+1; ++i) {
        long sidx = gkyl_ridx(slice_range, i);

        const double *waves = gkyl_array_cfetch(wv->waves, sidx);
        const double *s = gkyl_array_cfetch(wv->speeds, sidx);
        double *flux2 = gkyl_array_fetch(wv->flux2, sidx);

        for (int mw=0; mw<mwaves; ++mw)
          calc_second_order_flux(meqn, dtdx, s[mw], &waves[mw*meqn], flux2);
      }

      // add second correction flux to solution in each interior cell
      for (int i=update_range->lower[dir]; i<=update_range->upper[dir]; ++i) {
        idxl[dir] = i;
        long sidx = gkyl_ridx(slice_range, i);
        
        calc_second_order_update(meqn, dtdx,
          gkyl_array_fetch(qout, gkyl_range_idx(update_range, idxl)),
          gkyl_array_cfetch(wv->flux2, gkyl_ridx(slice_range, i)),
          gkyl_array_cfetch(wv->flux2, gkyl_ridx(slice_range, i+1))
        );
      }
    }
  }

  // compute allowable time-step from this update, but suggest only
  // bigger time-step; (Only way dt can reduce is if the update
  // fails. If the code comes here the update suceeded and so we
  // should not allow dt to reduce).
  double dt_suggested = dt*cfl/cfla;

  return (struct gkyl_wave_prop_status) {
    .success = 1,
    .dt_suggested = dt_suggested > dt ? dt_suggested : dt
  };
}

double
gkyl_wave_prop_max_dt(const gkyl_wave_prop *wv, const struct gkyl_range *update_range,
  const struct gkyl_array *qin)
{
  double max_dt = DBL_MAX;
  
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_range);
  while (gkyl_range_iter_next(&iter)) {

    for (int d=0; d<wv->num_up_dirs; ++d) {
      int dir = wv->update_dirs[d];
      double dx = wv->grid.dx[dir];

      const double *q = gkyl_array_cfetch(qin, gkyl_range_idx(update_range, iter.idx));
      double maxs = gkyl_wv_eqn_max_speed(wv->equation, dir, q);
      max_dt = fmin(max_dt, wv->cfl*dx/maxs);
    }
    
  }

  return max_dt;
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
