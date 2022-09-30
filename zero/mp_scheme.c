#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_mp_scheme.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <gkyl_wave_geom.h>

#include <float.h>
#include <math.h>

// type signature for function to do recovery
typedef void (*recovery_fn_t)(int meqn,
  const double *f3m, const double *f2m, const double *fm,
  const double *fp, const double *f2p, const double *f3p,
  double *outl, double *outr);

struct gkyl_mp_scheme {
  struct gkyl_rect_grid grid; // grid object
  int ndim; // number of dimensions
  int num_up_dirs; // number of update directions
  int update_dirs[GKYL_MAX_DIM]; // directions to update
  enum gkyl_mp_recon mp_recon; // base reconstruction to use

  double cfl; // CFL number  
  
  const struct gkyl_wv_eqn *equation; // equation object
  struct gkyl_wave_geom *geom; // geometry object

  recovery_fn_t recovery_fn; // function to do recovery
};

// Each of the recovery methods below take 6 cells, three to the left
// and three to the right, and return value at the left/right of the
// interface. Note that depending on the scheme, some of these values
// may be ignored.

static inline void
c2_recovery(int meqn,
  const double *f3m, const double *f2m, const double *fm,
  const double *fp, const double *f2p, const double *f3p,
  double *outl, double *outr)
{
  // c2 is symmetric 2nd order scheme, so outl and outr are same
  for (int m=0; m<meqn; ++m)
    outr[m] = outl[m] = fm[m]/2.0 + fp[m]/2.0;
}

static inline void
c4_recovery(int meqn,
  const double *f3m, const double *f2m, const double *fm,
  const double *fp, const double *f2p, const double *f3p,
  double *outl, double *outr)
{
  // c4 is symmetric 4th order scheme, so outl and outr are same
  for (int m=0; m<meqn; ++m)
    outr[m] = outl[m] = -f2m[m]/12.0 + 7.0*fm[m]/12.0 + 7.0*fp[m]/12.0 - f2p[m]/12.0;
}

static inline void
c6_recovery(int meqn,
  const double *f3m, const double *f2m, const double *fm,
  const double *fp, const double *f2p, const double *f3p,
  double *outl, double *outr)
{
  // c6 is symmetric 6th order scheme, so outl and outr are same
  for (int m=0; m<meqn; ++m)  
    outr[m] = outl[m] =
      37.0*fp[m]/60.0+37.0*fm[m]/60.0+f3p[m]/60.0+f3m[m]/60.0-2.0*f2p[m]/15.0-2.0*f2m[m]/15.0;
}

static inline void
u1_recovery(int meqn,
  const double *f3m, const double *f2m, const double *fm,
  const double *fp, const double *f2p, const double *f3p,
  double *outl, double *outr)
{
  // u1 is upwind-biased 1st order scheme
  for (int m=0; m<meqn; ++m) {
    outl[m] = fm[m];
    outr[m] = fp[m];
  }
}

static inline void
u3_recovery(int meqn,
  const double *f3m, const double *f2m, const double *fm,
  const double *fp, const double *f2p, const double *f3p,
  double *outl, double *outr)
{
  // u3 is upwind-biased 3rd order scheme
  for (int m=0; m<meqn; ++m) {
    outl[m] = -1.0/6.0*f2m[m] + 5.0/6.0*fm[m] + 1.0/3.0*fp[m];
    outr[m] = 1.0/3.0*fm[m] + 5.0/6.0*fp[m] - 1.0/6.0*f2p[m];
  }
}

static inline void
u5_recovery(int meqn,
  const double *f3m, const double *f2m, const double *fm,
  const double *fp, const double *f2p, const double *f3p,
  double *outl, double *outr)
{
  // u5 is upwind-biased 5th order scheme
  for (int m=0; m<meqn; ++m) {
    outl[m] = 1.0/30.0*f3m[m] - 13.0/60.0*f2m[m] + 47.0/60.0*fm[m] + 9.0/20.0*fp[m] - 1.0/20.0*f2p[m];
    outr[m] = -1.0/20.0*f2m[m] + 9.0/20.0*fm[m] + 47.0/60.0*fp[m] - 13.0/60.0*f2p[m] + 1.0/30.0*f3p[m];
  }
}

gkyl_mp_scheme*
gkyl_mp_scheme_new(const struct gkyl_mp_scheme_inp *mpinp)
{
  struct gkyl_mp_scheme *mp = gkyl_malloc(sizeof(*mp));

  mp->grid = *(mpinp->grid);
  mp->ndim = mp->grid.ndim;
  
  mp->num_up_dirs = mpinp->num_up_dirs;
  for (int i=0; i<mpinp->num_up_dirs; ++i)
    mp->update_dirs[i] = mpinp->update_dirs[i];

  mp->mp_recon = mpinp->mp_recon;
  mp->cfl = mpinp->cfl;
  
  mp->equation = gkyl_wv_eqn_acquire(mpinp->equation);
  mp->geom = gkyl_wave_geom_acquire(mpinp->geom);

  switch (mpinp->mp_recon) {
    case GKYL_MP_C2:
      mp->recovery_fn = c2_recovery;
      break;
    case GKYL_MP_C4:
      mp->recovery_fn = c4_recovery;
      break;
    case GKYL_MP_C6:
      mp->recovery_fn = c6_recovery;
      break;
    case GKYL_MP_U1:
      mp->recovery_fn = u1_recovery;
      break;
    case GKYL_MP_U3:
      mp->recovery_fn = u3_recovery;
      break;
    case GKYL_MP_U5:
      mp->recovery_fn = u5_recovery;
      break;
  }

  return mp;
}

static inline long
get_offset(int dir, int loc, const struct gkyl_range *range)
{
  int idx[GKYL_MAX_CDIM] = { 0, 0, 0 };
  idx[dir] = loc;
  return gkyl_range_offset(range, idx);
}

void
gkyl_mp_scheme_advance(gkyl_mp_scheme *mp,
  const struct gkyl_range *update_range, const struct gkyl_array *qin,
  struct gkyl_array *qrec_l, struct gkyl_array *qrec_r,
  struct gkyl_array *amdq, struct gkyl_array *apdq,
  struct gkyl_array *cflrate, struct gkyl_array *rhs)
{
  int ndim = update_range->ndim;
  int meqn = mp->equation->num_equations;

  // labels for three cells to left, three cells to right of edge:
  enum { I3M, I2M, IM, IP, I2P, I3P };

  gkyl_array_clear_range(rhs, 0.0, *update_range);
  
  // outer loop is over direction: the RHS is updated direction
  // by direction
  for (int d=0; d<mp->num_up_dirs; ++d) {
    int dir = mp->update_dirs[d];

    double dx = mp->grid.dx[dir];

    // compute index offsets of cells on left/right of edge
    long offsets[6];
    offsets[IP]  = get_offset(dir, 0, update_range);
    offsets[I2P] = get_offset(dir, 1, update_range);
    offsets[I3P] = get_offset(dir, 2, update_range);
    offsets[IM]  = get_offset(dir, -1, update_range);
    offsets[I2M] = get_offset(dir, -2, update_range);
    offsets[I3M] = get_offset(dir, -3, update_range);

    const double *qrec_in[6]; // pointers to cells attached to edge

    // create range that includes one extra layer on the upper size
    int upper[GKYL_MAX_CDIM] = { 0 };
    for (int d=0; d<update_range->ndim; ++d) upper[d] = update_range->upper[d];
    upper[dir] += 1;
    struct gkyl_range update_range_ext;
    gkyl_range_init(&update_range_ext, update_range->ndim, update_range->lower, upper);

    struct gkyl_range_iter iter;
    // loop over all cells and recover left/right edge values in
    // direction 'dir'. Note this is effectively a loop over the
    // edges, and includes upper most edge also
    gkyl_range_iter_init(&iter, &update_range_ext);
    while (gkyl_range_iter_next(&iter)) {
      // Note: Edge is between cells IM and IP
      
      long loc = gkyl_range_idx(update_range, iter.idx);

      // attach pointers to cells for recovery
      for (int i=0; i<6; ++i)
        qrec_in[i] = gkyl_array_cfetch(qin, loc+offsets[i]);

      // qr_l is left of edge (right edge of left cell), qr_r right of
      // edge (left edge of right cell)
      double *qr_l = gkyl_array_fetch(qrec_r, loc+offsets[IM]);
      double *qr_r = gkyl_array_fetch(qrec_l, loc+offsets[IP]);

      // recover variables
      mp->recovery_fn(meqn, qrec_in[I3M], qrec_in[I2M], qrec_in[IM],
        qrec_in[IP], qrec_in[I2P], qrec_in[I3P],
        qr_l, qr_r);

      // TODO: Apply MP limiter

      double *amdq_p = gkyl_array_fetch(amdq, loc+offsets[IM]);
      double *apdq_p = gkyl_array_fetch(apdq, loc+offsets[IP]);
      // compute fluctuations: note the equation system must not use
      // either the waves or speeds
      gkyl_wv_eqn_qfluct(mp->equation, GKYL_WV_HIGH_ORDER_FLUX,
        qr_l, qr_r, 0, 0, amdq_p, apdq_p);
    }

    double deltaf[meqn];
    // Update RHS with contribution from flux jumps. Note this loop is
    // over interior cells
    gkyl_range_iter_init(&iter, update_range);
    while (gkyl_range_iter_next(&iter)) {
      long loc = gkyl_range_idx(update_range, iter.idx);

      double amax = gkyl_wv_eqn_flux_jump(mp->equation,
        gkyl_array_cfetch(qrec_l, loc), gkyl_array_cfetch(qrec_r, loc),
        deltaf);

      const double *amdq_p = gkyl_array_cfetch(amdq, loc);
      const double *apdq_p = gkyl_array_cfetch(apdq, loc);

      double *rhs_p = gkyl_array_fetch(rhs, loc);
      for (int m=0; m<meqn; ++m)
        rhs_p[m] += -deltaf[m]/dx - (apdq_p[m]+amdq_p[m])/dx;

      double *cflrate_p = gkyl_array_fetch(cflrate, loc);
      cflrate_p[0] += amax/dx;
    }
  }
}

double
gkyl_mp_scheme_max_dt(const gkyl_mp_scheme *mp, const struct gkyl_range *update_range,
  const struct gkyl_array *qin)
{
  double max_dt = DBL_MAX;
  
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_range);
  while (gkyl_range_iter_next(&iter)) {

    for (int d=0; d<mp->num_up_dirs; ++d) {
      int dir = mp->update_dirs[d];
      double dx = mp->grid.dx[dir];

      const double *q = gkyl_array_cfetch(qin, gkyl_range_idx(update_range, iter.idx));
      double maxs = gkyl_wv_eqn_max_speed(mp->equation, q);
      max_dt = fmin(max_dt, mp->cfl*dx/maxs);
    }
    
  }

  return max_dt;  
}

void
gkyl_mp_scheme_release(gkyl_mp_scheme* mp)
{
  gkyl_wv_eqn_release(mp->equation);
  gkyl_wave_geom_release(mp->geom);
  
  gkyl_free(mp);
}
