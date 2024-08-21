#include <math.h>
#include <stdbool.h>

#include <gkyl_alloc.h>
#include <gkyl_wv_advect.h>

struct wv_advect {
  struct gkyl_wv_eqn eqn; // base object
  double c;               // advection speed
};

static void advect_free(const struct gkyl_ref_count *ref) {
  struct gkyl_wv_eqn *base = container_of(ref, struct gkyl_wv_eqn, ref_count);
  struct wv_advect *advect = container_of(base, struct wv_advect, eqn);
  gkyl_free(advect);
}

static inline void
cons_to_riem(const struct gkyl_wv_eqn *eqn,
  const double *qstate, const double *qin, double *wout)
{
  wout[0] = qin[0];
}
static inline void
riem_to_cons(const struct gkyl_wv_eqn *eqn,
  const double *qstate, const double *win, double *qout)
{
  qout[0] = win[0];
}

static inline void
rot_to_local(const double *tau1, const double *tau2, const double *norm,
  const double *GKYL_RESTRICT qglobal, double *GKYL_RESTRICT qlocal)
{
  qlocal[0] = qglobal[0];
}

static inline void
rot_to_global(const double *tau1, const double *tau2, const double *norm,
  const double *GKYL_RESTRICT qlocal, double *GKYL_RESTRICT qglobal)
{
  qglobal[0] = qlocal[0];
}

static double
wave(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *delta, const double *ql, const double *qr, 
  double *waves, double *s)
{
  const struct wv_advect *advect = container_of(eqn, struct wv_advect, eqn);

  double c = advect->c;

  waves[0] = delta[0];
  s[0] = c;

  return fabs(c);
}

static void
qfluct(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *ql, const double *qr, const double *waves, const double *s,
  double *amdq, double *apdq)
{
  const struct wv_advect *advect = container_of(eqn, struct wv_advect, eqn);

  double c = advect->c;

  if (c > 0) {
    apdq[0] = c * (qr[0]-ql[0]);
    amdq[0] = 0.0;
  } else {
    apdq[0] = 0.0;
    amdq[0] = c * (qr[0]-ql[0]);
  }
}

static double
flux_jump(const struct gkyl_wv_eqn *eqn, const double *ql, const double *qr, double *flux_jump)
{
  const struct wv_advect *advect = container_of(eqn, struct wv_advect, eqn);

  double c = advect->c;

  flux_jump[0] = c * (qr[0] - ql[0]);

  return fabs(c);
}

static bool
check_inv(const struct gkyl_wv_eqn *eqn, const double *q)
{
  return true; // no negative states in advection
}

static double
max_speed(const struct gkyl_wv_eqn *eqn, const double *q)
{
  const struct wv_advect *advect = container_of(eqn, struct wv_advect, eqn);
  return fabs(advect->c);
}

static inline void
advect_source(const struct gkyl_wv_eqn* eqn, const double* qin, double* sout)
{
  sout[0] = 0.0;
}

struct gkyl_wv_eqn *gkyl_wv_advect_new(double c) {
  struct wv_advect *advect = gkyl_malloc(sizeof(struct wv_advect));

  advect->eqn.type = GKYL_EQN_MAXWELL;
  advect->eqn.num_equations = 1;
  advect->eqn.num_waves = 1;
  advect->eqn.num_diag = 1;

  advect->c = c;

  advect->eqn.waves_func = wave;
  advect->eqn.qfluct_func = qfluct;

  advect->eqn.flux_jump = flux_jump;
  advect->eqn.check_inv_func = check_inv;
  advect->eqn.max_speed_func = max_speed;
  advect->eqn.rotate_to_local_func = rot_to_local;
  advect->eqn.rotate_to_global_func = rot_to_global;

  advect->eqn.cons_to_riem = cons_to_riem;
  advect->eqn.riem_to_cons = riem_to_cons;

  advect->eqn.cons_to_diag = gkyl_default_cons_to_diag;

  advect->eqn.source_func = advect_source;

  advect->eqn.ref_count = gkyl_ref_count_init(advect_free);

  return &advect->eqn;
}
