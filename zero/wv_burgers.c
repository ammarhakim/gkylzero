#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_wv_burgers.h>

struct wv_burgers {
  struct gkyl_wv_eqn eqn; // base object
};

static void
burgers_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_wv_eqn *base = container_of(ref, struct gkyl_wv_eqn, ref_count);
  struct wv_burgers *burgers = container_of(base, struct wv_burgers, eqn);
  gkyl_free(burgers);
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

// Waves and speeds using Roe averaging
static double
wave_roe(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *delta, const double *ql, const double *qr, 
  double *waves, double *s)
{
  const struct wv_burgers *burgers = container_of(eqn, struct wv_burgers, eqn);

  waves[0] = delta[0];
  s[0] = 0.5*(ql[0]+qr[0]);
  
  return s[0];
}

static void
qfluct_roe(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *ql, const double *qr, const double *waves, const double *s,
  double *amdq, double *apdq)
{
  amdq[0] = apdq[0] = 0.0;
  if (s[0] < 0)
    amdq[0] = s[0]*waves[0];
  else
    apdq[0] = s[0]*waves[0];

  if ((ql[0]<0.0) && (0.0<qr[0])) {
    amdq[0] = -0.5*ql[0]*ql[0];
    apdq[0] = 0.5*qr[0]*qr[0];
  }  
}

static void
ffluct_roe(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *ql, const double *qr, const double *waves, const double *s,
  double *amdq, double *apdq)
{
  amdq[0] = apdq[0] = 0.0;

  if (fabs(s[0]) < 1e-15) {
    amdq[0] = apdq[0] = 0.5*waves[0];
  }
  else {
    if (s[0] < 0)
      amdq[0] = waves[0];
    else
      apdq[0] = waves[0];
  }
  if ((ql[0]<0.0) && (0.0<qr[0])) {
    amdq[0] = 0.0 -0.5*ql[0]*ql[0];
    apdq[0] = 0.5*qr[0]*qr[0] - 0.0;
  }
}

static double
flux_jump(const struct gkyl_wv_eqn *eqn, const double *ql, const double *qr, double *flux_jump)
{
  flux_jump[0] = 0.5*(qr[0]*qr[0]-ql[0]*ql[0]);
  return fmax(fabs(ql[0]), fabs(qr[0]));
}

static bool
check_inv(const struct gkyl_wv_eqn *eqn, const double *q)
{
  return true; // no negative states in Burgers
}

static double
max_speed(const struct gkyl_wv_eqn *eqn, const double *q)
{
  const struct wv_burgers *burgers = container_of(eqn, struct wv_burgers, eqn);
  return fabs(q[0]);
}

static inline void
burgers_source(const struct gkyl_wv_eqn* eqn, const double* qin, double* sout)
{
  sout[0] = 0.0;
}

struct gkyl_wv_eqn*
gkyl_wv_burgers_new(void)
{
  struct wv_burgers *burgers = gkyl_malloc(sizeof(struct wv_burgers));

  burgers->eqn.type = GKYL_EQN_BURGERS;
  burgers->eqn.num_equations = 1;
  burgers->eqn.num_waves = 1;
  burgers->eqn.num_diag = 1;
  
  burgers->eqn.waves_func = wave_roe;
  burgers->eqn.qfluct_func = qfluct_roe;
  burgers->eqn.ffluct_func = ffluct_roe;

  burgers->eqn.flux_jump = flux_jump;
  burgers->eqn.check_inv_func = check_inv;
  burgers->eqn.max_speed_func = max_speed;
  burgers->eqn.rotate_to_local_func = rot_to_local;
  burgers->eqn.rotate_to_global_func = rot_to_global;

  burgers->eqn.cons_to_riem = cons_to_riem;
  burgers->eqn.riem_to_cons = riem_to_cons;

  burgers->eqn.cons_to_diag = gkyl_default_cons_to_diag;

  burgers->eqn.source_func = burgers_source;

  burgers->eqn.ref_count = gkyl_ref_count_init(burgers_free);

  return &burgers->eqn;
}
