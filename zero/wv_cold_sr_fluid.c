#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_wv_coldfluid.h>
#include <gkyl_wv_cold_sr_fluid.h>

#define NUX 1
#define NUY 2
#define NUZ 3

#define SQ(x) ((x)*(x))

struct wv_cold_sr_fluid {
  struct gkyl_wv_sr_eqn eqn; // base object
};

void
cold_sr_fluid_flux(const double q[5], double flux[5])
{
  // TODO: 1. Assumes normalized by c?, 2. Select either Vx, Vy, Vz?
  double Vx = q[NUX]/sqrt(q[0]^2 + (q[NUX]^2 + q[NUY]^2 + q[NUZ]^2)); // Vx = NUx/(N^2 + NU^2/c^2) 
  flux[0] = q[0]*Vx; // N*Vx
  flux[NUX] = q[NUX]*Vx; // N*Ux*Vx
  flux[NUY] = q[NUY]*Vx; // N*Uy*Vx
  flux[NUZ] = q[NUZ]*Vx; // N*Uz*Vx
}

static inline void
cold_sr_fluid_cons_to_diag(const struct gkyl_wv_sr_eqn *eqn,
  const double *qin, double *diag)
{
  // density and moment is copied as-is
  for (int i=0; i<4; ++i) diag[i] = qin[i];
  double ke = sqrt(qin[0]^2 + (qin[1]^2 + qin[2]^2 + qin[3]^2))/qin[0] - 1.0; // Ke-sr density (gamma-1)
  diag[4] = ke;
}

static void
cold_sr_fluid_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_wv_sr_eqn *base = container_of(ref, struct gkyl_wv_sr_eqn, ref_count);
  struct wv_cold_sr_fluid *cold_sr_fluid = container_of(base, struct wv_cold_sr_fluid, eqn);
  gkyl_free(cold_sr_fluid);
}

static inline void
cons_to_riem_sr(const struct gkyl_wv_sr_eqn *eqn,
  const double *qstate, const double *qin, double *wout)
{
  // TODO: this should use proper L matrix
  for (int i=0; i<4; ++i)
    wout[i] = qin[i];
}
static inline void
riem_to_cons_sr(const struct gkyl_wv_sr_eqn *eqn,
  const double *qstate, const double *win, double *qout)
{
  // TODO: this should use proper L matrix
  for (int i=0; i<4; ++i)
    qout[i] = win[i];
}

// Waves and speeds using Roe averaging, correction term
static double
wave_roe_sr(const struct gkyl_wv_sr_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *delta, const double *ql, const double *qr, double *waves, double *s)
{
  double f[4];
  double ur = qr[NUX]/qr[0], ul = ql[NUX]/ql[0];

  double *wv = 0;

  // corrective term
  if ((ul < 0) && (0 < ur)) { // vacuum intermediate state will be formed
    cold_sr_fluid_flux(ql, f);
    wv = &waves[0];
    for(int m=0; m<4; ++m) wv[m] = -f[m];
    s[0] = ul;

    cold_sr_fluid_flux(qr, f);
    wv = &waves[4];
    for(int m=0; m<4; ++m) wv[m] = f[m];
    s[1] = ur;
  }
  else {
    // no vacuum state
    double rl = ql[0];
    double rr = qr[0];
    // compute Roe averaged speed
    double Vx_hat = compute_sr_roe_avereged_velocity(ql,qr,uav);
            
    if(uav<0) {
      wv = &waves[0];
      for(int m=0; m<4; ++m)
        wv[m] = delta[m];

      wv = &waves[4];
      for(int m=0; m<4; ++m)
        wv[m] = 0.0;
    }
    else {
      wv = &waves[0];
      for(int m=0; m<4; ++m)
        wv[m] = 0;

      wv = &waves[4];
      for(int m=0; m<4; ++m)
        wv[m] = delta[m];
    }
    s[0] = uav;
    s[1] = uav;
  }

  return fmax(fabs(s[0]), fabs(s[1]));
}

static void
qfluct_roe(const struct gkyl_wv_sr_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *ql, const double *qr, const double *waves, const double *s,
  double *amdq, double *apdq)
{
  int meqn = 4, mwaves = 2;
  
  for (int m=0; m<meqn; ++m) {
    amdq[m] = 0.0; apdq[m] = 0.0;

    for (int mw=0; mw<mwaves; ++mw) {
      const double *wv = &waves[mw*meqn];
      
      if (s[mw] < 0.0)
        amdq[m] += s[mw]*wv[m];
      else
        apdq[m] += s[mw]*wv[m];
    }
  }
}

static void
ffluct_roe(const struct gkyl_wv_sr_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *ql, const double *qr, const double *waves, const double *s,
  double *amdq, double *apdq)
{
  int meqn = 4, mwaves = 2;
  
  for (int m=0; m<meqn; ++m) {
    amdq[m] = 0.0; apdq[m] = 0.0;

    for (int mw=0; mw<mwaves; ++mw) {
      const double *wv = &waves[mw*meqn];
      
      if (s[mw] < 0.0) {
        amdq[m] += wv[m];
      }
      else if (s[mw] > 0.0) {
        apdq[m] += wv[m];
      }
      else {
        amdq[m] += 0.5*wv[m];
        apdq[m] += 0.5*wv[m];
      }
    }
  }
}

static double
flux_jump(const struct gkyl_wv_sr_eqn *eqn, const double *ql, const double *qr, double *flux_jump)
{
  double fr[4], fl[4];
  cold_sr_fluid_flux(ql, fl);
  cold_sr_fluid_flux(qr, fr);

  for (int m=0; m<4; ++m) flux_jump[m] = fr[m]-fl[m];

  double amaxl = ql[NUX]/ql[0];
  double amaxr = qr[NUX]/qr[0];

  return fmax(amaxl, amaxr);
}

static bool
check_inv(const struct gkyl_wv_sr_eqn *eqn, const double *q)
{
  return q[0] > 0.0;
}

static double
max_speed(const struct gkyl_wv_sr_eqn *eqn, const double *q)
{
  const struct wv_cold_sr_fluid *cold_sr_fluid = container_of(eqn, struct wv_cold_sr_fluid, eqn);
  return fabs(q[NUX]/q[0]);
}

struct gkyl_wv_sr_eqn*
gkyl_wv_cold_sr_fluid_new(void)
{
  struct wv_cold_sr_fluid *cold_sr_fluid = gkyl_malloc(sizeof(struct wv_cold_sr_fluid));

  cold_sr_fluid->eqn.type = GKYL_EQN_cold_sr_flUID;
  cold_sr_fluid->eqn.num_equations = 4;
  cold_sr_fluid->eqn.num_waves = 2;
  cold_sr_fluid->eqn.num_diag = 5; // KE is final component
  
  cold_sr_fluid->eqn.waves_func = wave_roe_sr;
  cold_sr_fluid->eqn.qfluct_func = qfluct_roe;
  cold_sr_fluid->eqn.ffluct_func = ffluct_roe;
  cold_sr_fluid->eqn.flux_jump = flux_jump;
  
  cold_sr_fluid->eqn.check_inv_func = check_inv;
  cold_sr_fluid->eqn.max_speed_func = max_speed;

  cold_sr_fluid->eqn.rotate_to_local_func = rot_to_local;
  cold_sr_fluid->eqn.rotate_to_global_func = rot_to_global;

  cold_sr_fluid->eqn.cons_to_riem_sr = cons_to_riem_sr;
  cold_sr_fluid->eqn.riem_to_cons_sr = riem_to_cons_sr;

  cold_sr_fluid->eqn.cons_to_diag = cold_sr_fluid_cons_to_diag;

  cold_sr_fluid->eqn.ref_count = gkyl_ref_count_init(cold_sr_fluid_free);

  return &cold_sr_fluid->eqn;
}
