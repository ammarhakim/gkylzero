#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_wv_cold_sr_fluid.h>

#define NUX 1
#define NUY 2
#define NUZ 3

#define SQ(x) ((x)*(x))

struct wv_cold_sr_fluid {
  struct gkyl_wv_eqn eqn; // base object
};

void
cold_sr_fluid_flux(const double q[5], double flux[5])
{
  // Vx = NUx/(N^2 + NU^2/c^2) 
  const double c = 299792458.0;
  double Vx = q[NUX]/sqrt(q[0]*q[0] + (q[NUX]*q[NUX] + q[NUY]*q[NUY] + q[NUZ]*q[NUZ])/(c*c)); 
  flux[0] = q[0]*Vx; // N*Vx
  flux[NUX] = q[NUX]*Vx; // N*Ux*Vx
  flux[NUY] = q[NUY]*Vx; // N*Uy*Vx
  flux[NUZ] = q[NUZ]*Vx; // N*Uz*Vx
}

static inline void
cold_sr_fluid_cons_to_diag(const struct gkyl_wv_eqn *eqn,
  const double *qin, double *diag)
{
  // density and moment is copied as-is
  for (int i=0; i<4; ++i) diag[i] = qin[i];
  // Ke-sr density (gamma-1)
  const double c = 299792458.0;
  double ke = sqrt(qin[0]*qin[0] + (qin[1]*qin[1] + qin[2]*qin[2] + qin[3]*qin[3])/(c*c))/qin[0] - 1.0; 
  diag[4] = ke;
}

static void
cold_sr_fluid_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_wv_eqn *base = container_of(ref, struct gkyl_wv_eqn, ref_count);
  struct wv_cold_sr_fluid *cold_sr_fluid = container_of(base, struct wv_cold_sr_fluid, eqn);
  gkyl_free(cold_sr_fluid);
}

static inline void
cons_to_riem_sr(const struct gkyl_wv_eqn *eqn,
  const double *qstate, const double *qin, double *wout)
{
  // TODO: this should use proper L matrix
  for (int i=0; i<4; ++i)
    wout[i] = qin[i];
}
static inline void
riem_to_cons_sr(const struct gkyl_wv_eqn *eqn,
  const double *qstate, const double *win, double *qout)
{
  // TODO: this should use proper L matrix
  for (int i=0; i<4; ++i)
    qout[i] = win[i];
}

static inline void
rot_to_local(const double *tau1, const double *tau2, const double *norm,
  const double *GKYL_RESTRICT qglobal, double *GKYL_RESTRICT qlocal)
{
  qlocal[0] = qglobal[0];
  qlocal[1] = qglobal[1]*norm[0] + qglobal[2]*norm[1] + qglobal[3]*norm[2];
  qlocal[2] = qglobal[1]*tau1[0] + qglobal[2]*tau1[1] + qglobal[3]*tau1[2];
  qlocal[3] = qglobal[1]*tau2[0] + qglobal[2]*tau2[1] + qglobal[3]*tau2[2];
}

static inline void
rot_to_global(const double *tau1, const double *tau2, const double *norm,
  const double *GKYL_RESTRICT qlocal, double *GKYL_RESTRICT qglobal)
{
  qglobal[0] = qlocal[0];
  qglobal[1] = qlocal[1]*norm[0] + qlocal[2]*tau1[0] + qlocal[3]*tau2[0];
  qglobal[2] = qlocal[1]*norm[1] + qlocal[2]*tau1[1] + qlocal[3]*tau2[1];
  qglobal[3] = qlocal[1]*norm[2] + qlocal[2]*tau1[2] + qlocal[3]*tau2[2];
}


// Waves and speeds using Roe averaging
// corrective terms
static double
wave_roe_sr(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *delta, const double *ql, const double *qr, double *waves, double *s)
{
  double f[4];
  const double c = 299792458.0;
  double vl = ql[NUX]/sqrt(ql[0]*ql[0] + (ql[NUX]*ql[NUX] + ql[NUY]*ql[NUY] + ql[NUZ]*ql[NUZ])/(c*c));
  double vr = qr[NUX]/sqrt(qr[0]*qr[0] + (qr[NUX]*qr[NUX] + qr[NUY]*qr[NUY] + qr[NUZ]*qr[NUZ])/(c*c)); 

  double *wv = 0;

  // corrective term
  if ((vl < 0.0) && (0.0 < vr)) { // vacuum intermediate state will be formed
    cold_sr_fluid_flux(ql, f);
    wv = &waves[0];
    for(int m=0; m<4; ++m) wv[m] = -f[m];
    s[0] = vl;

    cold_sr_fluid_flux(qr, f);
    wv = &waves[4];
    for(int m=0; m<4; ++m) wv[m] = f[m];
    s[1] = vr;
  }
  else {
    // no vacuum state
    double rl = ql[0];
    double rr = qr[0];
    double u_diag[3];
    // compute Roe averaged speed
    double vel = compute_sr_roe_averaged_velocity_via_ridders(ql,qr,c,u_diag);
    //printf("vL: %1.16f, vs %1.16f, vR: %1.16f\n",vl,vel,vr);
            
    if(vel<0) {
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
        wv[m] = 0.0;

      wv = &waves[4];
      for(int m=0; m<4; ++m)
        wv[m] = delta[m];
    }
    s[0] = vel;
    s[1] = vel;
  }

  return fmax(fabs(s[0]), fabs(s[1]));
}

static void
qfluct_roe(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
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

// First order solve
static void
ffluct_roe(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
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
flux_jump_sr(const struct gkyl_wv_eqn *eqn, const double *ql, const double *qr, double *flux_jump_sr)
{
  double fr[4], fl[4];
  cold_sr_fluid_flux(ql, fl);
  cold_sr_fluid_flux(qr, fr);

  for (int m=0; m<4; ++m) flux_jump_sr[m] = fr[m]-fl[m];

 // Vn = NUn/(N^2 + NU^2/c^2)
  const double c = 299792458.0;
  double amaxl = ql[NUX]/sqrt(ql[0]*ql[0] + (ql[NUX]*ql[NUX] + ql[NUY]*ql[NUY] + ql[NUZ]*ql[NUZ])/(c*c));
  double amaxr = qr[NUX]/sqrt(qr[0]*qr[0] + (qr[NUX]*qr[NUX] + qr[NUY]*qr[NUY] + qr[NUZ]*qr[NUZ])/(c*c)); 

  return fmax(amaxl, amaxr);
}

static bool
check_inv(const struct gkyl_wv_eqn *eqn, const double *q)
{
  return q[0] > 0.0;
}

static double
max_speed_sr(const struct gkyl_wv_eqn *eqn, const double *q)
{
  const struct wv_cold_sr_fluid *cold_sr_fluid = container_of(eqn, struct wv_cold_sr_fluid, eqn);
  const double c = 299792458.0;
  return fabs(q[NUX]/sqrt(q[0]*q[0] + (q[NUX]*q[NUX] + q[NUY]*q[NUY] + q[NUZ]*q[NUZ])/(c*c)));
}

struct gkyl_wv_eqn*
gkyl_wv_cold_sr_fluid_new(void)
{
  struct wv_cold_sr_fluid *cold_sr_fluid = gkyl_malloc(sizeof(struct wv_cold_sr_fluid));

  cold_sr_fluid->eqn.type = GKYL_EQN_COLDFLUID_SR;
  cold_sr_fluid->eqn.num_equations = 4;
  cold_sr_fluid->eqn.num_waves = 2;
  cold_sr_fluid->eqn.num_diag = 5; // KE is final component
  
  cold_sr_fluid->eqn.waves_func = wave_roe_sr;
  cold_sr_fluid->eqn.qfluct_func = qfluct_roe;
  cold_sr_fluid->eqn.ffluct_func = ffluct_roe;
  cold_sr_fluid->eqn.flux_jump = flux_jump_sr;
  
  cold_sr_fluid->eqn.check_inv_func = check_inv;
  cold_sr_fluid->eqn.max_speed_func = max_speed_sr;

  cold_sr_fluid->eqn.rotate_to_local_func = rot_to_local;
  cold_sr_fluid->eqn.rotate_to_global_func = rot_to_global;

  cold_sr_fluid->eqn.cons_to_riem = cons_to_riem_sr;
  cold_sr_fluid->eqn.riem_to_cons = riem_to_cons_sr;

  cold_sr_fluid->eqn.cons_to_diag = cold_sr_fluid_cons_to_diag;

  cold_sr_fluid->eqn.ref_count = gkyl_ref_count_init(cold_sr_fluid_free);

  return &cold_sr_fluid->eqn;
}
