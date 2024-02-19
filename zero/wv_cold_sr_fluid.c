#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_wv_cold_sr_fluid.h>

#define NUX 1
#define NUY 2
#define NUZ 3

#define SQ(x) ((x)*(x))

struct wv_cold_sr_fluid {
  struct gkyl_wv_eqn eqn; // base object
  double clight;
};

// relativistic case
struct cold_sr {
  double rhol, rhor;
  double vl[3], vr[3];
};

static inline void
cold_sr_fluid_flux(const double q[4], double *flux)
{
  const double c = 299792458.0; 
  double Vx = q[NUX]/sqrt(q[0]*q[0] + (q[NUX]*q[NUX] + q[NUY]*q[NUY] + q[NUZ]*q[NUZ])/(c*c)); 
  flux[0] =  q[0]*Vx; // N*Vx
  flux[NUX] =  q[NUX]*Vx; // N*Ux*Vx
  flux[NUY] =  q[NUY]*Vx; // N*Uy*Vx
  flux[NUZ] =  q[NUZ]*Vx; // N*Uz*Vx
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

  const double c = 299792458.0;
  double *wv = 0;

  // isolate left and right states
  double rhor = qr[0];
  if (qr[0] < 0.0) {
    rhor = 1.0;
  }
  double urx = qr[1]/(rhor);
  double ury = qr[2]/(rhor);
  double urz = qr[3]/(rhor);
  double rhol = ql[0];
  // Fix density if negative:
  if (ql[0] < 0.0) {
    rhol = 1.0;
  }
  double ulx = ql[1]/(rhol);
  double uly = ql[2]/(rhol);
  double ulz = ql[3]/(rhol);

  // compute the constants:
  double gammal = sqrt(1.0 + (ulx*ulx + uly*uly + ulz*ulz)/(c*c));
  double gammar = sqrt(1.0 + (urx*urx + ury*ury + urz*urz)/(c*c));
  double vlx = ulx/gammal;
  double vrx = urx/gammar;
  double vly = uly/gammal;
  double vry = ury/gammar;
  double vlz = ulz/gammal;
  double vrz = urz/gammar;

  // Primative rho
  double rhol_prim = rhol/(gammal);
  double rhor_prim = rhor/(gammar);

  double k = (sqrt(rhol_prim) + sqrt(rhor_prim))/(2.0*c);
  double w0 = (sqrt(rhol_prim)*gammal + sqrt(rhor_prim)*gammar)/(2.0);
  double w1 = (sqrt(rhol_prim)*gammal*vlx/c + sqrt(rhor_prim)*gammar*vrx/c)/(2.0);
  double w2 = (sqrt(rhol_prim)*gammal*vly/c + sqrt(rhor_prim)*gammar*vry/c)/(2.0);
  double w3 = (sqrt(rhol_prim)*gammal*vlz/c + sqrt(rhor_prim)*gammar*vrz/c)/(2.0);

  // TEMP: for prim waves, use_conserved_var = 0;
  // Decide between primative or conserved jumps
  bool use_conserved_var = 1;

  if (use_conserved_var){
    // Assign the jump in the state vector, d = (d0,d1,d2,d3) 
    double d0 = qr[0] - ql[0]; 
    double d1 = qr[1] - ql[1]; 
    double d2 = qr[2] - ql[2]; 
    double d3 = qr[3] - ql[3]; 

    // Wave 1: eigenvalue is (c*w1)/w0;
    wv = &waves[0];
    wv[0] =  (c * (d1 * c * c * k * k * k - d0 * c * c * k * k * w1 + d1 * k * w0 * w0 - d1 * k * w1 * w1 - 2 * d2 * k * w1 * w2 - 2 * d3 * k * w1 * w3 + d1 * k * w2 * w2 + d1 * k * w3 * w3 - d0 * w0 * w0 * w1 + d0 * w1 * w1 * w1 + d0 * w1 * w2 * w2 + d0 * w1 * w3 * w3)) / (w0 * (c * c * k * k - w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3)); 
    wv[1] =  (2 * c * w1 * (d1 * c * c * k * k - d0 * w1 * c * c * k + d1 * w2 * w2 - d2 * w1 * w2 + d1 * w3 * w3 - d3 * w1 * w3)) / (w0 * (c * c * k * k - w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3));
    wv[2] =  (c * (d2 * c * c * k * k * w1 + d1 * c * c * k * k * w2 - 2 * d0 * c * c * k * w1 * w2 - d2 * w0 * w0 * w1 + d1 * w0 * w0 * w2 + d2 * w1 * w1 * w1 - d1 * w1 * w1 * w2 - d2 * w1 * w2 * w2 - 2 * d3 * w1 * w2 * w3 + d2 * w1 * w3 * w3 + d1 * w2 * w2 * w2 + d1 * w2 * w3 * w3)) / (w0 * (c * c * k * k - w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3));
    wv[3] =  (c * (d3 * c * c * k * k * w1 + d1 * c * c * k * k * w3 - 2 * d0 * c * c * k * w1 * w3 - d3 * w0 * w0 * w1 + d1 * w0 * w0 * w3 + d3 * w1 * w1 * w1 - d1 * w1 * w1 * w3 + d3 * w1 * w2 * w2 - 2 * d2 * w1 * w2 * w3 - d3 * w1 * w3 * w3 + d1 * w2 * w2 * w3 + d1 * w3 * w3 * w3)) / (w0 * (c * c * k * k - w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3));
    for (int i=0; i<4; ++i) if (isnan(wv[i]))  wv[i] = 0;
    for (int i=0; i<4; ++i) if (!(isfinite(wv[i])))  wv[i] = 0;
    s[0] = (c*w1)/w0;
    if (isnan(s[0])) s[0] = 0;
    if (!(isfinite(s[0]))) s[0] = 0;

    // Wave 2: eigenvalue is 2*w0*w1/( k*k + w0*w0 + w1*w1 + w2*w2 + w3*w3 );
    wv = &waves[4];
    wv[0] =  -(2 * c *  k * w0 * (d1 * c * c * k * k - 2 * d0 * c * c * k * w1 + d1 * w0 * w0 - d1 * w1 * w1 - 2 * d2 * w1 * w2 - 2 * d3 * w1 * w3 + d1 * w2 * w2 + d1 * w3 * w3)) / ((c * c * k * k - w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3) * (c * c * k * k + w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3));
    wv[1] =  -(2 * c * w0 * w1 * (d1 * c * c * k * k - 2 * d0 * c * c * k * w1 + d1 * w0 * w0 - d1 * w1 * w1 - 2 * d2 * w1 * w2 - 2 * d3 * w1 * w3 + d1 * w2 * w2 + d1 * w3 * w3)) / ((c * c * k * k - w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3) * (c * c * k * k + w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3));
    wv[2] =  -(2 * c * w0 * w2 * (d1 * c * c * k * k - 2 * d0 * c * c * k * w1 + d1 * w0 * w0 - d1 * w1 * w1 - 2 * d2 * w1 * w2 - 2 * d3 * w1 * w3 + d1 * w2 * w2 + d1 * w3 * w3)) / ((c * c * k * k - w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3) * (c * c * k * k + w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3));
    wv[3] =  -(2 * c * w0 * w3 * (d1 * c * c * k * k - 2 * d0 * c * c * k * w1 + d1 * w0 * w0 - d1 * w1 * w1 - 2 * d2 * w1 * w2 - 2 * d3 * w1 * w3 + d1 * w2 * w2 + d1 * w3 * w3)) / ((c * c * k * k - w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3) * (c * c * k * k + w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3));
    for (int i=0; i<4; ++i) if (isnan(wv[i])) wv[i] = 0;
    for (int i=0; i<4; ++i) if (!(isfinite(wv[i]))) wv[i] = 0;
    s[1] =  2.0*c*w0*w1/( c*c*k*k + w0*w0 + w1*w1 + w2*w2 + w3*w3 ); 
    if (isnan(s[1])) s[1] = 0;
    if (!(isfinite(s[1]))) s[1] = 0;


  // Use Primative variables
  } else {

    // Alternate primitive jump formulation, requires 3 waves instead of 2!

    // Assign the jump in the PRIMATIVE state vector, d = (d0,d1,d2,d3) 
    double d0 = rhor - rhol; 
    double d1 = urx - ulx; 
    double d2 = ury - uly; 
    double d3 = urz - ulz; 

    // Compute some useful terms
    double kstar = sqrt(rhol_prim)*sqrt(rhor_prim)/(c*c);
    double a_sqrt = sqrt((c * c * c * c * c * c * k * k * k * k * kstar * kstar) +
                    (2 * c * c * c * c * k * k * kstar * kstar * w0 * w0) +
                    (c * c * kstar * kstar * w0 * w0 * w0 * w0) -
                    (2 * c * k * kstar * w0 * w1 * w1) -
                    (2 * c * k * kstar * w0 * w2 * w2) -
                    (2 * c * k * kstar * w0 * w3 * w3) +
                    (k * k * w0 * w0));
    double sigma = c*c*k*k + w0*w0 + w1*w1 + w2*w2 + w3*w3; 
    double epsilon = -c*c*k*k - w0*w0 + w1*w1 + w2*w2 + w3*w3; 

    //Norm of fwave1
    double c1_wv1 = c * (kstar * c * c * c * k * k + kstar * c * w0 * w0 + k * w0 + a_sqrt);
    double c2_wv1 = (kstar*c*c*c*k*k + kstar*c*w0*w0 - k*w0 + a_sqrt);
    double c3_wv1 = (2*c2_wv1*d1*w1*w1 - c2_wv1*d1*sigma + 2*c2_wv1*d2*w1*w2 + 2*c2_wv1*d3*w1*w3 + 2*d0*epsilon*w0*w1);
    double nf1 = (c1_wv1*c3_wv1/(2*a_sqrt*epsilon*k*sigma));

    //Norm of fwave2
    double c1_wv2 = kstar * c * c * c * k * k + kstar * c * w0 * w0 + k * w0 - a_sqrt;
    double c2_wv2 = -kstar * c * c * c * k * k - kstar * c * w0 * w0 + k * w0 + a_sqrt;
    double c3_wv2 = 2*c2_wv2*d1*w1*w1 - c2_wv2*d1*sigma + 2*c2_wv2*d2*w1*w2 + 2*c2_wv2*d3*w1*w3 - 2*d0*epsilon*w0*w1;
    double nf2 = c*c1_wv2*c3_wv2/(2*a_sqrt*epsilon*k*sigma);

    // Move from eigenvalues to eigenvectors 
    double c_norm = w0*c*kstar/k; 

    // Set the waves
    wv = &waves[0];
    wv[0] = (c2_wv2*nf1)/(2*w0);
    wv[1] = nf1*w1;
    wv[2] = nf1*w2;
    wv[3] = nf1*w3;
    s[0] = (c*w1*(kstar*c*c*c*k*k + kstar*c*w0*w0 + k*w0 + a_sqrt))/(k*sigma);

    wv = &waves[4];
    wv[0] = -(c2_wv1*nf2)/(2*w0);
    wv[1] = nf2*w1;
    wv[2] = nf2*w2;
    wv[3] = nf2*w3;
    s[1] = (c*w1*(kstar*c*c*c*k*k + kstar*c*w0*w0 + k*w0 - a_sqrt))/(k*sigma);

    wv = &waves[8];
    wv[0] = 0;
    wv[1] = (2 * c * c * kstar * w1 * (d1 * w2 * w2 - d2 * w1 * w2 + d1 * w3 * w3 - d3 * w1 * w3)) / (epsilon * k);
    wv[2] = -(c * c * kstar * (d2 * w0 * w0 * w1 - d2 * w1 * w1 * w1 + d2 * w1 * w2 * w2 - d2 * w1 * w3 * w3 - d1 * w2 * (-2 * w1 * w1 + sigma) + 2 * d3 * w1 * w2 * w3 + c * c * d2 * k * k * w1)) / (epsilon * k);
    wv[3] = -(c * c * kstar * (d3 * w0 * w0 * w1 - d3 * w1 * w1 * w1 - d3 * w1 * w2 * w2 + d3 * w1 * w3 * w3 - d1 * w3 * (-2 * w1 * w1 + sigma) + 2 * d2 * w1 * w2 * w3 + c * c * d3 * k * k * w1)) / (epsilon * k);
    s[2] = (c*c*kstar*w1)/k;

    // Iterate over waves, correcting for overflow/underflow
    for (int i=0; i<3; ++i){
      s[i] = s[i]/c_norm; // Very important, converts lambda->V, by moving const to F
      if (isnan(s[i]))  s[i] = 0;
      if (!(isfinite(s[i])))  s[i] = 0;
    }
    for (int i=0; i<12; ++i){
      waves[i] = waves[i]/c_norm; // Very important, converts lambda->V, by moving const to F
      if (isnan(waves[i])) waves[i] = 0;
      if (!(isfinite(waves[i])))  waves[i] = 0;
    }

  }


  return fmax(fabs(s[0]), fabs(s[1]));
}

static void
qfluct_roe(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *ql, const double *qr, const double *waves, const double *s,
  double *amdq, double *apdq)
{
  //Q-Waves will not work with this system (L,R Eigenvectors are not unique)
}

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

  const double c = 299792458.0;
  double amaxl =  ql[NUX]/sqrt(ql[0]*ql[0] + (ql[NUX]*ql[NUX] + ql[NUY]*ql[NUY] + ql[NUZ]*ql[NUZ])/(c*c));
  double amaxr =  qr[NUX]/sqrt(qr[0]*qr[0] + (qr[NUX]*qr[NUX] + qr[NUY]*qr[NUY] + qr[NUZ]*qr[NUZ])/(c*c)); 

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
