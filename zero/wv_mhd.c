#include <math.h>
#include <float.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_moment_prim_mhd.h>
#include <gkyl_wv_mhd.h>

// Make indexing cleaner and clearer
#define DN (0)
#define MX (1)
#define MY (2)
#define MZ (3)
#define ER (4)
#define BX (5)
#define BY (6)
#define BZ (7)
#define PSI_GLM (8)

#define sq(x) ((x)*(x))

static inline void
rot_to_local_rect(const double *tau1, const double *tau2, const double *norm,
  const double *GKYL_RESTRICT qglobal, double *GKYL_RESTRICT qlocal)
{
  // Mass density is a scalar
  qlocal[0] = qglobal[0];
  // Rotate momentum to local coordinates
  qlocal[1] = qglobal[1]*norm[0] + qglobal[2]*norm[1] + qglobal[3]*norm[2];
  qlocal[2] = qglobal[1]*tau1[0] + qglobal[2]*tau1[1] + qglobal[3]*tau1[2];
  qlocal[3] = qglobal[1]*tau2[0] + qglobal[2]*tau2[1] + qglobal[3]*tau2[2];
  // Total energy is a scalar  
  qlocal[4] = qglobal[4];
  // Rotate B to local coordinates
  qlocal[5] = qglobal[5]*norm[0] + qglobal[6]*norm[1] + qglobal[7]*norm[2];
  qlocal[6] = qglobal[5]*tau1[0] + qglobal[6]*tau1[1] + qglobal[7]*tau1[2];
  qlocal[7] = qglobal[5]*tau2[0] + qglobal[6]*tau2[1] + qglobal[7]*tau2[2];
}

static inline void
rot_to_global_rect(const double *tau1, const double *tau2, const double *norm,
  const double *GKYL_RESTRICT qlocal, double *GKYL_RESTRICT qglobal)
{
  // Mass density is a scalar
  qglobal[0] = qlocal[0];
  // Rotate momentum back to global coordinates
  qglobal[1] = qlocal[1]*norm[0] + qlocal[2]*tau1[0] + qlocal[3]*tau2[0];
  qglobal[2] = qlocal[1]*norm[1] + qlocal[2]*tau1[1] + qlocal[3]*tau2[1];
  qglobal[3] = qlocal[1]*norm[2] + qlocal[2]*tau1[2] + qlocal[3]*tau2[2];
  // Total energy is a scalar  
  qglobal[4] = qlocal[4];
  // Rotate B back to global coordinates
  qglobal[5] = qlocal[5]*norm[0] + qlocal[6]*tau1[0] + qlocal[7]*tau2[0];
  qglobal[6] = qlocal[5]*norm[1] + qlocal[6]*tau1[1] + qlocal[7]*tau2[1];
  qglobal[7] = qlocal[5]*norm[2] + qlocal[6]*tau1[2] + qlocal[7]*tau2[2];
}

static inline void
rot_to_local_rect_glm(const double *tau1, const double *tau2, const double *norm,
  const double *GKYL_RESTRICT qglobal, double *GKYL_RESTRICT qlocal)
{
  rot_to_local_rect(tau1, tau2, norm, qglobal, qlocal);
  qlocal[8] = qglobal[8];
}

static inline void
rot_to_global_rect_glm(const double *tau1, const double *tau2, const double *norm,
  const double *GKYL_RESTRICT qlocal, double *GKYL_RESTRICT qglobal)
{
  rot_to_global_rect(tau1, tau2, norm, qlocal, qglobal);
  qglobal[8] = qlocal[8];
}

static void
mhd_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_wv_eqn *base = container_of(ref, struct gkyl_wv_eqn, ref_count);
  struct wv_mhd *mhd = container_of(base, struct wv_mhd, eqn);
  gkyl_free(mhd);
}

// Computing waves and waves speeds from Roe linearization.
// Following Cargo & Gallice 1997 section 4.2.
static double
wave_roe(const struct gkyl_wv_eqn *eqn,
  const double *dQ, const double *ql, const double *qr, double *waves, double *ev)
{
  const struct wv_mhd *mhd = container_of(eqn, struct wv_mhd, eqn);
  double gamma = mhd->gas_gamma;

  //////////////////////////////////////////////////////////////////////////////
  // STEP 1: COMPUTE PRIMITIVE VARIABLES                                      //
  //////////////////////////////////////////////////////////////////////////////
  double ul = ql[MX] / ql[DN], ur = qr[MX] / qr[DN];
  double vl = ql[MY] / ql[DN], vr = qr[MY] / qr[DN];
  double wl = ql[MZ] / ql[DN], wr = qr[MZ] / qr[DN];
  double pl = gkyl_mhd_pressure(gamma, ql);
  double pr = gkyl_mhd_pressure(gamma, qr);
  double pbl = 0.5 * (sq(ql[BX]) + sq(ql[BY]) + sq(ql[BZ]));
  double pbr = 0.5 * (sq(qr[BX]) + sq(qr[BY]) + sq(qr[BZ]));
  // total enthalpy (density) in CG97 eq. 2.2
  double Hl = (ql[ER] + pl + pbl) / ql[DN];
  double Hr = (qr[ER] + pr + pbr) / qr[DN];

  //////////////////////////////////////////////////////////////////////////////
  // STEP 2: COMPUTE ROE AVERAGES OF PRIMITIVE VARIABLES                      //
  //////////////////////////////////////////////////////////////////////////////
  double srrhol = sqrt(ql[DN]);
  double srrhor = sqrt(qr[DN]);
  double sl = srrhol / (srrhol + srrhor);
  double sr = srrhor / (srrhol + srrhor);

  double rho = srrhol * srrhor;
  double u = sl * ul + sr * ur;
  double v = sl * vl + sr * vr;
  double w = sl * wl + sr * wr;
  double H = sl * Hl + sr * Hr;  // total enthalpy
  double Bx = sr * ql[BX] + sl * qr[BX];
  double By = sr * ql[BY] + sl * qr[BY];
  double Bz = sr * ql[BZ] + sl * qr[BZ];

  //////////////////////////////////////////////////////////////////////////////
  // STEP 3: COMPUTE CHARACTERASTIC WAVE SPEEDS AND OTHER USEFUL QUANTITIES   //
  //////////////////////////////////////////////////////////////////////////////
  // CG97 eq. 4.12
  double X = (sq(dQ[BY]) + sq(dQ[BZ])) / (2*sq(srrhol+srrhor));

  // CG97 eq 4.17, wave speeds
  double ca2 = Bx*Bx/rho; // for alfven speed due to normal B field
  double b2 = (Bx*Bx+By*By+Bz*Bz) / rho; // for alfven speed due to full B field
  double v2 = u*u+v*v+w*w;
  double Hgas = H - b2;  // enthalpy of the gas
  double a2 = (2-gamma)*X + (gamma-1)*(Hgas-0.5*v2);  // for sound speed

  double astar2 = a2 + b2;
  double cf2 = (astar2 + sqrt(sq(astar2)-4*a2*ca2)) / 2;  // fast wave speed
  double cs2 = (astar2 - sqrt(sq(astar2)-4*a2*ca2)) / 2;  // slow wave speed

  double a = sqrt(a2);  // sound speed
  double ca = sqrt(ca2);  // alfven speed due to normal B field
  double cf = sqrt(cf2);  // fast magnetosonic speed
  double cs = sqrt(cs2);  // slow magnetosonic speed

  // CG97 eq 4.21, S, beta, and alpha
  double S = Bx >=0? 1: -1;

  double Bt = sqrt(By*By+Bz*Bz);
  double betay, betaz;
  if (Bt>0) { // TODO compare wtih a tiny number
    betay = By / Bt;
    betaz = Bz / Bt;
  } else {
    betay = betaz = 1/sqrt(2);
  }

  double alphaf2 = (a2 - cs2) / (cf2 - cs2);
  double alphas2 = (cf2 - a2) / (cf2 - cs2);
  double alphaf = sqrt(alphaf2);
  double alphas = sqrt(alphas2);

  //////////////////////////////////////////////////////////////////////////////
  // STEP 4: COMPUTE WAVE SPEEDS AND FLUCTUATIONS DUE TO ROE LINEARIZATION    //
  //                                                                          //
  // For each eigensolution of the Roe linearizing matrix, compute:           //
  // 0. Eigenvalue, which is the wave speed. The results are stored in ev     //
  //     array and sorted by their values, that is,                           //
  //    u-c_fast<u-c_alfven<u-c_slow<u<u+c_slow<u+c_alfven<u+c_fast           //
  // 1. Projection coefficient of jumps in conservative state on to the       //
  //    left-eigenvector according to CG97 eq. 4.20.                          //
  // 2. The wave (fluctuations) due to this mode.                             //
  //////////////////////////////////////////////////////////////////////////////
  // CG97 eq. 4.20 gave projection coefficients in primitive varialbles.
  // Alternatively, one may compute the full left-eigenvectors and dot it with
  // dQ. The left-eigenvector of the Roe matrix can be found in Stone et al.
  // (2008), Appendix B2.
  double drho = dQ[DN];
  double du = ur - ul;
  double dv = vr - vl;
  double dw = wr - wl;
  double dp = pr - pl;
  double dBy = dQ[BY];
  double dBz = dQ[BZ];

  const int meqns = eqn->num_equations;
  const int mwaves = eqn->num_waves;
  for (int i=0; i<meqns*mwaves; ++i) waves[i] = 0.0;
  double *wv;
  double eta[mwaves];

  /////////////////////////////
  // Fast magnetosonic waves //
  /////////////////////////////

  // For projection coefficients in eq. 4.20 for magnetosonic waves
  double t1 = alphaf * (X*drho + dp);
  double t2 = alphas * cs * rho * S * (betay*dv + betaz*dw);
  double t3 = alphaf * cf * rho * du;
  double t4 = alphas * sqrt(rho) * a * (betay*dBy + betaz*dBz);
  // For right eigenvectors for magnetosonic waves in eq. 4.19
  double t5 = alphas * cs * S * (v*betay + w*betaz);
  double t6 = alphas * a * Bt / sqrt(rho);

  ev[0] = u - cf;
  // Coefficients of projection onto the left-eigenvector
  eta[0] = (t1 + t2 - t3 + t4) / 2; // eq. 4.20
  eta[0] /= a2; // merge a2 in right-eigenvector in 4.19 into coefficents
  // Fluctuations due to this mode
  wv = &waves[0*meqns];
  wv[DN] = eta[0]*alphaf;
  wv[MX] = eta[0]*(alphaf*(u - cf));
  wv[MY] = eta[0]*(alphaf*v + alphas*cs*betay*S);
  wv[MZ] = eta[0]*(alphaf*w + alphas*cs*betaz*S);
  wv[ER] = eta[0]*(alphaf*(Hgas-u*cf)+t5+t6); // -t6 in CG97 typo? +t6 in SG08
  wv[BY] = eta[0]*alphas*a*betay/sqrt(rho);
  wv[BZ] = eta[0]*alphas*a*betaz/sqrt(rho);

  ev[6] = u + cf;
  // Coefficients of projection onto the left-eigenvector
  eta[6] = (t1 - t2 + t3 + t4) / 2; // eq. 4.20
  eta[6] /= a2; // merge a2 in right-eigenvector in 4.19 into coefficents
  // Fluctuations due to this mode
  wv = &waves[6*meqns];
  wv[DN] = eta[6]*alphaf;
  wv[MX] = eta[6]*(alphaf*(u + cf));
  wv[MY] = eta[6]*(alphaf*v - alphas*cs*betay*S);
  wv[MZ] = eta[6]*(alphaf*w - alphas*cs*betaz*S);
  wv[ER] = eta[6]*(alphaf*(Hgas+u*cf)-t5+t6); // -t6 in CG97 typo? +t6 in SG08
  wv[BY] = eta[6]*alphas*a*betay/sqrt(rho);
  wv[BZ] = eta[6]*alphas*a*betaz/sqrt(rho);

  //////////////////
  // Alfven waves //
  //////////////////

  // For projection coefficients in eq. 4.20 for alfven waves
  double t7 = betay*dw - betaz*dv;
  double t8 = (S/sqrt(rho)) * (betay*dBz-betaz*dBy);

  ev[1] = u - ca;
  // Coefficients of projection onto the left-eigenvector
  eta[1] = (t7 + t8) / 2; // eg. 4.20
  eta[1] *= rho; // merge rho in right-eigenvector in 4.18 into coefficents
  wv = &waves[1*meqns];
  // Fluctuations due to this mode
  wv[MY] = -eta[1]*betaz;
  wv[MZ] = eta[1]*betay;
  wv[ER] = -eta[1]*(v*betaz - w*betay);
  wv[BY] = -eta[1]*S*betaz/sqrt(rho);
  wv[BZ] = eta[1]*S*betay/sqrt(rho);

  ev[5] = u + ca;
  // Coefficients of projection onto the left-eigenvector
  eta[5] = (-t7 + t8) / 2; // eq. 4.20
  eta[5] *= rho; // merge rho in right-eigenvector in 4.18 into coefficents
  wv = &waves[5*meqns];
  // Fluctuations due to this mode
  wv[MY] = eta[5]*betaz;
  wv[MZ] = -eta[5]*betay;
  wv[ER] = +eta[5]*(v*betaz - w*betay);
  wv[BY] = -eta[5]*S*betaz/sqrt(rho);
  wv[BZ] = eta[5]*S*betay/sqrt(rho);

  /////////////////////////////
  // Slow magnetosonic waves //
  /////////////////////////////

  // For projection coefficients in eq. 4.20 for magnetosonic waves
  // Compared to the fast wave coefficients, alphas/alphaf & cs/cf are switched
  t1 = alphas * (X*drho + dp);
  t2 = alphaf * cf * rho * S * (betay*dv + betaz*dw);
  t3 = alphas * cs * rho * du;
  t4 = alphaf * sqrt(rho) * a * (betay*dBy + betaz*dBz);
  // for right eigenvectors for magnetosonic waves in eq. 4.19
  t5 = alphaf * cf * S * (v*betay + w*betaz);
  t6 = alphaf * a * Bt / sqrt(rho);

  ev[2] = u - cs;
  // Coefficients of projection onto the left-eigenvector
  eta[2] = (t1 - t2 - t3 - t4) / 2; // eq. 4.20
  eta[2] /= a2; // merge a2 in right-eigenvector in 4.19 into coefficents
  // Fluctuations due to this mode
  wv = &waves[2*meqns];
  wv[DN] = eta[2]*alphas;
  wv[MX] = eta[2]*(alphas*(u - cs));
  wv[MY] = eta[2]*(alphas*v - alphaf*cf*betay*S);
  wv[MZ] = eta[2]*(alphas*w - alphaf*cf*betaz*S);
  wv[ER] = eta[2]*(alphas*(Hgas-u*cs) - t5 - t6);
  wv[BY] = -eta[2]*alphaf*a*betay/sqrt(rho);
  wv[BZ] = -eta[2]*alphaf*a*betaz/sqrt(rho);

  ev[4] = u + cs;
  // Coefficients of projection onto the left-eigenvector
  eta[4] = (t1 + t2 + t3 - t4) / 2; // eq. 4.20
  eta[4] /= a2; // merge a2 in right-eigenvector in 4.19 into coefficents
  // Fluctuations due to this mode
  wv = &waves[4*meqns];
  wv[DN] = eta[4]*alphas;
  wv[MX] = eta[4]*(alphas*(u + cs));
  wv[MY] = eta[4]*(alphas*v + alphaf*cf*betay*S);
  wv[MZ] = eta[4]*(alphas*w + alphaf*cf*betaz*S);
  wv[ER] = eta[4]*(alphas*(Hgas+u*cs) + t5 - t6);
  wv[BY] = -eta[4]*alphaf*a*betay/sqrt(rho);
  wv[BZ] = -eta[4]*alphaf*a*betaz/sqrt(rho);

  //////////////////
  // Enropy wave  //
  //////////////////
  ev[3] = u;
  // Coefficients of projection onto the left-eigenvector
  eta[3] = (a2 - X) * drho - dp; // eq. 4.20
  eta[3] /= a2; // merge a2 in right-eigenvector in 4.18 into coefficents
  // Fluctuations due to this mode
  wv = &waves[3*meqns];
  wv[DN] = eta[3];
  wv[MX] = eta[3]*u;
  wv[MY] = eta[3]*v;
  wv[MZ] = eta[3]*w;
  wv[ER] = eta[3]*(v2/2 + X*(gamma-2)/(gamma-1));

  /////////////////
  // div(B) wave //
  /////////////////
  // This wave exists in the eight-wave scheme. In this implementation, it is
  // incorporated into the last wave, the entropy wave, since both have eigen
  // value ev=u. This wave only advects jump in Bx at the speed u.
  if (mhd->divergence_constraint == GKYL_MHD_DIVB_EIGHT_WAVES)
    wv[BX] = dQ[BX];

  double max_speed =  fabs(u) + cf;

  // GLM divB and psi waves; XXX set ch properly.
  if (mhd->divergence_constraint == GKYL_MHD_DIVB_GLM)
  {
    double ch;
    if (mhd->glm_ch > 0)
      ch = mhd->glm_ch;
    else
      // candidates:
      // 1. Dedner (2002) sec 4: cfl * min(dx, dy, dz) / dt
      // 2. Derigs (2018) eq 3.30: max_speed_global - max_global(|u|, |v|, |w|)
      // 3. local version of 2; needs to store ch for the source step
      ch = max_speed;

    ev[7] = -ch;
    eta[7] = 0.5 * (-dQ[BX]*ch+dQ[PSI_GLM]);
    wv = &waves[7*meqns];
    wv[BX] = -eta[7]/ch;
    wv[PSI_GLM] = eta[7];

    ev[8] = ch;
    eta[8] = 0.5 * (dQ[BX]*ch+dQ[PSI_GLM]);
    wv = &waves[8*meqns];
    wv[BX] = eta[8]/ch;
    wv[PSI_GLM] = eta[8];

    max_speed = max_speed > ch ? max_speed : ch;
  }

  return max_speed;
}

static void
qfluct_roe(const struct gkyl_wv_eqn *eqn,
  const double *ql, const double *qr, const double *waves, const double *s, double *amdq,
  double *apdq)
{
  int meqn = eqn->num_equations;
  for (int i=0; i<meqn; ++i) {
    amdq[i] = fmin(0.0, s[0]) * waves[i];
    apdq[i] = fmax(0.0, s[0]) * waves[i];
    for (int mw=1; mw<eqn->num_waves; ++mw) {
      amdq[i] += fmin(0.0, s[mw]) * waves[meqn*mw+i];
      apdq[i] += fmax(0.0, s[mw]) * waves[meqn*mw+i];
    }
  }
}

// Computing waves and waves speeds from Lax fluxes
static double
wave_lax(const struct gkyl_wv_eqn *eqn,
  const double *dQ, const double *ql, const double *qr, double *waves, double *ev)
{
  const struct wv_mhd *mhd = container_of(eqn, struct wv_mhd, eqn);
  double gas_gamma = mhd->gas_gamma;  

  double sl = gkyl_mhd_max_abs_speed(gas_gamma, ql);
  double sr = gkyl_mhd_max_abs_speed(gas_gamma, qr);

  double *wv = &waves[0]; // single wave
  for (int i=0; i<mhd->eqn.num_equations; ++i)  wv[i] = dQ[i];
  
  ev[0] = 0.5*(sl+sr);

  return ev[0];
}

static void
qfluct_lax(const struct gkyl_wv_eqn *eqn,
  const double *ql, const double *qr, const double *waves, const double *s,
  double *amdq, double *apdq)
{
  const struct wv_mhd *mhd = container_of(eqn, struct wv_mhd, eqn);
  double gas_gamma = mhd->gas_gamma;  

  double sl = gkyl_mhd_max_abs_speed(gas_gamma, ql);
  double sr = gkyl_mhd_max_abs_speed(gas_gamma, qr);
  double amax = fmax(sl, sr);

  double fl[10], fr[10];
  
  if (mhd->divergence_constraint == GKYL_MHD_DIVB_GLM) {
    double ch = amax; // Is this right?
    gkyl_glm_mhd_flux(gas_gamma, ch, ql, fl);
    gkyl_glm_mhd_flux(gas_gamma, ch, qr, fr);
  } else {
    gkyl_mhd_flux(gas_gamma, ql, fl);
    gkyl_mhd_flux(gas_gamma, qr, fr);
  }

  for (int i=0; i<mhd->eqn.num_equations; ++i) {
    amdq[i] = 0.5*(fr[i]-fl[i] - amax*(qr[i]-ql[i]));
    apdq[i] = 0.5*(fr[i]-fl[i] + amax*(qr[i]-ql[i]));
  }
}

static double
wave(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *delta, const double *ql, const double *qr, double *waves, double *s)
{
  if (type == GKYL_WV_HIGH_ORDER_FLUX)
    return wave_roe(eqn, delta, ql, qr, waves, s);
  else
    return wave_lax(eqn, delta, ql, qr, waves, s);

  return 0.0; // can't happen
}

static void
qfluct(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *ql, const double *qr, const double *waves, const double *s,
  double *amdq, double *apdq)
{
  if (type == GKYL_WV_HIGH_ORDER_FLUX)
    return qfluct_roe(eqn, ql, qr, waves, s, amdq, apdq);
  else
    return qfluct_lax(eqn, ql, qr, waves, s, amdq, apdq);
}

static bool
check_inv(const struct gkyl_wv_eqn *eqn, const double *q)
{
  const struct wv_mhd *mhd = container_of(eqn, struct wv_mhd, eqn);
  
  if (q[0] < 0.0)
    return false;

  double pr = gkyl_mhd_pressure(mhd->gas_gamma, q);
  if (pr < 0.0)
    return false;

  return true;
}

static double
max_speed(const struct gkyl_wv_eqn *eqn, const double *q)
{
  const struct wv_mhd *mhd = container_of(eqn, struct wv_mhd, eqn);
  return gkyl_mhd_max_abs_speed(mhd->gas_gamma, q);
}

struct gkyl_wv_eqn*
gkyl_wv_mhd_new(const struct wv_mhd_inp *inp)
{
  struct wv_mhd *mhd = gkyl_malloc(sizeof(struct wv_mhd));

  mhd->eqn.type = GKYL_EQN_MHD;
  mhd->gas_gamma = inp->gas_gamma;
  mhd->eqn.waves_func = wave;
  mhd->eqn.qfluct_func = qfluct;
  
  mhd->eqn.check_inv_func = check_inv;
  mhd->eqn.max_speed_func = max_speed;
  
  mhd->eqn.rotate_to_local_func = rot_to_local_rect;
  mhd->eqn.rotate_to_global_func = rot_to_global_rect;

  mhd->divergence_constraint = inp->divergence_constraint;
  switch (inp->divergence_constraint) {
    case GKYL_MHD_DIVB_NONE:
      mhd->eqn.num_equations = 8;
      mhd->eqn.num_waves = 7;
      break;

    case GKYL_MHD_DIVB_EIGHT_WAVES:
      mhd->eqn.num_equations = 8;
      mhd->eqn.num_waves = 7;  // will merge the entropy wave and divB wave
      break;

    case GKYL_MHD_DIVB_GLM:
      mhd->eqn.num_equations = 9;
      mhd->eqn.num_waves = 9;
      mhd->eqn.rotate_to_local_func = rot_to_local_rect_glm;
      mhd->eqn.rotate_to_global_func = rot_to_global_rect_glm;
      mhd->glm_ch = inp->glm_ch;
      mhd->glm_cp = inp->glm_cp;
      break;
  }

  mhd->eqn.ref_count = gkyl_ref_count_init(mhd_free);

  return &mhd->eqn;
}

double
gkyl_wv_mhd_gas_gamma(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_mhd *mhd = container_of(eqn, struct wv_mhd, eqn);
  return mhd->gas_gamma;
}
