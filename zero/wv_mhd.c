#include <math.h>
#include <float.h>

#include <gkyl_alloc.h>
#include <gkyl_prim_mhd.h>
#include <gkyl_wv_mhd.h>

static const int dir_shuffle[][6] = {
  {1, 2, 3, 5, 6, 7},
  {2, 3, 1, 6, 7, 5},
  {3, 1, 2, 7, 5, 6}
};

// Make indexing cleaner with the dir_shuffle
#define DN (0)
#define MX (idx[0])
#define MY (idx[1])
#define MZ (idx[2])
#define ER (4)
#define BX (idx[3])
#define BY (idx[4])
#define BZ (idx[5])

#define sq(x) ((x)*(x))

static inline void
rot_to_local_rect(int dir, const double *tau1, const double *tau2, const double *norm,
  const double *GKYL_RESTRICT qglobal, double *GKYL_RESTRICT qlocal)
{
  const int *idx = dir_shuffle[dir];
  qlocal[0] = qglobal[0];
  qlocal[1] = qglobal[MX];
  qlocal[2] = qglobal[MY];
  qlocal[3] = qglobal[MZ];
  qlocal[4] = qglobal[4];
  qlocal[5] = qglobal[BX];
  qlocal[6] = qglobal[BY];
  qlocal[7] = qglobal[BZ];
}

static inline void
rot_to_global_rect(int dir, const double *tau1, const double *tau2, const double *norm,
  const double *GKYL_RESTRICT qlocal, double *GKYL_RESTRICT qglobal)
{
  const int *idx = dir_shuffle[dir];  
  qglobal[0] = qlocal[0];
  qglobal[MX] = qlocal[1];
  qglobal[MY] = qlocal[2];
  qglobal[MZ] = qlocal[3];
  qglobal[4] = qlocal[4];
  qglobal[BX] = qlocal[5];
  qglobal[BY] = qlocal[6];
  qglobal[BZ] = qlocal[7];
}

struct wv_mhd {
  struct gkyl_wv_eqn eqn; // base object
  double gas_gamma; // gas adiabatic constant
};

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
wave_roe(const struct gkyl_wv_eqn *eqn, int dir, const double *dQ,
  const double *ql, const double *qr, double *waves, double *ev)
{
  const struct wv_mhd *mhd = container_of(eqn, struct wv_mhd, eqn);
  const int *idx = dir_shuffle[0]; // no need for shuffle
  double gamma = mhd->gas_gamma;

  // Compute primitive states
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

  // Compute Roe averages of primitive states
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

  // CG97 eq. 4.12; FIXME: dBx is included as needed by 8-wave scheme
  double X = (sq(dQ[BX]) + sq(dQ[BY]) + sq(dQ[BZ])) / (2*sq(srrhol+srrhor));

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
  if (Bt>0) {
    betay = By / Bt;
    betaz = Bz / Bt;
  } else {
    betay = betaz = 1/sqrt(2);
  }

  double alphaf2 = (a2 - cs2) / (cf2 - cs2);
  double alphas2 = (cf2 - a2) / (cf2 - cs2);
  double alphaf = sqrt(alphaf2);
  double alphas = sqrt(alphas2);

  // For each eigensolution of the Roe linearizing matrix, compute:
  // 0. Eigenvalue, which is the wave speed. The results are stored in ev array
  //    and sorted by their values, that is,
  //    u-c_fast < u-c_alfven < u-c_slow < u < u+c_slow < u+c_alfven < u+c_fast
  // 1. Projection coefficient of jumps in conservative state on to the
  //    left-eigenvector according to CG97 eq. 4.20.
  // 2. The wave (fluctuations) due to this mode.
  double drho = dQ[DN];
  double du = ur - ul;
  double dv = vr - vl;
  double dw = wr - wl;
  double dp = pr - pl;
  double dBy = dQ[BY];
  double dBz = dQ[BZ];

  const int meqns = 8;
  for (int i=0; i<meqns*7; ++i) waves[i] = 0.0;
  double *wv;
  double eta[7];

  /////////////////////////////////////////////////////////////
  // Fast magnetosonic modes

  // For projection coefficients in eq. 4.20 for magnetosonic waves
  double t1 = alphaf * (X*drho + dp);
  double t2 = alphas * cs * rho * S * (betay*dv + betaz*dw);
  double t3 = alphaf * cf * rho * du;
  double t4 = alphas * sqrt(rho) * a * (betay*dBy + betaz*dBz);
  // for right eigenvectors for magnetosonic waves in eq. 4.19
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

  /////////////////////////////////////////////////////////////
  // Alfven waves

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

  /////////////////////////////////////////////////////////////
  // Slow magnetosonic waves

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

  /////////////////////////////////////////////////////////////
  // Enropy wave; TODO: merge divB wave
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

  /* printf("\n"); */
  /* printf("u=%f, a=%f, ca=%f, cs=%f, cf=%f\n", u, a, ca, cs, cf); */
  /* printf("alphaf=%f, alphas=%f, betay=%f, betaz=%f\n", */
  /*        alphaf, alphas, betay, betaz); */
  /* for(int i=0;i<7;++i) */
  /*   printf("[%d], eta=%10g, ev=%10g\n", i, eta[i], ev[i]); */

  return fabs(u) + cf;
}

static void
qfluct_roe(const struct gkyl_wv_eqn *eqn, int dir, const double *ql,
  const double *qr, const double *waves, const double *s, double *amdq,
  double *apdq)
{
  const double *w0 = &waves[0*8], *w1 = &waves[1*8], *w2 = &waves[2*8];
  const double *w3 = &waves[3*8], *w4 = &waves[4*8], *w5 = &waves[5*8];
  const double *w6 = &waves[6*8];

  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]), s2m = fmin(0.0, s[2]);
  double s3m = fmin(0.0, s[3]), s4m = fmin(0.0, s[4]), s5m = fmin(0.0, s[5]);
  double s6m = fmin(0.0, s[6]);

  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]), s2p = fmax(0.0, s[2]);
  double s3p = fmax(0.0, s[3]), s4p = fmax(0.0, s[4]), s5p = fmax(0.0, s[5]);
  double s6p = fmax(0.0, s[6]);

  for (int i=0; i<8; ++i) {
    amdq[i] = s0m*w0[i] + s1m*w1[i] + s2m*w2[i] + s3m*w3[i] + s4m*w4[i]
      + s5m*w5[i] + s5m*w6[i];
    apdq[i] = s0p*w0[i] + s1p*w1[i] + s2p*w2[i] + s3p*w3[i] + s4p*w4[i]
      + s5p*w5[i] + s6p*w6[i];
  }
}

static double
max_speed(const struct gkyl_wv_eqn *eqn, int dir, const double *q)
{
  const struct wv_mhd *mhd = container_of(eqn, struct wv_mhd, eqn);
  return gkyl_mhd_max_abs_speed(dir, mhd->gas_gamma, q);
}

struct gkyl_wv_eqn*
gkyl_wv_mhd_new(double gas_gamma)
{
  struct wv_mhd *mhd = gkyl_malloc(sizeof(struct wv_mhd));

  mhd->eqn.type = GKYL_EQN_MHD;
  mhd->eqn.num_equations = 8;
  mhd->eqn.num_waves = 7;
  mhd->gas_gamma = gas_gamma;
  mhd->eqn.waves_func = wave_roe;
  mhd->eqn.qfluct_func = qfluct_roe;
  mhd->eqn.max_speed_func = max_speed;
  mhd->eqn.rotate_to_local_func = rot_to_local_rect;
  mhd->eqn.rotate_to_global_func = rot_to_global_rect;

  mhd->eqn.ref_count = (struct gkyl_ref_count) { mhd_free, 1 };

  return &mhd->eqn;
}

double
gkyl_wv_mhd_gas_gamma(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_mhd *mhd = container_of(eqn, struct wv_mhd, eqn);
  return mhd->gas_gamma;
}
