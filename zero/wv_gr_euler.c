#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
//#include <gkyl_wv_gr_euler.h>
#include <gkyl_wv_gr_euler_priv.h>

static inline void
gkyl_gr_euler_prim_vars(double gas_gamma, const double q[5], double v[5])
{
  double D = q[0];
  double momx = q[1];
  double momy = q[2];
  double momz = q[3];
  double Etot = q[4];

  double C = D / sqrt(((Etot + D) * (Etot + D)) - ((momx * momx) + (momy * momy) + (momz * momz)));
  double C0 = (D + Etot) / sqrt(((Etot + D) * (Etot + D)) - ((momx * momx) + (momy * momy) + (momz * momz)));
  if (((Etot + D) * (Etot + D)) - ((momx * momx) + (momy * momy) + (momz * momz)) < pow(10.0, -8.0)) {
    C = D / sqrt(pow(10.0, -8.0));
    C0 = (D + Etot) / sqrt(pow(10.0, -8.0));
  }

  double alpha0 = -1.0 / (gas_gamma * gas_gamma);
  double alpha1 = -2.0 * C * ((gas_gamma - 1.0) / (gas_gamma * gas_gamma));
  double alpha2 = ((gas_gamma - 2.0) / gas_gamma) * ((C0 * C0) - 1.0) + 1.0 - (C * C) * ((gas_gamma - 1.0) / gas_gamma) * ((gas_gamma - 1.0) / gas_gamma);
  double alpha4 = (C0 * C0) - 1.0;
  double eta = 2.0 * C *((gas_gamma - 1.0) / gas_gamma);

  double guess = 1.0;
  int iter = 0;

  while (iter < 1000) {
    double poly = (alpha4 * (guess * guess * guess) * (guess - eta)) + (alpha2 * (guess * guess)) + (alpha1 * guess) + alpha0;
    double poly_der = alpha1 + (2.0 * alpha2 * guess) + (4.0 * alpha4 * (guess * guess * guess)) - (3.0 * eta * alpha4 * (guess * guess));

    double guess_new = guess - (poly / poly_der);

    if (fabs(guess - guess_new) < pow(10.0, -8.0)) {
      iter = 1000;
    }
    else {
      iter += 1;
      guess = guess_new;
    }
  }

  double W = 0.5 * C0 * guess * (1.0 + sqrt(1.0 + (4.0 * ((gas_gamma - 1.0) / gas_gamma) * ((1.0 - (C * guess)) / ((C0 * C0) * (guess * guess))))));
  double h = 1.0 / (C * guess);

  v[0] = D / W;
  v[1] = momx / (D * h * (W * W));
  v[2] = momy / (D * h * (W * W));
  v[3] = momz / (D * h * (W * W));
  v[4] = (D * h * (W * W)) - D - Etot;

  if (v[0] < pow(10.0, -8.0)) {
    v[0] = pow(10.0, -8.0);
  }
  if (v[4] < pow(10.0, -8.0)) {
    v[4] = pow(10.0, -8.0);
  }
}

static inline double
gkyl_gr_euler_max_abs_speed(double gas_gamma, const double q[5])
{
  double v[5] = { 0.0 };
  gkyl_gr_euler_prim_vars(gas_gamma, q, v);
  double rho = v[0];
  double p = v[4];

  double num = (gas_gamma * p) / rho;
  double den = 1.0 * ((p / rho) * (gas_gamma / (gas_gamma - 1.0)));
  double cs = sqrt(num / den);

  double v_sq = sqrt((v[1] * v[1]) + (v[2] * v[2]) + (v[3] * v[3]));
  return fabs(v_sq) + cs;
}

static void
gkyl_gr_euler_flux(double gas_gamma, const double q[5], double flux[5])
{
  double v[5] = { 0.0 };
  gkyl_gr_euler_prim_vars(gas_gamma, q, v);
  double rho =  v[0];
  double vx = v[1];
  double vy = v[2];
  double vz = v[3];
  double p = v[4];

  double W = 1.0 / (sqrt(1.0 - ((vx * vx) + (vy * vy) + (vz * vz))));
  double h = 1.0 + ((p / rho) * (gas_gamma / (gas_gamma - 1.0)));

  flux[0] = rho * W * vx;
  flux[1] = rho * h * (W * W) * (vx * vx) + p;
  flux[2] = rho * h * (W * W) * (vx * vy);
  flux[3] = rho * h * (W * W) * (vx * vz);
  flux[4] = ((rho * h * (W * W)) - p - (rho * W)) * vx + (p * vx);
}

static inline void
cons_to_riem(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* qin, double* wout)
{
  // TODO: This should use a proper L matrix.
  for (int i = 0; i < 5; i++) {
    wout[i] = qin[i];
  }
}

static inline void
riem_to_cons(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* win, double* qout)
{
  // TODO: This should use a proper L matrix.
  for (int i = 0; i < 5; i++) {
    qout[i] = win[i];
  }
}

static void
gr_euler_wall(double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  for (int i = 0; i < 5; i++) {
    ghost[i] = skin[i];
  }

  ghost[1] = -ghost[1];
}

static void
gr_euler_no_slip(double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  for (int i = 1; i < 4; i++) {
    ghost[i] = -skin[i];
  }

  ghost[0] = skin[0];
  ghost[4] = skin[4];
}

static inline void
rot_to_local(const double* tau1, const double* tau2, const double* norm, const double* GKYL_RESTRICT qglobal, double* GKYL_RESTRICT qlocal)
{
  qlocal[0] = qglobal[0];
  qlocal[1] = (qglobal[1] * norm[0]) + (qglobal[2] * norm[1]) + (qglobal[3] * norm[2]);
  qlocal[2] = (qglobal[1] * tau1[0]) + (qglobal[2] * tau1[1]) + (qglobal[3] * tau1[2]);
  qlocal[3] = (qglobal[1] * tau2[0]) + (qglobal[2] * tau2[1]) + (qglobal[3] * tau2[2]);
  qlocal[4] = qglobal[4];
}

static inline void
rot_to_global(const double* tau1, const double* tau2, const double* norm, const double* GKYL_RESTRICT qlocal, double* GKYL_RESTRICT qglobal)
{
  qglobal[0] = qlocal[0];
  qglobal[1] = (qlocal[1] * norm[0]) + (qlocal[2] * tau1[0]) + (qlocal[3] * tau2[0]);
  qglobal[2] = (qlocal[1] * norm[1]) + (qlocal[2] * tau1[1]) + (qlocal[3] * tau2[1]);
  qglobal[3] = (qlocal[1] * norm[2]) + (qlocal[2] * tau1[2]) + (qlocal[3] * tau2[2]);
  qglobal[4] = qlocal[4];
}