#include <math.h>
#include <float.h>

#include <gkyl_moment_prim_mhd.h>

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

double
gkyl_mhd_fast_speed(double gas_gamma, const double q[8])
{
  double u1 = q[MX] / q[DN];
  double u2 = q[MY] / q[DN];
  double u3 = q[MZ] / q[DN];
  double k = q[DN] * (u1*u1 + u2*u2 + u3*u3) / 2; // bulk kinetic energy
  double BX_sq = sq(q[BX]);
  double B_sq = BX_sq + sq(q[BY]) + sq(q[BZ]);
  double pb = B_sq / 2; // magnetic pressure
  double p = (gas_gamma-1) * (q[ER] - k - pb); // plasma pressure

  double a_sq = gas_gamma * p / q[DN]; // sound speed
  double ca_sq = B_sq / q[DN];  // Alfven speed
  double ca1_sq = BX_sq / q[DN];  // Alfven speed due to normal B field
  // fast speed
  double cf = sqrt(a_sq+ca_sq + sqrt(sq(a_sq + ca_sq) - 4*a_sq*ca1_sq)) / 2;

  return cf;
}

double
gkyl_mhd_max_abs_speed(double gas_gamma, const double q[8])
{
  double cf = gkyl_mhd_fast_speed(gas_gamma, q);
  double u = sqrt( sq(q[1]) + sq(q[2]) + sq(q[3]) ) / q[0];

  return u + cf;
}

double
gkyl_mhd_max_abs_speed_roe(const double gamma, const double *ql, const double *qr)
{
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
  double X = (sq(qr[BX]-ql[BX]) + sq(qr[BY]-ql[BY]) + sq(qr[BZ]-qr[BZ])) / (2*sq(srrhol+srrhor));

  //////////////////////////////////////////////////////////////////////////////
  // STEP 3: COMPUTE CHARACTERASTIC WAVE SPEEDS AND OTHER USEFUL QUANTITIES   //
  //////////////////////////////////////////////////////////////////////////////
  double ca2 = Bx*Bx/rho; // for alfven speed due to normal B field
  double b2 = (Bx*Bx+By*By+Bz*Bz) / rho; // for alfven speed due to full B field
  double v2 = u*u+v*v+w*w;
  double Hgas = H - b2;  // enthalpy of the gas
  double a2 = (2-gamma)*X + (gamma-1)*(Hgas-0.5*v2);  // for sound speed

  double astar2 = a2 + b2;
  double cf2 = (astar2 + sqrt(sq(astar2)-4*a2*ca2)) / 2;  // fast wave speed

  return sqrt(cf2) + sqrt(v2);
}

void
gkyl_mhd_flux(double gas_gamma, const double q[8], double flux[8])
{
  double u1 = q[MX] / q[DN];
  double u2 = q[MY] / q[DN];
  double u3 = q[MZ] / q[DN];
  double k = q[DN] * (u1*u1 + u2*u2 + u3*u3) / 2;  // bulk kinetic energy
  double pb = (sq(q[BX]) + sq(q[BY]) + sq(q[BZ])) / 2; // magnetic pressure
  double p = (gas_gamma-1) * (q[ER] - k - pb); // plasma pressure

  flux[DN] = q[MX];
  flux[MX] = u1*q[MX] - q[BX]*q[BX] + p + pb;
  flux[MY] = u1*q[MY] - q[BX]*q[BY];
  flux[MZ] = u1*q[MZ] - q[BX]*q[BZ];
  flux[ER] = u1*(q[ER]+p+pb) - q[BX]*(u1*q[BX]+u2*q[BY]+u3*q[BZ]);
  flux[BX] = 0.0;
  flux[BY] = u1*q[BY] - u2*q[BX];
  flux[BZ] = u1*q[BZ] - u3*q[BX];
}

void
gkyl_glm_mhd_flux(double gas_gamma, double ch, const double q[9], double flux[9])
{
  gkyl_mhd_flux(gas_gamma, q, flux);

  flux[BX] = q[PSI_GLM];
  flux[PSI_GLM] = ch*ch*q[BX];
}
