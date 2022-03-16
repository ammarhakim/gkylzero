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
gkyl_mhd_max_abs_speed(double gas_gamma, const double q[8])
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

  return fabs(u1) + cf;
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
