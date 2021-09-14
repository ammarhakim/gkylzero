#include <math.h>
#include <float.h>

#include <gkyl_prim_mhd.h>

static const int dir_shuffle[][6] = {
  {1, 2, 3, 5, 6, 7},
  {2, 3, 1, 6, 7, 5},
  {3, 1, 2, 7, 5, 6}
};

#define RHO (0)
#define M1 (d[0])
#define M2 (d[1])
#define M3 (d[2])
#define ER (4)
#define B1 (d[3])
#define B2 (d[4])
#define B3 (d[5])

#define sq(x) ((x)*(x))

double
gkyl_mhd_max_abs_speed(int dir, double gas_gamma, const double q[8])
{
  int const *const d = dir_shuffle[dir];

  double u1 = q[M1] / q[RHO];
  double u2 = q[M2] / q[RHO];
  double u3 = q[M3] / q[RHO];
  double k = q[RHO] * (u1*u1 + u2*u2 + u3*u3) / 2; // bulk kinetic energy
  double B1_sq = sq(q[B1]);
  double B_sq = B1_sq + sq(q[B2]) + sq(q[B3]);
  double pb = B_sq / 2; // magnetic pressure
  double p = (gas_gamma-1) * (q[ER] - k - pb); // plasma pressure

  double a_sq = gas_gamma * p / q[RHO]; // sound speed
  double ca_sq = B_sq / q[RHO];  // Alfven speed
  double ca1_sq = B1_sq / q[RHO];  // Alfven speed due to normal B field
  // fast speed
  double cf = sqrt(a_sq+ca_sq + sqrt(sq(a_sq + ca_sq) - 4*a_sq*ca1_sq)) / 2;

  return fabs(u1) + cf;
}

void
gkyl_mhd_flux(int dir, double gas_gamma, const double q[8], double flux[8])
{
  int const *const d = dir_shuffle[dir];

  double u1 = q[M1] / q[RHO];
  double u2 = q[M2] / q[RHO];
  double u3 = q[M3] / q[RHO];
  double k = q[RHO] * (u1*u1 + u2*u2 + u3*u3) / 2;  // bulk kinetic energy
  double pb = (sq(q[B1]) + sq(q[B2]) + sq(q[B3])) / 2; // magnetic pressure
  double p = (gas_gamma-1) * (q[ER] - k - pb); // plasma pressure

  flux[RHO] = q[M1];
  flux[M1] = u1*q[M1] - q[B1]*q[B1] + p + pb;
  flux[M2] = u1*q[M2] - q[B1]*q[B2];
  flux[M3] = u1*q[M3] - q[B1]*q[B3];
  flux[ER] = u1*(q[ER]+p+pb) - q[B1]*(u1*q[B1]+u2*q[B2]+u3*q[B3]);
  flux[B1] = 0;
  flux[B2] = u1*q[B2] - u2*q[B1];
  flux[B3] = u1*q[B3] - u3*q[B1];
}
