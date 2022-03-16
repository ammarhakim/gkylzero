#include <math.h>

#include <gkyl_moment_prim_euler.h>

#define RHOU 1
#define RHOV 2
#define RHOW 3

double
gkyl_euler_max_abs_speed(double gas_gamma, const double q[5])
{
  double pr = gkyl_euler_pressure(gas_gamma, q), u = q[RHOU]/q[0];
  return fabs(u) + sqrt(gas_gamma*pr/q[0]);
}

void
gkyl_euler_flux(double gas_gamma, const double q[5], double flux[5])
{
  double pr = gkyl_euler_pressure(gas_gamma, q), u = q[RHOU]/q[0];
  flux[0] = q[RHOU]; // rho*u
  flux[RHOU] = q[RHOU]*u + pr; // rho*u*u + pr
  flux[RHOV] = q[RHOV]*u; // rho*v*u
  flux[RHOW] = q[RHOW]*u; // rho*w*u
  flux[4] = (q[4]+pr)*u; // (E+p)*u
}
