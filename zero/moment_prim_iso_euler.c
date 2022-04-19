#include <math.h>

#include <gkyl_moment_prim_iso_euler.h>

#define RHOU 1
#define RHOV 2
#define RHOW 3

double
gkyl_iso_euler_max_abs_speed(double vt, const double q[4])
{
  double u = q[RHOU]/q[0];
  return fabs(u) + vt;
}

void
gkyl_iso_euler_flux(double vt, const double q[4], double flux[4])
{
  double u = q[RHOU]/q[0];
  flux[0] = q[RHOU]; // rho*u
  flux[RHOU] = q[RHOU]*u + q[0]*vt*vt; // rho*(u*u + vt*vt)
  flux[RHOV] = q[RHOV]*u; // rho*v*u
  flux[RHOW] = q[RHOW]*u; // rho*w*u
}
