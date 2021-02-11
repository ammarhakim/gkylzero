#include <gkyl_euler_prim.h>

#include <math.h>

static const int dir_shuffle[][4] = {
  { 0, 1, 2, 3},
  { 0, 2, 3, 1},
  { 0, 3, 1, 2}
};

#define RHOU d[1]
#define RHOV d[2]
#define RHOW d[3]

void
gkyl_euler_flux(int dir, double gas_gamma, const double q[5], double flux[5])
{
  const int *d = dir_shuffle[dir];
  double pr = gkyl_euler_pressure(gas_gamma, q), u = q[RHOU]/q[0];
  flux[0] = q[RHOU]; // rho*u
  flux[RHOU] = q[RHOU]*u + pr; // rho*u*u + pr
  flux[RHOV] = q[RHOV]*u; // rho*v*u
  flux[RHOW] = q[RHOW]*u; // rho*w*u
  flux[4] = (q[4]+pr)*u; // (E+p)*u
}

double
gkyl_euler_max_abs_speed(int dir, double gas_gamma, const double q[5])
{
  const int *d = dir_shuffle[dir];
  double pr = gkyl_euler_pressure(gas_gamma, q), u = q[RHOU]/q[0];
  return fabs(u) + sqrt(gas_gamma*pr/q[0]);
}
