#include <math.h>

#include <gkyl_prim_euler.h>

static const int dir_shuffle[][3] = {
  {1, 2, 3},
  {2, 3, 1},
  {3, 1, 2}
};

#define RHOU d[0]
#define RHOV d[1]
#define RHOW d[2]

double
gkyl_euler_max_abs_speed(int dir, double gas_gamma, const double q[5])
{
  int const *const d = dir_shuffle[dir];
  double pr = gkyl_euler_pressure(gas_gamma, q), u = q[RHOU]/q[0];
  return fabs(u) + sqrt(gas_gamma*pr/q[0]);
}

void
gkyl_euler_flux(int dir, double gas_gamma, const double q[5], double flux[5])
{
  int const *const d = dir_shuffle[dir];
  double pr = gkyl_euler_pressure(gas_gamma, q), u = q[RHOU]/q[0];
  flux[0] = q[RHOU]; // rho*u
  flux[RHOU] = q[RHOU]*u + pr; // rho*u*u + pr
  flux[RHOV] = q[RHOV]*u; // rho*v*u
  flux[RHOW] = q[RHOW]*u; // rho*w*u
  flux[4] = (q[4]+pr)*u; // (E+p)*u
}
