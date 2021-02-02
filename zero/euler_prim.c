#include <gkyl_euler_prim.h>

enum euler_cv { RHOU = 1, RHOV = 2, RHOW = 3 };

static const int dir_shuffle[][4] = {
  { 0, 1, 2, 3},
  { 0, 2, 3, 1},
  { 0, 3, 1, 2}
};

void
gkyl_euler_flux(int dir, double gas_gamma, const double *q, double *flux)
{
  const int *d = dir_shuffle[dir];
  
#define RHOU d[1]
#define RHOV d[2]
#define RHOW d[3]  
  
  double pr = gkyl_euler_pressure(gas_gamma, q), u = q[RHOU]/q[0];
  flux[0] = q[RHOU]; // rho*u
  flux[RHOU] = q[RHOU]*u + pr; // rho*u*u + pr
  flux[RHOV] = q[RHOV]*u; // rho*v*u
  flux[RHOW] = q[RHOW]*u; // rho*w*u
  flux[4] = (q[5]+pr)*u; // (E+p)*u
}
