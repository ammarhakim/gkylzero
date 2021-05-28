#include <math.h>

#include <gkyl_prim_iso_euler.h>

static const int dir_shuffle[][3] = {
  {1, 2, 3},
  {2, 3, 1},
  {3, 1, 2}
};

#define RHOU d[0]
#define RHOV d[1]
#define RHOW d[2]

double
gkyl_iso_euler_max_abs_speed(int dir, double vt, const double q[4])
{
  int const *const d = dir_shuffle[dir];
  double u = q[RHOU]/q[0];
  return fabs(u) + vt;
}

void
gkyl_iso_euler_flux(int dir, double vt, const double q[4], double flux[4])
{
  int const *const d = dir_shuffle[dir];
  double u = q[RHOU]/q[0];
  flux[0] = q[RHOU]; // rho*u
  flux[RHOU] = q[RHOU]*u + q[0]*vt*vt; // rho*(u*u + vt*vt)
  flux[RHOV] = q[RHOV]*u; // rho*v*u
  flux[RHOW] = q[RHOW]*u; // rho*w*u
}
