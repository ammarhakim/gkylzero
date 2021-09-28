#include <math.h>

#include <gkyl_prim_ten_moment.h>

// Make indexing cleaner with the dir_shuffle
#define RHOU 1
#define RHOV 2
#define RHOW 3

#define PXX 4
#define PXY 5
#define PXZ 6
#define PYY 7
#define PYZ 8
#define PZZ 9

double
gkyl_ten_moment_max_abs_speed(const double q[10])
{
  double u = q[RHOU]/q[0];
  double p11 = q[PXX] - q[0]*u*u;
  return fabs(u) + sqrt(3.0*p11/q[0]);
}

void
gkyl_ten_moment_flux(const double q[10], double flux[10])
{
  double v[10];
  gkyl_ten_moment_primitive(q, v);

  flux[0] = q[RHOU]; // rho*u
  flux[RHOU] = q[PXX]; // Pxx
  flux[RHOV] = q[PXY]; // Pxy
  flux[RHOW] = q[PXZ]; // Pxz
  flux[PXX] = v[0]*v[RHOU]*v[RHOU]*v[RHOU] + 3*v[RHOU]*v[PXX]; // rho u^3 + 3*u*Pxx
  flux[PXY] = v[0]*v[RHOU]*v[RHOU]*v[RHOV] + 2*v[RHOU]*v[PXY] + v[RHOV]*v[PXX]; // rho*u^2*v + 2*u*Pxy + v*Pxx
  flux[PXZ] = v[0]*v[RHOU]*v[RHOU]*v[RHOW] + 2*v[RHOU]*v[PXZ] + v[RHOW]*v[PXX]; // rho*u^2*w + 2*u*Pxz + w*Pxx
  flux[PYY] = v[0]*v[RHOU]*v[RHOV]*v[RHOV] + 2*v[RHOV]*v[PXY] + v[RHOU]*v[PYY]; // rho*u*v^2 + 2*v*Pxy + u*Pyy
  flux[PYZ] = v[0]*v[RHOU]*v[RHOV]*v[RHOW] + v[RHOU]*v[PYZ] + v[RHOV]*v[PXZ] + v[RHOW]*v[PXY]; // rho*u*v*w + u*Pyz + v*Pxz + w*Pxy
  flux[PZZ] = v[0]*v[RHOU]*v[RHOW]*v[RHOW] + 2*v[RHOW]*v[PXZ] + v[RHOU]*v[PZZ]; // rho*u*w^2 + 2*w*Pxz + u*Pzz
}
