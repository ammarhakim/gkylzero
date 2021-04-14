#include <math.h>

#include <gkyl_maxwell_prim.h>

static const int dir_shuffle[][6] = {
  {0, 1, 2, 3, 4, 5},
  {1, 2, 0, 4, 5, 3},
  {2, 0, 1, 5, 3, 4}
};

// Make indexing cleaner with the dir_shuffle
#define EX d[0]
#define EY d[1]
#define EZ d[2]
#define BX d[3]
#define BY d[4]
#define BZ d[5]

double
gkyl_maxwell_max_abs_speed(int dir, double c, double e_fact, double b_fact, const double q[8])
{
  return c;
}

void
gkyl_maxwell_flux(int dir, double c, double e_fact, double b_fact, const double q[8], double flux[8])
{
  const int *d = dir_shuffle[dir];
  double c2 = c*c;

  flux[EX] = e_fact*c2*q[6]; // e_fact*c^2*phi
  flux[EY] = c2*q[BZ]; // c^2*Bz
  flux[EZ] = -c2*q[BY]; // -c^2*By
  flux[BX] = b_fact*q[7]; // b_fact*psi
  flux[BY] = -q[EZ]; // -Ez
  flux[BZ] = q[EY]; // Ey
  flux[6] = e_fact*q[EX]; // e_fact*Ex
  flux[7] = b_fact*c2*q[BX]; // b_fact*c^2*Bx
}
