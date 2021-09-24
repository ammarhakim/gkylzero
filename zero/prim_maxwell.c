#include <math.h>

#include <gkyl_prim_maxwell.h>

// Make indexing cleaner with the dir_shuffle
#define EX 0
#define EY 1
#define EZ 2
#define BX 3
#define BY 4
#define BZ 5

double
gkyl_maxwell_max_abs_speed(double c, double e_fact, double b_fact, const double q[8])
{
  return c;
}

void
gkyl_maxwell_flux(double c, double e_fact, double b_fact, const double q[8], double flux[8])
{
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
