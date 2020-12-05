#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gkyl_util.h>

void
gkyl_exit(const char* msg)
{
  fprintf(stderr, "Error: %s", msg);
  exit(EXIT_FAILURE);
}

int
gkyl_compare_flt(float a, float b, float eps)
{
  float absa = fabs(a), absb = fabs(b), diff = fabs(a-b);

  if (a == b) return 1;
  if (a == 0 || b == 0 || (absa+absb < FLT_MIN)) return diff < eps*FLT_MIN;
  return diff/fminf(absa+absb, FLT_MAX) < eps;
}

int
gkyl_compare_dbl(double a, double b, double eps)
{
  double absa = fabs(a), absb = fabs(b), diff = fabs(a-b);

  if (a == b) return 1;
  if (a == 0 || b == 0 || (absa+absb < DBL_MIN)) return diff < eps*DBL_MIN;
  return diff/fmin(absa+absb, DBL_MAX) < eps;
}
