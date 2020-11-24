#include <gkyl_util.h>

/**
 * Print error message to stderr and exit.
 *
 * @param msg Error message.
 */
void
gkyl_exit(const char* msg)
{
  fprintf(stderr, "Error: %s", msg);
  exit(EXIT_FAILURE);
}

/**
 * Compares two float numbers 'a' and 'b' to check if they are
 * sufficiently close by, where 'eps' is the relative tolerance.
 */
int
gkyl_compare_float(float a, float b, float eps)
{
  float absa = fabs(a), absb = fabs(b), diff = fabs(a-b);

  if (a == b) return 1;
  if (a == 0 || b == 0 || (absa+absb < FLT_MIN)) return diff < eps*FLT_MIN;
  return diff/fminf(absa+absb, FLT_MAX) < eps;
}

/**
 * Compares two double numbers 'a' and 'b' to check if they are
 * sufficiently close by, where 'eps' is the relative tolerance.
 */
int
gkyl_compare_double(double a, double b, double eps)
{
  double absa = fabs(a), absb = fabs(b), diff = fabs(a-b);

  if (a == b) return 1;
  if (a == 0 || b == 0 || (absa+absb < DBL_MIN)) return diff < eps*DBL_MIN;
  return diff/fmin(absa+absb, DBL_MAX) < eps;
}
