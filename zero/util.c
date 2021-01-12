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

int
gkyl_wall_clock(struct timespec *ts)
{
  return clock_gettime(CLOCK_REALTIME, ts);
}

struct timespec
gkyl_time_diff(struct timespec *start, struct timespec *end)
{
  struct timespec tm;
  if ((end->tv_nsec-start->tv_nsec)<0) {
    tm.tv_sec = end->tv_sec-start->tv_sec-1;
    tm.tv_nsec = 1000000000+end->tv_nsec-start->tv_nsec;
  }
  else {
    tm.tv_sec = end->tv_sec-start->tv_sec;
    tm.tv_nsec = end->tv_nsec-start->tv_nsec;
  }
  return tm;  
}
   
double
gkyl_time_sec(struct timespec tm)
{
  return tm.tv_sec + 1e-9*tm.tv_nsec;
}
