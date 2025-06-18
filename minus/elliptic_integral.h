#pragma once

#include <math.h>       // required for fabs(), fabsl(), sqrtl(), and M_PI_2
#include <float.h>      // required for LDBL_EPSILON, DBL_MAX

double Complete_Elliptic_Integral_First_Kind(char arg, double x);
double Complete_Elliptic_Integral_Second_Kind(char arg, double x);