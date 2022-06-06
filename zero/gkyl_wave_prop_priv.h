#pragma once

#include <assert.h>
#include <float.h>
#include <math.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_util.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wave_geom_priv.h>
#include <gkyl_wave_prop.h>

GKYL_CU_DH
static inline double
fmax3(double a, double b, double c)
{
  return fmax(fmax(a,b),c);
}

GKYL_CU_DH
static inline double
fmin3(double a, double b, double c)
{
  return fmin(fmin(a,b),c);
}

// limiter function
GKYL_CU_DH
static inline double
limiter_function(double r, enum gkyl_wave_limiter limiter)
{
  double theta = 0.0;
  switch (limiter) {
    case GKYL_NO_LIMITER:
      theta = 1.0;
      break;
    
    case GKYL_MIN_MOD:
      theta = fmax(0., fmin(1., r));
      break;

    case GKYL_SUPERBEE:
      theta = fmax3(0.0, fmin(1., 2*r), fmin(2.0, r));
      break;

    case GKYL_VAN_LEER:
      theta = (r+fabs(r))/(1+fabs(r));
      break;

    case GKYL_MONOTONIZED_CENTERED:
      theta = fmax(0.0, fmin3((1.0+r)/2, 2, 2*r));
      break;

    case GKYL_BEAM_WARMING:
      theta = r;
      break;

    case GKYL_ZERO:
      theta = 0;
      break;
  }
  return theta;
}

GKYL_CU_DH
static inline void
calc_jump(int n, const double *ql, const double *qr, double * GKYL_RESTRICT jump)
{
  for (int d=0; d<n; ++d) jump[d] = qr[d]-ql[d];
}

GKYL_CU_DH
static inline void
calc_first_order_update(int meqn, double dtdx,
  double * GKYL_RESTRICT ql, double * GKYL_RESTRICT qr, const double *amdq, const double *apdq)
{
  for (int i=0; i<meqn; ++i) {
    qr[i] = qr[i] - dtdx*apdq[i];
    ql[i] = ql[i] - dtdx*amdq[i];
  }
}

GKYL_CU_DH
static inline double
calc_cfla(int mwaves, double cfla, double dtdx, const double *s)
{
  double c = cfla;
  for (int i=0; i<mwaves; ++i)
    c = fmax(c, dtdx*fabs(s[i]));
  return c;
}

GKYL_CU_DH
static inline double
wave_dot_prod(int meqn, const double * GKYL_RESTRICT wa, const double * GKYL_RESTRICT wb)
{
  double dot = 0.0;
  for (int i=0; i<meqn; ++i) dot += wa[i]*wb[i];
  return dot;
}

GKYL_CU_DH
static inline void
wave_rescale(int meqn, double fact, double *w)
{
  for (int i=0; i<meqn; ++i) w[i] *= fact; 
}

GKYL_CU_DH
static inline void
calc_second_order_flux(int meqn, double dtdx, double s,
  const double *waves, double * GKYL_RESTRICT flux2)
{
  double sfact = 0.5*fabs(s)*(1-fabs(s)*dtdx);
  for (int i=0; i<meqn; ++i)
    flux2[i] += sfact*waves[i];
}

GKYL_CU_DH
static inline void
calc_second_order_update(int meqn, double dtdx, double * GKYL_RESTRICT qout,
  const double *fl, const double *fr)
{
  for (int i=0; i<meqn; ++i)
    qout[i] += -dtdx*(fr[i]-fl[i]);
}

