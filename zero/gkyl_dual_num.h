#pragma once

#include <math.h>
#include <gkyl_util.h>

// Dual or hyperreal number
struct gkyl_dn { double x[2]; };

// Construct new dual numbers
GKYL_CU_DH
static inline struct gkyl_dn
gdn_new(double x0, double x1)
{
  return (struct gkyl_dn) { x0, x1 };
}

GKYL_CU_DH
static inline struct gkyl_dn
gdn_new0(double x0)
{
  return (struct gkyl_dn) { x0, 0.0 };
}

GKYL_CU_DH
static inline struct gkyl_dn
gdn_new1(double x0)
{
  return (struct gkyl_dn) { x0, 1.0 };
}

// Basic operators:
// 
// neg, add, sadd, sub, ssub, mul, smul, div, sdiv

GKYL_CU_DH
static inline struct gkyl_dn
gdn_neg(struct gkyl_dn d1)
{
  return (struct gkyl_dn) { -d1.x[0], -d1.x[1] };
}

GKYL_CU_DH
static inline struct gkyl_dn
gdn_add(struct gkyl_dn d1, struct gkyl_dn d2)
{
  return (struct gkyl_dn) { d1.x[0]+d2.x[0], d1.x[1]+d2.x[1] };
}

GKYL_CU_DH
static inline struct gkyl_dn
gdn_sadd(double s, struct gkyl_dn d1)
{
  return (struct gkyl_dn) { s+d1.x[0], d1.x[1] };
}

GKYL_CU_DH
static inline struct gkyl_dn
gdn_sub(struct gkyl_dn d1, struct gkyl_dn d2)
{
  return (struct gkyl_dn) { d1.x[0]-d2.x[0], d1.x[1]-d2.x[1] };
}

GKYL_CU_DH
static inline struct gkyl_dn
gdn_ssub(double s, struct gkyl_dn d1)
{
  return (struct gkyl_dn) { s-d1.x[0], -d1.x[1] };
}

GKYL_CU_DH
static inline struct gkyl_dn
gdn_mul(struct gkyl_dn d1, struct gkyl_dn d2)
{
  return (struct gkyl_dn) {
      d1.x[0]*d2.x[0],
      d1.x[0]*d2.x[1] + d1.x[1]*d2.x[0]    
  };
}

GKYL_CU_DH
static inline struct gkyl_dn
gdn_smul(double s, struct gkyl_dn d1)
{
  return (struct gkyl_dn) { s*d1.x[0], s*d1.x[1] };
}

GKYL_CU_DH
static inline struct gkyl_dn
gdn_inv(struct gkyl_dn d1)
{
  return (struct gkyl_dn) {
    1.0/d1.x[0],
    -d1.x[1]/(d1.x[0]*d1.x[0])
  };
}

GKYL_CU_DH
static inline struct gkyl_dn
gdn_div(struct gkyl_dn d1, struct gkyl_dn d2)
{
  return gdn_mul(d1, gdn_inv(d2));
}

GKYL_CU_DH
static inline struct gkyl_dn
gdn_sdiv(double s, struct gkyl_dn d1)
{
  return gdn_smul(s, gdn_inv(d1));
}

// Functions
// 
// sq, cube, npow, sqrt

GKYL_CU_DH
static inline struct gkyl_dn
gdn_sq(struct gkyl_dn d1)
{
  return (struct gkyl_dn) { d1.x[0]*d1.x[0], 2.0*d1.x[0]*d1.x[1] };
}

GKYL_CU_DH
static inline struct gkyl_dn
gdn_cube(struct gkyl_dn d1)
{
  return gdn_mul(d1, gdn_sq(d1));
}

GKYL_CU_DH
static inline struct gkyl_dn
gdn_npow(struct gkyl_dn d1, int n)
{
  double pn1 = pow(d1.x[0], n-1);
  return (struct gkyl_dn) { pn1*d1.x[0], n*pn1*d1.x[1] };
}

GKYL_CU_DH
static inline struct gkyl_dn
gdn_sqrt(struct gkyl_dn d1)
{
  double sx = sqrt(d1.x[0]);
  return (struct gkyl_dn) { sx, d1.x[1]*0.5/sx  };
}
