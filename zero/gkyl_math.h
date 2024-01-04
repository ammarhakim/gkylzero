#pragma once

#include <math.h>

// Three-component vector
struct  gkyl_vec3 { double x[3]; };

// new vector with zeros 
static inline struct gkyl_vec3
gkyl_vec3_zeros()
{
  return (struct gkyl_vec3) { .x = { 0.0, 0.0, 0.0} };
}

// new vector with x,y,z components
static inline struct gkyl_vec3
gkyl_vec3_new(double x, double y, double z)
{
  return (struct gkyl_vec3) { .x = { x, y, z} };
}

// a+b
static inline struct gkyl_vec3
gkyl_vec3_add(struct gkyl_vec3 a, struct gkyl_vec3 b)
{
  return (struct gkyl_vec3) { .x = { a.x[0]+b.x[0], a.x[1]+b.x[1], a.x[2]+b.x[2] } };
}

// a-b
static inline struct gkyl_vec3
gkyl_vec3_sub(struct gkyl_vec3 a, struct gkyl_vec3 b)
{
  return (struct gkyl_vec3) { .x = { a.x[0]-b.x[0], a.x[1]-b.x[1], a.x[2]-b.x[2] } };
}

// |a|
static inline double
gkyl_vec3_len(struct gkyl_vec3 a)
{
  return sqrt( a.x[0]*a.x[0] + a.x[1]*a.x[1] + a.x[2]*a.x[2] );
}

// normalize a
static inline struct gkyl_vec3
gkyl_vec3_norm(struct gkyl_vec3 a)
{
  double len = gkyl_vec3_len(a);
  return (struct gkyl_vec3) { .x = { a.x[0]/len, a.x[1]/len, a.x[2]/len } };
}

// a \dot b
static inline double
gkyl_vec3_dot(struct gkyl_vec3 a, struct gkyl_vec3 b)
{
  return a.x[0]*b.x[0] + a.x[1]*b.x[1] + a.x[2]*b.x[2];
}

// a \times b
static inline struct gkyl_vec3
gkyl_vec3_cross(struct gkyl_vec3 a, struct gkyl_vec3 b)
{
  return (struct gkyl_vec3) { .x = {
      a.x[1]*b.x[2]-a.x[2]*b.x[1],
      a.x[2]*b.x[0]-a.x[0]*b.x[2],
      a.x[0]*b.x[1]-a.x[1]*b.x[0]
    }
  };
}

// minmod(a,b,c, ...) returns min of all parameters are positive, max
// of all parameters if all are negative, and zero otherwise.

// minmod(x,y)
static inline double
gkyl_minmod_2(double x, double y)
{
  if (x>0 && y>0)
    return fmin(x,y);
  if (x<0 && y<0)
    return fmax(x,y);
  return 0.0;
}

// minmod(x,y)
static inline double
gkyl_minmod_3(double x, double y, double z)
{
  if (x>0 && y>0 && z>0)
    return fmin(fmin(x,y),z);
  if (x<0 && y<0 && z<0)
    return fmax(fmax(x,y),z);
  return 0.0;
}

// minmod(x,y,z,w)
static inline double
gkyl_minmod_4(double x, double y, double z, double w)
{
  if (x>0 && y>0 && z>0 && w>0)
    return fmin(fmin(x,y),fmin(z,w));
  if (x<0 && y<0 && z<0 && w<0)
    return fmax(fmax(x,y),fmax(z,w));
  return 0.0;
}

// median(x,y,z) returns the value between two others
static inline double
gkyl_median(double x, double y, double z)
{
  return x + gkyl_minmod_2(y-x,z-x);
}

// min(x,y,z)
static inline double
gkyl_min_3(double x, double y, double z)
{
  return fmin(x,fmin(y,z));
}

// max(x,y,z)
static inline double
gkyl_max_3(double x, double y, double z)
{
  return fmax(x,fmax(y,z));
}

// Result from numerical quadrature
struct gkyl_qr_res {
  double res; // result of quadrature
  double error; // error estimate
  int status; // 0 - success.
  int nevals; // number of function evaluations
  int nlevels; // number of levels
};

// polynomial order for low-order polynomial
enum gkyl_lo_poly_order {
  GKYL_LO_POLY_2 = 2,
  GKYL_LO_POLY_3,
  GKYL_LO_POLY_4,  
};  

// Polynomial roots for low-order polynomials (upto 4th order)
struct gkyl_lo_poly_roots {
  double rpart[4], impart[4]; // real and imaginary part of roots
  double err[4]; // estimated error for each root
  int niter; // number of iterations
};

// Polynomial roots for arbitrary order polynomials: allocate an free
// using new and release methods as usual
struct gkyl_poly_roots {
  int poly_order; // polynomial order
  double *rpart, *impart; // real and imaginary part of roots
  double *err; // estimated error for each root
  int niter;   // number of iterations
  void *work; // some memory needed internally. Do not muck or access!
};  

/**
 * Integrate func(x,ctx) in interval [a,b], using maximum levels
 * "n". Typical value of levels is 7. The status flag of the return
 * struct is set to 1 if quadrature failed due to number of levels
 * exceeded.
 *
 * @param func Function to integrate
 * @param ctx Context to pass to function
 * @param a Lower limit of interval
 * @param b Upper limit of interval
 * @param n Number of levels
 * @param eps Error to compute integral
 * @return Result of quadrature
 */
struct gkyl_qr_res gkyl_dbl_exp(double (*func)(double, void *), void *ctx,
  double a, double b, int n, double eps);

/**
 * Compute single real root of a function using Ridders' method. See
 * IEEE Tran. Circuit and Systems, vol CAS-26 No 11, Pg 976 1976. The
 * root must be inside the interval specified, and there must only be
 * * one such root in the interval.
 *
 * @param func Function to integrate
 * @param ctx Context to pass to function
 * @param x1 Lower limit of interval
 * @param x2 Upper limit of interval
 * @param f1 Function at x1
 * @param f2 Function at x2
 * @param max_iter Maximum number of iterations
 * @param eps Desired error bound
 * @return Result of root finder
 */
struct gkyl_qr_res gkyl_ridders(double (*func)(double,void*), void *ctx,
  double x1, double x2, double f1, double f2, int max_iter, double eps);

/**
 * Compute all simple roots of the low-order polynomial with monomial
 * coefficients given by @a coeff. The leading coefficient is always *
 * assumed be 1.0 and so the coeff[i] give the coefficient for the *
 * monomial x^i. For example for the quartic:
 *
 * p(x) = x^4 + 2 x^3 + 4
 *
 * coeff[4] = { 4.0, 0.0, 0.0, 2.0 };
 *
 * @param order Order of polynomial
 * @param coeff Coefficients of the polynomial
 * @return Roots of the polynomial
 */
struct gkyl_lo_poly_roots gkyl_calc_lo_poly_roots(enum gkyl_lo_poly_order order,
  double coeff[4]);

/**
 * Allocate memory to store polynomial roots.
 *
 * @param poly_order Polynomial order
 * @return newly allocated memory. Free using release method.
 */
struct gkyl_poly_roots* gkyl_poly_roots_new(int poly_order);

/**
 * Compute all simple roots of polymomial with real coefficients. The
 * leading coefficient must be 1.0 and the coeff[i] is the coefficient
 * of x^i.
 *
 * @param pr On output, contains all the roots. Must be preallocated
 * @param coeff Coefficients of monomials
 */
void gkyl_calc_poly_roots(struct gkyl_poly_roots *pr, const double *coeff);

/**
 * Release memory for use in polynomial root finder.
 *
 * @param pr Memory to release.
 */
void gkyl_poly_roots_release(struct gkyl_poly_roots *pr);
