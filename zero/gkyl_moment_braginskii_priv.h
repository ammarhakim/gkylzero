#pragma once
#include <math.h>

// Private header, not for direct use in user code

// lower case l/u denote lower/upper with respect to edge

// Calculate symmetrized gradient 1D
// Based on Günter, Lackner, & Tichmann 2005 JCP
static inline double
calc_sym_grad_1D(double dx, double a_l, double a_u)
{
  return (a_u - a_l)/dx;
}

// Calculate symmetrized gradients 2D
static inline double
calc_sym_gradx_2D(double dx, double a_ll, double a_lu, double a_ul, double a_uu)
{
  return (a_ul + a_uu - a_ll - a_lu)/(2*dx);
}

static inline double
calc_sym_grady_2D(double dy, double a_ll, double a_lu, double a_ul, double a_uu)
{
  return (a_lu + a_uu - a_ll - a_ul)/(2*dy);
}

// Calculate symmetrized gradients 3D
static inline double
calc_sym_gradx_3D(double dx, double a_lll, double a_llu, double a_lul, double a_luu, double a_ull, double a_ulu, double a_uul, double a_uuu)
{
  return (a_ull + a_ulu + a_uul + a_uuu - a_lll - a_llu - a_lul - a_luu)/(4*dx);
}

static inline double
calc_sym_grady_3D(double dy, double a_lll, double a_llu, double a_lul, double a_luu, double a_ull, double a_ulu, double a_uul, double a_uuu)
{
  return (a_lul + a_luu + a_uul + a_uuu - a_lll - a_llu - a_ull - a_ulu)/(4*dy);
}

static inline double
calc_sym_gradz_3D(double dz, double a_lll, double a_llu, double a_lul, double a_luu, double a_ull, double a_ulu, double a_uul, double a_uuu)
{
  return (a_llu + a_luu + a_ulu + a_uuu - a_lll - a_lul - a_ull - a_uul)/(4*dz);
}

// In 1D, computes quantity at cell edge of two-cell interface
static inline double
calc_arithm_avg_1D(double a_l, double a_u)
{
  return 0.5*(a_l + a_u);
}

static inline double
calc_harmonic_avg_1D(double a_l, double a_u)
{
  return 1.0/(0.5/a_l + 0.5/a_u);
}

// In 2D, computes quantity at cell corner of four-cell interface
static inline double
calc_arithm_avg_2D(double a_ll, double a_lu, double a_ul, double a_uu)
{
  return 0.25*(a_ll + a_lu + a_ul + a_uu);
}

static inline double
calc_harmonic_avg_2D(double a_ll, double a_lu, double a_ul, double a_uu)
{
  return 1.0/(0.25/a_ll + 0.25/a_lu + 0.25/a_ul + 0.25/a_uu);
}

// In 3D, computes quantity at cell corner of eight-cell interface
static inline double
calc_arithm_avg_3D(double a_lll, double a_llu, double a_lul, double a_luu, double a_ull, double a_ulu, double a_uul, double a_uuu)
{
  return 0.125*(a_lll + a_llu + a_lul + a_luu + a_ull + a_ulu + a_uul + a_uuu);
}

static inline double
calc_harmonic_avg_3D(double a_lll, double a_llu, double a_lul, double a_luu, double a_ull, double a_ulu, double a_uul, double a_uuu)
{
  return 1.0/(0.125/a_lll + 0.125/a_llu + 0.125/a_lul + 0.125/a_luu + 0.125/a_ull + 0.125/a_ulu + 0.125/a_uul + 0.125/a_uuu);
}

// Calculate the collision time based on the species' parameters
// Note: assumes the electron-ion collision frequency so sqrt(2) may be missing
//       coulomb_log considered constant, rho is mass density, temp is temperature
static inline double
calc_tau(double coulomb_log, double epsilon0, double charge, double mass, double rho, double temp)
{
  double a = M_PI*mass*temp;
  return 6.0*sqrt(2.0*a*a*a)*epsilon0*epsilon0/(coulomb_log*charge*charge*charge*charge*rho);
}

// Calculate rate of strain tensor
// In 1D, computes a single tensor at cell edge of two-cell interface
static void
calc_ros_1D(double dx, double u_l[3], double u_u[3], double w[6])
{
  double gradx_ux = calc_sym_grad_1D(dx, u_l[0], u_u[0]);
  double gradx_uy = calc_sym_grad_1D(dx, u_l[1], u_u[1]);
  double gradx_uz = calc_sym_grad_1D(dx, u_l[2], u_u[2]);

  w[0] = 4.0/3.0*gradx_ux;
  w[1] = gradx_uy;
  w[2] = gradx_uz;
  w[3] = -2.0/3.0*gradx_ux;
  w[4] = 0.0;
  w[5] = -2.0/3.0*gradx_ux;
}

// Calculate rate of strain tensor
// In 2D, computes a tensor in one corner of four-cell interface (down left, down right, up left, & up right)
static void
calc_ros_2D(double dx, double dy, double u_ll[3], double u_lu[3], double u_ul[3], double u_uu[3], double w[6])
{
  double gradx_ux = calc_sym_gradx_2D(dx, u_ll[0], u_lu[0], u_ul[0], u_uu[0]);
  double gradx_uy = calc_sym_gradx_2D(dx, u_ll[1], u_lu[1], u_ul[1], u_uu[1]);
  double gradx_uz = calc_sym_gradx_2D(dx, u_ll[2], u_lu[2], u_ul[2], u_uu[2]);

  double grady_ux = calc_sym_grady_2D(dy, u_ll[0], u_lu[0], u_ul[0], u_uu[0]);
  double grady_uy = calc_sym_grady_2D(dy, u_ll[1], u_lu[1], u_ul[1], u_uu[1]);
  double grady_uz = calc_sym_grady_2D(dy, u_ll[2], u_lu[2], u_ul[2], u_uu[2]);

  double divu = gradx_ux + grady_uy;
  w[0] = 2.0*gradx_ux - 2.0/3.0*divu;
  w[1] = gradx_uy + grady_ux;
  w[2] = gradx_uz;
  w[3] = 2.0*grady_uy - 2.0/3.0*divu;
  w[4] = grady_uz;
  w[5] = -2.0/3.0*divu;
}
