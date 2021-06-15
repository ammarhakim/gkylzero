#pragma once

// Private header, not for direct use in user code

// lower case l/u denote lower/upper with respect to edge

// Calculate symmetrized gradient 1D
// Based on Günter, Lackner, & Tichmann 2005 JCP
static void
calc_sym_grad_1D(double dx, double a_l, double a_u, double gradx_a)
{
  gradx_a = (a_u - a_l)/dx;
}

// Calculate symmetrized gradient 2D
static void
calc_sym_grad_2D(double dx, double dy, double a_ll, double a_lu, double a_ul, double a_uu, double gradx_a, double grady_a)
{
  gradx_a = (a_ul + a_uu - a_ll - a_lu)/(2*dx);
  grady_a = (a_lu + a_uu - a_ll - a_ul)/(2*dy);
}

// Calculate rate of strain tensor
// In 1D, computes a single tensor at cell edge of two-cell interface
static void
calc_ros_1D(double dx, double u_l[3], double u_u[3], double w[6])
{
  double gradx_ux, gradx_uy, gradx_uz;
  calc_sym_grad_1D(dx, u_l[0], u_u[0], gradx_ux);
  calc_sym_grad_1D(dx, u_l[1], u_u[1], gradx_uy);
  calc_sym_grad_1D(dx, u_l[2], u_u[2], gradx_uz);

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
  double gradx_ux, gradx_uy, gradx_uz;
  double grady_ux, grady_uy, grady_uz;
  calc_sym_grad_2D(dx, dy, u_ll[0], u_lu[0], u_ul[0], u_uu[0], gradx_ux, grady_ux);
  calc_sym_grad_2D(dx, dy, u_ll[1], u_lu[1], u_ul[1], u_uu[1], gradx_uy, grady_uy);
  calc_sym_grad_2D(dx, dy, u_ll[2], u_lu[2], u_ul[2], u_uu[2], gradx_uz, grady_uz);

  double divu = gradx_ux + grady_uy;
  w[0] = 2.0*gradx_ux - 2.0/3.0*divu;
  w[1] = gradx_uy + grady_ux;
  w[2] = gradx_uz;
  w[3] = 2.0*grady_uy - 2.0/3.0*divu;
  w[4] = grady_uz;
  w[5] = -2.0/3.0*divu;
}

// Calculate arithmetic average
// In 1D, computes quantity at cell edge of two-cell interface
static void
calc_arithm_avg_1D(double a_l, double a_u, double avg)
{
  avg = 0.5*(a_l + a_u);
}

// In 2D, computes quantity at cell corner of four-cell interface (down left, down right, up left, & up right)
static void
calc_arithm_avg_2D(double a_ll, double a_lu, double a_ul, double a_uu, double avg)
{
  avg = 0.25*(a_ll + a_lu + a_ul + a_uu);
}

// Calculate harmonic average
// In 1D, computes quantity at cell edge of two-cell interface
static void
calc_harmonic_avg_1D(double a_l, double a_u, double avg)
{
  avg = 1.0/(0.5/a_l + 0.5/a_u);
}

// In 2D, computes quantity at cell corner of four-cell interface (down left, down right, up left, & up right)
static void
calc_harmonic_avg_2D(double a_ll, double a_lu, double a_ul, double a_uu, double avg)
{
  avg = 1.0/(0.25/a_ll + 0.25/a_lu + 0.25/a_ul + 0.25/a_uu);
}