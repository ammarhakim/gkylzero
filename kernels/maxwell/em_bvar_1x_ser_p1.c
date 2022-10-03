#include <gkyl_maxwell_kernels.h> 
#include <gkyl_basis_ser_1x_p1_exp_sq.h> 
#include <gkyl_basis_ser_1x_p1_inv.h> 
#include <gkyl_basis_ser_1x_p1_sqrt.h> 
GKYL_CU_DH void em_bvar_1x_ser_p1(const double *em, double* GKYL_RESTRICT bvar) 
{ 
  // em:   Input electromagnetic fields. 
  // bvar: b_i = B_i/|B| (first 3 components), b_i b_j = B_i B_j/|B|^2 (last 6 components). 
 
  const double *B_x = &em[6]; 
  const double *B_y = &em[8]; 
  const double *B_z = &em[10]; 
 
  double *bx = &bvar[0]; 
  double *by = &bvar[2]; 
  double *bz = &bvar[4]; 
  double *bxbx = &bvar[6]; 
  double *bxby = &bvar[8]; 
  double *bxbz = &bvar[10]; 
  double *byby = &bvar[12]; 
  double *bybz = &bvar[14]; 
  double *bzbz = &bvar[16]; 
 
  // Calculate |B|^2 and get expansion of 1/|B|^2. 
  double B_x_sq[2] = {0.0}; 
  ser_1x_p1_exp_sq(B_x, B_x_sq); 
 
  double B_y_sq[2] = {0.0}; 
  ser_1x_p1_exp_sq(B_y, B_y_sq); 
 
  double B_z_sq[2] = {0.0}; 
  ser_1x_p1_exp_sq(B_z, B_z_sq); 
 
  double magB_sq[2] = {0.0}; 

  magB_sq[0] = B_x_sq[0] + B_y_sq[0] + B_z_sq[0]; 
  magB_sq[1] = B_x_sq[1] + B_y_sq[1] + B_z_sq[1]; 
 
  double magB_sq_inv[2] = {0.0}; 

  ser_1x_p1_inv(magB_sq, magB_sq_inv); 
  // Calculate expansions of B_i B_j/|B|^2, which can be calculated free of aliasing errors. 
  bxbx[0] = 0.7071067811865475*B_x_sq[1]*magB_sq_inv[1]+0.7071067811865475*B_x_sq[0]*magB_sq_inv[0]; 
  bxbx[1] = 0.7071067811865475*B_x_sq[0]*magB_sq_inv[1]+0.7071067811865475*magB_sq_inv[0]*B_x_sq[1]; 

  bxby[0] = 0.5*B_x[0]*B_y[1]*magB_sq_inv[1]+0.5*B_y[0]*B_x[1]*magB_sq_inv[1]+0.5*magB_sq_inv[0]*B_x[1]*B_y[1]+0.5*B_x[0]*B_y[0]*magB_sq_inv[0]; 
  bxby[1] = 0.9*B_x[1]*B_y[1]*magB_sq_inv[1]+0.5*B_x[0]*B_y[0]*magB_sq_inv[1]+0.5*B_x[0]*magB_sq_inv[0]*B_y[1]+0.5*B_y[0]*magB_sq_inv[0]*B_x[1]; 

  bxbz[0] = 0.5*B_x[0]*B_z[1]*magB_sq_inv[1]+0.5*B_z[0]*B_x[1]*magB_sq_inv[1]+0.5*magB_sq_inv[0]*B_x[1]*B_z[1]+0.5*B_x[0]*B_z[0]*magB_sq_inv[0]; 
  bxbz[1] = 0.9*B_x[1]*B_z[1]*magB_sq_inv[1]+0.5*B_x[0]*B_z[0]*magB_sq_inv[1]+0.5*B_x[0]*magB_sq_inv[0]*B_z[1]+0.5*B_z[0]*magB_sq_inv[0]*B_x[1]; 

  byby[0] = 0.7071067811865475*B_y_sq[1]*magB_sq_inv[1]+0.7071067811865475*B_y_sq[0]*magB_sq_inv[0]; 
  byby[1] = 0.7071067811865475*B_y_sq[0]*magB_sq_inv[1]+0.7071067811865475*magB_sq_inv[0]*B_y_sq[1]; 

  bybz[0] = 0.5*B_y[0]*B_z[1]*magB_sq_inv[1]+0.5*B_z[0]*B_y[1]*magB_sq_inv[1]+0.5*magB_sq_inv[0]*B_y[1]*B_z[1]+0.5*B_y[0]*B_z[0]*magB_sq_inv[0]; 
  bybz[1] = 0.9*B_y[1]*B_z[1]*magB_sq_inv[1]+0.5*B_y[0]*B_z[0]*magB_sq_inv[1]+0.5*B_y[0]*magB_sq_inv[0]*B_z[1]+0.5*B_z[0]*magB_sq_inv[0]*B_y[1]; 

  bzbz[0] = 0.7071067811865475*B_z_sq[1]*magB_sq_inv[1]+0.7071067811865475*B_z_sq[0]*magB_sq_inv[0]; 
  bzbz[1] = 0.7071067811865475*B_z_sq[0]*magB_sq_inv[1]+0.7071067811865475*magB_sq_inv[0]*B_z_sq[1]; 

  // Calculate b_i = B_i/|B| by taking square root of B_i^2/|B|^2 at quadrature points. 
  ser_1x_p1_sqrt(bxbx, bx); 
  ser_1x_p1_sqrt(byby, by); 
  ser_1x_p1_sqrt(bzbz, bz); 
} 
 
