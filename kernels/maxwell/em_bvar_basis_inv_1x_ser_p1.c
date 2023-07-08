#include <gkyl_maxwell_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
#include <gkyl_basis_ser_1x_p1_inv.h> 
#include <gkyl_basis_ser_1x_p1_sqrt_with_sign.h> 
GKYL_CU_DH void em_bvar_basis_inv_1x_ser_p1(const double *em, double* bvar) 
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
 
  // Calculate B_i B_j. 
  double B_x_sq[2] = {0.0}; 
  binop_mul_1d_ser_p1(B_x, B_x, B_x_sq); 
 
  double B_y_sq[2] = {0.0}; 
  binop_mul_1d_ser_p1(B_y, B_y, B_y_sq); 
 
  double B_z_sq[2] = {0.0}; 
  binop_mul_1d_ser_p1(B_z, B_z, B_z_sq); 
 
  double B_x_B_y[2] = {0.0}; 
  binop_mul_1d_ser_p1(B_x, B_y, B_x_B_y); 
 
  double B_x_B_z[2] = {0.0}; 
  binop_mul_1d_ser_p1(B_x, B_z, B_x_B_z); 
 
  double B_y_B_z[2] = {0.0}; 
  binop_mul_1d_ser_p1(B_y, B_z, B_y_B_z); 
 
  double magB2[2] = {0.0}; 

  magB2[0] = B_z_sq[0]+B_y_sq[0]+B_x_sq[0]; 
  magB2[1] = B_z_sq[1]+B_y_sq[1]+B_x_sq[1]; 

  int cell_avg = 0;
  // Check if |B|^2 < 0 at control points. 
  if (0.7071067811865475*magB2[0]-1.224744871391589*magB2[1] < 0.0) cell_avg = 1; 
  if (1.224744871391589*magB2[1]+0.7071067811865475*magB2[0] < 0.0) cell_avg = 1; 
  double magB2_inv[2] = {0.0}; 

  if (cell_avg) { 
  // If |B|^2 < 0 at control points, only use cell average to get 1/|B|^2. 
  magB2_inv[0] = 2.0/magB2[0]; 
  } else { 
  ser_1x_p1_inv(magB2, magB2_inv); 
  } 
  // Calculate expansions of B_i B_j/|B|^2, which can be calculated free of aliasing errors. 
  binop_mul_1d_ser_p1(magB2_inv, B_x_sq, bxbx); 
  binop_mul_1d_ser_p1(magB2_inv, B_y_sq, byby); 
  binop_mul_1d_ser_p1(magB2_inv, B_z_sq, bzbz); 
  binop_mul_1d_ser_p1(magB2_inv, B_x_B_y, bxby); 
  binop_mul_1d_ser_p1(magB2_inv, B_x_B_z, bxbz); 
  binop_mul_1d_ser_p1(magB2_inv, B_y_B_z, bybz); 
 
  int cell_avg_bb = 0;
  if (0.7071067811865475*bxbx[0]-0.7071067811865475*bxbx[1] < 0.0) cell_avg_bb = 1; 
  if (0.7071067811865475*byby[0]-0.7071067811865475*byby[1] < 0.0) cell_avg_bb = 1; 
  if (0.7071067811865475*bzbz[0]-0.7071067811865475*bzbz[1] < 0.0) cell_avg_bb = 1; 
  if (0.7071067811865475*bxbx[1]+0.7071067811865475*bxbx[0] < 0.0) cell_avg_bb = 1; 
  if (0.7071067811865475*byby[1]+0.7071067811865475*byby[0] < 0.0) cell_avg_bb = 1; 
  if (0.7071067811865475*bzbz[1]+0.7071067811865475*bzbz[0] < 0.0) cell_avg_bb = 1; 
  if (cell_avg_bb || cell_avg) { 
    bxbx[1] = 0.0; 
    bxby[1] = 0.0; 
    bxbz[1] = 0.0; 
    byby[1] = 0.0; 
    bybz[1] = 0.0; 
    bzbz[1] = 0.0; 
  } 
  // Calculate b_i = B_i/|B| by taking square root of B_i^2/|B|^2 at quadrature points. 
  // Uses the sign of B_i at quadrature points to get the correct sign of b_i. 
  // Also checks if B_i^2/|B|^2 < 0.0 at quadrature points and zeros out the value there. 
  ser_1x_p1_sqrt_with_sign(B_x, bxbx, bx); 
  ser_1x_p1_sqrt_with_sign(B_y, byby, by); 
  ser_1x_p1_sqrt_with_sign(B_z, bzbz, bz); 
} 
 
