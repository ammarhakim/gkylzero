#include <gkyl_maxwell_kernels.h> 
#include <gkyl_basis_ser_2x_p1_exp_sq.h> 
#include <gkyl_basis_ser_2x_p1_inv.h> 
#include <gkyl_basis_ser_2x_p1_sqrt_with_sign.h> 
GKYL_CU_DH void em_bvar_2x_ser_p1(const double *em, double* bvar) 
{ 
  // em:   Input electromagnetic fields. 
  // bvar: b_i = B_i/|B| (first 3 components), b_i b_j = B_i B_j/|B|^2 (last 6 components). 
 
  const double *B_x = &em[12]; 
  const double *B_y = &em[16]; 
  const double *B_z = &em[20]; 
 
  double *bx = &bvar[0]; 
  double *by = &bvar[4]; 
  double *bz = &bvar[8]; 
  double *bxbx = &bvar[12]; 
  double *bxby = &bvar[16]; 
  double *bxbz = &bvar[20]; 
  double *byby = &bvar[24]; 
  double *bybz = &bvar[28]; 
  double *bzbz = &bvar[32]; 
 
  // Calculate |B|^2 and get expansion of 1/|B|^2. 
  double B_x_sq[4] = {0.0}; 
  ser_2x_p1_exp_sq(B_x, B_x_sq); 
 
  double B_y_sq[4] = {0.0}; 
  ser_2x_p1_exp_sq(B_y, B_y_sq); 
 
  double B_z_sq[4] = {0.0}; 
  ser_2x_p1_exp_sq(B_z, B_z_sq); 
 
  double magB_sq[4] = {0.0}; 

  double B_x_B_y[4] = {0.0}; 
  B_x_B_y[0] = 0.5*B_x[3]*B_y[3]+0.5*B_x[2]*B_y[2]+0.5*B_x[1]*B_y[1]+0.5*B_x[0]*B_y[0]; 
  B_x_B_y[1] = 0.5*B_x[2]*B_y[3]+0.5*B_y[2]*B_x[3]+0.5*B_x[0]*B_y[1]+0.5*B_y[0]*B_x[1]; 
  B_x_B_y[2] = 0.5*B_x[1]*B_y[3]+0.5*B_y[1]*B_x[3]+0.5*B_x[0]*B_y[2]+0.5*B_y[0]*B_x[2]; 
  B_x_B_y[3] = 0.5*B_x[0]*B_y[3]+0.5*B_y[0]*B_x[3]+0.5*B_x[1]*B_y[2]+0.5*B_y[1]*B_x[2]; 

  double B_x_B_z[4] = {0.0}; 
  B_x_B_z[0] = 0.5*B_x[3]*B_z[3]+0.5*B_x[2]*B_z[2]+0.5*B_x[1]*B_z[1]+0.5*B_x[0]*B_z[0]; 
  B_x_B_z[1] = 0.5*B_x[2]*B_z[3]+0.5*B_z[2]*B_x[3]+0.5*B_x[0]*B_z[1]+0.5*B_z[0]*B_x[1]; 
  B_x_B_z[2] = 0.5*B_x[1]*B_z[3]+0.5*B_z[1]*B_x[3]+0.5*B_x[0]*B_z[2]+0.5*B_z[0]*B_x[2]; 
  B_x_B_z[3] = 0.5*B_x[0]*B_z[3]+0.5*B_z[0]*B_x[3]+0.5*B_x[1]*B_z[2]+0.5*B_z[1]*B_x[2]; 

  double B_y_B_z[4] = {0.0}; 
  B_y_B_z[0] = 0.5*B_y[3]*B_z[3]+0.5*B_y[2]*B_z[2]+0.5*B_y[1]*B_z[1]+0.5*B_y[0]*B_z[0]; 
  B_y_B_z[1] = 0.5*B_y[2]*B_z[3]+0.5*B_z[2]*B_y[3]+0.5*B_y[0]*B_z[1]+0.5*B_z[0]*B_y[1]; 
  B_y_B_z[2] = 0.5*B_y[1]*B_z[3]+0.5*B_z[1]*B_y[3]+0.5*B_y[0]*B_z[2]+0.5*B_z[0]*B_y[2]; 
  B_y_B_z[3] = 0.5*B_y[0]*B_z[3]+0.5*B_z[0]*B_y[3]+0.5*B_y[1]*B_z[2]+0.5*B_z[1]*B_y[2]; 

  magB_sq[0] = B_z_sq[0]+B_y_sq[0]+B_x_sq[0]; 
  magB_sq[1] = B_z_sq[1]+B_y_sq[1]+B_x_sq[1]; 
  magB_sq[2] = B_z_sq[2]+B_y_sq[2]+B_x_sq[2]; 
  magB_sq[3] = B_z_sq[3]+B_y_sq[3]+B_x_sq[3]; 

  bool notCellAvg = true;
  // Check if |B|^2 < 0 at control points. 
  if (notCellAvg && (1.5*magB_sq[3]-0.8660254037844386*magB_sq[2]-0.8660254037844386*magB_sq[1]+0.5*magB_sq[0] < 0)) notCellAvg = false; 
  if (notCellAvg && ((-1.5*magB_sq[3])-0.8660254037844386*magB_sq[2]+0.8660254037844386*magB_sq[1]+0.5*magB_sq[0] < 0)) notCellAvg = false; 
  if (notCellAvg && ((-1.5*magB_sq[3])+0.8660254037844386*magB_sq[2]-0.8660254037844386*magB_sq[1]+0.5*magB_sq[0] < 0)) notCellAvg = false; 
  if (notCellAvg && (1.5*magB_sq[3]+0.8660254037844386*magB_sq[2]+0.8660254037844386*magB_sq[1]+0.5*magB_sq[0] < 0)) notCellAvg = false; 
  double magB_sq_inv[4] = {0.0}; 

  // Calculate expansions of B_i B_j/|B|^2, which can be calculated free of aliasing errors. 
  if (notCellAvg) { 
  ser_2x_p1_inv(magB_sq, magB_sq_inv); 
  bxbx[0] = 0.5*B_x_sq[3]*magB_sq_inv[3]+0.5*B_x_sq[2]*magB_sq_inv[2]+0.5*B_x_sq[1]*magB_sq_inv[1]+0.5*B_x_sq[0]*magB_sq_inv[0]; 
  bxbx[1] = 0.5*B_x_sq[2]*magB_sq_inv[3]+0.5*magB_sq_inv[2]*B_x_sq[3]+0.5*B_x_sq[0]*magB_sq_inv[1]+0.5*magB_sq_inv[0]*B_x_sq[1]; 
  bxbx[2] = 0.5*B_x_sq[1]*magB_sq_inv[3]+0.5*magB_sq_inv[1]*B_x_sq[3]+0.5*B_x_sq[0]*magB_sq_inv[2]+0.5*magB_sq_inv[0]*B_x_sq[2]; 
  bxbx[3] = 0.5*B_x_sq[0]*magB_sq_inv[3]+0.5*magB_sq_inv[0]*B_x_sq[3]+0.5*B_x_sq[1]*magB_sq_inv[2]+0.5*magB_sq_inv[1]*B_x_sq[2]; 

  bxby[0] = 0.5*B_x_B_y[3]*magB_sq_inv[3]+0.5*B_x_B_y[2]*magB_sq_inv[2]+0.5*B_x_B_y[1]*magB_sq_inv[1]+0.5*B_x_B_y[0]*magB_sq_inv[0]; 
  bxby[1] = 0.5*B_x_B_y[2]*magB_sq_inv[3]+0.5*magB_sq_inv[2]*B_x_B_y[3]+0.5*B_x_B_y[0]*magB_sq_inv[1]+0.5*magB_sq_inv[0]*B_x_B_y[1]; 
  bxby[2] = 0.5*B_x_B_y[1]*magB_sq_inv[3]+0.5*magB_sq_inv[1]*B_x_B_y[3]+0.5*B_x_B_y[0]*magB_sq_inv[2]+0.5*magB_sq_inv[0]*B_x_B_y[2]; 
  bxby[3] = 0.5*B_x_B_y[0]*magB_sq_inv[3]+0.5*magB_sq_inv[0]*B_x_B_y[3]+0.5*B_x_B_y[1]*magB_sq_inv[2]+0.5*magB_sq_inv[1]*B_x_B_y[2]; 

  bxbz[0] = 0.5*B_x_B_z[3]*magB_sq_inv[3]+0.5*B_x_B_z[2]*magB_sq_inv[2]+0.5*B_x_B_z[1]*magB_sq_inv[1]+0.5*B_x_B_z[0]*magB_sq_inv[0]; 
  bxbz[1] = 0.5*B_x_B_z[2]*magB_sq_inv[3]+0.5*magB_sq_inv[2]*B_x_B_z[3]+0.5*B_x_B_z[0]*magB_sq_inv[1]+0.5*magB_sq_inv[0]*B_x_B_z[1]; 
  bxbz[2] = 0.5*B_x_B_z[1]*magB_sq_inv[3]+0.5*magB_sq_inv[1]*B_x_B_z[3]+0.5*B_x_B_z[0]*magB_sq_inv[2]+0.5*magB_sq_inv[0]*B_x_B_z[2]; 
  bxbz[3] = 0.5*B_x_B_z[0]*magB_sq_inv[3]+0.5*magB_sq_inv[0]*B_x_B_z[3]+0.5*B_x_B_z[1]*magB_sq_inv[2]+0.5*magB_sq_inv[1]*B_x_B_z[2]; 

  byby[0] = 0.5*B_y_sq[3]*magB_sq_inv[3]+0.5*B_y_sq[2]*magB_sq_inv[2]+0.5*B_y_sq[1]*magB_sq_inv[1]+0.5*B_y_sq[0]*magB_sq_inv[0]; 
  byby[1] = 0.5*B_y_sq[2]*magB_sq_inv[3]+0.5*magB_sq_inv[2]*B_y_sq[3]+0.5*B_y_sq[0]*magB_sq_inv[1]+0.5*magB_sq_inv[0]*B_y_sq[1]; 
  byby[2] = 0.5*B_y_sq[1]*magB_sq_inv[3]+0.5*magB_sq_inv[1]*B_y_sq[3]+0.5*B_y_sq[0]*magB_sq_inv[2]+0.5*magB_sq_inv[0]*B_y_sq[2]; 
  byby[3] = 0.5*B_y_sq[0]*magB_sq_inv[3]+0.5*magB_sq_inv[0]*B_y_sq[3]+0.5*B_y_sq[1]*magB_sq_inv[2]+0.5*magB_sq_inv[1]*B_y_sq[2]; 

  bybz[0] = 0.5*B_y_B_z[3]*magB_sq_inv[3]+0.5*B_y_B_z[2]*magB_sq_inv[2]+0.5*B_y_B_z[1]*magB_sq_inv[1]+0.5*B_y_B_z[0]*magB_sq_inv[0]; 
  bybz[1] = 0.5*B_y_B_z[2]*magB_sq_inv[3]+0.5*magB_sq_inv[2]*B_y_B_z[3]+0.5*B_y_B_z[0]*magB_sq_inv[1]+0.5*magB_sq_inv[0]*B_y_B_z[1]; 
  bybz[2] = 0.5*B_y_B_z[1]*magB_sq_inv[3]+0.5*magB_sq_inv[1]*B_y_B_z[3]+0.5*B_y_B_z[0]*magB_sq_inv[2]+0.5*magB_sq_inv[0]*B_y_B_z[2]; 
  bybz[3] = 0.5*B_y_B_z[0]*magB_sq_inv[3]+0.5*magB_sq_inv[0]*B_y_B_z[3]+0.5*B_y_B_z[1]*magB_sq_inv[2]+0.5*magB_sq_inv[1]*B_y_B_z[2]; 

  bzbz[0] = 0.5*B_z_sq[3]*magB_sq_inv[3]+0.5*B_z_sq[2]*magB_sq_inv[2]+0.5*B_z_sq[1]*magB_sq_inv[1]+0.5*B_z_sq[0]*magB_sq_inv[0]; 
  bzbz[1] = 0.5*B_z_sq[2]*magB_sq_inv[3]+0.5*magB_sq_inv[2]*B_z_sq[3]+0.5*B_z_sq[0]*magB_sq_inv[1]+0.5*magB_sq_inv[0]*B_z_sq[1]; 
  bzbz[2] = 0.5*B_z_sq[1]*magB_sq_inv[3]+0.5*magB_sq_inv[1]*B_z_sq[3]+0.5*B_z_sq[0]*magB_sq_inv[2]+0.5*magB_sq_inv[0]*B_z_sq[2]; 
  bzbz[3] = 0.5*B_z_sq[0]*magB_sq_inv[3]+0.5*magB_sq_inv[0]*B_z_sq[3]+0.5*B_z_sq[1]*magB_sq_inv[2]+0.5*magB_sq_inv[1]*B_z_sq[2]; 

  } else { 
  // If |B|^2 at control points, only use cell average to get 1/|B|^2. 
  magB_sq_inv[0] = 4.0/magB_sq[0]; 
  bxbx[0] = 0.5*B_x_sq[3]*magB_sq_inv[3]+0.5*B_x_sq[2]*magB_sq_inv[2]+0.5*B_x_sq[1]*magB_sq_inv[1]+0.5*B_x_sq[0]*magB_sq_inv[0]; 
  bxby[0] = 0.5*B_x_B_y[3]*magB_sq_inv[3]+0.5*B_x_B_y[2]*magB_sq_inv[2]+0.5*B_x_B_y[1]*magB_sq_inv[1]+0.5*B_x_B_y[0]*magB_sq_inv[0]; 
  bxbz[0] = 0.5*B_x_B_z[3]*magB_sq_inv[3]+0.5*B_x_B_z[2]*magB_sq_inv[2]+0.5*B_x_B_z[1]*magB_sq_inv[1]+0.5*B_x_B_z[0]*magB_sq_inv[0]; 
  byby[0] = 0.5*B_y_sq[3]*magB_sq_inv[3]+0.5*B_y_sq[2]*magB_sq_inv[2]+0.5*B_y_sq[1]*magB_sq_inv[1]+0.5*B_y_sq[0]*magB_sq_inv[0]; 
  bybz[0] = 0.5*B_y_B_z[3]*magB_sq_inv[3]+0.5*B_y_B_z[2]*magB_sq_inv[2]+0.5*B_y_B_z[1]*magB_sq_inv[1]+0.5*B_y_B_z[0]*magB_sq_inv[0]; 
  bzbz[0] = 0.5*B_z_sq[3]*magB_sq_inv[3]+0.5*B_z_sq[2]*magB_sq_inv[2]+0.5*B_z_sq[1]*magB_sq_inv[1]+0.5*B_z_sq[0]*magB_sq_inv[0]; 

  bxbx[1] = 0.0; 
  bxby[1] = 0.0; 
  bxbz[1] = 0.0; 
  byby[1] = 0.0; 
  bybz[1] = 0.0; 
  bzbz[1] = 0.0; 

  bxbx[2] = 0.0; 
  bxby[2] = 0.0; 
  bxbz[2] = 0.0; 
  byby[2] = 0.0; 
  bybz[2] = 0.0; 
  bzbz[2] = 0.0; 

  bxbx[3] = 0.0; 
  bxby[3] = 0.0; 
  bxbz[3] = 0.0; 
  byby[3] = 0.0; 
  bybz[3] = 0.0; 
  bzbz[3] = 0.0; 

  } 
  // Calculate b_i = B_i/|B| by taking square root of B_i^2/|B|^2 at quadrature points. 
  // Uses the sign of B_i at quadrature points to get the correct sign of b_i. 
  // Also checks if B_i^2/|B|^2 < 0.0 at quadrature points and zeros out the value there. 
  ser_2x_p1_sqrt_with_sign(B_x, bxbx, bx); 
  ser_2x_p1_sqrt_with_sign(B_y, byby, by); 
  ser_2x_p1_sqrt_with_sign(B_z, bzbz, bz); 
} 
 
