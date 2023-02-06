#include <gkyl_maxwell_kernels.h> 
#include <gkyl_basis_ser_1x_p2_exp_sq.h> 
#include <gkyl_basis_ser_1x_p2_inv.h> 
#include <gkyl_basis_ser_1x_p2_sqrt_with_sign.h> 
GKYL_CU_DH void em_bvar_1x_ser_p2(const double *em, double* bvar) 
{ 
  // em:   Input electromagnetic fields. 
  // bvar: b_i = B_i/|B| (first 3 components), b_i b_j = B_i B_j/|B|^2 (last 6 components). 
 
  const double *B_x = &em[9]; 
  const double *B_y = &em[12]; 
  const double *B_z = &em[15]; 
 
  double *bx = &bvar[0]; 
  double *by = &bvar[3]; 
  double *bz = &bvar[6]; 
  double *bxbx = &bvar[9]; 
  double *bxby = &bvar[12]; 
  double *bxbz = &bvar[15]; 
  double *byby = &bvar[18]; 
  double *bybz = &bvar[21]; 
  double *bzbz = &bvar[24]; 
 
  // Calculate |B|^2 and get expansion of 1/|B|^2. 
  double B_x_sq[3] = {0.0}; 
  ser_1x_p2_exp_sq(B_x, B_x_sq); 
 
  double B_y_sq[3] = {0.0}; 
  ser_1x_p2_exp_sq(B_y, B_y_sq); 
 
  double B_z_sq[3] = {0.0}; 
  ser_1x_p2_exp_sq(B_z, B_z_sq); 
 
  double magB_sq[3] = {0.0}; 

  double B_x_B_y[3] = {0.0}; 
  B_x_B_y[0] = 0.7071067811865475*B_x[2]*B_y[2]+0.7071067811865475*B_x[1]*B_y[1]+0.7071067811865475*B_x[0]*B_y[0]; 
  B_x_B_y[1] = 0.6324555320336759*B_x[1]*B_y[2]+0.6324555320336759*B_y[1]*B_x[2]+0.7071067811865475*B_x[0]*B_y[1]+0.7071067811865475*B_y[0]*B_x[1]; 
  B_x_B_y[2] = 0.4517539514526256*B_x[2]*B_y[2]+0.7071067811865475*B_x[0]*B_y[2]+0.7071067811865475*B_y[0]*B_x[2]+0.6324555320336759*B_x[1]*B_y[1]; 

  double B_x_B_z[3] = {0.0}; 
  B_x_B_z[0] = 0.7071067811865475*B_x[2]*B_z[2]+0.7071067811865475*B_x[1]*B_z[1]+0.7071067811865475*B_x[0]*B_z[0]; 
  B_x_B_z[1] = 0.6324555320336759*B_x[1]*B_z[2]+0.6324555320336759*B_z[1]*B_x[2]+0.7071067811865475*B_x[0]*B_z[1]+0.7071067811865475*B_z[0]*B_x[1]; 
  B_x_B_z[2] = 0.4517539514526256*B_x[2]*B_z[2]+0.7071067811865475*B_x[0]*B_z[2]+0.7071067811865475*B_z[0]*B_x[2]+0.6324555320336759*B_x[1]*B_z[1]; 

  double B_y_B_z[3] = {0.0}; 
  B_y_B_z[0] = 0.7071067811865475*B_y[2]*B_z[2]+0.7071067811865475*B_y[1]*B_z[1]+0.7071067811865475*B_y[0]*B_z[0]; 
  B_y_B_z[1] = 0.6324555320336759*B_y[1]*B_z[2]+0.6324555320336759*B_z[1]*B_y[2]+0.7071067811865475*B_y[0]*B_z[1]+0.7071067811865475*B_z[0]*B_y[1]; 
  B_y_B_z[2] = 0.4517539514526256*B_y[2]*B_z[2]+0.7071067811865475*B_y[0]*B_z[2]+0.7071067811865475*B_z[0]*B_y[2]+0.6324555320336759*B_y[1]*B_z[1]; 

  magB_sq[0] = B_z_sq[0]+B_y_sq[0]+B_x_sq[0]; 
  magB_sq[1] = B_z_sq[1]+B_y_sq[1]+B_x_sq[1]; 
  magB_sq[2] = B_z_sq[2]+B_y_sq[2]+B_x_sq[2]; 

  bool notCellAvg = true;
  // Check if either Bx^2, By^2, or Bz^2 < 0 at control points (Gauss-Legendre quadrature points). 
  // Or if |B|^2 < 0 at cell corners. 
  if (notCellAvg && (0.6324555320336759*B_x_sq[2]-0.9486832980505137*B_x_sq[1]+0.7071067811865475*B_x_sq[0] < 0)) notCellAvg = false; 
  if (notCellAvg && (0.6324555320336759*B_y_sq[2]-0.9486832980505137*B_y_sq[1]+0.7071067811865475*B_y_sq[0] < 0)) notCellAvg = false; 
  if (notCellAvg && (0.6324555320336759*B_z_sq[2]-0.9486832980505137*B_z_sq[1]+0.7071067811865475*B_z_sq[0] < 0)) notCellAvg = false; 
  if (notCellAvg && (0.7071067811865475*B_x_sq[0]-0.7905694150420947*B_x_sq[2] < 0)) notCellAvg = false; 
  if (notCellAvg && (0.7071067811865475*B_y_sq[0]-0.7905694150420947*B_y_sq[2] < 0)) notCellAvg = false; 
  if (notCellAvg && (0.7071067811865475*B_z_sq[0]-0.7905694150420947*B_z_sq[2] < 0)) notCellAvg = false; 
  if (notCellAvg && (0.6324555320336759*B_x_sq[2]+0.9486832980505137*B_x_sq[1]+0.7071067811865475*B_x_sq[0] < 0)) notCellAvg = false; 
  if (notCellAvg && (0.6324555320336759*B_y_sq[2]+0.9486832980505137*B_y_sq[1]+0.7071067811865475*B_y_sq[0] < 0)) notCellAvg = false; 
  if (notCellAvg && (0.6324555320336759*B_z_sq[2]+0.9486832980505137*B_z_sq[1]+0.7071067811865475*B_z_sq[0] < 0)) notCellAvg = false; 
  if (notCellAvg && (1.58113883008419*magB_sq[2]-1.224744871391589*magB_sq[1]+0.7071067811865475*magB_sq[0] < 0)) notCellAvg = false; 
  if (notCellAvg && (1.58113883008419*magB_sq[2]+1.224744871391589*magB_sq[1]+0.7071067811865475*magB_sq[0] < 0)) notCellAvg = false; 
  double magB_sq_inv[3] = {0.0}; 

  // Calculate expansions of B_i B_j/|B|^2, which can be calculated free of aliasing errors. 
  if (notCellAvg) { 
  ser_1x_p2_inv(magB_sq, magB_sq_inv); 
  bxbx[0] = 0.7071067811865475*B_x_sq[2]*magB_sq_inv[2]+0.7071067811865475*B_x_sq[1]*magB_sq_inv[1]+0.7071067811865475*B_x_sq[0]*magB_sq_inv[0]; 
  bxbx[1] = 0.6324555320336759*B_x_sq[1]*magB_sq_inv[2]+0.6324555320336759*magB_sq_inv[1]*B_x_sq[2]+0.7071067811865475*B_x_sq[0]*magB_sq_inv[1]+0.7071067811865475*magB_sq_inv[0]*B_x_sq[1]; 
  bxbx[2] = 0.4517539514526256*B_x_sq[2]*magB_sq_inv[2]+0.7071067811865475*B_x_sq[0]*magB_sq_inv[2]+0.7071067811865475*magB_sq_inv[0]*B_x_sq[2]+0.6324555320336759*B_x_sq[1]*magB_sq_inv[1]; 

  bxby[0] = 0.7071067811865475*B_x_B_y[2]*magB_sq_inv[2]+0.7071067811865475*B_x_B_y[1]*magB_sq_inv[1]+0.7071067811865475*B_x_B_y[0]*magB_sq_inv[0]; 
  bxby[1] = 0.6324555320336759*B_x_B_y[1]*magB_sq_inv[2]+0.6324555320336759*magB_sq_inv[1]*B_x_B_y[2]+0.7071067811865475*B_x_B_y[0]*magB_sq_inv[1]+0.7071067811865475*magB_sq_inv[0]*B_x_B_y[1]; 
  bxby[2] = 0.4517539514526256*B_x_B_y[2]*magB_sq_inv[2]+0.7071067811865475*B_x_B_y[0]*magB_sq_inv[2]+0.7071067811865475*magB_sq_inv[0]*B_x_B_y[2]+0.6324555320336759*B_x_B_y[1]*magB_sq_inv[1]; 

  bxbz[0] = 0.7071067811865475*B_x_B_z[2]*magB_sq_inv[2]+0.7071067811865475*B_x_B_z[1]*magB_sq_inv[1]+0.7071067811865475*B_x_B_z[0]*magB_sq_inv[0]; 
  bxbz[1] = 0.6324555320336759*B_x_B_z[1]*magB_sq_inv[2]+0.6324555320336759*magB_sq_inv[1]*B_x_B_z[2]+0.7071067811865475*B_x_B_z[0]*magB_sq_inv[1]+0.7071067811865475*magB_sq_inv[0]*B_x_B_z[1]; 
  bxbz[2] = 0.4517539514526256*B_x_B_z[2]*magB_sq_inv[2]+0.7071067811865475*B_x_B_z[0]*magB_sq_inv[2]+0.7071067811865475*magB_sq_inv[0]*B_x_B_z[2]+0.6324555320336759*B_x_B_z[1]*magB_sq_inv[1]; 

  byby[0] = 0.7071067811865475*B_y_sq[2]*magB_sq_inv[2]+0.7071067811865475*B_y_sq[1]*magB_sq_inv[1]+0.7071067811865475*B_y_sq[0]*magB_sq_inv[0]; 
  byby[1] = 0.6324555320336759*B_y_sq[1]*magB_sq_inv[2]+0.6324555320336759*magB_sq_inv[1]*B_y_sq[2]+0.7071067811865475*B_y_sq[0]*magB_sq_inv[1]+0.7071067811865475*magB_sq_inv[0]*B_y_sq[1]; 
  byby[2] = 0.4517539514526256*B_y_sq[2]*magB_sq_inv[2]+0.7071067811865475*B_y_sq[0]*magB_sq_inv[2]+0.7071067811865475*magB_sq_inv[0]*B_y_sq[2]+0.6324555320336759*B_y_sq[1]*magB_sq_inv[1]; 

  bybz[0] = 0.7071067811865475*B_y_B_z[2]*magB_sq_inv[2]+0.7071067811865475*B_y_B_z[1]*magB_sq_inv[1]+0.7071067811865475*B_y_B_z[0]*magB_sq_inv[0]; 
  bybz[1] = 0.6324555320336759*B_y_B_z[1]*magB_sq_inv[2]+0.6324555320336759*magB_sq_inv[1]*B_y_B_z[2]+0.7071067811865475*B_y_B_z[0]*magB_sq_inv[1]+0.7071067811865475*magB_sq_inv[0]*B_y_B_z[1]; 
  bybz[2] = 0.4517539514526256*B_y_B_z[2]*magB_sq_inv[2]+0.7071067811865475*B_y_B_z[0]*magB_sq_inv[2]+0.7071067811865475*magB_sq_inv[0]*B_y_B_z[2]+0.6324555320336759*B_y_B_z[1]*magB_sq_inv[1]; 

  bzbz[0] = 0.7071067811865475*B_z_sq[2]*magB_sq_inv[2]+0.7071067811865475*B_z_sq[1]*magB_sq_inv[1]+0.7071067811865475*B_z_sq[0]*magB_sq_inv[0]; 
  bzbz[1] = 0.6324555320336759*B_z_sq[1]*magB_sq_inv[2]+0.6324555320336759*magB_sq_inv[1]*B_z_sq[2]+0.7071067811865475*B_z_sq[0]*magB_sq_inv[1]+0.7071067811865475*magB_sq_inv[0]*B_z_sq[1]; 
  bzbz[2] = 0.4517539514526256*B_z_sq[2]*magB_sq_inv[2]+0.7071067811865475*B_z_sq[0]*magB_sq_inv[2]+0.7071067811865475*magB_sq_inv[0]*B_z_sq[2]+0.6324555320336759*B_z_sq[1]*magB_sq_inv[1]; 

  } else { 
  // If either Bx^2, By^2, or Bz^2 < 0 at control points, only use cell average to get 1/|B|^2. 
  magB_sq_inv[0] = 2.0/magB_sq[0]; 
  bxbx[0] = 0.7071067811865475*B_x_sq[2]*magB_sq_inv[2]+0.7071067811865475*B_x_sq[1]*magB_sq_inv[1]+0.7071067811865475*B_x_sq[0]*magB_sq_inv[0]; 
  bxby[0] = 0.7071067811865475*B_x_B_y[2]*magB_sq_inv[2]+0.7071067811865475*B_x_B_y[1]*magB_sq_inv[1]+0.7071067811865475*B_x_B_y[0]*magB_sq_inv[0]; 
  bxbz[0] = 0.7071067811865475*B_x_B_z[2]*magB_sq_inv[2]+0.7071067811865475*B_x_B_z[1]*magB_sq_inv[1]+0.7071067811865475*B_x_B_z[0]*magB_sq_inv[0]; 
  byby[0] = 0.7071067811865475*B_y_sq[2]*magB_sq_inv[2]+0.7071067811865475*B_y_sq[1]*magB_sq_inv[1]+0.7071067811865475*B_y_sq[0]*magB_sq_inv[0]; 
  bybz[0] = 0.7071067811865475*B_y_B_z[2]*magB_sq_inv[2]+0.7071067811865475*B_y_B_z[1]*magB_sq_inv[1]+0.7071067811865475*B_y_B_z[0]*magB_sq_inv[0]; 
  bzbz[0] = 0.7071067811865475*B_z_sq[2]*magB_sq_inv[2]+0.7071067811865475*B_z_sq[1]*magB_sq_inv[1]+0.7071067811865475*B_z_sq[0]*magB_sq_inv[0]; 

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

  } 
  // Calculate b_i = B_i/|B| by taking square root of B_i^2/|B|^2 at quadrature points. 
  ser_1x_p2_sqrt_with_sign(bxbx, bx); 
  ser_1x_p2_sqrt_with_sign(byby, by); 
  ser_1x_p2_sqrt_with_sign(bzbz, bz); 
} 
 
