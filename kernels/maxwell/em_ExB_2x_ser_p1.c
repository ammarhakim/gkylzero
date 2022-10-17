#include <gkyl_maxwell_kernels.h> 
#include <gkyl_basis_ser_2x_p1_exp_sq.h> 
#include <gkyl_basis_ser_2x_p1_inv.h> 
GKYL_CU_DH void em_ExB_2x_ser_p1(const double *em, double* GKYL_RESTRICT ExB) 
{ 
  // em:  Input electromagnetic fields. 
  // ExB: E x B velocity = E x B/|B|^2. 
 
  const double *E_x = &em[0]; 
  const double *E_y = &em[4]; 
  const double *E_z = &em[8]; 
  const double *B_x = &em[12]; 
  const double *B_y = &em[16]; 
  const double *B_z = &em[20]; 
 
  double *ExB_x = &ExB[0]; 
  double *ExB_y = &ExB[4]; 
  double *ExB_z = &ExB[8]; 
 
  // Calculate |B|^2 and get expansion of 1/|B|^2. 
  double B_x_sq[4] = {0.0}; 
  ser_2x_p1_exp_sq(B_x, B_x_sq); 
 
  double B_y_sq[4] = {0.0}; 
  ser_2x_p1_exp_sq(B_y, B_y_sq); 
 
  double B_z_sq[4] = {0.0}; 
  ser_2x_p1_exp_sq(B_z, B_z_sq); 
 
  double magB_sq[4] = {0.0}; 

  double num_ExB_x[4] = {0.0}; 
  num_ExB_x[0] = (-0.5*B_y[3]*E_z[3])+0.5*B_z[3]*E_y[3]-0.5*B_y[2]*E_z[2]+0.5*B_z[2]*E_y[2]-0.5*B_y[1]*E_z[1]+0.5*B_z[1]*E_y[1]-0.5*B_y[0]*E_z[0]+0.5*B_z[0]*E_y[0]; 
  num_ExB_x[1] = (-0.5*B_y[2]*E_z[3])+0.5*B_z[2]*E_y[3]+0.5*E_y[2]*B_z[3]-0.5*E_z[2]*B_y[3]-0.5*B_y[0]*E_z[1]+0.5*B_z[0]*E_y[1]+0.5*E_y[0]*B_z[1]-0.5*E_z[0]*B_y[1]; 
  num_ExB_x[2] = (-0.5*B_y[1]*E_z[3])+0.5*B_z[1]*E_y[3]+0.5*E_y[1]*B_z[3]-0.5*E_z[1]*B_y[3]-0.5*B_y[0]*E_z[2]+0.5*B_z[0]*E_y[2]+0.5*E_y[0]*B_z[2]-0.5*E_z[0]*B_y[2]; 
  num_ExB_x[3] = (-0.5*B_y[0]*E_z[3])+0.5*B_z[0]*E_y[3]+0.5*E_y[0]*B_z[3]-0.5*E_z[0]*B_y[3]-0.5*B_y[1]*E_z[2]+0.5*B_z[1]*E_y[2]+0.5*E_y[1]*B_z[2]-0.5*E_z[1]*B_y[2]; 

  double num_ExB_y[4] = {0.0}; 
  num_ExB_y[0] = 0.5*B_x[3]*E_z[3]-0.5*B_z[3]*E_x[3]+0.5*B_x[2]*E_z[2]-0.5*B_z[2]*E_x[2]+0.5*B_x[1]*E_z[1]-0.5*B_z[1]*E_x[1]+0.5*B_x[0]*E_z[0]-0.5*B_z[0]*E_x[0]; 
  num_ExB_y[1] = 0.5*B_x[2]*E_z[3]-0.5*B_z[2]*E_x[3]-0.5*E_x[2]*B_z[3]+0.5*E_z[2]*B_x[3]+0.5*B_x[0]*E_z[1]-0.5*B_z[0]*E_x[1]-0.5*E_x[0]*B_z[1]+0.5*E_z[0]*B_x[1]; 
  num_ExB_y[2] = 0.5*B_x[1]*E_z[3]-0.5*B_z[1]*E_x[3]-0.5*E_x[1]*B_z[3]+0.5*E_z[1]*B_x[3]+0.5*B_x[0]*E_z[2]-0.5*B_z[0]*E_x[2]-0.5*E_x[0]*B_z[2]+0.5*E_z[0]*B_x[2]; 
  num_ExB_y[3] = 0.5*B_x[0]*E_z[3]-0.5*B_z[0]*E_x[3]-0.5*E_x[0]*B_z[3]+0.5*E_z[0]*B_x[3]+0.5*B_x[1]*E_z[2]-0.5*B_z[1]*E_x[2]-0.5*E_x[1]*B_z[2]+0.5*E_z[1]*B_x[2]; 

  double num_ExB_z[4] = {0.0}; 
  num_ExB_z[0] = (-0.5*B_x[3]*E_y[3])+0.5*B_y[3]*E_x[3]-0.5*B_x[2]*E_y[2]+0.5*B_y[2]*E_x[2]-0.5*B_x[1]*E_y[1]+0.5*B_y[1]*E_x[1]-0.5*B_x[0]*E_y[0]+0.5*B_y[0]*E_x[0]; 
  num_ExB_z[1] = (-0.5*B_x[2]*E_y[3])+0.5*B_y[2]*E_x[3]+0.5*E_x[2]*B_y[3]-0.5*E_y[2]*B_x[3]-0.5*B_x[0]*E_y[1]+0.5*B_y[0]*E_x[1]+0.5*E_x[0]*B_y[1]-0.5*E_y[0]*B_x[1]; 
  num_ExB_z[2] = (-0.5*B_x[1]*E_y[3])+0.5*B_y[1]*E_x[3]+0.5*E_x[1]*B_y[3]-0.5*E_y[1]*B_x[3]-0.5*B_x[0]*E_y[2]+0.5*B_y[0]*E_x[2]+0.5*E_x[0]*B_y[2]-0.5*E_y[0]*B_x[2]; 
  num_ExB_z[3] = (-0.5*B_x[0]*E_y[3])+0.5*B_y[0]*E_x[3]+0.5*E_x[0]*B_y[3]-0.5*E_y[0]*B_x[3]-0.5*B_x[1]*E_y[2]+0.5*B_y[1]*E_x[2]+0.5*E_x[1]*B_y[2]-0.5*E_y[1]*B_x[2]; 

  magB_sq[0] = B_z_sq[0]+B_y_sq[0]+B_x_sq[0]; 
  magB_sq[1] = B_z_sq[1]+B_y_sq[1]+B_x_sq[1]; 
  magB_sq[2] = B_z_sq[2]+B_y_sq[2]+B_x_sq[2]; 
  magB_sq[3] = B_z_sq[3]+B_y_sq[3]+B_x_sq[3]; 

  bool notCellAvg = true;
  if (notCellAvg && (1.5*magB_sq[3]-0.8660254037844386*magB_sq[2]-0.8660254037844386*magB_sq[1]+0.5*magB_sq[0] < 0)) notCellAvg = false; 
  if (notCellAvg && ((-1.5*magB_sq[3])-0.8660254037844386*magB_sq[2]+0.8660254037844386*magB_sq[1]+0.5*magB_sq[0] < 0)) notCellAvg = false; 
  if (notCellAvg && ((-1.5*magB_sq[3])+0.8660254037844386*magB_sq[2]-0.8660254037844386*magB_sq[1]+0.5*magB_sq[0] < 0)) notCellAvg = false; 
  if (notCellAvg && (1.5*magB_sq[3]+0.8660254037844386*magB_sq[2]+0.8660254037844386*magB_sq[1]+0.5*magB_sq[0] < 0)) notCellAvg = false; 
  double magB_sq_inv[4] = {0.0}; 

  if (notCellAvg) { 
  ser_2x_p1_inv(magB_sq, magB_sq_inv); 
  ExB_x[0] = 0.5*magB_sq_inv[3]*num_ExB_x[3]+0.5*magB_sq_inv[2]*num_ExB_x[2]+0.5*magB_sq_inv[1]*num_ExB_x[1]+0.5*magB_sq_inv[0]*num_ExB_x[0]; 
  ExB_x[1] = 0.5*magB_sq_inv[2]*num_ExB_x[3]+0.5*num_ExB_x[2]*magB_sq_inv[3]+0.5*magB_sq_inv[0]*num_ExB_x[1]+0.5*num_ExB_x[0]*magB_sq_inv[1]; 
  ExB_x[2] = 0.5*magB_sq_inv[1]*num_ExB_x[3]+0.5*num_ExB_x[1]*magB_sq_inv[3]+0.5*magB_sq_inv[0]*num_ExB_x[2]+0.5*num_ExB_x[0]*magB_sq_inv[2]; 
  ExB_x[3] = 0.5*magB_sq_inv[0]*num_ExB_x[3]+0.5*num_ExB_x[0]*magB_sq_inv[3]+0.5*magB_sq_inv[1]*num_ExB_x[2]+0.5*num_ExB_x[1]*magB_sq_inv[2]; 

  ExB_y[0] = 0.5*magB_sq_inv[3]*num_ExB_y[3]+0.5*magB_sq_inv[2]*num_ExB_y[2]+0.5*magB_sq_inv[1]*num_ExB_y[1]+0.5*magB_sq_inv[0]*num_ExB_y[0]; 
  ExB_y[1] = 0.5*magB_sq_inv[2]*num_ExB_y[3]+0.5*num_ExB_y[2]*magB_sq_inv[3]+0.5*magB_sq_inv[0]*num_ExB_y[1]+0.5*num_ExB_y[0]*magB_sq_inv[1]; 
  ExB_y[2] = 0.5*magB_sq_inv[1]*num_ExB_y[3]+0.5*num_ExB_y[1]*magB_sq_inv[3]+0.5*magB_sq_inv[0]*num_ExB_y[2]+0.5*num_ExB_y[0]*magB_sq_inv[2]; 
  ExB_y[3] = 0.5*magB_sq_inv[0]*num_ExB_y[3]+0.5*num_ExB_y[0]*magB_sq_inv[3]+0.5*magB_sq_inv[1]*num_ExB_y[2]+0.5*num_ExB_y[1]*magB_sq_inv[2]; 

  ExB_z[0] = 0.5*magB_sq_inv[3]*num_ExB_z[3]+0.5*magB_sq_inv[2]*num_ExB_z[2]+0.5*magB_sq_inv[1]*num_ExB_z[1]+0.5*magB_sq_inv[0]*num_ExB_z[0]; 
  ExB_z[1] = 0.5*magB_sq_inv[2]*num_ExB_z[3]+0.5*num_ExB_z[2]*magB_sq_inv[3]+0.5*magB_sq_inv[0]*num_ExB_z[1]+0.5*num_ExB_z[0]*magB_sq_inv[1]; 
  ExB_z[2] = 0.5*magB_sq_inv[1]*num_ExB_z[3]+0.5*num_ExB_z[1]*magB_sq_inv[3]+0.5*magB_sq_inv[0]*num_ExB_z[2]+0.5*num_ExB_z[0]*magB_sq_inv[2]; 
  ExB_z[3] = 0.5*magB_sq_inv[0]*num_ExB_z[3]+0.5*num_ExB_z[0]*magB_sq_inv[3]+0.5*magB_sq_inv[1]*num_ExB_z[2]+0.5*num_ExB_z[1]*magB_sq_inv[2]; 

  } else { 
  magB_sq_inv[0] = 4.0/magB_sq[0]; 
  ExB_x[0] = 0.5*magB_sq_inv[3]*num_ExB_x[3]+0.5*magB_sq_inv[2]*num_ExB_x[2]+0.5*magB_sq_inv[1]*num_ExB_x[1]+0.5*magB_sq_inv[0]*num_ExB_x[0]; 
  ExB_y[0] = 0.5*magB_sq_inv[3]*num_ExB_y[3]+0.5*magB_sq_inv[2]*num_ExB_y[2]+0.5*magB_sq_inv[1]*num_ExB_y[1]+0.5*magB_sq_inv[0]*num_ExB_y[0]; 
  ExB_z[0] = 0.5*magB_sq_inv[3]*num_ExB_z[3]+0.5*magB_sq_inv[2]*num_ExB_z[2]+0.5*magB_sq_inv[1]*num_ExB_z[1]+0.5*magB_sq_inv[0]*num_ExB_z[0]; 

  } 
} 
 
