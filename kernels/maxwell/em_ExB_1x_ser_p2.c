#include <gkyl_maxwell_kernels.h> 
#include <gkyl_basis_ser_1x_p2_exp_sq.h> 
#include <gkyl_basis_ser_1x_p2_inv.h> 
GKYL_CU_DH void em_ExB_1x_ser_p2(const double *em, double* GKYL_RESTRICT ExB) 
{ 
  // em:  Input electromagnetic fields. 
  // ExB: E x B velocity = E x B/|B|^2. 
 
  const double *E_x = &em[0]; 
  const double *E_y = &em[3]; 
  const double *E_z = &em[6]; 
  const double *B_x = &em[9]; 
  const double *B_y = &em[12]; 
  const double *B_z = &em[15]; 
 
  double *ExB_x = &ExB[0]; 
  double *ExB_y = &ExB[3]; 
  double *ExB_z = &ExB[6]; 
 
  // Calculate |B|^2 and get expansion of 1/|B|^2. 
  double B_x_sq[3] = {0.0}; 
  ser_1x_p2_exp_sq(B_x, B_x_sq); 
 
  double B_y_sq[3] = {0.0}; 
  ser_1x_p2_exp_sq(B_y, B_y_sq); 
 
  double B_z_sq[3] = {0.0}; 
  ser_1x_p2_exp_sq(B_z, B_z_sq); 
 
  double magB_sq[3] = {0.0}; 

  double num_ExB_x[3] = {0.0}; 
  num_ExB_x[0] = (-0.7071067811865475*B_y[2]*E_z[2])+0.7071067811865475*B_z[2]*E_y[2]-0.7071067811865475*B_y[1]*E_z[1]+0.7071067811865475*B_z[1]*E_y[1]-0.7071067811865475*B_y[0]*E_z[0]+0.7071067811865475*B_z[0]*E_y[0]; 
  num_ExB_x[1] = (-0.6324555320336759*B_y[1]*E_z[2])+0.6324555320336759*B_z[1]*E_y[2]+0.6324555320336759*E_y[1]*B_z[2]-0.6324555320336759*E_z[1]*B_y[2]-0.7071067811865475*B_y[0]*E_z[1]+0.7071067811865475*B_z[0]*E_y[1]+0.7071067811865475*E_y[0]*B_z[1]-0.7071067811865475*E_z[0]*B_y[1]; 
  num_ExB_x[2] = (-0.4517539514526256*B_y[2]*E_z[2])-0.7071067811865475*B_y[0]*E_z[2]+0.4517539514526256*B_z[2]*E_y[2]+0.7071067811865475*B_z[0]*E_y[2]+0.7071067811865475*E_y[0]*B_z[2]-0.7071067811865475*E_z[0]*B_y[2]-0.6324555320336759*B_y[1]*E_z[1]+0.6324555320336759*B_z[1]*E_y[1]; 

  double num_ExB_y[3] = {0.0}; 
  num_ExB_y[0] = 0.7071067811865475*B_x[2]*E_z[2]-0.7071067811865475*B_z[2]*E_x[2]+0.7071067811865475*B_x[1]*E_z[1]-0.7071067811865475*B_z[1]*E_x[1]+0.7071067811865475*B_x[0]*E_z[0]-0.7071067811865475*B_z[0]*E_x[0]; 
  num_ExB_y[1] = 0.6324555320336759*B_x[1]*E_z[2]-0.6324555320336759*B_z[1]*E_x[2]-0.6324555320336759*E_x[1]*B_z[2]+0.6324555320336759*E_z[1]*B_x[2]+0.7071067811865475*B_x[0]*E_z[1]-0.7071067811865475*B_z[0]*E_x[1]-0.7071067811865475*E_x[0]*B_z[1]+0.7071067811865475*E_z[0]*B_x[1]; 
  num_ExB_y[2] = 0.4517539514526256*B_x[2]*E_z[2]+0.7071067811865475*B_x[0]*E_z[2]-0.4517539514526256*B_z[2]*E_x[2]-0.7071067811865475*B_z[0]*E_x[2]-0.7071067811865475*E_x[0]*B_z[2]+0.7071067811865475*E_z[0]*B_x[2]+0.6324555320336759*B_x[1]*E_z[1]-0.6324555320336759*B_z[1]*E_x[1]; 

  double num_ExB_z[3] = {0.0}; 
  num_ExB_z[0] = (-0.7071067811865475*B_x[2]*E_y[2])+0.7071067811865475*B_y[2]*E_x[2]-0.7071067811865475*B_x[1]*E_y[1]+0.7071067811865475*B_y[1]*E_x[1]-0.7071067811865475*B_x[0]*E_y[0]+0.7071067811865475*B_y[0]*E_x[0]; 
  num_ExB_z[1] = (-0.6324555320336759*B_x[1]*E_y[2])+0.6324555320336759*B_y[1]*E_x[2]+0.6324555320336759*E_x[1]*B_y[2]-0.6324555320336759*E_y[1]*B_x[2]-0.7071067811865475*B_x[0]*E_y[1]+0.7071067811865475*B_y[0]*E_x[1]+0.7071067811865475*E_x[0]*B_y[1]-0.7071067811865475*E_y[0]*B_x[1]; 
  num_ExB_z[2] = (-0.4517539514526256*B_x[2]*E_y[2])-0.7071067811865475*B_x[0]*E_y[2]+0.4517539514526256*B_y[2]*E_x[2]+0.7071067811865475*B_y[0]*E_x[2]+0.7071067811865475*E_x[0]*B_y[2]-0.7071067811865475*E_y[0]*B_x[2]-0.6324555320336759*B_x[1]*E_y[1]+0.6324555320336759*B_y[1]*E_x[1]; 

  magB_sq[0] = B_z_sq[0]+B_y_sq[0]+B_x_sq[0]; 
  magB_sq[1] = B_z_sq[1]+B_y_sq[1]+B_x_sq[1]; 
  magB_sq[2] = B_z_sq[2]+B_y_sq[2]+B_x_sq[2]; 

  bool notCellAvg = true;
  if (notCellAvg && (1.58113883008419*magB_sq[2]-1.224744871391589*magB_sq[1]+0.7071067811865475*magB_sq[0] < 0)) notCellAvg = false; 
  if (notCellAvg && (1.58113883008419*magB_sq[2]+1.224744871391589*magB_sq[1]+0.7071067811865475*magB_sq[0] < 0)) notCellAvg = false; 
  double magB_sq_inv[3] = {0.0}; 

  if (notCellAvg) { 
  ser_1x_p2_inv(magB_sq, magB_sq_inv); 
  ExB_x[0] = 0.7071067811865475*magB_sq_inv[2]*num_ExB_x[2]+0.7071067811865475*magB_sq_inv[1]*num_ExB_x[1]+0.7071067811865475*magB_sq_inv[0]*num_ExB_x[0]; 
  ExB_x[1] = 0.6324555320336759*magB_sq_inv[1]*num_ExB_x[2]+0.6324555320336759*num_ExB_x[1]*magB_sq_inv[2]+0.7071067811865475*magB_sq_inv[0]*num_ExB_x[1]+0.7071067811865475*num_ExB_x[0]*magB_sq_inv[1]; 
  ExB_x[2] = 0.4517539514526256*magB_sq_inv[2]*num_ExB_x[2]+0.7071067811865475*magB_sq_inv[0]*num_ExB_x[2]+0.7071067811865475*num_ExB_x[0]*magB_sq_inv[2]+0.6324555320336759*magB_sq_inv[1]*num_ExB_x[1]; 

  ExB_y[0] = 0.7071067811865475*magB_sq_inv[2]*num_ExB_y[2]+0.7071067811865475*magB_sq_inv[1]*num_ExB_y[1]+0.7071067811865475*magB_sq_inv[0]*num_ExB_y[0]; 
  ExB_y[1] = 0.6324555320336759*magB_sq_inv[1]*num_ExB_y[2]+0.6324555320336759*num_ExB_y[1]*magB_sq_inv[2]+0.7071067811865475*magB_sq_inv[0]*num_ExB_y[1]+0.7071067811865475*num_ExB_y[0]*magB_sq_inv[1]; 
  ExB_y[2] = 0.4517539514526256*magB_sq_inv[2]*num_ExB_y[2]+0.7071067811865475*magB_sq_inv[0]*num_ExB_y[2]+0.7071067811865475*num_ExB_y[0]*magB_sq_inv[2]+0.6324555320336759*magB_sq_inv[1]*num_ExB_y[1]; 

  ExB_z[0] = 0.7071067811865475*magB_sq_inv[2]*num_ExB_z[2]+0.7071067811865475*magB_sq_inv[1]*num_ExB_z[1]+0.7071067811865475*magB_sq_inv[0]*num_ExB_z[0]; 
  ExB_z[1] = 0.6324555320336759*magB_sq_inv[1]*num_ExB_z[2]+0.6324555320336759*num_ExB_z[1]*magB_sq_inv[2]+0.7071067811865475*magB_sq_inv[0]*num_ExB_z[1]+0.7071067811865475*num_ExB_z[0]*magB_sq_inv[1]; 
  ExB_z[2] = 0.4517539514526256*magB_sq_inv[2]*num_ExB_z[2]+0.7071067811865475*magB_sq_inv[0]*num_ExB_z[2]+0.7071067811865475*num_ExB_z[0]*magB_sq_inv[2]+0.6324555320336759*magB_sq_inv[1]*num_ExB_z[1]; 

  } else { 
  magB_sq_inv[0] = 2.0/magB_sq[0]; 
  ExB_x[0] = 0.7071067811865475*magB_sq_inv[2]*num_ExB_x[2]+0.7071067811865475*magB_sq_inv[1]*num_ExB_x[1]+0.7071067811865475*magB_sq_inv[0]*num_ExB_x[0]; 
  ExB_y[0] = 0.7071067811865475*magB_sq_inv[2]*num_ExB_y[2]+0.7071067811865475*magB_sq_inv[1]*num_ExB_y[1]+0.7071067811865475*magB_sq_inv[0]*num_ExB_y[0]; 
  ExB_z[0] = 0.7071067811865475*magB_sq_inv[2]*num_ExB_z[2]+0.7071067811865475*magB_sq_inv[1]*num_ExB_z[1]+0.7071067811865475*magB_sq_inv[0]*num_ExB_z[0]; 

  } 
} 
 
