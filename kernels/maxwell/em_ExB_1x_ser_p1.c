#include <gkyl_maxwell_kernels.h> 
#include <gkyl_basis_ser_1x_p1_exp_sq.h> 
#include <gkyl_basis_ser_1x_p1_inv.h> 
GKYL_CU_DH void em_ExB_1x_ser_p1(const double *em, double* GKYL_RESTRICT ExB) 
{ 
  // em:  Input electromagnetic fields. 
  // ExB: E x B velocity = E x B/|B|^2. 
 
  const double *E_x = &em[0]; 
  const double *E_y = &em[2]; 
  const double *E_z = &em[4]; 
  const double *B_x = &em[6]; 
  const double *B_y = &em[8]; 
  const double *B_z = &em[10]; 
 
  double *ExB_x = &ExB[0]; 
  double *ExB_y = &ExB[2]; 
  double *ExB_z = &ExB[4]; 
 
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
  ExB_x[0] = (-0.5*B_y[0]*E_z[1]*magB_sq_inv[1])+0.5*B_z[0]*E_y[1]*magB_sq_inv[1]+0.5*E_y[0]*B_z[1]*magB_sq_inv[1]-0.5*E_z[0]*B_y[1]*magB_sq_inv[1]-0.5*magB_sq_inv[0]*B_y[1]*E_z[1]+0.5*magB_sq_inv[0]*B_z[1]*E_y[1]-0.5*B_y[0]*E_z[0]*magB_sq_inv[0]+0.5*B_z[0]*E_y[0]*magB_sq_inv[0]; 
  ExB_x[1] = (-0.9*B_y[1]*E_z[1]*magB_sq_inv[1])+0.9*B_z[1]*E_y[1]*magB_sq_inv[1]-0.5*B_y[0]*E_z[0]*magB_sq_inv[1]+0.5*B_z[0]*E_y[0]*magB_sq_inv[1]-0.5*B_y[0]*magB_sq_inv[0]*E_z[1]+0.5*B_z[0]*magB_sq_inv[0]*E_y[1]+0.5*E_y[0]*magB_sq_inv[0]*B_z[1]-0.5*E_z[0]*magB_sq_inv[0]*B_y[1]; 

  ExB_y[0] = 0.5*B_x[0]*E_z[1]*magB_sq_inv[1]-0.5*B_z[0]*E_x[1]*magB_sq_inv[1]-0.5*E_x[0]*B_z[1]*magB_sq_inv[1]+0.5*E_z[0]*B_x[1]*magB_sq_inv[1]+0.5*magB_sq_inv[0]*B_x[1]*E_z[1]-0.5*magB_sq_inv[0]*B_z[1]*E_x[1]+0.5*B_x[0]*E_z[0]*magB_sq_inv[0]-0.5*B_z[0]*E_x[0]*magB_sq_inv[0]; 
  ExB_y[1] = 0.9*B_x[1]*E_z[1]*magB_sq_inv[1]-0.9*B_z[1]*E_x[1]*magB_sq_inv[1]+0.5*B_x[0]*E_z[0]*magB_sq_inv[1]-0.5*B_z[0]*E_x[0]*magB_sq_inv[1]+0.5*B_x[0]*magB_sq_inv[0]*E_z[1]-0.5*B_z[0]*magB_sq_inv[0]*E_x[1]-0.5*E_x[0]*magB_sq_inv[0]*B_z[1]+0.5*E_z[0]*magB_sq_inv[0]*B_x[1]; 

  ExB_z[0] = (-0.5*B_x[0]*E_y[1]*magB_sq_inv[1])+0.5*B_y[0]*E_x[1]*magB_sq_inv[1]+0.5*E_x[0]*B_y[1]*magB_sq_inv[1]-0.5*E_y[0]*B_x[1]*magB_sq_inv[1]-0.5*magB_sq_inv[0]*B_x[1]*E_y[1]+0.5*magB_sq_inv[0]*B_y[1]*E_x[1]-0.5*B_x[0]*E_y[0]*magB_sq_inv[0]+0.5*B_y[0]*E_x[0]*magB_sq_inv[0]; 
  ExB_z[1] = (-0.9*B_x[1]*E_y[1]*magB_sq_inv[1])+0.9*B_y[1]*E_x[1]*magB_sq_inv[1]-0.5*B_x[0]*E_y[0]*magB_sq_inv[1]+0.5*B_y[0]*E_x[0]*magB_sq_inv[1]-0.5*B_x[0]*magB_sq_inv[0]*E_y[1]+0.5*B_y[0]*magB_sq_inv[0]*E_x[1]+0.5*E_x[0]*magB_sq_inv[0]*B_y[1]-0.5*E_y[0]*magB_sq_inv[0]*B_x[1]; 

} 
 
