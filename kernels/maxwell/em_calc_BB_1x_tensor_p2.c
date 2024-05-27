#include <gkyl_maxwell_kernels.h> 
#include <gkyl_basis_tensor_1x_2p_exp_sq.h> 
GKYL_CU_DH void em_calc_BB_1x_tensor_p2(const double *em, double* GKYL_RESTRICT out) 
{ 
  // em:  Input electromagnetic fields. 
  // out: Output B_i B_j tensor. 
 
  double *B_x_sq  = &out[0]; 
  double *B_x_B_y = &out[3]; 
  double *B_x_B_z = &out[6]; 
  double *B_y_sq  = &out[9]; 
  double *B_y_B_z = &out[12]; 
  double *B_z_sq  = &out[15]; 
 
  const double *B_x = &em[6]; 
  const double *B_y = &em[8]; 
  const double *B_z = &em[10]; 
 
  // Calculate B_i B_i. 
  tensor_1x_2p_exp_sq(B_x, B_x_sq); 
  tensor_1x_2p_exp_sq(B_y, B_y_sq); 
  tensor_1x_2p_exp_sq(B_z, B_z_sq); 
 
  // Calculate B_i B_j. 
  B_x_B_y[0] = 0.7071067811865475*B_x[1]*B_y[1]+0.7071067811865475*B_x[0]*B_y[0]; 
  B_x_B_y[1] = 0.7071067811865475*B_x[0]*B_y[1]+0.7071067811865475*B_y[0]*B_x[1]; 
  B_x_B_y[2] = 0.6324555320336759*B_x[1]*B_y[1]; 
  B_x_B_z[0] = 0.7071067811865475*B_x[1]*B_z[1]+0.7071067811865475*B_x[0]*B_z[0]; 
  B_x_B_z[1] = 0.7071067811865475*B_x[0]*B_z[1]+0.7071067811865475*B_z[0]*B_x[1]; 
  B_x_B_z[2] = 0.6324555320336759*B_x[1]*B_z[1]; 
  B_y_B_z[0] = 0.7071067811865475*B_y[1]*B_z[1]+0.7071067811865475*B_y[0]*B_z[0]; 
  B_y_B_z[1] = 0.7071067811865475*B_y[0]*B_z[1]+0.7071067811865475*B_z[0]*B_y[1]; 
  B_y_B_z[2] = 0.6324555320336759*B_y[1]*B_z[1]; 
 
} 
 
