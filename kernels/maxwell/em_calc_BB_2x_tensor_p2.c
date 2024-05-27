#include <gkyl_maxwell_kernels.h> 
#include <gkyl_basis_tensor_2x_2p_exp_sq.h> 
GKYL_CU_DH void em_calc_BB_2x_tensor_p2(const double *em, double* GKYL_RESTRICT out) 
{ 
  // em:  Input electromagnetic fields. 
  // out: Output B_i B_j tensor. 
 
  double *B_x_sq  = &out[0]; 
  double *B_x_B_y = &out[9]; 
  double *B_x_B_z = &out[18]; 
  double *B_y_sq  = &out[27]; 
  double *B_y_B_z = &out[36]; 
  double *B_z_sq  = &out[45]; 
 
  const double *B_x = &em[12]; 
  const double *B_y = &em[16]; 
  const double *B_z = &em[20]; 
 
  // Calculate B_i B_i. 
  tensor_2x_2p_exp_sq(B_x, B_x_sq); 
  tensor_2x_2p_exp_sq(B_y, B_y_sq); 
  tensor_2x_2p_exp_sq(B_z, B_z_sq); 
 
  // Calculate B_i B_j. 
  B_x_B_y[0] = 0.5*B_x[3]*B_y[3]+0.5*B_x[2]*B_y[2]+0.5*B_x[1]*B_y[1]+0.5*B_x[0]*B_y[0]; 
  B_x_B_y[1] = 0.5*B_x[2]*B_y[3]+0.5*B_y[2]*B_x[3]+0.5*B_x[0]*B_y[1]+0.5*B_y[0]*B_x[1]; 
  B_x_B_y[2] = 0.5*B_x[1]*B_y[3]+0.5*B_y[1]*B_x[3]+0.5*B_x[0]*B_y[2]+0.5*B_y[0]*B_x[2]; 
  B_x_B_y[3] = 0.5*B_x[0]*B_y[3]+0.5*B_y[0]*B_x[3]+0.5*B_x[1]*B_y[2]+0.5*B_y[1]*B_x[2]; 
  B_x_B_y[4] = 0.4472135954999579*B_x[3]*B_y[3]+0.4472135954999579*B_x[1]*B_y[1]; 
  B_x_B_y[5] = 0.4472135954999579*B_x[3]*B_y[3]+0.4472135954999579*B_x[2]*B_y[2]; 
  B_x_B_y[6] = 0.447213595499958*B_x[1]*B_y[3]+0.447213595499958*B_y[1]*B_x[3]; 
  B_x_B_y[7] = 0.447213595499958*B_x[2]*B_y[3]+0.447213595499958*B_y[2]*B_x[3]; 
  B_x_B_y[8] = 0.4*B_x[3]*B_y[3]; 
  B_x_B_z[0] = 0.5*B_x[3]*B_z[3]+0.5*B_x[2]*B_z[2]+0.5*B_x[1]*B_z[1]+0.5*B_x[0]*B_z[0]; 
  B_x_B_z[1] = 0.5*B_x[2]*B_z[3]+0.5*B_z[2]*B_x[3]+0.5*B_x[0]*B_z[1]+0.5*B_z[0]*B_x[1]; 
  B_x_B_z[2] = 0.5*B_x[1]*B_z[3]+0.5*B_z[1]*B_x[3]+0.5*B_x[0]*B_z[2]+0.5*B_z[0]*B_x[2]; 
  B_x_B_z[3] = 0.5*B_x[0]*B_z[3]+0.5*B_z[0]*B_x[3]+0.5*B_x[1]*B_z[2]+0.5*B_z[1]*B_x[2]; 
  B_x_B_z[4] = 0.4472135954999579*B_x[3]*B_z[3]+0.4472135954999579*B_x[1]*B_z[1]; 
  B_x_B_z[5] = 0.4472135954999579*B_x[3]*B_z[3]+0.4472135954999579*B_x[2]*B_z[2]; 
  B_x_B_z[6] = 0.447213595499958*B_x[1]*B_z[3]+0.447213595499958*B_z[1]*B_x[3]; 
  B_x_B_z[7] = 0.447213595499958*B_x[2]*B_z[3]+0.447213595499958*B_z[2]*B_x[3]; 
  B_x_B_z[8] = 0.4*B_x[3]*B_z[3]; 
  B_y_B_z[0] = 0.5*B_y[3]*B_z[3]+0.5*B_y[2]*B_z[2]+0.5*B_y[1]*B_z[1]+0.5*B_y[0]*B_z[0]; 
  B_y_B_z[1] = 0.5*B_y[2]*B_z[3]+0.5*B_z[2]*B_y[3]+0.5*B_y[0]*B_z[1]+0.5*B_z[0]*B_y[1]; 
  B_y_B_z[2] = 0.5*B_y[1]*B_z[3]+0.5*B_z[1]*B_y[3]+0.5*B_y[0]*B_z[2]+0.5*B_z[0]*B_y[2]; 
  B_y_B_z[3] = 0.5*B_y[0]*B_z[3]+0.5*B_z[0]*B_y[3]+0.5*B_y[1]*B_z[2]+0.5*B_z[1]*B_y[2]; 
  B_y_B_z[4] = 0.4472135954999579*B_y[3]*B_z[3]+0.4472135954999579*B_y[1]*B_z[1]; 
  B_y_B_z[5] = 0.4472135954999579*B_y[3]*B_z[3]+0.4472135954999579*B_y[2]*B_z[2]; 
  B_y_B_z[6] = 0.447213595499958*B_y[1]*B_z[3]+0.447213595499958*B_z[1]*B_y[3]; 
  B_y_B_z[7] = 0.447213595499958*B_y[2]*B_z[3]+0.447213595499958*B_z[2]*B_y[3]; 
  B_y_B_z[8] = 0.4*B_y[3]*B_z[3]; 
 
} 
 
