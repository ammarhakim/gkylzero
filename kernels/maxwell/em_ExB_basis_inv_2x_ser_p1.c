#include <gkyl_maxwell_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
#include <gkyl_basis_ser_2x_p1_inv.h> 
GKYL_CU_DH void em_ExB_basis_inv_2x_ser_p1(const double *em, double* ExB) 
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
 
  // Calculate B_x^2, B_y^2, and B_z^2. 
  double B_x_sq[4] = {0.0}; 
  binop_mul_2d_ser_p1(B_x, B_x, B_x_sq); 
 
  double B_y_sq[4] = {0.0}; 
  binop_mul_2d_ser_p1(B_y, B_y, B_y_sq); 
 
  double B_z_sq[4] = {0.0}; 
  binop_mul_2d_ser_p1(B_z, B_z, B_z_sq); 
 
  double magB2[4] = {0.0}; 

  magB2[0] = B_z_sq[0]+B_y_sq[0]+B_x_sq[0]; 
  magB2[1] = B_z_sq[1]+B_y_sq[1]+B_x_sq[1]; 
  magB2[2] = B_z_sq[2]+B_y_sq[2]+B_x_sq[2]; 
  magB2[3] = B_z_sq[3]+B_y_sq[3]+B_x_sq[3]; 

  int cell_avg = 0;
  // Check if |B|^2 < 0 at control points. 
  if (1.5*magB2[3]-0.8660254037844386*magB2[2]-0.8660254037844386*magB2[1]+0.5*magB2[0] < 0.0) cell_avg = 1; 
  if ((-1.5*magB2[3])-0.8660254037844386*magB2[2]+0.8660254037844386*magB2[1]+0.5*magB2[0] < 0.0) cell_avg = 1; 
  if ((-1.5*magB2[3])+0.8660254037844386*magB2[2]-0.8660254037844386*magB2[1]+0.5*magB2[0] < 0.0) cell_avg = 1; 
  if (1.5*magB2[3]+0.8660254037844386*magB2[2]+0.8660254037844386*magB2[1]+0.5*magB2[0] < 0.0) cell_avg = 1; 
  double magB2_inv[4] = {0.0}; 

  if (cell_avg) { 
  // If |B|^2 < 0 at control points, only use cell average to get 1/|B|^2. 
  magB2_inv[0] = 4.0/magB2[0]; 
  } else { 
  ser_2x_p1_inv(magB2, magB2_inv); 
  } 
  // Calculate E_i B_j. 
  double E_x_B_y[4] = {0.0}; 
  binop_mul_2d_ser_p1(E_x, B_y, E_x_B_y); 
 
  double E_x_B_z[4] = {0.0}; 
  binop_mul_2d_ser_p1(E_x, B_z, E_x_B_z); 
 
  double E_y_B_x[4] = {0.0}; 
  binop_mul_2d_ser_p1(E_y, B_x, E_y_B_x); 
 
  double E_y_B_z[4] = {0.0}; 
  binop_mul_2d_ser_p1(E_y, B_z, E_y_B_z); 
 
  double E_z_B_x[4] = {0.0}; 
  binop_mul_2d_ser_p1(E_z, B_x, E_z_B_x); 
 
  double E_z_B_y[4] = {0.0}; 
  binop_mul_2d_ser_p1(E_z, B_y, E_z_B_y); 
 
  double num_ExB_x[4] = {0.0}; 
  double num_ExB_y[4] = {0.0}; 
  double num_ExB_z[4] = {0.0}; 
  num_ExB_x[0] = E_y_B_z[0] - E_z_B_y[0]; 
  num_ExB_y[0] = E_z_B_x[0] - E_x_B_z[0]; 
  num_ExB_z[0] = E_x_B_y[0] - E_y_B_x[0]; 
  num_ExB_x[1] = E_y_B_z[1] - E_z_B_y[1]; 
  num_ExB_y[1] = E_z_B_x[1] - E_x_B_z[1]; 
  num_ExB_z[1] = E_x_B_y[1] - E_y_B_x[1]; 
  num_ExB_x[2] = E_y_B_z[2] - E_z_B_y[2]; 
  num_ExB_y[2] = E_z_B_x[2] - E_x_B_z[2]; 
  num_ExB_z[2] = E_x_B_y[2] - E_y_B_x[2]; 
  num_ExB_x[3] = E_y_B_z[3] - E_z_B_y[3]; 
  num_ExB_y[3] = E_z_B_x[3] - E_x_B_z[3]; 
  num_ExB_z[3] = E_x_B_y[3] - E_y_B_x[3]; 
 
  if (cell_avg) { 
    num_ExB_x[1] = 0.0; 
    num_ExB_y[1] = 0.0; 
    num_ExB_z[1] = 0.0; 
    num_ExB_x[2] = 0.0; 
    num_ExB_y[2] = 0.0; 
    num_ExB_z[2] = 0.0; 
    num_ExB_x[3] = 0.0; 
    num_ExB_y[3] = 0.0; 
    num_ExB_z[3] = 0.0; 
  } 
  // Calculate expansions of E x B|B|^2, which can be calculated free of aliasing errors. 
  binop_mul_2d_ser_p1(magB2_inv, num_ExB_x, ExB_x); 
  binop_mul_2d_ser_p1(magB2_inv, num_ExB_y, ExB_y); 
  binop_mul_2d_ser_p1(magB2_inv, num_ExB_z, ExB_z); 
 
} 
 
