#include <gkyl_mat.h> 
#include <gkyl_maxwell_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
GKYL_CU_DH void em_copy_ExB_1x_ser_p1(struct gkyl_mat *x, const double *em, int* cell_avg_magB2, double* ExB) 
{ 
  // x:   Input solution vector for 1/|B|^2. 
  // em:  Input electromagnetic fields. 
  // cell_avg_magB2: Input flag for cell average if 1/|B|^2 only used cell averages. 
  // ExB: E x B velocity = E x B/|B|^2. 
 
  double magB2_inv[2] = {0.0}; 

  magB2_inv[0] = gkyl_mat_get(x,0,0); 
  magB2_inv[1] = gkyl_mat_get(x,1,0); 

  const double *E_x = &em[0]; 
  const double *E_y = &em[2]; 
  const double *E_z = &em[4]; 
  const double *B_x = &em[6]; 
  const double *B_y = &em[8]; 
  const double *B_z = &em[10]; 
 
  double *ExB_x = &ExB[0]; 
  double *ExB_y = &ExB[2]; 
  double *ExB_z = &ExB[4]; 
 
  // Calculate E_i B_j. 
  double E_x_B_y[2] = {0.0}; 
  binop_mul_1d_ser_p1(E_x, B_y, E_x_B_y); 
 
  double E_x_B_z[2] = {0.0}; 
  binop_mul_1d_ser_p1(E_x, B_z, E_x_B_z); 
 
  double E_y_B_x[2] = {0.0}; 
  binop_mul_1d_ser_p1(E_y, B_x, E_y_B_x); 
 
  double E_y_B_z[2] = {0.0}; 
  binop_mul_1d_ser_p1(E_y, B_z, E_y_B_z); 
 
  double E_z_B_x[2] = {0.0}; 
  binop_mul_1d_ser_p1(E_z, B_x, E_z_B_x); 
 
  double E_z_B_y[2] = {0.0}; 
  binop_mul_1d_ser_p1(E_z, B_y, E_z_B_y); 
 
  double num_ExB_x[2] = {0.0}; 
  double num_ExB_y[2] = {0.0}; 
  double num_ExB_z[2] = {0.0}; 
  num_ExB_x[0] = E_y_B_z[0] - E_z_B_y[0]; 
  num_ExB_y[0] = E_z_B_x[0] - E_x_B_z[0]; 
  num_ExB_z[0] = E_x_B_y[0] - E_y_B_x[0]; 
  num_ExB_x[1] = E_y_B_z[1] - E_z_B_y[1]; 
  num_ExB_y[1] = E_z_B_x[1] - E_x_B_z[1]; 
  num_ExB_z[1] = E_x_B_y[1] - E_y_B_x[1]; 
 
  if (cell_avg_magB2[0]) { 
    num_ExB_x[1] = 0.0; 
    num_ExB_y[1] = 0.0; 
    num_ExB_z[1] = 0.0; 
  } 
  // Calculate expansions of E x B|B|^2, which can be calculated free of aliasing errors. 
  binop_mul_1d_ser_p1(magB2_inv, num_ExB_x, ExB_x); 
  binop_mul_1d_ser_p1(magB2_inv, num_ExB_y, ExB_y); 
  binop_mul_1d_ser_p1(magB2_inv, num_ExB_z, ExB_z); 
 
} 
 
