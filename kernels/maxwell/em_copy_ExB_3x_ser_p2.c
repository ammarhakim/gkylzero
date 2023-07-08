#include <gkyl_mat.h> 
#include <gkyl_maxwell_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
GKYL_CU_DH void em_copy_ExB_3x_ser_p2(struct gkyl_mat *x, const double *em, int* cell_avg_magB2, double* ExB) 
{ 
  // x:   Input solution vector for 1/|B|^2. 
  // em:  Input electromagnetic fields. 
  // cell_avg_magB2: Input flag for cell average if 1/|B|^2 only used cell averages. 
  // ExB: E x B velocity = E x B/|B|^2. 
 
  double magB2_inv[20] = {0.0}; 

  magB2_inv[0] = gkyl_mat_get(x,0,0); 
  magB2_inv[1] = gkyl_mat_get(x,1,0); 
  magB2_inv[2] = gkyl_mat_get(x,2,0); 
  magB2_inv[3] = gkyl_mat_get(x,3,0); 
  magB2_inv[4] = gkyl_mat_get(x,4,0); 
  magB2_inv[5] = gkyl_mat_get(x,5,0); 
  magB2_inv[6] = gkyl_mat_get(x,6,0); 
  magB2_inv[7] = gkyl_mat_get(x,7,0); 
  magB2_inv[8] = gkyl_mat_get(x,8,0); 
  magB2_inv[9] = gkyl_mat_get(x,9,0); 
  magB2_inv[10] = gkyl_mat_get(x,10,0); 
  magB2_inv[11] = gkyl_mat_get(x,11,0); 
  magB2_inv[12] = gkyl_mat_get(x,12,0); 
  magB2_inv[13] = gkyl_mat_get(x,13,0); 
  magB2_inv[14] = gkyl_mat_get(x,14,0); 
  magB2_inv[15] = gkyl_mat_get(x,15,0); 
  magB2_inv[16] = gkyl_mat_get(x,16,0); 
  magB2_inv[17] = gkyl_mat_get(x,17,0); 
  magB2_inv[18] = gkyl_mat_get(x,18,0); 
  magB2_inv[19] = gkyl_mat_get(x,19,0); 

  const double *E_x = &em[0]; 
  const double *E_y = &em[20]; 
  const double *E_z = &em[40]; 
  const double *B_x = &em[60]; 
  const double *B_y = &em[80]; 
  const double *B_z = &em[100]; 
 
  double *ExB_x = &ExB[0]; 
  double *ExB_y = &ExB[20]; 
  double *ExB_z = &ExB[40]; 
 
  // Calculate E_i B_j. 
  double E_x_B_y[20] = {0.0}; 
  binop_mul_3d_ser_p2(E_x, B_y, E_x_B_y); 
 
  double E_x_B_z[20] = {0.0}; 
  binop_mul_3d_ser_p2(E_x, B_z, E_x_B_z); 
 
  double E_y_B_x[20] = {0.0}; 
  binop_mul_3d_ser_p2(E_y, B_x, E_y_B_x); 
 
  double E_y_B_z[20] = {0.0}; 
  binop_mul_3d_ser_p2(E_y, B_z, E_y_B_z); 
 
  double E_z_B_x[20] = {0.0}; 
  binop_mul_3d_ser_p2(E_z, B_x, E_z_B_x); 
 
  double E_z_B_y[20] = {0.0}; 
  binop_mul_3d_ser_p2(E_z, B_y, E_z_B_y); 
 
  double num_ExB_x[20] = {0.0}; 
  double num_ExB_y[20] = {0.0}; 
  double num_ExB_z[20] = {0.0}; 
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
  num_ExB_x[4] = E_y_B_z[4] - E_z_B_y[4]; 
  num_ExB_y[4] = E_z_B_x[4] - E_x_B_z[4]; 
  num_ExB_z[4] = E_x_B_y[4] - E_y_B_x[4]; 
  num_ExB_x[5] = E_y_B_z[5] - E_z_B_y[5]; 
  num_ExB_y[5] = E_z_B_x[5] - E_x_B_z[5]; 
  num_ExB_z[5] = E_x_B_y[5] - E_y_B_x[5]; 
  num_ExB_x[6] = E_y_B_z[6] - E_z_B_y[6]; 
  num_ExB_y[6] = E_z_B_x[6] - E_x_B_z[6]; 
  num_ExB_z[6] = E_x_B_y[6] - E_y_B_x[6]; 
  num_ExB_x[7] = E_y_B_z[7] - E_z_B_y[7]; 
  num_ExB_y[7] = E_z_B_x[7] - E_x_B_z[7]; 
  num_ExB_z[7] = E_x_B_y[7] - E_y_B_x[7]; 
  num_ExB_x[8] = E_y_B_z[8] - E_z_B_y[8]; 
  num_ExB_y[8] = E_z_B_x[8] - E_x_B_z[8]; 
  num_ExB_z[8] = E_x_B_y[8] - E_y_B_x[8]; 
  num_ExB_x[9] = E_y_B_z[9] - E_z_B_y[9]; 
  num_ExB_y[9] = E_z_B_x[9] - E_x_B_z[9]; 
  num_ExB_z[9] = E_x_B_y[9] - E_y_B_x[9]; 
  num_ExB_x[10] = E_y_B_z[10] - E_z_B_y[10]; 
  num_ExB_y[10] = E_z_B_x[10] - E_x_B_z[10]; 
  num_ExB_z[10] = E_x_B_y[10] - E_y_B_x[10]; 
  num_ExB_x[11] = E_y_B_z[11] - E_z_B_y[11]; 
  num_ExB_y[11] = E_z_B_x[11] - E_x_B_z[11]; 
  num_ExB_z[11] = E_x_B_y[11] - E_y_B_x[11]; 
  num_ExB_x[12] = E_y_B_z[12] - E_z_B_y[12]; 
  num_ExB_y[12] = E_z_B_x[12] - E_x_B_z[12]; 
  num_ExB_z[12] = E_x_B_y[12] - E_y_B_x[12]; 
  num_ExB_x[13] = E_y_B_z[13] - E_z_B_y[13]; 
  num_ExB_y[13] = E_z_B_x[13] - E_x_B_z[13]; 
  num_ExB_z[13] = E_x_B_y[13] - E_y_B_x[13]; 
  num_ExB_x[14] = E_y_B_z[14] - E_z_B_y[14]; 
  num_ExB_y[14] = E_z_B_x[14] - E_x_B_z[14]; 
  num_ExB_z[14] = E_x_B_y[14] - E_y_B_x[14]; 
  num_ExB_x[15] = E_y_B_z[15] - E_z_B_y[15]; 
  num_ExB_y[15] = E_z_B_x[15] - E_x_B_z[15]; 
  num_ExB_z[15] = E_x_B_y[15] - E_y_B_x[15]; 
  num_ExB_x[16] = E_y_B_z[16] - E_z_B_y[16]; 
  num_ExB_y[16] = E_z_B_x[16] - E_x_B_z[16]; 
  num_ExB_z[16] = E_x_B_y[16] - E_y_B_x[16]; 
  num_ExB_x[17] = E_y_B_z[17] - E_z_B_y[17]; 
  num_ExB_y[17] = E_z_B_x[17] - E_x_B_z[17]; 
  num_ExB_z[17] = E_x_B_y[17] - E_y_B_x[17]; 
  num_ExB_x[18] = E_y_B_z[18] - E_z_B_y[18]; 
  num_ExB_y[18] = E_z_B_x[18] - E_x_B_z[18]; 
  num_ExB_z[18] = E_x_B_y[18] - E_y_B_x[18]; 
  num_ExB_x[19] = E_y_B_z[19] - E_z_B_y[19]; 
  num_ExB_y[19] = E_z_B_x[19] - E_x_B_z[19]; 
  num_ExB_z[19] = E_x_B_y[19] - E_y_B_x[19]; 
 
  if (cell_avg_magB2[0]) { 
    num_ExB_x[1] = 0.0; 
    num_ExB_y[1] = 0.0; 
    num_ExB_z[1] = 0.0; 
    num_ExB_x[2] = 0.0; 
    num_ExB_y[2] = 0.0; 
    num_ExB_z[2] = 0.0; 
    num_ExB_x[3] = 0.0; 
    num_ExB_y[3] = 0.0; 
    num_ExB_z[3] = 0.0; 
    num_ExB_x[4] = 0.0; 
    num_ExB_y[4] = 0.0; 
    num_ExB_z[4] = 0.0; 
    num_ExB_x[5] = 0.0; 
    num_ExB_y[5] = 0.0; 
    num_ExB_z[5] = 0.0; 
    num_ExB_x[6] = 0.0; 
    num_ExB_y[6] = 0.0; 
    num_ExB_z[6] = 0.0; 
    num_ExB_x[7] = 0.0; 
    num_ExB_y[7] = 0.0; 
    num_ExB_z[7] = 0.0; 
    num_ExB_x[8] = 0.0; 
    num_ExB_y[8] = 0.0; 
    num_ExB_z[8] = 0.0; 
    num_ExB_x[9] = 0.0; 
    num_ExB_y[9] = 0.0; 
    num_ExB_z[9] = 0.0; 
    num_ExB_x[10] = 0.0; 
    num_ExB_y[10] = 0.0; 
    num_ExB_z[10] = 0.0; 
    num_ExB_x[11] = 0.0; 
    num_ExB_y[11] = 0.0; 
    num_ExB_z[11] = 0.0; 
    num_ExB_x[12] = 0.0; 
    num_ExB_y[12] = 0.0; 
    num_ExB_z[12] = 0.0; 
    num_ExB_x[13] = 0.0; 
    num_ExB_y[13] = 0.0; 
    num_ExB_z[13] = 0.0; 
    num_ExB_x[14] = 0.0; 
    num_ExB_y[14] = 0.0; 
    num_ExB_z[14] = 0.0; 
    num_ExB_x[15] = 0.0; 
    num_ExB_y[15] = 0.0; 
    num_ExB_z[15] = 0.0; 
    num_ExB_x[16] = 0.0; 
    num_ExB_y[16] = 0.0; 
    num_ExB_z[16] = 0.0; 
    num_ExB_x[17] = 0.0; 
    num_ExB_y[17] = 0.0; 
    num_ExB_z[17] = 0.0; 
    num_ExB_x[18] = 0.0; 
    num_ExB_y[18] = 0.0; 
    num_ExB_z[18] = 0.0; 
    num_ExB_x[19] = 0.0; 
    num_ExB_y[19] = 0.0; 
    num_ExB_z[19] = 0.0; 
  } 
  // Calculate expansions of E x B|B|^2, which can be calculated free of aliasing errors. 
  binop_mul_3d_ser_p2(magB2_inv, num_ExB_x, ExB_x); 
  binop_mul_3d_ser_p2(magB2_inv, num_ExB_y, ExB_y); 
  binop_mul_3d_ser_p2(magB2_inv, num_ExB_z, ExB_z); 
 
} 
 
