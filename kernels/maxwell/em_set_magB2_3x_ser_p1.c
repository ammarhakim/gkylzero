#include <gkyl_mat.h> 
#include <gkyl_maxwell_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
GKYL_CU_DH int em_set_magB2_3x_ser_p1(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *em) 
{ 
  // A:     preallocated LHS matrix. 
  // rhs:   preallocated RHS vector. 
  // em:    Input electromagnetic fields. 
 
  const double *B_x = &em[24]; 
  const double *B_y = &em[32]; 
  const double *B_z = &em[40]; 
 
  // Calculate |B|^2 and set matrix to solve for 1/|B|^2. 
  double B_x_sq[8] = {0.0}; 
  binop_mul_3d_ser_p1(B_x, B_x, B_x_sq); 
 
  double B_y_sq[8] = {0.0}; 
  binop_mul_3d_ser_p1(B_y, B_y, B_y_sq); 
 
  double B_z_sq[8] = {0.0}; 
  binop_mul_3d_ser_p1(B_z, B_z, B_z_sq); 
 
  double magB2[8] = {0.0}; 

  magB2[0] = B_z_sq[0]+B_y_sq[0]+B_x_sq[0]; 
  magB2[1] = B_z_sq[1]+B_y_sq[1]+B_x_sq[1]; 
  magB2[2] = B_z_sq[2]+B_y_sq[2]+B_x_sq[2]; 
  magB2[3] = B_z_sq[3]+B_y_sq[3]+B_x_sq[3]; 
  magB2[4] = B_z_sq[4]+B_y_sq[4]+B_x_sq[4]; 
  magB2[5] = B_z_sq[5]+B_y_sq[5]+B_x_sq[5]; 
  magB2[6] = B_z_sq[6]+B_y_sq[6]+B_x_sq[6]; 
  magB2[7] = B_z_sq[7]+B_y_sq[7]+B_x_sq[7]; 

  int cell_avg = 0;
  // Check if |B|^2 < 0 at control points. 
  if ((-1.837117307087383*magB2[7])+1.060660171779821*magB2[6]+1.060660171779821*magB2[5]+1.060660171779821*magB2[4]-0.6123724356957944*magB2[3]-0.6123724356957944*magB2[2]-0.6123724356957944*magB2[1]+0.3535533905932737*magB2[0] < 0.0) cell_avg = 1; 
  if (1.837117307087383*magB2[7]+1.060660171779821*magB2[6]-1.060660171779821*magB2[5]-1.060660171779821*magB2[4]-0.6123724356957944*magB2[3]-0.6123724356957944*magB2[2]+0.6123724356957944*magB2[1]+0.3535533905932737*magB2[0] < 0.0) cell_avg = 1; 
  if (1.837117307087383*magB2[7]-1.060660171779821*magB2[6]+1.060660171779821*magB2[5]-1.060660171779821*magB2[4]-0.6123724356957944*magB2[3]+0.6123724356957944*magB2[2]-0.6123724356957944*magB2[1]+0.3535533905932737*magB2[0] < 0.0) cell_avg = 1; 
  if ((-1.837117307087383*magB2[7])-1.060660171779821*magB2[6]-1.060660171779821*magB2[5]+1.060660171779821*magB2[4]-0.6123724356957944*magB2[3]+0.6123724356957944*magB2[2]+0.6123724356957944*magB2[1]+0.3535533905932737*magB2[0] < 0.0) cell_avg = 1; 
  if (1.837117307087383*magB2[7]-1.060660171779821*magB2[6]-1.060660171779821*magB2[5]+1.060660171779821*magB2[4]+0.6123724356957944*magB2[3]-0.6123724356957944*magB2[2]-0.6123724356957944*magB2[1]+0.3535533905932737*magB2[0] < 0.0) cell_avg = 1; 
  if ((-1.837117307087383*magB2[7])-1.060660171779821*magB2[6]+1.060660171779821*magB2[5]-1.060660171779821*magB2[4]+0.6123724356957944*magB2[3]-0.6123724356957944*magB2[2]+0.6123724356957944*magB2[1]+0.3535533905932737*magB2[0] < 0.0) cell_avg = 1; 
  if ((-1.837117307087383*magB2[7])+1.060660171779821*magB2[6]-1.060660171779821*magB2[5]-1.060660171779821*magB2[4]+0.6123724356957944*magB2[3]+0.6123724356957944*magB2[2]-0.6123724356957944*magB2[1]+0.3535533905932737*magB2[0] < 0.0) cell_avg = 1; 
  if (1.837117307087383*magB2[7]+1.060660171779821*magB2[6]+1.060660171779821*magB2[5]+1.060660171779821*magB2[4]+0.6123724356957944*magB2[3]+0.6123724356957944*magB2[2]+0.6123724356957944*magB2[1]+0.3535533905932737*magB2[0] < 0.0) cell_avg = 1; 
  if (cell_avg) { 
    magB2[1] = 0.0; 
    magB2[2] = 0.0; 
    magB2[3] = 0.0; 
    magB2[4] = 0.0; 
    magB2[5] = 0.0; 
    magB2[6] = 0.0; 
    magB2[7] = 0.0; 
  } 
  gkyl_mat_set(rhs,0,0,2.828427124746191); 
  gkyl_mat_set(A,0,0,0.3535533905932737*magB2[0]); 
  gkyl_mat_set(A,0,1,0.3535533905932737*magB2[1]); 
  gkyl_mat_set(A,0,2,0.3535533905932737*magB2[2]); 
  gkyl_mat_set(A,0,3,0.3535533905932737*magB2[3]); 
  gkyl_mat_set(A,0,4,0.3535533905932737*magB2[4]); 
  gkyl_mat_set(A,0,5,0.3535533905932737*magB2[5]); 
  gkyl_mat_set(A,0,6,0.3535533905932737*magB2[6]); 
  gkyl_mat_set(A,0,7,0.3535533905932737*magB2[7]); 
  gkyl_mat_set(A,1,0,0.3535533905932737*magB2[1]); 
  gkyl_mat_set(A,1,1,0.3535533905932737*magB2[0]); 
  gkyl_mat_set(A,1,2,0.3535533905932737*magB2[4]); 
  gkyl_mat_set(A,1,3,0.3535533905932737*magB2[5]); 
  gkyl_mat_set(A,1,4,0.3535533905932737*magB2[2]); 
  gkyl_mat_set(A,1,5,0.3535533905932737*magB2[3]); 
  gkyl_mat_set(A,1,6,0.3535533905932737*magB2[7]); 
  gkyl_mat_set(A,1,7,0.3535533905932737*magB2[6]); 
  gkyl_mat_set(A,2,0,0.3535533905932737*magB2[2]); 
  gkyl_mat_set(A,2,1,0.3535533905932737*magB2[4]); 
  gkyl_mat_set(A,2,2,0.3535533905932737*magB2[0]); 
  gkyl_mat_set(A,2,3,0.3535533905932737*magB2[6]); 
  gkyl_mat_set(A,2,4,0.3535533905932737*magB2[1]); 
  gkyl_mat_set(A,2,5,0.3535533905932737*magB2[7]); 
  gkyl_mat_set(A,2,6,0.3535533905932737*magB2[3]); 
  gkyl_mat_set(A,2,7,0.3535533905932737*magB2[5]); 
  gkyl_mat_set(A,3,0,0.3535533905932737*magB2[3]); 
  gkyl_mat_set(A,3,1,0.3535533905932737*magB2[5]); 
  gkyl_mat_set(A,3,2,0.3535533905932737*magB2[6]); 
  gkyl_mat_set(A,3,3,0.3535533905932737*magB2[0]); 
  gkyl_mat_set(A,3,4,0.3535533905932737*magB2[7]); 
  gkyl_mat_set(A,3,5,0.3535533905932737*magB2[1]); 
  gkyl_mat_set(A,3,6,0.3535533905932737*magB2[2]); 
  gkyl_mat_set(A,3,7,0.3535533905932737*magB2[4]); 
  gkyl_mat_set(A,4,0,0.3535533905932737*magB2[4]); 
  gkyl_mat_set(A,4,1,0.3535533905932737*magB2[2]); 
  gkyl_mat_set(A,4,2,0.3535533905932737*magB2[1]); 
  gkyl_mat_set(A,4,3,0.3535533905932737*magB2[7]); 
  gkyl_mat_set(A,4,4,0.3535533905932737*magB2[0]); 
  gkyl_mat_set(A,4,5,0.3535533905932737*magB2[6]); 
  gkyl_mat_set(A,4,6,0.3535533905932737*magB2[5]); 
  gkyl_mat_set(A,4,7,0.3535533905932737*magB2[3]); 
  gkyl_mat_set(A,5,0,0.3535533905932737*magB2[5]); 
  gkyl_mat_set(A,5,1,0.3535533905932737*magB2[3]); 
  gkyl_mat_set(A,5,2,0.3535533905932737*magB2[7]); 
  gkyl_mat_set(A,5,3,0.3535533905932737*magB2[1]); 
  gkyl_mat_set(A,5,4,0.3535533905932737*magB2[6]); 
  gkyl_mat_set(A,5,5,0.3535533905932737*magB2[0]); 
  gkyl_mat_set(A,5,6,0.3535533905932737*magB2[4]); 
  gkyl_mat_set(A,5,7,0.3535533905932737*magB2[2]); 
  gkyl_mat_set(A,6,0,0.3535533905932737*magB2[6]); 
  gkyl_mat_set(A,6,1,0.3535533905932737*magB2[7]); 
  gkyl_mat_set(A,6,2,0.3535533905932737*magB2[3]); 
  gkyl_mat_set(A,6,3,0.3535533905932737*magB2[2]); 
  gkyl_mat_set(A,6,4,0.3535533905932737*magB2[5]); 
  gkyl_mat_set(A,6,5,0.3535533905932737*magB2[4]); 
  gkyl_mat_set(A,6,6,0.3535533905932737*magB2[0]); 
  gkyl_mat_set(A,6,7,0.3535533905932737*magB2[1]); 
  gkyl_mat_set(A,7,0,0.3535533905932737*magB2[7]); 
  gkyl_mat_set(A,7,1,0.3535533905932737*magB2[6]); 
  gkyl_mat_set(A,7,2,0.3535533905932737*magB2[5]); 
  gkyl_mat_set(A,7,3,0.3535533905932737*magB2[4]); 
  gkyl_mat_set(A,7,4,0.3535533905932737*magB2[3]); 
  gkyl_mat_set(A,7,5,0.3535533905932737*magB2[2]); 
  gkyl_mat_set(A,7,6,0.3535533905932737*magB2[1]); 
  gkyl_mat_set(A,7,7,0.3535533905932737*magB2[0]); 
  return cell_avg;
} 
