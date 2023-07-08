#include <gkyl_mat.h> 
#include <gkyl_maxwell_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
GKYL_CU_DH int em_set_magB2_1x_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *em) 
{ 
  // A:     preallocated LHS matrix. 
  // rhs:   preallocated RHS vector. 
  // em:    Input electromagnetic fields. 
 
  const double *B_x = &em[9]; 
  const double *B_y = &em[12]; 
  const double *B_z = &em[15]; 
 
  // Calculate |B|^2 and set matrix to solve for 1/|B|^2. 
  double B_x_sq[3] = {0.0}; 
  binop_mul_1d_ser_p2(B_x, B_x, B_x_sq); 
 
  double B_y_sq[3] = {0.0}; 
  binop_mul_1d_ser_p2(B_y, B_y, B_y_sq); 
 
  double B_z_sq[3] = {0.0}; 
  binop_mul_1d_ser_p2(B_z, B_z, B_z_sq); 
 
  double magB2[3] = {0.0}; 

  magB2[0] = B_z_sq[0]+B_y_sq[0]+B_x_sq[0]; 
  magB2[1] = B_z_sq[1]+B_y_sq[1]+B_x_sq[1]; 
  magB2[2] = B_z_sq[2]+B_y_sq[2]+B_x_sq[2]; 

  int cell_avg = 0;
  // Check if |B|^2 < 0 at control points. 
  if (1.58113883008419*magB2[2]-1.224744871391589*magB2[1]+0.7071067811865475*magB2[0] < 0.0) cell_avg = 1; 
  if (0.7071067811865475*magB2[0]-0.7905694150420947*magB2[2] < 0.0) cell_avg = 1; 
  if (1.58113883008419*magB2[2]+1.224744871391589*magB2[1]+0.7071067811865475*magB2[0] < 0.0) cell_avg = 1; 
  if (cell_avg) { 
    magB2[1] = 0.0; 
    magB2[2] = 0.0; 
  } 
  gkyl_mat_set(rhs,0,0,1.414213562373095); 
  gkyl_mat_set(A,0,0,0.7071067811865475*magB2[0]); 
  gkyl_mat_set(A,0,1,0.7071067811865475*magB2[1]); 
  gkyl_mat_set(A,0,2,0.7071067811865475*magB2[2]); 
  gkyl_mat_set(A,1,0,0.7071067811865475*magB2[1]); 
  gkyl_mat_set(A,1,1,0.6324555320336759*magB2[2]+0.7071067811865475*magB2[0]); 
  gkyl_mat_set(A,1,2,0.6324555320336759*magB2[1]); 
  gkyl_mat_set(A,2,0,0.7071067811865475*magB2[2]); 
  gkyl_mat_set(A,2,1,0.6324555320336759*magB2[1]); 
  gkyl_mat_set(A,2,2,0.4517539514526256*magB2[2]+0.7071067811865475*magB2[0]); 
  return cell_avg;
} 
