#include <gkyl_mat.h> 
#include <gkyl_maxwell_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
#include <gkyl_basis_ser_1x_p1_sqrt_with_sign.h> 
GKYL_CU_DH void em_copy_bvar_1x_ser_p1(struct gkyl_mat *x, const double *em, int* cell_avg_magB2, double* bvar) 
{ 
  // x:              Input solution vector for 1/|B|^2. 
  // em:             Input electromagnetic fields. 
  // cell_avg_magB2: Input flag for cell average if 1/|B|^2 only used cell averages. Also is adjusted if bb < 0 at control points. 
  // bvar:           Output b_i = B_i/|B| (first 3 components), b_i b_j = B_i B_j/|B|^2 (last 6 components). 
 
  double magB2_inv[2] = {0.0}; 

  magB2_inv[0] = gkyl_mat_get(x,0,0); 
  magB2_inv[1] = gkyl_mat_get(x,1,0); 

  double *bx = &bvar[0]; 
  double *by = &bvar[2]; 
  double *bz = &bvar[4]; 
  double *bxbx = &bvar[6]; 
  double *bxby = &bvar[8]; 
  double *bxbz = &bvar[10]; 
  double *byby = &bvar[12]; 
  double *bybz = &bvar[14]; 
  double *bzbz = &bvar[16]; 
 
  const double *B_x = &em[6]; 
  const double *B_y = &em[8]; 
  const double *B_z = &em[10]; 
 
  // Calculate B_i B_j. 
  double B_x_sq[2] = {0.0}; 
  binop_mul_1d_ser_p1(B_x, B_x, B_x_sq); 
 
  double B_y_sq[2] = {0.0}; 
  binop_mul_1d_ser_p1(B_y, B_y, B_y_sq); 
 
  double B_z_sq[2] = {0.0}; 
  binop_mul_1d_ser_p1(B_z, B_z, B_z_sq); 
 
  double B_x_B_y[2] = {0.0}; 
  binop_mul_1d_ser_p1(B_x, B_y, B_x_B_y); 
 
  double B_x_B_z[2] = {0.0}; 
  binop_mul_1d_ser_p1(B_x, B_z, B_x_B_z); 
 
  double B_y_B_z[2] = {0.0}; 
  binop_mul_1d_ser_p1(B_y, B_z, B_y_B_z); 
 
  // Calculate expansions of B_i B_j/|B|^2, which can be calculated free of aliasing errors. 
  binop_mul_1d_ser_p1(magB2_inv, B_x_sq, bxbx); 
  binop_mul_1d_ser_p1(magB2_inv, B_y_sq, byby); 
  binop_mul_1d_ser_p1(magB2_inv, B_z_sq, bzbz); 
  binop_mul_1d_ser_p1(magB2_inv, B_x_B_y, bxby); 
  binop_mul_1d_ser_p1(magB2_inv, B_x_B_z, bxbz); 
  binop_mul_1d_ser_p1(magB2_inv, B_y_B_z, bybz); 
 
  int cell_avg = 0;
  if (0.7071067811865475*bxbx[0]-0.7071067811865475*bxbx[1] < 0.0) cell_avg = 1; 
  if (0.7071067811865475*byby[0]-0.7071067811865475*byby[1] < 0.0) cell_avg = 1; 
  if (0.7071067811865475*bzbz[0]-0.7071067811865475*bzbz[1] < 0.0) cell_avg = 1; 
  if (0.7071067811865475*bxbx[1]+0.7071067811865475*bxbx[0] < 0.0) cell_avg = 1; 
  if (0.7071067811865475*byby[1]+0.7071067811865475*byby[0] < 0.0) cell_avg = 1; 
  if (0.7071067811865475*bzbz[1]+0.7071067811865475*bzbz[0] < 0.0) cell_avg = 1; 
  if (cell_avg || cell_avg_magB2[0]) { 
    bxbx[1] = 0.0; 
    bxby[1] = 0.0; 
    bxbz[1] = 0.0; 
    byby[1] = 0.0; 
    bybz[1] = 0.0; 
    bzbz[1] = 0.0; 
  // If bxbx, byby, or bzbz < 0.0 at the quadrature points, 
  // set cell_avg_magB2 to be true if it was not true before. 
  cell_avg_magB2[0] = 1; 
  } 
  // Calculate b_i = B_i/|B| by taking square root of B_i^2/|B|^2 at quadrature points. 
  // Uses the sign of B_i at quadrature points to get the correct sign of b_i. 
  // Also checks if B_i^2/|B|^2 < 0.0 at quadrature points and zeros out the value there. 
  ser_1x_p1_sqrt_with_sign(B_x, bxbx, bx); 
  ser_1x_p1_sqrt_with_sign(B_y, byby, by); 
  ser_1x_p1_sqrt_with_sign(B_z, bzbz, bz); 
} 
 
