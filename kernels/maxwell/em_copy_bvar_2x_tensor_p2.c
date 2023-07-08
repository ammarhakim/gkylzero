#include <gkyl_mat.h> 
#include <gkyl_maxwell_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
#include <gkyl_basis_tensor_2x_p2_sqrt_with_sign.h> 
GKYL_CU_DH void em_copy_bvar_2x_tensor_p2(struct gkyl_mat *x, const double *em, int* cell_avg_magB2, double* bvar) 
{ 
  // x:              Input solution vector for 1/|B|^2. 
  // em:             Input electromagnetic fields. 
  // cell_avg_magB2: Input flag for cell average if 1/|B|^2 only used cell averages. Also is adjusted if bb < 0 at control points. 
  // bvar:           Output b_i = B_i/|B| (first 3 components), b_i b_j = B_i B_j/|B|^2 (last 6 components). 
 
  double magB2_inv[9] = {0.0}; 

  magB2_inv[0] = gkyl_mat_get(x,0,0); 
  magB2_inv[1] = gkyl_mat_get(x,1,0); 
  magB2_inv[2] = gkyl_mat_get(x,2,0); 
  magB2_inv[3] = gkyl_mat_get(x,3,0); 
  magB2_inv[4] = gkyl_mat_get(x,4,0); 
  magB2_inv[5] = gkyl_mat_get(x,5,0); 
  magB2_inv[6] = gkyl_mat_get(x,6,0); 
  magB2_inv[7] = gkyl_mat_get(x,7,0); 
  magB2_inv[8] = gkyl_mat_get(x,8,0); 

  double *bx = &bvar[0]; 
  double *by = &bvar[9]; 
  double *bz = &bvar[18]; 
  double *bxbx = &bvar[27]; 
  double *bxby = &bvar[36]; 
  double *bxbz = &bvar[45]; 
  double *byby = &bvar[54]; 
  double *bybz = &bvar[63]; 
  double *bzbz = &bvar[72]; 
 
  const double *B_x = &em[27]; 
  const double *B_y = &em[36]; 
  const double *B_z = &em[45]; 
 
  // Calculate B_i B_j. 
  double B_x_sq[9] = {0.0}; 
  binop_mul_2d_tensor_p2(B_x, B_x, B_x_sq); 
 
  double B_y_sq[9] = {0.0}; 
  binop_mul_2d_tensor_p2(B_y, B_y, B_y_sq); 
 
  double B_z_sq[9] = {0.0}; 
  binop_mul_2d_tensor_p2(B_z, B_z, B_z_sq); 
 
  double B_x_B_y[9] = {0.0}; 
  binop_mul_2d_tensor_p2(B_x, B_y, B_x_B_y); 
 
  double B_x_B_z[9] = {0.0}; 
  binop_mul_2d_tensor_p2(B_x, B_z, B_x_B_z); 
 
  double B_y_B_z[9] = {0.0}; 
  binop_mul_2d_tensor_p2(B_y, B_z, B_y_B_z); 
 
  // Calculate expansions of B_i B_j/|B|^2, which can be calculated free of aliasing errors. 
  binop_mul_2d_tensor_p2(magB2_inv, B_x_sq, bxbx); 
  binop_mul_2d_tensor_p2(magB2_inv, B_y_sq, byby); 
  binop_mul_2d_tensor_p2(magB2_inv, B_z_sq, bzbz); 
  binop_mul_2d_tensor_p2(magB2_inv, B_x_B_y, bxby); 
  binop_mul_2d_tensor_p2(magB2_inv, B_x_B_z, bxbz); 
  binop_mul_2d_tensor_p2(magB2_inv, B_y_B_z, bybz); 
 
  int cell_avg = 0;
  if (0.4*bxbx[8]-0.5999999999999995*bxbx[7]-0.5999999999999999*bxbx[6]+0.4472135954999579*bxbx[5]+0.4472135954999579*bxbx[4]+0.9*bxbx[3]-0.6708203932499369*bxbx[2]-0.6708203932499369*bxbx[1]+0.5*bxbx[0] < 0.0) cell_avg = 1; 
  if (0.4*byby[8]-0.5999999999999995*byby[7]-0.5999999999999999*byby[6]+0.4472135954999579*byby[5]+0.4472135954999579*byby[4]+0.9*byby[3]-0.6708203932499369*byby[2]-0.6708203932499369*byby[1]+0.5*byby[0] < 0.0) cell_avg = 1; 
  if (0.4*bzbz[8]-0.5999999999999995*bzbz[7]-0.5999999999999999*bzbz[6]+0.4472135954999579*bzbz[5]+0.4472135954999579*bzbz[4]+0.9*bzbz[3]-0.6708203932499369*bzbz[2]-0.6708203932499369*bzbz[1]+0.5*bzbz[0] < 0.0) cell_avg = 1; 
  if ((-0.5*bxbx[8])+0.75*bxbx[7]-0.5590169943749475*bxbx[5]+0.4472135954999579*bxbx[4]-0.6708203932499369*bxbx[1]+0.5*bxbx[0] < 0.0) cell_avg = 1; 
  if ((-0.5*byby[8])+0.75*byby[7]-0.5590169943749475*byby[5]+0.4472135954999579*byby[4]-0.6708203932499369*byby[1]+0.5*byby[0] < 0.0) cell_avg = 1; 
  if ((-0.5*bzbz[8])+0.75*bzbz[7]-0.5590169943749475*bzbz[5]+0.4472135954999579*bzbz[4]-0.6708203932499369*bzbz[1]+0.5*bzbz[0] < 0.0) cell_avg = 1; 
  if (0.4*bxbx[8]-0.5999999999999995*bxbx[7]+0.5999999999999999*bxbx[6]+0.4472135954999579*bxbx[5]+0.4472135954999579*bxbx[4]-0.9*bxbx[3]+0.6708203932499369*bxbx[2]-0.6708203932499369*bxbx[1]+0.5*bxbx[0] < 0.0) cell_avg = 1; 
  if (0.4*byby[8]-0.5999999999999995*byby[7]+0.5999999999999999*byby[6]+0.4472135954999579*byby[5]+0.4472135954999579*byby[4]-0.9*byby[3]+0.6708203932499369*byby[2]-0.6708203932499369*byby[1]+0.5*byby[0] < 0.0) cell_avg = 1; 
  if (0.4*bzbz[8]-0.5999999999999995*bzbz[7]+0.5999999999999999*bzbz[6]+0.4472135954999579*bzbz[5]+0.4472135954999579*bzbz[4]-0.9*bzbz[3]+0.6708203932499369*bzbz[2]-0.6708203932499369*bzbz[1]+0.5*bzbz[0] < 0.0) cell_avg = 1; 
  if ((-0.5*bxbx[8])+0.75*bxbx[6]+0.4472135954999579*bxbx[5]-0.5590169943749475*bxbx[4]-0.6708203932499369*bxbx[2]+0.5*bxbx[0] < 0.0) cell_avg = 1; 
  if ((-0.5*byby[8])+0.75*byby[6]+0.4472135954999579*byby[5]-0.5590169943749475*byby[4]-0.6708203932499369*byby[2]+0.5*byby[0] < 0.0) cell_avg = 1; 
  if ((-0.5*bzbz[8])+0.75*bzbz[6]+0.4472135954999579*bzbz[5]-0.5590169943749475*bzbz[4]-0.6708203932499369*bzbz[2]+0.5*bzbz[0] < 0.0) cell_avg = 1; 
  if (0.625*bxbx[8]-0.5590169943749475*bxbx[5]-0.5590169943749475*bxbx[4]+0.5*bxbx[0] < 0.0) cell_avg = 1; 
  if (0.625*byby[8]-0.5590169943749475*byby[5]-0.5590169943749475*byby[4]+0.5*byby[0] < 0.0) cell_avg = 1; 
  if (0.625*bzbz[8]-0.5590169943749475*bzbz[5]-0.5590169943749475*bzbz[4]+0.5*bzbz[0] < 0.0) cell_avg = 1; 
  if ((-0.5*bxbx[8])-0.75*bxbx[6]+0.4472135954999579*bxbx[5]-0.5590169943749475*bxbx[4]+0.6708203932499369*bxbx[2]+0.5*bxbx[0] < 0.0) cell_avg = 1; 
  if ((-0.5*byby[8])-0.75*byby[6]+0.4472135954999579*byby[5]-0.5590169943749475*byby[4]+0.6708203932499369*byby[2]+0.5*byby[0] < 0.0) cell_avg = 1; 
  if ((-0.5*bzbz[8])-0.75*bzbz[6]+0.4472135954999579*bzbz[5]-0.5590169943749475*bzbz[4]+0.6708203932499369*bzbz[2]+0.5*bzbz[0] < 0.0) cell_avg = 1; 
  if (0.4*bxbx[8]+0.5999999999999995*bxbx[7]-0.5999999999999999*bxbx[6]+0.4472135954999579*bxbx[5]+0.4472135954999579*bxbx[4]-0.9*bxbx[3]-0.6708203932499369*bxbx[2]+0.6708203932499369*bxbx[1]+0.5*bxbx[0] < 0.0) cell_avg = 1; 
  if (0.4*byby[8]+0.5999999999999995*byby[7]-0.5999999999999999*byby[6]+0.4472135954999579*byby[5]+0.4472135954999579*byby[4]-0.9*byby[3]-0.6708203932499369*byby[2]+0.6708203932499369*byby[1]+0.5*byby[0] < 0.0) cell_avg = 1; 
  if (0.4*bzbz[8]+0.5999999999999995*bzbz[7]-0.5999999999999999*bzbz[6]+0.4472135954999579*bzbz[5]+0.4472135954999579*bzbz[4]-0.9*bzbz[3]-0.6708203932499369*bzbz[2]+0.6708203932499369*bzbz[1]+0.5*bzbz[0] < 0.0) cell_avg = 1; 
  if ((-0.5*bxbx[8])-0.75*bxbx[7]-0.5590169943749475*bxbx[5]+0.4472135954999579*bxbx[4]+0.6708203932499369*bxbx[1]+0.5*bxbx[0] < 0.0) cell_avg = 1; 
  if ((-0.5*byby[8])-0.75*byby[7]-0.5590169943749475*byby[5]+0.4472135954999579*byby[4]+0.6708203932499369*byby[1]+0.5*byby[0] < 0.0) cell_avg = 1; 
  if ((-0.5*bzbz[8])-0.75*bzbz[7]-0.5590169943749475*bzbz[5]+0.4472135954999579*bzbz[4]+0.6708203932499369*bzbz[1]+0.5*bzbz[0] < 0.0) cell_avg = 1; 
  if (0.4*bxbx[8]+0.5999999999999995*bxbx[7]+0.5999999999999999*bxbx[6]+0.4472135954999579*bxbx[5]+0.4472135954999579*bxbx[4]+0.9*bxbx[3]+0.6708203932499369*bxbx[2]+0.6708203932499369*bxbx[1]+0.5*bxbx[0] < 0.0) cell_avg = 1; 
  if (0.4*byby[8]+0.5999999999999995*byby[7]+0.5999999999999999*byby[6]+0.4472135954999579*byby[5]+0.4472135954999579*byby[4]+0.9*byby[3]+0.6708203932499369*byby[2]+0.6708203932499369*byby[1]+0.5*byby[0] < 0.0) cell_avg = 1; 
  if (0.4*bzbz[8]+0.5999999999999995*bzbz[7]+0.5999999999999999*bzbz[6]+0.4472135954999579*bzbz[5]+0.4472135954999579*bzbz[4]+0.9*bzbz[3]+0.6708203932499369*bzbz[2]+0.6708203932499369*bzbz[1]+0.5*bzbz[0] < 0.0) cell_avg = 1; 
  if (cell_avg || cell_avg_magB2[0]) { 
    bxbx[1] = 0.0; 
    bxby[1] = 0.0; 
    bxbz[1] = 0.0; 
    byby[1] = 0.0; 
    bybz[1] = 0.0; 
    bzbz[1] = 0.0; 
    bxbx[2] = 0.0; 
    bxby[2] = 0.0; 
    bxbz[2] = 0.0; 
    byby[2] = 0.0; 
    bybz[2] = 0.0; 
    bzbz[2] = 0.0; 
    bxbx[3] = 0.0; 
    bxby[3] = 0.0; 
    bxbz[3] = 0.0; 
    byby[3] = 0.0; 
    bybz[3] = 0.0; 
    bzbz[3] = 0.0; 
    bxbx[4] = 0.0; 
    bxby[4] = 0.0; 
    bxbz[4] = 0.0; 
    byby[4] = 0.0; 
    bybz[4] = 0.0; 
    bzbz[4] = 0.0; 
    bxbx[5] = 0.0; 
    bxby[5] = 0.0; 
    bxbz[5] = 0.0; 
    byby[5] = 0.0; 
    bybz[5] = 0.0; 
    bzbz[5] = 0.0; 
    bxbx[6] = 0.0; 
    bxby[6] = 0.0; 
    bxbz[6] = 0.0; 
    byby[6] = 0.0; 
    bybz[6] = 0.0; 
    bzbz[6] = 0.0; 
    bxbx[7] = 0.0; 
    bxby[7] = 0.0; 
    bxbz[7] = 0.0; 
    byby[7] = 0.0; 
    bybz[7] = 0.0; 
    bzbz[7] = 0.0; 
    bxbx[8] = 0.0; 
    bxby[8] = 0.0; 
    bxbz[8] = 0.0; 
    byby[8] = 0.0; 
    bybz[8] = 0.0; 
    bzbz[8] = 0.0; 
  // If bxbx, byby, or bzbz < 0.0 at the quadrature points, 
  // set cell_avg_magB2 to be true if it was not true before. 
  cell_avg_magB2[0] = 1; 
  } 
  // Calculate b_i = B_i/|B| by taking square root of B_i^2/|B|^2 at quadrature points. 
  // Uses the sign of B_i at quadrature points to get the correct sign of b_i. 
  // Also checks if B_i^2/|B|^2 < 0.0 at quadrature points and zeros out the value there. 
  tensor_2x_p2_sqrt_with_sign(B_x, bxbx, bx); 
  tensor_2x_p2_sqrt_with_sign(B_y, byby, by); 
  tensor_2x_p2_sqrt_with_sign(B_z, bzbz, bz); 
} 
 
