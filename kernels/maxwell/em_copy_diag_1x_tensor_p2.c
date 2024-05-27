#include <gkyl_mat.h> 
#include <gkyl_maxwell_kernels.h> 
#include <gkyl_basis_tensor_1x_2p_sqrt_with_sign.h> 
GKYL_CU_DH void em_copy_diag_1x_tensor_p2(int count, struct gkyl_nmat *x, const double *em, int* cell_avg_bb, 
    double* GKYL_RESTRICT out) 
{ 
  // count:       Integer to indicate which matrix being fetched. 
  // x:           Input solution vector. 
  // em:          Input electromagnetic fields. 
  // cell_avg_bb: Output flag for cell average if bb only used cell averages. 
  // out:         Output volume expansion of diagnostic EM variables. 
 
  struct gkyl_mat x_bxbx = gkyl_nmat_get(x, count); 
  struct gkyl_mat x_bxby = gkyl_nmat_get(x, count+1); 
  struct gkyl_mat x_bxbz = gkyl_nmat_get(x, count+2); 
  struct gkyl_mat x_byby = gkyl_nmat_get(x, count+3); 
  struct gkyl_mat x_bybz = gkyl_nmat_get(x, count+4); 
  struct gkyl_mat x_bzbz = gkyl_nmat_get(x, count+5); 
  struct gkyl_mat x_ExBx = gkyl_nmat_get(x, count+6); 
  struct gkyl_mat x_ExBy = gkyl_nmat_get(x, count+7); 
  struct gkyl_mat x_ExBz = gkyl_nmat_get(x, count+8); 
  double *bx = &out[0]; 
  double *by = &out[3]; 
  double *bz = &out[6]; 
  double *bxbx = &out[9]; 
  double *bxby = &out[12]; 
  double *bxbz = &out[15]; 
  double *byby = &out[18]; 
  double *bybz = &out[21]; 
  double *bzbz = &out[24]; 
  double *ExBx = &out[27]; 
  double *ExBy = &out[30]; 
  double *ExBz = &out[33]; 
 
  bxbx[0] = gkyl_mat_get(&x_bxbx,0,0); 
  bxby[0] = gkyl_mat_get(&x_bxby,0,0); 
  bxbz[0] = gkyl_mat_get(&x_bxbz,0,0); 
  byby[0] = gkyl_mat_get(&x_byby,0,0); 
  bybz[0] = gkyl_mat_get(&x_bybz,0,0); 
  bzbz[0] = gkyl_mat_get(&x_bzbz,0,0); 
  ExBx[0] = gkyl_mat_get(&x_ExBx,0,0); 
  ExBy[0] = gkyl_mat_get(&x_ExBy,0,0); 
  ExBz[0] = gkyl_mat_get(&x_ExBz,0,0); 
  bxbx[1] = gkyl_mat_get(&x_bxbx,1,0); 
  bxby[1] = gkyl_mat_get(&x_bxby,1,0); 
  bxbz[1] = gkyl_mat_get(&x_bxbz,1,0); 
  byby[1] = gkyl_mat_get(&x_byby,1,0); 
  bybz[1] = gkyl_mat_get(&x_bybz,1,0); 
  bzbz[1] = gkyl_mat_get(&x_bzbz,1,0); 
  ExBx[1] = gkyl_mat_get(&x_ExBx,1,0); 
  ExBy[1] = gkyl_mat_get(&x_ExBy,1,0); 
  ExBz[1] = gkyl_mat_get(&x_ExBz,1,0); 
  bxbx[2] = gkyl_mat_get(&x_bxbx,2,0); 
  bxby[2] = gkyl_mat_get(&x_bxby,2,0); 
  bxbz[2] = gkyl_mat_get(&x_bxbz,2,0); 
  byby[2] = gkyl_mat_get(&x_byby,2,0); 
  bybz[2] = gkyl_mat_get(&x_bybz,2,0); 
  bzbz[2] = gkyl_mat_get(&x_bzbz,2,0); 
  ExBx[2] = gkyl_mat_get(&x_ExBx,2,0); 
  ExBy[2] = gkyl_mat_get(&x_ExBy,2,0); 
  ExBz[2] = gkyl_mat_get(&x_ExBz,2,0); 
  const double *B_x = &em[6]; 
  const double *B_y = &em[8]; 
  const double *B_z = &em[10]; 
 
  int cell_avg_bxbx = 0;
  int cell_avg_byby = 0;
  int cell_avg_bzbz = 0;
  if (0.6324555320336759*bxbx[2]-0.9486832980505137*bxbx[1]+0.7071067811865475*bxbx[0] < 0.0) cell_avg_bxbx = 1; 
  if (0.6324555320336759*byby[2]-0.9486832980505137*byby[1]+0.7071067811865475*byby[0] < 0.0) cell_avg_byby = 1; 
  if (0.6324555320336759*bzbz[2]-0.9486832980505137*bzbz[1]+0.7071067811865475*bzbz[0] < 0.0) cell_avg_bzbz = 1; 
  if (0.7071067811865475*bxbx[0]-0.7905694150420947*bxbx[2] < 0.0) cell_avg_bxbx = 1; 
  if (0.7071067811865475*byby[0]-0.7905694150420947*byby[2] < 0.0) cell_avg_byby = 1; 
  if (0.7071067811865475*bzbz[0]-0.7905694150420947*bzbz[2] < 0.0) cell_avg_bzbz = 1; 
  if (0.6324555320336759*bxbx[2]+0.9486832980505137*bxbx[1]+0.7071067811865475*bxbx[0] < 0.0) cell_avg_bxbx = 1; 
  if (0.6324555320336759*byby[2]+0.9486832980505137*byby[1]+0.7071067811865475*byby[0] < 0.0) cell_avg_byby = 1; 
  if (0.6324555320336759*bzbz[2]+0.9486832980505137*bzbz[1]+0.7071067811865475*bzbz[0] < 0.0) cell_avg_bzbz = 1; 
  if (cell_avg_bxbx || cell_avg_byby || cell_avg_bzbz) { 
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
    cell_avg_bb[0] = cell_avg_bxbx; 
    cell_avg_bb[1] = cell_avg_byby; 
    cell_avg_bb[2] = cell_avg_bzbz; 
  } 
  // Calculate b_i = B_i/|B| by taking square root of B_i^2/|B|^2 at quadrature points. 
  // Uses the sign of B_i at quadrature points to get the correct sign of b_i. 
  // Also checks if B_i^2/|B|^2 < 0.0 at quadrature points and zeros out the value there. 
  tensor_1x_2p_sqrt_with_sign(B_x, bxbx, bx); 
  tensor_1x_2p_sqrt_with_sign(B_y, byby, by); 
  tensor_1x_2p_sqrt_with_sign(B_z, bzbz, bz); 
 
} 
