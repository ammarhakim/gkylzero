#include <gkyl_mat.h> 
#include <gkyl_maxwell_kernels.h> 
#include <gkyl_basis_tensor_2x_2p_sqrt_with_sign.h> 
GKYL_CU_DH void em_copy_diag_2x_tensor_p2(int count, struct gkyl_nmat *x, const double *em, int* cell_avg_bb, 
    double* GKYL_RESTRICT em_vars_diag, double* GKYL_RESTRICT bvar, double* GKYL_RESTRICT bvar_surf) 
{ 
  // count:        Integer to indicate which matrix being fetched. 
  // x:            Input solution vector. 
  // em:           Input electromagnetic fields. 
  // cell_avg_bb:  Output flag for cell average if bb only used cell averages. 
  // em_vars_diag: Output volume expansion of diagnostic EM variables [bb, ExB/|B|^2]. 
  // bvar:         Output volume expansion of b_i = B_i/|B| (3 components). 
  // bvar_surf:    Output surface expansion of b_i = B_i/|B|. 
  //               [bx_xl, bx_xr, by_yl, by_yr, bz_zl, bz_zr]. 
 
  struct gkyl_mat x_bxbx = gkyl_nmat_get(x, count); 
  struct gkyl_mat x_bxby = gkyl_nmat_get(x, count+1); 
  struct gkyl_mat x_bxbz = gkyl_nmat_get(x, count+2); 
  struct gkyl_mat x_byby = gkyl_nmat_get(x, count+3); 
  struct gkyl_mat x_bybz = gkyl_nmat_get(x, count+4); 
  struct gkyl_mat x_bzbz = gkyl_nmat_get(x, count+5); 
  struct gkyl_mat x_ExBx = gkyl_nmat_get(x, count+6); 
  struct gkyl_mat x_ExBy = gkyl_nmat_get(x, count+7); 
  struct gkyl_mat x_ExBz = gkyl_nmat_get(x, count+8); 
  double *bx = &bvar[0]; 
  double *by = &bvar[4]; 
  double *bz = &bvar[8]; 
 
  double *bxbx = &em_vars_diag[0]; 
  double *bxby = &em_vars_diag[9]; 
  double *bxbz = &em_vars_diag[18]; 
  double *byby = &em_vars_diag[27]; 
  double *bybz = &em_vars_diag[36]; 
  double *bzbz = &em_vars_diag[45]; 
  double *ExBx = &em_vars_diag[54]; 
  double *ExBy = &em_vars_diag[63]; 
  double *ExBz = &em_vars_diag[72]; 
 
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
  bxbx[3] = gkyl_mat_get(&x_bxbx,3,0); 
  bxby[3] = gkyl_mat_get(&x_bxby,3,0); 
  bxbz[3] = gkyl_mat_get(&x_bxbz,3,0); 
  byby[3] = gkyl_mat_get(&x_byby,3,0); 
  bybz[3] = gkyl_mat_get(&x_bybz,3,0); 
  bzbz[3] = gkyl_mat_get(&x_bzbz,3,0); 
  ExBx[3] = gkyl_mat_get(&x_ExBx,3,0); 
  ExBy[3] = gkyl_mat_get(&x_ExBy,3,0); 
  ExBz[3] = gkyl_mat_get(&x_ExBz,3,0); 
  bxbx[4] = gkyl_mat_get(&x_bxbx,4,0); 
  bxby[4] = gkyl_mat_get(&x_bxby,4,0); 
  bxbz[4] = gkyl_mat_get(&x_bxbz,4,0); 
  byby[4] = gkyl_mat_get(&x_byby,4,0); 
  bybz[4] = gkyl_mat_get(&x_bybz,4,0); 
  bzbz[4] = gkyl_mat_get(&x_bzbz,4,0); 
  ExBx[4] = gkyl_mat_get(&x_ExBx,4,0); 
  ExBy[4] = gkyl_mat_get(&x_ExBy,4,0); 
  ExBz[4] = gkyl_mat_get(&x_ExBz,4,0); 
  bxbx[5] = gkyl_mat_get(&x_bxbx,5,0); 
  bxby[5] = gkyl_mat_get(&x_bxby,5,0); 
  bxbz[5] = gkyl_mat_get(&x_bxbz,5,0); 
  byby[5] = gkyl_mat_get(&x_byby,5,0); 
  bybz[5] = gkyl_mat_get(&x_bybz,5,0); 
  bzbz[5] = gkyl_mat_get(&x_bzbz,5,0); 
  ExBx[5] = gkyl_mat_get(&x_ExBx,5,0); 
  ExBy[5] = gkyl_mat_get(&x_ExBy,5,0); 
  ExBz[5] = gkyl_mat_get(&x_ExBz,5,0); 
  bxbx[6] = gkyl_mat_get(&x_bxbx,6,0); 
  bxby[6] = gkyl_mat_get(&x_bxby,6,0); 
  bxbz[6] = gkyl_mat_get(&x_bxbz,6,0); 
  byby[6] = gkyl_mat_get(&x_byby,6,0); 
  bybz[6] = gkyl_mat_get(&x_bybz,6,0); 
  bzbz[6] = gkyl_mat_get(&x_bzbz,6,0); 
  ExBx[6] = gkyl_mat_get(&x_ExBx,6,0); 
  ExBy[6] = gkyl_mat_get(&x_ExBy,6,0); 
  ExBz[6] = gkyl_mat_get(&x_ExBz,6,0); 
  bxbx[7] = gkyl_mat_get(&x_bxbx,7,0); 
  bxby[7] = gkyl_mat_get(&x_bxby,7,0); 
  bxbz[7] = gkyl_mat_get(&x_bxbz,7,0); 
  byby[7] = gkyl_mat_get(&x_byby,7,0); 
  bybz[7] = gkyl_mat_get(&x_bybz,7,0); 
  bzbz[7] = gkyl_mat_get(&x_bzbz,7,0); 
  ExBx[7] = gkyl_mat_get(&x_ExBx,7,0); 
  ExBy[7] = gkyl_mat_get(&x_ExBy,7,0); 
  ExBz[7] = gkyl_mat_get(&x_ExBz,7,0); 
  bxbx[8] = gkyl_mat_get(&x_bxbx,8,0); 
  bxby[8] = gkyl_mat_get(&x_bxby,8,0); 
  bxbz[8] = gkyl_mat_get(&x_bxbz,8,0); 
  byby[8] = gkyl_mat_get(&x_byby,8,0); 
  bybz[8] = gkyl_mat_get(&x_bybz,8,0); 
  bzbz[8] = gkyl_mat_get(&x_bzbz,8,0); 
  ExBx[8] = gkyl_mat_get(&x_ExBx,8,0); 
  ExBy[8] = gkyl_mat_get(&x_ExBy,8,0); 
  ExBz[8] = gkyl_mat_get(&x_ExBz,8,0); 
  const double *B_x = &em[12]; 
  const double *B_y = &em[16]; 
  const double *B_z = &em[20]; 
 
  int cell_avg_bxbx = 0;
  int cell_avg_byby = 0;
  int cell_avg_bzbz = 0;
  if (0.5*bxbx[3]-0.5*bxbx[2]-0.5*bxbx[1]+0.5*bxbx[0] < 0.0) cell_avg_bxbx = 1; 
  if (0.5*byby[3]-0.5*byby[2]-0.5*byby[1]+0.5*byby[0] < 0.0) cell_avg_byby = 1; 
  if (0.5*bzbz[3]-0.5*bzbz[2]-0.5*bzbz[1]+0.5*bzbz[0] < 0.0) cell_avg_bzbz = 1; 
  if ((-0.5*bxbx[3])+0.5*bxbx[2]-0.5*bxbx[1]+0.5*bxbx[0] < 0.0) cell_avg_bxbx = 1; 
  if ((-0.5*byby[3])+0.5*byby[2]-0.5*byby[1]+0.5*byby[0] < 0.0) cell_avg_byby = 1; 
  if ((-0.5*bzbz[3])+0.5*bzbz[2]-0.5*bzbz[1]+0.5*bzbz[0] < 0.0) cell_avg_bzbz = 1; 
  if ((-0.5*bxbx[3])-0.5*bxbx[2]+0.5*bxbx[1]+0.5*bxbx[0] < 0.0) cell_avg_bxbx = 1; 
  if ((-0.5*byby[3])-0.5*byby[2]+0.5*byby[1]+0.5*byby[0] < 0.0) cell_avg_byby = 1; 
  if ((-0.5*bzbz[3])-0.5*bzbz[2]+0.5*bzbz[1]+0.5*bzbz[0] < 0.0) cell_avg_bzbz = 1; 
  if (0.5*bxbx[3]+0.5*bxbx[2]+0.5*bxbx[1]+0.5*bxbx[0] < 0.0) cell_avg_bxbx = 1; 
  if (0.5*byby[3]+0.5*byby[2]+0.5*byby[1]+0.5*byby[0] < 0.0) cell_avg_byby = 1; 
  if (0.5*bzbz[3]+0.5*bzbz[2]+0.5*bzbz[1]+0.5*bzbz[0] < 0.0) cell_avg_bzbz = 1; 
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
    cell_avg_bb[0] = cell_avg_bxbx; 
    cell_avg_bb[1] = cell_avg_byby; 
    cell_avg_bb[2] = cell_avg_bzbz; 
  } 
  // Calculate b_i = B_i/|B| by taking square root of B_i^2/|B|^2 at quadrature points. 
  // Uses the sign of B_i at quadrature points to get the correct sign of b_i. 
  tensor_2x_2p_sqrt_with_sign(B_x, bxbx, bx); 
  tensor_2x_2p_sqrt_with_sign(B_y, byby, by); 
  tensor_2x_2p_sqrt_with_sign(B_z, bzbz, bz); 
 
  double *bx_xl = &bvar_surf[0]; 
  double *bx_xr = &bvar_surf[2]; 
 
  bx_xl[0] = 0.7071067811865475*bx[0]-1.224744871391589*bx[1]; 
  bx_xl[1] = 0.7071067811865475*bx[2]-1.224744871391589*bx[3]; 
 
  bx_xr[0] = 1.224744871391589*bx[1]+0.7071067811865475*bx[0]; 
  bx_xr[1] = 1.224744871391589*bx[3]+0.7071067811865475*bx[2]; 
 
  double *by_yl = &bvar_surf[4]; 
  double *by_yr = &bvar_surf[6]; 
 
  by_yl[0] = 0.7071067811865475*by[0]-1.224744871391589*by[2]; 
  by_yl[1] = 0.7071067811865475*by[1]-1.224744871391589*by[3]; 
 
  by_yr[0] = 1.224744871391589*by[2]+0.7071067811865475*by[0]; 
  by_yr[1] = 1.224744871391589*by[3]+0.7071067811865475*by[1]; 
 
} 
 
