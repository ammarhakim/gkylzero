#include <gkyl_mat.h> 
#include <gkyl_maxwell_kernels.h> 
#include <gkyl_basis_tensor_1x_2p_sqrt_with_sign.h> 
GKYL_CU_DH void em_copy_bvar_1x_tensor_p2(int count, struct gkyl_nmat *x, const double *em, int* cell_avg_bb, 
    double* GKYL_RESTRICT bvar, double* GKYL_RESTRICT bvar_surf) 
{ 
  // count:       Integer to indicate which matrix being fetched. 
  // x:           Input solution vector. 
  // em:          Input electromagnetic fields. 
  // cell_avg_bb: Output flag for cell average if bb only used cell averages. 
  // bvar:        Output volume expansion of b_i = B_i/|B| (first 3 components), b_i b_j = B_i B_j/|B|^2 (last 6 components). 
  // bvar_surf:   Output surface expansion of b_i = B_i/|B|. 
  //              [bx_xl, bx_xr, by_yl, by_yr, bz_zl, bz_zr]. 
 
  struct gkyl_mat x_bxbx = gkyl_nmat_get(x, count); 
  struct gkyl_mat x_bxby = gkyl_nmat_get(x, count+1); 
  struct gkyl_mat x_bxbz = gkyl_nmat_get(x, count+2); 
  struct gkyl_mat x_byby = gkyl_nmat_get(x, count+3); 
  struct gkyl_mat x_bybz = gkyl_nmat_get(x, count+4); 
  struct gkyl_mat x_bzbz = gkyl_nmat_get(x, count+5); 
  double *bx = &bvar[0]; 
  double *by = &bvar[3]; 
  double *bz = &bvar[6]; 
  double *bxbx = &bvar[9]; 
  double *bxby = &bvar[12]; 
  double *bxbz = &bvar[15]; 
  double *byby = &bvar[18]; 
  double *bybz = &bvar[21]; 
  double *bzbz = &bvar[24]; 
 
  bxbx[0] = gkyl_mat_get(&x_bxbx,0,0); 
  bxby[0] = gkyl_mat_get(&x_bxby,0,0); 
  bxbz[0] = gkyl_mat_get(&x_bxbz,0,0); 
  byby[0] = gkyl_mat_get(&x_byby,0,0); 
  bybz[0] = gkyl_mat_get(&x_bybz,0,0); 
  bzbz[0] = gkyl_mat_get(&x_bzbz,0,0); 
  bxbx[1] = gkyl_mat_get(&x_bxbx,1,0); 
  bxby[1] = gkyl_mat_get(&x_bxby,1,0); 
  bxbz[1] = gkyl_mat_get(&x_bxbz,1,0); 
  byby[1] = gkyl_mat_get(&x_byby,1,0); 
  bybz[1] = gkyl_mat_get(&x_bybz,1,0); 
  bzbz[1] = gkyl_mat_get(&x_bzbz,1,0); 
  bxbx[2] = gkyl_mat_get(&x_bxbx,2,0); 
  bxby[2] = gkyl_mat_get(&x_bxby,2,0); 
  bxbz[2] = gkyl_mat_get(&x_bxbz,2,0); 
  byby[2] = gkyl_mat_get(&x_byby,2,0); 
  bybz[2] = gkyl_mat_get(&x_bybz,2,0); 
  bzbz[2] = gkyl_mat_get(&x_bzbz,2,0); 
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
 
  double *bx_xl = &bvar_surf[0]; 
  double *bx_xr = &bvar_surf[1]; 
 
  bx_xl[0] = 1.58113883008419*bx[2]-1.224744871391589*bx[1]+0.7071067811865475*bx[0]; 
  bx_xr[0] = 1.58113883008419*bx[2]+1.224744871391589*bx[1]+0.7071067811865475*bx[0]; 
 
} 
 