#include <gkyl_mat.h> 
#include <gkyl_maxwell_kernels.h> 
#include <gkyl_basis_ser_2x_p1_sqrt_with_sign.h> 
GKYL_CU_DH void em_surf_copy_bvar_3x_ser_p1(int count, struct gkyl_nmat *x, const double *em, int* cell_avg_magB2_surf, double* GKYL_RESTRICT bvar_surf) 
{ 
  // count:               Integer to indicate which matrix being fetched. 
  // x:                   Input solution vector. 
  // em:                  Input electromagnetic fields. 
  // cell_avg_magB2_surf: Output flag for cell average if 1/|B|^2 at a surface only used cell averages. 
  // bvar_surf:           Output magnetic field unit tensor and unit vector at surfaces. 
  //                      [bx_xl, bx_xr, bxbx_xl, bxbx_xr, bxby_xl, bxby_xr, bxbz_xl, bxbz_xr, 
  //                       by_yl, by_yr, byby_yl, byby_yr, bxby_yl, bxby_yr, bybz_yl, bybz_yr, 
  //                       bz_zl, bz_zr, bzbz_zl, bzbz_zr, bxbz_zl, bxbz_zr, bybz_zl, bybz_zr] 
 
  struct gkyl_mat x_bxbx_xl = gkyl_nmat_get(x, count); 
  struct gkyl_mat x_bxbx_xr = gkyl_nmat_get(x, count+1); 
  struct gkyl_mat x_bxby_xl = gkyl_nmat_get(x, count+2); 
  struct gkyl_mat x_bxby_xr = gkyl_nmat_get(x, count+3); 
  struct gkyl_mat x_bxbz_xl = gkyl_nmat_get(x, count+4); 
  struct gkyl_mat x_bxbz_xr = gkyl_nmat_get(x, count+5); 
 
  double *bx_xl = &bvar_surf[0]; 
  double *bx_xr = &bvar_surf[4]; 
  double *bxbx_xl = &bvar_surf[8]; 
  double *bxbx_xr = &bvar_surf[12]; 
  double *bxby_xl = &bvar_surf[16]; 
  double *bxby_xr = &bvar_surf[20]; 
  double *bxbz_xl = &bvar_surf[24]; 
  double *bxbz_xr = &bvar_surf[28]; 
  int *cell_avg_magB2_xl = &cell_avg_magB2_surf[0]; 
  int *cell_avg_magB2_xr = &cell_avg_magB2_surf[1]; 
 
  struct gkyl_mat x_byby_yl = gkyl_nmat_get(x, count+6); 
  struct gkyl_mat x_byby_yr = gkyl_nmat_get(x, count+7); 
  struct gkyl_mat x_bxby_yl = gkyl_nmat_get(x, count+8); 
  struct gkyl_mat x_bxby_yr = gkyl_nmat_get(x, count+9); 
  struct gkyl_mat x_bybz_yl = gkyl_nmat_get(x, count+10); 
  struct gkyl_mat x_bybz_yr = gkyl_nmat_get(x, count+11); 
 
  double *by_yl = &bvar_surf[32]; 
  double *by_yr = &bvar_surf[36]; 
  double *byby_yl = &bvar_surf[40]; 
  double *byby_yr = &bvar_surf[44]; 
  double *bxby_yl = &bvar_surf[48]; 
  double *bxby_yr = &bvar_surf[52]; 
  double *bybz_yl = &bvar_surf[56]; 
  double *bybz_yr = &bvar_surf[60]; 
  int *cell_avg_magB2_yl = &cell_avg_magB2_surf[2]; 
  int *cell_avg_magB2_yr = &cell_avg_magB2_surf[3]; 
 
  struct gkyl_mat x_bzbz_zl = gkyl_nmat_get(x, count+12); 
  struct gkyl_mat x_bzbz_zr = gkyl_nmat_get(x, count+13); 
  struct gkyl_mat x_bxbz_zl = gkyl_nmat_get(x, count+14); 
  struct gkyl_mat x_bxbz_zr = gkyl_nmat_get(x, count+15); 
  struct gkyl_mat x_bybz_zl = gkyl_nmat_get(x, count+16); 
  struct gkyl_mat x_bybz_zr = gkyl_nmat_get(x, count+17); 
 
  double *bz_zl = &bvar_surf[64]; 
  double *bz_zr = &bvar_surf[68]; 
  double *bzbz_zl = &bvar_surf[72]; 
  double *bzbz_zr = &bvar_surf[76]; 
  double *bxbz_zl = &bvar_surf[80]; 
  double *bxbz_zr = &bvar_surf[84]; 
  double *bybz_zl = &bvar_surf[88]; 
  double *bybz_zr = &bvar_surf[92]; 
  int *cell_avg_magB2_zl = &cell_avg_magB2_surf[4]; 
  int *cell_avg_magB2_zr = &cell_avg_magB2_surf[5]; 
 
  bxbx_xl[0] = gkyl_mat_get(&x_bxbx_xl,0,0); 
  bxbx_xr[0] = gkyl_mat_get(&x_bxbx_xr,0,0); 
  bxby_xl[0] = gkyl_mat_get(&x_bxby_xl,0,0); 
  bxby_xr[0] = gkyl_mat_get(&x_bxby_xr,0,0); 
  bxbz_xl[0] = gkyl_mat_get(&x_bxbz_xl,0,0); 
  bxbz_xr[0] = gkyl_mat_get(&x_bxbz_xr,0,0); 
 
  byby_yl[0] = gkyl_mat_get(&x_byby_yl,0,0); 
  byby_yr[0] = gkyl_mat_get(&x_byby_yr,0,0); 
  bxby_yl[0] = gkyl_mat_get(&x_bxby_yl,0,0); 
  bxby_yr[0] = gkyl_mat_get(&x_bxby_yr,0,0); 
  bybz_yl[0] = gkyl_mat_get(&x_bybz_yl,0,0); 
  bybz_yr[0] = gkyl_mat_get(&x_bybz_yr,0,0); 
 
  bzbz_zl[0] = gkyl_mat_get(&x_bzbz_zl,0,0); 
  bzbz_zr[0] = gkyl_mat_get(&x_bzbz_zr,0,0); 
  bxbz_zl[0] = gkyl_mat_get(&x_bxbz_zl,0,0); 
  bxbz_zr[0] = gkyl_mat_get(&x_bxbz_zr,0,0); 
  bybz_zl[0] = gkyl_mat_get(&x_bybz_zl,0,0); 
  bybz_zr[0] = gkyl_mat_get(&x_bybz_zr,0,0); 
 
  bxbx_xl[1] = gkyl_mat_get(&x_bxbx_xl,1,0); 
  bxbx_xr[1] = gkyl_mat_get(&x_bxbx_xr,1,0); 
  bxby_xl[1] = gkyl_mat_get(&x_bxby_xl,1,0); 
  bxby_xr[1] = gkyl_mat_get(&x_bxby_xr,1,0); 
  bxbz_xl[1] = gkyl_mat_get(&x_bxbz_xl,1,0); 
  bxbz_xr[1] = gkyl_mat_get(&x_bxbz_xr,1,0); 
 
  byby_yl[1] = gkyl_mat_get(&x_byby_yl,1,0); 
  byby_yr[1] = gkyl_mat_get(&x_byby_yr,1,0); 
  bxby_yl[1] = gkyl_mat_get(&x_bxby_yl,1,0); 
  bxby_yr[1] = gkyl_mat_get(&x_bxby_yr,1,0); 
  bybz_yl[1] = gkyl_mat_get(&x_bybz_yl,1,0); 
  bybz_yr[1] = gkyl_mat_get(&x_bybz_yr,1,0); 
 
  bzbz_zl[1] = gkyl_mat_get(&x_bzbz_zl,1,0); 
  bzbz_zr[1] = gkyl_mat_get(&x_bzbz_zr,1,0); 
  bxbz_zl[1] = gkyl_mat_get(&x_bxbz_zl,1,0); 
  bxbz_zr[1] = gkyl_mat_get(&x_bxbz_zr,1,0); 
  bybz_zl[1] = gkyl_mat_get(&x_bybz_zl,1,0); 
  bybz_zr[1] = gkyl_mat_get(&x_bybz_zr,1,0); 
 
  bxbx_xl[2] = gkyl_mat_get(&x_bxbx_xl,2,0); 
  bxbx_xr[2] = gkyl_mat_get(&x_bxbx_xr,2,0); 
  bxby_xl[2] = gkyl_mat_get(&x_bxby_xl,2,0); 
  bxby_xr[2] = gkyl_mat_get(&x_bxby_xr,2,0); 
  bxbz_xl[2] = gkyl_mat_get(&x_bxbz_xl,2,0); 
  bxbz_xr[2] = gkyl_mat_get(&x_bxbz_xr,2,0); 
 
  byby_yl[2] = gkyl_mat_get(&x_byby_yl,2,0); 
  byby_yr[2] = gkyl_mat_get(&x_byby_yr,2,0); 
  bxby_yl[2] = gkyl_mat_get(&x_bxby_yl,2,0); 
  bxby_yr[2] = gkyl_mat_get(&x_bxby_yr,2,0); 
  bybz_yl[2] = gkyl_mat_get(&x_bybz_yl,2,0); 
  bybz_yr[2] = gkyl_mat_get(&x_bybz_yr,2,0); 
 
  bzbz_zl[2] = gkyl_mat_get(&x_bzbz_zl,2,0); 
  bzbz_zr[2] = gkyl_mat_get(&x_bzbz_zr,2,0); 
  bxbz_zl[2] = gkyl_mat_get(&x_bxbz_zl,2,0); 
  bxbz_zr[2] = gkyl_mat_get(&x_bxbz_zr,2,0); 
  bybz_zl[2] = gkyl_mat_get(&x_bybz_zl,2,0); 
  bybz_zr[2] = gkyl_mat_get(&x_bybz_zr,2,0); 
 
  bxbx_xl[3] = gkyl_mat_get(&x_bxbx_xl,3,0); 
  bxbx_xr[3] = gkyl_mat_get(&x_bxbx_xr,3,0); 
  bxby_xl[3] = gkyl_mat_get(&x_bxby_xl,3,0); 
  bxby_xr[3] = gkyl_mat_get(&x_bxby_xr,3,0); 
  bxbz_xl[3] = gkyl_mat_get(&x_bxbz_xl,3,0); 
  bxbz_xr[3] = gkyl_mat_get(&x_bxbz_xr,3,0); 
 
  byby_yl[3] = gkyl_mat_get(&x_byby_yl,3,0); 
  byby_yr[3] = gkyl_mat_get(&x_byby_yr,3,0); 
  bxby_yl[3] = gkyl_mat_get(&x_bxby_yl,3,0); 
  bxby_yr[3] = gkyl_mat_get(&x_bxby_yr,3,0); 
  bybz_yl[3] = gkyl_mat_get(&x_bybz_yl,3,0); 
  bybz_yr[3] = gkyl_mat_get(&x_bybz_yr,3,0); 
 
  bzbz_zl[3] = gkyl_mat_get(&x_bzbz_zl,3,0); 
  bzbz_zr[3] = gkyl_mat_get(&x_bzbz_zr,3,0); 
  bxbz_zl[3] = gkyl_mat_get(&x_bxbz_zl,3,0); 
  bxbz_zr[3] = gkyl_mat_get(&x_bxbz_zr,3,0); 
  bybz_zl[3] = gkyl_mat_get(&x_bybz_zl,3,0); 
  bybz_zr[3] = gkyl_mat_get(&x_bybz_zr,3,0); 
 
  const double *B_x = &em[24]; 
  const double *B_y = &em[32]; 
  const double *B_z = &em[40]; 
 
  int cell_avg_xl = 0;
  int cell_avg_xr = 0;
  if (0.5*bxbx_xl[3]-0.5*bxbx_xl[2]-0.5*bxbx_xl[1]+0.5*bxbx_xl[0] < 0.0) cell_avg_xl = 1; 
  if (0.5*bxbx_xr[3]-0.5*bxbx_xr[2]-0.5*bxbx_xr[1]+0.5*bxbx_xr[0] < 0.0) cell_avg_xr = 1; 
  if ((-0.5*bxbx_xl[3])+0.5*bxbx_xl[2]-0.5*bxbx_xl[1]+0.5*bxbx_xl[0] < 0.0) cell_avg_xl = 1; 
  if ((-0.5*bxbx_xr[3])+0.5*bxbx_xr[2]-0.5*bxbx_xr[1]+0.5*bxbx_xr[0] < 0.0) cell_avg_xr = 1; 
  if ((-0.5*bxbx_xl[3])-0.5*bxbx_xl[2]+0.5*bxbx_xl[1]+0.5*bxbx_xl[0] < 0.0) cell_avg_xl = 1; 
  if ((-0.5*bxbx_xr[3])-0.5*bxbx_xr[2]+0.5*bxbx_xr[1]+0.5*bxbx_xr[0] < 0.0) cell_avg_xr = 1; 
  if (0.5*bxbx_xl[3]+0.5*bxbx_xl[2]+0.5*bxbx_xl[1]+0.5*bxbx_xl[0] < 0.0) cell_avg_xl = 1; 
  if (0.5*bxbx_xr[3]+0.5*bxbx_xr[2]+0.5*bxbx_xr[1]+0.5*bxbx_xr[0] < 0.0) cell_avg_xr = 1; 
 
  if (cell_avg_xl || cell_avg_magB2_xl[0]) { 
    bxbx_xl[1] = 0.0; 
    bxby_xl[1] = 0.0; 
    bxbz_xl[1] = 0.0; 
    bxbx_xl[2] = 0.0; 
    bxby_xl[2] = 0.0; 
    bxbz_xl[2] = 0.0; 
    bxbx_xl[3] = 0.0; 
    bxby_xl[3] = 0.0; 
    bxbz_xl[3] = 0.0; 
    // If bxbx, bxby, or bxbz < 0.0 at the lower x surface quadrature points, 
    // set cell_avg_magB2_xl to be true in case it was not true before. 
    cell_avg_magB2_xl[0] = 1; 
  } 
 
  if (cell_avg_xr || cell_avg_magB2_xr[0]) { 
    bxbx_xr[1] = 0.0; 
    bxby_xr[1] = 0.0; 
    bxbz_xr[1] = 0.0; 
    bxbx_xr[2] = 0.0; 
    bxby_xr[2] = 0.0; 
    bxbz_xr[2] = 0.0; 
    bxbx_xr[3] = 0.0; 
    bxby_xr[3] = 0.0; 
    bxbz_xr[3] = 0.0; 
    // If bxbx, bxby, or bxbz < 0.0 at the upper x surface quadrature points, 
    // set cell_avg_magB2_xr to be true in case it was not true before. 
    cell_avg_magB2_xr[0] = 1; 
  } 
 
  int cell_avg_yl = 0;
  int cell_avg_yr = 0;
  if (0.5*byby_yl[3]-0.5*byby_yl[2]-0.5*byby_yl[1]+0.5*byby_yl[0] < 0.0) cell_avg_yl = 1; 
  if (0.5*byby_yr[3]-0.5*byby_yr[2]-0.5*byby_yr[1]+0.5*byby_yr[0] < 0.0) cell_avg_yr = 1; 
  if ((-0.5*byby_yl[3])+0.5*byby_yl[2]-0.5*byby_yl[1]+0.5*byby_yl[0] < 0.0) cell_avg_yl = 1; 
  if ((-0.5*byby_yr[3])+0.5*byby_yr[2]-0.5*byby_yr[1]+0.5*byby_yr[0] < 0.0) cell_avg_yr = 1; 
  if ((-0.5*byby_yl[3])-0.5*byby_yl[2]+0.5*byby_yl[1]+0.5*byby_yl[0] < 0.0) cell_avg_yl = 1; 
  if ((-0.5*byby_yr[3])-0.5*byby_yr[2]+0.5*byby_yr[1]+0.5*byby_yr[0] < 0.0) cell_avg_yr = 1; 
  if (0.5*byby_yl[3]+0.5*byby_yl[2]+0.5*byby_yl[1]+0.5*byby_yl[0] < 0.0) cell_avg_yl = 1; 
  if (0.5*byby_yr[3]+0.5*byby_yr[2]+0.5*byby_yr[1]+0.5*byby_yr[0] < 0.0) cell_avg_yr = 1; 
 
  if (cell_avg_yl || cell_avg_magB2_yl[0]) { 
    byby_yl[1] = 0.0; 
    bxby_yl[1] = 0.0; 
    bybz_yl[1] = 0.0; 
    byby_yl[2] = 0.0; 
    bxby_yl[2] = 0.0; 
    bybz_yl[2] = 0.0; 
    byby_yl[3] = 0.0; 
    bxby_yl[3] = 0.0; 
    bybz_yl[3] = 0.0; 
    // If byby, bxby, or bybz < 0.0 at the lower y surface quadrature points, 
    // set cell_avg_magB2_yl to be true in case it was not true before. 
    cell_avg_magB2_yl[0] = 1; 
  } 
 
  if (cell_avg_yr || cell_avg_magB2_yr[0]) { 
    byby_yr[1] = 0.0; 
    bxby_yr[1] = 0.0; 
    bybz_yr[1] = 0.0; 
    byby_yr[2] = 0.0; 
    bxby_yr[2] = 0.0; 
    bybz_yr[2] = 0.0; 
    byby_yr[3] = 0.0; 
    bxby_yr[3] = 0.0; 
    bybz_yr[3] = 0.0; 
    // If byby, bxby, or bybz < 0.0 at the upper y surface quadrature points, 
    // set cell_avg_magB2_yr to be true in case it was not true before. 
    cell_avg_magB2_yr[0] = 1; 
  } 
 
  int cell_avg_zl = 0;
  int cell_avg_zr = 0;
  if (0.5*bzbz_zl[3]-0.5*bzbz_zl[2]-0.5*bzbz_zl[1]+0.5*bzbz_zl[0] < 0.0) cell_avg_zl = 1; 
  if (0.5*bzbz_zr[3]-0.5*bzbz_zr[2]-0.5*bzbz_zr[1]+0.5*bzbz_zr[0] < 0.0) cell_avg_zr = 1; 
  if ((-0.5*bzbz_zl[3])+0.5*bzbz_zl[2]-0.5*bzbz_zl[1]+0.5*bzbz_zl[0] < 0.0) cell_avg_zl = 1; 
  if ((-0.5*bzbz_zr[3])+0.5*bzbz_zr[2]-0.5*bzbz_zr[1]+0.5*bzbz_zr[0] < 0.0) cell_avg_zr = 1; 
  if ((-0.5*bzbz_zl[3])-0.5*bzbz_zl[2]+0.5*bzbz_zl[1]+0.5*bzbz_zl[0] < 0.0) cell_avg_zl = 1; 
  if ((-0.5*bzbz_zr[3])-0.5*bzbz_zr[2]+0.5*bzbz_zr[1]+0.5*bzbz_zr[0] < 0.0) cell_avg_zr = 1; 
  if (0.5*bzbz_zl[3]+0.5*bzbz_zl[2]+0.5*bzbz_zl[1]+0.5*bzbz_zl[0] < 0.0) cell_avg_zl = 1; 
  if (0.5*bzbz_zr[3]+0.5*bzbz_zr[2]+0.5*bzbz_zr[1]+0.5*bzbz_zr[0] < 0.0) cell_avg_zr = 1; 
 
  if (cell_avg_zl || cell_avg_magB2_zl[0]) { 
    bzbz_zl[1] = 0.0; 
    bxbz_zl[1] = 0.0; 
    bybz_zl[1] = 0.0; 
    bzbz_zl[2] = 0.0; 
    bxbz_zl[2] = 0.0; 
    bybz_zl[2] = 0.0; 
    bzbz_zl[3] = 0.0; 
    bxbz_zl[3] = 0.0; 
    bybz_zl[3] = 0.0; 
    // If bzbz, bxbz, or bybz < 0.0 at the lower z surface quadrature points, 
    // set cell_avg_magB2_zl to be true in case it was not true before. 
    cell_avg_magB2_zl[0] = 1; 
  } 
 
  if (cell_avg_zr || cell_avg_magB2_zr[0]) { 
    bzbz_zr[1] = 0.0; 
    bxbz_zr[1] = 0.0; 
    bybz_zr[1] = 0.0; 
    bzbz_zr[2] = 0.0; 
    bxbz_zr[2] = 0.0; 
    bybz_zr[2] = 0.0; 
    bzbz_zr[3] = 0.0; 
    bxbz_zr[3] = 0.0; 
    bybz_zr[3] = 0.0; 
    // If bzbz, bxbz, or bybz < 0.0 at the upper z surface quadrature points, 
    // set cell_avg_magB2_zr to be true in case it was not true before. 
    cell_avg_magB2_zr[0] = 1; 
  } 
 
  // Calculate b_i = B_i/|B| by taking square root of B_i^2/|B|^2 at quadrature points. 
  // Uses the sign of B_i at quadrature points *on the surface* to get the correct sign of b_i *on the surface*. 
  // Note: positivity check already happened, so only uses cell average if needed to avoid imaginary values of b_hat. 
  double B_x_xl[4] = {0.0}; 
  double B_x_xr[4] = {0.0}; 
 
  B_x_xl[0] = 0.7071067811865475*B_x[0]-1.224744871391589*B_x[1]; 
  B_x_xl[1] = 0.7071067811865475*B_x[2]-1.224744871391589*B_x[4]; 
  B_x_xl[2] = 0.7071067811865475*B_x[3]-1.224744871391589*B_x[5]; 
  B_x_xl[3] = 0.7071067811865475*B_x[6]-1.224744871391589*B_x[7]; 
  B_x_xr[0] = 1.224744871391589*B_x[1]+0.7071067811865475*B_x[0]; 
  B_x_xr[1] = 1.224744871391589*B_x[4]+0.7071067811865475*B_x[2]; 
  B_x_xr[2] = 1.224744871391589*B_x[5]+0.7071067811865475*B_x[3]; 
  B_x_xr[3] = 1.224744871391589*B_x[7]+0.7071067811865475*B_x[6]; 
  ser_2x_p1_sqrt_with_sign(B_x_xl, bxbx_xl, bx_xl); 
  ser_2x_p1_sqrt_with_sign(B_x_xr, bxbx_xr, bx_xr); 
 
  double B_y_yl[4] = {0.0}; 
  double B_y_yr[4] = {0.0}; 
 
  B_y_yl[0] = 0.7071067811865475*B_y[0]-1.224744871391589*B_y[2]; 
  B_y_yl[1] = 0.7071067811865475*B_y[1]-1.224744871391589*B_y[4]; 
  B_y_yl[2] = 0.7071067811865475*B_y[3]-1.224744871391589*B_y[6]; 
  B_y_yl[3] = 0.7071067811865475*B_y[5]-1.224744871391589*B_y[7]; 
  B_y_yr[0] = 1.224744871391589*B_y[2]+0.7071067811865475*B_y[0]; 
  B_y_yr[1] = 1.224744871391589*B_y[4]+0.7071067811865475*B_y[1]; 
  B_y_yr[2] = 1.224744871391589*B_y[6]+0.7071067811865475*B_y[3]; 
  B_y_yr[3] = 1.224744871391589*B_y[7]+0.7071067811865475*B_y[5]; 
  ser_2x_p1_sqrt_with_sign(B_y_yl, byby_yl, by_yl); 
  ser_2x_p1_sqrt_with_sign(B_y_yr, byby_yr, by_yr); 
 
  double B_z_zl[4] = {0.0}; 
  double B_z_zr[4] = {0.0}; 
 
  B_z_zl[0] = 0.7071067811865475*B_z[0]-1.224744871391589*B_z[3]; 
  B_z_zl[1] = 0.7071067811865475*B_z[1]-1.224744871391589*B_z[5]; 
  B_z_zl[2] = 0.7071067811865475*B_z[2]-1.224744871391589*B_z[6]; 
  B_z_zl[3] = 0.7071067811865475*B_z[4]-1.224744871391589*B_z[7]; 
  B_z_zr[0] = 1.224744871391589*B_z[3]+0.7071067811865475*B_z[0]; 
  B_z_zr[1] = 1.224744871391589*B_z[5]+0.7071067811865475*B_z[1]; 
  B_z_zr[2] = 1.224744871391589*B_z[6]+0.7071067811865475*B_z[2]; 
  B_z_zr[3] = 1.224744871391589*B_z[7]+0.7071067811865475*B_z[4]; 
  ser_2x_p1_sqrt_with_sign(B_z_zl, bzbz_zl, bz_zl); 
  ser_2x_p1_sqrt_with_sign(B_z_zl, bzbz_zl, bz_zl); 
 
} 
 
