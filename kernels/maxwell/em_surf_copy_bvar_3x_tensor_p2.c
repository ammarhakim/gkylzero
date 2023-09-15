#include <gkyl_mat.h> 
#include <gkyl_maxwell_kernels.h> 
#include <gkyl_basis_tensor_2x_p2_sqrt_with_sign.h> 
GKYL_CU_DH void em_surf_copy_bvar_3x_tensor_p2(int count, struct gkyl_nmat *x, const double *em, int* cell_avg_magB2_surf, double* GKYL_RESTRICT bvar_surf) 
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
  double *bx_xr = &bvar_surf[9]; 
  double *bxbx_xl = &bvar_surf[18]; 
  double *bxbx_xr = &bvar_surf[27]; 
  double *bxby_xl = &bvar_surf[36]; 
  double *bxby_xr = &bvar_surf[45]; 
  double *bxbz_xl = &bvar_surf[54]; 
  double *bxbz_xr = &bvar_surf[63]; 
  int *cell_avg_magB2_xl = &cell_avg_magB2_surf[0]; 
  int *cell_avg_magB2_xr = &cell_avg_magB2_surf[1]; 
 
  struct gkyl_mat x_byby_yl = gkyl_nmat_get(x, count+6); 
  struct gkyl_mat x_byby_yr = gkyl_nmat_get(x, count+7); 
  struct gkyl_mat x_bxby_yl = gkyl_nmat_get(x, count+8); 
  struct gkyl_mat x_bxby_yr = gkyl_nmat_get(x, count+9); 
  struct gkyl_mat x_bybz_yl = gkyl_nmat_get(x, count+10); 
  struct gkyl_mat x_bybz_yr = gkyl_nmat_get(x, count+11); 
 
  double *by_yl = &bvar_surf[72]; 
  double *by_yr = &bvar_surf[81]; 
  double *byby_yl = &bvar_surf[90]; 
  double *byby_yr = &bvar_surf[99]; 
  double *bxby_yl = &bvar_surf[108]; 
  double *bxby_yr = &bvar_surf[117]; 
  double *bybz_yl = &bvar_surf[126]; 
  double *bybz_yr = &bvar_surf[135]; 
  int *cell_avg_magB2_yl = &cell_avg_magB2_surf[2]; 
  int *cell_avg_magB2_yr = &cell_avg_magB2_surf[3]; 
 
  struct gkyl_mat x_bzbz_zl = gkyl_nmat_get(x, count+12); 
  struct gkyl_mat x_bzbz_zr = gkyl_nmat_get(x, count+13); 
  struct gkyl_mat x_bxbz_zl = gkyl_nmat_get(x, count+14); 
  struct gkyl_mat x_bxbz_zr = gkyl_nmat_get(x, count+15); 
  struct gkyl_mat x_bybz_zl = gkyl_nmat_get(x, count+16); 
  struct gkyl_mat x_bybz_zr = gkyl_nmat_get(x, count+17); 
 
  double *bz_zl = &bvar_surf[144]; 
  double *bz_zr = &bvar_surf[153]; 
  double *bzbz_zl = &bvar_surf[162]; 
  double *bzbz_zr = &bvar_surf[171]; 
  double *bxbz_zl = &bvar_surf[180]; 
  double *bxbz_zr = &bvar_surf[189]; 
  double *bybz_zl = &bvar_surf[198]; 
  double *bybz_zr = &bvar_surf[207]; 
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
 
  bxbx_xl[4] = gkyl_mat_get(&x_bxbx_xl,4,0); 
  bxbx_xr[4] = gkyl_mat_get(&x_bxbx_xr,4,0); 
  bxby_xl[4] = gkyl_mat_get(&x_bxby_xl,4,0); 
  bxby_xr[4] = gkyl_mat_get(&x_bxby_xr,4,0); 
  bxbz_xl[4] = gkyl_mat_get(&x_bxbz_xl,4,0); 
  bxbz_xr[4] = gkyl_mat_get(&x_bxbz_xr,4,0); 
 
  byby_yl[4] = gkyl_mat_get(&x_byby_yl,4,0); 
  byby_yr[4] = gkyl_mat_get(&x_byby_yr,4,0); 
  bxby_yl[4] = gkyl_mat_get(&x_bxby_yl,4,0); 
  bxby_yr[4] = gkyl_mat_get(&x_bxby_yr,4,0); 
  bybz_yl[4] = gkyl_mat_get(&x_bybz_yl,4,0); 
  bybz_yr[4] = gkyl_mat_get(&x_bybz_yr,4,0); 
 
  bzbz_zl[4] = gkyl_mat_get(&x_bzbz_zl,4,0); 
  bzbz_zr[4] = gkyl_mat_get(&x_bzbz_zr,4,0); 
  bxbz_zl[4] = gkyl_mat_get(&x_bxbz_zl,4,0); 
  bxbz_zr[4] = gkyl_mat_get(&x_bxbz_zr,4,0); 
  bybz_zl[4] = gkyl_mat_get(&x_bybz_zl,4,0); 
  bybz_zr[4] = gkyl_mat_get(&x_bybz_zr,4,0); 
 
  bxbx_xl[5] = gkyl_mat_get(&x_bxbx_xl,5,0); 
  bxbx_xr[5] = gkyl_mat_get(&x_bxbx_xr,5,0); 
  bxby_xl[5] = gkyl_mat_get(&x_bxby_xl,5,0); 
  bxby_xr[5] = gkyl_mat_get(&x_bxby_xr,5,0); 
  bxbz_xl[5] = gkyl_mat_get(&x_bxbz_xl,5,0); 
  bxbz_xr[5] = gkyl_mat_get(&x_bxbz_xr,5,0); 
 
  byby_yl[5] = gkyl_mat_get(&x_byby_yl,5,0); 
  byby_yr[5] = gkyl_mat_get(&x_byby_yr,5,0); 
  bxby_yl[5] = gkyl_mat_get(&x_bxby_yl,5,0); 
  bxby_yr[5] = gkyl_mat_get(&x_bxby_yr,5,0); 
  bybz_yl[5] = gkyl_mat_get(&x_bybz_yl,5,0); 
  bybz_yr[5] = gkyl_mat_get(&x_bybz_yr,5,0); 
 
  bzbz_zl[5] = gkyl_mat_get(&x_bzbz_zl,5,0); 
  bzbz_zr[5] = gkyl_mat_get(&x_bzbz_zr,5,0); 
  bxbz_zl[5] = gkyl_mat_get(&x_bxbz_zl,5,0); 
  bxbz_zr[5] = gkyl_mat_get(&x_bxbz_zr,5,0); 
  bybz_zl[5] = gkyl_mat_get(&x_bybz_zl,5,0); 
  bybz_zr[5] = gkyl_mat_get(&x_bybz_zr,5,0); 
 
  bxbx_xl[6] = gkyl_mat_get(&x_bxbx_xl,6,0); 
  bxbx_xr[6] = gkyl_mat_get(&x_bxbx_xr,6,0); 
  bxby_xl[6] = gkyl_mat_get(&x_bxby_xl,6,0); 
  bxby_xr[6] = gkyl_mat_get(&x_bxby_xr,6,0); 
  bxbz_xl[6] = gkyl_mat_get(&x_bxbz_xl,6,0); 
  bxbz_xr[6] = gkyl_mat_get(&x_bxbz_xr,6,0); 
 
  byby_yl[6] = gkyl_mat_get(&x_byby_yl,6,0); 
  byby_yr[6] = gkyl_mat_get(&x_byby_yr,6,0); 
  bxby_yl[6] = gkyl_mat_get(&x_bxby_yl,6,0); 
  bxby_yr[6] = gkyl_mat_get(&x_bxby_yr,6,0); 
  bybz_yl[6] = gkyl_mat_get(&x_bybz_yl,6,0); 
  bybz_yr[6] = gkyl_mat_get(&x_bybz_yr,6,0); 
 
  bzbz_zl[6] = gkyl_mat_get(&x_bzbz_zl,6,0); 
  bzbz_zr[6] = gkyl_mat_get(&x_bzbz_zr,6,0); 
  bxbz_zl[6] = gkyl_mat_get(&x_bxbz_zl,6,0); 
  bxbz_zr[6] = gkyl_mat_get(&x_bxbz_zr,6,0); 
  bybz_zl[6] = gkyl_mat_get(&x_bybz_zl,6,0); 
  bybz_zr[6] = gkyl_mat_get(&x_bybz_zr,6,0); 
 
  bxbx_xl[7] = gkyl_mat_get(&x_bxbx_xl,7,0); 
  bxbx_xr[7] = gkyl_mat_get(&x_bxbx_xr,7,0); 
  bxby_xl[7] = gkyl_mat_get(&x_bxby_xl,7,0); 
  bxby_xr[7] = gkyl_mat_get(&x_bxby_xr,7,0); 
  bxbz_xl[7] = gkyl_mat_get(&x_bxbz_xl,7,0); 
  bxbz_xr[7] = gkyl_mat_get(&x_bxbz_xr,7,0); 
 
  byby_yl[7] = gkyl_mat_get(&x_byby_yl,7,0); 
  byby_yr[7] = gkyl_mat_get(&x_byby_yr,7,0); 
  bxby_yl[7] = gkyl_mat_get(&x_bxby_yl,7,0); 
  bxby_yr[7] = gkyl_mat_get(&x_bxby_yr,7,0); 
  bybz_yl[7] = gkyl_mat_get(&x_bybz_yl,7,0); 
  bybz_yr[7] = gkyl_mat_get(&x_bybz_yr,7,0); 
 
  bzbz_zl[7] = gkyl_mat_get(&x_bzbz_zl,7,0); 
  bzbz_zr[7] = gkyl_mat_get(&x_bzbz_zr,7,0); 
  bxbz_zl[7] = gkyl_mat_get(&x_bxbz_zl,7,0); 
  bxbz_zr[7] = gkyl_mat_get(&x_bxbz_zr,7,0); 
  bybz_zl[7] = gkyl_mat_get(&x_bybz_zl,7,0); 
  bybz_zr[7] = gkyl_mat_get(&x_bybz_zr,7,0); 
 
  bxbx_xl[8] = gkyl_mat_get(&x_bxbx_xl,8,0); 
  bxbx_xr[8] = gkyl_mat_get(&x_bxbx_xr,8,0); 
  bxby_xl[8] = gkyl_mat_get(&x_bxby_xl,8,0); 
  bxby_xr[8] = gkyl_mat_get(&x_bxby_xr,8,0); 
  bxbz_xl[8] = gkyl_mat_get(&x_bxbz_xl,8,0); 
  bxbz_xr[8] = gkyl_mat_get(&x_bxbz_xr,8,0); 
 
  byby_yl[8] = gkyl_mat_get(&x_byby_yl,8,0); 
  byby_yr[8] = gkyl_mat_get(&x_byby_yr,8,0); 
  bxby_yl[8] = gkyl_mat_get(&x_bxby_yl,8,0); 
  bxby_yr[8] = gkyl_mat_get(&x_bxby_yr,8,0); 
  bybz_yl[8] = gkyl_mat_get(&x_bybz_yl,8,0); 
  bybz_yr[8] = gkyl_mat_get(&x_bybz_yr,8,0); 
 
  bzbz_zl[8] = gkyl_mat_get(&x_bzbz_zl,8,0); 
  bzbz_zr[8] = gkyl_mat_get(&x_bzbz_zr,8,0); 
  bxbz_zl[8] = gkyl_mat_get(&x_bxbz_zl,8,0); 
  bxbz_zr[8] = gkyl_mat_get(&x_bxbz_zr,8,0); 
  bybz_zl[8] = gkyl_mat_get(&x_bybz_zl,8,0); 
  bybz_zr[8] = gkyl_mat_get(&x_bybz_zr,8,0); 
 
  const double *B_x = &em[81]; 
  const double *B_y = &em[108]; 
  const double *B_z = &em[135]; 
 
  int cell_avg_xl = 0;
  int cell_avg_xr = 0;
  if (0.4*bxbx_xl[8]-0.5999999999999995*bxbx_xl[7]-0.5999999999999999*bxbx_xl[6]+0.4472135954999579*bxbx_xl[5]+0.4472135954999579*bxbx_xl[4]+0.9*bxbx_xl[3]-0.6708203932499369*bxbx_xl[2]-0.6708203932499369*bxbx_xl[1]+0.5*bxbx_xl[0] < 0.0) cell_avg_xl = 1; 
  if (0.4*bxbx_xr[8]-0.5999999999999995*bxbx_xr[7]-0.5999999999999999*bxbx_xr[6]+0.4472135954999579*bxbx_xr[5]+0.4472135954999579*bxbx_xr[4]+0.9*bxbx_xr[3]-0.6708203932499369*bxbx_xr[2]-0.6708203932499369*bxbx_xr[1]+0.5*bxbx_xr[0] < 0.0) cell_avg_xr = 1; 
  if ((-0.5*bxbx_xl[8])+0.75*bxbx_xl[7]-0.5590169943749475*bxbx_xl[5]+0.4472135954999579*bxbx_xl[4]-0.6708203932499369*bxbx_xl[1]+0.5*bxbx_xl[0] < 0.0) cell_avg_xl = 1; 
  if ((-0.5*bxbx_xr[8])+0.75*bxbx_xr[7]-0.5590169943749475*bxbx_xr[5]+0.4472135954999579*bxbx_xr[4]-0.6708203932499369*bxbx_xr[1]+0.5*bxbx_xr[0] < 0.0) cell_avg_xr = 1; 
  if (0.4*bxbx_xl[8]-0.5999999999999995*bxbx_xl[7]+0.5999999999999999*bxbx_xl[6]+0.4472135954999579*bxbx_xl[5]+0.4472135954999579*bxbx_xl[4]-0.9*bxbx_xl[3]+0.6708203932499369*bxbx_xl[2]-0.6708203932499369*bxbx_xl[1]+0.5*bxbx_xl[0] < 0.0) cell_avg_xl = 1; 
  if (0.4*bxbx_xr[8]-0.5999999999999995*bxbx_xr[7]+0.5999999999999999*bxbx_xr[6]+0.4472135954999579*bxbx_xr[5]+0.4472135954999579*bxbx_xr[4]-0.9*bxbx_xr[3]+0.6708203932499369*bxbx_xr[2]-0.6708203932499369*bxbx_xr[1]+0.5*bxbx_xr[0] < 0.0) cell_avg_xr = 1; 
  if ((-0.5*bxbx_xl[8])+0.75*bxbx_xl[6]+0.4472135954999579*bxbx_xl[5]-0.5590169943749475*bxbx_xl[4]-0.6708203932499369*bxbx_xl[2]+0.5*bxbx_xl[0] < 0.0) cell_avg_xl = 1; 
  if ((-0.5*bxbx_xr[8])+0.75*bxbx_xr[6]+0.4472135954999579*bxbx_xr[5]-0.5590169943749475*bxbx_xr[4]-0.6708203932499369*bxbx_xr[2]+0.5*bxbx_xr[0] < 0.0) cell_avg_xr = 1; 
  if (0.625*bxbx_xl[8]-0.5590169943749475*bxbx_xl[5]-0.5590169943749475*bxbx_xl[4]+0.5*bxbx_xl[0] < 0.0) cell_avg_xl = 1; 
  if (0.625*bxbx_xr[8]-0.5590169943749475*bxbx_xr[5]-0.5590169943749475*bxbx_xr[4]+0.5*bxbx_xr[0] < 0.0) cell_avg_xr = 1; 
  if ((-0.5*bxbx_xl[8])-0.75*bxbx_xl[6]+0.4472135954999579*bxbx_xl[5]-0.5590169943749475*bxbx_xl[4]+0.6708203932499369*bxbx_xl[2]+0.5*bxbx_xl[0] < 0.0) cell_avg_xl = 1; 
  if ((-0.5*bxbx_xr[8])-0.75*bxbx_xr[6]+0.4472135954999579*bxbx_xr[5]-0.5590169943749475*bxbx_xr[4]+0.6708203932499369*bxbx_xr[2]+0.5*bxbx_xr[0] < 0.0) cell_avg_xr = 1; 
  if (0.4*bxbx_xl[8]+0.5999999999999995*bxbx_xl[7]-0.5999999999999999*bxbx_xl[6]+0.4472135954999579*bxbx_xl[5]+0.4472135954999579*bxbx_xl[4]-0.9*bxbx_xl[3]-0.6708203932499369*bxbx_xl[2]+0.6708203932499369*bxbx_xl[1]+0.5*bxbx_xl[0] < 0.0) cell_avg_xl = 1; 
  if (0.4*bxbx_xr[8]+0.5999999999999995*bxbx_xr[7]-0.5999999999999999*bxbx_xr[6]+0.4472135954999579*bxbx_xr[5]+0.4472135954999579*bxbx_xr[4]-0.9*bxbx_xr[3]-0.6708203932499369*bxbx_xr[2]+0.6708203932499369*bxbx_xr[1]+0.5*bxbx_xr[0] < 0.0) cell_avg_xr = 1; 
  if ((-0.5*bxbx_xl[8])-0.75*bxbx_xl[7]-0.5590169943749475*bxbx_xl[5]+0.4472135954999579*bxbx_xl[4]+0.6708203932499369*bxbx_xl[1]+0.5*bxbx_xl[0] < 0.0) cell_avg_xl = 1; 
  if ((-0.5*bxbx_xr[8])-0.75*bxbx_xr[7]-0.5590169943749475*bxbx_xr[5]+0.4472135954999579*bxbx_xr[4]+0.6708203932499369*bxbx_xr[1]+0.5*bxbx_xr[0] < 0.0) cell_avg_xr = 1; 
  if (0.4*bxbx_xl[8]+0.5999999999999995*bxbx_xl[7]+0.5999999999999999*bxbx_xl[6]+0.4472135954999579*bxbx_xl[5]+0.4472135954999579*bxbx_xl[4]+0.9*bxbx_xl[3]+0.6708203932499369*bxbx_xl[2]+0.6708203932499369*bxbx_xl[1]+0.5*bxbx_xl[0] < 0.0) cell_avg_xl = 1; 
  if (0.4*bxbx_xr[8]+0.5999999999999995*bxbx_xr[7]+0.5999999999999999*bxbx_xr[6]+0.4472135954999579*bxbx_xr[5]+0.4472135954999579*bxbx_xr[4]+0.9*bxbx_xr[3]+0.6708203932499369*bxbx_xr[2]+0.6708203932499369*bxbx_xr[1]+0.5*bxbx_xr[0] < 0.0) cell_avg_xr = 1; 
 
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
    bxbx_xl[4] = 0.0; 
    bxby_xl[4] = 0.0; 
    bxbz_xl[4] = 0.0; 
    bxbx_xl[5] = 0.0; 
    bxby_xl[5] = 0.0; 
    bxbz_xl[5] = 0.0; 
    bxbx_xl[6] = 0.0; 
    bxby_xl[6] = 0.0; 
    bxbz_xl[6] = 0.0; 
    bxbx_xl[7] = 0.0; 
    bxby_xl[7] = 0.0; 
    bxbz_xl[7] = 0.0; 
    bxbx_xl[8] = 0.0; 
    bxby_xl[8] = 0.0; 
    bxbz_xl[8] = 0.0; 
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
    bxbx_xr[4] = 0.0; 
    bxby_xr[4] = 0.0; 
    bxbz_xr[4] = 0.0; 
    bxbx_xr[5] = 0.0; 
    bxby_xr[5] = 0.0; 
    bxbz_xr[5] = 0.0; 
    bxbx_xr[6] = 0.0; 
    bxby_xr[6] = 0.0; 
    bxbz_xr[6] = 0.0; 
    bxbx_xr[7] = 0.0; 
    bxby_xr[7] = 0.0; 
    bxbz_xr[7] = 0.0; 
    bxbx_xr[8] = 0.0; 
    bxby_xr[8] = 0.0; 
    bxbz_xr[8] = 0.0; 
    // If bxbx, bxby, or bxbz < 0.0 at the upper x surface quadrature points, 
    // set cell_avg_magB2_xr to be true in case it was not true before. 
    cell_avg_magB2_xr[0] = 1; 
  } 
 
  int cell_avg_yl = 0;
  int cell_avg_yr = 0;
  if (0.4*byby_yl[8]-0.5999999999999995*byby_yl[7]-0.5999999999999999*byby_yl[6]+0.4472135954999579*byby_yl[5]+0.4472135954999579*byby_yl[4]+0.9*byby_yl[3]-0.6708203932499369*byby_yl[2]-0.6708203932499369*byby_yl[1]+0.5*byby_yl[0] < 0.0) cell_avg_yl = 1; 
  if (0.4*byby_yr[8]-0.5999999999999995*byby_yr[7]-0.5999999999999999*byby_yr[6]+0.4472135954999579*byby_yr[5]+0.4472135954999579*byby_yr[4]+0.9*byby_yr[3]-0.6708203932499369*byby_yr[2]-0.6708203932499369*byby_yr[1]+0.5*byby_yr[0] < 0.0) cell_avg_yr = 1; 
  if ((-0.5*byby_yl[8])+0.75*byby_yl[7]-0.5590169943749475*byby_yl[5]+0.4472135954999579*byby_yl[4]-0.6708203932499369*byby_yl[1]+0.5*byby_yl[0] < 0.0) cell_avg_yl = 1; 
  if ((-0.5*byby_yr[8])+0.75*byby_yr[7]-0.5590169943749475*byby_yr[5]+0.4472135954999579*byby_yr[4]-0.6708203932499369*byby_yr[1]+0.5*byby_yr[0] < 0.0) cell_avg_yr = 1; 
  if (0.4*byby_yl[8]-0.5999999999999995*byby_yl[7]+0.5999999999999999*byby_yl[6]+0.4472135954999579*byby_yl[5]+0.4472135954999579*byby_yl[4]-0.9*byby_yl[3]+0.6708203932499369*byby_yl[2]-0.6708203932499369*byby_yl[1]+0.5*byby_yl[0] < 0.0) cell_avg_yl = 1; 
  if (0.4*byby_yr[8]-0.5999999999999995*byby_yr[7]+0.5999999999999999*byby_yr[6]+0.4472135954999579*byby_yr[5]+0.4472135954999579*byby_yr[4]-0.9*byby_yr[3]+0.6708203932499369*byby_yr[2]-0.6708203932499369*byby_yr[1]+0.5*byby_yr[0] < 0.0) cell_avg_yr = 1; 
  if ((-0.5*byby_yl[8])+0.75*byby_yl[6]+0.4472135954999579*byby_yl[5]-0.5590169943749475*byby_yl[4]-0.6708203932499369*byby_yl[2]+0.5*byby_yl[0] < 0.0) cell_avg_yl = 1; 
  if ((-0.5*byby_yr[8])+0.75*byby_yr[6]+0.4472135954999579*byby_yr[5]-0.5590169943749475*byby_yr[4]-0.6708203932499369*byby_yr[2]+0.5*byby_yr[0] < 0.0) cell_avg_yr = 1; 
  if (0.625*byby_yl[8]-0.5590169943749475*byby_yl[5]-0.5590169943749475*byby_yl[4]+0.5*byby_yl[0] < 0.0) cell_avg_yl = 1; 
  if (0.625*byby_yr[8]-0.5590169943749475*byby_yr[5]-0.5590169943749475*byby_yr[4]+0.5*byby_yr[0] < 0.0) cell_avg_yr = 1; 
  if ((-0.5*byby_yl[8])-0.75*byby_yl[6]+0.4472135954999579*byby_yl[5]-0.5590169943749475*byby_yl[4]+0.6708203932499369*byby_yl[2]+0.5*byby_yl[0] < 0.0) cell_avg_yl = 1; 
  if ((-0.5*byby_yr[8])-0.75*byby_yr[6]+0.4472135954999579*byby_yr[5]-0.5590169943749475*byby_yr[4]+0.6708203932499369*byby_yr[2]+0.5*byby_yr[0] < 0.0) cell_avg_yr = 1; 
  if (0.4*byby_yl[8]+0.5999999999999995*byby_yl[7]-0.5999999999999999*byby_yl[6]+0.4472135954999579*byby_yl[5]+0.4472135954999579*byby_yl[4]-0.9*byby_yl[3]-0.6708203932499369*byby_yl[2]+0.6708203932499369*byby_yl[1]+0.5*byby_yl[0] < 0.0) cell_avg_yl = 1; 
  if (0.4*byby_yr[8]+0.5999999999999995*byby_yr[7]-0.5999999999999999*byby_yr[6]+0.4472135954999579*byby_yr[5]+0.4472135954999579*byby_yr[4]-0.9*byby_yr[3]-0.6708203932499369*byby_yr[2]+0.6708203932499369*byby_yr[1]+0.5*byby_yr[0] < 0.0) cell_avg_yr = 1; 
  if ((-0.5*byby_yl[8])-0.75*byby_yl[7]-0.5590169943749475*byby_yl[5]+0.4472135954999579*byby_yl[4]+0.6708203932499369*byby_yl[1]+0.5*byby_yl[0] < 0.0) cell_avg_yl = 1; 
  if ((-0.5*byby_yr[8])-0.75*byby_yr[7]-0.5590169943749475*byby_yr[5]+0.4472135954999579*byby_yr[4]+0.6708203932499369*byby_yr[1]+0.5*byby_yr[0] < 0.0) cell_avg_yr = 1; 
  if (0.4*byby_yl[8]+0.5999999999999995*byby_yl[7]+0.5999999999999999*byby_yl[6]+0.4472135954999579*byby_yl[5]+0.4472135954999579*byby_yl[4]+0.9*byby_yl[3]+0.6708203932499369*byby_yl[2]+0.6708203932499369*byby_yl[1]+0.5*byby_yl[0] < 0.0) cell_avg_yl = 1; 
  if (0.4*byby_yr[8]+0.5999999999999995*byby_yr[7]+0.5999999999999999*byby_yr[6]+0.4472135954999579*byby_yr[5]+0.4472135954999579*byby_yr[4]+0.9*byby_yr[3]+0.6708203932499369*byby_yr[2]+0.6708203932499369*byby_yr[1]+0.5*byby_yr[0] < 0.0) cell_avg_yr = 1; 
 
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
    byby_yl[4] = 0.0; 
    bxby_yl[4] = 0.0; 
    bybz_yl[4] = 0.0; 
    byby_yl[5] = 0.0; 
    bxby_yl[5] = 0.0; 
    bybz_yl[5] = 0.0; 
    byby_yl[6] = 0.0; 
    bxby_yl[6] = 0.0; 
    bybz_yl[6] = 0.0; 
    byby_yl[7] = 0.0; 
    bxby_yl[7] = 0.0; 
    bybz_yl[7] = 0.0; 
    byby_yl[8] = 0.0; 
    bxby_yl[8] = 0.0; 
    bybz_yl[8] = 0.0; 
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
    byby_yr[4] = 0.0; 
    bxby_yr[4] = 0.0; 
    bybz_yr[4] = 0.0; 
    byby_yr[5] = 0.0; 
    bxby_yr[5] = 0.0; 
    bybz_yr[5] = 0.0; 
    byby_yr[6] = 0.0; 
    bxby_yr[6] = 0.0; 
    bybz_yr[6] = 0.0; 
    byby_yr[7] = 0.0; 
    bxby_yr[7] = 0.0; 
    bybz_yr[7] = 0.0; 
    byby_yr[8] = 0.0; 
    bxby_yr[8] = 0.0; 
    bybz_yr[8] = 0.0; 
    // If byby, bxby, or bybz < 0.0 at the upper y surface quadrature points, 
    // set cell_avg_magB2_yr to be true in case it was not true before. 
    cell_avg_magB2_yr[0] = 1; 
  } 
 
  int cell_avg_zl = 0;
  int cell_avg_zr = 0;
  if (0.4*bzbz_zl[8]-0.5999999999999995*bzbz_zl[7]-0.5999999999999999*bzbz_zl[6]+0.4472135954999579*bzbz_zl[5]+0.4472135954999579*bzbz_zl[4]+0.9*bzbz_zl[3]-0.6708203932499369*bzbz_zl[2]-0.6708203932499369*bzbz_zl[1]+0.5*bzbz_zl[0] < 0.0) cell_avg_zl = 1; 
  if (0.4*bzbz_zr[8]-0.5999999999999995*bzbz_zr[7]-0.5999999999999999*bzbz_zr[6]+0.4472135954999579*bzbz_zr[5]+0.4472135954999579*bzbz_zr[4]+0.9*bzbz_zr[3]-0.6708203932499369*bzbz_zr[2]-0.6708203932499369*bzbz_zr[1]+0.5*bzbz_zr[0] < 0.0) cell_avg_zr = 1; 
  if ((-0.5*bzbz_zl[8])+0.75*bzbz_zl[7]-0.5590169943749475*bzbz_zl[5]+0.4472135954999579*bzbz_zl[4]-0.6708203932499369*bzbz_zl[1]+0.5*bzbz_zl[0] < 0.0) cell_avg_zl = 1; 
  if ((-0.5*bzbz_zr[8])+0.75*bzbz_zr[7]-0.5590169943749475*bzbz_zr[5]+0.4472135954999579*bzbz_zr[4]-0.6708203932499369*bzbz_zr[1]+0.5*bzbz_zr[0] < 0.0) cell_avg_zr = 1; 
  if (0.4*bzbz_zl[8]-0.5999999999999995*bzbz_zl[7]+0.5999999999999999*bzbz_zl[6]+0.4472135954999579*bzbz_zl[5]+0.4472135954999579*bzbz_zl[4]-0.9*bzbz_zl[3]+0.6708203932499369*bzbz_zl[2]-0.6708203932499369*bzbz_zl[1]+0.5*bzbz_zl[0] < 0.0) cell_avg_zl = 1; 
  if (0.4*bzbz_zr[8]-0.5999999999999995*bzbz_zr[7]+0.5999999999999999*bzbz_zr[6]+0.4472135954999579*bzbz_zr[5]+0.4472135954999579*bzbz_zr[4]-0.9*bzbz_zr[3]+0.6708203932499369*bzbz_zr[2]-0.6708203932499369*bzbz_zr[1]+0.5*bzbz_zr[0] < 0.0) cell_avg_zr = 1; 
  if ((-0.5*bzbz_zl[8])+0.75*bzbz_zl[6]+0.4472135954999579*bzbz_zl[5]-0.5590169943749475*bzbz_zl[4]-0.6708203932499369*bzbz_zl[2]+0.5*bzbz_zl[0] < 0.0) cell_avg_zl = 1; 
  if ((-0.5*bzbz_zr[8])+0.75*bzbz_zr[6]+0.4472135954999579*bzbz_zr[5]-0.5590169943749475*bzbz_zr[4]-0.6708203932499369*bzbz_zr[2]+0.5*bzbz_zr[0] < 0.0) cell_avg_zr = 1; 
  if (0.625*bzbz_zl[8]-0.5590169943749475*bzbz_zl[5]-0.5590169943749475*bzbz_zl[4]+0.5*bzbz_zl[0] < 0.0) cell_avg_zl = 1; 
  if (0.625*bzbz_zr[8]-0.5590169943749475*bzbz_zr[5]-0.5590169943749475*bzbz_zr[4]+0.5*bzbz_zr[0] < 0.0) cell_avg_zr = 1; 
  if ((-0.5*bzbz_zl[8])-0.75*bzbz_zl[6]+0.4472135954999579*bzbz_zl[5]-0.5590169943749475*bzbz_zl[4]+0.6708203932499369*bzbz_zl[2]+0.5*bzbz_zl[0] < 0.0) cell_avg_zl = 1; 
  if ((-0.5*bzbz_zr[8])-0.75*bzbz_zr[6]+0.4472135954999579*bzbz_zr[5]-0.5590169943749475*bzbz_zr[4]+0.6708203932499369*bzbz_zr[2]+0.5*bzbz_zr[0] < 0.0) cell_avg_zr = 1; 
  if (0.4*bzbz_zl[8]+0.5999999999999995*bzbz_zl[7]-0.5999999999999999*bzbz_zl[6]+0.4472135954999579*bzbz_zl[5]+0.4472135954999579*bzbz_zl[4]-0.9*bzbz_zl[3]-0.6708203932499369*bzbz_zl[2]+0.6708203932499369*bzbz_zl[1]+0.5*bzbz_zl[0] < 0.0) cell_avg_zl = 1; 
  if (0.4*bzbz_zr[8]+0.5999999999999995*bzbz_zr[7]-0.5999999999999999*bzbz_zr[6]+0.4472135954999579*bzbz_zr[5]+0.4472135954999579*bzbz_zr[4]-0.9*bzbz_zr[3]-0.6708203932499369*bzbz_zr[2]+0.6708203932499369*bzbz_zr[1]+0.5*bzbz_zr[0] < 0.0) cell_avg_zr = 1; 
  if ((-0.5*bzbz_zl[8])-0.75*bzbz_zl[7]-0.5590169943749475*bzbz_zl[5]+0.4472135954999579*bzbz_zl[4]+0.6708203932499369*bzbz_zl[1]+0.5*bzbz_zl[0] < 0.0) cell_avg_zl = 1; 
  if ((-0.5*bzbz_zr[8])-0.75*bzbz_zr[7]-0.5590169943749475*bzbz_zr[5]+0.4472135954999579*bzbz_zr[4]+0.6708203932499369*bzbz_zr[1]+0.5*bzbz_zr[0] < 0.0) cell_avg_zr = 1; 
  if (0.4*bzbz_zl[8]+0.5999999999999995*bzbz_zl[7]+0.5999999999999999*bzbz_zl[6]+0.4472135954999579*bzbz_zl[5]+0.4472135954999579*bzbz_zl[4]+0.9*bzbz_zl[3]+0.6708203932499369*bzbz_zl[2]+0.6708203932499369*bzbz_zl[1]+0.5*bzbz_zl[0] < 0.0) cell_avg_zl = 1; 
  if (0.4*bzbz_zr[8]+0.5999999999999995*bzbz_zr[7]+0.5999999999999999*bzbz_zr[6]+0.4472135954999579*bzbz_zr[5]+0.4472135954999579*bzbz_zr[4]+0.9*bzbz_zr[3]+0.6708203932499369*bzbz_zr[2]+0.6708203932499369*bzbz_zr[1]+0.5*bzbz_zr[0] < 0.0) cell_avg_zr = 1; 
 
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
    bzbz_zl[4] = 0.0; 
    bxbz_zl[4] = 0.0; 
    bybz_zl[4] = 0.0; 
    bzbz_zl[5] = 0.0; 
    bxbz_zl[5] = 0.0; 
    bybz_zl[5] = 0.0; 
    bzbz_zl[6] = 0.0; 
    bxbz_zl[6] = 0.0; 
    bybz_zl[6] = 0.0; 
    bzbz_zl[7] = 0.0; 
    bxbz_zl[7] = 0.0; 
    bybz_zl[7] = 0.0; 
    bzbz_zl[8] = 0.0; 
    bxbz_zl[8] = 0.0; 
    bybz_zl[8] = 0.0; 
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
    bzbz_zr[4] = 0.0; 
    bxbz_zr[4] = 0.0; 
    bybz_zr[4] = 0.0; 
    bzbz_zr[5] = 0.0; 
    bxbz_zr[5] = 0.0; 
    bybz_zr[5] = 0.0; 
    bzbz_zr[6] = 0.0; 
    bxbz_zr[6] = 0.0; 
    bybz_zr[6] = 0.0; 
    bzbz_zr[7] = 0.0; 
    bxbz_zr[7] = 0.0; 
    bybz_zr[7] = 0.0; 
    bzbz_zr[8] = 0.0; 
    bxbz_zr[8] = 0.0; 
    bybz_zr[8] = 0.0; 
    // If bzbz, bxbz, or bybz < 0.0 at the upper z surface quadrature points, 
    // set cell_avg_magB2_zr to be true in case it was not true before. 
    cell_avg_magB2_zr[0] = 1; 
  } 
 
  // Calculate b_i = B_i/|B| by taking square root of B_i^2/|B|^2 at quadrature points. 
  // Uses the sign of B_i at quadrature points *on the surface* to get the correct sign of b_i *on the surface*. 
  // Note: positivity check already happened, so only uses cell average if needed to avoid imaginary values of b_hat. 
  double B_x_xl[9] = {0.0}; 
  double B_x_xr[9] = {0.0}; 
 
  B_x_xl[0] = 1.58113883008419*B_x[7]-1.224744871391589*B_x[1]+0.7071067811865475*B_x[0]; 
  B_x_xl[1] = 1.58113883008419*B_x[11]-1.224744871391589*B_x[4]+0.7071067811865475*B_x[2]; 
  B_x_xl[2] = 1.58113883008419*B_x[13]-1.224744871391589*B_x[5]+0.7071067811865475*B_x[3]; 
  B_x_xl[3] = 1.58113883008419*B_x[17]-1.224744871391589*B_x[10]+0.7071067811865475*B_x[6]; 
  B_x_xl[4] = 1.58113883008419*B_x[20]-1.224744871391589*B_x[12]+0.7071067811865475*B_x[8]; 
  B_x_xl[5] = 1.58113883008419*B_x[21]-1.224744871391589*B_x[15]+0.7071067811865475*B_x[9]; 
  B_x_xl[6] = 1.58113883008419*B_x[23]-1.224744871391589*B_x[18]+0.7071067811865475*B_x[14]; 
  B_x_xl[7] = 1.58113883008419*B_x[24]-1.224744871391589*B_x[19]+0.7071067811865475*B_x[16]; 
  B_x_xl[8] = 1.58113883008419*B_x[26]-1.224744871391589*B_x[25]+0.7071067811865475*B_x[22]; 
  B_x_xr[0] = 1.58113883008419*B_x[7]+1.224744871391589*B_x[1]+0.7071067811865475*B_x[0]; 
  B_x_xr[1] = 1.58113883008419*B_x[11]+1.224744871391589*B_x[4]+0.7071067811865475*B_x[2]; 
  B_x_xr[2] = 1.58113883008419*B_x[13]+1.224744871391589*B_x[5]+0.7071067811865475*B_x[3]; 
  B_x_xr[3] = 1.58113883008419*B_x[17]+1.224744871391589*B_x[10]+0.7071067811865475*B_x[6]; 
  B_x_xr[4] = 1.58113883008419*B_x[20]+1.224744871391589*B_x[12]+0.7071067811865475*B_x[8]; 
  B_x_xr[5] = 1.58113883008419*B_x[21]+1.224744871391589*B_x[15]+0.7071067811865475*B_x[9]; 
  B_x_xr[6] = 1.58113883008419*B_x[23]+1.224744871391589*B_x[18]+0.7071067811865475*B_x[14]; 
  B_x_xr[7] = 1.58113883008419*B_x[24]+1.224744871391589*B_x[19]+0.7071067811865475*B_x[16]; 
  B_x_xr[8] = 1.58113883008419*B_x[26]+1.224744871391589*B_x[25]+0.7071067811865475*B_x[22]; 
  tensor_2x_p2_sqrt_with_sign(B_x_xl, bxbx_xl, bx_xl); 
  tensor_2x_p2_sqrt_with_sign(B_x_xr, bxbx_xr, bx_xr); 
 
  double B_y_yl[9] = {0.0}; 
  double B_y_yr[9] = {0.0}; 
 
  B_y_yl[0] = 1.58113883008419*B_y[8]-1.224744871391589*B_y[2]+0.7071067811865475*B_y[0]; 
  B_y_yl[1] = 1.58113883008419*B_y[12]-1.224744871391589*B_y[4]+0.7071067811865475*B_y[1]; 
  B_y_yl[2] = 1.58113883008419*B_y[14]-1.224744871391589*B_y[6]+0.7071067811865475*B_y[3]; 
  B_y_yl[3] = 1.58113883008419*B_y[18]-1.224744871391589*B_y[10]+0.7071067811865475*B_y[5]; 
  B_y_yl[4] = 1.58113883008419*B_y[20]-1.224744871391589*B_y[11]+0.7071067811865475*B_y[7]; 
  B_y_yl[5] = 1.58113883008419*B_y[22]-1.224744871391589*B_y[16]+0.7071067811865475*B_y[9]; 
  B_y_yl[6] = 1.58113883008419*B_y[23]-1.224744871391589*B_y[17]+0.7071067811865475*B_y[13]; 
  B_y_yl[7] = 1.58113883008419*B_y[25]-1.224744871391589*B_y[19]+0.7071067811865475*B_y[15]; 
  B_y_yl[8] = 1.58113883008419*B_y[26]-1.224744871391589*B_y[24]+0.7071067811865475*B_y[21]; 
  B_y_yr[0] = 1.58113883008419*B_y[8]+1.224744871391589*B_y[2]+0.7071067811865475*B_y[0]; 
  B_y_yr[1] = 1.58113883008419*B_y[12]+1.224744871391589*B_y[4]+0.7071067811865475*B_y[1]; 
  B_y_yr[2] = 1.58113883008419*B_y[14]+1.224744871391589*B_y[6]+0.7071067811865475*B_y[3]; 
  B_y_yr[3] = 1.58113883008419*B_y[18]+1.224744871391589*B_y[10]+0.7071067811865475*B_y[5]; 
  B_y_yr[4] = 1.58113883008419*B_y[20]+1.224744871391589*B_y[11]+0.7071067811865475*B_y[7]; 
  B_y_yr[5] = 1.58113883008419*B_y[22]+1.224744871391589*B_y[16]+0.7071067811865475*B_y[9]; 
  B_y_yr[6] = 1.58113883008419*B_y[23]+1.224744871391589*B_y[17]+0.7071067811865475*B_y[13]; 
  B_y_yr[7] = 1.58113883008419*B_y[25]+1.224744871391589*B_y[19]+0.7071067811865475*B_y[15]; 
  B_y_yr[8] = 1.58113883008419*B_y[26]+1.224744871391589*B_y[24]+0.7071067811865475*B_y[21]; 
  tensor_2x_p2_sqrt_with_sign(B_y_yl, byby_yl, by_yl); 
  tensor_2x_p2_sqrt_with_sign(B_y_yr, byby_yr, by_yr); 
 
  double B_z_zl[9] = {0.0}; 
  double B_z_zr[9] = {0.0}; 
 
  B_z_zl[0] = 1.58113883008419*B_z[9]-1.224744871391589*B_z[3]+0.7071067811865475*B_z[0]; 
  B_z_zl[1] = 1.58113883008419*B_z[15]-1.224744871391589*B_z[5]+0.7071067811865475*B_z[1]; 
  B_z_zl[2] = 1.58113883008419*B_z[16]-1.224744871391589*B_z[6]+0.7071067811865475*B_z[2]; 
  B_z_zl[3] = 1.58113883008419*B_z[19]-1.224744871391589*B_z[10]+0.7071067811865475*B_z[4]; 
  B_z_zl[4] = 1.58113883008419*B_z[21]-1.224744871391589*B_z[13]+0.7071067811865475*B_z[7]; 
  B_z_zl[5] = 1.58113883008419*B_z[22]-1.224744871391589*B_z[14]+0.7071067811865475*B_z[8]; 
  B_z_zl[6] = 1.58113883008419*B_z[24]-1.224744871391589*B_z[17]+0.7071067811865475*B_z[11]; 
  B_z_zl[7] = 1.58113883008419*B_z[25]-1.224744871391589*B_z[18]+0.7071067811865475*B_z[12]; 
  B_z_zl[8] = 1.58113883008419*B_z[26]-1.224744871391589*B_z[23]+0.7071067811865475*B_z[20]; 
  B_z_zr[0] = 1.58113883008419*B_z[9]+1.224744871391589*B_z[3]+0.7071067811865475*B_z[0]; 
  B_z_zr[1] = 1.58113883008419*B_z[15]+1.224744871391589*B_z[5]+0.7071067811865475*B_z[1]; 
  B_z_zr[2] = 1.58113883008419*B_z[16]+1.224744871391589*B_z[6]+0.7071067811865475*B_z[2]; 
  B_z_zr[3] = 1.58113883008419*B_z[19]+1.224744871391589*B_z[10]+0.7071067811865475*B_z[4]; 
  B_z_zr[4] = 1.58113883008419*B_z[21]+1.224744871391589*B_z[13]+0.7071067811865475*B_z[7]; 
  B_z_zr[5] = 1.58113883008419*B_z[22]+1.224744871391589*B_z[14]+0.7071067811865475*B_z[8]; 
  B_z_zr[6] = 1.58113883008419*B_z[24]+1.224744871391589*B_z[17]+0.7071067811865475*B_z[11]; 
  B_z_zr[7] = 1.58113883008419*B_z[25]+1.224744871391589*B_z[18]+0.7071067811865475*B_z[12]; 
  B_z_zr[8] = 1.58113883008419*B_z[26]+1.224744871391589*B_z[23]+0.7071067811865475*B_z[20]; 
  tensor_2x_p2_sqrt_with_sign(B_z_zl, bzbz_zl, bz_zl); 
  tensor_2x_p2_sqrt_with_sign(B_z_zl, bzbz_zl, bz_zl); 
 
} 
 
