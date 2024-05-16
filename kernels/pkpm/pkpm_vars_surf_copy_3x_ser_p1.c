#include <gkyl_mat.h> 
#include <gkyl_euler_pkpm_kernels.h> 
GKYL_CU_DH void pkpm_vars_surf_copy_3x_ser_p1(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT prim_surf) 
{ 
  // count: integer to indicate which matrix being fetched. 
  // x:     Input solution vector. 
  // prim_surf: Primitive variables at surfaces.
 
  struct gkyl_mat x_ux_xl = gkyl_nmat_get(x, count); 
  struct gkyl_mat x_ux_xr = gkyl_nmat_get(x, count+1); 
  struct gkyl_mat x_uy_xl = gkyl_nmat_get(x, count+2); 
  struct gkyl_mat x_uy_xr = gkyl_nmat_get(x, count+3); 
  struct gkyl_mat x_uz_xl = gkyl_nmat_get(x, count+4); 
  struct gkyl_mat x_uz_xr = gkyl_nmat_get(x, count+5); 
  struct gkyl_mat x_Txx_xl = gkyl_nmat_get(x, count+6); 
  struct gkyl_mat x_Txx_xr = gkyl_nmat_get(x, count+7); 
 
  double *ux_xl = &prim_surf[0]; 
  double *ux_xr = &prim_surf[4]; 
  double *uy_xl = &prim_surf[8]; 
  double *uy_xr = &prim_surf[12]; 
  double *uz_xl = &prim_surf[16]; 
  double *uz_xr = &prim_surf[20]; 
  double *Txx_xl = &prim_surf[24]; 
  double *Txx_xr = &prim_surf[28]; 
 
  struct gkyl_mat x_ux_yl = gkyl_nmat_get(x, count+8); 
  struct gkyl_mat x_ux_yr = gkyl_nmat_get(x, count+9); 
  struct gkyl_mat x_uy_yl = gkyl_nmat_get(x, count+10); 
  struct gkyl_mat x_uy_yr = gkyl_nmat_get(x, count+11); 
  struct gkyl_mat x_uz_yl = gkyl_nmat_get(x, count+12); 
  struct gkyl_mat x_uz_yr = gkyl_nmat_get(x, count+13); 
  struct gkyl_mat x_Tyy_yl = gkyl_nmat_get(x, count+14); 
  struct gkyl_mat x_Tyy_yr = gkyl_nmat_get(x, count+15); 
 
  double *ux_yl = &prim_surf[32]; 
  double *ux_yr = &prim_surf[36]; 
  double *uy_yl = &prim_surf[40]; 
  double *uy_yr = &prim_surf[44]; 
  double *uz_yl = &prim_surf[48]; 
  double *uz_yr = &prim_surf[52]; 
  double *Tyy_yl = &prim_surf[56]; 
  double *Tyy_yr = &prim_surf[60]; 
 
  struct gkyl_mat x_ux_zl = gkyl_nmat_get(x, count+16); 
  struct gkyl_mat x_ux_zr = gkyl_nmat_get(x, count+17); 
  struct gkyl_mat x_uy_zl = gkyl_nmat_get(x, count+18); 
  struct gkyl_mat x_uy_zr = gkyl_nmat_get(x, count+19); 
  struct gkyl_mat x_uz_zl = gkyl_nmat_get(x, count+20); 
  struct gkyl_mat x_uz_zr = gkyl_nmat_get(x, count+21); 
  struct gkyl_mat x_Tzz_zl = gkyl_nmat_get(x, count+22); 
  struct gkyl_mat x_Tzz_zr = gkyl_nmat_get(x, count+23); 
 
  double *ux_zl = &prim_surf[64]; 
  double *ux_zr = &prim_surf[68]; 
  double *uy_zl = &prim_surf[72]; 
  double *uy_zr = &prim_surf[76]; 
  double *uz_zl = &prim_surf[80]; 
  double *uz_zr = &prim_surf[84]; 
  double *Tzz_zl = &prim_surf[88]; 
  double *Tzz_zr = &prim_surf[92]; 
 
  ux_xl[0] = gkyl_mat_get(&x_ux_xl,0,0); 
  ux_xr[0] = gkyl_mat_get(&x_ux_xr,0,0); 
  uy_xl[0] = gkyl_mat_get(&x_uy_xl,0,0); 
  uy_xr[0] = gkyl_mat_get(&x_uy_xr,0,0); 
  uz_xl[0] = gkyl_mat_get(&x_uz_xl,0,0); 
  uz_xr[0] = gkyl_mat_get(&x_uz_xr,0,0); 
  Txx_xl[0] = gkyl_mat_get(&x_Txx_xl,0,0); 
  Txx_xr[0] = gkyl_mat_get(&x_Txx_xr,0,0); 
 
  ux_yl[0] = gkyl_mat_get(&x_ux_yl,0,0); 
  ux_yr[0] = gkyl_mat_get(&x_ux_yr,0,0); 
  uy_yl[0] = gkyl_mat_get(&x_uy_yl,0,0); 
  uy_yr[0] = gkyl_mat_get(&x_uy_yr,0,0); 
  uz_yl[0] = gkyl_mat_get(&x_uz_yl,0,0); 
  uz_yr[0] = gkyl_mat_get(&x_uz_yr,0,0); 
  Tyy_yl[0] = gkyl_mat_get(&x_Tyy_yl,0,0); 
  Tyy_yr[0] = gkyl_mat_get(&x_Tyy_yr,0,0); 
 
  ux_zl[0] = gkyl_mat_get(&x_ux_zl,0,0); 
  ux_zr[0] = gkyl_mat_get(&x_ux_zr,0,0); 
  uy_zl[0] = gkyl_mat_get(&x_uy_zl,0,0); 
  uy_zr[0] = gkyl_mat_get(&x_uy_zr,0,0); 
  uz_zl[0] = gkyl_mat_get(&x_uz_zl,0,0); 
  uz_zr[0] = gkyl_mat_get(&x_uz_zr,0,0); 
  Tzz_zl[0] = gkyl_mat_get(&x_Tzz_zl,0,0); 
  Tzz_zr[0] = gkyl_mat_get(&x_Tzz_zr,0,0); 
 
  ux_xl[1] = gkyl_mat_get(&x_ux_xl,1,0); 
  ux_xr[1] = gkyl_mat_get(&x_ux_xr,1,0); 
  uy_xl[1] = gkyl_mat_get(&x_uy_xl,1,0); 
  uy_xr[1] = gkyl_mat_get(&x_uy_xr,1,0); 
  uz_xl[1] = gkyl_mat_get(&x_uz_xl,1,0); 
  uz_xr[1] = gkyl_mat_get(&x_uz_xr,1,0); 
  Txx_xl[1] = gkyl_mat_get(&x_Txx_xl,1,0); 
  Txx_xr[1] = gkyl_mat_get(&x_Txx_xr,1,0); 
 
  ux_yl[1] = gkyl_mat_get(&x_ux_yl,1,0); 
  ux_yr[1] = gkyl_mat_get(&x_ux_yr,1,0); 
  uy_yl[1] = gkyl_mat_get(&x_uy_yl,1,0); 
  uy_yr[1] = gkyl_mat_get(&x_uy_yr,1,0); 
  uz_yl[1] = gkyl_mat_get(&x_uz_yl,1,0); 
  uz_yr[1] = gkyl_mat_get(&x_uz_yr,1,0); 
  Tyy_yl[1] = gkyl_mat_get(&x_Tyy_yl,1,0); 
  Tyy_yr[1] = gkyl_mat_get(&x_Tyy_yr,1,0); 
 
  ux_zl[1] = gkyl_mat_get(&x_ux_zl,1,0); 
  ux_zr[1] = gkyl_mat_get(&x_ux_zr,1,0); 
  uy_zl[1] = gkyl_mat_get(&x_uy_zl,1,0); 
  uy_zr[1] = gkyl_mat_get(&x_uy_zr,1,0); 
  uz_zl[1] = gkyl_mat_get(&x_uz_zl,1,0); 
  uz_zr[1] = gkyl_mat_get(&x_uz_zr,1,0); 
  Tzz_zl[1] = gkyl_mat_get(&x_Tzz_zl,1,0); 
  Tzz_zr[1] = gkyl_mat_get(&x_Tzz_zr,1,0); 
 
  ux_xl[2] = gkyl_mat_get(&x_ux_xl,2,0); 
  ux_xr[2] = gkyl_mat_get(&x_ux_xr,2,0); 
  uy_xl[2] = gkyl_mat_get(&x_uy_xl,2,0); 
  uy_xr[2] = gkyl_mat_get(&x_uy_xr,2,0); 
  uz_xl[2] = gkyl_mat_get(&x_uz_xl,2,0); 
  uz_xr[2] = gkyl_mat_get(&x_uz_xr,2,0); 
  Txx_xl[2] = gkyl_mat_get(&x_Txx_xl,2,0); 
  Txx_xr[2] = gkyl_mat_get(&x_Txx_xr,2,0); 
 
  ux_yl[2] = gkyl_mat_get(&x_ux_yl,2,0); 
  ux_yr[2] = gkyl_mat_get(&x_ux_yr,2,0); 
  uy_yl[2] = gkyl_mat_get(&x_uy_yl,2,0); 
  uy_yr[2] = gkyl_mat_get(&x_uy_yr,2,0); 
  uz_yl[2] = gkyl_mat_get(&x_uz_yl,2,0); 
  uz_yr[2] = gkyl_mat_get(&x_uz_yr,2,0); 
  Tyy_yl[2] = gkyl_mat_get(&x_Tyy_yl,2,0); 
  Tyy_yr[2] = gkyl_mat_get(&x_Tyy_yr,2,0); 
 
  ux_zl[2] = gkyl_mat_get(&x_ux_zl,2,0); 
  ux_zr[2] = gkyl_mat_get(&x_ux_zr,2,0); 
  uy_zl[2] = gkyl_mat_get(&x_uy_zl,2,0); 
  uy_zr[2] = gkyl_mat_get(&x_uy_zr,2,0); 
  uz_zl[2] = gkyl_mat_get(&x_uz_zl,2,0); 
  uz_zr[2] = gkyl_mat_get(&x_uz_zr,2,0); 
  Tzz_zl[2] = gkyl_mat_get(&x_Tzz_zl,2,0); 
  Tzz_zr[2] = gkyl_mat_get(&x_Tzz_zr,2,0); 
 
  ux_xl[3] = gkyl_mat_get(&x_ux_xl,3,0); 
  ux_xr[3] = gkyl_mat_get(&x_ux_xr,3,0); 
  uy_xl[3] = gkyl_mat_get(&x_uy_xl,3,0); 
  uy_xr[3] = gkyl_mat_get(&x_uy_xr,3,0); 
  uz_xl[3] = gkyl_mat_get(&x_uz_xl,3,0); 
  uz_xr[3] = gkyl_mat_get(&x_uz_xr,3,0); 
  Txx_xl[3] = gkyl_mat_get(&x_Txx_xl,3,0); 
  Txx_xr[3] = gkyl_mat_get(&x_Txx_xr,3,0); 
 
  ux_yl[3] = gkyl_mat_get(&x_ux_yl,3,0); 
  ux_yr[3] = gkyl_mat_get(&x_ux_yr,3,0); 
  uy_yl[3] = gkyl_mat_get(&x_uy_yl,3,0); 
  uy_yr[3] = gkyl_mat_get(&x_uy_yr,3,0); 
  uz_yl[3] = gkyl_mat_get(&x_uz_yl,3,0); 
  uz_yr[3] = gkyl_mat_get(&x_uz_yr,3,0); 
  Tyy_yl[3] = gkyl_mat_get(&x_Tyy_yl,3,0); 
  Tyy_yr[3] = gkyl_mat_get(&x_Tyy_yr,3,0); 
 
  ux_zl[3] = gkyl_mat_get(&x_ux_zl,3,0); 
  ux_zr[3] = gkyl_mat_get(&x_ux_zr,3,0); 
  uy_zl[3] = gkyl_mat_get(&x_uy_zl,3,0); 
  uy_zr[3] = gkyl_mat_get(&x_uy_zr,3,0); 
  uz_zl[3] = gkyl_mat_get(&x_uz_zl,3,0); 
  uz_zr[3] = gkyl_mat_get(&x_uz_zr,3,0); 
  Tzz_zl[3] = gkyl_mat_get(&x_Tzz_zl,3,0); 
  Tzz_zr[3] = gkyl_mat_get(&x_Tzz_zr,3,0); 
 
} 
 
