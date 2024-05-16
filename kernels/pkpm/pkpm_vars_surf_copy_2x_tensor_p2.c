#include <gkyl_mat.h> 
#include <gkyl_euler_pkpm_kernels.h> 
GKYL_CU_DH void pkpm_vars_surf_copy_2x_tensor_p2(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT prim_surf) 
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
  double *ux_xr = &prim_surf[3]; 
  double *uy_xl = &prim_surf[6]; 
  double *uy_xr = &prim_surf[9]; 
  double *uz_xl = &prim_surf[12]; 
  double *uz_xr = &prim_surf[15]; 
  double *Txx_xl = &prim_surf[18]; 
  double *Txx_xr = &prim_surf[21]; 
 
  struct gkyl_mat x_ux_yl = gkyl_nmat_get(x, count+8); 
  struct gkyl_mat x_ux_yr = gkyl_nmat_get(x, count+9); 
  struct gkyl_mat x_uy_yl = gkyl_nmat_get(x, count+10); 
  struct gkyl_mat x_uy_yr = gkyl_nmat_get(x, count+11); 
  struct gkyl_mat x_uz_yl = gkyl_nmat_get(x, count+12); 
  struct gkyl_mat x_uz_yr = gkyl_nmat_get(x, count+13); 
  struct gkyl_mat x_Tyy_yl = gkyl_nmat_get(x, count+14); 
  struct gkyl_mat x_Tyy_yr = gkyl_nmat_get(x, count+15); 
 
  double *ux_yl = &prim_surf[24]; 
  double *ux_yr = &prim_surf[27]; 
  double *uy_yl = &prim_surf[30]; 
  double *uy_yr = &prim_surf[33]; 
  double *uz_yl = &prim_surf[36]; 
  double *uz_yr = &prim_surf[39]; 
  double *Tyy_yl = &prim_surf[42]; 
  double *Tyy_yr = &prim_surf[45]; 
 
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
 
} 
 
