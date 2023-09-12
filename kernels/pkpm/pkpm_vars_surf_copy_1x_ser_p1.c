#include <gkyl_mat.h> 
#include <gkyl_euler_pkpm_kernels.h> 
GKYL_CU_DH void pkpm_vars_surf_copy_1x_ser_p1(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT prim_surf) 
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
  double *ux_xr = &prim_surf[1]; 
  double *uy_xl = &prim_surf[2]; 
  double *uy_xr = &prim_surf[3]; 
  double *uz_xl = &prim_surf[4]; 
  double *uz_xr = &prim_surf[5]; 
  double *Txx_xl = &prim_surf[6]; 
  double *Txx_xr = &prim_surf[7]; 
 
  ux_xl[0] = gkyl_mat_get(&x_ux_xl,0,0); 
  ux_xr[0] = gkyl_mat_get(&x_ux_xr,0,0); 
  uy_xl[0] = gkyl_mat_get(&x_uy_xl,0,0); 
  uy_xr[0] = gkyl_mat_get(&x_uy_xr,0,0); 
  uz_xl[0] = gkyl_mat_get(&x_uz_xl,0,0); 
  uz_xr[0] = gkyl_mat_get(&x_uz_xr,0,0); 
  Txx_xl[0] = gkyl_mat_get(&x_Txx_xl,0,0); 
  Txx_xr[0] = gkyl_mat_get(&x_Txx_xr,0,0); 
 
} 
 
