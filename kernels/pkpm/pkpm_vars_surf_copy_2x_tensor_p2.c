#include <gkyl_mat.h> 
#include <gkyl_euler_pkpm_kernels.h> 
GKYL_CU_DH void pkpm_vars_surf_copy_2x_tensor_p2(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT prim_surf) 
{ 
  // count: integer to indicate which matrix being fetched. 
  // x:     Input solution vector. 
  // prim_surf: Primitive variables at surfaces.
 
  struct gkyl_mat x_Txx_xl = gkyl_nmat_get(x, count); 
  struct gkyl_mat x_Txx_xr = gkyl_nmat_get(x, count+1); 
 
  double *Txx_xl = &prim_surf[0]; 
  double *Txx_xr = &prim_surf[3]; 
 
  struct gkyl_mat x_Tyy_yl = gkyl_nmat_get(x, count+2); 
  struct gkyl_mat x_Tyy_yr = gkyl_nmat_get(x, count+3); 
 
  double *Tyy_yl = &prim_surf[6]; 
  double *Tyy_yr = &prim_surf[9]; 
 
  Txx_xl[0] = gkyl_mat_get(&x_Txx_xl,0,0); 
  Txx_xr[0] = gkyl_mat_get(&x_Txx_xr,0,0); 
 
  Tyy_yl[0] = gkyl_mat_get(&x_Tyy_yl,0,0); 
  Tyy_yr[0] = gkyl_mat_get(&x_Tyy_yr,0,0); 
 
  Txx_xl[1] = gkyl_mat_get(&x_Txx_xl,1,0); 
  Txx_xr[1] = gkyl_mat_get(&x_Txx_xr,1,0); 
 
  Tyy_yl[1] = gkyl_mat_get(&x_Tyy_yl,1,0); 
  Tyy_yr[1] = gkyl_mat_get(&x_Tyy_yr,1,0); 
 
  Txx_xl[2] = gkyl_mat_get(&x_Txx_xl,2,0); 
  Txx_xr[2] = gkyl_mat_get(&x_Txx_xr,2,0); 
 
  Tyy_yl[2] = gkyl_mat_get(&x_Tyy_yl,2,0); 
  Tyy_yr[2] = gkyl_mat_get(&x_Tyy_yr,2,0); 
 
} 
 
