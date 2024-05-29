#include <gkyl_mat.h> 
#include <gkyl_euler_pkpm_kernels.h> 
GKYL_CU_DH void pkpm_vars_surf_copy_3x_tensor_p2(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT prim_surf) 
{ 
  // count: integer to indicate which matrix being fetched. 
  // x:     Input solution vector. 
  // prim_surf: Primitive variables at surfaces.
 
  struct gkyl_mat x_Txx_xl = gkyl_nmat_get(x, count); 
  struct gkyl_mat x_Txx_xr = gkyl_nmat_get(x, count+1); 
 
  double *Txx_xl = &prim_surf[0]; 
  double *Txx_xr = &prim_surf[9]; 
 
  struct gkyl_mat x_Tyy_yl = gkyl_nmat_get(x, count+2); 
  struct gkyl_mat x_Tyy_yr = gkyl_nmat_get(x, count+3); 
 
  double *Tyy_yl = &prim_surf[18]; 
  double *Tyy_yr = &prim_surf[27]; 
 
  struct gkyl_mat x_Tzz_zl = gkyl_nmat_get(x, count+4); 
  struct gkyl_mat x_Tzz_zr = gkyl_nmat_get(x, count+5); 
 
  double *Tzz_zl = &prim_surf[36]; 
  double *Tzz_zr = &prim_surf[45]; 
 
  Txx_xl[0] = gkyl_mat_get(&x_Txx_xl,0,0); 
  Txx_xr[0] = gkyl_mat_get(&x_Txx_xr,0,0); 
 
  Tyy_yl[0] = gkyl_mat_get(&x_Tyy_yl,0,0); 
  Tyy_yr[0] = gkyl_mat_get(&x_Tyy_yr,0,0); 
 
  Tzz_zl[0] = gkyl_mat_get(&x_Tzz_zl,0,0); 
  Tzz_zr[0] = gkyl_mat_get(&x_Tzz_zr,0,0); 
 
  Txx_xl[1] = gkyl_mat_get(&x_Txx_xl,1,0); 
  Txx_xr[1] = gkyl_mat_get(&x_Txx_xr,1,0); 
 
  Tyy_yl[1] = gkyl_mat_get(&x_Tyy_yl,1,0); 
  Tyy_yr[1] = gkyl_mat_get(&x_Tyy_yr,1,0); 
 
  Tzz_zl[1] = gkyl_mat_get(&x_Tzz_zl,1,0); 
  Tzz_zr[1] = gkyl_mat_get(&x_Tzz_zr,1,0); 
 
  Txx_xl[2] = gkyl_mat_get(&x_Txx_xl,2,0); 
  Txx_xr[2] = gkyl_mat_get(&x_Txx_xr,2,0); 
 
  Tyy_yl[2] = gkyl_mat_get(&x_Tyy_yl,2,0); 
  Tyy_yr[2] = gkyl_mat_get(&x_Tyy_yr,2,0); 
 
  Tzz_zl[2] = gkyl_mat_get(&x_Tzz_zl,2,0); 
  Tzz_zr[2] = gkyl_mat_get(&x_Tzz_zr,2,0); 
 
  Txx_xl[3] = gkyl_mat_get(&x_Txx_xl,3,0); 
  Txx_xr[3] = gkyl_mat_get(&x_Txx_xr,3,0); 
 
  Tyy_yl[3] = gkyl_mat_get(&x_Tyy_yl,3,0); 
  Tyy_yr[3] = gkyl_mat_get(&x_Tyy_yr,3,0); 
 
  Tzz_zl[3] = gkyl_mat_get(&x_Tzz_zl,3,0); 
  Tzz_zr[3] = gkyl_mat_get(&x_Tzz_zr,3,0); 
 
  Txx_xl[4] = gkyl_mat_get(&x_Txx_xl,4,0); 
  Txx_xr[4] = gkyl_mat_get(&x_Txx_xr,4,0); 
 
  Tyy_yl[4] = gkyl_mat_get(&x_Tyy_yl,4,0); 
  Tyy_yr[4] = gkyl_mat_get(&x_Tyy_yr,4,0); 
 
  Tzz_zl[4] = gkyl_mat_get(&x_Tzz_zl,4,0); 
  Tzz_zr[4] = gkyl_mat_get(&x_Tzz_zr,4,0); 
 
  Txx_xl[5] = gkyl_mat_get(&x_Txx_xl,5,0); 
  Txx_xr[5] = gkyl_mat_get(&x_Txx_xr,5,0); 
 
  Tyy_yl[5] = gkyl_mat_get(&x_Tyy_yl,5,0); 
  Tyy_yr[5] = gkyl_mat_get(&x_Tyy_yr,5,0); 
 
  Tzz_zl[5] = gkyl_mat_get(&x_Tzz_zl,5,0); 
  Tzz_zr[5] = gkyl_mat_get(&x_Tzz_zr,5,0); 
 
  Txx_xl[6] = gkyl_mat_get(&x_Txx_xl,6,0); 
  Txx_xr[6] = gkyl_mat_get(&x_Txx_xr,6,0); 
 
  Tyy_yl[6] = gkyl_mat_get(&x_Tyy_yl,6,0); 
  Tyy_yr[6] = gkyl_mat_get(&x_Tyy_yr,6,0); 
 
  Tzz_zl[6] = gkyl_mat_get(&x_Tzz_zl,6,0); 
  Tzz_zr[6] = gkyl_mat_get(&x_Tzz_zr,6,0); 
 
  Txx_xl[7] = gkyl_mat_get(&x_Txx_xl,7,0); 
  Txx_xr[7] = gkyl_mat_get(&x_Txx_xr,7,0); 
 
  Tyy_yl[7] = gkyl_mat_get(&x_Tyy_yl,7,0); 
  Tyy_yr[7] = gkyl_mat_get(&x_Tyy_yr,7,0); 
 
  Tzz_zl[7] = gkyl_mat_get(&x_Tzz_zl,7,0); 
  Tzz_zr[7] = gkyl_mat_get(&x_Tzz_zr,7,0); 
 
  Txx_xl[8] = gkyl_mat_get(&x_Txx_xl,8,0); 
  Txx_xr[8] = gkyl_mat_get(&x_Txx_xr,8,0); 
 
  Tyy_yl[8] = gkyl_mat_get(&x_Tyy_yl,8,0); 
  Tyy_yr[8] = gkyl_mat_get(&x_Tyy_yr,8,0); 
 
  Tzz_zl[8] = gkyl_mat_get(&x_Tzz_zl,8,0); 
  Tzz_zr[8] = gkyl_mat_get(&x_Tzz_zr,8,0); 
 
} 
 
