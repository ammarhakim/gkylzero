#include <gkyl_mat.h> 
#include <gkyl_euler_pkpm_kernels.h> 
GKYL_CU_DH void pkpm_vars_copy_1x_tensor_p2(int count, struct gkyl_nmat *x, 
    double* GKYL_RESTRICT pkpm_prim) 
{ 
  // count:     integer to indicate which matrix being fetched. 
  // x:         Input solution vector. 
  // pkpm_prim: Output volume expansion of pkpm pressure force variables: 
  //            [1/rho div(p_par b), T_perp/m, m/T_perp, 3*Txx/m, 3*Tyy/m, 3*Tzz/m]. 
 
  struct gkyl_mat x_pkpm_div_ppar_over_rho = gkyl_nmat_get(x, count); 
  struct gkyl_mat x_T_perp_over_m = gkyl_nmat_get(x, count+1); 
  struct gkyl_mat x_T_perp_over_m_inv = gkyl_nmat_get(x, count+2); 
  struct gkyl_mat x_Txx = gkyl_nmat_get(x, count+3); 
  struct gkyl_mat x_Tyy = gkyl_nmat_get(x, count+4); 
  struct gkyl_mat x_Tzz = gkyl_nmat_get(x, count+5); 
  double *pkpm_div_ppar_over_rho = &pkpm_prim[0]; 
  double *T_perp_over_m = &pkpm_prim[3]; 
  double *T_perp_over_m_inv = &pkpm_prim[6]; 
  double *Txx = &pkpm_prim[9]; 
  double *Tyy = &pkpm_prim[12]; 
  double *Tzz = &pkpm_prim[15]; 

  pkpm_div_ppar_over_rho[0] = gkyl_mat_get(&x_pkpm_div_ppar_over_rho,0,0); 
  T_perp_over_m[0] = gkyl_mat_get(&x_T_perp_over_m,0,0); 
  T_perp_over_m_inv[0] = gkyl_mat_get(&x_T_perp_over_m_inv,0,0); 
  Txx[0] = 3.0*gkyl_mat_get(&x_Txx,0,0); 
  Tyy[0] = 3.0*gkyl_mat_get(&x_Tyy,0,0); 
  Tzz[0] = 3.0*gkyl_mat_get(&x_Tzz,0,0); 
  pkpm_div_ppar_over_rho[1] = gkyl_mat_get(&x_pkpm_div_ppar_over_rho,1,0); 
  T_perp_over_m[1] = gkyl_mat_get(&x_T_perp_over_m,1,0); 
  T_perp_over_m_inv[1] = gkyl_mat_get(&x_T_perp_over_m_inv,1,0); 
  Txx[1] = 3.0*gkyl_mat_get(&x_Txx,1,0); 
  Tyy[1] = 3.0*gkyl_mat_get(&x_Tyy,1,0); 
  Tzz[1] = 3.0*gkyl_mat_get(&x_Tzz,1,0); 
  pkpm_div_ppar_over_rho[2] = gkyl_mat_get(&x_pkpm_div_ppar_over_rho,2,0); 
  T_perp_over_m[2] = gkyl_mat_get(&x_T_perp_over_m,2,0); 
  T_perp_over_m_inv[2] = gkyl_mat_get(&x_T_perp_over_m_inv,2,0); 
  Txx[2] = 3.0*gkyl_mat_get(&x_Txx,2,0); 
  Tyy[2] = 3.0*gkyl_mat_get(&x_Tyy,2,0); 
  Tzz[2] = 3.0*gkyl_mat_get(&x_Tzz,2,0); 

} 
 
