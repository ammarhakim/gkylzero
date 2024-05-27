#include <gkyl_mat.h> 
#include <gkyl_euler_pkpm_kernels.h> 
GKYL_CU_DH void pkpm_vars_copy_2x_tensor_p2(int count, struct gkyl_nmat *x, 
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
  double *T_perp_over_m = &pkpm_prim[9]; 
  double *T_perp_over_m_inv = &pkpm_prim[18]; 
  double *Txx = &pkpm_prim[27]; 
  double *Tyy = &pkpm_prim[36]; 
  double *Tzz = &pkpm_prim[45]; 

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
  pkpm_div_ppar_over_rho[3] = gkyl_mat_get(&x_pkpm_div_ppar_over_rho,3,0); 
  T_perp_over_m[3] = gkyl_mat_get(&x_T_perp_over_m,3,0); 
  T_perp_over_m_inv[3] = gkyl_mat_get(&x_T_perp_over_m_inv,3,0); 
  Txx[3] = 3.0*gkyl_mat_get(&x_Txx,3,0); 
  Tyy[3] = 3.0*gkyl_mat_get(&x_Tyy,3,0); 
  Tzz[3] = 3.0*gkyl_mat_get(&x_Tzz,3,0); 
  pkpm_div_ppar_over_rho[4] = gkyl_mat_get(&x_pkpm_div_ppar_over_rho,4,0); 
  T_perp_over_m[4] = gkyl_mat_get(&x_T_perp_over_m,4,0); 
  T_perp_over_m_inv[4] = gkyl_mat_get(&x_T_perp_over_m_inv,4,0); 
  Txx[4] = 3.0*gkyl_mat_get(&x_Txx,4,0); 
  Tyy[4] = 3.0*gkyl_mat_get(&x_Tyy,4,0); 
  Tzz[4] = 3.0*gkyl_mat_get(&x_Tzz,4,0); 
  pkpm_div_ppar_over_rho[5] = gkyl_mat_get(&x_pkpm_div_ppar_over_rho,5,0); 
  T_perp_over_m[5] = gkyl_mat_get(&x_T_perp_over_m,5,0); 
  T_perp_over_m_inv[5] = gkyl_mat_get(&x_T_perp_over_m_inv,5,0); 
  Txx[5] = 3.0*gkyl_mat_get(&x_Txx,5,0); 
  Tyy[5] = 3.0*gkyl_mat_get(&x_Tyy,5,0); 
  Tzz[5] = 3.0*gkyl_mat_get(&x_Tzz,5,0); 
  pkpm_div_ppar_over_rho[6] = gkyl_mat_get(&x_pkpm_div_ppar_over_rho,6,0); 
  T_perp_over_m[6] = gkyl_mat_get(&x_T_perp_over_m,6,0); 
  T_perp_over_m_inv[6] = gkyl_mat_get(&x_T_perp_over_m_inv,6,0); 
  Txx[6] = 3.0*gkyl_mat_get(&x_Txx,6,0); 
  Tyy[6] = 3.0*gkyl_mat_get(&x_Tyy,6,0); 
  Tzz[6] = 3.0*gkyl_mat_get(&x_Tzz,6,0); 
  pkpm_div_ppar_over_rho[7] = gkyl_mat_get(&x_pkpm_div_ppar_over_rho,7,0); 
  T_perp_over_m[7] = gkyl_mat_get(&x_T_perp_over_m,7,0); 
  T_perp_over_m_inv[7] = gkyl_mat_get(&x_T_perp_over_m_inv,7,0); 
  Txx[7] = 3.0*gkyl_mat_get(&x_Txx,7,0); 
  Tyy[7] = 3.0*gkyl_mat_get(&x_Tyy,7,0); 
  Tzz[7] = 3.0*gkyl_mat_get(&x_Tzz,7,0); 
  pkpm_div_ppar_over_rho[8] = gkyl_mat_get(&x_pkpm_div_ppar_over_rho,8,0); 
  T_perp_over_m[8] = gkyl_mat_get(&x_T_perp_over_m,8,0); 
  T_perp_over_m_inv[8] = gkyl_mat_get(&x_T_perp_over_m_inv,8,0); 
  Txx[8] = 3.0*gkyl_mat_get(&x_Txx,8,0); 
  Tyy[8] = 3.0*gkyl_mat_get(&x_Tyy,8,0); 
  Tzz[8] = 3.0*gkyl_mat_get(&x_Tzz,8,0); 

} 
 
