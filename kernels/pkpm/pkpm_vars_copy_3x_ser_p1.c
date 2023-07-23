#include <gkyl_mat.h> 
#include <gkyl_euler_pkpm_kernels.h> 
GKYL_CU_DH void pkpm_vars_copy_3x_ser_p1(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT prim) 
{ 
  // count: integer to indicate which matrix being fetched. 
  // x:     Input solution vector. 
  // prim:  [ux, uy, uz, 3*Txx/m, 3*Tyy/m, 3*Tzz/m, 1/rho div(p_par b), T_perp/m, m/T_perp].
 
  struct gkyl_mat x_ux = gkyl_nmat_get(x, count); 
  struct gkyl_mat x_uy = gkyl_nmat_get(x, count+1); 
  struct gkyl_mat x_uz = gkyl_nmat_get(x, count+2); 
  struct gkyl_mat x_Txx = gkyl_nmat_get(x, count+3); 
  struct gkyl_mat x_Tyy = gkyl_nmat_get(x, count+4); 
  struct gkyl_mat x_Tzz = gkyl_nmat_get(x, count+5); 
  struct gkyl_mat x_pkpm_div_ppar = gkyl_nmat_get(x, count+6); 
  struct gkyl_mat x_T_perp_over_m = gkyl_nmat_get(x, count+7); 
  struct gkyl_mat x_T_perp_over_m_inv = gkyl_nmat_get(x, count+8); 
  double *ux = &prim[0]; 
  double *uy = &prim[8]; 
  double *uz = &prim[16]; 
  double *Txx = &prim[24]; 
  double *Tyy = &prim[32]; 
  double *Tzz = &prim[40]; 
  double *p_force = &prim[48]; 
  double *T_perp_over_m = &prim[56]; 
  double *T_perp_over_m_inv = &prim[64]; 

  ux[0] = gkyl_mat_get(&x_ux,0,0); 
  uy[0] = gkyl_mat_get(&x_uy,0,0); 
  uz[0] = gkyl_mat_get(&x_uz,0,0); 
  Txx[0] = 3.0*gkyl_mat_get(&x_Txx,0,0); 
  Tyy[0] = 3.0*gkyl_mat_get(&x_Tyy,0,0); 
  Tzz[0] = 3.0*gkyl_mat_get(&x_Tzz,0,0); 
  p_force[0] = gkyl_mat_get(&x_pkpm_div_ppar,0,0); 
  T_perp_over_m[0] = gkyl_mat_get(&x_T_perp_over_m,0,0); 
  T_perp_over_m_inv[0] = gkyl_mat_get(&x_T_perp_over_m_inv,0,0); 
  ux[1] = gkyl_mat_get(&x_ux,1,0); 
  uy[1] = gkyl_mat_get(&x_uy,1,0); 
  uz[1] = gkyl_mat_get(&x_uz,1,0); 
  Txx[1] = 3.0*gkyl_mat_get(&x_Txx,1,0); 
  Tyy[1] = 3.0*gkyl_mat_get(&x_Tyy,1,0); 
  Tzz[1] = 3.0*gkyl_mat_get(&x_Tzz,1,0); 
  p_force[1] = gkyl_mat_get(&x_pkpm_div_ppar,1,0); 
  T_perp_over_m[1] = gkyl_mat_get(&x_T_perp_over_m,1,0); 
  T_perp_over_m_inv[1] = gkyl_mat_get(&x_T_perp_over_m_inv,1,0); 
  ux[2] = gkyl_mat_get(&x_ux,2,0); 
  uy[2] = gkyl_mat_get(&x_uy,2,0); 
  uz[2] = gkyl_mat_get(&x_uz,2,0); 
  Txx[2] = 3.0*gkyl_mat_get(&x_Txx,2,0); 
  Tyy[2] = 3.0*gkyl_mat_get(&x_Tyy,2,0); 
  Tzz[2] = 3.0*gkyl_mat_get(&x_Tzz,2,0); 
  p_force[2] = gkyl_mat_get(&x_pkpm_div_ppar,2,0); 
  T_perp_over_m[2] = gkyl_mat_get(&x_T_perp_over_m,2,0); 
  T_perp_over_m_inv[2] = gkyl_mat_get(&x_T_perp_over_m_inv,2,0); 
  ux[3] = gkyl_mat_get(&x_ux,3,0); 
  uy[3] = gkyl_mat_get(&x_uy,3,0); 
  uz[3] = gkyl_mat_get(&x_uz,3,0); 
  Txx[3] = 3.0*gkyl_mat_get(&x_Txx,3,0); 
  Tyy[3] = 3.0*gkyl_mat_get(&x_Tyy,3,0); 
  Tzz[3] = 3.0*gkyl_mat_get(&x_Tzz,3,0); 
  p_force[3] = gkyl_mat_get(&x_pkpm_div_ppar,3,0); 
  T_perp_over_m[3] = gkyl_mat_get(&x_T_perp_over_m,3,0); 
  T_perp_over_m_inv[3] = gkyl_mat_get(&x_T_perp_over_m_inv,3,0); 
  ux[4] = gkyl_mat_get(&x_ux,4,0); 
  uy[4] = gkyl_mat_get(&x_uy,4,0); 
  uz[4] = gkyl_mat_get(&x_uz,4,0); 
  Txx[4] = 3.0*gkyl_mat_get(&x_Txx,4,0); 
  Tyy[4] = 3.0*gkyl_mat_get(&x_Tyy,4,0); 
  Tzz[4] = 3.0*gkyl_mat_get(&x_Tzz,4,0); 
  p_force[4] = gkyl_mat_get(&x_pkpm_div_ppar,4,0); 
  T_perp_over_m[4] = gkyl_mat_get(&x_T_perp_over_m,4,0); 
  T_perp_over_m_inv[4] = gkyl_mat_get(&x_T_perp_over_m_inv,4,0); 
  ux[5] = gkyl_mat_get(&x_ux,5,0); 
  uy[5] = gkyl_mat_get(&x_uy,5,0); 
  uz[5] = gkyl_mat_get(&x_uz,5,0); 
  Txx[5] = 3.0*gkyl_mat_get(&x_Txx,5,0); 
  Tyy[5] = 3.0*gkyl_mat_get(&x_Tyy,5,0); 
  Tzz[5] = 3.0*gkyl_mat_get(&x_Tzz,5,0); 
  p_force[5] = gkyl_mat_get(&x_pkpm_div_ppar,5,0); 
  T_perp_over_m[5] = gkyl_mat_get(&x_T_perp_over_m,5,0); 
  T_perp_over_m_inv[5] = gkyl_mat_get(&x_T_perp_over_m_inv,5,0); 
  ux[6] = gkyl_mat_get(&x_ux,6,0); 
  uy[6] = gkyl_mat_get(&x_uy,6,0); 
  uz[6] = gkyl_mat_get(&x_uz,6,0); 
  Txx[6] = 3.0*gkyl_mat_get(&x_Txx,6,0); 
  Tyy[6] = 3.0*gkyl_mat_get(&x_Tyy,6,0); 
  Tzz[6] = 3.0*gkyl_mat_get(&x_Tzz,6,0); 
  p_force[6] = gkyl_mat_get(&x_pkpm_div_ppar,6,0); 
  T_perp_over_m[6] = gkyl_mat_get(&x_T_perp_over_m,6,0); 
  T_perp_over_m_inv[6] = gkyl_mat_get(&x_T_perp_over_m_inv,6,0); 
  ux[7] = gkyl_mat_get(&x_ux,7,0); 
  uy[7] = gkyl_mat_get(&x_uy,7,0); 
  uz[7] = gkyl_mat_get(&x_uz,7,0); 
  Txx[7] = 3.0*gkyl_mat_get(&x_Txx,7,0); 
  Tyy[7] = 3.0*gkyl_mat_get(&x_Tyy,7,0); 
  Tzz[7] = 3.0*gkyl_mat_get(&x_Tzz,7,0); 
  p_force[7] = gkyl_mat_get(&x_pkpm_div_ppar,7,0); 
  T_perp_over_m[7] = gkyl_mat_get(&x_T_perp_over_m,7,0); 
  T_perp_over_m_inv[7] = gkyl_mat_get(&x_T_perp_over_m_inv,7,0); 
} 
 
