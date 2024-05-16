#include <gkyl_mat.h> 
#include <gkyl_euler_pkpm_kernels.h> 
GKYL_CU_DH void pkpm_vars_copy_1x_ser_p2(int count, struct gkyl_nmat *x, 
    double* GKYL_RESTRICT prim) 
{ 
  // count:     integer to indicate which matrix being fetched. 
  // x:         Input solution vector. 
  // prim:      Output volume expansion of primitive variables: 
  //            [ux, uy, uz, 1/rho div(p_par b), T_perp/m, m/T_perp]. 
 
  struct gkyl_mat x_ux = gkyl_nmat_get(x, count); 
  struct gkyl_mat x_uy = gkyl_nmat_get(x, count+1); 
  struct gkyl_mat x_uz = gkyl_nmat_get(x, count+2); 
  struct gkyl_mat x_pkpm_div_ppar = gkyl_nmat_get(x, count+3); 
  struct gkyl_mat x_T_perp_over_m = gkyl_nmat_get(x, count+4); 
  struct gkyl_mat x_T_perp_over_m_inv = gkyl_nmat_get(x, count+5); 
  double *ux = &prim[0]; 
  double *uy = &prim[3]; 
  double *uz = &prim[6]; 
  double *p_force = &prim[9]; 
  double *T_perp_over_m = &prim[12]; 
  double *T_perp_over_m_inv = &prim[15]; 

  ux[0] = gkyl_mat_get(&x_ux,0,0); 
  uy[0] = gkyl_mat_get(&x_uy,0,0); 
  uz[0] = gkyl_mat_get(&x_uz,0,0); 
  p_force[0] = gkyl_mat_get(&x_pkpm_div_ppar,0,0); 
  T_perp_over_m[0] = gkyl_mat_get(&x_T_perp_over_m,0,0); 
  T_perp_over_m_inv[0] = gkyl_mat_get(&x_T_perp_over_m_inv,0,0); 
  ux[1] = gkyl_mat_get(&x_ux,1,0); 
  uy[1] = gkyl_mat_get(&x_uy,1,0); 
  uz[1] = gkyl_mat_get(&x_uz,1,0); 
  p_force[1] = gkyl_mat_get(&x_pkpm_div_ppar,1,0); 
  T_perp_over_m[1] = gkyl_mat_get(&x_T_perp_over_m,1,0); 
  T_perp_over_m_inv[1] = gkyl_mat_get(&x_T_perp_over_m_inv,1,0); 
  ux[2] = gkyl_mat_get(&x_ux,2,0); 
  uy[2] = gkyl_mat_get(&x_uy,2,0); 
  uz[2] = gkyl_mat_get(&x_uz,2,0); 
  p_force[2] = gkyl_mat_get(&x_pkpm_div_ppar,2,0); 
  T_perp_over_m[2] = gkyl_mat_get(&x_T_perp_over_m,2,0); 
  T_perp_over_m_inv[2] = gkyl_mat_get(&x_T_perp_over_m_inv,2,0); 

} 
 
