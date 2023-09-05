#include <gkyl_mat.h> 
#include <gkyl_euler_pkpm_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
GKYL_CU_DH void pkpm_vars_surf_set_1x_ser_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  const double *p_ij_surf, const int *cell_avg_prim) 
{ 
  // count:            integer to indicate which matrix being fetched. 
  // A:                preallocated LHS matrix. 
  // rhs:              preallocated RHS vector. 
  // vlasov_pkpm_moms: Input [rho, p_parallel, p_perp], Moments computed from kinetic equation in pkpm model.
  // euler_pkpm:       Input [rho ux, rho uy, rho uz], Fluid state vector.
  // p_ij_surf:        Input surface expansion of p_ij = (p_par - p_perp) b_i b_j + p_perp g_ij.
  //                   [Pxx_xl, Pxx_xr, Pxy_xl, Pxy_xr, Pxz_xl, Pxz_xr, 
  //                    Pxy_yl, Pxy_yr, Pyy_yl, Pyy_yr, Pyz_yl, Pyz_yr, 
  //                    Pxz_zl, Pxz_zr, Pyz_zl, Pyz_zr, Pzz_zl, Pzz_zr] 
  // cell_avg_prim:    Boolean array to determine if we only use cell averages when computing surface expansions.

  // For poly_order = 1, we can analytically invert the matrix and just store the solution 
  struct gkyl_mat rhs_ux_l = gkyl_nmat_get(rhs, count); 
  struct gkyl_mat rhs_ux_r = gkyl_nmat_get(rhs, count+1); 
  struct gkyl_mat rhs_uy_l = gkyl_nmat_get(rhs, count+2); 
  struct gkyl_mat rhs_uy_r = gkyl_nmat_get(rhs, count+3); 
  struct gkyl_mat rhs_uz_l = gkyl_nmat_get(rhs, count+4); 
  struct gkyl_mat rhs_uz_r = gkyl_nmat_get(rhs, count+5); 
  struct gkyl_mat rhs_Txx_l = gkyl_nmat_get(rhs, count+6); 
  struct gkyl_mat rhs_Txx_r = gkyl_nmat_get(rhs, count+7); 
  // Clear rhs for each component of primitive variables being solved for 
  gkyl_mat_clear(&rhs_ux_l, 0.0); 
  gkyl_mat_clear(&rhs_ux_r, 0.0); 
  gkyl_mat_clear(&rhs_uy_l, 0.0); 
  gkyl_mat_clear(&rhs_uy_r, 0.0); 
  gkyl_mat_clear(&rhs_uz_l, 0.0); 
  gkyl_mat_clear(&rhs_uz_r, 0.0); 
  gkyl_mat_clear(&rhs_Txx_l, 0.0); 
  gkyl_mat_clear(&rhs_Txx_r, 0.0); 
  const double *rhoux = &euler_pkpm[0]; 
  const double *rhouy = &euler_pkpm[3]; 
  const double *rhouz = &euler_pkpm[6]; 
  const double *rho = &vlasov_pkpm_moms[0]; 
  const double *Pxx_xl = &p_ij_surf[0]; 
  const double *Pxx_xr = &p_ij_surf[1]; 
  double ux_l = 0.0; 
  double ux_r = 0.0; 
  double uy_l = 0.0; 
  double uy_r = 0.0; 
  double uz_l = 0.0; 
  double uz_r = 0.0; 
  double Txx_l = 0.0; 
  double Txx_r = 0.0; 
  if (cell_avg_prim[0]) { 
  // If rho or p_perp < 0 at control points, only use cell average. 
  ux_l = rhoux[0]/rho[0]; 
  ux_r = rhoux[0]/rho[0]; 
  uy_l = rhouy[0]/rho[0]; 
  uy_r = rhouy[0]/rho[0]; 
  uz_l = rhouz[0]/rho[0]; 
  uz_r = rhouz[0]/rho[0]; 
  Txx_l = Pxx_xl[0]/rho[0]; 
  Txx_r = Pxx_xr[0]/rho[0]; 
  } else { 
  double rhoux_l = 1.58113883008419*rhoux[2]-1.224744871391589*rhoux[1]+0.7071067811865475*rhoux[0]; 
  double rhoux_r = 1.58113883008419*rhoux[2]+1.224744871391589*rhoux[1]+0.7071067811865475*rhoux[0]; 
  double rhouy_l = 1.58113883008419*rhouy[2]-1.224744871391589*rhouy[1]+0.7071067811865475*rhouy[0]; 
  double rhouy_r = 1.58113883008419*rhouy[2]+1.224744871391589*rhouy[1]+0.7071067811865475*rhouy[0]; 
  double rhouz_l = 1.58113883008419*rhouz[2]-1.224744871391589*rhouz[1]+0.7071067811865475*rhouz[0]; 
  double rhouz_r = 1.58113883008419*rhouz[2]+1.224744871391589*rhouz[1]+0.7071067811865475*rhouz[0]; 
  double rho_l = 1.58113883008419*rho[2]-1.224744871391589*rho[1]+0.7071067811865475*rho[0]; 
  double rho_r = 1.58113883008419*rho[2]+1.224744871391589*rho[1]+0.7071067811865475*rho[0]; 
  double Pxx_l = Pxx_xl[0]; 
  double Pxx_r = Pxx_xr[0]; 
  ux_l = rhoux_l/rho_l; 
  ux_r = rhoux_r/rho_r; 
  uy_l = rhouy_l/rho_l; 
  uy_r = rhouy_r/rho_r; 
  uz_l = rhouz_l/rho_l; 
  uz_r = rhouz_r/rho_r; 
  Txx_l = Pxx_xl[0]/rho_l; 
  Txx_r = Pxx_xr[0]/rho_r; 
  } 
 
  gkyl_mat_set(&rhs_ux_l,0,0,ux_l); 
  gkyl_mat_set(&rhs_ux_r,0,0,ux_r); 
  gkyl_mat_set(&rhs_uy_l,0,0,uy_l); 
  gkyl_mat_set(&rhs_uy_r,0,0,uy_r); 
  gkyl_mat_set(&rhs_uz_l,0,0,uz_l); 
  gkyl_mat_set(&rhs_uz_r,0,0,uz_r); 
  gkyl_mat_set(&rhs_Txx_l,0,0,3.0*Txx_l); 
  gkyl_mat_set(&rhs_Txx_r,0,0,3.0*Txx_r); 
} 
