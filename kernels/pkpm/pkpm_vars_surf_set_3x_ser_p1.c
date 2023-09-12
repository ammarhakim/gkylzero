#include <gkyl_mat.h> 
#include <gkyl_euler_pkpm_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
#include <gkyl_basis_ser_2x_p1_inv.h> 
GKYL_CU_DH void pkpm_vars_surf_set_3x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  const double *p_ij_surf, const int *cell_avg_prim) 
{ 
  // count:            integer to indicate which matrix being fetched. 
  // A:                preallocated LHS matrix. 
  // rhs:              preallocated RHS vector. 
  // vlasov_pkpm_moms: Input [rho, p_parallel, p_perp], Moments computed from kinetic equation in pkpm model.
  // euler_pkpm:       Input [rho ux, rho uy, rho uz], Fluid input state vector.
  // p_ij_surf:        Input surface expansion of p_ij = (p_par - p_perp) b_i b_j + p_perp g_ij.
  //                   [Pxx_xl, Pxx_xr, Pxy_xl, Pxy_xr, Pxz_xl, Pxz_xr, 
  //                    Pxy_yl, Pxy_yr, Pyy_yl, Pyy_yr, Pyz_yl, Pyz_yr, 
  //                    Pxz_zl, Pxz_zr, Pyz_zl, Pyz_zr, Pzz_zl, Pzz_zr] 
  // cell_avg_prim:    Boolean array to determine if we only use cell averages when computing surface expansions.

  // For poly_order = 1, we can analytically invert the matrix and just store the solution 
  struct gkyl_mat rhs_ux_xl = gkyl_nmat_get(rhs, count); 
  struct gkyl_mat rhs_ux_xr = gkyl_nmat_get(rhs, count+1); 
  struct gkyl_mat rhs_uy_xl = gkyl_nmat_get(rhs, count+2); 
  struct gkyl_mat rhs_uy_xr = gkyl_nmat_get(rhs, count+3); 
  struct gkyl_mat rhs_uz_xl = gkyl_nmat_get(rhs, count+4); 
  struct gkyl_mat rhs_uz_xr = gkyl_nmat_get(rhs, count+5); 
  struct gkyl_mat rhs_Txx_xl = gkyl_nmat_get(rhs, count+6); 
  struct gkyl_mat rhs_Txx_xr = gkyl_nmat_get(rhs, count+7); 
  struct gkyl_mat rhs_ux_yl = gkyl_nmat_get(rhs, count+8); 
  struct gkyl_mat rhs_ux_yr = gkyl_nmat_get(rhs, count+9); 
  struct gkyl_mat rhs_uy_yl = gkyl_nmat_get(rhs, count+10); 
  struct gkyl_mat rhs_uy_yr = gkyl_nmat_get(rhs, count+11); 
  struct gkyl_mat rhs_uz_yl = gkyl_nmat_get(rhs, count+12); 
  struct gkyl_mat rhs_uz_yr = gkyl_nmat_get(rhs, count+13); 
  struct gkyl_mat rhs_Tyy_yl = gkyl_nmat_get(rhs, count+14); 
  struct gkyl_mat rhs_Tyy_yr = gkyl_nmat_get(rhs, count+15); 
  struct gkyl_mat rhs_ux_zl = gkyl_nmat_get(rhs, count+16); 
  struct gkyl_mat rhs_ux_zr = gkyl_nmat_get(rhs, count+17); 
  struct gkyl_mat rhs_uy_zl = gkyl_nmat_get(rhs, count+18); 
  struct gkyl_mat rhs_uy_zr = gkyl_nmat_get(rhs, count+19); 
  struct gkyl_mat rhs_uz_zl = gkyl_nmat_get(rhs, count+20); 
  struct gkyl_mat rhs_uz_zr = gkyl_nmat_get(rhs, count+21); 
  struct gkyl_mat rhs_Tzz_zl = gkyl_nmat_get(rhs, count+22); 
  struct gkyl_mat rhs_Tzz_zr = gkyl_nmat_get(rhs, count+23); 
  // Clear rhs for each component of primitive variables being solved for 
  gkyl_mat_clear(&rhs_ux_xl, 0.0); 
  gkyl_mat_clear(&rhs_ux_xr, 0.0); 
  gkyl_mat_clear(&rhs_uy_xl, 0.0); 
  gkyl_mat_clear(&rhs_uy_xr, 0.0); 
  gkyl_mat_clear(&rhs_uz_xl, 0.0); 
  gkyl_mat_clear(&rhs_uz_xr, 0.0); 
  gkyl_mat_clear(&rhs_Txx_xl, 0.0); 
  gkyl_mat_clear(&rhs_Txx_xr, 0.0); 
  gkyl_mat_clear(&rhs_ux_yl, 0.0); 
  gkyl_mat_clear(&rhs_ux_yr, 0.0); 
  gkyl_mat_clear(&rhs_uy_yl, 0.0); 
  gkyl_mat_clear(&rhs_uy_yr, 0.0); 
  gkyl_mat_clear(&rhs_uz_yl, 0.0); 
  gkyl_mat_clear(&rhs_uz_yr, 0.0); 
  gkyl_mat_clear(&rhs_Tyy_yl, 0.0); 
  gkyl_mat_clear(&rhs_Tyy_yr, 0.0); 
  gkyl_mat_clear(&rhs_ux_zl, 0.0); 
  gkyl_mat_clear(&rhs_ux_zr, 0.0); 
  gkyl_mat_clear(&rhs_uy_zl, 0.0); 
  gkyl_mat_clear(&rhs_uy_zr, 0.0); 
  gkyl_mat_clear(&rhs_uz_zl, 0.0); 
  gkyl_mat_clear(&rhs_uz_zr, 0.0); 
  gkyl_mat_clear(&rhs_Tzz_zl, 0.0); 
  gkyl_mat_clear(&rhs_Tzz_zr, 0.0); 
  const double *rhoux = &euler_pkpm[0]; 
  const double *rhouy = &euler_pkpm[8]; 
  const double *rhouz = &euler_pkpm[16]; 
  const double *rho = &vlasov_pkpm_moms[0]; 
  const double *Pxx_xl = &p_ij_surf[0]; 
  const double *Pxx_xr = &p_ij_surf[4]; 
  const double *Pyy_yl = &p_ij_surf[32]; 
  const double *Pyy_yr = &p_ij_surf[36]; 
  const double *Pzz_zl = &p_ij_surf[64]; 
  const double *Pzz_zr = &p_ij_surf[68]; 
  double ux_xl[4] = {0.0}; 
  double ux_xr[4] = {0.0}; 
  double uy_xl[4] = {0.0}; 
  double uy_xr[4] = {0.0}; 
  double uz_xl[4] = {0.0}; 
  double uz_xr[4] = {0.0}; 
  double Txx_xl[4] = {0.0}; 
  double Txx_xr[4] = {0.0}; 
  double ux_yl[4] = {0.0}; 
  double ux_yr[4] = {0.0}; 
  double uy_yl[4] = {0.0}; 
  double uy_yr[4] = {0.0}; 
  double uz_yl[4] = {0.0}; 
  double uz_yr[4] = {0.0}; 
  double Tyy_yl[4] = {0.0}; 
  double Tyy_yr[4] = {0.0}; 
  double ux_zl[4] = {0.0}; 
  double ux_zr[4] = {0.0}; 
  double uy_zl[4] = {0.0}; 
  double uy_zr[4] = {0.0}; 
  double uz_zl[4] = {0.0}; 
  double uz_zr[4] = {0.0}; 
  double Tzz_zl[4] = {0.0}; 
  double Tzz_zr[4] = {0.0}; 
  if (cell_avg_prim[0]) { 
  // If rho or p_perp < 0 at control points, only use cell average. 
  ux_xl[0] = rhoux[0]/rho[0]; 
  ux_xr[0] = rhoux[0]/rho[0]; 
  uy_xl[0] = rhouy[0]/rho[0]; 
  uy_xr[0] = rhouy[0]/rho[0]; 
  uz_xl[0] = rhouz[0]/rho[0]; 
  uz_xr[0] = rhouz[0]/rho[0]; 
  Txx_xl[0] = Pxx_xl[0]/rho[0]; 
  Txx_xr[0] = Pxx_xr[0]/rho[0]; 
  ux_yl[0] = rhoux[0]/rho[0]; 
  ux_yr[0] = rhoux[0]/rho[0]; 
  uy_yl[0] = rhouy[0]/rho[0]; 
  uy_yr[0] = rhouy[0]/rho[0]; 
  uz_yl[0] = rhouz[0]/rho[0]; 
  uz_yr[0] = rhouz[0]/rho[0]; 
  Tyy_yl[0] = Pyy_yl[0]/rho[0]; 
  Tyy_yr[0] = Pyy_yr[0]/rho[0]; 
  ux_zl[0] = rhoux[0]/rho[0]; 
  ux_zr[0] = rhoux[0]/rho[0]; 
  uy_zl[0] = rhouy[0]/rho[0]; 
  uy_zr[0] = rhouy[0]/rho[0]; 
  uz_zl[0] = rhouz[0]/rho[0]; 
  uz_zr[0] = rhouz[0]/rho[0]; 
  Tzz_zl[0] = Pzz_zl[0]/rho[0]; 
  Tzz_zr[0] = Pzz_zr[0]/rho[0]; 
  } else { 
  double rhoux_xl[4] = {0.0}; 
  double rhoux_xr[4] = {0.0}; 
  double rhouy_xl[4] = {0.0}; 
  double rhouy_xr[4] = {0.0}; 
  double rhouz_xl[4] = {0.0}; 
  double rhouz_xr[4] = {0.0}; 
  double rho_xl[4] = {0.0}; 
  double rho_xr[4] = {0.0}; 
  double rho_inv_xl[4] = {0.0}; 
  double rho_inv_xr[4] = {0.0}; 
 
  rhoux_xl[0] = 0.7071067811865475*rhoux[0]-1.224744871391589*rhoux[1]; 
  rhoux_xl[1] = 0.7071067811865475*rhoux[2]-1.224744871391589*rhoux[4]; 
  rhoux_xl[2] = 0.7071067811865475*rhoux[3]-1.224744871391589*rhoux[5]; 
  rhoux_xl[3] = 0.7071067811865475*rhoux[6]-1.224744871391589*rhoux[7]; 
  rhoux_xr[0] = 1.224744871391589*rhoux[1]+0.7071067811865475*rhoux[0]; 
  rhoux_xr[1] = 1.224744871391589*rhoux[4]+0.7071067811865475*rhoux[2]; 
  rhoux_xr[2] = 1.224744871391589*rhoux[5]+0.7071067811865475*rhoux[3]; 
  rhoux_xr[3] = 1.224744871391589*rhoux[7]+0.7071067811865475*rhoux[6]; 
  rhouy_xl[0] = 0.7071067811865475*rhouy[0]-1.224744871391589*rhouy[1]; 
  rhouy_xl[1] = 0.7071067811865475*rhouy[2]-1.224744871391589*rhouy[4]; 
  rhouy_xl[2] = 0.7071067811865475*rhouy[3]-1.224744871391589*rhouy[5]; 
  rhouy_xl[3] = 0.7071067811865475*rhouy[6]-1.224744871391589*rhouy[7]; 
  rhouy_xr[0] = 1.224744871391589*rhouy[1]+0.7071067811865475*rhouy[0]; 
  rhouy_xr[1] = 1.224744871391589*rhouy[4]+0.7071067811865475*rhouy[2]; 
  rhouy_xr[2] = 1.224744871391589*rhouy[5]+0.7071067811865475*rhouy[3]; 
  rhouy_xr[3] = 1.224744871391589*rhouy[7]+0.7071067811865475*rhouy[6]; 
  rhouz_xl[0] = 0.7071067811865475*rhouz[0]-1.224744871391589*rhouz[1]; 
  rhouz_xl[1] = 0.7071067811865475*rhouz[2]-1.224744871391589*rhouz[4]; 
  rhouz_xl[2] = 0.7071067811865475*rhouz[3]-1.224744871391589*rhouz[5]; 
  rhouz_xl[3] = 0.7071067811865475*rhouz[6]-1.224744871391589*rhouz[7]; 
  rhouz_xr[0] = 1.224744871391589*rhouz[1]+0.7071067811865475*rhouz[0]; 
  rhouz_xr[1] = 1.224744871391589*rhouz[4]+0.7071067811865475*rhouz[2]; 
  rhouz_xr[2] = 1.224744871391589*rhouz[5]+0.7071067811865475*rhouz[3]; 
  rhouz_xr[3] = 1.224744871391589*rhouz[7]+0.7071067811865475*rhouz[6]; 
  rho_xl[0] = 0.7071067811865475*rho[0]-1.224744871391589*rho[1]; 
  rho_xl[1] = 0.7071067811865475*rho[2]-1.224744871391589*rho[4]; 
  rho_xl[2] = 0.7071067811865475*rho[3]-1.224744871391589*rho[5]; 
  rho_xl[3] = 0.7071067811865475*rho[6]-1.224744871391589*rho[7]; 
  rho_xr[0] = 1.224744871391589*rho[1]+0.7071067811865475*rho[0]; 
  rho_xr[1] = 1.224744871391589*rho[4]+0.7071067811865475*rho[2]; 
  rho_xr[2] = 1.224744871391589*rho[5]+0.7071067811865475*rho[3]; 
  rho_xr[3] = 1.224744871391589*rho[7]+0.7071067811865475*rho[6]; 
  ser_2x_p1_inv(rho_xl, rho_inv_xl); 
  ser_2x_p1_inv(rho_xr, rho_inv_xr); 
  binop_mul_2d_ser_p1(rho_inv_xl, rhoux_xl, ux_xl); 
  binop_mul_2d_ser_p1(rho_inv_xr, rhoux_xr, ux_xr); 
  binop_mul_2d_ser_p1(rho_inv_xl, rhouy_xl, uy_xl); 
  binop_mul_2d_ser_p1(rho_inv_xr, rhouy_xr, uy_xr); 
  binop_mul_2d_ser_p1(rho_inv_xl, rhouz_xl, uz_xl); 
  binop_mul_2d_ser_p1(rho_inv_xr, rhouz_xr, uz_xr); 
  binop_mul_2d_ser_p1(rho_inv_xl, Pxx_xl, Txx_xl); 
  binop_mul_2d_ser_p1(rho_inv_xr, Pxx_xr, Txx_xr); 
 
  double rhoux_yl[4] = {0.0}; 
  double rhoux_yr[4] = {0.0}; 
  double rhouy_yl[4] = {0.0}; 
  double rhouy_yr[4] = {0.0}; 
  double rhouz_yl[4] = {0.0}; 
  double rhouz_yr[4] = {0.0}; 
  double rho_yl[4] = {0.0}; 
  double rho_yr[4] = {0.0}; 
  double rho_inv_yl[4] = {0.0}; 
  double rho_inv_yr[4] = {0.0}; 
 
  rhoux_yl[0] = 0.7071067811865475*rhoux[0]-1.224744871391589*rhoux[2]; 
  rhoux_yl[1] = 0.7071067811865475*rhoux[1]-1.224744871391589*rhoux[4]; 
  rhoux_yl[2] = 0.7071067811865475*rhoux[3]-1.224744871391589*rhoux[6]; 
  rhoux_yl[3] = 0.7071067811865475*rhoux[5]-1.224744871391589*rhoux[7]; 
  rhoux_yr[0] = 1.224744871391589*rhoux[2]+0.7071067811865475*rhoux[0]; 
  rhoux_yr[1] = 1.224744871391589*rhoux[4]+0.7071067811865475*rhoux[1]; 
  rhoux_yr[2] = 1.224744871391589*rhoux[6]+0.7071067811865475*rhoux[3]; 
  rhoux_yr[3] = 1.224744871391589*rhoux[7]+0.7071067811865475*rhoux[5]; 
  rhouy_yl[0] = 0.7071067811865475*rhouy[0]-1.224744871391589*rhouy[2]; 
  rhouy_yl[1] = 0.7071067811865475*rhouy[1]-1.224744871391589*rhouy[4]; 
  rhouy_yl[2] = 0.7071067811865475*rhouy[3]-1.224744871391589*rhouy[6]; 
  rhouy_yl[3] = 0.7071067811865475*rhouy[5]-1.224744871391589*rhouy[7]; 
  rhouy_yr[0] = 1.224744871391589*rhouy[2]+0.7071067811865475*rhouy[0]; 
  rhouy_yr[1] = 1.224744871391589*rhouy[4]+0.7071067811865475*rhouy[1]; 
  rhouy_yr[2] = 1.224744871391589*rhouy[6]+0.7071067811865475*rhouy[3]; 
  rhouy_yr[3] = 1.224744871391589*rhouy[7]+0.7071067811865475*rhouy[5]; 
  rhouz_yl[0] = 0.7071067811865475*rhouz[0]-1.224744871391589*rhouz[2]; 
  rhouz_yl[1] = 0.7071067811865475*rhouz[1]-1.224744871391589*rhouz[4]; 
  rhouz_yl[2] = 0.7071067811865475*rhouz[3]-1.224744871391589*rhouz[6]; 
  rhouz_yl[3] = 0.7071067811865475*rhouz[5]-1.224744871391589*rhouz[7]; 
  rhouz_yr[0] = 1.224744871391589*rhouz[2]+0.7071067811865475*rhouz[0]; 
  rhouz_yr[1] = 1.224744871391589*rhouz[4]+0.7071067811865475*rhouz[1]; 
  rhouz_yr[2] = 1.224744871391589*rhouz[6]+0.7071067811865475*rhouz[3]; 
  rhouz_yr[3] = 1.224744871391589*rhouz[7]+0.7071067811865475*rhouz[5]; 
  rho_yl[0] = 0.7071067811865475*rho[0]-1.224744871391589*rho[2]; 
  rho_yl[1] = 0.7071067811865475*rho[1]-1.224744871391589*rho[4]; 
  rho_yl[2] = 0.7071067811865475*rho[3]-1.224744871391589*rho[6]; 
  rho_yl[3] = 0.7071067811865475*rho[5]-1.224744871391589*rho[7]; 
  rho_yr[0] = 1.224744871391589*rho[2]+0.7071067811865475*rho[0]; 
  rho_yr[1] = 1.224744871391589*rho[4]+0.7071067811865475*rho[1]; 
  rho_yr[2] = 1.224744871391589*rho[6]+0.7071067811865475*rho[3]; 
  rho_yr[3] = 1.224744871391589*rho[7]+0.7071067811865475*rho[5]; 
  ser_2x_p1_inv(rho_yl, rho_inv_yl); 
  ser_2x_p1_inv(rho_yr, rho_inv_yr); 
  binop_mul_2d_ser_p1(rho_inv_yl, rhoux_yl, ux_yl); 
  binop_mul_2d_ser_p1(rho_inv_yr, rhoux_yr, ux_yr); 
  binop_mul_2d_ser_p1(rho_inv_yl, rhouy_yl, uy_yl); 
  binop_mul_2d_ser_p1(rho_inv_yr, rhouy_yr, uy_yr); 
  binop_mul_2d_ser_p1(rho_inv_yl, rhouz_yl, uz_yl); 
  binop_mul_2d_ser_p1(rho_inv_yr, rhouz_yr, uz_yr); 
  binop_mul_2d_ser_p1(rho_inv_yl, Pyy_yl, Tyy_yl); 
  binop_mul_2d_ser_p1(rho_inv_yr, Pyy_yr, Tyy_yr); 
 
  double rhoux_zl[4] = {0.0}; 
  double rhoux_zr[4] = {0.0}; 
  double rhouy_zl[4] = {0.0}; 
  double rhouy_zr[4] = {0.0}; 
  double rhouz_zl[4] = {0.0}; 
  double rhouz_zr[4] = {0.0}; 
  double rho_zl[4] = {0.0}; 
  double rho_zr[4] = {0.0}; 
  double rho_inv_zl[4] = {0.0}; 
  double rho_inv_zr[4] = {0.0}; 
 
  rhoux_zl[0] = 0.7071067811865475*rhoux[0]-1.224744871391589*rhoux[3]; 
  rhoux_zl[1] = 0.7071067811865475*rhoux[1]-1.224744871391589*rhoux[5]; 
  rhoux_zl[2] = 0.7071067811865475*rhoux[2]-1.224744871391589*rhoux[6]; 
  rhoux_zl[3] = 0.7071067811865475*rhoux[4]-1.224744871391589*rhoux[7]; 
  rhoux_zr[0] = 1.224744871391589*rhoux[3]+0.7071067811865475*rhoux[0]; 
  rhoux_zr[1] = 1.224744871391589*rhoux[5]+0.7071067811865475*rhoux[1]; 
  rhoux_zr[2] = 1.224744871391589*rhoux[6]+0.7071067811865475*rhoux[2]; 
  rhoux_zr[3] = 1.224744871391589*rhoux[7]+0.7071067811865475*rhoux[4]; 
  rhouy_zl[0] = 0.7071067811865475*rhouy[0]-1.224744871391589*rhouy[3]; 
  rhouy_zl[1] = 0.7071067811865475*rhouy[1]-1.224744871391589*rhouy[5]; 
  rhouy_zl[2] = 0.7071067811865475*rhouy[2]-1.224744871391589*rhouy[6]; 
  rhouy_zl[3] = 0.7071067811865475*rhouy[4]-1.224744871391589*rhouy[7]; 
  rhouy_zr[0] = 1.224744871391589*rhouy[3]+0.7071067811865475*rhouy[0]; 
  rhouy_zr[1] = 1.224744871391589*rhouy[5]+0.7071067811865475*rhouy[1]; 
  rhouy_zr[2] = 1.224744871391589*rhouy[6]+0.7071067811865475*rhouy[2]; 
  rhouy_zr[3] = 1.224744871391589*rhouy[7]+0.7071067811865475*rhouy[4]; 
  rhouz_zl[0] = 0.7071067811865475*rhouz[0]-1.224744871391589*rhouz[3]; 
  rhouz_zl[1] = 0.7071067811865475*rhouz[1]-1.224744871391589*rhouz[5]; 
  rhouz_zl[2] = 0.7071067811865475*rhouz[2]-1.224744871391589*rhouz[6]; 
  rhouz_zl[3] = 0.7071067811865475*rhouz[4]-1.224744871391589*rhouz[7]; 
  rhouz_zr[0] = 1.224744871391589*rhouz[3]+0.7071067811865475*rhouz[0]; 
  rhouz_zr[1] = 1.224744871391589*rhouz[5]+0.7071067811865475*rhouz[1]; 
  rhouz_zr[2] = 1.224744871391589*rhouz[6]+0.7071067811865475*rhouz[2]; 
  rhouz_zr[3] = 1.224744871391589*rhouz[7]+0.7071067811865475*rhouz[4]; 
  rho_zl[0] = 0.7071067811865475*rho[0]-1.224744871391589*rho[3]; 
  rho_zl[1] = 0.7071067811865475*rho[1]-1.224744871391589*rho[5]; 
  rho_zl[2] = 0.7071067811865475*rho[2]-1.224744871391589*rho[6]; 
  rho_zl[3] = 0.7071067811865475*rho[4]-1.224744871391589*rho[7]; 
  rho_zr[0] = 1.224744871391589*rho[3]+0.7071067811865475*rho[0]; 
  rho_zr[1] = 1.224744871391589*rho[5]+0.7071067811865475*rho[1]; 
  rho_zr[2] = 1.224744871391589*rho[6]+0.7071067811865475*rho[2]; 
  rho_zr[3] = 1.224744871391589*rho[7]+0.7071067811865475*rho[4]; 
  ser_2x_p1_inv(rho_zl, rho_inv_zl); 
  ser_2x_p1_inv(rho_zr, rho_inv_zr); 
  binop_mul_2d_ser_p1(rho_inv_zl, rhoux_zl, ux_zl); 
  binop_mul_2d_ser_p1(rho_inv_zr, rhoux_zr, ux_zr); 
  binop_mul_2d_ser_p1(rho_inv_zl, rhouy_zl, uy_zl); 
  binop_mul_2d_ser_p1(rho_inv_zr, rhouy_zr, uy_zr); 
  binop_mul_2d_ser_p1(rho_inv_zl, rhouz_zl, uz_zl); 
  binop_mul_2d_ser_p1(rho_inv_zr, rhouz_zr, uz_zr); 
  binop_mul_2d_ser_p1(rho_inv_zl, Pzz_zl, Tzz_zl); 
  binop_mul_2d_ser_p1(rho_inv_zr, Pzz_zr, Tzz_zr); 
 
  } 
 
  gkyl_mat_set(&rhs_ux_xl,0,0,ux_xl[0]); 
  gkyl_mat_set(&rhs_ux_xr,0,0,ux_xr[0]); 
  gkyl_mat_set(&rhs_uy_xl,0,0,uy_xl[0]); 
  gkyl_mat_set(&rhs_uy_xr,0,0,uy_xr[0]); 
  gkyl_mat_set(&rhs_uz_xl,0,0,uz_xl[0]); 
  gkyl_mat_set(&rhs_uz_xr,0,0,uz_xr[0]); 
  gkyl_mat_set(&rhs_Txx_xl,0,0,3.0*Txx_xl[0]); 
  gkyl_mat_set(&rhs_Txx_xr,0,0,3.0*Txx_xr[0]); 
 
  gkyl_mat_set(&rhs_ux_yl,0,0,ux_yl[0]); 
  gkyl_mat_set(&rhs_ux_yr,0,0,ux_yr[0]); 
  gkyl_mat_set(&rhs_uy_yl,0,0,uy_yl[0]); 
  gkyl_mat_set(&rhs_uy_yr,0,0,uy_yr[0]); 
  gkyl_mat_set(&rhs_uz_yl,0,0,uz_yl[0]); 
  gkyl_mat_set(&rhs_uz_yr,0,0,uz_yr[0]); 
  gkyl_mat_set(&rhs_Tyy_yl,0,0,3.0*Tyy_yl[0]); 
  gkyl_mat_set(&rhs_Tyy_yr,0,0,3.0*Tyy_yr[0]); 
 
  gkyl_mat_set(&rhs_ux_zl,0,0,ux_zl[0]); 
  gkyl_mat_set(&rhs_ux_zr,0,0,ux_zr[0]); 
  gkyl_mat_set(&rhs_uy_zl,0,0,uy_zl[0]); 
  gkyl_mat_set(&rhs_uy_zr,0,0,uy_zr[0]); 
  gkyl_mat_set(&rhs_uz_zl,0,0,uz_zl[0]); 
  gkyl_mat_set(&rhs_uz_zr,0,0,uz_zr[0]); 
  gkyl_mat_set(&rhs_Tzz_zl,0,0,3.0*Tzz_zl[0]); 
  gkyl_mat_set(&rhs_Tzz_zr,0,0,3.0*Tzz_zr[0]); 
 
  gkyl_mat_set(&rhs_ux_xl,1,0,ux_xl[1]); 
  gkyl_mat_set(&rhs_ux_xr,1,0,ux_xr[1]); 
  gkyl_mat_set(&rhs_uy_xl,1,0,uy_xl[1]); 
  gkyl_mat_set(&rhs_uy_xr,1,0,uy_xr[1]); 
  gkyl_mat_set(&rhs_uz_xl,1,0,uz_xl[1]); 
  gkyl_mat_set(&rhs_uz_xr,1,0,uz_xr[1]); 
  gkyl_mat_set(&rhs_Txx_xl,1,0,3.0*Txx_xl[1]); 
  gkyl_mat_set(&rhs_Txx_xr,1,0,3.0*Txx_xr[1]); 
 
  gkyl_mat_set(&rhs_ux_yl,1,0,ux_yl[1]); 
  gkyl_mat_set(&rhs_ux_yr,1,0,ux_yr[1]); 
  gkyl_mat_set(&rhs_uy_yl,1,0,uy_yl[1]); 
  gkyl_mat_set(&rhs_uy_yr,1,0,uy_yr[1]); 
  gkyl_mat_set(&rhs_uz_yl,1,0,uz_yl[1]); 
  gkyl_mat_set(&rhs_uz_yr,1,0,uz_yr[1]); 
  gkyl_mat_set(&rhs_Tyy_yl,1,0,3.0*Tyy_yl[1]); 
  gkyl_mat_set(&rhs_Tyy_yr,1,0,3.0*Tyy_yr[1]); 
 
  gkyl_mat_set(&rhs_ux_zl,1,0,ux_zl[1]); 
  gkyl_mat_set(&rhs_ux_zr,1,0,ux_zr[1]); 
  gkyl_mat_set(&rhs_uy_zl,1,0,uy_zl[1]); 
  gkyl_mat_set(&rhs_uy_zr,1,0,uy_zr[1]); 
  gkyl_mat_set(&rhs_uz_zl,1,0,uz_zl[1]); 
  gkyl_mat_set(&rhs_uz_zr,1,0,uz_zr[1]); 
  gkyl_mat_set(&rhs_Tzz_zl,1,0,3.0*Tzz_zl[1]); 
  gkyl_mat_set(&rhs_Tzz_zr,1,0,3.0*Tzz_zr[1]); 
 
  gkyl_mat_set(&rhs_ux_xl,2,0,ux_xl[2]); 
  gkyl_mat_set(&rhs_ux_xr,2,0,ux_xr[2]); 
  gkyl_mat_set(&rhs_uy_xl,2,0,uy_xl[2]); 
  gkyl_mat_set(&rhs_uy_xr,2,0,uy_xr[2]); 
  gkyl_mat_set(&rhs_uz_xl,2,0,uz_xl[2]); 
  gkyl_mat_set(&rhs_uz_xr,2,0,uz_xr[2]); 
  gkyl_mat_set(&rhs_Txx_xl,2,0,3.0*Txx_xl[2]); 
  gkyl_mat_set(&rhs_Txx_xr,2,0,3.0*Txx_xr[2]); 
 
  gkyl_mat_set(&rhs_ux_yl,2,0,ux_yl[2]); 
  gkyl_mat_set(&rhs_ux_yr,2,0,ux_yr[2]); 
  gkyl_mat_set(&rhs_uy_yl,2,0,uy_yl[2]); 
  gkyl_mat_set(&rhs_uy_yr,2,0,uy_yr[2]); 
  gkyl_mat_set(&rhs_uz_yl,2,0,uz_yl[2]); 
  gkyl_mat_set(&rhs_uz_yr,2,0,uz_yr[2]); 
  gkyl_mat_set(&rhs_Tyy_yl,2,0,3.0*Tyy_yl[2]); 
  gkyl_mat_set(&rhs_Tyy_yr,2,0,3.0*Tyy_yr[2]); 
 
  gkyl_mat_set(&rhs_ux_zl,2,0,ux_zl[2]); 
  gkyl_mat_set(&rhs_ux_zr,2,0,ux_zr[2]); 
  gkyl_mat_set(&rhs_uy_zl,2,0,uy_zl[2]); 
  gkyl_mat_set(&rhs_uy_zr,2,0,uy_zr[2]); 
  gkyl_mat_set(&rhs_uz_zl,2,0,uz_zl[2]); 
  gkyl_mat_set(&rhs_uz_zr,2,0,uz_zr[2]); 
  gkyl_mat_set(&rhs_Tzz_zl,2,0,3.0*Tzz_zl[2]); 
  gkyl_mat_set(&rhs_Tzz_zr,2,0,3.0*Tzz_zr[2]); 
 
  gkyl_mat_set(&rhs_ux_xl,3,0,ux_xl[3]); 
  gkyl_mat_set(&rhs_ux_xr,3,0,ux_xr[3]); 
  gkyl_mat_set(&rhs_uy_xl,3,0,uy_xl[3]); 
  gkyl_mat_set(&rhs_uy_xr,3,0,uy_xr[3]); 
  gkyl_mat_set(&rhs_uz_xl,3,0,uz_xl[3]); 
  gkyl_mat_set(&rhs_uz_xr,3,0,uz_xr[3]); 
  gkyl_mat_set(&rhs_Txx_xl,3,0,3.0*Txx_xl[3]); 
  gkyl_mat_set(&rhs_Txx_xr,3,0,3.0*Txx_xr[3]); 
 
  gkyl_mat_set(&rhs_ux_yl,3,0,ux_yl[3]); 
  gkyl_mat_set(&rhs_ux_yr,3,0,ux_yr[3]); 
  gkyl_mat_set(&rhs_uy_yl,3,0,uy_yl[3]); 
  gkyl_mat_set(&rhs_uy_yr,3,0,uy_yr[3]); 
  gkyl_mat_set(&rhs_uz_yl,3,0,uz_yl[3]); 
  gkyl_mat_set(&rhs_uz_yr,3,0,uz_yr[3]); 
  gkyl_mat_set(&rhs_Tyy_yl,3,0,3.0*Tyy_yl[3]); 
  gkyl_mat_set(&rhs_Tyy_yr,3,0,3.0*Tyy_yr[3]); 
 
  gkyl_mat_set(&rhs_ux_zl,3,0,ux_zl[3]); 
  gkyl_mat_set(&rhs_ux_zr,3,0,ux_zr[3]); 
  gkyl_mat_set(&rhs_uy_zl,3,0,uy_zl[3]); 
  gkyl_mat_set(&rhs_uy_zr,3,0,uy_zr[3]); 
  gkyl_mat_set(&rhs_uz_zl,3,0,uz_zl[3]); 
  gkyl_mat_set(&rhs_uz_zr,3,0,uz_zr[3]); 
  gkyl_mat_set(&rhs_Tzz_zl,3,0,3.0*Tzz_zl[3]); 
  gkyl_mat_set(&rhs_Tzz_zr,3,0,3.0*Tzz_zr[3]); 
 
} 
