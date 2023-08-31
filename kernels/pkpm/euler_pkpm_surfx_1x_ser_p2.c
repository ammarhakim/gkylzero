#include <gkyl_euler_pkpm_kernels.h> 
#include <gkyl_basis_ser_1x_p2_surfx1_eval_quad.h> 
GKYL_CU_DH double euler_pkpm_surfx_1x_ser_p2(const double *w, const double *dxv, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r,
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r,
    const double *p_ij_surf_l, const double *p_ij_surf_c, const double *p_ij_surf_r,
    const double *euler_pkpm_l, const double *euler_pkpm_c, const double *euler_pkpm_r, 
    const double *pkpm_lax, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:                Cell-center coordinates.
  // dxv[NDIM]:              Cell spacing.
  // vlasov_pkpm_moms_l/c/r: Input pkpm moments in left/center/right cells.
  // prim_surf_l/c/r:        Input surface primitive variables [u_i, 3*T_ii/m] in left/center/right cells in each direction.
  // p_ij_surf_l/c/r:        Input surface expansion of p_ij = (p_par - p_perp) b_i b_j + p_perp g_ij.
  //                         [Pxx_xl, Pxx_xr, Pxy_xl, Pxy_xr, Pxz_xl, Pxz_xr, 
  //                          Pxy_yl, Pxy_yr, Pyy_yl, Pyy_yr, Pyz_yl, Pyz_yr, 
  //                          Pxz_zl, Pxz_zr, Pyz_zl, Pyz_zr, Pzz_zl, Pzz_zr] 
  // euler_pkpm_l/c/r:       Input [rho ux, rho uy, rho uz], Fluid input state vector in left/center/right cells.
  // pkpm_lax:               Surface expansion of pkpm Lax penalization: lambda_i = |u_i| + sqrt(3.0*T_ii/m).
  // out: Incremented output.

  const double dx1 = 2.0/dxv[0]; 

  const double *rhoux_l = &euler_pkpm_l[0]; 
  const double *rhouy_l = &euler_pkpm_l[3]; 
  const double *rhouz_l = &euler_pkpm_l[6]; 

  const double *rhoux_c = &euler_pkpm_c[0]; 
  const double *rhouy_c = &euler_pkpm_c[3]; 
  const double *rhouz_c = &euler_pkpm_c[6]; 

  const double *rhoux_r = &euler_pkpm_r[0]; 
  const double *rhouy_r = &euler_pkpm_r[3]; 
  const double *rhouz_r = &euler_pkpm_r[6]; 

  const double *rho_l = &vlasov_pkpm_moms_l[0]; 
  const double *rho_c = &vlasov_pkpm_moms_c[0]; 
  const double *rho_r = &vlasov_pkpm_moms_r[0]; 

  const double *ux_surf_lr = &prim_surf_l[1]; 
  const double *uy_surf_lr = &prim_surf_l[3]; 
  const double *uz_surf_lr = &prim_surf_l[5]; 

  const double *ux_surf_cl = &prim_surf_c[0]; 
  const double *uy_surf_cl = &prim_surf_c[2]; 
  const double *uz_surf_cl = &prim_surf_c[4]; 

  const double *ux_surf_cr = &prim_surf_c[1]; 
  const double *uy_surf_cr = &prim_surf_c[3]; 
  const double *uz_surf_cr = &prim_surf_c[5]; 

  const double *ux_surf_rl = &prim_surf_r[0]; 
  const double *uy_surf_rl = &prim_surf_r[2]; 
  const double *uz_surf_rl = &prim_surf_r[4]; 

  const double *Pxx_surf_lr = &p_ij_surf_l[1]; 
  const double *Pxy_surf_lr = &p_ij_surf_l[3]; 
  const double *Pxz_surf_lr = &p_ij_surf_l[5]; 

  const double *Pxx_surf_cl = &p_ij_surf_c[0]; 
  const double *Pxy_surf_cl = &p_ij_surf_c[2]; 
  const double *Pxz_surf_cl = &p_ij_surf_c[4]; 

  const double *Pxx_surf_cr = &p_ij_surf_c[1]; 
  const double *Pxy_surf_cr = &p_ij_surf_c[3]; 
  const double *Pxz_surf_cr = &p_ij_surf_c[5]; 

  const double *Pxx_surf_rl = &p_ij_surf_r[0]; 
  const double *Pxy_surf_rl = &p_ij_surf_r[2]; 
  const double *Pxz_surf_rl = &p_ij_surf_r[4]; 

  const double *pkpm_lax_l = &pkpm_lax[0]; 
  const double *pkpm_lax_r = &pkpm_lax[1]; 

  double *outrhou0 = &out[0]; 
  double *outrhou1 = &out[3]; 
  double *outrhou2 = &out[6]; 

  double Ghat_rhoux_l = 0.0; 
  double Ghat_rhoux_r = 0.0; 
  double Ghat_rhouy_l = 0.0; 
  double Ghat_rhouy_r = 0.0; 
  double Ghat_rhouz_l = 0.0; 
  double Ghat_rhouz_r = 0.0; 
  double uxl_r = ux_surf_lr[0]; 
  double uxc_l = ux_surf_cl[0]; 
  double uxc_r = ux_surf_cr[0]; 
  double uxr_l = ux_surf_rl[0]; 

  double Pxxl_r = Pxx_surf_lr[0]; 
  double Pxxc_l = Pxx_surf_cl[0]; 
  double Pxxc_r = Pxx_surf_cr[0]; 
  double Pxxr_l = Pxx_surf_rl[0]; 

  double uxl_r_sq = uxl_r*uxl_r; 
  double uxc_l_sq = uxc_l*uxc_l; 
  double uxc_r_sq = uxc_r*uxc_r; 
  double uxr_l_sq = uxr_l*uxr_l; 

  double uyl_r = uy_surf_lr[0]; 
  double uyc_l = uy_surf_cl[0]; 
  double uyc_r = uy_surf_cr[0]; 
  double uyr_l = uy_surf_rl[0]; 

  double Pxyl_r = Pxy_surf_lr[0]; 
  double Pxyc_l = Pxy_surf_cl[0]; 
  double Pxyc_r = Pxy_surf_cr[0]; 
  double Pxyr_l = Pxy_surf_rl[0]; 

  double uzl_r = uz_surf_lr[0]; 
  double uzc_l = uz_surf_cl[0]; 
  double uzc_r = uz_surf_cr[0]; 
  double uzr_l = uz_surf_rl[0]; 

  double Pxzl_r = Pxz_surf_lr[0]; 
  double Pxzc_l = Pxz_surf_cl[0]; 
  double Pxzc_r = Pxz_surf_cr[0]; 
  double Pxzr_l = Pxz_surf_rl[0]; 

  double max_speed_l = pkpm_lax_l[0]; 
  double max_speed_r = pkpm_lax_r[0]; 

  Ghat_rhoux_l = 0.1976423537605236*rho_l[2]*uxl_r_sq+0.1976423537605236*rho_c[2]*uxl_r_sq+0.1530931089239486*rho_l[1]*uxl_r_sq-0.1530931089239486*rho_c[1]*uxl_r_sq+0.0883883476483184*rho_l[0]*uxl_r_sq+0.0883883476483184*rho_c[0]*uxl_r_sq+0.3952847075210473*rho_l[2]*uxc_l*uxl_r+0.3952847075210473*rho_c[2]*uxc_l*uxl_r+0.3061862178478971*rho_l[1]*uxc_l*uxl_r-0.3061862178478971*rho_c[1]*uxc_l*uxl_r+0.1767766952966368*rho_l[0]*uxc_l*uxl_r+0.1767766952966368*rho_c[0]*uxc_l*uxl_r+0.1976423537605236*rho_l[2]*uxc_l_sq+0.1976423537605236*rho_c[2]*uxc_l_sq+0.1530931089239486*rho_l[1]*uxc_l_sq-0.1530931089239486*rho_c[1]*uxc_l_sq+0.0883883476483184*rho_l[0]*uxc_l_sq+0.0883883476483184*rho_c[0]*uxc_l_sq+0.7905694150420947*rhoux_l[2]*max_speed_l-0.7905694150420947*rhoux_c[2]*max_speed_l+0.6123724356957944*rhoux_l[1]*max_speed_l+0.6123724356957944*rhoux_c[1]*max_speed_l+0.3535533905932737*rhoux_l[0]*max_speed_l-0.3535533905932737*rhoux_c[0]*max_speed_l+0.5*Pxxl_r+0.5*Pxxc_l; 
  Ghat_rhouy_l = 0.1976423537605237*rho_l[2]*uxl_r*uyl_r+0.1976423537605237*rho_c[2]*uxl_r*uyl_r+0.1530931089239486*rho_l[1]*uxl_r*uyl_r-0.1530931089239486*rho_c[1]*uxl_r*uyl_r+0.08838834764831843*rho_l[0]*uxl_r*uyl_r+0.08838834764831843*rho_c[0]*uxl_r*uyl_r+0.1976423537605237*rho_l[2]*uxc_l*uyl_r+0.1976423537605237*rho_c[2]*uxc_l*uyl_r+0.1530931089239486*rho_l[1]*uxc_l*uyl_r-0.1530931089239486*rho_c[1]*uxc_l*uyl_r+0.08838834764831843*rho_l[0]*uxc_l*uyl_r+0.08838834764831843*rho_c[0]*uxc_l*uyl_r+0.1976423537605237*rho_l[2]*uxl_r*uyc_l+0.1976423537605237*rho_c[2]*uxl_r*uyc_l+0.1530931089239486*rho_l[1]*uxl_r*uyc_l-0.1530931089239486*rho_c[1]*uxl_r*uyc_l+0.08838834764831843*rho_l[0]*uxl_r*uyc_l+0.08838834764831843*rho_c[0]*uxl_r*uyc_l+0.1976423537605237*rho_l[2]*uxc_l*uyc_l+0.1976423537605237*rho_c[2]*uxc_l*uyc_l+0.1530931089239486*rho_l[1]*uxc_l*uyc_l-0.1530931089239486*rho_c[1]*uxc_l*uyc_l+0.08838834764831843*rho_l[0]*uxc_l*uyc_l+0.08838834764831843*rho_c[0]*uxc_l*uyc_l+0.7905694150420948*rhouy_l[2]*max_speed_l-0.7905694150420948*rhouy_c[2]*max_speed_l+0.6123724356957945*rhouy_l[1]*max_speed_l+0.6123724356957945*rhouy_c[1]*max_speed_l+0.3535533905932737*rhouy_l[0]*max_speed_l-0.3535533905932737*rhouy_c[0]*max_speed_l+0.5*Pxyl_r+0.5*Pxyc_l; 
  Ghat_rhouz_l = 0.1976423537605237*rho_l[2]*uxl_r*uzl_r+0.1976423537605237*rho_c[2]*uxl_r*uzl_r+0.1530931089239486*rho_l[1]*uxl_r*uzl_r-0.1530931089239486*rho_c[1]*uxl_r*uzl_r+0.08838834764831843*rho_l[0]*uxl_r*uzl_r+0.08838834764831843*rho_c[0]*uxl_r*uzl_r+0.1976423537605237*rho_l[2]*uxc_l*uzl_r+0.1976423537605237*rho_c[2]*uxc_l*uzl_r+0.1530931089239486*rho_l[1]*uxc_l*uzl_r-0.1530931089239486*rho_c[1]*uxc_l*uzl_r+0.08838834764831843*rho_l[0]*uxc_l*uzl_r+0.08838834764831843*rho_c[0]*uxc_l*uzl_r+0.1976423537605237*rho_l[2]*uxl_r*uzc_l+0.1976423537605237*rho_c[2]*uxl_r*uzc_l+0.1530931089239486*rho_l[1]*uxl_r*uzc_l-0.1530931089239486*rho_c[1]*uxl_r*uzc_l+0.08838834764831843*rho_l[0]*uxl_r*uzc_l+0.08838834764831843*rho_c[0]*uxl_r*uzc_l+0.1976423537605237*rho_l[2]*uxc_l*uzc_l+0.1976423537605237*rho_c[2]*uxc_l*uzc_l+0.1530931089239486*rho_l[1]*uxc_l*uzc_l-0.1530931089239486*rho_c[1]*uxc_l*uzc_l+0.08838834764831843*rho_l[0]*uxc_l*uzc_l+0.08838834764831843*rho_c[0]*uxc_l*uzc_l+0.7905694150420948*rhouz_l[2]*max_speed_l-0.7905694150420948*rhouz_c[2]*max_speed_l+0.6123724356957945*rhouz_l[1]*max_speed_l+0.6123724356957945*rhouz_c[1]*max_speed_l+0.3535533905932737*rhouz_l[0]*max_speed_l-0.3535533905932737*rhouz_c[0]*max_speed_l+0.5*Pxzl_r+0.5*Pxzc_l; 
  Ghat_rhoux_r = 0.1976423537605236*rho_r[2]*uxr_l_sq+0.1976423537605236*rho_c[2]*uxr_l_sq-0.1530931089239486*rho_r[1]*uxr_l_sq+0.1530931089239486*rho_c[1]*uxr_l_sq+0.0883883476483184*rho_r[0]*uxr_l_sq+0.0883883476483184*rho_c[0]*uxr_l_sq+0.3952847075210473*rho_r[2]*uxc_r*uxr_l+0.3952847075210473*rho_c[2]*uxc_r*uxr_l-0.3061862178478971*rho_r[1]*uxc_r*uxr_l+0.3061862178478971*rho_c[1]*uxc_r*uxr_l+0.1767766952966368*rho_r[0]*uxc_r*uxr_l+0.1767766952966368*rho_c[0]*uxc_r*uxr_l+0.1976423537605236*rho_r[2]*uxc_r_sq+0.1976423537605236*rho_c[2]*uxc_r_sq-0.1530931089239486*rho_r[1]*uxc_r_sq+0.1530931089239486*rho_c[1]*uxc_r_sq+0.0883883476483184*rho_r[0]*uxc_r_sq+0.0883883476483184*rho_c[0]*uxc_r_sq-0.7905694150420947*rhoux_r[2]*max_speed_r+0.7905694150420947*rhoux_c[2]*max_speed_r+0.6123724356957944*rhoux_r[1]*max_speed_r+0.6123724356957944*rhoux_c[1]*max_speed_r-0.3535533905932737*rhoux_r[0]*max_speed_r+0.3535533905932737*rhoux_c[0]*max_speed_r+0.5*Pxxr_l+0.5*Pxxc_r; 
  Ghat_rhouy_r = 0.1976423537605237*rho_r[2]*uxr_l*uyr_l+0.1976423537605237*rho_c[2]*uxr_l*uyr_l-0.1530931089239486*rho_r[1]*uxr_l*uyr_l+0.1530931089239486*rho_c[1]*uxr_l*uyr_l+0.08838834764831843*rho_r[0]*uxr_l*uyr_l+0.08838834764831843*rho_c[0]*uxr_l*uyr_l+0.1976423537605237*rho_r[2]*uxc_r*uyr_l+0.1976423537605237*rho_c[2]*uxc_r*uyr_l-0.1530931089239486*rho_r[1]*uxc_r*uyr_l+0.1530931089239486*rho_c[1]*uxc_r*uyr_l+0.08838834764831843*rho_r[0]*uxc_r*uyr_l+0.08838834764831843*rho_c[0]*uxc_r*uyr_l+0.1976423537605237*rho_r[2]*uxr_l*uyc_r+0.1976423537605237*rho_c[2]*uxr_l*uyc_r-0.1530931089239486*rho_r[1]*uxr_l*uyc_r+0.1530931089239486*rho_c[1]*uxr_l*uyc_r+0.08838834764831843*rho_r[0]*uxr_l*uyc_r+0.08838834764831843*rho_c[0]*uxr_l*uyc_r+0.1976423537605237*rho_r[2]*uxc_r*uyc_r+0.1976423537605237*rho_c[2]*uxc_r*uyc_r-0.1530931089239486*rho_r[1]*uxc_r*uyc_r+0.1530931089239486*rho_c[1]*uxc_r*uyc_r+0.08838834764831843*rho_r[0]*uxc_r*uyc_r+0.08838834764831843*rho_c[0]*uxc_r*uyc_r-0.7905694150420948*rhouy_r[2]*max_speed_r+0.7905694150420948*rhouy_c[2]*max_speed_r+0.6123724356957945*rhouy_r[1]*max_speed_r+0.6123724356957945*rhouy_c[1]*max_speed_r-0.3535533905932737*rhouy_r[0]*max_speed_r+0.3535533905932737*rhouy_c[0]*max_speed_r+0.5*Pxyr_l+0.5*Pxyc_r; 
  Ghat_rhouz_r = 0.1976423537605237*rho_r[2]*uxr_l*uzr_l+0.1976423537605237*rho_c[2]*uxr_l*uzr_l-0.1530931089239486*rho_r[1]*uxr_l*uzr_l+0.1530931089239486*rho_c[1]*uxr_l*uzr_l+0.08838834764831843*rho_r[0]*uxr_l*uzr_l+0.08838834764831843*rho_c[0]*uxr_l*uzr_l+0.1976423537605237*rho_r[2]*uxc_r*uzr_l+0.1976423537605237*rho_c[2]*uxc_r*uzr_l-0.1530931089239486*rho_r[1]*uxc_r*uzr_l+0.1530931089239486*rho_c[1]*uxc_r*uzr_l+0.08838834764831843*rho_r[0]*uxc_r*uzr_l+0.08838834764831843*rho_c[0]*uxc_r*uzr_l+0.1976423537605237*rho_r[2]*uxr_l*uzc_r+0.1976423537605237*rho_c[2]*uxr_l*uzc_r-0.1530931089239486*rho_r[1]*uxr_l*uzc_r+0.1530931089239486*rho_c[1]*uxr_l*uzc_r+0.08838834764831843*rho_r[0]*uxr_l*uzc_r+0.08838834764831843*rho_c[0]*uxr_l*uzc_r+0.1976423537605237*rho_r[2]*uxc_r*uzc_r+0.1976423537605237*rho_c[2]*uxc_r*uzc_r-0.1530931089239486*rho_r[1]*uxc_r*uzc_r+0.1530931089239486*rho_c[1]*uxc_r*uzc_r+0.08838834764831843*rho_r[0]*uxc_r*uzc_r+0.08838834764831843*rho_c[0]*uxc_r*uzc_r-0.7905694150420948*rhouz_r[2]*max_speed_r+0.7905694150420948*rhouz_c[2]*max_speed_r+0.6123724356957945*rhouz_r[1]*max_speed_r+0.6123724356957945*rhouz_c[1]*max_speed_r-0.3535533905932737*rhouz_r[0]*max_speed_r+0.3535533905932737*rhouz_c[0]*max_speed_r+0.5*Pxzr_l+0.5*Pxzc_r; 

  outrhou0[0] += (0.7071067811865475*Ghat_rhoux_l-0.7071067811865475*Ghat_rhoux_r)*dx1; 
  outrhou0[1] += -1.224744871391589*(Ghat_rhoux_r+Ghat_rhoux_l)*dx1; 
  outrhou0[2] += (1.58113883008419*Ghat_rhoux_l-1.58113883008419*Ghat_rhoux_r)*dx1; 

  outrhou1[0] += (0.7071067811865475*Ghat_rhouy_l-0.7071067811865475*Ghat_rhouy_r)*dx1; 
  outrhou1[1] += -1.224744871391589*(Ghat_rhouy_r+Ghat_rhouy_l)*dx1; 
  outrhou1[2] += (1.58113883008419*Ghat_rhouy_l-1.58113883008419*Ghat_rhouy_r)*dx1; 

  outrhou2[0] += (0.7071067811865475*Ghat_rhouz_l-0.7071067811865475*Ghat_rhouz_r)*dx1; 
  outrhou2[1] += -1.224744871391589*(Ghat_rhouz_r+Ghat_rhouz_l)*dx1; 
  outrhou2[2] += (1.58113883008419*Ghat_rhouz_l-1.58113883008419*Ghat_rhouz_r)*dx1; 

  return 0.;

} 
