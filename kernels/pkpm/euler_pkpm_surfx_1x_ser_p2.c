#include <gkyl_euler_pkpm_kernels.h> 
#include <gkyl_basis_ser_1x_p2_surfx1_eval_quad.h> 
GKYL_CU_DH double euler_pkpm_surfx_1x_ser_p2(const double *w, const double *dxv, 
    const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom_l, const struct gkyl_wave_cell_geom *geom_r, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r,
    const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r,
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r,
    const double *euler_pkpm_l, const double *euler_pkpm_c, const double *euler_pkpm_r, 
    const double *pkpm_lax, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:                Cell-center coordinates.
  // dxv[NDIM]:              Cell spacing.
  // wv_eqn:                 Wave equation for computing fluctuations at the interface for upwinding.
  // geom_l:                 Geometry for the left surface update.
  // geom_r:                 Geometry for the right surface update.
  // vlasov_pkpm_moms_l/c/r: Input pkpm moments in left/center/right cells.
  // prim_surf_l/c/r:        Input surface primitive variables [u_i, 3*T_ii/m] in left/center/right cells in each direction.
  // p_ij_l/c/r:             Input volume expansion of p_ij = (p_par - p_perp) b_i b_j + p_perp g_ij in left/center/right cells.
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

  const double *Pxx_l = &p_ij_l[0]; 
  const double *Pxy_l = &p_ij_l[3]; 
  const double *Pxz_l = &p_ij_l[6]; 
  const double *Pyy_l = &p_ij_l[9]; 
  const double *Pyz_l = &p_ij_l[12]; 
  const double *Pzz_l = &p_ij_l[15]; 

  const double *Pxx_c = &p_ij_c[0]; 
  const double *Pxy_c = &p_ij_c[3]; 
  const double *Pxz_c = &p_ij_c[6]; 
  const double *Pyy_c = &p_ij_c[9]; 
  const double *Pyz_c = &p_ij_c[12]; 
  const double *Pzz_c = &p_ij_c[15]; 

  const double *Pxx_r = &p_ij_r[0]; 
  const double *Pxy_r = &p_ij_r[3]; 
  const double *Pxz_r = &p_ij_r[6]; 
  const double *Pyy_r = &p_ij_r[9]; 
  const double *Pyz_r = &p_ij_r[12]; 
  const double *Pzz_r = &p_ij_r[15]; 

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
  double q_lr[10] = {0.0}; 
  double q_cl[10] = {0.0}; 
  double q_cr[10] = {0.0}; 
  double q_rl[10] = {0.0}; 
  q_lr[0] = 1.58113883008419*rho_l[2]+1.224744871391589*rho_l[1]+0.7071067811865475*rho_l[0]; 
  q_lr[1] = 1.58113883008419*rhoux_l[2]+1.224744871391589*rhoux_l[1]+0.7071067811865475*rhoux_l[0]; 
  q_lr[2] = 1.58113883008419*rhouy_l[2]+1.224744871391589*rhouy_l[1]+0.7071067811865475*rhouy_l[0]; 
  q_lr[3] = 1.58113883008419*rhouz_l[2]+1.224744871391589*rhouz_l[1]+0.7071067811865475*rhouz_l[0]; 
  q_lr[4] = 1.58113883008419*Pxx_l[2]+1.224744871391589*Pxx_l[1]+0.7071067811865475*Pxx_l[0] + q_lr[1]*q_lr[1]/q_lr[0]; 
  q_lr[5] = 1.58113883008419*Pxy_l[2]+1.224744871391589*Pxy_l[1]+0.7071067811865475*Pxy_l[0] + q_lr[1]*q_lr[2]/q_lr[0]; 
  q_lr[6] = 1.58113883008419*Pxz_l[2]+1.224744871391589*Pxz_l[1]+0.7071067811865475*Pxz_l[0] + q_lr[1]*q_lr[3]/q_lr[0]; 
  q_lr[7] = 1.58113883008419*Pyy_l[2]+1.224744871391589*Pyy_l[1]+0.7071067811865475*Pyy_l[0] + q_lr[2]*q_lr[2]/q_lr[0]; 
  q_lr[8] = 1.58113883008419*Pyz_l[2]+1.224744871391589*Pyz_l[1]+0.7071067811865475*Pyz_l[0] + q_lr[2]*q_lr[3]/q_lr[0]; 
  q_lr[9] = 1.58113883008419*Pzz_l[2]+1.224744871391589*Pzz_l[1]+0.7071067811865475*Pzz_l[0] + q_lr[3]*q_lr[3]/q_lr[0]; 
  q_cl[0] = 1.58113883008419*rho_c[2]-1.224744871391589*rho_c[1]+0.7071067811865475*rho_c[0]; 
  q_cl[1] = 1.58113883008419*rhoux_c[2]-1.224744871391589*rhoux_c[1]+0.7071067811865475*rhoux_c[0]; 
  q_cl[2] = 1.58113883008419*rhouy_c[2]-1.224744871391589*rhouy_c[1]+0.7071067811865475*rhouy_c[0]; 
  q_cl[3] = 1.58113883008419*rhouz_c[2]-1.224744871391589*rhouz_c[1]+0.7071067811865475*rhouz_c[0]; 
  q_cl[4] = 1.58113883008419*Pxx_c[2]-1.224744871391589*Pxx_c[1]+0.7071067811865475*Pxx_c[0] + q_cl[1]*q_cl[1]/q_cl[0]; 
  q_cl[5] = 1.58113883008419*Pxy_c[2]-1.224744871391589*Pxy_c[1]+0.7071067811865475*Pxy_c[0] + q_cl[1]*q_cl[2]/q_cl[0]; 
  q_cl[6] = 1.58113883008419*Pxz_c[2]-1.224744871391589*Pxz_c[1]+0.7071067811865475*Pxz_c[0] + q_cl[1]*q_cl[3]/q_cl[0]; 
  q_cl[7] = 1.58113883008419*Pyy_c[2]-1.224744871391589*Pyy_c[1]+0.7071067811865475*Pyy_c[0] + q_cl[2]*q_cl[2]/q_cl[0]; 
  q_cl[8] = 1.58113883008419*Pyz_c[2]-1.224744871391589*Pyz_c[1]+0.7071067811865475*Pyz_c[0] + q_cl[2]*q_cl[3]/q_cl[0]; 
  q_cl[9] = 1.58113883008419*Pzz_c[2]-1.224744871391589*Pzz_c[1]+0.7071067811865475*Pzz_c[0] + q_cl[3]*q_cl[3]/q_cl[0]; 
  q_cr[0] = 1.58113883008419*rho_c[2]+1.224744871391589*rho_c[1]+0.7071067811865475*rho_c[0]; 
  q_cr[1] = 1.58113883008419*rhoux_c[2]+1.224744871391589*rhoux_c[1]+0.7071067811865475*rhoux_c[0]; 
  q_cr[2] = 1.58113883008419*rhouy_c[2]+1.224744871391589*rhouy_c[1]+0.7071067811865475*rhouy_c[0]; 
  q_cr[3] = 1.58113883008419*rhouz_c[2]+1.224744871391589*rhouz_c[1]+0.7071067811865475*rhouz_c[0]; 
  q_cr[4] = 1.58113883008419*Pxx_c[2]+1.224744871391589*Pxx_c[1]+0.7071067811865475*Pxx_c[0] + q_cr[1]*q_cr[1]/q_cr[0]; 
  q_cr[5] = 1.58113883008419*Pxy_c[2]+1.224744871391589*Pxy_c[1]+0.7071067811865475*Pxy_c[0] + q_cr[1]*q_cr[2]/q_cr[0]; 
  q_cr[6] = 1.58113883008419*Pxz_c[2]+1.224744871391589*Pxz_c[1]+0.7071067811865475*Pxz_c[0] + q_cr[1]*q_cr[3]/q_cr[0]; 
  q_cr[7] = 1.58113883008419*Pyy_c[2]+1.224744871391589*Pyy_c[1]+0.7071067811865475*Pyy_c[0] + q_cr[2]*q_cr[2]/q_cr[0]; 
  q_cr[8] = 1.58113883008419*Pyz_c[2]+1.224744871391589*Pyz_c[1]+0.7071067811865475*Pyz_c[0] + q_cr[2]*q_cr[3]/q_cr[0]; 
  q_cr[9] = 1.58113883008419*Pzz_c[2]+1.224744871391589*Pzz_c[1]+0.7071067811865475*Pzz_c[0] + q_cr[3]*q_cr[3]/q_cr[0]; 
  q_rl[0] = 1.58113883008419*rho_r[2]-1.224744871391589*rho_r[1]+0.7071067811865475*rho_r[0]; 
  q_rl[1] = 1.58113883008419*rhoux_r[2]-1.224744871391589*rhoux_r[1]+0.7071067811865475*rhoux_r[0]; 
  q_rl[2] = 1.58113883008419*rhouy_r[2]-1.224744871391589*rhouy_r[1]+0.7071067811865475*rhouy_r[0]; 
  q_rl[3] = 1.58113883008419*rhouz_r[2]-1.224744871391589*rhouz_r[1]+0.7071067811865475*rhouz_r[0]; 
  q_rl[4] = 1.58113883008419*Pxx_r[2]-1.224744871391589*Pxx_r[1]+0.7071067811865475*Pxx_r[0] + q_rl[1]*q_rl[1]/q_rl[0]; 
  q_rl[5] = 1.58113883008419*Pxy_r[2]-1.224744871391589*Pxy_r[1]+0.7071067811865475*Pxy_r[0] + q_rl[1]*q_rl[2]/q_rl[0]; 
  q_rl[6] = 1.58113883008419*Pxz_r[2]-1.224744871391589*Pxz_r[1]+0.7071067811865475*Pxz_r[0] + q_rl[1]*q_rl[3]/q_rl[0]; 
  q_rl[7] = 1.58113883008419*Pyy_r[2]-1.224744871391589*Pyy_r[1]+0.7071067811865475*Pyy_r[0] + q_rl[2]*q_rl[2]/q_rl[0]; 
  q_rl[8] = 1.58113883008419*Pyz_r[2]-1.224744871391589*Pyz_r[1]+0.7071067811865475*Pyz_r[0] + q_rl[2]*q_rl[3]/q_rl[0]; 
  q_rl[9] = 1.58113883008419*Pzz_r[2]-1.224744871391589*Pzz_r[1]+0.7071067811865475*Pzz_r[0] + q_rl[3]*q_rl[3]/q_rl[0]; 

  double q_lr_local[10] = {0.0}; 
  double q_cl_local[10] = {0.0}; 
  double q_cr_local[10] = {0.0}; 
  double q_rl_local[10] = {0.0}; 
  gkyl_wv_eqn_rotate_to_local(wv_eqn, geom_l->tau1[0], geom_l->tau2[0], geom_l->norm[0], q_lr, q_lr_local); 
  gkyl_wv_eqn_rotate_to_local(wv_eqn, geom_l->tau1[0], geom_l->tau2[0], geom_l->norm[0], q_cl, q_cl_local); 
  gkyl_wv_eqn_rotate_to_local(wv_eqn, geom_r->tau1[0], geom_r->tau2[0], geom_r->norm[0], q_cr, q_cr_local); 
  gkyl_wv_eqn_rotate_to_local(wv_eqn, geom_r->tau1[0], geom_r->tau2[0], geom_r->norm[0], q_rl, q_rl_local); 

  double delta_l[10] = {0.0}; 
  double delta_r[10] = {0.0}; 
  delta_l[0] = q_cl_local[0] - q_lr_local[0]; 
  delta_l[1] = q_cl_local[1] - q_lr_local[1]; 
  delta_l[2] = q_cl_local[2] - q_lr_local[2]; 
  delta_l[3] = q_cl_local[3] - q_lr_local[3]; 
  delta_l[4] = q_cl_local[4] - q_lr_local[4]; 
  delta_l[5] = q_cl_local[5] - q_lr_local[5]; 
  delta_l[6] = q_cl_local[6] - q_lr_local[6]; 
  delta_l[7] = q_cl_local[7] - q_lr_local[7]; 
  delta_l[8] = q_cl_local[8] - q_lr_local[8]; 
  delta_l[9] = q_cl_local[9] - q_lr_local[9]; 
  delta_r[0] = q_rl_local[0] - q_cr_local[0]; 
  delta_r[1] = q_rl_local[1] - q_cr_local[1]; 
  delta_r[2] = q_rl_local[2] - q_cr_local[2]; 
  delta_r[3] = q_rl_local[3] - q_cr_local[3]; 
  delta_r[4] = q_rl_local[4] - q_cr_local[4]; 
  delta_r[5] = q_rl_local[5] - q_cr_local[5]; 
  delta_r[6] = q_rl_local[6] - q_cr_local[6]; 
  delta_r[7] = q_rl_local[7] - q_cr_local[7]; 
  delta_r[8] = q_rl_local[8] - q_cr_local[8]; 
  delta_r[9] = q_rl_local[9] - q_cr_local[9]; 

  double waves_l[50] = {0.0}; 
  double waves_r[50] = {0.0}; 
  double speeds_l[5] = {0.0}; 
  double speeds_r[5] = {0.0}; 
  double my_max_speed_l = gkyl_wv_eqn_waves(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, delta_l, q_lr_local, q_cl_local, waves_l, speeds_l); 
  double my_max_speed_r = gkyl_wv_eqn_waves(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, delta_r, q_cr_local, q_rl_local, waves_r, speeds_r); 

  double lenr_l = geom_l->lenr[0]; 
  speeds_l[0] *= lenr_l; 
  speeds_l[1] *= lenr_l; 
  speeds_l[2] *= lenr_l; 
  speeds_l[3] *= lenr_l; 
  speeds_l[4] *= lenr_l; 
  double lenr_r = geom_r->lenr[0]; 
  speeds_r[0] *= lenr_r; 
  speeds_r[1] *= lenr_r; 
  speeds_r[2] *= lenr_r; 
  speeds_r[3] *= lenr_r; 
  speeds_r[4] *= lenr_r; 

  double amdq_l_local[10] = {0.0}; 
  double apdq_l_local[10] = {0.0}; 
  double amdq_r_local[10] = {0.0}; 
  double apdq_r_local[10] = {0.0}; 
  gkyl_wv_eqn_qfluct(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, q_lr_local, q_cl_local, waves_l, speeds_l, amdq_l_local, apdq_l_local); 
  gkyl_wv_eqn_qfluct(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, q_cr_local, q_rl_local, waves_r, speeds_r, amdq_r_local, apdq_r_local); 

  double amdq_l[10] = {0.0}; 
  double apdq_l[10] = {0.0}; 
  double amdq_r[10] = {0.0}; 
  double apdq_r[10] = {0.0}; 
  gkyl_wv_eqn_rotate_to_global(wv_eqn, geom_l->tau1[0], geom_l->tau2[0], geom_l->norm[0], amdq_l_local, amdq_l); 
  gkyl_wv_eqn_rotate_to_global(wv_eqn, geom_l->tau1[0], geom_l->tau2[0], geom_l->norm[0], apdq_l_local, apdq_l); 
  gkyl_wv_eqn_rotate_to_global(wv_eqn, geom_r->tau1[0], geom_r->tau2[0], geom_r->norm[0], amdq_r_local, amdq_r); 
  gkyl_wv_eqn_rotate_to_global(wv_eqn, geom_r->tau1[0], geom_r->tau2[0], geom_r->norm[0], apdq_r_local, apdq_r); 

  double uxl_r = ux_surf_lr[0]; 
  double uxc_l = ux_surf_cl[0]; 
  double uxc_r = ux_surf_cr[0]; 
  double uxr_l = ux_surf_rl[0]; 

  double uxl_r_sq = uxl_r*uxl_r; 
  double uxc_l_sq = uxc_l*uxc_l; 
  double uxc_r_sq = uxc_r*uxc_r; 
  double uxr_l_sq = uxr_l*uxr_l; 

  double uyl_r = uy_surf_lr[0]; 
  double uyc_l = uy_surf_cl[0]; 
  double uyc_r = uy_surf_cr[0]; 
  double uyr_l = uy_surf_rl[0]; 

  double uzl_r = uz_surf_lr[0]; 
  double uzc_l = uz_surf_cl[0]; 
  double uzc_r = uz_surf_cr[0]; 
  double uzr_l = uz_surf_rl[0]; 

  double max_speed_l = pkpm_lax_l[0]; 
  double max_speed_r = pkpm_lax_r[0]; 

  Ghat_rhoux_l = 0.1976423537605236*rho_l[2]*uxl_r_sq+0.1976423537605236*rho_c[2]*uxl_r_sq+0.1530931089239486*rho_l[1]*uxl_r_sq-0.1530931089239486*rho_c[1]*uxl_r_sq+0.0883883476483184*rho_l[0]*uxl_r_sq+0.0883883476483184*rho_c[0]*uxl_r_sq+0.3952847075210473*rho_l[2]*uxc_l*uxl_r+0.3952847075210473*rho_c[2]*uxc_l*uxl_r+0.3061862178478971*rho_l[1]*uxc_l*uxl_r-0.3061862178478971*rho_c[1]*uxc_l*uxl_r+0.1767766952966368*rho_l[0]*uxc_l*uxl_r+0.1767766952966368*rho_c[0]*uxc_l*uxl_r+0.3952847075210473*rho_l[2]*max_speed_l*uxl_r-0.3952847075210473*rho_c[2]*max_speed_l*uxl_r+0.3061862178478971*rho_l[1]*max_speed_l*uxl_r+0.3061862178478971*rho_c[1]*max_speed_l*uxl_r+0.1767766952966368*rho_l[0]*max_speed_l*uxl_r-0.1767766952966368*rho_c[0]*max_speed_l*uxl_r+0.1976423537605236*rho_l[2]*uxc_l_sq+0.1976423537605236*rho_c[2]*uxc_l_sq+0.1530931089239486*rho_l[1]*uxc_l_sq-0.1530931089239486*rho_c[1]*uxc_l_sq+0.0883883476483184*rho_l[0]*uxc_l_sq+0.0883883476483184*rho_c[0]*uxc_l_sq+0.3952847075210473*rho_l[2]*max_speed_l*uxc_l-0.3952847075210473*rho_c[2]*max_speed_l*uxc_l+0.3061862178478971*rho_l[1]*max_speed_l*uxc_l+0.3061862178478971*rho_c[1]*max_speed_l*uxc_l+0.1767766952966368*rho_l[0]*max_speed_l*uxc_l-0.1767766952966368*rho_c[0]*max_speed_l*uxc_l+0.7905694150420947*Pxx_l[2]+0.7905694150420947*Pxx_c[2]-0.5*apdq_l[1]+0.5*amdq_l[1]+0.6123724356957944*Pxx_l[1]-0.6123724356957944*Pxx_c[1]+0.3535533905932737*Pxx_l[0]+0.3535533905932737*Pxx_c[0]; 
  Ghat_rhouy_l = 0.1976423537605237*rho_l[2]*uxl_r*uyl_r+0.1976423537605237*rho_c[2]*uxl_r*uyl_r+0.1530931089239486*rho_l[1]*uxl_r*uyl_r-0.1530931089239486*rho_c[1]*uxl_r*uyl_r+0.08838834764831843*rho_l[0]*uxl_r*uyl_r+0.08838834764831843*rho_c[0]*uxl_r*uyl_r+0.1976423537605237*rho_l[2]*uxc_l*uyl_r+0.1976423537605237*rho_c[2]*uxc_l*uyl_r+0.1530931089239486*rho_l[1]*uxc_l*uyl_r-0.1530931089239486*rho_c[1]*uxc_l*uyl_r+0.08838834764831843*rho_l[0]*uxc_l*uyl_r+0.08838834764831843*rho_c[0]*uxc_l*uyl_r+0.3952847075210474*rho_l[2]*max_speed_l*uyl_r-0.3952847075210474*rho_c[2]*max_speed_l*uyl_r+0.3061862178478972*rho_l[1]*max_speed_l*uyl_r+0.3061862178478972*rho_c[1]*max_speed_l*uyl_r+0.1767766952966369*rho_l[0]*max_speed_l*uyl_r-0.1767766952966369*rho_c[0]*max_speed_l*uyl_r+0.1976423537605237*rho_l[2]*uxl_r*uyc_l+0.1976423537605237*rho_c[2]*uxl_r*uyc_l+0.1530931089239486*rho_l[1]*uxl_r*uyc_l-0.1530931089239486*rho_c[1]*uxl_r*uyc_l+0.08838834764831843*rho_l[0]*uxl_r*uyc_l+0.08838834764831843*rho_c[0]*uxl_r*uyc_l+0.1976423537605237*rho_l[2]*uxc_l*uyc_l+0.1976423537605237*rho_c[2]*uxc_l*uyc_l+0.1530931089239486*rho_l[1]*uxc_l*uyc_l-0.1530931089239486*rho_c[1]*uxc_l*uyc_l+0.08838834764831843*rho_l[0]*uxc_l*uyc_l+0.08838834764831843*rho_c[0]*uxc_l*uyc_l+0.3952847075210474*rho_l[2]*max_speed_l*uyc_l-0.3952847075210474*rho_c[2]*max_speed_l*uyc_l+0.3061862178478972*rho_l[1]*max_speed_l*uyc_l+0.3061862178478972*rho_c[1]*max_speed_l*uyc_l+0.1767766952966369*rho_l[0]*max_speed_l*uyc_l-0.1767766952966369*rho_c[0]*max_speed_l*uyc_l-0.5*apdq_l[2]+0.5*amdq_l[2]+0.7905694150420948*Pxy_l[2]+0.7905694150420948*Pxy_c[2]+0.6123724356957945*Pxy_l[1]-0.6123724356957945*Pxy_c[1]+0.3535533905932737*Pxy_l[0]+0.3535533905932737*Pxy_c[0]; 
  Ghat_rhouz_l = 0.1976423537605237*rho_l[2]*uxl_r*uzl_r+0.1976423537605237*rho_c[2]*uxl_r*uzl_r+0.1530931089239486*rho_l[1]*uxl_r*uzl_r-0.1530931089239486*rho_c[1]*uxl_r*uzl_r+0.08838834764831843*rho_l[0]*uxl_r*uzl_r+0.08838834764831843*rho_c[0]*uxl_r*uzl_r+0.1976423537605237*rho_l[2]*uxc_l*uzl_r+0.1976423537605237*rho_c[2]*uxc_l*uzl_r+0.1530931089239486*rho_l[1]*uxc_l*uzl_r-0.1530931089239486*rho_c[1]*uxc_l*uzl_r+0.08838834764831843*rho_l[0]*uxc_l*uzl_r+0.08838834764831843*rho_c[0]*uxc_l*uzl_r+0.3952847075210474*rho_l[2]*max_speed_l*uzl_r-0.3952847075210474*rho_c[2]*max_speed_l*uzl_r+0.3061862178478972*rho_l[1]*max_speed_l*uzl_r+0.3061862178478972*rho_c[1]*max_speed_l*uzl_r+0.1767766952966369*rho_l[0]*max_speed_l*uzl_r-0.1767766952966369*rho_c[0]*max_speed_l*uzl_r+0.1976423537605237*rho_l[2]*uxl_r*uzc_l+0.1976423537605237*rho_c[2]*uxl_r*uzc_l+0.1530931089239486*rho_l[1]*uxl_r*uzc_l-0.1530931089239486*rho_c[1]*uxl_r*uzc_l+0.08838834764831843*rho_l[0]*uxl_r*uzc_l+0.08838834764831843*rho_c[0]*uxl_r*uzc_l+0.1976423537605237*rho_l[2]*uxc_l*uzc_l+0.1976423537605237*rho_c[2]*uxc_l*uzc_l+0.1530931089239486*rho_l[1]*uxc_l*uzc_l-0.1530931089239486*rho_c[1]*uxc_l*uzc_l+0.08838834764831843*rho_l[0]*uxc_l*uzc_l+0.08838834764831843*rho_c[0]*uxc_l*uzc_l+0.3952847075210474*rho_l[2]*max_speed_l*uzc_l-0.3952847075210474*rho_c[2]*max_speed_l*uzc_l+0.3061862178478972*rho_l[1]*max_speed_l*uzc_l+0.3061862178478972*rho_c[1]*max_speed_l*uzc_l+0.1767766952966369*rho_l[0]*max_speed_l*uzc_l-0.1767766952966369*rho_c[0]*max_speed_l*uzc_l-0.5*apdq_l[3]+0.5*amdq_l[3]+0.7905694150420948*Pxz_l[2]+0.7905694150420948*Pxz_c[2]+0.6123724356957945*Pxz_l[1]-0.6123724356957945*Pxz_c[1]+0.3535533905932737*Pxz_l[0]+0.3535533905932737*Pxz_c[0]; 
  Ghat_rhoux_r = 0.1976423537605236*rho_r[2]*uxr_l_sq+0.1976423537605236*rho_c[2]*uxr_l_sq-0.1530931089239486*rho_r[1]*uxr_l_sq+0.1530931089239486*rho_c[1]*uxr_l_sq+0.0883883476483184*rho_r[0]*uxr_l_sq+0.0883883476483184*rho_c[0]*uxr_l_sq+0.3952847075210473*rho_r[2]*uxc_r*uxr_l+0.3952847075210473*rho_c[2]*uxc_r*uxr_l-0.3061862178478971*rho_r[1]*uxc_r*uxr_l+0.3061862178478971*rho_c[1]*uxc_r*uxr_l+0.1767766952966368*rho_r[0]*uxc_r*uxr_l+0.1767766952966368*rho_c[0]*uxc_r*uxr_l-0.3952847075210473*rho_r[2]*max_speed_r*uxr_l+0.3952847075210473*rho_c[2]*max_speed_r*uxr_l+0.3061862178478971*rho_r[1]*max_speed_r*uxr_l+0.3061862178478971*rho_c[1]*max_speed_r*uxr_l-0.1767766952966368*rho_r[0]*max_speed_r*uxr_l+0.1767766952966368*rho_c[0]*max_speed_r*uxr_l+0.1976423537605236*rho_r[2]*uxc_r_sq+0.1976423537605236*rho_c[2]*uxc_r_sq-0.1530931089239486*rho_r[1]*uxc_r_sq+0.1530931089239486*rho_c[1]*uxc_r_sq+0.0883883476483184*rho_r[0]*uxc_r_sq+0.0883883476483184*rho_c[0]*uxc_r_sq-0.3952847075210473*rho_r[2]*max_speed_r*uxc_r+0.3952847075210473*rho_c[2]*max_speed_r*uxc_r+0.3061862178478971*rho_r[1]*max_speed_r*uxc_r+0.3061862178478971*rho_c[1]*max_speed_r*uxc_r-0.1767766952966368*rho_r[0]*max_speed_r*uxc_r+0.1767766952966368*rho_c[0]*max_speed_r*uxc_r+0.7905694150420947*Pxx_r[2]+0.7905694150420947*Pxx_c[2]-0.5*apdq_r[1]+0.5*amdq_r[1]-0.6123724356957944*Pxx_r[1]+0.6123724356957944*Pxx_c[1]+0.3535533905932737*Pxx_r[0]+0.3535533905932737*Pxx_c[0]; 
  Ghat_rhouy_r = 0.1976423537605237*rho_r[2]*uxr_l*uyr_l+0.1976423537605237*rho_c[2]*uxr_l*uyr_l-0.1530931089239486*rho_r[1]*uxr_l*uyr_l+0.1530931089239486*rho_c[1]*uxr_l*uyr_l+0.08838834764831843*rho_r[0]*uxr_l*uyr_l+0.08838834764831843*rho_c[0]*uxr_l*uyr_l+0.1976423537605237*rho_r[2]*uxc_r*uyr_l+0.1976423537605237*rho_c[2]*uxc_r*uyr_l-0.1530931089239486*rho_r[1]*uxc_r*uyr_l+0.1530931089239486*rho_c[1]*uxc_r*uyr_l+0.08838834764831843*rho_r[0]*uxc_r*uyr_l+0.08838834764831843*rho_c[0]*uxc_r*uyr_l-0.3952847075210474*rho_r[2]*max_speed_r*uyr_l+0.3952847075210474*rho_c[2]*max_speed_r*uyr_l+0.3061862178478972*rho_r[1]*max_speed_r*uyr_l+0.3061862178478972*rho_c[1]*max_speed_r*uyr_l-0.1767766952966369*rho_r[0]*max_speed_r*uyr_l+0.1767766952966369*rho_c[0]*max_speed_r*uyr_l+0.1976423537605237*rho_r[2]*uxr_l*uyc_r+0.1976423537605237*rho_c[2]*uxr_l*uyc_r-0.1530931089239486*rho_r[1]*uxr_l*uyc_r+0.1530931089239486*rho_c[1]*uxr_l*uyc_r+0.08838834764831843*rho_r[0]*uxr_l*uyc_r+0.08838834764831843*rho_c[0]*uxr_l*uyc_r+0.1976423537605237*rho_r[2]*uxc_r*uyc_r+0.1976423537605237*rho_c[2]*uxc_r*uyc_r-0.1530931089239486*rho_r[1]*uxc_r*uyc_r+0.1530931089239486*rho_c[1]*uxc_r*uyc_r+0.08838834764831843*rho_r[0]*uxc_r*uyc_r+0.08838834764831843*rho_c[0]*uxc_r*uyc_r-0.3952847075210474*rho_r[2]*max_speed_r*uyc_r+0.3952847075210474*rho_c[2]*max_speed_r*uyc_r+0.3061862178478972*rho_r[1]*max_speed_r*uyc_r+0.3061862178478972*rho_c[1]*max_speed_r*uyc_r-0.1767766952966369*rho_r[0]*max_speed_r*uyc_r+0.1767766952966369*rho_c[0]*max_speed_r*uyc_r-0.5*apdq_r[2]+0.5*amdq_r[2]+0.7905694150420948*Pxy_r[2]+0.7905694150420948*Pxy_c[2]-0.6123724356957945*Pxy_r[1]+0.6123724356957945*Pxy_c[1]+0.3535533905932737*Pxy_r[0]+0.3535533905932737*Pxy_c[0]; 
  Ghat_rhouz_r = 0.1976423537605237*rho_r[2]*uxr_l*uzr_l+0.1976423537605237*rho_c[2]*uxr_l*uzr_l-0.1530931089239486*rho_r[1]*uxr_l*uzr_l+0.1530931089239486*rho_c[1]*uxr_l*uzr_l+0.08838834764831843*rho_r[0]*uxr_l*uzr_l+0.08838834764831843*rho_c[0]*uxr_l*uzr_l+0.1976423537605237*rho_r[2]*uxc_r*uzr_l+0.1976423537605237*rho_c[2]*uxc_r*uzr_l-0.1530931089239486*rho_r[1]*uxc_r*uzr_l+0.1530931089239486*rho_c[1]*uxc_r*uzr_l+0.08838834764831843*rho_r[0]*uxc_r*uzr_l+0.08838834764831843*rho_c[0]*uxc_r*uzr_l-0.3952847075210474*rho_r[2]*max_speed_r*uzr_l+0.3952847075210474*rho_c[2]*max_speed_r*uzr_l+0.3061862178478972*rho_r[1]*max_speed_r*uzr_l+0.3061862178478972*rho_c[1]*max_speed_r*uzr_l-0.1767766952966369*rho_r[0]*max_speed_r*uzr_l+0.1767766952966369*rho_c[0]*max_speed_r*uzr_l+0.1976423537605237*rho_r[2]*uxr_l*uzc_r+0.1976423537605237*rho_c[2]*uxr_l*uzc_r-0.1530931089239486*rho_r[1]*uxr_l*uzc_r+0.1530931089239486*rho_c[1]*uxr_l*uzc_r+0.08838834764831843*rho_r[0]*uxr_l*uzc_r+0.08838834764831843*rho_c[0]*uxr_l*uzc_r+0.1976423537605237*rho_r[2]*uxc_r*uzc_r+0.1976423537605237*rho_c[2]*uxc_r*uzc_r-0.1530931089239486*rho_r[1]*uxc_r*uzc_r+0.1530931089239486*rho_c[1]*uxc_r*uzc_r+0.08838834764831843*rho_r[0]*uxc_r*uzc_r+0.08838834764831843*rho_c[0]*uxc_r*uzc_r-0.3952847075210474*rho_r[2]*max_speed_r*uzc_r+0.3952847075210474*rho_c[2]*max_speed_r*uzc_r+0.3061862178478972*rho_r[1]*max_speed_r*uzc_r+0.3061862178478972*rho_c[1]*max_speed_r*uzc_r-0.1767766952966369*rho_r[0]*max_speed_r*uzc_r+0.1767766952966369*rho_c[0]*max_speed_r*uzc_r-0.5*apdq_r[3]+0.5*amdq_r[3]+0.7905694150420948*Pxz_r[2]+0.7905694150420948*Pxz_c[2]-0.6123724356957945*Pxz_r[1]+0.6123724356957945*Pxz_c[1]+0.3535533905932737*Pxz_r[0]+0.3535533905932737*Pxz_c[0]; 

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
