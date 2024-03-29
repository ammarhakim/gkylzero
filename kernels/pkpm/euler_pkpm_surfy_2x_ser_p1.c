#include <gkyl_euler_pkpm_kernels.h> 
#include <gkyl_basis_ser_2x_p1_surfx2_eval_quad.h> 
#include <gkyl_basis_ser_2x_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double euler_pkpm_surfy_2x_ser_p1(const double *w, const double *dxv,
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

  const double dx1 = 2.0/dxv[1]; 

  const double *rhoux_l = &euler_pkpm_l[0]; 
  const double *rhouy_l = &euler_pkpm_l[4]; 
  const double *rhouz_l = &euler_pkpm_l[8]; 

  const double *rhoux_c = &euler_pkpm_c[0]; 
  const double *rhouy_c = &euler_pkpm_c[4]; 
  const double *rhouz_c = &euler_pkpm_c[8]; 

  const double *rhoux_r = &euler_pkpm_r[0]; 
  const double *rhouy_r = &euler_pkpm_r[4]; 
  const double *rhouz_r = &euler_pkpm_r[8]; 

  const double *rho_l = &vlasov_pkpm_moms_l[0]; 
  const double *rho_c = &vlasov_pkpm_moms_c[0]; 
  const double *rho_r = &vlasov_pkpm_moms_r[0]; 

  const double *ux_surf_lr = &prim_surf_l[18]; 
  const double *uy_surf_lr = &prim_surf_l[22]; 
  const double *uz_surf_lr = &prim_surf_l[26]; 

  const double *ux_surf_cl = &prim_surf_c[16]; 
  const double *uy_surf_cl = &prim_surf_c[20]; 
  const double *uz_surf_cl = &prim_surf_c[24]; 

  const double *ux_surf_cr = &prim_surf_c[18]; 
  const double *uy_surf_cr = &prim_surf_c[22]; 
  const double *uz_surf_cr = &prim_surf_c[26]; 

  const double *ux_surf_rl = &prim_surf_r[16]; 
  const double *uy_surf_rl = &prim_surf_r[20]; 
  const double *uz_surf_rl = &prim_surf_r[24]; 

  const double *P_surfx_lr = &p_ij_surf_l[14]; 
  const double *P_surfy_lr = &p_ij_surf_l[18]; 
  const double *P_surfz_lr = &p_ij_surf_l[22]; 

  const double *P_surfx_cl = &p_ij_surf_c[12]; 
  const double *P_surfy_cl = &p_ij_surf_c[16]; 
  const double *P_surfz_cl = &p_ij_surf_c[20]; 

  const double *P_surfx_cr = &p_ij_surf_c[14]; 
  const double *P_surfy_cr = &p_ij_surf_c[18]; 
  const double *P_surfz_cr = &p_ij_surf_c[22]; 

  const double *P_surfx_rl = &p_ij_surf_r[12]; 
  const double *P_surfy_rl = &p_ij_surf_r[16]; 
  const double *P_surfz_rl = &p_ij_surf_r[20]; 

  const double *pkpm_lax_l = &pkpm_lax[4]; 
  const double *pkpm_lax_r = &pkpm_lax[6]; 

  double *outrhou0 = &out[0]; 
  double *outrhou1 = &out[4]; 
  double *outrhou2 = &out[8]; 

  double flux_rho_l[2] = {0.0}; 
  double flux_rho_r[2] = {0.0}; 

  double avg_p_ij_x_l[2] = {0.0}; 
  double avg_p_ij_x_r[2] = {0.0}; 
  double avg_p_ij_y_l[2] = {0.0}; 
  double avg_p_ij_y_r[2] = {0.0}; 
  double avg_p_ij_z_l[2] = {0.0}; 
  double avg_p_ij_z_r[2] = {0.0}; 

  double jump_rhoux_l[2] = {0.0}; 
  double jump_rhoux_r[2] = {0.0}; 
  double jump_rhouy_l[2] = {0.0}; 
  double jump_rhouy_r[2] = {0.0}; 
  double jump_rhouz_l[2] = {0.0}; 
  double jump_rhouz_r[2] = {0.0}; 

  avg_p_ij_x_l[0] = 0.5*(P_surfx_lr[0] + P_surfx_cl[0]); 
  avg_p_ij_x_r[0] = 0.5*(P_surfx_cr[0] + P_surfx_rl[0]); 
  avg_p_ij_y_l[0] = 0.5*(P_surfy_lr[0] + P_surfy_cl[0]); 
  avg_p_ij_y_r[0] = 0.5*(P_surfy_cr[0] + P_surfy_rl[0]); 
  avg_p_ij_z_l[0] = 0.5*(P_surfz_lr[0] + P_surfz_cl[0]); 
  avg_p_ij_z_r[0] = 0.5*(P_surfz_cr[0] + P_surfz_rl[0]); 

  avg_p_ij_x_l[1] = 0.5*(P_surfx_lr[1] + P_surfx_cl[1]); 
  avg_p_ij_x_r[1] = 0.5*(P_surfx_cr[1] + P_surfx_rl[1]); 
  avg_p_ij_y_l[1] = 0.5*(P_surfy_lr[1] + P_surfy_cl[1]); 
  avg_p_ij_y_r[1] = 0.5*(P_surfy_cr[1] + P_surfy_rl[1]); 
  avg_p_ij_z_l[1] = 0.5*(P_surfz_lr[1] + P_surfz_cl[1]); 
  avg_p_ij_z_r[1] = 0.5*(P_surfz_cr[1] + P_surfz_rl[1]); 

  flux_rho_l[0] = 0.2165063509461096*uy_surf_lr[1]*rho_l[3]+0.2165063509461096*uy_surf_cl[1]*rho_l[3]-0.2165063509461096*uy_surf_lr[1]*rho_c[3]-0.2165063509461096*uy_surf_cl[1]*rho_c[3]+0.2165063509461096*uy_surf_lr[0]*rho_l[2]+0.2165063509461096*uy_surf_cl[0]*rho_l[2]-0.2165063509461096*uy_surf_lr[0]*rho_c[2]-0.2165063509461096*uy_surf_cl[0]*rho_c[2]+0.125*rho_l[1]*uy_surf_lr[1]+0.125*rho_c[1]*uy_surf_lr[1]+0.125*rho_l[1]*uy_surf_cl[1]+0.125*rho_c[1]*uy_surf_cl[1]+0.125*rho_l[0]*uy_surf_lr[0]+0.125*rho_c[0]*uy_surf_lr[0]+0.125*rho_l[0]*uy_surf_cl[0]+0.125*rho_c[0]*uy_surf_cl[0]; 
  flux_rho_l[1] = 0.2165063509461096*uy_surf_lr[0]*rho_l[3]+0.2165063509461096*uy_surf_cl[0]*rho_l[3]-0.2165063509461096*uy_surf_lr[0]*rho_c[3]-0.2165063509461096*uy_surf_cl[0]*rho_c[3]+0.2165063509461096*uy_surf_lr[1]*rho_l[2]+0.2165063509461096*uy_surf_cl[1]*rho_l[2]-0.2165063509461096*uy_surf_lr[1]*rho_c[2]-0.2165063509461096*uy_surf_cl[1]*rho_c[2]+0.125*rho_l[0]*uy_surf_lr[1]+0.125*rho_c[0]*uy_surf_lr[1]+0.125*rho_l[0]*uy_surf_cl[1]+0.125*rho_c[0]*uy_surf_cl[1]+0.125*uy_surf_lr[0]*rho_l[1]+0.125*uy_surf_cl[0]*rho_l[1]+0.125*uy_surf_lr[0]*rho_c[1]+0.125*uy_surf_cl[0]*rho_c[1]; 

  flux_rho_r[0] = (-0.2165063509461096*uy_surf_rl[1]*rho_r[3])-0.2165063509461096*uy_surf_cr[1]*rho_r[3]+0.2165063509461096*uy_surf_rl[1]*rho_c[3]+0.2165063509461096*uy_surf_cr[1]*rho_c[3]-0.2165063509461096*uy_surf_rl[0]*rho_r[2]-0.2165063509461096*uy_surf_cr[0]*rho_r[2]+0.2165063509461096*uy_surf_rl[0]*rho_c[2]+0.2165063509461096*uy_surf_cr[0]*rho_c[2]+0.125*rho_r[1]*uy_surf_rl[1]+0.125*rho_c[1]*uy_surf_rl[1]+0.125*rho_r[1]*uy_surf_cr[1]+0.125*rho_c[1]*uy_surf_cr[1]+0.125*rho_r[0]*uy_surf_rl[0]+0.125*rho_c[0]*uy_surf_rl[0]+0.125*rho_r[0]*uy_surf_cr[0]+0.125*rho_c[0]*uy_surf_cr[0]; 
  flux_rho_r[1] = (-0.2165063509461096*uy_surf_rl[0]*rho_r[3])-0.2165063509461096*uy_surf_cr[0]*rho_r[3]+0.2165063509461096*uy_surf_rl[0]*rho_c[3]+0.2165063509461096*uy_surf_cr[0]*rho_c[3]-0.2165063509461096*uy_surf_rl[1]*rho_r[2]-0.2165063509461096*uy_surf_cr[1]*rho_r[2]+0.2165063509461096*uy_surf_rl[1]*rho_c[2]+0.2165063509461096*uy_surf_cr[1]*rho_c[2]+0.125*rho_r[0]*uy_surf_rl[1]+0.125*rho_c[0]*uy_surf_rl[1]+0.125*rho_r[0]*uy_surf_cr[1]+0.125*rho_c[0]*uy_surf_cr[1]+0.125*uy_surf_rl[0]*rho_r[1]+0.125*uy_surf_cr[0]*rho_r[1]+0.125*uy_surf_rl[0]*rho_c[1]+0.125*uy_surf_cr[0]*rho_c[1]; 

  jump_rhoux_l[0] = (-0.6123724356957944*rhoux_l[2])-0.6123724356957944*rhoux_c[2]-0.3535533905932737*rhoux_l[0]+0.3535533905932737*rhoux_c[0]; 
  jump_rhoux_l[1] = (-0.6123724356957944*rhoux_l[3])-0.6123724356957944*rhoux_c[3]-0.3535533905932737*rhoux_l[1]+0.3535533905932737*rhoux_c[1]; 

  jump_rhoux_r[0] = (-0.6123724356957944*rhoux_r[2])-0.6123724356957944*rhoux_c[2]+0.3535533905932737*rhoux_r[0]-0.3535533905932737*rhoux_c[0]; 
  jump_rhoux_r[1] = (-0.6123724356957944*rhoux_r[3])-0.6123724356957944*rhoux_c[3]+0.3535533905932737*rhoux_r[1]-0.3535533905932737*rhoux_c[1]; 

  jump_rhouy_l[0] = (-0.6123724356957944*rhouy_l[2])-0.6123724356957944*rhouy_c[2]-0.3535533905932737*rhouy_l[0]+0.3535533905932737*rhouy_c[0]; 
  jump_rhouy_l[1] = (-0.6123724356957944*rhouy_l[3])-0.6123724356957944*rhouy_c[3]-0.3535533905932737*rhouy_l[1]+0.3535533905932737*rhouy_c[1]; 

  jump_rhouy_r[0] = (-0.6123724356957944*rhouy_r[2])-0.6123724356957944*rhouy_c[2]+0.3535533905932737*rhouy_r[0]-0.3535533905932737*rhouy_c[0]; 
  jump_rhouy_r[1] = (-0.6123724356957944*rhouy_r[3])-0.6123724356957944*rhouy_c[3]+0.3535533905932737*rhouy_r[1]-0.3535533905932737*rhouy_c[1]; 

  jump_rhouz_l[0] = (-0.6123724356957944*rhouz_l[2])-0.6123724356957944*rhouz_c[2]-0.3535533905932737*rhouz_l[0]+0.3535533905932737*rhouz_c[0]; 
  jump_rhouz_l[1] = (-0.6123724356957944*rhouz_l[3])-0.6123724356957944*rhouz_c[3]-0.3535533905932737*rhouz_l[1]+0.3535533905932737*rhouz_c[1]; 

  jump_rhouz_r[0] = (-0.6123724356957944*rhouz_r[2])-0.6123724356957944*rhouz_c[2]+0.3535533905932737*rhouz_r[0]-0.3535533905932737*rhouz_c[0]; 
  jump_rhouz_r[1] = (-0.6123724356957944*rhouz_r[3])-0.6123724356957944*rhouz_c[3]+0.3535533905932737*rhouz_r[1]-0.3535533905932737*rhouz_c[1]; 

  outrhou0[0] += ((-0.25*flux_rho_r[1]*ux_surf_rl[1])+0.25*flux_rho_l[1]*ux_surf_lr[1]-0.25*flux_rho_r[1]*ux_surf_cr[1]+0.25*flux_rho_l[1]*ux_surf_cl[1]+0.5*jump_rhoux_r[1]*pkpm_lax_r[1]-0.5*jump_rhoux_l[1]*pkpm_lax_l[1]-0.25*flux_rho_r[0]*ux_surf_rl[0]+0.25*flux_rho_l[0]*ux_surf_lr[0]-0.25*flux_rho_r[0]*ux_surf_cr[0]+0.25*flux_rho_l[0]*ux_surf_cl[0]+0.5*jump_rhoux_r[0]*pkpm_lax_r[0]-0.5*jump_rhoux_l[0]*pkpm_lax_l[0]-0.7071067811865475*avg_p_ij_x_r[0]+0.7071067811865475*avg_p_ij_x_l[0])*dx1; 
  outrhou0[1] += ((-0.25*flux_rho_r[0]*ux_surf_rl[1])+0.25*flux_rho_l[0]*ux_surf_lr[1]-0.25*flux_rho_r[0]*ux_surf_cr[1]+0.25*flux_rho_l[0]*ux_surf_cl[1]+0.5*jump_rhoux_r[0]*pkpm_lax_r[1]-0.5*jump_rhoux_l[0]*pkpm_lax_l[1]+0.5*pkpm_lax_r[0]*jump_rhoux_r[1]-0.5*pkpm_lax_l[0]*jump_rhoux_l[1]-0.25*(ux_surf_rl[0]+ux_surf_cr[0])*flux_rho_r[1]+0.25*(ux_surf_lr[0]+ux_surf_cl[0])*flux_rho_l[1]-0.7071067811865475*avg_p_ij_x_r[1]+0.7071067811865475*avg_p_ij_x_l[1])*dx1; 
  outrhou0[2] += ((-0.4330127018922193*(flux_rho_r[1]*ux_surf_rl[1]+flux_rho_l[1]*ux_surf_lr[1]+flux_rho_r[1]*ux_surf_cr[1]+flux_rho_l[1]*ux_surf_cl[1]))+0.8660254037844386*(jump_rhoux_r[1]*pkpm_lax_r[1]+jump_rhoux_l[1]*pkpm_lax_l[1])-0.4330127018922193*(flux_rho_r[0]*ux_surf_rl[0]+flux_rho_l[0]*ux_surf_lr[0]+flux_rho_r[0]*ux_surf_cr[0]+flux_rho_l[0]*ux_surf_cl[0])+0.8660254037844386*(jump_rhoux_r[0]*pkpm_lax_r[0]+jump_rhoux_l[0]*pkpm_lax_l[0])-1.224744871391589*(avg_p_ij_x_r[0]+avg_p_ij_x_l[0]))*dx1; 
  outrhou0[3] += ((-0.4330127018922193*(flux_rho_r[0]*ux_surf_rl[1]+flux_rho_l[0]*ux_surf_lr[1]+flux_rho_r[0]*ux_surf_cr[1]+flux_rho_l[0]*ux_surf_cl[1]))+0.8660254037844386*(jump_rhoux_r[0]*pkpm_lax_r[1]+jump_rhoux_l[0]*pkpm_lax_l[1]+pkpm_lax_r[0]*jump_rhoux_r[1]+pkpm_lax_l[0]*jump_rhoux_l[1])-0.4330127018922193*((ux_surf_rl[0]+ux_surf_cr[0])*flux_rho_r[1]+(ux_surf_lr[0]+ux_surf_cl[0])*flux_rho_l[1])-1.224744871391589*(avg_p_ij_x_r[1]+avg_p_ij_x_l[1]))*dx1; 

  outrhou1[0] += ((-0.25*flux_rho_r[1]*uy_surf_rl[1])+0.25*flux_rho_l[1]*uy_surf_lr[1]-0.25*flux_rho_r[1]*uy_surf_cr[1]+0.25*flux_rho_l[1]*uy_surf_cl[1]+0.5*jump_rhouy_r[1]*pkpm_lax_r[1]-0.5*jump_rhouy_l[1]*pkpm_lax_l[1]-0.25*flux_rho_r[0]*uy_surf_rl[0]+0.25*flux_rho_l[0]*uy_surf_lr[0]-0.25*flux_rho_r[0]*uy_surf_cr[0]+0.25*flux_rho_l[0]*uy_surf_cl[0]+0.5*jump_rhouy_r[0]*pkpm_lax_r[0]-0.5*jump_rhouy_l[0]*pkpm_lax_l[0]-0.7071067811865475*avg_p_ij_y_r[0]+0.7071067811865475*avg_p_ij_y_l[0])*dx1; 
  outrhou1[1] += ((-0.25*flux_rho_r[0]*uy_surf_rl[1])+0.25*flux_rho_l[0]*uy_surf_lr[1]-0.25*flux_rho_r[0]*uy_surf_cr[1]+0.25*flux_rho_l[0]*uy_surf_cl[1]+0.5*jump_rhouy_r[0]*pkpm_lax_r[1]-0.5*jump_rhouy_l[0]*pkpm_lax_l[1]+0.5*pkpm_lax_r[0]*jump_rhouy_r[1]-0.5*pkpm_lax_l[0]*jump_rhouy_l[1]-0.25*(uy_surf_rl[0]+uy_surf_cr[0])*flux_rho_r[1]+0.25*(uy_surf_lr[0]+uy_surf_cl[0])*flux_rho_l[1]-0.7071067811865475*avg_p_ij_y_r[1]+0.7071067811865475*avg_p_ij_y_l[1])*dx1; 
  outrhou1[2] += ((-0.4330127018922193*(flux_rho_r[1]*uy_surf_rl[1]+flux_rho_l[1]*uy_surf_lr[1]+flux_rho_r[1]*uy_surf_cr[1]+flux_rho_l[1]*uy_surf_cl[1]))+0.8660254037844386*(jump_rhouy_r[1]*pkpm_lax_r[1]+jump_rhouy_l[1]*pkpm_lax_l[1])-0.4330127018922193*(flux_rho_r[0]*uy_surf_rl[0]+flux_rho_l[0]*uy_surf_lr[0]+flux_rho_r[0]*uy_surf_cr[0]+flux_rho_l[0]*uy_surf_cl[0])+0.8660254037844386*(jump_rhouy_r[0]*pkpm_lax_r[0]+jump_rhouy_l[0]*pkpm_lax_l[0])-1.224744871391589*(avg_p_ij_y_r[0]+avg_p_ij_y_l[0]))*dx1; 
  outrhou1[3] += ((-0.4330127018922193*(flux_rho_r[0]*uy_surf_rl[1]+flux_rho_l[0]*uy_surf_lr[1]+flux_rho_r[0]*uy_surf_cr[1]+flux_rho_l[0]*uy_surf_cl[1]))+0.8660254037844386*(jump_rhouy_r[0]*pkpm_lax_r[1]+jump_rhouy_l[0]*pkpm_lax_l[1]+pkpm_lax_r[0]*jump_rhouy_r[1]+pkpm_lax_l[0]*jump_rhouy_l[1])-0.4330127018922193*((uy_surf_rl[0]+uy_surf_cr[0])*flux_rho_r[1]+(uy_surf_lr[0]+uy_surf_cl[0])*flux_rho_l[1])-1.224744871391589*(avg_p_ij_y_r[1]+avg_p_ij_y_l[1]))*dx1; 

  outrhou2[0] += ((-0.25*flux_rho_r[1]*uz_surf_rl[1])+0.25*flux_rho_l[1]*uz_surf_lr[1]-0.25*flux_rho_r[1]*uz_surf_cr[1]+0.25*flux_rho_l[1]*uz_surf_cl[1]+0.5*jump_rhouz_r[1]*pkpm_lax_r[1]-0.5*jump_rhouz_l[1]*pkpm_lax_l[1]-0.25*flux_rho_r[0]*uz_surf_rl[0]+0.25*flux_rho_l[0]*uz_surf_lr[0]-0.25*flux_rho_r[0]*uz_surf_cr[0]+0.25*flux_rho_l[0]*uz_surf_cl[0]+0.5*jump_rhouz_r[0]*pkpm_lax_r[0]-0.5*jump_rhouz_l[0]*pkpm_lax_l[0]-0.7071067811865475*avg_p_ij_z_r[0]+0.7071067811865475*avg_p_ij_z_l[0])*dx1; 
  outrhou2[1] += ((-0.25*flux_rho_r[0]*uz_surf_rl[1])+0.25*flux_rho_l[0]*uz_surf_lr[1]-0.25*flux_rho_r[0]*uz_surf_cr[1]+0.25*flux_rho_l[0]*uz_surf_cl[1]+0.5*jump_rhouz_r[0]*pkpm_lax_r[1]-0.5*jump_rhouz_l[0]*pkpm_lax_l[1]+0.5*pkpm_lax_r[0]*jump_rhouz_r[1]-0.5*pkpm_lax_l[0]*jump_rhouz_l[1]-0.25*(uz_surf_rl[0]+uz_surf_cr[0])*flux_rho_r[1]+0.25*(uz_surf_lr[0]+uz_surf_cl[0])*flux_rho_l[1]-0.7071067811865475*avg_p_ij_z_r[1]+0.7071067811865475*avg_p_ij_z_l[1])*dx1; 
  outrhou2[2] += ((-0.4330127018922193*(flux_rho_r[1]*uz_surf_rl[1]+flux_rho_l[1]*uz_surf_lr[1]+flux_rho_r[1]*uz_surf_cr[1]+flux_rho_l[1]*uz_surf_cl[1]))+0.8660254037844386*(jump_rhouz_r[1]*pkpm_lax_r[1]+jump_rhouz_l[1]*pkpm_lax_l[1])-0.4330127018922193*(flux_rho_r[0]*uz_surf_rl[0]+flux_rho_l[0]*uz_surf_lr[0]+flux_rho_r[0]*uz_surf_cr[0]+flux_rho_l[0]*uz_surf_cl[0])+0.8660254037844386*(jump_rhouz_r[0]*pkpm_lax_r[0]+jump_rhouz_l[0]*pkpm_lax_l[0])-1.224744871391589*(avg_p_ij_z_r[0]+avg_p_ij_z_l[0]))*dx1; 
  outrhou2[3] += ((-0.4330127018922193*(flux_rho_r[0]*uz_surf_rl[1]+flux_rho_l[0]*uz_surf_lr[1]+flux_rho_r[0]*uz_surf_cr[1]+flux_rho_l[0]*uz_surf_cl[1]))+0.8660254037844386*(jump_rhouz_r[0]*pkpm_lax_r[1]+jump_rhouz_l[0]*pkpm_lax_l[1]+pkpm_lax_r[0]*jump_rhouz_r[1]+pkpm_lax_l[0]*jump_rhouz_l[1])-0.4330127018922193*((uz_surf_rl[0]+uz_surf_cr[0])*flux_rho_r[1]+(uz_surf_lr[0]+uz_surf_cl[0])*flux_rho_l[1])-1.224744871391589*(avg_p_ij_z_r[1]+avg_p_ij_z_l[1]))*dx1; 

  return 0.;

} 
