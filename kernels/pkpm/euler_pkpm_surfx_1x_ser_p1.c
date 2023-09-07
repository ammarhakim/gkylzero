#include <gkyl_euler_pkpm_kernels.h> 
#include <gkyl_basis_ser_1x_p1_surfx1_eval_quad.h> 
GKYL_CU_DH double euler_pkpm_surfx_1x_ser_p1(const double *w, const double *dxv, 
  const double *vlasov_pkpm_momsl, const double *vlasov_pkpm_momsc, const double *vlasov_pkpm_momsr,
  const double *u_il, const double *u_ic, const double *u_ir,
  const double *T_ijl, const double *T_ijc, const double *T_ijr,
  const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:                       Cell-center coordinates.
  // dxv[NDIM]:                     Cell spacing.
  // vlasov_pkpm_momsl/vlasov_pkpm_momsc/vlasov_pkpm_momsr: Input pkpm moments in left/center/right cells.
  // u_il/u_ic/u_ir:                Input bulk velocity (ux,uy,uz) in left/center/right cells.
  // T_ijl/T_ijc/T_ijr:             Input Temperature tensor/mass (for penalization) in left/center/right cells.
  // statevecl/statevecc/statevecr: [rho ux, rho uy, rho uz], Fluid input state vector in left/center/right cells.
  // out: Incremented output.

  const double dx1 = 2.0/dxv[0]; 

  const double *rhoux_l = &statevecl[0]; 
  const double *rhouy_l = &statevecl[2]; 
  const double *rhouz_l = &statevecl[4]; 

  const double *rhoux_c = &statevecc[0]; 
  const double *rhouy_c = &statevecc[2]; 
  const double *rhouz_c = &statevecc[4]; 

  const double *rhoux_r = &statevecr[0]; 
  const double *rhouy_r = &statevecr[2]; 
  const double *rhouz_r = &statevecr[4]; 

  const double *rho_l = &vlasov_pkpm_momsl[0]; 
  const double *rho_c = &vlasov_pkpm_momsc[0]; 
  const double *rho_r = &vlasov_pkpm_momsr[0]; 

  const double *ux_l = &u_il[0]; 
  const double *ux_c = &u_ic[0]; 
  const double *ux_r = &u_ir[0]; 

  const double *uy_l = &u_il[2]; 
  const double *uy_c = &u_ic[2]; 
  const double *uy_r = &u_ir[2]; 

  const double *uz_l = &u_il[4]; 
  const double *uz_c = &u_ic[4]; 
  const double *uz_r = &u_ir[4]; 

  const double *vth_sql = &T_ijl[0]; 
  const double *vth_sqc = &T_ijc[0]; 
  const double *vth_sqr = &T_ijr[0]; 

  double *outrhou0 = &out[0]; 
  double *outrhou1 = &out[2]; 
  double *outrhou2 = &out[4]; 

  double Ghat_rhoux_l = 0.0; 
  double Ghat_rhoux_r = 0.0; 
  double Ghat_rhouy_l = 0.0; 
  double Ghat_rhouy_r = 0.0; 
  double Ghat_rhouz_l = 0.0; 
  double Ghat_rhouz_r = 0.0; 
  double uxl_r = ser_1x_p1_surfx1_eval_quad_node_0_r(ux_l); 
  double uxc_l = ser_1x_p1_surfx1_eval_quad_node_0_l(ux_c); 
  double uxc_r = ser_1x_p1_surfx1_eval_quad_node_0_r(ux_c); 
  double uxr_l = ser_1x_p1_surfx1_eval_quad_node_0_l(ux_r); 

  double uxl_r_sq = uxl_r*uxl_r; 
  double uxc_l_sq = uxc_l*uxc_l; 
  double uxc_r_sq = uxc_r*uxc_r; 
  double uxr_l_sq = uxr_l*uxr_l; 

  double uyl_r = ser_1x_p1_surfx1_eval_quad_node_0_r(uy_l); 
  double uyc_l = ser_1x_p1_surfx1_eval_quad_node_0_l(uy_c); 
  double uyc_r = ser_1x_p1_surfx1_eval_quad_node_0_r(uy_c); 
  double uyr_l = ser_1x_p1_surfx1_eval_quad_node_0_l(uy_r); 

  double uzl_r = ser_1x_p1_surfx1_eval_quad_node_0_r(uz_l); 
  double uzc_l = ser_1x_p1_surfx1_eval_quad_node_0_l(uz_c); 
  double uzc_r = ser_1x_p1_surfx1_eval_quad_node_0_r(uz_c); 
  double uzr_l = ser_1x_p1_surfx1_eval_quad_node_0_l(uz_r); 

  double vth_sq_l_r = ser_1x_p1_surfx1_eval_quad_node_0_r(vth_sql); 
  double vth_sq_c_l = ser_1x_p1_surfx1_eval_quad_node_0_l(vth_sqc); 
  double vth_sq_c_r = ser_1x_p1_surfx1_eval_quad_node_0_r(vth_sqc); 
  double vth_sq_r_l = ser_1x_p1_surfx1_eval_quad_node_0_l(vth_sqr); 

  double ux_max_l = fmax(fabs(uxl_r), fabs(uxc_l)); 
  double ux_max_r = fmax(fabs(uxc_r), fabs(uxr_l)); 
  double vth_max_l = fmax(sqrt(fabs(vth_sq_l_r)), sqrt(fabs(vth_sq_c_l))); 
  double vth_max_r = fmax(sqrt(fabs(vth_sq_c_r)), sqrt(fabs(vth_sq_r_l))); 
  double max_speed_l = ux_max_l + vth_max_l; 
  double max_speed_r = ux_max_r + vth_max_r; 

  Ghat_rhoux_l = 0.1530931089239486*rho_l[1]*uxl_r_sq-0.1530931089239486*rho_c[1]*uxl_r_sq+0.0883883476483184*rho_l[0]*uxl_r_sq+0.0883883476483184*rho_c[0]*uxl_r_sq+0.3061862178478971*rho_l[1]*uxc_l*uxl_r-0.3061862178478971*rho_c[1]*uxc_l*uxl_r+0.1767766952966368*rho_l[0]*uxc_l*uxl_r+0.1767766952966368*rho_c[0]*uxc_l*uxl_r+0.1530931089239486*rho_l[1]*uxc_l_sq-0.1530931089239486*rho_c[1]*uxc_l_sq+0.0883883476483184*rho_l[0]*uxc_l_sq+0.0883883476483184*rho_c[0]*uxc_l_sq+0.6123724356957944*rhoux_l[1]*max_speed_l+0.6123724356957944*rhoux_c[1]*max_speed_l+0.3535533905932737*rhoux_l[0]*max_speed_l-0.3535533905932737*rhoux_c[0]*max_speed_l; 
  Ghat_rhouy_l = 0.1530931089239486*rho_l[1]*uxl_r*uyl_r-0.1530931089239486*rho_c[1]*uxl_r*uyl_r+0.08838834764831843*rho_l[0]*uxl_r*uyl_r+0.08838834764831843*rho_c[0]*uxl_r*uyl_r+0.1530931089239486*rho_l[1]*uxc_l*uyl_r-0.1530931089239486*rho_c[1]*uxc_l*uyl_r+0.08838834764831843*rho_l[0]*uxc_l*uyl_r+0.08838834764831843*rho_c[0]*uxc_l*uyl_r+0.1530931089239486*rho_l[1]*uxl_r*uyc_l-0.1530931089239486*rho_c[1]*uxl_r*uyc_l+0.08838834764831843*rho_l[0]*uxl_r*uyc_l+0.08838834764831843*rho_c[0]*uxl_r*uyc_l+0.1530931089239486*rho_l[1]*uxc_l*uyc_l-0.1530931089239486*rho_c[1]*uxc_l*uyc_l+0.08838834764831843*rho_l[0]*uxc_l*uyc_l+0.08838834764831843*rho_c[0]*uxc_l*uyc_l+0.6123724356957945*rhouy_l[1]*max_speed_l+0.6123724356957945*rhouy_c[1]*max_speed_l+0.3535533905932737*rhouy_l[0]*max_speed_l-0.3535533905932737*rhouy_c[0]*max_speed_l; 
  Ghat_rhouz_l = 0.1530931089239486*rho_l[1]*uxl_r*uzl_r-0.1530931089239486*rho_c[1]*uxl_r*uzl_r+0.08838834764831843*rho_l[0]*uxl_r*uzl_r+0.08838834764831843*rho_c[0]*uxl_r*uzl_r+0.1530931089239486*rho_l[1]*uxc_l*uzl_r-0.1530931089239486*rho_c[1]*uxc_l*uzl_r+0.08838834764831843*rho_l[0]*uxc_l*uzl_r+0.08838834764831843*rho_c[0]*uxc_l*uzl_r+0.1530931089239486*rho_l[1]*uxl_r*uzc_l-0.1530931089239486*rho_c[1]*uxl_r*uzc_l+0.08838834764831843*rho_l[0]*uxl_r*uzc_l+0.08838834764831843*rho_c[0]*uxl_r*uzc_l+0.1530931089239486*rho_l[1]*uxc_l*uzc_l-0.1530931089239486*rho_c[1]*uxc_l*uzc_l+0.08838834764831843*rho_l[0]*uxc_l*uzc_l+0.08838834764831843*rho_c[0]*uxc_l*uzc_l+0.6123724356957945*rhouz_l[1]*max_speed_l+0.6123724356957945*rhouz_c[1]*max_speed_l+0.3535533905932737*rhouz_l[0]*max_speed_l-0.3535533905932737*rhouz_c[0]*max_speed_l; 
  Ghat_rhoux_r = (-0.1530931089239486*rho_r[1]*uxr_l_sq)+0.1530931089239486*rho_c[1]*uxr_l_sq+0.0883883476483184*rho_r[0]*uxr_l_sq+0.0883883476483184*rho_c[0]*uxr_l_sq-0.3061862178478971*rho_r[1]*uxc_r*uxr_l+0.3061862178478971*rho_c[1]*uxc_r*uxr_l+0.1767766952966368*rho_r[0]*uxc_r*uxr_l+0.1767766952966368*rho_c[0]*uxc_r*uxr_l-0.1530931089239486*rho_r[1]*uxc_r_sq+0.1530931089239486*rho_c[1]*uxc_r_sq+0.0883883476483184*rho_r[0]*uxc_r_sq+0.0883883476483184*rho_c[0]*uxc_r_sq+0.6123724356957944*rhoux_r[1]*max_speed_r+0.6123724356957944*rhoux_c[1]*max_speed_r-0.3535533905932737*rhoux_r[0]*max_speed_r+0.3535533905932737*rhoux_c[0]*max_speed_r; 
  Ghat_rhouy_r = (-0.1530931089239486*rho_r[1]*uxr_l*uyr_l)+0.1530931089239486*rho_c[1]*uxr_l*uyr_l+0.08838834764831843*rho_r[0]*uxr_l*uyr_l+0.08838834764831843*rho_c[0]*uxr_l*uyr_l-0.1530931089239486*rho_r[1]*uxc_r*uyr_l+0.1530931089239486*rho_c[1]*uxc_r*uyr_l+0.08838834764831843*rho_r[0]*uxc_r*uyr_l+0.08838834764831843*rho_c[0]*uxc_r*uyr_l-0.1530931089239486*rho_r[1]*uxr_l*uyc_r+0.1530931089239486*rho_c[1]*uxr_l*uyc_r+0.08838834764831843*rho_r[0]*uxr_l*uyc_r+0.08838834764831843*rho_c[0]*uxr_l*uyc_r-0.1530931089239486*rho_r[1]*uxc_r*uyc_r+0.1530931089239486*rho_c[1]*uxc_r*uyc_r+0.08838834764831843*rho_r[0]*uxc_r*uyc_r+0.08838834764831843*rho_c[0]*uxc_r*uyc_r+0.6123724356957945*rhouy_r[1]*max_speed_r+0.6123724356957945*rhouy_c[1]*max_speed_r-0.3535533905932737*rhouy_r[0]*max_speed_r+0.3535533905932737*rhouy_c[0]*max_speed_r; 
  Ghat_rhouz_r = (-0.1530931089239486*rho_r[1]*uxr_l*uzr_l)+0.1530931089239486*rho_c[1]*uxr_l*uzr_l+0.08838834764831843*rho_r[0]*uxr_l*uzr_l+0.08838834764831843*rho_c[0]*uxr_l*uzr_l-0.1530931089239486*rho_r[1]*uxc_r*uzr_l+0.1530931089239486*rho_c[1]*uxc_r*uzr_l+0.08838834764831843*rho_r[0]*uxc_r*uzr_l+0.08838834764831843*rho_c[0]*uxc_r*uzr_l-0.1530931089239486*rho_r[1]*uxr_l*uzc_r+0.1530931089239486*rho_c[1]*uxr_l*uzc_r+0.08838834764831843*rho_r[0]*uxr_l*uzc_r+0.08838834764831843*rho_c[0]*uxr_l*uzc_r-0.1530931089239486*rho_r[1]*uxc_r*uzc_r+0.1530931089239486*rho_c[1]*uxc_r*uzc_r+0.08838834764831843*rho_r[0]*uxc_r*uzc_r+0.08838834764831843*rho_c[0]*uxc_r*uzc_r+0.6123724356957945*rhouz_r[1]*max_speed_r+0.6123724356957945*rhouz_c[1]*max_speed_r-0.3535533905932737*rhouz_r[0]*max_speed_r+0.3535533905932737*rhouz_c[0]*max_speed_r; 

  outrhou0[0] += 0.7071067811865475*Ghat_rhoux_l*dx1-0.7071067811865475*Ghat_rhoux_r*dx1; 
  outrhou0[1] += (-1.224744871391589*Ghat_rhoux_r*dx1)-1.224744871391589*Ghat_rhoux_l*dx1; 

  outrhou1[0] += 0.7071067811865475*Ghat_rhouy_l*dx1-0.7071067811865475*Ghat_rhouy_r*dx1; 
  outrhou1[1] += (-1.224744871391589*Ghat_rhouy_r*dx1)-1.224744871391589*Ghat_rhouy_l*dx1; 

  outrhou2[0] += 0.7071067811865475*Ghat_rhouz_l*dx1-0.7071067811865475*Ghat_rhouz_r*dx1; 
  outrhou2[1] += (-1.224744871391589*Ghat_rhouz_r*dx1)-1.224744871391589*Ghat_rhouz_l*dx1; 

  return 0.;

} 
