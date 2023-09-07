#include <gkyl_euler_pkpm_kernels.h> 
#include <gkyl_basis_ser_2x_p1_surfx2_eval_quad.h> 
#include <gkyl_basis_ser_2x_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double euler_pkpm_surfy_2x_ser_p1(const double *w, const double *dxv,
  const double *vlasov_pkpm_momsl, const double *vlasov_pkpm_momsc, const double *vlasov_pkpm_momsr,
  const double *u_il, const double *u_ic, const double *u_ir,
  const double *T_ijl, const double *T_ijc, const double *T_ijr,
  const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:                       Cell-center coordinates.
  // dxv[NDIM]:                     Cell spacing.
  // u_il/u_ic/u_ir:                Input bulk velocity (ux,uy,uz) in left/center/right cells.
  // T_ijl/T_ijc/T_ijr:             Input Temperature tensor/mass (for penalization) in left/center/right cells.
  // statevecl/statevecc/statevecr: [rho ux, rho uy, rho uz, E_perp], Fluid input state vector in left/center/right cells.
  // out: Incremented output.

  const double dx1 = 2.0/dxv[1]; 

  const double *rhoux_l = &statevecl[0]; 
  const double *rhouy_l = &statevecl[4]; 
  const double *rhouz_l = &statevecl[8]; 

  const double *rhoux_c = &statevecc[0]; 
  const double *rhouy_c = &statevecc[4]; 
  const double *rhouz_c = &statevecc[8]; 

  const double *rhoux_r = &statevecr[0]; 
  const double *rhouy_r = &statevecr[4]; 
  const double *rhouz_r = &statevecr[8]; 

  const double *rho_l = &vlasov_pkpm_momsl[0]; 
  const double *rho_c = &vlasov_pkpm_momsc[0]; 
  const double *rho_r = &vlasov_pkpm_momsr[0]; 

  const double *ux_l = &u_il[0]; 
  const double *ux_c = &u_ic[0]; 
  const double *ux_r = &u_ir[0]; 

  const double *uy_l = &u_il[4]; 
  const double *uy_c = &u_ic[4]; 
  const double *uy_r = &u_ir[4]; 

  const double *uz_l = &u_il[8]; 
  const double *uz_c = &u_ic[8]; 
  const double *uz_r = &u_ir[8]; 

  // Get another pointer for u in direction of update for ease of flux calculation. 
  const double *u_l = &u_il[4]; 
  const double *u_c = &u_ic[4]; 
  const double *u_r = &u_ir[4]; 

  // Get thermal velocity in direction of update for penalization vth^2 = 3.0*T_ii/m. 
  const double *vth_sql = &T_ijl[12]; 
  const double *vth_sqc = &T_ijc[12]; 
  const double *vth_sqr = &T_ijr[12]; 

  double *outrhou0 = &out[0]; 
  double *outrhou1 = &out[4]; 
  double *outrhou2 = &out[8]; 

  double flux_rho_l[2] = {0.0}; 
  double flux_rho_r[2] = {0.0}; 
  double avg_ux_l[2] = {0.0}; 
  double avg_ux_r[2] = {0.0}; 
  double avg_uy_l[2] = {0.0}; 
  double avg_uy_r[2] = {0.0}; 
  double avg_uz_l[2] = {0.0}; 
  double avg_uz_r[2] = {0.0}; 
  double jump_rhoux_l[2] = {0.0}; 
  double jump_rhoux_r[2] = {0.0}; 
  double jump_rhouy_l[2] = {0.0}; 
  double jump_rhouy_r[2] = {0.0}; 
  double jump_rhouz_l[2] = {0.0}; 
  double jump_rhouz_r[2] = {0.0}; 
  double max_speed_quad_l[2] = {0.0}; 
  double max_speed_quad_r[2] = {0.0}; 
  double max_speed_modal_l[2] = {0.0}; 
  double max_speed_modal_r[2] = {0.0}; 
  flux_rho_l[0] = 0.2651650429449552*rho_l[3]*u_l[3]-0.2651650429449552*rho_c[3]*u_l[3]+0.1530931089239486*rho_l[1]*u_l[3]+0.1530931089239486*rho_c[1]*u_l[3]-0.2651650429449552*rho_l[3]*u_c[3]+0.2651650429449552*rho_c[3]*u_c[3]-0.1530931089239486*rho_l[1]*u_c[3]-0.1530931089239486*rho_c[1]*u_c[3]+0.1530931089239486*u_l[1]*rho_l[3]+0.1530931089239486*u_c[1]*rho_l[3]-0.1530931089239486*u_l[1]*rho_c[3]-0.1530931089239486*u_c[1]*rho_c[3]+0.2651650429449552*rho_l[2]*u_l[2]-0.2651650429449552*rho_c[2]*u_l[2]+0.1530931089239486*rho_l[0]*u_l[2]+0.1530931089239486*rho_c[0]*u_l[2]-0.2651650429449552*rho_l[2]*u_c[2]+0.2651650429449552*rho_c[2]*u_c[2]-0.1530931089239486*rho_l[0]*u_c[2]-0.1530931089239486*rho_c[0]*u_c[2]+0.1530931089239486*u_l[0]*rho_l[2]+0.1530931089239486*u_c[0]*rho_l[2]-0.1530931089239486*u_l[0]*rho_c[2]-0.1530931089239486*u_c[0]*rho_c[2]+0.0883883476483184*rho_l[1]*u_l[1]+0.0883883476483184*rho_c[1]*u_l[1]+0.0883883476483184*rho_l[1]*u_c[1]+0.0883883476483184*rho_c[1]*u_c[1]+0.0883883476483184*rho_l[0]*u_l[0]+0.0883883476483184*rho_c[0]*u_l[0]+0.0883883476483184*rho_l[0]*u_c[0]+0.0883883476483184*rho_c[0]*u_c[0]; 
  flux_rho_l[1] = 0.2651650429449552*rho_l[2]*u_l[3]-0.2651650429449552*rho_c[2]*u_l[3]+0.1530931089239486*rho_l[0]*u_l[3]+0.1530931089239486*rho_c[0]*u_l[3]-0.2651650429449552*rho_l[2]*u_c[3]+0.2651650429449552*rho_c[2]*u_c[3]-0.1530931089239486*rho_l[0]*u_c[3]-0.1530931089239486*rho_c[0]*u_c[3]+0.2651650429449552*u_l[2]*rho_l[3]-0.2651650429449552*u_c[2]*rho_l[3]+0.1530931089239486*u_l[0]*rho_l[3]+0.1530931089239486*u_c[0]*rho_l[3]-0.2651650429449552*u_l[2]*rho_c[3]+0.2651650429449552*u_c[2]*rho_c[3]-0.1530931089239486*u_l[0]*rho_c[3]-0.1530931089239486*u_c[0]*rho_c[3]+0.1530931089239486*rho_l[1]*u_l[2]+0.1530931089239486*rho_c[1]*u_l[2]-0.1530931089239486*rho_l[1]*u_c[2]-0.1530931089239486*rho_c[1]*u_c[2]+0.1530931089239486*u_l[1]*rho_l[2]+0.1530931089239486*u_c[1]*rho_l[2]-0.1530931089239486*u_l[1]*rho_c[2]-0.1530931089239486*u_c[1]*rho_c[2]+0.0883883476483184*rho_l[0]*u_l[1]+0.0883883476483184*rho_c[0]*u_l[1]+0.0883883476483184*rho_l[0]*u_c[1]+0.0883883476483184*rho_c[0]*u_c[1]+0.0883883476483184*u_l[0]*rho_l[1]+0.0883883476483184*u_c[0]*rho_l[1]+0.0883883476483184*u_l[0]*rho_c[1]+0.0883883476483184*u_c[0]*rho_c[1]; 

  flux_rho_r[0] = 0.2651650429449552*rho_r[3]*u_r[3]-0.2651650429449552*rho_c[3]*u_r[3]-0.1530931089239486*rho_r[1]*u_r[3]-0.1530931089239486*rho_c[1]*u_r[3]-0.2651650429449552*rho_r[3]*u_c[3]+0.2651650429449552*rho_c[3]*u_c[3]+0.1530931089239486*rho_r[1]*u_c[3]+0.1530931089239486*rho_c[1]*u_c[3]-0.1530931089239486*u_r[1]*rho_r[3]-0.1530931089239486*u_c[1]*rho_r[3]+0.1530931089239486*u_r[1]*rho_c[3]+0.1530931089239486*u_c[1]*rho_c[3]+0.2651650429449552*rho_r[2]*u_r[2]-0.2651650429449552*rho_c[2]*u_r[2]-0.1530931089239486*rho_r[0]*u_r[2]-0.1530931089239486*rho_c[0]*u_r[2]-0.2651650429449552*rho_r[2]*u_c[2]+0.2651650429449552*rho_c[2]*u_c[2]+0.1530931089239486*rho_r[0]*u_c[2]+0.1530931089239486*rho_c[0]*u_c[2]-0.1530931089239486*u_r[0]*rho_r[2]-0.1530931089239486*u_c[0]*rho_r[2]+0.1530931089239486*u_r[0]*rho_c[2]+0.1530931089239486*u_c[0]*rho_c[2]+0.0883883476483184*rho_r[1]*u_r[1]+0.0883883476483184*rho_c[1]*u_r[1]+0.0883883476483184*rho_r[1]*u_c[1]+0.0883883476483184*rho_c[1]*u_c[1]+0.0883883476483184*rho_r[0]*u_r[0]+0.0883883476483184*rho_c[0]*u_r[0]+0.0883883476483184*rho_r[0]*u_c[0]+0.0883883476483184*rho_c[0]*u_c[0]; 
  flux_rho_r[1] = 0.2651650429449552*rho_r[2]*u_r[3]-0.2651650429449552*rho_c[2]*u_r[3]-0.1530931089239486*rho_r[0]*u_r[3]-0.1530931089239486*rho_c[0]*u_r[3]-0.2651650429449552*rho_r[2]*u_c[3]+0.2651650429449552*rho_c[2]*u_c[3]+0.1530931089239486*rho_r[0]*u_c[3]+0.1530931089239486*rho_c[0]*u_c[3]+0.2651650429449552*u_r[2]*rho_r[3]-0.2651650429449552*u_c[2]*rho_r[3]-0.1530931089239486*u_r[0]*rho_r[3]-0.1530931089239486*u_c[0]*rho_r[3]-0.2651650429449552*u_r[2]*rho_c[3]+0.2651650429449552*u_c[2]*rho_c[3]+0.1530931089239486*u_r[0]*rho_c[3]+0.1530931089239486*u_c[0]*rho_c[3]-0.1530931089239486*rho_r[1]*u_r[2]-0.1530931089239486*rho_c[1]*u_r[2]+0.1530931089239486*rho_r[1]*u_c[2]+0.1530931089239486*rho_c[1]*u_c[2]-0.1530931089239486*u_r[1]*rho_r[2]-0.1530931089239486*u_c[1]*rho_r[2]+0.1530931089239486*u_r[1]*rho_c[2]+0.1530931089239486*u_c[1]*rho_c[2]+0.0883883476483184*rho_r[0]*u_r[1]+0.0883883476483184*rho_c[0]*u_r[1]+0.0883883476483184*rho_r[0]*u_c[1]+0.0883883476483184*rho_c[0]*u_c[1]+0.0883883476483184*u_r[0]*rho_r[1]+0.0883883476483184*u_c[0]*rho_r[1]+0.0883883476483184*u_r[0]*rho_c[1]+0.0883883476483184*u_c[0]*rho_c[1]; 

  avg_ux_l[0] = 0.6123724356957944*ux_l[2]-0.6123724356957944*ux_c[2]+0.3535533905932737*ux_l[0]+0.3535533905932737*ux_c[0]; 
  avg_ux_l[1] = 0.6123724356957944*ux_l[3]-0.6123724356957944*ux_c[3]+0.3535533905932737*ux_l[1]+0.3535533905932737*ux_c[1]; 

  avg_ux_r[0] = (-0.6123724356957944*ux_r[2])+0.6123724356957944*ux_c[2]+0.3535533905932737*ux_r[0]+0.3535533905932737*ux_c[0]; 
  avg_ux_r[1] = (-0.6123724356957944*ux_r[3])+0.6123724356957944*ux_c[3]+0.3535533905932737*ux_r[1]+0.3535533905932737*ux_c[1]; 

  avg_uy_l[0] = 0.6123724356957944*uy_l[2]-0.6123724356957944*uy_c[2]+0.3535533905932737*uy_l[0]+0.3535533905932737*uy_c[0]; 
  avg_uy_l[1] = 0.6123724356957944*uy_l[3]-0.6123724356957944*uy_c[3]+0.3535533905932737*uy_l[1]+0.3535533905932737*uy_c[1]; 

  avg_uy_r[0] = (-0.6123724356957944*uy_r[2])+0.6123724356957944*uy_c[2]+0.3535533905932737*uy_r[0]+0.3535533905932737*uy_c[0]; 
  avg_uy_r[1] = (-0.6123724356957944*uy_r[3])+0.6123724356957944*uy_c[3]+0.3535533905932737*uy_r[1]+0.3535533905932737*uy_c[1]; 

  avg_uz_l[0] = 0.6123724356957944*uz_l[2]-0.6123724356957944*uz_c[2]+0.3535533905932737*uz_l[0]+0.3535533905932737*uz_c[0]; 
  avg_uz_l[1] = 0.6123724356957944*uz_l[3]-0.6123724356957944*uz_c[3]+0.3535533905932737*uz_l[1]+0.3535533905932737*uz_c[1]; 

  avg_uz_r[0] = (-0.6123724356957944*uz_r[2])+0.6123724356957944*uz_c[2]+0.3535533905932737*uz_r[0]+0.3535533905932737*uz_c[0]; 
  avg_uz_r[1] = (-0.6123724356957944*uz_r[3])+0.6123724356957944*uz_c[3]+0.3535533905932737*uz_r[1]+0.3535533905932737*uz_c[1]; 

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

  double ul_r = 0.0; 
  double uc_l = 0.0; 
  double uc_r = 0.0; 
  double ur_l = 0.0; 
  double uQuad_l = 0.0; 
  double uQuad_r = 0.0; 

  double vth_sql_r = 0.0; 
  double vth_sqc_l = 0.0; 
  double vth_sqc_r = 0.0; 
  double vth_sqr_l = 0.0; 
  double vthQuad_l = 0.0; 
  double vthQuad_r = 0.0; 

  ul_r = ser_2x_p1_surfx2_eval_quad_node_0_r(u_l); 
  uc_l = ser_2x_p1_surfx2_eval_quad_node_0_l(u_c); 
  uc_r = ser_2x_p1_surfx2_eval_quad_node_0_r(u_c); 
  ur_l = ser_2x_p1_surfx2_eval_quad_node_0_l(u_r); 
  uQuad_l = fmax(fabs(ul_r), fabs(uc_l)); 
  uQuad_r = fmax(fabs(uc_r), fabs(ur_l)); 
  vth_sql_r = ser_2x_p1_surfx2_eval_quad_node_0_r(vth_sql); 
  vth_sqc_l = ser_2x_p1_surfx2_eval_quad_node_0_l(vth_sqc); 
  vth_sqc_r = ser_2x_p1_surfx2_eval_quad_node_0_r(vth_sqc); 
  vth_sqr_l = ser_2x_p1_surfx2_eval_quad_node_0_l(vth_sqr); 
  vthQuad_l = fmax(sqrt(fabs(vth_sql_r)), sqrt(fabs(vth_sqc_l))); 
  vthQuad_r = fmax(sqrt(fabs(vth_sqc_r)), sqrt(fabs(vth_sqr_l))); 
  max_speed_quad_l[0] = uQuad_l + vthQuad_l; 
  max_speed_quad_r[0] = uQuad_r + vthQuad_r; 
  ul_r = ser_2x_p1_surfx2_eval_quad_node_1_r(u_l); 
  uc_l = ser_2x_p1_surfx2_eval_quad_node_1_l(u_c); 
  uc_r = ser_2x_p1_surfx2_eval_quad_node_1_r(u_c); 
  ur_l = ser_2x_p1_surfx2_eval_quad_node_1_l(u_r); 
  uQuad_l = fmax(fabs(ul_r), fabs(uc_l)); 
  uQuad_r = fmax(fabs(uc_r), fabs(ur_l)); 
  vth_sql_r = ser_2x_p1_surfx2_eval_quad_node_1_r(vth_sql); 
  vth_sqc_l = ser_2x_p1_surfx2_eval_quad_node_1_l(vth_sqc); 
  vth_sqc_r = ser_2x_p1_surfx2_eval_quad_node_1_r(vth_sqc); 
  vth_sqr_l = ser_2x_p1_surfx2_eval_quad_node_1_l(vth_sqr); 
  vthQuad_l = fmax(sqrt(fabs(vth_sql_r)), sqrt(fabs(vth_sqc_l))); 
  vthQuad_r = fmax(sqrt(fabs(vth_sqc_r)), sqrt(fabs(vth_sqr_l))); 
  max_speed_quad_l[1] = uQuad_l + vthQuad_l; 
  max_speed_quad_r[1] = uQuad_r + vthQuad_r; 
  ser_2x_p1_upwind_quad_to_modal(max_speed_quad_l, max_speed_modal_l); 
  ser_2x_p1_upwind_quad_to_modal(max_speed_quad_r, max_speed_modal_r); 
  outrhou0[0] += 0.5*jump_rhoux_r[1]*max_speed_modal_r[1]*dx1-0.5*jump_rhoux_l[1]*max_speed_modal_l[1]*dx1-0.5*avg_ux_r[1]*flux_rho_r[1]*dx1+0.5*avg_ux_l[1]*flux_rho_l[1]*dx1+0.5*jump_rhoux_r[0]*max_speed_modal_r[0]*dx1-0.5*jump_rhoux_l[0]*max_speed_modal_l[0]*dx1-0.5*avg_ux_r[0]*flux_rho_r[0]*dx1+0.5*avg_ux_l[0]*flux_rho_l[0]*dx1; 
  outrhou0[1] += 0.5*jump_rhoux_r[0]*max_speed_modal_r[1]*dx1-0.5*jump_rhoux_l[0]*max_speed_modal_l[1]*dx1+0.5*max_speed_modal_r[0]*jump_rhoux_r[1]*dx1-0.5*max_speed_modal_l[0]*jump_rhoux_l[1]*dx1-0.5*avg_ux_r[0]*flux_rho_r[1]*dx1+0.5*avg_ux_l[0]*flux_rho_l[1]*dx1-0.5*flux_rho_r[0]*avg_ux_r[1]*dx1+0.5*flux_rho_l[0]*avg_ux_l[1]*dx1; 
  outrhou0[2] += 0.8660254037844386*jump_rhoux_r[1]*max_speed_modal_r[1]*dx1+0.8660254037844386*jump_rhoux_l[1]*max_speed_modal_l[1]*dx1-0.8660254037844386*avg_ux_r[1]*flux_rho_r[1]*dx1-0.8660254037844386*avg_ux_l[1]*flux_rho_l[1]*dx1+0.8660254037844386*jump_rhoux_r[0]*max_speed_modal_r[0]*dx1+0.8660254037844386*jump_rhoux_l[0]*max_speed_modal_l[0]*dx1-0.8660254037844386*avg_ux_r[0]*flux_rho_r[0]*dx1-0.8660254037844386*avg_ux_l[0]*flux_rho_l[0]*dx1; 
  outrhou0[3] += 0.8660254037844386*jump_rhoux_r[0]*max_speed_modal_r[1]*dx1+0.8660254037844386*jump_rhoux_l[0]*max_speed_modal_l[1]*dx1+0.8660254037844386*max_speed_modal_r[0]*jump_rhoux_r[1]*dx1+0.8660254037844386*max_speed_modal_l[0]*jump_rhoux_l[1]*dx1-0.8660254037844386*avg_ux_r[0]*flux_rho_r[1]*dx1-0.8660254037844386*avg_ux_l[0]*flux_rho_l[1]*dx1-0.8660254037844386*flux_rho_r[0]*avg_ux_r[1]*dx1-0.8660254037844386*flux_rho_l[0]*avg_ux_l[1]*dx1; 

  outrhou1[0] += 0.5*jump_rhouy_r[1]*max_speed_modal_r[1]*dx1-0.5*jump_rhouy_l[1]*max_speed_modal_l[1]*dx1-0.5*avg_uy_r[1]*flux_rho_r[1]*dx1+0.5*avg_uy_l[1]*flux_rho_l[1]*dx1+0.5*jump_rhouy_r[0]*max_speed_modal_r[0]*dx1-0.5*jump_rhouy_l[0]*max_speed_modal_l[0]*dx1-0.5*avg_uy_r[0]*flux_rho_r[0]*dx1+0.5*avg_uy_l[0]*flux_rho_l[0]*dx1; 
  outrhou1[1] += 0.5*jump_rhouy_r[0]*max_speed_modal_r[1]*dx1-0.5*jump_rhouy_l[0]*max_speed_modal_l[1]*dx1+0.5*max_speed_modal_r[0]*jump_rhouy_r[1]*dx1-0.5*max_speed_modal_l[0]*jump_rhouy_l[1]*dx1-0.5*avg_uy_r[0]*flux_rho_r[1]*dx1+0.5*avg_uy_l[0]*flux_rho_l[1]*dx1-0.5*flux_rho_r[0]*avg_uy_r[1]*dx1+0.5*flux_rho_l[0]*avg_uy_l[1]*dx1; 
  outrhou1[2] += 0.8660254037844386*jump_rhouy_r[1]*max_speed_modal_r[1]*dx1+0.8660254037844386*jump_rhouy_l[1]*max_speed_modal_l[1]*dx1-0.8660254037844386*avg_uy_r[1]*flux_rho_r[1]*dx1-0.8660254037844386*avg_uy_l[1]*flux_rho_l[1]*dx1+0.8660254037844386*jump_rhouy_r[0]*max_speed_modal_r[0]*dx1+0.8660254037844386*jump_rhouy_l[0]*max_speed_modal_l[0]*dx1-0.8660254037844386*avg_uy_r[0]*flux_rho_r[0]*dx1-0.8660254037844386*avg_uy_l[0]*flux_rho_l[0]*dx1; 
  outrhou1[3] += 0.8660254037844386*jump_rhouy_r[0]*max_speed_modal_r[1]*dx1+0.8660254037844386*jump_rhouy_l[0]*max_speed_modal_l[1]*dx1+0.8660254037844386*max_speed_modal_r[0]*jump_rhouy_r[1]*dx1+0.8660254037844386*max_speed_modal_l[0]*jump_rhouy_l[1]*dx1-0.8660254037844386*avg_uy_r[0]*flux_rho_r[1]*dx1-0.8660254037844386*avg_uy_l[0]*flux_rho_l[1]*dx1-0.8660254037844386*flux_rho_r[0]*avg_uy_r[1]*dx1-0.8660254037844386*flux_rho_l[0]*avg_uy_l[1]*dx1; 

  outrhou2[0] += 0.5*jump_rhouz_r[1]*max_speed_modal_r[1]*dx1-0.5*jump_rhouz_l[1]*max_speed_modal_l[1]*dx1-0.5*avg_uz_r[1]*flux_rho_r[1]*dx1+0.5*avg_uz_l[1]*flux_rho_l[1]*dx1+0.5*jump_rhouz_r[0]*max_speed_modal_r[0]*dx1-0.5*jump_rhouz_l[0]*max_speed_modal_l[0]*dx1-0.5*avg_uz_r[0]*flux_rho_r[0]*dx1+0.5*avg_uz_l[0]*flux_rho_l[0]*dx1; 
  outrhou2[1] += 0.5*jump_rhouz_r[0]*max_speed_modal_r[1]*dx1-0.5*jump_rhouz_l[0]*max_speed_modal_l[1]*dx1+0.5*max_speed_modal_r[0]*jump_rhouz_r[1]*dx1-0.5*max_speed_modal_l[0]*jump_rhouz_l[1]*dx1-0.5*avg_uz_r[0]*flux_rho_r[1]*dx1+0.5*avg_uz_l[0]*flux_rho_l[1]*dx1-0.5*flux_rho_r[0]*avg_uz_r[1]*dx1+0.5*flux_rho_l[0]*avg_uz_l[1]*dx1; 
  outrhou2[2] += 0.8660254037844386*jump_rhouz_r[1]*max_speed_modal_r[1]*dx1+0.8660254037844386*jump_rhouz_l[1]*max_speed_modal_l[1]*dx1-0.8660254037844386*avg_uz_r[1]*flux_rho_r[1]*dx1-0.8660254037844386*avg_uz_l[1]*flux_rho_l[1]*dx1+0.8660254037844386*jump_rhouz_r[0]*max_speed_modal_r[0]*dx1+0.8660254037844386*jump_rhouz_l[0]*max_speed_modal_l[0]*dx1-0.8660254037844386*avg_uz_r[0]*flux_rho_r[0]*dx1-0.8660254037844386*avg_uz_l[0]*flux_rho_l[0]*dx1; 
  outrhou2[3] += 0.8660254037844386*jump_rhouz_r[0]*max_speed_modal_r[1]*dx1+0.8660254037844386*jump_rhouz_l[0]*max_speed_modal_l[1]*dx1+0.8660254037844386*max_speed_modal_r[0]*jump_rhouz_r[1]*dx1+0.8660254037844386*max_speed_modal_l[0]*jump_rhouz_l[1]*dx1-0.8660254037844386*avg_uz_r[0]*flux_rho_r[1]*dx1-0.8660254037844386*avg_uz_l[0]*flux_rho_l[1]*dx1-0.8660254037844386*flux_rho_r[0]*avg_uz_r[1]*dx1-0.8660254037844386*flux_rho_l[0]*avg_uz_l[1]*dx1; 

  return 0.;

} 
