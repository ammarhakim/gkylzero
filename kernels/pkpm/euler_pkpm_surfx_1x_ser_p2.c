#include <gkyl_euler_pkpm_kernels.h> 
#include <gkyl_basis_ser_1x_p2_surfx1_eval_quad.h> 
GKYL_CU_DH void euler_pkpm_surfx_1x_ser_p2(const double *w, const double *dxv, 
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
  const double *rhouy_l = &statevecl[3]; 
  const double *rhouz_l = &statevecl[6]; 

  const double *rhoux_c = &statevecc[0]; 
  const double *rhouy_c = &statevecc[3]; 
  const double *rhouz_c = &statevecc[6]; 

  const double *rhoux_r = &statevecr[0]; 
  const double *rhouy_r = &statevecr[3]; 
  const double *rhouz_r = &statevecr[6]; 

  const double *rho_l = &vlasov_pkpm_momsl[0]; 
  const double *rho_c = &vlasov_pkpm_momsc[0]; 
  const double *rho_r = &vlasov_pkpm_momsr[0]; 

  const double *ux_l = &u_il[0]; 
  const double *ux_c = &u_ic[0]; 
  const double *ux_r = &u_ir[0]; 

  const double *uy_l = &u_il[3]; 
  const double *uy_c = &u_ic[3]; 
  const double *uy_r = &u_ir[3]; 

  const double *uz_l = &u_il[6]; 
  const double *uz_c = &u_ic[6]; 
  const double *uz_r = &u_ir[6]; 

  const double *vth_sql = &T_ijl[0]; 
  const double *vth_sqc = &T_ijc[0]; 
  const double *vth_sqr = &T_ijr[0]; 

  double *outrhou0 = &out[0]; 
  double *outrhou1 = &out[3]; 
  double *outrhou2 = &out[6]; 

  double uxrec_l = 0.3458741190809163*ux_l[2]+0.3458741190809163*ux_c[2]+0.4975526040028326*ux_l[1]-0.4975526040028326*ux_c[1]+0.3535533905932737*ux_l[0]+0.3535533905932737*ux_c[0]; 
  double uxrec_r = 0.3458741190809163*ux_r[2]+0.3458741190809163*ux_c[2]-0.4975526040028326*ux_r[1]+0.4975526040028326*ux_c[1]+0.3535533905932737*ux_r[0]+0.3535533905932737*ux_c[0]; 
  double uxrec_l2 = uxrec_l*uxrec_l; 
  double uxrec_r2 = uxrec_r*uxrec_r; 

  double uyrec_l = 0.3458741190809163*uy_l[2]+0.3458741190809163*uy_c[2]+0.4975526040028326*uy_l[1]-0.4975526040028326*uy_c[1]+0.3535533905932737*uy_l[0]+0.3535533905932737*uy_c[0]; 
  double uyrec_r = 0.3458741190809163*uy_r[2]+0.3458741190809163*uy_c[2]-0.4975526040028326*uy_r[1]+0.4975526040028326*uy_c[1]+0.3535533905932737*uy_r[0]+0.3535533905932737*uy_c[0]; 

  double uzrec_l = 0.3458741190809163*uz_l[2]+0.3458741190809163*uz_c[2]+0.4975526040028326*uz_l[1]-0.4975526040028326*uz_c[1]+0.3535533905932737*uz_l[0]+0.3535533905932737*uz_c[0]; 
  double uzrec_r = 0.3458741190809163*uz_r[2]+0.3458741190809163*uz_c[2]-0.4975526040028326*uz_r[1]+0.4975526040028326*uz_c[1]+0.3535533905932737*uz_r[0]+0.3535533905932737*uz_c[0]; 

  double Ghat_rhoux_l = 0.0; 
  double Ghat_rhoux_r = 0.0; 
  double Ghat_rhouy_l = 0.0; 
  double Ghat_rhouy_r = 0.0; 
  double Ghat_rhouz_l = 0.0; 
  double Ghat_rhouz_r = 0.0; 
  double ux_l_r = ser_1x_p2_surfx1_eval_quad_node_0_r(ux_l); 
  double ux_c_l = ser_1x_p2_surfx1_eval_quad_node_0_l(ux_c); 
  double ux_c_r = ser_1x_p2_surfx1_eval_quad_node_0_r(ux_c); 
  double ux_r_l = ser_1x_p2_surfx1_eval_quad_node_0_l(ux_r); 
  double vth_sq_l_r = ser_1x_p2_surfx1_eval_quad_node_0_r(vth_sql); 
  double vth_sq_c_l = ser_1x_p2_surfx1_eval_quad_node_0_l(vth_sqc); 
  double vth_sq_c_r = ser_1x_p2_surfx1_eval_quad_node_0_r(vth_sqc); 
  double vth_sq_r_l = ser_1x_p2_surfx1_eval_quad_node_0_l(vth_sqr); 

  double ux_max_l = fmax(fabs(ux_l_r), fabs(ux_c_l)); 
  double ux_max_r = fmax(fabs(ux_c_r), fabs(ux_r_l)); 
  double vth_max_l = fmax(sqrt(fabs(vth_sq_l_r)), sqrt(fabs(vth_sq_c_l))); 
  double vth_max_r = fmax(sqrt(fabs(vth_sq_c_r)), sqrt(fabs(vth_sq_r_l))); 
  double max_speed_l = ux_max_l + vth_max_l; 
  double max_speed_r = ux_max_r + vth_max_r; 

  if (uxrec_l > 0) { 
  Ghat_rhoux_l = (1.58113883008419*rho_l[2]+1.224744871391589*rho_l[1]+0.7071067811865475*rho_l[0])*uxrec_l2+0.7905694150420948*rhoux_l[2]*max_speed_l-0.7905694150420948*rhoux_c[2]*max_speed_l+0.6123724356957945*rhoux_l[1]*max_speed_l+0.6123724356957945*rhoux_c[1]*max_speed_l+0.3535533905932737*rhoux_l[0]*max_speed_l-0.3535533905932737*rhoux_c[0]*max_speed_l; 
  Ghat_rhouy_l = (1.58113883008419*rho_l[2]+1.224744871391589*rho_l[1]+0.7071067811865475*rho_l[0])*uxrec_l*uyrec_l+0.7905694150420948*rhouy_l[2]*max_speed_l-0.7905694150420948*rhouy_c[2]*max_speed_l+0.6123724356957945*rhouy_l[1]*max_speed_l+0.6123724356957945*rhouy_c[1]*max_speed_l+0.3535533905932737*rhouy_l[0]*max_speed_l-0.3535533905932737*rhouy_c[0]*max_speed_l; 
  Ghat_rhouz_l = (1.58113883008419*rho_l[2]+1.224744871391589*rho_l[1]+0.7071067811865475*rho_l[0])*uxrec_l*uzrec_l+0.7905694150420948*rhouz_l[2]*max_speed_l-0.7905694150420948*rhouz_c[2]*max_speed_l+0.6123724356957945*rhouz_l[1]*max_speed_l+0.6123724356957945*rhouz_c[1]*max_speed_l+0.3535533905932737*rhouz_l[0]*max_speed_l-0.3535533905932737*rhouz_c[0]*max_speed_l; 
  } else { 
  Ghat_rhoux_l = (1.58113883008419*rho_c[2]-1.224744871391589*rho_c[1]+0.7071067811865475*rho_c[0])*uxrec_l2+0.7905694150420948*rhoux_l[2]*max_speed_l-0.7905694150420948*rhoux_c[2]*max_speed_l+0.6123724356957945*rhoux_l[1]*max_speed_l+0.6123724356957945*rhoux_c[1]*max_speed_l+0.3535533905932737*rhoux_l[0]*max_speed_l-0.3535533905932737*rhoux_c[0]*max_speed_l; 
  Ghat_rhouy_l = (1.58113883008419*rho_c[2]-1.224744871391589*rho_c[1]+0.7071067811865475*rho_c[0])*uxrec_l*uyrec_l+0.7905694150420948*rhouy_l[2]*max_speed_l-0.7905694150420948*rhouy_c[2]*max_speed_l+0.6123724356957945*rhouy_l[1]*max_speed_l+0.6123724356957945*rhouy_c[1]*max_speed_l+0.3535533905932737*rhouy_l[0]*max_speed_l-0.3535533905932737*rhouy_c[0]*max_speed_l; 
  Ghat_rhouz_l = (1.58113883008419*rho_c[2]-1.224744871391589*rho_c[1]+0.7071067811865475*rho_c[0])*uxrec_l*uzrec_l+0.7905694150420948*rhouz_l[2]*max_speed_l-0.7905694150420948*rhouz_c[2]*max_speed_l+0.6123724356957945*rhouz_l[1]*max_speed_l+0.6123724356957945*rhouz_c[1]*max_speed_l+0.3535533905932737*rhouz_l[0]*max_speed_l-0.3535533905932737*rhouz_c[0]*max_speed_l; 
  } 
  if (uxrec_r > 0) { 
  Ghat_rhoux_r = (1.58113883008419*rho_c[2]+1.224744871391589*rho_c[1]+0.7071067811865475*rho_c[0])*uxrec_r2-0.7905694150420948*rhoux_r[2]*max_speed_r+0.7905694150420948*rhoux_c[2]*max_speed_r+0.6123724356957945*rhoux_r[1]*max_speed_r+0.6123724356957945*rhoux_c[1]*max_speed_r-0.3535533905932737*rhoux_r[0]*max_speed_r+0.3535533905932737*rhoux_c[0]*max_speed_r; 
  Ghat_rhouy_r = (1.58113883008419*rho_c[2]+1.224744871391589*rho_c[1]+0.7071067811865475*rho_c[0])*uxrec_r*uyrec_r-0.7905694150420948*rhouy_r[2]*max_speed_r+0.7905694150420948*rhouy_c[2]*max_speed_r+0.6123724356957945*rhouy_r[1]*max_speed_r+0.6123724356957945*rhouy_c[1]*max_speed_r-0.3535533905932737*rhouy_r[0]*max_speed_r+0.3535533905932737*rhouy_c[0]*max_speed_r; 
  Ghat_rhouz_r = (1.58113883008419*rho_c[2]+1.224744871391589*rho_c[1]+0.7071067811865475*rho_c[0])*uxrec_r*uzrec_r-0.7905694150420948*rhouz_r[2]*max_speed_r+0.7905694150420948*rhouz_c[2]*max_speed_r+0.6123724356957945*rhouz_r[1]*max_speed_r+0.6123724356957945*rhouz_c[1]*max_speed_r-0.3535533905932737*rhouz_r[0]*max_speed_r+0.3535533905932737*rhouz_c[0]*max_speed_r; 
  } else { 
  Ghat_rhoux_r = (1.58113883008419*rho_r[2]-1.224744871391589*rho_r[1]+0.7071067811865475*rho_r[0])*uxrec_r2-0.7905694150420948*rhoux_r[2]*max_speed_r+0.7905694150420948*rhoux_c[2]*max_speed_r+0.6123724356957945*rhoux_r[1]*max_speed_r+0.6123724356957945*rhoux_c[1]*max_speed_r-0.3535533905932737*rhoux_r[0]*max_speed_r+0.3535533905932737*rhoux_c[0]*max_speed_r; 
  Ghat_rhouy_r = (1.58113883008419*rho_r[2]-1.224744871391589*rho_r[1]+0.7071067811865475*rho_r[0])*uxrec_r*uyrec_r-0.7905694150420948*rhouy_r[2]*max_speed_r+0.7905694150420948*rhouy_c[2]*max_speed_r+0.6123724356957945*rhouy_r[1]*max_speed_r+0.6123724356957945*rhouy_c[1]*max_speed_r-0.3535533905932737*rhouy_r[0]*max_speed_r+0.3535533905932737*rhouy_c[0]*max_speed_r; 
  Ghat_rhouz_r = (1.58113883008419*rho_r[2]-1.224744871391589*rho_r[1]+0.7071067811865475*rho_r[0])*uxrec_r*uzrec_r-0.7905694150420948*rhouz_r[2]*max_speed_r+0.7905694150420948*rhouz_c[2]*max_speed_r+0.6123724356957945*rhouz_r[1]*max_speed_r+0.6123724356957945*rhouz_c[1]*max_speed_r-0.3535533905932737*rhouz_r[0]*max_speed_r+0.3535533905932737*rhouz_c[0]*max_speed_r; 
  } 

  outrhou0[0] += 0.7071067811865475*Ghat_rhoux_l*dx1-0.7071067811865475*Ghat_rhoux_r*dx1; 
  outrhou0[1] += (-1.224744871391589*Ghat_rhoux_r*dx1)-1.224744871391589*Ghat_rhoux_l*dx1; 
  outrhou0[2] += 1.58113883008419*Ghat_rhoux_l*dx1-1.58113883008419*Ghat_rhoux_r*dx1; 

  outrhou1[0] += 0.7071067811865475*Ghat_rhouy_l*dx1-0.7071067811865475*Ghat_rhouy_r*dx1; 
  outrhou1[1] += (-1.224744871391589*Ghat_rhouy_r*dx1)-1.224744871391589*Ghat_rhouy_l*dx1; 
  outrhou1[2] += 1.58113883008419*Ghat_rhouy_l*dx1-1.58113883008419*Ghat_rhouy_r*dx1; 

  outrhou2[0] += 0.7071067811865475*Ghat_rhouz_l*dx1-0.7071067811865475*Ghat_rhouz_r*dx1; 
  outrhou2[1] += (-1.224744871391589*Ghat_rhouz_r*dx1)-1.224744871391589*Ghat_rhouz_l*dx1; 
  outrhou2[2] += 1.58113883008419*Ghat_rhouz_l*dx1-1.58113883008419*Ghat_rhouz_r*dx1; 

} 
