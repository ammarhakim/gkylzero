#include <gkyl_euler_pkpm_kernels.h> 
#include <gkyl_basis_ser_3x_p1_surfx3_eval_quad.h> 
#include <gkyl_basis_ser_3x_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double euler_pkpm_surfz_3x_ser_p1(const double *w, const double *dxv,
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

  const double dx1 = 2.0/dxv[2]; 

  const double *rhoux_l = &statevecl[0]; 
  const double *rhouy_l = &statevecl[8]; 
  const double *rhouz_l = &statevecl[16]; 

  const double *rhoux_c = &statevecc[0]; 
  const double *rhouy_c = &statevecc[8]; 
  const double *rhouz_c = &statevecc[16]; 

  const double *rhoux_r = &statevecr[0]; 
  const double *rhouy_r = &statevecr[8]; 
  const double *rhouz_r = &statevecr[16]; 

  const double *rho_l = &vlasov_pkpm_momsl[0]; 
  const double *rho_c = &vlasov_pkpm_momsc[0]; 
  const double *rho_r = &vlasov_pkpm_momsr[0]; 

  const double *ux_l = &u_il[0]; 
  const double *ux_c = &u_ic[0]; 
  const double *ux_r = &u_ir[0]; 

  const double *uy_l = &u_il[8]; 
  const double *uy_c = &u_ic[8]; 
  const double *uy_r = &u_ir[8]; 

  const double *uz_l = &u_il[16]; 
  const double *uz_c = &u_ic[16]; 
  const double *uz_r = &u_ir[16]; 

  // Get another pointer for u in direction of update for ease of flux calculation. 
  const double *u_l = &u_il[16]; 
  const double *u_c = &u_ic[16]; 
  const double *u_r = &u_ir[16]; 

  // Get thermal velocity in direction of update for penalization vth^2 = 3.0*T_ii/m. 
  const double *vth_sql = &T_ijl[40]; 
  const double *vth_sqc = &T_ijc[40]; 
  const double *vth_sqr = &T_ijr[40]; 

  double *outrhou0 = &out[0]; 
  double *outrhou1 = &out[8]; 
  double *outrhou2 = &out[16]; 

  double flux_rho_l[4] = {0.0}; 
  double flux_rho_r[4] = {0.0}; 
  double avg_ux_l[4] = {0.0}; 
  double avg_ux_r[4] = {0.0}; 
  double avg_uy_l[4] = {0.0}; 
  double avg_uy_r[4] = {0.0}; 
  double avg_uz_l[4] = {0.0}; 
  double avg_uz_r[4] = {0.0}; 
  double jump_rhoux_l[4] = {0.0}; 
  double jump_rhoux_r[4] = {0.0}; 
  double jump_rhouy_l[4] = {0.0}; 
  double jump_rhouy_r[4] = {0.0}; 
  double jump_rhouz_l[4] = {0.0}; 
  double jump_rhouz_r[4] = {0.0}; 
  double max_speed_quad_l[4] = {0.0}; 
  double max_speed_quad_r[4] = {0.0}; 
  double max_speed_modal_l[4] = {0.0}; 
  double max_speed_modal_r[4] = {0.0}; 
  flux_rho_l[0] = 0.1875*rho_l[7]*u_l[7]-0.1875*rho_c[7]*u_l[7]+0.1082531754730548*rho_l[4]*u_l[7]+0.1082531754730548*rho_c[4]*u_l[7]-0.1875*rho_l[7]*u_c[7]+0.1875*rho_c[7]*u_c[7]-0.1082531754730548*rho_l[4]*u_c[7]-0.1082531754730548*rho_c[4]*u_c[7]+0.1082531754730548*u_l[4]*rho_l[7]+0.1082531754730548*u_c[4]*rho_l[7]-0.1082531754730548*u_l[4]*rho_c[7]-0.1082531754730548*u_c[4]*rho_c[7]+0.1875*rho_l[6]*u_l[6]-0.1875*rho_c[6]*u_l[6]+0.1082531754730548*rho_l[2]*u_l[6]+0.1082531754730548*rho_c[2]*u_l[6]-0.1875*rho_l[6]*u_c[6]+0.1875*rho_c[6]*u_c[6]-0.1082531754730548*rho_l[2]*u_c[6]-0.1082531754730548*rho_c[2]*u_c[6]+0.1082531754730548*u_l[2]*rho_l[6]+0.1082531754730548*u_c[2]*rho_l[6]-0.1082531754730548*u_l[2]*rho_c[6]-0.1082531754730548*u_c[2]*rho_c[6]+0.1875*rho_l[5]*u_l[5]-0.1875*rho_c[5]*u_l[5]+0.1082531754730548*rho_l[1]*u_l[5]+0.1082531754730548*rho_c[1]*u_l[5]-0.1875*rho_l[5]*u_c[5]+0.1875*rho_c[5]*u_c[5]-0.1082531754730548*rho_l[1]*u_c[5]-0.1082531754730548*rho_c[1]*u_c[5]+0.1082531754730548*u_l[1]*rho_l[5]+0.1082531754730548*u_c[1]*rho_l[5]-0.1082531754730548*u_l[1]*rho_c[5]-0.1082531754730548*u_c[1]*rho_c[5]+0.0625*rho_l[4]*u_l[4]+0.0625*rho_c[4]*u_l[4]+0.0625*rho_l[4]*u_c[4]+0.0625*rho_c[4]*u_c[4]+0.1875*rho_l[3]*u_l[3]-0.1875*rho_c[3]*u_l[3]+0.1082531754730548*rho_l[0]*u_l[3]+0.1082531754730548*rho_c[0]*u_l[3]-0.1875*rho_l[3]*u_c[3]+0.1875*rho_c[3]*u_c[3]-0.1082531754730548*rho_l[0]*u_c[3]-0.1082531754730548*rho_c[0]*u_c[3]+0.1082531754730548*u_l[0]*rho_l[3]+0.1082531754730548*u_c[0]*rho_l[3]-0.1082531754730548*u_l[0]*rho_c[3]-0.1082531754730548*u_c[0]*rho_c[3]+0.0625*rho_l[2]*u_l[2]+0.0625*rho_c[2]*u_l[2]+0.0625*rho_l[2]*u_c[2]+0.0625*rho_c[2]*u_c[2]+0.0625*rho_l[1]*u_l[1]+0.0625*rho_c[1]*u_l[1]+0.0625*rho_l[1]*u_c[1]+0.0625*rho_c[1]*u_c[1]+0.0625*rho_l[0]*u_l[0]+0.0625*rho_c[0]*u_l[0]+0.0625*rho_l[0]*u_c[0]+0.0625*rho_c[0]*u_c[0]; 
  flux_rho_l[1] = 0.1875*rho_l[6]*u_l[7]-0.1875*rho_c[6]*u_l[7]+0.1082531754730548*rho_l[2]*u_l[7]+0.1082531754730548*rho_c[2]*u_l[7]-0.1875*rho_l[6]*u_c[7]+0.1875*rho_c[6]*u_c[7]-0.1082531754730548*rho_l[2]*u_c[7]-0.1082531754730548*rho_c[2]*u_c[7]+0.1875*u_l[6]*rho_l[7]-0.1875*u_c[6]*rho_l[7]+0.1082531754730548*u_l[2]*rho_l[7]+0.1082531754730548*u_c[2]*rho_l[7]-0.1875*u_l[6]*rho_c[7]+0.1875*u_c[6]*rho_c[7]-0.1082531754730548*u_l[2]*rho_c[7]-0.1082531754730548*u_c[2]*rho_c[7]+0.1082531754730548*rho_l[4]*u_l[6]+0.1082531754730548*rho_c[4]*u_l[6]-0.1082531754730548*rho_l[4]*u_c[6]-0.1082531754730548*rho_c[4]*u_c[6]+0.1082531754730548*u_l[4]*rho_l[6]+0.1082531754730548*u_c[4]*rho_l[6]-0.1082531754730548*u_l[4]*rho_c[6]-0.1082531754730548*u_c[4]*rho_c[6]+0.1875*rho_l[3]*u_l[5]-0.1875*rho_c[3]*u_l[5]+0.1082531754730548*rho_l[0]*u_l[5]+0.1082531754730548*rho_c[0]*u_l[5]-0.1875*rho_l[3]*u_c[5]+0.1875*rho_c[3]*u_c[5]-0.1082531754730548*rho_l[0]*u_c[5]-0.1082531754730548*rho_c[0]*u_c[5]+0.1875*u_l[3]*rho_l[5]-0.1875*u_c[3]*rho_l[5]+0.1082531754730548*u_l[0]*rho_l[5]+0.1082531754730548*u_c[0]*rho_l[5]-0.1875*u_l[3]*rho_c[5]+0.1875*u_c[3]*rho_c[5]-0.1082531754730548*u_l[0]*rho_c[5]-0.1082531754730548*u_c[0]*rho_c[5]+0.0625*rho_l[2]*u_l[4]+0.0625*rho_c[2]*u_l[4]+0.0625*rho_l[2]*u_c[4]+0.0625*rho_c[2]*u_c[4]+0.0625*u_l[2]*rho_l[4]+0.0625*u_c[2]*rho_l[4]+0.0625*u_l[2]*rho_c[4]+0.0625*u_c[2]*rho_c[4]+0.1082531754730548*rho_l[1]*u_l[3]+0.1082531754730548*rho_c[1]*u_l[3]-0.1082531754730548*rho_l[1]*u_c[3]-0.1082531754730548*rho_c[1]*u_c[3]+0.1082531754730548*u_l[1]*rho_l[3]+0.1082531754730548*u_c[1]*rho_l[3]-0.1082531754730548*u_l[1]*rho_c[3]-0.1082531754730548*u_c[1]*rho_c[3]+0.0625*rho_l[0]*u_l[1]+0.0625*rho_c[0]*u_l[1]+0.0625*rho_l[0]*u_c[1]+0.0625*rho_c[0]*u_c[1]+0.0625*u_l[0]*rho_l[1]+0.0625*u_c[0]*rho_l[1]+0.0625*u_l[0]*rho_c[1]+0.0625*u_c[0]*rho_c[1]; 
  flux_rho_l[2] = 0.1875*rho_l[5]*u_l[7]-0.1875*rho_c[5]*u_l[7]+0.1082531754730548*rho_l[1]*u_l[7]+0.1082531754730548*rho_c[1]*u_l[7]-0.1875*rho_l[5]*u_c[7]+0.1875*rho_c[5]*u_c[7]-0.1082531754730548*rho_l[1]*u_c[7]-0.1082531754730548*rho_c[1]*u_c[7]+0.1875*u_l[5]*rho_l[7]-0.1875*u_c[5]*rho_l[7]+0.1082531754730548*u_l[1]*rho_l[7]+0.1082531754730548*u_c[1]*rho_l[7]-0.1875*u_l[5]*rho_c[7]+0.1875*u_c[5]*rho_c[7]-0.1082531754730548*u_l[1]*rho_c[7]-0.1082531754730548*u_c[1]*rho_c[7]+0.1875*rho_l[3]*u_l[6]-0.1875*rho_c[3]*u_l[6]+0.1082531754730548*rho_l[0]*u_l[6]+0.1082531754730548*rho_c[0]*u_l[6]-0.1875*rho_l[3]*u_c[6]+0.1875*rho_c[3]*u_c[6]-0.1082531754730548*rho_l[0]*u_c[6]-0.1082531754730548*rho_c[0]*u_c[6]+0.1875*u_l[3]*rho_l[6]-0.1875*u_c[3]*rho_l[6]+0.1082531754730548*u_l[0]*rho_l[6]+0.1082531754730548*u_c[0]*rho_l[6]-0.1875*u_l[3]*rho_c[6]+0.1875*u_c[3]*rho_c[6]-0.1082531754730548*u_l[0]*rho_c[6]-0.1082531754730548*u_c[0]*rho_c[6]+0.1082531754730548*rho_l[4]*u_l[5]+0.1082531754730548*rho_c[4]*u_l[5]-0.1082531754730548*rho_l[4]*u_c[5]-0.1082531754730548*rho_c[4]*u_c[5]+0.1082531754730548*u_l[4]*rho_l[5]+0.1082531754730548*u_c[4]*rho_l[5]-0.1082531754730548*u_l[4]*rho_c[5]-0.1082531754730548*u_c[4]*rho_c[5]+0.0625*rho_l[1]*u_l[4]+0.0625*rho_c[1]*u_l[4]+0.0625*rho_l[1]*u_c[4]+0.0625*rho_c[1]*u_c[4]+0.0625*u_l[1]*rho_l[4]+0.0625*u_c[1]*rho_l[4]+0.0625*u_l[1]*rho_c[4]+0.0625*u_c[1]*rho_c[4]+0.1082531754730548*rho_l[2]*u_l[3]+0.1082531754730548*rho_c[2]*u_l[3]-0.1082531754730548*rho_l[2]*u_c[3]-0.1082531754730548*rho_c[2]*u_c[3]+0.1082531754730548*u_l[2]*rho_l[3]+0.1082531754730548*u_c[2]*rho_l[3]-0.1082531754730548*u_l[2]*rho_c[3]-0.1082531754730548*u_c[2]*rho_c[3]+0.0625*rho_l[0]*u_l[2]+0.0625*rho_c[0]*u_l[2]+0.0625*rho_l[0]*u_c[2]+0.0625*rho_c[0]*u_c[2]+0.0625*u_l[0]*rho_l[2]+0.0625*u_c[0]*rho_l[2]+0.0625*u_l[0]*rho_c[2]+0.0625*u_c[0]*rho_c[2]; 
  flux_rho_l[3] = 0.1875*rho_l[3]*u_l[7]-0.1875*rho_c[3]*u_l[7]+0.1082531754730548*rho_l[0]*u_l[7]+0.1082531754730548*rho_c[0]*u_l[7]-0.1875*rho_l[3]*u_c[7]+0.1875*rho_c[3]*u_c[7]-0.1082531754730548*rho_l[0]*u_c[7]-0.1082531754730548*rho_c[0]*u_c[7]+0.1875*u_l[3]*rho_l[7]-0.1875*u_c[3]*rho_l[7]+0.1082531754730548*u_l[0]*rho_l[7]+0.1082531754730548*u_c[0]*rho_l[7]-0.1875*u_l[3]*rho_c[7]+0.1875*u_c[3]*rho_c[7]-0.1082531754730548*u_l[0]*rho_c[7]-0.1082531754730548*u_c[0]*rho_c[7]+0.1875*rho_l[5]*u_l[6]-0.1875*rho_c[5]*u_l[6]+0.1082531754730548*rho_l[1]*u_l[6]+0.1082531754730548*rho_c[1]*u_l[6]-0.1875*rho_l[5]*u_c[6]+0.1875*rho_c[5]*u_c[6]-0.1082531754730548*rho_l[1]*u_c[6]-0.1082531754730548*rho_c[1]*u_c[6]+0.1875*u_l[5]*rho_l[6]-0.1875*u_c[5]*rho_l[6]+0.1082531754730548*u_l[1]*rho_l[6]+0.1082531754730548*u_c[1]*rho_l[6]-0.1875*u_l[5]*rho_c[6]+0.1875*u_c[5]*rho_c[6]-0.1082531754730548*u_l[1]*rho_c[6]-0.1082531754730548*u_c[1]*rho_c[6]+0.1082531754730548*rho_l[2]*u_l[5]+0.1082531754730548*rho_c[2]*u_l[5]-0.1082531754730548*rho_l[2]*u_c[5]-0.1082531754730548*rho_c[2]*u_c[5]+0.1082531754730548*u_l[2]*rho_l[5]+0.1082531754730548*u_c[2]*rho_l[5]-0.1082531754730548*u_l[2]*rho_c[5]-0.1082531754730548*u_c[2]*rho_c[5]+0.1082531754730548*rho_l[3]*u_l[4]-0.1082531754730548*rho_c[3]*u_l[4]+0.0625*rho_l[0]*u_l[4]+0.0625*rho_c[0]*u_l[4]+0.1082531754730548*rho_l[3]*u_c[4]-0.1082531754730548*rho_c[3]*u_c[4]+0.0625*rho_l[0]*u_c[4]+0.0625*rho_c[0]*u_c[4]+0.1082531754730548*u_l[3]*rho_l[4]-0.1082531754730548*u_c[3]*rho_l[4]+0.0625*u_l[0]*rho_l[4]+0.0625*u_c[0]*rho_l[4]+0.1082531754730548*u_l[3]*rho_c[4]-0.1082531754730548*u_c[3]*rho_c[4]+0.0625*u_l[0]*rho_c[4]+0.0625*u_c[0]*rho_c[4]+0.0625*rho_l[1]*u_l[2]+0.0625*rho_c[1]*u_l[2]+0.0625*rho_l[1]*u_c[2]+0.0625*rho_c[1]*u_c[2]+0.0625*u_l[1]*rho_l[2]+0.0625*u_c[1]*rho_l[2]+0.0625*u_l[1]*rho_c[2]+0.0625*u_c[1]*rho_c[2]; 

  flux_rho_r[0] = 0.1875*rho_r[7]*u_r[7]-0.1875*rho_c[7]*u_r[7]-0.1082531754730548*rho_r[4]*u_r[7]-0.1082531754730548*rho_c[4]*u_r[7]-0.1875*rho_r[7]*u_c[7]+0.1875*rho_c[7]*u_c[7]+0.1082531754730548*rho_r[4]*u_c[7]+0.1082531754730548*rho_c[4]*u_c[7]-0.1082531754730548*u_r[4]*rho_r[7]-0.1082531754730548*u_c[4]*rho_r[7]+0.1082531754730548*u_r[4]*rho_c[7]+0.1082531754730548*u_c[4]*rho_c[7]+0.1875*rho_r[6]*u_r[6]-0.1875*rho_c[6]*u_r[6]-0.1082531754730548*rho_r[2]*u_r[6]-0.1082531754730548*rho_c[2]*u_r[6]-0.1875*rho_r[6]*u_c[6]+0.1875*rho_c[6]*u_c[6]+0.1082531754730548*rho_r[2]*u_c[6]+0.1082531754730548*rho_c[2]*u_c[6]-0.1082531754730548*u_r[2]*rho_r[6]-0.1082531754730548*u_c[2]*rho_r[6]+0.1082531754730548*u_r[2]*rho_c[6]+0.1082531754730548*u_c[2]*rho_c[6]+0.1875*rho_r[5]*u_r[5]-0.1875*rho_c[5]*u_r[5]-0.1082531754730548*rho_r[1]*u_r[5]-0.1082531754730548*rho_c[1]*u_r[5]-0.1875*rho_r[5]*u_c[5]+0.1875*rho_c[5]*u_c[5]+0.1082531754730548*rho_r[1]*u_c[5]+0.1082531754730548*rho_c[1]*u_c[5]-0.1082531754730548*u_r[1]*rho_r[5]-0.1082531754730548*u_c[1]*rho_r[5]+0.1082531754730548*u_r[1]*rho_c[5]+0.1082531754730548*u_c[1]*rho_c[5]+0.0625*rho_r[4]*u_r[4]+0.0625*rho_c[4]*u_r[4]+0.0625*rho_r[4]*u_c[4]+0.0625*rho_c[4]*u_c[4]+0.1875*rho_r[3]*u_r[3]-0.1875*rho_c[3]*u_r[3]-0.1082531754730548*rho_r[0]*u_r[3]-0.1082531754730548*rho_c[0]*u_r[3]-0.1875*rho_r[3]*u_c[3]+0.1875*rho_c[3]*u_c[3]+0.1082531754730548*rho_r[0]*u_c[3]+0.1082531754730548*rho_c[0]*u_c[3]-0.1082531754730548*u_r[0]*rho_r[3]-0.1082531754730548*u_c[0]*rho_r[3]+0.1082531754730548*u_r[0]*rho_c[3]+0.1082531754730548*u_c[0]*rho_c[3]+0.0625*rho_r[2]*u_r[2]+0.0625*rho_c[2]*u_r[2]+0.0625*rho_r[2]*u_c[2]+0.0625*rho_c[2]*u_c[2]+0.0625*rho_r[1]*u_r[1]+0.0625*rho_c[1]*u_r[1]+0.0625*rho_r[1]*u_c[1]+0.0625*rho_c[1]*u_c[1]+0.0625*rho_r[0]*u_r[0]+0.0625*rho_c[0]*u_r[0]+0.0625*rho_r[0]*u_c[0]+0.0625*rho_c[0]*u_c[0]; 
  flux_rho_r[1] = 0.1875*rho_r[6]*u_r[7]-0.1875*rho_c[6]*u_r[7]-0.1082531754730548*rho_r[2]*u_r[7]-0.1082531754730548*rho_c[2]*u_r[7]-0.1875*rho_r[6]*u_c[7]+0.1875*rho_c[6]*u_c[7]+0.1082531754730548*rho_r[2]*u_c[7]+0.1082531754730548*rho_c[2]*u_c[7]+0.1875*u_r[6]*rho_r[7]-0.1875*u_c[6]*rho_r[7]-0.1082531754730548*u_r[2]*rho_r[7]-0.1082531754730548*u_c[2]*rho_r[7]-0.1875*u_r[6]*rho_c[7]+0.1875*u_c[6]*rho_c[7]+0.1082531754730548*u_r[2]*rho_c[7]+0.1082531754730548*u_c[2]*rho_c[7]-0.1082531754730548*rho_r[4]*u_r[6]-0.1082531754730548*rho_c[4]*u_r[6]+0.1082531754730548*rho_r[4]*u_c[6]+0.1082531754730548*rho_c[4]*u_c[6]-0.1082531754730548*u_r[4]*rho_r[6]-0.1082531754730548*u_c[4]*rho_r[6]+0.1082531754730548*u_r[4]*rho_c[6]+0.1082531754730548*u_c[4]*rho_c[6]+0.1875*rho_r[3]*u_r[5]-0.1875*rho_c[3]*u_r[5]-0.1082531754730548*rho_r[0]*u_r[5]-0.1082531754730548*rho_c[0]*u_r[5]-0.1875*rho_r[3]*u_c[5]+0.1875*rho_c[3]*u_c[5]+0.1082531754730548*rho_r[0]*u_c[5]+0.1082531754730548*rho_c[0]*u_c[5]+0.1875*u_r[3]*rho_r[5]-0.1875*u_c[3]*rho_r[5]-0.1082531754730548*u_r[0]*rho_r[5]-0.1082531754730548*u_c[0]*rho_r[5]-0.1875*u_r[3]*rho_c[5]+0.1875*u_c[3]*rho_c[5]+0.1082531754730548*u_r[0]*rho_c[5]+0.1082531754730548*u_c[0]*rho_c[5]+0.0625*rho_r[2]*u_r[4]+0.0625*rho_c[2]*u_r[4]+0.0625*rho_r[2]*u_c[4]+0.0625*rho_c[2]*u_c[4]+0.0625*u_r[2]*rho_r[4]+0.0625*u_c[2]*rho_r[4]+0.0625*u_r[2]*rho_c[4]+0.0625*u_c[2]*rho_c[4]-0.1082531754730548*rho_r[1]*u_r[3]-0.1082531754730548*rho_c[1]*u_r[3]+0.1082531754730548*rho_r[1]*u_c[3]+0.1082531754730548*rho_c[1]*u_c[3]-0.1082531754730548*u_r[1]*rho_r[3]-0.1082531754730548*u_c[1]*rho_r[3]+0.1082531754730548*u_r[1]*rho_c[3]+0.1082531754730548*u_c[1]*rho_c[3]+0.0625*rho_r[0]*u_r[1]+0.0625*rho_c[0]*u_r[1]+0.0625*rho_r[0]*u_c[1]+0.0625*rho_c[0]*u_c[1]+0.0625*u_r[0]*rho_r[1]+0.0625*u_c[0]*rho_r[1]+0.0625*u_r[0]*rho_c[1]+0.0625*u_c[0]*rho_c[1]; 
  flux_rho_r[2] = 0.1875*rho_r[5]*u_r[7]-0.1875*rho_c[5]*u_r[7]-0.1082531754730548*rho_r[1]*u_r[7]-0.1082531754730548*rho_c[1]*u_r[7]-0.1875*rho_r[5]*u_c[7]+0.1875*rho_c[5]*u_c[7]+0.1082531754730548*rho_r[1]*u_c[7]+0.1082531754730548*rho_c[1]*u_c[7]+0.1875*u_r[5]*rho_r[7]-0.1875*u_c[5]*rho_r[7]-0.1082531754730548*u_r[1]*rho_r[7]-0.1082531754730548*u_c[1]*rho_r[7]-0.1875*u_r[5]*rho_c[7]+0.1875*u_c[5]*rho_c[7]+0.1082531754730548*u_r[1]*rho_c[7]+0.1082531754730548*u_c[1]*rho_c[7]+0.1875*rho_r[3]*u_r[6]-0.1875*rho_c[3]*u_r[6]-0.1082531754730548*rho_r[0]*u_r[6]-0.1082531754730548*rho_c[0]*u_r[6]-0.1875*rho_r[3]*u_c[6]+0.1875*rho_c[3]*u_c[6]+0.1082531754730548*rho_r[0]*u_c[6]+0.1082531754730548*rho_c[0]*u_c[6]+0.1875*u_r[3]*rho_r[6]-0.1875*u_c[3]*rho_r[6]-0.1082531754730548*u_r[0]*rho_r[6]-0.1082531754730548*u_c[0]*rho_r[6]-0.1875*u_r[3]*rho_c[6]+0.1875*u_c[3]*rho_c[6]+0.1082531754730548*u_r[0]*rho_c[6]+0.1082531754730548*u_c[0]*rho_c[6]-0.1082531754730548*rho_r[4]*u_r[5]-0.1082531754730548*rho_c[4]*u_r[5]+0.1082531754730548*rho_r[4]*u_c[5]+0.1082531754730548*rho_c[4]*u_c[5]-0.1082531754730548*u_r[4]*rho_r[5]-0.1082531754730548*u_c[4]*rho_r[5]+0.1082531754730548*u_r[4]*rho_c[5]+0.1082531754730548*u_c[4]*rho_c[5]+0.0625*rho_r[1]*u_r[4]+0.0625*rho_c[1]*u_r[4]+0.0625*rho_r[1]*u_c[4]+0.0625*rho_c[1]*u_c[4]+0.0625*u_r[1]*rho_r[4]+0.0625*u_c[1]*rho_r[4]+0.0625*u_r[1]*rho_c[4]+0.0625*u_c[1]*rho_c[4]-0.1082531754730548*rho_r[2]*u_r[3]-0.1082531754730548*rho_c[2]*u_r[3]+0.1082531754730548*rho_r[2]*u_c[3]+0.1082531754730548*rho_c[2]*u_c[3]-0.1082531754730548*u_r[2]*rho_r[3]-0.1082531754730548*u_c[2]*rho_r[3]+0.1082531754730548*u_r[2]*rho_c[3]+0.1082531754730548*u_c[2]*rho_c[3]+0.0625*rho_r[0]*u_r[2]+0.0625*rho_c[0]*u_r[2]+0.0625*rho_r[0]*u_c[2]+0.0625*rho_c[0]*u_c[2]+0.0625*u_r[0]*rho_r[2]+0.0625*u_c[0]*rho_r[2]+0.0625*u_r[0]*rho_c[2]+0.0625*u_c[0]*rho_c[2]; 
  flux_rho_r[3] = 0.1875*rho_r[3]*u_r[7]-0.1875*rho_c[3]*u_r[7]-0.1082531754730548*rho_r[0]*u_r[7]-0.1082531754730548*rho_c[0]*u_r[7]-0.1875*rho_r[3]*u_c[7]+0.1875*rho_c[3]*u_c[7]+0.1082531754730548*rho_r[0]*u_c[7]+0.1082531754730548*rho_c[0]*u_c[7]+0.1875*u_r[3]*rho_r[7]-0.1875*u_c[3]*rho_r[7]-0.1082531754730548*u_r[0]*rho_r[7]-0.1082531754730548*u_c[0]*rho_r[7]-0.1875*u_r[3]*rho_c[7]+0.1875*u_c[3]*rho_c[7]+0.1082531754730548*u_r[0]*rho_c[7]+0.1082531754730548*u_c[0]*rho_c[7]+0.1875*rho_r[5]*u_r[6]-0.1875*rho_c[5]*u_r[6]-0.1082531754730548*rho_r[1]*u_r[6]-0.1082531754730548*rho_c[1]*u_r[6]-0.1875*rho_r[5]*u_c[6]+0.1875*rho_c[5]*u_c[6]+0.1082531754730548*rho_r[1]*u_c[6]+0.1082531754730548*rho_c[1]*u_c[6]+0.1875*u_r[5]*rho_r[6]-0.1875*u_c[5]*rho_r[6]-0.1082531754730548*u_r[1]*rho_r[6]-0.1082531754730548*u_c[1]*rho_r[6]-0.1875*u_r[5]*rho_c[6]+0.1875*u_c[5]*rho_c[6]+0.1082531754730548*u_r[1]*rho_c[6]+0.1082531754730548*u_c[1]*rho_c[6]-0.1082531754730548*rho_r[2]*u_r[5]-0.1082531754730548*rho_c[2]*u_r[5]+0.1082531754730548*rho_r[2]*u_c[5]+0.1082531754730548*rho_c[2]*u_c[5]-0.1082531754730548*u_r[2]*rho_r[5]-0.1082531754730548*u_c[2]*rho_r[5]+0.1082531754730548*u_r[2]*rho_c[5]+0.1082531754730548*u_c[2]*rho_c[5]-0.1082531754730548*rho_r[3]*u_r[4]+0.1082531754730548*rho_c[3]*u_r[4]+0.0625*rho_r[0]*u_r[4]+0.0625*rho_c[0]*u_r[4]-0.1082531754730548*rho_r[3]*u_c[4]+0.1082531754730548*rho_c[3]*u_c[4]+0.0625*rho_r[0]*u_c[4]+0.0625*rho_c[0]*u_c[4]-0.1082531754730548*u_r[3]*rho_r[4]+0.1082531754730548*u_c[3]*rho_r[4]+0.0625*u_r[0]*rho_r[4]+0.0625*u_c[0]*rho_r[4]-0.1082531754730548*u_r[3]*rho_c[4]+0.1082531754730548*u_c[3]*rho_c[4]+0.0625*u_r[0]*rho_c[4]+0.0625*u_c[0]*rho_c[4]+0.0625*rho_r[1]*u_r[2]+0.0625*rho_c[1]*u_r[2]+0.0625*rho_r[1]*u_c[2]+0.0625*rho_c[1]*u_c[2]+0.0625*u_r[1]*rho_r[2]+0.0625*u_c[1]*rho_r[2]+0.0625*u_r[1]*rho_c[2]+0.0625*u_c[1]*rho_c[2]; 

  avg_ux_l[0] = 0.6123724356957944*ux_l[3]-0.6123724356957944*ux_c[3]+0.3535533905932737*ux_l[0]+0.3535533905932737*ux_c[0]; 
  avg_ux_l[1] = 0.6123724356957944*ux_l[5]-0.6123724356957944*ux_c[5]+0.3535533905932737*ux_l[1]+0.3535533905932737*ux_c[1]; 
  avg_ux_l[2] = 0.6123724356957944*ux_l[6]-0.6123724356957944*ux_c[6]+0.3535533905932737*ux_l[2]+0.3535533905932737*ux_c[2]; 
  avg_ux_l[3] = 0.6123724356957944*ux_l[7]-0.6123724356957944*ux_c[7]+0.3535533905932737*ux_l[4]+0.3535533905932737*ux_c[4]; 

  avg_ux_r[0] = (-0.6123724356957944*ux_r[3])+0.6123724356957944*ux_c[3]+0.3535533905932737*ux_r[0]+0.3535533905932737*ux_c[0]; 
  avg_ux_r[1] = (-0.6123724356957944*ux_r[5])+0.6123724356957944*ux_c[5]+0.3535533905932737*ux_r[1]+0.3535533905932737*ux_c[1]; 
  avg_ux_r[2] = (-0.6123724356957944*ux_r[6])+0.6123724356957944*ux_c[6]+0.3535533905932737*ux_r[2]+0.3535533905932737*ux_c[2]; 
  avg_ux_r[3] = (-0.6123724356957944*ux_r[7])+0.6123724356957944*ux_c[7]+0.3535533905932737*ux_r[4]+0.3535533905932737*ux_c[4]; 

  avg_uy_l[0] = 0.6123724356957944*uy_l[3]-0.6123724356957944*uy_c[3]+0.3535533905932737*uy_l[0]+0.3535533905932737*uy_c[0]; 
  avg_uy_l[1] = 0.6123724356957944*uy_l[5]-0.6123724356957944*uy_c[5]+0.3535533905932737*uy_l[1]+0.3535533905932737*uy_c[1]; 
  avg_uy_l[2] = 0.6123724356957944*uy_l[6]-0.6123724356957944*uy_c[6]+0.3535533905932737*uy_l[2]+0.3535533905932737*uy_c[2]; 
  avg_uy_l[3] = 0.6123724356957944*uy_l[7]-0.6123724356957944*uy_c[7]+0.3535533905932737*uy_l[4]+0.3535533905932737*uy_c[4]; 

  avg_uy_r[0] = (-0.6123724356957944*uy_r[3])+0.6123724356957944*uy_c[3]+0.3535533905932737*uy_r[0]+0.3535533905932737*uy_c[0]; 
  avg_uy_r[1] = (-0.6123724356957944*uy_r[5])+0.6123724356957944*uy_c[5]+0.3535533905932737*uy_r[1]+0.3535533905932737*uy_c[1]; 
  avg_uy_r[2] = (-0.6123724356957944*uy_r[6])+0.6123724356957944*uy_c[6]+0.3535533905932737*uy_r[2]+0.3535533905932737*uy_c[2]; 
  avg_uy_r[3] = (-0.6123724356957944*uy_r[7])+0.6123724356957944*uy_c[7]+0.3535533905932737*uy_r[4]+0.3535533905932737*uy_c[4]; 

  avg_uz_l[0] = 0.6123724356957944*uz_l[3]-0.6123724356957944*uz_c[3]+0.3535533905932737*uz_l[0]+0.3535533905932737*uz_c[0]; 
  avg_uz_l[1] = 0.6123724356957944*uz_l[5]-0.6123724356957944*uz_c[5]+0.3535533905932737*uz_l[1]+0.3535533905932737*uz_c[1]; 
  avg_uz_l[2] = 0.6123724356957944*uz_l[6]-0.6123724356957944*uz_c[6]+0.3535533905932737*uz_l[2]+0.3535533905932737*uz_c[2]; 
  avg_uz_l[3] = 0.6123724356957944*uz_l[7]-0.6123724356957944*uz_c[7]+0.3535533905932737*uz_l[4]+0.3535533905932737*uz_c[4]; 

  avg_uz_r[0] = (-0.6123724356957944*uz_r[3])+0.6123724356957944*uz_c[3]+0.3535533905932737*uz_r[0]+0.3535533905932737*uz_c[0]; 
  avg_uz_r[1] = (-0.6123724356957944*uz_r[5])+0.6123724356957944*uz_c[5]+0.3535533905932737*uz_r[1]+0.3535533905932737*uz_c[1]; 
  avg_uz_r[2] = (-0.6123724356957944*uz_r[6])+0.6123724356957944*uz_c[6]+0.3535533905932737*uz_r[2]+0.3535533905932737*uz_c[2]; 
  avg_uz_r[3] = (-0.6123724356957944*uz_r[7])+0.6123724356957944*uz_c[7]+0.3535533905932737*uz_r[4]+0.3535533905932737*uz_c[4]; 

  jump_rhoux_l[0] = (-0.6123724356957944*rhoux_l[3])-0.6123724356957944*rhoux_c[3]-0.3535533905932737*rhoux_l[0]+0.3535533905932737*rhoux_c[0]; 
  jump_rhoux_l[1] = (-0.6123724356957944*rhoux_l[5])-0.6123724356957944*rhoux_c[5]-0.3535533905932737*rhoux_l[1]+0.3535533905932737*rhoux_c[1]; 
  jump_rhoux_l[2] = (-0.6123724356957944*rhoux_l[6])-0.6123724356957944*rhoux_c[6]-0.3535533905932737*rhoux_l[2]+0.3535533905932737*rhoux_c[2]; 
  jump_rhoux_l[3] = (-0.6123724356957944*rhoux_l[7])-0.6123724356957944*rhoux_c[7]-0.3535533905932737*rhoux_l[4]+0.3535533905932737*rhoux_c[4]; 

  jump_rhoux_r[0] = (-0.6123724356957944*rhoux_r[3])-0.6123724356957944*rhoux_c[3]+0.3535533905932737*rhoux_r[0]-0.3535533905932737*rhoux_c[0]; 
  jump_rhoux_r[1] = (-0.6123724356957944*rhoux_r[5])-0.6123724356957944*rhoux_c[5]+0.3535533905932737*rhoux_r[1]-0.3535533905932737*rhoux_c[1]; 
  jump_rhoux_r[2] = (-0.6123724356957944*rhoux_r[6])-0.6123724356957944*rhoux_c[6]+0.3535533905932737*rhoux_r[2]-0.3535533905932737*rhoux_c[2]; 
  jump_rhoux_r[3] = (-0.6123724356957944*rhoux_r[7])-0.6123724356957944*rhoux_c[7]+0.3535533905932737*rhoux_r[4]-0.3535533905932737*rhoux_c[4]; 

  jump_rhouy_l[0] = (-0.6123724356957944*rhouy_l[3])-0.6123724356957944*rhouy_c[3]-0.3535533905932737*rhouy_l[0]+0.3535533905932737*rhouy_c[0]; 
  jump_rhouy_l[1] = (-0.6123724356957944*rhouy_l[5])-0.6123724356957944*rhouy_c[5]-0.3535533905932737*rhouy_l[1]+0.3535533905932737*rhouy_c[1]; 
  jump_rhouy_l[2] = (-0.6123724356957944*rhouy_l[6])-0.6123724356957944*rhouy_c[6]-0.3535533905932737*rhouy_l[2]+0.3535533905932737*rhouy_c[2]; 
  jump_rhouy_l[3] = (-0.6123724356957944*rhouy_l[7])-0.6123724356957944*rhouy_c[7]-0.3535533905932737*rhouy_l[4]+0.3535533905932737*rhouy_c[4]; 

  jump_rhouy_r[0] = (-0.6123724356957944*rhouy_r[3])-0.6123724356957944*rhouy_c[3]+0.3535533905932737*rhouy_r[0]-0.3535533905932737*rhouy_c[0]; 
  jump_rhouy_r[1] = (-0.6123724356957944*rhouy_r[5])-0.6123724356957944*rhouy_c[5]+0.3535533905932737*rhouy_r[1]-0.3535533905932737*rhouy_c[1]; 
  jump_rhouy_r[2] = (-0.6123724356957944*rhouy_r[6])-0.6123724356957944*rhouy_c[6]+0.3535533905932737*rhouy_r[2]-0.3535533905932737*rhouy_c[2]; 
  jump_rhouy_r[3] = (-0.6123724356957944*rhouy_r[7])-0.6123724356957944*rhouy_c[7]+0.3535533905932737*rhouy_r[4]-0.3535533905932737*rhouy_c[4]; 

  jump_rhouz_l[0] = (-0.6123724356957944*rhouz_l[3])-0.6123724356957944*rhouz_c[3]-0.3535533905932737*rhouz_l[0]+0.3535533905932737*rhouz_c[0]; 
  jump_rhouz_l[1] = (-0.6123724356957944*rhouz_l[5])-0.6123724356957944*rhouz_c[5]-0.3535533905932737*rhouz_l[1]+0.3535533905932737*rhouz_c[1]; 
  jump_rhouz_l[2] = (-0.6123724356957944*rhouz_l[6])-0.6123724356957944*rhouz_c[6]-0.3535533905932737*rhouz_l[2]+0.3535533905932737*rhouz_c[2]; 
  jump_rhouz_l[3] = (-0.6123724356957944*rhouz_l[7])-0.6123724356957944*rhouz_c[7]-0.3535533905932737*rhouz_l[4]+0.3535533905932737*rhouz_c[4]; 

  jump_rhouz_r[0] = (-0.6123724356957944*rhouz_r[3])-0.6123724356957944*rhouz_c[3]+0.3535533905932737*rhouz_r[0]-0.3535533905932737*rhouz_c[0]; 
  jump_rhouz_r[1] = (-0.6123724356957944*rhouz_r[5])-0.6123724356957944*rhouz_c[5]+0.3535533905932737*rhouz_r[1]-0.3535533905932737*rhouz_c[1]; 
  jump_rhouz_r[2] = (-0.6123724356957944*rhouz_r[6])-0.6123724356957944*rhouz_c[6]+0.3535533905932737*rhouz_r[2]-0.3535533905932737*rhouz_c[2]; 
  jump_rhouz_r[3] = (-0.6123724356957944*rhouz_r[7])-0.6123724356957944*rhouz_c[7]+0.3535533905932737*rhouz_r[4]-0.3535533905932737*rhouz_c[4]; 

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

  ul_r = ser_3x_p1_surfx3_eval_quad_node_0_r(u_l); 
  uc_l = ser_3x_p1_surfx3_eval_quad_node_0_l(u_c); 
  uc_r = ser_3x_p1_surfx3_eval_quad_node_0_r(u_c); 
  ur_l = ser_3x_p1_surfx3_eval_quad_node_0_l(u_r); 
  uQuad_l = fmax(fabs(ul_r), fabs(uc_l)); 
  uQuad_r = fmax(fabs(uc_r), fabs(ur_l)); 
  vth_sql_r = ser_3x_p1_surfx3_eval_quad_node_0_r(vth_sql); 
  vth_sqc_l = ser_3x_p1_surfx3_eval_quad_node_0_l(vth_sqc); 
  vth_sqc_r = ser_3x_p1_surfx3_eval_quad_node_0_r(vth_sqc); 
  vth_sqr_l = ser_3x_p1_surfx3_eval_quad_node_0_l(vth_sqr); 
  vthQuad_l = fmax(sqrt(fabs(vth_sql_r)), sqrt(fabs(vth_sqc_l))); 
  vthQuad_r = fmax(sqrt(fabs(vth_sqc_r)), sqrt(fabs(vth_sqr_l))); 
  max_speed_quad_l[0] = uQuad_l + vthQuad_l; 
  max_speed_quad_r[0] = uQuad_r + vthQuad_r; 
  ul_r = ser_3x_p1_surfx3_eval_quad_node_1_r(u_l); 
  uc_l = ser_3x_p1_surfx3_eval_quad_node_1_l(u_c); 
  uc_r = ser_3x_p1_surfx3_eval_quad_node_1_r(u_c); 
  ur_l = ser_3x_p1_surfx3_eval_quad_node_1_l(u_r); 
  uQuad_l = fmax(fabs(ul_r), fabs(uc_l)); 
  uQuad_r = fmax(fabs(uc_r), fabs(ur_l)); 
  vth_sql_r = ser_3x_p1_surfx3_eval_quad_node_1_r(vth_sql); 
  vth_sqc_l = ser_3x_p1_surfx3_eval_quad_node_1_l(vth_sqc); 
  vth_sqc_r = ser_3x_p1_surfx3_eval_quad_node_1_r(vth_sqc); 
  vth_sqr_l = ser_3x_p1_surfx3_eval_quad_node_1_l(vth_sqr); 
  vthQuad_l = fmax(sqrt(fabs(vth_sql_r)), sqrt(fabs(vth_sqc_l))); 
  vthQuad_r = fmax(sqrt(fabs(vth_sqc_r)), sqrt(fabs(vth_sqr_l))); 
  max_speed_quad_l[1] = uQuad_l + vthQuad_l; 
  max_speed_quad_r[1] = uQuad_r + vthQuad_r; 
  ul_r = ser_3x_p1_surfx3_eval_quad_node_2_r(u_l); 
  uc_l = ser_3x_p1_surfx3_eval_quad_node_2_l(u_c); 
  uc_r = ser_3x_p1_surfx3_eval_quad_node_2_r(u_c); 
  ur_l = ser_3x_p1_surfx3_eval_quad_node_2_l(u_r); 
  uQuad_l = fmax(fabs(ul_r), fabs(uc_l)); 
  uQuad_r = fmax(fabs(uc_r), fabs(ur_l)); 
  vth_sql_r = ser_3x_p1_surfx3_eval_quad_node_2_r(vth_sql); 
  vth_sqc_l = ser_3x_p1_surfx3_eval_quad_node_2_l(vth_sqc); 
  vth_sqc_r = ser_3x_p1_surfx3_eval_quad_node_2_r(vth_sqc); 
  vth_sqr_l = ser_3x_p1_surfx3_eval_quad_node_2_l(vth_sqr); 
  vthQuad_l = fmax(sqrt(fabs(vth_sql_r)), sqrt(fabs(vth_sqc_l))); 
  vthQuad_r = fmax(sqrt(fabs(vth_sqc_r)), sqrt(fabs(vth_sqr_l))); 
  max_speed_quad_l[2] = uQuad_l + vthQuad_l; 
  max_speed_quad_r[2] = uQuad_r + vthQuad_r; 
  ul_r = ser_3x_p1_surfx3_eval_quad_node_3_r(u_l); 
  uc_l = ser_3x_p1_surfx3_eval_quad_node_3_l(u_c); 
  uc_r = ser_3x_p1_surfx3_eval_quad_node_3_r(u_c); 
  ur_l = ser_3x_p1_surfx3_eval_quad_node_3_l(u_r); 
  uQuad_l = fmax(fabs(ul_r), fabs(uc_l)); 
  uQuad_r = fmax(fabs(uc_r), fabs(ur_l)); 
  vth_sql_r = ser_3x_p1_surfx3_eval_quad_node_3_r(vth_sql); 
  vth_sqc_l = ser_3x_p1_surfx3_eval_quad_node_3_l(vth_sqc); 
  vth_sqc_r = ser_3x_p1_surfx3_eval_quad_node_3_r(vth_sqc); 
  vth_sqr_l = ser_3x_p1_surfx3_eval_quad_node_3_l(vth_sqr); 
  vthQuad_l = fmax(sqrt(fabs(vth_sql_r)), sqrt(fabs(vth_sqc_l))); 
  vthQuad_r = fmax(sqrt(fabs(vth_sqc_r)), sqrt(fabs(vth_sqr_l))); 
  max_speed_quad_l[3] = uQuad_l + vthQuad_l; 
  max_speed_quad_r[3] = uQuad_r + vthQuad_r; 
  ser_3x_p1_upwind_quad_to_modal(max_speed_quad_l, max_speed_modal_l); 
  ser_3x_p1_upwind_quad_to_modal(max_speed_quad_r, max_speed_modal_r); 
  outrhou0[0] += 0.3535533905932737*jump_rhoux_r[3]*max_speed_modal_r[3]*dx1-0.3535533905932737*jump_rhoux_l[3]*max_speed_modal_l[3]*dx1-0.3535533905932737*avg_ux_r[3]*flux_rho_r[3]*dx1+0.3535533905932737*avg_ux_l[3]*flux_rho_l[3]*dx1+0.3535533905932737*jump_rhoux_r[2]*max_speed_modal_r[2]*dx1-0.3535533905932737*jump_rhoux_l[2]*max_speed_modal_l[2]*dx1-0.3535533905932737*avg_ux_r[2]*flux_rho_r[2]*dx1+0.3535533905932737*avg_ux_l[2]*flux_rho_l[2]*dx1+0.3535533905932737*jump_rhoux_r[1]*max_speed_modal_r[1]*dx1-0.3535533905932737*jump_rhoux_l[1]*max_speed_modal_l[1]*dx1-0.3535533905932737*avg_ux_r[1]*flux_rho_r[1]*dx1+0.3535533905932737*avg_ux_l[1]*flux_rho_l[1]*dx1+0.3535533905932737*jump_rhoux_r[0]*max_speed_modal_r[0]*dx1-0.3535533905932737*jump_rhoux_l[0]*max_speed_modal_l[0]*dx1-0.3535533905932737*avg_ux_r[0]*flux_rho_r[0]*dx1+0.3535533905932737*avg_ux_l[0]*flux_rho_l[0]*dx1; 
  outrhou0[1] += 0.3535533905932737*jump_rhoux_r[2]*max_speed_modal_r[3]*dx1-0.3535533905932737*jump_rhoux_l[2]*max_speed_modal_l[3]*dx1+0.3535533905932737*max_speed_modal_r[2]*jump_rhoux_r[3]*dx1-0.3535533905932737*max_speed_modal_l[2]*jump_rhoux_l[3]*dx1-0.3535533905932737*avg_ux_r[2]*flux_rho_r[3]*dx1+0.3535533905932737*avg_ux_l[2]*flux_rho_l[3]*dx1-0.3535533905932737*flux_rho_r[2]*avg_ux_r[3]*dx1+0.3535533905932737*flux_rho_l[2]*avg_ux_l[3]*dx1+0.3535533905932737*jump_rhoux_r[0]*max_speed_modal_r[1]*dx1-0.3535533905932737*jump_rhoux_l[0]*max_speed_modal_l[1]*dx1+0.3535533905932737*max_speed_modal_r[0]*jump_rhoux_r[1]*dx1-0.3535533905932737*max_speed_modal_l[0]*jump_rhoux_l[1]*dx1-0.3535533905932737*avg_ux_r[0]*flux_rho_r[1]*dx1+0.3535533905932737*avg_ux_l[0]*flux_rho_l[1]*dx1-0.3535533905932737*flux_rho_r[0]*avg_ux_r[1]*dx1+0.3535533905932737*flux_rho_l[0]*avg_ux_l[1]*dx1; 
  outrhou0[2] += 0.3535533905932737*jump_rhoux_r[1]*max_speed_modal_r[3]*dx1-0.3535533905932737*jump_rhoux_l[1]*max_speed_modal_l[3]*dx1+0.3535533905932737*max_speed_modal_r[1]*jump_rhoux_r[3]*dx1-0.3535533905932737*max_speed_modal_l[1]*jump_rhoux_l[3]*dx1-0.3535533905932737*avg_ux_r[1]*flux_rho_r[3]*dx1+0.3535533905932737*avg_ux_l[1]*flux_rho_l[3]*dx1-0.3535533905932737*flux_rho_r[1]*avg_ux_r[3]*dx1+0.3535533905932737*flux_rho_l[1]*avg_ux_l[3]*dx1+0.3535533905932737*jump_rhoux_r[0]*max_speed_modal_r[2]*dx1-0.3535533905932737*jump_rhoux_l[0]*max_speed_modal_l[2]*dx1+0.3535533905932737*max_speed_modal_r[0]*jump_rhoux_r[2]*dx1-0.3535533905932737*max_speed_modal_l[0]*jump_rhoux_l[2]*dx1-0.3535533905932737*avg_ux_r[0]*flux_rho_r[2]*dx1+0.3535533905932737*avg_ux_l[0]*flux_rho_l[2]*dx1-0.3535533905932737*flux_rho_r[0]*avg_ux_r[2]*dx1+0.3535533905932737*flux_rho_l[0]*avg_ux_l[2]*dx1; 
  outrhou0[3] += 0.6123724356957944*jump_rhoux_r[3]*max_speed_modal_r[3]*dx1+0.6123724356957944*jump_rhoux_l[3]*max_speed_modal_l[3]*dx1-0.6123724356957944*avg_ux_r[3]*flux_rho_r[3]*dx1-0.6123724356957944*avg_ux_l[3]*flux_rho_l[3]*dx1+0.6123724356957944*jump_rhoux_r[2]*max_speed_modal_r[2]*dx1+0.6123724356957944*jump_rhoux_l[2]*max_speed_modal_l[2]*dx1-0.6123724356957944*avg_ux_r[2]*flux_rho_r[2]*dx1-0.6123724356957944*avg_ux_l[2]*flux_rho_l[2]*dx1+0.6123724356957944*jump_rhoux_r[1]*max_speed_modal_r[1]*dx1+0.6123724356957944*jump_rhoux_l[1]*max_speed_modal_l[1]*dx1-0.6123724356957944*avg_ux_r[1]*flux_rho_r[1]*dx1-0.6123724356957944*avg_ux_l[1]*flux_rho_l[1]*dx1+0.6123724356957944*jump_rhoux_r[0]*max_speed_modal_r[0]*dx1+0.6123724356957944*jump_rhoux_l[0]*max_speed_modal_l[0]*dx1-0.6123724356957944*avg_ux_r[0]*flux_rho_r[0]*dx1-0.6123724356957944*avg_ux_l[0]*flux_rho_l[0]*dx1; 
  outrhou0[4] += 0.3535533905932737*jump_rhoux_r[0]*max_speed_modal_r[3]*dx1-0.3535533905932737*jump_rhoux_l[0]*max_speed_modal_l[3]*dx1+0.3535533905932737*max_speed_modal_r[0]*jump_rhoux_r[3]*dx1-0.3535533905932737*max_speed_modal_l[0]*jump_rhoux_l[3]*dx1-0.3535533905932737*avg_ux_r[0]*flux_rho_r[3]*dx1+0.3535533905932737*avg_ux_l[0]*flux_rho_l[3]*dx1-0.3535533905932737*flux_rho_r[0]*avg_ux_r[3]*dx1+0.3535533905932737*flux_rho_l[0]*avg_ux_l[3]*dx1+0.3535533905932737*jump_rhoux_r[1]*max_speed_modal_r[2]*dx1-0.3535533905932737*jump_rhoux_l[1]*max_speed_modal_l[2]*dx1+0.3535533905932737*max_speed_modal_r[1]*jump_rhoux_r[2]*dx1-0.3535533905932737*max_speed_modal_l[1]*jump_rhoux_l[2]*dx1-0.3535533905932737*avg_ux_r[1]*flux_rho_r[2]*dx1+0.3535533905932737*avg_ux_l[1]*flux_rho_l[2]*dx1-0.3535533905932737*flux_rho_r[1]*avg_ux_r[2]*dx1+0.3535533905932737*flux_rho_l[1]*avg_ux_l[2]*dx1; 
  outrhou0[5] += 0.6123724356957944*jump_rhoux_r[2]*max_speed_modal_r[3]*dx1+0.6123724356957944*jump_rhoux_l[2]*max_speed_modal_l[3]*dx1+0.6123724356957944*max_speed_modal_r[2]*jump_rhoux_r[3]*dx1+0.6123724356957944*max_speed_modal_l[2]*jump_rhoux_l[3]*dx1-0.6123724356957944*avg_ux_r[2]*flux_rho_r[3]*dx1-0.6123724356957944*avg_ux_l[2]*flux_rho_l[3]*dx1-0.6123724356957944*flux_rho_r[2]*avg_ux_r[3]*dx1-0.6123724356957944*flux_rho_l[2]*avg_ux_l[3]*dx1+0.6123724356957944*jump_rhoux_r[0]*max_speed_modal_r[1]*dx1+0.6123724356957944*jump_rhoux_l[0]*max_speed_modal_l[1]*dx1+0.6123724356957944*max_speed_modal_r[0]*jump_rhoux_r[1]*dx1+0.6123724356957944*max_speed_modal_l[0]*jump_rhoux_l[1]*dx1-0.6123724356957944*avg_ux_r[0]*flux_rho_r[1]*dx1-0.6123724356957944*avg_ux_l[0]*flux_rho_l[1]*dx1-0.6123724356957944*flux_rho_r[0]*avg_ux_r[1]*dx1-0.6123724356957944*flux_rho_l[0]*avg_ux_l[1]*dx1; 
  outrhou0[6] += 0.6123724356957944*jump_rhoux_r[1]*max_speed_modal_r[3]*dx1+0.6123724356957944*jump_rhoux_l[1]*max_speed_modal_l[3]*dx1+0.6123724356957944*max_speed_modal_r[1]*jump_rhoux_r[3]*dx1+0.6123724356957944*max_speed_modal_l[1]*jump_rhoux_l[3]*dx1-0.6123724356957944*avg_ux_r[1]*flux_rho_r[3]*dx1-0.6123724356957944*avg_ux_l[1]*flux_rho_l[3]*dx1-0.6123724356957944*flux_rho_r[1]*avg_ux_r[3]*dx1-0.6123724356957944*flux_rho_l[1]*avg_ux_l[3]*dx1+0.6123724356957944*jump_rhoux_r[0]*max_speed_modal_r[2]*dx1+0.6123724356957944*jump_rhoux_l[0]*max_speed_modal_l[2]*dx1+0.6123724356957944*max_speed_modal_r[0]*jump_rhoux_r[2]*dx1+0.6123724356957944*max_speed_modal_l[0]*jump_rhoux_l[2]*dx1-0.6123724356957944*avg_ux_r[0]*flux_rho_r[2]*dx1-0.6123724356957944*avg_ux_l[0]*flux_rho_l[2]*dx1-0.6123724356957944*flux_rho_r[0]*avg_ux_r[2]*dx1-0.6123724356957944*flux_rho_l[0]*avg_ux_l[2]*dx1; 
  outrhou0[7] += 0.6123724356957944*jump_rhoux_r[0]*max_speed_modal_r[3]*dx1+0.6123724356957944*jump_rhoux_l[0]*max_speed_modal_l[3]*dx1+0.6123724356957944*max_speed_modal_r[0]*jump_rhoux_r[3]*dx1+0.6123724356957944*max_speed_modal_l[0]*jump_rhoux_l[3]*dx1-0.6123724356957944*avg_ux_r[0]*flux_rho_r[3]*dx1-0.6123724356957944*avg_ux_l[0]*flux_rho_l[3]*dx1-0.6123724356957944*flux_rho_r[0]*avg_ux_r[3]*dx1-0.6123724356957944*flux_rho_l[0]*avg_ux_l[3]*dx1+0.6123724356957944*jump_rhoux_r[1]*max_speed_modal_r[2]*dx1+0.6123724356957944*jump_rhoux_l[1]*max_speed_modal_l[2]*dx1+0.6123724356957944*max_speed_modal_r[1]*jump_rhoux_r[2]*dx1+0.6123724356957944*max_speed_modal_l[1]*jump_rhoux_l[2]*dx1-0.6123724356957944*avg_ux_r[1]*flux_rho_r[2]*dx1-0.6123724356957944*avg_ux_l[1]*flux_rho_l[2]*dx1-0.6123724356957944*flux_rho_r[1]*avg_ux_r[2]*dx1-0.6123724356957944*flux_rho_l[1]*avg_ux_l[2]*dx1; 

  outrhou1[0] += 0.3535533905932737*jump_rhouy_r[3]*max_speed_modal_r[3]*dx1-0.3535533905932737*jump_rhouy_l[3]*max_speed_modal_l[3]*dx1-0.3535533905932737*avg_uy_r[3]*flux_rho_r[3]*dx1+0.3535533905932737*avg_uy_l[3]*flux_rho_l[3]*dx1+0.3535533905932737*jump_rhouy_r[2]*max_speed_modal_r[2]*dx1-0.3535533905932737*jump_rhouy_l[2]*max_speed_modal_l[2]*dx1-0.3535533905932737*avg_uy_r[2]*flux_rho_r[2]*dx1+0.3535533905932737*avg_uy_l[2]*flux_rho_l[2]*dx1+0.3535533905932737*jump_rhouy_r[1]*max_speed_modal_r[1]*dx1-0.3535533905932737*jump_rhouy_l[1]*max_speed_modal_l[1]*dx1-0.3535533905932737*avg_uy_r[1]*flux_rho_r[1]*dx1+0.3535533905932737*avg_uy_l[1]*flux_rho_l[1]*dx1+0.3535533905932737*jump_rhouy_r[0]*max_speed_modal_r[0]*dx1-0.3535533905932737*jump_rhouy_l[0]*max_speed_modal_l[0]*dx1-0.3535533905932737*avg_uy_r[0]*flux_rho_r[0]*dx1+0.3535533905932737*avg_uy_l[0]*flux_rho_l[0]*dx1; 
  outrhou1[1] += 0.3535533905932737*jump_rhouy_r[2]*max_speed_modal_r[3]*dx1-0.3535533905932737*jump_rhouy_l[2]*max_speed_modal_l[3]*dx1+0.3535533905932737*max_speed_modal_r[2]*jump_rhouy_r[3]*dx1-0.3535533905932737*max_speed_modal_l[2]*jump_rhouy_l[3]*dx1-0.3535533905932737*avg_uy_r[2]*flux_rho_r[3]*dx1+0.3535533905932737*avg_uy_l[2]*flux_rho_l[3]*dx1-0.3535533905932737*flux_rho_r[2]*avg_uy_r[3]*dx1+0.3535533905932737*flux_rho_l[2]*avg_uy_l[3]*dx1+0.3535533905932737*jump_rhouy_r[0]*max_speed_modal_r[1]*dx1-0.3535533905932737*jump_rhouy_l[0]*max_speed_modal_l[1]*dx1+0.3535533905932737*max_speed_modal_r[0]*jump_rhouy_r[1]*dx1-0.3535533905932737*max_speed_modal_l[0]*jump_rhouy_l[1]*dx1-0.3535533905932737*avg_uy_r[0]*flux_rho_r[1]*dx1+0.3535533905932737*avg_uy_l[0]*flux_rho_l[1]*dx1-0.3535533905932737*flux_rho_r[0]*avg_uy_r[1]*dx1+0.3535533905932737*flux_rho_l[0]*avg_uy_l[1]*dx1; 
  outrhou1[2] += 0.3535533905932737*jump_rhouy_r[1]*max_speed_modal_r[3]*dx1-0.3535533905932737*jump_rhouy_l[1]*max_speed_modal_l[3]*dx1+0.3535533905932737*max_speed_modal_r[1]*jump_rhouy_r[3]*dx1-0.3535533905932737*max_speed_modal_l[1]*jump_rhouy_l[3]*dx1-0.3535533905932737*avg_uy_r[1]*flux_rho_r[3]*dx1+0.3535533905932737*avg_uy_l[1]*flux_rho_l[3]*dx1-0.3535533905932737*flux_rho_r[1]*avg_uy_r[3]*dx1+0.3535533905932737*flux_rho_l[1]*avg_uy_l[3]*dx1+0.3535533905932737*jump_rhouy_r[0]*max_speed_modal_r[2]*dx1-0.3535533905932737*jump_rhouy_l[0]*max_speed_modal_l[2]*dx1+0.3535533905932737*max_speed_modal_r[0]*jump_rhouy_r[2]*dx1-0.3535533905932737*max_speed_modal_l[0]*jump_rhouy_l[2]*dx1-0.3535533905932737*avg_uy_r[0]*flux_rho_r[2]*dx1+0.3535533905932737*avg_uy_l[0]*flux_rho_l[2]*dx1-0.3535533905932737*flux_rho_r[0]*avg_uy_r[2]*dx1+0.3535533905932737*flux_rho_l[0]*avg_uy_l[2]*dx1; 
  outrhou1[3] += 0.6123724356957944*jump_rhouy_r[3]*max_speed_modal_r[3]*dx1+0.6123724356957944*jump_rhouy_l[3]*max_speed_modal_l[3]*dx1-0.6123724356957944*avg_uy_r[3]*flux_rho_r[3]*dx1-0.6123724356957944*avg_uy_l[3]*flux_rho_l[3]*dx1+0.6123724356957944*jump_rhouy_r[2]*max_speed_modal_r[2]*dx1+0.6123724356957944*jump_rhouy_l[2]*max_speed_modal_l[2]*dx1-0.6123724356957944*avg_uy_r[2]*flux_rho_r[2]*dx1-0.6123724356957944*avg_uy_l[2]*flux_rho_l[2]*dx1+0.6123724356957944*jump_rhouy_r[1]*max_speed_modal_r[1]*dx1+0.6123724356957944*jump_rhouy_l[1]*max_speed_modal_l[1]*dx1-0.6123724356957944*avg_uy_r[1]*flux_rho_r[1]*dx1-0.6123724356957944*avg_uy_l[1]*flux_rho_l[1]*dx1+0.6123724356957944*jump_rhouy_r[0]*max_speed_modal_r[0]*dx1+0.6123724356957944*jump_rhouy_l[0]*max_speed_modal_l[0]*dx1-0.6123724356957944*avg_uy_r[0]*flux_rho_r[0]*dx1-0.6123724356957944*avg_uy_l[0]*flux_rho_l[0]*dx1; 
  outrhou1[4] += 0.3535533905932737*jump_rhouy_r[0]*max_speed_modal_r[3]*dx1-0.3535533905932737*jump_rhouy_l[0]*max_speed_modal_l[3]*dx1+0.3535533905932737*max_speed_modal_r[0]*jump_rhouy_r[3]*dx1-0.3535533905932737*max_speed_modal_l[0]*jump_rhouy_l[3]*dx1-0.3535533905932737*avg_uy_r[0]*flux_rho_r[3]*dx1+0.3535533905932737*avg_uy_l[0]*flux_rho_l[3]*dx1-0.3535533905932737*flux_rho_r[0]*avg_uy_r[3]*dx1+0.3535533905932737*flux_rho_l[0]*avg_uy_l[3]*dx1+0.3535533905932737*jump_rhouy_r[1]*max_speed_modal_r[2]*dx1-0.3535533905932737*jump_rhouy_l[1]*max_speed_modal_l[2]*dx1+0.3535533905932737*max_speed_modal_r[1]*jump_rhouy_r[2]*dx1-0.3535533905932737*max_speed_modal_l[1]*jump_rhouy_l[2]*dx1-0.3535533905932737*avg_uy_r[1]*flux_rho_r[2]*dx1+0.3535533905932737*avg_uy_l[1]*flux_rho_l[2]*dx1-0.3535533905932737*flux_rho_r[1]*avg_uy_r[2]*dx1+0.3535533905932737*flux_rho_l[1]*avg_uy_l[2]*dx1; 
  outrhou1[5] += 0.6123724356957944*jump_rhouy_r[2]*max_speed_modal_r[3]*dx1+0.6123724356957944*jump_rhouy_l[2]*max_speed_modal_l[3]*dx1+0.6123724356957944*max_speed_modal_r[2]*jump_rhouy_r[3]*dx1+0.6123724356957944*max_speed_modal_l[2]*jump_rhouy_l[3]*dx1-0.6123724356957944*avg_uy_r[2]*flux_rho_r[3]*dx1-0.6123724356957944*avg_uy_l[2]*flux_rho_l[3]*dx1-0.6123724356957944*flux_rho_r[2]*avg_uy_r[3]*dx1-0.6123724356957944*flux_rho_l[2]*avg_uy_l[3]*dx1+0.6123724356957944*jump_rhouy_r[0]*max_speed_modal_r[1]*dx1+0.6123724356957944*jump_rhouy_l[0]*max_speed_modal_l[1]*dx1+0.6123724356957944*max_speed_modal_r[0]*jump_rhouy_r[1]*dx1+0.6123724356957944*max_speed_modal_l[0]*jump_rhouy_l[1]*dx1-0.6123724356957944*avg_uy_r[0]*flux_rho_r[1]*dx1-0.6123724356957944*avg_uy_l[0]*flux_rho_l[1]*dx1-0.6123724356957944*flux_rho_r[0]*avg_uy_r[1]*dx1-0.6123724356957944*flux_rho_l[0]*avg_uy_l[1]*dx1; 
  outrhou1[6] += 0.6123724356957944*jump_rhouy_r[1]*max_speed_modal_r[3]*dx1+0.6123724356957944*jump_rhouy_l[1]*max_speed_modal_l[3]*dx1+0.6123724356957944*max_speed_modal_r[1]*jump_rhouy_r[3]*dx1+0.6123724356957944*max_speed_modal_l[1]*jump_rhouy_l[3]*dx1-0.6123724356957944*avg_uy_r[1]*flux_rho_r[3]*dx1-0.6123724356957944*avg_uy_l[1]*flux_rho_l[3]*dx1-0.6123724356957944*flux_rho_r[1]*avg_uy_r[3]*dx1-0.6123724356957944*flux_rho_l[1]*avg_uy_l[3]*dx1+0.6123724356957944*jump_rhouy_r[0]*max_speed_modal_r[2]*dx1+0.6123724356957944*jump_rhouy_l[0]*max_speed_modal_l[2]*dx1+0.6123724356957944*max_speed_modal_r[0]*jump_rhouy_r[2]*dx1+0.6123724356957944*max_speed_modal_l[0]*jump_rhouy_l[2]*dx1-0.6123724356957944*avg_uy_r[0]*flux_rho_r[2]*dx1-0.6123724356957944*avg_uy_l[0]*flux_rho_l[2]*dx1-0.6123724356957944*flux_rho_r[0]*avg_uy_r[2]*dx1-0.6123724356957944*flux_rho_l[0]*avg_uy_l[2]*dx1; 
  outrhou1[7] += 0.6123724356957944*jump_rhouy_r[0]*max_speed_modal_r[3]*dx1+0.6123724356957944*jump_rhouy_l[0]*max_speed_modal_l[3]*dx1+0.6123724356957944*max_speed_modal_r[0]*jump_rhouy_r[3]*dx1+0.6123724356957944*max_speed_modal_l[0]*jump_rhouy_l[3]*dx1-0.6123724356957944*avg_uy_r[0]*flux_rho_r[3]*dx1-0.6123724356957944*avg_uy_l[0]*flux_rho_l[3]*dx1-0.6123724356957944*flux_rho_r[0]*avg_uy_r[3]*dx1-0.6123724356957944*flux_rho_l[0]*avg_uy_l[3]*dx1+0.6123724356957944*jump_rhouy_r[1]*max_speed_modal_r[2]*dx1+0.6123724356957944*jump_rhouy_l[1]*max_speed_modal_l[2]*dx1+0.6123724356957944*max_speed_modal_r[1]*jump_rhouy_r[2]*dx1+0.6123724356957944*max_speed_modal_l[1]*jump_rhouy_l[2]*dx1-0.6123724356957944*avg_uy_r[1]*flux_rho_r[2]*dx1-0.6123724356957944*avg_uy_l[1]*flux_rho_l[2]*dx1-0.6123724356957944*flux_rho_r[1]*avg_uy_r[2]*dx1-0.6123724356957944*flux_rho_l[1]*avg_uy_l[2]*dx1; 

  outrhou2[0] += 0.3535533905932737*jump_rhouz_r[3]*max_speed_modal_r[3]*dx1-0.3535533905932737*jump_rhouz_l[3]*max_speed_modal_l[3]*dx1-0.3535533905932737*avg_uz_r[3]*flux_rho_r[3]*dx1+0.3535533905932737*avg_uz_l[3]*flux_rho_l[3]*dx1+0.3535533905932737*jump_rhouz_r[2]*max_speed_modal_r[2]*dx1-0.3535533905932737*jump_rhouz_l[2]*max_speed_modal_l[2]*dx1-0.3535533905932737*avg_uz_r[2]*flux_rho_r[2]*dx1+0.3535533905932737*avg_uz_l[2]*flux_rho_l[2]*dx1+0.3535533905932737*jump_rhouz_r[1]*max_speed_modal_r[1]*dx1-0.3535533905932737*jump_rhouz_l[1]*max_speed_modal_l[1]*dx1-0.3535533905932737*avg_uz_r[1]*flux_rho_r[1]*dx1+0.3535533905932737*avg_uz_l[1]*flux_rho_l[1]*dx1+0.3535533905932737*jump_rhouz_r[0]*max_speed_modal_r[0]*dx1-0.3535533905932737*jump_rhouz_l[0]*max_speed_modal_l[0]*dx1-0.3535533905932737*avg_uz_r[0]*flux_rho_r[0]*dx1+0.3535533905932737*avg_uz_l[0]*flux_rho_l[0]*dx1; 
  outrhou2[1] += 0.3535533905932737*jump_rhouz_r[2]*max_speed_modal_r[3]*dx1-0.3535533905932737*jump_rhouz_l[2]*max_speed_modal_l[3]*dx1+0.3535533905932737*max_speed_modal_r[2]*jump_rhouz_r[3]*dx1-0.3535533905932737*max_speed_modal_l[2]*jump_rhouz_l[3]*dx1-0.3535533905932737*avg_uz_r[2]*flux_rho_r[3]*dx1+0.3535533905932737*avg_uz_l[2]*flux_rho_l[3]*dx1-0.3535533905932737*flux_rho_r[2]*avg_uz_r[3]*dx1+0.3535533905932737*flux_rho_l[2]*avg_uz_l[3]*dx1+0.3535533905932737*jump_rhouz_r[0]*max_speed_modal_r[1]*dx1-0.3535533905932737*jump_rhouz_l[0]*max_speed_modal_l[1]*dx1+0.3535533905932737*max_speed_modal_r[0]*jump_rhouz_r[1]*dx1-0.3535533905932737*max_speed_modal_l[0]*jump_rhouz_l[1]*dx1-0.3535533905932737*avg_uz_r[0]*flux_rho_r[1]*dx1+0.3535533905932737*avg_uz_l[0]*flux_rho_l[1]*dx1-0.3535533905932737*flux_rho_r[0]*avg_uz_r[1]*dx1+0.3535533905932737*flux_rho_l[0]*avg_uz_l[1]*dx1; 
  outrhou2[2] += 0.3535533905932737*jump_rhouz_r[1]*max_speed_modal_r[3]*dx1-0.3535533905932737*jump_rhouz_l[1]*max_speed_modal_l[3]*dx1+0.3535533905932737*max_speed_modal_r[1]*jump_rhouz_r[3]*dx1-0.3535533905932737*max_speed_modal_l[1]*jump_rhouz_l[3]*dx1-0.3535533905932737*avg_uz_r[1]*flux_rho_r[3]*dx1+0.3535533905932737*avg_uz_l[1]*flux_rho_l[3]*dx1-0.3535533905932737*flux_rho_r[1]*avg_uz_r[3]*dx1+0.3535533905932737*flux_rho_l[1]*avg_uz_l[3]*dx1+0.3535533905932737*jump_rhouz_r[0]*max_speed_modal_r[2]*dx1-0.3535533905932737*jump_rhouz_l[0]*max_speed_modal_l[2]*dx1+0.3535533905932737*max_speed_modal_r[0]*jump_rhouz_r[2]*dx1-0.3535533905932737*max_speed_modal_l[0]*jump_rhouz_l[2]*dx1-0.3535533905932737*avg_uz_r[0]*flux_rho_r[2]*dx1+0.3535533905932737*avg_uz_l[0]*flux_rho_l[2]*dx1-0.3535533905932737*flux_rho_r[0]*avg_uz_r[2]*dx1+0.3535533905932737*flux_rho_l[0]*avg_uz_l[2]*dx1; 
  outrhou2[3] += 0.6123724356957944*jump_rhouz_r[3]*max_speed_modal_r[3]*dx1+0.6123724356957944*jump_rhouz_l[3]*max_speed_modal_l[3]*dx1-0.6123724356957944*avg_uz_r[3]*flux_rho_r[3]*dx1-0.6123724356957944*avg_uz_l[3]*flux_rho_l[3]*dx1+0.6123724356957944*jump_rhouz_r[2]*max_speed_modal_r[2]*dx1+0.6123724356957944*jump_rhouz_l[2]*max_speed_modal_l[2]*dx1-0.6123724356957944*avg_uz_r[2]*flux_rho_r[2]*dx1-0.6123724356957944*avg_uz_l[2]*flux_rho_l[2]*dx1+0.6123724356957944*jump_rhouz_r[1]*max_speed_modal_r[1]*dx1+0.6123724356957944*jump_rhouz_l[1]*max_speed_modal_l[1]*dx1-0.6123724356957944*avg_uz_r[1]*flux_rho_r[1]*dx1-0.6123724356957944*avg_uz_l[1]*flux_rho_l[1]*dx1+0.6123724356957944*jump_rhouz_r[0]*max_speed_modal_r[0]*dx1+0.6123724356957944*jump_rhouz_l[0]*max_speed_modal_l[0]*dx1-0.6123724356957944*avg_uz_r[0]*flux_rho_r[0]*dx1-0.6123724356957944*avg_uz_l[0]*flux_rho_l[0]*dx1; 
  outrhou2[4] += 0.3535533905932737*jump_rhouz_r[0]*max_speed_modal_r[3]*dx1-0.3535533905932737*jump_rhouz_l[0]*max_speed_modal_l[3]*dx1+0.3535533905932737*max_speed_modal_r[0]*jump_rhouz_r[3]*dx1-0.3535533905932737*max_speed_modal_l[0]*jump_rhouz_l[3]*dx1-0.3535533905932737*avg_uz_r[0]*flux_rho_r[3]*dx1+0.3535533905932737*avg_uz_l[0]*flux_rho_l[3]*dx1-0.3535533905932737*flux_rho_r[0]*avg_uz_r[3]*dx1+0.3535533905932737*flux_rho_l[0]*avg_uz_l[3]*dx1+0.3535533905932737*jump_rhouz_r[1]*max_speed_modal_r[2]*dx1-0.3535533905932737*jump_rhouz_l[1]*max_speed_modal_l[2]*dx1+0.3535533905932737*max_speed_modal_r[1]*jump_rhouz_r[2]*dx1-0.3535533905932737*max_speed_modal_l[1]*jump_rhouz_l[2]*dx1-0.3535533905932737*avg_uz_r[1]*flux_rho_r[2]*dx1+0.3535533905932737*avg_uz_l[1]*flux_rho_l[2]*dx1-0.3535533905932737*flux_rho_r[1]*avg_uz_r[2]*dx1+0.3535533905932737*flux_rho_l[1]*avg_uz_l[2]*dx1; 
  outrhou2[5] += 0.6123724356957944*jump_rhouz_r[2]*max_speed_modal_r[3]*dx1+0.6123724356957944*jump_rhouz_l[2]*max_speed_modal_l[3]*dx1+0.6123724356957944*max_speed_modal_r[2]*jump_rhouz_r[3]*dx1+0.6123724356957944*max_speed_modal_l[2]*jump_rhouz_l[3]*dx1-0.6123724356957944*avg_uz_r[2]*flux_rho_r[3]*dx1-0.6123724356957944*avg_uz_l[2]*flux_rho_l[3]*dx1-0.6123724356957944*flux_rho_r[2]*avg_uz_r[3]*dx1-0.6123724356957944*flux_rho_l[2]*avg_uz_l[3]*dx1+0.6123724356957944*jump_rhouz_r[0]*max_speed_modal_r[1]*dx1+0.6123724356957944*jump_rhouz_l[0]*max_speed_modal_l[1]*dx1+0.6123724356957944*max_speed_modal_r[0]*jump_rhouz_r[1]*dx1+0.6123724356957944*max_speed_modal_l[0]*jump_rhouz_l[1]*dx1-0.6123724356957944*avg_uz_r[0]*flux_rho_r[1]*dx1-0.6123724356957944*avg_uz_l[0]*flux_rho_l[1]*dx1-0.6123724356957944*flux_rho_r[0]*avg_uz_r[1]*dx1-0.6123724356957944*flux_rho_l[0]*avg_uz_l[1]*dx1; 
  outrhou2[6] += 0.6123724356957944*jump_rhouz_r[1]*max_speed_modal_r[3]*dx1+0.6123724356957944*jump_rhouz_l[1]*max_speed_modal_l[3]*dx1+0.6123724356957944*max_speed_modal_r[1]*jump_rhouz_r[3]*dx1+0.6123724356957944*max_speed_modal_l[1]*jump_rhouz_l[3]*dx1-0.6123724356957944*avg_uz_r[1]*flux_rho_r[3]*dx1-0.6123724356957944*avg_uz_l[1]*flux_rho_l[3]*dx1-0.6123724356957944*flux_rho_r[1]*avg_uz_r[3]*dx1-0.6123724356957944*flux_rho_l[1]*avg_uz_l[3]*dx1+0.6123724356957944*jump_rhouz_r[0]*max_speed_modal_r[2]*dx1+0.6123724356957944*jump_rhouz_l[0]*max_speed_modal_l[2]*dx1+0.6123724356957944*max_speed_modal_r[0]*jump_rhouz_r[2]*dx1+0.6123724356957944*max_speed_modal_l[0]*jump_rhouz_l[2]*dx1-0.6123724356957944*avg_uz_r[0]*flux_rho_r[2]*dx1-0.6123724356957944*avg_uz_l[0]*flux_rho_l[2]*dx1-0.6123724356957944*flux_rho_r[0]*avg_uz_r[2]*dx1-0.6123724356957944*flux_rho_l[0]*avg_uz_l[2]*dx1; 
  outrhou2[7] += 0.6123724356957944*jump_rhouz_r[0]*max_speed_modal_r[3]*dx1+0.6123724356957944*jump_rhouz_l[0]*max_speed_modal_l[3]*dx1+0.6123724356957944*max_speed_modal_r[0]*jump_rhouz_r[3]*dx1+0.6123724356957944*max_speed_modal_l[0]*jump_rhouz_l[3]*dx1-0.6123724356957944*avg_uz_r[0]*flux_rho_r[3]*dx1-0.6123724356957944*avg_uz_l[0]*flux_rho_l[3]*dx1-0.6123724356957944*flux_rho_r[0]*avg_uz_r[3]*dx1-0.6123724356957944*flux_rho_l[0]*avg_uz_l[3]*dx1+0.6123724356957944*jump_rhouz_r[1]*max_speed_modal_r[2]*dx1+0.6123724356957944*jump_rhouz_l[1]*max_speed_modal_l[2]*dx1+0.6123724356957944*max_speed_modal_r[1]*jump_rhouz_r[2]*dx1+0.6123724356957944*max_speed_modal_l[1]*jump_rhouz_l[2]*dx1-0.6123724356957944*avg_uz_r[1]*flux_rho_r[2]*dx1-0.6123724356957944*avg_uz_l[1]*flux_rho_l[2]*dx1-0.6123724356957944*flux_rho_r[1]*avg_uz_r[2]*dx1-0.6123724356957944*flux_rho_l[1]*avg_uz_l[2]*dx1; 

  return 0.;

} 
