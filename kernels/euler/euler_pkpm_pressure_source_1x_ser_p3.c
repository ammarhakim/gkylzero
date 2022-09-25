#include <gkyl_euler_kernels.h> 
GKYL_CU_DH void euler_pkpm_pressure_source_1x_ser_p3(const double *nu, const double *nu_vth_sq, const double *vlasov_pkpm_moms, const double *u_i, const double *euler_pkpm, double* GKYL_RESTRICT out) 
{ 
  // nu: Collisionality.
  // nu_vth_sq: nu*vth^2 = nu*T/m.
  // vlasov_pkpm_moms: [rho, p_parallel, q_parallel], Moments computed from kinetic equation in pkpm model.
  // u_i: Input flow velocity.
  // euler_pkpm: Input fluid variables.
  // out: Output increment
  const double *rho = &vlasov_pkpm_moms[0]; 
  const double *p_perp = &euler_pkpm[12]; 

  double *outp_perp = &out[12]; 

  outp_perp[0] += 0.7071067811865475*nu_vth_sq[3]*rho[3]-0.7071067811865475*nu[3]*p_perp[3]+0.7071067811865475*nu_vth_sq[2]*rho[2]-0.7071067811865475*nu[2]*p_perp[2]+0.7071067811865475*nu_vth_sq[1]*rho[1]-0.7071067811865475*nu[1]*p_perp[1]+0.7071067811865475*nu_vth_sq[0]*rho[0]-0.7071067811865475*nu[0]*p_perp[0]; 
  outp_perp[1] += 0.6210590034081186*nu_vth_sq[2]*rho[3]-0.6210590034081186*nu[2]*p_perp[3]+0.6210590034081186*rho[2]*nu_vth_sq[3]-0.6210590034081186*p_perp[2]*nu[3]+0.6324555320336759*nu_vth_sq[1]*rho[2]-0.6324555320336759*nu[1]*p_perp[2]+0.6324555320336759*rho[1]*nu_vth_sq[2]-0.6324555320336759*p_perp[1]*nu[2]+0.7071067811865475*nu_vth_sq[0]*rho[1]-0.7071067811865475*nu[0]*p_perp[1]+0.7071067811865475*rho[0]*nu_vth_sq[1]-0.7071067811865475*p_perp[0]*nu[1]; 
  outp_perp[2] += 0.421637021355784*nu_vth_sq[3]*rho[3]+0.6210590034081186*nu_vth_sq[1]*rho[3]-0.421637021355784*nu[3]*p_perp[3]-0.6210590034081186*nu[1]*p_perp[3]+0.6210590034081186*rho[1]*nu_vth_sq[3]-0.6210590034081186*p_perp[1]*nu[3]+0.4517539514526256*nu_vth_sq[2]*rho[2]+0.7071067811865475*nu_vth_sq[0]*rho[2]-0.4517539514526256*nu[2]*p_perp[2]-0.7071067811865475*nu[0]*p_perp[2]+0.7071067811865475*rho[0]*nu_vth_sq[2]-0.7071067811865475*p_perp[0]*nu[2]+0.6324555320336759*nu_vth_sq[1]*rho[1]-0.6324555320336759*nu[1]*p_perp[1]; 
  outp_perp[3] += 0.421637021355784*nu_vth_sq[2]*rho[3]+0.7071067811865475*nu_vth_sq[0]*rho[3]-0.421637021355784*nu[2]*p_perp[3]-0.7071067811865475*nu[0]*p_perp[3]+0.421637021355784*rho[2]*nu_vth_sq[3]+0.7071067811865475*rho[0]*nu_vth_sq[3]-0.421637021355784*p_perp[2]*nu[3]-0.7071067811865475*p_perp[0]*nu[3]+0.6210590034081186*nu_vth_sq[1]*rho[2]-0.6210590034081186*nu[1]*p_perp[2]+0.6210590034081186*rho[1]*nu_vth_sq[2]-0.6210590034081186*p_perp[1]*nu[2]; 
} 
