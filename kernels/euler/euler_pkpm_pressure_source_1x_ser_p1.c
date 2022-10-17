#include <gkyl_euler_kernels.h> 
GKYL_CU_DH void euler_pkpm_pressure_source_1x_ser_p1(const double *qmem, const double *nu, const double *nu_vth_sq, const double *vlasov_pkpm_moms, const double *u_i, const double *euler_pkpm, double* GKYL_RESTRICT out) 
{ 
  // qmem: q/m*EM fields.
  // nu: Collisionality.
  // nu_vth_sq: nu*vth^2 = nu*T/m.
  // vlasov_pkpm_moms: [rho, p_parallel, q_parallel], Moments computed from kinetic equation in pkpm model.
  // u_i: Input flow velocity.
  // euler_pkpm: Input fluid variables.
  // out: Output increment
  const double *rho = &vlasov_pkpm_moms[0]; 
  const double *Ex = &qmem[0]; 
  const double *Ey = &qmem[2]; 
  const double *Ez = &qmem[4]; 
  const double *p_perp = &euler_pkpm[6]; 

  double *outrhoux = &out[0]; 
  double *outrhouy = &out[2]; 
  double *outrhouz = &out[4]; 
  double *outp_perp = &out[6]; 

  outrhoux[0] += 0.7071067811865475*Ex[1]*rho[1]+0.7071067811865475*Ex[0]*rho[0]; 
  outrhoux[1] += 0.7071067811865475*Ex[0]*rho[1]+0.7071067811865475*rho[0]*Ex[1]; 

  outrhouy[0] += 0.7071067811865475*Ey[1]*rho[1]+0.7071067811865475*Ey[0]*rho[0]; 
  outrhouy[1] += 0.7071067811865475*Ey[0]*rho[1]+0.7071067811865475*rho[0]*Ey[1]; 

  outrhouz[0] += 0.7071067811865475*Ez[1]*rho[1]+0.7071067811865475*Ez[0]*rho[0]; 
  outrhouz[1] += 0.7071067811865475*Ez[0]*rho[1]+0.7071067811865475*rho[0]*Ez[1]; 

  outp_perp[0] += 0.7071067811865475*nu_vth_sq[1]*rho[1]-0.7071067811865475*nu[1]*p_perp[1]+0.7071067811865475*nu_vth_sq[0]*rho[0]-0.7071067811865475*nu[0]*p_perp[0]; 
  outp_perp[1] += 0.7071067811865475*nu_vth_sq[0]*rho[1]-0.7071067811865475*nu[0]*p_perp[1]+0.7071067811865475*rho[0]*nu_vth_sq[1]-0.7071067811865475*p_perp[0]*nu[1]; 

} 
