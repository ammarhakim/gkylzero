#include <gkyl_euler_pkpm_kernels.h> 
GKYL_CU_DH void pkpm_vars_integrated_2x_tensor_p2(const double *vlasov_pkpm_moms, 
  const double* pkpm_u, double* GKYL_RESTRICT int_pkpm_vars) 
{ 
  // vlasov_pkpm_moms: [rho, p_parallel, p_perp], Moments computed from kinetic equation in pkpm model.
  // pkpm_u:           [ux, uy, uz], Fluid input state vector.
  // int_pkpm_vars:    Output integrated variables.

  const double *rho = &vlasov_pkpm_moms[0]; 
  const double *p_parallel = &vlasov_pkpm_moms[9]; 
  const double *p_perp = &vlasov_pkpm_moms[18]; 

  const double *ux = &pkpm_u[0]; 
  const double *uy = &pkpm_u[4]; 
  const double *uz = &pkpm_u[8]; 

  double rhoux[4] = {0.0}; 
  double rhouy[4] = {0.0}; 
  double rhouz[4] = {0.0}; 
  rhoux[0] = 0.5*rho[3]*ux[3]+0.5*rho[2]*ux[2]+0.5*rho[1]*ux[1]+0.5*rho[0]*ux[0]; 
  rhoux[1] = 0.447213595499958*ux[3]*rho[6]+0.4472135954999579*ux[1]*rho[4]+0.5*rho[2]*ux[3]+0.5*ux[2]*rho[3]+0.5*rho[0]*ux[1]+0.5*ux[0]*rho[1]; 
  rhoux[2] = 0.447213595499958*ux[3]*rho[7]+0.4472135954999579*ux[2]*rho[5]+0.5*rho[1]*ux[3]+0.5*ux[1]*rho[3]+0.5*rho[0]*ux[2]+0.5*ux[0]*rho[2]; 
  rhoux[3] = 0.4*ux[3]*rho[8]+0.447213595499958*ux[2]*rho[7]+0.447213595499958*ux[1]*rho[6]+0.4472135954999579*ux[3]*rho[5]+0.4472135954999579*ux[3]*rho[4]+0.5*rho[0]*ux[3]+0.5*ux[0]*rho[3]+0.5*rho[1]*ux[2]+0.5*ux[1]*rho[2]; 
  rhouy[0] = 0.5*rho[3]*uy[3]+0.5*rho[2]*uy[2]+0.5*rho[1]*uy[1]+0.5*rho[0]*uy[0]; 
  rhouy[1] = 0.447213595499958*uy[3]*rho[6]+0.4472135954999579*uy[1]*rho[4]+0.5*rho[2]*uy[3]+0.5*uy[2]*rho[3]+0.5*rho[0]*uy[1]+0.5*uy[0]*rho[1]; 
  rhouy[2] = 0.447213595499958*uy[3]*rho[7]+0.4472135954999579*uy[2]*rho[5]+0.5*rho[1]*uy[3]+0.5*uy[1]*rho[3]+0.5*rho[0]*uy[2]+0.5*uy[0]*rho[2]; 
  rhouy[3] = 0.4*uy[3]*rho[8]+0.447213595499958*uy[2]*rho[7]+0.447213595499958*uy[1]*rho[6]+0.4472135954999579*uy[3]*rho[5]+0.4472135954999579*uy[3]*rho[4]+0.5*rho[0]*uy[3]+0.5*uy[0]*rho[3]+0.5*rho[1]*uy[2]+0.5*uy[1]*rho[2]; 
  rhouz[0] = 0.5*rho[3]*uz[3]+0.5*rho[2]*uz[2]+0.5*rho[1]*uz[1]+0.5*rho[0]*uz[0]; 
  rhouz[1] = 0.447213595499958*uz[3]*rho[6]+0.4472135954999579*uz[1]*rho[4]+0.5*rho[2]*uz[3]+0.5*uz[2]*rho[3]+0.5*rho[0]*uz[1]+0.5*uz[0]*rho[1]; 
  rhouz[2] = 0.447213595499958*uz[3]*rho[7]+0.4472135954999579*uz[2]*rho[5]+0.5*rho[1]*uz[3]+0.5*uz[1]*rho[3]+0.5*rho[0]*uz[2]+0.5*uz[0]*rho[2]; 
  rhouz[3] = 0.4*uz[3]*rho[8]+0.447213595499958*uz[2]*rho[7]+0.447213595499958*uz[1]*rho[6]+0.4472135954999579*uz[3]*rho[5]+0.4472135954999579*uz[3]*rho[4]+0.5*rho[0]*uz[3]+0.5*uz[0]*rho[3]+0.5*rho[1]*uz[2]+0.5*uz[1]*rho[2]; 

  // Order of integrated variables is (rho, rhoux, rhouy, rhouz, rho ux^2, rho uy^2, rho uz^2, p_parallel, p_perp) 
  int_pkpm_vars[0] += 2.0*rho[0]; 
  int_pkpm_vars[1] += 2.0*rhoux[0]; 
  int_pkpm_vars[2] += 2.0*rhouy[0]; 
  int_pkpm_vars[3] += 2.0*rhouz[0]; 
  int_pkpm_vars[4] += rhoux[3]*ux[3]+rhoux[2]*ux[2]+rhoux[1]*ux[1]+rhoux[0]*ux[0]; 
  int_pkpm_vars[5] += rhouy[3]*uy[3]+rhouy[2]*uy[2]+rhouy[1]*uy[1]+rhouy[0]*uy[0]; 
  int_pkpm_vars[6] += rhouz[3]*uz[3]+rhouz[2]*uz[2]+rhouz[1]*uz[1]+rhouz[0]*uz[0]; 
  int_pkpm_vars[7] += 2.0*p_parallel[0]; 
  int_pkpm_vars[8] += 2.0*p_perp[0]; 
} 
