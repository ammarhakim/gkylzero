#include <gkyl_euler_pkpm_kernels.h> 
GKYL_CU_DH void pkpm_vars_integrated_1x_tensor_p2(const double *vlasov_pkpm_moms, 
  const double* pkpm_u, double* GKYL_RESTRICT int_pkpm_vars) 
{ 
  // vlasov_pkpm_moms: [rho, p_parallel, p_perp], Moments computed from kinetic equation in pkpm model.
  // pkpm_u:           [ux, uy, uz], Fluid input state vector.
  // int_pkpm_vars:    Output integrated variables.

  const double *rho = &vlasov_pkpm_moms[0]; 
  const double *p_parallel = &vlasov_pkpm_moms[3]; 
  const double *p_perp = &vlasov_pkpm_moms[6]; 

  const double *ux = &pkpm_u[0]; 
  const double *uy = &pkpm_u[2]; 
  const double *uz = &pkpm_u[4]; 

  double rhoux[2] = {0.0}; 
  double rhouy[2] = {0.0}; 
  double rhouz[2] = {0.0}; 
  rhoux[0] = 0.7071067811865475*rho[1]*ux[1]+0.7071067811865475*rho[0]*ux[0]; 
  rhoux[1] = 0.6324555320336759*ux[1]*rho[2]+0.7071067811865475*rho[0]*ux[1]+0.7071067811865475*ux[0]*rho[1]; 
  rhouy[0] = 0.7071067811865475*rho[1]*uy[1]+0.7071067811865475*rho[0]*uy[0]; 
  rhouy[1] = 0.6324555320336759*uy[1]*rho[2]+0.7071067811865475*rho[0]*uy[1]+0.7071067811865475*uy[0]*rho[1]; 
  rhouz[0] = 0.7071067811865475*rho[1]*uz[1]+0.7071067811865475*rho[0]*uz[0]; 
  rhouz[1] = 0.6324555320336759*uz[1]*rho[2]+0.7071067811865475*rho[0]*uz[1]+0.7071067811865475*uz[0]*rho[1]; 

  // Order of integrated variables is (rho, rhoux, rhouy, rhouz, rho ux^2, rho uy^2, rho uz^2, p_parallel, p_perp) 
  int_pkpm_vars[0] += 1.414213562373095*rho[0]; 
  int_pkpm_vars[1] += 1.414213562373095*rhoux[0]; 
  int_pkpm_vars[2] += 1.414213562373095*rhouy[0]; 
  int_pkpm_vars[3] += 1.414213562373095*rhouz[0]; 
  int_pkpm_vars[4] += rhoux[1]*ux[1]+rhoux[0]*ux[0]; 
  int_pkpm_vars[5] += rhouy[1]*uy[1]+rhouy[0]*uy[0]; 
  int_pkpm_vars[6] += rhouz[1]*uz[1]+rhouz[0]*uz[0]; 
  int_pkpm_vars[7] += 1.414213562373095*p_parallel[0]; 
  int_pkpm_vars[8] += 1.414213562373095*p_perp[0]; 
} 
