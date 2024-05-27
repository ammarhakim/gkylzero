#include <gkyl_euler_pkpm_kernels.h> 
#include <gkyl_basis_tensor_1x_2p_exp_sq.h> 
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

  // Calculate u^2. 
  double ux_sq[3] = {0.0}; 
  double uy_sq[3] = {0.0}; 
  double uz_sq[3] = {0.0}; 
  tensor_1x_2p_exp_sq(ux, ux_sq); 
  tensor_1x_2p_exp_sq(uy, uy_sq); 
  tensor_1x_2p_exp_sq(uz, uz_sq); 
 
  // Order of integrated variables is (rho, rhoux, rhouy, rhouz, rho ux^2, rho uy^2, rho uz^2, p_parallel, p_perp) 
  int_pkpm_vars[0] += 1.414213562373095*rho[0]; 
  int_pkpm_vars[1] += rho[1]*ux[1]+rho[0]*ux[0]; 
  int_pkpm_vars[2] += rho[1]*uy[1]+rho[0]*uy[0]; 
  int_pkpm_vars[3] += rho[1]*uz[1]+rho[0]*uz[0]; 
  int_pkpm_vars[4] += rho[2]*ux_sq[2]+rho[1]*ux_sq[1]+rho[0]*ux_sq[0]; 
  int_pkpm_vars[5] += rho[2]*uy_sq[2]+rho[1]*uy_sq[1]+rho[0]*uy_sq[0]; 
  int_pkpm_vars[6] += rho[2]*uz_sq[2]+rho[1]*uz_sq[1]+rho[0]*uz_sq[0]; 
  int_pkpm_vars[7] += 1.414213562373095*p_parallel[0]; 
  int_pkpm_vars[8] += 1.414213562373095*p_perp[0]; 
} 
