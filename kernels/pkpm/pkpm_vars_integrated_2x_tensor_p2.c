#include <gkyl_euler_pkpm_kernels.h> 
#include <gkyl_basis_tensor_2x_2p_exp_sq.h> 
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

  // Calculate u^2. 
  double ux_sq[9] = {0.0}; 
  double uy_sq[9] = {0.0}; 
  double uz_sq[9] = {0.0}; 
  tensor_2x_2p_exp_sq(ux, ux_sq); 
  tensor_2x_2p_exp_sq(uy, uy_sq); 
  tensor_2x_2p_exp_sq(uz, uz_sq); 
 
  // Order of integrated variables is (rho, rhoux, rhouy, rhouz, rho ux^2, rho uy^2, rho uz^2, p_parallel, p_perp) 
  int_pkpm_vars[0] += 2.0*rho[0]; 
  int_pkpm_vars[1] += rho[3]*ux[3]+rho[2]*ux[2]+rho[1]*ux[1]+rho[0]*ux[0]; 
  int_pkpm_vars[2] += rho[3]*uy[3]+rho[2]*uy[2]+rho[1]*uy[1]+rho[0]*uy[0]; 
  int_pkpm_vars[3] += rho[3]*uz[3]+rho[2]*uz[2]+rho[1]*uz[1]+rho[0]*uz[0]; 
  int_pkpm_vars[4] += rho[8]*ux_sq[8]+rho[7]*ux_sq[7]+rho[6]*ux_sq[6]+rho[5]*ux_sq[5]+rho[4]*ux_sq[4]+rho[3]*ux_sq[3]+rho[2]*ux_sq[2]+rho[1]*ux_sq[1]+rho[0]*ux_sq[0]; 
  int_pkpm_vars[5] += rho[8]*uy_sq[8]+rho[7]*uy_sq[7]+rho[6]*uy_sq[6]+rho[5]*uy_sq[5]+rho[4]*uy_sq[4]+rho[3]*uy_sq[3]+rho[2]*uy_sq[2]+rho[1]*uy_sq[1]+rho[0]*uy_sq[0]; 
  int_pkpm_vars[6] += rho[8]*uz_sq[8]+rho[7]*uz_sq[7]+rho[6]*uz_sq[6]+rho[5]*uz_sq[5]+rho[4]*uz_sq[4]+rho[3]*uz_sq[3]+rho[2]*uz_sq[2]+rho[1]*uz_sq[1]+rho[0]*uz_sq[0]; 
  int_pkpm_vars[7] += 2.0*p_parallel[0]; 
  int_pkpm_vars[8] += 2.0*p_perp[0]; 
} 
