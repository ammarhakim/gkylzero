#include <gkyl_euler_pkpm_kernels.h> 
GKYL_CU_DH void euler_pkpm_int_vars_1x_ser_p2(const double *vlasov_pkpm_moms, const double *statevec, 
  const double* u_i, double* int_pkpm_vars) 
{ 
  // vlasov_pkpm_moms:  [rho, p_parallel, p_perp], Moments computed from kinetic equation in pkpm model.
  // statevec:          [rho ux, rho uy, rho uz], Fluid input state vector.
  // u_i:               [ux, uy, uz] Flow velocity.
  // int_pkpm_vars:     Output integrated variables.

  const double *rho = &vlasov_pkpm_moms[0]; 
  const double *p_parallel = &vlasov_pkpm_moms[3]; 
  const double *p_perp = &vlasov_pkpm_moms[6]; 

  const double *rhoux = &statevec[0]; 
  const double *rhouy = &statevec[3]; 
  const double *rhouz = &statevec[6]; 

  const double *ux = &u_i[0]; 
  const double *uy = &u_i[3]; 
  const double *uz = &u_i[6]; 

  double int_p_par = 1.414213562373095*p_parallel[0]; 
  double int_p_perp = 1.414213562373095*p_perp[0]; 
  double int_rhoux2 = rhoux[2]*ux[2]+rhoux[1]*ux[1]+rhoux[0]*ux[0]; 
  double int_rhouy2 = rhouy[2]*uy[2]+rhouy[1]*uy[1]+rhouy[0]*uy[0]; 
  double int_rhouz2 = rhouz[2]*uz[2]+rhouz[1]*uz[1]+rhouz[0]*uz[0]; 
  // Order of integrated variables is (rho, p_parallel, p_perp, rho ux^2, rho uy^2, rho uz^2, 3/2 p, 1/2 rho u^2, E = 3/2 p + 1/2 rho u^2) 
  int_pkpm_vars[0] += 1.414213562373095*rho[0]; 
  int_pkpm_vars[1] += int_p_par; 
  int_pkpm_vars[2] += int_p_perp; 
  int_pkpm_vars[3] += int_rhoux2; 
  int_pkpm_vars[4] += int_rhouy2; 
  int_pkpm_vars[5] += int_rhouz2; 
  int_pkpm_vars[6] += 0.5*int_p_par + int_p_perp; 
  int_pkpm_vars[7] += 0.5*(int_rhoux2 + int_rhouy2 + int_rhouz2); 
  int_pkpm_vars[8] += 0.5*(int_p_par + 2.0*int_p_perp + int_rhoux2 + int_rhouy2 + int_rhouz2); 
} 
