#include <gkyl_euler_pkpm_kernels.h> 
GKYL_CU_DH void pkpm_vars_integrated_1x_ser_p1(const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  const double* prim, double* GKYL_RESTRICT int_pkpm_vars) 
{ 
  // vlasov_pkpm_moms: [rho, p_parallel, p_perp], Moments computed from kinetic equation in pkpm model.
  // euler_pkpm:       [rho ux, rho uy, rho uz], Fluid input state vector.
  // prim:             [ux, uy, uz, 1/rho div(p_par b), T_perp/m, m/T_perp].
  // int_pkpm_vars:    Output integrated variables.

  const double *rho = &vlasov_pkpm_moms[0]; 
  const double *p_parallel = &vlasov_pkpm_moms[2]; 
  const double *p_perp = &vlasov_pkpm_moms[4]; 

  const double *rhoux = &euler_pkpm[0]; 
  const double *rhouy = &euler_pkpm[2]; 
  const double *rhouz = &euler_pkpm[4]; 

  const double *ux = &prim[0]; 
  const double *uy = &prim[2]; 
  const double *uz = &prim[4]; 

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
