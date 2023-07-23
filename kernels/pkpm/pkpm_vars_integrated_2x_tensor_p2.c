#include <gkyl_euler_pkpm_kernels.h> 
GKYL_CU_DH void pkpm_vars_integrated_2x_tensor_p2(const double *vlasov_pkpm_moms, const double *statevec, 
  const double* pkpm_prim, double* int_pkpm_vars) 
{ 
  // vlasov_pkpm_moms: [rho, p_parallel, p_perp], Moments computed from kinetic equation in pkpm model.
  // statevec:         [rho ux, rho uy, rho uz], Fluid input state vector.
  // pkpm_prim:        [ux, uy, uz, pxx, pxy, pxz, pyy, pyz, pzz, T_perp/m = p_perp/rho, m/T_perp = rho/p_perp].
  // int_pkpm_vars:    Output integrated variables.

  const double *rho = &vlasov_pkpm_moms[0]; 
  const double *p_parallel = &vlasov_pkpm_moms[9]; 
  const double *p_perp = &vlasov_pkpm_moms[18]; 

  const double *rhoux = &statevec[0]; 
  const double *rhouy = &statevec[9]; 
  const double *rhouz = &statevec[18]; 

  const double *ux = &pkpm_prim[0]; 
  const double *uy = &pkpm_prim[9]; 
  const double *uz = &pkpm_prim[18]; 

  double int_p_par = 2.0*p_parallel[0]; 
  double int_p_perp = 2.0*p_perp[0]; 
  double int_rhoux2 = rhoux[8]*ux[8]+rhoux[7]*ux[7]+rhoux[6]*ux[6]+rhoux[5]*ux[5]+rhoux[4]*ux[4]+rhoux[3]*ux[3]+rhoux[2]*ux[2]+rhoux[1]*ux[1]+rhoux[0]*ux[0]; 
  double int_rhouy2 = rhouy[8]*uy[8]+rhouy[7]*uy[7]+rhouy[6]*uy[6]+rhouy[5]*uy[5]+rhouy[4]*uy[4]+rhouy[3]*uy[3]+rhouy[2]*uy[2]+rhouy[1]*uy[1]+rhouy[0]*uy[0]; 
  double int_rhouz2 = rhouz[8]*uz[8]+rhouz[7]*uz[7]+rhouz[6]*uz[6]+rhouz[5]*uz[5]+rhouz[4]*uz[4]+rhouz[3]*uz[3]+rhouz[2]*uz[2]+rhouz[1]*uz[1]+rhouz[0]*uz[0]; 
  // Order of integrated variables is (rho, p_parallel, p_perp, rho ux^2, rho uy^2, rho uz^2) 
  int_pkpm_vars[0] += 2.0*rho[0]; 
  int_pkpm_vars[1] += int_p_par; 
  int_pkpm_vars[2] += int_p_perp; 
  int_pkpm_vars[3] += int_rhoux2; 
  int_pkpm_vars[4] += int_rhouy2; 
  int_pkpm_vars[5] += int_rhouz2; 
} 
