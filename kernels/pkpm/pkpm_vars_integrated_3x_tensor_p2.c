#include <gkyl_euler_pkpm_kernels.h> 
#include <gkyl_basis_tensor_3x_2p_exp_sq.h> 
GKYL_CU_DH void pkpm_vars_integrated_3x_tensor_p2(const double *vlasov_pkpm_moms, 
  const double* pkpm_u, double* GKYL_RESTRICT int_pkpm_vars) 
{ 
  // vlasov_pkpm_moms: [rho, p_parallel, p_perp], Moments computed from kinetic equation in pkpm model.
  // pkpm_u:           [ux, uy, uz], Fluid input state vector.
  // int_pkpm_vars:    Output integrated variables.

  const double *rho = &vlasov_pkpm_moms[0]; 
  const double *p_parallel = &vlasov_pkpm_moms[27]; 
  const double *p_perp = &vlasov_pkpm_moms[54]; 

  const double *ux = &pkpm_u[0]; 
  const double *uy = &pkpm_u[8]; 
  const double *uz = &pkpm_u[16]; 

  // Calculate u^2. 
  double ux_sq[27] = {0.0}; 
  double uy_sq[27] = {0.0}; 
  double uz_sq[27] = {0.0}; 
  tensor_3x_2p_exp_sq(ux, ux_sq); 
  tensor_3x_2p_exp_sq(uy, uy_sq); 
  tensor_3x_2p_exp_sq(uz, uz_sq); 
 
  // Order of integrated variables is (rho, rhoux, rhouy, rhouz, rho ux^2, rho uy^2, rho uz^2, p_parallel, p_perp) 
  int_pkpm_vars[0] += 2.828427124746191*rho[0]; 
  int_pkpm_vars[1] += ux[7]*rho[10]+rho[6]*ux[6]+rho[5]*ux[5]+rho[4]*ux[4]+rho[3]*ux[3]+rho[2]*ux[2]+rho[1]*ux[1]+rho[0]*ux[0]; 
  int_pkpm_vars[2] += uy[7]*rho[10]+rho[6]*uy[6]+rho[5]*uy[5]+rho[4]*uy[4]+rho[3]*uy[3]+rho[2]*uy[2]+rho[1]*uy[1]+rho[0]*uy[0]; 
  int_pkpm_vars[3] += uz[7]*rho[10]+rho[6]*uz[6]+rho[5]*uz[5]+rho[4]*uz[4]+rho[3]*uz[3]+rho[2]*uz[2]+rho[1]*uz[1]+rho[0]*uz[0]; 
  int_pkpm_vars[4] += rho[26]*ux_sq[26]+rho[25]*ux_sq[25]+rho[24]*ux_sq[24]+rho[23]*ux_sq[23]+rho[22]*ux_sq[22]+rho[21]*ux_sq[21]+rho[20]*ux_sq[20]+rho[19]*ux_sq[19]+rho[18]*ux_sq[18]+rho[17]*ux_sq[17]+rho[16]*ux_sq[16]+rho[15]*ux_sq[15]+rho[14]*ux_sq[14]+rho[13]*ux_sq[13]+rho[12]*ux_sq[12]+rho[11]*ux_sq[11]+rho[10]*ux_sq[10]+rho[9]*ux_sq[9]+rho[8]*ux_sq[8]+rho[7]*ux_sq[7]+rho[6]*ux_sq[6]+rho[5]*ux_sq[5]+rho[4]*ux_sq[4]+rho[3]*ux_sq[3]+rho[2]*ux_sq[2]+rho[1]*ux_sq[1]+rho[0]*ux_sq[0]; 
  int_pkpm_vars[5] += rho[26]*uy_sq[26]+rho[25]*uy_sq[25]+rho[24]*uy_sq[24]+rho[23]*uy_sq[23]+rho[22]*uy_sq[22]+rho[21]*uy_sq[21]+rho[20]*uy_sq[20]+rho[19]*uy_sq[19]+rho[18]*uy_sq[18]+rho[17]*uy_sq[17]+rho[16]*uy_sq[16]+rho[15]*uy_sq[15]+rho[14]*uy_sq[14]+rho[13]*uy_sq[13]+rho[12]*uy_sq[12]+rho[11]*uy_sq[11]+rho[10]*uy_sq[10]+rho[9]*uy_sq[9]+rho[8]*uy_sq[8]+rho[7]*uy_sq[7]+rho[6]*uy_sq[6]+rho[5]*uy_sq[5]+rho[4]*uy_sq[4]+rho[3]*uy_sq[3]+rho[2]*uy_sq[2]+rho[1]*uy_sq[1]+rho[0]*uy_sq[0]; 
  int_pkpm_vars[6] += rho[26]*uz_sq[26]+rho[25]*uz_sq[25]+rho[24]*uz_sq[24]+rho[23]*uz_sq[23]+rho[22]*uz_sq[22]+rho[21]*uz_sq[21]+rho[20]*uz_sq[20]+rho[19]*uz_sq[19]+rho[18]*uz_sq[18]+rho[17]*uz_sq[17]+rho[16]*uz_sq[16]+rho[15]*uz_sq[15]+rho[14]*uz_sq[14]+rho[13]*uz_sq[13]+rho[12]*uz_sq[12]+rho[11]*uz_sq[11]+rho[10]*uz_sq[10]+rho[9]*uz_sq[9]+rho[8]*uz_sq[8]+rho[7]*uz_sq[7]+rho[6]*uz_sq[6]+rho[5]*uz_sq[5]+rho[4]*uz_sq[4]+rho[3]*uz_sq[3]+rho[2]*uz_sq[2]+rho[1]*uz_sq[1]+rho[0]*uz_sq[0]; 
  int_pkpm_vars[7] += 2.828427124746191*p_parallel[0]; 
  int_pkpm_vars[8] += 2.828427124746191*p_perp[0]; 
} 
