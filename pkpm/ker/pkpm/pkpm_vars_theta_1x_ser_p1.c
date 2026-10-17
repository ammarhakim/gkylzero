#include <gkyl_pkpm_vars_kernels.h> 
GKYL_CU_DH int pkpm_vars_theta_1x_ser_p1(const double *vlasov_pkpm_moms, double* GKYL_RESTRICT theta) 
{ 
  // vlasov_pkpm_moms: Input [rho, p_parallel, p_perp, M1], moments computed from the kinetic equation. 
  // theta:            Output T_perp/m = p_perp/rho, evaluated at the 2 Gauss-Legendre nodes 
  //                   and interpolated to the p=1 basis. 
  // Returns the number of nodes where rho <= 0 (theta is set to zero at such nodes). 
 
  const double *rho = &vlasov_pkpm_moms[0]; 
  const double *p_perp = &vlasov_pkpm_moms[4]; 
 
  double theta_n[2] = {0.0}; 
  int num_null_rho = 0; 
 
  { 
    const double rho_q = 0.7071067811865475*rho[0]-0.7071067811865475*rho[1]; 
    const double p_perp_q = 0.7071067811865475*p_perp[0]-0.7071067811865475*p_perp[1]; 
    if (rho_q > 0.0) { 
      theta_n[0] = p_perp_q/rho_q; 
    } else { 
      num_null_rho += 1; 
    } 
  } 
  { 
    const double rho_q = 0.7071067811865475*rho[1]+0.7071067811865475*rho[0]; 
    const double p_perp_q = 0.7071067811865475*p_perp[1]+0.7071067811865475*p_perp[0]; 
    if (rho_q > 0.0) { 
      theta_n[1] = p_perp_q/rho_q; 
    } else { 
      num_null_rho += 1; 
    } 
  } 
 
  theta[0] = 0.7071067811865475*theta_n[1]+0.7071067811865475*theta_n[0]; 
  theta[1] = 0.7071067811865475*theta_n[1]-0.7071067811865475*theta_n[0]; 
 
  return num_null_rho; 
} 
