#include <gkyl_pkpm_vars_kernels.h> 
GKYL_CU_DH int pkpm_vars_theta_2x_ser_p1(const double *vlasov_pkpm_moms, double* GKYL_RESTRICT theta) 
{ 
  // vlasov_pkpm_moms: Input [rho, p_parallel, p_perp, M1], moments computed from the kinetic equation. 
  // theta:            Output T_perp/m = p_perp/rho, evaluated at the 4 Gauss-Legendre nodes 
  //                   and interpolated to the p=1 basis. 
  // Returns the number of nodes where rho <= 0 (theta is set to zero at such nodes). 
 
  const double *rho = &vlasov_pkpm_moms[0]; 
  const double *p_perp = &vlasov_pkpm_moms[8]; 
 
  double theta_n[4] = {0.0}; 
  int num_null_rho = 0; 
 
  { 
    const double rho_q = 0.5*rho[3]-0.5*rho[2]-0.5*rho[1]+0.5*rho[0]; 
    const double p_perp_q = 0.5*p_perp[3]-0.5*p_perp[2]-0.5*p_perp[1]+0.5*p_perp[0]; 
    if (rho_q > 0.0) { 
      theta_n[0] = p_perp_q/rho_q; 
    } else { 
      num_null_rho += 1; 
    } 
  } 
  { 
    const double rho_q = -(0.5*rho[3])+0.5*rho[2]-0.5*rho[1]+0.5*rho[0]; 
    const double p_perp_q = -(0.5*p_perp[3])+0.5*p_perp[2]-0.5*p_perp[1]+0.5*p_perp[0]; 
    if (rho_q > 0.0) { 
      theta_n[1] = p_perp_q/rho_q; 
    } else { 
      num_null_rho += 1; 
    } 
  } 
  { 
    const double rho_q = -(0.5*rho[3])-0.5*rho[2]+0.5*rho[1]+0.5*rho[0]; 
    const double p_perp_q = -(0.5*p_perp[3])-0.5*p_perp[2]+0.5*p_perp[1]+0.5*p_perp[0]; 
    if (rho_q > 0.0) { 
      theta_n[2] = p_perp_q/rho_q; 
    } else { 
      num_null_rho += 1; 
    } 
  } 
  { 
    const double rho_q = 0.5*rho[3]+0.5*rho[2]+0.5*rho[1]+0.5*rho[0]; 
    const double p_perp_q = 0.5*p_perp[3]+0.5*p_perp[2]+0.5*p_perp[1]+0.5*p_perp[0]; 
    if (rho_q > 0.0) { 
      theta_n[3] = p_perp_q/rho_q; 
    } else { 
      num_null_rho += 1; 
    } 
  } 
 
  theta[0] = 0.5*theta_n[3]+0.5*theta_n[2]+0.5*theta_n[1]+0.5*theta_n[0]; 
  theta[1] = 0.5*theta_n[3]+0.5*theta_n[2]-0.5*theta_n[1]-0.5*theta_n[0]; 
  theta[2] = 0.5*theta_n[3]-0.5*theta_n[2]+0.5*theta_n[1]-0.5*theta_n[0]; 
  theta[3] = 0.5*theta_n[3]-0.5*theta_n[2]-0.5*theta_n[1]+0.5*theta_n[0]; 
 
  return num_null_rho; 
} 
