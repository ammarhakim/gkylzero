#include <gkyl_pkpm_vars_kernels.h> 
GKYL_CU_DH int pkpm_vars_p_force_closure_2x_ser_p1(const double *vlasov_pkpm_moms, const double *rhs_moms, double* GKYL_RESTRICT a_P) 
{ 
  // vlasov_pkpm_moms: Input [rho, p_parallel, p_perp, M1], moments of the kinetic equation (stage state). 
  // rhs_moms:         Input [rho, p_parallel, p_perp, M1] moments of the stage-A right-hand side of F_0. 
  // a_P:              Output pressure force closing the first moment: rho a_P = -M1(rhs_A) weakly, 
  //                   evaluated at the 4 Gauss-Legendre nodes and interpolated to the p=1 basis. 
  // Returns the number of nodes where rho <= 0 (a_P is set to zero in the cell in that case). 
 
  const double *rho = &vlasov_pkpm_moms[0]; 
  const double *m1_rhs = &rhs_moms[12]; 
 
  double a_P_n[4] = {0.0}; 
  int num_null_rho = 0; 
 
  { 
    const double rho_q = 0.5*rho[3]-0.5*rho[2]-0.5*rho[1]+0.5*rho[0]; 
    const double m1_q = 0.5*m1_rhs[3]-0.5*m1_rhs[2]-0.5*m1_rhs[1]+0.5*m1_rhs[0]; 
    if (rho_q > 0.0) { 
      a_P_n[0] = -m1_q/rho_q; 
    } else { 
      num_null_rho += 1; 
    } 
  } 
  { 
    const double rho_q = -(0.5*rho[3])+0.5*rho[2]-0.5*rho[1]+0.5*rho[0]; 
    const double m1_q = -(0.5*m1_rhs[3])+0.5*m1_rhs[2]-0.5*m1_rhs[1]+0.5*m1_rhs[0]; 
    if (rho_q > 0.0) { 
      a_P_n[1] = -m1_q/rho_q; 
    } else { 
      num_null_rho += 1; 
    } 
  } 
  { 
    const double rho_q = -(0.5*rho[3])-0.5*rho[2]+0.5*rho[1]+0.5*rho[0]; 
    const double m1_q = -(0.5*m1_rhs[3])-0.5*m1_rhs[2]+0.5*m1_rhs[1]+0.5*m1_rhs[0]; 
    if (rho_q > 0.0) { 
      a_P_n[2] = -m1_q/rho_q; 
    } else { 
      num_null_rho += 1; 
    } 
  } 
  { 
    const double rho_q = 0.5*rho[3]+0.5*rho[2]+0.5*rho[1]+0.5*rho[0]; 
    const double m1_q = 0.5*m1_rhs[3]+0.5*m1_rhs[2]+0.5*m1_rhs[1]+0.5*m1_rhs[0]; 
    if (rho_q > 0.0) { 
      a_P_n[3] = -m1_q/rho_q; 
    } else { 
      num_null_rho += 1; 
    } 
  } 
 
  if (num_null_rho > 0) { 
    a_P[0] = 0.0; 
    a_P[1] = 0.0; 
    a_P[2] = 0.0; 
    a_P[3] = 0.0; 
    return num_null_rho; 
  } 
 
  a_P[0] = 0.5*a_P_n[3]+0.5*a_P_n[2]+0.5*a_P_n[1]+0.5*a_P_n[0]; 
  a_P[1] = 0.5*a_P_n[3]+0.5*a_P_n[2]-0.5*a_P_n[1]-0.5*a_P_n[0]; 
  a_P[2] = 0.5*a_P_n[3]-0.5*a_P_n[2]+0.5*a_P_n[1]-0.5*a_P_n[0]; 
  a_P[3] = 0.5*a_P_n[3]-0.5*a_P_n[2]-0.5*a_P_n[1]+0.5*a_P_n[0]; 
 
  return 0; 
} 
