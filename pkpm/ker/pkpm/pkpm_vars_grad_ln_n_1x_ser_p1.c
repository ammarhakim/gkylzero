#include <gkyl_pkpm_vars_kernels.h> 
GKYL_CU_DH int pkpm_vars_grad_ln_n_1x_ser_p1(const double *rho, const double *grad_rho, double* GKYL_RESTRICT grad_ln_n) 
{ 
  // rho:       Input mass density. 
  // grad_rho:  Input weak gradient of [rho, 0, 0] (pkpm_vars_grad_u_* layout, cdim x 3 expansions). 
  // grad_ln_n: Output [d_x ln n, ...] (cdim expansions): grad(rho) divided by the cell average of rho. 
  // Returns 1 if the cell average of rho is not positive (the quotient is zero then), 0 otherwise. 
 
  const double *grho_x = &grad_rho[0]; 
  double *gln_x = &grad_ln_n[0]; 
 
  const double rho_avg = 0.7071067811865475*rho[0]; 
  const double inv_rho = rho_avg > 0.0 ? 1.0/rho_avg : 0.0; 
 
  gln_x[0] = inv_rho*grho_x[0]; 
  gln_x[1] = inv_rho*grho_x[1]; 
 
  return rho_avg > 0.0 ? 0 : 1; 
} 
