#include <gkyl_euler_kernels.h> 
#include <gkyl_basis_ser_1x_p2_inv.h> 
GKYL_CU_DH void euler_pkpm_p_force_1x_ser_p2(const double *bvar, const double *div_p, const double *vlasov_pkpm_moms, double* GKYL_RESTRICT p_force) 
{ 
  // bvar:             magnetic field unit vector (nine components; first three components, b_i, other six components, b_i b_j.) 
  // div_p:            Volume expansion of div(p).
  // vlasov_pkpm_moms: [rho, p_parallel, q_parallel], Moments computed from kinetic equation in pkpm model.
  // p_force:          total pressure force = 1/rho (b . div(P) + p_perp div(b)) for Euler PKPM.

  const double *bx = &bvar[0]; 
  const double *by = &bvar[3]; 
  const double *bz = &bvar[6]; 
  const double *div_p_x = &div_p[0]; 
  const double *div_p_y = &div_p[3]; 
  const double *div_p_z = &div_p[6]; 
  const double *rho = &vlasov_pkpm_moms[0]; 

  double b_div_p[3] = {0.0}; 
  double rho_inv[3] = {0.0}; 
  ser_1x_p2_inv(rho, rho_inv); 
  b_div_p[0] = 0.7071067811865475*bz[2]*div_p_z[2]+0.7071067811865475*by[2]*div_p_y[2]+0.7071067811865475*bx[2]*div_p_x[2]+0.7071067811865475*bz[1]*div_p_z[1]+0.7071067811865475*by[1]*div_p_y[1]+0.7071067811865475*bx[1]*div_p_x[1]+0.7071067811865475*bz[0]*div_p_z[0]+0.7071067811865475*by[0]*div_p_y[0]+0.7071067811865475*bx[0]*div_p_x[0]; 
  b_div_p[1] = 0.6324555320336759*bz[1]*div_p_z[2]+0.6324555320336759*by[1]*div_p_y[2]+0.6324555320336759*bx[1]*div_p_x[2]+0.6324555320336759*div_p_z[1]*bz[2]+0.6324555320336759*div_p_y[1]*by[2]+0.6324555320336759*div_p_x[1]*bx[2]+0.7071067811865475*bz[0]*div_p_z[1]+0.7071067811865475*by[0]*div_p_y[1]+0.7071067811865475*bx[0]*div_p_x[1]+0.7071067811865475*div_p_z[0]*bz[1]+0.7071067811865475*div_p_y[0]*by[1]+0.7071067811865475*div_p_x[0]*bx[1]; 
  b_div_p[2] = 0.4517539514526256*bz[2]*div_p_z[2]+0.7071067811865475*bz[0]*div_p_z[2]+0.4517539514526256*by[2]*div_p_y[2]+0.7071067811865475*by[0]*div_p_y[2]+0.4517539514526256*bx[2]*div_p_x[2]+0.7071067811865475*bx[0]*div_p_x[2]+0.7071067811865475*div_p_z[0]*bz[2]+0.7071067811865475*div_p_y[0]*by[2]+0.7071067811865475*div_p_x[0]*bx[2]+0.6324555320336759*bz[1]*div_p_z[1]+0.6324555320336759*by[1]*div_p_y[1]+0.6324555320336759*bx[1]*div_p_x[1]; 

  p_force[0] = 0.7071067811865475*b_div_p[2]*rho_inv[2]+0.7071067811865475*b_div_p[1]*rho_inv[1]+0.7071067811865475*b_div_p[0]*rho_inv[0]; 
  p_force[1] = 0.6324555320336759*b_div_p[1]*rho_inv[2]+0.6324555320336759*rho_inv[1]*b_div_p[2]+0.7071067811865475*b_div_p[0]*rho_inv[1]+0.7071067811865475*rho_inv[0]*b_div_p[1]; 
  p_force[2] = 0.4517539514526256*b_div_p[2]*rho_inv[2]+0.7071067811865475*b_div_p[0]*rho_inv[2]+0.7071067811865475*rho_inv[0]*b_div_p[2]+0.6324555320336759*b_div_p[1]*rho_inv[1]; 

} 
