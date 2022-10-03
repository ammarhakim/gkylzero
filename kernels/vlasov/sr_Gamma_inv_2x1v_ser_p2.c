#include <gkyl_sr_Gamma_kernels.h> 
#include <gkyl_basis_ser_2x_p2_exp_sq.h> 
#include <gkyl_basis_ser_2x_p2_sqrt.h> 
GKYL_CU_DH void sr_Gamma_inv_2x1v_ser_p2(const double *V, double* GKYL_RESTRICT Gamma_inv) 
{ 
  // V:     Input velocity. 
  // Gamma: Gamma = 1/sqrt(1 - V^2/c^2). 
 
  const double *V_0 = &V[0]; 
  double V_0_sq[8] = {0.0}; 
  ser_2x_p2_exp_sq(V_0, V_0_sq); 
 
  double Gamma2_inv[8] = {0.0}; 
 
  Gamma2_inv[0] = 2.0-1.0*V_0_sq[0]; 
  Gamma2_inv[1] = -1.0*V_0_sq[1]; 
  Gamma2_inv[2] = -1.0*V_0_sq[2]; 
  Gamma2_inv[3] = -1.0*V_0_sq[3]; 
  Gamma2_inv[4] = -1.0*V_0_sq[4]; 
  Gamma2_inv[5] = -1.0*V_0_sq[5]; 
  Gamma2_inv[6] = -1.0*V_0_sq[6]; 
  Gamma2_inv[7] = -1.0*V_0_sq[7]; 

  ser_2x_p2_sqrt(Gamma2_inv, Gamma_inv); 
} 
 
