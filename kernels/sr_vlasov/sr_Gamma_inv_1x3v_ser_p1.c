#include <gkyl_sr_Gamma_kernels.h> 
#include <gkyl_basis_ser_1x_p1_exp_sq.h> 
#include <gkyl_basis_ser_1x_p1_sqrt.h> 
GKYL_CU_DH void sr_Gamma_inv_1x3v_ser_p1(const double *V, double* GKYL_RESTRICT Gamma_inv) 
{ 
  // V:     Input velocity. 
  // Gamma: Gamma = 1/sqrt(1 - V^2/c^2). 
 
  const double *V_0 = &V[0]; 
  double V_0_sq[2] = {0.0}; 
  ser_1x_p1_exp_sq(V_0, V_0_sq); 
 
  const double *V_1 = &V[2]; 
  double V_1_sq[2] = {0.0}; 
  ser_1x_p1_exp_sq(V_1, V_1_sq); 
 
  const double *V_2 = &V[4]; 
  double V_2_sq[2] = {0.0}; 
  ser_1x_p1_exp_sq(V_2, V_2_sq); 
 
  double V_sq_avg = (V_0_sq[0]+V_1_sq[0]+V_2_sq[0])/(1.414213562373095); 
  double V_sq_avg_lo = (V_0[0]*V_0[0] + V_1[0]*V_1[0] + V_2[0]*V_2[0])/(2.0); 
 
  double Gamma2_inv[2] = {0.0}; 
  Gamma2_inv[0] = (-1.0*V_2_sq[0])-1.0*V_1_sq[0]-1.0*V_0_sq[0]+1.414213562373095; 
  Gamma2_inv[1] = (-1.0*V_2_sq[1])-1.0*V_1_sq[1]-1.0*V_0_sq[1]; 

  // Check if cell average of Gamma^{-2} = 1 - V^2/c^2 < 0. 
  if (V_sq_avg >= 1.0) { 
    Gamma2_inv[0] = 1.414213562373095*(1.0 - V_sq_avg_lo); 
    Gamma2_inv[1] = 0.0; 
  } 
 
  int cell_avg = 0;
  // Check if Gamma^{-2} = 1 - V^2/c^2 < 0 at control points. 
  if (0.7071067811865475*Gamma2_inv[0]-1.224744871391589*Gamma2_inv[1] < 0.0) cell_avg = 1; 
  if (1.224744871391589*Gamma2_inv[1]+0.7071067811865475*Gamma2_inv[0] < 0.0) cell_avg = 1; 
 
  if (cell_avg) { 
    Gamma2_inv[1] = 0.0; 
    ser_1x_p1_sqrt(Gamma2_inv, Gamma_inv); 
  } else { 
    ser_1x_p1_sqrt(Gamma2_inv, Gamma_inv); 
  } 
} 
 
