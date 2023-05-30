#include <gkyl_sr_Gamma_kernels.h> 
#include <gkyl_basis_ser_1x_p1_exp_sq.h> 
#include <gkyl_basis_ser_1x_p1_inv.h> 
#include <gkyl_basis_ser_1x_p1_sqrt.h> 
GKYL_CU_DH void sr_Gamma_1x2v_ser_p1(const double *V, double* GKYL_RESTRICT Gamma) 
{ 
  // V:     Input velocity. 
  // Gamma: Gamma = 1/sqrt(1 - V^2/c^2). 
 
  const double *V_0 = &V[0]; 
  double V_0_sq[2] = {0.0}; 
  ser_1x_p1_exp_sq(V_0, V_0_sq); 
 
  const double *V_1 = &V[2]; 
  double V_1_sq[2] = {0.0}; 
  ser_1x_p1_exp_sq(V_1, V_1_sq); 
 
  double Gamma2_inv[2] = {0.0}; 
  double Gamma2[2] = {0.0}; 
 
  Gamma2_inv[0] = (-1.0*V_1_sq[0])-1.0*V_0_sq[0]+1.414213562373095; 
  Gamma2_inv[1] = (-1.0*V_1_sq[1])-1.0*V_0_sq[1]; 

  bool notCellAvg = true;
  if (notCellAvg && (0.7071067811865475*Gamma2_inv[0]-1.224744871391589*Gamma2_inv[1] < 0)) notCellAvg = false; 
  if (notCellAvg && (1.224744871391589*Gamma2_inv[1]+0.7071067811865475*Gamma2_inv[0] < 0)) notCellAvg = false; 
 
  if (notCellAvg) { 
  ser_1x_p1_inv(Gamma2_inv, Gamma2); 
  } else { 
  Gamma2[0] = 2.0/Gamma2_inv[0]; 
  Gamma2[1] = 0.0; 
  } 
  ser_1x_p1_sqrt(Gamma2, Gamma); 
} 
 
