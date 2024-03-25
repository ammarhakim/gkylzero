#include <gkyl_sr_Gamma_kernels.h> 
#include <gkyl_basis_ser_1x_p2_exp_sq.h> 
#include <gkyl_basis_ser_1x_p2_inv.h> 
GKYL_CU_DH void sr_Gamma2_1x3v_ser_p2(const double *V, double* GKYL_RESTRICT Gamma2) 
{ 
  // V:      Input velocity. 
  // Gamma2: Gamma^2 = 1/(1 - V^2/c^2). 
 
  const double *V_0 = &V[0]; 
  double V_0_sq[3] = {0.0}; 
  ser_1x_p2_exp_sq(V_0, V_0_sq); 
 
  const double *V_1 = &V[3]; 
  double V_1_sq[3] = {0.0}; 
  ser_1x_p2_exp_sq(V_1, V_1_sq); 
 
  const double *V_2 = &V[6]; 
  double V_2_sq[3] = {0.0}; 
  ser_1x_p2_exp_sq(V_2, V_2_sq); 
 
  double V_sq_avg = (V_0_sq[0]+V_1_sq[0]+V_2_sq[0])/(1.414213562373095); 
  double V_sq_avg_lo = (V_0[0]*V_0[0] + V_1[0]*V_1[0] + V_2[0]*V_2[0])/(2.0); 
 
  double Gamma2_inv[3] = {0.0}; 
  Gamma2_inv[0] = (-1.0*V_2_sq[0])-1.0*V_1_sq[0]-1.0*V_0_sq[0]+1.414213562373095; 
  Gamma2_inv[1] = (-1.0*V_2_sq[1])-1.0*V_1_sq[1]-1.0*V_0_sq[1]; 
  Gamma2_inv[2] = (-1.0*V_2_sq[2])-1.0*V_1_sq[2]-1.0*V_0_sq[2]; 

  // Check if cell average of Gamma^{-2} = 1 - V^2/c^2 < 0. 
  if (V_sq_avg >= 1.0) { 
    Gamma2_inv[0] = 1.414213562373095*(1.0 - V_sq_avg_lo); 
    Gamma2_inv[1] = 0.0; 
    Gamma2_inv[2] = 0.0; 
  } 
 
  int cell_avg = 0;
  // Check if Gamma^{-2} = 1 - V^2/c^2 < 0 at control points. 
  if (1.58113883008419*Gamma2_inv[2]-1.224744871391589*Gamma2_inv[1]+0.7071067811865475*Gamma2_inv[0] < 0.0) cell_avg = 1; 
  if (0.7071067811865475*Gamma2_inv[0]-0.7905694150420947*Gamma2_inv[2] < 0.0) cell_avg = 1; 
  if (1.58113883008419*Gamma2_inv[2]+1.224744871391589*Gamma2_inv[1]+0.7071067811865475*Gamma2_inv[0] < 0.0) cell_avg = 1; 
 
  if (cell_avg) { 
    Gamma2[0] = 2.0/Gamma2_inv[0]; 
    Gamma2[1] = 0.0; 
    Gamma2[2] = 0.0; 
  } else { 
    ser_1x_p2_inv(Gamma2_inv, Gamma2); 
  } 
} 
 
