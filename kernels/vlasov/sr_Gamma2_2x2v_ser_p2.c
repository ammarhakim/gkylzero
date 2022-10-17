#include <gkyl_sr_Gamma_kernels.h> 
#include <gkyl_basis_ser_2x_p2_exp_sq.h> 
#include <gkyl_basis_ser_2x_p2_inv.h> 
GKYL_CU_DH void sr_Gamma2_2x2v_ser_p2(const double *V, double* GKYL_RESTRICT Gamma2) 
{ 
  // V:      Input velocity. 
  // Gamma2: Gamma^2 = 1/(1 - V^2/c^2). 
 
  const double *V_0 = &V[0]; 
  double V_0_sq[8] = {0.0}; 
  ser_2x_p2_exp_sq(V_0, V_0_sq); 
 
  const double *V_1 = &V[8]; 
  double V_1_sq[8] = {0.0}; 
  ser_2x_p2_exp_sq(V_1, V_1_sq); 
 
  double Gamma2_inv[8] = {0.0}; 
 
  Gamma2_inv[0] = (-1.0*V_1_sq[0])-1.0*V_0_sq[0]+2.0; 
  Gamma2_inv[1] = (-1.0*V_1_sq[1])-1.0*V_0_sq[1]; 
  Gamma2_inv[2] = (-1.0*V_1_sq[2])-1.0*V_0_sq[2]; 
  Gamma2_inv[3] = (-1.0*V_1_sq[3])-1.0*V_0_sq[3]; 
  Gamma2_inv[4] = (-1.0*V_1_sq[4])-1.0*V_0_sq[4]; 
  Gamma2_inv[5] = (-1.0*V_1_sq[5])-1.0*V_0_sq[5]; 
  Gamma2_inv[6] = (-1.0*V_1_sq[6])-1.0*V_0_sq[6]; 
  Gamma2_inv[7] = (-1.0*V_1_sq[7])-1.0*V_0_sq[7]; 

  bool notCellAvg = true;
  if (notCellAvg && ((-1.936491673103709*Gamma2_inv[7])-1.936491673103709*Gamma2_inv[6]+1.118033988749895*Gamma2_inv[5]+1.118033988749895*Gamma2_inv[4]+1.5*Gamma2_inv[3]-0.8660254037844386*Gamma2_inv[2]-0.8660254037844386*Gamma2_inv[1]+0.5*Gamma2_inv[0] < 0)) notCellAvg = false; 
  if (notCellAvg && (1.936491673103709*Gamma2_inv[7]-1.936491673103709*Gamma2_inv[6]+1.118033988749895*Gamma2_inv[5]+1.118033988749895*Gamma2_inv[4]-1.5*Gamma2_inv[3]-0.8660254037844386*Gamma2_inv[2]+0.8660254037844386*Gamma2_inv[1]+0.5*Gamma2_inv[0] < 0)) notCellAvg = false; 
  if (notCellAvg && ((-1.936491673103709*Gamma2_inv[7])+1.936491673103709*Gamma2_inv[6]+1.118033988749895*Gamma2_inv[5]+1.118033988749895*Gamma2_inv[4]-1.5*Gamma2_inv[3]+0.8660254037844386*Gamma2_inv[2]-0.8660254037844386*Gamma2_inv[1]+0.5*Gamma2_inv[0] < 0)) notCellAvg = false; 
  if (notCellAvg && (1.936491673103709*Gamma2_inv[7]+1.936491673103709*Gamma2_inv[6]+1.118033988749895*Gamma2_inv[5]+1.118033988749895*Gamma2_inv[4]+1.5*Gamma2_inv[3]+0.8660254037844386*Gamma2_inv[2]+0.8660254037844386*Gamma2_inv[1]+0.5*Gamma2_inv[0] < 0)) notCellAvg = false; 
 
  if (notCellAvg) { 
  ser_2x_p2_inv(Gamma2_inv, Gamma2); 
  } else { 
  Gamma2[0] = 4.0/Gamma2_inv[0]; 
  Gamma2[1] = 0.0; 
  Gamma2[2] = 0.0; 
  Gamma2[3] = 0.0; 
  Gamma2[4] = 0.0; 
  Gamma2[5] = 0.0; 
  Gamma2[6] = 0.0; 
  Gamma2[7] = 0.0; 
  } 
} 
 
