#include <gkyl_sr_Gamma_kernels.h> 
#include <gkyl_basis_ser_3x_p1_exp_sq.h> 
#include <gkyl_basis_ser_3x_p1_inv.h> 
#include <gkyl_basis_ser_3x_p1_sqrt.h> 
GKYL_CU_DH void sr_Gamma_3x1v_ser_p1(const double *V, double* GKYL_RESTRICT Gamma) 
{ 
  // V:     Input velocity. 
  // Gamma: Gamma = 1/sqrt(1 - V^2/c^2). 
 
  const double *V_0 = &V[0]; 
  double V_0_sq[8] = {0.0}; 
  ser_3x_p1_exp_sq(V_0, V_0_sq); 
 
  double Gamma2_inv[8] = {0.0}; 
  double Gamma2[8] = {0.0}; 
 
  Gamma2_inv[0] = 2.828427124746191-1.0*V_0_sq[0]; 
  Gamma2_inv[1] = -1.0*V_0_sq[1]; 
  Gamma2_inv[2] = -1.0*V_0_sq[2]; 
  Gamma2_inv[3] = -1.0*V_0_sq[3]; 
  Gamma2_inv[4] = -1.0*V_0_sq[4]; 
  Gamma2_inv[5] = -1.0*V_0_sq[5]; 
  Gamma2_inv[6] = -1.0*V_0_sq[6]; 
  Gamma2_inv[7] = -1.0*V_0_sq[7]; 

  bool notCellAvg = true;
  if (notCellAvg && ((-1.837117307087383*Gamma2_inv[7])+1.060660171779821*Gamma2_inv[6]+1.060660171779821*Gamma2_inv[5]+1.060660171779821*Gamma2_inv[4]-0.6123724356957944*Gamma2_inv[3]-0.6123724356957944*Gamma2_inv[2]-0.6123724356957944*Gamma2_inv[1]+0.3535533905932737*Gamma2_inv[0] < 0)) notCellAvg = false; 
  if (notCellAvg && (1.837117307087383*Gamma2_inv[7]+1.060660171779821*Gamma2_inv[6]-1.060660171779821*Gamma2_inv[5]-1.060660171779821*Gamma2_inv[4]-0.6123724356957944*Gamma2_inv[3]-0.6123724356957944*Gamma2_inv[2]+0.6123724356957944*Gamma2_inv[1]+0.3535533905932737*Gamma2_inv[0] < 0)) notCellAvg = false; 
  if (notCellAvg && (1.837117307087383*Gamma2_inv[7]-1.060660171779821*Gamma2_inv[6]+1.060660171779821*Gamma2_inv[5]-1.060660171779821*Gamma2_inv[4]-0.6123724356957944*Gamma2_inv[3]+0.6123724356957944*Gamma2_inv[2]-0.6123724356957944*Gamma2_inv[1]+0.3535533905932737*Gamma2_inv[0] < 0)) notCellAvg = false; 
  if (notCellAvg && ((-1.837117307087383*Gamma2_inv[7])-1.060660171779821*Gamma2_inv[6]-1.060660171779821*Gamma2_inv[5]+1.060660171779821*Gamma2_inv[4]-0.6123724356957944*Gamma2_inv[3]+0.6123724356957944*Gamma2_inv[2]+0.6123724356957944*Gamma2_inv[1]+0.3535533905932737*Gamma2_inv[0] < 0)) notCellAvg = false; 
  if (notCellAvg && (1.837117307087383*Gamma2_inv[7]-1.060660171779821*Gamma2_inv[6]-1.060660171779821*Gamma2_inv[5]+1.060660171779821*Gamma2_inv[4]+0.6123724356957944*Gamma2_inv[3]-0.6123724356957944*Gamma2_inv[2]-0.6123724356957944*Gamma2_inv[1]+0.3535533905932737*Gamma2_inv[0] < 0)) notCellAvg = false; 
  if (notCellAvg && ((-1.837117307087383*Gamma2_inv[7])-1.060660171779821*Gamma2_inv[6]+1.060660171779821*Gamma2_inv[5]-1.060660171779821*Gamma2_inv[4]+0.6123724356957944*Gamma2_inv[3]-0.6123724356957944*Gamma2_inv[2]+0.6123724356957944*Gamma2_inv[1]+0.3535533905932737*Gamma2_inv[0] < 0)) notCellAvg = false; 
  if (notCellAvg && ((-1.837117307087383*Gamma2_inv[7])+1.060660171779821*Gamma2_inv[6]-1.060660171779821*Gamma2_inv[5]-1.060660171779821*Gamma2_inv[4]+0.6123724356957944*Gamma2_inv[3]+0.6123724356957944*Gamma2_inv[2]-0.6123724356957944*Gamma2_inv[1]+0.3535533905932737*Gamma2_inv[0] < 0)) notCellAvg = false; 
  if (notCellAvg && (1.837117307087383*Gamma2_inv[7]+1.060660171779821*Gamma2_inv[6]+1.060660171779821*Gamma2_inv[5]+1.060660171779821*Gamma2_inv[4]+0.6123724356957944*Gamma2_inv[3]+0.6123724356957944*Gamma2_inv[2]+0.6123724356957944*Gamma2_inv[1]+0.3535533905932737*Gamma2_inv[0] < 0)) notCellAvg = false; 
 
  if (notCellAvg) { 
  ser_3x_p1_inv(Gamma2_inv, Gamma2); 
  } else { 
  Gamma2[0] = 8.0/Gamma2_inv[0]; 
  Gamma2[1] = 0.0; 
  Gamma2[2] = 0.0; 
  Gamma2[3] = 0.0; 
  Gamma2[4] = 0.0; 
  Gamma2[5] = 0.0; 
  Gamma2[6] = 0.0; 
  Gamma2[7] = 0.0; 
  } 
  ser_3x_p1_sqrt(Gamma2, Gamma); 
} 
 
