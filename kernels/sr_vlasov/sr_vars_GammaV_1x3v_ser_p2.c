#include <gkyl_sr_Gamma_kernels.h> 
#include <gkyl_basis_ser_1x_p2_exp_sq.h> 
#include <gkyl_basis_ser_1x_p2_sqrt.h> 
GKYL_CU_DH void sr_vars_GammaV_1x3v_ser_p2(const double *u_i, double* GKYL_RESTRICT u_i_sq, double* GKYL_RESTRICT GammaV, double* GKYL_RESTRICT GammaV_sq) 
{ 
  // u_i:       Input spatial components of bulk four-velocity = GammaV*V_drift. 
  // u_i_sq:    Output squared spatial components of bulk four-velocity = u_i^2. 
  // GammaV:    Output bulk four-velocity Lorentz factor = sqrt(1 + |u_i|^2). 
  // GammaV_sq: Output squared bulk four-velocity Lorentz factor = 1 + |u_i|^2. 
 
  const double *V_0 = &u_i[0]; 
  double *V_0_sq = &u_i_sq[0]; 
  ser_1x_p2_exp_sq(V_0, V_0_sq); 
 
  const double *V_1 = &u_i[3]; 
  double *V_1_sq = &u_i_sq[3]; 
  ser_1x_p2_exp_sq(V_1, V_1_sq); 
 
  const double *V_2 = &u_i[6]; 
  double *V_2_sq = &u_i_sq[6]; 
  ser_1x_p2_exp_sq(V_2, V_2_sq); 
 
  GammaV_sq[0] = V_2_sq[0]+V_1_sq[0]+V_0_sq[0]+1.414213562373095; 
  GammaV_sq[1] = V_2_sq[1]+V_1_sq[1]+V_0_sq[1]; 
  GammaV_sq[2] = V_2_sq[2]+V_1_sq[2]+V_0_sq[2]; 

  ser_1x_p2_sqrt(GammaV_sq, GammaV); 
} 
 
