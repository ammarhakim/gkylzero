#include <gkyl_sr_Gamma_kernels.h> 
#include <gkyl_basis_ser_1x_p1_exp_sq.h> 
#include <gkyl_basis_ser_1x_p1_sqrt.h> 
GKYL_CU_DH void sr_vars_GammaV_1x2v_ser_p1(const double *u_i, double* GKYL_RESTRICT u_i_sq, double* GKYL_RESTRICT GammaV, double* GKYL_RESTRICT GammaV_sq) 
{ 
  // u_i:       Input spatial components of bulk four-velocity = GammaV*V_drift. 
  // u_i_sq:    Output squared spatial components of bulk four-velocity = u_i^2. 
  // GammaV:    Output bulk four-velocity Lorentz factor = sqrt(1 + |u_i|^2). 
  // GammaV_sq: Output squared bulk four-velocity Lorentz factor = 1 + |u_i|^2. 
 
  const double *V_0 = &u_i[0]; 
  double *V_0_sq = &u_i_sq[0]; 
  ser_1x_p1_exp_sq(V_0, V_0_sq); 
 
  const double *V_1 = &u_i[2]; 
  double *V_1_sq = &u_i_sq[2]; 
  ser_1x_p1_exp_sq(V_1, V_1_sq); 
 
  GammaV_sq[0] = V_1_sq[0]+V_0_sq[0]+1.414213562373095; 
  GammaV_sq[1] = V_1_sq[1]+V_0_sq[1]; 

  ser_1x_p1_sqrt(GammaV_sq, GammaV); 
} 
 
