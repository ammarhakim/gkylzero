#include <gkyl_euler_pkpm_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
GKYL_CU_DH void pkpm_vars_p_force_3x_tensor_p2(const double *prim_c, const double *div_b, 
  double* GKYL_RESTRICT pkpm_accel) 
{ 
  // prim_c:     Input volume expansion of primitive variables in center cell. 
  //             [1/rho div(p_par b), T_perp/m, m/T_perp, 3*Txx/m, 3*Tyy/m, 3*Tzz/m]. 
  // div_b:      Input volume expansion of div(b) in center cell. 
  // pkpm_accel: Output volume expansion of pkpm acceleration variables. 

  const double *pkpm_div_ppar = &prim_c[0]; 
  const double *T_perp_over_m = &prim_c[27]; 

  double *p_perp_div_b = &pkpm_accel[0]; 
  double *p_force = &pkpm_accel[54]; 

  binop_mul_3d_tensor_p2(T_perp_over_m, div_b, p_perp_div_b); 

  p_force[0] += pkpm_div_ppar[0]-1.0*p_perp_div_b[0]; 
  p_force[1] += pkpm_div_ppar[1]-1.0*p_perp_div_b[1]; 
  p_force[2] += pkpm_div_ppar[2]-1.0*p_perp_div_b[2]; 
  p_force[3] += pkpm_div_ppar[3]-1.0*p_perp_div_b[3]; 
  p_force[4] += pkpm_div_ppar[4]-1.0*p_perp_div_b[4]; 
  p_force[5] += pkpm_div_ppar[5]-1.0*p_perp_div_b[5]; 
  p_force[6] += pkpm_div_ppar[6]-1.0*p_perp_div_b[6]; 
  p_force[7] += pkpm_div_ppar[7]-1.0*p_perp_div_b[7]; 
  p_force[8] += pkpm_div_ppar[8]-1.0*p_perp_div_b[8]; 
  p_force[9] += pkpm_div_ppar[9]-1.0*p_perp_div_b[9]; 
  p_force[10] += pkpm_div_ppar[10]-1.0*p_perp_div_b[10]; 
  p_force[11] += pkpm_div_ppar[11]-1.0*p_perp_div_b[11]; 
  p_force[12] += pkpm_div_ppar[12]-1.0*p_perp_div_b[12]; 
  p_force[13] += pkpm_div_ppar[13]-1.0*p_perp_div_b[13]; 
  p_force[14] += pkpm_div_ppar[14]-1.0*p_perp_div_b[14]; 
  p_force[15] += pkpm_div_ppar[15]-1.0*p_perp_div_b[15]; 
  p_force[16] += pkpm_div_ppar[16]-1.0*p_perp_div_b[16]; 
  p_force[17] += pkpm_div_ppar[17]-1.0*p_perp_div_b[17]; 
  p_force[18] += pkpm_div_ppar[18]-1.0*p_perp_div_b[18]; 
  p_force[19] += pkpm_div_ppar[19]-1.0*p_perp_div_b[19]; 
  p_force[20] += pkpm_div_ppar[20]-1.0*p_perp_div_b[20]; 
  p_force[21] += pkpm_div_ppar[21]-1.0*p_perp_div_b[21]; 
  p_force[22] += pkpm_div_ppar[22]-1.0*p_perp_div_b[22]; 
  p_force[23] += pkpm_div_ppar[23]-1.0*p_perp_div_b[23]; 
  p_force[24] += pkpm_div_ppar[24]-1.0*p_perp_div_b[24]; 
  p_force[25] += pkpm_div_ppar[25]-1.0*p_perp_div_b[25]; 
  p_force[26] += pkpm_div_ppar[26]-1.0*p_perp_div_b[26]; 

} 
