#include <gkyl_canonical_pb_kernels.h>  
#include <gkyl_binop_mul_ser.h> 
GKYL_CU_DH void canonical_pb_vars_m1i_contra_to_cov_1x1v_tensor_p2(const double *h_ij, const double *v_j, const double *nv_i, double* GKYL_RESTRICT v_j_cov, double* GKYL_RESTRICT nv_i_cov) 
{ 
  // h_ij:         Input volume expansion of the covariant components of the metric tensor.
  //                   [Hxx, Hxy, Hxz, 
  //                    Hxy, Hyy, Hyz, 
  //                    Hxz, Hyz, Hzz] 
  // v_j:              Input volume expansion of V_drift (contravariant components).
  //                   [vx, vy, vz] 
  // nv_i:              Input volume expansion of M1i = N*Vdrift (contravariant components).
  //                   [nvx, nvy, nvz] 
  // v_j_cov:           Input volume expansion of V_drift (covariant components).
  //                   [vx, vy, vz] 
  // nv_i_cov:          Input volume expansion of M1i = N*Vdrift (covariant components).
  //                   [nvx, nvy, nvz] 

  const double *NVx = &nv_i[0]; 

  const double *Vx = &v_j[0]; 

  double *NVx_cov = &nv_i_cov[0]; 

  double *Vx_cov = &v_j_cov[0]; 

  const double *Hxx = &h_ij[0]; 

  // h_{ij}v^j 
  double Hxx_Vx[3] = {0.0}; 
  binop_mul_1d_ser_p2(Hxx, Vx, Hxx_Vx); 
 
  // h_{ij}Jnv^j 
  double Hxx_NVx[3] = {0.0}; 
  binop_mul_1d_ser_p2(Hxx, NVx, Hxx_NVx); 
 
  // u_i_cov = h_{ij}v^j 
  Vx_cov[0] = Hxx_Vx[0]; 
  Vx_cov[1] = Hxx_Vx[1]; 
  Vx_cov[2] = Hxx_Vx[2]; 
 
  // Jnu_i_cov = h_{ij}Jnv^j 
  NVx_cov[0] = Hxx_NVx[0]; 
  NVx_cov[1] = Hxx_NVx[1]; 
  NVx_cov[2] = Hxx_NVx[2]; 
 
} 
