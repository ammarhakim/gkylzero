#include <gkyl_canonical_pb_kernels.h>  
#include <gkyl_binop_mul_ser.h> 
GKYL_CU_DH void canonical_pb_vars_pressure_1x_ser_p1(const double *h_ij_inv, const double *M2_ij, const double *v_j, const double *nv_i, double* GKYL_RESTRICT d_Jv_P) 
{ 
  // h_ij_inv:         Input volume expansion of the inverse metric tensor.
  //                   [Hxx, Hxy, Hxz, 
  //                    Hxy, Hyy, Hyz, 
  //                    Hxz, Hyz, Hzz] 
  // M2_ij:             Input volume expansion of the M2_ij moment.
  //                   [M2xx, M2xy, M2xz, 
  //                    M2xy, M2yy, M2yz, 
  //                    M2xz, M2yz, M2zz] 
  // v_j:              Input volume expansion of V_drift.
  //                   [vx, vy, vz] 
  // nv_i:              Input volume expansion of M1i = N*Vdrift.
  //                   [nvx, nvy, nvz] 
  // d_Jv_P:            Output volume expansion of d*P*Jv = h^{ij}*M2_{ij} - n*h^{ij}*u_i*u_j .

  const double *NVx = &nv_i[0]; 

  const double *Vx = &v_j[0]; 

  const double *M2xx = &M2_ij[0]; 

  const double *Hxx = &h_ij_inv[0]; 

  // h^{ij}M2_ij 
  double Hxx_M2xx[2] = {0.0}; 
  binop_mul_1d_ser_p1(Hxx, M2xx, Hxx_M2xx); 
 
  // h^{ij}*nv_i*v_j 
  double Hxx_M1x[2] = {0.0}; 
  double Hxx_M1x_Vx[2] = {0.0}; 
  binop_mul_1d_ser_p1(Hxx, NVx, Hxx_M1x); 
  binop_mul_1d_ser_p1(Hxx_M1x, Vx, Hxx_M1x_Vx); 
 
  d_Jv_P[0] = 0.0; 
  d_Jv_P[0] +=  Hxx_M2xx[0] - Hxx_M1x_Vx[0]; 
  d_Jv_P[1] = 0.0; 
  d_Jv_P[1] +=  Hxx_M2xx[1] - Hxx_M1x_Vx[1]; 
 
} 
