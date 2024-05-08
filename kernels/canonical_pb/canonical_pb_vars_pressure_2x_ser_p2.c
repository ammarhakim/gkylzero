#include <gkyl_canonical_pb_kernels.h>  
#include <gkyl_binop_mul_ser.h> 
GKYL_CU_DH void canonical_pb_vars_pressure_2x_ser_p2(const double *h_ij_inv, const double *M2_ij, const double *v_j, const double *nv_i, double* GKYL_RESTRICT d_Jv_P) 
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
  const double *NVy = &nv_i[8]; 

  const double *Vx = &v_j[0]; 
  const double *Vy = &v_j[8]; 

  const double *M2xx = &M2_ij[0]; 
  const double *M2xy = &M2_ij[8]; 
  const double *M2yy = &M2_ij[16]; 

  const double *Hxx = &h_ij_inv[0]; 
  const double *Hxy = &h_ij_inv[8]; 
  const double *Hyy = &h_ij_inv[16]; 

  // h^{ij}M2_ij 
  double Hxx_M2xx[8] = {0.0}; 
  binop_mul_2d_ser_p2(Hxx, M2xx, Hxx_M2xx); 
 
  double Hxy_M2xy[8] = {0.0}; 
  binop_mul_2d_ser_p2(Hxy, M2xy, Hxy_M2xy); 
 
  double Hyy_M2yy[8] = {0.0}; 
  binop_mul_2d_ser_p2(Hyy, M2yy, Hyy_M2yy); 
 
  // h^{ij}*nv_i*v_j 
  double Hxx_M1x[8] = {0.0}; 
  double Hxx_M1x_Vx[8] = {0.0}; 
  binop_mul_2d_ser_p2(Hxx, NVx, Hxx_M1x); 
  binop_mul_2d_ser_p2(Hxx_M1x, Vx, Hxx_M1x_Vx); 
 
  double Hxy_M1x[8] = {0.0}; 
  double Hxy_M1x_Vy[8] = {0.0}; 
  binop_mul_2d_ser_p2(Hxy, NVx, Hxy_M1x); 
  binop_mul_2d_ser_p2(Hxy_M1x, Vy, Hxy_M1x_Vy); 
 
  double Hyy_M1y[8] = {0.0}; 
  double Hyy_M1y_Vy[8] = {0.0}; 
  binop_mul_2d_ser_p2(Hyy, NVy, Hyy_M1y); 
  binop_mul_2d_ser_p2(Hyy_M1y, Vy, Hyy_M1y_Vy); 
 
  d_Jv_P[0] = 0.0; 
  d_Jv_P[0] +=  Hxx_M2xx[0] - Hxx_M1x_Vx[0]; 
  d_Jv_P[0] += (Hxy_M2xy[0] - Hxy_M1x_Vy[0])*2.0; 
  d_Jv_P[0] +=  Hyy_M2yy[0] - Hyy_M1y_Vy[0]; 
  d_Jv_P[1] = 0.0; 
  d_Jv_P[1] +=  Hxx_M2xx[1] - Hxx_M1x_Vx[1]; 
  d_Jv_P[1] += (Hxy_M2xy[1] - Hxy_M1x_Vy[1])*2.0; 
  d_Jv_P[1] +=  Hyy_M2yy[1] - Hyy_M1y_Vy[1]; 
  d_Jv_P[2] = 0.0; 
  d_Jv_P[2] +=  Hxx_M2xx[2] - Hxx_M1x_Vx[2]; 
  d_Jv_P[2] += (Hxy_M2xy[2] - Hxy_M1x_Vy[2])*2.0; 
  d_Jv_P[2] +=  Hyy_M2yy[2] - Hyy_M1y_Vy[2]; 
  d_Jv_P[3] = 0.0; 
  d_Jv_P[3] +=  Hxx_M2xx[3] - Hxx_M1x_Vx[3]; 
  d_Jv_P[3] += (Hxy_M2xy[3] - Hxy_M1x_Vy[3])*2.0; 
  d_Jv_P[3] +=  Hyy_M2yy[3] - Hyy_M1y_Vy[3]; 
  d_Jv_P[4] = 0.0; 
  d_Jv_P[4] +=  Hxx_M2xx[4] - Hxx_M1x_Vx[4]; 
  d_Jv_P[4] += (Hxy_M2xy[4] - Hxy_M1x_Vy[4])*2.0; 
  d_Jv_P[4] +=  Hyy_M2yy[4] - Hyy_M1y_Vy[4]; 
  d_Jv_P[5] = 0.0; 
  d_Jv_P[5] +=  Hxx_M2xx[5] - Hxx_M1x_Vx[5]; 
  d_Jv_P[5] += (Hxy_M2xy[5] - Hxy_M1x_Vy[5])*2.0; 
  d_Jv_P[5] +=  Hyy_M2yy[5] - Hyy_M1y_Vy[5]; 
  d_Jv_P[6] = 0.0; 
  d_Jv_P[6] +=  Hxx_M2xx[6] - Hxx_M1x_Vx[6]; 
  d_Jv_P[6] += (Hxy_M2xy[6] - Hxy_M1x_Vy[6])*2.0; 
  d_Jv_P[6] +=  Hyy_M2yy[6] - Hyy_M1y_Vy[6]; 
  d_Jv_P[7] = 0.0; 
  d_Jv_P[7] +=  Hxx_M2xx[7] - Hxx_M1x_Vx[7]; 
  d_Jv_P[7] += (Hxy_M2xy[7] - Hxy_M1x_Vy[7])*2.0; 
  d_Jv_P[7] +=  Hyy_M2yy[7] - Hyy_M1y_Vy[7]; 
 
} 
