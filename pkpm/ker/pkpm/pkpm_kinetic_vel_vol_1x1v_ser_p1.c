#include <gkyl_pkpm_kinetic_kernels.h> 
GKYL_CU_DH double pkpm_kinetic_vel_vol_1x1v_ser_p1(const double *w, const double *dxv,
  const double *a_P, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates. 
  // dxv[NDIM]: Cell spacing. 
  // a_P:       Pressure force closing the first moment. 
  // f:         Input distribution functions [F_0, G]. 
  // out:       Incremented output [F_0, G]. 
  const double dv1 = 2.0/dxv[1]; 
  const double *F_0 = &f[0]; 
  const double *G_1 = &f[6]; 
  double *out_F_0 = &out[0]; 
  double *out_G_1 = &out[6]; 
 
  out_F_0[2] += (1.224744871391589*F_0[1]*a_P[1]+1.224744871391589*F_0[0]*a_P[0])*dv1; 
  out_F_0[3] += (1.224744871391589*F_0[0]*a_P[1]+1.224744871391589*a_P[0]*F_0[1])*dv1; 
  out_F_0[4] += (2.7386127875258306*a_P[1]*F_0[3]+2.7386127875258306*a_P[0]*F_0[2])*dv1; 
  out_F_0[5] += (2.7386127875258306*a_P[0]*F_0[3]+2.7386127875258306*a_P[1]*F_0[2])*dv1; 
  out_G_1[2] += (1.224744871391589*G_1[1]*a_P[1]+1.224744871391589*G_1[0]*a_P[0])*dv1; 
  out_G_1[3] += (1.224744871391589*G_1[0]*a_P[1]+1.224744871391589*a_P[0]*G_1[1])*dv1; 
  out_G_1[4] += (2.7386127875258306*a_P[1]*G_1[3]+2.7386127875258306*a_P[0]*G_1[2])*dv1; 
  out_G_1[5] += (2.7386127875258306*a_P[0]*G_1[3]+2.7386127875258306*a_P[1]*G_1[2])*dv1; 
 
  return 0.0; 
} 
