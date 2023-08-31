#include <gkyl_dg_diffusion_kernels.h> 
GKYL_CU_DH double dg_diffusion_pkpm_surfx_1x_ser_p2(const double* w, const double* dx, double D, 
  const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: Cell-center coordinates
  // dxv[NDIM]: Cell spacing
  // D: Diffusion coefficient in the center cell
  // ql: Input field in the left cell
  // qc: Input field in the center cell
  // qr: Input field in the right cell
  // out: Incremented output

  const double dx1 = 2.0/dx[0]; 
  const double J = pow(dx1, 2.0);

  const double *q0l = &ql[0]; 
  const double *q0c = &qc[0]; 
  const double *q0r = &qr[0]; 
  double *out0 = &out[0]; 

  out0[0] += J*D*(0.6708203932499369*q0r[2]+0.6708203932499369*q0l[2]-1.341640786499874*q0c[2]-1.190784930203603*q0r[1]+1.190784930203603*q0l[1]+0.9375*q0r[0]+0.9375*q0l[0]-1.875*q0c[0]); 
  out0[1] += J*D*(0.7382874503707888*q0r[2]-0.7382874503707888*q0l[2]-1.453125*q0r[1]-1.453125*q0l[1]-5.34375*q0c[1]+1.190784930203603*q0r[0]-1.190784930203603*q0l[0]); 
  out0[2] += J*D*((-0.140625*q0r[2])-0.140625*q0l[2]-6.28125*q0c[2]-0.3025768239224545*q0r[1]+0.3025768239224545*q0l[1]+0.4192627457812106*q0r[0]+0.4192627457812106*q0l[0]-0.8385254915624212*q0c[0]); 

  const double *q1l = &ql[3]; 
  const double *q1c = &qc[3]; 
  const double *q1r = &qr[3]; 
  double *out1 = &out[3]; 

  out1[0] += J*D*(0.6708203932499369*q1r[2]+0.6708203932499369*q1l[2]-1.341640786499874*q1c[2]-1.190784930203603*q1r[1]+1.190784930203603*q1l[1]+0.9375*q1r[0]+0.9375*q1l[0]-1.875*q1c[0]); 
  out1[1] += J*D*(0.7382874503707888*q1r[2]-0.7382874503707888*q1l[2]-1.453125*q1r[1]-1.453125*q1l[1]-5.34375*q1c[1]+1.190784930203603*q1r[0]-1.190784930203603*q1l[0]); 
  out1[2] += J*D*((-0.140625*q1r[2])-0.140625*q1l[2]-6.28125*q1c[2]-0.3025768239224545*q1r[1]+0.3025768239224545*q1l[1]+0.4192627457812106*q1r[0]+0.4192627457812106*q1l[0]-0.8385254915624212*q1c[0]); 

  const double *q2l = &ql[6]; 
  const double *q2c = &qc[6]; 
  const double *q2r = &qr[6]; 
  double *out2 = &out[6]; 

  out2[0] += J*D*(0.6708203932499369*q2r[2]+0.6708203932499369*q2l[2]-1.341640786499874*q2c[2]-1.190784930203603*q2r[1]+1.190784930203603*q2l[1]+0.9375*q2r[0]+0.9375*q2l[0]-1.875*q2c[0]); 
  out2[1] += J*D*(0.7382874503707888*q2r[2]-0.7382874503707888*q2l[2]-1.453125*q2r[1]-1.453125*q2l[1]-5.34375*q2c[1]+1.190784930203603*q2r[0]-1.190784930203603*q2l[0]); 
  out2[2] += J*D*((-0.140625*q2r[2])-0.140625*q2l[2]-6.28125*q2c[2]-0.3025768239224545*q2r[1]+0.3025768239224545*q2l[1]+0.4192627457812106*q2r[0]+0.4192627457812106*q2l[0]-0.8385254915624212*q2c[0]); 

  return 0.;

} 
