#include <gkyl_dg_diffusion_kernels.h> 
GKYL_CU_DH void dg_diffusion_surfx_1x_ser_p2(const double* w, const double* dx, double D, 
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
  double *out0= &out[0]; 

  out0[0] += J*D*(0.6708203932499369*q0r[2]+0.6708203932499369*q0l[2]-1.341640786499874*q0c[2]-1.190784930203603*q0r[1]+1.190784930203603*q0l[1]+0.9375*q0r[0]+0.9375*q0l[0]-1.875*q0c[0]); 
  out0[1] += J*D*(0.7382874503707888*q0r[2]-0.7382874503707888*q0l[2]-1.453125*q0r[1]-1.453125*q0l[1]-5.34375*q0c[1]+1.190784930203603*q0r[0]-1.190784930203603*q0l[0]); 
  out0[2] += J*D*((-0.140625*q0r[2])-0.140625*q0l[2]-6.28125*q0c[2]-0.3025768239224545*q0r[1]+0.3025768239224545*q0l[1]+0.4192627457812106*q0r[0]+0.4192627457812106*q0l[0]-0.8385254915624212*q0c[0]); 

} 
