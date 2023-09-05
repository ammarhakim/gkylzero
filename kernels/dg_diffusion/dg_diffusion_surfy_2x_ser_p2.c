#include <gkyl_dg_diffusion_kernels.h> 
GKYL_CU_DH double dg_diffusion_surfy_2x_ser_p2(const double* w, const double* dx, double D, 
  const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: Cell-center coordinates
  // dxv[NDIM]: Cell spacing
  // D: Diffusion coefficient in the center cell
  // ql: Input field in the left cell
  // qc: Input field in the center cell
  // qr: Input field in the right cell
  // out: Incremented output

  const double dx1 = 2.0/dx[1]; 
  const double J = pow(dx1, 2.0);

  const double *q0l = &ql[0]; 
  const double *q0c = &qc[0]; 
  const double *q0r = &qr[0]; 
  double *out0 = &out[0]; 

  out0[0] += J*D*(0.6708203932499369*q0r[5]+0.6708203932499369*q0l[5]-1.341640786499874*q0c[5]-1.190784930203603*q0r[2]+1.190784930203603*q0l[2]+0.9375*q0r[0]+0.9375*q0l[0]-1.875*q0c[0]); 
  out0[1] += J*D*(0.6708203932499369*q0r[7]+0.6708203932499369*q0l[7]-1.341640786499874*q0c[7]-1.190784930203603*q0r[3]+1.190784930203603*q0l[3]+0.9375*q0r[1]+0.9375*q0l[1]-1.875*q0c[1]); 
  out0[2] += J*D*(0.7382874503707888*q0r[5]-0.7382874503707888*q0l[5]-1.453125*q0r[2]-1.453125*q0l[2]-5.34375*q0c[2]+1.190784930203603*q0r[0]-1.190784930203603*q0l[0]); 
  out0[3] += J*D*(0.7382874503707888*q0r[7]-0.7382874503707888*q0l[7]-1.453125*q0r[3]-1.453125*q0l[3]-5.34375*q0c[3]+1.190784930203603*q0r[1]-1.190784930203603*q0l[1]); 
  out0[4] += J*D*((-1.190784930203603*q0r[6])+1.190784930203603*q0l[6]+0.9375*q0r[4]+0.9375*q0l[4]-1.875*q0c[4]); 
  out0[5] += J*D*((-0.140625*q0r[5])-0.140625*q0l[5]-6.28125*q0c[5]-0.3025768239224545*q0r[2]+0.3025768239224545*q0l[2]+0.4192627457812106*q0r[0]+0.4192627457812106*q0l[0]-0.8385254915624212*q0c[0]); 
  out0[6] += J*D*((-1.453125*q0r[6])-1.453125*q0l[6]-5.34375*q0c[6]+1.190784930203603*q0r[4]-1.190784930203603*q0l[4]); 
  out0[7] += J*D*((-0.140625*q0r[7])-0.140625*q0l[7]-6.28125*q0c[7]-0.3025768239224544*q0r[3]+0.3025768239224544*q0l[3]+0.4192627457812105*q0r[1]+0.4192627457812105*q0l[1]-0.8385254915624211*q0c[1]); 

  return 0.;

} 
