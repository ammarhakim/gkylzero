#include <gkyl_dg_diffusion_kernels.h> 
GKYL_CU_DH void dg_diffusion_pkpm_surfx_1x_ser_p1(const double* w, const double* dx, double D, 
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

  out0[0] += J*D*((-0.5412658773652741*q0r[1])+0.5412658773652741*q0l[1]+0.5625*q0r[0]+0.5625*q0l[0]-1.125*q0c[0]); 
  out0[1] += J*D*((-0.4375*q0r[1])-0.4375*q0l[1]-2.875*q0c[1]+0.5412658773652741*q0r[0]-0.5412658773652741*q0l[0]); 

  const double *q1l = &ql[2]; 
  const double *q1c = &qc[2]; 
  const double *q1r = &qr[2]; 
  double *out1= &out[2]; 

  out1[0] += J*D*((-0.5412658773652741*q1r[1])+0.5412658773652741*q1l[1]+0.5625*q1r[0]+0.5625*q1l[0]-1.125*q1c[0]); 
  out1[1] += J*D*((-0.4375*q1r[1])-0.4375*q1l[1]-2.875*q1c[1]+0.5412658773652741*q1r[0]-0.5412658773652741*q1l[0]); 

  const double *q2l = &ql[4]; 
  const double *q2c = &qc[4]; 
  const double *q2r = &qr[4]; 
  double *out2= &out[4]; 

  out2[0] += J*D*((-0.5412658773652741*q2r[1])+0.5412658773652741*q2l[1]+0.5625*q2r[0]+0.5625*q2l[0]-1.125*q2c[0]); 
  out2[1] += J*D*((-0.4375*q2r[1])-0.4375*q2l[1]-2.875*q2c[1]+0.5412658773652741*q2r[0]-0.5412658773652741*q2l[0]); 

} 
