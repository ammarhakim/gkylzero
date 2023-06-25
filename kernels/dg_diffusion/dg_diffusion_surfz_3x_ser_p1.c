#include <gkyl_dg_diffusion_kernels.h> 
GKYL_CU_DH void dg_diffusion_surfz_3x_ser_p1(const double* w, const double* dx, double D, 
  const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: Cell-center coordinates
  // dxv[NDIM]: Cell spacing
  // D: Diffusion coefficient in the center cell
  // ql: Input field in the left cell
  // qc: Input field in the center cell
  // qr: Input field in the right cell
  // out: Incremented output

  const double dx1 = 2.0/dx[2]; 
  const double J = pow(dx1, 2.0);

  const double *q0l = &ql[0]; 
  const double *q0c = &qc[0]; 
  const double *q0r = &qr[0]; 
  double *out0= &out[0]; 

  out0[0] += J*D*((-0.5412658773652741*q0r[3])+0.5412658773652741*q0l[3]+0.5625*q0r[0]+0.5625*q0l[0]-1.125*q0c[0]); 
  out0[1] += J*D*((-0.5412658773652741*q0r[5])+0.5412658773652741*q0l[5]+0.5625*q0r[1]+0.5625*q0l[1]-1.125*q0c[1]); 
  out0[2] += J*D*((-0.5412658773652741*q0r[6])+0.5412658773652741*q0l[6]+0.5625*q0r[2]+0.5625*q0l[2]-1.125*q0c[2]); 
  out0[3] += J*D*((-0.4375*q0r[3])-0.4375*q0l[3]-2.875*q0c[3]+0.5412658773652741*q0r[0]-0.5412658773652741*q0l[0]); 
  out0[4] += J*D*((-0.5412658773652741*q0r[7])+0.5412658773652741*q0l[7]+0.5625*q0r[4]+0.5625*q0l[4]-1.125*q0c[4]); 
  out0[5] += J*D*((-0.4375*q0r[5])-0.4375*q0l[5]-2.875*q0c[5]+0.5412658773652741*q0r[1]-0.5412658773652741*q0l[1]); 
  out0[6] += J*D*((-0.4375*q0r[6])-0.4375*q0l[6]-2.875*q0c[6]+0.5412658773652741*q0r[2]-0.5412658773652741*q0l[2]); 
  out0[7] += J*D*((-0.4375*q0r[7])-0.4375*q0l[7]-2.875*q0c[7]+0.5412658773652741*q0r[4]-0.5412658773652741*q0l[4]); 

} 
