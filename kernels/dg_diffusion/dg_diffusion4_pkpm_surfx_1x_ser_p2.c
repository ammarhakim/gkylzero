#include <gkyl_dg_diffusion_kernels.h> 
GKYL_CU_DH void dg_diffusion4_pkpm_surfx_1x_ser_p2(const double* w, const double* dx, double D, 
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
  const double J = -1.0*pow(dx1, 4.0);

  const double *q0l = &ql[0]; 
  const double *q0c = &qc[0]; 
  const double *q0r = &qr[0]; 
  double *out0= &out[0]; 

  out0[0] += J*D*((-6.708203932499369*q0r[2])-6.708203932499369*q0l[2]+13.41640786499874*q0c[2]+8.11898816047911*q0r[1]-8.11898816047911*q0l[1]-4.6875*q0r[0]-4.6875*q0l[0]+9.375*q0c[0]); 
  out0[1] += J*D*((-9.077304717673634*q0r[2])+9.077304717673634*q0l[2]+12.65625*q0r[1]+12.65625*q0l[1]+30.9375*q0c[1]-8.11898816047911*q0r[0]+8.11898816047911*q0l[0]); 
  out0[2] += J*D*((-0.65625*q0r[2])-0.65625*q0l[2]+40.6875*q0c[2]+4.720198453190289*q0r[1]-4.720198453190289*q0l[1]-4.192627457812106*q0r[0]-4.192627457812106*q0l[0]+8.385254915624213*q0c[0]); 

  const double *q1l = &ql[3]; 
  const double *q1c = &qc[3]; 
  const double *q1r = &qr[3]; 
  double *out1= &out[3]; 

  out1[0] += J*D*((-6.708203932499369*q1r[2])-6.708203932499369*q1l[2]+13.41640786499874*q1c[2]+8.11898816047911*q1r[1]-8.11898816047911*q1l[1]-4.6875*q1r[0]-4.6875*q1l[0]+9.375*q1c[0]); 
  out1[1] += J*D*((-9.077304717673634*q1r[2])+9.077304717673634*q1l[2]+12.65625*q1r[1]+12.65625*q1l[1]+30.9375*q1c[1]-8.11898816047911*q1r[0]+8.11898816047911*q1l[0]); 
  out1[2] += J*D*((-0.65625*q1r[2])-0.65625*q1l[2]+40.6875*q1c[2]+4.720198453190289*q1r[1]-4.720198453190289*q1l[1]-4.192627457812106*q1r[0]-4.192627457812106*q1l[0]+8.385254915624213*q1c[0]); 

  const double *q2l = &ql[6]; 
  const double *q2c = &qc[6]; 
  const double *q2r = &qr[6]; 
  double *out2= &out[6]; 

  out2[0] += J*D*((-6.708203932499369*q2r[2])-6.708203932499369*q2l[2]+13.41640786499874*q2c[2]+8.11898816047911*q2r[1]-8.11898816047911*q2l[1]-4.6875*q2r[0]-4.6875*q2l[0]+9.375*q2c[0]); 
  out2[1] += J*D*((-9.077304717673634*q2r[2])+9.077304717673634*q2l[2]+12.65625*q2r[1]+12.65625*q2l[1]+30.9375*q2c[1]-8.11898816047911*q2r[0]+8.11898816047911*q2l[0]); 
  out2[2] += J*D*((-0.65625*q2r[2])-0.65625*q2l[2]+40.6875*q2c[2]+4.720198453190289*q2r[1]-4.720198453190289*q2l[1]-4.192627457812106*q2r[0]-4.192627457812106*q2l[0]+8.385254915624213*q2c[0]); 

} 
