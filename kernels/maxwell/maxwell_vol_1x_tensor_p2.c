#include <gkyl_maxwell_kernels.h> 
GKYL_CU_DH double maxwell_vol_1x_tensor_p2(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *q, double* GKYL_RESTRICT out) 
{ 
  const double c2 = meq->c*meq->c, chi = meq->chi, gamma = meq->gamma; 
  const double c2chi = c2*chi, c2gamma = c2*gamma; 
 
  const double *ex = &q[0]; 
  const double *ey = &q[3]; 
  const double *ez = &q[6]; 
  const double *bx = &q[9]; 
  const double *by = &q[12]; 
  const double *bz = &q[15]; 
  const double *ph = &q[18]; 
  const double *ps = &q[21]; 
 
  double *outEx = &out[0]; 
  double *outEy = &out[3]; 
  double *outEz = &out[6]; 
  double *outBx = &out[9]; 
  double *outBy = &out[12]; 
  double *outBz = &out[15]; 
  double *outPh = &out[18]; 
  double *outPs = &out[21]; 
 
  double dx0 = 2.0/dx[0]; 

  outEx[1] += 1.732050807568877*ph[0]*c2chi*dx0; 
  outEx[2] += 3.872983346207417*ph[1]*c2chi*dx0; 

  outEy[1] += 1.732050807568877*bz[0]*c2*dx0; 
  outEy[2] += 3.872983346207417*bz[1]*c2*dx0; 

  outEz[1] += -1.732050807568877*by[0]*c2*dx0; 
  outEz[2] += -3.872983346207417*by[1]*c2*dx0; 

  outBx[1] += 1.732050807568877*ps[0]*dx0*gamma; 
  outBx[2] += 3.872983346207417*ps[1]*dx0*gamma; 

  outBy[1] += -1.732050807568877*ez[0]*dx0; 
  outBy[2] += -3.872983346207417*ez[1]*dx0; 

  outBz[1] += 1.732050807568877*ey[0]*dx0; 
  outBz[2] += 3.872983346207417*ey[1]*dx0; 

  outPh[1] += 1.732050807568877*ex[0]*chi*dx0; 
  outPh[2] += 3.872983346207417*ex[1]*chi*dx0; 

  outPs[1] += 1.732050807568877*bx[0]*c2gamma*dx0; 
  outPs[2] += 3.872983346207417*bx[1]*c2gamma*dx0; 

  double cflFreq = 0.0; 
  cflFreq += meq->c/dx[0]; 
  return 5.0*cflFreq; 
} 
