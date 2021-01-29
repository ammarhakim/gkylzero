#include <gkyl_maxwell_kernels.h> 
gkyl_real maxwell_vol_1x_ser_p1(const gkyl_maxwell_inp *meq, const gkyl_real *w, const gkyl_real *dx, const gkyl_real *q, gkyl_real* restrict out) 
{ 
  const gkyl_real c2 = meq->c*meq->c, chi = meq->chi, gamma = meq->gamma; 
  const gkyl_real c2chi = c2*chi, c2gamma = c2*gamma; 
 
  const gkyl_real *ex = &q[0]; 
  const gkyl_real *ey = &q[2]; 
  const gkyl_real *ez = &q[4]; 
  const gkyl_real *bx = &q[6]; 
  const gkyl_real *by = &q[8]; 
  const gkyl_real *bz = &q[10]; 
  const gkyl_real *ph = &q[12]; 
  const gkyl_real *ps = &q[14]; 
 
  gkyl_real *outEx = &out[0]; 
  gkyl_real *outEy = &out[2]; 
  gkyl_real *outEz = &out[4]; 
  gkyl_real *outBx = &out[6]; 
  gkyl_real *outBy = &out[8]; 
  gkyl_real *outBz = &out[10]; 
  gkyl_real *outPh = &out[12]; 
  gkyl_real *outPs = &out[14]; 
 
  gkyl_real dx0 = 2.0/dx[0]; 

  outEx[1] += 1.732050807568877*ph[0]*c2chi*dx0; 

  outEy[1] += 1.732050807568877*bz[0]*c2*dx0; 

  outEz[1] += -1.732050807568877*by[0]*c2*dx0; 

  outBx[1] += 1.732050807568877*ps[0]*dx0*gamma; 

  outBy[1] += -1.732050807568877*ez[0]*dx0; 

  outBz[1] += 1.732050807568877*ey[0]*dx0; 

  outPh[1] += 1.732050807568877*ex[0]*chi*dx0; 

  outPs[1] += 1.732050807568877*bx[0]*c2gamma*dx0; 

  gkyl_real cflFreq = 0.0; 
  cflFreq += meq->c/dx[0]; 
  return cflFreq; 
} 
