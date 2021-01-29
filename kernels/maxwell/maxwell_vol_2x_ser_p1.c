#include <gkyl_maxwell_kernels.h> 
gkyl_real maxwell_vol_2x_ser_p1(const gkyl_maxwell_inp *meq, const gkyl_real *w, const gkyl_real *dx, const gkyl_real *q, gkyl_real* restrict out) 
{ 
  const gkyl_real c2 = meq->c*meq->c, chi = meq->chi, gamma = meq->gamma; 
  const gkyl_real c2chi = c2*chi, c2gamma = c2*gamma; 
 
  const gkyl_real *ex = &q[0]; 
  const gkyl_real *ey = &q[4]; 
  const gkyl_real *ez = &q[8]; 
  const gkyl_real *bx = &q[12]; 
  const gkyl_real *by = &q[16]; 
  const gkyl_real *bz = &q[20]; 
  const gkyl_real *ph = &q[24]; 
  const gkyl_real *ps = &q[28]; 
 
  gkyl_real *outEx = &out[0]; 
  gkyl_real *outEy = &out[4]; 
  gkyl_real *outEz = &out[8]; 
  gkyl_real *outBx = &out[12]; 
  gkyl_real *outBy = &out[16]; 
  gkyl_real *outBz = &out[20]; 
  gkyl_real *outPh = &out[24]; 
  gkyl_real *outPs = &out[28]; 
 
  gkyl_real dx0 = 2.0/dx[0]; 
  gkyl_real dx1 = 2.0/dx[1]; 

  outEx[1] += 1.732050807568877*ph[0]*c2chi*dx0; 
  outEx[2] += -1.732050807568877*bz[0]*c2*dx1; 
  outEx[3] += 1.732050807568877*ph[2]*c2chi*dx0-1.732050807568877*bz[1]*c2*dx1; 

  outEy[1] += 1.732050807568877*bz[0]*c2*dx0; 
  outEy[2] += 1.732050807568877*ph[0]*c2chi*dx1; 
  outEy[3] += 1.732050807568877*ph[1]*c2chi*dx1+1.732050807568877*bz[2]*c2*dx0; 

  outEz[1] += -1.732050807568877*by[0]*c2*dx0; 
  outEz[2] += 1.732050807568877*bx[0]*c2*dx1; 
  outEz[3] += 1.732050807568877*bx[1]*c2*dx1-1.732050807568877*by[2]*c2*dx0; 

  outBx[1] += 1.732050807568877*ps[0]*dx0*gamma; 
  outBx[2] += 1.732050807568877*ez[0]*dx1; 
  outBx[3] += 1.732050807568877*ps[2]*dx0*gamma+1.732050807568877*ez[1]*dx1; 

  outBy[1] += -1.732050807568877*ez[0]*dx0; 
  outBy[2] += 1.732050807568877*ps[0]*dx1*gamma; 
  outBy[3] += 1.732050807568877*ps[1]*dx1*gamma-1.732050807568877*ez[2]*dx0; 

  outBz[1] += 1.732050807568877*ey[0]*dx0; 
  outBz[2] += -1.732050807568877*ex[0]*dx1; 
  outBz[3] += 1.732050807568877*ey[2]*dx0-1.732050807568877*ex[1]*dx1; 

  outPh[1] += 1.732050807568877*ex[0]*chi*dx0; 
  outPh[2] += 1.732050807568877*ey[0]*chi*dx1; 
  outPh[3] += 1.732050807568877*ey[1]*chi*dx1+1.732050807568877*ex[2]*chi*dx0; 

  outPs[1] += 1.732050807568877*bx[0]*c2gamma*dx0; 
  outPs[2] += 1.732050807568877*by[0]*c2gamma*dx1; 
  outPs[3] += 1.732050807568877*by[1]*c2gamma*dx1+1.732050807568877*bx[2]*c2gamma*dx0; 

  gkyl_real cflFreq = 0.0; 
  cflFreq += meq->c/dx[0]; 
  cflFreq += meq->c/dx[1]; 
  return cflFreq; 
} 
