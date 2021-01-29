#include <gkyl_maxwell_kernels.h> 
gkyl_real maxwell_vol_2x_ser_p2(const gkyl_maxwell_inp *meq, const gkyl_real *w, const gkyl_real *dx, const gkyl_real *q, gkyl_real* restrict out) 
{ 
  const gkyl_real c2 = meq->c*meq->c, chi = meq->chi, gamma = meq->gamma; 
  const gkyl_real c2chi = c2*chi, c2gamma = c2*gamma; 
 
  const gkyl_real *ex = &q[0]; 
  const gkyl_real *ey = &q[8]; 
  const gkyl_real *ez = &q[16]; 
  const gkyl_real *bx = &q[24]; 
  const gkyl_real *by = &q[32]; 
  const gkyl_real *bz = &q[40]; 
  const gkyl_real *ph = &q[48]; 
  const gkyl_real *ps = &q[56]; 
 
  gkyl_real *outEx = &out[0]; 
  gkyl_real *outEy = &out[8]; 
  gkyl_real *outEz = &out[16]; 
  gkyl_real *outBx = &out[24]; 
  gkyl_real *outBy = &out[32]; 
  gkyl_real *outBz = &out[40]; 
  gkyl_real *outPh = &out[48]; 
  gkyl_real *outPs = &out[56]; 
 
  gkyl_real dx0 = 2.0/dx[0]; 
  gkyl_real dx1 = 2.0/dx[1]; 

  outEx[1] += 1.732050807568877*ph[0]*c2chi*dx0; 
  outEx[2] += -1.732050807568877*bz[0]*c2*dx1; 
  outEx[3] += 1.732050807568877*ph[2]*c2chi*dx0-1.732050807568877*bz[1]*c2*dx1; 
  outEx[4] += 3.872983346207417*ph[1]*c2chi*dx0; 
  outEx[5] += -3.872983346207417*bz[2]*c2*dx1; 
  outEx[6] += 3.872983346207417*ph[3]*c2chi*dx0-1.732050807568877*bz[4]*c2*dx1; 
  outEx[7] += 1.732050807568877*ph[5]*c2chi*dx0-3.872983346207417*bz[3]*c2*dx1; 

  outEy[1] += 1.732050807568877*bz[0]*c2*dx0; 
  outEy[2] += 1.732050807568877*ph[0]*c2chi*dx1; 
  outEy[3] += 1.732050807568877*ph[1]*c2chi*dx1+1.732050807568877*bz[2]*c2*dx0; 
  outEy[4] += 3.872983346207417*bz[1]*c2*dx0; 
  outEy[5] += 3.872983346207417*ph[2]*c2chi*dx1; 
  outEy[6] += 1.732050807568877*ph[4]*c2chi*dx1+3.872983346207417*bz[3]*c2*dx0; 
  outEy[7] += 3.872983346207417*ph[3]*c2chi*dx1+1.732050807568877*bz[5]*c2*dx0; 

  outEz[1] += -1.732050807568877*by[0]*c2*dx0; 
  outEz[2] += 1.732050807568877*bx[0]*c2*dx1; 
  outEz[3] += 1.732050807568877*bx[1]*c2*dx1-1.732050807568877*by[2]*c2*dx0; 
  outEz[4] += -3.872983346207417*by[1]*c2*dx0; 
  outEz[5] += 3.872983346207417*bx[2]*c2*dx1; 
  outEz[6] += 1.732050807568877*bx[4]*c2*dx1-3.872983346207417*by[3]*c2*dx0; 
  outEz[7] += 3.872983346207417*bx[3]*c2*dx1-1.732050807568877*by[5]*c2*dx0; 

  outBx[1] += 1.732050807568877*ps[0]*dx0*gamma; 
  outBx[2] += 1.732050807568877*ez[0]*dx1; 
  outBx[3] += 1.732050807568877*ps[2]*dx0*gamma+1.732050807568877*ez[1]*dx1; 
  outBx[4] += 3.872983346207417*ps[1]*dx0*gamma; 
  outBx[5] += 3.872983346207417*ez[2]*dx1; 
  outBx[6] += 3.872983346207417*ps[3]*dx0*gamma+1.732050807568877*ez[4]*dx1; 
  outBx[7] += 1.732050807568877*ps[5]*dx0*gamma+3.872983346207417*ez[3]*dx1; 

  outBy[1] += -1.732050807568877*ez[0]*dx0; 
  outBy[2] += 1.732050807568877*ps[0]*dx1*gamma; 
  outBy[3] += 1.732050807568877*ps[1]*dx1*gamma-1.732050807568877*ez[2]*dx0; 
  outBy[4] += -3.872983346207417*ez[1]*dx0; 
  outBy[5] += 3.872983346207417*ps[2]*dx1*gamma; 
  outBy[6] += 1.732050807568877*ps[4]*dx1*gamma-3.872983346207417*ez[3]*dx0; 
  outBy[7] += 3.872983346207417*ps[3]*dx1*gamma-1.732050807568877*ez[5]*dx0; 

  outBz[1] += 1.732050807568877*ey[0]*dx0; 
  outBz[2] += -1.732050807568877*ex[0]*dx1; 
  outBz[3] += 1.732050807568877*ey[2]*dx0-1.732050807568877*ex[1]*dx1; 
  outBz[4] += 3.872983346207417*ey[1]*dx0; 
  outBz[5] += -3.872983346207417*ex[2]*dx1; 
  outBz[6] += 3.872983346207417*ey[3]*dx0-1.732050807568877*ex[4]*dx1; 
  outBz[7] += 1.732050807568877*ey[5]*dx0-3.872983346207417*ex[3]*dx1; 

  outPh[1] += 1.732050807568877*ex[0]*chi*dx0; 
  outPh[2] += 1.732050807568877*ey[0]*chi*dx1; 
  outPh[3] += 1.732050807568877*ey[1]*chi*dx1+1.732050807568877*ex[2]*chi*dx0; 
  outPh[4] += 3.872983346207417*ex[1]*chi*dx0; 
  outPh[5] += 3.872983346207417*ey[2]*chi*dx1; 
  outPh[6] += 1.732050807568877*ey[4]*chi*dx1+3.872983346207417*ex[3]*chi*dx0; 
  outPh[7] += 3.872983346207417*ey[3]*chi*dx1+1.732050807568877*ex[5]*chi*dx0; 

  outPs[1] += 1.732050807568877*bx[0]*c2gamma*dx0; 
  outPs[2] += 1.732050807568877*by[0]*c2gamma*dx1; 
  outPs[3] += 1.732050807568877*by[1]*c2gamma*dx1+1.732050807568877*bx[2]*c2gamma*dx0; 
  outPs[4] += 3.872983346207417*bx[1]*c2gamma*dx0; 
  outPs[5] += 3.872983346207417*by[2]*c2gamma*dx1; 
  outPs[6] += 1.732050807568877*by[4]*c2gamma*dx1+3.872983346207417*bx[3]*c2gamma*dx0; 
  outPs[7] += 3.872983346207417*by[3]*c2gamma*dx1+1.732050807568877*bx[5]*c2gamma*dx0; 

  gkyl_real cflFreq = 0.0; 
  cflFreq += meq->c/dx[0]; 
  cflFreq += meq->c/dx[1]; 
  return cflFreq; 
} 
