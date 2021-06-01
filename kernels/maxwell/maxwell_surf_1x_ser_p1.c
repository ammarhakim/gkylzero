#include <gkyl_maxwell_kernels.h> 
double maxwell_surfx_1x_ser_p1(const gkyl_maxwell_inp *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double* GKYL_RESTRICT outl, double* GKYL_RESTRICT outr) 
{ 
  const double c = meq->c, chi = meq->chi, gamma = meq->gamma; 
  const double c2 = c*c; 
  const double c2chi = c2*chi, c2gamma = c2*gamma; 
 
  const double dxl1 = 2.0/dxl[0]; 
  const double dxr1 = 2.0/dxr[0]; 
  const double *exl = &ql[0]; 
  const double *eyl = &ql[2]; 
  const double *ezl = &ql[4]; 
  const double *bxl = &ql[6]; 
  const double *byl = &ql[8]; 
  const double *bzl = &ql[10]; 
  const double *phl = &ql[12]; 
  const double *psl = &ql[14]; 
 
  double *outExl = &outl[0]; 
  double *outEyl = &outl[2]; 
  double *outEzl = &outl[4]; 
  double *outBxl = &outl[6]; 
  double *outByl = &outl[8]; 
  double *outBzl = &outl[10]; 
  double *outPhl = &outl[12]; 
  double *outPsl = &outl[14]; 
 
  const double *exr = &qr[0]; 
  const double *eyr = &qr[2]; 
  const double *ezr = &qr[4]; 
  const double *bxr = &qr[6]; 
  const double *byr = &qr[8]; 
  const double *bzr = &qr[10]; 
  const double *phr = &qr[12]; 
  const double *psr = &qr[14]; 
 
  double *outExr = &outr[0]; 
  double *outEyr = &outr[2]; 
  double *outEzr = &outr[4]; 
  double *outBxr = &outr[6]; 
  double *outByr = &outr[8]; 
  double *outBzr = &outr[10]; 
  double *outPhr = &outr[12]; 
  double *outPsr = &outr[14]; 
 
  double incr[2]; 
 
  incr[0] = (0.4330127018922193*(exr[1]+exl[1])-0.25*exr[0]+0.25*exl[0])*c*chi+((-0.4330127018922193*phr[1])+0.4330127018922193*phl[1]+0.25*(phr[0]+phl[0]))*c2chi; 
  incr[1] = ((-0.75*(exr[1]+exl[1]))+0.4330127018922193*exr[0]-0.4330127018922193*exl[0])*c*chi+(0.75*phr[1]-0.75*phl[1]-0.4330127018922193*(phr[0]+phl[0]))*c2chi; 

  outExr[0] += incr[0]*dxr1; 
  outExr[1] += incr[1]*dxr1; 

  outExl[0] += -1.0*incr[0]*dxl1; 
  outExl[1] += incr[1]*dxl1; 

  incr[0] = (0.4330127018922193*(eyr[1]+eyl[1])-0.25*eyr[0]+0.25*eyl[0])*tau+((-0.4330127018922193*bzr[1])+0.4330127018922193*bzl[1]+0.25*(bzr[0]+bzl[0]))*c2; 
  incr[1] = ((-0.75*(eyr[1]+eyl[1]))+0.4330127018922193*eyr[0]-0.4330127018922193*eyl[0])*tau+(0.75*bzr[1]-0.75*bzl[1]-0.4330127018922193*(bzr[0]+bzl[0]))*c2; 

  outEyr[0] += incr[0]*dxr1; 
  outEyr[1] += incr[1]*dxr1; 

  outEyl[0] += -1.0*incr[0]*dxl1; 
  outEyl[1] += incr[1]*dxl1; 

  incr[0] = (0.4330127018922193*(ezr[1]+ezl[1])-0.25*ezr[0]+0.25*ezl[0])*tau+(0.4330127018922193*byr[1]-0.4330127018922193*byl[1]-0.25*(byr[0]+byl[0]))*c2; 
  incr[1] = ((-0.75*(ezr[1]+ezl[1]))+0.4330127018922193*ezr[0]-0.4330127018922193*ezl[0])*tau+((-0.75*byr[1])+0.75*byl[1]+0.4330127018922193*(byr[0]+byl[0]))*c2; 

  outEzr[0] += incr[0]*dxr1; 
  outEzr[1] += incr[1]*dxr1; 

  outEzl[0] += -1.0*incr[0]*dxl1; 
  outEzl[1] += incr[1]*dxl1; 

  incr[0] = ((0.4330127018922193*(bxr[1]+bxl[1])-0.25*bxr[0]+0.25*bxl[0])*c-0.4330127018922193*psr[1]+0.4330127018922193*psl[1]+0.25*(psr[0]+psl[0]))*gamma; 
  incr[1] = (((-0.75*(bxr[1]+bxl[1]))+0.4330127018922193*bxr[0]-0.4330127018922193*bxl[0])*c+0.75*psr[1]-0.75*psl[1]-0.4330127018922193*(psr[0]+psl[0]))*gamma; 

  outBxr[0] += incr[0]*dxr1; 
  outBxr[1] += incr[1]*dxr1; 

  outBxl[0] += -1.0*incr[0]*dxl1; 
  outBxl[1] += incr[1]*dxl1; 

  incr[0] = ((0.4330127018922193*(byr[1]+byl[1])-0.25*byr[0]+0.25*byl[0])*c2)/tau+0.4330127018922193*ezr[1]-0.4330127018922193*ezl[1]-0.25*(ezr[0]+ezl[0]); 
  incr[1] = (((-0.75*(byr[1]+byl[1]))+0.4330127018922193*byr[0]-0.4330127018922193*byl[0])*c2)/tau-0.75*ezr[1]+0.75*ezl[1]+0.4330127018922193*(ezr[0]+ezl[0]); 

  outByr[0] += incr[0]*dxr1; 
  outByr[1] += incr[1]*dxr1; 

  outByl[0] += -1.0*incr[0]*dxl1; 
  outByl[1] += incr[1]*dxl1; 

  incr[0] = ((0.4330127018922193*(bzr[1]+bzl[1])-0.25*bzr[0]+0.25*bzl[0])*c2)/tau-0.4330127018922193*eyr[1]+0.4330127018922193*eyl[1]+0.25*(eyr[0]+eyl[0]); 
  incr[1] = (((-0.75*(bzr[1]+bzl[1]))+0.4330127018922193*bzr[0]-0.4330127018922193*bzl[0])*c2)/tau+0.75*eyr[1]-0.75*eyl[1]-0.4330127018922193*(eyr[0]+eyl[0]); 

  outBzr[0] += incr[0]*dxr1; 
  outBzr[1] += incr[1]*dxr1; 

  outBzl[0] += -1.0*incr[0]*dxl1; 
  outBzl[1] += incr[1]*dxl1; 

  incr[0] = ((0.4330127018922193*(phr[1]+phl[1])-0.25*phr[0]+0.25*phl[0])*c-0.4330127018922193*exr[1]+0.4330127018922193*exl[1]+0.25*(exr[0]+exl[0]))*chi; 
  incr[1] = (((-0.75*(phr[1]+phl[1]))+0.4330127018922193*phr[0]-0.4330127018922193*phl[0])*c+0.75*exr[1]-0.75*exl[1]-0.4330127018922193*(exr[0]+exl[0]))*chi; 

  outPhr[0] += incr[0]*dxr1; 
  outPhr[1] += incr[1]*dxr1; 

  outPhl[0] += -1.0*incr[0]*dxl1; 
  outPhl[1] += incr[1]*dxl1; 

  incr[0] = (0.4330127018922193*(psr[1]+psl[1])-0.25*psr[0]+0.25*psl[0])*c*gamma+((-0.4330127018922193*bxr[1])+0.4330127018922193*bxl[1]+0.25*(bxr[0]+bxl[0]))*c2gamma; 
  incr[1] = ((-0.75*(psr[1]+psl[1]))+0.4330127018922193*psr[0]-0.4330127018922193*psl[0])*c*gamma+(0.75*bxr[1]-0.75*bxl[1]-0.4330127018922193*(bxr[0]+bxl[0]))*c2gamma; 

  outPsr[0] += incr[0]*dxr1; 
  outPsr[1] += incr[1]*dxr1; 

  outPsl[0] += -1.0*incr[0]*dxl1; 
  outPsl[1] += incr[1]*dxl1; 

  return fmax(c, tau); 
} 
