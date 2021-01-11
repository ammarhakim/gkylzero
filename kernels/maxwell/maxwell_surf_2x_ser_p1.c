#include <gkyl_maxwell_kernels.h> 
double maxwell_surfx_2x_ser_p1(const gkyl_maxwell_inp *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double* restrict outl, double* restrict outr) 
{ 
  const double c = meq->c, chi = meq->chi, gamma = meq->gamma; 
  const double c2 = c*c; 
  const double c2chi = c2*chi, c2gamma = c2*gamma; 
 
  const double dxl1 = 2.0/dxl[0]; 
  const double dxr1 = 2.0/dxr[0]; 
  const double *exl = &ql[0]; 
  const double *eyl = &ql[4]; 
  const double *ezl = &ql[8]; 
  const double *bxl = &ql[12]; 
  const double *byl = &ql[16]; 
  const double *bzl = &ql[20]; 
  const double *phl = &ql[24]; 
  const double *psl = &ql[28]; 
 
  double *outExl = &outl[0]; 
  double *outEyl = &outl[4]; 
  double *outEzl = &outl[8]; 
  double *outBxl = &outl[12]; 
  double *outByl = &outl[16]; 
  double *outBzl = &outl[20]; 
  double *outPhl = &outl[24]; 
  double *outPsl = &outl[28]; 
 
  const double *exr = &qr[0]; 
  const double *eyr = &qr[4]; 
  const double *ezr = &qr[8]; 
  const double *bxr = &qr[12]; 
  const double *byr = &qr[16]; 
  const double *bzr = &qr[20]; 
  const double *phr = &qr[24]; 
  const double *psr = &qr[28]; 
 
  double *outExr = &outr[0]; 
  double *outEyr = &outr[4]; 
  double *outEzr = &outr[8]; 
  double *outBxr = &outr[12]; 
  double *outByr = &outr[16]; 
  double *outBzr = &outr[20]; 
  double *outPhr = &outr[24]; 
  double *outPsr = &outr[28]; 
 
  double incr[4]; 
 
  incr[0] = (0.4330127018922193*(exr[1]+exl[1])-0.25*exr[0]+0.25*exl[0])*c*chi+((-0.4330127018922193*phr[1])+0.4330127018922193*phl[1]+0.25*(phr[0]+phl[0]))*c2chi; 
  incr[1] = ((-0.75*(exr[1]+exl[1]))+0.4330127018922193*exr[0]-0.4330127018922193*exl[0])*c*chi+(0.75*phr[1]-0.75*phl[1]-0.4330127018922193*(phr[0]+phl[0]))*c2chi; 
  incr[2] = (0.4330127018922193*(exr[3]+exl[3])-0.25*exr[2]+0.25*exl[2])*c*chi+((-0.4330127018922193*phr[3])+0.4330127018922193*phl[3]+0.25*(phr[2]+phl[2]))*c2chi; 
  incr[3] = ((-0.75*(exr[3]+exl[3]))+0.4330127018922193*exr[2]-0.4330127018922193*exl[2])*c*chi+(0.75*phr[3]-0.75*phl[3]-0.4330127018922193*(phr[2]+phl[2]))*c2chi; 

  outExr[0] += incr[0]*dxr1; 
  outExr[1] += incr[1]*dxr1; 
  outExr[2] += incr[2]*dxr1; 
  outExr[3] += incr[3]*dxr1; 

  outExl[0] += -1.0*incr[0]*dxl1; 
  outExl[1] += incr[1]*dxl1; 
  outExl[2] += -1.0*incr[2]*dxl1; 
  outExl[3] += incr[3]*dxl1; 

  incr[0] = (0.4330127018922193*(eyr[1]+eyl[1])-0.25*eyr[0]+0.25*eyl[0])*tau+((-0.4330127018922193*bzr[1])+0.4330127018922193*bzl[1]+0.25*(bzr[0]+bzl[0]))*c2; 
  incr[1] = ((-0.75*(eyr[1]+eyl[1]))+0.4330127018922193*eyr[0]-0.4330127018922193*eyl[0])*tau+(0.75*bzr[1]-0.75*bzl[1]-0.4330127018922193*(bzr[0]+bzl[0]))*c2; 
  incr[2] = (0.4330127018922193*(eyr[3]+eyl[3])-0.25*eyr[2]+0.25*eyl[2])*tau+((-0.4330127018922193*bzr[3])+0.4330127018922193*bzl[3]+0.25*(bzr[2]+bzl[2]))*c2; 
  incr[3] = ((-0.75*(eyr[3]+eyl[3]))+0.4330127018922193*eyr[2]-0.4330127018922193*eyl[2])*tau+(0.75*bzr[3]-0.75*bzl[3]-0.4330127018922193*(bzr[2]+bzl[2]))*c2; 

  outEyr[0] += incr[0]*dxr1; 
  outEyr[1] += incr[1]*dxr1; 
  outEyr[2] += incr[2]*dxr1; 
  outEyr[3] += incr[3]*dxr1; 

  outEyl[0] += -1.0*incr[0]*dxl1; 
  outEyl[1] += incr[1]*dxl1; 
  outEyl[2] += -1.0*incr[2]*dxl1; 
  outEyl[3] += incr[3]*dxl1; 

  incr[0] = (0.4330127018922193*(ezr[1]+ezl[1])-0.25*ezr[0]+0.25*ezl[0])*tau+(0.4330127018922193*byr[1]-0.4330127018922193*byl[1]-0.25*(byr[0]+byl[0]))*c2; 
  incr[1] = ((-0.75*(ezr[1]+ezl[1]))+0.4330127018922193*ezr[0]-0.4330127018922193*ezl[0])*tau+((-0.75*byr[1])+0.75*byl[1]+0.4330127018922193*(byr[0]+byl[0]))*c2; 
  incr[2] = (0.4330127018922193*(ezr[3]+ezl[3])-0.25*ezr[2]+0.25*ezl[2])*tau+(0.4330127018922193*byr[3]-0.4330127018922193*byl[3]-0.25*(byr[2]+byl[2]))*c2; 
  incr[3] = ((-0.75*(ezr[3]+ezl[3]))+0.4330127018922193*ezr[2]-0.4330127018922193*ezl[2])*tau+((-0.75*byr[3])+0.75*byl[3]+0.4330127018922193*(byr[2]+byl[2]))*c2; 

  outEzr[0] += incr[0]*dxr1; 
  outEzr[1] += incr[1]*dxr1; 
  outEzr[2] += incr[2]*dxr1; 
  outEzr[3] += incr[3]*dxr1; 

  outEzl[0] += -1.0*incr[0]*dxl1; 
  outEzl[1] += incr[1]*dxl1; 
  outEzl[2] += -1.0*incr[2]*dxl1; 
  outEzl[3] += incr[3]*dxl1; 

  incr[0] = ((0.4330127018922193*(bxr[1]+bxl[1])-0.25*bxr[0]+0.25*bxl[0])*c-0.4330127018922193*psr[1]+0.4330127018922193*psl[1]+0.25*(psr[0]+psl[0]))*gamma; 
  incr[1] = (((-0.75*(bxr[1]+bxl[1]))+0.4330127018922193*bxr[0]-0.4330127018922193*bxl[0])*c+0.75*psr[1]-0.75*psl[1]-0.4330127018922193*(psr[0]+psl[0]))*gamma; 
  incr[2] = ((0.4330127018922193*(bxr[3]+bxl[3])-0.25*bxr[2]+0.25*bxl[2])*c-0.4330127018922193*psr[3]+0.4330127018922193*psl[3]+0.25*(psr[2]+psl[2]))*gamma; 
  incr[3] = (((-0.75*(bxr[3]+bxl[3]))+0.4330127018922193*bxr[2]-0.4330127018922193*bxl[2])*c+0.75*psr[3]-0.75*psl[3]-0.4330127018922193*(psr[2]+psl[2]))*gamma; 

  outBxr[0] += incr[0]*dxr1; 
  outBxr[1] += incr[1]*dxr1; 
  outBxr[2] += incr[2]*dxr1; 
  outBxr[3] += incr[3]*dxr1; 

  outBxl[0] += -1.0*incr[0]*dxl1; 
  outBxl[1] += incr[1]*dxl1; 
  outBxl[2] += -1.0*incr[2]*dxl1; 
  outBxl[3] += incr[3]*dxl1; 

  incr[0] = ((0.4330127018922193*(byr[1]+byl[1])-0.25*byr[0]+0.25*byl[0])*c2)/tau+0.4330127018922193*ezr[1]-0.4330127018922193*ezl[1]-0.25*(ezr[0]+ezl[0]); 
  incr[1] = (((-0.75*(byr[1]+byl[1]))+0.4330127018922193*byr[0]-0.4330127018922193*byl[0])*c2)/tau-0.75*ezr[1]+0.75*ezl[1]+0.4330127018922193*(ezr[0]+ezl[0]); 
  incr[2] = ((0.4330127018922193*(byr[3]+byl[3])-0.25*byr[2]+0.25*byl[2])*c2)/tau+0.4330127018922193*ezr[3]-0.4330127018922193*ezl[3]-0.25*(ezr[2]+ezl[2]); 
  incr[3] = (((-0.75*(byr[3]+byl[3]))+0.4330127018922193*byr[2]-0.4330127018922193*byl[2])*c2)/tau-0.75*ezr[3]+0.75*ezl[3]+0.4330127018922193*(ezr[2]+ezl[2]); 

  outByr[0] += incr[0]*dxr1; 
  outByr[1] += incr[1]*dxr1; 
  outByr[2] += incr[2]*dxr1; 
  outByr[3] += incr[3]*dxr1; 

  outByl[0] += -1.0*incr[0]*dxl1; 
  outByl[1] += incr[1]*dxl1; 
  outByl[2] += -1.0*incr[2]*dxl1; 
  outByl[3] += incr[3]*dxl1; 

  incr[0] = ((0.4330127018922193*(bzr[1]+bzl[1])-0.25*bzr[0]+0.25*bzl[0])*c2)/tau-0.4330127018922193*eyr[1]+0.4330127018922193*eyl[1]+0.25*(eyr[0]+eyl[0]); 
  incr[1] = (((-0.75*(bzr[1]+bzl[1]))+0.4330127018922193*bzr[0]-0.4330127018922193*bzl[0])*c2)/tau+0.75*eyr[1]-0.75*eyl[1]-0.4330127018922193*(eyr[0]+eyl[0]); 
  incr[2] = ((0.4330127018922193*(bzr[3]+bzl[3])-0.25*bzr[2]+0.25*bzl[2])*c2)/tau-0.4330127018922193*eyr[3]+0.4330127018922193*eyl[3]+0.25*(eyr[2]+eyl[2]); 
  incr[3] = (((-0.75*(bzr[3]+bzl[3]))+0.4330127018922193*bzr[2]-0.4330127018922193*bzl[2])*c2)/tau+0.75*eyr[3]-0.75*eyl[3]-0.4330127018922193*(eyr[2]+eyl[2]); 

  outBzr[0] += incr[0]*dxr1; 
  outBzr[1] += incr[1]*dxr1; 
  outBzr[2] += incr[2]*dxr1; 
  outBzr[3] += incr[3]*dxr1; 

  outBzl[0] += -1.0*incr[0]*dxl1; 
  outBzl[1] += incr[1]*dxl1; 
  outBzl[2] += -1.0*incr[2]*dxl1; 
  outBzl[3] += incr[3]*dxl1; 

  incr[0] = ((0.4330127018922193*(phr[1]+phl[1])-0.25*phr[0]+0.25*phl[0])*c-0.4330127018922193*exr[1]+0.4330127018922193*exl[1]+0.25*(exr[0]+exl[0]))*chi; 
  incr[1] = (((-0.75*(phr[1]+phl[1]))+0.4330127018922193*phr[0]-0.4330127018922193*phl[0])*c+0.75*exr[1]-0.75*exl[1]-0.4330127018922193*(exr[0]+exl[0]))*chi; 
  incr[2] = ((0.4330127018922193*(phr[3]+phl[3])-0.25*phr[2]+0.25*phl[2])*c-0.4330127018922193*exr[3]+0.4330127018922193*exl[3]+0.25*(exr[2]+exl[2]))*chi; 
  incr[3] = (((-0.75*(phr[3]+phl[3]))+0.4330127018922193*phr[2]-0.4330127018922193*phl[2])*c+0.75*exr[3]-0.75*exl[3]-0.4330127018922193*(exr[2]+exl[2]))*chi; 

  outPhr[0] += incr[0]*dxr1; 
  outPhr[1] += incr[1]*dxr1; 
  outPhr[2] += incr[2]*dxr1; 
  outPhr[3] += incr[3]*dxr1; 

  outPhl[0] += -1.0*incr[0]*dxl1; 
  outPhl[1] += incr[1]*dxl1; 
  outPhl[2] += -1.0*incr[2]*dxl1; 
  outPhl[3] += incr[3]*dxl1; 

  incr[0] = (0.4330127018922193*(psr[1]+psl[1])-0.25*psr[0]+0.25*psl[0])*c*gamma+((-0.4330127018922193*bxr[1])+0.4330127018922193*bxl[1]+0.25*(bxr[0]+bxl[0]))*c2gamma; 
  incr[1] = ((-0.75*(psr[1]+psl[1]))+0.4330127018922193*psr[0]-0.4330127018922193*psl[0])*c*gamma+(0.75*bxr[1]-0.75*bxl[1]-0.4330127018922193*(bxr[0]+bxl[0]))*c2gamma; 
  incr[2] = (0.4330127018922193*(psr[3]+psl[3])-0.25*psr[2]+0.25*psl[2])*c*gamma+((-0.4330127018922193*bxr[3])+0.4330127018922193*bxl[3]+0.25*(bxr[2]+bxl[2]))*c2gamma; 
  incr[3] = ((-0.75*(psr[3]+psl[3]))+0.4330127018922193*psr[2]-0.4330127018922193*psl[2])*c*gamma+(0.75*bxr[3]-0.75*bxl[3]-0.4330127018922193*(bxr[2]+bxl[2]))*c2gamma; 

  outPsr[0] += incr[0]*dxr1; 
  outPsr[1] += incr[1]*dxr1; 
  outPsr[2] += incr[2]*dxr1; 
  outPsr[3] += incr[3]*dxr1; 

  outPsl[0] += -1.0*incr[0]*dxl1; 
  outPsl[1] += incr[1]*dxl1; 
  outPsl[2] += -1.0*incr[2]*dxl1; 
  outPsl[3] += incr[3]*dxl1; 

  return fmax(c, tau); 
} 
double maxwell_surfy_2x_ser_p1(const gkyl_maxwell_inp *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double* restrict outl, double* restrict outr) 
{ 
  const double c = meq->c, chi = meq->chi, gamma = meq->gamma; 
  const double c2 = c*c; 
  const double c2chi = c2*chi, c2gamma = c2*gamma; 
 
  const double dxl1 = 2.0/dxl[1]; 
  const double dxr1 = 2.0/dxr[1]; 
  const double *exl = &ql[0]; 
  const double *eyl = &ql[4]; 
  const double *ezl = &ql[8]; 
  const double *bxl = &ql[12]; 
  const double *byl = &ql[16]; 
  const double *bzl = &ql[20]; 
  const double *phl = &ql[24]; 
  const double *psl = &ql[28]; 
 
  double *outExl = &outl[0]; 
  double *outEyl = &outl[4]; 
  double *outEzl = &outl[8]; 
  double *outBxl = &outl[12]; 
  double *outByl = &outl[16]; 
  double *outBzl = &outl[20]; 
  double *outPhl = &outl[24]; 
  double *outPsl = &outl[28]; 
 
  const double *exr = &qr[0]; 
  const double *eyr = &qr[4]; 
  const double *ezr = &qr[8]; 
  const double *bxr = &qr[12]; 
  const double *byr = &qr[16]; 
  const double *bzr = &qr[20]; 
  const double *phr = &qr[24]; 
  const double *psr = &qr[28]; 
 
  double *outExr = &outr[0]; 
  double *outEyr = &outr[4]; 
  double *outEzr = &outr[8]; 
  double *outBxr = &outr[12]; 
  double *outByr = &outr[16]; 
  double *outBzr = &outr[20]; 
  double *outPhr = &outr[24]; 
  double *outPsr = &outr[28]; 
 
  double incr[4]; 
 
  incr[0] = (0.4330127018922193*(exr[2]+exl[2])-0.25*exr[0]+0.25*exl[0])*tau+(0.4330127018922193*bzr[2]-0.4330127018922193*bzl[2]-0.25*(bzr[0]+bzl[0]))*c2; 
  incr[1] = (0.4330127018922193*(exr[3]+exl[3])-0.25*exr[1]+0.25*exl[1])*tau+(0.4330127018922193*bzr[3]-0.4330127018922193*bzl[3]-0.25*(bzr[1]+bzl[1]))*c2; 
  incr[2] = ((-0.75*(exr[2]+exl[2]))+0.4330127018922193*exr[0]-0.4330127018922193*exl[0])*tau+((-0.75*bzr[2])+0.75*bzl[2]+0.4330127018922193*(bzr[0]+bzl[0]))*c2; 
  incr[3] = ((-0.75*(exr[3]+exl[3]))+0.4330127018922193*exr[1]-0.4330127018922193*exl[1])*tau+((-0.75*bzr[3])+0.75*bzl[3]+0.4330127018922193*(bzr[1]+bzl[1]))*c2; 

  outExr[0] += incr[0]*dxr1; 
  outExr[1] += incr[1]*dxr1; 
  outExr[2] += incr[2]*dxr1; 
  outExr[3] += incr[3]*dxr1; 

  outExl[0] += -1.0*incr[0]*dxl1; 
  outExl[1] += -1.0*incr[1]*dxl1; 
  outExl[2] += incr[2]*dxl1; 
  outExl[3] += incr[3]*dxl1; 

  incr[0] = (0.4330127018922193*(eyr[2]+eyl[2])-0.25*eyr[0]+0.25*eyl[0])*c*chi+((-0.4330127018922193*phr[2])+0.4330127018922193*phl[2]+0.25*(phr[0]+phl[0]))*c2chi; 
  incr[1] = (0.4330127018922193*(eyr[3]+eyl[3])-0.25*eyr[1]+0.25*eyl[1])*c*chi+((-0.4330127018922193*phr[3])+0.4330127018922193*phl[3]+0.25*(phr[1]+phl[1]))*c2chi; 
  incr[2] = ((-0.75*(eyr[2]+eyl[2]))+0.4330127018922193*eyr[0]-0.4330127018922193*eyl[0])*c*chi+(0.75*phr[2]-0.75*phl[2]-0.4330127018922193*(phr[0]+phl[0]))*c2chi; 
  incr[3] = ((-0.75*(eyr[3]+eyl[3]))+0.4330127018922193*eyr[1]-0.4330127018922193*eyl[1])*c*chi+(0.75*phr[3]-0.75*phl[3]-0.4330127018922193*(phr[1]+phl[1]))*c2chi; 

  outEyr[0] += incr[0]*dxr1; 
  outEyr[1] += incr[1]*dxr1; 
  outEyr[2] += incr[2]*dxr1; 
  outEyr[3] += incr[3]*dxr1; 

  outEyl[0] += -1.0*incr[0]*dxl1; 
  outEyl[1] += -1.0*incr[1]*dxl1; 
  outEyl[2] += incr[2]*dxl1; 
  outEyl[3] += incr[3]*dxl1; 

  incr[0] = (0.4330127018922193*(ezr[2]+ezl[2])-0.25*ezr[0]+0.25*ezl[0])*tau+((-0.4330127018922193*bxr[2])+0.4330127018922193*bxl[2]+0.25*(bxr[0]+bxl[0]))*c2; 
  incr[1] = (0.4330127018922193*(ezr[3]+ezl[3])-0.25*ezr[1]+0.25*ezl[1])*tau+((-0.4330127018922193*bxr[3])+0.4330127018922193*bxl[3]+0.25*(bxr[1]+bxl[1]))*c2; 
  incr[2] = ((-0.75*(ezr[2]+ezl[2]))+0.4330127018922193*ezr[0]-0.4330127018922193*ezl[0])*tau+(0.75*bxr[2]-0.75*bxl[2]-0.4330127018922193*(bxr[0]+bxl[0]))*c2; 
  incr[3] = ((-0.75*(ezr[3]+ezl[3]))+0.4330127018922193*ezr[1]-0.4330127018922193*ezl[1])*tau+(0.75*bxr[3]-0.75*bxl[3]-0.4330127018922193*(bxr[1]+bxl[1]))*c2; 

  outEzr[0] += incr[0]*dxr1; 
  outEzr[1] += incr[1]*dxr1; 
  outEzr[2] += incr[2]*dxr1; 
  outEzr[3] += incr[3]*dxr1; 

  outEzl[0] += -1.0*incr[0]*dxl1; 
  outEzl[1] += -1.0*incr[1]*dxl1; 
  outEzl[2] += incr[2]*dxl1; 
  outEzl[3] += incr[3]*dxl1; 

  incr[0] = ((0.4330127018922193*(bxr[2]+bxl[2])-0.25*bxr[0]+0.25*bxl[0])*c2)/tau-0.4330127018922193*ezr[2]+0.4330127018922193*ezl[2]+0.25*(ezr[0]+ezl[0]); 
  incr[1] = ((0.4330127018922193*(bxr[3]+bxl[3])-0.25*bxr[1]+0.25*bxl[1])*c2)/tau-0.4330127018922193*ezr[3]+0.4330127018922193*ezl[3]+0.25*(ezr[1]+ezl[1]); 
  incr[2] = (((-0.75*(bxr[2]+bxl[2]))+0.4330127018922193*bxr[0]-0.4330127018922193*bxl[0])*c2)/tau+0.75*ezr[2]-0.75*ezl[2]-0.4330127018922193*(ezr[0]+ezl[0]); 
  incr[3] = (((-0.75*(bxr[3]+bxl[3]))+0.4330127018922193*bxr[1]-0.4330127018922193*bxl[1])*c2)/tau+0.75*ezr[3]-0.75*ezl[3]-0.4330127018922193*(ezr[1]+ezl[1]); 

  outBxr[0] += incr[0]*dxr1; 
  outBxr[1] += incr[1]*dxr1; 
  outBxr[2] += incr[2]*dxr1; 
  outBxr[3] += incr[3]*dxr1; 

  outBxl[0] += -1.0*incr[0]*dxl1; 
  outBxl[1] += -1.0*incr[1]*dxl1; 
  outBxl[2] += incr[2]*dxl1; 
  outBxl[3] += incr[3]*dxl1; 

  incr[0] = ((0.4330127018922193*(byr[2]+byl[2])-0.25*byr[0]+0.25*byl[0])*c-0.4330127018922193*psr[2]+0.4330127018922193*psl[2]+0.25*(psr[0]+psl[0]))*gamma; 
  incr[1] = ((0.4330127018922193*(byr[3]+byl[3])-0.25*byr[1]+0.25*byl[1])*c-0.4330127018922193*psr[3]+0.4330127018922193*psl[3]+0.25*(psr[1]+psl[1]))*gamma; 
  incr[2] = (((-0.75*(byr[2]+byl[2]))+0.4330127018922193*byr[0]-0.4330127018922193*byl[0])*c+0.75*psr[2]-0.75*psl[2]-0.4330127018922193*(psr[0]+psl[0]))*gamma; 
  incr[3] = (((-0.75*(byr[3]+byl[3]))+0.4330127018922193*byr[1]-0.4330127018922193*byl[1])*c+0.75*psr[3]-0.75*psl[3]-0.4330127018922193*(psr[1]+psl[1]))*gamma; 

  outByr[0] += incr[0]*dxr1; 
  outByr[1] += incr[1]*dxr1; 
  outByr[2] += incr[2]*dxr1; 
  outByr[3] += incr[3]*dxr1; 

  outByl[0] += -1.0*incr[0]*dxl1; 
  outByl[1] += -1.0*incr[1]*dxl1; 
  outByl[2] += incr[2]*dxl1; 
  outByl[3] += incr[3]*dxl1; 

  incr[0] = ((0.4330127018922193*(bzr[2]+bzl[2])-0.25*bzr[0]+0.25*bzl[0])*c2)/tau+0.4330127018922193*exr[2]-0.4330127018922193*exl[2]-0.25*(exr[0]+exl[0]); 
  incr[1] = ((0.4330127018922193*(bzr[3]+bzl[3])-0.25*bzr[1]+0.25*bzl[1])*c2)/tau+0.4330127018922193*exr[3]-0.4330127018922193*exl[3]-0.25*(exr[1]+exl[1]); 
  incr[2] = (((-0.75*(bzr[2]+bzl[2]))+0.4330127018922193*bzr[0]-0.4330127018922193*bzl[0])*c2)/tau-0.75*exr[2]+0.75*exl[2]+0.4330127018922193*(exr[0]+exl[0]); 
  incr[3] = (((-0.75*(bzr[3]+bzl[3]))+0.4330127018922193*bzr[1]-0.4330127018922193*bzl[1])*c2)/tau-0.75*exr[3]+0.75*exl[3]+0.4330127018922193*(exr[1]+exl[1]); 

  outBzr[0] += incr[0]*dxr1; 
  outBzr[1] += incr[1]*dxr1; 
  outBzr[2] += incr[2]*dxr1; 
  outBzr[3] += incr[3]*dxr1; 

  outBzl[0] += -1.0*incr[0]*dxl1; 
  outBzl[1] += -1.0*incr[1]*dxl1; 
  outBzl[2] += incr[2]*dxl1; 
  outBzl[3] += incr[3]*dxl1; 

  incr[0] = ((0.4330127018922193*(phr[2]+phl[2])-0.25*phr[0]+0.25*phl[0])*c-0.4330127018922193*eyr[2]+0.4330127018922193*eyl[2]+0.25*(eyr[0]+eyl[0]))*chi; 
  incr[1] = ((0.4330127018922193*(phr[3]+phl[3])-0.25*phr[1]+0.25*phl[1])*c-0.4330127018922193*eyr[3]+0.4330127018922193*eyl[3]+0.25*(eyr[1]+eyl[1]))*chi; 
  incr[2] = (((-0.75*(phr[2]+phl[2]))+0.4330127018922193*phr[0]-0.4330127018922193*phl[0])*c+0.75*eyr[2]-0.75*eyl[2]-0.4330127018922193*(eyr[0]+eyl[0]))*chi; 
  incr[3] = (((-0.75*(phr[3]+phl[3]))+0.4330127018922193*phr[1]-0.4330127018922193*phl[1])*c+0.75*eyr[3]-0.75*eyl[3]-0.4330127018922193*(eyr[1]+eyl[1]))*chi; 

  outPhr[0] += incr[0]*dxr1; 
  outPhr[1] += incr[1]*dxr1; 
  outPhr[2] += incr[2]*dxr1; 
  outPhr[3] += incr[3]*dxr1; 

  outPhl[0] += -1.0*incr[0]*dxl1; 
  outPhl[1] += -1.0*incr[1]*dxl1; 
  outPhl[2] += incr[2]*dxl1; 
  outPhl[3] += incr[3]*dxl1; 

  incr[0] = (0.4330127018922193*(psr[2]+psl[2])-0.25*psr[0]+0.25*psl[0])*c*gamma+((-0.4330127018922193*byr[2])+0.4330127018922193*byl[2]+0.25*(byr[0]+byl[0]))*c2gamma; 
  incr[1] = (0.4330127018922193*(psr[3]+psl[3])-0.25*psr[1]+0.25*psl[1])*c*gamma+((-0.4330127018922193*byr[3])+0.4330127018922193*byl[3]+0.25*(byr[1]+byl[1]))*c2gamma; 
  incr[2] = ((-0.75*(psr[2]+psl[2]))+0.4330127018922193*psr[0]-0.4330127018922193*psl[0])*c*gamma+(0.75*byr[2]-0.75*byl[2]-0.4330127018922193*(byr[0]+byl[0]))*c2gamma; 
  incr[3] = ((-0.75*(psr[3]+psl[3]))+0.4330127018922193*psr[1]-0.4330127018922193*psl[1])*c*gamma+(0.75*byr[3]-0.75*byl[3]-0.4330127018922193*(byr[1]+byl[1]))*c2gamma; 

  outPsr[0] += incr[0]*dxr1; 
  outPsr[1] += incr[1]*dxr1; 
  outPsr[2] += incr[2]*dxr1; 
  outPsr[3] += incr[3]*dxr1; 

  outPsl[0] += -1.0*incr[0]*dxl1; 
  outPsl[1] += -1.0*incr[1]*dxl1; 
  outPsl[2] += incr[2]*dxl1; 
  outPsl[3] += incr[3]*dxl1; 

  return fmax(c, tau); 
} 
