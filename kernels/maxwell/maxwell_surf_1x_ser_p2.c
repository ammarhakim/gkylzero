#include <gkyl_maxwell_kernels.h> 
gkyl_real maxwell_surfx_1x_ser_p2(const gkyl_maxwell_inp *meq, const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxl, const gkyl_real *dxr, const gkyl_real tau, const gkyl_real *ql, const gkyl_real *qr, gkyl_real* restrict outl, gkyl_real* restrict outr) 
{ 
  const gkyl_real c = meq->c, chi = meq->chi, gamma = meq->gamma; 
  const gkyl_real c2 = c*c; 
  const gkyl_real c2chi = c2*chi, c2gamma = c2*gamma; 
 
  const gkyl_real dxl1 = 2.0/dxl[0]; 
  const gkyl_real dxr1 = 2.0/dxr[0]; 
  const gkyl_real *exl = &ql[0]; 
  const gkyl_real *eyl = &ql[3]; 
  const gkyl_real *ezl = &ql[6]; 
  const gkyl_real *bxl = &ql[9]; 
  const gkyl_real *byl = &ql[12]; 
  const gkyl_real *bzl = &ql[15]; 
  const gkyl_real *phl = &ql[18]; 
  const gkyl_real *psl = &ql[21]; 
 
  gkyl_real *outExl = &outl[0]; 
  gkyl_real *outEyl = &outl[3]; 
  gkyl_real *outEzl = &outl[6]; 
  gkyl_real *outBxl = &outl[9]; 
  gkyl_real *outByl = &outl[12]; 
  gkyl_real *outBzl = &outl[15]; 
  gkyl_real *outPhl = &outl[18]; 
  gkyl_real *outPsl = &outl[21]; 
 
  const gkyl_real *exr = &qr[0]; 
  const gkyl_real *eyr = &qr[3]; 
  const gkyl_real *ezr = &qr[6]; 
  const gkyl_real *bxr = &qr[9]; 
  const gkyl_real *byr = &qr[12]; 
  const gkyl_real *bzr = &qr[15]; 
  const gkyl_real *phr = &qr[18]; 
  const gkyl_real *psr = &qr[21]; 
 
  gkyl_real *outExr = &outr[0]; 
  gkyl_real *outEyr = &outr[3]; 
  gkyl_real *outEzr = &outr[6]; 
  gkyl_real *outBxr = &outr[9]; 
  gkyl_real *outByr = &outr[12]; 
  gkyl_real *outBzr = &outr[15]; 
  gkyl_real *outPhr = &outr[18]; 
  gkyl_real *outPsr = &outr[21]; 
 
  gkyl_real incr[3]; 
 
  incr[0] = ((-0.5590169943749475*exr[2])+0.5590169943749475*exl[2]+0.4330127018922193*(exr[1]+exl[1])-0.25*exr[0]+0.25*exl[0])*c*chi+(0.5590169943749475*(phr[2]+phl[2])-0.4330127018922193*phr[1]+0.4330127018922193*phl[1]+0.25*(phr[0]+phl[0]))*c2chi; 
  incr[1] = (0.9682458365518543*exr[2]-0.9682458365518543*exl[2]-0.75*(exr[1]+exl[1])+0.4330127018922193*exr[0]-0.4330127018922193*exl[0])*c*chi+((-0.9682458365518543*(phr[2]+phl[2]))+0.75*phr[1]-0.75*phl[1]-0.4330127018922193*(phr[0]+phl[0]))*c2chi; 
  incr[2] = ((-1.25*exr[2])+1.25*exl[2]+0.9682458365518543*(exr[1]+exl[1])-0.5590169943749475*exr[0]+0.5590169943749475*exl[0])*c*chi+(1.25*(phr[2]+phl[2])-0.9682458365518543*phr[1]+0.9682458365518543*phl[1]+0.5590169943749475*(phr[0]+phl[0]))*c2chi; 

  outExr[0] += incr[0]*dxr1; 
  outExr[1] += incr[1]*dxr1; 
  outExr[2] += incr[2]*dxr1; 

  outExl[0] += -1.0*incr[0]*dxl1; 
  outExl[1] += incr[1]*dxl1; 
  outExl[2] += -1.0*incr[2]*dxl1; 

  incr[0] = ((-0.5590169943749475*eyr[2])+0.5590169943749475*eyl[2]+0.4330127018922193*(eyr[1]+eyl[1])-0.25*eyr[0]+0.25*eyl[0])*tau+(0.5590169943749475*(bzr[2]+bzl[2])-0.4330127018922193*bzr[1]+0.4330127018922193*bzl[1]+0.25*(bzr[0]+bzl[0]))*c2; 
  incr[1] = (0.9682458365518543*eyr[2]-0.9682458365518543*eyl[2]-0.75*(eyr[1]+eyl[1])+0.4330127018922193*eyr[0]-0.4330127018922193*eyl[0])*tau+((-0.9682458365518543*(bzr[2]+bzl[2]))+0.75*bzr[1]-0.75*bzl[1]-0.4330127018922193*(bzr[0]+bzl[0]))*c2; 
  incr[2] = ((-1.25*eyr[2])+1.25*eyl[2]+0.9682458365518543*(eyr[1]+eyl[1])-0.5590169943749475*eyr[0]+0.5590169943749475*eyl[0])*tau+(1.25*(bzr[2]+bzl[2])-0.9682458365518543*bzr[1]+0.9682458365518543*bzl[1]+0.5590169943749475*(bzr[0]+bzl[0]))*c2; 

  outEyr[0] += incr[0]*dxr1; 
  outEyr[1] += incr[1]*dxr1; 
  outEyr[2] += incr[2]*dxr1; 

  outEyl[0] += -1.0*incr[0]*dxl1; 
  outEyl[1] += incr[1]*dxl1; 
  outEyl[2] += -1.0*incr[2]*dxl1; 

  incr[0] = ((-0.5590169943749475*ezr[2])+0.5590169943749475*ezl[2]+0.4330127018922193*(ezr[1]+ezl[1])-0.25*ezr[0]+0.25*ezl[0])*tau+((-0.5590169943749475*(byr[2]+byl[2]))+0.4330127018922193*byr[1]-0.4330127018922193*byl[1]-0.25*(byr[0]+byl[0]))*c2; 
  incr[1] = (0.9682458365518543*ezr[2]-0.9682458365518543*ezl[2]-0.75*(ezr[1]+ezl[1])+0.4330127018922193*ezr[0]-0.4330127018922193*ezl[0])*tau+(0.9682458365518543*(byr[2]+byl[2])-0.75*byr[1]+0.75*byl[1]+0.4330127018922193*(byr[0]+byl[0]))*c2; 
  incr[2] = ((-1.25*ezr[2])+1.25*ezl[2]+0.9682458365518543*(ezr[1]+ezl[1])-0.5590169943749475*ezr[0]+0.5590169943749475*ezl[0])*tau+((-1.25*(byr[2]+byl[2]))+0.9682458365518543*byr[1]-0.9682458365518543*byl[1]-0.5590169943749475*(byr[0]+byl[0]))*c2; 

  outEzr[0] += incr[0]*dxr1; 
  outEzr[1] += incr[1]*dxr1; 
  outEzr[2] += incr[2]*dxr1; 

  outEzl[0] += -1.0*incr[0]*dxl1; 
  outEzl[1] += incr[1]*dxl1; 
  outEzl[2] += -1.0*incr[2]*dxl1; 

  incr[0] = (((-0.5590169943749475*bxr[2])+0.5590169943749475*bxl[2]+0.4330127018922193*(bxr[1]+bxl[1])-0.25*bxr[0]+0.25*bxl[0])*c+0.5590169943749475*(psr[2]+psl[2])-0.4330127018922193*psr[1]+0.4330127018922193*psl[1]+0.25*(psr[0]+psl[0]))*gamma; 
  incr[1] = ((0.9682458365518543*bxr[2]-0.9682458365518543*bxl[2]-0.75*(bxr[1]+bxl[1])+0.4330127018922193*bxr[0]-0.4330127018922193*bxl[0])*c-0.9682458365518543*(psr[2]+psl[2])+0.75*psr[1]-0.75*psl[1]-0.4330127018922193*(psr[0]+psl[0]))*gamma; 
  incr[2] = (((-1.25*bxr[2])+1.25*bxl[2]+0.9682458365518543*(bxr[1]+bxl[1])-0.5590169943749475*bxr[0]+0.5590169943749475*bxl[0])*c+1.25*(psr[2]+psl[2])-0.9682458365518543*psr[1]+0.9682458365518543*psl[1]+0.5590169943749475*(psr[0]+psl[0]))*gamma; 

  outBxr[0] += incr[0]*dxr1; 
  outBxr[1] += incr[1]*dxr1; 
  outBxr[2] += incr[2]*dxr1; 

  outBxl[0] += -1.0*incr[0]*dxl1; 
  outBxl[1] += incr[1]*dxl1; 
  outBxl[2] += -1.0*incr[2]*dxl1; 

  incr[0] = (((-0.5590169943749475*byr[2])+0.5590169943749475*byl[2]+0.4330127018922193*(byr[1]+byl[1])-0.25*byr[0]+0.25*byl[0])*c2)/tau-0.5590169943749475*(ezr[2]+ezl[2])+0.4330127018922193*ezr[1]-0.4330127018922193*ezl[1]-0.25*(ezr[0]+ezl[0]); 
  incr[1] = ((0.9682458365518543*byr[2]-0.9682458365518543*byl[2]-0.75*(byr[1]+byl[1])+0.4330127018922193*byr[0]-0.4330127018922193*byl[0])*c2)/tau+0.9682458365518543*(ezr[2]+ezl[2])-0.75*ezr[1]+0.75*ezl[1]+0.4330127018922193*(ezr[0]+ezl[0]); 
  incr[2] = (((-1.25*byr[2])+1.25*byl[2]+0.9682458365518543*(byr[1]+byl[1])-0.5590169943749475*byr[0]+0.5590169943749475*byl[0])*c2)/tau-1.25*(ezr[2]+ezl[2])+0.9682458365518543*ezr[1]-0.9682458365518543*ezl[1]-0.5590169943749475*(ezr[0]+ezl[0]); 

  outByr[0] += incr[0]*dxr1; 
  outByr[1] += incr[1]*dxr1; 
  outByr[2] += incr[2]*dxr1; 

  outByl[0] += -1.0*incr[0]*dxl1; 
  outByl[1] += incr[1]*dxl1; 
  outByl[2] += -1.0*incr[2]*dxl1; 

  incr[0] = (((-0.5590169943749475*bzr[2])+0.5590169943749475*bzl[2]+0.4330127018922193*(bzr[1]+bzl[1])-0.25*bzr[0]+0.25*bzl[0])*c2)/tau+0.5590169943749475*(eyr[2]+eyl[2])-0.4330127018922193*eyr[1]+0.4330127018922193*eyl[1]+0.25*(eyr[0]+eyl[0]); 
  incr[1] = ((0.9682458365518543*bzr[2]-0.9682458365518543*bzl[2]-0.75*(bzr[1]+bzl[1])+0.4330127018922193*bzr[0]-0.4330127018922193*bzl[0])*c2)/tau-0.9682458365518543*(eyr[2]+eyl[2])+0.75*eyr[1]-0.75*eyl[1]-0.4330127018922193*(eyr[0]+eyl[0]); 
  incr[2] = (((-1.25*bzr[2])+1.25*bzl[2]+0.9682458365518543*(bzr[1]+bzl[1])-0.5590169943749475*bzr[0]+0.5590169943749475*bzl[0])*c2)/tau+1.25*(eyr[2]+eyl[2])-0.9682458365518543*eyr[1]+0.9682458365518543*eyl[1]+0.5590169943749475*(eyr[0]+eyl[0]); 

  outBzr[0] += incr[0]*dxr1; 
  outBzr[1] += incr[1]*dxr1; 
  outBzr[2] += incr[2]*dxr1; 

  outBzl[0] += -1.0*incr[0]*dxl1; 
  outBzl[1] += incr[1]*dxl1; 
  outBzl[2] += -1.0*incr[2]*dxl1; 

  incr[0] = (((-0.5590169943749475*phr[2])+0.5590169943749475*phl[2]+0.4330127018922193*(phr[1]+phl[1])-0.25*phr[0]+0.25*phl[0])*c+0.5590169943749475*(exr[2]+exl[2])-0.4330127018922193*exr[1]+0.4330127018922193*exl[1]+0.25*(exr[0]+exl[0]))*chi; 
  incr[1] = ((0.9682458365518543*phr[2]-0.9682458365518543*phl[2]-0.75*(phr[1]+phl[1])+0.4330127018922193*phr[0]-0.4330127018922193*phl[0])*c-0.9682458365518543*(exr[2]+exl[2])+0.75*exr[1]-0.75*exl[1]-0.4330127018922193*(exr[0]+exl[0]))*chi; 
  incr[2] = (((-1.25*phr[2])+1.25*phl[2]+0.9682458365518543*(phr[1]+phl[1])-0.5590169943749475*phr[0]+0.5590169943749475*phl[0])*c+1.25*(exr[2]+exl[2])-0.9682458365518543*exr[1]+0.9682458365518543*exl[1]+0.5590169943749475*(exr[0]+exl[0]))*chi; 

  outPhr[0] += incr[0]*dxr1; 
  outPhr[1] += incr[1]*dxr1; 
  outPhr[2] += incr[2]*dxr1; 

  outPhl[0] += -1.0*incr[0]*dxl1; 
  outPhl[1] += incr[1]*dxl1; 
  outPhl[2] += -1.0*incr[2]*dxl1; 

  incr[0] = ((-0.5590169943749475*psr[2])+0.5590169943749475*psl[2]+0.4330127018922193*(psr[1]+psl[1])-0.25*psr[0]+0.25*psl[0])*c*gamma+(0.5590169943749475*(bxr[2]+bxl[2])-0.4330127018922193*bxr[1]+0.4330127018922193*bxl[1]+0.25*(bxr[0]+bxl[0]))*c2gamma; 
  incr[1] = (0.9682458365518543*psr[2]-0.9682458365518543*psl[2]-0.75*(psr[1]+psl[1])+0.4330127018922193*psr[0]-0.4330127018922193*psl[0])*c*gamma+((-0.9682458365518543*(bxr[2]+bxl[2]))+0.75*bxr[1]-0.75*bxl[1]-0.4330127018922193*(bxr[0]+bxl[0]))*c2gamma; 
  incr[2] = ((-1.25*psr[2])+1.25*psl[2]+0.9682458365518543*(psr[1]+psl[1])-0.5590169943749475*psr[0]+0.5590169943749475*psl[0])*c*gamma+(1.25*(bxr[2]+bxl[2])-0.9682458365518543*bxr[1]+0.9682458365518543*bxl[1]+0.5590169943749475*(bxr[0]+bxl[0]))*c2gamma; 

  outPsr[0] += incr[0]*dxr1; 
  outPsr[1] += incr[1]*dxr1; 
  outPsr[2] += incr[2]*dxr1; 

  outPsl[0] += -1.0*incr[0]*dxl1; 
  outPsl[1] += incr[1]*dxl1; 
  outPsl[2] += -1.0*incr[2]*dxl1; 

  return fmax(c, tau); 
} 
