#include <gkyl_maxwell_kernels.h> 
GKYL_CU_DH double maxwell_surfy_2x_ser_p3(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{ 
  const double c = meq->c, chi = meq->chi, gamma = meq->gamma; 
  const double c2 = c*c; 
  const double c2chi = c2*chi, c2gamma = c2*gamma; 
 
  const double dx1 = 2.0/dx[1]; 
  const double *exl = &ql[0]; 
  const double *eyl = &ql[12]; 
  const double *ezl = &ql[24]; 
  const double *bxl = &ql[36]; 
  const double *byl = &ql[48]; 
  const double *bzl = &ql[60]; 
  const double *phl = &ql[72]; 
  const double *psl = &ql[84]; 
 
  const double *exc = &qc[0]; 
  const double *eyc = &qc[12]; 
  const double *ezc = &qc[24]; 
  const double *bxc = &qc[36]; 
  const double *byc = &qc[48]; 
  const double *bzc = &qc[60]; 
  const double *phc = &qc[72]; 
  const double *psc = &qc[84]; 
 
  const double *exr = &qr[0]; 
  const double *eyr = &qr[12]; 
  const double *ezr = &qr[24]; 
  const double *bxr = &qr[36]; 
  const double *byr = &qr[48]; 
  const double *bzr = &qr[60]; 
  const double *phr = &qr[72]; 
  const double *psr = &qr[84]; 
 
  double *outEx = &out[0]; 
  double *outEy = &out[12]; 
  double *outEz = &out[24]; 
  double *outBx = &out[36]; 
  double *outBy = &out[48]; 
  double *outBz = &out[60]; 
  double *outPh = &out[72]; 
  double *outPs = &out[84]; 
 
  double incr_l[12]; 
 
  double incr_r[12]; 
 
  incr_l[0] = ((-0.6614378277661477*bzl[9])+0.6614378277661477*bzc[9]-0.5590169943749475*(bzl[5]+bzc[5])-0.4330127018922193*bzl[2]+0.4330127018922193*bzc[2]-0.25*(bzl[0]+bzc[0]))*c2; 
  incr_l[1] = ((-0.6614378277661477*bzl[11])+0.6614378277661477*bzc[11]-0.5590169943749476*(bzl[7]+bzc[7])-0.4330127018922193*bzl[3]+0.4330127018922193*bzc[3]-0.25*(bzl[1]+bzc[1]))*c2; 
  incr_l[2] = (1.14564392373896*bzl[9]-1.14564392373896*bzc[9]+0.9682458365518543*(bzl[5]+bzc[5])+0.75*bzl[2]-0.75*bzc[2]+0.4330127018922193*(bzl[0]+bzc[0]))*c2; 
  incr_l[3] = (1.14564392373896*bzl[11]-1.14564392373896*bzc[11]+0.9682458365518543*(bzl[7]+bzc[7])+0.75*bzl[3]-0.75*bzc[3]+0.4330127018922193*(bzl[1]+bzc[1]))*c2; 
  incr_l[4] = ((-0.4330127018922194*bzl[6])+0.4330127018922194*bzc[6]-0.25*(bzl[4]+bzc[4]))*c2; 
  incr_l[5] = ((-1.479019945774904*bzl[9])+1.479019945774904*bzc[9]-1.25*(bzl[5]+bzc[5])-0.9682458365518543*bzl[2]+0.9682458365518543*bzc[2]-0.5590169943749475*(bzl[0]+bzc[0]))*c2; 
  incr_l[6] = (0.75*bzl[6]-0.75*bzc[6]+0.4330127018922194*(bzl[4]+bzc[4]))*c2; 
  incr_l[7] = ((-1.479019945774904*bzl[11])+1.479019945774904*bzc[11]-1.25*(bzl[7]+bzc[7])-0.9682458365518543*bzl[3]+0.9682458365518543*bzc[3]-0.5590169943749476*(bzl[1]+bzc[1]))*c2; 
  incr_l[8] = ((-0.4330127018922193*bzl[10])+0.4330127018922193*bzc[10]-0.25*(bzl[8]+bzc[8]))*c2; 
  incr_l[9] = (1.75*bzl[9]-1.75*bzc[9]+1.479019945774904*(bzl[5]+bzc[5])+1.14564392373896*bzl[2]-1.14564392373896*bzc[2]+0.6614378277661477*(bzl[0]+bzc[0]))*c2; 
  incr_l[10] = (0.75*bzl[10]-0.75*bzc[10]+0.4330127018922193*(bzl[8]+bzc[8]))*c2; 
  incr_l[11] = (1.75*bzl[11]-1.75*bzc[11]+1.479019945774904*(bzl[7]+bzc[7])+1.14564392373896*bzl[3]-1.14564392373896*bzc[3]+0.6614378277661477*(bzl[1]+bzc[1]))*c2; 

  incr_r[0] = ((-0.6614378277661477*bzr[9])+0.6614378277661477*bzc[9]+0.5590169943749475*(bzr[5]+bzc[5])-0.4330127018922193*bzr[2]+0.4330127018922193*bzc[2]+0.25*(bzr[0]+bzc[0]))*c2; 
  incr_r[1] = ((-0.6614378277661477*bzr[11])+0.6614378277661477*bzc[11]+0.5590169943749476*(bzr[7]+bzc[7])-0.4330127018922193*bzr[3]+0.4330127018922193*bzc[3]+0.25*(bzr[1]+bzc[1]))*c2; 
  incr_r[2] = ((-1.14564392373896*bzr[9])+1.14564392373896*bzc[9]+0.9682458365518543*(bzr[5]+bzc[5])-0.75*bzr[2]+0.75*bzc[2]+0.4330127018922193*(bzr[0]+bzc[0]))*c2; 
  incr_r[3] = ((-1.14564392373896*bzr[11])+1.14564392373896*bzc[11]+0.9682458365518543*(bzr[7]+bzc[7])-0.75*bzr[3]+0.75*bzc[3]+0.4330127018922193*(bzr[1]+bzc[1]))*c2; 
  incr_r[4] = ((-0.4330127018922194*bzr[6])+0.4330127018922194*bzc[6]+0.25*(bzr[4]+bzc[4]))*c2; 
  incr_r[5] = ((-1.479019945774904*bzr[9])+1.479019945774904*bzc[9]+1.25*(bzr[5]+bzc[5])-0.9682458365518543*bzr[2]+0.9682458365518543*bzc[2]+0.5590169943749475*(bzr[0]+bzc[0]))*c2; 
  incr_r[6] = ((-0.75*bzr[6])+0.75*bzc[6]+0.4330127018922194*(bzr[4]+bzc[4]))*c2; 
  incr_r[7] = ((-1.479019945774904*bzr[11])+1.479019945774904*bzc[11]+1.25*(bzr[7]+bzc[7])-0.9682458365518543*bzr[3]+0.9682458365518543*bzc[3]+0.5590169943749476*(bzr[1]+bzc[1]))*c2; 
  incr_r[8] = ((-0.4330127018922193*bzr[10])+0.4330127018922193*bzc[10]+0.25*(bzr[8]+bzc[8]))*c2; 
  incr_r[9] = ((-1.75*bzr[9])+1.75*bzc[9]+1.479019945774904*(bzr[5]+bzc[5])-1.14564392373896*bzr[2]+1.14564392373896*bzc[2]+0.6614378277661477*(bzr[0]+bzc[0]))*c2; 
  incr_r[10] = ((-0.75*bzr[10])+0.75*bzc[10]+0.4330127018922193*(bzr[8]+bzc[8]))*c2; 
  incr_r[11] = ((-1.75*bzr[11])+1.75*bzc[11]+1.479019945774904*(bzr[7]+bzc[7])-1.14564392373896*bzr[3]+1.14564392373896*bzc[3]+0.6614378277661477*(bzr[1]+bzc[1]))*c2; 

  outEx[0] += (incr_r[0]+incr_l[0])*dx1; 
  outEx[1] += (incr_r[1]+incr_l[1])*dx1; 
  outEx[2] += (incr_r[2]+incr_l[2])*dx1; 
  outEx[3] += (incr_r[3]+incr_l[3])*dx1; 
  outEx[4] += (incr_r[4]+incr_l[4])*dx1; 
  outEx[5] += (incr_r[5]+incr_l[5])*dx1; 
  outEx[6] += (incr_r[6]+incr_l[6])*dx1; 
  outEx[7] += (incr_r[7]+incr_l[7])*dx1; 
  outEx[8] += (incr_r[8]+incr_l[8])*dx1; 
  outEx[9] += (incr_r[9]+incr_l[9])*dx1; 
  outEx[10] += (incr_r[10]+incr_l[10])*dx1; 
  outEx[11] += (incr_r[11]+incr_l[11])*dx1; 

  incr_l[0] = (0.6614378277661477*phl[9]-0.6614378277661477*phc[9]+0.5590169943749475*(phl[5]+phc[5])+0.4330127018922193*phl[2]-0.4330127018922193*phc[2]+0.25*(phl[0]+phc[0]))*c2chi; 
  incr_l[1] = (0.6614378277661477*phl[11]-0.6614378277661477*phc[11]+0.5590169943749476*(phl[7]+phc[7])+0.4330127018922193*phl[3]-0.4330127018922193*phc[3]+0.25*(phl[1]+phc[1]))*c2chi; 
  incr_l[2] = ((-1.14564392373896*phl[9])+1.14564392373896*phc[9]-0.9682458365518543*(phl[5]+phc[5])-0.75*phl[2]+0.75*phc[2]-0.4330127018922193*(phl[0]+phc[0]))*c2chi; 
  incr_l[3] = ((-1.14564392373896*phl[11])+1.14564392373896*phc[11]-0.9682458365518543*(phl[7]+phc[7])-0.75*phl[3]+0.75*phc[3]-0.4330127018922193*(phl[1]+phc[1]))*c2chi; 
  incr_l[4] = (0.4330127018922194*phl[6]-0.4330127018922194*phc[6]+0.25*(phl[4]+phc[4]))*c2chi; 
  incr_l[5] = (1.479019945774904*phl[9]-1.479019945774904*phc[9]+1.25*(phl[5]+phc[5])+0.9682458365518543*phl[2]-0.9682458365518543*phc[2]+0.5590169943749475*(phl[0]+phc[0]))*c2chi; 
  incr_l[6] = ((-0.75*phl[6])+0.75*phc[6]-0.4330127018922194*(phl[4]+phc[4]))*c2chi; 
  incr_l[7] = (1.479019945774904*phl[11]-1.479019945774904*phc[11]+1.25*(phl[7]+phc[7])+0.9682458365518543*phl[3]-0.9682458365518543*phc[3]+0.5590169943749476*(phl[1]+phc[1]))*c2chi; 
  incr_l[8] = (0.4330127018922193*phl[10]-0.4330127018922193*phc[10]+0.25*(phl[8]+phc[8]))*c2chi; 
  incr_l[9] = ((-1.75*phl[9])+1.75*phc[9]-1.479019945774904*(phl[5]+phc[5])-1.14564392373896*phl[2]+1.14564392373896*phc[2]-0.6614378277661477*(phl[0]+phc[0]))*c2chi; 
  incr_l[10] = ((-0.75*phl[10])+0.75*phc[10]-0.4330127018922193*(phl[8]+phc[8]))*c2chi; 
  incr_l[11] = ((-1.75*phl[11])+1.75*phc[11]-1.479019945774904*(phl[7]+phc[7])-1.14564392373896*phl[3]+1.14564392373896*phc[3]-0.6614378277661477*(phl[1]+phc[1]))*c2chi; 

  incr_r[0] = (0.6614378277661477*phr[9]-0.6614378277661477*phc[9]-0.5590169943749475*(phr[5]+phc[5])+0.4330127018922193*phr[2]-0.4330127018922193*phc[2]-0.25*(phr[0]+phc[0]))*c2chi; 
  incr_r[1] = (0.6614378277661477*phr[11]-0.6614378277661477*phc[11]-0.5590169943749476*(phr[7]+phc[7])+0.4330127018922193*phr[3]-0.4330127018922193*phc[3]-0.25*(phr[1]+phc[1]))*c2chi; 
  incr_r[2] = (1.14564392373896*phr[9]-1.14564392373896*phc[9]-0.9682458365518543*(phr[5]+phc[5])+0.75*phr[2]-0.75*phc[2]-0.4330127018922193*(phr[0]+phc[0]))*c2chi; 
  incr_r[3] = (1.14564392373896*phr[11]-1.14564392373896*phc[11]-0.9682458365518543*(phr[7]+phc[7])+0.75*phr[3]-0.75*phc[3]-0.4330127018922193*(phr[1]+phc[1]))*c2chi; 
  incr_r[4] = (0.4330127018922194*phr[6]-0.4330127018922194*phc[6]-0.25*(phr[4]+phc[4]))*c2chi; 
  incr_r[5] = (1.479019945774904*phr[9]-1.479019945774904*phc[9]-1.25*(phr[5]+phc[5])+0.9682458365518543*phr[2]-0.9682458365518543*phc[2]-0.5590169943749475*(phr[0]+phc[0]))*c2chi; 
  incr_r[6] = (0.75*phr[6]-0.75*phc[6]-0.4330127018922194*(phr[4]+phc[4]))*c2chi; 
  incr_r[7] = (1.479019945774904*phr[11]-1.479019945774904*phc[11]-1.25*(phr[7]+phc[7])+0.9682458365518543*phr[3]-0.9682458365518543*phc[3]-0.5590169943749476*(phr[1]+phc[1]))*c2chi; 
  incr_r[8] = (0.4330127018922193*phr[10]-0.4330127018922193*phc[10]-0.25*(phr[8]+phc[8]))*c2chi; 
  incr_r[9] = (1.75*phr[9]-1.75*phc[9]-1.479019945774904*(phr[5]+phc[5])+1.14564392373896*phr[2]-1.14564392373896*phc[2]-0.6614378277661477*(phr[0]+phc[0]))*c2chi; 
  incr_r[10] = (0.75*phr[10]-0.75*phc[10]-0.4330127018922193*(phr[8]+phc[8]))*c2chi; 
  incr_r[11] = (1.75*phr[11]-1.75*phc[11]-1.479019945774904*(phr[7]+phc[7])+1.14564392373896*phr[3]-1.14564392373896*phc[3]-0.6614378277661477*(phr[1]+phc[1]))*c2chi; 

  outEy[0] += (incr_r[0]+incr_l[0])*dx1; 
  outEy[1] += (incr_r[1]+incr_l[1])*dx1; 
  outEy[2] += (incr_r[2]+incr_l[2])*dx1; 
  outEy[3] += (incr_r[3]+incr_l[3])*dx1; 
  outEy[4] += (incr_r[4]+incr_l[4])*dx1; 
  outEy[5] += (incr_r[5]+incr_l[5])*dx1; 
  outEy[6] += (incr_r[6]+incr_l[6])*dx1; 
  outEy[7] += (incr_r[7]+incr_l[7])*dx1; 
  outEy[8] += (incr_r[8]+incr_l[8])*dx1; 
  outEy[9] += (incr_r[9]+incr_l[9])*dx1; 
  outEy[10] += (incr_r[10]+incr_l[10])*dx1; 
  outEy[11] += (incr_r[11]+incr_l[11])*dx1; 

  incr_l[0] = (0.6614378277661477*bxl[9]-0.6614378277661477*bxc[9]+0.5590169943749475*(bxl[5]+bxc[5])+0.4330127018922193*bxl[2]-0.4330127018922193*bxc[2]+0.25*(bxl[0]+bxc[0]))*c2; 
  incr_l[1] = (0.6614378277661477*bxl[11]-0.6614378277661477*bxc[11]+0.5590169943749476*(bxl[7]+bxc[7])+0.4330127018922193*bxl[3]-0.4330127018922193*bxc[3]+0.25*(bxl[1]+bxc[1]))*c2; 
  incr_l[2] = ((-1.14564392373896*bxl[9])+1.14564392373896*bxc[9]-0.9682458365518543*(bxl[5]+bxc[5])-0.75*bxl[2]+0.75*bxc[2]-0.4330127018922193*(bxl[0]+bxc[0]))*c2; 
  incr_l[3] = ((-1.14564392373896*bxl[11])+1.14564392373896*bxc[11]-0.9682458365518543*(bxl[7]+bxc[7])-0.75*bxl[3]+0.75*bxc[3]-0.4330127018922193*(bxl[1]+bxc[1]))*c2; 
  incr_l[4] = (0.4330127018922194*bxl[6]-0.4330127018922194*bxc[6]+0.25*(bxl[4]+bxc[4]))*c2; 
  incr_l[5] = (1.479019945774904*bxl[9]-1.479019945774904*bxc[9]+1.25*(bxl[5]+bxc[5])+0.9682458365518543*bxl[2]-0.9682458365518543*bxc[2]+0.5590169943749475*(bxl[0]+bxc[0]))*c2; 
  incr_l[6] = ((-0.75*bxl[6])+0.75*bxc[6]-0.4330127018922194*(bxl[4]+bxc[4]))*c2; 
  incr_l[7] = (1.479019945774904*bxl[11]-1.479019945774904*bxc[11]+1.25*(bxl[7]+bxc[7])+0.9682458365518543*bxl[3]-0.9682458365518543*bxc[3]+0.5590169943749476*(bxl[1]+bxc[1]))*c2; 
  incr_l[8] = (0.4330127018922193*bxl[10]-0.4330127018922193*bxc[10]+0.25*(bxl[8]+bxc[8]))*c2; 
  incr_l[9] = ((-1.75*bxl[9])+1.75*bxc[9]-1.479019945774904*(bxl[5]+bxc[5])-1.14564392373896*bxl[2]+1.14564392373896*bxc[2]-0.6614378277661477*(bxl[0]+bxc[0]))*c2; 
  incr_l[10] = ((-0.75*bxl[10])+0.75*bxc[10]-0.4330127018922193*(bxl[8]+bxc[8]))*c2; 
  incr_l[11] = ((-1.75*bxl[11])+1.75*bxc[11]-1.479019945774904*(bxl[7]+bxc[7])-1.14564392373896*bxl[3]+1.14564392373896*bxc[3]-0.6614378277661477*(bxl[1]+bxc[1]))*c2; 

  incr_r[0] = (0.6614378277661477*bxr[9]-0.6614378277661477*bxc[9]-0.5590169943749475*(bxr[5]+bxc[5])+0.4330127018922193*bxr[2]-0.4330127018922193*bxc[2]-0.25*(bxr[0]+bxc[0]))*c2; 
  incr_r[1] = (0.6614378277661477*bxr[11]-0.6614378277661477*bxc[11]-0.5590169943749476*(bxr[7]+bxc[7])+0.4330127018922193*bxr[3]-0.4330127018922193*bxc[3]-0.25*(bxr[1]+bxc[1]))*c2; 
  incr_r[2] = (1.14564392373896*bxr[9]-1.14564392373896*bxc[9]-0.9682458365518543*(bxr[5]+bxc[5])+0.75*bxr[2]-0.75*bxc[2]-0.4330127018922193*(bxr[0]+bxc[0]))*c2; 
  incr_r[3] = (1.14564392373896*bxr[11]-1.14564392373896*bxc[11]-0.9682458365518543*(bxr[7]+bxc[7])+0.75*bxr[3]-0.75*bxc[3]-0.4330127018922193*(bxr[1]+bxc[1]))*c2; 
  incr_r[4] = (0.4330127018922194*bxr[6]-0.4330127018922194*bxc[6]-0.25*(bxr[4]+bxc[4]))*c2; 
  incr_r[5] = (1.479019945774904*bxr[9]-1.479019945774904*bxc[9]-1.25*(bxr[5]+bxc[5])+0.9682458365518543*bxr[2]-0.9682458365518543*bxc[2]-0.5590169943749475*(bxr[0]+bxc[0]))*c2; 
  incr_r[6] = (0.75*bxr[6]-0.75*bxc[6]-0.4330127018922194*(bxr[4]+bxc[4]))*c2; 
  incr_r[7] = (1.479019945774904*bxr[11]-1.479019945774904*bxc[11]-1.25*(bxr[7]+bxc[7])+0.9682458365518543*bxr[3]-0.9682458365518543*bxc[3]-0.5590169943749476*(bxr[1]+bxc[1]))*c2; 
  incr_r[8] = (0.4330127018922193*bxr[10]-0.4330127018922193*bxc[10]-0.25*(bxr[8]+bxc[8]))*c2; 
  incr_r[9] = (1.75*bxr[9]-1.75*bxc[9]-1.479019945774904*(bxr[5]+bxc[5])+1.14564392373896*bxr[2]-1.14564392373896*bxc[2]-0.6614378277661477*(bxr[0]+bxc[0]))*c2; 
  incr_r[10] = (0.75*bxr[10]-0.75*bxc[10]-0.4330127018922193*(bxr[8]+bxc[8]))*c2; 
  incr_r[11] = (1.75*bxr[11]-1.75*bxc[11]-1.479019945774904*(bxr[7]+bxc[7])+1.14564392373896*bxr[3]-1.14564392373896*bxc[3]-0.6614378277661477*(bxr[1]+bxc[1]))*c2; 

  outEz[0] += (incr_r[0]+incr_l[0])*dx1; 
  outEz[1] += (incr_r[1]+incr_l[1])*dx1; 
  outEz[2] += (incr_r[2]+incr_l[2])*dx1; 
  outEz[3] += (incr_r[3]+incr_l[3])*dx1; 
  outEz[4] += (incr_r[4]+incr_l[4])*dx1; 
  outEz[5] += (incr_r[5]+incr_l[5])*dx1; 
  outEz[6] += (incr_r[6]+incr_l[6])*dx1; 
  outEz[7] += (incr_r[7]+incr_l[7])*dx1; 
  outEz[8] += (incr_r[8]+incr_l[8])*dx1; 
  outEz[9] += (incr_r[9]+incr_l[9])*dx1; 
  outEz[10] += (incr_r[10]+incr_l[10])*dx1; 
  outEz[11] += (incr_r[11]+incr_l[11])*dx1; 

  incr_l[0] = 0.6614378277661477*ezl[9]-0.6614378277661477*ezc[9]+0.5590169943749475*(ezl[5]+ezc[5])+0.4330127018922193*ezl[2]-0.4330127018922193*ezc[2]+0.25*(ezl[0]+ezc[0]); 
  incr_l[1] = 0.6614378277661477*ezl[11]-0.6614378277661477*ezc[11]+0.5590169943749476*(ezl[7]+ezc[7])+0.4330127018922193*ezl[3]-0.4330127018922193*ezc[3]+0.25*(ezl[1]+ezc[1]); 
  incr_l[2] = (-1.14564392373896*ezl[9])+1.14564392373896*ezc[9]-0.9682458365518543*(ezl[5]+ezc[5])-0.75*ezl[2]+0.75*ezc[2]-0.4330127018922193*(ezl[0]+ezc[0]); 
  incr_l[3] = (-1.14564392373896*ezl[11])+1.14564392373896*ezc[11]-0.9682458365518543*(ezl[7]+ezc[7])-0.75*ezl[3]+0.75*ezc[3]-0.4330127018922193*(ezl[1]+ezc[1]); 
  incr_l[4] = 0.4330127018922194*ezl[6]-0.4330127018922194*ezc[6]+0.25*(ezl[4]+ezc[4]); 
  incr_l[5] = 1.479019945774904*ezl[9]-1.479019945774904*ezc[9]+1.25*(ezl[5]+ezc[5])+0.9682458365518543*ezl[2]-0.9682458365518543*ezc[2]+0.5590169943749475*(ezl[0]+ezc[0]); 
  incr_l[6] = (-0.75*ezl[6])+0.75*ezc[6]-0.4330127018922194*(ezl[4]+ezc[4]); 
  incr_l[7] = 1.479019945774904*ezl[11]-1.479019945774904*ezc[11]+1.25*(ezl[7]+ezc[7])+0.9682458365518543*ezl[3]-0.9682458365518543*ezc[3]+0.5590169943749476*(ezl[1]+ezc[1]); 
  incr_l[8] = 0.4330127018922193*ezl[10]-0.4330127018922193*ezc[10]+0.25*(ezl[8]+ezc[8]); 
  incr_l[9] = (-1.75*ezl[9])+1.75*ezc[9]-1.479019945774904*(ezl[5]+ezc[5])-1.14564392373896*ezl[2]+1.14564392373896*ezc[2]-0.6614378277661477*(ezl[0]+ezc[0]); 
  incr_l[10] = (-0.75*ezl[10])+0.75*ezc[10]-0.4330127018922193*(ezl[8]+ezc[8]); 
  incr_l[11] = (-1.75*ezl[11])+1.75*ezc[11]-1.479019945774904*(ezl[7]+ezc[7])-1.14564392373896*ezl[3]+1.14564392373896*ezc[3]-0.6614378277661477*(ezl[1]+ezc[1]); 

  incr_r[0] = 0.6614378277661477*ezr[9]-0.6614378277661477*ezc[9]-0.5590169943749475*(ezr[5]+ezc[5])+0.4330127018922193*ezr[2]-0.4330127018922193*ezc[2]-0.25*(ezr[0]+ezc[0]); 
  incr_r[1] = 0.6614378277661477*ezr[11]-0.6614378277661477*ezc[11]-0.5590169943749476*(ezr[7]+ezc[7])+0.4330127018922193*ezr[3]-0.4330127018922193*ezc[3]-0.25*(ezr[1]+ezc[1]); 
  incr_r[2] = 1.14564392373896*ezr[9]-1.14564392373896*ezc[9]-0.9682458365518543*(ezr[5]+ezc[5])+0.75*ezr[2]-0.75*ezc[2]-0.4330127018922193*(ezr[0]+ezc[0]); 
  incr_r[3] = 1.14564392373896*ezr[11]-1.14564392373896*ezc[11]-0.9682458365518543*(ezr[7]+ezc[7])+0.75*ezr[3]-0.75*ezc[3]-0.4330127018922193*(ezr[1]+ezc[1]); 
  incr_r[4] = 0.4330127018922194*ezr[6]-0.4330127018922194*ezc[6]-0.25*(ezr[4]+ezc[4]); 
  incr_r[5] = 1.479019945774904*ezr[9]-1.479019945774904*ezc[9]-1.25*(ezr[5]+ezc[5])+0.9682458365518543*ezr[2]-0.9682458365518543*ezc[2]-0.5590169943749475*(ezr[0]+ezc[0]); 
  incr_r[6] = 0.75*ezr[6]-0.75*ezc[6]-0.4330127018922194*(ezr[4]+ezc[4]); 
  incr_r[7] = 1.479019945774904*ezr[11]-1.479019945774904*ezc[11]-1.25*(ezr[7]+ezc[7])+0.9682458365518543*ezr[3]-0.9682458365518543*ezc[3]-0.5590169943749476*(ezr[1]+ezc[1]); 
  incr_r[8] = 0.4330127018922193*ezr[10]-0.4330127018922193*ezc[10]-0.25*(ezr[8]+ezc[8]); 
  incr_r[9] = 1.75*ezr[9]-1.75*ezc[9]-1.479019945774904*(ezr[5]+ezc[5])+1.14564392373896*ezr[2]-1.14564392373896*ezc[2]-0.6614378277661477*(ezr[0]+ezc[0]); 
  incr_r[10] = 0.75*ezr[10]-0.75*ezc[10]-0.4330127018922193*(ezr[8]+ezc[8]); 
  incr_r[11] = 1.75*ezr[11]-1.75*ezc[11]-1.479019945774904*(ezr[7]+ezc[7])+1.14564392373896*ezr[3]-1.14564392373896*ezc[3]-0.6614378277661477*(ezr[1]+ezc[1]); 

  outBx[0] += (incr_r[0]+incr_l[0])*dx1; 
  outBx[1] += (incr_r[1]+incr_l[1])*dx1; 
  outBx[2] += (incr_r[2]+incr_l[2])*dx1; 
  outBx[3] += (incr_r[3]+incr_l[3])*dx1; 
  outBx[4] += (incr_r[4]+incr_l[4])*dx1; 
  outBx[5] += (incr_r[5]+incr_l[5])*dx1; 
  outBx[6] += (incr_r[6]+incr_l[6])*dx1; 
  outBx[7] += (incr_r[7]+incr_l[7])*dx1; 
  outBx[8] += (incr_r[8]+incr_l[8])*dx1; 
  outBx[9] += (incr_r[9]+incr_l[9])*dx1; 
  outBx[10] += (incr_r[10]+incr_l[10])*dx1; 
  outBx[11] += (incr_r[11]+incr_l[11])*dx1; 

  incr_l[0] = (0.6614378277661477*psl[9]-0.6614378277661477*psc[9]+0.5590169943749475*(psl[5]+psc[5])+0.4330127018922193*psl[2]-0.4330127018922193*psc[2]+0.25*(psl[0]+psc[0]))*gamma; 
  incr_l[1] = (0.6614378277661477*psl[11]-0.6614378277661477*psc[11]+0.5590169943749476*(psl[7]+psc[7])+0.4330127018922193*psl[3]-0.4330127018922193*psc[3]+0.25*(psl[1]+psc[1]))*gamma; 
  incr_l[2] = ((-1.14564392373896*psl[9])+1.14564392373896*psc[9]-0.9682458365518543*(psl[5]+psc[5])-0.75*psl[2]+0.75*psc[2]-0.4330127018922193*(psl[0]+psc[0]))*gamma; 
  incr_l[3] = ((-1.14564392373896*psl[11])+1.14564392373896*psc[11]-0.9682458365518543*(psl[7]+psc[7])-0.75*psl[3]+0.75*psc[3]-0.4330127018922193*(psl[1]+psc[1]))*gamma; 
  incr_l[4] = (0.4330127018922194*psl[6]-0.4330127018922194*psc[6]+0.25*(psl[4]+psc[4]))*gamma; 
  incr_l[5] = (1.479019945774904*psl[9]-1.479019945774904*psc[9]+1.25*(psl[5]+psc[5])+0.9682458365518543*psl[2]-0.9682458365518543*psc[2]+0.5590169943749475*(psl[0]+psc[0]))*gamma; 
  incr_l[6] = ((-0.75*psl[6])+0.75*psc[6]-0.4330127018922194*(psl[4]+psc[4]))*gamma; 
  incr_l[7] = (1.479019945774904*psl[11]-1.479019945774904*psc[11]+1.25*(psl[7]+psc[7])+0.9682458365518543*psl[3]-0.9682458365518543*psc[3]+0.5590169943749476*(psl[1]+psc[1]))*gamma; 
  incr_l[8] = (0.4330127018922193*psl[10]-0.4330127018922193*psc[10]+0.25*(psl[8]+psc[8]))*gamma; 
  incr_l[9] = ((-1.75*psl[9])+1.75*psc[9]-1.479019945774904*(psl[5]+psc[5])-1.14564392373896*psl[2]+1.14564392373896*psc[2]-0.6614378277661477*(psl[0]+psc[0]))*gamma; 
  incr_l[10] = ((-0.75*psl[10])+0.75*psc[10]-0.4330127018922193*(psl[8]+psc[8]))*gamma; 
  incr_l[11] = ((-1.75*psl[11])+1.75*psc[11]-1.479019945774904*(psl[7]+psc[7])-1.14564392373896*psl[3]+1.14564392373896*psc[3]-0.6614378277661477*(psl[1]+psc[1]))*gamma; 

  incr_r[0] = (0.6614378277661477*psr[9]-0.6614378277661477*psc[9]-0.5590169943749475*(psr[5]+psc[5])+0.4330127018922193*psr[2]-0.4330127018922193*psc[2]-0.25*(psr[0]+psc[0]))*gamma; 
  incr_r[1] = (0.6614378277661477*psr[11]-0.6614378277661477*psc[11]-0.5590169943749476*(psr[7]+psc[7])+0.4330127018922193*psr[3]-0.4330127018922193*psc[3]-0.25*(psr[1]+psc[1]))*gamma; 
  incr_r[2] = (1.14564392373896*psr[9]-1.14564392373896*psc[9]-0.9682458365518543*(psr[5]+psc[5])+0.75*psr[2]-0.75*psc[2]-0.4330127018922193*(psr[0]+psc[0]))*gamma; 
  incr_r[3] = (1.14564392373896*psr[11]-1.14564392373896*psc[11]-0.9682458365518543*(psr[7]+psc[7])+0.75*psr[3]-0.75*psc[3]-0.4330127018922193*(psr[1]+psc[1]))*gamma; 
  incr_r[4] = (0.4330127018922194*psr[6]-0.4330127018922194*psc[6]-0.25*(psr[4]+psc[4]))*gamma; 
  incr_r[5] = (1.479019945774904*psr[9]-1.479019945774904*psc[9]-1.25*(psr[5]+psc[5])+0.9682458365518543*psr[2]-0.9682458365518543*psc[2]-0.5590169943749475*(psr[0]+psc[0]))*gamma; 
  incr_r[6] = (0.75*psr[6]-0.75*psc[6]-0.4330127018922194*(psr[4]+psc[4]))*gamma; 
  incr_r[7] = (1.479019945774904*psr[11]-1.479019945774904*psc[11]-1.25*(psr[7]+psc[7])+0.9682458365518543*psr[3]-0.9682458365518543*psc[3]-0.5590169943749476*(psr[1]+psc[1]))*gamma; 
  incr_r[8] = (0.4330127018922193*psr[10]-0.4330127018922193*psc[10]-0.25*(psr[8]+psc[8]))*gamma; 
  incr_r[9] = (1.75*psr[9]-1.75*psc[9]-1.479019945774904*(psr[5]+psc[5])+1.14564392373896*psr[2]-1.14564392373896*psc[2]-0.6614378277661477*(psr[0]+psc[0]))*gamma; 
  incr_r[10] = (0.75*psr[10]-0.75*psc[10]-0.4330127018922193*(psr[8]+psc[8]))*gamma; 
  incr_r[11] = (1.75*psr[11]-1.75*psc[11]-1.479019945774904*(psr[7]+psc[7])+1.14564392373896*psr[3]-1.14564392373896*psc[3]-0.6614378277661477*(psr[1]+psc[1]))*gamma; 

  outBy[0] += (incr_r[0]+incr_l[0])*dx1; 
  outBy[1] += (incr_r[1]+incr_l[1])*dx1; 
  outBy[2] += (incr_r[2]+incr_l[2])*dx1; 
  outBy[3] += (incr_r[3]+incr_l[3])*dx1; 
  outBy[4] += (incr_r[4]+incr_l[4])*dx1; 
  outBy[5] += (incr_r[5]+incr_l[5])*dx1; 
  outBy[6] += (incr_r[6]+incr_l[6])*dx1; 
  outBy[7] += (incr_r[7]+incr_l[7])*dx1; 
  outBy[8] += (incr_r[8]+incr_l[8])*dx1; 
  outBy[9] += (incr_r[9]+incr_l[9])*dx1; 
  outBy[10] += (incr_r[10]+incr_l[10])*dx1; 
  outBy[11] += (incr_r[11]+incr_l[11])*dx1; 

  incr_l[0] = (-0.6614378277661477*exl[9])+0.6614378277661477*exc[9]-0.5590169943749475*(exl[5]+exc[5])-0.4330127018922193*exl[2]+0.4330127018922193*exc[2]-0.25*(exl[0]+exc[0]); 
  incr_l[1] = (-0.6614378277661477*exl[11])+0.6614378277661477*exc[11]-0.5590169943749476*(exl[7]+exc[7])-0.4330127018922193*exl[3]+0.4330127018922193*exc[3]-0.25*(exl[1]+exc[1]); 
  incr_l[2] = 1.14564392373896*exl[9]-1.14564392373896*exc[9]+0.9682458365518543*(exl[5]+exc[5])+0.75*exl[2]-0.75*exc[2]+0.4330127018922193*(exl[0]+exc[0]); 
  incr_l[3] = 1.14564392373896*exl[11]-1.14564392373896*exc[11]+0.9682458365518543*(exl[7]+exc[7])+0.75*exl[3]-0.75*exc[3]+0.4330127018922193*(exl[1]+exc[1]); 
  incr_l[4] = (-0.4330127018922194*exl[6])+0.4330127018922194*exc[6]-0.25*(exl[4]+exc[4]); 
  incr_l[5] = (-1.479019945774904*exl[9])+1.479019945774904*exc[9]-1.25*(exl[5]+exc[5])-0.9682458365518543*exl[2]+0.9682458365518543*exc[2]-0.5590169943749475*(exl[0]+exc[0]); 
  incr_l[6] = 0.75*exl[6]-0.75*exc[6]+0.4330127018922194*(exl[4]+exc[4]); 
  incr_l[7] = (-1.479019945774904*exl[11])+1.479019945774904*exc[11]-1.25*(exl[7]+exc[7])-0.9682458365518543*exl[3]+0.9682458365518543*exc[3]-0.5590169943749476*(exl[1]+exc[1]); 
  incr_l[8] = (-0.4330127018922193*exl[10])+0.4330127018922193*exc[10]-0.25*(exl[8]+exc[8]); 
  incr_l[9] = 1.75*exl[9]-1.75*exc[9]+1.479019945774904*(exl[5]+exc[5])+1.14564392373896*exl[2]-1.14564392373896*exc[2]+0.6614378277661477*(exl[0]+exc[0]); 
  incr_l[10] = 0.75*exl[10]-0.75*exc[10]+0.4330127018922193*(exl[8]+exc[8]); 
  incr_l[11] = 1.75*exl[11]-1.75*exc[11]+1.479019945774904*(exl[7]+exc[7])+1.14564392373896*exl[3]-1.14564392373896*exc[3]+0.6614378277661477*(exl[1]+exc[1]); 

  incr_r[0] = (-0.6614378277661477*exr[9])+0.6614378277661477*exc[9]+0.5590169943749475*(exr[5]+exc[5])-0.4330127018922193*exr[2]+0.4330127018922193*exc[2]+0.25*(exr[0]+exc[0]); 
  incr_r[1] = (-0.6614378277661477*exr[11])+0.6614378277661477*exc[11]+0.5590169943749476*(exr[7]+exc[7])-0.4330127018922193*exr[3]+0.4330127018922193*exc[3]+0.25*(exr[1]+exc[1]); 
  incr_r[2] = (-1.14564392373896*exr[9])+1.14564392373896*exc[9]+0.9682458365518543*(exr[5]+exc[5])-0.75*exr[2]+0.75*exc[2]+0.4330127018922193*(exr[0]+exc[0]); 
  incr_r[3] = (-1.14564392373896*exr[11])+1.14564392373896*exc[11]+0.9682458365518543*(exr[7]+exc[7])-0.75*exr[3]+0.75*exc[3]+0.4330127018922193*(exr[1]+exc[1]); 
  incr_r[4] = (-0.4330127018922194*exr[6])+0.4330127018922194*exc[6]+0.25*(exr[4]+exc[4]); 
  incr_r[5] = (-1.479019945774904*exr[9])+1.479019945774904*exc[9]+1.25*(exr[5]+exc[5])-0.9682458365518543*exr[2]+0.9682458365518543*exc[2]+0.5590169943749475*(exr[0]+exc[0]); 
  incr_r[6] = (-0.75*exr[6])+0.75*exc[6]+0.4330127018922194*(exr[4]+exc[4]); 
  incr_r[7] = (-1.479019945774904*exr[11])+1.479019945774904*exc[11]+1.25*(exr[7]+exc[7])-0.9682458365518543*exr[3]+0.9682458365518543*exc[3]+0.5590169943749476*(exr[1]+exc[1]); 
  incr_r[8] = (-0.4330127018922193*exr[10])+0.4330127018922193*exc[10]+0.25*(exr[8]+exc[8]); 
  incr_r[9] = (-1.75*exr[9])+1.75*exc[9]+1.479019945774904*(exr[5]+exc[5])-1.14564392373896*exr[2]+1.14564392373896*exc[2]+0.6614378277661477*(exr[0]+exc[0]); 
  incr_r[10] = (-0.75*exr[10])+0.75*exc[10]+0.4330127018922193*(exr[8]+exc[8]); 
  incr_r[11] = (-1.75*exr[11])+1.75*exc[11]+1.479019945774904*(exr[7]+exc[7])-1.14564392373896*exr[3]+1.14564392373896*exc[3]+0.6614378277661477*(exr[1]+exc[1]); 

  outBz[0] += (incr_r[0]+incr_l[0])*dx1; 
  outBz[1] += (incr_r[1]+incr_l[1])*dx1; 
  outBz[2] += (incr_r[2]+incr_l[2])*dx1; 
  outBz[3] += (incr_r[3]+incr_l[3])*dx1; 
  outBz[4] += (incr_r[4]+incr_l[4])*dx1; 
  outBz[5] += (incr_r[5]+incr_l[5])*dx1; 
  outBz[6] += (incr_r[6]+incr_l[6])*dx1; 
  outBz[7] += (incr_r[7]+incr_l[7])*dx1; 
  outBz[8] += (incr_r[8]+incr_l[8])*dx1; 
  outBz[9] += (incr_r[9]+incr_l[9])*dx1; 
  outBz[10] += (incr_r[10]+incr_l[10])*dx1; 
  outBz[11] += (incr_r[11]+incr_l[11])*dx1; 

  incr_l[0] = (0.6614378277661477*eyl[9]-0.6614378277661477*eyc[9]+0.5590169943749475*(eyl[5]+eyc[5])+0.4330127018922193*eyl[2]-0.4330127018922193*eyc[2]+0.25*(eyl[0]+eyc[0]))*chi; 
  incr_l[1] = (0.6614378277661477*eyl[11]-0.6614378277661477*eyc[11]+0.5590169943749476*(eyl[7]+eyc[7])+0.4330127018922193*eyl[3]-0.4330127018922193*eyc[3]+0.25*(eyl[1]+eyc[1]))*chi; 
  incr_l[2] = ((-1.14564392373896*eyl[9])+1.14564392373896*eyc[9]-0.9682458365518543*(eyl[5]+eyc[5])-0.75*eyl[2]+0.75*eyc[2]-0.4330127018922193*(eyl[0]+eyc[0]))*chi; 
  incr_l[3] = ((-1.14564392373896*eyl[11])+1.14564392373896*eyc[11]-0.9682458365518543*(eyl[7]+eyc[7])-0.75*eyl[3]+0.75*eyc[3]-0.4330127018922193*(eyl[1]+eyc[1]))*chi; 
  incr_l[4] = (0.4330127018922194*eyl[6]-0.4330127018922194*eyc[6]+0.25*(eyl[4]+eyc[4]))*chi; 
  incr_l[5] = (1.479019945774904*eyl[9]-1.479019945774904*eyc[9]+1.25*(eyl[5]+eyc[5])+0.9682458365518543*eyl[2]-0.9682458365518543*eyc[2]+0.5590169943749475*(eyl[0]+eyc[0]))*chi; 
  incr_l[6] = ((-0.75*eyl[6])+0.75*eyc[6]-0.4330127018922194*(eyl[4]+eyc[4]))*chi; 
  incr_l[7] = (1.479019945774904*eyl[11]-1.479019945774904*eyc[11]+1.25*(eyl[7]+eyc[7])+0.9682458365518543*eyl[3]-0.9682458365518543*eyc[3]+0.5590169943749476*(eyl[1]+eyc[1]))*chi; 
  incr_l[8] = (0.4330127018922193*eyl[10]-0.4330127018922193*eyc[10]+0.25*(eyl[8]+eyc[8]))*chi; 
  incr_l[9] = ((-1.75*eyl[9])+1.75*eyc[9]-1.479019945774904*(eyl[5]+eyc[5])-1.14564392373896*eyl[2]+1.14564392373896*eyc[2]-0.6614378277661477*(eyl[0]+eyc[0]))*chi; 
  incr_l[10] = ((-0.75*eyl[10])+0.75*eyc[10]-0.4330127018922193*(eyl[8]+eyc[8]))*chi; 
  incr_l[11] = ((-1.75*eyl[11])+1.75*eyc[11]-1.479019945774904*(eyl[7]+eyc[7])-1.14564392373896*eyl[3]+1.14564392373896*eyc[3]-0.6614378277661477*(eyl[1]+eyc[1]))*chi; 

  incr_r[0] = (0.6614378277661477*eyr[9]-0.6614378277661477*eyc[9]-0.5590169943749475*(eyr[5]+eyc[5])+0.4330127018922193*eyr[2]-0.4330127018922193*eyc[2]-0.25*(eyr[0]+eyc[0]))*chi; 
  incr_r[1] = (0.6614378277661477*eyr[11]-0.6614378277661477*eyc[11]-0.5590169943749476*(eyr[7]+eyc[7])+0.4330127018922193*eyr[3]-0.4330127018922193*eyc[3]-0.25*(eyr[1]+eyc[1]))*chi; 
  incr_r[2] = (1.14564392373896*eyr[9]-1.14564392373896*eyc[9]-0.9682458365518543*(eyr[5]+eyc[5])+0.75*eyr[2]-0.75*eyc[2]-0.4330127018922193*(eyr[0]+eyc[0]))*chi; 
  incr_r[3] = (1.14564392373896*eyr[11]-1.14564392373896*eyc[11]-0.9682458365518543*(eyr[7]+eyc[7])+0.75*eyr[3]-0.75*eyc[3]-0.4330127018922193*(eyr[1]+eyc[1]))*chi; 
  incr_r[4] = (0.4330127018922194*eyr[6]-0.4330127018922194*eyc[6]-0.25*(eyr[4]+eyc[4]))*chi; 
  incr_r[5] = (1.479019945774904*eyr[9]-1.479019945774904*eyc[9]-1.25*(eyr[5]+eyc[5])+0.9682458365518543*eyr[2]-0.9682458365518543*eyc[2]-0.5590169943749475*(eyr[0]+eyc[0]))*chi; 
  incr_r[6] = (0.75*eyr[6]-0.75*eyc[6]-0.4330127018922194*(eyr[4]+eyc[4]))*chi; 
  incr_r[7] = (1.479019945774904*eyr[11]-1.479019945774904*eyc[11]-1.25*(eyr[7]+eyc[7])+0.9682458365518543*eyr[3]-0.9682458365518543*eyc[3]-0.5590169943749476*(eyr[1]+eyc[1]))*chi; 
  incr_r[8] = (0.4330127018922193*eyr[10]-0.4330127018922193*eyc[10]-0.25*(eyr[8]+eyc[8]))*chi; 
  incr_r[9] = (1.75*eyr[9]-1.75*eyc[9]-1.479019945774904*(eyr[5]+eyc[5])+1.14564392373896*eyr[2]-1.14564392373896*eyc[2]-0.6614378277661477*(eyr[0]+eyc[0]))*chi; 
  incr_r[10] = (0.75*eyr[10]-0.75*eyc[10]-0.4330127018922193*(eyr[8]+eyc[8]))*chi; 
  incr_r[11] = (1.75*eyr[11]-1.75*eyc[11]-1.479019945774904*(eyr[7]+eyc[7])+1.14564392373896*eyr[3]-1.14564392373896*eyc[3]-0.6614378277661477*(eyr[1]+eyc[1]))*chi; 

  outPh[0] += (incr_r[0]+incr_l[0])*dx1; 
  outPh[1] += (incr_r[1]+incr_l[1])*dx1; 
  outPh[2] += (incr_r[2]+incr_l[2])*dx1; 
  outPh[3] += (incr_r[3]+incr_l[3])*dx1; 
  outPh[4] += (incr_r[4]+incr_l[4])*dx1; 
  outPh[5] += (incr_r[5]+incr_l[5])*dx1; 
  outPh[6] += (incr_r[6]+incr_l[6])*dx1; 
  outPh[7] += (incr_r[7]+incr_l[7])*dx1; 
  outPh[8] += (incr_r[8]+incr_l[8])*dx1; 
  outPh[9] += (incr_r[9]+incr_l[9])*dx1; 
  outPh[10] += (incr_r[10]+incr_l[10])*dx1; 
  outPh[11] += (incr_r[11]+incr_l[11])*dx1; 

  incr_l[0] = (0.6614378277661477*byl[9]-0.6614378277661477*byc[9]+0.5590169943749475*(byl[5]+byc[5])+0.4330127018922193*byl[2]-0.4330127018922193*byc[2]+0.25*(byl[0]+byc[0]))*c2gamma; 
  incr_l[1] = (0.6614378277661477*byl[11]-0.6614378277661477*byc[11]+0.5590169943749476*(byl[7]+byc[7])+0.4330127018922193*byl[3]-0.4330127018922193*byc[3]+0.25*(byl[1]+byc[1]))*c2gamma; 
  incr_l[2] = ((-1.14564392373896*byl[9])+1.14564392373896*byc[9]-0.9682458365518543*(byl[5]+byc[5])-0.75*byl[2]+0.75*byc[2]-0.4330127018922193*(byl[0]+byc[0]))*c2gamma; 
  incr_l[3] = ((-1.14564392373896*byl[11])+1.14564392373896*byc[11]-0.9682458365518543*(byl[7]+byc[7])-0.75*byl[3]+0.75*byc[3]-0.4330127018922193*(byl[1]+byc[1]))*c2gamma; 
  incr_l[4] = (0.4330127018922194*byl[6]-0.4330127018922194*byc[6]+0.25*(byl[4]+byc[4]))*c2gamma; 
  incr_l[5] = (1.479019945774904*byl[9]-1.479019945774904*byc[9]+1.25*(byl[5]+byc[5])+0.9682458365518543*byl[2]-0.9682458365518543*byc[2]+0.5590169943749475*(byl[0]+byc[0]))*c2gamma; 
  incr_l[6] = ((-0.75*byl[6])+0.75*byc[6]-0.4330127018922194*(byl[4]+byc[4]))*c2gamma; 
  incr_l[7] = (1.479019945774904*byl[11]-1.479019945774904*byc[11]+1.25*(byl[7]+byc[7])+0.9682458365518543*byl[3]-0.9682458365518543*byc[3]+0.5590169943749476*(byl[1]+byc[1]))*c2gamma; 
  incr_l[8] = (0.4330127018922193*byl[10]-0.4330127018922193*byc[10]+0.25*(byl[8]+byc[8]))*c2gamma; 
  incr_l[9] = ((-1.75*byl[9])+1.75*byc[9]-1.479019945774904*(byl[5]+byc[5])-1.14564392373896*byl[2]+1.14564392373896*byc[2]-0.6614378277661477*(byl[0]+byc[0]))*c2gamma; 
  incr_l[10] = ((-0.75*byl[10])+0.75*byc[10]-0.4330127018922193*(byl[8]+byc[8]))*c2gamma; 
  incr_l[11] = ((-1.75*byl[11])+1.75*byc[11]-1.479019945774904*(byl[7]+byc[7])-1.14564392373896*byl[3]+1.14564392373896*byc[3]-0.6614378277661477*(byl[1]+byc[1]))*c2gamma; 

  incr_r[0] = (0.6614378277661477*byr[9]-0.6614378277661477*byc[9]-0.5590169943749475*(byr[5]+byc[5])+0.4330127018922193*byr[2]-0.4330127018922193*byc[2]-0.25*(byr[0]+byc[0]))*c2gamma; 
  incr_r[1] = (0.6614378277661477*byr[11]-0.6614378277661477*byc[11]-0.5590169943749476*(byr[7]+byc[7])+0.4330127018922193*byr[3]-0.4330127018922193*byc[3]-0.25*(byr[1]+byc[1]))*c2gamma; 
  incr_r[2] = (1.14564392373896*byr[9]-1.14564392373896*byc[9]-0.9682458365518543*(byr[5]+byc[5])+0.75*byr[2]-0.75*byc[2]-0.4330127018922193*(byr[0]+byc[0]))*c2gamma; 
  incr_r[3] = (1.14564392373896*byr[11]-1.14564392373896*byc[11]-0.9682458365518543*(byr[7]+byc[7])+0.75*byr[3]-0.75*byc[3]-0.4330127018922193*(byr[1]+byc[1]))*c2gamma; 
  incr_r[4] = (0.4330127018922194*byr[6]-0.4330127018922194*byc[6]-0.25*(byr[4]+byc[4]))*c2gamma; 
  incr_r[5] = (1.479019945774904*byr[9]-1.479019945774904*byc[9]-1.25*(byr[5]+byc[5])+0.9682458365518543*byr[2]-0.9682458365518543*byc[2]-0.5590169943749475*(byr[0]+byc[0]))*c2gamma; 
  incr_r[6] = (0.75*byr[6]-0.75*byc[6]-0.4330127018922194*(byr[4]+byc[4]))*c2gamma; 
  incr_r[7] = (1.479019945774904*byr[11]-1.479019945774904*byc[11]-1.25*(byr[7]+byc[7])+0.9682458365518543*byr[3]-0.9682458365518543*byc[3]-0.5590169943749476*(byr[1]+byc[1]))*c2gamma; 
  incr_r[8] = (0.4330127018922193*byr[10]-0.4330127018922193*byc[10]-0.25*(byr[8]+byc[8]))*c2gamma; 
  incr_r[9] = (1.75*byr[9]-1.75*byc[9]-1.479019945774904*(byr[5]+byc[5])+1.14564392373896*byr[2]-1.14564392373896*byc[2]-0.6614378277661477*(byr[0]+byc[0]))*c2gamma; 
  incr_r[10] = (0.75*byr[10]-0.75*byc[10]-0.4330127018922193*(byr[8]+byc[8]))*c2gamma; 
  incr_r[11] = (1.75*byr[11]-1.75*byc[11]-1.479019945774904*(byr[7]+byc[7])+1.14564392373896*byr[3]-1.14564392373896*byc[3]-0.6614378277661477*(byr[1]+byc[1]))*c2gamma; 

  outPs[0] += (incr_r[0]+incr_l[0])*dx1; 
  outPs[1] += (incr_r[1]+incr_l[1])*dx1; 
  outPs[2] += (incr_r[2]+incr_l[2])*dx1; 
  outPs[3] += (incr_r[3]+incr_l[3])*dx1; 
  outPs[4] += (incr_r[4]+incr_l[4])*dx1; 
  outPs[5] += (incr_r[5]+incr_l[5])*dx1; 
  outPs[6] += (incr_r[6]+incr_l[6])*dx1; 
  outPs[7] += (incr_r[7]+incr_l[7])*dx1; 
  outPs[8] += (incr_r[8]+incr_l[8])*dx1; 
  outPs[9] += (incr_r[9]+incr_l[9])*dx1; 
  outPs[10] += (incr_r[10]+incr_l[10])*dx1; 
  outPs[11] += (incr_r[11]+incr_l[11])*dx1; 

  return 0.;

} 
