#include <gkyl_vlasov_kernels.h> 
void vlasov_stream_surfx_1x3v_ser_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells.
  double rdxl2 = 2.0/dxvl[0]; 
  double rdxr2 = 2.0/dxvr[0]; 

  double incr[16]; 

  if (wr[1]>0) { 
  incr[0] = dxvl[1]*(0.25*fl[5]+0.1443375672974065*fl[2])+(0.8660254037844386*fl[1]+0.5*fl[0])*wl[1]; 
  incr[1] = dxvl[1]*((-0.4330127018922193*fl[5])-0.25*fl[2])+((-1.5*fl[1])-0.8660254037844386*fl[0])*wl[1]; 
  incr[2] = wl[1]*(0.8660254037844386*fl[5]+0.5*fl[2])+dxvl[1]*(0.25*fl[1]+0.1443375672974065*fl[0]); 
  incr[3] = dxvl[1]*(0.25*fl[11]+0.1443375672974065*fl[7])+wl[1]*(0.8660254037844386*fl[6]+0.5*fl[3]); 
  incr[4] = dxvl[1]*(0.25*fl[12]+0.1443375672974065*fl[9])+wl[1]*(0.8660254037844386*fl[8]+0.5*fl[4]); 
  incr[5] = wl[1]*((-1.5*fl[5])-0.8660254037844386*fl[2])+dxvl[1]*((-0.4330127018922193*fl[1])-0.25*fl[0]); 
  incr[6] = dxvl[1]*((-0.4330127018922193*fl[11])-0.25*fl[7])+wl[1]*((-1.5*fl[6])-0.8660254037844386*fl[3]); 
  incr[7] = wl[1]*(0.8660254037844386*fl[11]+0.5*fl[7])+dxvl[1]*(0.25*fl[6]+0.1443375672974065*fl[3]); 
  incr[8] = dxvl[1]*((-0.4330127018922193*fl[12])-0.25*fl[9])+wl[1]*((-1.5*fl[8])-0.8660254037844386*fl[4]); 
  incr[9] = wl[1]*(0.8660254037844386*fl[12]+0.5*fl[9])+dxvl[1]*(0.25*fl[8]+0.1443375672974065*fl[4]); 
  incr[10] = dxvl[1]*(0.25*fl[15]+0.1443375672974065*fl[14])+wl[1]*(0.8660254037844386*fl[13]+0.5*fl[10]); 
  incr[11] = wl[1]*((-1.5*fl[11])-0.8660254037844386*fl[7])+dxvl[1]*((-0.4330127018922193*fl[6])-0.25*fl[3]); 
  incr[12] = wl[1]*((-1.5*fl[12])-0.8660254037844386*fl[9])+dxvl[1]*((-0.4330127018922193*fl[8])-0.25*fl[4]); 
  incr[13] = dxvl[1]*((-0.4330127018922193*fl[15])-0.25*fl[14])+wl[1]*((-1.5*fl[13])-0.8660254037844386*fl[10]); 
  incr[14] = wl[1]*(0.8660254037844386*fl[15]+0.5*fl[14])+dxvl[1]*(0.25*fl[13]+0.1443375672974065*fl[10]); 
  incr[15] = wl[1]*((-1.5*fl[15])-0.8660254037844386*fl[14])+dxvl[1]*((-0.4330127018922193*fl[13])-0.25*fl[10]); 

  outr[0] += incr[0]*rdxr2; 
  outr[1] += incr[1]*rdxr2; 
  outr[2] += incr[2]*rdxr2; 
  outr[3] += incr[3]*rdxr2; 
  outr[4] += incr[4]*rdxr2; 
  outr[5] += incr[5]*rdxr2; 
  outr[6] += incr[6]*rdxr2; 
  outr[7] += incr[7]*rdxr2; 
  outr[8] += incr[8]*rdxr2; 
  outr[9] += incr[9]*rdxr2; 
  outr[10] += incr[10]*rdxr2; 
  outr[11] += incr[11]*rdxr2; 
  outr[12] += incr[12]*rdxr2; 
  outr[13] += incr[13]*rdxr2; 
  outr[14] += incr[14]*rdxr2; 
  outr[15] += incr[15]*rdxr2; 

  outl[0] += -1.0*incr[0]*rdxl2; 
  outl[1] += incr[1]*rdxl2; 
  outl[2] += -1.0*incr[2]*rdxl2; 
  outl[3] += -1.0*incr[3]*rdxl2; 
  outl[4] += -1.0*incr[4]*rdxl2; 
  outl[5] += incr[5]*rdxl2; 
  outl[6] += incr[6]*rdxl2; 
  outl[7] += -1.0*incr[7]*rdxl2; 
  outl[8] += incr[8]*rdxl2; 
  outl[9] += -1.0*incr[9]*rdxl2; 
  outl[10] += -1.0*incr[10]*rdxl2; 
  outl[11] += incr[11]*rdxl2; 
  outl[12] += incr[12]*rdxl2; 
  outl[13] += incr[13]*rdxl2; 
  outl[14] += -1.0*incr[14]*rdxl2; 
  outl[15] += incr[15]*rdxl2; 
  } else { 
  incr[0] = dxvr[1]*(0.1443375672974065*fr[2]-0.25*fr[5])+(0.5*fr[0]-0.8660254037844386*fr[1])*wr[1]; 
  incr[1] = dxvr[1]*(0.4330127018922193*fr[5]-0.25*fr[2])+(1.5*fr[1]-0.8660254037844386*fr[0])*wr[1]; 
  incr[2] = wr[1]*(0.5*fr[2]-0.8660254037844386*fr[5])+dxvr[1]*(0.1443375672974065*fr[0]-0.25*fr[1]); 
  incr[3] = dxvr[1]*(0.1443375672974065*fr[7]-0.25*fr[11])+wr[1]*(0.5*fr[3]-0.8660254037844386*fr[6]); 
  incr[4] = dxvr[1]*(0.1443375672974065*fr[9]-0.25*fr[12])+wr[1]*(0.5*fr[4]-0.8660254037844386*fr[8]); 
  incr[5] = wr[1]*(1.5*fr[5]-0.8660254037844386*fr[2])+dxvr[1]*(0.4330127018922193*fr[1]-0.25*fr[0]); 
  incr[6] = dxvr[1]*(0.4330127018922193*fr[11]-0.25*fr[7])+wr[1]*(1.5*fr[6]-0.8660254037844386*fr[3]); 
  incr[7] = wr[1]*(0.5*fr[7]-0.8660254037844386*fr[11])+dxvr[1]*(0.1443375672974065*fr[3]-0.25*fr[6]); 
  incr[8] = dxvr[1]*(0.4330127018922193*fr[12]-0.25*fr[9])+wr[1]*(1.5*fr[8]-0.8660254037844386*fr[4]); 
  incr[9] = wr[1]*(0.5*fr[9]-0.8660254037844386*fr[12])+dxvr[1]*(0.1443375672974065*fr[4]-0.25*fr[8]); 
  incr[10] = dxvr[1]*(0.1443375672974065*fr[14]-0.25*fr[15])+wr[1]*(0.5*fr[10]-0.8660254037844386*fr[13]); 
  incr[11] = wr[1]*(1.5*fr[11]-0.8660254037844386*fr[7])+dxvr[1]*(0.4330127018922193*fr[6]-0.25*fr[3]); 
  incr[12] = wr[1]*(1.5*fr[12]-0.8660254037844386*fr[9])+dxvr[1]*(0.4330127018922193*fr[8]-0.25*fr[4]); 
  incr[13] = dxvr[1]*(0.4330127018922193*fr[15]-0.25*fr[14])+wr[1]*(1.5*fr[13]-0.8660254037844386*fr[10]); 
  incr[14] = wr[1]*(0.5*fr[14]-0.8660254037844386*fr[15])+dxvr[1]*(0.1443375672974065*fr[10]-0.25*fr[13]); 
  incr[15] = wr[1]*(1.5*fr[15]-0.8660254037844386*fr[14])+dxvr[1]*(0.4330127018922193*fr[13]-0.25*fr[10]); 

  outr[0] += incr[0]*rdxr2; 
  outr[1] += incr[1]*rdxr2; 
  outr[2] += incr[2]*rdxr2; 
  outr[3] += incr[3]*rdxr2; 
  outr[4] += incr[4]*rdxr2; 
  outr[5] += incr[5]*rdxr2; 
  outr[6] += incr[6]*rdxr2; 
  outr[7] += incr[7]*rdxr2; 
  outr[8] += incr[8]*rdxr2; 
  outr[9] += incr[9]*rdxr2; 
  outr[10] += incr[10]*rdxr2; 
  outr[11] += incr[11]*rdxr2; 
  outr[12] += incr[12]*rdxr2; 
  outr[13] += incr[13]*rdxr2; 
  outr[14] += incr[14]*rdxr2; 
  outr[15] += incr[15]*rdxr2; 

  outl[0] += -1.0*incr[0]*rdxl2; 
  outl[1] += incr[1]*rdxl2; 
  outl[2] += -1.0*incr[2]*rdxl2; 
  outl[3] += -1.0*incr[3]*rdxl2; 
  outl[4] += -1.0*incr[4]*rdxl2; 
  outl[5] += incr[5]*rdxl2; 
  outl[6] += incr[6]*rdxl2; 
  outl[7] += -1.0*incr[7]*rdxl2; 
  outl[8] += incr[8]*rdxl2; 
  outl[9] += -1.0*incr[9]*rdxl2; 
  outl[10] += -1.0*incr[10]*rdxl2; 
  outl[11] += incr[11]*rdxl2; 
  outl[12] += incr[12]*rdxl2; 
  outl[13] += incr[13]*rdxl2; 
  outl[14] += -1.0*incr[14]*rdxl2; 
  outl[15] += incr[15]*rdxl2; 
  } 
} 
