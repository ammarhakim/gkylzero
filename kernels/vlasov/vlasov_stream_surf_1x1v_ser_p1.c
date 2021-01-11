#include <gkyl_vlasov_kernels.h> 
void vlasov_stream_surfx_1x1v_ser_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells.
  double rdxl2 = 2.0/dxvl[0]; 
  double rdxr2 = 2.0/dxvr[0]; 

  double incr[4]; 

  if (wr[1]>0) { 
  incr[0] = dxvl[1]*(0.25*fl[3]+0.1443375672974065*fl[2])+(0.8660254037844386*fl[1]+0.5*fl[0])*wl[1]; 
  incr[1] = dxvl[1]*((-0.4330127018922193*fl[3])-0.25*fl[2])+((-1.5*fl[1])-0.8660254037844386*fl[0])*wl[1]; 
  incr[2] = wl[1]*(0.8660254037844386*fl[3]+0.5*fl[2])+dxvl[1]*(0.25*fl[1]+0.1443375672974065*fl[0]); 
  incr[3] = wl[1]*((-1.5*fl[3])-0.8660254037844386*fl[2])+dxvl[1]*((-0.4330127018922193*fl[1])-0.25*fl[0]); 

  outr[0] += incr[0]*rdxr2; 
  outr[1] += incr[1]*rdxr2; 
  outr[2] += incr[2]*rdxr2; 
  outr[3] += incr[3]*rdxr2; 

  outl[0] += -1.0*incr[0]*rdxl2; 
  outl[1] += incr[1]*rdxl2; 
  outl[2] += -1.0*incr[2]*rdxl2; 
  outl[3] += incr[3]*rdxl2; 
  } else { 
  incr[0] = dxvr[1]*(0.1443375672974065*fr[2]-0.25*fr[3])+(0.5*fr[0]-0.8660254037844386*fr[1])*wr[1]; 
  incr[1] = dxvr[1]*(0.4330127018922193*fr[3]-0.25*fr[2])+(1.5*fr[1]-0.8660254037844386*fr[0])*wr[1]; 
  incr[2] = wr[1]*(0.5*fr[2]-0.8660254037844386*fr[3])+dxvr[1]*(0.1443375672974065*fr[0]-0.25*fr[1]); 
  incr[3] = wr[1]*(1.5*fr[3]-0.8660254037844386*fr[2])+dxvr[1]*(0.4330127018922193*fr[1]-0.25*fr[0]); 

  outr[0] += incr[0]*rdxr2; 
  outr[1] += incr[1]*rdxr2; 
  outr[2] += incr[2]*rdxr2; 
  outr[3] += incr[3]*rdxr2; 

  outl[0] += -1.0*incr[0]*rdxl2; 
  outl[1] += incr[1]*rdxl2; 
  outl[2] += -1.0*incr[2]*rdxl2; 
  outl[3] += incr[3]*rdxl2; 
  } 
} 
