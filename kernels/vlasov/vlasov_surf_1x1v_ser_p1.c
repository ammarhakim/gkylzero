#include <gkyl_vlasov_kernels.h> 
void vlasov_surfx_1x1v_ser_p1(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells.
  gkyl_real rdxl2 = 2.0/dxvl[0]; 
  gkyl_real rdxr2 = 2.0/dxvr[0]; 

  gkyl_real incr[4]; 

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
gkyl_real vlasov_surfvx_1x1v_ser_p1(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real amax, const gkyl_real *qmem, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // amax:      amax in global lax flux.
  // qmem:      q/m*EM fields.
  // fl/fr:     Distribution function in left/right cells 
  // outl/outr: output distribution function in left/right cells 
  // returns abs(amid) for use in determining amax in cfl and global lax flux 
  gkyl_real dv10l = 2/dxvl[1]; 
  gkyl_real dv10r = 2/dxvr[1]; 

  const gkyl_real dv1 = dxvr[1], wv1 = wr[1]; 
  const gkyl_real *E0 = &qmem[0]; 

  gkyl_real Ghat[2]; 
  gkyl_real favg[2]; 
  gkyl_real alpha[2]; 

  favg[0] = (-1.224744871391589*fr[2])+1.224744871391589*fl[2]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0]; 
  favg[1] = (-1.224744871391589*fr[3])+1.224744871391589*fl[3]+0.7071067811865475*fr[1]+0.7071067811865475*fl[1]; 

  alpha[0] = E0[0]; 
  alpha[1] = E0[1]; 

  gkyl_real amid = 0.7071067811865475*alpha[0]; 

  Ghat[0] = 0.3535533905932737*((1.732050807568877*(fr[2]+fl[2])-1.0*fr[0]+fl[0])*amax+alpha[1]*favg[1]+alpha[0]*favg[0]); 
  Ghat[1] = 0.3535533905932737*((1.732050807568877*(fr[3]+fl[3])-1.0*fr[1]+fl[1])*amax+alpha[0]*favg[1]+favg[0]*alpha[1]); 

  outr[0] += 0.7071067811865475*Ghat[0]*dv10r; 
  outr[1] += 0.7071067811865475*Ghat[1]*dv10r; 
  outr[2] += -1.224744871391589*Ghat[0]*dv10r; 
  outr[3] += -1.224744871391589*Ghat[1]*dv10r; 

  outl[0] += -0.7071067811865475*Ghat[0]*dv10l; 
  outl[1] += -0.7071067811865475*Ghat[1]*dv10l; 
  outl[2] += -1.224744871391589*Ghat[0]*dv10l; 
  outl[3] += -1.224744871391589*Ghat[1]*dv10l; 

  return fabs(amid); 
} 
