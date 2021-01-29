#include <gkyl_vlasov_kernels.h> 
void vlasov_surfx_1x1v_ser_p2(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells.
  gkyl_real rdxl2 = 2.0/dxvl[0]; 
  gkyl_real rdxr2 = 2.0/dxvr[0]; 

  gkyl_real incr[8]; 

  if (wr[1]>0) { 
  incr[0] = dxvl[1]*(0.3227486121839514*fl[6]+0.25*fl[3]+0.1443375672974065*fl[2])+wl[1]*(1.118033988749895*fl[4]+0.8660254037844386*fl[1]+0.5*fl[0]); 
  incr[1] = dxvl[1]*((-0.5590169943749476*fl[6])-0.4330127018922193*fl[3]-0.25*fl[2])+wl[1]*((-1.936491673103709*fl[4])-1.5*fl[1]-0.8660254037844386*fl[0]); 
  incr[2] = dxvl[1]*(0.223606797749979*fl[7]+0.1290994448735806*fl[5]+0.3227486121839515*fl[4]+0.25*fl[1]+0.1443375672974065*fl[0])+wl[1]*(1.118033988749895*fl[6]+0.8660254037844386*fl[3]+0.5*fl[2]); 
  incr[3] = dxvl[1]*((-0.3872983346207417*fl[7])-0.223606797749979*fl[5]-0.5590169943749475*fl[4]-0.4330127018922193*fl[1]-0.25*fl[0])+wl[1]*((-1.936491673103709*fl[6])-1.5*fl[3]-0.8660254037844386*fl[2]); 
  incr[4] = dxvl[1]*(0.7216878364870323*fl[6]+0.5590169943749475*fl[3]+0.3227486121839515*fl[2])+wl[1]*(2.5*fl[4]+1.936491673103709*fl[1]+1.118033988749895*fl[0]); 
  incr[5] = wl[1]*(0.8660254037844387*fl[7]+0.5*fl[5])+dxvl[1]*(0.2886751345948129*fl[6]+0.223606797749979*fl[3]+0.1290994448735806*fl[2]); 
  incr[6] = dxvl[1]*(0.5*fl[7]+0.2886751345948129*fl[5]+0.7216878364870323*fl[4]+0.5590169943749476*fl[1]+0.3227486121839514*fl[0])+wl[1]*(2.5*fl[6]+1.936491673103709*fl[3]+1.118033988749895*fl[2]); 
  incr[7] = wl[1]*((-1.5*fl[7])-0.8660254037844387*fl[5])+dxvl[1]*((-0.5*fl[6])-0.3872983346207417*fl[3]-0.223606797749979*fl[2]); 

  outr[0] += incr[0]*rdxr2; 
  outr[1] += incr[1]*rdxr2; 
  outr[2] += incr[2]*rdxr2; 
  outr[3] += incr[3]*rdxr2; 
  outr[4] += incr[4]*rdxr2; 
  outr[5] += incr[5]*rdxr2; 
  outr[6] += incr[6]*rdxr2; 
  outr[7] += incr[7]*rdxr2; 

  outl[0] += -1.0*incr[0]*rdxl2; 
  outl[1] += incr[1]*rdxl2; 
  outl[2] += -1.0*incr[2]*rdxl2; 
  outl[3] += incr[3]*rdxl2; 
  outl[4] += -1.0*incr[4]*rdxl2; 
  outl[5] += -1.0*incr[5]*rdxl2; 
  outl[6] += -1.0*incr[6]*rdxl2; 
  outl[7] += incr[7]*rdxl2; 
  } else { 
  incr[0] = dxvr[1]*(0.3227486121839514*fr[6]-0.25*fr[3]+0.1443375672974065*fr[2])+wr[1]*(1.118033988749895*fr[4]-0.8660254037844386*fr[1]+0.5*fr[0]); 
  incr[1] = dxvr[1]*((-0.5590169943749476*fr[6])+0.4330127018922193*fr[3]-0.25*fr[2])+wr[1]*((-1.936491673103709*fr[4])+1.5*fr[1]-0.8660254037844386*fr[0]); 
  incr[2] = dxvr[1]*((-0.223606797749979*fr[7])+0.1290994448735806*fr[5]+0.3227486121839515*fr[4]-0.25*fr[1]+0.1443375672974065*fr[0])+wr[1]*(1.118033988749895*fr[6]-0.8660254037844386*fr[3]+0.5*fr[2]); 
  incr[3] = dxvr[1]*(0.3872983346207417*fr[7]-0.223606797749979*fr[5]-0.5590169943749475*fr[4]+0.4330127018922193*fr[1]-0.25*fr[0])+wr[1]*((-1.936491673103709*fr[6])+1.5*fr[3]-0.8660254037844386*fr[2]); 
  incr[4] = dxvr[1]*(0.7216878364870323*fr[6]-0.5590169943749475*fr[3]+0.3227486121839515*fr[2])+wr[1]*(2.5*fr[4]-1.936491673103709*fr[1]+1.118033988749895*fr[0]); 
  incr[5] = wr[1]*(0.5*fr[5]-0.8660254037844387*fr[7])+dxvr[1]*(0.2886751345948129*fr[6]-0.223606797749979*fr[3]+0.1290994448735806*fr[2]); 
  incr[6] = dxvr[1]*((-0.5*fr[7])+0.2886751345948129*fr[5]+0.7216878364870323*fr[4]-0.5590169943749476*fr[1]+0.3227486121839514*fr[0])+wr[1]*(2.5*fr[6]-1.936491673103709*fr[3]+1.118033988749895*fr[2]); 
  incr[7] = wr[1]*(1.5*fr[7]-0.8660254037844387*fr[5])+dxvr[1]*((-0.5*fr[6])+0.3872983346207417*fr[3]-0.223606797749979*fr[2]); 

  outr[0] += incr[0]*rdxr2; 
  outr[1] += incr[1]*rdxr2; 
  outr[2] += incr[2]*rdxr2; 
  outr[3] += incr[3]*rdxr2; 
  outr[4] += incr[4]*rdxr2; 
  outr[5] += incr[5]*rdxr2; 
  outr[6] += incr[6]*rdxr2; 
  outr[7] += incr[7]*rdxr2; 

  outl[0] += -1.0*incr[0]*rdxl2; 
  outl[1] += incr[1]*rdxl2; 
  outl[2] += -1.0*incr[2]*rdxl2; 
  outl[3] += incr[3]*rdxl2; 
  outl[4] += -1.0*incr[4]*rdxl2; 
  outl[5] += -1.0*incr[5]*rdxl2; 
  outl[6] += -1.0*incr[6]*rdxl2; 
  outl[7] += incr[7]*rdxl2; 
  } 
} 
gkyl_real vlasov_surfvx_1x1v_ser_p2(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real amax, const gkyl_real *qmem, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr) 
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

  gkyl_real Ghat[3]; 
  gkyl_real favg[3]; 
  gkyl_real alpha[3]; 

  favg[0] = 1.58113883008419*fr[5]+1.58113883008419*fl[5]-1.224744871391589*fr[2]+1.224744871391589*fl[2]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0]; 
  favg[1] = 1.58113883008419*fr[7]+1.58113883008419*fl[7]-1.224744871391589*fr[3]+1.224744871391589*fl[3]+0.7071067811865475*fr[1]+0.7071067811865475*fl[1]; 
  favg[2] = (-1.224744871391589*fr[6])+1.224744871391589*fl[6]+0.7071067811865475*fr[4]+0.7071067811865475*fl[4]; 

  alpha[0] = E0[0]; 
  alpha[1] = E0[1]; 
  alpha[2] = E0[2]; 

  gkyl_real amid = 0.7071067811865475*alpha[0]-0.7905694150420947*alpha[2]; 

  Ghat[0] = 0.3535533905932737*(alpha[2]*favg[2]+alpha[1]*favg[1]+alpha[0]*favg[0])-0.3535533905932737*(2.23606797749979*fr[5]-1.0*(2.23606797749979*fl[5]+1.732050807568877*(fr[2]+fl[2]))+fr[0]-1.0*fl[0])*amax; 
  Ghat[1] = 0.07071067811865475*(4.47213595499958*(alpha[1]*favg[2]+favg[1]*alpha[2])+5.0*(alpha[0]*favg[1]+favg[0]*alpha[1]))-0.1178511301977579*(1.732050807568877*(3.872983346207417*fr[7]-1.0*(3.872983346207417*fl[7]+3.0*(fr[3]+fl[3])))+3.0*(fr[1]-1.0*fl[1]))*amax; 
  Ghat[2] = 0.07071067811865475*(8.660254037844387*(fr[6]+fl[6])+5.0*(fl[4]-1.0*fr[4]))*amax+0.01010152544552211*((22.3606797749979*alpha[2]+35.0*alpha[0])*favg[2]+35.0*favg[0]*alpha[2]+31.30495168499706*alpha[1]*favg[1]); 

  outr[0] += 0.7071067811865475*Ghat[0]*dv10r; 
  outr[1] += 0.7071067811865475*Ghat[1]*dv10r; 
  outr[2] += -1.224744871391589*Ghat[0]*dv10r; 
  outr[3] += -1.224744871391589*Ghat[1]*dv10r; 
  outr[4] += 0.7071067811865475*Ghat[2]*dv10r; 
  outr[5] += 1.58113883008419*Ghat[0]*dv10r; 
  outr[6] += -1.224744871391589*Ghat[2]*dv10r; 
  outr[7] += 1.58113883008419*Ghat[1]*dv10r; 

  outl[0] += -0.7071067811865475*Ghat[0]*dv10l; 
  outl[1] += -0.7071067811865475*Ghat[1]*dv10l; 
  outl[2] += -1.224744871391589*Ghat[0]*dv10l; 
  outl[3] += -1.224744871391589*Ghat[1]*dv10l; 
  outl[4] += -0.7071067811865475*Ghat[2]*dv10l; 
  outl[5] += -1.58113883008419*Ghat[0]*dv10l; 
  outl[6] += -1.224744871391589*Ghat[2]*dv10l; 
  outl[7] += -1.58113883008419*Ghat[1]*dv10l; 

  return fabs(amid); 
} 
