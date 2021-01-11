#include <gkyl_vlasov_kernels.h> 
void vlasov_surfx_1x2v_ser_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells.
  double rdxl2 = 2.0/dxvl[0]; 
  double rdxr2 = 2.0/dxvr[0]; 

  double incr[8]; 

  if (wr[1]>0) { 
  incr[0] = dxvl[1]*(0.25*fl[4]+0.1443375672974065*fl[2])+(0.8660254037844386*fl[1]+0.5*fl[0])*wl[1]; 
  incr[1] = dxvl[1]*((-0.4330127018922193*fl[4])-0.25*fl[2])+((-1.5*fl[1])-0.8660254037844386*fl[0])*wl[1]; 
  incr[2] = wl[1]*(0.8660254037844386*fl[4]+0.5*fl[2])+dxvl[1]*(0.25*fl[1]+0.1443375672974065*fl[0]); 
  incr[3] = dxvl[1]*(0.25*fl[7]+0.1443375672974065*fl[6])+wl[1]*(0.8660254037844386*fl[5]+0.5*fl[3]); 
  incr[4] = wl[1]*((-1.5*fl[4])-0.8660254037844386*fl[2])+dxvl[1]*((-0.4330127018922193*fl[1])-0.25*fl[0]); 
  incr[5] = dxvl[1]*((-0.4330127018922193*fl[7])-0.25*fl[6])+wl[1]*((-1.5*fl[5])-0.8660254037844386*fl[3]); 
  incr[6] = wl[1]*(0.8660254037844386*fl[7]+0.5*fl[6])+dxvl[1]*(0.25*fl[5]+0.1443375672974065*fl[3]); 
  incr[7] = wl[1]*((-1.5*fl[7])-0.8660254037844386*fl[6])+dxvl[1]*((-0.4330127018922193*fl[5])-0.25*fl[3]); 

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
  outl[3] += -1.0*incr[3]*rdxl2; 
  outl[4] += incr[4]*rdxl2; 
  outl[5] += incr[5]*rdxl2; 
  outl[6] += -1.0*incr[6]*rdxl2; 
  outl[7] += incr[7]*rdxl2; 
  } else { 
  incr[0] = dxvr[1]*(0.1443375672974065*fr[2]-0.25*fr[4])+(0.5*fr[0]-0.8660254037844386*fr[1])*wr[1]; 
  incr[1] = dxvr[1]*(0.4330127018922193*fr[4]-0.25*fr[2])+(1.5*fr[1]-0.8660254037844386*fr[0])*wr[1]; 
  incr[2] = wr[1]*(0.5*fr[2]-0.8660254037844386*fr[4])+dxvr[1]*(0.1443375672974065*fr[0]-0.25*fr[1]); 
  incr[3] = dxvr[1]*(0.1443375672974065*fr[6]-0.25*fr[7])+wr[1]*(0.5*fr[3]-0.8660254037844386*fr[5]); 
  incr[4] = wr[1]*(1.5*fr[4]-0.8660254037844386*fr[2])+dxvr[1]*(0.4330127018922193*fr[1]-0.25*fr[0]); 
  incr[5] = dxvr[1]*(0.4330127018922193*fr[7]-0.25*fr[6])+wr[1]*(1.5*fr[5]-0.8660254037844386*fr[3]); 
  incr[6] = wr[1]*(0.5*fr[6]-0.8660254037844386*fr[7])+dxvr[1]*(0.1443375672974065*fr[3]-0.25*fr[5]); 
  incr[7] = wr[1]*(1.5*fr[7]-0.8660254037844386*fr[6])+dxvr[1]*(0.4330127018922193*fr[5]-0.25*fr[3]); 

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
  outl[3] += -1.0*incr[3]*rdxl2; 
  outl[4] += incr[4]*rdxl2; 
  outl[5] += incr[5]*rdxl2; 
  outl[6] += -1.0*incr[6]*rdxl2; 
  outl[7] += incr[7]*rdxl2; 
  } 
} 
double vlasov_surfvx_1x2v_ser_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *qmem, const double *fl, const double *fr, double* restrict outl, double* restrict outr) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // amax:      amax in global lax flux.
  // qmem:      q/m*EM fields.
  // fl/fr:     Distribution function in left/right cells 
  // outl/outr: output distribution function in left/right cells 
  // returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv10l = 2/dxvl[1]; 
  double dv10r = 2/dxvr[1]; 

  const double dv1 = dxvr[1], wv1 = wr[1]; 
  const double dv2 = dxvr[2], wv2 = wr[2]; 
  const double *E0 = &qmem[0]; 
  const double *B2 = &qmem[10]; 

  double Ghat[4]; 
  double favg[4]; 
  double alpha[4]; 

  favg[0] = (-1.224744871391589*fr[2])+1.224744871391589*fl[2]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0]; 
  favg[1] = (-1.224744871391589*fr[4])+1.224744871391589*fl[4]+0.7071067811865475*fr[1]+0.7071067811865475*fl[1]; 
  favg[2] = (-1.224744871391589*fr[6])+1.224744871391589*fl[6]+0.7071067811865475*fr[3]+0.7071067811865475*fl[3]; 
  favg[3] = (-1.224744871391589*fr[7])+1.224744871391589*fl[7]+0.7071067811865475*fr[5]+0.7071067811865475*fl[5]; 

  alpha[0] = 1.414213562373095*(B2[0]*wv2+E0[0]); 
  alpha[1] = 1.414213562373095*(B2[1]*wv2+E0[1]); 
  alpha[2] = 0.408248290463863*B2[0]*dv2; 
  alpha[3] = 0.408248290463863*B2[1]*dv2; 

  double amid = 0.5*alpha[0]; 

  Ghat[0] = 0.3535533905932737*(1.732050807568877*(fr[2]+fl[2])-1.0*fr[0]+fl[0])*amax+0.25*(alpha[3]*favg[3]+alpha[2]*favg[2]+alpha[1]*favg[1]+alpha[0]*favg[0]); 
  Ghat[1] = 0.3535533905932737*(1.732050807568877*(fr[4]+fl[4])-1.0*fr[1]+fl[1])*amax+0.25*(alpha[2]*favg[3]+favg[2]*alpha[3]+alpha[0]*favg[1]+favg[0]*alpha[1]); 
  Ghat[2] = 0.3535533905932737*(1.732050807568877*(fr[6]+fl[6])-1.0*fr[3]+fl[3])*amax+0.25*(alpha[1]*favg[3]+favg[1]*alpha[3]+alpha[0]*favg[2]+favg[0]*alpha[2]); 
  Ghat[3] = 0.3535533905932737*(1.732050807568877*(fr[7]+fl[7])-1.0*fr[5]+fl[5])*amax+0.25*(alpha[0]*favg[3]+favg[0]*alpha[3]+alpha[1]*favg[2]+favg[1]*alpha[2]); 

  outr[0] += 0.7071067811865475*Ghat[0]*dv10r; 
  outr[1] += 0.7071067811865475*Ghat[1]*dv10r; 
  outr[2] += -1.224744871391589*Ghat[0]*dv10r; 
  outr[3] += 0.7071067811865475*Ghat[2]*dv10r; 
  outr[4] += -1.224744871391589*Ghat[1]*dv10r; 
  outr[5] += 0.7071067811865475*Ghat[3]*dv10r; 
  outr[6] += -1.224744871391589*Ghat[2]*dv10r; 
  outr[7] += -1.224744871391589*Ghat[3]*dv10r; 

  outl[0] += -0.7071067811865475*Ghat[0]*dv10l; 
  outl[1] += -0.7071067811865475*Ghat[1]*dv10l; 
  outl[2] += -1.224744871391589*Ghat[0]*dv10l; 
  outl[3] += -0.7071067811865475*Ghat[2]*dv10l; 
  outl[4] += -1.224744871391589*Ghat[1]*dv10l; 
  outl[5] += -0.7071067811865475*Ghat[3]*dv10l; 
  outl[6] += -1.224744871391589*Ghat[2]*dv10l; 
  outl[7] += -1.224744871391589*Ghat[3]*dv10l; 

  return fabs(amid); 
} 
double vlasov_surfvy_1x2v_ser_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *qmem, const double *fl, const double *fr, double* restrict outl, double* restrict outr) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // amax:      amax in global lax flux.
  // qmem:      q/m*EM fields.
  // fl/fr:     Distribution function in left/right cells 
  // outl/outr: output distribution function in left/right cells 
  // returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv11l = 2/dxvl[2]; 
  double dv11r = 2/dxvr[2]; 

  const double dv1 = dxvr[1], wv1 = wr[1]; 
  const double dv2 = dxvr[2], wv2 = wr[2]; 
  const double *E1 = &qmem[2]; 
  const double *B2 = &qmem[10]; 

  double Ghat[4]; 
  double favg[4]; 
  double alpha[4]; 

  favg[0] = (-1.224744871391589*fr[3])+1.224744871391589*fl[3]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0]; 
  favg[1] = (-1.224744871391589*fr[5])+1.224744871391589*fl[5]+0.7071067811865475*fr[1]+0.7071067811865475*fl[1]; 
  favg[2] = (-1.224744871391589*fr[6])+1.224744871391589*fl[6]+0.7071067811865475*fr[2]+0.7071067811865475*fl[2]; 
  favg[3] = (-1.224744871391589*fr[7])+1.224744871391589*fl[7]+0.7071067811865475*fr[4]+0.7071067811865475*fl[4]; 

  alpha[0] = 1.414213562373095*E1[0]-1.414213562373095*B2[0]*wv1; 
  alpha[1] = 1.414213562373095*E1[1]-1.414213562373095*B2[1]*wv1; 
  alpha[2] = -0.408248290463863*B2[0]*dv1; 
  alpha[3] = -0.408248290463863*B2[1]*dv1; 

  double amid = 0.5*alpha[0]; 

  Ghat[0] = 0.3535533905932737*(1.732050807568877*(fr[3]+fl[3])-1.0*fr[0]+fl[0])*amax+0.25*(alpha[3]*favg[3]+alpha[2]*favg[2]+alpha[1]*favg[1]+alpha[0]*favg[0]); 
  Ghat[1] = 0.3535533905932737*(1.732050807568877*(fr[5]+fl[5])-1.0*fr[1]+fl[1])*amax+0.25*(alpha[2]*favg[3]+favg[2]*alpha[3]+alpha[0]*favg[1]+favg[0]*alpha[1]); 
  Ghat[2] = 0.3535533905932737*(1.732050807568877*(fr[6]+fl[6])-1.0*fr[2]+fl[2])*amax+0.25*(alpha[1]*favg[3]+favg[1]*alpha[3]+alpha[0]*favg[2]+favg[0]*alpha[2]); 
  Ghat[3] = 0.3535533905932737*(1.732050807568877*(fr[7]+fl[7])-1.0*fr[4]+fl[4])*amax+0.25*(alpha[0]*favg[3]+favg[0]*alpha[3]+alpha[1]*favg[2]+favg[1]*alpha[2]); 

  outr[0] += 0.7071067811865475*Ghat[0]*dv11r; 
  outr[1] += 0.7071067811865475*Ghat[1]*dv11r; 
  outr[2] += 0.7071067811865475*Ghat[2]*dv11r; 
  outr[3] += -1.224744871391589*Ghat[0]*dv11r; 
  outr[4] += 0.7071067811865475*Ghat[3]*dv11r; 
  outr[5] += -1.224744871391589*Ghat[1]*dv11r; 
  outr[6] += -1.224744871391589*Ghat[2]*dv11r; 
  outr[7] += -1.224744871391589*Ghat[3]*dv11r; 

  outl[0] += -0.7071067811865475*Ghat[0]*dv11l; 
  outl[1] += -0.7071067811865475*Ghat[1]*dv11l; 
  outl[2] += -0.7071067811865475*Ghat[2]*dv11l; 
  outl[3] += -1.224744871391589*Ghat[0]*dv11l; 
  outl[4] += -0.7071067811865475*Ghat[3]*dv11l; 
  outl[5] += -1.224744871391589*Ghat[1]*dv11l; 
  outl[6] += -1.224744871391589*Ghat[2]*dv11l; 
  outl[7] += -1.224744871391589*Ghat[3]*dv11l; 

  return fabs(amid); 
} 
