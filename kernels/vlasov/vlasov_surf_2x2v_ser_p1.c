#include <gkyl_vlasov_kernels.h> 
void vlasov_surfx_2x2v_ser_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells.
  double rdxl2 = 2.0/dxvl[0]; 
  double rdxr2 = 2.0/dxvr[0]; 

  double incr[16]; 

  if (wr[2]>0) { 
  incr[0] = dxvl[2]*(0.25*fl[6]+0.1443375672974065*fl[3])+(0.8660254037844386*fl[1]+0.5*fl[0])*wl[2]; 
  incr[1] = dxvl[2]*((-0.4330127018922193*fl[6])-0.25*fl[3])+((-1.5*fl[1])-0.8660254037844386*fl[0])*wl[2]; 
  incr[2] = dxvl[2]*(0.25*fl[11]+0.1443375672974065*fl[7])+wl[2]*(0.8660254037844386*fl[5]+0.5*fl[2]); 
  incr[3] = wl[2]*(0.8660254037844386*fl[6]+0.5*fl[3])+(0.25*fl[1]+0.1443375672974065*fl[0])*dxvl[2]; 
  incr[4] = dxvl[2]*(0.25*fl[13]+0.1443375672974065*fl[10])+wl[2]*(0.8660254037844386*fl[8]+0.5*fl[4]); 
  incr[5] = dxvl[2]*((-0.4330127018922193*fl[11])-0.25*fl[7])+wl[2]*((-1.5*fl[5])-0.8660254037844386*fl[2]); 
  incr[6] = wl[2]*((-1.5*fl[6])-0.8660254037844386*fl[3])+((-0.4330127018922193*fl[1])-0.25*fl[0])*dxvl[2]; 
  incr[7] = wl[2]*(0.8660254037844386*fl[11]+0.5*fl[7])+dxvl[2]*(0.25*fl[5]+0.1443375672974065*fl[2]); 
  incr[8] = dxvl[2]*((-0.4330127018922193*fl[13])-0.25*fl[10])+wl[2]*((-1.5*fl[8])-0.8660254037844386*fl[4]); 
  incr[9] = dxvl[2]*(0.25*fl[15]+0.1443375672974065*fl[14])+wl[2]*(0.8660254037844386*fl[12]+0.5*fl[9]); 
  incr[10] = wl[2]*(0.8660254037844386*fl[13]+0.5*fl[10])+dxvl[2]*(0.25*fl[8]+0.1443375672974065*fl[4]); 
  incr[11] = wl[2]*((-1.5*fl[11])-0.8660254037844386*fl[7])+dxvl[2]*((-0.4330127018922193*fl[5])-0.25*fl[2]); 
  incr[12] = dxvl[2]*((-0.4330127018922193*fl[15])-0.25*fl[14])+wl[2]*((-1.5*fl[12])-0.8660254037844386*fl[9]); 
  incr[13] = wl[2]*((-1.5*fl[13])-0.8660254037844386*fl[10])+dxvl[2]*((-0.4330127018922193*fl[8])-0.25*fl[4]); 
  incr[14] = wl[2]*(0.8660254037844386*fl[15]+0.5*fl[14])+dxvl[2]*(0.25*fl[12]+0.1443375672974065*fl[9]); 
  incr[15] = wl[2]*((-1.5*fl[15])-0.8660254037844386*fl[14])+dxvl[2]*((-0.4330127018922193*fl[12])-0.25*fl[9]); 

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
  incr[0] = dxvr[2]*(0.1443375672974065*fr[3]-0.25*fr[6])+(0.5*fr[0]-0.8660254037844386*fr[1])*wr[2]; 
  incr[1] = dxvr[2]*(0.4330127018922193*fr[6]-0.25*fr[3])+(1.5*fr[1]-0.8660254037844386*fr[0])*wr[2]; 
  incr[2] = dxvr[2]*(0.1443375672974065*fr[7]-0.25*fr[11])+wr[2]*(0.5*fr[2]-0.8660254037844386*fr[5]); 
  incr[3] = wr[2]*(0.5*fr[3]-0.8660254037844386*fr[6])+(0.1443375672974065*fr[0]-0.25*fr[1])*dxvr[2]; 
  incr[4] = dxvr[2]*(0.1443375672974065*fr[10]-0.25*fr[13])+wr[2]*(0.5*fr[4]-0.8660254037844386*fr[8]); 
  incr[5] = dxvr[2]*(0.4330127018922193*fr[11]-0.25*fr[7])+wr[2]*(1.5*fr[5]-0.8660254037844386*fr[2]); 
  incr[6] = wr[2]*(1.5*fr[6]-0.8660254037844386*fr[3])+(0.4330127018922193*fr[1]-0.25*fr[0])*dxvr[2]; 
  incr[7] = wr[2]*(0.5*fr[7]-0.8660254037844386*fr[11])+dxvr[2]*(0.1443375672974065*fr[2]-0.25*fr[5]); 
  incr[8] = dxvr[2]*(0.4330127018922193*fr[13]-0.25*fr[10])+wr[2]*(1.5*fr[8]-0.8660254037844386*fr[4]); 
  incr[9] = dxvr[2]*(0.1443375672974065*fr[14]-0.25*fr[15])+wr[2]*(0.5*fr[9]-0.8660254037844386*fr[12]); 
  incr[10] = wr[2]*(0.5*fr[10]-0.8660254037844386*fr[13])+dxvr[2]*(0.1443375672974065*fr[4]-0.25*fr[8]); 
  incr[11] = wr[2]*(1.5*fr[11]-0.8660254037844386*fr[7])+dxvr[2]*(0.4330127018922193*fr[5]-0.25*fr[2]); 
  incr[12] = dxvr[2]*(0.4330127018922193*fr[15]-0.25*fr[14])+wr[2]*(1.5*fr[12]-0.8660254037844386*fr[9]); 
  incr[13] = wr[2]*(1.5*fr[13]-0.8660254037844386*fr[10])+dxvr[2]*(0.4330127018922193*fr[8]-0.25*fr[4]); 
  incr[14] = wr[2]*(0.5*fr[14]-0.8660254037844386*fr[15])+dxvr[2]*(0.1443375672974065*fr[9]-0.25*fr[12]); 
  incr[15] = wr[2]*(1.5*fr[15]-0.8660254037844386*fr[14])+dxvr[2]*(0.4330127018922193*fr[12]-0.25*fr[9]); 

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
void vlasov_surfy_2x2v_ser_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells.
  double rdxl2 = 2.0/dxvl[1]; 
  double rdxr2 = 2.0/dxvr[1]; 

  double incr[16]; 

  if (wr[3]>0) { 
  incr[0] = dxvl[3]*(0.25*fl[9]+0.1443375672974065*fl[4])+(0.8660254037844386*fl[2]+0.5*fl[0])*wl[3]; 
  incr[1] = dxvl[3]*(0.25*fl[12]+0.1443375672974065*fl[8])+wl[3]*(0.8660254037844386*fl[5]+0.5*fl[1]); 
  incr[2] = dxvl[3]*((-0.4330127018922193*fl[9])-0.25*fl[4])+((-1.5*fl[2])-0.8660254037844386*fl[0])*wl[3]; 
  incr[3] = dxvl[3]*(0.25*fl[14]+0.1443375672974065*fl[10])+wl[3]*(0.8660254037844386*fl[7]+0.5*fl[3]); 
  incr[4] = wl[3]*(0.8660254037844386*fl[9]+0.5*fl[4])+(0.25*fl[2]+0.1443375672974065*fl[0])*dxvl[3]; 
  incr[5] = dxvl[3]*((-0.4330127018922193*fl[12])-0.25*fl[8])+wl[3]*((-1.5*fl[5])-0.8660254037844386*fl[1]); 
  incr[6] = dxvl[3]*(0.25*fl[15]+0.1443375672974065*fl[13])+wl[3]*(0.8660254037844386*fl[11]+0.5*fl[6]); 
  incr[7] = dxvl[3]*((-0.4330127018922193*fl[14])-0.25*fl[10])+wl[3]*((-1.5*fl[7])-0.8660254037844386*fl[3]); 
  incr[8] = wl[3]*(0.8660254037844386*fl[12]+0.5*fl[8])+dxvl[3]*(0.25*fl[5]+0.1443375672974065*fl[1]); 
  incr[9] = wl[3]*((-1.5*fl[9])-0.8660254037844386*fl[4])+((-0.4330127018922193*fl[2])-0.25*fl[0])*dxvl[3]; 
  incr[10] = wl[3]*(0.8660254037844386*fl[14]+0.5*fl[10])+dxvl[3]*(0.25*fl[7]+0.1443375672974065*fl[3]); 
  incr[11] = dxvl[3]*((-0.4330127018922193*fl[15])-0.25*fl[13])+wl[3]*((-1.5*fl[11])-0.8660254037844386*fl[6]); 
  incr[12] = wl[3]*((-1.5*fl[12])-0.8660254037844386*fl[8])+dxvl[3]*((-0.4330127018922193*fl[5])-0.25*fl[1]); 
  incr[13] = wl[3]*(0.8660254037844386*fl[15]+0.5*fl[13])+dxvl[3]*(0.25*fl[11]+0.1443375672974065*fl[6]); 
  incr[14] = wl[3]*((-1.5*fl[14])-0.8660254037844386*fl[10])+dxvl[3]*((-0.4330127018922193*fl[7])-0.25*fl[3]); 
  incr[15] = wl[3]*((-1.5*fl[15])-0.8660254037844386*fl[13])+dxvl[3]*((-0.4330127018922193*fl[11])-0.25*fl[6]); 

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
  outl[1] += -1.0*incr[1]*rdxl2; 
  outl[2] += incr[2]*rdxl2; 
  outl[3] += -1.0*incr[3]*rdxl2; 
  outl[4] += -1.0*incr[4]*rdxl2; 
  outl[5] += incr[5]*rdxl2; 
  outl[6] += -1.0*incr[6]*rdxl2; 
  outl[7] += incr[7]*rdxl2; 
  outl[8] += -1.0*incr[8]*rdxl2; 
  outl[9] += incr[9]*rdxl2; 
  outl[10] += -1.0*incr[10]*rdxl2; 
  outl[11] += incr[11]*rdxl2; 
  outl[12] += incr[12]*rdxl2; 
  outl[13] += -1.0*incr[13]*rdxl2; 
  outl[14] += incr[14]*rdxl2; 
  outl[15] += incr[15]*rdxl2; 
  } else { 
  incr[0] = dxvr[3]*(0.1443375672974065*fr[4]-0.25*fr[9])+(0.5*fr[0]-0.8660254037844386*fr[2])*wr[3]; 
  incr[1] = dxvr[3]*(0.1443375672974065*fr[8]-0.25*fr[12])+wr[3]*(0.5*fr[1]-0.8660254037844386*fr[5]); 
  incr[2] = dxvr[3]*(0.4330127018922193*fr[9]-0.25*fr[4])+(1.5*fr[2]-0.8660254037844386*fr[0])*wr[3]; 
  incr[3] = dxvr[3]*(0.1443375672974065*fr[10]-0.25*fr[14])+wr[3]*(0.5*fr[3]-0.8660254037844386*fr[7]); 
  incr[4] = wr[3]*(0.5*fr[4]-0.8660254037844386*fr[9])+(0.1443375672974065*fr[0]-0.25*fr[2])*dxvr[3]; 
  incr[5] = dxvr[3]*(0.4330127018922193*fr[12]-0.25*fr[8])+wr[3]*(1.5*fr[5]-0.8660254037844386*fr[1]); 
  incr[6] = dxvr[3]*(0.1443375672974065*fr[13]-0.25*fr[15])+wr[3]*(0.5*fr[6]-0.8660254037844386*fr[11]); 
  incr[7] = dxvr[3]*(0.4330127018922193*fr[14]-0.25*fr[10])+wr[3]*(1.5*fr[7]-0.8660254037844386*fr[3]); 
  incr[8] = wr[3]*(0.5*fr[8]-0.8660254037844386*fr[12])+dxvr[3]*(0.1443375672974065*fr[1]-0.25*fr[5]); 
  incr[9] = wr[3]*(1.5*fr[9]-0.8660254037844386*fr[4])+(0.4330127018922193*fr[2]-0.25*fr[0])*dxvr[3]; 
  incr[10] = wr[3]*(0.5*fr[10]-0.8660254037844386*fr[14])+dxvr[3]*(0.1443375672974065*fr[3]-0.25*fr[7]); 
  incr[11] = dxvr[3]*(0.4330127018922193*fr[15]-0.25*fr[13])+wr[3]*(1.5*fr[11]-0.8660254037844386*fr[6]); 
  incr[12] = wr[3]*(1.5*fr[12]-0.8660254037844386*fr[8])+dxvr[3]*(0.4330127018922193*fr[5]-0.25*fr[1]); 
  incr[13] = wr[3]*(0.5*fr[13]-0.8660254037844386*fr[15])+dxvr[3]*(0.1443375672974065*fr[6]-0.25*fr[11]); 
  incr[14] = wr[3]*(1.5*fr[14]-0.8660254037844386*fr[10])+dxvr[3]*(0.4330127018922193*fr[7]-0.25*fr[3]); 
  incr[15] = wr[3]*(1.5*fr[15]-0.8660254037844386*fr[13])+dxvr[3]*(0.4330127018922193*fr[11]-0.25*fr[6]); 

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
  outl[1] += -1.0*incr[1]*rdxl2; 
  outl[2] += incr[2]*rdxl2; 
  outl[3] += -1.0*incr[3]*rdxl2; 
  outl[4] += -1.0*incr[4]*rdxl2; 
  outl[5] += incr[5]*rdxl2; 
  outl[6] += -1.0*incr[6]*rdxl2; 
  outl[7] += incr[7]*rdxl2; 
  outl[8] += -1.0*incr[8]*rdxl2; 
  outl[9] += incr[9]*rdxl2; 
  outl[10] += -1.0*incr[10]*rdxl2; 
  outl[11] += incr[11]*rdxl2; 
  outl[12] += incr[12]*rdxl2; 
  outl[13] += -1.0*incr[13]*rdxl2; 
  outl[14] += incr[14]*rdxl2; 
  outl[15] += incr[15]*rdxl2; 
  } 
} 
double vlasov_surfvx_2x2v_ser_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *qmem, const double *fl, const double *fr, double* restrict outl, double* restrict outr) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // amax:      amax in global lax flux.
  // qmem:      q/m*EM fields.
  // fl/fr:     Distribution function in left/right cells 
  // outl/outr: output distribution function in left/right cells 
  // returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv10l = 2/dxvl[2]; 
  double dv10r = 2/dxvr[2]; 

  const double dv1 = dxvr[2], wv1 = wr[2]; 
  const double dv2 = dxvr[3], wv2 = wr[3]; 
  const double *E0 = &qmem[0]; 
  const double *B2 = &qmem[20]; 

  double Ghat[8]; 
  double favg[8]; 
  double alpha[8]; 

  favg[0] = (-1.224744871391589*fr[3])+1.224744871391589*fl[3]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0]; 
  favg[1] = (-1.224744871391589*fr[6])+1.224744871391589*fl[6]+0.7071067811865475*fr[1]+0.7071067811865475*fl[1]; 
  favg[2] = (-1.224744871391589*fr[7])+1.224744871391589*fl[7]+0.7071067811865475*fr[2]+0.7071067811865475*fl[2]; 
  favg[3] = (-1.224744871391589*fr[10])+1.224744871391589*fl[10]+0.7071067811865475*fr[4]+0.7071067811865475*fl[4]; 
  favg[4] = (-1.224744871391589*fr[11])+1.224744871391589*fl[11]+0.7071067811865475*fr[5]+0.7071067811865475*fl[5]; 
  favg[5] = (-1.224744871391589*fr[13])+1.224744871391589*fl[13]+0.7071067811865475*fr[8]+0.7071067811865475*fl[8]; 
  favg[6] = (-1.224744871391589*fr[14])+1.224744871391589*fl[14]+0.7071067811865475*fr[9]+0.7071067811865475*fl[9]; 
  favg[7] = (-1.224744871391589*fr[15])+1.224744871391589*fl[15]+0.7071067811865475*fr[12]+0.7071067811865475*fl[12]; 

  alpha[0] = 1.414213562373095*(B2[0]*wv2+E0[0]); 
  alpha[1] = 1.414213562373095*(B2[1]*wv2+E0[1]); 
  alpha[2] = 1.414213562373095*(B2[2]*wv2+E0[2]); 
  alpha[3] = 0.408248290463863*B2[0]*dv2; 
  alpha[4] = 1.414213562373095*(B2[3]*wv2+E0[3]); 
  alpha[5] = 0.408248290463863*B2[1]*dv2; 
  alpha[6] = 0.408248290463863*B2[2]*dv2; 
  alpha[7] = 0.408248290463863*B2[3]*dv2; 

  double amid = 0.3535533905932737*alpha[0]; 

  Ghat[0] = 0.3535533905932737*(1.732050807568877*(fr[3]+fl[3])-1.0*fr[0]+fl[0])*amax+0.1767766952966368*(alpha[7]*favg[7]+alpha[6]*favg[6]+alpha[5]*favg[5]+alpha[4]*favg[4]+alpha[3]*favg[3]+alpha[2]*favg[2]+alpha[1]*favg[1]+alpha[0]*favg[0]); 
  Ghat[1] = 0.3535533905932737*(1.732050807568877*(fr[6]+fl[6])-1.0*fr[1]+fl[1])*amax+0.1767766952966368*(alpha[6]*favg[7]+favg[6]*alpha[7]+alpha[3]*favg[5]+favg[3]*alpha[5]+alpha[2]*favg[4]+favg[2]*alpha[4]+alpha[0]*favg[1]+favg[0]*alpha[1]); 
  Ghat[2] = 0.3535533905932737*(1.732050807568877*(fr[7]+fl[7])-1.0*fr[2]+fl[2])*amax+0.1767766952966368*(alpha[5]*favg[7]+favg[5]*alpha[7]+alpha[3]*favg[6]+favg[3]*alpha[6]+alpha[1]*favg[4]+favg[1]*alpha[4]+alpha[0]*favg[2]+favg[0]*alpha[2]); 
  Ghat[3] = 0.3535533905932737*(1.732050807568877*(fr[10]+fl[10])-1.0*fr[4]+fl[4])*amax+0.1767766952966368*(alpha[4]*favg[7]+favg[4]*alpha[7]+alpha[2]*favg[6]+favg[2]*alpha[6]+alpha[1]*favg[5]+favg[1]*alpha[5]+alpha[0]*favg[3]+favg[0]*alpha[3]); 
  Ghat[4] = 0.3535533905932737*(1.732050807568877*(fr[11]+fl[11])-1.0*fr[5]+fl[5])*amax+0.1767766952966368*(alpha[3]*favg[7]+favg[3]*alpha[7]+alpha[5]*favg[6]+favg[5]*alpha[6]+alpha[0]*favg[4]+favg[0]*alpha[4]+alpha[1]*favg[2]+favg[1]*alpha[2]); 
  Ghat[5] = 0.3535533905932737*(1.732050807568877*(fr[13]+fl[13])-1.0*fr[8]+fl[8])*amax+0.1767766952966368*(alpha[2]*favg[7]+favg[2]*alpha[7]+alpha[4]*favg[6]+favg[4]*alpha[6]+alpha[0]*favg[5]+favg[0]*alpha[5]+alpha[1]*favg[3]+favg[1]*alpha[3]); 
  Ghat[6] = 0.3535533905932737*(1.732050807568877*(fr[14]+fl[14])-1.0*fr[9]+fl[9])*amax+0.1767766952966368*(alpha[1]*favg[7]+favg[1]*alpha[7]+alpha[0]*favg[6]+favg[0]*alpha[6]+alpha[4]*favg[5]+favg[4]*alpha[5]+alpha[2]*favg[3]+favg[2]*alpha[3]); 
  Ghat[7] = 0.3535533905932737*(1.732050807568877*(fr[15]+fl[15])-1.0*fr[12]+fl[12])*amax+0.1767766952966368*(alpha[0]*favg[7]+favg[0]*alpha[7]+alpha[1]*favg[6]+favg[1]*alpha[6]+alpha[2]*favg[5]+favg[2]*alpha[5]+alpha[3]*favg[4]+favg[3]*alpha[4]); 

  outr[0] += 0.7071067811865475*Ghat[0]*dv10r; 
  outr[1] += 0.7071067811865475*Ghat[1]*dv10r; 
  outr[2] += 0.7071067811865475*Ghat[2]*dv10r; 
  outr[3] += -1.224744871391589*Ghat[0]*dv10r; 
  outr[4] += 0.7071067811865475*Ghat[3]*dv10r; 
  outr[5] += 0.7071067811865475*Ghat[4]*dv10r; 
  outr[6] += -1.224744871391589*Ghat[1]*dv10r; 
  outr[7] += -1.224744871391589*Ghat[2]*dv10r; 
  outr[8] += 0.7071067811865475*Ghat[5]*dv10r; 
  outr[9] += 0.7071067811865475*Ghat[6]*dv10r; 
  outr[10] += -1.224744871391589*Ghat[3]*dv10r; 
  outr[11] += -1.224744871391589*Ghat[4]*dv10r; 
  outr[12] += 0.7071067811865475*Ghat[7]*dv10r; 
  outr[13] += -1.224744871391589*Ghat[5]*dv10r; 
  outr[14] += -1.224744871391589*Ghat[6]*dv10r; 
  outr[15] += -1.224744871391589*Ghat[7]*dv10r; 

  outl[0] += -0.7071067811865475*Ghat[0]*dv10l; 
  outl[1] += -0.7071067811865475*Ghat[1]*dv10l; 
  outl[2] += -0.7071067811865475*Ghat[2]*dv10l; 
  outl[3] += -1.224744871391589*Ghat[0]*dv10l; 
  outl[4] += -0.7071067811865475*Ghat[3]*dv10l; 
  outl[5] += -0.7071067811865475*Ghat[4]*dv10l; 
  outl[6] += -1.224744871391589*Ghat[1]*dv10l; 
  outl[7] += -1.224744871391589*Ghat[2]*dv10l; 
  outl[8] += -0.7071067811865475*Ghat[5]*dv10l; 
  outl[9] += -0.7071067811865475*Ghat[6]*dv10l; 
  outl[10] += -1.224744871391589*Ghat[3]*dv10l; 
  outl[11] += -1.224744871391589*Ghat[4]*dv10l; 
  outl[12] += -0.7071067811865475*Ghat[7]*dv10l; 
  outl[13] += -1.224744871391589*Ghat[5]*dv10l; 
  outl[14] += -1.224744871391589*Ghat[6]*dv10l; 
  outl[15] += -1.224744871391589*Ghat[7]*dv10l; 

  return fabs(amid); 
} 
double vlasov_surfvy_2x2v_ser_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *qmem, const double *fl, const double *fr, double* restrict outl, double* restrict outr) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // amax:      amax in global lax flux.
  // qmem:      q/m*EM fields.
  // fl/fr:     Distribution function in left/right cells 
  // outl/outr: output distribution function in left/right cells 
  // returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv11l = 2/dxvl[3]; 
  double dv11r = 2/dxvr[3]; 

  const double dv1 = dxvr[2], wv1 = wr[2]; 
  const double dv2 = dxvr[3], wv2 = wr[3]; 
  const double *E1 = &qmem[4]; 
  const double *B2 = &qmem[20]; 

  double Ghat[8]; 
  double favg[8]; 
  double alpha[8]; 

  favg[0] = (-1.224744871391589*fr[4])+1.224744871391589*fl[4]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0]; 
  favg[1] = (-1.224744871391589*fr[8])+1.224744871391589*fl[8]+0.7071067811865475*fr[1]+0.7071067811865475*fl[1]; 
  favg[2] = (-1.224744871391589*fr[9])+1.224744871391589*fl[9]+0.7071067811865475*fr[2]+0.7071067811865475*fl[2]; 
  favg[3] = (-1.224744871391589*fr[10])+1.224744871391589*fl[10]+0.7071067811865475*fr[3]+0.7071067811865475*fl[3]; 
  favg[4] = (-1.224744871391589*fr[12])+1.224744871391589*fl[12]+0.7071067811865475*fr[5]+0.7071067811865475*fl[5]; 
  favg[5] = (-1.224744871391589*fr[13])+1.224744871391589*fl[13]+0.7071067811865475*fr[6]+0.7071067811865475*fl[6]; 
  favg[6] = (-1.224744871391589*fr[14])+1.224744871391589*fl[14]+0.7071067811865475*fr[7]+0.7071067811865475*fl[7]; 
  favg[7] = (-1.224744871391589*fr[15])+1.224744871391589*fl[15]+0.7071067811865475*fr[11]+0.7071067811865475*fl[11]; 

  alpha[0] = 1.414213562373095*E1[0]-1.414213562373095*B2[0]*wv1; 
  alpha[1] = 1.414213562373095*E1[1]-1.414213562373095*B2[1]*wv1; 
  alpha[2] = 1.414213562373095*E1[2]-1.414213562373095*B2[2]*wv1; 
  alpha[3] = -0.408248290463863*B2[0]*dv1; 
  alpha[4] = 1.414213562373095*E1[3]-1.414213562373095*B2[3]*wv1; 
  alpha[5] = -0.408248290463863*B2[1]*dv1; 
  alpha[6] = -0.408248290463863*B2[2]*dv1; 
  alpha[7] = -0.408248290463863*B2[3]*dv1; 

  double amid = 0.3535533905932737*alpha[0]; 

  Ghat[0] = 0.3535533905932737*(1.732050807568877*(fr[4]+fl[4])-1.0*fr[0]+fl[0])*amax+0.1767766952966368*(alpha[7]*favg[7]+alpha[6]*favg[6]+alpha[5]*favg[5]+alpha[4]*favg[4]+alpha[3]*favg[3]+alpha[2]*favg[2]+alpha[1]*favg[1]+alpha[0]*favg[0]); 
  Ghat[1] = 0.3535533905932737*(1.732050807568877*(fr[8]+fl[8])-1.0*fr[1]+fl[1])*amax+0.1767766952966368*(alpha[6]*favg[7]+favg[6]*alpha[7]+alpha[3]*favg[5]+favg[3]*alpha[5]+alpha[2]*favg[4]+favg[2]*alpha[4]+alpha[0]*favg[1]+favg[0]*alpha[1]); 
  Ghat[2] = 0.3535533905932737*(1.732050807568877*(fr[9]+fl[9])-1.0*fr[2]+fl[2])*amax+0.1767766952966368*(alpha[5]*favg[7]+favg[5]*alpha[7]+alpha[3]*favg[6]+favg[3]*alpha[6]+alpha[1]*favg[4]+favg[1]*alpha[4]+alpha[0]*favg[2]+favg[0]*alpha[2]); 
  Ghat[3] = 0.3535533905932737*(1.732050807568877*(fr[10]+fl[10])-1.0*fr[3]+fl[3])*amax+0.1767766952966368*(alpha[4]*favg[7]+favg[4]*alpha[7]+alpha[2]*favg[6]+favg[2]*alpha[6]+alpha[1]*favg[5]+favg[1]*alpha[5]+alpha[0]*favg[3]+favg[0]*alpha[3]); 
  Ghat[4] = 0.3535533905932737*(1.732050807568877*(fr[12]+fl[12])-1.0*fr[5]+fl[5])*amax+0.1767766952966368*(alpha[3]*favg[7]+favg[3]*alpha[7]+alpha[5]*favg[6]+favg[5]*alpha[6]+alpha[0]*favg[4]+favg[0]*alpha[4]+alpha[1]*favg[2]+favg[1]*alpha[2]); 
  Ghat[5] = 0.3535533905932737*(1.732050807568877*(fr[13]+fl[13])-1.0*fr[6]+fl[6])*amax+0.1767766952966368*(alpha[2]*favg[7]+favg[2]*alpha[7]+alpha[4]*favg[6]+favg[4]*alpha[6]+alpha[0]*favg[5]+favg[0]*alpha[5]+alpha[1]*favg[3]+favg[1]*alpha[3]); 
  Ghat[6] = 0.3535533905932737*(1.732050807568877*(fr[14]+fl[14])-1.0*fr[7]+fl[7])*amax+0.1767766952966368*(alpha[1]*favg[7]+favg[1]*alpha[7]+alpha[0]*favg[6]+favg[0]*alpha[6]+alpha[4]*favg[5]+favg[4]*alpha[5]+alpha[2]*favg[3]+favg[2]*alpha[3]); 
  Ghat[7] = 0.3535533905932737*(1.732050807568877*(fr[15]+fl[15])-1.0*fr[11]+fl[11])*amax+0.1767766952966368*(alpha[0]*favg[7]+favg[0]*alpha[7]+alpha[1]*favg[6]+favg[1]*alpha[6]+alpha[2]*favg[5]+favg[2]*alpha[5]+alpha[3]*favg[4]+favg[3]*alpha[4]); 

  outr[0] += 0.7071067811865475*Ghat[0]*dv11r; 
  outr[1] += 0.7071067811865475*Ghat[1]*dv11r; 
  outr[2] += 0.7071067811865475*Ghat[2]*dv11r; 
  outr[3] += 0.7071067811865475*Ghat[3]*dv11r; 
  outr[4] += -1.224744871391589*Ghat[0]*dv11r; 
  outr[5] += 0.7071067811865475*Ghat[4]*dv11r; 
  outr[6] += 0.7071067811865475*Ghat[5]*dv11r; 
  outr[7] += 0.7071067811865475*Ghat[6]*dv11r; 
  outr[8] += -1.224744871391589*Ghat[1]*dv11r; 
  outr[9] += -1.224744871391589*Ghat[2]*dv11r; 
  outr[10] += -1.224744871391589*Ghat[3]*dv11r; 
  outr[11] += 0.7071067811865475*Ghat[7]*dv11r; 
  outr[12] += -1.224744871391589*Ghat[4]*dv11r; 
  outr[13] += -1.224744871391589*Ghat[5]*dv11r; 
  outr[14] += -1.224744871391589*Ghat[6]*dv11r; 
  outr[15] += -1.224744871391589*Ghat[7]*dv11r; 

  outl[0] += -0.7071067811865475*Ghat[0]*dv11l; 
  outl[1] += -0.7071067811865475*Ghat[1]*dv11l; 
  outl[2] += -0.7071067811865475*Ghat[2]*dv11l; 
  outl[3] += -0.7071067811865475*Ghat[3]*dv11l; 
  outl[4] += -1.224744871391589*Ghat[0]*dv11l; 
  outl[5] += -0.7071067811865475*Ghat[4]*dv11l; 
  outl[6] += -0.7071067811865475*Ghat[5]*dv11l; 
  outl[7] += -0.7071067811865475*Ghat[6]*dv11l; 
  outl[8] += -1.224744871391589*Ghat[1]*dv11l; 
  outl[9] += -1.224744871391589*Ghat[2]*dv11l; 
  outl[10] += -1.224744871391589*Ghat[3]*dv11l; 
  outl[11] += -0.7071067811865475*Ghat[7]*dv11l; 
  outl[12] += -1.224744871391589*Ghat[4]*dv11l; 
  outl[13] += -1.224744871391589*Ghat[5]*dv11l; 
  outl[14] += -1.224744871391589*Ghat[6]*dv11l; 
  outl[15] += -1.224744871391589*Ghat[7]*dv11l; 

  return fabs(amid); 
} 
