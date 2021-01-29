#include <gkyl_vlasov_kernels.h> 
void vlasov_surfx_2x3v_ser_p1(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells.
  gkyl_real rdxl2 = 2.0/dxvl[0]; 
  gkyl_real rdxr2 = 2.0/dxvr[0]; 

  gkyl_real incr[32]; 

  if (wr[2]>0) { 
  incr[0] = dxvl[2]*(0.25*fl[7]+0.1443375672974065*fl[3])+(0.8660254037844386*fl[1]+0.5*fl[0])*wl[2]; 
  incr[1] = dxvl[2]*((-0.4330127018922193*fl[7])-0.25*fl[3])+((-1.5*fl[1])-0.8660254037844386*fl[0])*wl[2]; 
  incr[2] = dxvl[2]*(0.25*fl[16]+0.1443375672974065*fl[8])+wl[2]*(0.8660254037844386*fl[6]+0.5*fl[2]); 
  incr[3] = wl[2]*(0.8660254037844386*fl[7]+0.5*fl[3])+(0.25*fl[1]+0.1443375672974065*fl[0])*dxvl[2]; 
  incr[4] = dxvl[2]*(0.25*fl[18]+0.1443375672974065*fl[11])+wl[2]*(0.8660254037844386*fl[9]+0.5*fl[4]); 
  incr[5] = dxvl[2]*(0.25*fl[21]+0.1443375672974065*fl[14])+wl[2]*(0.8660254037844386*fl[12]+0.5*fl[5]); 
  incr[6] = dxvl[2]*((-0.4330127018922193*fl[16])-0.25*fl[8])+wl[2]*((-1.5*fl[6])-0.8660254037844386*fl[2]); 
  incr[7] = wl[2]*((-1.5*fl[7])-0.8660254037844386*fl[3])+((-0.4330127018922193*fl[1])-0.25*fl[0])*dxvl[2]; 
  incr[8] = wl[2]*(0.8660254037844386*fl[16]+0.5*fl[8])+dxvl[2]*(0.25*fl[6]+0.1443375672974065*fl[2]); 
  incr[9] = dxvl[2]*((-0.4330127018922193*fl[18])-0.25*fl[11])+wl[2]*((-1.5*fl[9])-0.8660254037844386*fl[4]); 
  incr[10] = dxvl[2]*(0.25*fl[26]+0.1443375672974065*fl[19])+wl[2]*(0.8660254037844386*fl[17]+0.5*fl[10]); 
  incr[11] = wl[2]*(0.8660254037844386*fl[18]+0.5*fl[11])+dxvl[2]*(0.25*fl[9]+0.1443375672974065*fl[4]); 
  incr[12] = dxvl[2]*((-0.4330127018922193*fl[21])-0.25*fl[14])+wl[2]*((-1.5*fl[12])-0.8660254037844386*fl[5]); 
  incr[13] = dxvl[2]*(0.25*fl[27]+0.1443375672974065*fl[22])+wl[2]*(0.8660254037844386*fl[20]+0.5*fl[13]); 
  incr[14] = wl[2]*(0.8660254037844386*fl[21]+0.5*fl[14])+dxvl[2]*(0.25*fl[12]+0.1443375672974065*fl[5]); 
  incr[15] = dxvl[2]*(0.25*fl[29]+0.1443375672974065*fl[25])+wl[2]*(0.8660254037844386*fl[23]+0.5*fl[15]); 
  incr[16] = wl[2]*((-1.5*fl[16])-0.8660254037844386*fl[8])+dxvl[2]*((-0.4330127018922193*fl[6])-0.25*fl[2]); 
  incr[17] = dxvl[2]*((-0.4330127018922193*fl[26])-0.25*fl[19])+wl[2]*((-1.5*fl[17])-0.8660254037844386*fl[10]); 
  incr[18] = wl[2]*((-1.5*fl[18])-0.8660254037844386*fl[11])+dxvl[2]*((-0.4330127018922193*fl[9])-0.25*fl[4]); 
  incr[19] = wl[2]*(0.8660254037844386*fl[26]+0.5*fl[19])+dxvl[2]*(0.25*fl[17]+0.1443375672974065*fl[10]); 
  incr[20] = dxvl[2]*((-0.4330127018922193*fl[27])-0.25*fl[22])+wl[2]*((-1.5*fl[20])-0.8660254037844386*fl[13]); 
  incr[21] = wl[2]*((-1.5*fl[21])-0.8660254037844386*fl[14])+dxvl[2]*((-0.4330127018922193*fl[12])-0.25*fl[5]); 
  incr[22] = wl[2]*(0.8660254037844386*fl[27]+0.5*fl[22])+dxvl[2]*(0.25*fl[20]+0.1443375672974065*fl[13]); 
  incr[23] = dxvl[2]*((-0.4330127018922193*fl[29])-0.25*fl[25])+wl[2]*((-1.5*fl[23])-0.8660254037844386*fl[15]); 
  incr[24] = dxvl[2]*(0.25*fl[31]+0.1443375672974065*fl[30])+wl[2]*(0.8660254037844386*fl[28]+0.5*fl[24]); 
  incr[25] = wl[2]*(0.8660254037844386*fl[29]+0.5*fl[25])+dxvl[2]*(0.25*fl[23]+0.1443375672974065*fl[15]); 
  incr[26] = wl[2]*((-1.5*fl[26])-0.8660254037844386*fl[19])+dxvl[2]*((-0.4330127018922193*fl[17])-0.25*fl[10]); 
  incr[27] = wl[2]*((-1.5*fl[27])-0.8660254037844386*fl[22])+dxvl[2]*((-0.4330127018922193*fl[20])-0.25*fl[13]); 
  incr[28] = dxvl[2]*((-0.4330127018922193*fl[31])-0.25*fl[30])+wl[2]*((-1.5*fl[28])-0.8660254037844386*fl[24]); 
  incr[29] = wl[2]*((-1.5*fl[29])-0.8660254037844386*fl[25])+dxvl[2]*((-0.4330127018922193*fl[23])-0.25*fl[15]); 
  incr[30] = wl[2]*(0.8660254037844386*fl[31]+0.5*fl[30])+dxvl[2]*(0.25*fl[28]+0.1443375672974065*fl[24]); 
  incr[31] = wl[2]*((-1.5*fl[31])-0.8660254037844386*fl[30])+dxvl[2]*((-0.4330127018922193*fl[28])-0.25*fl[24]); 

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
  outr[16] += incr[16]*rdxr2; 
  outr[17] += incr[17]*rdxr2; 
  outr[18] += incr[18]*rdxr2; 
  outr[19] += incr[19]*rdxr2; 
  outr[20] += incr[20]*rdxr2; 
  outr[21] += incr[21]*rdxr2; 
  outr[22] += incr[22]*rdxr2; 
  outr[23] += incr[23]*rdxr2; 
  outr[24] += incr[24]*rdxr2; 
  outr[25] += incr[25]*rdxr2; 
  outr[26] += incr[26]*rdxr2; 
  outr[27] += incr[27]*rdxr2; 
  outr[28] += incr[28]*rdxr2; 
  outr[29] += incr[29]*rdxr2; 
  outr[30] += incr[30]*rdxr2; 
  outr[31] += incr[31]*rdxr2; 

  outl[0] += -1.0*incr[0]*rdxl2; 
  outl[1] += incr[1]*rdxl2; 
  outl[2] += -1.0*incr[2]*rdxl2; 
  outl[3] += -1.0*incr[3]*rdxl2; 
  outl[4] += -1.0*incr[4]*rdxl2; 
  outl[5] += -1.0*incr[5]*rdxl2; 
  outl[6] += incr[6]*rdxl2; 
  outl[7] += incr[7]*rdxl2; 
  outl[8] += -1.0*incr[8]*rdxl2; 
  outl[9] += incr[9]*rdxl2; 
  outl[10] += -1.0*incr[10]*rdxl2; 
  outl[11] += -1.0*incr[11]*rdxl2; 
  outl[12] += incr[12]*rdxl2; 
  outl[13] += -1.0*incr[13]*rdxl2; 
  outl[14] += -1.0*incr[14]*rdxl2; 
  outl[15] += -1.0*incr[15]*rdxl2; 
  outl[16] += incr[16]*rdxl2; 
  outl[17] += incr[17]*rdxl2; 
  outl[18] += incr[18]*rdxl2; 
  outl[19] += -1.0*incr[19]*rdxl2; 
  outl[20] += incr[20]*rdxl2; 
  outl[21] += incr[21]*rdxl2; 
  outl[22] += -1.0*incr[22]*rdxl2; 
  outl[23] += incr[23]*rdxl2; 
  outl[24] += -1.0*incr[24]*rdxl2; 
  outl[25] += -1.0*incr[25]*rdxl2; 
  outl[26] += incr[26]*rdxl2; 
  outl[27] += incr[27]*rdxl2; 
  outl[28] += incr[28]*rdxl2; 
  outl[29] += incr[29]*rdxl2; 
  outl[30] += -1.0*incr[30]*rdxl2; 
  outl[31] += incr[31]*rdxl2; 
  } else { 
  incr[0] = dxvr[2]*(0.1443375672974065*fr[3]-0.25*fr[7])+(0.5*fr[0]-0.8660254037844386*fr[1])*wr[2]; 
  incr[1] = dxvr[2]*(0.4330127018922193*fr[7]-0.25*fr[3])+(1.5*fr[1]-0.8660254037844386*fr[0])*wr[2]; 
  incr[2] = dxvr[2]*(0.1443375672974065*fr[8]-0.25*fr[16])+wr[2]*(0.5*fr[2]-0.8660254037844386*fr[6]); 
  incr[3] = wr[2]*(0.5*fr[3]-0.8660254037844386*fr[7])+(0.1443375672974065*fr[0]-0.25*fr[1])*dxvr[2]; 
  incr[4] = dxvr[2]*(0.1443375672974065*fr[11]-0.25*fr[18])+wr[2]*(0.5*fr[4]-0.8660254037844386*fr[9]); 
  incr[5] = dxvr[2]*(0.1443375672974065*fr[14]-0.25*fr[21])+wr[2]*(0.5*fr[5]-0.8660254037844386*fr[12]); 
  incr[6] = dxvr[2]*(0.4330127018922193*fr[16]-0.25*fr[8])+wr[2]*(1.5*fr[6]-0.8660254037844386*fr[2]); 
  incr[7] = wr[2]*(1.5*fr[7]-0.8660254037844386*fr[3])+(0.4330127018922193*fr[1]-0.25*fr[0])*dxvr[2]; 
  incr[8] = wr[2]*(0.5*fr[8]-0.8660254037844386*fr[16])+dxvr[2]*(0.1443375672974065*fr[2]-0.25*fr[6]); 
  incr[9] = dxvr[2]*(0.4330127018922193*fr[18]-0.25*fr[11])+wr[2]*(1.5*fr[9]-0.8660254037844386*fr[4]); 
  incr[10] = dxvr[2]*(0.1443375672974065*fr[19]-0.25*fr[26])+wr[2]*(0.5*fr[10]-0.8660254037844386*fr[17]); 
  incr[11] = wr[2]*(0.5*fr[11]-0.8660254037844386*fr[18])+dxvr[2]*(0.1443375672974065*fr[4]-0.25*fr[9]); 
  incr[12] = dxvr[2]*(0.4330127018922193*fr[21]-0.25*fr[14])+wr[2]*(1.5*fr[12]-0.8660254037844386*fr[5]); 
  incr[13] = dxvr[2]*(0.1443375672974065*fr[22]-0.25*fr[27])+wr[2]*(0.5*fr[13]-0.8660254037844386*fr[20]); 
  incr[14] = wr[2]*(0.5*fr[14]-0.8660254037844386*fr[21])+dxvr[2]*(0.1443375672974065*fr[5]-0.25*fr[12]); 
  incr[15] = dxvr[2]*(0.1443375672974065*fr[25]-0.25*fr[29])+wr[2]*(0.5*fr[15]-0.8660254037844386*fr[23]); 
  incr[16] = wr[2]*(1.5*fr[16]-0.8660254037844386*fr[8])+dxvr[2]*(0.4330127018922193*fr[6]-0.25*fr[2]); 
  incr[17] = dxvr[2]*(0.4330127018922193*fr[26]-0.25*fr[19])+wr[2]*(1.5*fr[17]-0.8660254037844386*fr[10]); 
  incr[18] = wr[2]*(1.5*fr[18]-0.8660254037844386*fr[11])+dxvr[2]*(0.4330127018922193*fr[9]-0.25*fr[4]); 
  incr[19] = wr[2]*(0.5*fr[19]-0.8660254037844386*fr[26])+dxvr[2]*(0.1443375672974065*fr[10]-0.25*fr[17]); 
  incr[20] = dxvr[2]*(0.4330127018922193*fr[27]-0.25*fr[22])+wr[2]*(1.5*fr[20]-0.8660254037844386*fr[13]); 
  incr[21] = wr[2]*(1.5*fr[21]-0.8660254037844386*fr[14])+dxvr[2]*(0.4330127018922193*fr[12]-0.25*fr[5]); 
  incr[22] = wr[2]*(0.5*fr[22]-0.8660254037844386*fr[27])+dxvr[2]*(0.1443375672974065*fr[13]-0.25*fr[20]); 
  incr[23] = dxvr[2]*(0.4330127018922193*fr[29]-0.25*fr[25])+wr[2]*(1.5*fr[23]-0.8660254037844386*fr[15]); 
  incr[24] = dxvr[2]*(0.1443375672974065*fr[30]-0.25*fr[31])+wr[2]*(0.5*fr[24]-0.8660254037844386*fr[28]); 
  incr[25] = wr[2]*(0.5*fr[25]-0.8660254037844386*fr[29])+dxvr[2]*(0.1443375672974065*fr[15]-0.25*fr[23]); 
  incr[26] = wr[2]*(1.5*fr[26]-0.8660254037844386*fr[19])+dxvr[2]*(0.4330127018922193*fr[17]-0.25*fr[10]); 
  incr[27] = wr[2]*(1.5*fr[27]-0.8660254037844386*fr[22])+dxvr[2]*(0.4330127018922193*fr[20]-0.25*fr[13]); 
  incr[28] = dxvr[2]*(0.4330127018922193*fr[31]-0.25*fr[30])+wr[2]*(1.5*fr[28]-0.8660254037844386*fr[24]); 
  incr[29] = wr[2]*(1.5*fr[29]-0.8660254037844386*fr[25])+dxvr[2]*(0.4330127018922193*fr[23]-0.25*fr[15]); 
  incr[30] = wr[2]*(0.5*fr[30]-0.8660254037844386*fr[31])+dxvr[2]*(0.1443375672974065*fr[24]-0.25*fr[28]); 
  incr[31] = wr[2]*(1.5*fr[31]-0.8660254037844386*fr[30])+dxvr[2]*(0.4330127018922193*fr[28]-0.25*fr[24]); 

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
  outr[16] += incr[16]*rdxr2; 
  outr[17] += incr[17]*rdxr2; 
  outr[18] += incr[18]*rdxr2; 
  outr[19] += incr[19]*rdxr2; 
  outr[20] += incr[20]*rdxr2; 
  outr[21] += incr[21]*rdxr2; 
  outr[22] += incr[22]*rdxr2; 
  outr[23] += incr[23]*rdxr2; 
  outr[24] += incr[24]*rdxr2; 
  outr[25] += incr[25]*rdxr2; 
  outr[26] += incr[26]*rdxr2; 
  outr[27] += incr[27]*rdxr2; 
  outr[28] += incr[28]*rdxr2; 
  outr[29] += incr[29]*rdxr2; 
  outr[30] += incr[30]*rdxr2; 
  outr[31] += incr[31]*rdxr2; 

  outl[0] += -1.0*incr[0]*rdxl2; 
  outl[1] += incr[1]*rdxl2; 
  outl[2] += -1.0*incr[2]*rdxl2; 
  outl[3] += -1.0*incr[3]*rdxl2; 
  outl[4] += -1.0*incr[4]*rdxl2; 
  outl[5] += -1.0*incr[5]*rdxl2; 
  outl[6] += incr[6]*rdxl2; 
  outl[7] += incr[7]*rdxl2; 
  outl[8] += -1.0*incr[8]*rdxl2; 
  outl[9] += incr[9]*rdxl2; 
  outl[10] += -1.0*incr[10]*rdxl2; 
  outl[11] += -1.0*incr[11]*rdxl2; 
  outl[12] += incr[12]*rdxl2; 
  outl[13] += -1.0*incr[13]*rdxl2; 
  outl[14] += -1.0*incr[14]*rdxl2; 
  outl[15] += -1.0*incr[15]*rdxl2; 
  outl[16] += incr[16]*rdxl2; 
  outl[17] += incr[17]*rdxl2; 
  outl[18] += incr[18]*rdxl2; 
  outl[19] += -1.0*incr[19]*rdxl2; 
  outl[20] += incr[20]*rdxl2; 
  outl[21] += incr[21]*rdxl2; 
  outl[22] += -1.0*incr[22]*rdxl2; 
  outl[23] += incr[23]*rdxl2; 
  outl[24] += -1.0*incr[24]*rdxl2; 
  outl[25] += -1.0*incr[25]*rdxl2; 
  outl[26] += incr[26]*rdxl2; 
  outl[27] += incr[27]*rdxl2; 
  outl[28] += incr[28]*rdxl2; 
  outl[29] += incr[29]*rdxl2; 
  outl[30] += -1.0*incr[30]*rdxl2; 
  outl[31] += incr[31]*rdxl2; 
  } 
} 
void vlasov_surfy_2x3v_ser_p1(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells.
  gkyl_real rdxl2 = 2.0/dxvl[1]; 
  gkyl_real rdxr2 = 2.0/dxvr[1]; 

  gkyl_real incr[32]; 

  if (wr[3]>0) { 
  incr[0] = dxvl[3]*(0.25*fl[10]+0.1443375672974065*fl[4])+(0.8660254037844386*fl[2]+0.5*fl[0])*wl[3]; 
  incr[1] = dxvl[3]*(0.25*fl[17]+0.1443375672974065*fl[9])+wl[3]*(0.8660254037844386*fl[6]+0.5*fl[1]); 
  incr[2] = dxvl[3]*((-0.4330127018922193*fl[10])-0.25*fl[4])+((-1.5*fl[2])-0.8660254037844386*fl[0])*wl[3]; 
  incr[3] = dxvl[3]*(0.25*fl[19]+0.1443375672974065*fl[11])+wl[3]*(0.8660254037844386*fl[8]+0.5*fl[3]); 
  incr[4] = wl[3]*(0.8660254037844386*fl[10]+0.5*fl[4])+(0.25*fl[2]+0.1443375672974065*fl[0])*dxvl[3]; 
  incr[5] = dxvl[3]*(0.25*fl[24]+0.1443375672974065*fl[15])+wl[3]*(0.8660254037844386*fl[13]+0.5*fl[5]); 
  incr[6] = dxvl[3]*((-0.4330127018922193*fl[17])-0.25*fl[9])+wl[3]*((-1.5*fl[6])-0.8660254037844386*fl[1]); 
  incr[7] = dxvl[3]*(0.25*fl[26]+0.1443375672974065*fl[18])+wl[3]*(0.8660254037844386*fl[16]+0.5*fl[7]); 
  incr[8] = dxvl[3]*((-0.4330127018922193*fl[19])-0.25*fl[11])+wl[3]*((-1.5*fl[8])-0.8660254037844386*fl[3]); 
  incr[9] = wl[3]*(0.8660254037844386*fl[17]+0.5*fl[9])+dxvl[3]*(0.25*fl[6]+0.1443375672974065*fl[1]); 
  incr[10] = wl[3]*((-1.5*fl[10])-0.8660254037844386*fl[4])+((-0.4330127018922193*fl[2])-0.25*fl[0])*dxvl[3]; 
  incr[11] = wl[3]*(0.8660254037844386*fl[19]+0.5*fl[11])+dxvl[3]*(0.25*fl[8]+0.1443375672974065*fl[3]); 
  incr[12] = dxvl[3]*(0.25*fl[28]+0.1443375672974065*fl[23])+wl[3]*(0.8660254037844386*fl[20]+0.5*fl[12]); 
  incr[13] = dxvl[3]*((-0.4330127018922193*fl[24])-0.25*fl[15])+wl[3]*((-1.5*fl[13])-0.8660254037844386*fl[5]); 
  incr[14] = dxvl[3]*(0.25*fl[30]+0.1443375672974065*fl[25])+wl[3]*(0.8660254037844386*fl[22]+0.5*fl[14]); 
  incr[15] = wl[3]*(0.8660254037844386*fl[24]+0.5*fl[15])+dxvl[3]*(0.25*fl[13]+0.1443375672974065*fl[5]); 
  incr[16] = dxvl[3]*((-0.4330127018922193*fl[26])-0.25*fl[18])+wl[3]*((-1.5*fl[16])-0.8660254037844386*fl[7]); 
  incr[17] = wl[3]*((-1.5*fl[17])-0.8660254037844386*fl[9])+dxvl[3]*((-0.4330127018922193*fl[6])-0.25*fl[1]); 
  incr[18] = wl[3]*(0.8660254037844386*fl[26]+0.5*fl[18])+dxvl[3]*(0.25*fl[16]+0.1443375672974065*fl[7]); 
  incr[19] = wl[3]*((-1.5*fl[19])-0.8660254037844386*fl[11])+dxvl[3]*((-0.4330127018922193*fl[8])-0.25*fl[3]); 
  incr[20] = dxvl[3]*((-0.4330127018922193*fl[28])-0.25*fl[23])+wl[3]*((-1.5*fl[20])-0.8660254037844386*fl[12]); 
  incr[21] = dxvl[3]*(0.25*fl[31]+0.1443375672974065*fl[29])+wl[3]*(0.8660254037844386*fl[27]+0.5*fl[21]); 
  incr[22] = dxvl[3]*((-0.4330127018922193*fl[30])-0.25*fl[25])+wl[3]*((-1.5*fl[22])-0.8660254037844386*fl[14]); 
  incr[23] = wl[3]*(0.8660254037844386*fl[28]+0.5*fl[23])+dxvl[3]*(0.25*fl[20]+0.1443375672974065*fl[12]); 
  incr[24] = wl[3]*((-1.5*fl[24])-0.8660254037844386*fl[15])+dxvl[3]*((-0.4330127018922193*fl[13])-0.25*fl[5]); 
  incr[25] = wl[3]*(0.8660254037844386*fl[30]+0.5*fl[25])+dxvl[3]*(0.25*fl[22]+0.1443375672974065*fl[14]); 
  incr[26] = wl[3]*((-1.5*fl[26])-0.8660254037844386*fl[18])+dxvl[3]*((-0.4330127018922193*fl[16])-0.25*fl[7]); 
  incr[27] = dxvl[3]*((-0.4330127018922193*fl[31])-0.25*fl[29])+wl[3]*((-1.5*fl[27])-0.8660254037844386*fl[21]); 
  incr[28] = wl[3]*((-1.5*fl[28])-0.8660254037844386*fl[23])+dxvl[3]*((-0.4330127018922193*fl[20])-0.25*fl[12]); 
  incr[29] = wl[3]*(0.8660254037844386*fl[31]+0.5*fl[29])+dxvl[3]*(0.25*fl[27]+0.1443375672974065*fl[21]); 
  incr[30] = wl[3]*((-1.5*fl[30])-0.8660254037844386*fl[25])+dxvl[3]*((-0.4330127018922193*fl[22])-0.25*fl[14]); 
  incr[31] = wl[3]*((-1.5*fl[31])-0.8660254037844386*fl[29])+dxvl[3]*((-0.4330127018922193*fl[27])-0.25*fl[21]); 

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
  outr[16] += incr[16]*rdxr2; 
  outr[17] += incr[17]*rdxr2; 
  outr[18] += incr[18]*rdxr2; 
  outr[19] += incr[19]*rdxr2; 
  outr[20] += incr[20]*rdxr2; 
  outr[21] += incr[21]*rdxr2; 
  outr[22] += incr[22]*rdxr2; 
  outr[23] += incr[23]*rdxr2; 
  outr[24] += incr[24]*rdxr2; 
  outr[25] += incr[25]*rdxr2; 
  outr[26] += incr[26]*rdxr2; 
  outr[27] += incr[27]*rdxr2; 
  outr[28] += incr[28]*rdxr2; 
  outr[29] += incr[29]*rdxr2; 
  outr[30] += incr[30]*rdxr2; 
  outr[31] += incr[31]*rdxr2; 

  outl[0] += -1.0*incr[0]*rdxl2; 
  outl[1] += -1.0*incr[1]*rdxl2; 
  outl[2] += incr[2]*rdxl2; 
  outl[3] += -1.0*incr[3]*rdxl2; 
  outl[4] += -1.0*incr[4]*rdxl2; 
  outl[5] += -1.0*incr[5]*rdxl2; 
  outl[6] += incr[6]*rdxl2; 
  outl[7] += -1.0*incr[7]*rdxl2; 
  outl[8] += incr[8]*rdxl2; 
  outl[9] += -1.0*incr[9]*rdxl2; 
  outl[10] += incr[10]*rdxl2; 
  outl[11] += -1.0*incr[11]*rdxl2; 
  outl[12] += -1.0*incr[12]*rdxl2; 
  outl[13] += incr[13]*rdxl2; 
  outl[14] += -1.0*incr[14]*rdxl2; 
  outl[15] += -1.0*incr[15]*rdxl2; 
  outl[16] += incr[16]*rdxl2; 
  outl[17] += incr[17]*rdxl2; 
  outl[18] += -1.0*incr[18]*rdxl2; 
  outl[19] += incr[19]*rdxl2; 
  outl[20] += incr[20]*rdxl2; 
  outl[21] += -1.0*incr[21]*rdxl2; 
  outl[22] += incr[22]*rdxl2; 
  outl[23] += -1.0*incr[23]*rdxl2; 
  outl[24] += incr[24]*rdxl2; 
  outl[25] += -1.0*incr[25]*rdxl2; 
  outl[26] += incr[26]*rdxl2; 
  outl[27] += incr[27]*rdxl2; 
  outl[28] += incr[28]*rdxl2; 
  outl[29] += -1.0*incr[29]*rdxl2; 
  outl[30] += incr[30]*rdxl2; 
  outl[31] += incr[31]*rdxl2; 
  } else { 
  incr[0] = dxvr[3]*(0.1443375672974065*fr[4]-0.25*fr[10])+(0.5*fr[0]-0.8660254037844386*fr[2])*wr[3]; 
  incr[1] = dxvr[3]*(0.1443375672974065*fr[9]-0.25*fr[17])+wr[3]*(0.5*fr[1]-0.8660254037844386*fr[6]); 
  incr[2] = dxvr[3]*(0.4330127018922193*fr[10]-0.25*fr[4])+(1.5*fr[2]-0.8660254037844386*fr[0])*wr[3]; 
  incr[3] = dxvr[3]*(0.1443375672974065*fr[11]-0.25*fr[19])+wr[3]*(0.5*fr[3]-0.8660254037844386*fr[8]); 
  incr[4] = wr[3]*(0.5*fr[4]-0.8660254037844386*fr[10])+(0.1443375672974065*fr[0]-0.25*fr[2])*dxvr[3]; 
  incr[5] = dxvr[3]*(0.1443375672974065*fr[15]-0.25*fr[24])+wr[3]*(0.5*fr[5]-0.8660254037844386*fr[13]); 
  incr[6] = dxvr[3]*(0.4330127018922193*fr[17]-0.25*fr[9])+wr[3]*(1.5*fr[6]-0.8660254037844386*fr[1]); 
  incr[7] = dxvr[3]*(0.1443375672974065*fr[18]-0.25*fr[26])+wr[3]*(0.5*fr[7]-0.8660254037844386*fr[16]); 
  incr[8] = dxvr[3]*(0.4330127018922193*fr[19]-0.25*fr[11])+wr[3]*(1.5*fr[8]-0.8660254037844386*fr[3]); 
  incr[9] = wr[3]*(0.5*fr[9]-0.8660254037844386*fr[17])+dxvr[3]*(0.1443375672974065*fr[1]-0.25*fr[6]); 
  incr[10] = wr[3]*(1.5*fr[10]-0.8660254037844386*fr[4])+(0.4330127018922193*fr[2]-0.25*fr[0])*dxvr[3]; 
  incr[11] = wr[3]*(0.5*fr[11]-0.8660254037844386*fr[19])+dxvr[3]*(0.1443375672974065*fr[3]-0.25*fr[8]); 
  incr[12] = dxvr[3]*(0.1443375672974065*fr[23]-0.25*fr[28])+wr[3]*(0.5*fr[12]-0.8660254037844386*fr[20]); 
  incr[13] = dxvr[3]*(0.4330127018922193*fr[24]-0.25*fr[15])+wr[3]*(1.5*fr[13]-0.8660254037844386*fr[5]); 
  incr[14] = dxvr[3]*(0.1443375672974065*fr[25]-0.25*fr[30])+wr[3]*(0.5*fr[14]-0.8660254037844386*fr[22]); 
  incr[15] = wr[3]*(0.5*fr[15]-0.8660254037844386*fr[24])+dxvr[3]*(0.1443375672974065*fr[5]-0.25*fr[13]); 
  incr[16] = dxvr[3]*(0.4330127018922193*fr[26]-0.25*fr[18])+wr[3]*(1.5*fr[16]-0.8660254037844386*fr[7]); 
  incr[17] = wr[3]*(1.5*fr[17]-0.8660254037844386*fr[9])+dxvr[3]*(0.4330127018922193*fr[6]-0.25*fr[1]); 
  incr[18] = wr[3]*(0.5*fr[18]-0.8660254037844386*fr[26])+dxvr[3]*(0.1443375672974065*fr[7]-0.25*fr[16]); 
  incr[19] = wr[3]*(1.5*fr[19]-0.8660254037844386*fr[11])+dxvr[3]*(0.4330127018922193*fr[8]-0.25*fr[3]); 
  incr[20] = dxvr[3]*(0.4330127018922193*fr[28]-0.25*fr[23])+wr[3]*(1.5*fr[20]-0.8660254037844386*fr[12]); 
  incr[21] = dxvr[3]*(0.1443375672974065*fr[29]-0.25*fr[31])+wr[3]*(0.5*fr[21]-0.8660254037844386*fr[27]); 
  incr[22] = dxvr[3]*(0.4330127018922193*fr[30]-0.25*fr[25])+wr[3]*(1.5*fr[22]-0.8660254037844386*fr[14]); 
  incr[23] = wr[3]*(0.5*fr[23]-0.8660254037844386*fr[28])+dxvr[3]*(0.1443375672974065*fr[12]-0.25*fr[20]); 
  incr[24] = wr[3]*(1.5*fr[24]-0.8660254037844386*fr[15])+dxvr[3]*(0.4330127018922193*fr[13]-0.25*fr[5]); 
  incr[25] = wr[3]*(0.5*fr[25]-0.8660254037844386*fr[30])+dxvr[3]*(0.1443375672974065*fr[14]-0.25*fr[22]); 
  incr[26] = wr[3]*(1.5*fr[26]-0.8660254037844386*fr[18])+dxvr[3]*(0.4330127018922193*fr[16]-0.25*fr[7]); 
  incr[27] = dxvr[3]*(0.4330127018922193*fr[31]-0.25*fr[29])+wr[3]*(1.5*fr[27]-0.8660254037844386*fr[21]); 
  incr[28] = wr[3]*(1.5*fr[28]-0.8660254037844386*fr[23])+dxvr[3]*(0.4330127018922193*fr[20]-0.25*fr[12]); 
  incr[29] = wr[3]*(0.5*fr[29]-0.8660254037844386*fr[31])+dxvr[3]*(0.1443375672974065*fr[21]-0.25*fr[27]); 
  incr[30] = wr[3]*(1.5*fr[30]-0.8660254037844386*fr[25])+dxvr[3]*(0.4330127018922193*fr[22]-0.25*fr[14]); 
  incr[31] = wr[3]*(1.5*fr[31]-0.8660254037844386*fr[29])+dxvr[3]*(0.4330127018922193*fr[27]-0.25*fr[21]); 

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
  outr[16] += incr[16]*rdxr2; 
  outr[17] += incr[17]*rdxr2; 
  outr[18] += incr[18]*rdxr2; 
  outr[19] += incr[19]*rdxr2; 
  outr[20] += incr[20]*rdxr2; 
  outr[21] += incr[21]*rdxr2; 
  outr[22] += incr[22]*rdxr2; 
  outr[23] += incr[23]*rdxr2; 
  outr[24] += incr[24]*rdxr2; 
  outr[25] += incr[25]*rdxr2; 
  outr[26] += incr[26]*rdxr2; 
  outr[27] += incr[27]*rdxr2; 
  outr[28] += incr[28]*rdxr2; 
  outr[29] += incr[29]*rdxr2; 
  outr[30] += incr[30]*rdxr2; 
  outr[31] += incr[31]*rdxr2; 

  outl[0] += -1.0*incr[0]*rdxl2; 
  outl[1] += -1.0*incr[1]*rdxl2; 
  outl[2] += incr[2]*rdxl2; 
  outl[3] += -1.0*incr[3]*rdxl2; 
  outl[4] += -1.0*incr[4]*rdxl2; 
  outl[5] += -1.0*incr[5]*rdxl2; 
  outl[6] += incr[6]*rdxl2; 
  outl[7] += -1.0*incr[7]*rdxl2; 
  outl[8] += incr[8]*rdxl2; 
  outl[9] += -1.0*incr[9]*rdxl2; 
  outl[10] += incr[10]*rdxl2; 
  outl[11] += -1.0*incr[11]*rdxl2; 
  outl[12] += -1.0*incr[12]*rdxl2; 
  outl[13] += incr[13]*rdxl2; 
  outl[14] += -1.0*incr[14]*rdxl2; 
  outl[15] += -1.0*incr[15]*rdxl2; 
  outl[16] += incr[16]*rdxl2; 
  outl[17] += incr[17]*rdxl2; 
  outl[18] += -1.0*incr[18]*rdxl2; 
  outl[19] += incr[19]*rdxl2; 
  outl[20] += incr[20]*rdxl2; 
  outl[21] += -1.0*incr[21]*rdxl2; 
  outl[22] += incr[22]*rdxl2; 
  outl[23] += -1.0*incr[23]*rdxl2; 
  outl[24] += incr[24]*rdxl2; 
  outl[25] += -1.0*incr[25]*rdxl2; 
  outl[26] += incr[26]*rdxl2; 
  outl[27] += incr[27]*rdxl2; 
  outl[28] += incr[28]*rdxl2; 
  outl[29] += -1.0*incr[29]*rdxl2; 
  outl[30] += incr[30]*rdxl2; 
  outl[31] += incr[31]*rdxl2; 
  } 
} 
gkyl_real vlasov_surfvx_2x3v_ser_p1(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real amax, const gkyl_real *qmem, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // amax:      amax in global lax flux.
  // qmem:      q/m*EM fields.
  // fl/fr:     Distribution function in left/right cells 
  // outl/outr: output distribution function in left/right cells 
  // returns abs(amid) for use in determining amax in cfl and global lax flux 
  gkyl_real dv10l = 2/dxvl[2]; 
  gkyl_real dv10r = 2/dxvr[2]; 

  const gkyl_real dv1 = dxvr[2], wv1 = wr[2]; 
  const gkyl_real dv2 = dxvr[3], wv2 = wr[3]; 
  const gkyl_real dv3 = dxvr[4], wv3 = wr[4]; 
  const gkyl_real *E0 = &qmem[0]; 
  const gkyl_real *B0 = &qmem[12]; 
  const gkyl_real *B1 = &qmem[16]; 
  const gkyl_real *B2 = &qmem[20]; 

  gkyl_real Ghat[16]; 
  gkyl_real favg[16]; 
  gkyl_real alpha[16]; 

  favg[0] = (-1.224744871391589*fr[3])+1.224744871391589*fl[3]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0]; 
  favg[1] = (-1.224744871391589*fr[7])+1.224744871391589*fl[7]+0.7071067811865475*fr[1]+0.7071067811865475*fl[1]; 
  favg[2] = (-1.224744871391589*fr[8])+1.224744871391589*fl[8]+0.7071067811865475*fr[2]+0.7071067811865475*fl[2]; 
  favg[3] = (-1.224744871391589*fr[11])+1.224744871391589*fl[11]+0.7071067811865475*fr[4]+0.7071067811865475*fl[4]; 
  favg[4] = (-1.224744871391589*fr[14])+1.224744871391589*fl[14]+0.7071067811865475*fr[5]+0.7071067811865475*fl[5]; 
  favg[5] = (-1.224744871391589*fr[16])+1.224744871391589*fl[16]+0.7071067811865475*fr[6]+0.7071067811865475*fl[6]; 
  favg[6] = (-1.224744871391589*fr[18])+1.224744871391589*fl[18]+0.7071067811865475*fr[9]+0.7071067811865475*fl[9]; 
  favg[7] = (-1.224744871391589*fr[19])+1.224744871391589*fl[19]+0.7071067811865475*fr[10]+0.7071067811865475*fl[10]; 
  favg[8] = (-1.224744871391589*fr[21])+1.224744871391589*fl[21]+0.7071067811865475*fr[12]+0.7071067811865475*fl[12]; 
  favg[9] = (-1.224744871391589*fr[22])+1.224744871391589*fl[22]+0.7071067811865475*fr[13]+0.7071067811865475*fl[13]; 
  favg[10] = (-1.224744871391589*fr[25])+1.224744871391589*fl[25]+0.7071067811865475*fr[15]+0.7071067811865475*fl[15]; 
  favg[11] = (-1.224744871391589*fr[26])+1.224744871391589*fl[26]+0.7071067811865475*fr[17]+0.7071067811865475*fl[17]; 
  favg[12] = (-1.224744871391589*fr[27])+1.224744871391589*fl[27]+0.7071067811865475*fr[20]+0.7071067811865475*fl[20]; 
  favg[13] = (-1.224744871391589*fr[29])+1.224744871391589*fl[29]+0.7071067811865475*fr[23]+0.7071067811865475*fl[23]; 
  favg[14] = (-1.224744871391589*fr[30])+1.224744871391589*fl[30]+0.7071067811865475*fr[24]+0.7071067811865475*fl[24]; 
  favg[15] = (-1.224744871391589*fr[31])+1.224744871391589*fl[31]+0.7071067811865475*fr[28]+0.7071067811865475*fl[28]; 

  alpha[0] = 2.0*(B2[0]*wv2+E0[0])-2.0*B1[0]*wv3; 
  alpha[1] = 2.0*(B2[1]*wv2+E0[1])-2.0*B1[1]*wv3; 
  alpha[2] = 2.0*(B2[2]*wv2+E0[2])-2.0*B1[2]*wv3; 
  alpha[3] = 0.5773502691896258*B2[0]*dv2; 
  alpha[4] = -0.5773502691896258*B1[0]*dv3; 
  alpha[5] = 2.0*(B2[3]*wv2+E0[3])-2.0*B1[3]*wv3; 
  alpha[6] = 0.5773502691896258*B2[1]*dv2; 
  alpha[7] = 0.5773502691896258*B2[2]*dv2; 
  alpha[8] = -0.5773502691896258*B1[1]*dv3; 
  alpha[9] = -0.5773502691896258*B1[2]*dv3; 
  alpha[11] = 0.5773502691896258*B2[3]*dv2; 
  alpha[12] = -0.5773502691896258*B1[3]*dv3; 

  gkyl_real amid = 0.25*alpha[0]; 

  Ghat[0] = 0.3535533905932737*(1.732050807568877*(fr[3]+fl[3])-1.0*fr[0]+fl[0])*amax+0.125*(alpha[12]*favg[12]+alpha[11]*favg[11]+alpha[9]*favg[9]+alpha[8]*favg[8]+alpha[7]*favg[7]+alpha[6]*favg[6]+alpha[5]*favg[5]+alpha[4]*favg[4]+alpha[3]*favg[3]+alpha[2]*favg[2]+alpha[1]*favg[1]+alpha[0]*favg[0]); 
  Ghat[1] = 0.3535533905932737*(1.732050807568877*(fr[7]+fl[7])-1.0*fr[1]+fl[1])*amax+0.125*(alpha[9]*favg[12]+favg[9]*alpha[12]+alpha[7]*favg[11]+favg[7]*alpha[11]+alpha[4]*favg[8]+favg[4]*alpha[8]+alpha[3]*favg[6]+favg[3]*alpha[6]+alpha[2]*favg[5]+favg[2]*alpha[5]+alpha[0]*favg[1]+favg[0]*alpha[1]); 
  Ghat[2] = 0.3535533905932737*(1.732050807568877*(fr[8]+fl[8])-1.0*fr[2]+fl[2])*amax+0.125*(alpha[8]*favg[12]+favg[8]*alpha[12]+alpha[6]*favg[11]+favg[6]*alpha[11]+alpha[4]*favg[9]+favg[4]*alpha[9]+alpha[3]*favg[7]+favg[3]*alpha[7]+alpha[1]*favg[5]+favg[1]*alpha[5]+alpha[0]*favg[2]+favg[0]*alpha[2]); 
  Ghat[3] = 0.3535533905932737*(1.732050807568877*(fr[11]+fl[11])-1.0*fr[4]+fl[4])*amax+0.125*(alpha[12]*favg[15]+alpha[9]*favg[14]+alpha[8]*favg[13]+alpha[5]*favg[11]+favg[5]*alpha[11]+alpha[4]*favg[10]+alpha[2]*favg[7]+favg[2]*alpha[7]+alpha[1]*favg[6]+favg[1]*alpha[6]+alpha[0]*favg[3]+favg[0]*alpha[3]); 
  Ghat[4] = 0.3535533905932737*(1.732050807568877*(fr[14]+fl[14])-1.0*fr[5]+fl[5])*amax+0.125*(alpha[11]*favg[15]+alpha[7]*favg[14]+alpha[6]*favg[13]+alpha[5]*favg[12]+favg[5]*alpha[12]+alpha[3]*favg[10]+alpha[2]*favg[9]+favg[2]*alpha[9]+alpha[1]*favg[8]+favg[1]*alpha[8]+alpha[0]*favg[4]+favg[0]*alpha[4]); 
  Ghat[5] = 0.3535533905932737*(1.732050807568877*(fr[16]+fl[16])-1.0*fr[6]+fl[6])*amax+0.125*(alpha[4]*favg[12]+favg[4]*alpha[12]+alpha[3]*favg[11]+favg[3]*alpha[11]+alpha[8]*favg[9]+favg[8]*alpha[9]+alpha[6]*favg[7]+favg[6]*alpha[7]+alpha[0]*favg[5]+favg[0]*alpha[5]+alpha[1]*favg[2]+favg[1]*alpha[2]); 
  Ghat[6] = 0.3535533905932737*(1.732050807568877*(fr[18]+fl[18])-1.0*fr[9]+fl[9])*amax+0.125*(alpha[9]*favg[15]+alpha[12]*favg[14]+alpha[4]*favg[13]+alpha[2]*favg[11]+favg[2]*alpha[11]+alpha[8]*favg[10]+alpha[5]*favg[7]+favg[5]*alpha[7]+alpha[0]*favg[6]+favg[0]*alpha[6]+alpha[1]*favg[3]+favg[1]*alpha[3]); 
  Ghat[7] = 0.3535533905932737*(1.732050807568877*(fr[19]+fl[19])-1.0*fr[10]+fl[10])*amax+0.125*(alpha[8]*favg[15]+alpha[4]*favg[14]+alpha[12]*favg[13]+alpha[1]*favg[11]+favg[1]*alpha[11]+alpha[9]*favg[10]+alpha[0]*favg[7]+favg[0]*alpha[7]+alpha[5]*favg[6]+favg[5]*alpha[6]+alpha[2]*favg[3]+favg[2]*alpha[3]); 
  Ghat[8] = 0.3535533905932737*(1.732050807568877*(fr[21]+fl[21])-1.0*fr[12]+fl[12])*amax+0.125*(alpha[7]*favg[15]+alpha[11]*favg[14]+alpha[3]*favg[13]+alpha[2]*favg[12]+favg[2]*alpha[12]+alpha[6]*favg[10]+alpha[5]*favg[9]+favg[5]*alpha[9]+alpha[0]*favg[8]+favg[0]*alpha[8]+alpha[1]*favg[4]+favg[1]*alpha[4]); 
  Ghat[9] = 0.3535533905932737*(1.732050807568877*(fr[22]+fl[22])-1.0*fr[13]+fl[13])*amax+0.125*(alpha[6]*favg[15]+alpha[3]*favg[14]+alpha[11]*favg[13]+alpha[1]*favg[12]+favg[1]*alpha[12]+alpha[7]*favg[10]+alpha[0]*favg[9]+favg[0]*alpha[9]+alpha[5]*favg[8]+favg[5]*alpha[8]+alpha[2]*favg[4]+favg[2]*alpha[4]); 
  Ghat[10] = 0.3535533905932737*(1.732050807568877*(fr[25]+fl[25])-1.0*fr[15]+fl[15])*amax+0.125*(alpha[5]*favg[15]+alpha[2]*favg[14]+alpha[1]*favg[13]+alpha[11]*favg[12]+favg[11]*alpha[12]+alpha[0]*favg[10]+alpha[7]*favg[9]+favg[7]*alpha[9]+alpha[6]*favg[8]+favg[6]*alpha[8]+alpha[3]*favg[4]+favg[3]*alpha[4]); 
  Ghat[11] = 0.3535533905932737*(1.732050807568877*(fr[26]+fl[26])-1.0*fr[17]+fl[17])*amax+0.125*(alpha[4]*favg[15]+alpha[8]*favg[14]+alpha[9]*favg[13]+favg[10]*alpha[12]+alpha[0]*favg[11]+favg[0]*alpha[11]+alpha[1]*favg[7]+favg[1]*alpha[7]+alpha[2]*favg[6]+favg[2]*alpha[6]+alpha[3]*favg[5]+favg[3]*alpha[5]); 
  Ghat[12] = 0.3535533905932737*(1.732050807568877*(fr[27]+fl[27])-1.0*fr[20]+fl[20])*amax+0.125*(alpha[3]*favg[15]+alpha[6]*favg[14]+alpha[7]*favg[13]+alpha[0]*favg[12]+favg[0]*alpha[12]+favg[10]*alpha[11]+alpha[1]*favg[9]+favg[1]*alpha[9]+alpha[2]*favg[8]+favg[2]*alpha[8]+alpha[4]*favg[5]+favg[4]*alpha[5]); 
  Ghat[13] = 0.3535533905932737*(1.732050807568877*(fr[29]+fl[29])-1.0*fr[23]+fl[23])*amax+0.125*(alpha[2]*favg[15]+alpha[5]*favg[14]+alpha[0]*favg[13]+alpha[7]*favg[12]+favg[7]*alpha[12]+alpha[9]*favg[11]+favg[9]*alpha[11]+alpha[1]*favg[10]+alpha[3]*favg[8]+favg[3]*alpha[8]+alpha[4]*favg[6]+favg[4]*alpha[6]); 
  Ghat[14] = 0.3535533905932737*(1.732050807568877*(fr[30]+fl[30])-1.0*fr[24]+fl[24])*amax+0.125*(alpha[1]*favg[15]+alpha[0]*favg[14]+alpha[5]*favg[13]+alpha[6]*favg[12]+favg[6]*alpha[12]+alpha[8]*favg[11]+favg[8]*alpha[11]+alpha[2]*favg[10]+alpha[3]*favg[9]+favg[3]*alpha[9]+alpha[4]*favg[7]+favg[4]*alpha[7]); 
  Ghat[15] = 0.3535533905932737*(1.732050807568877*(fr[31]+fl[31])-1.0*fr[28]+fl[28])*amax+0.125*(alpha[0]*favg[15]+alpha[1]*favg[14]+alpha[2]*favg[13]+alpha[3]*favg[12]+favg[3]*alpha[12]+alpha[4]*favg[11]+favg[4]*alpha[11]+alpha[5]*favg[10]+alpha[6]*favg[9]+favg[6]*alpha[9]+alpha[7]*favg[8]+favg[7]*alpha[8]); 

  outr[0] += 0.7071067811865475*Ghat[0]*dv10r; 
  outr[1] += 0.7071067811865475*Ghat[1]*dv10r; 
  outr[2] += 0.7071067811865475*Ghat[2]*dv10r; 
  outr[3] += -1.224744871391589*Ghat[0]*dv10r; 
  outr[4] += 0.7071067811865475*Ghat[3]*dv10r; 
  outr[5] += 0.7071067811865475*Ghat[4]*dv10r; 
  outr[6] += 0.7071067811865475*Ghat[5]*dv10r; 
  outr[7] += -1.224744871391589*Ghat[1]*dv10r; 
  outr[8] += -1.224744871391589*Ghat[2]*dv10r; 
  outr[9] += 0.7071067811865475*Ghat[6]*dv10r; 
  outr[10] += 0.7071067811865475*Ghat[7]*dv10r; 
  outr[11] += -1.224744871391589*Ghat[3]*dv10r; 
  outr[12] += 0.7071067811865475*Ghat[8]*dv10r; 
  outr[13] += 0.7071067811865475*Ghat[9]*dv10r; 
  outr[14] += -1.224744871391589*Ghat[4]*dv10r; 
  outr[15] += 0.7071067811865475*Ghat[10]*dv10r; 
  outr[16] += -1.224744871391589*Ghat[5]*dv10r; 
  outr[17] += 0.7071067811865475*Ghat[11]*dv10r; 
  outr[18] += -1.224744871391589*Ghat[6]*dv10r; 
  outr[19] += -1.224744871391589*Ghat[7]*dv10r; 
  outr[20] += 0.7071067811865475*Ghat[12]*dv10r; 
  outr[21] += -1.224744871391589*Ghat[8]*dv10r; 
  outr[22] += -1.224744871391589*Ghat[9]*dv10r; 
  outr[23] += 0.7071067811865475*Ghat[13]*dv10r; 
  outr[24] += 0.7071067811865475*Ghat[14]*dv10r; 
  outr[25] += -1.224744871391589*Ghat[10]*dv10r; 
  outr[26] += -1.224744871391589*Ghat[11]*dv10r; 
  outr[27] += -1.224744871391589*Ghat[12]*dv10r; 
  outr[28] += 0.7071067811865475*Ghat[15]*dv10r; 
  outr[29] += -1.224744871391589*Ghat[13]*dv10r; 
  outr[30] += -1.224744871391589*Ghat[14]*dv10r; 
  outr[31] += -1.224744871391589*Ghat[15]*dv10r; 

  outl[0] += -0.7071067811865475*Ghat[0]*dv10l; 
  outl[1] += -0.7071067811865475*Ghat[1]*dv10l; 
  outl[2] += -0.7071067811865475*Ghat[2]*dv10l; 
  outl[3] += -1.224744871391589*Ghat[0]*dv10l; 
  outl[4] += -0.7071067811865475*Ghat[3]*dv10l; 
  outl[5] += -0.7071067811865475*Ghat[4]*dv10l; 
  outl[6] += -0.7071067811865475*Ghat[5]*dv10l; 
  outl[7] += -1.224744871391589*Ghat[1]*dv10l; 
  outl[8] += -1.224744871391589*Ghat[2]*dv10l; 
  outl[9] += -0.7071067811865475*Ghat[6]*dv10l; 
  outl[10] += -0.7071067811865475*Ghat[7]*dv10l; 
  outl[11] += -1.224744871391589*Ghat[3]*dv10l; 
  outl[12] += -0.7071067811865475*Ghat[8]*dv10l; 
  outl[13] += -0.7071067811865475*Ghat[9]*dv10l; 
  outl[14] += -1.224744871391589*Ghat[4]*dv10l; 
  outl[15] += -0.7071067811865475*Ghat[10]*dv10l; 
  outl[16] += -1.224744871391589*Ghat[5]*dv10l; 
  outl[17] += -0.7071067811865475*Ghat[11]*dv10l; 
  outl[18] += -1.224744871391589*Ghat[6]*dv10l; 
  outl[19] += -1.224744871391589*Ghat[7]*dv10l; 
  outl[20] += -0.7071067811865475*Ghat[12]*dv10l; 
  outl[21] += -1.224744871391589*Ghat[8]*dv10l; 
  outl[22] += -1.224744871391589*Ghat[9]*dv10l; 
  outl[23] += -0.7071067811865475*Ghat[13]*dv10l; 
  outl[24] += -0.7071067811865475*Ghat[14]*dv10l; 
  outl[25] += -1.224744871391589*Ghat[10]*dv10l; 
  outl[26] += -1.224744871391589*Ghat[11]*dv10l; 
  outl[27] += -1.224744871391589*Ghat[12]*dv10l; 
  outl[28] += -0.7071067811865475*Ghat[15]*dv10l; 
  outl[29] += -1.224744871391589*Ghat[13]*dv10l; 
  outl[30] += -1.224744871391589*Ghat[14]*dv10l; 
  outl[31] += -1.224744871391589*Ghat[15]*dv10l; 

  return fabs(amid); 
} 
gkyl_real vlasov_surfvy_2x3v_ser_p1(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real amax, const gkyl_real *qmem, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // amax:      amax in global lax flux.
  // qmem:      q/m*EM fields.
  // fl/fr:     Distribution function in left/right cells 
  // outl/outr: output distribution function in left/right cells 
  // returns abs(amid) for use in determining amax in cfl and global lax flux 
  gkyl_real dv11l = 2/dxvl[3]; 
  gkyl_real dv11r = 2/dxvr[3]; 

  const gkyl_real dv1 = dxvr[2], wv1 = wr[2]; 
  const gkyl_real dv2 = dxvr[3], wv2 = wr[3]; 
  const gkyl_real dv3 = dxvr[4], wv3 = wr[4]; 
  const gkyl_real *E1 = &qmem[4]; 
  const gkyl_real *B0 = &qmem[12]; 
  const gkyl_real *B1 = &qmem[16]; 
  const gkyl_real *B2 = &qmem[20]; 

  gkyl_real Ghat[16]; 
  gkyl_real favg[16]; 
  gkyl_real alpha[16]; 

  favg[0] = (-1.224744871391589*fr[4])+1.224744871391589*fl[4]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0]; 
  favg[1] = (-1.224744871391589*fr[9])+1.224744871391589*fl[9]+0.7071067811865475*fr[1]+0.7071067811865475*fl[1]; 
  favg[2] = (-1.224744871391589*fr[10])+1.224744871391589*fl[10]+0.7071067811865475*fr[2]+0.7071067811865475*fl[2]; 
  favg[3] = (-1.224744871391589*fr[11])+1.224744871391589*fl[11]+0.7071067811865475*fr[3]+0.7071067811865475*fl[3]; 
  favg[4] = (-1.224744871391589*fr[15])+1.224744871391589*fl[15]+0.7071067811865475*fr[5]+0.7071067811865475*fl[5]; 
  favg[5] = (-1.224744871391589*fr[17])+1.224744871391589*fl[17]+0.7071067811865475*fr[6]+0.7071067811865475*fl[6]; 
  favg[6] = (-1.224744871391589*fr[18])+1.224744871391589*fl[18]+0.7071067811865475*fr[7]+0.7071067811865475*fl[7]; 
  favg[7] = (-1.224744871391589*fr[19])+1.224744871391589*fl[19]+0.7071067811865475*fr[8]+0.7071067811865475*fl[8]; 
  favg[8] = (-1.224744871391589*fr[23])+1.224744871391589*fl[23]+0.7071067811865475*fr[12]+0.7071067811865475*fl[12]; 
  favg[9] = (-1.224744871391589*fr[24])+1.224744871391589*fl[24]+0.7071067811865475*fr[13]+0.7071067811865475*fl[13]; 
  favg[10] = (-1.224744871391589*fr[25])+1.224744871391589*fl[25]+0.7071067811865475*fr[14]+0.7071067811865475*fl[14]; 
  favg[11] = (-1.224744871391589*fr[26])+1.224744871391589*fl[26]+0.7071067811865475*fr[16]+0.7071067811865475*fl[16]; 
  favg[12] = (-1.224744871391589*fr[28])+1.224744871391589*fl[28]+0.7071067811865475*fr[20]+0.7071067811865475*fl[20]; 
  favg[13] = (-1.224744871391589*fr[29])+1.224744871391589*fl[29]+0.7071067811865475*fr[21]+0.7071067811865475*fl[21]; 
  favg[14] = (-1.224744871391589*fr[30])+1.224744871391589*fl[30]+0.7071067811865475*fr[22]+0.7071067811865475*fl[22]; 
  favg[15] = (-1.224744871391589*fr[31])+1.224744871391589*fl[31]+0.7071067811865475*fr[27]+0.7071067811865475*fl[27]; 

  alpha[0] = 2.0*B0[0]*wv3-2.0*B2[0]*wv1+2.0*E1[0]; 
  alpha[1] = 2.0*B0[1]*wv3-2.0*B2[1]*wv1+2.0*E1[1]; 
  alpha[2] = 2.0*B0[2]*wv3-2.0*B2[2]*wv1+2.0*E1[2]; 
  alpha[3] = -0.5773502691896258*B2[0]*dv1; 
  alpha[4] = 0.5773502691896258*B0[0]*dv3; 
  alpha[5] = 2.0*B0[3]*wv3-2.0*B2[3]*wv1+2.0*E1[3]; 
  alpha[6] = -0.5773502691896258*B2[1]*dv1; 
  alpha[7] = -0.5773502691896258*B2[2]*dv1; 
  alpha[8] = 0.5773502691896258*B0[1]*dv3; 
  alpha[9] = 0.5773502691896258*B0[2]*dv3; 
  alpha[11] = -0.5773502691896258*B2[3]*dv1; 
  alpha[12] = 0.5773502691896258*B0[3]*dv3; 

  gkyl_real amid = 0.25*alpha[0]; 

  Ghat[0] = 0.3535533905932737*(1.732050807568877*(fr[4]+fl[4])-1.0*fr[0]+fl[0])*amax+0.125*(alpha[12]*favg[12]+alpha[11]*favg[11]+alpha[9]*favg[9]+alpha[8]*favg[8]+alpha[7]*favg[7]+alpha[6]*favg[6]+alpha[5]*favg[5]+alpha[4]*favg[4]+alpha[3]*favg[3]+alpha[2]*favg[2]+alpha[1]*favg[1]+alpha[0]*favg[0]); 
  Ghat[1] = 0.3535533905932737*(1.732050807568877*(fr[9]+fl[9])-1.0*fr[1]+fl[1])*amax+0.125*(alpha[9]*favg[12]+favg[9]*alpha[12]+alpha[7]*favg[11]+favg[7]*alpha[11]+alpha[4]*favg[8]+favg[4]*alpha[8]+alpha[3]*favg[6]+favg[3]*alpha[6]+alpha[2]*favg[5]+favg[2]*alpha[5]+alpha[0]*favg[1]+favg[0]*alpha[1]); 
  Ghat[2] = 0.3535533905932737*(1.732050807568877*(fr[10]+fl[10])-1.0*fr[2]+fl[2])*amax+0.125*(alpha[8]*favg[12]+favg[8]*alpha[12]+alpha[6]*favg[11]+favg[6]*alpha[11]+alpha[4]*favg[9]+favg[4]*alpha[9]+alpha[3]*favg[7]+favg[3]*alpha[7]+alpha[1]*favg[5]+favg[1]*alpha[5]+alpha[0]*favg[2]+favg[0]*alpha[2]); 
  Ghat[3] = 0.3535533905932737*(1.732050807568877*(fr[11]+fl[11])-1.0*fr[3]+fl[3])*amax+0.125*(alpha[12]*favg[15]+alpha[9]*favg[14]+alpha[8]*favg[13]+alpha[5]*favg[11]+favg[5]*alpha[11]+alpha[4]*favg[10]+alpha[2]*favg[7]+favg[2]*alpha[7]+alpha[1]*favg[6]+favg[1]*alpha[6]+alpha[0]*favg[3]+favg[0]*alpha[3]); 
  Ghat[4] = 0.3535533905932737*(1.732050807568877*(fr[15]+fl[15])-1.0*fr[5]+fl[5])*amax+0.125*(alpha[11]*favg[15]+alpha[7]*favg[14]+alpha[6]*favg[13]+alpha[5]*favg[12]+favg[5]*alpha[12]+alpha[3]*favg[10]+alpha[2]*favg[9]+favg[2]*alpha[9]+alpha[1]*favg[8]+favg[1]*alpha[8]+alpha[0]*favg[4]+favg[0]*alpha[4]); 
  Ghat[5] = 0.3535533905932737*(1.732050807568877*(fr[17]+fl[17])-1.0*fr[6]+fl[6])*amax+0.125*(alpha[4]*favg[12]+favg[4]*alpha[12]+alpha[3]*favg[11]+favg[3]*alpha[11]+alpha[8]*favg[9]+favg[8]*alpha[9]+alpha[6]*favg[7]+favg[6]*alpha[7]+alpha[0]*favg[5]+favg[0]*alpha[5]+alpha[1]*favg[2]+favg[1]*alpha[2]); 
  Ghat[6] = 0.3535533905932737*(1.732050807568877*(fr[18]+fl[18])-1.0*fr[7]+fl[7])*amax+0.125*(alpha[9]*favg[15]+alpha[12]*favg[14]+alpha[4]*favg[13]+alpha[2]*favg[11]+favg[2]*alpha[11]+alpha[8]*favg[10]+alpha[5]*favg[7]+favg[5]*alpha[7]+alpha[0]*favg[6]+favg[0]*alpha[6]+alpha[1]*favg[3]+favg[1]*alpha[3]); 
  Ghat[7] = 0.3535533905932737*(1.732050807568877*(fr[19]+fl[19])-1.0*fr[8]+fl[8])*amax+0.125*(alpha[8]*favg[15]+alpha[4]*favg[14]+alpha[12]*favg[13]+alpha[1]*favg[11]+favg[1]*alpha[11]+alpha[9]*favg[10]+alpha[0]*favg[7]+favg[0]*alpha[7]+alpha[5]*favg[6]+favg[5]*alpha[6]+alpha[2]*favg[3]+favg[2]*alpha[3]); 
  Ghat[8] = 0.3535533905932737*(1.732050807568877*(fr[23]+fl[23])-1.0*fr[12]+fl[12])*amax+0.125*(alpha[7]*favg[15]+alpha[11]*favg[14]+alpha[3]*favg[13]+alpha[2]*favg[12]+favg[2]*alpha[12]+alpha[6]*favg[10]+alpha[5]*favg[9]+favg[5]*alpha[9]+alpha[0]*favg[8]+favg[0]*alpha[8]+alpha[1]*favg[4]+favg[1]*alpha[4]); 
  Ghat[9] = 0.3535533905932737*(1.732050807568877*(fr[24]+fl[24])-1.0*fr[13]+fl[13])*amax+0.125*(alpha[6]*favg[15]+alpha[3]*favg[14]+alpha[11]*favg[13]+alpha[1]*favg[12]+favg[1]*alpha[12]+alpha[7]*favg[10]+alpha[0]*favg[9]+favg[0]*alpha[9]+alpha[5]*favg[8]+favg[5]*alpha[8]+alpha[2]*favg[4]+favg[2]*alpha[4]); 
  Ghat[10] = 0.3535533905932737*(1.732050807568877*(fr[25]+fl[25])-1.0*fr[14]+fl[14])*amax+0.125*(alpha[5]*favg[15]+alpha[2]*favg[14]+alpha[1]*favg[13]+alpha[11]*favg[12]+favg[11]*alpha[12]+alpha[0]*favg[10]+alpha[7]*favg[9]+favg[7]*alpha[9]+alpha[6]*favg[8]+favg[6]*alpha[8]+alpha[3]*favg[4]+favg[3]*alpha[4]); 
  Ghat[11] = 0.3535533905932737*(1.732050807568877*(fr[26]+fl[26])-1.0*fr[16]+fl[16])*amax+0.125*(alpha[4]*favg[15]+alpha[8]*favg[14]+alpha[9]*favg[13]+favg[10]*alpha[12]+alpha[0]*favg[11]+favg[0]*alpha[11]+alpha[1]*favg[7]+favg[1]*alpha[7]+alpha[2]*favg[6]+favg[2]*alpha[6]+alpha[3]*favg[5]+favg[3]*alpha[5]); 
  Ghat[12] = 0.3535533905932737*(1.732050807568877*(fr[28]+fl[28])-1.0*fr[20]+fl[20])*amax+0.125*(alpha[3]*favg[15]+alpha[6]*favg[14]+alpha[7]*favg[13]+alpha[0]*favg[12]+favg[0]*alpha[12]+favg[10]*alpha[11]+alpha[1]*favg[9]+favg[1]*alpha[9]+alpha[2]*favg[8]+favg[2]*alpha[8]+alpha[4]*favg[5]+favg[4]*alpha[5]); 
  Ghat[13] = 0.3535533905932737*(1.732050807568877*(fr[29]+fl[29])-1.0*fr[21]+fl[21])*amax+0.125*(alpha[2]*favg[15]+alpha[5]*favg[14]+alpha[0]*favg[13]+alpha[7]*favg[12]+favg[7]*alpha[12]+alpha[9]*favg[11]+favg[9]*alpha[11]+alpha[1]*favg[10]+alpha[3]*favg[8]+favg[3]*alpha[8]+alpha[4]*favg[6]+favg[4]*alpha[6]); 
  Ghat[14] = 0.3535533905932737*(1.732050807568877*(fr[30]+fl[30])-1.0*fr[22]+fl[22])*amax+0.125*(alpha[1]*favg[15]+alpha[0]*favg[14]+alpha[5]*favg[13]+alpha[6]*favg[12]+favg[6]*alpha[12]+alpha[8]*favg[11]+favg[8]*alpha[11]+alpha[2]*favg[10]+alpha[3]*favg[9]+favg[3]*alpha[9]+alpha[4]*favg[7]+favg[4]*alpha[7]); 
  Ghat[15] = 0.3535533905932737*(1.732050807568877*(fr[31]+fl[31])-1.0*fr[27]+fl[27])*amax+0.125*(alpha[0]*favg[15]+alpha[1]*favg[14]+alpha[2]*favg[13]+alpha[3]*favg[12]+favg[3]*alpha[12]+alpha[4]*favg[11]+favg[4]*alpha[11]+alpha[5]*favg[10]+alpha[6]*favg[9]+favg[6]*alpha[9]+alpha[7]*favg[8]+favg[7]*alpha[8]); 

  outr[0] += 0.7071067811865475*Ghat[0]*dv11r; 
  outr[1] += 0.7071067811865475*Ghat[1]*dv11r; 
  outr[2] += 0.7071067811865475*Ghat[2]*dv11r; 
  outr[3] += 0.7071067811865475*Ghat[3]*dv11r; 
  outr[4] += -1.224744871391589*Ghat[0]*dv11r; 
  outr[5] += 0.7071067811865475*Ghat[4]*dv11r; 
  outr[6] += 0.7071067811865475*Ghat[5]*dv11r; 
  outr[7] += 0.7071067811865475*Ghat[6]*dv11r; 
  outr[8] += 0.7071067811865475*Ghat[7]*dv11r; 
  outr[9] += -1.224744871391589*Ghat[1]*dv11r; 
  outr[10] += -1.224744871391589*Ghat[2]*dv11r; 
  outr[11] += -1.224744871391589*Ghat[3]*dv11r; 
  outr[12] += 0.7071067811865475*Ghat[8]*dv11r; 
  outr[13] += 0.7071067811865475*Ghat[9]*dv11r; 
  outr[14] += 0.7071067811865475*Ghat[10]*dv11r; 
  outr[15] += -1.224744871391589*Ghat[4]*dv11r; 
  outr[16] += 0.7071067811865475*Ghat[11]*dv11r; 
  outr[17] += -1.224744871391589*Ghat[5]*dv11r; 
  outr[18] += -1.224744871391589*Ghat[6]*dv11r; 
  outr[19] += -1.224744871391589*Ghat[7]*dv11r; 
  outr[20] += 0.7071067811865475*Ghat[12]*dv11r; 
  outr[21] += 0.7071067811865475*Ghat[13]*dv11r; 
  outr[22] += 0.7071067811865475*Ghat[14]*dv11r; 
  outr[23] += -1.224744871391589*Ghat[8]*dv11r; 
  outr[24] += -1.224744871391589*Ghat[9]*dv11r; 
  outr[25] += -1.224744871391589*Ghat[10]*dv11r; 
  outr[26] += -1.224744871391589*Ghat[11]*dv11r; 
  outr[27] += 0.7071067811865475*Ghat[15]*dv11r; 
  outr[28] += -1.224744871391589*Ghat[12]*dv11r; 
  outr[29] += -1.224744871391589*Ghat[13]*dv11r; 
  outr[30] += -1.224744871391589*Ghat[14]*dv11r; 
  outr[31] += -1.224744871391589*Ghat[15]*dv11r; 

  outl[0] += -0.7071067811865475*Ghat[0]*dv11l; 
  outl[1] += -0.7071067811865475*Ghat[1]*dv11l; 
  outl[2] += -0.7071067811865475*Ghat[2]*dv11l; 
  outl[3] += -0.7071067811865475*Ghat[3]*dv11l; 
  outl[4] += -1.224744871391589*Ghat[0]*dv11l; 
  outl[5] += -0.7071067811865475*Ghat[4]*dv11l; 
  outl[6] += -0.7071067811865475*Ghat[5]*dv11l; 
  outl[7] += -0.7071067811865475*Ghat[6]*dv11l; 
  outl[8] += -0.7071067811865475*Ghat[7]*dv11l; 
  outl[9] += -1.224744871391589*Ghat[1]*dv11l; 
  outl[10] += -1.224744871391589*Ghat[2]*dv11l; 
  outl[11] += -1.224744871391589*Ghat[3]*dv11l; 
  outl[12] += -0.7071067811865475*Ghat[8]*dv11l; 
  outl[13] += -0.7071067811865475*Ghat[9]*dv11l; 
  outl[14] += -0.7071067811865475*Ghat[10]*dv11l; 
  outl[15] += -1.224744871391589*Ghat[4]*dv11l; 
  outl[16] += -0.7071067811865475*Ghat[11]*dv11l; 
  outl[17] += -1.224744871391589*Ghat[5]*dv11l; 
  outl[18] += -1.224744871391589*Ghat[6]*dv11l; 
  outl[19] += -1.224744871391589*Ghat[7]*dv11l; 
  outl[20] += -0.7071067811865475*Ghat[12]*dv11l; 
  outl[21] += -0.7071067811865475*Ghat[13]*dv11l; 
  outl[22] += -0.7071067811865475*Ghat[14]*dv11l; 
  outl[23] += -1.224744871391589*Ghat[8]*dv11l; 
  outl[24] += -1.224744871391589*Ghat[9]*dv11l; 
  outl[25] += -1.224744871391589*Ghat[10]*dv11l; 
  outl[26] += -1.224744871391589*Ghat[11]*dv11l; 
  outl[27] += -0.7071067811865475*Ghat[15]*dv11l; 
  outl[28] += -1.224744871391589*Ghat[12]*dv11l; 
  outl[29] += -1.224744871391589*Ghat[13]*dv11l; 
  outl[30] += -1.224744871391589*Ghat[14]*dv11l; 
  outl[31] += -1.224744871391589*Ghat[15]*dv11l; 

  return fabs(amid); 
} 
gkyl_real vlasov_surfvz_2x3v_ser_p1(const gkyl_real *wl, const gkyl_real *wr, const gkyl_real *dxvl, const gkyl_real *dxvr, const gkyl_real amax, const gkyl_real *qmem, const gkyl_real *fl, const gkyl_real *fr, gkyl_real* restrict outl, gkyl_real* restrict outr) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // amax:      amax in global lax flux.
  // qmem:      q/m*EM fields.
  // fl/fr:     Distribution function in left/right cells 
  // outl/outr: output distribution function in left/right cells 
  // returns abs(amid) for use in determining amax in cfl and global lax flux 
  gkyl_real dv12l = 2/dxvl[4]; 
  gkyl_real dv12r = 2/dxvr[4]; 

  const gkyl_real dv1 = dxvr[2], wv1 = wr[2]; 
  const gkyl_real dv2 = dxvr[3], wv2 = wr[3]; 
  const gkyl_real dv3 = dxvr[4], wv3 = wr[4]; 
  const gkyl_real *E2 = &qmem[8]; 
  const gkyl_real *B0 = &qmem[12]; 
  const gkyl_real *B1 = &qmem[16]; 
  const gkyl_real *B2 = &qmem[20]; 

  gkyl_real Ghat[16]; 
  gkyl_real favg[16]; 
  gkyl_real alpha[16]; 

  favg[0] = (-1.224744871391589*fr[5])+1.224744871391589*fl[5]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0]; 
  favg[1] = (-1.224744871391589*fr[12])+1.224744871391589*fl[12]+0.7071067811865475*fr[1]+0.7071067811865475*fl[1]; 
  favg[2] = (-1.224744871391589*fr[13])+1.224744871391589*fl[13]+0.7071067811865475*fr[2]+0.7071067811865475*fl[2]; 
  favg[3] = (-1.224744871391589*fr[14])+1.224744871391589*fl[14]+0.7071067811865475*fr[3]+0.7071067811865475*fl[3]; 
  favg[4] = (-1.224744871391589*fr[15])+1.224744871391589*fl[15]+0.7071067811865475*fr[4]+0.7071067811865475*fl[4]; 
  favg[5] = (-1.224744871391589*fr[20])+1.224744871391589*fl[20]+0.7071067811865475*fr[6]+0.7071067811865475*fl[6]; 
  favg[6] = (-1.224744871391589*fr[21])+1.224744871391589*fl[21]+0.7071067811865475*fr[7]+0.7071067811865475*fl[7]; 
  favg[7] = (-1.224744871391589*fr[22])+1.224744871391589*fl[22]+0.7071067811865475*fr[8]+0.7071067811865475*fl[8]; 
  favg[8] = (-1.224744871391589*fr[23])+1.224744871391589*fl[23]+0.7071067811865475*fr[9]+0.7071067811865475*fl[9]; 
  favg[9] = (-1.224744871391589*fr[24])+1.224744871391589*fl[24]+0.7071067811865475*fr[10]+0.7071067811865475*fl[10]; 
  favg[10] = (-1.224744871391589*fr[25])+1.224744871391589*fl[25]+0.7071067811865475*fr[11]+0.7071067811865475*fl[11]; 
  favg[11] = (-1.224744871391589*fr[27])+1.224744871391589*fl[27]+0.7071067811865475*fr[16]+0.7071067811865475*fl[16]; 
  favg[12] = (-1.224744871391589*fr[28])+1.224744871391589*fl[28]+0.7071067811865475*fr[17]+0.7071067811865475*fl[17]; 
  favg[13] = (-1.224744871391589*fr[29])+1.224744871391589*fl[29]+0.7071067811865475*fr[18]+0.7071067811865475*fl[18]; 
  favg[14] = (-1.224744871391589*fr[30])+1.224744871391589*fl[30]+0.7071067811865475*fr[19]+0.7071067811865475*fl[19]; 
  favg[15] = (-1.224744871391589*fr[31])+1.224744871391589*fl[31]+0.7071067811865475*fr[26]+0.7071067811865475*fl[26]; 

  alpha[0] = 2.0*(B1[0]*wv1+E2[0])-2.0*B0[0]*wv2; 
  alpha[1] = 2.0*(B1[1]*wv1+E2[1])-2.0*B0[1]*wv2; 
  alpha[2] = 2.0*(B1[2]*wv1+E2[2])-2.0*B0[2]*wv2; 
  alpha[3] = 0.5773502691896258*B1[0]*dv1; 
  alpha[4] = -0.5773502691896258*B0[0]*dv2; 
  alpha[5] = 2.0*(B1[3]*wv1+E2[3])-2.0*B0[3]*wv2; 
  alpha[6] = 0.5773502691896258*B1[1]*dv1; 
  alpha[7] = 0.5773502691896258*B1[2]*dv1; 
  alpha[8] = -0.5773502691896258*B0[1]*dv2; 
  alpha[9] = -0.5773502691896258*B0[2]*dv2; 
  alpha[11] = 0.5773502691896258*B1[3]*dv1; 
  alpha[12] = -0.5773502691896258*B0[3]*dv2; 

  gkyl_real amid = 0.25*alpha[0]; 

  Ghat[0] = 0.3535533905932737*(1.732050807568877*(fr[5]+fl[5])-1.0*fr[0]+fl[0])*amax+0.125*(alpha[12]*favg[12]+alpha[11]*favg[11]+alpha[9]*favg[9]+alpha[8]*favg[8]+alpha[7]*favg[7]+alpha[6]*favg[6]+alpha[5]*favg[5]+alpha[4]*favg[4]+alpha[3]*favg[3]+alpha[2]*favg[2]+alpha[1]*favg[1]+alpha[0]*favg[0]); 
  Ghat[1] = 0.3535533905932737*(1.732050807568877*(fr[12]+fl[12])-1.0*fr[1]+fl[1])*amax+0.125*(alpha[9]*favg[12]+favg[9]*alpha[12]+alpha[7]*favg[11]+favg[7]*alpha[11]+alpha[4]*favg[8]+favg[4]*alpha[8]+alpha[3]*favg[6]+favg[3]*alpha[6]+alpha[2]*favg[5]+favg[2]*alpha[5]+alpha[0]*favg[1]+favg[0]*alpha[1]); 
  Ghat[2] = 0.3535533905932737*(1.732050807568877*(fr[13]+fl[13])-1.0*fr[2]+fl[2])*amax+0.125*(alpha[8]*favg[12]+favg[8]*alpha[12]+alpha[6]*favg[11]+favg[6]*alpha[11]+alpha[4]*favg[9]+favg[4]*alpha[9]+alpha[3]*favg[7]+favg[3]*alpha[7]+alpha[1]*favg[5]+favg[1]*alpha[5]+alpha[0]*favg[2]+favg[0]*alpha[2]); 
  Ghat[3] = 0.3535533905932737*(1.732050807568877*(fr[14]+fl[14])-1.0*fr[3]+fl[3])*amax+0.125*(alpha[12]*favg[15]+alpha[9]*favg[14]+alpha[8]*favg[13]+alpha[5]*favg[11]+favg[5]*alpha[11]+alpha[4]*favg[10]+alpha[2]*favg[7]+favg[2]*alpha[7]+alpha[1]*favg[6]+favg[1]*alpha[6]+alpha[0]*favg[3]+favg[0]*alpha[3]); 
  Ghat[4] = 0.3535533905932737*(1.732050807568877*(fr[15]+fl[15])-1.0*fr[4]+fl[4])*amax+0.125*(alpha[11]*favg[15]+alpha[7]*favg[14]+alpha[6]*favg[13]+alpha[5]*favg[12]+favg[5]*alpha[12]+alpha[3]*favg[10]+alpha[2]*favg[9]+favg[2]*alpha[9]+alpha[1]*favg[8]+favg[1]*alpha[8]+alpha[0]*favg[4]+favg[0]*alpha[4]); 
  Ghat[5] = 0.3535533905932737*(1.732050807568877*(fr[20]+fl[20])-1.0*fr[6]+fl[6])*amax+0.125*(alpha[4]*favg[12]+favg[4]*alpha[12]+alpha[3]*favg[11]+favg[3]*alpha[11]+alpha[8]*favg[9]+favg[8]*alpha[9]+alpha[6]*favg[7]+favg[6]*alpha[7]+alpha[0]*favg[5]+favg[0]*alpha[5]+alpha[1]*favg[2]+favg[1]*alpha[2]); 
  Ghat[6] = 0.3535533905932737*(1.732050807568877*(fr[21]+fl[21])-1.0*fr[7]+fl[7])*amax+0.125*(alpha[9]*favg[15]+alpha[12]*favg[14]+alpha[4]*favg[13]+alpha[2]*favg[11]+favg[2]*alpha[11]+alpha[8]*favg[10]+alpha[5]*favg[7]+favg[5]*alpha[7]+alpha[0]*favg[6]+favg[0]*alpha[6]+alpha[1]*favg[3]+favg[1]*alpha[3]); 
  Ghat[7] = 0.3535533905932737*(1.732050807568877*(fr[22]+fl[22])-1.0*fr[8]+fl[8])*amax+0.125*(alpha[8]*favg[15]+alpha[4]*favg[14]+alpha[12]*favg[13]+alpha[1]*favg[11]+favg[1]*alpha[11]+alpha[9]*favg[10]+alpha[0]*favg[7]+favg[0]*alpha[7]+alpha[5]*favg[6]+favg[5]*alpha[6]+alpha[2]*favg[3]+favg[2]*alpha[3]); 
  Ghat[8] = 0.3535533905932737*(1.732050807568877*(fr[23]+fl[23])-1.0*fr[9]+fl[9])*amax+0.125*(alpha[7]*favg[15]+alpha[11]*favg[14]+alpha[3]*favg[13]+alpha[2]*favg[12]+favg[2]*alpha[12]+alpha[6]*favg[10]+alpha[5]*favg[9]+favg[5]*alpha[9]+alpha[0]*favg[8]+favg[0]*alpha[8]+alpha[1]*favg[4]+favg[1]*alpha[4]); 
  Ghat[9] = 0.3535533905932737*(1.732050807568877*(fr[24]+fl[24])-1.0*fr[10]+fl[10])*amax+0.125*(alpha[6]*favg[15]+alpha[3]*favg[14]+alpha[11]*favg[13]+alpha[1]*favg[12]+favg[1]*alpha[12]+alpha[7]*favg[10]+alpha[0]*favg[9]+favg[0]*alpha[9]+alpha[5]*favg[8]+favg[5]*alpha[8]+alpha[2]*favg[4]+favg[2]*alpha[4]); 
  Ghat[10] = 0.3535533905932737*(1.732050807568877*(fr[25]+fl[25])-1.0*fr[11]+fl[11])*amax+0.125*(alpha[5]*favg[15]+alpha[2]*favg[14]+alpha[1]*favg[13]+alpha[11]*favg[12]+favg[11]*alpha[12]+alpha[0]*favg[10]+alpha[7]*favg[9]+favg[7]*alpha[9]+alpha[6]*favg[8]+favg[6]*alpha[8]+alpha[3]*favg[4]+favg[3]*alpha[4]); 
  Ghat[11] = 0.3535533905932737*(1.732050807568877*(fr[27]+fl[27])-1.0*fr[16]+fl[16])*amax+0.125*(alpha[4]*favg[15]+alpha[8]*favg[14]+alpha[9]*favg[13]+favg[10]*alpha[12]+alpha[0]*favg[11]+favg[0]*alpha[11]+alpha[1]*favg[7]+favg[1]*alpha[7]+alpha[2]*favg[6]+favg[2]*alpha[6]+alpha[3]*favg[5]+favg[3]*alpha[5]); 
  Ghat[12] = 0.3535533905932737*(1.732050807568877*(fr[28]+fl[28])-1.0*fr[17]+fl[17])*amax+0.125*(alpha[3]*favg[15]+alpha[6]*favg[14]+alpha[7]*favg[13]+alpha[0]*favg[12]+favg[0]*alpha[12]+favg[10]*alpha[11]+alpha[1]*favg[9]+favg[1]*alpha[9]+alpha[2]*favg[8]+favg[2]*alpha[8]+alpha[4]*favg[5]+favg[4]*alpha[5]); 
  Ghat[13] = 0.3535533905932737*(1.732050807568877*(fr[29]+fl[29])-1.0*fr[18]+fl[18])*amax+0.125*(alpha[2]*favg[15]+alpha[5]*favg[14]+alpha[0]*favg[13]+alpha[7]*favg[12]+favg[7]*alpha[12]+alpha[9]*favg[11]+favg[9]*alpha[11]+alpha[1]*favg[10]+alpha[3]*favg[8]+favg[3]*alpha[8]+alpha[4]*favg[6]+favg[4]*alpha[6]); 
  Ghat[14] = 0.3535533905932737*(1.732050807568877*(fr[30]+fl[30])-1.0*fr[19]+fl[19])*amax+0.125*(alpha[1]*favg[15]+alpha[0]*favg[14]+alpha[5]*favg[13]+alpha[6]*favg[12]+favg[6]*alpha[12]+alpha[8]*favg[11]+favg[8]*alpha[11]+alpha[2]*favg[10]+alpha[3]*favg[9]+favg[3]*alpha[9]+alpha[4]*favg[7]+favg[4]*alpha[7]); 
  Ghat[15] = 0.3535533905932737*(1.732050807568877*(fr[31]+fl[31])-1.0*fr[26]+fl[26])*amax+0.125*(alpha[0]*favg[15]+alpha[1]*favg[14]+alpha[2]*favg[13]+alpha[3]*favg[12]+favg[3]*alpha[12]+alpha[4]*favg[11]+favg[4]*alpha[11]+alpha[5]*favg[10]+alpha[6]*favg[9]+favg[6]*alpha[9]+alpha[7]*favg[8]+favg[7]*alpha[8]); 

  outr[0] += 0.7071067811865475*Ghat[0]*dv12r; 
  outr[1] += 0.7071067811865475*Ghat[1]*dv12r; 
  outr[2] += 0.7071067811865475*Ghat[2]*dv12r; 
  outr[3] += 0.7071067811865475*Ghat[3]*dv12r; 
  outr[4] += 0.7071067811865475*Ghat[4]*dv12r; 
  outr[5] += -1.224744871391589*Ghat[0]*dv12r; 
  outr[6] += 0.7071067811865475*Ghat[5]*dv12r; 
  outr[7] += 0.7071067811865475*Ghat[6]*dv12r; 
  outr[8] += 0.7071067811865475*Ghat[7]*dv12r; 
  outr[9] += 0.7071067811865475*Ghat[8]*dv12r; 
  outr[10] += 0.7071067811865475*Ghat[9]*dv12r; 
  outr[11] += 0.7071067811865475*Ghat[10]*dv12r; 
  outr[12] += -1.224744871391589*Ghat[1]*dv12r; 
  outr[13] += -1.224744871391589*Ghat[2]*dv12r; 
  outr[14] += -1.224744871391589*Ghat[3]*dv12r; 
  outr[15] += -1.224744871391589*Ghat[4]*dv12r; 
  outr[16] += 0.7071067811865475*Ghat[11]*dv12r; 
  outr[17] += 0.7071067811865475*Ghat[12]*dv12r; 
  outr[18] += 0.7071067811865475*Ghat[13]*dv12r; 
  outr[19] += 0.7071067811865475*Ghat[14]*dv12r; 
  outr[20] += -1.224744871391589*Ghat[5]*dv12r; 
  outr[21] += -1.224744871391589*Ghat[6]*dv12r; 
  outr[22] += -1.224744871391589*Ghat[7]*dv12r; 
  outr[23] += -1.224744871391589*Ghat[8]*dv12r; 
  outr[24] += -1.224744871391589*Ghat[9]*dv12r; 
  outr[25] += -1.224744871391589*Ghat[10]*dv12r; 
  outr[26] += 0.7071067811865475*Ghat[15]*dv12r; 
  outr[27] += -1.224744871391589*Ghat[11]*dv12r; 
  outr[28] += -1.224744871391589*Ghat[12]*dv12r; 
  outr[29] += -1.224744871391589*Ghat[13]*dv12r; 
  outr[30] += -1.224744871391589*Ghat[14]*dv12r; 
  outr[31] += -1.224744871391589*Ghat[15]*dv12r; 

  outl[0] += -0.7071067811865475*Ghat[0]*dv12l; 
  outl[1] += -0.7071067811865475*Ghat[1]*dv12l; 
  outl[2] += -0.7071067811865475*Ghat[2]*dv12l; 
  outl[3] += -0.7071067811865475*Ghat[3]*dv12l; 
  outl[4] += -0.7071067811865475*Ghat[4]*dv12l; 
  outl[5] += -1.224744871391589*Ghat[0]*dv12l; 
  outl[6] += -0.7071067811865475*Ghat[5]*dv12l; 
  outl[7] += -0.7071067811865475*Ghat[6]*dv12l; 
  outl[8] += -0.7071067811865475*Ghat[7]*dv12l; 
  outl[9] += -0.7071067811865475*Ghat[8]*dv12l; 
  outl[10] += -0.7071067811865475*Ghat[9]*dv12l; 
  outl[11] += -0.7071067811865475*Ghat[10]*dv12l; 
  outl[12] += -1.224744871391589*Ghat[1]*dv12l; 
  outl[13] += -1.224744871391589*Ghat[2]*dv12l; 
  outl[14] += -1.224744871391589*Ghat[3]*dv12l; 
  outl[15] += -1.224744871391589*Ghat[4]*dv12l; 
  outl[16] += -0.7071067811865475*Ghat[11]*dv12l; 
  outl[17] += -0.7071067811865475*Ghat[12]*dv12l; 
  outl[18] += -0.7071067811865475*Ghat[13]*dv12l; 
  outl[19] += -0.7071067811865475*Ghat[14]*dv12l; 
  outl[20] += -1.224744871391589*Ghat[5]*dv12l; 
  outl[21] += -1.224744871391589*Ghat[6]*dv12l; 
  outl[22] += -1.224744871391589*Ghat[7]*dv12l; 
  outl[23] += -1.224744871391589*Ghat[8]*dv12l; 
  outl[24] += -1.224744871391589*Ghat[9]*dv12l; 
  outl[25] += -1.224744871391589*Ghat[10]*dv12l; 
  outl[26] += -0.7071067811865475*Ghat[15]*dv12l; 
  outl[27] += -1.224744871391589*Ghat[11]*dv12l; 
  outl[28] += -1.224744871391589*Ghat[12]*dv12l; 
  outl[29] += -1.224744871391589*Ghat[13]*dv12l; 
  outl[30] += -1.224744871391589*Ghat[14]*dv12l; 
  outl[31] += -1.224744871391589*Ghat[15]*dv12l; 

  return fabs(amid); 
} 
