#include <gkyl_fpo_vlasov_kernels.h> 
GKYL_CU_DH void fpo_diff_coeff_recov_vy_2x3v_ser_p1(const double *dxv, const double *G_l, const double *G_c, const double *G_r, double *diff_coeff) { 
  // dxv[NDIM]: Cell spacing. 
  // G_l,c,r:   Input potential in left/center/right cells in recovery direction.

  const double dv1 = 2.0/dxv[3]; 
  const double dv1_sq = dv1*dv1; 
  double *diff_coeff_xx = &diff_coeff[0]; 
  double *diff_coeff_yy = &diff_coeff[32]; 
  double *diff_coeff_zz = &diff_coeff[64]; 
  diff_coeff_yy[0] = (-0.5412658773652741*G_r[4]*dv1_sq)+0.5412658773652741*G_l[4]*dv1_sq+0.5625*G_r[0]*dv1_sq+0.5625*G_l[0]*dv1_sq-1.125*G_c[0]*dv1_sq; 
  diff_coeff_yy[1] = (-0.5412658773652741*G_r[9]*dv1_sq)+0.5412658773652741*G_l[9]*dv1_sq+0.5625*G_r[1]*dv1_sq+0.5625*G_l[1]*dv1_sq-1.125*G_c[1]*dv1_sq; 
  diff_coeff_yy[2] = (-0.5412658773652741*G_r[10]*dv1_sq)+0.5412658773652741*G_l[10]*dv1_sq+0.5625*G_r[2]*dv1_sq+0.5625*G_l[2]*dv1_sq-1.125*G_c[2]*dv1_sq; 
  diff_coeff_yy[3] = (-0.5412658773652741*G_r[11]*dv1_sq)+0.5412658773652741*G_l[11]*dv1_sq+0.5625*G_r[3]*dv1_sq+0.5625*G_l[3]*dv1_sq-1.125*G_c[3]*dv1_sq; 
  diff_coeff_yy[4] = (-0.4375*G_r[4]*dv1_sq)-0.4375*G_l[4]*dv1_sq-2.875*G_c[4]*dv1_sq+0.5412658773652741*G_r[0]*dv1_sq-0.5412658773652741*G_l[0]*dv1_sq; 
  diff_coeff_yy[5] = (-0.5412658773652741*G_r[15]*dv1_sq)+0.5412658773652741*G_l[15]*dv1_sq+0.5625*G_r[5]*dv1_sq+0.5625*G_l[5]*dv1_sq-1.125*G_c[5]*dv1_sq; 
  diff_coeff_yy[6] = (-0.5412658773652741*G_r[17]*dv1_sq)+0.5412658773652741*G_l[17]*dv1_sq+0.5625*G_r[6]*dv1_sq+0.5625*G_l[6]*dv1_sq-1.125*G_c[6]*dv1_sq; 
  diff_coeff_yy[7] = (-0.5412658773652741*G_r[18]*dv1_sq)+0.5412658773652741*G_l[18]*dv1_sq+0.5625*G_r[7]*dv1_sq+0.5625*G_l[7]*dv1_sq-1.125*G_c[7]*dv1_sq; 
  diff_coeff_yy[8] = (-0.5412658773652741*G_r[19]*dv1_sq)+0.5412658773652741*G_l[19]*dv1_sq+0.5625*G_r[8]*dv1_sq+0.5625*G_l[8]*dv1_sq-1.125*G_c[8]*dv1_sq; 
  diff_coeff_yy[9] = (-0.4375*G_r[9]*dv1_sq)-0.4375*G_l[9]*dv1_sq-2.875*G_c[9]*dv1_sq+0.5412658773652741*G_r[1]*dv1_sq-0.5412658773652741*G_l[1]*dv1_sq; 
  diff_coeff_yy[10] = (-0.4375*G_r[10]*dv1_sq)-0.4375*G_l[10]*dv1_sq-2.875*G_c[10]*dv1_sq+0.5412658773652741*G_r[2]*dv1_sq-0.5412658773652741*G_l[2]*dv1_sq; 
  diff_coeff_yy[11] = (-0.4375*G_r[11]*dv1_sq)-0.4375*G_l[11]*dv1_sq-2.875*G_c[11]*dv1_sq+0.5412658773652741*G_r[3]*dv1_sq-0.5412658773652741*G_l[3]*dv1_sq; 
  diff_coeff_yy[12] = (-0.5412658773652741*G_r[23]*dv1_sq)+0.5412658773652741*G_l[23]*dv1_sq+0.5625*G_r[12]*dv1_sq+0.5625*G_l[12]*dv1_sq-1.125*G_c[12]*dv1_sq; 
  diff_coeff_yy[13] = (-0.5412658773652741*G_r[24]*dv1_sq)+0.5412658773652741*G_l[24]*dv1_sq+0.5625*G_r[13]*dv1_sq+0.5625*G_l[13]*dv1_sq-1.125*G_c[13]*dv1_sq; 
  diff_coeff_yy[14] = (-0.5412658773652741*G_r[25]*dv1_sq)+0.5412658773652741*G_l[25]*dv1_sq+0.5625*G_r[14]*dv1_sq+0.5625*G_l[14]*dv1_sq-1.125*G_c[14]*dv1_sq; 
  diff_coeff_yy[15] = (-0.4375*G_r[15]*dv1_sq)-0.4375*G_l[15]*dv1_sq-2.875*G_c[15]*dv1_sq+0.5412658773652741*G_r[5]*dv1_sq-0.5412658773652741*G_l[5]*dv1_sq; 
  diff_coeff_yy[16] = (-0.5412658773652741*G_r[26]*dv1_sq)+0.5412658773652741*G_l[26]*dv1_sq+0.5625*G_r[16]*dv1_sq+0.5625*G_l[16]*dv1_sq-1.125*G_c[16]*dv1_sq; 
  diff_coeff_yy[17] = (-0.4375*G_r[17]*dv1_sq)-0.4375*G_l[17]*dv1_sq-2.875*G_c[17]*dv1_sq+0.5412658773652741*G_r[6]*dv1_sq-0.5412658773652741*G_l[6]*dv1_sq; 
  diff_coeff_yy[18] = (-0.4375*G_r[18]*dv1_sq)-0.4375*G_l[18]*dv1_sq-2.875*G_c[18]*dv1_sq+0.5412658773652741*G_r[7]*dv1_sq-0.5412658773652741*G_l[7]*dv1_sq; 
  diff_coeff_yy[19] = (-0.4375*G_r[19]*dv1_sq)-0.4375*G_l[19]*dv1_sq-2.875*G_c[19]*dv1_sq+0.5412658773652741*G_r[8]*dv1_sq-0.5412658773652741*G_l[8]*dv1_sq; 
  diff_coeff_yy[20] = (-0.5412658773652741*G_r[28]*dv1_sq)+0.5412658773652741*G_l[28]*dv1_sq+0.5625*G_r[20]*dv1_sq+0.5625*G_l[20]*dv1_sq-1.125*G_c[20]*dv1_sq; 
  diff_coeff_yy[21] = (-0.5412658773652741*G_r[29]*dv1_sq)+0.5412658773652741*G_l[29]*dv1_sq+0.5625*G_r[21]*dv1_sq+0.5625*G_l[21]*dv1_sq-1.125*G_c[21]*dv1_sq; 
  diff_coeff_yy[22] = (-0.5412658773652741*G_r[30]*dv1_sq)+0.5412658773652741*G_l[30]*dv1_sq+0.5625*G_r[22]*dv1_sq+0.5625*G_l[22]*dv1_sq-1.125*G_c[22]*dv1_sq; 
  diff_coeff_yy[23] = (-0.4375*G_r[23]*dv1_sq)-0.4375*G_l[23]*dv1_sq-2.875*G_c[23]*dv1_sq+0.5412658773652741*G_r[12]*dv1_sq-0.5412658773652741*G_l[12]*dv1_sq; 
  diff_coeff_yy[24] = (-0.4375*G_r[24]*dv1_sq)-0.4375*G_l[24]*dv1_sq-2.875*G_c[24]*dv1_sq+0.5412658773652741*G_r[13]*dv1_sq-0.5412658773652741*G_l[13]*dv1_sq; 
  diff_coeff_yy[25] = (-0.4375*G_r[25]*dv1_sq)-0.4375*G_l[25]*dv1_sq-2.875*G_c[25]*dv1_sq+0.5412658773652741*G_r[14]*dv1_sq-0.5412658773652741*G_l[14]*dv1_sq; 
  diff_coeff_yy[26] = (-0.4375*G_r[26]*dv1_sq)-0.4375*G_l[26]*dv1_sq-2.875*G_c[26]*dv1_sq+0.5412658773652741*G_r[16]*dv1_sq-0.5412658773652741*G_l[16]*dv1_sq; 
  diff_coeff_yy[27] = (-0.5412658773652741*G_r[31]*dv1_sq)+0.5412658773652741*G_l[31]*dv1_sq+0.5625*G_r[27]*dv1_sq+0.5625*G_l[27]*dv1_sq-1.125*G_c[27]*dv1_sq; 
  diff_coeff_yy[28] = (-0.4375*G_r[28]*dv1_sq)-0.4375*G_l[28]*dv1_sq-2.875*G_c[28]*dv1_sq+0.5412658773652741*G_r[20]*dv1_sq-0.5412658773652741*G_l[20]*dv1_sq; 
  diff_coeff_yy[29] = (-0.4375*G_r[29]*dv1_sq)-0.4375*G_l[29]*dv1_sq-2.875*G_c[29]*dv1_sq+0.5412658773652741*G_r[21]*dv1_sq-0.5412658773652741*G_l[21]*dv1_sq; 
  diff_coeff_yy[30] = (-0.4375*G_r[30]*dv1_sq)-0.4375*G_l[30]*dv1_sq-2.875*G_c[30]*dv1_sq+0.5412658773652741*G_r[22]*dv1_sq-0.5412658773652741*G_l[22]*dv1_sq; 
  diff_coeff_yy[31] = (-0.4375*G_r[31]*dv1_sq)-0.4375*G_l[31]*dv1_sq-2.875*G_c[31]*dv1_sq+0.5412658773652741*G_r[27]*dv1_sq-0.5412658773652741*G_l[27]*dv1_sq; 
} 

