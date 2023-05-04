#include <gkyl_fpo_vlasov_kernels.h> 
GKYL_CU_DH void fpo_diff_coeff_recov_vz_1x3v_ser_p1(const double *dxv, const double *G_l, const double *G_c, const double *G_r, double *diff_coeff) { 
  // dxv[NDIM]: Cell spacing. 
  // G_l,c,r:   Input potential in left/center/right cells in recovery direction.

  const double dv1 = 2.0/dxv[3]; 
  const double dv1_sq = dv1*dv1; 
  double *diff_coeff_xx = &diff_coeff[0]; 
  double *diff_coeff_yy = &diff_coeff[16]; 
  double *diff_coeff_zz = &diff_coeff[32]; 
  diff_coeff_zz[0] = (-0.5412658773652741*G_r[4]*dv1_sq)+0.5412658773652741*G_l[4]*dv1_sq+0.5625*G_r[0]*dv1_sq+0.5625*G_l[0]*dv1_sq-1.125*G_c[0]*dv1_sq; 
  diff_coeff_zz[1] = (-0.5412658773652741*G_r[8]*dv1_sq)+0.5412658773652741*G_l[8]*dv1_sq+0.5625*G_r[1]*dv1_sq+0.5625*G_l[1]*dv1_sq-1.125*G_c[1]*dv1_sq; 
  diff_coeff_zz[2] = (-0.5412658773652741*G_r[9]*dv1_sq)+0.5412658773652741*G_l[9]*dv1_sq+0.5625*G_r[2]*dv1_sq+0.5625*G_l[2]*dv1_sq-1.125*G_c[2]*dv1_sq; 
  diff_coeff_zz[3] = (-0.5412658773652741*G_r[10]*dv1_sq)+0.5412658773652741*G_l[10]*dv1_sq+0.5625*G_r[3]*dv1_sq+0.5625*G_l[3]*dv1_sq-1.125*G_c[3]*dv1_sq; 
  diff_coeff_zz[4] = (-0.4375*G_r[4]*dv1_sq)-0.4375*G_l[4]*dv1_sq-2.875*G_c[4]*dv1_sq+0.5412658773652741*G_r[0]*dv1_sq-0.5412658773652741*G_l[0]*dv1_sq; 
  diff_coeff_zz[5] = (-0.5412658773652741*G_r[12]*dv1_sq)+0.5412658773652741*G_l[12]*dv1_sq+0.5625*G_r[5]*dv1_sq+0.5625*G_l[5]*dv1_sq-1.125*G_c[5]*dv1_sq; 
  diff_coeff_zz[6] = (-0.5412658773652741*G_r[13]*dv1_sq)+0.5412658773652741*G_l[13]*dv1_sq+0.5625*G_r[6]*dv1_sq+0.5625*G_l[6]*dv1_sq-1.125*G_c[6]*dv1_sq; 
  diff_coeff_zz[7] = (-0.5412658773652741*G_r[14]*dv1_sq)+0.5412658773652741*G_l[14]*dv1_sq+0.5625*G_r[7]*dv1_sq+0.5625*G_l[7]*dv1_sq-1.125*G_c[7]*dv1_sq; 
  diff_coeff_zz[8] = (-0.4375*G_r[8]*dv1_sq)-0.4375*G_l[8]*dv1_sq-2.875*G_c[8]*dv1_sq+0.5412658773652741*G_r[1]*dv1_sq-0.5412658773652741*G_l[1]*dv1_sq; 
  diff_coeff_zz[9] = (-0.4375*G_r[9]*dv1_sq)-0.4375*G_l[9]*dv1_sq-2.875*G_c[9]*dv1_sq+0.5412658773652741*G_r[2]*dv1_sq-0.5412658773652741*G_l[2]*dv1_sq; 
  diff_coeff_zz[10] = (-0.4375*G_r[10]*dv1_sq)-0.4375*G_l[10]*dv1_sq-2.875*G_c[10]*dv1_sq+0.5412658773652741*G_r[3]*dv1_sq-0.5412658773652741*G_l[3]*dv1_sq; 
  diff_coeff_zz[11] = (-0.5412658773652741*G_r[15]*dv1_sq)+0.5412658773652741*G_l[15]*dv1_sq+0.5625*G_r[11]*dv1_sq+0.5625*G_l[11]*dv1_sq-1.125*G_c[11]*dv1_sq; 
  diff_coeff_zz[12] = (-0.4375*G_r[12]*dv1_sq)-0.4375*G_l[12]*dv1_sq-2.875*G_c[12]*dv1_sq+0.5412658773652741*G_r[5]*dv1_sq-0.5412658773652741*G_l[5]*dv1_sq; 
  diff_coeff_zz[13] = (-0.4375*G_r[13]*dv1_sq)-0.4375*G_l[13]*dv1_sq-2.875*G_c[13]*dv1_sq+0.5412658773652741*G_r[6]*dv1_sq-0.5412658773652741*G_l[6]*dv1_sq; 
  diff_coeff_zz[14] = (-0.4375*G_r[14]*dv1_sq)-0.4375*G_l[14]*dv1_sq-2.875*G_c[14]*dv1_sq+0.5412658773652741*G_r[7]*dv1_sq-0.5412658773652741*G_l[7]*dv1_sq; 
  diff_coeff_zz[15] = (-0.4375*G_r[15]*dv1_sq)-0.4375*G_l[15]*dv1_sq-2.875*G_c[15]*dv1_sq+0.5412658773652741*G_r[11]*dv1_sq-0.5412658773652741*G_l[11]*dv1_sq; 
} 

