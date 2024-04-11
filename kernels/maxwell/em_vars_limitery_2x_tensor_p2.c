#include <gkyl_maxwell_kernels.h> 
GKYL_CU_DH void em_vars_limitery_2x_tensor_p2(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, 
  double *ql, double *qc, double *qr) 
{ 
  // limiter_fac: Factor for relationship between cell slopes and cell average differences (by default: 1/sqrt(3))
  // wv_eqn:      Wave equation for computing waves for limiting characteristics
  // ql/c/r:      [Ex, Ey, Ez, Bx, By, Bz, phi, psi], Input state vector in left/center/right cells.

  double *exl = &ql[0]; 
  double *eyl = &ql[9]; 
  double *ezl = &ql[18]; 
  double *bxl = &ql[27]; 
  double *byl = &ql[36]; 
  double *bzl = &ql[45]; 
  double *phl = &ql[54]; 
  double *psl = &ql[63]; 
 
  double *exc = &qc[0]; 
  double *eyc = &qc[9]; 
  double *ezc = &qc[18]; 
  double *bxc = &qc[27]; 
  double *byc = &qc[36]; 
  double *bzc = &qc[45]; 
  double *phc = &qc[54]; 
  double *psc = &qc[63]; 
 
  double *exr = &qr[0]; 
  double *eyr = &qr[9]; 
  double *ezr = &qr[18]; 
  double *bxr = &qr[27]; 
  double *byr = &qr[36]; 
  double *bzr = &qr[45]; 
  double *phr = &qr[54]; 
  double *psr = &qr[63]; 

  double em_avg_l[8], em_avg_c[8], em_avg_r[8] = {0.0};
  double delta_l[8], delta_c[8], delta_r[8] = {0.0};
  double waves_slope_l[48], waves_slope_c[48], waves_slope_r[48] = {0.0};
  double speeds[6];

  em_avg_l[0] = exl[0];
  em_avg_l[1] = eyl[0];
  em_avg_l[2] = ezl[0];
  em_avg_l[3] = bxl[0];
  em_avg_l[4] = byl[0];
  em_avg_l[5] = bzl[0];
  em_avg_l[6] = phl[0];
  em_avg_l[7] = psl[0];

  em_avg_c[0] = exc[0];
  em_avg_c[1] = eyc[0];
  em_avg_c[2] = ezc[0];
  em_avg_c[3] = bxc[0];
  em_avg_c[4] = byc[0];
  em_avg_c[5] = bzc[0];
  em_avg_c[6] = phc[0];
  em_avg_c[7] = psc[0];

  em_avg_r[0] = exr[0];
  em_avg_r[1] = eyr[0];
  em_avg_r[2] = ezr[0];
  em_avg_r[3] = bxr[0];
  em_avg_r[4] = byr[0];
  em_avg_r[5] = bzr[0];
  em_avg_r[6] = phr[0];
  em_avg_r[7] = psr[0];

  delta_l[0] = limiter_fac*(em_avg_c[0] - em_avg_l[0]);
  delta_l[1] = limiter_fac*(em_avg_c[1] - em_avg_l[1]);
  delta_l[2] = limiter_fac*(em_avg_c[2] - em_avg_l[2]);
  delta_l[3] = limiter_fac*(em_avg_c[3] - em_avg_l[3]);
  delta_l[4] = limiter_fac*(em_avg_c[4] - em_avg_l[4]);
  delta_l[5] = limiter_fac*(em_avg_c[5] - em_avg_l[5]);
  delta_l[6] = limiter_fac*(em_avg_c[6] - em_avg_l[6]);
  delta_l[7] = limiter_fac*(em_avg_c[7] - em_avg_l[7]);

  delta_c[0] = exc[2];
  delta_c[1] = eyc[2];
  delta_c[2] = ezc[2];
  delta_c[3] = bxc[2];
  delta_c[4] = byc[2];
  delta_c[5] = bzc[2];
  delta_c[6] = phc[2];
  delta_c[7] = psc[2];

  delta_r[0] = limiter_fac*(em_avg_r[0] - em_avg_c[0]);
  delta_r[1] = limiter_fac*(em_avg_r[1] - em_avg_c[1]);
  delta_r[2] = limiter_fac*(em_avg_r[2] - em_avg_c[2]);
  delta_r[3] = limiter_fac*(em_avg_r[3] - em_avg_c[3]);
  delta_r[4] = limiter_fac*(em_avg_r[4] - em_avg_c[4]);
  delta_r[5] = limiter_fac*(em_avg_r[5] - em_avg_c[5]);
  delta_r[6] = limiter_fac*(em_avg_r[6] - em_avg_c[6]);
  delta_r[7] = limiter_fac*(em_avg_r[7] - em_avg_c[7]);

  double my_max_speed_l = gkyl_wv_eqn_waves(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, delta_l,
    em_avg_c, em_avg_c, waves_slope_l, speeds);
  double my_max_speed_c = gkyl_wv_eqn_waves(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, delta_c,
    em_avg_c, em_avg_c, waves_slope_c, speeds);
  double my_max_speed_r = gkyl_wv_eqn_waves(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, delta_r,
    em_avg_c, em_avg_c, waves_slope_r, speeds);

  // limit wave components, accumulating it to slope
  double mm[48] = {0.0};
  double slope[8] = {0.0};
  mm[0] = minmod(waves_slope_c[0], waves_slope_l[0], waves_slope_r[0]);
  mm[1] = minmod(waves_slope_c[1], waves_slope_l[1], waves_slope_r[1]);
  mm[2] = minmod(waves_slope_c[2], waves_slope_l[2], waves_slope_r[2]);
  mm[3] = minmod(waves_slope_c[3], waves_slope_l[3], waves_slope_r[3]);
  mm[4] = minmod(waves_slope_c[4], waves_slope_l[4], waves_slope_r[4]);
  mm[5] = minmod(waves_slope_c[5], waves_slope_l[5], waves_slope_r[5]);
  mm[6] = minmod(waves_slope_c[6], waves_slope_l[6], waves_slope_r[6]);
  mm[7] = minmod(waves_slope_c[7], waves_slope_l[7], waves_slope_r[7]);

  slope[0] += mm[0];
  slope[1] += mm[1];
  slope[2] += mm[2];
  slope[3] += mm[3];
  slope[4] += mm[4];
  slope[5] += mm[5];
  slope[6] += mm[6];
  slope[7] += mm[7];

  mm[8] = minmod(waves_slope_c[8], waves_slope_l[8], waves_slope_r[8]);
  mm[9] = minmod(waves_slope_c[9], waves_slope_l[9], waves_slope_r[9]);
  mm[10] = minmod(waves_slope_c[10], waves_slope_l[10], waves_slope_r[10]);
  mm[11] = minmod(waves_slope_c[11], waves_slope_l[11], waves_slope_r[11]);
  mm[12] = minmod(waves_slope_c[12], waves_slope_l[12], waves_slope_r[12]);
  mm[13] = minmod(waves_slope_c[13], waves_slope_l[13], waves_slope_r[13]);
  mm[14] = minmod(waves_slope_c[14], waves_slope_l[14], waves_slope_r[14]);
  mm[15] = minmod(waves_slope_c[15], waves_slope_l[15], waves_slope_r[15]);

  slope[0] += mm[8];
  slope[1] += mm[9];
  slope[2] += mm[10];
  slope[3] += mm[11];
  slope[4] += mm[12];
  slope[5] += mm[13];
  slope[6] += mm[14];
  slope[7] += mm[15];

  mm[16] = minmod(waves_slope_c[16], waves_slope_l[16], waves_slope_r[16]);
  mm[17] = minmod(waves_slope_c[17], waves_slope_l[17], waves_slope_r[17]);
  mm[18] = minmod(waves_slope_c[18], waves_slope_l[18], waves_slope_r[18]);
  mm[19] = minmod(waves_slope_c[19], waves_slope_l[19], waves_slope_r[19]);
  mm[20] = minmod(waves_slope_c[20], waves_slope_l[20], waves_slope_r[20]);
  mm[21] = minmod(waves_slope_c[21], waves_slope_l[21], waves_slope_r[21]);
  mm[22] = minmod(waves_slope_c[22], waves_slope_l[22], waves_slope_r[22]);
  mm[23] = minmod(waves_slope_c[23], waves_slope_l[23], waves_slope_r[23]);

  slope[0] += mm[16];
  slope[1] += mm[17];
  slope[2] += mm[18];
  slope[3] += mm[19];
  slope[4] += mm[20];
  slope[5] += mm[21];
  slope[6] += mm[22];
  slope[7] += mm[23];

  mm[24] = minmod(waves_slope_c[24], waves_slope_l[24], waves_slope_r[24]);
  mm[25] = minmod(waves_slope_c[25], waves_slope_l[25], waves_slope_r[25]);
  mm[26] = minmod(waves_slope_c[26], waves_slope_l[26], waves_slope_r[26]);
  mm[27] = minmod(waves_slope_c[27], waves_slope_l[27], waves_slope_r[27]);
  mm[28] = minmod(waves_slope_c[28], waves_slope_l[28], waves_slope_r[28]);
  mm[29] = minmod(waves_slope_c[29], waves_slope_l[29], waves_slope_r[29]);
  mm[30] = minmod(waves_slope_c[30], waves_slope_l[30], waves_slope_r[30]);
  mm[31] = minmod(waves_slope_c[31], waves_slope_l[31], waves_slope_r[31]);

  slope[0] += mm[24];
  slope[1] += mm[25];
  slope[2] += mm[26];
  slope[3] += mm[27];
  slope[4] += mm[28];
  slope[5] += mm[29];
  slope[6] += mm[30];
  slope[7] += mm[31];

  mm[32] = minmod(waves_slope_c[32], waves_slope_l[32], waves_slope_r[32]);
  mm[33] = minmod(waves_slope_c[33], waves_slope_l[33], waves_slope_r[33]);
  mm[34] = minmod(waves_slope_c[34], waves_slope_l[34], waves_slope_r[34]);
  mm[35] = minmod(waves_slope_c[35], waves_slope_l[35], waves_slope_r[35]);
  mm[36] = minmod(waves_slope_c[36], waves_slope_l[36], waves_slope_r[36]);
  mm[37] = minmod(waves_slope_c[37], waves_slope_l[37], waves_slope_r[37]);
  mm[38] = minmod(waves_slope_c[38], waves_slope_l[38], waves_slope_r[38]);
  mm[39] = minmod(waves_slope_c[39], waves_slope_l[39], waves_slope_r[39]);

  slope[0] += mm[32];
  slope[1] += mm[33];
  slope[2] += mm[34];
  slope[3] += mm[35];
  slope[4] += mm[36];
  slope[5] += mm[37];
  slope[6] += mm[38];
  slope[7] += mm[39];

  mm[40] = minmod(waves_slope_c[40], waves_slope_l[40], waves_slope_r[40]);
  mm[41] = minmod(waves_slope_c[41], waves_slope_l[41], waves_slope_r[41]);
  mm[42] = minmod(waves_slope_c[42], waves_slope_l[42], waves_slope_r[42]);
  mm[43] = minmod(waves_slope_c[43], waves_slope_l[43], waves_slope_r[43]);
  mm[44] = minmod(waves_slope_c[44], waves_slope_l[44], waves_slope_r[44]);
  mm[45] = minmod(waves_slope_c[45], waves_slope_l[45], waves_slope_r[45]);
  mm[46] = minmod(waves_slope_c[46], waves_slope_l[46], waves_slope_r[46]);
  mm[47] = minmod(waves_slope_c[47], waves_slope_l[47], waves_slope_r[47]);

  slope[0] += mm[40];
  slope[1] += mm[41];
  slope[2] += mm[42];
  slope[3] += mm[43];
  slope[4] += mm[44];
  slope[5] += mm[45];
  slope[6] += mm[46];
  slope[7] += mm[47];

  exc[2] = slope[0];
  eyc[2] = slope[1];
  ezc[2] = slope[2];
  bxc[2] = slope[3];
  byc[2] = slope[4];
  bzc[2] = slope[5];
  phc[2] = slope[6];
  psc[2] = slope[7];

  // Zero out x*y cross coefficient and all higher order coefficients if limiter applied
  for (int i=0; i<6; ++i) {
    if (mm[i*8] != waves_slope_c[i*8]) {
      exc[3] = 0.0;
      exc[4] = 0.0;
      exc[5] = 0.0;
      exc[6] = 0.0;
      exc[7] = 0.0;
      exc[8] = 0.0;
    }
    if (mm[i*8+1] != waves_slope_c[i*8+1]) {
      eyc[3] = 0.0;
      eyc[4] = 0.0;
      eyc[5] = 0.0;
      eyc[6] = 0.0;
      eyc[7] = 0.0;
      eyc[8] = 0.0;
    }
    if (mm[i*8+2] != waves_slope_c[i*8+2]) {
      ezc[3] = 0.0;
      ezc[4] = 0.0;
      ezc[5] = 0.0;
      ezc[6] = 0.0;
      ezc[7] = 0.0;
      ezc[8] = 0.0;
    }
    if (mm[i*8+3] != waves_slope_c[i*8+3]) {
      bxc[3] = 0.0;
      bxc[4] = 0.0;
      bxc[5] = 0.0;
      bxc[6] = 0.0;
      bxc[7] = 0.0;
      bxc[8] = 0.0;
    }
    if (mm[i*8+4] != waves_slope_c[i*8+4]) {
      byc[3] = 0.0;
      byc[4] = 0.0;
      byc[5] = 0.0;
      byc[6] = 0.0;
      byc[7] = 0.0;
      byc[8] = 0.0;
    }
    if (mm[i*8+5] != waves_slope_c[i*8+5]) {
      bzc[3] = 0.0;
      bzc[4] = 0.0;
      bzc[5] = 0.0;
      bzc[6] = 0.0;
      bzc[7] = 0.0;
      bzc[8] = 0.0;
    }
    if (mm[i*8+6] != waves_slope_c[i*8+6]) {
      phc[3] = 0.0;
      phc[4] = 0.0;
      phc[5] = 0.0;
      phc[6] = 0.0;
      phc[7] = 0.0;
      phc[8] = 0.0;
    }
    if (mm[i*8+7] != waves_slope_c[i*8+7]) {
      psc[3] = 0.0;
      psc[4] = 0.0;
      psc[5] = 0.0;
      psc[6] = 0.0;
      psc[7] = 0.0;
      psc[8] = 0.0;
    }
  }    
}