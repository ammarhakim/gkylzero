#include <gkyl_pkpm_fluid_kernels.h> 
#include <gkyl_pkpm_harm_coeffs.h> 
#include <gkyl_pkpm_flux_surf_tables_1x1v_p1.h> 
GKYL_CU_DH void pkpm_fluid_flux_accum_x_1x1v_ser_p1_g(int tid, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT G_l, double* GKYL_RESTRICT G_r) 
{ 
  if (tid >= 3) return; 
  const int a = tid/3; 
  const int j = tid - a*3; 
  double g_l = 0.0; 
  double g_r = 0.0; 
  for (int q = vst_pkpm_1x1v_p1_ph_x0_aoff[a]; q < vst_pkpm_1x1v_p1_ph_x0_aoff[a+1]; ++q) { 
    const int k = vst_pkpm_1x1v_p1_ph_x0_aks[q]; 
    g_l += vst_pkpm_1x1v_p1_ph_x0_Wl[k*3 + j]*f_l[k]; 
    g_r += vst_pkpm_1x1v_p1_ph_x0_Wr[k*3 + j]*f_r[k]; 
  } 
  G_l[j*1 + a] = g_l; 
  G_r[j*1 + a] = g_r; 
} 

GKYL_CU_DH double pkpm_fluid_flux_accum_x_1x1v_ser_p1(const double *w, const double *dxv,
  const double *jacob_pos_l, const double *jacob_pos_r, const double *b_l, const double *b_r,
  const double *u_l, const double *u_r, const double *xi_l, const double *xi_r,
  const double *theta_l, const double *theta_r, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT I) 
{ 
  // w[NDIM]:   Cell-center coordinates of the phase-space cell (r side). 
  // dxv[NDIM]: Cell spacing of the phase-space cell. 
  // jacob_pos_l/r: Position-map Jacobian per configuration direction (cdim blocks of 2 equal 
  //            entries) of the left/right cell. 
  // b_l/r:     Magnetic field unit vector [bx, by, bz] of the left/right cell. 
  // u_l/r:     Flow velocity [ux, uy, uz] of the left/right cell. 
  // xi_l/r:    Thermal speed [xi, xi^2, p] of the left/right cell (the first expansion is used). 
  // theta_l/r: T_perp/T of the left/right cell. 
  // f_l/r:     Jacobian-weighted normalized distributions J [F, G, M0_perp, N_perp] of the left/right 
  //            phase-space cell. 
  // I:         Accumulated nodal moments of the kinetic HLL flux at the 1 transverse face nodes: 
  //            I[k*1 + n], k = 0..6 for [int Phi_F, int w Phi_F, int w^2 Phi_F, int Phi_G, int Phi_M0 (x, y, z)] 
  //            of the left contribution Phi^l and k = 7..13 the same of the right one; the traces are 
  //            divided by J_x of their cell, so the moments carry the face Jacobian J/J_x. 
  // Returns the largest width S_R - S_L of the speed bounds over the nodes (for the CFL rate). 
  // Gauss-Legendre weights in w times 1/2 (the w integral of a cell is dv sum_j wg[j] (.)_j) 
  const double wg[3] = {0.2777777777777778, 0.4444444444444444, 0.2777777777777778}; 
  double alpha_max = 0.0; 
  const double jdinv_l = 1.0/jacob_pos_l[0], jdinv_r = 1.0/jacob_pos_r[0]; 
  const double wv = w[1], dv = dxv[1]; 
  // w of the nodes: wv + wn[j] dv 
  const double wn[3] = {-0.38729833462074165, 0.0, 0.38729833462074165}; 
  const double *b0_l = &b_l[0], *b0_r = &b_r[0]; 
  const double *b1_l = &b_l[2], *b1_r = &b_r[2]; 
  const double *b2_l = &b_l[4], *b2_r = &b_r[4]; 
  const double *ud_l = &u_l[0], *ud_r = &u_r[0]; 
  const double *xd_l = &xi_l[0], *xd_r = &xi_r[0]; 
 
  // Traces of b, u_x, xi and theta at the transverse nodes: the left cell at its upper 
  // edge, the current cell at its lower edge. 
  double bq0_l[1], bq0_r[1]; 
  bq0_l[0] = 1.224744871391589*b0_l[1]+0.7071067811865475*b0_l[0]; 
  bq0_r[0] = 0.7071067811865475*b0_r[0]-1.224744871391589*b0_r[1]; 
  double bq1_l[1], bq1_r[1]; 
  bq1_l[0] = 1.224744871391589*b1_l[1]+0.7071067811865475*b1_l[0]; 
  bq1_r[0] = 0.7071067811865475*b1_r[0]-1.224744871391589*b1_r[1]; 
  double bq2_l[1], bq2_r[1]; 
  bq2_l[0] = 1.224744871391589*b2_l[1]+0.7071067811865475*b2_l[0]; 
  bq2_r[0] = 0.7071067811865475*b2_r[0]-1.224744871391589*b2_r[1]; 
  double uq_l[1], uq_r[1], xq_l[1], xq_r[1], thq_l[1], thq_r[1]; 
  uq_l[0] = 1.224744871391589*ud_l[1]+0.7071067811865475*ud_l[0]; 
  uq_r[0] = 0.7071067811865475*ud_r[0]-1.224744871391589*ud_r[1]; 
  xq_l[0] = 1.224744871391589*xd_l[1]+0.7071067811865475*xd_l[0]; 
  xq_r[0] = 0.7071067811865475*xd_r[0]-1.224744871391589*xd_r[1]; 
  thq_l[0] = 1.224744871391589*theta_l[1]+0.7071067811865475*theta_l[0]; 
  thq_r[0] = 0.7071067811865475*theta_r[0]-1.224744871391589*theta_r[1]; 
 
  // Trace factors of the eight blocks, scaled by 1/J_x of their cell. 
  double GU_l[8][3], GU_r[8][3]; 
  for (int c = 0; c < 8; ++c) { 
    for (int tid = 0; tid < 3; ++tid) pkpm_fluid_flux_accum_x_1x1v_ser_p1_g(tid, &f_l[c*6], &f_r[c*6], GU_l[c], GU_r[c]); 
    for (int k = 0; k < 3; ++k) { GU_l[c][k] *= jdinv_l; GU_r[c][k] *= jdinv_r; } 
  } 
 
  for (int i = 0; i < 1; ++i) { 
    const double bn_l[3] = {bq0_l[i], bq1_l[i], bq2_l[i]}, bn_r[3] = {bq0_r[i], bq1_r[i], bq2_r[i]}; 
    for (int j = 0; j < 3; ++j) { 
      const double wq = wv + wn[j]*dv; 
      double X_l[8], X_r[8]; 
      for (int c = 0; c < 8; ++c) { 
        double x_l = 0.0, x_r = 0.0; 
        for (int a = 0; a < 1; ++a) { 
          x_l += vst_pkpm_1x1v_p1_ph_x0_Cm[i*1 + a]*GU_l[c][j*1 + a]; 
          x_r += vst_pkpm_1x1v_p1_ph_x0_Cm[i*1 + a]*GU_r[c][j*1 + a]; 
        } 
        X_l[c] = x_l; X_r[c] = x_r; 
      } 
      double phi_l[8], phi_r[8], sb[2]; 
      pkpm_harm_conf_flux(0, X_l, X_r, xq_l[i]*bn_l[0]*wq + uq_l[i], xq_r[i]*bn_r[0]*wq + uq_r[i], 
        xq_l[i], xq_r[i], bn_l, bn_r, thq_l[i], thq_r[i], phi_l, phi_r, sb); 
      const double wgt = wg[j]*dv; 
      I[0 + i] += wgt*phi_l[0]; 
      I[1 + i] += wgt*wq*phi_l[0]; 
      I[2 + i] += wgt*wq*wq*phi_l[0]; 
      I[3 + i] += wgt*phi_l[1]; 
      for (int a = 0; a < 3; ++a) I[(4 + a)*1 + i] += wgt*phi_l[2 + a]; 
      I[7 + i] += wgt*phi_r[0]; 
      I[8 + i] += wgt*wq*phi_r[0]; 
      I[9 + i] += wgt*wq*wq*phi_r[0]; 
      I[10 + i] += wgt*phi_r[1]; 
      for (int a = 0; a < 3; ++a) I[(11 + a)*1 + i] += wgt*phi_r[2 + a]; 
      alpha_max = fmax(alpha_max, sb[1] - sb[0]); 
    } 
  } 
  return alpha_max; 
} 
