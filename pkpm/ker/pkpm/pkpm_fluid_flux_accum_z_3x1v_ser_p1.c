#include <gkyl_pkpm_fluid_kernels.h> 
#include <gkyl_pkpm_harm_coeffs.h> 
#include <gkyl_pkpm_flux_surf_tables_3x1v_p1.h> 
GKYL_CU_DH void pkpm_fluid_flux_accum_z_3x1v_ser_p1_g(int tid, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT G_l, double* GKYL_RESTRICT G_r) 
{ 
  if (tid >= 12) return; 
  const int a = tid/3; 
  const int j = tid - a*3; 
  double g_l = 0.0; 
  double g_r = 0.0; 
  for (int q = vst_pkpm_3x1v_p1_ph_x2_aoff[a]; q < vst_pkpm_3x1v_p1_ph_x2_aoff[a+1]; ++q) { 
    const int k = vst_pkpm_3x1v_p1_ph_x2_aks[q]; 
    g_l += vst_pkpm_3x1v_p1_ph_x2_Wl[k*3 + j]*f_l[k]; 
    g_r += vst_pkpm_3x1v_p1_ph_x2_Wr[k*3 + j]*f_r[k]; 
  } 
  G_l[j*4 + a] = g_l; 
  G_r[j*4 + a] = g_r; 
} 

GKYL_CU_DH double pkpm_fluid_flux_accum_z_3x1v_ser_p1(const double *w, const double *dxv,
  const double *jacob_pos_l, const double *jacob_pos_r, const double *fc, const double *f_l,
  const double *f_r, double* GKYL_RESTRICT I) 
{ 
  // w[NDIM]:   Cell-center coordinates of the phase-space cell (r side). 
  // dxv[NDIM]: Cell spacing of the phase-space cell. 
  // jacob_pos_l/r: Position-map Jacobian per configuration direction (cdim blocks of 2 equal 
  //            entries) of the left/right cell. 
  // fc:        Face coefficients of the face (pkpm_conf_face_coef_z_3x_ser_p1). 
  // f_l/r:     Jacobian-weighted normalized distributions J [F, G, M0_perp, N_perp] of the left/right 
  //            phase-space cell. 
  // I:         Accumulated nodal moments of the kinetic HLL flux at the 25 transverse face nodes: 
  //            I[k*25 + n], k = 0..6 for [int Phi_F, int w Phi_F, int w^2 Phi_F, int Phi_G, int Phi_M0 (x, y, z)] 
  //            of the left contribution Phi^l and k = 7..13 the same of the right one; the traces are 
  //            divided by J_z of their cell, so the moments carry the face Jacobian J/J_z. 
  // Returns the largest width S_R - S_L of the speed bounds over the nodes (for the CFL rate). 
  // Gauss-Legendre weights in w times 1/2 (the w integral of a cell is dv sum_j wg[j] (.)_j) 
  const double wg[3] = {0.2777777777777778, 0.4444444444444444, 0.2777777777777778}; 
  double alpha_max = 0.0; 
  const double jdinv_l = 1.0/jacob_pos_l[4], jdinv_r = 1.0/jacob_pos_r[4]; 
  const double wv = w[3], dv = dxv[3]; 
  // w of the nodes: wv + wn[j] dv 
  const double wn[3] = {-0.38729833462074165, 0.0, 0.38729833462074165}; 
 
  // Trace factors of the eight blocks, scaled by 1/J_z of their cell. 
  double GU_l[8][12], GU_r[8][12]; 
  for (int c = 0; c < 8; ++c) { 
    for (int tid = 0; tid < 12; ++tid) pkpm_fluid_flux_accum_z_3x1v_ser_p1_g(tid, &f_l[c*24], &f_r[c*24], GU_l[c], GU_r[c]); 
    for (int k = 0; k < 12; ++k) { GU_l[c][k] *= jdinv_l; GU_r[c][k] *= jdinv_r; } 
  } 
 
  for (int i = 0; i < 25; ++i) { 
    for (int j = 0; j < 3; ++j) { 
      const double wq = wv + wn[j]*dv; 
      double X_l[8], X_r[8]; 
      for (int c = 0; c < 8; ++c) { 
        double x_l = 0.0, x_r = 0.0; 
        for (int a = 0; a < 4; ++a) { 
          x_l += vst_pkpm_3x1v_p1_ph_x2_Cm[i*4 + a]*GU_l[c][j*4 + a]; 
          x_r += vst_pkpm_3x1v_p1_ph_x2_Cm[i*4 + a]*GU_r[c][j*4 + a]; 
        } 
        X_l[c] = x_l; X_r[c] = x_r; 
      } 
      double phi_l[8], phi_r[8], sb[2]; 
      pkpm_harm_conf_flux(2, X_l, X_r, wq, &fc[i], 25, phi_l, phi_r, sb); 
      const double wgt = wg[j]*dv; 
      I[0 + i] += wgt*phi_l[0]; 
      I[25 + i] += wgt*wq*phi_l[0]; 
      I[50 + i] += wgt*wq*wq*phi_l[0]; 
      I[75 + i] += wgt*phi_l[1]; 
      for (int a = 0; a < 3; ++a) I[(4 + a)*25 + i] += wgt*phi_l[2 + a]; 
      I[175 + i] += wgt*phi_r[0]; 
      I[200 + i] += wgt*wq*phi_r[0]; 
      I[225 + i] += wgt*wq*wq*phi_r[0]; 
      I[250 + i] += wgt*phi_r[1]; 
      for (int a = 0; a < 3; ++a) I[(11 + a)*25 + i] += wgt*phi_r[2 + a]; 
      alpha_max = fmax(alpha_max, sb[1] - sb[0]); 
    } 
  } 
  return alpha_max; 
} 
