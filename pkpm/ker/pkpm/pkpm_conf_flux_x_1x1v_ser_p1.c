#include <gkyl_pkpm_kinetic_kernels.h> 
#include <gkyl_pkpm_harm_coeffs.h> 
#include <gkyl_pkpm_flux_surf_tables_1x1v_p1.h> 
GKYL_CU_DH void pkpm_conf_flux_x_1x1v_ser_p1_g(int tid, const double *f_l, const double *f_r,
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

GKYL_CU_DH void pkpm_conf_flux_x_1x1v_ser_p1_prj(int tid, const double *Fhat_nodal, double* GKYL_RESTRICT flux) 
{ 
  if (tid >= 3) return; 
  const int b = tid; 
  double t[1]; 
  for (int i = 0; i < 1; ++i) t[i] = 0.0; 
  for (int j = 0; j < 3; ++j) { 
    const double w = vst_pkpm_1x1v_p1_prj_x0_Vw[j*3 + b]; 
    for (int i = 0; i < 1; ++i) t[i] += w*Fhat_nodal[i*3 + j]; 
  } 
  for (int q = vst_pkpm_1x1v_p1_prj_x0_boff[b]; q < vst_pkpm_1x1v_p1_prj_x0_boff[b+1]; ++q) { 
    const int k = vst_pkpm_1x1v_p1_prj_x0_bks[q]; 
    const int a = vst_pkpm_1x1v_p1_prj_x0_kamap[k]; 
    double g = 0.0; 
    for (int i = 0; i < 1; ++i) g += vst_pkpm_1x1v_p1_prj_x0_Cw[i*1 + a]*t[i]; 
    flux[0 + k] = g; 
  } 
} 

GKYL_CU_DH double pkpm_conf_flux_x_1x1v_ser_p1(const double *w, const double *dxv,
  const double *jacob_pos_l, const double *jacob_pos_r, const double *fc, const double *f_l,
  const double *f_r, double* GKYL_RESTRICT flux) 
{ 
  // w[NDIM]:   Cell-center coordinates (current cell, computational). 
  // dxv[NDIM]: Cell spacing (current cell, computational). 
  // jacob_pos_l/r: Position-map Jacobian per configuration direction (cdim blocks of 2 equal 
  //            entries) of the left/current cell; the total Jacobian J is their product. 
  // fc:        Face coefficients of the face (pkpm_conf_face_coef_x_1x_ser_p1). 
  // f_l/f_r:   Jacobian-weighted normalized distributions J [F, G, M0_perp, N_perp] in the left/current cell. 
  // flux:      Output modal surface fluxes of the eight evolved blocks on the lower face of the 
  //            current cell in direction x: the physical HLL flux of the physical traces 
  //            (J f)/J times the face Jacobian J/J_x, i.e. the traces scaled by 1/J_x of their cell. 
  // Returns the stable time-step estimate from the largest speed bound over the nodes. 
  const double dx1 = 2.0/dxv[0]; 
  double Fhat[8][3]; 
  double alpha_max = 0.0; 
  const double jdinv_l = 1.0/jacob_pos_l[0], jdinv_r = 1.0/jacob_pos_r[0]; 
  const double wv = w[1], dv = dxv[1]; 
  // w of the nodes: wv + wn[j] dv 
  const double wn[3] = {-0.38729833462074165, 0.0, 0.38729833462074165}; 
 
  // Trace factors of the eight blocks, scaled by 1/J_x of their cell. 
  double GU_l[8][3], GU_r[8][3]; 
  for (int c = 0; c < 8; ++c) { 
    for (int tid = 0; tid < 3; ++tid) pkpm_conf_flux_x_1x1v_ser_p1_g(tid, &f_l[c*6], &f_r[c*6], GU_l[c], GU_r[c]); 
    for (int k = 0; k < 3; ++k) { GU_l[c][k] *= jdinv_l; GU_r[c][k] *= jdinv_r; } 
  } 
 
  for (int i = 0; i < 1; ++i) { 
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
      pkpm_harm_conf_flux(0, X_l, X_r, wq, &fc[i], 1, phi_l, phi_r, sb); 
      for (int c = 0; c < 8; ++c) Fhat[c][i*3 + j] = phi_l[c] + phi_r[c]; 
      alpha_max = fmax(alpha_max, fmax(-sb[0], sb[1])); 
    } 
  } 
 
  for (int c = 0; c < 8; ++c) { 
    for (int tid = 0; tid < 3; ++tid) pkpm_conf_flux_x_1x1v_ser_p1_prj(tid, Fhat[c], &flux[c*3]); 
  } 
  return 1.5*dx1*jdinv_r*alpha_max; 
} 
