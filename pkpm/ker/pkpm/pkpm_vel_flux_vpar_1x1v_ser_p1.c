#include <gkyl_pkpm_kinetic_kernels.h> 
#include <gkyl_pkpm_flux_surf_tables_1x1v_p1.h> 
#include <gkyl_pkpm_harm_coeffs.h> 
GKYL_CU_DH void pkpm_vel_flux_vpar_1x1v_ser_p1_g(int tid, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT G_l, double* GKYL_RESTRICT G_r) 
{ 
  if (tid >= 2) return; 
  const int a = tid/1; 
  const int j = tid - a*1; 
  double g_l = 0.0; 
  double g_r = 0.0; 
  for (int q = vst_pkpm_1x1v_p1_ph_v0_aoff[a]; q < vst_pkpm_1x1v_p1_ph_v0_aoff[a+1]; ++q) { 
    const int k = vst_pkpm_1x1v_p1_ph_v0_aks[q]; 
    g_l += vst_pkpm_1x1v_p1_ph_v0_Wl[k*1 + j]*f_l[k]; 
    g_r += vst_pkpm_1x1v_p1_ph_v0_Wr[k*1 + j]*f_r[k]; 
  } 
  G_l[j*2 + a] = g_l; 
  G_r[j*2 + a] = g_r; 
} 

GKYL_CU_DH void pkpm_vel_flux_vpar_1x1v_ser_p1_prj(int tid, const double *Fhat_nodal, double* GKYL_RESTRICT flux) 
{ 
  if (tid >= 1) return; 
  const int b = tid; 
  double t[5]; 
  for (int i = 0; i < 5; ++i) t[i] = 0.0; 
  for (int j = 0; j < 1; ++j) { 
    const double w = vst_pkpm_1x1v_p1_prj_v0_Vw[j*1 + b]; 
    for (int i = 0; i < 5; ++i) t[i] += w*Fhat_nodal[i*1 + j]; 
  } 
  for (int q = vst_pkpm_1x1v_p1_prj_v0_boff[b]; q < vst_pkpm_1x1v_p1_prj_v0_boff[b+1]; ++q) { 
    const int k = vst_pkpm_1x1v_p1_prj_v0_bks[q]; 
    const int a = vst_pkpm_1x1v_p1_prj_v0_kamap[k]; 
    double g = 0.0; 
    for (int i = 0; i < 5; ++i) g += vst_pkpm_1x1v_p1_prj_v0_Cw[i*2 + a]*t[i]; 
    flux[0 + k] = g; 
  } 
} 

GKYL_CU_DH double pkpm_vel_flux_vpar_1x1v_ser_p1(const double *w, const double *dxv, const double *nv,
  const double *f_l, const double *f_r, double* GKYL_RESTRICT flux) 
{ 
  // w[NDIM]:   Cell-center coordinates (current cell). 
  // dxv[NDIM]: Cell spacing (current cell). 
  // nv:        Node values [5 nodes][PKPM_HS_NODE_NV] of the w flux of the configuration cell. 
  // f_l/f_r:   Distributions [F, G, M0_perp, N_perp] in the cell below/current cell. 
  // flux:      Output modal surface fluxes [8 blocks] on the lower w face. 
  // Returns the stable time-step estimate over the face nodes. 
  const double dv1 = 2.0/dxv[1]; 
  const double v_face = w[1] - 0.5*dxv[1]; 
 
  double GU_l[8][2], GU_r[8][2]; 
  for (int tid = 0; tid < 2; ++tid) { 
    for (int c = 0; c < 8; ++c) { 
      pkpm_vel_flux_vpar_1x1v_ser_p1_g(tid, &f_l[c*6], &f_r[c*6], GU_l[c], GU_r[c]); 
    } 
  } 
 
  double Fhat[8][5]; 
  double cfl_max = 0.0; 
  for (int i = 0; i < 5; ++i) { 
    // Nodal traces of the eight blocks on both sides of the face. 
    double U_l[8] = {0.0}, U_r[8] = {0.0}, Fn[8]; 
    for (int c = 0; c < 8; ++c) { 
      for (int a = 0; a < 2; ++a) { 
        const double cm = vst_pkpm_1x1v_p1_ph_v0_Cm[i*2 + a]; 
        U_l[c] += cm*GU_l[c][a]; U_r[c] += cm*GU_r[c][a]; 
      } 
    } 
    cfl_max = fmax(cfl_max, pkpm_harm_node_wflux(&nv[i*PKPM_HS_NODE_NV], v_face, U_l, U_r, Fn)); 
    for (int c = 0; c < 8; ++c) Fhat[c][i] = Fn[c]; 
  } 
 
  for (int tid = 0; tid < 5; ++tid) { 
    for (int c = 0; c < 8; ++c) pkpm_vel_flux_vpar_1x1v_ser_p1_prj(tid, Fhat[c], &flux[c*2]); 
  } 
  return 2.5*dv1*cfl_max; 
} 
