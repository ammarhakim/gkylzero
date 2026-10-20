#include <gkyl_pkpm_kinetic_kernels.h> 
#include <gkyl_pkpm_flux_surf_tables_2x1v_p1.h> 
#include <gkyl_pkpm_harm_coeffs.h> 
GKYL_CU_DH void pkpm_vel_flux_vpar_2x1v_ser_p1_g(int tid, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT G_l, double* GKYL_RESTRICT G_r) 
{ 
  if (tid >= 4) return; 
  const int a = tid/1; 
  const int j = tid - a*1; 
  double g_l = 0.0; 
  double g_r = 0.0; 
  for (int q = vst_pkpm_2x1v_p1_ph_v0_aoff[a]; q < vst_pkpm_2x1v_p1_ph_v0_aoff[a+1]; ++q) { 
    const int k = vst_pkpm_2x1v_p1_ph_v0_aks[q]; 
    g_l += vst_pkpm_2x1v_p1_ph_v0_Wl[k*1 + j]*f_l[k]; 
    g_r += vst_pkpm_2x1v_p1_ph_v0_Wr[k*1 + j]*f_r[k]; 
  } 
  G_l[j*4 + a] = g_l; 
  G_r[j*4 + a] = g_r; 
} 

GKYL_CU_DH void pkpm_vel_flux_vpar_2x1v_ser_p1_prj(int tid, const double *Fhat_nodal, double* GKYL_RESTRICT flux) 
{ 
  if (tid >= 1) return; 
  const int b = tid; 
  double t[25]; 
  for (int i = 0; i < 25; ++i) t[i] = 0.0; 
  for (int j = 0; j < 1; ++j) { 
    const double w = vst_pkpm_2x1v_p1_prj_v0_Vw[j*1 + b]; 
    for (int i = 0; i < 25; ++i) t[i] += w*Fhat_nodal[i*1 + j]; 
  } 
  for (int q = vst_pkpm_2x1v_p1_prj_v0_boff[b]; q < vst_pkpm_2x1v_p1_prj_v0_boff[b+1]; ++q) { 
    const int k = vst_pkpm_2x1v_p1_prj_v0_bks[q]; 
    const int a = vst_pkpm_2x1v_p1_prj_v0_kamap[k]; 
    double g = 0.0; 
    for (int i = 0; i < 25; ++i) g += vst_pkpm_2x1v_p1_prj_v0_Cw[i*4 + a]*t[i]; 
    flux[0 + k] = g; 
  } 
} 

GKYL_CU_DH double pkpm_vel_flux_vpar_2x1v_ser_p1(const double *w, const double *dxv, const double *nv,
  const double *f_l, const double *f_r, double* GKYL_RESTRICT flux) 
{ 
  // w[NDIM]:   Cell-center coordinates (current cell). 
  // dxv[NDIM]: Cell spacing (current cell). 
  // nv:        Node values [25 nodes][PKPM_HS_NODE_NV] of the w flux of the configuration cell. 
  // f_l/f_r:   Distributions [F, G, M0_perp, N_perp] in the cell below/current cell. 
  // flux:      Output modal surface fluxes [8 blocks] on the lower w face. 
  // Returns the stable time-step estimate over the face nodes. 
  const double dv1 = 2.0/dxv[2]; 
  const double v_face = w[2] - 0.5*dxv[2]; 
 
  double GU_l[8][4], GU_r[8][4]; 
  for (int tid = 0; tid < 4; ++tid) { 
    for (int c = 0; c < 8; ++c) { 
      pkpm_vel_flux_vpar_2x1v_ser_p1_g(tid, &f_l[c*12], &f_r[c*12], GU_l[c], GU_r[c]); 
    } 
  } 
 
  double Fhat[8][25]; 
  double cfl_max = 0.0; 
  for (int i = 0; i < 25; ++i) { 
    // Nodal traces of the eight blocks on both sides of the face. 
    double U_l[8] = {0.0}, U_r[8] = {0.0}, Fn[8]; 
    for (int c = 0; c < 8; ++c) { 
      for (int a = 0; a < 4; ++a) { 
        const double cm = vst_pkpm_2x1v_p1_ph_v0_Cm[i*4 + a]; 
        U_l[c] += cm*GU_l[c][a]; U_r[c] += cm*GU_r[c][a]; 
      } 
    } 
    cfl_max = fmax(cfl_max, pkpm_harm_node_wflux(&nv[i*PKPM_HS_NODE_NV], v_face, U_l, U_r, Fn)); 
    for (int c = 0; c < 8; ++c) Fhat[c][i] = Fn[c]; 
  } 
 
  for (int tid = 0; tid < 25; ++tid) { 
    for (int c = 0; c < 8; ++c) pkpm_vel_flux_vpar_2x1v_ser_p1_prj(tid, Fhat[c], &flux[c*4]); 
  } 
  return 2.5*dv1*cfl_max; 
} 
