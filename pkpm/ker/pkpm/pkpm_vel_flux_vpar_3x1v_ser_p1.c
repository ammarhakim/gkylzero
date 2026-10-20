#include <gkyl_pkpm_kinetic_kernels.h> 
#include <gkyl_pkpm_flux_surf_tables_3x1v_p1.h> 
#include <gkyl_pkpm_harm_coeffs.h> 
GKYL_CU_DH void pkpm_vel_flux_vpar_3x1v_ser_p1_g(int tid, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT G_l, double* GKYL_RESTRICT G_r) 
{ 
  if (tid >= 8) return; 
  const int a = tid/1; 
  const int j = tid - a*1; 
  double g_l = 0.0; 
  double g_r = 0.0; 
  for (int q = vst_pkpm_3x1v_p1_ph_v0_aoff[a]; q < vst_pkpm_3x1v_p1_ph_v0_aoff[a+1]; ++q) { 
    const int k = vst_pkpm_3x1v_p1_ph_v0_aks[q]; 
    g_l += vst_pkpm_3x1v_p1_ph_v0_Wl[k*1 + j]*f_l[k]; 
    g_r += vst_pkpm_3x1v_p1_ph_v0_Wr[k*1 + j]*f_r[k]; 
  } 
  G_l[j*8 + a] = g_l; 
  G_r[j*8 + a] = g_r; 
} 

GKYL_CU_DH void pkpm_vel_flux_vpar_3x1v_ser_p1_prj(int tid, const double *Fhat_nodal, double* GKYL_RESTRICT flux) 
{ 
  if (tid >= 1) return; 
  const int b = tid; 
  double t[125]; 
  for (int i = 0; i < 125; ++i) t[i] = 0.0; 
  for (int j = 0; j < 1; ++j) { 
    const double w = vst_pkpm_3x1v_p1_prj_v0_Vw[j*1 + b]; 
    for (int i = 0; i < 125; ++i) t[i] += w*Fhat_nodal[i*1 + j]; 
  } 
  for (int q = vst_pkpm_3x1v_p1_prj_v0_boff[b]; q < vst_pkpm_3x1v_p1_prj_v0_boff[b+1]; ++q) { 
    const int k = vst_pkpm_3x1v_p1_prj_v0_bks[q]; 
    const int a = vst_pkpm_3x1v_p1_prj_v0_kamap[k]; 
    double g = 0.0; 
    for (int i = 0; i < 125; ++i) g += vst_pkpm_3x1v_p1_prj_v0_Cw[i*8 + a]*t[i]; 
    flux[0 + k] = g; 
  } 
} 

GKYL_CU_DH double pkpm_vel_flux_vpar_3x1v_ser_p1(const double *w, const double *dxv,
  const double *bvar, const double *grad_b, const double *u, const double *grad_u,
  const double *xi, const double *grad_xi, const double *theta, const double *dtB_over_magB,
  const double *alpha_0, const double *c_h, const double *f_l, const double *f_r, double* GKYL_RESTRICT flux) 
{ 
  // w[NDIM]:       Cell-center coordinates (current cell). 
  // dxv[NDIM]:     Cell spacing (current cell). 
  // bvar:          Magnetic field unit vector [bx, by, bz]. 
  // grad_b:        Weak gradient of b, [d_i bx, d_i by, d_i bz] per direction i. 
  // u:             Flow velocity [ux, uy, uz]. 
  // grad_u:        Weak gradient of u, [d_i ux, d_i uy, d_i uz] per direction i. 
  // xi:            Thermal speed [xi, xi^2] (the first expansion is used). 
  // grad_xi:       Weak gradient of [xi, 0, 0] (pkpm_vars_grad_u_* layout). 
  // theta:         T_perp/T. 
  // dtB_over_magB: d_t B/|B| of the stage. 
  // alpha_0:       Pressure force a_P/xi closing the first moment. 
  // c_h:           Heating rate closing the scalar second moment. 
  // f_l/f_r:       Distributions [F, G, M0_perp, N_perp] in the cell below/current cell. 
  // flux:          Output modal surface fluxes [8 blocks] on the lower w face. 
  // Returns the stable time-step estimate over the face nodes. 
  const double dv1 = 2.0/dxv[3]; 
  const double v_face = w[3] - 0.5*dxv[3]; 
  const double *bx = &bvar[0]; 
  const double *by = &bvar[8]; 
  const double *bz = &bvar[16]; 
  const double *ux = &u[0]; 
  const double *uy = &u[8]; 
  const double *uz = &u[16]; 
 
  double GU_l[8][8], GU_r[8][8]; 
  for (int tid = 0; tid < 8; ++tid) { 
    for (int c = 0; c < 8; ++c) { 
      pkpm_vel_flux_vpar_3x1v_ser_p1_g(tid, &f_l[c*24], &f_r[c*24], GU_l[c], GU_r[c]); 
    } 
  } 
 
  double Fhat[8][125]; 
  double cfl_max = 0.0; 
  for (int i = 0; i < 125; ++i) { 
    // Configuration-space coefficients at node i. 
    double b_q[3] = {0.0}, u_q[3] = {0.0}, d_q[3] = {0.0}, theta_q = 0.0, a0_q = 0.0, ch_q = 0.0, xi_q = 0.0; 
    double gu_q[9] = {0.0}, gb_q[9] = {0.0}, gx_q[3] = {0.0}; 
    for (int a = 0; a < 8; ++a) { 
      const double ev = vst_pkpm_3x1v_p1_conf_ev[i*8 + a]; 
      b_q[0] += ev*bx[a]; b_q[1] += ev*by[a]; b_q[2] += ev*bz[a]; 
      u_q[0] += ev*ux[a]; u_q[1] += ev*uy[a]; u_q[2] += ev*uz[a]; 
      for (int j = 0; j < 3; ++j) d_q[j] += ev*dtB_over_magB[j*8 + a]; 
      theta_q += ev*theta[a]; 
      a0_q += ev*alpha_0[a]; 
      ch_q += ev*c_h[a]; 
      xi_q += ev*xi[a]; 
      for (int d = 0; d < 3; ++d) { 
        gx_q[d] += ev*grad_xi[(d*3)*8 + a]; 
        for (int j = 0; j < 3; ++j) { 
          gu_q[d*3 + j] += ev*grad_u[(d*3 + j)*8 + a]; 
          gb_q[d*3 + j] += ev*grad_b[(d*3 + j)*8 + a]; 
        } 
      } 
    } 
    struct pkpm_harm_coeffs hc; 
    pkpm_harm_coeffs_set(&hc, 3, b_q, u_q, xi_q, theta_q, d_q, gu_q, gb_q, gx_q, 0, 0.0, 0); 
    const double alpha = a0_q + (ch_q - hc.K)*v_face - hc.bg*v_face*v_face; 
    // Nodal traces of the eight blocks on both sides of the face. 
    double U_l[8] = {0.0}, U_r[8] = {0.0}; 
    for (int c = 0; c < 8; ++c) { 
      for (int a = 0; a < 8; ++a) { 
        const double cm = vst_pkpm_3x1v_p1_ph_v0_Cm[i*8 + a]; 
        U_l[c] += cm*GU_l[c][a]; U_r[c] += cm*GU_r[c][a]; 
      } 
    } 
    double T_l[8], T_r[8], lam_lo, lam_hi; 
    pkpm_harm_coupling_flux(&hc, v_face, U_l, T_l); 
    pkpm_harm_coupling_flux(&hc, v_face, U_r, T_r); 
    pkpm_harm_speed_bounds(&hc, v_face, &lam_lo, &lam_hi); 
    const double sl = alpha + lam_lo, sr = alpha + lam_hi; 
    if (sl >= 0.0) { 
      for (int c = 0; c < 8; ++c) Fhat[c][i] = alpha*U_l[c] + T_l[c]; 
    } else if (sr <= 0.0) { 
      for (int c = 0; c < 8; ++c) Fhat[c][i] = alpha*U_r[c] + T_r[c]; 
    } else { 
      for (int c = 0; c < 8; ++c) 
        Fhat[c][i] = (sr*(alpha*U_l[c] + T_l[c]) - sl*(alpha*U_r[c] + T_r[c]) + sl*sr*(U_r[c] - U_l[c]))/(sr - sl); 
    } 
    cfl_max = fmax(cfl_max, fmax(fabs(sl), fabs(sr))); 
  } 
 
  for (int tid = 0; tid < 125; ++tid) { 
    for (int c = 0; c < 8; ++c) pkpm_vel_flux_vpar_3x1v_ser_p1_prj(tid, Fhat[c], &flux[c*8]); 
  } 
  return 2.5*dv1*cfl_max; 
} 
