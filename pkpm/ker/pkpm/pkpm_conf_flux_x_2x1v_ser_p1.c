#include <gkyl_pkpm_kinetic_kernels.h> 
#include <gkyl_pkpm_harm_coeffs.h> 
#include <gkyl_pkpm_flux_surf_tables_2x1v_p1.h> 
GKYL_CU_DH void pkpm_conf_flux_x_2x1v_ser_p1_g(int tid, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT G_l, double* GKYL_RESTRICT G_r) 
{ 
  if (tid >= 6) return; 
  const int a = tid/3; 
  const int j = tid - a*3; 
  double g_l = 0.0; 
  double g_r = 0.0; 
  for (int q = vst_pkpm_2x1v_p1_ph_x0_aoff[a]; q < vst_pkpm_2x1v_p1_ph_x0_aoff[a+1]; ++q) { 
    const int k = vst_pkpm_2x1v_p1_ph_x0_aks[q]; 
    g_l += vst_pkpm_2x1v_p1_ph_x0_Wl[k*3 + j]*f_l[k]; 
    g_r += vst_pkpm_2x1v_p1_ph_x0_Wr[k*3 + j]*f_r[k]; 
  } 
  G_l[j*2 + a] = g_l; 
  G_r[j*2 + a] = g_r; 
} 

GKYL_CU_DH void pkpm_conf_flux_x_2x1v_ser_p1_prj(int tid, const double *Fhat_nodal, double* GKYL_RESTRICT flux) 
{ 
  if (tid >= 3) return; 
  const int b = tid; 
  double t[5]; 
  for (int i = 0; i < 5; ++i) t[i] = 0.0; 
  for (int j = 0; j < 3; ++j) { 
    const double w = vst_pkpm_2x1v_p1_prj_x0_Vw[j*3 + b]; 
    for (int i = 0; i < 5; ++i) t[i] += w*Fhat_nodal[i*3 + j]; 
  } 
  for (int q = vst_pkpm_2x1v_p1_prj_x0_boff[b]; q < vst_pkpm_2x1v_p1_prj_x0_boff[b+1]; ++q) { 
    const int k = vst_pkpm_2x1v_p1_prj_x0_bks[q]; 
    const int a = vst_pkpm_2x1v_p1_prj_x0_kamap[k]; 
    double g = 0.0; 
    for (int i = 0; i < 5; ++i) g += vst_pkpm_2x1v_p1_prj_x0_Cw[i*2 + a]*t[i]; 
    flux[0 + k] = g; 
  } 
} 

GKYL_CU_DH double pkpm_conf_flux_x_2x1v_ser_p1(const double *w, const double *dxv,
  const double *jacob_pos_l, const double *jacob_pos_r, const double *b_l, const double *b_r,
  const double *u_l, const double *u_r, const double *xi_l, const double *xi_r,
  const double *theta_l, const double *theta_r, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT flux) 
{ 
  // w[NDIM]:   Cell-center coordinates (current cell, computational). 
  // dxv[NDIM]: Cell spacing (current cell, computational). 
  // jacob_pos_l/r: Position-map Jacobian per configuration direction (cdim blocks of 2 equal 
  //            entries) of the left/current cell; the total Jacobian J is their product. 
  // b_l/b_r:   Magnetic field unit vector [bx, by, bz] in the left/current cell. 
  // u_l/u_r:   Flow velocity [ux, uy, uz] in the left/current cell. 
  // xi_l/xi_r: Thermal speed [xi, xi^2] in the left/current cell (the first expansion is used). 
  // theta_l/r: T_perp/T in the left/current cell. 
  // f_l/f_r:   Jacobian-weighted normalized distributions J [F, G, M0_perp, N_perp] in the left/current cell. 
  // flux:      Output modal surface fluxes of the eight evolved blocks on the lower face of the 
  //            current cell in direction x: the physical HLL flux of the physical traces 
  //            (J f)/J times the face Jacobian J/J_x, i.e. the traces scaled by 1/J_x of their cell. 
  // Returns the stable time-step estimate from the largest speed bound over the nodes. 
  const double dx1 = 2.0/dxv[0]; 
  double Fhat[8][15]; 
  double alpha_max = 0.0; 
  const double jdinv_l = 1.0/jacob_pos_l[0], jdinv_r = 1.0/jacob_pos_r[0]; 
  const double wv = w[2], dv = dxv[2]; 
  // w of the nodes: wv + wn[j] dv 
  const double wn[3] = {-0.38729833462074165, 0.0, 0.38729833462074165}; 
  const double *b0_l = &b_l[0], *b0_r = &b_r[0]; 
  const double *b1_l = &b_l[4], *b1_r = &b_r[4]; 
  const double *b2_l = &b_l[8], *b2_r = &b_r[8]; 
  const double *ud_l = &u_l[0], *ud_r = &u_r[0]; 
  const double *xd_l = &xi_l[0], *xd_r = &xi_r[0]; 
 
  // Traces of b, u_x, xi and theta at the transverse nodes: the left cell at its upper 
  // edge, the current cell at its lower edge. 
  double bq0_l[5], bq0_r[5]; 
  bq0_l[0] = -(1.359269768907996*b0_l[3])-0.7847747669803519*b0_l[2]+0.8660254037844386*b0_l[1]+0.5*b0_l[0]; 
  bq0_l[1] = -(0.8077039651585247*b0_l[3])-0.4663281017098023*b0_l[2]+0.8660254037844386*b0_l[1]+0.5*b0_l[0]; 
  bq0_l[2] = 0.8660254037844386*b0_l[1]+0.5*b0_l[0]; 
  bq0_l[3] = 0.8077039651585247*b0_l[3]+0.4663281017098023*b0_l[2]+0.8660254037844386*b0_l[1]+0.5*b0_l[0]; 
  bq0_l[4] = 1.359269768907996*b0_l[3]+0.7847747669803519*b0_l[2]+0.8660254037844386*b0_l[1]+0.5*b0_l[0]; 
  bq0_r[0] = 1.359269768907996*b0_r[3]-0.7847747669803519*b0_r[2]-0.8660254037844386*b0_r[1]+0.5*b0_r[0]; 
  bq0_r[1] = 0.8077039651585247*b0_r[3]-0.4663281017098023*b0_r[2]-0.8660254037844386*b0_r[1]+0.5*b0_r[0]; 
  bq0_r[2] = 0.5*b0_r[0]-0.8660254037844386*b0_r[1]; 
  bq0_r[3] = -(0.8077039651585247*b0_r[3])+0.4663281017098023*b0_r[2]-0.8660254037844386*b0_r[1]+0.5*b0_r[0]; 
  bq0_r[4] = -(1.359269768907996*b0_r[3])+0.7847747669803519*b0_r[2]-0.8660254037844386*b0_r[1]+0.5*b0_r[0]; 
  double bq1_l[5], bq1_r[5]; 
  bq1_l[0] = -(1.359269768907996*b1_l[3])-0.7847747669803519*b1_l[2]+0.8660254037844386*b1_l[1]+0.5*b1_l[0]; 
  bq1_l[1] = -(0.8077039651585247*b1_l[3])-0.4663281017098023*b1_l[2]+0.8660254037844386*b1_l[1]+0.5*b1_l[0]; 
  bq1_l[2] = 0.8660254037844386*b1_l[1]+0.5*b1_l[0]; 
  bq1_l[3] = 0.8077039651585247*b1_l[3]+0.4663281017098023*b1_l[2]+0.8660254037844386*b1_l[1]+0.5*b1_l[0]; 
  bq1_l[4] = 1.359269768907996*b1_l[3]+0.7847747669803519*b1_l[2]+0.8660254037844386*b1_l[1]+0.5*b1_l[0]; 
  bq1_r[0] = 1.359269768907996*b1_r[3]-0.7847747669803519*b1_r[2]-0.8660254037844386*b1_r[1]+0.5*b1_r[0]; 
  bq1_r[1] = 0.8077039651585247*b1_r[3]-0.4663281017098023*b1_r[2]-0.8660254037844386*b1_r[1]+0.5*b1_r[0]; 
  bq1_r[2] = 0.5*b1_r[0]-0.8660254037844386*b1_r[1]; 
  bq1_r[3] = -(0.8077039651585247*b1_r[3])+0.4663281017098023*b1_r[2]-0.8660254037844386*b1_r[1]+0.5*b1_r[0]; 
  bq1_r[4] = -(1.359269768907996*b1_r[3])+0.7847747669803519*b1_r[2]-0.8660254037844386*b1_r[1]+0.5*b1_r[0]; 
  double bq2_l[5], bq2_r[5]; 
  bq2_l[0] = -(1.359269768907996*b2_l[3])-0.7847747669803519*b2_l[2]+0.8660254037844386*b2_l[1]+0.5*b2_l[0]; 
  bq2_l[1] = -(0.8077039651585247*b2_l[3])-0.4663281017098023*b2_l[2]+0.8660254037844386*b2_l[1]+0.5*b2_l[0]; 
  bq2_l[2] = 0.8660254037844386*b2_l[1]+0.5*b2_l[0]; 
  bq2_l[3] = 0.8077039651585247*b2_l[3]+0.4663281017098023*b2_l[2]+0.8660254037844386*b2_l[1]+0.5*b2_l[0]; 
  bq2_l[4] = 1.359269768907996*b2_l[3]+0.7847747669803519*b2_l[2]+0.8660254037844386*b2_l[1]+0.5*b2_l[0]; 
  bq2_r[0] = 1.359269768907996*b2_r[3]-0.7847747669803519*b2_r[2]-0.8660254037844386*b2_r[1]+0.5*b2_r[0]; 
  bq2_r[1] = 0.8077039651585247*b2_r[3]-0.4663281017098023*b2_r[2]-0.8660254037844386*b2_r[1]+0.5*b2_r[0]; 
  bq2_r[2] = 0.5*b2_r[0]-0.8660254037844386*b2_r[1]; 
  bq2_r[3] = -(0.8077039651585247*b2_r[3])+0.4663281017098023*b2_r[2]-0.8660254037844386*b2_r[1]+0.5*b2_r[0]; 
  bq2_r[4] = -(1.359269768907996*b2_r[3])+0.7847747669803519*b2_r[2]-0.8660254037844386*b2_r[1]+0.5*b2_r[0]; 
  double uq_l[5], uq_r[5], xq_l[5], xq_r[5], thq_l[5], thq_r[5]; 
  uq_l[0] = -(1.359269768907996*ud_l[3])-0.7847747669803519*ud_l[2]+0.8660254037844386*ud_l[1]+0.5*ud_l[0]; 
  uq_l[1] = -(0.8077039651585247*ud_l[3])-0.4663281017098023*ud_l[2]+0.8660254037844386*ud_l[1]+0.5*ud_l[0]; 
  uq_l[2] = 0.8660254037844386*ud_l[1]+0.5*ud_l[0]; 
  uq_l[3] = 0.8077039651585247*ud_l[3]+0.4663281017098023*ud_l[2]+0.8660254037844386*ud_l[1]+0.5*ud_l[0]; 
  uq_l[4] = 1.359269768907996*ud_l[3]+0.7847747669803519*ud_l[2]+0.8660254037844386*ud_l[1]+0.5*ud_l[0]; 
  uq_r[0] = 1.359269768907996*ud_r[3]-0.7847747669803519*ud_r[2]-0.8660254037844386*ud_r[1]+0.5*ud_r[0]; 
  uq_r[1] = 0.8077039651585247*ud_r[3]-0.4663281017098023*ud_r[2]-0.8660254037844386*ud_r[1]+0.5*ud_r[0]; 
  uq_r[2] = 0.5*ud_r[0]-0.8660254037844386*ud_r[1]; 
  uq_r[3] = -(0.8077039651585247*ud_r[3])+0.4663281017098023*ud_r[2]-0.8660254037844386*ud_r[1]+0.5*ud_r[0]; 
  uq_r[4] = -(1.359269768907996*ud_r[3])+0.7847747669803519*ud_r[2]-0.8660254037844386*ud_r[1]+0.5*ud_r[0]; 
  xq_l[0] = -(1.359269768907996*xd_l[3])-0.7847747669803519*xd_l[2]+0.8660254037844386*xd_l[1]+0.5*xd_l[0]; 
  xq_l[1] = -(0.8077039651585247*xd_l[3])-0.4663281017098023*xd_l[2]+0.8660254037844386*xd_l[1]+0.5*xd_l[0]; 
  xq_l[2] = 0.8660254037844386*xd_l[1]+0.5*xd_l[0]; 
  xq_l[3] = 0.8077039651585247*xd_l[3]+0.4663281017098023*xd_l[2]+0.8660254037844386*xd_l[1]+0.5*xd_l[0]; 
  xq_l[4] = 1.359269768907996*xd_l[3]+0.7847747669803519*xd_l[2]+0.8660254037844386*xd_l[1]+0.5*xd_l[0]; 
  xq_r[0] = 1.359269768907996*xd_r[3]-0.7847747669803519*xd_r[2]-0.8660254037844386*xd_r[1]+0.5*xd_r[0]; 
  xq_r[1] = 0.8077039651585247*xd_r[3]-0.4663281017098023*xd_r[2]-0.8660254037844386*xd_r[1]+0.5*xd_r[0]; 
  xq_r[2] = 0.5*xd_r[0]-0.8660254037844386*xd_r[1]; 
  xq_r[3] = -(0.8077039651585247*xd_r[3])+0.4663281017098023*xd_r[2]-0.8660254037844386*xd_r[1]+0.5*xd_r[0]; 
  xq_r[4] = -(1.359269768907996*xd_r[3])+0.7847747669803519*xd_r[2]-0.8660254037844386*xd_r[1]+0.5*xd_r[0]; 
  thq_l[0] = -(1.359269768907996*theta_l[3])-0.7847747669803519*theta_l[2]+0.8660254037844386*theta_l[1]+0.5*theta_l[0]; 
  thq_l[1] = -(0.8077039651585247*theta_l[3])-0.4663281017098023*theta_l[2]+0.8660254037844386*theta_l[1]+0.5*theta_l[0]; 
  thq_l[2] = 0.8660254037844386*theta_l[1]+0.5*theta_l[0]; 
  thq_l[3] = 0.8077039651585247*theta_l[3]+0.4663281017098023*theta_l[2]+0.8660254037844386*theta_l[1]+0.5*theta_l[0]; 
  thq_l[4] = 1.359269768907996*theta_l[3]+0.7847747669803519*theta_l[2]+0.8660254037844386*theta_l[1]+0.5*theta_l[0]; 
  thq_r[0] = 1.359269768907996*theta_r[3]-0.7847747669803519*theta_r[2]-0.8660254037844386*theta_r[1]+0.5*theta_r[0]; 
  thq_r[1] = 0.8077039651585247*theta_r[3]-0.4663281017098023*theta_r[2]-0.8660254037844386*theta_r[1]+0.5*theta_r[0]; 
  thq_r[2] = 0.5*theta_r[0]-0.8660254037844386*theta_r[1]; 
  thq_r[3] = -(0.8077039651585247*theta_r[3])+0.4663281017098023*theta_r[2]-0.8660254037844386*theta_r[1]+0.5*theta_r[0]; 
  thq_r[4] = -(1.359269768907996*theta_r[3])+0.7847747669803519*theta_r[2]-0.8660254037844386*theta_r[1]+0.5*theta_r[0]; 
 
  // Trace factors of the eight blocks, scaled by 1/J_x of their cell. 
  double GU_l[8][6], GU_r[8][6]; 
  for (int c = 0; c < 8; ++c) { 
    for (int tid = 0; tid < 6; ++tid) pkpm_conf_flux_x_2x1v_ser_p1_g(tid, &f_l[c*12], &f_r[c*12], GU_l[c], GU_r[c]); 
    for (int k = 0; k < 6; ++k) { GU_l[c][k] *= jdinv_l; GU_r[c][k] *= jdinv_r; } 
  } 
 
  for (int i = 0; i < 5; ++i) { 
    const double bn_l[3] = {bq0_l[i], bq1_l[i], bq2_l[i]}, bn_r[3] = {bq0_r[i], bq1_r[i], bq2_r[i]}; 
    for (int j = 0; j < 3; ++j) { 
      const double wq = wv + wn[j]*dv; 
      double X_l[8], X_r[8]; 
      for (int c = 0; c < 8; ++c) { 
        double x_l = 0.0, x_r = 0.0; 
        for (int a = 0; a < 2; ++a) { 
          x_l += vst_pkpm_2x1v_p1_ph_x0_Cm[i*2 + a]*GU_l[c][j*2 + a]; 
          x_r += vst_pkpm_2x1v_p1_ph_x0_Cm[i*2 + a]*GU_r[c][j*2 + a]; 
        } 
        X_l[c] = x_l; X_r[c] = x_r; 
      } 
      double phi_l[8], phi_r[8], sb[2]; 
      pkpm_harm_conf_flux(0, X_l, X_r, xq_l[i]*bn_l[0]*wq + uq_l[i], xq_r[i]*bn_r[0]*wq + uq_r[i], 
        xq_l[i], xq_r[i], bn_l, bn_r, thq_l[i], thq_r[i], phi_l, phi_r, sb); 
      for (int c = 0; c < 8; ++c) Fhat[c][i*3 + j] = phi_l[c] + phi_r[c]; 
      alpha_max = fmax(alpha_max, fmax(-sb[0], sb[1])); 
    } 
  } 
 
  for (int c = 0; c < 8; ++c) { 
    for (int tid = 0; tid < 3; ++tid) pkpm_conf_flux_x_2x1v_ser_p1_prj(tid, Fhat[c], &flux[c*6]); 
  } 
  return 1.5*dx1*jdinv_r*alpha_max; 
} 
