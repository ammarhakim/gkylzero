#include <gkyl_pkpm_kinetic_kernels.h> 
#include <gkyl_pkpm_flux_surf_tables_1x1v_p1.h> 
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
  double t[3]; 
  for (int i = 0; i < 3; ++i) t[i] = 0.0; 
  for (int j = 0; j < 1; ++j) { 
    const double w = vst_pkpm_1x1v_p1_prj_v0_Vw[j*1 + b]; 
    for (int i = 0; i < 3; ++i) t[i] += w*Fhat_nodal[i*1 + j]; 
  } 
  for (int q = vst_pkpm_1x1v_p1_prj_v0_boff[b]; q < vst_pkpm_1x1v_p1_prj_v0_boff[b+1]; ++q) { 
    const int k = vst_pkpm_1x1v_p1_prj_v0_bks[q]; 
    const int a = vst_pkpm_1x1v_p1_prj_v0_kamap[k]; 
    double g = 0.0; 
    for (int i = 0; i < 3; ++i) g += vst_pkpm_1x1v_p1_prj_v0_Cw[i*2 + a]*t[i]; 
    flux[0 + k] = g; 
  } 
} 

GKYL_CU_DH double pkpm_vel_flux_vpar_1x1v_ser_p1(const double *w, const double *dxv,
  const double *bvar, const double *grad_b, const double *u, const double *grad_u, const double *theta,
  const double *a_P, const double *f_l, const double *f_r, double* GKYL_RESTRICT flux) 
{ 
  // w[NDIM]:   Cell-center coordinates (current cell). 
  // dxv[NDIM]: Cell spacing (current cell). 
  // bvar:      Magnetic field unit vector [bx, by, bz]. 
  // grad_b:    Weak gradient of b, [d_i bx, d_i by, d_i bz] per direction i. 
  // u:         Flow velocity [ux, uy, uz]. 
  // grad_u:    Weak gradient of u, [d_i ux, d_i uy, d_i uz] per direction i. 
  // theta:     T_perp/m. 
  // a_P:       Pressure force closing the first moment. 
  // f_l/f_r:   Distribution functions [F_0, G] in the cell below/current cell. 
  // flux:      Output modal surface fluxes [F_0 flux, G flux] on the lower v_par face. 
  // Returns the stable time-step estimate over the face nodes. 
  const double dx10 = 2.0/dxv[0]; 
  const double dv1 = 2.0/dxv[1]; 
  const double v_face = w[1] - 0.5*dxv[1]; 
  const double *bx = &bvar[0]; 
  const double *by = &bvar[2]; 
  const double *bz = &bvar[4]; 
  const double *ux = &u[0]; 
  const double *uy = &u[2]; 
  const double *uz = &u[4]; 
 
  double GF_l[2], GF_r[2], GG_l[2], GG_r[2]; 
  for (int tid = 0; tid < 2; ++tid) { 
    pkpm_vel_flux_vpar_1x1v_ser_p1_g(tid, &f_l[0], &f_r[0], GF_l, GF_r); 
    pkpm_vel_flux_vpar_1x1v_ser_p1_g(tid, &f_l[6], &f_r[6], GG_l, GG_r); 
  } 
 
  double Fhat_F[3], Fhat_G[3]; 
  double cfl_max = 0.0; 
  for (int i = 0; i < 3; ++i) { 
    // Configuration-space coefficients at node i. 
    double b_q[3] = {0.0}, u_q[3] = {0.0}, gradu_q[3][3] = {{0.0}}, div_b = 0.0, theta_q = 0.0, aP_q = 0.0; 
    for (int a = 0; a < 2; ++a) { 
      const double ev = vst_pkpm_1x1v_p1_conf_ev[i*2 + a]; 
      b_q[0] += ev*bx[a]; b_q[1] += ev*by[a]; b_q[2] += ev*bz[a]; 
      u_q[0] += ev*ux[a]; u_q[1] += ev*uy[a]; u_q[2] += ev*uz[a]; 
      theta_q += ev*theta[a]; 
      aP_q += ev*a_P[a]; 
      div_b += ev*grad_b[0 + a]; 
      gradu_q[0][0] += ev*grad_u[0 + a]; 
      gradu_q[0][1] += ev*grad_u[2 + a]; 
      gradu_q[0][2] += ev*grad_u[4 + a]; 
    } 
    double K = 0.0; 
    for (int d = 0; d < 1; ++d) { 
      for (int j = 0; j < 3; ++j) K += b_q[d]*b_q[j]*gradu_q[d][j]; 
    } 
    const double alpha = aP_q - v_face*K; 
    const double beta = div_b; 
    // Nodal traces of F and G on both sides of the face. 
    double F_l = 0.0, F_r = 0.0, G_l = 0.0, G_r = 0.0; 
    for (int a = 0; a < 2; ++a) { 
      const double cm = vst_pkpm_1x1v_p1_ph_v0_Cm[i*2 + a]; 
      F_l += cm*GF_l[a]; F_r += cm*GF_r[a]; G_l += cm*GG_l[a]; G_r += cm*GG_r[a]; 
    } 
    const double H_l = 4.0*theta_q*G_l - 2.0*theta_q*theta_q*F_l; 
    const double H_r = 4.0*theta_q*G_r - 2.0*theta_q*theta_q*F_r; 
    const double ap = fmax(alpha, 0.0), am = fmin(alpha, 0.0); 
    const double bp = fmax(beta, 0.0), bm = fmin(beta, 0.0); 
    Fhat_F[i] = ap*F_l + am*F_r + bp*G_l + bm*G_r; 
    Fhat_G[i] = ap*G_l + am*G_r + bp*H_l + bm*H_r; 
    cfl_max = fmax(cfl_max, fabs(alpha) + 3.414213562373095*fabs(theta_q*beta)); 
  } 
 
  for (int tid = 0; tid < 3; ++tid) { 
    pkpm_vel_flux_vpar_1x1v_ser_p1_prj(tid, Fhat_F, &flux[0]); 
    pkpm_vel_flux_vpar_1x1v_ser_p1_prj(tid, Fhat_G, &flux[2]); 
  } 
  return 2.5*dv1*cfl_max; 
} 
