#include <gkyl_pkpm_kinetic_kernels.h> 
#include <gkyl_pkpm_flux_surf_tables_2x1v_p1.h> 
GKYL_CU_DH void pkpm_conf_flux_y_2x1v_ser_p1_g(int tid, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT G_l, double* GKYL_RESTRICT G_r) 
{ 
  if (tid >= 6) return; 
  const int a = tid/3; 
  const int j = tid - a*3; 
  double g_l = 0.0; 
  double g_r = 0.0; 
  for (int q = vst_pkpm_2x1v_p1_ph_x1_aoff[a]; q < vst_pkpm_2x1v_p1_ph_x1_aoff[a+1]; ++q) { 
    const int k = vst_pkpm_2x1v_p1_ph_x1_aks[q]; 
    g_l += vst_pkpm_2x1v_p1_ph_x1_Wl[k*3 + j]*f_l[k]; 
    g_r += vst_pkpm_2x1v_p1_ph_x1_Wr[k*3 + j]*f_r[k]; 
  } 
  G_l[j*2 + a] = g_l; 
  G_r[j*2 + a] = g_r; 
} 

GKYL_CU_DH double pkpm_conf_flux_y_2x1v_ser_p1_node(int i, int j, double a_l, double a_r,
  const double *G_l, const double *G_r, double* GKYL_RESTRICT Fhat_nodal) 
{ 
  double f_l_quad = 0.0; 
  double f_r_quad = 0.0; 
  for (int a = 0; a < 2; ++a) { 
    f_l_quad += vst_pkpm_2x1v_p1_ph_x1_Cm[i*2 + a]*G_l[j*2 + a]; 
    f_r_quad += vst_pkpm_2x1v_p1_ph_x1_Cm[i*2 + a]*G_r[j*2 + a]; 
  } 
  const int n = i*3 + j; 
  // HLL with S_L = min(a_l, a_r, 0), S_R = max(a_l, a_r, 0). 
  const double s_l = fmin(fmin(a_l, a_r), 0.0); 
  const double s_r = fmax(fmax(a_l, a_r), 0.0); 
  const double ds = s_r - s_l; 
  Fhat_nodal[n] = ds > 0.0 ? (s_r*a_l*f_l_quad - s_l*a_r*f_r_quad + s_l*s_r*(f_r_quad - f_l_quad))/ds : 0.0; 
  return fmax(-s_l, s_r); 
} 

GKYL_CU_DH void pkpm_conf_flux_y_2x1v_ser_p1_prj(int tid, const double *Fhat_nodal, double* GKYL_RESTRICT flux) 
{ 
  if (tid >= 3) return; 
  const int b = tid; 
  double t[3]; 
  for (int i = 0; i < 3; ++i) t[i] = 0.0; 
  for (int j = 0; j < 3; ++j) { 
    const double w = vst_pkpm_2x1v_p1_prj_x1_Vw[j*3 + b]; 
    for (int i = 0; i < 3; ++i) t[i] += w*Fhat_nodal[i*3 + j]; 
  } 
  for (int q = vst_pkpm_2x1v_p1_prj_x1_boff[b]; q < vst_pkpm_2x1v_p1_prj_x1_boff[b+1]; ++q) { 
    const int k = vst_pkpm_2x1v_p1_prj_x1_bks[q]; 
    const int a = vst_pkpm_2x1v_p1_prj_x1_kamap[k]; 
    double g = 0.0; 
    for (int i = 0; i < 3; ++i) g += vst_pkpm_2x1v_p1_prj_x1_Cw[i*2 + a]*t[i]; 
    flux[0 + k] = g; 
  } 
} 

GKYL_CU_DH double pkpm_conf_flux_y_2x1v_ser_p1(const double *w, const double *dxv,
  const double *jacob_pos_l, const double *jacob_pos_r,
  const double *b_l, const double *b_r, const double *u_l, const double *u_r,
  const double *xi_l, const double *xi_r, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT flux) 
{ 
  // w[NDIM]:   Cell-center coordinates (current cell, computational). 
  // dxv[NDIM]: Cell spacing (current cell, computational). 
  // jacob_pos_l/r: Position-map Jacobian per configuration direction (cdim blocks of 2 equal 
  //            entries) of the left/current cell; the total Jacobian J is their product. 
  // b_l/b_r:   Magnetic field unit vector [bx, by, bz] in the left/current cell. 
  // u_l/u_r:   Flow velocity [ux, uy, uz] in the left/current cell. 
  // xi_l/xi_r: Thermal speed [xi, xi^2] in the left/current cell (the first expansion is used). 
  // f_l/f_r:   Jacobian-weighted normalized distribution functions J [F, G] in the left/current cell. 
  // flux:      Output modal surface fluxes of the evolved J [F, G] on the lower face of the 
  //            current cell in direction y: the physical HLL flux of the physical traces 
  //            (J f)/J times the face Jacobian J/J_y, i.e. the traces scaled by 1/J_y of their cell. 
  // Returns the stable time-step estimate from the largest speed bound over the nodes. 
  const double dx1 = 2.0/dxv[1]; 
  const double jdinv_l = 1.0/jacob_pos_l[2], jdinv_r = 1.0/jacob_pos_r[2]; 
  const double wv = w[2], dv = dxv[2]; 
  const double *bd_l = &b_l[4]; 
  const double *bd_r = &b_r[4]; 
  const double *ud_l = &u_l[4]; 
  const double *ud_r = &u_r[4]; 
  const double *xd_l = &xi_l[0]; 
  const double *xd_r = &xi_r[0]; 
 
  double bl_q[3] = {0.0}; 
  double ul_q[3] = {0.0}; 
  double br_q[3] = {0.0}; 
  double ur_q[3] = {0.0}; 
  double xl_q[3] = {0.0}; 
  double xr_q[3] = {0.0}; 
  xl_q[0] = -(1.1618950038622249*xd_l[3])+0.8660254037844386*xd_l[2]-0.6708203932499369*xd_l[1]+0.5*xd_l[0]; 
  xl_q[1] = 0.8660254037844386*xd_l[2]+0.5*xd_l[0]; 
  xl_q[2] = 1.1618950038622249*xd_l[3]+0.8660254037844386*xd_l[2]+0.6708203932499369*xd_l[1]+0.5*xd_l[0]; 
  xr_q[0] = 1.1618950038622249*xd_r[3]-0.8660254037844386*xd_r[2]-0.6708203932499369*xd_r[1]+0.5*xd_r[0]; 
  xr_q[1] = 0.5*xd_r[0]-0.8660254037844386*xd_r[2]; 
  xr_q[2] = -(1.1618950038622249*xd_r[3])-0.8660254037844386*xd_r[2]+0.6708203932499369*xd_r[1]+0.5*xd_r[0]; 
  bl_q[0] = -(1.1618950038622249*bd_l[3])+0.8660254037844386*bd_l[2]-0.6708203932499369*bd_l[1]+0.5*bd_l[0]; 
  bl_q[1] = 0.8660254037844386*bd_l[2]+0.5*bd_l[0]; 
  bl_q[2] = 1.1618950038622249*bd_l[3]+0.8660254037844386*bd_l[2]+0.6708203932499369*bd_l[1]+0.5*bd_l[0]; 
  ul_q[0] = -(1.1618950038622249*ud_l[3])+0.8660254037844386*ud_l[2]-0.6708203932499369*ud_l[1]+0.5*ud_l[0]; 
  ul_q[1] = 0.8660254037844386*ud_l[2]+0.5*ud_l[0]; 
  ul_q[2] = 1.1618950038622249*ud_l[3]+0.8660254037844386*ud_l[2]+0.6708203932499369*ud_l[1]+0.5*ud_l[0]; 
  br_q[0] = 1.1618950038622249*bd_r[3]-0.8660254037844386*bd_r[2]-0.6708203932499369*bd_r[1]+0.5*bd_r[0]; 
  br_q[1] = 0.5*bd_r[0]-0.8660254037844386*bd_r[2]; 
  br_q[2] = -(1.1618950038622249*bd_r[3])-0.8660254037844386*bd_r[2]+0.6708203932499369*bd_r[1]+0.5*bd_r[0]; 
  ur_q[0] = 1.1618950038622249*ud_r[3]-0.8660254037844386*ud_r[2]-0.6708203932499369*ud_r[1]+0.5*ud_r[0]; 
  ur_q[1] = 0.5*ud_r[0]-0.8660254037844386*ud_r[2]; 
  ur_q[2] = -(1.1618950038622249*ud_r[3])-0.8660254037844386*ud_r[2]+0.6708203932499369*ud_r[1]+0.5*ud_r[0]; 
 
  double a_l[9] = {0.0}; 
  double a_r[9] = {0.0}; 
  for (int i = 0; i < 3; ++i) { 
    a_l[i*3 + 0] = xl_q[i]*bl_q[i]*(wv + -0.38729833462074165*dv) + ul_q[i]; 
    a_r[i*3 + 0] = xr_q[i]*br_q[i]*(wv + -0.38729833462074165*dv) + ur_q[i]; 
    a_l[i*3 + 1] = xl_q[i]*bl_q[i]*(wv + 0.0*dv) + ul_q[i]; 
    a_r[i*3 + 1] = xr_q[i]*br_q[i]*(wv + 0.0*dv) + ur_q[i]; 
    a_l[i*3 + 2] = xl_q[i]*bl_q[i]*(wv + 0.38729833462074165*dv) + ul_q[i]; 
    a_r[i*3 + 2] = xr_q[i]*br_q[i]*(wv + 0.38729833462074165*dv) + ur_q[i]; 
  } 
 
  double G_l[6]; 
  double G_r[6]; 
  double Fhat_nodal[9]; 
  double alpha_max = 0.0; 
  for (int c = 0; c < 2; ++c) { 
    for (int tid = 0; tid < 6; ++tid) pkpm_conf_flux_y_2x1v_ser_p1_g(tid, &f_l[c*12], &f_r[c*12], G_l, G_r); 
    for (int k = 0; k < 6; ++k) { G_l[k] *= jdinv_l; G_r[k] *= jdinv_r; } 
    for (int i = 0; i < 3; ++i) { 
      for (int j = 0; j < 3; ++j) { 
        alpha_max = fmax(alpha_max, pkpm_conf_flux_y_2x1v_ser_p1_node(i, j, a_l[i*3 + j], a_r[i*3 + j], G_l, G_r, Fhat_nodal)); 
      } 
    } 
    for (int tid = 0; tid < 9; ++tid) pkpm_conf_flux_y_2x1v_ser_p1_prj(tid, Fhat_nodal, &flux[c*6]); 
  } 
  return 1.5*dx1*jdinv_r*alpha_max; 
} 
