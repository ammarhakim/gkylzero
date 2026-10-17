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

GKYL_CU_DH double pkpm_conf_flux_y_2x1v_ser_p1_node(int i, int j, double a_l, double a_r, double c_l, double c_r,
  const double *G_l, const double *G_r, double* GKYL_RESTRICT Fhat_nodal) 
{ 
  double f_l_quad = 0.0; 
  double f_r_quad = 0.0; 
  for (int a = 0; a < 2; ++a) { 
    f_l_quad += vst_pkpm_2x1v_p1_ph_x1_Cm[i*2 + a]*G_l[j*2 + a]; 
    f_r_quad += vst_pkpm_2x1v_p1_ph_x1_Cm[i*2 + a]*G_r[j*2 + a]; 
  } 
  const int n = i*3 + j; 
  // HLL with S_L = min(a_l - c_l, a_r - c_r, 0), S_R = max(a_l + c_l, a_r + c_r, 0), where c is 
  // the acoustic speed of the perpendicular pressure, which the streaming speed does not carry. 
  const double s_l = fmin(fmin(a_l - c_l, a_r - c_r), 0.0); 
  const double s_r = fmax(fmax(a_l + c_l, a_r + c_r), 0.0); 
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
  const double *b_l, const double *b_r, const double *u_l, const double *u_r,
  const double *theta_l, const double *theta_r, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT flux, double* GKYL_RESTRICT mass_flux) 
{ 
  // w[NDIM]:   Cell-center coordinates (current cell). 
  // dxv[NDIM]: Cell spacing (current cell). 
  // b_l/b_r:   Magnetic field unit vector [bx, by, bz] in the left/current cell. 
  // u_l/u_r:   Flow velocity [ux, uy, uz] in the left/current cell. 
  // theta_l/theta_r: p_perp/rho in the left/current cell (acoustic floor of the HLL 
  //            speeds where b is perpendicular to the face); NULL disables the floor. 
  // f_l/f_r:   Distribution functions [F_0, G] in the left/current cell. 
  // flux:      Output modal surface fluxes on the lower face of the current cell, 
  //            [F_0 flux, G flux] in direction y. 
  // mass_flux: Output velocity integral over this cell of the F_0 flux, on the 
  //            configuration face basis (this cell's share of the kinetic mass flux). 
  // Returns the stable time-step estimate from the largest speed bound over the nodes. 
  const double dx1 = 2.0/dxv[1]; 
  const double wv = w[2], dv = dxv[2]; 
  const double *bd_l = &b_l[4]; 
  const double *bd_r = &b_r[4]; 
  const double *ud_l = &u_l[4]; 
  const double *ud_r = &u_r[4]; 
 
  double bl_q[3] = {0.0}; 
  double ul_q[3] = {0.0}; 
  double br_q[3] = {0.0}; 
  double ur_q[3] = {0.0}; 
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
    a_l[i*3 + 0] = bl_q[i]*(wv + -0.38729833462074165*dv) + ul_q[i]; 
    a_r[i*3 + 0] = br_q[i]*(wv + -0.38729833462074165*dv) + ur_q[i]; 
    a_l[i*3 + 1] = bl_q[i]*(wv + 0.0*dv) + ul_q[i]; 
    a_r[i*3 + 1] = br_q[i]*(wv + 0.0*dv) + ur_q[i]; 
    a_l[i*3 + 2] = bl_q[i]*(wv + 0.38729833462074165*dv) + ul_q[i]; 
    a_r[i*3 + 2] = br_q[i]*(wv + 0.38729833462074165*dv) + ur_q[i]; 
  } 
 
  double cl_q[3] = {0.0}; 
  double cr_q[3] = {0.0}; 
  if (theta_l && theta_r) { 
    double tl_q[3] = {0.0}; 
    double tr_q[3] = {0.0}; 
    tl_q[0] = -(1.1618950038622249*theta_l[3])+0.8660254037844386*theta_l[2]-0.6708203932499369*theta_l[1]+0.5*theta_l[0]; 
    tl_q[1] = 0.8660254037844386*theta_l[2]+0.5*theta_l[0]; 
    tl_q[2] = 1.1618950038622249*theta_l[3]+0.8660254037844386*theta_l[2]+0.6708203932499369*theta_l[1]+0.5*theta_l[0]; 
    tr_q[0] = 1.1618950038622249*theta_r[3]-0.8660254037844386*theta_r[2]-0.6708203932499369*theta_r[1]+0.5*theta_r[0]; 
    tr_q[1] = 0.5*theta_r[0]-0.8660254037844386*theta_r[2]; 
    tr_q[2] = -(1.1618950038622249*theta_r[3])-0.8660254037844386*theta_r[2]+0.6708203932499369*theta_r[1]+0.5*theta_r[0]; 
    for (int i = 0; i < 3; ++i) { 
      cl_q[i] = sqrt(3.0*fmax(tl_q[i], 0.0)*fmax(1.0 - bl_q[i]*bl_q[i], 0.0)); 
      cr_q[i] = sqrt(3.0*fmax(tr_q[i], 0.0)*fmax(1.0 - br_q[i]*br_q[i], 0.0)); 
    } 
  } 
 
  double G_l[6]; 
  double G_r[6]; 
  double Fhat_nodal[9]; 
  double alpha_max = 0.0; 
  for (int c = 0; c < 2; ++c) { 
    for (int tid = 0; tid < 6; ++tid) pkpm_conf_flux_y_2x1v_ser_p1_g(tid, &f_l[c*12], &f_r[c*12], G_l, G_r); 
    for (int i = 0; i < 3; ++i) { 
      for (int j = 0; j < 3; ++j) { 
        alpha_max = fmax(alpha_max, pkpm_conf_flux_y_2x1v_ser_p1_node(i, j, a_l[i*3 + j], a_r[i*3 + j], cl_q[i], cr_q[i], G_l, G_r, Fhat_nodal)); 
      } 
    } 
    for (int tid = 0; tid < 9; ++tid) pkpm_conf_flux_y_2x1v_ser_p1_prj(tid, Fhat_nodal, &flux[c*6]); 
    if (c == 0) { 
      // Velocity integral of the F_0 flux over this cell at the face nodes, on the face basis. 
      double mq[3] = {0.0}; 
      for (int i = 0; i < 3; ++i) { 
        mq[i] = 0.5*dv*(0.5555555555555556*Fhat_nodal[i*3 + 0] + 0.8888888888888888*Fhat_nodal[i*3 + 1] + 0.5555555555555556*Fhat_nodal[i*3 + 2]); 
      } 
      mass_flux[0] = 0.39283710065919303*mq[2]+0.6285393610547091*mq[1]+0.39283710065919303*mq[0]; 
      mass_flux[1] = 0.5270462766947298*mq[2]-0.5270462766947298*mq[0]; 
    } 
  } 
  return 1.5*dx1*alpha_max; 
} 

GKYL_CU_DH void pkpm_conf_flux_corr_y_2x1v_ser_p1(double mass, const double *moms_l, const double *moms_r,
  const double *roe_mass_flux, const double *kin_mass_flux, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT flux) 
{ 
  // mass:          Species mass (moms and roe_mass_flux are mass densities and fluxes). 
  // moms_l/r:      [rho, p_parallel, p_perp, M1] in the left/current cell. 
  // roe_mass_flux: Mass flux of the momentum solve through the lower face, face modal. 
  // kin_mass_flux: Velocity integral of the stored F_0 flux through the face, face modal. 
  // f_l/f_r:       Distribution functions [F_0, G] in the left/current cell. 
  // flux:          Stored modal fluxes [F_0 flux, G flux] of the lower face in direction y, 
  //                incremented so that the F_0 flux integrates to roe_mass_flux/mass. 
  const double *rho_l = &moms_l[0]; 
  const double *rho_r = &moms_r[0]; 
  double rl_q[3] = {0.0}; 
  double rr_q[3] = {0.0}; 
  rl_q[0] = -(1.1618950038622249*rho_l[3])+0.8660254037844386*rho_l[2]-0.6708203932499369*rho_l[1]+0.5*rho_l[0]; 
  rl_q[1] = 0.8660254037844386*rho_l[2]+0.5*rho_l[0]; 
  rl_q[2] = 1.1618950038622249*rho_l[3]+0.8660254037844386*rho_l[2]+0.6708203932499369*rho_l[1]+0.5*rho_l[0]; 
  rr_q[0] = 1.1618950038622249*rho_r[3]-0.8660254037844386*rho_r[2]-0.6708203932499369*rho_r[1]+0.5*rho_r[0]; 
  rr_q[1] = 0.5*rho_r[0]-0.8660254037844386*rho_r[2]; 
  rr_q[2] = -(1.1618950038622249*rho_r[3])-0.8660254037844386*rho_r[2]+0.6708203932499369*rho_r[1]+0.5*rho_r[0]; 
  // Face-modal defect (number flux) and its value at the face nodes. 
  double dm[2] = {0.0}; 
  for (int k = 0; k < 2; ++k) dm[k] = roe_mass_flux[k]/mass - kin_mass_flux[k]; 
  double dq[3] = {0.0}; 
  dq[0] = 0.7071067811865475*dm[0]-0.9486832980505137*dm[1]; 
  dq[1] = 0.7071067811865475*dm[0]; 
  dq[2] = 0.9486832980505137*dm[1]+0.7071067811865475*dm[0]; 
  // Defect per unit face-averaged number density, applied only where the correction adds 
  // dissipation to the density (dq (n_r - n_l) <= 0) and the density is positive. 
  double sh[3] = {0.0}; 
  for (int i = 0; i < 3; ++i) { 
    const double nbar = 0.5*(rl_q[i] + rr_q[i])/mass; 
    sh[i] = (nbar > 0.0 && dq[i]*(rr_q[i] - rl_q[i]) <= 0.0) ? dq[i]/nbar : 0.0; 
  } 
 
  double G_l[6]; 
  double G_r[6]; 
  double corr_nodal[9]; 
  double tmp[6]; 
  for (int c = 0; c < 2; ++c) { 
    for (int tid = 0; tid < 6; ++tid) pkpm_conf_flux_y_2x1v_ser_p1_g(tid, &f_l[c*12], &f_r[c*12], G_l, G_r); 
    for (int i = 0; i < 3; ++i) { 
      for (int j = 0; j < 3; ++j) { 
        double fl = 0.0, fr = 0.0; 
        for (int a = 0; a < 2; ++a) { 
          fl += vst_pkpm_2x1v_p1_ph_x1_Cm[i*2 + a]*G_l[j*2 + a]; 
          fr += vst_pkpm_2x1v_p1_ph_x1_Cm[i*2 + a]*G_r[j*2 + a]; 
        } 
        corr_nodal[i*3 + j] = 0.5*sh[i]*(fl + fr); 
      } 
    } 
    for (int k = 0; k < 6; ++k) tmp[k] = 0.0; 
    for (int tid = 0; tid < 9; ++tid) pkpm_conf_flux_y_2x1v_ser_p1_prj(tid, corr_nodal, tmp); 
    for (int k = 0; k < 6; ++k) flux[c*6 + k] += tmp[k]; 
  } 
} 
