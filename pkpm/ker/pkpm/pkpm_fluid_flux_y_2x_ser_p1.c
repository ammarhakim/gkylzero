#include <gkyl_pkpm_fluid_kernels.h> 
GKYL_CU_DH double pkpm_fluid_flux_y_2x_ser_p1(
  const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom,
  const double *w, const double *dxv,
  const double *moms_l, const double *moms_r, const double *bvar_l, const double *bvar_r,
  const double *u_l, const double *u_r, const double *euler_l, const double *euler_r,
  double* GKYL_RESTRICT flux, double* GKYL_RESTRICT mass_flux) 
{ 
  // wv_eqn:      Ten-moment wave equation for the Roe fluctuations. 
  // geom:        Wave geometry of the current cell (normal and tangents of its lower face). 
  // w/dxv:       Cell center and spacing of the current cell. 
  // moms_l/r:    [rho, p_parallel, p_perp, M1] in the left/current cell. 
  // bvar_l/r:    Magnetic field unit vector [bx, by, bz] in the left/current cell. 
  // u_l/r:       Flow velocity [ux, uy, uz] in the left/current cell. 
  // euler_l/r:   Momentum [rho ux, rho uy, rho uz] in the left/current cell. 
  // flux:        Output momentum flux [F_x, F_y, F_z] on the lower face, surface modal. 
  // mass_flux:   Output mass flux of the same Riemann solve on the lower face, surface modal. 
  // Returns the stable time-step estimate from the normal speed bound at the nodes. 
  const double dx1 = 2.0/dxv[1]; 
  const double *rho_l = &moms_l[0]; 
  const double *ppar_l = &moms_l[4]; 
  const double *pperp_l = &moms_l[8]; 
  const double *rho_r = &moms_r[0]; 
  const double *ppar_r = &moms_r[4]; 
  const double *pperp_r = &moms_r[8]; 
  const double *bx_l = &bvar_l[0]; 
  const double *bx_r = &bvar_r[0]; 
  const double *ux_l = &u_l[0]; 
  const double *ux_r = &u_r[0]; 
  const double *Mx_l = &euler_l[0]; 
  const double *Mx_r = &euler_r[0]; 
  const double *by_l = &bvar_l[4]; 
  const double *by_r = &bvar_r[4]; 
  const double *uy_l = &u_l[4]; 
  const double *uy_r = &u_r[4]; 
  const double *My_l = &euler_l[4]; 
  const double *My_r = &euler_r[4]; 
  const double *bz_l = &bvar_l[8]; 
  const double *bz_r = &bvar_r[8]; 
  const double *uz_l = &u_l[8]; 
  const double *uz_r = &u_r[8]; 
  const double *Mz_l = &euler_l[8]; 
  const double *Mz_r = &euler_r[8]; 
 
  double rho_lq[3] = {0.0}; 
  double rho_rq[3] = {0.0}; 
  double ppar_lq[3] = {0.0}; 
  double ppar_rq[3] = {0.0}; 
  double pperp_lq[3] = {0.0}; 
  double pperp_rq[3] = {0.0}; 
  double bx_lq[3] = {0.0}; 
  double bx_rq[3] = {0.0}; 
  double by_lq[3] = {0.0}; 
  double by_rq[3] = {0.0}; 
  double bz_lq[3] = {0.0}; 
  double bz_rq[3] = {0.0}; 
  double ux_lq[3] = {0.0}; 
  double ux_rq[3] = {0.0}; 
  double uy_lq[3] = {0.0}; 
  double uy_rq[3] = {0.0}; 
  double uz_lq[3] = {0.0}; 
  double uz_rq[3] = {0.0}; 
  double Mx_lq[3] = {0.0}; 
  double Mx_rq[3] = {0.0}; 
  double My_lq[3] = {0.0}; 
  double My_rq[3] = {0.0}; 
  double Mz_lq[3] = {0.0}; 
  double Mz_rq[3] = {0.0}; 
  rho_lq[0] = -(1.1618950038622249*rho_l[3])+0.8660254037844386*rho_l[2]-0.6708203932499369*rho_l[1]+0.5*rho_l[0]; 
  rho_lq[1] = 0.8660254037844386*rho_l[2]+0.5*rho_l[0]; 
  rho_lq[2] = 1.1618950038622249*rho_l[3]+0.8660254037844386*rho_l[2]+0.6708203932499369*rho_l[1]+0.5*rho_l[0]; 
  rho_rq[0] = 1.1618950038622249*rho_r[3]-0.8660254037844386*rho_r[2]-0.6708203932499369*rho_r[1]+0.5*rho_r[0]; 
  rho_rq[1] = 0.5*rho_r[0]-0.8660254037844386*rho_r[2]; 
  rho_rq[2] = -(1.1618950038622249*rho_r[3])-0.8660254037844386*rho_r[2]+0.6708203932499369*rho_r[1]+0.5*rho_r[0]; 
  ppar_lq[0] = -(1.1618950038622249*ppar_l[3])+0.8660254037844386*ppar_l[2]-0.6708203932499369*ppar_l[1]+0.5*ppar_l[0]; 
  ppar_lq[1] = 0.8660254037844386*ppar_l[2]+0.5*ppar_l[0]; 
  ppar_lq[2] = 1.1618950038622249*ppar_l[3]+0.8660254037844386*ppar_l[2]+0.6708203932499369*ppar_l[1]+0.5*ppar_l[0]; 
  ppar_rq[0] = 1.1618950038622249*ppar_r[3]-0.8660254037844386*ppar_r[2]-0.6708203932499369*ppar_r[1]+0.5*ppar_r[0]; 
  ppar_rq[1] = 0.5*ppar_r[0]-0.8660254037844386*ppar_r[2]; 
  ppar_rq[2] = -(1.1618950038622249*ppar_r[3])-0.8660254037844386*ppar_r[2]+0.6708203932499369*ppar_r[1]+0.5*ppar_r[0]; 
  pperp_lq[0] = -(1.1618950038622249*pperp_l[3])+0.8660254037844386*pperp_l[2]-0.6708203932499369*pperp_l[1]+0.5*pperp_l[0]; 
  pperp_lq[1] = 0.8660254037844386*pperp_l[2]+0.5*pperp_l[0]; 
  pperp_lq[2] = 1.1618950038622249*pperp_l[3]+0.8660254037844386*pperp_l[2]+0.6708203932499369*pperp_l[1]+0.5*pperp_l[0]; 
  pperp_rq[0] = 1.1618950038622249*pperp_r[3]-0.8660254037844386*pperp_r[2]-0.6708203932499369*pperp_r[1]+0.5*pperp_r[0]; 
  pperp_rq[1] = 0.5*pperp_r[0]-0.8660254037844386*pperp_r[2]; 
  pperp_rq[2] = -(1.1618950038622249*pperp_r[3])-0.8660254037844386*pperp_r[2]+0.6708203932499369*pperp_r[1]+0.5*pperp_r[0]; 
  bx_lq[0] = -(1.1618950038622249*bx_l[3])+0.8660254037844386*bx_l[2]-0.6708203932499369*bx_l[1]+0.5*bx_l[0]; 
  bx_lq[1] = 0.8660254037844386*bx_l[2]+0.5*bx_l[0]; 
  bx_lq[2] = 1.1618950038622249*bx_l[3]+0.8660254037844386*bx_l[2]+0.6708203932499369*bx_l[1]+0.5*bx_l[0]; 
  bx_rq[0] = 1.1618950038622249*bx_r[3]-0.8660254037844386*bx_r[2]-0.6708203932499369*bx_r[1]+0.5*bx_r[0]; 
  bx_rq[1] = 0.5*bx_r[0]-0.8660254037844386*bx_r[2]; 
  bx_rq[2] = -(1.1618950038622249*bx_r[3])-0.8660254037844386*bx_r[2]+0.6708203932499369*bx_r[1]+0.5*bx_r[0]; 
  by_lq[0] = -(1.1618950038622249*by_l[3])+0.8660254037844386*by_l[2]-0.6708203932499369*by_l[1]+0.5*by_l[0]; 
  by_lq[1] = 0.8660254037844386*by_l[2]+0.5*by_l[0]; 
  by_lq[2] = 1.1618950038622249*by_l[3]+0.8660254037844386*by_l[2]+0.6708203932499369*by_l[1]+0.5*by_l[0]; 
  by_rq[0] = 1.1618950038622249*by_r[3]-0.8660254037844386*by_r[2]-0.6708203932499369*by_r[1]+0.5*by_r[0]; 
  by_rq[1] = 0.5*by_r[0]-0.8660254037844386*by_r[2]; 
  by_rq[2] = -(1.1618950038622249*by_r[3])-0.8660254037844386*by_r[2]+0.6708203932499369*by_r[1]+0.5*by_r[0]; 
  bz_lq[0] = -(1.1618950038622249*bz_l[3])+0.8660254037844386*bz_l[2]-0.6708203932499369*bz_l[1]+0.5*bz_l[0]; 
  bz_lq[1] = 0.8660254037844386*bz_l[2]+0.5*bz_l[0]; 
  bz_lq[2] = 1.1618950038622249*bz_l[3]+0.8660254037844386*bz_l[2]+0.6708203932499369*bz_l[1]+0.5*bz_l[0]; 
  bz_rq[0] = 1.1618950038622249*bz_r[3]-0.8660254037844386*bz_r[2]-0.6708203932499369*bz_r[1]+0.5*bz_r[0]; 
  bz_rq[1] = 0.5*bz_r[0]-0.8660254037844386*bz_r[2]; 
  bz_rq[2] = -(1.1618950038622249*bz_r[3])-0.8660254037844386*bz_r[2]+0.6708203932499369*bz_r[1]+0.5*bz_r[0]; 
  ux_lq[0] = -(1.1618950038622249*ux_l[3])+0.8660254037844386*ux_l[2]-0.6708203932499369*ux_l[1]+0.5*ux_l[0]; 
  ux_lq[1] = 0.8660254037844386*ux_l[2]+0.5*ux_l[0]; 
  ux_lq[2] = 1.1618950038622249*ux_l[3]+0.8660254037844386*ux_l[2]+0.6708203932499369*ux_l[1]+0.5*ux_l[0]; 
  ux_rq[0] = 1.1618950038622249*ux_r[3]-0.8660254037844386*ux_r[2]-0.6708203932499369*ux_r[1]+0.5*ux_r[0]; 
  ux_rq[1] = 0.5*ux_r[0]-0.8660254037844386*ux_r[2]; 
  ux_rq[2] = -(1.1618950038622249*ux_r[3])-0.8660254037844386*ux_r[2]+0.6708203932499369*ux_r[1]+0.5*ux_r[0]; 
  uy_lq[0] = -(1.1618950038622249*uy_l[3])+0.8660254037844386*uy_l[2]-0.6708203932499369*uy_l[1]+0.5*uy_l[0]; 
  uy_lq[1] = 0.8660254037844386*uy_l[2]+0.5*uy_l[0]; 
  uy_lq[2] = 1.1618950038622249*uy_l[3]+0.8660254037844386*uy_l[2]+0.6708203932499369*uy_l[1]+0.5*uy_l[0]; 
  uy_rq[0] = 1.1618950038622249*uy_r[3]-0.8660254037844386*uy_r[2]-0.6708203932499369*uy_r[1]+0.5*uy_r[0]; 
  uy_rq[1] = 0.5*uy_r[0]-0.8660254037844386*uy_r[2]; 
  uy_rq[2] = -(1.1618950038622249*uy_r[3])-0.8660254037844386*uy_r[2]+0.6708203932499369*uy_r[1]+0.5*uy_r[0]; 
  uz_lq[0] = -(1.1618950038622249*uz_l[3])+0.8660254037844386*uz_l[2]-0.6708203932499369*uz_l[1]+0.5*uz_l[0]; 
  uz_lq[1] = 0.8660254037844386*uz_l[2]+0.5*uz_l[0]; 
  uz_lq[2] = 1.1618950038622249*uz_l[3]+0.8660254037844386*uz_l[2]+0.6708203932499369*uz_l[1]+0.5*uz_l[0]; 
  uz_rq[0] = 1.1618950038622249*uz_r[3]-0.8660254037844386*uz_r[2]-0.6708203932499369*uz_r[1]+0.5*uz_r[0]; 
  uz_rq[1] = 0.5*uz_r[0]-0.8660254037844386*uz_r[2]; 
  uz_rq[2] = -(1.1618950038622249*uz_r[3])-0.8660254037844386*uz_r[2]+0.6708203932499369*uz_r[1]+0.5*uz_r[0]; 
  Mx_lq[0] = -(1.1618950038622249*Mx_l[3])+0.8660254037844386*Mx_l[2]-0.6708203932499369*Mx_l[1]+0.5*Mx_l[0]; 
  Mx_lq[1] = 0.8660254037844386*Mx_l[2]+0.5*Mx_l[0]; 
  Mx_lq[2] = 1.1618950038622249*Mx_l[3]+0.8660254037844386*Mx_l[2]+0.6708203932499369*Mx_l[1]+0.5*Mx_l[0]; 
  Mx_rq[0] = 1.1618950038622249*Mx_r[3]-0.8660254037844386*Mx_r[2]-0.6708203932499369*Mx_r[1]+0.5*Mx_r[0]; 
  Mx_rq[1] = 0.5*Mx_r[0]-0.8660254037844386*Mx_r[2]; 
  Mx_rq[2] = -(1.1618950038622249*Mx_r[3])-0.8660254037844386*Mx_r[2]+0.6708203932499369*Mx_r[1]+0.5*Mx_r[0]; 
  My_lq[0] = -(1.1618950038622249*My_l[3])+0.8660254037844386*My_l[2]-0.6708203932499369*My_l[1]+0.5*My_l[0]; 
  My_lq[1] = 0.8660254037844386*My_l[2]+0.5*My_l[0]; 
  My_lq[2] = 1.1618950038622249*My_l[3]+0.8660254037844386*My_l[2]+0.6708203932499369*My_l[1]+0.5*My_l[0]; 
  My_rq[0] = 1.1618950038622249*My_r[3]-0.8660254037844386*My_r[2]-0.6708203932499369*My_r[1]+0.5*My_r[0]; 
  My_rq[1] = 0.5*My_r[0]-0.8660254037844386*My_r[2]; 
  My_rq[2] = -(1.1618950038622249*My_r[3])-0.8660254037844386*My_r[2]+0.6708203932499369*My_r[1]+0.5*My_r[0]; 
  Mz_lq[0] = -(1.1618950038622249*Mz_l[3])+0.8660254037844386*Mz_l[2]-0.6708203932499369*Mz_l[1]+0.5*Mz_l[0]; 
  Mz_lq[1] = 0.8660254037844386*Mz_l[2]+0.5*Mz_l[0]; 
  Mz_lq[2] = 1.1618950038622249*Mz_l[3]+0.8660254037844386*Mz_l[2]+0.6708203932499369*Mz_l[1]+0.5*Mz_l[0]; 
  Mz_rq[0] = 1.1618950038622249*Mz_r[3]-0.8660254037844386*Mz_r[2]-0.6708203932499369*Mz_r[1]+0.5*Mz_r[0]; 
  Mz_rq[1] = 0.5*Mz_r[0]-0.8660254037844386*Mz_r[2]; 
  Mz_rq[2] = -(1.1618950038622249*Mz_r[3])-0.8660254037844386*Mz_r[2]+0.6708203932499369*Mz_r[1]+0.5*Mz_r[0]; 
  double F_q[3][3] = {{0.0}}; 
  double F0_q[3] = {0.0}; 
  double lambda_max = 0.0; 
  double q_l[10] = {0.0}, q_r[10] = {0.0}, q_l_local[10] = {0.0}, q_r_local[10] = {0.0}; 
  double delta[10] = {0.0}, waves[50] = {0.0}, speeds[5] = {0.0}; 
  double amdq_local[10] = {0.0}, apdq_local[10] = {0.0}, amdq[10] = {0.0}, apdq[10] = {0.0}; 
  for (int n = 0; n < 3; ++n) { 
    // Full pressure tensor from the primitive factors at the node. 
    const double b_l[3] = {bx_lq[n], by_lq[n], bz_lq[n]}, b_r[3] = {bx_rq[n], by_rq[n], bz_rq[n]}; 
    const double u_lv[3] = {ux_lq[n], uy_lq[n], uz_lq[n]}, u_rv[3] = {ux_rq[n], uy_rq[n], uz_rq[n]}; 
    const double M_lv[3] = {Mx_lq[n], My_lq[n], Mz_lq[n]}, M_rv[3] = {Mx_rq[n], My_rq[n], Mz_rq[n]}; 
    const double dp_l = ppar_lq[n] - pperp_lq[n], dp_r = ppar_rq[n] - pperp_rq[n]; 
    double P_l[3][3], P_r[3][3]; 
    for (int i = 0; i < 3; ++i) { 
      for (int j = 0; j < 3; ++j) { 
        P_l[i][j] = dp_l*b_l[i]*b_l[j] + (i == j ? pperp_lq[n] : 0.0); 
        P_r[i][j] = dp_r*b_r[i]*b_r[j] + (i == j ? pperp_rq[n] : 0.0); 
      } 
    } 
    const int d = 1; 
    // Admissibility of both traces, decided before any division or square root. 
    const bool ok_l = rho_lq[n] > 0.0 && P_l[d][d] > 0.0; 
    const bool ok_r = rho_rq[n] > 0.0 && P_r[d][d] > 0.0; 
    const double c_l = ok_l ? sqrt(3.0*P_l[d][d]/rho_lq[n]) : 0.0; 
    const double c_r = ok_r ? sqrt(3.0*P_r[d][d]/rho_rq[n]) : 0.0; 
    lambda_max = fmax(lambda_max, fmax(fabs(u_lv[d]) + c_l, fabs(u_rv[d]) + c_r)); 
    double F_lv[3], F_rv[3]; 
    for (int i = 0; i < 3; ++i) { 
      F_lv[i] = M_lv[i]*u_lv[d] + P_l[i][d]; 
      F_rv[i] = M_rv[i]*u_rv[d] + P_r[i][d]; 
    } 
    if (ok_l && ok_r) { 
      // Ten-moment Roe flux of the nodal states. 
      q_l[0] = rho_lq[n]; q_r[0] = rho_rq[n]; 
      for (int i = 0; i < 3; ++i) { q_l[1+i] = M_lv[i]; q_r[1+i] = M_rv[i]; } 
      q_l[4] = P_l[0][0] + M_lv[0]*M_lv[0]/rho_lq[n]; q_r[4] = P_r[0][0] + M_rv[0]*M_rv[0]/rho_rq[n]; 
      q_l[5] = P_l[0][1] + M_lv[0]*M_lv[1]/rho_lq[n]; q_r[5] = P_r[0][1] + M_rv[0]*M_rv[1]/rho_rq[n]; 
      q_l[6] = P_l[0][2] + M_lv[0]*M_lv[2]/rho_lq[n]; q_r[6] = P_r[0][2] + M_rv[0]*M_rv[2]/rho_rq[n]; 
      q_l[7] = P_l[1][1] + M_lv[1]*M_lv[1]/rho_lq[n]; q_r[7] = P_r[1][1] + M_rv[1]*M_rv[1]/rho_rq[n]; 
      q_l[8] = P_l[1][2] + M_lv[1]*M_lv[2]/rho_lq[n]; q_r[8] = P_r[1][2] + M_rv[1]*M_rv[2]/rho_rq[n]; 
      q_l[9] = P_l[2][2] + M_lv[2]*M_lv[2]/rho_lq[n]; q_r[9] = P_r[2][2] + M_rv[2]*M_rv[2]/rho_rq[n]; 
      rot_to_local(wv_eqn, geom->tau1[d], geom->tau2[d], geom->norm[d], q_l, q_l_local); 
      rot_to_local(wv_eqn, geom->tau1[d], geom->tau2[d], geom->norm[d], q_r, q_r_local); 
      for (int i = 0; i < 10; ++i) delta[i] = q_r_local[i] - q_l_local[i]; 
      wave_roe(wv_eqn, delta, q_l_local, q_r_local, waves, speeds); 
      for (int i = 0; i < 5; ++i) speeds[i] *= geom->lenr[d]; 
      qfluct_roe(wv_eqn, q_l_local, q_r_local, waves, speeds, amdq_local, apdq_local); 
      rot_to_global(wv_eqn, geom->tau1[d], geom->tau2[d], geom->norm[d], amdq_local, amdq); 
      rot_to_global(wv_eqn, geom->tau1[d], geom->tau2[d], geom->norm[d], apdq_local, apdq); 
      for (int i = 0; i < 3; ++i) F_q[i][n] = 0.5*(F_lv[i] + F_rv[i]) - 0.5*(apdq[1+i] - amdq[1+i]); 
      F0_q[n] = 0.5*(M_lv[d] + M_rv[d]) - 0.5*(apdq[0] - amdq[0]); 
    } else { 
      // HLL flux of the momentum with the normal speed bounds (c = 0 on an inadmissible trace). 
      const double S_L = fmin(u_lv[d] - c_l, u_rv[d] - c_r); 
      const double S_R = fmax(u_lv[d] + c_l, u_rv[d] + c_r); 
      for (int i = 0; i < 3; ++i) { 
        if (S_L >= 0.0) F_q[i][n] = F_lv[i]; 
        else if (S_R <= 0.0) F_q[i][n] = F_rv[i]; 
        else F_q[i][n] = (S_R*F_lv[i] - S_L*F_rv[i] + S_L*S_R*(M_rv[i] - M_lv[i]))/(S_R - S_L); 
      } 
      if (S_L >= 0.0) F0_q[n] = M_lv[d]; 
      else if (S_R <= 0.0) F0_q[n] = M_rv[d]; 
      else F0_q[n] = (S_R*M_lv[d] - S_L*M_rv[d] + S_L*S_R*(rho_rq[n] - rho_lq[n]))/(S_R - S_L); 
    } 
  } 
 
  // Projection of the nodal fluxes onto the surface basis. 
  flux[0] = 0.39283710065919303*F_q[0][2]+0.6285393610547091*F_q[0][1]+0.39283710065919303*F_q[0][0]; 
  flux[1] = 0.5270462766947298*F_q[0][2]-0.5270462766947298*F_q[0][0]; 
  flux[2] = 0.39283710065919303*F_q[1][2]+0.6285393610547091*F_q[1][1]+0.39283710065919303*F_q[1][0]; 
  flux[3] = 0.5270462766947298*F_q[1][2]-0.5270462766947298*F_q[1][0]; 
  flux[4] = 0.39283710065919303*F_q[2][2]+0.6285393610547091*F_q[2][1]+0.39283710065919303*F_q[2][0]; 
  flux[5] = 0.5270462766947298*F_q[2][2]-0.5270462766947298*F_q[2][0]; 
  mass_flux[0] = 0.39283710065919303*F0_q[2]+0.6285393610547091*F0_q[1]+0.39283710065919303*F0_q[0]; 
  mass_flux[1] = 0.5270462766947298*F0_q[2]-0.5270462766947298*F0_q[0]; 
  return 1.5*dx1*lambda_max; 
} 
