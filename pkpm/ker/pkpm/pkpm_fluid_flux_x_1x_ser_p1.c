#include <gkyl_pkpm_fluid_kernels.h> 
GKYL_CU_DH double pkpm_fluid_flux_x_1x_ser_p1(
  const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom,
  const double *w, const double *dxv,
  const double *moms_l, const double *moms_r, const double *bvar_l, const double *bvar_r,
  const double *u_l, const double *u_r, const double *fluid_l, const double *fluid_r,
  double* GKYL_RESTRICT flux) 
{ 
  // wv_eqn:      Ten-moment wave equation (unused; kept for the limiter interface). 
  // geom:        Wave geometry of the current cell (unused). 
  // w/dxv:       Cell center and spacing of the current cell. 
  // moms_l/r:    Kinetic coupling moments [p_par, p_perp, q_par, q_perp] in the left/current cell. 
  // bvar_l/r:    Magnetic field unit vector [bx, by, bz] in the left/current cell. 
  // u_l/r:       Flow velocity [ux, uy, uz] in the left/current cell. 
  // fluid_l/r:   Evolved fluid state [rho, rho ux, rho uy, rho uz, E] in the left/current cell. 
  // flux:        Output fluxes [F_rho, F_x, F_y, F_z, F_E] on the lower face, surface modal. 
  // Returns the stable time-step estimate from the normal speed bound at the nodes. 
  const double dx1 = 2.0/dxv[0]; 
  const double *rho_l = &fluid_l[0]; 
  const double *E_l = &fluid_l[8]; 
  const double *ppar_l = &moms_l[0]; 
  const double *pperp_l = &moms_l[2]; 
  const double *qpar_l = &moms_l[4]; 
  const double *qperp_l = &moms_l[6]; 
  const double *rho_r = &fluid_r[0]; 
  const double *E_r = &fluid_r[8]; 
  const double *ppar_r = &moms_r[0]; 
  const double *pperp_r = &moms_r[2]; 
  const double *qpar_r = &moms_r[4]; 
  const double *qperp_r = &moms_r[6]; 
  const double *bx_l = &bvar_l[0]; 
  const double *bx_r = &bvar_r[0]; 
  const double *ux_l = &u_l[0]; 
  const double *ux_r = &u_r[0]; 
  const double *Mx_l = &fluid_l[2]; 
  const double *Mx_r = &fluid_r[2]; 
  const double *by_l = &bvar_l[2]; 
  const double *by_r = &bvar_r[2]; 
  const double *uy_l = &u_l[2]; 
  const double *uy_r = &u_r[2]; 
  const double *My_l = &fluid_l[4]; 
  const double *My_r = &fluid_r[4]; 
  const double *bz_l = &bvar_l[4]; 
  const double *bz_r = &bvar_r[4]; 
  const double *uz_l = &u_l[4]; 
  const double *uz_r = &u_r[4]; 
  const double *Mz_l = &fluid_l[6]; 
  const double *Mz_r = &fluid_r[6]; 
 
  double rho_lq[1] = {0.0}; 
  double rho_rq[1] = {0.0}; 
  double E_lq[1] = {0.0}; 
  double E_rq[1] = {0.0}; 
  double ppar_lq[1] = {0.0}; 
  double ppar_rq[1] = {0.0}; 
  double pperp_lq[1] = {0.0}; 
  double pperp_rq[1] = {0.0}; 
  double qpar_lq[1] = {0.0}; 
  double qpar_rq[1] = {0.0}; 
  double qperp_lq[1] = {0.0}; 
  double qperp_rq[1] = {0.0}; 
  double bx_lq[1] = {0.0}; 
  double bx_rq[1] = {0.0}; 
  double by_lq[1] = {0.0}; 
  double by_rq[1] = {0.0}; 
  double bz_lq[1] = {0.0}; 
  double bz_rq[1] = {0.0}; 
  double ux_lq[1] = {0.0}; 
  double ux_rq[1] = {0.0}; 
  double uy_lq[1] = {0.0}; 
  double uy_rq[1] = {0.0}; 
  double uz_lq[1] = {0.0}; 
  double uz_rq[1] = {0.0}; 
  double Mx_lq[1] = {0.0}; 
  double Mx_rq[1] = {0.0}; 
  double My_lq[1] = {0.0}; 
  double My_rq[1] = {0.0}; 
  double Mz_lq[1] = {0.0}; 
  double Mz_rq[1] = {0.0}; 
  rho_lq[0] = 1.224744871391589*rho_l[1]+0.7071067811865475*rho_l[0]; 
  rho_rq[0] = 0.7071067811865475*rho_r[0]-1.224744871391589*rho_r[1]; 
  E_lq[0] = 1.224744871391589*E_l[1]+0.7071067811865475*E_l[0]; 
  E_rq[0] = 0.7071067811865475*E_r[0]-1.224744871391589*E_r[1]; 
  ppar_lq[0] = 1.224744871391589*ppar_l[1]+0.7071067811865475*ppar_l[0]; 
  ppar_rq[0] = 0.7071067811865475*ppar_r[0]-1.224744871391589*ppar_r[1]; 
  pperp_lq[0] = 1.224744871391589*pperp_l[1]+0.7071067811865475*pperp_l[0]; 
  pperp_rq[0] = 0.7071067811865475*pperp_r[0]-1.224744871391589*pperp_r[1]; 
  qpar_lq[0] = 1.224744871391589*qpar_l[1]+0.7071067811865475*qpar_l[0]; 
  qpar_rq[0] = 0.7071067811865475*qpar_r[0]-1.224744871391589*qpar_r[1]; 
  qperp_lq[0] = 1.224744871391589*qperp_l[1]+0.7071067811865475*qperp_l[0]; 
  qperp_rq[0] = 0.7071067811865475*qperp_r[0]-1.224744871391589*qperp_r[1]; 
  bx_lq[0] = 1.224744871391589*bx_l[1]+0.7071067811865475*bx_l[0]; 
  bx_rq[0] = 0.7071067811865475*bx_r[0]-1.224744871391589*bx_r[1]; 
  by_lq[0] = 1.224744871391589*by_l[1]+0.7071067811865475*by_l[0]; 
  by_rq[0] = 0.7071067811865475*by_r[0]-1.224744871391589*by_r[1]; 
  bz_lq[0] = 1.224744871391589*bz_l[1]+0.7071067811865475*bz_l[0]; 
  bz_rq[0] = 0.7071067811865475*bz_r[0]-1.224744871391589*bz_r[1]; 
  ux_lq[0] = 1.224744871391589*ux_l[1]+0.7071067811865475*ux_l[0]; 
  ux_rq[0] = 0.7071067811865475*ux_r[0]-1.224744871391589*ux_r[1]; 
  uy_lq[0] = 1.224744871391589*uy_l[1]+0.7071067811865475*uy_l[0]; 
  uy_rq[0] = 0.7071067811865475*uy_r[0]-1.224744871391589*uy_r[1]; 
  uz_lq[0] = 1.224744871391589*uz_l[1]+0.7071067811865475*uz_l[0]; 
  uz_rq[0] = 0.7071067811865475*uz_r[0]-1.224744871391589*uz_r[1]; 
  Mx_lq[0] = 1.224744871391589*Mx_l[1]+0.7071067811865475*Mx_l[0]; 
  Mx_rq[0] = 0.7071067811865475*Mx_r[0]-1.224744871391589*Mx_r[1]; 
  My_lq[0] = 1.224744871391589*My_l[1]+0.7071067811865475*My_l[0]; 
  My_rq[0] = 0.7071067811865475*My_r[0]-1.224744871391589*My_r[1]; 
  Mz_lq[0] = 1.224744871391589*Mz_l[1]+0.7071067811865475*Mz_l[0]; 
  Mz_rq[0] = 0.7071067811865475*Mz_r[0]-1.224744871391589*Mz_r[1]; 
  double F_q[5][1] = {{0.0}}; 
  double lambda_max = 0.0; 
  for (int n = 0; n < 1; ++n) { 
    // Pressure tensor and heat flux from the primitive factors at the node. 
    const double b_l[3] = {bx_lq[n], by_lq[n], bz_lq[n]}, b_r[3] = {bx_rq[n], by_rq[n], bz_rq[n]}; 
    const double u_lv[3] = {ux_lq[n], uy_lq[n], uz_lq[n]}, u_rv[3] = {ux_rq[n], uy_rq[n], uz_rq[n]}; 
    const double M_lv[3] = {Mx_lq[n], My_lq[n], Mz_lq[n]}, M_rv[3] = {Mx_rq[n], My_rq[n], Mz_rq[n]}; 
    const double dp_l = ppar_lq[n] - pperp_lq[n], dp_r = ppar_rq[n] - pperp_rq[n]; 
    const double pperp_l = pperp_lq[n], pperp_r = pperp_rq[n]; 
    double P_l[3][3], P_r[3][3]; 
    for (int i = 0; i < 3; ++i) { 
      for (int j = 0; j < 3; ++j) { 
        P_l[i][j] = dp_l*b_l[i]*b_l[j] + (i == j ? pperp_l : 0.0); 
        P_r[i][j] = dp_r*b_r[i]*b_r[j] + (i == j ? pperp_r : 0.0); 
      } 
    } 
    const int d = 0; 
    const double qn_l = 0.5*b_l[d]*(qpar_lq[n] + 2.0*qperp_lq[n]); 
    const double qn_r = 0.5*b_r[d]*(qpar_rq[n] + 2.0*qperp_rq[n]); 
    // Admissibility of both traces, decided before any division or square root. 
    const bool ok_l = rho_lq[n] > 0.0 && P_l[d][d] > 0.0; 
    const bool ok_r = rho_rq[n] > 0.0 && P_r[d][d] > 0.0; 
    const double c_l = ok_l ? sqrt((3.0 + 2.0*b_l[d]*b_l[d])*P_l[d][d]/rho_lq[n]) : 0.0; 
    const double c_r = ok_r ? sqrt((3.0 + 2.0*b_r[d]*b_r[d])*P_r[d][d]/rho_rq[n]) : 0.0; 
    lambda_max = fmax(lambda_max, fmax(fabs(u_lv[d]) + c_l, fabs(u_rv[d]) + c_r)); 
    double F_lv[5], F_rv[5], q_lv[5], q_rv[5]; 
    double Pu_l = 0.0, Pu_r = 0.0; 
    for (int i = 0; i < 3; ++i) { 
      F_lv[1+i] = M_lv[i]*u_lv[d] + P_l[i][d]; 
      F_rv[1+i] = M_rv[i]*u_rv[d] + P_r[i][d]; 
      q_lv[1+i] = M_lv[i]; q_rv[1+i] = M_rv[i]; 
      Pu_l += P_l[i][d]*u_lv[i]; Pu_r += P_r[i][d]*u_rv[i]; 
    } 
    F_lv[0] = M_lv[d]; F_rv[0] = M_rv[d]; q_lv[0] = rho_lq[n]; q_rv[0] = rho_rq[n]; 
    F_lv[4] = u_lv[d]*E_lq[n] + Pu_l + qn_l; F_rv[4] = u_rv[d]*E_rq[n] + Pu_r + qn_r; q_lv[4] = E_lq[n]; q_rv[4] = E_rq[n]; 
    // HLL flux with the normal speed bounds (c = 0 on an inadmissible trace). 
    const double S_L = fmin(u_lv[d] - c_l, u_rv[d] - c_r); 
    const double S_R = fmax(u_lv[d] + c_l, u_rv[d] + c_r); 
    for (int i = 0; i < 5; ++i) { 
      if (S_L >= 0.0) F_q[i][n] = F_lv[i]; 
      else if (S_R <= 0.0) F_q[i][n] = F_rv[i]; 
      else F_q[i][n] = (S_R*F_lv[i] - S_L*F_rv[i] + S_L*S_R*(q_rv[i] - q_lv[i]))/(S_R - S_L); 
    } 
  } 
 
  // Projection of the nodal fluxes onto the surface basis. 
  flux[0] = F_q[0][0]; 
  flux[1] = F_q[1][0]; 
  flux[2] = F_q[2][0]; 
  flux[3] = F_q[3][0]; 
  flux[4] = F_q[4][0]; 
  return 1.5*dx1*lambda_max; 
} 
