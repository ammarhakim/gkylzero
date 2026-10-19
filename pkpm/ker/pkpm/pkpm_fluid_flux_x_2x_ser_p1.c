#include <gkyl_pkpm_fluid_kernels.h> 
GKYL_CU_DH double pkpm_fluid_flux_x_2x_ser_p1(const double *jacob_pos, const double *bvar_l, const double *bvar_r,
  const double *u_l, const double *u_r, const double *xi_l, const double *xi_r, const double *theta_l,
  const double *theta_r, const double *fluid_l, const double *fluid_r, const double *I, double* GKYL_RESTRICT flux) 
{ 
  // jacob_pos: Position-map Jacobian per configuration direction (cdim blocks of 2 equal entries) of 
  //            either cell (the face Jacobian J/J_x is the same on both sides). 
  // bvar_l/r:  Magnetic field unit vector [bx, by, bz] of the left/right cell. 
  // u_l/r:     Flow velocity [ux, uy, uz] of the left/right cell. 
  // xi_l/r:    Thermal speed [xi, xi^2, p] of the left/right cell (the first expansion is used: 
  //            xi^2 is the square of its trace, the xi of the streaming speed). 
  // theta_l/r: T_perp/T = M0(G) of the left/right cell (p_perp = rho xi^2 theta). 
  // fluid_l/r: PHYSICAL fluid state [rho, rho ux, rho uy, rho uz, E] of the left/right cell (Jacobian divided out). 
  // I:         Nodal moments [I_0, I_1, I_2, J]^l, [I_0, I_1, I_2, J]^r of the kinetic HLL flux at the 
  //            3 face nodes (pkpm_fluid_flux_accum_x_2x1v_ser_p1 summed over the velocity cells); 
  //            they carry the face Jacobian J/J_x. 
  // flux:      Output modal surface fluxes [F_rho, F_x, F_y, F_z, F_E] of the lower face of the right 
  //            cell for the evolved J-weighted fluid. 
  // Returns the largest perpendicular-sound bound |u_n| + xi sqrt(2 theta (1 - b_n^2)) over the 
  //            nodes and traces (the perpendicular pressure is a central flux; its waves need 
  //            their own CFL bound, the parallel streaming has the kinetic one). 
  double cfl_max = 0.0; 
  // face Jacobian J/J_x of the perpendicular-pressure terms (the moments carry it already) 
  const double jperp = jacob_pos[2]; 
  const double *rho_l = &fluid_l[0]; 
  const double *q1_l = &fluid_l[4]; 
  const double *q2_l = &fluid_l[8]; 
  const double *q3_l = &fluid_l[12]; 
  const double *q4_l = &fluid_l[16]; 
  const double *xd_l = &xi_l[0]; 
  const double *thd_l = &theta_l[0]; 
  const double *b0_l = &bvar_l[0]; 
  const double *u0_l = &u_l[0]; 
  const double *b1_l = &bvar_l[4]; 
  const double *u1_l = &u_l[4]; 
  const double *b2_l = &bvar_l[8]; 
  const double *u2_l = &u_l[8]; 
  const double *rho_r = &fluid_r[0]; 
  const double *q1_r = &fluid_r[4]; 
  const double *q2_r = &fluid_r[8]; 
  const double *q3_r = &fluid_r[12]; 
  const double *q4_r = &fluid_r[16]; 
  const double *xd_r = &xi_r[0]; 
  const double *thd_r = &theta_r[0]; 
  const double *b0_r = &bvar_r[0]; 
  const double *u0_r = &u_r[0]; 
  const double *b1_r = &bvar_r[4]; 
  const double *u1_r = &u_r[4]; 
  const double *b2_r = &bvar_r[8]; 
  const double *u2_r = &u_r[8]; 
 
  double rq_l[3] = {0.0}; 
  double xq_l[3] = {0.0}; 
  double thq_l[3] = {0.0}; 
  double qq1_l[3] = {0.0}; 
  double qq2_l[3] = {0.0}; 
  double qq3_l[3] = {0.0}; 
  double qq4_l[3] = {0.0}; 
  double bq0_l[3] = {0.0}; 
  double uq0_l[3] = {0.0}; 
  double bq1_l[3] = {0.0}; 
  double uq1_l[3] = {0.0}; 
  double bq2_l[3] = {0.0}; 
  double uq2_l[3] = {0.0}; 
  rq_l[0] = -(1.1618950038622249*rho_l[3])-0.6708203932499369*rho_l[2]+0.8660254037844386*rho_l[1]+0.5*rho_l[0]; 
  rq_l[1] = 0.8660254037844386*rho_l[1]+0.5*rho_l[0]; 
  rq_l[2] = 1.1618950038622249*rho_l[3]+0.6708203932499369*rho_l[2]+0.8660254037844386*rho_l[1]+0.5*rho_l[0]; 
  xq_l[0] = -(1.1618950038622249*xd_l[3])-0.6708203932499369*xd_l[2]+0.8660254037844386*xd_l[1]+0.5*xd_l[0]; 
  xq_l[1] = 0.8660254037844386*xd_l[1]+0.5*xd_l[0]; 
  xq_l[2] = 1.1618950038622249*xd_l[3]+0.6708203932499369*xd_l[2]+0.8660254037844386*xd_l[1]+0.5*xd_l[0]; 
  thq_l[0] = -(1.1618950038622249*thd_l[3])-0.6708203932499369*thd_l[2]+0.8660254037844386*thd_l[1]+0.5*thd_l[0]; 
  thq_l[1] = 0.8660254037844386*thd_l[1]+0.5*thd_l[0]; 
  thq_l[2] = 1.1618950038622249*thd_l[3]+0.6708203932499369*thd_l[2]+0.8660254037844386*thd_l[1]+0.5*thd_l[0]; 
  qq1_l[0] = -(1.1618950038622249*q1_l[3])-0.6708203932499369*q1_l[2]+0.8660254037844386*q1_l[1]+0.5*q1_l[0]; 
  qq1_l[1] = 0.8660254037844386*q1_l[1]+0.5*q1_l[0]; 
  qq1_l[2] = 1.1618950038622249*q1_l[3]+0.6708203932499369*q1_l[2]+0.8660254037844386*q1_l[1]+0.5*q1_l[0]; 
  qq2_l[0] = -(1.1618950038622249*q2_l[3])-0.6708203932499369*q2_l[2]+0.8660254037844386*q2_l[1]+0.5*q2_l[0]; 
  qq2_l[1] = 0.8660254037844386*q2_l[1]+0.5*q2_l[0]; 
  qq2_l[2] = 1.1618950038622249*q2_l[3]+0.6708203932499369*q2_l[2]+0.8660254037844386*q2_l[1]+0.5*q2_l[0]; 
  qq3_l[0] = -(1.1618950038622249*q3_l[3])-0.6708203932499369*q3_l[2]+0.8660254037844386*q3_l[1]+0.5*q3_l[0]; 
  qq3_l[1] = 0.8660254037844386*q3_l[1]+0.5*q3_l[0]; 
  qq3_l[2] = 1.1618950038622249*q3_l[3]+0.6708203932499369*q3_l[2]+0.8660254037844386*q3_l[1]+0.5*q3_l[0]; 
  qq4_l[0] = -(1.1618950038622249*q4_l[3])-0.6708203932499369*q4_l[2]+0.8660254037844386*q4_l[1]+0.5*q4_l[0]; 
  qq4_l[1] = 0.8660254037844386*q4_l[1]+0.5*q4_l[0]; 
  qq4_l[2] = 1.1618950038622249*q4_l[3]+0.6708203932499369*q4_l[2]+0.8660254037844386*q4_l[1]+0.5*q4_l[0]; 
  bq0_l[0] = -(1.1618950038622249*b0_l[3])-0.6708203932499369*b0_l[2]+0.8660254037844386*b0_l[1]+0.5*b0_l[0]; 
  bq0_l[1] = 0.8660254037844386*b0_l[1]+0.5*b0_l[0]; 
  bq0_l[2] = 1.1618950038622249*b0_l[3]+0.6708203932499369*b0_l[2]+0.8660254037844386*b0_l[1]+0.5*b0_l[0]; 
  uq0_l[0] = -(1.1618950038622249*u0_l[3])-0.6708203932499369*u0_l[2]+0.8660254037844386*u0_l[1]+0.5*u0_l[0]; 
  uq0_l[1] = 0.8660254037844386*u0_l[1]+0.5*u0_l[0]; 
  uq0_l[2] = 1.1618950038622249*u0_l[3]+0.6708203932499369*u0_l[2]+0.8660254037844386*u0_l[1]+0.5*u0_l[0]; 
  bq1_l[0] = -(1.1618950038622249*b1_l[3])-0.6708203932499369*b1_l[2]+0.8660254037844386*b1_l[1]+0.5*b1_l[0]; 
  bq1_l[1] = 0.8660254037844386*b1_l[1]+0.5*b1_l[0]; 
  bq1_l[2] = 1.1618950038622249*b1_l[3]+0.6708203932499369*b1_l[2]+0.8660254037844386*b1_l[1]+0.5*b1_l[0]; 
  uq1_l[0] = -(1.1618950038622249*u1_l[3])-0.6708203932499369*u1_l[2]+0.8660254037844386*u1_l[1]+0.5*u1_l[0]; 
  uq1_l[1] = 0.8660254037844386*u1_l[1]+0.5*u1_l[0]; 
  uq1_l[2] = 1.1618950038622249*u1_l[3]+0.6708203932499369*u1_l[2]+0.8660254037844386*u1_l[1]+0.5*u1_l[0]; 
  bq2_l[0] = -(1.1618950038622249*b2_l[3])-0.6708203932499369*b2_l[2]+0.8660254037844386*b2_l[1]+0.5*b2_l[0]; 
  bq2_l[1] = 0.8660254037844386*b2_l[1]+0.5*b2_l[0]; 
  bq2_l[2] = 1.1618950038622249*b2_l[3]+0.6708203932499369*b2_l[2]+0.8660254037844386*b2_l[1]+0.5*b2_l[0]; 
  uq2_l[0] = -(1.1618950038622249*u2_l[3])-0.6708203932499369*u2_l[2]+0.8660254037844386*u2_l[1]+0.5*u2_l[0]; 
  uq2_l[1] = 0.8660254037844386*u2_l[1]+0.5*u2_l[0]; 
  uq2_l[2] = 1.1618950038622249*u2_l[3]+0.6708203932499369*u2_l[2]+0.8660254037844386*u2_l[1]+0.5*u2_l[0]; 
  double rq_r[3] = {0.0}; 
  double xq_r[3] = {0.0}; 
  double thq_r[3] = {0.0}; 
  double qq1_r[3] = {0.0}; 
  double qq2_r[3] = {0.0}; 
  double qq3_r[3] = {0.0}; 
  double qq4_r[3] = {0.0}; 
  double bq0_r[3] = {0.0}; 
  double uq0_r[3] = {0.0}; 
  double bq1_r[3] = {0.0}; 
  double uq1_r[3] = {0.0}; 
  double bq2_r[3] = {0.0}; 
  double uq2_r[3] = {0.0}; 
  rq_r[0] = 1.1618950038622249*rho_r[3]-0.6708203932499369*rho_r[2]-0.8660254037844386*rho_r[1]+0.5*rho_r[0]; 
  rq_r[1] = 0.5*rho_r[0]-0.8660254037844386*rho_r[1]; 
  rq_r[2] = -(1.1618950038622249*rho_r[3])+0.6708203932499369*rho_r[2]-0.8660254037844386*rho_r[1]+0.5*rho_r[0]; 
  xq_r[0] = 1.1618950038622249*xd_r[3]-0.6708203932499369*xd_r[2]-0.8660254037844386*xd_r[1]+0.5*xd_r[0]; 
  xq_r[1] = 0.5*xd_r[0]-0.8660254037844386*xd_r[1]; 
  xq_r[2] = -(1.1618950038622249*xd_r[3])+0.6708203932499369*xd_r[2]-0.8660254037844386*xd_r[1]+0.5*xd_r[0]; 
  thq_r[0] = 1.1618950038622249*thd_r[3]-0.6708203932499369*thd_r[2]-0.8660254037844386*thd_r[1]+0.5*thd_r[0]; 
  thq_r[1] = 0.5*thd_r[0]-0.8660254037844386*thd_r[1]; 
  thq_r[2] = -(1.1618950038622249*thd_r[3])+0.6708203932499369*thd_r[2]-0.8660254037844386*thd_r[1]+0.5*thd_r[0]; 
  qq1_r[0] = 1.1618950038622249*q1_r[3]-0.6708203932499369*q1_r[2]-0.8660254037844386*q1_r[1]+0.5*q1_r[0]; 
  qq1_r[1] = 0.5*q1_r[0]-0.8660254037844386*q1_r[1]; 
  qq1_r[2] = -(1.1618950038622249*q1_r[3])+0.6708203932499369*q1_r[2]-0.8660254037844386*q1_r[1]+0.5*q1_r[0]; 
  qq2_r[0] = 1.1618950038622249*q2_r[3]-0.6708203932499369*q2_r[2]-0.8660254037844386*q2_r[1]+0.5*q2_r[0]; 
  qq2_r[1] = 0.5*q2_r[0]-0.8660254037844386*q2_r[1]; 
  qq2_r[2] = -(1.1618950038622249*q2_r[3])+0.6708203932499369*q2_r[2]-0.8660254037844386*q2_r[1]+0.5*q2_r[0]; 
  qq3_r[0] = 1.1618950038622249*q3_r[3]-0.6708203932499369*q3_r[2]-0.8660254037844386*q3_r[1]+0.5*q3_r[0]; 
  qq3_r[1] = 0.5*q3_r[0]-0.8660254037844386*q3_r[1]; 
  qq3_r[2] = -(1.1618950038622249*q3_r[3])+0.6708203932499369*q3_r[2]-0.8660254037844386*q3_r[1]+0.5*q3_r[0]; 
  qq4_r[0] = 1.1618950038622249*q4_r[3]-0.6708203932499369*q4_r[2]-0.8660254037844386*q4_r[1]+0.5*q4_r[0]; 
  qq4_r[1] = 0.5*q4_r[0]-0.8660254037844386*q4_r[1]; 
  qq4_r[2] = -(1.1618950038622249*q4_r[3])+0.6708203932499369*q4_r[2]-0.8660254037844386*q4_r[1]+0.5*q4_r[0]; 
  bq0_r[0] = 1.1618950038622249*b0_r[3]-0.6708203932499369*b0_r[2]-0.8660254037844386*b0_r[1]+0.5*b0_r[0]; 
  bq0_r[1] = 0.5*b0_r[0]-0.8660254037844386*b0_r[1]; 
  bq0_r[2] = -(1.1618950038622249*b0_r[3])+0.6708203932499369*b0_r[2]-0.8660254037844386*b0_r[1]+0.5*b0_r[0]; 
  uq0_r[0] = 1.1618950038622249*u0_r[3]-0.6708203932499369*u0_r[2]-0.8660254037844386*u0_r[1]+0.5*u0_r[0]; 
  uq0_r[1] = 0.5*u0_r[0]-0.8660254037844386*u0_r[1]; 
  uq0_r[2] = -(1.1618950038622249*u0_r[3])+0.6708203932499369*u0_r[2]-0.8660254037844386*u0_r[1]+0.5*u0_r[0]; 
  bq1_r[0] = 1.1618950038622249*b1_r[3]-0.6708203932499369*b1_r[2]-0.8660254037844386*b1_r[1]+0.5*b1_r[0]; 
  bq1_r[1] = 0.5*b1_r[0]-0.8660254037844386*b1_r[1]; 
  bq1_r[2] = -(1.1618950038622249*b1_r[3])+0.6708203932499369*b1_r[2]-0.8660254037844386*b1_r[1]+0.5*b1_r[0]; 
  uq1_r[0] = 1.1618950038622249*u1_r[3]-0.6708203932499369*u1_r[2]-0.8660254037844386*u1_r[1]+0.5*u1_r[0]; 
  uq1_r[1] = 0.5*u1_r[0]-0.8660254037844386*u1_r[1]; 
  uq1_r[2] = -(1.1618950038622249*u1_r[3])+0.6708203932499369*u1_r[2]-0.8660254037844386*u1_r[1]+0.5*u1_r[0]; 
  bq2_r[0] = 1.1618950038622249*b2_r[3]-0.6708203932499369*b2_r[2]-0.8660254037844386*b2_r[1]+0.5*b2_r[0]; 
  bq2_r[1] = 0.5*b2_r[0]-0.8660254037844386*b2_r[1]; 
  bq2_r[2] = -(1.1618950038622249*b2_r[3])+0.6708203932499369*b2_r[2]-0.8660254037844386*b2_r[1]+0.5*b2_r[0]; 
  uq2_r[0] = 1.1618950038622249*u2_r[3]-0.6708203932499369*u2_r[2]-0.8660254037844386*u2_r[1]+0.5*u2_r[0]; 
  uq2_r[1] = 0.5*u2_r[0]-0.8660254037844386*u2_r[1]; 
  uq2_r[2] = -(1.1618950038622249*u2_r[3])+0.6708203932499369*u2_r[2]-0.8660254037844386*u2_r[1]+0.5*u2_r[0]; 
 
  // Nodal fluxes [F_rho, F_x, F_y, F_z, F_E]: the two one-sided contributions. 
  double Fq[15] = {0.0}; 
  for (int n = 0; n < 3; ++n) { 
    { 
      const double I0 = I[0 + n], I1 = I[3 + n], I2 = I[6 + n], J0 = I[9 + n]; 
      const double rq = rq_l[n], xq = xq_l[n], x2q = xq*xq; 
      const double pperp = 0.5*jperp*rq*x2q*thq_l[n]; // half: the central average of the two traces 
      cfl_max = fmax(cfl_max, fabs(uq0_l[n]) + xq*sqrt(fmax(2.0*thq_l[n]*(1.0 - bq0_l[n]*bq0_l[n]), 0.0))); 
      const double bn = bq0_l[n]; 
      const double bu = bq0_l[n]*uq0_l[n] + bq1_l[n]*uq1_l[n] + bq2_l[n]*uq2_l[n]; 
      const double uu = uq0_l[n]*uq0_l[n] + uq1_l[n]*uq1_l[n] + uq2_l[n]*uq2_l[n]; 
      Fq[n] += rq*I0; 
      Fq[3 + n] += rq*(uq0_l[n]*I0 + xq*bq0_l[n]*I1) + pperp*(1.0 - bq0_l[n]*bn); 
      Fq[6 + n] += rq*(uq1_l[n]*I0 + xq*bq1_l[n]*I1) + pperp*(0.0 - bq1_l[n]*bn); 
      Fq[9 + n] += rq*(uq2_l[n]*I0 + xq*bq2_l[n]*I1) + pperp*(0.0 - bq2_l[n]*bn); 
      Fq[12 + n] += 0.5*rq*(uu*I0 + 2.0*xq*bu*I1 + x2q*I2) + rq*x2q*J0 + pperp*(uq0_l[n] - bu*bn); 
    } 
    { 
      const double I0 = I[12 + n], I1 = I[15 + n], I2 = I[18 + n], J0 = I[21 + n]; 
      const double rq = rq_r[n], xq = xq_r[n], x2q = xq*xq; 
      const double pperp = 0.5*jperp*rq*x2q*thq_r[n]; // half: the central average of the two traces 
      cfl_max = fmax(cfl_max, fabs(uq0_r[n]) + xq*sqrt(fmax(2.0*thq_r[n]*(1.0 - bq0_r[n]*bq0_r[n]), 0.0))); 
      const double bn = bq0_r[n]; 
      const double bu = bq0_r[n]*uq0_r[n] + bq1_r[n]*uq1_r[n] + bq2_r[n]*uq2_r[n]; 
      const double uu = uq0_r[n]*uq0_r[n] + uq1_r[n]*uq1_r[n] + uq2_r[n]*uq2_r[n]; 
      Fq[n] += rq*I0; 
      Fq[3 + n] += rq*(uq0_r[n]*I0 + xq*bq0_r[n]*I1) + pperp*(1.0 - bq0_r[n]*bn); 
      Fq[6 + n] += rq*(uq1_r[n]*I0 + xq*bq1_r[n]*I1) + pperp*(0.0 - bq1_r[n]*bn); 
      Fq[9 + n] += rq*(uq2_r[n]*I0 + xq*bq2_r[n]*I1) + pperp*(0.0 - bq2_r[n]*bn); 
      Fq[12 + n] += 0.5*rq*(uu*I0 + 2.0*xq*bu*I1 + x2q*I2) + rq*x2q*J0 + pperp*(uq0_r[n] - bu*bn); 
    } 
  } 
  // Rusanov dissipation of the perpendicular sound waves of the central pressure flux. 
  for (int n = 0; n < 3; ++n) { 
    const double s_l = xq_l[n]*sqrt(fmax(2.0*thq_l[n]*(1.0 - bq0_l[n]*bq0_l[n]), 0.0)); 
    const double s_r = xq_r[n]*sqrt(fmax(2.0*thq_r[n]*(1.0 - bq0_r[n]*bq0_r[n]), 0.0)); 
    const double s = 0.5*jperp*fmax(s_l, s_r); 
    Fq[n] -= s*(rq_r[n] - rq_l[n]); 
    Fq[3 + n] -= s*(qq1_r[n] - qq1_l[n]); 
    Fq[6 + n] -= s*(qq2_r[n] - qq2_l[n]); 
    Fq[9 + n] -= s*(qq3_r[n] - qq3_l[n]); 
    Fq[12 + n] -= s*(qq4_r[n] - qq4_l[n]); 
  } 
 
  // Projection of the nodal fluxes onto the surface basis. 
  flux[0] = 0.39283710065919303*Fq[0] + 0.6285393610547091*Fq[1] + 0.39283710065919303*Fq[2]; 
  flux[1] = -0.5270462766947298*Fq[0] + 0.0*Fq[1] + 0.5270462766947298*Fq[2]; 
  flux[2] = 0.39283710065919303*Fq[3] + 0.6285393610547091*Fq[4] + 0.39283710065919303*Fq[5]; 
  flux[3] = -0.5270462766947298*Fq[3] + 0.0*Fq[4] + 0.5270462766947298*Fq[5]; 
  flux[4] = 0.39283710065919303*Fq[6] + 0.6285393610547091*Fq[7] + 0.39283710065919303*Fq[8]; 
  flux[5] = -0.5270462766947298*Fq[6] + 0.0*Fq[7] + 0.5270462766947298*Fq[8]; 
  flux[6] = 0.39283710065919303*Fq[9] + 0.6285393610547091*Fq[10] + 0.39283710065919303*Fq[11]; 
  flux[7] = -0.5270462766947298*Fq[9] + 0.0*Fq[10] + 0.5270462766947298*Fq[11]; 
  flux[8] = 0.39283710065919303*Fq[12] + 0.6285393610547091*Fq[13] + 0.39283710065919303*Fq[14]; 
  flux[9] = -0.5270462766947298*Fq[12] + 0.0*Fq[13] + 0.5270462766947298*Fq[14]; 
  return cfl_max; 
} 
