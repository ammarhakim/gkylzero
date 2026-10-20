#include <gkyl_pkpm_fluid_kernels.h> 
GKYL_CU_DH void pkpm_fluid_flux_x_1x_ser_p1(const double *bvar_l, const double *bvar_r,
  const double *u_l, const double *u_r, const double *xi_l, const double *xi_r,
  const double *fluid_l, const double *fluid_r, const double *I, double* GKYL_RESTRICT flux) 
{ 
  // bvar_l/r:  Magnetic field unit vector [bx, by, bz] of the left/right cell. 
  // u_l/r:     Flow velocity [ux, uy, uz] of the left/right cell. 
  // xi_l/r:    Thermal speed [xi, xi^2, p] of the left/right cell (the first expansion is used: 
  //            xi^2 is the square of its trace, the xi of the streaming speed). 
  // fluid_l/r: PHYSICAL fluid state [rho, rho ux, rho uy, rho uz, E] of the left/right cell (Jacobian divided out). 
  // I:         Nodal moments [I_0, I_1, I_2, I_G, I_M (3)]^l, then ^r, of the kinetic HLL flux at the 
  //            1 face nodes (pkpm_fluid_flux_accum_x_1x1v_ser_p1 summed over the velocity cells); 
  //            they carry the face Jacobian J/J_x. 
  // flux:      Output modal surface fluxes [F_rho, F_x, F_y, F_z, F_E] of the lower face of the right 
  //            cell for the evolved J-weighted fluid. 
  const double *rho_l = &fluid_l[0]; 
  const double *xd_l = &xi_l[0]; 
  const double *b0_l = &bvar_l[0]; 
  const double *u0_l = &u_l[0]; 
  const double *b1_l = &bvar_l[2]; 
  const double *u1_l = &u_l[2]; 
  const double *b2_l = &bvar_l[4]; 
  const double *u2_l = &u_l[4]; 
  const double *rho_r = &fluid_r[0]; 
  const double *xd_r = &xi_r[0]; 
  const double *b0_r = &bvar_r[0]; 
  const double *u0_r = &u_r[0]; 
  const double *b1_r = &bvar_r[2]; 
  const double *u1_r = &u_r[2]; 
  const double *b2_r = &bvar_r[4]; 
  const double *u2_r = &u_r[4]; 
 
  double rq_l[1] = {0.0}; 
  double xq_l[1] = {0.0}; 
  double bq0_l[1] = {0.0}; 
  double uq0_l[1] = {0.0}; 
  double bq1_l[1] = {0.0}; 
  double uq1_l[1] = {0.0}; 
  double bq2_l[1] = {0.0}; 
  double uq2_l[1] = {0.0}; 
  rq_l[0] = 1.224744871391589*rho_l[1]+0.7071067811865475*rho_l[0]; 
  xq_l[0] = 1.224744871391589*xd_l[1]+0.7071067811865475*xd_l[0]; 
  bq0_l[0] = 1.224744871391589*b0_l[1]+0.7071067811865475*b0_l[0]; 
  uq0_l[0] = 1.224744871391589*u0_l[1]+0.7071067811865475*u0_l[0]; 
  bq1_l[0] = 1.224744871391589*b1_l[1]+0.7071067811865475*b1_l[0]; 
  uq1_l[0] = 1.224744871391589*u1_l[1]+0.7071067811865475*u1_l[0]; 
  bq2_l[0] = 1.224744871391589*b2_l[1]+0.7071067811865475*b2_l[0]; 
  uq2_l[0] = 1.224744871391589*u2_l[1]+0.7071067811865475*u2_l[0]; 
  double rq_r[1] = {0.0}; 
  double xq_r[1] = {0.0}; 
  double bq0_r[1] = {0.0}; 
  double uq0_r[1] = {0.0}; 
  double bq1_r[1] = {0.0}; 
  double uq1_r[1] = {0.0}; 
  double bq2_r[1] = {0.0}; 
  double uq2_r[1] = {0.0}; 
  rq_r[0] = 0.7071067811865475*rho_r[0]-1.224744871391589*rho_r[1]; 
  xq_r[0] = 0.7071067811865475*xd_r[0]-1.224744871391589*xd_r[1]; 
  bq0_r[0] = 0.7071067811865475*b0_r[0]-1.224744871391589*b0_r[1]; 
  uq0_r[0] = 0.7071067811865475*u0_r[0]-1.224744871391589*u0_r[1]; 
  bq1_r[0] = 0.7071067811865475*b1_r[0]-1.224744871391589*b1_r[1]; 
  uq1_r[0] = 0.7071067811865475*u1_r[0]-1.224744871391589*u1_r[1]; 
  bq2_r[0] = 0.7071067811865475*b2_r[0]-1.224744871391589*b2_r[1]; 
  uq2_r[0] = 0.7071067811865475*u2_r[0]-1.224744871391589*u2_r[1]; 
 
  // Nodal fluxes [F_rho, F_x, F_y, F_z, F_E]: the moment maps of the two one-sided contributions. 
  double Fq[5] = {0.0}; 
  for (int n = 0; n < 1; ++n) { 
    { 
      const double I0 = I[0 + n], I1 = I[1 + n], I2 = I[2 + n], J0 = I[3 + n]; 
      const double IM[3] = {I[4 + n], I[5 + n], I[6 + n]}; 
      const double rq = rq_l[n], xq = xq_l[n], x2q = xq*xq; 
      const double bu = bq0_l[n]*uq0_l[n] + bq1_l[n]*uq1_l[n] + bq2_l[n]*uq2_l[n]; 
      const double uu = uq0_l[n]*uq0_l[n] + uq1_l[n]*uq1_l[n] + uq2_l[n]*uq2_l[n]; 
      const double uIM = uq0_l[n]*IM[0] + uq1_l[n]*IM[1] + uq2_l[n]*IM[2]; 
      Fq[n] += rq*I0; 
      Fq[1 + n] += rq*(uq0_l[n]*I0 + xq*bq0_l[n]*I1) + rq*xq*IM[0]; 
      Fq[2 + n] += rq*(uq1_l[n]*I0 + xq*bq1_l[n]*I1) + rq*xq*IM[1]; 
      Fq[3 + n] += rq*(uq2_l[n]*I0 + xq*bq2_l[n]*I1) + rq*xq*IM[2]; 
      Fq[4 + n] += 0.5*rq*(uu*I0 + 2.0*xq*bu*I1 + x2q*I2) + rq*x2q*J0 + rq*xq*uIM; 
    } 
    { 
      const double I0 = I[7 + n], I1 = I[8 + n], I2 = I[9 + n], J0 = I[10 + n]; 
      const double IM[3] = {I[11 + n], I[12 + n], I[13 + n]}; 
      const double rq = rq_r[n], xq = xq_r[n], x2q = xq*xq; 
      const double bu = bq0_r[n]*uq0_r[n] + bq1_r[n]*uq1_r[n] + bq2_r[n]*uq2_r[n]; 
      const double uu = uq0_r[n]*uq0_r[n] + uq1_r[n]*uq1_r[n] + uq2_r[n]*uq2_r[n]; 
      const double uIM = uq0_r[n]*IM[0] + uq1_r[n]*IM[1] + uq2_r[n]*IM[2]; 
      Fq[n] += rq*I0; 
      Fq[1 + n] += rq*(uq0_r[n]*I0 + xq*bq0_r[n]*I1) + rq*xq*IM[0]; 
      Fq[2 + n] += rq*(uq1_r[n]*I0 + xq*bq1_r[n]*I1) + rq*xq*IM[1]; 
      Fq[3 + n] += rq*(uq2_r[n]*I0 + xq*bq2_r[n]*I1) + rq*xq*IM[2]; 
      Fq[4 + n] += 0.5*rq*(uu*I0 + 2.0*xq*bu*I1 + x2q*I2) + rq*x2q*J0 + rq*xq*uIM; 
    } 
  } 
 
  // Projection of the nodal fluxes onto the surface basis. 
  flux[0] = 1.0*Fq[0]; 
  flux[1] = 1.0*Fq[1]; 
  flux[2] = 1.0*Fq[2]; 
  flux[3] = 1.0*Fq[3]; 
  flux[4] = 1.0*Fq[4]; 
} 
