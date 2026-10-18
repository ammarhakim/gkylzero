#include <gkyl_pkpm_fluid_kernels.h> 
GKYL_CU_DH void pkpm_fluid_flux_absorb_y_2x_ser_p1(int edge, const double *bvar, const double *u,
  const double *xi, const double *fluid, const double *I, double* GKYL_RESTRICT flux) 
{ 
  // edge:  -1 for the lower boundary (face at z_y = -1 of the skin cell), +1 for the upper. 
  // bvar:  Magnetic field unit vector [bx, by, bz] of the skin cell. 
  // u:     Flow velocity [ux, uy, uz] of the skin cell. 
  // xi:    Thermal speed [xi, xi^2, p] of the skin cell. 
  // fluid: Fluid state [rho, rho ux, rho uy, rho uz, E] of the skin cell. 
  // I:     Nodal moments [I_0, I_1, I_2, J_0] of the kinetic upwind flux at the 3 face nodes 
  //        (pkpm_fluid_flux_absorb_accum_y_2x1v_ser_p1 summed over the velocity cells). 
  // flux:  Output modal surface fluxes [F_rho, F_x, F_y, F_z, F_E] of the face. 
  const double ze = edge < 0 ? -1.0 : 1.0; 
  const double *rho = &fluid[0]; 
  const double *xd = &xi[0]; 
  const double *x2d = &xi[4]; 
  const double *b0 = &bvar[0]; 
  const double *u0 = &u[0]; 
  const double *b1 = &bvar[4]; 
  const double *u1 = &u[4]; 
  const double *b2 = &bvar[8]; 
  const double *u2 = &u[8]; 
 
  double rq[3] = {0.0}; 
  double xq[3] = {0.0}; 
  double x2q[3] = {0.0}; 
  double bq0[3] = {0.0}; 
  double uq0[3] = {0.0}; 
  double bq1[3] = {0.0}; 
  double uq1[3] = {0.0}; 
  double bq2[3] = {0.0}; 
  double uq2[3] = {0.0}; 
  rq[0] = -(1.1618950038622249*rho[3]*ze)+0.8660254037844386*rho[2]*ze-0.6708203932499369*rho[1]+0.5*rho[0]; 
  rq[1] = 0.8660254037844386*rho[2]*ze+0.5*rho[0]; 
  rq[2] = 1.1618950038622249*rho[3]*ze+0.8660254037844386*rho[2]*ze+0.6708203932499369*rho[1]+0.5*rho[0]; 
  xq[0] = -(1.1618950038622249*xd[3]*ze)+0.8660254037844386*xd[2]*ze-0.6708203932499369*xd[1]+0.5*xd[0]; 
  xq[1] = 0.8660254037844386*xd[2]*ze+0.5*xd[0]; 
  xq[2] = 1.1618950038622249*xd[3]*ze+0.8660254037844386*xd[2]*ze+0.6708203932499369*xd[1]+0.5*xd[0]; 
  x2q[0] = -(1.1618950038622249*x2d[3]*ze)+0.8660254037844386*x2d[2]*ze-0.6708203932499369*x2d[1]+0.5*x2d[0]; 
  x2q[1] = 0.8660254037844386*x2d[2]*ze+0.5*x2d[0]; 
  x2q[2] = 1.1618950038622249*x2d[3]*ze+0.8660254037844386*x2d[2]*ze+0.6708203932499369*x2d[1]+0.5*x2d[0]; 
  bq0[0] = -(1.1618950038622249*b0[3]*ze)+0.8660254037844386*b0[2]*ze-0.6708203932499369*b0[1]+0.5*b0[0]; 
  bq0[1] = 0.8660254037844386*b0[2]*ze+0.5*b0[0]; 
  bq0[2] = 1.1618950038622249*b0[3]*ze+0.8660254037844386*b0[2]*ze+0.6708203932499369*b0[1]+0.5*b0[0]; 
  uq0[0] = -(1.1618950038622249*u0[3]*ze)+0.8660254037844386*u0[2]*ze-0.6708203932499369*u0[1]+0.5*u0[0]; 
  uq0[1] = 0.8660254037844386*u0[2]*ze+0.5*u0[0]; 
  uq0[2] = 1.1618950038622249*u0[3]*ze+0.8660254037844386*u0[2]*ze+0.6708203932499369*u0[1]+0.5*u0[0]; 
  bq1[0] = -(1.1618950038622249*b1[3]*ze)+0.8660254037844386*b1[2]*ze-0.6708203932499369*b1[1]+0.5*b1[0]; 
  bq1[1] = 0.8660254037844386*b1[2]*ze+0.5*b1[0]; 
  bq1[2] = 1.1618950038622249*b1[3]*ze+0.8660254037844386*b1[2]*ze+0.6708203932499369*b1[1]+0.5*b1[0]; 
  uq1[0] = -(1.1618950038622249*u1[3]*ze)+0.8660254037844386*u1[2]*ze-0.6708203932499369*u1[1]+0.5*u1[0]; 
  uq1[1] = 0.8660254037844386*u1[2]*ze+0.5*u1[0]; 
  uq1[2] = 1.1618950038622249*u1[3]*ze+0.8660254037844386*u1[2]*ze+0.6708203932499369*u1[1]+0.5*u1[0]; 
  bq2[0] = -(1.1618950038622249*b2[3]*ze)+0.8660254037844386*b2[2]*ze-0.6708203932499369*b2[1]+0.5*b2[0]; 
  bq2[1] = 0.8660254037844386*b2[2]*ze+0.5*b2[0]; 
  bq2[2] = 1.1618950038622249*b2[3]*ze+0.8660254037844386*b2[2]*ze+0.6708203932499369*b2[1]+0.5*b2[0]; 
  uq2[0] = -(1.1618950038622249*u2[3]*ze)+0.8660254037844386*u2[2]*ze-0.6708203932499369*u2[1]+0.5*u2[0]; 
  uq2[1] = 0.8660254037844386*u2[2]*ze+0.5*u2[0]; 
  uq2[2] = 1.1618950038622249*u2[3]*ze+0.8660254037844386*u2[2]*ze+0.6708203932499369*u2[1]+0.5*u2[0]; 
 
  // Nodal fluxes [F_rho, F_x, F_y, F_z, F_E]. 
  double Fq[15] = {0.0}; 
  for (int n = 0; n < 3; ++n) { 
    const double I0 = I[n], I1 = I[3 + n], I2 = I[6 + n], J0 = I[9 + n]; 
    const double bn = bq1[n]; 
    const double bu = bq0[n]*uq0[n] + bq1[n]*uq1[n] + bq2[n]*uq2[n]; 
    const double uu = uq0[n]*uq0[n] + uq1[n]*uq1[n] + uq2[n]*uq2[n]; 
    Fq[n] = rq[n]*I0; 
    Fq[3 + n] = rq[n]*(uq0[n]*I0 + xq[n]*bq0[n]*I1) + rq[n]*x2q[n]*(0.0 - bq0[n]*bn)*J0; 
    Fq[6 + n] = rq[n]*(uq1[n]*I0 + xq[n]*bq1[n]*I1) + rq[n]*x2q[n]*(1.0 - bq1[n]*bn)*J0; 
    Fq[9 + n] = rq[n]*(uq2[n]*I0 + xq[n]*bq2[n]*I1) + rq[n]*x2q[n]*(0.0 - bq2[n]*bn)*J0; 
    Fq[12 + n] = 0.5*rq[n]*(uu*I0 + 2.0*xq[n]*bu*I1 + x2q[n]*I2) + rq[n]*x2q[n]*J0; 
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
} 
