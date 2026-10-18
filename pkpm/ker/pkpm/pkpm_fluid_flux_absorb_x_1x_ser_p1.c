#include <gkyl_pkpm_fluid_kernels.h> 
GKYL_CU_DH void pkpm_fluid_flux_absorb_x_1x_ser_p1(int edge, const double *bvar, const double *u,
  const double *xi, const double *fluid, const double *I, double* GKYL_RESTRICT flux) 
{ 
  // edge:  -1 for the lower boundary (face at z_x = -1 of the skin cell), +1 for the upper. 
  // bvar:  Magnetic field unit vector [bx, by, bz] of the skin cell. 
  // u:     Flow velocity [ux, uy, uz] of the skin cell. 
  // xi:    Thermal speed [xi, xi^2, p] of the skin cell. 
  // fluid: Fluid state [rho, rho ux, rho uy, rho uz, E] of the skin cell. 
  // I:     Nodal moments [I_0, I_1, I_2, J_0] of the kinetic upwind flux at the 1 face nodes 
  //        (pkpm_fluid_flux_absorb_accum_x_1x1v_ser_p1 summed over the velocity cells). 
  // flux:  Output modal surface fluxes [F_rho, F_x, F_y, F_z, F_E] of the face. 
  const double ze = edge < 0 ? -1.0 : 1.0; 
  const double *rho = &fluid[0]; 
  const double *xd = &xi[0]; 
  const double *x2d = &xi[2]; 
  const double *b0 = &bvar[0]; 
  const double *u0 = &u[0]; 
  const double *b1 = &bvar[2]; 
  const double *u1 = &u[2]; 
  const double *b2 = &bvar[4]; 
  const double *u2 = &u[4]; 
 
  double rq[1] = {0.0}; 
  double xq[1] = {0.0}; 
  double x2q[1] = {0.0}; 
  double bq0[1] = {0.0}; 
  double uq0[1] = {0.0}; 
  double bq1[1] = {0.0}; 
  double uq1[1] = {0.0}; 
  double bq2[1] = {0.0}; 
  double uq2[1] = {0.0}; 
  rq[0] = 1.224744871391589*rho[1]*ze+0.7071067811865475*rho[0]; 
  xq[0] = 1.224744871391589*xd[1]*ze+0.7071067811865475*xd[0]; 
  x2q[0] = 1.224744871391589*x2d[1]*ze+0.7071067811865475*x2d[0]; 
  bq0[0] = 1.224744871391589*b0[1]*ze+0.7071067811865475*b0[0]; 
  uq0[0] = 1.224744871391589*u0[1]*ze+0.7071067811865475*u0[0]; 
  bq1[0] = 1.224744871391589*b1[1]*ze+0.7071067811865475*b1[0]; 
  uq1[0] = 1.224744871391589*u1[1]*ze+0.7071067811865475*u1[0]; 
  bq2[0] = 1.224744871391589*b2[1]*ze+0.7071067811865475*b2[0]; 
  uq2[0] = 1.224744871391589*u2[1]*ze+0.7071067811865475*u2[0]; 
 
  // Nodal fluxes [F_rho, F_x, F_y, F_z, F_E]. 
  double Fq[5] = {0.0}; 
  for (int n = 0; n < 1; ++n) { 
    const double I0 = I[n], I1 = I[1 + n], I2 = I[2 + n], J0 = I[3 + n]; 
    const double bn = bq0[n]; 
    const double bu = bq0[n]*uq0[n] + bq1[n]*uq1[n] + bq2[n]*uq2[n]; 
    const double uu = uq0[n]*uq0[n] + uq1[n]*uq1[n] + uq2[n]*uq2[n]; 
    Fq[n] = rq[n]*I0; 
    Fq[1 + n] = rq[n]*(uq0[n]*I0 + xq[n]*bq0[n]*I1) + rq[n]*x2q[n]*(1.0 - bq0[n]*bn)*J0; 
    Fq[2 + n] = rq[n]*(uq1[n]*I0 + xq[n]*bq1[n]*I1) + rq[n]*x2q[n]*(0.0 - bq1[n]*bn)*J0; 
    Fq[3 + n] = rq[n]*(uq2[n]*I0 + xq[n]*bq2[n]*I1) + rq[n]*x2q[n]*(0.0 - bq2[n]*bn)*J0; 
    Fq[4 + n] = 0.5*rq[n]*(uu*I0 + 2.0*xq[n]*bu*I1 + x2q[n]*I2) + rq[n]*x2q[n]*J0; 
  } 
 
  // Projection of the nodal fluxes onto the surface basis. 
  flux[0] = 1.0*Fq[0]; 
  flux[1] = 1.0*Fq[1]; 
  flux[2] = 1.0*Fq[2]; 
  flux[3] = 1.0*Fq[3]; 
  flux[4] = 1.0*Fq[4]; 
} 
