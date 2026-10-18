#include <gkyl_pkpm_fluid_kernels.h> 
GKYL_CU_DH void pkpm_fluid_flux_absorb_accum_x_1x1v_ser_p1(int edge, const double *w, const double *dxv,
  const double *bvar, const double *u, const double *xi, const double *f,
  double* GKYL_RESTRICT I) 
{ 
  // edge:      -1 for the lower boundary (face at z_x = -1 of the skin cell), +1 for the upper. 
  // w[NDIM]:   Cell-center coordinates of the phase-space cell. 
  // dxv[NDIM]: Cell spacing of the phase-space cell. 
  // bvar:      Magnetic field unit vector [bx, by, bz] of the skin cell. 
  // u:         Flow velocity [ux, uy, uz] of the skin cell. 
  // xi:        Thermal speed [xi, xi^2, p] of the skin cell (the first expansion is used). 
  // f:         Normalized distribution functions [F, G] of the phase-space cell. 
  // I:         Accumulated nodal moments of the upwind flux at the 1 transverse face nodes: 
  //            I[k*1 + n], k = 0..3 for [int Phi_F dw, int w Phi_F dw, int w^2 Phi_F dw, int Phi_G dw]. 
  const double ze = edge < 0 ? -1.0 : 1.0; 
  const double wv = w[1], dv = dxv[1]; 
  const double *bd = &bvar[0]; 
  const double *ud = &u[0]; 
  const double *xd = &xi[0]; 
  const double *F = &f[0]; 
  const double *G = &f[6]; 
  const double vz[3] = {-0.7745966692414833, 0.0, 0.7745966692414833}; 
  const double vw[3] = {0.5555555555555556, 0.8888888888888888, 0.5555555555555556}; 
 
  double bq[1] = {0.0}; 
  double uq[1] = {0.0}; 
  double xq[1] = {0.0}; 
  bq[0] = 1.224744871391589*bd[1]*ze+0.7071067811865475*bd[0]; 
  uq[0] = 1.224744871391589*ud[1]*ze+0.7071067811865475*ud[0]; 
  xq[0] = 1.224744871391589*xd[1]*ze+0.7071067811865475*xd[0]; 
 
  double Fq[3] = {0.0}; 
  double Gq[3] = {0.0}; 
  Fq[0] = 0.7745966692414833*F[5]*ze-1.1618950038622249*F[3]*ze+0.8660254037844386*F[1]*ze+0.4472135954999579*F[4]-0.6708203932499369*F[2]+0.5*F[0]; 
  Gq[0] = 0.7745966692414833*G[5]*ze-1.1618950038622249*G[3]*ze+0.8660254037844386*G[1]*ze+0.4472135954999579*G[4]-0.6708203932499369*G[2]+0.5*G[0]; 
  Fq[1] = -(0.9682458365518543*F[5]*ze)+0.8660254037844386*F[1]*ze-0.5590169943749475*F[4]+0.5*F[0]; 
  Gq[1] = -(0.9682458365518543*G[5]*ze)+0.8660254037844386*G[1]*ze-0.5590169943749475*G[4]+0.5*G[0]; 
  Fq[2] = 0.7745966692414833*F[5]*ze+1.1618950038622249*F[3]*ze+0.8660254037844386*F[1]*ze+0.4472135954999579*F[4]+0.6708203932499369*F[2]+0.5*F[0]; 
  Gq[2] = 0.7745966692414833*G[5]*ze+1.1618950038622249*G[3]*ze+0.8660254037844386*G[1]*ze+0.4472135954999579*G[4]+0.6708203932499369*G[2]+0.5*G[0]; 
 
  for (int n = 0; n < 1; ++n) { 
    for (int j = 0; j < 3; ++j) { 
      const double wq = wv + 0.5*dv*vz[j]; 
      const double a = xq[n]*bq[n]*wq + uq[n]; 
      const double phi = (edge < 0 ? fmin(a, 0.0) : fmax(a, 0.0))*0.5*dv*vw[j]; 
      const double fq = Fq[n*3 + j], gq = Gq[n*3 + j]; 
      I[0*1 + n] += phi*fq; 
      I[1*1 + n] += wq*phi*fq; 
      I[2*1 + n] += wq*wq*phi*fq; 
      I[3*1 + n] += phi*gq; 
    } 
  } 
} 
