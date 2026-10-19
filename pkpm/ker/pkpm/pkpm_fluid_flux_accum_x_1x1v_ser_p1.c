#include <gkyl_pkpm_fluid_kernels.h> 
GKYL_CU_DH double pkpm_fluid_flux_accum_x_1x1v_ser_p1(const double *w, const double *dxv,
  const double *jacob_pos_l, const double *jacob_pos_r, const double *bvar_l, const double *bvar_r,
  const double *u_l, const double *u_r, const double *xi_l, const double *xi_r,
  const double *f_l, const double *f_r, double* GKYL_RESTRICT I) 
{ 
  // w[NDIM]:   Cell-center coordinates of the phase-space cell (r side). 
  // dxv[NDIM]: Cell spacing of the phase-space cell. 
  // jacob_pos_l/r: Position-map Jacobian per configuration direction (cdim blocks of 2 equal 
  //            entries) of the left/right cell. 
  // bvar_l/r:  Magnetic field unit vector [bx, by, bz] of the left/right cell. 
  // u_l/r:     Flow velocity [ux, uy, uz] of the left/right cell. 
  // xi_l/r:    Thermal speed [xi, xi^2, p] of the left/right cell (the first expansion is used). 
  // f_l/r:     Jacobian-weighted normalized distribution functions J [F, G] of the left/right phase-space cell. 
  // I:         Accumulated nodal moments of the kinetic HLL flux at the 1 transverse face nodes: 
  //            I[k*1 + n], k = 0..3 for [int A_l F_l dw, int w A_l F_l dw, int w^2 A_l F_l dw, int A_l G_l dw] 
  //            and k = 4..7 the same with the right trace and weight A_r; the traces are divided by 
  //            J_x of their cell, so the moments carry the face Jacobian J/J_x. 
  // Returns the largest speed bound max(-S_L, S_R) over the nodes (for the CFL rate). 
  const double wv = w[1], dv = dxv[1]; 
  double alpha_max = 0.0; 
  const double jdinv_l = 1.0/jacob_pos_l[0], jdinv_r = 1.0/jacob_pos_r[0]; 
  const double *bdl = &bvar_l[0]; 
  const double *bdr = &bvar_r[0]; 
  const double *udl = &u_l[0]; 
  const double *udr = &u_r[0]; 
  const double *xdl = &xi_l[0]; 
  const double *xdr = &xi_r[0]; 
  const double *Fl = &f_l[0]; 
  const double *Gl = &f_l[6]; 
  const double *Fr = &f_r[0]; 
  const double *Gr = &f_r[6]; 
  const double vz[3] = {-0.7745966692414833, 0.0, 0.7745966692414833}; 
  const double vw[3] = {0.5555555555555556, 0.8888888888888888, 0.5555555555555556}; 
 
  double bql[1] = {0.0}; 
  double uql[1] = {0.0}; 
  double xql[1] = {0.0}; 
  double bqr[1] = {0.0}; 
  double uqr[1] = {0.0}; 
  double xqr[1] = {0.0}; 
  bql[0] = 1.224744871391589*bdl[1]+0.7071067811865475*bdl[0]; 
  uql[0] = 1.224744871391589*udl[1]+0.7071067811865475*udl[0]; 
  xql[0] = 1.224744871391589*xdl[1]+0.7071067811865475*xdl[0]; 
  bqr[0] = 0.7071067811865475*bdr[0]-1.224744871391589*bdr[1]; 
  uqr[0] = 0.7071067811865475*udr[0]-1.224744871391589*udr[1]; 
  xqr[0] = 0.7071067811865475*xdr[0]-1.224744871391589*xdr[1]; 
 
  double Flq[3] = {0.0}; 
  double Glq[3] = {0.0}; 
  double Frq[3] = {0.0}; 
  double Grq[3] = {0.0}; 
  Flq[0] = 0.7745966692414833*Fl[5]+0.4472135954999579*Fl[4]-1.1618950038622249*Fl[3]-0.6708203932499369*Fl[2]+0.8660254037844386*Fl[1]+0.5*Fl[0]; 
  Glq[0] = 0.7745966692414833*Gl[5]+0.4472135954999579*Gl[4]-1.1618950038622249*Gl[3]-0.6708203932499369*Gl[2]+0.8660254037844386*Gl[1]+0.5*Gl[0]; 
  Frq[0] = -(0.7745966692414833*Fr[5])+0.4472135954999579*Fr[4]+1.1618950038622249*Fr[3]-0.6708203932499369*Fr[2]-0.8660254037844386*Fr[1]+0.5*Fr[0]; 
  Grq[0] = -(0.7745966692414833*Gr[5])+0.4472135954999579*Gr[4]+1.1618950038622249*Gr[3]-0.6708203932499369*Gr[2]-0.8660254037844386*Gr[1]+0.5*Gr[0]; 
  Flq[1] = -(0.9682458365518543*Fl[5])-0.5590169943749475*Fl[4]+0.8660254037844386*Fl[1]+0.5*Fl[0]; 
  Glq[1] = -(0.9682458365518543*Gl[5])-0.5590169943749475*Gl[4]+0.8660254037844386*Gl[1]+0.5*Gl[0]; 
  Frq[1] = 0.9682458365518543*Fr[5]-0.5590169943749475*Fr[4]-0.8660254037844386*Fr[1]+0.5*Fr[0]; 
  Grq[1] = 0.9682458365518543*Gr[5]-0.5590169943749475*Gr[4]-0.8660254037844386*Gr[1]+0.5*Gr[0]; 
  Flq[2] = 0.7745966692414833*Fl[5]+0.4472135954999579*Fl[4]+1.1618950038622249*Fl[3]+0.6708203932499369*Fl[2]+0.8660254037844386*Fl[1]+0.5*Fl[0]; 
  Glq[2] = 0.7745966692414833*Gl[5]+0.4472135954999579*Gl[4]+1.1618950038622249*Gl[3]+0.6708203932499369*Gl[2]+0.8660254037844386*Gl[1]+0.5*Gl[0]; 
  Frq[2] = -(0.7745966692414833*Fr[5])+0.4472135954999579*Fr[4]-1.1618950038622249*Fr[3]+0.6708203932499369*Fr[2]-0.8660254037844386*Fr[1]+0.5*Fr[0]; 
  Grq[2] = -(0.7745966692414833*Gr[5])+0.4472135954999579*Gr[4]-1.1618950038622249*Gr[3]+0.6708203932499369*Gr[2]-0.8660254037844386*Gr[1]+0.5*Gr[0]; 
 
  for (int n = 0; n < 1; ++n) { 
    for (int j = 0; j < 3; ++j) { 
      const double wq = wv + 0.5*dv*vz[j]; 
      const double a_l = xql[n]*bql[n]*wq + uql[n]; 
      const double a_r = xqr[n]*bqr[n]*wq + uqr[n]; 
      // HLL with S_L = min(a_l, a_r, 0), S_R = max(a_l, a_r, 0): Phi = A_l F_l + A_r F_r. 
      const double s_l = fmin(fmin(a_l, a_r), 0.0); 
      const double s_r = fmax(fmax(a_l, a_r), 0.0); 
      const double ds = s_r - s_l; 
      alpha_max = fmax(alpha_max, ds); 
      const double A_l = ds > 0.0 ? (s_r*a_l - s_l*s_r)/ds : 0.0; 
      const double A_r = ds > 0.0 ? (s_l*s_r - s_l*a_r)/ds : 0.0; 
      const double phi_l = A_l*0.5*dv*vw[j]*jdinv_l; 
      const double phi_r = A_r*0.5*dv*vw[j]*jdinv_r; 
      const double flq = Flq[n*3 + j], glq = Glq[n*3 + j]; 
      const double frq = Frq[n*3 + j], grq = Grq[n*3 + j]; 
      I[0*1 + n] += phi_l*flq; 
      I[1*1 + n] += wq*phi_l*flq; 
      I[2*1 + n] += wq*wq*phi_l*flq; 
      I[3*1 + n] += phi_l*glq; 
      I[4*1 + n] += phi_r*frq; 
      I[5*1 + n] += wq*phi_r*frq; 
      I[6*1 + n] += wq*wq*phi_r*frq; 
      I[7*1 + n] += phi_r*grq; 
    } 
  } 
  return alpha_max; 
} 
