#include <gkyl_pkpm_fluid_kernels.h> 
GKYL_CU_DH double pkpm_fluid_vol_1x_ser_p1(const double *w, const double *dxv,
  const double *moms, const double *bvar, const double *u, const double *euler,
  double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates. 
  // dxv[NDIM]: Cell spacing. 
  // moms:      [rho, p_parallel, p_perp, M1], moments of the kinetic equation. 
  // bvar:      Magnetic field unit vector [bx, by, bz]. 
  // u:         Flow velocity [ux, uy, uz]. 
  // euler:     Momentum [rho ux, rho uy, rho uz]. 
  // out:       Incremented momentum increment. 
  const double dx10 = 2.0/dxv[0]; 
  const double *ppar = &moms[2]; 
  const double *pperp = &moms[4]; 
  const double *bx = &bvar[0]; 
  const double *ux = &u[0]; 
  const double *Mx = &euler[0]; 
  double *out0 = &out[0]; 
  const double *by = &bvar[2]; 
  const double *uy = &u[2]; 
  const double *My = &euler[2]; 
  double *out1 = &out[2]; 
  const double *bz = &bvar[4]; 
  const double *uz = &u[4]; 
  const double *Mz = &euler[4]; 
  double *out2 = &out[4]; 
 
  out0[1] += (-(0.8660254037844386*pperp[0]*pow(bx[1],2))+0.8660254037844386*ppar[0]*pow(bx[1],2)-0.8660254037844386*pperp[0]*pow(bx[0],2)+0.8660254037844386*ppar[0]*pow(bx[0],2)+1.224744871391589*Mx[1]*ux[1]-1.7320508075688772*bx[0]*bx[1]*pperp[1]+1.7320508075688772*bx[0]*bx[1]*ppar[1]+1.224744871391589*Mx[0]*ux[0]+1.7320508075688772*pperp[0])*dx10; 
  out1[1] += (1.224744871391589*My[1]*ux[1]-0.8660254037844386*bx[0]*by[1]*pperp[1]-0.8660254037844386*by[0]*bx[1]*pperp[1]+0.8660254037844386*bx[0]*by[1]*ppar[1]+0.8660254037844386*by[0]*bx[1]*ppar[1]-0.8660254037844386*pperp[0]*bx[1]*by[1]+0.8660254037844386*ppar[0]*bx[1]*by[1]+1.224744871391589*My[0]*ux[0]-0.8660254037844386*bx[0]*by[0]*pperp[0]+0.8660254037844386*bx[0]*by[0]*ppar[0])*dx10; 
  out2[1] += (1.224744871391589*Mz[1]*ux[1]-0.8660254037844386*bx[0]*bz[1]*pperp[1]-0.8660254037844386*bz[0]*bx[1]*pperp[1]+0.8660254037844386*bx[0]*bz[1]*ppar[1]+0.8660254037844386*bz[0]*bx[1]*ppar[1]-0.8660254037844386*pperp[0]*bx[1]*bz[1]+0.8660254037844386*ppar[0]*bx[1]*bz[1]+1.224744871391589*Mz[0]*ux[0]-0.8660254037844386*bx[0]*bz[0]*pperp[0]+0.8660254037844386*bx[0]*bz[0]*ppar[0])*dx10; 
 
  double cflFreq_mid = 0.0; 
  cflFreq_mid += 1.5*dx10*fabs(0.7071067811865475*ux[0]); 
  return cflFreq_mid; 
} 
