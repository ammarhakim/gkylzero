#include <gkyl_pkpm_vars_kernels.h> 
GKYL_CU_DH void pkpm_vars_pressure_proj_1x_ser_p1(const double *vlasov_pkpm_moms, const double *bvar, double* GKYL_RESTRICT p_ij) 
{ 
  // vlasov_pkpm_moms: [rho, p_parallel, p_perp, M1], moments of the kinetic equation. 
  // bvar:             Magnetic field unit vector [bx, by, bz]. 
  // p_ij:             Output p=1 projection of p_perp delta_ij + (p_par - p_perp) b_i b_j 
  //                   [Pxx, Pxy, Pxz, Pyy, Pyz, Pzz] (limiter scratch and output only). 
  const double *ppar = &vlasov_pkpm_moms[2]; 
  const double *pperp = &vlasov_pkpm_moms[4]; 
  const double *bx = &bvar[0]; 
  const double *by = &bvar[2]; 
  const double *bz = &bvar[4]; 
  double *Pxx = &p_ij[0]; 
  double *Pxy = &p_ij[2]; 
  double *Pxz = &p_ij[4]; 
  double *Pyy = &p_ij[6]; 
  double *Pyz = &p_ij[8]; 
  double *Pzz = &p_ij[10]; 
 
  Pxx[0] = -(0.5*pperp[0]*pow(bx[1],2))+0.5*ppar[0]*pow(bx[1],2)-0.5*pperp[0]*pow(bx[0],2)+0.5*ppar[0]*pow(bx[0],2)-1.0*bx[0]*bx[1]*pperp[1]+bx[0]*bx[1]*ppar[1]+pperp[0]; 
  Pxx[1] = -(0.9*pperp[1]*pow(bx[1],2))+0.9*ppar[1]*pow(bx[1],2)-0.5*pperp[1]*pow(bx[0],2)+0.5*ppar[1]*pow(bx[0],2)+pperp[1]-1.0*bx[0]*pperp[0]*bx[1]+bx[0]*ppar[0]*bx[1]; 
  Pxy[0] = -(0.5*bx[0]*by[1]*pperp[1])-0.5*by[0]*bx[1]*pperp[1]+0.5*bx[0]*by[1]*ppar[1]+0.5*by[0]*bx[1]*ppar[1]-0.5*pperp[0]*bx[1]*by[1]+0.5*ppar[0]*bx[1]*by[1]-0.5*bx[0]*by[0]*pperp[0]+0.5*bx[0]*by[0]*ppar[0]; 
  Pxy[1] = -(0.9*bx[1]*by[1]*pperp[1])-0.5*bx[0]*by[0]*pperp[1]+0.9*bx[1]*by[1]*ppar[1]+0.5*bx[0]*by[0]*ppar[1]-0.5*bx[0]*pperp[0]*by[1]+0.5*bx[0]*ppar[0]*by[1]-0.5*by[0]*pperp[0]*bx[1]+0.5*by[0]*ppar[0]*bx[1]; 
  Pxz[0] = -(0.5*bx[0]*bz[1]*pperp[1])-0.5*bz[0]*bx[1]*pperp[1]+0.5*bx[0]*bz[1]*ppar[1]+0.5*bz[0]*bx[1]*ppar[1]-0.5*pperp[0]*bx[1]*bz[1]+0.5*ppar[0]*bx[1]*bz[1]-0.5*bx[0]*bz[0]*pperp[0]+0.5*bx[0]*bz[0]*ppar[0]; 
  Pxz[1] = -(0.9*bx[1]*bz[1]*pperp[1])-0.5*bx[0]*bz[0]*pperp[1]+0.9*bx[1]*bz[1]*ppar[1]+0.5*bx[0]*bz[0]*ppar[1]-0.5*bx[0]*pperp[0]*bz[1]+0.5*bx[0]*ppar[0]*bz[1]-0.5*bz[0]*pperp[0]*bx[1]+0.5*bz[0]*ppar[0]*bx[1]; 
  Pyy[0] = -(0.5*pperp[0]*pow(by[1],2))+0.5*ppar[0]*pow(by[1],2)-0.5*pperp[0]*pow(by[0],2)+0.5*ppar[0]*pow(by[0],2)-1.0*by[0]*by[1]*pperp[1]+by[0]*by[1]*ppar[1]+pperp[0]; 
  Pyy[1] = -(0.9*pperp[1]*pow(by[1],2))+0.9*ppar[1]*pow(by[1],2)-0.5*pperp[1]*pow(by[0],2)+0.5*ppar[1]*pow(by[0],2)+pperp[1]-1.0*by[0]*pperp[0]*by[1]+by[0]*ppar[0]*by[1]; 
  Pyz[0] = -(0.5*by[0]*bz[1]*pperp[1])-0.5*bz[0]*by[1]*pperp[1]+0.5*by[0]*bz[1]*ppar[1]+0.5*bz[0]*by[1]*ppar[1]-0.5*pperp[0]*by[1]*bz[1]+0.5*ppar[0]*by[1]*bz[1]-0.5*by[0]*bz[0]*pperp[0]+0.5*by[0]*bz[0]*ppar[0]; 
  Pyz[1] = -(0.9*by[1]*bz[1]*pperp[1])-0.5*by[0]*bz[0]*pperp[1]+0.9*by[1]*bz[1]*ppar[1]+0.5*by[0]*bz[0]*ppar[1]-0.5*by[0]*pperp[0]*bz[1]+0.5*by[0]*ppar[0]*bz[1]-0.5*bz[0]*pperp[0]*by[1]+0.5*bz[0]*ppar[0]*by[1]; 
  Pzz[0] = -(0.5*pperp[0]*pow(bz[1],2))+0.5*ppar[0]*pow(bz[1],2)-0.5*pperp[0]*pow(bz[0],2)+0.5*ppar[0]*pow(bz[0],2)-1.0*bz[0]*bz[1]*pperp[1]+bz[0]*bz[1]*ppar[1]+pperp[0]; 
  Pzz[1] = -(0.9*pperp[1]*pow(bz[1],2))+0.9*ppar[1]*pow(bz[1],2)-0.5*pperp[1]*pow(bz[0],2)+0.5*ppar[1]*pow(bz[0],2)+pperp[1]-1.0*bz[0]*pperp[0]*bz[1]+bz[0]*ppar[0]*bz[1]; 
} 
