#include <gkyl_pkpm_fluid_kernels.h> 
GKYL_CU_DH double pkpm_fluid_vol_1x_ser_p1(const double *w, const double *dxv,
  const double *vlasov_pkpm_moms, const double *bvar, const double *u, const double *fluid,
  double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:          Cell-center coordinates. 
  // dxv[NDIM]:        Cell spacing. 
  // vlasov_pkpm_moms: Kinetic coupling moments [p_par, p_perp, q_par, q_perp]. 
  // bvar:             Magnetic field unit vector [bx, by, bz]. 
  // u:                Flow velocity [ux, uy, uz]. 
  // fluid:            Evolved fluid state [rho, rho ux, rho uy, rho uz, E]. 
  // out:              Incremented increment of [rho, rho ux, rho uy, rho uz, E]. 
  const double dx10 = 2.0/dxv[0]; 
  const double *ppar = &vlasov_pkpm_moms[0]; 
  const double *pperp = &vlasov_pkpm_moms[2]; 
  const double *qpar = &vlasov_pkpm_moms[4]; 
  const double *qperp = &vlasov_pkpm_moms[6]; 
  const double *E = &fluid[8]; 
  const double *bx = &bvar[0]; 
  const double *ux = &u[0]; 
  const double *Mx = &fluid[2]; 
  double *out1 = &out[2]; 
  const double *by = &bvar[2]; 
  const double *uy = &u[2]; 
  const double *My = &fluid[4]; 
  double *out2 = &out[4]; 
  const double *bz = &bvar[4]; 
  const double *uz = &u[4]; 
  const double *Mz = &fluid[6]; 
  double *out3 = &out[6]; 
  double *out0 = &out[0]; 
  double *out4 = &out[8]; 
 
  out0[1] += 1.7320508075688772*Mx[0]*dx10; 
  out1[1] += (-(0.8660254037844386*pperp[0]*pow(bx[1],2))+0.8660254037844386*ppar[0]*pow(bx[1],2)-0.8660254037844386*pperp[0]*pow(bx[0],2)+0.8660254037844386*ppar[0]*pow(bx[0],2)+1.224744871391589*Mx[1]*ux[1]-1.7320508075688772*bx[0]*bx[1]*pperp[1]+1.7320508075688772*bx[0]*bx[1]*ppar[1]+1.224744871391589*Mx[0]*ux[0]+1.7320508075688772*pperp[0])*dx10; 
  out2[1] += (1.224744871391589*My[1]*ux[1]-0.8660254037844386*bx[0]*by[1]*pperp[1]-0.8660254037844386*by[0]*bx[1]*pperp[1]+0.8660254037844386*bx[0]*by[1]*ppar[1]+0.8660254037844386*by[0]*bx[1]*ppar[1]-0.8660254037844386*pperp[0]*bx[1]*by[1]+0.8660254037844386*ppar[0]*bx[1]*by[1]+1.224744871391589*My[0]*ux[0]-0.8660254037844386*bx[0]*by[0]*pperp[0]+0.8660254037844386*bx[0]*by[0]*ppar[0])*dx10; 
  out3[1] += (1.224744871391589*Mz[1]*ux[1]-0.8660254037844386*bx[0]*bz[1]*pperp[1]-0.8660254037844386*bz[0]*bx[1]*pperp[1]+0.8660254037844386*bx[0]*bz[1]*ppar[1]+0.8660254037844386*bz[0]*bx[1]*ppar[1]-0.8660254037844386*pperp[0]*bx[1]*bz[1]+0.8660254037844386*ppar[0]*bx[1]*bz[1]+1.224744871391589*Mz[0]*ux[0]-0.8660254037844386*bx[0]*bz[0]*pperp[0]+0.8660254037844386*bx[0]*bz[0]*ppar[0])*dx10; 
  out4[1] += (-(1.1022703842524295*pperp[1]*ux[1]*pow(bx[1],2))+1.1022703842524295*ppar[1]*ux[1]*pow(bx[1],2)-0.6123724356957944*pperp[0]*ux[0]*pow(bx[1],2)+0.6123724356957944*ppar[0]*ux[0]*pow(bx[1],2)-0.6123724356957944*pperp[1]*ux[1]*pow(bx[0],2)+0.6123724356957944*ppar[1]*ux[1]*pow(bx[0],2)-0.6123724356957944*pperp[0]*ux[0]*pow(bx[0],2)+0.6123724356957944*ppar[0]*ux[0]*pow(bx[0],2)-1.1022703842524295*bx[1]*bz[1]*pperp[1]*uz[1]-0.6123724356957944*bx[0]*bz[0]*pperp[1]*uz[1]+1.1022703842524295*bx[1]*bz[1]*ppar[1]*uz[1]+0.6123724356957944*bx[0]*bz[0]*ppar[1]*uz[1]-0.6123724356957944*bx[0]*pperp[0]*bz[1]*uz[1]+0.6123724356957944*bx[0]*ppar[0]*bz[1]*uz[1]-0.6123724356957944*bz[0]*pperp[0]*bx[1]*uz[1]+0.6123724356957944*bz[0]*ppar[0]*bx[1]*uz[1]-1.1022703842524295*bx[1]*by[1]*pperp[1]*uy[1]-0.6123724356957944*bx[0]*by[0]*pperp[1]*uy[1]+1.1022703842524295*bx[1]*by[1]*ppar[1]*uy[1]+0.6123724356957944*bx[0]*by[0]*ppar[1]*uy[1]-0.6123724356957944*bx[0]*pperp[0]*by[1]*uy[1]+0.6123724356957944*bx[0]*ppar[0]*by[1]*uy[1]-0.6123724356957944*by[0]*pperp[0]*bx[1]*uy[1]+0.6123724356957944*by[0]*ppar[0]*bx[1]*uy[1]+1.224744871391589*pperp[1]*ux[1]-1.224744871391589*bx[0]*pperp[0]*bx[1]*ux[1]+1.224744871391589*bx[0]*ppar[0]*bx[1]*ux[1]+1.224744871391589*E[1]*ux[1]+1.224744871391589*bx[1]*qperp[1]+0.6123724356957944*bx[1]*qpar[1]-0.6123724356957944*bx[0]*uz[0]*bz[1]*pperp[1]-0.6123724356957944*bx[0]*uy[0]*by[1]*pperp[1]-0.6123724356957944*bz[0]*uz[0]*bx[1]*pperp[1]-0.6123724356957944*by[0]*uy[0]*bx[1]*pperp[1]-1.224744871391589*bx[0]*ux[0]*bx[1]*pperp[1]+0.6123724356957944*bx[0]*uz[0]*bz[1]*ppar[1]+0.6123724356957944*bx[0]*uy[0]*by[1]*ppar[1]+0.6123724356957944*bz[0]*uz[0]*bx[1]*ppar[1]+0.6123724356957944*by[0]*uy[0]*bx[1]*ppar[1]+1.224744871391589*bx[0]*ux[0]*bx[1]*ppar[1]-0.6123724356957944*pperp[0]*uz[0]*bx[1]*bz[1]+0.6123724356957944*ppar[0]*uz[0]*bx[1]*bz[1]-0.6123724356957944*pperp[0]*uy[0]*bx[1]*by[1]+0.6123724356957944*ppar[0]*uy[0]*bx[1]*by[1]-0.6123724356957944*bx[0]*bz[0]*pperp[0]*uz[0]+0.6123724356957944*bx[0]*bz[0]*ppar[0]*uz[0]-0.6123724356957944*bx[0]*by[0]*pperp[0]*uy[0]+0.6123724356957944*bx[0]*by[0]*ppar[0]*uy[0]+1.224744871391589*pperp[0]*ux[0]+1.224744871391589*E[0]*ux[0]+1.224744871391589*bx[0]*qperp[0]+0.6123724356957944*bx[0]*qpar[0])*dx10; 
 
  double cflFreq_mid = 0.0; 
  cflFreq_mid += 1.5*dx10*fabs(0.7071067811865475*ux[0]); 
  return cflFreq_mid; 
} 
