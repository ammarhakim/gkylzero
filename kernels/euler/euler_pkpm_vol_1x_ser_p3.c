#include <gkyl_euler_kernels.h> 
GKYL_CU_DH double euler_pkpm_vol_1x_ser_p3(const double *w, const double *dxv, const double *u_i, const double *p_ij, const double *statevec, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // u_i:       bulk flow velocity (ux, uy, uz).
  // p_ij:      pressure tensor (P_xx, P_xy, P_xz, P_yy, P_yz, P_zz).
  // statevec: [rho ux, rho uy, rho uz, p_perp], Fluid input state vector.
  // out: Incremented output.

  double dx10 = 2./dxv[0]; 

  const double *rhoux = &statevec[0]; 
  const double *rhouy = &statevec[4]; 
  const double *rhouz = &statevec[8]; 
  const double *p_perp = &statevec[12]; 

  const double *ux = &u_i[0]; 
  const double *uy = &u_i[4]; 
  const double *uz = &u_i[8]; 

  const double *Pxx = &p_ij[0]; 
  const double *Pxy = &p_ij[4]; 
  const double *Pxz = &p_ij[8]; 
  const double *Pyy = &p_ij[12]; 
  const double *Pyz = &p_ij[16]; 
  const double *Pzz = &p_ij[20]; 

  double *outrhoux = &out[0]; 
  double *outrhouy = &out[4]; 
  double *outrhouz = &out[8]; 
  double *outp_perp = &out[12]; 

  double cflFreq_mid = 0.0; 
  cflFreq_mid += 0.5*7.0*dx10*(fabs(0.7071067811865475*ux[0]-0.7905694150420947*ux[2])); 

  outrhoux[1] += 1.224744871391589*rhoux[3]*ux[3]*dx10+1.224744871391589*rhoux[2]*ux[2]*dx10+1.224744871391589*rhoux[1]*ux[1]*dx10+1.224744871391589*rhoux[0]*ux[0]*dx10+1.732050807568877*Pxx[0]*dx10; 
  outrhoux[2] += 2.405351177211819*rhoux[2]*ux[3]*dx10+2.405351177211819*ux[2]*rhoux[3]*dx10+2.449489742783178*rhoux[1]*ux[2]*dx10+2.449489742783178*ux[1]*rhoux[2]*dx10+2.738612787525831*rhoux[0]*ux[1]*dx10+2.738612787525831*ux[0]*rhoux[1]*dx10+3.872983346207417*Pxx[1]*dx10; 
  outrhoux[3] += 4.365266951236265*rhoux[3]*ux[3]*dx10+3.674234614174766*rhoux[1]*ux[3]*dx10+3.674234614174766*ux[1]*rhoux[3]*dx10+4.543441112511214*rhoux[2]*ux[2]*dx10+4.183300132670378*rhoux[0]*ux[2]*dx10+4.183300132670378*ux[0]*rhoux[2]*dx10+5.916079783099617*Pxx[2]*dx10+5.612486080160912*rhoux[1]*ux[1]*dx10+1.870828693386971*rhoux[0]*ux[0]*dx10+2.645751311064591*Pxx[0]*dx10; 

  outrhouy[1] += 1.224744871391589*rhouy[3]*ux[3]*dx10+1.224744871391589*rhouy[2]*ux[2]*dx10+1.224744871391589*rhouy[1]*ux[1]*dx10+1.224744871391589*rhouy[0]*ux[0]*dx10+1.732050807568877*Pxy[0]*dx10; 
  outrhouy[2] += 2.405351177211819*rhouy[2]*ux[3]*dx10+2.405351177211819*ux[2]*rhouy[3]*dx10+2.449489742783178*rhouy[1]*ux[2]*dx10+2.449489742783178*ux[1]*rhouy[2]*dx10+2.738612787525831*rhouy[0]*ux[1]*dx10+2.738612787525831*ux[0]*rhouy[1]*dx10+3.872983346207417*Pxy[1]*dx10; 
  outrhouy[3] += 4.365266951236265*rhouy[3]*ux[3]*dx10+3.674234614174766*rhouy[1]*ux[3]*dx10+3.674234614174766*ux[1]*rhouy[3]*dx10+4.543441112511214*rhouy[2]*ux[2]*dx10+4.183300132670378*rhouy[0]*ux[2]*dx10+4.183300132670378*ux[0]*rhouy[2]*dx10+5.916079783099617*Pxy[2]*dx10+5.612486080160912*rhouy[1]*ux[1]*dx10+1.870828693386971*rhouy[0]*ux[0]*dx10+2.645751311064591*Pxy[0]*dx10; 

  outrhouz[1] += 1.224744871391589*rhouz[3]*ux[3]*dx10+1.224744871391589*rhouz[2]*ux[2]*dx10+1.224744871391589*rhouz[1]*ux[1]*dx10+1.224744871391589*rhouz[0]*ux[0]*dx10+1.732050807568877*Pxz[0]*dx10; 
  outrhouz[2] += 2.405351177211819*rhouz[2]*ux[3]*dx10+2.405351177211819*ux[2]*rhouz[3]*dx10+2.449489742783178*rhouz[1]*ux[2]*dx10+2.449489742783178*ux[1]*rhouz[2]*dx10+2.738612787525831*rhouz[0]*ux[1]*dx10+2.738612787525831*ux[0]*rhouz[1]*dx10+3.872983346207417*Pxz[1]*dx10; 
  outrhouz[3] += 4.365266951236265*rhouz[3]*ux[3]*dx10+3.674234614174766*rhouz[1]*ux[3]*dx10+3.674234614174766*ux[1]*rhouz[3]*dx10+4.543441112511214*rhouz[2]*ux[2]*dx10+4.183300132670378*rhouz[0]*ux[2]*dx10+4.183300132670378*ux[0]*rhouz[2]*dx10+5.916079783099617*Pxz[2]*dx10+5.612486080160912*rhouz[1]*ux[1]*dx10+1.870828693386971*rhouz[0]*ux[0]*dx10+2.645751311064591*Pxz[0]*dx10; 

  outp_perp[1] += 1.224744871391589*p_perp[3]*ux[3]*dx10+1.224744871391589*p_perp[2]*ux[2]*dx10+1.224744871391589*p_perp[1]*ux[1]*dx10+1.224744871391589*p_perp[0]*ux[0]*dx10; 
  outp_perp[2] += 2.405351177211819*p_perp[2]*ux[3]*dx10+2.405351177211819*ux[2]*p_perp[3]*dx10+2.449489742783178*p_perp[1]*ux[2]*dx10+2.449489742783178*ux[1]*p_perp[2]*dx10+2.738612787525831*p_perp[0]*ux[1]*dx10+2.738612787525831*ux[0]*p_perp[1]*dx10; 
  outp_perp[3] += 4.365266951236265*p_perp[3]*ux[3]*dx10+3.674234614174766*p_perp[1]*ux[3]*dx10+3.674234614174766*ux[1]*p_perp[3]*dx10+4.543441112511214*p_perp[2]*ux[2]*dx10+4.183300132670378*p_perp[0]*ux[2]*dx10+4.183300132670378*ux[0]*p_perp[2]*dx10+5.612486080160912*p_perp[1]*ux[1]*dx10+1.870828693386971*p_perp[0]*ux[0]*dx10; 

  return cflFreq_mid; 
} 
