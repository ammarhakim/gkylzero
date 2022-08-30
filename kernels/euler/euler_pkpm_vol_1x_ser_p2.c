#include <gkyl_euler_kernels.h> 
GKYL_CU_DH double euler_pkpm_vol_1x_ser_p2(const double *w, const double *dxv, const double *u_i, const double *p_ij, const double *vlasov_pkpm_moms, const double *statevec, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // uvar: [ux, uy, uz], Fluid flow.
  // ppar: Fluid parallel pressure (computed from kinetic equation).
  // pperp: Fluid perpendicular pressure.
  // qpar: Fluid parallel heat flux (computed from kinetic equation).
  // statevec: [rho ux, rho uy, rho uz, energy], Fluid input state vector.
  // out: Incremented output.

  double dx10 = 2./dxv[0]; 

  const double *rhou0 = &statevec[0]; 
  const double *rhou1 = &statevec[3]; 
  const double *rhou2 = &statevec[6]; 
  const double *energy = &statevec[9]; 

  const double *ux = &u_i[0]; 
  const double *uy = &u_i[3]; 
  const double *uz = &u_i[6]; 

  const double *Pxx = &p_ij[0]; 
  const double *Pxy = &p_ij[3]; 
  const double *Pxz = &p_ij[6]; 
  const double *Pyy = &p_ij[9]; 
  const double *Pyz = &p_ij[12]; 
  const double *Pzz = &p_ij[15]; 

  const double *qx = &vlasov_pkpm_moms[6]; 

  double *outrhou0 = &out[0]; 
  double *outrhou1 = &out[3]; 
  double *outrhou2 = &out[6]; 
  double *outenergy = &out[9]; 

  double cflFreq_mid = 0.0; 
  cflFreq_mid += 0.5*5.0*dx10*(fabs(0.7071067811865475*ux[0]-0.7905694150420947*ux[2])); 

  outrhou0[1] += 1.224744871391589*rhou0[2]*ux[2]*dx10+1.224744871391589*rhou0[1]*ux[1]*dx10+1.224744871391589*rhou0[0]*ux[0]*dx10+1.732050807568877*Pxx[0]*dx10; 
  outrhou0[2] += 2.449489742783178*rhou0[1]*ux[2]*dx10+2.449489742783178*ux[1]*rhou0[2]*dx10+2.738612787525831*rhou0[0]*ux[1]*dx10+2.738612787525831*ux[0]*rhou0[1]*dx10+3.872983346207417*Pxx[1]*dx10; 

  outrhou1[1] += 1.224744871391589*rhou1[2]*ux[2]*dx10+1.224744871391589*rhou1[1]*ux[1]*dx10+1.224744871391589*rhou1[0]*ux[0]*dx10+1.732050807568877*Pxy[0]*dx10; 
  outrhou1[2] += 2.449489742783178*rhou1[1]*ux[2]*dx10+2.449489742783178*ux[1]*rhou1[2]*dx10+2.738612787525831*rhou1[0]*ux[1]*dx10+2.738612787525831*ux[0]*rhou1[1]*dx10+3.872983346207417*Pxy[1]*dx10; 

  outrhou2[1] += 1.224744871391589*rhou2[2]*ux[2]*dx10+1.224744871391589*rhou2[1]*ux[1]*dx10+1.224744871391589*rhou2[0]*ux[0]*dx10+1.732050807568877*Pxz[0]*dx10; 
  outrhou2[2] += 2.449489742783178*rhou2[1]*ux[2]*dx10+2.449489742783178*ux[1]*rhou2[2]*dx10+2.738612787525831*rhou2[0]*ux[1]*dx10+2.738612787525831*ux[0]*rhou2[1]*dx10+3.872983346207417*Pxz[1]*dx10; 

  outenergy[1] += 1.224744871391589*Pxz[2]*uz[2]*dx10+1.224744871391589*Pxy[2]*uy[2]*dx10+1.224744871391589*energy[2]*ux[2]*dx10+1.224744871391589*Pxx[2]*ux[2]*dx10+1.224744871391589*Pxz[1]*uz[1]*dx10+1.224744871391589*Pxy[1]*uy[1]*dx10+1.224744871391589*energy[1]*ux[1]*dx10+1.224744871391589*Pxx[1]*ux[1]*dx10+1.224744871391589*Pxz[0]*uz[0]*dx10+1.224744871391589*Pxy[0]*uy[0]*dx10+1.224744871391589*energy[0]*ux[0]*dx10+1.224744871391589*Pxx[0]*ux[0]*dx10+1.732050807568877*qx[0]*dx10; 
  outenergy[2] += 2.449489742783178*Pxz[1]*uz[2]*dx10+2.449489742783178*Pxy[1]*uy[2]*dx10+2.449489742783178*energy[1]*ux[2]*dx10+2.449489742783178*Pxx[1]*ux[2]*dx10+2.449489742783178*ux[1]*energy[2]*dx10+2.449489742783178*uz[1]*Pxz[2]*dx10+2.449489742783178*uy[1]*Pxy[2]*dx10+2.449489742783178*ux[1]*Pxx[2]*dx10+2.738612787525831*Pxz[0]*uz[1]*dx10+2.738612787525831*Pxy[0]*uy[1]*dx10+2.738612787525831*energy[0]*ux[1]*dx10+2.738612787525831*Pxx[0]*ux[1]*dx10+3.872983346207417*qx[1]*dx10+2.738612787525831*ux[0]*energy[1]*dx10+2.738612787525831*uz[0]*Pxz[1]*dx10+2.738612787525831*uy[0]*Pxy[1]*dx10+2.738612787525831*ux[0]*Pxx[1]*dx10; 

  return cflFreq_mid; 
} 
