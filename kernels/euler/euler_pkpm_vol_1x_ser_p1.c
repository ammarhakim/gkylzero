#include <gkyl_euler_kernels.h> 
GKYL_CU_DH double euler_pkpm_vol_1x_ser_p1(const double *w, const double *dxv, const double *u_i, const double *p_ij, const double *vlasov_pkpm_moms, const double *statevec, double* GKYL_RESTRICT out) 
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
  const double *rhou1 = &statevec[2]; 
  const double *rhou2 = &statevec[4]; 
  const double *energy = &statevec[6]; 

  const double *ux = &u_i[0]; 
  const double *uy = &u_i[2]; 
  const double *uz = &u_i[4]; 

  const double *Pxx = &p_ij[0]; 
  const double *Pxy = &p_ij[2]; 
  const double *Pxz = &p_ij[4]; 
  const double *Pyy = &p_ij[6]; 
  const double *Pyz = &p_ij[8]; 
  const double *Pzz = &p_ij[10]; 

  const double *qx = &vlasov_pkpm_moms[4]; 

  double *outrhou0 = &out[0]; 
  double *outrhou1 = &out[2]; 
  double *outrhou2 = &out[4]; 
  double *outenergy = &out[6]; 

  double cflFreq_mid = 0.0; 
  cflFreq_mid += 0.5*3.0*dx10*(fabs(0.7071067811865475*ux[0])); 

  outrhou0[1] += 1.224744871391589*rhou0[1]*ux[1]*dx10+1.224744871391589*rhou0[0]*ux[0]*dx10+1.732050807568877*Pxx[0]*dx10; 

  outrhou1[1] += 1.224744871391589*rhou1[1]*ux[1]*dx10+1.224744871391589*rhou1[0]*ux[0]*dx10+1.732050807568877*Pxy[0]*dx10; 

  outrhou2[1] += 1.224744871391589*rhou2[1]*ux[1]*dx10+1.224744871391589*rhou2[0]*ux[0]*dx10+1.732050807568877*Pxz[0]*dx10; 

  outenergy[1] += 1.224744871391589*Pxz[1]*uz[1]*dx10+1.224744871391589*Pxy[1]*uy[1]*dx10+1.224744871391589*energy[1]*ux[1]*dx10+1.224744871391589*Pxx[1]*ux[1]*dx10+1.224744871391589*Pxz[0]*uz[0]*dx10+1.224744871391589*Pxy[0]*uy[0]*dx10+1.224744871391589*energy[0]*ux[0]*dx10+1.224744871391589*Pxx[0]*ux[0]*dx10+1.732050807568877*qx[0]*dx10; 

  return cflFreq_mid; 
} 
