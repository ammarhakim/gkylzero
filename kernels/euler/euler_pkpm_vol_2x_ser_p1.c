#include <gkyl_euler_kernels.h> 
GKYL_CU_DH double euler_pkpm_vol_2x_ser_p1(const double *w, const double *dxv, const double *u_i, const double *p_ij, const double *vlasov_pkpm_moms, const double *statevec, double* GKYL_RESTRICT out) 
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
  double dx11 = 2./dxv[1]; 

  const double *rhou0 = &statevec[0]; 
  const double *rhou1 = &statevec[4]; 
  const double *rhou2 = &statevec[8]; 
  const double *energy = &statevec[12]; 

  const double *ux = &u_i[0]; 
  const double *uy = &u_i[4]; 
  const double *uz = &u_i[8]; 

  const double *Pxx = &p_ij[0]; 
  const double *Pxy = &p_ij[4]; 
  const double *Pxz = &p_ij[8]; 
  const double *Pyy = &p_ij[12]; 
  const double *Pyz = &p_ij[16]; 
  const double *Pzz = &p_ij[20]; 

  const double *qx = &vlasov_pkpm_moms[8]; 
  const double *qy = &vlasov_pkpm_moms[12]; 

  double *outrhou0 = &out[0]; 
  double *outrhou1 = &out[4]; 
  double *outrhou2 = &out[8]; 
  double *outenergy = &out[12]; 

  double alpha_mid = 0.0; 
  alpha_mid += 0.5*dx10*(fabs(0.5*ux[0])); 
  alpha_mid += 0.5*dx11*(fabs(0.5*uy[0])); 

  outrhou0[1] += 0.8660254037844386*rhou0[3]*ux[3]*dx10+0.8660254037844386*rhou0[2]*ux[2]*dx10+0.8660254037844386*rhou0[1]*ux[1]*dx10+0.8660254037844386*rhou0[0]*ux[0]*dx10+1.732050807568877*Pxx[0]*dx10; 
  outrhou0[2] += 0.8660254037844386*rhou0[3]*uy[3]*dx11+0.8660254037844386*rhou0[2]*uy[2]*dx11+0.8660254037844386*rhou0[1]*uy[1]*dx11+0.8660254037844386*rhou0[0]*uy[0]*dx11+1.732050807568877*Pxy[0]*dx11; 
  outrhou0[3] += 0.8660254037844386*rhou0[2]*uy[3]*dx11+0.8660254037844386*uy[2]*rhou0[3]*dx11+0.8660254037844386*rhou0[0]*uy[1]*dx11+0.8660254037844386*uy[0]*rhou0[1]*dx11+1.732050807568877*Pxy[1]*dx11+0.8660254037844386*rhou0[1]*ux[3]*dx10+0.8660254037844386*ux[1]*rhou0[3]*dx10+0.8660254037844386*rhou0[0]*ux[2]*dx10+0.8660254037844386*ux[0]*rhou0[2]*dx10+1.732050807568877*Pxx[2]*dx10; 

  outrhou1[1] += 0.8660254037844386*rhou1[3]*ux[3]*dx10+0.8660254037844386*rhou1[2]*ux[2]*dx10+0.8660254037844386*rhou1[1]*ux[1]*dx10+0.8660254037844386*rhou1[0]*ux[0]*dx10+1.732050807568877*Pxy[0]*dx10; 
  outrhou1[2] += 0.8660254037844386*rhou1[3]*uy[3]*dx11+0.8660254037844386*rhou1[2]*uy[2]*dx11+0.8660254037844386*rhou1[1]*uy[1]*dx11+0.8660254037844386*rhou1[0]*uy[0]*dx11+1.732050807568877*Pyy[0]*dx11; 
  outrhou1[3] += 0.8660254037844386*rhou1[2]*uy[3]*dx11+0.8660254037844386*uy[2]*rhou1[3]*dx11+0.8660254037844386*rhou1[0]*uy[1]*dx11+0.8660254037844386*uy[0]*rhou1[1]*dx11+1.732050807568877*Pyy[1]*dx11+0.8660254037844386*rhou1[1]*ux[3]*dx10+0.8660254037844386*ux[1]*rhou1[3]*dx10+0.8660254037844386*rhou1[0]*ux[2]*dx10+0.8660254037844386*ux[0]*rhou1[2]*dx10+1.732050807568877*Pxy[2]*dx10; 

  outrhou2[1] += 0.8660254037844386*rhou2[3]*ux[3]*dx10+0.8660254037844386*rhou2[2]*ux[2]*dx10+0.8660254037844386*rhou2[1]*ux[1]*dx10+0.8660254037844386*rhou2[0]*ux[0]*dx10+1.732050807568877*Pxz[0]*dx10; 
  outrhou2[2] += 0.8660254037844386*rhou2[3]*uy[3]*dx11+0.8660254037844386*rhou2[2]*uy[2]*dx11+0.8660254037844386*rhou2[1]*uy[1]*dx11+0.8660254037844386*rhou2[0]*uy[0]*dx11+1.732050807568877*Pyz[0]*dx11; 
  outrhou2[3] += 0.8660254037844386*rhou2[2]*uy[3]*dx11+0.8660254037844386*uy[2]*rhou2[3]*dx11+0.8660254037844386*rhou2[0]*uy[1]*dx11+0.8660254037844386*uy[0]*rhou2[1]*dx11+1.732050807568877*Pyz[1]*dx11+0.8660254037844386*rhou2[1]*ux[3]*dx10+0.8660254037844386*ux[1]*rhou2[3]*dx10+0.8660254037844386*rhou2[0]*ux[2]*dx10+0.8660254037844386*ux[0]*rhou2[2]*dx10+1.732050807568877*Pxz[2]*dx10; 

  outenergy[1] += 0.8660254037844386*Pxz[3]*uz[3]*dx10+0.8660254037844386*Pxy[3]*uy[3]*dx10+0.8660254037844386*energy[3]*ux[3]*dx10+0.8660254037844386*Pxx[3]*ux[3]*dx10+0.8660254037844386*Pxz[2]*uz[2]*dx10+0.8660254037844386*Pxy[2]*uy[2]*dx10+0.8660254037844386*energy[2]*ux[2]*dx10+0.8660254037844386*Pxx[2]*ux[2]*dx10+0.8660254037844386*Pxz[1]*uz[1]*dx10+0.8660254037844386*Pxy[1]*uy[1]*dx10+0.8660254037844386*energy[1]*ux[1]*dx10+0.8660254037844386*Pxx[1]*ux[1]*dx10+0.8660254037844386*Pxz[0]*uz[0]*dx10+0.8660254037844386*Pxy[0]*uy[0]*dx10+0.8660254037844386*energy[0]*ux[0]*dx10+0.8660254037844386*Pxx[0]*ux[0]*dx10+1.732050807568877*qx[0]*dx10; 
  outenergy[2] += 0.8660254037844386*Pyz[3]*uz[3]*dx11+0.8660254037844386*energy[3]*uy[3]*dx11+0.8660254037844386*Pyy[3]*uy[3]*dx11+0.8660254037844386*Pxy[3]*ux[3]*dx11+0.8660254037844386*Pyz[2]*uz[2]*dx11+0.8660254037844386*energy[2]*uy[2]*dx11+0.8660254037844386*Pyy[2]*uy[2]*dx11+0.8660254037844386*Pxy[2]*ux[2]*dx11+0.8660254037844386*Pyz[1]*uz[1]*dx11+0.8660254037844386*energy[1]*uy[1]*dx11+0.8660254037844386*Pyy[1]*uy[1]*dx11+0.8660254037844386*Pxy[1]*ux[1]*dx11+0.8660254037844386*Pyz[0]*uz[0]*dx11+0.8660254037844386*energy[0]*uy[0]*dx11+0.8660254037844386*Pyy[0]*uy[0]*dx11+0.8660254037844386*Pxy[0]*ux[0]*dx11+1.732050807568877*qy[0]*dx11; 
  outenergy[3] += 0.8660254037844386*Pyz[2]*uz[3]*dx11+0.8660254037844386*energy[2]*uy[3]*dx11+0.8660254037844386*Pyy[2]*uy[3]*dx11+0.8660254037844386*Pxy[2]*ux[3]*dx11+0.8660254037844386*uy[2]*energy[3]*dx11+0.8660254037844386*uz[2]*Pyz[3]*dx11+0.8660254037844386*uy[2]*Pyy[3]*dx11+0.8660254037844386*ux[2]*Pxy[3]*dx11+0.8660254037844386*Pyz[0]*uz[1]*dx11+0.8660254037844386*energy[0]*uy[1]*dx11+0.8660254037844386*Pyy[0]*uy[1]*dx11+0.8660254037844386*Pxy[0]*ux[1]*dx11+1.732050807568877*qy[1]*dx11+0.8660254037844386*uy[0]*energy[1]*dx11+0.8660254037844386*uz[0]*Pyz[1]*dx11+0.8660254037844386*uy[0]*Pyy[1]*dx11+0.8660254037844386*ux[0]*Pxy[1]*dx11+0.8660254037844386*Pxz[1]*uz[3]*dx10+0.8660254037844386*Pxy[1]*uy[3]*dx10+0.8660254037844386*energy[1]*ux[3]*dx10+0.8660254037844386*Pxx[1]*ux[3]*dx10+0.8660254037844386*ux[1]*energy[3]*dx10+0.8660254037844386*uz[1]*Pxz[3]*dx10+0.8660254037844386*uy[1]*Pxy[3]*dx10+0.8660254037844386*ux[1]*Pxx[3]*dx10+0.8660254037844386*Pxz[0]*uz[2]*dx10+0.8660254037844386*Pxy[0]*uy[2]*dx10+0.8660254037844386*energy[0]*ux[2]*dx10+0.8660254037844386*Pxx[0]*ux[2]*dx10+1.732050807568877*qx[2]*dx10+0.8660254037844386*ux[0]*energy[2]*dx10+0.8660254037844386*uz[0]*Pxz[2]*dx10+0.8660254037844386*uy[0]*Pxy[2]*dx10+0.8660254037844386*ux[0]*Pxx[2]*dx10; 

  return alpha_mid; 
} 
