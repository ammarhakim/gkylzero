#include <gkyl_euler_kernels.h> 
GKYL_CU_DH double euler_pkpm_vol_1x_ser_p1(const double *w, const double *dxv, const double *u_i, const double *p_ij, const double *statevec, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // u_i:       bulk flow velocity (ux, uy, uz).
  // p_ij:      pressure tensor (P_xx, P_xy, P_xz, P_yy, P_yz, P_zz).
  // statevec: [rho ux, rho uy, rho uz, p_perp], Fluid input state vector.
  // out: Incremented output.

  double dx10 = 2./dxv[0]; 

  const double *rhoux = &statevec[0]; 
  const double *rhouy = &statevec[2]; 
  const double *rhouz = &statevec[4]; 
  const double *p_perp = &statevec[6]; 

  const double *ux = &u_i[0]; 
  const double *uy = &u_i[2]; 
  const double *uz = &u_i[4]; 

  const double *Pxx = &p_ij[0]; 
  const double *Pxy = &p_ij[2]; 
  const double *Pxz = &p_ij[4]; 
  const double *Pyy = &p_ij[6]; 
  const double *Pyz = &p_ij[8]; 
  const double *Pzz = &p_ij[10]; 

  double *outrhoux = &out[0]; 
  double *outrhouy = &out[2]; 
  double *outrhouz = &out[4]; 
  double *outp_perp = &out[6]; 

  double cflFreq_mid = 0.0; 
  cflFreq_mid += 0.5*3.0*dx10*(fabs(0.7071067811865475*ux[0])); 

  outrhoux[1] += 1.224744871391589*rhoux[1]*ux[1]*dx10+1.224744871391589*rhoux[0]*ux[0]*dx10+1.732050807568877*Pxx[0]*dx10; 

  outrhouy[1] += 1.224744871391589*rhouy[1]*ux[1]*dx10+1.224744871391589*rhouy[0]*ux[0]*dx10+1.732050807568877*Pxy[0]*dx10; 

  outrhouz[1] += 1.224744871391589*rhouz[1]*ux[1]*dx10+1.224744871391589*rhouz[0]*ux[0]*dx10+1.732050807568877*Pxz[0]*dx10; 

  outp_perp[1] += 1.224744871391589*p_perp[1]*ux[1]*dx10+1.224744871391589*p_perp[0]*ux[0]*dx10; 

  return cflFreq_mid; 
} 
