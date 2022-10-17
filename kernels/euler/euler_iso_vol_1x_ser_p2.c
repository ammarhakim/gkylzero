#include <gkyl_euler_kernels.h> 
GKYL_CU_DH double euler_iso_vol_1x_ser_p2(const double *w, const double *dxv, const double vth, const double *uvar, const double *statevec, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // vth: Thermal velocity.
  // uvar: [ux, uy, uz], Fluid flow.
  // statevec: [rho, rho ux, rho uy, rho uz].
  // out: Incremented output.

  const double *rho = &statevec[0]; 
  const double *rhou0 = &statevec[3]; 
  const double *rhou1 = &statevec[6]; 
  const double *rhou2 = &statevec[9]; 
  const double *uvar0 = &uvar[0]; 
  const double *uvar1 = &uvar[3]; 
  const double *uvar2 = &uvar[6]; 
  double *outrho = &out[0]; 
  double *outrhou0 = &out[3]; 
  double *outrhou1 = &out[6]; 
  double *outrhou2 = &out[9]; 
  double dx10 = 2./dxv[0]; 

  double cflFreq_mid = 0.0; 
  double vthsq = vth*vth; 
  cflFreq_mid += 0.5*5.0*dx10*(fabs(0.7071067811865475*uvar0[0]-0.7905694150420947*uvar0[2])+vth); 

  outrho[1] += 1.224744871391589*rho[2]*rhou0[2]*dx10+1.224744871391589*rho[1]*rhou0[1]*dx10+1.224744871391589*rho[0]*rhou0[0]*dx10; 
  outrho[2] += 2.449489742783178*rho[1]*rhou0[2]*dx10+2.449489742783178*rhou0[1]*rho[2]*dx10+2.738612787525831*rho[0]*rhou0[1]*dx10+2.738612787525831*rhou0[0]*rho[1]*dx10; 

  outrhou0[1] += 1.732050807568877*rho[0]*dx10*vthsq+1.224744871391589*rhou0[2]*uvar0[2]*dx10+1.224744871391589*rhou0[1]*uvar0[1]*dx10+1.224744871391589*rhou0[0]*uvar0[0]*dx10; 
  outrhou0[2] += 3.872983346207417*rho[1]*dx10*vthsq+2.449489742783178*rhou0[1]*uvar0[2]*dx10+2.449489742783178*uvar0[1]*rhou0[2]*dx10+2.738612787525831*rhou0[0]*uvar0[1]*dx10+2.738612787525831*uvar0[0]*rhou0[1]*dx10; 

  outrhou1[1] += 1.224744871391589*rhou1[2]*uvar0[2]*dx10+1.224744871391589*rhou1[1]*uvar0[1]*dx10+1.224744871391589*rhou1[0]*uvar0[0]*dx10; 
  outrhou1[2] += 2.449489742783178*rhou1[1]*uvar0[2]*dx10+2.449489742783178*uvar0[1]*rhou1[2]*dx10+2.738612787525831*rhou1[0]*uvar0[1]*dx10+2.738612787525831*uvar0[0]*rhou1[1]*dx10; 

  outrhou2[1] += 1.224744871391589*rhou2[2]*uvar0[2]*dx10+1.224744871391589*rhou2[1]*uvar0[1]*dx10+1.224744871391589*rhou2[0]*uvar0[0]*dx10; 
  outrhou2[2] += 2.449489742783178*rhou2[1]*uvar0[2]*dx10+2.449489742783178*uvar0[1]*rhou2[2]*dx10+2.738612787525831*rhou2[0]*uvar0[1]*dx10+2.738612787525831*uvar0[0]*rhou2[1]*dx10; 

  return cflFreq_mid; 
} 
