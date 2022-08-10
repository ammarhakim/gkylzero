#include <gkyl_isoeuler_kernels.h> 
GKYL_CU_DH double isoeuler_vol_2x_ser_p2(const double *w, const double *dxv, const double vth, const double *uvar, const double *statevec, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // statevec: [rho, rho ux, rho uy, rho uz].
  // uvar: [ux, uy, uz].
  // out: Incremented output.

  const double *rho = &statevec[0]; 
  const double *rhou0 = &statevec[8]; 
  const double *rhou1 = &statevec[16]; 
  const double *rhou2 = &statevec[24]; 
  const double *uvar0 = &uvar[0]; 
  const double *uvar1 = &uvar[8]; 
  const double *uvar2 = &uvar[16]; 
  double *outrho = &out[0]; 
  double *outrhou0 = &out[8]; 
  double *outrhou1 = &out[16]; 
  double *outrhou2 = &out[24]; 
  double dx10 = 2./dxv[0]; 
  double dx11 = 2./dxv[1]; 

  double vthsq = vth*vth; 
  outrho[0] += rho[7]*rhou1[7]*dx11+rho[6]*rhou1[6]*dx11+rho[5]*rhou1[5]*dx11+rho[4]*rhou1[4]*dx11+rho[3]*rhou1[3]*dx11+rho[2]*rhou1[2]*dx11+rho[1]*rhou1[1]*dx11+rho[0]*rhou1[0]*dx11+rho[7]*rhou0[7]*dx10+rho[6]*rhou0[6]*dx10+rho[5]*rhou0[5]*dx10+rho[4]*rhou0[4]*dx10+rho[3]*rhou0[3]*dx10+rho[2]*rhou0[2]*dx10+rho[1]*rhou0[1]*dx10+rho[0]*rhou0[0]*dx10; 

  outrhou0[0] += 1.732050807568877*rho[1]*vthsq+rhou0[7]*uvar1[7]*dx11+rhou0[6]*uvar1[6]*dx11+rhou0[5]*uvar1[5]*dx11+rhou0[4]*uvar1[4]*dx11+rhou0[3]*uvar1[3]*dx11+rhou0[2]*uvar1[2]*dx11+rhou0[1]*uvar1[1]*dx11+rhou0[0]*uvar1[0]*dx11+rhou0[7]*uvar0[7]*dx10+rhou0[6]*uvar0[6]*dx10+rhou0[5]*uvar0[5]*dx10+rhou0[4]*uvar0[4]*dx10+rhou0[3]*uvar0[3]*dx10+rhou0[2]*uvar0[2]*dx10+rhou0[1]*uvar0[1]*dx10+rhou0[0]*uvar0[0]*dx10; 
  outrhou0[1] += 3.872983346207417*rho[4]*vthsq+rhou0[7]*uvar1[7]*dx11+rhou0[6]*uvar1[6]*dx11+rhou0[5]*uvar1[5]*dx11+rhou0[4]*uvar1[4]*dx11+rhou0[3]*uvar1[3]*dx11+rhou0[2]*uvar1[2]*dx11+rhou0[1]*uvar1[1]*dx11+rhou0[0]*uvar1[0]*dx11+rhou0[7]*uvar0[7]*dx10+rhou0[6]*uvar0[6]*dx10+rhou0[5]*uvar0[5]*dx10+rhou0[4]*uvar0[4]*dx10+rhou0[3]*uvar0[3]*dx10+rhou0[2]*uvar0[2]*dx10+rhou0[1]*uvar0[1]*dx10+rhou0[0]*uvar0[0]*dx10; 
  outrhou0[2] += 1.732050807568877*rho[3]*vthsq+rhou0[7]*uvar1[7]*dx11+rhou0[6]*uvar1[6]*dx11+rhou0[5]*uvar1[5]*dx11+rhou0[4]*uvar1[4]*dx11+rhou0[3]*uvar1[3]*dx11+rhou0[2]*uvar1[2]*dx11+rhou0[1]*uvar1[1]*dx11+rhou0[0]*uvar1[0]*dx11+rhou0[7]*uvar0[7]*dx10+rhou0[6]*uvar0[6]*dx10+rhou0[5]*uvar0[5]*dx10+rhou0[4]*uvar0[4]*dx10+rhou0[3]*uvar0[3]*dx10+rhou0[2]*uvar0[2]*dx10+rhou0[1]*uvar0[1]*dx10+rhou0[0]*uvar0[0]*dx10; 
  outrhou0[3] += 3.872983346207417*rho[6]*vthsq+rhou0[7]*uvar1[7]*dx11+rhou0[6]*uvar1[6]*dx11+rhou0[5]*uvar1[5]*dx11+rhou0[4]*uvar1[4]*dx11+rhou0[3]*uvar1[3]*dx11+rhou0[2]*uvar1[2]*dx11+rhou0[1]*uvar1[1]*dx11+rhou0[0]*uvar1[0]*dx11+rhou0[7]*uvar0[7]*dx10+rhou0[6]*uvar0[6]*dx10+rhou0[5]*uvar0[5]*dx10+rhou0[4]*uvar0[4]*dx10+rhou0[3]*uvar0[3]*dx10+rhou0[2]*uvar0[2]*dx10+rhou0[1]*uvar0[1]*dx10+rhou0[0]*uvar0[0]*dx10; 
  outrhou0[4] += rhou0[7]*uvar1[7]*dx11+rhou0[6]*uvar1[6]*dx11+rhou0[5]*uvar1[5]*dx11+rhou0[4]*uvar1[4]*dx11+rhou0[3]*uvar1[3]*dx11+rhou0[2]*uvar1[2]*dx11+rhou0[1]*uvar1[1]*dx11+rhou0[0]*uvar1[0]*dx11+rhou0[7]*uvar0[7]*dx10+rhou0[6]*uvar0[6]*dx10+rhou0[5]*uvar0[5]*dx10+rhou0[4]*uvar0[4]*dx10+rhou0[3]*uvar0[3]*dx10+rhou0[2]*uvar0[2]*dx10+rhou0[1]*uvar0[1]*dx10+rhou0[0]*uvar0[0]*dx10; 
  outrhou0[5] += 1.732050807568877*rho[7]*vthsq+rhou0[7]*uvar1[7]*dx11+rhou0[6]*uvar1[6]*dx11+rhou0[5]*uvar1[5]*dx11+rhou0[4]*uvar1[4]*dx11+rhou0[3]*uvar1[3]*dx11+rhou0[2]*uvar1[2]*dx11+rhou0[1]*uvar1[1]*dx11+rhou0[0]*uvar1[0]*dx11+rhou0[7]*uvar0[7]*dx10+rhou0[6]*uvar0[6]*dx10+rhou0[5]*uvar0[5]*dx10+rhou0[4]*uvar0[4]*dx10+rhou0[3]*uvar0[3]*dx10+rhou0[2]*uvar0[2]*dx10+rhou0[1]*uvar0[1]*dx10+rhou0[0]*uvar0[0]*dx10; 
  outrhou0[6] += rhou0[7]*uvar1[7]*dx11+rhou0[6]*uvar1[6]*dx11+rhou0[5]*uvar1[5]*dx11+rhou0[4]*uvar1[4]*dx11+rhou0[3]*uvar1[3]*dx11+rhou0[2]*uvar1[2]*dx11+rhou0[1]*uvar1[1]*dx11+rhou0[0]*uvar1[0]*dx11+rhou0[7]*uvar0[7]*dx10+rhou0[6]*uvar0[6]*dx10+rhou0[5]*uvar0[5]*dx10+rhou0[4]*uvar0[4]*dx10+rhou0[3]*uvar0[3]*dx10+rhou0[2]*uvar0[2]*dx10+rhou0[1]*uvar0[1]*dx10+rhou0[0]*uvar0[0]*dx10; 
  outrhou0[7] += rhou0[7]*uvar1[7]*dx11+rhou0[6]*uvar1[6]*dx11+rhou0[5]*uvar1[5]*dx11+rhou0[4]*uvar1[4]*dx11+rhou0[3]*uvar1[3]*dx11+rhou0[2]*uvar1[2]*dx11+rhou0[1]*uvar1[1]*dx11+rhou0[0]*uvar1[0]*dx11+rhou0[7]*uvar0[7]*dx10+rhou0[6]*uvar0[6]*dx10+rhou0[5]*uvar0[5]*dx10+rhou0[4]*uvar0[4]*dx10+rhou0[3]*uvar0[3]*dx10+rhou0[2]*uvar0[2]*dx10+rhou0[1]*uvar0[1]*dx10+rhou0[0]*uvar0[0]*dx10; 

  outrhou1[0] += 1.732050807568877*rho[2]*vthsq+rhou1[7]*uvar1[7]*dx11+rhou1[6]*uvar1[6]*dx11+rhou1[5]*uvar1[5]*dx11+rhou1[4]*uvar1[4]*dx11+rhou1[3]*uvar1[3]*dx11+rhou1[2]*uvar1[2]*dx11+rhou1[1]*uvar1[1]*dx11+rhou1[0]*uvar1[0]*dx11+rhou1[7]*uvar0[7]*dx10+rhou1[6]*uvar0[6]*dx10+rhou1[5]*uvar0[5]*dx10+rhou1[4]*uvar0[4]*dx10+rhou1[3]*uvar0[3]*dx10+rhou1[2]*uvar0[2]*dx10+rhou1[1]*uvar0[1]*dx10+rhou1[0]*uvar0[0]*dx10; 
  outrhou1[1] += 1.732050807568877*rho[3]*vthsq+rhou1[7]*uvar1[7]*dx11+rhou1[6]*uvar1[6]*dx11+rhou1[5]*uvar1[5]*dx11+rhou1[4]*uvar1[4]*dx11+rhou1[3]*uvar1[3]*dx11+rhou1[2]*uvar1[2]*dx11+rhou1[1]*uvar1[1]*dx11+rhou1[0]*uvar1[0]*dx11+rhou1[7]*uvar0[7]*dx10+rhou1[6]*uvar0[6]*dx10+rhou1[5]*uvar0[5]*dx10+rhou1[4]*uvar0[4]*dx10+rhou1[3]*uvar0[3]*dx10+rhou1[2]*uvar0[2]*dx10+rhou1[1]*uvar0[1]*dx10+rhou1[0]*uvar0[0]*dx10; 
  outrhou1[2] += 3.872983346207417*rho[5]*vthsq+rhou1[7]*uvar1[7]*dx11+rhou1[6]*uvar1[6]*dx11+rhou1[5]*uvar1[5]*dx11+rhou1[4]*uvar1[4]*dx11+rhou1[3]*uvar1[3]*dx11+rhou1[2]*uvar1[2]*dx11+rhou1[1]*uvar1[1]*dx11+rhou1[0]*uvar1[0]*dx11+rhou1[7]*uvar0[7]*dx10+rhou1[6]*uvar0[6]*dx10+rhou1[5]*uvar0[5]*dx10+rhou1[4]*uvar0[4]*dx10+rhou1[3]*uvar0[3]*dx10+rhou1[2]*uvar0[2]*dx10+rhou1[1]*uvar0[1]*dx10+rhou1[0]*uvar0[0]*dx10; 
  outrhou1[3] += 3.872983346207417*rho[7]*vthsq+rhou1[7]*uvar1[7]*dx11+rhou1[6]*uvar1[6]*dx11+rhou1[5]*uvar1[5]*dx11+rhou1[4]*uvar1[4]*dx11+rhou1[3]*uvar1[3]*dx11+rhou1[2]*uvar1[2]*dx11+rhou1[1]*uvar1[1]*dx11+rhou1[0]*uvar1[0]*dx11+rhou1[7]*uvar0[7]*dx10+rhou1[6]*uvar0[6]*dx10+rhou1[5]*uvar0[5]*dx10+rhou1[4]*uvar0[4]*dx10+rhou1[3]*uvar0[3]*dx10+rhou1[2]*uvar0[2]*dx10+rhou1[1]*uvar0[1]*dx10+rhou1[0]*uvar0[0]*dx10; 
  outrhou1[4] += 1.732050807568877*rho[6]*vthsq+rhou1[7]*uvar1[7]*dx11+rhou1[6]*uvar1[6]*dx11+rhou1[5]*uvar1[5]*dx11+rhou1[4]*uvar1[4]*dx11+rhou1[3]*uvar1[3]*dx11+rhou1[2]*uvar1[2]*dx11+rhou1[1]*uvar1[1]*dx11+rhou1[0]*uvar1[0]*dx11+rhou1[7]*uvar0[7]*dx10+rhou1[6]*uvar0[6]*dx10+rhou1[5]*uvar0[5]*dx10+rhou1[4]*uvar0[4]*dx10+rhou1[3]*uvar0[3]*dx10+rhou1[2]*uvar0[2]*dx10+rhou1[1]*uvar0[1]*dx10+rhou1[0]*uvar0[0]*dx10; 
  outrhou1[5] += rhou1[7]*uvar1[7]*dx11+rhou1[6]*uvar1[6]*dx11+rhou1[5]*uvar1[5]*dx11+rhou1[4]*uvar1[4]*dx11+rhou1[3]*uvar1[3]*dx11+rhou1[2]*uvar1[2]*dx11+rhou1[1]*uvar1[1]*dx11+rhou1[0]*uvar1[0]*dx11+rhou1[7]*uvar0[7]*dx10+rhou1[6]*uvar0[6]*dx10+rhou1[5]*uvar0[5]*dx10+rhou1[4]*uvar0[4]*dx10+rhou1[3]*uvar0[3]*dx10+rhou1[2]*uvar0[2]*dx10+rhou1[1]*uvar0[1]*dx10+rhou1[0]*uvar0[0]*dx10; 
  outrhou1[6] += rhou1[7]*uvar1[7]*dx11+rhou1[6]*uvar1[6]*dx11+rhou1[5]*uvar1[5]*dx11+rhou1[4]*uvar1[4]*dx11+rhou1[3]*uvar1[3]*dx11+rhou1[2]*uvar1[2]*dx11+rhou1[1]*uvar1[1]*dx11+rhou1[0]*uvar1[0]*dx11+rhou1[7]*uvar0[7]*dx10+rhou1[6]*uvar0[6]*dx10+rhou1[5]*uvar0[5]*dx10+rhou1[4]*uvar0[4]*dx10+rhou1[3]*uvar0[3]*dx10+rhou1[2]*uvar0[2]*dx10+rhou1[1]*uvar0[1]*dx10+rhou1[0]*uvar0[0]*dx10; 
  outrhou1[7] += rhou1[7]*uvar1[7]*dx11+rhou1[6]*uvar1[6]*dx11+rhou1[5]*uvar1[5]*dx11+rhou1[4]*uvar1[4]*dx11+rhou1[3]*uvar1[3]*dx11+rhou1[2]*uvar1[2]*dx11+rhou1[1]*uvar1[1]*dx11+rhou1[0]*uvar1[0]*dx11+rhou1[7]*uvar0[7]*dx10+rhou1[6]*uvar0[6]*dx10+rhou1[5]*uvar0[5]*dx10+rhou1[4]*uvar0[4]*dx10+rhou1[3]*uvar0[3]*dx10+rhou1[2]*uvar0[2]*dx10+rhou1[1]*uvar0[1]*dx10+rhou1[0]*uvar0[0]*dx10; 

  return 0.; 
} 
