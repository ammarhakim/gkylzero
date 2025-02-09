#include <gkyl_canonical_pb_kernels.h> 
void canonical_pb_fluid_hasegawa_wakatani_source_2x_ser_p2(const double *dxv, double alpha, double kappa, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  // dxv[NDIM]: cell length.
  // alpha: Adiabaticity parameter for adiabatic coupling of vorticity and density (zero for Hasegawa-Mima).
  // kappa: Constant density gradient scale length.
  // phi: Potential in fluid system given by grad^2 phi = f where f is (one of) the evolved quantities.
  // f: input state vector in center cell.
  // out: output increment in center cell.
  double dy = 2.0/dxv[1]; 
  const double *n = &f[8]; 
  double *out_zeta = &out[0]; 
  double *out_n = &out[8]; 
  out_zeta[0] += (phi[0]-1.0*n[0])*alpha; 
  out_zeta[1] += (phi[1]-1.0*n[1])*alpha; 
  out_zeta[2] += (phi[2]-1.0*n[2])*alpha; 
  out_zeta[3] += (phi[3]-1.0*n[3])*alpha; 
  out_zeta[4] += (phi[4]-1.0*n[4])*alpha; 
  out_zeta[5] += (phi[5]-1.0*n[5])*alpha; 
  out_zeta[6] += (phi[6]-1.0*n[6])*alpha; 
  out_zeta[7] += (phi[7]-1.0*n[7])*alpha; 
  out_n[0] += (phi[0]-1.0*n[0])*alpha-1.7320508075688772*phi[2]*dy*kappa; 
  out_n[1] += (phi[1]-1.0*n[1])*alpha-1.7320508075688772*phi[3]*dy*kappa; 
  out_n[2] += (phi[2]-1.0*n[2])*alpha-3.872983346207417*phi[5]*dy*kappa; 
  out_n[3] += (phi[3]-1.0*n[3])*alpha-3.872983346207417*phi[7]*dy*kappa; 
  out_n[4] += -(0.2*(8.660254037844387*phi[6]*dy*kappa+(5.0*n[4]-5.0*phi[4])*alpha)); 
  out_n[5] += (phi[5]-1.0*n[5])*alpha; 
  out_n[6] += (phi[6]-1.0*n[6])*alpha; 
  out_n[7] += (phi[7]-1.0*n[7])*alpha; 
} 
