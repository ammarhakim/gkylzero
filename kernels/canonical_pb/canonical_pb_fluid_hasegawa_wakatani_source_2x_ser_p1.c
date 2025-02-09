#include <gkyl_canonical_pb_kernels.h> 
void canonical_pb_fluid_hasegawa_wakatani_source_2x_ser_p1(const double *dxv, double alpha, double kappa, const double *background_n_gradient, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  // dxv[NDIM]: cell length.
  // alpha: Adiabaticity parameter for adiabatic coupling of vorticity and density (zero for Hasegawa-Mima).
  // kappa: Constant density gradient scale length.
  // background_n_gradient: Background density gradient.
  // phi: Potential in fluid system given by grad^2 phi = f where f is (one of) the evolved quantities.
  // f: input state vector in center cell.
  // out: output increment in center cell.
  double dx = 2.0/dxv[0]; 
  double dy = 2.0/dxv[1]; 
  const double *n = &f[4]; 
  double *out_zeta = &out[0]; 
  double *out_n = &out[4]; 
  out_zeta[0] += (phi[0]-1.0*n[0])*alpha; 
  out_zeta[1] += (phi[1]-1.0*n[1])*alpha; 
  out_zeta[2] += (phi[2]-1.0*n[2])*alpha; 
  out_zeta[3] += (phi[3]-1.0*n[3])*alpha; 
  out_n[0] += -(0.5*((3.0*background_n_gradient[1]*phi[2]-3.0*phi[1]*background_n_gradient[2])*dx*dy+(2.0*n[0]-2.0*phi[0])*alpha)); 
  out_n[1] += -(0.5*((3.0*background_n_gradient[1]*phi[3]-3.0*phi[1]*background_n_gradient[3])*dx*dy+(2.0*n[1]-2.0*phi[1])*alpha)); 
  out_n[2] += 0.5*((3.0*background_n_gradient[2]*phi[3]-3.0*phi[2]*background_n_gradient[3])*dx*dy+(2.0*phi[2]-2.0*n[2])*alpha); 
  out_n[3] += (phi[3]-1.0*n[3])*alpha; 
} 
