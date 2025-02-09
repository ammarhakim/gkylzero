#include <gkyl_canonical_pb_kernels.h> 
void canonical_pb_fluid_hasegawa_mima_source_2x_ser_p1(const double *dxv, double alpha, double kappa, const double *background_n_gradient, const double *phi, const double *f, double* GKYL_RESTRICT out) 
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
  out[0] += -(0.5*(3.0*background_n_gradient[1]*phi[2]-3.0*phi[1]*background_n_gradient[2])*dx*dy); 
  out[1] += -(0.5*(3.0*background_n_gradient[1]*phi[3]-3.0*phi[1]*background_n_gradient[3])*dx*dy); 
  out[2] += 0.5*(3.0*background_n_gradient[2]*phi[3]-3.0*phi[2]*background_n_gradient[3])*dx*dy; 
} 
