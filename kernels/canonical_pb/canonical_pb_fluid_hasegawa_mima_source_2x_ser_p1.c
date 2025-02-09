#include <gkyl_canonical_pb_kernels.h> 
void canonical_pb_fluid_hasegawa_mima_source_2x_ser_p1(const double *dxv, double alpha, double kappa, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  // dxv[NDIM]: cell length.
  // alpha: Adiabaticity parameter for adiabatic coupling of vorticity and density (zero for Hasegawa-Mima).
  // kappa: Constant density gradient scale length.
  // phi: Potential in fluid system given by grad^2 phi = f where f is (one of) the evolved quantities.
  // f: input state vector in center cell.
  // out: output increment in center cell.
  double dy = 2.0/dxv[1]; 
  out[0] += 1.7320508075688772*phi[2]*dy*kappa; 
  out[1] += 1.7320508075688772*phi[3]*dy*kappa; 
} 
