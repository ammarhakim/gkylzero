#include <gkyl_canonical_pb_kernels.h> 
void canonical_pb_fluid_hasegawa_mima_source_2x_ser_p2(const double *dxv, double alpha, double kappa, const double *background_n_gradient, const double *phi, const double *f, double* GKYL_RESTRICT out) 
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
  out[0] += -(7.5*background_n_gradient[6]*phi[7]*dx*dy)+7.5*phi[6]*background_n_gradient[7]*dx*dy-3.3541019662496847*background_n_gradient[3]*phi[5]*dx*dy+3.3541019662496847*phi[3]*background_n_gradient[5]*dx*dy+3.3541019662496847*background_n_gradient[3]*phi[4]*dx*dy-3.3541019662496847*phi[3]*background_n_gradient[4]*dx*dy-1.5*background_n_gradient[1]*phi[2]*dx*dy+1.5*phi[1]*background_n_gradient[2]*dx*dy; 
  out[1] += -(3.3541019662496843*background_n_gradient[3]*phi[7]*dx*dy)+3.3541019662496843*phi[3]*background_n_gradient[7]*dx*dy+7.500000000000001*background_n_gradient[5]*phi[6]*dx*dy-3.0*background_n_gradient[4]*phi[6]*dx*dy-7.500000000000001*phi[5]*background_n_gradient[6]*dx*dy+3.0*phi[4]*background_n_gradient[6]*dx*dy+3.3541019662496847*background_n_gradient[2]*phi[4]*dx*dy-3.3541019662496847*phi[2]*background_n_gradient[4]*dx*dy-1.5*background_n_gradient[1]*phi[3]*dx*dy+1.5*phi[1]*background_n_gradient[3]*dx*dy; 
  out[2] += 3.0*background_n_gradient[5]*phi[7]*dx*dy-7.500000000000001*background_n_gradient[4]*phi[7]*dx*dy-3.0*phi[5]*background_n_gradient[7]*dx*dy+7.500000000000001*phi[4]*background_n_gradient[7]*dx*dy+3.3541019662496843*background_n_gradient[3]*phi[6]*dx*dy-3.3541019662496843*phi[3]*background_n_gradient[6]*dx*dy-3.3541019662496847*background_n_gradient[1]*phi[5]*dx*dy+3.3541019662496847*phi[1]*background_n_gradient[5]*dx*dy+1.5*background_n_gradient[2]*phi[3]*dx*dy-1.5*phi[2]*background_n_gradient[3]*dx*dy; 
  out[3] += -(3.3541019662496843*background_n_gradient[1]*phi[7]*dx*dy)+3.3541019662496843*phi[1]*background_n_gradient[7]*dx*dy+3.3541019662496843*background_n_gradient[2]*phi[6]*dx*dy-3.3541019662496843*phi[2]*background_n_gradient[6]*dx*dy-7.5*background_n_gradient[4]*phi[5]*dx*dy+7.5*phi[4]*background_n_gradient[5]*dx*dy; 
  out[4] += -(6.708203932499369*background_n_gradient[6]*phi[7]*dx*dy)+6.708203932499369*phi[6]*background_n_gradient[7]*dx*dy-1.5*background_n_gradient[1]*phi[6]*dx*dy+1.5*phi[1]*background_n_gradient[6]*dx*dy+3.0*background_n_gradient[3]*phi[4]*dx*dy-3.0*phi[3]*background_n_gradient[4]*dx*dy; 
  out[5] += -(6.708203932499369*background_n_gradient[6]*phi[7]*dx*dy)+1.5*background_n_gradient[2]*phi[7]*dx*dy+6.708203932499369*phi[6]*background_n_gradient[7]*dx*dy-1.5*phi[2]*background_n_gradient[7]*dx*dy-3.0*background_n_gradient[3]*phi[5]*dx*dy+3.0*phi[3]*background_n_gradient[5]*dx*dy; 
  out[6] += -(6.708203932499369*background_n_gradient[4]*phi[7]*dx*dy)+6.708203932499369*phi[4]*background_n_gradient[7]*dx*dy+1.5*background_n_gradient[3]*phi[6]*dx*dy-1.5*phi[3]*background_n_gradient[6]*dx*dy; 
  out[7] += -(1.5*background_n_gradient[3]*phi[7]*dx*dy)+1.5*phi[3]*background_n_gradient[7]*dx*dy+6.708203932499369*background_n_gradient[5]*phi[6]*dx*dy-6.708203932499369*phi[5]*background_n_gradient[6]*dx*dy; 
} 
