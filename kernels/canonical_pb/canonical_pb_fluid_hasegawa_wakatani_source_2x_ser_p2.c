#include <gkyl_canonical_pb_kernels.h> 
void canonical_pb_fluid_hasegawa_wakatani_source_2x_ser_p2(const double *dxv, double alpha, double kappa, const double *background_n_gradient, const double *phi, const double *f, double* GKYL_RESTRICT out) 
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
  out_n[0] += -(0.5*((15.0*background_n_gradient[6]*phi[7]-15.0*phi[6]*background_n_gradient[7]+6.708203932499369*background_n_gradient[3]*phi[5]-6.708203932499369*phi[3]*background_n_gradient[5]-6.708203932499369*background_n_gradient[3]*phi[4]+6.708203932499369*phi[3]*background_n_gradient[4]+3.0*background_n_gradient[1]*phi[2]-3.0*phi[1]*background_n_gradient[2])*dx*dy+(2.0*n[0]-2.0*phi[0])*alpha)); 
  out_n[1] += -(0.1*((33.54101966249684*background_n_gradient[3]*phi[7]-33.54101966249684*phi[3]*background_n_gradient[7]+(30.000000000000004*background_n_gradient[4]-75.00000000000001*background_n_gradient[5])*phi[6]+(75.00000000000001*phi[5]-30.000000000000004*phi[4])*background_n_gradient[6]-33.54101966249685*background_n_gradient[2]*phi[4]+33.54101966249685*phi[2]*background_n_gradient[4]+15.0*background_n_gradient[1]*phi[3]-15.0*phi[1]*background_n_gradient[3])*dx*dy+(10.0*n[1]-10.0*phi[1])*alpha)); 
  out_n[2] += 0.1*(((30.000000000000004*background_n_gradient[5]-75.00000000000001*background_n_gradient[4])*phi[7]+(75.00000000000001*phi[4]-30.000000000000004*phi[5])*background_n_gradient[7]+33.54101966249684*background_n_gradient[3]*phi[6]-33.54101966249684*phi[3]*background_n_gradient[6]-33.54101966249685*background_n_gradient[1]*phi[5]+33.54101966249685*phi[1]*background_n_gradient[5]+15.0*background_n_gradient[2]*phi[3]-15.0*phi[2]*background_n_gradient[3])*dx*dy+(10.0*phi[2]-10.0*n[2])*alpha); 
  out_n[3] += -(0.5*((6.7082039324993685*background_n_gradient[1]*phi[7]-6.7082039324993685*phi[1]*background_n_gradient[7]-6.7082039324993685*background_n_gradient[2]*phi[6]+6.7082039324993685*phi[2]*background_n_gradient[6]+15.0*background_n_gradient[4]*phi[5]-15.0*phi[4]*background_n_gradient[5])*dx*dy+(2.0*n[3]-2.0*phi[3])*alpha)); 
  out_n[4] += -(0.1*((67.0820393249937*background_n_gradient[6]*phi[7]-67.0820393249937*phi[6]*background_n_gradient[7]+15.000000000000002*background_n_gradient[1]*phi[6]-15.000000000000002*phi[1]*background_n_gradient[6]-30.0*background_n_gradient[3]*phi[4]+30.0*phi[3]*background_n_gradient[4])*dx*dy+(10.0*n[4]-10.0*phi[4])*alpha)); 
  out_n[5] += -(0.1*(((67.0820393249937*background_n_gradient[6]-15.000000000000002*background_n_gradient[2])*phi[7]+(15.000000000000002*phi[2]-67.0820393249937*phi[6])*background_n_gradient[7]+30.0*background_n_gradient[3]*phi[5]-30.0*phi[3]*background_n_gradient[5])*dx*dy+(10.0*n[5]-10.0*phi[5])*alpha)); 
  out_n[6] += -(0.5*((13.416407864998739*background_n_gradient[4]*phi[7]-13.416407864998739*phi[4]*background_n_gradient[7]-3.0*background_n_gradient[3]*phi[6]+3.0*phi[3]*background_n_gradient[6])*dx*dy+(2.0*n[6]-2.0*phi[6])*alpha)); 
  out_n[7] += -(0.5*((3.0*background_n_gradient[3]*phi[7]-3.0*phi[3]*background_n_gradient[7]-13.416407864998739*background_n_gradient[5]*phi[6]+13.416407864998739*phi[5]*background_n_gradient[6])*dx*dy+(2.0*n[7]-2.0*phi[7])*alpha)); 
} 
