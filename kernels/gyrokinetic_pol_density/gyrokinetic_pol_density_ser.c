#include <gkyl_gyrokinetic_pol_density_kernels.h>

void gkyl_gyrokinetic_pol_density_2x_ser_p1(const double *dx, const double *epsilon, const double *phi, double *out) 
{ 
  // dx: cell lengths.
  // epsilon: polarization weight, sum_s jacobgeo*(n_s*m_s/B^2)*g^{ij}.
  // phi: electrostatic potential.
  // out: polarization density.

  const double *epsxx = &epsilon[0];

  double rdx00 = 4.0/(dx[0]*dx[0]);

  out[0] = -(13.74772708486752*epsxx[3]*phi[11]*rdx00)-13.747727084867519*epsxx[1]*phi[8]*rdx00-3.3541019662496843*epsxx[2]*phi[6]*rdx00-3.3541019662496847*epsxx[0]*phi[4]*rdx00-1.5*epsxx[3]*phi[3]*rdx00-1.5*epsxx[1]*phi[1]*rdx00; 
  out[1] = -(11.4564392373896*epsxx[2]*phi[11]*rdx00)-11.4564392373896*epsxx[0]*phi[8]*rdx00-6.708203932499368*epsxx[3]*phi[6]*rdx00-6.708203932499369*epsxx[1]*phi[4]*rdx00; 
  out[2] = -(12.296340919151517*epsxx[3]*phi[13]*rdx00)-13.74772708486752*epsxx[1]*phi[11]*rdx00-3.0*epsxx[2]*phi[10]*rdx00-13.747727084867519*epsxx[3]*phi[8]*rdx00-1.3416407864998738*epsxx[3]*phi[7]*rdx00-3.3541019662496843*epsxx[0]*phi[6]*rdx00-3.3541019662496847*epsxx[2]*phi[4]*rdx00-1.5*epsxx[1]*phi[3]*rdx00-1.5*phi[1]*epsxx[3]*rdx00; 
  out[3] = -(10.246950765959596*epsxx[2]*phi[13]*rdx00)-11.4564392373896*epsxx[0]*phi[11]*rdx00-6.0*epsxx[3]*phi[10]*rdx00-11.4564392373896*epsxx[2]*phi[8]*rdx00-6.708203932499368*epsxx[1]*phi[6]*rdx00-6.708203932499369*epsxx[3]*phi[4]*rdx00; 

}
