#include <math.h>
#include <gkyl_flash_compute_cfl_kernels.h>

double flash_compute_cfl_vlasov_pkpm_surf_accel_1x1v_ser_p1(const double *w, const double *dxv, 
  const double *pkpm_accel_vars) 
{
  // w[NDIM]:         Cell-center coordinates. 
  // dxv[NDIM]:       Cell spacing. 
  const double dv1par = 2.0/dxv[1]; 
  const double dvpar = dxv[1], wvpar = w[1]; 

  const double *bb_grad_u = &pkpm_accel_vars[2]; 
  const double *p_force = &pkpm_accel_vars[4]; 

  double alphaSurf_l[2] = {0.0}; 
  alphaSurf_l[0] = -(1.0*bb_grad_u[0]*wvpar)+0.5*bb_grad_u[0]*dvpar+p_force[0]; 
  alphaSurf_l[1] = -(1.0*bb_grad_u[1]*wvpar)+0.5*bb_grad_u[1]*dvpar+p_force[1]; 

  double alphaSurf_r[2] = {0.0}; 
  alphaSurf_r[0] = -(1.0*bb_grad_u[0]*wvpar)-0.5*bb_grad_u[0]*dvpar+p_force[0]; 
  alphaSurf_r[1] = -(1.0*bb_grad_u[1]*wvpar)-0.5*bb_grad_u[1]*dvpar+p_force[1]; 

  // get stable timestep of alpha_v = 1/rho (div(p_par b) - p_perp div(b)) - v_par bb : grad(u) 
  // from the quadrature point evaluation needed to compute upwinded distribution functions 
  double cflFreq = 0.0;
  double alphaOrd_l, alphaOrd_r = 0.0;

  alphaOrd_l = 0.7071067811865475*alphaSurf_l[0]-0.7071067811865475*alphaSurf_l[1];
  alphaOrd_r = 0.7071067811865475*alphaSurf_r[0]-0.7071067811865475*alphaSurf_r[1];
  cflFreq = fmax(cflFreq, fmax(fabs(alphaOrd_l), fabs(alphaOrd_r)));

  alphaOrd_l = 0.7071067811865475*(alphaSurf_l[1]+alphaSurf_l[0]);
  alphaOrd_r = 0.7071067811865475*(alphaSurf_r[1]+alphaSurf_r[0]);
  cflFreq = fmax(cflFreq, fmax(fabs(alphaOrd_l), fabs(alphaOrd_r)));

  return 2.5*dv1par*cflFreq;
} 
