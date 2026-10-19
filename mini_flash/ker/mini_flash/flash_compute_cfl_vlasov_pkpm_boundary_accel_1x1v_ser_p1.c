#include <math.h>
#include <gkyl_flash_compute_cfl_kernels.h>

double flash_compute_cfl_vlasov_pkpm_boundary_accel_1x1v_ser_p1(const double *w, const double *dxv, 
  const double *div_b, const double *pkpm_accel_vars, const int edge) 
{
  // w[NDIM]:         Cell-center coordinates. 
  // dxv[NDIM]:       Cell spacing. 
  // edge:            Determines if the update is for the left edge (-1) or right edge (+1). 
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

  double alphaSurf[2] = {0.0}; 

  // get stable timestep of alpha_v = 1/rho (div(p_par b) - p_perp div(b)) - v_par bb : grad(u) 
  // from the quadrature point evaluation needed to compute upwinded distribution functions 
  double cflFreq = 0.0;
  double alphaOrd = 0.0;

  if (edge == -1) { 
    alphaOrd = 0.7071067811865475*alphaSurf_l[0]-0.7071067811865475*alphaSurf_l[1];
    cflFreq = fmax(cflFreq, fabs(alphaOrd));

    alphaOrd = 0.7071067811865475*(alphaSurf_l[1]+alphaSurf_l[0]);
    cflFreq = fmax(cflFreq, fabs(alphaOrd));

  } else { 
    alphaOrd = 0.7071067811865475*alphaSurf_r[0]-0.7071067811865475*alphaSurf_r[1];
    cflFreq = fmax(cflFreq, fabs(alphaOrd));

    alphaOrd = 0.7071067811865475*(alphaSurf_r[1]+alphaSurf_r[0]);
    cflFreq = fmax(cflFreq, fabs(alphaOrd));

  }
  return 2.5*dv1par*cflFreq;
} 
