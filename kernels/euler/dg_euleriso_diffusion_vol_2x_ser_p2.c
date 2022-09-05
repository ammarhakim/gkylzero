#include <dg_euleriso_diffusion_kernels.h>

GKYL_CU_DH double
dg_euleriso_diffusion_vol_2x_ser_p2(const double* w, const double* dx,
  const double* D, const double* uvar, const double* statevec, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dxv[NDIM]: Cell spacing
  // D: Diffusion tensor (rho diffusion, rhoux diffusion, rhouyz diffusion, rhouz diffusion)- advection eq (rho) has no diffusion, first pOrder values should be zero
  // uvar: unaliased form of rhoui/ rho
  // q: Input state [rho rhoux rhouy rhouz]
  // out: Incremented output

  const double retval;
  const double *D1 = &D[0]; 
  const double *D2 = &D[8]; 
  const double *D3 = &D[16]; 
  const double *D4 = &D[24]; 
  double mu = D2[0]; 
  double dx10 = 2./dx[0]; 
  double dx11 = 2./dx[1]; 
  double cfl = mu*9*(dx10*dx10);
  return cfl;
}
