#include <math.h>
#include <gkyl_flash_compute_cfl_kernels.h>

double flash_compute_cfl_lbo_pkpm_vol_1x1v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double* nuPrimMomsSum) 
{
  // w[NDIM]:         Cell-center coordinates. 
  // dxv[NDIM]:       Cell spacing. 
  // nuSum:         Collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum: Sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 

  const double dvpar = dxv[1], wvpar = w[1]; 
  const double dv1par = 2.0/dxv[1]; 

  const double dv1parSq4 = 4.0/(dxv[1]*dxv[1]); 
  const double *sumNuUPar = &nuPrimMomsSum[0]; 
  const double *nuVtSqSum = &nuPrimMomsSum[2];

  double alphaDrag[6]; 
  alphaDrag[0] = dv1par*(1.4142135623730951*sumNuUPar[0]-1.4142135623730951*nuSum[0]*wvpar); 
  alphaDrag[1] = dv1par*(1.4142135623730951*sumNuUPar[1]-1.4142135623730951*nuSum[1]*wvpar); 
  alphaDrag[2] = -(0.408248290463863*nuSum[0]*dv1par*dvpar); 
  alphaDrag[3] = -(0.408248290463863*nuSum[1]*dv1par*dvpar); 

  double cflFreq_drag = fabs(1.25*alphaDrag[0]); 
  double cflFreq_diff = fabs(6.363961030678928*nuVtSqSum[0]*dv1parSq4); 

  return fabs(cflFreq_drag) + fabs(cflFreq_diff); 
} 
