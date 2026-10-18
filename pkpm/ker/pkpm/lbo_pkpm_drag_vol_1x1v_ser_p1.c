#include <gkyl_lbo_pkpm_kernels.h> 
GKYL_CU_DH double lbo_pkpm_drag_vol_1x1v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:       Cell-center coordinates. 
  // dxv[NDIM]:     Cell spacing. 
  // nuSum:         Collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum: Sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 
  // f:             Input distribution functions [F_0, T_perp/m G_1 = T_perp/m (F_0 - F_1)].
  // out:           Incremented output distribution functions. 
  const double rdvpar = 2.0/dxv[1]; 
  const double dvpar = dxv[1], wvpar = w[1]; 
  const double *F_0 = &f[0]; 
  const double *G_1 = &f[6]; 
  double *out_F_0 = &out[0]; 
  double *out_G_1 = &out[6]; 

  const double *sumNuUPar = &nuPrimMomsSum[0]; 

  double alphaDrag[6]; 
  alphaDrag[0] = rdvpar*(1.4142135623730951*sumNuUPar[0]-1.4142135623730951*nuSum[0]*wvpar); 
  alphaDrag[1] = rdvpar*(1.4142135623730951*sumNuUPar[1]-1.4142135623730951*nuSum[1]*wvpar); 
  alphaDrag[2] = -(0.408248290463863*nuSum[0]*dvpar*rdvpar); 
  alphaDrag[3] = -(0.408248290463863*nuSum[1]*dvpar*rdvpar); 

  const double *sumNuVtSq = &nuPrimMomsSum[2]; 
  out_F_0[2] += 0.8660254037844386*F_0[3]*alphaDrag[3]+0.8660254037844386*F_0[2]*alphaDrag[2]+0.8660254037844386*F_0[1]*alphaDrag[1]+0.8660254037844386*F_0[0]*alphaDrag[0]; 
  out_F_0[3] += 0.8660254037844386*F_0[2]*alphaDrag[3]+0.8660254037844386*alphaDrag[2]*F_0[3]+0.8660254037844386*F_0[0]*alphaDrag[1]+0.8660254037844386*alphaDrag[0]*F_0[1]; 
  out_F_0[4] += 1.7320508075688774*alphaDrag[3]*F_0[5]+1.7320508075688772*alphaDrag[2]*F_0[4]+1.9364916731037085*F_0[1]*alphaDrag[3]+1.9364916731037085*alphaDrag[1]*F_0[3]+1.9364916731037085*F_0[0]*alphaDrag[2]+1.9364916731037085*alphaDrag[0]*F_0[2]; 
  out_F_0[5] += 1.7320508075688772*alphaDrag[2]*F_0[5]+1.7320508075688774*alphaDrag[3]*F_0[4]+1.9364916731037085*F_0[0]*alphaDrag[3]+1.9364916731037085*alphaDrag[0]*F_0[3]+1.9364916731037085*F_0[1]*alphaDrag[2]+1.9364916731037085*alphaDrag[1]*F_0[2]; 
  out_G_1[0] += 1.4142135623730951*F_0[1]*sumNuVtSq[1]-1.4142135623730951*G_1[1]*nuSum[1]+1.4142135623730951*F_0[0]*sumNuVtSq[0]-1.4142135623730951*G_1[0]*nuSum[0]; 
  out_G_1[1] += 1.4142135623730951*F_0[0]*sumNuVtSq[1]-1.4142135623730951*G_1[0]*nuSum[1]-1.4142135623730951*nuSum[0]*G_1[1]+1.4142135623730951*sumNuVtSq[0]*F_0[1]; 
  out_G_1[2] += 0.8660254037844386*G_1[3]*alphaDrag[3]-1.4142135623730951*nuSum[1]*G_1[3]+1.4142135623730951*sumNuVtSq[1]*F_0[3]+0.8660254037844386*G_1[2]*alphaDrag[2]-1.4142135623730951*nuSum[0]*G_1[2]+1.4142135623730951*sumNuVtSq[0]*F_0[2]+0.8660254037844386*G_1[1]*alphaDrag[1]+0.8660254037844386*G_1[0]*alphaDrag[0]; 
  out_G_1[3] += 0.8660254037844386*G_1[2]*alphaDrag[3]+0.8660254037844386*alphaDrag[2]*G_1[3]-1.4142135623730951*nuSum[0]*G_1[3]+1.4142135623730951*sumNuVtSq[0]*F_0[3]-1.4142135623730951*nuSum[1]*G_1[2]+1.4142135623730951*sumNuVtSq[1]*F_0[2]+0.8660254037844386*G_1[0]*alphaDrag[1]+0.8660254037844386*alphaDrag[0]*G_1[1]; 
  out_G_1[4] += 1.7320508075688774*alphaDrag[3]*G_1[5]-1.4142135623730951*nuSum[1]*G_1[5]+1.4142135623730951*sumNuVtSq[1]*F_0[5]+1.7320508075688772*alphaDrag[2]*G_1[4]-1.4142135623730951*nuSum[0]*G_1[4]+1.4142135623730951*sumNuVtSq[0]*F_0[4]+1.9364916731037085*G_1[1]*alphaDrag[3]+1.9364916731037085*alphaDrag[1]*G_1[3]+1.9364916731037085*G_1[0]*alphaDrag[2]+1.9364916731037085*alphaDrag[0]*G_1[2]; 
  out_G_1[5] += 1.7320508075688772*alphaDrag[2]*G_1[5]-1.4142135623730951*nuSum[0]*G_1[5]+1.4142135623730951*sumNuVtSq[0]*F_0[5]+1.7320508075688774*alphaDrag[3]*G_1[4]-1.4142135623730951*nuSum[1]*G_1[4]+1.4142135623730951*sumNuVtSq[1]*F_0[4]+1.9364916731037085*G_1[0]*alphaDrag[3]+1.9364916731037085*alphaDrag[0]*G_1[3]+1.9364916731037085*G_1[1]*alphaDrag[2]+1.9364916731037085*alphaDrag[1]*G_1[2]; 

  return fabs(1.25*alphaDrag[0]) + 2.0*fabs(0.7071067811865475*nuSum[0]); 

} 
