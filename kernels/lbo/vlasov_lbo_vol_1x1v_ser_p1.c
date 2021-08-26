#include <gkyl_vlasov_lbo_kernels.h> 
GKYL_CU_DH double vlasov_lbo_vol_1x1v_ser_p1(const double *w, const double *dxv, const double nuSum, const double *nuUSum, const double *nuVtSqSum, const double *f, double GKYL_RESTRICT *out) 
{ 
  // w[2]:      Cell-center coordinates. 
  // dxv[2]:    Cell spacing. 
  // nuSum:     collisionalities added (self and cross species collisionalities). 
  // nuUSum:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum: sum of thermal speeds squared time their respective collisionalities. 
  // f:         Input distribution function.
  // out:       Incremented output 
  const double rdvx2 = 2.0/dxv[1]; 
  const double rdvxSq4 = 4.0/(dxv[1]*dxv[1]); 

  double alphaDrag[4]; 
  // Expand rdv2*(nu*vx-nuUSumx) in phase basis.
  alphaDrag[0] = (1.414213562373095*nuUSum[0]-2.0*w[1]*nuSum)*rdvx2; 
  alphaDrag[1] = 1.414213562373095*nuUSum[1]*rdvx2; 
  alphaDrag[2] = -0.5773502691896258*dxv[1]*nuSum*rdvx2; 

  // Put together updates due to drag and diffusion terms.
  out[2] += 0.8660254037844386*(alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[3] += 0.8660254037844386*(alphaDrag[2]*f[3]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 

  return fabs(0.25*alphaDrag[0])+fabs(0.9428090415820636*nuVtSqSum[0]*rdvxSq4); 

} 
