#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH double lbo_gyrokinetic_drag_vol_1x1v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[2]:      cell-center coordinates. 
  // dxv[2]:    cell spacing. 
  // m_:        species mass.
  // bmag_inv:  1/(magnetic field magnitude). 
  // nuSum:     collisionalities added (self and cross species collisionalities). 
  // nuUSum:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum: sum of thermal speeds squared time their respective collisionalities. 
  // f:         input distribution function.
  // out:       incremented output 
  double rdv2[1]; 
  rdv2[0]   = 2.0/dxv[1]; 

  double alphaDrag[6]; 
  // Expand rdv2*(nu*vpar-nuUparSum) in phase basis.
  alphaDrag[0] = rdv2[0]*(1.414213562373095*nuUSum[0]-1.414213562373095*nuSum[0]*w[1]); 
  alphaDrag[1] = rdv2[0]*(1.414213562373095*nuUSum[1]-1.414213562373095*nuSum[1]*w[1]); 
  alphaDrag[2] = -0.408248290463863*nuSum[0]*rdv2[0]*dxv[1]; 
  alphaDrag[3] = -0.408248290463863*rdv2[0]*dxv[1]*nuSum[1]; 

  out[2] += 0.8660254037844386*(alphaDrag[3]*f[3]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[3] += 0.8660254037844386*(alphaDrag[2]*f[3]+f[2]*alphaDrag[3]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[4] += 1.732050807568877*(alphaDrag[3]*f[5]+alphaDrag[2]*f[4])+1.936491673103709*(alphaDrag[1]*f[3]+f[1]*alphaDrag[3]+alphaDrag[0]*f[2]+f[0]*alphaDrag[2]); 
  out[5] += 1.732050807568877*(alphaDrag[2]*f[5]+alphaDrag[3]*f[4])+1.936491673103709*(alphaDrag[0]*f[3]+f[0]*alphaDrag[3]+alphaDrag[1]*f[2]+f[1]*alphaDrag[2]); 

  return fabs(0.25*alphaDrag[0]); 

} 
