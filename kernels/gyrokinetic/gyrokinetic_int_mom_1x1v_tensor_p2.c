#include <gkyl_mom_gyrokinetic_kernels.h> 
GKYL_CU_DH void gyrokinetic_int_mom_1x1v_tensor_p2(const double *w, const double *dxv, const int *idx, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += volFact*(2.0*f[0]*wx1+0.5773502691896258*f[2]*dv1); 
  out[2] += volFact*(2.0*f[0]*wx1_sq+1.154700538379252*f[2]*dv1*wx1+0.149071198499986*f[5]*dv1_sq+0.1666666666666667*f[0]*dv1_sq); 
} 