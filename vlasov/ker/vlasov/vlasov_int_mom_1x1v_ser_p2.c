#include <gkyl_mom_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_int_five_moments_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[0]*dxv[1]*0.25; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += volFact*(2.0*f[0]*wx1+0.5773502691896258*f[2]*dv1); 
  out[2] += volFact*(2.0*f[0]*wx1_sq+1.1547005383792517*f[2]*dv1*wx1+0.14907119849998596*f[5]*dv1_sq+0.16666666666666666*f[0]*dv1_sq); 
} 
