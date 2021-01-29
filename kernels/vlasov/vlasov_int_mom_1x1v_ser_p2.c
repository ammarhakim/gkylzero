#include <gkyl_vlasov_mom_kernels.h> 
void vlasov_int_mom_1x1v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out) 
{ 
  const gkyl_real volFact = dxv[0]*dxv[1]*0.25; 
  const gkyl_real wx1 = w[1], dv1 = dxv[1]; 
  const gkyl_real wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += volFact*(2.0*f[0]*wx1+0.5773502691896258*f[2]*dv1); 
  out[2] += volFact*(2.0*f[0]*wx1_sq+1.154700538379252*f[2]*dv1*wx1+0.149071198499986*f[5]*dv1_sq+0.1666666666666667*f[0]*dv1_sq); 
} 
