#include <gkyl_vlasov_mom_kernels.h> 
void vlasov_int_mom_1x2v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out) 
{ 
  const gkyl_real volFact = dxv[0]*dxv[1]*dxv[2]*0.125; 
  const gkyl_real wx1 = w[1], dv1 = dxv[1]; 
  const gkyl_real wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const gkyl_real wx2 = w[2], dv2 = dxv[2]; 
  const gkyl_real wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
 
  out[0] += 2.828427124746191*f[0]*volFact; 
  out[1] += volFact*(2.828427124746191*f[0]*wx1+0.8164965809277261*f[2]*dv1); 
  out[2] += volFact*(2.828427124746191*f[0]*wx2+0.8164965809277261*f[3]*dv2); 
  out[3] += volFact*(2.828427124746191*f[0]*wx2_sq+1.632993161855453*f[3]*dv2*wx2+2.828427124746191*f[0]*wx1_sq+1.632993161855453*f[2]*dv1*wx1+0.2357022603955158*f[0]*dv2_sq+0.2357022603955158*f[0]*dv1_sq); 
} 
