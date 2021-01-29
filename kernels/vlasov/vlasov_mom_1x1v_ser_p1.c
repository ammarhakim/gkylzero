#include <gkyl_vlasov_mom_kernels.h> 
void vlasov_M0_1x1v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out) 
{ 
  const gkyl_real volFact = dxv[1]/2; 
  out[0] += 1.414213562373095*f[0]*volFact; 
  out[1] += 1.414213562373095*f[1]*volFact; 
} 
void vlasov_M1i_1x1v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out) 
{ 
  const gkyl_real volFact = dxv[1]/2; 
  const gkyl_real wx1 = w[1], dv1 = dxv[1]; 
  out[0] += volFact*(1.414213562373095*f[0]*wx1+0.408248290463863*f[2]*dv1); 
  out[1] += volFact*(1.414213562373095*f[1]*wx1+0.408248290463863*f[3]*dv1); 
} 
void vlasov_M2_1x1v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out) 
{ 
  const gkyl_real volFact = dxv[1]/2; 
  const gkyl_real wx1 = w[1], dv1 = dxv[1]; 
  const gkyl_real wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  out[0] += volFact*(1.414213562373095*f[0]*wx1_sq+0.8164965809277261*f[2]*dv1*wx1+0.1178511301977579*f[0]*dv1_sq); 
  out[1] += volFact*(1.414213562373095*f[1]*wx1_sq+0.8164965809277261*f[3]*dv1*wx1+0.1178511301977579*f[1]*dv1_sq); 
} 
void vlasov_FiveMoments_1x1v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict outM0, gkyl_real* restrict outM1i, gkyl_real* restrict outM2) 
{ 
  const gkyl_real volFact = dxv[1]/2; 
  const gkyl_real wx1 = w[1], dv1 = dxv[1]; 
  const gkyl_real wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  gkyl_real tempM0[2], tempM1i[2]; 

  tempM0[0] = 1.414213562373095*f[0]*volFact; 
  tempM0[1] = 1.414213562373095*f[1]*volFact; 

  tempM1i[0] = tempM0[0]*wx1+0.408248290463863*f[2]*dv1*volFact; 
  tempM1i[1] = tempM0[1]*wx1+0.408248290463863*f[3]*dv1*volFact; 

  outM0[0] += tempM0[0]; 
  outM0[1] += tempM0[1]; 
  outM1i[0] += tempM1i[0]; 
  outM1i[1] += tempM1i[1]; 
  outM2[0] += (-1.0*tempM0[0]*wx1_sq)+2.0*tempM1i[0]*wx1+0.1178511301977579*f[0]*dv1_sq*volFact; 
  outM2[1] += (-1.0*tempM0[1]*wx1_sq)+2.0*tempM1i[1]*wx1+0.1178511301977579*f[1]*dv1_sq*volFact; 
} 
void vlasov_M2ij_1x1v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out) 
{ 
  const gkyl_real volFact = dxv[1]/2; 
  const gkyl_real wx1 = w[1], dv1 = dxv[1]; 
  const gkyl_real wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  out[0] += volFact*(1.414213562373095*f[0]*wx1_sq+0.8164965809277261*f[2]*dv1*wx1+0.1178511301977579*f[0]*dv1_sq); 
  out[1] += volFact*(1.414213562373095*f[1]*wx1_sq+0.8164965809277261*f[3]*dv1*wx1+0.1178511301977579*f[1]*dv1_sq); 
} 
void vlasov_M3i_1x1v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out) 
{ 
  const gkyl_real volFact = dxv[1]/2; 
  const gkyl_real wx1 = w[1], dv1 = dxv[1]; 
  const gkyl_real wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const gkyl_real wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  out[0] += volFact*(1.414213562373095*f[0]*wx1*wx1_sq+1.224744871391589*f[2]*dv1*wx1_sq+0.3535533905932737*f[0]*dv1_sq*wx1+0.06123724356957942*f[2]*dv1*dv1_sq); 
  out[1] += volFact*(1.414213562373095*f[1]*wx1*wx1_sq+1.224744871391589*f[3]*dv1*wx1_sq+0.3535533905932737*f[1]*dv1_sq*wx1+0.06123724356957942*f[3]*dv1*dv1_sq); 
} 
void vlasov_M3ijk_1x1v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out) 
{ 
  const gkyl_real volFact = dxv[1]/2; 
  const gkyl_real wx1 = w[1], dv1 = dxv[1]; 
  const gkyl_real wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const gkyl_real wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  out[0] += volFact*(1.414213562373095*f[0]*wx1*wx1_sq+1.224744871391589*f[2]*dv1*wx1_sq+0.3535533905932737*f[0]*dv1_sq*wx1+0.06123724356957942*f[2]*dv1*dv1_sq); 
  out[1] += volFact*(1.414213562373095*f[1]*wx1*wx1_sq+1.224744871391589*f[3]*dv1*wx1_sq+0.3535533905932737*f[1]*dv1_sq*wx1+0.06123724356957942*f[3]*dv1*dv1_sq); 
} 
