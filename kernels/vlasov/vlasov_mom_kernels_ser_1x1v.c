#include <gkyl_vlasov_mom_kernels.h>

void vlasov_mom_1x1v_m0_p1(const double *w, const double *dxv, const int *idx, const double *restrict f, double* restrict out) 
{ 
  const double volFact = dxv[1]/2; 
  out[0] += 1.414213562373095*f[0]*volFact; 
  out[1] += 1.414213562373095*f[1]*volFact; 
} 
void vlasov_mom_1x1v_m0_p2(const double *w, const double *dxv, const int *idx, const double *restrict f, double* restrict out) 
{ 
  const double volFact = dxv[1]/2; 
  out[0] += 1.414213562373095*f[0]*volFact; 
  out[1] += 1.414213562373095*f[1]*volFact; 
  out[2] += 1.414213562373095*f[4]*volFact; 
} 
void vlasov_mom_1x1v_m1i_p1(const double *w, const double *dxv, const int *idx, const double *restrict f, double* restrict out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  out[0] += volFact*(1.414213562373095*f[0]*wx1+0.408248290463863*f[2]*dv1); 
  out[1] += volFact*(1.414213562373095*f[1]*wx1+0.408248290463863*f[3]*dv1); 
} 
void vlasov_mom_1x1v_m1i_p2(const double *w, const double *dxv, const int *idx, const double *restrict f, double* restrict out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  out[0] += volFact*(1.414213562373095*f[0]*wx1+0.408248290463863*f[2]*dv1); 
  out[1] += volFact*(1.414213562373095*f[1]*wx1+0.408248290463863*f[3]*dv1); 
  out[2] += volFact*(1.414213562373095*f[4]*wx1+0.408248290463863*f[6]*dv1); 
} 
void vlasov_mom_1x1v_m2ij_p1(const double *w, const double *dxv, const int *idx, const double *restrict f, double* restrict out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  out[0] += volFact*(1.414213562373095*f[0]*wx1_sq+0.8164965809277261*f[2]*dv1*wx1+0.1178511301977579*f[0]*dv1_sq); 
  out[1] += volFact*(1.414213562373095*f[1]*wx1_sq+0.8164965809277261*f[3]*dv1*wx1+0.1178511301977579*f[1]*dv1_sq); 
} 
void vlasov_mom_1x1v_m2ij_p2(const double *w, const double *dxv, const int *idx, const double *restrict f, double* restrict out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  out[0] += volFact*(1.414213562373095*f[0]*wx1_sq+0.8164965809277261*f[2]*dv1*wx1+0.105409255338946*f[5]*dv1_sq+0.1178511301977579*f[0]*dv1_sq); 
  out[1] += volFact*(1.414213562373095*f[1]*wx1_sq+0.8164965809277261*f[3]*dv1*wx1+0.105409255338946*f[7]*dv1_sq+0.1178511301977579*f[1]*dv1_sq); 
  out[2] += volFact*(1.414213562373095*f[4]*wx1_sq+0.816496580927726*f[6]*dv1*wx1+0.1178511301977579*f[4]*dv1_sq); 
} 
void vlasov_mom_1x1v_m2_p1(const double *w, const double *dxv, const int *idx, const double *restrict f, double* restrict out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  out[0] += volFact*(1.414213562373095*f[0]*wx1_sq+0.8164965809277261*f[2]*dv1*wx1+0.1178511301977579*f[0]*dv1_sq); 
  out[1] += volFact*(1.414213562373095*f[1]*wx1_sq+0.8164965809277261*f[3]*dv1*wx1+0.1178511301977579*f[1]*dv1_sq); 
} 
void vlasov_mom_1x1v_m2_p2(const double *w, const double *dxv, const int *idx, const double *restrict f, double* restrict out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  out[0] += volFact*(1.414213562373095*f[0]*wx1_sq+0.8164965809277261*f[2]*dv1*wx1+0.105409255338946*f[5]*dv1_sq+0.1178511301977579*f[0]*dv1_sq); 
  out[1] += volFact*(1.414213562373095*f[1]*wx1_sq+0.8164965809277261*f[3]*dv1*wx1+0.105409255338946*f[7]*dv1_sq+0.1178511301977579*f[1]*dv1_sq); 
  out[2] += volFact*(1.414213562373095*f[4]*wx1_sq+0.816496580927726*f[6]*dv1*wx1+0.1178511301977579*f[4]*dv1_sq); 
} 
void vlasov_mom_1x1v_m3i_p1(const double *w, const double *dxv, const int *idx, const double *restrict f, double* restrict out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  out[0] += volFact*(1.414213562373095*f[0]*wx1*wx1_sq+1.224744871391589*f[2]*dv1*wx1_sq+0.3535533905932737*f[0]*dv1_sq*wx1+0.06123724356957942*f[2]*dv1*dv1_sq); 
  out[1] += volFact*(1.414213562373095*f[1]*wx1*wx1_sq+1.224744871391589*f[3]*dv1*wx1_sq+0.3535533905932737*f[1]*dv1_sq*wx1+0.06123724356957942*f[3]*dv1*dv1_sq); 
} 
void vlasov_mom_1x1v_m3i_p2(const double *w, const double *dxv, const int *idx, const double *restrict f, double* restrict out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  out[0] += volFact*(1.414213562373095*f[0]*wx1*wx1_sq+1.224744871391589*f[2]*dv1*wx1_sq+0.3162277660168379*f[5]*dv1_sq*wx1+0.3535533905932737*f[0]*dv1_sq*wx1+0.06123724356957942*f[2]*dv1*dv1_sq); 
  out[1] += volFact*(1.414213562373095*f[1]*wx1*wx1_sq+1.224744871391589*f[3]*dv1*wx1_sq+0.3162277660168379*f[7]*dv1_sq*wx1+0.3535533905932737*f[1]*dv1_sq*wx1+0.06123724356957942*f[3]*dv1*dv1_sq); 
  out[2] += volFact*(1.414213562373095*f[4]*wx1*wx1_sq+1.224744871391589*f[6]*dv1*wx1_sq+0.3535533905932737*f[4]*dv1_sq*wx1+0.06123724356957942*f[6]*dv1*dv1_sq); 
} 
