#include <gkyl_mom_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_M0_1x3v_tensor_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  out[0] += 2.828427124746191*f[0]*volFact; 
  out[1] += 2.828427124746191*f[1]*volFact; 
} 
GKYL_CU_DH void vlasov_M1i_1x3v_tensor_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx3 = w[3], dv3 = dxv[3]; 
  out[0] += volFact*(2.828427124746191*f[0]*wx1+0.8164965809277261*f[2]*dv1); 
  out[1] += volFact*(2.828427124746191*f[1]*wx1+0.8164965809277261*f[5]*dv1); 
  out[2] += volFact*(2.828427124746191*f[0]*wx2+0.8164965809277261*f[3]*dv2); 
  out[3] += volFact*(2.828427124746191*f[1]*wx2+0.8164965809277261*f[6]*dv2); 
  out[4] += volFact*(2.828427124746191*f[0]*wx3+0.8164965809277261*f[4]*dv3); 
  out[5] += volFact*(2.828427124746191*f[1]*wx3+0.8164965809277261*f[8]*dv3); 
} 
GKYL_CU_DH void vlasov_M2_1x3v_tensor_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx3 = w[3], dv3 = dxv[3]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
  out[0] += volFact*(2.828427124746191*f[0]*wx3_sq+1.632993161855453*f[4]*dv3*wx3+2.828427124746191*f[0]*wx2_sq+1.632993161855453*f[3]*dv2*wx2+2.828427124746191*f[0]*wx1_sq+1.632993161855453*f[2]*dv1*wx1+0.2357022603955158*f[0]*dv3_sq+0.2357022603955158*f[0]*dv2_sq+0.2357022603955158*f[0]*dv1_sq); 
  out[1] += volFact*(2.828427124746191*f[1]*wx3_sq+1.632993161855453*f[8]*dv3*wx3+2.828427124746191*f[1]*wx2_sq+1.632993161855453*f[6]*dv2*wx2+2.828427124746191*f[1]*wx1_sq+1.632993161855453*f[5]*dv1*wx1+0.2357022603955158*f[1]*dv3_sq+0.2357022603955158*f[1]*dv2_sq+0.2357022603955158*f[1]*dv1_sq); 
} 
GKYL_CU_DH void vlasov_M2ij_1x3v_tensor_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx3 = w[3], dv3 = dxv[3]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
  out[0] += volFact*(2.828427124746191*f[0]*wx1_sq+1.632993161855453*f[2]*dv1*wx1+0.2357022603955158*f[0]*dv1_sq); 
  out[1] += volFact*(2.828427124746191*f[1]*wx1_sq+1.632993161855453*f[5]*dv1*wx1+0.2357022603955158*f[1]*dv1_sq); 
  out[2] += volFact*(2.828427124746191*f[0]*wx1*wx2+0.8164965809277261*f[2]*dv1*wx2+0.8164965809277261*f[3]*dv2*wx1+0.2357022603955158*f[7]*dv1*dv2); 
  out[3] += volFact*(2.828427124746191*f[1]*wx1*wx2+0.8164965809277261*f[5]*dv1*wx2+0.8164965809277261*f[6]*dv2*wx1+0.2357022603955158*f[11]*dv1*dv2); 
  out[4] += volFact*(2.828427124746191*f[0]*wx1*wx3+0.8164965809277261*f[2]*dv1*wx3+0.8164965809277261*f[4]*dv3*wx1+0.2357022603955158*f[9]*dv1*dv3); 
  out[5] += volFact*(2.828427124746191*f[1]*wx1*wx3+0.8164965809277261*f[5]*dv1*wx3+0.8164965809277261*f[8]*dv3*wx1+0.2357022603955158*f[12]*dv1*dv3); 
  out[6] += volFact*(2.828427124746191*f[0]*wx2_sq+1.632993161855453*f[3]*dv2*wx2+0.2357022603955158*f[0]*dv2_sq); 
  out[7] += volFact*(2.828427124746191*f[1]*wx2_sq+1.632993161855453*f[6]*dv2*wx2+0.2357022603955158*f[1]*dv2_sq); 
  out[8] += volFact*(2.828427124746191*f[0]*wx2*wx3+0.8164965809277261*f[3]*dv2*wx3+0.8164965809277261*f[4]*dv3*wx2+0.2357022603955158*f[10]*dv2*dv3); 
  out[9] += volFact*(2.828427124746191*f[1]*wx2*wx3+0.8164965809277261*f[6]*dv2*wx3+0.8164965809277261*f[8]*dv3*wx2+0.2357022603955158*f[13]*dv2*dv3); 
  out[10] += volFact*(2.828427124746191*f[0]*wx3_sq+1.632993161855453*f[4]*dv3*wx3+0.2357022603955158*f[0]*dv3_sq); 
  out[11] += volFact*(2.828427124746191*f[1]*wx3_sq+1.632993161855453*f[8]*dv3*wx3+0.2357022603955158*f[1]*dv3_sq); 
} 
GKYL_CU_DH void vlasov_M3i_1x3v_tensor_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx2_cu = wx2*wx2*wx2, dv2_cu = dv2*dv2*dv2; 
  const double wx3 = w[3], dv3 = dxv[3]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
  const double wx3_cu = wx3*wx3*wx3, dv3_cu = dv3*dv3*dv3; 
  out[0] += volFact*(2.828427124746191*f[0]*wx1*wx3_sq+0.8164965809277261*f[2]*dv1*wx3_sq+1.632993161855453*f[4]*dv3*wx1*wx3+0.4714045207910317*f[9]*dv1*dv3*wx3+2.828427124746191*f[0]*wx1*wx2_sq+0.8164965809277261*f[2]*dv1*wx2_sq+1.632993161855453*f[3]*dv2*wx1*wx2+0.4714045207910317*f[7]*dv1*dv2*wx2+2.828427124746191*f[0]*wx1*wx1_sq+2.449489742783178*f[2]*dv1*wx1_sq+0.2357022603955158*f[0]*dv3_sq*wx1+0.2357022603955158*f[0]*dv2_sq*wx1+0.7071067811865475*f[0]*dv1_sq*wx1+0.06804138174397717*f[2]*dv1*dv3_sq+0.06804138174397717*f[2]*dv1*dv2_sq+0.1224744871391589*f[2]*dv1*dv1_sq); 
  out[1] += volFact*(2.828427124746191*f[1]*wx1*wx3_sq+0.8164965809277261*f[5]*dv1*wx3_sq+1.632993161855453*f[8]*dv3*wx1*wx3+0.4714045207910317*f[12]*dv1*dv3*wx3+2.828427124746191*f[1]*wx1*wx2_sq+0.8164965809277261*f[5]*dv1*wx2_sq+1.632993161855453*f[6]*dv2*wx1*wx2+0.4714045207910317*f[11]*dv1*dv2*wx2+2.828427124746191*f[1]*wx1*wx1_sq+2.449489742783178*f[5]*dv1*wx1_sq+0.2357022603955158*f[1]*dv3_sq*wx1+0.2357022603955158*f[1]*dv2_sq*wx1+0.7071067811865475*f[1]*dv1_sq*wx1+0.06804138174397717*f[5]*dv1*dv3_sq+0.06804138174397717*f[5]*dv1*dv2_sq+0.1224744871391589*f[5]*dv1*dv1_sq); 
  out[2] += volFact*(2.828427124746191*f[0]*wx2*wx3_sq+0.8164965809277261*f[3]*dv2*wx3_sq+1.632993161855453*f[4]*dv3*wx2*wx3+0.4714045207910317*f[10]*dv2*dv3*wx3+2.828427124746191*f[0]*wx2*wx2_sq+2.449489742783178*f[3]*dv2*wx2_sq+2.828427124746191*f[0]*wx1_sq*wx2+1.632993161855453*f[2]*dv1*wx1*wx2+0.2357022603955158*f[0]*dv3_sq*wx2+0.7071067811865475*f[0]*dv2_sq*wx2+0.2357022603955158*f[0]*dv1_sq*wx2+0.8164965809277261*f[3]*dv2*wx1_sq+0.4714045207910317*f[7]*dv1*dv2*wx1+0.06804138174397717*f[3]*dv2*dv3_sq+0.1224744871391589*f[3]*dv2*dv2_sq+0.06804138174397717*f[3]*dv1_sq*dv2); 
  out[3] += volFact*(2.828427124746191*f[1]*wx2*wx3_sq+0.8164965809277261*f[6]*dv2*wx3_sq+1.632993161855453*f[8]*dv3*wx2*wx3+0.4714045207910317*f[13]*dv2*dv3*wx3+2.828427124746191*f[1]*wx2*wx2_sq+2.449489742783178*f[6]*dv2*wx2_sq+2.828427124746191*f[1]*wx1_sq*wx2+1.632993161855453*f[5]*dv1*wx1*wx2+0.2357022603955158*f[1]*dv3_sq*wx2+0.7071067811865475*f[1]*dv2_sq*wx2+0.2357022603955158*f[1]*dv1_sq*wx2+0.8164965809277261*f[6]*dv2*wx1_sq+0.4714045207910317*f[11]*dv1*dv2*wx1+0.06804138174397717*f[6]*dv2*dv3_sq+0.1224744871391589*f[6]*dv2*dv2_sq+0.06804138174397717*f[6]*dv1_sq*dv2); 
  out[4] += volFact*(2.828427124746191*f[0]*wx3*wx3_sq+2.449489742783178*f[4]*dv3*wx3_sq+2.828427124746191*f[0]*wx2_sq*wx3+1.632993161855453*f[3]*dv2*wx2*wx3+2.828427124746191*f[0]*wx1_sq*wx3+1.632993161855453*f[2]*dv1*wx1*wx3+0.7071067811865475*f[0]*dv3_sq*wx3+0.2357022603955158*f[0]*dv2_sq*wx3+0.2357022603955158*f[0]*dv1_sq*wx3+0.8164965809277261*f[4]*dv3*wx2_sq+0.4714045207910317*f[10]*dv2*dv3*wx2+0.8164965809277261*f[4]*dv3*wx1_sq+0.4714045207910317*f[9]*dv1*dv3*wx1+0.1224744871391589*f[4]*dv3*dv3_sq+0.06804138174397717*f[4]*dv2_sq*dv3+0.06804138174397717*f[4]*dv1_sq*dv3); 
  out[5] += volFact*(2.828427124746191*f[1]*wx3*wx3_sq+2.449489742783178*f[8]*dv3*wx3_sq+2.828427124746191*f[1]*wx2_sq*wx3+1.632993161855453*f[6]*dv2*wx2*wx3+2.828427124746191*f[1]*wx1_sq*wx3+1.632993161855453*f[5]*dv1*wx1*wx3+0.7071067811865475*f[1]*dv3_sq*wx3+0.2357022603955158*f[1]*dv2_sq*wx3+0.2357022603955158*f[1]*dv1_sq*wx3+0.8164965809277261*f[8]*dv3*wx2_sq+0.4714045207910317*f[13]*dv2*dv3*wx2+0.8164965809277261*f[8]*dv3*wx1_sq+0.4714045207910317*f[12]*dv1*dv3*wx1+0.1224744871391589*f[8]*dv3*dv3_sq+0.06804138174397717*f[8]*dv2_sq*dv3+0.06804138174397717*f[8]*dv1_sq*dv3); 
} 
GKYL_CU_DH void vlasov_M3ijk_1x3v_tensor_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx2_cu = wx2*wx2*wx2, dv2_cu = dv2*dv2*dv2; 
  const double wx3 = w[3], dv3 = dxv[3]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
  const double wx3_cu = wx3*wx3*wx3, dv3_cu = dv3*dv3*dv3; 
  out[0] += volFact*(2.828427124746191*f[0]*wx1*wx1_sq+2.449489742783178*f[2]*dv1*wx1_sq+0.7071067811865475*f[0]*dv1_sq*wx1+0.1224744871391589*f[2]*dv1*dv1_sq); 
  out[1] += volFact*(2.828427124746191*f[1]*wx1*wx1_sq+2.449489742783178*f[5]*dv1*wx1_sq+0.7071067811865475*f[1]*dv1_sq*wx1+0.1224744871391589*f[5]*dv1*dv1_sq); 
  out[2] += volFact*(2.828427124746191*f[0]*wx1_sq*wx2+1.632993161855453*f[2]*dv1*wx1*wx2+0.2357022603955158*f[0]*dv1_sq*wx2+0.8164965809277261*f[3]*dv2*wx1_sq+0.4714045207910317*f[7]*dv1*dv2*wx1+0.06804138174397717*f[3]*dv1_sq*dv2); 
  out[3] += volFact*(2.828427124746191*f[1]*wx1_sq*wx2+1.632993161855453*f[5]*dv1*wx1*wx2+0.2357022603955158*f[1]*dv1_sq*wx2+0.8164965809277261*f[6]*dv2*wx1_sq+0.4714045207910317*f[11]*dv1*dv2*wx1+0.06804138174397717*f[6]*dv1_sq*dv2); 
  out[4] += volFact*(2.828427124746191*f[0]*wx1_sq*wx3+1.632993161855453*f[2]*dv1*wx1*wx3+0.2357022603955158*f[0]*dv1_sq*wx3+0.8164965809277261*f[4]*dv3*wx1_sq+0.4714045207910317*f[9]*dv1*dv3*wx1+0.06804138174397717*f[4]*dv1_sq*dv3); 
  out[5] += volFact*(2.828427124746191*f[1]*wx1_sq*wx3+1.632993161855453*f[5]*dv1*wx1*wx3+0.2357022603955158*f[1]*dv1_sq*wx3+0.8164965809277261*f[8]*dv3*wx1_sq+0.4714045207910317*f[12]*dv1*dv3*wx1+0.06804138174397717*f[8]*dv1_sq*dv3); 
  out[6] += volFact*(2.828427124746191*f[0]*wx1*wx2_sq+0.8164965809277261*f[2]*dv1*wx2_sq+1.632993161855453*f[3]*dv2*wx1*wx2+0.4714045207910317*f[7]*dv1*dv2*wx2+0.2357022603955158*f[0]*dv2_sq*wx1+0.06804138174397717*f[2]*dv1*dv2_sq); 
  out[7] += volFact*(2.828427124746191*f[1]*wx1*wx2_sq+0.8164965809277261*f[5]*dv1*wx2_sq+1.632993161855453*f[6]*dv2*wx1*wx2+0.4714045207910317*f[11]*dv1*dv2*wx2+0.2357022603955158*f[1]*dv2_sq*wx1+0.06804138174397717*f[5]*dv1*dv2_sq); 
  out[8] += volFact*(2.828427124746191*f[0]*wx1*wx2*wx3+0.8164965809277261*f[2]*dv1*wx2*wx3+0.8164965809277261*f[3]*dv2*wx1*wx3+0.2357022603955158*f[7]*dv1*dv2*wx3+0.8164965809277261*f[4]*dv3*wx1*wx2+0.2357022603955158*f[9]*dv1*dv3*wx2+0.2357022603955158*f[10]*dv2*dv3*wx1+0.06804138174397717*f[14]*dv1*dv2*dv3); 
  out[9] += volFact*(2.828427124746191*f[1]*wx1*wx2*wx3+0.8164965809277261*f[5]*dv1*wx2*wx3+0.8164965809277261*f[6]*dv2*wx1*wx3+0.2357022603955158*f[11]*dv1*dv2*wx3+0.8164965809277261*f[8]*dv3*wx1*wx2+0.2357022603955158*f[12]*dv1*dv3*wx2+0.2357022603955158*f[13]*dv2*dv3*wx1+0.06804138174397717*f[15]*dv1*dv2*dv3); 
  out[10] += volFact*(2.828427124746191*f[0]*wx1*wx3_sq+0.8164965809277261*f[2]*dv1*wx3_sq+1.632993161855453*f[4]*dv3*wx1*wx3+0.4714045207910317*f[9]*dv1*dv3*wx3+0.2357022603955158*f[0]*dv3_sq*wx1+0.06804138174397717*f[2]*dv1*dv3_sq); 
  out[11] += volFact*(2.828427124746191*f[1]*wx1*wx3_sq+0.8164965809277261*f[5]*dv1*wx3_sq+1.632993161855453*f[8]*dv3*wx1*wx3+0.4714045207910317*f[12]*dv1*dv3*wx3+0.2357022603955158*f[1]*dv3_sq*wx1+0.06804138174397717*f[5]*dv1*dv3_sq); 
  out[12] += volFact*(2.828427124746191*f[0]*wx2*wx2_sq+2.449489742783178*f[3]*dv2*wx2_sq+0.7071067811865475*f[0]*dv2_sq*wx2+0.1224744871391589*f[3]*dv2*dv2_sq); 
  out[13] += volFact*(2.828427124746191*f[1]*wx2*wx2_sq+2.449489742783178*f[6]*dv2*wx2_sq+0.7071067811865475*f[1]*dv2_sq*wx2+0.1224744871391589*f[6]*dv2*dv2_sq); 
  out[14] += volFact*(2.828427124746191*f[0]*wx2_sq*wx3+1.632993161855453*f[3]*dv2*wx2*wx3+0.2357022603955158*f[0]*dv2_sq*wx3+0.8164965809277261*f[4]*dv3*wx2_sq+0.4714045207910317*f[10]*dv2*dv3*wx2+0.06804138174397717*f[4]*dv2_sq*dv3); 
  out[15] += volFact*(2.828427124746191*f[1]*wx2_sq*wx3+1.632993161855453*f[6]*dv2*wx2*wx3+0.2357022603955158*f[1]*dv2_sq*wx3+0.8164965809277261*f[8]*dv3*wx2_sq+0.4714045207910317*f[13]*dv2*dv3*wx2+0.06804138174397717*f[8]*dv2_sq*dv3); 
  out[16] += volFact*(2.828427124746191*f[0]*wx2*wx3_sq+0.8164965809277261*f[3]*dv2*wx3_sq+1.632993161855453*f[4]*dv3*wx2*wx3+0.4714045207910317*f[10]*dv2*dv3*wx3+0.2357022603955158*f[0]*dv3_sq*wx2+0.06804138174397717*f[3]*dv2*dv3_sq); 
  out[17] += volFact*(2.828427124746191*f[1]*wx2*wx3_sq+0.8164965809277261*f[6]*dv2*wx3_sq+1.632993161855453*f[8]*dv3*wx2*wx3+0.4714045207910317*f[13]*dv2*dv3*wx3+0.2357022603955158*f[1]*dv3_sq*wx2+0.06804138174397717*f[6]*dv2*dv3_sq); 
  out[18] += volFact*(2.828427124746191*f[0]*wx3*wx3_sq+2.449489742783178*f[4]*dv3*wx3_sq+0.7071067811865475*f[0]*dv3_sq*wx3+0.1224744871391589*f[4]*dv3*dv3_sq); 
  out[19] += volFact*(2.828427124746191*f[1]*wx3*wx3_sq+2.449489742783178*f[8]*dv3*wx3_sq+0.7071067811865475*f[1]*dv3_sq*wx3+0.1224744871391589*f[8]*dv3*dv3_sq); 
} 
GKYL_CU_DH void vlasov_five_moments_1x3v_tensor_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx3 = w[3], dv3 = dxv[3]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
  double tempM0[2], tempM1i[6]; 

  tempM0[0] = 2.828427124746191*f[0]*volFact; 
  tempM0[1] = 2.828427124746191*f[1]*volFact; 

  tempM1i[0] = tempM0[0]*wx1+0.8164965809277261*f[2]*dv1*volFact; 
  tempM1i[1] = tempM0[1]*wx1+0.8164965809277261*f[5]*dv1*volFact; 
  tempM1i[2] = tempM0[0]*wx2+0.8164965809277261*f[3]*dv2*volFact; 
  tempM1i[3] = tempM0[1]*wx2+0.8164965809277261*f[6]*dv2*volFact; 
  tempM1i[4] = tempM0[0]*wx3+0.8164965809277261*f[4]*dv3*volFact; 
  tempM1i[5] = tempM0[1]*wx3+0.8164965809277261*f[8]*dv3*volFact; 

  out[0] += tempM0[0]; 
  out[1] += tempM0[1]; 
  out[2] += tempM1i[0]; 
  out[3] += tempM1i[1]; 
  out[4] += tempM1i[2]; 
  out[5] += tempM1i[3]; 
  out[6] += tempM1i[4]; 
  out[7] += tempM1i[5]; 
  out[8] += tempM0[0]*((-1.0*wx3_sq)-1.0*wx2_sq-1.0*wx1_sq)+2.0*tempM1i[4]*wx3+2.0*tempM1i[2]*wx2+2.0*tempM1i[0]*wx1+(0.2357022603955158*f[0]*dv3_sq+0.2357022603955158*f[0]*dv2_sq+0.2357022603955158*f[0]*dv1_sq)*volFact; 
  out[9] += tempM0[1]*((-1.0*wx3_sq)-1.0*wx2_sq-1.0*wx1_sq)+2.0*tempM1i[5]*wx3+2.0*tempM1i[3]*wx2+2.0*tempM1i[1]*wx1+(0.2357022603955158*f[1]*dv3_sq+0.2357022603955158*f[1]*dv2_sq+0.2357022603955158*f[1]*dv1_sq)*volFact; 
} 
