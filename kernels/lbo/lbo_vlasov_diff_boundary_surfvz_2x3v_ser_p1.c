#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH void lbo_vlasov_diff_boundary_surfvz_2x3v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[5]:         Cell-center coordinates. 
  // dxv[5]:       Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[12]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[4]: sum of thermal speeds squared time their respective collisionalities. 
  // fSkin/Edge:    Distribution function in cells 
  // out:           Incremented distribution function in cell 
  double rdvSq4 = 4.0/(dxv[4]*dxv[4]); 

  double temp_diff[32] = {0.0}; 
  double temp_edge[32] = {0.0}; 
  double diff_incr[32] = {0.0}; 
  double edge_incr[32] = {0.0}; 
  double vol_incr[32] = {0.0}; 


  if (edge == -1) { 

  temp_diff[0] = (-0.5412658773652741*fSkin[5])-0.5412658773652741*fEdge[5]-0.5625*fSkin[0]+0.5625*fEdge[0]; 
  temp_diff[1] = (-0.5412658773652741*fSkin[12])-0.5412658773652741*fEdge[12]-0.5625*fSkin[1]+0.5625*fEdge[1]; 
  temp_diff[2] = (-0.5412658773652741*fSkin[13])-0.5412658773652741*fEdge[13]-0.5625*fSkin[2]+0.5625*fEdge[2]; 
  temp_diff[3] = (-0.5412658773652741*fSkin[14])-0.5412658773652741*fEdge[14]-0.5625*fSkin[3]+0.5625*fEdge[3]; 
  temp_diff[4] = (-0.5412658773652741*fSkin[15])-0.5412658773652741*fEdge[15]-0.5625*fSkin[4]+0.5625*fEdge[4]; 
  temp_diff[5] = (-1.4375*fSkin[5])-0.4375*fEdge[5]-1.407291281149713*fSkin[0]+0.5412658773652741*fEdge[0]; 
  temp_diff[6] = (-0.5412658773652741*fSkin[20])-0.5412658773652741*fEdge[20]-0.5625*fSkin[6]+0.5625*fEdge[6]; 
  temp_diff[7] = (-0.5412658773652741*fSkin[21])-0.5412658773652741*fEdge[21]-0.5625*fSkin[7]+0.5625*fEdge[7]; 
  temp_diff[8] = (-0.5412658773652741*fSkin[22])-0.5412658773652741*fEdge[22]-0.5625*fSkin[8]+0.5625*fEdge[8]; 
  temp_diff[9] = (-0.5412658773652741*fSkin[23])-0.5412658773652741*fEdge[23]-0.5625*fSkin[9]+0.5625*fEdge[9]; 
  temp_diff[10] = (-0.5412658773652741*fSkin[24])-0.5412658773652741*fEdge[24]-0.5625*fSkin[10]+0.5625*fEdge[10]; 
  temp_diff[11] = (-0.5412658773652741*fSkin[25])-0.5412658773652741*fEdge[25]-0.5625*fSkin[11]+0.5625*fEdge[11]; 
  temp_diff[12] = (-1.4375*fSkin[12])-0.4375*fEdge[12]-1.407291281149713*fSkin[1]+0.5412658773652741*fEdge[1]; 
  temp_diff[13] = (-1.4375*fSkin[13])-0.4375*fEdge[13]-1.407291281149713*fSkin[2]+0.5412658773652741*fEdge[2]; 
  temp_diff[14] = (-1.4375*fSkin[14])-0.4375*fEdge[14]-1.407291281149713*fSkin[3]+0.5412658773652741*fEdge[3]; 
  temp_diff[15] = (-1.4375*fSkin[15])-0.4375*fEdge[15]-1.407291281149713*fSkin[4]+0.5412658773652741*fEdge[4]; 
  temp_diff[16] = (-0.5412658773652741*fSkin[27])-0.5412658773652741*fEdge[27]-0.5625*fSkin[16]+0.5625*fEdge[16]; 
  temp_diff[17] = (-0.5412658773652741*fSkin[28])-0.5412658773652741*fEdge[28]-0.5625*fSkin[17]+0.5625*fEdge[17]; 
  temp_diff[18] = (-0.5412658773652741*fSkin[29])-0.5412658773652741*fEdge[29]-0.5625*fSkin[18]+0.5625*fEdge[18]; 
  temp_diff[19] = (-0.5412658773652741*fSkin[30])-0.5412658773652741*fEdge[30]-0.5625*fSkin[19]+0.5625*fEdge[19]; 
  temp_diff[20] = (-1.4375*fSkin[20])-0.4375*fEdge[20]-1.407291281149713*fSkin[6]+0.5412658773652741*fEdge[6]; 
  temp_diff[21] = (-1.4375*fSkin[21])-0.4375*fEdge[21]-1.407291281149713*fSkin[7]+0.5412658773652741*fEdge[7]; 
  temp_diff[22] = (-1.4375*fSkin[22])-0.4375*fEdge[22]-1.407291281149713*fSkin[8]+0.5412658773652741*fEdge[8]; 
  temp_diff[23] = (-1.4375*fSkin[23])-0.4375*fEdge[23]-1.407291281149713*fSkin[9]+0.5412658773652741*fEdge[9]; 
  temp_diff[24] = (-1.4375*fSkin[24])-0.4375*fEdge[24]-1.407291281149713*fSkin[10]+0.5412658773652741*fEdge[10]; 
  temp_diff[25] = (-1.4375*fSkin[25])-0.4375*fEdge[25]-1.407291281149713*fSkin[11]+0.5412658773652741*fEdge[11]; 
  temp_diff[26] = (-0.5412658773652741*fSkin[31])-0.5412658773652741*fEdge[31]-0.5625*fSkin[26]+0.5625*fEdge[26]; 
  temp_diff[27] = (-1.4375*fSkin[27])-0.4375*fEdge[27]-1.407291281149713*fSkin[16]+0.5412658773652741*fEdge[16]; 
  temp_diff[28] = (-1.4375*fSkin[28])-0.4375*fEdge[28]-1.407291281149713*fSkin[17]+0.5412658773652741*fEdge[17]; 
  temp_diff[29] = (-1.4375*fSkin[29])-0.4375*fEdge[29]-1.407291281149713*fSkin[18]+0.5412658773652741*fEdge[18]; 
  temp_diff[30] = (-1.4375*fSkin[30])-0.4375*fEdge[30]-1.407291281149713*fSkin[19]+0.5412658773652741*fEdge[19]; 
  temp_diff[31] = (-1.4375*fSkin[31])-0.4375*fEdge[31]-1.407291281149713*fSkin[26]+0.5412658773652741*fEdge[26]; 

  temp_edge[5] = 0.8660254037844386*fSkin[0]-1.5*fSkin[5]; 
  temp_edge[12] = 0.8660254037844386*fSkin[1]-1.5*fSkin[12]; 
  temp_edge[13] = 0.8660254037844386*fSkin[2]-1.5*fSkin[13]; 
  temp_edge[14] = 0.8660254037844386*fSkin[3]-1.5*fSkin[14]; 
  temp_edge[15] = 0.8660254037844386*fSkin[4]-1.5*fSkin[15]; 
  temp_edge[20] = 0.8660254037844386*fSkin[6]-1.5*fSkin[20]; 
  temp_edge[21] = 0.8660254037844386*fSkin[7]-1.5*fSkin[21]; 
  temp_edge[22] = 0.8660254037844386*fSkin[8]-1.5*fSkin[22]; 
  temp_edge[23] = 0.8660254037844386*fSkin[9]-1.5*fSkin[23]; 
  temp_edge[24] = 0.8660254037844386*fSkin[10]-1.5*fSkin[24]; 
  temp_edge[25] = 0.8660254037844386*fSkin[11]-1.5*fSkin[25]; 
  temp_edge[27] = 0.8660254037844386*fSkin[16]-1.5*fSkin[27]; 
  temp_edge[28] = 0.8660254037844386*fSkin[17]-1.5*fSkin[28]; 
  temp_edge[29] = 0.8660254037844386*fSkin[18]-1.5*fSkin[29]; 
  temp_edge[30] = 0.8660254037844386*fSkin[19]-1.5*fSkin[30]; 
  temp_edge[31] = 0.8660254037844386*fSkin[26]-1.5*fSkin[31]; 

  diff_incr[0] = 0.5*nuVtSqSum[3]*temp_diff[6]+0.5*nuVtSqSum[2]*temp_diff[2]+0.5*nuVtSqSum[1]*temp_diff[1]+0.5*nuVtSqSum[0]*temp_diff[0]; 
  diff_incr[1] = 0.5*nuVtSqSum[2]*temp_diff[6]+0.5*temp_diff[2]*nuVtSqSum[3]+0.5*nuVtSqSum[0]*temp_diff[1]+0.5*temp_diff[0]*nuVtSqSum[1]; 
  diff_incr[2] = 0.5*nuVtSqSum[1]*temp_diff[6]+0.5*temp_diff[1]*nuVtSqSum[3]+0.5*nuVtSqSum[0]*temp_diff[2]+0.5*temp_diff[0]*nuVtSqSum[2]; 
  diff_incr[3] = 0.5*nuVtSqSum[3]*temp_diff[16]+0.5*nuVtSqSum[2]*temp_diff[8]+0.5*nuVtSqSum[1]*temp_diff[7]+0.5*nuVtSqSum[0]*temp_diff[3]; 
  diff_incr[4] = 0.5*nuVtSqSum[3]*temp_diff[17]+0.5*nuVtSqSum[2]*temp_diff[10]+0.5*nuVtSqSum[1]*temp_diff[9]+0.5*nuVtSqSum[0]*temp_diff[4]; 
  diff_incr[5] = 0.5*nuVtSqSum[3]*temp_diff[20]+0.5*nuVtSqSum[2]*temp_diff[13]+0.5*nuVtSqSum[1]*temp_diff[12]+0.5*nuVtSqSum[0]*temp_diff[5]; 
  diff_incr[6] = 0.5*nuVtSqSum[0]*temp_diff[6]+0.5*temp_diff[0]*nuVtSqSum[3]+0.5*nuVtSqSum[1]*temp_diff[2]+0.5*temp_diff[1]*nuVtSqSum[2]; 
  diff_incr[7] = 0.5*nuVtSqSum[2]*temp_diff[16]+0.5*nuVtSqSum[3]*temp_diff[8]+0.5*nuVtSqSum[0]*temp_diff[7]+0.5*nuVtSqSum[1]*temp_diff[3]; 
  diff_incr[8] = 0.5*nuVtSqSum[1]*temp_diff[16]+0.5*nuVtSqSum[0]*temp_diff[8]+0.5*nuVtSqSum[3]*temp_diff[7]+0.5*nuVtSqSum[2]*temp_diff[3]; 
  diff_incr[9] = 0.5*nuVtSqSum[2]*temp_diff[17]+0.5*nuVtSqSum[3]*temp_diff[10]+0.5*nuVtSqSum[0]*temp_diff[9]+0.5*nuVtSqSum[1]*temp_diff[4]; 
  diff_incr[10] = 0.5*nuVtSqSum[1]*temp_diff[17]+0.5*nuVtSqSum[0]*temp_diff[10]+0.5*nuVtSqSum[3]*temp_diff[9]+0.5*nuVtSqSum[2]*temp_diff[4]; 
  diff_incr[11] = 0.5*nuVtSqSum[3]*temp_diff[26]+0.5*nuVtSqSum[2]*temp_diff[19]+0.5*nuVtSqSum[1]*temp_diff[18]+0.5*nuVtSqSum[0]*temp_diff[11]; 
  diff_incr[12] = 0.5*nuVtSqSum[2]*temp_diff[20]+0.5*nuVtSqSum[3]*temp_diff[13]+0.5*nuVtSqSum[0]*temp_diff[12]+0.5*nuVtSqSum[1]*temp_diff[5]; 
  diff_incr[13] = 0.5*nuVtSqSum[1]*temp_diff[20]+0.5*nuVtSqSum[0]*temp_diff[13]+0.5*nuVtSqSum[3]*temp_diff[12]+0.5*nuVtSqSum[2]*temp_diff[5]; 
  diff_incr[14] = 0.5*nuVtSqSum[3]*temp_diff[27]+0.5*nuVtSqSum[2]*temp_diff[22]+0.5*nuVtSqSum[1]*temp_diff[21]+0.5*nuVtSqSum[0]*temp_diff[14]; 
  diff_incr[15] = 0.5*nuVtSqSum[3]*temp_diff[28]+0.5*nuVtSqSum[2]*temp_diff[24]+0.5*nuVtSqSum[1]*temp_diff[23]+0.5*nuVtSqSum[0]*temp_diff[15]; 
  diff_incr[16] = 0.5*nuVtSqSum[0]*temp_diff[16]+0.5*nuVtSqSum[1]*temp_diff[8]+0.5*nuVtSqSum[2]*temp_diff[7]+0.5*nuVtSqSum[3]*temp_diff[3]; 
  diff_incr[17] = 0.5*nuVtSqSum[0]*temp_diff[17]+0.5*nuVtSqSum[1]*temp_diff[10]+0.5*nuVtSqSum[2]*temp_diff[9]+0.5*nuVtSqSum[3]*temp_diff[4]; 
  diff_incr[18] = 0.5*nuVtSqSum[2]*temp_diff[26]+0.5*nuVtSqSum[3]*temp_diff[19]+0.5*nuVtSqSum[0]*temp_diff[18]+0.5*nuVtSqSum[1]*temp_diff[11]; 
  diff_incr[19] = 0.5*nuVtSqSum[1]*temp_diff[26]+0.5*nuVtSqSum[0]*temp_diff[19]+0.5*nuVtSqSum[3]*temp_diff[18]+0.5*nuVtSqSum[2]*temp_diff[11]; 
  diff_incr[20] = 0.5*nuVtSqSum[0]*temp_diff[20]+0.5*nuVtSqSum[1]*temp_diff[13]+0.5*nuVtSqSum[2]*temp_diff[12]+0.5*nuVtSqSum[3]*temp_diff[5]; 
  diff_incr[21] = 0.5*nuVtSqSum[2]*temp_diff[27]+0.5*nuVtSqSum[3]*temp_diff[22]+0.5*nuVtSqSum[0]*temp_diff[21]+0.5*nuVtSqSum[1]*temp_diff[14]; 
  diff_incr[22] = 0.5*nuVtSqSum[1]*temp_diff[27]+0.5*nuVtSqSum[0]*temp_diff[22]+0.5*nuVtSqSum[3]*temp_diff[21]+0.5*nuVtSqSum[2]*temp_diff[14]; 
  diff_incr[23] = 0.5*nuVtSqSum[2]*temp_diff[28]+0.5*nuVtSqSum[3]*temp_diff[24]+0.5*nuVtSqSum[0]*temp_diff[23]+0.5*nuVtSqSum[1]*temp_diff[15]; 
  diff_incr[24] = 0.5*nuVtSqSum[1]*temp_diff[28]+0.5*nuVtSqSum[0]*temp_diff[24]+0.5*nuVtSqSum[3]*temp_diff[23]+0.5*nuVtSqSum[2]*temp_diff[15]; 
  diff_incr[25] = 0.5*nuVtSqSum[3]*temp_diff[31]+0.5*nuVtSqSum[2]*temp_diff[30]+0.5*nuVtSqSum[1]*temp_diff[29]+0.5*nuVtSqSum[0]*temp_diff[25]; 
  diff_incr[26] = 0.5*nuVtSqSum[0]*temp_diff[26]+0.5*nuVtSqSum[1]*temp_diff[19]+0.5*nuVtSqSum[2]*temp_diff[18]+0.5*nuVtSqSum[3]*temp_diff[11]; 
  diff_incr[27] = 0.5*nuVtSqSum[0]*temp_diff[27]+0.5*nuVtSqSum[1]*temp_diff[22]+0.5*nuVtSqSum[2]*temp_diff[21]+0.5*nuVtSqSum[3]*temp_diff[14]; 
  diff_incr[28] = 0.5*nuVtSqSum[0]*temp_diff[28]+0.5*nuVtSqSum[1]*temp_diff[24]+0.5*nuVtSqSum[2]*temp_diff[23]+0.5*nuVtSqSum[3]*temp_diff[15]; 
  diff_incr[29] = 0.5*nuVtSqSum[2]*temp_diff[31]+0.5*nuVtSqSum[3]*temp_diff[30]+0.5*nuVtSqSum[0]*temp_diff[29]+0.5*nuVtSqSum[1]*temp_diff[25]; 
  diff_incr[30] = 0.5*nuVtSqSum[1]*temp_diff[31]+0.5*nuVtSqSum[0]*temp_diff[30]+0.5*nuVtSqSum[3]*temp_diff[29]+0.5*nuVtSqSum[2]*temp_diff[25]; 
  diff_incr[31] = 0.5*nuVtSqSum[0]*temp_diff[31]+0.5*nuVtSqSum[1]*temp_diff[30]+0.5*nuVtSqSum[2]*temp_diff[29]+0.5*nuVtSqSum[3]*temp_diff[25]; 

  edge_incr[0] = 0.5*nuVtSqSum[3]*temp_edge[6]+0.5*nuVtSqSum[2]*temp_edge[2]+0.5*nuVtSqSum[1]*temp_edge[1]+0.5*nuVtSqSum[0]*temp_edge[0]; 
  edge_incr[1] = 0.5*nuVtSqSum[2]*temp_edge[6]+0.5*temp_edge[2]*nuVtSqSum[3]+0.5*nuVtSqSum[0]*temp_edge[1]+0.5*temp_edge[0]*nuVtSqSum[1]; 
  edge_incr[2] = 0.5*nuVtSqSum[1]*temp_edge[6]+0.5*temp_edge[1]*nuVtSqSum[3]+0.5*nuVtSqSum[0]*temp_edge[2]+0.5*temp_edge[0]*nuVtSqSum[2]; 
  edge_incr[3] = 0.5*nuVtSqSum[3]*temp_edge[16]+0.5*nuVtSqSum[2]*temp_edge[8]+0.5*nuVtSqSum[1]*temp_edge[7]+0.5*nuVtSqSum[0]*temp_edge[3]; 
  edge_incr[4] = 0.5*nuVtSqSum[3]*temp_edge[17]+0.5*nuVtSqSum[2]*temp_edge[10]+0.5*nuVtSqSum[1]*temp_edge[9]+0.5*nuVtSqSum[0]*temp_edge[4]; 
  edge_incr[5] = 0.5*nuVtSqSum[3]*temp_edge[20]+0.5*nuVtSqSum[2]*temp_edge[13]+0.5*nuVtSqSum[1]*temp_edge[12]+0.5*nuVtSqSum[0]*temp_edge[5]; 
  edge_incr[6] = 0.5*nuVtSqSum[0]*temp_edge[6]+0.5*temp_edge[0]*nuVtSqSum[3]+0.5*nuVtSqSum[1]*temp_edge[2]+0.5*temp_edge[1]*nuVtSqSum[2]; 
  edge_incr[7] = 0.5*nuVtSqSum[2]*temp_edge[16]+0.5*nuVtSqSum[3]*temp_edge[8]+0.5*nuVtSqSum[0]*temp_edge[7]+0.5*nuVtSqSum[1]*temp_edge[3]; 
  edge_incr[8] = 0.5*nuVtSqSum[1]*temp_edge[16]+0.5*nuVtSqSum[0]*temp_edge[8]+0.5*nuVtSqSum[3]*temp_edge[7]+0.5*nuVtSqSum[2]*temp_edge[3]; 
  edge_incr[9] = 0.5*nuVtSqSum[2]*temp_edge[17]+0.5*nuVtSqSum[3]*temp_edge[10]+0.5*nuVtSqSum[0]*temp_edge[9]+0.5*nuVtSqSum[1]*temp_edge[4]; 
  edge_incr[10] = 0.5*nuVtSqSum[1]*temp_edge[17]+0.5*nuVtSqSum[0]*temp_edge[10]+0.5*nuVtSqSum[3]*temp_edge[9]+0.5*nuVtSqSum[2]*temp_edge[4]; 
  edge_incr[11] = 0.5*nuVtSqSum[3]*temp_edge[26]+0.5*nuVtSqSum[2]*temp_edge[19]+0.5*nuVtSqSum[1]*temp_edge[18]+0.5*nuVtSqSum[0]*temp_edge[11]; 
  edge_incr[12] = 0.5*nuVtSqSum[2]*temp_edge[20]+0.5*nuVtSqSum[3]*temp_edge[13]+0.5*nuVtSqSum[0]*temp_edge[12]+0.5*nuVtSqSum[1]*temp_edge[5]; 
  edge_incr[13] = 0.5*nuVtSqSum[1]*temp_edge[20]+0.5*nuVtSqSum[0]*temp_edge[13]+0.5*nuVtSqSum[3]*temp_edge[12]+0.5*nuVtSqSum[2]*temp_edge[5]; 
  edge_incr[14] = 0.5*nuVtSqSum[3]*temp_edge[27]+0.5*nuVtSqSum[2]*temp_edge[22]+0.5*nuVtSqSum[1]*temp_edge[21]+0.5*nuVtSqSum[0]*temp_edge[14]; 
  edge_incr[15] = 0.5*nuVtSqSum[3]*temp_edge[28]+0.5*nuVtSqSum[2]*temp_edge[24]+0.5*nuVtSqSum[1]*temp_edge[23]+0.5*nuVtSqSum[0]*temp_edge[15]; 
  edge_incr[16] = 0.5*nuVtSqSum[0]*temp_edge[16]+0.5*nuVtSqSum[1]*temp_edge[8]+0.5*nuVtSqSum[2]*temp_edge[7]+0.5*nuVtSqSum[3]*temp_edge[3]; 
  edge_incr[17] = 0.5*nuVtSqSum[0]*temp_edge[17]+0.5*nuVtSqSum[1]*temp_edge[10]+0.5*nuVtSqSum[2]*temp_edge[9]+0.5*nuVtSqSum[3]*temp_edge[4]; 
  edge_incr[18] = 0.5*nuVtSqSum[2]*temp_edge[26]+0.5*nuVtSqSum[3]*temp_edge[19]+0.5*nuVtSqSum[0]*temp_edge[18]+0.5*nuVtSqSum[1]*temp_edge[11]; 
  edge_incr[19] = 0.5*nuVtSqSum[1]*temp_edge[26]+0.5*nuVtSqSum[0]*temp_edge[19]+0.5*nuVtSqSum[3]*temp_edge[18]+0.5*nuVtSqSum[2]*temp_edge[11]; 
  edge_incr[20] = 0.5*nuVtSqSum[0]*temp_edge[20]+0.5*nuVtSqSum[1]*temp_edge[13]+0.5*nuVtSqSum[2]*temp_edge[12]+0.5*nuVtSqSum[3]*temp_edge[5]; 
  edge_incr[21] = 0.5*nuVtSqSum[2]*temp_edge[27]+0.5*nuVtSqSum[3]*temp_edge[22]+0.5*nuVtSqSum[0]*temp_edge[21]+0.5*nuVtSqSum[1]*temp_edge[14]; 
  edge_incr[22] = 0.5*nuVtSqSum[1]*temp_edge[27]+0.5*nuVtSqSum[0]*temp_edge[22]+0.5*nuVtSqSum[3]*temp_edge[21]+0.5*nuVtSqSum[2]*temp_edge[14]; 
  edge_incr[23] = 0.5*nuVtSqSum[2]*temp_edge[28]+0.5*nuVtSqSum[3]*temp_edge[24]+0.5*nuVtSqSum[0]*temp_edge[23]+0.5*nuVtSqSum[1]*temp_edge[15]; 
  edge_incr[24] = 0.5*nuVtSqSum[1]*temp_edge[28]+0.5*nuVtSqSum[0]*temp_edge[24]+0.5*nuVtSqSum[3]*temp_edge[23]+0.5*nuVtSqSum[2]*temp_edge[15]; 
  edge_incr[25] = 0.5*nuVtSqSum[3]*temp_edge[31]+0.5*nuVtSqSum[2]*temp_edge[30]+0.5*nuVtSqSum[1]*temp_edge[29]+0.5*nuVtSqSum[0]*temp_edge[25]; 
  edge_incr[26] = 0.5*nuVtSqSum[0]*temp_edge[26]+0.5*nuVtSqSum[1]*temp_edge[19]+0.5*nuVtSqSum[2]*temp_edge[18]+0.5*nuVtSqSum[3]*temp_edge[11]; 
  edge_incr[27] = 0.5*nuVtSqSum[0]*temp_edge[27]+0.5*nuVtSqSum[1]*temp_edge[22]+0.5*nuVtSqSum[2]*temp_edge[21]+0.5*nuVtSqSum[3]*temp_edge[14]; 
  edge_incr[28] = 0.5*nuVtSqSum[0]*temp_edge[28]+0.5*nuVtSqSum[1]*temp_edge[24]+0.5*nuVtSqSum[2]*temp_edge[23]+0.5*nuVtSqSum[3]*temp_edge[15]; 
  edge_incr[29] = 0.5*nuVtSqSum[2]*temp_edge[31]+0.5*nuVtSqSum[3]*temp_edge[30]+0.5*nuVtSqSum[0]*temp_edge[29]+0.5*nuVtSqSum[1]*temp_edge[25]; 
  edge_incr[30] = 0.5*nuVtSqSum[1]*temp_edge[31]+0.5*nuVtSqSum[0]*temp_edge[30]+0.5*nuVtSqSum[3]*temp_edge[29]+0.5*nuVtSqSum[2]*temp_edge[25]; 
  edge_incr[31] = 0.5*nuVtSqSum[0]*temp_edge[31]+0.5*nuVtSqSum[1]*temp_edge[30]+0.5*nuVtSqSum[2]*temp_edge[29]+0.5*nuVtSqSum[3]*temp_edge[25]; 


  } else { 

  temp_diff[0] = 0.5412658773652741*fSkin[5]+0.5412658773652741*fEdge[5]-0.5625*fSkin[0]+0.5625*fEdge[0]; 
  temp_diff[1] = 0.5412658773652741*fSkin[12]+0.5412658773652741*fEdge[12]-0.5625*fSkin[1]+0.5625*fEdge[1]; 
  temp_diff[2] = 0.5412658773652741*fSkin[13]+0.5412658773652741*fEdge[13]-0.5625*fSkin[2]+0.5625*fEdge[2]; 
  temp_diff[3] = 0.5412658773652741*fSkin[14]+0.5412658773652741*fEdge[14]-0.5625*fSkin[3]+0.5625*fEdge[3]; 
  temp_diff[4] = 0.5412658773652741*fSkin[15]+0.5412658773652741*fEdge[15]-0.5625*fSkin[4]+0.5625*fEdge[4]; 
  temp_diff[5] = (-1.4375*fSkin[5])-0.4375*fEdge[5]+1.407291281149713*fSkin[0]-0.5412658773652741*fEdge[0]; 
  temp_diff[6] = 0.5412658773652741*fSkin[20]+0.5412658773652741*fEdge[20]-0.5625*fSkin[6]+0.5625*fEdge[6]; 
  temp_diff[7] = 0.5412658773652741*fSkin[21]+0.5412658773652741*fEdge[21]-0.5625*fSkin[7]+0.5625*fEdge[7]; 
  temp_diff[8] = 0.5412658773652741*fSkin[22]+0.5412658773652741*fEdge[22]-0.5625*fSkin[8]+0.5625*fEdge[8]; 
  temp_diff[9] = 0.5412658773652741*fSkin[23]+0.5412658773652741*fEdge[23]-0.5625*fSkin[9]+0.5625*fEdge[9]; 
  temp_diff[10] = 0.5412658773652741*fSkin[24]+0.5412658773652741*fEdge[24]-0.5625*fSkin[10]+0.5625*fEdge[10]; 
  temp_diff[11] = 0.5412658773652741*fSkin[25]+0.5412658773652741*fEdge[25]-0.5625*fSkin[11]+0.5625*fEdge[11]; 
  temp_diff[12] = (-1.4375*fSkin[12])-0.4375*fEdge[12]+1.407291281149713*fSkin[1]-0.5412658773652741*fEdge[1]; 
  temp_diff[13] = (-1.4375*fSkin[13])-0.4375*fEdge[13]+1.407291281149713*fSkin[2]-0.5412658773652741*fEdge[2]; 
  temp_diff[14] = (-1.4375*fSkin[14])-0.4375*fEdge[14]+1.407291281149713*fSkin[3]-0.5412658773652741*fEdge[3]; 
  temp_diff[15] = (-1.4375*fSkin[15])-0.4375*fEdge[15]+1.407291281149713*fSkin[4]-0.5412658773652741*fEdge[4]; 
  temp_diff[16] = 0.5412658773652741*fSkin[27]+0.5412658773652741*fEdge[27]-0.5625*fSkin[16]+0.5625*fEdge[16]; 
  temp_diff[17] = 0.5412658773652741*fSkin[28]+0.5412658773652741*fEdge[28]-0.5625*fSkin[17]+0.5625*fEdge[17]; 
  temp_diff[18] = 0.5412658773652741*fSkin[29]+0.5412658773652741*fEdge[29]-0.5625*fSkin[18]+0.5625*fEdge[18]; 
  temp_diff[19] = 0.5412658773652741*fSkin[30]+0.5412658773652741*fEdge[30]-0.5625*fSkin[19]+0.5625*fEdge[19]; 
  temp_diff[20] = (-1.4375*fSkin[20])-0.4375*fEdge[20]+1.407291281149713*fSkin[6]-0.5412658773652741*fEdge[6]; 
  temp_diff[21] = (-1.4375*fSkin[21])-0.4375*fEdge[21]+1.407291281149713*fSkin[7]-0.5412658773652741*fEdge[7]; 
  temp_diff[22] = (-1.4375*fSkin[22])-0.4375*fEdge[22]+1.407291281149713*fSkin[8]-0.5412658773652741*fEdge[8]; 
  temp_diff[23] = (-1.4375*fSkin[23])-0.4375*fEdge[23]+1.407291281149713*fSkin[9]-0.5412658773652741*fEdge[9]; 
  temp_diff[24] = (-1.4375*fSkin[24])-0.4375*fEdge[24]+1.407291281149713*fSkin[10]-0.5412658773652741*fEdge[10]; 
  temp_diff[25] = (-1.4375*fSkin[25])-0.4375*fEdge[25]+1.407291281149713*fSkin[11]-0.5412658773652741*fEdge[11]; 
  temp_diff[26] = 0.5412658773652741*fSkin[31]+0.5412658773652741*fEdge[31]-0.5625*fSkin[26]+0.5625*fEdge[26]; 
  temp_diff[27] = (-1.4375*fSkin[27])-0.4375*fEdge[27]+1.407291281149713*fSkin[16]-0.5412658773652741*fEdge[16]; 
  temp_diff[28] = (-1.4375*fSkin[28])-0.4375*fEdge[28]+1.407291281149713*fSkin[17]-0.5412658773652741*fEdge[17]; 
  temp_diff[29] = (-1.4375*fSkin[29])-0.4375*fEdge[29]+1.407291281149713*fSkin[18]-0.5412658773652741*fEdge[18]; 
  temp_diff[30] = (-1.4375*fSkin[30])-0.4375*fEdge[30]+1.407291281149713*fSkin[19]-0.5412658773652741*fEdge[19]; 
  temp_diff[31] = (-1.4375*fSkin[31])-0.4375*fEdge[31]+1.407291281149713*fSkin[26]-0.5412658773652741*fEdge[26]; 

  temp_edge[5] = (-1.5*fSkin[5])-0.8660254037844386*fSkin[0]; 
  temp_edge[12] = (-1.5*fSkin[12])-0.8660254037844386*fSkin[1]; 
  temp_edge[13] = (-1.5*fSkin[13])-0.8660254037844386*fSkin[2]; 
  temp_edge[14] = (-1.5*fSkin[14])-0.8660254037844386*fSkin[3]; 
  temp_edge[15] = (-1.5*fSkin[15])-0.8660254037844386*fSkin[4]; 
  temp_edge[20] = (-1.5*fSkin[20])-0.8660254037844386*fSkin[6]; 
  temp_edge[21] = (-1.5*fSkin[21])-0.8660254037844386*fSkin[7]; 
  temp_edge[22] = (-1.5*fSkin[22])-0.8660254037844386*fSkin[8]; 
  temp_edge[23] = (-1.5*fSkin[23])-0.8660254037844386*fSkin[9]; 
  temp_edge[24] = (-1.5*fSkin[24])-0.8660254037844386*fSkin[10]; 
  temp_edge[25] = (-1.5*fSkin[25])-0.8660254037844386*fSkin[11]; 
  temp_edge[27] = (-1.5*fSkin[27])-0.8660254037844386*fSkin[16]; 
  temp_edge[28] = (-1.5*fSkin[28])-0.8660254037844386*fSkin[17]; 
  temp_edge[29] = (-1.5*fSkin[29])-0.8660254037844386*fSkin[18]; 
  temp_edge[30] = (-1.5*fSkin[30])-0.8660254037844386*fSkin[19]; 
  temp_edge[31] = (-1.5*fSkin[31])-0.8660254037844386*fSkin[26]; 

  diff_incr[0] = 0.5*nuVtSqSum[3]*temp_diff[6]+0.5*nuVtSqSum[2]*temp_diff[2]+0.5*nuVtSqSum[1]*temp_diff[1]+0.5*nuVtSqSum[0]*temp_diff[0]; 
  diff_incr[1] = 0.5*nuVtSqSum[2]*temp_diff[6]+0.5*temp_diff[2]*nuVtSqSum[3]+0.5*nuVtSqSum[0]*temp_diff[1]+0.5*temp_diff[0]*nuVtSqSum[1]; 
  diff_incr[2] = 0.5*nuVtSqSum[1]*temp_diff[6]+0.5*temp_diff[1]*nuVtSqSum[3]+0.5*nuVtSqSum[0]*temp_diff[2]+0.5*temp_diff[0]*nuVtSqSum[2]; 
  diff_incr[3] = 0.5*nuVtSqSum[3]*temp_diff[16]+0.5*nuVtSqSum[2]*temp_diff[8]+0.5*nuVtSqSum[1]*temp_diff[7]+0.5*nuVtSqSum[0]*temp_diff[3]; 
  diff_incr[4] = 0.5*nuVtSqSum[3]*temp_diff[17]+0.5*nuVtSqSum[2]*temp_diff[10]+0.5*nuVtSqSum[1]*temp_diff[9]+0.5*nuVtSqSum[0]*temp_diff[4]; 
  diff_incr[5] = 0.5*nuVtSqSum[3]*temp_diff[20]+0.5*nuVtSqSum[2]*temp_diff[13]+0.5*nuVtSqSum[1]*temp_diff[12]+0.5*nuVtSqSum[0]*temp_diff[5]; 
  diff_incr[6] = 0.5*nuVtSqSum[0]*temp_diff[6]+0.5*temp_diff[0]*nuVtSqSum[3]+0.5*nuVtSqSum[1]*temp_diff[2]+0.5*temp_diff[1]*nuVtSqSum[2]; 
  diff_incr[7] = 0.5*nuVtSqSum[2]*temp_diff[16]+0.5*nuVtSqSum[3]*temp_diff[8]+0.5*nuVtSqSum[0]*temp_diff[7]+0.5*nuVtSqSum[1]*temp_diff[3]; 
  diff_incr[8] = 0.5*nuVtSqSum[1]*temp_diff[16]+0.5*nuVtSqSum[0]*temp_diff[8]+0.5*nuVtSqSum[3]*temp_diff[7]+0.5*nuVtSqSum[2]*temp_diff[3]; 
  diff_incr[9] = 0.5*nuVtSqSum[2]*temp_diff[17]+0.5*nuVtSqSum[3]*temp_diff[10]+0.5*nuVtSqSum[0]*temp_diff[9]+0.5*nuVtSqSum[1]*temp_diff[4]; 
  diff_incr[10] = 0.5*nuVtSqSum[1]*temp_diff[17]+0.5*nuVtSqSum[0]*temp_diff[10]+0.5*nuVtSqSum[3]*temp_diff[9]+0.5*nuVtSqSum[2]*temp_diff[4]; 
  diff_incr[11] = 0.5*nuVtSqSum[3]*temp_diff[26]+0.5*nuVtSqSum[2]*temp_diff[19]+0.5*nuVtSqSum[1]*temp_diff[18]+0.5*nuVtSqSum[0]*temp_diff[11]; 
  diff_incr[12] = 0.5*nuVtSqSum[2]*temp_diff[20]+0.5*nuVtSqSum[3]*temp_diff[13]+0.5*nuVtSqSum[0]*temp_diff[12]+0.5*nuVtSqSum[1]*temp_diff[5]; 
  diff_incr[13] = 0.5*nuVtSqSum[1]*temp_diff[20]+0.5*nuVtSqSum[0]*temp_diff[13]+0.5*nuVtSqSum[3]*temp_diff[12]+0.5*nuVtSqSum[2]*temp_diff[5]; 
  diff_incr[14] = 0.5*nuVtSqSum[3]*temp_diff[27]+0.5*nuVtSqSum[2]*temp_diff[22]+0.5*nuVtSqSum[1]*temp_diff[21]+0.5*nuVtSqSum[0]*temp_diff[14]; 
  diff_incr[15] = 0.5*nuVtSqSum[3]*temp_diff[28]+0.5*nuVtSqSum[2]*temp_diff[24]+0.5*nuVtSqSum[1]*temp_diff[23]+0.5*nuVtSqSum[0]*temp_diff[15]; 
  diff_incr[16] = 0.5*nuVtSqSum[0]*temp_diff[16]+0.5*nuVtSqSum[1]*temp_diff[8]+0.5*nuVtSqSum[2]*temp_diff[7]+0.5*nuVtSqSum[3]*temp_diff[3]; 
  diff_incr[17] = 0.5*nuVtSqSum[0]*temp_diff[17]+0.5*nuVtSqSum[1]*temp_diff[10]+0.5*nuVtSqSum[2]*temp_diff[9]+0.5*nuVtSqSum[3]*temp_diff[4]; 
  diff_incr[18] = 0.5*nuVtSqSum[2]*temp_diff[26]+0.5*nuVtSqSum[3]*temp_diff[19]+0.5*nuVtSqSum[0]*temp_diff[18]+0.5*nuVtSqSum[1]*temp_diff[11]; 
  diff_incr[19] = 0.5*nuVtSqSum[1]*temp_diff[26]+0.5*nuVtSqSum[0]*temp_diff[19]+0.5*nuVtSqSum[3]*temp_diff[18]+0.5*nuVtSqSum[2]*temp_diff[11]; 
  diff_incr[20] = 0.5*nuVtSqSum[0]*temp_diff[20]+0.5*nuVtSqSum[1]*temp_diff[13]+0.5*nuVtSqSum[2]*temp_diff[12]+0.5*nuVtSqSum[3]*temp_diff[5]; 
  diff_incr[21] = 0.5*nuVtSqSum[2]*temp_diff[27]+0.5*nuVtSqSum[3]*temp_diff[22]+0.5*nuVtSqSum[0]*temp_diff[21]+0.5*nuVtSqSum[1]*temp_diff[14]; 
  diff_incr[22] = 0.5*nuVtSqSum[1]*temp_diff[27]+0.5*nuVtSqSum[0]*temp_diff[22]+0.5*nuVtSqSum[3]*temp_diff[21]+0.5*nuVtSqSum[2]*temp_diff[14]; 
  diff_incr[23] = 0.5*nuVtSqSum[2]*temp_diff[28]+0.5*nuVtSqSum[3]*temp_diff[24]+0.5*nuVtSqSum[0]*temp_diff[23]+0.5*nuVtSqSum[1]*temp_diff[15]; 
  diff_incr[24] = 0.5*nuVtSqSum[1]*temp_diff[28]+0.5*nuVtSqSum[0]*temp_diff[24]+0.5*nuVtSqSum[3]*temp_diff[23]+0.5*nuVtSqSum[2]*temp_diff[15]; 
  diff_incr[25] = 0.5*nuVtSqSum[3]*temp_diff[31]+0.5*nuVtSqSum[2]*temp_diff[30]+0.5*nuVtSqSum[1]*temp_diff[29]+0.5*nuVtSqSum[0]*temp_diff[25]; 
  diff_incr[26] = 0.5*nuVtSqSum[0]*temp_diff[26]+0.5*nuVtSqSum[1]*temp_diff[19]+0.5*nuVtSqSum[2]*temp_diff[18]+0.5*nuVtSqSum[3]*temp_diff[11]; 
  diff_incr[27] = 0.5*nuVtSqSum[0]*temp_diff[27]+0.5*nuVtSqSum[1]*temp_diff[22]+0.5*nuVtSqSum[2]*temp_diff[21]+0.5*nuVtSqSum[3]*temp_diff[14]; 
  diff_incr[28] = 0.5*nuVtSqSum[0]*temp_diff[28]+0.5*nuVtSqSum[1]*temp_diff[24]+0.5*nuVtSqSum[2]*temp_diff[23]+0.5*nuVtSqSum[3]*temp_diff[15]; 
  diff_incr[29] = 0.5*nuVtSqSum[2]*temp_diff[31]+0.5*nuVtSqSum[3]*temp_diff[30]+0.5*nuVtSqSum[0]*temp_diff[29]+0.5*nuVtSqSum[1]*temp_diff[25]; 
  diff_incr[30] = 0.5*nuVtSqSum[1]*temp_diff[31]+0.5*nuVtSqSum[0]*temp_diff[30]+0.5*nuVtSqSum[3]*temp_diff[29]+0.5*nuVtSqSum[2]*temp_diff[25]; 
  diff_incr[31] = 0.5*nuVtSqSum[0]*temp_diff[31]+0.5*nuVtSqSum[1]*temp_diff[30]+0.5*nuVtSqSum[2]*temp_diff[29]+0.5*nuVtSqSum[3]*temp_diff[25]; 

  edge_incr[0] = 0.5*nuVtSqSum[3]*temp_edge[6]+0.5*nuVtSqSum[2]*temp_edge[2]+0.5*nuVtSqSum[1]*temp_edge[1]+0.5*nuVtSqSum[0]*temp_edge[0]; 
  edge_incr[1] = 0.5*nuVtSqSum[2]*temp_edge[6]+0.5*temp_edge[2]*nuVtSqSum[3]+0.5*nuVtSqSum[0]*temp_edge[1]+0.5*temp_edge[0]*nuVtSqSum[1]; 
  edge_incr[2] = 0.5*nuVtSqSum[1]*temp_edge[6]+0.5*temp_edge[1]*nuVtSqSum[3]+0.5*nuVtSqSum[0]*temp_edge[2]+0.5*temp_edge[0]*nuVtSqSum[2]; 
  edge_incr[3] = 0.5*nuVtSqSum[3]*temp_edge[16]+0.5*nuVtSqSum[2]*temp_edge[8]+0.5*nuVtSqSum[1]*temp_edge[7]+0.5*nuVtSqSum[0]*temp_edge[3]; 
  edge_incr[4] = 0.5*nuVtSqSum[3]*temp_edge[17]+0.5*nuVtSqSum[2]*temp_edge[10]+0.5*nuVtSqSum[1]*temp_edge[9]+0.5*nuVtSqSum[0]*temp_edge[4]; 
  edge_incr[5] = 0.5*nuVtSqSum[3]*temp_edge[20]+0.5*nuVtSqSum[2]*temp_edge[13]+0.5*nuVtSqSum[1]*temp_edge[12]+0.5*nuVtSqSum[0]*temp_edge[5]; 
  edge_incr[6] = 0.5*nuVtSqSum[0]*temp_edge[6]+0.5*temp_edge[0]*nuVtSqSum[3]+0.5*nuVtSqSum[1]*temp_edge[2]+0.5*temp_edge[1]*nuVtSqSum[2]; 
  edge_incr[7] = 0.5*nuVtSqSum[2]*temp_edge[16]+0.5*nuVtSqSum[3]*temp_edge[8]+0.5*nuVtSqSum[0]*temp_edge[7]+0.5*nuVtSqSum[1]*temp_edge[3]; 
  edge_incr[8] = 0.5*nuVtSqSum[1]*temp_edge[16]+0.5*nuVtSqSum[0]*temp_edge[8]+0.5*nuVtSqSum[3]*temp_edge[7]+0.5*nuVtSqSum[2]*temp_edge[3]; 
  edge_incr[9] = 0.5*nuVtSqSum[2]*temp_edge[17]+0.5*nuVtSqSum[3]*temp_edge[10]+0.5*nuVtSqSum[0]*temp_edge[9]+0.5*nuVtSqSum[1]*temp_edge[4]; 
  edge_incr[10] = 0.5*nuVtSqSum[1]*temp_edge[17]+0.5*nuVtSqSum[0]*temp_edge[10]+0.5*nuVtSqSum[3]*temp_edge[9]+0.5*nuVtSqSum[2]*temp_edge[4]; 
  edge_incr[11] = 0.5*nuVtSqSum[3]*temp_edge[26]+0.5*nuVtSqSum[2]*temp_edge[19]+0.5*nuVtSqSum[1]*temp_edge[18]+0.5*nuVtSqSum[0]*temp_edge[11]; 
  edge_incr[12] = 0.5*nuVtSqSum[2]*temp_edge[20]+0.5*nuVtSqSum[3]*temp_edge[13]+0.5*nuVtSqSum[0]*temp_edge[12]+0.5*nuVtSqSum[1]*temp_edge[5]; 
  edge_incr[13] = 0.5*nuVtSqSum[1]*temp_edge[20]+0.5*nuVtSqSum[0]*temp_edge[13]+0.5*nuVtSqSum[3]*temp_edge[12]+0.5*nuVtSqSum[2]*temp_edge[5]; 
  edge_incr[14] = 0.5*nuVtSqSum[3]*temp_edge[27]+0.5*nuVtSqSum[2]*temp_edge[22]+0.5*nuVtSqSum[1]*temp_edge[21]+0.5*nuVtSqSum[0]*temp_edge[14]; 
  edge_incr[15] = 0.5*nuVtSqSum[3]*temp_edge[28]+0.5*nuVtSqSum[2]*temp_edge[24]+0.5*nuVtSqSum[1]*temp_edge[23]+0.5*nuVtSqSum[0]*temp_edge[15]; 
  edge_incr[16] = 0.5*nuVtSqSum[0]*temp_edge[16]+0.5*nuVtSqSum[1]*temp_edge[8]+0.5*nuVtSqSum[2]*temp_edge[7]+0.5*nuVtSqSum[3]*temp_edge[3]; 
  edge_incr[17] = 0.5*nuVtSqSum[0]*temp_edge[17]+0.5*nuVtSqSum[1]*temp_edge[10]+0.5*nuVtSqSum[2]*temp_edge[9]+0.5*nuVtSqSum[3]*temp_edge[4]; 
  edge_incr[18] = 0.5*nuVtSqSum[2]*temp_edge[26]+0.5*nuVtSqSum[3]*temp_edge[19]+0.5*nuVtSqSum[0]*temp_edge[18]+0.5*nuVtSqSum[1]*temp_edge[11]; 
  edge_incr[19] = 0.5*nuVtSqSum[1]*temp_edge[26]+0.5*nuVtSqSum[0]*temp_edge[19]+0.5*nuVtSqSum[3]*temp_edge[18]+0.5*nuVtSqSum[2]*temp_edge[11]; 
  edge_incr[20] = 0.5*nuVtSqSum[0]*temp_edge[20]+0.5*nuVtSqSum[1]*temp_edge[13]+0.5*nuVtSqSum[2]*temp_edge[12]+0.5*nuVtSqSum[3]*temp_edge[5]; 
  edge_incr[21] = 0.5*nuVtSqSum[2]*temp_edge[27]+0.5*nuVtSqSum[3]*temp_edge[22]+0.5*nuVtSqSum[0]*temp_edge[21]+0.5*nuVtSqSum[1]*temp_edge[14]; 
  edge_incr[22] = 0.5*nuVtSqSum[1]*temp_edge[27]+0.5*nuVtSqSum[0]*temp_edge[22]+0.5*nuVtSqSum[3]*temp_edge[21]+0.5*nuVtSqSum[2]*temp_edge[14]; 
  edge_incr[23] = 0.5*nuVtSqSum[2]*temp_edge[28]+0.5*nuVtSqSum[3]*temp_edge[24]+0.5*nuVtSqSum[0]*temp_edge[23]+0.5*nuVtSqSum[1]*temp_edge[15]; 
  edge_incr[24] = 0.5*nuVtSqSum[1]*temp_edge[28]+0.5*nuVtSqSum[0]*temp_edge[24]+0.5*nuVtSqSum[3]*temp_edge[23]+0.5*nuVtSqSum[2]*temp_edge[15]; 
  edge_incr[25] = 0.5*nuVtSqSum[3]*temp_edge[31]+0.5*nuVtSqSum[2]*temp_edge[30]+0.5*nuVtSqSum[1]*temp_edge[29]+0.5*nuVtSqSum[0]*temp_edge[25]; 
  edge_incr[26] = 0.5*nuVtSqSum[0]*temp_edge[26]+0.5*nuVtSqSum[1]*temp_edge[19]+0.5*nuVtSqSum[2]*temp_edge[18]+0.5*nuVtSqSum[3]*temp_edge[11]; 
  edge_incr[27] = 0.5*nuVtSqSum[0]*temp_edge[27]+0.5*nuVtSqSum[1]*temp_edge[22]+0.5*nuVtSqSum[2]*temp_edge[21]+0.5*nuVtSqSum[3]*temp_edge[14]; 
  edge_incr[28] = 0.5*nuVtSqSum[0]*temp_edge[28]+0.5*nuVtSqSum[1]*temp_edge[24]+0.5*nuVtSqSum[2]*temp_edge[23]+0.5*nuVtSqSum[3]*temp_edge[15]; 
  edge_incr[29] = 0.5*nuVtSqSum[2]*temp_edge[31]+0.5*nuVtSqSum[3]*temp_edge[30]+0.5*nuVtSqSum[0]*temp_edge[29]+0.5*nuVtSqSum[1]*temp_edge[25]; 
  edge_incr[30] = 0.5*nuVtSqSum[1]*temp_edge[31]+0.5*nuVtSqSum[0]*temp_edge[30]+0.5*nuVtSqSum[3]*temp_edge[29]+0.5*nuVtSqSum[2]*temp_edge[25]; 
  edge_incr[31] = 0.5*nuVtSqSum[0]*temp_edge[31]+0.5*nuVtSqSum[1]*temp_edge[30]+0.5*nuVtSqSum[2]*temp_edge[29]+0.5*nuVtSqSum[3]*temp_edge[25]; 

  } 

  out[0] += edge_incr[0]*rdvSq4+diff_incr[0]*rdvSq4+vol_incr[0]; 
  out[1] += edge_incr[1]*rdvSq4+diff_incr[1]*rdvSq4+vol_incr[1]; 
  out[2] += edge_incr[2]*rdvSq4+diff_incr[2]*rdvSq4+vol_incr[2]; 
  out[3] += edge_incr[3]*rdvSq4+diff_incr[3]*rdvSq4+vol_incr[3]; 
  out[4] += edge_incr[4]*rdvSq4+diff_incr[4]*rdvSq4+vol_incr[4]; 
  out[5] += edge_incr[5]*rdvSq4+diff_incr[5]*rdvSq4+vol_incr[5]; 
  out[6] += edge_incr[6]*rdvSq4+diff_incr[6]*rdvSq4+vol_incr[6]; 
  out[7] += edge_incr[7]*rdvSq4+diff_incr[7]*rdvSq4+vol_incr[7]; 
  out[8] += edge_incr[8]*rdvSq4+diff_incr[8]*rdvSq4+vol_incr[8]; 
  out[9] += edge_incr[9]*rdvSq4+diff_incr[9]*rdvSq4+vol_incr[9]; 
  out[10] += edge_incr[10]*rdvSq4+diff_incr[10]*rdvSq4+vol_incr[10]; 
  out[11] += edge_incr[11]*rdvSq4+diff_incr[11]*rdvSq4+vol_incr[11]; 
  out[12] += edge_incr[12]*rdvSq4+diff_incr[12]*rdvSq4+vol_incr[12]; 
  out[13] += edge_incr[13]*rdvSq4+diff_incr[13]*rdvSq4+vol_incr[13]; 
  out[14] += edge_incr[14]*rdvSq4+diff_incr[14]*rdvSq4+vol_incr[14]; 
  out[15] += edge_incr[15]*rdvSq4+diff_incr[15]*rdvSq4+vol_incr[15]; 
  out[16] += edge_incr[16]*rdvSq4+diff_incr[16]*rdvSq4+vol_incr[16]; 
  out[17] += edge_incr[17]*rdvSq4+diff_incr[17]*rdvSq4+vol_incr[17]; 
  out[18] += edge_incr[18]*rdvSq4+diff_incr[18]*rdvSq4+vol_incr[18]; 
  out[19] += edge_incr[19]*rdvSq4+diff_incr[19]*rdvSq4+vol_incr[19]; 
  out[20] += edge_incr[20]*rdvSq4+diff_incr[20]*rdvSq4+vol_incr[20]; 
  out[21] += edge_incr[21]*rdvSq4+diff_incr[21]*rdvSq4+vol_incr[21]; 
  out[22] += edge_incr[22]*rdvSq4+diff_incr[22]*rdvSq4+vol_incr[22]; 
  out[23] += edge_incr[23]*rdvSq4+diff_incr[23]*rdvSq4+vol_incr[23]; 
  out[24] += edge_incr[24]*rdvSq4+diff_incr[24]*rdvSq4+vol_incr[24]; 
  out[25] += edge_incr[25]*rdvSq4+diff_incr[25]*rdvSq4+vol_incr[25]; 
  out[26] += edge_incr[26]*rdvSq4+diff_incr[26]*rdvSq4+vol_incr[26]; 
  out[27] += edge_incr[27]*rdvSq4+diff_incr[27]*rdvSq4+vol_incr[27]; 
  out[28] += edge_incr[28]*rdvSq4+diff_incr[28]*rdvSq4+vol_incr[28]; 
  out[29] += edge_incr[29]*rdvSq4+diff_incr[29]*rdvSq4+vol_incr[29]; 
  out[30] += edge_incr[30]*rdvSq4+diff_incr[30]*rdvSq4+vol_incr[30]; 
  out[31] += edge_incr[31]*rdvSq4+diff_incr[31]*rdvSq4+vol_incr[31]; 
} 
