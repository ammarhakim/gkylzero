#include <gkyl_vlasov_lbo_kernels.h> 
GKYL_CU_DH void vlasov_lbo_diff_boundary_surfvy_1x2v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[3]:         Cell-center coordinates. 
  // dxv[3]:       Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[4]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[2]: sum of thermal speeds squared time their respective collisionalities. 
  // fSkin/Edge:    Distribution function in cells 
  // out:           Incremented distribution function in cell 
  double rdvSq4 = 4.0/(dxv[2]*dxv[2]); 

  double temp_diff[8] = {0.0}; 
  double temp_edge[8] = {0.0}; 
  double diff_incr[8] = {0.0}; 
  double edge_incr[8] = {0.0}; 
  double vol_incr[8] = {0.0}; 


  if (edge == -1) { 

  temp_diff[0] = (-0.5412658773652741*fSkin[3])-0.5412658773652741*fEdge[3]-0.5625*fSkin[0]+0.5625*fEdge[0]; 
  temp_diff[1] = (-0.5412658773652741*fSkin[5])-0.5412658773652741*fEdge[5]-0.5625*fSkin[1]+0.5625*fEdge[1]; 
  temp_diff[2] = (-0.5412658773652741*fSkin[6])-0.5412658773652741*fEdge[6]-0.5625*fSkin[2]+0.5625*fEdge[2]; 
  temp_diff[3] = (-1.4375*fSkin[3])-0.4375*fEdge[3]-1.407291281149713*fSkin[0]+0.5412658773652741*fEdge[0]; 
  temp_diff[4] = (-0.5412658773652741*fSkin[7])-0.5412658773652741*fEdge[7]-0.5625*fSkin[4]+0.5625*fEdge[4]; 
  temp_diff[5] = (-1.4375*fSkin[5])-0.4375*fEdge[5]-1.407291281149713*fSkin[1]+0.5412658773652741*fEdge[1]; 
  temp_diff[6] = (-1.4375*fSkin[6])-0.4375*fEdge[6]-1.407291281149713*fSkin[2]+0.5412658773652741*fEdge[2]; 
  temp_diff[7] = (-1.4375*fSkin[7])-0.4375*fEdge[7]-1.407291281149713*fSkin[4]+0.5412658773652741*fEdge[4]; 

  temp_edge[3] = 0.8660254037844386*fSkin[0]-1.5*fSkin[3]; 
  temp_edge[5] = 0.8660254037844386*fSkin[1]-1.5*fSkin[5]; 
  temp_edge[6] = 0.8660254037844386*fSkin[2]-1.5*fSkin[6]; 
  temp_edge[7] = 0.8660254037844386*fSkin[4]-1.5*fSkin[7]; 

  diff_incr[0] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[1]+0.7071067811865475*nuVtSqSum[0]*temp_diff[0]; 
  diff_incr[1] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[1]+0.7071067811865475*temp_diff[0]*nuVtSqSum[1]; 
  diff_incr[2] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[4]+0.7071067811865475*nuVtSqSum[0]*temp_diff[2]; 
  diff_incr[3] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[5]+0.7071067811865475*nuVtSqSum[0]*temp_diff[3]; 
  diff_incr[4] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[4]+0.7071067811865475*nuVtSqSum[1]*temp_diff[2]; 
  diff_incr[5] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[5]+0.7071067811865475*nuVtSqSum[1]*temp_diff[3]; 
  diff_incr[6] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[7]+0.7071067811865475*nuVtSqSum[0]*temp_diff[6]; 
  diff_incr[7] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[7]+0.7071067811865475*nuVtSqSum[1]*temp_diff[6]; 

  edge_incr[0] = 0.7071067811865475*nuVtSqSum[1]*temp_edge[1]+0.7071067811865475*nuVtSqSum[0]*temp_edge[0]; 
  edge_incr[1] = 0.7071067811865475*nuVtSqSum[0]*temp_edge[1]+0.7071067811865475*temp_edge[0]*nuVtSqSum[1]; 
  edge_incr[2] = 0.7071067811865475*nuVtSqSum[1]*temp_edge[4]+0.7071067811865475*nuVtSqSum[0]*temp_edge[2]; 
  edge_incr[3] = 0.7071067811865475*nuVtSqSum[1]*temp_edge[5]+0.7071067811865475*nuVtSqSum[0]*temp_edge[3]; 
  edge_incr[4] = 0.7071067811865475*nuVtSqSum[0]*temp_edge[4]+0.7071067811865475*nuVtSqSum[1]*temp_edge[2]; 
  edge_incr[5] = 0.7071067811865475*nuVtSqSum[0]*temp_edge[5]+0.7071067811865475*nuVtSqSum[1]*temp_edge[3]; 
  edge_incr[6] = 0.7071067811865475*nuVtSqSum[1]*temp_edge[7]+0.7071067811865475*nuVtSqSum[0]*temp_edge[6]; 
  edge_incr[7] = 0.7071067811865475*nuVtSqSum[0]*temp_edge[7]+0.7071067811865475*nuVtSqSum[1]*temp_edge[6]; 


  } else { 

  temp_diff[0] = 0.5412658773652741*fSkin[3]+0.5412658773652741*fEdge[3]-0.5625*fSkin[0]+0.5625*fEdge[0]; 
  temp_diff[1] = 0.5412658773652741*fSkin[5]+0.5412658773652741*fEdge[5]-0.5625*fSkin[1]+0.5625*fEdge[1]; 
  temp_diff[2] = 0.5412658773652741*fSkin[6]+0.5412658773652741*fEdge[6]-0.5625*fSkin[2]+0.5625*fEdge[2]; 
  temp_diff[3] = (-1.4375*fSkin[3])-0.4375*fEdge[3]+1.407291281149713*fSkin[0]-0.5412658773652741*fEdge[0]; 
  temp_diff[4] = 0.5412658773652741*fSkin[7]+0.5412658773652741*fEdge[7]-0.5625*fSkin[4]+0.5625*fEdge[4]; 
  temp_diff[5] = (-1.4375*fSkin[5])-0.4375*fEdge[5]+1.407291281149713*fSkin[1]-0.5412658773652741*fEdge[1]; 
  temp_diff[6] = (-1.4375*fSkin[6])-0.4375*fEdge[6]+1.407291281149713*fSkin[2]-0.5412658773652741*fEdge[2]; 
  temp_diff[7] = (-1.4375*fSkin[7])-0.4375*fEdge[7]+1.407291281149713*fSkin[4]-0.5412658773652741*fEdge[4]; 

  temp_edge[3] = (-1.5*fSkin[3])-0.8660254037844386*fSkin[0]; 
  temp_edge[5] = (-1.5*fSkin[5])-0.8660254037844386*fSkin[1]; 
  temp_edge[6] = (-1.5*fSkin[6])-0.8660254037844386*fSkin[2]; 
  temp_edge[7] = (-1.5*fSkin[7])-0.8660254037844386*fSkin[4]; 

  diff_incr[0] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[1]+0.7071067811865475*nuVtSqSum[0]*temp_diff[0]; 
  diff_incr[1] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[1]+0.7071067811865475*temp_diff[0]*nuVtSqSum[1]; 
  diff_incr[2] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[4]+0.7071067811865475*nuVtSqSum[0]*temp_diff[2]; 
  diff_incr[3] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[5]+0.7071067811865475*nuVtSqSum[0]*temp_diff[3]; 
  diff_incr[4] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[4]+0.7071067811865475*nuVtSqSum[1]*temp_diff[2]; 
  diff_incr[5] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[5]+0.7071067811865475*nuVtSqSum[1]*temp_diff[3]; 
  diff_incr[6] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[7]+0.7071067811865475*nuVtSqSum[0]*temp_diff[6]; 
  diff_incr[7] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[7]+0.7071067811865475*nuVtSqSum[1]*temp_diff[6]; 

  edge_incr[0] = 0.7071067811865475*nuVtSqSum[1]*temp_edge[1]+0.7071067811865475*nuVtSqSum[0]*temp_edge[0]; 
  edge_incr[1] = 0.7071067811865475*nuVtSqSum[0]*temp_edge[1]+0.7071067811865475*temp_edge[0]*nuVtSqSum[1]; 
  edge_incr[2] = 0.7071067811865475*nuVtSqSum[1]*temp_edge[4]+0.7071067811865475*nuVtSqSum[0]*temp_edge[2]; 
  edge_incr[3] = 0.7071067811865475*nuVtSqSum[1]*temp_edge[5]+0.7071067811865475*nuVtSqSum[0]*temp_edge[3]; 
  edge_incr[4] = 0.7071067811865475*nuVtSqSum[0]*temp_edge[4]+0.7071067811865475*nuVtSqSum[1]*temp_edge[2]; 
  edge_incr[5] = 0.7071067811865475*nuVtSqSum[0]*temp_edge[5]+0.7071067811865475*nuVtSqSum[1]*temp_edge[3]; 
  edge_incr[6] = 0.7071067811865475*nuVtSqSum[1]*temp_edge[7]+0.7071067811865475*nuVtSqSum[0]*temp_edge[6]; 
  edge_incr[7] = 0.7071067811865475*nuVtSqSum[0]*temp_edge[7]+0.7071067811865475*nuVtSqSum[1]*temp_edge[6]; 

  } 

  out[0] += edge_incr[0]*rdvSq4+diff_incr[0]*rdvSq4+vol_incr[0]; 
  out[1] += edge_incr[1]*rdvSq4+diff_incr[1]*rdvSq4+vol_incr[1]; 
  out[2] += edge_incr[2]*rdvSq4+diff_incr[2]*rdvSq4+vol_incr[2]; 
  out[3] += edge_incr[3]*rdvSq4+diff_incr[3]*rdvSq4+vol_incr[3]; 
  out[4] += edge_incr[4]*rdvSq4+diff_incr[4]*rdvSq4+vol_incr[4]; 
  out[5] += edge_incr[5]*rdvSq4+diff_incr[5]*rdvSq4+vol_incr[5]; 
  out[6] += edge_incr[6]*rdvSq4+diff_incr[6]*rdvSq4+vol_incr[6]; 
  out[7] += edge_incr[7]*rdvSq4+diff_incr[7]*rdvSq4+vol_incr[7]; 
} 
