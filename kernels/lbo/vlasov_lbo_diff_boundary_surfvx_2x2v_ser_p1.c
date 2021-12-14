#include <gkyl_vlasov_lbo_kernels.h> 
GKYL_CU_DH void vlasov_lbo_diff_boundary_surfvx_2x2v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[4]:         Cell-center coordinates. 
  // dxv[4]:       Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[8]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[4]: sum of thermal speeds squared time their respective collisionalities. 
  // fSkin/Edge:    Distribution function in cells 
  // out:           Incremented distribution function in cell 
  double rdvSq4 = 4.0/(dxv[2]*dxv[2]); 

  double temp_diff[16] = {0.0}; 
  double temp_edge[16] = {0.0}; 
  double diff_incr[16] = {0.0}; 
  double edge_incr[16] = {0.0}; 
  double vol_incr[16] = {0.0}; 


  if (edge == -1) { 

  temp_diff[0] = (-0.5412658773652741*fSkin[3])-0.5412658773652741*fEdge[3]-0.5625*fSkin[0]+0.5625*fEdge[0]; 
  temp_diff[1] = (-0.5412658773652741*fSkin[6])-0.5412658773652741*fEdge[6]-0.5625*fSkin[1]+0.5625*fEdge[1]; 
  temp_diff[2] = (-0.5412658773652741*fSkin[7])-0.5412658773652741*fEdge[7]-0.5625*fSkin[2]+0.5625*fEdge[2]; 
  temp_diff[3] = (-1.4375*fSkin[3])-0.4375*fEdge[3]-1.407291281149713*fSkin[0]+0.5412658773652741*fEdge[0]; 
  temp_diff[4] = (-0.5412658773652741*fSkin[10])-0.5412658773652741*fEdge[10]-0.5625*fSkin[4]+0.5625*fEdge[4]; 
  temp_diff[5] = (-0.5412658773652741*fSkin[11])-0.5412658773652741*fEdge[11]-0.5625*fSkin[5]+0.5625*fEdge[5]; 
  temp_diff[6] = (-1.4375*fSkin[6])-0.4375*fEdge[6]-1.407291281149713*fSkin[1]+0.5412658773652741*fEdge[1]; 
  temp_diff[7] = (-1.4375*fSkin[7])-0.4375*fEdge[7]-1.407291281149713*fSkin[2]+0.5412658773652741*fEdge[2]; 
  temp_diff[8] = (-0.5412658773652741*fSkin[13])-0.5412658773652741*fEdge[13]-0.5625*fSkin[8]+0.5625*fEdge[8]; 
  temp_diff[9] = (-0.5412658773652741*fSkin[14])-0.5412658773652741*fEdge[14]-0.5625*fSkin[9]+0.5625*fEdge[9]; 
  temp_diff[10] = (-1.4375*fSkin[10])-0.4375*fEdge[10]-1.407291281149713*fSkin[4]+0.5412658773652741*fEdge[4]; 
  temp_diff[11] = (-1.4375*fSkin[11])-0.4375*fEdge[11]-1.407291281149713*fSkin[5]+0.5412658773652741*fEdge[5]; 
  temp_diff[12] = (-0.5412658773652741*fSkin[15])-0.5412658773652741*fEdge[15]-0.5625*fSkin[12]+0.5625*fEdge[12]; 
  temp_diff[13] = (-1.4375*fSkin[13])-0.4375*fEdge[13]-1.407291281149713*fSkin[8]+0.5412658773652741*fEdge[8]; 
  temp_diff[14] = (-1.4375*fSkin[14])-0.4375*fEdge[14]-1.407291281149713*fSkin[9]+0.5412658773652741*fEdge[9]; 
  temp_diff[15] = (-1.4375*fSkin[15])-0.4375*fEdge[15]-1.407291281149713*fSkin[12]+0.5412658773652741*fEdge[12]; 

  temp_edge[3] = 0.8660254037844386*fSkin[0]-1.5*fSkin[3]; 
  temp_edge[6] = 0.8660254037844386*fSkin[1]-1.5*fSkin[6]; 
  temp_edge[7] = 0.8660254037844386*fSkin[2]-1.5*fSkin[7]; 
  temp_edge[10] = 0.8660254037844386*fSkin[4]-1.5*fSkin[10]; 
  temp_edge[11] = 0.8660254037844386*fSkin[5]-1.5*fSkin[11]; 
  temp_edge[13] = 0.8660254037844386*fSkin[8]-1.5*fSkin[13]; 
  temp_edge[14] = 0.8660254037844386*fSkin[9]-1.5*fSkin[14]; 
  temp_edge[15] = 0.8660254037844386*fSkin[12]-1.5*fSkin[15]; 

  diff_incr[0] = 0.5*nuVtSqSum[3]*temp_diff[5]+0.5*nuVtSqSum[2]*temp_diff[2]+0.5*nuVtSqSum[1]*temp_diff[1]+0.5*nuVtSqSum[0]*temp_diff[0]; 
  diff_incr[1] = 0.5*nuVtSqSum[2]*temp_diff[5]+0.5*temp_diff[2]*nuVtSqSum[3]+0.5*nuVtSqSum[0]*temp_diff[1]+0.5*temp_diff[0]*nuVtSqSum[1]; 
  diff_incr[2] = 0.5*nuVtSqSum[1]*temp_diff[5]+0.5*temp_diff[1]*nuVtSqSum[3]+0.5*nuVtSqSum[0]*temp_diff[2]+0.5*temp_diff[0]*nuVtSqSum[2]; 
  diff_incr[3] = 0.5*nuVtSqSum[3]*temp_diff[11]+0.5*nuVtSqSum[2]*temp_diff[7]+0.5*nuVtSqSum[1]*temp_diff[6]+0.5*nuVtSqSum[0]*temp_diff[3]; 
  diff_incr[4] = 0.5*nuVtSqSum[3]*temp_diff[12]+0.5*nuVtSqSum[2]*temp_diff[9]+0.5*nuVtSqSum[1]*temp_diff[8]+0.5*nuVtSqSum[0]*temp_diff[4]; 
  diff_incr[5] = 0.5*nuVtSqSum[0]*temp_diff[5]+0.5*temp_diff[0]*nuVtSqSum[3]+0.5*nuVtSqSum[1]*temp_diff[2]+0.5*temp_diff[1]*nuVtSqSum[2]; 
  diff_incr[6] = 0.5*nuVtSqSum[2]*temp_diff[11]+0.5*nuVtSqSum[3]*temp_diff[7]+0.5*nuVtSqSum[0]*temp_diff[6]+0.5*nuVtSqSum[1]*temp_diff[3]; 
  diff_incr[7] = 0.5*nuVtSqSum[1]*temp_diff[11]+0.5*nuVtSqSum[0]*temp_diff[7]+0.5*nuVtSqSum[3]*temp_diff[6]+0.5*nuVtSqSum[2]*temp_diff[3]; 
  diff_incr[8] = 0.5*nuVtSqSum[2]*temp_diff[12]+0.5*nuVtSqSum[3]*temp_diff[9]+0.5*nuVtSqSum[0]*temp_diff[8]+0.5*nuVtSqSum[1]*temp_diff[4]; 
  diff_incr[9] = 0.5*nuVtSqSum[1]*temp_diff[12]+0.5*nuVtSqSum[0]*temp_diff[9]+0.5*nuVtSqSum[3]*temp_diff[8]+0.5*nuVtSqSum[2]*temp_diff[4]; 
  diff_incr[10] = 0.5*nuVtSqSum[3]*temp_diff[15]+0.5*nuVtSqSum[2]*temp_diff[14]+0.5*nuVtSqSum[1]*temp_diff[13]+0.5*nuVtSqSum[0]*temp_diff[10]; 
  diff_incr[11] = 0.5*nuVtSqSum[0]*temp_diff[11]+0.5*nuVtSqSum[1]*temp_diff[7]+0.5*nuVtSqSum[2]*temp_diff[6]+0.5*nuVtSqSum[3]*temp_diff[3]; 
  diff_incr[12] = 0.5*nuVtSqSum[0]*temp_diff[12]+0.5*nuVtSqSum[1]*temp_diff[9]+0.5*nuVtSqSum[2]*temp_diff[8]+0.5*nuVtSqSum[3]*temp_diff[4]; 
  diff_incr[13] = 0.5*nuVtSqSum[2]*temp_diff[15]+0.5*nuVtSqSum[3]*temp_diff[14]+0.5*nuVtSqSum[0]*temp_diff[13]+0.5*nuVtSqSum[1]*temp_diff[10]; 
  diff_incr[14] = 0.5*nuVtSqSum[1]*temp_diff[15]+0.5*nuVtSqSum[0]*temp_diff[14]+0.5*nuVtSqSum[3]*temp_diff[13]+0.5*nuVtSqSum[2]*temp_diff[10]; 
  diff_incr[15] = 0.5*nuVtSqSum[0]*temp_diff[15]+0.5*nuVtSqSum[1]*temp_diff[14]+0.5*nuVtSqSum[2]*temp_diff[13]+0.5*nuVtSqSum[3]*temp_diff[10]; 

  edge_incr[0] = 0.5*nuVtSqSum[3]*temp_edge[5]+0.5*nuVtSqSum[2]*temp_edge[2]+0.5*nuVtSqSum[1]*temp_edge[1]+0.5*nuVtSqSum[0]*temp_edge[0]; 
  edge_incr[1] = 0.5*nuVtSqSum[2]*temp_edge[5]+0.5*temp_edge[2]*nuVtSqSum[3]+0.5*nuVtSqSum[0]*temp_edge[1]+0.5*temp_edge[0]*nuVtSqSum[1]; 
  edge_incr[2] = 0.5*nuVtSqSum[1]*temp_edge[5]+0.5*temp_edge[1]*nuVtSqSum[3]+0.5*nuVtSqSum[0]*temp_edge[2]+0.5*temp_edge[0]*nuVtSqSum[2]; 
  edge_incr[3] = 0.5*nuVtSqSum[3]*temp_edge[11]+0.5*nuVtSqSum[2]*temp_edge[7]+0.5*nuVtSqSum[1]*temp_edge[6]+0.5*nuVtSqSum[0]*temp_edge[3]; 
  edge_incr[4] = 0.5*nuVtSqSum[3]*temp_edge[12]+0.5*nuVtSqSum[2]*temp_edge[9]+0.5*nuVtSqSum[1]*temp_edge[8]+0.5*nuVtSqSum[0]*temp_edge[4]; 
  edge_incr[5] = 0.5*nuVtSqSum[0]*temp_edge[5]+0.5*temp_edge[0]*nuVtSqSum[3]+0.5*nuVtSqSum[1]*temp_edge[2]+0.5*temp_edge[1]*nuVtSqSum[2]; 
  edge_incr[6] = 0.5*nuVtSqSum[2]*temp_edge[11]+0.5*nuVtSqSum[3]*temp_edge[7]+0.5*nuVtSqSum[0]*temp_edge[6]+0.5*nuVtSqSum[1]*temp_edge[3]; 
  edge_incr[7] = 0.5*nuVtSqSum[1]*temp_edge[11]+0.5*nuVtSqSum[0]*temp_edge[7]+0.5*nuVtSqSum[3]*temp_edge[6]+0.5*nuVtSqSum[2]*temp_edge[3]; 
  edge_incr[8] = 0.5*nuVtSqSum[2]*temp_edge[12]+0.5*nuVtSqSum[3]*temp_edge[9]+0.5*nuVtSqSum[0]*temp_edge[8]+0.5*nuVtSqSum[1]*temp_edge[4]; 
  edge_incr[9] = 0.5*nuVtSqSum[1]*temp_edge[12]+0.5*nuVtSqSum[0]*temp_edge[9]+0.5*nuVtSqSum[3]*temp_edge[8]+0.5*nuVtSqSum[2]*temp_edge[4]; 
  edge_incr[10] = 0.5*nuVtSqSum[3]*temp_edge[15]+0.5*nuVtSqSum[2]*temp_edge[14]+0.5*nuVtSqSum[1]*temp_edge[13]+0.5*nuVtSqSum[0]*temp_edge[10]; 
  edge_incr[11] = 0.5*nuVtSqSum[0]*temp_edge[11]+0.5*nuVtSqSum[1]*temp_edge[7]+0.5*nuVtSqSum[2]*temp_edge[6]+0.5*nuVtSqSum[3]*temp_edge[3]; 
  edge_incr[12] = 0.5*nuVtSqSum[0]*temp_edge[12]+0.5*nuVtSqSum[1]*temp_edge[9]+0.5*nuVtSqSum[2]*temp_edge[8]+0.5*nuVtSqSum[3]*temp_edge[4]; 
  edge_incr[13] = 0.5*nuVtSqSum[2]*temp_edge[15]+0.5*nuVtSqSum[3]*temp_edge[14]+0.5*nuVtSqSum[0]*temp_edge[13]+0.5*nuVtSqSum[1]*temp_edge[10]; 
  edge_incr[14] = 0.5*nuVtSqSum[1]*temp_edge[15]+0.5*nuVtSqSum[0]*temp_edge[14]+0.5*nuVtSqSum[3]*temp_edge[13]+0.5*nuVtSqSum[2]*temp_edge[10]; 
  edge_incr[15] = 0.5*nuVtSqSum[0]*temp_edge[15]+0.5*nuVtSqSum[1]*temp_edge[14]+0.5*nuVtSqSum[2]*temp_edge[13]+0.5*nuVtSqSum[3]*temp_edge[10]; 


  } else { 

  temp_diff[0] = 0.5412658773652741*fSkin[3]+0.5412658773652741*fEdge[3]-0.5625*fSkin[0]+0.5625*fEdge[0]; 
  temp_diff[1] = 0.5412658773652741*fSkin[6]+0.5412658773652741*fEdge[6]-0.5625*fSkin[1]+0.5625*fEdge[1]; 
  temp_diff[2] = 0.5412658773652741*fSkin[7]+0.5412658773652741*fEdge[7]-0.5625*fSkin[2]+0.5625*fEdge[2]; 
  temp_diff[3] = (-1.4375*fSkin[3])-0.4375*fEdge[3]+1.407291281149713*fSkin[0]-0.5412658773652741*fEdge[0]; 
  temp_diff[4] = 0.5412658773652741*fSkin[10]+0.5412658773652741*fEdge[10]-0.5625*fSkin[4]+0.5625*fEdge[4]; 
  temp_diff[5] = 0.5412658773652741*fSkin[11]+0.5412658773652741*fEdge[11]-0.5625*fSkin[5]+0.5625*fEdge[5]; 
  temp_diff[6] = (-1.4375*fSkin[6])-0.4375*fEdge[6]+1.407291281149713*fSkin[1]-0.5412658773652741*fEdge[1]; 
  temp_diff[7] = (-1.4375*fSkin[7])-0.4375*fEdge[7]+1.407291281149713*fSkin[2]-0.5412658773652741*fEdge[2]; 
  temp_diff[8] = 0.5412658773652741*fSkin[13]+0.5412658773652741*fEdge[13]-0.5625*fSkin[8]+0.5625*fEdge[8]; 
  temp_diff[9] = 0.5412658773652741*fSkin[14]+0.5412658773652741*fEdge[14]-0.5625*fSkin[9]+0.5625*fEdge[9]; 
  temp_diff[10] = (-1.4375*fSkin[10])-0.4375*fEdge[10]+1.407291281149713*fSkin[4]-0.5412658773652741*fEdge[4]; 
  temp_diff[11] = (-1.4375*fSkin[11])-0.4375*fEdge[11]+1.407291281149713*fSkin[5]-0.5412658773652741*fEdge[5]; 
  temp_diff[12] = 0.5412658773652741*fSkin[15]+0.5412658773652741*fEdge[15]-0.5625*fSkin[12]+0.5625*fEdge[12]; 
  temp_diff[13] = (-1.4375*fSkin[13])-0.4375*fEdge[13]+1.407291281149713*fSkin[8]-0.5412658773652741*fEdge[8]; 
  temp_diff[14] = (-1.4375*fSkin[14])-0.4375*fEdge[14]+1.407291281149713*fSkin[9]-0.5412658773652741*fEdge[9]; 
  temp_diff[15] = (-1.4375*fSkin[15])-0.4375*fEdge[15]+1.407291281149713*fSkin[12]-0.5412658773652741*fEdge[12]; 

  temp_edge[3] = (-1.5*fSkin[3])-0.8660254037844386*fSkin[0]; 
  temp_edge[6] = (-1.5*fSkin[6])-0.8660254037844386*fSkin[1]; 
  temp_edge[7] = (-1.5*fSkin[7])-0.8660254037844386*fSkin[2]; 
  temp_edge[10] = (-1.5*fSkin[10])-0.8660254037844386*fSkin[4]; 
  temp_edge[11] = (-1.5*fSkin[11])-0.8660254037844386*fSkin[5]; 
  temp_edge[13] = (-1.5*fSkin[13])-0.8660254037844386*fSkin[8]; 
  temp_edge[14] = (-1.5*fSkin[14])-0.8660254037844386*fSkin[9]; 
  temp_edge[15] = (-1.5*fSkin[15])-0.8660254037844386*fSkin[12]; 

  diff_incr[0] = 0.5*nuVtSqSum[3]*temp_diff[5]+0.5*nuVtSqSum[2]*temp_diff[2]+0.5*nuVtSqSum[1]*temp_diff[1]+0.5*nuVtSqSum[0]*temp_diff[0]; 
  diff_incr[1] = 0.5*nuVtSqSum[2]*temp_diff[5]+0.5*temp_diff[2]*nuVtSqSum[3]+0.5*nuVtSqSum[0]*temp_diff[1]+0.5*temp_diff[0]*nuVtSqSum[1]; 
  diff_incr[2] = 0.5*nuVtSqSum[1]*temp_diff[5]+0.5*temp_diff[1]*nuVtSqSum[3]+0.5*nuVtSqSum[0]*temp_diff[2]+0.5*temp_diff[0]*nuVtSqSum[2]; 
  diff_incr[3] = 0.5*nuVtSqSum[3]*temp_diff[11]+0.5*nuVtSqSum[2]*temp_diff[7]+0.5*nuVtSqSum[1]*temp_diff[6]+0.5*nuVtSqSum[0]*temp_diff[3]; 
  diff_incr[4] = 0.5*nuVtSqSum[3]*temp_diff[12]+0.5*nuVtSqSum[2]*temp_diff[9]+0.5*nuVtSqSum[1]*temp_diff[8]+0.5*nuVtSqSum[0]*temp_diff[4]; 
  diff_incr[5] = 0.5*nuVtSqSum[0]*temp_diff[5]+0.5*temp_diff[0]*nuVtSqSum[3]+0.5*nuVtSqSum[1]*temp_diff[2]+0.5*temp_diff[1]*nuVtSqSum[2]; 
  diff_incr[6] = 0.5*nuVtSqSum[2]*temp_diff[11]+0.5*nuVtSqSum[3]*temp_diff[7]+0.5*nuVtSqSum[0]*temp_diff[6]+0.5*nuVtSqSum[1]*temp_diff[3]; 
  diff_incr[7] = 0.5*nuVtSqSum[1]*temp_diff[11]+0.5*nuVtSqSum[0]*temp_diff[7]+0.5*nuVtSqSum[3]*temp_diff[6]+0.5*nuVtSqSum[2]*temp_diff[3]; 
  diff_incr[8] = 0.5*nuVtSqSum[2]*temp_diff[12]+0.5*nuVtSqSum[3]*temp_diff[9]+0.5*nuVtSqSum[0]*temp_diff[8]+0.5*nuVtSqSum[1]*temp_diff[4]; 
  diff_incr[9] = 0.5*nuVtSqSum[1]*temp_diff[12]+0.5*nuVtSqSum[0]*temp_diff[9]+0.5*nuVtSqSum[3]*temp_diff[8]+0.5*nuVtSqSum[2]*temp_diff[4]; 
  diff_incr[10] = 0.5*nuVtSqSum[3]*temp_diff[15]+0.5*nuVtSqSum[2]*temp_diff[14]+0.5*nuVtSqSum[1]*temp_diff[13]+0.5*nuVtSqSum[0]*temp_diff[10]; 
  diff_incr[11] = 0.5*nuVtSqSum[0]*temp_diff[11]+0.5*nuVtSqSum[1]*temp_diff[7]+0.5*nuVtSqSum[2]*temp_diff[6]+0.5*nuVtSqSum[3]*temp_diff[3]; 
  diff_incr[12] = 0.5*nuVtSqSum[0]*temp_diff[12]+0.5*nuVtSqSum[1]*temp_diff[9]+0.5*nuVtSqSum[2]*temp_diff[8]+0.5*nuVtSqSum[3]*temp_diff[4]; 
  diff_incr[13] = 0.5*nuVtSqSum[2]*temp_diff[15]+0.5*nuVtSqSum[3]*temp_diff[14]+0.5*nuVtSqSum[0]*temp_diff[13]+0.5*nuVtSqSum[1]*temp_diff[10]; 
  diff_incr[14] = 0.5*nuVtSqSum[1]*temp_diff[15]+0.5*nuVtSqSum[0]*temp_diff[14]+0.5*nuVtSqSum[3]*temp_diff[13]+0.5*nuVtSqSum[2]*temp_diff[10]; 
  diff_incr[15] = 0.5*nuVtSqSum[0]*temp_diff[15]+0.5*nuVtSqSum[1]*temp_diff[14]+0.5*nuVtSqSum[2]*temp_diff[13]+0.5*nuVtSqSum[3]*temp_diff[10]; 

  edge_incr[0] = 0.5*nuVtSqSum[3]*temp_edge[5]+0.5*nuVtSqSum[2]*temp_edge[2]+0.5*nuVtSqSum[1]*temp_edge[1]+0.5*nuVtSqSum[0]*temp_edge[0]; 
  edge_incr[1] = 0.5*nuVtSqSum[2]*temp_edge[5]+0.5*temp_edge[2]*nuVtSqSum[3]+0.5*nuVtSqSum[0]*temp_edge[1]+0.5*temp_edge[0]*nuVtSqSum[1]; 
  edge_incr[2] = 0.5*nuVtSqSum[1]*temp_edge[5]+0.5*temp_edge[1]*nuVtSqSum[3]+0.5*nuVtSqSum[0]*temp_edge[2]+0.5*temp_edge[0]*nuVtSqSum[2]; 
  edge_incr[3] = 0.5*nuVtSqSum[3]*temp_edge[11]+0.5*nuVtSqSum[2]*temp_edge[7]+0.5*nuVtSqSum[1]*temp_edge[6]+0.5*nuVtSqSum[0]*temp_edge[3]; 
  edge_incr[4] = 0.5*nuVtSqSum[3]*temp_edge[12]+0.5*nuVtSqSum[2]*temp_edge[9]+0.5*nuVtSqSum[1]*temp_edge[8]+0.5*nuVtSqSum[0]*temp_edge[4]; 
  edge_incr[5] = 0.5*nuVtSqSum[0]*temp_edge[5]+0.5*temp_edge[0]*nuVtSqSum[3]+0.5*nuVtSqSum[1]*temp_edge[2]+0.5*temp_edge[1]*nuVtSqSum[2]; 
  edge_incr[6] = 0.5*nuVtSqSum[2]*temp_edge[11]+0.5*nuVtSqSum[3]*temp_edge[7]+0.5*nuVtSqSum[0]*temp_edge[6]+0.5*nuVtSqSum[1]*temp_edge[3]; 
  edge_incr[7] = 0.5*nuVtSqSum[1]*temp_edge[11]+0.5*nuVtSqSum[0]*temp_edge[7]+0.5*nuVtSqSum[3]*temp_edge[6]+0.5*nuVtSqSum[2]*temp_edge[3]; 
  edge_incr[8] = 0.5*nuVtSqSum[2]*temp_edge[12]+0.5*nuVtSqSum[3]*temp_edge[9]+0.5*nuVtSqSum[0]*temp_edge[8]+0.5*nuVtSqSum[1]*temp_edge[4]; 
  edge_incr[9] = 0.5*nuVtSqSum[1]*temp_edge[12]+0.5*nuVtSqSum[0]*temp_edge[9]+0.5*nuVtSqSum[3]*temp_edge[8]+0.5*nuVtSqSum[2]*temp_edge[4]; 
  edge_incr[10] = 0.5*nuVtSqSum[3]*temp_edge[15]+0.5*nuVtSqSum[2]*temp_edge[14]+0.5*nuVtSqSum[1]*temp_edge[13]+0.5*nuVtSqSum[0]*temp_edge[10]; 
  edge_incr[11] = 0.5*nuVtSqSum[0]*temp_edge[11]+0.5*nuVtSqSum[1]*temp_edge[7]+0.5*nuVtSqSum[2]*temp_edge[6]+0.5*nuVtSqSum[3]*temp_edge[3]; 
  edge_incr[12] = 0.5*nuVtSqSum[0]*temp_edge[12]+0.5*nuVtSqSum[1]*temp_edge[9]+0.5*nuVtSqSum[2]*temp_edge[8]+0.5*nuVtSqSum[3]*temp_edge[4]; 
  edge_incr[13] = 0.5*nuVtSqSum[2]*temp_edge[15]+0.5*nuVtSqSum[3]*temp_edge[14]+0.5*nuVtSqSum[0]*temp_edge[13]+0.5*nuVtSqSum[1]*temp_edge[10]; 
  edge_incr[14] = 0.5*nuVtSqSum[1]*temp_edge[15]+0.5*nuVtSqSum[0]*temp_edge[14]+0.5*nuVtSqSum[3]*temp_edge[13]+0.5*nuVtSqSum[2]*temp_edge[10]; 
  edge_incr[15] = 0.5*nuVtSqSum[0]*temp_edge[15]+0.5*nuVtSqSum[1]*temp_edge[14]+0.5*nuVtSqSum[2]*temp_edge[13]+0.5*nuVtSqSum[3]*temp_edge[10]; 

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
} 
