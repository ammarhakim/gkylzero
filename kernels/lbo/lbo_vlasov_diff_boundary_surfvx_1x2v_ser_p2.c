#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH void lbo_vlasov_diff_boundary_surfvx_1x2v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[3]:         Cell-center coordinates. 
  // dxv[3]:       Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[6]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[3]: sum of thermal speeds squared time their respective collisionalities. 
  // fSkin/Edge:    Distribution function in cells 
  // out:           Incremented distribution function in cell 
  double rdvSq4 = 4.0/(dxv[1]*dxv[1]); 

  double facDiff[3]; 
  // Expand nuVtSqSum in phase basis.
  facDiff[0] = nuVtSqSum[0]; 
  facDiff[1] = nuVtSqSum[1]; 
  facDiff[2] = nuVtSqSum[2]; 

  double temp_diff[20] = {0.0}; 
  double temp_edge[20] = {0.0}; 
  double diff_incr[20] = {0.0}; 
  double edge_incr[20] = {0.0}; 
  double vol_incr[20] = {0.0}; 

  vol_incr[8] = 4.743416490252569*facDiff[2]*fSkin[7]*rdvSq4+4.743416490252569*fSkin[1]*facDiff[1]*rdvSq4+4.743416490252569*fSkin[0]*facDiff[0]*rdvSq4; 
  vol_incr[12] = 4.242640687119286*facDiff[1]*fSkin[7]*rdvSq4+4.242640687119286*fSkin[1]*facDiff[2]*rdvSq4+4.743416490252569*fSkin[0]*facDiff[1]*rdvSq4+4.743416490252569*facDiff[0]*fSkin[1]*rdvSq4; 
  vol_incr[14] = 4.743416490252569*facDiff[2]*fSkin[13]*rdvSq4+4.743416490252569*facDiff[1]*fSkin[5]*rdvSq4+4.743416490252569*facDiff[0]*fSkin[3]*rdvSq4; 
  vol_incr[18] = 4.242640687119286*facDiff[1]*fSkin[13]*rdvSq4+4.242640687119286*facDiff[2]*fSkin[5]*rdvSq4+4.743416490252569*facDiff[0]*fSkin[5]*rdvSq4+4.743416490252569*facDiff[1]*fSkin[3]*rdvSq4; 

  if (edge == -1) { 

  temp_diff[0] = (-0.6708203932499369*fSkin[8])+0.6708203932499369*fEdge[8]-1.190784930203603*fSkin[2]-1.190784930203603*fEdge[2]-0.9375*fSkin[0]+0.9375*fEdge[0]; 
  temp_diff[1] = (-0.6708203932499369*fSkin[12])+0.6708203932499369*fEdge[12]-1.190784930203603*fSkin[4]-1.190784930203603*fEdge[4]-0.9375*fSkin[1]+0.9375*fEdge[1]; 
  temp_diff[2] = (-1.585502557353661*fSkin[8])+0.7382874503707888*fEdge[8]-2.671875*fSkin[2]-1.453125*fEdge[2]-2.056810333988042*fSkin[0]+1.190784930203603*fEdge[0]; 
  temp_diff[3] = (-0.6708203932499369*fSkin[14])+0.6708203932499369*fEdge[14]-1.190784930203603*fSkin[6]-1.190784930203603*fEdge[6]-0.9375*fSkin[3]+0.9375*fEdge[3]; 
  temp_diff[4] = (-1.585502557353661*fSkin[12])+0.7382874503707888*fEdge[12]-2.671875*fSkin[4]-1.453125*fEdge[4]-2.056810333988042*fSkin[1]+1.190784930203603*fEdge[1]; 
  temp_diff[5] = (-0.6708203932499369*fSkin[18])+0.6708203932499369*fEdge[18]-1.190784930203603*fSkin[10]-1.190784930203603*fEdge[10]-0.9375*fSkin[5]+0.9375*fEdge[5]; 
  temp_diff[6] = (-1.585502557353661*fSkin[14])+0.7382874503707888*fEdge[14]-2.671875*fSkin[6]-1.453125*fEdge[6]-2.056810333988042*fSkin[3]+1.190784930203603*fEdge[3]; 
  temp_diff[7] = (-1.190784930203603*fSkin[11])-1.190784930203603*fEdge[11]-0.9375*fSkin[7]+0.9375*fEdge[7]; 
  temp_diff[8] = (-3.140625*fSkin[8])-0.140625*fEdge[8]-5.022775277112744*fSkin[2]-0.3025768239224545*fEdge[2]-3.773364712030896*fSkin[0]+0.4192627457812106*fEdge[0]; 
  temp_diff[9] = (-1.190784930203603*fSkin[16])-1.190784930203603*fEdge[16]-0.9375*fSkin[9]+0.9375*fEdge[9]; 
  temp_diff[10] = (-1.585502557353661*fSkin[18])+0.7382874503707888*fEdge[18]-2.671875*fSkin[10]-1.453125*fEdge[10]-2.056810333988042*fSkin[5]+1.190784930203603*fEdge[5]; 
  temp_diff[11] = (-2.671875*fSkin[11])-1.453125*fEdge[11]-2.056810333988042*fSkin[7]+1.190784930203603*fEdge[7]; 
  temp_diff[12] = (-3.140625*fSkin[12])-0.140625*fEdge[12]-5.022775277112744*fSkin[4]-0.3025768239224544*fEdge[4]-3.773364712030894*fSkin[1]+0.4192627457812105*fEdge[1]; 
  temp_diff[13] = (-1.190784930203603*fSkin[17])-1.190784930203603*fEdge[17]-0.9375*fSkin[13]+0.9375*fEdge[13]; 
  temp_diff[14] = (-3.140625*fSkin[14])-0.140625*fEdge[14]-5.022775277112744*fSkin[6]-0.3025768239224544*fEdge[6]-3.773364712030894*fSkin[3]+0.4192627457812105*fEdge[3]; 
  temp_diff[15] = (-1.190784930203603*fSkin[19])-1.190784930203603*fEdge[19]-0.9375*fSkin[15]+0.9375*fEdge[15]; 
  temp_diff[16] = (-2.671875*fSkin[16])-1.453125*fEdge[16]-2.056810333988042*fSkin[9]+1.190784930203603*fEdge[9]; 
  temp_diff[17] = (-2.671875*fSkin[17])-1.453125*fEdge[17]-2.056810333988042*fSkin[13]+1.190784930203603*fEdge[13]; 
  temp_diff[18] = (-3.140625*fSkin[18])-0.140625*fEdge[18]-5.022775277112744*fSkin[10]-0.3025768239224545*fEdge[10]-3.773364712030896*fSkin[5]+0.4192627457812106*fEdge[5]; 
  temp_diff[19] = (-2.671875*fSkin[19])-1.453125*fEdge[19]-2.056810333988042*fSkin[15]+1.190784930203603*fEdge[15]; 

  temp_edge[2] = 1.936491673103709*fSkin[8]-1.5*fSkin[2]+0.8660254037844386*fSkin[0]; 
  temp_edge[4] = 1.936491673103709*fSkin[12]-1.5*fSkin[4]+0.8660254037844386*fSkin[1]; 
  temp_edge[6] = 1.936491673103709*fSkin[14]-1.5*fSkin[6]+0.8660254037844386*fSkin[3]; 
  temp_edge[8] = (-7.5*fSkin[8])+5.809475019311125*fSkin[2]-3.354101966249685*fSkin[0]; 
  temp_edge[10] = 1.936491673103709*fSkin[18]-1.5*fSkin[10]+0.8660254037844386*fSkin[5]; 
  temp_edge[11] = 0.8660254037844387*fSkin[7]-1.5*fSkin[11]; 
  temp_edge[12] = (-7.5*fSkin[12])+5.809475019311126*fSkin[4]-3.354101966249684*fSkin[1]; 
  temp_edge[14] = (-7.5*fSkin[14])+5.809475019311126*fSkin[6]-3.354101966249684*fSkin[3]; 
  temp_edge[16] = 0.8660254037844387*fSkin[9]-1.5*fSkin[16]; 
  temp_edge[17] = 0.8660254037844387*fSkin[13]-1.5*fSkin[17]; 
  temp_edge[18] = (-7.5*fSkin[18])+5.809475019311125*fSkin[10]-3.354101966249685*fSkin[5]; 
  temp_edge[19] = 0.8660254037844387*fSkin[15]-1.5*fSkin[19]; 

  diff_incr[0] = 0.7071067811865475*nuVtSqSum[2]*temp_diff[7]+0.7071067811865475*nuVtSqSum[1]*temp_diff[1]+0.7071067811865475*nuVtSqSum[0]*temp_diff[0]; 
  diff_incr[1] = 0.6324555320336759*nuVtSqSum[1]*temp_diff[7]+0.6324555320336759*temp_diff[1]*nuVtSqSum[2]+0.7071067811865475*nuVtSqSum[0]*temp_diff[1]+0.7071067811865475*temp_diff[0]*nuVtSqSum[1]; 
  diff_incr[2] = 0.7071067811865475*nuVtSqSum[2]*temp_diff[11]+0.7071067811865475*nuVtSqSum[1]*temp_diff[4]+0.7071067811865475*nuVtSqSum[0]*temp_diff[2]; 
  diff_incr[3] = 0.7071067811865475*nuVtSqSum[2]*temp_diff[13]+0.7071067811865475*nuVtSqSum[1]*temp_diff[5]+0.7071067811865475*nuVtSqSum[0]*temp_diff[3]; 
  diff_incr[4] = 0.632455532033676*nuVtSqSum[1]*temp_diff[11]+0.6324555320336759*nuVtSqSum[2]*temp_diff[4]+0.7071067811865475*nuVtSqSum[0]*temp_diff[4]+0.7071067811865475*nuVtSqSum[1]*temp_diff[2]; 
  diff_incr[5] = 0.632455532033676*nuVtSqSum[1]*temp_diff[13]+0.6324555320336759*nuVtSqSum[2]*temp_diff[5]+0.7071067811865475*nuVtSqSum[0]*temp_diff[5]+0.7071067811865475*nuVtSqSum[1]*temp_diff[3]; 
  diff_incr[6] = 0.7071067811865475*nuVtSqSum[2]*temp_diff[17]+0.7071067811865475*nuVtSqSum[1]*temp_diff[10]+0.7071067811865475*nuVtSqSum[0]*temp_diff[6]; 
  diff_incr[7] = 0.4517539514526256*nuVtSqSum[2]*temp_diff[7]+0.7071067811865475*nuVtSqSum[0]*temp_diff[7]+0.7071067811865475*temp_diff[0]*nuVtSqSum[2]+0.6324555320336759*nuVtSqSum[1]*temp_diff[1]; 
  diff_incr[8] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[12]+0.7071067811865475*nuVtSqSum[0]*temp_diff[8]; 
  diff_incr[9] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[15]+0.7071067811865475*nuVtSqSum[0]*temp_diff[9]; 
  diff_incr[10] = 0.6324555320336759*nuVtSqSum[1]*temp_diff[17]+0.6324555320336759*nuVtSqSum[2]*temp_diff[10]+0.7071067811865475*nuVtSqSum[0]*temp_diff[10]+0.7071067811865475*nuVtSqSum[1]*temp_diff[6]; 
  diff_incr[11] = 0.4517539514526256*nuVtSqSum[2]*temp_diff[11]+0.7071067811865475*nuVtSqSum[0]*temp_diff[11]+0.632455532033676*nuVtSqSum[1]*temp_diff[4]+0.7071067811865475*nuVtSqSum[2]*temp_diff[2]; 
  diff_incr[12] = 0.6324555320336759*nuVtSqSum[2]*temp_diff[12]+0.7071067811865475*nuVtSqSum[0]*temp_diff[12]+0.7071067811865475*nuVtSqSum[1]*temp_diff[8]; 
  diff_incr[13] = 0.4517539514526256*nuVtSqSum[2]*temp_diff[13]+0.7071067811865475*nuVtSqSum[0]*temp_diff[13]+0.632455532033676*nuVtSqSum[1]*temp_diff[5]+0.7071067811865475*nuVtSqSum[2]*temp_diff[3]; 
  diff_incr[14] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[18]+0.7071067811865475*nuVtSqSum[0]*temp_diff[14]; 
  diff_incr[15] = 0.6324555320336759*nuVtSqSum[2]*temp_diff[15]+0.7071067811865475*nuVtSqSum[0]*temp_diff[15]+0.7071067811865475*nuVtSqSum[1]*temp_diff[9]; 
  diff_incr[16] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[19]+0.7071067811865475*nuVtSqSum[0]*temp_diff[16]; 
  diff_incr[17] = 0.4517539514526256*nuVtSqSum[2]*temp_diff[17]+0.7071067811865475*nuVtSqSum[0]*temp_diff[17]+0.6324555320336759*nuVtSqSum[1]*temp_diff[10]+0.7071067811865475*nuVtSqSum[2]*temp_diff[6]; 
  diff_incr[18] = 0.6324555320336759*nuVtSqSum[2]*temp_diff[18]+0.7071067811865475*nuVtSqSum[0]*temp_diff[18]+0.7071067811865475*nuVtSqSum[1]*temp_diff[14]; 
  diff_incr[19] = 0.6324555320336759*nuVtSqSum[2]*temp_diff[19]+0.7071067811865475*nuVtSqSum[0]*temp_diff[19]+0.7071067811865475*nuVtSqSum[1]*temp_diff[16]; 

  edge_incr[0] = 0.7071067811865475*nuVtSqSum[2]*temp_edge[7]+0.7071067811865475*nuVtSqSum[1]*temp_edge[1]+0.7071067811865475*nuVtSqSum[0]*temp_edge[0]; 
  edge_incr[1] = 0.6324555320336759*nuVtSqSum[1]*temp_edge[7]+0.6324555320336759*temp_edge[1]*nuVtSqSum[2]+0.7071067811865475*nuVtSqSum[0]*temp_edge[1]+0.7071067811865475*temp_edge[0]*nuVtSqSum[1]; 
  edge_incr[2] = 0.7071067811865475*nuVtSqSum[2]*temp_edge[11]+0.7071067811865475*nuVtSqSum[1]*temp_edge[4]+0.7071067811865475*nuVtSqSum[0]*temp_edge[2]; 
  edge_incr[3] = 0.7071067811865475*nuVtSqSum[2]*temp_edge[13]+0.7071067811865475*nuVtSqSum[1]*temp_edge[5]+0.7071067811865475*nuVtSqSum[0]*temp_edge[3]; 
  edge_incr[4] = 0.632455532033676*nuVtSqSum[1]*temp_edge[11]+0.6324555320336759*nuVtSqSum[2]*temp_edge[4]+0.7071067811865475*nuVtSqSum[0]*temp_edge[4]+0.7071067811865475*nuVtSqSum[1]*temp_edge[2]; 
  edge_incr[5] = 0.632455532033676*nuVtSqSum[1]*temp_edge[13]+0.6324555320336759*nuVtSqSum[2]*temp_edge[5]+0.7071067811865475*nuVtSqSum[0]*temp_edge[5]+0.7071067811865475*nuVtSqSum[1]*temp_edge[3]; 
  edge_incr[6] = 0.7071067811865475*nuVtSqSum[2]*temp_edge[17]+0.7071067811865475*nuVtSqSum[1]*temp_edge[10]+0.7071067811865475*nuVtSqSum[0]*temp_edge[6]; 
  edge_incr[7] = 0.4517539514526256*nuVtSqSum[2]*temp_edge[7]+0.7071067811865475*nuVtSqSum[0]*temp_edge[7]+0.7071067811865475*temp_edge[0]*nuVtSqSum[2]+0.6324555320336759*nuVtSqSum[1]*temp_edge[1]; 
  edge_incr[8] = 0.7071067811865475*nuVtSqSum[1]*temp_edge[12]+0.7071067811865475*nuVtSqSum[0]*temp_edge[8]; 
  edge_incr[9] = 0.7071067811865475*nuVtSqSum[1]*temp_edge[15]+0.7071067811865475*nuVtSqSum[0]*temp_edge[9]; 
  edge_incr[10] = 0.6324555320336759*nuVtSqSum[1]*temp_edge[17]+0.6324555320336759*nuVtSqSum[2]*temp_edge[10]+0.7071067811865475*nuVtSqSum[0]*temp_edge[10]+0.7071067811865475*nuVtSqSum[1]*temp_edge[6]; 
  edge_incr[11] = 0.4517539514526256*nuVtSqSum[2]*temp_edge[11]+0.7071067811865475*nuVtSqSum[0]*temp_edge[11]+0.632455532033676*nuVtSqSum[1]*temp_edge[4]+0.7071067811865475*nuVtSqSum[2]*temp_edge[2]; 
  edge_incr[12] = 0.6324555320336759*nuVtSqSum[2]*temp_edge[12]+0.7071067811865475*nuVtSqSum[0]*temp_edge[12]+0.7071067811865475*nuVtSqSum[1]*temp_edge[8]; 
  edge_incr[13] = 0.4517539514526256*nuVtSqSum[2]*temp_edge[13]+0.7071067811865475*nuVtSqSum[0]*temp_edge[13]+0.632455532033676*nuVtSqSum[1]*temp_edge[5]+0.7071067811865475*nuVtSqSum[2]*temp_edge[3]; 
  edge_incr[14] = 0.7071067811865475*nuVtSqSum[1]*temp_edge[18]+0.7071067811865475*nuVtSqSum[0]*temp_edge[14]; 
  edge_incr[15] = 0.6324555320336759*nuVtSqSum[2]*temp_edge[15]+0.7071067811865475*nuVtSqSum[0]*temp_edge[15]+0.7071067811865475*nuVtSqSum[1]*temp_edge[9]; 
  edge_incr[16] = 0.7071067811865475*nuVtSqSum[1]*temp_edge[19]+0.7071067811865475*nuVtSqSum[0]*temp_edge[16]; 
  edge_incr[17] = 0.4517539514526256*nuVtSqSum[2]*temp_edge[17]+0.7071067811865475*nuVtSqSum[0]*temp_edge[17]+0.6324555320336759*nuVtSqSum[1]*temp_edge[10]+0.7071067811865475*nuVtSqSum[2]*temp_edge[6]; 
  edge_incr[18] = 0.6324555320336759*nuVtSqSum[2]*temp_edge[18]+0.7071067811865475*nuVtSqSum[0]*temp_edge[18]+0.7071067811865475*nuVtSqSum[1]*temp_edge[14]; 
  edge_incr[19] = 0.6324555320336759*nuVtSqSum[2]*temp_edge[19]+0.7071067811865475*nuVtSqSum[0]*temp_edge[19]+0.7071067811865475*nuVtSqSum[1]*temp_edge[16]; 


  } else { 

  temp_diff[0] = (-0.6708203932499369*fSkin[8])+0.6708203932499369*fEdge[8]+1.190784930203603*fSkin[2]+1.190784930203603*fEdge[2]-0.9375*fSkin[0]+0.9375*fEdge[0]; 
  temp_diff[1] = (-0.6708203932499369*fSkin[12])+0.6708203932499369*fEdge[12]+1.190784930203603*fSkin[4]+1.190784930203603*fEdge[4]-0.9375*fSkin[1]+0.9375*fEdge[1]; 
  temp_diff[2] = 1.585502557353661*fSkin[8]-0.7382874503707888*fEdge[8]-2.671875*fSkin[2]-1.453125*fEdge[2]+2.056810333988042*fSkin[0]-1.190784930203603*fEdge[0]; 
  temp_diff[3] = (-0.6708203932499369*fSkin[14])+0.6708203932499369*fEdge[14]+1.190784930203603*fSkin[6]+1.190784930203603*fEdge[6]-0.9375*fSkin[3]+0.9375*fEdge[3]; 
  temp_diff[4] = 1.585502557353661*fSkin[12]-0.7382874503707888*fEdge[12]-2.671875*fSkin[4]-1.453125*fEdge[4]+2.056810333988042*fSkin[1]-1.190784930203603*fEdge[1]; 
  temp_diff[5] = (-0.6708203932499369*fSkin[18])+0.6708203932499369*fEdge[18]+1.190784930203603*fSkin[10]+1.190784930203603*fEdge[10]-0.9375*fSkin[5]+0.9375*fEdge[5]; 
  temp_diff[6] = 1.585502557353661*fSkin[14]-0.7382874503707888*fEdge[14]-2.671875*fSkin[6]-1.453125*fEdge[6]+2.056810333988042*fSkin[3]-1.190784930203603*fEdge[3]; 
  temp_diff[7] = 1.190784930203603*fSkin[11]+1.190784930203603*fEdge[11]-0.9375*fSkin[7]+0.9375*fEdge[7]; 
  temp_diff[8] = (-3.140625*fSkin[8])-0.140625*fEdge[8]+5.022775277112744*fSkin[2]+0.3025768239224545*fEdge[2]-3.773364712030896*fSkin[0]+0.4192627457812106*fEdge[0]; 
  temp_diff[9] = 1.190784930203603*fSkin[16]+1.190784930203603*fEdge[16]-0.9375*fSkin[9]+0.9375*fEdge[9]; 
  temp_diff[10] = 1.585502557353661*fSkin[18]-0.7382874503707888*fEdge[18]-2.671875*fSkin[10]-1.453125*fEdge[10]+2.056810333988042*fSkin[5]-1.190784930203603*fEdge[5]; 
  temp_diff[11] = (-2.671875*fSkin[11])-1.453125*fEdge[11]+2.056810333988042*fSkin[7]-1.190784930203603*fEdge[7]; 
  temp_diff[12] = (-3.140625*fSkin[12])-0.140625*fEdge[12]+5.022775277112744*fSkin[4]+0.3025768239224544*fEdge[4]-3.773364712030894*fSkin[1]+0.4192627457812105*fEdge[1]; 
  temp_diff[13] = 1.190784930203603*fSkin[17]+1.190784930203603*fEdge[17]-0.9375*fSkin[13]+0.9375*fEdge[13]; 
  temp_diff[14] = (-3.140625*fSkin[14])-0.140625*fEdge[14]+5.022775277112744*fSkin[6]+0.3025768239224544*fEdge[6]-3.773364712030894*fSkin[3]+0.4192627457812105*fEdge[3]; 
  temp_diff[15] = 1.190784930203603*fSkin[19]+1.190784930203603*fEdge[19]-0.9375*fSkin[15]+0.9375*fEdge[15]; 
  temp_diff[16] = (-2.671875*fSkin[16])-1.453125*fEdge[16]+2.056810333988042*fSkin[9]-1.190784930203603*fEdge[9]; 
  temp_diff[17] = (-2.671875*fSkin[17])-1.453125*fEdge[17]+2.056810333988042*fSkin[13]-1.190784930203603*fEdge[13]; 
  temp_diff[18] = (-3.140625*fSkin[18])-0.140625*fEdge[18]+5.022775277112744*fSkin[10]+0.3025768239224545*fEdge[10]-3.773364712030896*fSkin[5]+0.4192627457812106*fEdge[5]; 
  temp_diff[19] = (-2.671875*fSkin[19])-1.453125*fEdge[19]+2.056810333988042*fSkin[15]-1.190784930203603*fEdge[15]; 

  temp_edge[2] = (-1.936491673103709*fSkin[8])-1.5*fSkin[2]-0.8660254037844386*fSkin[0]; 
  temp_edge[4] = (-1.936491673103709*fSkin[12])-1.5*fSkin[4]-0.8660254037844386*fSkin[1]; 
  temp_edge[6] = (-1.936491673103709*fSkin[14])-1.5*fSkin[6]-0.8660254037844386*fSkin[3]; 
  temp_edge[8] = (-7.5*fSkin[8])-5.809475019311125*fSkin[2]-3.354101966249685*fSkin[0]; 
  temp_edge[10] = (-1.936491673103709*fSkin[18])-1.5*fSkin[10]-0.8660254037844386*fSkin[5]; 
  temp_edge[11] = (-1.5*fSkin[11])-0.8660254037844387*fSkin[7]; 
  temp_edge[12] = (-7.5*fSkin[12])-5.809475019311126*fSkin[4]-3.354101966249684*fSkin[1]; 
  temp_edge[14] = (-7.5*fSkin[14])-5.809475019311126*fSkin[6]-3.354101966249684*fSkin[3]; 
  temp_edge[16] = (-1.5*fSkin[16])-0.8660254037844387*fSkin[9]; 
  temp_edge[17] = (-1.5*fSkin[17])-0.8660254037844387*fSkin[13]; 
  temp_edge[18] = (-7.5*fSkin[18])-5.809475019311125*fSkin[10]-3.354101966249685*fSkin[5]; 
  temp_edge[19] = (-1.5*fSkin[19])-0.8660254037844387*fSkin[15]; 

  diff_incr[0] = 0.7071067811865475*nuVtSqSum[2]*temp_diff[7]+0.7071067811865475*nuVtSqSum[1]*temp_diff[1]+0.7071067811865475*nuVtSqSum[0]*temp_diff[0]; 
  diff_incr[1] = 0.6324555320336759*nuVtSqSum[1]*temp_diff[7]+0.6324555320336759*temp_diff[1]*nuVtSqSum[2]+0.7071067811865475*nuVtSqSum[0]*temp_diff[1]+0.7071067811865475*temp_diff[0]*nuVtSqSum[1]; 
  diff_incr[2] = 0.7071067811865475*nuVtSqSum[2]*temp_diff[11]+0.7071067811865475*nuVtSqSum[1]*temp_diff[4]+0.7071067811865475*nuVtSqSum[0]*temp_diff[2]; 
  diff_incr[3] = 0.7071067811865475*nuVtSqSum[2]*temp_diff[13]+0.7071067811865475*nuVtSqSum[1]*temp_diff[5]+0.7071067811865475*nuVtSqSum[0]*temp_diff[3]; 
  diff_incr[4] = 0.632455532033676*nuVtSqSum[1]*temp_diff[11]+0.6324555320336759*nuVtSqSum[2]*temp_diff[4]+0.7071067811865475*nuVtSqSum[0]*temp_diff[4]+0.7071067811865475*nuVtSqSum[1]*temp_diff[2]; 
  diff_incr[5] = 0.632455532033676*nuVtSqSum[1]*temp_diff[13]+0.6324555320336759*nuVtSqSum[2]*temp_diff[5]+0.7071067811865475*nuVtSqSum[0]*temp_diff[5]+0.7071067811865475*nuVtSqSum[1]*temp_diff[3]; 
  diff_incr[6] = 0.7071067811865475*nuVtSqSum[2]*temp_diff[17]+0.7071067811865475*nuVtSqSum[1]*temp_diff[10]+0.7071067811865475*nuVtSqSum[0]*temp_diff[6]; 
  diff_incr[7] = 0.4517539514526256*nuVtSqSum[2]*temp_diff[7]+0.7071067811865475*nuVtSqSum[0]*temp_diff[7]+0.7071067811865475*temp_diff[0]*nuVtSqSum[2]+0.6324555320336759*nuVtSqSum[1]*temp_diff[1]; 
  diff_incr[8] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[12]+0.7071067811865475*nuVtSqSum[0]*temp_diff[8]; 
  diff_incr[9] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[15]+0.7071067811865475*nuVtSqSum[0]*temp_diff[9]; 
  diff_incr[10] = 0.6324555320336759*nuVtSqSum[1]*temp_diff[17]+0.6324555320336759*nuVtSqSum[2]*temp_diff[10]+0.7071067811865475*nuVtSqSum[0]*temp_diff[10]+0.7071067811865475*nuVtSqSum[1]*temp_diff[6]; 
  diff_incr[11] = 0.4517539514526256*nuVtSqSum[2]*temp_diff[11]+0.7071067811865475*nuVtSqSum[0]*temp_diff[11]+0.632455532033676*nuVtSqSum[1]*temp_diff[4]+0.7071067811865475*nuVtSqSum[2]*temp_diff[2]; 
  diff_incr[12] = 0.6324555320336759*nuVtSqSum[2]*temp_diff[12]+0.7071067811865475*nuVtSqSum[0]*temp_diff[12]+0.7071067811865475*nuVtSqSum[1]*temp_diff[8]; 
  diff_incr[13] = 0.4517539514526256*nuVtSqSum[2]*temp_diff[13]+0.7071067811865475*nuVtSqSum[0]*temp_diff[13]+0.632455532033676*nuVtSqSum[1]*temp_diff[5]+0.7071067811865475*nuVtSqSum[2]*temp_diff[3]; 
  diff_incr[14] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[18]+0.7071067811865475*nuVtSqSum[0]*temp_diff[14]; 
  diff_incr[15] = 0.6324555320336759*nuVtSqSum[2]*temp_diff[15]+0.7071067811865475*nuVtSqSum[0]*temp_diff[15]+0.7071067811865475*nuVtSqSum[1]*temp_diff[9]; 
  diff_incr[16] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[19]+0.7071067811865475*nuVtSqSum[0]*temp_diff[16]; 
  diff_incr[17] = 0.4517539514526256*nuVtSqSum[2]*temp_diff[17]+0.7071067811865475*nuVtSqSum[0]*temp_diff[17]+0.6324555320336759*nuVtSqSum[1]*temp_diff[10]+0.7071067811865475*nuVtSqSum[2]*temp_diff[6]; 
  diff_incr[18] = 0.6324555320336759*nuVtSqSum[2]*temp_diff[18]+0.7071067811865475*nuVtSqSum[0]*temp_diff[18]+0.7071067811865475*nuVtSqSum[1]*temp_diff[14]; 
  diff_incr[19] = 0.6324555320336759*nuVtSqSum[2]*temp_diff[19]+0.7071067811865475*nuVtSqSum[0]*temp_diff[19]+0.7071067811865475*nuVtSqSum[1]*temp_diff[16]; 

  edge_incr[0] = 0.7071067811865475*nuVtSqSum[2]*temp_edge[7]+0.7071067811865475*nuVtSqSum[1]*temp_edge[1]+0.7071067811865475*nuVtSqSum[0]*temp_edge[0]; 
  edge_incr[1] = 0.6324555320336759*nuVtSqSum[1]*temp_edge[7]+0.6324555320336759*temp_edge[1]*nuVtSqSum[2]+0.7071067811865475*nuVtSqSum[0]*temp_edge[1]+0.7071067811865475*temp_edge[0]*nuVtSqSum[1]; 
  edge_incr[2] = 0.7071067811865475*nuVtSqSum[2]*temp_edge[11]+0.7071067811865475*nuVtSqSum[1]*temp_edge[4]+0.7071067811865475*nuVtSqSum[0]*temp_edge[2]; 
  edge_incr[3] = 0.7071067811865475*nuVtSqSum[2]*temp_edge[13]+0.7071067811865475*nuVtSqSum[1]*temp_edge[5]+0.7071067811865475*nuVtSqSum[0]*temp_edge[3]; 
  edge_incr[4] = 0.632455532033676*nuVtSqSum[1]*temp_edge[11]+0.6324555320336759*nuVtSqSum[2]*temp_edge[4]+0.7071067811865475*nuVtSqSum[0]*temp_edge[4]+0.7071067811865475*nuVtSqSum[1]*temp_edge[2]; 
  edge_incr[5] = 0.632455532033676*nuVtSqSum[1]*temp_edge[13]+0.6324555320336759*nuVtSqSum[2]*temp_edge[5]+0.7071067811865475*nuVtSqSum[0]*temp_edge[5]+0.7071067811865475*nuVtSqSum[1]*temp_edge[3]; 
  edge_incr[6] = 0.7071067811865475*nuVtSqSum[2]*temp_edge[17]+0.7071067811865475*nuVtSqSum[1]*temp_edge[10]+0.7071067811865475*nuVtSqSum[0]*temp_edge[6]; 
  edge_incr[7] = 0.4517539514526256*nuVtSqSum[2]*temp_edge[7]+0.7071067811865475*nuVtSqSum[0]*temp_edge[7]+0.7071067811865475*temp_edge[0]*nuVtSqSum[2]+0.6324555320336759*nuVtSqSum[1]*temp_edge[1]; 
  edge_incr[8] = 0.7071067811865475*nuVtSqSum[1]*temp_edge[12]+0.7071067811865475*nuVtSqSum[0]*temp_edge[8]; 
  edge_incr[9] = 0.7071067811865475*nuVtSqSum[1]*temp_edge[15]+0.7071067811865475*nuVtSqSum[0]*temp_edge[9]; 
  edge_incr[10] = 0.6324555320336759*nuVtSqSum[1]*temp_edge[17]+0.6324555320336759*nuVtSqSum[2]*temp_edge[10]+0.7071067811865475*nuVtSqSum[0]*temp_edge[10]+0.7071067811865475*nuVtSqSum[1]*temp_edge[6]; 
  edge_incr[11] = 0.4517539514526256*nuVtSqSum[2]*temp_edge[11]+0.7071067811865475*nuVtSqSum[0]*temp_edge[11]+0.632455532033676*nuVtSqSum[1]*temp_edge[4]+0.7071067811865475*nuVtSqSum[2]*temp_edge[2]; 
  edge_incr[12] = 0.6324555320336759*nuVtSqSum[2]*temp_edge[12]+0.7071067811865475*nuVtSqSum[0]*temp_edge[12]+0.7071067811865475*nuVtSqSum[1]*temp_edge[8]; 
  edge_incr[13] = 0.4517539514526256*nuVtSqSum[2]*temp_edge[13]+0.7071067811865475*nuVtSqSum[0]*temp_edge[13]+0.632455532033676*nuVtSqSum[1]*temp_edge[5]+0.7071067811865475*nuVtSqSum[2]*temp_edge[3]; 
  edge_incr[14] = 0.7071067811865475*nuVtSqSum[1]*temp_edge[18]+0.7071067811865475*nuVtSqSum[0]*temp_edge[14]; 
  edge_incr[15] = 0.6324555320336759*nuVtSqSum[2]*temp_edge[15]+0.7071067811865475*nuVtSqSum[0]*temp_edge[15]+0.7071067811865475*nuVtSqSum[1]*temp_edge[9]; 
  edge_incr[16] = 0.7071067811865475*nuVtSqSum[1]*temp_edge[19]+0.7071067811865475*nuVtSqSum[0]*temp_edge[16]; 
  edge_incr[17] = 0.4517539514526256*nuVtSqSum[2]*temp_edge[17]+0.7071067811865475*nuVtSqSum[0]*temp_edge[17]+0.6324555320336759*nuVtSqSum[1]*temp_edge[10]+0.7071067811865475*nuVtSqSum[2]*temp_edge[6]; 
  edge_incr[18] = 0.6324555320336759*nuVtSqSum[2]*temp_edge[18]+0.7071067811865475*nuVtSqSum[0]*temp_edge[18]+0.7071067811865475*nuVtSqSum[1]*temp_edge[14]; 
  edge_incr[19] = 0.6324555320336759*nuVtSqSum[2]*temp_edge[19]+0.7071067811865475*nuVtSqSum[0]*temp_edge[19]+0.7071067811865475*nuVtSqSum[1]*temp_edge[16]; 

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
} 
