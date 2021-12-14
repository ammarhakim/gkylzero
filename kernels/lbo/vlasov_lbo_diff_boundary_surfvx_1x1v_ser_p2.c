#include <gkyl_vlasov_lbo_kernels.h> 
GKYL_CU_DH void vlasov_lbo_diff_boundary_surfvx_1x1v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[2]:         Cell-center coordinates. 
  // dxv[2]:       Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[3]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[3]: sum of thermal speeds squared time their respective collisionalities. 
  // fSkin/Edge:    Distribution function in cells 
  // out:           Incremented distribution function in cell 
  double rdvSq4 = 4.0/(dxv[1]*dxv[1]); 

  double facDiff[3]; 
  // Expand nuVtSqSum in phase basis.
  facDiff[0] = nuVtSqSum[0]; 
  facDiff[1] = nuVtSqSum[1]; 
  facDiff[2] = nuVtSqSum[2]; 

  double temp_diff[8] = {0.0}; 
  double temp_edge[8] = {0.0}; 
  double diff_incr[8] = {0.0}; 
  double edge_incr[8] = {0.0}; 
  double vol_incr[8] = {0.0}; 

  vol_incr[5] = 4.743416490252569*facDiff[2]*fSkin[4]*rdvSq4+4.743416490252569*fSkin[1]*facDiff[1]*rdvSq4+4.743416490252569*fSkin[0]*facDiff[0]*rdvSq4; 
  vol_incr[7] = 4.242640687119286*facDiff[1]*fSkin[4]*rdvSq4+4.242640687119286*fSkin[1]*facDiff[2]*rdvSq4+4.743416490252569*fSkin[0]*facDiff[1]*rdvSq4+4.743416490252569*facDiff[0]*fSkin[1]*rdvSq4; 

  if (edge == -1) { 

  temp_diff[0] = (-0.6708203932499369*fSkin[5])+0.6708203932499369*fEdge[5]-1.190784930203603*fSkin[2]-1.190784930203603*fEdge[2]-0.9375*fSkin[0]+0.9375*fEdge[0]; 
  temp_diff[1] = (-0.6708203932499369*fSkin[7])+0.6708203932499369*fEdge[7]-1.190784930203603*fSkin[3]-1.190784930203603*fEdge[3]-0.9375*fSkin[1]+0.9375*fEdge[1]; 
  temp_diff[2] = (-1.585502557353661*fSkin[5])+0.7382874503707888*fEdge[5]-2.671875*fSkin[2]-1.453125*fEdge[2]-2.056810333988042*fSkin[0]+1.190784930203603*fEdge[0]; 
  temp_diff[3] = (-1.585502557353661*fSkin[7])+0.7382874503707888*fEdge[7]-2.671875*fSkin[3]-1.453125*fEdge[3]-2.056810333988042*fSkin[1]+1.190784930203603*fEdge[1]; 
  temp_diff[4] = (-1.190784930203603*fSkin[6])-1.190784930203603*fEdge[6]-0.9375*fSkin[4]+0.9375*fEdge[4]; 
  temp_diff[5] = (-3.140625*fSkin[5])-0.140625*fEdge[5]-5.022775277112744*fSkin[2]-0.3025768239224545*fEdge[2]-3.773364712030896*fSkin[0]+0.4192627457812106*fEdge[0]; 
  temp_diff[6] = (-2.671875*fSkin[6])-1.453125*fEdge[6]-2.056810333988042*fSkin[4]+1.190784930203603*fEdge[4]; 
  temp_diff[7] = (-3.140625*fSkin[7])-0.140625*fEdge[7]-5.022775277112744*fSkin[3]-0.3025768239224544*fEdge[3]-3.773364712030894*fSkin[1]+0.4192627457812105*fEdge[1]; 

  temp_edge[2] = 1.936491673103709*fSkin[5]-1.5*fSkin[2]+0.8660254037844386*fSkin[0]; 
  temp_edge[3] = 1.936491673103709*fSkin[7]-1.5*fSkin[3]+0.8660254037844386*fSkin[1]; 
  temp_edge[5] = (-7.5*fSkin[5])+5.809475019311125*fSkin[2]-3.354101966249685*fSkin[0]; 
  temp_edge[6] = 0.8660254037844387*fSkin[4]-1.5*fSkin[6]; 
  temp_edge[7] = (-7.5*fSkin[7])+5.809475019311126*fSkin[3]-3.354101966249684*fSkin[1]; 

  diff_incr[0] = 0.7071067811865475*nuVtSqSum[2]*temp_diff[4]+0.7071067811865475*nuVtSqSum[1]*temp_diff[1]+0.7071067811865475*nuVtSqSum[0]*temp_diff[0]; 
  diff_incr[1] = 0.6324555320336759*nuVtSqSum[1]*temp_diff[4]+0.6324555320336759*temp_diff[1]*nuVtSqSum[2]+0.7071067811865475*nuVtSqSum[0]*temp_diff[1]+0.7071067811865475*temp_diff[0]*nuVtSqSum[1]; 
  diff_incr[2] = 0.7071067811865475*nuVtSqSum[2]*temp_diff[6]+0.7071067811865475*nuVtSqSum[1]*temp_diff[3]+0.7071067811865475*nuVtSqSum[0]*temp_diff[2]; 
  diff_incr[3] = 0.632455532033676*nuVtSqSum[1]*temp_diff[6]+0.6324555320336759*nuVtSqSum[2]*temp_diff[3]+0.7071067811865475*nuVtSqSum[0]*temp_diff[3]+0.7071067811865475*nuVtSqSum[1]*temp_diff[2]; 
  diff_incr[4] = 0.4517539514526256*nuVtSqSum[2]*temp_diff[4]+0.7071067811865475*nuVtSqSum[0]*temp_diff[4]+0.7071067811865475*temp_diff[0]*nuVtSqSum[2]+0.6324555320336759*nuVtSqSum[1]*temp_diff[1]; 
  diff_incr[5] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[7]+0.7071067811865475*nuVtSqSum[0]*temp_diff[5]; 
  diff_incr[6] = 0.4517539514526256*nuVtSqSum[2]*temp_diff[6]+0.7071067811865475*nuVtSqSum[0]*temp_diff[6]+0.632455532033676*nuVtSqSum[1]*temp_diff[3]+0.7071067811865475*nuVtSqSum[2]*temp_diff[2]; 
  diff_incr[7] = 0.6324555320336759*nuVtSqSum[2]*temp_diff[7]+0.7071067811865475*nuVtSqSum[0]*temp_diff[7]+0.7071067811865475*nuVtSqSum[1]*temp_diff[5]; 

  edge_incr[0] = 0.7071067811865475*nuVtSqSum[2]*temp_edge[4]+0.7071067811865475*nuVtSqSum[1]*temp_edge[1]+0.7071067811865475*nuVtSqSum[0]*temp_edge[0]; 
  edge_incr[1] = 0.6324555320336759*nuVtSqSum[1]*temp_edge[4]+0.6324555320336759*temp_edge[1]*nuVtSqSum[2]+0.7071067811865475*nuVtSqSum[0]*temp_edge[1]+0.7071067811865475*temp_edge[0]*nuVtSqSum[1]; 
  edge_incr[2] = 0.7071067811865475*nuVtSqSum[2]*temp_edge[6]+0.7071067811865475*nuVtSqSum[1]*temp_edge[3]+0.7071067811865475*nuVtSqSum[0]*temp_edge[2]; 
  edge_incr[3] = 0.632455532033676*nuVtSqSum[1]*temp_edge[6]+0.6324555320336759*nuVtSqSum[2]*temp_edge[3]+0.7071067811865475*nuVtSqSum[0]*temp_edge[3]+0.7071067811865475*nuVtSqSum[1]*temp_edge[2]; 
  edge_incr[4] = 0.4517539514526256*nuVtSqSum[2]*temp_edge[4]+0.7071067811865475*nuVtSqSum[0]*temp_edge[4]+0.7071067811865475*temp_edge[0]*nuVtSqSum[2]+0.6324555320336759*nuVtSqSum[1]*temp_edge[1]; 
  edge_incr[5] = 0.7071067811865475*nuVtSqSum[1]*temp_edge[7]+0.7071067811865475*nuVtSqSum[0]*temp_edge[5]; 
  edge_incr[6] = 0.4517539514526256*nuVtSqSum[2]*temp_edge[6]+0.7071067811865475*nuVtSqSum[0]*temp_edge[6]+0.632455532033676*nuVtSqSum[1]*temp_edge[3]+0.7071067811865475*nuVtSqSum[2]*temp_edge[2]; 
  edge_incr[7] = 0.6324555320336759*nuVtSqSum[2]*temp_edge[7]+0.7071067811865475*nuVtSqSum[0]*temp_edge[7]+0.7071067811865475*nuVtSqSum[1]*temp_edge[5]; 


  } else { 

  temp_diff[0] = (-0.6708203932499369*fSkin[5])+0.6708203932499369*fEdge[5]-1.190784930203603*fSkin[2]-1.190784930203603*fEdge[2]-0.9375*fSkin[0]+0.9375*fEdge[0]; 
  temp_diff[1] = (-0.6708203932499369*fSkin[7])+0.6708203932499369*fEdge[7]-1.190784930203603*fSkin[3]-1.190784930203603*fEdge[3]-0.9375*fSkin[1]+0.9375*fEdge[1]; 
  temp_diff[2] = (-1.585502557353661*fSkin[5])+0.7382874503707888*fEdge[5]-2.671875*fSkin[2]-1.453125*fEdge[2]-2.056810333988042*fSkin[0]+1.190784930203603*fEdge[0]; 
  temp_diff[3] = (-1.585502557353661*fSkin[7])+0.7382874503707888*fEdge[7]-2.671875*fSkin[3]-1.453125*fEdge[3]-2.056810333988042*fSkin[1]+1.190784930203603*fEdge[1]; 
  temp_diff[4] = (-1.190784930203603*fSkin[6])-1.190784930203603*fEdge[6]-0.9375*fSkin[4]+0.9375*fEdge[4]; 
  temp_diff[5] = (-3.140625*fSkin[5])-0.140625*fEdge[5]-5.022775277112744*fSkin[2]-0.3025768239224545*fEdge[2]-3.773364712030896*fSkin[0]+0.4192627457812106*fEdge[0]; 
  temp_diff[6] = (-2.671875*fSkin[6])-1.453125*fEdge[6]-2.056810333988042*fSkin[4]+1.190784930203603*fEdge[4]; 
  temp_diff[7] = (-3.140625*fSkin[7])-0.140625*fEdge[7]-5.022775277112744*fSkin[3]-0.3025768239224544*fEdge[3]-3.773364712030894*fSkin[1]+0.4192627457812105*fEdge[1]; 

  temp_edge[2] = 1.936491673103709*fSkin[5]-1.5*fSkin[2]+0.8660254037844386*fSkin[0]; 
  temp_edge[3] = 1.936491673103709*fSkin[7]-1.5*fSkin[3]+0.8660254037844386*fSkin[1]; 
  temp_edge[5] = (-7.5*fSkin[5])+5.809475019311125*fSkin[2]-3.354101966249685*fSkin[0]; 
  temp_edge[6] = 0.8660254037844387*fSkin[4]-1.5*fSkin[6]; 
  temp_edge[7] = (-7.5*fSkin[7])+5.809475019311126*fSkin[3]-3.354101966249684*fSkin[1]; 

  diff_incr[0] = 0.7071067811865475*nuVtSqSum[2]*temp_diff[4]+0.7071067811865475*nuVtSqSum[1]*temp_diff[1]+0.7071067811865475*nuVtSqSum[0]*temp_diff[0]; 
  diff_incr[1] = 0.6324555320336759*nuVtSqSum[1]*temp_diff[4]+0.6324555320336759*temp_diff[1]*nuVtSqSum[2]+0.7071067811865475*nuVtSqSum[0]*temp_diff[1]+0.7071067811865475*temp_diff[0]*nuVtSqSum[1]; 
  diff_incr[2] = 0.7071067811865475*nuVtSqSum[2]*temp_diff[6]+0.7071067811865475*nuVtSqSum[1]*temp_diff[3]+0.7071067811865475*nuVtSqSum[0]*temp_diff[2]; 
  diff_incr[3] = 0.632455532033676*nuVtSqSum[1]*temp_diff[6]+0.6324555320336759*nuVtSqSum[2]*temp_diff[3]+0.7071067811865475*nuVtSqSum[0]*temp_diff[3]+0.7071067811865475*nuVtSqSum[1]*temp_diff[2]; 
  diff_incr[4] = 0.4517539514526256*nuVtSqSum[2]*temp_diff[4]+0.7071067811865475*nuVtSqSum[0]*temp_diff[4]+0.7071067811865475*temp_diff[0]*nuVtSqSum[2]+0.6324555320336759*nuVtSqSum[1]*temp_diff[1]; 
  diff_incr[5] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[7]+0.7071067811865475*nuVtSqSum[0]*temp_diff[5]; 
  diff_incr[6] = 0.4517539514526256*nuVtSqSum[2]*temp_diff[6]+0.7071067811865475*nuVtSqSum[0]*temp_diff[6]+0.632455532033676*nuVtSqSum[1]*temp_diff[3]+0.7071067811865475*nuVtSqSum[2]*temp_diff[2]; 
  diff_incr[7] = 0.6324555320336759*nuVtSqSum[2]*temp_diff[7]+0.7071067811865475*nuVtSqSum[0]*temp_diff[7]+0.7071067811865475*nuVtSqSum[1]*temp_diff[5]; 

  edge_incr[0] = 0.7071067811865475*nuVtSqSum[2]*temp_edge[4]+0.7071067811865475*nuVtSqSum[1]*temp_edge[1]+0.7071067811865475*nuVtSqSum[0]*temp_edge[0]; 
  edge_incr[1] = 0.6324555320336759*nuVtSqSum[1]*temp_edge[4]+0.6324555320336759*temp_edge[1]*nuVtSqSum[2]+0.7071067811865475*nuVtSqSum[0]*temp_edge[1]+0.7071067811865475*temp_edge[0]*nuVtSqSum[1]; 
  edge_incr[2] = 0.7071067811865475*nuVtSqSum[2]*temp_edge[6]+0.7071067811865475*nuVtSqSum[1]*temp_edge[3]+0.7071067811865475*nuVtSqSum[0]*temp_edge[2]; 
  edge_incr[3] = 0.632455532033676*nuVtSqSum[1]*temp_edge[6]+0.6324555320336759*nuVtSqSum[2]*temp_edge[3]+0.7071067811865475*nuVtSqSum[0]*temp_edge[3]+0.7071067811865475*nuVtSqSum[1]*temp_edge[2]; 
  edge_incr[4] = 0.4517539514526256*nuVtSqSum[2]*temp_edge[4]+0.7071067811865475*nuVtSqSum[0]*temp_edge[4]+0.7071067811865475*temp_edge[0]*nuVtSqSum[2]+0.6324555320336759*nuVtSqSum[1]*temp_edge[1]; 
  edge_incr[5] = 0.7071067811865475*nuVtSqSum[1]*temp_edge[7]+0.7071067811865475*nuVtSqSum[0]*temp_edge[5]; 
  edge_incr[6] = 0.4517539514526256*nuVtSqSum[2]*temp_edge[6]+0.7071067811865475*nuVtSqSum[0]*temp_edge[6]+0.632455532033676*nuVtSqSum[1]*temp_edge[3]+0.7071067811865475*nuVtSqSum[2]*temp_edge[2]; 
  edge_incr[7] = 0.6324555320336759*nuVtSqSum[2]*temp_edge[7]+0.7071067811865475*nuVtSqSum[0]*temp_edge[7]+0.7071067811865475*nuVtSqSum[1]*temp_edge[5]; 

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
