#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH void lbo_vlasov_diff_boundary_surfvx_1x2v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[3]:         Cell-center coordinates. 
  // dxv[3]:       Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[4]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[2]: sum of thermal speeds squared time their respective collisionalities. 
  // fSkin/Edge:    Distribution function in cells 
  // out:           Incremented distribution function in cell 
  double rdvSq4 = 4.0/(dxv[1]*dxv[1]); 

  double vol_incr[16] = {0.0}; 
  vol_incr[8] = 4.743416490252569*fSkin[1]*nuVtSqSum[1]*rdvSq4+4.743416490252569*fSkin[0]*nuVtSqSum[0]*rdvSq4; 
  vol_incr[9] = 4.743416490252569*fSkin[0]*nuVtSqSum[1]*rdvSq4+4.743416490252569*nuVtSqSum[0]*fSkin[1]*rdvSq4; 
  vol_incr[10] = 4.743416490252569*nuVtSqSum[1]*fSkin[5]*rdvSq4+4.743416490252569*nuVtSqSum[0]*fSkin[3]*rdvSq4; 
  vol_incr[11] = 4.743416490252569*nuVtSqSum[0]*fSkin[5]*rdvSq4+4.743416490252569*nuVtSqSum[1]*fSkin[3]*rdvSq4; 

  double temp_diff[16] = {0.0}; 
  double temp_edge[16] = {0.0}; 
  double diff_incr[16] = {0.0}; 
  double edge_incr[16] = {0.0}; 

  if (edge == -1) { 

  temp_diff[0] = (-0.6708203932499369*fSkin[8])+0.6708203932499369*fEdge[8]-1.190784930203603*fSkin[2]-1.190784930203603*fEdge[2]-0.9375*fSkin[0]+0.9375*fEdge[0]; 
  temp_diff[1] = (-0.6708203932499369*fSkin[9])+0.6708203932499369*fEdge[9]-1.190784930203603*fSkin[4]-1.190784930203603*fEdge[4]-0.9375*fSkin[1]+0.9375*fEdge[1]; 
  temp_diff[2] = (-1.585502557353661*fSkin[8])+0.7382874503707888*fEdge[8]-2.671875*fSkin[2]-1.453125*fEdge[2]-2.056810333988042*fSkin[0]+1.190784930203603*fEdge[0]; 
  temp_diff[3] = (-0.6708203932499369*fSkin[10])+0.6708203932499369*fEdge[10]-1.190784930203603*fSkin[6]-1.190784930203603*fEdge[6]-0.9375*fSkin[3]+0.9375*fEdge[3]; 
  temp_diff[4] = (-1.585502557353661*fSkin[9])+0.7382874503707888*fEdge[9]-2.671875*fSkin[4]-1.453125*fEdge[4]-2.056810333988042*fSkin[1]+1.190784930203603*fEdge[1]; 
  temp_diff[5] = (-0.6708203932499369*fSkin[11])+0.6708203932499369*fEdge[11]-1.190784930203603*fSkin[7]-1.190784930203603*fEdge[7]-0.9375*fSkin[5]+0.9375*fEdge[5]; 
  temp_diff[6] = (-1.585502557353661*fSkin[10])+0.7382874503707888*fEdge[10]-2.671875*fSkin[6]-1.453125*fEdge[6]-2.056810333988042*fSkin[3]+1.190784930203603*fEdge[3]; 
  temp_diff[7] = (-1.585502557353661*fSkin[11])+0.7382874503707888*fEdge[11]-2.671875*fSkin[7]-1.453125*fEdge[7]-2.056810333988042*fSkin[5]+1.190784930203603*fEdge[5]; 
  temp_diff[8] = (-3.140625*fSkin[8])-0.140625*fEdge[8]-5.022775277112744*fSkin[2]-0.3025768239224545*fEdge[2]-3.773364712030896*fSkin[0]+0.4192627457812106*fEdge[0]; 
  temp_diff[9] = (-3.140625*fSkin[9])-0.140625*fEdge[9]-5.022775277112744*fSkin[4]-0.3025768239224544*fEdge[4]-3.773364712030894*fSkin[1]+0.4192627457812105*fEdge[1]; 
  temp_diff[10] = (-3.140625*fSkin[10])-0.140625*fEdge[10]-5.022775277112744*fSkin[6]-0.3025768239224544*fEdge[6]-3.773364712030894*fSkin[3]+0.4192627457812105*fEdge[3]; 
  temp_diff[11] = (-3.140625*fSkin[11])-0.140625*fEdge[11]-5.022775277112744*fSkin[7]-0.3025768239224545*fEdge[7]-3.773364712030896*fSkin[5]+0.4192627457812106*fEdge[5]; 
  temp_diff[12] = (-1.190784930203603*fSkin[14])-1.190784930203603*fEdge[14]-0.9375*fSkin[12]+0.9375*fEdge[12]; 
  temp_diff[13] = (-1.190784930203603*fSkin[15])-1.190784930203603*fEdge[15]-0.9375*fSkin[13]+0.9375*fEdge[13]; 
  temp_diff[14] = (-2.671875*fSkin[14])-1.453125*fEdge[14]-2.056810333988042*fSkin[12]+1.190784930203603*fEdge[12]; 
  temp_diff[15] = (-2.671875*fSkin[15])-1.453125*fEdge[15]-2.056810333988042*fSkin[13]+1.190784930203603*fEdge[13]; 

  temp_edge[2] = 1.936491673103709*fSkin[8]-1.5*fSkin[2]+0.8660254037844386*fSkin[0]; 
  temp_edge[4] = 1.936491673103709*fSkin[9]-1.5*fSkin[4]+0.8660254037844386*fSkin[1]; 
  temp_edge[6] = 1.936491673103709*fSkin[10]-1.5*fSkin[6]+0.8660254037844386*fSkin[3]; 
  temp_edge[7] = 1.936491673103709*fSkin[11]-1.5*fSkin[7]+0.8660254037844386*fSkin[5]; 
  temp_edge[8] = (-7.5*fSkin[8])+5.809475019311125*fSkin[2]-3.354101966249685*fSkin[0]; 
  temp_edge[9] = (-7.5*fSkin[9])+5.809475019311126*fSkin[4]-3.354101966249684*fSkin[1]; 
  temp_edge[10] = (-7.5*fSkin[10])+5.809475019311126*fSkin[6]-3.354101966249684*fSkin[3]; 
  temp_edge[11] = (-7.5*fSkin[11])+5.809475019311125*fSkin[7]-3.354101966249685*fSkin[5]; 
  temp_edge[14] = 0.8660254037844387*fSkin[12]-1.5*fSkin[14]; 
  temp_edge[15] = 0.8660254037844387*fSkin[13]-1.5*fSkin[15]; 

  diff_incr[0] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[1]+0.7071067811865475*nuVtSqSum[0]*temp_diff[0]; 
  diff_incr[1] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[1]+0.7071067811865475*temp_diff[0]*nuVtSqSum[1]; 
  diff_incr[2] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[4]+0.7071067811865475*nuVtSqSum[0]*temp_diff[2]; 
  diff_incr[3] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[5]+0.7071067811865475*nuVtSqSum[0]*temp_diff[3]; 
  diff_incr[4] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[4]+0.7071067811865475*nuVtSqSum[1]*temp_diff[2]; 
  diff_incr[5] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[5]+0.7071067811865475*nuVtSqSum[1]*temp_diff[3]; 
  diff_incr[6] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[7]+0.7071067811865475*nuVtSqSum[0]*temp_diff[6]; 
  diff_incr[7] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[7]+0.7071067811865475*nuVtSqSum[1]*temp_diff[6]; 
  diff_incr[8] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[9]+0.7071067811865475*nuVtSqSum[0]*temp_diff[8]; 
  diff_incr[9] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[9]+0.7071067811865475*nuVtSqSum[1]*temp_diff[8]; 
  diff_incr[10] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[11]+0.7071067811865475*nuVtSqSum[0]*temp_diff[10]; 
  diff_incr[11] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[11]+0.7071067811865475*nuVtSqSum[1]*temp_diff[10]; 
  diff_incr[12] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[13]+0.7071067811865475*nuVtSqSum[0]*temp_diff[12]; 
  diff_incr[13] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[13]+0.7071067811865475*nuVtSqSum[1]*temp_diff[12]; 
  diff_incr[14] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[15]+0.7071067811865475*nuVtSqSum[0]*temp_diff[14]; 
  diff_incr[15] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[15]+0.7071067811865475*nuVtSqSum[1]*temp_diff[14]; 

  edge_incr[0] = 0.7071067811865475*nuVtSqSum[1]*temp_edge[1]+0.7071067811865475*nuVtSqSum[0]*temp_edge[0]; 
  edge_incr[1] = 0.7071067811865475*nuVtSqSum[0]*temp_edge[1]+0.7071067811865475*temp_edge[0]*nuVtSqSum[1]; 
  edge_incr[2] = 0.7071067811865475*nuVtSqSum[1]*temp_edge[4]+0.7071067811865475*nuVtSqSum[0]*temp_edge[2]; 
  edge_incr[3] = 0.7071067811865475*nuVtSqSum[1]*temp_edge[5]+0.7071067811865475*nuVtSqSum[0]*temp_edge[3]; 
  edge_incr[4] = 0.7071067811865475*nuVtSqSum[0]*temp_edge[4]+0.7071067811865475*nuVtSqSum[1]*temp_edge[2]; 
  edge_incr[5] = 0.7071067811865475*nuVtSqSum[0]*temp_edge[5]+0.7071067811865475*nuVtSqSum[1]*temp_edge[3]; 
  edge_incr[6] = 0.7071067811865475*nuVtSqSum[1]*temp_edge[7]+0.7071067811865475*nuVtSqSum[0]*temp_edge[6]; 
  edge_incr[7] = 0.7071067811865475*nuVtSqSum[0]*temp_edge[7]+0.7071067811865475*nuVtSqSum[1]*temp_edge[6]; 
  edge_incr[8] = 0.7071067811865475*nuVtSqSum[1]*temp_edge[9]+0.7071067811865475*nuVtSqSum[0]*temp_edge[8]; 
  edge_incr[9] = 0.7071067811865475*nuVtSqSum[0]*temp_edge[9]+0.7071067811865475*nuVtSqSum[1]*temp_edge[8]; 
  edge_incr[10] = 0.7071067811865475*nuVtSqSum[1]*temp_edge[11]+0.7071067811865475*nuVtSqSum[0]*temp_edge[10]; 
  edge_incr[11] = 0.7071067811865475*nuVtSqSum[0]*temp_edge[11]+0.7071067811865475*nuVtSqSum[1]*temp_edge[10]; 
  edge_incr[12] = 0.7071067811865475*nuVtSqSum[1]*temp_edge[13]+0.7071067811865475*nuVtSqSum[0]*temp_edge[12]; 
  edge_incr[13] = 0.7071067811865475*nuVtSqSum[0]*temp_edge[13]+0.7071067811865475*nuVtSqSum[1]*temp_edge[12]; 
  edge_incr[14] = 0.7071067811865475*nuVtSqSum[1]*temp_edge[15]+0.7071067811865475*nuVtSqSum[0]*temp_edge[14]; 
  edge_incr[15] = 0.7071067811865475*nuVtSqSum[0]*temp_edge[15]+0.7071067811865475*nuVtSqSum[1]*temp_edge[14]; 


  } else { 

  temp_diff[0] = (-0.6708203932499369*fSkin[8])+0.6708203932499369*fEdge[8]+1.190784930203603*fSkin[2]+1.190784930203603*fEdge[2]-0.9375*fSkin[0]+0.9375*fEdge[0]; 
  temp_diff[1] = (-0.6708203932499369*fSkin[9])+0.6708203932499369*fEdge[9]+1.190784930203603*fSkin[4]+1.190784930203603*fEdge[4]-0.9375*fSkin[1]+0.9375*fEdge[1]; 
  temp_diff[2] = 1.585502557353661*fSkin[8]-0.7382874503707888*fEdge[8]-2.671875*fSkin[2]-1.453125*fEdge[2]+2.056810333988042*fSkin[0]-1.190784930203603*fEdge[0]; 
  temp_diff[3] = (-0.6708203932499369*fSkin[10])+0.6708203932499369*fEdge[10]+1.190784930203603*fSkin[6]+1.190784930203603*fEdge[6]-0.9375*fSkin[3]+0.9375*fEdge[3]; 
  temp_diff[4] = 1.585502557353661*fSkin[9]-0.7382874503707888*fEdge[9]-2.671875*fSkin[4]-1.453125*fEdge[4]+2.056810333988042*fSkin[1]-1.190784930203603*fEdge[1]; 
  temp_diff[5] = (-0.6708203932499369*fSkin[11])+0.6708203932499369*fEdge[11]+1.190784930203603*fSkin[7]+1.190784930203603*fEdge[7]-0.9375*fSkin[5]+0.9375*fEdge[5]; 
  temp_diff[6] = 1.585502557353661*fSkin[10]-0.7382874503707888*fEdge[10]-2.671875*fSkin[6]-1.453125*fEdge[6]+2.056810333988042*fSkin[3]-1.190784930203603*fEdge[3]; 
  temp_diff[7] = 1.585502557353661*fSkin[11]-0.7382874503707888*fEdge[11]-2.671875*fSkin[7]-1.453125*fEdge[7]+2.056810333988042*fSkin[5]-1.190784930203603*fEdge[5]; 
  temp_diff[8] = (-3.140625*fSkin[8])-0.140625*fEdge[8]+5.022775277112744*fSkin[2]+0.3025768239224545*fEdge[2]-3.773364712030896*fSkin[0]+0.4192627457812106*fEdge[0]; 
  temp_diff[9] = (-3.140625*fSkin[9])-0.140625*fEdge[9]+5.022775277112744*fSkin[4]+0.3025768239224544*fEdge[4]-3.773364712030894*fSkin[1]+0.4192627457812105*fEdge[1]; 
  temp_diff[10] = (-3.140625*fSkin[10])-0.140625*fEdge[10]+5.022775277112744*fSkin[6]+0.3025768239224544*fEdge[6]-3.773364712030894*fSkin[3]+0.4192627457812105*fEdge[3]; 
  temp_diff[11] = (-3.140625*fSkin[11])-0.140625*fEdge[11]+5.022775277112744*fSkin[7]+0.3025768239224545*fEdge[7]-3.773364712030896*fSkin[5]+0.4192627457812106*fEdge[5]; 
  temp_diff[12] = 1.190784930203603*fSkin[14]+1.190784930203603*fEdge[14]-0.9375*fSkin[12]+0.9375*fEdge[12]; 
  temp_diff[13] = 1.190784930203603*fSkin[15]+1.190784930203603*fEdge[15]-0.9375*fSkin[13]+0.9375*fEdge[13]; 
  temp_diff[14] = (-2.671875*fSkin[14])-1.453125*fEdge[14]+2.056810333988042*fSkin[12]-1.190784930203603*fEdge[12]; 
  temp_diff[15] = (-2.671875*fSkin[15])-1.453125*fEdge[15]+2.056810333988042*fSkin[13]-1.190784930203603*fEdge[13]; 

  temp_edge[2] = (-1.936491673103709*fSkin[8])-1.5*fSkin[2]-0.8660254037844386*fSkin[0]; 
  temp_edge[4] = (-1.936491673103709*fSkin[9])-1.5*fSkin[4]-0.8660254037844386*fSkin[1]; 
  temp_edge[6] = (-1.936491673103709*fSkin[10])-1.5*fSkin[6]-0.8660254037844386*fSkin[3]; 
  temp_edge[7] = (-1.936491673103709*fSkin[11])-1.5*fSkin[7]-0.8660254037844386*fSkin[5]; 
  temp_edge[8] = (-7.5*fSkin[8])-5.809475019311125*fSkin[2]-3.354101966249685*fSkin[0]; 
  temp_edge[9] = (-7.5*fSkin[9])-5.809475019311126*fSkin[4]-3.354101966249684*fSkin[1]; 
  temp_edge[10] = (-7.5*fSkin[10])-5.809475019311126*fSkin[6]-3.354101966249684*fSkin[3]; 
  temp_edge[11] = (-7.5*fSkin[11])-5.809475019311125*fSkin[7]-3.354101966249685*fSkin[5]; 
  temp_edge[14] = (-1.5*fSkin[14])-0.8660254037844387*fSkin[12]; 
  temp_edge[15] = (-1.5*fSkin[15])-0.8660254037844387*fSkin[13]; 

  diff_incr[0] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[1]+0.7071067811865475*nuVtSqSum[0]*temp_diff[0]; 
  diff_incr[1] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[1]+0.7071067811865475*temp_diff[0]*nuVtSqSum[1]; 
  diff_incr[2] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[4]+0.7071067811865475*nuVtSqSum[0]*temp_diff[2]; 
  diff_incr[3] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[5]+0.7071067811865475*nuVtSqSum[0]*temp_diff[3]; 
  diff_incr[4] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[4]+0.7071067811865475*nuVtSqSum[1]*temp_diff[2]; 
  diff_incr[5] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[5]+0.7071067811865475*nuVtSqSum[1]*temp_diff[3]; 
  diff_incr[6] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[7]+0.7071067811865475*nuVtSqSum[0]*temp_diff[6]; 
  diff_incr[7] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[7]+0.7071067811865475*nuVtSqSum[1]*temp_diff[6]; 
  diff_incr[8] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[9]+0.7071067811865475*nuVtSqSum[0]*temp_diff[8]; 
  diff_incr[9] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[9]+0.7071067811865475*nuVtSqSum[1]*temp_diff[8]; 
  diff_incr[10] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[11]+0.7071067811865475*nuVtSqSum[0]*temp_diff[10]; 
  diff_incr[11] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[11]+0.7071067811865475*nuVtSqSum[1]*temp_diff[10]; 
  diff_incr[12] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[13]+0.7071067811865475*nuVtSqSum[0]*temp_diff[12]; 
  diff_incr[13] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[13]+0.7071067811865475*nuVtSqSum[1]*temp_diff[12]; 
  diff_incr[14] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[15]+0.7071067811865475*nuVtSqSum[0]*temp_diff[14]; 
  diff_incr[15] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[15]+0.7071067811865475*nuVtSqSum[1]*temp_diff[14]; 

  edge_incr[0] = 0.7071067811865475*nuVtSqSum[1]*temp_edge[1]+0.7071067811865475*nuVtSqSum[0]*temp_edge[0]; 
  edge_incr[1] = 0.7071067811865475*nuVtSqSum[0]*temp_edge[1]+0.7071067811865475*temp_edge[0]*nuVtSqSum[1]; 
  edge_incr[2] = 0.7071067811865475*nuVtSqSum[1]*temp_edge[4]+0.7071067811865475*nuVtSqSum[0]*temp_edge[2]; 
  edge_incr[3] = 0.7071067811865475*nuVtSqSum[1]*temp_edge[5]+0.7071067811865475*nuVtSqSum[0]*temp_edge[3]; 
  edge_incr[4] = 0.7071067811865475*nuVtSqSum[0]*temp_edge[4]+0.7071067811865475*nuVtSqSum[1]*temp_edge[2]; 
  edge_incr[5] = 0.7071067811865475*nuVtSqSum[0]*temp_edge[5]+0.7071067811865475*nuVtSqSum[1]*temp_edge[3]; 
  edge_incr[6] = 0.7071067811865475*nuVtSqSum[1]*temp_edge[7]+0.7071067811865475*nuVtSqSum[0]*temp_edge[6]; 
  edge_incr[7] = 0.7071067811865475*nuVtSqSum[0]*temp_edge[7]+0.7071067811865475*nuVtSqSum[1]*temp_edge[6]; 
  edge_incr[8] = 0.7071067811865475*nuVtSqSum[1]*temp_edge[9]+0.7071067811865475*nuVtSqSum[0]*temp_edge[8]; 
  edge_incr[9] = 0.7071067811865475*nuVtSqSum[0]*temp_edge[9]+0.7071067811865475*nuVtSqSum[1]*temp_edge[8]; 
  edge_incr[10] = 0.7071067811865475*nuVtSqSum[1]*temp_edge[11]+0.7071067811865475*nuVtSqSum[0]*temp_edge[10]; 
  edge_incr[11] = 0.7071067811865475*nuVtSqSum[0]*temp_edge[11]+0.7071067811865475*nuVtSqSum[1]*temp_edge[10]; 
  edge_incr[12] = 0.7071067811865475*nuVtSqSum[1]*temp_edge[13]+0.7071067811865475*nuVtSqSum[0]*temp_edge[12]; 
  edge_incr[13] = 0.7071067811865475*nuVtSqSum[0]*temp_edge[13]+0.7071067811865475*nuVtSqSum[1]*temp_edge[12]; 
  edge_incr[14] = 0.7071067811865475*nuVtSqSum[1]*temp_edge[15]+0.7071067811865475*nuVtSqSum[0]*temp_edge[14]; 
  edge_incr[15] = 0.7071067811865475*nuVtSqSum[0]*temp_edge[15]+0.7071067811865475*nuVtSqSum[1]*temp_edge[14]; 

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
