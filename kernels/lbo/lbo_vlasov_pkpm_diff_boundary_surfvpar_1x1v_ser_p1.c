#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH void lbo_vlasov_pkpm_diff_boundary_surfvpar_1x1v_ser_p1(const double *w, const double *dxv, const double *nuVtSq, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[2]:         Cell-center coordinates. 
  // dxv[2]:       Cell spacing. 
  // nuVtSq[2]:    Thermal speeds squared times collisionality. 
  // fSkin/Edge:   Distribution function in cells 
  // out:          Incremented distribution function in cell 
  double rdvSq4 = 4.0/(dxv[1]*dxv[1]); 

  double vol_incr[6] = {0.0}; 
  vol_incr[4] = 4.743416490252569*fSkin[1]*nuVtSq[1]*rdvSq4+4.743416490252569*fSkin[0]*nuVtSq[0]*rdvSq4; 
  vol_incr[5] = 4.743416490252569*fSkin[0]*nuVtSq[1]*rdvSq4+4.743416490252569*nuVtSq[0]*fSkin[1]*rdvSq4; 

  double temp_diff[6] = {0.0}; 
  double temp_edge[6] = {0.0}; 
  double diff_incr[6] = {0.0}; 
  double edge_incr[6] = {0.0}; 

  if (edge == -1) { 

  temp_diff[0] = (-0.6708203932499369*fSkin[4])+0.6708203932499369*fEdge[4]-1.190784930203603*fSkin[2]-1.190784930203603*fEdge[2]-0.9375*fSkin[0]+0.9375*fEdge[0]; 
  temp_diff[1] = (-0.6708203932499369*fSkin[5])+0.6708203932499369*fEdge[5]-1.190784930203603*fSkin[3]-1.190784930203603*fEdge[3]-0.9375*fSkin[1]+0.9375*fEdge[1]; 
  temp_diff[2] = (-1.585502557353661*fSkin[4])+0.7382874503707888*fEdge[4]-2.671875*fSkin[2]-1.453125*fEdge[2]-2.056810333988042*fSkin[0]+1.190784930203603*fEdge[0]; 
  temp_diff[3] = (-1.585502557353661*fSkin[5])+0.7382874503707888*fEdge[5]-2.671875*fSkin[3]-1.453125*fEdge[3]-2.056810333988042*fSkin[1]+1.190784930203603*fEdge[1]; 
  temp_diff[4] = (-3.140625*fSkin[4])-0.140625*fEdge[4]-5.022775277112744*fSkin[2]-0.3025768239224545*fEdge[2]-3.773364712030896*fSkin[0]+0.4192627457812106*fEdge[0]; 
  temp_diff[5] = (-3.140625*fSkin[5])-0.140625*fEdge[5]-5.022775277112744*fSkin[3]-0.3025768239224544*fEdge[3]-3.773364712030894*fSkin[1]+0.4192627457812105*fEdge[1]; 

  temp_edge[2] = 1.936491673103709*fSkin[4]-1.5*fSkin[2]+0.8660254037844386*fSkin[0]; 
  temp_edge[3] = 1.936491673103709*fSkin[5]-1.5*fSkin[3]+0.8660254037844386*fSkin[1]; 
  temp_edge[4] = (-7.5*fSkin[4])+5.809475019311125*fSkin[2]-3.354101966249685*fSkin[0]; 
  temp_edge[5] = (-7.5*fSkin[5])+5.809475019311126*fSkin[3]-3.354101966249684*fSkin[1]; 

  diff_incr[0] = 0.7071067811865475*nuVtSq[1]*temp_diff[1]+0.7071067811865475*nuVtSq[0]*temp_diff[0]; 
  diff_incr[1] = 0.7071067811865475*nuVtSq[0]*temp_diff[1]+0.7071067811865475*temp_diff[0]*nuVtSq[1]; 
  diff_incr[2] = 0.7071067811865475*nuVtSq[1]*temp_diff[3]+0.7071067811865475*nuVtSq[0]*temp_diff[2]; 
  diff_incr[3] = 0.7071067811865475*nuVtSq[0]*temp_diff[3]+0.7071067811865475*nuVtSq[1]*temp_diff[2]; 
  diff_incr[4] = 0.7071067811865475*nuVtSq[1]*temp_diff[5]+0.7071067811865475*nuVtSq[0]*temp_diff[4]; 
  diff_incr[5] = 0.7071067811865475*nuVtSq[0]*temp_diff[5]+0.7071067811865475*nuVtSq[1]*temp_diff[4]; 

  edge_incr[0] = 0.7071067811865475*nuVtSq[1]*temp_edge[1]+0.7071067811865475*nuVtSq[0]*temp_edge[0]; 
  edge_incr[1] = 0.7071067811865475*nuVtSq[0]*temp_edge[1]+0.7071067811865475*temp_edge[0]*nuVtSq[1]; 
  edge_incr[2] = 0.7071067811865475*nuVtSq[1]*temp_edge[3]+0.7071067811865475*nuVtSq[0]*temp_edge[2]; 
  edge_incr[3] = 0.7071067811865475*nuVtSq[0]*temp_edge[3]+0.7071067811865475*nuVtSq[1]*temp_edge[2]; 
  edge_incr[4] = 0.7071067811865475*nuVtSq[1]*temp_edge[5]+0.7071067811865475*nuVtSq[0]*temp_edge[4]; 
  edge_incr[5] = 0.7071067811865475*nuVtSq[0]*temp_edge[5]+0.7071067811865475*nuVtSq[1]*temp_edge[4]; 


  } else { 

  temp_diff[0] = (-0.6708203932499369*fSkin[4])+0.6708203932499369*fEdge[4]+1.190784930203603*fSkin[2]+1.190784930203603*fEdge[2]-0.9375*fSkin[0]+0.9375*fEdge[0]; 
  temp_diff[1] = (-0.6708203932499369*fSkin[5])+0.6708203932499369*fEdge[5]+1.190784930203603*fSkin[3]+1.190784930203603*fEdge[3]-0.9375*fSkin[1]+0.9375*fEdge[1]; 
  temp_diff[2] = 1.585502557353661*fSkin[4]-0.7382874503707888*fEdge[4]-2.671875*fSkin[2]-1.453125*fEdge[2]+2.056810333988042*fSkin[0]-1.190784930203603*fEdge[0]; 
  temp_diff[3] = 1.585502557353661*fSkin[5]-0.7382874503707888*fEdge[5]-2.671875*fSkin[3]-1.453125*fEdge[3]+2.056810333988042*fSkin[1]-1.190784930203603*fEdge[1]; 
  temp_diff[4] = (-3.140625*fSkin[4])-0.140625*fEdge[4]+5.022775277112744*fSkin[2]+0.3025768239224545*fEdge[2]-3.773364712030896*fSkin[0]+0.4192627457812106*fEdge[0]; 
  temp_diff[5] = (-3.140625*fSkin[5])-0.140625*fEdge[5]+5.022775277112744*fSkin[3]+0.3025768239224544*fEdge[3]-3.773364712030894*fSkin[1]+0.4192627457812105*fEdge[1]; 

  temp_edge[2] = (-1.936491673103709*fSkin[4])-1.5*fSkin[2]-0.8660254037844386*fSkin[0]; 
  temp_edge[3] = (-1.936491673103709*fSkin[5])-1.5*fSkin[3]-0.8660254037844386*fSkin[1]; 
  temp_edge[4] = (-7.5*fSkin[4])-5.809475019311125*fSkin[2]-3.354101966249685*fSkin[0]; 
  temp_edge[5] = (-7.5*fSkin[5])-5.809475019311126*fSkin[3]-3.354101966249684*fSkin[1]; 

  diff_incr[0] = 0.7071067811865475*nuVtSq[1]*temp_diff[1]+0.7071067811865475*nuVtSq[0]*temp_diff[0]; 
  diff_incr[1] = 0.7071067811865475*nuVtSq[0]*temp_diff[1]+0.7071067811865475*temp_diff[0]*nuVtSq[1]; 
  diff_incr[2] = 0.7071067811865475*nuVtSq[1]*temp_diff[3]+0.7071067811865475*nuVtSq[0]*temp_diff[2]; 
  diff_incr[3] = 0.7071067811865475*nuVtSq[0]*temp_diff[3]+0.7071067811865475*nuVtSq[1]*temp_diff[2]; 
  diff_incr[4] = 0.7071067811865475*nuVtSq[1]*temp_diff[5]+0.7071067811865475*nuVtSq[0]*temp_diff[4]; 
  diff_incr[5] = 0.7071067811865475*nuVtSq[0]*temp_diff[5]+0.7071067811865475*nuVtSq[1]*temp_diff[4]; 

  edge_incr[0] = 0.7071067811865475*nuVtSq[1]*temp_edge[1]+0.7071067811865475*nuVtSq[0]*temp_edge[0]; 
  edge_incr[1] = 0.7071067811865475*nuVtSq[0]*temp_edge[1]+0.7071067811865475*temp_edge[0]*nuVtSq[1]; 
  edge_incr[2] = 0.7071067811865475*nuVtSq[1]*temp_edge[3]+0.7071067811865475*nuVtSq[0]*temp_edge[2]; 
  edge_incr[3] = 0.7071067811865475*nuVtSq[0]*temp_edge[3]+0.7071067811865475*nuVtSq[1]*temp_edge[2]; 
  edge_incr[4] = 0.7071067811865475*nuVtSq[1]*temp_edge[5]+0.7071067811865475*nuVtSq[0]*temp_edge[4]; 
  edge_incr[5] = 0.7071067811865475*nuVtSq[0]*temp_edge[5]+0.7071067811865475*nuVtSq[1]*temp_edge[4]; 

  } 

  out[0] += edge_incr[0]*rdvSq4+diff_incr[0]*rdvSq4+vol_incr[0]; 
  out[1] += edge_incr[1]*rdvSq4+diff_incr[1]*rdvSq4+vol_incr[1]; 
  out[2] += edge_incr[2]*rdvSq4+diff_incr[2]*rdvSq4+vol_incr[2]; 
  out[3] += edge_incr[3]*rdvSq4+diff_incr[3]*rdvSq4+vol_incr[3]; 
  out[4] += edge_incr[4]*rdvSq4+diff_incr[4]*rdvSq4+vol_incr[4]; 
  out[5] += edge_incr[5]*rdvSq4+diff_incr[5]*rdvSq4+vol_incr[5]; 
} 
