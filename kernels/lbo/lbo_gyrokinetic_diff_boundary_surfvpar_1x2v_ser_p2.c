#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH void lbo_gyrokinetic_diff_boundary_surfvpar_1x2v_ser_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  // w[3]:         Cell-center coordinates. 
  // dxv[3]:       Cell spacing. 
  // m_:           species mass.
  // bmag_inv:     1/(magnetic field magnitude). 
  // nuSum:        collisionalities added (self and cross species collisionalities). 
  // nuUSum[6]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[3]: sum of thermal speeds squared time their respective collisionalities. 
  // fskin/edge:   Distribution function in cells 
  // out:          Incremented distribution function in cell 
  double rdvSq4 = 4.0/(dxv[1]*dxv[1]); 

  double vol_incr[20] = {0.0};
  vol_incr[8] = 4.743416490252569*nuVtSqSum[2]*fskin[7]*rdvSq4+4.743416490252569*fskin[1]*nuVtSqSum[1]*rdvSq4+4.743416490252569*fskin[0]*nuVtSqSum[0]*rdvSq4; 
  vol_incr[12] = 4.242640687119286*nuVtSqSum[1]*fskin[7]*rdvSq4+4.242640687119286*fskin[1]*nuVtSqSum[2]*rdvSq4+4.743416490252569*fskin[0]*nuVtSqSum[1]*rdvSq4+4.743416490252569*nuVtSqSum[0]*fskin[1]*rdvSq4; 
  vol_incr[14] = 4.743416490252569*nuVtSqSum[2]*fskin[13]*rdvSq4+4.743416490252569*nuVtSqSum[1]*fskin[5]*rdvSq4+4.743416490252569*nuVtSqSum[0]*fskin[3]*rdvSq4; 
  vol_incr[18] = 4.242640687119286*nuVtSqSum[1]*fskin[13]*rdvSq4+4.242640687119286*nuVtSqSum[2]*fskin[5]*rdvSq4+4.743416490252569*nuVtSqSum[0]*fskin[5]*rdvSq4+4.743416490252569*nuVtSqSum[1]*fskin[3]*rdvSq4; 

  double temp_diff[20] = {0.0}; 
  double temp_edge[20] = {0.0}; 
  double diff_incr[20] = {0.0}; 
  double edge_incr[20] = {0.0}; 

  if (edge == -1) { 

  temp_diff[0] = (-0.6708203932499369*fskin[8])+0.6708203932499369*fedge[8]-1.190784930203603*fskin[2]-1.190784930203603*fedge[2]-0.9375*fskin[0]+0.9375*fedge[0]; 
  temp_diff[1] = (-0.6708203932499369*fskin[12])+0.6708203932499369*fedge[12]-1.190784930203603*fskin[4]-1.190784930203603*fedge[4]-0.9375*fskin[1]+0.9375*fedge[1]; 
  temp_diff[2] = (-1.585502557353661*fskin[8])+0.7382874503707888*fedge[8]-2.671875*fskin[2]-1.453125*fedge[2]-2.056810333988042*fskin[0]+1.190784930203603*fedge[0]; 
  temp_diff[3] = (-0.6708203932499369*fskin[14])+0.6708203932499369*fedge[14]-1.190784930203603*fskin[6]-1.190784930203603*fedge[6]-0.9375*fskin[3]+0.9375*fedge[3]; 
  temp_diff[4] = (-1.585502557353661*fskin[12])+0.7382874503707888*fedge[12]-2.671875*fskin[4]-1.453125*fedge[4]-2.056810333988042*fskin[1]+1.190784930203603*fedge[1]; 
  temp_diff[5] = (-0.6708203932499369*fskin[18])+0.6708203932499369*fedge[18]-1.190784930203603*fskin[10]-1.190784930203603*fedge[10]-0.9375*fskin[5]+0.9375*fedge[5]; 
  temp_diff[6] = (-1.585502557353661*fskin[14])+0.7382874503707888*fedge[14]-2.671875*fskin[6]-1.453125*fedge[6]-2.056810333988042*fskin[3]+1.190784930203603*fedge[3]; 
  temp_diff[7] = (-1.190784930203603*fskin[11])-1.190784930203603*fedge[11]-0.9375*fskin[7]+0.9375*fedge[7]; 
  temp_diff[8] = (-3.140625*fskin[8])-0.140625*fedge[8]-5.022775277112744*fskin[2]-0.3025768239224545*fedge[2]-3.773364712030896*fskin[0]+0.4192627457812106*fedge[0]; 
  temp_diff[9] = (-1.190784930203603*fskin[16])-1.190784930203603*fedge[16]-0.9375*fskin[9]+0.9375*fedge[9]; 
  temp_diff[10] = (-1.585502557353661*fskin[18])+0.7382874503707888*fedge[18]-2.671875*fskin[10]-1.453125*fedge[10]-2.056810333988042*fskin[5]+1.190784930203603*fedge[5]; 
  temp_diff[11] = (-2.671875*fskin[11])-1.453125*fedge[11]-2.056810333988042*fskin[7]+1.190784930203603*fedge[7]; 
  temp_diff[12] = (-3.140625*fskin[12])-0.140625*fedge[12]-5.022775277112744*fskin[4]-0.3025768239224544*fedge[4]-3.773364712030894*fskin[1]+0.4192627457812105*fedge[1]; 
  temp_diff[13] = (-1.190784930203603*fskin[17])-1.190784930203603*fedge[17]-0.9375*fskin[13]+0.9375*fedge[13]; 
  temp_diff[14] = (-3.140625*fskin[14])-0.140625*fedge[14]-5.022775277112744*fskin[6]-0.3025768239224544*fedge[6]-3.773364712030894*fskin[3]+0.4192627457812105*fedge[3]; 
  temp_diff[15] = (-1.190784930203603*fskin[19])-1.190784930203603*fedge[19]-0.9375*fskin[15]+0.9375*fedge[15]; 
  temp_diff[16] = (-2.671875*fskin[16])-1.453125*fedge[16]-2.056810333988042*fskin[9]+1.190784930203603*fedge[9]; 
  temp_diff[17] = (-2.671875*fskin[17])-1.453125*fedge[17]-2.056810333988042*fskin[13]+1.190784930203603*fedge[13]; 
  temp_diff[18] = (-3.140625*fskin[18])-0.140625*fedge[18]-5.022775277112744*fskin[10]-0.3025768239224545*fedge[10]-3.773364712030896*fskin[5]+0.4192627457812106*fedge[5]; 
  temp_diff[19] = (-2.671875*fskin[19])-1.453125*fedge[19]-2.056810333988042*fskin[15]+1.190784930203603*fedge[15]; 

  temp_edge[2] = 1.936491673103709*fskin[8]-1.5*fskin[2]+0.8660254037844386*fskin[0]; 
  temp_edge[4] = 1.936491673103709*fskin[12]-1.5*fskin[4]+0.8660254037844386*fskin[1]; 
  temp_edge[6] = 1.936491673103709*fskin[14]-1.5*fskin[6]+0.8660254037844386*fskin[3]; 
  temp_edge[8] = (-7.5*fskin[8])+5.809475019311125*fskin[2]-3.354101966249685*fskin[0]; 
  temp_edge[10] = 1.936491673103709*fskin[18]-1.5*fskin[10]+0.8660254037844386*fskin[5]; 
  temp_edge[11] = 0.8660254037844387*fskin[7]-1.5*fskin[11]; 
  temp_edge[12] = (-7.5*fskin[12])+5.809475019311126*fskin[4]-3.354101966249684*fskin[1]; 
  temp_edge[14] = (-7.5*fskin[14])+5.809475019311126*fskin[6]-3.354101966249684*fskin[3]; 
  temp_edge[16] = 0.8660254037844387*fskin[9]-1.5*fskin[16]; 
  temp_edge[17] = 0.8660254037844387*fskin[13]-1.5*fskin[17]; 
  temp_edge[18] = (-7.5*fskin[18])+5.809475019311125*fskin[10]-3.354101966249685*fskin[5]; 
  temp_edge[19] = 0.8660254037844387*fskin[15]-1.5*fskin[19]; 

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

  temp_diff[0] = (-0.6708203932499369*fskin[8])+0.6708203932499369*fedge[8]+1.190784930203603*fskin[2]+1.190784930203603*fedge[2]-0.9375*fskin[0]+0.9375*fedge[0]; 
  temp_diff[1] = (-0.6708203932499369*fskin[12])+0.6708203932499369*fedge[12]+1.190784930203603*fskin[4]+1.190784930203603*fedge[4]-0.9375*fskin[1]+0.9375*fedge[1]; 
  temp_diff[2] = 1.585502557353661*fskin[8]-0.7382874503707888*fedge[8]-2.671875*fskin[2]-1.453125*fedge[2]+2.056810333988042*fskin[0]-1.190784930203603*fedge[0]; 
  temp_diff[3] = (-0.6708203932499369*fskin[14])+0.6708203932499369*fedge[14]+1.190784930203603*fskin[6]+1.190784930203603*fedge[6]-0.9375*fskin[3]+0.9375*fedge[3]; 
  temp_diff[4] = 1.585502557353661*fskin[12]-0.7382874503707888*fedge[12]-2.671875*fskin[4]-1.453125*fedge[4]+2.056810333988042*fskin[1]-1.190784930203603*fedge[1]; 
  temp_diff[5] = (-0.6708203932499369*fskin[18])+0.6708203932499369*fedge[18]+1.190784930203603*fskin[10]+1.190784930203603*fedge[10]-0.9375*fskin[5]+0.9375*fedge[5]; 
  temp_diff[6] = 1.585502557353661*fskin[14]-0.7382874503707888*fedge[14]-2.671875*fskin[6]-1.453125*fedge[6]+2.056810333988042*fskin[3]-1.190784930203603*fedge[3]; 
  temp_diff[7] = 1.190784930203603*fskin[11]+1.190784930203603*fedge[11]-0.9375*fskin[7]+0.9375*fedge[7]; 
  temp_diff[8] = (-3.140625*fskin[8])-0.140625*fedge[8]+5.022775277112744*fskin[2]+0.3025768239224545*fedge[2]-3.773364712030896*fskin[0]+0.4192627457812106*fedge[0]; 
  temp_diff[9] = 1.190784930203603*fskin[16]+1.190784930203603*fedge[16]-0.9375*fskin[9]+0.9375*fedge[9]; 
  temp_diff[10] = 1.585502557353661*fskin[18]-0.7382874503707888*fedge[18]-2.671875*fskin[10]-1.453125*fedge[10]+2.056810333988042*fskin[5]-1.190784930203603*fedge[5]; 
  temp_diff[11] = (-2.671875*fskin[11])-1.453125*fedge[11]+2.056810333988042*fskin[7]-1.190784930203603*fedge[7]; 
  temp_diff[12] = (-3.140625*fskin[12])-0.140625*fedge[12]+5.022775277112744*fskin[4]+0.3025768239224544*fedge[4]-3.773364712030894*fskin[1]+0.4192627457812105*fedge[1]; 
  temp_diff[13] = 1.190784930203603*fskin[17]+1.190784930203603*fedge[17]-0.9375*fskin[13]+0.9375*fedge[13]; 
  temp_diff[14] = (-3.140625*fskin[14])-0.140625*fedge[14]+5.022775277112744*fskin[6]+0.3025768239224544*fedge[6]-3.773364712030894*fskin[3]+0.4192627457812105*fedge[3]; 
  temp_diff[15] = 1.190784930203603*fskin[19]+1.190784930203603*fedge[19]-0.9375*fskin[15]+0.9375*fedge[15]; 
  temp_diff[16] = (-2.671875*fskin[16])-1.453125*fedge[16]+2.056810333988042*fskin[9]-1.190784930203603*fedge[9]; 
  temp_diff[17] = (-2.671875*fskin[17])-1.453125*fedge[17]+2.056810333988042*fskin[13]-1.190784930203603*fedge[13]; 
  temp_diff[18] = (-3.140625*fskin[18])-0.140625*fedge[18]+5.022775277112744*fskin[10]+0.3025768239224545*fedge[10]-3.773364712030896*fskin[5]+0.4192627457812106*fedge[5]; 
  temp_diff[19] = (-2.671875*fskin[19])-1.453125*fedge[19]+2.056810333988042*fskin[15]-1.190784930203603*fedge[15]; 

  temp_edge[2] = (-1.936491673103709*fskin[8])-1.5*fskin[2]-0.8660254037844386*fskin[0]; 
  temp_edge[4] = (-1.936491673103709*fskin[12])-1.5*fskin[4]-0.8660254037844386*fskin[1]; 
  temp_edge[6] = (-1.936491673103709*fskin[14])-1.5*fskin[6]-0.8660254037844386*fskin[3]; 
  temp_edge[8] = (-7.5*fskin[8])-5.809475019311125*fskin[2]-3.354101966249685*fskin[0]; 
  temp_edge[10] = (-1.936491673103709*fskin[18])-1.5*fskin[10]-0.8660254037844386*fskin[5]; 
  temp_edge[11] = (-1.5*fskin[11])-0.8660254037844387*fskin[7]; 
  temp_edge[12] = (-7.5*fskin[12])-5.809475019311126*fskin[4]-3.354101966249684*fskin[1]; 
  temp_edge[14] = (-7.5*fskin[14])-5.809475019311126*fskin[6]-3.354101966249684*fskin[3]; 
  temp_edge[16] = (-1.5*fskin[16])-0.8660254037844387*fskin[9]; 
  temp_edge[17] = (-1.5*fskin[17])-0.8660254037844387*fskin[13]; 
  temp_edge[18] = (-7.5*fskin[18])-5.809475019311125*fskin[10]-3.354101966249685*fskin[5]; 
  temp_edge[19] = (-1.5*fskin[19])-0.8660254037844387*fskin[15]; 

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
