#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH void lbo_gyrokinetic_diff_boundary_surfmu_1x2v_ser_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
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
  double rdvSq4 = 4.0/(dxv[2]*dxv[2]); 

  double facDiff[3]; 
  // Expand diffusion coefficient in conf basis.
  facDiff[0] = 1.414213562373095*(bmag_inv[2]*nuVtSqSum[2]+bmag_inv[1]*nuVtSqSum[1]+bmag_inv[0]*nuVtSqSum[0])*m_; 
  facDiff[1] = (1.264911064067352*(bmag_inv[1]*nuVtSqSum[2]+nuVtSqSum[1]*bmag_inv[2])+1.414213562373095*(bmag_inv[0]*nuVtSqSum[1]+nuVtSqSum[0]*bmag_inv[1]))*m_; 
  facDiff[2] = (0.9035079029052515*bmag_inv[2]*nuVtSqSum[2]+1.414213562373095*(bmag_inv[0]*nuVtSqSum[2]+nuVtSqSum[0]*bmag_inv[2])+1.264911064067352*bmag_inv[1]*nuVtSqSum[1])*m_; 

  double vol_incr[20] = {0.0};
  vol_incr[3] = 0.6123724356957944*dxv[2]*facDiff[2]*fskin[7]*rdvSq4+0.6123724356957944*facDiff[1]*fskin[1]*dxv[2]*rdvSq4+0.6123724356957944*facDiff[0]*fskin[0]*dxv[2]*rdvSq4; 
  vol_incr[5] = 0.5477225575051661*facDiff[1]*dxv[2]*fskin[7]*rdvSq4+0.5477225575051661*fskin[1]*dxv[2]*facDiff[2]*rdvSq4+0.6123724356957944*facDiff[0]*fskin[1]*dxv[2]*rdvSq4+0.6123724356957944*fskin[0]*facDiff[1]*dxv[2]*rdvSq4; 
  vol_incr[6] = 0.6123724356957944*dxv[2]*facDiff[2]*fskin[11]*rdvSq4+0.6123724356957944*facDiff[1]*dxv[2]*fskin[4]*rdvSq4+0.6123724356957944*facDiff[0]*dxv[2]*fskin[2]*rdvSq4; 
  vol_incr[9] = 2.738612787525831*dxv[2]*facDiff[2]*fskin[13]*rdvSq4+4.743416490252569*facDiff[2]*w[2]*fskin[7]*rdvSq4+2.738612787525831*facDiff[1]*dxv[2]*fskin[5]*rdvSq4+2.738612787525831*facDiff[0]*dxv[2]*fskin[3]*rdvSq4+4.743416490252569*facDiff[1]*fskin[1]*w[2]*rdvSq4+4.743416490252569*facDiff[0]*fskin[0]*w[2]*rdvSq4; 
  vol_incr[10] = 0.5477225575051661*facDiff[1]*dxv[2]*fskin[11]*rdvSq4+0.5477225575051661*dxv[2]*facDiff[2]*fskin[4]*rdvSq4+0.6123724356957944*facDiff[0]*dxv[2]*fskin[4]*rdvSq4+0.6123724356957944*facDiff[1]*dxv[2]*fskin[2]*rdvSq4; 
  vol_incr[13] = 0.3912303982179757*dxv[2]*facDiff[2]*fskin[7]*rdvSq4+0.6123724356957944*facDiff[0]*dxv[2]*fskin[7]*rdvSq4+0.6123724356957944*fskin[0]*dxv[2]*facDiff[2]*rdvSq4+0.5477225575051661*facDiff[1]*fskin[1]*dxv[2]*rdvSq4; 
  vol_incr[14] = 0.6123724356957944*facDiff[1]*dxv[2]*fskin[12]*rdvSq4+0.6123724356957944*facDiff[0]*dxv[2]*fskin[8]*rdvSq4; 
  vol_incr[15] = 2.449489742783178*facDiff[1]*dxv[2]*fskin[13]*rdvSq4+4.242640687119286*facDiff[1]*w[2]*fskin[7]*rdvSq4+2.449489742783178*dxv[2]*facDiff[2]*fskin[5]*rdvSq4+2.738612787525831*facDiff[0]*dxv[2]*fskin[5]*rdvSq4+2.738612787525831*facDiff[1]*dxv[2]*fskin[3]*rdvSq4+4.242640687119286*fskin[1]*facDiff[2]*w[2]*rdvSq4+4.743416490252569*facDiff[0]*fskin[1]*w[2]*rdvSq4+4.743416490252569*fskin[0]*facDiff[1]*w[2]*rdvSq4; 
  vol_incr[16] = 2.738612787525831*dxv[2]*facDiff[2]*fskin[17]*rdvSq4+4.743416490252569*facDiff[2]*w[2]*fskin[11]*rdvSq4+2.738612787525831*facDiff[1]*dxv[2]*fskin[10]*rdvSq4+2.738612787525831*facDiff[0]*dxv[2]*fskin[6]*rdvSq4+4.743416490252569*facDiff[1]*w[2]*fskin[4]*rdvSq4+4.743416490252569*facDiff[0]*fskin[2]*w[2]*rdvSq4; 
  vol_incr[17] = 0.3912303982179757*dxv[2]*facDiff[2]*fskin[11]*rdvSq4+0.6123724356957944*facDiff[0]*dxv[2]*fskin[11]*rdvSq4+0.5477225575051661*facDiff[1]*dxv[2]*fskin[4]*rdvSq4+0.6123724356957944*dxv[2]*facDiff[2]*fskin[2]*rdvSq4; 
  vol_incr[18] = 0.5477225575051661*dxv[2]*facDiff[2]*fskin[12]*rdvSq4+0.6123724356957944*facDiff[0]*dxv[2]*fskin[12]*rdvSq4+0.6123724356957944*facDiff[1]*dxv[2]*fskin[8]*rdvSq4; 
  vol_incr[19] = 2.449489742783178*facDiff[1]*dxv[2]*fskin[17]*rdvSq4+4.242640687119286*facDiff[1]*w[2]*fskin[11]*rdvSq4+2.449489742783178*dxv[2]*facDiff[2]*fskin[10]*rdvSq4+2.738612787525831*facDiff[0]*dxv[2]*fskin[10]*rdvSq4+2.738612787525831*facDiff[1]*dxv[2]*fskin[6]*rdvSq4+4.242640687119286*facDiff[2]*w[2]*fskin[4]*rdvSq4+4.743416490252569*facDiff[0]*w[2]*fskin[4]*rdvSq4+4.743416490252569*facDiff[1]*fskin[2]*w[2]*rdvSq4; 

  double temp_diff[20] = {0.0}; 
  double temp_edge[20] = {0.0}; 
  double diff_incr[20] = {0.0}; 
  double edge_incr[20] = {0.0}; 

  if (edge == -1) { 

  temp_diff[0] = (-0.6708203932499369*w[2]*fskin[9])-0.3354101966249685*dxv[2]*fskin[9]+0.6708203932499369*w[2]*fedge[9]+0.3354101966249685*dxv[2]*fedge[9]-1.190784930203603*w[2]*fskin[3]-0.5953924651018015*dxv[2]*fskin[3]-1.190784930203603*w[2]*fedge[3]-0.5953924651018015*dxv[2]*fedge[3]-0.9375*fskin[0]*w[2]+0.9375*fedge[0]*w[2]-0.46875*fskin[0]*dxv[2]+0.46875*fedge[0]*dxv[2]; 
  temp_diff[1] = (-0.6708203932499369*w[2]*fskin[15])-0.3354101966249685*dxv[2]*fskin[15]+0.6708203932499369*w[2]*fedge[15]+0.3354101966249685*dxv[2]*fedge[15]-1.190784930203603*w[2]*fskin[5]-0.5953924651018015*dxv[2]*fskin[5]-1.190784930203603*w[2]*fedge[5]-0.5953924651018015*dxv[2]*fedge[5]-0.9375*fskin[1]*w[2]+0.9375*fedge[1]*w[2]-0.46875*fskin[1]*dxv[2]+0.46875*fedge[1]*dxv[2]; 
  temp_diff[2] = (-0.6708203932499369*w[2]*fskin[16])-0.3354101966249685*dxv[2]*fskin[16]+0.6708203932499369*w[2]*fedge[16]+0.3354101966249685*dxv[2]*fedge[16]-1.190784930203603*w[2]*fskin[6]-0.5953924651018015*dxv[2]*fskin[6]-1.190784930203603*w[2]*fedge[6]-0.5953924651018015*dxv[2]*fedge[6]-0.9375*fskin[2]*w[2]+0.9375*fedge[2]*w[2]-0.46875*dxv[2]*fskin[2]+0.46875*dxv[2]*fedge[2]; 
  temp_diff[3] = (-1.585502557353661*w[2]*fskin[9])-0.7927512786768306*dxv[2]*fskin[9]+0.7382874503707888*w[2]*fedge[9]+0.3691437251853944*dxv[2]*fedge[9]-2.671875*w[2]*fskin[3]-1.3359375*dxv[2]*fskin[3]-1.453125*w[2]*fedge[3]-0.7265625*dxv[2]*fedge[3]-2.056810333988042*fskin[0]*w[2]+1.190784930203603*fedge[0]*w[2]-1.028405166994021*fskin[0]*dxv[2]+0.5953924651018015*fedge[0]*dxv[2]; 
  temp_diff[4] = (-0.6708203932499369*w[2]*fskin[19])-0.3354101966249685*dxv[2]*fskin[19]+0.6708203932499369*w[2]*fedge[19]+0.3354101966249685*dxv[2]*fedge[19]-1.190784930203603*w[2]*fskin[10]-0.5953924651018015*dxv[2]*fskin[10]-1.190784930203603*w[2]*fedge[10]-0.5953924651018015*dxv[2]*fedge[10]-0.9375*w[2]*fskin[4]-0.46875*dxv[2]*fskin[4]+0.9375*w[2]*fedge[4]+0.46875*dxv[2]*fedge[4]; 
  temp_diff[5] = (-1.585502557353661*w[2]*fskin[15])-0.7927512786768306*dxv[2]*fskin[15]+0.7382874503707888*w[2]*fedge[15]+0.3691437251853944*dxv[2]*fedge[15]-2.671875*w[2]*fskin[5]-1.3359375*dxv[2]*fskin[5]-1.453125*w[2]*fedge[5]-0.7265625*dxv[2]*fedge[5]-2.056810333988042*fskin[1]*w[2]+1.190784930203603*fedge[1]*w[2]-1.028405166994021*fskin[1]*dxv[2]+0.5953924651018015*fedge[1]*dxv[2]; 
  temp_diff[6] = (-1.585502557353661*w[2]*fskin[16])-0.7927512786768306*dxv[2]*fskin[16]+0.7382874503707888*w[2]*fedge[16]+0.3691437251853944*dxv[2]*fedge[16]-2.671875*w[2]*fskin[6]-1.3359375*dxv[2]*fskin[6]-1.453125*w[2]*fedge[6]-0.7265625*dxv[2]*fedge[6]-2.056810333988042*fskin[2]*w[2]+1.190784930203603*fedge[2]*w[2]-1.028405166994021*dxv[2]*fskin[2]+0.5953924651018015*dxv[2]*fedge[2]; 
  temp_diff[7] = (-1.190784930203603*w[2]*fskin[13])-0.5953924651018015*dxv[2]*fskin[13]-1.190784930203603*w[2]*fedge[13]-0.5953924651018015*dxv[2]*fedge[13]-0.9375*w[2]*fskin[7]-0.46875*dxv[2]*fskin[7]+0.9375*w[2]*fedge[7]+0.46875*dxv[2]*fedge[7]; 
  temp_diff[8] = (-1.190784930203603*w[2]*fskin[14])-0.5953924651018015*dxv[2]*fskin[14]-1.190784930203603*w[2]*fedge[14]-0.5953924651018015*dxv[2]*fedge[14]-0.9375*w[2]*fskin[8]-0.46875*dxv[2]*fskin[8]+0.9375*w[2]*fedge[8]+0.46875*dxv[2]*fedge[8]; 
  temp_diff[9] = (-3.140625*w[2]*fskin[9])-1.5703125*dxv[2]*fskin[9]-0.140625*w[2]*fedge[9]-0.0703125*dxv[2]*fedge[9]-5.022775277112744*w[2]*fskin[3]-2.511387638556372*dxv[2]*fskin[3]-0.3025768239224545*w[2]*fedge[3]-0.1512884119612272*dxv[2]*fedge[3]-3.773364712030896*fskin[0]*w[2]+0.4192627457812106*fedge[0]*w[2]-1.886682356015448*fskin[0]*dxv[2]+0.2096313728906053*fedge[0]*dxv[2]; 
  temp_diff[10] = (-1.585502557353661*w[2]*fskin[19])-0.7927512786768306*dxv[2]*fskin[19]+0.7382874503707888*w[2]*fedge[19]+0.3691437251853944*dxv[2]*fedge[19]-2.671875*w[2]*fskin[10]-1.3359375*dxv[2]*fskin[10]-1.453125*w[2]*fedge[10]-0.7265625*dxv[2]*fedge[10]-2.056810333988042*w[2]*fskin[4]-1.028405166994021*dxv[2]*fskin[4]+1.190784930203603*w[2]*fedge[4]+0.5953924651018015*dxv[2]*fedge[4]; 
  temp_diff[11] = (-1.190784930203603*w[2]*fskin[17])-0.5953924651018015*dxv[2]*fskin[17]-1.190784930203603*w[2]*fedge[17]-0.5953924651018015*dxv[2]*fedge[17]-0.9375*w[2]*fskin[11]-0.46875*dxv[2]*fskin[11]+0.9375*w[2]*fedge[11]+0.46875*dxv[2]*fedge[11]; 
  temp_diff[12] = (-1.190784930203603*w[2]*fskin[18])-0.5953924651018015*dxv[2]*fskin[18]-1.190784930203603*w[2]*fedge[18]-0.5953924651018015*dxv[2]*fedge[18]-0.9375*w[2]*fskin[12]-0.46875*dxv[2]*fskin[12]+0.9375*w[2]*fedge[12]+0.46875*dxv[2]*fedge[12]; 
  temp_diff[13] = (-2.671875*w[2]*fskin[13])-1.3359375*dxv[2]*fskin[13]-1.453125*w[2]*fedge[13]-0.7265625*dxv[2]*fedge[13]-2.056810333988042*w[2]*fskin[7]-1.028405166994021*dxv[2]*fskin[7]+1.190784930203603*w[2]*fedge[7]+0.5953924651018015*dxv[2]*fedge[7]; 
  temp_diff[14] = (-2.671875*w[2]*fskin[14])-1.3359375*dxv[2]*fskin[14]-1.453125*w[2]*fedge[14]-0.7265625*dxv[2]*fedge[14]-2.056810333988042*w[2]*fskin[8]-1.028405166994021*dxv[2]*fskin[8]+1.190784930203603*w[2]*fedge[8]+0.5953924651018015*dxv[2]*fedge[8]; 
  temp_diff[15] = (-3.140625*w[2]*fskin[15])-1.5703125*dxv[2]*fskin[15]-0.140625*w[2]*fedge[15]-0.0703125*dxv[2]*fedge[15]-5.022775277112744*w[2]*fskin[5]-2.511387638556372*dxv[2]*fskin[5]-0.3025768239224544*w[2]*fedge[5]-0.1512884119612272*dxv[2]*fedge[5]-3.773364712030894*fskin[1]*w[2]+0.4192627457812105*fedge[1]*w[2]-1.886682356015447*fskin[1]*dxv[2]+0.2096313728906053*fedge[1]*dxv[2]; 
  temp_diff[16] = (-3.140625*w[2]*fskin[16])-1.5703125*dxv[2]*fskin[16]-0.140625*w[2]*fedge[16]-0.0703125*dxv[2]*fedge[16]-5.022775277112744*w[2]*fskin[6]-2.511387638556372*dxv[2]*fskin[6]-0.3025768239224544*w[2]*fedge[6]-0.1512884119612272*dxv[2]*fedge[6]-3.773364712030894*fskin[2]*w[2]+0.4192627457812105*fedge[2]*w[2]-1.886682356015447*dxv[2]*fskin[2]+0.2096313728906053*dxv[2]*fedge[2]; 
  temp_diff[17] = (-2.671875*w[2]*fskin[17])-1.3359375*dxv[2]*fskin[17]-1.453125*w[2]*fedge[17]-0.7265625*dxv[2]*fedge[17]-2.056810333988042*w[2]*fskin[11]-1.028405166994021*dxv[2]*fskin[11]+1.190784930203603*w[2]*fedge[11]+0.5953924651018015*dxv[2]*fedge[11]; 
  temp_diff[18] = (-2.671875*w[2]*fskin[18])-1.3359375*dxv[2]*fskin[18]-1.453125*w[2]*fedge[18]-0.7265625*dxv[2]*fedge[18]-2.056810333988042*w[2]*fskin[12]-1.028405166994021*dxv[2]*fskin[12]+1.190784930203603*w[2]*fedge[12]+0.5953924651018015*dxv[2]*fedge[12]; 
  temp_diff[19] = (-3.140625*w[2]*fskin[19])-1.5703125*dxv[2]*fskin[19]-0.140625*w[2]*fedge[19]-0.0703125*dxv[2]*fedge[19]-5.022775277112744*w[2]*fskin[10]-2.511387638556372*dxv[2]*fskin[10]-0.3025768239224545*w[2]*fedge[10]-0.1512884119612272*dxv[2]*fedge[10]-3.773364712030896*w[2]*fskin[4]-1.886682356015448*dxv[2]*fskin[4]+0.4192627457812106*w[2]*fedge[4]+0.2096313728906053*dxv[2]*fedge[4]; 

  temp_edge[3] = 1.936491673103709*w[2]*fskin[9]-0.9682458365518543*dxv[2]*fskin[9]-1.5*w[2]*fskin[3]+0.75*dxv[2]*fskin[3]+0.8660254037844386*fskin[0]*w[2]-0.4330127018922193*fskin[0]*dxv[2]; 
  temp_edge[5] = 1.936491673103709*w[2]*fskin[15]-0.9682458365518543*dxv[2]*fskin[15]-1.5*w[2]*fskin[5]+0.75*dxv[2]*fskin[5]+0.8660254037844386*fskin[1]*w[2]-0.4330127018922193*fskin[1]*dxv[2]; 
  temp_edge[6] = 1.936491673103709*w[2]*fskin[16]-0.9682458365518543*dxv[2]*fskin[16]-1.5*w[2]*fskin[6]+0.75*dxv[2]*fskin[6]+0.8660254037844386*fskin[2]*w[2]-0.4330127018922193*dxv[2]*fskin[2]; 
  temp_edge[9] = (-7.5*w[2]*fskin[9])+3.75*dxv[2]*fskin[9]+5.809475019311125*w[2]*fskin[3]-2.904737509655563*dxv[2]*fskin[3]-3.354101966249685*fskin[0]*w[2]+1.677050983124842*fskin[0]*dxv[2]; 
  temp_edge[10] = 1.936491673103709*w[2]*fskin[19]-0.9682458365518543*dxv[2]*fskin[19]-1.5*w[2]*fskin[10]+0.75*dxv[2]*fskin[10]+0.8660254037844386*w[2]*fskin[4]-0.4330127018922193*dxv[2]*fskin[4]; 
  temp_edge[13] = (-1.5*w[2]*fskin[13])+0.75*dxv[2]*fskin[13]+0.8660254037844387*w[2]*fskin[7]-0.4330127018922194*dxv[2]*fskin[7]; 
  temp_edge[14] = (-1.5*w[2]*fskin[14])+0.75*dxv[2]*fskin[14]+0.8660254037844387*w[2]*fskin[8]-0.4330127018922194*dxv[2]*fskin[8]; 
  temp_edge[15] = (-7.5*w[2]*fskin[15])+3.75*dxv[2]*fskin[15]+5.809475019311126*w[2]*fskin[5]-2.904737509655563*dxv[2]*fskin[5]-3.354101966249684*fskin[1]*w[2]+1.677050983124842*fskin[1]*dxv[2]; 
  temp_edge[16] = (-7.5*w[2]*fskin[16])+3.75*dxv[2]*fskin[16]+5.809475019311126*w[2]*fskin[6]-2.904737509655563*dxv[2]*fskin[6]-3.354101966249684*fskin[2]*w[2]+1.677050983124842*dxv[2]*fskin[2]; 
  temp_edge[17] = (-1.5*w[2]*fskin[17])+0.75*dxv[2]*fskin[17]+0.8660254037844387*w[2]*fskin[11]-0.4330127018922194*dxv[2]*fskin[11]; 
  temp_edge[18] = (-1.5*w[2]*fskin[18])+0.75*dxv[2]*fskin[18]+0.8660254037844387*w[2]*fskin[12]-0.4330127018922194*dxv[2]*fskin[12]; 
  temp_edge[19] = (-7.5*w[2]*fskin[19])+3.75*dxv[2]*fskin[19]+5.809475019311125*w[2]*fskin[10]-2.904737509655563*dxv[2]*fskin[10]-3.354101966249685*w[2]*fskin[4]+1.677050983124842*dxv[2]*fskin[4]; 

  diff_incr[0] = 0.7071067811865475*facDiff[2]*temp_diff[7]+0.7071067811865475*facDiff[1]*temp_diff[1]+0.7071067811865475*facDiff[0]*temp_diff[0]; 
  diff_incr[1] = 0.6324555320336759*facDiff[1]*temp_diff[7]+0.6324555320336759*temp_diff[1]*facDiff[2]+0.7071067811865475*facDiff[0]*temp_diff[1]+0.7071067811865475*temp_diff[0]*facDiff[1]; 
  diff_incr[2] = 0.7071067811865475*facDiff[2]*temp_diff[11]+0.7071067811865475*facDiff[1]*temp_diff[4]+0.7071067811865475*facDiff[0]*temp_diff[2]; 
  diff_incr[3] = 0.7071067811865475*facDiff[2]*temp_diff[13]+0.7071067811865475*facDiff[1]*temp_diff[5]+0.7071067811865475*facDiff[0]*temp_diff[3]; 
  diff_incr[4] = 0.632455532033676*facDiff[1]*temp_diff[11]+0.6324555320336759*facDiff[2]*temp_diff[4]+0.7071067811865475*facDiff[0]*temp_diff[4]+0.7071067811865475*facDiff[1]*temp_diff[2]; 
  diff_incr[5] = 0.632455532033676*facDiff[1]*temp_diff[13]+0.6324555320336759*facDiff[2]*temp_diff[5]+0.7071067811865475*facDiff[0]*temp_diff[5]+0.7071067811865475*facDiff[1]*temp_diff[3]; 
  diff_incr[6] = 0.7071067811865475*facDiff[2]*temp_diff[17]+0.7071067811865475*facDiff[1]*temp_diff[10]+0.7071067811865475*facDiff[0]*temp_diff[6]; 
  diff_incr[7] = 0.4517539514526256*facDiff[2]*temp_diff[7]+0.7071067811865475*facDiff[0]*temp_diff[7]+0.7071067811865475*temp_diff[0]*facDiff[2]+0.6324555320336759*facDiff[1]*temp_diff[1]; 
  diff_incr[8] = 0.7071067811865475*facDiff[1]*temp_diff[12]+0.7071067811865475*facDiff[0]*temp_diff[8]; 
  diff_incr[9] = 0.7071067811865475*facDiff[1]*temp_diff[15]+0.7071067811865475*facDiff[0]*temp_diff[9]; 
  diff_incr[10] = 0.6324555320336759*facDiff[1]*temp_diff[17]+0.6324555320336759*facDiff[2]*temp_diff[10]+0.7071067811865475*facDiff[0]*temp_diff[10]+0.7071067811865475*facDiff[1]*temp_diff[6]; 
  diff_incr[11] = 0.4517539514526256*facDiff[2]*temp_diff[11]+0.7071067811865475*facDiff[0]*temp_diff[11]+0.632455532033676*facDiff[1]*temp_diff[4]+0.7071067811865475*facDiff[2]*temp_diff[2]; 
  diff_incr[12] = 0.6324555320336759*facDiff[2]*temp_diff[12]+0.7071067811865475*facDiff[0]*temp_diff[12]+0.7071067811865475*facDiff[1]*temp_diff[8]; 
  diff_incr[13] = 0.4517539514526256*facDiff[2]*temp_diff[13]+0.7071067811865475*facDiff[0]*temp_diff[13]+0.632455532033676*facDiff[1]*temp_diff[5]+0.7071067811865475*facDiff[2]*temp_diff[3]; 
  diff_incr[14] = 0.7071067811865475*facDiff[1]*temp_diff[18]+0.7071067811865475*facDiff[0]*temp_diff[14]; 
  diff_incr[15] = 0.6324555320336759*facDiff[2]*temp_diff[15]+0.7071067811865475*facDiff[0]*temp_diff[15]+0.7071067811865475*facDiff[1]*temp_diff[9]; 
  diff_incr[16] = 0.7071067811865475*facDiff[1]*temp_diff[19]+0.7071067811865475*facDiff[0]*temp_diff[16]; 
  diff_incr[17] = 0.4517539514526256*facDiff[2]*temp_diff[17]+0.7071067811865475*facDiff[0]*temp_diff[17]+0.6324555320336759*facDiff[1]*temp_diff[10]+0.7071067811865475*facDiff[2]*temp_diff[6]; 
  diff_incr[18] = 0.6324555320336759*facDiff[2]*temp_diff[18]+0.7071067811865475*facDiff[0]*temp_diff[18]+0.7071067811865475*facDiff[1]*temp_diff[14]; 
  diff_incr[19] = 0.6324555320336759*facDiff[2]*temp_diff[19]+0.7071067811865475*facDiff[0]*temp_diff[19]+0.7071067811865475*facDiff[1]*temp_diff[16]; 

  edge_incr[0] = 0.7071067811865475*facDiff[2]*temp_edge[7]+0.7071067811865475*facDiff[1]*temp_edge[1]+0.7071067811865475*facDiff[0]*temp_edge[0]; 
  edge_incr[1] = 0.6324555320336759*facDiff[1]*temp_edge[7]+0.6324555320336759*temp_edge[1]*facDiff[2]+0.7071067811865475*facDiff[0]*temp_edge[1]+0.7071067811865475*temp_edge[0]*facDiff[1]; 
  edge_incr[2] = 0.7071067811865475*facDiff[2]*temp_edge[11]+0.7071067811865475*facDiff[1]*temp_edge[4]+0.7071067811865475*facDiff[0]*temp_edge[2]; 
  edge_incr[3] = 0.7071067811865475*facDiff[2]*temp_edge[13]+0.7071067811865475*facDiff[1]*temp_edge[5]+0.7071067811865475*facDiff[0]*temp_edge[3]; 
  edge_incr[4] = 0.632455532033676*facDiff[1]*temp_edge[11]+0.6324555320336759*facDiff[2]*temp_edge[4]+0.7071067811865475*facDiff[0]*temp_edge[4]+0.7071067811865475*facDiff[1]*temp_edge[2]; 
  edge_incr[5] = 0.632455532033676*facDiff[1]*temp_edge[13]+0.6324555320336759*facDiff[2]*temp_edge[5]+0.7071067811865475*facDiff[0]*temp_edge[5]+0.7071067811865475*facDiff[1]*temp_edge[3]; 
  edge_incr[6] = 0.7071067811865475*facDiff[2]*temp_edge[17]+0.7071067811865475*facDiff[1]*temp_edge[10]+0.7071067811865475*facDiff[0]*temp_edge[6]; 
  edge_incr[7] = 0.4517539514526256*facDiff[2]*temp_edge[7]+0.7071067811865475*facDiff[0]*temp_edge[7]+0.7071067811865475*temp_edge[0]*facDiff[2]+0.6324555320336759*facDiff[1]*temp_edge[1]; 
  edge_incr[8] = 0.7071067811865475*facDiff[1]*temp_edge[12]+0.7071067811865475*facDiff[0]*temp_edge[8]; 
  edge_incr[9] = 0.7071067811865475*facDiff[1]*temp_edge[15]+0.7071067811865475*facDiff[0]*temp_edge[9]; 
  edge_incr[10] = 0.6324555320336759*facDiff[1]*temp_edge[17]+0.6324555320336759*facDiff[2]*temp_edge[10]+0.7071067811865475*facDiff[0]*temp_edge[10]+0.7071067811865475*facDiff[1]*temp_edge[6]; 
  edge_incr[11] = 0.4517539514526256*facDiff[2]*temp_edge[11]+0.7071067811865475*facDiff[0]*temp_edge[11]+0.632455532033676*facDiff[1]*temp_edge[4]+0.7071067811865475*facDiff[2]*temp_edge[2]; 
  edge_incr[12] = 0.6324555320336759*facDiff[2]*temp_edge[12]+0.7071067811865475*facDiff[0]*temp_edge[12]+0.7071067811865475*facDiff[1]*temp_edge[8]; 
  edge_incr[13] = 0.4517539514526256*facDiff[2]*temp_edge[13]+0.7071067811865475*facDiff[0]*temp_edge[13]+0.632455532033676*facDiff[1]*temp_edge[5]+0.7071067811865475*facDiff[2]*temp_edge[3]; 
  edge_incr[14] = 0.7071067811865475*facDiff[1]*temp_edge[18]+0.7071067811865475*facDiff[0]*temp_edge[14]; 
  edge_incr[15] = 0.6324555320336759*facDiff[2]*temp_edge[15]+0.7071067811865475*facDiff[0]*temp_edge[15]+0.7071067811865475*facDiff[1]*temp_edge[9]; 
  edge_incr[16] = 0.7071067811865475*facDiff[1]*temp_edge[19]+0.7071067811865475*facDiff[0]*temp_edge[16]; 
  edge_incr[17] = 0.4517539514526256*facDiff[2]*temp_edge[17]+0.7071067811865475*facDiff[0]*temp_edge[17]+0.6324555320336759*facDiff[1]*temp_edge[10]+0.7071067811865475*facDiff[2]*temp_edge[6]; 
  edge_incr[18] = 0.6324555320336759*facDiff[2]*temp_edge[18]+0.7071067811865475*facDiff[0]*temp_edge[18]+0.7071067811865475*facDiff[1]*temp_edge[14]; 
  edge_incr[19] = 0.6324555320336759*facDiff[2]*temp_edge[19]+0.7071067811865475*facDiff[0]*temp_edge[19]+0.7071067811865475*facDiff[1]*temp_edge[16]; 


  } else { 

  temp_diff[0] = (-0.6708203932499369*w[2]*fskin[9])+0.3354101966249685*dxv[2]*fskin[9]+0.6708203932499369*w[2]*fedge[9]-0.3354101966249685*dxv[2]*fedge[9]+1.190784930203603*w[2]*fskin[3]-0.5953924651018015*dxv[2]*fskin[3]+1.190784930203603*w[2]*fedge[3]-0.5953924651018015*dxv[2]*fedge[3]-0.9375*fskin[0]*w[2]+0.9375*fedge[0]*w[2]+0.46875*fskin[0]*dxv[2]-0.46875*fedge[0]*dxv[2]; 
  temp_diff[1] = (-0.6708203932499369*w[2]*fskin[15])+0.3354101966249685*dxv[2]*fskin[15]+0.6708203932499369*w[2]*fedge[15]-0.3354101966249685*dxv[2]*fedge[15]+1.190784930203603*w[2]*fskin[5]-0.5953924651018015*dxv[2]*fskin[5]+1.190784930203603*w[2]*fedge[5]-0.5953924651018015*dxv[2]*fedge[5]-0.9375*fskin[1]*w[2]+0.9375*fedge[1]*w[2]+0.46875*fskin[1]*dxv[2]-0.46875*fedge[1]*dxv[2]; 
  temp_diff[2] = (-0.6708203932499369*w[2]*fskin[16])+0.3354101966249685*dxv[2]*fskin[16]+0.6708203932499369*w[2]*fedge[16]-0.3354101966249685*dxv[2]*fedge[16]+1.190784930203603*w[2]*fskin[6]-0.5953924651018015*dxv[2]*fskin[6]+1.190784930203603*w[2]*fedge[6]-0.5953924651018015*dxv[2]*fedge[6]-0.9375*fskin[2]*w[2]+0.9375*fedge[2]*w[2]+0.46875*dxv[2]*fskin[2]-0.46875*dxv[2]*fedge[2]; 
  temp_diff[3] = 1.585502557353661*w[2]*fskin[9]-0.7927512786768306*dxv[2]*fskin[9]-0.7382874503707888*w[2]*fedge[9]+0.3691437251853944*dxv[2]*fedge[9]-2.671875*w[2]*fskin[3]+1.3359375*dxv[2]*fskin[3]-1.453125*w[2]*fedge[3]+0.7265625*dxv[2]*fedge[3]+2.056810333988042*fskin[0]*w[2]-1.190784930203603*fedge[0]*w[2]-1.028405166994021*fskin[0]*dxv[2]+0.5953924651018015*fedge[0]*dxv[2]; 
  temp_diff[4] = (-0.6708203932499369*w[2]*fskin[19])+0.3354101966249685*dxv[2]*fskin[19]+0.6708203932499369*w[2]*fedge[19]-0.3354101966249685*dxv[2]*fedge[19]+1.190784930203603*w[2]*fskin[10]-0.5953924651018015*dxv[2]*fskin[10]+1.190784930203603*w[2]*fedge[10]-0.5953924651018015*dxv[2]*fedge[10]-0.9375*w[2]*fskin[4]+0.46875*dxv[2]*fskin[4]+0.9375*w[2]*fedge[4]-0.46875*dxv[2]*fedge[4]; 
  temp_diff[5] = 1.585502557353661*w[2]*fskin[15]-0.7927512786768306*dxv[2]*fskin[15]-0.7382874503707888*w[2]*fedge[15]+0.3691437251853944*dxv[2]*fedge[15]-2.671875*w[2]*fskin[5]+1.3359375*dxv[2]*fskin[5]-1.453125*w[2]*fedge[5]+0.7265625*dxv[2]*fedge[5]+2.056810333988042*fskin[1]*w[2]-1.190784930203603*fedge[1]*w[2]-1.028405166994021*fskin[1]*dxv[2]+0.5953924651018015*fedge[1]*dxv[2]; 
  temp_diff[6] = 1.585502557353661*w[2]*fskin[16]-0.7927512786768306*dxv[2]*fskin[16]-0.7382874503707888*w[2]*fedge[16]+0.3691437251853944*dxv[2]*fedge[16]-2.671875*w[2]*fskin[6]+1.3359375*dxv[2]*fskin[6]-1.453125*w[2]*fedge[6]+0.7265625*dxv[2]*fedge[6]+2.056810333988042*fskin[2]*w[2]-1.190784930203603*fedge[2]*w[2]-1.028405166994021*dxv[2]*fskin[2]+0.5953924651018015*dxv[2]*fedge[2]; 
  temp_diff[7] = 1.190784930203603*w[2]*fskin[13]-0.5953924651018015*dxv[2]*fskin[13]+1.190784930203603*w[2]*fedge[13]-0.5953924651018015*dxv[2]*fedge[13]-0.9375*w[2]*fskin[7]+0.46875*dxv[2]*fskin[7]+0.9375*w[2]*fedge[7]-0.46875*dxv[2]*fedge[7]; 
  temp_diff[8] = 1.190784930203603*w[2]*fskin[14]-0.5953924651018015*dxv[2]*fskin[14]+1.190784930203603*w[2]*fedge[14]-0.5953924651018015*dxv[2]*fedge[14]-0.9375*w[2]*fskin[8]+0.46875*dxv[2]*fskin[8]+0.9375*w[2]*fedge[8]-0.46875*dxv[2]*fedge[8]; 
  temp_diff[9] = (-3.140625*w[2]*fskin[9])+1.5703125*dxv[2]*fskin[9]-0.140625*w[2]*fedge[9]+0.0703125*dxv[2]*fedge[9]+5.022775277112744*w[2]*fskin[3]-2.511387638556372*dxv[2]*fskin[3]+0.3025768239224545*w[2]*fedge[3]-0.1512884119612272*dxv[2]*fedge[3]-3.773364712030896*fskin[0]*w[2]+0.4192627457812106*fedge[0]*w[2]+1.886682356015448*fskin[0]*dxv[2]-0.2096313728906053*fedge[0]*dxv[2]; 
  temp_diff[10] = 1.585502557353661*w[2]*fskin[19]-0.7927512786768306*dxv[2]*fskin[19]-0.7382874503707888*w[2]*fedge[19]+0.3691437251853944*dxv[2]*fedge[19]-2.671875*w[2]*fskin[10]+1.3359375*dxv[2]*fskin[10]-1.453125*w[2]*fedge[10]+0.7265625*dxv[2]*fedge[10]+2.056810333988042*w[2]*fskin[4]-1.028405166994021*dxv[2]*fskin[4]-1.190784930203603*w[2]*fedge[4]+0.5953924651018015*dxv[2]*fedge[4]; 
  temp_diff[11] = 1.190784930203603*w[2]*fskin[17]-0.5953924651018015*dxv[2]*fskin[17]+1.190784930203603*w[2]*fedge[17]-0.5953924651018015*dxv[2]*fedge[17]-0.9375*w[2]*fskin[11]+0.46875*dxv[2]*fskin[11]+0.9375*w[2]*fedge[11]-0.46875*dxv[2]*fedge[11]; 
  temp_diff[12] = 1.190784930203603*w[2]*fskin[18]-0.5953924651018015*dxv[2]*fskin[18]+1.190784930203603*w[2]*fedge[18]-0.5953924651018015*dxv[2]*fedge[18]-0.9375*w[2]*fskin[12]+0.46875*dxv[2]*fskin[12]+0.9375*w[2]*fedge[12]-0.46875*dxv[2]*fedge[12]; 
  temp_diff[13] = (-2.671875*w[2]*fskin[13])+1.3359375*dxv[2]*fskin[13]-1.453125*w[2]*fedge[13]+0.7265625*dxv[2]*fedge[13]+2.056810333988042*w[2]*fskin[7]-1.028405166994021*dxv[2]*fskin[7]-1.190784930203603*w[2]*fedge[7]+0.5953924651018015*dxv[2]*fedge[7]; 
  temp_diff[14] = (-2.671875*w[2]*fskin[14])+1.3359375*dxv[2]*fskin[14]-1.453125*w[2]*fedge[14]+0.7265625*dxv[2]*fedge[14]+2.056810333988042*w[2]*fskin[8]-1.028405166994021*dxv[2]*fskin[8]-1.190784930203603*w[2]*fedge[8]+0.5953924651018015*dxv[2]*fedge[8]; 
  temp_diff[15] = (-3.140625*w[2]*fskin[15])+1.5703125*dxv[2]*fskin[15]-0.140625*w[2]*fedge[15]+0.0703125*dxv[2]*fedge[15]+5.022775277112744*w[2]*fskin[5]-2.511387638556372*dxv[2]*fskin[5]+0.3025768239224544*w[2]*fedge[5]-0.1512884119612272*dxv[2]*fedge[5]-3.773364712030894*fskin[1]*w[2]+0.4192627457812105*fedge[1]*w[2]+1.886682356015447*fskin[1]*dxv[2]-0.2096313728906053*fedge[1]*dxv[2]; 
  temp_diff[16] = (-3.140625*w[2]*fskin[16])+1.5703125*dxv[2]*fskin[16]-0.140625*w[2]*fedge[16]+0.0703125*dxv[2]*fedge[16]+5.022775277112744*w[2]*fskin[6]-2.511387638556372*dxv[2]*fskin[6]+0.3025768239224544*w[2]*fedge[6]-0.1512884119612272*dxv[2]*fedge[6]-3.773364712030894*fskin[2]*w[2]+0.4192627457812105*fedge[2]*w[2]+1.886682356015447*dxv[2]*fskin[2]-0.2096313728906053*dxv[2]*fedge[2]; 
  temp_diff[17] = (-2.671875*w[2]*fskin[17])+1.3359375*dxv[2]*fskin[17]-1.453125*w[2]*fedge[17]+0.7265625*dxv[2]*fedge[17]+2.056810333988042*w[2]*fskin[11]-1.028405166994021*dxv[2]*fskin[11]-1.190784930203603*w[2]*fedge[11]+0.5953924651018015*dxv[2]*fedge[11]; 
  temp_diff[18] = (-2.671875*w[2]*fskin[18])+1.3359375*dxv[2]*fskin[18]-1.453125*w[2]*fedge[18]+0.7265625*dxv[2]*fedge[18]+2.056810333988042*w[2]*fskin[12]-1.028405166994021*dxv[2]*fskin[12]-1.190784930203603*w[2]*fedge[12]+0.5953924651018015*dxv[2]*fedge[12]; 
  temp_diff[19] = (-3.140625*w[2]*fskin[19])+1.5703125*dxv[2]*fskin[19]-0.140625*w[2]*fedge[19]+0.0703125*dxv[2]*fedge[19]+5.022775277112744*w[2]*fskin[10]-2.511387638556372*dxv[2]*fskin[10]+0.3025768239224545*w[2]*fedge[10]-0.1512884119612272*dxv[2]*fedge[10]-3.773364712030896*w[2]*fskin[4]+1.886682356015448*dxv[2]*fskin[4]+0.4192627457812106*w[2]*fedge[4]-0.2096313728906053*dxv[2]*fedge[4]; 

  temp_edge[3] = (-1.936491673103709*w[2]*fskin[9])-0.9682458365518543*dxv[2]*fskin[9]-1.5*w[2]*fskin[3]-0.75*dxv[2]*fskin[3]-0.8660254037844386*fskin[0]*w[2]-0.4330127018922193*fskin[0]*dxv[2]; 
  temp_edge[5] = (-1.936491673103709*w[2]*fskin[15])-0.9682458365518543*dxv[2]*fskin[15]-1.5*w[2]*fskin[5]-0.75*dxv[2]*fskin[5]-0.8660254037844386*fskin[1]*w[2]-0.4330127018922193*fskin[1]*dxv[2]; 
  temp_edge[6] = (-1.936491673103709*w[2]*fskin[16])-0.9682458365518543*dxv[2]*fskin[16]-1.5*w[2]*fskin[6]-0.75*dxv[2]*fskin[6]-0.8660254037844386*fskin[2]*w[2]-0.4330127018922193*dxv[2]*fskin[2]; 
  temp_edge[9] = (-7.5*w[2]*fskin[9])-3.75*dxv[2]*fskin[9]-5.809475019311125*w[2]*fskin[3]-2.904737509655563*dxv[2]*fskin[3]-3.354101966249685*fskin[0]*w[2]-1.677050983124842*fskin[0]*dxv[2]; 
  temp_edge[10] = (-1.936491673103709*w[2]*fskin[19])-0.9682458365518543*dxv[2]*fskin[19]-1.5*w[2]*fskin[10]-0.75*dxv[2]*fskin[10]-0.8660254037844386*w[2]*fskin[4]-0.4330127018922193*dxv[2]*fskin[4]; 
  temp_edge[13] = (-1.5*w[2]*fskin[13])-0.75*dxv[2]*fskin[13]-0.8660254037844387*w[2]*fskin[7]-0.4330127018922194*dxv[2]*fskin[7]; 
  temp_edge[14] = (-1.5*w[2]*fskin[14])-0.75*dxv[2]*fskin[14]-0.8660254037844387*w[2]*fskin[8]-0.4330127018922194*dxv[2]*fskin[8]; 
  temp_edge[15] = (-7.5*w[2]*fskin[15])-3.75*dxv[2]*fskin[15]-5.809475019311126*w[2]*fskin[5]-2.904737509655563*dxv[2]*fskin[5]-3.354101966249684*fskin[1]*w[2]-1.677050983124842*fskin[1]*dxv[2]; 
  temp_edge[16] = (-7.5*w[2]*fskin[16])-3.75*dxv[2]*fskin[16]-5.809475019311126*w[2]*fskin[6]-2.904737509655563*dxv[2]*fskin[6]-3.354101966249684*fskin[2]*w[2]-1.677050983124842*dxv[2]*fskin[2]; 
  temp_edge[17] = (-1.5*w[2]*fskin[17])-0.75*dxv[2]*fskin[17]-0.8660254037844387*w[2]*fskin[11]-0.4330127018922194*dxv[2]*fskin[11]; 
  temp_edge[18] = (-1.5*w[2]*fskin[18])-0.75*dxv[2]*fskin[18]-0.8660254037844387*w[2]*fskin[12]-0.4330127018922194*dxv[2]*fskin[12]; 
  temp_edge[19] = (-7.5*w[2]*fskin[19])-3.75*dxv[2]*fskin[19]-5.809475019311125*w[2]*fskin[10]-2.904737509655563*dxv[2]*fskin[10]-3.354101966249685*w[2]*fskin[4]-1.677050983124842*dxv[2]*fskin[4]; 

  diff_incr[0] = 0.7071067811865475*facDiff[2]*temp_diff[7]+0.7071067811865475*facDiff[1]*temp_diff[1]+0.7071067811865475*facDiff[0]*temp_diff[0]; 
  diff_incr[1] = 0.6324555320336759*facDiff[1]*temp_diff[7]+0.6324555320336759*temp_diff[1]*facDiff[2]+0.7071067811865475*facDiff[0]*temp_diff[1]+0.7071067811865475*temp_diff[0]*facDiff[1]; 
  diff_incr[2] = 0.7071067811865475*facDiff[2]*temp_diff[11]+0.7071067811865475*facDiff[1]*temp_diff[4]+0.7071067811865475*facDiff[0]*temp_diff[2]; 
  diff_incr[3] = 0.7071067811865475*facDiff[2]*temp_diff[13]+0.7071067811865475*facDiff[1]*temp_diff[5]+0.7071067811865475*facDiff[0]*temp_diff[3]; 
  diff_incr[4] = 0.632455532033676*facDiff[1]*temp_diff[11]+0.6324555320336759*facDiff[2]*temp_diff[4]+0.7071067811865475*facDiff[0]*temp_diff[4]+0.7071067811865475*facDiff[1]*temp_diff[2]; 
  diff_incr[5] = 0.632455532033676*facDiff[1]*temp_diff[13]+0.6324555320336759*facDiff[2]*temp_diff[5]+0.7071067811865475*facDiff[0]*temp_diff[5]+0.7071067811865475*facDiff[1]*temp_diff[3]; 
  diff_incr[6] = 0.7071067811865475*facDiff[2]*temp_diff[17]+0.7071067811865475*facDiff[1]*temp_diff[10]+0.7071067811865475*facDiff[0]*temp_diff[6]; 
  diff_incr[7] = 0.4517539514526256*facDiff[2]*temp_diff[7]+0.7071067811865475*facDiff[0]*temp_diff[7]+0.7071067811865475*temp_diff[0]*facDiff[2]+0.6324555320336759*facDiff[1]*temp_diff[1]; 
  diff_incr[8] = 0.7071067811865475*facDiff[1]*temp_diff[12]+0.7071067811865475*facDiff[0]*temp_diff[8]; 
  diff_incr[9] = 0.7071067811865475*facDiff[1]*temp_diff[15]+0.7071067811865475*facDiff[0]*temp_diff[9]; 
  diff_incr[10] = 0.6324555320336759*facDiff[1]*temp_diff[17]+0.6324555320336759*facDiff[2]*temp_diff[10]+0.7071067811865475*facDiff[0]*temp_diff[10]+0.7071067811865475*facDiff[1]*temp_diff[6]; 
  diff_incr[11] = 0.4517539514526256*facDiff[2]*temp_diff[11]+0.7071067811865475*facDiff[0]*temp_diff[11]+0.632455532033676*facDiff[1]*temp_diff[4]+0.7071067811865475*facDiff[2]*temp_diff[2]; 
  diff_incr[12] = 0.6324555320336759*facDiff[2]*temp_diff[12]+0.7071067811865475*facDiff[0]*temp_diff[12]+0.7071067811865475*facDiff[1]*temp_diff[8]; 
  diff_incr[13] = 0.4517539514526256*facDiff[2]*temp_diff[13]+0.7071067811865475*facDiff[0]*temp_diff[13]+0.632455532033676*facDiff[1]*temp_diff[5]+0.7071067811865475*facDiff[2]*temp_diff[3]; 
  diff_incr[14] = 0.7071067811865475*facDiff[1]*temp_diff[18]+0.7071067811865475*facDiff[0]*temp_diff[14]; 
  diff_incr[15] = 0.6324555320336759*facDiff[2]*temp_diff[15]+0.7071067811865475*facDiff[0]*temp_diff[15]+0.7071067811865475*facDiff[1]*temp_diff[9]; 
  diff_incr[16] = 0.7071067811865475*facDiff[1]*temp_diff[19]+0.7071067811865475*facDiff[0]*temp_diff[16]; 
  diff_incr[17] = 0.4517539514526256*facDiff[2]*temp_diff[17]+0.7071067811865475*facDiff[0]*temp_diff[17]+0.6324555320336759*facDiff[1]*temp_diff[10]+0.7071067811865475*facDiff[2]*temp_diff[6]; 
  diff_incr[18] = 0.6324555320336759*facDiff[2]*temp_diff[18]+0.7071067811865475*facDiff[0]*temp_diff[18]+0.7071067811865475*facDiff[1]*temp_diff[14]; 
  diff_incr[19] = 0.6324555320336759*facDiff[2]*temp_diff[19]+0.7071067811865475*facDiff[0]*temp_diff[19]+0.7071067811865475*facDiff[1]*temp_diff[16]; 

  edge_incr[0] = 0.7071067811865475*facDiff[2]*temp_edge[7]+0.7071067811865475*facDiff[1]*temp_edge[1]+0.7071067811865475*facDiff[0]*temp_edge[0]; 
  edge_incr[1] = 0.6324555320336759*facDiff[1]*temp_edge[7]+0.6324555320336759*temp_edge[1]*facDiff[2]+0.7071067811865475*facDiff[0]*temp_edge[1]+0.7071067811865475*temp_edge[0]*facDiff[1]; 
  edge_incr[2] = 0.7071067811865475*facDiff[2]*temp_edge[11]+0.7071067811865475*facDiff[1]*temp_edge[4]+0.7071067811865475*facDiff[0]*temp_edge[2]; 
  edge_incr[3] = 0.7071067811865475*facDiff[2]*temp_edge[13]+0.7071067811865475*facDiff[1]*temp_edge[5]+0.7071067811865475*facDiff[0]*temp_edge[3]; 
  edge_incr[4] = 0.632455532033676*facDiff[1]*temp_edge[11]+0.6324555320336759*facDiff[2]*temp_edge[4]+0.7071067811865475*facDiff[0]*temp_edge[4]+0.7071067811865475*facDiff[1]*temp_edge[2]; 
  edge_incr[5] = 0.632455532033676*facDiff[1]*temp_edge[13]+0.6324555320336759*facDiff[2]*temp_edge[5]+0.7071067811865475*facDiff[0]*temp_edge[5]+0.7071067811865475*facDiff[1]*temp_edge[3]; 
  edge_incr[6] = 0.7071067811865475*facDiff[2]*temp_edge[17]+0.7071067811865475*facDiff[1]*temp_edge[10]+0.7071067811865475*facDiff[0]*temp_edge[6]; 
  edge_incr[7] = 0.4517539514526256*facDiff[2]*temp_edge[7]+0.7071067811865475*facDiff[0]*temp_edge[7]+0.7071067811865475*temp_edge[0]*facDiff[2]+0.6324555320336759*facDiff[1]*temp_edge[1]; 
  edge_incr[8] = 0.7071067811865475*facDiff[1]*temp_edge[12]+0.7071067811865475*facDiff[0]*temp_edge[8]; 
  edge_incr[9] = 0.7071067811865475*facDiff[1]*temp_edge[15]+0.7071067811865475*facDiff[0]*temp_edge[9]; 
  edge_incr[10] = 0.6324555320336759*facDiff[1]*temp_edge[17]+0.6324555320336759*facDiff[2]*temp_edge[10]+0.7071067811865475*facDiff[0]*temp_edge[10]+0.7071067811865475*facDiff[1]*temp_edge[6]; 
  edge_incr[11] = 0.4517539514526256*facDiff[2]*temp_edge[11]+0.7071067811865475*facDiff[0]*temp_edge[11]+0.632455532033676*facDiff[1]*temp_edge[4]+0.7071067811865475*facDiff[2]*temp_edge[2]; 
  edge_incr[12] = 0.6324555320336759*facDiff[2]*temp_edge[12]+0.7071067811865475*facDiff[0]*temp_edge[12]+0.7071067811865475*facDiff[1]*temp_edge[8]; 
  edge_incr[13] = 0.4517539514526256*facDiff[2]*temp_edge[13]+0.7071067811865475*facDiff[0]*temp_edge[13]+0.632455532033676*facDiff[1]*temp_edge[5]+0.7071067811865475*facDiff[2]*temp_edge[3]; 
  edge_incr[14] = 0.7071067811865475*facDiff[1]*temp_edge[18]+0.7071067811865475*facDiff[0]*temp_edge[14]; 
  edge_incr[15] = 0.6324555320336759*facDiff[2]*temp_edge[15]+0.7071067811865475*facDiff[0]*temp_edge[15]+0.7071067811865475*facDiff[1]*temp_edge[9]; 
  edge_incr[16] = 0.7071067811865475*facDiff[1]*temp_edge[19]+0.7071067811865475*facDiff[0]*temp_edge[16]; 
  edge_incr[17] = 0.4517539514526256*facDiff[2]*temp_edge[17]+0.7071067811865475*facDiff[0]*temp_edge[17]+0.6324555320336759*facDiff[1]*temp_edge[10]+0.7071067811865475*facDiff[2]*temp_edge[6]; 
  edge_incr[18] = 0.6324555320336759*facDiff[2]*temp_edge[18]+0.7071067811865475*facDiff[0]*temp_edge[18]+0.7071067811865475*facDiff[1]*temp_edge[14]; 
  edge_incr[19] = 0.6324555320336759*facDiff[2]*temp_edge[19]+0.7071067811865475*facDiff[0]*temp_edge[19]+0.7071067811865475*facDiff[1]*temp_edge[16]; 

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
