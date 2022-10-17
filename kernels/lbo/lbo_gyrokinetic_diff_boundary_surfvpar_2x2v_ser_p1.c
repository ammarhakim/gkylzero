#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH void lbo_gyrokinetic_diff_boundary_surfvpar_2x2v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  // w[4]: Cell-center coordinates. 
  // dxv[4]: Cell spacing. 
  // m_: species mass.
  // bmag_inv: 1/(magnetic field magnitude). 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[8]: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // fskin/edge: Distribution function in cells 
  // out: Incremented distribution function in cell 

  const double *nuVtSqSum = &nuPrimMomsSum[4];

  double rdvSq4 = 4.0/(dxv[2]*dxv[2]); 

  double vol_incr[24] = {0.0};
  vol_incr[16] = 3.354101966249685*nuVtSqSum[3]*fskin[5]*rdvSq4+3.354101966249685*fskin[2]*nuVtSqSum[2]*rdvSq4+3.354101966249685*fskin[1]*nuVtSqSum[1]*rdvSq4+3.354101966249685*fskin[0]*nuVtSqSum[0]*rdvSq4; 
  vol_incr[17] = 3.354101966249684*nuVtSqSum[2]*fskin[5]*rdvSq4+3.354101966249684*fskin[2]*nuVtSqSum[3]*rdvSq4+3.354101966249684*fskin[0]*nuVtSqSum[1]*rdvSq4+3.354101966249684*nuVtSqSum[0]*fskin[1]*rdvSq4; 
  vol_incr[18] = 3.354101966249684*nuVtSqSum[1]*fskin[5]*rdvSq4+3.354101966249684*fskin[1]*nuVtSqSum[3]*rdvSq4+3.354101966249684*fskin[0]*nuVtSqSum[2]*rdvSq4+3.354101966249684*nuVtSqSum[0]*fskin[2]*rdvSq4; 
  vol_incr[19] = 3.354101966249684*nuVtSqSum[3]*fskin[12]*rdvSq4+3.354101966249684*nuVtSqSum[2]*fskin[9]*rdvSq4+3.354101966249684*nuVtSqSum[1]*fskin[8]*rdvSq4+3.354101966249684*nuVtSqSum[0]*fskin[4]*rdvSq4; 
  vol_incr[20] = 3.354101966249685*nuVtSqSum[0]*fskin[5]*rdvSq4+3.354101966249685*fskin[0]*nuVtSqSum[3]*rdvSq4+3.354101966249685*fskin[1]*nuVtSqSum[2]*rdvSq4+3.354101966249685*nuVtSqSum[1]*fskin[2]*rdvSq4; 
  vol_incr[21] = 3.354101966249685*nuVtSqSum[2]*fskin[12]*rdvSq4+3.354101966249685*nuVtSqSum[3]*fskin[9]*rdvSq4+3.354101966249685*nuVtSqSum[0]*fskin[8]*rdvSq4+3.354101966249685*nuVtSqSum[1]*fskin[4]*rdvSq4; 
  vol_incr[22] = 3.354101966249685*nuVtSqSum[1]*fskin[12]*rdvSq4+3.354101966249685*nuVtSqSum[0]*fskin[9]*rdvSq4+3.354101966249685*nuVtSqSum[3]*fskin[8]*rdvSq4+3.354101966249685*nuVtSqSum[2]*fskin[4]*rdvSq4; 
  vol_incr[23] = 3.354101966249684*nuVtSqSum[0]*fskin[12]*rdvSq4+3.354101966249684*nuVtSqSum[1]*fskin[9]*rdvSq4+3.354101966249684*nuVtSqSum[2]*fskin[8]*rdvSq4+3.354101966249684*nuVtSqSum[3]*fskin[4]*rdvSq4; 

  double temp_diff[24] = {0.0}; 
  double temp_edge[24] = {0.0}; 
  double diff_incr[24] = {0.0}; 
  double edge_incr[24] = {0.0}; 

  if (edge == -1) { 

  temp_diff[0] = (-0.6708203932499369*fskin[16])+0.6708203932499369*fedge[16]-1.190784930203603*fskin[3]-1.190784930203603*fedge[3]-0.9375*fskin[0]+0.9375*fedge[0]; 
  temp_diff[1] = (-0.6708203932499369*fskin[17])+0.6708203932499369*fedge[17]-1.190784930203603*fskin[6]-1.190784930203603*fedge[6]-0.9375*fskin[1]+0.9375*fedge[1]; 
  temp_diff[2] = (-0.6708203932499369*fskin[18])+0.6708203932499369*fedge[18]-1.190784930203603*fskin[7]-1.190784930203603*fedge[7]-0.9375*fskin[2]+0.9375*fedge[2]; 
  temp_diff[3] = (-1.585502557353661*fskin[16])+0.7382874503707888*fedge[16]-2.671875*fskin[3]-1.453125*fedge[3]-2.056810333988042*fskin[0]+1.190784930203603*fedge[0]; 
  temp_diff[4] = (-0.6708203932499369*fskin[19])+0.6708203932499369*fedge[19]-1.190784930203603*fskin[10]-1.190784930203603*fedge[10]-0.9375*fskin[4]+0.9375*fedge[4]; 
  temp_diff[5] = (-0.6708203932499369*fskin[20])+0.6708203932499369*fedge[20]-1.190784930203603*fskin[11]-1.190784930203603*fedge[11]-0.9375*fskin[5]+0.9375*fedge[5]; 
  temp_diff[6] = (-1.585502557353661*fskin[17])+0.7382874503707888*fedge[17]-2.671875*fskin[6]-1.453125*fedge[6]-2.056810333988042*fskin[1]+1.190784930203603*fedge[1]; 
  temp_diff[7] = (-1.585502557353661*fskin[18])+0.7382874503707888*fedge[18]-2.671875*fskin[7]-1.453125*fedge[7]-2.056810333988042*fskin[2]+1.190784930203603*fedge[2]; 
  temp_diff[8] = (-0.6708203932499369*fskin[21])+0.6708203932499369*fedge[21]-1.190784930203603*fskin[13]-1.190784930203603*fedge[13]-0.9375*fskin[8]+0.9375*fedge[8]; 
  temp_diff[9] = (-0.6708203932499369*fskin[22])+0.6708203932499369*fedge[22]-1.190784930203603*fskin[14]-1.190784930203603*fedge[14]-0.9375*fskin[9]+0.9375*fedge[9]; 
  temp_diff[10] = (-1.585502557353661*fskin[19])+0.7382874503707888*fedge[19]-2.671875*fskin[10]-1.453125*fedge[10]-2.056810333988042*fskin[4]+1.190784930203603*fedge[4]; 
  temp_diff[11] = (-1.585502557353661*fskin[20])+0.7382874503707888*fedge[20]-2.671875*fskin[11]-1.453125*fedge[11]-2.056810333988042*fskin[5]+1.190784930203603*fedge[5]; 
  temp_diff[12] = (-0.6708203932499369*fskin[23])+0.6708203932499369*fedge[23]-1.190784930203603*fskin[15]-1.190784930203603*fedge[15]-0.9375*fskin[12]+0.9375*fedge[12]; 
  temp_diff[13] = (-1.585502557353661*fskin[21])+0.7382874503707888*fedge[21]-2.671875*fskin[13]-1.453125*fedge[13]-2.056810333988042*fskin[8]+1.190784930203603*fedge[8]; 
  temp_diff[14] = (-1.585502557353661*fskin[22])+0.7382874503707888*fedge[22]-2.671875*fskin[14]-1.453125*fedge[14]-2.056810333988042*fskin[9]+1.190784930203603*fedge[9]; 
  temp_diff[15] = (-1.585502557353661*fskin[23])+0.7382874503707888*fedge[23]-2.671875*fskin[15]-1.453125*fedge[15]-2.056810333988042*fskin[12]+1.190784930203603*fedge[12]; 
  temp_diff[16] = (-3.140625*fskin[16])-0.140625*fedge[16]-5.022775277112744*fskin[3]-0.3025768239224545*fedge[3]-3.773364712030896*fskin[0]+0.4192627457812106*fedge[0]; 
  temp_diff[17] = (-3.140625*fskin[17])-0.140625*fedge[17]-5.022775277112744*fskin[6]-0.3025768239224544*fedge[6]-3.773364712030894*fskin[1]+0.4192627457812105*fedge[1]; 
  temp_diff[18] = (-3.140625*fskin[18])-0.140625*fedge[18]-5.022775277112744*fskin[7]-0.3025768239224544*fedge[7]-3.773364712030894*fskin[2]+0.4192627457812105*fedge[2]; 
  temp_diff[19] = (-3.140625*fskin[19])-0.140625*fedge[19]-5.022775277112744*fskin[10]-0.3025768239224544*fedge[10]-3.773364712030894*fskin[4]+0.4192627457812105*fedge[4]; 
  temp_diff[20] = (-3.140625*fskin[20])-0.140625*fedge[20]-5.022775277112744*fskin[11]-0.3025768239224545*fedge[11]-3.773364712030896*fskin[5]+0.4192627457812106*fedge[5]; 
  temp_diff[21] = (-3.140625*fskin[21])-0.140625*fedge[21]-5.022775277112744*fskin[13]-0.3025768239224545*fedge[13]-3.773364712030896*fskin[8]+0.4192627457812106*fedge[8]; 
  temp_diff[22] = (-3.140625*fskin[22])-0.140625*fedge[22]-5.022775277112744*fskin[14]-0.3025768239224545*fedge[14]-3.773364712030896*fskin[9]+0.4192627457812106*fedge[9]; 
  temp_diff[23] = (-3.140625*fskin[23])-0.140625*fedge[23]-5.022775277112744*fskin[15]-0.3025768239224544*fedge[15]-3.773364712030894*fskin[12]+0.4192627457812105*fedge[12]; 

  temp_edge[3] = 1.936491673103709*fskin[16]-1.5*fskin[3]+0.8660254037844386*fskin[0]; 
  temp_edge[6] = 1.936491673103709*fskin[17]-1.5*fskin[6]+0.8660254037844386*fskin[1]; 
  temp_edge[7] = 1.936491673103709*fskin[18]-1.5*fskin[7]+0.8660254037844386*fskin[2]; 
  temp_edge[10] = 1.936491673103709*fskin[19]-1.5*fskin[10]+0.8660254037844386*fskin[4]; 
  temp_edge[11] = 1.936491673103709*fskin[20]-1.5*fskin[11]+0.8660254037844386*fskin[5]; 
  temp_edge[13] = 1.936491673103709*fskin[21]-1.5*fskin[13]+0.8660254037844386*fskin[8]; 
  temp_edge[14] = 1.936491673103709*fskin[22]-1.5*fskin[14]+0.8660254037844386*fskin[9]; 
  temp_edge[15] = 1.936491673103709*fskin[23]-1.5*fskin[15]+0.8660254037844386*fskin[12]; 
  temp_edge[16] = (-7.5*fskin[16])+5.809475019311125*fskin[3]-3.354101966249685*fskin[0]; 
  temp_edge[17] = (-7.5*fskin[17])+5.809475019311126*fskin[6]-3.354101966249684*fskin[1]; 
  temp_edge[18] = (-7.5*fskin[18])+5.809475019311126*fskin[7]-3.354101966249684*fskin[2]; 
  temp_edge[19] = (-7.5*fskin[19])+5.809475019311126*fskin[10]-3.354101966249684*fskin[4]; 
  temp_edge[20] = (-7.5*fskin[20])+5.809475019311125*fskin[11]-3.354101966249685*fskin[5]; 
  temp_edge[21] = (-7.5*fskin[21])+5.809475019311125*fskin[13]-3.354101966249685*fskin[8]; 
  temp_edge[22] = (-7.5*fskin[22])+5.809475019311125*fskin[14]-3.354101966249685*fskin[9]; 
  temp_edge[23] = (-7.5*fskin[23])+5.809475019311126*fskin[15]-3.354101966249684*fskin[12]; 

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
  diff_incr[16] = 0.5*nuVtSqSum[3]*temp_diff[20]+0.5000000000000001*nuVtSqSum[2]*temp_diff[18]+0.5000000000000001*nuVtSqSum[1]*temp_diff[17]+0.5*nuVtSqSum[0]*temp_diff[16]; 
  diff_incr[17] = 0.5000000000000001*nuVtSqSum[2]*temp_diff[20]+0.5*nuVtSqSum[3]*temp_diff[18]+0.5*nuVtSqSum[0]*temp_diff[17]+0.5000000000000001*nuVtSqSum[1]*temp_diff[16]; 
  diff_incr[18] = 0.5000000000000001*nuVtSqSum[1]*temp_diff[20]+0.5*nuVtSqSum[0]*temp_diff[18]+0.5*nuVtSqSum[3]*temp_diff[17]+0.5000000000000001*nuVtSqSum[2]*temp_diff[16]; 
  diff_incr[19] = 0.5*nuVtSqSum[3]*temp_diff[23]+0.5000000000000001*nuVtSqSum[2]*temp_diff[22]+0.5000000000000001*nuVtSqSum[1]*temp_diff[21]+0.5*nuVtSqSum[0]*temp_diff[19]; 
  diff_incr[20] = 0.5*nuVtSqSum[0]*temp_diff[20]+0.5000000000000001*nuVtSqSum[1]*temp_diff[18]+0.5000000000000001*nuVtSqSum[2]*temp_diff[17]+0.5*nuVtSqSum[3]*temp_diff[16]; 
  diff_incr[21] = 0.5000000000000001*nuVtSqSum[2]*temp_diff[23]+0.5*nuVtSqSum[3]*temp_diff[22]+0.5*nuVtSqSum[0]*temp_diff[21]+0.5000000000000001*nuVtSqSum[1]*temp_diff[19]; 
  diff_incr[22] = 0.5000000000000001*nuVtSqSum[1]*temp_diff[23]+0.5*nuVtSqSum[0]*temp_diff[22]+0.5*nuVtSqSum[3]*temp_diff[21]+0.5000000000000001*nuVtSqSum[2]*temp_diff[19]; 
  diff_incr[23] = 0.5*nuVtSqSum[0]*temp_diff[23]+0.5000000000000001*nuVtSqSum[1]*temp_diff[22]+0.5000000000000001*nuVtSqSum[2]*temp_diff[21]+0.5*nuVtSqSum[3]*temp_diff[19]; 

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
  edge_incr[16] = 0.5*nuVtSqSum[3]*temp_edge[20]+0.5000000000000001*nuVtSqSum[2]*temp_edge[18]+0.5000000000000001*nuVtSqSum[1]*temp_edge[17]+0.5*nuVtSqSum[0]*temp_edge[16]; 
  edge_incr[17] = 0.5000000000000001*nuVtSqSum[2]*temp_edge[20]+0.5*nuVtSqSum[3]*temp_edge[18]+0.5*nuVtSqSum[0]*temp_edge[17]+0.5000000000000001*nuVtSqSum[1]*temp_edge[16]; 
  edge_incr[18] = 0.5000000000000001*nuVtSqSum[1]*temp_edge[20]+0.5*nuVtSqSum[0]*temp_edge[18]+0.5*nuVtSqSum[3]*temp_edge[17]+0.5000000000000001*nuVtSqSum[2]*temp_edge[16]; 
  edge_incr[19] = 0.5*nuVtSqSum[3]*temp_edge[23]+0.5000000000000001*nuVtSqSum[2]*temp_edge[22]+0.5000000000000001*nuVtSqSum[1]*temp_edge[21]+0.5*nuVtSqSum[0]*temp_edge[19]; 
  edge_incr[20] = 0.5*nuVtSqSum[0]*temp_edge[20]+0.5000000000000001*nuVtSqSum[1]*temp_edge[18]+0.5000000000000001*nuVtSqSum[2]*temp_edge[17]+0.5*nuVtSqSum[3]*temp_edge[16]; 
  edge_incr[21] = 0.5000000000000001*nuVtSqSum[2]*temp_edge[23]+0.5*nuVtSqSum[3]*temp_edge[22]+0.5*nuVtSqSum[0]*temp_edge[21]+0.5000000000000001*nuVtSqSum[1]*temp_edge[19]; 
  edge_incr[22] = 0.5000000000000001*nuVtSqSum[1]*temp_edge[23]+0.5*nuVtSqSum[0]*temp_edge[22]+0.5*nuVtSqSum[3]*temp_edge[21]+0.5000000000000001*nuVtSqSum[2]*temp_edge[19]; 
  edge_incr[23] = 0.5*nuVtSqSum[0]*temp_edge[23]+0.5000000000000001*nuVtSqSum[1]*temp_edge[22]+0.5000000000000001*nuVtSqSum[2]*temp_edge[21]+0.5*nuVtSqSum[3]*temp_edge[19]; 


  } else { 

  temp_diff[0] = (-0.6708203932499369*fskin[16])+0.6708203932499369*fedge[16]+1.190784930203603*fskin[3]+1.190784930203603*fedge[3]-0.9375*fskin[0]+0.9375*fedge[0]; 
  temp_diff[1] = (-0.6708203932499369*fskin[17])+0.6708203932499369*fedge[17]+1.190784930203603*fskin[6]+1.190784930203603*fedge[6]-0.9375*fskin[1]+0.9375*fedge[1]; 
  temp_diff[2] = (-0.6708203932499369*fskin[18])+0.6708203932499369*fedge[18]+1.190784930203603*fskin[7]+1.190784930203603*fedge[7]-0.9375*fskin[2]+0.9375*fedge[2]; 
  temp_diff[3] = 1.585502557353661*fskin[16]-0.7382874503707888*fedge[16]-2.671875*fskin[3]-1.453125*fedge[3]+2.056810333988042*fskin[0]-1.190784930203603*fedge[0]; 
  temp_diff[4] = (-0.6708203932499369*fskin[19])+0.6708203932499369*fedge[19]+1.190784930203603*fskin[10]+1.190784930203603*fedge[10]-0.9375*fskin[4]+0.9375*fedge[4]; 
  temp_diff[5] = (-0.6708203932499369*fskin[20])+0.6708203932499369*fedge[20]+1.190784930203603*fskin[11]+1.190784930203603*fedge[11]-0.9375*fskin[5]+0.9375*fedge[5]; 
  temp_diff[6] = 1.585502557353661*fskin[17]-0.7382874503707888*fedge[17]-2.671875*fskin[6]-1.453125*fedge[6]+2.056810333988042*fskin[1]-1.190784930203603*fedge[1]; 
  temp_diff[7] = 1.585502557353661*fskin[18]-0.7382874503707888*fedge[18]-2.671875*fskin[7]-1.453125*fedge[7]+2.056810333988042*fskin[2]-1.190784930203603*fedge[2]; 
  temp_diff[8] = (-0.6708203932499369*fskin[21])+0.6708203932499369*fedge[21]+1.190784930203603*fskin[13]+1.190784930203603*fedge[13]-0.9375*fskin[8]+0.9375*fedge[8]; 
  temp_diff[9] = (-0.6708203932499369*fskin[22])+0.6708203932499369*fedge[22]+1.190784930203603*fskin[14]+1.190784930203603*fedge[14]-0.9375*fskin[9]+0.9375*fedge[9]; 
  temp_diff[10] = 1.585502557353661*fskin[19]-0.7382874503707888*fedge[19]-2.671875*fskin[10]-1.453125*fedge[10]+2.056810333988042*fskin[4]-1.190784930203603*fedge[4]; 
  temp_diff[11] = 1.585502557353661*fskin[20]-0.7382874503707888*fedge[20]-2.671875*fskin[11]-1.453125*fedge[11]+2.056810333988042*fskin[5]-1.190784930203603*fedge[5]; 
  temp_diff[12] = (-0.6708203932499369*fskin[23])+0.6708203932499369*fedge[23]+1.190784930203603*fskin[15]+1.190784930203603*fedge[15]-0.9375*fskin[12]+0.9375*fedge[12]; 
  temp_diff[13] = 1.585502557353661*fskin[21]-0.7382874503707888*fedge[21]-2.671875*fskin[13]-1.453125*fedge[13]+2.056810333988042*fskin[8]-1.190784930203603*fedge[8]; 
  temp_diff[14] = 1.585502557353661*fskin[22]-0.7382874503707888*fedge[22]-2.671875*fskin[14]-1.453125*fedge[14]+2.056810333988042*fskin[9]-1.190784930203603*fedge[9]; 
  temp_diff[15] = 1.585502557353661*fskin[23]-0.7382874503707888*fedge[23]-2.671875*fskin[15]-1.453125*fedge[15]+2.056810333988042*fskin[12]-1.190784930203603*fedge[12]; 
  temp_diff[16] = (-3.140625*fskin[16])-0.140625*fedge[16]+5.022775277112744*fskin[3]+0.3025768239224545*fedge[3]-3.773364712030896*fskin[0]+0.4192627457812106*fedge[0]; 
  temp_diff[17] = (-3.140625*fskin[17])-0.140625*fedge[17]+5.022775277112744*fskin[6]+0.3025768239224544*fedge[6]-3.773364712030894*fskin[1]+0.4192627457812105*fedge[1]; 
  temp_diff[18] = (-3.140625*fskin[18])-0.140625*fedge[18]+5.022775277112744*fskin[7]+0.3025768239224544*fedge[7]-3.773364712030894*fskin[2]+0.4192627457812105*fedge[2]; 
  temp_diff[19] = (-3.140625*fskin[19])-0.140625*fedge[19]+5.022775277112744*fskin[10]+0.3025768239224544*fedge[10]-3.773364712030894*fskin[4]+0.4192627457812105*fedge[4]; 
  temp_diff[20] = (-3.140625*fskin[20])-0.140625*fedge[20]+5.022775277112744*fskin[11]+0.3025768239224545*fedge[11]-3.773364712030896*fskin[5]+0.4192627457812106*fedge[5]; 
  temp_diff[21] = (-3.140625*fskin[21])-0.140625*fedge[21]+5.022775277112744*fskin[13]+0.3025768239224545*fedge[13]-3.773364712030896*fskin[8]+0.4192627457812106*fedge[8]; 
  temp_diff[22] = (-3.140625*fskin[22])-0.140625*fedge[22]+5.022775277112744*fskin[14]+0.3025768239224545*fedge[14]-3.773364712030896*fskin[9]+0.4192627457812106*fedge[9]; 
  temp_diff[23] = (-3.140625*fskin[23])-0.140625*fedge[23]+5.022775277112744*fskin[15]+0.3025768239224544*fedge[15]-3.773364712030894*fskin[12]+0.4192627457812105*fedge[12]; 

  temp_edge[3] = (-1.936491673103709*fskin[16])-1.5*fskin[3]-0.8660254037844386*fskin[0]; 
  temp_edge[6] = (-1.936491673103709*fskin[17])-1.5*fskin[6]-0.8660254037844386*fskin[1]; 
  temp_edge[7] = (-1.936491673103709*fskin[18])-1.5*fskin[7]-0.8660254037844386*fskin[2]; 
  temp_edge[10] = (-1.936491673103709*fskin[19])-1.5*fskin[10]-0.8660254037844386*fskin[4]; 
  temp_edge[11] = (-1.936491673103709*fskin[20])-1.5*fskin[11]-0.8660254037844386*fskin[5]; 
  temp_edge[13] = (-1.936491673103709*fskin[21])-1.5*fskin[13]-0.8660254037844386*fskin[8]; 
  temp_edge[14] = (-1.936491673103709*fskin[22])-1.5*fskin[14]-0.8660254037844386*fskin[9]; 
  temp_edge[15] = (-1.936491673103709*fskin[23])-1.5*fskin[15]-0.8660254037844386*fskin[12]; 
  temp_edge[16] = (-7.5*fskin[16])-5.809475019311125*fskin[3]-3.354101966249685*fskin[0]; 
  temp_edge[17] = (-7.5*fskin[17])-5.809475019311126*fskin[6]-3.354101966249684*fskin[1]; 
  temp_edge[18] = (-7.5*fskin[18])-5.809475019311126*fskin[7]-3.354101966249684*fskin[2]; 
  temp_edge[19] = (-7.5*fskin[19])-5.809475019311126*fskin[10]-3.354101966249684*fskin[4]; 
  temp_edge[20] = (-7.5*fskin[20])-5.809475019311125*fskin[11]-3.354101966249685*fskin[5]; 
  temp_edge[21] = (-7.5*fskin[21])-5.809475019311125*fskin[13]-3.354101966249685*fskin[8]; 
  temp_edge[22] = (-7.5*fskin[22])-5.809475019311125*fskin[14]-3.354101966249685*fskin[9]; 
  temp_edge[23] = (-7.5*fskin[23])-5.809475019311126*fskin[15]-3.354101966249684*fskin[12]; 

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
  diff_incr[16] = 0.5*nuVtSqSum[3]*temp_diff[20]+0.5000000000000001*nuVtSqSum[2]*temp_diff[18]+0.5000000000000001*nuVtSqSum[1]*temp_diff[17]+0.5*nuVtSqSum[0]*temp_diff[16]; 
  diff_incr[17] = 0.5000000000000001*nuVtSqSum[2]*temp_diff[20]+0.5*nuVtSqSum[3]*temp_diff[18]+0.5*nuVtSqSum[0]*temp_diff[17]+0.5000000000000001*nuVtSqSum[1]*temp_diff[16]; 
  diff_incr[18] = 0.5000000000000001*nuVtSqSum[1]*temp_diff[20]+0.5*nuVtSqSum[0]*temp_diff[18]+0.5*nuVtSqSum[3]*temp_diff[17]+0.5000000000000001*nuVtSqSum[2]*temp_diff[16]; 
  diff_incr[19] = 0.5*nuVtSqSum[3]*temp_diff[23]+0.5000000000000001*nuVtSqSum[2]*temp_diff[22]+0.5000000000000001*nuVtSqSum[1]*temp_diff[21]+0.5*nuVtSqSum[0]*temp_diff[19]; 
  diff_incr[20] = 0.5*nuVtSqSum[0]*temp_diff[20]+0.5000000000000001*nuVtSqSum[1]*temp_diff[18]+0.5000000000000001*nuVtSqSum[2]*temp_diff[17]+0.5*nuVtSqSum[3]*temp_diff[16]; 
  diff_incr[21] = 0.5000000000000001*nuVtSqSum[2]*temp_diff[23]+0.5*nuVtSqSum[3]*temp_diff[22]+0.5*nuVtSqSum[0]*temp_diff[21]+0.5000000000000001*nuVtSqSum[1]*temp_diff[19]; 
  diff_incr[22] = 0.5000000000000001*nuVtSqSum[1]*temp_diff[23]+0.5*nuVtSqSum[0]*temp_diff[22]+0.5*nuVtSqSum[3]*temp_diff[21]+0.5000000000000001*nuVtSqSum[2]*temp_diff[19]; 
  diff_incr[23] = 0.5*nuVtSqSum[0]*temp_diff[23]+0.5000000000000001*nuVtSqSum[1]*temp_diff[22]+0.5000000000000001*nuVtSqSum[2]*temp_diff[21]+0.5*nuVtSqSum[3]*temp_diff[19]; 

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
  edge_incr[16] = 0.5*nuVtSqSum[3]*temp_edge[20]+0.5000000000000001*nuVtSqSum[2]*temp_edge[18]+0.5000000000000001*nuVtSqSum[1]*temp_edge[17]+0.5*nuVtSqSum[0]*temp_edge[16]; 
  edge_incr[17] = 0.5000000000000001*nuVtSqSum[2]*temp_edge[20]+0.5*nuVtSqSum[3]*temp_edge[18]+0.5*nuVtSqSum[0]*temp_edge[17]+0.5000000000000001*nuVtSqSum[1]*temp_edge[16]; 
  edge_incr[18] = 0.5000000000000001*nuVtSqSum[1]*temp_edge[20]+0.5*nuVtSqSum[0]*temp_edge[18]+0.5*nuVtSqSum[3]*temp_edge[17]+0.5000000000000001*nuVtSqSum[2]*temp_edge[16]; 
  edge_incr[19] = 0.5*nuVtSqSum[3]*temp_edge[23]+0.5000000000000001*nuVtSqSum[2]*temp_edge[22]+0.5000000000000001*nuVtSqSum[1]*temp_edge[21]+0.5*nuVtSqSum[0]*temp_edge[19]; 
  edge_incr[20] = 0.5*nuVtSqSum[0]*temp_edge[20]+0.5000000000000001*nuVtSqSum[1]*temp_edge[18]+0.5000000000000001*nuVtSqSum[2]*temp_edge[17]+0.5*nuVtSqSum[3]*temp_edge[16]; 
  edge_incr[21] = 0.5000000000000001*nuVtSqSum[2]*temp_edge[23]+0.5*nuVtSqSum[3]*temp_edge[22]+0.5*nuVtSqSum[0]*temp_edge[21]+0.5000000000000001*nuVtSqSum[1]*temp_edge[19]; 
  edge_incr[22] = 0.5000000000000001*nuVtSqSum[1]*temp_edge[23]+0.5*nuVtSqSum[0]*temp_edge[22]+0.5*nuVtSqSum[3]*temp_edge[21]+0.5000000000000001*nuVtSqSum[2]*temp_edge[19]; 
  edge_incr[23] = 0.5*nuVtSqSum[0]*temp_edge[23]+0.5000000000000001*nuVtSqSum[1]*temp_edge[22]+0.5000000000000001*nuVtSqSum[2]*temp_edge[21]+0.5*nuVtSqSum[3]*temp_edge[19]; 

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
} 
