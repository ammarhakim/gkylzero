#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH void lbo_gyrokinetic_diff_boundary_surfmu_1x2v_tensor_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  // w[3]: Cell-center coordinates. 
  // dxv[3]: Cell spacing. 
  // m_: species mass.
  // bmag_inv: 1/(magnetic field magnitude). 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[6]: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // fskin/edge: Distribution function in cells 
  // out: Incremented distribution function in cell 

  const double *nuVtSqSum = &nuPrimMomsSum[3];

  double rdvSq4 = 4.0/(dxv[2]*dxv[2]); 

  double facDiff[3]; 
  // Expand diffusion coefficient in conf basis.
  facDiff[0] = 1.414213562373095*(bmag_inv[2]*nuVtSqSum[2]+bmag_inv[1]*nuVtSqSum[1]+bmag_inv[0]*nuVtSqSum[0])*m_; 
  facDiff[1] = (1.264911064067352*(bmag_inv[1]*nuVtSqSum[2]+nuVtSqSum[1]*bmag_inv[2])+1.414213562373095*(bmag_inv[0]*nuVtSqSum[1]+nuVtSqSum[0]*bmag_inv[1]))*m_; 
  facDiff[2] = (0.9035079029052515*bmag_inv[2]*nuVtSqSum[2]+1.414213562373095*(bmag_inv[0]*nuVtSqSum[2]+nuVtSqSum[0]*bmag_inv[2])+1.264911064067352*bmag_inv[1]*nuVtSqSum[1])*m_; 

  double vol_incr[27] = {0.0};
  vol_incr[3] = 0.6123724356957944*dxv[2]*facDiff[2]*fskin[7]+0.6123724356957944*facDiff[1]*fskin[1]*dxv[2]+0.6123724356957944*facDiff[0]*fskin[0]*dxv[2]; 
  vol_incr[5] = 0.5477225575051661*facDiff[1]*dxv[2]*fskin[7]+0.5477225575051661*fskin[1]*dxv[2]*facDiff[2]+0.6123724356957944*facDiff[0]*fskin[1]*dxv[2]+0.6123724356957944*fskin[0]*facDiff[1]*dxv[2]; 
  vol_incr[6] = 0.6123724356957944*dxv[2]*facDiff[2]*fskin[11]+0.6123724356957944*facDiff[1]*dxv[2]*fskin[4]+0.6123724356957944*facDiff[0]*dxv[2]*fskin[2]; 
  vol_incr[9] = 2.738612787525831*dxv[2]*facDiff[2]*fskin[13]+4.743416490252569*facDiff[2]*w[2]*fskin[7]+2.738612787525831*facDiff[1]*dxv[2]*fskin[5]+2.738612787525831*facDiff[0]*dxv[2]*fskin[3]+4.743416490252569*facDiff[1]*fskin[1]*w[2]+4.743416490252569*facDiff[0]*fskin[0]*w[2]; 
  vol_incr[10] = 0.5477225575051661*facDiff[1]*dxv[2]*fskin[11]+0.5477225575051661*dxv[2]*facDiff[2]*fskin[4]+0.6123724356957944*facDiff[0]*dxv[2]*fskin[4]+0.6123724356957944*facDiff[1]*dxv[2]*fskin[2]; 
  vol_incr[13] = 0.3912303982179757*dxv[2]*facDiff[2]*fskin[7]+0.6123724356957944*facDiff[0]*dxv[2]*fskin[7]+0.6123724356957944*fskin[0]*dxv[2]*facDiff[2]+0.5477225575051661*facDiff[1]*fskin[1]*dxv[2]; 
  vol_incr[14] = 0.6123724356957944*dxv[2]*facDiff[2]*fskin[20]+0.6123724356957944*facDiff[1]*dxv[2]*fskin[12]+0.6123724356957944*facDiff[0]*dxv[2]*fskin[8]; 
  vol_incr[15] = 2.449489742783178*facDiff[1]*dxv[2]*fskin[13]+4.242640687119286*facDiff[1]*w[2]*fskin[7]+2.449489742783178*dxv[2]*facDiff[2]*fskin[5]+2.738612787525831*facDiff[0]*dxv[2]*fskin[5]+2.738612787525831*facDiff[1]*dxv[2]*fskin[3]+4.242640687119286*fskin[1]*facDiff[2]*w[2]+4.743416490252569*facDiff[0]*fskin[1]*w[2]+4.743416490252569*fskin[0]*facDiff[1]*w[2]; 
  vol_incr[16] = 2.738612787525831*dxv[2]*facDiff[2]*fskin[17]+4.743416490252569*facDiff[2]*w[2]*fskin[11]+2.738612787525831*facDiff[1]*dxv[2]*fskin[10]+2.738612787525831*facDiff[0]*dxv[2]*fskin[6]+4.743416490252569*facDiff[1]*w[2]*fskin[4]+4.743416490252569*facDiff[0]*fskin[2]*w[2]; 
  vol_incr[17] = 0.3912303982179757*dxv[2]*facDiff[2]*fskin[11]+0.6123724356957944*facDiff[0]*dxv[2]*fskin[11]+0.5477225575051661*facDiff[1]*dxv[2]*fskin[4]+0.6123724356957944*dxv[2]*facDiff[2]*fskin[2]; 
  vol_incr[18] = 0.5477225575051661*facDiff[1]*dxv[2]*fskin[20]+0.5477225575051661*dxv[2]*facDiff[2]*fskin[12]+0.6123724356957944*facDiff[0]*dxv[2]*fskin[12]+0.6123724356957944*facDiff[1]*dxv[2]*fskin[8]; 
  vol_incr[19] = 2.449489742783178*facDiff[1]*dxv[2]*fskin[17]+4.242640687119286*facDiff[1]*w[2]*fskin[11]+2.449489742783178*dxv[2]*facDiff[2]*fskin[10]+2.738612787525831*facDiff[0]*dxv[2]*fskin[10]+2.738612787525831*facDiff[1]*dxv[2]*fskin[6]+4.242640687119286*facDiff[2]*w[2]*fskin[4]+4.743416490252569*facDiff[0]*w[2]*fskin[4]+4.743416490252569*facDiff[1]*fskin[2]*w[2]; 
  vol_incr[21] = 1.749635530559413*dxv[2]*facDiff[2]*fskin[13]+2.738612787525831*facDiff[0]*dxv[2]*fskin[13]+3.030457633656632*facDiff[2]*w[2]*fskin[7]+4.743416490252569*facDiff[0]*w[2]*fskin[7]+2.449489742783178*facDiff[1]*dxv[2]*fskin[5]+2.738612787525831*dxv[2]*facDiff[2]*fskin[3]+4.743416490252569*fskin[0]*facDiff[2]*w[2]+4.242640687119286*facDiff[1]*fskin[1]*w[2]; 
  vol_incr[22] = 2.738612787525831*dxv[2]*facDiff[2]*fskin[23]+4.743416490252569*facDiff[2]*w[2]*fskin[20]+2.738612787525831*facDiff[1]*dxv[2]*fskin[18]+2.738612787525831*facDiff[0]*dxv[2]*fskin[14]+4.743416490252569*facDiff[1]*w[2]*fskin[12]+4.743416490252569*facDiff[0]*w[2]*fskin[8]; 
  vol_incr[23] = 0.3912303982179757*dxv[2]*facDiff[2]*fskin[20]+0.6123724356957944*facDiff[0]*dxv[2]*fskin[20]+0.5477225575051661*facDiff[1]*dxv[2]*fskin[12]+0.6123724356957944*dxv[2]*facDiff[2]*fskin[8]; 
  vol_incr[24] = 1.749635530559413*dxv[2]*facDiff[2]*fskin[17]+2.738612787525831*facDiff[0]*dxv[2]*fskin[17]+3.030457633656633*facDiff[2]*w[2]*fskin[11]+4.743416490252569*facDiff[0]*w[2]*fskin[11]+2.449489742783178*facDiff[1]*dxv[2]*fskin[10]+2.738612787525831*dxv[2]*facDiff[2]*fskin[6]+4.242640687119286*facDiff[1]*w[2]*fskin[4]+4.743416490252569*facDiff[2]*fskin[2]*w[2]; 
  vol_incr[25] = 2.449489742783178*facDiff[1]*dxv[2]*fskin[23]+4.242640687119286*facDiff[1]*w[2]*fskin[20]+2.449489742783178*dxv[2]*facDiff[2]*fskin[18]+2.738612787525831*facDiff[0]*dxv[2]*fskin[18]+2.738612787525831*facDiff[1]*dxv[2]*fskin[14]+4.242640687119286*facDiff[2]*w[2]*fskin[12]+4.743416490252569*facDiff[0]*w[2]*fskin[12]+4.743416490252569*facDiff[1]*w[2]*fskin[8]; 
  vol_incr[26] = 1.749635530559413*dxv[2]*facDiff[2]*fskin[23]+2.738612787525831*facDiff[0]*dxv[2]*fskin[23]+3.030457633656632*facDiff[2]*w[2]*fskin[20]+4.743416490252569*facDiff[0]*w[2]*fskin[20]+2.449489742783178*facDiff[1]*dxv[2]*fskin[18]+2.738612787525831*dxv[2]*facDiff[2]*fskin[14]+4.242640687119286*facDiff[1]*w[2]*fskin[12]+4.743416490252569*facDiff[2]*w[2]*fskin[8]; 

  double edgeSurf_incr[27] = {0.0}; 
  double boundSurf_incr[27] = {0.0}; 

  if (edge == -1) { 

  double edgeSurf[27] = {0.0}; 
  edgeSurf[0] = (-0.6708203932499369*w[2]*fskin[9])-0.3354101966249685*dxv[2]*fskin[9]+0.6708203932499369*w[2]*fedge[9]+0.3354101966249685*dxv[2]*fedge[9]-1.190784930203603*w[2]*fskin[3]-0.5953924651018015*dxv[2]*fskin[3]-1.190784930203603*w[2]*fedge[3]-0.5953924651018015*dxv[2]*fedge[3]-0.9375*fskin[0]*w[2]+0.9375*fedge[0]*w[2]-0.46875*fskin[0]*dxv[2]+0.46875*fedge[0]*dxv[2]; 
  edgeSurf[1] = (-0.6708203932499369*w[2]*fskin[15])-0.3354101966249685*dxv[2]*fskin[15]+0.6708203932499369*w[2]*fedge[15]+0.3354101966249685*dxv[2]*fedge[15]-1.190784930203603*w[2]*fskin[5]-0.5953924651018015*dxv[2]*fskin[5]-1.190784930203603*w[2]*fedge[5]-0.5953924651018015*dxv[2]*fedge[5]-0.9375*fskin[1]*w[2]+0.9375*fedge[1]*w[2]-0.46875*fskin[1]*dxv[2]+0.46875*fedge[1]*dxv[2]; 
  edgeSurf[2] = (-0.6708203932499369*w[2]*fskin[16])-0.3354101966249685*dxv[2]*fskin[16]+0.6708203932499369*w[2]*fedge[16]+0.3354101966249685*dxv[2]*fedge[16]-1.190784930203603*w[2]*fskin[6]-0.5953924651018015*dxv[2]*fskin[6]-1.190784930203603*w[2]*fedge[6]-0.5953924651018015*dxv[2]*fedge[6]-0.9375*fskin[2]*w[2]+0.9375*fedge[2]*w[2]-0.46875*dxv[2]*fskin[2]+0.46875*dxv[2]*fedge[2]; 
  edgeSurf[3] = (-1.585502557353661*w[2]*fskin[9])-0.7927512786768306*dxv[2]*fskin[9]+0.7382874503707888*w[2]*fedge[9]+0.3691437251853944*dxv[2]*fedge[9]-2.671875*w[2]*fskin[3]-1.3359375*dxv[2]*fskin[3]-1.453125*w[2]*fedge[3]-0.7265625*dxv[2]*fedge[3]-2.056810333988042*fskin[0]*w[2]+1.190784930203603*fedge[0]*w[2]-1.028405166994021*fskin[0]*dxv[2]+0.5953924651018015*fedge[0]*dxv[2]; 
  edgeSurf[4] = (-0.6708203932499369*w[2]*fskin[19])-0.3354101966249685*dxv[2]*fskin[19]+0.6708203932499369*w[2]*fedge[19]+0.3354101966249685*dxv[2]*fedge[19]-1.190784930203603*w[2]*fskin[10]-0.5953924651018015*dxv[2]*fskin[10]-1.190784930203603*w[2]*fedge[10]-0.5953924651018015*dxv[2]*fedge[10]-0.9375*w[2]*fskin[4]-0.46875*dxv[2]*fskin[4]+0.9375*w[2]*fedge[4]+0.46875*dxv[2]*fedge[4]; 
  edgeSurf[5] = (-1.585502557353661*w[2]*fskin[15])-0.7927512786768306*dxv[2]*fskin[15]+0.7382874503707888*w[2]*fedge[15]+0.3691437251853944*dxv[2]*fedge[15]-2.671875*w[2]*fskin[5]-1.3359375*dxv[2]*fskin[5]-1.453125*w[2]*fedge[5]-0.7265625*dxv[2]*fedge[5]-2.056810333988042*fskin[1]*w[2]+1.190784930203603*fedge[1]*w[2]-1.028405166994021*fskin[1]*dxv[2]+0.5953924651018015*fedge[1]*dxv[2]; 
  edgeSurf[6] = (-1.585502557353661*w[2]*fskin[16])-0.7927512786768306*dxv[2]*fskin[16]+0.7382874503707888*w[2]*fedge[16]+0.3691437251853944*dxv[2]*fedge[16]-2.671875*w[2]*fskin[6]-1.3359375*dxv[2]*fskin[6]-1.453125*w[2]*fedge[6]-0.7265625*dxv[2]*fedge[6]-2.056810333988042*fskin[2]*w[2]+1.190784930203603*fedge[2]*w[2]-1.028405166994021*dxv[2]*fskin[2]+0.5953924651018015*dxv[2]*fedge[2]; 
  edgeSurf[7] = (-0.6708203932499369*w[2]*fskin[21])-0.3354101966249685*dxv[2]*fskin[21]+0.6708203932499369*w[2]*fedge[21]+0.3354101966249685*dxv[2]*fedge[21]-1.190784930203603*w[2]*fskin[13]-0.5953924651018015*dxv[2]*fskin[13]-1.190784930203603*w[2]*fedge[13]-0.5953924651018015*dxv[2]*fedge[13]-0.9375*w[2]*fskin[7]-0.46875*dxv[2]*fskin[7]+0.9375*w[2]*fedge[7]+0.46875*dxv[2]*fedge[7]; 
  edgeSurf[8] = (-0.6708203932499369*w[2]*fskin[22])-0.3354101966249685*dxv[2]*fskin[22]+0.6708203932499369*w[2]*fedge[22]+0.3354101966249685*dxv[2]*fedge[22]-1.190784930203603*w[2]*fskin[14]-0.5953924651018015*dxv[2]*fskin[14]-1.190784930203603*w[2]*fedge[14]-0.5953924651018015*dxv[2]*fedge[14]-0.9375*w[2]*fskin[8]-0.46875*dxv[2]*fskin[8]+0.9375*w[2]*fedge[8]+0.46875*dxv[2]*fedge[8]; 
  edgeSurf[9] = (-3.140625*w[2]*fskin[9])-1.5703125*dxv[2]*fskin[9]-0.140625*w[2]*fedge[9]-0.0703125*dxv[2]*fedge[9]-5.022775277112744*w[2]*fskin[3]-2.511387638556372*dxv[2]*fskin[3]-0.3025768239224545*w[2]*fedge[3]-0.1512884119612272*dxv[2]*fedge[3]-3.773364712030896*fskin[0]*w[2]+0.4192627457812106*fedge[0]*w[2]-1.886682356015448*fskin[0]*dxv[2]+0.2096313728906053*fedge[0]*dxv[2]; 
  edgeSurf[10] = (-1.585502557353661*w[2]*fskin[19])-0.7927512786768306*dxv[2]*fskin[19]+0.7382874503707888*w[2]*fedge[19]+0.3691437251853944*dxv[2]*fedge[19]-2.671875*w[2]*fskin[10]-1.3359375*dxv[2]*fskin[10]-1.453125*w[2]*fedge[10]-0.7265625*dxv[2]*fedge[10]-2.056810333988042*w[2]*fskin[4]-1.028405166994021*dxv[2]*fskin[4]+1.190784930203603*w[2]*fedge[4]+0.5953924651018015*dxv[2]*fedge[4]; 
  edgeSurf[11] = (-0.6708203932499369*w[2]*fskin[24])-0.3354101966249685*dxv[2]*fskin[24]+0.6708203932499369*w[2]*fedge[24]+0.3354101966249685*dxv[2]*fedge[24]-1.190784930203603*w[2]*fskin[17]-0.5953924651018015*dxv[2]*fskin[17]-1.190784930203603*w[2]*fedge[17]-0.5953924651018015*dxv[2]*fedge[17]-0.9375*w[2]*fskin[11]-0.46875*dxv[2]*fskin[11]+0.9375*w[2]*fedge[11]+0.46875*dxv[2]*fedge[11]; 
  edgeSurf[12] = (-0.6708203932499369*w[2]*fskin[25])-0.3354101966249685*dxv[2]*fskin[25]+0.6708203932499369*w[2]*fedge[25]+0.3354101966249685*dxv[2]*fedge[25]-1.190784930203603*w[2]*fskin[18]-0.5953924651018015*dxv[2]*fskin[18]-1.190784930203603*w[2]*fedge[18]-0.5953924651018015*dxv[2]*fedge[18]-0.9375*w[2]*fskin[12]-0.46875*dxv[2]*fskin[12]+0.9375*w[2]*fedge[12]+0.46875*dxv[2]*fedge[12]; 
  edgeSurf[13] = (-1.585502557353661*w[2]*fskin[21])-0.7927512786768306*dxv[2]*fskin[21]+0.7382874503707888*w[2]*fedge[21]+0.3691437251853944*dxv[2]*fedge[21]-2.671875*w[2]*fskin[13]-1.3359375*dxv[2]*fskin[13]-1.453125*w[2]*fedge[13]-0.7265625*dxv[2]*fedge[13]-2.056810333988042*w[2]*fskin[7]-1.028405166994021*dxv[2]*fskin[7]+1.190784930203603*w[2]*fedge[7]+0.5953924651018015*dxv[2]*fedge[7]; 
  edgeSurf[14] = (-1.585502557353661*w[2]*fskin[22])-0.7927512786768306*dxv[2]*fskin[22]+0.7382874503707888*w[2]*fedge[22]+0.3691437251853944*dxv[2]*fedge[22]-2.671875*w[2]*fskin[14]-1.3359375*dxv[2]*fskin[14]-1.453125*w[2]*fedge[14]-0.7265625*dxv[2]*fedge[14]-2.056810333988042*w[2]*fskin[8]-1.028405166994021*dxv[2]*fskin[8]+1.190784930203603*w[2]*fedge[8]+0.5953924651018015*dxv[2]*fedge[8]; 
  edgeSurf[15] = (-3.140625*w[2]*fskin[15])-1.5703125*dxv[2]*fskin[15]-0.140625*w[2]*fedge[15]-0.0703125*dxv[2]*fedge[15]-5.022775277112744*w[2]*fskin[5]-2.511387638556372*dxv[2]*fskin[5]-0.3025768239224544*w[2]*fedge[5]-0.1512884119612272*dxv[2]*fedge[5]-3.773364712030894*fskin[1]*w[2]+0.4192627457812105*fedge[1]*w[2]-1.886682356015447*fskin[1]*dxv[2]+0.2096313728906053*fedge[1]*dxv[2]; 
  edgeSurf[16] = (-3.140625*w[2]*fskin[16])-1.5703125*dxv[2]*fskin[16]-0.140625*w[2]*fedge[16]-0.0703125*dxv[2]*fedge[16]-5.022775277112744*w[2]*fskin[6]-2.511387638556372*dxv[2]*fskin[6]-0.3025768239224544*w[2]*fedge[6]-0.1512884119612272*dxv[2]*fedge[6]-3.773364712030894*fskin[2]*w[2]+0.4192627457812105*fedge[2]*w[2]-1.886682356015447*dxv[2]*fskin[2]+0.2096313728906053*dxv[2]*fedge[2]; 
  edgeSurf[17] = (-1.585502557353661*w[2]*fskin[24])-0.7927512786768306*dxv[2]*fskin[24]+0.7382874503707888*w[2]*fedge[24]+0.3691437251853944*dxv[2]*fedge[24]-2.671875*w[2]*fskin[17]-1.3359375*dxv[2]*fskin[17]-1.453125*w[2]*fedge[17]-0.7265625*dxv[2]*fedge[17]-2.056810333988042*w[2]*fskin[11]-1.028405166994021*dxv[2]*fskin[11]+1.190784930203603*w[2]*fedge[11]+0.5953924651018015*dxv[2]*fedge[11]; 
  edgeSurf[18] = (-1.585502557353661*w[2]*fskin[25])-0.7927512786768306*dxv[2]*fskin[25]+0.7382874503707888*w[2]*fedge[25]+0.3691437251853944*dxv[2]*fedge[25]-2.671875*w[2]*fskin[18]-1.3359375*dxv[2]*fskin[18]-1.453125*w[2]*fedge[18]-0.7265625*dxv[2]*fedge[18]-2.056810333988042*w[2]*fskin[12]-1.028405166994021*dxv[2]*fskin[12]+1.190784930203603*w[2]*fedge[12]+0.5953924651018015*dxv[2]*fedge[12]; 
  edgeSurf[19] = (-3.140625*w[2]*fskin[19])-1.5703125*dxv[2]*fskin[19]-0.140625*w[2]*fedge[19]-0.0703125*dxv[2]*fedge[19]-5.022775277112744*w[2]*fskin[10]-2.511387638556372*dxv[2]*fskin[10]-0.3025768239224545*w[2]*fedge[10]-0.1512884119612272*dxv[2]*fedge[10]-3.773364712030896*w[2]*fskin[4]-1.886682356015448*dxv[2]*fskin[4]+0.4192627457812106*w[2]*fedge[4]+0.2096313728906053*dxv[2]*fedge[4]; 
  edgeSurf[20] = (-0.6708203932499369*w[2]*fskin[26])-0.3354101966249685*dxv[2]*fskin[26]+0.6708203932499369*w[2]*fedge[26]+0.3354101966249685*dxv[2]*fedge[26]-1.190784930203603*w[2]*fskin[23]-0.5953924651018015*dxv[2]*fskin[23]-1.190784930203603*w[2]*fedge[23]-0.5953924651018015*dxv[2]*fedge[23]-0.9375*w[2]*fskin[20]-0.46875*dxv[2]*fskin[20]+0.9375*w[2]*fedge[20]+0.46875*dxv[2]*fedge[20]; 
  edgeSurf[21] = (-3.140625*w[2]*fskin[21])-1.5703125*dxv[2]*fskin[21]-0.140625*w[2]*fedge[21]-0.0703125*dxv[2]*fedge[21]-5.022775277112744*w[2]*fskin[13]-2.511387638556372*dxv[2]*fskin[13]-0.3025768239224544*w[2]*fedge[13]-0.1512884119612272*dxv[2]*fedge[13]-3.773364712030896*w[2]*fskin[7]-1.886682356015448*dxv[2]*fskin[7]+0.4192627457812106*w[2]*fedge[7]+0.2096313728906053*dxv[2]*fedge[7]; 
  edgeSurf[22] = (-3.140625*w[2]*fskin[22])-1.5703125*dxv[2]*fskin[22]-0.140625*w[2]*fedge[22]-0.0703125*dxv[2]*fedge[22]-5.022775277112744*w[2]*fskin[14]-2.511387638556372*dxv[2]*fskin[14]-0.3025768239224544*w[2]*fedge[14]-0.1512884119612272*dxv[2]*fedge[14]-3.773364712030896*w[2]*fskin[8]-1.886682356015448*dxv[2]*fskin[8]+0.4192627457812106*w[2]*fedge[8]+0.2096313728906053*dxv[2]*fedge[8]; 
  edgeSurf[23] = (-1.585502557353661*w[2]*fskin[26])-0.7927512786768306*dxv[2]*fskin[26]+0.7382874503707888*w[2]*fedge[26]+0.3691437251853944*dxv[2]*fedge[26]-2.671875*w[2]*fskin[23]-1.3359375*dxv[2]*fskin[23]-1.453125*w[2]*fedge[23]-0.7265625*dxv[2]*fedge[23]-2.056810333988042*w[2]*fskin[20]-1.028405166994021*dxv[2]*fskin[20]+1.190784930203603*w[2]*fedge[20]+0.5953924651018015*dxv[2]*fedge[20]; 
  edgeSurf[24] = (-3.140625*w[2]*fskin[24])-1.5703125*dxv[2]*fskin[24]-0.140625*w[2]*fedge[24]-0.0703125*dxv[2]*fedge[24]-5.022775277112744*w[2]*fskin[17]-2.511387638556372*dxv[2]*fskin[17]-0.3025768239224545*w[2]*fedge[17]-0.1512884119612272*dxv[2]*fedge[17]-3.773364712030894*w[2]*fskin[11]-1.886682356015447*dxv[2]*fskin[11]+0.4192627457812105*w[2]*fedge[11]+0.2096313728906053*dxv[2]*fedge[11]; 
  edgeSurf[25] = (-3.140625*w[2]*fskin[25])-1.5703125*dxv[2]*fskin[25]-0.140625*w[2]*fedge[25]-0.0703125*dxv[2]*fedge[25]-5.022775277112744*w[2]*fskin[18]-2.511387638556372*dxv[2]*fskin[18]-0.3025768239224545*w[2]*fedge[18]-0.1512884119612272*dxv[2]*fedge[18]-3.773364712030894*w[2]*fskin[12]-1.886682356015447*dxv[2]*fskin[12]+0.4192627457812105*w[2]*fedge[12]+0.2096313728906053*dxv[2]*fedge[12]; 
  edgeSurf[26] = (-3.140625*w[2]*fskin[26])-1.5703125*dxv[2]*fskin[26]-0.140625*w[2]*fedge[26]-0.0703125*dxv[2]*fedge[26]-5.022775277112744*w[2]*fskin[23]-2.511387638556372*dxv[2]*fskin[23]-0.3025768239224545*w[2]*fedge[23]-0.1512884119612272*dxv[2]*fedge[23]-3.773364712030896*w[2]*fskin[20]-1.886682356015448*dxv[2]*fskin[20]+0.4192627457812106*w[2]*fedge[20]+0.2096313728906053*dxv[2]*fedge[20]; 

  double boundSurf[27] = {0.0}; 
  boundSurf[3] = 1.936491673103709*w[2]*fskin[9]-0.9682458365518543*dxv[2]*fskin[9]-1.5*w[2]*fskin[3]+0.75*dxv[2]*fskin[3]+0.8660254037844386*fskin[0]*w[2]-0.4330127018922193*fskin[0]*dxv[2]; 
  boundSurf[5] = 1.936491673103709*w[2]*fskin[15]-0.9682458365518543*dxv[2]*fskin[15]-1.5*w[2]*fskin[5]+0.75*dxv[2]*fskin[5]+0.8660254037844386*fskin[1]*w[2]-0.4330127018922193*fskin[1]*dxv[2]; 
  boundSurf[6] = 1.936491673103709*w[2]*fskin[16]-0.9682458365518543*dxv[2]*fskin[16]-1.5*w[2]*fskin[6]+0.75*dxv[2]*fskin[6]+0.8660254037844386*fskin[2]*w[2]-0.4330127018922193*dxv[2]*fskin[2]; 
  boundSurf[9] = (-7.5*w[2]*fskin[9])+3.75*dxv[2]*fskin[9]+5.809475019311125*w[2]*fskin[3]-2.904737509655563*dxv[2]*fskin[3]-3.354101966249685*fskin[0]*w[2]+1.677050983124842*fskin[0]*dxv[2]; 
  boundSurf[10] = 1.936491673103709*w[2]*fskin[19]-0.9682458365518543*dxv[2]*fskin[19]-1.5*w[2]*fskin[10]+0.75*dxv[2]*fskin[10]+0.8660254037844386*w[2]*fskin[4]-0.4330127018922193*dxv[2]*fskin[4]; 
  boundSurf[13] = 1.936491673103709*w[2]*fskin[21]-0.9682458365518543*dxv[2]*fskin[21]-1.5*w[2]*fskin[13]+0.75*dxv[2]*fskin[13]+0.8660254037844387*w[2]*fskin[7]-0.4330127018922194*dxv[2]*fskin[7]; 
  boundSurf[14] = 1.936491673103709*w[2]*fskin[22]-0.9682458365518543*dxv[2]*fskin[22]-1.5*w[2]*fskin[14]+0.75*dxv[2]*fskin[14]+0.8660254037844387*w[2]*fskin[8]-0.4330127018922194*dxv[2]*fskin[8]; 
  boundSurf[15] = (-7.5*w[2]*fskin[15])+3.75*dxv[2]*fskin[15]+5.809475019311126*w[2]*fskin[5]-2.904737509655563*dxv[2]*fskin[5]-3.354101966249684*fskin[1]*w[2]+1.677050983124842*fskin[1]*dxv[2]; 
  boundSurf[16] = (-7.5*w[2]*fskin[16])+3.75*dxv[2]*fskin[16]+5.809475019311126*w[2]*fskin[6]-2.904737509655563*dxv[2]*fskin[6]-3.354101966249684*fskin[2]*w[2]+1.677050983124842*dxv[2]*fskin[2]; 
  boundSurf[17] = 1.936491673103709*w[2]*fskin[24]-0.9682458365518543*dxv[2]*fskin[24]-1.5*w[2]*fskin[17]+0.75*dxv[2]*fskin[17]+0.8660254037844387*w[2]*fskin[11]-0.4330127018922194*dxv[2]*fskin[11]; 
  boundSurf[18] = 1.936491673103709*w[2]*fskin[25]-0.9682458365518543*dxv[2]*fskin[25]-1.5*w[2]*fskin[18]+0.75*dxv[2]*fskin[18]+0.8660254037844387*w[2]*fskin[12]-0.4330127018922194*dxv[2]*fskin[12]; 
  boundSurf[19] = (-7.5*w[2]*fskin[19])+3.75*dxv[2]*fskin[19]+5.809475019311125*w[2]*fskin[10]-2.904737509655563*dxv[2]*fskin[10]-3.354101966249685*w[2]*fskin[4]+1.677050983124842*dxv[2]*fskin[4]; 
  boundSurf[21] = (-7.5*w[2]*fskin[21])+3.75*dxv[2]*fskin[21]+5.809475019311126*w[2]*fskin[13]-2.904737509655563*dxv[2]*fskin[13]-3.354101966249685*w[2]*fskin[7]+1.677050983124842*dxv[2]*fskin[7]; 
  boundSurf[22] = (-7.5*w[2]*fskin[22])+3.75*dxv[2]*fskin[22]+5.809475019311126*w[2]*fskin[14]-2.904737509655563*dxv[2]*fskin[14]-3.354101966249685*w[2]*fskin[8]+1.677050983124842*dxv[2]*fskin[8]; 
  boundSurf[23] = 1.936491673103709*w[2]*fskin[26]-0.9682458365518543*dxv[2]*fskin[26]-1.5*w[2]*fskin[23]+0.75*dxv[2]*fskin[23]+0.8660254037844386*w[2]*fskin[20]-0.4330127018922193*dxv[2]*fskin[20]; 
  boundSurf[24] = (-7.5*w[2]*fskin[24])+3.75*dxv[2]*fskin[24]+5.809475019311125*w[2]*fskin[17]-2.904737509655563*dxv[2]*fskin[17]-3.354101966249684*w[2]*fskin[11]+1.677050983124842*dxv[2]*fskin[11]; 
  boundSurf[25] = (-7.5*w[2]*fskin[25])+3.75*dxv[2]*fskin[25]+5.809475019311125*w[2]*fskin[18]-2.904737509655563*dxv[2]*fskin[18]-3.354101966249684*w[2]*fskin[12]+1.677050983124842*dxv[2]*fskin[12]; 
  boundSurf[26] = (-7.5*w[2]*fskin[26])+3.75*dxv[2]*fskin[26]+5.809475019311125*w[2]*fskin[23]-2.904737509655563*dxv[2]*fskin[23]-3.354101966249685*w[2]*fskin[20]+1.677050983124842*dxv[2]*fskin[20]; 

  edgeSurf_incr[0] = 0.7071067811865475*facDiff[2]*edgeSurf[7]+0.7071067811865475*edgeSurf[1]*facDiff[1]+0.7071067811865475*edgeSurf[0]*facDiff[0]; 
  edgeSurf_incr[1] = 0.6324555320336759*facDiff[1]*edgeSurf[7]+0.6324555320336759*edgeSurf[1]*facDiff[2]+0.7071067811865475*edgeSurf[0]*facDiff[1]+0.7071067811865475*facDiff[0]*edgeSurf[1]; 
  edgeSurf_incr[2] = 0.7071067811865475*facDiff[2]*edgeSurf[11]+0.7071067811865475*facDiff[1]*edgeSurf[4]+0.7071067811865475*facDiff[0]*edgeSurf[2]; 
  edgeSurf_incr[3] = 0.7071067811865475*facDiff[2]*edgeSurf[13]+0.7071067811865475*facDiff[1]*edgeSurf[5]+0.7071067811865475*facDiff[0]*edgeSurf[3]; 
  edgeSurf_incr[4] = 0.632455532033676*facDiff[1]*edgeSurf[11]+0.6324555320336759*facDiff[2]*edgeSurf[4]+0.7071067811865475*facDiff[0]*edgeSurf[4]+0.7071067811865475*facDiff[1]*edgeSurf[2]; 
  edgeSurf_incr[5] = 0.632455532033676*facDiff[1]*edgeSurf[13]+0.6324555320336759*facDiff[2]*edgeSurf[5]+0.7071067811865475*facDiff[0]*edgeSurf[5]+0.7071067811865475*facDiff[1]*edgeSurf[3]; 
  edgeSurf_incr[6] = 0.7071067811865475*facDiff[2]*edgeSurf[17]+0.7071067811865475*facDiff[1]*edgeSurf[10]+0.7071067811865475*facDiff[0]*edgeSurf[6]; 
  edgeSurf_incr[7] = 0.4517539514526256*facDiff[2]*edgeSurf[7]+0.7071067811865475*facDiff[0]*edgeSurf[7]+0.7071067811865475*edgeSurf[0]*facDiff[2]+0.6324555320336759*edgeSurf[1]*facDiff[1]; 
  edgeSurf_incr[8] = 0.7071067811865475*facDiff[2]*edgeSurf[20]+0.7071067811865475*facDiff[1]*edgeSurf[12]+0.7071067811865475*facDiff[0]*edgeSurf[8]; 
  edgeSurf_incr[9] = 0.7071067811865475*facDiff[2]*edgeSurf[21]+0.7071067811865475*facDiff[1]*edgeSurf[15]+0.7071067811865475*facDiff[0]*edgeSurf[9]; 
  edgeSurf_incr[10] = 0.6324555320336759*facDiff[1]*edgeSurf[17]+0.6324555320336759*facDiff[2]*edgeSurf[10]+0.7071067811865475*facDiff[0]*edgeSurf[10]+0.7071067811865475*facDiff[1]*edgeSurf[6]; 
  edgeSurf_incr[11] = 0.4517539514526256*facDiff[2]*edgeSurf[11]+0.7071067811865475*facDiff[0]*edgeSurf[11]+0.632455532033676*facDiff[1]*edgeSurf[4]+0.7071067811865475*edgeSurf[2]*facDiff[2]; 
  edgeSurf_incr[12] = 0.632455532033676*facDiff[1]*edgeSurf[20]+0.6324555320336759*facDiff[2]*edgeSurf[12]+0.7071067811865475*facDiff[0]*edgeSurf[12]+0.7071067811865475*facDiff[1]*edgeSurf[8]; 
  edgeSurf_incr[13] = 0.4517539514526256*facDiff[2]*edgeSurf[13]+0.7071067811865475*facDiff[0]*edgeSurf[13]+0.632455532033676*facDiff[1]*edgeSurf[5]+0.7071067811865475*facDiff[2]*edgeSurf[3]; 
  edgeSurf_incr[14] = 0.7071067811865475*facDiff[2]*edgeSurf[23]+0.7071067811865475*facDiff[1]*edgeSurf[18]+0.7071067811865475*facDiff[0]*edgeSurf[14]; 
  edgeSurf_incr[15] = 0.632455532033676*facDiff[1]*edgeSurf[21]+0.6324555320336759*facDiff[2]*edgeSurf[15]+0.7071067811865475*facDiff[0]*edgeSurf[15]+0.7071067811865475*facDiff[1]*edgeSurf[9]; 
  edgeSurf_incr[16] = 0.7071067811865475*facDiff[2]*edgeSurf[24]+0.7071067811865475*facDiff[1]*edgeSurf[19]+0.7071067811865475*facDiff[0]*edgeSurf[16]; 
  edgeSurf_incr[17] = 0.4517539514526256*facDiff[2]*edgeSurf[17]+0.7071067811865475*facDiff[0]*edgeSurf[17]+0.6324555320336759*facDiff[1]*edgeSurf[10]+0.7071067811865475*facDiff[2]*edgeSurf[6]; 
  edgeSurf_incr[18] = 0.6324555320336759*facDiff[1]*edgeSurf[23]+0.6324555320336759*facDiff[2]*edgeSurf[18]+0.7071067811865475*facDiff[0]*edgeSurf[18]+0.7071067811865475*facDiff[1]*edgeSurf[14]; 
  edgeSurf_incr[19] = 0.6324555320336759*facDiff[1]*edgeSurf[24]+0.6324555320336759*facDiff[2]*edgeSurf[19]+0.7071067811865475*facDiff[0]*edgeSurf[19]+0.7071067811865475*facDiff[1]*edgeSurf[16]; 
  edgeSurf_incr[20] = 0.4517539514526256*facDiff[2]*edgeSurf[20]+0.7071067811865475*facDiff[0]*edgeSurf[20]+0.632455532033676*facDiff[1]*edgeSurf[12]+0.7071067811865475*facDiff[2]*edgeSurf[8]; 
  edgeSurf_incr[21] = 0.4517539514526256*facDiff[2]*edgeSurf[21]+0.7071067811865475*facDiff[0]*edgeSurf[21]+0.632455532033676*facDiff[1]*edgeSurf[15]+0.7071067811865475*facDiff[2]*edgeSurf[9]; 
  edgeSurf_incr[22] = 0.7071067811865475*facDiff[2]*edgeSurf[26]+0.7071067811865475*facDiff[1]*edgeSurf[25]+0.7071067811865475*facDiff[0]*edgeSurf[22]; 
  edgeSurf_incr[23] = 0.4517539514526256*facDiff[2]*edgeSurf[23]+0.7071067811865475*facDiff[0]*edgeSurf[23]+0.6324555320336759*facDiff[1]*edgeSurf[18]+0.7071067811865475*facDiff[2]*edgeSurf[14]; 
  edgeSurf_incr[24] = 0.4517539514526256*facDiff[2]*edgeSurf[24]+0.7071067811865475*facDiff[0]*edgeSurf[24]+0.6324555320336759*facDiff[1]*edgeSurf[19]+0.7071067811865475*facDiff[2]*edgeSurf[16]; 
  edgeSurf_incr[25] = 0.6324555320336759*facDiff[1]*edgeSurf[26]+0.6324555320336759*facDiff[2]*edgeSurf[25]+0.7071067811865475*facDiff[0]*edgeSurf[25]+0.7071067811865475*facDiff[1]*edgeSurf[22]; 
  edgeSurf_incr[26] = 0.4517539514526256*facDiff[2]*edgeSurf[26]+0.7071067811865475*facDiff[0]*edgeSurf[26]+0.6324555320336759*facDiff[1]*edgeSurf[25]+0.7071067811865475*facDiff[2]*edgeSurf[22]; 

  boundSurf_incr[0] = 0.7071067811865475*facDiff[2]*boundSurf[7]+0.7071067811865475*boundSurf[1]*facDiff[1]+0.7071067811865475*boundSurf[0]*facDiff[0]; 
  boundSurf_incr[1] = 0.6324555320336759*facDiff[1]*boundSurf[7]+0.6324555320336759*boundSurf[1]*facDiff[2]+0.7071067811865475*boundSurf[0]*facDiff[1]+0.7071067811865475*facDiff[0]*boundSurf[1]; 
  boundSurf_incr[2] = 0.7071067811865475*facDiff[2]*boundSurf[11]+0.7071067811865475*facDiff[1]*boundSurf[4]+0.7071067811865475*facDiff[0]*boundSurf[2]; 
  boundSurf_incr[3] = 0.7071067811865475*facDiff[2]*boundSurf[13]+0.7071067811865475*facDiff[1]*boundSurf[5]+0.7071067811865475*facDiff[0]*boundSurf[3]; 
  boundSurf_incr[4] = 0.632455532033676*facDiff[1]*boundSurf[11]+0.6324555320336759*facDiff[2]*boundSurf[4]+0.7071067811865475*facDiff[0]*boundSurf[4]+0.7071067811865475*facDiff[1]*boundSurf[2]; 
  boundSurf_incr[5] = 0.632455532033676*facDiff[1]*boundSurf[13]+0.6324555320336759*facDiff[2]*boundSurf[5]+0.7071067811865475*facDiff[0]*boundSurf[5]+0.7071067811865475*facDiff[1]*boundSurf[3]; 
  boundSurf_incr[6] = 0.7071067811865475*facDiff[2]*boundSurf[17]+0.7071067811865475*facDiff[1]*boundSurf[10]+0.7071067811865475*facDiff[0]*boundSurf[6]; 
  boundSurf_incr[7] = 0.4517539514526256*facDiff[2]*boundSurf[7]+0.7071067811865475*facDiff[0]*boundSurf[7]+0.7071067811865475*boundSurf[0]*facDiff[2]+0.6324555320336759*boundSurf[1]*facDiff[1]; 
  boundSurf_incr[8] = 0.7071067811865475*facDiff[2]*boundSurf[20]+0.7071067811865475*facDiff[1]*boundSurf[12]+0.7071067811865475*facDiff[0]*boundSurf[8]; 
  boundSurf_incr[9] = 0.7071067811865475*facDiff[2]*boundSurf[21]+0.7071067811865475*facDiff[1]*boundSurf[15]+0.7071067811865475*facDiff[0]*boundSurf[9]; 
  boundSurf_incr[10] = 0.6324555320336759*facDiff[1]*boundSurf[17]+0.6324555320336759*facDiff[2]*boundSurf[10]+0.7071067811865475*facDiff[0]*boundSurf[10]+0.7071067811865475*facDiff[1]*boundSurf[6]; 
  boundSurf_incr[11] = 0.4517539514526256*facDiff[2]*boundSurf[11]+0.7071067811865475*facDiff[0]*boundSurf[11]+0.632455532033676*facDiff[1]*boundSurf[4]+0.7071067811865475*boundSurf[2]*facDiff[2]; 
  boundSurf_incr[12] = 0.632455532033676*facDiff[1]*boundSurf[20]+0.6324555320336759*facDiff[2]*boundSurf[12]+0.7071067811865475*facDiff[0]*boundSurf[12]+0.7071067811865475*facDiff[1]*boundSurf[8]; 
  boundSurf_incr[13] = 0.4517539514526256*facDiff[2]*boundSurf[13]+0.7071067811865475*facDiff[0]*boundSurf[13]+0.632455532033676*facDiff[1]*boundSurf[5]+0.7071067811865475*facDiff[2]*boundSurf[3]; 
  boundSurf_incr[14] = 0.7071067811865475*facDiff[2]*boundSurf[23]+0.7071067811865475*facDiff[1]*boundSurf[18]+0.7071067811865475*facDiff[0]*boundSurf[14]; 
  boundSurf_incr[15] = 0.632455532033676*facDiff[1]*boundSurf[21]+0.6324555320336759*facDiff[2]*boundSurf[15]+0.7071067811865475*facDiff[0]*boundSurf[15]+0.7071067811865475*facDiff[1]*boundSurf[9]; 
  boundSurf_incr[16] = 0.7071067811865475*facDiff[2]*boundSurf[24]+0.7071067811865475*facDiff[1]*boundSurf[19]+0.7071067811865475*facDiff[0]*boundSurf[16]; 
  boundSurf_incr[17] = 0.4517539514526256*facDiff[2]*boundSurf[17]+0.7071067811865475*facDiff[0]*boundSurf[17]+0.6324555320336759*facDiff[1]*boundSurf[10]+0.7071067811865475*facDiff[2]*boundSurf[6]; 
  boundSurf_incr[18] = 0.6324555320336759*facDiff[1]*boundSurf[23]+0.6324555320336759*facDiff[2]*boundSurf[18]+0.7071067811865475*facDiff[0]*boundSurf[18]+0.7071067811865475*facDiff[1]*boundSurf[14]; 
  boundSurf_incr[19] = 0.6324555320336759*facDiff[1]*boundSurf[24]+0.6324555320336759*facDiff[2]*boundSurf[19]+0.7071067811865475*facDiff[0]*boundSurf[19]+0.7071067811865475*facDiff[1]*boundSurf[16]; 
  boundSurf_incr[20] = 0.4517539514526256*facDiff[2]*boundSurf[20]+0.7071067811865475*facDiff[0]*boundSurf[20]+0.632455532033676*facDiff[1]*boundSurf[12]+0.7071067811865475*facDiff[2]*boundSurf[8]; 
  boundSurf_incr[21] = 0.4517539514526256*facDiff[2]*boundSurf[21]+0.7071067811865475*facDiff[0]*boundSurf[21]+0.632455532033676*facDiff[1]*boundSurf[15]+0.7071067811865475*facDiff[2]*boundSurf[9]; 
  boundSurf_incr[22] = 0.7071067811865475*facDiff[2]*boundSurf[26]+0.7071067811865475*facDiff[1]*boundSurf[25]+0.7071067811865475*facDiff[0]*boundSurf[22]; 
  boundSurf_incr[23] = 0.4517539514526256*facDiff[2]*boundSurf[23]+0.7071067811865475*facDiff[0]*boundSurf[23]+0.6324555320336759*facDiff[1]*boundSurf[18]+0.7071067811865475*facDiff[2]*boundSurf[14]; 
  boundSurf_incr[24] = 0.4517539514526256*facDiff[2]*boundSurf[24]+0.7071067811865475*facDiff[0]*boundSurf[24]+0.6324555320336759*facDiff[1]*boundSurf[19]+0.7071067811865475*facDiff[2]*boundSurf[16]; 
  boundSurf_incr[25] = 0.6324555320336759*facDiff[1]*boundSurf[26]+0.6324555320336759*facDiff[2]*boundSurf[25]+0.7071067811865475*facDiff[0]*boundSurf[25]+0.7071067811865475*facDiff[1]*boundSurf[22]; 
  boundSurf_incr[26] = 0.4517539514526256*facDiff[2]*boundSurf[26]+0.7071067811865475*facDiff[0]*boundSurf[26]+0.6324555320336759*facDiff[1]*boundSurf[25]+0.7071067811865475*facDiff[2]*boundSurf[22]; 


  } else { 

  double edgeSurf[27] = {0.0}; 
  edgeSurf[0] = (-0.6708203932499369*w[2]*fskin[9])+0.3354101966249685*dxv[2]*fskin[9]+0.6708203932499369*w[2]*fedge[9]-0.3354101966249685*dxv[2]*fedge[9]+1.190784930203603*w[2]*fskin[3]-0.5953924651018015*dxv[2]*fskin[3]+1.190784930203603*w[2]*fedge[3]-0.5953924651018015*dxv[2]*fedge[3]-0.9375*fskin[0]*w[2]+0.9375*fedge[0]*w[2]+0.46875*fskin[0]*dxv[2]-0.46875*fedge[0]*dxv[2]; 
  edgeSurf[1] = (-0.6708203932499369*w[2]*fskin[15])+0.3354101966249685*dxv[2]*fskin[15]+0.6708203932499369*w[2]*fedge[15]-0.3354101966249685*dxv[2]*fedge[15]+1.190784930203603*w[2]*fskin[5]-0.5953924651018015*dxv[2]*fskin[5]+1.190784930203603*w[2]*fedge[5]-0.5953924651018015*dxv[2]*fedge[5]-0.9375*fskin[1]*w[2]+0.9375*fedge[1]*w[2]+0.46875*fskin[1]*dxv[2]-0.46875*fedge[1]*dxv[2]; 
  edgeSurf[2] = (-0.6708203932499369*w[2]*fskin[16])+0.3354101966249685*dxv[2]*fskin[16]+0.6708203932499369*w[2]*fedge[16]-0.3354101966249685*dxv[2]*fedge[16]+1.190784930203603*w[2]*fskin[6]-0.5953924651018015*dxv[2]*fskin[6]+1.190784930203603*w[2]*fedge[6]-0.5953924651018015*dxv[2]*fedge[6]-0.9375*fskin[2]*w[2]+0.9375*fedge[2]*w[2]+0.46875*dxv[2]*fskin[2]-0.46875*dxv[2]*fedge[2]; 
  edgeSurf[3] = 1.585502557353661*w[2]*fskin[9]-0.7927512786768306*dxv[2]*fskin[9]-0.7382874503707888*w[2]*fedge[9]+0.3691437251853944*dxv[2]*fedge[9]-2.671875*w[2]*fskin[3]+1.3359375*dxv[2]*fskin[3]-1.453125*w[2]*fedge[3]+0.7265625*dxv[2]*fedge[3]+2.056810333988042*fskin[0]*w[2]-1.190784930203603*fedge[0]*w[2]-1.028405166994021*fskin[0]*dxv[2]+0.5953924651018015*fedge[0]*dxv[2]; 
  edgeSurf[4] = (-0.6708203932499369*w[2]*fskin[19])+0.3354101966249685*dxv[2]*fskin[19]+0.6708203932499369*w[2]*fedge[19]-0.3354101966249685*dxv[2]*fedge[19]+1.190784930203603*w[2]*fskin[10]-0.5953924651018015*dxv[2]*fskin[10]+1.190784930203603*w[2]*fedge[10]-0.5953924651018015*dxv[2]*fedge[10]-0.9375*w[2]*fskin[4]+0.46875*dxv[2]*fskin[4]+0.9375*w[2]*fedge[4]-0.46875*dxv[2]*fedge[4]; 
  edgeSurf[5] = 1.585502557353661*w[2]*fskin[15]-0.7927512786768306*dxv[2]*fskin[15]-0.7382874503707888*w[2]*fedge[15]+0.3691437251853944*dxv[2]*fedge[15]-2.671875*w[2]*fskin[5]+1.3359375*dxv[2]*fskin[5]-1.453125*w[2]*fedge[5]+0.7265625*dxv[2]*fedge[5]+2.056810333988042*fskin[1]*w[2]-1.190784930203603*fedge[1]*w[2]-1.028405166994021*fskin[1]*dxv[2]+0.5953924651018015*fedge[1]*dxv[2]; 
  edgeSurf[6] = 1.585502557353661*w[2]*fskin[16]-0.7927512786768306*dxv[2]*fskin[16]-0.7382874503707888*w[2]*fedge[16]+0.3691437251853944*dxv[2]*fedge[16]-2.671875*w[2]*fskin[6]+1.3359375*dxv[2]*fskin[6]-1.453125*w[2]*fedge[6]+0.7265625*dxv[2]*fedge[6]+2.056810333988042*fskin[2]*w[2]-1.190784930203603*fedge[2]*w[2]-1.028405166994021*dxv[2]*fskin[2]+0.5953924651018015*dxv[2]*fedge[2]; 
  edgeSurf[7] = (-0.6708203932499369*w[2]*fskin[21])+0.3354101966249685*dxv[2]*fskin[21]+0.6708203932499369*w[2]*fedge[21]-0.3354101966249685*dxv[2]*fedge[21]+1.190784930203603*w[2]*fskin[13]-0.5953924651018015*dxv[2]*fskin[13]+1.190784930203603*w[2]*fedge[13]-0.5953924651018015*dxv[2]*fedge[13]-0.9375*w[2]*fskin[7]+0.46875*dxv[2]*fskin[7]+0.9375*w[2]*fedge[7]-0.46875*dxv[2]*fedge[7]; 
  edgeSurf[8] = (-0.6708203932499369*w[2]*fskin[22])+0.3354101966249685*dxv[2]*fskin[22]+0.6708203932499369*w[2]*fedge[22]-0.3354101966249685*dxv[2]*fedge[22]+1.190784930203603*w[2]*fskin[14]-0.5953924651018015*dxv[2]*fskin[14]+1.190784930203603*w[2]*fedge[14]-0.5953924651018015*dxv[2]*fedge[14]-0.9375*w[2]*fskin[8]+0.46875*dxv[2]*fskin[8]+0.9375*w[2]*fedge[8]-0.46875*dxv[2]*fedge[8]; 
  edgeSurf[9] = (-3.140625*w[2]*fskin[9])+1.5703125*dxv[2]*fskin[9]-0.140625*w[2]*fedge[9]+0.0703125*dxv[2]*fedge[9]+5.022775277112744*w[2]*fskin[3]-2.511387638556372*dxv[2]*fskin[3]+0.3025768239224545*w[2]*fedge[3]-0.1512884119612272*dxv[2]*fedge[3]-3.773364712030896*fskin[0]*w[2]+0.4192627457812106*fedge[0]*w[2]+1.886682356015448*fskin[0]*dxv[2]-0.2096313728906053*fedge[0]*dxv[2]; 
  edgeSurf[10] = 1.585502557353661*w[2]*fskin[19]-0.7927512786768306*dxv[2]*fskin[19]-0.7382874503707888*w[2]*fedge[19]+0.3691437251853944*dxv[2]*fedge[19]-2.671875*w[2]*fskin[10]+1.3359375*dxv[2]*fskin[10]-1.453125*w[2]*fedge[10]+0.7265625*dxv[2]*fedge[10]+2.056810333988042*w[2]*fskin[4]-1.028405166994021*dxv[2]*fskin[4]-1.190784930203603*w[2]*fedge[4]+0.5953924651018015*dxv[2]*fedge[4]; 
  edgeSurf[11] = (-0.6708203932499369*w[2]*fskin[24])+0.3354101966249685*dxv[2]*fskin[24]+0.6708203932499369*w[2]*fedge[24]-0.3354101966249685*dxv[2]*fedge[24]+1.190784930203603*w[2]*fskin[17]-0.5953924651018015*dxv[2]*fskin[17]+1.190784930203603*w[2]*fedge[17]-0.5953924651018015*dxv[2]*fedge[17]-0.9375*w[2]*fskin[11]+0.46875*dxv[2]*fskin[11]+0.9375*w[2]*fedge[11]-0.46875*dxv[2]*fedge[11]; 
  edgeSurf[12] = (-0.6708203932499369*w[2]*fskin[25])+0.3354101966249685*dxv[2]*fskin[25]+0.6708203932499369*w[2]*fedge[25]-0.3354101966249685*dxv[2]*fedge[25]+1.190784930203603*w[2]*fskin[18]-0.5953924651018015*dxv[2]*fskin[18]+1.190784930203603*w[2]*fedge[18]-0.5953924651018015*dxv[2]*fedge[18]-0.9375*w[2]*fskin[12]+0.46875*dxv[2]*fskin[12]+0.9375*w[2]*fedge[12]-0.46875*dxv[2]*fedge[12]; 
  edgeSurf[13] = 1.585502557353661*w[2]*fskin[21]-0.7927512786768306*dxv[2]*fskin[21]-0.7382874503707888*w[2]*fedge[21]+0.3691437251853944*dxv[2]*fedge[21]-2.671875*w[2]*fskin[13]+1.3359375*dxv[2]*fskin[13]-1.453125*w[2]*fedge[13]+0.7265625*dxv[2]*fedge[13]+2.056810333988042*w[2]*fskin[7]-1.028405166994021*dxv[2]*fskin[7]-1.190784930203603*w[2]*fedge[7]+0.5953924651018015*dxv[2]*fedge[7]; 
  edgeSurf[14] = 1.585502557353661*w[2]*fskin[22]-0.7927512786768306*dxv[2]*fskin[22]-0.7382874503707888*w[2]*fedge[22]+0.3691437251853944*dxv[2]*fedge[22]-2.671875*w[2]*fskin[14]+1.3359375*dxv[2]*fskin[14]-1.453125*w[2]*fedge[14]+0.7265625*dxv[2]*fedge[14]+2.056810333988042*w[2]*fskin[8]-1.028405166994021*dxv[2]*fskin[8]-1.190784930203603*w[2]*fedge[8]+0.5953924651018015*dxv[2]*fedge[8]; 
  edgeSurf[15] = (-3.140625*w[2]*fskin[15])+1.5703125*dxv[2]*fskin[15]-0.140625*w[2]*fedge[15]+0.0703125*dxv[2]*fedge[15]+5.022775277112744*w[2]*fskin[5]-2.511387638556372*dxv[2]*fskin[5]+0.3025768239224544*w[2]*fedge[5]-0.1512884119612272*dxv[2]*fedge[5]-3.773364712030894*fskin[1]*w[2]+0.4192627457812105*fedge[1]*w[2]+1.886682356015447*fskin[1]*dxv[2]-0.2096313728906053*fedge[1]*dxv[2]; 
  edgeSurf[16] = (-3.140625*w[2]*fskin[16])+1.5703125*dxv[2]*fskin[16]-0.140625*w[2]*fedge[16]+0.0703125*dxv[2]*fedge[16]+5.022775277112744*w[2]*fskin[6]-2.511387638556372*dxv[2]*fskin[6]+0.3025768239224544*w[2]*fedge[6]-0.1512884119612272*dxv[2]*fedge[6]-3.773364712030894*fskin[2]*w[2]+0.4192627457812105*fedge[2]*w[2]+1.886682356015447*dxv[2]*fskin[2]-0.2096313728906053*dxv[2]*fedge[2]; 
  edgeSurf[17] = 1.585502557353661*w[2]*fskin[24]-0.7927512786768306*dxv[2]*fskin[24]-0.7382874503707888*w[2]*fedge[24]+0.3691437251853944*dxv[2]*fedge[24]-2.671875*w[2]*fskin[17]+1.3359375*dxv[2]*fskin[17]-1.453125*w[2]*fedge[17]+0.7265625*dxv[2]*fedge[17]+2.056810333988042*w[2]*fskin[11]-1.028405166994021*dxv[2]*fskin[11]-1.190784930203603*w[2]*fedge[11]+0.5953924651018015*dxv[2]*fedge[11]; 
  edgeSurf[18] = 1.585502557353661*w[2]*fskin[25]-0.7927512786768306*dxv[2]*fskin[25]-0.7382874503707888*w[2]*fedge[25]+0.3691437251853944*dxv[2]*fedge[25]-2.671875*w[2]*fskin[18]+1.3359375*dxv[2]*fskin[18]-1.453125*w[2]*fedge[18]+0.7265625*dxv[2]*fedge[18]+2.056810333988042*w[2]*fskin[12]-1.028405166994021*dxv[2]*fskin[12]-1.190784930203603*w[2]*fedge[12]+0.5953924651018015*dxv[2]*fedge[12]; 
  edgeSurf[19] = (-3.140625*w[2]*fskin[19])+1.5703125*dxv[2]*fskin[19]-0.140625*w[2]*fedge[19]+0.0703125*dxv[2]*fedge[19]+5.022775277112744*w[2]*fskin[10]-2.511387638556372*dxv[2]*fskin[10]+0.3025768239224545*w[2]*fedge[10]-0.1512884119612272*dxv[2]*fedge[10]-3.773364712030896*w[2]*fskin[4]+1.886682356015448*dxv[2]*fskin[4]+0.4192627457812106*w[2]*fedge[4]-0.2096313728906053*dxv[2]*fedge[4]; 
  edgeSurf[20] = (-0.6708203932499369*w[2]*fskin[26])+0.3354101966249685*dxv[2]*fskin[26]+0.6708203932499369*w[2]*fedge[26]-0.3354101966249685*dxv[2]*fedge[26]+1.190784930203603*w[2]*fskin[23]-0.5953924651018015*dxv[2]*fskin[23]+1.190784930203603*w[2]*fedge[23]-0.5953924651018015*dxv[2]*fedge[23]-0.9375*w[2]*fskin[20]+0.46875*dxv[2]*fskin[20]+0.9375*w[2]*fedge[20]-0.46875*dxv[2]*fedge[20]; 
  edgeSurf[21] = (-3.140625*w[2]*fskin[21])+1.5703125*dxv[2]*fskin[21]-0.140625*w[2]*fedge[21]+0.0703125*dxv[2]*fedge[21]+5.022775277112744*w[2]*fskin[13]-2.511387638556372*dxv[2]*fskin[13]+0.3025768239224544*w[2]*fedge[13]-0.1512884119612272*dxv[2]*fedge[13]-3.773364712030896*w[2]*fskin[7]+1.886682356015448*dxv[2]*fskin[7]+0.4192627457812106*w[2]*fedge[7]-0.2096313728906053*dxv[2]*fedge[7]; 
  edgeSurf[22] = (-3.140625*w[2]*fskin[22])+1.5703125*dxv[2]*fskin[22]-0.140625*w[2]*fedge[22]+0.0703125*dxv[2]*fedge[22]+5.022775277112744*w[2]*fskin[14]-2.511387638556372*dxv[2]*fskin[14]+0.3025768239224544*w[2]*fedge[14]-0.1512884119612272*dxv[2]*fedge[14]-3.773364712030896*w[2]*fskin[8]+1.886682356015448*dxv[2]*fskin[8]+0.4192627457812106*w[2]*fedge[8]-0.2096313728906053*dxv[2]*fedge[8]; 
  edgeSurf[23] = 1.585502557353661*w[2]*fskin[26]-0.7927512786768306*dxv[2]*fskin[26]-0.7382874503707888*w[2]*fedge[26]+0.3691437251853944*dxv[2]*fedge[26]-2.671875*w[2]*fskin[23]+1.3359375*dxv[2]*fskin[23]-1.453125*w[2]*fedge[23]+0.7265625*dxv[2]*fedge[23]+2.056810333988042*w[2]*fskin[20]-1.028405166994021*dxv[2]*fskin[20]-1.190784930203603*w[2]*fedge[20]+0.5953924651018015*dxv[2]*fedge[20]; 
  edgeSurf[24] = (-3.140625*w[2]*fskin[24])+1.5703125*dxv[2]*fskin[24]-0.140625*w[2]*fedge[24]+0.0703125*dxv[2]*fedge[24]+5.022775277112744*w[2]*fskin[17]-2.511387638556372*dxv[2]*fskin[17]+0.3025768239224545*w[2]*fedge[17]-0.1512884119612272*dxv[2]*fedge[17]-3.773364712030894*w[2]*fskin[11]+1.886682356015447*dxv[2]*fskin[11]+0.4192627457812105*w[2]*fedge[11]-0.2096313728906053*dxv[2]*fedge[11]; 
  edgeSurf[25] = (-3.140625*w[2]*fskin[25])+1.5703125*dxv[2]*fskin[25]-0.140625*w[2]*fedge[25]+0.0703125*dxv[2]*fedge[25]+5.022775277112744*w[2]*fskin[18]-2.511387638556372*dxv[2]*fskin[18]+0.3025768239224545*w[2]*fedge[18]-0.1512884119612272*dxv[2]*fedge[18]-3.773364712030894*w[2]*fskin[12]+1.886682356015447*dxv[2]*fskin[12]+0.4192627457812105*w[2]*fedge[12]-0.2096313728906053*dxv[2]*fedge[12]; 
  edgeSurf[26] = (-3.140625*w[2]*fskin[26])+1.5703125*dxv[2]*fskin[26]-0.140625*w[2]*fedge[26]+0.0703125*dxv[2]*fedge[26]+5.022775277112744*w[2]*fskin[23]-2.511387638556372*dxv[2]*fskin[23]+0.3025768239224545*w[2]*fedge[23]-0.1512884119612272*dxv[2]*fedge[23]-3.773364712030896*w[2]*fskin[20]+1.886682356015448*dxv[2]*fskin[20]+0.4192627457812106*w[2]*fedge[20]-0.2096313728906053*dxv[2]*fedge[20]; 

  double boundSurf[27] = {0.0}; 
  boundSurf[3] = (-1.936491673103709*w[2]*fskin[9])-0.9682458365518543*dxv[2]*fskin[9]-1.5*w[2]*fskin[3]-0.75*dxv[2]*fskin[3]-0.8660254037844386*fskin[0]*w[2]-0.4330127018922193*fskin[0]*dxv[2]; 
  boundSurf[5] = (-1.936491673103709*w[2]*fskin[15])-0.9682458365518543*dxv[2]*fskin[15]-1.5*w[2]*fskin[5]-0.75*dxv[2]*fskin[5]-0.8660254037844386*fskin[1]*w[2]-0.4330127018922193*fskin[1]*dxv[2]; 
  boundSurf[6] = (-1.936491673103709*w[2]*fskin[16])-0.9682458365518543*dxv[2]*fskin[16]-1.5*w[2]*fskin[6]-0.75*dxv[2]*fskin[6]-0.8660254037844386*fskin[2]*w[2]-0.4330127018922193*dxv[2]*fskin[2]; 
  boundSurf[9] = (-7.5*w[2]*fskin[9])-3.75*dxv[2]*fskin[9]-5.809475019311125*w[2]*fskin[3]-2.904737509655563*dxv[2]*fskin[3]-3.354101966249685*fskin[0]*w[2]-1.677050983124842*fskin[0]*dxv[2]; 
  boundSurf[10] = (-1.936491673103709*w[2]*fskin[19])-0.9682458365518543*dxv[2]*fskin[19]-1.5*w[2]*fskin[10]-0.75*dxv[2]*fskin[10]-0.8660254037844386*w[2]*fskin[4]-0.4330127018922193*dxv[2]*fskin[4]; 
  boundSurf[13] = (-1.936491673103709*w[2]*fskin[21])-0.9682458365518543*dxv[2]*fskin[21]-1.5*w[2]*fskin[13]-0.75*dxv[2]*fskin[13]-0.8660254037844387*w[2]*fskin[7]-0.4330127018922194*dxv[2]*fskin[7]; 
  boundSurf[14] = (-1.936491673103709*w[2]*fskin[22])-0.9682458365518543*dxv[2]*fskin[22]-1.5*w[2]*fskin[14]-0.75*dxv[2]*fskin[14]-0.8660254037844387*w[2]*fskin[8]-0.4330127018922194*dxv[2]*fskin[8]; 
  boundSurf[15] = (-7.5*w[2]*fskin[15])-3.75*dxv[2]*fskin[15]-5.809475019311126*w[2]*fskin[5]-2.904737509655563*dxv[2]*fskin[5]-3.354101966249684*fskin[1]*w[2]-1.677050983124842*fskin[1]*dxv[2]; 
  boundSurf[16] = (-7.5*w[2]*fskin[16])-3.75*dxv[2]*fskin[16]-5.809475019311126*w[2]*fskin[6]-2.904737509655563*dxv[2]*fskin[6]-3.354101966249684*fskin[2]*w[2]-1.677050983124842*dxv[2]*fskin[2]; 
  boundSurf[17] = (-1.936491673103709*w[2]*fskin[24])-0.9682458365518543*dxv[2]*fskin[24]-1.5*w[2]*fskin[17]-0.75*dxv[2]*fskin[17]-0.8660254037844387*w[2]*fskin[11]-0.4330127018922194*dxv[2]*fskin[11]; 
  boundSurf[18] = (-1.936491673103709*w[2]*fskin[25])-0.9682458365518543*dxv[2]*fskin[25]-1.5*w[2]*fskin[18]-0.75*dxv[2]*fskin[18]-0.8660254037844387*w[2]*fskin[12]-0.4330127018922194*dxv[2]*fskin[12]; 
  boundSurf[19] = (-7.5*w[2]*fskin[19])-3.75*dxv[2]*fskin[19]-5.809475019311125*w[2]*fskin[10]-2.904737509655563*dxv[2]*fskin[10]-3.354101966249685*w[2]*fskin[4]-1.677050983124842*dxv[2]*fskin[4]; 
  boundSurf[21] = (-7.5*w[2]*fskin[21])-3.75*dxv[2]*fskin[21]-5.809475019311126*w[2]*fskin[13]-2.904737509655563*dxv[2]*fskin[13]-3.354101966249685*w[2]*fskin[7]-1.677050983124842*dxv[2]*fskin[7]; 
  boundSurf[22] = (-7.5*w[2]*fskin[22])-3.75*dxv[2]*fskin[22]-5.809475019311126*w[2]*fskin[14]-2.904737509655563*dxv[2]*fskin[14]-3.354101966249685*w[2]*fskin[8]-1.677050983124842*dxv[2]*fskin[8]; 
  boundSurf[23] = (-1.936491673103709*w[2]*fskin[26])-0.9682458365518543*dxv[2]*fskin[26]-1.5*w[2]*fskin[23]-0.75*dxv[2]*fskin[23]-0.8660254037844386*w[2]*fskin[20]-0.4330127018922193*dxv[2]*fskin[20]; 
  boundSurf[24] = (-7.5*w[2]*fskin[24])-3.75*dxv[2]*fskin[24]-5.809475019311125*w[2]*fskin[17]-2.904737509655563*dxv[2]*fskin[17]-3.354101966249684*w[2]*fskin[11]-1.677050983124842*dxv[2]*fskin[11]; 
  boundSurf[25] = (-7.5*w[2]*fskin[25])-3.75*dxv[2]*fskin[25]-5.809475019311125*w[2]*fskin[18]-2.904737509655563*dxv[2]*fskin[18]-3.354101966249684*w[2]*fskin[12]-1.677050983124842*dxv[2]*fskin[12]; 
  boundSurf[26] = (-7.5*w[2]*fskin[26])-3.75*dxv[2]*fskin[26]-5.809475019311125*w[2]*fskin[23]-2.904737509655563*dxv[2]*fskin[23]-3.354101966249685*w[2]*fskin[20]-1.677050983124842*dxv[2]*fskin[20]; 

  edgeSurf_incr[0] = 0.7071067811865475*facDiff[2]*edgeSurf[7]+0.7071067811865475*edgeSurf[1]*facDiff[1]+0.7071067811865475*edgeSurf[0]*facDiff[0]; 
  edgeSurf_incr[1] = 0.6324555320336759*facDiff[1]*edgeSurf[7]+0.6324555320336759*edgeSurf[1]*facDiff[2]+0.7071067811865475*edgeSurf[0]*facDiff[1]+0.7071067811865475*facDiff[0]*edgeSurf[1]; 
  edgeSurf_incr[2] = 0.7071067811865475*facDiff[2]*edgeSurf[11]+0.7071067811865475*facDiff[1]*edgeSurf[4]+0.7071067811865475*facDiff[0]*edgeSurf[2]; 
  edgeSurf_incr[3] = 0.7071067811865475*facDiff[2]*edgeSurf[13]+0.7071067811865475*facDiff[1]*edgeSurf[5]+0.7071067811865475*facDiff[0]*edgeSurf[3]; 
  edgeSurf_incr[4] = 0.632455532033676*facDiff[1]*edgeSurf[11]+0.6324555320336759*facDiff[2]*edgeSurf[4]+0.7071067811865475*facDiff[0]*edgeSurf[4]+0.7071067811865475*facDiff[1]*edgeSurf[2]; 
  edgeSurf_incr[5] = 0.632455532033676*facDiff[1]*edgeSurf[13]+0.6324555320336759*facDiff[2]*edgeSurf[5]+0.7071067811865475*facDiff[0]*edgeSurf[5]+0.7071067811865475*facDiff[1]*edgeSurf[3]; 
  edgeSurf_incr[6] = 0.7071067811865475*facDiff[2]*edgeSurf[17]+0.7071067811865475*facDiff[1]*edgeSurf[10]+0.7071067811865475*facDiff[0]*edgeSurf[6]; 
  edgeSurf_incr[7] = 0.4517539514526256*facDiff[2]*edgeSurf[7]+0.7071067811865475*facDiff[0]*edgeSurf[7]+0.7071067811865475*edgeSurf[0]*facDiff[2]+0.6324555320336759*edgeSurf[1]*facDiff[1]; 
  edgeSurf_incr[8] = 0.7071067811865475*facDiff[2]*edgeSurf[20]+0.7071067811865475*facDiff[1]*edgeSurf[12]+0.7071067811865475*facDiff[0]*edgeSurf[8]; 
  edgeSurf_incr[9] = 0.7071067811865475*facDiff[2]*edgeSurf[21]+0.7071067811865475*facDiff[1]*edgeSurf[15]+0.7071067811865475*facDiff[0]*edgeSurf[9]; 
  edgeSurf_incr[10] = 0.6324555320336759*facDiff[1]*edgeSurf[17]+0.6324555320336759*facDiff[2]*edgeSurf[10]+0.7071067811865475*facDiff[0]*edgeSurf[10]+0.7071067811865475*facDiff[1]*edgeSurf[6]; 
  edgeSurf_incr[11] = 0.4517539514526256*facDiff[2]*edgeSurf[11]+0.7071067811865475*facDiff[0]*edgeSurf[11]+0.632455532033676*facDiff[1]*edgeSurf[4]+0.7071067811865475*edgeSurf[2]*facDiff[2]; 
  edgeSurf_incr[12] = 0.632455532033676*facDiff[1]*edgeSurf[20]+0.6324555320336759*facDiff[2]*edgeSurf[12]+0.7071067811865475*facDiff[0]*edgeSurf[12]+0.7071067811865475*facDiff[1]*edgeSurf[8]; 
  edgeSurf_incr[13] = 0.4517539514526256*facDiff[2]*edgeSurf[13]+0.7071067811865475*facDiff[0]*edgeSurf[13]+0.632455532033676*facDiff[1]*edgeSurf[5]+0.7071067811865475*facDiff[2]*edgeSurf[3]; 
  edgeSurf_incr[14] = 0.7071067811865475*facDiff[2]*edgeSurf[23]+0.7071067811865475*facDiff[1]*edgeSurf[18]+0.7071067811865475*facDiff[0]*edgeSurf[14]; 
  edgeSurf_incr[15] = 0.632455532033676*facDiff[1]*edgeSurf[21]+0.6324555320336759*facDiff[2]*edgeSurf[15]+0.7071067811865475*facDiff[0]*edgeSurf[15]+0.7071067811865475*facDiff[1]*edgeSurf[9]; 
  edgeSurf_incr[16] = 0.7071067811865475*facDiff[2]*edgeSurf[24]+0.7071067811865475*facDiff[1]*edgeSurf[19]+0.7071067811865475*facDiff[0]*edgeSurf[16]; 
  edgeSurf_incr[17] = 0.4517539514526256*facDiff[2]*edgeSurf[17]+0.7071067811865475*facDiff[0]*edgeSurf[17]+0.6324555320336759*facDiff[1]*edgeSurf[10]+0.7071067811865475*facDiff[2]*edgeSurf[6]; 
  edgeSurf_incr[18] = 0.6324555320336759*facDiff[1]*edgeSurf[23]+0.6324555320336759*facDiff[2]*edgeSurf[18]+0.7071067811865475*facDiff[0]*edgeSurf[18]+0.7071067811865475*facDiff[1]*edgeSurf[14]; 
  edgeSurf_incr[19] = 0.6324555320336759*facDiff[1]*edgeSurf[24]+0.6324555320336759*facDiff[2]*edgeSurf[19]+0.7071067811865475*facDiff[0]*edgeSurf[19]+0.7071067811865475*facDiff[1]*edgeSurf[16]; 
  edgeSurf_incr[20] = 0.4517539514526256*facDiff[2]*edgeSurf[20]+0.7071067811865475*facDiff[0]*edgeSurf[20]+0.632455532033676*facDiff[1]*edgeSurf[12]+0.7071067811865475*facDiff[2]*edgeSurf[8]; 
  edgeSurf_incr[21] = 0.4517539514526256*facDiff[2]*edgeSurf[21]+0.7071067811865475*facDiff[0]*edgeSurf[21]+0.632455532033676*facDiff[1]*edgeSurf[15]+0.7071067811865475*facDiff[2]*edgeSurf[9]; 
  edgeSurf_incr[22] = 0.7071067811865475*facDiff[2]*edgeSurf[26]+0.7071067811865475*facDiff[1]*edgeSurf[25]+0.7071067811865475*facDiff[0]*edgeSurf[22]; 
  edgeSurf_incr[23] = 0.4517539514526256*facDiff[2]*edgeSurf[23]+0.7071067811865475*facDiff[0]*edgeSurf[23]+0.6324555320336759*facDiff[1]*edgeSurf[18]+0.7071067811865475*facDiff[2]*edgeSurf[14]; 
  edgeSurf_incr[24] = 0.4517539514526256*facDiff[2]*edgeSurf[24]+0.7071067811865475*facDiff[0]*edgeSurf[24]+0.6324555320336759*facDiff[1]*edgeSurf[19]+0.7071067811865475*facDiff[2]*edgeSurf[16]; 
  edgeSurf_incr[25] = 0.6324555320336759*facDiff[1]*edgeSurf[26]+0.6324555320336759*facDiff[2]*edgeSurf[25]+0.7071067811865475*facDiff[0]*edgeSurf[25]+0.7071067811865475*facDiff[1]*edgeSurf[22]; 
  edgeSurf_incr[26] = 0.4517539514526256*facDiff[2]*edgeSurf[26]+0.7071067811865475*facDiff[0]*edgeSurf[26]+0.6324555320336759*facDiff[1]*edgeSurf[25]+0.7071067811865475*facDiff[2]*edgeSurf[22]; 

  boundSurf_incr[0] = 0.7071067811865475*facDiff[2]*boundSurf[7]+0.7071067811865475*boundSurf[1]*facDiff[1]+0.7071067811865475*boundSurf[0]*facDiff[0]; 
  boundSurf_incr[1] = 0.6324555320336759*facDiff[1]*boundSurf[7]+0.6324555320336759*boundSurf[1]*facDiff[2]+0.7071067811865475*boundSurf[0]*facDiff[1]+0.7071067811865475*facDiff[0]*boundSurf[1]; 
  boundSurf_incr[2] = 0.7071067811865475*facDiff[2]*boundSurf[11]+0.7071067811865475*facDiff[1]*boundSurf[4]+0.7071067811865475*facDiff[0]*boundSurf[2]; 
  boundSurf_incr[3] = 0.7071067811865475*facDiff[2]*boundSurf[13]+0.7071067811865475*facDiff[1]*boundSurf[5]+0.7071067811865475*facDiff[0]*boundSurf[3]; 
  boundSurf_incr[4] = 0.632455532033676*facDiff[1]*boundSurf[11]+0.6324555320336759*facDiff[2]*boundSurf[4]+0.7071067811865475*facDiff[0]*boundSurf[4]+0.7071067811865475*facDiff[1]*boundSurf[2]; 
  boundSurf_incr[5] = 0.632455532033676*facDiff[1]*boundSurf[13]+0.6324555320336759*facDiff[2]*boundSurf[5]+0.7071067811865475*facDiff[0]*boundSurf[5]+0.7071067811865475*facDiff[1]*boundSurf[3]; 
  boundSurf_incr[6] = 0.7071067811865475*facDiff[2]*boundSurf[17]+0.7071067811865475*facDiff[1]*boundSurf[10]+0.7071067811865475*facDiff[0]*boundSurf[6]; 
  boundSurf_incr[7] = 0.4517539514526256*facDiff[2]*boundSurf[7]+0.7071067811865475*facDiff[0]*boundSurf[7]+0.7071067811865475*boundSurf[0]*facDiff[2]+0.6324555320336759*boundSurf[1]*facDiff[1]; 
  boundSurf_incr[8] = 0.7071067811865475*facDiff[2]*boundSurf[20]+0.7071067811865475*facDiff[1]*boundSurf[12]+0.7071067811865475*facDiff[0]*boundSurf[8]; 
  boundSurf_incr[9] = 0.7071067811865475*facDiff[2]*boundSurf[21]+0.7071067811865475*facDiff[1]*boundSurf[15]+0.7071067811865475*facDiff[0]*boundSurf[9]; 
  boundSurf_incr[10] = 0.6324555320336759*facDiff[1]*boundSurf[17]+0.6324555320336759*facDiff[2]*boundSurf[10]+0.7071067811865475*facDiff[0]*boundSurf[10]+0.7071067811865475*facDiff[1]*boundSurf[6]; 
  boundSurf_incr[11] = 0.4517539514526256*facDiff[2]*boundSurf[11]+0.7071067811865475*facDiff[0]*boundSurf[11]+0.632455532033676*facDiff[1]*boundSurf[4]+0.7071067811865475*boundSurf[2]*facDiff[2]; 
  boundSurf_incr[12] = 0.632455532033676*facDiff[1]*boundSurf[20]+0.6324555320336759*facDiff[2]*boundSurf[12]+0.7071067811865475*facDiff[0]*boundSurf[12]+0.7071067811865475*facDiff[1]*boundSurf[8]; 
  boundSurf_incr[13] = 0.4517539514526256*facDiff[2]*boundSurf[13]+0.7071067811865475*facDiff[0]*boundSurf[13]+0.632455532033676*facDiff[1]*boundSurf[5]+0.7071067811865475*facDiff[2]*boundSurf[3]; 
  boundSurf_incr[14] = 0.7071067811865475*facDiff[2]*boundSurf[23]+0.7071067811865475*facDiff[1]*boundSurf[18]+0.7071067811865475*facDiff[0]*boundSurf[14]; 
  boundSurf_incr[15] = 0.632455532033676*facDiff[1]*boundSurf[21]+0.6324555320336759*facDiff[2]*boundSurf[15]+0.7071067811865475*facDiff[0]*boundSurf[15]+0.7071067811865475*facDiff[1]*boundSurf[9]; 
  boundSurf_incr[16] = 0.7071067811865475*facDiff[2]*boundSurf[24]+0.7071067811865475*facDiff[1]*boundSurf[19]+0.7071067811865475*facDiff[0]*boundSurf[16]; 
  boundSurf_incr[17] = 0.4517539514526256*facDiff[2]*boundSurf[17]+0.7071067811865475*facDiff[0]*boundSurf[17]+0.6324555320336759*facDiff[1]*boundSurf[10]+0.7071067811865475*facDiff[2]*boundSurf[6]; 
  boundSurf_incr[18] = 0.6324555320336759*facDiff[1]*boundSurf[23]+0.6324555320336759*facDiff[2]*boundSurf[18]+0.7071067811865475*facDiff[0]*boundSurf[18]+0.7071067811865475*facDiff[1]*boundSurf[14]; 
  boundSurf_incr[19] = 0.6324555320336759*facDiff[1]*boundSurf[24]+0.6324555320336759*facDiff[2]*boundSurf[19]+0.7071067811865475*facDiff[0]*boundSurf[19]+0.7071067811865475*facDiff[1]*boundSurf[16]; 
  boundSurf_incr[20] = 0.4517539514526256*facDiff[2]*boundSurf[20]+0.7071067811865475*facDiff[0]*boundSurf[20]+0.632455532033676*facDiff[1]*boundSurf[12]+0.7071067811865475*facDiff[2]*boundSurf[8]; 
  boundSurf_incr[21] = 0.4517539514526256*facDiff[2]*boundSurf[21]+0.7071067811865475*facDiff[0]*boundSurf[21]+0.632455532033676*facDiff[1]*boundSurf[15]+0.7071067811865475*facDiff[2]*boundSurf[9]; 
  boundSurf_incr[22] = 0.7071067811865475*facDiff[2]*boundSurf[26]+0.7071067811865475*facDiff[1]*boundSurf[25]+0.7071067811865475*facDiff[0]*boundSurf[22]; 
  boundSurf_incr[23] = 0.4517539514526256*facDiff[2]*boundSurf[23]+0.7071067811865475*facDiff[0]*boundSurf[23]+0.6324555320336759*facDiff[1]*boundSurf[18]+0.7071067811865475*facDiff[2]*boundSurf[14]; 
  boundSurf_incr[24] = 0.4517539514526256*facDiff[2]*boundSurf[24]+0.7071067811865475*facDiff[0]*boundSurf[24]+0.6324555320336759*facDiff[1]*boundSurf[19]+0.7071067811865475*facDiff[2]*boundSurf[16]; 
  boundSurf_incr[25] = 0.6324555320336759*facDiff[1]*boundSurf[26]+0.6324555320336759*facDiff[2]*boundSurf[25]+0.7071067811865475*facDiff[0]*boundSurf[25]+0.7071067811865475*facDiff[1]*boundSurf[22]; 
  boundSurf_incr[26] = 0.4517539514526256*facDiff[2]*boundSurf[26]+0.7071067811865475*facDiff[0]*boundSurf[26]+0.6324555320336759*facDiff[1]*boundSurf[25]+0.7071067811865475*facDiff[2]*boundSurf[22]; 

  } 

  out[0] += (vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*rdvSq4; 
  out[1] += (vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*rdvSq4; 
  out[2] += (vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*rdvSq4; 
  out[3] += (vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*rdvSq4; 
  out[4] += (vol_incr[4]+edgeSurf_incr[4]+boundSurf_incr[4])*rdvSq4; 
  out[5] += (vol_incr[5]+edgeSurf_incr[5]+boundSurf_incr[5])*rdvSq4; 
  out[6] += (vol_incr[6]+edgeSurf_incr[6]+boundSurf_incr[6])*rdvSq4; 
  out[7] += (vol_incr[7]+edgeSurf_incr[7]+boundSurf_incr[7])*rdvSq4; 
  out[8] += (vol_incr[8]+edgeSurf_incr[8]+boundSurf_incr[8])*rdvSq4; 
  out[9] += (vol_incr[9]+edgeSurf_incr[9]+boundSurf_incr[9])*rdvSq4; 
  out[10] += (vol_incr[10]+edgeSurf_incr[10]+boundSurf_incr[10])*rdvSq4; 
  out[11] += (vol_incr[11]+edgeSurf_incr[11]+boundSurf_incr[11])*rdvSq4; 
  out[12] += (vol_incr[12]+edgeSurf_incr[12]+boundSurf_incr[12])*rdvSq4; 
  out[13] += (vol_incr[13]+edgeSurf_incr[13]+boundSurf_incr[13])*rdvSq4; 
  out[14] += (vol_incr[14]+edgeSurf_incr[14]+boundSurf_incr[14])*rdvSq4; 
  out[15] += (vol_incr[15]+edgeSurf_incr[15]+boundSurf_incr[15])*rdvSq4; 
  out[16] += (vol_incr[16]+edgeSurf_incr[16]+boundSurf_incr[16])*rdvSq4; 
  out[17] += (vol_incr[17]+edgeSurf_incr[17]+boundSurf_incr[17])*rdvSq4; 
  out[18] += (vol_incr[18]+edgeSurf_incr[18]+boundSurf_incr[18])*rdvSq4; 
  out[19] += (vol_incr[19]+edgeSurf_incr[19]+boundSurf_incr[19])*rdvSq4; 
  out[20] += (vol_incr[20]+edgeSurf_incr[20]+boundSurf_incr[20])*rdvSq4; 
  out[21] += (vol_incr[21]+edgeSurf_incr[21]+boundSurf_incr[21])*rdvSq4; 
  out[22] += (vol_incr[22]+edgeSurf_incr[22]+boundSurf_incr[22])*rdvSq4; 
  out[23] += (vol_incr[23]+edgeSurf_incr[23]+boundSurf_incr[23])*rdvSq4; 
  out[24] += (vol_incr[24]+edgeSurf_incr[24]+boundSurf_incr[24])*rdvSq4; 
  out[25] += (vol_incr[25]+edgeSurf_incr[25]+boundSurf_incr[25])*rdvSq4; 
  out[26] += (vol_incr[26]+edgeSurf_incr[26]+boundSurf_incr[26])*rdvSq4; 
} 
