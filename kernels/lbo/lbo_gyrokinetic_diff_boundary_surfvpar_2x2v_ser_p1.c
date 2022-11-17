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
  vol_incr[16] = 3.354101966249685*nuVtSqSum[3]*fskin[5]+3.354101966249685*fskin[2]*nuVtSqSum[2]+3.354101966249685*fskin[1]*nuVtSqSum[1]+3.354101966249685*fskin[0]*nuVtSqSum[0]; 
  vol_incr[17] = 3.354101966249684*nuVtSqSum[2]*fskin[5]+3.354101966249684*fskin[2]*nuVtSqSum[3]+3.354101966249684*fskin[0]*nuVtSqSum[1]+3.354101966249684*nuVtSqSum[0]*fskin[1]; 
  vol_incr[18] = 3.354101966249684*nuVtSqSum[1]*fskin[5]+3.354101966249684*fskin[1]*nuVtSqSum[3]+3.354101966249684*fskin[0]*nuVtSqSum[2]+3.354101966249684*nuVtSqSum[0]*fskin[2]; 
  vol_incr[19] = 3.354101966249684*nuVtSqSum[3]*fskin[12]+3.354101966249684*nuVtSqSum[2]*fskin[9]+3.354101966249684*nuVtSqSum[1]*fskin[8]+3.354101966249684*nuVtSqSum[0]*fskin[4]; 
  vol_incr[20] = 3.354101966249685*nuVtSqSum[0]*fskin[5]+3.354101966249685*fskin[0]*nuVtSqSum[3]+3.354101966249685*fskin[1]*nuVtSqSum[2]+3.354101966249685*nuVtSqSum[1]*fskin[2]; 
  vol_incr[21] = 3.354101966249685*nuVtSqSum[2]*fskin[12]+3.354101966249685*nuVtSqSum[3]*fskin[9]+3.354101966249685*nuVtSqSum[0]*fskin[8]+3.354101966249685*nuVtSqSum[1]*fskin[4]; 
  vol_incr[22] = 3.354101966249685*nuVtSqSum[1]*fskin[12]+3.354101966249685*nuVtSqSum[0]*fskin[9]+3.354101966249685*nuVtSqSum[3]*fskin[8]+3.354101966249685*nuVtSqSum[2]*fskin[4]; 
  vol_incr[23] = 3.354101966249684*nuVtSqSum[0]*fskin[12]+3.354101966249684*nuVtSqSum[1]*fskin[9]+3.354101966249684*nuVtSqSum[2]*fskin[8]+3.354101966249684*nuVtSqSum[3]*fskin[4]; 

  double edgeSurf_incr[24] = {0.0}; 
  double boundSurf_incr[24] = {0.0}; 

  if (edge == -1) { 

  double edgeSurf[24] = {0.0}; 
  edgeSurf[0] = (-0.6708203932499369*fskin[16])+0.6708203932499369*fedge[16]-1.190784930203603*fskin[3]-1.190784930203603*fedge[3]-0.9375*fskin[0]+0.9375*fedge[0]; 
  edgeSurf[1] = (-0.6708203932499369*fskin[17])+0.6708203932499369*fedge[17]-1.190784930203603*fskin[6]-1.190784930203603*fedge[6]-0.9375*fskin[1]+0.9375*fedge[1]; 
  edgeSurf[2] = (-0.6708203932499369*fskin[18])+0.6708203932499369*fedge[18]-1.190784930203603*fskin[7]-1.190784930203603*fedge[7]-0.9375*fskin[2]+0.9375*fedge[2]; 
  edgeSurf[3] = (-1.585502557353661*fskin[16])+0.7382874503707888*fedge[16]-2.671875*fskin[3]-1.453125*fedge[3]-2.056810333988042*fskin[0]+1.190784930203603*fedge[0]; 
  edgeSurf[4] = (-0.6708203932499369*fskin[19])+0.6708203932499369*fedge[19]-1.190784930203603*fskin[10]-1.190784930203603*fedge[10]-0.9375*fskin[4]+0.9375*fedge[4]; 
  edgeSurf[5] = (-0.6708203932499369*fskin[20])+0.6708203932499369*fedge[20]-1.190784930203603*fskin[11]-1.190784930203603*fedge[11]-0.9375*fskin[5]+0.9375*fedge[5]; 
  edgeSurf[6] = (-1.585502557353661*fskin[17])+0.7382874503707888*fedge[17]-2.671875*fskin[6]-1.453125*fedge[6]-2.056810333988042*fskin[1]+1.190784930203603*fedge[1]; 
  edgeSurf[7] = (-1.585502557353661*fskin[18])+0.7382874503707888*fedge[18]-2.671875*fskin[7]-1.453125*fedge[7]-2.056810333988042*fskin[2]+1.190784930203603*fedge[2]; 
  edgeSurf[8] = (-0.6708203932499369*fskin[21])+0.6708203932499369*fedge[21]-1.190784930203603*fskin[13]-1.190784930203603*fedge[13]-0.9375*fskin[8]+0.9375*fedge[8]; 
  edgeSurf[9] = (-0.6708203932499369*fskin[22])+0.6708203932499369*fedge[22]-1.190784930203603*fskin[14]-1.190784930203603*fedge[14]-0.9375*fskin[9]+0.9375*fedge[9]; 
  edgeSurf[10] = (-1.585502557353661*fskin[19])+0.7382874503707888*fedge[19]-2.671875*fskin[10]-1.453125*fedge[10]-2.056810333988042*fskin[4]+1.190784930203603*fedge[4]; 
  edgeSurf[11] = (-1.585502557353661*fskin[20])+0.7382874503707888*fedge[20]-2.671875*fskin[11]-1.453125*fedge[11]-2.056810333988042*fskin[5]+1.190784930203603*fedge[5]; 
  edgeSurf[12] = (-0.6708203932499369*fskin[23])+0.6708203932499369*fedge[23]-1.190784930203603*fskin[15]-1.190784930203603*fedge[15]-0.9375*fskin[12]+0.9375*fedge[12]; 
  edgeSurf[13] = (-1.585502557353661*fskin[21])+0.7382874503707888*fedge[21]-2.671875*fskin[13]-1.453125*fedge[13]-2.056810333988042*fskin[8]+1.190784930203603*fedge[8]; 
  edgeSurf[14] = (-1.585502557353661*fskin[22])+0.7382874503707888*fedge[22]-2.671875*fskin[14]-1.453125*fedge[14]-2.056810333988042*fskin[9]+1.190784930203603*fedge[9]; 
  edgeSurf[15] = (-1.585502557353661*fskin[23])+0.7382874503707888*fedge[23]-2.671875*fskin[15]-1.453125*fedge[15]-2.056810333988042*fskin[12]+1.190784930203603*fedge[12]; 
  edgeSurf[16] = (-3.140625*fskin[16])-0.140625*fedge[16]-5.022775277112744*fskin[3]-0.3025768239224545*fedge[3]-3.773364712030896*fskin[0]+0.4192627457812106*fedge[0]; 
  edgeSurf[17] = (-3.140625*fskin[17])-0.140625*fedge[17]-5.022775277112744*fskin[6]-0.3025768239224544*fedge[6]-3.773364712030894*fskin[1]+0.4192627457812105*fedge[1]; 
  edgeSurf[18] = (-3.140625*fskin[18])-0.140625*fedge[18]-5.022775277112744*fskin[7]-0.3025768239224544*fedge[7]-3.773364712030894*fskin[2]+0.4192627457812105*fedge[2]; 
  edgeSurf[19] = (-3.140625*fskin[19])-0.140625*fedge[19]-5.022775277112744*fskin[10]-0.3025768239224544*fedge[10]-3.773364712030894*fskin[4]+0.4192627457812105*fedge[4]; 
  edgeSurf[20] = (-3.140625*fskin[20])-0.140625*fedge[20]-5.022775277112744*fskin[11]-0.3025768239224545*fedge[11]-3.773364712030896*fskin[5]+0.4192627457812106*fedge[5]; 
  edgeSurf[21] = (-3.140625*fskin[21])-0.140625*fedge[21]-5.022775277112744*fskin[13]-0.3025768239224545*fedge[13]-3.773364712030896*fskin[8]+0.4192627457812106*fedge[8]; 
  edgeSurf[22] = (-3.140625*fskin[22])-0.140625*fedge[22]-5.022775277112744*fskin[14]-0.3025768239224545*fedge[14]-3.773364712030896*fskin[9]+0.4192627457812106*fedge[9]; 
  edgeSurf[23] = (-3.140625*fskin[23])-0.140625*fedge[23]-5.022775277112744*fskin[15]-0.3025768239224544*fedge[15]-3.773364712030894*fskin[12]+0.4192627457812105*fedge[12]; 

  double boundSurf[24] = {0.0}; 
  boundSurf[3] = 1.936491673103709*fskin[16]-1.5*fskin[3]+0.8660254037844386*fskin[0]; 
  boundSurf[6] = 1.936491673103709*fskin[17]-1.5*fskin[6]+0.8660254037844386*fskin[1]; 
  boundSurf[7] = 1.936491673103709*fskin[18]-1.5*fskin[7]+0.8660254037844386*fskin[2]; 
  boundSurf[10] = 1.936491673103709*fskin[19]-1.5*fskin[10]+0.8660254037844386*fskin[4]; 
  boundSurf[11] = 1.936491673103709*fskin[20]-1.5*fskin[11]+0.8660254037844386*fskin[5]; 
  boundSurf[13] = 1.936491673103709*fskin[21]-1.5*fskin[13]+0.8660254037844386*fskin[8]; 
  boundSurf[14] = 1.936491673103709*fskin[22]-1.5*fskin[14]+0.8660254037844386*fskin[9]; 
  boundSurf[15] = 1.936491673103709*fskin[23]-1.5*fskin[15]+0.8660254037844386*fskin[12]; 
  boundSurf[16] = (-7.5*fskin[16])+5.809475019311125*fskin[3]-3.354101966249685*fskin[0]; 
  boundSurf[17] = (-7.5*fskin[17])+5.809475019311126*fskin[6]-3.354101966249684*fskin[1]; 
  boundSurf[18] = (-7.5*fskin[18])+5.809475019311126*fskin[7]-3.354101966249684*fskin[2]; 
  boundSurf[19] = (-7.5*fskin[19])+5.809475019311126*fskin[10]-3.354101966249684*fskin[4]; 
  boundSurf[20] = (-7.5*fskin[20])+5.809475019311125*fskin[11]-3.354101966249685*fskin[5]; 
  boundSurf[21] = (-7.5*fskin[21])+5.809475019311125*fskin[13]-3.354101966249685*fskin[8]; 
  boundSurf[22] = (-7.5*fskin[22])+5.809475019311125*fskin[14]-3.354101966249685*fskin[9]; 
  boundSurf[23] = (-7.5*fskin[23])+5.809475019311126*fskin[15]-3.354101966249684*fskin[12]; 

  edgeSurf_incr[0] = 0.5*nuVtSqSum[3]*edgeSurf[5]+0.5*edgeSurf[2]*nuVtSqSum[2]+0.5*edgeSurf[1]*nuVtSqSum[1]+0.5*edgeSurf[0]*nuVtSqSum[0]; 
  edgeSurf_incr[1] = 0.5*nuVtSqSum[2]*edgeSurf[5]+0.5*edgeSurf[2]*nuVtSqSum[3]+0.5*edgeSurf[0]*nuVtSqSum[1]+0.5*nuVtSqSum[0]*edgeSurf[1]; 
  edgeSurf_incr[2] = 0.5*nuVtSqSum[1]*edgeSurf[5]+0.5*edgeSurf[1]*nuVtSqSum[3]+0.5*edgeSurf[0]*nuVtSqSum[2]+0.5*nuVtSqSum[0]*edgeSurf[2]; 
  edgeSurf_incr[3] = 0.5*nuVtSqSum[3]*edgeSurf[11]+0.5*nuVtSqSum[2]*edgeSurf[7]+0.5*nuVtSqSum[1]*edgeSurf[6]+0.5*nuVtSqSum[0]*edgeSurf[3]; 
  edgeSurf_incr[4] = 0.5*nuVtSqSum[3]*edgeSurf[12]+0.5*nuVtSqSum[2]*edgeSurf[9]+0.5*nuVtSqSum[1]*edgeSurf[8]+0.5*nuVtSqSum[0]*edgeSurf[4]; 
  edgeSurf_incr[5] = 0.5*nuVtSqSum[0]*edgeSurf[5]+0.5*edgeSurf[0]*nuVtSqSum[3]+0.5*edgeSurf[1]*nuVtSqSum[2]+0.5*nuVtSqSum[1]*edgeSurf[2]; 
  edgeSurf_incr[6] = 0.5*nuVtSqSum[2]*edgeSurf[11]+0.5*nuVtSqSum[3]*edgeSurf[7]+0.5*nuVtSqSum[0]*edgeSurf[6]+0.5*nuVtSqSum[1]*edgeSurf[3]; 
  edgeSurf_incr[7] = 0.5*nuVtSqSum[1]*edgeSurf[11]+0.5*nuVtSqSum[0]*edgeSurf[7]+0.5*nuVtSqSum[3]*edgeSurf[6]+0.5*nuVtSqSum[2]*edgeSurf[3]; 
  edgeSurf_incr[8] = 0.5*nuVtSqSum[2]*edgeSurf[12]+0.5*nuVtSqSum[3]*edgeSurf[9]+0.5*nuVtSqSum[0]*edgeSurf[8]+0.5*nuVtSqSum[1]*edgeSurf[4]; 
  edgeSurf_incr[9] = 0.5*nuVtSqSum[1]*edgeSurf[12]+0.5*nuVtSqSum[0]*edgeSurf[9]+0.5*nuVtSqSum[3]*edgeSurf[8]+0.5*nuVtSqSum[2]*edgeSurf[4]; 
  edgeSurf_incr[10] = 0.5*nuVtSqSum[3]*edgeSurf[15]+0.5*nuVtSqSum[2]*edgeSurf[14]+0.5*nuVtSqSum[1]*edgeSurf[13]+0.5*nuVtSqSum[0]*edgeSurf[10]; 
  edgeSurf_incr[11] = 0.5*nuVtSqSum[0]*edgeSurf[11]+0.5*nuVtSqSum[1]*edgeSurf[7]+0.5*nuVtSqSum[2]*edgeSurf[6]+0.5*edgeSurf[3]*nuVtSqSum[3]; 
  edgeSurf_incr[12] = 0.5*nuVtSqSum[0]*edgeSurf[12]+0.5*nuVtSqSum[1]*edgeSurf[9]+0.5*nuVtSqSum[2]*edgeSurf[8]+0.5*nuVtSqSum[3]*edgeSurf[4]; 
  edgeSurf_incr[13] = 0.5*nuVtSqSum[2]*edgeSurf[15]+0.5*nuVtSqSum[3]*edgeSurf[14]+0.5*nuVtSqSum[0]*edgeSurf[13]+0.5*nuVtSqSum[1]*edgeSurf[10]; 
  edgeSurf_incr[14] = 0.5*nuVtSqSum[1]*edgeSurf[15]+0.5*nuVtSqSum[0]*edgeSurf[14]+0.5*nuVtSqSum[3]*edgeSurf[13]+0.5*nuVtSqSum[2]*edgeSurf[10]; 
  edgeSurf_incr[15] = 0.5*nuVtSqSum[0]*edgeSurf[15]+0.5*nuVtSqSum[1]*edgeSurf[14]+0.5*nuVtSqSum[2]*edgeSurf[13]+0.5*nuVtSqSum[3]*edgeSurf[10]; 
  edgeSurf_incr[16] = 0.5*nuVtSqSum[3]*edgeSurf[20]+0.5000000000000001*nuVtSqSum[2]*edgeSurf[18]+0.5000000000000001*nuVtSqSum[1]*edgeSurf[17]+0.5*nuVtSqSum[0]*edgeSurf[16]; 
  edgeSurf_incr[17] = 0.5000000000000001*nuVtSqSum[2]*edgeSurf[20]+0.5*nuVtSqSum[3]*edgeSurf[18]+0.5*nuVtSqSum[0]*edgeSurf[17]+0.5000000000000001*nuVtSqSum[1]*edgeSurf[16]; 
  edgeSurf_incr[18] = 0.5000000000000001*nuVtSqSum[1]*edgeSurf[20]+0.5*nuVtSqSum[0]*edgeSurf[18]+0.5*nuVtSqSum[3]*edgeSurf[17]+0.5000000000000001*nuVtSqSum[2]*edgeSurf[16]; 
  edgeSurf_incr[19] = 0.5*nuVtSqSum[3]*edgeSurf[23]+0.5000000000000001*nuVtSqSum[2]*edgeSurf[22]+0.5000000000000001*nuVtSqSum[1]*edgeSurf[21]+0.5*nuVtSqSum[0]*edgeSurf[19]; 
  edgeSurf_incr[20] = 0.5*nuVtSqSum[0]*edgeSurf[20]+0.5000000000000001*nuVtSqSum[1]*edgeSurf[18]+0.5000000000000001*nuVtSqSum[2]*edgeSurf[17]+0.5*nuVtSqSum[3]*edgeSurf[16]; 
  edgeSurf_incr[21] = 0.5000000000000001*nuVtSqSum[2]*edgeSurf[23]+0.5*nuVtSqSum[3]*edgeSurf[22]+0.5*nuVtSqSum[0]*edgeSurf[21]+0.5000000000000001*nuVtSqSum[1]*edgeSurf[19]; 
  edgeSurf_incr[22] = 0.5000000000000001*nuVtSqSum[1]*edgeSurf[23]+0.5*nuVtSqSum[0]*edgeSurf[22]+0.5*nuVtSqSum[3]*edgeSurf[21]+0.5000000000000001*nuVtSqSum[2]*edgeSurf[19]; 
  edgeSurf_incr[23] = 0.5*nuVtSqSum[0]*edgeSurf[23]+0.5000000000000001*nuVtSqSum[1]*edgeSurf[22]+0.5000000000000001*nuVtSqSum[2]*edgeSurf[21]+0.5*nuVtSqSum[3]*edgeSurf[19]; 

  boundSurf_incr[0] = 0.5*nuVtSqSum[3]*boundSurf[5]+0.5*boundSurf[2]*nuVtSqSum[2]+0.5*boundSurf[1]*nuVtSqSum[1]+0.5*boundSurf[0]*nuVtSqSum[0]; 
  boundSurf_incr[1] = 0.5*nuVtSqSum[2]*boundSurf[5]+0.5*boundSurf[2]*nuVtSqSum[3]+0.5*boundSurf[0]*nuVtSqSum[1]+0.5*nuVtSqSum[0]*boundSurf[1]; 
  boundSurf_incr[2] = 0.5*nuVtSqSum[1]*boundSurf[5]+0.5*boundSurf[1]*nuVtSqSum[3]+0.5*boundSurf[0]*nuVtSqSum[2]+0.5*nuVtSqSum[0]*boundSurf[2]; 
  boundSurf_incr[3] = 0.5*nuVtSqSum[3]*boundSurf[11]+0.5*nuVtSqSum[2]*boundSurf[7]+0.5*nuVtSqSum[1]*boundSurf[6]+0.5*nuVtSqSum[0]*boundSurf[3]; 
  boundSurf_incr[4] = 0.5*nuVtSqSum[3]*boundSurf[12]+0.5*nuVtSqSum[2]*boundSurf[9]+0.5*nuVtSqSum[1]*boundSurf[8]+0.5*nuVtSqSum[0]*boundSurf[4]; 
  boundSurf_incr[5] = 0.5*nuVtSqSum[0]*boundSurf[5]+0.5*boundSurf[0]*nuVtSqSum[3]+0.5*boundSurf[1]*nuVtSqSum[2]+0.5*nuVtSqSum[1]*boundSurf[2]; 
  boundSurf_incr[6] = 0.5*nuVtSqSum[2]*boundSurf[11]+0.5*nuVtSqSum[3]*boundSurf[7]+0.5*nuVtSqSum[0]*boundSurf[6]+0.5*nuVtSqSum[1]*boundSurf[3]; 
  boundSurf_incr[7] = 0.5*nuVtSqSum[1]*boundSurf[11]+0.5*nuVtSqSum[0]*boundSurf[7]+0.5*nuVtSqSum[3]*boundSurf[6]+0.5*nuVtSqSum[2]*boundSurf[3]; 
  boundSurf_incr[8] = 0.5*nuVtSqSum[2]*boundSurf[12]+0.5*nuVtSqSum[3]*boundSurf[9]+0.5*nuVtSqSum[0]*boundSurf[8]+0.5*nuVtSqSum[1]*boundSurf[4]; 
  boundSurf_incr[9] = 0.5*nuVtSqSum[1]*boundSurf[12]+0.5*nuVtSqSum[0]*boundSurf[9]+0.5*nuVtSqSum[3]*boundSurf[8]+0.5*nuVtSqSum[2]*boundSurf[4]; 
  boundSurf_incr[10] = 0.5*nuVtSqSum[3]*boundSurf[15]+0.5*nuVtSqSum[2]*boundSurf[14]+0.5*nuVtSqSum[1]*boundSurf[13]+0.5*nuVtSqSum[0]*boundSurf[10]; 
  boundSurf_incr[11] = 0.5*nuVtSqSum[0]*boundSurf[11]+0.5*nuVtSqSum[1]*boundSurf[7]+0.5*nuVtSqSum[2]*boundSurf[6]+0.5*boundSurf[3]*nuVtSqSum[3]; 
  boundSurf_incr[12] = 0.5*nuVtSqSum[0]*boundSurf[12]+0.5*nuVtSqSum[1]*boundSurf[9]+0.5*nuVtSqSum[2]*boundSurf[8]+0.5*nuVtSqSum[3]*boundSurf[4]; 
  boundSurf_incr[13] = 0.5*nuVtSqSum[2]*boundSurf[15]+0.5*nuVtSqSum[3]*boundSurf[14]+0.5*nuVtSqSum[0]*boundSurf[13]+0.5*nuVtSqSum[1]*boundSurf[10]; 
  boundSurf_incr[14] = 0.5*nuVtSqSum[1]*boundSurf[15]+0.5*nuVtSqSum[0]*boundSurf[14]+0.5*nuVtSqSum[3]*boundSurf[13]+0.5*nuVtSqSum[2]*boundSurf[10]; 
  boundSurf_incr[15] = 0.5*nuVtSqSum[0]*boundSurf[15]+0.5*nuVtSqSum[1]*boundSurf[14]+0.5*nuVtSqSum[2]*boundSurf[13]+0.5*nuVtSqSum[3]*boundSurf[10]; 
  boundSurf_incr[16] = 0.5*nuVtSqSum[3]*boundSurf[20]+0.5000000000000001*nuVtSqSum[2]*boundSurf[18]+0.5000000000000001*nuVtSqSum[1]*boundSurf[17]+0.5*nuVtSqSum[0]*boundSurf[16]; 
  boundSurf_incr[17] = 0.5000000000000001*nuVtSqSum[2]*boundSurf[20]+0.5*nuVtSqSum[3]*boundSurf[18]+0.5*nuVtSqSum[0]*boundSurf[17]+0.5000000000000001*nuVtSqSum[1]*boundSurf[16]; 
  boundSurf_incr[18] = 0.5000000000000001*nuVtSqSum[1]*boundSurf[20]+0.5*nuVtSqSum[0]*boundSurf[18]+0.5*nuVtSqSum[3]*boundSurf[17]+0.5000000000000001*nuVtSqSum[2]*boundSurf[16]; 
  boundSurf_incr[19] = 0.5*nuVtSqSum[3]*boundSurf[23]+0.5000000000000001*nuVtSqSum[2]*boundSurf[22]+0.5000000000000001*nuVtSqSum[1]*boundSurf[21]+0.5*nuVtSqSum[0]*boundSurf[19]; 
  boundSurf_incr[20] = 0.5*nuVtSqSum[0]*boundSurf[20]+0.5000000000000001*nuVtSqSum[1]*boundSurf[18]+0.5000000000000001*nuVtSqSum[2]*boundSurf[17]+0.5*nuVtSqSum[3]*boundSurf[16]; 
  boundSurf_incr[21] = 0.5000000000000001*nuVtSqSum[2]*boundSurf[23]+0.5*nuVtSqSum[3]*boundSurf[22]+0.5*nuVtSqSum[0]*boundSurf[21]+0.5000000000000001*nuVtSqSum[1]*boundSurf[19]; 
  boundSurf_incr[22] = 0.5000000000000001*nuVtSqSum[1]*boundSurf[23]+0.5*nuVtSqSum[0]*boundSurf[22]+0.5*nuVtSqSum[3]*boundSurf[21]+0.5000000000000001*nuVtSqSum[2]*boundSurf[19]; 
  boundSurf_incr[23] = 0.5*nuVtSqSum[0]*boundSurf[23]+0.5000000000000001*nuVtSqSum[1]*boundSurf[22]+0.5000000000000001*nuVtSqSum[2]*boundSurf[21]+0.5*nuVtSqSum[3]*boundSurf[19]; 


  } else { 

  double edgeSurf[24] = {0.0}; 
  edgeSurf[0] = (-0.6708203932499369*fskin[16])+0.6708203932499369*fedge[16]+1.190784930203603*fskin[3]+1.190784930203603*fedge[3]-0.9375*fskin[0]+0.9375*fedge[0]; 
  edgeSurf[1] = (-0.6708203932499369*fskin[17])+0.6708203932499369*fedge[17]+1.190784930203603*fskin[6]+1.190784930203603*fedge[6]-0.9375*fskin[1]+0.9375*fedge[1]; 
  edgeSurf[2] = (-0.6708203932499369*fskin[18])+0.6708203932499369*fedge[18]+1.190784930203603*fskin[7]+1.190784930203603*fedge[7]-0.9375*fskin[2]+0.9375*fedge[2]; 
  edgeSurf[3] = 1.585502557353661*fskin[16]-0.7382874503707888*fedge[16]-2.671875*fskin[3]-1.453125*fedge[3]+2.056810333988042*fskin[0]-1.190784930203603*fedge[0]; 
  edgeSurf[4] = (-0.6708203932499369*fskin[19])+0.6708203932499369*fedge[19]+1.190784930203603*fskin[10]+1.190784930203603*fedge[10]-0.9375*fskin[4]+0.9375*fedge[4]; 
  edgeSurf[5] = (-0.6708203932499369*fskin[20])+0.6708203932499369*fedge[20]+1.190784930203603*fskin[11]+1.190784930203603*fedge[11]-0.9375*fskin[5]+0.9375*fedge[5]; 
  edgeSurf[6] = 1.585502557353661*fskin[17]-0.7382874503707888*fedge[17]-2.671875*fskin[6]-1.453125*fedge[6]+2.056810333988042*fskin[1]-1.190784930203603*fedge[1]; 
  edgeSurf[7] = 1.585502557353661*fskin[18]-0.7382874503707888*fedge[18]-2.671875*fskin[7]-1.453125*fedge[7]+2.056810333988042*fskin[2]-1.190784930203603*fedge[2]; 
  edgeSurf[8] = (-0.6708203932499369*fskin[21])+0.6708203932499369*fedge[21]+1.190784930203603*fskin[13]+1.190784930203603*fedge[13]-0.9375*fskin[8]+0.9375*fedge[8]; 
  edgeSurf[9] = (-0.6708203932499369*fskin[22])+0.6708203932499369*fedge[22]+1.190784930203603*fskin[14]+1.190784930203603*fedge[14]-0.9375*fskin[9]+0.9375*fedge[9]; 
  edgeSurf[10] = 1.585502557353661*fskin[19]-0.7382874503707888*fedge[19]-2.671875*fskin[10]-1.453125*fedge[10]+2.056810333988042*fskin[4]-1.190784930203603*fedge[4]; 
  edgeSurf[11] = 1.585502557353661*fskin[20]-0.7382874503707888*fedge[20]-2.671875*fskin[11]-1.453125*fedge[11]+2.056810333988042*fskin[5]-1.190784930203603*fedge[5]; 
  edgeSurf[12] = (-0.6708203932499369*fskin[23])+0.6708203932499369*fedge[23]+1.190784930203603*fskin[15]+1.190784930203603*fedge[15]-0.9375*fskin[12]+0.9375*fedge[12]; 
  edgeSurf[13] = 1.585502557353661*fskin[21]-0.7382874503707888*fedge[21]-2.671875*fskin[13]-1.453125*fedge[13]+2.056810333988042*fskin[8]-1.190784930203603*fedge[8]; 
  edgeSurf[14] = 1.585502557353661*fskin[22]-0.7382874503707888*fedge[22]-2.671875*fskin[14]-1.453125*fedge[14]+2.056810333988042*fskin[9]-1.190784930203603*fedge[9]; 
  edgeSurf[15] = 1.585502557353661*fskin[23]-0.7382874503707888*fedge[23]-2.671875*fskin[15]-1.453125*fedge[15]+2.056810333988042*fskin[12]-1.190784930203603*fedge[12]; 
  edgeSurf[16] = (-3.140625*fskin[16])-0.140625*fedge[16]+5.022775277112744*fskin[3]+0.3025768239224545*fedge[3]-3.773364712030896*fskin[0]+0.4192627457812106*fedge[0]; 
  edgeSurf[17] = (-3.140625*fskin[17])-0.140625*fedge[17]+5.022775277112744*fskin[6]+0.3025768239224544*fedge[6]-3.773364712030894*fskin[1]+0.4192627457812105*fedge[1]; 
  edgeSurf[18] = (-3.140625*fskin[18])-0.140625*fedge[18]+5.022775277112744*fskin[7]+0.3025768239224544*fedge[7]-3.773364712030894*fskin[2]+0.4192627457812105*fedge[2]; 
  edgeSurf[19] = (-3.140625*fskin[19])-0.140625*fedge[19]+5.022775277112744*fskin[10]+0.3025768239224544*fedge[10]-3.773364712030894*fskin[4]+0.4192627457812105*fedge[4]; 
  edgeSurf[20] = (-3.140625*fskin[20])-0.140625*fedge[20]+5.022775277112744*fskin[11]+0.3025768239224545*fedge[11]-3.773364712030896*fskin[5]+0.4192627457812106*fedge[5]; 
  edgeSurf[21] = (-3.140625*fskin[21])-0.140625*fedge[21]+5.022775277112744*fskin[13]+0.3025768239224545*fedge[13]-3.773364712030896*fskin[8]+0.4192627457812106*fedge[8]; 
  edgeSurf[22] = (-3.140625*fskin[22])-0.140625*fedge[22]+5.022775277112744*fskin[14]+0.3025768239224545*fedge[14]-3.773364712030896*fskin[9]+0.4192627457812106*fedge[9]; 
  edgeSurf[23] = (-3.140625*fskin[23])-0.140625*fedge[23]+5.022775277112744*fskin[15]+0.3025768239224544*fedge[15]-3.773364712030894*fskin[12]+0.4192627457812105*fedge[12]; 

  double boundSurf[24] = {0.0}; 
  boundSurf[3] = (-1.936491673103709*fskin[16])-1.5*fskin[3]-0.8660254037844386*fskin[0]; 
  boundSurf[6] = (-1.936491673103709*fskin[17])-1.5*fskin[6]-0.8660254037844386*fskin[1]; 
  boundSurf[7] = (-1.936491673103709*fskin[18])-1.5*fskin[7]-0.8660254037844386*fskin[2]; 
  boundSurf[10] = (-1.936491673103709*fskin[19])-1.5*fskin[10]-0.8660254037844386*fskin[4]; 
  boundSurf[11] = (-1.936491673103709*fskin[20])-1.5*fskin[11]-0.8660254037844386*fskin[5]; 
  boundSurf[13] = (-1.936491673103709*fskin[21])-1.5*fskin[13]-0.8660254037844386*fskin[8]; 
  boundSurf[14] = (-1.936491673103709*fskin[22])-1.5*fskin[14]-0.8660254037844386*fskin[9]; 
  boundSurf[15] = (-1.936491673103709*fskin[23])-1.5*fskin[15]-0.8660254037844386*fskin[12]; 
  boundSurf[16] = (-7.5*fskin[16])-5.809475019311125*fskin[3]-3.354101966249685*fskin[0]; 
  boundSurf[17] = (-7.5*fskin[17])-5.809475019311126*fskin[6]-3.354101966249684*fskin[1]; 
  boundSurf[18] = (-7.5*fskin[18])-5.809475019311126*fskin[7]-3.354101966249684*fskin[2]; 
  boundSurf[19] = (-7.5*fskin[19])-5.809475019311126*fskin[10]-3.354101966249684*fskin[4]; 
  boundSurf[20] = (-7.5*fskin[20])-5.809475019311125*fskin[11]-3.354101966249685*fskin[5]; 
  boundSurf[21] = (-7.5*fskin[21])-5.809475019311125*fskin[13]-3.354101966249685*fskin[8]; 
  boundSurf[22] = (-7.5*fskin[22])-5.809475019311125*fskin[14]-3.354101966249685*fskin[9]; 
  boundSurf[23] = (-7.5*fskin[23])-5.809475019311126*fskin[15]-3.354101966249684*fskin[12]; 

  edgeSurf_incr[0] = 0.5*nuVtSqSum[3]*edgeSurf[5]+0.5*edgeSurf[2]*nuVtSqSum[2]+0.5*edgeSurf[1]*nuVtSqSum[1]+0.5*edgeSurf[0]*nuVtSqSum[0]; 
  edgeSurf_incr[1] = 0.5*nuVtSqSum[2]*edgeSurf[5]+0.5*edgeSurf[2]*nuVtSqSum[3]+0.5*edgeSurf[0]*nuVtSqSum[1]+0.5*nuVtSqSum[0]*edgeSurf[1]; 
  edgeSurf_incr[2] = 0.5*nuVtSqSum[1]*edgeSurf[5]+0.5*edgeSurf[1]*nuVtSqSum[3]+0.5*edgeSurf[0]*nuVtSqSum[2]+0.5*nuVtSqSum[0]*edgeSurf[2]; 
  edgeSurf_incr[3] = 0.5*nuVtSqSum[3]*edgeSurf[11]+0.5*nuVtSqSum[2]*edgeSurf[7]+0.5*nuVtSqSum[1]*edgeSurf[6]+0.5*nuVtSqSum[0]*edgeSurf[3]; 
  edgeSurf_incr[4] = 0.5*nuVtSqSum[3]*edgeSurf[12]+0.5*nuVtSqSum[2]*edgeSurf[9]+0.5*nuVtSqSum[1]*edgeSurf[8]+0.5*nuVtSqSum[0]*edgeSurf[4]; 
  edgeSurf_incr[5] = 0.5*nuVtSqSum[0]*edgeSurf[5]+0.5*edgeSurf[0]*nuVtSqSum[3]+0.5*edgeSurf[1]*nuVtSqSum[2]+0.5*nuVtSqSum[1]*edgeSurf[2]; 
  edgeSurf_incr[6] = 0.5*nuVtSqSum[2]*edgeSurf[11]+0.5*nuVtSqSum[3]*edgeSurf[7]+0.5*nuVtSqSum[0]*edgeSurf[6]+0.5*nuVtSqSum[1]*edgeSurf[3]; 
  edgeSurf_incr[7] = 0.5*nuVtSqSum[1]*edgeSurf[11]+0.5*nuVtSqSum[0]*edgeSurf[7]+0.5*nuVtSqSum[3]*edgeSurf[6]+0.5*nuVtSqSum[2]*edgeSurf[3]; 
  edgeSurf_incr[8] = 0.5*nuVtSqSum[2]*edgeSurf[12]+0.5*nuVtSqSum[3]*edgeSurf[9]+0.5*nuVtSqSum[0]*edgeSurf[8]+0.5*nuVtSqSum[1]*edgeSurf[4]; 
  edgeSurf_incr[9] = 0.5*nuVtSqSum[1]*edgeSurf[12]+0.5*nuVtSqSum[0]*edgeSurf[9]+0.5*nuVtSqSum[3]*edgeSurf[8]+0.5*nuVtSqSum[2]*edgeSurf[4]; 
  edgeSurf_incr[10] = 0.5*nuVtSqSum[3]*edgeSurf[15]+0.5*nuVtSqSum[2]*edgeSurf[14]+0.5*nuVtSqSum[1]*edgeSurf[13]+0.5*nuVtSqSum[0]*edgeSurf[10]; 
  edgeSurf_incr[11] = 0.5*nuVtSqSum[0]*edgeSurf[11]+0.5*nuVtSqSum[1]*edgeSurf[7]+0.5*nuVtSqSum[2]*edgeSurf[6]+0.5*edgeSurf[3]*nuVtSqSum[3]; 
  edgeSurf_incr[12] = 0.5*nuVtSqSum[0]*edgeSurf[12]+0.5*nuVtSqSum[1]*edgeSurf[9]+0.5*nuVtSqSum[2]*edgeSurf[8]+0.5*nuVtSqSum[3]*edgeSurf[4]; 
  edgeSurf_incr[13] = 0.5*nuVtSqSum[2]*edgeSurf[15]+0.5*nuVtSqSum[3]*edgeSurf[14]+0.5*nuVtSqSum[0]*edgeSurf[13]+0.5*nuVtSqSum[1]*edgeSurf[10]; 
  edgeSurf_incr[14] = 0.5*nuVtSqSum[1]*edgeSurf[15]+0.5*nuVtSqSum[0]*edgeSurf[14]+0.5*nuVtSqSum[3]*edgeSurf[13]+0.5*nuVtSqSum[2]*edgeSurf[10]; 
  edgeSurf_incr[15] = 0.5*nuVtSqSum[0]*edgeSurf[15]+0.5*nuVtSqSum[1]*edgeSurf[14]+0.5*nuVtSqSum[2]*edgeSurf[13]+0.5*nuVtSqSum[3]*edgeSurf[10]; 
  edgeSurf_incr[16] = 0.5*nuVtSqSum[3]*edgeSurf[20]+0.5000000000000001*nuVtSqSum[2]*edgeSurf[18]+0.5000000000000001*nuVtSqSum[1]*edgeSurf[17]+0.5*nuVtSqSum[0]*edgeSurf[16]; 
  edgeSurf_incr[17] = 0.5000000000000001*nuVtSqSum[2]*edgeSurf[20]+0.5*nuVtSqSum[3]*edgeSurf[18]+0.5*nuVtSqSum[0]*edgeSurf[17]+0.5000000000000001*nuVtSqSum[1]*edgeSurf[16]; 
  edgeSurf_incr[18] = 0.5000000000000001*nuVtSqSum[1]*edgeSurf[20]+0.5*nuVtSqSum[0]*edgeSurf[18]+0.5*nuVtSqSum[3]*edgeSurf[17]+0.5000000000000001*nuVtSqSum[2]*edgeSurf[16]; 
  edgeSurf_incr[19] = 0.5*nuVtSqSum[3]*edgeSurf[23]+0.5000000000000001*nuVtSqSum[2]*edgeSurf[22]+0.5000000000000001*nuVtSqSum[1]*edgeSurf[21]+0.5*nuVtSqSum[0]*edgeSurf[19]; 
  edgeSurf_incr[20] = 0.5*nuVtSqSum[0]*edgeSurf[20]+0.5000000000000001*nuVtSqSum[1]*edgeSurf[18]+0.5000000000000001*nuVtSqSum[2]*edgeSurf[17]+0.5*nuVtSqSum[3]*edgeSurf[16]; 
  edgeSurf_incr[21] = 0.5000000000000001*nuVtSqSum[2]*edgeSurf[23]+0.5*nuVtSqSum[3]*edgeSurf[22]+0.5*nuVtSqSum[0]*edgeSurf[21]+0.5000000000000001*nuVtSqSum[1]*edgeSurf[19]; 
  edgeSurf_incr[22] = 0.5000000000000001*nuVtSqSum[1]*edgeSurf[23]+0.5*nuVtSqSum[0]*edgeSurf[22]+0.5*nuVtSqSum[3]*edgeSurf[21]+0.5000000000000001*nuVtSqSum[2]*edgeSurf[19]; 
  edgeSurf_incr[23] = 0.5*nuVtSqSum[0]*edgeSurf[23]+0.5000000000000001*nuVtSqSum[1]*edgeSurf[22]+0.5000000000000001*nuVtSqSum[2]*edgeSurf[21]+0.5*nuVtSqSum[3]*edgeSurf[19]; 

  boundSurf_incr[0] = 0.5*nuVtSqSum[3]*boundSurf[5]+0.5*boundSurf[2]*nuVtSqSum[2]+0.5*boundSurf[1]*nuVtSqSum[1]+0.5*boundSurf[0]*nuVtSqSum[0]; 
  boundSurf_incr[1] = 0.5*nuVtSqSum[2]*boundSurf[5]+0.5*boundSurf[2]*nuVtSqSum[3]+0.5*boundSurf[0]*nuVtSqSum[1]+0.5*nuVtSqSum[0]*boundSurf[1]; 
  boundSurf_incr[2] = 0.5*nuVtSqSum[1]*boundSurf[5]+0.5*boundSurf[1]*nuVtSqSum[3]+0.5*boundSurf[0]*nuVtSqSum[2]+0.5*nuVtSqSum[0]*boundSurf[2]; 
  boundSurf_incr[3] = 0.5*nuVtSqSum[3]*boundSurf[11]+0.5*nuVtSqSum[2]*boundSurf[7]+0.5*nuVtSqSum[1]*boundSurf[6]+0.5*nuVtSqSum[0]*boundSurf[3]; 
  boundSurf_incr[4] = 0.5*nuVtSqSum[3]*boundSurf[12]+0.5*nuVtSqSum[2]*boundSurf[9]+0.5*nuVtSqSum[1]*boundSurf[8]+0.5*nuVtSqSum[0]*boundSurf[4]; 
  boundSurf_incr[5] = 0.5*nuVtSqSum[0]*boundSurf[5]+0.5*boundSurf[0]*nuVtSqSum[3]+0.5*boundSurf[1]*nuVtSqSum[2]+0.5*nuVtSqSum[1]*boundSurf[2]; 
  boundSurf_incr[6] = 0.5*nuVtSqSum[2]*boundSurf[11]+0.5*nuVtSqSum[3]*boundSurf[7]+0.5*nuVtSqSum[0]*boundSurf[6]+0.5*nuVtSqSum[1]*boundSurf[3]; 
  boundSurf_incr[7] = 0.5*nuVtSqSum[1]*boundSurf[11]+0.5*nuVtSqSum[0]*boundSurf[7]+0.5*nuVtSqSum[3]*boundSurf[6]+0.5*nuVtSqSum[2]*boundSurf[3]; 
  boundSurf_incr[8] = 0.5*nuVtSqSum[2]*boundSurf[12]+0.5*nuVtSqSum[3]*boundSurf[9]+0.5*nuVtSqSum[0]*boundSurf[8]+0.5*nuVtSqSum[1]*boundSurf[4]; 
  boundSurf_incr[9] = 0.5*nuVtSqSum[1]*boundSurf[12]+0.5*nuVtSqSum[0]*boundSurf[9]+0.5*nuVtSqSum[3]*boundSurf[8]+0.5*nuVtSqSum[2]*boundSurf[4]; 
  boundSurf_incr[10] = 0.5*nuVtSqSum[3]*boundSurf[15]+0.5*nuVtSqSum[2]*boundSurf[14]+0.5*nuVtSqSum[1]*boundSurf[13]+0.5*nuVtSqSum[0]*boundSurf[10]; 
  boundSurf_incr[11] = 0.5*nuVtSqSum[0]*boundSurf[11]+0.5*nuVtSqSum[1]*boundSurf[7]+0.5*nuVtSqSum[2]*boundSurf[6]+0.5*boundSurf[3]*nuVtSqSum[3]; 
  boundSurf_incr[12] = 0.5*nuVtSqSum[0]*boundSurf[12]+0.5*nuVtSqSum[1]*boundSurf[9]+0.5*nuVtSqSum[2]*boundSurf[8]+0.5*nuVtSqSum[3]*boundSurf[4]; 
  boundSurf_incr[13] = 0.5*nuVtSqSum[2]*boundSurf[15]+0.5*nuVtSqSum[3]*boundSurf[14]+0.5*nuVtSqSum[0]*boundSurf[13]+0.5*nuVtSqSum[1]*boundSurf[10]; 
  boundSurf_incr[14] = 0.5*nuVtSqSum[1]*boundSurf[15]+0.5*nuVtSqSum[0]*boundSurf[14]+0.5*nuVtSqSum[3]*boundSurf[13]+0.5*nuVtSqSum[2]*boundSurf[10]; 
  boundSurf_incr[15] = 0.5*nuVtSqSum[0]*boundSurf[15]+0.5*nuVtSqSum[1]*boundSurf[14]+0.5*nuVtSqSum[2]*boundSurf[13]+0.5*nuVtSqSum[3]*boundSurf[10]; 
  boundSurf_incr[16] = 0.5*nuVtSqSum[3]*boundSurf[20]+0.5000000000000001*nuVtSqSum[2]*boundSurf[18]+0.5000000000000001*nuVtSqSum[1]*boundSurf[17]+0.5*nuVtSqSum[0]*boundSurf[16]; 
  boundSurf_incr[17] = 0.5000000000000001*nuVtSqSum[2]*boundSurf[20]+0.5*nuVtSqSum[3]*boundSurf[18]+0.5*nuVtSqSum[0]*boundSurf[17]+0.5000000000000001*nuVtSqSum[1]*boundSurf[16]; 
  boundSurf_incr[18] = 0.5000000000000001*nuVtSqSum[1]*boundSurf[20]+0.5*nuVtSqSum[0]*boundSurf[18]+0.5*nuVtSqSum[3]*boundSurf[17]+0.5000000000000001*nuVtSqSum[2]*boundSurf[16]; 
  boundSurf_incr[19] = 0.5*nuVtSqSum[3]*boundSurf[23]+0.5000000000000001*nuVtSqSum[2]*boundSurf[22]+0.5000000000000001*nuVtSqSum[1]*boundSurf[21]+0.5*nuVtSqSum[0]*boundSurf[19]; 
  boundSurf_incr[20] = 0.5*nuVtSqSum[0]*boundSurf[20]+0.5000000000000001*nuVtSqSum[1]*boundSurf[18]+0.5000000000000001*nuVtSqSum[2]*boundSurf[17]+0.5*nuVtSqSum[3]*boundSurf[16]; 
  boundSurf_incr[21] = 0.5000000000000001*nuVtSqSum[2]*boundSurf[23]+0.5*nuVtSqSum[3]*boundSurf[22]+0.5*nuVtSqSum[0]*boundSurf[21]+0.5000000000000001*nuVtSqSum[1]*boundSurf[19]; 
  boundSurf_incr[22] = 0.5000000000000001*nuVtSqSum[1]*boundSurf[23]+0.5*nuVtSqSum[0]*boundSurf[22]+0.5*nuVtSqSum[3]*boundSurf[21]+0.5000000000000001*nuVtSqSum[2]*boundSurf[19]; 
  boundSurf_incr[23] = 0.5*nuVtSqSum[0]*boundSurf[23]+0.5000000000000001*nuVtSqSum[1]*boundSurf[22]+0.5000000000000001*nuVtSqSum[2]*boundSurf[21]+0.5*nuVtSqSum[3]*boundSurf[19]; 

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
} 
