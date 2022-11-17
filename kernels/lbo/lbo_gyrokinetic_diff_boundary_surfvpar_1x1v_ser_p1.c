#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH void lbo_gyrokinetic_diff_boundary_surfvpar_1x1v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  // w[2]: Cell-center coordinates. 
  // dxv[2]: Cell spacing. 
  // m_: species mass.
  // bmag_inv: 1/(magnetic field magnitude). 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[4]: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // fskin/edge: Distribution function in cells 
  // out: Incremented distribution function in cell 

  const double *nuVtSqSum = &nuPrimMomsSum[2];

  double rdvSq4 = 4.0/(dxv[1]*dxv[1]); 

  double vol_incr[6] = {0.0};
  vol_incr[4] = 4.743416490252569*fskin[1]*nuVtSqSum[1]+4.743416490252569*fskin[0]*nuVtSqSum[0]; 
  vol_incr[5] = 4.743416490252569*fskin[0]*nuVtSqSum[1]+4.743416490252569*nuVtSqSum[0]*fskin[1]; 

  double edgeSurf_incr[6] = {0.0}; 
  double boundSurf_incr[6] = {0.0}; 

  if (edge == -1) { 

  double edgeSurf[6] = {0.0}; 
  edgeSurf[0] = (-0.6708203932499369*fskin[4])+0.6708203932499369*fedge[4]-1.190784930203603*fskin[2]-1.190784930203603*fedge[2]-0.9375*fskin[0]+0.9375*fedge[0]; 
  edgeSurf[1] = (-0.6708203932499369*fskin[5])+0.6708203932499369*fedge[5]-1.190784930203603*fskin[3]-1.190784930203603*fedge[3]-0.9375*fskin[1]+0.9375*fedge[1]; 
  edgeSurf[2] = (-1.585502557353661*fskin[4])+0.7382874503707888*fedge[4]-2.671875*fskin[2]-1.453125*fedge[2]-2.056810333988042*fskin[0]+1.190784930203603*fedge[0]; 
  edgeSurf[3] = (-1.585502557353661*fskin[5])+0.7382874503707888*fedge[5]-2.671875*fskin[3]-1.453125*fedge[3]-2.056810333988042*fskin[1]+1.190784930203603*fedge[1]; 
  edgeSurf[4] = (-3.140625*fskin[4])-0.140625*fedge[4]-5.022775277112744*fskin[2]-0.3025768239224545*fedge[2]-3.773364712030896*fskin[0]+0.4192627457812106*fedge[0]; 
  edgeSurf[5] = (-3.140625*fskin[5])-0.140625*fedge[5]-5.022775277112744*fskin[3]-0.3025768239224544*fedge[3]-3.773364712030894*fskin[1]+0.4192627457812105*fedge[1]; 

  double boundSurf[6] = {0.0}; 
  boundSurf[2] = 1.936491673103709*fskin[4]-1.5*fskin[2]+0.8660254037844386*fskin[0]; 
  boundSurf[3] = 1.936491673103709*fskin[5]-1.5*fskin[3]+0.8660254037844386*fskin[1]; 
  boundSurf[4] = (-7.5*fskin[4])+5.809475019311125*fskin[2]-3.354101966249685*fskin[0]; 
  boundSurf[5] = (-7.5*fskin[5])+5.809475019311126*fskin[3]-3.354101966249684*fskin[1]; 

  edgeSurf_incr[0] = 0.7071067811865475*edgeSurf[1]*nuVtSqSum[1]+0.7071067811865475*edgeSurf[0]*nuVtSqSum[0]; 
  edgeSurf_incr[1] = 0.7071067811865475*edgeSurf[0]*nuVtSqSum[1]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[1]; 
  edgeSurf_incr[2] = 0.7071067811865475*nuVtSqSum[1]*edgeSurf[3]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[2]; 
  edgeSurf_incr[3] = 0.7071067811865475*nuVtSqSum[0]*edgeSurf[3]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[2]; 
  edgeSurf_incr[4] = 0.7071067811865475*nuVtSqSum[1]*edgeSurf[5]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[4]; 
  edgeSurf_incr[5] = 0.7071067811865475*nuVtSqSum[0]*edgeSurf[5]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[4]; 

  boundSurf_incr[0] = 0.7071067811865475*boundSurf[1]*nuVtSqSum[1]+0.7071067811865475*boundSurf[0]*nuVtSqSum[0]; 
  boundSurf_incr[1] = 0.7071067811865475*boundSurf[0]*nuVtSqSum[1]+0.7071067811865475*nuVtSqSum[0]*boundSurf[1]; 
  boundSurf_incr[2] = 0.7071067811865475*nuVtSqSum[1]*boundSurf[3]+0.7071067811865475*nuVtSqSum[0]*boundSurf[2]; 
  boundSurf_incr[3] = 0.7071067811865475*nuVtSqSum[0]*boundSurf[3]+0.7071067811865475*nuVtSqSum[1]*boundSurf[2]; 
  boundSurf_incr[4] = 0.7071067811865475*nuVtSqSum[1]*boundSurf[5]+0.7071067811865475*nuVtSqSum[0]*boundSurf[4]; 
  boundSurf_incr[5] = 0.7071067811865475*nuVtSqSum[0]*boundSurf[5]+0.7071067811865475*nuVtSqSum[1]*boundSurf[4]; 


  } else { 

  double edgeSurf[6] = {0.0}; 
  edgeSurf[0] = (-0.6708203932499369*fskin[4])+0.6708203932499369*fedge[4]+1.190784930203603*fskin[2]+1.190784930203603*fedge[2]-0.9375*fskin[0]+0.9375*fedge[0]; 
  edgeSurf[1] = (-0.6708203932499369*fskin[5])+0.6708203932499369*fedge[5]+1.190784930203603*fskin[3]+1.190784930203603*fedge[3]-0.9375*fskin[1]+0.9375*fedge[1]; 
  edgeSurf[2] = 1.585502557353661*fskin[4]-0.7382874503707888*fedge[4]-2.671875*fskin[2]-1.453125*fedge[2]+2.056810333988042*fskin[0]-1.190784930203603*fedge[0]; 
  edgeSurf[3] = 1.585502557353661*fskin[5]-0.7382874503707888*fedge[5]-2.671875*fskin[3]-1.453125*fedge[3]+2.056810333988042*fskin[1]-1.190784930203603*fedge[1]; 
  edgeSurf[4] = (-3.140625*fskin[4])-0.140625*fedge[4]+5.022775277112744*fskin[2]+0.3025768239224545*fedge[2]-3.773364712030896*fskin[0]+0.4192627457812106*fedge[0]; 
  edgeSurf[5] = (-3.140625*fskin[5])-0.140625*fedge[5]+5.022775277112744*fskin[3]+0.3025768239224544*fedge[3]-3.773364712030894*fskin[1]+0.4192627457812105*fedge[1]; 

  double boundSurf[6] = {0.0}; 
  boundSurf[2] = (-1.936491673103709*fskin[4])-1.5*fskin[2]-0.8660254037844386*fskin[0]; 
  boundSurf[3] = (-1.936491673103709*fskin[5])-1.5*fskin[3]-0.8660254037844386*fskin[1]; 
  boundSurf[4] = (-7.5*fskin[4])-5.809475019311125*fskin[2]-3.354101966249685*fskin[0]; 
  boundSurf[5] = (-7.5*fskin[5])-5.809475019311126*fskin[3]-3.354101966249684*fskin[1]; 

  edgeSurf_incr[0] = 0.7071067811865475*edgeSurf[1]*nuVtSqSum[1]+0.7071067811865475*edgeSurf[0]*nuVtSqSum[0]; 
  edgeSurf_incr[1] = 0.7071067811865475*edgeSurf[0]*nuVtSqSum[1]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[1]; 
  edgeSurf_incr[2] = 0.7071067811865475*nuVtSqSum[1]*edgeSurf[3]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[2]; 
  edgeSurf_incr[3] = 0.7071067811865475*nuVtSqSum[0]*edgeSurf[3]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[2]; 
  edgeSurf_incr[4] = 0.7071067811865475*nuVtSqSum[1]*edgeSurf[5]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[4]; 
  edgeSurf_incr[5] = 0.7071067811865475*nuVtSqSum[0]*edgeSurf[5]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[4]; 

  boundSurf_incr[0] = 0.7071067811865475*boundSurf[1]*nuVtSqSum[1]+0.7071067811865475*boundSurf[0]*nuVtSqSum[0]; 
  boundSurf_incr[1] = 0.7071067811865475*boundSurf[0]*nuVtSqSum[1]+0.7071067811865475*nuVtSqSum[0]*boundSurf[1]; 
  boundSurf_incr[2] = 0.7071067811865475*nuVtSqSum[1]*boundSurf[3]+0.7071067811865475*nuVtSqSum[0]*boundSurf[2]; 
  boundSurf_incr[3] = 0.7071067811865475*nuVtSqSum[0]*boundSurf[3]+0.7071067811865475*nuVtSqSum[1]*boundSurf[2]; 
  boundSurf_incr[4] = 0.7071067811865475*nuVtSqSum[1]*boundSurf[5]+0.7071067811865475*nuVtSqSum[0]*boundSurf[4]; 
  boundSurf_incr[5] = 0.7071067811865475*nuVtSqSum[0]*boundSurf[5]+0.7071067811865475*nuVtSqSum[1]*boundSurf[4]; 

  } 

  out[0] += (vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*rdvSq4; 
  out[1] += (vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*rdvSq4; 
  out[2] += (vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*rdvSq4; 
  out[3] += (vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*rdvSq4; 
  out[4] += (vol_incr[4]+edgeSurf_incr[4]+boundSurf_incr[4])*rdvSq4; 
  out[5] += (vol_incr[5]+edgeSurf_incr[5]+boundSurf_incr[5])*rdvSq4; 
} 
