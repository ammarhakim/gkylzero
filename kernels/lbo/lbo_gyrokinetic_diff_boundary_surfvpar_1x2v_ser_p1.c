#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH double lbo_gyrokinetic_diff_boundary_surfvpar_1x2v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  // w[3]: Cell-center coordinates. 
  // dxv[3]: Cell spacing. 
  // m_: species mass.
  // bmag_inv: 1/(magnetic field magnitude). 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[4]: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // fskin/edge: Distribution function in cells 
  // out: Incremented distribution function in cell 

  const double *nuVtSqSum = &nuPrimMomsSum[2];

  double rdvSq4 = 4.0/(dxv[1]*dxv[1]); 

  double bFacFskin[12] = {0.0}; 
  bFacFskin[8] = 6.708203932499369*fskin[0]; 
  bFacFskin[9] = 6.708203932499369*fskin[1]; 
  bFacFskin[10] = 6.708203932499369*fskin[3]; 
  bFacFskin[11] = 6.708203932499369*fskin[5]; 

  double vol_incr[12] = {0.0};
  vol_incr[8] = 0.105409255338946*(6.708203932499369*nuVtSqSum[1]*bFacFskin[9]+6.708203932499369*nuVtSqSum[0]*bFacFskin[8]); 
  vol_incr[9] = 0.7071067811865474*(nuVtSqSum[0]*bFacFskin[9]+nuVtSqSum[1]*bFacFskin[8]); 
  vol_incr[10] = 0.7071067811865474*(nuVtSqSum[1]*bFacFskin[11]+nuVtSqSum[0]*bFacFskin[10]); 
  vol_incr[11] = 0.105409255338946*(6.708203932499369*nuVtSqSum[0]*bFacFskin[11]+6.708203932499369*nuVtSqSum[1]*bFacFskin[10]); 

  double edgeSurf_incr[12] = {0.0}; 
  double boundSurf_incr[12] = {0.0}; 

  if (edge == -1) { 

  double edgeSurf[12] = {0.0}; 
  edgeSurf[0] = (-0.6708203932499369*fskin[8])+0.6708203932499369*fedge[8]-1.190784930203603*fskin[2]-1.190784930203603*fedge[2]-0.9375*fskin[0]+0.9375*fedge[0]; 
  edgeSurf[1] = (-0.6708203932499369*fskin[9])+0.6708203932499369*fedge[9]-1.190784930203603*fskin[4]-1.190784930203603*fedge[4]-0.9375*fskin[1]+0.9375*fedge[1]; 
  edgeSurf[2] = (-1.585502557353661*fskin[8])+0.7382874503707888*fedge[8]-2.671875*fskin[2]-1.453125*fedge[2]-2.056810333988042*fskin[0]+1.190784930203603*fedge[0]; 
  edgeSurf[3] = (-0.6708203932499369*fskin[10])+0.6708203932499369*fedge[10]-1.190784930203603*fskin[6]-1.190784930203603*fedge[6]-0.9375*fskin[3]+0.9375*fedge[3]; 
  edgeSurf[4] = (-1.585502557353661*fskin[9])+0.7382874503707888*fedge[9]-2.671875*fskin[4]-1.453125*fedge[4]-2.056810333988042*fskin[1]+1.190784930203603*fedge[1]; 
  edgeSurf[5] = (-0.6708203932499369*fskin[11])+0.6708203932499369*fedge[11]-1.190784930203603*fskin[7]-1.190784930203603*fedge[7]-0.9375*fskin[5]+0.9375*fedge[5]; 
  edgeSurf[6] = (-1.585502557353661*fskin[10])+0.7382874503707888*fedge[10]-2.671875*fskin[6]-1.453125*fedge[6]-2.056810333988042*fskin[3]+1.190784930203603*fedge[3]; 
  edgeSurf[7] = (-1.585502557353661*fskin[11])+0.7382874503707888*fedge[11]-2.671875*fskin[7]-1.453125*fedge[7]-2.056810333988042*fskin[5]+1.190784930203603*fedge[5]; 
  edgeSurf[8] = (-3.140625*fskin[8])-0.140625*fedge[8]-5.022775277112744*fskin[2]-0.3025768239224545*fedge[2]-3.773364712030896*fskin[0]+0.4192627457812106*fedge[0]; 
  edgeSurf[9] = (-3.140625*fskin[9])-0.140625*fedge[9]-5.022775277112744*fskin[4]-0.3025768239224544*fedge[4]-3.773364712030894*fskin[1]+0.4192627457812105*fedge[1]; 
  edgeSurf[10] = (-3.140625*fskin[10])-0.140625*fedge[10]-5.022775277112744*fskin[6]-0.3025768239224544*fedge[6]-3.773364712030894*fskin[3]+0.4192627457812105*fedge[3]; 
  edgeSurf[11] = (-3.140625*fskin[11])-0.140625*fedge[11]-5.022775277112744*fskin[7]-0.3025768239224545*fedge[7]-3.773364712030896*fskin[5]+0.4192627457812106*fedge[5]; 

  double boundSurf[12] = {0.0}; 
  boundSurf[2] = 1.936491673103709*fskin[8]-1.5*fskin[2]+0.8660254037844386*fskin[0]; 
  boundSurf[4] = 1.936491673103709*fskin[9]-1.5*fskin[4]+0.8660254037844386*fskin[1]; 
  boundSurf[6] = 1.936491673103709*fskin[10]-1.5*fskin[6]+0.8660254037844386*fskin[3]; 
  boundSurf[7] = 1.936491673103709*fskin[11]-1.5*fskin[7]+0.8660254037844386*fskin[5]; 
  boundSurf[8] = (-7.5*fskin[8])+5.809475019311125*fskin[2]-3.354101966249685*fskin[0]; 
  boundSurf[9] = (-7.5*fskin[9])+5.809475019311126*fskin[4]-3.354101966249684*fskin[1]; 
  boundSurf[10] = (-7.5*fskin[10])+5.809475019311126*fskin[6]-3.354101966249684*fskin[3]; 
  boundSurf[11] = (-7.5*fskin[11])+5.809475019311125*fskin[7]-3.354101966249685*fskin[5]; 

  edgeSurf_incr[0] = 0.7071067811865475*(edgeSurf[1]*nuVtSqSum[1]+edgeSurf[0]*nuVtSqSum[0]); 
  edgeSurf_incr[1] = 0.7071067811865475*(edgeSurf[0]*nuVtSqSum[1]+nuVtSqSum[0]*edgeSurf[1]); 
  edgeSurf_incr[2] = 0.7071067811865475*(nuVtSqSum[1]*edgeSurf[4]+nuVtSqSum[0]*edgeSurf[2]); 
  edgeSurf_incr[3] = 0.7071067811865475*(nuVtSqSum[1]*edgeSurf[5]+nuVtSqSum[0]*edgeSurf[3]); 
  edgeSurf_incr[4] = 0.7071067811865475*(nuVtSqSum[0]*edgeSurf[4]+nuVtSqSum[1]*edgeSurf[2]); 
  edgeSurf_incr[5] = 0.7071067811865475*(nuVtSqSum[0]*edgeSurf[5]+nuVtSqSum[1]*edgeSurf[3]); 
  edgeSurf_incr[6] = 0.7071067811865475*(nuVtSqSum[1]*edgeSurf[7]+nuVtSqSum[0]*edgeSurf[6]); 
  edgeSurf_incr[7] = 0.7071067811865475*(nuVtSqSum[0]*edgeSurf[7]+nuVtSqSum[1]*edgeSurf[6]); 
  edgeSurf_incr[8] = 0.04714045207910316*(15.0*nuVtSqSum[1]*edgeSurf[9]+15.0*nuVtSqSum[0]*edgeSurf[8]); 
  edgeSurf_incr[9] = 0.04714045207910316*(15.0*nuVtSqSum[0]*edgeSurf[9]+15.0*nuVtSqSum[1]*edgeSurf[8]); 
  edgeSurf_incr[10] = 0.04714045207910316*(15.0*nuVtSqSum[1]*edgeSurf[11]+15.0*nuVtSqSum[0]*edgeSurf[10]); 
  edgeSurf_incr[11] = 0.04714045207910316*(15.0*nuVtSqSum[0]*edgeSurf[11]+15.0*nuVtSqSum[1]*edgeSurf[10]); 

  boundSurf_incr[0] = 0.7071067811865475*(boundSurf[1]*nuVtSqSum[1]+boundSurf[0]*nuVtSqSum[0]); 
  boundSurf_incr[1] = 0.7071067811865475*(boundSurf[0]*nuVtSqSum[1]+nuVtSqSum[0]*boundSurf[1]); 
  boundSurf_incr[2] = 0.7071067811865475*(nuVtSqSum[1]*boundSurf[4]+nuVtSqSum[0]*boundSurf[2]); 
  boundSurf_incr[3] = 0.7071067811865475*(nuVtSqSum[1]*boundSurf[5]+nuVtSqSum[0]*boundSurf[3]); 
  boundSurf_incr[4] = 0.7071067811865475*(nuVtSqSum[0]*boundSurf[4]+nuVtSqSum[1]*boundSurf[2]); 
  boundSurf_incr[5] = 0.7071067811865475*(nuVtSqSum[0]*boundSurf[5]+nuVtSqSum[1]*boundSurf[3]); 
  boundSurf_incr[6] = 0.7071067811865475*(nuVtSqSum[1]*boundSurf[7]+nuVtSqSum[0]*boundSurf[6]); 
  boundSurf_incr[7] = 0.7071067811865475*(nuVtSqSum[0]*boundSurf[7]+nuVtSqSum[1]*boundSurf[6]); 
  boundSurf_incr[8] = 0.04714045207910316*(15.0*nuVtSqSum[1]*boundSurf[9]+15.0*nuVtSqSum[0]*boundSurf[8]); 
  boundSurf_incr[9] = 0.04714045207910316*(15.0*nuVtSqSum[0]*boundSurf[9]+15.0*nuVtSqSum[1]*boundSurf[8]); 
  boundSurf_incr[10] = 0.04714045207910316*(15.0*nuVtSqSum[1]*boundSurf[11]+15.0*nuVtSqSum[0]*boundSurf[10]); 
  boundSurf_incr[11] = 0.04714045207910316*(15.0*nuVtSqSum[0]*boundSurf[11]+15.0*nuVtSqSum[1]*boundSurf[10]); 


  } else { 

  double edgeSurf[12] = {0.0}; 
  edgeSurf[0] = -0.0125*(53.66563145999496*fskin[8]-53.66563145999496*fedge[8]-95.26279441628824*(fskin[2]+fedge[2])+75.0*fskin[0]-75.0*fedge[0]); 
  edgeSurf[1] = -0.0125*(53.66563145999495*fskin[9]-53.66563145999495*fedge[9]-95.26279441628824*(fskin[4]+fedge[4])+75.0*fskin[1]-75.0*fedge[1]); 
  edgeSurf[2] = 0.003125*(507.3608183531716*fskin[8]-236.2519841186524*fedge[8]-855.0*fskin[2]-465.0*fedge[2]+658.1793068761733*fskin[0]-381.051177665153*fedge[0]); 
  edgeSurf[3] = -0.0125*(53.66563145999495*fskin[10]-53.66563145999495*fedge[10]-95.26279441628824*(fskin[6]+fedge[6])+75.0*fskin[3]-75.0*fedge[3]); 
  edgeSurf[4] = 0.003125*(507.3608183531716*fskin[9]-236.2519841186524*fedge[9]-855.0*fskin[4]-465.0*fedge[4]+658.1793068761733*fskin[1]-381.051177665153*fedge[1]); 
  edgeSurf[5] = -0.0125*(53.66563145999496*fskin[11]-53.66563145999496*fedge[11]-95.26279441628824*(fskin[7]+fedge[7])+75.0*fskin[5]-75.0*fedge[5]); 
  edgeSurf[6] = 0.003125*(507.3608183531716*fskin[10]-236.2519841186524*fedge[10]-855.0*fskin[6]-465.0*fedge[6]+658.1793068761733*fskin[3]-381.051177665153*fedge[3]); 
  edgeSurf[7] = 0.003125*(507.3608183531716*fskin[11]-236.2519841186524*fedge[11]-855.0*fskin[7]-465.0*fedge[7]+658.1793068761733*fskin[5]-381.051177665153*fedge[5]); 
  edgeSurf[8] = -0.015625*(201.0*fskin[8]+9.0*fedge[8]-321.4576177352156*fskin[2]-19.36491673103709*fedge[2]+241.4953415699773*fskin[0]-26.83281572999748*fedge[0]); 
  edgeSurf[9] = -0.015625*(201.0*fskin[9]+9.0*fedge[9]-321.4576177352156*fskin[4]-19.36491673103708*fedge[4]+241.4953415699772*fskin[1]-26.83281572999747*fedge[1]); 
  edgeSurf[10] = -0.015625*(201.0*fskin[10]+9.0*fedge[10]-321.4576177352156*fskin[6]-19.36491673103708*fedge[6]+241.4953415699772*fskin[3]-26.83281572999747*fedge[3]); 
  edgeSurf[11] = -0.015625*(201.0*fskin[11]+9.0*fedge[11]-321.4576177352156*fskin[7]-19.36491673103709*fedge[7]+241.4953415699773*fskin[5]-26.83281572999748*fedge[5]); 

  double boundSurf[12] = {0.0}; 
  boundSurf[2] = -0.5*(3.872983346207417*fskin[8]+3.0*fskin[2]+1.732050807568877*fskin[0]); 
  boundSurf[4] = -0.5*(3.872983346207417*fskin[9]+3.0*fskin[4]+1.732050807568877*fskin[1]); 
  boundSurf[6] = -0.5*(3.872983346207417*fskin[10]+3.0*fskin[6]+1.732050807568877*fskin[3]); 
  boundSurf[7] = -0.5*(3.872983346207417*fskin[11]+3.0*fskin[7]+1.732050807568877*fskin[5]); 
  boundSurf[8] = -0.5*(15.0*fskin[8]+11.61895003862225*fskin[2]+6.708203932499369*fskin[0]); 
  boundSurf[9] = -0.5*(15.0*fskin[9]+11.61895003862225*fskin[4]+6.708203932499369*fskin[1]); 
  boundSurf[10] = -0.5*(15.0*fskin[10]+11.61895003862225*fskin[6]+6.708203932499369*fskin[3]); 
  boundSurf[11] = -0.5*(15.0*fskin[11]+11.61895003862225*fskin[7]+6.708203932499369*fskin[5]); 

  edgeSurf_incr[0] = 0.7071067811865475*(edgeSurf[1]*nuVtSqSum[1]+edgeSurf[0]*nuVtSqSum[0]); 
  edgeSurf_incr[1] = 0.7071067811865475*(edgeSurf[0]*nuVtSqSum[1]+nuVtSqSum[0]*edgeSurf[1]); 
  edgeSurf_incr[2] = 0.7071067811865475*(nuVtSqSum[1]*edgeSurf[4]+nuVtSqSum[0]*edgeSurf[2]); 
  edgeSurf_incr[3] = 0.7071067811865475*(nuVtSqSum[1]*edgeSurf[5]+nuVtSqSum[0]*edgeSurf[3]); 
  edgeSurf_incr[4] = 0.7071067811865475*(nuVtSqSum[0]*edgeSurf[4]+nuVtSqSum[1]*edgeSurf[2]); 
  edgeSurf_incr[5] = 0.7071067811865475*(nuVtSqSum[0]*edgeSurf[5]+nuVtSqSum[1]*edgeSurf[3]); 
  edgeSurf_incr[6] = 0.7071067811865475*(nuVtSqSum[1]*edgeSurf[7]+nuVtSqSum[0]*edgeSurf[6]); 
  edgeSurf_incr[7] = 0.7071067811865475*(nuVtSqSum[0]*edgeSurf[7]+nuVtSqSum[1]*edgeSurf[6]); 
  edgeSurf_incr[8] = 0.04714045207910316*(15.0*nuVtSqSum[1]*edgeSurf[9]+15.0*nuVtSqSum[0]*edgeSurf[8]); 
  edgeSurf_incr[9] = 0.04714045207910316*(15.0*nuVtSqSum[0]*edgeSurf[9]+15.0*nuVtSqSum[1]*edgeSurf[8]); 
  edgeSurf_incr[10] = 0.04714045207910316*(15.0*nuVtSqSum[1]*edgeSurf[11]+15.0*nuVtSqSum[0]*edgeSurf[10]); 
  edgeSurf_incr[11] = 0.04714045207910316*(15.0*nuVtSqSum[0]*edgeSurf[11]+15.0*nuVtSqSum[1]*edgeSurf[10]); 

  boundSurf_incr[0] = 0.7071067811865475*(boundSurf[1]*nuVtSqSum[1]+boundSurf[0]*nuVtSqSum[0]); 
  boundSurf_incr[1] = 0.7071067811865475*(boundSurf[0]*nuVtSqSum[1]+nuVtSqSum[0]*boundSurf[1]); 
  boundSurf_incr[2] = 0.7071067811865475*(nuVtSqSum[1]*boundSurf[4]+nuVtSqSum[0]*boundSurf[2]); 
  boundSurf_incr[3] = 0.7071067811865475*(nuVtSqSum[1]*boundSurf[5]+nuVtSqSum[0]*boundSurf[3]); 
  boundSurf_incr[4] = 0.7071067811865475*(nuVtSqSum[0]*boundSurf[4]+nuVtSqSum[1]*boundSurf[2]); 
  boundSurf_incr[5] = 0.7071067811865475*(nuVtSqSum[0]*boundSurf[5]+nuVtSqSum[1]*boundSurf[3]); 
  boundSurf_incr[6] = 0.7071067811865475*(nuVtSqSum[1]*boundSurf[7]+nuVtSqSum[0]*boundSurf[6]); 
  boundSurf_incr[7] = 0.7071067811865475*(nuVtSqSum[0]*boundSurf[7]+nuVtSqSum[1]*boundSurf[6]); 
  boundSurf_incr[8] = 0.04714045207910316*(15.0*nuVtSqSum[1]*boundSurf[9]+15.0*nuVtSqSum[0]*boundSurf[8]); 
  boundSurf_incr[9] = 0.04714045207910316*(15.0*nuVtSqSum[0]*boundSurf[9]+15.0*nuVtSqSum[1]*boundSurf[8]); 
  boundSurf_incr[10] = 0.04714045207910316*(15.0*nuVtSqSum[1]*boundSurf[11]+15.0*nuVtSqSum[0]*boundSurf[10]); 
  boundSurf_incr[11] = 0.04714045207910316*(15.0*nuVtSqSum[0]*boundSurf[11]+15.0*nuVtSqSum[1]*boundSurf[10]); 

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
  return 0.;

} 
