#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH void lbo_gyrokinetic_diff_boundary_surfvpar_1x2v_tensor_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
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

  double rdvSq4 = 4.0/(dxv[1]*dxv[1]); 

  double bFacFskin[27] = {0.0}; 
  bFacFskin[8] = 6.708203932499369*fskin[0]; 
  bFacFskin[12] = 6.708203932499369*fskin[1]; 
  bFacFskin[14] = 6.708203932499369*fskin[3]; 
  bFacFskin[18] = 6.708203932499369*fskin[5]; 
  bFacFskin[20] = 6.708203932499369*fskin[7]; 
  bFacFskin[22] = 6.708203932499369*fskin[9]; 
  bFacFskin[23] = 6.708203932499369*fskin[13]; 
  bFacFskin[25] = 6.708203932499369*fskin[15]; 
  bFacFskin[26] = 6.708203932499369*fskin[21]; 

  double vol_incr[27] = {0.0};
  vol_incr[8] = 0.105409255338946*(6.708203932499369*nuVtSqSum[2]*bFacFskin[20]+6.708203932499369*nuVtSqSum[1]*bFacFskin[12]+6.708203932499369*nuVtSqSum[0]*bFacFskin[8]); 
  vol_incr[12] = 0.03651483716701107*(17.32050807568877*nuVtSqSum[1]*bFacFskin[20]+(17.32050807568877*nuVtSqSum[2]+19.36491673103708*nuVtSqSum[0])*bFacFskin[12]+19.36491673103709*nuVtSqSum[1]*bFacFskin[8]); 
  vol_incr[14] = 0.7071067811865474*(nuVtSqSum[2]*bFacFskin[23]+nuVtSqSum[1]*bFacFskin[18]+nuVtSqSum[0]*bFacFskin[14]); 
  vol_incr[18] = 0.105409255338946*(6.0*nuVtSqSum[1]*bFacFskin[23]+(6.0*nuVtSqSum[2]+6.708203932499369*nuVtSqSum[0])*bFacFskin[18]+6.708203932499369*nuVtSqSum[1]*bFacFskin[14]); 
  vol_incr[20] = 0.006734350297014738*((67.0820393249937*nuVtSqSum[2]+105.0*nuVtSqSum[0])*bFacFskin[20]+93.91485505499116*nuVtSqSum[1]*bFacFskin[12]+105.0*nuVtSqSum[2]*bFacFskin[8]); 
  vol_incr[22] = 0.7071067811865475*(nuVtSqSum[2]*bFacFskin[26]+nuVtSqSum[1]*bFacFskin[25]+nuVtSqSum[0]*bFacFskin[22]); 
  vol_incr[23] = 0.01166423687039609*((38.72983346207418*nuVtSqSum[2]+60.6217782649107*nuVtSqSum[0])*bFacFskin[23]+54.22176684690384*nuVtSqSum[1]*bFacFskin[18]+60.62177826491071*nuVtSqSum[2]*bFacFskin[14]); 
  vol_incr[25] = 0.1414213562373095*(4.47213595499958*nuVtSqSum[1]*bFacFskin[26]+(4.47213595499958*nuVtSqSum[2]+5.0*nuVtSqSum[0])*bFacFskin[25]+5.0*nuVtSqSum[1]*bFacFskin[22]); 
  vol_incr[26] = 0.04517539514526255*((10.0*nuVtSqSum[2]+15.65247584249853*nuVtSqSum[0])*bFacFskin[26]+14.0*nuVtSqSum[1]*bFacFskin[25]+15.65247584249853*nuVtSqSum[2]*bFacFskin[22]); 

  double edgeSurf_incr[27] = {0.0}; 
  double boundSurf_incr[27] = {0.0}; 

  if (edge == -1) { 

  double edgeSurf[27] = {0.0}; 
  edgeSurf[0] = (-0.6708203932499369*fskin[8])+0.6708203932499369*fedge[8]-1.190784930203603*fskin[2]-1.190784930203603*fedge[2]-0.9375*fskin[0]+0.9375*fedge[0]; 
  edgeSurf[1] = (-0.6708203932499369*fskin[12])+0.6708203932499369*fedge[12]-1.190784930203603*fskin[4]-1.190784930203603*fedge[4]-0.9375*fskin[1]+0.9375*fedge[1]; 
  edgeSurf[2] = (-1.585502557353661*fskin[8])+0.7382874503707888*fedge[8]-2.671875*fskin[2]-1.453125*fedge[2]-2.056810333988042*fskin[0]+1.190784930203603*fedge[0]; 
  edgeSurf[3] = (-0.6708203932499369*fskin[14])+0.6708203932499369*fedge[14]-1.190784930203603*fskin[6]-1.190784930203603*fedge[6]-0.9375*fskin[3]+0.9375*fedge[3]; 
  edgeSurf[4] = (-1.585502557353661*fskin[12])+0.7382874503707888*fedge[12]-2.671875*fskin[4]-1.453125*fedge[4]-2.056810333988042*fskin[1]+1.190784930203603*fedge[1]; 
  edgeSurf[5] = (-0.6708203932499369*fskin[18])+0.6708203932499369*fedge[18]-1.190784930203603*fskin[10]-1.190784930203603*fedge[10]-0.9375*fskin[5]+0.9375*fedge[5]; 
  edgeSurf[6] = (-1.585502557353661*fskin[14])+0.7382874503707888*fedge[14]-2.671875*fskin[6]-1.453125*fedge[6]-2.056810333988042*fskin[3]+1.190784930203603*fedge[3]; 
  edgeSurf[7] = (-0.6708203932499369*fskin[20])+0.6708203932499369*fedge[20]-1.190784930203603*fskin[11]-1.190784930203603*fedge[11]-0.9375*fskin[7]+0.9375*fedge[7]; 
  edgeSurf[8] = (-3.140625*fskin[8])-0.140625*fedge[8]-5.022775277112744*fskin[2]-0.3025768239224545*fedge[2]-3.773364712030896*fskin[0]+0.4192627457812106*fedge[0]; 
  edgeSurf[9] = (-0.6708203932499369*fskin[22])+0.6708203932499369*fedge[22]-1.190784930203603*fskin[16]-1.190784930203603*fedge[16]-0.9375*fskin[9]+0.9375*fedge[9]; 
  edgeSurf[10] = (-1.585502557353661*fskin[18])+0.7382874503707888*fedge[18]-2.671875*fskin[10]-1.453125*fedge[10]-2.056810333988042*fskin[5]+1.190784930203603*fedge[5]; 
  edgeSurf[11] = (-1.585502557353661*fskin[20])+0.7382874503707888*fedge[20]-2.671875*fskin[11]-1.453125*fedge[11]-2.056810333988042*fskin[7]+1.190784930203603*fedge[7]; 
  edgeSurf[12] = (-3.140625*fskin[12])-0.140625*fedge[12]-5.022775277112744*fskin[4]-0.3025768239224544*fedge[4]-3.773364712030894*fskin[1]+0.4192627457812105*fedge[1]; 
  edgeSurf[13] = (-0.6708203932499369*fskin[23])+0.6708203932499369*fedge[23]-1.190784930203603*fskin[17]-1.190784930203603*fedge[17]-0.9375*fskin[13]+0.9375*fedge[13]; 
  edgeSurf[14] = (-3.140625*fskin[14])-0.140625*fedge[14]-5.022775277112744*fskin[6]-0.3025768239224544*fedge[6]-3.773364712030894*fskin[3]+0.4192627457812105*fedge[3]; 
  edgeSurf[15] = (-0.6708203932499369*fskin[25])+0.6708203932499369*fedge[25]-1.190784930203603*fskin[19]-1.190784930203603*fedge[19]-0.9375*fskin[15]+0.9375*fedge[15]; 
  edgeSurf[16] = (-1.585502557353661*fskin[22])+0.7382874503707888*fedge[22]-2.671875*fskin[16]-1.453125*fedge[16]-2.056810333988042*fskin[9]+1.190784930203603*fedge[9]; 
  edgeSurf[17] = (-1.585502557353661*fskin[23])+0.7382874503707888*fedge[23]-2.671875*fskin[17]-1.453125*fedge[17]-2.056810333988042*fskin[13]+1.190784930203603*fedge[13]; 
  edgeSurf[18] = (-3.140625*fskin[18])-0.140625*fedge[18]-5.022775277112744*fskin[10]-0.3025768239224545*fedge[10]-3.773364712030896*fskin[5]+0.4192627457812106*fedge[5]; 
  edgeSurf[19] = (-1.585502557353661*fskin[25])+0.7382874503707888*fedge[25]-2.671875*fskin[19]-1.453125*fedge[19]-2.056810333988042*fskin[15]+1.190784930203603*fedge[15]; 
  edgeSurf[20] = (-3.140625*fskin[20])-0.140625*fedge[20]-5.022775277112744*fskin[11]-0.3025768239224544*fedge[11]-3.773364712030896*fskin[7]+0.4192627457812106*fedge[7]; 
  edgeSurf[21] = (-0.6708203932499369*fskin[26])+0.6708203932499369*fedge[26]-1.190784930203603*fskin[24]-1.190784930203603*fedge[24]-0.9375*fskin[21]+0.9375*fedge[21]; 
  edgeSurf[22] = (-3.140625*fskin[22])-0.140625*fedge[22]-5.022775277112744*fskin[16]-0.3025768239224544*fedge[16]-3.773364712030896*fskin[9]+0.4192627457812106*fedge[9]; 
  edgeSurf[23] = (-3.140625*fskin[23])-0.140625*fedge[23]-5.022775277112744*fskin[17]-0.3025768239224545*fedge[17]-3.773364712030894*fskin[13]+0.4192627457812105*fedge[13]; 
  edgeSurf[24] = (-1.585502557353661*fskin[26])+0.7382874503707888*fedge[26]-2.671875*fskin[24]-1.453125*fedge[24]-2.056810333988042*fskin[21]+1.190784930203603*fedge[21]; 
  edgeSurf[25] = (-3.140625*fskin[25])-0.140625*fedge[25]-5.022775277112744*fskin[19]-0.3025768239224545*fedge[19]-3.773364712030894*fskin[15]+0.4192627457812105*fedge[15]; 
  edgeSurf[26] = (-3.140625*fskin[26])-0.140625*fedge[26]-5.022775277112744*fskin[24]-0.3025768239224545*fedge[24]-3.773364712030896*fskin[21]+0.4192627457812106*fedge[21]; 

  double boundSurf[27] = {0.0}; 
  boundSurf[2] = 1.936491673103709*fskin[8]-1.5*fskin[2]+0.8660254037844386*fskin[0]; 
  boundSurf[4] = 1.936491673103709*fskin[12]-1.5*fskin[4]+0.8660254037844386*fskin[1]; 
  boundSurf[6] = 1.936491673103709*fskin[14]-1.5*fskin[6]+0.8660254037844386*fskin[3]; 
  boundSurf[8] = (-7.5*fskin[8])+5.809475019311125*fskin[2]-3.354101966249685*fskin[0]; 
  boundSurf[10] = 1.936491673103709*fskin[18]-1.5*fskin[10]+0.8660254037844386*fskin[5]; 
  boundSurf[11] = 1.936491673103709*fskin[20]-1.5*fskin[11]+0.8660254037844387*fskin[7]; 
  boundSurf[12] = (-7.5*fskin[12])+5.809475019311126*fskin[4]-3.354101966249684*fskin[1]; 
  boundSurf[14] = (-7.5*fskin[14])+5.809475019311126*fskin[6]-3.354101966249684*fskin[3]; 
  boundSurf[16] = 1.936491673103709*fskin[22]-1.5*fskin[16]+0.8660254037844387*fskin[9]; 
  boundSurf[17] = 1.936491673103709*fskin[23]-1.5*fskin[17]+0.8660254037844387*fskin[13]; 
  boundSurf[18] = (-7.5*fskin[18])+5.809475019311125*fskin[10]-3.354101966249685*fskin[5]; 
  boundSurf[19] = 1.936491673103709*fskin[25]-1.5*fskin[19]+0.8660254037844387*fskin[15]; 
  boundSurf[20] = (-7.5*fskin[20])+5.809475019311126*fskin[11]-3.354101966249685*fskin[7]; 
  boundSurf[22] = (-7.5*fskin[22])+5.809475019311126*fskin[16]-3.354101966249685*fskin[9]; 
  boundSurf[23] = (-7.5*fskin[23])+5.809475019311125*fskin[17]-3.354101966249684*fskin[13]; 
  boundSurf[24] = 1.936491673103709*fskin[26]-1.5*fskin[24]+0.8660254037844386*fskin[21]; 
  boundSurf[25] = (-7.5*fskin[25])+5.809475019311125*fskin[19]-3.354101966249684*fskin[15]; 
  boundSurf[26] = (-7.5*fskin[26])+5.809475019311125*fskin[24]-3.354101966249685*fskin[21]; 

  edgeSurf_incr[0] = 0.7071067811865475*(nuVtSqSum[2]*edgeSurf[7]+edgeSurf[1]*nuVtSqSum[1]+edgeSurf[0]*nuVtSqSum[0]); 
  edgeSurf_incr[1] = 0.1414213562373095*(4.47213595499958*(nuVtSqSum[1]*edgeSurf[7]+edgeSurf[1]*nuVtSqSum[2])+5.0*(edgeSurf[0]*nuVtSqSum[1]+nuVtSqSum[0]*edgeSurf[1])); 
  edgeSurf_incr[2] = 0.04714045207910316*(15.0*nuVtSqSum[2]*edgeSurf[11]+15.0*(nuVtSqSum[1]*edgeSurf[4]+nuVtSqSum[0]*edgeSurf[2])); 
  edgeSurf_incr[3] = 0.04714045207910316*(15.0*nuVtSqSum[2]*edgeSurf[13]+15.0*(nuVtSqSum[1]*edgeSurf[5]+nuVtSqSum[0]*edgeSurf[3])); 
  edgeSurf_incr[4] = 0.04714045207910316*(13.41640786499874*nuVtSqSum[1]*edgeSurf[11]+(13.41640786499874*nuVtSqSum[2]+15.0*nuVtSqSum[0])*edgeSurf[4]+15.0*nuVtSqSum[1]*edgeSurf[2]); 
  edgeSurf_incr[5] = 0.04714045207910316*(13.41640786499874*nuVtSqSum[1]*edgeSurf[13]+(13.41640786499874*nuVtSqSum[2]+15.0*nuVtSqSum[0])*edgeSurf[5]+15.0*nuVtSqSum[1]*edgeSurf[3]); 
  edgeSurf_incr[6] = 0.7071067811865475*(nuVtSqSum[2]*edgeSurf[17]+nuVtSqSum[1]*edgeSurf[10]+nuVtSqSum[0]*edgeSurf[6]); 
  edgeSurf_incr[7] = 0.02020305089104421*((22.3606797749979*nuVtSqSum[2]+35.0*nuVtSqSum[0])*edgeSurf[7]+35.0*edgeSurf[0]*nuVtSqSum[2]+31.30495168499706*edgeSurf[1]*nuVtSqSum[1]); 
  edgeSurf_incr[8] = 0.04714045207910316*(15.0*nuVtSqSum[2]*edgeSurf[20]+15.0*nuVtSqSum[1]*edgeSurf[12]+15.0*nuVtSqSum[0]*edgeSurf[8]); 
  edgeSurf_incr[9] = 0.04714045207910316*(15.0*nuVtSqSum[2]*edgeSurf[21]+15.0*nuVtSqSum[1]*edgeSurf[15]+15.0*nuVtSqSum[0]*edgeSurf[9]); 
  edgeSurf_incr[10] = 0.1414213562373095*(4.47213595499958*nuVtSqSum[1]*edgeSurf[17]+(4.47213595499958*nuVtSqSum[2]+5.0*nuVtSqSum[0])*edgeSurf[10]+5.0*nuVtSqSum[1]*edgeSurf[6]); 
  edgeSurf_incr[11] = 0.006734350297014738*((67.0820393249937*nuVtSqSum[2]+105.0*nuVtSqSum[0])*edgeSurf[11]+93.91485505499116*nuVtSqSum[1]*edgeSurf[4]+105.0*edgeSurf[2]*nuVtSqSum[2]); 
  edgeSurf_incr[12] = 0.04714045207910316*(13.41640786499874*nuVtSqSum[1]*edgeSurf[20]+(13.41640786499874*nuVtSqSum[2]+15.0*nuVtSqSum[0])*edgeSurf[12]+15.0*nuVtSqSum[1]*edgeSurf[8]); 
  edgeSurf_incr[13] = 0.006734350297014738*((67.0820393249937*nuVtSqSum[2]+105.0*nuVtSqSum[0])*edgeSurf[13]+93.91485505499116*nuVtSqSum[1]*edgeSurf[5]+105.0*nuVtSqSum[2]*edgeSurf[3]); 
  edgeSurf_incr[14] = 0.04714045207910316*(15.0*(nuVtSqSum[2]*edgeSurf[23]+nuVtSqSum[1]*edgeSurf[18])+15.0*nuVtSqSum[0]*edgeSurf[14]); 
  edgeSurf_incr[15] = 0.04714045207910316*(13.41640786499874*nuVtSqSum[1]*edgeSurf[21]+(13.41640786499874*nuVtSqSum[2]+15.0*nuVtSqSum[0])*edgeSurf[15]+15.0*nuVtSqSum[1]*edgeSurf[9]); 
  edgeSurf_incr[16] = 0.04714045207910316*(15.0*(nuVtSqSum[2]*edgeSurf[24]+nuVtSqSum[1]*edgeSurf[19])+15.0*nuVtSqSum[0]*edgeSurf[16]); 
  edgeSurf_incr[17] = 0.02020305089104421*((22.3606797749979*nuVtSqSum[2]+35.0*nuVtSqSum[0])*edgeSurf[17]+31.30495168499706*nuVtSqSum[1]*edgeSurf[10]+35.0*nuVtSqSum[2]*edgeSurf[6]); 
  edgeSurf_incr[18] = 0.04714045207910316*(13.41640786499874*nuVtSqSum[1]*edgeSurf[23]+(13.41640786499874*nuVtSqSum[2]+15.0*nuVtSqSum[0])*edgeSurf[18]+15.0*nuVtSqSum[1]*edgeSurf[14]); 
  edgeSurf_incr[19] = 0.04714045207910316*(13.41640786499874*nuVtSqSum[1]*edgeSurf[24]+(13.41640786499874*nuVtSqSum[2]+15.0*nuVtSqSum[0])*edgeSurf[19]+15.0*nuVtSqSum[1]*edgeSurf[16]); 
  edgeSurf_incr[20] = 0.006734350297014738*((67.0820393249937*nuVtSqSum[2]+105.0*nuVtSqSum[0])*edgeSurf[20]+93.91485505499116*nuVtSqSum[1]*edgeSurf[12]+105.0*nuVtSqSum[2]*edgeSurf[8]); 
  edgeSurf_incr[21] = 0.006734350297014738*((67.0820393249937*nuVtSqSum[2]+105.0*nuVtSqSum[0])*edgeSurf[21]+93.91485505499116*nuVtSqSum[1]*edgeSurf[15]+105.0*nuVtSqSum[2]*edgeSurf[9]); 
  edgeSurf_incr[22] = 0.7071067811865475*(nuVtSqSum[2]*edgeSurf[26]+nuVtSqSum[1]*edgeSurf[25]+nuVtSqSum[0]*edgeSurf[22]); 
  edgeSurf_incr[23] = 0.006734350297014738*((67.0820393249937*nuVtSqSum[2]+105.0*nuVtSqSum[0])*edgeSurf[23]+93.91485505499116*nuVtSqSum[1]*edgeSurf[18]+105.0*nuVtSqSum[2]*edgeSurf[14]); 
  edgeSurf_incr[24] = 0.006734350297014738*((67.0820393249937*nuVtSqSum[2]+105.0*nuVtSqSum[0])*edgeSurf[24]+93.91485505499116*nuVtSqSum[1]*edgeSurf[19]+105.0*nuVtSqSum[2]*edgeSurf[16]); 
  edgeSurf_incr[25] = 0.1414213562373095*(4.47213595499958*nuVtSqSum[1]*edgeSurf[26]+(4.47213595499958*nuVtSqSum[2]+5.0*nuVtSqSum[0])*edgeSurf[25]+5.0*nuVtSqSum[1]*edgeSurf[22]); 
  edgeSurf_incr[26] = 0.02020305089104421*((22.3606797749979*nuVtSqSum[2]+35.0*nuVtSqSum[0])*edgeSurf[26]+31.30495168499706*nuVtSqSum[1]*edgeSurf[25]+35.0*nuVtSqSum[2]*edgeSurf[22]); 

  boundSurf_incr[0] = 0.7071067811865475*(nuVtSqSum[2]*boundSurf[7]+boundSurf[1]*nuVtSqSum[1]+boundSurf[0]*nuVtSqSum[0]); 
  boundSurf_incr[1] = 0.1414213562373095*(4.47213595499958*(nuVtSqSum[1]*boundSurf[7]+boundSurf[1]*nuVtSqSum[2])+5.0*(boundSurf[0]*nuVtSqSum[1]+nuVtSqSum[0]*boundSurf[1])); 
  boundSurf_incr[2] = 0.04714045207910316*(15.0*nuVtSqSum[2]*boundSurf[11]+15.0*(nuVtSqSum[1]*boundSurf[4]+nuVtSqSum[0]*boundSurf[2])); 
  boundSurf_incr[3] = 0.04714045207910316*(15.0*nuVtSqSum[2]*boundSurf[13]+15.0*(nuVtSqSum[1]*boundSurf[5]+nuVtSqSum[0]*boundSurf[3])); 
  boundSurf_incr[4] = 0.04714045207910316*(13.41640786499874*nuVtSqSum[1]*boundSurf[11]+(13.41640786499874*nuVtSqSum[2]+15.0*nuVtSqSum[0])*boundSurf[4]+15.0*nuVtSqSum[1]*boundSurf[2]); 
  boundSurf_incr[5] = 0.04714045207910316*(13.41640786499874*nuVtSqSum[1]*boundSurf[13]+(13.41640786499874*nuVtSqSum[2]+15.0*nuVtSqSum[0])*boundSurf[5]+15.0*nuVtSqSum[1]*boundSurf[3]); 
  boundSurf_incr[6] = 0.7071067811865475*(nuVtSqSum[2]*boundSurf[17]+nuVtSqSum[1]*boundSurf[10]+nuVtSqSum[0]*boundSurf[6]); 
  boundSurf_incr[7] = 0.02020305089104421*((22.3606797749979*nuVtSqSum[2]+35.0*nuVtSqSum[0])*boundSurf[7]+35.0*boundSurf[0]*nuVtSqSum[2]+31.30495168499706*boundSurf[1]*nuVtSqSum[1]); 
  boundSurf_incr[8] = 0.04714045207910316*(15.0*nuVtSqSum[2]*boundSurf[20]+15.0*nuVtSqSum[1]*boundSurf[12]+15.0*nuVtSqSum[0]*boundSurf[8]); 
  boundSurf_incr[9] = 0.04714045207910316*(15.0*nuVtSqSum[2]*boundSurf[21]+15.0*nuVtSqSum[1]*boundSurf[15]+15.0*nuVtSqSum[0]*boundSurf[9]); 
  boundSurf_incr[10] = 0.1414213562373095*(4.47213595499958*nuVtSqSum[1]*boundSurf[17]+(4.47213595499958*nuVtSqSum[2]+5.0*nuVtSqSum[0])*boundSurf[10]+5.0*nuVtSqSum[1]*boundSurf[6]); 
  boundSurf_incr[11] = 0.006734350297014738*((67.0820393249937*nuVtSqSum[2]+105.0*nuVtSqSum[0])*boundSurf[11]+93.91485505499116*nuVtSqSum[1]*boundSurf[4]+105.0*boundSurf[2]*nuVtSqSum[2]); 
  boundSurf_incr[12] = 0.04714045207910316*(13.41640786499874*nuVtSqSum[1]*boundSurf[20]+(13.41640786499874*nuVtSqSum[2]+15.0*nuVtSqSum[0])*boundSurf[12]+15.0*nuVtSqSum[1]*boundSurf[8]); 
  boundSurf_incr[13] = 0.006734350297014738*((67.0820393249937*nuVtSqSum[2]+105.0*nuVtSqSum[0])*boundSurf[13]+93.91485505499116*nuVtSqSum[1]*boundSurf[5]+105.0*nuVtSqSum[2]*boundSurf[3]); 
  boundSurf_incr[14] = 0.04714045207910316*(15.0*(nuVtSqSum[2]*boundSurf[23]+nuVtSqSum[1]*boundSurf[18])+15.0*nuVtSqSum[0]*boundSurf[14]); 
  boundSurf_incr[15] = 0.04714045207910316*(13.41640786499874*nuVtSqSum[1]*boundSurf[21]+(13.41640786499874*nuVtSqSum[2]+15.0*nuVtSqSum[0])*boundSurf[15]+15.0*nuVtSqSum[1]*boundSurf[9]); 
  boundSurf_incr[16] = 0.04714045207910316*(15.0*(nuVtSqSum[2]*boundSurf[24]+nuVtSqSum[1]*boundSurf[19])+15.0*nuVtSqSum[0]*boundSurf[16]); 
  boundSurf_incr[17] = 0.02020305089104421*((22.3606797749979*nuVtSqSum[2]+35.0*nuVtSqSum[0])*boundSurf[17]+31.30495168499706*nuVtSqSum[1]*boundSurf[10]+35.0*nuVtSqSum[2]*boundSurf[6]); 
  boundSurf_incr[18] = 0.04714045207910316*(13.41640786499874*nuVtSqSum[1]*boundSurf[23]+(13.41640786499874*nuVtSqSum[2]+15.0*nuVtSqSum[0])*boundSurf[18]+15.0*nuVtSqSum[1]*boundSurf[14]); 
  boundSurf_incr[19] = 0.04714045207910316*(13.41640786499874*nuVtSqSum[1]*boundSurf[24]+(13.41640786499874*nuVtSqSum[2]+15.0*nuVtSqSum[0])*boundSurf[19]+15.0*nuVtSqSum[1]*boundSurf[16]); 
  boundSurf_incr[20] = 0.006734350297014738*((67.0820393249937*nuVtSqSum[2]+105.0*nuVtSqSum[0])*boundSurf[20]+93.91485505499116*nuVtSqSum[1]*boundSurf[12]+105.0*nuVtSqSum[2]*boundSurf[8]); 
  boundSurf_incr[21] = 0.006734350297014738*((67.0820393249937*nuVtSqSum[2]+105.0*nuVtSqSum[0])*boundSurf[21]+93.91485505499116*nuVtSqSum[1]*boundSurf[15]+105.0*nuVtSqSum[2]*boundSurf[9]); 
  boundSurf_incr[22] = 0.7071067811865475*(nuVtSqSum[2]*boundSurf[26]+nuVtSqSum[1]*boundSurf[25]+nuVtSqSum[0]*boundSurf[22]); 
  boundSurf_incr[23] = 0.006734350297014738*((67.0820393249937*nuVtSqSum[2]+105.0*nuVtSqSum[0])*boundSurf[23]+93.91485505499116*nuVtSqSum[1]*boundSurf[18]+105.0*nuVtSqSum[2]*boundSurf[14]); 
  boundSurf_incr[24] = 0.006734350297014738*((67.0820393249937*nuVtSqSum[2]+105.0*nuVtSqSum[0])*boundSurf[24]+93.91485505499116*nuVtSqSum[1]*boundSurf[19]+105.0*nuVtSqSum[2]*boundSurf[16]); 
  boundSurf_incr[25] = 0.1414213562373095*(4.47213595499958*nuVtSqSum[1]*boundSurf[26]+(4.47213595499958*nuVtSqSum[2]+5.0*nuVtSqSum[0])*boundSurf[25]+5.0*nuVtSqSum[1]*boundSurf[22]); 
  boundSurf_incr[26] = 0.02020305089104421*((22.3606797749979*nuVtSqSum[2]+35.0*nuVtSqSum[0])*boundSurf[26]+31.30495168499706*nuVtSqSum[1]*boundSurf[25]+35.0*nuVtSqSum[2]*boundSurf[22]); 


  } else { 

  double edgeSurf[27] = {0.0}; 
  edgeSurf[0] = -0.0125*(53.66563145999496*fskin[8]-53.66563145999496*fedge[8]-95.26279441628824*(fskin[2]+fedge[2])+75.0*fskin[0]-75.0*fedge[0]); 
  edgeSurf[1] = -0.0125*(53.66563145999495*fskin[12]-53.66563145999495*fedge[12]-95.26279441628824*(fskin[4]+fedge[4])+75.0*fskin[1]-75.0*fedge[1]); 
  edgeSurf[2] = 0.003125*(507.3608183531716*fskin[8]-236.2519841186524*fedge[8]-855.0*fskin[2]-465.0*fedge[2]+658.1793068761733*fskin[0]-381.051177665153*fedge[0]); 
  edgeSurf[3] = -0.0125*(53.66563145999495*fskin[14]-53.66563145999495*fedge[14]-95.26279441628824*(fskin[6]+fedge[6])+75.0*fskin[3]-75.0*fedge[3]); 
  edgeSurf[4] = 0.003125*(507.3608183531716*fskin[12]-236.2519841186524*fedge[12]-855.0*fskin[4]-465.0*fedge[4]+658.1793068761733*fskin[1]-381.051177665153*fedge[1]); 
  edgeSurf[5] = -0.0125*(53.66563145999496*fskin[18]-53.66563145999496*fedge[18]-95.26279441628824*(fskin[10]+fedge[10])+75.0*fskin[5]-75.0*fedge[5]); 
  edgeSurf[6] = 0.003125*(507.3608183531716*fskin[14]-236.2519841186524*fedge[14]-855.0*fskin[6]-465.0*fedge[6]+658.1793068761733*fskin[3]-381.051177665153*fedge[3]); 
  edgeSurf[7] = -0.0125*(53.66563145999496*fskin[20]-53.66563145999496*fedge[20]-95.26279441628826*(fskin[11]+fedge[11])+75.0*fskin[7]-75.0*fedge[7]); 
  edgeSurf[8] = -0.015625*(201.0*fskin[8]+9.0*fedge[8]-321.4576177352156*fskin[2]-19.36491673103709*fedge[2]+241.4953415699773*fskin[0]-26.83281572999748*fedge[0]); 
  edgeSurf[9] = -0.0125*(53.66563145999496*fskin[22]-53.66563145999496*fedge[22]-95.26279441628826*(fskin[16]+fedge[16])+75.0*fskin[9]-75.0*fedge[9]); 
  edgeSurf[10] = 0.003125*(507.3608183531716*fskin[18]-236.2519841186524*fedge[18]-855.0*fskin[10]-465.0*fedge[10]+658.1793068761733*fskin[5]-381.051177665153*fedge[5]); 
  edgeSurf[11] = 0.003125*(507.3608183531716*fskin[20]-236.2519841186524*fedge[20]-855.0*fskin[11]-465.0*fedge[11]+658.1793068761734*fskin[7]-381.051177665153*fedge[7]); 
  edgeSurf[12] = -0.015625*(201.0*fskin[12]+9.0*fedge[12]-321.4576177352156*fskin[4]-19.36491673103708*fedge[4]+241.4953415699772*fskin[1]-26.83281572999747*fedge[1]); 
  edgeSurf[13] = -0.0125*(53.66563145999495*fskin[23]-53.66563145999495*fedge[23]-95.26279441628826*(fskin[17]+fedge[17])+75.0*fskin[13]-75.0*fedge[13]); 
  edgeSurf[14] = -0.015625*(201.0*fskin[14]+9.0*fedge[14]-321.4576177352156*fskin[6]-19.36491673103708*fedge[6]+241.4953415699772*fskin[3]-26.83281572999747*fedge[3]); 
  edgeSurf[15] = -0.0125*(53.66563145999495*fskin[25]-53.66563145999495*fedge[25]-95.26279441628826*(fskin[19]+fedge[19])+75.0*fskin[15]-75.0*fedge[15]); 
  edgeSurf[16] = 0.003125*(507.3608183531716*fskin[22]-236.2519841186524*fedge[22]-855.0*fskin[16]-465.0*fedge[16]+658.1793068761734*fskin[9]-381.051177665153*fedge[9]); 
  edgeSurf[17] = 0.003125*(507.3608183531716*fskin[23]-236.2519841186524*fedge[23]-855.0*fskin[17]-465.0*fedge[17]+658.1793068761734*fskin[13]-381.051177665153*fedge[13]); 
  edgeSurf[18] = -0.015625*(201.0*fskin[18]+9.0*fedge[18]-321.4576177352156*fskin[10]-19.36491673103709*fedge[10]+241.4953415699773*fskin[5]-26.83281572999748*fedge[5]); 
  edgeSurf[19] = 0.003125*(507.3608183531716*fskin[25]-236.2519841186524*fedge[25]-855.0*fskin[19]-465.0*fedge[19]+658.1793068761734*fskin[15]-381.051177665153*fedge[15]); 
  edgeSurf[20] = -0.015625*(201.0*fskin[20]+9.0*fedge[20]-321.4576177352156*fskin[11]-19.36491673103708*fedge[11]+241.4953415699773*fskin[7]-26.83281572999748*fedge[7]); 
  edgeSurf[21] = -0.0125*(53.66563145999496*fskin[26]-53.66563145999496*fedge[26]-95.26279441628824*(fskin[24]+fedge[24])+75.0*fskin[21]-75.0*fedge[21]); 
  edgeSurf[22] = -0.015625*(201.0*fskin[22]+9.0*fedge[22]-321.4576177352156*fskin[16]-19.36491673103708*fedge[16]+241.4953415699773*fskin[9]-26.83281572999748*fedge[9]); 
  edgeSurf[23] = -0.015625*(201.0*fskin[23]+9.0*fedge[23]-321.4576177352156*fskin[17]-19.36491673103709*fedge[17]+241.4953415699772*fskin[13]-26.83281572999747*fedge[13]); 
  edgeSurf[24] = 0.003125*(507.3608183531716*fskin[26]-236.2519841186524*fedge[26]-855.0*fskin[24]-465.0*fedge[24]+658.1793068761733*fskin[21]-381.051177665153*fedge[21]); 
  edgeSurf[25] = -0.015625*(201.0*fskin[25]+9.0*fedge[25]-321.4576177352156*fskin[19]-19.36491673103709*fedge[19]+241.4953415699772*fskin[15]-26.83281572999747*fedge[15]); 
  edgeSurf[26] = -0.015625*(201.0*fskin[26]+9.0*fedge[26]-321.4576177352156*fskin[24]-19.36491673103709*fedge[24]+241.4953415699773*fskin[21]-26.83281572999748*fedge[21]); 

  double boundSurf[27] = {0.0}; 
  boundSurf[2] = -0.5*(3.872983346207417*fskin[8]+3.0*fskin[2]+1.732050807568877*fskin[0]); 
  boundSurf[4] = -0.5*(3.872983346207417*fskin[12]+3.0*fskin[4]+1.732050807568877*fskin[1]); 
  boundSurf[6] = -0.5*(3.872983346207417*fskin[14]+3.0*fskin[6]+1.732050807568877*fskin[3]); 
  boundSurf[8] = -0.5*(15.0*fskin[8]+11.61895003862225*fskin[2]+6.708203932499369*fskin[0]); 
  boundSurf[10] = -0.5*(3.872983346207417*fskin[18]+3.0*fskin[10]+1.732050807568877*fskin[5]); 
  boundSurf[11] = -0.1*(19.36491673103708*fskin[20]+15.0*fskin[11]+8.660254037844387*fskin[7]); 
  boundSurf[12] = -0.5*(15.0*fskin[12]+11.61895003862225*fskin[4]+6.708203932499369*fskin[1]); 
  boundSurf[14] = -0.5*(15.0*fskin[14]+11.61895003862225*fskin[6]+6.708203932499369*fskin[3]); 
  boundSurf[16] = -0.1*(19.36491673103708*fskin[22]+15.0*fskin[16]+8.660254037844387*fskin[9]); 
  boundSurf[17] = -0.1*(19.36491673103709*fskin[23]+15.0*fskin[17]+8.660254037844387*fskin[13]); 
  boundSurf[18] = -0.5*(15.0*fskin[18]+11.61895003862225*fskin[10]+6.708203932499369*fskin[5]); 
  boundSurf[19] = -0.1*(19.36491673103709*fskin[25]+15.0*fskin[19]+8.660254037844387*fskin[15]); 
  boundSurf[20] = -0.5*(15.0*fskin[20]+11.61895003862225*fskin[11]+6.708203932499369*fskin[7]); 
  boundSurf[22] = -0.5*(15.0*fskin[22]+11.61895003862225*fskin[16]+6.708203932499369*fskin[9]); 
  boundSurf[23] = -0.5*(15.0*fskin[23]+11.61895003862225*fskin[17]+6.708203932499369*fskin[13]); 
  boundSurf[24] = -0.5*(3.872983346207417*fskin[26]+3.0*fskin[24]+1.732050807568877*fskin[21]); 
  boundSurf[25] = -0.5*(15.0*fskin[25]+11.61895003862225*fskin[19]+6.708203932499369*fskin[15]); 
  boundSurf[26] = -0.5*(15.0*fskin[26]+11.61895003862225*fskin[24]+6.708203932499369*fskin[21]); 

  edgeSurf_incr[0] = 0.7071067811865475*(nuVtSqSum[2]*edgeSurf[7]+edgeSurf[1]*nuVtSqSum[1]+edgeSurf[0]*nuVtSqSum[0]); 
  edgeSurf_incr[1] = 0.1414213562373095*(4.47213595499958*(nuVtSqSum[1]*edgeSurf[7]+edgeSurf[1]*nuVtSqSum[2])+5.0*(edgeSurf[0]*nuVtSqSum[1]+nuVtSqSum[0]*edgeSurf[1])); 
  edgeSurf_incr[2] = 0.04714045207910316*(15.0*nuVtSqSum[2]*edgeSurf[11]+15.0*(nuVtSqSum[1]*edgeSurf[4]+nuVtSqSum[0]*edgeSurf[2])); 
  edgeSurf_incr[3] = 0.04714045207910316*(15.0*nuVtSqSum[2]*edgeSurf[13]+15.0*(nuVtSqSum[1]*edgeSurf[5]+nuVtSqSum[0]*edgeSurf[3])); 
  edgeSurf_incr[4] = 0.04714045207910316*(13.41640786499874*nuVtSqSum[1]*edgeSurf[11]+(13.41640786499874*nuVtSqSum[2]+15.0*nuVtSqSum[0])*edgeSurf[4]+15.0*nuVtSqSum[1]*edgeSurf[2]); 
  edgeSurf_incr[5] = 0.04714045207910316*(13.41640786499874*nuVtSqSum[1]*edgeSurf[13]+(13.41640786499874*nuVtSqSum[2]+15.0*nuVtSqSum[0])*edgeSurf[5]+15.0*nuVtSqSum[1]*edgeSurf[3]); 
  edgeSurf_incr[6] = 0.7071067811865475*(nuVtSqSum[2]*edgeSurf[17]+nuVtSqSum[1]*edgeSurf[10]+nuVtSqSum[0]*edgeSurf[6]); 
  edgeSurf_incr[7] = 0.02020305089104421*((22.3606797749979*nuVtSqSum[2]+35.0*nuVtSqSum[0])*edgeSurf[7]+35.0*edgeSurf[0]*nuVtSqSum[2]+31.30495168499706*edgeSurf[1]*nuVtSqSum[1]); 
  edgeSurf_incr[8] = 0.04714045207910316*(15.0*nuVtSqSum[2]*edgeSurf[20]+15.0*nuVtSqSum[1]*edgeSurf[12]+15.0*nuVtSqSum[0]*edgeSurf[8]); 
  edgeSurf_incr[9] = 0.04714045207910316*(15.0*nuVtSqSum[2]*edgeSurf[21]+15.0*nuVtSqSum[1]*edgeSurf[15]+15.0*nuVtSqSum[0]*edgeSurf[9]); 
  edgeSurf_incr[10] = 0.1414213562373095*(4.47213595499958*nuVtSqSum[1]*edgeSurf[17]+(4.47213595499958*nuVtSqSum[2]+5.0*nuVtSqSum[0])*edgeSurf[10]+5.0*nuVtSqSum[1]*edgeSurf[6]); 
  edgeSurf_incr[11] = 0.006734350297014738*((67.0820393249937*nuVtSqSum[2]+105.0*nuVtSqSum[0])*edgeSurf[11]+93.91485505499116*nuVtSqSum[1]*edgeSurf[4]+105.0*edgeSurf[2]*nuVtSqSum[2]); 
  edgeSurf_incr[12] = 0.04714045207910316*(13.41640786499874*nuVtSqSum[1]*edgeSurf[20]+(13.41640786499874*nuVtSqSum[2]+15.0*nuVtSqSum[0])*edgeSurf[12]+15.0*nuVtSqSum[1]*edgeSurf[8]); 
  edgeSurf_incr[13] = 0.006734350297014738*((67.0820393249937*nuVtSqSum[2]+105.0*nuVtSqSum[0])*edgeSurf[13]+93.91485505499116*nuVtSqSum[1]*edgeSurf[5]+105.0*nuVtSqSum[2]*edgeSurf[3]); 
  edgeSurf_incr[14] = 0.04714045207910316*(15.0*(nuVtSqSum[2]*edgeSurf[23]+nuVtSqSum[1]*edgeSurf[18])+15.0*nuVtSqSum[0]*edgeSurf[14]); 
  edgeSurf_incr[15] = 0.04714045207910316*(13.41640786499874*nuVtSqSum[1]*edgeSurf[21]+(13.41640786499874*nuVtSqSum[2]+15.0*nuVtSqSum[0])*edgeSurf[15]+15.0*nuVtSqSum[1]*edgeSurf[9]); 
  edgeSurf_incr[16] = 0.04714045207910316*(15.0*(nuVtSqSum[2]*edgeSurf[24]+nuVtSqSum[1]*edgeSurf[19])+15.0*nuVtSqSum[0]*edgeSurf[16]); 
  edgeSurf_incr[17] = 0.02020305089104421*((22.3606797749979*nuVtSqSum[2]+35.0*nuVtSqSum[0])*edgeSurf[17]+31.30495168499706*nuVtSqSum[1]*edgeSurf[10]+35.0*nuVtSqSum[2]*edgeSurf[6]); 
  edgeSurf_incr[18] = 0.04714045207910316*(13.41640786499874*nuVtSqSum[1]*edgeSurf[23]+(13.41640786499874*nuVtSqSum[2]+15.0*nuVtSqSum[0])*edgeSurf[18]+15.0*nuVtSqSum[1]*edgeSurf[14]); 
  edgeSurf_incr[19] = 0.04714045207910316*(13.41640786499874*nuVtSqSum[1]*edgeSurf[24]+(13.41640786499874*nuVtSqSum[2]+15.0*nuVtSqSum[0])*edgeSurf[19]+15.0*nuVtSqSum[1]*edgeSurf[16]); 
  edgeSurf_incr[20] = 0.006734350297014738*((67.0820393249937*nuVtSqSum[2]+105.0*nuVtSqSum[0])*edgeSurf[20]+93.91485505499116*nuVtSqSum[1]*edgeSurf[12]+105.0*nuVtSqSum[2]*edgeSurf[8]); 
  edgeSurf_incr[21] = 0.006734350297014738*((67.0820393249937*nuVtSqSum[2]+105.0*nuVtSqSum[0])*edgeSurf[21]+93.91485505499116*nuVtSqSum[1]*edgeSurf[15]+105.0*nuVtSqSum[2]*edgeSurf[9]); 
  edgeSurf_incr[22] = 0.7071067811865475*(nuVtSqSum[2]*edgeSurf[26]+nuVtSqSum[1]*edgeSurf[25]+nuVtSqSum[0]*edgeSurf[22]); 
  edgeSurf_incr[23] = 0.006734350297014738*((67.0820393249937*nuVtSqSum[2]+105.0*nuVtSqSum[0])*edgeSurf[23]+93.91485505499116*nuVtSqSum[1]*edgeSurf[18]+105.0*nuVtSqSum[2]*edgeSurf[14]); 
  edgeSurf_incr[24] = 0.006734350297014738*((67.0820393249937*nuVtSqSum[2]+105.0*nuVtSqSum[0])*edgeSurf[24]+93.91485505499116*nuVtSqSum[1]*edgeSurf[19]+105.0*nuVtSqSum[2]*edgeSurf[16]); 
  edgeSurf_incr[25] = 0.1414213562373095*(4.47213595499958*nuVtSqSum[1]*edgeSurf[26]+(4.47213595499958*nuVtSqSum[2]+5.0*nuVtSqSum[0])*edgeSurf[25]+5.0*nuVtSqSum[1]*edgeSurf[22]); 
  edgeSurf_incr[26] = 0.02020305089104421*((22.3606797749979*nuVtSqSum[2]+35.0*nuVtSqSum[0])*edgeSurf[26]+31.30495168499706*nuVtSqSum[1]*edgeSurf[25]+35.0*nuVtSqSum[2]*edgeSurf[22]); 

  boundSurf_incr[0] = 0.7071067811865475*(nuVtSqSum[2]*boundSurf[7]+boundSurf[1]*nuVtSqSum[1]+boundSurf[0]*nuVtSqSum[0]); 
  boundSurf_incr[1] = 0.1414213562373095*(4.47213595499958*(nuVtSqSum[1]*boundSurf[7]+boundSurf[1]*nuVtSqSum[2])+5.0*(boundSurf[0]*nuVtSqSum[1]+nuVtSqSum[0]*boundSurf[1])); 
  boundSurf_incr[2] = 0.04714045207910316*(15.0*nuVtSqSum[2]*boundSurf[11]+15.0*(nuVtSqSum[1]*boundSurf[4]+nuVtSqSum[0]*boundSurf[2])); 
  boundSurf_incr[3] = 0.04714045207910316*(15.0*nuVtSqSum[2]*boundSurf[13]+15.0*(nuVtSqSum[1]*boundSurf[5]+nuVtSqSum[0]*boundSurf[3])); 
  boundSurf_incr[4] = 0.04714045207910316*(13.41640786499874*nuVtSqSum[1]*boundSurf[11]+(13.41640786499874*nuVtSqSum[2]+15.0*nuVtSqSum[0])*boundSurf[4]+15.0*nuVtSqSum[1]*boundSurf[2]); 
  boundSurf_incr[5] = 0.04714045207910316*(13.41640786499874*nuVtSqSum[1]*boundSurf[13]+(13.41640786499874*nuVtSqSum[2]+15.0*nuVtSqSum[0])*boundSurf[5]+15.0*nuVtSqSum[1]*boundSurf[3]); 
  boundSurf_incr[6] = 0.7071067811865475*(nuVtSqSum[2]*boundSurf[17]+nuVtSqSum[1]*boundSurf[10]+nuVtSqSum[0]*boundSurf[6]); 
  boundSurf_incr[7] = 0.02020305089104421*((22.3606797749979*nuVtSqSum[2]+35.0*nuVtSqSum[0])*boundSurf[7]+35.0*boundSurf[0]*nuVtSqSum[2]+31.30495168499706*boundSurf[1]*nuVtSqSum[1]); 
  boundSurf_incr[8] = 0.04714045207910316*(15.0*nuVtSqSum[2]*boundSurf[20]+15.0*nuVtSqSum[1]*boundSurf[12]+15.0*nuVtSqSum[0]*boundSurf[8]); 
  boundSurf_incr[9] = 0.04714045207910316*(15.0*nuVtSqSum[2]*boundSurf[21]+15.0*nuVtSqSum[1]*boundSurf[15]+15.0*nuVtSqSum[0]*boundSurf[9]); 
  boundSurf_incr[10] = 0.1414213562373095*(4.47213595499958*nuVtSqSum[1]*boundSurf[17]+(4.47213595499958*nuVtSqSum[2]+5.0*nuVtSqSum[0])*boundSurf[10]+5.0*nuVtSqSum[1]*boundSurf[6]); 
  boundSurf_incr[11] = 0.006734350297014738*((67.0820393249937*nuVtSqSum[2]+105.0*nuVtSqSum[0])*boundSurf[11]+93.91485505499116*nuVtSqSum[1]*boundSurf[4]+105.0*boundSurf[2]*nuVtSqSum[2]); 
  boundSurf_incr[12] = 0.04714045207910316*(13.41640786499874*nuVtSqSum[1]*boundSurf[20]+(13.41640786499874*nuVtSqSum[2]+15.0*nuVtSqSum[0])*boundSurf[12]+15.0*nuVtSqSum[1]*boundSurf[8]); 
  boundSurf_incr[13] = 0.006734350297014738*((67.0820393249937*nuVtSqSum[2]+105.0*nuVtSqSum[0])*boundSurf[13]+93.91485505499116*nuVtSqSum[1]*boundSurf[5]+105.0*nuVtSqSum[2]*boundSurf[3]); 
  boundSurf_incr[14] = 0.04714045207910316*(15.0*(nuVtSqSum[2]*boundSurf[23]+nuVtSqSum[1]*boundSurf[18])+15.0*nuVtSqSum[0]*boundSurf[14]); 
  boundSurf_incr[15] = 0.04714045207910316*(13.41640786499874*nuVtSqSum[1]*boundSurf[21]+(13.41640786499874*nuVtSqSum[2]+15.0*nuVtSqSum[0])*boundSurf[15]+15.0*nuVtSqSum[1]*boundSurf[9]); 
  boundSurf_incr[16] = 0.04714045207910316*(15.0*(nuVtSqSum[2]*boundSurf[24]+nuVtSqSum[1]*boundSurf[19])+15.0*nuVtSqSum[0]*boundSurf[16]); 
  boundSurf_incr[17] = 0.02020305089104421*((22.3606797749979*nuVtSqSum[2]+35.0*nuVtSqSum[0])*boundSurf[17]+31.30495168499706*nuVtSqSum[1]*boundSurf[10]+35.0*nuVtSqSum[2]*boundSurf[6]); 
  boundSurf_incr[18] = 0.04714045207910316*(13.41640786499874*nuVtSqSum[1]*boundSurf[23]+(13.41640786499874*nuVtSqSum[2]+15.0*nuVtSqSum[0])*boundSurf[18]+15.0*nuVtSqSum[1]*boundSurf[14]); 
  boundSurf_incr[19] = 0.04714045207910316*(13.41640786499874*nuVtSqSum[1]*boundSurf[24]+(13.41640786499874*nuVtSqSum[2]+15.0*nuVtSqSum[0])*boundSurf[19]+15.0*nuVtSqSum[1]*boundSurf[16]); 
  boundSurf_incr[20] = 0.006734350297014738*((67.0820393249937*nuVtSqSum[2]+105.0*nuVtSqSum[0])*boundSurf[20]+93.91485505499116*nuVtSqSum[1]*boundSurf[12]+105.0*nuVtSqSum[2]*boundSurf[8]); 
  boundSurf_incr[21] = 0.006734350297014738*((67.0820393249937*nuVtSqSum[2]+105.0*nuVtSqSum[0])*boundSurf[21]+93.91485505499116*nuVtSqSum[1]*boundSurf[15]+105.0*nuVtSqSum[2]*boundSurf[9]); 
  boundSurf_incr[22] = 0.7071067811865475*(nuVtSqSum[2]*boundSurf[26]+nuVtSqSum[1]*boundSurf[25]+nuVtSqSum[0]*boundSurf[22]); 
  boundSurf_incr[23] = 0.006734350297014738*((67.0820393249937*nuVtSqSum[2]+105.0*nuVtSqSum[0])*boundSurf[23]+93.91485505499116*nuVtSqSum[1]*boundSurf[18]+105.0*nuVtSqSum[2]*boundSurf[14]); 
  boundSurf_incr[24] = 0.006734350297014738*((67.0820393249937*nuVtSqSum[2]+105.0*nuVtSqSum[0])*boundSurf[24]+93.91485505499116*nuVtSqSum[1]*boundSurf[19]+105.0*nuVtSqSum[2]*boundSurf[16]); 
  boundSurf_incr[25] = 0.1414213562373095*(4.47213595499958*nuVtSqSum[1]*boundSurf[26]+(4.47213595499958*nuVtSqSum[2]+5.0*nuVtSqSum[0])*boundSurf[25]+5.0*nuVtSqSum[1]*boundSurf[22]); 
  boundSurf_incr[26] = 0.02020305089104421*((22.3606797749979*nuVtSqSum[2]+35.0*nuVtSqSum[0])*boundSurf[26]+31.30495168499706*nuVtSqSum[1]*boundSurf[25]+35.0*nuVtSqSum[2]*boundSurf[22]); 

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
