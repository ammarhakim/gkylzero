#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH double lbo_gyrokinetic_diff_surfvpar_1x2v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // m_: species mass.
  // bmag_inv: 1/(magnetic field magnitude). 
  // w[3]: cell-center coordinates. 
  // dxv[3]: cell spacing. 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[4]: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // fl/fc/fr: distribution function in cells 
  // out: incremented distribution function in cell 

  const double *nuVtSqSum = &nuPrimMomsSum[2];

  double rdvSq4 = 4.0/(dxv[1]*dxv[1]); 
  double incr[12] = {0.0}; 

  double f_xx[12] = {0.0}; 
  f_xx[0] = 0.6708203932499369*fr[8]+0.6708203932499369*fl[8]-1.341640786499874*fc[8]-1.190784930203603*fr[2]+1.190784930203603*fl[2]+0.9375*fr[0]+0.9375*fl[0]-1.875*fc[0]; 
  f_xx[1] = 0.6708203932499369*fr[9]+0.6708203932499369*fl[9]-1.341640786499874*fc[9]-1.190784930203603*fr[4]+1.190784930203603*fl[4]+0.9375*fr[1]+0.9375*fl[1]-1.875*fc[1]; 
  f_xx[2] = 0.7382874503707888*fr[8]-0.7382874503707888*fl[8]-1.453125*fr[2]-1.453125*fl[2]-5.34375*fc[2]+1.190784930203603*fr[0]-1.190784930203603*fl[0]; 
  f_xx[3] = 0.6708203932499369*fr[10]+0.6708203932499369*fl[10]-1.341640786499874*fc[10]-1.190784930203603*fr[6]+1.190784930203603*fl[6]+0.9375*fr[3]+0.9375*fl[3]-1.875*fc[3]; 
  f_xx[4] = 0.7382874503707888*fr[9]-0.7382874503707888*fl[9]-1.453125*fr[4]-1.453125*fl[4]-5.34375*fc[4]+1.190784930203603*fr[1]-1.190784930203603*fl[1]; 
  f_xx[5] = 0.6708203932499369*fr[11]+0.6708203932499369*fl[11]-1.341640786499874*fc[11]-1.190784930203603*fr[7]+1.190784930203603*fl[7]+0.9375*fr[5]+0.9375*fl[5]-1.875*fc[5]; 
  f_xx[6] = 0.7382874503707888*fr[10]-0.7382874503707888*fl[10]-1.453125*fr[6]-1.453125*fl[6]-5.34375*fc[6]+1.190784930203603*fr[3]-1.190784930203603*fl[3]; 
  f_xx[7] = 0.7382874503707888*fr[11]-0.7382874503707888*fl[11]-1.453125*fr[7]-1.453125*fl[7]-5.34375*fc[7]+1.190784930203603*fr[5]-1.190784930203603*fl[5]; 
  f_xx[8] = (-0.140625*fr[8])-0.140625*fl[8]-6.28125*fc[8]-0.3025768239224545*fr[2]+0.3025768239224545*fl[2]+0.4192627457812106*fr[0]+0.4192627457812106*fl[0]-0.8385254915624212*fc[0]; 
  f_xx[9] = (-0.140625*fr[9])-0.140625*fl[9]-6.28125*fc[9]-0.3025768239224544*fr[4]+0.3025768239224544*fl[4]+0.4192627457812105*fr[1]+0.4192627457812105*fl[1]-0.8385254915624211*fc[1]; 
  f_xx[10] = (-0.140625*fr[10])-0.140625*fl[10]-6.28125*fc[10]-0.3025768239224544*fr[6]+0.3025768239224544*fl[6]+0.4192627457812105*fr[3]+0.4192627457812105*fl[3]-0.8385254915624211*fc[3]; 
  f_xx[11] = (-0.140625*fr[11])-0.140625*fl[11]-6.28125*fc[11]-0.3025768239224545*fr[7]+0.3025768239224545*fl[7]+0.4192627457812106*fr[5]+0.4192627457812106*fl[5]-0.8385254915624212*fc[5]; 

  incr[0] = 0.7071067811865475*f_xx[1]*nuVtSqSum[1]+0.7071067811865475*f_xx[0]*nuVtSqSum[0]; 
  incr[1] = 0.7071067811865475*f_xx[0]*nuVtSqSum[1]+0.7071067811865475*nuVtSqSum[0]*f_xx[1]; 
  incr[2] = 0.7071067811865475*nuVtSqSum[1]*f_xx[4]+0.7071067811865475*nuVtSqSum[0]*f_xx[2]; 
  incr[3] = 0.7071067811865475*nuVtSqSum[1]*f_xx[5]+0.7071067811865475*nuVtSqSum[0]*f_xx[3]; 
  incr[4] = 0.7071067811865475*nuVtSqSum[0]*f_xx[4]+0.7071067811865475*nuVtSqSum[1]*f_xx[2]; 
  incr[5] = 0.7071067811865475*nuVtSqSum[0]*f_xx[5]+0.7071067811865475*nuVtSqSum[1]*f_xx[3]; 
  incr[6] = 0.7071067811865475*nuVtSqSum[1]*f_xx[7]+0.7071067811865475*nuVtSqSum[0]*f_xx[6]; 
  incr[7] = 0.7071067811865475*nuVtSqSum[0]*f_xx[7]+0.7071067811865475*nuVtSqSum[1]*f_xx[6]; 
  incr[8] = 0.7071067811865475*nuVtSqSum[1]*f_xx[9]+0.7071067811865475*nuVtSqSum[0]*f_xx[8]; 
  incr[9] = 0.7071067811865475*nuVtSqSum[0]*f_xx[9]+0.7071067811865475*nuVtSqSum[1]*f_xx[8]; 
  incr[10] = 0.7071067811865475*nuVtSqSum[1]*f_xx[11]+0.7071067811865475*nuVtSqSum[0]*f_xx[10]; 
  incr[11] = 0.7071067811865475*nuVtSqSum[0]*f_xx[11]+0.7071067811865475*nuVtSqSum[1]*f_xx[10]; 

  out[0] += incr[0]*rdvSq4; 
  out[1] += incr[1]*rdvSq4; 
  out[2] += incr[2]*rdvSq4; 
  out[3] += incr[3]*rdvSq4; 
  out[4] += incr[4]*rdvSq4; 
  out[5] += incr[5]*rdvSq4; 
  out[6] += incr[6]*rdvSq4; 
  out[7] += incr[7]*rdvSq4; 
  out[8] += incr[8]*rdvSq4; 
  out[9] += incr[9]*rdvSq4; 
  out[10] += incr[10]*rdvSq4; 
  out[11] += incr[11]*rdvSq4; 

  return 0.;

} 
