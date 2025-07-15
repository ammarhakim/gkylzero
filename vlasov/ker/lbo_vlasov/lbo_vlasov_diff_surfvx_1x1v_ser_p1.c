#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH double lbo_vlasov_diff_surfvx_1x1v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[2]: cell-center coordinates. 
  // dxv[2]: cell spacing. 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[4]: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // fl/fc/fr: distribution function in cells 
  // out: incremented distribution function in cell 

  const double *nuVtSqSum = &nuPrimMomsSum[2];

  double rdvSq4 = 4.0/(dxv[1]*dxv[1]); 
  double incr[6]; 

  double f_xx[6] = {0.0}; 
  f_xx[0] = 0.6708203932499369*fr[4]+0.6708203932499369*fl[4]-1.341640786499874*fc[4]-1.190784930203603*fr[2]+1.190784930203603*fl[2]+0.9375*fr[0]+0.9375*fl[0]-1.875*fc[0]; 
  f_xx[1] = 0.6708203932499369*fr[5]+0.6708203932499369*fl[5]-1.341640786499874*fc[5]-1.190784930203603*fr[3]+1.190784930203603*fl[3]+0.9375*fr[1]+0.9375*fl[1]-1.875*fc[1]; 
  f_xx[2] = 0.7382874503707888*fr[4]-0.7382874503707888*fl[4]-1.453125*fr[2]-1.453125*fl[2]-5.34375*fc[2]+1.190784930203603*fr[0]-1.190784930203603*fl[0]; 
  f_xx[3] = 0.7382874503707888*fr[5]-0.7382874503707888*fl[5]-1.453125*fr[3]-1.453125*fl[3]-5.34375*fc[3]+1.190784930203603*fr[1]-1.190784930203603*fl[1]; 
  f_xx[4] = (-0.140625*fr[4])-0.140625*fl[4]-6.28125*fc[4]-0.3025768239224545*fr[2]+0.3025768239224545*fl[2]+0.4192627457812106*fr[0]+0.4192627457812106*fl[0]-0.8385254915624212*fc[0]; 
  f_xx[5] = (-0.140625*fr[5])-0.140625*fl[5]-6.28125*fc[5]-0.3025768239224544*fr[3]+0.3025768239224544*fl[3]+0.4192627457812105*fr[1]+0.4192627457812105*fl[1]-0.8385254915624211*fc[1]; 

  incr[0] = 0.7071067811865475*f_xx[1]*nuVtSqSum[1]+0.7071067811865475*f_xx[0]*nuVtSqSum[0]; 
  incr[1] = 0.7071067811865475*f_xx[0]*nuVtSqSum[1]+0.7071067811865475*nuVtSqSum[0]*f_xx[1]; 
  incr[2] = 0.7071067811865475*nuVtSqSum[1]*f_xx[3]+0.7071067811865475*nuVtSqSum[0]*f_xx[2]; 
  incr[3] = 0.7071067811865475*nuVtSqSum[0]*f_xx[3]+0.7071067811865475*nuVtSqSum[1]*f_xx[2]; 
  incr[4] = 0.7071067811865475*nuVtSqSum[1]*f_xx[5]+0.7071067811865475*nuVtSqSum[0]*f_xx[4]; 
  incr[5] = 0.7071067811865475*nuVtSqSum[0]*f_xx[5]+0.7071067811865475*nuVtSqSum[1]*f_xx[4]; 

  out[0] += incr[0]*rdvSq4; 
  out[1] += incr[1]*rdvSq4; 
  out[2] += incr[2]*rdvSq4; 
  out[3] += incr[3]*rdvSq4; 
  out[4] += incr[4]*rdvSq4; 
  out[5] += incr[5]*rdvSq4; 

  return 0.;

} 
