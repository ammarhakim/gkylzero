#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH void lbo_gyrokinetic_diff_surfvpar_1x1v_ser_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // m_:            species mass.
  // bmag_inv:      1/(magnetic field magnitude). 
  // w[2]:         cell-center coordinates. 
  // dxv[2]:       cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[3]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[3]: sum of thermal speeds squared time their respective collisionalities. 
  // fl/fc/fr:      distribution function in cells 
  // out:           incremented distribution function in cell 
  double rdvSq4 = 4.0/(dxv[1]*dxv[1]); 
  double temp_diff[8] = {0.0}; 
  double diff_incr[8] = {0.0}; 

  temp_diff[0] = 0.6708203932499369*fr[5]+0.6708203932499369*fl[5]-1.341640786499874*fc[5]-1.190784930203603*fr[2]+1.190784930203603*fl[2]+0.9375*fr[0]+0.9375*fl[0]-1.875*fc[0]; 
  temp_diff[1] = 0.6708203932499369*fr[7]+0.6708203932499369*fl[7]-1.341640786499874*fc[7]-1.190784930203603*fr[3]+1.190784930203603*fl[3]+0.9375*fr[1]+0.9375*fl[1]-1.875*fc[1]; 
  temp_diff[2] = 0.7382874503707888*fr[5]-0.7382874503707888*fl[5]-1.453125*fr[2]-1.453125*fl[2]-5.34375*fc[2]+1.190784930203603*fr[0]-1.190784930203603*fl[0]; 
  temp_diff[3] = 0.7382874503707888*fr[7]-0.7382874503707888*fl[7]-1.453125*fr[3]-1.453125*fl[3]-5.34375*fc[3]+1.190784930203603*fr[1]-1.190784930203603*fl[1]; 
  temp_diff[4] = (-1.190784930203603*fr[6])+1.190784930203603*fl[6]+0.9375*fr[4]+0.9375*fl[4]-1.875*fc[4]; 
  temp_diff[5] = (-0.140625*fr[5])-0.140625*fl[5]-6.28125*fc[5]-0.3025768239224545*fr[2]+0.3025768239224545*fl[2]+0.4192627457812106*fr[0]+0.4192627457812106*fl[0]-0.8385254915624212*fc[0]; 
  temp_diff[6] = (-1.453125*fr[6])-1.453125*fl[6]-5.34375*fc[6]+1.190784930203603*fr[4]-1.190784930203603*fl[4]; 
  temp_diff[7] = (-0.140625*fr[7])-0.140625*fl[7]-6.28125*fc[7]-0.3025768239224544*fr[3]+0.3025768239224544*fl[3]+0.4192627457812105*fr[1]+0.4192627457812105*fl[1]-0.8385254915624211*fc[1]; 

  diff_incr[0] = 0.7071067811865475*nuVtSqSum[2]*temp_diff[4]+0.7071067811865475*nuVtSqSum[1]*temp_diff[1]+0.7071067811865475*nuVtSqSum[0]*temp_diff[0]; 
  diff_incr[1] = 0.6324555320336759*nuVtSqSum[1]*temp_diff[4]+0.6324555320336759*temp_diff[1]*nuVtSqSum[2]+0.7071067811865475*nuVtSqSum[0]*temp_diff[1]+0.7071067811865475*temp_diff[0]*nuVtSqSum[1]; 
  diff_incr[2] = 0.7071067811865475*nuVtSqSum[2]*temp_diff[6]+0.7071067811865475*nuVtSqSum[1]*temp_diff[3]+0.7071067811865475*nuVtSqSum[0]*temp_diff[2]; 
  diff_incr[3] = 0.632455532033676*nuVtSqSum[1]*temp_diff[6]+0.6324555320336759*nuVtSqSum[2]*temp_diff[3]+0.7071067811865475*nuVtSqSum[0]*temp_diff[3]+0.7071067811865475*nuVtSqSum[1]*temp_diff[2]; 
  diff_incr[4] = 0.4517539514526256*nuVtSqSum[2]*temp_diff[4]+0.7071067811865475*nuVtSqSum[0]*temp_diff[4]+0.7071067811865475*temp_diff[0]*nuVtSqSum[2]+0.6324555320336759*nuVtSqSum[1]*temp_diff[1]; 
  diff_incr[5] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[7]+0.7071067811865475*nuVtSqSum[0]*temp_diff[5]; 
  diff_incr[6] = 0.4517539514526256*nuVtSqSum[2]*temp_diff[6]+0.7071067811865475*nuVtSqSum[0]*temp_diff[6]+0.632455532033676*nuVtSqSum[1]*temp_diff[3]+0.7071067811865475*nuVtSqSum[2]*temp_diff[2]; 
  diff_incr[7] = 0.6324555320336759*nuVtSqSum[2]*temp_diff[7]+0.7071067811865475*nuVtSqSum[0]*temp_diff[7]+0.7071067811865475*nuVtSqSum[1]*temp_diff[5]; 

  out[0] += diff_incr[0]*rdvSq4; 
  out[1] += diff_incr[1]*rdvSq4; 
  out[2] += diff_incr[2]*rdvSq4; 
  out[3] += diff_incr[3]*rdvSq4; 
  out[4] += diff_incr[4]*rdvSq4; 
  out[5] += diff_incr[5]*rdvSq4; 
  out[6] += diff_incr[6]*rdvSq4; 
  out[7] += diff_incr[7]*rdvSq4; 
} 
