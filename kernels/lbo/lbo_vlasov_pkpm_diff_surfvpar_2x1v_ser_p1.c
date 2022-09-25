#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH void lbo_vlasov_pkpm_diff_surfvpar_2x1v_ser_p1(const double *w, const double *dxv, const double *nuVtSq, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[3]:         cell-center coordinates. 
  // dxv[3]:       cell spacing. 
  // nuVtSqSum[4]: sum of thermal speeds squared time their respective collisionalities. 
  // fl/fc/fr:      distribution function in cells 
  // out:           incremented distribution function in cell 
  double rdvSq4 = 4.0/(dxv[2]*dxv[2]); 
  double temp_diff[12] = {0.0}; 
  double diff_incr[12] = {0.0}; 

  temp_diff[0] = 0.6708203932499369*fr[8]+0.6708203932499369*fl[8]-1.341640786499874*fc[8]-1.190784930203603*fr[3]+1.190784930203603*fl[3]+0.9375*fr[0]+0.9375*fl[0]-1.875*fc[0]; 
  temp_diff[1] = 0.6708203932499369*fr[9]+0.6708203932499369*fl[9]-1.341640786499874*fc[9]-1.190784930203603*fr[5]+1.190784930203603*fl[5]+0.9375*fr[1]+0.9375*fl[1]-1.875*fc[1]; 
  temp_diff[2] = 0.6708203932499369*fr[10]+0.6708203932499369*fl[10]-1.341640786499874*fc[10]-1.190784930203603*fr[6]+1.190784930203603*fl[6]+0.9375*fr[2]+0.9375*fl[2]-1.875*fc[2]; 
  temp_diff[3] = 0.7382874503707888*fr[8]-0.7382874503707888*fl[8]-1.453125*fr[3]-1.453125*fl[3]-5.34375*fc[3]+1.190784930203603*fr[0]-1.190784930203603*fl[0]; 
  temp_diff[4] = 0.6708203932499369*fr[11]+0.6708203932499369*fl[11]-1.341640786499874*fc[11]-1.190784930203603*fr[7]+1.190784930203603*fl[7]+0.9375*fr[4]+0.9375*fl[4]-1.875*fc[4]; 
  temp_diff[5] = 0.7382874503707888*fr[9]-0.7382874503707888*fl[9]-1.453125*fr[5]-1.453125*fl[5]-5.34375*fc[5]+1.190784930203603*fr[1]-1.190784930203603*fl[1]; 
  temp_diff[6] = 0.7382874503707888*fr[10]-0.7382874503707888*fl[10]-1.453125*fr[6]-1.453125*fl[6]-5.34375*fc[6]+1.190784930203603*fr[2]-1.190784930203603*fl[2]; 
  temp_diff[7] = 0.7382874503707888*fr[11]-0.7382874503707888*fl[11]-1.453125*fr[7]-1.453125*fl[7]-5.34375*fc[7]+1.190784930203603*fr[4]-1.190784930203603*fl[4]; 
  temp_diff[8] = (-0.140625*fr[8])-0.140625*fl[8]-6.28125*fc[8]-0.3025768239224545*fr[3]+0.3025768239224545*fl[3]+0.4192627457812106*fr[0]+0.4192627457812106*fl[0]-0.8385254915624212*fc[0]; 
  temp_diff[9] = (-0.140625*fr[9])-0.140625*fl[9]-6.28125*fc[9]-0.3025768239224544*fr[5]+0.3025768239224544*fl[5]+0.4192627457812105*fr[1]+0.4192627457812105*fl[1]-0.8385254915624211*fc[1]; 
  temp_diff[10] = (-0.140625*fr[10])-0.140625*fl[10]-6.28125*fc[10]-0.3025768239224544*fr[6]+0.3025768239224544*fl[6]+0.4192627457812105*fr[2]+0.4192627457812105*fl[2]-0.8385254915624211*fc[2]; 
  temp_diff[11] = (-0.140625*fr[11])-0.140625*fl[11]-6.28125*fc[11]-0.3025768239224545*fr[7]+0.3025768239224545*fl[7]+0.4192627457812106*fr[4]+0.4192627457812106*fl[4]-0.8385254915624212*fc[4]; 

  diff_incr[0] = 0.5*nuVtSq[3]*temp_diff[4]+0.5*nuVtSq[2]*temp_diff[2]+0.5*nuVtSq[1]*temp_diff[1]+0.5*nuVtSq[0]*temp_diff[0]; 
  diff_incr[1] = 0.5*nuVtSq[2]*temp_diff[4]+0.5*temp_diff[2]*nuVtSq[3]+0.5*nuVtSq[0]*temp_diff[1]+0.5*temp_diff[0]*nuVtSq[1]; 
  diff_incr[2] = 0.5*nuVtSq[1]*temp_diff[4]+0.5*temp_diff[1]*nuVtSq[3]+0.5*nuVtSq[0]*temp_diff[2]+0.5*temp_diff[0]*nuVtSq[2]; 
  diff_incr[3] = 0.5*nuVtSq[3]*temp_diff[7]+0.5*nuVtSq[2]*temp_diff[6]+0.5*nuVtSq[1]*temp_diff[5]+0.5*nuVtSq[0]*temp_diff[3]; 
  diff_incr[4] = 0.5*nuVtSq[0]*temp_diff[4]+0.5*temp_diff[0]*nuVtSq[3]+0.5*nuVtSq[1]*temp_diff[2]+0.5*temp_diff[1]*nuVtSq[2]; 
  diff_incr[5] = 0.5*nuVtSq[2]*temp_diff[7]+0.5*nuVtSq[3]*temp_diff[6]+0.5*nuVtSq[0]*temp_diff[5]+0.5*nuVtSq[1]*temp_diff[3]; 
  diff_incr[6] = 0.5*nuVtSq[1]*temp_diff[7]+0.5*nuVtSq[0]*temp_diff[6]+0.5*nuVtSq[3]*temp_diff[5]+0.5*nuVtSq[2]*temp_diff[3]; 
  diff_incr[7] = 0.5*nuVtSq[0]*temp_diff[7]+0.5*nuVtSq[1]*temp_diff[6]+0.5*nuVtSq[2]*temp_diff[5]+0.5*nuVtSq[3]*temp_diff[3]; 
  diff_incr[8] = 0.5*nuVtSq[3]*temp_diff[11]+0.5000000000000001*nuVtSq[2]*temp_diff[10]+0.5000000000000001*nuVtSq[1]*temp_diff[9]+0.5*nuVtSq[0]*temp_diff[8]; 
  diff_incr[9] = 0.5000000000000001*nuVtSq[2]*temp_diff[11]+0.5*nuVtSq[3]*temp_diff[10]+0.5*nuVtSq[0]*temp_diff[9]+0.5000000000000001*nuVtSq[1]*temp_diff[8]; 
  diff_incr[10] = 0.5000000000000001*nuVtSq[1]*temp_diff[11]+0.5*nuVtSq[0]*temp_diff[10]+0.5*nuVtSq[3]*temp_diff[9]+0.5000000000000001*nuVtSq[2]*temp_diff[8]; 
  diff_incr[11] = 0.5*nuVtSq[0]*temp_diff[11]+0.5000000000000001*nuVtSq[1]*temp_diff[10]+0.5000000000000001*nuVtSq[2]*temp_diff[9]+0.5*nuVtSq[3]*temp_diff[8]; 

  out[0] += diff_incr[0]*rdvSq4; 
  out[1] += diff_incr[1]*rdvSq4; 
  out[2] += diff_incr[2]*rdvSq4; 
  out[3] += diff_incr[3]*rdvSq4; 
  out[4] += diff_incr[4]*rdvSq4; 
  out[5] += diff_incr[5]*rdvSq4; 
  out[6] += diff_incr[6]*rdvSq4; 
  out[7] += diff_incr[7]*rdvSq4; 
  out[8] += diff_incr[8]*rdvSq4; 
  out[9] += diff_incr[9]*rdvSq4; 
  out[10] += diff_incr[10]*rdvSq4; 
  out[11] += diff_incr[11]*rdvSq4; 
} 
