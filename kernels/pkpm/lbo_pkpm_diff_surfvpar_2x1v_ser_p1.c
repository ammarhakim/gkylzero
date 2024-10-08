#include <gkyl_lbo_pkpm_kernels.h> 
GKYL_CU_DH double lbo_pkpm_diff_surfvpar_2x1v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:       Cell-center coordinates. 
  // dxv[NDIM]:     Cell spacing. 
  // nuSum:         Collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum: Sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 
  // fl/fc/fr:      Input Distribution function [F_0, T_perp G = T_perp (F_1 - F_0)] in left/center/right cells 
  // out:           Incremented distribution functions in center cell. 
  double rdvSq4 = 4.0/(dxv[2]*dxv[2]); 
  const double *nuVtSqSum = &nuPrimMomsSum[4];

  const double *F_0l = &fl[0]; 
  const double *G_1l = &fl[12]; 
  const double *F_0c = &fc[0]; 
  const double *G_1c = &fc[12]; 
  const double *F_0r = &fr[0]; 
  const double *G_1r = &fr[12]; 
  double *out_F_0 = &out[0]; 
  double *out_G_1 = &out[12]; 
  double incr_F_0[12] = {0.0}; 
  double incr_G_1[12] = {0.0}; 

  double F_0_xx[12] = {0.0}; 
  double G_1_xx[12] = {0.0}; 
  F_0_xx[0] = 0.6708203932499369*F_0r[8]+0.6708203932499369*F_0l[8]-1.341640786499874*F_0c[8]-1.190784930203603*F_0r[3]+1.190784930203603*F_0l[3]+0.9375*F_0r[0]+0.9375*F_0l[0]-1.875*F_0c[0]; 
  F_0_xx[1] = 0.6708203932499369*F_0r[9]+0.6708203932499369*F_0l[9]-1.341640786499874*F_0c[9]-1.190784930203603*F_0r[5]+1.190784930203603*F_0l[5]+0.9375*F_0r[1]+0.9375*F_0l[1]-1.875*F_0c[1]; 
  F_0_xx[2] = 0.6708203932499369*F_0r[10]+0.6708203932499369*F_0l[10]-1.341640786499874*F_0c[10]-1.190784930203603*F_0r[6]+1.190784930203603*F_0l[6]+0.9375*F_0r[2]+0.9375*F_0l[2]-1.875*F_0c[2]; 
  F_0_xx[3] = 0.7382874503707888*F_0r[8]-0.7382874503707888*F_0l[8]-1.453125*F_0r[3]-1.453125*F_0l[3]-5.34375*F_0c[3]+1.190784930203603*F_0r[0]-1.190784930203603*F_0l[0]; 
  F_0_xx[4] = 0.6708203932499369*F_0r[11]+0.6708203932499369*F_0l[11]-1.341640786499874*F_0c[11]-1.190784930203603*F_0r[7]+1.190784930203603*F_0l[7]+0.9375*F_0r[4]+0.9375*F_0l[4]-1.875*F_0c[4]; 
  F_0_xx[5] = 0.7382874503707888*F_0r[9]-0.7382874503707888*F_0l[9]-1.453125*F_0r[5]-1.453125*F_0l[5]-5.34375*F_0c[5]+1.190784930203603*F_0r[1]-1.190784930203603*F_0l[1]; 
  F_0_xx[6] = 0.7382874503707888*F_0r[10]-0.7382874503707888*F_0l[10]-1.453125*F_0r[6]-1.453125*F_0l[6]-5.34375*F_0c[6]+1.190784930203603*F_0r[2]-1.190784930203603*F_0l[2]; 
  F_0_xx[7] = 0.7382874503707888*F_0r[11]-0.7382874503707888*F_0l[11]-1.453125*F_0r[7]-1.453125*F_0l[7]-5.34375*F_0c[7]+1.190784930203603*F_0r[4]-1.190784930203603*F_0l[4]; 
  F_0_xx[8] = (-0.140625*F_0r[8])-0.140625*F_0l[8]-6.28125*F_0c[8]-0.3025768239224545*F_0r[3]+0.3025768239224545*F_0l[3]+0.4192627457812106*F_0r[0]+0.4192627457812106*F_0l[0]-0.8385254915624212*F_0c[0]; 
  F_0_xx[9] = (-0.140625*F_0r[9])-0.140625*F_0l[9]-6.28125*F_0c[9]-0.3025768239224544*F_0r[5]+0.3025768239224544*F_0l[5]+0.4192627457812105*F_0r[1]+0.4192627457812105*F_0l[1]-0.8385254915624211*F_0c[1]; 
  F_0_xx[10] = (-0.140625*F_0r[10])-0.140625*F_0l[10]-6.28125*F_0c[10]-0.3025768239224544*F_0r[6]+0.3025768239224544*F_0l[6]+0.4192627457812105*F_0r[2]+0.4192627457812105*F_0l[2]-0.8385254915624211*F_0c[2]; 
  F_0_xx[11] = (-0.140625*F_0r[11])-0.140625*F_0l[11]-6.28125*F_0c[11]-0.3025768239224545*F_0r[7]+0.3025768239224545*F_0l[7]+0.4192627457812106*F_0r[4]+0.4192627457812106*F_0l[4]-0.8385254915624212*F_0c[4]; 
  G_1_xx[0] = 0.6708203932499369*G_1r[8]+0.6708203932499369*G_1l[8]-1.341640786499874*G_1c[8]-1.190784930203603*G_1r[3]+1.190784930203603*G_1l[3]+0.9375*G_1r[0]+0.9375*G_1l[0]-1.875*G_1c[0]; 
  G_1_xx[1] = 0.6708203932499369*G_1r[9]+0.6708203932499369*G_1l[9]-1.341640786499874*G_1c[9]-1.190784930203603*G_1r[5]+1.190784930203603*G_1l[5]+0.9375*G_1r[1]+0.9375*G_1l[1]-1.875*G_1c[1]; 
  G_1_xx[2] = 0.6708203932499369*G_1r[10]+0.6708203932499369*G_1l[10]-1.341640786499874*G_1c[10]-1.190784930203603*G_1r[6]+1.190784930203603*G_1l[6]+0.9375*G_1r[2]+0.9375*G_1l[2]-1.875*G_1c[2]; 
  G_1_xx[3] = 0.7382874503707888*G_1r[8]-0.7382874503707888*G_1l[8]-1.453125*G_1r[3]-1.453125*G_1l[3]-5.34375*G_1c[3]+1.190784930203603*G_1r[0]-1.190784930203603*G_1l[0]; 
  G_1_xx[4] = 0.6708203932499369*G_1r[11]+0.6708203932499369*G_1l[11]-1.341640786499874*G_1c[11]-1.190784930203603*G_1r[7]+1.190784930203603*G_1l[7]+0.9375*G_1r[4]+0.9375*G_1l[4]-1.875*G_1c[4]; 
  G_1_xx[5] = 0.7382874503707888*G_1r[9]-0.7382874503707888*G_1l[9]-1.453125*G_1r[5]-1.453125*G_1l[5]-5.34375*G_1c[5]+1.190784930203603*G_1r[1]-1.190784930203603*G_1l[1]; 
  G_1_xx[6] = 0.7382874503707888*G_1r[10]-0.7382874503707888*G_1l[10]-1.453125*G_1r[6]-1.453125*G_1l[6]-5.34375*G_1c[6]+1.190784930203603*G_1r[2]-1.190784930203603*G_1l[2]; 
  G_1_xx[7] = 0.7382874503707888*G_1r[11]-0.7382874503707888*G_1l[11]-1.453125*G_1r[7]-1.453125*G_1l[7]-5.34375*G_1c[7]+1.190784930203603*G_1r[4]-1.190784930203603*G_1l[4]; 
  G_1_xx[8] = (-0.140625*G_1r[8])-0.140625*G_1l[8]-6.28125*G_1c[8]-0.3025768239224545*G_1r[3]+0.3025768239224545*G_1l[3]+0.4192627457812106*G_1r[0]+0.4192627457812106*G_1l[0]-0.8385254915624212*G_1c[0]; 
  G_1_xx[9] = (-0.140625*G_1r[9])-0.140625*G_1l[9]-6.28125*G_1c[9]-0.3025768239224544*G_1r[5]+0.3025768239224544*G_1l[5]+0.4192627457812105*G_1r[1]+0.4192627457812105*G_1l[1]-0.8385254915624211*G_1c[1]; 
  G_1_xx[10] = (-0.140625*G_1r[10])-0.140625*G_1l[10]-6.28125*G_1c[10]-0.3025768239224544*G_1r[6]+0.3025768239224544*G_1l[6]+0.4192627457812105*G_1r[2]+0.4192627457812105*G_1l[2]-0.8385254915624211*G_1c[2]; 
  G_1_xx[11] = (-0.140625*G_1r[11])-0.140625*G_1l[11]-6.28125*G_1c[11]-0.3025768239224545*G_1r[7]+0.3025768239224545*G_1l[7]+0.4192627457812106*G_1r[4]+0.4192627457812106*G_1l[4]-0.8385254915624212*G_1c[4]; 

  incr_F_0[0] = 0.5*nuVtSqSum[3]*F_0_xx[4]+0.5*F_0_xx[2]*nuVtSqSum[2]+0.5*F_0_xx[1]*nuVtSqSum[1]+0.5*F_0_xx[0]*nuVtSqSum[0]; 
  incr_F_0[1] = 0.5*nuVtSqSum[2]*F_0_xx[4]+0.5*F_0_xx[2]*nuVtSqSum[3]+0.5*F_0_xx[0]*nuVtSqSum[1]+0.5*nuVtSqSum[0]*F_0_xx[1]; 
  incr_F_0[2] = 0.5*nuVtSqSum[1]*F_0_xx[4]+0.5*F_0_xx[1]*nuVtSqSum[3]+0.5*F_0_xx[0]*nuVtSqSum[2]+0.5*nuVtSqSum[0]*F_0_xx[2]; 
  incr_F_0[3] = 0.5*nuVtSqSum[3]*F_0_xx[7]+0.5*nuVtSqSum[2]*F_0_xx[6]+0.5*nuVtSqSum[1]*F_0_xx[5]+0.5*nuVtSqSum[0]*F_0_xx[3]; 
  incr_F_0[4] = 0.5*nuVtSqSum[0]*F_0_xx[4]+0.5*F_0_xx[0]*nuVtSqSum[3]+0.5*F_0_xx[1]*nuVtSqSum[2]+0.5*nuVtSqSum[1]*F_0_xx[2]; 
  incr_F_0[5] = 0.5*nuVtSqSum[2]*F_0_xx[7]+0.5*nuVtSqSum[3]*F_0_xx[6]+0.5*nuVtSqSum[0]*F_0_xx[5]+0.5*nuVtSqSum[1]*F_0_xx[3]; 
  incr_F_0[6] = 0.5*nuVtSqSum[1]*F_0_xx[7]+0.5*nuVtSqSum[0]*F_0_xx[6]+0.5*nuVtSqSum[3]*F_0_xx[5]+0.5*nuVtSqSum[2]*F_0_xx[3]; 
  incr_F_0[7] = 0.5*nuVtSqSum[0]*F_0_xx[7]+0.5*nuVtSqSum[1]*F_0_xx[6]+0.5*nuVtSqSum[2]*F_0_xx[5]+0.5*F_0_xx[3]*nuVtSqSum[3]; 
  incr_F_0[8] = 0.5*nuVtSqSum[3]*F_0_xx[11]+0.5000000000000001*nuVtSqSum[2]*F_0_xx[10]+0.5000000000000001*nuVtSqSum[1]*F_0_xx[9]+0.5*nuVtSqSum[0]*F_0_xx[8]; 
  incr_F_0[9] = 0.5000000000000001*nuVtSqSum[2]*F_0_xx[11]+0.5*nuVtSqSum[3]*F_0_xx[10]+0.5*nuVtSqSum[0]*F_0_xx[9]+0.5000000000000001*nuVtSqSum[1]*F_0_xx[8]; 
  incr_F_0[10] = 0.5000000000000001*nuVtSqSum[1]*F_0_xx[11]+0.5*nuVtSqSum[0]*F_0_xx[10]+0.5*nuVtSqSum[3]*F_0_xx[9]+0.5000000000000001*nuVtSqSum[2]*F_0_xx[8]; 
  incr_F_0[11] = 0.5*nuVtSqSum[0]*F_0_xx[11]+0.5000000000000001*nuVtSqSum[1]*F_0_xx[10]+0.5000000000000001*nuVtSqSum[2]*F_0_xx[9]+0.5*nuVtSqSum[3]*F_0_xx[8]; 
  incr_G_1[0] = 0.5*nuVtSqSum[3]*G_1_xx[4]+0.5*G_1_xx[2]*nuVtSqSum[2]+0.5*G_1_xx[1]*nuVtSqSum[1]+0.5*G_1_xx[0]*nuVtSqSum[0]; 
  incr_G_1[1] = 0.5*nuVtSqSum[2]*G_1_xx[4]+0.5*G_1_xx[2]*nuVtSqSum[3]+0.5*G_1_xx[0]*nuVtSqSum[1]+0.5*nuVtSqSum[0]*G_1_xx[1]; 
  incr_G_1[2] = 0.5*nuVtSqSum[1]*G_1_xx[4]+0.5*G_1_xx[1]*nuVtSqSum[3]+0.5*G_1_xx[0]*nuVtSqSum[2]+0.5*nuVtSqSum[0]*G_1_xx[2]; 
  incr_G_1[3] = 0.5*nuVtSqSum[3]*G_1_xx[7]+0.5*nuVtSqSum[2]*G_1_xx[6]+0.5*nuVtSqSum[1]*G_1_xx[5]+0.5*nuVtSqSum[0]*G_1_xx[3]; 
  incr_G_1[4] = 0.5*nuVtSqSum[0]*G_1_xx[4]+0.5*G_1_xx[0]*nuVtSqSum[3]+0.5*G_1_xx[1]*nuVtSqSum[2]+0.5*nuVtSqSum[1]*G_1_xx[2]; 
  incr_G_1[5] = 0.5*nuVtSqSum[2]*G_1_xx[7]+0.5*nuVtSqSum[3]*G_1_xx[6]+0.5*nuVtSqSum[0]*G_1_xx[5]+0.5*nuVtSqSum[1]*G_1_xx[3]; 
  incr_G_1[6] = 0.5*nuVtSqSum[1]*G_1_xx[7]+0.5*nuVtSqSum[0]*G_1_xx[6]+0.5*nuVtSqSum[3]*G_1_xx[5]+0.5*nuVtSqSum[2]*G_1_xx[3]; 
  incr_G_1[7] = 0.5*nuVtSqSum[0]*G_1_xx[7]+0.5*nuVtSqSum[1]*G_1_xx[6]+0.5*nuVtSqSum[2]*G_1_xx[5]+0.5*G_1_xx[3]*nuVtSqSum[3]; 
  incr_G_1[8] = 0.5*nuVtSqSum[3]*G_1_xx[11]+0.5000000000000001*nuVtSqSum[2]*G_1_xx[10]+0.5000000000000001*nuVtSqSum[1]*G_1_xx[9]+0.5*nuVtSqSum[0]*G_1_xx[8]; 
  incr_G_1[9] = 0.5000000000000001*nuVtSqSum[2]*G_1_xx[11]+0.5*nuVtSqSum[3]*G_1_xx[10]+0.5*nuVtSqSum[0]*G_1_xx[9]+0.5000000000000001*nuVtSqSum[1]*G_1_xx[8]; 
  incr_G_1[10] = 0.5000000000000001*nuVtSqSum[1]*G_1_xx[11]+0.5*nuVtSqSum[0]*G_1_xx[10]+0.5*nuVtSqSum[3]*G_1_xx[9]+0.5000000000000001*nuVtSqSum[2]*G_1_xx[8]; 
  incr_G_1[11] = 0.5*nuVtSqSum[0]*G_1_xx[11]+0.5000000000000001*nuVtSqSum[1]*G_1_xx[10]+0.5000000000000001*nuVtSqSum[2]*G_1_xx[9]+0.5*nuVtSqSum[3]*G_1_xx[8]; 

  out_F_0[0] += incr_F_0[0]*rdvSq4; 
  out_F_0[1] += incr_F_0[1]*rdvSq4; 
  out_F_0[2] += incr_F_0[2]*rdvSq4; 
  out_F_0[3] += incr_F_0[3]*rdvSq4; 
  out_F_0[4] += incr_F_0[4]*rdvSq4; 
  out_F_0[5] += incr_F_0[5]*rdvSq4; 
  out_F_0[6] += incr_F_0[6]*rdvSq4; 
  out_F_0[7] += incr_F_0[7]*rdvSq4; 
  out_F_0[8] += incr_F_0[8]*rdvSq4; 
  out_F_0[9] += incr_F_0[9]*rdvSq4; 
  out_F_0[10] += incr_F_0[10]*rdvSq4; 
  out_F_0[11] += incr_F_0[11]*rdvSq4; 
  out_G_1[0] += incr_G_1[0]*rdvSq4; 
  out_G_1[1] += incr_G_1[1]*rdvSq4; 
  out_G_1[2] += incr_G_1[2]*rdvSq4; 
  out_G_1[3] += incr_G_1[3]*rdvSq4; 
  out_G_1[4] += incr_G_1[4]*rdvSq4; 
  out_G_1[5] += incr_G_1[5]*rdvSq4; 
  out_G_1[6] += incr_G_1[6]*rdvSq4; 
  out_G_1[7] += incr_G_1[7]*rdvSq4; 
  out_G_1[8] += incr_G_1[8]*rdvSq4; 
  out_G_1[9] += incr_G_1[9]*rdvSq4; 
  out_G_1[10] += incr_G_1[10]*rdvSq4; 
  out_G_1[11] += incr_G_1[11]*rdvSq4; 

  return 0.;

} 
