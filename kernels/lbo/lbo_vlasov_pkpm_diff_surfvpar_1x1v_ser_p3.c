#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH void lbo_vlasov_pkpm_diff_surfvpar_1x1v_ser_p3(const double *w, const double *dxv, const double *nuVtSq, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[2]:         cell-center coordinates. 
  // dxv[2]:       cell spacing. 
  // nuVtSqSum[4]: sum of thermal speeds squared time their respective collisionalities. 
  // fl/fc/fr:      distribution function in cells 
  // out:           incremented distribution function in cell 
  double rdvSq4 = 4.0/(dxv[1]*dxv[1]); 
  double temp_diff[12] = {0.0}; 
  double diff_incr[12] = {0.0}; 

  temp_diff[0] = (-0.6821077598838398*fr[9])+0.6821077598838398*fl[9]+1.554766015605323*fr[5]+1.554766015605323*fl[5]-3.109532031210645*fc[5]-1.935025511580855*fr[2]+1.935025511580855*fl[2]+1.3671875*fr[0]+1.3671875*fl[0]-2.734375*fc[0]; 
  temp_diff[1] = (-0.6821077598838398*fr[11])+0.6821077598838398*fl[11]+1.554766015605323*fr[7]+1.554766015605323*fl[7]-3.109532031210645*fc[7]-1.935025511580855*fr[3]+1.935025511580855*fl[3]+1.3671875*fr[1]+1.3671875*fl[1]-2.734375*fc[1]; 
  temp_diff[2] = (-0.8541184610018139*fr[9])-0.8541184610018139*fl[9]-3.017544263419582*fc[9]+2.087780085064936*fr[5]-2.087780085064936*fl[5]-2.6953125*fr[2]-2.6953125*fl[2]-8.015625*fc[2]+1.935025511580855*fr[0]-1.935025511580855*fl[0]; 
  temp_diff[3] = (-0.8541184610018139*fr[11])-0.8541184610018139*fl[11]-3.017544263419582*fc[11]+2.087780085064936*fr[7]-2.087780085064936*fl[7]-2.6953125*fr[3]-2.6953125*fl[3]-8.015625*fc[3]+1.935025511580855*fr[1]-1.935025511580855*fl[1]; 
  temp_diff[4] = (-1.935025511580855*fr[6])+1.935025511580855*fl[6]+1.3671875*fr[4]+1.3671875*fl[4]-2.734375*fc[4]; 
  temp_diff[5] = (-0.2575079369875949*fr[9])+0.2575079369875949*fl[9]+1.1328125*fr[5]+1.1328125*fl[5]-11.640625*fc[5]-1.785203261142481*fr[2]+1.785203261142481*fl[2]+1.380073204863151*fr[0]+1.380073204863151*fl[0]-2.760146409726303*fc[0]; 
  temp_diff[6] = (-2.6953125*fr[6])-2.6953125*fl[6]-8.015625*fc[6]+1.935025511580855*fr[4]-1.935025511580855*fl[4]; 
  temp_diff[7] = (-0.2575079369875949*fr[11])+0.2575079369875949*fl[11]+1.1328125*fr[7]+1.1328125*fl[7]-11.640625*fc[7]-1.785203261142481*fr[3]+1.785203261142481*fl[3]+1.380073204863152*fr[1]+1.380073204863152*fl[1]-2.760146409726303*fc[1]; 
  temp_diff[8] = (-1.935025511580855*fr[10])+1.935025511580855*fl[10]+1.3671875*fr[8]+1.3671875*fl[8]-2.734375*fc[8]; 
  temp_diff[9] = 1.1953125*fr[9]+1.1953125*fl[9]-9.609375*fc[9]-1.432800572469438*fr[5]+1.432800572469438*fl[5]+0.8950343154210625*fr[2]+0.8950343154210625*fl[2]+0.6444247071031648*fc[2]-0.351388846000766*fr[0]+0.351388846000766*fl[0]; 
  temp_diff[10] = (-2.6953125*fr[10])-2.6953125*fl[10]-8.015625*fc[10]+1.935025511580855*fr[8]-1.935025511580855*fl[8]; 
  temp_diff[11] = 1.1953125*fr[11]+1.1953125*fl[11]-9.609375*fc[11]-1.432800572469438*fr[7]+1.432800572469438*fl[7]+0.8950343154210625*fr[3]+0.8950343154210625*fl[3]+0.644424707103165*fc[3]-0.351388846000766*fr[1]+0.351388846000766*fl[1]; 

  diff_incr[0] = 0.7071067811865475*nuVtSq[3]*temp_diff[8]+0.7071067811865475*nuVtSq[2]*temp_diff[4]+0.7071067811865475*nuVtSq[1]*temp_diff[1]+0.7071067811865475*nuVtSq[0]*temp_diff[0]; 
  diff_incr[1] = 0.6210590034081186*nuVtSq[2]*temp_diff[8]+0.6210590034081186*nuVtSq[3]*temp_diff[4]+0.6324555320336759*nuVtSq[1]*temp_diff[4]+0.6324555320336759*temp_diff[1]*nuVtSq[2]+0.7071067811865475*nuVtSq[0]*temp_diff[1]+0.7071067811865475*temp_diff[0]*nuVtSq[1]; 
  diff_incr[2] = 0.7071067811865474*nuVtSq[3]*temp_diff[10]+0.7071067811865475*nuVtSq[2]*temp_diff[6]+0.7071067811865475*nuVtSq[1]*temp_diff[3]+0.7071067811865475*nuVtSq[0]*temp_diff[2]; 
  diff_incr[3] = 0.6210590034081187*nuVtSq[2]*temp_diff[10]+0.6210590034081187*nuVtSq[3]*temp_diff[6]+0.632455532033676*nuVtSq[1]*temp_diff[6]+0.6324555320336759*nuVtSq[2]*temp_diff[3]+0.7071067811865475*nuVtSq[0]*temp_diff[3]+0.7071067811865475*nuVtSq[1]*temp_diff[2]; 
  diff_incr[4] = 0.421637021355784*nuVtSq[3]*temp_diff[8]+0.6210590034081186*nuVtSq[1]*temp_diff[8]+0.4517539514526256*nuVtSq[2]*temp_diff[4]+0.7071067811865475*nuVtSq[0]*temp_diff[4]+0.6210590034081186*temp_diff[1]*nuVtSq[3]+0.7071067811865475*temp_diff[0]*nuVtSq[2]+0.6324555320336759*nuVtSq[1]*temp_diff[1]; 
  diff_incr[5] = 0.7071067811865475*nuVtSq[1]*temp_diff[7]+0.7071067811865475*nuVtSq[0]*temp_diff[5]; 
  diff_incr[6] = 0.4216370213557839*nuVtSq[3]*temp_diff[10]+0.6210590034081187*nuVtSq[1]*temp_diff[10]+0.4517539514526256*nuVtSq[2]*temp_diff[6]+0.7071067811865475*nuVtSq[0]*temp_diff[6]+0.6210590034081187*nuVtSq[3]*temp_diff[3]+0.632455532033676*nuVtSq[1]*temp_diff[3]+0.7071067811865475*nuVtSq[2]*temp_diff[2]; 
  diff_incr[7] = 0.6324555320336759*nuVtSq[2]*temp_diff[7]+0.7071067811865475*nuVtSq[0]*temp_diff[7]+0.7071067811865475*nuVtSq[1]*temp_diff[5]; 
  diff_incr[8] = 0.421637021355784*nuVtSq[2]*temp_diff[8]+0.7071067811865475*nuVtSq[0]*temp_diff[8]+0.421637021355784*nuVtSq[3]*temp_diff[4]+0.6210590034081186*nuVtSq[1]*temp_diff[4]+0.7071067811865475*temp_diff[0]*nuVtSq[3]+0.6210590034081186*temp_diff[1]*nuVtSq[2]; 
  diff_incr[9] = 0.7071067811865474*nuVtSq[1]*temp_diff[11]+0.7071067811865475*nuVtSq[0]*temp_diff[9]; 
  diff_incr[10] = 0.421637021355784*nuVtSq[2]*temp_diff[10]+0.7071067811865475*nuVtSq[0]*temp_diff[10]+0.4216370213557839*nuVtSq[3]*temp_diff[6]+0.6210590034081187*nuVtSq[1]*temp_diff[6]+0.6210590034081187*nuVtSq[2]*temp_diff[3]+0.7071067811865474*temp_diff[2]*nuVtSq[3]; 
  diff_incr[11] = 0.6324555320336759*nuVtSq[2]*temp_diff[11]+0.7071067811865475*nuVtSq[0]*temp_diff[11]+0.7071067811865474*nuVtSq[1]*temp_diff[9]; 

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
