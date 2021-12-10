#include <gkyl_vlasov_lbo_kernels.h> 
GKYL_CU_DH void vlasov_lbo_diff_surfvx_1x2v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[3]:         cell-center coordinates. 
  // dxv[3]:       cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[4]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[2]: sum of thermal speeds squared time their respective collisionalities. 
  // fl/fc/fr:      distribution function in cells 
  // out:           incremented distribution function in cell 
  double rdvSq4 = 4.0/(dxv[1]*dxv[1]); 
  double diff_incr[8] = {0.0}; 

  diff_incr[0] = (-0.3827327723098713*nuVtSqSum[1]*fr[4])+0.3827327723098713*nuVtSqSum[1]*fl[4]-0.3827327723098713*nuVtSqSum[0]*fr[2]+0.3827327723098713*nuVtSqSum[0]*fl[2]+0.3977475644174328*fr[1]*nuVtSqSum[1]+0.3977475644174328*fl[1]*nuVtSqSum[1]-0.7954951288348656*fc[1]*nuVtSqSum[1]+0.3977475644174328*fr[0]*nuVtSqSum[0]+0.3977475644174328*fl[0]*nuVtSqSum[0]-0.7954951288348656*fc[0]*nuVtSqSum[0]; 
  diff_incr[1] = (-0.3827327723098713*nuVtSqSum[0]*fr[4])+0.3827327723098713*nuVtSqSum[0]*fl[4]-0.3827327723098713*nuVtSqSum[1]*fr[2]+0.3827327723098713*nuVtSqSum[1]*fl[2]+0.3977475644174328*fr[0]*nuVtSqSum[1]+0.3977475644174328*fl[0]*nuVtSqSum[1]-0.7954951288348656*fc[0]*nuVtSqSum[1]+0.3977475644174328*nuVtSqSum[0]*fr[1]+0.3977475644174328*nuVtSqSum[0]*fl[1]-0.7954951288348656*nuVtSqSum[0]*fc[1]; 
  diff_incr[2] = (-0.3093592167691144*nuVtSqSum[1]*fr[4])-0.3093592167691144*nuVtSqSum[1]*fl[4]-2.032931995911323*nuVtSqSum[1]*fc[4]-0.3093592167691144*nuVtSqSum[0]*fr[2]-0.3093592167691144*nuVtSqSum[0]*fl[2]-2.032931995911323*nuVtSqSum[0]*fc[2]+0.3827327723098713*fr[1]*nuVtSqSum[1]-0.3827327723098713*fl[1]*nuVtSqSum[1]+0.3827327723098713*fr[0]*nuVtSqSum[0]-0.3827327723098713*fl[0]*nuVtSqSum[0]; 
  diff_incr[3] = (-0.3827327723098713*nuVtSqSum[1]*fr[7])+0.3827327723098713*nuVtSqSum[1]*fl[7]-0.3827327723098713*nuVtSqSum[0]*fr[6]+0.3827327723098713*nuVtSqSum[0]*fl[6]+0.3977475644174328*nuVtSqSum[1]*fr[5]+0.3977475644174328*nuVtSqSum[1]*fl[5]-0.7954951288348656*nuVtSqSum[1]*fc[5]+0.3977475644174328*nuVtSqSum[0]*fr[3]+0.3977475644174328*nuVtSqSum[0]*fl[3]-0.7954951288348656*nuVtSqSum[0]*fc[3]; 
  diff_incr[4] = (-0.3093592167691144*nuVtSqSum[0]*fr[4])-0.3093592167691144*nuVtSqSum[0]*fl[4]-2.032931995911323*nuVtSqSum[0]*fc[4]-0.3093592167691144*nuVtSqSum[1]*fr[2]-0.3093592167691144*nuVtSqSum[1]*fl[2]-2.032931995911323*nuVtSqSum[1]*fc[2]+0.3827327723098713*fr[0]*nuVtSqSum[1]-0.3827327723098713*fl[0]*nuVtSqSum[1]+0.3827327723098713*nuVtSqSum[0]*fr[1]-0.3827327723098713*nuVtSqSum[0]*fl[1]; 
  diff_incr[5] = (-0.3827327723098713*nuVtSqSum[0]*fr[7])+0.3827327723098713*nuVtSqSum[0]*fl[7]-0.3827327723098713*nuVtSqSum[1]*fr[6]+0.3827327723098713*nuVtSqSum[1]*fl[6]+0.3977475644174328*nuVtSqSum[0]*fr[5]+0.3977475644174328*nuVtSqSum[0]*fl[5]-0.7954951288348656*nuVtSqSum[0]*fc[5]+0.3977475644174328*nuVtSqSum[1]*fr[3]+0.3977475644174328*nuVtSqSum[1]*fl[3]-0.7954951288348656*nuVtSqSum[1]*fc[3]; 
  diff_incr[6] = (-0.3093592167691144*nuVtSqSum[1]*fr[7])-0.3093592167691144*nuVtSqSum[1]*fl[7]-2.032931995911323*nuVtSqSum[1]*fc[7]-0.3093592167691144*nuVtSqSum[0]*fr[6]-0.3093592167691144*nuVtSqSum[0]*fl[6]-2.032931995911323*nuVtSqSum[0]*fc[6]+0.3827327723098713*nuVtSqSum[1]*fr[5]-0.3827327723098713*nuVtSqSum[1]*fl[5]+0.3827327723098713*nuVtSqSum[0]*fr[3]-0.3827327723098713*nuVtSqSum[0]*fl[3]; 
  diff_incr[7] = (-0.3093592167691144*nuVtSqSum[0]*fr[7])-0.3093592167691144*nuVtSqSum[0]*fl[7]-2.032931995911323*nuVtSqSum[0]*fc[7]-0.3093592167691144*nuVtSqSum[1]*fr[6]-0.3093592167691144*nuVtSqSum[1]*fl[6]-2.032931995911323*nuVtSqSum[1]*fc[6]+0.3827327723098713*nuVtSqSum[0]*fr[5]-0.3827327723098713*nuVtSqSum[0]*fl[5]+0.3827327723098713*nuVtSqSum[1]*fr[3]-0.3827327723098713*nuVtSqSum[1]*fl[3]; 

  out[0] += diff_incr[0]*rdvSq4; 
  out[1] += diff_incr[1]*rdvSq4; 
  out[2] += diff_incr[2]*rdvSq4; 
  out[3] += diff_incr[3]*rdvSq4; 
  out[4] += diff_incr[4]*rdvSq4; 
  out[5] += diff_incr[5]*rdvSq4; 
  out[6] += diff_incr[6]*rdvSq4; 
  out[7] += diff_incr[7]*rdvSq4; 
} 
