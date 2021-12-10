#include <gkyl_vlasov_lbo_kernels.h> 
GKYL_CU_DH void vlasov_lbo_diff_surfvx_1x1v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[2]:         cell-center coordinates. 
  // dxv[2]:       cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[2]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[2]: sum of thermal speeds squared time their respective collisionalities. 
  // fl/fc/fr:      distribution function in cells 
  // out:           incremented distribution function in cell 
  double rdvSq4 = 4.0/(dxv[1]*dxv[1]); 
  double diff_incr[4] = {0.0}; 

  diff_incr[0] = (-0.3827327723098713*nuVtSqSum[1]*fr[3])+0.3827327723098713*nuVtSqSum[1]*fl[3]-0.3827327723098713*nuVtSqSum[0]*fr[2]+0.3827327723098713*nuVtSqSum[0]*fl[2]+0.3977475644174328*fr[1]*nuVtSqSum[1]+0.3977475644174328*fl[1]*nuVtSqSum[1]-0.7954951288348656*fc[1]*nuVtSqSum[1]+0.3977475644174328*fr[0]*nuVtSqSum[0]+0.3977475644174328*fl[0]*nuVtSqSum[0]-0.7954951288348656*fc[0]*nuVtSqSum[0]; 
  diff_incr[1] = (-0.3827327723098713*nuVtSqSum[0]*fr[3])+0.3827327723098713*nuVtSqSum[0]*fl[3]-0.3827327723098713*nuVtSqSum[1]*fr[2]+0.3827327723098713*nuVtSqSum[1]*fl[2]+0.3977475644174328*fr[0]*nuVtSqSum[1]+0.3977475644174328*fl[0]*nuVtSqSum[1]-0.7954951288348656*fc[0]*nuVtSqSum[1]+0.3977475644174328*nuVtSqSum[0]*fr[1]+0.3977475644174328*nuVtSqSum[0]*fl[1]-0.7954951288348656*nuVtSqSum[0]*fc[1]; 
  diff_incr[2] = (-0.3093592167691144*nuVtSqSum[1]*fr[3])-0.3093592167691144*nuVtSqSum[1]*fl[3]-2.032931995911323*nuVtSqSum[1]*fc[3]-0.3093592167691144*nuVtSqSum[0]*fr[2]-0.3093592167691144*nuVtSqSum[0]*fl[2]-2.032931995911323*nuVtSqSum[0]*fc[2]+0.3827327723098713*fr[1]*nuVtSqSum[1]-0.3827327723098713*fl[1]*nuVtSqSum[1]+0.3827327723098713*fr[0]*nuVtSqSum[0]-0.3827327723098713*fl[0]*nuVtSqSum[0]; 
  diff_incr[3] = (-0.3093592167691144*nuVtSqSum[0]*fr[3])-0.3093592167691144*nuVtSqSum[0]*fl[3]-2.032931995911323*nuVtSqSum[0]*fc[3]-0.3093592167691144*nuVtSqSum[1]*fr[2]-0.3093592167691144*nuVtSqSum[1]*fl[2]-2.032931995911323*nuVtSqSum[1]*fc[2]+0.3827327723098713*fr[0]*nuVtSqSum[1]-0.3827327723098713*fl[0]*nuVtSqSum[1]+0.3827327723098713*nuVtSqSum[0]*fr[1]-0.3827327723098713*nuVtSqSum[0]*fl[1]; 

  out[0] += diff_incr[0]*rdvSq4; 
  out[1] += diff_incr[1]*rdvSq4; 
  out[2] += diff_incr[2]*rdvSq4; 
  out[3] += diff_incr[3]*rdvSq4; 
} 
