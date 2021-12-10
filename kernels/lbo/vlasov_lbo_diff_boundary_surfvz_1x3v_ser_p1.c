#include <gkyl_vlasov_lbo_kernels.h> 
GKYL_CU_DH void vlasov_lbo_diff_boundary_surfvz_1x3v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[4]:         Cell-center coordinates. 
  // dxv[4]:       Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[6]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[2]: sum of thermal speeds squared time their respective collisionalities. 
  // fSkin/Edge:    Distribution function in cells 
  // out:           Incremented distribution function in cell 
  double rdvSq4 = 4.0/(dxv[3]*dxv[3]); 

  double diff_incr[16] = {0.0}; 

  if (edge == -1) { 

  diff_incr[0] = 0.2296396633859228*nuVtSqSum[1]*fSkin[8]-0.3827327723098713*nuVtSqSum[1]*fEdge[8]+0.2296396633859228*nuVtSqSum[0]*fSkin[4]-0.3827327723098713*nuVtSqSum[0]*fEdge[4]-0.7513009550107064*fSkin[1]*nuVtSqSum[1]+0.3977475644174328*fEdge[1]*nuVtSqSum[1]-0.7513009550107064*fSkin[0]*nuVtSqSum[0]+0.3977475644174328*fEdge[0]*nuVtSqSum[0]; 
  diff_incr[1] = 0.2296396633859228*nuVtSqSum[0]*fSkin[8]-0.3827327723098713*nuVtSqSum[0]*fEdge[8]+0.2296396633859228*nuVtSqSum[1]*fSkin[4]-0.3827327723098713*nuVtSqSum[1]*fEdge[4]-0.7513009550107064*fSkin[0]*nuVtSqSum[1]+0.3977475644174328*fEdge[0]*nuVtSqSum[1]-0.7513009550107064*nuVtSqSum[0]*fSkin[1]+0.3977475644174328*nuVtSqSum[0]*fEdge[1]; 
  diff_incr[2] = 0.2296396633859228*nuVtSqSum[1]*fSkin[12]-0.3827327723098713*nuVtSqSum[1]*fEdge[12]+0.2296396633859228*nuVtSqSum[0]*fSkin[9]-0.3827327723098713*nuVtSqSum[0]*fEdge[9]-0.7513009550107064*nuVtSqSum[1]*fSkin[5]+0.3977475644174328*nuVtSqSum[1]*fEdge[5]-0.7513009550107064*nuVtSqSum[0]*fSkin[2]+0.3977475644174328*nuVtSqSum[0]*fEdge[2]; 
  diff_incr[3] = 0.2296396633859228*nuVtSqSum[1]*fSkin[13]-0.3827327723098713*nuVtSqSum[1]*fEdge[13]+0.2296396633859228*nuVtSqSum[0]*fSkin[10]-0.3827327723098713*nuVtSqSum[0]*fEdge[10]-0.7513009550107064*nuVtSqSum[1]*fSkin[6]+0.3977475644174328*nuVtSqSum[1]*fEdge[6]-0.7513009550107064*nuVtSqSum[0]*fSkin[3]+0.3977475644174328*nuVtSqSum[0]*fEdge[3]; 
  diff_incr[4] = (-2.077126169735482*nuVtSqSum[1]*fSkin[8])-0.3093592167691144*nuVtSqSum[1]*fEdge[8]-2.077126169735482*nuVtSqSum[0]*fSkin[4]-0.3093592167691144*nuVtSqSum[0]*fEdge[4]-0.3827327723098713*fSkin[1]*nuVtSqSum[1]+0.3827327723098713*fEdge[1]*nuVtSqSum[1]-0.3827327723098713*fSkin[0]*nuVtSqSum[0]+0.3827327723098713*fEdge[0]*nuVtSqSum[0]; 
  diff_incr[5] = 0.2296396633859228*nuVtSqSum[0]*fSkin[12]-0.3827327723098713*nuVtSqSum[0]*fEdge[12]+0.2296396633859228*nuVtSqSum[1]*fSkin[9]-0.3827327723098713*nuVtSqSum[1]*fEdge[9]-0.7513009550107064*nuVtSqSum[0]*fSkin[5]+0.3977475644174328*nuVtSqSum[0]*fEdge[5]-0.7513009550107064*nuVtSqSum[1]*fSkin[2]+0.3977475644174328*nuVtSqSum[1]*fEdge[2]; 
  diff_incr[6] = 0.2296396633859228*nuVtSqSum[0]*fSkin[13]-0.3827327723098713*nuVtSqSum[0]*fEdge[13]+0.2296396633859228*nuVtSqSum[1]*fSkin[10]-0.3827327723098713*nuVtSqSum[1]*fEdge[10]-0.7513009550107064*nuVtSqSum[0]*fSkin[6]+0.3977475644174328*nuVtSqSum[0]*fEdge[6]-0.7513009550107064*nuVtSqSum[1]*fSkin[3]+0.3977475644174328*nuVtSqSum[1]*fEdge[3]; 
  diff_incr[7] = 0.2296396633859228*nuVtSqSum[1]*fSkin[15]-0.3827327723098713*nuVtSqSum[1]*fEdge[15]+0.2296396633859228*nuVtSqSum[0]*fSkin[14]-0.3827327723098713*nuVtSqSum[0]*fEdge[14]-0.7513009550107064*nuVtSqSum[1]*fSkin[11]+0.3977475644174328*nuVtSqSum[1]*fEdge[11]-0.7513009550107064*nuVtSqSum[0]*fSkin[7]+0.3977475644174328*nuVtSqSum[0]*fEdge[7]; 
  diff_incr[8] = (-2.077126169735482*nuVtSqSum[0]*fSkin[8])-0.3093592167691144*nuVtSqSum[0]*fEdge[8]-2.077126169735482*nuVtSqSum[1]*fSkin[4]-0.3093592167691144*nuVtSqSum[1]*fEdge[4]-0.3827327723098713*fSkin[0]*nuVtSqSum[1]+0.3827327723098713*fEdge[0]*nuVtSqSum[1]-0.3827327723098713*nuVtSqSum[0]*fSkin[1]+0.3827327723098713*nuVtSqSum[0]*fEdge[1]; 
  diff_incr[9] = (-2.077126169735482*nuVtSqSum[1]*fSkin[12])-0.3093592167691144*nuVtSqSum[1]*fEdge[12]-2.077126169735482*nuVtSqSum[0]*fSkin[9]-0.3093592167691144*nuVtSqSum[0]*fEdge[9]-0.3827327723098713*nuVtSqSum[1]*fSkin[5]+0.3827327723098713*nuVtSqSum[1]*fEdge[5]-0.3827327723098713*nuVtSqSum[0]*fSkin[2]+0.3827327723098713*nuVtSqSum[0]*fEdge[2]; 
  diff_incr[10] = (-2.077126169735482*nuVtSqSum[1]*fSkin[13])-0.3093592167691144*nuVtSqSum[1]*fEdge[13]-2.077126169735482*nuVtSqSum[0]*fSkin[10]-0.3093592167691144*nuVtSqSum[0]*fEdge[10]-0.3827327723098713*nuVtSqSum[1]*fSkin[6]+0.3827327723098713*nuVtSqSum[1]*fEdge[6]-0.3827327723098713*nuVtSqSum[0]*fSkin[3]+0.3827327723098713*nuVtSqSum[0]*fEdge[3]; 
  diff_incr[11] = 0.2296396633859228*nuVtSqSum[0]*fSkin[15]-0.3827327723098713*nuVtSqSum[0]*fEdge[15]+0.2296396633859228*nuVtSqSum[1]*fSkin[14]-0.3827327723098713*nuVtSqSum[1]*fEdge[14]-0.7513009550107064*nuVtSqSum[0]*fSkin[11]+0.3977475644174328*nuVtSqSum[0]*fEdge[11]-0.7513009550107064*nuVtSqSum[1]*fSkin[7]+0.3977475644174328*nuVtSqSum[1]*fEdge[7]; 
  diff_incr[12] = (-2.077126169735482*nuVtSqSum[0]*fSkin[12])-0.3093592167691144*nuVtSqSum[0]*fEdge[12]-2.077126169735482*nuVtSqSum[1]*fSkin[9]-0.3093592167691144*nuVtSqSum[1]*fEdge[9]-0.3827327723098713*nuVtSqSum[0]*fSkin[5]+0.3827327723098713*nuVtSqSum[0]*fEdge[5]-0.3827327723098713*nuVtSqSum[1]*fSkin[2]+0.3827327723098713*nuVtSqSum[1]*fEdge[2]; 
  diff_incr[13] = (-2.077126169735482*nuVtSqSum[0]*fSkin[13])-0.3093592167691144*nuVtSqSum[0]*fEdge[13]-2.077126169735482*nuVtSqSum[1]*fSkin[10]-0.3093592167691144*nuVtSqSum[1]*fEdge[10]-0.3827327723098713*nuVtSqSum[0]*fSkin[6]+0.3827327723098713*nuVtSqSum[0]*fEdge[6]-0.3827327723098713*nuVtSqSum[1]*fSkin[3]+0.3827327723098713*nuVtSqSum[1]*fEdge[3]; 
  diff_incr[14] = (-2.077126169735482*nuVtSqSum[1]*fSkin[15])-0.3093592167691144*nuVtSqSum[1]*fEdge[15]-2.077126169735482*nuVtSqSum[0]*fSkin[14]-0.3093592167691144*nuVtSqSum[0]*fEdge[14]-0.3827327723098713*nuVtSqSum[1]*fSkin[11]+0.3827327723098713*nuVtSqSum[1]*fEdge[11]-0.3827327723098713*nuVtSqSum[0]*fSkin[7]+0.3827327723098713*nuVtSqSum[0]*fEdge[7]; 
  diff_incr[15] = (-2.077126169735482*nuVtSqSum[0]*fSkin[15])-0.3093592167691144*nuVtSqSum[0]*fEdge[15]-2.077126169735482*nuVtSqSum[1]*fSkin[14]-0.3093592167691144*nuVtSqSum[1]*fEdge[14]-0.3827327723098713*nuVtSqSum[0]*fSkin[11]+0.3827327723098713*nuVtSqSum[0]*fEdge[11]-0.3827327723098713*nuVtSqSum[1]*fSkin[7]+0.3827327723098713*nuVtSqSum[1]*fEdge[7]; 


  } else { 

  diff_incr[0] = 0.9951052080056654*nuVtSqSum[1]*fSkin[8]+0.3827327723098713*nuVtSqSum[1]*fEdge[8]+0.9951052080056654*nuVtSqSum[0]*fSkin[4]+0.3827327723098713*nuVtSqSum[0]*fEdge[4]-0.0441941738241592*fSkin[1]*nuVtSqSum[1]+0.3977475644174328*fEdge[1]*nuVtSqSum[1]-0.0441941738241592*fSkin[0]*nuVtSqSum[0]+0.3977475644174328*fEdge[0]*nuVtSqSum[0]; 
  diff_incr[1] = 0.9951052080056654*nuVtSqSum[0]*fSkin[8]+0.3827327723098713*nuVtSqSum[0]*fEdge[8]+0.9951052080056654*nuVtSqSum[1]*fSkin[4]+0.3827327723098713*nuVtSqSum[1]*fEdge[4]-0.0441941738241592*fSkin[0]*nuVtSqSum[1]+0.3977475644174328*fEdge[0]*nuVtSqSum[1]-0.0441941738241592*nuVtSqSum[0]*fSkin[1]+0.3977475644174328*nuVtSqSum[0]*fEdge[1]; 
  diff_incr[2] = 0.9951052080056654*nuVtSqSum[1]*fSkin[12]+0.3827327723098713*nuVtSqSum[1]*fEdge[12]+0.9951052080056654*nuVtSqSum[0]*fSkin[9]+0.3827327723098713*nuVtSqSum[0]*fEdge[9]-0.0441941738241592*nuVtSqSum[1]*fSkin[5]+0.3977475644174328*nuVtSqSum[1]*fEdge[5]-0.0441941738241592*nuVtSqSum[0]*fSkin[2]+0.3977475644174328*nuVtSqSum[0]*fEdge[2]; 
  diff_incr[3] = 0.9951052080056654*nuVtSqSum[1]*fSkin[13]+0.3827327723098713*nuVtSqSum[1]*fEdge[13]+0.9951052080056654*nuVtSqSum[0]*fSkin[10]+0.3827327723098713*nuVtSqSum[0]*fEdge[10]-0.0441941738241592*nuVtSqSum[1]*fSkin[6]+0.3977475644174328*nuVtSqSum[1]*fEdge[6]-0.0441941738241592*nuVtSqSum[0]*fSkin[3]+0.3977475644174328*nuVtSqSum[0]*fEdge[3]; 
  diff_incr[4] = 0.0441941738241592*nuVtSqSum[1]*fSkin[8]-0.3093592167691144*nuVtSqSum[1]*fEdge[8]+0.0441941738241592*nuVtSqSum[0]*fSkin[4]-0.3093592167691144*nuVtSqSum[0]*fEdge[4]+1.60747764370146*fSkin[1]*nuVtSqSum[1]-0.3827327723098713*fEdge[1]*nuVtSqSum[1]+1.60747764370146*fSkin[0]*nuVtSqSum[0]-0.3827327723098713*fEdge[0]*nuVtSqSum[0]; 
  diff_incr[5] = 0.9951052080056654*nuVtSqSum[0]*fSkin[12]+0.3827327723098713*nuVtSqSum[0]*fEdge[12]+0.9951052080056654*nuVtSqSum[1]*fSkin[9]+0.3827327723098713*nuVtSqSum[1]*fEdge[9]-0.0441941738241592*nuVtSqSum[0]*fSkin[5]+0.3977475644174328*nuVtSqSum[0]*fEdge[5]-0.0441941738241592*nuVtSqSum[1]*fSkin[2]+0.3977475644174328*nuVtSqSum[1]*fEdge[2]; 
  diff_incr[6] = 0.9951052080056654*nuVtSqSum[0]*fSkin[13]+0.3827327723098713*nuVtSqSum[0]*fEdge[13]+0.9951052080056654*nuVtSqSum[1]*fSkin[10]+0.3827327723098713*nuVtSqSum[1]*fEdge[10]-0.0441941738241592*nuVtSqSum[0]*fSkin[6]+0.3977475644174328*nuVtSqSum[0]*fEdge[6]-0.0441941738241592*nuVtSqSum[1]*fSkin[3]+0.3977475644174328*nuVtSqSum[1]*fEdge[3]; 
  diff_incr[7] = 0.9951052080056654*nuVtSqSum[1]*fSkin[15]+0.3827327723098713*nuVtSqSum[1]*fEdge[15]+0.9951052080056654*nuVtSqSum[0]*fSkin[14]+0.3827327723098713*nuVtSqSum[0]*fEdge[14]-0.0441941738241592*nuVtSqSum[1]*fSkin[11]+0.3977475644174328*nuVtSqSum[1]*fEdge[11]-0.0441941738241592*nuVtSqSum[0]*fSkin[7]+0.3977475644174328*nuVtSqSum[0]*fEdge[7]; 
  diff_incr[8] = 0.0441941738241592*nuVtSqSum[0]*fSkin[8]-0.3093592167691144*nuVtSqSum[0]*fEdge[8]+0.0441941738241592*nuVtSqSum[1]*fSkin[4]-0.3093592167691144*nuVtSqSum[1]*fEdge[4]+1.60747764370146*fSkin[0]*nuVtSqSum[1]-0.3827327723098713*fEdge[0]*nuVtSqSum[1]+1.60747764370146*nuVtSqSum[0]*fSkin[1]-0.3827327723098713*nuVtSqSum[0]*fEdge[1]; 
  diff_incr[9] = 0.0441941738241592*nuVtSqSum[1]*fSkin[12]-0.3093592167691144*nuVtSqSum[1]*fEdge[12]+0.0441941738241592*nuVtSqSum[0]*fSkin[9]-0.3093592167691144*nuVtSqSum[0]*fEdge[9]+1.60747764370146*nuVtSqSum[1]*fSkin[5]-0.3827327723098713*nuVtSqSum[1]*fEdge[5]+1.60747764370146*nuVtSqSum[0]*fSkin[2]-0.3827327723098713*nuVtSqSum[0]*fEdge[2]; 
  diff_incr[10] = 0.0441941738241592*nuVtSqSum[1]*fSkin[13]-0.3093592167691144*nuVtSqSum[1]*fEdge[13]+0.0441941738241592*nuVtSqSum[0]*fSkin[10]-0.3093592167691144*nuVtSqSum[0]*fEdge[10]+1.60747764370146*nuVtSqSum[1]*fSkin[6]-0.3827327723098713*nuVtSqSum[1]*fEdge[6]+1.60747764370146*nuVtSqSum[0]*fSkin[3]-0.3827327723098713*nuVtSqSum[0]*fEdge[3]; 
  diff_incr[11] = 0.9951052080056654*nuVtSqSum[0]*fSkin[15]+0.3827327723098713*nuVtSqSum[0]*fEdge[15]+0.9951052080056654*nuVtSqSum[1]*fSkin[14]+0.3827327723098713*nuVtSqSum[1]*fEdge[14]-0.0441941738241592*nuVtSqSum[0]*fSkin[11]+0.3977475644174328*nuVtSqSum[0]*fEdge[11]-0.0441941738241592*nuVtSqSum[1]*fSkin[7]+0.3977475644174328*nuVtSqSum[1]*fEdge[7]; 
  diff_incr[12] = 0.0441941738241592*nuVtSqSum[0]*fSkin[12]-0.3093592167691144*nuVtSqSum[0]*fEdge[12]+0.0441941738241592*nuVtSqSum[1]*fSkin[9]-0.3093592167691144*nuVtSqSum[1]*fEdge[9]+1.60747764370146*nuVtSqSum[0]*fSkin[5]-0.3827327723098713*nuVtSqSum[0]*fEdge[5]+1.60747764370146*nuVtSqSum[1]*fSkin[2]-0.3827327723098713*nuVtSqSum[1]*fEdge[2]; 
  diff_incr[13] = 0.0441941738241592*nuVtSqSum[0]*fSkin[13]-0.3093592167691144*nuVtSqSum[0]*fEdge[13]+0.0441941738241592*nuVtSqSum[1]*fSkin[10]-0.3093592167691144*nuVtSqSum[1]*fEdge[10]+1.60747764370146*nuVtSqSum[0]*fSkin[6]-0.3827327723098713*nuVtSqSum[0]*fEdge[6]+1.60747764370146*nuVtSqSum[1]*fSkin[3]-0.3827327723098713*nuVtSqSum[1]*fEdge[3]; 
  diff_incr[14] = 0.0441941738241592*nuVtSqSum[1]*fSkin[15]-0.3093592167691144*nuVtSqSum[1]*fEdge[15]+0.0441941738241592*nuVtSqSum[0]*fSkin[14]-0.3093592167691144*nuVtSqSum[0]*fEdge[14]+1.60747764370146*nuVtSqSum[1]*fSkin[11]-0.3827327723098713*nuVtSqSum[1]*fEdge[11]+1.60747764370146*nuVtSqSum[0]*fSkin[7]-0.3827327723098713*nuVtSqSum[0]*fEdge[7]; 
  diff_incr[15] = 0.0441941738241592*nuVtSqSum[0]*fSkin[15]-0.3093592167691144*nuVtSqSum[0]*fEdge[15]+0.0441941738241592*nuVtSqSum[1]*fSkin[14]-0.3093592167691144*nuVtSqSum[1]*fEdge[14]+1.60747764370146*nuVtSqSum[0]*fSkin[11]-0.3827327723098713*nuVtSqSum[0]*fEdge[11]+1.60747764370146*nuVtSqSum[1]*fSkin[7]-0.3827327723098713*nuVtSqSum[1]*fEdge[7]; 

  } 

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
  out[12] += diff_incr[12]*rdvSq4; 
  out[13] += diff_incr[13]*rdvSq4; 
  out[14] += diff_incr[14]*rdvSq4; 
  out[15] += diff_incr[15]*rdvSq4; 
} 
