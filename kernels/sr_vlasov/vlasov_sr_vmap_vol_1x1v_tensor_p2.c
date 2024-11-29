#include <gkyl_vlasov_sr_kernels.h> 
GKYL_CU_DH double vlasov_sr_vmap_vol_1x1v_tensor_p2(const double *w, const double *dxv, const double *jacob_vel_inv, const double *gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:             Cell-center coordinates.
  // dxv[NDIM]:           Cell spacing.
  // jacob_vel_inv[VDIM]: Inverse velocity space Jacobian in each direction.
  // gamma:               Particle Lorentz boost factor sqrt(1 + p^2).
  // qmem:                q/m*EM fields.
  // f:                   Input distribution function.
  // out:                 Incremented output.
  const double dx10 = 2.0/dxv[0]; 
  const double dv10 = 2.0/dxv[1]; 
  const double *E0 = &qmem[0]; 
  const double *jacob_vel_inv0 = &jacob_vel_inv[0]; 
  double p0_over_gamma[3] = {0.0}; 
  p0_over_gamma[0] = 2.738612787525831*jacob_vel_inv0[1]*gamma[2]*dv10+1.224744871391589*jacob_vel_inv0[0]*gamma[1]*dv10; 
  p0_over_gamma[1] = 2.449489742783178*jacob_vel_inv0[2]*gamma[2]*dv10+2.738612787525831*jacob_vel_inv0[0]*gamma[2]*dv10+1.224744871391589*jacob_vel_inv0[1]*gamma[1]*dv10; 
  p0_over_gamma[2] = 2.449489742783178*jacob_vel_inv0[1]*gamma[2]*dv10+1.224744871391589*gamma[1]*jacob_vel_inv0[2]*dv10; 

  double cflFreq_mid = 0.0; 
  double alpha_vdim[9] = {0.0}; 

  cflFreq_mid += 5.0*fabs((0.3535533905932737*p0_over_gamma[0]-0.3952847075210473*p0_over_gamma[2])*dx10); 

  out[1] += 1.224744871391589*(p0_over_gamma[2]*f[5]+p0_over_gamma[1]*f[2]+f[0]*p0_over_gamma[0])*dx10; 
  out[3] += (1.095445115010332*(p0_over_gamma[1]*f[5]+f[2]*p0_over_gamma[2])+1.224744871391589*(p0_over_gamma[0]*f[2]+f[0]*p0_over_gamma[1]))*dx10; 
  out[4] += 2.738612787525831*(p0_over_gamma[2]*f[7]+p0_over_gamma[1]*f[3]+p0_over_gamma[0]*f[1])*dx10; 
  out[6] += (2.449489742783178*(p0_over_gamma[1]*f[7]+p0_over_gamma[2]*f[3])+2.738612787525831*(p0_over_gamma[0]*f[3]+f[1]*p0_over_gamma[1]))*dx10; 
  out[7] += (0.7824607964359517*p0_over_gamma[2]*f[5]+1.224744871391589*(p0_over_gamma[0]*f[5]+f[0]*p0_over_gamma[2])+1.095445115010332*p0_over_gamma[1]*f[2])*dx10; 
  out[8] += ((1.749635530559413*p0_over_gamma[2]+2.738612787525831*p0_over_gamma[0])*f[7]+2.449489742783178*p0_over_gamma[1]*f[3]+2.738612787525831*f[1]*p0_over_gamma[2])*dx10; 

  alpha_vdim[0] = E0[0]*jacob_vel_inv0[0]*dv10; 
  alpha_vdim[1] = jacob_vel_inv0[0]*E0[1]*dv10; 
  alpha_vdim[2] = E0[0]*jacob_vel_inv0[1]*dv10; 
  alpha_vdim[3] = E0[1]*jacob_vel_inv0[1]*dv10; 
  alpha_vdim[4] = jacob_vel_inv0[0]*E0[2]*dv10; 
  alpha_vdim[5] = E0[0]*jacob_vel_inv0[2]*dv10; 
  alpha_vdim[6] = jacob_vel_inv0[1]*E0[2]*dv10; 
  alpha_vdim[7] = E0[1]*jacob_vel_inv0[2]*dv10; 
  alpha_vdim[8] = E0[2]*jacob_vel_inv0[2]*dv10; 
  cflFreq_mid += 5.0*fabs(0.3125*alpha_vdim[8]-0.2795084971874737*(alpha_vdim[5]+alpha_vdim[4])+0.25*alpha_vdim[0]); 

  out[2] += 0.8660254037844386*(alpha_vdim[8]*f[8]+alpha_vdim[7]*f[7]+alpha_vdim[6]*f[6]+alpha_vdim[5]*f[5]+alpha_vdim[4]*f[4]+alpha_vdim[3]*f[3]+alpha_vdim[2]*f[2]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[3] += 0.7745966692414834*(alpha_vdim[7]*f[8]+f[7]*alpha_vdim[8])+0.8660254037844387*(alpha_vdim[5]*f[7]+f[5]*alpha_vdim[7])+0.7745966692414834*(alpha_vdim[3]*f[6]+f[3]*alpha_vdim[6])+0.7745966692414833*(alpha_vdim[1]*f[4]+f[1]*alpha_vdim[4])+0.8660254037844386*(alpha_vdim[2]*f[3]+f[2]*alpha_vdim[3]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[5] += 1.732050807568877*(alpha_vdim[6]*f[8]+f[6]*alpha_vdim[8]+alpha_vdim[3]*f[7]+f[3]*alpha_vdim[7])+1.936491673103709*(alpha_vdim[4]*f[6]+f[4]*alpha_vdim[6])+1.732050807568877*(alpha_vdim[2]*f[5]+f[2]*alpha_vdim[5])+1.936491673103709*(alpha_vdim[1]*f[3]+f[1]*alpha_vdim[3]+alpha_vdim[0]*f[2]+f[0]*alpha_vdim[2]); 
  out[6] += 0.5532833351724881*alpha_vdim[8]*f[8]+0.8660254037844387*(alpha_vdim[5]*f[8]+f[5]*alpha_vdim[8])+0.7745966692414834*alpha_vdim[7]*f[7]+0.5532833351724881*alpha_vdim[6]*f[6]+0.8660254037844386*(alpha_vdim[2]*f[6]+f[2]*alpha_vdim[6])+0.5532833351724881*alpha_vdim[4]*f[4]+0.8660254037844387*(alpha_vdim[0]*f[4]+f[0]*alpha_vdim[4])+0.7745966692414834*(alpha_vdim[3]*f[3]+alpha_vdim[1]*f[1]); 
  out[7] += 1.549193338482967*(alpha_vdim[3]*f[8]+f[3]*alpha_vdim[8])+(1.549193338482967*alpha_vdim[6]+1.732050807568877*alpha_vdim[2])*f[7]+1.549193338482967*f[6]*alpha_vdim[7]+1.732050807568877*(f[2]*alpha_vdim[7]+alpha_vdim[1]*f[6]+f[1]*alpha_vdim[6])+1.732050807568877*(alpha_vdim[3]*f[5]+f[3]*alpha_vdim[5]+alpha_vdim[3]*f[4]+f[3]*alpha_vdim[4])+1.936491673103709*(alpha_vdim[0]*f[3]+f[0]*alpha_vdim[3]+alpha_vdim[1]*f[2]+f[1]*alpha_vdim[2]); 
  out[8] += (1.106566670344976*alpha_vdim[6]+1.732050807568877*alpha_vdim[2])*f[8]+(1.106566670344976*f[6]+1.732050807568877*f[2])*alpha_vdim[8]+1.549193338482967*(alpha_vdim[3]*f[7]+f[3]*alpha_vdim[7])+(1.732050807568877*alpha_vdim[5]+1.237179148263484*alpha_vdim[4]+1.936491673103709*alpha_vdim[0])*f[6]+(1.732050807568877*f[5]+1.237179148263484*f[4])*alpha_vdim[6]+1.936491673103709*(f[0]*alpha_vdim[6]+alpha_vdim[2]*f[4]+f[2]*alpha_vdim[4])+1.732050807568877*(alpha_vdim[1]*f[3]+f[1]*alpha_vdim[3]); 

  return cflFreq_mid; 
} 
