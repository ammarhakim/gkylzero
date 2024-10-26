#include <gkyl_vlasov_sr_kernels.h> 
GKYL_CU_DH double vlasov_sr_vmap_stream_vol_1x1v_tensor_p2(const double *w, const double *dxv, const double *jacob_vel_inv, const double *gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:             Cell-center coordinates.
  // dxv[NDIM]:           Cell spacing.
  // jacob_vel_inv[VDIM]: Inverse velocity space Jacobian in each direction.
  // gamma:               Particle Lorentz boost factor sqrt(1 + p^2).
  // qmem:                q/m*EM fields (unused in streaming-only update).
  // f:                   Input distribution function.
  // out:                 Incremented output.
  const double dx10 = 2.0/dxv[0]; 
  const double dv10 = 2.0/dxv[1]; 
  const double *jacob_vel_inv0 = &jacob_vel_inv[0]; 
  double p0_over_gamma[3] = {0.0}; 
  p0_over_gamma[0] = 2.738612787525831*jacob_vel_inv0[1]*gamma[2]*dv10+1.224744871391589*jacob_vel_inv0[0]*gamma[1]*dv10; 
  p0_over_gamma[1] = 2.449489742783178*jacob_vel_inv0[2]*gamma[2]*dv10+2.738612787525831*jacob_vel_inv0[0]*gamma[2]*dv10+1.224744871391589*jacob_vel_inv0[1]*gamma[1]*dv10; 
  p0_over_gamma[2] = 2.449489742783178*jacob_vel_inv0[1]*gamma[2]*dv10+1.224744871391589*gamma[1]*jacob_vel_inv0[2]*dv10; 

  double cflFreq_mid = 0.0; 
  double alpha_cdim[9] = {0.0}; 
  alpha_cdim[0] = 1.414213562373095*p0_over_gamma[0]*dx10; 
  alpha_cdim[2] = 1.414213562373095*p0_over_gamma[1]*dx10; 
  alpha_cdim[5] = 1.414213562373095*p0_over_gamma[2]*dx10; 
  cflFreq_mid += fabs(0.25*alpha_cdim[0]-0.2795084971874737*alpha_cdim[5]); 

  out[1] += 0.8660254037844386*(alpha_cdim[5]*f[5]+alpha_cdim[2]*f[2]+alpha_cdim[0]*f[0]); 
  out[3] += 0.7745966692414833*(alpha_cdim[2]*f[5]+f[2]*alpha_cdim[5])+0.8660254037844386*(alpha_cdim[0]*f[2]+f[0]*alpha_cdim[2]); 
  out[4] += 1.936491673103709*(alpha_cdim[5]*f[7]+alpha_cdim[2]*f[3]+alpha_cdim[0]*f[1]); 
  out[6] += 1.732050807568877*(alpha_cdim[2]*f[7]+f[3]*alpha_cdim[5])+1.936491673103709*(alpha_cdim[0]*f[3]+f[1]*alpha_cdim[2]); 
  out[7] += 0.5532833351724881*alpha_cdim[5]*f[5]+0.8660254037844386*(alpha_cdim[0]*f[5]+f[0]*alpha_cdim[5])+0.7745966692414833*alpha_cdim[2]*f[2]; 
  out[8] += 1.237179148263484*alpha_cdim[5]*f[7]+1.936491673103709*(alpha_cdim[0]*f[7]+f[1]*alpha_cdim[5])+1.732050807568877*alpha_cdim[2]*f[3]; 

  return 5.0*cflFreq_mid; 
} 
