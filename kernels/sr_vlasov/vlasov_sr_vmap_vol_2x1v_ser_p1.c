#include <gkyl_vlasov_sr_kernels.h> 
GKYL_CU_DH double vlasov_sr_vmap_vol_2x1v_ser_p1(const double *w, const double *dxv, const double *jacob_vel_inv, const double *gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:             Cell-center coordinates.
  // dxv[NDIM]:           Cell spacing.
  // jacob_vel_inv[VDIM]: Inverse velocity space Jacobian in each direction.
  // gamma:               Particle Lorentz boost factor sqrt(1 + p^2).
  // qmem:                q/m*EM fields.
  // f:                   Input distribution function.
  // out:                 Incremented output.
  const double dx10 = 2.0/dxv[0]; 
  const double dx11 = 2.0/dxv[1]; 
  const double dv10 = 2.0/dxv[2]; 
  const double *E0 = &qmem[0]; 
  const double *jacob_vel_inv0 = &jacob_vel_inv[0]; 
  double p0_over_gamma[3] = {0.0}; 
  p0_over_gamma[0] = 2.738612787525831*jacob_vel_inv0[1]*gamma[2]*dv10+1.224744871391589*jacob_vel_inv0[0]*gamma[1]*dv10; 
  p0_over_gamma[1] = 2.449489742783178*jacob_vel_inv0[2]*gamma[2]*dv10+2.738612787525831*jacob_vel_inv0[0]*gamma[2]*dv10+1.224744871391589*jacob_vel_inv0[1]*gamma[1]*dv10; 
  p0_over_gamma[2] = 2.449489742783178*jacob_vel_inv0[1]*gamma[2]*dv10+1.224744871391589*gamma[1]*jacob_vel_inv0[2]*dv10; 

  double cflFreq_mid = 0.0; 
  double alpha_vdim[12] = {0.0}; 

  cflFreq_mid += 3.0*fabs((0.3535533905932737*p0_over_gamma[0]-0.3952847075210473*p0_over_gamma[2])*dx10); 

  out[1] += 1.224744871391589*(p0_over_gamma[2]*f[7]+p0_over_gamma[1]*f[3]+f[0]*p0_over_gamma[0])*dx10; 
  out[4] += 1.224744871391589*(p0_over_gamma[2]*f[10]+p0_over_gamma[1]*f[6]+p0_over_gamma[0]*f[2])*dx10; 
  out[5] += (1.095445115010332*(p0_over_gamma[1]*f[7]+p0_over_gamma[2]*f[3])+1.224744871391589*(p0_over_gamma[0]*f[3]+f[0]*p0_over_gamma[1]))*dx10; 
  out[8] += (1.095445115010332*(p0_over_gamma[1]*f[10]+p0_over_gamma[2]*f[6])+1.224744871391589*(p0_over_gamma[0]*f[6]+p0_over_gamma[1]*f[2]))*dx10; 
  out[9] += ((0.7824607964359517*p0_over_gamma[2]+1.224744871391589*p0_over_gamma[0])*f[7]+1.095445115010332*p0_over_gamma[1]*f[3]+1.224744871391589*f[0]*p0_over_gamma[2])*dx10; 
  out[11] += ((0.7824607964359517*p0_over_gamma[2]+1.224744871391589*p0_over_gamma[0])*f[10]+1.095445115010332*p0_over_gamma[1]*f[6]+1.224744871391589*f[2]*p0_over_gamma[2])*dx10; 

  alpha_vdim[0] = E0[0]*jacob_vel_inv0[0]*dv10; 
  alpha_vdim[1] = jacob_vel_inv0[0]*E0[1]*dv10; 
  alpha_vdim[2] = jacob_vel_inv0[0]*E0[2]*dv10; 
  alpha_vdim[3] = E0[0]*jacob_vel_inv0[1]*dv10; 
  alpha_vdim[4] = jacob_vel_inv0[0]*E0[3]*dv10; 
  alpha_vdim[5] = E0[1]*jacob_vel_inv0[1]*dv10; 
  alpha_vdim[6] = jacob_vel_inv0[1]*E0[2]*dv10; 
  alpha_vdim[7] = E0[0]*jacob_vel_inv0[2]*dv10; 
  alpha_vdim[8] = jacob_vel_inv0[1]*E0[3]*dv10; 
  alpha_vdim[9] = E0[1]*jacob_vel_inv0[2]*dv10; 
  alpha_vdim[10] = E0[2]*jacob_vel_inv0[2]*dv10; 
  alpha_vdim[11] = jacob_vel_inv0[2]*E0[3]*dv10; 
  cflFreq_mid += 5.0*fabs(0.1767766952966368*alpha_vdim[0]-0.1976423537605236*alpha_vdim[7]); 

  out[3] += 0.6123724356957944*(alpha_vdim[11]*f[11]+alpha_vdim[10]*f[10]+alpha_vdim[9]*f[9]+alpha_vdim[8]*f[8]+alpha_vdim[7]*f[7]+alpha_vdim[6]*f[6]+alpha_vdim[5]*f[5]+alpha_vdim[4]*f[4]+alpha_vdim[3]*f[3]+alpha_vdim[2]*f[2]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[5] += 0.6123724356957944*(alpha_vdim[10]*f[11]+f[10]*alpha_vdim[11]+alpha_vdim[7]*f[9]+f[7]*alpha_vdim[9]+alpha_vdim[6]*f[8]+f[6]*alpha_vdim[8]+alpha_vdim[3]*f[5]+f[3]*alpha_vdim[5]+alpha_vdim[2]*f[4]+f[2]*alpha_vdim[4]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[6] += 0.6123724356957944*(alpha_vdim[9]*f[11]+f[9]*alpha_vdim[11]+alpha_vdim[7]*f[10]+f[7]*alpha_vdim[10]+alpha_vdim[5]*f[8]+f[5]*alpha_vdim[8]+alpha_vdim[3]*f[6]+f[3]*alpha_vdim[6]+alpha_vdim[1]*f[4]+f[1]*alpha_vdim[4]+alpha_vdim[0]*f[2]+f[0]*alpha_vdim[2]); 
  out[7] += 1.224744871391589*(alpha_vdim[8]*f[11]+f[8]*alpha_vdim[11]+alpha_vdim[6]*f[10]+f[6]*alpha_vdim[10]+alpha_vdim[5]*f[9]+f[5]*alpha_vdim[9])+1.369306393762915*(alpha_vdim[4]*f[8]+f[4]*alpha_vdim[8])+1.224744871391589*(alpha_vdim[3]*f[7]+f[3]*alpha_vdim[7])+1.369306393762915*(alpha_vdim[2]*f[6]+f[2]*alpha_vdim[6]+alpha_vdim[1]*f[5]+f[1]*alpha_vdim[5]+alpha_vdim[0]*f[3]+f[0]*alpha_vdim[3]); 
  out[8] += 0.6123724356957944*(alpha_vdim[7]*f[11]+f[7]*alpha_vdim[11]+alpha_vdim[9]*f[10]+f[9]*alpha_vdim[10]+alpha_vdim[3]*f[8]+f[3]*alpha_vdim[8]+alpha_vdim[5]*f[6]+f[5]*alpha_vdim[6]+alpha_vdim[0]*f[4]+f[0]*alpha_vdim[4]+alpha_vdim[1]*f[2]+f[1]*alpha_vdim[2]); 
  out[9] += 1.224744871391589*(alpha_vdim[6]*f[11]+f[6]*alpha_vdim[11]+alpha_vdim[8]*f[10]+f[8]*alpha_vdim[10]+alpha_vdim[3]*f[9]+f[3]*alpha_vdim[9])+1.369306393762915*(alpha_vdim[2]*f[8]+f[2]*alpha_vdim[8])+1.224744871391589*(alpha_vdim[5]*f[7]+f[5]*alpha_vdim[7])+1.369306393762915*(alpha_vdim[4]*f[6]+f[4]*alpha_vdim[6]+alpha_vdim[0]*f[5]+f[0]*alpha_vdim[5]+alpha_vdim[1]*f[3]+f[1]*alpha_vdim[3]); 
  out[10] += 1.224744871391589*(alpha_vdim[5]*f[11]+f[5]*alpha_vdim[11]+alpha_vdim[3]*f[10]+f[3]*alpha_vdim[10]+alpha_vdim[8]*f[9]+f[8]*alpha_vdim[9])+1.369306393762915*(alpha_vdim[1]*f[8]+f[1]*alpha_vdim[8])+1.224744871391589*(alpha_vdim[6]*f[7]+f[6]*alpha_vdim[7])+1.369306393762915*(alpha_vdim[0]*f[6]+f[0]*alpha_vdim[6]+alpha_vdim[4]*f[5]+f[4]*alpha_vdim[5]+alpha_vdim[2]*f[3]+f[2]*alpha_vdim[3]); 
  out[11] += 1.224744871391589*(alpha_vdim[3]*f[11]+f[3]*alpha_vdim[11]+alpha_vdim[5]*f[10]+f[5]*alpha_vdim[10]+alpha_vdim[6]*f[9]+f[6]*alpha_vdim[9])+(1.224744871391589*alpha_vdim[7]+1.369306393762915*alpha_vdim[0])*f[8]+1.224744871391589*f[7]*alpha_vdim[8]+1.369306393762915*(f[0]*alpha_vdim[8]+alpha_vdim[1]*f[6]+f[1]*alpha_vdim[6]+alpha_vdim[2]*f[5]+f[2]*alpha_vdim[5]+alpha_vdim[3]*f[4]+f[3]*alpha_vdim[4]); 

  return cflFreq_mid; 
} 
