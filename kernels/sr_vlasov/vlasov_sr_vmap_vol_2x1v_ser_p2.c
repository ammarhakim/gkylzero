#include <gkyl_vlasov_sr_kernels.h> 
GKYL_CU_DH double vlasov_sr_vmap_vol_2x1v_ser_p2(const double *w, const double *dxv, const double *jacob_vel_inv, const double *gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out) 
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
  double alpha_vdim[20] = {0.0}; 

  cflFreq_mid += 5.0*fabs((0.3535533905932737*p0_over_gamma[0]-0.3952847075210473*p0_over_gamma[2])*dx10); 

  out[1] += 1.224744871391589*(p0_over_gamma[2]*f[9]+p0_over_gamma[1]*f[3]+f[0]*p0_over_gamma[0])*dx10; 
  out[4] += 1.224744871391589*(p0_over_gamma[2]*f[16]+p0_over_gamma[1]*f[6]+p0_over_gamma[0]*f[2])*dx10; 
  out[5] += (1.095445115010332*(p0_over_gamma[1]*f[9]+p0_over_gamma[2]*f[3])+1.224744871391589*(p0_over_gamma[0]*f[3]+f[0]*p0_over_gamma[1]))*dx10; 
  out[7] += 2.738612787525831*(p0_over_gamma[2]*f[15]+p0_over_gamma[1]*f[5]+p0_over_gamma[0]*f[1])*dx10; 
  out[10] += (1.095445115010332*(p0_over_gamma[1]*f[16]+p0_over_gamma[2]*f[6])+1.224744871391589*(p0_over_gamma[0]*f[6]+p0_over_gamma[1]*f[2]))*dx10; 
  out[11] += 2.738612787525831*(p0_over_gamma[2]*f[19]+p0_over_gamma[1]*f[10]+p0_over_gamma[0]*f[4])*dx10; 
  out[12] += 1.224744871391589*(p0_over_gamma[1]*f[14]+p0_over_gamma[0]*f[8])*dx10; 
  out[13] += (2.449489742783178*(p0_over_gamma[1]*f[15]+p0_over_gamma[2]*f[5])+2.738612787525831*(p0_over_gamma[0]*f[5]+f[1]*p0_over_gamma[1]))*dx10; 
  out[15] += ((0.7824607964359517*p0_over_gamma[2]+1.224744871391589*p0_over_gamma[0])*f[9]+1.095445115010332*p0_over_gamma[1]*f[3]+1.224744871391589*f[0]*p0_over_gamma[2])*dx10; 
  out[17] += (2.449489742783178*(p0_over_gamma[1]*f[19]+p0_over_gamma[2]*f[10])+2.738612787525831*(p0_over_gamma[0]*f[10]+p0_over_gamma[1]*f[4]))*dx10; 
  out[18] += (1.095445115010332*p0_over_gamma[2]*f[14]+1.224744871391589*(p0_over_gamma[0]*f[14]+p0_over_gamma[1]*f[8]))*dx10; 
  out[19] += ((0.7824607964359517*p0_over_gamma[2]+1.224744871391589*p0_over_gamma[0])*f[16]+1.095445115010332*p0_over_gamma[1]*f[6]+1.224744871391589*f[2]*p0_over_gamma[2])*dx10; 

  alpha_vdim[0] = E0[0]*jacob_vel_inv0[0]*dv10; 
  alpha_vdim[1] = jacob_vel_inv0[0]*E0[1]*dv10; 
  alpha_vdim[2] = jacob_vel_inv0[0]*E0[2]*dv10; 
  alpha_vdim[3] = E0[0]*jacob_vel_inv0[1]*dv10; 
  alpha_vdim[4] = jacob_vel_inv0[0]*E0[3]*dv10; 
  alpha_vdim[5] = E0[1]*jacob_vel_inv0[1]*dv10; 
  alpha_vdim[6] = jacob_vel_inv0[1]*E0[2]*dv10; 
  alpha_vdim[7] = jacob_vel_inv0[0]*E0[4]*dv10; 
  alpha_vdim[8] = jacob_vel_inv0[0]*E0[5]*dv10; 
  alpha_vdim[9] = E0[0]*jacob_vel_inv0[2]*dv10; 
  alpha_vdim[10] = jacob_vel_inv0[1]*E0[3]*dv10; 
  alpha_vdim[11] = jacob_vel_inv0[0]*E0[6]*dv10; 
  alpha_vdim[12] = jacob_vel_inv0[0]*E0[7]*dv10; 
  alpha_vdim[13] = jacob_vel_inv0[1]*E0[4]*dv10; 
  alpha_vdim[14] = jacob_vel_inv0[1]*E0[5]*dv10; 
  alpha_vdim[15] = E0[1]*jacob_vel_inv0[2]*dv10; 
  alpha_vdim[16] = E0[2]*jacob_vel_inv0[2]*dv10; 
  alpha_vdim[17] = jacob_vel_inv0[1]*E0[6]*dv10; 
  alpha_vdim[18] = jacob_vel_inv0[1]*E0[7]*dv10; 
  alpha_vdim[19] = jacob_vel_inv0[2]*E0[3]*dv10; 
  cflFreq_mid += 5.0*fabs(0.1767766952966368*alpha_vdim[0]-0.1976423537605236*(alpha_vdim[9]+alpha_vdim[8]+alpha_vdim[7])); 

  out[3] += 0.6123724356957944*(alpha_vdim[19]*f[19]+alpha_vdim[18]*f[18]+alpha_vdim[17]*f[17]+alpha_vdim[16]*f[16]+alpha_vdim[15]*f[15]+alpha_vdim[14]*f[14]+alpha_vdim[13]*f[13]+alpha_vdim[12]*f[12]+alpha_vdim[11]*f[11]+alpha_vdim[10]*f[10]+alpha_vdim[9]*f[9]+alpha_vdim[8]*f[8]+alpha_vdim[7]*f[7]+alpha_vdim[6]*f[6]+alpha_vdim[5]*f[5]+alpha_vdim[4]*f[4]+alpha_vdim[3]*f[3]+alpha_vdim[2]*f[2]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[5] += 0.6123724356957944*(alpha_vdim[16]*f[19]+f[16]*alpha_vdim[19]+alpha_vdim[14]*f[18]+f[14]*alpha_vdim[18])+0.5477225575051661*(alpha_vdim[10]*f[17]+f[10]*alpha_vdim[17])+0.6123724356957944*(alpha_vdim[9]*f[15]+f[9]*alpha_vdim[15])+0.5477225575051661*(alpha_vdim[5]*f[13]+f[5]*alpha_vdim[13])+0.6123724356957944*(alpha_vdim[8]*f[12]+f[8]*alpha_vdim[12])+0.5477225575051661*(alpha_vdim[4]*f[11]+f[4]*alpha_vdim[11])+0.6123724356957944*(alpha_vdim[6]*f[10]+f[6]*alpha_vdim[10])+0.5477225575051661*(alpha_vdim[1]*f[7]+f[1]*alpha_vdim[7])+0.6123724356957944*(alpha_vdim[3]*f[5]+f[3]*alpha_vdim[5]+alpha_vdim[2]*f[4]+f[2]*alpha_vdim[4]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[6] += 0.6123724356957944*(alpha_vdim[15]*f[19]+f[15]*alpha_vdim[19])+0.5477225575051661*(alpha_vdim[10]*f[18]+f[10]*alpha_vdim[18])+0.6123724356957944*(alpha_vdim[13]*f[17]+f[13]*alpha_vdim[17]+alpha_vdim[9]*f[16]+f[9]*alpha_vdim[16])+0.5477225575051661*(alpha_vdim[6]*f[14]+f[6]*alpha_vdim[14]+alpha_vdim[4]*f[12]+f[4]*alpha_vdim[12])+0.6123724356957944*(alpha_vdim[7]*f[11]+f[7]*alpha_vdim[11]+alpha_vdim[5]*f[10]+f[5]*alpha_vdim[10])+0.5477225575051661*(alpha_vdim[2]*f[8]+f[2]*alpha_vdim[8])+0.6123724356957944*(alpha_vdim[3]*f[6]+f[3]*alpha_vdim[6]+alpha_vdim[1]*f[4]+f[1]*alpha_vdim[4]+alpha_vdim[0]*f[2]+f[0]*alpha_vdim[2]); 
  out[9] += 1.224744871391589*(alpha_vdim[10]*f[19]+f[10]*alpha_vdim[19])+1.369306393762915*(alpha_vdim[12]*f[18]+f[12]*alpha_vdim[18]+alpha_vdim[11]*f[17]+f[11]*alpha_vdim[17])+1.224744871391589*(alpha_vdim[6]*f[16]+f[6]*alpha_vdim[16]+alpha_vdim[5]*f[15]+f[5]*alpha_vdim[15])+1.369306393762915*(alpha_vdim[8]*f[14]+f[8]*alpha_vdim[14]+alpha_vdim[7]*f[13]+f[7]*alpha_vdim[13]+alpha_vdim[4]*f[10]+f[4]*alpha_vdim[10])+1.224744871391589*(alpha_vdim[3]*f[9]+f[3]*alpha_vdim[9])+1.369306393762915*(alpha_vdim[2]*f[6]+f[2]*alpha_vdim[6]+alpha_vdim[1]*f[5]+f[1]*alpha_vdim[5]+alpha_vdim[0]*f[3]+f[0]*alpha_vdim[3]); 
  out[10] += 0.6123724356957944*(alpha_vdim[9]*f[19]+f[9]*alpha_vdim[19])+(0.4898979485566357*alpha_vdim[17]+0.5477225575051661*alpha_vdim[6])*f[18]+0.4898979485566357*f[17]*alpha_vdim[18]+0.5477225575051661*(f[6]*alpha_vdim[18]+alpha_vdim[5]*f[17]+f[5]*alpha_vdim[17])+0.6123724356957944*(alpha_vdim[15]*f[16]+f[15]*alpha_vdim[16])+0.5477225575051661*(alpha_vdim[10]*f[14]+f[10]*alpha_vdim[14]+alpha_vdim[10]*f[13]+f[10]*alpha_vdim[13])+(0.4898979485566357*alpha_vdim[11]+0.5477225575051661*alpha_vdim[2])*f[12]+0.4898979485566357*f[11]*alpha_vdim[12]+0.5477225575051661*(f[2]*alpha_vdim[12]+alpha_vdim[1]*f[11]+f[1]*alpha_vdim[11])+0.6123724356957944*(alpha_vdim[3]*f[10]+f[3]*alpha_vdim[10])+0.5477225575051661*(alpha_vdim[4]*f[8]+f[4]*alpha_vdim[8]+alpha_vdim[4]*f[7]+f[4]*alpha_vdim[7])+0.6123724356957944*(alpha_vdim[5]*f[6]+f[5]*alpha_vdim[6]+alpha_vdim[0]*f[4]+f[0]*alpha_vdim[4]+alpha_vdim[1]*f[2]+f[1]*alpha_vdim[2]); 
  out[13] += 0.5477225575051661*(alpha_vdim[19]*f[19]+alpha_vdim[18]*f[18])+0.3912303982179757*alpha_vdim[17]*f[17]+0.6123724356957944*(alpha_vdim[6]*f[17]+f[6]*alpha_vdim[17])+0.5477225575051661*alpha_vdim[15]*f[15]+0.3912303982179757*alpha_vdim[13]*f[13]+0.6123724356957944*(alpha_vdim[3]*f[13]+f[3]*alpha_vdim[13])+0.5477225575051661*alpha_vdim[12]*f[12]+0.3912303982179757*alpha_vdim[11]*f[11]+0.6123724356957944*(alpha_vdim[2]*f[11]+f[2]*alpha_vdim[11])+0.5477225575051661*alpha_vdim[10]*f[10]+0.3912303982179757*alpha_vdim[7]*f[7]+0.6123724356957944*(alpha_vdim[0]*f[7]+f[0]*alpha_vdim[7])+0.5477225575051661*(alpha_vdim[5]*f[5]+alpha_vdim[4]*f[4]+alpha_vdim[1]*f[1]); 
  out[14] += 0.5477225575051661*alpha_vdim[19]*f[19]+0.3912303982179757*alpha_vdim[18]*f[18]+0.6123724356957944*(alpha_vdim[5]*f[18]+f[5]*alpha_vdim[18])+0.5477225575051661*(alpha_vdim[17]*f[17]+alpha_vdim[16]*f[16])+0.3912303982179757*alpha_vdim[14]*f[14]+0.6123724356957944*(alpha_vdim[3]*f[14]+f[3]*alpha_vdim[14])+0.3912303982179757*alpha_vdim[12]*f[12]+0.6123724356957944*(alpha_vdim[1]*f[12]+f[1]*alpha_vdim[12])+0.5477225575051661*(alpha_vdim[11]*f[11]+alpha_vdim[10]*f[10])+0.3912303982179757*alpha_vdim[8]*f[8]+0.6123724356957944*(alpha_vdim[0]*f[8]+f[0]*alpha_vdim[8])+0.5477225575051661*(alpha_vdim[6]*f[6]+alpha_vdim[4]*f[4]+alpha_vdim[2]*f[2]); 
  out[15] += (1.095445115010332*alpha_vdim[17]+1.224744871391589*alpha_vdim[6])*f[19]+(1.095445115010332*f[17]+1.224744871391589*f[6])*alpha_vdim[19]+1.369306393762915*(alpha_vdim[8]*f[18]+f[8]*alpha_vdim[18])+1.224744871391589*(alpha_vdim[4]*f[17]+f[4]*alpha_vdim[17]+alpha_vdim[10]*f[16]+f[10]*alpha_vdim[16])+(1.095445115010332*alpha_vdim[13]+1.224744871391589*alpha_vdim[3])*f[15]+(1.095445115010332*f[13]+1.224744871391589*f[3])*alpha_vdim[15]+1.369306393762915*(alpha_vdim[12]*f[14]+f[12]*alpha_vdim[14])+1.224744871391589*(alpha_vdim[1]*f[13]+f[1]*alpha_vdim[13]+alpha_vdim[10]*f[11]+f[10]*alpha_vdim[11])+1.369306393762915*(alpha_vdim[2]*f[10]+f[2]*alpha_vdim[10])+1.224744871391589*(alpha_vdim[5]*f[9]+f[5]*alpha_vdim[9]+alpha_vdim[5]*f[7]+f[5]*alpha_vdim[7])+1.369306393762915*(alpha_vdim[4]*f[6]+f[4]*alpha_vdim[6]+alpha_vdim[0]*f[5]+f[0]*alpha_vdim[5]+alpha_vdim[1]*f[3]+f[1]*alpha_vdim[3]); 
  out[16] += (1.095445115010332*alpha_vdim[18]+1.224744871391589*alpha_vdim[5])*f[19]+1.095445115010332*f[18]*alpha_vdim[19]+1.224744871391589*(f[5]*alpha_vdim[19]+alpha_vdim[4]*f[18]+f[4]*alpha_vdim[18])+1.369306393762915*(alpha_vdim[7]*f[17]+f[7]*alpha_vdim[17])+(1.095445115010332*alpha_vdim[14]+1.224744871391589*alpha_vdim[3])*f[16]+1.095445115010332*f[14]*alpha_vdim[16]+1.224744871391589*(f[3]*alpha_vdim[16]+alpha_vdim[10]*f[15]+f[10]*alpha_vdim[15]+alpha_vdim[2]*f[14]+f[2]*alpha_vdim[14])+1.369306393762915*(alpha_vdim[11]*f[13]+f[11]*alpha_vdim[13])+1.224744871391589*(alpha_vdim[10]*f[12]+f[10]*alpha_vdim[12])+1.369306393762915*(alpha_vdim[1]*f[10]+f[1]*alpha_vdim[10])+1.224744871391589*(alpha_vdim[6]*f[9]+f[6]*alpha_vdim[9]+alpha_vdim[6]*f[8]+f[6]*alpha_vdim[8])+1.369306393762915*(alpha_vdim[0]*f[6]+f[0]*alpha_vdim[6]+alpha_vdim[4]*f[5]+f[4]*alpha_vdim[5]+alpha_vdim[2]*f[3]+f[2]*alpha_vdim[3]); 
  out[17] += 0.5477225575051661*(alpha_vdim[15]*f[19]+f[15]*alpha_vdim[19])+0.4898979485566357*(alpha_vdim[10]*f[18]+f[10]*alpha_vdim[18])+(0.5477225575051661*alpha_vdim[14]+0.3912303982179757*alpha_vdim[13]+0.6123724356957944*alpha_vdim[3])*f[17]+(0.5477225575051661*f[14]+0.3912303982179757*f[13])*alpha_vdim[17]+0.6123724356957944*(f[3]*alpha_vdim[17]+alpha_vdim[6]*f[13]+f[6]*alpha_vdim[13])+0.4898979485566356*(alpha_vdim[4]*f[12]+f[4]*alpha_vdim[12])+(0.5477225575051661*alpha_vdim[8]+0.3912303982179757*alpha_vdim[7]+0.6123724356957944*alpha_vdim[0])*f[11]+(0.5477225575051661*f[8]+0.3912303982179757*f[7]+0.6123724356957944*f[0])*alpha_vdim[11]+0.5477225575051661*(alpha_vdim[5]*f[10]+f[5]*alpha_vdim[10])+0.6123724356957944*(alpha_vdim[2]*f[7]+f[2]*alpha_vdim[7])+0.5477225575051661*(alpha_vdim[1]*f[4]+f[1]*alpha_vdim[4]); 
  out[18] += 0.5477225575051661*(alpha_vdim[16]*f[19]+f[16]*alpha_vdim[19])+(0.3912303982179757*alpha_vdim[14]+0.5477225575051661*alpha_vdim[13]+0.6123724356957944*alpha_vdim[3])*f[18]+(0.3912303982179757*f[14]+0.5477225575051661*f[13]+0.6123724356957944*f[3])*alpha_vdim[18]+0.4898979485566357*(alpha_vdim[10]*f[17]+f[10]*alpha_vdim[17])+0.6123724356957944*(alpha_vdim[5]*f[14]+f[5]*alpha_vdim[14])+(0.3912303982179757*alpha_vdim[8]+0.5477225575051661*alpha_vdim[7]+0.6123724356957944*alpha_vdim[0])*f[12]+(0.3912303982179757*f[8]+0.5477225575051661*f[7]+0.6123724356957944*f[0])*alpha_vdim[12]+0.4898979485566356*(alpha_vdim[4]*f[11]+f[4]*alpha_vdim[11])+0.5477225575051661*(alpha_vdim[6]*f[10]+f[6]*alpha_vdim[10])+0.6123724356957944*(alpha_vdim[1]*f[8]+f[1]*alpha_vdim[8])+0.5477225575051661*(alpha_vdim[2]*f[4]+f[2]*alpha_vdim[4]); 
  out[19] += (1.095445115010332*(alpha_vdim[14]+alpha_vdim[13])+1.224744871391589*alpha_vdim[3])*f[19]+(1.095445115010332*(f[14]+f[13])+1.224744871391589*f[3])*alpha_vdim[19]+(1.095445115010332*(alpha_vdim[16]+alpha_vdim[11])+1.224744871391589*alpha_vdim[2])*f[18]+(1.095445115010332*(f[16]+f[11])+1.224744871391589*f[2])*alpha_vdim[18]+(1.095445115010332*(alpha_vdim[15]+alpha_vdim[12])+1.224744871391589*alpha_vdim[1])*f[17]+1.095445115010332*(f[15]+f[12])*alpha_vdim[17]+1.224744871391589*(f[1]*alpha_vdim[17]+alpha_vdim[5]*f[16]+f[5]*alpha_vdim[16]+alpha_vdim[6]*f[15]+f[6]*alpha_vdim[15]+alpha_vdim[4]*f[14]+f[4]*alpha_vdim[14]+alpha_vdim[4]*f[13]+f[4]*alpha_vdim[13]+alpha_vdim[6]*f[12]+f[6]*alpha_vdim[12]+alpha_vdim[5]*f[11]+f[5]*alpha_vdim[11])+(1.224744871391589*(alpha_vdim[9]+alpha_vdim[8]+alpha_vdim[7])+1.369306393762915*alpha_vdim[0])*f[10]+1.224744871391589*(f[9]+f[8]+f[7])*alpha_vdim[10]+1.369306393762915*(f[0]*alpha_vdim[10]+alpha_vdim[1]*f[6]+f[1]*alpha_vdim[6]+alpha_vdim[2]*f[5]+f[2]*alpha_vdim[5]+alpha_vdim[3]*f[4]+f[3]*alpha_vdim[4]); 

  return cflFreq_mid; 
} 
