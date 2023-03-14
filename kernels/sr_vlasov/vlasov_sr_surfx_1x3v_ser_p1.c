#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_sr_surfx_1x3v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // p_over_gamma: p/gamma (velocity).
  // fl/fc/fr:  Input Distribution function in left/center/right cells.
  // out:       Incremented distribution function in center cell.
  const double *p0_over_gamma = &p_over_gamma[0]; 
  double alpha[8] = {0.0}; 
  alpha[0] = p0_over_gamma[0]; 
  alpha[1] = p0_over_gamma[1]; 
  alpha[2] = p0_over_gamma[2]; 
  alpha[3] = p0_over_gamma[3]; 
  alpha[4] = p0_over_gamma[4]; 
  alpha[5] = p0_over_gamma[5]; 
  alpha[6] = p0_over_gamma[6]; 
  alpha[7] = p0_over_gamma[7]; 

  const double dx10 = 2/dxv[0]; 
  const double dv = dxv[1], wv = w[1]; 
  double Ghat_r[8]; 
  double Ghat_l[8]; 
  if (wv>0) { 

  Ghat_r[0] = alpha[7]*(0.4330127018922193*fc[15]+0.25*fc[14])+0.4330127018922193*(alpha[6]*fc[13]+alpha[5]*fc[12]+alpha[4]*fc[11])+0.25*(alpha[6]*fc[10]+alpha[5]*fc[9])+0.4330127018922193*alpha[3]*fc[8]+0.25*alpha[4]*fc[7]+0.4330127018922193*(alpha[2]*fc[6]+alpha[1]*fc[5])+0.25*(alpha[3]*fc[4]+alpha[2]*fc[3]+alpha[1]*fc[2])+alpha[0]*(0.4330127018922193*fc[1]+0.25*fc[0]); 
  Ghat_r[1] = alpha[6]*(0.4330127018922193*fc[15]+0.25*fc[14])+0.4330127018922193*(alpha[7]*fc[13]+alpha[3]*fc[12]+alpha[2]*fc[11])+0.25*(alpha[7]*fc[10]+alpha[3]*fc[9])+0.4330127018922193*alpha[5]*fc[8]+0.25*alpha[2]*fc[7]+0.4330127018922193*(alpha[4]*fc[6]+alpha[0]*fc[5])+0.25*(fc[4]*alpha[5]+fc[3]*alpha[4]+alpha[0]*fc[2])+alpha[1]*(0.4330127018922193*fc[1]+0.25*fc[0]); 
  Ghat_r[2] = alpha[5]*(0.4330127018922193*fc[15]+0.25*fc[14])+0.4330127018922193*(alpha[3]*fc[13]+alpha[7]*fc[12]+alpha[1]*fc[11])+0.25*(alpha[3]*fc[10]+alpha[7]*fc[9])+0.4330127018922193*alpha[6]*fc[8]+0.25*alpha[1]*fc[7]+0.4330127018922193*alpha[0]*fc[6]+0.25*fc[4]*alpha[6]+0.4330127018922193*alpha[4]*fc[5]+0.25*(fc[2]*alpha[4]+alpha[0]*fc[3])+(0.4330127018922193*fc[1]+0.25*fc[0])*alpha[2]; 
  Ghat_r[3] = alpha[4]*(0.4330127018922193*fc[15]+0.25*fc[14])+0.4330127018922193*(alpha[2]*fc[13]+alpha[1]*fc[12]+alpha[7]*fc[11])+0.25*(alpha[2]*fc[10]+alpha[1]*fc[9])+0.4330127018922193*alpha[0]*fc[8]+0.25*alpha[7]*fc[7]+alpha[6]*(0.4330127018922193*fc[6]+0.25*fc[3])+0.4330127018922193*alpha[5]*fc[5]+0.25*(fc[2]*alpha[5]+alpha[0]*fc[4])+(0.4330127018922193*fc[1]+0.25*fc[0])*alpha[3]; 
  Ghat_r[4] = alpha[3]*(0.4330127018922193*fc[15]+0.25*fc[14])+0.4330127018922193*(alpha[5]*fc[13]+alpha[6]*fc[12]+alpha[0]*fc[11])+0.25*(alpha[5]*fc[10]+alpha[6]*fc[9])+0.4330127018922193*alpha[7]*fc[8]+0.25*(alpha[0]*fc[7]+fc[4]*alpha[7])+0.4330127018922193*(alpha[1]*fc[6]+alpha[2]*fc[5]+fc[1]*alpha[4])+0.25*(fc[0]*alpha[4]+alpha[1]*fc[3]+alpha[2]*fc[2]); 
  Ghat_r[5] = alpha[2]*(0.4330127018922193*fc[15]+0.25*fc[14])+0.4330127018922193*(alpha[4]*fc[13]+alpha[0]*fc[12]+alpha[6]*fc[11])+0.25*(alpha[4]*fc[10]+alpha[0]*fc[9])+0.4330127018922193*alpha[1]*fc[8]+0.25*alpha[6]*fc[7]+(0.4330127018922193*fc[6]+0.25*fc[3])*alpha[7]+0.4330127018922193*(alpha[3]*fc[5]+fc[1]*alpha[5])+0.25*(fc[0]*alpha[5]+alpha[1]*fc[4]+fc[2]*alpha[3]); 
  Ghat_r[6] = alpha[1]*(0.4330127018922193*fc[15]+0.25*fc[14])+0.4330127018922193*(alpha[0]*fc[13]+alpha[4]*fc[12]+alpha[5]*fc[11])+0.25*(alpha[0]*fc[10]+alpha[4]*fc[9])+0.4330127018922193*alpha[2]*fc[8]+0.25*alpha[5]*fc[7]+(0.4330127018922193*fc[5]+0.25*fc[2])*alpha[7]+0.4330127018922193*(alpha[3]*fc[6]+fc[1]*alpha[6])+0.25*(fc[0]*alpha[6]+alpha[2]*fc[4]+alpha[3]*fc[3]); 
  Ghat_r[7] = alpha[0]*(0.4330127018922193*fc[15]+0.25*fc[14])+0.4330127018922193*(alpha[1]*fc[13]+alpha[2]*fc[12]+alpha[3]*fc[11])+0.25*(alpha[1]*fc[10]+alpha[2]*fc[9])+0.4330127018922193*alpha[4]*fc[8]+0.25*alpha[3]*fc[7]+(0.4330127018922193*fc[1]+0.25*fc[0])*alpha[7]+0.4330127018922193*(alpha[5]*fc[6]+fc[5]*alpha[6])+0.25*(fc[2]*alpha[6]+fc[3]*alpha[5]+alpha[4]*fc[4]); 

  Ghat_l[0] = alpha[7]*(0.4330127018922193*fl[15]+0.25*fl[14])+0.4330127018922193*(alpha[6]*fl[13]+alpha[5]*fl[12]+alpha[4]*fl[11])+0.25*(alpha[6]*fl[10]+alpha[5]*fl[9])+0.4330127018922193*alpha[3]*fl[8]+0.25*alpha[4]*fl[7]+0.4330127018922193*(alpha[2]*fl[6]+alpha[1]*fl[5])+0.25*(alpha[3]*fl[4]+alpha[2]*fl[3]+alpha[1]*fl[2])+alpha[0]*(0.4330127018922193*fl[1]+0.25*fl[0]); 
  Ghat_l[1] = alpha[6]*(0.4330127018922193*fl[15]+0.25*fl[14])+0.4330127018922193*(alpha[7]*fl[13]+alpha[3]*fl[12]+alpha[2]*fl[11])+0.25*(alpha[7]*fl[10]+alpha[3]*fl[9])+0.4330127018922193*alpha[5]*fl[8]+0.25*alpha[2]*fl[7]+0.4330127018922193*(alpha[4]*fl[6]+alpha[0]*fl[5])+0.25*(fl[4]*alpha[5]+fl[3]*alpha[4]+alpha[0]*fl[2])+alpha[1]*(0.4330127018922193*fl[1]+0.25*fl[0]); 
  Ghat_l[2] = alpha[5]*(0.4330127018922193*fl[15]+0.25*fl[14])+0.4330127018922193*(alpha[3]*fl[13]+alpha[7]*fl[12]+alpha[1]*fl[11])+0.25*(alpha[3]*fl[10]+alpha[7]*fl[9])+0.4330127018922193*alpha[6]*fl[8]+0.25*alpha[1]*fl[7]+0.4330127018922193*alpha[0]*fl[6]+0.25*fl[4]*alpha[6]+0.4330127018922193*alpha[4]*fl[5]+0.25*(fl[2]*alpha[4]+alpha[0]*fl[3])+(0.4330127018922193*fl[1]+0.25*fl[0])*alpha[2]; 
  Ghat_l[3] = alpha[4]*(0.4330127018922193*fl[15]+0.25*fl[14])+0.4330127018922193*(alpha[2]*fl[13]+alpha[1]*fl[12]+alpha[7]*fl[11])+0.25*(alpha[2]*fl[10]+alpha[1]*fl[9])+0.4330127018922193*alpha[0]*fl[8]+0.25*alpha[7]*fl[7]+alpha[6]*(0.4330127018922193*fl[6]+0.25*fl[3])+0.4330127018922193*alpha[5]*fl[5]+0.25*(fl[2]*alpha[5]+alpha[0]*fl[4])+(0.4330127018922193*fl[1]+0.25*fl[0])*alpha[3]; 
  Ghat_l[4] = alpha[3]*(0.4330127018922193*fl[15]+0.25*fl[14])+0.4330127018922193*(alpha[5]*fl[13]+alpha[6]*fl[12]+alpha[0]*fl[11])+0.25*(alpha[5]*fl[10]+alpha[6]*fl[9])+0.4330127018922193*alpha[7]*fl[8]+0.25*(alpha[0]*fl[7]+fl[4]*alpha[7])+0.4330127018922193*(alpha[1]*fl[6]+alpha[2]*fl[5]+fl[1]*alpha[4])+0.25*(fl[0]*alpha[4]+alpha[1]*fl[3]+alpha[2]*fl[2]); 
  Ghat_l[5] = alpha[2]*(0.4330127018922193*fl[15]+0.25*fl[14])+0.4330127018922193*(alpha[4]*fl[13]+alpha[0]*fl[12]+alpha[6]*fl[11])+0.25*(alpha[4]*fl[10]+alpha[0]*fl[9])+0.4330127018922193*alpha[1]*fl[8]+0.25*alpha[6]*fl[7]+(0.4330127018922193*fl[6]+0.25*fl[3])*alpha[7]+0.4330127018922193*(alpha[3]*fl[5]+fl[1]*alpha[5])+0.25*(fl[0]*alpha[5]+alpha[1]*fl[4]+fl[2]*alpha[3]); 
  Ghat_l[6] = alpha[1]*(0.4330127018922193*fl[15]+0.25*fl[14])+0.4330127018922193*(alpha[0]*fl[13]+alpha[4]*fl[12]+alpha[5]*fl[11])+0.25*(alpha[0]*fl[10]+alpha[4]*fl[9])+0.4330127018922193*alpha[2]*fl[8]+0.25*alpha[5]*fl[7]+(0.4330127018922193*fl[5]+0.25*fl[2])*alpha[7]+0.4330127018922193*(alpha[3]*fl[6]+fl[1]*alpha[6])+0.25*(fl[0]*alpha[6]+alpha[2]*fl[4]+alpha[3]*fl[3]); 
  Ghat_l[7] = alpha[0]*(0.4330127018922193*fl[15]+0.25*fl[14])+0.4330127018922193*(alpha[1]*fl[13]+alpha[2]*fl[12]+alpha[3]*fl[11])+0.25*(alpha[1]*fl[10]+alpha[2]*fl[9])+0.4330127018922193*alpha[4]*fl[8]+0.25*alpha[3]*fl[7]+(0.4330127018922193*fl[1]+0.25*fl[0])*alpha[7]+0.4330127018922193*(alpha[5]*fl[6]+fl[5]*alpha[6])+0.25*(fl[2]*alpha[6]+fl[3]*alpha[5]+alpha[4]*fl[4]); 

  } else { 

  Ghat_r[0] = -0.25*(alpha[7]*(1.732050807568877*fr[15]-1.0*fr[14])+1.732050807568877*(alpha[6]*fr[13]+alpha[5]*fr[12]+alpha[4]*fr[11])-1.0*(alpha[6]*fr[10]+alpha[5]*fr[9])+1.732050807568877*alpha[3]*fr[8]-1.0*alpha[4]*fr[7]+1.732050807568877*(alpha[2]*fr[6]+alpha[1]*fr[5])-1.0*(alpha[3]*fr[4]+alpha[2]*fr[3]+alpha[1]*fr[2])+alpha[0]*(1.732050807568877*fr[1]-1.0*fr[0])); 
  Ghat_r[1] = -0.25*(alpha[6]*(1.732050807568877*fr[15]-1.0*fr[14])+1.732050807568877*(alpha[7]*fr[13]+alpha[3]*fr[12]+alpha[2]*fr[11])-1.0*(alpha[7]*fr[10]+alpha[3]*fr[9])+1.732050807568877*alpha[5]*fr[8]-1.0*alpha[2]*fr[7]+1.732050807568877*(alpha[4]*fr[6]+alpha[0]*fr[5])-1.0*(fr[4]*alpha[5]+fr[3]*alpha[4]+alpha[0]*fr[2])+alpha[1]*(1.732050807568877*fr[1]-1.0*fr[0])); 
  Ghat_r[2] = -0.25*(alpha[5]*(1.732050807568877*fr[15]-1.0*fr[14])+1.732050807568877*(alpha[3]*fr[13]+alpha[7]*fr[12]+alpha[1]*fr[11])-1.0*(alpha[3]*fr[10]+alpha[7]*fr[9])+1.732050807568877*alpha[6]*fr[8]-1.0*alpha[1]*fr[7]+1.732050807568877*alpha[0]*fr[6]-1.0*fr[4]*alpha[6]+1.732050807568877*alpha[4]*fr[5]-1.0*(fr[2]*alpha[4]+alpha[0]*fr[3])+(1.732050807568877*fr[1]-1.0*fr[0])*alpha[2]); 
  Ghat_r[3] = -0.25*(alpha[4]*(1.732050807568877*fr[15]-1.0*fr[14])+1.732050807568877*(alpha[2]*fr[13]+alpha[1]*fr[12]+alpha[7]*fr[11])-1.0*(alpha[2]*fr[10]+alpha[1]*fr[9])+1.732050807568877*alpha[0]*fr[8]-1.0*alpha[7]*fr[7]+alpha[6]*(1.732050807568877*fr[6]-1.0*fr[3])+1.732050807568877*alpha[5]*fr[5]-1.0*(fr[2]*alpha[5]+alpha[0]*fr[4])+(1.732050807568877*fr[1]-1.0*fr[0])*alpha[3]); 
  Ghat_r[4] = -0.25*(alpha[3]*(1.732050807568877*fr[15]-1.0*fr[14])+1.732050807568877*(alpha[5]*fr[13]+alpha[6]*fr[12]+alpha[0]*fr[11])-1.0*(alpha[5]*fr[10]+alpha[6]*fr[9])+1.732050807568877*alpha[7]*fr[8]-1.0*(alpha[0]*fr[7]+fr[4]*alpha[7])+1.732050807568877*(alpha[1]*fr[6]+alpha[2]*fr[5])+(1.732050807568877*fr[1]-1.0*fr[0])*alpha[4]-1.0*(alpha[1]*fr[3]+alpha[2]*fr[2])); 
  Ghat_r[5] = -0.25*(alpha[2]*(1.732050807568877*fr[15]-1.0*fr[14])+1.732050807568877*(alpha[4]*fr[13]+alpha[0]*fr[12]+alpha[6]*fr[11])-1.0*(alpha[4]*fr[10]+alpha[0]*fr[9])+1.732050807568877*alpha[1]*fr[8]-1.0*alpha[6]*fr[7]+(1.732050807568877*fr[6]-1.0*fr[3])*alpha[7]+1.732050807568877*alpha[3]*fr[5]+(1.732050807568877*fr[1]-1.0*fr[0])*alpha[5]-1.0*(alpha[1]*fr[4]+fr[2]*alpha[3])); 
  Ghat_r[6] = -0.25*(alpha[1]*(1.732050807568877*fr[15]-1.0*fr[14])+1.732050807568877*(alpha[0]*fr[13]+alpha[4]*fr[12]+alpha[5]*fr[11])-1.0*(alpha[0]*fr[10]+alpha[4]*fr[9])+1.732050807568877*alpha[2]*fr[8]-1.0*alpha[5]*fr[7]+(1.732050807568877*fr[5]-1.0*fr[2])*alpha[7]+1.732050807568877*alpha[3]*fr[6]+(1.732050807568877*fr[1]-1.0*fr[0])*alpha[6]-1.0*(alpha[2]*fr[4]+alpha[3]*fr[3])); 
  Ghat_r[7] = -0.25*(alpha[0]*(1.732050807568877*fr[15]-1.0*fr[14])+1.732050807568877*(alpha[1]*fr[13]+alpha[2]*fr[12]+alpha[3]*fr[11])-1.0*(alpha[1]*fr[10]+alpha[2]*fr[9])+1.732050807568877*alpha[4]*fr[8]-1.0*alpha[3]*fr[7]+(1.732050807568877*fr[1]-1.0*fr[0])*alpha[7]+1.732050807568877*alpha[5]*fr[6]+(1.732050807568877*fr[5]-1.0*fr[2])*alpha[6]-1.0*(fr[3]*alpha[5]+alpha[4]*fr[4])); 

  Ghat_l[0] = -0.25*(alpha[7]*(1.732050807568877*fc[15]-1.0*fc[14])+1.732050807568877*(alpha[6]*fc[13]+alpha[5]*fc[12]+alpha[4]*fc[11])-1.0*(alpha[6]*fc[10]+alpha[5]*fc[9])+1.732050807568877*alpha[3]*fc[8]-1.0*alpha[4]*fc[7]+1.732050807568877*(alpha[2]*fc[6]+alpha[1]*fc[5])-1.0*(alpha[3]*fc[4]+alpha[2]*fc[3]+alpha[1]*fc[2])+alpha[0]*(1.732050807568877*fc[1]-1.0*fc[0])); 
  Ghat_l[1] = -0.25*(alpha[6]*(1.732050807568877*fc[15]-1.0*fc[14])+1.732050807568877*(alpha[7]*fc[13]+alpha[3]*fc[12]+alpha[2]*fc[11])-1.0*(alpha[7]*fc[10]+alpha[3]*fc[9])+1.732050807568877*alpha[5]*fc[8]-1.0*alpha[2]*fc[7]+1.732050807568877*(alpha[4]*fc[6]+alpha[0]*fc[5])-1.0*(fc[4]*alpha[5]+fc[3]*alpha[4]+alpha[0]*fc[2])+alpha[1]*(1.732050807568877*fc[1]-1.0*fc[0])); 
  Ghat_l[2] = -0.25*(alpha[5]*(1.732050807568877*fc[15]-1.0*fc[14])+1.732050807568877*(alpha[3]*fc[13]+alpha[7]*fc[12]+alpha[1]*fc[11])-1.0*(alpha[3]*fc[10]+alpha[7]*fc[9])+1.732050807568877*alpha[6]*fc[8]-1.0*alpha[1]*fc[7]+1.732050807568877*alpha[0]*fc[6]-1.0*fc[4]*alpha[6]+1.732050807568877*alpha[4]*fc[5]-1.0*(fc[2]*alpha[4]+alpha[0]*fc[3])+(1.732050807568877*fc[1]-1.0*fc[0])*alpha[2]); 
  Ghat_l[3] = -0.25*(alpha[4]*(1.732050807568877*fc[15]-1.0*fc[14])+1.732050807568877*(alpha[2]*fc[13]+alpha[1]*fc[12]+alpha[7]*fc[11])-1.0*(alpha[2]*fc[10]+alpha[1]*fc[9])+1.732050807568877*alpha[0]*fc[8]-1.0*alpha[7]*fc[7]+alpha[6]*(1.732050807568877*fc[6]-1.0*fc[3])+1.732050807568877*alpha[5]*fc[5]-1.0*(fc[2]*alpha[5]+alpha[0]*fc[4])+(1.732050807568877*fc[1]-1.0*fc[0])*alpha[3]); 
  Ghat_l[4] = -0.25*(alpha[3]*(1.732050807568877*fc[15]-1.0*fc[14])+1.732050807568877*(alpha[5]*fc[13]+alpha[6]*fc[12]+alpha[0]*fc[11])-1.0*(alpha[5]*fc[10]+alpha[6]*fc[9])+1.732050807568877*alpha[7]*fc[8]-1.0*(alpha[0]*fc[7]+fc[4]*alpha[7])+1.732050807568877*(alpha[1]*fc[6]+alpha[2]*fc[5])+(1.732050807568877*fc[1]-1.0*fc[0])*alpha[4]-1.0*(alpha[1]*fc[3]+alpha[2]*fc[2])); 
  Ghat_l[5] = -0.25*(alpha[2]*(1.732050807568877*fc[15]-1.0*fc[14])+1.732050807568877*(alpha[4]*fc[13]+alpha[0]*fc[12]+alpha[6]*fc[11])-1.0*(alpha[4]*fc[10]+alpha[0]*fc[9])+1.732050807568877*alpha[1]*fc[8]-1.0*alpha[6]*fc[7]+(1.732050807568877*fc[6]-1.0*fc[3])*alpha[7]+1.732050807568877*alpha[3]*fc[5]+(1.732050807568877*fc[1]-1.0*fc[0])*alpha[5]-1.0*(alpha[1]*fc[4]+fc[2]*alpha[3])); 
  Ghat_l[6] = -0.25*(alpha[1]*(1.732050807568877*fc[15]-1.0*fc[14])+1.732050807568877*(alpha[0]*fc[13]+alpha[4]*fc[12]+alpha[5]*fc[11])-1.0*(alpha[0]*fc[10]+alpha[4]*fc[9])+1.732050807568877*alpha[2]*fc[8]-1.0*alpha[5]*fc[7]+(1.732050807568877*fc[5]-1.0*fc[2])*alpha[7]+1.732050807568877*alpha[3]*fc[6]+(1.732050807568877*fc[1]-1.0*fc[0])*alpha[6]-1.0*(alpha[2]*fc[4]+alpha[3]*fc[3])); 
  Ghat_l[7] = -0.25*(alpha[0]*(1.732050807568877*fc[15]-1.0*fc[14])+1.732050807568877*(alpha[1]*fc[13]+alpha[2]*fc[12]+alpha[3]*fc[11])-1.0*(alpha[1]*fc[10]+alpha[2]*fc[9])+1.732050807568877*alpha[4]*fc[8]-1.0*alpha[3]*fc[7]+(1.732050807568877*fc[1]-1.0*fc[0])*alpha[7]+1.732050807568877*alpha[5]*fc[6]+(1.732050807568877*fc[5]-1.0*fc[2])*alpha[6]-1.0*(fc[3]*alpha[5]+alpha[4]*fc[4])); 

  } 
  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dx10; 
  out[1] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dx10; 
  out[2] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dx10; 
  out[3] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dx10; 
  out[4] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dx10; 
  out[5] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dx10; 
  out[6] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dx10; 
  out[7] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*dx10; 
  out[8] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dx10; 
  out[9] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*dx10; 
  out[10] += (0.7071067811865475*Ghat_l[6]-0.7071067811865475*Ghat_r[6])*dx10; 
  out[11] += -1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dx10; 
  out[12] += -1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dx10; 
  out[13] += -1.224744871391589*(Ghat_r[6]+Ghat_l[6])*dx10; 
  out[14] += (0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*dx10; 
  out[15] += -1.224744871391589*(Ghat_r[7]+Ghat_l[7])*dx10; 
} 
