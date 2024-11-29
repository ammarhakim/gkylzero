#include <gkyl_vlasov_sr_kernels.h> 
GKYL_CU_DH double vlasov_sr_vmap_surfy_2x2v_ser_p1(const double *w, const double *dxv, const double *jacob_vel_inv, const double *gamma, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:             Cell-center coordinates.
  // dxv[NDIM]:           Cell spacing.
  // jacob_vel_inv[VDIM]: Inverse velocity space Jacobian in each direction.
  // gamma:               Particle Lorentz boost factor sqrt(1 + p^2).
  // fl/fc/fr:            Input Distribution function in left/center/right cells.
  // out:                 Incremented distribution function in center cell.
  const double dx11 = 2.0/dxv[1]; 
  const double dv11 = 2.0/dxv[3]; 
  const double wv = w[3]; 
  const double *jacob_vel_inv_dir = &jacob_vel_inv[3]; 

  double p_over_gamma[8] = {0.0}; 
  p_over_gamma[0] = 2.738612787525831*jacob_vel_inv_dir[1]*gamma[5]*dv11+1.224744871391589*jacob_vel_inv_dir[0]*gamma[2]*dv11; 
  p_over_gamma[1] = 2.738612787525831*jacob_vel_inv_dir[1]*gamma[7]*dv11+1.224744871391589*jacob_vel_inv_dir[0]*gamma[3]*dv11; 
  p_over_gamma[2] = 2.449489742783178*jacob_vel_inv_dir[2]*gamma[5]*dv11+2.738612787525831*jacob_vel_inv_dir[0]*gamma[5]*dv11+1.224744871391589*jacob_vel_inv_dir[1]*gamma[2]*dv11; 
  p_over_gamma[3] = 2.449489742783178*jacob_vel_inv_dir[2]*gamma[7]*dv11+2.738612787525831*jacob_vel_inv_dir[0]*gamma[7]*dv11+1.224744871391589*jacob_vel_inv_dir[1]*gamma[3]*dv11; 
  p_over_gamma[4] = 1.224744871391589*jacob_vel_inv_dir[0]*gamma[6]*dv11; 
  p_over_gamma[5] = 2.449489742783178*jacob_vel_inv_dir[1]*gamma[5]*dv11+1.224744871391589*jacob_vel_inv_dir[2]*gamma[2]*dv11; 
  p_over_gamma[6] = 1.224744871391589*jacob_vel_inv_dir[1]*gamma[6]*dv11; 
  p_over_gamma[7] = 2.449489742783178*jacob_vel_inv_dir[1]*gamma[7]*dv11+1.224744871391589*jacob_vel_inv_dir[2]*gamma[3]*dv11; 
  double fUpwind_l[16]; 
  double fUpwind_r[16]; 
  if (wv>0) { 

  fUpwind_l[0] = 1.224744871391589*fl[2]+0.7071067811865475*fl[0]; 
  fUpwind_l[1] = 1.224744871391589*fl[5]+0.7071067811865475*fl[1]; 
  fUpwind_l[2] = 1.224744871391589*fl[7]+0.7071067811865475*fl[3]; 
  fUpwind_l[3] = 1.224744871391589*fl[9]+0.7071067811865475*fl[4]; 
  fUpwind_l[4] = 1.224744871391589*fl[11]+0.7071067811865475*fl[6]; 
  fUpwind_l[5] = 1.224744871391589*fl[12]+0.7071067811865475*fl[8]; 
  fUpwind_l[6] = 1.224744871391589*fl[14]+0.7071067811865475*fl[10]; 
  fUpwind_l[7] = 1.224744871391589*fl[15]+0.7071067811865475*fl[13]; 
  fUpwind_l[8] = 1.224744871391589*fl[18]+0.7071067811865475*fl[16]; 
  fUpwind_l[9] = 1.224744871391589*fl[20]+0.7071067811865475*fl[17]; 
  fUpwind_l[10] = 1.224744871391589*fl[22]+0.7071067811865475*fl[19]; 
  fUpwind_l[11] = 1.224744871391589*fl[23]+0.7071067811865475*fl[21]; 
  fUpwind_l[12] = 1.224744871391589*fl[26]+0.7071067811865475*fl[24]; 
  fUpwind_l[13] = 1.224744871391589*fl[28]+0.7071067811865475*fl[25]; 
  fUpwind_l[14] = 1.224744871391589*fl[30]+0.7071067811865475*fl[27]; 
  fUpwind_l[15] = 1.224744871391589*fl[31]+0.7071067811865475*fl[29]; 

  fUpwind_r[0] = 1.224744871391589*fc[2]+0.7071067811865475*fc[0]; 
  fUpwind_r[1] = 1.224744871391589*fc[5]+0.7071067811865475*fc[1]; 
  fUpwind_r[2] = 1.224744871391589*fc[7]+0.7071067811865475*fc[3]; 
  fUpwind_r[3] = 1.224744871391589*fc[9]+0.7071067811865475*fc[4]; 
  fUpwind_r[4] = 1.224744871391589*fc[11]+0.7071067811865475*fc[6]; 
  fUpwind_r[5] = 1.224744871391589*fc[12]+0.7071067811865475*fc[8]; 
  fUpwind_r[6] = 1.224744871391589*fc[14]+0.7071067811865475*fc[10]; 
  fUpwind_r[7] = 1.224744871391589*fc[15]+0.7071067811865475*fc[13]; 
  fUpwind_r[8] = 1.224744871391589*fc[18]+0.7071067811865475*fc[16]; 
  fUpwind_r[9] = 1.224744871391589*fc[20]+0.7071067811865475*fc[17]; 
  fUpwind_r[10] = 1.224744871391589*fc[22]+0.7071067811865475*fc[19]; 
  fUpwind_r[11] = 1.224744871391589*fc[23]+0.7071067811865475*fc[21]; 
  fUpwind_r[12] = 1.224744871391589*fc[26]+0.7071067811865475*fc[24]; 
  fUpwind_r[13] = 1.224744871391589*fc[28]+0.7071067811865475*fc[25]; 
  fUpwind_r[14] = 1.224744871391589*fc[30]+0.7071067811865475*fc[27]; 
  fUpwind_r[15] = 1.224744871391589*fc[31]+0.7071067811865475*fc[29]; 

  } else { 

  fUpwind_l[0] = 0.7071067811865475*fc[0]-1.224744871391589*fc[2]; 
  fUpwind_l[1] = 0.7071067811865475*fc[1]-1.224744871391589*fc[5]; 
  fUpwind_l[2] = 0.7071067811865475*fc[3]-1.224744871391589*fc[7]; 
  fUpwind_l[3] = 0.7071067811865475*fc[4]-1.224744871391589*fc[9]; 
  fUpwind_l[4] = 0.7071067811865475*fc[6]-1.224744871391589*fc[11]; 
  fUpwind_l[5] = 0.7071067811865475*fc[8]-1.224744871391589*fc[12]; 
  fUpwind_l[6] = 0.7071067811865475*fc[10]-1.224744871391589*fc[14]; 
  fUpwind_l[7] = 0.7071067811865475*fc[13]-1.224744871391589*fc[15]; 
  fUpwind_l[8] = 0.7071067811865475*fc[16]-1.224744871391589*fc[18]; 
  fUpwind_l[9] = 0.7071067811865475*fc[17]-1.224744871391589*fc[20]; 
  fUpwind_l[10] = 0.7071067811865475*fc[19]-1.224744871391589*fc[22]; 
  fUpwind_l[11] = 0.7071067811865475*fc[21]-1.224744871391589*fc[23]; 
  fUpwind_l[12] = 0.7071067811865475*fc[24]-1.224744871391589*fc[26]; 
  fUpwind_l[13] = 0.7071067811865475*fc[25]-1.224744871391589*fc[28]; 
  fUpwind_l[14] = 0.7071067811865475*fc[27]-1.224744871391589*fc[30]; 
  fUpwind_l[15] = 0.7071067811865475*fc[29]-1.224744871391589*fc[31]; 

  fUpwind_r[0] = 0.7071067811865475*fr[0]-1.224744871391589*fr[2]; 
  fUpwind_r[1] = 0.7071067811865475*fr[1]-1.224744871391589*fr[5]; 
  fUpwind_r[2] = 0.7071067811865475*fr[3]-1.224744871391589*fr[7]; 
  fUpwind_r[3] = 0.7071067811865475*fr[4]-1.224744871391589*fr[9]; 
  fUpwind_r[4] = 0.7071067811865475*fr[6]-1.224744871391589*fr[11]; 
  fUpwind_r[5] = 0.7071067811865475*fr[8]-1.224744871391589*fr[12]; 
  fUpwind_r[6] = 0.7071067811865475*fr[10]-1.224744871391589*fr[14]; 
  fUpwind_r[7] = 0.7071067811865475*fr[13]-1.224744871391589*fr[15]; 
  fUpwind_r[8] = 0.7071067811865475*fr[16]-1.224744871391589*fr[18]; 
  fUpwind_r[9] = 0.7071067811865475*fr[17]-1.224744871391589*fr[20]; 
  fUpwind_r[10] = 0.7071067811865475*fr[19]-1.224744871391589*fr[22]; 
  fUpwind_r[11] = 0.7071067811865475*fr[21]-1.224744871391589*fr[23]; 
  fUpwind_r[12] = 0.7071067811865475*fr[24]-1.224744871391589*fr[26]; 
  fUpwind_r[13] = 0.7071067811865475*fr[25]-1.224744871391589*fr[28]; 
  fUpwind_r[14] = 0.7071067811865475*fr[27]-1.224744871391589*fr[30]; 
  fUpwind_r[15] = 0.7071067811865475*fr[29]-1.224744871391589*fr[31]; 

  } 
  double Ghat_l[16]; 
  double Ghat_r[16]; 
  Ghat_l[0] = 0.5*(p_over_gamma[7]*fUpwind_l[14]+p_over_gamma[5]*fUpwind_l[12]+p_over_gamma[6]*fUpwind_l[10]+p_over_gamma[4]*fUpwind_l[8]+p_over_gamma[3]*fUpwind_l[6]+p_over_gamma[2]*fUpwind_l[3]+p_over_gamma[1]*fUpwind_l[2]+fUpwind_l[0]*p_over_gamma[0]); 
  Ghat_l[1] = 0.5000000000000001*(p_over_gamma[7]*fUpwind_l[15]+p_over_gamma[5]*fUpwind_l[13]+p_over_gamma[6]*fUpwind_l[11]+p_over_gamma[4]*fUpwind_l[9])+0.5*(p_over_gamma[3]*fUpwind_l[7]+p_over_gamma[2]*fUpwind_l[5]+p_over_gamma[1]*fUpwind_l[4]+p_over_gamma[0]*fUpwind_l[1]); 
  Ghat_l[2] = 0.5000000000000001*(p_over_gamma[5]*fUpwind_l[14]+p_over_gamma[7]*fUpwind_l[12])+0.447213595499958*p_over_gamma[3]*fUpwind_l[10]+0.4472135954999579*p_over_gamma[1]*fUpwind_l[8]+fUpwind_l[6]*(0.447213595499958*p_over_gamma[6]+0.5*p_over_gamma[2])+0.4472135954999579*fUpwind_l[2]*p_over_gamma[4]+0.5*(fUpwind_l[3]*p_over_gamma[3]+p_over_gamma[0]*fUpwind_l[2]+fUpwind_l[0]*p_over_gamma[1]); 
  Ghat_l[3] = 0.447213595499958*p_over_gamma[3]*fUpwind_l[14]+0.4472135954999579*p_over_gamma[2]*fUpwind_l[12]+0.5000000000000001*(p_over_gamma[4]*fUpwind_l[10]+p_over_gamma[6]*fUpwind_l[8])+fUpwind_l[6]*(0.447213595499958*p_over_gamma[7]+0.5*p_over_gamma[1])+0.4472135954999579*fUpwind_l[3]*p_over_gamma[5]+0.5*(fUpwind_l[2]*p_over_gamma[3]+p_over_gamma[0]*fUpwind_l[3]+fUpwind_l[0]*p_over_gamma[2]); 
  Ghat_l[4] = 0.5*(p_over_gamma[5]*fUpwind_l[15]+p_over_gamma[7]*fUpwind_l[13])+0.4472135954999579*p_over_gamma[3]*fUpwind_l[11]+0.447213595499958*(p_over_gamma[1]*fUpwind_l[9]+p_over_gamma[6]*fUpwind_l[7])+0.5*(p_over_gamma[2]*fUpwind_l[7]+p_over_gamma[3]*fUpwind_l[5])+0.4472135954999579*fUpwind_l[4]*p_over_gamma[4]+0.5*(p_over_gamma[0]*fUpwind_l[4]+fUpwind_l[1]*p_over_gamma[1]); 
  Ghat_l[5] = 0.4472135954999579*p_over_gamma[3]*fUpwind_l[15]+0.447213595499958*p_over_gamma[2]*fUpwind_l[13]+0.5*(p_over_gamma[4]*fUpwind_l[11]+p_over_gamma[6]*fUpwind_l[9])+fUpwind_l[7]*(0.447213595499958*p_over_gamma[7]+0.5*p_over_gamma[1])+0.4472135954999579*fUpwind_l[5]*p_over_gamma[5]+0.5*(p_over_gamma[0]*fUpwind_l[5]+p_over_gamma[3]*fUpwind_l[4]+fUpwind_l[1]*p_over_gamma[2]); 
  Ghat_l[6] = (0.4*p_over_gamma[6]+0.447213595499958*p_over_gamma[2])*fUpwind_l[14]+0.4472135954999579*p_over_gamma[3]*fUpwind_l[12]+(0.4*p_over_gamma[7]+0.447213595499958*p_over_gamma[1])*fUpwind_l[10]+0.4472135954999579*p_over_gamma[3]*fUpwind_l[8]+0.447213595499958*(fUpwind_l[3]*p_over_gamma[7]+fUpwind_l[2]*p_over_gamma[6])+0.4472135954999579*(p_over_gamma[5]+p_over_gamma[4])*fUpwind_l[6]+0.5*(p_over_gamma[0]*fUpwind_l[6]+fUpwind_l[0]*p_over_gamma[3]+p_over_gamma[1]*fUpwind_l[3]+fUpwind_l[2]*p_over_gamma[2]); 
  Ghat_l[7] = (0.4*p_over_gamma[6]+0.4472135954999579*p_over_gamma[2])*fUpwind_l[15]+0.447213595499958*p_over_gamma[3]*fUpwind_l[13]+(0.4*p_over_gamma[7]+0.4472135954999579*p_over_gamma[1])*fUpwind_l[11]+0.447213595499958*(p_over_gamma[3]*fUpwind_l[9]+fUpwind_l[5]*p_over_gamma[7])+(0.4472135954999579*(p_over_gamma[5]+p_over_gamma[4])+0.5*p_over_gamma[0])*fUpwind_l[7]+0.447213595499958*fUpwind_l[4]*p_over_gamma[6]+0.5*(p_over_gamma[1]*fUpwind_l[5]+p_over_gamma[2]*fUpwind_l[4]+fUpwind_l[1]*p_over_gamma[3]); 
  Ghat_l[8] = 0.4472135954999579*p_over_gamma[7]*fUpwind_l[14]+(0.31943828249997*p_over_gamma[6]+0.5000000000000001*p_over_gamma[2])*fUpwind_l[10]+(0.31943828249997*p_over_gamma[4]+0.5*p_over_gamma[0])*fUpwind_l[8]+0.5000000000000001*fUpwind_l[3]*p_over_gamma[6]+0.4472135954999579*p_over_gamma[3]*fUpwind_l[6]+0.5*fUpwind_l[0]*p_over_gamma[4]+0.4472135954999579*p_over_gamma[1]*fUpwind_l[2]; 
  Ghat_l[9] = 0.4472135954999579*p_over_gamma[7]*fUpwind_l[15]+(0.31943828249997*p_over_gamma[6]+0.5000000000000001*p_over_gamma[2])*fUpwind_l[11]+(0.31943828249997*p_over_gamma[4]+0.5*p_over_gamma[0])*fUpwind_l[9]+0.447213595499958*p_over_gamma[3]*fUpwind_l[7]+0.5*fUpwind_l[5]*p_over_gamma[6]+0.5000000000000001*fUpwind_l[1]*p_over_gamma[4]+0.447213595499958*p_over_gamma[1]*fUpwind_l[4]; 
  Ghat_l[10] = 0.4*p_over_gamma[3]*fUpwind_l[14]+0.4472135954999579*p_over_gamma[6]*fUpwind_l[12]+(0.4472135954999579*p_over_gamma[5]+0.31943828249997*p_over_gamma[4]+0.5*p_over_gamma[0])*fUpwind_l[10]+(0.31943828249997*p_over_gamma[6]+0.5000000000000001*p_over_gamma[2])*fUpwind_l[8]+0.4*fUpwind_l[6]*p_over_gamma[7]+0.5*fUpwind_l[0]*p_over_gamma[6]+0.447213595499958*p_over_gamma[1]*fUpwind_l[6]+0.5000000000000001*fUpwind_l[3]*p_over_gamma[4]+0.447213595499958*fUpwind_l[2]*p_over_gamma[3]; 
  Ghat_l[11] = 0.4*p_over_gamma[3]*fUpwind_l[15]+0.4472135954999579*p_over_gamma[6]*fUpwind_l[13]+(0.4472135954999579*p_over_gamma[5]+0.31943828249997*p_over_gamma[4]+0.5*p_over_gamma[0])*fUpwind_l[11]+(0.31943828249997*p_over_gamma[6]+0.5000000000000001*p_over_gamma[2])*fUpwind_l[9]+fUpwind_l[7]*(0.4*p_over_gamma[7]+0.4472135954999579*p_over_gamma[1])+0.5000000000000001*fUpwind_l[1]*p_over_gamma[6]+0.5*p_over_gamma[4]*fUpwind_l[5]+0.4472135954999579*p_over_gamma[3]*fUpwind_l[4]; 
  Ghat_l[12] = (0.31943828249997*p_over_gamma[7]+0.5000000000000001*p_over_gamma[1])*fUpwind_l[14]+(0.31943828249997*p_over_gamma[5]+0.5*p_over_gamma[0])*fUpwind_l[12]+0.4472135954999579*p_over_gamma[6]*fUpwind_l[10]+0.5000000000000001*fUpwind_l[2]*p_over_gamma[7]+0.4472135954999579*p_over_gamma[3]*fUpwind_l[6]+0.5*fUpwind_l[0]*p_over_gamma[5]+0.4472135954999579*p_over_gamma[2]*fUpwind_l[3]; 
  Ghat_l[13] = (0.31943828249997*p_over_gamma[7]+0.5000000000000001*p_over_gamma[1])*fUpwind_l[15]+(0.31943828249997*p_over_gamma[5]+0.5*p_over_gamma[0])*fUpwind_l[13]+0.4472135954999579*p_over_gamma[6]*fUpwind_l[11]+0.5*fUpwind_l[4]*p_over_gamma[7]+0.447213595499958*p_over_gamma[3]*fUpwind_l[7]+0.5000000000000001*fUpwind_l[1]*p_over_gamma[5]+0.447213595499958*p_over_gamma[2]*fUpwind_l[5]; 
  Ghat_l[14] = (0.31943828249997*p_over_gamma[5]+0.4472135954999579*p_over_gamma[4]+0.5*p_over_gamma[0])*fUpwind_l[14]+(0.31943828249997*p_over_gamma[7]+0.5000000000000001*p_over_gamma[1])*fUpwind_l[12]+0.4*p_over_gamma[3]*fUpwind_l[10]+p_over_gamma[7]*(0.4472135954999579*fUpwind_l[8]+0.5*fUpwind_l[0])+fUpwind_l[6]*(0.4*p_over_gamma[6]+0.447213595499958*p_over_gamma[2])+0.5000000000000001*fUpwind_l[2]*p_over_gamma[5]+0.447213595499958*fUpwind_l[3]*p_over_gamma[3]; 
  Ghat_l[15] = (0.31943828249997*p_over_gamma[5]+0.4472135954999579*p_over_gamma[4]+0.5*p_over_gamma[0])*fUpwind_l[15]+(0.31943828249997*p_over_gamma[7]+0.5000000000000001*p_over_gamma[1])*fUpwind_l[13]+0.4*p_over_gamma[3]*fUpwind_l[11]+p_over_gamma[7]*(0.4472135954999579*fUpwind_l[9]+0.5000000000000001*fUpwind_l[1])+(0.4*p_over_gamma[6]+0.4472135954999579*p_over_gamma[2])*fUpwind_l[7]+0.5*fUpwind_l[4]*p_over_gamma[5]+0.4472135954999579*p_over_gamma[3]*fUpwind_l[5]; 

  Ghat_r[0] = 0.5*(p_over_gamma[7]*fUpwind_r[14]+p_over_gamma[5]*fUpwind_r[12]+p_over_gamma[6]*fUpwind_r[10]+p_over_gamma[4]*fUpwind_r[8]+p_over_gamma[3]*fUpwind_r[6]+p_over_gamma[2]*fUpwind_r[3]+p_over_gamma[1]*fUpwind_r[2]+fUpwind_r[0]*p_over_gamma[0]); 
  Ghat_r[1] = 0.5000000000000001*(p_over_gamma[7]*fUpwind_r[15]+p_over_gamma[5]*fUpwind_r[13]+p_over_gamma[6]*fUpwind_r[11]+p_over_gamma[4]*fUpwind_r[9])+0.5*(p_over_gamma[3]*fUpwind_r[7]+p_over_gamma[2]*fUpwind_r[5]+p_over_gamma[1]*fUpwind_r[4]+p_over_gamma[0]*fUpwind_r[1]); 
  Ghat_r[2] = 0.5000000000000001*(p_over_gamma[5]*fUpwind_r[14]+p_over_gamma[7]*fUpwind_r[12])+0.447213595499958*p_over_gamma[3]*fUpwind_r[10]+0.4472135954999579*p_over_gamma[1]*fUpwind_r[8]+fUpwind_r[6]*(0.447213595499958*p_over_gamma[6]+0.5*p_over_gamma[2])+0.4472135954999579*fUpwind_r[2]*p_over_gamma[4]+0.5*(fUpwind_r[3]*p_over_gamma[3]+p_over_gamma[0]*fUpwind_r[2]+fUpwind_r[0]*p_over_gamma[1]); 
  Ghat_r[3] = 0.447213595499958*p_over_gamma[3]*fUpwind_r[14]+0.4472135954999579*p_over_gamma[2]*fUpwind_r[12]+0.5000000000000001*(p_over_gamma[4]*fUpwind_r[10]+p_over_gamma[6]*fUpwind_r[8])+fUpwind_r[6]*(0.447213595499958*p_over_gamma[7]+0.5*p_over_gamma[1])+0.4472135954999579*fUpwind_r[3]*p_over_gamma[5]+0.5*(fUpwind_r[2]*p_over_gamma[3]+p_over_gamma[0]*fUpwind_r[3]+fUpwind_r[0]*p_over_gamma[2]); 
  Ghat_r[4] = 0.5*(p_over_gamma[5]*fUpwind_r[15]+p_over_gamma[7]*fUpwind_r[13])+0.4472135954999579*p_over_gamma[3]*fUpwind_r[11]+0.447213595499958*(p_over_gamma[1]*fUpwind_r[9]+p_over_gamma[6]*fUpwind_r[7])+0.5*(p_over_gamma[2]*fUpwind_r[7]+p_over_gamma[3]*fUpwind_r[5])+0.4472135954999579*fUpwind_r[4]*p_over_gamma[4]+0.5*(p_over_gamma[0]*fUpwind_r[4]+fUpwind_r[1]*p_over_gamma[1]); 
  Ghat_r[5] = 0.4472135954999579*p_over_gamma[3]*fUpwind_r[15]+0.447213595499958*p_over_gamma[2]*fUpwind_r[13]+0.5*(p_over_gamma[4]*fUpwind_r[11]+p_over_gamma[6]*fUpwind_r[9])+fUpwind_r[7]*(0.447213595499958*p_over_gamma[7]+0.5*p_over_gamma[1])+0.4472135954999579*fUpwind_r[5]*p_over_gamma[5]+0.5*(p_over_gamma[0]*fUpwind_r[5]+p_over_gamma[3]*fUpwind_r[4]+fUpwind_r[1]*p_over_gamma[2]); 
  Ghat_r[6] = (0.4*p_over_gamma[6]+0.447213595499958*p_over_gamma[2])*fUpwind_r[14]+0.4472135954999579*p_over_gamma[3]*fUpwind_r[12]+(0.4*p_over_gamma[7]+0.447213595499958*p_over_gamma[1])*fUpwind_r[10]+0.4472135954999579*p_over_gamma[3]*fUpwind_r[8]+0.447213595499958*(fUpwind_r[3]*p_over_gamma[7]+fUpwind_r[2]*p_over_gamma[6])+0.4472135954999579*(p_over_gamma[5]+p_over_gamma[4])*fUpwind_r[6]+0.5*(p_over_gamma[0]*fUpwind_r[6]+fUpwind_r[0]*p_over_gamma[3]+p_over_gamma[1]*fUpwind_r[3]+fUpwind_r[2]*p_over_gamma[2]); 
  Ghat_r[7] = (0.4*p_over_gamma[6]+0.4472135954999579*p_over_gamma[2])*fUpwind_r[15]+0.447213595499958*p_over_gamma[3]*fUpwind_r[13]+(0.4*p_over_gamma[7]+0.4472135954999579*p_over_gamma[1])*fUpwind_r[11]+0.447213595499958*(p_over_gamma[3]*fUpwind_r[9]+fUpwind_r[5]*p_over_gamma[7])+(0.4472135954999579*(p_over_gamma[5]+p_over_gamma[4])+0.5*p_over_gamma[0])*fUpwind_r[7]+0.447213595499958*fUpwind_r[4]*p_over_gamma[6]+0.5*(p_over_gamma[1]*fUpwind_r[5]+p_over_gamma[2]*fUpwind_r[4]+fUpwind_r[1]*p_over_gamma[3]); 
  Ghat_r[8] = 0.4472135954999579*p_over_gamma[7]*fUpwind_r[14]+(0.31943828249997*p_over_gamma[6]+0.5000000000000001*p_over_gamma[2])*fUpwind_r[10]+(0.31943828249997*p_over_gamma[4]+0.5*p_over_gamma[0])*fUpwind_r[8]+0.5000000000000001*fUpwind_r[3]*p_over_gamma[6]+0.4472135954999579*p_over_gamma[3]*fUpwind_r[6]+0.5*fUpwind_r[0]*p_over_gamma[4]+0.4472135954999579*p_over_gamma[1]*fUpwind_r[2]; 
  Ghat_r[9] = 0.4472135954999579*p_over_gamma[7]*fUpwind_r[15]+(0.31943828249997*p_over_gamma[6]+0.5000000000000001*p_over_gamma[2])*fUpwind_r[11]+(0.31943828249997*p_over_gamma[4]+0.5*p_over_gamma[0])*fUpwind_r[9]+0.447213595499958*p_over_gamma[3]*fUpwind_r[7]+0.5*fUpwind_r[5]*p_over_gamma[6]+0.5000000000000001*fUpwind_r[1]*p_over_gamma[4]+0.447213595499958*p_over_gamma[1]*fUpwind_r[4]; 
  Ghat_r[10] = 0.4*p_over_gamma[3]*fUpwind_r[14]+0.4472135954999579*p_over_gamma[6]*fUpwind_r[12]+(0.4472135954999579*p_over_gamma[5]+0.31943828249997*p_over_gamma[4]+0.5*p_over_gamma[0])*fUpwind_r[10]+(0.31943828249997*p_over_gamma[6]+0.5000000000000001*p_over_gamma[2])*fUpwind_r[8]+0.4*fUpwind_r[6]*p_over_gamma[7]+0.5*fUpwind_r[0]*p_over_gamma[6]+0.447213595499958*p_over_gamma[1]*fUpwind_r[6]+0.5000000000000001*fUpwind_r[3]*p_over_gamma[4]+0.447213595499958*fUpwind_r[2]*p_over_gamma[3]; 
  Ghat_r[11] = 0.4*p_over_gamma[3]*fUpwind_r[15]+0.4472135954999579*p_over_gamma[6]*fUpwind_r[13]+(0.4472135954999579*p_over_gamma[5]+0.31943828249997*p_over_gamma[4]+0.5*p_over_gamma[0])*fUpwind_r[11]+(0.31943828249997*p_over_gamma[6]+0.5000000000000001*p_over_gamma[2])*fUpwind_r[9]+fUpwind_r[7]*(0.4*p_over_gamma[7]+0.4472135954999579*p_over_gamma[1])+0.5000000000000001*fUpwind_r[1]*p_over_gamma[6]+0.5*p_over_gamma[4]*fUpwind_r[5]+0.4472135954999579*p_over_gamma[3]*fUpwind_r[4]; 
  Ghat_r[12] = (0.31943828249997*p_over_gamma[7]+0.5000000000000001*p_over_gamma[1])*fUpwind_r[14]+(0.31943828249997*p_over_gamma[5]+0.5*p_over_gamma[0])*fUpwind_r[12]+0.4472135954999579*p_over_gamma[6]*fUpwind_r[10]+0.5000000000000001*fUpwind_r[2]*p_over_gamma[7]+0.4472135954999579*p_over_gamma[3]*fUpwind_r[6]+0.5*fUpwind_r[0]*p_over_gamma[5]+0.4472135954999579*p_over_gamma[2]*fUpwind_r[3]; 
  Ghat_r[13] = (0.31943828249997*p_over_gamma[7]+0.5000000000000001*p_over_gamma[1])*fUpwind_r[15]+(0.31943828249997*p_over_gamma[5]+0.5*p_over_gamma[0])*fUpwind_r[13]+0.4472135954999579*p_over_gamma[6]*fUpwind_r[11]+0.5*fUpwind_r[4]*p_over_gamma[7]+0.447213595499958*p_over_gamma[3]*fUpwind_r[7]+0.5000000000000001*fUpwind_r[1]*p_over_gamma[5]+0.447213595499958*p_over_gamma[2]*fUpwind_r[5]; 
  Ghat_r[14] = (0.31943828249997*p_over_gamma[5]+0.4472135954999579*p_over_gamma[4]+0.5*p_over_gamma[0])*fUpwind_r[14]+(0.31943828249997*p_over_gamma[7]+0.5000000000000001*p_over_gamma[1])*fUpwind_r[12]+0.4*p_over_gamma[3]*fUpwind_r[10]+p_over_gamma[7]*(0.4472135954999579*fUpwind_r[8]+0.5*fUpwind_r[0])+fUpwind_r[6]*(0.4*p_over_gamma[6]+0.447213595499958*p_over_gamma[2])+0.5000000000000001*fUpwind_r[2]*p_over_gamma[5]+0.447213595499958*fUpwind_r[3]*p_over_gamma[3]; 
  Ghat_r[15] = (0.31943828249997*p_over_gamma[5]+0.4472135954999579*p_over_gamma[4]+0.5*p_over_gamma[0])*fUpwind_r[15]+(0.31943828249997*p_over_gamma[7]+0.5000000000000001*p_over_gamma[1])*fUpwind_r[13]+0.4*p_over_gamma[3]*fUpwind_r[11]+p_over_gamma[7]*(0.4472135954999579*fUpwind_r[9]+0.5000000000000001*fUpwind_r[1])+(0.4*p_over_gamma[6]+0.4472135954999579*p_over_gamma[2])*fUpwind_r[7]+0.5*fUpwind_r[4]*p_over_gamma[5]+0.4472135954999579*p_over_gamma[3]*fUpwind_r[5]; 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dx11; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dx11; 
  out[2] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dx11; 
  out[3] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dx11; 
  out[4] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dx11; 
  out[5] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dx11; 
  out[6] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*dx11; 
  out[7] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dx11; 
  out[8] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*dx11; 
  out[9] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dx11; 
  out[10] += (0.7071067811865475*Ghat_l[6]-0.7071067811865475*Ghat_r[6])*dx11; 
  out[11] += -1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dx11; 
  out[12] += -1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dx11; 
  out[13] += (0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*dx11; 
  out[14] += -1.224744871391589*(Ghat_r[6]+Ghat_l[6])*dx11; 
  out[15] += -1.224744871391589*(Ghat_r[7]+Ghat_l[7])*dx11; 
  out[16] += (0.7071067811865475*Ghat_l[8]-0.7071067811865475*Ghat_r[8])*dx11; 
  out[17] += (0.7071067811865475*Ghat_l[9]-0.7071067811865475*Ghat_r[9])*dx11; 
  out[18] += -1.224744871391589*(Ghat_r[8]+Ghat_l[8])*dx11; 
  out[19] += (0.7071067811865475*Ghat_l[10]-0.7071067811865475*Ghat_r[10])*dx11; 
  out[20] += -1.224744871391589*(Ghat_r[9]+Ghat_l[9])*dx11; 
  out[21] += (0.7071067811865475*Ghat_l[11]-0.7071067811865475*Ghat_r[11])*dx11; 
  out[22] += -1.224744871391589*(Ghat_r[10]+Ghat_l[10])*dx11; 
  out[23] += -1.224744871391589*(Ghat_r[11]+Ghat_l[11])*dx11; 
  out[24] += (0.7071067811865475*Ghat_l[12]-0.7071067811865475*Ghat_r[12])*dx11; 
  out[25] += (0.7071067811865475*Ghat_l[13]-0.7071067811865475*Ghat_r[13])*dx11; 
  out[26] += -1.224744871391589*(Ghat_r[12]+Ghat_l[12])*dx11; 
  out[27] += (0.7071067811865475*Ghat_l[14]-0.7071067811865475*Ghat_r[14])*dx11; 
  out[28] += -1.224744871391589*(Ghat_r[13]+Ghat_l[13])*dx11; 
  out[29] += (0.7071067811865475*Ghat_l[15]-0.7071067811865475*Ghat_r[15])*dx11; 
  out[30] += -1.224744871391589*(Ghat_r[14]+Ghat_l[14])*dx11; 
  out[31] += -1.224744871391589*(Ghat_r[15]+Ghat_l[15])*dx11; 

  return 0.;

} 
