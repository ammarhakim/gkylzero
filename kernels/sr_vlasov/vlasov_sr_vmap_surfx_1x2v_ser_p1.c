#include <gkyl_vlasov_sr_kernels.h> 
GKYL_CU_DH double vlasov_sr_vmap_surfx_1x2v_ser_p1(const double *w, const double *dxv, const double *jacob_vel_inv, const double *gamma, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:             Cell-center coordinates.
  // dxv[NDIM]:           Cell spacing.
  // jacob_vel_inv[VDIM]: Inverse velocity space Jacobian in each direction.
  // gamma:               Particle Lorentz boost factor sqrt(1 + p^2).
  // fl/fc/fr:            Input Distribution function in left/center/right cells.
  // out:                 Incremented distribution function in center cell.
  const double dx10 = 2.0/dxv[0]; 
  const double dv10 = 2.0/dxv[1]; 
  const double wv = w[1]; 
  const double *jacob_vel_inv_dir = &jacob_vel_inv[0]; 

  double p_over_gamma[8] = {0.0}; 
  p_over_gamma[0] = 2.738612787525831*jacob_vel_inv_dir[1]*gamma[4]*dv10+1.224744871391589*jacob_vel_inv_dir[0]*gamma[1]*dv10; 
  p_over_gamma[1] = 2.449489742783178*jacob_vel_inv_dir[2]*gamma[4]*dv10+2.738612787525831*jacob_vel_inv_dir[0]*gamma[4]*dv10+1.224744871391589*jacob_vel_inv_dir[1]*gamma[1]*dv10; 
  p_over_gamma[2] = 2.738612787525831*jacob_vel_inv_dir[1]*gamma[6]*dv10+1.224744871391589*jacob_vel_inv_dir[0]*gamma[3]*dv10; 
  p_over_gamma[3] = 2.449489742783178*jacob_vel_inv_dir[2]*gamma[6]*dv10+2.738612787525831*jacob_vel_inv_dir[0]*gamma[6]*dv10+1.224744871391589*jacob_vel_inv_dir[1]*gamma[3]*dv10; 
  p_over_gamma[4] = 2.449489742783178*jacob_vel_inv_dir[1]*gamma[4]*dv10+1.224744871391589*gamma[1]*jacob_vel_inv_dir[2]*dv10; 
  p_over_gamma[5] = 1.224744871391589*jacob_vel_inv_dir[0]*gamma[7]*dv10; 
  p_over_gamma[6] = 2.449489742783178*jacob_vel_inv_dir[1]*gamma[6]*dv10+1.224744871391589*jacob_vel_inv_dir[2]*gamma[3]*dv10; 
  p_over_gamma[7] = 1.224744871391589*jacob_vel_inv_dir[1]*gamma[7]*dv10; 
  double fUpwind_l[8]; 
  double fUpwind_r[8]; 
  if (wv>0) { 

  fUpwind_l[0] = 1.224744871391589*fl[1]+0.7071067811865475*fl[0]; 
  fUpwind_l[1] = 1.224744871391589*fl[4]+0.7071067811865475*fl[2]; 
  fUpwind_l[2] = 1.224744871391589*fl[5]+0.7071067811865475*fl[3]; 
  fUpwind_l[3] = 1.224744871391589*fl[7]+0.7071067811865475*fl[6]; 
  fUpwind_l[4] = 1.224744871391589*fl[9]+0.7071067811865475*fl[8]; 
  fUpwind_l[5] = 1.224744871391589*fl[13]+0.7071067811865475*fl[12]; 
  fUpwind_l[6] = 1.224744871391589*fl[11]+0.7071067811865475*fl[10]; 
  fUpwind_l[7] = 1.224744871391589*fl[15]+0.7071067811865475*fl[14]; 

  fUpwind_r[0] = 1.224744871391589*fc[1]+0.7071067811865475*fc[0]; 
  fUpwind_r[1] = 1.224744871391589*fc[4]+0.7071067811865475*fc[2]; 
  fUpwind_r[2] = 1.224744871391589*fc[5]+0.7071067811865475*fc[3]; 
  fUpwind_r[3] = 1.224744871391589*fc[7]+0.7071067811865475*fc[6]; 
  fUpwind_r[4] = 1.224744871391589*fc[9]+0.7071067811865475*fc[8]; 
  fUpwind_r[5] = 1.224744871391589*fc[13]+0.7071067811865475*fc[12]; 
  fUpwind_r[6] = 1.224744871391589*fc[11]+0.7071067811865475*fc[10]; 
  fUpwind_r[7] = 1.224744871391589*fc[15]+0.7071067811865475*fc[14]; 

  } else { 

  fUpwind_l[0] = 0.7071067811865475*fc[0]-1.224744871391589*fc[1]; 
  fUpwind_l[1] = 0.7071067811865475*fc[2]-1.224744871391589*fc[4]; 
  fUpwind_l[2] = 0.7071067811865475*fc[3]-1.224744871391589*fc[5]; 
  fUpwind_l[3] = 0.7071067811865475*fc[6]-1.224744871391589*fc[7]; 
  fUpwind_l[4] = 0.7071067811865475*fc[8]-1.224744871391589*fc[9]; 
  fUpwind_l[5] = 0.7071067811865475*fc[12]-1.224744871391589*fc[13]; 
  fUpwind_l[6] = 0.7071067811865475*fc[10]-1.224744871391589*fc[11]; 
  fUpwind_l[7] = 0.7071067811865475*fc[14]-1.224744871391589*fc[15]; 

  fUpwind_r[0] = 0.7071067811865475*fr[0]-1.224744871391589*fr[1]; 
  fUpwind_r[1] = 0.7071067811865475*fr[2]-1.224744871391589*fr[4]; 
  fUpwind_r[2] = 0.7071067811865475*fr[3]-1.224744871391589*fr[5]; 
  fUpwind_r[3] = 0.7071067811865475*fr[6]-1.224744871391589*fr[7]; 
  fUpwind_r[4] = 0.7071067811865475*fr[8]-1.224744871391589*fr[9]; 
  fUpwind_r[5] = 0.7071067811865475*fr[12]-1.224744871391589*fr[13]; 
  fUpwind_r[6] = 0.7071067811865475*fr[10]-1.224744871391589*fr[11]; 
  fUpwind_r[7] = 0.7071067811865475*fr[14]-1.224744871391589*fr[15]; 

  } 
  double Ghat_l[8]; 
  double Ghat_r[8]; 
  Ghat_l[0] = 0.5*(fUpwind_l[7]*p_over_gamma[7]+fUpwind_l[6]*p_over_gamma[6]+fUpwind_l[5]*p_over_gamma[5]+fUpwind_l[4]*p_over_gamma[4]+fUpwind_l[3]*p_over_gamma[3]+fUpwind_l[2]*p_over_gamma[2]+fUpwind_l[1]*p_over_gamma[1]+fUpwind_l[0]*p_over_gamma[0]); 
  Ghat_l[1] = 0.5000000000000001*(fUpwind_l[5]*p_over_gamma[7]+p_over_gamma[5]*fUpwind_l[7])+0.447213595499958*(fUpwind_l[3]*p_over_gamma[6]+p_over_gamma[3]*fUpwind_l[6])+0.4472135954999579*(fUpwind_l[1]*p_over_gamma[4]+p_over_gamma[1]*fUpwind_l[4])+0.5*(fUpwind_l[2]*p_over_gamma[3]+p_over_gamma[2]*fUpwind_l[3]+fUpwind_l[0]*p_over_gamma[1]+p_over_gamma[0]*fUpwind_l[1]); 
  Ghat_l[2] = 0.447213595499958*(fUpwind_l[3]*p_over_gamma[7]+p_over_gamma[3]*fUpwind_l[7])+0.5000000000000001*(fUpwind_l[4]*p_over_gamma[6]+p_over_gamma[4]*fUpwind_l[6])+0.4472135954999579*(fUpwind_l[2]*p_over_gamma[5]+p_over_gamma[2]*fUpwind_l[5])+0.5*(fUpwind_l[1]*p_over_gamma[3]+p_over_gamma[1]*fUpwind_l[3]+fUpwind_l[0]*p_over_gamma[2]+p_over_gamma[0]*fUpwind_l[2]); 
  Ghat_l[3] = (0.4*fUpwind_l[6]+0.447213595499958*fUpwind_l[2])*p_over_gamma[7]+0.4*p_over_gamma[6]*fUpwind_l[7]+0.447213595499958*(p_over_gamma[2]*fUpwind_l[7]+fUpwind_l[1]*p_over_gamma[6]+p_over_gamma[1]*fUpwind_l[6])+0.4472135954999579*(fUpwind_l[3]*p_over_gamma[5]+p_over_gamma[3]*fUpwind_l[5]+fUpwind_l[3]*p_over_gamma[4]+p_over_gamma[3]*fUpwind_l[4])+0.5*(fUpwind_l[0]*p_over_gamma[3]+p_over_gamma[0]*fUpwind_l[3]+fUpwind_l[1]*p_over_gamma[2]+p_over_gamma[1]*fUpwind_l[2]); 
  Ghat_l[4] = 0.4472135954999579*fUpwind_l[7]*p_over_gamma[7]+0.31943828249997*fUpwind_l[6]*p_over_gamma[6]+0.5000000000000001*(fUpwind_l[2]*p_over_gamma[6]+p_over_gamma[2]*fUpwind_l[6])+0.31943828249997*fUpwind_l[4]*p_over_gamma[4]+0.5*(fUpwind_l[0]*p_over_gamma[4]+p_over_gamma[0]*fUpwind_l[4])+0.4472135954999579*(fUpwind_l[3]*p_over_gamma[3]+fUpwind_l[1]*p_over_gamma[1]); 
  Ghat_l[5] = 0.31943828249997*fUpwind_l[7]*p_over_gamma[7]+0.5000000000000001*(fUpwind_l[1]*p_over_gamma[7]+p_over_gamma[1]*fUpwind_l[7])+0.4472135954999579*fUpwind_l[6]*p_over_gamma[6]+0.31943828249997*fUpwind_l[5]*p_over_gamma[5]+0.5*(fUpwind_l[0]*p_over_gamma[5]+p_over_gamma[0]*fUpwind_l[5])+0.4472135954999579*(fUpwind_l[3]*p_over_gamma[3]+fUpwind_l[2]*p_over_gamma[2]); 
  Ghat_l[6] = 0.4*(fUpwind_l[3]*p_over_gamma[7]+p_over_gamma[3]*fUpwind_l[7])+(0.4472135954999579*fUpwind_l[5]+0.31943828249997*fUpwind_l[4]+0.5*fUpwind_l[0])*p_over_gamma[6]+(0.4472135954999579*p_over_gamma[5]+0.31943828249997*p_over_gamma[4]+0.5*p_over_gamma[0])*fUpwind_l[6]+0.5000000000000001*(fUpwind_l[2]*p_over_gamma[4]+p_over_gamma[2]*fUpwind_l[4])+0.447213595499958*(fUpwind_l[1]*p_over_gamma[3]+p_over_gamma[1]*fUpwind_l[3]); 
  Ghat_l[7] = (0.31943828249997*fUpwind_l[5]+0.4472135954999579*fUpwind_l[4]+0.5*fUpwind_l[0])*p_over_gamma[7]+(0.31943828249997*p_over_gamma[5]+0.4472135954999579*p_over_gamma[4]+0.5*p_over_gamma[0])*fUpwind_l[7]+0.4*(fUpwind_l[3]*p_over_gamma[6]+p_over_gamma[3]*fUpwind_l[6])+0.5000000000000001*(fUpwind_l[1]*p_over_gamma[5]+p_over_gamma[1]*fUpwind_l[5])+0.447213595499958*(fUpwind_l[2]*p_over_gamma[3]+p_over_gamma[2]*fUpwind_l[3]); 

  Ghat_r[0] = 0.5*(fUpwind_r[7]*p_over_gamma[7]+fUpwind_r[6]*p_over_gamma[6]+fUpwind_r[5]*p_over_gamma[5]+fUpwind_r[4]*p_over_gamma[4]+fUpwind_r[3]*p_over_gamma[3]+fUpwind_r[2]*p_over_gamma[2]+fUpwind_r[1]*p_over_gamma[1]+fUpwind_r[0]*p_over_gamma[0]); 
  Ghat_r[1] = 0.5000000000000001*(fUpwind_r[5]*p_over_gamma[7]+p_over_gamma[5]*fUpwind_r[7])+0.447213595499958*(fUpwind_r[3]*p_over_gamma[6]+p_over_gamma[3]*fUpwind_r[6])+0.4472135954999579*(fUpwind_r[1]*p_over_gamma[4]+p_over_gamma[1]*fUpwind_r[4])+0.5*(fUpwind_r[2]*p_over_gamma[3]+p_over_gamma[2]*fUpwind_r[3]+fUpwind_r[0]*p_over_gamma[1]+p_over_gamma[0]*fUpwind_r[1]); 
  Ghat_r[2] = 0.447213595499958*(fUpwind_r[3]*p_over_gamma[7]+p_over_gamma[3]*fUpwind_r[7])+0.5000000000000001*(fUpwind_r[4]*p_over_gamma[6]+p_over_gamma[4]*fUpwind_r[6])+0.4472135954999579*(fUpwind_r[2]*p_over_gamma[5]+p_over_gamma[2]*fUpwind_r[5])+0.5*(fUpwind_r[1]*p_over_gamma[3]+p_over_gamma[1]*fUpwind_r[3]+fUpwind_r[0]*p_over_gamma[2]+p_over_gamma[0]*fUpwind_r[2]); 
  Ghat_r[3] = (0.4*fUpwind_r[6]+0.447213595499958*fUpwind_r[2])*p_over_gamma[7]+0.4*p_over_gamma[6]*fUpwind_r[7]+0.447213595499958*(p_over_gamma[2]*fUpwind_r[7]+fUpwind_r[1]*p_over_gamma[6]+p_over_gamma[1]*fUpwind_r[6])+0.4472135954999579*(fUpwind_r[3]*p_over_gamma[5]+p_over_gamma[3]*fUpwind_r[5]+fUpwind_r[3]*p_over_gamma[4]+p_over_gamma[3]*fUpwind_r[4])+0.5*(fUpwind_r[0]*p_over_gamma[3]+p_over_gamma[0]*fUpwind_r[3]+fUpwind_r[1]*p_over_gamma[2]+p_over_gamma[1]*fUpwind_r[2]); 
  Ghat_r[4] = 0.4472135954999579*fUpwind_r[7]*p_over_gamma[7]+0.31943828249997*fUpwind_r[6]*p_over_gamma[6]+0.5000000000000001*(fUpwind_r[2]*p_over_gamma[6]+p_over_gamma[2]*fUpwind_r[6])+0.31943828249997*fUpwind_r[4]*p_over_gamma[4]+0.5*(fUpwind_r[0]*p_over_gamma[4]+p_over_gamma[0]*fUpwind_r[4])+0.4472135954999579*(fUpwind_r[3]*p_over_gamma[3]+fUpwind_r[1]*p_over_gamma[1]); 
  Ghat_r[5] = 0.31943828249997*fUpwind_r[7]*p_over_gamma[7]+0.5000000000000001*(fUpwind_r[1]*p_over_gamma[7]+p_over_gamma[1]*fUpwind_r[7])+0.4472135954999579*fUpwind_r[6]*p_over_gamma[6]+0.31943828249997*fUpwind_r[5]*p_over_gamma[5]+0.5*(fUpwind_r[0]*p_over_gamma[5]+p_over_gamma[0]*fUpwind_r[5])+0.4472135954999579*(fUpwind_r[3]*p_over_gamma[3]+fUpwind_r[2]*p_over_gamma[2]); 
  Ghat_r[6] = 0.4*(fUpwind_r[3]*p_over_gamma[7]+p_over_gamma[3]*fUpwind_r[7])+(0.4472135954999579*fUpwind_r[5]+0.31943828249997*fUpwind_r[4]+0.5*fUpwind_r[0])*p_over_gamma[6]+(0.4472135954999579*p_over_gamma[5]+0.31943828249997*p_over_gamma[4]+0.5*p_over_gamma[0])*fUpwind_r[6]+0.5000000000000001*(fUpwind_r[2]*p_over_gamma[4]+p_over_gamma[2]*fUpwind_r[4])+0.447213595499958*(fUpwind_r[1]*p_over_gamma[3]+p_over_gamma[1]*fUpwind_r[3]); 
  Ghat_r[7] = (0.31943828249997*fUpwind_r[5]+0.4472135954999579*fUpwind_r[4]+0.5*fUpwind_r[0])*p_over_gamma[7]+(0.31943828249997*p_over_gamma[5]+0.4472135954999579*p_over_gamma[4]+0.5*p_over_gamma[0])*fUpwind_r[7]+0.4*(fUpwind_r[3]*p_over_gamma[6]+p_over_gamma[3]*fUpwind_r[6])+0.5000000000000001*(fUpwind_r[1]*p_over_gamma[5]+p_over_gamma[1]*fUpwind_r[5])+0.447213595499958*(fUpwind_r[2]*p_over_gamma[3]+p_over_gamma[2]*fUpwind_r[3]); 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dx10; 
  out[1] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dx10; 
  out[2] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dx10; 
  out[3] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dx10; 
  out[4] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dx10; 
  out[5] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dx10; 
  out[6] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dx10; 
  out[7] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dx10; 
  out[8] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*dx10; 
  out[9] += -1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dx10; 
  out[10] += (0.7071067811865475*Ghat_l[6]-0.7071067811865475*Ghat_r[6])*dx10; 
  out[11] += -1.224744871391589*(Ghat_r[6]+Ghat_l[6])*dx10; 
  out[12] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*dx10; 
  out[13] += -1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dx10; 
  out[14] += (0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*dx10; 
  out[15] += -1.224744871391589*(Ghat_r[7]+Ghat_l[7])*dx10; 

  return 0.;

} 
