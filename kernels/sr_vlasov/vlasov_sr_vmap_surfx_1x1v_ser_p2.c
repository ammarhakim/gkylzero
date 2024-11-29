#include <gkyl_vlasov_sr_kernels.h> 
GKYL_CU_DH double vlasov_sr_vmap_surfx_1x1v_ser_p2(const double *w, const double *dxv, const double *jacob_vel_inv, const double *gamma, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
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

  double p_over_gamma[3] = {0.0}; 
  p_over_gamma[0] = 2.738612787525831*jacob_vel_inv_dir[1]*gamma[2]*dv10+1.224744871391589*jacob_vel_inv_dir[0]*gamma[1]*dv10; 
  p_over_gamma[1] = 2.449489742783178*jacob_vel_inv_dir[2]*gamma[2]*dv10+2.738612787525831*jacob_vel_inv_dir[0]*gamma[2]*dv10+1.224744871391589*jacob_vel_inv_dir[1]*gamma[1]*dv10; 
  p_over_gamma[2] = 2.449489742783178*jacob_vel_inv_dir[1]*gamma[2]*dv10+1.224744871391589*gamma[1]*jacob_vel_inv_dir[2]*dv10; 
  double fUpwind_l[3]; 
  double fUpwind_r[3]; 
  if (wv>0) { 

  fUpwind_l[0] = 1.58113883008419*fl[4]+1.224744871391589*fl[1]+0.7071067811865475*fl[0]; 
  fUpwind_l[1] = 1.58113883008419*fl[6]+1.224744871391589*fl[3]+0.7071067811865475*fl[2]; 
  fUpwind_l[2] = 1.224744871391589*fl[7]+0.7071067811865475*fl[5]; 

  fUpwind_r[0] = 1.58113883008419*fc[4]+1.224744871391589*fc[1]+0.7071067811865475*fc[0]; 
  fUpwind_r[1] = 1.58113883008419*fc[6]+1.224744871391589*fc[3]+0.7071067811865475*fc[2]; 
  fUpwind_r[2] = 1.224744871391589*fc[7]+0.7071067811865475*fc[5]; 

  } else { 

  fUpwind_l[0] = 1.58113883008419*fc[4]-1.224744871391589*fc[1]+0.7071067811865475*fc[0]; 
  fUpwind_l[1] = 1.58113883008419*fc[6]-1.224744871391589*fc[3]+0.7071067811865475*fc[2]; 
  fUpwind_l[2] = 0.7071067811865475*fc[5]-1.224744871391589*fc[7]; 

  fUpwind_r[0] = 1.58113883008419*fr[4]-1.224744871391589*fr[1]+0.7071067811865475*fr[0]; 
  fUpwind_r[1] = 1.58113883008419*fr[6]-1.224744871391589*fr[3]+0.7071067811865475*fr[2]; 
  fUpwind_r[2] = 0.7071067811865475*fr[5]-1.224744871391589*fr[7]; 

  } 
  double Ghat_l[3]; 
  double Ghat_r[3]; 
  Ghat_l[0] = 0.7071067811865475*(fUpwind_l[2]*p_over_gamma[2]+fUpwind_l[1]*p_over_gamma[1]+fUpwind_l[0]*p_over_gamma[0]); 
  Ghat_l[1] = 0.6324555320336759*(fUpwind_l[1]*p_over_gamma[2]+p_over_gamma[1]*fUpwind_l[2])+0.7071067811865475*(fUpwind_l[0]*p_over_gamma[1]+p_over_gamma[0]*fUpwind_l[1]); 
  Ghat_l[2] = 0.4517539514526256*fUpwind_l[2]*p_over_gamma[2]+0.7071067811865475*(fUpwind_l[0]*p_over_gamma[2]+p_over_gamma[0]*fUpwind_l[2])+0.6324555320336759*fUpwind_l[1]*p_over_gamma[1]; 

  Ghat_r[0] = 0.7071067811865475*(fUpwind_r[2]*p_over_gamma[2]+fUpwind_r[1]*p_over_gamma[1]+fUpwind_r[0]*p_over_gamma[0]); 
  Ghat_r[1] = 0.6324555320336759*(fUpwind_r[1]*p_over_gamma[2]+p_over_gamma[1]*fUpwind_r[2])+0.7071067811865475*(fUpwind_r[0]*p_over_gamma[1]+p_over_gamma[0]*fUpwind_r[1]); 
  Ghat_r[2] = 0.4517539514526256*fUpwind_r[2]*p_over_gamma[2]+0.7071067811865475*(fUpwind_r[0]*p_over_gamma[2]+p_over_gamma[0]*fUpwind_r[2])+0.6324555320336759*fUpwind_r[1]*p_over_gamma[1]; 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dx10; 
  out[1] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dx10; 
  out[2] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dx10; 
  out[3] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dx10; 
  out[4] += (1.58113883008419*Ghat_l[0]-1.58113883008419*Ghat_r[0])*dx10; 
  out[5] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dx10; 
  out[6] += (1.58113883008419*Ghat_l[1]-1.58113883008419*Ghat_r[1])*dx10; 
  out[7] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dx10; 

  return 0.;

} 
