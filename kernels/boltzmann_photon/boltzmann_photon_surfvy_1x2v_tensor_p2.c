#include <gkyl_boltzmann_photon_kernels.h> 
GKYL_CU_DH double boltzmann_photon_surfvy_1x2v_tensor_p2(const double *w, const double *dxv,  double light_speed, double rho_curv, 
  const double *jacob_vel_inv, const double *kpar_abs, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:             Cell-center coordinates.
  // dxv[NDIM]:           Cell spacing.
  // light_speed:         Speed of light.
  // rho_curv:            Curvature of the magnetic field.
  // jacob_vel_inv[VDIM]: Inverse velocity space Jacobian in each direction.
  // kpar_abs:            Continuous expansion of |kpar| on the grid.
  // fl/fc/fr:            Input photon distribution function in left/center/right cells.
  // out:                 Incremented photon distribution function in center cell.
  const double dv10 = 2.0/dxv[1]; 
  const double dv11 = 2.0/dxv[2]; 
  const double *jacob_vel_inv_dir = &jacob_vel_inv[3]; 
  double alpha_r[9] = {0.0}; 
  double alpha_l[9] = {0.0}; 
  alpha_r[0] = 2.23606797749979*kpar_abs[0]*jacob_vel_inv_dir[2]+1.732050807568877*kpar_abs[0]*jacob_vel_inv_dir[1]+jacob_vel_inv_dir[0]*kpar_abs[0]; 
  alpha_r[2] = 2.23606797749979*kpar_abs[1]*jacob_vel_inv_dir[2]+1.732050807568877*jacob_vel_inv_dir[1]*kpar_abs[1]+jacob_vel_inv_dir[0]*kpar_abs[1]; 
  alpha_r[5] = 2.23606797749979*jacob_vel_inv_dir[2]*kpar_abs[2]+1.732050807568877*jacob_vel_inv_dir[1]*kpar_abs[2]+jacob_vel_inv_dir[0]*kpar_abs[2]; 

  alpha_l[0] = 2.23606797749979*kpar_abs[0]*jacob_vel_inv_dir[2]-1.732050807568877*kpar_abs[0]*jacob_vel_inv_dir[1]+jacob_vel_inv_dir[0]*kpar_abs[0]; 
  alpha_l[2] = 2.23606797749979*kpar_abs[1]*jacob_vel_inv_dir[2]-1.732050807568877*jacob_vel_inv_dir[1]*kpar_abs[1]+jacob_vel_inv_dir[0]*kpar_abs[1]; 
  alpha_l[5] = 2.23606797749979*jacob_vel_inv_dir[2]*kpar_abs[2]-1.732050807568877*jacob_vel_inv_dir[1]*kpar_abs[2]+jacob_vel_inv_dir[0]*kpar_abs[2]; 

  double Ghat_r[9] = {0.0}; 
  double Ghat_l[9] = {0.0}; 
  Ghat_l[0] = 0.7905694150420947*alpha_l[5]*fl[22]+0.7905694150420948*alpha_l[2]*fl[16]+0.6123724356957944*alpha_l[5]*fl[14]+0.7905694150420947*alpha_l[0]*fl[9]+0.3535533905932737*alpha_l[5]*fl[8]+0.6123724356957944*alpha_l[2]*fl[6]+0.6123724356957944*alpha_l[0]*fl[3]+0.3535533905932737*alpha_l[2]*fl[2]+0.3535533905932737*alpha_l[0]*fl[0]; 
  Ghat_l[1] = 0.7905694150420947*alpha_l[5]*fl[25]+0.7905694150420947*alpha_l[2]*fl[19]+0.6123724356957944*alpha_l[5]*fl[18]+0.7905694150420948*alpha_l[0]*fl[15]+0.3535533905932737*alpha_l[5]*fl[12]+0.6123724356957944*alpha_l[2]*fl[10]+0.6123724356957944*alpha_l[0]*fl[5]+0.3535533905932737*alpha_l[2]*fl[4]+0.3535533905932737*alpha_l[0]*fl[1]; 
  Ghat_l[2] = 0.7071067811865475*alpha_l[2]*fl[22]+0.7071067811865475*alpha_l[5]*fl[16]+0.7905694150420948*alpha_l[0]*fl[16]+0.5477225575051661*alpha_l[2]*fl[14]+0.7905694150420947*alpha_l[2]*fl[9]+0.3162277660168379*alpha_l[2]*fl[8]+0.5477225575051661*alpha_l[5]*fl[6]+0.6123724356957944*alpha_l[0]*fl[6]+0.3162277660168379*fl[2]*alpha_l[5]+0.6123724356957944*alpha_l[2]*fl[3]+0.3535533905932737*alpha_l[0]*fl[2]+0.3535533905932737*fl[0]*alpha_l[2]; 
  Ghat_l[3] = 0.7071067811865475*alpha_l[2]*fl[25]+0.7071067811865475*alpha_l[5]*fl[19]+0.7905694150420947*alpha_l[0]*fl[19]+0.5477225575051661*alpha_l[2]*fl[18]+0.7905694150420948*alpha_l[2]*fl[15]+0.3162277660168379*alpha_l[2]*fl[12]+0.5477225575051661*alpha_l[5]*fl[10]+0.6123724356957944*alpha_l[0]*fl[10]+0.6123724356957944*alpha_l[2]*fl[5]+0.3162277660168379*fl[4]*alpha_l[5]+0.3535533905932737*alpha_l[0]*fl[4]+0.3535533905932737*fl[1]*alpha_l[2]; 
  Ghat_l[4] = 0.7905694150420947*alpha_l[5]*fl[26]+0.7905694150420947*alpha_l[2]*fl[24]+0.6123724356957944*alpha_l[5]*fl[23]+0.7905694150420947*alpha_l[0]*fl[21]+0.3535533905932737*alpha_l[5]*fl[20]+0.6123724356957944*alpha_l[2]*fl[17]+0.6123724356957944*alpha_l[0]*fl[13]+0.3535533905932737*alpha_l[2]*fl[11]+0.3535533905932737*alpha_l[0]*fl[7]; 
  Ghat_l[5] = 0.5050762722761053*alpha_l[5]*fl[22]+0.7905694150420947*alpha_l[0]*fl[22]+0.7071067811865475*alpha_l[2]*fl[16]+0.3912303982179757*alpha_l[5]*fl[14]+0.6123724356957944*alpha_l[0]*fl[14]+0.7905694150420947*alpha_l[5]*fl[9]+0.2258769757263128*alpha_l[5]*fl[8]+0.3535533905932737*alpha_l[0]*fl[8]+0.5477225575051661*alpha_l[2]*fl[6]+0.6123724356957944*fl[3]*alpha_l[5]+0.3535533905932737*fl[0]*alpha_l[5]+0.3162277660168379*alpha_l[2]*fl[2]; 
  Ghat_l[6] = 0.7071067811865475*alpha_l[2]*fl[26]+0.7071067811865475*alpha_l[5]*fl[24]+0.7905694150420948*alpha_l[0]*fl[24]+0.5477225575051661*alpha_l[2]*fl[23]+0.7905694150420948*alpha_l[2]*fl[21]+0.3162277660168379*alpha_l[2]*fl[20]+0.5477225575051661*alpha_l[5]*fl[17]+0.6123724356957944*alpha_l[0]*fl[17]+0.6123724356957944*alpha_l[2]*fl[13]+0.3162277660168379*alpha_l[5]*fl[11]+0.3535533905932737*alpha_l[0]*fl[11]+0.3535533905932737*alpha_l[2]*fl[7]; 
  Ghat_l[7] = 0.5050762722761054*alpha_l[5]*fl[25]+0.7905694150420948*alpha_l[0]*fl[25]+0.7071067811865475*alpha_l[2]*fl[19]+0.3912303982179757*alpha_l[5]*fl[18]+0.6123724356957944*alpha_l[0]*fl[18]+0.7905694150420947*alpha_l[5]*fl[15]+0.2258769757263128*alpha_l[5]*fl[12]+0.3535533905932737*alpha_l[0]*fl[12]+0.5477225575051661*alpha_l[2]*fl[10]+0.6123724356957944*alpha_l[5]*fl[5]+0.3535533905932737*fl[1]*alpha_l[5]+0.3162277660168379*alpha_l[2]*fl[4]; 
  Ghat_l[8] = 0.5050762722761053*alpha_l[5]*fl[26]+0.7905694150420947*alpha_l[0]*fl[26]+0.7071067811865475*alpha_l[2]*fl[24]+0.3912303982179757*alpha_l[5]*fl[23]+0.6123724356957944*alpha_l[0]*fl[23]+0.7905694150420947*alpha_l[5]*fl[21]+0.2258769757263128*alpha_l[5]*fl[20]+0.3535533905932737*alpha_l[0]*fl[20]+0.5477225575051661*alpha_l[2]*fl[17]+0.6123724356957944*alpha_l[5]*fl[13]+0.3162277660168379*alpha_l[2]*fl[11]+0.3535533905932737*alpha_l[5]*fl[7]; 

  Ghat_r[0] = 0.7905694150420947*alpha_r[5]*fc[22]+0.7905694150420948*alpha_r[2]*fc[16]+0.6123724356957944*alpha_r[5]*fc[14]+0.7905694150420947*alpha_r[0]*fc[9]+0.3535533905932737*alpha_r[5]*fc[8]+0.6123724356957944*alpha_r[2]*fc[6]+0.6123724356957944*alpha_r[0]*fc[3]+0.3535533905932737*alpha_r[2]*fc[2]+0.3535533905932737*alpha_r[0]*fc[0]; 
  Ghat_r[1] = 0.7905694150420947*alpha_r[5]*fc[25]+0.7905694150420947*alpha_r[2]*fc[19]+0.6123724356957944*alpha_r[5]*fc[18]+0.7905694150420948*alpha_r[0]*fc[15]+0.3535533905932737*alpha_r[5]*fc[12]+0.6123724356957944*alpha_r[2]*fc[10]+0.6123724356957944*alpha_r[0]*fc[5]+0.3535533905932737*alpha_r[2]*fc[4]+0.3535533905932737*alpha_r[0]*fc[1]; 
  Ghat_r[2] = 0.7071067811865475*alpha_r[2]*fc[22]+0.7071067811865475*alpha_r[5]*fc[16]+0.7905694150420948*alpha_r[0]*fc[16]+0.5477225575051661*alpha_r[2]*fc[14]+0.7905694150420947*alpha_r[2]*fc[9]+0.3162277660168379*alpha_r[2]*fc[8]+0.5477225575051661*alpha_r[5]*fc[6]+0.6123724356957944*alpha_r[0]*fc[6]+0.3162277660168379*fc[2]*alpha_r[5]+0.6123724356957944*alpha_r[2]*fc[3]+0.3535533905932737*alpha_r[0]*fc[2]+0.3535533905932737*fc[0]*alpha_r[2]; 
  Ghat_r[3] = 0.7071067811865475*alpha_r[2]*fc[25]+0.7071067811865475*alpha_r[5]*fc[19]+0.7905694150420947*alpha_r[0]*fc[19]+0.5477225575051661*alpha_r[2]*fc[18]+0.7905694150420948*alpha_r[2]*fc[15]+0.3162277660168379*alpha_r[2]*fc[12]+0.5477225575051661*alpha_r[5]*fc[10]+0.6123724356957944*alpha_r[0]*fc[10]+0.6123724356957944*alpha_r[2]*fc[5]+0.3162277660168379*fc[4]*alpha_r[5]+0.3535533905932737*alpha_r[0]*fc[4]+0.3535533905932737*fc[1]*alpha_r[2]; 
  Ghat_r[4] = 0.7905694150420947*alpha_r[5]*fc[26]+0.7905694150420947*alpha_r[2]*fc[24]+0.6123724356957944*alpha_r[5]*fc[23]+0.7905694150420947*alpha_r[0]*fc[21]+0.3535533905932737*alpha_r[5]*fc[20]+0.6123724356957944*alpha_r[2]*fc[17]+0.6123724356957944*alpha_r[0]*fc[13]+0.3535533905932737*alpha_r[2]*fc[11]+0.3535533905932737*alpha_r[0]*fc[7]; 
  Ghat_r[5] = 0.5050762722761053*alpha_r[5]*fc[22]+0.7905694150420947*alpha_r[0]*fc[22]+0.7071067811865475*alpha_r[2]*fc[16]+0.3912303982179757*alpha_r[5]*fc[14]+0.6123724356957944*alpha_r[0]*fc[14]+0.7905694150420947*alpha_r[5]*fc[9]+0.2258769757263128*alpha_r[5]*fc[8]+0.3535533905932737*alpha_r[0]*fc[8]+0.5477225575051661*alpha_r[2]*fc[6]+0.6123724356957944*fc[3]*alpha_r[5]+0.3535533905932737*fc[0]*alpha_r[5]+0.3162277660168379*alpha_r[2]*fc[2]; 
  Ghat_r[6] = 0.7071067811865475*alpha_r[2]*fc[26]+0.7071067811865475*alpha_r[5]*fc[24]+0.7905694150420948*alpha_r[0]*fc[24]+0.5477225575051661*alpha_r[2]*fc[23]+0.7905694150420948*alpha_r[2]*fc[21]+0.3162277660168379*alpha_r[2]*fc[20]+0.5477225575051661*alpha_r[5]*fc[17]+0.6123724356957944*alpha_r[0]*fc[17]+0.6123724356957944*alpha_r[2]*fc[13]+0.3162277660168379*alpha_r[5]*fc[11]+0.3535533905932737*alpha_r[0]*fc[11]+0.3535533905932737*alpha_r[2]*fc[7]; 
  Ghat_r[7] = 0.5050762722761054*alpha_r[5]*fc[25]+0.7905694150420948*alpha_r[0]*fc[25]+0.7071067811865475*alpha_r[2]*fc[19]+0.3912303982179757*alpha_r[5]*fc[18]+0.6123724356957944*alpha_r[0]*fc[18]+0.7905694150420947*alpha_r[5]*fc[15]+0.2258769757263128*alpha_r[5]*fc[12]+0.3535533905932737*alpha_r[0]*fc[12]+0.5477225575051661*alpha_r[2]*fc[10]+0.6123724356957944*alpha_r[5]*fc[5]+0.3535533905932737*fc[1]*alpha_r[5]+0.3162277660168379*alpha_r[2]*fc[4]; 
  Ghat_r[8] = 0.5050762722761053*alpha_r[5]*fc[26]+0.7905694150420947*alpha_r[0]*fc[26]+0.7071067811865475*alpha_r[2]*fc[24]+0.3912303982179757*alpha_r[5]*fc[23]+0.6123724356957944*alpha_r[0]*fc[23]+0.7905694150420947*alpha_r[5]*fc[21]+0.2258769757263128*alpha_r[5]*fc[20]+0.3535533905932737*alpha_r[0]*fc[20]+0.5477225575051661*alpha_r[2]*fc[17]+0.6123724356957944*alpha_r[5]*fc[13]+0.3162277660168379*alpha_r[2]*fc[11]+0.3535533905932737*alpha_r[5]*fc[7]; 

  out[0] += ((0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv11*light_speed)/rho_curv; 
  out[1] += ((0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv11*light_speed)/rho_curv; 
  out[2] += ((0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dv11*light_speed)/rho_curv; 
  out[3] += -(1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv11*light_speed)/rho_curv; 
  out[4] += ((0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dv11*light_speed)/rho_curv; 
  out[5] += -(1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv11*light_speed)/rho_curv; 
  out[6] += -(1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dv11*light_speed)/rho_curv; 
  out[7] += ((0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*dv11*light_speed)/rho_curv; 
  out[8] += ((0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*dv11*light_speed)/rho_curv; 
  out[9] += ((1.58113883008419*Ghat_l[0]-1.58113883008419*Ghat_r[0])*dv11*light_speed)/rho_curv; 
  out[10] += -(1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dv11*light_speed)/rho_curv; 
  out[11] += ((0.7071067811865475*Ghat_l[6]-0.7071067811865475*Ghat_r[6])*dv11*light_speed)/rho_curv; 
  out[12] += ((0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*dv11*light_speed)/rho_curv; 
  out[13] += -(1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dv11*light_speed)/rho_curv; 
  out[14] += -(1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dv11*light_speed)/rho_curv; 
  out[15] += ((1.58113883008419*Ghat_l[1]-1.58113883008419*Ghat_r[1])*dv11*light_speed)/rho_curv; 
  out[16] += ((1.58113883008419*Ghat_l[2]-1.58113883008419*Ghat_r[2])*dv11*light_speed)/rho_curv; 
  out[17] += -(1.224744871391589*(Ghat_r[6]+Ghat_l[6])*dv11*light_speed)/rho_curv; 
  out[18] += -(1.224744871391589*(Ghat_r[7]+Ghat_l[7])*dv11*light_speed)/rho_curv; 
  out[19] += ((1.58113883008419*Ghat_l[3]-1.58113883008419*Ghat_r[3])*dv11*light_speed)/rho_curv; 
  out[20] += ((0.7071067811865475*Ghat_l[8]-0.7071067811865475*Ghat_r[8])*dv11*light_speed)/rho_curv; 
  out[21] += ((1.58113883008419*Ghat_l[4]-1.58113883008419*Ghat_r[4])*dv11*light_speed)/rho_curv; 
  out[22] += ((1.58113883008419*Ghat_l[5]-1.58113883008419*Ghat_r[5])*dv11*light_speed)/rho_curv; 
  out[23] += -(1.224744871391589*(Ghat_r[8]+Ghat_l[8])*dv11*light_speed)/rho_curv; 
  out[24] += ((1.58113883008419*Ghat_l[6]-1.58113883008419*Ghat_r[6])*dv11*light_speed)/rho_curv; 
  out[25] += ((1.58113883008419*Ghat_l[7]-1.58113883008419*Ghat_r[7])*dv11*light_speed)/rho_curv; 
  out[26] += ((1.58113883008419*Ghat_l[8]-1.58113883008419*Ghat_r[8])*dv11*light_speed)/rho_curv; 

  double cflFreq = light_speed/rho_curv*fmax(fabs(alpha_l[0]), fabs(alpha_r[0])); 
  return 1.25*dv11*cflFreq; 

} 
