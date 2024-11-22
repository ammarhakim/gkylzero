#include <gkyl_boltzmann_photon_kernels.h> 
GKYL_CU_DH double boltzmann_photon_boundary_surfvy_1x2v_ser_p1(const double *w, const double *dxv,  double light_speed, double rho_curv, 
  const double *jacob_vel_inv, const double *kpar_abs, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:             Cell-center coordinates.
  // dxv[NDIM]:           Cell spacing.
  // light_speed:         Speed of light.
  // rho_curv:            Curvature of the magnetic field.
  // jacob_vel_inv[VDIM]: Inverse velocity space Jacobian in each direction.
  // kpar_abs:            Continuous expansion of |kpar| on the grid.
  // edge:                Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge:         Input photon distribution function in skin cell/last edge cell 
  // out:                 Output photon distribution function in skin cell 
  const double dv10 = 2.0/dxv[1]; 
  const double dv11 = 2.0/dxv[2]; 
  const double *jacob_vel_inv_dir = &jacob_vel_inv[3]; 
  double alpha[6] = {0.0}; 
  double Ghat[6] = {0.0}; 
  if (edge == -1) { 

  alpha[0] = 2.23606797749979*kpar_abs[0]*jacob_vel_inv_dir[2]+1.732050807568877*kpar_abs[0]*jacob_vel_inv_dir[1]+jacob_vel_inv_dir[0]*kpar_abs[0]; 
  alpha[2] = 2.23606797749979*kpar_abs[1]*jacob_vel_inv_dir[2]+1.732050807568877*jacob_vel_inv_dir[1]*kpar_abs[1]+jacob_vel_inv_dir[0]*kpar_abs[1]; 
  alpha[4] = 2.23606797749979*jacob_vel_inv_dir[2]*kpar_abs[2]+1.732050807568877*jacob_vel_inv_dir[1]*kpar_abs[2]+jacob_vel_inv_dir[0]*kpar_abs[2]; 

  Ghat[0] = 0.7905694150420948*alpha[2]*fSkin[14]+0.7905694150420947*alpha[0]*fSkin[12]+0.6123724356957944*alpha[4]*fSkin[10]+0.3535533905932737*alpha[4]*fSkin[8]+0.6123724356957944*alpha[2]*fSkin[6]+0.6123724356957944*alpha[0]*fSkin[3]+0.3535533905932737*alpha[2]*fSkin[2]+0.3535533905932737*alpha[0]*fSkin[0]; 
  Ghat[1] = 0.7905694150420947*alpha[2]*fSkin[15]+0.7905694150420948*alpha[0]*fSkin[13]+0.6123724356957944*alpha[4]*fSkin[11]+0.3535533905932737*alpha[4]*fSkin[9]+0.6123724356957944*alpha[2]*fSkin[7]+0.6123724356957944*alpha[0]*fSkin[5]+0.3535533905932737*alpha[2]*fSkin[4]+0.3535533905932737*alpha[0]*fSkin[1]; 
  Ghat[2] = 0.7071067811865475*alpha[4]*fSkin[14]+0.7905694150420948*alpha[0]*fSkin[14]+0.7905694150420947*alpha[2]*fSkin[12]+0.5477225575051661*alpha[2]*fSkin[10]+0.3162277660168379*alpha[2]*fSkin[8]+0.5477225575051661*alpha[4]*fSkin[6]+0.6123724356957944*alpha[0]*fSkin[6]+0.3162277660168379*fSkin[2]*alpha[4]+0.6123724356957944*alpha[2]*fSkin[3]+0.3535533905932737*alpha[0]*fSkin[2]+0.3535533905932737*fSkin[0]*alpha[2]; 
  Ghat[3] = 0.7071067811865475*alpha[4]*fSkin[15]+0.7905694150420947*alpha[0]*fSkin[15]+0.7905694150420948*alpha[2]*fSkin[13]+0.5477225575051661*alpha[2]*fSkin[11]+0.3162277660168379*alpha[2]*fSkin[9]+0.5477225575051661*alpha[4]*fSkin[7]+0.6123724356957944*alpha[0]*fSkin[7]+0.6123724356957944*alpha[2]*fSkin[5]+0.3162277660168379*alpha[4]*fSkin[4]+0.3535533905932737*alpha[0]*fSkin[4]+0.3535533905932737*fSkin[1]*alpha[2]; 
  Ghat[4] = 0.7071067811865475*alpha[2]*fSkin[14]+0.7905694150420947*alpha[4]*fSkin[12]+0.3912303982179757*alpha[4]*fSkin[10]+0.6123724356957944*alpha[0]*fSkin[10]+0.2258769757263128*alpha[4]*fSkin[8]+0.3535533905932737*alpha[0]*fSkin[8]+0.5477225575051661*alpha[2]*fSkin[6]+0.6123724356957944*fSkin[3]*alpha[4]+0.3535533905932737*fSkin[0]*alpha[4]+0.3162277660168379*alpha[2]*fSkin[2]; 
  Ghat[5] = 0.7071067811865475*alpha[2]*fSkin[15]+0.7905694150420947*alpha[4]*fSkin[13]+0.3912303982179757*alpha[4]*fSkin[11]+0.6123724356957944*alpha[0]*fSkin[11]+0.2258769757263128*alpha[4]*fSkin[9]+0.3535533905932737*alpha[0]*fSkin[9]+0.5477225575051661*alpha[2]*fSkin[7]+0.6123724356957944*alpha[4]*fSkin[5]+0.3162277660168379*alpha[2]*fSkin[4]+0.3535533905932737*fSkin[1]*alpha[4]; 

  out[0] += -(0.7071067811865475*Ghat[0]*dv11*light_speed)/rho_curv; 
  out[1] += -(0.7071067811865475*Ghat[1]*dv11*light_speed)/rho_curv; 
  out[2] += -(0.7071067811865475*Ghat[2]*dv11*light_speed)/rho_curv; 
  out[3] += -(1.224744871391589*Ghat[0]*dv11*light_speed)/rho_curv; 
  out[4] += -(0.7071067811865475*Ghat[3]*dv11*light_speed)/rho_curv; 
  out[5] += -(1.224744871391589*Ghat[1]*dv11*light_speed)/rho_curv; 
  out[6] += -(1.224744871391589*Ghat[2]*dv11*light_speed)/rho_curv; 
  out[7] += -(1.224744871391589*Ghat[3]*dv11*light_speed)/rho_curv; 
  out[8] += -(0.7071067811865475*Ghat[4]*dv11*light_speed)/rho_curv; 
  out[9] += -(0.7071067811865475*Ghat[5]*dv11*light_speed)/rho_curv; 
  out[10] += -(1.224744871391589*Ghat[4]*dv11*light_speed)/rho_curv; 
  out[11] += -(1.224744871391589*Ghat[5]*dv11*light_speed)/rho_curv; 
  out[12] += -(1.58113883008419*Ghat[0]*dv11*light_speed)/rho_curv; 
  out[13] += -(1.58113883008419*Ghat[1]*dv11*light_speed)/rho_curv; 
  out[14] += -(1.58113883008419*Ghat[2]*dv11*light_speed)/rho_curv; 
  out[15] += -(1.58113883008419*Ghat[3]*dv11*light_speed)/rho_curv; 

    } 
    else { 
  alpha[0] = 2.23606797749979*kpar_abs[0]*jacob_vel_inv_dir[2]-1.732050807568877*kpar_abs[0]*jacob_vel_inv_dir[1]+jacob_vel_inv_dir[0]*kpar_abs[0]; 
  alpha[2] = 2.23606797749979*kpar_abs[1]*jacob_vel_inv_dir[2]-1.732050807568877*jacob_vel_inv_dir[1]*kpar_abs[1]+jacob_vel_inv_dir[0]*kpar_abs[1]; 
  alpha[4] = 2.23606797749979*jacob_vel_inv_dir[2]*kpar_abs[2]-1.732050807568877*jacob_vel_inv_dir[1]*kpar_abs[2]+jacob_vel_inv_dir[0]*kpar_abs[2]; 

  Ghat[0] = 0.7905694150420948*alpha[2]*fEdge[14]+0.7905694150420947*alpha[0]*fEdge[12]+0.6123724356957944*alpha[4]*fEdge[10]+0.3535533905932737*alpha[4]*fEdge[8]+0.6123724356957944*alpha[2]*fEdge[6]+0.6123724356957944*alpha[0]*fEdge[3]+0.3535533905932737*alpha[2]*fEdge[2]+0.3535533905932737*alpha[0]*fEdge[0]; 
  Ghat[1] = 0.7905694150420947*alpha[2]*fEdge[15]+0.7905694150420948*alpha[0]*fEdge[13]+0.6123724356957944*alpha[4]*fEdge[11]+0.3535533905932737*alpha[4]*fEdge[9]+0.6123724356957944*alpha[2]*fEdge[7]+0.6123724356957944*alpha[0]*fEdge[5]+0.3535533905932737*alpha[2]*fEdge[4]+0.3535533905932737*alpha[0]*fEdge[1]; 
  Ghat[2] = 0.7071067811865475*alpha[4]*fEdge[14]+0.7905694150420948*alpha[0]*fEdge[14]+0.7905694150420947*alpha[2]*fEdge[12]+0.5477225575051661*alpha[2]*fEdge[10]+0.3162277660168379*alpha[2]*fEdge[8]+0.5477225575051661*alpha[4]*fEdge[6]+0.6123724356957944*alpha[0]*fEdge[6]+0.3162277660168379*fEdge[2]*alpha[4]+0.6123724356957944*alpha[2]*fEdge[3]+0.3535533905932737*alpha[0]*fEdge[2]+0.3535533905932737*fEdge[0]*alpha[2]; 
  Ghat[3] = 0.7071067811865475*alpha[4]*fEdge[15]+0.7905694150420947*alpha[0]*fEdge[15]+0.7905694150420948*alpha[2]*fEdge[13]+0.5477225575051661*alpha[2]*fEdge[11]+0.3162277660168379*alpha[2]*fEdge[9]+0.5477225575051661*alpha[4]*fEdge[7]+0.6123724356957944*alpha[0]*fEdge[7]+0.6123724356957944*alpha[2]*fEdge[5]+0.3162277660168379*alpha[4]*fEdge[4]+0.3535533905932737*alpha[0]*fEdge[4]+0.3535533905932737*fEdge[1]*alpha[2]; 
  Ghat[4] = 0.7071067811865475*alpha[2]*fEdge[14]+0.7905694150420947*alpha[4]*fEdge[12]+0.3912303982179757*alpha[4]*fEdge[10]+0.6123724356957944*alpha[0]*fEdge[10]+0.2258769757263128*alpha[4]*fEdge[8]+0.3535533905932737*alpha[0]*fEdge[8]+0.5477225575051661*alpha[2]*fEdge[6]+0.6123724356957944*fEdge[3]*alpha[4]+0.3535533905932737*fEdge[0]*alpha[4]+0.3162277660168379*alpha[2]*fEdge[2]; 
  Ghat[5] = 0.7071067811865475*alpha[2]*fEdge[15]+0.7905694150420947*alpha[4]*fEdge[13]+0.3912303982179757*alpha[4]*fEdge[11]+0.6123724356957944*alpha[0]*fEdge[11]+0.2258769757263128*alpha[4]*fEdge[9]+0.3535533905932737*alpha[0]*fEdge[9]+0.5477225575051661*alpha[2]*fEdge[7]+0.6123724356957944*alpha[4]*fEdge[5]+0.3162277660168379*alpha[2]*fEdge[4]+0.3535533905932737*fEdge[1]*alpha[4]; 

  out[0] += (0.7071067811865475*Ghat[0]*dv11*light_speed)/rho_curv; 
  out[1] += (0.7071067811865475*Ghat[1]*dv11*light_speed)/rho_curv; 
  out[2] += (0.7071067811865475*Ghat[2]*dv11*light_speed)/rho_curv; 
  out[3] += -(1.224744871391589*Ghat[0]*dv11*light_speed)/rho_curv; 
  out[4] += (0.7071067811865475*Ghat[3]*dv11*light_speed)/rho_curv; 
  out[5] += -(1.224744871391589*Ghat[1]*dv11*light_speed)/rho_curv; 
  out[6] += -(1.224744871391589*Ghat[2]*dv11*light_speed)/rho_curv; 
  out[7] += -(1.224744871391589*Ghat[3]*dv11*light_speed)/rho_curv; 
  out[8] += (0.7071067811865475*Ghat[4]*dv11*light_speed)/rho_curv; 
  out[9] += (0.7071067811865475*Ghat[5]*dv11*light_speed)/rho_curv; 
  out[10] += -(1.224744871391589*Ghat[4]*dv11*light_speed)/rho_curv; 
  out[11] += -(1.224744871391589*Ghat[5]*dv11*light_speed)/rho_curv; 
  out[12] += (1.58113883008419*Ghat[0]*dv11*light_speed)/rho_curv; 
  out[13] += (1.58113883008419*Ghat[1]*dv11*light_speed)/rho_curv; 
  out[14] += (1.58113883008419*Ghat[2]*dv11*light_speed)/rho_curv; 
  out[15] += (1.58113883008419*Ghat[3]*dv11*light_speed)/rho_curv; 

  } 
  double cflFreq = light_speed/rho_curv*alpha[0]; 
  return 0.75*dv11*cflFreq; 

} 
