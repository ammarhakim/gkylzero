#include <gkyl_boltzmann_photon_kernels.h> 
GKYL_CU_DH double boltzmann_photon_boundary_surfx_1x2v_tensor_p2(const double *w, const double *dxv,  double light_speed, double rho_curv, 
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
  const double dx10 = 2.0/dxv[0]; 
  const double wv = w[1]; 

  double sign_kpar = 1.0; 
  if (wv < 0.0) sign_kpar = -1.0; 
  double Ghat[9]; 
  if (edge == -1) { 
    if (wv>0) { 
      Ghat[0] = 1.58113883008419*fSkin[7]+1.224744871391589*fSkin[1]+0.7071067811865475*fSkin[0]; 
      Ghat[1] = 1.58113883008419*fSkin[11]+1.224744871391589*fSkin[4]+0.7071067811865475*fSkin[2]; 
      Ghat[2] = 1.58113883008419*fSkin[13]+1.224744871391589*fSkin[5]+0.7071067811865475*fSkin[3]; 
      Ghat[3] = 1.58113883008419*fSkin[17]+1.224744871391589*fSkin[10]+0.7071067811865475*fSkin[6]; 
      Ghat[4] = 1.58113883008419*fSkin[20]+1.224744871391589*fSkin[12]+0.7071067811865475*fSkin[8]; 
      Ghat[5] = 1.58113883008419*fSkin[21]+1.224744871391589*fSkin[15]+0.7071067811865475*fSkin[9]; 
      Ghat[6] = 1.58113883008419*fSkin[23]+1.224744871391589*fSkin[18]+0.7071067811865475*fSkin[14]; 
      Ghat[7] = 1.58113883008419*fSkin[24]+1.224744871391589*fSkin[19]+0.7071067811865475*fSkin[16]; 
      Ghat[8] = 1.58113883008419*fSkin[26]+1.224744871391589*fSkin[25]+0.7071067811865475*fSkin[22]; 
    } 
    else { 
      Ghat[0] = 1.58113883008419*fEdge[7]-1.224744871391589*fEdge[1]+0.7071067811865475*fEdge[0]; 
      Ghat[1] = 1.58113883008419*fEdge[11]-1.224744871391589*fEdge[4]+0.7071067811865475*fEdge[2]; 
      Ghat[2] = 1.58113883008419*fEdge[13]-1.224744871391589*fEdge[5]+0.7071067811865475*fEdge[3]; 
      Ghat[3] = 1.58113883008419*fEdge[17]-1.224744871391589*fEdge[10]+0.7071067811865475*fEdge[6]; 
      Ghat[4] = 1.58113883008419*fEdge[20]-1.224744871391589*fEdge[12]+0.7071067811865475*fEdge[8]; 
      Ghat[5] = 1.58113883008419*fEdge[21]-1.224744871391589*fEdge[15]+0.7071067811865475*fEdge[9]; 
      Ghat[6] = 1.58113883008419*fEdge[23]-1.224744871391589*fEdge[18]+0.7071067811865475*fEdge[14]; 
      Ghat[7] = 1.58113883008419*fEdge[24]-1.224744871391589*fEdge[19]+0.7071067811865475*fEdge[16]; 
      Ghat[8] = 1.58113883008419*fEdge[26]-1.224744871391589*fEdge[25]+0.7071067811865475*fEdge[22]; 
    } 
    out[0] += light_speed*sign_kpar*dx10*(-0.7071067811865475*Ghat[0]); 
    out[1] += light_speed*sign_kpar*dx10*(-1.224744871391589*Ghat[0]); 
    out[2] += light_speed*sign_kpar*dx10*(-0.7071067811865475*Ghat[1]); 
    out[3] += light_speed*sign_kpar*dx10*(-0.7071067811865475*Ghat[2]); 
    out[4] += light_speed*sign_kpar*dx10*(-1.224744871391589*Ghat[1]); 
    out[5] += light_speed*sign_kpar*dx10*(-1.224744871391589*Ghat[2]); 
    out[6] += light_speed*sign_kpar*dx10*(-0.7071067811865475*Ghat[3]); 
    out[7] += light_speed*sign_kpar*dx10*(-1.58113883008419*Ghat[0]); 
    out[8] += light_speed*sign_kpar*dx10*(-0.7071067811865475*Ghat[4]); 
  } 
  else { 
    if (wv>0) { 
      Ghat[0] = 1.58113883008419*fEdge[7]+1.224744871391589*fEdge[1]+0.7071067811865475*fEdge[0]; 
      Ghat[1] = 1.58113883008419*fEdge[11]+1.224744871391589*fEdge[4]+0.7071067811865475*fEdge[2]; 
      Ghat[2] = 1.58113883008419*fEdge[13]+1.224744871391589*fEdge[5]+0.7071067811865475*fEdge[3]; 
      Ghat[3] = 1.58113883008419*fEdge[17]+1.224744871391589*fEdge[10]+0.7071067811865475*fEdge[6]; 
      Ghat[4] = 1.58113883008419*fEdge[20]+1.224744871391589*fEdge[12]+0.7071067811865475*fEdge[8]; 
      Ghat[5] = 1.58113883008419*fEdge[21]+1.224744871391589*fEdge[15]+0.7071067811865475*fEdge[9]; 
      Ghat[6] = 1.58113883008419*fEdge[23]+1.224744871391589*fEdge[18]+0.7071067811865475*fEdge[14]; 
      Ghat[7] = 1.58113883008419*fEdge[24]+1.224744871391589*fEdge[19]+0.7071067811865475*fEdge[16]; 
      Ghat[8] = 1.58113883008419*fEdge[26]+1.224744871391589*fEdge[25]+0.7071067811865475*fEdge[22]; 
    } 
    else { 
      Ghat[0] = 1.58113883008419*fSkin[7]-1.224744871391589*fSkin[1]+0.7071067811865475*fSkin[0]; 
      Ghat[1] = 1.58113883008419*fSkin[11]-1.224744871391589*fSkin[4]+0.7071067811865475*fSkin[2]; 
      Ghat[2] = 1.58113883008419*fSkin[13]-1.224744871391589*fSkin[5]+0.7071067811865475*fSkin[3]; 
      Ghat[3] = 1.58113883008419*fSkin[17]-1.224744871391589*fSkin[10]+0.7071067811865475*fSkin[6]; 
      Ghat[4] = 1.58113883008419*fSkin[20]-1.224744871391589*fSkin[12]+0.7071067811865475*fSkin[8]; 
      Ghat[5] = 1.58113883008419*fSkin[21]-1.224744871391589*fSkin[15]+0.7071067811865475*fSkin[9]; 
      Ghat[6] = 1.58113883008419*fSkin[23]-1.224744871391589*fSkin[18]+0.7071067811865475*fSkin[14]; 
      Ghat[7] = 1.58113883008419*fSkin[24]-1.224744871391589*fSkin[19]+0.7071067811865475*fSkin[16]; 
      Ghat[8] = 1.58113883008419*fSkin[26]-1.224744871391589*fSkin[25]+0.7071067811865475*fSkin[22]; 
    } 
    out[0] += light_speed*sign_kpar*dx10*(0.7071067811865475*Ghat[0]); 
    out[1] += light_speed*sign_kpar*dx10*(-1.224744871391589*Ghat[0]); 
    out[2] += light_speed*sign_kpar*dx10*(0.7071067811865475*Ghat[1]); 
    out[3] += light_speed*sign_kpar*dx10*(0.7071067811865475*Ghat[2]); 
    out[4] += light_speed*sign_kpar*dx10*(-1.224744871391589*Ghat[1]); 
    out[5] += light_speed*sign_kpar*dx10*(-1.224744871391589*Ghat[2]); 
    out[6] += light_speed*sign_kpar*dx10*(0.7071067811865475*Ghat[3]); 
    out[7] += light_speed*sign_kpar*dx10*(1.58113883008419*Ghat[0]); 
    out[8] += light_speed*sign_kpar*dx10*(0.7071067811865475*Ghat[4]); 
  } 
  double cflFreq = light_speed; 
  return 2.5*dx10*cflFreq; 

} 
