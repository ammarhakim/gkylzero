#include <gkyl_boltzmann_photon_kernels.h> 
GKYL_CU_DH double boltzmann_photon_boundary_surfx_1x2v_ser_p1(const double *w, const double *dxv,  double light_speed, double rho_curv, 
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
  double Ghat[8]; 
  if (edge == -1) { 
    if (wv>0) { 
      Ghat[0] = 1.224744871391589*fSkin[1]+0.7071067811865475*fSkin[0]; 
      Ghat[1] = 1.224744871391589*fSkin[4]+0.7071067811865475*fSkin[2]; 
      Ghat[2] = 1.224744871391589*fSkin[5]+0.7071067811865475*fSkin[3]; 
      Ghat[3] = 1.224744871391589*fSkin[7]+0.7071067811865475*fSkin[6]; 
      Ghat[4] = 1.224744871391589*fSkin[9]+0.7071067811865475*fSkin[8]; 
      Ghat[5] = 1.224744871391589*fSkin[13]+0.7071067811865475*fSkin[12]; 
      Ghat[6] = 1.224744871391589*fSkin[11]+0.7071067811865475*fSkin[10]; 
      Ghat[7] = 1.224744871391589*fSkin[15]+0.7071067811865475*fSkin[14]; 
    } 
    else { 
      Ghat[0] = 0.7071067811865475*fEdge[0]-1.224744871391589*fEdge[1]; 
      Ghat[1] = 0.7071067811865475*fEdge[2]-1.224744871391589*fEdge[4]; 
      Ghat[2] = 0.7071067811865475*fEdge[3]-1.224744871391589*fEdge[5]; 
      Ghat[3] = 0.7071067811865475*fEdge[6]-1.224744871391589*fEdge[7]; 
      Ghat[4] = 0.7071067811865475*fEdge[8]-1.224744871391589*fEdge[9]; 
      Ghat[5] = 0.7071067811865475*fEdge[12]-1.224744871391589*fEdge[13]; 
      Ghat[6] = 0.7071067811865475*fEdge[10]-1.224744871391589*fEdge[11]; 
      Ghat[7] = 0.7071067811865475*fEdge[14]-1.224744871391589*fEdge[15]; 
    } 
    out[0] += light_speed*sign_kpar*dx10*(-0.7071067811865475*Ghat[0]); 
    out[1] += light_speed*sign_kpar*dx10*(-1.224744871391589*Ghat[0]); 
    out[2] += light_speed*sign_kpar*dx10*(-0.7071067811865475*Ghat[1]); 
    out[3] += light_speed*sign_kpar*dx10*(-0.7071067811865475*Ghat[2]); 
    out[4] += light_speed*sign_kpar*dx10*(-1.224744871391589*Ghat[1]); 
    out[5] += light_speed*sign_kpar*dx10*(-1.224744871391589*Ghat[2]); 
    out[6] += light_speed*sign_kpar*dx10*(-0.7071067811865475*Ghat[3]); 
    out[7] += light_speed*sign_kpar*dx10*(-1.224744871391589*Ghat[3]); 
  } 
  else { 
    if (wv>0) { 
      Ghat[0] = 1.224744871391589*fEdge[1]+0.7071067811865475*fEdge[0]; 
      Ghat[1] = 1.224744871391589*fEdge[4]+0.7071067811865475*fEdge[2]; 
      Ghat[2] = 1.224744871391589*fEdge[5]+0.7071067811865475*fEdge[3]; 
      Ghat[3] = 1.224744871391589*fEdge[7]+0.7071067811865475*fEdge[6]; 
      Ghat[4] = 1.224744871391589*fEdge[9]+0.7071067811865475*fEdge[8]; 
      Ghat[5] = 1.224744871391589*fEdge[13]+0.7071067811865475*fEdge[12]; 
      Ghat[6] = 1.224744871391589*fEdge[11]+0.7071067811865475*fEdge[10]; 
      Ghat[7] = 1.224744871391589*fEdge[15]+0.7071067811865475*fEdge[14]; 
    } 
    else { 
      Ghat[0] = 0.7071067811865475*fSkin[0]-1.224744871391589*fSkin[1]; 
      Ghat[1] = 0.7071067811865475*fSkin[2]-1.224744871391589*fSkin[4]; 
      Ghat[2] = 0.7071067811865475*fSkin[3]-1.224744871391589*fSkin[5]; 
      Ghat[3] = 0.7071067811865475*fSkin[6]-1.224744871391589*fSkin[7]; 
      Ghat[4] = 0.7071067811865475*fSkin[8]-1.224744871391589*fSkin[9]; 
      Ghat[5] = 0.7071067811865475*fSkin[12]-1.224744871391589*fSkin[13]; 
      Ghat[6] = 0.7071067811865475*fSkin[10]-1.224744871391589*fSkin[11]; 
      Ghat[7] = 0.7071067811865475*fSkin[14]-1.224744871391589*fSkin[15]; 
    } 
    out[0] += light_speed*sign_kpar*dx10*(0.7071067811865475*Ghat[0]); 
    out[1] += light_speed*sign_kpar*dx10*(-1.224744871391589*Ghat[0]); 
    out[2] += light_speed*sign_kpar*dx10*(0.7071067811865475*Ghat[1]); 
    out[3] += light_speed*sign_kpar*dx10*(0.7071067811865475*Ghat[2]); 
    out[4] += light_speed*sign_kpar*dx10*(-1.224744871391589*Ghat[1]); 
    out[5] += light_speed*sign_kpar*dx10*(-1.224744871391589*Ghat[2]); 
    out[6] += light_speed*sign_kpar*dx10*(0.7071067811865475*Ghat[3]); 
    out[7] += light_speed*sign_kpar*dx10*(-1.224744871391589*Ghat[3]); 
  } 
  double cflFreq = light_speed; 
  return 1.5*dx10*cflFreq; 

} 
