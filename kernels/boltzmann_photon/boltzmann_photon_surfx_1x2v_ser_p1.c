#include <gkyl_boltzmann_photon_kernels.h> 
GKYL_CU_DH double boltzmann_photon_surfx_1x2v_ser_p1(const double *w, const double *dxv,  double light_speed, double rho_curv, 
  const double *jacob_vel_inv, const double *kpar_abs, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:     Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // light_speed: Speed of light.
  // rho_cuv:     Curvature of the magnetic field.
  // jacob_vel_inv[VDIM]: Inverse velocity space Jacobian in each direction.
  // kpar_abs:    Continuous expansion of |kpar| on the grid.
  // fl/fc/fr:    Input photon distribution function in left/center/right cells.
  // out:         Incremented photon distribution function in center cell.
  const double dx10 = 2.0/dxv[0]; 
  const double dv10 = 2.0/dxv[1]; 
  const double wv = w[1]; 

  double sign_kpar = 1.0; 
  if (wv < 0.0) sign_kpar = -1.0; 
  double Ghat_r[8]; 
  double Ghat_l[8]; 
  if (wv>0) { 

  Ghat_r[0] = 1.224744871391589*fc[1]+0.7071067811865475*fc[0]; 
  Ghat_r[1] = 1.224744871391589*fc[4]+0.7071067811865475*fc[2]; 
  Ghat_r[2] = 1.224744871391589*fc[5]+0.7071067811865475*fc[3]; 
  Ghat_r[3] = 1.224744871391589*fc[7]+0.7071067811865475*fc[6]; 
  Ghat_r[4] = 1.224744871391589*fc[9]+0.7071067811865475*fc[8]; 
  Ghat_r[5] = 1.224744871391589*fc[13]+0.7071067811865475*fc[12]; 
  Ghat_r[6] = 1.224744871391589*fc[11]+0.7071067811865475*fc[10]; 
  Ghat_r[7] = 1.224744871391589*fc[15]+0.7071067811865475*fc[14]; 

  Ghat_l[0] = 1.224744871391589*fl[1]+0.7071067811865475*fl[0]; 
  Ghat_l[1] = 1.224744871391589*fl[4]+0.7071067811865475*fl[2]; 
  Ghat_l[2] = 1.224744871391589*fl[5]+0.7071067811865475*fl[3]; 
  Ghat_l[3] = 1.224744871391589*fl[7]+0.7071067811865475*fl[6]; 
  Ghat_l[4] = 1.224744871391589*fl[9]+0.7071067811865475*fl[8]; 
  Ghat_l[5] = 1.224744871391589*fl[13]+0.7071067811865475*fl[12]; 
  Ghat_l[6] = 1.224744871391589*fl[11]+0.7071067811865475*fl[10]; 
  Ghat_l[7] = 1.224744871391589*fl[15]+0.7071067811865475*fl[14]; 

  } else { 

  Ghat_r[0] = 0.7071067811865475*fr[0]-1.224744871391589*fr[1]; 
  Ghat_r[1] = 0.7071067811865475*fr[2]-1.224744871391589*fr[4]; 
  Ghat_r[2] = 0.7071067811865475*fr[3]-1.224744871391589*fr[5]; 
  Ghat_r[3] = 0.7071067811865475*fr[6]-1.224744871391589*fr[7]; 
  Ghat_r[4] = 0.7071067811865475*fr[8]-1.224744871391589*fr[9]; 
  Ghat_r[5] = 0.7071067811865475*fr[12]-1.224744871391589*fr[13]; 
  Ghat_r[6] = 0.7071067811865475*fr[10]-1.224744871391589*fr[11]; 
  Ghat_r[7] = 0.7071067811865475*fr[14]-1.224744871391589*fr[15]; 

  Ghat_l[0] = 0.7071067811865475*fc[0]-1.224744871391589*fc[1]; 
  Ghat_l[1] = 0.7071067811865475*fc[2]-1.224744871391589*fc[4]; 
  Ghat_l[2] = 0.7071067811865475*fc[3]-1.224744871391589*fc[5]; 
  Ghat_l[3] = 0.7071067811865475*fc[6]-1.224744871391589*fc[7]; 
  Ghat_l[4] = 0.7071067811865475*fc[8]-1.224744871391589*fc[9]; 
  Ghat_l[5] = 0.7071067811865475*fc[12]-1.224744871391589*fc[13]; 
  Ghat_l[6] = 0.7071067811865475*fc[10]-1.224744871391589*fc[11]; 
  Ghat_l[7] = 0.7071067811865475*fc[14]-1.224744871391589*fc[15]; 

  } 
  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dx10*light_speed*sign_kpar; 
  out[1] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dx10*light_speed*sign_kpar; 
  out[2] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dx10*light_speed*sign_kpar; 
  out[3] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dx10*light_speed*sign_kpar; 
  out[4] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dx10*light_speed*sign_kpar; 
  out[5] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dx10*light_speed*sign_kpar; 
  out[6] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dx10*light_speed*sign_kpar; 
  out[7] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dx10*light_speed*sign_kpar; 
  out[8] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*dx10*light_speed*sign_kpar; 
  out[9] += -1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dx10*light_speed*sign_kpar; 
  out[10] += (0.7071067811865475*Ghat_l[6]-0.7071067811865475*Ghat_r[6])*dx10*light_speed*sign_kpar; 
  out[11] += -1.224744871391589*(Ghat_r[6]+Ghat_l[6])*dx10*light_speed*sign_kpar; 
  out[12] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*dx10*light_speed*sign_kpar; 
  out[13] += -1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dx10*light_speed*sign_kpar; 
  out[14] += (0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*dx10*light_speed*sign_kpar; 
  out[15] += -1.224744871391589*(Ghat_r[7]+Ghat_l[7])*dx10*light_speed*sign_kpar; 

  double cflFreq = light_speed; 
  return 1.5*dx10*cflFreq; 

} 
