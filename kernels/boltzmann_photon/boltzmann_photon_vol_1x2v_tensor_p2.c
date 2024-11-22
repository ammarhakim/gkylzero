#include <gkyl_boltzmann_photon_kernels.h> 
GKYL_CU_DH double boltzmann_photon_vol_1x2v_tensor_p2(const double *w, const double *dxv, double light_speed, double rho_curv, 
  const double *jacob_vel_inv, const double *kpar_abs, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:     Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // light_speed: Speed of light.
  // rho_curv:    Curvature of the magnetic field.
  // jacob_vel_inv[VDIM]: Inverse velocity space Jacobian in each direction.
  // kpar_abs:    Continuous expansion of |kpar| on the grid.
  // f:           Input photon distribution function.
  // out:         Incremented output.
  const double dx10 = 2.0/dxv[0]; 
  const double dv10 = 2.0/dxv[1]; 
  const double *jacob_vel_inv0 = &jacob_vel_inv[0]; 
  const double dv11 = 2.0/dxv[2]; 
  const double *jacob_vel_inv1 = &jacob_vel_inv[3]; 
  double incr_cdim[27] = {0.0}; 
  double incr_vdim[27] = {0.0}; 

  double sign_kpar = 1.0; 
  if (w[1] < 0.0) sign_kpar = -1.0; 
  incr_cdim[1] = light_speed*sign_kpar*dx10*(1.732050807568877*f[0]); 
  incr_cdim[4] = light_speed*sign_kpar*dx10*(1.732050807568877*f[2]); 
  incr_cdim[5] = light_speed*sign_kpar*dx10*(1.732050807568877*f[3]); 
  incr_cdim[7] = light_speed*sign_kpar*dx10*(3.872983346207417*f[1]); 
  incr_cdim[10] = light_speed*sign_kpar*dx10*(1.732050807568877*f[6]); 
  incr_cdim[11] = light_speed*sign_kpar*dx10*(3.872983346207417*f[4]); 
  incr_cdim[12] = light_speed*sign_kpar*dx10*(1.732050807568877*f[8]); 
  incr_cdim[13] = light_speed*sign_kpar*dx10*(3.872983346207417*f[5]); 
  incr_cdim[15] = light_speed*sign_kpar*dx10*(1.732050807568877*f[9]); 
  incr_cdim[17] = light_speed*sign_kpar*dx10*(3.872983346207417*f[10]); 
  incr_cdim[18] = light_speed*sign_kpar*dx10*(1.732050807568877*f[14]); 
  incr_cdim[19] = light_speed*sign_kpar*dx10*(1.732050807568877*f[16]); 
  incr_cdim[20] = light_speed*sign_kpar*dx10*(3.872983346207417*f[12]); 
  incr_cdim[21] = light_speed*sign_kpar*dx10*(3.872983346207417*f[15]); 
  incr_cdim[23] = light_speed*sign_kpar*dx10*(3.872983346207417*f[18]); 
  incr_cdim[24] = light_speed*sign_kpar*dx10*(3.872983346207417*f[19]); 
  incr_cdim[25] = light_speed*sign_kpar*dx10*(1.732050807568877*f[22]); 
  incr_cdim[26] = light_speed*sign_kpar*dx10*(3.872983346207417*f[25]); 

  incr_vdim[3] = light_speed/rho_curv*dv11*(0.8660254037844386*jacob_vel_inv1[2]*kpar_abs[2]*f[22]+0.8660254037844387*kpar_abs[1]*jacob_vel_inv1[2]*f[16]+0.8660254037844387*jacob_vel_inv1[1]*kpar_abs[2]*f[14]+0.8660254037844386*kpar_abs[0]*jacob_vel_inv1[2]*f[9]+0.8660254037844386*jacob_vel_inv1[0]*kpar_abs[2]*f[8]+0.8660254037844386*jacob_vel_inv1[1]*kpar_abs[1]*f[6]+0.8660254037844386*kpar_abs[0]*jacob_vel_inv1[1]*f[3]+0.8660254037844386*jacob_vel_inv1[0]*kpar_abs[1]*f[2]+0.8660254037844386*f[0]*jacob_vel_inv1[0]*kpar_abs[0]); 
  incr_vdim[5] = light_speed/rho_curv*dv11*(0.8660254037844386*jacob_vel_inv1[2]*kpar_abs[2]*f[25]+0.8660254037844386*kpar_abs[1]*jacob_vel_inv1[2]*f[19]+0.8660254037844386*jacob_vel_inv1[1]*kpar_abs[2]*f[18]+0.8660254037844387*kpar_abs[0]*jacob_vel_inv1[2]*f[15]+0.8660254037844387*jacob_vel_inv1[0]*kpar_abs[2]*f[12]+0.8660254037844386*jacob_vel_inv1[1]*kpar_abs[1]*f[10]+0.8660254037844386*kpar_abs[0]*jacob_vel_inv1[1]*f[5]+0.8660254037844386*jacob_vel_inv1[0]*kpar_abs[1]*f[4]+0.8660254037844386*jacob_vel_inv1[0]*kpar_abs[0]*f[1]); 
  incr_vdim[6] = light_speed/rho_curv*dv11*(0.7745966692414833*kpar_abs[1]*jacob_vel_inv1[2]*f[22]+0.7745966692414834*jacob_vel_inv1[2]*kpar_abs[2]*f[16]+0.8660254037844387*kpar_abs[0]*jacob_vel_inv1[2]*f[16]+0.7745966692414834*jacob_vel_inv1[1]*kpar_abs[1]*f[14]+0.8660254037844386*kpar_abs[1]*jacob_vel_inv1[2]*f[9]+0.7745966692414833*jacob_vel_inv1[0]*kpar_abs[1]*f[8]+0.7745966692414833*jacob_vel_inv1[1]*kpar_abs[2]*f[6]+0.8660254037844386*kpar_abs[0]*jacob_vel_inv1[1]*f[6]+0.8660254037844386*jacob_vel_inv1[1]*kpar_abs[1]*f[3]+0.7745966692414833*jacob_vel_inv1[0]*f[2]*kpar_abs[2]+0.8660254037844386*jacob_vel_inv1[0]*kpar_abs[0]*f[2]+0.8660254037844386*f[0]*jacob_vel_inv1[0]*kpar_abs[1]); 
  incr_vdim[9] = light_speed/rho_curv*dv11*(1.732050807568877*jacob_vel_inv1[1]*kpar_abs[2]*f[22]+1.732050807568877*jacob_vel_inv1[1]*kpar_abs[1]*f[16]+1.732050807568877*jacob_vel_inv1[2]*kpar_abs[2]*f[14]+1.936491673103709*jacob_vel_inv1[0]*kpar_abs[2]*f[14]+1.732050807568877*kpar_abs[0]*jacob_vel_inv1[1]*f[9]+1.936491673103709*jacob_vel_inv1[1]*kpar_abs[2]*f[8]+1.732050807568877*kpar_abs[1]*jacob_vel_inv1[2]*f[6]+1.936491673103709*jacob_vel_inv1[0]*kpar_abs[1]*f[6]+1.732050807568877*kpar_abs[0]*jacob_vel_inv1[2]*f[3]+1.936491673103709*jacob_vel_inv1[0]*kpar_abs[0]*f[3]+1.936491673103709*jacob_vel_inv1[1]*kpar_abs[1]*f[2]+1.936491673103709*f[0]*kpar_abs[0]*jacob_vel_inv1[1]); 
  incr_vdim[10] = light_speed/rho_curv*dv11*(0.7745966692414833*kpar_abs[1]*jacob_vel_inv1[2]*f[25]+0.7745966692414833*jacob_vel_inv1[2]*kpar_abs[2]*f[19]+0.8660254037844386*kpar_abs[0]*jacob_vel_inv1[2]*f[19]+0.7745966692414833*jacob_vel_inv1[1]*kpar_abs[1]*f[18]+0.8660254037844387*kpar_abs[1]*jacob_vel_inv1[2]*f[15]+0.7745966692414834*jacob_vel_inv1[0]*kpar_abs[1]*f[12]+0.7745966692414833*jacob_vel_inv1[1]*kpar_abs[2]*f[10]+0.8660254037844386*kpar_abs[0]*jacob_vel_inv1[1]*f[10]+0.8660254037844386*jacob_vel_inv1[1]*kpar_abs[1]*f[5]+0.7745966692414833*jacob_vel_inv1[0]*kpar_abs[2]*f[4]+0.8660254037844386*jacob_vel_inv1[0]*kpar_abs[0]*f[4]+0.8660254037844386*jacob_vel_inv1[0]*f[1]*kpar_abs[1]); 
  incr_vdim[13] = light_speed/rho_curv*dv11*(0.8660254037844387*jacob_vel_inv1[2]*kpar_abs[2]*f[26]+0.8660254037844387*kpar_abs[1]*jacob_vel_inv1[2]*f[24]+0.8660254037844387*jacob_vel_inv1[1]*kpar_abs[2]*f[23]+0.8660254037844387*kpar_abs[0]*jacob_vel_inv1[2]*f[21]+0.8660254037844387*jacob_vel_inv1[0]*kpar_abs[2]*f[20]+0.8660254037844387*jacob_vel_inv1[1]*kpar_abs[1]*f[17]+0.8660254037844386*kpar_abs[0]*jacob_vel_inv1[1]*f[13]+0.8660254037844386*jacob_vel_inv1[0]*kpar_abs[1]*f[11]+0.8660254037844387*jacob_vel_inv1[0]*kpar_abs[0]*f[7]); 
  incr_vdim[14] = light_speed/rho_curv*dv11*(0.5532833351724881*jacob_vel_inv1[2]*kpar_abs[2]*f[22]+0.8660254037844387*kpar_abs[0]*jacob_vel_inv1[2]*f[22]+0.7745966692414833*kpar_abs[1]*jacob_vel_inv1[2]*f[16]+0.5532833351724881*jacob_vel_inv1[1]*kpar_abs[2]*f[14]+0.8660254037844386*kpar_abs[0]*jacob_vel_inv1[1]*f[14]+0.8660254037844387*jacob_vel_inv1[2]*kpar_abs[2]*f[9]+0.5532833351724881*jacob_vel_inv1[0]*kpar_abs[2]*f[8]+0.8660254037844387*jacob_vel_inv1[0]*kpar_abs[0]*f[8]+0.7745966692414834*jacob_vel_inv1[1]*kpar_abs[1]*f[6]+0.8660254037844387*jacob_vel_inv1[1]*kpar_abs[2]*f[3]+0.8660254037844387*f[0]*jacob_vel_inv1[0]*kpar_abs[2]+0.7745966692414834*jacob_vel_inv1[0]*kpar_abs[1]*f[2]); 
  incr_vdim[15] = light_speed/rho_curv*dv11*(1.732050807568877*jacob_vel_inv1[1]*kpar_abs[2]*f[25]+1.732050807568877*jacob_vel_inv1[1]*kpar_abs[1]*f[19]+1.732050807568877*jacob_vel_inv1[2]*kpar_abs[2]*f[18]+1.936491673103709*jacob_vel_inv1[0]*kpar_abs[2]*f[18]+1.732050807568877*kpar_abs[0]*jacob_vel_inv1[1]*f[15]+1.936491673103709*jacob_vel_inv1[1]*kpar_abs[2]*f[12]+1.732050807568877*kpar_abs[1]*jacob_vel_inv1[2]*f[10]+1.936491673103709*jacob_vel_inv1[0]*kpar_abs[1]*f[10]+1.732050807568877*kpar_abs[0]*jacob_vel_inv1[2]*f[5]+1.936491673103709*jacob_vel_inv1[0]*kpar_abs[0]*f[5]+1.936491673103709*jacob_vel_inv1[1]*kpar_abs[1]*f[4]+1.936491673103709*kpar_abs[0]*f[1]*jacob_vel_inv1[1]); 
  incr_vdim[16] = light_speed/rho_curv*dv11*(1.549193338482967*jacob_vel_inv1[1]*kpar_abs[1]*f[22]+1.549193338482967*jacob_vel_inv1[1]*kpar_abs[2]*f[16]+1.732050807568877*kpar_abs[0]*jacob_vel_inv1[1]*f[16]+1.549193338482967*kpar_abs[1]*jacob_vel_inv1[2]*f[14]+1.732050807568877*jacob_vel_inv1[0]*kpar_abs[1]*f[14]+1.732050807568877*jacob_vel_inv1[1]*kpar_abs[1]*f[9]+1.732050807568877*jacob_vel_inv1[1]*kpar_abs[1]*f[8]+1.549193338482967*jacob_vel_inv1[2]*kpar_abs[2]*f[6]+1.732050807568877*jacob_vel_inv1[0]*kpar_abs[2]*f[6]+1.732050807568877*kpar_abs[0]*jacob_vel_inv1[2]*f[6]+1.936491673103709*jacob_vel_inv1[0]*kpar_abs[0]*f[6]+1.732050807568877*kpar_abs[1]*jacob_vel_inv1[2]*f[3]+1.936491673103709*jacob_vel_inv1[0]*kpar_abs[1]*f[3]+1.732050807568877*jacob_vel_inv1[1]*f[2]*kpar_abs[2]+1.936491673103709*kpar_abs[0]*jacob_vel_inv1[1]*f[2]+1.936491673103709*f[0]*jacob_vel_inv1[1]*kpar_abs[1]); 
  incr_vdim[17] = light_speed/rho_curv*dv11*(0.7745966692414833*kpar_abs[1]*jacob_vel_inv1[2]*f[26]+0.7745966692414833*jacob_vel_inv1[2]*kpar_abs[2]*f[24]+0.8660254037844386*kpar_abs[0]*jacob_vel_inv1[2]*f[24]+0.7745966692414833*jacob_vel_inv1[1]*kpar_abs[1]*f[23]+0.8660254037844386*kpar_abs[1]*jacob_vel_inv1[2]*f[21]+0.7745966692414833*jacob_vel_inv1[0]*kpar_abs[1]*f[20]+0.7745966692414833*jacob_vel_inv1[1]*kpar_abs[2]*f[17]+0.8660254037844386*kpar_abs[0]*jacob_vel_inv1[1]*f[17]+0.8660254037844387*jacob_vel_inv1[1]*kpar_abs[1]*f[13]+0.7745966692414834*jacob_vel_inv1[0]*kpar_abs[2]*f[11]+0.8660254037844387*jacob_vel_inv1[0]*kpar_abs[0]*f[11]+0.8660254037844386*jacob_vel_inv1[0]*kpar_abs[1]*f[7]); 
  incr_vdim[18] = light_speed/rho_curv*dv11*(0.5532833351724881*jacob_vel_inv1[2]*kpar_abs[2]*f[25]+0.8660254037844386*kpar_abs[0]*jacob_vel_inv1[2]*f[25]+0.7745966692414833*kpar_abs[1]*jacob_vel_inv1[2]*f[19]+0.5532833351724881*jacob_vel_inv1[1]*kpar_abs[2]*f[18]+0.8660254037844386*kpar_abs[0]*jacob_vel_inv1[1]*f[18]+0.8660254037844387*jacob_vel_inv1[2]*kpar_abs[2]*f[15]+0.5532833351724881*jacob_vel_inv1[0]*kpar_abs[2]*f[12]+0.8660254037844387*jacob_vel_inv1[0]*kpar_abs[0]*f[12]+0.7745966692414833*jacob_vel_inv1[1]*kpar_abs[1]*f[10]+0.8660254037844386*jacob_vel_inv1[1]*kpar_abs[2]*f[5]+0.7745966692414833*jacob_vel_inv1[0]*kpar_abs[1]*f[4]+0.8660254037844386*jacob_vel_inv1[0]*f[1]*kpar_abs[2]); 
  incr_vdim[19] = light_speed/rho_curv*dv11*(1.549193338482967*jacob_vel_inv1[1]*kpar_abs[1]*f[25]+1.549193338482967*jacob_vel_inv1[1]*kpar_abs[2]*f[19]+1.732050807568877*kpar_abs[0]*jacob_vel_inv1[1]*f[19]+1.549193338482967*kpar_abs[1]*jacob_vel_inv1[2]*f[18]+1.732050807568877*jacob_vel_inv1[0]*kpar_abs[1]*f[18]+1.732050807568877*jacob_vel_inv1[1]*kpar_abs[1]*f[15]+1.732050807568877*jacob_vel_inv1[1]*kpar_abs[1]*f[12]+1.549193338482967*jacob_vel_inv1[2]*kpar_abs[2]*f[10]+1.732050807568877*jacob_vel_inv1[0]*kpar_abs[2]*f[10]+1.732050807568877*kpar_abs[0]*jacob_vel_inv1[2]*f[10]+1.936491673103709*jacob_vel_inv1[0]*kpar_abs[0]*f[10]+1.732050807568877*kpar_abs[1]*jacob_vel_inv1[2]*f[5]+1.936491673103709*jacob_vel_inv1[0]*kpar_abs[1]*f[5]+1.732050807568877*jacob_vel_inv1[1]*kpar_abs[2]*f[4]+1.936491673103709*kpar_abs[0]*jacob_vel_inv1[1]*f[4]+1.936491673103709*f[1]*jacob_vel_inv1[1]*kpar_abs[1]); 
  incr_vdim[21] = light_speed/rho_curv*dv11*(1.732050807568877*jacob_vel_inv1[1]*kpar_abs[2]*f[26]+1.732050807568877*jacob_vel_inv1[1]*kpar_abs[1]*f[24]+1.732050807568877*jacob_vel_inv1[2]*kpar_abs[2]*f[23]+1.936491673103709*jacob_vel_inv1[0]*kpar_abs[2]*f[23]+1.732050807568877*kpar_abs[0]*jacob_vel_inv1[1]*f[21]+1.936491673103709*jacob_vel_inv1[1]*kpar_abs[2]*f[20]+1.732050807568877*kpar_abs[1]*jacob_vel_inv1[2]*f[17]+1.936491673103709*jacob_vel_inv1[0]*kpar_abs[1]*f[17]+1.732050807568877*kpar_abs[0]*jacob_vel_inv1[2]*f[13]+1.936491673103709*jacob_vel_inv1[0]*kpar_abs[0]*f[13]+1.936491673103709*jacob_vel_inv1[1]*kpar_abs[1]*f[11]+1.936491673103709*kpar_abs[0]*jacob_vel_inv1[1]*f[7]); 
  incr_vdim[22] = light_speed/rho_curv*dv11*(1.106566670344976*jacob_vel_inv1[1]*kpar_abs[2]*f[22]+1.732050807568877*kpar_abs[0]*jacob_vel_inv1[1]*f[22]+1.549193338482967*jacob_vel_inv1[1]*kpar_abs[1]*f[16]+1.106566670344976*jacob_vel_inv1[2]*kpar_abs[2]*f[14]+1.237179148263484*jacob_vel_inv1[0]*kpar_abs[2]*f[14]+1.732050807568877*kpar_abs[0]*jacob_vel_inv1[2]*f[14]+1.936491673103709*jacob_vel_inv1[0]*kpar_abs[0]*f[14]+1.732050807568877*jacob_vel_inv1[1]*kpar_abs[2]*f[9]+1.237179148263484*jacob_vel_inv1[1]*kpar_abs[2]*f[8]+1.936491673103709*kpar_abs[0]*jacob_vel_inv1[1]*f[8]+1.549193338482967*kpar_abs[1]*jacob_vel_inv1[2]*f[6]+1.732050807568877*jacob_vel_inv1[0]*kpar_abs[1]*f[6]+1.732050807568877*jacob_vel_inv1[2]*kpar_abs[2]*f[3]+1.936491673103709*jacob_vel_inv1[0]*kpar_abs[2]*f[3]+1.936491673103709*f[0]*jacob_vel_inv1[1]*kpar_abs[2]+1.732050807568877*jacob_vel_inv1[1]*kpar_abs[1]*f[2]); 
  incr_vdim[23] = light_speed/rho_curv*dv11*(0.5532833351724881*jacob_vel_inv1[2]*kpar_abs[2]*f[26]+0.8660254037844386*kpar_abs[0]*jacob_vel_inv1[2]*f[26]+0.7745966692414833*kpar_abs[1]*jacob_vel_inv1[2]*f[24]+0.5532833351724881*jacob_vel_inv1[1]*kpar_abs[2]*f[23]+0.8660254037844386*kpar_abs[0]*jacob_vel_inv1[1]*f[23]+0.8660254037844386*jacob_vel_inv1[2]*kpar_abs[2]*f[21]+0.5532833351724881*jacob_vel_inv1[0]*kpar_abs[2]*f[20]+0.8660254037844386*jacob_vel_inv1[0]*kpar_abs[0]*f[20]+0.7745966692414833*jacob_vel_inv1[1]*kpar_abs[1]*f[17]+0.8660254037844387*jacob_vel_inv1[1]*kpar_abs[2]*f[13]+0.7745966692414834*jacob_vel_inv1[0]*kpar_abs[1]*f[11]+0.8660254037844386*jacob_vel_inv1[0]*kpar_abs[2]*f[7]); 
  incr_vdim[24] = light_speed/rho_curv*dv11*(1.549193338482967*jacob_vel_inv1[1]*kpar_abs[1]*f[26]+1.549193338482967*jacob_vel_inv1[1]*kpar_abs[2]*f[24]+1.732050807568877*kpar_abs[0]*jacob_vel_inv1[1]*f[24]+1.549193338482967*kpar_abs[1]*jacob_vel_inv1[2]*f[23]+1.732050807568877*jacob_vel_inv1[0]*kpar_abs[1]*f[23]+1.732050807568877*jacob_vel_inv1[1]*kpar_abs[1]*f[21]+1.732050807568877*jacob_vel_inv1[1]*kpar_abs[1]*f[20]+1.549193338482967*jacob_vel_inv1[2]*kpar_abs[2]*f[17]+1.732050807568877*jacob_vel_inv1[0]*kpar_abs[2]*f[17]+1.732050807568877*kpar_abs[0]*jacob_vel_inv1[2]*f[17]+1.936491673103709*jacob_vel_inv1[0]*kpar_abs[0]*f[17]+1.732050807568877*kpar_abs[1]*jacob_vel_inv1[2]*f[13]+1.936491673103709*jacob_vel_inv1[0]*kpar_abs[1]*f[13]+1.732050807568877*jacob_vel_inv1[1]*kpar_abs[2]*f[11]+1.936491673103709*kpar_abs[0]*jacob_vel_inv1[1]*f[11]+1.936491673103709*jacob_vel_inv1[1]*kpar_abs[1]*f[7]); 
  incr_vdim[25] = light_speed/rho_curv*dv11*(1.106566670344976*jacob_vel_inv1[1]*kpar_abs[2]*f[25]+1.732050807568877*kpar_abs[0]*jacob_vel_inv1[1]*f[25]+1.549193338482967*jacob_vel_inv1[1]*kpar_abs[1]*f[19]+1.106566670344976*jacob_vel_inv1[2]*kpar_abs[2]*f[18]+1.237179148263484*jacob_vel_inv1[0]*kpar_abs[2]*f[18]+1.732050807568877*kpar_abs[0]*jacob_vel_inv1[2]*f[18]+1.936491673103709*jacob_vel_inv1[0]*kpar_abs[0]*f[18]+1.732050807568877*jacob_vel_inv1[1]*kpar_abs[2]*f[15]+1.237179148263484*jacob_vel_inv1[1]*kpar_abs[2]*f[12]+1.936491673103709*kpar_abs[0]*jacob_vel_inv1[1]*f[12]+1.549193338482967*kpar_abs[1]*jacob_vel_inv1[2]*f[10]+1.732050807568877*jacob_vel_inv1[0]*kpar_abs[1]*f[10]+1.732050807568877*jacob_vel_inv1[2]*kpar_abs[2]*f[5]+1.936491673103709*jacob_vel_inv1[0]*kpar_abs[2]*f[5]+1.732050807568877*jacob_vel_inv1[1]*kpar_abs[1]*f[4]+1.936491673103709*f[1]*jacob_vel_inv1[1]*kpar_abs[2]); 
  incr_vdim[26] = light_speed/rho_curv*dv11*(1.106566670344976*jacob_vel_inv1[1]*kpar_abs[2]*f[26]+1.732050807568877*kpar_abs[0]*jacob_vel_inv1[1]*f[26]+1.549193338482967*jacob_vel_inv1[1]*kpar_abs[1]*f[24]+1.106566670344976*jacob_vel_inv1[2]*kpar_abs[2]*f[23]+1.237179148263484*jacob_vel_inv1[0]*kpar_abs[2]*f[23]+1.732050807568877*kpar_abs[0]*jacob_vel_inv1[2]*f[23]+1.936491673103709*jacob_vel_inv1[0]*kpar_abs[0]*f[23]+1.732050807568877*jacob_vel_inv1[1]*kpar_abs[2]*f[21]+1.237179148263484*jacob_vel_inv1[1]*kpar_abs[2]*f[20]+1.936491673103709*kpar_abs[0]*jacob_vel_inv1[1]*f[20]+1.549193338482967*kpar_abs[1]*jacob_vel_inv1[2]*f[17]+1.732050807568877*jacob_vel_inv1[0]*kpar_abs[1]*f[17]+1.732050807568877*jacob_vel_inv1[2]*kpar_abs[2]*f[13]+1.936491673103709*jacob_vel_inv1[0]*kpar_abs[2]*f[13]+1.732050807568877*jacob_vel_inv1[1]*kpar_abs[1]*f[11]+1.936491673103709*jacob_vel_inv1[1]*kpar_abs[2]*f[7]); 

  out[0] += incr_cdim[0] + incr_vdim[0]; 
  out[1] += incr_cdim[1] + incr_vdim[1]; 
  out[2] += incr_cdim[2] + incr_vdim[2]; 
  out[3] += incr_cdim[3] + incr_vdim[3]; 
  out[4] += incr_cdim[4] + incr_vdim[4]; 
  out[5] += incr_cdim[5] + incr_vdim[5]; 
  out[6] += incr_cdim[6] + incr_vdim[6]; 
  out[7] += incr_cdim[7] + incr_vdim[7]; 
  out[8] += incr_cdim[8] + incr_vdim[8]; 
  out[9] += incr_cdim[9] + incr_vdim[9]; 
  out[10] += incr_cdim[10] + incr_vdim[10]; 
  out[11] += incr_cdim[11] + incr_vdim[11]; 
  out[12] += incr_cdim[12] + incr_vdim[12]; 
  out[13] += incr_cdim[13] + incr_vdim[13]; 
  out[14] += incr_cdim[14] + incr_vdim[14]; 
  out[15] += incr_cdim[15] + incr_vdim[15]; 
  out[16] += incr_cdim[16] + incr_vdim[16]; 
  out[17] += incr_cdim[17] + incr_vdim[17]; 
  out[18] += incr_cdim[18] + incr_vdim[18]; 
  out[19] += incr_cdim[19] + incr_vdim[19]; 
  out[20] += incr_cdim[20] + incr_vdim[20]; 
  out[21] += incr_cdim[21] + incr_vdim[21]; 
  out[22] += incr_cdim[22] + incr_vdim[22]; 
  out[23] += incr_cdim[23] + incr_vdim[23]; 
  out[24] += incr_cdim[24] + incr_vdim[24]; 
  out[25] += incr_cdim[25] + incr_vdim[25]; 
  out[26] += incr_cdim[26] + incr_vdim[26]; 

  return 0.0; 
} 
