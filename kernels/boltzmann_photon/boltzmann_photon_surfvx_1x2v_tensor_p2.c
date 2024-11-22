#include <gkyl_boltzmann_photon_kernels.h> 
GKYL_CU_DH double boltzmann_photon_surfvx_1x2v_tensor_p2(const double *w, const double *dxv,  double light_speed, double rho_curv, 
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
  const double *jacob_vel_inv_dir = &jacob_vel_inv[0]; 
  return 0.0; 

} 
