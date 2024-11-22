#include <gkyl_boltzmann_photon_kernels.h> 
GKYL_CU_DH double boltzmann_photon_boundary_surfvx_1x2v_ser_p1(const double *w, const double *dxv,  double light_speed, double rho_curv, 
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
  const double *jacob_vel_inv_dir = &jacob_vel_inv[0]; 
  return 0.0; 

} 
