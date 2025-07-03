#include <gkyl_vlasov_sr_kernels.h> 
GKYL_CU_DH double vlasov_sr_vmap_surfy_2x1v_ser_p2(const double *w, const double *dxv, const double *jacob_vel_inv, const double *gamma, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // jacob_vel_inv[VDIM]: Inverse velocity space Jacobian in each direction (unused in uniform grid simulations).
  // gamma:     Particle Lorentz boost factor sqrt(1 + p^2).
  // fl/fc/fr:  Input Distribution function in left/center/right cells.
  // out:       Incremented distribution function in center cell.
  return 0.;

} 
