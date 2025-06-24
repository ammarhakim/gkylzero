#include <gkyl_gyrokinetic_kernels.h>
#include <gkyl_basis_gkhyb_2x2v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double gyrokinetic_no_by_surfvpar_2x2v_ser_p1(const double *w, const double *dxv,
    const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r,
    const double *flux_surf_l, const double *flux_surf_r, 
    double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // vmap_prime_l,vmap_prime_c,vmap_prime_r: velocity space mapping derivative in left, center and right cells.
  // flux_surf_l: Surface expansion of phase space flux on the left.
  // flux_surf_r: Surface expansion of phase space flux on the right.
  // out: output increment in center cell.

  double rdvpar2 = 2.0/dxv[2];

  const double *GhatL = &flux_surf_l[24];
  const double *GhatR = &flux_surf_r[24];
  out[0] += (0.7071067811865475*GhatL[0]-0.7071067811865475*GhatR[0])*rdvpar2; 
  out[1] += (0.7071067811865475*GhatL[1]-0.7071067811865475*GhatR[1])*rdvpar2; 
  out[2] += (0.7071067811865475*GhatL[2]-0.7071067811865475*GhatR[2])*rdvpar2; 
  out[3] += ((-1.224744871391589*GhatR[0])-1.224744871391589*GhatL[0])*rdvpar2; 
  out[4] += (0.7071067811865475*GhatL[3]-0.7071067811865475*GhatR[3])*rdvpar2; 
  out[5] += (0.7071067811865475*GhatL[4]-0.7071067811865475*GhatR[4])*rdvpar2; 
  out[6] += ((-1.224744871391589*GhatR[1])-1.224744871391589*GhatL[1])*rdvpar2; 
  out[7] += ((-1.224744871391589*GhatR[2])-1.224744871391589*GhatL[2])*rdvpar2; 
  out[8] += (0.7071067811865475*GhatL[5]-0.7071067811865475*GhatR[5])*rdvpar2; 
  out[9] += (0.7071067811865475*GhatL[6]-0.7071067811865475*GhatR[6])*rdvpar2; 
  out[10] += ((-1.224744871391589*GhatR[3])-1.224744871391589*GhatL[3])*rdvpar2; 
  out[11] += ((-1.224744871391589*GhatR[4])-1.224744871391589*GhatL[4])*rdvpar2; 
  out[12] += (0.7071067811865475*GhatL[7]-0.7071067811865475*GhatR[7])*rdvpar2; 
  out[13] += ((-1.224744871391589*GhatR[5])-1.224744871391589*GhatL[5])*rdvpar2; 
  out[14] += ((-1.224744871391589*GhatR[6])-1.224744871391589*GhatL[6])*rdvpar2; 
  out[15] += ((-1.224744871391589*GhatR[7])-1.224744871391589*GhatL[7])*rdvpar2; 
  out[16] += (1.58113883008419*GhatL[0]-1.58113883008419*GhatR[0])*rdvpar2; 
  out[17] += (1.58113883008419*GhatL[1]-1.58113883008419*GhatR[1])*rdvpar2; 
  out[18] += (1.58113883008419*GhatL[2]-1.58113883008419*GhatR[2])*rdvpar2; 
  out[19] += (1.58113883008419*GhatL[3]-1.58113883008419*GhatR[3])*rdvpar2; 
  out[20] += (1.58113883008419*GhatL[4]-1.58113883008419*GhatR[4])*rdvpar2; 
  out[21] += (1.58113883008419*GhatL[5]-1.58113883008419*GhatR[5])*rdvpar2; 
  out[22] += (1.58113883008419*GhatL[6]-1.58113883008419*GhatR[6])*rdvpar2; 
  out[23] += (1.58113883008419*GhatL[7]-1.58113883008419*GhatR[7])*rdvpar2; 

  double vmap_prime_min = fmin(fmin(fabs(vmap_prime_l[0]),fabs(vmap_prime_c[0])),fabs(vmap_prime_r[0]));
  double cflFreq = fmax(fabs(flux_surf_l[0]/vmap_prime_min), fabs(flux_surf_r[0]/vmap_prime_min)); 
  return 0.8838834764831842*rdvpar2*cflFreq; 

} 
