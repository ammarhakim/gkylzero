#include <gkyl_gyrokinetic_kernels.h>
#include <gkyl_basis_gkhyb_2x2v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double gyrokinetic_surfx_2x2v_ser_p1(const double *w, const double *dxv,
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

  double rdx2 = 2.0/dxv[0];

  const double *GhatL = &flux_surf_l[0];
  const double *GhatR = &flux_surf_r[0];
  out[0] += (0.7071067811865475*GhatL[0]-0.7071067811865475*GhatR[0])*rdx2; 
  out[1] += ((-1.224744871391589*GhatR[0])-1.224744871391589*GhatL[0])*rdx2; 
  out[2] += (0.7071067811865475*GhatL[1]-0.7071067811865475*GhatR[1])*rdx2; 
  out[3] += (0.7071067811865475*GhatL[2]-0.7071067811865475*GhatR[2])*rdx2; 
  out[4] += (0.7071067811865475*GhatL[3]-0.7071067811865475*GhatR[3])*rdx2; 
  out[5] += ((-1.224744871391589*GhatR[1])-1.224744871391589*GhatL[1])*rdx2; 
  out[6] += ((-1.224744871391589*GhatR[2])-1.224744871391589*GhatL[2])*rdx2; 
  out[7] += (0.7071067811865475*GhatL[4]-0.7071067811865475*GhatR[4])*rdx2; 
  out[8] += ((-1.224744871391589*GhatR[3])-1.224744871391589*GhatL[3])*rdx2; 
  out[9] += (0.7071067811865475*GhatL[5]-0.7071067811865475*GhatR[5])*rdx2; 
  out[10] += (0.7071067811865475*GhatL[6]-0.7071067811865475*GhatR[6])*rdx2; 
  out[11] += ((-1.224744871391589*GhatR[4])-1.224744871391589*GhatL[4])*rdx2; 
  out[12] += ((-1.224744871391589*GhatR[5])-1.224744871391589*GhatL[5])*rdx2; 
  out[13] += ((-1.224744871391589*GhatR[6])-1.224744871391589*GhatL[6])*rdx2; 
  out[14] += (0.7071067811865475*GhatL[7]-0.7071067811865475*GhatR[7])*rdx2; 
  out[15] += ((-1.224744871391589*GhatR[7])-1.224744871391589*GhatL[7])*rdx2; 
  out[16] += (0.7071067811865475*GhatL[8]-0.7071067811865475*GhatR[8])*rdx2; 
  out[17] += ((-1.224744871391589*GhatR[8])-1.224744871391589*GhatL[8])*rdx2; 
  out[18] += (0.7071067811865475*GhatL[9]-0.7071067811865475*GhatR[9])*rdx2; 
  out[19] += (0.7071067811865475*GhatL[10]-0.7071067811865475*GhatR[10])*rdx2; 
  out[20] += ((-1.224744871391589*GhatR[9])-1.224744871391589*GhatL[9])*rdx2; 
  out[21] += ((-1.224744871391589*GhatR[10])-1.224744871391589*GhatL[10])*rdx2; 
  out[22] += (0.7071067811865475*GhatL[11]-0.7071067811865475*GhatR[11])*rdx2; 
  out[23] += ((-1.224744871391589*GhatR[11])-1.224744871391589*GhatL[11])*rdx2; 

  double cflFreq = fmax(fabs(flux_surf_l[0]), fabs(flux_surf_r[0])); 
  return 0.5303300858899105*rdx2*cflFreq; 

} 
