#include <gkyl_gyrokinetic_kernels.h>
#include <gkyl_basis_gkhyb_2x2v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double gyrokinetic_no_by_boundary_surfy_2x2v_ser_p1(const double *w, const double *dxv,
    const double *vmap_prime_edge, const double *vmap_prime_skin,
    const double *flux_surf_edge, const double *flux_surf_skin, 
    const int edge, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // vmap_prime_edge,vmap_prime_skin: velocity space mapping derivative in edge and skin cells.
  // flux_surf_edge: Surface expansion of phase space flux on the lower edges of the edge cell.
  // flux_surf_skin: Surface expansion of phase space flux on the lower edges of the skin cell.
  // edge: determines if the update is for the left edge (-1) or right edge (+1).
  // out: output increment in center cell.

  double rdz2 = 2.0/dxv[1];

  const double *GhatL = &alpha_surf_skin[12];
  const double *GhatR = &alpha_surf_edge[12];

  if (edge == -1) { 

  out[0] += -0.7071067811865475*GhatR[0]*rdz2; 
  out[1] += -0.7071067811865475*GhatR[1]*rdz2; 
  out[2] += -1.224744871391589*GhatR[0]*rdz2; 
  out[3] += -0.7071067811865475*GhatR[2]*rdz2; 
  out[4] += -0.7071067811865475*GhatR[3]*rdz2; 
  out[5] += -1.224744871391589*GhatR[1]*rdz2; 
  out[6] += -0.7071067811865475*GhatR[4]*rdz2; 
  out[7] += -1.224744871391589*GhatR[2]*rdz2; 
  out[8] += -0.7071067811865475*GhatR[5]*rdz2; 
  out[9] += -1.224744871391589*GhatR[3]*rdz2; 
  out[10] += -0.7071067811865475*GhatR[6]*rdz2; 
  out[11] += -1.224744871391589*GhatR[4]*rdz2; 
  out[12] += -1.224744871391589*GhatR[5]*rdz2; 
  out[13] += -0.7071067811865475*GhatR[7]*rdz2; 
  out[14] += -1.224744871391589*GhatR[6]*rdz2; 
  out[15] += -1.224744871391589*GhatR[7]*rdz2; 
  out[16] += -0.7071067811865475*GhatR[8]*rdz2; 
  out[17] += -0.7071067811865475*GhatR[9]*rdz2; 
  out[18] += -1.224744871391589*GhatR[8]*rdz2; 
  out[19] += -0.7071067811865475*GhatR[10]*rdz2; 
  out[20] += -1.224744871391589*GhatR[9]*rdz2; 
  out[21] += -0.7071067811865475*GhatR[11]*rdz2; 
  out[22] += -1.224744871391589*GhatR[10]*rdz2; 
  out[23] += -1.224744871391589*GhatR[11]*rdz2; 

  } else { 

  out[0] += 0.7071067811865475*GhatL[0]*rdz2; 
  out[1] += 0.7071067811865475*GhatL[1]*rdz2; 
  out[2] += -1.224744871391589*GhatL[0]*rdz2; 
  out[3] += 0.7071067811865475*GhatL[2]*rdz2; 
  out[4] += 0.7071067811865475*GhatL[3]*rdz2; 
  out[5] += -1.224744871391589*GhatL[1]*rdz2; 
  out[6] += 0.7071067811865475*GhatL[4]*rdz2; 
  out[7] += -1.224744871391589*GhatL[2]*rdz2; 
  out[8] += 0.7071067811865475*GhatL[5]*rdz2; 
  out[9] += -1.224744871391589*GhatL[3]*rdz2; 
  out[10] += 0.7071067811865475*GhatL[6]*rdz2; 
  out[11] += -1.224744871391589*GhatL[4]*rdz2; 
  out[12] += -1.224744871391589*GhatL[5]*rdz2; 
  out[13] += 0.7071067811865475*GhatL[7]*rdz2; 
  out[14] += -1.224744871391589*GhatL[6]*rdz2; 
  out[15] += -1.224744871391589*GhatL[7]*rdz2; 
  out[16] += 0.7071067811865475*GhatL[8]*rdz2; 
  out[17] += 0.7071067811865475*GhatL[9]*rdz2; 
  out[18] += -1.224744871391589*GhatL[8]*rdz2; 
  out[19] += 0.7071067811865475*GhatL[10]*rdz2; 
  out[20] += -1.224744871391589*GhatL[9]*rdz2; 
  out[21] += 0.7071067811865475*GhatL[11]*rdz2; 
  out[22] += -1.224744871391589*GhatL[10]*rdz2; 
  out[23] += -1.224744871391589*GhatL[11]*rdz2; 

  } 

  double cflFreq = fmax(fabs(alphaL[0]), fabs(alphaR[0])); 
  return 0.5303300858899105*rdz2*cflFreq; 

} 
