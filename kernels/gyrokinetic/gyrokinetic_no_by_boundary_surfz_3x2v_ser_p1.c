#include <gkyl_gyrokinetic_kernels.h>
#include <gkyl_basis_gkhyb_3x2v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double gyrokinetic_no_by_boundary_surfz_3x2v_ser_p1(const double *w, const double *dxv,
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

  double rdz2 = 2.0/dxv[2];

  const double *GhatL = &flux_surf_skin[48];
  const double *GhatR = &flux_surf_edge[48];

  if (edge == -1) { 

  out[0] += -0.7071067811865475*GhatR[0]*rdz2; 
  out[1] += -0.7071067811865475*GhatR[1]*rdz2; 
  out[2] += -0.7071067811865475*GhatR[2]*rdz2; 
  out[3] += -1.224744871391589*GhatR[0]*rdz2; 
  out[4] += -0.7071067811865475*GhatR[3]*rdz2; 
  out[5] += -0.7071067811865475*GhatR[4]*rdz2; 
  out[6] += -0.7071067811865475*GhatR[5]*rdz2; 
  out[7] += -1.224744871391589*GhatR[1]*rdz2; 
  out[8] += -1.224744871391589*GhatR[2]*rdz2; 
  out[9] += -0.7071067811865475*GhatR[6]*rdz2; 
  out[10] += -0.7071067811865475*GhatR[7]*rdz2; 
  out[11] += -1.224744871391589*GhatR[3]*rdz2; 
  out[12] += -0.7071067811865475*GhatR[8]*rdz2; 
  out[13] += -0.7071067811865475*GhatR[9]*rdz2; 
  out[14] += -1.224744871391589*GhatR[4]*rdz2; 
  out[15] += -0.7071067811865475*GhatR[10]*rdz2; 
  out[16] += -1.224744871391589*GhatR[5]*rdz2; 
  out[17] += -0.7071067811865475*GhatR[11]*rdz2; 
  out[18] += -1.224744871391589*GhatR[6]*rdz2; 
  out[19] += -1.224744871391589*GhatR[7]*rdz2; 
  out[20] += -0.7071067811865475*GhatR[12]*rdz2; 
  out[21] += -1.224744871391589*GhatR[8]*rdz2; 
  out[22] += -1.224744871391589*GhatR[9]*rdz2; 
  out[23] += -0.7071067811865475*GhatR[13]*rdz2; 
  out[24] += -0.7071067811865475*GhatR[14]*rdz2; 
  out[25] += -1.224744871391589*GhatR[10]*rdz2; 
  out[26] += -1.224744871391589*GhatR[11]*rdz2; 
  out[27] += -1.224744871391589*GhatR[12]*rdz2; 
  out[28] += -0.7071067811865475*GhatR[15]*rdz2; 
  out[29] += -1.224744871391589*GhatR[13]*rdz2; 
  out[30] += -1.224744871391589*GhatR[14]*rdz2; 
  out[31] += -1.224744871391589*GhatR[15]*rdz2; 
  out[32] += -0.7071067811865475*GhatR[16]*rdz2; 
  out[33] += -0.7071067811865475*GhatR[17]*rdz2; 
  out[34] += -0.7071067811865475*GhatR[18]*rdz2; 
  out[35] += -1.224744871391589*GhatR[16]*rdz2; 
  out[36] += -0.7071067811865475*GhatR[19]*rdz2; 
  out[37] += -0.7071067811865475*GhatR[20]*rdz2; 
  out[38] += -1.224744871391589*GhatR[17]*rdz2; 
  out[39] += -1.224744871391589*GhatR[18]*rdz2; 
  out[40] += -0.7071067811865475*GhatR[21]*rdz2; 
  out[41] += -0.7071067811865475*GhatR[22]*rdz2; 
  out[42] += -1.224744871391589*GhatR[19]*rdz2; 
  out[43] += -1.224744871391589*GhatR[20]*rdz2; 
  out[44] += -0.7071067811865475*GhatR[23]*rdz2; 
  out[45] += -1.224744871391589*GhatR[21]*rdz2; 
  out[46] += -1.224744871391589*GhatR[22]*rdz2; 
  out[47] += -1.224744871391589*GhatR[23]*rdz2; 

  } else { 

  out[0] += 0.7071067811865475*GhatL[0]*rdz2; 
  out[1] += 0.7071067811865475*GhatL[1]*rdz2; 
  out[2] += 0.7071067811865475*GhatL[2]*rdz2; 
  out[3] += -1.224744871391589*GhatL[0]*rdz2; 
  out[4] += 0.7071067811865475*GhatL[3]*rdz2; 
  out[5] += 0.7071067811865475*GhatL[4]*rdz2; 
  out[6] += 0.7071067811865475*GhatL[5]*rdz2; 
  out[7] += -1.224744871391589*GhatL[1]*rdz2; 
  out[8] += -1.224744871391589*GhatL[2]*rdz2; 
  out[9] += 0.7071067811865475*GhatL[6]*rdz2; 
  out[10] += 0.7071067811865475*GhatL[7]*rdz2; 
  out[11] += -1.224744871391589*GhatL[3]*rdz2; 
  out[12] += 0.7071067811865475*GhatL[8]*rdz2; 
  out[13] += 0.7071067811865475*GhatL[9]*rdz2; 
  out[14] += -1.224744871391589*GhatL[4]*rdz2; 
  out[15] += 0.7071067811865475*GhatL[10]*rdz2; 
  out[16] += -1.224744871391589*GhatL[5]*rdz2; 
  out[17] += 0.7071067811865475*GhatL[11]*rdz2; 
  out[18] += -1.224744871391589*GhatL[6]*rdz2; 
  out[19] += -1.224744871391589*GhatL[7]*rdz2; 
  out[20] += 0.7071067811865475*GhatL[12]*rdz2; 
  out[21] += -1.224744871391589*GhatL[8]*rdz2; 
  out[22] += -1.224744871391589*GhatL[9]*rdz2; 
  out[23] += 0.7071067811865475*GhatL[13]*rdz2; 
  out[24] += 0.7071067811865475*GhatL[14]*rdz2; 
  out[25] += -1.224744871391589*GhatL[10]*rdz2; 
  out[26] += -1.224744871391589*GhatL[11]*rdz2; 
  out[27] += -1.224744871391589*GhatL[12]*rdz2; 
  out[28] += 0.7071067811865475*GhatL[15]*rdz2; 
  out[29] += -1.224744871391589*GhatL[13]*rdz2; 
  out[30] += -1.224744871391589*GhatL[14]*rdz2; 
  out[31] += -1.224744871391589*GhatL[15]*rdz2; 
  out[32] += 0.7071067811865475*GhatL[16]*rdz2; 
  out[33] += 0.7071067811865475*GhatL[17]*rdz2; 
  out[34] += 0.7071067811865475*GhatL[18]*rdz2; 
  out[35] += -1.224744871391589*GhatL[16]*rdz2; 
  out[36] += 0.7071067811865475*GhatL[19]*rdz2; 
  out[37] += 0.7071067811865475*GhatL[20]*rdz2; 
  out[38] += -1.224744871391589*GhatL[17]*rdz2; 
  out[39] += -1.224744871391589*GhatL[18]*rdz2; 
  out[40] += 0.7071067811865475*GhatL[21]*rdz2; 
  out[41] += 0.7071067811865475*GhatL[22]*rdz2; 
  out[42] += -1.224744871391589*GhatL[19]*rdz2; 
  out[43] += -1.224744871391589*GhatL[20]*rdz2; 
  out[44] += 0.7071067811865475*GhatL[23]*rdz2; 
  out[45] += -1.224744871391589*GhatL[21]*rdz2; 
  out[46] += -1.224744871391589*GhatL[22]*rdz2; 
  out[47] += -1.224744871391589*GhatL[23]*rdz2; 

  } 

  return 0.0; 

} 
