#include <gkyl_gyrokinetic_kernels.h>
#include <gkyl_basis_gkhyb_3x2v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double gyrokinetic_boundary_surfvpar_3x2v_ser_p1(const double *w, const double *dxv,
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

  double rdvpar2 = 2.0/dxv[3];

  const double *GhatL = &flux_surf_skin[72];
  const double *GhatR = &flux_surf_edge[72];

  if (edge == -1) { 

  out[0] += -0.7071067811865475*GhatR[0]*rdvpar2; 
  out[1] += -0.7071067811865475*GhatR[1]*rdvpar2; 
  out[2] += -0.7071067811865475*GhatR[2]*rdvpar2; 
  out[3] += -0.7071067811865475*GhatR[3]*rdvpar2; 
  out[4] += -1.224744871391589*GhatR[0]*rdvpar2; 
  out[5] += -0.7071067811865475*GhatR[4]*rdvpar2; 
  out[6] += -0.7071067811865475*GhatR[5]*rdvpar2; 
  out[7] += -0.7071067811865475*GhatR[6]*rdvpar2; 
  out[8] += -0.7071067811865475*GhatR[7]*rdvpar2; 
  out[9] += -1.224744871391589*GhatR[1]*rdvpar2; 
  out[10] += -1.224744871391589*GhatR[2]*rdvpar2; 
  out[11] += -1.224744871391589*GhatR[3]*rdvpar2; 
  out[12] += -0.7071067811865475*GhatR[8]*rdvpar2; 
  out[13] += -0.7071067811865475*GhatR[9]*rdvpar2; 
  out[14] += -0.7071067811865475*GhatR[10]*rdvpar2; 
  out[15] += -1.224744871391589*GhatR[4]*rdvpar2; 
  out[16] += -0.7071067811865475*GhatR[11]*rdvpar2; 
  out[17] += -1.224744871391589*GhatR[5]*rdvpar2; 
  out[18] += -1.224744871391589*GhatR[6]*rdvpar2; 
  out[19] += -1.224744871391589*GhatR[7]*rdvpar2; 
  out[20] += -0.7071067811865475*GhatR[12]*rdvpar2; 
  out[21] += -0.7071067811865475*GhatR[13]*rdvpar2; 
  out[22] += -0.7071067811865475*GhatR[14]*rdvpar2; 
  out[23] += -1.224744871391589*GhatR[8]*rdvpar2; 
  out[24] += -1.224744871391589*GhatR[9]*rdvpar2; 
  out[25] += -1.224744871391589*GhatR[10]*rdvpar2; 
  out[26] += -1.224744871391589*GhatR[11]*rdvpar2; 
  out[27] += -0.7071067811865475*GhatR[15]*rdvpar2; 
  out[28] += -1.224744871391589*GhatR[12]*rdvpar2; 
  out[29] += -1.224744871391589*GhatR[13]*rdvpar2; 
  out[30] += -1.224744871391589*GhatR[14]*rdvpar2; 
  out[31] += -1.224744871391589*GhatR[15]*rdvpar2; 
  out[32] += -1.58113883008419*GhatR[0]*rdvpar2; 
  out[33] += -1.58113883008419*GhatR[1]*rdvpar2; 
  out[34] += -1.58113883008419*GhatR[2]*rdvpar2; 
  out[35] += -1.58113883008419*GhatR[3]*rdvpar2; 
  out[36] += -1.58113883008419*GhatR[4]*rdvpar2; 
  out[37] += -1.58113883008419*GhatR[5]*rdvpar2; 
  out[38] += -1.58113883008419*GhatR[6]*rdvpar2; 
  out[39] += -1.58113883008419*GhatR[7]*rdvpar2; 
  out[40] += -1.58113883008419*GhatR[8]*rdvpar2; 
  out[41] += -1.58113883008419*GhatR[9]*rdvpar2; 
  out[42] += -1.58113883008419*GhatR[10]*rdvpar2; 
  out[43] += -1.58113883008419*GhatR[11]*rdvpar2; 
  out[44] += -1.58113883008419*GhatR[12]*rdvpar2; 
  out[45] += -1.58113883008419*GhatR[13]*rdvpar2; 
  out[46] += -1.58113883008419*GhatR[14]*rdvpar2; 
  out[47] += -1.58113883008419*GhatR[15]*rdvpar2; 

  } else { 

  out[0] += 0.7071067811865475*GhatL[0]*rdvpar2; 
  out[1] += 0.7071067811865475*GhatL[1]*rdvpar2; 
  out[2] += 0.7071067811865475*GhatL[2]*rdvpar2; 
  out[3] += 0.7071067811865475*GhatL[3]*rdvpar2; 
  out[4] += -1.224744871391589*GhatL[0]*rdvpar2; 
  out[5] += 0.7071067811865475*GhatL[4]*rdvpar2; 
  out[6] += 0.7071067811865475*GhatL[5]*rdvpar2; 
  out[7] += 0.7071067811865475*GhatL[6]*rdvpar2; 
  out[8] += 0.7071067811865475*GhatL[7]*rdvpar2; 
  out[9] += -1.224744871391589*GhatL[1]*rdvpar2; 
  out[10] += -1.224744871391589*GhatL[2]*rdvpar2; 
  out[11] += -1.224744871391589*GhatL[3]*rdvpar2; 
  out[12] += 0.7071067811865475*GhatL[8]*rdvpar2; 
  out[13] += 0.7071067811865475*GhatL[9]*rdvpar2; 
  out[14] += 0.7071067811865475*GhatL[10]*rdvpar2; 
  out[15] += -1.224744871391589*GhatL[4]*rdvpar2; 
  out[16] += 0.7071067811865475*GhatL[11]*rdvpar2; 
  out[17] += -1.224744871391589*GhatL[5]*rdvpar2; 
  out[18] += -1.224744871391589*GhatL[6]*rdvpar2; 
  out[19] += -1.224744871391589*GhatL[7]*rdvpar2; 
  out[20] += 0.7071067811865475*GhatL[12]*rdvpar2; 
  out[21] += 0.7071067811865475*GhatL[13]*rdvpar2; 
  out[22] += 0.7071067811865475*GhatL[14]*rdvpar2; 
  out[23] += -1.224744871391589*GhatL[8]*rdvpar2; 
  out[24] += -1.224744871391589*GhatL[9]*rdvpar2; 
  out[25] += -1.224744871391589*GhatL[10]*rdvpar2; 
  out[26] += -1.224744871391589*GhatL[11]*rdvpar2; 
  out[27] += 0.7071067811865475*GhatL[15]*rdvpar2; 
  out[28] += -1.224744871391589*GhatL[12]*rdvpar2; 
  out[29] += -1.224744871391589*GhatL[13]*rdvpar2; 
  out[30] += -1.224744871391589*GhatL[14]*rdvpar2; 
  out[31] += -1.224744871391589*GhatL[15]*rdvpar2; 
  out[32] += 1.58113883008419*GhatL[0]*rdvpar2; 
  out[33] += 1.58113883008419*GhatL[1]*rdvpar2; 
  out[34] += 1.58113883008419*GhatL[2]*rdvpar2; 
  out[35] += 1.58113883008419*GhatL[3]*rdvpar2; 
  out[36] += 1.58113883008419*GhatL[4]*rdvpar2; 
  out[37] += 1.58113883008419*GhatL[5]*rdvpar2; 
  out[38] += 1.58113883008419*GhatL[6]*rdvpar2; 
  out[39] += 1.58113883008419*GhatL[7]*rdvpar2; 
  out[40] += 1.58113883008419*GhatL[8]*rdvpar2; 
  out[41] += 1.58113883008419*GhatL[9]*rdvpar2; 
  out[42] += 1.58113883008419*GhatL[10]*rdvpar2; 
  out[43] += 1.58113883008419*GhatL[11]*rdvpar2; 
  out[44] += 1.58113883008419*GhatL[12]*rdvpar2; 
  out[45] += 1.58113883008419*GhatL[13]*rdvpar2; 
  out[46] += 1.58113883008419*GhatL[14]*rdvpar2; 
  out[47] += 1.58113883008419*GhatL[15]*rdvpar2; 

  } 

  double vmap_prime_min = fmin(fabs(vmap_prime_edge[0]),fabs(vmap_prime_skin[0]));
  double cflFreq = fmax(fabs(flux_surf_edge[0]/vmap_prime_min), fabs(flux_surf_skin[0]/vmap_prime_min)); 
  return 0.625*rdvpar2*cflFreq; 

} 
