#include <gkyl_gyrokinetic_kernels.h>
#include <gkyl_basis_gkhyb_1x2v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double gyrokinetic_boundary_surfx_1x2v_ser_p1(const double *w, const double *dxv,
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

  double rdx2 = 2.0/dxv[0];

  const double *GhatL = &alpha_surf_skin[0];
  const double *GhatR = &alpha_surf_edge[0];

  if (edge == -1) { 

  out[0] += -0.7071067811865475*GhatR[0]*rdx2; 
  out[1] += -1.224744871391589*GhatR[0]*rdx2; 
  out[2] += -0.7071067811865475*GhatR[1]*rdx2; 
  out[3] += -0.7071067811865475*GhatR[2]*rdx2; 
  out[4] += -1.224744871391589*GhatR[1]*rdx2; 
  out[5] += -1.224744871391589*GhatR[2]*rdx2; 
  out[6] += -0.7071067811865475*GhatR[3]*rdx2; 
  out[7] += -1.224744871391589*GhatR[3]*rdx2; 
  out[8] += -0.7071067811865475*GhatR[4]*rdx2; 
  out[9] += -1.224744871391589*GhatR[4]*rdx2; 
  out[10] += -0.7071067811865475*GhatR[5]*rdx2; 
  out[11] += -1.224744871391589*GhatR[5]*rdx2; 

  } else { 

  out[0] += 0.7071067811865475*GhatL[0]*rdx2; 
  out[1] += -1.224744871391589*GhatL[0]*rdx2; 
  out[2] += 0.7071067811865475*GhatL[1]*rdx2; 
  out[3] += 0.7071067811865475*GhatL[2]*rdx2; 
  out[4] += -1.224744871391589*GhatL[1]*rdx2; 
  out[5] += -1.224744871391589*GhatL[2]*rdx2; 
  out[6] += 0.7071067811865475*GhatL[3]*rdx2; 
  out[7] += -1.224744871391589*GhatL[3]*rdx2; 
  out[8] += 0.7071067811865475*GhatL[4]*rdx2; 
  out[9] += -1.224744871391589*GhatL[4]*rdx2; 
  out[10] += 0.7071067811865475*GhatL[5]*rdx2; 
  out[11] += -1.224744871391589*GhatL[5]*rdx2; 

  } 

  double cflFreq = fmax(fabs(alphaL[0]), fabs(alphaR[0])); 
  return 0.75*rdx2*cflFreq; 

} 
