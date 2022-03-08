#include <gkyl_lbo_vlasov_kernels.h> 
#include <gkyl_basis_ser_2x_p1_surfx2_quad.h> 
#include <gkyl_basis_ser_2x_p1_upwind.h> 
GKYL_CU_DH void lbo_vlasov_drag_surfvx_1x1v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[2]:         cell-center coordinates. 
  // dxv[2]:       cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[2]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[2]: sum of thermal speeds squared time their respective collisionalities. 
  // fl/fc/fr:      distribution function in cells 
  // out:           incremented distribution function in cell 
  double rdv2 = 2.0/dxv[1]; 

  const double *sumNuUx = &nuUSum[0]; 

  double alphaDrSurf_l[2] = {0.0}; 
  alphaDrSurf_l[0] = nuSum[0]*w[1]-0.5*nuSum[0]*dxv[1]-1.0*sumNuUx[0]; 
  alphaDrSurf_l[1] = nuSum[1]*w[1]-1.0*sumNuUx[1]-0.5*dxv[1]*nuSum[1]; 

  double alphaDrSurf_r[2] = {0.0}; 
  alphaDrSurf_r[0] = nuSum[0]*w[1]+0.5*nuSum[0]*dxv[1]-1.0*sumNuUx[0]; 
  alphaDrSurf_r[1] = nuSum[1]*w[1]-1.0*sumNuUx[1]+0.5*dxv[1]*nuSum[1]; 

  double fUpwindQuad_l[2] = {0.0};
  double fUpwindQuad_r[2] = {0.0};
  double fUpwind_l[2] = {0.0};
  double fUpwind_r[2] = {0.0};
  double Gdrag_l[2] = {0.0}; 
  double Gdrag_r[2] = {0.0}; 

  if (alphaDrSurf_l[0]-alphaDrSurf_l[1] < 0) { 
    fUpwindQuad_l[0] = ser_2x_p1_surfx2_quad_0_r(fl); 
  } else { 
    fUpwindQuad_l[0] = ser_2x_p1_surfx2_quad_0_l(fc); 
  } 
  if (alphaDrSurf_r[0]-alphaDrSurf_r[1] < 0) { 
    fUpwindQuad_r[0] = ser_2x_p1_surfx2_quad_0_r(fc); 
  } else { 
    fUpwindQuad_r[0] = ser_2x_p1_surfx2_quad_0_l(fr); 
  } 
  if (alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[1] = ser_2x_p1_surfx2_quad_1_r(fl); 
  } else { 
    fUpwindQuad_l[1] = ser_2x_p1_surfx2_quad_1_l(fc); 
  } 
  if (alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[1] = ser_2x_p1_surfx2_quad_1_r(fc); 
  } else { 
    fUpwindQuad_r[1] = ser_2x_p1_surfx2_quad_1_l(fr); 
  } 

  // Project nodal basis back onto modal basis. 
  ser_2x_p1_upwind(fUpwindQuad_l, fUpwind_l); 
  ser_2x_p1_upwind(fUpwindQuad_r, fUpwind_r); 

  Gdrag_l[0] = 0.7071067811865475*alphaDrSurf_l[1]*fUpwind_l[1]+0.7071067811865475*alphaDrSurf_l[0]*fUpwind_l[0]; 
  Gdrag_l[1] = 0.7071067811865475*alphaDrSurf_l[0]*fUpwind_l[1]+0.7071067811865475*fUpwind_l[0]*alphaDrSurf_l[1]; 

  Gdrag_r[0] = 0.7071067811865475*alphaDrSurf_r[1]*fUpwind_r[1]+0.7071067811865475*alphaDrSurf_r[0]*fUpwind_r[0]; 
  Gdrag_r[1] = 0.7071067811865475*alphaDrSurf_r[0]*fUpwind_r[1]+0.7071067811865475*fUpwind_r[0]*alphaDrSurf_r[1]; 

  out[0] += 0.7071067811865475*Gdrag_r[0]*rdv2-0.7071067811865475*Gdrag_l[0]*rdv2; 
  out[1] += 0.7071067811865475*Gdrag_r[1]*rdv2-0.7071067811865475*Gdrag_l[1]*rdv2; 
  out[2] += 1.224744871391589*Gdrag_r[0]*rdv2+1.224744871391589*Gdrag_l[0]*rdv2; 
  out[3] += 1.224744871391589*Gdrag_r[1]*rdv2+1.224744871391589*Gdrag_l[1]*rdv2; 
} 
