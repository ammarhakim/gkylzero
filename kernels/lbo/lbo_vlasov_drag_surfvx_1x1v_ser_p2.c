#include <gkyl_lbo_vlasov_kernels.h> 
#include <gkyl_basis_ser_2x_p2_surfx2_quad.h> 
#include <gkyl_basis_ser_2x_p2_upwind.h> 
GKYL_CU_DH void lbo_vlasov_drag_surfvx_1x1v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[2]:         cell-center coordinates. 
  // dxv[2]:       cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[3]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[3]: sum of thermal speeds squared time their respective collisionalities. 
  // fl/fc/fr:      distribution function in cells 
  // out:           incremented distribution function in cell 
  double rdv2 = 2.0/dxv[1]; 

  const double *sumNuUx = &nuUSum[0]; 

  double alphaDrSurf_l[3] = {0.0}; 
  alphaDrSurf_l[0] = nuSum[0]*w[1]-0.5*nuSum[0]*dxv[1]-1.0*sumNuUx[0]; 
  alphaDrSurf_l[1] = nuSum[1]*w[1]-1.0*sumNuUx[1]-0.5*dxv[1]*nuSum[1]; 
  alphaDrSurf_l[2] = (-1.0*sumNuUx[2])+w[1]*nuSum[2]-0.5*dxv[1]*nuSum[2]; 

  double alphaDrSurf_r[3] = {0.0}; 
  alphaDrSurf_r[0] = nuSum[0]*w[1]+0.5*nuSum[0]*dxv[1]-1.0*sumNuUx[0]; 
  alphaDrSurf_r[1] = nuSum[1]*w[1]-1.0*sumNuUx[1]+0.5*dxv[1]*nuSum[1]; 
  alphaDrSurf_r[2] = (-1.0*sumNuUx[2])+w[1]*nuSum[2]+0.5*dxv[1]*nuSum[2]; 

  double fUpwindQuad_l[3] = {0.0};
  double fUpwindQuad_r[3] = {0.0};
  double fUpwind_l[3] = {0.0};
  double fUpwind_r[3] = {0.0};
  double Gdrag_l[3] = {0.0}; 
  double Gdrag_r[3] = {0.0}; 

  if (0.6324555320336759*alphaDrSurf_l[2]-0.9486832980505137*alphaDrSurf_l[1]+0.7071067811865475*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[0] = ser_2x_p2_surfx2_quad_0_r(fl); 
  } else { 
    fUpwindQuad_l[0] = ser_2x_p2_surfx2_quad_0_l(fc); 
  } 
  if (0.6324555320336759*alphaDrSurf_r[2]-0.9486832980505137*alphaDrSurf_r[1]+0.7071067811865475*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[0] = ser_2x_p2_surfx2_quad_0_r(fc); 
  } else { 
    fUpwindQuad_r[0] = ser_2x_p2_surfx2_quad_0_l(fr); 
  } 
  if (0.7071067811865475*alphaDrSurf_l[0]-0.7905694150420947*alphaDrSurf_l[2] < 0) { 
    fUpwindQuad_l[1] = ser_2x_p2_surfx2_quad_1_r(fl); 
  } else { 
    fUpwindQuad_l[1] = ser_2x_p2_surfx2_quad_1_l(fc); 
  } 
  if (0.7071067811865475*alphaDrSurf_r[0]-0.7905694150420947*alphaDrSurf_r[2] < 0) { 
    fUpwindQuad_r[1] = ser_2x_p2_surfx2_quad_1_r(fc); 
  } else { 
    fUpwindQuad_r[1] = ser_2x_p2_surfx2_quad_1_l(fr); 
  } 
  if (0.6324555320336759*alphaDrSurf_l[2]+0.9486832980505137*alphaDrSurf_l[1]+0.7071067811865475*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[2] = ser_2x_p2_surfx2_quad_2_r(fl); 
  } else { 
    fUpwindQuad_l[2] = ser_2x_p2_surfx2_quad_2_l(fc); 
  } 
  if (0.6324555320336759*alphaDrSurf_r[2]+0.9486832980505137*alphaDrSurf_r[1]+0.7071067811865475*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[2] = ser_2x_p2_surfx2_quad_2_r(fc); 
  } else { 
    fUpwindQuad_r[2] = ser_2x_p2_surfx2_quad_2_l(fr); 
  } 

  // Project nodal basis back onto modal basis. 
  ser_2x_p2_upwind(fUpwindQuad_l, fUpwind_l); 
  ser_2x_p2_upwind(fUpwindQuad_r, fUpwind_r); 

  Gdrag_l[0] = 0.7071067811865475*alphaDrSurf_l[2]*fUpwind_l[2]+0.7071067811865475*alphaDrSurf_l[1]*fUpwind_l[1]+0.7071067811865475*alphaDrSurf_l[0]*fUpwind_l[0]; 
  Gdrag_l[1] = 0.6324555320336759*alphaDrSurf_l[1]*fUpwind_l[2]+0.6324555320336759*fUpwind_l[1]*alphaDrSurf_l[2]+0.7071067811865475*alphaDrSurf_l[0]*fUpwind_l[1]+0.7071067811865475*fUpwind_l[0]*alphaDrSurf_l[1]; 
  Gdrag_l[2] = 0.4517539514526256*alphaDrSurf_l[2]*fUpwind_l[2]+0.7071067811865475*alphaDrSurf_l[0]*fUpwind_l[2]+0.7071067811865475*fUpwind_l[0]*alphaDrSurf_l[2]+0.6324555320336759*alphaDrSurf_l[1]*fUpwind_l[1]; 

  Gdrag_r[0] = 0.7071067811865475*alphaDrSurf_r[2]*fUpwind_r[2]+0.7071067811865475*alphaDrSurf_r[1]*fUpwind_r[1]+0.7071067811865475*alphaDrSurf_r[0]*fUpwind_r[0]; 
  Gdrag_r[1] = 0.6324555320336759*alphaDrSurf_r[1]*fUpwind_r[2]+0.6324555320336759*fUpwind_r[1]*alphaDrSurf_r[2]+0.7071067811865475*alphaDrSurf_r[0]*fUpwind_r[1]+0.7071067811865475*fUpwind_r[0]*alphaDrSurf_r[1]; 
  Gdrag_r[2] = 0.4517539514526256*alphaDrSurf_r[2]*fUpwind_r[2]+0.7071067811865475*alphaDrSurf_r[0]*fUpwind_r[2]+0.7071067811865475*fUpwind_r[0]*alphaDrSurf_r[2]+0.6324555320336759*alphaDrSurf_r[1]*fUpwind_r[1]; 

  out[0] += 0.7071067811865475*Gdrag_r[0]*rdv2-0.7071067811865475*Gdrag_l[0]*rdv2; 
  out[1] += 0.7071067811865475*Gdrag_r[1]*rdv2-0.7071067811865475*Gdrag_l[1]*rdv2; 
  out[2] += 1.224744871391589*Gdrag_r[0]*rdv2+1.224744871391589*Gdrag_l[0]*rdv2; 
  out[3] += 1.224744871391589*Gdrag_r[1]*rdv2+1.224744871391589*Gdrag_l[1]*rdv2; 
  out[4] += 0.7071067811865475*Gdrag_r[2]*rdv2-0.7071067811865475*Gdrag_l[2]*rdv2; 
  out[5] += 1.58113883008419*Gdrag_r[0]*rdv2-1.58113883008419*Gdrag_l[0]*rdv2; 
  out[6] += 1.224744871391589*Gdrag_r[2]*rdv2+1.224744871391589*Gdrag_l[2]*rdv2; 
  out[7] += 1.58113883008419*Gdrag_r[1]*rdv2-1.58113883008419*Gdrag_l[1]*rdv2; 
} 
