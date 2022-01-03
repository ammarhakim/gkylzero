#include <gkyl_vlasov_lbo_kernels.h> 
#include <gkyl_basis_ser_1x1v_p2_surfvx_quad.h> 
GKYL_CU_DH void vlasov_lbo_drag_boundary_surfvx_1x1v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[2]:         Cell-center coordinates. 
  // dxv[2]:       Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[3]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[3]: sum of thermal speeds squared time their respective collisionalities. 
  // fSkin/Edge:    Distribution function in cells 
  // out:           Incremented distribution function in cell 
  double rdv2 = 2.0/dxv[1]; 

  const double *sumNuUx = &nuUSum[0]; 

  double alphaDrSurf[3] = {0.0}; 
  double fUpwindQuad[3] = {0.0};
  double fUpwind[3] = {0.0};;
  double drag_incr[3] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = 0.5*(nuSum[0]*(2.0*w[1]+dxv[1])-2.0*sumNuUx[0]); 
  alphaDrSurf[1] = 0.5*(2.0*nuSum[1]*w[1]-2.0*sumNuUx[1]+dxv[1]*nuSum[1]); 
  alphaDrSurf[2] = -0.5*(2.0*sumNuUx[2]+((-2.0*w[1])-1.0*dxv[1])*nuSum[2]); 

  if (0.6324555320336759*alphaDrSurf[2]-0.9486832980505137*alphaDrSurf[1]+0.7071067811865475*alphaDrSurf[0] < 0) { 
    fUpwindQuad[0] = ser_1x1v_p2_surfvx_quad_0(1, fSkin); 
  } else { 
    fUpwindQuad[0] = ser_1x1v_p2_surfvx_quad_0(-1, fEdge); 
  } 
  if (0.7071067811865475*alphaDrSurf[0]-0.7905694150420947*alphaDrSurf[2] < 0) { 
    fUpwindQuad[1] = ser_1x1v_p2_surfvx_quad_1(1, fSkin); 
  } else { 
    fUpwindQuad[1] = ser_1x1v_p2_surfvx_quad_1(-1, fEdge); 
  } 
  if (0.6324555320336759*alphaDrSurf[2]+0.9486832980505137*alphaDrSurf[1]+0.7071067811865475*alphaDrSurf[0] < 0) { 
    fUpwindQuad[2] = ser_1x1v_p2_surfvx_quad_2(1, fSkin); 
  } else { 
    fUpwindQuad[2] = ser_1x1v_p2_surfvx_quad_2(-1, fEdge); 
  } 

  fUpwind[0] = 0.392837100659193*fUpwindQuad[2]+0.6285393610547091*fUpwindQuad[1]+0.392837100659193*fUpwindQuad[0]; 
  fUpwind[1] = 0.5270462766947298*fUpwindQuad[2]-0.5270462766947298*fUpwindQuad[0]; 
  fUpwind[2] = 0.3513641844631533*fUpwindQuad[2]-0.7027283689263066*fUpwindQuad[1]+0.3513641844631533*fUpwindQuad[0]; 

  drag_incr[0] = 0.7071067811865475*alphaDrSurf[2]*fUpwind[2]+0.7071067811865475*alphaDrSurf[1]*fUpwind[1]+0.7071067811865475*alphaDrSurf[0]*fUpwind[0]; 
  drag_incr[1] = 0.6324555320336759*alphaDrSurf[1]*fUpwind[2]+0.6324555320336759*fUpwind[1]*alphaDrSurf[2]+0.7071067811865475*alphaDrSurf[0]*fUpwind[1]+0.7071067811865475*fUpwind[0]*alphaDrSurf[1]; 
  drag_incr[2] = 0.4517539514526256*alphaDrSurf[2]*fUpwind[2]+0.7071067811865475*alphaDrSurf[0]*fUpwind[2]+0.7071067811865475*fUpwind[0]*alphaDrSurf[2]+0.6324555320336759*alphaDrSurf[1]*fUpwind[1]; 

  out[0] += 0.7071067811865475*drag_incr[0]*rdv2; 
  out[1] += 0.7071067811865475*drag_incr[1]*rdv2; 
  out[2] += 1.224744871391589*drag_incr[0]*rdv2; 
  out[3] += 1.224744871391589*drag_incr[1]*rdv2; 
  out[4] += 0.7071067811865475*drag_incr[2]*rdv2; 
  out[5] += 1.58113883008419*drag_incr[0]*rdv2; 
  out[6] += 1.224744871391589*drag_incr[2]*rdv2; 
  out[7] += 1.58113883008419*drag_incr[1]*rdv2; 

  } else { 

  alphaDrSurf[0] = 0.5*(nuSum[0]*(2.0*w[1]-1.0*dxv[1])-2.0*sumNuUx[0]); 
  alphaDrSurf[1] = 0.5*(2.0*nuSum[1]*w[1]-2.0*sumNuUx[1]-1.0*dxv[1]*nuSum[1]); 
  alphaDrSurf[2] = -0.5*(2.0*sumNuUx[2]+(dxv[1]-2.0*w[1])*nuSum[2]); 

  if (0.6324555320336759*alphaDrSurf[2]-0.9486832980505137*alphaDrSurf[1]+0.7071067811865475*alphaDrSurf[0] < 0) { 
    fUpwindQuad[0] = ser_1x1v_p2_surfvx_quad_0(1, fEdge); 
  } else { 
    fUpwindQuad[0] = ser_1x1v_p2_surfvx_quad_0(-1, fSkin); 
  } 
  if (0.7071067811865475*alphaDrSurf[0]-0.7905694150420947*alphaDrSurf[2] < 0) { 
    fUpwindQuad[1] = ser_1x1v_p2_surfvx_quad_1(1, fEdge); 
  } else { 
    fUpwindQuad[1] = ser_1x1v_p2_surfvx_quad_1(-1, fSkin); 
  } 
  if (0.6324555320336759*alphaDrSurf[2]+0.9486832980505137*alphaDrSurf[1]+0.7071067811865475*alphaDrSurf[0] < 0) { 
    fUpwindQuad[2] = ser_1x1v_p2_surfvx_quad_2(1, fEdge); 
  } else { 
    fUpwindQuad[2] = ser_1x1v_p2_surfvx_quad_2(-1, fSkin); 
  } 

  fUpwind[0] = 0.392837100659193*fUpwindQuad[2]+0.6285393610547091*fUpwindQuad[1]+0.392837100659193*fUpwindQuad[0]; 
  fUpwind[1] = 0.5270462766947298*fUpwindQuad[2]-0.5270462766947298*fUpwindQuad[0]; 
  fUpwind[2] = 0.3513641844631533*fUpwindQuad[2]-0.7027283689263066*fUpwindQuad[1]+0.3513641844631533*fUpwindQuad[0]; 

  drag_incr[0] = 0.7071067811865475*alphaDrSurf[2]*fUpwind[2]+0.7071067811865475*alphaDrSurf[1]*fUpwind[1]+0.7071067811865475*alphaDrSurf[0]*fUpwind[0]; 
  drag_incr[1] = 0.6324555320336759*alphaDrSurf[1]*fUpwind[2]+0.6324555320336759*fUpwind[1]*alphaDrSurf[2]+0.7071067811865475*alphaDrSurf[0]*fUpwind[1]+0.7071067811865475*fUpwind[0]*alphaDrSurf[1]; 
  drag_incr[2] = 0.4517539514526256*alphaDrSurf[2]*fUpwind[2]+0.7071067811865475*alphaDrSurf[0]*fUpwind[2]+0.7071067811865475*fUpwind[0]*alphaDrSurf[2]+0.6324555320336759*alphaDrSurf[1]*fUpwind[1]; 

  out[0] += -0.7071067811865475*drag_incr[0]*rdv2; 
  out[1] += -0.7071067811865475*drag_incr[1]*rdv2; 
  out[2] += 1.224744871391589*drag_incr[0]*rdv2; 
  out[3] += 1.224744871391589*drag_incr[1]*rdv2; 
  out[4] += -0.7071067811865475*drag_incr[2]*rdv2; 
  out[5] += -1.58113883008419*drag_incr[0]*rdv2; 
  out[6] += 1.224744871391589*drag_incr[2]*rdv2; 
  out[7] += -1.58113883008419*drag_incr[1]*rdv2; 

  } 
} 
