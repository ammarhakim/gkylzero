#include <gkyl_vlasov_lbo_kernels.h> 
#include <gkyl_basis_ser_1x1v_p1_surfvx_quad.h> 
GKYL_CU_DH void vlasov_lbo_boundary_surfvx_1x1v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[2]:         Cell-center coordinates. 
  // dxv[2]:       Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[2]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[2]: sum of thermal speeds squared time their respective collisionalities. 
  // fSkin/Edge:    Distribution function in cells 
  // out:           Incremented distribution function in cell 
  double rdv2 = 2.0/dxv[1]; 
  double rdvSq4 = 4.0/(dxv[1]*dxv[1]); 

  const double *sumNuUx = &nuUSum[0]; 

  double alphaDrSurf[2] = {0.0}; 
  double fUpwindQuad[2] = {0.0};
  double fUpwind[2] = {0.0};;
  double Ghat[2] = {0.0}; 
  double Gdiff[2] = {0.0}; 
  double Gdiff2[2] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = 0.5*(nuSum[0]*(2.0*w[1]+dxv[1])-2.0*sumNuUx[0]); 
  alphaDrSurf[1] = 0.5*(2.0*nuSum[1]*w[1]-2.0*sumNuUx[1]+dxv[1]*nuSum[1]); 

  if (alphaDrSurf[0]-alphaDrSurf[1] < 0) { 
    fUpwindQuad[0] = ser_1x1v_p1_surfvx_quad_0(1, fSkin); 
  } else { 

    fUpwindQuad[0] = ser_1x1v_p1_surfvx_quad_0(-1, fEdge); 
  } 
  if (alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[1] = ser_1x1v_p1_surfvx_quad_1(1, fSkin); 
  } else { 

    fUpwindQuad[1] = ser_1x1v_p1_surfvx_quad_1(-1, fEdge); 
  } 

  fUpwind[0] = 0.7071067811865475*(fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[1] = 0.7071067811865475*(fUpwindQuad[1]-1.0*fUpwindQuad[0]); 

  Gdiff2[0] = 0.2886751345948129*nuVtSqSum[1]*fSkin[3]-0.2886751345948129*nuVtSqSum[1]*fEdge[3]+0.2886751345948129*nuVtSqSum[0]*fSkin[2]-0.2886751345948129*nuVtSqSum[0]*fEdge[2]+0.25*fSkin[1]*nuVtSqSum[1]+0.25*fEdge[1]*nuVtSqSum[1]+0.25*fSkin[0]*nuVtSqSum[0]+0.25*fEdge[0]*nuVtSqSum[0]; 
  Gdiff2[1] = 0.2886751345948129*nuVtSqSum[0]*fSkin[3]-0.2886751345948129*nuVtSqSum[0]*fEdge[3]+0.2886751345948129*nuVtSqSum[1]*fSkin[2]-0.2886751345948129*nuVtSqSum[1]*fEdge[2]+0.25*fSkin[0]*nuVtSqSum[1]+0.25*fEdge[0]*nuVtSqSum[1]+0.25*nuVtSqSum[0]*fSkin[1]+0.25*nuVtSqSum[0]*fEdge[1]; 

  Gdiff[0] = (-0.5412658773652741*nuVtSqSum[1]*fSkin[3])-0.5412658773652741*nuVtSqSum[1]*fEdge[3]-0.5412658773652741*nuVtSqSum[0]*fSkin[2]-0.5412658773652741*nuVtSqSum[0]*fEdge[2]-0.5625*fSkin[1]*nuVtSqSum[1]+0.5625*fEdge[1]*nuVtSqSum[1]-0.5625*fSkin[0]*nuVtSqSum[0]+0.5625*fEdge[0]*nuVtSqSum[0]; 
  Gdiff[1] = (-0.5412658773652741*nuVtSqSum[0]*fSkin[3])-0.5412658773652741*nuVtSqSum[0]*fEdge[3]-0.5412658773652741*nuVtSqSum[1]*fSkin[2]-0.5412658773652741*nuVtSqSum[1]*fEdge[2]-0.5625*fSkin[0]*nuVtSqSum[1]+0.5625*fEdge[0]*nuVtSqSum[1]-0.5625*nuVtSqSum[0]*fSkin[1]+0.5625*nuVtSqSum[0]*fEdge[1]; 

  Ghat[0] = (-1.0*Gdiff[0]*rdv2)+0.7071067811865475*alphaDrSurf[1]*fUpwind[1]+0.7071067811865475*alphaDrSurf[0]*fUpwind[0]; 
  Ghat[1] = (-1.0*Gdiff[1]*rdv2)+0.7071067811865475*alphaDrSurf[0]*fUpwind[1]+0.7071067811865475*fUpwind[0]*alphaDrSurf[1]; 

  out[0] += -0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += -0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += 1.060660171779821*nuVtSqSum[1]*fSkin[3]*rdvSq4+1.060660171779821*nuVtSqSum[0]*fSkin[2]*rdvSq4-0.6123724356957944*fSkin[1]*nuVtSqSum[1]*rdvSq4-0.6123724356957944*fSkin[0]*nuVtSqSum[0]*rdvSq4+1.224744871391589*Gdiff2[0]*rdvSq4-1.224744871391589*Ghat[0]*rdv2; 
  out[3] += 1.060660171779821*nuVtSqSum[0]*fSkin[3]*rdvSq4+1.060660171779821*nuVtSqSum[1]*fSkin[2]*rdvSq4-0.6123724356957944*fSkin[0]*nuVtSqSum[1]*rdvSq4-0.6123724356957944*nuVtSqSum[0]*fSkin[1]*rdvSq4+1.224744871391589*Gdiff2[1]*rdvSq4-1.224744871391589*Ghat[1]*rdv2; 

  } else { 

  alphaDrSurf[0] = 0.5*(nuSum[0]*(2.0*w[1]+dxv[1])-2.0*sumNuUx[0]); 
  alphaDrSurf[1] = 0.5*(2.0*nuSum[1]*w[1]-2.0*sumNuUx[1]+dxv[1]*nuSum[1]); 

  if (alphaDrSurf[0]-alphaDrSurf[1] < 0) { 
    fUpwindQuad[0] = ser_1x1v_p1_surfvx_quad_0(1, fEdge); 
  } else { 
    fUpwindQuad[0] = ser_1x1v_p1_surfvx_quad_0(-1, fSkin); 
  } 
  if (alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[1] = ser_1x1v_p1_surfvx_quad_1(1, fEdge); 
  } else { 
    fUpwindQuad[1] = ser_1x1v_p1_surfvx_quad_1(-1, fSkin); 
  } 

  fUpwind[0] = 0.7071067811865475*(fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[1] = 0.7071067811865475*(fUpwindQuad[1]-1.0*fUpwindQuad[0]); 

  Gdiff2[0] = (-0.2886751345948129*nuVtSqSum[1]*fSkin[3])+0.2886751345948129*nuVtSqSum[1]*fEdge[3]-0.2886751345948129*nuVtSqSum[0]*fSkin[2]+0.2886751345948129*nuVtSqSum[0]*fEdge[2]+0.25*fSkin[1]*nuVtSqSum[1]+0.25*fEdge[1]*nuVtSqSum[1]+0.25*fSkin[0]*nuVtSqSum[0]+0.25*fEdge[0]*nuVtSqSum[0]; 
  Gdiff2[1] = (-0.2886751345948129*nuVtSqSum[0]*fSkin[3])+0.2886751345948129*nuVtSqSum[0]*fEdge[3]-0.2886751345948129*nuVtSqSum[1]*fSkin[2]+0.2886751345948129*nuVtSqSum[1]*fEdge[2]+0.25*fSkin[0]*nuVtSqSum[1]+0.25*fEdge[0]*nuVtSqSum[1]+0.25*nuVtSqSum[0]*fSkin[1]+0.25*nuVtSqSum[0]*fEdge[1]; 

  Gdiff[0] = (-0.5412658773652741*nuVtSqSum[1]*fSkin[3])-0.5412658773652741*nuVtSqSum[1]*fEdge[3]-0.5412658773652741*nuVtSqSum[0]*fSkin[2]-0.5412658773652741*nuVtSqSum[0]*fEdge[2]+0.5625*fSkin[1]*nuVtSqSum[1]-0.5625*fEdge[1]*nuVtSqSum[1]+0.5625*fSkin[0]*nuVtSqSum[0]-0.5625*fEdge[0]*nuVtSqSum[0]; 
  Gdiff[1] = (-0.5412658773652741*nuVtSqSum[0]*fSkin[3])-0.5412658773652741*nuVtSqSum[0]*fEdge[3]-0.5412658773652741*nuVtSqSum[1]*fSkin[2]-0.5412658773652741*nuVtSqSum[1]*fEdge[2]+0.5625*fSkin[0]*nuVtSqSum[1]-0.5625*fEdge[0]*nuVtSqSum[1]+0.5625*nuVtSqSum[0]*fSkin[1]-0.5625*nuVtSqSum[0]*fEdge[1]; 

  Ghat[0] = Gdiff[0]*rdv2+0.7071067811865475*alphaDrSurf[1]*fUpwind[1]+0.7071067811865475*alphaDrSurf[0]*fUpwind[0]; 
  Ghat[1] = Gdiff[1]*rdv2+0.7071067811865475*alphaDrSurf[0]*fUpwind[1]+0.7071067811865475*fUpwind[0]*alphaDrSurf[1]; 

  out[0] += 0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += 1.060660171779821*nuVtSqSum[1]*fSkin[3]*rdvSq4+1.060660171779821*nuVtSqSum[0]*fSkin[2]*rdvSq4+0.6123724356957944*fSkin[1]*nuVtSqSum[1]*rdvSq4+0.6123724356957944*fSkin[0]*nuVtSqSum[0]*rdvSq4-1.224744871391589*Gdiff2[0]*rdvSq4-1.224744871391589*Ghat[0]*rdv2; 
  out[3] += 1.060660171779821*nuVtSqSum[0]*fSkin[3]*rdvSq4+1.060660171779821*nuVtSqSum[1]*fSkin[2]*rdvSq4+0.6123724356957944*fSkin[0]*nuVtSqSum[1]*rdvSq4+0.6123724356957944*nuVtSqSum[0]*fSkin[1]*rdvSq4-1.224744871391589*Gdiff2[1]*rdvSq4-1.224744871391589*Ghat[1]*rdv2; 

  } 
} 
