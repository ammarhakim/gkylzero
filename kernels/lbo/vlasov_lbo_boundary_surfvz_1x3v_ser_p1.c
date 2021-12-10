#include <gkyl_vlasov_lbo_kernels.h> 
#include <gkyl_basis_ser_1x3v_p1_surfvz_quad.h> 
GKYL_CU_DH void vlasov_lbo_boundary_surfvz_1x3v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[4]:         Cell-center coordinates. 
  // dxv[4]:       Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[6]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[2]: sum of thermal speeds squared time their respective collisionalities. 
  // fSkin/Edge:    Distribution function in cells 
  // out:           Incremented distribution function in cell 
  double rdv2 = 2.0/dxv[3]; 

  const double *sumNuUz = &nuUSum[4]; 

  double alphaDrSurf[8] = {0.0}; 
  double fUpwindQuad[8] = {0.0};
  double fUpwind[8] = {0.0};;
  double drag_incr[8] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = nuSum[0]*(2.0*w[3]+dxv[3])-2.0*sumNuUz[0]; 
  alphaDrSurf[1] = nuSum[1]*(2.0*w[3]+dxv[3])-2.0*sumNuUz[1]; 

  if (alphaDrSurf[0]-alphaDrSurf[1] < 0) { 
    fUpwindQuad[0] = ser_1x3v_p1_surfvz_quad_0(1, fSkin); 
  } else { 

    fUpwindQuad[0] = ser_1x3v_p1_surfvz_quad_0(-1, fEdge); 
  } 
  if (alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[1] = ser_1x3v_p1_surfvz_quad_1(1, fSkin); 
  } else { 

    fUpwindQuad[1] = ser_1x3v_p1_surfvz_quad_1(-1, fEdge); 
  } 
  if (alphaDrSurf[0]-alphaDrSurf[1] < 0) { 
    fUpwindQuad[2] = ser_1x3v_p1_surfvz_quad_2(1, fSkin); 
  } else { 

    fUpwindQuad[2] = ser_1x3v_p1_surfvz_quad_2(-1, fEdge); 
  } 
  if (alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[3] = ser_1x3v_p1_surfvz_quad_3(1, fSkin); 
  } else { 

    fUpwindQuad[3] = ser_1x3v_p1_surfvz_quad_3(-1, fEdge); 
  } 
  if (alphaDrSurf[0]-alphaDrSurf[1] < 0) { 
    fUpwindQuad[4] = ser_1x3v_p1_surfvz_quad_4(1, fSkin); 
  } else { 

    fUpwindQuad[4] = ser_1x3v_p1_surfvz_quad_4(-1, fEdge); 
  } 
  if (alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[5] = ser_1x3v_p1_surfvz_quad_5(1, fSkin); 
  } else { 

    fUpwindQuad[5] = ser_1x3v_p1_surfvz_quad_5(-1, fEdge); 
  } 
  if (alphaDrSurf[0]-alphaDrSurf[1] < 0) { 
    fUpwindQuad[6] = ser_1x3v_p1_surfvz_quad_6(1, fSkin); 
  } else { 

    fUpwindQuad[6] = ser_1x3v_p1_surfvz_quad_6(-1, fEdge); 
  } 
  if (alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[7] = ser_1x3v_p1_surfvz_quad_7(1, fSkin); 
  } else { 

    fUpwindQuad[7] = ser_1x3v_p1_surfvz_quad_7(-1, fEdge); 
  } 

  fUpwind[0] = 0.3535533905932737*(fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[1] = 0.3535533905932737*(fUpwindQuad[7]-1.0*fUpwindQuad[6]+fUpwindQuad[5]-1.0*fUpwindQuad[4]+fUpwindQuad[3]-1.0*fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[2] = 0.3535533905932737*(fUpwindQuad[7]+fUpwindQuad[6]-1.0*(fUpwindQuad[5]+fUpwindQuad[4])+fUpwindQuad[3]+fUpwindQuad[2]-1.0*(fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[3] = 0.3535533905932737*(fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]-1.0*(fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[4] = 0.3535533905932737*(fUpwindQuad[7]-1.0*(fUpwindQuad[6]+fUpwindQuad[5])+fUpwindQuad[4]+fUpwindQuad[3]-1.0*(fUpwindQuad[2]+fUpwindQuad[1])+fUpwindQuad[0]); 
  fUpwind[5] = 0.3535533905932737*(fUpwindQuad[7]-1.0*fUpwindQuad[6]+fUpwindQuad[5]-1.0*(fUpwindQuad[4]+fUpwindQuad[3])+fUpwindQuad[2]-1.0*fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[6] = 0.3535533905932737*(fUpwindQuad[7]+fUpwindQuad[6]-1.0*(fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2])+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[7] = 0.3535533905932737*(fUpwindQuad[7]-1.0*(fUpwindQuad[6]+fUpwindQuad[5])+fUpwindQuad[4]-1.0*fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 

  drag_incr[0] = 0.3535533905932737*alphaDrSurf[1]*fUpwind[1]+0.3535533905932737*alphaDrSurf[0]*fUpwind[0]; 
  drag_incr[1] = 0.3535533905932737*alphaDrSurf[0]*fUpwind[1]+0.3535533905932737*fUpwind[0]*alphaDrSurf[1]; 
  drag_incr[2] = 0.3535533905932737*alphaDrSurf[1]*fUpwind[4]+0.3535533905932737*alphaDrSurf[0]*fUpwind[2]; 
  drag_incr[3] = 0.3535533905932737*alphaDrSurf[1]*fUpwind[5]+0.3535533905932737*alphaDrSurf[0]*fUpwind[3]; 
  drag_incr[4] = 0.3535533905932737*alphaDrSurf[0]*fUpwind[4]+0.3535533905932737*alphaDrSurf[1]*fUpwind[2]; 
  drag_incr[5] = 0.3535533905932737*alphaDrSurf[0]*fUpwind[5]+0.3535533905932737*alphaDrSurf[1]*fUpwind[3]; 
  drag_incr[6] = 0.3535533905932737*alphaDrSurf[1]*fUpwind[7]+0.3535533905932737*alphaDrSurf[0]*fUpwind[6]; 
  drag_incr[7] = 0.3535533905932737*alphaDrSurf[0]*fUpwind[7]+0.3535533905932737*alphaDrSurf[1]*fUpwind[6]; 

  out[0] += 0.7071067811865475*drag_incr[0]*rdv2; 
  out[1] += 0.7071067811865475*drag_incr[1]*rdv2; 
  out[2] += 0.7071067811865475*drag_incr[2]*rdv2; 
  out[3] += 0.7071067811865475*drag_incr[3]*rdv2; 
  out[4] += 1.224744871391589*drag_incr[0]*rdv2; 
  out[5] += 0.7071067811865475*drag_incr[4]*rdv2; 
  out[6] += 0.7071067811865475*drag_incr[5]*rdv2; 
  out[7] += 0.7071067811865475*drag_incr[6]*rdv2; 
  out[8] += 1.224744871391589*drag_incr[1]*rdv2; 
  out[9] += 1.224744871391589*drag_incr[2]*rdv2; 
  out[10] += 1.224744871391589*drag_incr[3]*rdv2; 
  out[11] += 0.7071067811865475*drag_incr[7]*rdv2; 
  out[12] += 1.224744871391589*drag_incr[4]*rdv2; 
  out[13] += 1.224744871391589*drag_incr[5]*rdv2; 
  out[14] += 1.224744871391589*drag_incr[6]*rdv2; 
  out[15] += 1.224744871391589*drag_incr[7]*rdv2; 

  } else { 

  alphaDrSurf[0] = nuSum[0]*(2.0*w[3]-1.0*dxv[3])-2.0*sumNuUz[0]; 
  alphaDrSurf[1] = nuSum[1]*(2.0*w[3]-1.0*dxv[3])-2.0*sumNuUz[1]; 

  if (alphaDrSurf[0]-alphaDrSurf[1] < 0) { 
    fUpwindQuad[0] = ser_1x3v_p1_surfvz_quad_0(1, fEdge); 
  } else { 
    fUpwindQuad[0] = ser_1x3v_p1_surfvz_quad_0(-1, fSkin); 
  } 
  if (alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[1] = ser_1x3v_p1_surfvz_quad_1(1, fEdge); 
  } else { 
    fUpwindQuad[1] = ser_1x3v_p1_surfvz_quad_1(-1, fSkin); 
  } 
  if (alphaDrSurf[0]-alphaDrSurf[1] < 0) { 
    fUpwindQuad[2] = ser_1x3v_p1_surfvz_quad_2(1, fEdge); 
  } else { 
    fUpwindQuad[2] = ser_1x3v_p1_surfvz_quad_2(-1, fSkin); 
  } 
  if (alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[3] = ser_1x3v_p1_surfvz_quad_3(1, fEdge); 
  } else { 
    fUpwindQuad[3] = ser_1x3v_p1_surfvz_quad_3(-1, fSkin); 
  } 
  if (alphaDrSurf[0]-alphaDrSurf[1] < 0) { 
    fUpwindQuad[4] = ser_1x3v_p1_surfvz_quad_4(1, fEdge); 
  } else { 
    fUpwindQuad[4] = ser_1x3v_p1_surfvz_quad_4(-1, fSkin); 
  } 
  if (alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[5] = ser_1x3v_p1_surfvz_quad_5(1, fEdge); 
  } else { 
    fUpwindQuad[5] = ser_1x3v_p1_surfvz_quad_5(-1, fSkin); 
  } 
  if (alphaDrSurf[0]-alphaDrSurf[1] < 0) { 
    fUpwindQuad[6] = ser_1x3v_p1_surfvz_quad_6(1, fEdge); 
  } else { 
    fUpwindQuad[6] = ser_1x3v_p1_surfvz_quad_6(-1, fSkin); 
  } 
  if (alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[7] = ser_1x3v_p1_surfvz_quad_7(1, fEdge); 
  } else { 
    fUpwindQuad[7] = ser_1x3v_p1_surfvz_quad_7(-1, fSkin); 
  } 

  fUpwind[0] = 0.3535533905932737*(fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[1] = 0.3535533905932737*(fUpwindQuad[7]-1.0*fUpwindQuad[6]+fUpwindQuad[5]-1.0*fUpwindQuad[4]+fUpwindQuad[3]-1.0*fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[2] = 0.3535533905932737*(fUpwindQuad[7]+fUpwindQuad[6]-1.0*(fUpwindQuad[5]+fUpwindQuad[4])+fUpwindQuad[3]+fUpwindQuad[2]-1.0*(fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[3] = 0.3535533905932737*(fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]-1.0*(fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[4] = 0.3535533905932737*(fUpwindQuad[7]-1.0*(fUpwindQuad[6]+fUpwindQuad[5])+fUpwindQuad[4]+fUpwindQuad[3]-1.0*(fUpwindQuad[2]+fUpwindQuad[1])+fUpwindQuad[0]); 
  fUpwind[5] = 0.3535533905932737*(fUpwindQuad[7]-1.0*fUpwindQuad[6]+fUpwindQuad[5]-1.0*(fUpwindQuad[4]+fUpwindQuad[3])+fUpwindQuad[2]-1.0*fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[6] = 0.3535533905932737*(fUpwindQuad[7]+fUpwindQuad[6]-1.0*(fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2])+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[7] = 0.3535533905932737*(fUpwindQuad[7]-1.0*(fUpwindQuad[6]+fUpwindQuad[5])+fUpwindQuad[4]-1.0*fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 

  drag_incr[0] = 0.3535533905932737*alphaDrSurf[1]*fUpwind[1]+0.3535533905932737*alphaDrSurf[0]*fUpwind[0]; 
  drag_incr[1] = 0.3535533905932737*alphaDrSurf[0]*fUpwind[1]+0.3535533905932737*fUpwind[0]*alphaDrSurf[1]; 
  drag_incr[2] = 0.3535533905932737*alphaDrSurf[1]*fUpwind[4]+0.3535533905932737*alphaDrSurf[0]*fUpwind[2]; 
  drag_incr[3] = 0.3535533905932737*alphaDrSurf[1]*fUpwind[5]+0.3535533905932737*alphaDrSurf[0]*fUpwind[3]; 
  drag_incr[4] = 0.3535533905932737*alphaDrSurf[0]*fUpwind[4]+0.3535533905932737*alphaDrSurf[1]*fUpwind[2]; 
  drag_incr[5] = 0.3535533905932737*alphaDrSurf[0]*fUpwind[5]+0.3535533905932737*alphaDrSurf[1]*fUpwind[3]; 
  drag_incr[6] = 0.3535533905932737*alphaDrSurf[1]*fUpwind[7]+0.3535533905932737*alphaDrSurf[0]*fUpwind[6]; 
  drag_incr[7] = 0.3535533905932737*alphaDrSurf[0]*fUpwind[7]+0.3535533905932737*alphaDrSurf[1]*fUpwind[6]; 

  out[0] += 0.7071067811865475*drag_incr[0]*rdv2; 
  out[1] += 0.7071067811865475*drag_incr[1]*rdv2; 
  out[2] += 0.7071067811865475*drag_incr[2]*rdv2; 
  out[3] += 0.7071067811865475*drag_incr[3]*rdv2; 
  out[4] += 1.224744871391589*drag_incr[0]*rdv2; 
  out[5] += 0.7071067811865475*drag_incr[4]*rdv2; 
  out[6] += 0.7071067811865475*drag_incr[5]*rdv2; 
  out[7] += 0.7071067811865475*drag_incr[6]*rdv2; 
  out[8] += 1.224744871391589*drag_incr[1]*rdv2; 
  out[9] += 1.224744871391589*drag_incr[2]*rdv2; 
  out[10] += 1.224744871391589*drag_incr[3]*rdv2; 
  out[11] += 0.7071067811865475*drag_incr[7]*rdv2; 
  out[12] += 1.224744871391589*drag_incr[4]*rdv2; 
  out[13] += 1.224744871391589*drag_incr[5]*rdv2; 
  out[14] += 1.224744871391589*drag_incr[6]*rdv2; 
  out[15] += 1.224744871391589*drag_incr[7]*rdv2; 

  } 
} 
