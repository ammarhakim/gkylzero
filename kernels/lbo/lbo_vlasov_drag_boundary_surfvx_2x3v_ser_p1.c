#include <gkyl_lbo_vlasov_kernels.h> 
#include <gkyl_basis_ser_2x3v_p1_surfvx_quad.h> 
GKYL_CU_DH void lbo_vlasov_drag_boundary_surfvx_2x3v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[5]:         Cell-center coordinates. 
  // dxv[5]:       Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[12]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[4]: sum of thermal speeds squared time their respective collisionalities. 
  // fSkin/Edge:    Distribution function in cells 
  // out:           Incremented distribution function in cell 
  double rdv2 = 2.0/dxv[2]; 

  const double *sumNuUx = &nuUSum[0]; 

  double alphaDrSurf[16] = {0.0}; 
  double fUpwindQuad[16] = {0.0};
  double fUpwind[16] = {0.0};;
  double drag_incr[16] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = nuSum[0]*(2.0*w[2]+dxv[2])-2.0*sumNuUx[0]; 
  alphaDrSurf[1] = nuSum[1]*(2.0*w[2]+dxv[2])-2.0*sumNuUx[1]; 
  alphaDrSurf[2] = 2.0*nuSum[2]*w[2]-2.0*sumNuUx[2]+dxv[2]*nuSum[2]; 
  alphaDrSurf[5] = (2.0*w[2]+dxv[2])*nuSum[3]-2.0*sumNuUx[3]; 

  if (alphaDrSurf[5]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[0] = ser_2x3v_p1_surfvx_quad_0(1, fSkin); 
  } else { 

    fUpwindQuad[0] = ser_2x3v_p1_surfvx_quad_0(-1, fEdge); 
  } 
  if ((-alphaDrSurf[5])-alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[1] = ser_2x3v_p1_surfvx_quad_1(1, fSkin); 
  } else { 

    fUpwindQuad[1] = ser_2x3v_p1_surfvx_quad_1(-1, fEdge); 
  } 
  if ((-alphaDrSurf[5])+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[2] = ser_2x3v_p1_surfvx_quad_2(1, fSkin); 
  } else { 

    fUpwindQuad[2] = ser_2x3v_p1_surfvx_quad_2(-1, fEdge); 
  } 
  if (alphaDrSurf[5]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[3] = ser_2x3v_p1_surfvx_quad_3(1, fSkin); 
  } else { 

    fUpwindQuad[3] = ser_2x3v_p1_surfvx_quad_3(-1, fEdge); 
  } 
  if (alphaDrSurf[5]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[4] = ser_2x3v_p1_surfvx_quad_4(1, fSkin); 
  } else { 

    fUpwindQuad[4] = ser_2x3v_p1_surfvx_quad_4(-1, fEdge); 
  } 
  if ((-alphaDrSurf[5])-alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[5] = ser_2x3v_p1_surfvx_quad_5(1, fSkin); 
  } else { 

    fUpwindQuad[5] = ser_2x3v_p1_surfvx_quad_5(-1, fEdge); 
  } 
  if ((-alphaDrSurf[5])+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[6] = ser_2x3v_p1_surfvx_quad_6(1, fSkin); 
  } else { 

    fUpwindQuad[6] = ser_2x3v_p1_surfvx_quad_6(-1, fEdge); 
  } 
  if (alphaDrSurf[5]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[7] = ser_2x3v_p1_surfvx_quad_7(1, fSkin); 
  } else { 

    fUpwindQuad[7] = ser_2x3v_p1_surfvx_quad_7(-1, fEdge); 
  } 
  if (alphaDrSurf[5]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[8] = ser_2x3v_p1_surfvx_quad_8(1, fSkin); 
  } else { 

    fUpwindQuad[8] = ser_2x3v_p1_surfvx_quad_8(-1, fEdge); 
  } 
  if ((-alphaDrSurf[5])-alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[9] = ser_2x3v_p1_surfvx_quad_9(1, fSkin); 
  } else { 

    fUpwindQuad[9] = ser_2x3v_p1_surfvx_quad_9(-1, fEdge); 
  } 
  if ((-alphaDrSurf[5])+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[10] = ser_2x3v_p1_surfvx_quad_10(1, fSkin); 
  } else { 

    fUpwindQuad[10] = ser_2x3v_p1_surfvx_quad_10(-1, fEdge); 
  } 
  if (alphaDrSurf[5]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[11] = ser_2x3v_p1_surfvx_quad_11(1, fSkin); 
  } else { 

    fUpwindQuad[11] = ser_2x3v_p1_surfvx_quad_11(-1, fEdge); 
  } 
  if (alphaDrSurf[5]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[12] = ser_2x3v_p1_surfvx_quad_12(1, fSkin); 
  } else { 

    fUpwindQuad[12] = ser_2x3v_p1_surfvx_quad_12(-1, fEdge); 
  } 
  if ((-alphaDrSurf[5])-alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[13] = ser_2x3v_p1_surfvx_quad_13(1, fSkin); 
  } else { 

    fUpwindQuad[13] = ser_2x3v_p1_surfvx_quad_13(-1, fEdge); 
  } 
  if ((-alphaDrSurf[5])+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[14] = ser_2x3v_p1_surfvx_quad_14(1, fSkin); 
  } else { 

    fUpwindQuad[14] = ser_2x3v_p1_surfvx_quad_14(-1, fEdge); 
  } 
  if (alphaDrSurf[5]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[15] = ser_2x3v_p1_surfvx_quad_15(1, fSkin); 
  } else { 

    fUpwindQuad[15] = ser_2x3v_p1_surfvx_quad_15(-1, fEdge); 
  } 

  fUpwind[0] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[1] = 0.25*(fUpwindQuad[15]-1.0*fUpwindQuad[14]+fUpwindQuad[13]-1.0*fUpwindQuad[12]+fUpwindQuad[11]-1.0*fUpwindQuad[10]+fUpwindQuad[9]-1.0*fUpwindQuad[8]+fUpwindQuad[7]-1.0*fUpwindQuad[6]+fUpwindQuad[5]-1.0*fUpwindQuad[4]+fUpwindQuad[3]-1.0*fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[2] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]-1.0*(fUpwindQuad[13]+fUpwindQuad[12])+fUpwindQuad[11]+fUpwindQuad[10]-1.0*(fUpwindQuad[9]+fUpwindQuad[8])+fUpwindQuad[7]+fUpwindQuad[6]-1.0*(fUpwindQuad[5]+fUpwindQuad[4])+fUpwindQuad[3]+fUpwindQuad[2]-1.0*(fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[3] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12]-1.0*(fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8])+fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]-1.0*(fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[4] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8]-1.0*(fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[5] = 0.25*(fUpwindQuad[15]-1.0*(fUpwindQuad[14]+fUpwindQuad[13])+fUpwindQuad[12]+fUpwindQuad[11]-1.0*(fUpwindQuad[10]+fUpwindQuad[9])+fUpwindQuad[8]+fUpwindQuad[7]-1.0*(fUpwindQuad[6]+fUpwindQuad[5])+fUpwindQuad[4]+fUpwindQuad[3]-1.0*(fUpwindQuad[2]+fUpwindQuad[1])+fUpwindQuad[0]); 
  fUpwind[6] = 0.25*(fUpwindQuad[15]-1.0*fUpwindQuad[14]+fUpwindQuad[13]-1.0*(fUpwindQuad[12]+fUpwindQuad[11])+fUpwindQuad[10]-1.0*fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]-1.0*fUpwindQuad[6]+fUpwindQuad[5]-1.0*(fUpwindQuad[4]+fUpwindQuad[3])+fUpwindQuad[2]-1.0*fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[7] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]-1.0*(fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10])+fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6]-1.0*(fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2])+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[8] = 0.25*(fUpwindQuad[15]-1.0*fUpwindQuad[14]+fUpwindQuad[13]-1.0*fUpwindQuad[12]+fUpwindQuad[11]-1.0*fUpwindQuad[10]+fUpwindQuad[9]-1.0*(fUpwindQuad[8]+fUpwindQuad[7])+fUpwindQuad[6]-1.0*fUpwindQuad[5]+fUpwindQuad[4]-1.0*fUpwindQuad[3]+fUpwindQuad[2]-1.0*fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[9] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]-1.0*(fUpwindQuad[13]+fUpwindQuad[12])+fUpwindQuad[11]+fUpwindQuad[10]-1.0*(fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6])+fUpwindQuad[5]+fUpwindQuad[4]-1.0*(fUpwindQuad[3]+fUpwindQuad[2])+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[10] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12]-1.0*(fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4])+fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[11] = 0.25*(fUpwindQuad[15]-1.0*(fUpwindQuad[14]+fUpwindQuad[13])+fUpwindQuad[12]-1.0*fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]-1.0*fUpwindQuad[8]+fUpwindQuad[7]-1.0*(fUpwindQuad[6]+fUpwindQuad[5])+fUpwindQuad[4]-1.0*fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[12] = 0.25*(fUpwindQuad[15]-1.0*(fUpwindQuad[14]+fUpwindQuad[13])+fUpwindQuad[12]+fUpwindQuad[11]-1.0*(fUpwindQuad[10]+fUpwindQuad[9])+fUpwindQuad[8]-1.0*fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]-1.0*(fUpwindQuad[4]+fUpwindQuad[3])+fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[13] = 0.25*(fUpwindQuad[15]-1.0*fUpwindQuad[14]+fUpwindQuad[13]-1.0*(fUpwindQuad[12]+fUpwindQuad[11])+fUpwindQuad[10]-1.0*fUpwindQuad[9]+fUpwindQuad[8]-1.0*fUpwindQuad[7]+fUpwindQuad[6]-1.0*fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]-1.0*fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[14] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]-1.0*(fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10])+fUpwindQuad[9]+fUpwindQuad[8]-1.0*(fUpwindQuad[7]+fUpwindQuad[6])+fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2]-1.0*(fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[15] = 0.25*(fUpwindQuad[15]-1.0*(fUpwindQuad[14]+fUpwindQuad[13])+fUpwindQuad[12]-1.0*fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]-1.0*(fUpwindQuad[8]+fUpwindQuad[7])+fUpwindQuad[6]+fUpwindQuad[5]-1.0*fUpwindQuad[4]+fUpwindQuad[3]-1.0*(fUpwindQuad[2]+fUpwindQuad[1])+fUpwindQuad[0]); 

  drag_incr[0] = 0.25*alphaDrSurf[5]*fUpwind[5]+0.25*alphaDrSurf[2]*fUpwind[2]+0.25*alphaDrSurf[1]*fUpwind[1]+0.25*alphaDrSurf[0]*fUpwind[0]; 
  drag_incr[1] = 0.25*alphaDrSurf[2]*fUpwind[5]+0.25*fUpwind[2]*alphaDrSurf[5]+0.25*alphaDrSurf[0]*fUpwind[1]+0.25*fUpwind[0]*alphaDrSurf[1]; 
  drag_incr[2] = 0.25*alphaDrSurf[1]*fUpwind[5]+0.25*fUpwind[1]*alphaDrSurf[5]+0.25*alphaDrSurf[0]*fUpwind[2]+0.25*fUpwind[0]*alphaDrSurf[2]; 
  drag_incr[3] = 0.25*alphaDrSurf[5]*fUpwind[11]+0.25*alphaDrSurf[2]*fUpwind[7]+0.25*alphaDrSurf[1]*fUpwind[6]+0.25*alphaDrSurf[0]*fUpwind[3]; 
  drag_incr[4] = 0.25*alphaDrSurf[5]*fUpwind[12]+0.25*alphaDrSurf[2]*fUpwind[9]+0.25*alphaDrSurf[1]*fUpwind[8]+0.25*alphaDrSurf[0]*fUpwind[4]; 
  drag_incr[5] = 0.25*alphaDrSurf[0]*fUpwind[5]+0.25*fUpwind[0]*alphaDrSurf[5]+0.25*alphaDrSurf[1]*fUpwind[2]+0.25*fUpwind[1]*alphaDrSurf[2]; 
  drag_incr[6] = 0.25*alphaDrSurf[2]*fUpwind[11]+0.25*alphaDrSurf[5]*fUpwind[7]+0.25*alphaDrSurf[0]*fUpwind[6]+0.25*alphaDrSurf[1]*fUpwind[3]; 
  drag_incr[7] = 0.25*alphaDrSurf[1]*fUpwind[11]+0.25*alphaDrSurf[0]*fUpwind[7]+0.25*alphaDrSurf[5]*fUpwind[6]+0.25*alphaDrSurf[2]*fUpwind[3]; 
  drag_incr[8] = 0.25*alphaDrSurf[2]*fUpwind[12]+0.25*alphaDrSurf[5]*fUpwind[9]+0.25*alphaDrSurf[0]*fUpwind[8]+0.25*alphaDrSurf[1]*fUpwind[4]; 
  drag_incr[9] = 0.25*alphaDrSurf[1]*fUpwind[12]+0.25*alphaDrSurf[0]*fUpwind[9]+0.25*alphaDrSurf[5]*fUpwind[8]+0.25*alphaDrSurf[2]*fUpwind[4]; 
  drag_incr[10] = 0.25*alphaDrSurf[5]*fUpwind[15]+0.25*alphaDrSurf[2]*fUpwind[14]+0.25*alphaDrSurf[1]*fUpwind[13]+0.25*alphaDrSurf[0]*fUpwind[10]; 
  drag_incr[11] = 0.25*alphaDrSurf[0]*fUpwind[11]+0.25*alphaDrSurf[1]*fUpwind[7]+0.25*alphaDrSurf[2]*fUpwind[6]+0.25*fUpwind[3]*alphaDrSurf[5]; 
  drag_incr[12] = 0.25*alphaDrSurf[0]*fUpwind[12]+0.25*alphaDrSurf[1]*fUpwind[9]+0.25*alphaDrSurf[2]*fUpwind[8]+0.25*fUpwind[4]*alphaDrSurf[5]; 
  drag_incr[13] = 0.25*alphaDrSurf[2]*fUpwind[15]+0.25*alphaDrSurf[5]*fUpwind[14]+0.25*alphaDrSurf[0]*fUpwind[13]+0.25*alphaDrSurf[1]*fUpwind[10]; 
  drag_incr[14] = 0.25*alphaDrSurf[1]*fUpwind[15]+0.25*alphaDrSurf[0]*fUpwind[14]+0.25*alphaDrSurf[5]*fUpwind[13]+0.25*alphaDrSurf[2]*fUpwind[10]; 
  drag_incr[15] = 0.25*alphaDrSurf[0]*fUpwind[15]+0.25*alphaDrSurf[1]*fUpwind[14]+0.25*alphaDrSurf[2]*fUpwind[13]+0.25*alphaDrSurf[5]*fUpwind[10]; 

  out[0] += 0.7071067811865475*drag_incr[0]*rdv2; 
  out[1] += 0.7071067811865475*drag_incr[1]*rdv2; 
  out[2] += 0.7071067811865475*drag_incr[2]*rdv2; 
  out[3] += 1.224744871391589*drag_incr[0]*rdv2; 
  out[4] += 0.7071067811865475*drag_incr[3]*rdv2; 
  out[5] += 0.7071067811865475*drag_incr[4]*rdv2; 
  out[6] += 0.7071067811865475*drag_incr[5]*rdv2; 
  out[7] += 1.224744871391589*drag_incr[1]*rdv2; 
  out[8] += 1.224744871391589*drag_incr[2]*rdv2; 
  out[9] += 0.7071067811865475*drag_incr[6]*rdv2; 
  out[10] += 0.7071067811865475*drag_incr[7]*rdv2; 
  out[11] += 1.224744871391589*drag_incr[3]*rdv2; 
  out[12] += 0.7071067811865475*drag_incr[8]*rdv2; 
  out[13] += 0.7071067811865475*drag_incr[9]*rdv2; 
  out[14] += 1.224744871391589*drag_incr[4]*rdv2; 
  out[15] += 0.7071067811865475*drag_incr[10]*rdv2; 
  out[16] += 1.224744871391589*drag_incr[5]*rdv2; 
  out[17] += 0.7071067811865475*drag_incr[11]*rdv2; 
  out[18] += 1.224744871391589*drag_incr[6]*rdv2; 
  out[19] += 1.224744871391589*drag_incr[7]*rdv2; 
  out[20] += 0.7071067811865475*drag_incr[12]*rdv2; 
  out[21] += 1.224744871391589*drag_incr[8]*rdv2; 
  out[22] += 1.224744871391589*drag_incr[9]*rdv2; 
  out[23] += 0.7071067811865475*drag_incr[13]*rdv2; 
  out[24] += 0.7071067811865475*drag_incr[14]*rdv2; 
  out[25] += 1.224744871391589*drag_incr[10]*rdv2; 
  out[26] += 1.224744871391589*drag_incr[11]*rdv2; 
  out[27] += 1.224744871391589*drag_incr[12]*rdv2; 
  out[28] += 0.7071067811865475*drag_incr[15]*rdv2; 
  out[29] += 1.224744871391589*drag_incr[13]*rdv2; 
  out[30] += 1.224744871391589*drag_incr[14]*rdv2; 
  out[31] += 1.224744871391589*drag_incr[15]*rdv2; 

  } else { 

  alphaDrSurf[0] = nuSum[0]*(2.0*w[2]-1.0*dxv[2])-2.0*sumNuUx[0]; 
  alphaDrSurf[1] = nuSum[1]*(2.0*w[2]-1.0*dxv[2])-2.0*sumNuUx[1]; 
  alphaDrSurf[2] = 2.0*nuSum[2]*w[2]-2.0*sumNuUx[2]-1.0*dxv[2]*nuSum[2]; 
  alphaDrSurf[5] = (2.0*w[2]-1.0*dxv[2])*nuSum[3]-2.0*sumNuUx[3]; 

  if (alphaDrSurf[5]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[0] = ser_2x3v_p1_surfvx_quad_0(1, fEdge); 
  } else { 
    fUpwindQuad[0] = ser_2x3v_p1_surfvx_quad_0(-1, fSkin); 
  } 
  if ((-alphaDrSurf[5])-alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[1] = ser_2x3v_p1_surfvx_quad_1(1, fEdge); 
  } else { 
    fUpwindQuad[1] = ser_2x3v_p1_surfvx_quad_1(-1, fSkin); 
  } 
  if ((-alphaDrSurf[5])+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[2] = ser_2x3v_p1_surfvx_quad_2(1, fEdge); 
  } else { 
    fUpwindQuad[2] = ser_2x3v_p1_surfvx_quad_2(-1, fSkin); 
  } 
  if (alphaDrSurf[5]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[3] = ser_2x3v_p1_surfvx_quad_3(1, fEdge); 
  } else { 
    fUpwindQuad[3] = ser_2x3v_p1_surfvx_quad_3(-1, fSkin); 
  } 
  if (alphaDrSurf[5]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[4] = ser_2x3v_p1_surfvx_quad_4(1, fEdge); 
  } else { 
    fUpwindQuad[4] = ser_2x3v_p1_surfvx_quad_4(-1, fSkin); 
  } 
  if ((-alphaDrSurf[5])-alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[5] = ser_2x3v_p1_surfvx_quad_5(1, fEdge); 
  } else { 
    fUpwindQuad[5] = ser_2x3v_p1_surfvx_quad_5(-1, fSkin); 
  } 
  if ((-alphaDrSurf[5])+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[6] = ser_2x3v_p1_surfvx_quad_6(1, fEdge); 
  } else { 
    fUpwindQuad[6] = ser_2x3v_p1_surfvx_quad_6(-1, fSkin); 
  } 
  if (alphaDrSurf[5]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[7] = ser_2x3v_p1_surfvx_quad_7(1, fEdge); 
  } else { 
    fUpwindQuad[7] = ser_2x3v_p1_surfvx_quad_7(-1, fSkin); 
  } 
  if (alphaDrSurf[5]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[8] = ser_2x3v_p1_surfvx_quad_8(1, fEdge); 
  } else { 
    fUpwindQuad[8] = ser_2x3v_p1_surfvx_quad_8(-1, fSkin); 
  } 
  if ((-alphaDrSurf[5])-alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[9] = ser_2x3v_p1_surfvx_quad_9(1, fEdge); 
  } else { 
    fUpwindQuad[9] = ser_2x3v_p1_surfvx_quad_9(-1, fSkin); 
  } 
  if ((-alphaDrSurf[5])+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[10] = ser_2x3v_p1_surfvx_quad_10(1, fEdge); 
  } else { 
    fUpwindQuad[10] = ser_2x3v_p1_surfvx_quad_10(-1, fSkin); 
  } 
  if (alphaDrSurf[5]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[11] = ser_2x3v_p1_surfvx_quad_11(1, fEdge); 
  } else { 
    fUpwindQuad[11] = ser_2x3v_p1_surfvx_quad_11(-1, fSkin); 
  } 
  if (alphaDrSurf[5]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[12] = ser_2x3v_p1_surfvx_quad_12(1, fEdge); 
  } else { 
    fUpwindQuad[12] = ser_2x3v_p1_surfvx_quad_12(-1, fSkin); 
  } 
  if ((-alphaDrSurf[5])-alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[13] = ser_2x3v_p1_surfvx_quad_13(1, fEdge); 
  } else { 
    fUpwindQuad[13] = ser_2x3v_p1_surfvx_quad_13(-1, fSkin); 
  } 
  if ((-alphaDrSurf[5])+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[14] = ser_2x3v_p1_surfvx_quad_14(1, fEdge); 
  } else { 
    fUpwindQuad[14] = ser_2x3v_p1_surfvx_quad_14(-1, fSkin); 
  } 
  if (alphaDrSurf[5]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[15] = ser_2x3v_p1_surfvx_quad_15(1, fEdge); 
  } else { 
    fUpwindQuad[15] = ser_2x3v_p1_surfvx_quad_15(-1, fSkin); 
  } 

  fUpwind[0] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[1] = 0.25*(fUpwindQuad[15]-1.0*fUpwindQuad[14]+fUpwindQuad[13]-1.0*fUpwindQuad[12]+fUpwindQuad[11]-1.0*fUpwindQuad[10]+fUpwindQuad[9]-1.0*fUpwindQuad[8]+fUpwindQuad[7]-1.0*fUpwindQuad[6]+fUpwindQuad[5]-1.0*fUpwindQuad[4]+fUpwindQuad[3]-1.0*fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[2] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]-1.0*(fUpwindQuad[13]+fUpwindQuad[12])+fUpwindQuad[11]+fUpwindQuad[10]-1.0*(fUpwindQuad[9]+fUpwindQuad[8])+fUpwindQuad[7]+fUpwindQuad[6]-1.0*(fUpwindQuad[5]+fUpwindQuad[4])+fUpwindQuad[3]+fUpwindQuad[2]-1.0*(fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[3] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12]-1.0*(fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8])+fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]-1.0*(fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[4] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8]-1.0*(fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[5] = 0.25*(fUpwindQuad[15]-1.0*(fUpwindQuad[14]+fUpwindQuad[13])+fUpwindQuad[12]+fUpwindQuad[11]-1.0*(fUpwindQuad[10]+fUpwindQuad[9])+fUpwindQuad[8]+fUpwindQuad[7]-1.0*(fUpwindQuad[6]+fUpwindQuad[5])+fUpwindQuad[4]+fUpwindQuad[3]-1.0*(fUpwindQuad[2]+fUpwindQuad[1])+fUpwindQuad[0]); 
  fUpwind[6] = 0.25*(fUpwindQuad[15]-1.0*fUpwindQuad[14]+fUpwindQuad[13]-1.0*(fUpwindQuad[12]+fUpwindQuad[11])+fUpwindQuad[10]-1.0*fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]-1.0*fUpwindQuad[6]+fUpwindQuad[5]-1.0*(fUpwindQuad[4]+fUpwindQuad[3])+fUpwindQuad[2]-1.0*fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[7] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]-1.0*(fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10])+fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6]-1.0*(fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2])+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[8] = 0.25*(fUpwindQuad[15]-1.0*fUpwindQuad[14]+fUpwindQuad[13]-1.0*fUpwindQuad[12]+fUpwindQuad[11]-1.0*fUpwindQuad[10]+fUpwindQuad[9]-1.0*(fUpwindQuad[8]+fUpwindQuad[7])+fUpwindQuad[6]-1.0*fUpwindQuad[5]+fUpwindQuad[4]-1.0*fUpwindQuad[3]+fUpwindQuad[2]-1.0*fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[9] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]-1.0*(fUpwindQuad[13]+fUpwindQuad[12])+fUpwindQuad[11]+fUpwindQuad[10]-1.0*(fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6])+fUpwindQuad[5]+fUpwindQuad[4]-1.0*(fUpwindQuad[3]+fUpwindQuad[2])+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[10] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12]-1.0*(fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4])+fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[11] = 0.25*(fUpwindQuad[15]-1.0*(fUpwindQuad[14]+fUpwindQuad[13])+fUpwindQuad[12]-1.0*fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]-1.0*fUpwindQuad[8]+fUpwindQuad[7]-1.0*(fUpwindQuad[6]+fUpwindQuad[5])+fUpwindQuad[4]-1.0*fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[12] = 0.25*(fUpwindQuad[15]-1.0*(fUpwindQuad[14]+fUpwindQuad[13])+fUpwindQuad[12]+fUpwindQuad[11]-1.0*(fUpwindQuad[10]+fUpwindQuad[9])+fUpwindQuad[8]-1.0*fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]-1.0*(fUpwindQuad[4]+fUpwindQuad[3])+fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[13] = 0.25*(fUpwindQuad[15]-1.0*fUpwindQuad[14]+fUpwindQuad[13]-1.0*(fUpwindQuad[12]+fUpwindQuad[11])+fUpwindQuad[10]-1.0*fUpwindQuad[9]+fUpwindQuad[8]-1.0*fUpwindQuad[7]+fUpwindQuad[6]-1.0*fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]-1.0*fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[14] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]-1.0*(fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10])+fUpwindQuad[9]+fUpwindQuad[8]-1.0*(fUpwindQuad[7]+fUpwindQuad[6])+fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2]-1.0*(fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[15] = 0.25*(fUpwindQuad[15]-1.0*(fUpwindQuad[14]+fUpwindQuad[13])+fUpwindQuad[12]-1.0*fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]-1.0*(fUpwindQuad[8]+fUpwindQuad[7])+fUpwindQuad[6]+fUpwindQuad[5]-1.0*fUpwindQuad[4]+fUpwindQuad[3]-1.0*(fUpwindQuad[2]+fUpwindQuad[1])+fUpwindQuad[0]); 

  drag_incr[0] = 0.25*alphaDrSurf[5]*fUpwind[5]+0.25*alphaDrSurf[2]*fUpwind[2]+0.25*alphaDrSurf[1]*fUpwind[1]+0.25*alphaDrSurf[0]*fUpwind[0]; 
  drag_incr[1] = 0.25*alphaDrSurf[2]*fUpwind[5]+0.25*fUpwind[2]*alphaDrSurf[5]+0.25*alphaDrSurf[0]*fUpwind[1]+0.25*fUpwind[0]*alphaDrSurf[1]; 
  drag_incr[2] = 0.25*alphaDrSurf[1]*fUpwind[5]+0.25*fUpwind[1]*alphaDrSurf[5]+0.25*alphaDrSurf[0]*fUpwind[2]+0.25*fUpwind[0]*alphaDrSurf[2]; 
  drag_incr[3] = 0.25*alphaDrSurf[5]*fUpwind[11]+0.25*alphaDrSurf[2]*fUpwind[7]+0.25*alphaDrSurf[1]*fUpwind[6]+0.25*alphaDrSurf[0]*fUpwind[3]; 
  drag_incr[4] = 0.25*alphaDrSurf[5]*fUpwind[12]+0.25*alphaDrSurf[2]*fUpwind[9]+0.25*alphaDrSurf[1]*fUpwind[8]+0.25*alphaDrSurf[0]*fUpwind[4]; 
  drag_incr[5] = 0.25*alphaDrSurf[0]*fUpwind[5]+0.25*fUpwind[0]*alphaDrSurf[5]+0.25*alphaDrSurf[1]*fUpwind[2]+0.25*fUpwind[1]*alphaDrSurf[2]; 
  drag_incr[6] = 0.25*alphaDrSurf[2]*fUpwind[11]+0.25*alphaDrSurf[5]*fUpwind[7]+0.25*alphaDrSurf[0]*fUpwind[6]+0.25*alphaDrSurf[1]*fUpwind[3]; 
  drag_incr[7] = 0.25*alphaDrSurf[1]*fUpwind[11]+0.25*alphaDrSurf[0]*fUpwind[7]+0.25*alphaDrSurf[5]*fUpwind[6]+0.25*alphaDrSurf[2]*fUpwind[3]; 
  drag_incr[8] = 0.25*alphaDrSurf[2]*fUpwind[12]+0.25*alphaDrSurf[5]*fUpwind[9]+0.25*alphaDrSurf[0]*fUpwind[8]+0.25*alphaDrSurf[1]*fUpwind[4]; 
  drag_incr[9] = 0.25*alphaDrSurf[1]*fUpwind[12]+0.25*alphaDrSurf[0]*fUpwind[9]+0.25*alphaDrSurf[5]*fUpwind[8]+0.25*alphaDrSurf[2]*fUpwind[4]; 
  drag_incr[10] = 0.25*alphaDrSurf[5]*fUpwind[15]+0.25*alphaDrSurf[2]*fUpwind[14]+0.25*alphaDrSurf[1]*fUpwind[13]+0.25*alphaDrSurf[0]*fUpwind[10]; 
  drag_incr[11] = 0.25*alphaDrSurf[0]*fUpwind[11]+0.25*alphaDrSurf[1]*fUpwind[7]+0.25*alphaDrSurf[2]*fUpwind[6]+0.25*fUpwind[3]*alphaDrSurf[5]; 
  drag_incr[12] = 0.25*alphaDrSurf[0]*fUpwind[12]+0.25*alphaDrSurf[1]*fUpwind[9]+0.25*alphaDrSurf[2]*fUpwind[8]+0.25*fUpwind[4]*alphaDrSurf[5]; 
  drag_incr[13] = 0.25*alphaDrSurf[2]*fUpwind[15]+0.25*alphaDrSurf[5]*fUpwind[14]+0.25*alphaDrSurf[0]*fUpwind[13]+0.25*alphaDrSurf[1]*fUpwind[10]; 
  drag_incr[14] = 0.25*alphaDrSurf[1]*fUpwind[15]+0.25*alphaDrSurf[0]*fUpwind[14]+0.25*alphaDrSurf[5]*fUpwind[13]+0.25*alphaDrSurf[2]*fUpwind[10]; 
  drag_incr[15] = 0.25*alphaDrSurf[0]*fUpwind[15]+0.25*alphaDrSurf[1]*fUpwind[14]+0.25*alphaDrSurf[2]*fUpwind[13]+0.25*alphaDrSurf[5]*fUpwind[10]; 

  out[0] += -0.7071067811865475*drag_incr[0]*rdv2; 
  out[1] += -0.7071067811865475*drag_incr[1]*rdv2; 
  out[2] += -0.7071067811865475*drag_incr[2]*rdv2; 
  out[3] += 1.224744871391589*drag_incr[0]*rdv2; 
  out[4] += -0.7071067811865475*drag_incr[3]*rdv2; 
  out[5] += -0.7071067811865475*drag_incr[4]*rdv2; 
  out[6] += -0.7071067811865475*drag_incr[5]*rdv2; 
  out[7] += 1.224744871391589*drag_incr[1]*rdv2; 
  out[8] += 1.224744871391589*drag_incr[2]*rdv2; 
  out[9] += -0.7071067811865475*drag_incr[6]*rdv2; 
  out[10] += -0.7071067811865475*drag_incr[7]*rdv2; 
  out[11] += 1.224744871391589*drag_incr[3]*rdv2; 
  out[12] += -0.7071067811865475*drag_incr[8]*rdv2; 
  out[13] += -0.7071067811865475*drag_incr[9]*rdv2; 
  out[14] += 1.224744871391589*drag_incr[4]*rdv2; 
  out[15] += -0.7071067811865475*drag_incr[10]*rdv2; 
  out[16] += 1.224744871391589*drag_incr[5]*rdv2; 
  out[17] += -0.7071067811865475*drag_incr[11]*rdv2; 
  out[18] += 1.224744871391589*drag_incr[6]*rdv2; 
  out[19] += 1.224744871391589*drag_incr[7]*rdv2; 
  out[20] += -0.7071067811865475*drag_incr[12]*rdv2; 
  out[21] += 1.224744871391589*drag_incr[8]*rdv2; 
  out[22] += 1.224744871391589*drag_incr[9]*rdv2; 
  out[23] += -0.7071067811865475*drag_incr[13]*rdv2; 
  out[24] += -0.7071067811865475*drag_incr[14]*rdv2; 
  out[25] += 1.224744871391589*drag_incr[10]*rdv2; 
  out[26] += 1.224744871391589*drag_incr[11]*rdv2; 
  out[27] += 1.224744871391589*drag_incr[12]*rdv2; 
  out[28] += -0.7071067811865475*drag_incr[15]*rdv2; 
  out[29] += 1.224744871391589*drag_incr[13]*rdv2; 
  out[30] += 1.224744871391589*drag_incr[14]*rdv2; 
  out[31] += 1.224744871391589*drag_incr[15]*rdv2; 

  } 
} 
