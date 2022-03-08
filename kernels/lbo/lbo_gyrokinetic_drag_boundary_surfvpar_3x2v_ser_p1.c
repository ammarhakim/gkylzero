#include <gkyl_lbo_gyrokinetic_kernels.h> 
#include <gkyl_basis_ser_5x_p1_surfx4_quad.h> 
#include <gkyl_basis_ser_5x_p1_upwind.h> 
GKYL_CU_DH void lbo_gyrokinetic_drag_boundary_surfvpar_3x2v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[5]:     cell-center coordinates. 
  // dxv[5]:   cell spacing. 
  // m_:        species mass.
  // bmag_inv:  1/(magnetic field magnitude). 
  // nuSum:     collisionalities added (self and cross species collisionalities). 
  // nuUSum[16]:sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum: sum of thermal speeds squared time their respective collisionalities. 
  // fSkin/Edge:    Distribution function in cells 
  // out:           Incremented distribution function in cell 
  double rdv2 = 2.0/dxv[3]; 


  double alphaDrSurf[16] = {0.0}; 
  double fUpwindQuad[16] = {0.0};
  double fUpwind[16] = {0.0};;
  double drag_incr[16] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = 0.5*(nuSum[0]*(2.828427124746191*w[3]+1.414213562373095*dxv[3])-2.828427124746191*nuUSum[0]); 
  alphaDrSurf[1] = 0.5*(nuSum[1]*(2.828427124746191*w[3]+1.414213562373095*dxv[3])-2.828427124746191*nuUSum[1]); 
  alphaDrSurf[2] = 0.5*(nuSum[2]*(2.828427124746191*w[3]+1.414213562373095*dxv[3])-2.828427124746191*nuUSum[2]); 
  alphaDrSurf[3] = 0.5*(2.828427124746191*nuSum[3]*w[3]-2.828427124746191*nuUSum[3]+1.414213562373095*dxv[3]*nuSum[3]); 
  alphaDrSurf[5] = -0.5*(2.828427124746191*nuUSum[4]+((-2.828427124746191*w[3])-1.414213562373095*dxv[3])*nuSum[4]); 
  alphaDrSurf[6] = -0.5*(2.828427124746191*nuUSum[5]+((-2.828427124746191*w[3])-1.414213562373095*dxv[3])*nuSum[5]); 
  alphaDrSurf[7] = -0.5*(2.828427124746191*nuUSum[6]+((-2.828427124746191*w[3])-1.414213562373095*dxv[3])*nuSum[6]); 
  alphaDrSurf[11] = -0.5*(2.828427124746191*nuUSum[7]+((-2.828427124746191*w[3])-1.414213562373095*dxv[3])*nuSum[7]); 

  if ((-alphaDrSurf[11])+alphaDrSurf[7]+alphaDrSurf[6]+alphaDrSurf[5]-alphaDrSurf[3]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[0] = ser_5x_p1_surfx4_quad_0_r(fSkin); 
    fUpwindQuad[8] = ser_5x_p1_surfx4_quad_8_r(fSkin); 
  } else { 
    fUpwindQuad[0] = ser_5x_p1_surfx4_quad_0_l(fEdge); 
    fUpwindQuad[8] = ser_5x_p1_surfx4_quad_8_l(fEdge); 
  } 
  if (alphaDrSurf[11]+alphaDrSurf[7]-alphaDrSurf[6]-alphaDrSurf[5]-alphaDrSurf[3]-alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[1] = ser_5x_p1_surfx4_quad_1_r(fSkin); 
    fUpwindQuad[9] = ser_5x_p1_surfx4_quad_9_r(fSkin); 
  } else { 
    fUpwindQuad[1] = ser_5x_p1_surfx4_quad_1_l(fEdge); 
    fUpwindQuad[9] = ser_5x_p1_surfx4_quad_9_l(fEdge); 
  } 
  if (alphaDrSurf[11]-alphaDrSurf[7]+alphaDrSurf[6]-alphaDrSurf[5]-alphaDrSurf[3]+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[2] = ser_5x_p1_surfx4_quad_2_r(fSkin); 
    fUpwindQuad[10] = ser_5x_p1_surfx4_quad_10_r(fSkin); 
  } else { 
    fUpwindQuad[2] = ser_5x_p1_surfx4_quad_2_l(fEdge); 
    fUpwindQuad[10] = ser_5x_p1_surfx4_quad_10_l(fEdge); 
  } 
  if ((-alphaDrSurf[11])-alphaDrSurf[7]-alphaDrSurf[6]+alphaDrSurf[5]-alphaDrSurf[3]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[3] = ser_5x_p1_surfx4_quad_3_r(fSkin); 
    fUpwindQuad[11] = ser_5x_p1_surfx4_quad_11_r(fSkin); 
  } else { 
    fUpwindQuad[3] = ser_5x_p1_surfx4_quad_3_l(fEdge); 
    fUpwindQuad[11] = ser_5x_p1_surfx4_quad_11_l(fEdge); 
  } 
  if (alphaDrSurf[11]-alphaDrSurf[7]-alphaDrSurf[6]+alphaDrSurf[5]+alphaDrSurf[3]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[4] = ser_5x_p1_surfx4_quad_4_r(fSkin); 
    fUpwindQuad[12] = ser_5x_p1_surfx4_quad_12_r(fSkin); 
  } else { 
    fUpwindQuad[4] = ser_5x_p1_surfx4_quad_4_l(fEdge); 
    fUpwindQuad[12] = ser_5x_p1_surfx4_quad_12_l(fEdge); 
  } 
  if ((-alphaDrSurf[11])-alphaDrSurf[7]+alphaDrSurf[6]-alphaDrSurf[5]+alphaDrSurf[3]-alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[5] = ser_5x_p1_surfx4_quad_5_r(fSkin); 
    fUpwindQuad[13] = ser_5x_p1_surfx4_quad_13_r(fSkin); 
  } else { 
    fUpwindQuad[5] = ser_5x_p1_surfx4_quad_5_l(fEdge); 
    fUpwindQuad[13] = ser_5x_p1_surfx4_quad_13_l(fEdge); 
  } 
  if ((-alphaDrSurf[11])+alphaDrSurf[7]-alphaDrSurf[6]-alphaDrSurf[5]+alphaDrSurf[3]+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[6] = ser_5x_p1_surfx4_quad_6_r(fSkin); 
    fUpwindQuad[14] = ser_5x_p1_surfx4_quad_14_r(fSkin); 
  } else { 
    fUpwindQuad[6] = ser_5x_p1_surfx4_quad_6_l(fEdge); 
    fUpwindQuad[14] = ser_5x_p1_surfx4_quad_14_l(fEdge); 
  } 
  if (alphaDrSurf[11]+alphaDrSurf[7]+alphaDrSurf[6]+alphaDrSurf[5]+alphaDrSurf[3]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[7] = ser_5x_p1_surfx4_quad_7_r(fSkin); 
    fUpwindQuad[15] = ser_5x_p1_surfx4_quad_15_r(fSkin); 
  } else { 
    fUpwindQuad[7] = ser_5x_p1_surfx4_quad_7_l(fEdge); 
    fUpwindQuad[15] = ser_5x_p1_surfx4_quad_15_l(fEdge); 
  } 

  // Project nodal basis back onto modal basis. 
  ser_5x_p1_upwind(fUpwindQuad, fUpwind); 

  drag_incr[0] = 0.25*alphaDrSurf[11]*fUpwind[11]+0.25*alphaDrSurf[7]*fUpwind[7]+0.25*alphaDrSurf[6]*fUpwind[6]+0.25*alphaDrSurf[5]*fUpwind[5]+0.25*alphaDrSurf[3]*fUpwind[3]+0.25*alphaDrSurf[2]*fUpwind[2]+0.25*alphaDrSurf[1]*fUpwind[1]+0.25*alphaDrSurf[0]*fUpwind[0]; 
  drag_incr[1] = 0.25*alphaDrSurf[7]*fUpwind[11]+0.25*fUpwind[7]*alphaDrSurf[11]+0.25*alphaDrSurf[3]*fUpwind[6]+0.25*fUpwind[3]*alphaDrSurf[6]+0.25*alphaDrSurf[2]*fUpwind[5]+0.25*fUpwind[2]*alphaDrSurf[5]+0.25*alphaDrSurf[0]*fUpwind[1]+0.25*fUpwind[0]*alphaDrSurf[1]; 
  drag_incr[2] = 0.25*alphaDrSurf[6]*fUpwind[11]+0.25*fUpwind[6]*alphaDrSurf[11]+0.25*alphaDrSurf[3]*fUpwind[7]+0.25*fUpwind[3]*alphaDrSurf[7]+0.25*alphaDrSurf[1]*fUpwind[5]+0.25*fUpwind[1]*alphaDrSurf[5]+0.25*alphaDrSurf[0]*fUpwind[2]+0.25*fUpwind[0]*alphaDrSurf[2]; 
  drag_incr[3] = 0.25*alphaDrSurf[5]*fUpwind[11]+0.25*fUpwind[5]*alphaDrSurf[11]+0.25*alphaDrSurf[2]*fUpwind[7]+0.25*fUpwind[2]*alphaDrSurf[7]+0.25*alphaDrSurf[1]*fUpwind[6]+0.25*fUpwind[1]*alphaDrSurf[6]+0.25*alphaDrSurf[0]*fUpwind[3]+0.25*fUpwind[0]*alphaDrSurf[3]; 
  drag_incr[4] = 0.25*alphaDrSurf[11]*fUpwind[15]+0.25*alphaDrSurf[7]*fUpwind[14]+0.25*alphaDrSurf[6]*fUpwind[13]+0.25*alphaDrSurf[5]*fUpwind[12]+0.25*alphaDrSurf[3]*fUpwind[10]+0.25*alphaDrSurf[2]*fUpwind[9]+0.25*alphaDrSurf[1]*fUpwind[8]+0.25*alphaDrSurf[0]*fUpwind[4]; 
  drag_incr[5] = 0.25*alphaDrSurf[3]*fUpwind[11]+0.25*fUpwind[3]*alphaDrSurf[11]+0.25*alphaDrSurf[6]*fUpwind[7]+0.25*fUpwind[6]*alphaDrSurf[7]+0.25*alphaDrSurf[0]*fUpwind[5]+0.25*fUpwind[0]*alphaDrSurf[5]+0.25*alphaDrSurf[1]*fUpwind[2]+0.25*fUpwind[1]*alphaDrSurf[2]; 
  drag_incr[6] = 0.25*alphaDrSurf[2]*fUpwind[11]+0.25*fUpwind[2]*alphaDrSurf[11]+0.25*alphaDrSurf[5]*fUpwind[7]+0.25*fUpwind[5]*alphaDrSurf[7]+0.25*alphaDrSurf[0]*fUpwind[6]+0.25*fUpwind[0]*alphaDrSurf[6]+0.25*alphaDrSurf[1]*fUpwind[3]+0.25*fUpwind[1]*alphaDrSurf[3]; 
  drag_incr[7] = 0.25*alphaDrSurf[1]*fUpwind[11]+0.25*fUpwind[1]*alphaDrSurf[11]+0.25*alphaDrSurf[0]*fUpwind[7]+0.25*fUpwind[0]*alphaDrSurf[7]+0.25*alphaDrSurf[5]*fUpwind[6]+0.25*fUpwind[5]*alphaDrSurf[6]+0.25*alphaDrSurf[2]*fUpwind[3]+0.25*fUpwind[2]*alphaDrSurf[3]; 
  drag_incr[8] = 0.25*alphaDrSurf[7]*fUpwind[15]+0.25*alphaDrSurf[11]*fUpwind[14]+0.25*alphaDrSurf[3]*fUpwind[13]+0.25*alphaDrSurf[2]*fUpwind[12]+0.25*alphaDrSurf[6]*fUpwind[10]+0.25*alphaDrSurf[5]*fUpwind[9]+0.25*alphaDrSurf[0]*fUpwind[8]+0.25*alphaDrSurf[1]*fUpwind[4]; 
  drag_incr[9] = 0.25*alphaDrSurf[6]*fUpwind[15]+0.25*alphaDrSurf[3]*fUpwind[14]+0.25*alphaDrSurf[11]*fUpwind[13]+0.25*alphaDrSurf[1]*fUpwind[12]+0.25*alphaDrSurf[7]*fUpwind[10]+0.25*alphaDrSurf[0]*fUpwind[9]+0.25*alphaDrSurf[5]*fUpwind[8]+0.25*alphaDrSurf[2]*fUpwind[4]; 
  drag_incr[10] = 0.25*alphaDrSurf[5]*fUpwind[15]+0.25*alphaDrSurf[2]*fUpwind[14]+0.25*alphaDrSurf[1]*fUpwind[13]+0.25*alphaDrSurf[11]*fUpwind[12]+0.25*alphaDrSurf[0]*fUpwind[10]+0.25*alphaDrSurf[7]*fUpwind[9]+0.25*alphaDrSurf[6]*fUpwind[8]+0.25*alphaDrSurf[3]*fUpwind[4]; 
  drag_incr[11] = 0.25*alphaDrSurf[0]*fUpwind[11]+0.25*fUpwind[0]*alphaDrSurf[11]+0.25*alphaDrSurf[1]*fUpwind[7]+0.25*fUpwind[1]*alphaDrSurf[7]+0.25*alphaDrSurf[2]*fUpwind[6]+0.25*fUpwind[2]*alphaDrSurf[6]+0.25*alphaDrSurf[3]*fUpwind[5]+0.25*fUpwind[3]*alphaDrSurf[5]; 
  drag_incr[12] = 0.25*alphaDrSurf[3]*fUpwind[15]+0.25*alphaDrSurf[6]*fUpwind[14]+0.25*alphaDrSurf[7]*fUpwind[13]+0.25*alphaDrSurf[0]*fUpwind[12]+0.25*fUpwind[10]*alphaDrSurf[11]+0.25*alphaDrSurf[1]*fUpwind[9]+0.25*alphaDrSurf[2]*fUpwind[8]+0.25*fUpwind[4]*alphaDrSurf[5]; 
  drag_incr[13] = 0.25*alphaDrSurf[2]*fUpwind[15]+0.25*alphaDrSurf[5]*fUpwind[14]+0.25*alphaDrSurf[0]*fUpwind[13]+0.25*alphaDrSurf[7]*fUpwind[12]+0.25*fUpwind[9]*alphaDrSurf[11]+0.25*alphaDrSurf[1]*fUpwind[10]+0.25*alphaDrSurf[3]*fUpwind[8]+0.25*fUpwind[4]*alphaDrSurf[6]; 
  drag_incr[14] = 0.25*alphaDrSurf[1]*fUpwind[15]+0.25*alphaDrSurf[0]*fUpwind[14]+0.25*alphaDrSurf[5]*fUpwind[13]+0.25*alphaDrSurf[6]*fUpwind[12]+0.25*fUpwind[8]*alphaDrSurf[11]+0.25*alphaDrSurf[2]*fUpwind[10]+0.25*alphaDrSurf[3]*fUpwind[9]+0.25*fUpwind[4]*alphaDrSurf[7]; 
  drag_incr[15] = 0.25*alphaDrSurf[0]*fUpwind[15]+0.25*alphaDrSurf[1]*fUpwind[14]+0.25*alphaDrSurf[2]*fUpwind[13]+0.25*alphaDrSurf[3]*fUpwind[12]+0.25*fUpwind[4]*alphaDrSurf[11]+0.25*alphaDrSurf[5]*fUpwind[10]+0.25*alphaDrSurf[6]*fUpwind[9]+0.25*alphaDrSurf[7]*fUpwind[8]; 

  out[0] += 0.7071067811865475*drag_incr[0]*rdv2; 
  out[1] += 0.7071067811865475*drag_incr[1]*rdv2; 
  out[2] += 0.7071067811865475*drag_incr[2]*rdv2; 
  out[3] += 0.7071067811865475*drag_incr[3]*rdv2; 
  out[4] += 1.224744871391589*drag_incr[0]*rdv2; 
  out[5] += 0.7071067811865475*drag_incr[4]*rdv2; 
  out[6] += 0.7071067811865475*drag_incr[5]*rdv2; 
  out[7] += 0.7071067811865475*drag_incr[6]*rdv2; 
  out[8] += 0.7071067811865475*drag_incr[7]*rdv2; 
  out[9] += 1.224744871391589*drag_incr[1]*rdv2; 
  out[10] += 1.224744871391589*drag_incr[2]*rdv2; 
  out[11] += 1.224744871391589*drag_incr[3]*rdv2; 
  out[12] += 0.7071067811865475*drag_incr[8]*rdv2; 
  out[13] += 0.7071067811865475*drag_incr[9]*rdv2; 
  out[14] += 0.7071067811865475*drag_incr[10]*rdv2; 
  out[15] += 1.224744871391589*drag_incr[4]*rdv2; 
  out[16] += 0.7071067811865475*drag_incr[11]*rdv2; 
  out[17] += 1.224744871391589*drag_incr[5]*rdv2; 
  out[18] += 1.224744871391589*drag_incr[6]*rdv2; 
  out[19] += 1.224744871391589*drag_incr[7]*rdv2; 
  out[20] += 0.7071067811865475*drag_incr[12]*rdv2; 
  out[21] += 0.7071067811865475*drag_incr[13]*rdv2; 
  out[22] += 0.7071067811865475*drag_incr[14]*rdv2; 
  out[23] += 1.224744871391589*drag_incr[8]*rdv2; 
  out[24] += 1.224744871391589*drag_incr[9]*rdv2; 
  out[25] += 1.224744871391589*drag_incr[10]*rdv2; 
  out[26] += 1.224744871391589*drag_incr[11]*rdv2; 
  out[27] += 0.7071067811865475*drag_incr[15]*rdv2; 
  out[28] += 1.224744871391589*drag_incr[12]*rdv2; 
  out[29] += 1.224744871391589*drag_incr[13]*rdv2; 
  out[30] += 1.224744871391589*drag_incr[14]*rdv2; 
  out[31] += 1.224744871391589*drag_incr[15]*rdv2; 

  } else { 

  alphaDrSurf[0] = 0.5*(nuSum[0]*(2.828427124746191*w[3]-1.414213562373095*dxv[3])-2.828427124746191*nuUSum[0]); 
  alphaDrSurf[1] = 0.5*(nuSum[1]*(2.828427124746191*w[3]-1.414213562373095*dxv[3])-2.828427124746191*nuUSum[1]); 
  alphaDrSurf[2] = 0.5*(nuSum[2]*(2.828427124746191*w[3]-1.414213562373095*dxv[3])-2.828427124746191*nuUSum[2]); 
  alphaDrSurf[3] = 0.5*(2.828427124746191*nuSum[3]*w[3]-2.828427124746191*nuUSum[3]-1.414213562373095*dxv[3]*nuSum[3]); 
  alphaDrSurf[5] = -0.5*(2.828427124746191*nuUSum[4]+(1.414213562373095*dxv[3]-2.828427124746191*w[3])*nuSum[4]); 
  alphaDrSurf[6] = -0.5*(2.828427124746191*nuUSum[5]+(1.414213562373095*dxv[3]-2.828427124746191*w[3])*nuSum[5]); 
  alphaDrSurf[7] = -0.5*(2.828427124746191*nuUSum[6]+(1.414213562373095*dxv[3]-2.828427124746191*w[3])*nuSum[6]); 
  alphaDrSurf[11] = -0.5*(2.828427124746191*nuUSum[7]+(1.414213562373095*dxv[3]-2.828427124746191*w[3])*nuSum[7]); 

  if ((-alphaDrSurf[11])+alphaDrSurf[7]+alphaDrSurf[6]+alphaDrSurf[5]-alphaDrSurf[3]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[0] = ser_5x_p1_surfx4_quad_0_r(fEdge); 
    fUpwindQuad[8] = ser_5x_p1_surfx4_quad_8_r(fEdge); 
  } else { 
    fUpwindQuad[0] = ser_5x_p1_surfx4_quad_0_l(fSkin); 
    fUpwindQuad[8] = ser_5x_p1_surfx4_quad_8_l(fSkin); 
  } 
  if (alphaDrSurf[11]+alphaDrSurf[7]-alphaDrSurf[6]-alphaDrSurf[5]-alphaDrSurf[3]-alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[1] = ser_5x_p1_surfx4_quad_1_r(fEdge); 
    fUpwindQuad[9] = ser_5x_p1_surfx4_quad_9_r(fEdge); 
  } else { 
    fUpwindQuad[1] = ser_5x_p1_surfx4_quad_1_l(fSkin); 
    fUpwindQuad[9] = ser_5x_p1_surfx4_quad_9_l(fSkin); 
  } 
  if (alphaDrSurf[11]-alphaDrSurf[7]+alphaDrSurf[6]-alphaDrSurf[5]-alphaDrSurf[3]+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[2] = ser_5x_p1_surfx4_quad_2_r(fEdge); 
    fUpwindQuad[10] = ser_5x_p1_surfx4_quad_10_r(fEdge); 
  } else { 
    fUpwindQuad[2] = ser_5x_p1_surfx4_quad_2_l(fSkin); 
    fUpwindQuad[10] = ser_5x_p1_surfx4_quad_10_l(fSkin); 
  } 
  if ((-alphaDrSurf[11])-alphaDrSurf[7]-alphaDrSurf[6]+alphaDrSurf[5]-alphaDrSurf[3]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[3] = ser_5x_p1_surfx4_quad_3_r(fEdge); 
    fUpwindQuad[11] = ser_5x_p1_surfx4_quad_11_r(fEdge); 
  } else { 
    fUpwindQuad[3] = ser_5x_p1_surfx4_quad_3_l(fSkin); 
    fUpwindQuad[11] = ser_5x_p1_surfx4_quad_11_l(fSkin); 
  } 
  if (alphaDrSurf[11]-alphaDrSurf[7]-alphaDrSurf[6]+alphaDrSurf[5]+alphaDrSurf[3]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[4] = ser_5x_p1_surfx4_quad_4_r(fEdge); 
    fUpwindQuad[12] = ser_5x_p1_surfx4_quad_12_r(fEdge); 
  } else { 
    fUpwindQuad[4] = ser_5x_p1_surfx4_quad_4_l(fSkin); 
    fUpwindQuad[12] = ser_5x_p1_surfx4_quad_12_l(fSkin); 
  } 
  if ((-alphaDrSurf[11])-alphaDrSurf[7]+alphaDrSurf[6]-alphaDrSurf[5]+alphaDrSurf[3]-alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[5] = ser_5x_p1_surfx4_quad_5_r(fEdge); 
    fUpwindQuad[13] = ser_5x_p1_surfx4_quad_13_r(fEdge); 
  } else { 
    fUpwindQuad[5] = ser_5x_p1_surfx4_quad_5_l(fSkin); 
    fUpwindQuad[13] = ser_5x_p1_surfx4_quad_13_l(fSkin); 
  } 
  if ((-alphaDrSurf[11])+alphaDrSurf[7]-alphaDrSurf[6]-alphaDrSurf[5]+alphaDrSurf[3]+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[6] = ser_5x_p1_surfx4_quad_6_r(fEdge); 
    fUpwindQuad[14] = ser_5x_p1_surfx4_quad_14_r(fEdge); 
  } else { 
    fUpwindQuad[6] = ser_5x_p1_surfx4_quad_6_l(fSkin); 
    fUpwindQuad[14] = ser_5x_p1_surfx4_quad_14_l(fSkin); 
  } 
  if (alphaDrSurf[11]+alphaDrSurf[7]+alphaDrSurf[6]+alphaDrSurf[5]+alphaDrSurf[3]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[7] = ser_5x_p1_surfx4_quad_7_r(fEdge); 
    fUpwindQuad[15] = ser_5x_p1_surfx4_quad_15_r(fEdge); 
  } else { 
    fUpwindQuad[7] = ser_5x_p1_surfx4_quad_7_l(fSkin); 
    fUpwindQuad[15] = ser_5x_p1_surfx4_quad_15_l(fSkin); 
  } 

  // Project nodal basis back onto modal basis. 
  ser_5x_p1_upwind(fUpwindQuad, fUpwind); 

  drag_incr[0] = 0.25*alphaDrSurf[11]*fUpwind[11]+0.25*alphaDrSurf[7]*fUpwind[7]+0.25*alphaDrSurf[6]*fUpwind[6]+0.25*alphaDrSurf[5]*fUpwind[5]+0.25*alphaDrSurf[3]*fUpwind[3]+0.25*alphaDrSurf[2]*fUpwind[2]+0.25*alphaDrSurf[1]*fUpwind[1]+0.25*alphaDrSurf[0]*fUpwind[0]; 
  drag_incr[1] = 0.25*alphaDrSurf[7]*fUpwind[11]+0.25*fUpwind[7]*alphaDrSurf[11]+0.25*alphaDrSurf[3]*fUpwind[6]+0.25*fUpwind[3]*alphaDrSurf[6]+0.25*alphaDrSurf[2]*fUpwind[5]+0.25*fUpwind[2]*alphaDrSurf[5]+0.25*alphaDrSurf[0]*fUpwind[1]+0.25*fUpwind[0]*alphaDrSurf[1]; 
  drag_incr[2] = 0.25*alphaDrSurf[6]*fUpwind[11]+0.25*fUpwind[6]*alphaDrSurf[11]+0.25*alphaDrSurf[3]*fUpwind[7]+0.25*fUpwind[3]*alphaDrSurf[7]+0.25*alphaDrSurf[1]*fUpwind[5]+0.25*fUpwind[1]*alphaDrSurf[5]+0.25*alphaDrSurf[0]*fUpwind[2]+0.25*fUpwind[0]*alphaDrSurf[2]; 
  drag_incr[3] = 0.25*alphaDrSurf[5]*fUpwind[11]+0.25*fUpwind[5]*alphaDrSurf[11]+0.25*alphaDrSurf[2]*fUpwind[7]+0.25*fUpwind[2]*alphaDrSurf[7]+0.25*alphaDrSurf[1]*fUpwind[6]+0.25*fUpwind[1]*alphaDrSurf[6]+0.25*alphaDrSurf[0]*fUpwind[3]+0.25*fUpwind[0]*alphaDrSurf[3]; 
  drag_incr[4] = 0.25*alphaDrSurf[11]*fUpwind[15]+0.25*alphaDrSurf[7]*fUpwind[14]+0.25*alphaDrSurf[6]*fUpwind[13]+0.25*alphaDrSurf[5]*fUpwind[12]+0.25*alphaDrSurf[3]*fUpwind[10]+0.25*alphaDrSurf[2]*fUpwind[9]+0.25*alphaDrSurf[1]*fUpwind[8]+0.25*alphaDrSurf[0]*fUpwind[4]; 
  drag_incr[5] = 0.25*alphaDrSurf[3]*fUpwind[11]+0.25*fUpwind[3]*alphaDrSurf[11]+0.25*alphaDrSurf[6]*fUpwind[7]+0.25*fUpwind[6]*alphaDrSurf[7]+0.25*alphaDrSurf[0]*fUpwind[5]+0.25*fUpwind[0]*alphaDrSurf[5]+0.25*alphaDrSurf[1]*fUpwind[2]+0.25*fUpwind[1]*alphaDrSurf[2]; 
  drag_incr[6] = 0.25*alphaDrSurf[2]*fUpwind[11]+0.25*fUpwind[2]*alphaDrSurf[11]+0.25*alphaDrSurf[5]*fUpwind[7]+0.25*fUpwind[5]*alphaDrSurf[7]+0.25*alphaDrSurf[0]*fUpwind[6]+0.25*fUpwind[0]*alphaDrSurf[6]+0.25*alphaDrSurf[1]*fUpwind[3]+0.25*fUpwind[1]*alphaDrSurf[3]; 
  drag_incr[7] = 0.25*alphaDrSurf[1]*fUpwind[11]+0.25*fUpwind[1]*alphaDrSurf[11]+0.25*alphaDrSurf[0]*fUpwind[7]+0.25*fUpwind[0]*alphaDrSurf[7]+0.25*alphaDrSurf[5]*fUpwind[6]+0.25*fUpwind[5]*alphaDrSurf[6]+0.25*alphaDrSurf[2]*fUpwind[3]+0.25*fUpwind[2]*alphaDrSurf[3]; 
  drag_incr[8] = 0.25*alphaDrSurf[7]*fUpwind[15]+0.25*alphaDrSurf[11]*fUpwind[14]+0.25*alphaDrSurf[3]*fUpwind[13]+0.25*alphaDrSurf[2]*fUpwind[12]+0.25*alphaDrSurf[6]*fUpwind[10]+0.25*alphaDrSurf[5]*fUpwind[9]+0.25*alphaDrSurf[0]*fUpwind[8]+0.25*alphaDrSurf[1]*fUpwind[4]; 
  drag_incr[9] = 0.25*alphaDrSurf[6]*fUpwind[15]+0.25*alphaDrSurf[3]*fUpwind[14]+0.25*alphaDrSurf[11]*fUpwind[13]+0.25*alphaDrSurf[1]*fUpwind[12]+0.25*alphaDrSurf[7]*fUpwind[10]+0.25*alphaDrSurf[0]*fUpwind[9]+0.25*alphaDrSurf[5]*fUpwind[8]+0.25*alphaDrSurf[2]*fUpwind[4]; 
  drag_incr[10] = 0.25*alphaDrSurf[5]*fUpwind[15]+0.25*alphaDrSurf[2]*fUpwind[14]+0.25*alphaDrSurf[1]*fUpwind[13]+0.25*alphaDrSurf[11]*fUpwind[12]+0.25*alphaDrSurf[0]*fUpwind[10]+0.25*alphaDrSurf[7]*fUpwind[9]+0.25*alphaDrSurf[6]*fUpwind[8]+0.25*alphaDrSurf[3]*fUpwind[4]; 
  drag_incr[11] = 0.25*alphaDrSurf[0]*fUpwind[11]+0.25*fUpwind[0]*alphaDrSurf[11]+0.25*alphaDrSurf[1]*fUpwind[7]+0.25*fUpwind[1]*alphaDrSurf[7]+0.25*alphaDrSurf[2]*fUpwind[6]+0.25*fUpwind[2]*alphaDrSurf[6]+0.25*alphaDrSurf[3]*fUpwind[5]+0.25*fUpwind[3]*alphaDrSurf[5]; 
  drag_incr[12] = 0.25*alphaDrSurf[3]*fUpwind[15]+0.25*alphaDrSurf[6]*fUpwind[14]+0.25*alphaDrSurf[7]*fUpwind[13]+0.25*alphaDrSurf[0]*fUpwind[12]+0.25*fUpwind[10]*alphaDrSurf[11]+0.25*alphaDrSurf[1]*fUpwind[9]+0.25*alphaDrSurf[2]*fUpwind[8]+0.25*fUpwind[4]*alphaDrSurf[5]; 
  drag_incr[13] = 0.25*alphaDrSurf[2]*fUpwind[15]+0.25*alphaDrSurf[5]*fUpwind[14]+0.25*alphaDrSurf[0]*fUpwind[13]+0.25*alphaDrSurf[7]*fUpwind[12]+0.25*fUpwind[9]*alphaDrSurf[11]+0.25*alphaDrSurf[1]*fUpwind[10]+0.25*alphaDrSurf[3]*fUpwind[8]+0.25*fUpwind[4]*alphaDrSurf[6]; 
  drag_incr[14] = 0.25*alphaDrSurf[1]*fUpwind[15]+0.25*alphaDrSurf[0]*fUpwind[14]+0.25*alphaDrSurf[5]*fUpwind[13]+0.25*alphaDrSurf[6]*fUpwind[12]+0.25*fUpwind[8]*alphaDrSurf[11]+0.25*alphaDrSurf[2]*fUpwind[10]+0.25*alphaDrSurf[3]*fUpwind[9]+0.25*fUpwind[4]*alphaDrSurf[7]; 
  drag_incr[15] = 0.25*alphaDrSurf[0]*fUpwind[15]+0.25*alphaDrSurf[1]*fUpwind[14]+0.25*alphaDrSurf[2]*fUpwind[13]+0.25*alphaDrSurf[3]*fUpwind[12]+0.25*fUpwind[4]*alphaDrSurf[11]+0.25*alphaDrSurf[5]*fUpwind[10]+0.25*alphaDrSurf[6]*fUpwind[9]+0.25*alphaDrSurf[7]*fUpwind[8]; 

  out[0] += -0.7071067811865475*drag_incr[0]*rdv2; 
  out[1] += -0.7071067811865475*drag_incr[1]*rdv2; 
  out[2] += -0.7071067811865475*drag_incr[2]*rdv2; 
  out[3] += -0.7071067811865475*drag_incr[3]*rdv2; 
  out[4] += 1.224744871391589*drag_incr[0]*rdv2; 
  out[5] += -0.7071067811865475*drag_incr[4]*rdv2; 
  out[6] += -0.7071067811865475*drag_incr[5]*rdv2; 
  out[7] += -0.7071067811865475*drag_incr[6]*rdv2; 
  out[8] += -0.7071067811865475*drag_incr[7]*rdv2; 
  out[9] += 1.224744871391589*drag_incr[1]*rdv2; 
  out[10] += 1.224744871391589*drag_incr[2]*rdv2; 
  out[11] += 1.224744871391589*drag_incr[3]*rdv2; 
  out[12] += -0.7071067811865475*drag_incr[8]*rdv2; 
  out[13] += -0.7071067811865475*drag_incr[9]*rdv2; 
  out[14] += -0.7071067811865475*drag_incr[10]*rdv2; 
  out[15] += 1.224744871391589*drag_incr[4]*rdv2; 
  out[16] += -0.7071067811865475*drag_incr[11]*rdv2; 
  out[17] += 1.224744871391589*drag_incr[5]*rdv2; 
  out[18] += 1.224744871391589*drag_incr[6]*rdv2; 
  out[19] += 1.224744871391589*drag_incr[7]*rdv2; 
  out[20] += -0.7071067811865475*drag_incr[12]*rdv2; 
  out[21] += -0.7071067811865475*drag_incr[13]*rdv2; 
  out[22] += -0.7071067811865475*drag_incr[14]*rdv2; 
  out[23] += 1.224744871391589*drag_incr[8]*rdv2; 
  out[24] += 1.224744871391589*drag_incr[9]*rdv2; 
  out[25] += 1.224744871391589*drag_incr[10]*rdv2; 
  out[26] += 1.224744871391589*drag_incr[11]*rdv2; 
  out[27] += -0.7071067811865475*drag_incr[15]*rdv2; 
  out[28] += 1.224744871391589*drag_incr[12]*rdv2; 
  out[29] += 1.224744871391589*drag_incr[13]*rdv2; 
  out[30] += 1.224744871391589*drag_incr[14]*rdv2; 
  out[31] += 1.224744871391589*drag_incr[15]*rdv2; 

  } 
} 
