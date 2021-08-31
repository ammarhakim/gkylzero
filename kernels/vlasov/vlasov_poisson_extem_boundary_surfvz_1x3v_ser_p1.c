#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_ser_1x3v_p1_surfvz_quad.h> 
GKYL_CU_DH void vlasov_poisson_extem_boundary_surfvz_1x3v_ser_p1(const double *w, const double *dxv, const double *fac_phi, const double *vecA, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w:           Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // fac_phi:     potential (scaled by appropriate factors).
  // vecA:        vector potential (scaled by appropriate factors). Unused in pure Vlasov-Poisson. 
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Output distribution function in skin cell 
  const double dv12 = 2/dxv[3]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv2 = dxv[2], wv2 = w[2]; 
  const double dv3 = dxv[3], wv3 = w[3]; 
  const double *phi = &fac_phi[0]; 
  const double dx10 = 2/dxv[0]; 
  const double *A0 = &vecA[0]; 
  const double *A1 = &vecA[2]; 
  const double *A2 = &vecA[4]; 
  double alpha[8] = {0.0}; 

  alpha[0] = -3.464101615137754*A2[1]*dx10*wv1; 
  alpha[2] = -1.0*A2[1]*dv1*dx10; 

  double fUpwindQuad[8] = {0.0};
  double fUpwind[8] = {0.0};
  double Ghat[8] = {0.0}; 

  if (edge == -1) { 

  if (alpha[0]-alpha[2] > 0) { 

    fUpwindQuad[0] = ser_1x3v_p1_surfvz_quad_0(1, fSkin); 
  } else { 

    fUpwindQuad[0] = ser_1x3v_p1_surfvz_quad_0(-1, fEdge); 
  } 
  if (alpha[0]-alpha[2] > 0) { 

    fUpwindQuad[1] = ser_1x3v_p1_surfvz_quad_1(1, fSkin); 
  } else { 

    fUpwindQuad[1] = ser_1x3v_p1_surfvz_quad_1(-1, fEdge); 
  } 
  if (alpha[2]+alpha[0] > 0) { 

    fUpwindQuad[2] = ser_1x3v_p1_surfvz_quad_2(1, fSkin); 
  } else { 

    fUpwindQuad[2] = ser_1x3v_p1_surfvz_quad_2(-1, fEdge); 
  } 
  if (alpha[2]+alpha[0] > 0) { 

    fUpwindQuad[3] = ser_1x3v_p1_surfvz_quad_3(1, fSkin); 
  } else { 

    fUpwindQuad[3] = ser_1x3v_p1_surfvz_quad_3(-1, fEdge); 
  } 
  if (alpha[0]-alpha[2] > 0) { 

    fUpwindQuad[4] = ser_1x3v_p1_surfvz_quad_4(1, fSkin); 
  } else { 

    fUpwindQuad[4] = ser_1x3v_p1_surfvz_quad_4(-1, fEdge); 
  } 
  if (alpha[0]-alpha[2] > 0) { 

    fUpwindQuad[5] = ser_1x3v_p1_surfvz_quad_5(1, fSkin); 
  } else { 

    fUpwindQuad[5] = ser_1x3v_p1_surfvz_quad_5(-1, fEdge); 
  } 
  if (alpha[2]+alpha[0] > 0) { 

    fUpwindQuad[6] = ser_1x3v_p1_surfvz_quad_6(1, fSkin); 
  } else { 

    fUpwindQuad[6] = ser_1x3v_p1_surfvz_quad_6(-1, fEdge); 
  } 
  if (alpha[2]+alpha[0] > 0) { 

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

  Ghat[0] += 0.3535533905932737*(alpha[2]*fUpwind[2]+alpha[0]*fUpwind[0]); 
  Ghat[1] += 0.3535533905932737*(alpha[2]*fUpwind[4]+alpha[0]*fUpwind[1]); 
  Ghat[2] += 0.3535533905932737*(alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2]); 
  Ghat[3] += 0.3535533905932737*(alpha[2]*fUpwind[6]+alpha[0]*fUpwind[3]); 
  Ghat[4] += 0.3535533905932737*(alpha[0]*fUpwind[4]+fUpwind[1]*alpha[2]); 
  Ghat[5] += 0.3535533905932737*(alpha[2]*fUpwind[7]+alpha[0]*fUpwind[5]); 
  Ghat[6] += 0.3535533905932737*(alpha[0]*fUpwind[6]+alpha[2]*fUpwind[3]); 
  Ghat[7] += 0.3535533905932737*(alpha[0]*fUpwind[7]+alpha[2]*fUpwind[5]); 

  out[0] += -0.7071067811865475*Ghat[0]*dv12; 
  out[1] += -0.7071067811865475*Ghat[1]*dv12; 
  out[2] += -0.7071067811865475*Ghat[2]*dv12; 
  out[3] += -0.7071067811865475*Ghat[3]*dv12; 
  out[4] += -1.224744871391589*Ghat[0]*dv12; 
  out[5] += -0.7071067811865475*Ghat[4]*dv12; 
  out[6] += -0.7071067811865475*Ghat[5]*dv12; 
  out[7] += -0.7071067811865475*Ghat[6]*dv12; 
  out[8] += -1.224744871391589*Ghat[1]*dv12; 
  out[9] += -1.224744871391589*Ghat[2]*dv12; 
  out[10] += -1.224744871391589*Ghat[3]*dv12; 
  out[11] += -0.7071067811865475*Ghat[7]*dv12; 
  out[12] += -1.224744871391589*Ghat[4]*dv12; 
  out[13] += -1.224744871391589*Ghat[5]*dv12; 
  out[14] += -1.224744871391589*Ghat[6]*dv12; 
  out[15] += -1.224744871391589*Ghat[7]*dv12; 

  } else { 

  if (alpha[0]-alpha[2] > 0) { 

    fUpwindQuad[0] = ser_1x3v_p1_surfvz_quad_0(1, fEdge); 
  } else { 

    fUpwindQuad[0] = ser_1x3v_p1_surfvz_quad_0(-1, fSkin); 
  } 
  if (alpha[0]-alpha[2] > 0) { 

    fUpwindQuad[1] = ser_1x3v_p1_surfvz_quad_1(1, fEdge); 
  } else { 

    fUpwindQuad[1] = ser_1x3v_p1_surfvz_quad_1(-1, fSkin); 
  } 
  if (alpha[2]+alpha[0] > 0) { 

    fUpwindQuad[2] = ser_1x3v_p1_surfvz_quad_2(1, fEdge); 
  } else { 

    fUpwindQuad[2] = ser_1x3v_p1_surfvz_quad_2(-1, fSkin); 
  } 
  if (alpha[2]+alpha[0] > 0) { 

    fUpwindQuad[3] = ser_1x3v_p1_surfvz_quad_3(1, fEdge); 
  } else { 

    fUpwindQuad[3] = ser_1x3v_p1_surfvz_quad_3(-1, fSkin); 
  } 
  if (alpha[0]-alpha[2] > 0) { 

    fUpwindQuad[4] = ser_1x3v_p1_surfvz_quad_4(1, fEdge); 
  } else { 

    fUpwindQuad[4] = ser_1x3v_p1_surfvz_quad_4(-1, fSkin); 
  } 
  if (alpha[0]-alpha[2] > 0) { 

    fUpwindQuad[5] = ser_1x3v_p1_surfvz_quad_5(1, fEdge); 
  } else { 

    fUpwindQuad[5] = ser_1x3v_p1_surfvz_quad_5(-1, fSkin); 
  } 
  if (alpha[2]+alpha[0] > 0) { 

    fUpwindQuad[6] = ser_1x3v_p1_surfvz_quad_6(1, fEdge); 
  } else { 

    fUpwindQuad[6] = ser_1x3v_p1_surfvz_quad_6(-1, fSkin); 
  } 
  if (alpha[2]+alpha[0] > 0) { 

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

  Ghat[0] += 0.3535533905932737*(alpha[2]*fUpwind[2]+alpha[0]*fUpwind[0]); 
  Ghat[1] += 0.3535533905932737*(alpha[2]*fUpwind[4]+alpha[0]*fUpwind[1]); 
  Ghat[2] += 0.3535533905932737*(alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2]); 
  Ghat[3] += 0.3535533905932737*(alpha[2]*fUpwind[6]+alpha[0]*fUpwind[3]); 
  Ghat[4] += 0.3535533905932737*(alpha[0]*fUpwind[4]+fUpwind[1]*alpha[2]); 
  Ghat[5] += 0.3535533905932737*(alpha[2]*fUpwind[7]+alpha[0]*fUpwind[5]); 
  Ghat[6] += 0.3535533905932737*(alpha[0]*fUpwind[6]+alpha[2]*fUpwind[3]); 
  Ghat[7] += 0.3535533905932737*(alpha[0]*fUpwind[7]+alpha[2]*fUpwind[5]); 

  out[0] += 0.7071067811865475*Ghat[0]*dv12; 
  out[1] += 0.7071067811865475*Ghat[1]*dv12; 
  out[2] += 0.7071067811865475*Ghat[2]*dv12; 
  out[3] += 0.7071067811865475*Ghat[3]*dv12; 
  out[4] += -1.224744871391589*Ghat[0]*dv12; 
  out[5] += 0.7071067811865475*Ghat[4]*dv12; 
  out[6] += 0.7071067811865475*Ghat[5]*dv12; 
  out[7] += 0.7071067811865475*Ghat[6]*dv12; 
  out[8] += -1.224744871391589*Ghat[1]*dv12; 
  out[9] += -1.224744871391589*Ghat[2]*dv12; 
  out[10] += -1.224744871391589*Ghat[3]*dv12; 
  out[11] += 0.7071067811865475*Ghat[7]*dv12; 
  out[12] += -1.224744871391589*Ghat[4]*dv12; 
  out[13] += -1.224744871391589*Ghat[5]*dv12; 
  out[14] += -1.224744871391589*Ghat[6]*dv12; 
  out[15] += -1.224744871391589*Ghat[7]*dv12; 

  } 
} 
