#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_ser_1x3v_p1_surfvz_quad.h> 
GKYL_CU_DH void vlasov_boundary_surfvz_1x3v_ser_p1(const double *w, const double *dxv, const double *qmem, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w:           Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // qmem:        q/m*EM fields.
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Output distribution function in skin cell 
  const double dv12 = 2/dxv[3]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv2 = dxv[2], wv2 = w[2]; 
  const double dv3 = dxv[3], wv3 = w[3]; 
  const double *E2 = &qmem[4]; 
  const double *B0 = &qmem[6]; 
  const double *B1 = &qmem[8]; 
  const double *B2 = &qmem[10]; 

  double alpha[8] = {0.0}; 

  alpha[0] = 2.0*(B1[0]*wv1+E2[0])-2.0*B0[0]*wv2; 
  alpha[1] = 2.0*(B1[1]*wv1+E2[1])-2.0*B0[1]*wv2; 
  alpha[2] = 0.5773502691896258*B1[0]*dv1; 
  alpha[3] = -0.5773502691896258*B0[0]*dv2; 
  alpha[4] = 0.5773502691896258*B1[1]*dv1; 
  alpha[5] = -0.5773502691896258*B0[1]*dv2; 

  double fUpwindQuad[8] = {0.0};
  double fUpwind[8] = {0.0};
  double Ghat[8] = {0.0}; 

  if (edge == -1) { 

  if (alpha[5]+alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[0] = ser_1x3v_p1_surfvz_quad_0(1, fSkin); 
  } else { 

    fUpwindQuad[0] = ser_1x3v_p1_surfvz_quad_0(-1, fEdge); 
  } 
  if ((-alpha[5])-alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[1] = ser_1x3v_p1_surfvz_quad_1(1, fSkin); 
  } else { 

    fUpwindQuad[1] = ser_1x3v_p1_surfvz_quad_1(-1, fEdge); 
  } 
  if (alpha[5]-alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[2] = ser_1x3v_p1_surfvz_quad_2(1, fSkin); 
  } else { 

    fUpwindQuad[2] = ser_1x3v_p1_surfvz_quad_2(-1, fEdge); 
  } 
  if ((-alpha[5])+alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[3] = ser_1x3v_p1_surfvz_quad_3(1, fSkin); 
  } else { 

    fUpwindQuad[3] = ser_1x3v_p1_surfvz_quad_3(-1, fEdge); 
  } 
  if ((-alpha[5])+alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[4] = ser_1x3v_p1_surfvz_quad_4(1, fSkin); 
  } else { 

    fUpwindQuad[4] = ser_1x3v_p1_surfvz_quad_4(-1, fEdge); 
  } 
  if (alpha[5]-alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[5] = ser_1x3v_p1_surfvz_quad_5(1, fSkin); 
  } else { 

    fUpwindQuad[5] = ser_1x3v_p1_surfvz_quad_5(-1, fEdge); 
  } 
  if ((-alpha[5])-alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[6] = ser_1x3v_p1_surfvz_quad_6(1, fSkin); 
  } else { 

    fUpwindQuad[6] = ser_1x3v_p1_surfvz_quad_6(-1, fEdge); 
  } 
  if (alpha[5]+alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

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

  Ghat[0] += 0.3535533905932737*(alpha[5]*fUpwind[5]+alpha[4]*fUpwind[4]+alpha[3]*fUpwind[3]+alpha[2]*fUpwind[2]+alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] += 0.3535533905932737*(alpha[3]*fUpwind[5]+fUpwind[3]*alpha[5]+alpha[2]*fUpwind[4]+fUpwind[2]*alpha[4]+alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]); 
  Ghat[2] += 0.3535533905932737*(alpha[5]*fUpwind[7]+alpha[3]*fUpwind[6]+alpha[1]*fUpwind[4]+fUpwind[1]*alpha[4]+alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2]); 
  Ghat[3] += 0.3535533905932737*(alpha[4]*fUpwind[7]+alpha[2]*fUpwind[6]+alpha[1]*fUpwind[5]+fUpwind[1]*alpha[5]+alpha[0]*fUpwind[3]+fUpwind[0]*alpha[3]); 
  Ghat[4] += 0.3535533905932737*(alpha[3]*fUpwind[7]+alpha[5]*fUpwind[6]+alpha[0]*fUpwind[4]+fUpwind[0]*alpha[4]+alpha[1]*fUpwind[2]+fUpwind[1]*alpha[2]); 
  Ghat[5] += 0.3535533905932737*(alpha[2]*fUpwind[7]+alpha[4]*fUpwind[6]+alpha[0]*fUpwind[5]+fUpwind[0]*alpha[5]+alpha[1]*fUpwind[3]+fUpwind[1]*alpha[3]); 
  Ghat[6] += 0.3535533905932737*(alpha[1]*fUpwind[7]+alpha[0]*fUpwind[6]+alpha[4]*fUpwind[5]+fUpwind[4]*alpha[5]+alpha[2]*fUpwind[3]+fUpwind[2]*alpha[3]); 
  Ghat[7] += 0.3535533905932737*(alpha[0]*fUpwind[7]+alpha[1]*fUpwind[6]+alpha[2]*fUpwind[5]+fUpwind[2]*alpha[5]+alpha[3]*fUpwind[4]+fUpwind[3]*alpha[4]); 

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

  if (alpha[5]+alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[0] = ser_1x3v_p1_surfvz_quad_0(1, fEdge); 
  } else { 

    fUpwindQuad[0] = ser_1x3v_p1_surfvz_quad_0(-1, fSkin); 
  } 
  if ((-alpha[5])-alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[1] = ser_1x3v_p1_surfvz_quad_1(1, fEdge); 
  } else { 

    fUpwindQuad[1] = ser_1x3v_p1_surfvz_quad_1(-1, fSkin); 
  } 
  if (alpha[5]-alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[2] = ser_1x3v_p1_surfvz_quad_2(1, fEdge); 
  } else { 

    fUpwindQuad[2] = ser_1x3v_p1_surfvz_quad_2(-1, fSkin); 
  } 
  if ((-alpha[5])+alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[3] = ser_1x3v_p1_surfvz_quad_3(1, fEdge); 
  } else { 

    fUpwindQuad[3] = ser_1x3v_p1_surfvz_quad_3(-1, fSkin); 
  } 
  if ((-alpha[5])+alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[4] = ser_1x3v_p1_surfvz_quad_4(1, fEdge); 
  } else { 

    fUpwindQuad[4] = ser_1x3v_p1_surfvz_quad_4(-1, fSkin); 
  } 
  if (alpha[5]-alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[5] = ser_1x3v_p1_surfvz_quad_5(1, fEdge); 
  } else { 

    fUpwindQuad[5] = ser_1x3v_p1_surfvz_quad_5(-1, fSkin); 
  } 
  if ((-alpha[5])-alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[6] = ser_1x3v_p1_surfvz_quad_6(1, fEdge); 
  } else { 

    fUpwindQuad[6] = ser_1x3v_p1_surfvz_quad_6(-1, fSkin); 
  } 
  if (alpha[5]+alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

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

  Ghat[0] += 0.3535533905932737*(alpha[5]*fUpwind[5]+alpha[4]*fUpwind[4]+alpha[3]*fUpwind[3]+alpha[2]*fUpwind[2]+alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] += 0.3535533905932737*(alpha[3]*fUpwind[5]+fUpwind[3]*alpha[5]+alpha[2]*fUpwind[4]+fUpwind[2]*alpha[4]+alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]); 
  Ghat[2] += 0.3535533905932737*(alpha[5]*fUpwind[7]+alpha[3]*fUpwind[6]+alpha[1]*fUpwind[4]+fUpwind[1]*alpha[4]+alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2]); 
  Ghat[3] += 0.3535533905932737*(alpha[4]*fUpwind[7]+alpha[2]*fUpwind[6]+alpha[1]*fUpwind[5]+fUpwind[1]*alpha[5]+alpha[0]*fUpwind[3]+fUpwind[0]*alpha[3]); 
  Ghat[4] += 0.3535533905932737*(alpha[3]*fUpwind[7]+alpha[5]*fUpwind[6]+alpha[0]*fUpwind[4]+fUpwind[0]*alpha[4]+alpha[1]*fUpwind[2]+fUpwind[1]*alpha[2]); 
  Ghat[5] += 0.3535533905932737*(alpha[2]*fUpwind[7]+alpha[4]*fUpwind[6]+alpha[0]*fUpwind[5]+fUpwind[0]*alpha[5]+alpha[1]*fUpwind[3]+fUpwind[1]*alpha[3]); 
  Ghat[6] += 0.3535533905932737*(alpha[1]*fUpwind[7]+alpha[0]*fUpwind[6]+alpha[4]*fUpwind[5]+fUpwind[4]*alpha[5]+alpha[2]*fUpwind[3]+fUpwind[2]*alpha[3]); 
  Ghat[7] += 0.3535533905932737*(alpha[0]*fUpwind[7]+alpha[1]*fUpwind[6]+alpha[2]*fUpwind[5]+fUpwind[2]*alpha[5]+alpha[3]*fUpwind[4]+fUpwind[3]*alpha[4]); 

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
