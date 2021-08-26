#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_ser_1x2v_p2_surfvx_quad.h> 
GKYL_CU_DH void vlasov_boundary_surfvx_1x2v_ser_p2(const double *w, const double *dxv, const double *qmem, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w:           Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // qmem:        q/m*EM fields.
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Output distribution function in skin cell 
  const double dv10 = 2/dxv[1]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv2 = dxv[2], wv2 = w[2]; 
  const double *E0 = &qmem[0]; 
  const double *B2 = &qmem[15]; 

  double alpha[8] = {0.0}; 

  alpha[0] = 1.414213562373095*(B2[0]*wv2+E0[0]); 
  alpha[1] = 1.414213562373095*(B2[1]*wv2+E0[1]); 
  alpha[2] = 0.408248290463863*B2[0]*dv2; 
  alpha[3] = 0.408248290463863*B2[1]*dv2; 
  alpha[4] = 1.414213562373095*(B2[2]*wv2+E0[2]); 
  alpha[6] = 0.408248290463863*B2[2]*dv2; 

  double fUpwindQuad[9] = {0.0};
  double fUpwind[8] = {0.0};
  double Ghat[8] = {0.0}; 

  if (edge == -1) { 

  if ((-0.5999999999999999*alpha[6])+0.4472135954999579*alpha[4]+0.9*alpha[3]-0.6708203932499369*(alpha[2]+alpha[1])+0.5*alpha[0] > 0) { 

    fUpwindQuad[0] = ser_1x2v_p2_surfvx_quad_0(1, fSkin); 
  } else { 

    fUpwindQuad[0] = ser_1x2v_p2_surfvx_quad_0(-1, fEdge); 
  } 
  if (0.75*alpha[6]-0.5590169943749475*alpha[4]-0.6708203932499369*alpha[2]+0.5*alpha[0] > 0) { 

    fUpwindQuad[1] = ser_1x2v_p2_surfvx_quad_1(1, fSkin); 
  } else { 

    fUpwindQuad[1] = ser_1x2v_p2_surfvx_quad_1(-1, fEdge); 
  } 
  if ((-0.5999999999999999*alpha[6])+0.4472135954999579*alpha[4]-0.9*alpha[3]-0.6708203932499369*alpha[2]+0.6708203932499369*alpha[1]+0.5*alpha[0] > 0) { 

    fUpwindQuad[2] = ser_1x2v_p2_surfvx_quad_2(1, fSkin); 
  } else { 

    fUpwindQuad[2] = ser_1x2v_p2_surfvx_quad_2(-1, fEdge); 
  } 
  if (0.4472135954999579*alpha[4]-0.6708203932499369*alpha[1]+0.5*alpha[0] > 0) { 

    fUpwindQuad[3] = ser_1x2v_p2_surfvx_quad_3(1, fSkin); 
  } else { 

    fUpwindQuad[3] = ser_1x2v_p2_surfvx_quad_3(-1, fEdge); 
  } 
  if (0.5*alpha[0]-0.5590169943749475*alpha[4] > 0) { 

    fUpwindQuad[4] = ser_1x2v_p2_surfvx_quad_4(1, fSkin); 
  } else { 

    fUpwindQuad[4] = ser_1x2v_p2_surfvx_quad_4(-1, fEdge); 
  } 
  if (0.4472135954999579*alpha[4]+0.6708203932499369*alpha[1]+0.5*alpha[0] > 0) { 

    fUpwindQuad[5] = ser_1x2v_p2_surfvx_quad_5(1, fSkin); 
  } else { 

    fUpwindQuad[5] = ser_1x2v_p2_surfvx_quad_5(-1, fEdge); 
  } 
  if (0.5999999999999999*alpha[6]+0.4472135954999579*alpha[4]-0.9*alpha[3]+0.6708203932499369*alpha[2]-0.6708203932499369*alpha[1]+0.5*alpha[0] > 0) { 

    fUpwindQuad[6] = ser_1x2v_p2_surfvx_quad_6(1, fSkin); 
  } else { 

    fUpwindQuad[6] = ser_1x2v_p2_surfvx_quad_6(-1, fEdge); 
  } 
  if ((-0.75*alpha[6])-0.5590169943749475*alpha[4]+0.6708203932499369*alpha[2]+0.5*alpha[0] > 0) { 

    fUpwindQuad[7] = ser_1x2v_p2_surfvx_quad_7(1, fSkin); 
  } else { 

    fUpwindQuad[7] = ser_1x2v_p2_surfvx_quad_7(-1, fEdge); 
  } 
  if (0.5999999999999999*alpha[6]+0.4472135954999579*alpha[4]+0.9*alpha[3]+0.6708203932499369*(alpha[2]+alpha[1])+0.5*alpha[0] > 0) { 

    fUpwindQuad[8] = ser_1x2v_p2_surfvx_quad_8(1, fSkin); 
  } else { 

    fUpwindQuad[8] = ser_1x2v_p2_surfvx_quad_8(-1, fEdge); 
  } 

  fUpwind[0] = 0.154320987654321*fUpwindQuad[8]+0.2469135802469136*fUpwindQuad[7]+0.154320987654321*fUpwindQuad[6]+0.2469135802469136*fUpwindQuad[5]+0.3950617283950617*fUpwindQuad[4]+0.2469135802469136*fUpwindQuad[3]+0.154320987654321*fUpwindQuad[2]+0.2469135802469136*fUpwindQuad[1]+0.154320987654321*fUpwindQuad[0]; 
  fUpwind[1] = 0.2070433312499806*fUpwindQuad[8]-0.2070433312499806*fUpwindQuad[6]+0.3312693299999688*fUpwindQuad[5]-0.3312693299999688*fUpwindQuad[3]+0.2070433312499806*fUpwindQuad[2]-0.2070433312499806*fUpwindQuad[0]; 
  fUpwind[2] = 0.2070433312499806*fUpwindQuad[8]+0.3312693299999688*fUpwindQuad[7]+0.2070433312499806*fUpwindQuad[6]-0.2070433312499806*fUpwindQuad[2]-0.3312693299999688*fUpwindQuad[1]-0.2070433312499806*fUpwindQuad[0]; 
  fUpwind[3] = 0.2777777777777778*fUpwindQuad[8]-0.2777777777777778*fUpwindQuad[6]-0.2777777777777778*fUpwindQuad[2]+0.2777777777777778*fUpwindQuad[0]; 
  fUpwind[4] = 0.138028887499987*fUpwindQuad[8]-0.2760577749999741*fUpwindQuad[7]+0.138028887499987*fUpwindQuad[6]+0.2208462199999792*fUpwindQuad[5]-0.4416924399999584*fUpwindQuad[4]+0.2208462199999792*fUpwindQuad[3]+0.138028887499987*fUpwindQuad[2]-0.2760577749999741*fUpwindQuad[1]+0.138028887499987*fUpwindQuad[0]; 
  fUpwind[5] = 0.138028887499987*fUpwindQuad[8]+0.2208462199999792*fUpwindQuad[7]+0.138028887499987*fUpwindQuad[6]-0.2760577749999741*fUpwindQuad[5]-0.4416924399999584*fUpwindQuad[4]-0.2760577749999741*fUpwindQuad[3]+0.138028887499987*fUpwindQuad[2]+0.2208462199999792*fUpwindQuad[1]+0.138028887499987*fUpwindQuad[0]; 
  fUpwind[6] = 0.1851851851851853*fUpwindQuad[8]-0.3703703703703705*fUpwindQuad[7]+0.1851851851851853*fUpwindQuad[6]-0.1851851851851853*fUpwindQuad[2]+0.3703703703703705*fUpwindQuad[1]-0.1851851851851853*fUpwindQuad[0]; 
  fUpwind[7] = 0.1851851851851853*fUpwindQuad[8]-0.1851851851851853*fUpwindQuad[6]-0.3703703703703705*fUpwindQuad[5]+0.3703703703703705*fUpwindQuad[3]+0.1851851851851853*fUpwindQuad[2]-0.1851851851851853*fUpwindQuad[0]; 

  Ghat[0] += 0.5*(alpha[6]*fUpwind[6]+alpha[4]*fUpwind[4]+alpha[3]*fUpwind[3]+alpha[2]*fUpwind[2]+alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] += 0.447213595499958*(alpha[3]*fUpwind[6]+fUpwind[3]*alpha[6])+0.4472135954999579*(alpha[1]*fUpwind[4]+fUpwind[1]*alpha[4])+0.5*(alpha[2]*fUpwind[3]+fUpwind[2]*alpha[3]+alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]); 
  Ghat[2] += 0.447213595499958*alpha[3]*fUpwind[7]+0.5000000000000001*(alpha[4]*fUpwind[6]+fUpwind[4]*alpha[6])+0.4472135954999579*alpha[2]*fUpwind[5]+0.5*(alpha[1]*fUpwind[3]+fUpwind[1]*alpha[3]+alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2]); 
  Ghat[3] += 0.4*alpha[6]*fUpwind[7]+0.447213595499958*(alpha[2]*fUpwind[7]+alpha[1]*fUpwind[6]+fUpwind[1]*alpha[6])+0.4472135954999579*(alpha[3]*(fUpwind[5]+fUpwind[4])+fUpwind[3]*alpha[4])+0.5*(alpha[0]*fUpwind[3]+fUpwind[0]*alpha[3]+alpha[1]*fUpwind[2]+fUpwind[1]*alpha[2]); 
  Ghat[4] += 0.31943828249997*alpha[6]*fUpwind[6]+0.5000000000000001*(alpha[2]*fUpwind[6]+fUpwind[2]*alpha[6])+0.31943828249997*alpha[4]*fUpwind[4]+0.5*(alpha[0]*fUpwind[4]+fUpwind[0]*alpha[4])+0.4472135954999579*(alpha[3]*fUpwind[3]+alpha[1]*fUpwind[1]); 
  Ghat[5] += 0.5000000000000001*alpha[1]*fUpwind[7]+0.4472135954999579*alpha[6]*fUpwind[6]+0.5*alpha[0]*fUpwind[5]+0.4472135954999579*(alpha[3]*fUpwind[3]+alpha[2]*fUpwind[2]); 
  Ghat[6] += 0.4*alpha[3]*fUpwind[7]+(0.31943828249997*alpha[4]+0.5*alpha[0])*fUpwind[6]+(0.4472135954999579*fUpwind[5]+0.31943828249997*fUpwind[4]+0.5*fUpwind[0])*alpha[6]+0.5000000000000001*(alpha[2]*fUpwind[4]+fUpwind[2]*alpha[4])+0.447213595499958*(alpha[1]*fUpwind[3]+fUpwind[1]*alpha[3]); 
  Ghat[7] += (0.4472135954999579*alpha[4]+0.5*alpha[0])*fUpwind[7]+0.4*(alpha[3]*fUpwind[6]+fUpwind[3]*alpha[6])+0.5000000000000001*alpha[1]*fUpwind[5]+0.447213595499958*(alpha[2]*fUpwind[3]+fUpwind[2]*alpha[3]); 

  out[0] += -0.7071067811865475*Ghat[0]*dv10; 
  out[1] += -0.7071067811865475*Ghat[1]*dv10; 
  out[2] += -1.224744871391589*Ghat[0]*dv10; 
  out[3] += -0.7071067811865475*Ghat[2]*dv10; 
  out[4] += -1.224744871391589*Ghat[1]*dv10; 
  out[5] += -0.7071067811865475*Ghat[3]*dv10; 
  out[6] += -1.224744871391589*Ghat[2]*dv10; 
  out[7] += -0.7071067811865475*Ghat[4]*dv10; 
  out[8] += -1.58113883008419*Ghat[0]*dv10; 
  out[9] += -0.7071067811865475*Ghat[5]*dv10; 
  out[10] += -1.224744871391589*Ghat[3]*dv10; 
  out[11] += -1.224744871391589*Ghat[4]*dv10; 
  out[12] += -1.58113883008419*Ghat[1]*dv10; 
  out[13] += -0.7071067811865475*Ghat[6]*dv10; 
  out[14] += -1.58113883008419*Ghat[2]*dv10; 
  out[15] += -0.7071067811865475*Ghat[7]*dv10; 
  out[16] += -1.224744871391589*Ghat[5]*dv10; 
  out[17] += -1.224744871391589*Ghat[6]*dv10; 
  out[18] += -1.58113883008419*Ghat[3]*dv10; 
  out[19] += -1.224744871391589*Ghat[7]*dv10; 

  } else { 

  if ((-0.5999999999999999*alpha[6])+0.4472135954999579*alpha[4]+0.9*alpha[3]-0.6708203932499369*(alpha[2]+alpha[1])+0.5*alpha[0] > 0) { 

    fUpwindQuad[0] = ser_1x2v_p2_surfvx_quad_0(1, fEdge); 
  } else { 

    fUpwindQuad[0] = ser_1x2v_p2_surfvx_quad_0(-1, fSkin); 
  } 
  if (0.75*alpha[6]-0.5590169943749475*alpha[4]-0.6708203932499369*alpha[2]+0.5*alpha[0] > 0) { 

    fUpwindQuad[1] = ser_1x2v_p2_surfvx_quad_1(1, fEdge); 
  } else { 

    fUpwindQuad[1] = ser_1x2v_p2_surfvx_quad_1(-1, fSkin); 
  } 
  if ((-0.5999999999999999*alpha[6])+0.4472135954999579*alpha[4]-0.9*alpha[3]-0.6708203932499369*alpha[2]+0.6708203932499369*alpha[1]+0.5*alpha[0] > 0) { 

    fUpwindQuad[2] = ser_1x2v_p2_surfvx_quad_2(1, fEdge); 
  } else { 

    fUpwindQuad[2] = ser_1x2v_p2_surfvx_quad_2(-1, fSkin); 
  } 
  if (0.4472135954999579*alpha[4]-0.6708203932499369*alpha[1]+0.5*alpha[0] > 0) { 

    fUpwindQuad[3] = ser_1x2v_p2_surfvx_quad_3(1, fEdge); 
  } else { 

    fUpwindQuad[3] = ser_1x2v_p2_surfvx_quad_3(-1, fSkin); 
  } 
  if (0.5*alpha[0]-0.5590169943749475*alpha[4] > 0) { 

    fUpwindQuad[4] = ser_1x2v_p2_surfvx_quad_4(1, fEdge); 
  } else { 

    fUpwindQuad[4] = ser_1x2v_p2_surfvx_quad_4(-1, fSkin); 
  } 
  if (0.4472135954999579*alpha[4]+0.6708203932499369*alpha[1]+0.5*alpha[0] > 0) { 

    fUpwindQuad[5] = ser_1x2v_p2_surfvx_quad_5(1, fEdge); 
  } else { 

    fUpwindQuad[5] = ser_1x2v_p2_surfvx_quad_5(-1, fSkin); 
  } 
  if (0.5999999999999999*alpha[6]+0.4472135954999579*alpha[4]-0.9*alpha[3]+0.6708203932499369*alpha[2]-0.6708203932499369*alpha[1]+0.5*alpha[0] > 0) { 

    fUpwindQuad[6] = ser_1x2v_p2_surfvx_quad_6(1, fEdge); 
  } else { 

    fUpwindQuad[6] = ser_1x2v_p2_surfvx_quad_6(-1, fSkin); 
  } 
  if ((-0.75*alpha[6])-0.5590169943749475*alpha[4]+0.6708203932499369*alpha[2]+0.5*alpha[0] > 0) { 

    fUpwindQuad[7] = ser_1x2v_p2_surfvx_quad_7(1, fEdge); 
  } else { 

    fUpwindQuad[7] = ser_1x2v_p2_surfvx_quad_7(-1, fSkin); 
  } 
  if (0.5999999999999999*alpha[6]+0.4472135954999579*alpha[4]+0.9*alpha[3]+0.6708203932499369*(alpha[2]+alpha[1])+0.5*alpha[0] > 0) { 

    fUpwindQuad[8] = ser_1x2v_p2_surfvx_quad_8(1, fEdge); 
  } else { 

    fUpwindQuad[8] = ser_1x2v_p2_surfvx_quad_8(-1, fSkin); 
  } 

  fUpwind[0] = 0.154320987654321*fUpwindQuad[8]+0.2469135802469136*fUpwindQuad[7]+0.154320987654321*fUpwindQuad[6]+0.2469135802469136*fUpwindQuad[5]+0.3950617283950617*fUpwindQuad[4]+0.2469135802469136*fUpwindQuad[3]+0.154320987654321*fUpwindQuad[2]+0.2469135802469136*fUpwindQuad[1]+0.154320987654321*fUpwindQuad[0]; 
  fUpwind[1] = 0.2070433312499806*fUpwindQuad[8]-0.2070433312499806*fUpwindQuad[6]+0.3312693299999688*fUpwindQuad[5]-0.3312693299999688*fUpwindQuad[3]+0.2070433312499806*fUpwindQuad[2]-0.2070433312499806*fUpwindQuad[0]; 
  fUpwind[2] = 0.2070433312499806*fUpwindQuad[8]+0.3312693299999688*fUpwindQuad[7]+0.2070433312499806*fUpwindQuad[6]-0.2070433312499806*fUpwindQuad[2]-0.3312693299999688*fUpwindQuad[1]-0.2070433312499806*fUpwindQuad[0]; 
  fUpwind[3] = 0.2777777777777778*fUpwindQuad[8]-0.2777777777777778*fUpwindQuad[6]-0.2777777777777778*fUpwindQuad[2]+0.2777777777777778*fUpwindQuad[0]; 
  fUpwind[4] = 0.138028887499987*fUpwindQuad[8]-0.2760577749999741*fUpwindQuad[7]+0.138028887499987*fUpwindQuad[6]+0.2208462199999792*fUpwindQuad[5]-0.4416924399999584*fUpwindQuad[4]+0.2208462199999792*fUpwindQuad[3]+0.138028887499987*fUpwindQuad[2]-0.2760577749999741*fUpwindQuad[1]+0.138028887499987*fUpwindQuad[0]; 
  fUpwind[5] = 0.138028887499987*fUpwindQuad[8]+0.2208462199999792*fUpwindQuad[7]+0.138028887499987*fUpwindQuad[6]-0.2760577749999741*fUpwindQuad[5]-0.4416924399999584*fUpwindQuad[4]-0.2760577749999741*fUpwindQuad[3]+0.138028887499987*fUpwindQuad[2]+0.2208462199999792*fUpwindQuad[1]+0.138028887499987*fUpwindQuad[0]; 
  fUpwind[6] = 0.1851851851851853*fUpwindQuad[8]-0.3703703703703705*fUpwindQuad[7]+0.1851851851851853*fUpwindQuad[6]-0.1851851851851853*fUpwindQuad[2]+0.3703703703703705*fUpwindQuad[1]-0.1851851851851853*fUpwindQuad[0]; 
  fUpwind[7] = 0.1851851851851853*fUpwindQuad[8]-0.1851851851851853*fUpwindQuad[6]-0.3703703703703705*fUpwindQuad[5]+0.3703703703703705*fUpwindQuad[3]+0.1851851851851853*fUpwindQuad[2]-0.1851851851851853*fUpwindQuad[0]; 

  Ghat[0] += 0.5*(alpha[6]*fUpwind[6]+alpha[4]*fUpwind[4]+alpha[3]*fUpwind[3]+alpha[2]*fUpwind[2]+alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] += 0.447213595499958*(alpha[3]*fUpwind[6]+fUpwind[3]*alpha[6])+0.4472135954999579*(alpha[1]*fUpwind[4]+fUpwind[1]*alpha[4])+0.5*(alpha[2]*fUpwind[3]+fUpwind[2]*alpha[3]+alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]); 
  Ghat[2] += 0.447213595499958*alpha[3]*fUpwind[7]+0.5000000000000001*(alpha[4]*fUpwind[6]+fUpwind[4]*alpha[6])+0.4472135954999579*alpha[2]*fUpwind[5]+0.5*(alpha[1]*fUpwind[3]+fUpwind[1]*alpha[3]+alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2]); 
  Ghat[3] += 0.4*alpha[6]*fUpwind[7]+0.447213595499958*(alpha[2]*fUpwind[7]+alpha[1]*fUpwind[6]+fUpwind[1]*alpha[6])+0.4472135954999579*(alpha[3]*(fUpwind[5]+fUpwind[4])+fUpwind[3]*alpha[4])+0.5*(alpha[0]*fUpwind[3]+fUpwind[0]*alpha[3]+alpha[1]*fUpwind[2]+fUpwind[1]*alpha[2]); 
  Ghat[4] += 0.31943828249997*alpha[6]*fUpwind[6]+0.5000000000000001*(alpha[2]*fUpwind[6]+fUpwind[2]*alpha[6])+0.31943828249997*alpha[4]*fUpwind[4]+0.5*(alpha[0]*fUpwind[4]+fUpwind[0]*alpha[4])+0.4472135954999579*(alpha[3]*fUpwind[3]+alpha[1]*fUpwind[1]); 
  Ghat[5] += 0.5000000000000001*alpha[1]*fUpwind[7]+0.4472135954999579*alpha[6]*fUpwind[6]+0.5*alpha[0]*fUpwind[5]+0.4472135954999579*(alpha[3]*fUpwind[3]+alpha[2]*fUpwind[2]); 
  Ghat[6] += 0.4*alpha[3]*fUpwind[7]+(0.31943828249997*alpha[4]+0.5*alpha[0])*fUpwind[6]+(0.4472135954999579*fUpwind[5]+0.31943828249997*fUpwind[4]+0.5*fUpwind[0])*alpha[6]+0.5000000000000001*(alpha[2]*fUpwind[4]+fUpwind[2]*alpha[4])+0.447213595499958*(alpha[1]*fUpwind[3]+fUpwind[1]*alpha[3]); 
  Ghat[7] += (0.4472135954999579*alpha[4]+0.5*alpha[0])*fUpwind[7]+0.4*(alpha[3]*fUpwind[6]+fUpwind[3]*alpha[6])+0.5000000000000001*alpha[1]*fUpwind[5]+0.447213595499958*(alpha[2]*fUpwind[3]+fUpwind[2]*alpha[3]); 

  out[0] += 0.7071067811865475*Ghat[0]*dv10; 
  out[1] += 0.7071067811865475*Ghat[1]*dv10; 
  out[2] += -1.224744871391589*Ghat[0]*dv10; 
  out[3] += 0.7071067811865475*Ghat[2]*dv10; 
  out[4] += -1.224744871391589*Ghat[1]*dv10; 
  out[5] += 0.7071067811865475*Ghat[3]*dv10; 
  out[6] += -1.224744871391589*Ghat[2]*dv10; 
  out[7] += 0.7071067811865475*Ghat[4]*dv10; 
  out[8] += 1.58113883008419*Ghat[0]*dv10; 
  out[9] += 0.7071067811865475*Ghat[5]*dv10; 
  out[10] += -1.224744871391589*Ghat[3]*dv10; 
  out[11] += -1.224744871391589*Ghat[4]*dv10; 
  out[12] += 1.58113883008419*Ghat[1]*dv10; 
  out[13] += 0.7071067811865475*Ghat[6]*dv10; 
  out[14] += 1.58113883008419*Ghat[2]*dv10; 
  out[15] += 0.7071067811865475*Ghat[7]*dv10; 
  out[16] += -1.224744871391589*Ghat[5]*dv10; 
  out[17] += -1.224744871391589*Ghat[6]*dv10; 
  out[18] += 1.58113883008419*Ghat[3]*dv10; 
  out[19] += -1.224744871391589*Ghat[7]*dv10; 

  } 
} 
