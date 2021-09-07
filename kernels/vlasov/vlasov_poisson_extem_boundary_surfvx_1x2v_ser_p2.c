#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_ser_1x2v_p2_surfvx_quad.h> 
GKYL_CU_DH void vlasov_poisson_extem_boundary_surfvx_1x2v_ser_p2(const double *w, const double *dxv, const double *fac_phi, const double *vecA, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w:           Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // fac_phi:     potential (scaled by appropriate factors).
  // vecA:        vector potential (scaled by appropriate factors). Unused in pure Vlasov-Poisson. 
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Output distribution function in skin cell 
  const double dv10 = 2/dxv[1]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv2 = dxv[2], wv2 = w[2]; 
  const double *phi = &fac_phi[0]; 
  const double dx10 = 2/dxv[0]; 
  const double *A0 = &vecA[0]; 
  const double *A1 = &vecA[3]; 
  double alpha[8] = {0.0}; 

  alpha[0] = 2.449489742783178*A1[1]*dx10*wv2-2.449489742783178*phi[1]*dx10; 
  alpha[1] = 5.477225575051662*A1[2]*dx10*wv2-5.477225575051662*phi[2]*dx10; 
  alpha[2] = 0.7071067811865475*A1[1]*dv2*dx10; 
  alpha[3] = 1.58113883008419*A1[2]*dv2*dx10; 

  double fUpwindQuad[9] = {0.0};
  double fUpwind[8] = {0.0};
  double Ghat[8] = {0.0}; 

  if (edge == -1) { 

  if (0.9*alpha[3]-0.6708203932499369*(alpha[2]+alpha[1])+0.5*alpha[0] > 0) { 

    fUpwindQuad[0] = ser_1x2v_p2_surfvx_quad_0(1, fSkin); 
  } else { 

    fUpwindQuad[0] = ser_1x2v_p2_surfvx_quad_0(-1, fEdge); 
  } 
  if (0.5*alpha[0]-0.6708203932499369*alpha[2] > 0) { 

    fUpwindQuad[1] = ser_1x2v_p2_surfvx_quad_1(1, fSkin); 
  } else { 

    fUpwindQuad[1] = ser_1x2v_p2_surfvx_quad_1(-1, fEdge); 
  } 
  if ((-0.9*alpha[3])-0.6708203932499369*alpha[2]+0.6708203932499369*alpha[1]+0.5*alpha[0] > 0) { 

    fUpwindQuad[2] = ser_1x2v_p2_surfvx_quad_2(1, fSkin); 
  } else { 

    fUpwindQuad[2] = ser_1x2v_p2_surfvx_quad_2(-1, fEdge); 
  } 
  if (0.5*alpha[0]-0.6708203932499369*alpha[1] > 0) { 

    fUpwindQuad[3] = ser_1x2v_p2_surfvx_quad_3(1, fSkin); 
  } else { 

    fUpwindQuad[3] = ser_1x2v_p2_surfvx_quad_3(-1, fEdge); 
  } 
  if (0.5*alpha[0] > 0) { 

    fUpwindQuad[4] = ser_1x2v_p2_surfvx_quad_4(1, fSkin); 
  } else { 

    fUpwindQuad[4] = ser_1x2v_p2_surfvx_quad_4(-1, fEdge); 
  } 
  if (0.6708203932499369*alpha[1]+0.5*alpha[0] > 0) { 

    fUpwindQuad[5] = ser_1x2v_p2_surfvx_quad_5(1, fSkin); 
  } else { 

    fUpwindQuad[5] = ser_1x2v_p2_surfvx_quad_5(-1, fEdge); 
  } 
  if ((-0.9*alpha[3])+0.6708203932499369*alpha[2]-0.6708203932499369*alpha[1]+0.5*alpha[0] > 0) { 

    fUpwindQuad[6] = ser_1x2v_p2_surfvx_quad_6(1, fSkin); 
  } else { 

    fUpwindQuad[6] = ser_1x2v_p2_surfvx_quad_6(-1, fEdge); 
  } 
  if (0.6708203932499369*alpha[2]+0.5*alpha[0] > 0) { 

    fUpwindQuad[7] = ser_1x2v_p2_surfvx_quad_7(1, fSkin); 
  } else { 

    fUpwindQuad[7] = ser_1x2v_p2_surfvx_quad_7(-1, fEdge); 
  } 
  if (0.9*alpha[3]+0.6708203932499369*(alpha[2]+alpha[1])+0.5*alpha[0] > 0) { 

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

  Ghat[0] += 0.5*(alpha[3]*fUpwind[3]+alpha[2]*fUpwind[2]+alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] += 0.447213595499958*alpha[3]*fUpwind[6]+0.4472135954999579*alpha[1]*fUpwind[4]+0.5*(alpha[2]*fUpwind[3]+fUpwind[2]*alpha[3]+alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]); 
  Ghat[2] += 0.447213595499958*alpha[3]*fUpwind[7]+0.4472135954999579*alpha[2]*fUpwind[5]+0.5*(alpha[1]*fUpwind[3]+fUpwind[1]*alpha[3]+alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2]); 
  Ghat[3] += 0.447213595499958*(alpha[2]*fUpwind[7]+alpha[1]*fUpwind[6])+0.4472135954999579*alpha[3]*(fUpwind[5]+fUpwind[4])+0.5*(alpha[0]*fUpwind[3]+fUpwind[0]*alpha[3]+alpha[1]*fUpwind[2]+fUpwind[1]*alpha[2]); 
  Ghat[4] += 0.5000000000000001*alpha[2]*fUpwind[6]+0.5*alpha[0]*fUpwind[4]+0.4472135954999579*(alpha[3]*fUpwind[3]+alpha[1]*fUpwind[1]); 
  Ghat[5] += 0.5000000000000001*alpha[1]*fUpwind[7]+0.5*alpha[0]*fUpwind[5]+0.4472135954999579*(alpha[3]*fUpwind[3]+alpha[2]*fUpwind[2]); 
  Ghat[6] += 0.4*alpha[3]*fUpwind[7]+0.5*alpha[0]*fUpwind[6]+0.5000000000000001*alpha[2]*fUpwind[4]+0.447213595499958*(alpha[1]*fUpwind[3]+fUpwind[1]*alpha[3]); 
  Ghat[7] += 0.5*alpha[0]*fUpwind[7]+0.4*alpha[3]*fUpwind[6]+0.5000000000000001*alpha[1]*fUpwind[5]+0.447213595499958*(alpha[2]*fUpwind[3]+fUpwind[2]*alpha[3]); 

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

  if (0.9*alpha[3]-0.6708203932499369*(alpha[2]+alpha[1])+0.5*alpha[0] > 0) { 

    fUpwindQuad[0] = ser_1x2v_p2_surfvx_quad_0(1, fEdge); 
  } else { 

    fUpwindQuad[0] = ser_1x2v_p2_surfvx_quad_0(-1, fSkin); 
  } 
  if (0.5*alpha[0]-0.6708203932499369*alpha[2] > 0) { 

    fUpwindQuad[1] = ser_1x2v_p2_surfvx_quad_1(1, fEdge); 
  } else { 

    fUpwindQuad[1] = ser_1x2v_p2_surfvx_quad_1(-1, fSkin); 
  } 
  if ((-0.9*alpha[3])-0.6708203932499369*alpha[2]+0.6708203932499369*alpha[1]+0.5*alpha[0] > 0) { 

    fUpwindQuad[2] = ser_1x2v_p2_surfvx_quad_2(1, fEdge); 
  } else { 

    fUpwindQuad[2] = ser_1x2v_p2_surfvx_quad_2(-1, fSkin); 
  } 
  if (0.5*alpha[0]-0.6708203932499369*alpha[1] > 0) { 

    fUpwindQuad[3] = ser_1x2v_p2_surfvx_quad_3(1, fEdge); 
  } else { 

    fUpwindQuad[3] = ser_1x2v_p2_surfvx_quad_3(-1, fSkin); 
  } 
  if (0.5*alpha[0] > 0) { 

    fUpwindQuad[4] = ser_1x2v_p2_surfvx_quad_4(1, fEdge); 
  } else { 

    fUpwindQuad[4] = ser_1x2v_p2_surfvx_quad_4(-1, fSkin); 
  } 
  if (0.6708203932499369*alpha[1]+0.5*alpha[0] > 0) { 

    fUpwindQuad[5] = ser_1x2v_p2_surfvx_quad_5(1, fEdge); 
  } else { 

    fUpwindQuad[5] = ser_1x2v_p2_surfvx_quad_5(-1, fSkin); 
  } 
  if ((-0.9*alpha[3])+0.6708203932499369*alpha[2]-0.6708203932499369*alpha[1]+0.5*alpha[0] > 0) { 

    fUpwindQuad[6] = ser_1x2v_p2_surfvx_quad_6(1, fEdge); 
  } else { 

    fUpwindQuad[6] = ser_1x2v_p2_surfvx_quad_6(-1, fSkin); 
  } 
  if (0.6708203932499369*alpha[2]+0.5*alpha[0] > 0) { 

    fUpwindQuad[7] = ser_1x2v_p2_surfvx_quad_7(1, fEdge); 
  } else { 

    fUpwindQuad[7] = ser_1x2v_p2_surfvx_quad_7(-1, fSkin); 
  } 
  if (0.9*alpha[3]+0.6708203932499369*(alpha[2]+alpha[1])+0.5*alpha[0] > 0) { 

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

  Ghat[0] += 0.5*(alpha[3]*fUpwind[3]+alpha[2]*fUpwind[2]+alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] += 0.447213595499958*alpha[3]*fUpwind[6]+0.4472135954999579*alpha[1]*fUpwind[4]+0.5*(alpha[2]*fUpwind[3]+fUpwind[2]*alpha[3]+alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]); 
  Ghat[2] += 0.447213595499958*alpha[3]*fUpwind[7]+0.4472135954999579*alpha[2]*fUpwind[5]+0.5*(alpha[1]*fUpwind[3]+fUpwind[1]*alpha[3]+alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2]); 
  Ghat[3] += 0.447213595499958*(alpha[2]*fUpwind[7]+alpha[1]*fUpwind[6])+0.4472135954999579*alpha[3]*(fUpwind[5]+fUpwind[4])+0.5*(alpha[0]*fUpwind[3]+fUpwind[0]*alpha[3]+alpha[1]*fUpwind[2]+fUpwind[1]*alpha[2]); 
  Ghat[4] += 0.5000000000000001*alpha[2]*fUpwind[6]+0.5*alpha[0]*fUpwind[4]+0.4472135954999579*(alpha[3]*fUpwind[3]+alpha[1]*fUpwind[1]); 
  Ghat[5] += 0.5000000000000001*alpha[1]*fUpwind[7]+0.5*alpha[0]*fUpwind[5]+0.4472135954999579*(alpha[3]*fUpwind[3]+alpha[2]*fUpwind[2]); 
  Ghat[6] += 0.4*alpha[3]*fUpwind[7]+0.5*alpha[0]*fUpwind[6]+0.5000000000000001*alpha[2]*fUpwind[4]+0.447213595499958*(alpha[1]*fUpwind[3]+fUpwind[1]*alpha[3]); 
  Ghat[7] += 0.5*alpha[0]*fUpwind[7]+0.4*alpha[3]*fUpwind[6]+0.5000000000000001*alpha[1]*fUpwind[5]+0.447213595499958*(alpha[2]*fUpwind[3]+fUpwind[2]*alpha[3]); 

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
