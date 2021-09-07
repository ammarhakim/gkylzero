#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_ser_3x3v_p1_surfvy_quad.h> 
GKYL_CU_DH void vlasov_poisson_extem_boundary_surfvy_3x3v_ser_p1(const double *w, const double *dxv, const double *fac_phi, const double *vecA, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w:           Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // fac_phi:     potential (scaled by appropriate factors).
  // vecA:        vector potential (scaled by appropriate factors). Unused in pure Vlasov-Poisson. 
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Output distribution function in skin cell 
  const double dv11 = 2/dxv[4]; 
  const double dv1 = dxv[3], wv1 = w[3]; 
  const double dv2 = dxv[4], wv2 = w[4]; 
  const double dv3 = dxv[5], wv3 = w[5]; 
  const double *phi = &fac_phi[0]; 
  const double dx10 = 2/dxv[0]; 
  const double dx11 = 2/dxv[1]; 
  const double dx12 = 2/dxv[2]; 
  const double *A0 = &vecA[0]; 
  const double *A1 = &vecA[8]; 
  const double *A2 = &vecA[16]; 
  double alpha[32] = {0.0}; 

  alpha[0] = (-3.464101615137754*A1[3]*dx12*wv3)+3.464101615137754*A2[2]*dx11*wv3+3.464101615137754*A0[2]*dx11*wv1-3.464101615137754*A1[1]*dx10*wv1-3.464101615137754*phi[2]*dx11; 
  alpha[1] = (-3.464101615137754*A1[5]*dx12*wv3)+3.464101615137754*A2[4]*dx11*wv3+3.464101615137754*A0[4]*dx11*wv1-3.464101615137754*phi[4]*dx11; 
  alpha[2] = (-3.464101615137754*A1[6]*dx12*wv3)-3.464101615137754*A1[4]*dx10*wv1; 
  alpha[3] = 3.464101615137754*A2[6]*dx11*wv3+3.464101615137754*A0[6]*dx11*wv1-3.464101615137754*A1[5]*dx10*wv1-3.464101615137754*phi[6]*dx11; 
  alpha[4] = A0[2]*dv1*dx11-1.0*A1[1]*dv1*dx10; 
  alpha[5] = A2[2]*dv3*dx11-1.0*A1[3]*dv3*dx12; 
  alpha[6] = -3.464101615137754*A1[7]*dx12*wv3; 
  alpha[7] = 3.464101615137754*A2[7]*dx11*wv3+3.464101615137754*A0[7]*dx11*wv1-3.464101615137754*phi[7]*dx11; 
  alpha[8] = -3.464101615137754*A1[7]*dx10*wv1; 
  alpha[9] = A0[4]*dv1*dx11; 
  alpha[10] = -1.0*A1[4]*dv1*dx10; 
  alpha[11] = A0[6]*dv1*dx11-1.0*A1[5]*dv1*dx10; 
  alpha[12] = A2[4]*dv3*dx11-1.0*A1[5]*dv3*dx12; 
  alpha[13] = -1.0*A1[6]*dv3*dx12; 
  alpha[14] = A2[6]*dv3*dx11; 
  alpha[18] = A0[7]*dv1*dx11; 
  alpha[19] = -1.0*A1[7]*dv1*dx10; 
  alpha[20] = -1.0*A1[7]*dv3*dx12; 
  alpha[21] = A2[7]*dv3*dx11; 

  double fUpwindQuad[32] = {0.0};
  double fUpwind[32] = {0.0};
  double Ghat[32] = {0.0}; 

  if (edge == -1) { 

  if ((-alpha[21])-alpha[20]-alpha[19]-alpha[18]+alpha[14]+alpha[13]+alpha[12]+alpha[11]+alpha[10]+alpha[9]+alpha[8]+alpha[7]+alpha[6]-alpha[5]-alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[0] = ser_3x3v_p1_surfvy_quad_0(1, fSkin); 
  } else { 

    fUpwindQuad[0] = ser_3x3v_p1_surfvy_quad_0(-1, fEdge); 
  } 
  if (alpha[21]+alpha[20]-alpha[19]+alpha[18]+alpha[14]+alpha[13]-alpha[12]+alpha[11]+alpha[10]-alpha[9]+alpha[8]-alpha[7]-alpha[6]-alpha[5]-alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[1] = ser_3x3v_p1_surfvy_quad_1(1, fSkin); 
  } else { 

    fUpwindQuad[1] = ser_3x3v_p1_surfvy_quad_1(-1, fEdge); 
  } 
  if ((-alpha[21])+alpha[20]+alpha[19]-alpha[18]+alpha[14]-alpha[13]+alpha[12]+alpha[11]-alpha[10]+alpha[9]-alpha[8]+alpha[7]-alpha[6]-alpha[5]-alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[2] = ser_3x3v_p1_surfvy_quad_2(1, fSkin); 
  } else { 

    fUpwindQuad[2] = ser_3x3v_p1_surfvy_quad_2(-1, fEdge); 
  } 
  if (alpha[21]-alpha[20]+alpha[19]+alpha[18]+alpha[14]-alpha[13]-alpha[12]+alpha[11]-alpha[10]-alpha[9]-alpha[8]-alpha[7]+alpha[6]-alpha[5]-alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[3] = ser_3x3v_p1_surfvy_quad_3(1, fSkin); 
  } else { 

    fUpwindQuad[3] = ser_3x3v_p1_surfvy_quad_3(-1, fEdge); 
  } 
  if (alpha[21]-alpha[20]+alpha[19]+alpha[18]-alpha[14]+alpha[13]+alpha[12]-alpha[11]+alpha[10]+alpha[9]-alpha[8]-alpha[7]+alpha[6]-alpha[5]-alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[4] = ser_3x3v_p1_surfvy_quad_4(1, fSkin); 
  } else { 

    fUpwindQuad[4] = ser_3x3v_p1_surfvy_quad_4(-1, fEdge); 
  } 
  if ((-alpha[21])+alpha[20]+alpha[19]-alpha[18]-alpha[14]+alpha[13]-alpha[12]-alpha[11]+alpha[10]-alpha[9]-alpha[8]+alpha[7]-alpha[6]-alpha[5]-alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[5] = ser_3x3v_p1_surfvy_quad_5(1, fSkin); 
  } else { 

    fUpwindQuad[5] = ser_3x3v_p1_surfvy_quad_5(-1, fEdge); 
  } 
  if (alpha[21]+alpha[20]-alpha[19]+alpha[18]-alpha[14]-alpha[13]+alpha[12]-alpha[11]-alpha[10]+alpha[9]+alpha[8]-alpha[7]-alpha[6]-alpha[5]-alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[6] = ser_3x3v_p1_surfvy_quad_6(1, fSkin); 
  } else { 

    fUpwindQuad[6] = ser_3x3v_p1_surfvy_quad_6(-1, fEdge); 
  } 
  if ((-alpha[21])-alpha[20]-alpha[19]-alpha[18]-alpha[14]-alpha[13]-alpha[12]-alpha[11]-alpha[10]-alpha[9]+alpha[8]+alpha[7]+alpha[6]-alpha[5]-alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[7] = ser_3x3v_p1_surfvy_quad_7(1, fSkin); 
  } else { 

    fUpwindQuad[7] = ser_3x3v_p1_surfvy_quad_7(-1, fEdge); 
  } 
  if ((-alpha[21])-alpha[20]+alpha[19]+alpha[18]+alpha[14]+alpha[13]+alpha[12]-alpha[11]-alpha[10]-alpha[9]+alpha[8]+alpha[7]+alpha[6]-alpha[5]+alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[8] = ser_3x3v_p1_surfvy_quad_8(1, fSkin); 
  } else { 

    fUpwindQuad[8] = ser_3x3v_p1_surfvy_quad_8(-1, fEdge); 
  } 
  if (alpha[21]+alpha[20]+alpha[19]-alpha[18]+alpha[14]+alpha[13]-alpha[12]-alpha[11]-alpha[10]+alpha[9]+alpha[8]-alpha[7]-alpha[6]-alpha[5]+alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[9] = ser_3x3v_p1_surfvy_quad_9(1, fSkin); 
  } else { 

    fUpwindQuad[9] = ser_3x3v_p1_surfvy_quad_9(-1, fEdge); 
  } 
  if ((-alpha[21])+alpha[20]-alpha[19]+alpha[18]+alpha[14]-alpha[13]+alpha[12]-alpha[11]+alpha[10]-alpha[9]-alpha[8]+alpha[7]-alpha[6]-alpha[5]+alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[10] = ser_3x3v_p1_surfvy_quad_10(1, fSkin); 
  } else { 

    fUpwindQuad[10] = ser_3x3v_p1_surfvy_quad_10(-1, fEdge); 
  } 
  if (alpha[21]-alpha[20]-alpha[19]-alpha[18]+alpha[14]-alpha[13]-alpha[12]-alpha[11]+alpha[10]+alpha[9]-alpha[8]-alpha[7]+alpha[6]-alpha[5]+alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[11] = ser_3x3v_p1_surfvy_quad_11(1, fSkin); 
  } else { 

    fUpwindQuad[11] = ser_3x3v_p1_surfvy_quad_11(-1, fEdge); 
  } 
  if (alpha[21]-alpha[20]-alpha[19]-alpha[18]-alpha[14]+alpha[13]+alpha[12]+alpha[11]-alpha[10]-alpha[9]-alpha[8]-alpha[7]+alpha[6]-alpha[5]+alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[12] = ser_3x3v_p1_surfvy_quad_12(1, fSkin); 
  } else { 

    fUpwindQuad[12] = ser_3x3v_p1_surfvy_quad_12(-1, fEdge); 
  } 
  if ((-alpha[21])+alpha[20]-alpha[19]+alpha[18]-alpha[14]+alpha[13]-alpha[12]+alpha[11]-alpha[10]+alpha[9]-alpha[8]+alpha[7]-alpha[6]-alpha[5]+alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[13] = ser_3x3v_p1_surfvy_quad_13(1, fSkin); 
  } else { 

    fUpwindQuad[13] = ser_3x3v_p1_surfvy_quad_13(-1, fEdge); 
  } 
  if (alpha[21]+alpha[20]+alpha[19]-alpha[18]-alpha[14]-alpha[13]+alpha[12]+alpha[11]+alpha[10]-alpha[9]+alpha[8]-alpha[7]-alpha[6]-alpha[5]+alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[14] = ser_3x3v_p1_surfvy_quad_14(1, fSkin); 
  } else { 

    fUpwindQuad[14] = ser_3x3v_p1_surfvy_quad_14(-1, fEdge); 
  } 
  if ((-alpha[21])-alpha[20]+alpha[19]+alpha[18]-alpha[14]-alpha[13]-alpha[12]+alpha[11]+alpha[10]+alpha[9]+alpha[8]+alpha[7]+alpha[6]-alpha[5]+alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[15] = ser_3x3v_p1_surfvy_quad_15(1, fSkin); 
  } else { 

    fUpwindQuad[15] = ser_3x3v_p1_surfvy_quad_15(-1, fEdge); 
  } 
  if (alpha[21]+alpha[20]-alpha[19]-alpha[18]-alpha[14]-alpha[13]-alpha[12]+alpha[11]+alpha[10]+alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5]-alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[16] = ser_3x3v_p1_surfvy_quad_16(1, fSkin); 
  } else { 

    fUpwindQuad[16] = ser_3x3v_p1_surfvy_quad_16(-1, fEdge); 
  } 
  if ((-alpha[21])-alpha[20]-alpha[19]+alpha[18]-alpha[14]-alpha[13]+alpha[12]+alpha[11]+alpha[10]-alpha[9]+alpha[8]-alpha[7]-alpha[6]+alpha[5]-alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[17] = ser_3x3v_p1_surfvy_quad_17(1, fSkin); 
  } else { 

    fUpwindQuad[17] = ser_3x3v_p1_surfvy_quad_17(-1, fEdge); 
  } 
  if (alpha[21]-alpha[20]+alpha[19]-alpha[18]-alpha[14]+alpha[13]-alpha[12]+alpha[11]-alpha[10]+alpha[9]-alpha[8]+alpha[7]-alpha[6]+alpha[5]-alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[18] = ser_3x3v_p1_surfvy_quad_18(1, fSkin); 
  } else { 

    fUpwindQuad[18] = ser_3x3v_p1_surfvy_quad_18(-1, fEdge); 
  } 
  if ((-alpha[21])+alpha[20]+alpha[19]+alpha[18]-alpha[14]+alpha[13]+alpha[12]+alpha[11]-alpha[10]-alpha[9]-alpha[8]-alpha[7]+alpha[6]+alpha[5]-alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[19] = ser_3x3v_p1_surfvy_quad_19(1, fSkin); 
  } else { 

    fUpwindQuad[19] = ser_3x3v_p1_surfvy_quad_19(-1, fEdge); 
  } 
  if ((-alpha[21])+alpha[20]+alpha[19]+alpha[18]+alpha[14]-alpha[13]-alpha[12]-alpha[11]+alpha[10]+alpha[9]-alpha[8]-alpha[7]+alpha[6]+alpha[5]-alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[20] = ser_3x3v_p1_surfvy_quad_20(1, fSkin); 
  } else { 

    fUpwindQuad[20] = ser_3x3v_p1_surfvy_quad_20(-1, fEdge); 
  } 
  if (alpha[21]-alpha[20]+alpha[19]-alpha[18]+alpha[14]-alpha[13]+alpha[12]-alpha[11]+alpha[10]-alpha[9]-alpha[8]+alpha[7]-alpha[6]+alpha[5]-alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[21] = ser_3x3v_p1_surfvy_quad_21(1, fSkin); 
  } else { 

    fUpwindQuad[21] = ser_3x3v_p1_surfvy_quad_21(-1, fEdge); 
  } 
  if ((-alpha[21])-alpha[20]-alpha[19]+alpha[18]+alpha[14]+alpha[13]-alpha[12]-alpha[11]-alpha[10]+alpha[9]+alpha[8]-alpha[7]-alpha[6]+alpha[5]-alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[22] = ser_3x3v_p1_surfvy_quad_22(1, fSkin); 
  } else { 

    fUpwindQuad[22] = ser_3x3v_p1_surfvy_quad_22(-1, fEdge); 
  } 
  if (alpha[21]+alpha[20]-alpha[19]-alpha[18]+alpha[14]+alpha[13]+alpha[12]-alpha[11]-alpha[10]-alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5]-alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[23] = ser_3x3v_p1_surfvy_quad_23(1, fSkin); 
  } else { 

    fUpwindQuad[23] = ser_3x3v_p1_surfvy_quad_23(-1, fEdge); 
  } 
  if (alpha[21]+alpha[20]+alpha[19]+alpha[18]-alpha[14]-alpha[13]-alpha[12]-alpha[11]-alpha[10]-alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5]+alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[24] = ser_3x3v_p1_surfvy_quad_24(1, fSkin); 
  } else { 

    fUpwindQuad[24] = ser_3x3v_p1_surfvy_quad_24(-1, fEdge); 
  } 
  if ((-alpha[21])-alpha[20]+alpha[19]-alpha[18]-alpha[14]-alpha[13]+alpha[12]-alpha[11]-alpha[10]+alpha[9]+alpha[8]-alpha[7]-alpha[6]+alpha[5]+alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[25] = ser_3x3v_p1_surfvy_quad_25(1, fSkin); 
  } else { 

    fUpwindQuad[25] = ser_3x3v_p1_surfvy_quad_25(-1, fEdge); 
  } 
  if (alpha[21]-alpha[20]-alpha[19]+alpha[18]-alpha[14]+alpha[13]-alpha[12]-alpha[11]+alpha[10]-alpha[9]-alpha[8]+alpha[7]-alpha[6]+alpha[5]+alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[26] = ser_3x3v_p1_surfvy_quad_26(1, fSkin); 
  } else { 

    fUpwindQuad[26] = ser_3x3v_p1_surfvy_quad_26(-1, fEdge); 
  } 
  if ((-alpha[21])+alpha[20]-alpha[19]-alpha[18]-alpha[14]+alpha[13]+alpha[12]-alpha[11]+alpha[10]+alpha[9]-alpha[8]-alpha[7]+alpha[6]+alpha[5]+alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[27] = ser_3x3v_p1_surfvy_quad_27(1, fSkin); 
  } else { 

    fUpwindQuad[27] = ser_3x3v_p1_surfvy_quad_27(-1, fEdge); 
  } 
  if ((-alpha[21])+alpha[20]-alpha[19]-alpha[18]+alpha[14]-alpha[13]-alpha[12]+alpha[11]-alpha[10]-alpha[9]-alpha[8]-alpha[7]+alpha[6]+alpha[5]+alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[28] = ser_3x3v_p1_surfvy_quad_28(1, fSkin); 
  } else { 

    fUpwindQuad[28] = ser_3x3v_p1_surfvy_quad_28(-1, fEdge); 
  } 
  if (alpha[21]-alpha[20]-alpha[19]+alpha[18]+alpha[14]-alpha[13]+alpha[12]+alpha[11]-alpha[10]+alpha[9]-alpha[8]+alpha[7]-alpha[6]+alpha[5]+alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[29] = ser_3x3v_p1_surfvy_quad_29(1, fSkin); 
  } else { 

    fUpwindQuad[29] = ser_3x3v_p1_surfvy_quad_29(-1, fEdge); 
  } 
  if ((-alpha[21])-alpha[20]+alpha[19]-alpha[18]+alpha[14]+alpha[13]-alpha[12]+alpha[11]+alpha[10]-alpha[9]+alpha[8]-alpha[7]-alpha[6]+alpha[5]+alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[30] = ser_3x3v_p1_surfvy_quad_30(1, fSkin); 
  } else { 

    fUpwindQuad[30] = ser_3x3v_p1_surfvy_quad_30(-1, fEdge); 
  } 
  if (alpha[21]+alpha[20]+alpha[19]+alpha[18]+alpha[14]+alpha[13]+alpha[12]+alpha[11]+alpha[10]+alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5]+alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[31] = ser_3x3v_p1_surfvy_quad_31(1, fSkin); 
  } else { 

    fUpwindQuad[31] = ser_3x3v_p1_surfvy_quad_31(-1, fEdge); 
  } 

  fUpwind[0] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]+fUpwindQuad[29]+fUpwindQuad[28]+fUpwindQuad[27]+fUpwindQuad[26]+fUpwindQuad[25]+fUpwindQuad[24]+fUpwindQuad[23]+fUpwindQuad[22]+fUpwindQuad[21]+fUpwindQuad[20]+fUpwindQuad[19]+fUpwindQuad[18]+fUpwindQuad[17]+fUpwindQuad[16]+fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[1] = 0.1767766952966368*(fUpwindQuad[31]-1.0*fUpwindQuad[30]+fUpwindQuad[29]-1.0*fUpwindQuad[28]+fUpwindQuad[27]-1.0*fUpwindQuad[26]+fUpwindQuad[25]-1.0*fUpwindQuad[24]+fUpwindQuad[23]-1.0*fUpwindQuad[22]+fUpwindQuad[21]-1.0*fUpwindQuad[20]+fUpwindQuad[19]-1.0*fUpwindQuad[18]+fUpwindQuad[17]-1.0*fUpwindQuad[16]+fUpwindQuad[15]-1.0*fUpwindQuad[14]+fUpwindQuad[13]-1.0*fUpwindQuad[12]+fUpwindQuad[11]-1.0*fUpwindQuad[10]+fUpwindQuad[9]-1.0*fUpwindQuad[8]+fUpwindQuad[7]-1.0*fUpwindQuad[6]+fUpwindQuad[5]-1.0*fUpwindQuad[4]+fUpwindQuad[3]-1.0*fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[2] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]-1.0*(fUpwindQuad[29]+fUpwindQuad[28])+fUpwindQuad[27]+fUpwindQuad[26]-1.0*(fUpwindQuad[25]+fUpwindQuad[24])+fUpwindQuad[23]+fUpwindQuad[22]-1.0*(fUpwindQuad[21]+fUpwindQuad[20])+fUpwindQuad[19]+fUpwindQuad[18]-1.0*(fUpwindQuad[17]+fUpwindQuad[16])+fUpwindQuad[15]+fUpwindQuad[14]-1.0*(fUpwindQuad[13]+fUpwindQuad[12])+fUpwindQuad[11]+fUpwindQuad[10]-1.0*(fUpwindQuad[9]+fUpwindQuad[8])+fUpwindQuad[7]+fUpwindQuad[6]-1.0*(fUpwindQuad[5]+fUpwindQuad[4])+fUpwindQuad[3]+fUpwindQuad[2]-1.0*(fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[3] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]+fUpwindQuad[29]+fUpwindQuad[28]-1.0*(fUpwindQuad[27]+fUpwindQuad[26]+fUpwindQuad[25]+fUpwindQuad[24])+fUpwindQuad[23]+fUpwindQuad[22]+fUpwindQuad[21]+fUpwindQuad[20]-1.0*(fUpwindQuad[19]+fUpwindQuad[18]+fUpwindQuad[17]+fUpwindQuad[16])+fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12]-1.0*(fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8])+fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]-1.0*(fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[4] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]+fUpwindQuad[29]+fUpwindQuad[28]+fUpwindQuad[27]+fUpwindQuad[26]+fUpwindQuad[25]+fUpwindQuad[24]-1.0*(fUpwindQuad[23]+fUpwindQuad[22]+fUpwindQuad[21]+fUpwindQuad[20]+fUpwindQuad[19]+fUpwindQuad[18]+fUpwindQuad[17]+fUpwindQuad[16])+fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8]-1.0*(fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[5] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]+fUpwindQuad[29]+fUpwindQuad[28]+fUpwindQuad[27]+fUpwindQuad[26]+fUpwindQuad[25]+fUpwindQuad[24]+fUpwindQuad[23]+fUpwindQuad[22]+fUpwindQuad[21]+fUpwindQuad[20]+fUpwindQuad[19]+fUpwindQuad[18]+fUpwindQuad[17]+fUpwindQuad[16]-1.0*(fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[6] = 0.1767766952966368*(fUpwindQuad[31]-1.0*(fUpwindQuad[30]+fUpwindQuad[29])+fUpwindQuad[28]+fUpwindQuad[27]-1.0*(fUpwindQuad[26]+fUpwindQuad[25])+fUpwindQuad[24]+fUpwindQuad[23]-1.0*(fUpwindQuad[22]+fUpwindQuad[21])+fUpwindQuad[20]+fUpwindQuad[19]-1.0*(fUpwindQuad[18]+fUpwindQuad[17])+fUpwindQuad[16]+fUpwindQuad[15]-1.0*(fUpwindQuad[14]+fUpwindQuad[13])+fUpwindQuad[12]+fUpwindQuad[11]-1.0*(fUpwindQuad[10]+fUpwindQuad[9])+fUpwindQuad[8]+fUpwindQuad[7]-1.0*(fUpwindQuad[6]+fUpwindQuad[5])+fUpwindQuad[4]+fUpwindQuad[3]-1.0*(fUpwindQuad[2]+fUpwindQuad[1])+fUpwindQuad[0]); 
  fUpwind[7] = 0.1767766952966368*(fUpwindQuad[31]-1.0*fUpwindQuad[30]+fUpwindQuad[29]-1.0*(fUpwindQuad[28]+fUpwindQuad[27])+fUpwindQuad[26]-1.0*fUpwindQuad[25]+fUpwindQuad[24]+fUpwindQuad[23]-1.0*fUpwindQuad[22]+fUpwindQuad[21]-1.0*(fUpwindQuad[20]+fUpwindQuad[19])+fUpwindQuad[18]-1.0*fUpwindQuad[17]+fUpwindQuad[16]+fUpwindQuad[15]-1.0*fUpwindQuad[14]+fUpwindQuad[13]-1.0*(fUpwindQuad[12]+fUpwindQuad[11])+fUpwindQuad[10]-1.0*fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]-1.0*fUpwindQuad[6]+fUpwindQuad[5]-1.0*(fUpwindQuad[4]+fUpwindQuad[3])+fUpwindQuad[2]-1.0*fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[8] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]-1.0*(fUpwindQuad[29]+fUpwindQuad[28]+fUpwindQuad[27]+fUpwindQuad[26])+fUpwindQuad[25]+fUpwindQuad[24]+fUpwindQuad[23]+fUpwindQuad[22]-1.0*(fUpwindQuad[21]+fUpwindQuad[20]+fUpwindQuad[19]+fUpwindQuad[18])+fUpwindQuad[17]+fUpwindQuad[16]+fUpwindQuad[15]+fUpwindQuad[14]-1.0*(fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10])+fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6]-1.0*(fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2])+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[9] = 0.1767766952966368*(fUpwindQuad[31]-1.0*fUpwindQuad[30]+fUpwindQuad[29]-1.0*fUpwindQuad[28]+fUpwindQuad[27]-1.0*fUpwindQuad[26]+fUpwindQuad[25]-1.0*(fUpwindQuad[24]+fUpwindQuad[23])+fUpwindQuad[22]-1.0*fUpwindQuad[21]+fUpwindQuad[20]-1.0*fUpwindQuad[19]+fUpwindQuad[18]-1.0*fUpwindQuad[17]+fUpwindQuad[16]+fUpwindQuad[15]-1.0*fUpwindQuad[14]+fUpwindQuad[13]-1.0*fUpwindQuad[12]+fUpwindQuad[11]-1.0*fUpwindQuad[10]+fUpwindQuad[9]-1.0*(fUpwindQuad[8]+fUpwindQuad[7])+fUpwindQuad[6]-1.0*fUpwindQuad[5]+fUpwindQuad[4]-1.0*fUpwindQuad[3]+fUpwindQuad[2]-1.0*fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[10] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]-1.0*(fUpwindQuad[29]+fUpwindQuad[28])+fUpwindQuad[27]+fUpwindQuad[26]-1.0*(fUpwindQuad[25]+fUpwindQuad[24]+fUpwindQuad[23]+fUpwindQuad[22])+fUpwindQuad[21]+fUpwindQuad[20]-1.0*(fUpwindQuad[19]+fUpwindQuad[18])+fUpwindQuad[17]+fUpwindQuad[16]+fUpwindQuad[15]+fUpwindQuad[14]-1.0*(fUpwindQuad[13]+fUpwindQuad[12])+fUpwindQuad[11]+fUpwindQuad[10]-1.0*(fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6])+fUpwindQuad[5]+fUpwindQuad[4]-1.0*(fUpwindQuad[3]+fUpwindQuad[2])+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[11] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]+fUpwindQuad[29]+fUpwindQuad[28]-1.0*(fUpwindQuad[27]+fUpwindQuad[26]+fUpwindQuad[25]+fUpwindQuad[24]+fUpwindQuad[23]+fUpwindQuad[22]+fUpwindQuad[21]+fUpwindQuad[20])+fUpwindQuad[19]+fUpwindQuad[18]+fUpwindQuad[17]+fUpwindQuad[16]+fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12]-1.0*(fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4])+fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[12] = 0.1767766952966368*(fUpwindQuad[31]-1.0*fUpwindQuad[30]+fUpwindQuad[29]-1.0*fUpwindQuad[28]+fUpwindQuad[27]-1.0*fUpwindQuad[26]+fUpwindQuad[25]-1.0*fUpwindQuad[24]+fUpwindQuad[23]-1.0*fUpwindQuad[22]+fUpwindQuad[21]-1.0*fUpwindQuad[20]+fUpwindQuad[19]-1.0*fUpwindQuad[18]+fUpwindQuad[17]-1.0*(fUpwindQuad[16]+fUpwindQuad[15])+fUpwindQuad[14]-1.0*fUpwindQuad[13]+fUpwindQuad[12]-1.0*fUpwindQuad[11]+fUpwindQuad[10]-1.0*fUpwindQuad[9]+fUpwindQuad[8]-1.0*fUpwindQuad[7]+fUpwindQuad[6]-1.0*fUpwindQuad[5]+fUpwindQuad[4]-1.0*fUpwindQuad[3]+fUpwindQuad[2]-1.0*fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[13] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]-1.0*(fUpwindQuad[29]+fUpwindQuad[28])+fUpwindQuad[27]+fUpwindQuad[26]-1.0*(fUpwindQuad[25]+fUpwindQuad[24])+fUpwindQuad[23]+fUpwindQuad[22]-1.0*(fUpwindQuad[21]+fUpwindQuad[20])+fUpwindQuad[19]+fUpwindQuad[18]-1.0*(fUpwindQuad[17]+fUpwindQuad[16]+fUpwindQuad[15]+fUpwindQuad[14])+fUpwindQuad[13]+fUpwindQuad[12]-1.0*(fUpwindQuad[11]+fUpwindQuad[10])+fUpwindQuad[9]+fUpwindQuad[8]-1.0*(fUpwindQuad[7]+fUpwindQuad[6])+fUpwindQuad[5]+fUpwindQuad[4]-1.0*(fUpwindQuad[3]+fUpwindQuad[2])+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[14] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]+fUpwindQuad[29]+fUpwindQuad[28]-1.0*(fUpwindQuad[27]+fUpwindQuad[26]+fUpwindQuad[25]+fUpwindQuad[24])+fUpwindQuad[23]+fUpwindQuad[22]+fUpwindQuad[21]+fUpwindQuad[20]-1.0*(fUpwindQuad[19]+fUpwindQuad[18]+fUpwindQuad[17]+fUpwindQuad[16]+fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12])+fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8]-1.0*(fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4])+fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[15] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]+fUpwindQuad[29]+fUpwindQuad[28]+fUpwindQuad[27]+fUpwindQuad[26]+fUpwindQuad[25]+fUpwindQuad[24]-1.0*(fUpwindQuad[23]+fUpwindQuad[22]+fUpwindQuad[21]+fUpwindQuad[20]+fUpwindQuad[19]+fUpwindQuad[18]+fUpwindQuad[17]+fUpwindQuad[16]+fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8])+fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[16] = 0.1767766952966368*(fUpwindQuad[31]-1.0*(fUpwindQuad[30]+fUpwindQuad[29])+fUpwindQuad[28]-1.0*fUpwindQuad[27]+fUpwindQuad[26]+fUpwindQuad[25]-1.0*fUpwindQuad[24]+fUpwindQuad[23]-1.0*(fUpwindQuad[22]+fUpwindQuad[21])+fUpwindQuad[20]-1.0*fUpwindQuad[19]+fUpwindQuad[18]+fUpwindQuad[17]-1.0*fUpwindQuad[16]+fUpwindQuad[15]-1.0*(fUpwindQuad[14]+fUpwindQuad[13])+fUpwindQuad[12]-1.0*fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]-1.0*fUpwindQuad[8]+fUpwindQuad[7]-1.0*(fUpwindQuad[6]+fUpwindQuad[5])+fUpwindQuad[4]-1.0*fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[17] = 0.1767766952966368*(fUpwindQuad[31]-1.0*(fUpwindQuad[30]+fUpwindQuad[29])+fUpwindQuad[28]+fUpwindQuad[27]-1.0*(fUpwindQuad[26]+fUpwindQuad[25])+fUpwindQuad[24]-1.0*fUpwindQuad[23]+fUpwindQuad[22]+fUpwindQuad[21]-1.0*(fUpwindQuad[20]+fUpwindQuad[19])+fUpwindQuad[18]+fUpwindQuad[17]-1.0*fUpwindQuad[16]+fUpwindQuad[15]-1.0*(fUpwindQuad[14]+fUpwindQuad[13])+fUpwindQuad[12]+fUpwindQuad[11]-1.0*(fUpwindQuad[10]+fUpwindQuad[9])+fUpwindQuad[8]-1.0*fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]-1.0*(fUpwindQuad[4]+fUpwindQuad[3])+fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[18] = 0.1767766952966368*(fUpwindQuad[31]-1.0*fUpwindQuad[30]+fUpwindQuad[29]-1.0*(fUpwindQuad[28]+fUpwindQuad[27])+fUpwindQuad[26]-1.0*fUpwindQuad[25]+fUpwindQuad[24]-1.0*fUpwindQuad[23]+fUpwindQuad[22]-1.0*fUpwindQuad[21]+fUpwindQuad[20]+fUpwindQuad[19]-1.0*fUpwindQuad[18]+fUpwindQuad[17]-1.0*fUpwindQuad[16]+fUpwindQuad[15]-1.0*fUpwindQuad[14]+fUpwindQuad[13]-1.0*(fUpwindQuad[12]+fUpwindQuad[11])+fUpwindQuad[10]-1.0*fUpwindQuad[9]+fUpwindQuad[8]-1.0*fUpwindQuad[7]+fUpwindQuad[6]-1.0*fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]-1.0*fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[19] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]-1.0*(fUpwindQuad[29]+fUpwindQuad[28]+fUpwindQuad[27]+fUpwindQuad[26])+fUpwindQuad[25]+fUpwindQuad[24]-1.0*(fUpwindQuad[23]+fUpwindQuad[22])+fUpwindQuad[21]+fUpwindQuad[20]+fUpwindQuad[19]+fUpwindQuad[18]-1.0*(fUpwindQuad[17]+fUpwindQuad[16])+fUpwindQuad[15]+fUpwindQuad[14]-1.0*(fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10])+fUpwindQuad[9]+fUpwindQuad[8]-1.0*(fUpwindQuad[7]+fUpwindQuad[6])+fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2]-1.0*(fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[20] = 0.1767766952966368*(fUpwindQuad[31]-1.0*(fUpwindQuad[30]+fUpwindQuad[29])+fUpwindQuad[28]+fUpwindQuad[27]-1.0*(fUpwindQuad[26]+fUpwindQuad[25])+fUpwindQuad[24]+fUpwindQuad[23]-1.0*(fUpwindQuad[22]+fUpwindQuad[21])+fUpwindQuad[20]+fUpwindQuad[19]-1.0*(fUpwindQuad[18]+fUpwindQuad[17])+fUpwindQuad[16]-1.0*fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]-1.0*(fUpwindQuad[12]+fUpwindQuad[11])+fUpwindQuad[10]+fUpwindQuad[9]-1.0*(fUpwindQuad[8]+fUpwindQuad[7])+fUpwindQuad[6]+fUpwindQuad[5]-1.0*(fUpwindQuad[4]+fUpwindQuad[3])+fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[21] = 0.1767766952966368*(fUpwindQuad[31]-1.0*fUpwindQuad[30]+fUpwindQuad[29]-1.0*(fUpwindQuad[28]+fUpwindQuad[27])+fUpwindQuad[26]-1.0*fUpwindQuad[25]+fUpwindQuad[24]+fUpwindQuad[23]-1.0*fUpwindQuad[22]+fUpwindQuad[21]-1.0*(fUpwindQuad[20]+fUpwindQuad[19])+fUpwindQuad[18]-1.0*fUpwindQuad[17]+fUpwindQuad[16]-1.0*fUpwindQuad[15]+fUpwindQuad[14]-1.0*fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]-1.0*fUpwindQuad[10]+fUpwindQuad[9]-1.0*(fUpwindQuad[8]+fUpwindQuad[7])+fUpwindQuad[6]-1.0*fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]-1.0*fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[22] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]-1.0*(fUpwindQuad[29]+fUpwindQuad[28]+fUpwindQuad[27]+fUpwindQuad[26])+fUpwindQuad[25]+fUpwindQuad[24]+fUpwindQuad[23]+fUpwindQuad[22]-1.0*(fUpwindQuad[21]+fUpwindQuad[20]+fUpwindQuad[19]+fUpwindQuad[18])+fUpwindQuad[17]+fUpwindQuad[16]-1.0*(fUpwindQuad[15]+fUpwindQuad[14])+fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10]-1.0*(fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6])+fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2]-1.0*(fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[23] = 0.1767766952966368*(fUpwindQuad[31]-1.0*fUpwindQuad[30]+fUpwindQuad[29]-1.0*fUpwindQuad[28]+fUpwindQuad[27]-1.0*fUpwindQuad[26]+fUpwindQuad[25]-1.0*(fUpwindQuad[24]+fUpwindQuad[23])+fUpwindQuad[22]-1.0*fUpwindQuad[21]+fUpwindQuad[20]-1.0*fUpwindQuad[19]+fUpwindQuad[18]-1.0*fUpwindQuad[17]+fUpwindQuad[16]-1.0*fUpwindQuad[15]+fUpwindQuad[14]-1.0*fUpwindQuad[13]+fUpwindQuad[12]-1.0*fUpwindQuad[11]+fUpwindQuad[10]-1.0*fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]-1.0*fUpwindQuad[6]+fUpwindQuad[5]-1.0*fUpwindQuad[4]+fUpwindQuad[3]-1.0*fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[24] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]-1.0*(fUpwindQuad[29]+fUpwindQuad[28])+fUpwindQuad[27]+fUpwindQuad[26]-1.0*(fUpwindQuad[25]+fUpwindQuad[24]+fUpwindQuad[23]+fUpwindQuad[22])+fUpwindQuad[21]+fUpwindQuad[20]-1.0*(fUpwindQuad[19]+fUpwindQuad[18])+fUpwindQuad[17]+fUpwindQuad[16]-1.0*(fUpwindQuad[15]+fUpwindQuad[14])+fUpwindQuad[13]+fUpwindQuad[12]-1.0*(fUpwindQuad[11]+fUpwindQuad[10])+fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6]-1.0*(fUpwindQuad[5]+fUpwindQuad[4])+fUpwindQuad[3]+fUpwindQuad[2]-1.0*(fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[25] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]+fUpwindQuad[29]+fUpwindQuad[28]-1.0*(fUpwindQuad[27]+fUpwindQuad[26]+fUpwindQuad[25]+fUpwindQuad[24]+fUpwindQuad[23]+fUpwindQuad[22]+fUpwindQuad[21]+fUpwindQuad[20])+fUpwindQuad[19]+fUpwindQuad[18]+fUpwindQuad[17]+fUpwindQuad[16]-1.0*(fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12])+fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]-1.0*(fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[26] = 0.1767766952966368*(fUpwindQuad[31]-1.0*(fUpwindQuad[30]+fUpwindQuad[29])+fUpwindQuad[28]-1.0*fUpwindQuad[27]+fUpwindQuad[26]+fUpwindQuad[25]-1.0*(fUpwindQuad[24]+fUpwindQuad[23])+fUpwindQuad[22]+fUpwindQuad[21]-1.0*fUpwindQuad[20]+fUpwindQuad[19]-1.0*(fUpwindQuad[18]+fUpwindQuad[17])+fUpwindQuad[16]+fUpwindQuad[15]-1.0*(fUpwindQuad[14]+fUpwindQuad[13])+fUpwindQuad[12]-1.0*fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]-1.0*(fUpwindQuad[8]+fUpwindQuad[7])+fUpwindQuad[6]+fUpwindQuad[5]-1.0*fUpwindQuad[4]+fUpwindQuad[3]-1.0*(fUpwindQuad[2]+fUpwindQuad[1])+fUpwindQuad[0]); 
  fUpwind[27] = 0.1767766952966368*(fUpwindQuad[31]-1.0*(fUpwindQuad[30]+fUpwindQuad[29])+fUpwindQuad[28]-1.0*fUpwindQuad[27]+fUpwindQuad[26]+fUpwindQuad[25]-1.0*fUpwindQuad[24]+fUpwindQuad[23]-1.0*(fUpwindQuad[22]+fUpwindQuad[21])+fUpwindQuad[20]-1.0*fUpwindQuad[19]+fUpwindQuad[18]+fUpwindQuad[17]-1.0*(fUpwindQuad[16]+fUpwindQuad[15])+fUpwindQuad[14]+fUpwindQuad[13]-1.0*fUpwindQuad[12]+fUpwindQuad[11]-1.0*(fUpwindQuad[10]+fUpwindQuad[9])+fUpwindQuad[8]-1.0*fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]-1.0*fUpwindQuad[4]+fUpwindQuad[3]-1.0*(fUpwindQuad[2]+fUpwindQuad[1])+fUpwindQuad[0]); 
  fUpwind[28] = 0.1767766952966368*(fUpwindQuad[31]-1.0*(fUpwindQuad[30]+fUpwindQuad[29])+fUpwindQuad[28]+fUpwindQuad[27]-1.0*(fUpwindQuad[26]+fUpwindQuad[25])+fUpwindQuad[24]-1.0*fUpwindQuad[23]+fUpwindQuad[22]+fUpwindQuad[21]-1.0*(fUpwindQuad[20]+fUpwindQuad[19])+fUpwindQuad[18]+fUpwindQuad[17]-1.0*(fUpwindQuad[16]+fUpwindQuad[15])+fUpwindQuad[14]+fUpwindQuad[13]-1.0*(fUpwindQuad[12]+fUpwindQuad[11])+fUpwindQuad[10]+fUpwindQuad[9]-1.0*fUpwindQuad[8]+fUpwindQuad[7]-1.0*(fUpwindQuad[6]+fUpwindQuad[5])+fUpwindQuad[4]+fUpwindQuad[3]-1.0*(fUpwindQuad[2]+fUpwindQuad[1])+fUpwindQuad[0]); 
  fUpwind[29] = 0.1767766952966368*(fUpwindQuad[31]-1.0*fUpwindQuad[30]+fUpwindQuad[29]-1.0*(fUpwindQuad[28]+fUpwindQuad[27])+fUpwindQuad[26]-1.0*fUpwindQuad[25]+fUpwindQuad[24]-1.0*fUpwindQuad[23]+fUpwindQuad[22]-1.0*fUpwindQuad[21]+fUpwindQuad[20]+fUpwindQuad[19]-1.0*fUpwindQuad[18]+fUpwindQuad[17]-1.0*(fUpwindQuad[16]+fUpwindQuad[15])+fUpwindQuad[14]-1.0*fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]-1.0*fUpwindQuad[10]+fUpwindQuad[9]-1.0*fUpwindQuad[8]+fUpwindQuad[7]-1.0*fUpwindQuad[6]+fUpwindQuad[5]-1.0*(fUpwindQuad[4]+fUpwindQuad[3])+fUpwindQuad[2]-1.0*fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[30] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]-1.0*(fUpwindQuad[29]+fUpwindQuad[28]+fUpwindQuad[27]+fUpwindQuad[26])+fUpwindQuad[25]+fUpwindQuad[24]-1.0*(fUpwindQuad[23]+fUpwindQuad[22])+fUpwindQuad[21]+fUpwindQuad[20]+fUpwindQuad[19]+fUpwindQuad[18]-1.0*(fUpwindQuad[17]+fUpwindQuad[16]+fUpwindQuad[15]+fUpwindQuad[14])+fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10]-1.0*(fUpwindQuad[9]+fUpwindQuad[8])+fUpwindQuad[7]+fUpwindQuad[6]-1.0*(fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2])+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[31] = 0.1767766952966368*(fUpwindQuad[31]-1.0*(fUpwindQuad[30]+fUpwindQuad[29])+fUpwindQuad[28]-1.0*fUpwindQuad[27]+fUpwindQuad[26]+fUpwindQuad[25]-1.0*(fUpwindQuad[24]+fUpwindQuad[23])+fUpwindQuad[22]+fUpwindQuad[21]-1.0*fUpwindQuad[20]+fUpwindQuad[19]-1.0*(fUpwindQuad[18]+fUpwindQuad[17])+fUpwindQuad[16]-1.0*fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]-1.0*fUpwindQuad[12]+fUpwindQuad[11]-1.0*(fUpwindQuad[10]+fUpwindQuad[9])+fUpwindQuad[8]+fUpwindQuad[7]-1.0*(fUpwindQuad[6]+fUpwindQuad[5])+fUpwindQuad[4]-1.0*fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 

  Ghat[0] += 0.1767766952966368*(alpha[21]*fUpwind[21]+alpha[20]*fUpwind[20]+alpha[19]*fUpwind[19]+alpha[18]*fUpwind[18]+alpha[14]*fUpwind[14]+alpha[13]*fUpwind[13]+alpha[12]*fUpwind[12]+alpha[11]*fUpwind[11]+alpha[10]*fUpwind[10]+alpha[9]*fUpwind[9]+alpha[8]*fUpwind[8]+alpha[7]*fUpwind[7]+alpha[6]*fUpwind[6]+alpha[5]*fUpwind[5]+alpha[4]*fUpwind[4]+alpha[3]*fUpwind[3]+alpha[2]*fUpwind[2]+alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] += 0.1767766952966368*(alpha[19]*fUpwind[26]+alpha[14]*fUpwind[21]+fUpwind[14]*alpha[21]+alpha[13]*fUpwind[20]+fUpwind[13]*alpha[20]+alpha[11]*fUpwind[18]+fUpwind[11]*alpha[18]+alpha[10]*fUpwind[17]+alpha[8]*fUpwind[16]+alpha[5]*fUpwind[12]+fUpwind[5]*alpha[12]+alpha[4]*fUpwind[9]+fUpwind[4]*alpha[9]+alpha[3]*fUpwind[7]+fUpwind[3]*alpha[7]+alpha[2]*fUpwind[6]+fUpwind[2]*alpha[6]+alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]); 
  Ghat[2] += 0.1767766952966368*(alpha[21]*fUpwind[27]+alpha[18]*fUpwind[26]+alpha[14]*fUpwind[22]+alpha[12]*fUpwind[20]+fUpwind[12]*alpha[20]+alpha[11]*fUpwind[19]+fUpwind[11]*alpha[19]+alpha[9]*fUpwind[17]+alpha[7]*fUpwind[16]+alpha[5]*fUpwind[13]+fUpwind[5]*alpha[13]+alpha[4]*fUpwind[10]+fUpwind[4]*alpha[10]+alpha[3]*fUpwind[8]+fUpwind[3]*alpha[8]+alpha[1]*fUpwind[6]+fUpwind[1]*alpha[6]+alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2]); 
  Ghat[3] += 0.1767766952966368*(alpha[20]*fUpwind[27]+alpha[13]*fUpwind[22]+alpha[12]*fUpwind[21]+fUpwind[12]*alpha[21]+alpha[10]*fUpwind[19]+fUpwind[10]*alpha[19]+alpha[9]*fUpwind[18]+fUpwind[9]*alpha[18]+alpha[6]*fUpwind[16]+alpha[5]*fUpwind[14]+fUpwind[5]*alpha[14]+alpha[4]*fUpwind[11]+fUpwind[4]*alpha[11]+alpha[2]*fUpwind[8]+fUpwind[2]*alpha[8]+alpha[1]*fUpwind[7]+fUpwind[1]*alpha[7]+alpha[0]*fUpwind[3]+fUpwind[0]*alpha[3]); 
  Ghat[4] += 0.1767766952966368*(alpha[21]*fUpwind[29]+alpha[20]*fUpwind[28]+alpha[14]*fUpwind[25]+alpha[13]*fUpwind[24]+alpha[12]*fUpwind[23]+alpha[8]*fUpwind[19]+fUpwind[8]*alpha[19]+alpha[7]*fUpwind[18]+fUpwind[7]*alpha[18]+alpha[6]*fUpwind[17]+alpha[5]*fUpwind[15]+alpha[3]*fUpwind[11]+fUpwind[3]*alpha[11]+alpha[2]*fUpwind[10]+fUpwind[2]*alpha[10]+alpha[1]*fUpwind[9]+fUpwind[1]*alpha[9]+alpha[0]*fUpwind[4]+fUpwind[0]*alpha[4]); 
  Ghat[5] += 0.1767766952966368*(alpha[19]*fUpwind[30]+alpha[18]*fUpwind[29]+alpha[11]*fUpwind[25]+alpha[10]*fUpwind[24]+alpha[9]*fUpwind[23]+alpha[8]*fUpwind[22]+alpha[7]*fUpwind[21]+fUpwind[7]*alpha[21]+alpha[6]*fUpwind[20]+fUpwind[6]*alpha[20]+alpha[4]*fUpwind[15]+alpha[3]*fUpwind[14]+fUpwind[3]*alpha[14]+alpha[2]*fUpwind[13]+fUpwind[2]*alpha[13]+alpha[1]*fUpwind[12]+fUpwind[1]*alpha[12]+alpha[0]*fUpwind[5]+fUpwind[0]*alpha[5]); 
  Ghat[6] += 0.1767766952966368*(alpha[14]*fUpwind[27]+alpha[11]*fUpwind[26]+alpha[21]*fUpwind[22]+alpha[5]*fUpwind[20]+fUpwind[5]*alpha[20]+alpha[18]*fUpwind[19]+fUpwind[18]*alpha[19]+alpha[4]*fUpwind[17]+alpha[3]*fUpwind[16]+alpha[12]*fUpwind[13]+fUpwind[12]*alpha[13]+alpha[9]*fUpwind[10]+fUpwind[9]*alpha[10]+alpha[7]*fUpwind[8]+fUpwind[7]*alpha[8]+alpha[0]*fUpwind[6]+fUpwind[0]*alpha[6]+alpha[1]*fUpwind[2]+fUpwind[1]*alpha[2]); 
  Ghat[7] += 0.1767766952966368*(alpha[13]*fUpwind[27]+alpha[10]*fUpwind[26]+alpha[20]*fUpwind[22]+alpha[5]*fUpwind[21]+fUpwind[5]*alpha[21]+fUpwind[17]*alpha[19]+alpha[4]*fUpwind[18]+fUpwind[4]*alpha[18]+alpha[2]*fUpwind[16]+alpha[12]*fUpwind[14]+fUpwind[12]*alpha[14]+alpha[9]*fUpwind[11]+fUpwind[9]*alpha[11]+alpha[6]*fUpwind[8]+fUpwind[6]*alpha[8]+alpha[0]*fUpwind[7]+fUpwind[0]*alpha[7]+alpha[1]*fUpwind[3]+fUpwind[1]*alpha[3]); 
  Ghat[8] += 0.1767766952966368*(alpha[12]*fUpwind[27]+alpha[9]*fUpwind[26]+alpha[5]*fUpwind[22]+alpha[20]*fUpwind[21]+fUpwind[20]*alpha[21]+alpha[4]*fUpwind[19]+fUpwind[4]*alpha[19]+fUpwind[17]*alpha[18]+alpha[1]*fUpwind[16]+alpha[13]*fUpwind[14]+fUpwind[13]*alpha[14]+alpha[10]*fUpwind[11]+fUpwind[10]*alpha[11]+alpha[0]*fUpwind[8]+fUpwind[0]*alpha[8]+alpha[6]*fUpwind[7]+fUpwind[6]*alpha[7]+alpha[2]*fUpwind[3]+fUpwind[2]*alpha[3]); 
  Ghat[9] += 0.1767766952966368*(alpha[14]*fUpwind[29]+alpha[13]*fUpwind[28]+alpha[8]*fUpwind[26]+alpha[21]*fUpwind[25]+alpha[20]*fUpwind[24]+alpha[5]*fUpwind[23]+fUpwind[16]*alpha[19]+alpha[3]*fUpwind[18]+fUpwind[3]*alpha[18]+alpha[2]*fUpwind[17]+alpha[12]*fUpwind[15]+alpha[7]*fUpwind[11]+fUpwind[7]*alpha[11]+alpha[6]*fUpwind[10]+fUpwind[6]*alpha[10]+alpha[0]*fUpwind[9]+fUpwind[0]*alpha[9]+alpha[1]*fUpwind[4]+fUpwind[1]*alpha[4]); 
  Ghat[10] += 0.1767766952966368*(alpha[21]*fUpwind[31]+alpha[14]*fUpwind[30]+alpha[12]*fUpwind[28]+alpha[7]*fUpwind[26]+alpha[5]*fUpwind[24]+alpha[20]*fUpwind[23]+alpha[3]*fUpwind[19]+fUpwind[3]*alpha[19]+fUpwind[16]*alpha[18]+alpha[1]*fUpwind[17]+alpha[13]*fUpwind[15]+alpha[8]*fUpwind[11]+fUpwind[8]*alpha[11]+alpha[0]*fUpwind[10]+fUpwind[0]*alpha[10]+alpha[6]*fUpwind[9]+fUpwind[6]*alpha[9]+alpha[2]*fUpwind[4]+fUpwind[2]*alpha[4]); 
  Ghat[11] += 0.1767766952966368*(alpha[20]*fUpwind[31]+alpha[13]*fUpwind[30]+alpha[12]*fUpwind[29]+alpha[6]*fUpwind[26]+alpha[5]*fUpwind[25]+alpha[21]*fUpwind[23]+alpha[2]*fUpwind[19]+fUpwind[2]*alpha[19]+alpha[1]*fUpwind[18]+fUpwind[1]*alpha[18]+alpha[14]*fUpwind[15]+alpha[0]*fUpwind[11]+fUpwind[0]*alpha[11]+alpha[8]*fUpwind[10]+fUpwind[8]*alpha[10]+alpha[7]*fUpwind[9]+fUpwind[7]*alpha[9]+alpha[3]*fUpwind[4]+fUpwind[3]*alpha[4]); 
  Ghat[12] += 0.1767766952966368*(alpha[19]*fUpwind[31]+alpha[11]*fUpwind[29]+alpha[10]*fUpwind[28]+alpha[8]*fUpwind[27]+alpha[18]*fUpwind[25]+alpha[4]*fUpwind[23]+alpha[3]*fUpwind[21]+fUpwind[3]*alpha[21]+alpha[2]*fUpwind[20]+fUpwind[2]*alpha[20]+alpha[9]*fUpwind[15]+alpha[7]*fUpwind[14]+fUpwind[7]*alpha[14]+alpha[6]*fUpwind[13]+fUpwind[6]*alpha[13]+alpha[0]*fUpwind[12]+fUpwind[0]*alpha[12]+alpha[1]*fUpwind[5]+fUpwind[1]*alpha[5]); 
  Ghat[13] += 0.1767766952966368*(alpha[18]*fUpwind[31]+alpha[11]*fUpwind[30]+alpha[9]*fUpwind[28]+alpha[7]*fUpwind[27]+alpha[19]*fUpwind[25]+alpha[4]*fUpwind[24]+alpha[3]*fUpwind[22]+fUpwind[16]*alpha[21]+alpha[1]*fUpwind[20]+fUpwind[1]*alpha[20]+alpha[10]*fUpwind[15]+alpha[8]*fUpwind[14]+fUpwind[8]*alpha[14]+alpha[0]*fUpwind[13]+fUpwind[0]*alpha[13]+alpha[6]*fUpwind[12]+fUpwind[6]*alpha[12]+alpha[2]*fUpwind[5]+fUpwind[2]*alpha[5]); 
  Ghat[14] += 0.1767766952966368*(alpha[10]*fUpwind[30]+alpha[9]*fUpwind[29]+alpha[6]*fUpwind[27]+alpha[4]*fUpwind[25]+alpha[19]*fUpwind[24]+alpha[18]*fUpwind[23]+alpha[2]*fUpwind[22]+alpha[1]*fUpwind[21]+fUpwind[1]*alpha[21]+fUpwind[16]*alpha[20]+alpha[11]*fUpwind[15]+alpha[0]*fUpwind[14]+fUpwind[0]*alpha[14]+alpha[8]*fUpwind[13]+fUpwind[8]*alpha[13]+alpha[7]*fUpwind[12]+fUpwind[7]*alpha[12]+alpha[3]*fUpwind[5]+fUpwind[3]*alpha[5]); 
  Ghat[15] += 0.1767766952966368*(alpha[8]*fUpwind[30]+alpha[7]*fUpwind[29]+alpha[6]*fUpwind[28]+alpha[3]*fUpwind[25]+alpha[2]*fUpwind[24]+alpha[1]*fUpwind[23]+alpha[19]*fUpwind[22]+alpha[18]*fUpwind[21]+fUpwind[18]*alpha[21]+fUpwind[17]*alpha[20]+alpha[0]*fUpwind[15]+alpha[11]*fUpwind[14]+fUpwind[11]*alpha[14]+alpha[10]*fUpwind[13]+fUpwind[10]*alpha[13]+alpha[9]*fUpwind[12]+fUpwind[9]*alpha[12]+alpha[4]*fUpwind[5]+fUpwind[4]*alpha[5]); 
  Ghat[16] += 0.1767766952966368*(alpha[5]*fUpwind[27]+alpha[4]*fUpwind[26]+alpha[12]*fUpwind[22]+alpha[13]*fUpwind[21]+fUpwind[13]*alpha[21]+alpha[14]*fUpwind[20]+fUpwind[14]*alpha[20]+alpha[9]*fUpwind[19]+fUpwind[9]*alpha[19]+alpha[10]*fUpwind[18]+fUpwind[10]*alpha[18]+alpha[11]*fUpwind[17]+alpha[0]*fUpwind[16]+alpha[1]*fUpwind[8]+fUpwind[1]*alpha[8]+alpha[2]*fUpwind[7]+fUpwind[2]*alpha[7]+alpha[3]*fUpwind[6]+fUpwind[3]*alpha[6]); 
  Ghat[17] += 0.1767766952966368*(alpha[14]*fUpwind[31]+alpha[21]*fUpwind[30]+alpha[5]*fUpwind[28]+alpha[3]*fUpwind[26]+alpha[12]*fUpwind[24]+alpha[13]*fUpwind[23]+fUpwind[15]*alpha[20]+alpha[7]*fUpwind[19]+fUpwind[7]*alpha[19]+alpha[8]*fUpwind[18]+fUpwind[8]*alpha[18]+alpha[0]*fUpwind[17]+alpha[11]*fUpwind[16]+alpha[1]*fUpwind[10]+fUpwind[1]*alpha[10]+alpha[2]*fUpwind[9]+fUpwind[2]*alpha[9]+alpha[4]*fUpwind[6]+fUpwind[4]*alpha[6]); 
  Ghat[18] += 0.1767766952966368*(alpha[13]*fUpwind[31]+alpha[20]*fUpwind[30]+alpha[5]*fUpwind[29]+alpha[2]*fUpwind[26]+alpha[12]*fUpwind[25]+alpha[14]*fUpwind[23]+fUpwind[15]*alpha[21]+alpha[6]*fUpwind[19]+fUpwind[6]*alpha[19]+alpha[0]*fUpwind[18]+fUpwind[0]*alpha[18]+alpha[8]*fUpwind[17]+alpha[10]*fUpwind[16]+alpha[1]*fUpwind[11]+fUpwind[1]*alpha[11]+alpha[3]*fUpwind[9]+fUpwind[3]*alpha[9]+alpha[4]*fUpwind[7]+fUpwind[4]*alpha[7]); 
  Ghat[19] += 0.1767766952966368*(alpha[12]*fUpwind[31]+alpha[5]*fUpwind[30]+alpha[20]*fUpwind[29]+alpha[21]*fUpwind[28]+alpha[1]*fUpwind[26]+alpha[13]*fUpwind[25]+alpha[14]*fUpwind[24]+alpha[0]*fUpwind[19]+fUpwind[0]*alpha[19]+alpha[6]*fUpwind[18]+fUpwind[6]*alpha[18]+alpha[7]*fUpwind[17]+alpha[9]*fUpwind[16]+alpha[2]*fUpwind[11]+fUpwind[2]*alpha[11]+alpha[3]*fUpwind[10]+fUpwind[3]*alpha[10]+alpha[4]*fUpwind[8]+fUpwind[4]*alpha[8]); 
  Ghat[20] += 0.1767766952966368*(alpha[11]*fUpwind[31]+alpha[18]*fUpwind[30]+alpha[19]*fUpwind[29]+alpha[4]*fUpwind[28]+alpha[3]*fUpwind[27]+alpha[9]*fUpwind[24]+alpha[10]*fUpwind[23]+alpha[7]*fUpwind[22]+alpha[8]*fUpwind[21]+fUpwind[8]*alpha[21]+alpha[0]*fUpwind[20]+fUpwind[0]*alpha[20]+alpha[14]*fUpwind[16]+alpha[1]*fUpwind[13]+fUpwind[1]*alpha[13]+alpha[2]*fUpwind[12]+fUpwind[2]*alpha[12]+alpha[5]*fUpwind[6]+fUpwind[5]*alpha[6]); 
  Ghat[21] += 0.1767766952966368*(alpha[10]*fUpwind[31]+alpha[4]*fUpwind[29]+alpha[19]*fUpwind[28]+alpha[2]*fUpwind[27]+alpha[9]*fUpwind[25]+alpha[11]*fUpwind[23]+alpha[6]*fUpwind[22]+alpha[0]*fUpwind[21]+fUpwind[0]*alpha[21]+alpha[8]*fUpwind[20]+fUpwind[8]*alpha[20]+fUpwind[15]*alpha[18]+alpha[13]*fUpwind[16]+alpha[1]*fUpwind[14]+fUpwind[1]*alpha[14]+alpha[3]*fUpwind[12]+fUpwind[3]*alpha[12]+alpha[5]*fUpwind[7]+fUpwind[5]*alpha[7]); 
  Ghat[22] += 0.1767766952966368*(alpha[9]*fUpwind[31]+alpha[4]*fUpwind[30]+alpha[18]*fUpwind[28]+alpha[1]*fUpwind[27]+alpha[10]*fUpwind[25]+alpha[11]*fUpwind[24]+alpha[0]*fUpwind[22]+alpha[6]*fUpwind[21]+fUpwind[6]*alpha[21]+alpha[7]*fUpwind[20]+fUpwind[7]*alpha[20]+fUpwind[15]*alpha[19]+alpha[12]*fUpwind[16]+alpha[2]*fUpwind[14]+fUpwind[2]*alpha[14]+alpha[3]*fUpwind[13]+fUpwind[3]*alpha[13]+alpha[5]*fUpwind[8]+fUpwind[5]*alpha[8]); 
  Ghat[23] += 0.1767766952966368*(alpha[8]*fUpwind[31]+alpha[3]*fUpwind[29]+alpha[2]*fUpwind[28]+alpha[19]*fUpwind[27]+alpha[7]*fUpwind[25]+alpha[6]*fUpwind[24]+alpha[0]*fUpwind[23]+alpha[11]*fUpwind[21]+fUpwind[11]*alpha[21]+alpha[10]*fUpwind[20]+fUpwind[10]*alpha[20]+alpha[14]*fUpwind[18]+fUpwind[14]*alpha[18]+alpha[13]*fUpwind[17]+alpha[1]*fUpwind[15]+alpha[4]*fUpwind[12]+fUpwind[4]*alpha[12]+alpha[5]*fUpwind[9]+fUpwind[5]*alpha[9]); 
  Ghat[24] += 0.1767766952966368*(alpha[7]*fUpwind[31]+alpha[3]*fUpwind[30]+alpha[1]*fUpwind[28]+alpha[18]*fUpwind[27]+alpha[21]*fUpwind[26]+alpha[8]*fUpwind[25]+alpha[0]*fUpwind[24]+alpha[6]*fUpwind[23]+alpha[11]*fUpwind[22]+alpha[9]*fUpwind[20]+fUpwind[9]*alpha[20]+alpha[14]*fUpwind[19]+fUpwind[14]*alpha[19]+alpha[12]*fUpwind[17]+alpha[2]*fUpwind[15]+alpha[4]*fUpwind[13]+fUpwind[4]*alpha[13]+alpha[5]*fUpwind[10]+fUpwind[5]*alpha[10]); 
  Ghat[25] += 0.1767766952966368*(alpha[6]*fUpwind[31]+alpha[2]*fUpwind[30]+alpha[1]*fUpwind[29]+alpha[20]*fUpwind[26]+alpha[0]*fUpwind[25]+alpha[8]*fUpwind[24]+alpha[7]*fUpwind[23]+alpha[10]*fUpwind[22]+alpha[9]*fUpwind[21]+fUpwind[9]*alpha[21]+alpha[13]*fUpwind[19]+fUpwind[13]*alpha[19]+alpha[12]*fUpwind[18]+fUpwind[12]*alpha[18]+alpha[3]*fUpwind[15]+alpha[4]*fUpwind[14]+fUpwind[4]*alpha[14]+alpha[5]*fUpwind[11]+fUpwind[5]*alpha[11]); 
  Ghat[26] += 0.1767766952966368*(alpha[5]*fUpwind[31]+alpha[12]*fUpwind[30]+alpha[13]*fUpwind[29]+alpha[14]*fUpwind[28]+alpha[0]*fUpwind[26]+alpha[20]*fUpwind[25]+alpha[21]*fUpwind[24]+alpha[1]*fUpwind[19]+fUpwind[1]*alpha[19]+alpha[2]*fUpwind[18]+fUpwind[2]*alpha[18]+alpha[3]*fUpwind[17]+alpha[4]*fUpwind[16]+alpha[6]*fUpwind[11]+fUpwind[6]*alpha[11]+alpha[7]*fUpwind[10]+fUpwind[7]*alpha[10]+alpha[8]*fUpwind[9]+fUpwind[8]*alpha[9]); 
  Ghat[27] += 0.1767766952966368*(alpha[4]*fUpwind[31]+alpha[9]*fUpwind[30]+alpha[10]*fUpwind[29]+alpha[11]*fUpwind[28]+alpha[0]*fUpwind[27]+alpha[18]*fUpwind[24]+alpha[19]*fUpwind[23]+alpha[1]*fUpwind[22]+alpha[2]*fUpwind[21]+fUpwind[2]*alpha[21]+alpha[3]*fUpwind[20]+fUpwind[3]*alpha[20]+alpha[5]*fUpwind[16]+alpha[6]*fUpwind[14]+fUpwind[6]*alpha[14]+alpha[7]*fUpwind[13]+fUpwind[7]*alpha[13]+alpha[8]*fUpwind[12]+fUpwind[8]*alpha[12]); 
  Ghat[28] += 0.1767766952966368*(alpha[3]*fUpwind[31]+alpha[7]*fUpwind[30]+alpha[8]*fUpwind[29]+alpha[0]*fUpwind[28]+alpha[11]*fUpwind[27]+alpha[14]*fUpwind[26]+alpha[1]*fUpwind[24]+alpha[2]*fUpwind[23]+alpha[18]*fUpwind[22]+alpha[19]*fUpwind[21]+fUpwind[19]*alpha[21]+alpha[4]*fUpwind[20]+fUpwind[4]*alpha[20]+alpha[5]*fUpwind[17]+alpha[6]*fUpwind[15]+alpha[9]*fUpwind[13]+fUpwind[9]*alpha[13]+alpha[10]*fUpwind[12]+fUpwind[10]*alpha[12]); 
  Ghat[29] += 0.1767766952966368*(alpha[2]*fUpwind[31]+alpha[6]*fUpwind[30]+alpha[0]*fUpwind[29]+alpha[8]*fUpwind[28]+alpha[10]*fUpwind[27]+alpha[13]*fUpwind[26]+alpha[1]*fUpwind[25]+alpha[3]*fUpwind[23]+alpha[4]*fUpwind[21]+fUpwind[4]*alpha[21]+alpha[19]*fUpwind[20]+fUpwind[19]*alpha[20]+alpha[5]*fUpwind[18]+fUpwind[5]*alpha[18]+alpha[7]*fUpwind[15]+alpha[9]*fUpwind[14]+fUpwind[9]*alpha[14]+alpha[11]*fUpwind[12]+fUpwind[11]*alpha[12]); 
  Ghat[30] += 0.1767766952966368*(alpha[1]*fUpwind[31]+alpha[0]*fUpwind[30]+alpha[6]*fUpwind[29]+alpha[7]*fUpwind[28]+alpha[9]*fUpwind[27]+alpha[12]*fUpwind[26]+alpha[2]*fUpwind[25]+alpha[3]*fUpwind[24]+alpha[4]*fUpwind[22]+fUpwind[17]*alpha[21]+alpha[18]*fUpwind[20]+fUpwind[18]*alpha[20]+alpha[5]*fUpwind[19]+fUpwind[5]*alpha[19]+alpha[8]*fUpwind[15]+alpha[10]*fUpwind[14]+fUpwind[10]*alpha[14]+alpha[11]*fUpwind[13]+fUpwind[11]*alpha[13]); 
  Ghat[31] += 0.1767766952966368*(alpha[0]*fUpwind[31]+alpha[1]*fUpwind[30]+alpha[2]*fUpwind[29]+alpha[3]*fUpwind[28]+alpha[4]*fUpwind[27]+alpha[5]*fUpwind[26]+alpha[6]*fUpwind[25]+alpha[7]*fUpwind[24]+alpha[8]*fUpwind[23]+alpha[9]*fUpwind[22]+alpha[10]*fUpwind[21]+fUpwind[10]*alpha[21]+alpha[11]*fUpwind[20]+fUpwind[11]*alpha[20]+alpha[12]*fUpwind[19]+fUpwind[12]*alpha[19]+alpha[13]*fUpwind[18]+fUpwind[13]*alpha[18]+alpha[14]*fUpwind[17]); 

  out[0] += -0.7071067811865475*Ghat[0]*dv11; 
  out[1] += -0.7071067811865475*Ghat[1]*dv11; 
  out[2] += -0.7071067811865475*Ghat[2]*dv11; 
  out[3] += -0.7071067811865475*Ghat[3]*dv11; 
  out[4] += -0.7071067811865475*Ghat[4]*dv11; 
  out[5] += -1.224744871391589*Ghat[0]*dv11; 
  out[6] += -0.7071067811865475*Ghat[5]*dv11; 
  out[7] += -0.7071067811865475*Ghat[6]*dv11; 
  out[8] += -0.7071067811865475*Ghat[7]*dv11; 
  out[9] += -0.7071067811865475*Ghat[8]*dv11; 
  out[10] += -0.7071067811865475*Ghat[9]*dv11; 
  out[11] += -0.7071067811865475*Ghat[10]*dv11; 
  out[12] += -0.7071067811865475*Ghat[11]*dv11; 
  out[13] += -1.224744871391589*Ghat[1]*dv11; 
  out[14] += -1.224744871391589*Ghat[2]*dv11; 
  out[15] += -1.224744871391589*Ghat[3]*dv11; 
  out[16] += -1.224744871391589*Ghat[4]*dv11; 
  out[17] += -0.7071067811865475*Ghat[12]*dv11; 
  out[18] += -0.7071067811865475*Ghat[13]*dv11; 
  out[19] += -0.7071067811865475*Ghat[14]*dv11; 
  out[20] += -0.7071067811865475*Ghat[15]*dv11; 
  out[21] += -1.224744871391589*Ghat[5]*dv11; 
  out[22] += -0.7071067811865475*Ghat[16]*dv11; 
  out[23] += -0.7071067811865475*Ghat[17]*dv11; 
  out[24] += -0.7071067811865475*Ghat[18]*dv11; 
  out[25] += -0.7071067811865475*Ghat[19]*dv11; 
  out[26] += -1.224744871391589*Ghat[6]*dv11; 
  out[27] += -1.224744871391589*Ghat[7]*dv11; 
  out[28] += -1.224744871391589*Ghat[8]*dv11; 
  out[29] += -1.224744871391589*Ghat[9]*dv11; 
  out[30] += -1.224744871391589*Ghat[10]*dv11; 
  out[31] += -1.224744871391589*Ghat[11]*dv11; 
  out[32] += -0.7071067811865475*Ghat[20]*dv11; 
  out[33] += -0.7071067811865475*Ghat[21]*dv11; 
  out[34] += -0.7071067811865475*Ghat[22]*dv11; 
  out[35] += -0.7071067811865475*Ghat[23]*dv11; 
  out[36] += -0.7071067811865475*Ghat[24]*dv11; 
  out[37] += -0.7071067811865475*Ghat[25]*dv11; 
  out[38] += -1.224744871391589*Ghat[12]*dv11; 
  out[39] += -1.224744871391589*Ghat[13]*dv11; 
  out[40] += -1.224744871391589*Ghat[14]*dv11; 
  out[41] += -1.224744871391589*Ghat[15]*dv11; 
  out[42] += -0.7071067811865475*Ghat[26]*dv11; 
  out[43] += -1.224744871391589*Ghat[16]*dv11; 
  out[44] += -1.224744871391589*Ghat[17]*dv11; 
  out[45] += -1.224744871391589*Ghat[18]*dv11; 
  out[46] += -1.224744871391589*Ghat[19]*dv11; 
  out[47] += -0.7071067811865475*Ghat[27]*dv11; 
  out[48] += -0.7071067811865475*Ghat[28]*dv11; 
  out[49] += -0.7071067811865475*Ghat[29]*dv11; 
  out[50] += -0.7071067811865475*Ghat[30]*dv11; 
  out[51] += -1.224744871391589*Ghat[20]*dv11; 
  out[52] += -1.224744871391589*Ghat[21]*dv11; 
  out[53] += -1.224744871391589*Ghat[22]*dv11; 
  out[54] += -1.224744871391589*Ghat[23]*dv11; 
  out[55] += -1.224744871391589*Ghat[24]*dv11; 
  out[56] += -1.224744871391589*Ghat[25]*dv11; 
  out[57] += -1.224744871391589*Ghat[26]*dv11; 
  out[58] += -0.7071067811865475*Ghat[31]*dv11; 
  out[59] += -1.224744871391589*Ghat[27]*dv11; 
  out[60] += -1.224744871391589*Ghat[28]*dv11; 
  out[61] += -1.224744871391589*Ghat[29]*dv11; 
  out[62] += -1.224744871391589*Ghat[30]*dv11; 
  out[63] += -1.224744871391589*Ghat[31]*dv11; 

  } else { 

  if ((-alpha[21])-alpha[20]-alpha[19]-alpha[18]+alpha[14]+alpha[13]+alpha[12]+alpha[11]+alpha[10]+alpha[9]+alpha[8]+alpha[7]+alpha[6]-alpha[5]-alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[0] = ser_3x3v_p1_surfvy_quad_0(1, fEdge); 
  } else { 

    fUpwindQuad[0] = ser_3x3v_p1_surfvy_quad_0(-1, fSkin); 
  } 
  if (alpha[21]+alpha[20]-alpha[19]+alpha[18]+alpha[14]+alpha[13]-alpha[12]+alpha[11]+alpha[10]-alpha[9]+alpha[8]-alpha[7]-alpha[6]-alpha[5]-alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[1] = ser_3x3v_p1_surfvy_quad_1(1, fEdge); 
  } else { 

    fUpwindQuad[1] = ser_3x3v_p1_surfvy_quad_1(-1, fSkin); 
  } 
  if ((-alpha[21])+alpha[20]+alpha[19]-alpha[18]+alpha[14]-alpha[13]+alpha[12]+alpha[11]-alpha[10]+alpha[9]-alpha[8]+alpha[7]-alpha[6]-alpha[5]-alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[2] = ser_3x3v_p1_surfvy_quad_2(1, fEdge); 
  } else { 

    fUpwindQuad[2] = ser_3x3v_p1_surfvy_quad_2(-1, fSkin); 
  } 
  if (alpha[21]-alpha[20]+alpha[19]+alpha[18]+alpha[14]-alpha[13]-alpha[12]+alpha[11]-alpha[10]-alpha[9]-alpha[8]-alpha[7]+alpha[6]-alpha[5]-alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[3] = ser_3x3v_p1_surfvy_quad_3(1, fEdge); 
  } else { 

    fUpwindQuad[3] = ser_3x3v_p1_surfvy_quad_3(-1, fSkin); 
  } 
  if (alpha[21]-alpha[20]+alpha[19]+alpha[18]-alpha[14]+alpha[13]+alpha[12]-alpha[11]+alpha[10]+alpha[9]-alpha[8]-alpha[7]+alpha[6]-alpha[5]-alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[4] = ser_3x3v_p1_surfvy_quad_4(1, fEdge); 
  } else { 

    fUpwindQuad[4] = ser_3x3v_p1_surfvy_quad_4(-1, fSkin); 
  } 
  if ((-alpha[21])+alpha[20]+alpha[19]-alpha[18]-alpha[14]+alpha[13]-alpha[12]-alpha[11]+alpha[10]-alpha[9]-alpha[8]+alpha[7]-alpha[6]-alpha[5]-alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[5] = ser_3x3v_p1_surfvy_quad_5(1, fEdge); 
  } else { 

    fUpwindQuad[5] = ser_3x3v_p1_surfvy_quad_5(-1, fSkin); 
  } 
  if (alpha[21]+alpha[20]-alpha[19]+alpha[18]-alpha[14]-alpha[13]+alpha[12]-alpha[11]-alpha[10]+alpha[9]+alpha[8]-alpha[7]-alpha[6]-alpha[5]-alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[6] = ser_3x3v_p1_surfvy_quad_6(1, fEdge); 
  } else { 

    fUpwindQuad[6] = ser_3x3v_p1_surfvy_quad_6(-1, fSkin); 
  } 
  if ((-alpha[21])-alpha[20]-alpha[19]-alpha[18]-alpha[14]-alpha[13]-alpha[12]-alpha[11]-alpha[10]-alpha[9]+alpha[8]+alpha[7]+alpha[6]-alpha[5]-alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[7] = ser_3x3v_p1_surfvy_quad_7(1, fEdge); 
  } else { 

    fUpwindQuad[7] = ser_3x3v_p1_surfvy_quad_7(-1, fSkin); 
  } 
  if ((-alpha[21])-alpha[20]+alpha[19]+alpha[18]+alpha[14]+alpha[13]+alpha[12]-alpha[11]-alpha[10]-alpha[9]+alpha[8]+alpha[7]+alpha[6]-alpha[5]+alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[8] = ser_3x3v_p1_surfvy_quad_8(1, fEdge); 
  } else { 

    fUpwindQuad[8] = ser_3x3v_p1_surfvy_quad_8(-1, fSkin); 
  } 
  if (alpha[21]+alpha[20]+alpha[19]-alpha[18]+alpha[14]+alpha[13]-alpha[12]-alpha[11]-alpha[10]+alpha[9]+alpha[8]-alpha[7]-alpha[6]-alpha[5]+alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[9] = ser_3x3v_p1_surfvy_quad_9(1, fEdge); 
  } else { 

    fUpwindQuad[9] = ser_3x3v_p1_surfvy_quad_9(-1, fSkin); 
  } 
  if ((-alpha[21])+alpha[20]-alpha[19]+alpha[18]+alpha[14]-alpha[13]+alpha[12]-alpha[11]+alpha[10]-alpha[9]-alpha[8]+alpha[7]-alpha[6]-alpha[5]+alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[10] = ser_3x3v_p1_surfvy_quad_10(1, fEdge); 
  } else { 

    fUpwindQuad[10] = ser_3x3v_p1_surfvy_quad_10(-1, fSkin); 
  } 
  if (alpha[21]-alpha[20]-alpha[19]-alpha[18]+alpha[14]-alpha[13]-alpha[12]-alpha[11]+alpha[10]+alpha[9]-alpha[8]-alpha[7]+alpha[6]-alpha[5]+alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[11] = ser_3x3v_p1_surfvy_quad_11(1, fEdge); 
  } else { 

    fUpwindQuad[11] = ser_3x3v_p1_surfvy_quad_11(-1, fSkin); 
  } 
  if (alpha[21]-alpha[20]-alpha[19]-alpha[18]-alpha[14]+alpha[13]+alpha[12]+alpha[11]-alpha[10]-alpha[9]-alpha[8]-alpha[7]+alpha[6]-alpha[5]+alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[12] = ser_3x3v_p1_surfvy_quad_12(1, fEdge); 
  } else { 

    fUpwindQuad[12] = ser_3x3v_p1_surfvy_quad_12(-1, fSkin); 
  } 
  if ((-alpha[21])+alpha[20]-alpha[19]+alpha[18]-alpha[14]+alpha[13]-alpha[12]+alpha[11]-alpha[10]+alpha[9]-alpha[8]+alpha[7]-alpha[6]-alpha[5]+alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[13] = ser_3x3v_p1_surfvy_quad_13(1, fEdge); 
  } else { 

    fUpwindQuad[13] = ser_3x3v_p1_surfvy_quad_13(-1, fSkin); 
  } 
  if (alpha[21]+alpha[20]+alpha[19]-alpha[18]-alpha[14]-alpha[13]+alpha[12]+alpha[11]+alpha[10]-alpha[9]+alpha[8]-alpha[7]-alpha[6]-alpha[5]+alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[14] = ser_3x3v_p1_surfvy_quad_14(1, fEdge); 
  } else { 

    fUpwindQuad[14] = ser_3x3v_p1_surfvy_quad_14(-1, fSkin); 
  } 
  if ((-alpha[21])-alpha[20]+alpha[19]+alpha[18]-alpha[14]-alpha[13]-alpha[12]+alpha[11]+alpha[10]+alpha[9]+alpha[8]+alpha[7]+alpha[6]-alpha[5]+alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[15] = ser_3x3v_p1_surfvy_quad_15(1, fEdge); 
  } else { 

    fUpwindQuad[15] = ser_3x3v_p1_surfvy_quad_15(-1, fSkin); 
  } 
  if (alpha[21]+alpha[20]-alpha[19]-alpha[18]-alpha[14]-alpha[13]-alpha[12]+alpha[11]+alpha[10]+alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5]-alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[16] = ser_3x3v_p1_surfvy_quad_16(1, fEdge); 
  } else { 

    fUpwindQuad[16] = ser_3x3v_p1_surfvy_quad_16(-1, fSkin); 
  } 
  if ((-alpha[21])-alpha[20]-alpha[19]+alpha[18]-alpha[14]-alpha[13]+alpha[12]+alpha[11]+alpha[10]-alpha[9]+alpha[8]-alpha[7]-alpha[6]+alpha[5]-alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[17] = ser_3x3v_p1_surfvy_quad_17(1, fEdge); 
  } else { 

    fUpwindQuad[17] = ser_3x3v_p1_surfvy_quad_17(-1, fSkin); 
  } 
  if (alpha[21]-alpha[20]+alpha[19]-alpha[18]-alpha[14]+alpha[13]-alpha[12]+alpha[11]-alpha[10]+alpha[9]-alpha[8]+alpha[7]-alpha[6]+alpha[5]-alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[18] = ser_3x3v_p1_surfvy_quad_18(1, fEdge); 
  } else { 

    fUpwindQuad[18] = ser_3x3v_p1_surfvy_quad_18(-1, fSkin); 
  } 
  if ((-alpha[21])+alpha[20]+alpha[19]+alpha[18]-alpha[14]+alpha[13]+alpha[12]+alpha[11]-alpha[10]-alpha[9]-alpha[8]-alpha[7]+alpha[6]+alpha[5]-alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[19] = ser_3x3v_p1_surfvy_quad_19(1, fEdge); 
  } else { 

    fUpwindQuad[19] = ser_3x3v_p1_surfvy_quad_19(-1, fSkin); 
  } 
  if ((-alpha[21])+alpha[20]+alpha[19]+alpha[18]+alpha[14]-alpha[13]-alpha[12]-alpha[11]+alpha[10]+alpha[9]-alpha[8]-alpha[7]+alpha[6]+alpha[5]-alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[20] = ser_3x3v_p1_surfvy_quad_20(1, fEdge); 
  } else { 

    fUpwindQuad[20] = ser_3x3v_p1_surfvy_quad_20(-1, fSkin); 
  } 
  if (alpha[21]-alpha[20]+alpha[19]-alpha[18]+alpha[14]-alpha[13]+alpha[12]-alpha[11]+alpha[10]-alpha[9]-alpha[8]+alpha[7]-alpha[6]+alpha[5]-alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[21] = ser_3x3v_p1_surfvy_quad_21(1, fEdge); 
  } else { 

    fUpwindQuad[21] = ser_3x3v_p1_surfvy_quad_21(-1, fSkin); 
  } 
  if ((-alpha[21])-alpha[20]-alpha[19]+alpha[18]+alpha[14]+alpha[13]-alpha[12]-alpha[11]-alpha[10]+alpha[9]+alpha[8]-alpha[7]-alpha[6]+alpha[5]-alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[22] = ser_3x3v_p1_surfvy_quad_22(1, fEdge); 
  } else { 

    fUpwindQuad[22] = ser_3x3v_p1_surfvy_quad_22(-1, fSkin); 
  } 
  if (alpha[21]+alpha[20]-alpha[19]-alpha[18]+alpha[14]+alpha[13]+alpha[12]-alpha[11]-alpha[10]-alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5]-alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[23] = ser_3x3v_p1_surfvy_quad_23(1, fEdge); 
  } else { 

    fUpwindQuad[23] = ser_3x3v_p1_surfvy_quad_23(-1, fSkin); 
  } 
  if (alpha[21]+alpha[20]+alpha[19]+alpha[18]-alpha[14]-alpha[13]-alpha[12]-alpha[11]-alpha[10]-alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5]+alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[24] = ser_3x3v_p1_surfvy_quad_24(1, fEdge); 
  } else { 

    fUpwindQuad[24] = ser_3x3v_p1_surfvy_quad_24(-1, fSkin); 
  } 
  if ((-alpha[21])-alpha[20]+alpha[19]-alpha[18]-alpha[14]-alpha[13]+alpha[12]-alpha[11]-alpha[10]+alpha[9]+alpha[8]-alpha[7]-alpha[6]+alpha[5]+alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[25] = ser_3x3v_p1_surfvy_quad_25(1, fEdge); 
  } else { 

    fUpwindQuad[25] = ser_3x3v_p1_surfvy_quad_25(-1, fSkin); 
  } 
  if (alpha[21]-alpha[20]-alpha[19]+alpha[18]-alpha[14]+alpha[13]-alpha[12]-alpha[11]+alpha[10]-alpha[9]-alpha[8]+alpha[7]-alpha[6]+alpha[5]+alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[26] = ser_3x3v_p1_surfvy_quad_26(1, fEdge); 
  } else { 

    fUpwindQuad[26] = ser_3x3v_p1_surfvy_quad_26(-1, fSkin); 
  } 
  if ((-alpha[21])+alpha[20]-alpha[19]-alpha[18]-alpha[14]+alpha[13]+alpha[12]-alpha[11]+alpha[10]+alpha[9]-alpha[8]-alpha[7]+alpha[6]+alpha[5]+alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[27] = ser_3x3v_p1_surfvy_quad_27(1, fEdge); 
  } else { 

    fUpwindQuad[27] = ser_3x3v_p1_surfvy_quad_27(-1, fSkin); 
  } 
  if ((-alpha[21])+alpha[20]-alpha[19]-alpha[18]+alpha[14]-alpha[13]-alpha[12]+alpha[11]-alpha[10]-alpha[9]-alpha[8]-alpha[7]+alpha[6]+alpha[5]+alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[28] = ser_3x3v_p1_surfvy_quad_28(1, fEdge); 
  } else { 

    fUpwindQuad[28] = ser_3x3v_p1_surfvy_quad_28(-1, fSkin); 
  } 
  if (alpha[21]-alpha[20]-alpha[19]+alpha[18]+alpha[14]-alpha[13]+alpha[12]+alpha[11]-alpha[10]+alpha[9]-alpha[8]+alpha[7]-alpha[6]+alpha[5]+alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[29] = ser_3x3v_p1_surfvy_quad_29(1, fEdge); 
  } else { 

    fUpwindQuad[29] = ser_3x3v_p1_surfvy_quad_29(-1, fSkin); 
  } 
  if ((-alpha[21])-alpha[20]+alpha[19]-alpha[18]+alpha[14]+alpha[13]-alpha[12]+alpha[11]+alpha[10]-alpha[9]+alpha[8]-alpha[7]-alpha[6]+alpha[5]+alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[30] = ser_3x3v_p1_surfvy_quad_30(1, fEdge); 
  } else { 

    fUpwindQuad[30] = ser_3x3v_p1_surfvy_quad_30(-1, fSkin); 
  } 
  if (alpha[21]+alpha[20]+alpha[19]+alpha[18]+alpha[14]+alpha[13]+alpha[12]+alpha[11]+alpha[10]+alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5]+alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[31] = ser_3x3v_p1_surfvy_quad_31(1, fEdge); 
  } else { 

    fUpwindQuad[31] = ser_3x3v_p1_surfvy_quad_31(-1, fSkin); 
  } 

  fUpwind[0] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]+fUpwindQuad[29]+fUpwindQuad[28]+fUpwindQuad[27]+fUpwindQuad[26]+fUpwindQuad[25]+fUpwindQuad[24]+fUpwindQuad[23]+fUpwindQuad[22]+fUpwindQuad[21]+fUpwindQuad[20]+fUpwindQuad[19]+fUpwindQuad[18]+fUpwindQuad[17]+fUpwindQuad[16]+fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[1] = 0.1767766952966368*(fUpwindQuad[31]-1.0*fUpwindQuad[30]+fUpwindQuad[29]-1.0*fUpwindQuad[28]+fUpwindQuad[27]-1.0*fUpwindQuad[26]+fUpwindQuad[25]-1.0*fUpwindQuad[24]+fUpwindQuad[23]-1.0*fUpwindQuad[22]+fUpwindQuad[21]-1.0*fUpwindQuad[20]+fUpwindQuad[19]-1.0*fUpwindQuad[18]+fUpwindQuad[17]-1.0*fUpwindQuad[16]+fUpwindQuad[15]-1.0*fUpwindQuad[14]+fUpwindQuad[13]-1.0*fUpwindQuad[12]+fUpwindQuad[11]-1.0*fUpwindQuad[10]+fUpwindQuad[9]-1.0*fUpwindQuad[8]+fUpwindQuad[7]-1.0*fUpwindQuad[6]+fUpwindQuad[5]-1.0*fUpwindQuad[4]+fUpwindQuad[3]-1.0*fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[2] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]-1.0*(fUpwindQuad[29]+fUpwindQuad[28])+fUpwindQuad[27]+fUpwindQuad[26]-1.0*(fUpwindQuad[25]+fUpwindQuad[24])+fUpwindQuad[23]+fUpwindQuad[22]-1.0*(fUpwindQuad[21]+fUpwindQuad[20])+fUpwindQuad[19]+fUpwindQuad[18]-1.0*(fUpwindQuad[17]+fUpwindQuad[16])+fUpwindQuad[15]+fUpwindQuad[14]-1.0*(fUpwindQuad[13]+fUpwindQuad[12])+fUpwindQuad[11]+fUpwindQuad[10]-1.0*(fUpwindQuad[9]+fUpwindQuad[8])+fUpwindQuad[7]+fUpwindQuad[6]-1.0*(fUpwindQuad[5]+fUpwindQuad[4])+fUpwindQuad[3]+fUpwindQuad[2]-1.0*(fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[3] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]+fUpwindQuad[29]+fUpwindQuad[28]-1.0*(fUpwindQuad[27]+fUpwindQuad[26]+fUpwindQuad[25]+fUpwindQuad[24])+fUpwindQuad[23]+fUpwindQuad[22]+fUpwindQuad[21]+fUpwindQuad[20]-1.0*(fUpwindQuad[19]+fUpwindQuad[18]+fUpwindQuad[17]+fUpwindQuad[16])+fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12]-1.0*(fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8])+fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]-1.0*(fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[4] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]+fUpwindQuad[29]+fUpwindQuad[28]+fUpwindQuad[27]+fUpwindQuad[26]+fUpwindQuad[25]+fUpwindQuad[24]-1.0*(fUpwindQuad[23]+fUpwindQuad[22]+fUpwindQuad[21]+fUpwindQuad[20]+fUpwindQuad[19]+fUpwindQuad[18]+fUpwindQuad[17]+fUpwindQuad[16])+fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8]-1.0*(fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[5] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]+fUpwindQuad[29]+fUpwindQuad[28]+fUpwindQuad[27]+fUpwindQuad[26]+fUpwindQuad[25]+fUpwindQuad[24]+fUpwindQuad[23]+fUpwindQuad[22]+fUpwindQuad[21]+fUpwindQuad[20]+fUpwindQuad[19]+fUpwindQuad[18]+fUpwindQuad[17]+fUpwindQuad[16]-1.0*(fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[6] = 0.1767766952966368*(fUpwindQuad[31]-1.0*(fUpwindQuad[30]+fUpwindQuad[29])+fUpwindQuad[28]+fUpwindQuad[27]-1.0*(fUpwindQuad[26]+fUpwindQuad[25])+fUpwindQuad[24]+fUpwindQuad[23]-1.0*(fUpwindQuad[22]+fUpwindQuad[21])+fUpwindQuad[20]+fUpwindQuad[19]-1.0*(fUpwindQuad[18]+fUpwindQuad[17])+fUpwindQuad[16]+fUpwindQuad[15]-1.0*(fUpwindQuad[14]+fUpwindQuad[13])+fUpwindQuad[12]+fUpwindQuad[11]-1.0*(fUpwindQuad[10]+fUpwindQuad[9])+fUpwindQuad[8]+fUpwindQuad[7]-1.0*(fUpwindQuad[6]+fUpwindQuad[5])+fUpwindQuad[4]+fUpwindQuad[3]-1.0*(fUpwindQuad[2]+fUpwindQuad[1])+fUpwindQuad[0]); 
  fUpwind[7] = 0.1767766952966368*(fUpwindQuad[31]-1.0*fUpwindQuad[30]+fUpwindQuad[29]-1.0*(fUpwindQuad[28]+fUpwindQuad[27])+fUpwindQuad[26]-1.0*fUpwindQuad[25]+fUpwindQuad[24]+fUpwindQuad[23]-1.0*fUpwindQuad[22]+fUpwindQuad[21]-1.0*(fUpwindQuad[20]+fUpwindQuad[19])+fUpwindQuad[18]-1.0*fUpwindQuad[17]+fUpwindQuad[16]+fUpwindQuad[15]-1.0*fUpwindQuad[14]+fUpwindQuad[13]-1.0*(fUpwindQuad[12]+fUpwindQuad[11])+fUpwindQuad[10]-1.0*fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]-1.0*fUpwindQuad[6]+fUpwindQuad[5]-1.0*(fUpwindQuad[4]+fUpwindQuad[3])+fUpwindQuad[2]-1.0*fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[8] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]-1.0*(fUpwindQuad[29]+fUpwindQuad[28]+fUpwindQuad[27]+fUpwindQuad[26])+fUpwindQuad[25]+fUpwindQuad[24]+fUpwindQuad[23]+fUpwindQuad[22]-1.0*(fUpwindQuad[21]+fUpwindQuad[20]+fUpwindQuad[19]+fUpwindQuad[18])+fUpwindQuad[17]+fUpwindQuad[16]+fUpwindQuad[15]+fUpwindQuad[14]-1.0*(fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10])+fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6]-1.0*(fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2])+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[9] = 0.1767766952966368*(fUpwindQuad[31]-1.0*fUpwindQuad[30]+fUpwindQuad[29]-1.0*fUpwindQuad[28]+fUpwindQuad[27]-1.0*fUpwindQuad[26]+fUpwindQuad[25]-1.0*(fUpwindQuad[24]+fUpwindQuad[23])+fUpwindQuad[22]-1.0*fUpwindQuad[21]+fUpwindQuad[20]-1.0*fUpwindQuad[19]+fUpwindQuad[18]-1.0*fUpwindQuad[17]+fUpwindQuad[16]+fUpwindQuad[15]-1.0*fUpwindQuad[14]+fUpwindQuad[13]-1.0*fUpwindQuad[12]+fUpwindQuad[11]-1.0*fUpwindQuad[10]+fUpwindQuad[9]-1.0*(fUpwindQuad[8]+fUpwindQuad[7])+fUpwindQuad[6]-1.0*fUpwindQuad[5]+fUpwindQuad[4]-1.0*fUpwindQuad[3]+fUpwindQuad[2]-1.0*fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[10] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]-1.0*(fUpwindQuad[29]+fUpwindQuad[28])+fUpwindQuad[27]+fUpwindQuad[26]-1.0*(fUpwindQuad[25]+fUpwindQuad[24]+fUpwindQuad[23]+fUpwindQuad[22])+fUpwindQuad[21]+fUpwindQuad[20]-1.0*(fUpwindQuad[19]+fUpwindQuad[18])+fUpwindQuad[17]+fUpwindQuad[16]+fUpwindQuad[15]+fUpwindQuad[14]-1.0*(fUpwindQuad[13]+fUpwindQuad[12])+fUpwindQuad[11]+fUpwindQuad[10]-1.0*(fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6])+fUpwindQuad[5]+fUpwindQuad[4]-1.0*(fUpwindQuad[3]+fUpwindQuad[2])+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[11] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]+fUpwindQuad[29]+fUpwindQuad[28]-1.0*(fUpwindQuad[27]+fUpwindQuad[26]+fUpwindQuad[25]+fUpwindQuad[24]+fUpwindQuad[23]+fUpwindQuad[22]+fUpwindQuad[21]+fUpwindQuad[20])+fUpwindQuad[19]+fUpwindQuad[18]+fUpwindQuad[17]+fUpwindQuad[16]+fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12]-1.0*(fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4])+fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[12] = 0.1767766952966368*(fUpwindQuad[31]-1.0*fUpwindQuad[30]+fUpwindQuad[29]-1.0*fUpwindQuad[28]+fUpwindQuad[27]-1.0*fUpwindQuad[26]+fUpwindQuad[25]-1.0*fUpwindQuad[24]+fUpwindQuad[23]-1.0*fUpwindQuad[22]+fUpwindQuad[21]-1.0*fUpwindQuad[20]+fUpwindQuad[19]-1.0*fUpwindQuad[18]+fUpwindQuad[17]-1.0*(fUpwindQuad[16]+fUpwindQuad[15])+fUpwindQuad[14]-1.0*fUpwindQuad[13]+fUpwindQuad[12]-1.0*fUpwindQuad[11]+fUpwindQuad[10]-1.0*fUpwindQuad[9]+fUpwindQuad[8]-1.0*fUpwindQuad[7]+fUpwindQuad[6]-1.0*fUpwindQuad[5]+fUpwindQuad[4]-1.0*fUpwindQuad[3]+fUpwindQuad[2]-1.0*fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[13] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]-1.0*(fUpwindQuad[29]+fUpwindQuad[28])+fUpwindQuad[27]+fUpwindQuad[26]-1.0*(fUpwindQuad[25]+fUpwindQuad[24])+fUpwindQuad[23]+fUpwindQuad[22]-1.0*(fUpwindQuad[21]+fUpwindQuad[20])+fUpwindQuad[19]+fUpwindQuad[18]-1.0*(fUpwindQuad[17]+fUpwindQuad[16]+fUpwindQuad[15]+fUpwindQuad[14])+fUpwindQuad[13]+fUpwindQuad[12]-1.0*(fUpwindQuad[11]+fUpwindQuad[10])+fUpwindQuad[9]+fUpwindQuad[8]-1.0*(fUpwindQuad[7]+fUpwindQuad[6])+fUpwindQuad[5]+fUpwindQuad[4]-1.0*(fUpwindQuad[3]+fUpwindQuad[2])+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[14] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]+fUpwindQuad[29]+fUpwindQuad[28]-1.0*(fUpwindQuad[27]+fUpwindQuad[26]+fUpwindQuad[25]+fUpwindQuad[24])+fUpwindQuad[23]+fUpwindQuad[22]+fUpwindQuad[21]+fUpwindQuad[20]-1.0*(fUpwindQuad[19]+fUpwindQuad[18]+fUpwindQuad[17]+fUpwindQuad[16]+fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12])+fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8]-1.0*(fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4])+fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[15] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]+fUpwindQuad[29]+fUpwindQuad[28]+fUpwindQuad[27]+fUpwindQuad[26]+fUpwindQuad[25]+fUpwindQuad[24]-1.0*(fUpwindQuad[23]+fUpwindQuad[22]+fUpwindQuad[21]+fUpwindQuad[20]+fUpwindQuad[19]+fUpwindQuad[18]+fUpwindQuad[17]+fUpwindQuad[16]+fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8])+fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[16] = 0.1767766952966368*(fUpwindQuad[31]-1.0*(fUpwindQuad[30]+fUpwindQuad[29])+fUpwindQuad[28]-1.0*fUpwindQuad[27]+fUpwindQuad[26]+fUpwindQuad[25]-1.0*fUpwindQuad[24]+fUpwindQuad[23]-1.0*(fUpwindQuad[22]+fUpwindQuad[21])+fUpwindQuad[20]-1.0*fUpwindQuad[19]+fUpwindQuad[18]+fUpwindQuad[17]-1.0*fUpwindQuad[16]+fUpwindQuad[15]-1.0*(fUpwindQuad[14]+fUpwindQuad[13])+fUpwindQuad[12]-1.0*fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]-1.0*fUpwindQuad[8]+fUpwindQuad[7]-1.0*(fUpwindQuad[6]+fUpwindQuad[5])+fUpwindQuad[4]-1.0*fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[17] = 0.1767766952966368*(fUpwindQuad[31]-1.0*(fUpwindQuad[30]+fUpwindQuad[29])+fUpwindQuad[28]+fUpwindQuad[27]-1.0*(fUpwindQuad[26]+fUpwindQuad[25])+fUpwindQuad[24]-1.0*fUpwindQuad[23]+fUpwindQuad[22]+fUpwindQuad[21]-1.0*(fUpwindQuad[20]+fUpwindQuad[19])+fUpwindQuad[18]+fUpwindQuad[17]-1.0*fUpwindQuad[16]+fUpwindQuad[15]-1.0*(fUpwindQuad[14]+fUpwindQuad[13])+fUpwindQuad[12]+fUpwindQuad[11]-1.0*(fUpwindQuad[10]+fUpwindQuad[9])+fUpwindQuad[8]-1.0*fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]-1.0*(fUpwindQuad[4]+fUpwindQuad[3])+fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[18] = 0.1767766952966368*(fUpwindQuad[31]-1.0*fUpwindQuad[30]+fUpwindQuad[29]-1.0*(fUpwindQuad[28]+fUpwindQuad[27])+fUpwindQuad[26]-1.0*fUpwindQuad[25]+fUpwindQuad[24]-1.0*fUpwindQuad[23]+fUpwindQuad[22]-1.0*fUpwindQuad[21]+fUpwindQuad[20]+fUpwindQuad[19]-1.0*fUpwindQuad[18]+fUpwindQuad[17]-1.0*fUpwindQuad[16]+fUpwindQuad[15]-1.0*fUpwindQuad[14]+fUpwindQuad[13]-1.0*(fUpwindQuad[12]+fUpwindQuad[11])+fUpwindQuad[10]-1.0*fUpwindQuad[9]+fUpwindQuad[8]-1.0*fUpwindQuad[7]+fUpwindQuad[6]-1.0*fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]-1.0*fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[19] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]-1.0*(fUpwindQuad[29]+fUpwindQuad[28]+fUpwindQuad[27]+fUpwindQuad[26])+fUpwindQuad[25]+fUpwindQuad[24]-1.0*(fUpwindQuad[23]+fUpwindQuad[22])+fUpwindQuad[21]+fUpwindQuad[20]+fUpwindQuad[19]+fUpwindQuad[18]-1.0*(fUpwindQuad[17]+fUpwindQuad[16])+fUpwindQuad[15]+fUpwindQuad[14]-1.0*(fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10])+fUpwindQuad[9]+fUpwindQuad[8]-1.0*(fUpwindQuad[7]+fUpwindQuad[6])+fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2]-1.0*(fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[20] = 0.1767766952966368*(fUpwindQuad[31]-1.0*(fUpwindQuad[30]+fUpwindQuad[29])+fUpwindQuad[28]+fUpwindQuad[27]-1.0*(fUpwindQuad[26]+fUpwindQuad[25])+fUpwindQuad[24]+fUpwindQuad[23]-1.0*(fUpwindQuad[22]+fUpwindQuad[21])+fUpwindQuad[20]+fUpwindQuad[19]-1.0*(fUpwindQuad[18]+fUpwindQuad[17])+fUpwindQuad[16]-1.0*fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]-1.0*(fUpwindQuad[12]+fUpwindQuad[11])+fUpwindQuad[10]+fUpwindQuad[9]-1.0*(fUpwindQuad[8]+fUpwindQuad[7])+fUpwindQuad[6]+fUpwindQuad[5]-1.0*(fUpwindQuad[4]+fUpwindQuad[3])+fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[21] = 0.1767766952966368*(fUpwindQuad[31]-1.0*fUpwindQuad[30]+fUpwindQuad[29]-1.0*(fUpwindQuad[28]+fUpwindQuad[27])+fUpwindQuad[26]-1.0*fUpwindQuad[25]+fUpwindQuad[24]+fUpwindQuad[23]-1.0*fUpwindQuad[22]+fUpwindQuad[21]-1.0*(fUpwindQuad[20]+fUpwindQuad[19])+fUpwindQuad[18]-1.0*fUpwindQuad[17]+fUpwindQuad[16]-1.0*fUpwindQuad[15]+fUpwindQuad[14]-1.0*fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]-1.0*fUpwindQuad[10]+fUpwindQuad[9]-1.0*(fUpwindQuad[8]+fUpwindQuad[7])+fUpwindQuad[6]-1.0*fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]-1.0*fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[22] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]-1.0*(fUpwindQuad[29]+fUpwindQuad[28]+fUpwindQuad[27]+fUpwindQuad[26])+fUpwindQuad[25]+fUpwindQuad[24]+fUpwindQuad[23]+fUpwindQuad[22]-1.0*(fUpwindQuad[21]+fUpwindQuad[20]+fUpwindQuad[19]+fUpwindQuad[18])+fUpwindQuad[17]+fUpwindQuad[16]-1.0*(fUpwindQuad[15]+fUpwindQuad[14])+fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10]-1.0*(fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6])+fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2]-1.0*(fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[23] = 0.1767766952966368*(fUpwindQuad[31]-1.0*fUpwindQuad[30]+fUpwindQuad[29]-1.0*fUpwindQuad[28]+fUpwindQuad[27]-1.0*fUpwindQuad[26]+fUpwindQuad[25]-1.0*(fUpwindQuad[24]+fUpwindQuad[23])+fUpwindQuad[22]-1.0*fUpwindQuad[21]+fUpwindQuad[20]-1.0*fUpwindQuad[19]+fUpwindQuad[18]-1.0*fUpwindQuad[17]+fUpwindQuad[16]-1.0*fUpwindQuad[15]+fUpwindQuad[14]-1.0*fUpwindQuad[13]+fUpwindQuad[12]-1.0*fUpwindQuad[11]+fUpwindQuad[10]-1.0*fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]-1.0*fUpwindQuad[6]+fUpwindQuad[5]-1.0*fUpwindQuad[4]+fUpwindQuad[3]-1.0*fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[24] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]-1.0*(fUpwindQuad[29]+fUpwindQuad[28])+fUpwindQuad[27]+fUpwindQuad[26]-1.0*(fUpwindQuad[25]+fUpwindQuad[24]+fUpwindQuad[23]+fUpwindQuad[22])+fUpwindQuad[21]+fUpwindQuad[20]-1.0*(fUpwindQuad[19]+fUpwindQuad[18])+fUpwindQuad[17]+fUpwindQuad[16]-1.0*(fUpwindQuad[15]+fUpwindQuad[14])+fUpwindQuad[13]+fUpwindQuad[12]-1.0*(fUpwindQuad[11]+fUpwindQuad[10])+fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6]-1.0*(fUpwindQuad[5]+fUpwindQuad[4])+fUpwindQuad[3]+fUpwindQuad[2]-1.0*(fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[25] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]+fUpwindQuad[29]+fUpwindQuad[28]-1.0*(fUpwindQuad[27]+fUpwindQuad[26]+fUpwindQuad[25]+fUpwindQuad[24]+fUpwindQuad[23]+fUpwindQuad[22]+fUpwindQuad[21]+fUpwindQuad[20])+fUpwindQuad[19]+fUpwindQuad[18]+fUpwindQuad[17]+fUpwindQuad[16]-1.0*(fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12])+fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]-1.0*(fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[26] = 0.1767766952966368*(fUpwindQuad[31]-1.0*(fUpwindQuad[30]+fUpwindQuad[29])+fUpwindQuad[28]-1.0*fUpwindQuad[27]+fUpwindQuad[26]+fUpwindQuad[25]-1.0*(fUpwindQuad[24]+fUpwindQuad[23])+fUpwindQuad[22]+fUpwindQuad[21]-1.0*fUpwindQuad[20]+fUpwindQuad[19]-1.0*(fUpwindQuad[18]+fUpwindQuad[17])+fUpwindQuad[16]+fUpwindQuad[15]-1.0*(fUpwindQuad[14]+fUpwindQuad[13])+fUpwindQuad[12]-1.0*fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]-1.0*(fUpwindQuad[8]+fUpwindQuad[7])+fUpwindQuad[6]+fUpwindQuad[5]-1.0*fUpwindQuad[4]+fUpwindQuad[3]-1.0*(fUpwindQuad[2]+fUpwindQuad[1])+fUpwindQuad[0]); 
  fUpwind[27] = 0.1767766952966368*(fUpwindQuad[31]-1.0*(fUpwindQuad[30]+fUpwindQuad[29])+fUpwindQuad[28]-1.0*fUpwindQuad[27]+fUpwindQuad[26]+fUpwindQuad[25]-1.0*fUpwindQuad[24]+fUpwindQuad[23]-1.0*(fUpwindQuad[22]+fUpwindQuad[21])+fUpwindQuad[20]-1.0*fUpwindQuad[19]+fUpwindQuad[18]+fUpwindQuad[17]-1.0*(fUpwindQuad[16]+fUpwindQuad[15])+fUpwindQuad[14]+fUpwindQuad[13]-1.0*fUpwindQuad[12]+fUpwindQuad[11]-1.0*(fUpwindQuad[10]+fUpwindQuad[9])+fUpwindQuad[8]-1.0*fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]-1.0*fUpwindQuad[4]+fUpwindQuad[3]-1.0*(fUpwindQuad[2]+fUpwindQuad[1])+fUpwindQuad[0]); 
  fUpwind[28] = 0.1767766952966368*(fUpwindQuad[31]-1.0*(fUpwindQuad[30]+fUpwindQuad[29])+fUpwindQuad[28]+fUpwindQuad[27]-1.0*(fUpwindQuad[26]+fUpwindQuad[25])+fUpwindQuad[24]-1.0*fUpwindQuad[23]+fUpwindQuad[22]+fUpwindQuad[21]-1.0*(fUpwindQuad[20]+fUpwindQuad[19])+fUpwindQuad[18]+fUpwindQuad[17]-1.0*(fUpwindQuad[16]+fUpwindQuad[15])+fUpwindQuad[14]+fUpwindQuad[13]-1.0*(fUpwindQuad[12]+fUpwindQuad[11])+fUpwindQuad[10]+fUpwindQuad[9]-1.0*fUpwindQuad[8]+fUpwindQuad[7]-1.0*(fUpwindQuad[6]+fUpwindQuad[5])+fUpwindQuad[4]+fUpwindQuad[3]-1.0*(fUpwindQuad[2]+fUpwindQuad[1])+fUpwindQuad[0]); 
  fUpwind[29] = 0.1767766952966368*(fUpwindQuad[31]-1.0*fUpwindQuad[30]+fUpwindQuad[29]-1.0*(fUpwindQuad[28]+fUpwindQuad[27])+fUpwindQuad[26]-1.0*fUpwindQuad[25]+fUpwindQuad[24]-1.0*fUpwindQuad[23]+fUpwindQuad[22]-1.0*fUpwindQuad[21]+fUpwindQuad[20]+fUpwindQuad[19]-1.0*fUpwindQuad[18]+fUpwindQuad[17]-1.0*(fUpwindQuad[16]+fUpwindQuad[15])+fUpwindQuad[14]-1.0*fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]-1.0*fUpwindQuad[10]+fUpwindQuad[9]-1.0*fUpwindQuad[8]+fUpwindQuad[7]-1.0*fUpwindQuad[6]+fUpwindQuad[5]-1.0*(fUpwindQuad[4]+fUpwindQuad[3])+fUpwindQuad[2]-1.0*fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[30] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]-1.0*(fUpwindQuad[29]+fUpwindQuad[28]+fUpwindQuad[27]+fUpwindQuad[26])+fUpwindQuad[25]+fUpwindQuad[24]-1.0*(fUpwindQuad[23]+fUpwindQuad[22])+fUpwindQuad[21]+fUpwindQuad[20]+fUpwindQuad[19]+fUpwindQuad[18]-1.0*(fUpwindQuad[17]+fUpwindQuad[16]+fUpwindQuad[15]+fUpwindQuad[14])+fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10]-1.0*(fUpwindQuad[9]+fUpwindQuad[8])+fUpwindQuad[7]+fUpwindQuad[6]-1.0*(fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2])+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[31] = 0.1767766952966368*(fUpwindQuad[31]-1.0*(fUpwindQuad[30]+fUpwindQuad[29])+fUpwindQuad[28]-1.0*fUpwindQuad[27]+fUpwindQuad[26]+fUpwindQuad[25]-1.0*(fUpwindQuad[24]+fUpwindQuad[23])+fUpwindQuad[22]+fUpwindQuad[21]-1.0*fUpwindQuad[20]+fUpwindQuad[19]-1.0*(fUpwindQuad[18]+fUpwindQuad[17])+fUpwindQuad[16]-1.0*fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]-1.0*fUpwindQuad[12]+fUpwindQuad[11]-1.0*(fUpwindQuad[10]+fUpwindQuad[9])+fUpwindQuad[8]+fUpwindQuad[7]-1.0*(fUpwindQuad[6]+fUpwindQuad[5])+fUpwindQuad[4]-1.0*fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 

  Ghat[0] += 0.1767766952966368*(alpha[21]*fUpwind[21]+alpha[20]*fUpwind[20]+alpha[19]*fUpwind[19]+alpha[18]*fUpwind[18]+alpha[14]*fUpwind[14]+alpha[13]*fUpwind[13]+alpha[12]*fUpwind[12]+alpha[11]*fUpwind[11]+alpha[10]*fUpwind[10]+alpha[9]*fUpwind[9]+alpha[8]*fUpwind[8]+alpha[7]*fUpwind[7]+alpha[6]*fUpwind[6]+alpha[5]*fUpwind[5]+alpha[4]*fUpwind[4]+alpha[3]*fUpwind[3]+alpha[2]*fUpwind[2]+alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] += 0.1767766952966368*(alpha[19]*fUpwind[26]+alpha[14]*fUpwind[21]+fUpwind[14]*alpha[21]+alpha[13]*fUpwind[20]+fUpwind[13]*alpha[20]+alpha[11]*fUpwind[18]+fUpwind[11]*alpha[18]+alpha[10]*fUpwind[17]+alpha[8]*fUpwind[16]+alpha[5]*fUpwind[12]+fUpwind[5]*alpha[12]+alpha[4]*fUpwind[9]+fUpwind[4]*alpha[9]+alpha[3]*fUpwind[7]+fUpwind[3]*alpha[7]+alpha[2]*fUpwind[6]+fUpwind[2]*alpha[6]+alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]); 
  Ghat[2] += 0.1767766952966368*(alpha[21]*fUpwind[27]+alpha[18]*fUpwind[26]+alpha[14]*fUpwind[22]+alpha[12]*fUpwind[20]+fUpwind[12]*alpha[20]+alpha[11]*fUpwind[19]+fUpwind[11]*alpha[19]+alpha[9]*fUpwind[17]+alpha[7]*fUpwind[16]+alpha[5]*fUpwind[13]+fUpwind[5]*alpha[13]+alpha[4]*fUpwind[10]+fUpwind[4]*alpha[10]+alpha[3]*fUpwind[8]+fUpwind[3]*alpha[8]+alpha[1]*fUpwind[6]+fUpwind[1]*alpha[6]+alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2]); 
  Ghat[3] += 0.1767766952966368*(alpha[20]*fUpwind[27]+alpha[13]*fUpwind[22]+alpha[12]*fUpwind[21]+fUpwind[12]*alpha[21]+alpha[10]*fUpwind[19]+fUpwind[10]*alpha[19]+alpha[9]*fUpwind[18]+fUpwind[9]*alpha[18]+alpha[6]*fUpwind[16]+alpha[5]*fUpwind[14]+fUpwind[5]*alpha[14]+alpha[4]*fUpwind[11]+fUpwind[4]*alpha[11]+alpha[2]*fUpwind[8]+fUpwind[2]*alpha[8]+alpha[1]*fUpwind[7]+fUpwind[1]*alpha[7]+alpha[0]*fUpwind[3]+fUpwind[0]*alpha[3]); 
  Ghat[4] += 0.1767766952966368*(alpha[21]*fUpwind[29]+alpha[20]*fUpwind[28]+alpha[14]*fUpwind[25]+alpha[13]*fUpwind[24]+alpha[12]*fUpwind[23]+alpha[8]*fUpwind[19]+fUpwind[8]*alpha[19]+alpha[7]*fUpwind[18]+fUpwind[7]*alpha[18]+alpha[6]*fUpwind[17]+alpha[5]*fUpwind[15]+alpha[3]*fUpwind[11]+fUpwind[3]*alpha[11]+alpha[2]*fUpwind[10]+fUpwind[2]*alpha[10]+alpha[1]*fUpwind[9]+fUpwind[1]*alpha[9]+alpha[0]*fUpwind[4]+fUpwind[0]*alpha[4]); 
  Ghat[5] += 0.1767766952966368*(alpha[19]*fUpwind[30]+alpha[18]*fUpwind[29]+alpha[11]*fUpwind[25]+alpha[10]*fUpwind[24]+alpha[9]*fUpwind[23]+alpha[8]*fUpwind[22]+alpha[7]*fUpwind[21]+fUpwind[7]*alpha[21]+alpha[6]*fUpwind[20]+fUpwind[6]*alpha[20]+alpha[4]*fUpwind[15]+alpha[3]*fUpwind[14]+fUpwind[3]*alpha[14]+alpha[2]*fUpwind[13]+fUpwind[2]*alpha[13]+alpha[1]*fUpwind[12]+fUpwind[1]*alpha[12]+alpha[0]*fUpwind[5]+fUpwind[0]*alpha[5]); 
  Ghat[6] += 0.1767766952966368*(alpha[14]*fUpwind[27]+alpha[11]*fUpwind[26]+alpha[21]*fUpwind[22]+alpha[5]*fUpwind[20]+fUpwind[5]*alpha[20]+alpha[18]*fUpwind[19]+fUpwind[18]*alpha[19]+alpha[4]*fUpwind[17]+alpha[3]*fUpwind[16]+alpha[12]*fUpwind[13]+fUpwind[12]*alpha[13]+alpha[9]*fUpwind[10]+fUpwind[9]*alpha[10]+alpha[7]*fUpwind[8]+fUpwind[7]*alpha[8]+alpha[0]*fUpwind[6]+fUpwind[0]*alpha[6]+alpha[1]*fUpwind[2]+fUpwind[1]*alpha[2]); 
  Ghat[7] += 0.1767766952966368*(alpha[13]*fUpwind[27]+alpha[10]*fUpwind[26]+alpha[20]*fUpwind[22]+alpha[5]*fUpwind[21]+fUpwind[5]*alpha[21]+fUpwind[17]*alpha[19]+alpha[4]*fUpwind[18]+fUpwind[4]*alpha[18]+alpha[2]*fUpwind[16]+alpha[12]*fUpwind[14]+fUpwind[12]*alpha[14]+alpha[9]*fUpwind[11]+fUpwind[9]*alpha[11]+alpha[6]*fUpwind[8]+fUpwind[6]*alpha[8]+alpha[0]*fUpwind[7]+fUpwind[0]*alpha[7]+alpha[1]*fUpwind[3]+fUpwind[1]*alpha[3]); 
  Ghat[8] += 0.1767766952966368*(alpha[12]*fUpwind[27]+alpha[9]*fUpwind[26]+alpha[5]*fUpwind[22]+alpha[20]*fUpwind[21]+fUpwind[20]*alpha[21]+alpha[4]*fUpwind[19]+fUpwind[4]*alpha[19]+fUpwind[17]*alpha[18]+alpha[1]*fUpwind[16]+alpha[13]*fUpwind[14]+fUpwind[13]*alpha[14]+alpha[10]*fUpwind[11]+fUpwind[10]*alpha[11]+alpha[0]*fUpwind[8]+fUpwind[0]*alpha[8]+alpha[6]*fUpwind[7]+fUpwind[6]*alpha[7]+alpha[2]*fUpwind[3]+fUpwind[2]*alpha[3]); 
  Ghat[9] += 0.1767766952966368*(alpha[14]*fUpwind[29]+alpha[13]*fUpwind[28]+alpha[8]*fUpwind[26]+alpha[21]*fUpwind[25]+alpha[20]*fUpwind[24]+alpha[5]*fUpwind[23]+fUpwind[16]*alpha[19]+alpha[3]*fUpwind[18]+fUpwind[3]*alpha[18]+alpha[2]*fUpwind[17]+alpha[12]*fUpwind[15]+alpha[7]*fUpwind[11]+fUpwind[7]*alpha[11]+alpha[6]*fUpwind[10]+fUpwind[6]*alpha[10]+alpha[0]*fUpwind[9]+fUpwind[0]*alpha[9]+alpha[1]*fUpwind[4]+fUpwind[1]*alpha[4]); 
  Ghat[10] += 0.1767766952966368*(alpha[21]*fUpwind[31]+alpha[14]*fUpwind[30]+alpha[12]*fUpwind[28]+alpha[7]*fUpwind[26]+alpha[5]*fUpwind[24]+alpha[20]*fUpwind[23]+alpha[3]*fUpwind[19]+fUpwind[3]*alpha[19]+fUpwind[16]*alpha[18]+alpha[1]*fUpwind[17]+alpha[13]*fUpwind[15]+alpha[8]*fUpwind[11]+fUpwind[8]*alpha[11]+alpha[0]*fUpwind[10]+fUpwind[0]*alpha[10]+alpha[6]*fUpwind[9]+fUpwind[6]*alpha[9]+alpha[2]*fUpwind[4]+fUpwind[2]*alpha[4]); 
  Ghat[11] += 0.1767766952966368*(alpha[20]*fUpwind[31]+alpha[13]*fUpwind[30]+alpha[12]*fUpwind[29]+alpha[6]*fUpwind[26]+alpha[5]*fUpwind[25]+alpha[21]*fUpwind[23]+alpha[2]*fUpwind[19]+fUpwind[2]*alpha[19]+alpha[1]*fUpwind[18]+fUpwind[1]*alpha[18]+alpha[14]*fUpwind[15]+alpha[0]*fUpwind[11]+fUpwind[0]*alpha[11]+alpha[8]*fUpwind[10]+fUpwind[8]*alpha[10]+alpha[7]*fUpwind[9]+fUpwind[7]*alpha[9]+alpha[3]*fUpwind[4]+fUpwind[3]*alpha[4]); 
  Ghat[12] += 0.1767766952966368*(alpha[19]*fUpwind[31]+alpha[11]*fUpwind[29]+alpha[10]*fUpwind[28]+alpha[8]*fUpwind[27]+alpha[18]*fUpwind[25]+alpha[4]*fUpwind[23]+alpha[3]*fUpwind[21]+fUpwind[3]*alpha[21]+alpha[2]*fUpwind[20]+fUpwind[2]*alpha[20]+alpha[9]*fUpwind[15]+alpha[7]*fUpwind[14]+fUpwind[7]*alpha[14]+alpha[6]*fUpwind[13]+fUpwind[6]*alpha[13]+alpha[0]*fUpwind[12]+fUpwind[0]*alpha[12]+alpha[1]*fUpwind[5]+fUpwind[1]*alpha[5]); 
  Ghat[13] += 0.1767766952966368*(alpha[18]*fUpwind[31]+alpha[11]*fUpwind[30]+alpha[9]*fUpwind[28]+alpha[7]*fUpwind[27]+alpha[19]*fUpwind[25]+alpha[4]*fUpwind[24]+alpha[3]*fUpwind[22]+fUpwind[16]*alpha[21]+alpha[1]*fUpwind[20]+fUpwind[1]*alpha[20]+alpha[10]*fUpwind[15]+alpha[8]*fUpwind[14]+fUpwind[8]*alpha[14]+alpha[0]*fUpwind[13]+fUpwind[0]*alpha[13]+alpha[6]*fUpwind[12]+fUpwind[6]*alpha[12]+alpha[2]*fUpwind[5]+fUpwind[2]*alpha[5]); 
  Ghat[14] += 0.1767766952966368*(alpha[10]*fUpwind[30]+alpha[9]*fUpwind[29]+alpha[6]*fUpwind[27]+alpha[4]*fUpwind[25]+alpha[19]*fUpwind[24]+alpha[18]*fUpwind[23]+alpha[2]*fUpwind[22]+alpha[1]*fUpwind[21]+fUpwind[1]*alpha[21]+fUpwind[16]*alpha[20]+alpha[11]*fUpwind[15]+alpha[0]*fUpwind[14]+fUpwind[0]*alpha[14]+alpha[8]*fUpwind[13]+fUpwind[8]*alpha[13]+alpha[7]*fUpwind[12]+fUpwind[7]*alpha[12]+alpha[3]*fUpwind[5]+fUpwind[3]*alpha[5]); 
  Ghat[15] += 0.1767766952966368*(alpha[8]*fUpwind[30]+alpha[7]*fUpwind[29]+alpha[6]*fUpwind[28]+alpha[3]*fUpwind[25]+alpha[2]*fUpwind[24]+alpha[1]*fUpwind[23]+alpha[19]*fUpwind[22]+alpha[18]*fUpwind[21]+fUpwind[18]*alpha[21]+fUpwind[17]*alpha[20]+alpha[0]*fUpwind[15]+alpha[11]*fUpwind[14]+fUpwind[11]*alpha[14]+alpha[10]*fUpwind[13]+fUpwind[10]*alpha[13]+alpha[9]*fUpwind[12]+fUpwind[9]*alpha[12]+alpha[4]*fUpwind[5]+fUpwind[4]*alpha[5]); 
  Ghat[16] += 0.1767766952966368*(alpha[5]*fUpwind[27]+alpha[4]*fUpwind[26]+alpha[12]*fUpwind[22]+alpha[13]*fUpwind[21]+fUpwind[13]*alpha[21]+alpha[14]*fUpwind[20]+fUpwind[14]*alpha[20]+alpha[9]*fUpwind[19]+fUpwind[9]*alpha[19]+alpha[10]*fUpwind[18]+fUpwind[10]*alpha[18]+alpha[11]*fUpwind[17]+alpha[0]*fUpwind[16]+alpha[1]*fUpwind[8]+fUpwind[1]*alpha[8]+alpha[2]*fUpwind[7]+fUpwind[2]*alpha[7]+alpha[3]*fUpwind[6]+fUpwind[3]*alpha[6]); 
  Ghat[17] += 0.1767766952966368*(alpha[14]*fUpwind[31]+alpha[21]*fUpwind[30]+alpha[5]*fUpwind[28]+alpha[3]*fUpwind[26]+alpha[12]*fUpwind[24]+alpha[13]*fUpwind[23]+fUpwind[15]*alpha[20]+alpha[7]*fUpwind[19]+fUpwind[7]*alpha[19]+alpha[8]*fUpwind[18]+fUpwind[8]*alpha[18]+alpha[0]*fUpwind[17]+alpha[11]*fUpwind[16]+alpha[1]*fUpwind[10]+fUpwind[1]*alpha[10]+alpha[2]*fUpwind[9]+fUpwind[2]*alpha[9]+alpha[4]*fUpwind[6]+fUpwind[4]*alpha[6]); 
  Ghat[18] += 0.1767766952966368*(alpha[13]*fUpwind[31]+alpha[20]*fUpwind[30]+alpha[5]*fUpwind[29]+alpha[2]*fUpwind[26]+alpha[12]*fUpwind[25]+alpha[14]*fUpwind[23]+fUpwind[15]*alpha[21]+alpha[6]*fUpwind[19]+fUpwind[6]*alpha[19]+alpha[0]*fUpwind[18]+fUpwind[0]*alpha[18]+alpha[8]*fUpwind[17]+alpha[10]*fUpwind[16]+alpha[1]*fUpwind[11]+fUpwind[1]*alpha[11]+alpha[3]*fUpwind[9]+fUpwind[3]*alpha[9]+alpha[4]*fUpwind[7]+fUpwind[4]*alpha[7]); 
  Ghat[19] += 0.1767766952966368*(alpha[12]*fUpwind[31]+alpha[5]*fUpwind[30]+alpha[20]*fUpwind[29]+alpha[21]*fUpwind[28]+alpha[1]*fUpwind[26]+alpha[13]*fUpwind[25]+alpha[14]*fUpwind[24]+alpha[0]*fUpwind[19]+fUpwind[0]*alpha[19]+alpha[6]*fUpwind[18]+fUpwind[6]*alpha[18]+alpha[7]*fUpwind[17]+alpha[9]*fUpwind[16]+alpha[2]*fUpwind[11]+fUpwind[2]*alpha[11]+alpha[3]*fUpwind[10]+fUpwind[3]*alpha[10]+alpha[4]*fUpwind[8]+fUpwind[4]*alpha[8]); 
  Ghat[20] += 0.1767766952966368*(alpha[11]*fUpwind[31]+alpha[18]*fUpwind[30]+alpha[19]*fUpwind[29]+alpha[4]*fUpwind[28]+alpha[3]*fUpwind[27]+alpha[9]*fUpwind[24]+alpha[10]*fUpwind[23]+alpha[7]*fUpwind[22]+alpha[8]*fUpwind[21]+fUpwind[8]*alpha[21]+alpha[0]*fUpwind[20]+fUpwind[0]*alpha[20]+alpha[14]*fUpwind[16]+alpha[1]*fUpwind[13]+fUpwind[1]*alpha[13]+alpha[2]*fUpwind[12]+fUpwind[2]*alpha[12]+alpha[5]*fUpwind[6]+fUpwind[5]*alpha[6]); 
  Ghat[21] += 0.1767766952966368*(alpha[10]*fUpwind[31]+alpha[4]*fUpwind[29]+alpha[19]*fUpwind[28]+alpha[2]*fUpwind[27]+alpha[9]*fUpwind[25]+alpha[11]*fUpwind[23]+alpha[6]*fUpwind[22]+alpha[0]*fUpwind[21]+fUpwind[0]*alpha[21]+alpha[8]*fUpwind[20]+fUpwind[8]*alpha[20]+fUpwind[15]*alpha[18]+alpha[13]*fUpwind[16]+alpha[1]*fUpwind[14]+fUpwind[1]*alpha[14]+alpha[3]*fUpwind[12]+fUpwind[3]*alpha[12]+alpha[5]*fUpwind[7]+fUpwind[5]*alpha[7]); 
  Ghat[22] += 0.1767766952966368*(alpha[9]*fUpwind[31]+alpha[4]*fUpwind[30]+alpha[18]*fUpwind[28]+alpha[1]*fUpwind[27]+alpha[10]*fUpwind[25]+alpha[11]*fUpwind[24]+alpha[0]*fUpwind[22]+alpha[6]*fUpwind[21]+fUpwind[6]*alpha[21]+alpha[7]*fUpwind[20]+fUpwind[7]*alpha[20]+fUpwind[15]*alpha[19]+alpha[12]*fUpwind[16]+alpha[2]*fUpwind[14]+fUpwind[2]*alpha[14]+alpha[3]*fUpwind[13]+fUpwind[3]*alpha[13]+alpha[5]*fUpwind[8]+fUpwind[5]*alpha[8]); 
  Ghat[23] += 0.1767766952966368*(alpha[8]*fUpwind[31]+alpha[3]*fUpwind[29]+alpha[2]*fUpwind[28]+alpha[19]*fUpwind[27]+alpha[7]*fUpwind[25]+alpha[6]*fUpwind[24]+alpha[0]*fUpwind[23]+alpha[11]*fUpwind[21]+fUpwind[11]*alpha[21]+alpha[10]*fUpwind[20]+fUpwind[10]*alpha[20]+alpha[14]*fUpwind[18]+fUpwind[14]*alpha[18]+alpha[13]*fUpwind[17]+alpha[1]*fUpwind[15]+alpha[4]*fUpwind[12]+fUpwind[4]*alpha[12]+alpha[5]*fUpwind[9]+fUpwind[5]*alpha[9]); 
  Ghat[24] += 0.1767766952966368*(alpha[7]*fUpwind[31]+alpha[3]*fUpwind[30]+alpha[1]*fUpwind[28]+alpha[18]*fUpwind[27]+alpha[21]*fUpwind[26]+alpha[8]*fUpwind[25]+alpha[0]*fUpwind[24]+alpha[6]*fUpwind[23]+alpha[11]*fUpwind[22]+alpha[9]*fUpwind[20]+fUpwind[9]*alpha[20]+alpha[14]*fUpwind[19]+fUpwind[14]*alpha[19]+alpha[12]*fUpwind[17]+alpha[2]*fUpwind[15]+alpha[4]*fUpwind[13]+fUpwind[4]*alpha[13]+alpha[5]*fUpwind[10]+fUpwind[5]*alpha[10]); 
  Ghat[25] += 0.1767766952966368*(alpha[6]*fUpwind[31]+alpha[2]*fUpwind[30]+alpha[1]*fUpwind[29]+alpha[20]*fUpwind[26]+alpha[0]*fUpwind[25]+alpha[8]*fUpwind[24]+alpha[7]*fUpwind[23]+alpha[10]*fUpwind[22]+alpha[9]*fUpwind[21]+fUpwind[9]*alpha[21]+alpha[13]*fUpwind[19]+fUpwind[13]*alpha[19]+alpha[12]*fUpwind[18]+fUpwind[12]*alpha[18]+alpha[3]*fUpwind[15]+alpha[4]*fUpwind[14]+fUpwind[4]*alpha[14]+alpha[5]*fUpwind[11]+fUpwind[5]*alpha[11]); 
  Ghat[26] += 0.1767766952966368*(alpha[5]*fUpwind[31]+alpha[12]*fUpwind[30]+alpha[13]*fUpwind[29]+alpha[14]*fUpwind[28]+alpha[0]*fUpwind[26]+alpha[20]*fUpwind[25]+alpha[21]*fUpwind[24]+alpha[1]*fUpwind[19]+fUpwind[1]*alpha[19]+alpha[2]*fUpwind[18]+fUpwind[2]*alpha[18]+alpha[3]*fUpwind[17]+alpha[4]*fUpwind[16]+alpha[6]*fUpwind[11]+fUpwind[6]*alpha[11]+alpha[7]*fUpwind[10]+fUpwind[7]*alpha[10]+alpha[8]*fUpwind[9]+fUpwind[8]*alpha[9]); 
  Ghat[27] += 0.1767766952966368*(alpha[4]*fUpwind[31]+alpha[9]*fUpwind[30]+alpha[10]*fUpwind[29]+alpha[11]*fUpwind[28]+alpha[0]*fUpwind[27]+alpha[18]*fUpwind[24]+alpha[19]*fUpwind[23]+alpha[1]*fUpwind[22]+alpha[2]*fUpwind[21]+fUpwind[2]*alpha[21]+alpha[3]*fUpwind[20]+fUpwind[3]*alpha[20]+alpha[5]*fUpwind[16]+alpha[6]*fUpwind[14]+fUpwind[6]*alpha[14]+alpha[7]*fUpwind[13]+fUpwind[7]*alpha[13]+alpha[8]*fUpwind[12]+fUpwind[8]*alpha[12]); 
  Ghat[28] += 0.1767766952966368*(alpha[3]*fUpwind[31]+alpha[7]*fUpwind[30]+alpha[8]*fUpwind[29]+alpha[0]*fUpwind[28]+alpha[11]*fUpwind[27]+alpha[14]*fUpwind[26]+alpha[1]*fUpwind[24]+alpha[2]*fUpwind[23]+alpha[18]*fUpwind[22]+alpha[19]*fUpwind[21]+fUpwind[19]*alpha[21]+alpha[4]*fUpwind[20]+fUpwind[4]*alpha[20]+alpha[5]*fUpwind[17]+alpha[6]*fUpwind[15]+alpha[9]*fUpwind[13]+fUpwind[9]*alpha[13]+alpha[10]*fUpwind[12]+fUpwind[10]*alpha[12]); 
  Ghat[29] += 0.1767766952966368*(alpha[2]*fUpwind[31]+alpha[6]*fUpwind[30]+alpha[0]*fUpwind[29]+alpha[8]*fUpwind[28]+alpha[10]*fUpwind[27]+alpha[13]*fUpwind[26]+alpha[1]*fUpwind[25]+alpha[3]*fUpwind[23]+alpha[4]*fUpwind[21]+fUpwind[4]*alpha[21]+alpha[19]*fUpwind[20]+fUpwind[19]*alpha[20]+alpha[5]*fUpwind[18]+fUpwind[5]*alpha[18]+alpha[7]*fUpwind[15]+alpha[9]*fUpwind[14]+fUpwind[9]*alpha[14]+alpha[11]*fUpwind[12]+fUpwind[11]*alpha[12]); 
  Ghat[30] += 0.1767766952966368*(alpha[1]*fUpwind[31]+alpha[0]*fUpwind[30]+alpha[6]*fUpwind[29]+alpha[7]*fUpwind[28]+alpha[9]*fUpwind[27]+alpha[12]*fUpwind[26]+alpha[2]*fUpwind[25]+alpha[3]*fUpwind[24]+alpha[4]*fUpwind[22]+fUpwind[17]*alpha[21]+alpha[18]*fUpwind[20]+fUpwind[18]*alpha[20]+alpha[5]*fUpwind[19]+fUpwind[5]*alpha[19]+alpha[8]*fUpwind[15]+alpha[10]*fUpwind[14]+fUpwind[10]*alpha[14]+alpha[11]*fUpwind[13]+fUpwind[11]*alpha[13]); 
  Ghat[31] += 0.1767766952966368*(alpha[0]*fUpwind[31]+alpha[1]*fUpwind[30]+alpha[2]*fUpwind[29]+alpha[3]*fUpwind[28]+alpha[4]*fUpwind[27]+alpha[5]*fUpwind[26]+alpha[6]*fUpwind[25]+alpha[7]*fUpwind[24]+alpha[8]*fUpwind[23]+alpha[9]*fUpwind[22]+alpha[10]*fUpwind[21]+fUpwind[10]*alpha[21]+alpha[11]*fUpwind[20]+fUpwind[11]*alpha[20]+alpha[12]*fUpwind[19]+fUpwind[12]*alpha[19]+alpha[13]*fUpwind[18]+fUpwind[13]*alpha[18]+alpha[14]*fUpwind[17]); 

  out[0] += 0.7071067811865475*Ghat[0]*dv11; 
  out[1] += 0.7071067811865475*Ghat[1]*dv11; 
  out[2] += 0.7071067811865475*Ghat[2]*dv11; 
  out[3] += 0.7071067811865475*Ghat[3]*dv11; 
  out[4] += 0.7071067811865475*Ghat[4]*dv11; 
  out[5] += -1.224744871391589*Ghat[0]*dv11; 
  out[6] += 0.7071067811865475*Ghat[5]*dv11; 
  out[7] += 0.7071067811865475*Ghat[6]*dv11; 
  out[8] += 0.7071067811865475*Ghat[7]*dv11; 
  out[9] += 0.7071067811865475*Ghat[8]*dv11; 
  out[10] += 0.7071067811865475*Ghat[9]*dv11; 
  out[11] += 0.7071067811865475*Ghat[10]*dv11; 
  out[12] += 0.7071067811865475*Ghat[11]*dv11; 
  out[13] += -1.224744871391589*Ghat[1]*dv11; 
  out[14] += -1.224744871391589*Ghat[2]*dv11; 
  out[15] += -1.224744871391589*Ghat[3]*dv11; 
  out[16] += -1.224744871391589*Ghat[4]*dv11; 
  out[17] += 0.7071067811865475*Ghat[12]*dv11; 
  out[18] += 0.7071067811865475*Ghat[13]*dv11; 
  out[19] += 0.7071067811865475*Ghat[14]*dv11; 
  out[20] += 0.7071067811865475*Ghat[15]*dv11; 
  out[21] += -1.224744871391589*Ghat[5]*dv11; 
  out[22] += 0.7071067811865475*Ghat[16]*dv11; 
  out[23] += 0.7071067811865475*Ghat[17]*dv11; 
  out[24] += 0.7071067811865475*Ghat[18]*dv11; 
  out[25] += 0.7071067811865475*Ghat[19]*dv11; 
  out[26] += -1.224744871391589*Ghat[6]*dv11; 
  out[27] += -1.224744871391589*Ghat[7]*dv11; 
  out[28] += -1.224744871391589*Ghat[8]*dv11; 
  out[29] += -1.224744871391589*Ghat[9]*dv11; 
  out[30] += -1.224744871391589*Ghat[10]*dv11; 
  out[31] += -1.224744871391589*Ghat[11]*dv11; 
  out[32] += 0.7071067811865475*Ghat[20]*dv11; 
  out[33] += 0.7071067811865475*Ghat[21]*dv11; 
  out[34] += 0.7071067811865475*Ghat[22]*dv11; 
  out[35] += 0.7071067811865475*Ghat[23]*dv11; 
  out[36] += 0.7071067811865475*Ghat[24]*dv11; 
  out[37] += 0.7071067811865475*Ghat[25]*dv11; 
  out[38] += -1.224744871391589*Ghat[12]*dv11; 
  out[39] += -1.224744871391589*Ghat[13]*dv11; 
  out[40] += -1.224744871391589*Ghat[14]*dv11; 
  out[41] += -1.224744871391589*Ghat[15]*dv11; 
  out[42] += 0.7071067811865475*Ghat[26]*dv11; 
  out[43] += -1.224744871391589*Ghat[16]*dv11; 
  out[44] += -1.224744871391589*Ghat[17]*dv11; 
  out[45] += -1.224744871391589*Ghat[18]*dv11; 
  out[46] += -1.224744871391589*Ghat[19]*dv11; 
  out[47] += 0.7071067811865475*Ghat[27]*dv11; 
  out[48] += 0.7071067811865475*Ghat[28]*dv11; 
  out[49] += 0.7071067811865475*Ghat[29]*dv11; 
  out[50] += 0.7071067811865475*Ghat[30]*dv11; 
  out[51] += -1.224744871391589*Ghat[20]*dv11; 
  out[52] += -1.224744871391589*Ghat[21]*dv11; 
  out[53] += -1.224744871391589*Ghat[22]*dv11; 
  out[54] += -1.224744871391589*Ghat[23]*dv11; 
  out[55] += -1.224744871391589*Ghat[24]*dv11; 
  out[56] += -1.224744871391589*Ghat[25]*dv11; 
  out[57] += -1.224744871391589*Ghat[26]*dv11; 
  out[58] += 0.7071067811865475*Ghat[31]*dv11; 
  out[59] += -1.224744871391589*Ghat[27]*dv11; 
  out[60] += -1.224744871391589*Ghat[28]*dv11; 
  out[61] += -1.224744871391589*Ghat[29]*dv11; 
  out[62] += -1.224744871391589*Ghat[30]*dv11; 
  out[63] += -1.224744871391589*Ghat[31]*dv11; 

  } 
} 
