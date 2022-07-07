#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_ser_6x_p1_surfx5_eval_quad.h> 
#include <gkyl_basis_ser_6x_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH void vlasov_boundary_surfvy_3x3v_ser_p1(const double *w, const double *dxv, const double *qmem, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w:           Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // qmem:        q/m*EM fields.
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Output distribution function in skin cell 
  const double dv11 = 2/dxv[4]; 
  const double dv1 = dxv[3], wv1 = w[3]; 
  const double dv2 = dxv[4], wv2 = w[4]; 
  const double dv3 = dxv[5], wv3 = w[5]; 
  const double *E1 = &qmem[8]; 
  const double *B0 = &qmem[24]; 
  const double *B1 = &qmem[32]; 
  const double *B2 = &qmem[40]; 

  double alpha[64] = {0.0}; 

  alpha[0] = 2.0*B0[0]*wv3-2.0*B2[0]*wv1+2.0*E1[0]; 
  alpha[1] = 2.0*B0[1]*wv3-2.0*B2[1]*wv1+2.0*E1[1]; 
  alpha[2] = 2.0*B0[2]*wv3-2.0*B2[2]*wv1+2.0*E1[2]; 
  alpha[3] = 2.0*B0[3]*wv3-2.0*B2[3]*wv1+2.0*E1[3]; 
  alpha[4] = -0.5773502691896258*B2[0]*dv1; 
  alpha[5] = 0.5773502691896258*B0[0]*dv3; 
  alpha[6] = 2.0*B0[4]*wv3-2.0*B2[4]*wv1+2.0*E1[4]; 
  alpha[7] = 2.0*B0[5]*wv3-2.0*B2[5]*wv1+2.0*E1[5]; 
  alpha[8] = 2.0*B0[6]*wv3-2.0*B2[6]*wv1+2.0*E1[6]; 
  alpha[9] = -0.5773502691896258*B2[1]*dv1; 
  alpha[10] = -0.5773502691896258*B2[2]*dv1; 
  alpha[11] = -0.5773502691896258*B2[3]*dv1; 
  alpha[12] = 0.5773502691896258*B0[1]*dv3; 
  alpha[13] = 0.5773502691896258*B0[2]*dv3; 
  alpha[14] = 0.5773502691896258*B0[3]*dv3; 
  alpha[16] = 2.0*B0[7]*wv3-2.0*B2[7]*wv1+2.0*E1[7]; 
  alpha[17] = -0.5773502691896258*B2[4]*dv1; 
  alpha[18] = -0.5773502691896258*B2[5]*dv1; 
  alpha[19] = -0.5773502691896258*B2[6]*dv1; 
  alpha[20] = 0.5773502691896258*B0[4]*dv3; 
  alpha[21] = 0.5773502691896258*B0[5]*dv3; 
  alpha[22] = 0.5773502691896258*B0[6]*dv3; 
  alpha[26] = -0.5773502691896258*B2[7]*dv1; 
  alpha[27] = 0.5773502691896258*B0[7]*dv3; 

  double fUpwindQuad[32] = {0.0};
  double fUpwind[64] = {0.0};
  double Ghat[64] = {0.0}; 

  if (edge == -1) { 

  if (alpha[27]+alpha[26]-alpha[22]-alpha[21]-alpha[20]-alpha[19]-alpha[18]-alpha[17]-alpha[16]+alpha[14]+alpha[13]+alpha[12]+alpha[11]+alpha[10]+alpha[9]+alpha[8]+alpha[7]+alpha[6]-alpha[5]-alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[0] = ser_6x_p1_surfx5_eval_quad_node_0_r(fSkin); 
  } else { 
    fUpwindQuad[0] = ser_6x_p1_surfx5_eval_quad_node_0_l(fEdge); 
  } 
  if ((-alpha[27])+alpha[26]+alpha[22]+alpha[21]+alpha[20]-alpha[19]-alpha[18]-alpha[17]-alpha[16]-alpha[14]-alpha[13]-alpha[12]+alpha[11]+alpha[10]+alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5]-alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[1] = ser_6x_p1_surfx5_eval_quad_node_1_r(fSkin); 
  } else { 
    fUpwindQuad[1] = ser_6x_p1_surfx5_eval_quad_node_1_l(fEdge); 
  } 
  if (alpha[27]-alpha[26]-alpha[22]-alpha[21]-alpha[20]+alpha[19]+alpha[18]+alpha[17]-alpha[16]+alpha[14]+alpha[13]+alpha[12]-alpha[11]-alpha[10]-alpha[9]+alpha[8]+alpha[7]+alpha[6]-alpha[5]+alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[2] = ser_6x_p1_surfx5_eval_quad_node_2_r(fSkin); 
  } else { 
    fUpwindQuad[2] = ser_6x_p1_surfx5_eval_quad_node_2_l(fEdge); 
  } 
  if ((-alpha[27])-alpha[26]+alpha[22]+alpha[21]+alpha[20]+alpha[19]+alpha[18]+alpha[17]-alpha[16]-alpha[14]-alpha[13]-alpha[12]-alpha[11]-alpha[10]-alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5]+alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[3] = ser_6x_p1_surfx5_eval_quad_node_3_r(fSkin); 
  } else { 
    fUpwindQuad[3] = ser_6x_p1_surfx5_eval_quad_node_3_l(fEdge); 
  } 
  if ((-alpha[27])-alpha[26]+alpha[22]+alpha[21]-alpha[20]+alpha[19]+alpha[18]-alpha[17]+alpha[16]-alpha[14]+alpha[13]+alpha[12]-alpha[11]+alpha[10]+alpha[9]-alpha[8]-alpha[7]+alpha[6]-alpha[5]-alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[4] = ser_6x_p1_surfx5_eval_quad_node_4_r(fSkin); 
  } else { 
    fUpwindQuad[4] = ser_6x_p1_surfx5_eval_quad_node_4_l(fEdge); 
  } 
  if (alpha[27]-alpha[26]-alpha[22]-alpha[21]+alpha[20]+alpha[19]+alpha[18]-alpha[17]+alpha[16]+alpha[14]-alpha[13]-alpha[12]-alpha[11]+alpha[10]+alpha[9]-alpha[8]-alpha[7]+alpha[6]+alpha[5]-alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[5] = ser_6x_p1_surfx5_eval_quad_node_5_r(fSkin); 
  } else { 
    fUpwindQuad[5] = ser_6x_p1_surfx5_eval_quad_node_5_l(fEdge); 
  } 
  if ((-alpha[27])+alpha[26]+alpha[22]+alpha[21]-alpha[20]-alpha[19]-alpha[18]+alpha[17]+alpha[16]-alpha[14]+alpha[13]+alpha[12]+alpha[11]-alpha[10]-alpha[9]-alpha[8]-alpha[7]+alpha[6]-alpha[5]+alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[6] = ser_6x_p1_surfx5_eval_quad_node_6_r(fSkin); 
  } else { 
    fUpwindQuad[6] = ser_6x_p1_surfx5_eval_quad_node_6_l(fEdge); 
  } 
  if (alpha[27]+alpha[26]-alpha[22]-alpha[21]+alpha[20]-alpha[19]-alpha[18]+alpha[17]+alpha[16]+alpha[14]-alpha[13]-alpha[12]+alpha[11]-alpha[10]-alpha[9]-alpha[8]-alpha[7]+alpha[6]+alpha[5]+alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[7] = ser_6x_p1_surfx5_eval_quad_node_7_r(fSkin); 
  } else { 
    fUpwindQuad[7] = ser_6x_p1_surfx5_eval_quad_node_7_l(fEdge); 
  } 
  if ((-alpha[27])-alpha[26]+alpha[22]-alpha[21]+alpha[20]+alpha[19]-alpha[18]+alpha[17]+alpha[16]+alpha[14]-alpha[13]+alpha[12]+alpha[11]-alpha[10]+alpha[9]-alpha[8]+alpha[7]-alpha[6]-alpha[5]-alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[8] = ser_6x_p1_surfx5_eval_quad_node_8_r(fSkin); 
  } else { 
    fUpwindQuad[8] = ser_6x_p1_surfx5_eval_quad_node_8_l(fEdge); 
  } 
  if (alpha[27]-alpha[26]-alpha[22]+alpha[21]-alpha[20]+alpha[19]-alpha[18]+alpha[17]+alpha[16]-alpha[14]+alpha[13]-alpha[12]+alpha[11]-alpha[10]+alpha[9]-alpha[8]+alpha[7]-alpha[6]+alpha[5]-alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[9] = ser_6x_p1_surfx5_eval_quad_node_9_r(fSkin); 
  } else { 
    fUpwindQuad[9] = ser_6x_p1_surfx5_eval_quad_node_9_l(fEdge); 
  } 
  if ((-alpha[27])+alpha[26]+alpha[22]-alpha[21]+alpha[20]-alpha[19]+alpha[18]-alpha[17]+alpha[16]+alpha[14]-alpha[13]+alpha[12]-alpha[11]+alpha[10]-alpha[9]-alpha[8]+alpha[7]-alpha[6]-alpha[5]+alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[10] = ser_6x_p1_surfx5_eval_quad_node_10_r(fSkin); 
  } else { 
    fUpwindQuad[10] = ser_6x_p1_surfx5_eval_quad_node_10_l(fEdge); 
  } 
  if (alpha[27]+alpha[26]-alpha[22]+alpha[21]-alpha[20]-alpha[19]+alpha[18]-alpha[17]+alpha[16]-alpha[14]+alpha[13]-alpha[12]-alpha[11]+alpha[10]-alpha[9]-alpha[8]+alpha[7]-alpha[6]+alpha[5]+alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[11] = ser_6x_p1_surfx5_eval_quad_node_11_r(fSkin); 
  } else { 
    fUpwindQuad[11] = ser_6x_p1_surfx5_eval_quad_node_11_l(fEdge); 
  } 
  if (alpha[27]+alpha[26]-alpha[22]+alpha[21]+alpha[20]-alpha[19]+alpha[18]+alpha[17]-alpha[16]-alpha[14]-alpha[13]+alpha[12]-alpha[11]-alpha[10]+alpha[9]+alpha[8]-alpha[7]-alpha[6]-alpha[5]-alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[12] = ser_6x_p1_surfx5_eval_quad_node_12_r(fSkin); 
  } else { 
    fUpwindQuad[12] = ser_6x_p1_surfx5_eval_quad_node_12_l(fEdge); 
  } 
  if ((-alpha[27])+alpha[26]+alpha[22]-alpha[21]-alpha[20]-alpha[19]+alpha[18]+alpha[17]-alpha[16]+alpha[14]+alpha[13]-alpha[12]-alpha[11]-alpha[10]+alpha[9]+alpha[8]-alpha[7]-alpha[6]+alpha[5]-alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[13] = ser_6x_p1_surfx5_eval_quad_node_13_r(fSkin); 
  } else { 
    fUpwindQuad[13] = ser_6x_p1_surfx5_eval_quad_node_13_l(fEdge); 
  } 
  if (alpha[27]-alpha[26]-alpha[22]+alpha[21]+alpha[20]+alpha[19]-alpha[18]-alpha[17]-alpha[16]-alpha[14]-alpha[13]+alpha[12]+alpha[11]+alpha[10]-alpha[9]+alpha[8]-alpha[7]-alpha[6]-alpha[5]+alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[14] = ser_6x_p1_surfx5_eval_quad_node_14_r(fSkin); 
  } else { 
    fUpwindQuad[14] = ser_6x_p1_surfx5_eval_quad_node_14_l(fEdge); 
  } 
  if ((-alpha[27])-alpha[26]+alpha[22]-alpha[21]-alpha[20]+alpha[19]-alpha[18]-alpha[17]-alpha[16]+alpha[14]+alpha[13]-alpha[12]+alpha[11]+alpha[10]-alpha[9]+alpha[8]-alpha[7]-alpha[6]+alpha[5]+alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[15] = ser_6x_p1_surfx5_eval_quad_node_15_r(fSkin); 
  } else { 
    fUpwindQuad[15] = ser_6x_p1_surfx5_eval_quad_node_15_l(fEdge); 
  } 
  if ((-alpha[27])-alpha[26]-alpha[22]+alpha[21]+alpha[20]-alpha[19]+alpha[18]+alpha[17]+alpha[16]+alpha[14]+alpha[13]-alpha[12]+alpha[11]+alpha[10]-alpha[9]+alpha[8]-alpha[7]-alpha[6]-alpha[5]-alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[16] = ser_6x_p1_surfx5_eval_quad_node_16_r(fSkin); 
  } else { 
    fUpwindQuad[16] = ser_6x_p1_surfx5_eval_quad_node_16_l(fEdge); 
  } 
  if (alpha[27]-alpha[26]+alpha[22]-alpha[21]-alpha[20]-alpha[19]+alpha[18]+alpha[17]+alpha[16]-alpha[14]-alpha[13]+alpha[12]+alpha[11]+alpha[10]-alpha[9]+alpha[8]-alpha[7]-alpha[6]+alpha[5]-alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[17] = ser_6x_p1_surfx5_eval_quad_node_17_r(fSkin); 
  } else { 
    fUpwindQuad[17] = ser_6x_p1_surfx5_eval_quad_node_17_l(fEdge); 
  } 
  if ((-alpha[27])+alpha[26]-alpha[22]+alpha[21]+alpha[20]+alpha[19]-alpha[18]-alpha[17]+alpha[16]+alpha[14]+alpha[13]-alpha[12]-alpha[11]-alpha[10]+alpha[9]+alpha[8]-alpha[7]-alpha[6]-alpha[5]+alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[18] = ser_6x_p1_surfx5_eval_quad_node_18_r(fSkin); 
  } else { 
    fUpwindQuad[18] = ser_6x_p1_surfx5_eval_quad_node_18_l(fEdge); 
  } 
  if (alpha[27]+alpha[26]+alpha[22]-alpha[21]-alpha[20]+alpha[19]-alpha[18]-alpha[17]+alpha[16]-alpha[14]-alpha[13]+alpha[12]-alpha[11]-alpha[10]+alpha[9]+alpha[8]-alpha[7]-alpha[6]+alpha[5]+alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[19] = ser_6x_p1_surfx5_eval_quad_node_19_r(fSkin); 
  } else { 
    fUpwindQuad[19] = ser_6x_p1_surfx5_eval_quad_node_19_l(fEdge); 
  } 
  if (alpha[27]+alpha[26]+alpha[22]-alpha[21]+alpha[20]+alpha[19]-alpha[18]+alpha[17]-alpha[16]-alpha[14]+alpha[13]-alpha[12]-alpha[11]+alpha[10]-alpha[9]-alpha[8]+alpha[7]-alpha[6]-alpha[5]-alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[20] = ser_6x_p1_surfx5_eval_quad_node_20_r(fSkin); 
  } else { 
    fUpwindQuad[20] = ser_6x_p1_surfx5_eval_quad_node_20_l(fEdge); 
  } 
  if ((-alpha[27])+alpha[26]-alpha[22]+alpha[21]-alpha[20]+alpha[19]-alpha[18]+alpha[17]-alpha[16]+alpha[14]-alpha[13]+alpha[12]-alpha[11]+alpha[10]-alpha[9]-alpha[8]+alpha[7]-alpha[6]+alpha[5]-alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[21] = ser_6x_p1_surfx5_eval_quad_node_21_r(fSkin); 
  } else { 
    fUpwindQuad[21] = ser_6x_p1_surfx5_eval_quad_node_21_l(fEdge); 
  } 
  if (alpha[27]-alpha[26]+alpha[22]-alpha[21]+alpha[20]-alpha[19]+alpha[18]-alpha[17]-alpha[16]-alpha[14]+alpha[13]-alpha[12]+alpha[11]-alpha[10]+alpha[9]-alpha[8]+alpha[7]-alpha[6]-alpha[5]+alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[22] = ser_6x_p1_surfx5_eval_quad_node_22_r(fSkin); 
  } else { 
    fUpwindQuad[22] = ser_6x_p1_surfx5_eval_quad_node_22_l(fEdge); 
  } 
  if ((-alpha[27])-alpha[26]-alpha[22]+alpha[21]-alpha[20]-alpha[19]+alpha[18]-alpha[17]-alpha[16]+alpha[14]-alpha[13]+alpha[12]+alpha[11]-alpha[10]+alpha[9]-alpha[8]+alpha[7]-alpha[6]+alpha[5]+alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[23] = ser_6x_p1_surfx5_eval_quad_node_23_r(fSkin); 
  } else { 
    fUpwindQuad[23] = ser_6x_p1_surfx5_eval_quad_node_23_l(fEdge); 
  } 
  if (alpha[27]+alpha[26]+alpha[22]+alpha[21]-alpha[20]+alpha[19]+alpha[18]-alpha[17]-alpha[16]+alpha[14]-alpha[13]-alpha[12]+alpha[11]-alpha[10]-alpha[9]-alpha[8]-alpha[7]+alpha[6]-alpha[5]-alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[24] = ser_6x_p1_surfx5_eval_quad_node_24_r(fSkin); 
  } else { 
    fUpwindQuad[24] = ser_6x_p1_surfx5_eval_quad_node_24_l(fEdge); 
  } 
  if ((-alpha[27])+alpha[26]-alpha[22]-alpha[21]+alpha[20]+alpha[19]+alpha[18]-alpha[17]-alpha[16]-alpha[14]+alpha[13]+alpha[12]+alpha[11]-alpha[10]-alpha[9]-alpha[8]-alpha[7]+alpha[6]+alpha[5]-alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[25] = ser_6x_p1_surfx5_eval_quad_node_25_r(fSkin); 
  } else { 
    fUpwindQuad[25] = ser_6x_p1_surfx5_eval_quad_node_25_l(fEdge); 
  } 
  if (alpha[27]-alpha[26]+alpha[22]+alpha[21]-alpha[20]-alpha[19]-alpha[18]+alpha[17]-alpha[16]+alpha[14]-alpha[13]-alpha[12]-alpha[11]+alpha[10]+alpha[9]-alpha[8]-alpha[7]+alpha[6]-alpha[5]+alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[26] = ser_6x_p1_surfx5_eval_quad_node_26_r(fSkin); 
  } else { 
    fUpwindQuad[26] = ser_6x_p1_surfx5_eval_quad_node_26_l(fEdge); 
  } 
  if ((-alpha[27])-alpha[26]-alpha[22]-alpha[21]+alpha[20]-alpha[19]-alpha[18]+alpha[17]-alpha[16]-alpha[14]+alpha[13]+alpha[12]-alpha[11]+alpha[10]+alpha[9]-alpha[8]-alpha[7]+alpha[6]+alpha[5]+alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[27] = ser_6x_p1_surfx5_eval_quad_node_27_r(fSkin); 
  } else { 
    fUpwindQuad[27] = ser_6x_p1_surfx5_eval_quad_node_27_l(fEdge); 
  } 
  if ((-alpha[27])-alpha[26]-alpha[22]-alpha[21]-alpha[20]-alpha[19]-alpha[18]-alpha[17]+alpha[16]-alpha[14]-alpha[13]-alpha[12]-alpha[11]-alpha[10]-alpha[9]+alpha[8]+alpha[7]+alpha[6]-alpha[5]-alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[28] = ser_6x_p1_surfx5_eval_quad_node_28_r(fSkin); 
  } else { 
    fUpwindQuad[28] = ser_6x_p1_surfx5_eval_quad_node_28_l(fEdge); 
  } 
  if (alpha[27]-alpha[26]+alpha[22]+alpha[21]+alpha[20]-alpha[19]-alpha[18]-alpha[17]+alpha[16]+alpha[14]+alpha[13]+alpha[12]-alpha[11]-alpha[10]-alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5]-alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[29] = ser_6x_p1_surfx5_eval_quad_node_29_r(fSkin); 
  } else { 
    fUpwindQuad[29] = ser_6x_p1_surfx5_eval_quad_node_29_l(fEdge); 
  } 
  if ((-alpha[27])+alpha[26]-alpha[22]-alpha[21]-alpha[20]+alpha[19]+alpha[18]+alpha[17]+alpha[16]-alpha[14]-alpha[13]-alpha[12]+alpha[11]+alpha[10]+alpha[9]+alpha[8]+alpha[7]+alpha[6]-alpha[5]+alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[30] = ser_6x_p1_surfx5_eval_quad_node_30_r(fSkin); 
  } else { 
    fUpwindQuad[30] = ser_6x_p1_surfx5_eval_quad_node_30_l(fEdge); 
  } 
  if (alpha[27]+alpha[26]+alpha[22]+alpha[21]+alpha[20]+alpha[19]+alpha[18]+alpha[17]+alpha[16]+alpha[14]+alpha[13]+alpha[12]+alpha[11]+alpha[10]+alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5]+alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[31] = ser_6x_p1_surfx5_eval_quad_node_31_r(fSkin); 
  } else { 
    fUpwindQuad[31] = ser_6x_p1_surfx5_eval_quad_node_31_l(fEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_6x_p1_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.1767766952966368*(alpha[27]*fUpwind[27]+alpha[26]*fUpwind[26]+alpha[22]*fUpwind[22]+alpha[21]*fUpwind[21]+alpha[20]*fUpwind[20]+alpha[19]*fUpwind[19]+alpha[18]*fUpwind[18]+alpha[17]*fUpwind[17]+alpha[16]*fUpwind[16]+alpha[14]*fUpwind[14]+alpha[13]*fUpwind[13]+alpha[12]*fUpwind[12]+alpha[11]*fUpwind[11]+alpha[10]*fUpwind[10]+alpha[9]*fUpwind[9]+alpha[8]*fUpwind[8]+alpha[7]*fUpwind[7]+alpha[6]*fUpwind[6]+alpha[5]*fUpwind[5]+alpha[4]*fUpwind[4]+alpha[3]*fUpwind[3]+alpha[2]*fUpwind[2]+alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] = 0.1767766952966368*(alpha[22]*fUpwind[27]+fUpwind[22]*alpha[27]+alpha[19]*fUpwind[26]+fUpwind[19]*alpha[26]+alpha[14]*fUpwind[21]+fUpwind[14]*alpha[21]+alpha[13]*fUpwind[20]+fUpwind[13]*alpha[20]+alpha[11]*fUpwind[18]+fUpwind[11]*alpha[18]+alpha[10]*fUpwind[17]+fUpwind[10]*alpha[17]+alpha[8]*fUpwind[16]+fUpwind[8]*alpha[16]+alpha[5]*fUpwind[12]+fUpwind[5]*alpha[12]+alpha[4]*fUpwind[9]+fUpwind[4]*alpha[9]+alpha[3]*fUpwind[7]+fUpwind[3]*alpha[7]+alpha[2]*fUpwind[6]+fUpwind[2]*alpha[6]+alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]); 
  Ghat[2] = 0.1767766952966368*(alpha[21]*fUpwind[27]+fUpwind[21]*alpha[27]+alpha[18]*fUpwind[26]+fUpwind[18]*alpha[26]+alpha[14]*fUpwind[22]+fUpwind[14]*alpha[22]+alpha[12]*fUpwind[20]+fUpwind[12]*alpha[20]+alpha[11]*fUpwind[19]+fUpwind[11]*alpha[19]+alpha[9]*fUpwind[17]+fUpwind[9]*alpha[17]+alpha[7]*fUpwind[16]+fUpwind[7]*alpha[16]+alpha[5]*fUpwind[13]+fUpwind[5]*alpha[13]+alpha[4]*fUpwind[10]+fUpwind[4]*alpha[10]+alpha[3]*fUpwind[8]+fUpwind[3]*alpha[8]+alpha[1]*fUpwind[6]+fUpwind[1]*alpha[6]+alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2]); 
  Ghat[3] = 0.1767766952966368*(alpha[20]*fUpwind[27]+fUpwind[20]*alpha[27]+alpha[17]*fUpwind[26]+fUpwind[17]*alpha[26]+alpha[13]*fUpwind[22]+fUpwind[13]*alpha[22]+alpha[12]*fUpwind[21]+fUpwind[12]*alpha[21]+alpha[10]*fUpwind[19]+fUpwind[10]*alpha[19]+alpha[9]*fUpwind[18]+fUpwind[9]*alpha[18]+alpha[6]*fUpwind[16]+fUpwind[6]*alpha[16]+alpha[5]*fUpwind[14]+fUpwind[5]*alpha[14]+alpha[4]*fUpwind[11]+fUpwind[4]*alpha[11]+alpha[2]*fUpwind[8]+fUpwind[2]*alpha[8]+alpha[1]*fUpwind[7]+fUpwind[1]*alpha[7]+alpha[0]*fUpwind[3]+fUpwind[0]*alpha[3]); 
  Ghat[4] = 0.01178511301977579*(13.41640786499874*alpha[26]*fUpwind[43]+13.41640786499874*(alpha[19]*fUpwind[39]+alpha[18]*fUpwind[38]+alpha[17]*fUpwind[37])+13.41640786499874*(alpha[11]*fUpwind[35]+alpha[10]*fUpwind[34]+alpha[9]*fUpwind[33])+13.41640786499874*alpha[4]*fUpwind[32]+15.0*(alpha[27]*fUpwind[31]+alpha[22]*fUpwind[30]+alpha[21]*fUpwind[29]+alpha[20]*fUpwind[28]+alpha[16]*fUpwind[26]+fUpwind[16]*alpha[26]+alpha[14]*fUpwind[25]+alpha[13]*fUpwind[24]+alpha[12]*fUpwind[23]+alpha[8]*fUpwind[19]+fUpwind[8]*alpha[19]+alpha[7]*fUpwind[18]+fUpwind[7]*alpha[18]+alpha[6]*fUpwind[17]+fUpwind[6]*alpha[17]+alpha[5]*fUpwind[15]+alpha[3]*fUpwind[11]+fUpwind[3]*alpha[11]+alpha[2]*fUpwind[10]+fUpwind[2]*alpha[10]+alpha[1]*fUpwind[9]+fUpwind[1]*alpha[9]+alpha[0]*fUpwind[4]+fUpwind[0]*alpha[4])); 
  Ghat[5] = 0.01178511301977579*(13.41640786499874*alpha[27]*fUpwind[59]+13.41640786499874*(alpha[22]*fUpwind[55]+alpha[21]*fUpwind[54]+alpha[20]*fUpwind[53])+13.41640786499874*(alpha[14]*fUpwind[51]+alpha[13]*fUpwind[50]+alpha[12]*fUpwind[49])+13.41640786499874*alpha[5]*fUpwind[48]+15.0*(alpha[26]*fUpwind[31]+alpha[19]*fUpwind[30]+alpha[18]*fUpwind[29]+alpha[17]*fUpwind[28]+alpha[16]*fUpwind[27]+fUpwind[16]*alpha[27]+alpha[11]*fUpwind[25]+alpha[10]*fUpwind[24]+alpha[9]*fUpwind[23]+alpha[8]*fUpwind[22]+fUpwind[8]*alpha[22]+alpha[7]*fUpwind[21]+fUpwind[7]*alpha[21]+alpha[6]*fUpwind[20]+fUpwind[6]*alpha[20]+alpha[4]*fUpwind[15]+alpha[3]*fUpwind[14]+fUpwind[3]*alpha[14]+alpha[2]*fUpwind[13]+fUpwind[2]*alpha[13]+alpha[1]*fUpwind[12]+fUpwind[1]*alpha[12]+alpha[0]*fUpwind[5]+fUpwind[0]*alpha[5])); 
  Ghat[6] = 0.1767766952966368*(alpha[14]*fUpwind[27]+fUpwind[14]*alpha[27]+alpha[11]*fUpwind[26]+fUpwind[11]*alpha[26]+alpha[21]*fUpwind[22]+fUpwind[21]*alpha[22]+alpha[5]*fUpwind[20]+fUpwind[5]*alpha[20]+alpha[18]*fUpwind[19]+fUpwind[18]*alpha[19]+alpha[4]*fUpwind[17]+fUpwind[4]*alpha[17]+alpha[3]*fUpwind[16]+fUpwind[3]*alpha[16]+alpha[12]*fUpwind[13]+fUpwind[12]*alpha[13]+alpha[9]*fUpwind[10]+fUpwind[9]*alpha[10]+alpha[7]*fUpwind[8]+fUpwind[7]*alpha[8]+alpha[0]*fUpwind[6]+fUpwind[0]*alpha[6]+alpha[1]*fUpwind[2]+fUpwind[1]*alpha[2]); 
  Ghat[7] = 0.1767766952966368*(alpha[13]*fUpwind[27]+fUpwind[13]*alpha[27]+alpha[10]*fUpwind[26]+fUpwind[10]*alpha[26]+alpha[20]*fUpwind[22]+fUpwind[20]*alpha[22]+alpha[5]*fUpwind[21]+fUpwind[5]*alpha[21]+alpha[17]*fUpwind[19]+fUpwind[17]*alpha[19]+alpha[4]*fUpwind[18]+fUpwind[4]*alpha[18]+alpha[2]*fUpwind[16]+fUpwind[2]*alpha[16]+alpha[12]*fUpwind[14]+fUpwind[12]*alpha[14]+alpha[9]*fUpwind[11]+fUpwind[9]*alpha[11]+alpha[6]*fUpwind[8]+fUpwind[6]*alpha[8]+alpha[0]*fUpwind[7]+fUpwind[0]*alpha[7]+alpha[1]*fUpwind[3]+fUpwind[1]*alpha[3]); 
  Ghat[8] = 0.1767766952966368*(alpha[12]*fUpwind[27]+fUpwind[12]*alpha[27]+alpha[9]*fUpwind[26]+fUpwind[9]*alpha[26]+alpha[5]*fUpwind[22]+fUpwind[5]*alpha[22]+alpha[20]*fUpwind[21]+fUpwind[20]*alpha[21]+alpha[4]*fUpwind[19]+fUpwind[4]*alpha[19]+alpha[17]*fUpwind[18]+fUpwind[17]*alpha[18]+alpha[1]*fUpwind[16]+fUpwind[1]*alpha[16]+alpha[13]*fUpwind[14]+fUpwind[13]*alpha[14]+alpha[10]*fUpwind[11]+fUpwind[10]*alpha[11]+alpha[0]*fUpwind[8]+fUpwind[0]*alpha[8]+alpha[6]*fUpwind[7]+fUpwind[6]*alpha[7]+alpha[2]*fUpwind[3]+fUpwind[2]*alpha[3]); 
  Ghat[9] = 0.01178511301977579*(13.41640786499874*alpha[19]*fUpwind[43]+13.41640786499874*(alpha[26]*fUpwind[39]+alpha[11]*fUpwind[38]+alpha[10]*fUpwind[37])+13.41640786499874*(alpha[18]*fUpwind[35]+alpha[17]*fUpwind[34]+alpha[4]*fUpwind[33])+13.41640786499874*alpha[9]*fUpwind[32]+15.0*(alpha[22]*fUpwind[31]+alpha[27]*fUpwind[30]+alpha[14]*fUpwind[29]+alpha[13]*fUpwind[28]+alpha[8]*fUpwind[26]+fUpwind[8]*alpha[26]+alpha[21]*fUpwind[25]+alpha[20]*fUpwind[24]+alpha[5]*fUpwind[23]+alpha[16]*fUpwind[19]+fUpwind[16]*alpha[19]+alpha[3]*fUpwind[18]+fUpwind[3]*alpha[18]+alpha[2]*fUpwind[17]+fUpwind[2]*alpha[17]+alpha[12]*fUpwind[15]+alpha[7]*fUpwind[11]+fUpwind[7]*alpha[11]+alpha[6]*fUpwind[10]+fUpwind[6]*alpha[10]+alpha[0]*fUpwind[9]+fUpwind[0]*alpha[9]+alpha[1]*fUpwind[4]+fUpwind[1]*alpha[4])); 
  Ghat[10] = 0.01178511301977579*(13.41640786499874*alpha[18]*fUpwind[43]+13.41640786499874*(alpha[11]*fUpwind[39]+alpha[26]*fUpwind[38]+alpha[9]*fUpwind[37])+13.41640786499874*(alpha[19]*fUpwind[35]+alpha[4]*fUpwind[34]+alpha[17]*fUpwind[33])+13.41640786499874*alpha[10]*fUpwind[32]+15.0*(alpha[21]*fUpwind[31]+alpha[14]*fUpwind[30]+alpha[27]*fUpwind[29]+alpha[12]*fUpwind[28]+alpha[7]*fUpwind[26]+fUpwind[7]*alpha[26]+alpha[22]*fUpwind[25]+alpha[5]*fUpwind[24]+alpha[20]*fUpwind[23]+alpha[3]*fUpwind[19]+fUpwind[3]*alpha[19]+alpha[16]*fUpwind[18]+fUpwind[16]*alpha[18]+alpha[1]*fUpwind[17]+fUpwind[1]*alpha[17]+alpha[13]*fUpwind[15]+alpha[8]*fUpwind[11]+fUpwind[8]*alpha[11]+alpha[0]*fUpwind[10]+fUpwind[0]*alpha[10]+alpha[6]*fUpwind[9]+fUpwind[6]*alpha[9]+alpha[2]*fUpwind[4]+fUpwind[2]*alpha[4])); 
  Ghat[11] = 0.01178511301977579*(13.41640786499874*alpha[17]*fUpwind[43]+13.41640786499874*(alpha[10]*fUpwind[39]+alpha[9]*fUpwind[38]+alpha[26]*fUpwind[37])+13.41640786499874*(alpha[4]*fUpwind[35]+alpha[19]*fUpwind[34]+alpha[18]*fUpwind[33])+13.41640786499874*alpha[11]*fUpwind[32]+15.0*(alpha[20]*fUpwind[31]+alpha[13]*fUpwind[30]+alpha[12]*fUpwind[29]+alpha[27]*fUpwind[28]+alpha[6]*fUpwind[26]+fUpwind[6]*alpha[26]+alpha[5]*fUpwind[25]+alpha[22]*fUpwind[24]+alpha[21]*fUpwind[23]+alpha[2]*fUpwind[19]+fUpwind[2]*alpha[19]+alpha[1]*fUpwind[18]+fUpwind[1]*alpha[18]+alpha[16]*fUpwind[17]+fUpwind[16]*alpha[17]+alpha[14]*fUpwind[15]+alpha[0]*fUpwind[11]+fUpwind[0]*alpha[11]+alpha[8]*fUpwind[10]+fUpwind[8]*alpha[10]+alpha[7]*fUpwind[9]+fUpwind[7]*alpha[9]+alpha[3]*fUpwind[4]+fUpwind[3]*alpha[4])); 
  Ghat[12] = 0.01178511301977579*(13.41640786499874*alpha[22]*fUpwind[59]+13.41640786499874*(alpha[27]*fUpwind[55]+alpha[14]*fUpwind[54]+alpha[13]*fUpwind[53])+13.41640786499874*(alpha[21]*fUpwind[51]+alpha[20]*fUpwind[50]+alpha[5]*fUpwind[49])+13.41640786499874*alpha[12]*fUpwind[48]+15.0*(alpha[19]*fUpwind[31]+alpha[26]*fUpwind[30]+alpha[11]*fUpwind[29]+alpha[10]*fUpwind[28]+alpha[8]*fUpwind[27]+fUpwind[8]*alpha[27]+alpha[18]*fUpwind[25]+alpha[17]*fUpwind[24]+alpha[4]*fUpwind[23]+alpha[16]*fUpwind[22]+fUpwind[16]*alpha[22]+alpha[3]*fUpwind[21]+fUpwind[3]*alpha[21]+alpha[2]*fUpwind[20]+fUpwind[2]*alpha[20]+alpha[9]*fUpwind[15]+alpha[7]*fUpwind[14]+fUpwind[7]*alpha[14]+alpha[6]*fUpwind[13]+fUpwind[6]*alpha[13]+alpha[0]*fUpwind[12]+fUpwind[0]*alpha[12]+alpha[1]*fUpwind[5]+fUpwind[1]*alpha[5])); 
  Ghat[13] = 0.01178511301977579*(13.41640786499874*alpha[21]*fUpwind[59]+13.41640786499874*(alpha[14]*fUpwind[55]+alpha[27]*fUpwind[54]+alpha[12]*fUpwind[53])+13.41640786499874*(alpha[22]*fUpwind[51]+alpha[5]*fUpwind[50]+alpha[20]*fUpwind[49])+13.41640786499874*alpha[13]*fUpwind[48]+15.0*(alpha[18]*fUpwind[31]+alpha[11]*fUpwind[30]+alpha[26]*fUpwind[29]+alpha[9]*fUpwind[28]+alpha[7]*fUpwind[27]+fUpwind[7]*alpha[27]+alpha[19]*fUpwind[25]+alpha[4]*fUpwind[24]+alpha[17]*fUpwind[23]+alpha[3]*fUpwind[22]+fUpwind[3]*alpha[22]+alpha[16]*fUpwind[21]+fUpwind[16]*alpha[21]+alpha[1]*fUpwind[20]+fUpwind[1]*alpha[20]+alpha[10]*fUpwind[15]+alpha[8]*fUpwind[14]+fUpwind[8]*alpha[14]+alpha[0]*fUpwind[13]+fUpwind[0]*alpha[13]+alpha[6]*fUpwind[12]+fUpwind[6]*alpha[12]+alpha[2]*fUpwind[5]+fUpwind[2]*alpha[5])); 
  Ghat[14] = 0.01178511301977579*(13.41640786499874*alpha[20]*fUpwind[59]+13.41640786499874*(alpha[13]*fUpwind[55]+alpha[12]*fUpwind[54]+alpha[27]*fUpwind[53])+13.41640786499874*(alpha[5]*fUpwind[51]+alpha[22]*fUpwind[50]+alpha[21]*fUpwind[49])+13.41640786499874*alpha[14]*fUpwind[48]+15.0*(alpha[17]*fUpwind[31]+alpha[10]*fUpwind[30]+alpha[9]*fUpwind[29]+alpha[26]*fUpwind[28]+alpha[6]*fUpwind[27]+fUpwind[6]*alpha[27]+alpha[4]*fUpwind[25]+alpha[19]*fUpwind[24]+alpha[18]*fUpwind[23]+alpha[2]*fUpwind[22]+fUpwind[2]*alpha[22]+alpha[1]*fUpwind[21]+fUpwind[1]*alpha[21]+alpha[16]*fUpwind[20]+fUpwind[16]*alpha[20]+alpha[11]*fUpwind[15]+alpha[0]*fUpwind[14]+fUpwind[0]*alpha[14]+alpha[8]*fUpwind[13]+fUpwind[8]*alpha[13]+alpha[7]*fUpwind[12]+fUpwind[7]*alpha[12]+alpha[3]*fUpwind[5]+fUpwind[3]*alpha[5])); 
  Ghat[15] = 0.01178511301977579*(13.41640786499874*alpha[27]*fUpwind[63]+13.41640786499874*(alpha[22]*fUpwind[62]+alpha[21]*fUpwind[61]+alpha[20]*fUpwind[60])+13.41640786499874*(alpha[14]*fUpwind[58]+alpha[13]*fUpwind[57]+alpha[12]*fUpwind[56])+13.41640786499874*alpha[5]*fUpwind[52]+13.41640786499874*alpha[26]*fUpwind[47]+13.41640786499874*(alpha[19]*fUpwind[46]+alpha[18]*fUpwind[45]+alpha[17]*fUpwind[44])+13.41640786499874*(alpha[11]*fUpwind[42]+alpha[10]*fUpwind[41]+alpha[9]*fUpwind[40])+13.41640786499874*alpha[4]*fUpwind[36]+15.0*(alpha[16]*fUpwind[31]+alpha[8]*fUpwind[30]+alpha[7]*fUpwind[29]+alpha[6]*fUpwind[28]+alpha[26]*fUpwind[27]+fUpwind[26]*alpha[27]+alpha[3]*fUpwind[25]+alpha[2]*fUpwind[24]+alpha[1]*fUpwind[23]+alpha[19]*fUpwind[22]+fUpwind[19]*alpha[22]+alpha[18]*fUpwind[21]+fUpwind[18]*alpha[21]+alpha[17]*fUpwind[20]+fUpwind[17]*alpha[20]+alpha[0]*fUpwind[15]+alpha[11]*fUpwind[14]+fUpwind[11]*alpha[14]+alpha[10]*fUpwind[13]+fUpwind[10]*alpha[13]+alpha[9]*fUpwind[12]+fUpwind[9]*alpha[12]+alpha[4]*fUpwind[5]+fUpwind[4]*alpha[5])); 
  Ghat[16] = 0.1767766952966368*(alpha[5]*fUpwind[27]+fUpwind[5]*alpha[27]+alpha[4]*fUpwind[26]+fUpwind[4]*alpha[26]+alpha[12]*fUpwind[22]+fUpwind[12]*alpha[22]+alpha[13]*fUpwind[21]+fUpwind[13]*alpha[21]+alpha[14]*fUpwind[20]+fUpwind[14]*alpha[20]+alpha[9]*fUpwind[19]+fUpwind[9]*alpha[19]+alpha[10]*fUpwind[18]+fUpwind[10]*alpha[18]+alpha[11]*fUpwind[17]+fUpwind[11]*alpha[17]+alpha[0]*fUpwind[16]+fUpwind[0]*alpha[16]+alpha[1]*fUpwind[8]+fUpwind[1]*alpha[8]+alpha[2]*fUpwind[7]+fUpwind[2]*alpha[7]+alpha[3]*fUpwind[6]+fUpwind[3]*alpha[6]); 
  Ghat[17] = 0.01178511301977579*(13.41640786499874*alpha[11]*fUpwind[43]+13.41640786499874*(alpha[18]*fUpwind[39]+alpha[19]*fUpwind[38]+alpha[4]*fUpwind[37])+13.41640786499874*(alpha[26]*fUpwind[35]+alpha[9]*fUpwind[34]+alpha[10]*fUpwind[33])+13.41640786499874*alpha[17]*fUpwind[32]+15.0*(alpha[14]*fUpwind[31]+alpha[21]*fUpwind[30]+alpha[22]*fUpwind[29]+alpha[5]*fUpwind[28]+fUpwind[25]*alpha[27]+alpha[3]*fUpwind[26]+fUpwind[3]*alpha[26]+alpha[12]*fUpwind[24]+alpha[13]*fUpwind[23]+fUpwind[15]*alpha[20]+alpha[7]*fUpwind[19]+fUpwind[7]*alpha[19]+alpha[8]*fUpwind[18]+fUpwind[8]*alpha[18]+alpha[0]*fUpwind[17]+fUpwind[0]*alpha[17]+alpha[11]*fUpwind[16]+fUpwind[11]*alpha[16]+alpha[1]*fUpwind[10]+fUpwind[1]*alpha[10]+alpha[2]*fUpwind[9]+fUpwind[2]*alpha[9]+alpha[4]*fUpwind[6]+fUpwind[4]*alpha[6])); 
  Ghat[18] = 0.01178511301977579*(13.41640786499874*alpha[10]*fUpwind[43]+13.41640786499874*(alpha[17]*fUpwind[39]+alpha[4]*fUpwind[38]+alpha[19]*fUpwind[37])+13.41640786499874*(alpha[9]*fUpwind[35]+alpha[26]*fUpwind[34]+alpha[11]*fUpwind[33])+13.41640786499874*alpha[18]*fUpwind[32]+15.0*(alpha[13]*fUpwind[31]+alpha[20]*fUpwind[30]+alpha[5]*fUpwind[29]+alpha[22]*fUpwind[28]+fUpwind[24]*alpha[27]+alpha[2]*fUpwind[26]+fUpwind[2]*alpha[26]+alpha[12]*fUpwind[25]+alpha[14]*fUpwind[23]+fUpwind[15]*alpha[21]+alpha[6]*fUpwind[19]+fUpwind[6]*alpha[19]+alpha[0]*fUpwind[18]+fUpwind[0]*alpha[18]+alpha[8]*fUpwind[17]+fUpwind[8]*alpha[17]+alpha[10]*fUpwind[16]+fUpwind[10]*alpha[16]+alpha[1]*fUpwind[11]+fUpwind[1]*alpha[11]+alpha[3]*fUpwind[9]+fUpwind[3]*alpha[9]+alpha[4]*fUpwind[7]+fUpwind[4]*alpha[7])); 
  Ghat[19] = 0.01178511301977579*(13.41640786499874*alpha[9]*fUpwind[43]+13.41640786499874*(alpha[4]*fUpwind[39]+alpha[17]*fUpwind[38]+alpha[18]*fUpwind[37])+13.41640786499874*(alpha[10]*fUpwind[35]+alpha[11]*fUpwind[34]+alpha[26]*fUpwind[33])+13.41640786499874*alpha[19]*fUpwind[32]+15.0*(alpha[12]*fUpwind[31]+alpha[5]*fUpwind[30]+alpha[20]*fUpwind[29]+alpha[21]*fUpwind[28]+fUpwind[23]*alpha[27]+alpha[1]*fUpwind[26]+fUpwind[1]*alpha[26]+alpha[13]*fUpwind[25]+alpha[14]*fUpwind[24]+fUpwind[15]*alpha[22]+alpha[0]*fUpwind[19]+fUpwind[0]*alpha[19]+alpha[6]*fUpwind[18]+fUpwind[6]*alpha[18]+alpha[7]*fUpwind[17]+fUpwind[7]*alpha[17]+alpha[9]*fUpwind[16]+fUpwind[9]*alpha[16]+alpha[2]*fUpwind[11]+fUpwind[2]*alpha[11]+alpha[3]*fUpwind[10]+fUpwind[3]*alpha[10]+alpha[4]*fUpwind[8]+fUpwind[4]*alpha[8])); 
  Ghat[20] = 0.01178511301977579*(13.41640786499874*alpha[14]*fUpwind[59]+13.41640786499874*(alpha[21]*fUpwind[55]+alpha[22]*fUpwind[54]+alpha[5]*fUpwind[53])+13.41640786499874*(alpha[27]*fUpwind[51]+alpha[12]*fUpwind[50]+alpha[13]*fUpwind[49])+13.41640786499874*alpha[20]*fUpwind[48]+15.0*(alpha[11]*fUpwind[31]+alpha[18]*fUpwind[30]+alpha[19]*fUpwind[29]+alpha[4]*fUpwind[28]+alpha[3]*fUpwind[27]+fUpwind[3]*alpha[27]+fUpwind[25]*alpha[26]+alpha[9]*fUpwind[24]+alpha[10]*fUpwind[23]+alpha[7]*fUpwind[22]+fUpwind[7]*alpha[22]+alpha[8]*fUpwind[21]+fUpwind[8]*alpha[21]+alpha[0]*fUpwind[20]+fUpwind[0]*alpha[20]+fUpwind[15]*alpha[17]+alpha[14]*fUpwind[16]+fUpwind[14]*alpha[16]+alpha[1]*fUpwind[13]+fUpwind[1]*alpha[13]+alpha[2]*fUpwind[12]+fUpwind[2]*alpha[12]+alpha[5]*fUpwind[6]+fUpwind[5]*alpha[6])); 
  Ghat[21] = 0.01178511301977579*(13.41640786499874*alpha[13]*fUpwind[59]+13.41640786499874*(alpha[20]*fUpwind[55]+alpha[5]*fUpwind[54]+alpha[22]*fUpwind[53])+13.41640786499874*(alpha[12]*fUpwind[51]+alpha[27]*fUpwind[50]+alpha[14]*fUpwind[49])+13.41640786499874*alpha[21]*fUpwind[48]+15.0*(alpha[10]*fUpwind[31]+alpha[17]*fUpwind[30]+alpha[4]*fUpwind[29]+alpha[19]*fUpwind[28]+alpha[2]*fUpwind[27]+fUpwind[2]*alpha[27]+fUpwind[24]*alpha[26]+alpha[9]*fUpwind[25]+alpha[11]*fUpwind[23]+alpha[6]*fUpwind[22]+fUpwind[6]*alpha[22]+alpha[0]*fUpwind[21]+fUpwind[0]*alpha[21]+alpha[8]*fUpwind[20]+fUpwind[8]*alpha[20]+fUpwind[15]*alpha[18]+alpha[13]*fUpwind[16]+fUpwind[13]*alpha[16]+alpha[1]*fUpwind[14]+fUpwind[1]*alpha[14]+alpha[3]*fUpwind[12]+fUpwind[3]*alpha[12]+alpha[5]*fUpwind[7]+fUpwind[5]*alpha[7])); 
  Ghat[22] = 0.01178511301977579*(13.41640786499874*alpha[12]*fUpwind[59]+13.41640786499874*(alpha[5]*fUpwind[55]+alpha[20]*fUpwind[54]+alpha[21]*fUpwind[53])+13.41640786499874*(alpha[13]*fUpwind[51]+alpha[14]*fUpwind[50]+alpha[27]*fUpwind[49])+13.41640786499874*alpha[22]*fUpwind[48]+15.0*(alpha[9]*fUpwind[31]+alpha[4]*fUpwind[30]+alpha[17]*fUpwind[29]+alpha[18]*fUpwind[28]+alpha[1]*fUpwind[27]+fUpwind[1]*alpha[27]+fUpwind[23]*alpha[26]+alpha[10]*fUpwind[25]+alpha[11]*fUpwind[24]+alpha[0]*fUpwind[22]+fUpwind[0]*alpha[22]+alpha[6]*fUpwind[21]+fUpwind[6]*alpha[21]+alpha[7]*fUpwind[20]+fUpwind[7]*alpha[20]+fUpwind[15]*alpha[19]+alpha[12]*fUpwind[16]+fUpwind[12]*alpha[16]+alpha[2]*fUpwind[14]+fUpwind[2]*alpha[14]+alpha[3]*fUpwind[13]+fUpwind[3]*alpha[13]+alpha[5]*fUpwind[8]+fUpwind[5]*alpha[8])); 
  Ghat[23] = 0.01178511301977579*(13.41640786499874*alpha[22]*fUpwind[63]+13.41640786499874*(alpha[27]*fUpwind[62]+alpha[14]*fUpwind[61]+alpha[13]*fUpwind[60])+13.41640786499874*(alpha[21]*fUpwind[58]+alpha[20]*fUpwind[57]+alpha[5]*fUpwind[56])+13.41640786499874*alpha[12]*fUpwind[52]+13.41640786499874*alpha[19]*fUpwind[47]+13.41640786499874*(alpha[26]*fUpwind[46]+alpha[11]*fUpwind[45]+alpha[10]*fUpwind[44])+13.41640786499874*(alpha[18]*fUpwind[42]+alpha[17]*fUpwind[41]+alpha[4]*fUpwind[40])+13.41640786499874*alpha[9]*fUpwind[36]+15.0*(alpha[8]*fUpwind[31]+alpha[16]*fUpwind[30]+alpha[3]*fUpwind[29]+alpha[2]*fUpwind[28]+alpha[19]*fUpwind[27]+fUpwind[19]*alpha[27]+alpha[22]*fUpwind[26]+fUpwind[22]*alpha[26]+alpha[7]*fUpwind[25]+alpha[6]*fUpwind[24]+alpha[0]*fUpwind[23]+alpha[11]*fUpwind[21]+fUpwind[11]*alpha[21]+alpha[10]*fUpwind[20]+fUpwind[10]*alpha[20]+alpha[14]*fUpwind[18]+fUpwind[14]*alpha[18]+alpha[13]*fUpwind[17]+fUpwind[13]*alpha[17]+alpha[1]*fUpwind[15]+alpha[4]*fUpwind[12]+fUpwind[4]*alpha[12]+alpha[5]*fUpwind[9]+fUpwind[5]*alpha[9])); 
  Ghat[24] = 0.01178511301977579*(13.41640786499874*alpha[21]*fUpwind[63]+13.41640786499874*(alpha[14]*fUpwind[62]+alpha[27]*fUpwind[61]+alpha[12]*fUpwind[60])+13.41640786499874*(alpha[22]*fUpwind[58]+alpha[5]*fUpwind[57]+alpha[20]*fUpwind[56])+13.41640786499874*alpha[13]*fUpwind[52]+13.41640786499874*alpha[18]*fUpwind[47]+13.41640786499874*(alpha[11]*fUpwind[46]+alpha[26]*fUpwind[45]+alpha[9]*fUpwind[44])+13.41640786499874*(alpha[19]*fUpwind[42]+alpha[4]*fUpwind[41]+alpha[17]*fUpwind[40])+13.41640786499874*alpha[10]*fUpwind[36]+15.0*(alpha[7]*fUpwind[31]+alpha[3]*fUpwind[30]+alpha[16]*fUpwind[29]+alpha[1]*fUpwind[28]+alpha[18]*fUpwind[27]+fUpwind[18]*alpha[27]+alpha[21]*fUpwind[26]+fUpwind[21]*alpha[26]+alpha[8]*fUpwind[25]+alpha[0]*fUpwind[24]+alpha[6]*fUpwind[23]+alpha[11]*fUpwind[22]+fUpwind[11]*alpha[22]+alpha[9]*fUpwind[20]+fUpwind[9]*alpha[20]+alpha[14]*fUpwind[19]+fUpwind[14]*alpha[19]+alpha[12]*fUpwind[17]+fUpwind[12]*alpha[17]+alpha[2]*fUpwind[15]+alpha[4]*fUpwind[13]+fUpwind[4]*alpha[13]+alpha[5]*fUpwind[10]+fUpwind[5]*alpha[10])); 
  Ghat[25] = 0.01178511301977579*(13.41640786499874*alpha[20]*fUpwind[63]+13.41640786499874*(alpha[13]*fUpwind[62]+alpha[12]*fUpwind[61]+alpha[27]*fUpwind[60])+13.41640786499874*(alpha[5]*fUpwind[58]+alpha[22]*fUpwind[57]+alpha[21]*fUpwind[56])+13.41640786499874*alpha[14]*fUpwind[52]+13.41640786499874*alpha[17]*fUpwind[47]+13.41640786499874*(alpha[10]*fUpwind[46]+alpha[9]*fUpwind[45]+alpha[26]*fUpwind[44])+13.41640786499874*(alpha[4]*fUpwind[42]+alpha[19]*fUpwind[41]+alpha[18]*fUpwind[40])+13.41640786499874*alpha[11]*fUpwind[36]+15.0*(alpha[6]*fUpwind[31]+alpha[2]*fUpwind[30]+alpha[1]*fUpwind[29]+alpha[16]*fUpwind[28]+alpha[17]*fUpwind[27]+fUpwind[17]*alpha[27]+alpha[20]*fUpwind[26]+fUpwind[20]*alpha[26]+alpha[0]*fUpwind[25]+alpha[8]*fUpwind[24]+alpha[7]*fUpwind[23]+alpha[10]*fUpwind[22]+fUpwind[10]*alpha[22]+alpha[9]*fUpwind[21]+fUpwind[9]*alpha[21]+alpha[13]*fUpwind[19]+fUpwind[13]*alpha[19]+alpha[12]*fUpwind[18]+fUpwind[12]*alpha[18]+alpha[3]*fUpwind[15]+alpha[4]*fUpwind[14]+fUpwind[4]*alpha[14]+alpha[5]*fUpwind[11]+fUpwind[5]*alpha[11])); 
  Ghat[26] = 0.01178511301977579*(13.41640786499874*alpha[4]*fUpwind[43]+13.41640786499874*(alpha[9]*fUpwind[39]+alpha[10]*fUpwind[38]+alpha[11]*fUpwind[37])+13.41640786499874*(alpha[17]*fUpwind[35]+alpha[18]*fUpwind[34]+alpha[19]*fUpwind[33])+13.41640786499874*alpha[26]*fUpwind[32]+15.0*(alpha[5]*fUpwind[31]+alpha[12]*fUpwind[30]+alpha[13]*fUpwind[29]+alpha[14]*fUpwind[28]+fUpwind[15]*alpha[27]+alpha[0]*fUpwind[26]+fUpwind[0]*alpha[26]+alpha[20]*fUpwind[25]+alpha[21]*fUpwind[24]+alpha[22]*fUpwind[23]+alpha[1]*fUpwind[19]+fUpwind[1]*alpha[19]+alpha[2]*fUpwind[18]+fUpwind[2]*alpha[18]+alpha[3]*fUpwind[17]+fUpwind[3]*alpha[17]+alpha[4]*fUpwind[16]+fUpwind[4]*alpha[16]+alpha[6]*fUpwind[11]+fUpwind[6]*alpha[11]+alpha[7]*fUpwind[10]+fUpwind[7]*alpha[10]+alpha[8]*fUpwind[9]+fUpwind[8]*alpha[9])); 
  Ghat[27] = 0.01178511301977579*(13.41640786499874*alpha[5]*fUpwind[59]+13.41640786499874*(alpha[12]*fUpwind[55]+alpha[13]*fUpwind[54]+alpha[14]*fUpwind[53])+13.41640786499874*(alpha[20]*fUpwind[51]+alpha[21]*fUpwind[50]+alpha[22]*fUpwind[49])+13.41640786499874*alpha[27]*fUpwind[48]+15.0*(alpha[4]*fUpwind[31]+alpha[9]*fUpwind[30]+alpha[10]*fUpwind[29]+alpha[11]*fUpwind[28]+alpha[0]*fUpwind[27]+fUpwind[0]*alpha[27]+fUpwind[15]*alpha[26]+alpha[17]*fUpwind[25]+alpha[18]*fUpwind[24]+alpha[19]*fUpwind[23]+alpha[1]*fUpwind[22]+fUpwind[1]*alpha[22]+alpha[2]*fUpwind[21]+fUpwind[2]*alpha[21]+alpha[3]*fUpwind[20]+fUpwind[3]*alpha[20]+alpha[5]*fUpwind[16]+fUpwind[5]*alpha[16]+alpha[6]*fUpwind[14]+fUpwind[6]*alpha[14]+alpha[7]*fUpwind[13]+fUpwind[7]*alpha[13]+alpha[8]*fUpwind[12]+fUpwind[8]*alpha[12])); 
  Ghat[28] = 0.01178511301977579*(13.41640786499874*alpha[14]*fUpwind[63]+13.41640786499874*(alpha[21]*fUpwind[62]+alpha[22]*fUpwind[61]+alpha[5]*fUpwind[60])+13.41640786499874*(alpha[27]*fUpwind[58]+alpha[12]*fUpwind[57]+alpha[13]*fUpwind[56])+13.41640786499874*alpha[20]*fUpwind[52]+13.41640786499874*alpha[11]*fUpwind[47]+13.41640786499874*(alpha[18]*fUpwind[46]+alpha[19]*fUpwind[45]+alpha[4]*fUpwind[44])+13.41640786499874*(alpha[26]*fUpwind[42]+alpha[9]*fUpwind[41]+alpha[10]*fUpwind[40])+13.41640786499874*alpha[17]*fUpwind[36]+15.0*(alpha[3]*fUpwind[31]+alpha[7]*fUpwind[30]+alpha[8]*fUpwind[29]+alpha[0]*fUpwind[28]+alpha[11]*fUpwind[27]+fUpwind[11]*alpha[27]+alpha[14]*fUpwind[26]+fUpwind[14]*alpha[26]+alpha[16]*fUpwind[25]+alpha[1]*fUpwind[24]+alpha[2]*fUpwind[23]+alpha[18]*fUpwind[22]+fUpwind[18]*alpha[22]+alpha[19]*fUpwind[21]+fUpwind[19]*alpha[21]+alpha[4]*fUpwind[20]+fUpwind[4]*alpha[20]+alpha[5]*fUpwind[17]+fUpwind[5]*alpha[17]+alpha[6]*fUpwind[15]+alpha[9]*fUpwind[13]+fUpwind[9]*alpha[13]+alpha[10]*fUpwind[12]+fUpwind[10]*alpha[12])); 
  Ghat[29] = 0.01178511301977579*(13.41640786499874*alpha[13]*fUpwind[63]+13.41640786499874*(alpha[20]*fUpwind[62]+alpha[5]*fUpwind[61]+alpha[22]*fUpwind[60])+13.41640786499874*(alpha[12]*fUpwind[58]+alpha[27]*fUpwind[57]+alpha[14]*fUpwind[56])+13.41640786499874*alpha[21]*fUpwind[52]+13.41640786499874*alpha[10]*fUpwind[47]+13.41640786499874*(alpha[17]*fUpwind[46]+alpha[4]*fUpwind[45]+alpha[19]*fUpwind[44])+13.41640786499874*(alpha[9]*fUpwind[42]+alpha[26]*fUpwind[41]+alpha[11]*fUpwind[40])+13.41640786499874*alpha[18]*fUpwind[36]+15.0*(alpha[2]*fUpwind[31]+alpha[6]*fUpwind[30]+alpha[0]*fUpwind[29]+alpha[8]*fUpwind[28]+alpha[10]*fUpwind[27]+fUpwind[10]*alpha[27]+alpha[13]*fUpwind[26]+fUpwind[13]*alpha[26]+alpha[1]*fUpwind[25]+alpha[16]*fUpwind[24]+alpha[3]*fUpwind[23]+alpha[17]*fUpwind[22]+fUpwind[17]*alpha[22]+alpha[4]*fUpwind[21]+fUpwind[4]*alpha[21]+alpha[19]*fUpwind[20]+fUpwind[19]*alpha[20]+alpha[5]*fUpwind[18]+fUpwind[5]*alpha[18]+alpha[7]*fUpwind[15]+alpha[9]*fUpwind[14]+fUpwind[9]*alpha[14]+alpha[11]*fUpwind[12]+fUpwind[11]*alpha[12])); 
  Ghat[30] = 0.01178511301977579*(13.41640786499874*alpha[12]*fUpwind[63]+13.41640786499874*(alpha[5]*fUpwind[62]+alpha[20]*fUpwind[61]+alpha[21]*fUpwind[60])+13.41640786499874*(alpha[13]*fUpwind[58]+alpha[14]*fUpwind[57]+alpha[27]*fUpwind[56])+13.41640786499874*alpha[22]*fUpwind[52]+13.41640786499874*alpha[9]*fUpwind[47]+13.41640786499874*(alpha[4]*fUpwind[46]+alpha[17]*fUpwind[45]+alpha[18]*fUpwind[44])+13.41640786499874*(alpha[10]*fUpwind[42]+alpha[11]*fUpwind[41]+alpha[26]*fUpwind[40])+13.41640786499874*alpha[19]*fUpwind[36]+15.0*(alpha[1]*fUpwind[31]+alpha[0]*fUpwind[30]+alpha[6]*fUpwind[29]+alpha[7]*fUpwind[28]+alpha[9]*fUpwind[27]+fUpwind[9]*alpha[27]+alpha[12]*fUpwind[26]+fUpwind[12]*alpha[26]+alpha[2]*fUpwind[25]+alpha[3]*fUpwind[24]+alpha[16]*fUpwind[23]+alpha[4]*fUpwind[22]+fUpwind[4]*alpha[22]+alpha[17]*fUpwind[21]+fUpwind[17]*alpha[21]+alpha[18]*fUpwind[20]+fUpwind[18]*alpha[20]+alpha[5]*fUpwind[19]+fUpwind[5]*alpha[19]+alpha[8]*fUpwind[15]+alpha[10]*fUpwind[14]+fUpwind[10]*alpha[14]+alpha[11]*fUpwind[13]+fUpwind[11]*alpha[13])); 
  Ghat[31] = 0.01178511301977579*(13.41640786499874*alpha[5]*fUpwind[63]+13.41640786499874*(alpha[12]*fUpwind[62]+alpha[13]*fUpwind[61]+alpha[14]*fUpwind[60])+13.41640786499874*(alpha[20]*fUpwind[58]+alpha[21]*fUpwind[57]+alpha[22]*fUpwind[56])+13.41640786499874*alpha[27]*fUpwind[52]+13.41640786499874*alpha[4]*fUpwind[47]+13.41640786499874*(alpha[9]*fUpwind[46]+alpha[10]*fUpwind[45]+alpha[11]*fUpwind[44])+13.41640786499874*(alpha[17]*fUpwind[42]+alpha[18]*fUpwind[41]+alpha[19]*fUpwind[40])+13.41640786499874*alpha[26]*fUpwind[36]+15.0*(alpha[0]*fUpwind[31]+alpha[1]*fUpwind[30]+alpha[2]*fUpwind[29]+alpha[3]*fUpwind[28]+alpha[4]*fUpwind[27]+fUpwind[4]*alpha[27]+alpha[5]*fUpwind[26]+fUpwind[5]*alpha[26]+alpha[6]*fUpwind[25]+alpha[7]*fUpwind[24]+alpha[8]*fUpwind[23]+alpha[9]*fUpwind[22]+fUpwind[9]*alpha[22]+alpha[10]*fUpwind[21]+fUpwind[10]*alpha[21]+alpha[11]*fUpwind[20]+fUpwind[11]*alpha[20]+alpha[12]*fUpwind[19]+fUpwind[12]*alpha[19]+alpha[13]*fUpwind[18]+fUpwind[13]*alpha[18]+alpha[14]*fUpwind[17]+fUpwind[14]*alpha[17]+fUpwind[15]*alpha[16])); 
  Ghat[32] = 0.01178511301977579*(15.0*alpha[27]*fUpwind[47]+15.0*(alpha[22]*fUpwind[46]+alpha[21]*fUpwind[45]+alpha[20]*fUpwind[44]+alpha[16]*fUpwind[43])+15.0*(alpha[14]*fUpwind[42]+alpha[13]*fUpwind[41]+alpha[12]*fUpwind[40]+alpha[8]*fUpwind[39]+alpha[7]*fUpwind[38]+alpha[6]*fUpwind[37])+15.0*(alpha[5]*fUpwind[36]+alpha[3]*fUpwind[35]+alpha[2]*fUpwind[34]+alpha[1]*fUpwind[33])+15.0*alpha[0]*fUpwind[32]+13.41640786499874*(alpha[26]*fUpwind[26]+alpha[19]*fUpwind[19]+alpha[18]*fUpwind[18]+alpha[17]*fUpwind[17]+alpha[11]*fUpwind[11]+alpha[10]*fUpwind[10]+alpha[9]*fUpwind[9]+alpha[4]*fUpwind[4])); 
  Ghat[33] = 0.01178511301977579*(15.0*alpha[22]*fUpwind[47]+15.0*(alpha[27]*fUpwind[46]+alpha[14]*fUpwind[45]+alpha[13]*fUpwind[44]+alpha[8]*fUpwind[43])+15.0*(alpha[21]*fUpwind[42]+alpha[20]*fUpwind[41]+alpha[5]*fUpwind[40]+alpha[16]*fUpwind[39]+alpha[3]*fUpwind[38]+alpha[2]*fUpwind[37])+15.0*(alpha[12]*fUpwind[36]+alpha[7]*fUpwind[35]+alpha[6]*fUpwind[34]+alpha[0]*fUpwind[33])+15.0*alpha[1]*fUpwind[32]+13.41640786499874*(alpha[19]*fUpwind[26]+fUpwind[19]*alpha[26]+alpha[11]*fUpwind[18]+fUpwind[11]*alpha[18]+alpha[10]*fUpwind[17]+fUpwind[10]*alpha[17]+alpha[4]*fUpwind[9]+fUpwind[4]*alpha[9])); 
  Ghat[34] = 0.01178511301977579*(15.0*alpha[21]*fUpwind[47]+15.0*(alpha[14]*fUpwind[46]+alpha[27]*fUpwind[45]+alpha[12]*fUpwind[44]+alpha[7]*fUpwind[43])+15.0*(alpha[22]*fUpwind[42]+alpha[5]*fUpwind[41]+alpha[20]*fUpwind[40]+alpha[3]*fUpwind[39]+alpha[16]*fUpwind[38]+alpha[1]*fUpwind[37])+15.0*(alpha[13]*fUpwind[36]+alpha[8]*fUpwind[35]+alpha[0]*fUpwind[34]+alpha[6]*fUpwind[33])+15.0*alpha[2]*fUpwind[32]+13.41640786499874*(alpha[18]*fUpwind[26]+fUpwind[18]*alpha[26]+alpha[11]*fUpwind[19]+fUpwind[11]*alpha[19]+alpha[9]*fUpwind[17]+fUpwind[9]*alpha[17]+alpha[4]*fUpwind[10]+fUpwind[4]*alpha[10])); 
  Ghat[35] = 0.01178511301977579*(15.0*alpha[20]*fUpwind[47]+15.0*(alpha[13]*fUpwind[46]+alpha[12]*fUpwind[45]+alpha[27]*fUpwind[44]+alpha[6]*fUpwind[43])+15.0*(alpha[5]*fUpwind[42]+alpha[22]*fUpwind[41]+alpha[21]*fUpwind[40]+alpha[2]*fUpwind[39]+alpha[1]*fUpwind[38]+alpha[16]*fUpwind[37])+15.0*(alpha[14]*fUpwind[36]+alpha[0]*fUpwind[35]+alpha[8]*fUpwind[34]+alpha[7]*fUpwind[33])+15.0*alpha[3]*fUpwind[32]+13.41640786499874*(alpha[17]*fUpwind[26]+fUpwind[17]*alpha[26]+alpha[10]*fUpwind[19]+fUpwind[10]*alpha[19]+alpha[9]*fUpwind[18]+fUpwind[9]*alpha[18]+alpha[4]*fUpwind[11]+fUpwind[4]*alpha[11])); 
  Ghat[36] = 0.01178511301977579*(15.0*alpha[16]*fUpwind[47]+15.0*(alpha[8]*fUpwind[46]+alpha[7]*fUpwind[45]+alpha[6]*fUpwind[44]+alpha[27]*fUpwind[43])+15.0*(alpha[3]*fUpwind[42]+alpha[2]*fUpwind[41]+alpha[1]*fUpwind[40]+alpha[22]*fUpwind[39]+alpha[21]*fUpwind[38]+alpha[20]*fUpwind[37])+15.0*(alpha[0]*fUpwind[36]+alpha[14]*fUpwind[35]+alpha[13]*fUpwind[34]+alpha[12]*fUpwind[33])+15.0*alpha[5]*fUpwind[32]+13.41640786499874*(alpha[26]*fUpwind[31]+alpha[19]*fUpwind[30]+alpha[18]*fUpwind[29]+alpha[17]*fUpwind[28]+alpha[11]*fUpwind[25]+alpha[10]*fUpwind[24]+alpha[9]*fUpwind[23]+alpha[4]*fUpwind[15])); 
  Ghat[37] = 0.01178511301977579*(15.0*alpha[14]*fUpwind[47]+15.0*(alpha[21]*fUpwind[46]+alpha[22]*fUpwind[45]+alpha[5]*fUpwind[44]+alpha[3]*fUpwind[43])+15.0*(alpha[27]*fUpwind[42]+alpha[12]*fUpwind[41]+alpha[13]*fUpwind[40]+alpha[7]*fUpwind[39]+alpha[8]*fUpwind[38]+alpha[0]*fUpwind[37])+15.0*(alpha[20]*fUpwind[36]+alpha[16]*fUpwind[35]+alpha[1]*fUpwind[34]+alpha[2]*fUpwind[33])+15.0*alpha[6]*fUpwind[32]+13.41640786499874*(alpha[11]*fUpwind[26]+fUpwind[11]*alpha[26]+alpha[18]*fUpwind[19]+fUpwind[18]*alpha[19]+alpha[4]*fUpwind[17]+fUpwind[4]*alpha[17]+alpha[9]*fUpwind[10]+fUpwind[9]*alpha[10])); 
  Ghat[38] = 0.01178511301977579*(15.0*alpha[13]*fUpwind[47]+15.0*(alpha[20]*fUpwind[46]+alpha[5]*fUpwind[45]+alpha[22]*fUpwind[44]+alpha[2]*fUpwind[43])+15.0*(alpha[12]*fUpwind[42]+alpha[27]*fUpwind[41]+alpha[14]*fUpwind[40]+alpha[6]*fUpwind[39]+alpha[0]*fUpwind[38]+alpha[8]*fUpwind[37])+15.0*(alpha[21]*fUpwind[36]+alpha[1]*fUpwind[35]+alpha[16]*fUpwind[34]+alpha[3]*fUpwind[33])+15.0*alpha[7]*fUpwind[32]+13.41640786499874*(alpha[10]*fUpwind[26]+fUpwind[10]*alpha[26]+alpha[17]*fUpwind[19]+fUpwind[17]*alpha[19]+alpha[4]*fUpwind[18]+fUpwind[4]*alpha[18]+alpha[9]*fUpwind[11]+fUpwind[9]*alpha[11])); 
  Ghat[39] = 0.01178511301977579*(15.0*alpha[12]*fUpwind[47]+15.0*(alpha[5]*fUpwind[46]+alpha[20]*fUpwind[45]+alpha[21]*fUpwind[44]+alpha[1]*fUpwind[43])+15.0*(alpha[13]*fUpwind[42]+alpha[14]*fUpwind[41]+alpha[27]*fUpwind[40]+alpha[0]*fUpwind[39]+alpha[6]*fUpwind[38]+alpha[7]*fUpwind[37])+15.0*(alpha[22]*fUpwind[36]+alpha[2]*fUpwind[35]+alpha[3]*fUpwind[34]+alpha[16]*fUpwind[33])+15.0*alpha[8]*fUpwind[32]+13.41640786499874*(alpha[9]*fUpwind[26]+fUpwind[9]*alpha[26]+alpha[4]*fUpwind[19]+fUpwind[4]*alpha[19]+alpha[17]*fUpwind[18]+fUpwind[17]*alpha[18]+alpha[10]*fUpwind[11]+fUpwind[10]*alpha[11])); 
  Ghat[40] = 0.01178511301977579*(15.0*alpha[8]*fUpwind[47]+15.0*(alpha[16]*fUpwind[46]+alpha[3]*fUpwind[45]+alpha[2]*fUpwind[44]+alpha[22]*fUpwind[43])+15.0*(alpha[7]*fUpwind[42]+alpha[6]*fUpwind[41]+alpha[0]*fUpwind[40]+alpha[27]*fUpwind[39]+alpha[14]*fUpwind[38]+alpha[13]*fUpwind[37])+15.0*(alpha[1]*fUpwind[36]+alpha[21]*fUpwind[35]+alpha[20]*fUpwind[34]+alpha[5]*fUpwind[33])+15.0*alpha[12]*fUpwind[32]+13.41640786499874*(alpha[19]*fUpwind[31]+alpha[26]*fUpwind[30]+alpha[11]*fUpwind[29]+alpha[10]*fUpwind[28]+alpha[18]*fUpwind[25]+alpha[17]*fUpwind[24]+alpha[4]*fUpwind[23]+alpha[9]*fUpwind[15])); 
  Ghat[41] = 0.01178511301977579*(15.0*alpha[7]*fUpwind[47]+15.0*(alpha[3]*fUpwind[46]+alpha[16]*fUpwind[45]+alpha[1]*fUpwind[44]+alpha[21]*fUpwind[43])+15.0*(alpha[8]*fUpwind[42]+alpha[0]*fUpwind[41]+alpha[6]*fUpwind[40]+alpha[14]*fUpwind[39]+alpha[27]*fUpwind[38]+alpha[12]*fUpwind[37])+15.0*(alpha[2]*fUpwind[36]+alpha[22]*fUpwind[35]+alpha[5]*fUpwind[34]+alpha[20]*fUpwind[33])+15.0*alpha[13]*fUpwind[32]+13.41640786499874*(alpha[18]*fUpwind[31]+alpha[11]*fUpwind[30]+alpha[26]*fUpwind[29]+alpha[9]*fUpwind[28]+alpha[19]*fUpwind[25]+alpha[4]*fUpwind[24]+alpha[17]*fUpwind[23]+alpha[10]*fUpwind[15])); 
  Ghat[42] = 0.01178511301977579*(15.0*alpha[6]*fUpwind[47]+15.0*(alpha[2]*fUpwind[46]+alpha[1]*fUpwind[45]+alpha[16]*fUpwind[44]+alpha[20]*fUpwind[43])+15.0*(alpha[0]*fUpwind[42]+alpha[8]*fUpwind[41]+alpha[7]*fUpwind[40]+alpha[13]*fUpwind[39]+alpha[12]*fUpwind[38]+alpha[27]*fUpwind[37])+15.0*(alpha[3]*fUpwind[36]+alpha[5]*fUpwind[35]+alpha[22]*fUpwind[34]+alpha[21]*fUpwind[33])+15.0*alpha[14]*fUpwind[32]+13.41640786499874*(alpha[17]*fUpwind[31]+alpha[10]*fUpwind[30]+alpha[9]*fUpwind[29]+alpha[26]*fUpwind[28]+alpha[4]*fUpwind[25]+alpha[19]*fUpwind[24]+alpha[18]*fUpwind[23]+alpha[11]*fUpwind[15])); 
  Ghat[43] = 0.01178511301977579*(15.0*alpha[5]*fUpwind[47]+15.0*(alpha[12]*fUpwind[46]+alpha[13]*fUpwind[45]+alpha[14]*fUpwind[44]+alpha[0]*fUpwind[43])+15.0*(alpha[20]*fUpwind[42]+alpha[21]*fUpwind[41]+alpha[22]*fUpwind[40]+alpha[1]*fUpwind[39]+alpha[2]*fUpwind[38]+alpha[3]*fUpwind[37])+15.0*(alpha[27]*fUpwind[36]+alpha[6]*fUpwind[35]+alpha[7]*fUpwind[34]+alpha[8]*fUpwind[33])+15.0*alpha[16]*fUpwind[32]+13.41640786499874*(alpha[4]*fUpwind[26]+fUpwind[4]*alpha[26]+alpha[9]*fUpwind[19]+fUpwind[9]*alpha[19]+alpha[10]*fUpwind[18]+fUpwind[10]*alpha[18]+alpha[11]*fUpwind[17]+fUpwind[11]*alpha[17])); 
  Ghat[44] = 0.01178511301977579*(15.0*alpha[3]*fUpwind[47]+15.0*(alpha[7]*fUpwind[46]+alpha[8]*fUpwind[45]+alpha[0]*fUpwind[44]+alpha[14]*fUpwind[43])+15.0*(alpha[16]*fUpwind[42]+alpha[1]*fUpwind[41]+alpha[2]*fUpwind[40]+alpha[21]*fUpwind[39]+alpha[22]*fUpwind[38]+alpha[5]*fUpwind[37])+15.0*(alpha[6]*fUpwind[36]+alpha[27]*fUpwind[35]+alpha[12]*fUpwind[34]+alpha[13]*fUpwind[33])+15.0*alpha[20]*fUpwind[32]+13.41640786499874*(alpha[11]*fUpwind[31]+alpha[18]*fUpwind[30]+alpha[19]*fUpwind[29]+alpha[4]*fUpwind[28]+fUpwind[25]*alpha[26]+alpha[9]*fUpwind[24]+alpha[10]*fUpwind[23]+fUpwind[15]*alpha[17])); 
  Ghat[45] = 0.01178511301977579*(15.0*alpha[2]*fUpwind[47]+15.0*(alpha[6]*fUpwind[46]+alpha[0]*fUpwind[45]+alpha[8]*fUpwind[44]+alpha[13]*fUpwind[43])+15.0*(alpha[1]*fUpwind[42]+alpha[16]*fUpwind[41]+alpha[3]*fUpwind[40]+alpha[20]*fUpwind[39]+alpha[5]*fUpwind[38]+alpha[22]*fUpwind[37])+15.0*(alpha[7]*fUpwind[36]+alpha[12]*fUpwind[35]+alpha[27]*fUpwind[34]+alpha[14]*fUpwind[33])+15.0*alpha[21]*fUpwind[32]+13.41640786499874*(alpha[10]*fUpwind[31]+alpha[17]*fUpwind[30]+alpha[4]*fUpwind[29]+alpha[19]*fUpwind[28]+fUpwind[24]*alpha[26]+alpha[9]*fUpwind[25]+alpha[11]*fUpwind[23]+fUpwind[15]*alpha[18])); 
  Ghat[46] = 0.01178511301977579*(15.0*alpha[1]*fUpwind[47]+15.0*(alpha[0]*fUpwind[46]+alpha[6]*fUpwind[45]+alpha[7]*fUpwind[44]+alpha[12]*fUpwind[43])+15.0*(alpha[2]*fUpwind[42]+alpha[3]*fUpwind[41]+alpha[16]*fUpwind[40]+alpha[5]*fUpwind[39]+alpha[20]*fUpwind[38]+alpha[21]*fUpwind[37])+15.0*(alpha[8]*fUpwind[36]+alpha[13]*fUpwind[35]+alpha[14]*fUpwind[34]+alpha[27]*fUpwind[33])+15.0*alpha[22]*fUpwind[32]+13.41640786499874*(alpha[9]*fUpwind[31]+alpha[4]*fUpwind[30]+alpha[17]*fUpwind[29]+alpha[18]*fUpwind[28]+fUpwind[23]*alpha[26]+alpha[10]*fUpwind[25]+alpha[11]*fUpwind[24]+fUpwind[15]*alpha[19])); 
  Ghat[47] = 0.01178511301977579*(15.0*alpha[0]*fUpwind[47]+15.0*(alpha[1]*fUpwind[46]+alpha[2]*fUpwind[45]+alpha[3]*fUpwind[44]+alpha[5]*fUpwind[43])+15.0*(alpha[6]*fUpwind[42]+alpha[7]*fUpwind[41]+alpha[8]*fUpwind[40]+alpha[12]*fUpwind[39]+alpha[13]*fUpwind[38]+alpha[14]*fUpwind[37])+15.0*(alpha[16]*fUpwind[36]+alpha[20]*fUpwind[35]+alpha[21]*fUpwind[34]+alpha[22]*fUpwind[33])+15.0*alpha[27]*fUpwind[32]+13.41640786499874*(alpha[4]*fUpwind[31]+alpha[9]*fUpwind[30]+alpha[10]*fUpwind[29]+alpha[11]*fUpwind[28]+fUpwind[15]*alpha[26]+alpha[17]*fUpwind[25]+alpha[18]*fUpwind[24]+alpha[19]*fUpwind[23])); 
  Ghat[48] = 0.01178511301977579*(15.0*alpha[26]*fUpwind[63]+15.0*(alpha[19]*fUpwind[62]+alpha[18]*fUpwind[61]+alpha[17]*fUpwind[60]+alpha[16]*fUpwind[59])+15.0*(alpha[11]*fUpwind[58]+alpha[10]*fUpwind[57]+alpha[9]*fUpwind[56]+alpha[8]*fUpwind[55]+alpha[7]*fUpwind[54]+alpha[6]*fUpwind[53])+15.0*(alpha[4]*fUpwind[52]+alpha[3]*fUpwind[51]+alpha[2]*fUpwind[50]+alpha[1]*fUpwind[49])+15.0*alpha[0]*fUpwind[48]+13.41640786499874*(alpha[27]*fUpwind[27]+alpha[22]*fUpwind[22]+alpha[21]*fUpwind[21]+alpha[20]*fUpwind[20]+alpha[14]*fUpwind[14]+alpha[13]*fUpwind[13]+alpha[12]*fUpwind[12]+alpha[5]*fUpwind[5])); 
  Ghat[49] = 0.01178511301977579*(15.0*alpha[19]*fUpwind[63]+15.0*(alpha[26]*fUpwind[62]+alpha[11]*fUpwind[61]+alpha[10]*fUpwind[60]+alpha[8]*fUpwind[59])+15.0*(alpha[18]*fUpwind[58]+alpha[17]*fUpwind[57]+alpha[4]*fUpwind[56]+alpha[16]*fUpwind[55]+alpha[3]*fUpwind[54]+alpha[2]*fUpwind[53])+15.0*(alpha[9]*fUpwind[52]+alpha[7]*fUpwind[51]+alpha[6]*fUpwind[50]+alpha[0]*fUpwind[49])+15.0*alpha[1]*fUpwind[48]+13.41640786499874*(alpha[22]*fUpwind[27]+fUpwind[22]*alpha[27]+alpha[14]*fUpwind[21]+fUpwind[14]*alpha[21]+alpha[13]*fUpwind[20]+fUpwind[13]*alpha[20]+alpha[5]*fUpwind[12]+fUpwind[5]*alpha[12])); 
  Ghat[50] = 0.01178511301977579*(15.0*alpha[18]*fUpwind[63]+15.0*(alpha[11]*fUpwind[62]+alpha[26]*fUpwind[61]+alpha[9]*fUpwind[60]+alpha[7]*fUpwind[59])+15.0*(alpha[19]*fUpwind[58]+alpha[4]*fUpwind[57]+alpha[17]*fUpwind[56]+alpha[3]*fUpwind[55]+alpha[16]*fUpwind[54]+alpha[1]*fUpwind[53])+15.0*(alpha[10]*fUpwind[52]+alpha[8]*fUpwind[51]+alpha[0]*fUpwind[50]+alpha[6]*fUpwind[49])+15.0*alpha[2]*fUpwind[48]+13.41640786499874*(alpha[21]*fUpwind[27]+fUpwind[21]*alpha[27]+alpha[14]*fUpwind[22]+fUpwind[14]*alpha[22]+alpha[12]*fUpwind[20]+fUpwind[12]*alpha[20]+alpha[5]*fUpwind[13]+fUpwind[5]*alpha[13])); 
  Ghat[51] = 0.01178511301977579*(15.0*alpha[17]*fUpwind[63]+15.0*(alpha[10]*fUpwind[62]+alpha[9]*fUpwind[61]+alpha[26]*fUpwind[60]+alpha[6]*fUpwind[59])+15.0*(alpha[4]*fUpwind[58]+alpha[19]*fUpwind[57]+alpha[18]*fUpwind[56]+alpha[2]*fUpwind[55]+alpha[1]*fUpwind[54]+alpha[16]*fUpwind[53])+15.0*(alpha[11]*fUpwind[52]+alpha[0]*fUpwind[51]+alpha[8]*fUpwind[50]+alpha[7]*fUpwind[49])+15.0*alpha[3]*fUpwind[48]+13.41640786499874*(alpha[20]*fUpwind[27]+fUpwind[20]*alpha[27]+alpha[13]*fUpwind[22]+fUpwind[13]*alpha[22]+alpha[12]*fUpwind[21]+fUpwind[12]*alpha[21]+alpha[5]*fUpwind[14]+fUpwind[5]*alpha[14])); 
  Ghat[52] = 0.01178511301977579*(15.0*alpha[16]*fUpwind[63]+15.0*(alpha[8]*fUpwind[62]+alpha[7]*fUpwind[61]+alpha[6]*fUpwind[60]+alpha[26]*fUpwind[59])+15.0*(alpha[3]*fUpwind[58]+alpha[2]*fUpwind[57]+alpha[1]*fUpwind[56]+alpha[19]*fUpwind[55]+alpha[18]*fUpwind[54]+alpha[17]*fUpwind[53])+15.0*(alpha[0]*fUpwind[52]+alpha[11]*fUpwind[51]+alpha[10]*fUpwind[50]+alpha[9]*fUpwind[49])+15.0*alpha[4]*fUpwind[48]+13.41640786499874*(alpha[27]*fUpwind[31]+alpha[22]*fUpwind[30]+alpha[21]*fUpwind[29]+alpha[20]*fUpwind[28]+alpha[14]*fUpwind[25]+alpha[13]*fUpwind[24]+alpha[12]*fUpwind[23]+alpha[5]*fUpwind[15])); 
  Ghat[53] = 0.01178511301977579*(15.0*alpha[11]*fUpwind[63]+15.0*(alpha[18]*fUpwind[62]+alpha[19]*fUpwind[61]+alpha[4]*fUpwind[60]+alpha[3]*fUpwind[59])+15.0*(alpha[26]*fUpwind[58]+alpha[9]*fUpwind[57]+alpha[10]*fUpwind[56]+alpha[7]*fUpwind[55]+alpha[8]*fUpwind[54]+alpha[0]*fUpwind[53])+15.0*(alpha[17]*fUpwind[52]+alpha[16]*fUpwind[51]+alpha[1]*fUpwind[50]+alpha[2]*fUpwind[49])+15.0*alpha[6]*fUpwind[48]+13.41640786499874*(alpha[14]*fUpwind[27]+fUpwind[14]*alpha[27]+alpha[21]*fUpwind[22]+fUpwind[21]*alpha[22]+alpha[5]*fUpwind[20]+fUpwind[5]*alpha[20]+alpha[12]*fUpwind[13]+fUpwind[12]*alpha[13])); 
  Ghat[54] = 0.01178511301977579*(15.0*alpha[10]*fUpwind[63]+15.0*(alpha[17]*fUpwind[62]+alpha[4]*fUpwind[61]+alpha[19]*fUpwind[60]+alpha[2]*fUpwind[59])+15.0*(alpha[9]*fUpwind[58]+alpha[26]*fUpwind[57]+alpha[11]*fUpwind[56]+alpha[6]*fUpwind[55]+alpha[0]*fUpwind[54]+alpha[8]*fUpwind[53])+15.0*(alpha[18]*fUpwind[52]+alpha[1]*fUpwind[51]+alpha[16]*fUpwind[50]+alpha[3]*fUpwind[49])+15.0*alpha[7]*fUpwind[48]+13.41640786499874*(alpha[13]*fUpwind[27]+fUpwind[13]*alpha[27]+alpha[20]*fUpwind[22]+fUpwind[20]*alpha[22]+alpha[5]*fUpwind[21]+fUpwind[5]*alpha[21]+alpha[12]*fUpwind[14]+fUpwind[12]*alpha[14])); 
  Ghat[55] = 0.01178511301977579*(15.0*alpha[9]*fUpwind[63]+15.0*(alpha[4]*fUpwind[62]+alpha[17]*fUpwind[61]+alpha[18]*fUpwind[60]+alpha[1]*fUpwind[59])+15.0*(alpha[10]*fUpwind[58]+alpha[11]*fUpwind[57]+alpha[26]*fUpwind[56]+alpha[0]*fUpwind[55]+alpha[6]*fUpwind[54]+alpha[7]*fUpwind[53])+15.0*(alpha[19]*fUpwind[52]+alpha[2]*fUpwind[51]+alpha[3]*fUpwind[50]+alpha[16]*fUpwind[49])+15.0*alpha[8]*fUpwind[48]+13.41640786499874*(alpha[12]*fUpwind[27]+fUpwind[12]*alpha[27]+alpha[5]*fUpwind[22]+fUpwind[5]*alpha[22]+alpha[20]*fUpwind[21]+fUpwind[20]*alpha[21]+alpha[13]*fUpwind[14]+fUpwind[13]*alpha[14])); 
  Ghat[56] = 0.01178511301977579*(15.0*alpha[8]*fUpwind[63]+15.0*(alpha[16]*fUpwind[62]+alpha[3]*fUpwind[61]+alpha[2]*fUpwind[60]+alpha[19]*fUpwind[59])+15.0*(alpha[7]*fUpwind[58]+alpha[6]*fUpwind[57]+alpha[0]*fUpwind[56]+alpha[26]*fUpwind[55]+alpha[11]*fUpwind[54]+alpha[10]*fUpwind[53])+15.0*(alpha[1]*fUpwind[52]+alpha[18]*fUpwind[51]+alpha[17]*fUpwind[50]+alpha[4]*fUpwind[49])+15.0*alpha[9]*fUpwind[48]+13.41640786499874*(alpha[22]*fUpwind[31]+alpha[27]*fUpwind[30]+alpha[14]*fUpwind[29]+alpha[13]*fUpwind[28]+alpha[21]*fUpwind[25]+alpha[20]*fUpwind[24]+alpha[5]*fUpwind[23]+alpha[12]*fUpwind[15])); 
  Ghat[57] = 0.01178511301977579*(15.0*alpha[7]*fUpwind[63]+15.0*(alpha[3]*fUpwind[62]+alpha[16]*fUpwind[61]+alpha[1]*fUpwind[60]+alpha[18]*fUpwind[59])+15.0*(alpha[8]*fUpwind[58]+alpha[0]*fUpwind[57]+alpha[6]*fUpwind[56]+alpha[11]*fUpwind[55]+alpha[26]*fUpwind[54]+alpha[9]*fUpwind[53])+15.0*(alpha[2]*fUpwind[52]+alpha[19]*fUpwind[51]+alpha[4]*fUpwind[50]+alpha[17]*fUpwind[49])+15.0*alpha[10]*fUpwind[48]+13.41640786499874*(alpha[21]*fUpwind[31]+alpha[14]*fUpwind[30]+alpha[27]*fUpwind[29]+alpha[12]*fUpwind[28]+alpha[22]*fUpwind[25]+alpha[5]*fUpwind[24]+alpha[20]*fUpwind[23]+alpha[13]*fUpwind[15])); 
  Ghat[58] = 0.01178511301977579*(15.0*alpha[6]*fUpwind[63]+15.0*(alpha[2]*fUpwind[62]+alpha[1]*fUpwind[61]+alpha[16]*fUpwind[60]+alpha[17]*fUpwind[59])+15.0*(alpha[0]*fUpwind[58]+alpha[8]*fUpwind[57]+alpha[7]*fUpwind[56]+alpha[10]*fUpwind[55]+alpha[9]*fUpwind[54]+alpha[26]*fUpwind[53])+15.0*(alpha[3]*fUpwind[52]+alpha[4]*fUpwind[51]+alpha[19]*fUpwind[50]+alpha[18]*fUpwind[49])+15.0*alpha[11]*fUpwind[48]+13.41640786499874*(alpha[20]*fUpwind[31]+alpha[13]*fUpwind[30]+alpha[12]*fUpwind[29]+alpha[27]*fUpwind[28]+alpha[5]*fUpwind[25]+alpha[22]*fUpwind[24]+alpha[21]*fUpwind[23]+alpha[14]*fUpwind[15])); 
  Ghat[59] = 0.01178511301977579*(15.0*alpha[4]*fUpwind[63]+15.0*(alpha[9]*fUpwind[62]+alpha[10]*fUpwind[61]+alpha[11]*fUpwind[60]+alpha[0]*fUpwind[59])+15.0*(alpha[17]*fUpwind[58]+alpha[18]*fUpwind[57]+alpha[19]*fUpwind[56]+alpha[1]*fUpwind[55]+alpha[2]*fUpwind[54]+alpha[3]*fUpwind[53])+15.0*(alpha[26]*fUpwind[52]+alpha[6]*fUpwind[51]+alpha[7]*fUpwind[50]+alpha[8]*fUpwind[49])+15.0*alpha[16]*fUpwind[48]+13.41640786499874*(alpha[5]*fUpwind[27]+fUpwind[5]*alpha[27]+alpha[12]*fUpwind[22]+fUpwind[12]*alpha[22]+alpha[13]*fUpwind[21]+fUpwind[13]*alpha[21]+alpha[14]*fUpwind[20]+fUpwind[14]*alpha[20])); 
  Ghat[60] = 0.01178511301977579*(15.0*alpha[3]*fUpwind[63]+15.0*(alpha[7]*fUpwind[62]+alpha[8]*fUpwind[61]+alpha[0]*fUpwind[60]+alpha[11]*fUpwind[59])+15.0*(alpha[16]*fUpwind[58]+alpha[1]*fUpwind[57]+alpha[2]*fUpwind[56]+alpha[18]*fUpwind[55]+alpha[19]*fUpwind[54]+alpha[4]*fUpwind[53])+15.0*(alpha[6]*fUpwind[52]+alpha[26]*fUpwind[51]+alpha[9]*fUpwind[50]+alpha[10]*fUpwind[49])+15.0*alpha[17]*fUpwind[48]+13.41640786499874*(alpha[14]*fUpwind[31]+alpha[21]*fUpwind[30]+alpha[22]*fUpwind[29]+alpha[5]*fUpwind[28]+fUpwind[25]*alpha[27]+alpha[12]*fUpwind[24]+alpha[13]*fUpwind[23]+fUpwind[15]*alpha[20])); 
  Ghat[61] = 0.01178511301977579*(15.0*alpha[2]*fUpwind[63]+15.0*(alpha[6]*fUpwind[62]+alpha[0]*fUpwind[61]+alpha[8]*fUpwind[60]+alpha[10]*fUpwind[59])+15.0*(alpha[1]*fUpwind[58]+alpha[16]*fUpwind[57]+alpha[3]*fUpwind[56]+alpha[17]*fUpwind[55]+alpha[4]*fUpwind[54]+alpha[19]*fUpwind[53])+15.0*(alpha[7]*fUpwind[52]+alpha[9]*fUpwind[51]+alpha[26]*fUpwind[50]+alpha[11]*fUpwind[49])+15.0*alpha[18]*fUpwind[48]+13.41640786499874*(alpha[13]*fUpwind[31]+alpha[20]*fUpwind[30]+alpha[5]*fUpwind[29]+alpha[22]*fUpwind[28]+fUpwind[24]*alpha[27]+alpha[12]*fUpwind[25]+alpha[14]*fUpwind[23]+fUpwind[15]*alpha[21])); 
  Ghat[62] = 0.01178511301977579*(15.0*alpha[1]*fUpwind[63]+15.0*(alpha[0]*fUpwind[62]+alpha[6]*fUpwind[61]+alpha[7]*fUpwind[60]+alpha[9]*fUpwind[59])+15.0*(alpha[2]*fUpwind[58]+alpha[3]*fUpwind[57]+alpha[16]*fUpwind[56]+alpha[4]*fUpwind[55]+alpha[17]*fUpwind[54]+alpha[18]*fUpwind[53])+15.0*(alpha[8]*fUpwind[52]+alpha[10]*fUpwind[51]+alpha[11]*fUpwind[50]+alpha[26]*fUpwind[49])+15.0*alpha[19]*fUpwind[48]+13.41640786499874*(alpha[12]*fUpwind[31]+alpha[5]*fUpwind[30]+alpha[20]*fUpwind[29]+alpha[21]*fUpwind[28]+fUpwind[23]*alpha[27]+alpha[13]*fUpwind[25]+alpha[14]*fUpwind[24]+fUpwind[15]*alpha[22])); 
  Ghat[63] = 0.01178511301977579*(15.0*alpha[0]*fUpwind[63]+15.0*(alpha[1]*fUpwind[62]+alpha[2]*fUpwind[61]+alpha[3]*fUpwind[60]+alpha[4]*fUpwind[59])+15.0*(alpha[6]*fUpwind[58]+alpha[7]*fUpwind[57]+alpha[8]*fUpwind[56]+alpha[9]*fUpwind[55]+alpha[10]*fUpwind[54]+alpha[11]*fUpwind[53])+15.0*(alpha[16]*fUpwind[52]+alpha[17]*fUpwind[51]+alpha[18]*fUpwind[50]+alpha[19]*fUpwind[49])+15.0*alpha[26]*fUpwind[48]+13.41640786499874*(alpha[5]*fUpwind[31]+alpha[12]*fUpwind[30]+alpha[13]*fUpwind[29]+alpha[14]*fUpwind[28]+fUpwind[15]*alpha[27]+alpha[20]*fUpwind[25]+alpha[21]*fUpwind[24]+alpha[22]*fUpwind[23])); 

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
  out[64] += -0.7071067811865475*Ghat[32]*dv11; 
  out[65] += -0.7071067811865475*Ghat[33]*dv11; 
  out[66] += -0.7071067811865475*Ghat[34]*dv11; 
  out[67] += -0.7071067811865475*Ghat[35]*dv11; 
  out[68] += -1.224744871391589*Ghat[32]*dv11; 
  out[69] += -0.7071067811865475*Ghat[36]*dv11; 
  out[70] += -0.7071067811865475*Ghat[37]*dv11; 
  out[71] += -0.7071067811865475*Ghat[38]*dv11; 
  out[72] += -0.7071067811865475*Ghat[39]*dv11; 
  out[73] += -1.224744871391589*Ghat[33]*dv11; 
  out[74] += -1.224744871391589*Ghat[34]*dv11; 
  out[75] += -1.224744871391589*Ghat[35]*dv11; 
  out[76] += -0.7071067811865475*Ghat[40]*dv11; 
  out[77] += -0.7071067811865475*Ghat[41]*dv11; 
  out[78] += -0.7071067811865475*Ghat[42]*dv11; 
  out[79] += -1.224744871391589*Ghat[36]*dv11; 
  out[80] += -0.7071067811865475*Ghat[43]*dv11; 
  out[81] += -1.224744871391589*Ghat[37]*dv11; 
  out[82] += -1.224744871391589*Ghat[38]*dv11; 
  out[83] += -1.224744871391589*Ghat[39]*dv11; 
  out[84] += -0.7071067811865475*Ghat[44]*dv11; 
  out[85] += -0.7071067811865475*Ghat[45]*dv11; 
  out[86] += -0.7071067811865475*Ghat[46]*dv11; 
  out[87] += -1.224744871391589*Ghat[40]*dv11; 
  out[88] += -1.224744871391589*Ghat[41]*dv11; 
  out[89] += -1.224744871391589*Ghat[42]*dv11; 
  out[90] += -1.224744871391589*Ghat[43]*dv11; 
  out[91] += -0.7071067811865475*Ghat[47]*dv11; 
  out[92] += -1.224744871391589*Ghat[44]*dv11; 
  out[93] += -1.224744871391589*Ghat[45]*dv11; 
  out[94] += -1.224744871391589*Ghat[46]*dv11; 
  out[95] += -1.224744871391589*Ghat[47]*dv11; 
  out[96] += -1.58113883008419*Ghat[0]*dv11; 
  out[97] += -1.58113883008419*Ghat[1]*dv11; 
  out[98] += -1.58113883008419*Ghat[2]*dv11; 
  out[99] += -1.58113883008419*Ghat[3]*dv11; 
  out[100] += -1.58113883008419*Ghat[4]*dv11; 
  out[101] += -1.58113883008419*Ghat[5]*dv11; 
  out[102] += -1.58113883008419*Ghat[6]*dv11; 
  out[103] += -1.58113883008419*Ghat[7]*dv11; 
  out[104] += -1.58113883008419*Ghat[8]*dv11; 
  out[105] += -1.58113883008419*Ghat[9]*dv11; 
  out[106] += -1.58113883008419*Ghat[10]*dv11; 
  out[107] += -1.58113883008419*Ghat[11]*dv11; 
  out[108] += -1.58113883008419*Ghat[12]*dv11; 
  out[109] += -1.58113883008419*Ghat[13]*dv11; 
  out[110] += -1.58113883008419*Ghat[14]*dv11; 
  out[111] += -1.58113883008419*Ghat[15]*dv11; 
  out[112] += -1.58113883008419*Ghat[16]*dv11; 
  out[113] += -1.58113883008419*Ghat[17]*dv11; 
  out[114] += -1.58113883008419*Ghat[18]*dv11; 
  out[115] += -1.58113883008419*Ghat[19]*dv11; 
  out[116] += -1.58113883008419*Ghat[20]*dv11; 
  out[117] += -1.58113883008419*Ghat[21]*dv11; 
  out[118] += -1.58113883008419*Ghat[22]*dv11; 
  out[119] += -1.58113883008419*Ghat[23]*dv11; 
  out[120] += -1.58113883008419*Ghat[24]*dv11; 
  out[121] += -1.58113883008419*Ghat[25]*dv11; 
  out[122] += -1.58113883008419*Ghat[26]*dv11; 
  out[123] += -1.58113883008419*Ghat[27]*dv11; 
  out[124] += -1.58113883008419*Ghat[28]*dv11; 
  out[125] += -1.58113883008419*Ghat[29]*dv11; 
  out[126] += -1.58113883008419*Ghat[30]*dv11; 
  out[127] += -1.58113883008419*Ghat[31]*dv11; 
  out[128] += -0.7071067811865475*Ghat[48]*dv11; 
  out[129] += -0.7071067811865475*Ghat[49]*dv11; 
  out[130] += -0.7071067811865475*Ghat[50]*dv11; 
  out[131] += -0.7071067811865475*Ghat[51]*dv11; 
  out[132] += -0.7071067811865475*Ghat[52]*dv11; 
  out[133] += -1.224744871391589*Ghat[48]*dv11; 
  out[134] += -0.7071067811865475*Ghat[53]*dv11; 
  out[135] += -0.7071067811865475*Ghat[54]*dv11; 
  out[136] += -0.7071067811865475*Ghat[55]*dv11; 
  out[137] += -0.7071067811865475*Ghat[56]*dv11; 
  out[138] += -0.7071067811865475*Ghat[57]*dv11; 
  out[139] += -0.7071067811865475*Ghat[58]*dv11; 
  out[140] += -1.224744871391589*Ghat[49]*dv11; 
  out[141] += -1.224744871391589*Ghat[50]*dv11; 
  out[142] += -1.224744871391589*Ghat[51]*dv11; 
  out[143] += -1.224744871391589*Ghat[52]*dv11; 
  out[144] += -0.7071067811865475*Ghat[59]*dv11; 
  out[145] += -0.7071067811865475*Ghat[60]*dv11; 
  out[146] += -0.7071067811865475*Ghat[61]*dv11; 
  out[147] += -0.7071067811865475*Ghat[62]*dv11; 
  out[148] += -1.224744871391589*Ghat[53]*dv11; 
  out[149] += -1.224744871391589*Ghat[54]*dv11; 
  out[150] += -1.224744871391589*Ghat[55]*dv11; 
  out[151] += -1.224744871391589*Ghat[56]*dv11; 
  out[152] += -1.224744871391589*Ghat[57]*dv11; 
  out[153] += -1.224744871391589*Ghat[58]*dv11; 
  out[154] += -0.7071067811865475*Ghat[63]*dv11; 
  out[155] += -1.224744871391589*Ghat[59]*dv11; 
  out[156] += -1.224744871391589*Ghat[60]*dv11; 
  out[157] += -1.224744871391589*Ghat[61]*dv11; 
  out[158] += -1.224744871391589*Ghat[62]*dv11; 
  out[159] += -1.224744871391589*Ghat[63]*dv11; 

  } else { 

  if (alpha[27]+alpha[26]-alpha[22]-alpha[21]-alpha[20]-alpha[19]-alpha[18]-alpha[17]-alpha[16]+alpha[14]+alpha[13]+alpha[12]+alpha[11]+alpha[10]+alpha[9]+alpha[8]+alpha[7]+alpha[6]-alpha[5]-alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[0] = ser_6x_p1_surfx5_eval_quad_node_0_r(fEdge); 
  } else { 
    fUpwindQuad[0] = ser_6x_p1_surfx5_eval_quad_node_0_l(fSkin); 
  } 
  if ((-alpha[27])+alpha[26]+alpha[22]+alpha[21]+alpha[20]-alpha[19]-alpha[18]-alpha[17]-alpha[16]-alpha[14]-alpha[13]-alpha[12]+alpha[11]+alpha[10]+alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5]-alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[1] = ser_6x_p1_surfx5_eval_quad_node_1_r(fEdge); 
  } else { 
    fUpwindQuad[1] = ser_6x_p1_surfx5_eval_quad_node_1_l(fSkin); 
  } 
  if (alpha[27]-alpha[26]-alpha[22]-alpha[21]-alpha[20]+alpha[19]+alpha[18]+alpha[17]-alpha[16]+alpha[14]+alpha[13]+alpha[12]-alpha[11]-alpha[10]-alpha[9]+alpha[8]+alpha[7]+alpha[6]-alpha[5]+alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[2] = ser_6x_p1_surfx5_eval_quad_node_2_r(fEdge); 
  } else { 
    fUpwindQuad[2] = ser_6x_p1_surfx5_eval_quad_node_2_l(fSkin); 
  } 
  if ((-alpha[27])-alpha[26]+alpha[22]+alpha[21]+alpha[20]+alpha[19]+alpha[18]+alpha[17]-alpha[16]-alpha[14]-alpha[13]-alpha[12]-alpha[11]-alpha[10]-alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5]+alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[3] = ser_6x_p1_surfx5_eval_quad_node_3_r(fEdge); 
  } else { 
    fUpwindQuad[3] = ser_6x_p1_surfx5_eval_quad_node_3_l(fSkin); 
  } 
  if ((-alpha[27])-alpha[26]+alpha[22]+alpha[21]-alpha[20]+alpha[19]+alpha[18]-alpha[17]+alpha[16]-alpha[14]+alpha[13]+alpha[12]-alpha[11]+alpha[10]+alpha[9]-alpha[8]-alpha[7]+alpha[6]-alpha[5]-alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[4] = ser_6x_p1_surfx5_eval_quad_node_4_r(fEdge); 
  } else { 
    fUpwindQuad[4] = ser_6x_p1_surfx5_eval_quad_node_4_l(fSkin); 
  } 
  if (alpha[27]-alpha[26]-alpha[22]-alpha[21]+alpha[20]+alpha[19]+alpha[18]-alpha[17]+alpha[16]+alpha[14]-alpha[13]-alpha[12]-alpha[11]+alpha[10]+alpha[9]-alpha[8]-alpha[7]+alpha[6]+alpha[5]-alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[5] = ser_6x_p1_surfx5_eval_quad_node_5_r(fEdge); 
  } else { 
    fUpwindQuad[5] = ser_6x_p1_surfx5_eval_quad_node_5_l(fSkin); 
  } 
  if ((-alpha[27])+alpha[26]+alpha[22]+alpha[21]-alpha[20]-alpha[19]-alpha[18]+alpha[17]+alpha[16]-alpha[14]+alpha[13]+alpha[12]+alpha[11]-alpha[10]-alpha[9]-alpha[8]-alpha[7]+alpha[6]-alpha[5]+alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[6] = ser_6x_p1_surfx5_eval_quad_node_6_r(fEdge); 
  } else { 
    fUpwindQuad[6] = ser_6x_p1_surfx5_eval_quad_node_6_l(fSkin); 
  } 
  if (alpha[27]+alpha[26]-alpha[22]-alpha[21]+alpha[20]-alpha[19]-alpha[18]+alpha[17]+alpha[16]+alpha[14]-alpha[13]-alpha[12]+alpha[11]-alpha[10]-alpha[9]-alpha[8]-alpha[7]+alpha[6]+alpha[5]+alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[7] = ser_6x_p1_surfx5_eval_quad_node_7_r(fEdge); 
  } else { 
    fUpwindQuad[7] = ser_6x_p1_surfx5_eval_quad_node_7_l(fSkin); 
  } 
  if ((-alpha[27])-alpha[26]+alpha[22]-alpha[21]+alpha[20]+alpha[19]-alpha[18]+alpha[17]+alpha[16]+alpha[14]-alpha[13]+alpha[12]+alpha[11]-alpha[10]+alpha[9]-alpha[8]+alpha[7]-alpha[6]-alpha[5]-alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[8] = ser_6x_p1_surfx5_eval_quad_node_8_r(fEdge); 
  } else { 
    fUpwindQuad[8] = ser_6x_p1_surfx5_eval_quad_node_8_l(fSkin); 
  } 
  if (alpha[27]-alpha[26]-alpha[22]+alpha[21]-alpha[20]+alpha[19]-alpha[18]+alpha[17]+alpha[16]-alpha[14]+alpha[13]-alpha[12]+alpha[11]-alpha[10]+alpha[9]-alpha[8]+alpha[7]-alpha[6]+alpha[5]-alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[9] = ser_6x_p1_surfx5_eval_quad_node_9_r(fEdge); 
  } else { 
    fUpwindQuad[9] = ser_6x_p1_surfx5_eval_quad_node_9_l(fSkin); 
  } 
  if ((-alpha[27])+alpha[26]+alpha[22]-alpha[21]+alpha[20]-alpha[19]+alpha[18]-alpha[17]+alpha[16]+alpha[14]-alpha[13]+alpha[12]-alpha[11]+alpha[10]-alpha[9]-alpha[8]+alpha[7]-alpha[6]-alpha[5]+alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[10] = ser_6x_p1_surfx5_eval_quad_node_10_r(fEdge); 
  } else { 
    fUpwindQuad[10] = ser_6x_p1_surfx5_eval_quad_node_10_l(fSkin); 
  } 
  if (alpha[27]+alpha[26]-alpha[22]+alpha[21]-alpha[20]-alpha[19]+alpha[18]-alpha[17]+alpha[16]-alpha[14]+alpha[13]-alpha[12]-alpha[11]+alpha[10]-alpha[9]-alpha[8]+alpha[7]-alpha[6]+alpha[5]+alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[11] = ser_6x_p1_surfx5_eval_quad_node_11_r(fEdge); 
  } else { 
    fUpwindQuad[11] = ser_6x_p1_surfx5_eval_quad_node_11_l(fSkin); 
  } 
  if (alpha[27]+alpha[26]-alpha[22]+alpha[21]+alpha[20]-alpha[19]+alpha[18]+alpha[17]-alpha[16]-alpha[14]-alpha[13]+alpha[12]-alpha[11]-alpha[10]+alpha[9]+alpha[8]-alpha[7]-alpha[6]-alpha[5]-alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[12] = ser_6x_p1_surfx5_eval_quad_node_12_r(fEdge); 
  } else { 
    fUpwindQuad[12] = ser_6x_p1_surfx5_eval_quad_node_12_l(fSkin); 
  } 
  if ((-alpha[27])+alpha[26]+alpha[22]-alpha[21]-alpha[20]-alpha[19]+alpha[18]+alpha[17]-alpha[16]+alpha[14]+alpha[13]-alpha[12]-alpha[11]-alpha[10]+alpha[9]+alpha[8]-alpha[7]-alpha[6]+alpha[5]-alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[13] = ser_6x_p1_surfx5_eval_quad_node_13_r(fEdge); 
  } else { 
    fUpwindQuad[13] = ser_6x_p1_surfx5_eval_quad_node_13_l(fSkin); 
  } 
  if (alpha[27]-alpha[26]-alpha[22]+alpha[21]+alpha[20]+alpha[19]-alpha[18]-alpha[17]-alpha[16]-alpha[14]-alpha[13]+alpha[12]+alpha[11]+alpha[10]-alpha[9]+alpha[8]-alpha[7]-alpha[6]-alpha[5]+alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[14] = ser_6x_p1_surfx5_eval_quad_node_14_r(fEdge); 
  } else { 
    fUpwindQuad[14] = ser_6x_p1_surfx5_eval_quad_node_14_l(fSkin); 
  } 
  if ((-alpha[27])-alpha[26]+alpha[22]-alpha[21]-alpha[20]+alpha[19]-alpha[18]-alpha[17]-alpha[16]+alpha[14]+alpha[13]-alpha[12]+alpha[11]+alpha[10]-alpha[9]+alpha[8]-alpha[7]-alpha[6]+alpha[5]+alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[15] = ser_6x_p1_surfx5_eval_quad_node_15_r(fEdge); 
  } else { 
    fUpwindQuad[15] = ser_6x_p1_surfx5_eval_quad_node_15_l(fSkin); 
  } 
  if ((-alpha[27])-alpha[26]-alpha[22]+alpha[21]+alpha[20]-alpha[19]+alpha[18]+alpha[17]+alpha[16]+alpha[14]+alpha[13]-alpha[12]+alpha[11]+alpha[10]-alpha[9]+alpha[8]-alpha[7]-alpha[6]-alpha[5]-alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[16] = ser_6x_p1_surfx5_eval_quad_node_16_r(fEdge); 
  } else { 
    fUpwindQuad[16] = ser_6x_p1_surfx5_eval_quad_node_16_l(fSkin); 
  } 
  if (alpha[27]-alpha[26]+alpha[22]-alpha[21]-alpha[20]-alpha[19]+alpha[18]+alpha[17]+alpha[16]-alpha[14]-alpha[13]+alpha[12]+alpha[11]+alpha[10]-alpha[9]+alpha[8]-alpha[7]-alpha[6]+alpha[5]-alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[17] = ser_6x_p1_surfx5_eval_quad_node_17_r(fEdge); 
  } else { 
    fUpwindQuad[17] = ser_6x_p1_surfx5_eval_quad_node_17_l(fSkin); 
  } 
  if ((-alpha[27])+alpha[26]-alpha[22]+alpha[21]+alpha[20]+alpha[19]-alpha[18]-alpha[17]+alpha[16]+alpha[14]+alpha[13]-alpha[12]-alpha[11]-alpha[10]+alpha[9]+alpha[8]-alpha[7]-alpha[6]-alpha[5]+alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[18] = ser_6x_p1_surfx5_eval_quad_node_18_r(fEdge); 
  } else { 
    fUpwindQuad[18] = ser_6x_p1_surfx5_eval_quad_node_18_l(fSkin); 
  } 
  if (alpha[27]+alpha[26]+alpha[22]-alpha[21]-alpha[20]+alpha[19]-alpha[18]-alpha[17]+alpha[16]-alpha[14]-alpha[13]+alpha[12]-alpha[11]-alpha[10]+alpha[9]+alpha[8]-alpha[7]-alpha[6]+alpha[5]+alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[19] = ser_6x_p1_surfx5_eval_quad_node_19_r(fEdge); 
  } else { 
    fUpwindQuad[19] = ser_6x_p1_surfx5_eval_quad_node_19_l(fSkin); 
  } 
  if (alpha[27]+alpha[26]+alpha[22]-alpha[21]+alpha[20]+alpha[19]-alpha[18]+alpha[17]-alpha[16]-alpha[14]+alpha[13]-alpha[12]-alpha[11]+alpha[10]-alpha[9]-alpha[8]+alpha[7]-alpha[6]-alpha[5]-alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[20] = ser_6x_p1_surfx5_eval_quad_node_20_r(fEdge); 
  } else { 
    fUpwindQuad[20] = ser_6x_p1_surfx5_eval_quad_node_20_l(fSkin); 
  } 
  if ((-alpha[27])+alpha[26]-alpha[22]+alpha[21]-alpha[20]+alpha[19]-alpha[18]+alpha[17]-alpha[16]+alpha[14]-alpha[13]+alpha[12]-alpha[11]+alpha[10]-alpha[9]-alpha[8]+alpha[7]-alpha[6]+alpha[5]-alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[21] = ser_6x_p1_surfx5_eval_quad_node_21_r(fEdge); 
  } else { 
    fUpwindQuad[21] = ser_6x_p1_surfx5_eval_quad_node_21_l(fSkin); 
  } 
  if (alpha[27]-alpha[26]+alpha[22]-alpha[21]+alpha[20]-alpha[19]+alpha[18]-alpha[17]-alpha[16]-alpha[14]+alpha[13]-alpha[12]+alpha[11]-alpha[10]+alpha[9]-alpha[8]+alpha[7]-alpha[6]-alpha[5]+alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[22] = ser_6x_p1_surfx5_eval_quad_node_22_r(fEdge); 
  } else { 
    fUpwindQuad[22] = ser_6x_p1_surfx5_eval_quad_node_22_l(fSkin); 
  } 
  if ((-alpha[27])-alpha[26]-alpha[22]+alpha[21]-alpha[20]-alpha[19]+alpha[18]-alpha[17]-alpha[16]+alpha[14]-alpha[13]+alpha[12]+alpha[11]-alpha[10]+alpha[9]-alpha[8]+alpha[7]-alpha[6]+alpha[5]+alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[23] = ser_6x_p1_surfx5_eval_quad_node_23_r(fEdge); 
  } else { 
    fUpwindQuad[23] = ser_6x_p1_surfx5_eval_quad_node_23_l(fSkin); 
  } 
  if (alpha[27]+alpha[26]+alpha[22]+alpha[21]-alpha[20]+alpha[19]+alpha[18]-alpha[17]-alpha[16]+alpha[14]-alpha[13]-alpha[12]+alpha[11]-alpha[10]-alpha[9]-alpha[8]-alpha[7]+alpha[6]-alpha[5]-alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[24] = ser_6x_p1_surfx5_eval_quad_node_24_r(fEdge); 
  } else { 
    fUpwindQuad[24] = ser_6x_p1_surfx5_eval_quad_node_24_l(fSkin); 
  } 
  if ((-alpha[27])+alpha[26]-alpha[22]-alpha[21]+alpha[20]+alpha[19]+alpha[18]-alpha[17]-alpha[16]-alpha[14]+alpha[13]+alpha[12]+alpha[11]-alpha[10]-alpha[9]-alpha[8]-alpha[7]+alpha[6]+alpha[5]-alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[25] = ser_6x_p1_surfx5_eval_quad_node_25_r(fEdge); 
  } else { 
    fUpwindQuad[25] = ser_6x_p1_surfx5_eval_quad_node_25_l(fSkin); 
  } 
  if (alpha[27]-alpha[26]+alpha[22]+alpha[21]-alpha[20]-alpha[19]-alpha[18]+alpha[17]-alpha[16]+alpha[14]-alpha[13]-alpha[12]-alpha[11]+alpha[10]+alpha[9]-alpha[8]-alpha[7]+alpha[6]-alpha[5]+alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[26] = ser_6x_p1_surfx5_eval_quad_node_26_r(fEdge); 
  } else { 
    fUpwindQuad[26] = ser_6x_p1_surfx5_eval_quad_node_26_l(fSkin); 
  } 
  if ((-alpha[27])-alpha[26]-alpha[22]-alpha[21]+alpha[20]-alpha[19]-alpha[18]+alpha[17]-alpha[16]-alpha[14]+alpha[13]+alpha[12]-alpha[11]+alpha[10]+alpha[9]-alpha[8]-alpha[7]+alpha[6]+alpha[5]+alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[27] = ser_6x_p1_surfx5_eval_quad_node_27_r(fEdge); 
  } else { 
    fUpwindQuad[27] = ser_6x_p1_surfx5_eval_quad_node_27_l(fSkin); 
  } 
  if ((-alpha[27])-alpha[26]-alpha[22]-alpha[21]-alpha[20]-alpha[19]-alpha[18]-alpha[17]+alpha[16]-alpha[14]-alpha[13]-alpha[12]-alpha[11]-alpha[10]-alpha[9]+alpha[8]+alpha[7]+alpha[6]-alpha[5]-alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[28] = ser_6x_p1_surfx5_eval_quad_node_28_r(fEdge); 
  } else { 
    fUpwindQuad[28] = ser_6x_p1_surfx5_eval_quad_node_28_l(fSkin); 
  } 
  if (alpha[27]-alpha[26]+alpha[22]+alpha[21]+alpha[20]-alpha[19]-alpha[18]-alpha[17]+alpha[16]+alpha[14]+alpha[13]+alpha[12]-alpha[11]-alpha[10]-alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5]-alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[29] = ser_6x_p1_surfx5_eval_quad_node_29_r(fEdge); 
  } else { 
    fUpwindQuad[29] = ser_6x_p1_surfx5_eval_quad_node_29_l(fSkin); 
  } 
  if ((-alpha[27])+alpha[26]-alpha[22]-alpha[21]-alpha[20]+alpha[19]+alpha[18]+alpha[17]+alpha[16]-alpha[14]-alpha[13]-alpha[12]+alpha[11]+alpha[10]+alpha[9]+alpha[8]+alpha[7]+alpha[6]-alpha[5]+alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[30] = ser_6x_p1_surfx5_eval_quad_node_30_r(fEdge); 
  } else { 
    fUpwindQuad[30] = ser_6x_p1_surfx5_eval_quad_node_30_l(fSkin); 
  } 
  if (alpha[27]+alpha[26]+alpha[22]+alpha[21]+alpha[20]+alpha[19]+alpha[18]+alpha[17]+alpha[16]+alpha[14]+alpha[13]+alpha[12]+alpha[11]+alpha[10]+alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5]+alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[31] = ser_6x_p1_surfx5_eval_quad_node_31_r(fEdge); 
  } else { 
    fUpwindQuad[31] = ser_6x_p1_surfx5_eval_quad_node_31_l(fSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_6x_p1_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.1767766952966368*(alpha[27]*fUpwind[27]+alpha[26]*fUpwind[26]+alpha[22]*fUpwind[22]+alpha[21]*fUpwind[21]+alpha[20]*fUpwind[20]+alpha[19]*fUpwind[19]+alpha[18]*fUpwind[18]+alpha[17]*fUpwind[17]+alpha[16]*fUpwind[16]+alpha[14]*fUpwind[14]+alpha[13]*fUpwind[13]+alpha[12]*fUpwind[12]+alpha[11]*fUpwind[11]+alpha[10]*fUpwind[10]+alpha[9]*fUpwind[9]+alpha[8]*fUpwind[8]+alpha[7]*fUpwind[7]+alpha[6]*fUpwind[6]+alpha[5]*fUpwind[5]+alpha[4]*fUpwind[4]+alpha[3]*fUpwind[3]+alpha[2]*fUpwind[2]+alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] = 0.1767766952966368*(alpha[22]*fUpwind[27]+fUpwind[22]*alpha[27]+alpha[19]*fUpwind[26]+fUpwind[19]*alpha[26]+alpha[14]*fUpwind[21]+fUpwind[14]*alpha[21]+alpha[13]*fUpwind[20]+fUpwind[13]*alpha[20]+alpha[11]*fUpwind[18]+fUpwind[11]*alpha[18]+alpha[10]*fUpwind[17]+fUpwind[10]*alpha[17]+alpha[8]*fUpwind[16]+fUpwind[8]*alpha[16]+alpha[5]*fUpwind[12]+fUpwind[5]*alpha[12]+alpha[4]*fUpwind[9]+fUpwind[4]*alpha[9]+alpha[3]*fUpwind[7]+fUpwind[3]*alpha[7]+alpha[2]*fUpwind[6]+fUpwind[2]*alpha[6]+alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]); 
  Ghat[2] = 0.1767766952966368*(alpha[21]*fUpwind[27]+fUpwind[21]*alpha[27]+alpha[18]*fUpwind[26]+fUpwind[18]*alpha[26]+alpha[14]*fUpwind[22]+fUpwind[14]*alpha[22]+alpha[12]*fUpwind[20]+fUpwind[12]*alpha[20]+alpha[11]*fUpwind[19]+fUpwind[11]*alpha[19]+alpha[9]*fUpwind[17]+fUpwind[9]*alpha[17]+alpha[7]*fUpwind[16]+fUpwind[7]*alpha[16]+alpha[5]*fUpwind[13]+fUpwind[5]*alpha[13]+alpha[4]*fUpwind[10]+fUpwind[4]*alpha[10]+alpha[3]*fUpwind[8]+fUpwind[3]*alpha[8]+alpha[1]*fUpwind[6]+fUpwind[1]*alpha[6]+alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2]); 
  Ghat[3] = 0.1767766952966368*(alpha[20]*fUpwind[27]+fUpwind[20]*alpha[27]+alpha[17]*fUpwind[26]+fUpwind[17]*alpha[26]+alpha[13]*fUpwind[22]+fUpwind[13]*alpha[22]+alpha[12]*fUpwind[21]+fUpwind[12]*alpha[21]+alpha[10]*fUpwind[19]+fUpwind[10]*alpha[19]+alpha[9]*fUpwind[18]+fUpwind[9]*alpha[18]+alpha[6]*fUpwind[16]+fUpwind[6]*alpha[16]+alpha[5]*fUpwind[14]+fUpwind[5]*alpha[14]+alpha[4]*fUpwind[11]+fUpwind[4]*alpha[11]+alpha[2]*fUpwind[8]+fUpwind[2]*alpha[8]+alpha[1]*fUpwind[7]+fUpwind[1]*alpha[7]+alpha[0]*fUpwind[3]+fUpwind[0]*alpha[3]); 
  Ghat[4] = 0.01178511301977579*(13.41640786499874*alpha[26]*fUpwind[43]+13.41640786499874*(alpha[19]*fUpwind[39]+alpha[18]*fUpwind[38]+alpha[17]*fUpwind[37])+13.41640786499874*(alpha[11]*fUpwind[35]+alpha[10]*fUpwind[34]+alpha[9]*fUpwind[33])+13.41640786499874*alpha[4]*fUpwind[32]+15.0*(alpha[27]*fUpwind[31]+alpha[22]*fUpwind[30]+alpha[21]*fUpwind[29]+alpha[20]*fUpwind[28]+alpha[16]*fUpwind[26]+fUpwind[16]*alpha[26]+alpha[14]*fUpwind[25]+alpha[13]*fUpwind[24]+alpha[12]*fUpwind[23]+alpha[8]*fUpwind[19]+fUpwind[8]*alpha[19]+alpha[7]*fUpwind[18]+fUpwind[7]*alpha[18]+alpha[6]*fUpwind[17]+fUpwind[6]*alpha[17]+alpha[5]*fUpwind[15]+alpha[3]*fUpwind[11]+fUpwind[3]*alpha[11]+alpha[2]*fUpwind[10]+fUpwind[2]*alpha[10]+alpha[1]*fUpwind[9]+fUpwind[1]*alpha[9]+alpha[0]*fUpwind[4]+fUpwind[0]*alpha[4])); 
  Ghat[5] = 0.01178511301977579*(13.41640786499874*alpha[27]*fUpwind[59]+13.41640786499874*(alpha[22]*fUpwind[55]+alpha[21]*fUpwind[54]+alpha[20]*fUpwind[53])+13.41640786499874*(alpha[14]*fUpwind[51]+alpha[13]*fUpwind[50]+alpha[12]*fUpwind[49])+13.41640786499874*alpha[5]*fUpwind[48]+15.0*(alpha[26]*fUpwind[31]+alpha[19]*fUpwind[30]+alpha[18]*fUpwind[29]+alpha[17]*fUpwind[28]+alpha[16]*fUpwind[27]+fUpwind[16]*alpha[27]+alpha[11]*fUpwind[25]+alpha[10]*fUpwind[24]+alpha[9]*fUpwind[23]+alpha[8]*fUpwind[22]+fUpwind[8]*alpha[22]+alpha[7]*fUpwind[21]+fUpwind[7]*alpha[21]+alpha[6]*fUpwind[20]+fUpwind[6]*alpha[20]+alpha[4]*fUpwind[15]+alpha[3]*fUpwind[14]+fUpwind[3]*alpha[14]+alpha[2]*fUpwind[13]+fUpwind[2]*alpha[13]+alpha[1]*fUpwind[12]+fUpwind[1]*alpha[12]+alpha[0]*fUpwind[5]+fUpwind[0]*alpha[5])); 
  Ghat[6] = 0.1767766952966368*(alpha[14]*fUpwind[27]+fUpwind[14]*alpha[27]+alpha[11]*fUpwind[26]+fUpwind[11]*alpha[26]+alpha[21]*fUpwind[22]+fUpwind[21]*alpha[22]+alpha[5]*fUpwind[20]+fUpwind[5]*alpha[20]+alpha[18]*fUpwind[19]+fUpwind[18]*alpha[19]+alpha[4]*fUpwind[17]+fUpwind[4]*alpha[17]+alpha[3]*fUpwind[16]+fUpwind[3]*alpha[16]+alpha[12]*fUpwind[13]+fUpwind[12]*alpha[13]+alpha[9]*fUpwind[10]+fUpwind[9]*alpha[10]+alpha[7]*fUpwind[8]+fUpwind[7]*alpha[8]+alpha[0]*fUpwind[6]+fUpwind[0]*alpha[6]+alpha[1]*fUpwind[2]+fUpwind[1]*alpha[2]); 
  Ghat[7] = 0.1767766952966368*(alpha[13]*fUpwind[27]+fUpwind[13]*alpha[27]+alpha[10]*fUpwind[26]+fUpwind[10]*alpha[26]+alpha[20]*fUpwind[22]+fUpwind[20]*alpha[22]+alpha[5]*fUpwind[21]+fUpwind[5]*alpha[21]+alpha[17]*fUpwind[19]+fUpwind[17]*alpha[19]+alpha[4]*fUpwind[18]+fUpwind[4]*alpha[18]+alpha[2]*fUpwind[16]+fUpwind[2]*alpha[16]+alpha[12]*fUpwind[14]+fUpwind[12]*alpha[14]+alpha[9]*fUpwind[11]+fUpwind[9]*alpha[11]+alpha[6]*fUpwind[8]+fUpwind[6]*alpha[8]+alpha[0]*fUpwind[7]+fUpwind[0]*alpha[7]+alpha[1]*fUpwind[3]+fUpwind[1]*alpha[3]); 
  Ghat[8] = 0.1767766952966368*(alpha[12]*fUpwind[27]+fUpwind[12]*alpha[27]+alpha[9]*fUpwind[26]+fUpwind[9]*alpha[26]+alpha[5]*fUpwind[22]+fUpwind[5]*alpha[22]+alpha[20]*fUpwind[21]+fUpwind[20]*alpha[21]+alpha[4]*fUpwind[19]+fUpwind[4]*alpha[19]+alpha[17]*fUpwind[18]+fUpwind[17]*alpha[18]+alpha[1]*fUpwind[16]+fUpwind[1]*alpha[16]+alpha[13]*fUpwind[14]+fUpwind[13]*alpha[14]+alpha[10]*fUpwind[11]+fUpwind[10]*alpha[11]+alpha[0]*fUpwind[8]+fUpwind[0]*alpha[8]+alpha[6]*fUpwind[7]+fUpwind[6]*alpha[7]+alpha[2]*fUpwind[3]+fUpwind[2]*alpha[3]); 
  Ghat[9] = 0.01178511301977579*(13.41640786499874*alpha[19]*fUpwind[43]+13.41640786499874*(alpha[26]*fUpwind[39]+alpha[11]*fUpwind[38]+alpha[10]*fUpwind[37])+13.41640786499874*(alpha[18]*fUpwind[35]+alpha[17]*fUpwind[34]+alpha[4]*fUpwind[33])+13.41640786499874*alpha[9]*fUpwind[32]+15.0*(alpha[22]*fUpwind[31]+alpha[27]*fUpwind[30]+alpha[14]*fUpwind[29]+alpha[13]*fUpwind[28]+alpha[8]*fUpwind[26]+fUpwind[8]*alpha[26]+alpha[21]*fUpwind[25]+alpha[20]*fUpwind[24]+alpha[5]*fUpwind[23]+alpha[16]*fUpwind[19]+fUpwind[16]*alpha[19]+alpha[3]*fUpwind[18]+fUpwind[3]*alpha[18]+alpha[2]*fUpwind[17]+fUpwind[2]*alpha[17]+alpha[12]*fUpwind[15]+alpha[7]*fUpwind[11]+fUpwind[7]*alpha[11]+alpha[6]*fUpwind[10]+fUpwind[6]*alpha[10]+alpha[0]*fUpwind[9]+fUpwind[0]*alpha[9]+alpha[1]*fUpwind[4]+fUpwind[1]*alpha[4])); 
  Ghat[10] = 0.01178511301977579*(13.41640786499874*alpha[18]*fUpwind[43]+13.41640786499874*(alpha[11]*fUpwind[39]+alpha[26]*fUpwind[38]+alpha[9]*fUpwind[37])+13.41640786499874*(alpha[19]*fUpwind[35]+alpha[4]*fUpwind[34]+alpha[17]*fUpwind[33])+13.41640786499874*alpha[10]*fUpwind[32]+15.0*(alpha[21]*fUpwind[31]+alpha[14]*fUpwind[30]+alpha[27]*fUpwind[29]+alpha[12]*fUpwind[28]+alpha[7]*fUpwind[26]+fUpwind[7]*alpha[26]+alpha[22]*fUpwind[25]+alpha[5]*fUpwind[24]+alpha[20]*fUpwind[23]+alpha[3]*fUpwind[19]+fUpwind[3]*alpha[19]+alpha[16]*fUpwind[18]+fUpwind[16]*alpha[18]+alpha[1]*fUpwind[17]+fUpwind[1]*alpha[17]+alpha[13]*fUpwind[15]+alpha[8]*fUpwind[11]+fUpwind[8]*alpha[11]+alpha[0]*fUpwind[10]+fUpwind[0]*alpha[10]+alpha[6]*fUpwind[9]+fUpwind[6]*alpha[9]+alpha[2]*fUpwind[4]+fUpwind[2]*alpha[4])); 
  Ghat[11] = 0.01178511301977579*(13.41640786499874*alpha[17]*fUpwind[43]+13.41640786499874*(alpha[10]*fUpwind[39]+alpha[9]*fUpwind[38]+alpha[26]*fUpwind[37])+13.41640786499874*(alpha[4]*fUpwind[35]+alpha[19]*fUpwind[34]+alpha[18]*fUpwind[33])+13.41640786499874*alpha[11]*fUpwind[32]+15.0*(alpha[20]*fUpwind[31]+alpha[13]*fUpwind[30]+alpha[12]*fUpwind[29]+alpha[27]*fUpwind[28]+alpha[6]*fUpwind[26]+fUpwind[6]*alpha[26]+alpha[5]*fUpwind[25]+alpha[22]*fUpwind[24]+alpha[21]*fUpwind[23]+alpha[2]*fUpwind[19]+fUpwind[2]*alpha[19]+alpha[1]*fUpwind[18]+fUpwind[1]*alpha[18]+alpha[16]*fUpwind[17]+fUpwind[16]*alpha[17]+alpha[14]*fUpwind[15]+alpha[0]*fUpwind[11]+fUpwind[0]*alpha[11]+alpha[8]*fUpwind[10]+fUpwind[8]*alpha[10]+alpha[7]*fUpwind[9]+fUpwind[7]*alpha[9]+alpha[3]*fUpwind[4]+fUpwind[3]*alpha[4])); 
  Ghat[12] = 0.01178511301977579*(13.41640786499874*alpha[22]*fUpwind[59]+13.41640786499874*(alpha[27]*fUpwind[55]+alpha[14]*fUpwind[54]+alpha[13]*fUpwind[53])+13.41640786499874*(alpha[21]*fUpwind[51]+alpha[20]*fUpwind[50]+alpha[5]*fUpwind[49])+13.41640786499874*alpha[12]*fUpwind[48]+15.0*(alpha[19]*fUpwind[31]+alpha[26]*fUpwind[30]+alpha[11]*fUpwind[29]+alpha[10]*fUpwind[28]+alpha[8]*fUpwind[27]+fUpwind[8]*alpha[27]+alpha[18]*fUpwind[25]+alpha[17]*fUpwind[24]+alpha[4]*fUpwind[23]+alpha[16]*fUpwind[22]+fUpwind[16]*alpha[22]+alpha[3]*fUpwind[21]+fUpwind[3]*alpha[21]+alpha[2]*fUpwind[20]+fUpwind[2]*alpha[20]+alpha[9]*fUpwind[15]+alpha[7]*fUpwind[14]+fUpwind[7]*alpha[14]+alpha[6]*fUpwind[13]+fUpwind[6]*alpha[13]+alpha[0]*fUpwind[12]+fUpwind[0]*alpha[12]+alpha[1]*fUpwind[5]+fUpwind[1]*alpha[5])); 
  Ghat[13] = 0.01178511301977579*(13.41640786499874*alpha[21]*fUpwind[59]+13.41640786499874*(alpha[14]*fUpwind[55]+alpha[27]*fUpwind[54]+alpha[12]*fUpwind[53])+13.41640786499874*(alpha[22]*fUpwind[51]+alpha[5]*fUpwind[50]+alpha[20]*fUpwind[49])+13.41640786499874*alpha[13]*fUpwind[48]+15.0*(alpha[18]*fUpwind[31]+alpha[11]*fUpwind[30]+alpha[26]*fUpwind[29]+alpha[9]*fUpwind[28]+alpha[7]*fUpwind[27]+fUpwind[7]*alpha[27]+alpha[19]*fUpwind[25]+alpha[4]*fUpwind[24]+alpha[17]*fUpwind[23]+alpha[3]*fUpwind[22]+fUpwind[3]*alpha[22]+alpha[16]*fUpwind[21]+fUpwind[16]*alpha[21]+alpha[1]*fUpwind[20]+fUpwind[1]*alpha[20]+alpha[10]*fUpwind[15]+alpha[8]*fUpwind[14]+fUpwind[8]*alpha[14]+alpha[0]*fUpwind[13]+fUpwind[0]*alpha[13]+alpha[6]*fUpwind[12]+fUpwind[6]*alpha[12]+alpha[2]*fUpwind[5]+fUpwind[2]*alpha[5])); 
  Ghat[14] = 0.01178511301977579*(13.41640786499874*alpha[20]*fUpwind[59]+13.41640786499874*(alpha[13]*fUpwind[55]+alpha[12]*fUpwind[54]+alpha[27]*fUpwind[53])+13.41640786499874*(alpha[5]*fUpwind[51]+alpha[22]*fUpwind[50]+alpha[21]*fUpwind[49])+13.41640786499874*alpha[14]*fUpwind[48]+15.0*(alpha[17]*fUpwind[31]+alpha[10]*fUpwind[30]+alpha[9]*fUpwind[29]+alpha[26]*fUpwind[28]+alpha[6]*fUpwind[27]+fUpwind[6]*alpha[27]+alpha[4]*fUpwind[25]+alpha[19]*fUpwind[24]+alpha[18]*fUpwind[23]+alpha[2]*fUpwind[22]+fUpwind[2]*alpha[22]+alpha[1]*fUpwind[21]+fUpwind[1]*alpha[21]+alpha[16]*fUpwind[20]+fUpwind[16]*alpha[20]+alpha[11]*fUpwind[15]+alpha[0]*fUpwind[14]+fUpwind[0]*alpha[14]+alpha[8]*fUpwind[13]+fUpwind[8]*alpha[13]+alpha[7]*fUpwind[12]+fUpwind[7]*alpha[12]+alpha[3]*fUpwind[5]+fUpwind[3]*alpha[5])); 
  Ghat[15] = 0.01178511301977579*(13.41640786499874*alpha[27]*fUpwind[63]+13.41640786499874*(alpha[22]*fUpwind[62]+alpha[21]*fUpwind[61]+alpha[20]*fUpwind[60])+13.41640786499874*(alpha[14]*fUpwind[58]+alpha[13]*fUpwind[57]+alpha[12]*fUpwind[56])+13.41640786499874*alpha[5]*fUpwind[52]+13.41640786499874*alpha[26]*fUpwind[47]+13.41640786499874*(alpha[19]*fUpwind[46]+alpha[18]*fUpwind[45]+alpha[17]*fUpwind[44])+13.41640786499874*(alpha[11]*fUpwind[42]+alpha[10]*fUpwind[41]+alpha[9]*fUpwind[40])+13.41640786499874*alpha[4]*fUpwind[36]+15.0*(alpha[16]*fUpwind[31]+alpha[8]*fUpwind[30]+alpha[7]*fUpwind[29]+alpha[6]*fUpwind[28]+alpha[26]*fUpwind[27]+fUpwind[26]*alpha[27]+alpha[3]*fUpwind[25]+alpha[2]*fUpwind[24]+alpha[1]*fUpwind[23]+alpha[19]*fUpwind[22]+fUpwind[19]*alpha[22]+alpha[18]*fUpwind[21]+fUpwind[18]*alpha[21]+alpha[17]*fUpwind[20]+fUpwind[17]*alpha[20]+alpha[0]*fUpwind[15]+alpha[11]*fUpwind[14]+fUpwind[11]*alpha[14]+alpha[10]*fUpwind[13]+fUpwind[10]*alpha[13]+alpha[9]*fUpwind[12]+fUpwind[9]*alpha[12]+alpha[4]*fUpwind[5]+fUpwind[4]*alpha[5])); 
  Ghat[16] = 0.1767766952966368*(alpha[5]*fUpwind[27]+fUpwind[5]*alpha[27]+alpha[4]*fUpwind[26]+fUpwind[4]*alpha[26]+alpha[12]*fUpwind[22]+fUpwind[12]*alpha[22]+alpha[13]*fUpwind[21]+fUpwind[13]*alpha[21]+alpha[14]*fUpwind[20]+fUpwind[14]*alpha[20]+alpha[9]*fUpwind[19]+fUpwind[9]*alpha[19]+alpha[10]*fUpwind[18]+fUpwind[10]*alpha[18]+alpha[11]*fUpwind[17]+fUpwind[11]*alpha[17]+alpha[0]*fUpwind[16]+fUpwind[0]*alpha[16]+alpha[1]*fUpwind[8]+fUpwind[1]*alpha[8]+alpha[2]*fUpwind[7]+fUpwind[2]*alpha[7]+alpha[3]*fUpwind[6]+fUpwind[3]*alpha[6]); 
  Ghat[17] = 0.01178511301977579*(13.41640786499874*alpha[11]*fUpwind[43]+13.41640786499874*(alpha[18]*fUpwind[39]+alpha[19]*fUpwind[38]+alpha[4]*fUpwind[37])+13.41640786499874*(alpha[26]*fUpwind[35]+alpha[9]*fUpwind[34]+alpha[10]*fUpwind[33])+13.41640786499874*alpha[17]*fUpwind[32]+15.0*(alpha[14]*fUpwind[31]+alpha[21]*fUpwind[30]+alpha[22]*fUpwind[29]+alpha[5]*fUpwind[28]+fUpwind[25]*alpha[27]+alpha[3]*fUpwind[26]+fUpwind[3]*alpha[26]+alpha[12]*fUpwind[24]+alpha[13]*fUpwind[23]+fUpwind[15]*alpha[20]+alpha[7]*fUpwind[19]+fUpwind[7]*alpha[19]+alpha[8]*fUpwind[18]+fUpwind[8]*alpha[18]+alpha[0]*fUpwind[17]+fUpwind[0]*alpha[17]+alpha[11]*fUpwind[16]+fUpwind[11]*alpha[16]+alpha[1]*fUpwind[10]+fUpwind[1]*alpha[10]+alpha[2]*fUpwind[9]+fUpwind[2]*alpha[9]+alpha[4]*fUpwind[6]+fUpwind[4]*alpha[6])); 
  Ghat[18] = 0.01178511301977579*(13.41640786499874*alpha[10]*fUpwind[43]+13.41640786499874*(alpha[17]*fUpwind[39]+alpha[4]*fUpwind[38]+alpha[19]*fUpwind[37])+13.41640786499874*(alpha[9]*fUpwind[35]+alpha[26]*fUpwind[34]+alpha[11]*fUpwind[33])+13.41640786499874*alpha[18]*fUpwind[32]+15.0*(alpha[13]*fUpwind[31]+alpha[20]*fUpwind[30]+alpha[5]*fUpwind[29]+alpha[22]*fUpwind[28]+fUpwind[24]*alpha[27]+alpha[2]*fUpwind[26]+fUpwind[2]*alpha[26]+alpha[12]*fUpwind[25]+alpha[14]*fUpwind[23]+fUpwind[15]*alpha[21]+alpha[6]*fUpwind[19]+fUpwind[6]*alpha[19]+alpha[0]*fUpwind[18]+fUpwind[0]*alpha[18]+alpha[8]*fUpwind[17]+fUpwind[8]*alpha[17]+alpha[10]*fUpwind[16]+fUpwind[10]*alpha[16]+alpha[1]*fUpwind[11]+fUpwind[1]*alpha[11]+alpha[3]*fUpwind[9]+fUpwind[3]*alpha[9]+alpha[4]*fUpwind[7]+fUpwind[4]*alpha[7])); 
  Ghat[19] = 0.01178511301977579*(13.41640786499874*alpha[9]*fUpwind[43]+13.41640786499874*(alpha[4]*fUpwind[39]+alpha[17]*fUpwind[38]+alpha[18]*fUpwind[37])+13.41640786499874*(alpha[10]*fUpwind[35]+alpha[11]*fUpwind[34]+alpha[26]*fUpwind[33])+13.41640786499874*alpha[19]*fUpwind[32]+15.0*(alpha[12]*fUpwind[31]+alpha[5]*fUpwind[30]+alpha[20]*fUpwind[29]+alpha[21]*fUpwind[28]+fUpwind[23]*alpha[27]+alpha[1]*fUpwind[26]+fUpwind[1]*alpha[26]+alpha[13]*fUpwind[25]+alpha[14]*fUpwind[24]+fUpwind[15]*alpha[22]+alpha[0]*fUpwind[19]+fUpwind[0]*alpha[19]+alpha[6]*fUpwind[18]+fUpwind[6]*alpha[18]+alpha[7]*fUpwind[17]+fUpwind[7]*alpha[17]+alpha[9]*fUpwind[16]+fUpwind[9]*alpha[16]+alpha[2]*fUpwind[11]+fUpwind[2]*alpha[11]+alpha[3]*fUpwind[10]+fUpwind[3]*alpha[10]+alpha[4]*fUpwind[8]+fUpwind[4]*alpha[8])); 
  Ghat[20] = 0.01178511301977579*(13.41640786499874*alpha[14]*fUpwind[59]+13.41640786499874*(alpha[21]*fUpwind[55]+alpha[22]*fUpwind[54]+alpha[5]*fUpwind[53])+13.41640786499874*(alpha[27]*fUpwind[51]+alpha[12]*fUpwind[50]+alpha[13]*fUpwind[49])+13.41640786499874*alpha[20]*fUpwind[48]+15.0*(alpha[11]*fUpwind[31]+alpha[18]*fUpwind[30]+alpha[19]*fUpwind[29]+alpha[4]*fUpwind[28]+alpha[3]*fUpwind[27]+fUpwind[3]*alpha[27]+fUpwind[25]*alpha[26]+alpha[9]*fUpwind[24]+alpha[10]*fUpwind[23]+alpha[7]*fUpwind[22]+fUpwind[7]*alpha[22]+alpha[8]*fUpwind[21]+fUpwind[8]*alpha[21]+alpha[0]*fUpwind[20]+fUpwind[0]*alpha[20]+fUpwind[15]*alpha[17]+alpha[14]*fUpwind[16]+fUpwind[14]*alpha[16]+alpha[1]*fUpwind[13]+fUpwind[1]*alpha[13]+alpha[2]*fUpwind[12]+fUpwind[2]*alpha[12]+alpha[5]*fUpwind[6]+fUpwind[5]*alpha[6])); 
  Ghat[21] = 0.01178511301977579*(13.41640786499874*alpha[13]*fUpwind[59]+13.41640786499874*(alpha[20]*fUpwind[55]+alpha[5]*fUpwind[54]+alpha[22]*fUpwind[53])+13.41640786499874*(alpha[12]*fUpwind[51]+alpha[27]*fUpwind[50]+alpha[14]*fUpwind[49])+13.41640786499874*alpha[21]*fUpwind[48]+15.0*(alpha[10]*fUpwind[31]+alpha[17]*fUpwind[30]+alpha[4]*fUpwind[29]+alpha[19]*fUpwind[28]+alpha[2]*fUpwind[27]+fUpwind[2]*alpha[27]+fUpwind[24]*alpha[26]+alpha[9]*fUpwind[25]+alpha[11]*fUpwind[23]+alpha[6]*fUpwind[22]+fUpwind[6]*alpha[22]+alpha[0]*fUpwind[21]+fUpwind[0]*alpha[21]+alpha[8]*fUpwind[20]+fUpwind[8]*alpha[20]+fUpwind[15]*alpha[18]+alpha[13]*fUpwind[16]+fUpwind[13]*alpha[16]+alpha[1]*fUpwind[14]+fUpwind[1]*alpha[14]+alpha[3]*fUpwind[12]+fUpwind[3]*alpha[12]+alpha[5]*fUpwind[7]+fUpwind[5]*alpha[7])); 
  Ghat[22] = 0.01178511301977579*(13.41640786499874*alpha[12]*fUpwind[59]+13.41640786499874*(alpha[5]*fUpwind[55]+alpha[20]*fUpwind[54]+alpha[21]*fUpwind[53])+13.41640786499874*(alpha[13]*fUpwind[51]+alpha[14]*fUpwind[50]+alpha[27]*fUpwind[49])+13.41640786499874*alpha[22]*fUpwind[48]+15.0*(alpha[9]*fUpwind[31]+alpha[4]*fUpwind[30]+alpha[17]*fUpwind[29]+alpha[18]*fUpwind[28]+alpha[1]*fUpwind[27]+fUpwind[1]*alpha[27]+fUpwind[23]*alpha[26]+alpha[10]*fUpwind[25]+alpha[11]*fUpwind[24]+alpha[0]*fUpwind[22]+fUpwind[0]*alpha[22]+alpha[6]*fUpwind[21]+fUpwind[6]*alpha[21]+alpha[7]*fUpwind[20]+fUpwind[7]*alpha[20]+fUpwind[15]*alpha[19]+alpha[12]*fUpwind[16]+fUpwind[12]*alpha[16]+alpha[2]*fUpwind[14]+fUpwind[2]*alpha[14]+alpha[3]*fUpwind[13]+fUpwind[3]*alpha[13]+alpha[5]*fUpwind[8]+fUpwind[5]*alpha[8])); 
  Ghat[23] = 0.01178511301977579*(13.41640786499874*alpha[22]*fUpwind[63]+13.41640786499874*(alpha[27]*fUpwind[62]+alpha[14]*fUpwind[61]+alpha[13]*fUpwind[60])+13.41640786499874*(alpha[21]*fUpwind[58]+alpha[20]*fUpwind[57]+alpha[5]*fUpwind[56])+13.41640786499874*alpha[12]*fUpwind[52]+13.41640786499874*alpha[19]*fUpwind[47]+13.41640786499874*(alpha[26]*fUpwind[46]+alpha[11]*fUpwind[45]+alpha[10]*fUpwind[44])+13.41640786499874*(alpha[18]*fUpwind[42]+alpha[17]*fUpwind[41]+alpha[4]*fUpwind[40])+13.41640786499874*alpha[9]*fUpwind[36]+15.0*(alpha[8]*fUpwind[31]+alpha[16]*fUpwind[30]+alpha[3]*fUpwind[29]+alpha[2]*fUpwind[28]+alpha[19]*fUpwind[27]+fUpwind[19]*alpha[27]+alpha[22]*fUpwind[26]+fUpwind[22]*alpha[26]+alpha[7]*fUpwind[25]+alpha[6]*fUpwind[24]+alpha[0]*fUpwind[23]+alpha[11]*fUpwind[21]+fUpwind[11]*alpha[21]+alpha[10]*fUpwind[20]+fUpwind[10]*alpha[20]+alpha[14]*fUpwind[18]+fUpwind[14]*alpha[18]+alpha[13]*fUpwind[17]+fUpwind[13]*alpha[17]+alpha[1]*fUpwind[15]+alpha[4]*fUpwind[12]+fUpwind[4]*alpha[12]+alpha[5]*fUpwind[9]+fUpwind[5]*alpha[9])); 
  Ghat[24] = 0.01178511301977579*(13.41640786499874*alpha[21]*fUpwind[63]+13.41640786499874*(alpha[14]*fUpwind[62]+alpha[27]*fUpwind[61]+alpha[12]*fUpwind[60])+13.41640786499874*(alpha[22]*fUpwind[58]+alpha[5]*fUpwind[57]+alpha[20]*fUpwind[56])+13.41640786499874*alpha[13]*fUpwind[52]+13.41640786499874*alpha[18]*fUpwind[47]+13.41640786499874*(alpha[11]*fUpwind[46]+alpha[26]*fUpwind[45]+alpha[9]*fUpwind[44])+13.41640786499874*(alpha[19]*fUpwind[42]+alpha[4]*fUpwind[41]+alpha[17]*fUpwind[40])+13.41640786499874*alpha[10]*fUpwind[36]+15.0*(alpha[7]*fUpwind[31]+alpha[3]*fUpwind[30]+alpha[16]*fUpwind[29]+alpha[1]*fUpwind[28]+alpha[18]*fUpwind[27]+fUpwind[18]*alpha[27]+alpha[21]*fUpwind[26]+fUpwind[21]*alpha[26]+alpha[8]*fUpwind[25]+alpha[0]*fUpwind[24]+alpha[6]*fUpwind[23]+alpha[11]*fUpwind[22]+fUpwind[11]*alpha[22]+alpha[9]*fUpwind[20]+fUpwind[9]*alpha[20]+alpha[14]*fUpwind[19]+fUpwind[14]*alpha[19]+alpha[12]*fUpwind[17]+fUpwind[12]*alpha[17]+alpha[2]*fUpwind[15]+alpha[4]*fUpwind[13]+fUpwind[4]*alpha[13]+alpha[5]*fUpwind[10]+fUpwind[5]*alpha[10])); 
  Ghat[25] = 0.01178511301977579*(13.41640786499874*alpha[20]*fUpwind[63]+13.41640786499874*(alpha[13]*fUpwind[62]+alpha[12]*fUpwind[61]+alpha[27]*fUpwind[60])+13.41640786499874*(alpha[5]*fUpwind[58]+alpha[22]*fUpwind[57]+alpha[21]*fUpwind[56])+13.41640786499874*alpha[14]*fUpwind[52]+13.41640786499874*alpha[17]*fUpwind[47]+13.41640786499874*(alpha[10]*fUpwind[46]+alpha[9]*fUpwind[45]+alpha[26]*fUpwind[44])+13.41640786499874*(alpha[4]*fUpwind[42]+alpha[19]*fUpwind[41]+alpha[18]*fUpwind[40])+13.41640786499874*alpha[11]*fUpwind[36]+15.0*(alpha[6]*fUpwind[31]+alpha[2]*fUpwind[30]+alpha[1]*fUpwind[29]+alpha[16]*fUpwind[28]+alpha[17]*fUpwind[27]+fUpwind[17]*alpha[27]+alpha[20]*fUpwind[26]+fUpwind[20]*alpha[26]+alpha[0]*fUpwind[25]+alpha[8]*fUpwind[24]+alpha[7]*fUpwind[23]+alpha[10]*fUpwind[22]+fUpwind[10]*alpha[22]+alpha[9]*fUpwind[21]+fUpwind[9]*alpha[21]+alpha[13]*fUpwind[19]+fUpwind[13]*alpha[19]+alpha[12]*fUpwind[18]+fUpwind[12]*alpha[18]+alpha[3]*fUpwind[15]+alpha[4]*fUpwind[14]+fUpwind[4]*alpha[14]+alpha[5]*fUpwind[11]+fUpwind[5]*alpha[11])); 
  Ghat[26] = 0.01178511301977579*(13.41640786499874*alpha[4]*fUpwind[43]+13.41640786499874*(alpha[9]*fUpwind[39]+alpha[10]*fUpwind[38]+alpha[11]*fUpwind[37])+13.41640786499874*(alpha[17]*fUpwind[35]+alpha[18]*fUpwind[34]+alpha[19]*fUpwind[33])+13.41640786499874*alpha[26]*fUpwind[32]+15.0*(alpha[5]*fUpwind[31]+alpha[12]*fUpwind[30]+alpha[13]*fUpwind[29]+alpha[14]*fUpwind[28]+fUpwind[15]*alpha[27]+alpha[0]*fUpwind[26]+fUpwind[0]*alpha[26]+alpha[20]*fUpwind[25]+alpha[21]*fUpwind[24]+alpha[22]*fUpwind[23]+alpha[1]*fUpwind[19]+fUpwind[1]*alpha[19]+alpha[2]*fUpwind[18]+fUpwind[2]*alpha[18]+alpha[3]*fUpwind[17]+fUpwind[3]*alpha[17]+alpha[4]*fUpwind[16]+fUpwind[4]*alpha[16]+alpha[6]*fUpwind[11]+fUpwind[6]*alpha[11]+alpha[7]*fUpwind[10]+fUpwind[7]*alpha[10]+alpha[8]*fUpwind[9]+fUpwind[8]*alpha[9])); 
  Ghat[27] = 0.01178511301977579*(13.41640786499874*alpha[5]*fUpwind[59]+13.41640786499874*(alpha[12]*fUpwind[55]+alpha[13]*fUpwind[54]+alpha[14]*fUpwind[53])+13.41640786499874*(alpha[20]*fUpwind[51]+alpha[21]*fUpwind[50]+alpha[22]*fUpwind[49])+13.41640786499874*alpha[27]*fUpwind[48]+15.0*(alpha[4]*fUpwind[31]+alpha[9]*fUpwind[30]+alpha[10]*fUpwind[29]+alpha[11]*fUpwind[28]+alpha[0]*fUpwind[27]+fUpwind[0]*alpha[27]+fUpwind[15]*alpha[26]+alpha[17]*fUpwind[25]+alpha[18]*fUpwind[24]+alpha[19]*fUpwind[23]+alpha[1]*fUpwind[22]+fUpwind[1]*alpha[22]+alpha[2]*fUpwind[21]+fUpwind[2]*alpha[21]+alpha[3]*fUpwind[20]+fUpwind[3]*alpha[20]+alpha[5]*fUpwind[16]+fUpwind[5]*alpha[16]+alpha[6]*fUpwind[14]+fUpwind[6]*alpha[14]+alpha[7]*fUpwind[13]+fUpwind[7]*alpha[13]+alpha[8]*fUpwind[12]+fUpwind[8]*alpha[12])); 
  Ghat[28] = 0.01178511301977579*(13.41640786499874*alpha[14]*fUpwind[63]+13.41640786499874*(alpha[21]*fUpwind[62]+alpha[22]*fUpwind[61]+alpha[5]*fUpwind[60])+13.41640786499874*(alpha[27]*fUpwind[58]+alpha[12]*fUpwind[57]+alpha[13]*fUpwind[56])+13.41640786499874*alpha[20]*fUpwind[52]+13.41640786499874*alpha[11]*fUpwind[47]+13.41640786499874*(alpha[18]*fUpwind[46]+alpha[19]*fUpwind[45]+alpha[4]*fUpwind[44])+13.41640786499874*(alpha[26]*fUpwind[42]+alpha[9]*fUpwind[41]+alpha[10]*fUpwind[40])+13.41640786499874*alpha[17]*fUpwind[36]+15.0*(alpha[3]*fUpwind[31]+alpha[7]*fUpwind[30]+alpha[8]*fUpwind[29]+alpha[0]*fUpwind[28]+alpha[11]*fUpwind[27]+fUpwind[11]*alpha[27]+alpha[14]*fUpwind[26]+fUpwind[14]*alpha[26]+alpha[16]*fUpwind[25]+alpha[1]*fUpwind[24]+alpha[2]*fUpwind[23]+alpha[18]*fUpwind[22]+fUpwind[18]*alpha[22]+alpha[19]*fUpwind[21]+fUpwind[19]*alpha[21]+alpha[4]*fUpwind[20]+fUpwind[4]*alpha[20]+alpha[5]*fUpwind[17]+fUpwind[5]*alpha[17]+alpha[6]*fUpwind[15]+alpha[9]*fUpwind[13]+fUpwind[9]*alpha[13]+alpha[10]*fUpwind[12]+fUpwind[10]*alpha[12])); 
  Ghat[29] = 0.01178511301977579*(13.41640786499874*alpha[13]*fUpwind[63]+13.41640786499874*(alpha[20]*fUpwind[62]+alpha[5]*fUpwind[61]+alpha[22]*fUpwind[60])+13.41640786499874*(alpha[12]*fUpwind[58]+alpha[27]*fUpwind[57]+alpha[14]*fUpwind[56])+13.41640786499874*alpha[21]*fUpwind[52]+13.41640786499874*alpha[10]*fUpwind[47]+13.41640786499874*(alpha[17]*fUpwind[46]+alpha[4]*fUpwind[45]+alpha[19]*fUpwind[44])+13.41640786499874*(alpha[9]*fUpwind[42]+alpha[26]*fUpwind[41]+alpha[11]*fUpwind[40])+13.41640786499874*alpha[18]*fUpwind[36]+15.0*(alpha[2]*fUpwind[31]+alpha[6]*fUpwind[30]+alpha[0]*fUpwind[29]+alpha[8]*fUpwind[28]+alpha[10]*fUpwind[27]+fUpwind[10]*alpha[27]+alpha[13]*fUpwind[26]+fUpwind[13]*alpha[26]+alpha[1]*fUpwind[25]+alpha[16]*fUpwind[24]+alpha[3]*fUpwind[23]+alpha[17]*fUpwind[22]+fUpwind[17]*alpha[22]+alpha[4]*fUpwind[21]+fUpwind[4]*alpha[21]+alpha[19]*fUpwind[20]+fUpwind[19]*alpha[20]+alpha[5]*fUpwind[18]+fUpwind[5]*alpha[18]+alpha[7]*fUpwind[15]+alpha[9]*fUpwind[14]+fUpwind[9]*alpha[14]+alpha[11]*fUpwind[12]+fUpwind[11]*alpha[12])); 
  Ghat[30] = 0.01178511301977579*(13.41640786499874*alpha[12]*fUpwind[63]+13.41640786499874*(alpha[5]*fUpwind[62]+alpha[20]*fUpwind[61]+alpha[21]*fUpwind[60])+13.41640786499874*(alpha[13]*fUpwind[58]+alpha[14]*fUpwind[57]+alpha[27]*fUpwind[56])+13.41640786499874*alpha[22]*fUpwind[52]+13.41640786499874*alpha[9]*fUpwind[47]+13.41640786499874*(alpha[4]*fUpwind[46]+alpha[17]*fUpwind[45]+alpha[18]*fUpwind[44])+13.41640786499874*(alpha[10]*fUpwind[42]+alpha[11]*fUpwind[41]+alpha[26]*fUpwind[40])+13.41640786499874*alpha[19]*fUpwind[36]+15.0*(alpha[1]*fUpwind[31]+alpha[0]*fUpwind[30]+alpha[6]*fUpwind[29]+alpha[7]*fUpwind[28]+alpha[9]*fUpwind[27]+fUpwind[9]*alpha[27]+alpha[12]*fUpwind[26]+fUpwind[12]*alpha[26]+alpha[2]*fUpwind[25]+alpha[3]*fUpwind[24]+alpha[16]*fUpwind[23]+alpha[4]*fUpwind[22]+fUpwind[4]*alpha[22]+alpha[17]*fUpwind[21]+fUpwind[17]*alpha[21]+alpha[18]*fUpwind[20]+fUpwind[18]*alpha[20]+alpha[5]*fUpwind[19]+fUpwind[5]*alpha[19]+alpha[8]*fUpwind[15]+alpha[10]*fUpwind[14]+fUpwind[10]*alpha[14]+alpha[11]*fUpwind[13]+fUpwind[11]*alpha[13])); 
  Ghat[31] = 0.01178511301977579*(13.41640786499874*alpha[5]*fUpwind[63]+13.41640786499874*(alpha[12]*fUpwind[62]+alpha[13]*fUpwind[61]+alpha[14]*fUpwind[60])+13.41640786499874*(alpha[20]*fUpwind[58]+alpha[21]*fUpwind[57]+alpha[22]*fUpwind[56])+13.41640786499874*alpha[27]*fUpwind[52]+13.41640786499874*alpha[4]*fUpwind[47]+13.41640786499874*(alpha[9]*fUpwind[46]+alpha[10]*fUpwind[45]+alpha[11]*fUpwind[44])+13.41640786499874*(alpha[17]*fUpwind[42]+alpha[18]*fUpwind[41]+alpha[19]*fUpwind[40])+13.41640786499874*alpha[26]*fUpwind[36]+15.0*(alpha[0]*fUpwind[31]+alpha[1]*fUpwind[30]+alpha[2]*fUpwind[29]+alpha[3]*fUpwind[28]+alpha[4]*fUpwind[27]+fUpwind[4]*alpha[27]+alpha[5]*fUpwind[26]+fUpwind[5]*alpha[26]+alpha[6]*fUpwind[25]+alpha[7]*fUpwind[24]+alpha[8]*fUpwind[23]+alpha[9]*fUpwind[22]+fUpwind[9]*alpha[22]+alpha[10]*fUpwind[21]+fUpwind[10]*alpha[21]+alpha[11]*fUpwind[20]+fUpwind[11]*alpha[20]+alpha[12]*fUpwind[19]+fUpwind[12]*alpha[19]+alpha[13]*fUpwind[18]+fUpwind[13]*alpha[18]+alpha[14]*fUpwind[17]+fUpwind[14]*alpha[17]+fUpwind[15]*alpha[16])); 
  Ghat[32] = 0.01178511301977579*(15.0*alpha[27]*fUpwind[47]+15.0*(alpha[22]*fUpwind[46]+alpha[21]*fUpwind[45]+alpha[20]*fUpwind[44]+alpha[16]*fUpwind[43])+15.0*(alpha[14]*fUpwind[42]+alpha[13]*fUpwind[41]+alpha[12]*fUpwind[40]+alpha[8]*fUpwind[39]+alpha[7]*fUpwind[38]+alpha[6]*fUpwind[37])+15.0*(alpha[5]*fUpwind[36]+alpha[3]*fUpwind[35]+alpha[2]*fUpwind[34]+alpha[1]*fUpwind[33])+15.0*alpha[0]*fUpwind[32]+13.41640786499874*(alpha[26]*fUpwind[26]+alpha[19]*fUpwind[19]+alpha[18]*fUpwind[18]+alpha[17]*fUpwind[17]+alpha[11]*fUpwind[11]+alpha[10]*fUpwind[10]+alpha[9]*fUpwind[9]+alpha[4]*fUpwind[4])); 
  Ghat[33] = 0.01178511301977579*(15.0*alpha[22]*fUpwind[47]+15.0*(alpha[27]*fUpwind[46]+alpha[14]*fUpwind[45]+alpha[13]*fUpwind[44]+alpha[8]*fUpwind[43])+15.0*(alpha[21]*fUpwind[42]+alpha[20]*fUpwind[41]+alpha[5]*fUpwind[40]+alpha[16]*fUpwind[39]+alpha[3]*fUpwind[38]+alpha[2]*fUpwind[37])+15.0*(alpha[12]*fUpwind[36]+alpha[7]*fUpwind[35]+alpha[6]*fUpwind[34]+alpha[0]*fUpwind[33])+15.0*alpha[1]*fUpwind[32]+13.41640786499874*(alpha[19]*fUpwind[26]+fUpwind[19]*alpha[26]+alpha[11]*fUpwind[18]+fUpwind[11]*alpha[18]+alpha[10]*fUpwind[17]+fUpwind[10]*alpha[17]+alpha[4]*fUpwind[9]+fUpwind[4]*alpha[9])); 
  Ghat[34] = 0.01178511301977579*(15.0*alpha[21]*fUpwind[47]+15.0*(alpha[14]*fUpwind[46]+alpha[27]*fUpwind[45]+alpha[12]*fUpwind[44]+alpha[7]*fUpwind[43])+15.0*(alpha[22]*fUpwind[42]+alpha[5]*fUpwind[41]+alpha[20]*fUpwind[40]+alpha[3]*fUpwind[39]+alpha[16]*fUpwind[38]+alpha[1]*fUpwind[37])+15.0*(alpha[13]*fUpwind[36]+alpha[8]*fUpwind[35]+alpha[0]*fUpwind[34]+alpha[6]*fUpwind[33])+15.0*alpha[2]*fUpwind[32]+13.41640786499874*(alpha[18]*fUpwind[26]+fUpwind[18]*alpha[26]+alpha[11]*fUpwind[19]+fUpwind[11]*alpha[19]+alpha[9]*fUpwind[17]+fUpwind[9]*alpha[17]+alpha[4]*fUpwind[10]+fUpwind[4]*alpha[10])); 
  Ghat[35] = 0.01178511301977579*(15.0*alpha[20]*fUpwind[47]+15.0*(alpha[13]*fUpwind[46]+alpha[12]*fUpwind[45]+alpha[27]*fUpwind[44]+alpha[6]*fUpwind[43])+15.0*(alpha[5]*fUpwind[42]+alpha[22]*fUpwind[41]+alpha[21]*fUpwind[40]+alpha[2]*fUpwind[39]+alpha[1]*fUpwind[38]+alpha[16]*fUpwind[37])+15.0*(alpha[14]*fUpwind[36]+alpha[0]*fUpwind[35]+alpha[8]*fUpwind[34]+alpha[7]*fUpwind[33])+15.0*alpha[3]*fUpwind[32]+13.41640786499874*(alpha[17]*fUpwind[26]+fUpwind[17]*alpha[26]+alpha[10]*fUpwind[19]+fUpwind[10]*alpha[19]+alpha[9]*fUpwind[18]+fUpwind[9]*alpha[18]+alpha[4]*fUpwind[11]+fUpwind[4]*alpha[11])); 
  Ghat[36] = 0.01178511301977579*(15.0*alpha[16]*fUpwind[47]+15.0*(alpha[8]*fUpwind[46]+alpha[7]*fUpwind[45]+alpha[6]*fUpwind[44]+alpha[27]*fUpwind[43])+15.0*(alpha[3]*fUpwind[42]+alpha[2]*fUpwind[41]+alpha[1]*fUpwind[40]+alpha[22]*fUpwind[39]+alpha[21]*fUpwind[38]+alpha[20]*fUpwind[37])+15.0*(alpha[0]*fUpwind[36]+alpha[14]*fUpwind[35]+alpha[13]*fUpwind[34]+alpha[12]*fUpwind[33])+15.0*alpha[5]*fUpwind[32]+13.41640786499874*(alpha[26]*fUpwind[31]+alpha[19]*fUpwind[30]+alpha[18]*fUpwind[29]+alpha[17]*fUpwind[28]+alpha[11]*fUpwind[25]+alpha[10]*fUpwind[24]+alpha[9]*fUpwind[23]+alpha[4]*fUpwind[15])); 
  Ghat[37] = 0.01178511301977579*(15.0*alpha[14]*fUpwind[47]+15.0*(alpha[21]*fUpwind[46]+alpha[22]*fUpwind[45]+alpha[5]*fUpwind[44]+alpha[3]*fUpwind[43])+15.0*(alpha[27]*fUpwind[42]+alpha[12]*fUpwind[41]+alpha[13]*fUpwind[40]+alpha[7]*fUpwind[39]+alpha[8]*fUpwind[38]+alpha[0]*fUpwind[37])+15.0*(alpha[20]*fUpwind[36]+alpha[16]*fUpwind[35]+alpha[1]*fUpwind[34]+alpha[2]*fUpwind[33])+15.0*alpha[6]*fUpwind[32]+13.41640786499874*(alpha[11]*fUpwind[26]+fUpwind[11]*alpha[26]+alpha[18]*fUpwind[19]+fUpwind[18]*alpha[19]+alpha[4]*fUpwind[17]+fUpwind[4]*alpha[17]+alpha[9]*fUpwind[10]+fUpwind[9]*alpha[10])); 
  Ghat[38] = 0.01178511301977579*(15.0*alpha[13]*fUpwind[47]+15.0*(alpha[20]*fUpwind[46]+alpha[5]*fUpwind[45]+alpha[22]*fUpwind[44]+alpha[2]*fUpwind[43])+15.0*(alpha[12]*fUpwind[42]+alpha[27]*fUpwind[41]+alpha[14]*fUpwind[40]+alpha[6]*fUpwind[39]+alpha[0]*fUpwind[38]+alpha[8]*fUpwind[37])+15.0*(alpha[21]*fUpwind[36]+alpha[1]*fUpwind[35]+alpha[16]*fUpwind[34]+alpha[3]*fUpwind[33])+15.0*alpha[7]*fUpwind[32]+13.41640786499874*(alpha[10]*fUpwind[26]+fUpwind[10]*alpha[26]+alpha[17]*fUpwind[19]+fUpwind[17]*alpha[19]+alpha[4]*fUpwind[18]+fUpwind[4]*alpha[18]+alpha[9]*fUpwind[11]+fUpwind[9]*alpha[11])); 
  Ghat[39] = 0.01178511301977579*(15.0*alpha[12]*fUpwind[47]+15.0*(alpha[5]*fUpwind[46]+alpha[20]*fUpwind[45]+alpha[21]*fUpwind[44]+alpha[1]*fUpwind[43])+15.0*(alpha[13]*fUpwind[42]+alpha[14]*fUpwind[41]+alpha[27]*fUpwind[40]+alpha[0]*fUpwind[39]+alpha[6]*fUpwind[38]+alpha[7]*fUpwind[37])+15.0*(alpha[22]*fUpwind[36]+alpha[2]*fUpwind[35]+alpha[3]*fUpwind[34]+alpha[16]*fUpwind[33])+15.0*alpha[8]*fUpwind[32]+13.41640786499874*(alpha[9]*fUpwind[26]+fUpwind[9]*alpha[26]+alpha[4]*fUpwind[19]+fUpwind[4]*alpha[19]+alpha[17]*fUpwind[18]+fUpwind[17]*alpha[18]+alpha[10]*fUpwind[11]+fUpwind[10]*alpha[11])); 
  Ghat[40] = 0.01178511301977579*(15.0*alpha[8]*fUpwind[47]+15.0*(alpha[16]*fUpwind[46]+alpha[3]*fUpwind[45]+alpha[2]*fUpwind[44]+alpha[22]*fUpwind[43])+15.0*(alpha[7]*fUpwind[42]+alpha[6]*fUpwind[41]+alpha[0]*fUpwind[40]+alpha[27]*fUpwind[39]+alpha[14]*fUpwind[38]+alpha[13]*fUpwind[37])+15.0*(alpha[1]*fUpwind[36]+alpha[21]*fUpwind[35]+alpha[20]*fUpwind[34]+alpha[5]*fUpwind[33])+15.0*alpha[12]*fUpwind[32]+13.41640786499874*(alpha[19]*fUpwind[31]+alpha[26]*fUpwind[30]+alpha[11]*fUpwind[29]+alpha[10]*fUpwind[28]+alpha[18]*fUpwind[25]+alpha[17]*fUpwind[24]+alpha[4]*fUpwind[23]+alpha[9]*fUpwind[15])); 
  Ghat[41] = 0.01178511301977579*(15.0*alpha[7]*fUpwind[47]+15.0*(alpha[3]*fUpwind[46]+alpha[16]*fUpwind[45]+alpha[1]*fUpwind[44]+alpha[21]*fUpwind[43])+15.0*(alpha[8]*fUpwind[42]+alpha[0]*fUpwind[41]+alpha[6]*fUpwind[40]+alpha[14]*fUpwind[39]+alpha[27]*fUpwind[38]+alpha[12]*fUpwind[37])+15.0*(alpha[2]*fUpwind[36]+alpha[22]*fUpwind[35]+alpha[5]*fUpwind[34]+alpha[20]*fUpwind[33])+15.0*alpha[13]*fUpwind[32]+13.41640786499874*(alpha[18]*fUpwind[31]+alpha[11]*fUpwind[30]+alpha[26]*fUpwind[29]+alpha[9]*fUpwind[28]+alpha[19]*fUpwind[25]+alpha[4]*fUpwind[24]+alpha[17]*fUpwind[23]+alpha[10]*fUpwind[15])); 
  Ghat[42] = 0.01178511301977579*(15.0*alpha[6]*fUpwind[47]+15.0*(alpha[2]*fUpwind[46]+alpha[1]*fUpwind[45]+alpha[16]*fUpwind[44]+alpha[20]*fUpwind[43])+15.0*(alpha[0]*fUpwind[42]+alpha[8]*fUpwind[41]+alpha[7]*fUpwind[40]+alpha[13]*fUpwind[39]+alpha[12]*fUpwind[38]+alpha[27]*fUpwind[37])+15.0*(alpha[3]*fUpwind[36]+alpha[5]*fUpwind[35]+alpha[22]*fUpwind[34]+alpha[21]*fUpwind[33])+15.0*alpha[14]*fUpwind[32]+13.41640786499874*(alpha[17]*fUpwind[31]+alpha[10]*fUpwind[30]+alpha[9]*fUpwind[29]+alpha[26]*fUpwind[28]+alpha[4]*fUpwind[25]+alpha[19]*fUpwind[24]+alpha[18]*fUpwind[23]+alpha[11]*fUpwind[15])); 
  Ghat[43] = 0.01178511301977579*(15.0*alpha[5]*fUpwind[47]+15.0*(alpha[12]*fUpwind[46]+alpha[13]*fUpwind[45]+alpha[14]*fUpwind[44]+alpha[0]*fUpwind[43])+15.0*(alpha[20]*fUpwind[42]+alpha[21]*fUpwind[41]+alpha[22]*fUpwind[40]+alpha[1]*fUpwind[39]+alpha[2]*fUpwind[38]+alpha[3]*fUpwind[37])+15.0*(alpha[27]*fUpwind[36]+alpha[6]*fUpwind[35]+alpha[7]*fUpwind[34]+alpha[8]*fUpwind[33])+15.0*alpha[16]*fUpwind[32]+13.41640786499874*(alpha[4]*fUpwind[26]+fUpwind[4]*alpha[26]+alpha[9]*fUpwind[19]+fUpwind[9]*alpha[19]+alpha[10]*fUpwind[18]+fUpwind[10]*alpha[18]+alpha[11]*fUpwind[17]+fUpwind[11]*alpha[17])); 
  Ghat[44] = 0.01178511301977579*(15.0*alpha[3]*fUpwind[47]+15.0*(alpha[7]*fUpwind[46]+alpha[8]*fUpwind[45]+alpha[0]*fUpwind[44]+alpha[14]*fUpwind[43])+15.0*(alpha[16]*fUpwind[42]+alpha[1]*fUpwind[41]+alpha[2]*fUpwind[40]+alpha[21]*fUpwind[39]+alpha[22]*fUpwind[38]+alpha[5]*fUpwind[37])+15.0*(alpha[6]*fUpwind[36]+alpha[27]*fUpwind[35]+alpha[12]*fUpwind[34]+alpha[13]*fUpwind[33])+15.0*alpha[20]*fUpwind[32]+13.41640786499874*(alpha[11]*fUpwind[31]+alpha[18]*fUpwind[30]+alpha[19]*fUpwind[29]+alpha[4]*fUpwind[28]+fUpwind[25]*alpha[26]+alpha[9]*fUpwind[24]+alpha[10]*fUpwind[23]+fUpwind[15]*alpha[17])); 
  Ghat[45] = 0.01178511301977579*(15.0*alpha[2]*fUpwind[47]+15.0*(alpha[6]*fUpwind[46]+alpha[0]*fUpwind[45]+alpha[8]*fUpwind[44]+alpha[13]*fUpwind[43])+15.0*(alpha[1]*fUpwind[42]+alpha[16]*fUpwind[41]+alpha[3]*fUpwind[40]+alpha[20]*fUpwind[39]+alpha[5]*fUpwind[38]+alpha[22]*fUpwind[37])+15.0*(alpha[7]*fUpwind[36]+alpha[12]*fUpwind[35]+alpha[27]*fUpwind[34]+alpha[14]*fUpwind[33])+15.0*alpha[21]*fUpwind[32]+13.41640786499874*(alpha[10]*fUpwind[31]+alpha[17]*fUpwind[30]+alpha[4]*fUpwind[29]+alpha[19]*fUpwind[28]+fUpwind[24]*alpha[26]+alpha[9]*fUpwind[25]+alpha[11]*fUpwind[23]+fUpwind[15]*alpha[18])); 
  Ghat[46] = 0.01178511301977579*(15.0*alpha[1]*fUpwind[47]+15.0*(alpha[0]*fUpwind[46]+alpha[6]*fUpwind[45]+alpha[7]*fUpwind[44]+alpha[12]*fUpwind[43])+15.0*(alpha[2]*fUpwind[42]+alpha[3]*fUpwind[41]+alpha[16]*fUpwind[40]+alpha[5]*fUpwind[39]+alpha[20]*fUpwind[38]+alpha[21]*fUpwind[37])+15.0*(alpha[8]*fUpwind[36]+alpha[13]*fUpwind[35]+alpha[14]*fUpwind[34]+alpha[27]*fUpwind[33])+15.0*alpha[22]*fUpwind[32]+13.41640786499874*(alpha[9]*fUpwind[31]+alpha[4]*fUpwind[30]+alpha[17]*fUpwind[29]+alpha[18]*fUpwind[28]+fUpwind[23]*alpha[26]+alpha[10]*fUpwind[25]+alpha[11]*fUpwind[24]+fUpwind[15]*alpha[19])); 
  Ghat[47] = 0.01178511301977579*(15.0*alpha[0]*fUpwind[47]+15.0*(alpha[1]*fUpwind[46]+alpha[2]*fUpwind[45]+alpha[3]*fUpwind[44]+alpha[5]*fUpwind[43])+15.0*(alpha[6]*fUpwind[42]+alpha[7]*fUpwind[41]+alpha[8]*fUpwind[40]+alpha[12]*fUpwind[39]+alpha[13]*fUpwind[38]+alpha[14]*fUpwind[37])+15.0*(alpha[16]*fUpwind[36]+alpha[20]*fUpwind[35]+alpha[21]*fUpwind[34]+alpha[22]*fUpwind[33])+15.0*alpha[27]*fUpwind[32]+13.41640786499874*(alpha[4]*fUpwind[31]+alpha[9]*fUpwind[30]+alpha[10]*fUpwind[29]+alpha[11]*fUpwind[28]+fUpwind[15]*alpha[26]+alpha[17]*fUpwind[25]+alpha[18]*fUpwind[24]+alpha[19]*fUpwind[23])); 
  Ghat[48] = 0.01178511301977579*(15.0*alpha[26]*fUpwind[63]+15.0*(alpha[19]*fUpwind[62]+alpha[18]*fUpwind[61]+alpha[17]*fUpwind[60]+alpha[16]*fUpwind[59])+15.0*(alpha[11]*fUpwind[58]+alpha[10]*fUpwind[57]+alpha[9]*fUpwind[56]+alpha[8]*fUpwind[55]+alpha[7]*fUpwind[54]+alpha[6]*fUpwind[53])+15.0*(alpha[4]*fUpwind[52]+alpha[3]*fUpwind[51]+alpha[2]*fUpwind[50]+alpha[1]*fUpwind[49])+15.0*alpha[0]*fUpwind[48]+13.41640786499874*(alpha[27]*fUpwind[27]+alpha[22]*fUpwind[22]+alpha[21]*fUpwind[21]+alpha[20]*fUpwind[20]+alpha[14]*fUpwind[14]+alpha[13]*fUpwind[13]+alpha[12]*fUpwind[12]+alpha[5]*fUpwind[5])); 
  Ghat[49] = 0.01178511301977579*(15.0*alpha[19]*fUpwind[63]+15.0*(alpha[26]*fUpwind[62]+alpha[11]*fUpwind[61]+alpha[10]*fUpwind[60]+alpha[8]*fUpwind[59])+15.0*(alpha[18]*fUpwind[58]+alpha[17]*fUpwind[57]+alpha[4]*fUpwind[56]+alpha[16]*fUpwind[55]+alpha[3]*fUpwind[54]+alpha[2]*fUpwind[53])+15.0*(alpha[9]*fUpwind[52]+alpha[7]*fUpwind[51]+alpha[6]*fUpwind[50]+alpha[0]*fUpwind[49])+15.0*alpha[1]*fUpwind[48]+13.41640786499874*(alpha[22]*fUpwind[27]+fUpwind[22]*alpha[27]+alpha[14]*fUpwind[21]+fUpwind[14]*alpha[21]+alpha[13]*fUpwind[20]+fUpwind[13]*alpha[20]+alpha[5]*fUpwind[12]+fUpwind[5]*alpha[12])); 
  Ghat[50] = 0.01178511301977579*(15.0*alpha[18]*fUpwind[63]+15.0*(alpha[11]*fUpwind[62]+alpha[26]*fUpwind[61]+alpha[9]*fUpwind[60]+alpha[7]*fUpwind[59])+15.0*(alpha[19]*fUpwind[58]+alpha[4]*fUpwind[57]+alpha[17]*fUpwind[56]+alpha[3]*fUpwind[55]+alpha[16]*fUpwind[54]+alpha[1]*fUpwind[53])+15.0*(alpha[10]*fUpwind[52]+alpha[8]*fUpwind[51]+alpha[0]*fUpwind[50]+alpha[6]*fUpwind[49])+15.0*alpha[2]*fUpwind[48]+13.41640786499874*(alpha[21]*fUpwind[27]+fUpwind[21]*alpha[27]+alpha[14]*fUpwind[22]+fUpwind[14]*alpha[22]+alpha[12]*fUpwind[20]+fUpwind[12]*alpha[20]+alpha[5]*fUpwind[13]+fUpwind[5]*alpha[13])); 
  Ghat[51] = 0.01178511301977579*(15.0*alpha[17]*fUpwind[63]+15.0*(alpha[10]*fUpwind[62]+alpha[9]*fUpwind[61]+alpha[26]*fUpwind[60]+alpha[6]*fUpwind[59])+15.0*(alpha[4]*fUpwind[58]+alpha[19]*fUpwind[57]+alpha[18]*fUpwind[56]+alpha[2]*fUpwind[55]+alpha[1]*fUpwind[54]+alpha[16]*fUpwind[53])+15.0*(alpha[11]*fUpwind[52]+alpha[0]*fUpwind[51]+alpha[8]*fUpwind[50]+alpha[7]*fUpwind[49])+15.0*alpha[3]*fUpwind[48]+13.41640786499874*(alpha[20]*fUpwind[27]+fUpwind[20]*alpha[27]+alpha[13]*fUpwind[22]+fUpwind[13]*alpha[22]+alpha[12]*fUpwind[21]+fUpwind[12]*alpha[21]+alpha[5]*fUpwind[14]+fUpwind[5]*alpha[14])); 
  Ghat[52] = 0.01178511301977579*(15.0*alpha[16]*fUpwind[63]+15.0*(alpha[8]*fUpwind[62]+alpha[7]*fUpwind[61]+alpha[6]*fUpwind[60]+alpha[26]*fUpwind[59])+15.0*(alpha[3]*fUpwind[58]+alpha[2]*fUpwind[57]+alpha[1]*fUpwind[56]+alpha[19]*fUpwind[55]+alpha[18]*fUpwind[54]+alpha[17]*fUpwind[53])+15.0*(alpha[0]*fUpwind[52]+alpha[11]*fUpwind[51]+alpha[10]*fUpwind[50]+alpha[9]*fUpwind[49])+15.0*alpha[4]*fUpwind[48]+13.41640786499874*(alpha[27]*fUpwind[31]+alpha[22]*fUpwind[30]+alpha[21]*fUpwind[29]+alpha[20]*fUpwind[28]+alpha[14]*fUpwind[25]+alpha[13]*fUpwind[24]+alpha[12]*fUpwind[23]+alpha[5]*fUpwind[15])); 
  Ghat[53] = 0.01178511301977579*(15.0*alpha[11]*fUpwind[63]+15.0*(alpha[18]*fUpwind[62]+alpha[19]*fUpwind[61]+alpha[4]*fUpwind[60]+alpha[3]*fUpwind[59])+15.0*(alpha[26]*fUpwind[58]+alpha[9]*fUpwind[57]+alpha[10]*fUpwind[56]+alpha[7]*fUpwind[55]+alpha[8]*fUpwind[54]+alpha[0]*fUpwind[53])+15.0*(alpha[17]*fUpwind[52]+alpha[16]*fUpwind[51]+alpha[1]*fUpwind[50]+alpha[2]*fUpwind[49])+15.0*alpha[6]*fUpwind[48]+13.41640786499874*(alpha[14]*fUpwind[27]+fUpwind[14]*alpha[27]+alpha[21]*fUpwind[22]+fUpwind[21]*alpha[22]+alpha[5]*fUpwind[20]+fUpwind[5]*alpha[20]+alpha[12]*fUpwind[13]+fUpwind[12]*alpha[13])); 
  Ghat[54] = 0.01178511301977579*(15.0*alpha[10]*fUpwind[63]+15.0*(alpha[17]*fUpwind[62]+alpha[4]*fUpwind[61]+alpha[19]*fUpwind[60]+alpha[2]*fUpwind[59])+15.0*(alpha[9]*fUpwind[58]+alpha[26]*fUpwind[57]+alpha[11]*fUpwind[56]+alpha[6]*fUpwind[55]+alpha[0]*fUpwind[54]+alpha[8]*fUpwind[53])+15.0*(alpha[18]*fUpwind[52]+alpha[1]*fUpwind[51]+alpha[16]*fUpwind[50]+alpha[3]*fUpwind[49])+15.0*alpha[7]*fUpwind[48]+13.41640786499874*(alpha[13]*fUpwind[27]+fUpwind[13]*alpha[27]+alpha[20]*fUpwind[22]+fUpwind[20]*alpha[22]+alpha[5]*fUpwind[21]+fUpwind[5]*alpha[21]+alpha[12]*fUpwind[14]+fUpwind[12]*alpha[14])); 
  Ghat[55] = 0.01178511301977579*(15.0*alpha[9]*fUpwind[63]+15.0*(alpha[4]*fUpwind[62]+alpha[17]*fUpwind[61]+alpha[18]*fUpwind[60]+alpha[1]*fUpwind[59])+15.0*(alpha[10]*fUpwind[58]+alpha[11]*fUpwind[57]+alpha[26]*fUpwind[56]+alpha[0]*fUpwind[55]+alpha[6]*fUpwind[54]+alpha[7]*fUpwind[53])+15.0*(alpha[19]*fUpwind[52]+alpha[2]*fUpwind[51]+alpha[3]*fUpwind[50]+alpha[16]*fUpwind[49])+15.0*alpha[8]*fUpwind[48]+13.41640786499874*(alpha[12]*fUpwind[27]+fUpwind[12]*alpha[27]+alpha[5]*fUpwind[22]+fUpwind[5]*alpha[22]+alpha[20]*fUpwind[21]+fUpwind[20]*alpha[21]+alpha[13]*fUpwind[14]+fUpwind[13]*alpha[14])); 
  Ghat[56] = 0.01178511301977579*(15.0*alpha[8]*fUpwind[63]+15.0*(alpha[16]*fUpwind[62]+alpha[3]*fUpwind[61]+alpha[2]*fUpwind[60]+alpha[19]*fUpwind[59])+15.0*(alpha[7]*fUpwind[58]+alpha[6]*fUpwind[57]+alpha[0]*fUpwind[56]+alpha[26]*fUpwind[55]+alpha[11]*fUpwind[54]+alpha[10]*fUpwind[53])+15.0*(alpha[1]*fUpwind[52]+alpha[18]*fUpwind[51]+alpha[17]*fUpwind[50]+alpha[4]*fUpwind[49])+15.0*alpha[9]*fUpwind[48]+13.41640786499874*(alpha[22]*fUpwind[31]+alpha[27]*fUpwind[30]+alpha[14]*fUpwind[29]+alpha[13]*fUpwind[28]+alpha[21]*fUpwind[25]+alpha[20]*fUpwind[24]+alpha[5]*fUpwind[23]+alpha[12]*fUpwind[15])); 
  Ghat[57] = 0.01178511301977579*(15.0*alpha[7]*fUpwind[63]+15.0*(alpha[3]*fUpwind[62]+alpha[16]*fUpwind[61]+alpha[1]*fUpwind[60]+alpha[18]*fUpwind[59])+15.0*(alpha[8]*fUpwind[58]+alpha[0]*fUpwind[57]+alpha[6]*fUpwind[56]+alpha[11]*fUpwind[55]+alpha[26]*fUpwind[54]+alpha[9]*fUpwind[53])+15.0*(alpha[2]*fUpwind[52]+alpha[19]*fUpwind[51]+alpha[4]*fUpwind[50]+alpha[17]*fUpwind[49])+15.0*alpha[10]*fUpwind[48]+13.41640786499874*(alpha[21]*fUpwind[31]+alpha[14]*fUpwind[30]+alpha[27]*fUpwind[29]+alpha[12]*fUpwind[28]+alpha[22]*fUpwind[25]+alpha[5]*fUpwind[24]+alpha[20]*fUpwind[23]+alpha[13]*fUpwind[15])); 
  Ghat[58] = 0.01178511301977579*(15.0*alpha[6]*fUpwind[63]+15.0*(alpha[2]*fUpwind[62]+alpha[1]*fUpwind[61]+alpha[16]*fUpwind[60]+alpha[17]*fUpwind[59])+15.0*(alpha[0]*fUpwind[58]+alpha[8]*fUpwind[57]+alpha[7]*fUpwind[56]+alpha[10]*fUpwind[55]+alpha[9]*fUpwind[54]+alpha[26]*fUpwind[53])+15.0*(alpha[3]*fUpwind[52]+alpha[4]*fUpwind[51]+alpha[19]*fUpwind[50]+alpha[18]*fUpwind[49])+15.0*alpha[11]*fUpwind[48]+13.41640786499874*(alpha[20]*fUpwind[31]+alpha[13]*fUpwind[30]+alpha[12]*fUpwind[29]+alpha[27]*fUpwind[28]+alpha[5]*fUpwind[25]+alpha[22]*fUpwind[24]+alpha[21]*fUpwind[23]+alpha[14]*fUpwind[15])); 
  Ghat[59] = 0.01178511301977579*(15.0*alpha[4]*fUpwind[63]+15.0*(alpha[9]*fUpwind[62]+alpha[10]*fUpwind[61]+alpha[11]*fUpwind[60]+alpha[0]*fUpwind[59])+15.0*(alpha[17]*fUpwind[58]+alpha[18]*fUpwind[57]+alpha[19]*fUpwind[56]+alpha[1]*fUpwind[55]+alpha[2]*fUpwind[54]+alpha[3]*fUpwind[53])+15.0*(alpha[26]*fUpwind[52]+alpha[6]*fUpwind[51]+alpha[7]*fUpwind[50]+alpha[8]*fUpwind[49])+15.0*alpha[16]*fUpwind[48]+13.41640786499874*(alpha[5]*fUpwind[27]+fUpwind[5]*alpha[27]+alpha[12]*fUpwind[22]+fUpwind[12]*alpha[22]+alpha[13]*fUpwind[21]+fUpwind[13]*alpha[21]+alpha[14]*fUpwind[20]+fUpwind[14]*alpha[20])); 
  Ghat[60] = 0.01178511301977579*(15.0*alpha[3]*fUpwind[63]+15.0*(alpha[7]*fUpwind[62]+alpha[8]*fUpwind[61]+alpha[0]*fUpwind[60]+alpha[11]*fUpwind[59])+15.0*(alpha[16]*fUpwind[58]+alpha[1]*fUpwind[57]+alpha[2]*fUpwind[56]+alpha[18]*fUpwind[55]+alpha[19]*fUpwind[54]+alpha[4]*fUpwind[53])+15.0*(alpha[6]*fUpwind[52]+alpha[26]*fUpwind[51]+alpha[9]*fUpwind[50]+alpha[10]*fUpwind[49])+15.0*alpha[17]*fUpwind[48]+13.41640786499874*(alpha[14]*fUpwind[31]+alpha[21]*fUpwind[30]+alpha[22]*fUpwind[29]+alpha[5]*fUpwind[28]+fUpwind[25]*alpha[27]+alpha[12]*fUpwind[24]+alpha[13]*fUpwind[23]+fUpwind[15]*alpha[20])); 
  Ghat[61] = 0.01178511301977579*(15.0*alpha[2]*fUpwind[63]+15.0*(alpha[6]*fUpwind[62]+alpha[0]*fUpwind[61]+alpha[8]*fUpwind[60]+alpha[10]*fUpwind[59])+15.0*(alpha[1]*fUpwind[58]+alpha[16]*fUpwind[57]+alpha[3]*fUpwind[56]+alpha[17]*fUpwind[55]+alpha[4]*fUpwind[54]+alpha[19]*fUpwind[53])+15.0*(alpha[7]*fUpwind[52]+alpha[9]*fUpwind[51]+alpha[26]*fUpwind[50]+alpha[11]*fUpwind[49])+15.0*alpha[18]*fUpwind[48]+13.41640786499874*(alpha[13]*fUpwind[31]+alpha[20]*fUpwind[30]+alpha[5]*fUpwind[29]+alpha[22]*fUpwind[28]+fUpwind[24]*alpha[27]+alpha[12]*fUpwind[25]+alpha[14]*fUpwind[23]+fUpwind[15]*alpha[21])); 
  Ghat[62] = 0.01178511301977579*(15.0*alpha[1]*fUpwind[63]+15.0*(alpha[0]*fUpwind[62]+alpha[6]*fUpwind[61]+alpha[7]*fUpwind[60]+alpha[9]*fUpwind[59])+15.0*(alpha[2]*fUpwind[58]+alpha[3]*fUpwind[57]+alpha[16]*fUpwind[56]+alpha[4]*fUpwind[55]+alpha[17]*fUpwind[54]+alpha[18]*fUpwind[53])+15.0*(alpha[8]*fUpwind[52]+alpha[10]*fUpwind[51]+alpha[11]*fUpwind[50]+alpha[26]*fUpwind[49])+15.0*alpha[19]*fUpwind[48]+13.41640786499874*(alpha[12]*fUpwind[31]+alpha[5]*fUpwind[30]+alpha[20]*fUpwind[29]+alpha[21]*fUpwind[28]+fUpwind[23]*alpha[27]+alpha[13]*fUpwind[25]+alpha[14]*fUpwind[24]+fUpwind[15]*alpha[22])); 
  Ghat[63] = 0.01178511301977579*(15.0*alpha[0]*fUpwind[63]+15.0*(alpha[1]*fUpwind[62]+alpha[2]*fUpwind[61]+alpha[3]*fUpwind[60]+alpha[4]*fUpwind[59])+15.0*(alpha[6]*fUpwind[58]+alpha[7]*fUpwind[57]+alpha[8]*fUpwind[56]+alpha[9]*fUpwind[55]+alpha[10]*fUpwind[54]+alpha[11]*fUpwind[53])+15.0*(alpha[16]*fUpwind[52]+alpha[17]*fUpwind[51]+alpha[18]*fUpwind[50]+alpha[19]*fUpwind[49])+15.0*alpha[26]*fUpwind[48]+13.41640786499874*(alpha[5]*fUpwind[31]+alpha[12]*fUpwind[30]+alpha[13]*fUpwind[29]+alpha[14]*fUpwind[28]+fUpwind[15]*alpha[27]+alpha[20]*fUpwind[25]+alpha[21]*fUpwind[24]+alpha[22]*fUpwind[23])); 

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
  out[64] += 0.7071067811865475*Ghat[32]*dv11; 
  out[65] += 0.7071067811865475*Ghat[33]*dv11; 
  out[66] += 0.7071067811865475*Ghat[34]*dv11; 
  out[67] += 0.7071067811865475*Ghat[35]*dv11; 
  out[68] += -1.224744871391589*Ghat[32]*dv11; 
  out[69] += 0.7071067811865475*Ghat[36]*dv11; 
  out[70] += 0.7071067811865475*Ghat[37]*dv11; 
  out[71] += 0.7071067811865475*Ghat[38]*dv11; 
  out[72] += 0.7071067811865475*Ghat[39]*dv11; 
  out[73] += -1.224744871391589*Ghat[33]*dv11; 
  out[74] += -1.224744871391589*Ghat[34]*dv11; 
  out[75] += -1.224744871391589*Ghat[35]*dv11; 
  out[76] += 0.7071067811865475*Ghat[40]*dv11; 
  out[77] += 0.7071067811865475*Ghat[41]*dv11; 
  out[78] += 0.7071067811865475*Ghat[42]*dv11; 
  out[79] += -1.224744871391589*Ghat[36]*dv11; 
  out[80] += 0.7071067811865475*Ghat[43]*dv11; 
  out[81] += -1.224744871391589*Ghat[37]*dv11; 
  out[82] += -1.224744871391589*Ghat[38]*dv11; 
  out[83] += -1.224744871391589*Ghat[39]*dv11; 
  out[84] += 0.7071067811865475*Ghat[44]*dv11; 
  out[85] += 0.7071067811865475*Ghat[45]*dv11; 
  out[86] += 0.7071067811865475*Ghat[46]*dv11; 
  out[87] += -1.224744871391589*Ghat[40]*dv11; 
  out[88] += -1.224744871391589*Ghat[41]*dv11; 
  out[89] += -1.224744871391589*Ghat[42]*dv11; 
  out[90] += -1.224744871391589*Ghat[43]*dv11; 
  out[91] += 0.7071067811865475*Ghat[47]*dv11; 
  out[92] += -1.224744871391589*Ghat[44]*dv11; 
  out[93] += -1.224744871391589*Ghat[45]*dv11; 
  out[94] += -1.224744871391589*Ghat[46]*dv11; 
  out[95] += -1.224744871391589*Ghat[47]*dv11; 
  out[96] += 1.58113883008419*Ghat[0]*dv11; 
  out[97] += 1.58113883008419*Ghat[1]*dv11; 
  out[98] += 1.58113883008419*Ghat[2]*dv11; 
  out[99] += 1.58113883008419*Ghat[3]*dv11; 
  out[100] += 1.58113883008419*Ghat[4]*dv11; 
  out[101] += 1.58113883008419*Ghat[5]*dv11; 
  out[102] += 1.58113883008419*Ghat[6]*dv11; 
  out[103] += 1.58113883008419*Ghat[7]*dv11; 
  out[104] += 1.58113883008419*Ghat[8]*dv11; 
  out[105] += 1.58113883008419*Ghat[9]*dv11; 
  out[106] += 1.58113883008419*Ghat[10]*dv11; 
  out[107] += 1.58113883008419*Ghat[11]*dv11; 
  out[108] += 1.58113883008419*Ghat[12]*dv11; 
  out[109] += 1.58113883008419*Ghat[13]*dv11; 
  out[110] += 1.58113883008419*Ghat[14]*dv11; 
  out[111] += 1.58113883008419*Ghat[15]*dv11; 
  out[112] += 1.58113883008419*Ghat[16]*dv11; 
  out[113] += 1.58113883008419*Ghat[17]*dv11; 
  out[114] += 1.58113883008419*Ghat[18]*dv11; 
  out[115] += 1.58113883008419*Ghat[19]*dv11; 
  out[116] += 1.58113883008419*Ghat[20]*dv11; 
  out[117] += 1.58113883008419*Ghat[21]*dv11; 
  out[118] += 1.58113883008419*Ghat[22]*dv11; 
  out[119] += 1.58113883008419*Ghat[23]*dv11; 
  out[120] += 1.58113883008419*Ghat[24]*dv11; 
  out[121] += 1.58113883008419*Ghat[25]*dv11; 
  out[122] += 1.58113883008419*Ghat[26]*dv11; 
  out[123] += 1.58113883008419*Ghat[27]*dv11; 
  out[124] += 1.58113883008419*Ghat[28]*dv11; 
  out[125] += 1.58113883008419*Ghat[29]*dv11; 
  out[126] += 1.58113883008419*Ghat[30]*dv11; 
  out[127] += 1.58113883008419*Ghat[31]*dv11; 
  out[128] += 0.7071067811865475*Ghat[48]*dv11; 
  out[129] += 0.7071067811865475*Ghat[49]*dv11; 
  out[130] += 0.7071067811865475*Ghat[50]*dv11; 
  out[131] += 0.7071067811865475*Ghat[51]*dv11; 
  out[132] += 0.7071067811865475*Ghat[52]*dv11; 
  out[133] += -1.224744871391589*Ghat[48]*dv11; 
  out[134] += 0.7071067811865475*Ghat[53]*dv11; 
  out[135] += 0.7071067811865475*Ghat[54]*dv11; 
  out[136] += 0.7071067811865475*Ghat[55]*dv11; 
  out[137] += 0.7071067811865475*Ghat[56]*dv11; 
  out[138] += 0.7071067811865475*Ghat[57]*dv11; 
  out[139] += 0.7071067811865475*Ghat[58]*dv11; 
  out[140] += -1.224744871391589*Ghat[49]*dv11; 
  out[141] += -1.224744871391589*Ghat[50]*dv11; 
  out[142] += -1.224744871391589*Ghat[51]*dv11; 
  out[143] += -1.224744871391589*Ghat[52]*dv11; 
  out[144] += 0.7071067811865475*Ghat[59]*dv11; 
  out[145] += 0.7071067811865475*Ghat[60]*dv11; 
  out[146] += 0.7071067811865475*Ghat[61]*dv11; 
  out[147] += 0.7071067811865475*Ghat[62]*dv11; 
  out[148] += -1.224744871391589*Ghat[53]*dv11; 
  out[149] += -1.224744871391589*Ghat[54]*dv11; 
  out[150] += -1.224744871391589*Ghat[55]*dv11; 
  out[151] += -1.224744871391589*Ghat[56]*dv11; 
  out[152] += -1.224744871391589*Ghat[57]*dv11; 
  out[153] += -1.224744871391589*Ghat[58]*dv11; 
  out[154] += 0.7071067811865475*Ghat[63]*dv11; 
  out[155] += -1.224744871391589*Ghat[59]*dv11; 
  out[156] += -1.224744871391589*Ghat[60]*dv11; 
  out[157] += -1.224744871391589*Ghat[61]*dv11; 
  out[158] += -1.224744871391589*Ghat[62]*dv11; 
  out[159] += -1.224744871391589*Ghat[63]*dv11; 

  } 
} 
