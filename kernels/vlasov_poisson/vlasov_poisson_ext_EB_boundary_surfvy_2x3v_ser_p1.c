#include <gkyl_vlasov_poisson_kernels.h> 
#include <gkyl_basis_hyb_2x3v_p1_surfx4_eval_quad.h> 
#include <gkyl_basis_hyb_2x3v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double vlasov_poisson_ext_EB_boundary_surfvy_2x3v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w:           Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // pots:        potentials phi_tot=phi+phi_ext and A_ext (scaled by q/m).
  // EBext:       external E and B fields (scaled by q/m).
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Output distribution function in skin cell 

  const double dv11 = 2/dxv[3]; 
  const double dx10 = 2/dxv[0]; 
  const double dx11 = 2/dxv[1]; 
  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv2 = dxv[3], wv2 = w[3]; 
  const double dv3 = dxv[4], wv3 = w[4]; 

  double alpha[32] = {0.0}; 

  const double *phi = &pots[0]; 

  const double *Ex = &EBext[0]; 
  const double *Ey = &EBext[4]; 
  const double *Ez = &EBext[8]; 
  const double *Bx = &EBext[12]; 
  const double *By = &EBext[16]; 
  const double *Bz = &EBext[20]; 

  alpha[0] = 2.0*Bx[0]*wv3-2.0*Bz[0]*wv1-3.4641016151377544*phi[2]*dx11+2.0*Ey[0]; 
  alpha[1] = 2.0*Bx[1]*wv3-2.0*Bz[1]*wv1-3.4641016151377544*phi[3]*dx11+2.0*Ey[1]; 
  alpha[2] = 2.0*Bx[2]*wv3-2.0*Bz[2]*wv1+2.0*Ey[2]; 
  alpha[3] = -(0.5773502691896258*Bz[0]*dv1); 
  alpha[4] = 0.5773502691896258*Bx[0]*dv3; 
  alpha[5] = 2.0*Bx[3]*wv3-2.0*Bz[3]*wv1+2.0*Ey[3]; 
  alpha[6] = -(0.5773502691896258*Bz[1]*dv1); 
  alpha[7] = -(0.5773502691896258*Bz[2]*dv1); 
  alpha[8] = 0.5773502691896258*Bx[1]*dv3; 
  alpha[9] = 0.5773502691896258*Bx[2]*dv3; 
  alpha[11] = -(0.5773502691896258*Bz[3]*dv1); 
  alpha[12] = 0.5773502691896258*Bx[3]*dv3; 

  double fUpwindQuad[36] = {0.0};
  double fUpwind[32] = {0.0};
  double Ghat[32] = {0.0}; 

  if (edge == -1) { 

  if (-(0.33541019662496846*(alpha[12]+alpha[11]))+0.33541019662496846*(alpha[9]+alpha[8]+alpha[7]+alpha[6])+0.25*alpha[5]-0.33541019662496846*(alpha[4]+alpha[3])-0.25*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[0] = hyb_2x3v_p1_surfx4_eval_quad_node_0_r(fSkin); 
  } else { 
    fUpwindQuad[0] = hyb_2x3v_p1_surfx4_eval_quad_node_0_l(fEdge); 
  } 
  if (-(0.33541019662496846*alpha[11])+0.33541019662496846*(alpha[7]+alpha[6])+0.25*alpha[5]-0.33541019662496846*alpha[3]-0.25*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[1] = hyb_2x3v_p1_surfx4_eval_quad_node_1_r(fSkin); 
  } else { 
    fUpwindQuad[1] = hyb_2x3v_p1_surfx4_eval_quad_node_1_l(fEdge); 
  } 
  if (0.33541019662496846*alpha[12]-0.33541019662496846*(alpha[11]+alpha[9]+alpha[8])+0.33541019662496846*(alpha[7]+alpha[6])+0.25*alpha[5]+0.33541019662496846*alpha[4]-0.33541019662496846*alpha[3]-0.25*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[2] = hyb_2x3v_p1_surfx4_eval_quad_node_2_r(fSkin); 
  } else { 
    fUpwindQuad[2] = hyb_2x3v_p1_surfx4_eval_quad_node_2_l(fEdge); 
  } 
  if (-(0.33541019662496846*alpha[12])+0.33541019662496846*(alpha[9]+alpha[8])+0.25*alpha[5]-0.33541019662496846*alpha[4]-0.25*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[3] = hyb_2x3v_p1_surfx4_eval_quad_node_3_r(fSkin); 
  } else { 
    fUpwindQuad[3] = hyb_2x3v_p1_surfx4_eval_quad_node_3_l(fEdge); 
  } 
  if (0.25*alpha[5]-0.25*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[4] = hyb_2x3v_p1_surfx4_eval_quad_node_4_r(fSkin); 
  } else { 
    fUpwindQuad[4] = hyb_2x3v_p1_surfx4_eval_quad_node_4_l(fEdge); 
  } 
  if (0.33541019662496846*alpha[12]-0.33541019662496846*(alpha[9]+alpha[8])+0.25*alpha[5]+0.33541019662496846*alpha[4]-0.25*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[5] = hyb_2x3v_p1_surfx4_eval_quad_node_5_r(fSkin); 
  } else { 
    fUpwindQuad[5] = hyb_2x3v_p1_surfx4_eval_quad_node_5_l(fEdge); 
  } 
  if (-(0.33541019662496846*alpha[12])+0.33541019662496846*(alpha[11]+alpha[9]+alpha[8])-0.33541019662496846*(alpha[7]+alpha[6])+0.25*alpha[5]-0.33541019662496846*alpha[4]+0.33541019662496846*alpha[3]-0.25*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[6] = hyb_2x3v_p1_surfx4_eval_quad_node_6_r(fSkin); 
  } else { 
    fUpwindQuad[6] = hyb_2x3v_p1_surfx4_eval_quad_node_6_l(fEdge); 
  } 
  if (0.33541019662496846*alpha[11]-0.33541019662496846*(alpha[7]+alpha[6])+0.25*alpha[5]+0.33541019662496846*alpha[3]-0.25*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[7] = hyb_2x3v_p1_surfx4_eval_quad_node_7_r(fSkin); 
  } else { 
    fUpwindQuad[7] = hyb_2x3v_p1_surfx4_eval_quad_node_7_l(fEdge); 
  } 
  if (0.33541019662496846*(alpha[12]+alpha[11])-0.33541019662496846*(alpha[9]+alpha[8]+alpha[7]+alpha[6])+0.25*alpha[5]+0.33541019662496846*(alpha[4]+alpha[3])-0.25*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[8] = hyb_2x3v_p1_surfx4_eval_quad_node_8_r(fSkin); 
  } else { 
    fUpwindQuad[8] = hyb_2x3v_p1_surfx4_eval_quad_node_8_l(fEdge); 
  } 
  if (0.33541019662496846*(alpha[12]+alpha[11])-0.33541019662496846*alpha[9]+0.33541019662496846*alpha[8]-0.33541019662496846*alpha[7]+0.33541019662496846*alpha[6]-0.25*alpha[5]-0.33541019662496846*(alpha[4]+alpha[3])+0.25*alpha[2]-0.25*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[9] = hyb_2x3v_p1_surfx4_eval_quad_node_9_r(fSkin); 
  } else { 
    fUpwindQuad[9] = hyb_2x3v_p1_surfx4_eval_quad_node_9_l(fEdge); 
  } 
  if (0.33541019662496846*alpha[11]-0.33541019662496846*alpha[7]+0.33541019662496846*alpha[6]-0.25*alpha[5]-0.33541019662496846*alpha[3]+0.25*alpha[2]-0.25*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[10] = hyb_2x3v_p1_surfx4_eval_quad_node_10_r(fSkin); 
  } else { 
    fUpwindQuad[10] = hyb_2x3v_p1_surfx4_eval_quad_node_10_l(fEdge); 
  } 
  if (-(0.33541019662496846*alpha[12])+0.33541019662496846*(alpha[11]+alpha[9])-0.33541019662496846*(alpha[8]+alpha[7])+0.33541019662496846*alpha[6]-0.25*alpha[5]+0.33541019662496846*alpha[4]-0.33541019662496846*alpha[3]+0.25*alpha[2]-0.25*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[11] = hyb_2x3v_p1_surfx4_eval_quad_node_11_r(fSkin); 
  } else { 
    fUpwindQuad[11] = hyb_2x3v_p1_surfx4_eval_quad_node_11_l(fEdge); 
  } 
  if (0.33541019662496846*alpha[12]-0.33541019662496846*alpha[9]+0.33541019662496846*alpha[8]-0.25*alpha[5]-0.33541019662496846*alpha[4]+0.25*alpha[2]-0.25*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[12] = hyb_2x3v_p1_surfx4_eval_quad_node_12_r(fSkin); 
  } else { 
    fUpwindQuad[12] = hyb_2x3v_p1_surfx4_eval_quad_node_12_l(fEdge); 
  } 
  if (-(0.25*alpha[5])+0.25*alpha[2]-0.25*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[13] = hyb_2x3v_p1_surfx4_eval_quad_node_13_r(fSkin); 
  } else { 
    fUpwindQuad[13] = hyb_2x3v_p1_surfx4_eval_quad_node_13_l(fEdge); 
  } 
  if (-(0.33541019662496846*alpha[12])+0.33541019662496846*alpha[9]-0.33541019662496846*alpha[8]-0.25*alpha[5]+0.33541019662496846*alpha[4]+0.25*alpha[2]-0.25*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[14] = hyb_2x3v_p1_surfx4_eval_quad_node_14_r(fSkin); 
  } else { 
    fUpwindQuad[14] = hyb_2x3v_p1_surfx4_eval_quad_node_14_l(fEdge); 
  } 
  if (0.33541019662496846*alpha[12]-0.33541019662496846*(alpha[11]+alpha[9])+0.33541019662496846*(alpha[8]+alpha[7])-0.33541019662496846*alpha[6]-0.25*alpha[5]-0.33541019662496846*alpha[4]+0.33541019662496846*alpha[3]+0.25*alpha[2]-0.25*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[15] = hyb_2x3v_p1_surfx4_eval_quad_node_15_r(fSkin); 
  } else { 
    fUpwindQuad[15] = hyb_2x3v_p1_surfx4_eval_quad_node_15_l(fEdge); 
  } 
  if (-(0.33541019662496846*alpha[11])+0.33541019662496846*alpha[7]-0.33541019662496846*alpha[6]-0.25*alpha[5]+0.33541019662496846*alpha[3]+0.25*alpha[2]-0.25*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[16] = hyb_2x3v_p1_surfx4_eval_quad_node_16_r(fSkin); 
  } else { 
    fUpwindQuad[16] = hyb_2x3v_p1_surfx4_eval_quad_node_16_l(fEdge); 
  } 
  if (-(0.33541019662496846*(alpha[12]+alpha[11]))+0.33541019662496846*alpha[9]-0.33541019662496846*alpha[8]+0.33541019662496846*alpha[7]-0.33541019662496846*alpha[6]-0.25*alpha[5]+0.33541019662496846*(alpha[4]+alpha[3])+0.25*alpha[2]-0.25*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[17] = hyb_2x3v_p1_surfx4_eval_quad_node_17_r(fSkin); 
  } else { 
    fUpwindQuad[17] = hyb_2x3v_p1_surfx4_eval_quad_node_17_l(fEdge); 
  } 
  if (0.33541019662496846*(alpha[12]+alpha[11]+alpha[9])-0.33541019662496846*alpha[8]+0.33541019662496846*alpha[7]-0.33541019662496846*alpha[6]-0.25*alpha[5]-0.33541019662496846*(alpha[4]+alpha[3])-0.25*alpha[2]+0.25*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[18] = hyb_2x3v_p1_surfx4_eval_quad_node_18_r(fSkin); 
  } else { 
    fUpwindQuad[18] = hyb_2x3v_p1_surfx4_eval_quad_node_18_l(fEdge); 
  } 
  if (0.33541019662496846*(alpha[11]+alpha[7])-0.33541019662496846*alpha[6]-0.25*alpha[5]-0.33541019662496846*alpha[3]-0.25*alpha[2]+0.25*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[19] = hyb_2x3v_p1_surfx4_eval_quad_node_19_r(fSkin); 
  } else { 
    fUpwindQuad[19] = hyb_2x3v_p1_surfx4_eval_quad_node_19_l(fEdge); 
  } 
  if (-(0.33541019662496846*alpha[12])+0.33541019662496846*alpha[11]-0.33541019662496846*alpha[9]+0.33541019662496846*(alpha[8]+alpha[7])-0.33541019662496846*alpha[6]-0.25*alpha[5]+0.33541019662496846*alpha[4]-0.33541019662496846*alpha[3]-0.25*alpha[2]+0.25*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[20] = hyb_2x3v_p1_surfx4_eval_quad_node_20_r(fSkin); 
  } else { 
    fUpwindQuad[20] = hyb_2x3v_p1_surfx4_eval_quad_node_20_l(fEdge); 
  } 
  if (0.33541019662496846*(alpha[12]+alpha[9])-0.33541019662496846*alpha[8]-0.25*alpha[5]-0.33541019662496846*alpha[4]-0.25*alpha[2]+0.25*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[21] = hyb_2x3v_p1_surfx4_eval_quad_node_21_r(fSkin); 
  } else { 
    fUpwindQuad[21] = hyb_2x3v_p1_surfx4_eval_quad_node_21_l(fEdge); 
  } 
  if (0.25*(alpha[1]+alpha[0])-0.25*(alpha[5]+alpha[2]) > 0) { 
    fUpwindQuad[22] = hyb_2x3v_p1_surfx4_eval_quad_node_22_r(fSkin); 
  } else { 
    fUpwindQuad[22] = hyb_2x3v_p1_surfx4_eval_quad_node_22_l(fEdge); 
  } 
  if (-(0.33541019662496846*(alpha[12]+alpha[9]))+0.33541019662496846*alpha[8]-0.25*alpha[5]+0.33541019662496846*alpha[4]-0.25*alpha[2]+0.25*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[23] = hyb_2x3v_p1_surfx4_eval_quad_node_23_r(fSkin); 
  } else { 
    fUpwindQuad[23] = hyb_2x3v_p1_surfx4_eval_quad_node_23_l(fEdge); 
  } 
  if (0.33541019662496846*alpha[12]-0.33541019662496846*alpha[11]+0.33541019662496846*alpha[9]-0.33541019662496846*(alpha[8]+alpha[7])+0.33541019662496846*alpha[6]-0.25*alpha[5]-0.33541019662496846*alpha[4]+0.33541019662496846*alpha[3]-0.25*alpha[2]+0.25*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[24] = hyb_2x3v_p1_surfx4_eval_quad_node_24_r(fSkin); 
  } else { 
    fUpwindQuad[24] = hyb_2x3v_p1_surfx4_eval_quad_node_24_l(fEdge); 
  } 
  if (-(0.33541019662496846*(alpha[11]+alpha[7]))+0.33541019662496846*alpha[6]-0.25*alpha[5]+0.33541019662496846*alpha[3]-0.25*alpha[2]+0.25*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[25] = hyb_2x3v_p1_surfx4_eval_quad_node_25_r(fSkin); 
  } else { 
    fUpwindQuad[25] = hyb_2x3v_p1_surfx4_eval_quad_node_25_l(fEdge); 
  } 
  if (-(0.33541019662496846*(alpha[12]+alpha[11]+alpha[9]))+0.33541019662496846*alpha[8]-0.33541019662496846*alpha[7]+0.33541019662496846*alpha[6]-0.25*alpha[5]+0.33541019662496846*(alpha[4]+alpha[3])-0.25*alpha[2]+0.25*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[26] = hyb_2x3v_p1_surfx4_eval_quad_node_26_r(fSkin); 
  } else { 
    fUpwindQuad[26] = hyb_2x3v_p1_surfx4_eval_quad_node_26_l(fEdge); 
  } 
  if (-(0.33541019662496846*(alpha[12]+alpha[11]+alpha[9]+alpha[8]+alpha[7]+alpha[6]))+0.25*alpha[5]-0.33541019662496846*(alpha[4]+alpha[3])+0.25*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[27] = hyb_2x3v_p1_surfx4_eval_quad_node_27_r(fSkin); 
  } else { 
    fUpwindQuad[27] = hyb_2x3v_p1_surfx4_eval_quad_node_27_l(fEdge); 
  } 
  if (-(0.33541019662496846*(alpha[11]+alpha[7]+alpha[6]))+0.25*alpha[5]-0.33541019662496846*alpha[3]+0.25*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[28] = hyb_2x3v_p1_surfx4_eval_quad_node_28_r(fSkin); 
  } else { 
    fUpwindQuad[28] = hyb_2x3v_p1_surfx4_eval_quad_node_28_l(fEdge); 
  } 
  if (0.33541019662496846*alpha[12]-0.33541019662496846*alpha[11]+0.33541019662496846*(alpha[9]+alpha[8])-0.33541019662496846*(alpha[7]+alpha[6])+0.25*alpha[5]+0.33541019662496846*alpha[4]-0.33541019662496846*alpha[3]+0.25*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[29] = hyb_2x3v_p1_surfx4_eval_quad_node_29_r(fSkin); 
  } else { 
    fUpwindQuad[29] = hyb_2x3v_p1_surfx4_eval_quad_node_29_l(fEdge); 
  } 
  if (-(0.33541019662496846*(alpha[12]+alpha[9]+alpha[8]))+0.25*alpha[5]-0.33541019662496846*alpha[4]+0.25*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[30] = hyb_2x3v_p1_surfx4_eval_quad_node_30_r(fSkin); 
  } else { 
    fUpwindQuad[30] = hyb_2x3v_p1_surfx4_eval_quad_node_30_l(fEdge); 
  } 
  if (0.25*(alpha[5]+alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[31] = hyb_2x3v_p1_surfx4_eval_quad_node_31_r(fSkin); 
  } else { 
    fUpwindQuad[31] = hyb_2x3v_p1_surfx4_eval_quad_node_31_l(fEdge); 
  } 
  if (0.33541019662496846*(alpha[12]+alpha[9]+alpha[8])+0.25*alpha[5]+0.33541019662496846*alpha[4]+0.25*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[32] = hyb_2x3v_p1_surfx4_eval_quad_node_32_r(fSkin); 
  } else { 
    fUpwindQuad[32] = hyb_2x3v_p1_surfx4_eval_quad_node_32_l(fEdge); 
  } 
  if (-(0.33541019662496846*alpha[12])+0.33541019662496846*alpha[11]-0.33541019662496846*(alpha[9]+alpha[8])+0.33541019662496846*(alpha[7]+alpha[6])+0.25*alpha[5]-0.33541019662496846*alpha[4]+0.33541019662496846*alpha[3]+0.25*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[33] = hyb_2x3v_p1_surfx4_eval_quad_node_33_r(fSkin); 
  } else { 
    fUpwindQuad[33] = hyb_2x3v_p1_surfx4_eval_quad_node_33_l(fEdge); 
  } 
  if (0.33541019662496846*(alpha[11]+alpha[7]+alpha[6])+0.25*alpha[5]+0.33541019662496846*alpha[3]+0.25*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[34] = hyb_2x3v_p1_surfx4_eval_quad_node_34_r(fSkin); 
  } else { 
    fUpwindQuad[34] = hyb_2x3v_p1_surfx4_eval_quad_node_34_l(fEdge); 
  } 
  if (0.33541019662496846*(alpha[12]+alpha[11]+alpha[9]+alpha[8]+alpha[7]+alpha[6])+0.25*alpha[5]+0.33541019662496846*(alpha[4]+alpha[3])+0.25*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[35] = hyb_2x3v_p1_surfx4_eval_quad_node_35_r(fSkin); 
  } else { 
    fUpwindQuad[35] = hyb_2x3v_p1_surfx4_eval_quad_node_35_l(fEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_2x3v_p1_vdir_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.25*alpha[12]*fUpwind[12]+0.25*alpha[11]*fUpwind[11]+0.25*alpha[9]*fUpwind[9]+0.25*alpha[8]*fUpwind[8]+0.25*alpha[7]*fUpwind[7]+0.25*alpha[6]*fUpwind[6]+0.25*alpha[5]*fUpwind[5]+0.25*alpha[4]*fUpwind[4]+0.25*alpha[3]*fUpwind[3]+0.25*alpha[2]*fUpwind[2]+0.25*alpha[1]*fUpwind[1]+0.25*alpha[0]*fUpwind[0]; 
  Ghat[1] = 0.25*alpha[9]*fUpwind[12]+0.25*fUpwind[9]*alpha[12]+0.25*alpha[7]*fUpwind[11]+0.25*fUpwind[7]*alpha[11]+0.25*alpha[4]*fUpwind[8]+0.25*fUpwind[4]*alpha[8]+0.25*alpha[3]*fUpwind[6]+0.25*fUpwind[3]*alpha[6]+0.25*alpha[2]*fUpwind[5]+0.25*fUpwind[2]*alpha[5]+0.25*alpha[0]*fUpwind[1]+0.25*fUpwind[0]*alpha[1]; 
  Ghat[2] = 0.25*alpha[8]*fUpwind[12]+0.25*fUpwind[8]*alpha[12]+0.25*alpha[6]*fUpwind[11]+0.25*fUpwind[6]*alpha[11]+0.25*alpha[4]*fUpwind[9]+0.25*fUpwind[4]*alpha[9]+0.25*alpha[3]*fUpwind[7]+0.25*fUpwind[3]*alpha[7]+0.25*alpha[1]*fUpwind[5]+0.25*fUpwind[1]*alpha[5]+0.25*alpha[0]*fUpwind[2]+0.25*fUpwind[0]*alpha[2]; 
  Ghat[3] = 0.22360679774997896*alpha[11]*fUpwind[20]+0.22360679774997902*alpha[7]*fUpwind[18]+0.22360679774997902*alpha[6]*fUpwind[17]+0.22360679774997896*alpha[3]*fUpwind[16]+0.25*alpha[12]*fUpwind[15]+0.25*alpha[9]*fUpwind[14]+0.25*alpha[8]*fUpwind[13]+0.25*alpha[5]*fUpwind[11]+0.25*fUpwind[5]*alpha[11]+0.25*alpha[4]*fUpwind[10]+0.25*alpha[2]*fUpwind[7]+0.25*fUpwind[2]*alpha[7]+0.25*alpha[1]*fUpwind[6]+0.25*fUpwind[1]*alpha[6]+0.25*alpha[0]*fUpwind[3]+0.25*fUpwind[0]*alpha[3]; 
  Ghat[4] = 0.22360679774997896*alpha[12]*fUpwind[28]+0.22360679774997902*alpha[9]*fUpwind[26]+0.22360679774997902*alpha[8]*fUpwind[25]+0.22360679774997896*alpha[4]*fUpwind[24]+0.25*alpha[11]*fUpwind[15]+0.25*alpha[7]*fUpwind[14]+0.25*alpha[6]*fUpwind[13]+0.25*alpha[5]*fUpwind[12]+0.25*fUpwind[5]*alpha[12]+0.25*alpha[3]*fUpwind[10]+0.25*alpha[2]*fUpwind[9]+0.25*fUpwind[2]*alpha[9]+0.25*alpha[1]*fUpwind[8]+0.25*fUpwind[1]*alpha[8]+0.25*alpha[0]*fUpwind[4]+0.25*fUpwind[0]*alpha[4]; 
  Ghat[5] = 0.25*alpha[4]*fUpwind[12]+0.25*fUpwind[4]*alpha[12]+0.25*alpha[3]*fUpwind[11]+0.25*fUpwind[3]*alpha[11]+0.25*alpha[8]*fUpwind[9]+0.25*fUpwind[8]*alpha[9]+0.25*alpha[6]*fUpwind[7]+0.25*fUpwind[6]*alpha[7]+0.25*alpha[0]*fUpwind[5]+0.25*fUpwind[0]*alpha[5]+0.25*alpha[1]*fUpwind[2]+0.25*fUpwind[1]*alpha[2]; 
  Ghat[6] = 0.22360679774997896*alpha[7]*fUpwind[20]+0.22360679774997902*alpha[11]*fUpwind[18]+0.22360679774997902*alpha[3]*fUpwind[17]+0.22360679774997896*alpha[6]*fUpwind[16]+0.25*alpha[9]*fUpwind[15]+0.25*alpha[12]*fUpwind[14]+0.25*alpha[4]*fUpwind[13]+0.25*alpha[2]*fUpwind[11]+0.25*fUpwind[2]*alpha[11]+0.25*alpha[8]*fUpwind[10]+0.25*alpha[5]*fUpwind[7]+0.25*fUpwind[5]*alpha[7]+0.25*alpha[0]*fUpwind[6]+0.25*fUpwind[0]*alpha[6]+0.25*alpha[1]*fUpwind[3]+0.25*fUpwind[1]*alpha[3]; 
  Ghat[7] = 0.22360679774997896*alpha[6]*fUpwind[20]+0.22360679774997902*alpha[3]*fUpwind[18]+0.22360679774997902*alpha[11]*fUpwind[17]+0.22360679774997896*alpha[7]*fUpwind[16]+0.25*alpha[8]*fUpwind[15]+0.25*alpha[4]*fUpwind[14]+0.25*alpha[12]*fUpwind[13]+0.25*alpha[1]*fUpwind[11]+0.25*fUpwind[1]*alpha[11]+0.25*alpha[9]*fUpwind[10]+0.25*alpha[0]*fUpwind[7]+0.25*fUpwind[0]*alpha[7]+0.25*alpha[5]*fUpwind[6]+0.25*fUpwind[5]*alpha[6]+0.25*alpha[2]*fUpwind[3]+0.25*fUpwind[2]*alpha[3]; 
  Ghat[8] = 0.22360679774997896*alpha[9]*fUpwind[28]+0.22360679774997902*alpha[12]*fUpwind[26]+0.22360679774997902*alpha[4]*fUpwind[25]+0.22360679774997896*alpha[8]*fUpwind[24]+0.25*alpha[7]*fUpwind[15]+0.25*alpha[11]*fUpwind[14]+0.25*alpha[3]*fUpwind[13]+0.25*alpha[2]*fUpwind[12]+0.25*fUpwind[2]*alpha[12]+0.25*alpha[6]*fUpwind[10]+0.25*alpha[5]*fUpwind[9]+0.25*fUpwind[5]*alpha[9]+0.25*alpha[0]*fUpwind[8]+0.25*fUpwind[0]*alpha[8]+0.25*alpha[1]*fUpwind[4]+0.25*fUpwind[1]*alpha[4]; 
  Ghat[9] = 0.22360679774997896*alpha[8]*fUpwind[28]+0.22360679774997902*alpha[4]*fUpwind[26]+0.22360679774997902*alpha[12]*fUpwind[25]+0.22360679774997896*alpha[9]*fUpwind[24]+0.25*alpha[6]*fUpwind[15]+0.25*alpha[3]*fUpwind[14]+0.25*alpha[11]*fUpwind[13]+0.25*alpha[1]*fUpwind[12]+0.25*fUpwind[1]*alpha[12]+0.25*alpha[7]*fUpwind[10]+0.25*alpha[0]*fUpwind[9]+0.25*fUpwind[0]*alpha[9]+0.25*alpha[5]*fUpwind[8]+0.25*fUpwind[5]*alpha[8]+0.25*alpha[2]*fUpwind[4]+0.25*fUpwind[2]*alpha[4]; 
  Ghat[10] = 0.22360679774997902*alpha[12]*fUpwind[31]+0.22360679774997896*alpha[9]*fUpwind[30]+0.22360679774997896*alpha[8]*fUpwind[29]+0.22360679774997902*alpha[4]*fUpwind[27]+0.22360679774997902*alpha[11]*fUpwind[23]+0.22360679774997896*alpha[7]*fUpwind[22]+0.22360679774997896*alpha[6]*fUpwind[21]+0.22360679774997902*alpha[3]*fUpwind[19]+0.25*alpha[5]*fUpwind[15]+0.25*alpha[2]*fUpwind[14]+0.25*alpha[1]*fUpwind[13]+0.25*alpha[11]*fUpwind[12]+0.25*fUpwind[11]*alpha[12]+0.25*alpha[0]*fUpwind[10]+0.25*alpha[7]*fUpwind[9]+0.25*fUpwind[7]*alpha[9]+0.25*alpha[6]*fUpwind[8]+0.25*fUpwind[6]*alpha[8]+0.25*alpha[3]*fUpwind[4]+0.25*fUpwind[3]*alpha[4]; 
  Ghat[11] = 0.22360679774997896*alpha[3]*fUpwind[20]+0.22360679774997902*alpha[6]*fUpwind[18]+0.22360679774997902*alpha[7]*fUpwind[17]+0.22360679774997896*alpha[11]*fUpwind[16]+0.25*alpha[4]*fUpwind[15]+0.25*alpha[8]*fUpwind[14]+0.25*alpha[9]*fUpwind[13]+0.25*fUpwind[10]*alpha[12]+0.25*alpha[0]*fUpwind[11]+0.25*fUpwind[0]*alpha[11]+0.25*alpha[1]*fUpwind[7]+0.25*fUpwind[1]*alpha[7]+0.25*alpha[2]*fUpwind[6]+0.25*fUpwind[2]*alpha[6]+0.25*alpha[3]*fUpwind[5]+0.25*fUpwind[3]*alpha[5]; 
  Ghat[12] = 0.22360679774997896*alpha[4]*fUpwind[28]+0.22360679774997902*alpha[8]*fUpwind[26]+0.22360679774997902*alpha[9]*fUpwind[25]+0.22360679774997896*alpha[12]*fUpwind[24]+0.25*alpha[3]*fUpwind[15]+0.25*alpha[6]*fUpwind[14]+0.25*alpha[7]*fUpwind[13]+0.25*alpha[0]*fUpwind[12]+0.25*fUpwind[0]*alpha[12]+0.25*fUpwind[10]*alpha[11]+0.25*alpha[1]*fUpwind[9]+0.25*fUpwind[1]*alpha[9]+0.25*alpha[2]*fUpwind[8]+0.25*fUpwind[2]*alpha[8]+0.25*alpha[4]*fUpwind[5]+0.25*fUpwind[4]*alpha[5]; 
  Ghat[13] = 0.22360679774997902*alpha[9]*fUpwind[31]+0.22360679774997896*alpha[12]*fUpwind[30]+0.22360679774997896*alpha[4]*fUpwind[29]+0.22360679774997902*alpha[8]*fUpwind[27]+0.22360679774997902*alpha[7]*fUpwind[23]+0.22360679774997896*alpha[11]*fUpwind[22]+0.22360679774997896*alpha[3]*fUpwind[21]+0.22360679774997902*alpha[6]*fUpwind[19]+0.25*alpha[2]*fUpwind[15]+0.25*alpha[5]*fUpwind[14]+0.25*alpha[0]*fUpwind[13]+0.25*alpha[7]*fUpwind[12]+0.25*fUpwind[7]*alpha[12]+0.25*alpha[9]*fUpwind[11]+0.25*fUpwind[9]*alpha[11]+0.25*alpha[1]*fUpwind[10]+0.25*alpha[3]*fUpwind[8]+0.25*fUpwind[3]*alpha[8]+0.25*alpha[4]*fUpwind[6]+0.25*fUpwind[4]*alpha[6]; 
  Ghat[14] = 0.22360679774997902*alpha[8]*fUpwind[31]+0.22360679774997896*alpha[4]*fUpwind[30]+0.22360679774997896*alpha[12]*fUpwind[29]+0.22360679774997902*alpha[9]*fUpwind[27]+0.22360679774997902*alpha[6]*fUpwind[23]+0.22360679774997896*alpha[3]*fUpwind[22]+0.22360679774997896*alpha[11]*fUpwind[21]+0.22360679774997902*alpha[7]*fUpwind[19]+0.25*alpha[1]*fUpwind[15]+0.25*alpha[0]*fUpwind[14]+0.25*alpha[5]*fUpwind[13]+0.25*alpha[6]*fUpwind[12]+0.25*fUpwind[6]*alpha[12]+0.25*alpha[8]*fUpwind[11]+0.25*fUpwind[8]*alpha[11]+0.25*alpha[2]*fUpwind[10]+0.25*alpha[3]*fUpwind[9]+0.25*fUpwind[3]*alpha[9]+0.25*alpha[4]*fUpwind[7]+0.25*fUpwind[4]*alpha[7]; 
  Ghat[15] = 0.22360679774997902*alpha[4]*fUpwind[31]+0.22360679774997896*alpha[8]*fUpwind[30]+0.22360679774997896*alpha[9]*fUpwind[29]+0.22360679774997902*alpha[12]*fUpwind[27]+0.22360679774997902*alpha[3]*fUpwind[23]+0.22360679774997896*alpha[6]*fUpwind[22]+0.22360679774997896*alpha[7]*fUpwind[21]+0.22360679774997902*alpha[11]*fUpwind[19]+0.25*alpha[0]*fUpwind[15]+0.25*alpha[1]*fUpwind[14]+0.25*alpha[2]*fUpwind[13]+0.25*alpha[3]*fUpwind[12]+0.25*fUpwind[3]*alpha[12]+0.25*alpha[4]*fUpwind[11]+0.25*fUpwind[4]*alpha[11]+0.25*alpha[5]*fUpwind[10]+0.25*alpha[6]*fUpwind[9]+0.25*fUpwind[6]*alpha[9]+0.25*alpha[7]*fUpwind[8]+0.25*fUpwind[7]*alpha[8]; 
  Ghat[16] = 0.25000000000000006*alpha[12]*fUpwind[23]+0.25*alpha[9]*fUpwind[22]+0.25*alpha[8]*fUpwind[21]+0.25*alpha[5]*fUpwind[20]+0.25000000000000006*alpha[4]*fUpwind[19]+0.25000000000000006*alpha[2]*fUpwind[18]+0.25000000000000006*alpha[1]*fUpwind[17]+0.25*alpha[0]*fUpwind[16]+0.22360679774997896*alpha[11]*fUpwind[11]+0.22360679774997896*alpha[7]*fUpwind[7]+0.22360679774997896*alpha[6]*fUpwind[6]+0.22360679774997896*alpha[3]*fUpwind[3]; 
  Ghat[17] = 0.25*alpha[9]*fUpwind[23]+0.25000000000000006*alpha[12]*fUpwind[22]+0.25000000000000006*alpha[4]*fUpwind[21]+0.25000000000000006*alpha[2]*fUpwind[20]+0.25*alpha[8]*fUpwind[19]+0.25*alpha[5]*fUpwind[18]+0.25*alpha[0]*fUpwind[17]+0.25000000000000006*alpha[1]*fUpwind[16]+0.22360679774997902*alpha[7]*fUpwind[11]+0.22360679774997902*fUpwind[7]*alpha[11]+0.22360679774997902*alpha[3]*fUpwind[6]+0.22360679774997902*fUpwind[3]*alpha[6]; 
  Ghat[18] = 0.25*alpha[8]*fUpwind[23]+0.25000000000000006*alpha[4]*fUpwind[22]+0.25000000000000006*alpha[12]*fUpwind[21]+0.25000000000000006*alpha[1]*fUpwind[20]+0.25*alpha[9]*fUpwind[19]+0.25*alpha[0]*fUpwind[18]+0.25*alpha[5]*fUpwind[17]+0.25000000000000006*alpha[2]*fUpwind[16]+0.22360679774997902*alpha[6]*fUpwind[11]+0.22360679774997902*fUpwind[6]*alpha[11]+0.22360679774997902*alpha[3]*fUpwind[7]+0.22360679774997902*fUpwind[3]*alpha[7]; 
  Ghat[19] = 0.25*alpha[5]*fUpwind[23]+0.25000000000000006*alpha[2]*fUpwind[22]+0.25000000000000006*alpha[1]*fUpwind[21]+0.25000000000000006*alpha[12]*fUpwind[20]+0.25*alpha[0]*fUpwind[19]+0.25*alpha[9]*fUpwind[18]+0.25*alpha[8]*fUpwind[17]+0.25000000000000006*alpha[4]*fUpwind[16]+0.22360679774997902*alpha[11]*fUpwind[15]+0.22360679774997902*alpha[7]*fUpwind[14]+0.22360679774997902*alpha[6]*fUpwind[13]+0.22360679774997902*alpha[3]*fUpwind[10]; 
  Ghat[20] = 0.25000000000000006*alpha[4]*fUpwind[23]+0.25*alpha[8]*fUpwind[22]+0.25*alpha[9]*fUpwind[21]+0.25*alpha[0]*fUpwind[20]+0.25000000000000006*alpha[12]*fUpwind[19]+0.25000000000000006*alpha[1]*fUpwind[18]+0.25000000000000006*alpha[2]*fUpwind[17]+0.25*alpha[5]*fUpwind[16]+0.22360679774997896*alpha[3]*fUpwind[11]+0.22360679774997896*fUpwind[3]*alpha[11]+0.22360679774997896*alpha[6]*fUpwind[7]+0.22360679774997896*fUpwind[6]*alpha[7]; 
  Ghat[21] = 0.25000000000000006*alpha[2]*fUpwind[23]+0.25*alpha[5]*fUpwind[22]+0.25*alpha[0]*fUpwind[21]+0.25*alpha[9]*fUpwind[20]+0.25000000000000006*alpha[1]*fUpwind[19]+0.25000000000000006*alpha[12]*fUpwind[18]+0.25000000000000006*alpha[4]*fUpwind[17]+0.25*alpha[8]*fUpwind[16]+0.22360679774997896*alpha[7]*fUpwind[15]+0.22360679774997896*alpha[11]*fUpwind[14]+0.22360679774997896*alpha[3]*fUpwind[13]+0.22360679774997896*alpha[6]*fUpwind[10]; 
  Ghat[22] = 0.25000000000000006*alpha[1]*fUpwind[23]+0.25*alpha[0]*fUpwind[22]+0.25*alpha[5]*fUpwind[21]+0.25*alpha[8]*fUpwind[20]+0.25000000000000006*alpha[2]*fUpwind[19]+0.25000000000000006*alpha[4]*fUpwind[18]+0.25000000000000006*alpha[12]*fUpwind[17]+0.25*alpha[9]*fUpwind[16]+0.22360679774997896*alpha[6]*fUpwind[15]+0.22360679774997896*alpha[3]*fUpwind[14]+0.22360679774997896*alpha[11]*fUpwind[13]+0.22360679774997896*alpha[7]*fUpwind[10]; 
  Ghat[23] = 0.25*alpha[0]*fUpwind[23]+0.25000000000000006*alpha[1]*fUpwind[22]+0.25000000000000006*alpha[2]*fUpwind[21]+0.25000000000000006*alpha[4]*fUpwind[20]+0.25*alpha[5]*fUpwind[19]+0.25*alpha[8]*fUpwind[18]+0.25*alpha[9]*fUpwind[17]+0.25000000000000006*alpha[12]*fUpwind[16]+0.22360679774997902*alpha[3]*fUpwind[15]+0.22360679774997902*alpha[6]*fUpwind[14]+0.22360679774997902*alpha[7]*fUpwind[13]+0.22360679774997902*fUpwind[10]*alpha[11]; 
  Ghat[24] = 0.25000000000000006*alpha[11]*fUpwind[31]+0.25*alpha[7]*fUpwind[30]+0.25*alpha[6]*fUpwind[29]+0.25*alpha[5]*fUpwind[28]+0.25000000000000006*alpha[3]*fUpwind[27]+0.25000000000000006*alpha[2]*fUpwind[26]+0.25000000000000006*alpha[1]*fUpwind[25]+0.25*alpha[0]*fUpwind[24]+0.22360679774997896*alpha[12]*fUpwind[12]+0.22360679774997896*alpha[9]*fUpwind[9]+0.22360679774997896*alpha[8]*fUpwind[8]+0.22360679774997896*alpha[4]*fUpwind[4]; 
  Ghat[25] = 0.25*alpha[7]*fUpwind[31]+0.25000000000000006*alpha[11]*fUpwind[30]+0.25000000000000006*alpha[3]*fUpwind[29]+0.25000000000000006*alpha[2]*fUpwind[28]+0.25*alpha[6]*fUpwind[27]+0.25*alpha[5]*fUpwind[26]+0.25*alpha[0]*fUpwind[25]+0.25000000000000006*alpha[1]*fUpwind[24]+0.22360679774997902*alpha[9]*fUpwind[12]+0.22360679774997902*fUpwind[9]*alpha[12]+0.22360679774997902*alpha[4]*fUpwind[8]+0.22360679774997902*fUpwind[4]*alpha[8]; 
  Ghat[26] = 0.25*alpha[6]*fUpwind[31]+0.25000000000000006*alpha[3]*fUpwind[30]+0.25000000000000006*alpha[11]*fUpwind[29]+0.25000000000000006*alpha[1]*fUpwind[28]+0.25*alpha[7]*fUpwind[27]+0.25*alpha[0]*fUpwind[26]+0.25*alpha[5]*fUpwind[25]+0.25000000000000006*alpha[2]*fUpwind[24]+0.22360679774997902*alpha[8]*fUpwind[12]+0.22360679774997902*fUpwind[8]*alpha[12]+0.22360679774997902*alpha[4]*fUpwind[9]+0.22360679774997902*fUpwind[4]*alpha[9]; 
  Ghat[27] = 0.25*alpha[5]*fUpwind[31]+0.25000000000000006*alpha[2]*fUpwind[30]+0.25000000000000006*alpha[1]*fUpwind[29]+0.25000000000000006*alpha[11]*fUpwind[28]+0.25*alpha[0]*fUpwind[27]+0.25*alpha[7]*fUpwind[26]+0.25*alpha[6]*fUpwind[25]+0.25000000000000006*alpha[3]*fUpwind[24]+0.22360679774997902*alpha[12]*fUpwind[15]+0.22360679774997902*alpha[9]*fUpwind[14]+0.22360679774997902*alpha[8]*fUpwind[13]+0.22360679774997902*alpha[4]*fUpwind[10]; 
  Ghat[28] = 0.25000000000000006*alpha[3]*fUpwind[31]+0.25*alpha[6]*fUpwind[30]+0.25*alpha[7]*fUpwind[29]+0.25*alpha[0]*fUpwind[28]+0.25000000000000006*alpha[11]*fUpwind[27]+0.25000000000000006*alpha[1]*fUpwind[26]+0.25000000000000006*alpha[2]*fUpwind[25]+0.25*alpha[5]*fUpwind[24]+0.22360679774997896*alpha[4]*fUpwind[12]+0.22360679774997896*fUpwind[4]*alpha[12]+0.22360679774997896*alpha[8]*fUpwind[9]+0.22360679774997896*fUpwind[8]*alpha[9]; 
  Ghat[29] = 0.25000000000000006*alpha[2]*fUpwind[31]+0.25*alpha[5]*fUpwind[30]+0.25*alpha[0]*fUpwind[29]+0.25*alpha[7]*fUpwind[28]+0.25000000000000006*alpha[1]*fUpwind[27]+0.25000000000000006*alpha[11]*fUpwind[26]+0.25000000000000006*alpha[3]*fUpwind[25]+0.25*alpha[6]*fUpwind[24]+0.22360679774997896*alpha[9]*fUpwind[15]+0.22360679774997896*alpha[12]*fUpwind[14]+0.22360679774997896*alpha[4]*fUpwind[13]+0.22360679774997896*alpha[8]*fUpwind[10]; 
  Ghat[30] = 0.25000000000000006*alpha[1]*fUpwind[31]+0.25*alpha[0]*fUpwind[30]+0.25*alpha[5]*fUpwind[29]+0.25*alpha[6]*fUpwind[28]+0.25000000000000006*alpha[2]*fUpwind[27]+0.25000000000000006*alpha[3]*fUpwind[26]+0.25000000000000006*alpha[11]*fUpwind[25]+0.25*alpha[7]*fUpwind[24]+0.22360679774997896*alpha[8]*fUpwind[15]+0.22360679774997896*alpha[4]*fUpwind[14]+0.22360679774997896*alpha[12]*fUpwind[13]+0.22360679774997896*alpha[9]*fUpwind[10]; 
  Ghat[31] = 0.25*alpha[0]*fUpwind[31]+0.25000000000000006*alpha[1]*fUpwind[30]+0.25000000000000006*alpha[2]*fUpwind[29]+0.25000000000000006*alpha[3]*fUpwind[28]+0.25*alpha[5]*fUpwind[27]+0.25*alpha[6]*fUpwind[26]+0.25*alpha[7]*fUpwind[25]+0.25000000000000006*alpha[11]*fUpwind[24]+0.22360679774997902*alpha[4]*fUpwind[15]+0.22360679774997902*alpha[8]*fUpwind[14]+0.22360679774997902*alpha[9]*fUpwind[13]+0.22360679774997902*fUpwind[10]*alpha[12]; 

  out[0] += -(0.7071067811865475*Ghat[0]*dv11); 
  out[1] += -(0.7071067811865475*Ghat[1]*dv11); 
  out[2] += -(0.7071067811865475*Ghat[2]*dv11); 
  out[3] += -(0.7071067811865475*Ghat[3]*dv11); 
  out[4] += -(1.224744871391589*Ghat[0]*dv11); 
  out[5] += -(0.7071067811865475*Ghat[4]*dv11); 
  out[6] += -(0.7071067811865475*Ghat[5]*dv11); 
  out[7] += -(0.7071067811865475*Ghat[6]*dv11); 
  out[8] += -(0.7071067811865475*Ghat[7]*dv11); 
  out[9] += -(1.224744871391589*Ghat[1]*dv11); 
  out[10] += -(1.224744871391589*Ghat[2]*dv11); 
  out[11] += -(1.224744871391589*Ghat[3]*dv11); 
  out[12] += -(0.7071067811865475*Ghat[8]*dv11); 
  out[13] += -(0.7071067811865475*Ghat[9]*dv11); 
  out[14] += -(0.7071067811865475*Ghat[10]*dv11); 
  out[15] += -(1.224744871391589*Ghat[4]*dv11); 
  out[16] += -(0.7071067811865475*Ghat[11]*dv11); 
  out[17] += -(1.224744871391589*Ghat[5]*dv11); 
  out[18] += -(1.224744871391589*Ghat[6]*dv11); 
  out[19] += -(1.224744871391589*Ghat[7]*dv11); 
  out[20] += -(0.7071067811865475*Ghat[12]*dv11); 
  out[21] += -(0.7071067811865475*Ghat[13]*dv11); 
  out[22] += -(0.7071067811865475*Ghat[14]*dv11); 
  out[23] += -(1.224744871391589*Ghat[8]*dv11); 
  out[24] += -(1.224744871391589*Ghat[9]*dv11); 
  out[25] += -(1.224744871391589*Ghat[10]*dv11); 
  out[26] += -(1.224744871391589*Ghat[11]*dv11); 
  out[27] += -(0.7071067811865475*Ghat[15]*dv11); 
  out[28] += -(1.224744871391589*Ghat[12]*dv11); 
  out[29] += -(1.224744871391589*Ghat[13]*dv11); 
  out[30] += -(1.224744871391589*Ghat[14]*dv11); 
  out[31] += -(1.224744871391589*Ghat[15]*dv11); 
  out[32] += -(0.7071067811865475*Ghat[16]*dv11); 
  out[33] += -(0.7071067811865475*Ghat[17]*dv11); 
  out[34] += -(0.7071067811865475*Ghat[18]*dv11); 
  out[35] += -(1.224744871391589*Ghat[16]*dv11); 
  out[36] += -(0.7071067811865475*Ghat[19]*dv11); 
  out[37] += -(0.7071067811865475*Ghat[20]*dv11); 
  out[38] += -(1.224744871391589*Ghat[17]*dv11); 
  out[39] += -(1.224744871391589*Ghat[18]*dv11); 
  out[40] += -(0.7071067811865475*Ghat[21]*dv11); 
  out[41] += -(0.7071067811865475*Ghat[22]*dv11); 
  out[42] += -(1.224744871391589*Ghat[19]*dv11); 
  out[43] += -(1.224744871391589*Ghat[20]*dv11); 
  out[44] += -(0.7071067811865475*Ghat[23]*dv11); 
  out[45] += -(1.224744871391589*Ghat[21]*dv11); 
  out[46] += -(1.224744871391589*Ghat[22]*dv11); 
  out[47] += -(1.224744871391589*Ghat[23]*dv11); 
  out[48] += -(1.5811388300841895*Ghat[0]*dv11); 
  out[49] += -(1.5811388300841898*Ghat[1]*dv11); 
  out[50] += -(1.5811388300841898*Ghat[2]*dv11); 
  out[51] += -(1.5811388300841898*Ghat[3]*dv11); 
  out[52] += -(1.5811388300841898*Ghat[4]*dv11); 
  out[53] += -(1.5811388300841895*Ghat[5]*dv11); 
  out[54] += -(1.5811388300841895*Ghat[6]*dv11); 
  out[55] += -(1.5811388300841895*Ghat[7]*dv11); 
  out[56] += -(1.5811388300841895*Ghat[8]*dv11); 
  out[57] += -(1.5811388300841895*Ghat[9]*dv11); 
  out[58] += -(1.5811388300841895*Ghat[10]*dv11); 
  out[59] += -(1.5811388300841898*Ghat[11]*dv11); 
  out[60] += -(1.5811388300841898*Ghat[12]*dv11); 
  out[61] += -(1.5811388300841898*Ghat[13]*dv11); 
  out[62] += -(1.5811388300841898*Ghat[14]*dv11); 
  out[63] += -(1.5811388300841895*Ghat[15]*dv11); 
  out[64] += -(0.7071067811865475*Ghat[24]*dv11); 
  out[65] += -(0.7071067811865475*Ghat[25]*dv11); 
  out[66] += -(0.7071067811865475*Ghat[26]*dv11); 
  out[67] += -(0.7071067811865475*Ghat[27]*dv11); 
  out[68] += -(1.224744871391589*Ghat[24]*dv11); 
  out[69] += -(0.7071067811865475*Ghat[28]*dv11); 
  out[70] += -(0.7071067811865475*Ghat[29]*dv11); 
  out[71] += -(0.7071067811865475*Ghat[30]*dv11); 
  out[72] += -(1.224744871391589*Ghat[25]*dv11); 
  out[73] += -(1.224744871391589*Ghat[26]*dv11); 
  out[74] += -(1.224744871391589*Ghat[27]*dv11); 
  out[75] += -(0.7071067811865475*Ghat[31]*dv11); 
  out[76] += -(1.224744871391589*Ghat[28]*dv11); 
  out[77] += -(1.224744871391589*Ghat[29]*dv11); 
  out[78] += -(1.224744871391589*Ghat[30]*dv11); 
  out[79] += -(1.224744871391589*Ghat[31]*dv11); 

  } else { 

  if (-(0.33541019662496846*(alpha[12]+alpha[11]))+0.33541019662496846*(alpha[9]+alpha[8]+alpha[7]+alpha[6])+0.25*alpha[5]-0.33541019662496846*(alpha[4]+alpha[3])-0.25*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[0] = hyb_2x3v_p1_surfx4_eval_quad_node_0_r(fEdge); 
  } else { 
    fUpwindQuad[0] = hyb_2x3v_p1_surfx4_eval_quad_node_0_l(fSkin); 
  } 
  if (-(0.33541019662496846*alpha[11])+0.33541019662496846*(alpha[7]+alpha[6])+0.25*alpha[5]-0.33541019662496846*alpha[3]-0.25*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[1] = hyb_2x3v_p1_surfx4_eval_quad_node_1_r(fEdge); 
  } else { 
    fUpwindQuad[1] = hyb_2x3v_p1_surfx4_eval_quad_node_1_l(fSkin); 
  } 
  if (0.33541019662496846*alpha[12]-0.33541019662496846*(alpha[11]+alpha[9]+alpha[8])+0.33541019662496846*(alpha[7]+alpha[6])+0.25*alpha[5]+0.33541019662496846*alpha[4]-0.33541019662496846*alpha[3]-0.25*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[2] = hyb_2x3v_p1_surfx4_eval_quad_node_2_r(fEdge); 
  } else { 
    fUpwindQuad[2] = hyb_2x3v_p1_surfx4_eval_quad_node_2_l(fSkin); 
  } 
  if (-(0.33541019662496846*alpha[12])+0.33541019662496846*(alpha[9]+alpha[8])+0.25*alpha[5]-0.33541019662496846*alpha[4]-0.25*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[3] = hyb_2x3v_p1_surfx4_eval_quad_node_3_r(fEdge); 
  } else { 
    fUpwindQuad[3] = hyb_2x3v_p1_surfx4_eval_quad_node_3_l(fSkin); 
  } 
  if (0.25*alpha[5]-0.25*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[4] = hyb_2x3v_p1_surfx4_eval_quad_node_4_r(fEdge); 
  } else { 
    fUpwindQuad[4] = hyb_2x3v_p1_surfx4_eval_quad_node_4_l(fSkin); 
  } 
  if (0.33541019662496846*alpha[12]-0.33541019662496846*(alpha[9]+alpha[8])+0.25*alpha[5]+0.33541019662496846*alpha[4]-0.25*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[5] = hyb_2x3v_p1_surfx4_eval_quad_node_5_r(fEdge); 
  } else { 
    fUpwindQuad[5] = hyb_2x3v_p1_surfx4_eval_quad_node_5_l(fSkin); 
  } 
  if (-(0.33541019662496846*alpha[12])+0.33541019662496846*(alpha[11]+alpha[9]+alpha[8])-0.33541019662496846*(alpha[7]+alpha[6])+0.25*alpha[5]-0.33541019662496846*alpha[4]+0.33541019662496846*alpha[3]-0.25*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[6] = hyb_2x3v_p1_surfx4_eval_quad_node_6_r(fEdge); 
  } else { 
    fUpwindQuad[6] = hyb_2x3v_p1_surfx4_eval_quad_node_6_l(fSkin); 
  } 
  if (0.33541019662496846*alpha[11]-0.33541019662496846*(alpha[7]+alpha[6])+0.25*alpha[5]+0.33541019662496846*alpha[3]-0.25*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[7] = hyb_2x3v_p1_surfx4_eval_quad_node_7_r(fEdge); 
  } else { 
    fUpwindQuad[7] = hyb_2x3v_p1_surfx4_eval_quad_node_7_l(fSkin); 
  } 
  if (0.33541019662496846*(alpha[12]+alpha[11])-0.33541019662496846*(alpha[9]+alpha[8]+alpha[7]+alpha[6])+0.25*alpha[5]+0.33541019662496846*(alpha[4]+alpha[3])-0.25*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad[8] = hyb_2x3v_p1_surfx4_eval_quad_node_8_r(fEdge); 
  } else { 
    fUpwindQuad[8] = hyb_2x3v_p1_surfx4_eval_quad_node_8_l(fSkin); 
  } 
  if (0.33541019662496846*(alpha[12]+alpha[11])-0.33541019662496846*alpha[9]+0.33541019662496846*alpha[8]-0.33541019662496846*alpha[7]+0.33541019662496846*alpha[6]-0.25*alpha[5]-0.33541019662496846*(alpha[4]+alpha[3])+0.25*alpha[2]-0.25*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[9] = hyb_2x3v_p1_surfx4_eval_quad_node_9_r(fEdge); 
  } else { 
    fUpwindQuad[9] = hyb_2x3v_p1_surfx4_eval_quad_node_9_l(fSkin); 
  } 
  if (0.33541019662496846*alpha[11]-0.33541019662496846*alpha[7]+0.33541019662496846*alpha[6]-0.25*alpha[5]-0.33541019662496846*alpha[3]+0.25*alpha[2]-0.25*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[10] = hyb_2x3v_p1_surfx4_eval_quad_node_10_r(fEdge); 
  } else { 
    fUpwindQuad[10] = hyb_2x3v_p1_surfx4_eval_quad_node_10_l(fSkin); 
  } 
  if (-(0.33541019662496846*alpha[12])+0.33541019662496846*(alpha[11]+alpha[9])-0.33541019662496846*(alpha[8]+alpha[7])+0.33541019662496846*alpha[6]-0.25*alpha[5]+0.33541019662496846*alpha[4]-0.33541019662496846*alpha[3]+0.25*alpha[2]-0.25*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[11] = hyb_2x3v_p1_surfx4_eval_quad_node_11_r(fEdge); 
  } else { 
    fUpwindQuad[11] = hyb_2x3v_p1_surfx4_eval_quad_node_11_l(fSkin); 
  } 
  if (0.33541019662496846*alpha[12]-0.33541019662496846*alpha[9]+0.33541019662496846*alpha[8]-0.25*alpha[5]-0.33541019662496846*alpha[4]+0.25*alpha[2]-0.25*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[12] = hyb_2x3v_p1_surfx4_eval_quad_node_12_r(fEdge); 
  } else { 
    fUpwindQuad[12] = hyb_2x3v_p1_surfx4_eval_quad_node_12_l(fSkin); 
  } 
  if (-(0.25*alpha[5])+0.25*alpha[2]-0.25*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[13] = hyb_2x3v_p1_surfx4_eval_quad_node_13_r(fEdge); 
  } else { 
    fUpwindQuad[13] = hyb_2x3v_p1_surfx4_eval_quad_node_13_l(fSkin); 
  } 
  if (-(0.33541019662496846*alpha[12])+0.33541019662496846*alpha[9]-0.33541019662496846*alpha[8]-0.25*alpha[5]+0.33541019662496846*alpha[4]+0.25*alpha[2]-0.25*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[14] = hyb_2x3v_p1_surfx4_eval_quad_node_14_r(fEdge); 
  } else { 
    fUpwindQuad[14] = hyb_2x3v_p1_surfx4_eval_quad_node_14_l(fSkin); 
  } 
  if (0.33541019662496846*alpha[12]-0.33541019662496846*(alpha[11]+alpha[9])+0.33541019662496846*(alpha[8]+alpha[7])-0.33541019662496846*alpha[6]-0.25*alpha[5]-0.33541019662496846*alpha[4]+0.33541019662496846*alpha[3]+0.25*alpha[2]-0.25*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[15] = hyb_2x3v_p1_surfx4_eval_quad_node_15_r(fEdge); 
  } else { 
    fUpwindQuad[15] = hyb_2x3v_p1_surfx4_eval_quad_node_15_l(fSkin); 
  } 
  if (-(0.33541019662496846*alpha[11])+0.33541019662496846*alpha[7]-0.33541019662496846*alpha[6]-0.25*alpha[5]+0.33541019662496846*alpha[3]+0.25*alpha[2]-0.25*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[16] = hyb_2x3v_p1_surfx4_eval_quad_node_16_r(fEdge); 
  } else { 
    fUpwindQuad[16] = hyb_2x3v_p1_surfx4_eval_quad_node_16_l(fSkin); 
  } 
  if (-(0.33541019662496846*(alpha[12]+alpha[11]))+0.33541019662496846*alpha[9]-0.33541019662496846*alpha[8]+0.33541019662496846*alpha[7]-0.33541019662496846*alpha[6]-0.25*alpha[5]+0.33541019662496846*(alpha[4]+alpha[3])+0.25*alpha[2]-0.25*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad[17] = hyb_2x3v_p1_surfx4_eval_quad_node_17_r(fEdge); 
  } else { 
    fUpwindQuad[17] = hyb_2x3v_p1_surfx4_eval_quad_node_17_l(fSkin); 
  } 
  if (0.33541019662496846*(alpha[12]+alpha[11]+alpha[9])-0.33541019662496846*alpha[8]+0.33541019662496846*alpha[7]-0.33541019662496846*alpha[6]-0.25*alpha[5]-0.33541019662496846*(alpha[4]+alpha[3])-0.25*alpha[2]+0.25*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[18] = hyb_2x3v_p1_surfx4_eval_quad_node_18_r(fEdge); 
  } else { 
    fUpwindQuad[18] = hyb_2x3v_p1_surfx4_eval_quad_node_18_l(fSkin); 
  } 
  if (0.33541019662496846*(alpha[11]+alpha[7])-0.33541019662496846*alpha[6]-0.25*alpha[5]-0.33541019662496846*alpha[3]-0.25*alpha[2]+0.25*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[19] = hyb_2x3v_p1_surfx4_eval_quad_node_19_r(fEdge); 
  } else { 
    fUpwindQuad[19] = hyb_2x3v_p1_surfx4_eval_quad_node_19_l(fSkin); 
  } 
  if (-(0.33541019662496846*alpha[12])+0.33541019662496846*alpha[11]-0.33541019662496846*alpha[9]+0.33541019662496846*(alpha[8]+alpha[7])-0.33541019662496846*alpha[6]-0.25*alpha[5]+0.33541019662496846*alpha[4]-0.33541019662496846*alpha[3]-0.25*alpha[2]+0.25*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[20] = hyb_2x3v_p1_surfx4_eval_quad_node_20_r(fEdge); 
  } else { 
    fUpwindQuad[20] = hyb_2x3v_p1_surfx4_eval_quad_node_20_l(fSkin); 
  } 
  if (0.33541019662496846*(alpha[12]+alpha[9])-0.33541019662496846*alpha[8]-0.25*alpha[5]-0.33541019662496846*alpha[4]-0.25*alpha[2]+0.25*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[21] = hyb_2x3v_p1_surfx4_eval_quad_node_21_r(fEdge); 
  } else { 
    fUpwindQuad[21] = hyb_2x3v_p1_surfx4_eval_quad_node_21_l(fSkin); 
  } 
  if (0.25*(alpha[1]+alpha[0])-0.25*(alpha[5]+alpha[2]) > 0) { 
    fUpwindQuad[22] = hyb_2x3v_p1_surfx4_eval_quad_node_22_r(fEdge); 
  } else { 
    fUpwindQuad[22] = hyb_2x3v_p1_surfx4_eval_quad_node_22_l(fSkin); 
  } 
  if (-(0.33541019662496846*(alpha[12]+alpha[9]))+0.33541019662496846*alpha[8]-0.25*alpha[5]+0.33541019662496846*alpha[4]-0.25*alpha[2]+0.25*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[23] = hyb_2x3v_p1_surfx4_eval_quad_node_23_r(fEdge); 
  } else { 
    fUpwindQuad[23] = hyb_2x3v_p1_surfx4_eval_quad_node_23_l(fSkin); 
  } 
  if (0.33541019662496846*alpha[12]-0.33541019662496846*alpha[11]+0.33541019662496846*alpha[9]-0.33541019662496846*(alpha[8]+alpha[7])+0.33541019662496846*alpha[6]-0.25*alpha[5]-0.33541019662496846*alpha[4]+0.33541019662496846*alpha[3]-0.25*alpha[2]+0.25*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[24] = hyb_2x3v_p1_surfx4_eval_quad_node_24_r(fEdge); 
  } else { 
    fUpwindQuad[24] = hyb_2x3v_p1_surfx4_eval_quad_node_24_l(fSkin); 
  } 
  if (-(0.33541019662496846*(alpha[11]+alpha[7]))+0.33541019662496846*alpha[6]-0.25*alpha[5]+0.33541019662496846*alpha[3]-0.25*alpha[2]+0.25*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[25] = hyb_2x3v_p1_surfx4_eval_quad_node_25_r(fEdge); 
  } else { 
    fUpwindQuad[25] = hyb_2x3v_p1_surfx4_eval_quad_node_25_l(fSkin); 
  } 
  if (-(0.33541019662496846*(alpha[12]+alpha[11]+alpha[9]))+0.33541019662496846*alpha[8]-0.33541019662496846*alpha[7]+0.33541019662496846*alpha[6]-0.25*alpha[5]+0.33541019662496846*(alpha[4]+alpha[3])-0.25*alpha[2]+0.25*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[26] = hyb_2x3v_p1_surfx4_eval_quad_node_26_r(fEdge); 
  } else { 
    fUpwindQuad[26] = hyb_2x3v_p1_surfx4_eval_quad_node_26_l(fSkin); 
  } 
  if (-(0.33541019662496846*(alpha[12]+alpha[11]+alpha[9]+alpha[8]+alpha[7]+alpha[6]))+0.25*alpha[5]-0.33541019662496846*(alpha[4]+alpha[3])+0.25*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[27] = hyb_2x3v_p1_surfx4_eval_quad_node_27_r(fEdge); 
  } else { 
    fUpwindQuad[27] = hyb_2x3v_p1_surfx4_eval_quad_node_27_l(fSkin); 
  } 
  if (-(0.33541019662496846*(alpha[11]+alpha[7]+alpha[6]))+0.25*alpha[5]-0.33541019662496846*alpha[3]+0.25*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[28] = hyb_2x3v_p1_surfx4_eval_quad_node_28_r(fEdge); 
  } else { 
    fUpwindQuad[28] = hyb_2x3v_p1_surfx4_eval_quad_node_28_l(fSkin); 
  } 
  if (0.33541019662496846*alpha[12]-0.33541019662496846*alpha[11]+0.33541019662496846*(alpha[9]+alpha[8])-0.33541019662496846*(alpha[7]+alpha[6])+0.25*alpha[5]+0.33541019662496846*alpha[4]-0.33541019662496846*alpha[3]+0.25*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[29] = hyb_2x3v_p1_surfx4_eval_quad_node_29_r(fEdge); 
  } else { 
    fUpwindQuad[29] = hyb_2x3v_p1_surfx4_eval_quad_node_29_l(fSkin); 
  } 
  if (-(0.33541019662496846*(alpha[12]+alpha[9]+alpha[8]))+0.25*alpha[5]-0.33541019662496846*alpha[4]+0.25*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[30] = hyb_2x3v_p1_surfx4_eval_quad_node_30_r(fEdge); 
  } else { 
    fUpwindQuad[30] = hyb_2x3v_p1_surfx4_eval_quad_node_30_l(fSkin); 
  } 
  if (0.25*(alpha[5]+alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[31] = hyb_2x3v_p1_surfx4_eval_quad_node_31_r(fEdge); 
  } else { 
    fUpwindQuad[31] = hyb_2x3v_p1_surfx4_eval_quad_node_31_l(fSkin); 
  } 
  if (0.33541019662496846*(alpha[12]+alpha[9]+alpha[8])+0.25*alpha[5]+0.33541019662496846*alpha[4]+0.25*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[32] = hyb_2x3v_p1_surfx4_eval_quad_node_32_r(fEdge); 
  } else { 
    fUpwindQuad[32] = hyb_2x3v_p1_surfx4_eval_quad_node_32_l(fSkin); 
  } 
  if (-(0.33541019662496846*alpha[12])+0.33541019662496846*alpha[11]-0.33541019662496846*(alpha[9]+alpha[8])+0.33541019662496846*(alpha[7]+alpha[6])+0.25*alpha[5]-0.33541019662496846*alpha[4]+0.33541019662496846*alpha[3]+0.25*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[33] = hyb_2x3v_p1_surfx4_eval_quad_node_33_r(fEdge); 
  } else { 
    fUpwindQuad[33] = hyb_2x3v_p1_surfx4_eval_quad_node_33_l(fSkin); 
  } 
  if (0.33541019662496846*(alpha[11]+alpha[7]+alpha[6])+0.25*alpha[5]+0.33541019662496846*alpha[3]+0.25*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[34] = hyb_2x3v_p1_surfx4_eval_quad_node_34_r(fEdge); 
  } else { 
    fUpwindQuad[34] = hyb_2x3v_p1_surfx4_eval_quad_node_34_l(fSkin); 
  } 
  if (0.33541019662496846*(alpha[12]+alpha[11]+alpha[9]+alpha[8]+alpha[7]+alpha[6])+0.25*alpha[5]+0.33541019662496846*(alpha[4]+alpha[3])+0.25*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[35] = hyb_2x3v_p1_surfx4_eval_quad_node_35_r(fEdge); 
  } else { 
    fUpwindQuad[35] = hyb_2x3v_p1_surfx4_eval_quad_node_35_l(fSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_2x3v_p1_vdir_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.25*alpha[12]*fUpwind[12]+0.25*alpha[11]*fUpwind[11]+0.25*alpha[9]*fUpwind[9]+0.25*alpha[8]*fUpwind[8]+0.25*alpha[7]*fUpwind[7]+0.25*alpha[6]*fUpwind[6]+0.25*alpha[5]*fUpwind[5]+0.25*alpha[4]*fUpwind[4]+0.25*alpha[3]*fUpwind[3]+0.25*alpha[2]*fUpwind[2]+0.25*alpha[1]*fUpwind[1]+0.25*alpha[0]*fUpwind[0]; 
  Ghat[1] = 0.25*alpha[9]*fUpwind[12]+0.25*fUpwind[9]*alpha[12]+0.25*alpha[7]*fUpwind[11]+0.25*fUpwind[7]*alpha[11]+0.25*alpha[4]*fUpwind[8]+0.25*fUpwind[4]*alpha[8]+0.25*alpha[3]*fUpwind[6]+0.25*fUpwind[3]*alpha[6]+0.25*alpha[2]*fUpwind[5]+0.25*fUpwind[2]*alpha[5]+0.25*alpha[0]*fUpwind[1]+0.25*fUpwind[0]*alpha[1]; 
  Ghat[2] = 0.25*alpha[8]*fUpwind[12]+0.25*fUpwind[8]*alpha[12]+0.25*alpha[6]*fUpwind[11]+0.25*fUpwind[6]*alpha[11]+0.25*alpha[4]*fUpwind[9]+0.25*fUpwind[4]*alpha[9]+0.25*alpha[3]*fUpwind[7]+0.25*fUpwind[3]*alpha[7]+0.25*alpha[1]*fUpwind[5]+0.25*fUpwind[1]*alpha[5]+0.25*alpha[0]*fUpwind[2]+0.25*fUpwind[0]*alpha[2]; 
  Ghat[3] = 0.22360679774997896*alpha[11]*fUpwind[20]+0.22360679774997902*alpha[7]*fUpwind[18]+0.22360679774997902*alpha[6]*fUpwind[17]+0.22360679774997896*alpha[3]*fUpwind[16]+0.25*alpha[12]*fUpwind[15]+0.25*alpha[9]*fUpwind[14]+0.25*alpha[8]*fUpwind[13]+0.25*alpha[5]*fUpwind[11]+0.25*fUpwind[5]*alpha[11]+0.25*alpha[4]*fUpwind[10]+0.25*alpha[2]*fUpwind[7]+0.25*fUpwind[2]*alpha[7]+0.25*alpha[1]*fUpwind[6]+0.25*fUpwind[1]*alpha[6]+0.25*alpha[0]*fUpwind[3]+0.25*fUpwind[0]*alpha[3]; 
  Ghat[4] = 0.22360679774997896*alpha[12]*fUpwind[28]+0.22360679774997902*alpha[9]*fUpwind[26]+0.22360679774997902*alpha[8]*fUpwind[25]+0.22360679774997896*alpha[4]*fUpwind[24]+0.25*alpha[11]*fUpwind[15]+0.25*alpha[7]*fUpwind[14]+0.25*alpha[6]*fUpwind[13]+0.25*alpha[5]*fUpwind[12]+0.25*fUpwind[5]*alpha[12]+0.25*alpha[3]*fUpwind[10]+0.25*alpha[2]*fUpwind[9]+0.25*fUpwind[2]*alpha[9]+0.25*alpha[1]*fUpwind[8]+0.25*fUpwind[1]*alpha[8]+0.25*alpha[0]*fUpwind[4]+0.25*fUpwind[0]*alpha[4]; 
  Ghat[5] = 0.25*alpha[4]*fUpwind[12]+0.25*fUpwind[4]*alpha[12]+0.25*alpha[3]*fUpwind[11]+0.25*fUpwind[3]*alpha[11]+0.25*alpha[8]*fUpwind[9]+0.25*fUpwind[8]*alpha[9]+0.25*alpha[6]*fUpwind[7]+0.25*fUpwind[6]*alpha[7]+0.25*alpha[0]*fUpwind[5]+0.25*fUpwind[0]*alpha[5]+0.25*alpha[1]*fUpwind[2]+0.25*fUpwind[1]*alpha[2]; 
  Ghat[6] = 0.22360679774997896*alpha[7]*fUpwind[20]+0.22360679774997902*alpha[11]*fUpwind[18]+0.22360679774997902*alpha[3]*fUpwind[17]+0.22360679774997896*alpha[6]*fUpwind[16]+0.25*alpha[9]*fUpwind[15]+0.25*alpha[12]*fUpwind[14]+0.25*alpha[4]*fUpwind[13]+0.25*alpha[2]*fUpwind[11]+0.25*fUpwind[2]*alpha[11]+0.25*alpha[8]*fUpwind[10]+0.25*alpha[5]*fUpwind[7]+0.25*fUpwind[5]*alpha[7]+0.25*alpha[0]*fUpwind[6]+0.25*fUpwind[0]*alpha[6]+0.25*alpha[1]*fUpwind[3]+0.25*fUpwind[1]*alpha[3]; 
  Ghat[7] = 0.22360679774997896*alpha[6]*fUpwind[20]+0.22360679774997902*alpha[3]*fUpwind[18]+0.22360679774997902*alpha[11]*fUpwind[17]+0.22360679774997896*alpha[7]*fUpwind[16]+0.25*alpha[8]*fUpwind[15]+0.25*alpha[4]*fUpwind[14]+0.25*alpha[12]*fUpwind[13]+0.25*alpha[1]*fUpwind[11]+0.25*fUpwind[1]*alpha[11]+0.25*alpha[9]*fUpwind[10]+0.25*alpha[0]*fUpwind[7]+0.25*fUpwind[0]*alpha[7]+0.25*alpha[5]*fUpwind[6]+0.25*fUpwind[5]*alpha[6]+0.25*alpha[2]*fUpwind[3]+0.25*fUpwind[2]*alpha[3]; 
  Ghat[8] = 0.22360679774997896*alpha[9]*fUpwind[28]+0.22360679774997902*alpha[12]*fUpwind[26]+0.22360679774997902*alpha[4]*fUpwind[25]+0.22360679774997896*alpha[8]*fUpwind[24]+0.25*alpha[7]*fUpwind[15]+0.25*alpha[11]*fUpwind[14]+0.25*alpha[3]*fUpwind[13]+0.25*alpha[2]*fUpwind[12]+0.25*fUpwind[2]*alpha[12]+0.25*alpha[6]*fUpwind[10]+0.25*alpha[5]*fUpwind[9]+0.25*fUpwind[5]*alpha[9]+0.25*alpha[0]*fUpwind[8]+0.25*fUpwind[0]*alpha[8]+0.25*alpha[1]*fUpwind[4]+0.25*fUpwind[1]*alpha[4]; 
  Ghat[9] = 0.22360679774997896*alpha[8]*fUpwind[28]+0.22360679774997902*alpha[4]*fUpwind[26]+0.22360679774997902*alpha[12]*fUpwind[25]+0.22360679774997896*alpha[9]*fUpwind[24]+0.25*alpha[6]*fUpwind[15]+0.25*alpha[3]*fUpwind[14]+0.25*alpha[11]*fUpwind[13]+0.25*alpha[1]*fUpwind[12]+0.25*fUpwind[1]*alpha[12]+0.25*alpha[7]*fUpwind[10]+0.25*alpha[0]*fUpwind[9]+0.25*fUpwind[0]*alpha[9]+0.25*alpha[5]*fUpwind[8]+0.25*fUpwind[5]*alpha[8]+0.25*alpha[2]*fUpwind[4]+0.25*fUpwind[2]*alpha[4]; 
  Ghat[10] = 0.22360679774997902*alpha[12]*fUpwind[31]+0.22360679774997896*alpha[9]*fUpwind[30]+0.22360679774997896*alpha[8]*fUpwind[29]+0.22360679774997902*alpha[4]*fUpwind[27]+0.22360679774997902*alpha[11]*fUpwind[23]+0.22360679774997896*alpha[7]*fUpwind[22]+0.22360679774997896*alpha[6]*fUpwind[21]+0.22360679774997902*alpha[3]*fUpwind[19]+0.25*alpha[5]*fUpwind[15]+0.25*alpha[2]*fUpwind[14]+0.25*alpha[1]*fUpwind[13]+0.25*alpha[11]*fUpwind[12]+0.25*fUpwind[11]*alpha[12]+0.25*alpha[0]*fUpwind[10]+0.25*alpha[7]*fUpwind[9]+0.25*fUpwind[7]*alpha[9]+0.25*alpha[6]*fUpwind[8]+0.25*fUpwind[6]*alpha[8]+0.25*alpha[3]*fUpwind[4]+0.25*fUpwind[3]*alpha[4]; 
  Ghat[11] = 0.22360679774997896*alpha[3]*fUpwind[20]+0.22360679774997902*alpha[6]*fUpwind[18]+0.22360679774997902*alpha[7]*fUpwind[17]+0.22360679774997896*alpha[11]*fUpwind[16]+0.25*alpha[4]*fUpwind[15]+0.25*alpha[8]*fUpwind[14]+0.25*alpha[9]*fUpwind[13]+0.25*fUpwind[10]*alpha[12]+0.25*alpha[0]*fUpwind[11]+0.25*fUpwind[0]*alpha[11]+0.25*alpha[1]*fUpwind[7]+0.25*fUpwind[1]*alpha[7]+0.25*alpha[2]*fUpwind[6]+0.25*fUpwind[2]*alpha[6]+0.25*alpha[3]*fUpwind[5]+0.25*fUpwind[3]*alpha[5]; 
  Ghat[12] = 0.22360679774997896*alpha[4]*fUpwind[28]+0.22360679774997902*alpha[8]*fUpwind[26]+0.22360679774997902*alpha[9]*fUpwind[25]+0.22360679774997896*alpha[12]*fUpwind[24]+0.25*alpha[3]*fUpwind[15]+0.25*alpha[6]*fUpwind[14]+0.25*alpha[7]*fUpwind[13]+0.25*alpha[0]*fUpwind[12]+0.25*fUpwind[0]*alpha[12]+0.25*fUpwind[10]*alpha[11]+0.25*alpha[1]*fUpwind[9]+0.25*fUpwind[1]*alpha[9]+0.25*alpha[2]*fUpwind[8]+0.25*fUpwind[2]*alpha[8]+0.25*alpha[4]*fUpwind[5]+0.25*fUpwind[4]*alpha[5]; 
  Ghat[13] = 0.22360679774997902*alpha[9]*fUpwind[31]+0.22360679774997896*alpha[12]*fUpwind[30]+0.22360679774997896*alpha[4]*fUpwind[29]+0.22360679774997902*alpha[8]*fUpwind[27]+0.22360679774997902*alpha[7]*fUpwind[23]+0.22360679774997896*alpha[11]*fUpwind[22]+0.22360679774997896*alpha[3]*fUpwind[21]+0.22360679774997902*alpha[6]*fUpwind[19]+0.25*alpha[2]*fUpwind[15]+0.25*alpha[5]*fUpwind[14]+0.25*alpha[0]*fUpwind[13]+0.25*alpha[7]*fUpwind[12]+0.25*fUpwind[7]*alpha[12]+0.25*alpha[9]*fUpwind[11]+0.25*fUpwind[9]*alpha[11]+0.25*alpha[1]*fUpwind[10]+0.25*alpha[3]*fUpwind[8]+0.25*fUpwind[3]*alpha[8]+0.25*alpha[4]*fUpwind[6]+0.25*fUpwind[4]*alpha[6]; 
  Ghat[14] = 0.22360679774997902*alpha[8]*fUpwind[31]+0.22360679774997896*alpha[4]*fUpwind[30]+0.22360679774997896*alpha[12]*fUpwind[29]+0.22360679774997902*alpha[9]*fUpwind[27]+0.22360679774997902*alpha[6]*fUpwind[23]+0.22360679774997896*alpha[3]*fUpwind[22]+0.22360679774997896*alpha[11]*fUpwind[21]+0.22360679774997902*alpha[7]*fUpwind[19]+0.25*alpha[1]*fUpwind[15]+0.25*alpha[0]*fUpwind[14]+0.25*alpha[5]*fUpwind[13]+0.25*alpha[6]*fUpwind[12]+0.25*fUpwind[6]*alpha[12]+0.25*alpha[8]*fUpwind[11]+0.25*fUpwind[8]*alpha[11]+0.25*alpha[2]*fUpwind[10]+0.25*alpha[3]*fUpwind[9]+0.25*fUpwind[3]*alpha[9]+0.25*alpha[4]*fUpwind[7]+0.25*fUpwind[4]*alpha[7]; 
  Ghat[15] = 0.22360679774997902*alpha[4]*fUpwind[31]+0.22360679774997896*alpha[8]*fUpwind[30]+0.22360679774997896*alpha[9]*fUpwind[29]+0.22360679774997902*alpha[12]*fUpwind[27]+0.22360679774997902*alpha[3]*fUpwind[23]+0.22360679774997896*alpha[6]*fUpwind[22]+0.22360679774997896*alpha[7]*fUpwind[21]+0.22360679774997902*alpha[11]*fUpwind[19]+0.25*alpha[0]*fUpwind[15]+0.25*alpha[1]*fUpwind[14]+0.25*alpha[2]*fUpwind[13]+0.25*alpha[3]*fUpwind[12]+0.25*fUpwind[3]*alpha[12]+0.25*alpha[4]*fUpwind[11]+0.25*fUpwind[4]*alpha[11]+0.25*alpha[5]*fUpwind[10]+0.25*alpha[6]*fUpwind[9]+0.25*fUpwind[6]*alpha[9]+0.25*alpha[7]*fUpwind[8]+0.25*fUpwind[7]*alpha[8]; 
  Ghat[16] = 0.25000000000000006*alpha[12]*fUpwind[23]+0.25*alpha[9]*fUpwind[22]+0.25*alpha[8]*fUpwind[21]+0.25*alpha[5]*fUpwind[20]+0.25000000000000006*alpha[4]*fUpwind[19]+0.25000000000000006*alpha[2]*fUpwind[18]+0.25000000000000006*alpha[1]*fUpwind[17]+0.25*alpha[0]*fUpwind[16]+0.22360679774997896*alpha[11]*fUpwind[11]+0.22360679774997896*alpha[7]*fUpwind[7]+0.22360679774997896*alpha[6]*fUpwind[6]+0.22360679774997896*alpha[3]*fUpwind[3]; 
  Ghat[17] = 0.25*alpha[9]*fUpwind[23]+0.25000000000000006*alpha[12]*fUpwind[22]+0.25000000000000006*alpha[4]*fUpwind[21]+0.25000000000000006*alpha[2]*fUpwind[20]+0.25*alpha[8]*fUpwind[19]+0.25*alpha[5]*fUpwind[18]+0.25*alpha[0]*fUpwind[17]+0.25000000000000006*alpha[1]*fUpwind[16]+0.22360679774997902*alpha[7]*fUpwind[11]+0.22360679774997902*fUpwind[7]*alpha[11]+0.22360679774997902*alpha[3]*fUpwind[6]+0.22360679774997902*fUpwind[3]*alpha[6]; 
  Ghat[18] = 0.25*alpha[8]*fUpwind[23]+0.25000000000000006*alpha[4]*fUpwind[22]+0.25000000000000006*alpha[12]*fUpwind[21]+0.25000000000000006*alpha[1]*fUpwind[20]+0.25*alpha[9]*fUpwind[19]+0.25*alpha[0]*fUpwind[18]+0.25*alpha[5]*fUpwind[17]+0.25000000000000006*alpha[2]*fUpwind[16]+0.22360679774997902*alpha[6]*fUpwind[11]+0.22360679774997902*fUpwind[6]*alpha[11]+0.22360679774997902*alpha[3]*fUpwind[7]+0.22360679774997902*fUpwind[3]*alpha[7]; 
  Ghat[19] = 0.25*alpha[5]*fUpwind[23]+0.25000000000000006*alpha[2]*fUpwind[22]+0.25000000000000006*alpha[1]*fUpwind[21]+0.25000000000000006*alpha[12]*fUpwind[20]+0.25*alpha[0]*fUpwind[19]+0.25*alpha[9]*fUpwind[18]+0.25*alpha[8]*fUpwind[17]+0.25000000000000006*alpha[4]*fUpwind[16]+0.22360679774997902*alpha[11]*fUpwind[15]+0.22360679774997902*alpha[7]*fUpwind[14]+0.22360679774997902*alpha[6]*fUpwind[13]+0.22360679774997902*alpha[3]*fUpwind[10]; 
  Ghat[20] = 0.25000000000000006*alpha[4]*fUpwind[23]+0.25*alpha[8]*fUpwind[22]+0.25*alpha[9]*fUpwind[21]+0.25*alpha[0]*fUpwind[20]+0.25000000000000006*alpha[12]*fUpwind[19]+0.25000000000000006*alpha[1]*fUpwind[18]+0.25000000000000006*alpha[2]*fUpwind[17]+0.25*alpha[5]*fUpwind[16]+0.22360679774997896*alpha[3]*fUpwind[11]+0.22360679774997896*fUpwind[3]*alpha[11]+0.22360679774997896*alpha[6]*fUpwind[7]+0.22360679774997896*fUpwind[6]*alpha[7]; 
  Ghat[21] = 0.25000000000000006*alpha[2]*fUpwind[23]+0.25*alpha[5]*fUpwind[22]+0.25*alpha[0]*fUpwind[21]+0.25*alpha[9]*fUpwind[20]+0.25000000000000006*alpha[1]*fUpwind[19]+0.25000000000000006*alpha[12]*fUpwind[18]+0.25000000000000006*alpha[4]*fUpwind[17]+0.25*alpha[8]*fUpwind[16]+0.22360679774997896*alpha[7]*fUpwind[15]+0.22360679774997896*alpha[11]*fUpwind[14]+0.22360679774997896*alpha[3]*fUpwind[13]+0.22360679774997896*alpha[6]*fUpwind[10]; 
  Ghat[22] = 0.25000000000000006*alpha[1]*fUpwind[23]+0.25*alpha[0]*fUpwind[22]+0.25*alpha[5]*fUpwind[21]+0.25*alpha[8]*fUpwind[20]+0.25000000000000006*alpha[2]*fUpwind[19]+0.25000000000000006*alpha[4]*fUpwind[18]+0.25000000000000006*alpha[12]*fUpwind[17]+0.25*alpha[9]*fUpwind[16]+0.22360679774997896*alpha[6]*fUpwind[15]+0.22360679774997896*alpha[3]*fUpwind[14]+0.22360679774997896*alpha[11]*fUpwind[13]+0.22360679774997896*alpha[7]*fUpwind[10]; 
  Ghat[23] = 0.25*alpha[0]*fUpwind[23]+0.25000000000000006*alpha[1]*fUpwind[22]+0.25000000000000006*alpha[2]*fUpwind[21]+0.25000000000000006*alpha[4]*fUpwind[20]+0.25*alpha[5]*fUpwind[19]+0.25*alpha[8]*fUpwind[18]+0.25*alpha[9]*fUpwind[17]+0.25000000000000006*alpha[12]*fUpwind[16]+0.22360679774997902*alpha[3]*fUpwind[15]+0.22360679774997902*alpha[6]*fUpwind[14]+0.22360679774997902*alpha[7]*fUpwind[13]+0.22360679774997902*fUpwind[10]*alpha[11]; 
  Ghat[24] = 0.25000000000000006*alpha[11]*fUpwind[31]+0.25*alpha[7]*fUpwind[30]+0.25*alpha[6]*fUpwind[29]+0.25*alpha[5]*fUpwind[28]+0.25000000000000006*alpha[3]*fUpwind[27]+0.25000000000000006*alpha[2]*fUpwind[26]+0.25000000000000006*alpha[1]*fUpwind[25]+0.25*alpha[0]*fUpwind[24]+0.22360679774997896*alpha[12]*fUpwind[12]+0.22360679774997896*alpha[9]*fUpwind[9]+0.22360679774997896*alpha[8]*fUpwind[8]+0.22360679774997896*alpha[4]*fUpwind[4]; 
  Ghat[25] = 0.25*alpha[7]*fUpwind[31]+0.25000000000000006*alpha[11]*fUpwind[30]+0.25000000000000006*alpha[3]*fUpwind[29]+0.25000000000000006*alpha[2]*fUpwind[28]+0.25*alpha[6]*fUpwind[27]+0.25*alpha[5]*fUpwind[26]+0.25*alpha[0]*fUpwind[25]+0.25000000000000006*alpha[1]*fUpwind[24]+0.22360679774997902*alpha[9]*fUpwind[12]+0.22360679774997902*fUpwind[9]*alpha[12]+0.22360679774997902*alpha[4]*fUpwind[8]+0.22360679774997902*fUpwind[4]*alpha[8]; 
  Ghat[26] = 0.25*alpha[6]*fUpwind[31]+0.25000000000000006*alpha[3]*fUpwind[30]+0.25000000000000006*alpha[11]*fUpwind[29]+0.25000000000000006*alpha[1]*fUpwind[28]+0.25*alpha[7]*fUpwind[27]+0.25*alpha[0]*fUpwind[26]+0.25*alpha[5]*fUpwind[25]+0.25000000000000006*alpha[2]*fUpwind[24]+0.22360679774997902*alpha[8]*fUpwind[12]+0.22360679774997902*fUpwind[8]*alpha[12]+0.22360679774997902*alpha[4]*fUpwind[9]+0.22360679774997902*fUpwind[4]*alpha[9]; 
  Ghat[27] = 0.25*alpha[5]*fUpwind[31]+0.25000000000000006*alpha[2]*fUpwind[30]+0.25000000000000006*alpha[1]*fUpwind[29]+0.25000000000000006*alpha[11]*fUpwind[28]+0.25*alpha[0]*fUpwind[27]+0.25*alpha[7]*fUpwind[26]+0.25*alpha[6]*fUpwind[25]+0.25000000000000006*alpha[3]*fUpwind[24]+0.22360679774997902*alpha[12]*fUpwind[15]+0.22360679774997902*alpha[9]*fUpwind[14]+0.22360679774997902*alpha[8]*fUpwind[13]+0.22360679774997902*alpha[4]*fUpwind[10]; 
  Ghat[28] = 0.25000000000000006*alpha[3]*fUpwind[31]+0.25*alpha[6]*fUpwind[30]+0.25*alpha[7]*fUpwind[29]+0.25*alpha[0]*fUpwind[28]+0.25000000000000006*alpha[11]*fUpwind[27]+0.25000000000000006*alpha[1]*fUpwind[26]+0.25000000000000006*alpha[2]*fUpwind[25]+0.25*alpha[5]*fUpwind[24]+0.22360679774997896*alpha[4]*fUpwind[12]+0.22360679774997896*fUpwind[4]*alpha[12]+0.22360679774997896*alpha[8]*fUpwind[9]+0.22360679774997896*fUpwind[8]*alpha[9]; 
  Ghat[29] = 0.25000000000000006*alpha[2]*fUpwind[31]+0.25*alpha[5]*fUpwind[30]+0.25*alpha[0]*fUpwind[29]+0.25*alpha[7]*fUpwind[28]+0.25000000000000006*alpha[1]*fUpwind[27]+0.25000000000000006*alpha[11]*fUpwind[26]+0.25000000000000006*alpha[3]*fUpwind[25]+0.25*alpha[6]*fUpwind[24]+0.22360679774997896*alpha[9]*fUpwind[15]+0.22360679774997896*alpha[12]*fUpwind[14]+0.22360679774997896*alpha[4]*fUpwind[13]+0.22360679774997896*alpha[8]*fUpwind[10]; 
  Ghat[30] = 0.25000000000000006*alpha[1]*fUpwind[31]+0.25*alpha[0]*fUpwind[30]+0.25*alpha[5]*fUpwind[29]+0.25*alpha[6]*fUpwind[28]+0.25000000000000006*alpha[2]*fUpwind[27]+0.25000000000000006*alpha[3]*fUpwind[26]+0.25000000000000006*alpha[11]*fUpwind[25]+0.25*alpha[7]*fUpwind[24]+0.22360679774997896*alpha[8]*fUpwind[15]+0.22360679774997896*alpha[4]*fUpwind[14]+0.22360679774997896*alpha[12]*fUpwind[13]+0.22360679774997896*alpha[9]*fUpwind[10]; 
  Ghat[31] = 0.25*alpha[0]*fUpwind[31]+0.25000000000000006*alpha[1]*fUpwind[30]+0.25000000000000006*alpha[2]*fUpwind[29]+0.25000000000000006*alpha[3]*fUpwind[28]+0.25*alpha[5]*fUpwind[27]+0.25*alpha[6]*fUpwind[26]+0.25*alpha[7]*fUpwind[25]+0.25000000000000006*alpha[11]*fUpwind[24]+0.22360679774997902*alpha[4]*fUpwind[15]+0.22360679774997902*alpha[8]*fUpwind[14]+0.22360679774997902*alpha[9]*fUpwind[13]+0.22360679774997902*fUpwind[10]*alpha[12]; 

  out[0] += 0.7071067811865475*Ghat[0]*dv11; 
  out[1] += 0.7071067811865475*Ghat[1]*dv11; 
  out[2] += 0.7071067811865475*Ghat[2]*dv11; 
  out[3] += 0.7071067811865475*Ghat[3]*dv11; 
  out[4] += -(1.224744871391589*Ghat[0]*dv11); 
  out[5] += 0.7071067811865475*Ghat[4]*dv11; 
  out[6] += 0.7071067811865475*Ghat[5]*dv11; 
  out[7] += 0.7071067811865475*Ghat[6]*dv11; 
  out[8] += 0.7071067811865475*Ghat[7]*dv11; 
  out[9] += -(1.224744871391589*Ghat[1]*dv11); 
  out[10] += -(1.224744871391589*Ghat[2]*dv11); 
  out[11] += -(1.224744871391589*Ghat[3]*dv11); 
  out[12] += 0.7071067811865475*Ghat[8]*dv11; 
  out[13] += 0.7071067811865475*Ghat[9]*dv11; 
  out[14] += 0.7071067811865475*Ghat[10]*dv11; 
  out[15] += -(1.224744871391589*Ghat[4]*dv11); 
  out[16] += 0.7071067811865475*Ghat[11]*dv11; 
  out[17] += -(1.224744871391589*Ghat[5]*dv11); 
  out[18] += -(1.224744871391589*Ghat[6]*dv11); 
  out[19] += -(1.224744871391589*Ghat[7]*dv11); 
  out[20] += 0.7071067811865475*Ghat[12]*dv11; 
  out[21] += 0.7071067811865475*Ghat[13]*dv11; 
  out[22] += 0.7071067811865475*Ghat[14]*dv11; 
  out[23] += -(1.224744871391589*Ghat[8]*dv11); 
  out[24] += -(1.224744871391589*Ghat[9]*dv11); 
  out[25] += -(1.224744871391589*Ghat[10]*dv11); 
  out[26] += -(1.224744871391589*Ghat[11]*dv11); 
  out[27] += 0.7071067811865475*Ghat[15]*dv11; 
  out[28] += -(1.224744871391589*Ghat[12]*dv11); 
  out[29] += -(1.224744871391589*Ghat[13]*dv11); 
  out[30] += -(1.224744871391589*Ghat[14]*dv11); 
  out[31] += -(1.224744871391589*Ghat[15]*dv11); 
  out[32] += 0.7071067811865475*Ghat[16]*dv11; 
  out[33] += 0.7071067811865475*Ghat[17]*dv11; 
  out[34] += 0.7071067811865475*Ghat[18]*dv11; 
  out[35] += -(1.224744871391589*Ghat[16]*dv11); 
  out[36] += 0.7071067811865475*Ghat[19]*dv11; 
  out[37] += 0.7071067811865475*Ghat[20]*dv11; 
  out[38] += -(1.224744871391589*Ghat[17]*dv11); 
  out[39] += -(1.224744871391589*Ghat[18]*dv11); 
  out[40] += 0.7071067811865475*Ghat[21]*dv11; 
  out[41] += 0.7071067811865475*Ghat[22]*dv11; 
  out[42] += -(1.224744871391589*Ghat[19]*dv11); 
  out[43] += -(1.224744871391589*Ghat[20]*dv11); 
  out[44] += 0.7071067811865475*Ghat[23]*dv11; 
  out[45] += -(1.224744871391589*Ghat[21]*dv11); 
  out[46] += -(1.224744871391589*Ghat[22]*dv11); 
  out[47] += -(1.224744871391589*Ghat[23]*dv11); 
  out[48] += 1.5811388300841895*Ghat[0]*dv11; 
  out[49] += 1.5811388300841898*Ghat[1]*dv11; 
  out[50] += 1.5811388300841898*Ghat[2]*dv11; 
  out[51] += 1.5811388300841898*Ghat[3]*dv11; 
  out[52] += 1.5811388300841898*Ghat[4]*dv11; 
  out[53] += 1.5811388300841895*Ghat[5]*dv11; 
  out[54] += 1.5811388300841895*Ghat[6]*dv11; 
  out[55] += 1.5811388300841895*Ghat[7]*dv11; 
  out[56] += 1.5811388300841895*Ghat[8]*dv11; 
  out[57] += 1.5811388300841895*Ghat[9]*dv11; 
  out[58] += 1.5811388300841895*Ghat[10]*dv11; 
  out[59] += 1.5811388300841898*Ghat[11]*dv11; 
  out[60] += 1.5811388300841898*Ghat[12]*dv11; 
  out[61] += 1.5811388300841898*Ghat[13]*dv11; 
  out[62] += 1.5811388300841898*Ghat[14]*dv11; 
  out[63] += 1.5811388300841895*Ghat[15]*dv11; 
  out[64] += 0.7071067811865475*Ghat[24]*dv11; 
  out[65] += 0.7071067811865475*Ghat[25]*dv11; 
  out[66] += 0.7071067811865475*Ghat[26]*dv11; 
  out[67] += 0.7071067811865475*Ghat[27]*dv11; 
  out[68] += -(1.224744871391589*Ghat[24]*dv11); 
  out[69] += 0.7071067811865475*Ghat[28]*dv11; 
  out[70] += 0.7071067811865475*Ghat[29]*dv11; 
  out[71] += 0.7071067811865475*Ghat[30]*dv11; 
  out[72] += -(1.224744871391589*Ghat[25]*dv11); 
  out[73] += -(1.224744871391589*Ghat[26]*dv11); 
  out[74] += -(1.224744871391589*Ghat[27]*dv11); 
  out[75] += 0.7071067811865475*Ghat[31]*dv11; 
  out[76] += -(1.224744871391589*Ghat[28]*dv11); 
  out[77] += -(1.224744871391589*Ghat[29]*dv11); 
  out[78] += -(1.224744871391589*Ghat[30]*dv11); 
  out[79] += -(1.224744871391589*Ghat[31]*dv11); 

  } 
  return 0.;

} 
