#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_ser_5x_p1_surfx3_quad.h> 
#include <gkyl_basis_ser_5x_p1_upwind.h> 
GKYL_CU_DH void vlasov_poisson_surfvx_2x3v_ser_p1(const double *w, const double *dxv, const double *fac_phi, const double *vecA, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // fac_phi:   potential (scaled by appropriate factors).
  // vecA:      vector potential (scaled by appropriate factors). Unused in pure Vlasov-Poisson. 
  // fl/fc/fr:  Input Distribution function in left/center/right cells 
  // out:       Output distribution function in center cell 
  const double dv10 = 2/dxv[2]; 
  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv2 = dxv[3], wv2 = w[3]; 
  const double dv3 = dxv[4], wv3 = w[4]; 
  const double *phi = &fac_phi[0]; 
  const double dx10 = 2/dxv[0]; 
  const double dx11 = 2/dxv[1]; 
  double alpha[16] = {0.0}; 

  alpha[0] = -3.464101615137754*phi[1]*dx10; 
  alpha[2] = -3.464101615137754*phi[3]*dx10; 

  double fUpwindQuad_l[16] = {0.0};
  double fUpwindQuad_r[16] = {0.0};
  double fUpwind_l[16] = {0.0};
  double fUpwind_r[16] = {0.0};
  double Ghat_l[16] = {0.0}; 
  double Ghat_r[16] = {0.0}; 

  if (alpha[0]-alpha[2] > 0) { 
    fUpwindQuad_l[0] = ser_5x_p1_surfx3_quad_0_r(fl); 
    fUpwindQuad_r[0] = ser_5x_p1_surfx3_quad_0_r(fc); 
    fUpwindQuad_l[4] = ser_5x_p1_surfx3_quad_4_r(fl); 
    fUpwindQuad_r[4] = ser_5x_p1_surfx3_quad_4_r(fc); 
    fUpwindQuad_l[8] = ser_5x_p1_surfx3_quad_8_r(fl); 
    fUpwindQuad_r[8] = ser_5x_p1_surfx3_quad_8_r(fc); 
    fUpwindQuad_l[12] = ser_5x_p1_surfx3_quad_12_r(fl); 
    fUpwindQuad_r[12] = ser_5x_p1_surfx3_quad_12_r(fc); 
  } else { 
    fUpwindQuad_l[0] = ser_5x_p1_surfx3_quad_0_l(fc); 
    fUpwindQuad_r[0] = ser_5x_p1_surfx3_quad_0_l(fr); 
    fUpwindQuad_l[4] = ser_5x_p1_surfx3_quad_4_l(fc); 
    fUpwindQuad_r[4] = ser_5x_p1_surfx3_quad_4_l(fr); 
    fUpwindQuad_l[8] = ser_5x_p1_surfx3_quad_8_l(fc); 
    fUpwindQuad_r[8] = ser_5x_p1_surfx3_quad_8_l(fr); 
    fUpwindQuad_l[12] = ser_5x_p1_surfx3_quad_12_l(fc); 
    fUpwindQuad_r[12] = ser_5x_p1_surfx3_quad_12_l(fr); 
  } 
  if (alpha[0]-alpha[2] > 0) { 
    fUpwindQuad_l[1] = ser_5x_p1_surfx3_quad_1_r(fl); 
    fUpwindQuad_r[1] = ser_5x_p1_surfx3_quad_1_r(fc); 
    fUpwindQuad_l[5] = ser_5x_p1_surfx3_quad_5_r(fl); 
    fUpwindQuad_r[5] = ser_5x_p1_surfx3_quad_5_r(fc); 
    fUpwindQuad_l[9] = ser_5x_p1_surfx3_quad_9_r(fl); 
    fUpwindQuad_r[9] = ser_5x_p1_surfx3_quad_9_r(fc); 
    fUpwindQuad_l[13] = ser_5x_p1_surfx3_quad_13_r(fl); 
    fUpwindQuad_r[13] = ser_5x_p1_surfx3_quad_13_r(fc); 
  } else { 
    fUpwindQuad_l[1] = ser_5x_p1_surfx3_quad_1_l(fc); 
    fUpwindQuad_r[1] = ser_5x_p1_surfx3_quad_1_l(fr); 
    fUpwindQuad_l[5] = ser_5x_p1_surfx3_quad_5_l(fc); 
    fUpwindQuad_r[5] = ser_5x_p1_surfx3_quad_5_l(fr); 
    fUpwindQuad_l[9] = ser_5x_p1_surfx3_quad_9_l(fc); 
    fUpwindQuad_r[9] = ser_5x_p1_surfx3_quad_9_l(fr); 
    fUpwindQuad_l[13] = ser_5x_p1_surfx3_quad_13_l(fc); 
    fUpwindQuad_r[13] = ser_5x_p1_surfx3_quad_13_l(fr); 
  } 
  if (alpha[2]+alpha[0] > 0) { 
    fUpwindQuad_l[2] = ser_5x_p1_surfx3_quad_2_r(fl); 
    fUpwindQuad_r[2] = ser_5x_p1_surfx3_quad_2_r(fc); 
    fUpwindQuad_l[6] = ser_5x_p1_surfx3_quad_6_r(fl); 
    fUpwindQuad_r[6] = ser_5x_p1_surfx3_quad_6_r(fc); 
    fUpwindQuad_l[10] = ser_5x_p1_surfx3_quad_10_r(fl); 
    fUpwindQuad_r[10] = ser_5x_p1_surfx3_quad_10_r(fc); 
    fUpwindQuad_l[14] = ser_5x_p1_surfx3_quad_14_r(fl); 
    fUpwindQuad_r[14] = ser_5x_p1_surfx3_quad_14_r(fc); 
  } else { 
    fUpwindQuad_l[2] = ser_5x_p1_surfx3_quad_2_l(fc); 
    fUpwindQuad_r[2] = ser_5x_p1_surfx3_quad_2_l(fr); 
    fUpwindQuad_l[6] = ser_5x_p1_surfx3_quad_6_l(fc); 
    fUpwindQuad_r[6] = ser_5x_p1_surfx3_quad_6_l(fr); 
    fUpwindQuad_l[10] = ser_5x_p1_surfx3_quad_10_l(fc); 
    fUpwindQuad_r[10] = ser_5x_p1_surfx3_quad_10_l(fr); 
    fUpwindQuad_l[14] = ser_5x_p1_surfx3_quad_14_l(fc); 
    fUpwindQuad_r[14] = ser_5x_p1_surfx3_quad_14_l(fr); 
  } 
  if (alpha[2]+alpha[0] > 0) { 
    fUpwindQuad_l[3] = ser_5x_p1_surfx3_quad_3_r(fl); 
    fUpwindQuad_r[3] = ser_5x_p1_surfx3_quad_3_r(fc); 
    fUpwindQuad_l[7] = ser_5x_p1_surfx3_quad_7_r(fl); 
    fUpwindQuad_r[7] = ser_5x_p1_surfx3_quad_7_r(fc); 
    fUpwindQuad_l[11] = ser_5x_p1_surfx3_quad_11_r(fl); 
    fUpwindQuad_r[11] = ser_5x_p1_surfx3_quad_11_r(fc); 
    fUpwindQuad_l[15] = ser_5x_p1_surfx3_quad_15_r(fl); 
    fUpwindQuad_r[15] = ser_5x_p1_surfx3_quad_15_r(fc); 
  } else { 
    fUpwindQuad_l[3] = ser_5x_p1_surfx3_quad_3_l(fc); 
    fUpwindQuad_r[3] = ser_5x_p1_surfx3_quad_3_l(fr); 
    fUpwindQuad_l[7] = ser_5x_p1_surfx3_quad_7_l(fc); 
    fUpwindQuad_r[7] = ser_5x_p1_surfx3_quad_7_l(fr); 
    fUpwindQuad_l[11] = ser_5x_p1_surfx3_quad_11_l(fc); 
    fUpwindQuad_r[11] = ser_5x_p1_surfx3_quad_11_l(fr); 
    fUpwindQuad_l[15] = ser_5x_p1_surfx3_quad_15_l(fc); 
    fUpwindQuad_r[15] = ser_5x_p1_surfx3_quad_15_l(fr); 
  } 

  // Project nodal basis back onto modal basis. 
  ser_5x_p1_upwind(fUpwindQuad_l, fUpwind_l); 
  ser_5x_p1_upwind(fUpwindQuad_r, fUpwind_r); 

  Ghat_l[0] += 0.25*(alpha[2]*fUpwind_l[2]+alpha[0]*fUpwind_l[0]); 
  Ghat_l[1] += 0.25*(alpha[2]*fUpwind_l[5]+alpha[0]*fUpwind_l[1]); 
  Ghat_l[2] += 0.25*(alpha[0]*fUpwind_l[2]+fUpwind_l[0]*alpha[2]); 
  Ghat_l[3] += 0.25*(alpha[2]*fUpwind_l[7]+alpha[0]*fUpwind_l[3]); 
  Ghat_l[4] += 0.25*(alpha[2]*fUpwind_l[9]+alpha[0]*fUpwind_l[4]); 
  Ghat_l[5] += 0.25*(alpha[0]*fUpwind_l[5]+fUpwind_l[1]*alpha[2]); 
  Ghat_l[6] += 0.25*(alpha[2]*fUpwind_l[11]+alpha[0]*fUpwind_l[6]); 
  Ghat_l[7] += 0.25*(alpha[0]*fUpwind_l[7]+alpha[2]*fUpwind_l[3]); 
  Ghat_l[8] += 0.25*(alpha[2]*fUpwind_l[12]+alpha[0]*fUpwind_l[8]); 
  Ghat_l[9] += 0.25*(alpha[0]*fUpwind_l[9]+alpha[2]*fUpwind_l[4]); 
  Ghat_l[10] += 0.25*(alpha[2]*fUpwind_l[14]+alpha[0]*fUpwind_l[10]); 
  Ghat_l[11] += 0.25*(alpha[0]*fUpwind_l[11]+alpha[2]*fUpwind_l[6]); 
  Ghat_l[12] += 0.25*(alpha[0]*fUpwind_l[12]+alpha[2]*fUpwind_l[8]); 
  Ghat_l[13] += 0.25*(alpha[2]*fUpwind_l[15]+alpha[0]*fUpwind_l[13]); 
  Ghat_l[14] += 0.25*(alpha[0]*fUpwind_l[14]+alpha[2]*fUpwind_l[10]); 
  Ghat_l[15] += 0.25*(alpha[0]*fUpwind_l[15]+alpha[2]*fUpwind_l[13]); 

  Ghat_r[0] += 0.25*(alpha[2]*fUpwind_r[2]+alpha[0]*fUpwind_r[0]); 
  Ghat_r[1] += 0.25*(alpha[2]*fUpwind_r[5]+alpha[0]*fUpwind_r[1]); 
  Ghat_r[2] += 0.25*(alpha[0]*fUpwind_r[2]+fUpwind_r[0]*alpha[2]); 
  Ghat_r[3] += 0.25*(alpha[2]*fUpwind_r[7]+alpha[0]*fUpwind_r[3]); 
  Ghat_r[4] += 0.25*(alpha[2]*fUpwind_r[9]+alpha[0]*fUpwind_r[4]); 
  Ghat_r[5] += 0.25*(alpha[0]*fUpwind_r[5]+fUpwind_r[1]*alpha[2]); 
  Ghat_r[6] += 0.25*(alpha[2]*fUpwind_r[11]+alpha[0]*fUpwind_r[6]); 
  Ghat_r[7] += 0.25*(alpha[0]*fUpwind_r[7]+alpha[2]*fUpwind_r[3]); 
  Ghat_r[8] += 0.25*(alpha[2]*fUpwind_r[12]+alpha[0]*fUpwind_r[8]); 
  Ghat_r[9] += 0.25*(alpha[0]*fUpwind_r[9]+alpha[2]*fUpwind_r[4]); 
  Ghat_r[10] += 0.25*(alpha[2]*fUpwind_r[14]+alpha[0]*fUpwind_r[10]); 
  Ghat_r[11] += 0.25*(alpha[0]*fUpwind_r[11]+alpha[2]*fUpwind_r[6]); 
  Ghat_r[12] += 0.25*(alpha[0]*fUpwind_r[12]+alpha[2]*fUpwind_r[8]); 
  Ghat_r[13] += 0.25*(alpha[2]*fUpwind_r[15]+alpha[0]*fUpwind_r[13]); 
  Ghat_r[14] += 0.25*(alpha[0]*fUpwind_r[14]+alpha[2]*fUpwind_r[10]); 
  Ghat_r[15] += 0.25*(alpha[0]*fUpwind_r[15]+alpha[2]*fUpwind_r[13]); 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv10; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv10; 
  out[2] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dv10; 
  out[3] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv10; 
  out[4] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dv10; 
  out[5] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*dv10; 
  out[6] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*dv10; 
  out[7] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv10; 
  out[8] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dv10; 
  out[9] += (0.7071067811865475*Ghat_l[6]-0.7071067811865475*Ghat_r[6])*dv10; 
  out[10] += (0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*dv10; 
  out[11] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dv10; 
  out[12] += (0.7071067811865475*Ghat_l[8]-0.7071067811865475*Ghat_r[8])*dv10; 
  out[13] += (0.7071067811865475*Ghat_l[9]-0.7071067811865475*Ghat_r[9])*dv10; 
  out[14] += -1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dv10; 
  out[15] += (0.7071067811865475*Ghat_l[10]-0.7071067811865475*Ghat_r[10])*dv10; 
  out[16] += -1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dv10; 
  out[17] += (0.7071067811865475*Ghat_l[11]-0.7071067811865475*Ghat_r[11])*dv10; 
  out[18] += -1.224744871391589*(Ghat_r[6]+Ghat_l[6])*dv10; 
  out[19] += -1.224744871391589*(Ghat_r[7]+Ghat_l[7])*dv10; 
  out[20] += (0.7071067811865475*Ghat_l[12]-0.7071067811865475*Ghat_r[12])*dv10; 
  out[21] += -1.224744871391589*(Ghat_r[8]+Ghat_l[8])*dv10; 
  out[22] += -1.224744871391589*(Ghat_r[9]+Ghat_l[9])*dv10; 
  out[23] += (0.7071067811865475*Ghat_l[13]-0.7071067811865475*Ghat_r[13])*dv10; 
  out[24] += (0.7071067811865475*Ghat_l[14]-0.7071067811865475*Ghat_r[14])*dv10; 
  out[25] += -1.224744871391589*(Ghat_r[10]+Ghat_l[10])*dv10; 
  out[26] += -1.224744871391589*(Ghat_r[11]+Ghat_l[11])*dv10; 
  out[27] += -1.224744871391589*(Ghat_r[12]+Ghat_l[12])*dv10; 
  out[28] += (0.7071067811865475*Ghat_l[15]-0.7071067811865475*Ghat_r[15])*dv10; 
  out[29] += -1.224744871391589*(Ghat_r[13]+Ghat_l[13])*dv10; 
  out[30] += -1.224744871391589*(Ghat_r[14]+Ghat_l[14])*dv10; 
  out[31] += -1.224744871391589*(Ghat_r[15]+Ghat_l[15])*dv10; 

} 
