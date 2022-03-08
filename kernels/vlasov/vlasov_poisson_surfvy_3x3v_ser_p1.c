#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_ser_6x_p1_surfx5_quad.h> 
#include <gkyl_basis_ser_6x_p1_upwind.h> 
GKYL_CU_DH void vlasov_poisson_surfvy_3x3v_ser_p1(const double *w, const double *dxv, const double *fac_phi, const double *vecA, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // fac_phi:   potential (scaled by appropriate factors).
  // vecA:      vector potential (scaled by appropriate factors). Unused in pure Vlasov-Poisson. 
  // fl/fc/fr:  Input Distribution function in left/center/right cells 
  // out:       Output distribution function in center cell 
  const double dv11 = 2/dxv[4]; 
  const double dv1 = dxv[3], wv1 = w[3]; 
  const double dv2 = dxv[4], wv2 = w[4]; 
  const double dv3 = dxv[5], wv3 = w[5]; 
  const double *phi = &fac_phi[0]; 
  const double dx10 = 2/dxv[0]; 
  const double dx11 = 2/dxv[1]; 
  const double dx12 = 2/dxv[2]; 
  double alpha[32] = {0.0}; 

  alpha[0] = -3.464101615137754*phi[2]*dx11; 
  alpha[1] = -3.464101615137754*phi[4]*dx11; 
  alpha[3] = -3.464101615137754*phi[6]*dx11; 
  alpha[7] = -3.464101615137754*phi[7]*dx11; 

  double fUpwindQuad_l[32] = {0.0};
  double fUpwindQuad_r[32] = {0.0};
  double fUpwind_l[32] = {0.0};
  double fUpwind_r[32] = {0.0};
  double Ghat_l[32] = {0.0}; 
  double Ghat_r[32] = {0.0}; 

  if (alpha[7]-alpha[3]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad_l[0] = ser_6x_p1_surfx5_quad_0_r(fl); 
    fUpwindQuad_r[0] = ser_6x_p1_surfx5_quad_0_r(fc); 
    fUpwindQuad_l[8] = ser_6x_p1_surfx5_quad_8_r(fl); 
    fUpwindQuad_r[8] = ser_6x_p1_surfx5_quad_8_r(fc); 
    fUpwindQuad_l[16] = ser_6x_p1_surfx5_quad_16_r(fl); 
    fUpwindQuad_r[16] = ser_6x_p1_surfx5_quad_16_r(fc); 
    fUpwindQuad_l[24] = ser_6x_p1_surfx5_quad_24_r(fl); 
    fUpwindQuad_r[24] = ser_6x_p1_surfx5_quad_24_r(fc); 
  } else { 
    fUpwindQuad_l[0] = ser_6x_p1_surfx5_quad_0_l(fc); 
    fUpwindQuad_r[0] = ser_6x_p1_surfx5_quad_0_l(fr); 
    fUpwindQuad_l[8] = ser_6x_p1_surfx5_quad_8_l(fc); 
    fUpwindQuad_r[8] = ser_6x_p1_surfx5_quad_8_l(fr); 
    fUpwindQuad_l[16] = ser_6x_p1_surfx5_quad_16_l(fc); 
    fUpwindQuad_r[16] = ser_6x_p1_surfx5_quad_16_l(fr); 
    fUpwindQuad_l[24] = ser_6x_p1_surfx5_quad_24_l(fc); 
    fUpwindQuad_r[24] = ser_6x_p1_surfx5_quad_24_l(fr); 
  } 
  if ((-alpha[7])-alpha[3]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad_l[1] = ser_6x_p1_surfx5_quad_1_r(fl); 
    fUpwindQuad_r[1] = ser_6x_p1_surfx5_quad_1_r(fc); 
    fUpwindQuad_l[9] = ser_6x_p1_surfx5_quad_9_r(fl); 
    fUpwindQuad_r[9] = ser_6x_p1_surfx5_quad_9_r(fc); 
    fUpwindQuad_l[17] = ser_6x_p1_surfx5_quad_17_r(fl); 
    fUpwindQuad_r[17] = ser_6x_p1_surfx5_quad_17_r(fc); 
    fUpwindQuad_l[25] = ser_6x_p1_surfx5_quad_25_r(fl); 
    fUpwindQuad_r[25] = ser_6x_p1_surfx5_quad_25_r(fc); 
  } else { 
    fUpwindQuad_l[1] = ser_6x_p1_surfx5_quad_1_l(fc); 
    fUpwindQuad_r[1] = ser_6x_p1_surfx5_quad_1_l(fr); 
    fUpwindQuad_l[9] = ser_6x_p1_surfx5_quad_9_l(fc); 
    fUpwindQuad_r[9] = ser_6x_p1_surfx5_quad_9_l(fr); 
    fUpwindQuad_l[17] = ser_6x_p1_surfx5_quad_17_l(fc); 
    fUpwindQuad_r[17] = ser_6x_p1_surfx5_quad_17_l(fr); 
    fUpwindQuad_l[25] = ser_6x_p1_surfx5_quad_25_l(fc); 
    fUpwindQuad_r[25] = ser_6x_p1_surfx5_quad_25_l(fr); 
  } 
  if (alpha[7]-alpha[3]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad_l[2] = ser_6x_p1_surfx5_quad_2_r(fl); 
    fUpwindQuad_r[2] = ser_6x_p1_surfx5_quad_2_r(fc); 
    fUpwindQuad_l[10] = ser_6x_p1_surfx5_quad_10_r(fl); 
    fUpwindQuad_r[10] = ser_6x_p1_surfx5_quad_10_r(fc); 
    fUpwindQuad_l[18] = ser_6x_p1_surfx5_quad_18_r(fl); 
    fUpwindQuad_r[18] = ser_6x_p1_surfx5_quad_18_r(fc); 
    fUpwindQuad_l[26] = ser_6x_p1_surfx5_quad_26_r(fl); 
    fUpwindQuad_r[26] = ser_6x_p1_surfx5_quad_26_r(fc); 
  } else { 
    fUpwindQuad_l[2] = ser_6x_p1_surfx5_quad_2_l(fc); 
    fUpwindQuad_r[2] = ser_6x_p1_surfx5_quad_2_l(fr); 
    fUpwindQuad_l[10] = ser_6x_p1_surfx5_quad_10_l(fc); 
    fUpwindQuad_r[10] = ser_6x_p1_surfx5_quad_10_l(fr); 
    fUpwindQuad_l[18] = ser_6x_p1_surfx5_quad_18_l(fc); 
    fUpwindQuad_r[18] = ser_6x_p1_surfx5_quad_18_l(fr); 
    fUpwindQuad_l[26] = ser_6x_p1_surfx5_quad_26_l(fc); 
    fUpwindQuad_r[26] = ser_6x_p1_surfx5_quad_26_l(fr); 
  } 
  if ((-alpha[7])-alpha[3]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad_l[3] = ser_6x_p1_surfx5_quad_3_r(fl); 
    fUpwindQuad_r[3] = ser_6x_p1_surfx5_quad_3_r(fc); 
    fUpwindQuad_l[11] = ser_6x_p1_surfx5_quad_11_r(fl); 
    fUpwindQuad_r[11] = ser_6x_p1_surfx5_quad_11_r(fc); 
    fUpwindQuad_l[19] = ser_6x_p1_surfx5_quad_19_r(fl); 
    fUpwindQuad_r[19] = ser_6x_p1_surfx5_quad_19_r(fc); 
    fUpwindQuad_l[27] = ser_6x_p1_surfx5_quad_27_r(fl); 
    fUpwindQuad_r[27] = ser_6x_p1_surfx5_quad_27_r(fc); 
  } else { 
    fUpwindQuad_l[3] = ser_6x_p1_surfx5_quad_3_l(fc); 
    fUpwindQuad_r[3] = ser_6x_p1_surfx5_quad_3_l(fr); 
    fUpwindQuad_l[11] = ser_6x_p1_surfx5_quad_11_l(fc); 
    fUpwindQuad_r[11] = ser_6x_p1_surfx5_quad_11_l(fr); 
    fUpwindQuad_l[19] = ser_6x_p1_surfx5_quad_19_l(fc); 
    fUpwindQuad_r[19] = ser_6x_p1_surfx5_quad_19_l(fr); 
    fUpwindQuad_l[27] = ser_6x_p1_surfx5_quad_27_l(fc); 
    fUpwindQuad_r[27] = ser_6x_p1_surfx5_quad_27_l(fr); 
  } 
  if ((-alpha[7])+alpha[3]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad_l[4] = ser_6x_p1_surfx5_quad_4_r(fl); 
    fUpwindQuad_r[4] = ser_6x_p1_surfx5_quad_4_r(fc); 
    fUpwindQuad_l[12] = ser_6x_p1_surfx5_quad_12_r(fl); 
    fUpwindQuad_r[12] = ser_6x_p1_surfx5_quad_12_r(fc); 
    fUpwindQuad_l[20] = ser_6x_p1_surfx5_quad_20_r(fl); 
    fUpwindQuad_r[20] = ser_6x_p1_surfx5_quad_20_r(fc); 
    fUpwindQuad_l[28] = ser_6x_p1_surfx5_quad_28_r(fl); 
    fUpwindQuad_r[28] = ser_6x_p1_surfx5_quad_28_r(fc); 
  } else { 
    fUpwindQuad_l[4] = ser_6x_p1_surfx5_quad_4_l(fc); 
    fUpwindQuad_r[4] = ser_6x_p1_surfx5_quad_4_l(fr); 
    fUpwindQuad_l[12] = ser_6x_p1_surfx5_quad_12_l(fc); 
    fUpwindQuad_r[12] = ser_6x_p1_surfx5_quad_12_l(fr); 
    fUpwindQuad_l[20] = ser_6x_p1_surfx5_quad_20_l(fc); 
    fUpwindQuad_r[20] = ser_6x_p1_surfx5_quad_20_l(fr); 
    fUpwindQuad_l[28] = ser_6x_p1_surfx5_quad_28_l(fc); 
    fUpwindQuad_r[28] = ser_6x_p1_surfx5_quad_28_l(fr); 
  } 
  if (alpha[7]+alpha[3]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad_l[5] = ser_6x_p1_surfx5_quad_5_r(fl); 
    fUpwindQuad_r[5] = ser_6x_p1_surfx5_quad_5_r(fc); 
    fUpwindQuad_l[13] = ser_6x_p1_surfx5_quad_13_r(fl); 
    fUpwindQuad_r[13] = ser_6x_p1_surfx5_quad_13_r(fc); 
    fUpwindQuad_l[21] = ser_6x_p1_surfx5_quad_21_r(fl); 
    fUpwindQuad_r[21] = ser_6x_p1_surfx5_quad_21_r(fc); 
    fUpwindQuad_l[29] = ser_6x_p1_surfx5_quad_29_r(fl); 
    fUpwindQuad_r[29] = ser_6x_p1_surfx5_quad_29_r(fc); 
  } else { 
    fUpwindQuad_l[5] = ser_6x_p1_surfx5_quad_5_l(fc); 
    fUpwindQuad_r[5] = ser_6x_p1_surfx5_quad_5_l(fr); 
    fUpwindQuad_l[13] = ser_6x_p1_surfx5_quad_13_l(fc); 
    fUpwindQuad_r[13] = ser_6x_p1_surfx5_quad_13_l(fr); 
    fUpwindQuad_l[21] = ser_6x_p1_surfx5_quad_21_l(fc); 
    fUpwindQuad_r[21] = ser_6x_p1_surfx5_quad_21_l(fr); 
    fUpwindQuad_l[29] = ser_6x_p1_surfx5_quad_29_l(fc); 
    fUpwindQuad_r[29] = ser_6x_p1_surfx5_quad_29_l(fr); 
  } 
  if ((-alpha[7])+alpha[3]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad_l[6] = ser_6x_p1_surfx5_quad_6_r(fl); 
    fUpwindQuad_r[6] = ser_6x_p1_surfx5_quad_6_r(fc); 
    fUpwindQuad_l[14] = ser_6x_p1_surfx5_quad_14_r(fl); 
    fUpwindQuad_r[14] = ser_6x_p1_surfx5_quad_14_r(fc); 
    fUpwindQuad_l[22] = ser_6x_p1_surfx5_quad_22_r(fl); 
    fUpwindQuad_r[22] = ser_6x_p1_surfx5_quad_22_r(fc); 
    fUpwindQuad_l[30] = ser_6x_p1_surfx5_quad_30_r(fl); 
    fUpwindQuad_r[30] = ser_6x_p1_surfx5_quad_30_r(fc); 
  } else { 
    fUpwindQuad_l[6] = ser_6x_p1_surfx5_quad_6_l(fc); 
    fUpwindQuad_r[6] = ser_6x_p1_surfx5_quad_6_l(fr); 
    fUpwindQuad_l[14] = ser_6x_p1_surfx5_quad_14_l(fc); 
    fUpwindQuad_r[14] = ser_6x_p1_surfx5_quad_14_l(fr); 
    fUpwindQuad_l[22] = ser_6x_p1_surfx5_quad_22_l(fc); 
    fUpwindQuad_r[22] = ser_6x_p1_surfx5_quad_22_l(fr); 
    fUpwindQuad_l[30] = ser_6x_p1_surfx5_quad_30_l(fc); 
    fUpwindQuad_r[30] = ser_6x_p1_surfx5_quad_30_l(fr); 
  } 
  if (alpha[7]+alpha[3]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad_l[7] = ser_6x_p1_surfx5_quad_7_r(fl); 
    fUpwindQuad_r[7] = ser_6x_p1_surfx5_quad_7_r(fc); 
    fUpwindQuad_l[15] = ser_6x_p1_surfx5_quad_15_r(fl); 
    fUpwindQuad_r[15] = ser_6x_p1_surfx5_quad_15_r(fc); 
    fUpwindQuad_l[23] = ser_6x_p1_surfx5_quad_23_r(fl); 
    fUpwindQuad_r[23] = ser_6x_p1_surfx5_quad_23_r(fc); 
    fUpwindQuad_l[31] = ser_6x_p1_surfx5_quad_31_r(fl); 
    fUpwindQuad_r[31] = ser_6x_p1_surfx5_quad_31_r(fc); 
  } else { 
    fUpwindQuad_l[7] = ser_6x_p1_surfx5_quad_7_l(fc); 
    fUpwindQuad_r[7] = ser_6x_p1_surfx5_quad_7_l(fr); 
    fUpwindQuad_l[15] = ser_6x_p1_surfx5_quad_15_l(fc); 
    fUpwindQuad_r[15] = ser_6x_p1_surfx5_quad_15_l(fr); 
    fUpwindQuad_l[23] = ser_6x_p1_surfx5_quad_23_l(fc); 
    fUpwindQuad_r[23] = ser_6x_p1_surfx5_quad_23_l(fr); 
    fUpwindQuad_l[31] = ser_6x_p1_surfx5_quad_31_l(fc); 
    fUpwindQuad_r[31] = ser_6x_p1_surfx5_quad_31_l(fr); 
  } 

  // Project nodal basis back onto modal basis. 
  ser_6x_p1_upwind(fUpwindQuad_l, fUpwind_l); 
  ser_6x_p1_upwind(fUpwindQuad_r, fUpwind_r); 

  Ghat_l[0] += 0.1767766952966368*(alpha[7]*fUpwind_l[7]+alpha[3]*fUpwind_l[3]+alpha[1]*fUpwind_l[1]+alpha[0]*fUpwind_l[0]); 
  Ghat_l[1] += 0.1767766952966368*(alpha[3]*fUpwind_l[7]+fUpwind_l[3]*alpha[7]+alpha[0]*fUpwind_l[1]+fUpwind_l[0]*alpha[1]); 
  Ghat_l[2] += 0.1767766952966368*(alpha[7]*fUpwind_l[16]+alpha[3]*fUpwind_l[8]+alpha[1]*fUpwind_l[6]+alpha[0]*fUpwind_l[2]); 
  Ghat_l[3] += 0.1767766952966368*(alpha[1]*fUpwind_l[7]+fUpwind_l[1]*alpha[7]+alpha[0]*fUpwind_l[3]+fUpwind_l[0]*alpha[3]); 
  Ghat_l[4] += 0.1767766952966368*(alpha[7]*fUpwind_l[18]+alpha[3]*fUpwind_l[11]+alpha[1]*fUpwind_l[9]+alpha[0]*fUpwind_l[4]); 
  Ghat_l[5] += 0.1767766952966368*(alpha[7]*fUpwind_l[21]+alpha[3]*fUpwind_l[14]+alpha[1]*fUpwind_l[12]+alpha[0]*fUpwind_l[5]); 
  Ghat_l[6] += 0.1767766952966368*(alpha[3]*fUpwind_l[16]+alpha[7]*fUpwind_l[8]+alpha[0]*fUpwind_l[6]+alpha[1]*fUpwind_l[2]); 
  Ghat_l[7] += 0.1767766952966368*(alpha[0]*fUpwind_l[7]+fUpwind_l[0]*alpha[7]+alpha[1]*fUpwind_l[3]+fUpwind_l[1]*alpha[3]); 
  Ghat_l[8] += 0.1767766952966368*(alpha[1]*fUpwind_l[16]+alpha[0]*fUpwind_l[8]+fUpwind_l[6]*alpha[7]+fUpwind_l[2]*alpha[3]); 
  Ghat_l[9] += 0.1767766952966368*(alpha[3]*fUpwind_l[18]+alpha[7]*fUpwind_l[11]+alpha[0]*fUpwind_l[9]+alpha[1]*fUpwind_l[4]); 
  Ghat_l[10] += 0.1767766952966368*(alpha[7]*fUpwind_l[26]+alpha[3]*fUpwind_l[19]+alpha[1]*fUpwind_l[17]+alpha[0]*fUpwind_l[10]); 
  Ghat_l[11] += 0.1767766952966368*(alpha[1]*fUpwind_l[18]+alpha[0]*fUpwind_l[11]+alpha[7]*fUpwind_l[9]+alpha[3]*fUpwind_l[4]); 
  Ghat_l[12] += 0.1767766952966368*(alpha[3]*fUpwind_l[21]+alpha[7]*fUpwind_l[14]+alpha[0]*fUpwind_l[12]+alpha[1]*fUpwind_l[5]); 
  Ghat_l[13] += 0.1767766952966368*(alpha[7]*fUpwind_l[27]+alpha[3]*fUpwind_l[22]+alpha[1]*fUpwind_l[20]+alpha[0]*fUpwind_l[13]); 
  Ghat_l[14] += 0.1767766952966368*(alpha[1]*fUpwind_l[21]+alpha[0]*fUpwind_l[14]+alpha[7]*fUpwind_l[12]+alpha[3]*fUpwind_l[5]); 
  Ghat_l[15] += 0.1767766952966368*(alpha[7]*fUpwind_l[29]+alpha[3]*fUpwind_l[25]+alpha[1]*fUpwind_l[23]+alpha[0]*fUpwind_l[15]); 
  Ghat_l[16] += 0.1767766952966368*(alpha[0]*fUpwind_l[16]+alpha[1]*fUpwind_l[8]+fUpwind_l[2]*alpha[7]+alpha[3]*fUpwind_l[6]); 
  Ghat_l[17] += 0.1767766952966368*(alpha[3]*fUpwind_l[26]+alpha[7]*fUpwind_l[19]+alpha[0]*fUpwind_l[17]+alpha[1]*fUpwind_l[10]); 
  Ghat_l[18] += 0.1767766952966368*(alpha[0]*fUpwind_l[18]+alpha[1]*fUpwind_l[11]+alpha[3]*fUpwind_l[9]+fUpwind_l[4]*alpha[7]); 
  Ghat_l[19] += 0.1767766952966368*(alpha[1]*fUpwind_l[26]+alpha[0]*fUpwind_l[19]+alpha[7]*fUpwind_l[17]+alpha[3]*fUpwind_l[10]); 
  Ghat_l[20] += 0.1767766952966368*(alpha[3]*fUpwind_l[27]+alpha[7]*fUpwind_l[22]+alpha[0]*fUpwind_l[20]+alpha[1]*fUpwind_l[13]); 
  Ghat_l[21] += 0.1767766952966368*(alpha[0]*fUpwind_l[21]+alpha[1]*fUpwind_l[14]+alpha[3]*fUpwind_l[12]+fUpwind_l[5]*alpha[7]); 
  Ghat_l[22] += 0.1767766952966368*(alpha[1]*fUpwind_l[27]+alpha[0]*fUpwind_l[22]+alpha[7]*fUpwind_l[20]+alpha[3]*fUpwind_l[13]); 
  Ghat_l[23] += 0.1767766952966368*(alpha[3]*fUpwind_l[29]+alpha[7]*fUpwind_l[25]+alpha[0]*fUpwind_l[23]+alpha[1]*fUpwind_l[15]); 
  Ghat_l[24] += 0.1767766952966368*(alpha[7]*fUpwind_l[31]+alpha[3]*fUpwind_l[30]+alpha[1]*fUpwind_l[28]+alpha[0]*fUpwind_l[24]); 
  Ghat_l[25] += 0.1767766952966368*(alpha[1]*fUpwind_l[29]+alpha[0]*fUpwind_l[25]+alpha[7]*fUpwind_l[23]+alpha[3]*fUpwind_l[15]); 
  Ghat_l[26] += 0.1767766952966368*(alpha[0]*fUpwind_l[26]+alpha[1]*fUpwind_l[19]+alpha[3]*fUpwind_l[17]+alpha[7]*fUpwind_l[10]); 
  Ghat_l[27] += 0.1767766952966368*(alpha[0]*fUpwind_l[27]+alpha[1]*fUpwind_l[22]+alpha[3]*fUpwind_l[20]+alpha[7]*fUpwind_l[13]); 
  Ghat_l[28] += 0.1767766952966368*(alpha[3]*fUpwind_l[31]+alpha[7]*fUpwind_l[30]+alpha[0]*fUpwind_l[28]+alpha[1]*fUpwind_l[24]); 
  Ghat_l[29] += 0.1767766952966368*(alpha[0]*fUpwind_l[29]+alpha[1]*fUpwind_l[25]+alpha[3]*fUpwind_l[23]+alpha[7]*fUpwind_l[15]); 
  Ghat_l[30] += 0.1767766952966368*(alpha[1]*fUpwind_l[31]+alpha[0]*fUpwind_l[30]+alpha[7]*fUpwind_l[28]+alpha[3]*fUpwind_l[24]); 
  Ghat_l[31] += 0.1767766952966368*(alpha[0]*fUpwind_l[31]+alpha[1]*fUpwind_l[30]+alpha[3]*fUpwind_l[28]+alpha[7]*fUpwind_l[24]); 

  Ghat_r[0] += 0.1767766952966368*(alpha[7]*fUpwind_r[7]+alpha[3]*fUpwind_r[3]+alpha[1]*fUpwind_r[1]+alpha[0]*fUpwind_r[0]); 
  Ghat_r[1] += 0.1767766952966368*(alpha[3]*fUpwind_r[7]+fUpwind_r[3]*alpha[7]+alpha[0]*fUpwind_r[1]+fUpwind_r[0]*alpha[1]); 
  Ghat_r[2] += 0.1767766952966368*(alpha[7]*fUpwind_r[16]+alpha[3]*fUpwind_r[8]+alpha[1]*fUpwind_r[6]+alpha[0]*fUpwind_r[2]); 
  Ghat_r[3] += 0.1767766952966368*(alpha[1]*fUpwind_r[7]+fUpwind_r[1]*alpha[7]+alpha[0]*fUpwind_r[3]+fUpwind_r[0]*alpha[3]); 
  Ghat_r[4] += 0.1767766952966368*(alpha[7]*fUpwind_r[18]+alpha[3]*fUpwind_r[11]+alpha[1]*fUpwind_r[9]+alpha[0]*fUpwind_r[4]); 
  Ghat_r[5] += 0.1767766952966368*(alpha[7]*fUpwind_r[21]+alpha[3]*fUpwind_r[14]+alpha[1]*fUpwind_r[12]+alpha[0]*fUpwind_r[5]); 
  Ghat_r[6] += 0.1767766952966368*(alpha[3]*fUpwind_r[16]+alpha[7]*fUpwind_r[8]+alpha[0]*fUpwind_r[6]+alpha[1]*fUpwind_r[2]); 
  Ghat_r[7] += 0.1767766952966368*(alpha[0]*fUpwind_r[7]+fUpwind_r[0]*alpha[7]+alpha[1]*fUpwind_r[3]+fUpwind_r[1]*alpha[3]); 
  Ghat_r[8] += 0.1767766952966368*(alpha[1]*fUpwind_r[16]+alpha[0]*fUpwind_r[8]+fUpwind_r[6]*alpha[7]+fUpwind_r[2]*alpha[3]); 
  Ghat_r[9] += 0.1767766952966368*(alpha[3]*fUpwind_r[18]+alpha[7]*fUpwind_r[11]+alpha[0]*fUpwind_r[9]+alpha[1]*fUpwind_r[4]); 
  Ghat_r[10] += 0.1767766952966368*(alpha[7]*fUpwind_r[26]+alpha[3]*fUpwind_r[19]+alpha[1]*fUpwind_r[17]+alpha[0]*fUpwind_r[10]); 
  Ghat_r[11] += 0.1767766952966368*(alpha[1]*fUpwind_r[18]+alpha[0]*fUpwind_r[11]+alpha[7]*fUpwind_r[9]+alpha[3]*fUpwind_r[4]); 
  Ghat_r[12] += 0.1767766952966368*(alpha[3]*fUpwind_r[21]+alpha[7]*fUpwind_r[14]+alpha[0]*fUpwind_r[12]+alpha[1]*fUpwind_r[5]); 
  Ghat_r[13] += 0.1767766952966368*(alpha[7]*fUpwind_r[27]+alpha[3]*fUpwind_r[22]+alpha[1]*fUpwind_r[20]+alpha[0]*fUpwind_r[13]); 
  Ghat_r[14] += 0.1767766952966368*(alpha[1]*fUpwind_r[21]+alpha[0]*fUpwind_r[14]+alpha[7]*fUpwind_r[12]+alpha[3]*fUpwind_r[5]); 
  Ghat_r[15] += 0.1767766952966368*(alpha[7]*fUpwind_r[29]+alpha[3]*fUpwind_r[25]+alpha[1]*fUpwind_r[23]+alpha[0]*fUpwind_r[15]); 
  Ghat_r[16] += 0.1767766952966368*(alpha[0]*fUpwind_r[16]+alpha[1]*fUpwind_r[8]+fUpwind_r[2]*alpha[7]+alpha[3]*fUpwind_r[6]); 
  Ghat_r[17] += 0.1767766952966368*(alpha[3]*fUpwind_r[26]+alpha[7]*fUpwind_r[19]+alpha[0]*fUpwind_r[17]+alpha[1]*fUpwind_r[10]); 
  Ghat_r[18] += 0.1767766952966368*(alpha[0]*fUpwind_r[18]+alpha[1]*fUpwind_r[11]+alpha[3]*fUpwind_r[9]+fUpwind_r[4]*alpha[7]); 
  Ghat_r[19] += 0.1767766952966368*(alpha[1]*fUpwind_r[26]+alpha[0]*fUpwind_r[19]+alpha[7]*fUpwind_r[17]+alpha[3]*fUpwind_r[10]); 
  Ghat_r[20] += 0.1767766952966368*(alpha[3]*fUpwind_r[27]+alpha[7]*fUpwind_r[22]+alpha[0]*fUpwind_r[20]+alpha[1]*fUpwind_r[13]); 
  Ghat_r[21] += 0.1767766952966368*(alpha[0]*fUpwind_r[21]+alpha[1]*fUpwind_r[14]+alpha[3]*fUpwind_r[12]+fUpwind_r[5]*alpha[7]); 
  Ghat_r[22] += 0.1767766952966368*(alpha[1]*fUpwind_r[27]+alpha[0]*fUpwind_r[22]+alpha[7]*fUpwind_r[20]+alpha[3]*fUpwind_r[13]); 
  Ghat_r[23] += 0.1767766952966368*(alpha[3]*fUpwind_r[29]+alpha[7]*fUpwind_r[25]+alpha[0]*fUpwind_r[23]+alpha[1]*fUpwind_r[15]); 
  Ghat_r[24] += 0.1767766952966368*(alpha[7]*fUpwind_r[31]+alpha[3]*fUpwind_r[30]+alpha[1]*fUpwind_r[28]+alpha[0]*fUpwind_r[24]); 
  Ghat_r[25] += 0.1767766952966368*(alpha[1]*fUpwind_r[29]+alpha[0]*fUpwind_r[25]+alpha[7]*fUpwind_r[23]+alpha[3]*fUpwind_r[15]); 
  Ghat_r[26] += 0.1767766952966368*(alpha[0]*fUpwind_r[26]+alpha[1]*fUpwind_r[19]+alpha[3]*fUpwind_r[17]+alpha[7]*fUpwind_r[10]); 
  Ghat_r[27] += 0.1767766952966368*(alpha[0]*fUpwind_r[27]+alpha[1]*fUpwind_r[22]+alpha[3]*fUpwind_r[20]+alpha[7]*fUpwind_r[13]); 
  Ghat_r[28] += 0.1767766952966368*(alpha[3]*fUpwind_r[31]+alpha[7]*fUpwind_r[30]+alpha[0]*fUpwind_r[28]+alpha[1]*fUpwind_r[24]); 
  Ghat_r[29] += 0.1767766952966368*(alpha[0]*fUpwind_r[29]+alpha[1]*fUpwind_r[25]+alpha[3]*fUpwind_r[23]+alpha[7]*fUpwind_r[15]); 
  Ghat_r[30] += 0.1767766952966368*(alpha[1]*fUpwind_r[31]+alpha[0]*fUpwind_r[30]+alpha[7]*fUpwind_r[28]+alpha[3]*fUpwind_r[24]); 
  Ghat_r[31] += 0.1767766952966368*(alpha[0]*fUpwind_r[31]+alpha[1]*fUpwind_r[30]+alpha[3]*fUpwind_r[28]+alpha[7]*fUpwind_r[24]); 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv11; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv11; 
  out[2] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dv11; 
  out[3] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dv11; 
  out[4] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*dv11; 
  out[5] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv11; 
  out[6] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*dv11; 
  out[7] += (0.7071067811865475*Ghat_l[6]-0.7071067811865475*Ghat_r[6])*dv11; 
  out[8] += (0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*dv11; 
  out[9] += (0.7071067811865475*Ghat_l[8]-0.7071067811865475*Ghat_r[8])*dv11; 
  out[10] += (0.7071067811865475*Ghat_l[9]-0.7071067811865475*Ghat_r[9])*dv11; 
  out[11] += (0.7071067811865475*Ghat_l[10]-0.7071067811865475*Ghat_r[10])*dv11; 
  out[12] += (0.7071067811865475*Ghat_l[11]-0.7071067811865475*Ghat_r[11])*dv11; 
  out[13] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv11; 
  out[14] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dv11; 
  out[15] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dv11; 
  out[16] += -1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dv11; 
  out[17] += (0.7071067811865475*Ghat_l[12]-0.7071067811865475*Ghat_r[12])*dv11; 
  out[18] += (0.7071067811865475*Ghat_l[13]-0.7071067811865475*Ghat_r[13])*dv11; 
  out[19] += (0.7071067811865475*Ghat_l[14]-0.7071067811865475*Ghat_r[14])*dv11; 
  out[20] += (0.7071067811865475*Ghat_l[15]-0.7071067811865475*Ghat_r[15])*dv11; 
  out[21] += -1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dv11; 
  out[22] += (0.7071067811865475*Ghat_l[16]-0.7071067811865475*Ghat_r[16])*dv11; 
  out[23] += (0.7071067811865475*Ghat_l[17]-0.7071067811865475*Ghat_r[17])*dv11; 
  out[24] += (0.7071067811865475*Ghat_l[18]-0.7071067811865475*Ghat_r[18])*dv11; 
  out[25] += (0.7071067811865475*Ghat_l[19]-0.7071067811865475*Ghat_r[19])*dv11; 
  out[26] += -1.224744871391589*(Ghat_r[6]+Ghat_l[6])*dv11; 
  out[27] += -1.224744871391589*(Ghat_r[7]+Ghat_l[7])*dv11; 
  out[28] += -1.224744871391589*(Ghat_r[8]+Ghat_l[8])*dv11; 
  out[29] += -1.224744871391589*(Ghat_r[9]+Ghat_l[9])*dv11; 
  out[30] += -1.224744871391589*(Ghat_r[10]+Ghat_l[10])*dv11; 
  out[31] += -1.224744871391589*(Ghat_r[11]+Ghat_l[11])*dv11; 
  out[32] += (0.7071067811865475*Ghat_l[20]-0.7071067811865475*Ghat_r[20])*dv11; 
  out[33] += (0.7071067811865475*Ghat_l[21]-0.7071067811865475*Ghat_r[21])*dv11; 
  out[34] += (0.7071067811865475*Ghat_l[22]-0.7071067811865475*Ghat_r[22])*dv11; 
  out[35] += (0.7071067811865475*Ghat_l[23]-0.7071067811865475*Ghat_r[23])*dv11; 
  out[36] += (0.7071067811865475*Ghat_l[24]-0.7071067811865475*Ghat_r[24])*dv11; 
  out[37] += (0.7071067811865475*Ghat_l[25]-0.7071067811865475*Ghat_r[25])*dv11; 
  out[38] += -1.224744871391589*(Ghat_r[12]+Ghat_l[12])*dv11; 
  out[39] += -1.224744871391589*(Ghat_r[13]+Ghat_l[13])*dv11; 
  out[40] += -1.224744871391589*(Ghat_r[14]+Ghat_l[14])*dv11; 
  out[41] += -1.224744871391589*(Ghat_r[15]+Ghat_l[15])*dv11; 
  out[42] += (0.7071067811865475*Ghat_l[26]-0.7071067811865475*Ghat_r[26])*dv11; 
  out[43] += -1.224744871391589*(Ghat_r[16]+Ghat_l[16])*dv11; 
  out[44] += -1.224744871391589*(Ghat_r[17]+Ghat_l[17])*dv11; 
  out[45] += -1.224744871391589*(Ghat_r[18]+Ghat_l[18])*dv11; 
  out[46] += -1.224744871391589*(Ghat_r[19]+Ghat_l[19])*dv11; 
  out[47] += (0.7071067811865475*Ghat_l[27]-0.7071067811865475*Ghat_r[27])*dv11; 
  out[48] += (0.7071067811865475*Ghat_l[28]-0.7071067811865475*Ghat_r[28])*dv11; 
  out[49] += (0.7071067811865475*Ghat_l[29]-0.7071067811865475*Ghat_r[29])*dv11; 
  out[50] += (0.7071067811865475*Ghat_l[30]-0.7071067811865475*Ghat_r[30])*dv11; 
  out[51] += -1.224744871391589*(Ghat_r[20]+Ghat_l[20])*dv11; 
  out[52] += -1.224744871391589*(Ghat_r[21]+Ghat_l[21])*dv11; 
  out[53] += -1.224744871391589*(Ghat_r[22]+Ghat_l[22])*dv11; 
  out[54] += -1.224744871391589*(Ghat_r[23]+Ghat_l[23])*dv11; 
  out[55] += -1.224744871391589*(Ghat_r[24]+Ghat_l[24])*dv11; 
  out[56] += -1.224744871391589*(Ghat_r[25]+Ghat_l[25])*dv11; 
  out[57] += -1.224744871391589*(Ghat_r[26]+Ghat_l[26])*dv11; 
  out[58] += (0.7071067811865475*Ghat_l[31]-0.7071067811865475*Ghat_r[31])*dv11; 
  out[59] += -1.224744871391589*(Ghat_r[27]+Ghat_l[27])*dv11; 
  out[60] += -1.224744871391589*(Ghat_r[28]+Ghat_l[28])*dv11; 
  out[61] += -1.224744871391589*(Ghat_r[29]+Ghat_l[29])*dv11; 
  out[62] += -1.224744871391589*(Ghat_r[30]+Ghat_l[30])*dv11; 
  out[63] += -1.224744871391589*(Ghat_r[31]+Ghat_l[31])*dv11; 

} 
