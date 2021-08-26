#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_ser_3x3v_p1_surfvx_quad.h> 
GKYL_CU_DH void vlasov_surfvx_3x3v_ser_p1(const double *w, const double *dxv, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // qmem:      q/m*EM fields.
  // fl/fc/fr:  Input Distribution function in left/center/right cells 
  // out:       Output distribution function in center cell 
  const double dv10 = 2/dxv[3]; 
  const double dv1 = dxv[3], wv1 = w[3]; 
  const double dv2 = dxv[4], wv2 = w[4]; 
  const double dv3 = dxv[5], wv3 = w[5]; 
  const double *E0 = &qmem[0]; 
  const double *B0 = &qmem[24]; 
  const double *B1 = &qmem[32]; 
  const double *B2 = &qmem[40]; 

  double alpha[32] = {0.0}; 

  alpha[0] = (-2.0*B1[0]*wv3)+2.0*B2[0]*wv2+2.0*E0[0]; 
  alpha[1] = (-2.0*B1[1]*wv3)+2.0*B2[1]*wv2+2.0*E0[1]; 
  alpha[2] = (-2.0*B1[2]*wv3)+2.0*B2[2]*wv2+2.0*E0[2]; 
  alpha[3] = (-2.0*B1[3]*wv3)+2.0*B2[3]*wv2+2.0*E0[3]; 
  alpha[4] = 0.5773502691896258*B2[0]*dv2; 
  alpha[5] = -0.5773502691896258*B1[0]*dv3; 
  alpha[6] = (-2.0*B1[4]*wv3)+2.0*B2[4]*wv2+2.0*E0[4]; 
  alpha[7] = (-2.0*B1[5]*wv3)+2.0*B2[5]*wv2+2.0*E0[5]; 
  alpha[8] = (-2.0*B1[6]*wv3)+2.0*B2[6]*wv2+2.0*E0[6]; 
  alpha[9] = 0.5773502691896258*B2[1]*dv2; 
  alpha[10] = 0.5773502691896258*B2[2]*dv2; 
  alpha[11] = 0.5773502691896258*B2[3]*dv2; 
  alpha[12] = -0.5773502691896258*B1[1]*dv3; 
  alpha[13] = -0.5773502691896258*B1[2]*dv3; 
  alpha[14] = -0.5773502691896258*B1[3]*dv3; 
  alpha[16] = (-2.0*B1[7]*wv3)+2.0*B2[7]*wv2+2.0*E0[7]; 
  alpha[17] = 0.5773502691896258*B2[4]*dv2; 
  alpha[18] = 0.5773502691896258*B2[5]*dv2; 
  alpha[19] = 0.5773502691896258*B2[6]*dv2; 
  alpha[20] = -0.5773502691896258*B1[4]*dv3; 
  alpha[21] = -0.5773502691896258*B1[5]*dv3; 
  alpha[22] = -0.5773502691896258*B1[6]*dv3; 
  alpha[26] = 0.5773502691896258*B2[7]*dv2; 
  alpha[27] = -0.5773502691896258*B1[7]*dv3; 

  double fUpwindQuad_l[32] = {0.0};
  double fUpwindQuad_r[32] = {0.0};
  double fUpwind_l[32] = {0.0};;
  double fUpwind_r[32] = {0.0};
  double Ghat_l[32] = {0.0}; 
  double Ghat_r[32] = {0.0}; 

  if (alpha[27]+alpha[26]-alpha[22]-alpha[21]-alpha[20]-alpha[19]-alpha[18]-alpha[17]-alpha[16]+alpha[14]+alpha[13]+alpha[12]+alpha[11]+alpha[10]+alpha[9]+alpha[8]+alpha[7]+alpha[6]-alpha[5]-alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[0] = ser_3x3v_p1_surfvx_quad_0(1, fl); 
    fUpwindQuad_r[0] = ser_3x3v_p1_surfvx_quad_0(1, fc); 
  } else { 

    fUpwindQuad_l[0] = ser_3x3v_p1_surfvx_quad_0(-1, fc); 
    fUpwindQuad_r[0] = ser_3x3v_p1_surfvx_quad_0(-1, fr); 
  } 
  if ((-alpha[27])-alpha[26]-alpha[22]+alpha[21]+alpha[20]-alpha[19]+alpha[18]+alpha[17]+alpha[16]+alpha[14]+alpha[13]-alpha[12]+alpha[11]+alpha[10]-alpha[9]+alpha[8]-alpha[7]-alpha[6]-alpha[5]-alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[1] = ser_3x3v_p1_surfvx_quad_1(1, fl); 
    fUpwindQuad_r[1] = ser_3x3v_p1_surfvx_quad_1(1, fc); 
  } else { 

    fUpwindQuad_l[1] = ser_3x3v_p1_surfvx_quad_1(-1, fc); 
    fUpwindQuad_r[1] = ser_3x3v_p1_surfvx_quad_1(-1, fr); 
  } 
  if ((-alpha[27])-alpha[26]+alpha[22]-alpha[21]+alpha[20]+alpha[19]-alpha[18]+alpha[17]+alpha[16]+alpha[14]-alpha[13]+alpha[12]+alpha[11]-alpha[10]+alpha[9]-alpha[8]+alpha[7]-alpha[6]-alpha[5]-alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[2] = ser_3x3v_p1_surfvx_quad_2(1, fl); 
    fUpwindQuad_r[2] = ser_3x3v_p1_surfvx_quad_2(1, fc); 
  } else { 

    fUpwindQuad_l[2] = ser_3x3v_p1_surfvx_quad_2(-1, fc); 
    fUpwindQuad_r[2] = ser_3x3v_p1_surfvx_quad_2(-1, fr); 
  } 
  if (alpha[27]+alpha[26]+alpha[22]+alpha[21]-alpha[20]+alpha[19]+alpha[18]-alpha[17]-alpha[16]+alpha[14]-alpha[13]-alpha[12]+alpha[11]-alpha[10]-alpha[9]-alpha[8]-alpha[7]+alpha[6]-alpha[5]-alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[3] = ser_3x3v_p1_surfvx_quad_3(1, fl); 
    fUpwindQuad_r[3] = ser_3x3v_p1_surfvx_quad_3(1, fc); 
  } else { 

    fUpwindQuad_l[3] = ser_3x3v_p1_surfvx_quad_3(-1, fc); 
    fUpwindQuad_r[3] = ser_3x3v_p1_surfvx_quad_3(-1, fr); 
  } 
  if ((-alpha[27])-alpha[26]+alpha[22]+alpha[21]-alpha[20]+alpha[19]+alpha[18]-alpha[17]+alpha[16]-alpha[14]+alpha[13]+alpha[12]-alpha[11]+alpha[10]+alpha[9]-alpha[8]-alpha[7]+alpha[6]-alpha[5]-alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[4] = ser_3x3v_p1_surfvx_quad_4(1, fl); 
    fUpwindQuad_r[4] = ser_3x3v_p1_surfvx_quad_4(1, fc); 
  } else { 

    fUpwindQuad_l[4] = ser_3x3v_p1_surfvx_quad_4(-1, fc); 
    fUpwindQuad_r[4] = ser_3x3v_p1_surfvx_quad_4(-1, fr); 
  } 
  if (alpha[27]+alpha[26]+alpha[22]-alpha[21]+alpha[20]+alpha[19]-alpha[18]+alpha[17]-alpha[16]-alpha[14]+alpha[13]-alpha[12]-alpha[11]+alpha[10]-alpha[9]-alpha[8]+alpha[7]-alpha[6]-alpha[5]-alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[5] = ser_3x3v_p1_surfvx_quad_5(1, fl); 
    fUpwindQuad_r[5] = ser_3x3v_p1_surfvx_quad_5(1, fc); 
  } else { 

    fUpwindQuad_l[5] = ser_3x3v_p1_surfvx_quad_5(-1, fc); 
    fUpwindQuad_r[5] = ser_3x3v_p1_surfvx_quad_5(-1, fr); 
  } 
  if (alpha[27]+alpha[26]-alpha[22]+alpha[21]+alpha[20]-alpha[19]+alpha[18]+alpha[17]-alpha[16]-alpha[14]-alpha[13]+alpha[12]-alpha[11]-alpha[10]+alpha[9]+alpha[8]-alpha[7]-alpha[6]-alpha[5]-alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[6] = ser_3x3v_p1_surfvx_quad_6(1, fl); 
    fUpwindQuad_r[6] = ser_3x3v_p1_surfvx_quad_6(1, fc); 
  } else { 

    fUpwindQuad_l[6] = ser_3x3v_p1_surfvx_quad_6(-1, fc); 
    fUpwindQuad_r[6] = ser_3x3v_p1_surfvx_quad_6(-1, fr); 
  } 
  if ((-alpha[27])-alpha[26]-alpha[22]-alpha[21]-alpha[20]-alpha[19]-alpha[18]-alpha[17]+alpha[16]-alpha[14]-alpha[13]-alpha[12]-alpha[11]-alpha[10]-alpha[9]+alpha[8]+alpha[7]+alpha[6]-alpha[5]-alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[7] = ser_3x3v_p1_surfvx_quad_7(1, fl); 
    fUpwindQuad_r[7] = ser_3x3v_p1_surfvx_quad_7(1, fc); 
  } else { 

    fUpwindQuad_l[7] = ser_3x3v_p1_surfvx_quad_7(-1, fc); 
    fUpwindQuad_r[7] = ser_3x3v_p1_surfvx_quad_7(-1, fr); 
  } 
  if (alpha[27]-alpha[26]-alpha[22]-alpha[21]-alpha[20]+alpha[19]+alpha[18]+alpha[17]-alpha[16]+alpha[14]+alpha[13]+alpha[12]-alpha[11]-alpha[10]-alpha[9]+alpha[8]+alpha[7]+alpha[6]-alpha[5]+alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[8] = ser_3x3v_p1_surfvx_quad_8(1, fl); 
    fUpwindQuad_r[8] = ser_3x3v_p1_surfvx_quad_8(1, fc); 
  } else { 

    fUpwindQuad_l[8] = ser_3x3v_p1_surfvx_quad_8(-1, fc); 
    fUpwindQuad_r[8] = ser_3x3v_p1_surfvx_quad_8(-1, fr); 
  } 
  if ((-alpha[27])+alpha[26]-alpha[22]+alpha[21]+alpha[20]+alpha[19]-alpha[18]-alpha[17]+alpha[16]+alpha[14]+alpha[13]-alpha[12]-alpha[11]-alpha[10]+alpha[9]+alpha[8]-alpha[7]-alpha[6]-alpha[5]+alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[9] = ser_3x3v_p1_surfvx_quad_9(1, fl); 
    fUpwindQuad_r[9] = ser_3x3v_p1_surfvx_quad_9(1, fc); 
  } else { 

    fUpwindQuad_l[9] = ser_3x3v_p1_surfvx_quad_9(-1, fc); 
    fUpwindQuad_r[9] = ser_3x3v_p1_surfvx_quad_9(-1, fr); 
  } 
  if ((-alpha[27])+alpha[26]+alpha[22]-alpha[21]+alpha[20]-alpha[19]+alpha[18]-alpha[17]+alpha[16]+alpha[14]-alpha[13]+alpha[12]-alpha[11]+alpha[10]-alpha[9]-alpha[8]+alpha[7]-alpha[6]-alpha[5]+alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[10] = ser_3x3v_p1_surfvx_quad_10(1, fl); 
    fUpwindQuad_r[10] = ser_3x3v_p1_surfvx_quad_10(1, fc); 
  } else { 

    fUpwindQuad_l[10] = ser_3x3v_p1_surfvx_quad_10(-1, fc); 
    fUpwindQuad_r[10] = ser_3x3v_p1_surfvx_quad_10(-1, fr); 
  } 
  if (alpha[27]-alpha[26]+alpha[22]+alpha[21]-alpha[20]-alpha[19]-alpha[18]+alpha[17]-alpha[16]+alpha[14]-alpha[13]-alpha[12]-alpha[11]+alpha[10]+alpha[9]-alpha[8]-alpha[7]+alpha[6]-alpha[5]+alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[11] = ser_3x3v_p1_surfvx_quad_11(1, fl); 
    fUpwindQuad_r[11] = ser_3x3v_p1_surfvx_quad_11(1, fc); 
  } else { 

    fUpwindQuad_l[11] = ser_3x3v_p1_surfvx_quad_11(-1, fc); 
    fUpwindQuad_r[11] = ser_3x3v_p1_surfvx_quad_11(-1, fr); 
  } 
  if ((-alpha[27])+alpha[26]+alpha[22]+alpha[21]-alpha[20]-alpha[19]-alpha[18]+alpha[17]+alpha[16]-alpha[14]+alpha[13]+alpha[12]+alpha[11]-alpha[10]-alpha[9]-alpha[8]-alpha[7]+alpha[6]-alpha[5]+alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[12] = ser_3x3v_p1_surfvx_quad_12(1, fl); 
    fUpwindQuad_r[12] = ser_3x3v_p1_surfvx_quad_12(1, fc); 
  } else { 

    fUpwindQuad_l[12] = ser_3x3v_p1_surfvx_quad_12(-1, fc); 
    fUpwindQuad_r[12] = ser_3x3v_p1_surfvx_quad_12(-1, fr); 
  } 
  if (alpha[27]-alpha[26]+alpha[22]-alpha[21]+alpha[20]-alpha[19]+alpha[18]-alpha[17]-alpha[16]-alpha[14]+alpha[13]-alpha[12]+alpha[11]-alpha[10]+alpha[9]-alpha[8]+alpha[7]-alpha[6]-alpha[5]+alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[13] = ser_3x3v_p1_surfvx_quad_13(1, fl); 
    fUpwindQuad_r[13] = ser_3x3v_p1_surfvx_quad_13(1, fc); 
  } else { 

    fUpwindQuad_l[13] = ser_3x3v_p1_surfvx_quad_13(-1, fc); 
    fUpwindQuad_r[13] = ser_3x3v_p1_surfvx_quad_13(-1, fr); 
  } 
  if (alpha[27]-alpha[26]-alpha[22]+alpha[21]+alpha[20]+alpha[19]-alpha[18]-alpha[17]-alpha[16]-alpha[14]-alpha[13]+alpha[12]+alpha[11]+alpha[10]-alpha[9]+alpha[8]-alpha[7]-alpha[6]-alpha[5]+alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[14] = ser_3x3v_p1_surfvx_quad_14(1, fl); 
    fUpwindQuad_r[14] = ser_3x3v_p1_surfvx_quad_14(1, fc); 
  } else { 

    fUpwindQuad_l[14] = ser_3x3v_p1_surfvx_quad_14(-1, fc); 
    fUpwindQuad_r[14] = ser_3x3v_p1_surfvx_quad_14(-1, fr); 
  } 
  if ((-alpha[27])+alpha[26]-alpha[22]-alpha[21]-alpha[20]+alpha[19]+alpha[18]+alpha[17]+alpha[16]-alpha[14]-alpha[13]-alpha[12]+alpha[11]+alpha[10]+alpha[9]+alpha[8]+alpha[7]+alpha[6]-alpha[5]+alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[15] = ser_3x3v_p1_surfvx_quad_15(1, fl); 
    fUpwindQuad_r[15] = ser_3x3v_p1_surfvx_quad_15(1, fc); 
  } else { 

    fUpwindQuad_l[15] = ser_3x3v_p1_surfvx_quad_15(-1, fc); 
    fUpwindQuad_r[15] = ser_3x3v_p1_surfvx_quad_15(-1, fr); 
  } 
  if ((-alpha[27])+alpha[26]+alpha[22]+alpha[21]+alpha[20]-alpha[19]-alpha[18]-alpha[17]-alpha[16]-alpha[14]-alpha[13]-alpha[12]+alpha[11]+alpha[10]+alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5]-alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[16] = ser_3x3v_p1_surfvx_quad_16(1, fl); 
    fUpwindQuad_r[16] = ser_3x3v_p1_surfvx_quad_16(1, fc); 
  } else { 

    fUpwindQuad_l[16] = ser_3x3v_p1_surfvx_quad_16(-1, fc); 
    fUpwindQuad_r[16] = ser_3x3v_p1_surfvx_quad_16(-1, fr); 
  } 
  if (alpha[27]-alpha[26]+alpha[22]-alpha[21]-alpha[20]-alpha[19]+alpha[18]+alpha[17]+alpha[16]-alpha[14]-alpha[13]+alpha[12]+alpha[11]+alpha[10]-alpha[9]+alpha[8]-alpha[7]-alpha[6]+alpha[5]-alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[17] = ser_3x3v_p1_surfvx_quad_17(1, fl); 
    fUpwindQuad_r[17] = ser_3x3v_p1_surfvx_quad_17(1, fc); 
  } else { 

    fUpwindQuad_l[17] = ser_3x3v_p1_surfvx_quad_17(-1, fc); 
    fUpwindQuad_r[17] = ser_3x3v_p1_surfvx_quad_17(-1, fr); 
  } 
  if (alpha[27]-alpha[26]-alpha[22]+alpha[21]-alpha[20]+alpha[19]-alpha[18]+alpha[17]+alpha[16]-alpha[14]+alpha[13]-alpha[12]+alpha[11]-alpha[10]+alpha[9]-alpha[8]+alpha[7]-alpha[6]+alpha[5]-alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[18] = ser_3x3v_p1_surfvx_quad_18(1, fl); 
    fUpwindQuad_r[18] = ser_3x3v_p1_surfvx_quad_18(1, fc); 
  } else { 

    fUpwindQuad_l[18] = ser_3x3v_p1_surfvx_quad_18(-1, fc); 
    fUpwindQuad_r[18] = ser_3x3v_p1_surfvx_quad_18(-1, fr); 
  } 
  if ((-alpha[27])+alpha[26]-alpha[22]-alpha[21]+alpha[20]+alpha[19]+alpha[18]-alpha[17]-alpha[16]-alpha[14]+alpha[13]+alpha[12]+alpha[11]-alpha[10]-alpha[9]-alpha[8]-alpha[7]+alpha[6]+alpha[5]-alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[19] = ser_3x3v_p1_surfvx_quad_19(1, fl); 
    fUpwindQuad_r[19] = ser_3x3v_p1_surfvx_quad_19(1, fc); 
  } else { 

    fUpwindQuad_l[19] = ser_3x3v_p1_surfvx_quad_19(-1, fc); 
    fUpwindQuad_r[19] = ser_3x3v_p1_surfvx_quad_19(-1, fr); 
  } 
  if (alpha[27]-alpha[26]-alpha[22]-alpha[21]+alpha[20]+alpha[19]+alpha[18]-alpha[17]+alpha[16]+alpha[14]-alpha[13]-alpha[12]-alpha[11]+alpha[10]+alpha[9]-alpha[8]-alpha[7]+alpha[6]+alpha[5]-alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[20] = ser_3x3v_p1_surfvx_quad_20(1, fl); 
    fUpwindQuad_r[20] = ser_3x3v_p1_surfvx_quad_20(1, fc); 
  } else { 

    fUpwindQuad_l[20] = ser_3x3v_p1_surfvx_quad_20(-1, fc); 
    fUpwindQuad_r[20] = ser_3x3v_p1_surfvx_quad_20(-1, fr); 
  } 
  if ((-alpha[27])+alpha[26]-alpha[22]+alpha[21]-alpha[20]+alpha[19]-alpha[18]+alpha[17]-alpha[16]+alpha[14]-alpha[13]+alpha[12]-alpha[11]+alpha[10]-alpha[9]-alpha[8]+alpha[7]-alpha[6]+alpha[5]-alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[21] = ser_3x3v_p1_surfvx_quad_21(1, fl); 
    fUpwindQuad_r[21] = ser_3x3v_p1_surfvx_quad_21(1, fc); 
  } else { 

    fUpwindQuad_l[21] = ser_3x3v_p1_surfvx_quad_21(-1, fc); 
    fUpwindQuad_r[21] = ser_3x3v_p1_surfvx_quad_21(-1, fr); 
  } 
  if ((-alpha[27])+alpha[26]+alpha[22]-alpha[21]-alpha[20]-alpha[19]+alpha[18]+alpha[17]-alpha[16]+alpha[14]+alpha[13]-alpha[12]-alpha[11]-alpha[10]+alpha[9]+alpha[8]-alpha[7]-alpha[6]+alpha[5]-alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[22] = ser_3x3v_p1_surfvx_quad_22(1, fl); 
    fUpwindQuad_r[22] = ser_3x3v_p1_surfvx_quad_22(1, fc); 
  } else { 

    fUpwindQuad_l[22] = ser_3x3v_p1_surfvx_quad_22(-1, fc); 
    fUpwindQuad_r[22] = ser_3x3v_p1_surfvx_quad_22(-1, fr); 
  } 
  if (alpha[27]-alpha[26]+alpha[22]+alpha[21]+alpha[20]-alpha[19]-alpha[18]-alpha[17]+alpha[16]+alpha[14]+alpha[13]+alpha[12]-alpha[11]-alpha[10]-alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5]-alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[23] = ser_3x3v_p1_surfvx_quad_23(1, fl); 
    fUpwindQuad_r[23] = ser_3x3v_p1_surfvx_quad_23(1, fc); 
  } else { 

    fUpwindQuad_l[23] = ser_3x3v_p1_surfvx_quad_23(-1, fc); 
    fUpwindQuad_r[23] = ser_3x3v_p1_surfvx_quad_23(-1, fr); 
  } 
  if ((-alpha[27])-alpha[26]+alpha[22]+alpha[21]+alpha[20]+alpha[19]+alpha[18]+alpha[17]-alpha[16]-alpha[14]-alpha[13]-alpha[12]-alpha[11]-alpha[10]-alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5]+alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[24] = ser_3x3v_p1_surfvx_quad_24(1, fl); 
    fUpwindQuad_r[24] = ser_3x3v_p1_surfvx_quad_24(1, fc); 
  } else { 

    fUpwindQuad_l[24] = ser_3x3v_p1_surfvx_quad_24(-1, fc); 
    fUpwindQuad_r[24] = ser_3x3v_p1_surfvx_quad_24(-1, fr); 
  } 
  if (alpha[27]+alpha[26]+alpha[22]-alpha[21]-alpha[20]+alpha[19]-alpha[18]-alpha[17]+alpha[16]-alpha[14]-alpha[13]+alpha[12]-alpha[11]-alpha[10]+alpha[9]+alpha[8]-alpha[7]-alpha[6]+alpha[5]+alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[25] = ser_3x3v_p1_surfvx_quad_25(1, fl); 
    fUpwindQuad_r[25] = ser_3x3v_p1_surfvx_quad_25(1, fc); 
  } else { 

    fUpwindQuad_l[25] = ser_3x3v_p1_surfvx_quad_25(-1, fc); 
    fUpwindQuad_r[25] = ser_3x3v_p1_surfvx_quad_25(-1, fr); 
  } 
  if (alpha[27]+alpha[26]-alpha[22]+alpha[21]-alpha[20]-alpha[19]+alpha[18]-alpha[17]+alpha[16]-alpha[14]+alpha[13]-alpha[12]-alpha[11]+alpha[10]-alpha[9]-alpha[8]+alpha[7]-alpha[6]+alpha[5]+alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[26] = ser_3x3v_p1_surfvx_quad_26(1, fl); 
    fUpwindQuad_r[26] = ser_3x3v_p1_surfvx_quad_26(1, fc); 
  } else { 

    fUpwindQuad_l[26] = ser_3x3v_p1_surfvx_quad_26(-1, fc); 
    fUpwindQuad_r[26] = ser_3x3v_p1_surfvx_quad_26(-1, fr); 
  } 
  if ((-alpha[27])-alpha[26]-alpha[22]-alpha[21]+alpha[20]-alpha[19]-alpha[18]+alpha[17]-alpha[16]-alpha[14]+alpha[13]+alpha[12]-alpha[11]+alpha[10]+alpha[9]-alpha[8]-alpha[7]+alpha[6]+alpha[5]+alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[27] = ser_3x3v_p1_surfvx_quad_27(1, fl); 
    fUpwindQuad_r[27] = ser_3x3v_p1_surfvx_quad_27(1, fc); 
  } else { 

    fUpwindQuad_l[27] = ser_3x3v_p1_surfvx_quad_27(-1, fc); 
    fUpwindQuad_r[27] = ser_3x3v_p1_surfvx_quad_27(-1, fr); 
  } 
  if (alpha[27]+alpha[26]-alpha[22]-alpha[21]+alpha[20]-alpha[19]-alpha[18]+alpha[17]+alpha[16]+alpha[14]-alpha[13]-alpha[12]+alpha[11]-alpha[10]-alpha[9]-alpha[8]-alpha[7]+alpha[6]+alpha[5]+alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[28] = ser_3x3v_p1_surfvx_quad_28(1, fl); 
    fUpwindQuad_r[28] = ser_3x3v_p1_surfvx_quad_28(1, fc); 
  } else { 

    fUpwindQuad_l[28] = ser_3x3v_p1_surfvx_quad_28(-1, fc); 
    fUpwindQuad_r[28] = ser_3x3v_p1_surfvx_quad_28(-1, fr); 
  } 
  if ((-alpha[27])-alpha[26]-alpha[22]+alpha[21]-alpha[20]-alpha[19]+alpha[18]-alpha[17]-alpha[16]+alpha[14]-alpha[13]+alpha[12]+alpha[11]-alpha[10]+alpha[9]-alpha[8]+alpha[7]-alpha[6]+alpha[5]+alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[29] = ser_3x3v_p1_surfvx_quad_29(1, fl); 
    fUpwindQuad_r[29] = ser_3x3v_p1_surfvx_quad_29(1, fc); 
  } else { 

    fUpwindQuad_l[29] = ser_3x3v_p1_surfvx_quad_29(-1, fc); 
    fUpwindQuad_r[29] = ser_3x3v_p1_surfvx_quad_29(-1, fr); 
  } 
  if ((-alpha[27])-alpha[26]+alpha[22]-alpha[21]-alpha[20]+alpha[19]-alpha[18]-alpha[17]-alpha[16]+alpha[14]+alpha[13]-alpha[12]+alpha[11]+alpha[10]-alpha[9]+alpha[8]-alpha[7]-alpha[6]+alpha[5]+alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[30] = ser_3x3v_p1_surfvx_quad_30(1, fl); 
    fUpwindQuad_r[30] = ser_3x3v_p1_surfvx_quad_30(1, fc); 
  } else { 

    fUpwindQuad_l[30] = ser_3x3v_p1_surfvx_quad_30(-1, fc); 
    fUpwindQuad_r[30] = ser_3x3v_p1_surfvx_quad_30(-1, fr); 
  } 
  if (alpha[27]+alpha[26]+alpha[22]+alpha[21]+alpha[20]+alpha[19]+alpha[18]+alpha[17]+alpha[16]+alpha[14]+alpha[13]+alpha[12]+alpha[11]+alpha[10]+alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5]+alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[31] = ser_3x3v_p1_surfvx_quad_31(1, fl); 
    fUpwindQuad_r[31] = ser_3x3v_p1_surfvx_quad_31(1, fc); 
  } else { 

    fUpwindQuad_l[31] = ser_3x3v_p1_surfvx_quad_31(-1, fc); 
    fUpwindQuad_r[31] = ser_3x3v_p1_surfvx_quad_31(-1, fr); 
  } 
  fUpwind_l[0] = 0.1767766952966368*(fUpwindQuad_l[31]+fUpwindQuad_l[30]+fUpwindQuad_l[29]+fUpwindQuad_l[28]+fUpwindQuad_l[27]+fUpwindQuad_l[26]+fUpwindQuad_l[25]+fUpwindQuad_l[24]+fUpwindQuad_l[23]+fUpwindQuad_l[22]+fUpwindQuad_l[21]+fUpwindQuad_l[20]+fUpwindQuad_l[19]+fUpwindQuad_l[18]+fUpwindQuad_l[17]+fUpwindQuad_l[16]+fUpwindQuad_l[15]+fUpwindQuad_l[14]+fUpwindQuad_l[13]+fUpwindQuad_l[12]+fUpwindQuad_l[11]+fUpwindQuad_l[10]+fUpwindQuad_l[9]+fUpwindQuad_l[8]+fUpwindQuad_l[7]+fUpwindQuad_l[6]+fUpwindQuad_l[5]+fUpwindQuad_l[4]+fUpwindQuad_l[3]+fUpwindQuad_l[2]+fUpwindQuad_l[1]+fUpwindQuad_l[0]); 
  fUpwind_l[1] = 0.1767766952966368*(fUpwindQuad_l[31]-1.0*fUpwindQuad_l[30]+fUpwindQuad_l[29]-1.0*fUpwindQuad_l[28]+fUpwindQuad_l[27]-1.0*fUpwindQuad_l[26]+fUpwindQuad_l[25]-1.0*fUpwindQuad_l[24]+fUpwindQuad_l[23]-1.0*fUpwindQuad_l[22]+fUpwindQuad_l[21]-1.0*fUpwindQuad_l[20]+fUpwindQuad_l[19]-1.0*fUpwindQuad_l[18]+fUpwindQuad_l[17]-1.0*fUpwindQuad_l[16]+fUpwindQuad_l[15]-1.0*fUpwindQuad_l[14]+fUpwindQuad_l[13]-1.0*fUpwindQuad_l[12]+fUpwindQuad_l[11]-1.0*fUpwindQuad_l[10]+fUpwindQuad_l[9]-1.0*fUpwindQuad_l[8]+fUpwindQuad_l[7]-1.0*fUpwindQuad_l[6]+fUpwindQuad_l[5]-1.0*fUpwindQuad_l[4]+fUpwindQuad_l[3]-1.0*fUpwindQuad_l[2]+fUpwindQuad_l[1]-1.0*fUpwindQuad_l[0]); 
  fUpwind_l[2] = 0.1767766952966368*(fUpwindQuad_l[31]+fUpwindQuad_l[30]-1.0*(fUpwindQuad_l[29]+fUpwindQuad_l[28])+fUpwindQuad_l[27]+fUpwindQuad_l[26]-1.0*(fUpwindQuad_l[25]+fUpwindQuad_l[24])+fUpwindQuad_l[23]+fUpwindQuad_l[22]-1.0*(fUpwindQuad_l[21]+fUpwindQuad_l[20])+fUpwindQuad_l[19]+fUpwindQuad_l[18]-1.0*(fUpwindQuad_l[17]+fUpwindQuad_l[16])+fUpwindQuad_l[15]+fUpwindQuad_l[14]-1.0*(fUpwindQuad_l[13]+fUpwindQuad_l[12])+fUpwindQuad_l[11]+fUpwindQuad_l[10]-1.0*(fUpwindQuad_l[9]+fUpwindQuad_l[8])+fUpwindQuad_l[7]+fUpwindQuad_l[6]-1.0*(fUpwindQuad_l[5]+fUpwindQuad_l[4])+fUpwindQuad_l[3]+fUpwindQuad_l[2]-1.0*(fUpwindQuad_l[1]+fUpwindQuad_l[0])); 
  fUpwind_l[3] = 0.1767766952966368*(fUpwindQuad_l[31]+fUpwindQuad_l[30]+fUpwindQuad_l[29]+fUpwindQuad_l[28]-1.0*(fUpwindQuad_l[27]+fUpwindQuad_l[26]+fUpwindQuad_l[25]+fUpwindQuad_l[24])+fUpwindQuad_l[23]+fUpwindQuad_l[22]+fUpwindQuad_l[21]+fUpwindQuad_l[20]-1.0*(fUpwindQuad_l[19]+fUpwindQuad_l[18]+fUpwindQuad_l[17]+fUpwindQuad_l[16])+fUpwindQuad_l[15]+fUpwindQuad_l[14]+fUpwindQuad_l[13]+fUpwindQuad_l[12]-1.0*(fUpwindQuad_l[11]+fUpwindQuad_l[10]+fUpwindQuad_l[9]+fUpwindQuad_l[8])+fUpwindQuad_l[7]+fUpwindQuad_l[6]+fUpwindQuad_l[5]+fUpwindQuad_l[4]-1.0*(fUpwindQuad_l[3]+fUpwindQuad_l[2]+fUpwindQuad_l[1]+fUpwindQuad_l[0])); 
  fUpwind_l[4] = 0.1767766952966368*(fUpwindQuad_l[31]+fUpwindQuad_l[30]+fUpwindQuad_l[29]+fUpwindQuad_l[28]+fUpwindQuad_l[27]+fUpwindQuad_l[26]+fUpwindQuad_l[25]+fUpwindQuad_l[24]-1.0*(fUpwindQuad_l[23]+fUpwindQuad_l[22]+fUpwindQuad_l[21]+fUpwindQuad_l[20]+fUpwindQuad_l[19]+fUpwindQuad_l[18]+fUpwindQuad_l[17]+fUpwindQuad_l[16])+fUpwindQuad_l[15]+fUpwindQuad_l[14]+fUpwindQuad_l[13]+fUpwindQuad_l[12]+fUpwindQuad_l[11]+fUpwindQuad_l[10]+fUpwindQuad_l[9]+fUpwindQuad_l[8]-1.0*(fUpwindQuad_l[7]+fUpwindQuad_l[6]+fUpwindQuad_l[5]+fUpwindQuad_l[4]+fUpwindQuad_l[3]+fUpwindQuad_l[2]+fUpwindQuad_l[1]+fUpwindQuad_l[0])); 
  fUpwind_l[5] = 0.1767766952966368*(fUpwindQuad_l[31]+fUpwindQuad_l[30]+fUpwindQuad_l[29]+fUpwindQuad_l[28]+fUpwindQuad_l[27]+fUpwindQuad_l[26]+fUpwindQuad_l[25]+fUpwindQuad_l[24]+fUpwindQuad_l[23]+fUpwindQuad_l[22]+fUpwindQuad_l[21]+fUpwindQuad_l[20]+fUpwindQuad_l[19]+fUpwindQuad_l[18]+fUpwindQuad_l[17]+fUpwindQuad_l[16]-1.0*(fUpwindQuad_l[15]+fUpwindQuad_l[14]+fUpwindQuad_l[13]+fUpwindQuad_l[12]+fUpwindQuad_l[11]+fUpwindQuad_l[10]+fUpwindQuad_l[9]+fUpwindQuad_l[8]+fUpwindQuad_l[7]+fUpwindQuad_l[6]+fUpwindQuad_l[5]+fUpwindQuad_l[4]+fUpwindQuad_l[3]+fUpwindQuad_l[2]+fUpwindQuad_l[1]+fUpwindQuad_l[0])); 
  fUpwind_l[6] = 0.1767766952966368*(fUpwindQuad_l[31]-1.0*(fUpwindQuad_l[30]+fUpwindQuad_l[29])+fUpwindQuad_l[28]+fUpwindQuad_l[27]-1.0*(fUpwindQuad_l[26]+fUpwindQuad_l[25])+fUpwindQuad_l[24]+fUpwindQuad_l[23]-1.0*(fUpwindQuad_l[22]+fUpwindQuad_l[21])+fUpwindQuad_l[20]+fUpwindQuad_l[19]-1.0*(fUpwindQuad_l[18]+fUpwindQuad_l[17])+fUpwindQuad_l[16]+fUpwindQuad_l[15]-1.0*(fUpwindQuad_l[14]+fUpwindQuad_l[13])+fUpwindQuad_l[12]+fUpwindQuad_l[11]-1.0*(fUpwindQuad_l[10]+fUpwindQuad_l[9])+fUpwindQuad_l[8]+fUpwindQuad_l[7]-1.0*(fUpwindQuad_l[6]+fUpwindQuad_l[5])+fUpwindQuad_l[4]+fUpwindQuad_l[3]-1.0*(fUpwindQuad_l[2]+fUpwindQuad_l[1])+fUpwindQuad_l[0]); 
  fUpwind_l[7] = 0.1767766952966368*(fUpwindQuad_l[31]-1.0*fUpwindQuad_l[30]+fUpwindQuad_l[29]-1.0*(fUpwindQuad_l[28]+fUpwindQuad_l[27])+fUpwindQuad_l[26]-1.0*fUpwindQuad_l[25]+fUpwindQuad_l[24]+fUpwindQuad_l[23]-1.0*fUpwindQuad_l[22]+fUpwindQuad_l[21]-1.0*(fUpwindQuad_l[20]+fUpwindQuad_l[19])+fUpwindQuad_l[18]-1.0*fUpwindQuad_l[17]+fUpwindQuad_l[16]+fUpwindQuad_l[15]-1.0*fUpwindQuad_l[14]+fUpwindQuad_l[13]-1.0*(fUpwindQuad_l[12]+fUpwindQuad_l[11])+fUpwindQuad_l[10]-1.0*fUpwindQuad_l[9]+fUpwindQuad_l[8]+fUpwindQuad_l[7]-1.0*fUpwindQuad_l[6]+fUpwindQuad_l[5]-1.0*(fUpwindQuad_l[4]+fUpwindQuad_l[3])+fUpwindQuad_l[2]-1.0*fUpwindQuad_l[1]+fUpwindQuad_l[0]); 
  fUpwind_l[8] = 0.1767766952966368*(fUpwindQuad_l[31]+fUpwindQuad_l[30]-1.0*(fUpwindQuad_l[29]+fUpwindQuad_l[28]+fUpwindQuad_l[27]+fUpwindQuad_l[26])+fUpwindQuad_l[25]+fUpwindQuad_l[24]+fUpwindQuad_l[23]+fUpwindQuad_l[22]-1.0*(fUpwindQuad_l[21]+fUpwindQuad_l[20]+fUpwindQuad_l[19]+fUpwindQuad_l[18])+fUpwindQuad_l[17]+fUpwindQuad_l[16]+fUpwindQuad_l[15]+fUpwindQuad_l[14]-1.0*(fUpwindQuad_l[13]+fUpwindQuad_l[12]+fUpwindQuad_l[11]+fUpwindQuad_l[10])+fUpwindQuad_l[9]+fUpwindQuad_l[8]+fUpwindQuad_l[7]+fUpwindQuad_l[6]-1.0*(fUpwindQuad_l[5]+fUpwindQuad_l[4]+fUpwindQuad_l[3]+fUpwindQuad_l[2])+fUpwindQuad_l[1]+fUpwindQuad_l[0]); 
  fUpwind_l[9] = 0.1767766952966368*(fUpwindQuad_l[31]-1.0*fUpwindQuad_l[30]+fUpwindQuad_l[29]-1.0*fUpwindQuad_l[28]+fUpwindQuad_l[27]-1.0*fUpwindQuad_l[26]+fUpwindQuad_l[25]-1.0*(fUpwindQuad_l[24]+fUpwindQuad_l[23])+fUpwindQuad_l[22]-1.0*fUpwindQuad_l[21]+fUpwindQuad_l[20]-1.0*fUpwindQuad_l[19]+fUpwindQuad_l[18]-1.0*fUpwindQuad_l[17]+fUpwindQuad_l[16]+fUpwindQuad_l[15]-1.0*fUpwindQuad_l[14]+fUpwindQuad_l[13]-1.0*fUpwindQuad_l[12]+fUpwindQuad_l[11]-1.0*fUpwindQuad_l[10]+fUpwindQuad_l[9]-1.0*(fUpwindQuad_l[8]+fUpwindQuad_l[7])+fUpwindQuad_l[6]-1.0*fUpwindQuad_l[5]+fUpwindQuad_l[4]-1.0*fUpwindQuad_l[3]+fUpwindQuad_l[2]-1.0*fUpwindQuad_l[1]+fUpwindQuad_l[0]); 
  fUpwind_l[10] = 0.1767766952966368*(fUpwindQuad_l[31]+fUpwindQuad_l[30]-1.0*(fUpwindQuad_l[29]+fUpwindQuad_l[28])+fUpwindQuad_l[27]+fUpwindQuad_l[26]-1.0*(fUpwindQuad_l[25]+fUpwindQuad_l[24]+fUpwindQuad_l[23]+fUpwindQuad_l[22])+fUpwindQuad_l[21]+fUpwindQuad_l[20]-1.0*(fUpwindQuad_l[19]+fUpwindQuad_l[18])+fUpwindQuad_l[17]+fUpwindQuad_l[16]+fUpwindQuad_l[15]+fUpwindQuad_l[14]-1.0*(fUpwindQuad_l[13]+fUpwindQuad_l[12])+fUpwindQuad_l[11]+fUpwindQuad_l[10]-1.0*(fUpwindQuad_l[9]+fUpwindQuad_l[8]+fUpwindQuad_l[7]+fUpwindQuad_l[6])+fUpwindQuad_l[5]+fUpwindQuad_l[4]-1.0*(fUpwindQuad_l[3]+fUpwindQuad_l[2])+fUpwindQuad_l[1]+fUpwindQuad_l[0]); 
  fUpwind_l[11] = 0.1767766952966368*(fUpwindQuad_l[31]+fUpwindQuad_l[30]+fUpwindQuad_l[29]+fUpwindQuad_l[28]-1.0*(fUpwindQuad_l[27]+fUpwindQuad_l[26]+fUpwindQuad_l[25]+fUpwindQuad_l[24]+fUpwindQuad_l[23]+fUpwindQuad_l[22]+fUpwindQuad_l[21]+fUpwindQuad_l[20])+fUpwindQuad_l[19]+fUpwindQuad_l[18]+fUpwindQuad_l[17]+fUpwindQuad_l[16]+fUpwindQuad_l[15]+fUpwindQuad_l[14]+fUpwindQuad_l[13]+fUpwindQuad_l[12]-1.0*(fUpwindQuad_l[11]+fUpwindQuad_l[10]+fUpwindQuad_l[9]+fUpwindQuad_l[8]+fUpwindQuad_l[7]+fUpwindQuad_l[6]+fUpwindQuad_l[5]+fUpwindQuad_l[4])+fUpwindQuad_l[3]+fUpwindQuad_l[2]+fUpwindQuad_l[1]+fUpwindQuad_l[0]); 
  fUpwind_l[12] = 0.1767766952966368*(fUpwindQuad_l[31]-1.0*fUpwindQuad_l[30]+fUpwindQuad_l[29]-1.0*fUpwindQuad_l[28]+fUpwindQuad_l[27]-1.0*fUpwindQuad_l[26]+fUpwindQuad_l[25]-1.0*fUpwindQuad_l[24]+fUpwindQuad_l[23]-1.0*fUpwindQuad_l[22]+fUpwindQuad_l[21]-1.0*fUpwindQuad_l[20]+fUpwindQuad_l[19]-1.0*fUpwindQuad_l[18]+fUpwindQuad_l[17]-1.0*(fUpwindQuad_l[16]+fUpwindQuad_l[15])+fUpwindQuad_l[14]-1.0*fUpwindQuad_l[13]+fUpwindQuad_l[12]-1.0*fUpwindQuad_l[11]+fUpwindQuad_l[10]-1.0*fUpwindQuad_l[9]+fUpwindQuad_l[8]-1.0*fUpwindQuad_l[7]+fUpwindQuad_l[6]-1.0*fUpwindQuad_l[5]+fUpwindQuad_l[4]-1.0*fUpwindQuad_l[3]+fUpwindQuad_l[2]-1.0*fUpwindQuad_l[1]+fUpwindQuad_l[0]); 
  fUpwind_l[13] = 0.1767766952966368*(fUpwindQuad_l[31]+fUpwindQuad_l[30]-1.0*(fUpwindQuad_l[29]+fUpwindQuad_l[28])+fUpwindQuad_l[27]+fUpwindQuad_l[26]-1.0*(fUpwindQuad_l[25]+fUpwindQuad_l[24])+fUpwindQuad_l[23]+fUpwindQuad_l[22]-1.0*(fUpwindQuad_l[21]+fUpwindQuad_l[20])+fUpwindQuad_l[19]+fUpwindQuad_l[18]-1.0*(fUpwindQuad_l[17]+fUpwindQuad_l[16]+fUpwindQuad_l[15]+fUpwindQuad_l[14])+fUpwindQuad_l[13]+fUpwindQuad_l[12]-1.0*(fUpwindQuad_l[11]+fUpwindQuad_l[10])+fUpwindQuad_l[9]+fUpwindQuad_l[8]-1.0*(fUpwindQuad_l[7]+fUpwindQuad_l[6])+fUpwindQuad_l[5]+fUpwindQuad_l[4]-1.0*(fUpwindQuad_l[3]+fUpwindQuad_l[2])+fUpwindQuad_l[1]+fUpwindQuad_l[0]); 
  fUpwind_l[14] = 0.1767766952966368*(fUpwindQuad_l[31]+fUpwindQuad_l[30]+fUpwindQuad_l[29]+fUpwindQuad_l[28]-1.0*(fUpwindQuad_l[27]+fUpwindQuad_l[26]+fUpwindQuad_l[25]+fUpwindQuad_l[24])+fUpwindQuad_l[23]+fUpwindQuad_l[22]+fUpwindQuad_l[21]+fUpwindQuad_l[20]-1.0*(fUpwindQuad_l[19]+fUpwindQuad_l[18]+fUpwindQuad_l[17]+fUpwindQuad_l[16]+fUpwindQuad_l[15]+fUpwindQuad_l[14]+fUpwindQuad_l[13]+fUpwindQuad_l[12])+fUpwindQuad_l[11]+fUpwindQuad_l[10]+fUpwindQuad_l[9]+fUpwindQuad_l[8]-1.0*(fUpwindQuad_l[7]+fUpwindQuad_l[6]+fUpwindQuad_l[5]+fUpwindQuad_l[4])+fUpwindQuad_l[3]+fUpwindQuad_l[2]+fUpwindQuad_l[1]+fUpwindQuad_l[0]); 
  fUpwind_l[15] = 0.1767766952966368*(fUpwindQuad_l[31]+fUpwindQuad_l[30]+fUpwindQuad_l[29]+fUpwindQuad_l[28]+fUpwindQuad_l[27]+fUpwindQuad_l[26]+fUpwindQuad_l[25]+fUpwindQuad_l[24]-1.0*(fUpwindQuad_l[23]+fUpwindQuad_l[22]+fUpwindQuad_l[21]+fUpwindQuad_l[20]+fUpwindQuad_l[19]+fUpwindQuad_l[18]+fUpwindQuad_l[17]+fUpwindQuad_l[16]+fUpwindQuad_l[15]+fUpwindQuad_l[14]+fUpwindQuad_l[13]+fUpwindQuad_l[12]+fUpwindQuad_l[11]+fUpwindQuad_l[10]+fUpwindQuad_l[9]+fUpwindQuad_l[8])+fUpwindQuad_l[7]+fUpwindQuad_l[6]+fUpwindQuad_l[5]+fUpwindQuad_l[4]+fUpwindQuad_l[3]+fUpwindQuad_l[2]+fUpwindQuad_l[1]+fUpwindQuad_l[0]); 
  fUpwind_l[16] = 0.1767766952966368*(fUpwindQuad_l[31]-1.0*(fUpwindQuad_l[30]+fUpwindQuad_l[29])+fUpwindQuad_l[28]-1.0*fUpwindQuad_l[27]+fUpwindQuad_l[26]+fUpwindQuad_l[25]-1.0*fUpwindQuad_l[24]+fUpwindQuad_l[23]-1.0*(fUpwindQuad_l[22]+fUpwindQuad_l[21])+fUpwindQuad_l[20]-1.0*fUpwindQuad_l[19]+fUpwindQuad_l[18]+fUpwindQuad_l[17]-1.0*fUpwindQuad_l[16]+fUpwindQuad_l[15]-1.0*(fUpwindQuad_l[14]+fUpwindQuad_l[13])+fUpwindQuad_l[12]-1.0*fUpwindQuad_l[11]+fUpwindQuad_l[10]+fUpwindQuad_l[9]-1.0*fUpwindQuad_l[8]+fUpwindQuad_l[7]-1.0*(fUpwindQuad_l[6]+fUpwindQuad_l[5])+fUpwindQuad_l[4]-1.0*fUpwindQuad_l[3]+fUpwindQuad_l[2]+fUpwindQuad_l[1]-1.0*fUpwindQuad_l[0]); 
  fUpwind_l[17] = 0.1767766952966368*(fUpwindQuad_l[31]-1.0*(fUpwindQuad_l[30]+fUpwindQuad_l[29])+fUpwindQuad_l[28]+fUpwindQuad_l[27]-1.0*(fUpwindQuad_l[26]+fUpwindQuad_l[25])+fUpwindQuad_l[24]-1.0*fUpwindQuad_l[23]+fUpwindQuad_l[22]+fUpwindQuad_l[21]-1.0*(fUpwindQuad_l[20]+fUpwindQuad_l[19])+fUpwindQuad_l[18]+fUpwindQuad_l[17]-1.0*fUpwindQuad_l[16]+fUpwindQuad_l[15]-1.0*(fUpwindQuad_l[14]+fUpwindQuad_l[13])+fUpwindQuad_l[12]+fUpwindQuad_l[11]-1.0*(fUpwindQuad_l[10]+fUpwindQuad_l[9])+fUpwindQuad_l[8]-1.0*fUpwindQuad_l[7]+fUpwindQuad_l[6]+fUpwindQuad_l[5]-1.0*(fUpwindQuad_l[4]+fUpwindQuad_l[3])+fUpwindQuad_l[2]+fUpwindQuad_l[1]-1.0*fUpwindQuad_l[0]); 
  fUpwind_l[18] = 0.1767766952966368*(fUpwindQuad_l[31]-1.0*fUpwindQuad_l[30]+fUpwindQuad_l[29]-1.0*(fUpwindQuad_l[28]+fUpwindQuad_l[27])+fUpwindQuad_l[26]-1.0*fUpwindQuad_l[25]+fUpwindQuad_l[24]-1.0*fUpwindQuad_l[23]+fUpwindQuad_l[22]-1.0*fUpwindQuad_l[21]+fUpwindQuad_l[20]+fUpwindQuad_l[19]-1.0*fUpwindQuad_l[18]+fUpwindQuad_l[17]-1.0*fUpwindQuad_l[16]+fUpwindQuad_l[15]-1.0*fUpwindQuad_l[14]+fUpwindQuad_l[13]-1.0*(fUpwindQuad_l[12]+fUpwindQuad_l[11])+fUpwindQuad_l[10]-1.0*fUpwindQuad_l[9]+fUpwindQuad_l[8]-1.0*fUpwindQuad_l[7]+fUpwindQuad_l[6]-1.0*fUpwindQuad_l[5]+fUpwindQuad_l[4]+fUpwindQuad_l[3]-1.0*fUpwindQuad_l[2]+fUpwindQuad_l[1]-1.0*fUpwindQuad_l[0]); 
  fUpwind_l[19] = 0.1767766952966368*(fUpwindQuad_l[31]+fUpwindQuad_l[30]-1.0*(fUpwindQuad_l[29]+fUpwindQuad_l[28]+fUpwindQuad_l[27]+fUpwindQuad_l[26])+fUpwindQuad_l[25]+fUpwindQuad_l[24]-1.0*(fUpwindQuad_l[23]+fUpwindQuad_l[22])+fUpwindQuad_l[21]+fUpwindQuad_l[20]+fUpwindQuad_l[19]+fUpwindQuad_l[18]-1.0*(fUpwindQuad_l[17]+fUpwindQuad_l[16])+fUpwindQuad_l[15]+fUpwindQuad_l[14]-1.0*(fUpwindQuad_l[13]+fUpwindQuad_l[12]+fUpwindQuad_l[11]+fUpwindQuad_l[10])+fUpwindQuad_l[9]+fUpwindQuad_l[8]-1.0*(fUpwindQuad_l[7]+fUpwindQuad_l[6])+fUpwindQuad_l[5]+fUpwindQuad_l[4]+fUpwindQuad_l[3]+fUpwindQuad_l[2]-1.0*(fUpwindQuad_l[1]+fUpwindQuad_l[0])); 
  fUpwind_l[20] = 0.1767766952966368*(fUpwindQuad_l[31]-1.0*(fUpwindQuad_l[30]+fUpwindQuad_l[29])+fUpwindQuad_l[28]+fUpwindQuad_l[27]-1.0*(fUpwindQuad_l[26]+fUpwindQuad_l[25])+fUpwindQuad_l[24]+fUpwindQuad_l[23]-1.0*(fUpwindQuad_l[22]+fUpwindQuad_l[21])+fUpwindQuad_l[20]+fUpwindQuad_l[19]-1.0*(fUpwindQuad_l[18]+fUpwindQuad_l[17])+fUpwindQuad_l[16]-1.0*fUpwindQuad_l[15]+fUpwindQuad_l[14]+fUpwindQuad_l[13]-1.0*(fUpwindQuad_l[12]+fUpwindQuad_l[11])+fUpwindQuad_l[10]+fUpwindQuad_l[9]-1.0*(fUpwindQuad_l[8]+fUpwindQuad_l[7])+fUpwindQuad_l[6]+fUpwindQuad_l[5]-1.0*(fUpwindQuad_l[4]+fUpwindQuad_l[3])+fUpwindQuad_l[2]+fUpwindQuad_l[1]-1.0*fUpwindQuad_l[0]); 
  fUpwind_l[21] = 0.1767766952966368*(fUpwindQuad_l[31]-1.0*fUpwindQuad_l[30]+fUpwindQuad_l[29]-1.0*(fUpwindQuad_l[28]+fUpwindQuad_l[27])+fUpwindQuad_l[26]-1.0*fUpwindQuad_l[25]+fUpwindQuad_l[24]+fUpwindQuad_l[23]-1.0*fUpwindQuad_l[22]+fUpwindQuad_l[21]-1.0*(fUpwindQuad_l[20]+fUpwindQuad_l[19])+fUpwindQuad_l[18]-1.0*fUpwindQuad_l[17]+fUpwindQuad_l[16]-1.0*fUpwindQuad_l[15]+fUpwindQuad_l[14]-1.0*fUpwindQuad_l[13]+fUpwindQuad_l[12]+fUpwindQuad_l[11]-1.0*fUpwindQuad_l[10]+fUpwindQuad_l[9]-1.0*(fUpwindQuad_l[8]+fUpwindQuad_l[7])+fUpwindQuad_l[6]-1.0*fUpwindQuad_l[5]+fUpwindQuad_l[4]+fUpwindQuad_l[3]-1.0*fUpwindQuad_l[2]+fUpwindQuad_l[1]-1.0*fUpwindQuad_l[0]); 
  fUpwind_l[22] = 0.1767766952966368*(fUpwindQuad_l[31]+fUpwindQuad_l[30]-1.0*(fUpwindQuad_l[29]+fUpwindQuad_l[28]+fUpwindQuad_l[27]+fUpwindQuad_l[26])+fUpwindQuad_l[25]+fUpwindQuad_l[24]+fUpwindQuad_l[23]+fUpwindQuad_l[22]-1.0*(fUpwindQuad_l[21]+fUpwindQuad_l[20]+fUpwindQuad_l[19]+fUpwindQuad_l[18])+fUpwindQuad_l[17]+fUpwindQuad_l[16]-1.0*(fUpwindQuad_l[15]+fUpwindQuad_l[14])+fUpwindQuad_l[13]+fUpwindQuad_l[12]+fUpwindQuad_l[11]+fUpwindQuad_l[10]-1.0*(fUpwindQuad_l[9]+fUpwindQuad_l[8]+fUpwindQuad_l[7]+fUpwindQuad_l[6])+fUpwindQuad_l[5]+fUpwindQuad_l[4]+fUpwindQuad_l[3]+fUpwindQuad_l[2]-1.0*(fUpwindQuad_l[1]+fUpwindQuad_l[0])); 
  fUpwind_l[23] = 0.1767766952966368*(fUpwindQuad_l[31]-1.0*fUpwindQuad_l[30]+fUpwindQuad_l[29]-1.0*fUpwindQuad_l[28]+fUpwindQuad_l[27]-1.0*fUpwindQuad_l[26]+fUpwindQuad_l[25]-1.0*(fUpwindQuad_l[24]+fUpwindQuad_l[23])+fUpwindQuad_l[22]-1.0*fUpwindQuad_l[21]+fUpwindQuad_l[20]-1.0*fUpwindQuad_l[19]+fUpwindQuad_l[18]-1.0*fUpwindQuad_l[17]+fUpwindQuad_l[16]-1.0*fUpwindQuad_l[15]+fUpwindQuad_l[14]-1.0*fUpwindQuad_l[13]+fUpwindQuad_l[12]-1.0*fUpwindQuad_l[11]+fUpwindQuad_l[10]-1.0*fUpwindQuad_l[9]+fUpwindQuad_l[8]+fUpwindQuad_l[7]-1.0*fUpwindQuad_l[6]+fUpwindQuad_l[5]-1.0*fUpwindQuad_l[4]+fUpwindQuad_l[3]-1.0*fUpwindQuad_l[2]+fUpwindQuad_l[1]-1.0*fUpwindQuad_l[0]); 
  fUpwind_l[24] = 0.1767766952966368*(fUpwindQuad_l[31]+fUpwindQuad_l[30]-1.0*(fUpwindQuad_l[29]+fUpwindQuad_l[28])+fUpwindQuad_l[27]+fUpwindQuad_l[26]-1.0*(fUpwindQuad_l[25]+fUpwindQuad_l[24]+fUpwindQuad_l[23]+fUpwindQuad_l[22])+fUpwindQuad_l[21]+fUpwindQuad_l[20]-1.0*(fUpwindQuad_l[19]+fUpwindQuad_l[18])+fUpwindQuad_l[17]+fUpwindQuad_l[16]-1.0*(fUpwindQuad_l[15]+fUpwindQuad_l[14])+fUpwindQuad_l[13]+fUpwindQuad_l[12]-1.0*(fUpwindQuad_l[11]+fUpwindQuad_l[10])+fUpwindQuad_l[9]+fUpwindQuad_l[8]+fUpwindQuad_l[7]+fUpwindQuad_l[6]-1.0*(fUpwindQuad_l[5]+fUpwindQuad_l[4])+fUpwindQuad_l[3]+fUpwindQuad_l[2]-1.0*(fUpwindQuad_l[1]+fUpwindQuad_l[0])); 
  fUpwind_l[25] = 0.1767766952966368*(fUpwindQuad_l[31]+fUpwindQuad_l[30]+fUpwindQuad_l[29]+fUpwindQuad_l[28]-1.0*(fUpwindQuad_l[27]+fUpwindQuad_l[26]+fUpwindQuad_l[25]+fUpwindQuad_l[24]+fUpwindQuad_l[23]+fUpwindQuad_l[22]+fUpwindQuad_l[21]+fUpwindQuad_l[20])+fUpwindQuad_l[19]+fUpwindQuad_l[18]+fUpwindQuad_l[17]+fUpwindQuad_l[16]-1.0*(fUpwindQuad_l[15]+fUpwindQuad_l[14]+fUpwindQuad_l[13]+fUpwindQuad_l[12])+fUpwindQuad_l[11]+fUpwindQuad_l[10]+fUpwindQuad_l[9]+fUpwindQuad_l[8]+fUpwindQuad_l[7]+fUpwindQuad_l[6]+fUpwindQuad_l[5]+fUpwindQuad_l[4]-1.0*(fUpwindQuad_l[3]+fUpwindQuad_l[2]+fUpwindQuad_l[1]+fUpwindQuad_l[0])); 
  fUpwind_l[26] = 0.1767766952966368*(fUpwindQuad_l[31]-1.0*(fUpwindQuad_l[30]+fUpwindQuad_l[29])+fUpwindQuad_l[28]-1.0*fUpwindQuad_l[27]+fUpwindQuad_l[26]+fUpwindQuad_l[25]-1.0*(fUpwindQuad_l[24]+fUpwindQuad_l[23])+fUpwindQuad_l[22]+fUpwindQuad_l[21]-1.0*fUpwindQuad_l[20]+fUpwindQuad_l[19]-1.0*(fUpwindQuad_l[18]+fUpwindQuad_l[17])+fUpwindQuad_l[16]+fUpwindQuad_l[15]-1.0*(fUpwindQuad_l[14]+fUpwindQuad_l[13])+fUpwindQuad_l[12]-1.0*fUpwindQuad_l[11]+fUpwindQuad_l[10]+fUpwindQuad_l[9]-1.0*(fUpwindQuad_l[8]+fUpwindQuad_l[7])+fUpwindQuad_l[6]+fUpwindQuad_l[5]-1.0*fUpwindQuad_l[4]+fUpwindQuad_l[3]-1.0*(fUpwindQuad_l[2]+fUpwindQuad_l[1])+fUpwindQuad_l[0]); 
  fUpwind_l[27] = 0.1767766952966368*(fUpwindQuad_l[31]-1.0*(fUpwindQuad_l[30]+fUpwindQuad_l[29])+fUpwindQuad_l[28]-1.0*fUpwindQuad_l[27]+fUpwindQuad_l[26]+fUpwindQuad_l[25]-1.0*fUpwindQuad_l[24]+fUpwindQuad_l[23]-1.0*(fUpwindQuad_l[22]+fUpwindQuad_l[21])+fUpwindQuad_l[20]-1.0*fUpwindQuad_l[19]+fUpwindQuad_l[18]+fUpwindQuad_l[17]-1.0*(fUpwindQuad_l[16]+fUpwindQuad_l[15])+fUpwindQuad_l[14]+fUpwindQuad_l[13]-1.0*fUpwindQuad_l[12]+fUpwindQuad_l[11]-1.0*(fUpwindQuad_l[10]+fUpwindQuad_l[9])+fUpwindQuad_l[8]-1.0*fUpwindQuad_l[7]+fUpwindQuad_l[6]+fUpwindQuad_l[5]-1.0*fUpwindQuad_l[4]+fUpwindQuad_l[3]-1.0*(fUpwindQuad_l[2]+fUpwindQuad_l[1])+fUpwindQuad_l[0]); 
  fUpwind_l[28] = 0.1767766952966368*(fUpwindQuad_l[31]-1.0*(fUpwindQuad_l[30]+fUpwindQuad_l[29])+fUpwindQuad_l[28]+fUpwindQuad_l[27]-1.0*(fUpwindQuad_l[26]+fUpwindQuad_l[25])+fUpwindQuad_l[24]-1.0*fUpwindQuad_l[23]+fUpwindQuad_l[22]+fUpwindQuad_l[21]-1.0*(fUpwindQuad_l[20]+fUpwindQuad_l[19])+fUpwindQuad_l[18]+fUpwindQuad_l[17]-1.0*(fUpwindQuad_l[16]+fUpwindQuad_l[15])+fUpwindQuad_l[14]+fUpwindQuad_l[13]-1.0*(fUpwindQuad_l[12]+fUpwindQuad_l[11])+fUpwindQuad_l[10]+fUpwindQuad_l[9]-1.0*fUpwindQuad_l[8]+fUpwindQuad_l[7]-1.0*(fUpwindQuad_l[6]+fUpwindQuad_l[5])+fUpwindQuad_l[4]+fUpwindQuad_l[3]-1.0*(fUpwindQuad_l[2]+fUpwindQuad_l[1])+fUpwindQuad_l[0]); 
  fUpwind_l[29] = 0.1767766952966368*(fUpwindQuad_l[31]-1.0*fUpwindQuad_l[30]+fUpwindQuad_l[29]-1.0*(fUpwindQuad_l[28]+fUpwindQuad_l[27])+fUpwindQuad_l[26]-1.0*fUpwindQuad_l[25]+fUpwindQuad_l[24]-1.0*fUpwindQuad_l[23]+fUpwindQuad_l[22]-1.0*fUpwindQuad_l[21]+fUpwindQuad_l[20]+fUpwindQuad_l[19]-1.0*fUpwindQuad_l[18]+fUpwindQuad_l[17]-1.0*(fUpwindQuad_l[16]+fUpwindQuad_l[15])+fUpwindQuad_l[14]-1.0*fUpwindQuad_l[13]+fUpwindQuad_l[12]+fUpwindQuad_l[11]-1.0*fUpwindQuad_l[10]+fUpwindQuad_l[9]-1.0*fUpwindQuad_l[8]+fUpwindQuad_l[7]-1.0*fUpwindQuad_l[6]+fUpwindQuad_l[5]-1.0*(fUpwindQuad_l[4]+fUpwindQuad_l[3])+fUpwindQuad_l[2]-1.0*fUpwindQuad_l[1]+fUpwindQuad_l[0]); 
  fUpwind_l[30] = 0.1767766952966368*(fUpwindQuad_l[31]+fUpwindQuad_l[30]-1.0*(fUpwindQuad_l[29]+fUpwindQuad_l[28]+fUpwindQuad_l[27]+fUpwindQuad_l[26])+fUpwindQuad_l[25]+fUpwindQuad_l[24]-1.0*(fUpwindQuad_l[23]+fUpwindQuad_l[22])+fUpwindQuad_l[21]+fUpwindQuad_l[20]+fUpwindQuad_l[19]+fUpwindQuad_l[18]-1.0*(fUpwindQuad_l[17]+fUpwindQuad_l[16]+fUpwindQuad_l[15]+fUpwindQuad_l[14])+fUpwindQuad_l[13]+fUpwindQuad_l[12]+fUpwindQuad_l[11]+fUpwindQuad_l[10]-1.0*(fUpwindQuad_l[9]+fUpwindQuad_l[8])+fUpwindQuad_l[7]+fUpwindQuad_l[6]-1.0*(fUpwindQuad_l[5]+fUpwindQuad_l[4]+fUpwindQuad_l[3]+fUpwindQuad_l[2])+fUpwindQuad_l[1]+fUpwindQuad_l[0]); 
  fUpwind_l[31] = 0.1767766952966368*(fUpwindQuad_l[31]-1.0*(fUpwindQuad_l[30]+fUpwindQuad_l[29])+fUpwindQuad_l[28]-1.0*fUpwindQuad_l[27]+fUpwindQuad_l[26]+fUpwindQuad_l[25]-1.0*(fUpwindQuad_l[24]+fUpwindQuad_l[23])+fUpwindQuad_l[22]+fUpwindQuad_l[21]-1.0*fUpwindQuad_l[20]+fUpwindQuad_l[19]-1.0*(fUpwindQuad_l[18]+fUpwindQuad_l[17])+fUpwindQuad_l[16]-1.0*fUpwindQuad_l[15]+fUpwindQuad_l[14]+fUpwindQuad_l[13]-1.0*fUpwindQuad_l[12]+fUpwindQuad_l[11]-1.0*(fUpwindQuad_l[10]+fUpwindQuad_l[9])+fUpwindQuad_l[8]+fUpwindQuad_l[7]-1.0*(fUpwindQuad_l[6]+fUpwindQuad_l[5])+fUpwindQuad_l[4]-1.0*fUpwindQuad_l[3]+fUpwindQuad_l[2]+fUpwindQuad_l[1]-1.0*fUpwindQuad_l[0]); 

  fUpwind_r[0] = 0.1767766952966368*(fUpwindQuad_r[31]+fUpwindQuad_r[30]+fUpwindQuad_r[29]+fUpwindQuad_r[28]+fUpwindQuad_r[27]+fUpwindQuad_r[26]+fUpwindQuad_r[25]+fUpwindQuad_r[24]+fUpwindQuad_r[23]+fUpwindQuad_r[22]+fUpwindQuad_r[21]+fUpwindQuad_r[20]+fUpwindQuad_r[19]+fUpwindQuad_r[18]+fUpwindQuad_r[17]+fUpwindQuad_r[16]+fUpwindQuad_r[15]+fUpwindQuad_r[14]+fUpwindQuad_r[13]+fUpwindQuad_r[12]+fUpwindQuad_r[11]+fUpwindQuad_r[10]+fUpwindQuad_r[9]+fUpwindQuad_r[8]+fUpwindQuad_r[7]+fUpwindQuad_r[6]+fUpwindQuad_r[5]+fUpwindQuad_r[4]+fUpwindQuad_r[3]+fUpwindQuad_r[2]+fUpwindQuad_r[1]+fUpwindQuad_r[0]); 
  fUpwind_r[1] = 0.1767766952966368*(fUpwindQuad_r[31]-1.0*fUpwindQuad_r[30]+fUpwindQuad_r[29]-1.0*fUpwindQuad_r[28]+fUpwindQuad_r[27]-1.0*fUpwindQuad_r[26]+fUpwindQuad_r[25]-1.0*fUpwindQuad_r[24]+fUpwindQuad_r[23]-1.0*fUpwindQuad_r[22]+fUpwindQuad_r[21]-1.0*fUpwindQuad_r[20]+fUpwindQuad_r[19]-1.0*fUpwindQuad_r[18]+fUpwindQuad_r[17]-1.0*fUpwindQuad_r[16]+fUpwindQuad_r[15]-1.0*fUpwindQuad_r[14]+fUpwindQuad_r[13]-1.0*fUpwindQuad_r[12]+fUpwindQuad_r[11]-1.0*fUpwindQuad_r[10]+fUpwindQuad_r[9]-1.0*fUpwindQuad_r[8]+fUpwindQuad_r[7]-1.0*fUpwindQuad_r[6]+fUpwindQuad_r[5]-1.0*fUpwindQuad_r[4]+fUpwindQuad_r[3]-1.0*fUpwindQuad_r[2]+fUpwindQuad_r[1]-1.0*fUpwindQuad_r[0]); 
  fUpwind_r[2] = 0.1767766952966368*(fUpwindQuad_r[31]+fUpwindQuad_r[30]-1.0*(fUpwindQuad_r[29]+fUpwindQuad_r[28])+fUpwindQuad_r[27]+fUpwindQuad_r[26]-1.0*(fUpwindQuad_r[25]+fUpwindQuad_r[24])+fUpwindQuad_r[23]+fUpwindQuad_r[22]-1.0*(fUpwindQuad_r[21]+fUpwindQuad_r[20])+fUpwindQuad_r[19]+fUpwindQuad_r[18]-1.0*(fUpwindQuad_r[17]+fUpwindQuad_r[16])+fUpwindQuad_r[15]+fUpwindQuad_r[14]-1.0*(fUpwindQuad_r[13]+fUpwindQuad_r[12])+fUpwindQuad_r[11]+fUpwindQuad_r[10]-1.0*(fUpwindQuad_r[9]+fUpwindQuad_r[8])+fUpwindQuad_r[7]+fUpwindQuad_r[6]-1.0*(fUpwindQuad_r[5]+fUpwindQuad_r[4])+fUpwindQuad_r[3]+fUpwindQuad_r[2]-1.0*(fUpwindQuad_r[1]+fUpwindQuad_r[0])); 
  fUpwind_r[3] = 0.1767766952966368*(fUpwindQuad_r[31]+fUpwindQuad_r[30]+fUpwindQuad_r[29]+fUpwindQuad_r[28]-1.0*(fUpwindQuad_r[27]+fUpwindQuad_r[26]+fUpwindQuad_r[25]+fUpwindQuad_r[24])+fUpwindQuad_r[23]+fUpwindQuad_r[22]+fUpwindQuad_r[21]+fUpwindQuad_r[20]-1.0*(fUpwindQuad_r[19]+fUpwindQuad_r[18]+fUpwindQuad_r[17]+fUpwindQuad_r[16])+fUpwindQuad_r[15]+fUpwindQuad_r[14]+fUpwindQuad_r[13]+fUpwindQuad_r[12]-1.0*(fUpwindQuad_r[11]+fUpwindQuad_r[10]+fUpwindQuad_r[9]+fUpwindQuad_r[8])+fUpwindQuad_r[7]+fUpwindQuad_r[6]+fUpwindQuad_r[5]+fUpwindQuad_r[4]-1.0*(fUpwindQuad_r[3]+fUpwindQuad_r[2]+fUpwindQuad_r[1]+fUpwindQuad_r[0])); 
  fUpwind_r[4] = 0.1767766952966368*(fUpwindQuad_r[31]+fUpwindQuad_r[30]+fUpwindQuad_r[29]+fUpwindQuad_r[28]+fUpwindQuad_r[27]+fUpwindQuad_r[26]+fUpwindQuad_r[25]+fUpwindQuad_r[24]-1.0*(fUpwindQuad_r[23]+fUpwindQuad_r[22]+fUpwindQuad_r[21]+fUpwindQuad_r[20]+fUpwindQuad_r[19]+fUpwindQuad_r[18]+fUpwindQuad_r[17]+fUpwindQuad_r[16])+fUpwindQuad_r[15]+fUpwindQuad_r[14]+fUpwindQuad_r[13]+fUpwindQuad_r[12]+fUpwindQuad_r[11]+fUpwindQuad_r[10]+fUpwindQuad_r[9]+fUpwindQuad_r[8]-1.0*(fUpwindQuad_r[7]+fUpwindQuad_r[6]+fUpwindQuad_r[5]+fUpwindQuad_r[4]+fUpwindQuad_r[3]+fUpwindQuad_r[2]+fUpwindQuad_r[1]+fUpwindQuad_r[0])); 
  fUpwind_r[5] = 0.1767766952966368*(fUpwindQuad_r[31]+fUpwindQuad_r[30]+fUpwindQuad_r[29]+fUpwindQuad_r[28]+fUpwindQuad_r[27]+fUpwindQuad_r[26]+fUpwindQuad_r[25]+fUpwindQuad_r[24]+fUpwindQuad_r[23]+fUpwindQuad_r[22]+fUpwindQuad_r[21]+fUpwindQuad_r[20]+fUpwindQuad_r[19]+fUpwindQuad_r[18]+fUpwindQuad_r[17]+fUpwindQuad_r[16]-1.0*(fUpwindQuad_r[15]+fUpwindQuad_r[14]+fUpwindQuad_r[13]+fUpwindQuad_r[12]+fUpwindQuad_r[11]+fUpwindQuad_r[10]+fUpwindQuad_r[9]+fUpwindQuad_r[8]+fUpwindQuad_r[7]+fUpwindQuad_r[6]+fUpwindQuad_r[5]+fUpwindQuad_r[4]+fUpwindQuad_r[3]+fUpwindQuad_r[2]+fUpwindQuad_r[1]+fUpwindQuad_r[0])); 
  fUpwind_r[6] = 0.1767766952966368*(fUpwindQuad_r[31]-1.0*(fUpwindQuad_r[30]+fUpwindQuad_r[29])+fUpwindQuad_r[28]+fUpwindQuad_r[27]-1.0*(fUpwindQuad_r[26]+fUpwindQuad_r[25])+fUpwindQuad_r[24]+fUpwindQuad_r[23]-1.0*(fUpwindQuad_r[22]+fUpwindQuad_r[21])+fUpwindQuad_r[20]+fUpwindQuad_r[19]-1.0*(fUpwindQuad_r[18]+fUpwindQuad_r[17])+fUpwindQuad_r[16]+fUpwindQuad_r[15]-1.0*(fUpwindQuad_r[14]+fUpwindQuad_r[13])+fUpwindQuad_r[12]+fUpwindQuad_r[11]-1.0*(fUpwindQuad_r[10]+fUpwindQuad_r[9])+fUpwindQuad_r[8]+fUpwindQuad_r[7]-1.0*(fUpwindQuad_r[6]+fUpwindQuad_r[5])+fUpwindQuad_r[4]+fUpwindQuad_r[3]-1.0*(fUpwindQuad_r[2]+fUpwindQuad_r[1])+fUpwindQuad_r[0]); 
  fUpwind_r[7] = 0.1767766952966368*(fUpwindQuad_r[31]-1.0*fUpwindQuad_r[30]+fUpwindQuad_r[29]-1.0*(fUpwindQuad_r[28]+fUpwindQuad_r[27])+fUpwindQuad_r[26]-1.0*fUpwindQuad_r[25]+fUpwindQuad_r[24]+fUpwindQuad_r[23]-1.0*fUpwindQuad_r[22]+fUpwindQuad_r[21]-1.0*(fUpwindQuad_r[20]+fUpwindQuad_r[19])+fUpwindQuad_r[18]-1.0*fUpwindQuad_r[17]+fUpwindQuad_r[16]+fUpwindQuad_r[15]-1.0*fUpwindQuad_r[14]+fUpwindQuad_r[13]-1.0*(fUpwindQuad_r[12]+fUpwindQuad_r[11])+fUpwindQuad_r[10]-1.0*fUpwindQuad_r[9]+fUpwindQuad_r[8]+fUpwindQuad_r[7]-1.0*fUpwindQuad_r[6]+fUpwindQuad_r[5]-1.0*(fUpwindQuad_r[4]+fUpwindQuad_r[3])+fUpwindQuad_r[2]-1.0*fUpwindQuad_r[1]+fUpwindQuad_r[0]); 
  fUpwind_r[8] = 0.1767766952966368*(fUpwindQuad_r[31]+fUpwindQuad_r[30]-1.0*(fUpwindQuad_r[29]+fUpwindQuad_r[28]+fUpwindQuad_r[27]+fUpwindQuad_r[26])+fUpwindQuad_r[25]+fUpwindQuad_r[24]+fUpwindQuad_r[23]+fUpwindQuad_r[22]-1.0*(fUpwindQuad_r[21]+fUpwindQuad_r[20]+fUpwindQuad_r[19]+fUpwindQuad_r[18])+fUpwindQuad_r[17]+fUpwindQuad_r[16]+fUpwindQuad_r[15]+fUpwindQuad_r[14]-1.0*(fUpwindQuad_r[13]+fUpwindQuad_r[12]+fUpwindQuad_r[11]+fUpwindQuad_r[10])+fUpwindQuad_r[9]+fUpwindQuad_r[8]+fUpwindQuad_r[7]+fUpwindQuad_r[6]-1.0*(fUpwindQuad_r[5]+fUpwindQuad_r[4]+fUpwindQuad_r[3]+fUpwindQuad_r[2])+fUpwindQuad_r[1]+fUpwindQuad_r[0]); 
  fUpwind_r[9] = 0.1767766952966368*(fUpwindQuad_r[31]-1.0*fUpwindQuad_r[30]+fUpwindQuad_r[29]-1.0*fUpwindQuad_r[28]+fUpwindQuad_r[27]-1.0*fUpwindQuad_r[26]+fUpwindQuad_r[25]-1.0*(fUpwindQuad_r[24]+fUpwindQuad_r[23])+fUpwindQuad_r[22]-1.0*fUpwindQuad_r[21]+fUpwindQuad_r[20]-1.0*fUpwindQuad_r[19]+fUpwindQuad_r[18]-1.0*fUpwindQuad_r[17]+fUpwindQuad_r[16]+fUpwindQuad_r[15]-1.0*fUpwindQuad_r[14]+fUpwindQuad_r[13]-1.0*fUpwindQuad_r[12]+fUpwindQuad_r[11]-1.0*fUpwindQuad_r[10]+fUpwindQuad_r[9]-1.0*(fUpwindQuad_r[8]+fUpwindQuad_r[7])+fUpwindQuad_r[6]-1.0*fUpwindQuad_r[5]+fUpwindQuad_r[4]-1.0*fUpwindQuad_r[3]+fUpwindQuad_r[2]-1.0*fUpwindQuad_r[1]+fUpwindQuad_r[0]); 
  fUpwind_r[10] = 0.1767766952966368*(fUpwindQuad_r[31]+fUpwindQuad_r[30]-1.0*(fUpwindQuad_r[29]+fUpwindQuad_r[28])+fUpwindQuad_r[27]+fUpwindQuad_r[26]-1.0*(fUpwindQuad_r[25]+fUpwindQuad_r[24]+fUpwindQuad_r[23]+fUpwindQuad_r[22])+fUpwindQuad_r[21]+fUpwindQuad_r[20]-1.0*(fUpwindQuad_r[19]+fUpwindQuad_r[18])+fUpwindQuad_r[17]+fUpwindQuad_r[16]+fUpwindQuad_r[15]+fUpwindQuad_r[14]-1.0*(fUpwindQuad_r[13]+fUpwindQuad_r[12])+fUpwindQuad_r[11]+fUpwindQuad_r[10]-1.0*(fUpwindQuad_r[9]+fUpwindQuad_r[8]+fUpwindQuad_r[7]+fUpwindQuad_r[6])+fUpwindQuad_r[5]+fUpwindQuad_r[4]-1.0*(fUpwindQuad_r[3]+fUpwindQuad_r[2])+fUpwindQuad_r[1]+fUpwindQuad_r[0]); 
  fUpwind_r[11] = 0.1767766952966368*(fUpwindQuad_r[31]+fUpwindQuad_r[30]+fUpwindQuad_r[29]+fUpwindQuad_r[28]-1.0*(fUpwindQuad_r[27]+fUpwindQuad_r[26]+fUpwindQuad_r[25]+fUpwindQuad_r[24]+fUpwindQuad_r[23]+fUpwindQuad_r[22]+fUpwindQuad_r[21]+fUpwindQuad_r[20])+fUpwindQuad_r[19]+fUpwindQuad_r[18]+fUpwindQuad_r[17]+fUpwindQuad_r[16]+fUpwindQuad_r[15]+fUpwindQuad_r[14]+fUpwindQuad_r[13]+fUpwindQuad_r[12]-1.0*(fUpwindQuad_r[11]+fUpwindQuad_r[10]+fUpwindQuad_r[9]+fUpwindQuad_r[8]+fUpwindQuad_r[7]+fUpwindQuad_r[6]+fUpwindQuad_r[5]+fUpwindQuad_r[4])+fUpwindQuad_r[3]+fUpwindQuad_r[2]+fUpwindQuad_r[1]+fUpwindQuad_r[0]); 
  fUpwind_r[12] = 0.1767766952966368*(fUpwindQuad_r[31]-1.0*fUpwindQuad_r[30]+fUpwindQuad_r[29]-1.0*fUpwindQuad_r[28]+fUpwindQuad_r[27]-1.0*fUpwindQuad_r[26]+fUpwindQuad_r[25]-1.0*fUpwindQuad_r[24]+fUpwindQuad_r[23]-1.0*fUpwindQuad_r[22]+fUpwindQuad_r[21]-1.0*fUpwindQuad_r[20]+fUpwindQuad_r[19]-1.0*fUpwindQuad_r[18]+fUpwindQuad_r[17]-1.0*(fUpwindQuad_r[16]+fUpwindQuad_r[15])+fUpwindQuad_r[14]-1.0*fUpwindQuad_r[13]+fUpwindQuad_r[12]-1.0*fUpwindQuad_r[11]+fUpwindQuad_r[10]-1.0*fUpwindQuad_r[9]+fUpwindQuad_r[8]-1.0*fUpwindQuad_r[7]+fUpwindQuad_r[6]-1.0*fUpwindQuad_r[5]+fUpwindQuad_r[4]-1.0*fUpwindQuad_r[3]+fUpwindQuad_r[2]-1.0*fUpwindQuad_r[1]+fUpwindQuad_r[0]); 
  fUpwind_r[13] = 0.1767766952966368*(fUpwindQuad_r[31]+fUpwindQuad_r[30]-1.0*(fUpwindQuad_r[29]+fUpwindQuad_r[28])+fUpwindQuad_r[27]+fUpwindQuad_r[26]-1.0*(fUpwindQuad_r[25]+fUpwindQuad_r[24])+fUpwindQuad_r[23]+fUpwindQuad_r[22]-1.0*(fUpwindQuad_r[21]+fUpwindQuad_r[20])+fUpwindQuad_r[19]+fUpwindQuad_r[18]-1.0*(fUpwindQuad_r[17]+fUpwindQuad_r[16]+fUpwindQuad_r[15]+fUpwindQuad_r[14])+fUpwindQuad_r[13]+fUpwindQuad_r[12]-1.0*(fUpwindQuad_r[11]+fUpwindQuad_r[10])+fUpwindQuad_r[9]+fUpwindQuad_r[8]-1.0*(fUpwindQuad_r[7]+fUpwindQuad_r[6])+fUpwindQuad_r[5]+fUpwindQuad_r[4]-1.0*(fUpwindQuad_r[3]+fUpwindQuad_r[2])+fUpwindQuad_r[1]+fUpwindQuad_r[0]); 
  fUpwind_r[14] = 0.1767766952966368*(fUpwindQuad_r[31]+fUpwindQuad_r[30]+fUpwindQuad_r[29]+fUpwindQuad_r[28]-1.0*(fUpwindQuad_r[27]+fUpwindQuad_r[26]+fUpwindQuad_r[25]+fUpwindQuad_r[24])+fUpwindQuad_r[23]+fUpwindQuad_r[22]+fUpwindQuad_r[21]+fUpwindQuad_r[20]-1.0*(fUpwindQuad_r[19]+fUpwindQuad_r[18]+fUpwindQuad_r[17]+fUpwindQuad_r[16]+fUpwindQuad_r[15]+fUpwindQuad_r[14]+fUpwindQuad_r[13]+fUpwindQuad_r[12])+fUpwindQuad_r[11]+fUpwindQuad_r[10]+fUpwindQuad_r[9]+fUpwindQuad_r[8]-1.0*(fUpwindQuad_r[7]+fUpwindQuad_r[6]+fUpwindQuad_r[5]+fUpwindQuad_r[4])+fUpwindQuad_r[3]+fUpwindQuad_r[2]+fUpwindQuad_r[1]+fUpwindQuad_r[0]); 
  fUpwind_r[15] = 0.1767766952966368*(fUpwindQuad_r[31]+fUpwindQuad_r[30]+fUpwindQuad_r[29]+fUpwindQuad_r[28]+fUpwindQuad_r[27]+fUpwindQuad_r[26]+fUpwindQuad_r[25]+fUpwindQuad_r[24]-1.0*(fUpwindQuad_r[23]+fUpwindQuad_r[22]+fUpwindQuad_r[21]+fUpwindQuad_r[20]+fUpwindQuad_r[19]+fUpwindQuad_r[18]+fUpwindQuad_r[17]+fUpwindQuad_r[16]+fUpwindQuad_r[15]+fUpwindQuad_r[14]+fUpwindQuad_r[13]+fUpwindQuad_r[12]+fUpwindQuad_r[11]+fUpwindQuad_r[10]+fUpwindQuad_r[9]+fUpwindQuad_r[8])+fUpwindQuad_r[7]+fUpwindQuad_r[6]+fUpwindQuad_r[5]+fUpwindQuad_r[4]+fUpwindQuad_r[3]+fUpwindQuad_r[2]+fUpwindQuad_r[1]+fUpwindQuad_r[0]); 
  fUpwind_r[16] = 0.1767766952966368*(fUpwindQuad_r[31]-1.0*(fUpwindQuad_r[30]+fUpwindQuad_r[29])+fUpwindQuad_r[28]-1.0*fUpwindQuad_r[27]+fUpwindQuad_r[26]+fUpwindQuad_r[25]-1.0*fUpwindQuad_r[24]+fUpwindQuad_r[23]-1.0*(fUpwindQuad_r[22]+fUpwindQuad_r[21])+fUpwindQuad_r[20]-1.0*fUpwindQuad_r[19]+fUpwindQuad_r[18]+fUpwindQuad_r[17]-1.0*fUpwindQuad_r[16]+fUpwindQuad_r[15]-1.0*(fUpwindQuad_r[14]+fUpwindQuad_r[13])+fUpwindQuad_r[12]-1.0*fUpwindQuad_r[11]+fUpwindQuad_r[10]+fUpwindQuad_r[9]-1.0*fUpwindQuad_r[8]+fUpwindQuad_r[7]-1.0*(fUpwindQuad_r[6]+fUpwindQuad_r[5])+fUpwindQuad_r[4]-1.0*fUpwindQuad_r[3]+fUpwindQuad_r[2]+fUpwindQuad_r[1]-1.0*fUpwindQuad_r[0]); 
  fUpwind_r[17] = 0.1767766952966368*(fUpwindQuad_r[31]-1.0*(fUpwindQuad_r[30]+fUpwindQuad_r[29])+fUpwindQuad_r[28]+fUpwindQuad_r[27]-1.0*(fUpwindQuad_r[26]+fUpwindQuad_r[25])+fUpwindQuad_r[24]-1.0*fUpwindQuad_r[23]+fUpwindQuad_r[22]+fUpwindQuad_r[21]-1.0*(fUpwindQuad_r[20]+fUpwindQuad_r[19])+fUpwindQuad_r[18]+fUpwindQuad_r[17]-1.0*fUpwindQuad_r[16]+fUpwindQuad_r[15]-1.0*(fUpwindQuad_r[14]+fUpwindQuad_r[13])+fUpwindQuad_r[12]+fUpwindQuad_r[11]-1.0*(fUpwindQuad_r[10]+fUpwindQuad_r[9])+fUpwindQuad_r[8]-1.0*fUpwindQuad_r[7]+fUpwindQuad_r[6]+fUpwindQuad_r[5]-1.0*(fUpwindQuad_r[4]+fUpwindQuad_r[3])+fUpwindQuad_r[2]+fUpwindQuad_r[1]-1.0*fUpwindQuad_r[0]); 
  fUpwind_r[18] = 0.1767766952966368*(fUpwindQuad_r[31]-1.0*fUpwindQuad_r[30]+fUpwindQuad_r[29]-1.0*(fUpwindQuad_r[28]+fUpwindQuad_r[27])+fUpwindQuad_r[26]-1.0*fUpwindQuad_r[25]+fUpwindQuad_r[24]-1.0*fUpwindQuad_r[23]+fUpwindQuad_r[22]-1.0*fUpwindQuad_r[21]+fUpwindQuad_r[20]+fUpwindQuad_r[19]-1.0*fUpwindQuad_r[18]+fUpwindQuad_r[17]-1.0*fUpwindQuad_r[16]+fUpwindQuad_r[15]-1.0*fUpwindQuad_r[14]+fUpwindQuad_r[13]-1.0*(fUpwindQuad_r[12]+fUpwindQuad_r[11])+fUpwindQuad_r[10]-1.0*fUpwindQuad_r[9]+fUpwindQuad_r[8]-1.0*fUpwindQuad_r[7]+fUpwindQuad_r[6]-1.0*fUpwindQuad_r[5]+fUpwindQuad_r[4]+fUpwindQuad_r[3]-1.0*fUpwindQuad_r[2]+fUpwindQuad_r[1]-1.0*fUpwindQuad_r[0]); 
  fUpwind_r[19] = 0.1767766952966368*(fUpwindQuad_r[31]+fUpwindQuad_r[30]-1.0*(fUpwindQuad_r[29]+fUpwindQuad_r[28]+fUpwindQuad_r[27]+fUpwindQuad_r[26])+fUpwindQuad_r[25]+fUpwindQuad_r[24]-1.0*(fUpwindQuad_r[23]+fUpwindQuad_r[22])+fUpwindQuad_r[21]+fUpwindQuad_r[20]+fUpwindQuad_r[19]+fUpwindQuad_r[18]-1.0*(fUpwindQuad_r[17]+fUpwindQuad_r[16])+fUpwindQuad_r[15]+fUpwindQuad_r[14]-1.0*(fUpwindQuad_r[13]+fUpwindQuad_r[12]+fUpwindQuad_r[11]+fUpwindQuad_r[10])+fUpwindQuad_r[9]+fUpwindQuad_r[8]-1.0*(fUpwindQuad_r[7]+fUpwindQuad_r[6])+fUpwindQuad_r[5]+fUpwindQuad_r[4]+fUpwindQuad_r[3]+fUpwindQuad_r[2]-1.0*(fUpwindQuad_r[1]+fUpwindQuad_r[0])); 
  fUpwind_r[20] = 0.1767766952966368*(fUpwindQuad_r[31]-1.0*(fUpwindQuad_r[30]+fUpwindQuad_r[29])+fUpwindQuad_r[28]+fUpwindQuad_r[27]-1.0*(fUpwindQuad_r[26]+fUpwindQuad_r[25])+fUpwindQuad_r[24]+fUpwindQuad_r[23]-1.0*(fUpwindQuad_r[22]+fUpwindQuad_r[21])+fUpwindQuad_r[20]+fUpwindQuad_r[19]-1.0*(fUpwindQuad_r[18]+fUpwindQuad_r[17])+fUpwindQuad_r[16]-1.0*fUpwindQuad_r[15]+fUpwindQuad_r[14]+fUpwindQuad_r[13]-1.0*(fUpwindQuad_r[12]+fUpwindQuad_r[11])+fUpwindQuad_r[10]+fUpwindQuad_r[9]-1.0*(fUpwindQuad_r[8]+fUpwindQuad_r[7])+fUpwindQuad_r[6]+fUpwindQuad_r[5]-1.0*(fUpwindQuad_r[4]+fUpwindQuad_r[3])+fUpwindQuad_r[2]+fUpwindQuad_r[1]-1.0*fUpwindQuad_r[0]); 
  fUpwind_r[21] = 0.1767766952966368*(fUpwindQuad_r[31]-1.0*fUpwindQuad_r[30]+fUpwindQuad_r[29]-1.0*(fUpwindQuad_r[28]+fUpwindQuad_r[27])+fUpwindQuad_r[26]-1.0*fUpwindQuad_r[25]+fUpwindQuad_r[24]+fUpwindQuad_r[23]-1.0*fUpwindQuad_r[22]+fUpwindQuad_r[21]-1.0*(fUpwindQuad_r[20]+fUpwindQuad_r[19])+fUpwindQuad_r[18]-1.0*fUpwindQuad_r[17]+fUpwindQuad_r[16]-1.0*fUpwindQuad_r[15]+fUpwindQuad_r[14]-1.0*fUpwindQuad_r[13]+fUpwindQuad_r[12]+fUpwindQuad_r[11]-1.0*fUpwindQuad_r[10]+fUpwindQuad_r[9]-1.0*(fUpwindQuad_r[8]+fUpwindQuad_r[7])+fUpwindQuad_r[6]-1.0*fUpwindQuad_r[5]+fUpwindQuad_r[4]+fUpwindQuad_r[3]-1.0*fUpwindQuad_r[2]+fUpwindQuad_r[1]-1.0*fUpwindQuad_r[0]); 
  fUpwind_r[22] = 0.1767766952966368*(fUpwindQuad_r[31]+fUpwindQuad_r[30]-1.0*(fUpwindQuad_r[29]+fUpwindQuad_r[28]+fUpwindQuad_r[27]+fUpwindQuad_r[26])+fUpwindQuad_r[25]+fUpwindQuad_r[24]+fUpwindQuad_r[23]+fUpwindQuad_r[22]-1.0*(fUpwindQuad_r[21]+fUpwindQuad_r[20]+fUpwindQuad_r[19]+fUpwindQuad_r[18])+fUpwindQuad_r[17]+fUpwindQuad_r[16]-1.0*(fUpwindQuad_r[15]+fUpwindQuad_r[14])+fUpwindQuad_r[13]+fUpwindQuad_r[12]+fUpwindQuad_r[11]+fUpwindQuad_r[10]-1.0*(fUpwindQuad_r[9]+fUpwindQuad_r[8]+fUpwindQuad_r[7]+fUpwindQuad_r[6])+fUpwindQuad_r[5]+fUpwindQuad_r[4]+fUpwindQuad_r[3]+fUpwindQuad_r[2]-1.0*(fUpwindQuad_r[1]+fUpwindQuad_r[0])); 
  fUpwind_r[23] = 0.1767766952966368*(fUpwindQuad_r[31]-1.0*fUpwindQuad_r[30]+fUpwindQuad_r[29]-1.0*fUpwindQuad_r[28]+fUpwindQuad_r[27]-1.0*fUpwindQuad_r[26]+fUpwindQuad_r[25]-1.0*(fUpwindQuad_r[24]+fUpwindQuad_r[23])+fUpwindQuad_r[22]-1.0*fUpwindQuad_r[21]+fUpwindQuad_r[20]-1.0*fUpwindQuad_r[19]+fUpwindQuad_r[18]-1.0*fUpwindQuad_r[17]+fUpwindQuad_r[16]-1.0*fUpwindQuad_r[15]+fUpwindQuad_r[14]-1.0*fUpwindQuad_r[13]+fUpwindQuad_r[12]-1.0*fUpwindQuad_r[11]+fUpwindQuad_r[10]-1.0*fUpwindQuad_r[9]+fUpwindQuad_r[8]+fUpwindQuad_r[7]-1.0*fUpwindQuad_r[6]+fUpwindQuad_r[5]-1.0*fUpwindQuad_r[4]+fUpwindQuad_r[3]-1.0*fUpwindQuad_r[2]+fUpwindQuad_r[1]-1.0*fUpwindQuad_r[0]); 
  fUpwind_r[24] = 0.1767766952966368*(fUpwindQuad_r[31]+fUpwindQuad_r[30]-1.0*(fUpwindQuad_r[29]+fUpwindQuad_r[28])+fUpwindQuad_r[27]+fUpwindQuad_r[26]-1.0*(fUpwindQuad_r[25]+fUpwindQuad_r[24]+fUpwindQuad_r[23]+fUpwindQuad_r[22])+fUpwindQuad_r[21]+fUpwindQuad_r[20]-1.0*(fUpwindQuad_r[19]+fUpwindQuad_r[18])+fUpwindQuad_r[17]+fUpwindQuad_r[16]-1.0*(fUpwindQuad_r[15]+fUpwindQuad_r[14])+fUpwindQuad_r[13]+fUpwindQuad_r[12]-1.0*(fUpwindQuad_r[11]+fUpwindQuad_r[10])+fUpwindQuad_r[9]+fUpwindQuad_r[8]+fUpwindQuad_r[7]+fUpwindQuad_r[6]-1.0*(fUpwindQuad_r[5]+fUpwindQuad_r[4])+fUpwindQuad_r[3]+fUpwindQuad_r[2]-1.0*(fUpwindQuad_r[1]+fUpwindQuad_r[0])); 
  fUpwind_r[25] = 0.1767766952966368*(fUpwindQuad_r[31]+fUpwindQuad_r[30]+fUpwindQuad_r[29]+fUpwindQuad_r[28]-1.0*(fUpwindQuad_r[27]+fUpwindQuad_r[26]+fUpwindQuad_r[25]+fUpwindQuad_r[24]+fUpwindQuad_r[23]+fUpwindQuad_r[22]+fUpwindQuad_r[21]+fUpwindQuad_r[20])+fUpwindQuad_r[19]+fUpwindQuad_r[18]+fUpwindQuad_r[17]+fUpwindQuad_r[16]-1.0*(fUpwindQuad_r[15]+fUpwindQuad_r[14]+fUpwindQuad_r[13]+fUpwindQuad_r[12])+fUpwindQuad_r[11]+fUpwindQuad_r[10]+fUpwindQuad_r[9]+fUpwindQuad_r[8]+fUpwindQuad_r[7]+fUpwindQuad_r[6]+fUpwindQuad_r[5]+fUpwindQuad_r[4]-1.0*(fUpwindQuad_r[3]+fUpwindQuad_r[2]+fUpwindQuad_r[1]+fUpwindQuad_r[0])); 
  fUpwind_r[26] = 0.1767766952966368*(fUpwindQuad_r[31]-1.0*(fUpwindQuad_r[30]+fUpwindQuad_r[29])+fUpwindQuad_r[28]-1.0*fUpwindQuad_r[27]+fUpwindQuad_r[26]+fUpwindQuad_r[25]-1.0*(fUpwindQuad_r[24]+fUpwindQuad_r[23])+fUpwindQuad_r[22]+fUpwindQuad_r[21]-1.0*fUpwindQuad_r[20]+fUpwindQuad_r[19]-1.0*(fUpwindQuad_r[18]+fUpwindQuad_r[17])+fUpwindQuad_r[16]+fUpwindQuad_r[15]-1.0*(fUpwindQuad_r[14]+fUpwindQuad_r[13])+fUpwindQuad_r[12]-1.0*fUpwindQuad_r[11]+fUpwindQuad_r[10]+fUpwindQuad_r[9]-1.0*(fUpwindQuad_r[8]+fUpwindQuad_r[7])+fUpwindQuad_r[6]+fUpwindQuad_r[5]-1.0*fUpwindQuad_r[4]+fUpwindQuad_r[3]-1.0*(fUpwindQuad_r[2]+fUpwindQuad_r[1])+fUpwindQuad_r[0]); 
  fUpwind_r[27] = 0.1767766952966368*(fUpwindQuad_r[31]-1.0*(fUpwindQuad_r[30]+fUpwindQuad_r[29])+fUpwindQuad_r[28]-1.0*fUpwindQuad_r[27]+fUpwindQuad_r[26]+fUpwindQuad_r[25]-1.0*fUpwindQuad_r[24]+fUpwindQuad_r[23]-1.0*(fUpwindQuad_r[22]+fUpwindQuad_r[21])+fUpwindQuad_r[20]-1.0*fUpwindQuad_r[19]+fUpwindQuad_r[18]+fUpwindQuad_r[17]-1.0*(fUpwindQuad_r[16]+fUpwindQuad_r[15])+fUpwindQuad_r[14]+fUpwindQuad_r[13]-1.0*fUpwindQuad_r[12]+fUpwindQuad_r[11]-1.0*(fUpwindQuad_r[10]+fUpwindQuad_r[9])+fUpwindQuad_r[8]-1.0*fUpwindQuad_r[7]+fUpwindQuad_r[6]+fUpwindQuad_r[5]-1.0*fUpwindQuad_r[4]+fUpwindQuad_r[3]-1.0*(fUpwindQuad_r[2]+fUpwindQuad_r[1])+fUpwindQuad_r[0]); 
  fUpwind_r[28] = 0.1767766952966368*(fUpwindQuad_r[31]-1.0*(fUpwindQuad_r[30]+fUpwindQuad_r[29])+fUpwindQuad_r[28]+fUpwindQuad_r[27]-1.0*(fUpwindQuad_r[26]+fUpwindQuad_r[25])+fUpwindQuad_r[24]-1.0*fUpwindQuad_r[23]+fUpwindQuad_r[22]+fUpwindQuad_r[21]-1.0*(fUpwindQuad_r[20]+fUpwindQuad_r[19])+fUpwindQuad_r[18]+fUpwindQuad_r[17]-1.0*(fUpwindQuad_r[16]+fUpwindQuad_r[15])+fUpwindQuad_r[14]+fUpwindQuad_r[13]-1.0*(fUpwindQuad_r[12]+fUpwindQuad_r[11])+fUpwindQuad_r[10]+fUpwindQuad_r[9]-1.0*fUpwindQuad_r[8]+fUpwindQuad_r[7]-1.0*(fUpwindQuad_r[6]+fUpwindQuad_r[5])+fUpwindQuad_r[4]+fUpwindQuad_r[3]-1.0*(fUpwindQuad_r[2]+fUpwindQuad_r[1])+fUpwindQuad_r[0]); 
  fUpwind_r[29] = 0.1767766952966368*(fUpwindQuad_r[31]-1.0*fUpwindQuad_r[30]+fUpwindQuad_r[29]-1.0*(fUpwindQuad_r[28]+fUpwindQuad_r[27])+fUpwindQuad_r[26]-1.0*fUpwindQuad_r[25]+fUpwindQuad_r[24]-1.0*fUpwindQuad_r[23]+fUpwindQuad_r[22]-1.0*fUpwindQuad_r[21]+fUpwindQuad_r[20]+fUpwindQuad_r[19]-1.0*fUpwindQuad_r[18]+fUpwindQuad_r[17]-1.0*(fUpwindQuad_r[16]+fUpwindQuad_r[15])+fUpwindQuad_r[14]-1.0*fUpwindQuad_r[13]+fUpwindQuad_r[12]+fUpwindQuad_r[11]-1.0*fUpwindQuad_r[10]+fUpwindQuad_r[9]-1.0*fUpwindQuad_r[8]+fUpwindQuad_r[7]-1.0*fUpwindQuad_r[6]+fUpwindQuad_r[5]-1.0*(fUpwindQuad_r[4]+fUpwindQuad_r[3])+fUpwindQuad_r[2]-1.0*fUpwindQuad_r[1]+fUpwindQuad_r[0]); 
  fUpwind_r[30] = 0.1767766952966368*(fUpwindQuad_r[31]+fUpwindQuad_r[30]-1.0*(fUpwindQuad_r[29]+fUpwindQuad_r[28]+fUpwindQuad_r[27]+fUpwindQuad_r[26])+fUpwindQuad_r[25]+fUpwindQuad_r[24]-1.0*(fUpwindQuad_r[23]+fUpwindQuad_r[22])+fUpwindQuad_r[21]+fUpwindQuad_r[20]+fUpwindQuad_r[19]+fUpwindQuad_r[18]-1.0*(fUpwindQuad_r[17]+fUpwindQuad_r[16]+fUpwindQuad_r[15]+fUpwindQuad_r[14])+fUpwindQuad_r[13]+fUpwindQuad_r[12]+fUpwindQuad_r[11]+fUpwindQuad_r[10]-1.0*(fUpwindQuad_r[9]+fUpwindQuad_r[8])+fUpwindQuad_r[7]+fUpwindQuad_r[6]-1.0*(fUpwindQuad_r[5]+fUpwindQuad_r[4]+fUpwindQuad_r[3]+fUpwindQuad_r[2])+fUpwindQuad_r[1]+fUpwindQuad_r[0]); 
  fUpwind_r[31] = 0.1767766952966368*(fUpwindQuad_r[31]-1.0*(fUpwindQuad_r[30]+fUpwindQuad_r[29])+fUpwindQuad_r[28]-1.0*fUpwindQuad_r[27]+fUpwindQuad_r[26]+fUpwindQuad_r[25]-1.0*(fUpwindQuad_r[24]+fUpwindQuad_r[23])+fUpwindQuad_r[22]+fUpwindQuad_r[21]-1.0*fUpwindQuad_r[20]+fUpwindQuad_r[19]-1.0*(fUpwindQuad_r[18]+fUpwindQuad_r[17])+fUpwindQuad_r[16]-1.0*fUpwindQuad_r[15]+fUpwindQuad_r[14]+fUpwindQuad_r[13]-1.0*fUpwindQuad_r[12]+fUpwindQuad_r[11]-1.0*(fUpwindQuad_r[10]+fUpwindQuad_r[9])+fUpwindQuad_r[8]+fUpwindQuad_r[7]-1.0*(fUpwindQuad_r[6]+fUpwindQuad_r[5])+fUpwindQuad_r[4]-1.0*fUpwindQuad_r[3]+fUpwindQuad_r[2]+fUpwindQuad_r[1]-1.0*fUpwindQuad_r[0]); 

  Ghat_l[0] += 0.1767766952966368*(alpha[27]*fUpwind_l[27]+alpha[26]*fUpwind_l[26]+alpha[22]*fUpwind_l[22]+alpha[21]*fUpwind_l[21]+alpha[20]*fUpwind_l[20]+alpha[19]*fUpwind_l[19]+alpha[18]*fUpwind_l[18]+alpha[17]*fUpwind_l[17]+alpha[16]*fUpwind_l[16]+alpha[14]*fUpwind_l[14]+alpha[13]*fUpwind_l[13]+alpha[12]*fUpwind_l[12]+alpha[11]*fUpwind_l[11]+alpha[10]*fUpwind_l[10]+alpha[9]*fUpwind_l[9]+alpha[8]*fUpwind_l[8]+alpha[7]*fUpwind_l[7]+alpha[6]*fUpwind_l[6]+alpha[5]*fUpwind_l[5]+alpha[4]*fUpwind_l[4]+alpha[3]*fUpwind_l[3]+alpha[2]*fUpwind_l[2]+alpha[1]*fUpwind_l[1]+alpha[0]*fUpwind_l[0]); 
  Ghat_l[1] += 0.1767766952966368*(alpha[22]*fUpwind_l[27]+fUpwind_l[22]*alpha[27]+alpha[19]*fUpwind_l[26]+fUpwind_l[19]*alpha[26]+alpha[14]*fUpwind_l[21]+fUpwind_l[14]*alpha[21]+alpha[13]*fUpwind_l[20]+fUpwind_l[13]*alpha[20]+alpha[11]*fUpwind_l[18]+fUpwind_l[11]*alpha[18]+alpha[10]*fUpwind_l[17]+fUpwind_l[10]*alpha[17]+alpha[8]*fUpwind_l[16]+fUpwind_l[8]*alpha[16]+alpha[5]*fUpwind_l[12]+fUpwind_l[5]*alpha[12]+alpha[4]*fUpwind_l[9]+fUpwind_l[4]*alpha[9]+alpha[3]*fUpwind_l[7]+fUpwind_l[3]*alpha[7]+alpha[2]*fUpwind_l[6]+fUpwind_l[2]*alpha[6]+alpha[0]*fUpwind_l[1]+fUpwind_l[0]*alpha[1]); 
  Ghat_l[2] += 0.1767766952966368*(alpha[21]*fUpwind_l[27]+fUpwind_l[21]*alpha[27]+alpha[18]*fUpwind_l[26]+fUpwind_l[18]*alpha[26]+alpha[14]*fUpwind_l[22]+fUpwind_l[14]*alpha[22]+alpha[12]*fUpwind_l[20]+fUpwind_l[12]*alpha[20]+alpha[11]*fUpwind_l[19]+fUpwind_l[11]*alpha[19]+alpha[9]*fUpwind_l[17]+fUpwind_l[9]*alpha[17]+alpha[7]*fUpwind_l[16]+fUpwind_l[7]*alpha[16]+alpha[5]*fUpwind_l[13]+fUpwind_l[5]*alpha[13]+alpha[4]*fUpwind_l[10]+fUpwind_l[4]*alpha[10]+alpha[3]*fUpwind_l[8]+fUpwind_l[3]*alpha[8]+alpha[1]*fUpwind_l[6]+fUpwind_l[1]*alpha[6]+alpha[0]*fUpwind_l[2]+fUpwind_l[0]*alpha[2]); 
  Ghat_l[3] += 0.1767766952966368*(alpha[20]*fUpwind_l[27]+fUpwind_l[20]*alpha[27]+alpha[17]*fUpwind_l[26]+fUpwind_l[17]*alpha[26]+alpha[13]*fUpwind_l[22]+fUpwind_l[13]*alpha[22]+alpha[12]*fUpwind_l[21]+fUpwind_l[12]*alpha[21]+alpha[10]*fUpwind_l[19]+fUpwind_l[10]*alpha[19]+alpha[9]*fUpwind_l[18]+fUpwind_l[9]*alpha[18]+alpha[6]*fUpwind_l[16]+fUpwind_l[6]*alpha[16]+alpha[5]*fUpwind_l[14]+fUpwind_l[5]*alpha[14]+alpha[4]*fUpwind_l[11]+fUpwind_l[4]*alpha[11]+alpha[2]*fUpwind_l[8]+fUpwind_l[2]*alpha[8]+alpha[1]*fUpwind_l[7]+fUpwind_l[1]*alpha[7]+alpha[0]*fUpwind_l[3]+fUpwind_l[0]*alpha[3]); 
  Ghat_l[4] += 0.1767766952966368*(alpha[27]*fUpwind_l[31]+alpha[22]*fUpwind_l[30]+alpha[21]*fUpwind_l[29]+alpha[20]*fUpwind_l[28]+alpha[16]*fUpwind_l[26]+fUpwind_l[16]*alpha[26]+alpha[14]*fUpwind_l[25]+alpha[13]*fUpwind_l[24]+alpha[12]*fUpwind_l[23]+alpha[8]*fUpwind_l[19]+fUpwind_l[8]*alpha[19]+alpha[7]*fUpwind_l[18]+fUpwind_l[7]*alpha[18]+alpha[6]*fUpwind_l[17]+fUpwind_l[6]*alpha[17]+alpha[5]*fUpwind_l[15]+alpha[3]*fUpwind_l[11]+fUpwind_l[3]*alpha[11]+alpha[2]*fUpwind_l[10]+fUpwind_l[2]*alpha[10]+alpha[1]*fUpwind_l[9]+fUpwind_l[1]*alpha[9]+alpha[0]*fUpwind_l[4]+fUpwind_l[0]*alpha[4]); 
  Ghat_l[5] += 0.1767766952966368*(alpha[26]*fUpwind_l[31]+alpha[19]*fUpwind_l[30]+alpha[18]*fUpwind_l[29]+alpha[17]*fUpwind_l[28]+alpha[16]*fUpwind_l[27]+fUpwind_l[16]*alpha[27]+alpha[11]*fUpwind_l[25]+alpha[10]*fUpwind_l[24]+alpha[9]*fUpwind_l[23]+alpha[8]*fUpwind_l[22]+fUpwind_l[8]*alpha[22]+alpha[7]*fUpwind_l[21]+fUpwind_l[7]*alpha[21]+alpha[6]*fUpwind_l[20]+fUpwind_l[6]*alpha[20]+alpha[4]*fUpwind_l[15]+alpha[3]*fUpwind_l[14]+fUpwind_l[3]*alpha[14]+alpha[2]*fUpwind_l[13]+fUpwind_l[2]*alpha[13]+alpha[1]*fUpwind_l[12]+fUpwind_l[1]*alpha[12]+alpha[0]*fUpwind_l[5]+fUpwind_l[0]*alpha[5]); 
  Ghat_l[6] += 0.1767766952966368*(alpha[14]*fUpwind_l[27]+fUpwind_l[14]*alpha[27]+alpha[11]*fUpwind_l[26]+fUpwind_l[11]*alpha[26]+alpha[21]*fUpwind_l[22]+fUpwind_l[21]*alpha[22]+alpha[5]*fUpwind_l[20]+fUpwind_l[5]*alpha[20]+alpha[18]*fUpwind_l[19]+fUpwind_l[18]*alpha[19]+alpha[4]*fUpwind_l[17]+fUpwind_l[4]*alpha[17]+alpha[3]*fUpwind_l[16]+fUpwind_l[3]*alpha[16]+alpha[12]*fUpwind_l[13]+fUpwind_l[12]*alpha[13]+alpha[9]*fUpwind_l[10]+fUpwind_l[9]*alpha[10]+alpha[7]*fUpwind_l[8]+fUpwind_l[7]*alpha[8]+alpha[0]*fUpwind_l[6]+fUpwind_l[0]*alpha[6]+alpha[1]*fUpwind_l[2]+fUpwind_l[1]*alpha[2]); 
  Ghat_l[7] += 0.1767766952966368*(alpha[13]*fUpwind_l[27]+fUpwind_l[13]*alpha[27]+alpha[10]*fUpwind_l[26]+fUpwind_l[10]*alpha[26]+alpha[20]*fUpwind_l[22]+fUpwind_l[20]*alpha[22]+alpha[5]*fUpwind_l[21]+fUpwind_l[5]*alpha[21]+alpha[17]*fUpwind_l[19]+fUpwind_l[17]*alpha[19]+alpha[4]*fUpwind_l[18]+fUpwind_l[4]*alpha[18]+alpha[2]*fUpwind_l[16]+fUpwind_l[2]*alpha[16]+alpha[12]*fUpwind_l[14]+fUpwind_l[12]*alpha[14]+alpha[9]*fUpwind_l[11]+fUpwind_l[9]*alpha[11]+alpha[6]*fUpwind_l[8]+fUpwind_l[6]*alpha[8]+alpha[0]*fUpwind_l[7]+fUpwind_l[0]*alpha[7]+alpha[1]*fUpwind_l[3]+fUpwind_l[1]*alpha[3]); 
  Ghat_l[8] += 0.1767766952966368*(alpha[12]*fUpwind_l[27]+fUpwind_l[12]*alpha[27]+alpha[9]*fUpwind_l[26]+fUpwind_l[9]*alpha[26]+alpha[5]*fUpwind_l[22]+fUpwind_l[5]*alpha[22]+alpha[20]*fUpwind_l[21]+fUpwind_l[20]*alpha[21]+alpha[4]*fUpwind_l[19]+fUpwind_l[4]*alpha[19]+alpha[17]*fUpwind_l[18]+fUpwind_l[17]*alpha[18]+alpha[1]*fUpwind_l[16]+fUpwind_l[1]*alpha[16]+alpha[13]*fUpwind_l[14]+fUpwind_l[13]*alpha[14]+alpha[10]*fUpwind_l[11]+fUpwind_l[10]*alpha[11]+alpha[0]*fUpwind_l[8]+fUpwind_l[0]*alpha[8]+alpha[6]*fUpwind_l[7]+fUpwind_l[6]*alpha[7]+alpha[2]*fUpwind_l[3]+fUpwind_l[2]*alpha[3]); 
  Ghat_l[9] += 0.1767766952966368*(alpha[22]*fUpwind_l[31]+alpha[27]*fUpwind_l[30]+alpha[14]*fUpwind_l[29]+alpha[13]*fUpwind_l[28]+alpha[8]*fUpwind_l[26]+fUpwind_l[8]*alpha[26]+alpha[21]*fUpwind_l[25]+alpha[20]*fUpwind_l[24]+alpha[5]*fUpwind_l[23]+alpha[16]*fUpwind_l[19]+fUpwind_l[16]*alpha[19]+alpha[3]*fUpwind_l[18]+fUpwind_l[3]*alpha[18]+alpha[2]*fUpwind_l[17]+fUpwind_l[2]*alpha[17]+alpha[12]*fUpwind_l[15]+alpha[7]*fUpwind_l[11]+fUpwind_l[7]*alpha[11]+alpha[6]*fUpwind_l[10]+fUpwind_l[6]*alpha[10]+alpha[0]*fUpwind_l[9]+fUpwind_l[0]*alpha[9]+alpha[1]*fUpwind_l[4]+fUpwind_l[1]*alpha[4]); 
  Ghat_l[10] += 0.1767766952966368*(alpha[21]*fUpwind_l[31]+alpha[14]*fUpwind_l[30]+alpha[27]*fUpwind_l[29]+alpha[12]*fUpwind_l[28]+alpha[7]*fUpwind_l[26]+fUpwind_l[7]*alpha[26]+alpha[22]*fUpwind_l[25]+alpha[5]*fUpwind_l[24]+alpha[20]*fUpwind_l[23]+alpha[3]*fUpwind_l[19]+fUpwind_l[3]*alpha[19]+alpha[16]*fUpwind_l[18]+fUpwind_l[16]*alpha[18]+alpha[1]*fUpwind_l[17]+fUpwind_l[1]*alpha[17]+alpha[13]*fUpwind_l[15]+alpha[8]*fUpwind_l[11]+fUpwind_l[8]*alpha[11]+alpha[0]*fUpwind_l[10]+fUpwind_l[0]*alpha[10]+alpha[6]*fUpwind_l[9]+fUpwind_l[6]*alpha[9]+alpha[2]*fUpwind_l[4]+fUpwind_l[2]*alpha[4]); 
  Ghat_l[11] += 0.1767766952966368*(alpha[20]*fUpwind_l[31]+alpha[13]*fUpwind_l[30]+alpha[12]*fUpwind_l[29]+alpha[27]*fUpwind_l[28]+alpha[6]*fUpwind_l[26]+fUpwind_l[6]*alpha[26]+alpha[5]*fUpwind_l[25]+alpha[22]*fUpwind_l[24]+alpha[21]*fUpwind_l[23]+alpha[2]*fUpwind_l[19]+fUpwind_l[2]*alpha[19]+alpha[1]*fUpwind_l[18]+fUpwind_l[1]*alpha[18]+alpha[16]*fUpwind_l[17]+fUpwind_l[16]*alpha[17]+alpha[14]*fUpwind_l[15]+alpha[0]*fUpwind_l[11]+fUpwind_l[0]*alpha[11]+alpha[8]*fUpwind_l[10]+fUpwind_l[8]*alpha[10]+alpha[7]*fUpwind_l[9]+fUpwind_l[7]*alpha[9]+alpha[3]*fUpwind_l[4]+fUpwind_l[3]*alpha[4]); 
  Ghat_l[12] += 0.1767766952966368*(alpha[19]*fUpwind_l[31]+alpha[26]*fUpwind_l[30]+alpha[11]*fUpwind_l[29]+alpha[10]*fUpwind_l[28]+alpha[8]*fUpwind_l[27]+fUpwind_l[8]*alpha[27]+alpha[18]*fUpwind_l[25]+alpha[17]*fUpwind_l[24]+alpha[4]*fUpwind_l[23]+alpha[16]*fUpwind_l[22]+fUpwind_l[16]*alpha[22]+alpha[3]*fUpwind_l[21]+fUpwind_l[3]*alpha[21]+alpha[2]*fUpwind_l[20]+fUpwind_l[2]*alpha[20]+alpha[9]*fUpwind_l[15]+alpha[7]*fUpwind_l[14]+fUpwind_l[7]*alpha[14]+alpha[6]*fUpwind_l[13]+fUpwind_l[6]*alpha[13]+alpha[0]*fUpwind_l[12]+fUpwind_l[0]*alpha[12]+alpha[1]*fUpwind_l[5]+fUpwind_l[1]*alpha[5]); 
  Ghat_l[13] += 0.1767766952966368*(alpha[18]*fUpwind_l[31]+alpha[11]*fUpwind_l[30]+alpha[26]*fUpwind_l[29]+alpha[9]*fUpwind_l[28]+alpha[7]*fUpwind_l[27]+fUpwind_l[7]*alpha[27]+alpha[19]*fUpwind_l[25]+alpha[4]*fUpwind_l[24]+alpha[17]*fUpwind_l[23]+alpha[3]*fUpwind_l[22]+fUpwind_l[3]*alpha[22]+alpha[16]*fUpwind_l[21]+fUpwind_l[16]*alpha[21]+alpha[1]*fUpwind_l[20]+fUpwind_l[1]*alpha[20]+alpha[10]*fUpwind_l[15]+alpha[8]*fUpwind_l[14]+fUpwind_l[8]*alpha[14]+alpha[0]*fUpwind_l[13]+fUpwind_l[0]*alpha[13]+alpha[6]*fUpwind_l[12]+fUpwind_l[6]*alpha[12]+alpha[2]*fUpwind_l[5]+fUpwind_l[2]*alpha[5]); 
  Ghat_l[14] += 0.1767766952966368*(alpha[17]*fUpwind_l[31]+alpha[10]*fUpwind_l[30]+alpha[9]*fUpwind_l[29]+alpha[26]*fUpwind_l[28]+alpha[6]*fUpwind_l[27]+fUpwind_l[6]*alpha[27]+alpha[4]*fUpwind_l[25]+alpha[19]*fUpwind_l[24]+alpha[18]*fUpwind_l[23]+alpha[2]*fUpwind_l[22]+fUpwind_l[2]*alpha[22]+alpha[1]*fUpwind_l[21]+fUpwind_l[1]*alpha[21]+alpha[16]*fUpwind_l[20]+fUpwind_l[16]*alpha[20]+alpha[11]*fUpwind_l[15]+alpha[0]*fUpwind_l[14]+fUpwind_l[0]*alpha[14]+alpha[8]*fUpwind_l[13]+fUpwind_l[8]*alpha[13]+alpha[7]*fUpwind_l[12]+fUpwind_l[7]*alpha[12]+alpha[3]*fUpwind_l[5]+fUpwind_l[3]*alpha[5]); 
  Ghat_l[15] += 0.1767766952966368*(alpha[16]*fUpwind_l[31]+alpha[8]*fUpwind_l[30]+alpha[7]*fUpwind_l[29]+alpha[6]*fUpwind_l[28]+alpha[26]*fUpwind_l[27]+fUpwind_l[26]*alpha[27]+alpha[3]*fUpwind_l[25]+alpha[2]*fUpwind_l[24]+alpha[1]*fUpwind_l[23]+alpha[19]*fUpwind_l[22]+fUpwind_l[19]*alpha[22]+alpha[18]*fUpwind_l[21]+fUpwind_l[18]*alpha[21]+alpha[17]*fUpwind_l[20]+fUpwind_l[17]*alpha[20]+alpha[0]*fUpwind_l[15]+alpha[11]*fUpwind_l[14]+fUpwind_l[11]*alpha[14]+alpha[10]*fUpwind_l[13]+fUpwind_l[10]*alpha[13]+alpha[9]*fUpwind_l[12]+fUpwind_l[9]*alpha[12]+alpha[4]*fUpwind_l[5]+fUpwind_l[4]*alpha[5]); 
  Ghat_l[16] += 0.1767766952966368*(alpha[5]*fUpwind_l[27]+fUpwind_l[5]*alpha[27]+alpha[4]*fUpwind_l[26]+fUpwind_l[4]*alpha[26]+alpha[12]*fUpwind_l[22]+fUpwind_l[12]*alpha[22]+alpha[13]*fUpwind_l[21]+fUpwind_l[13]*alpha[21]+alpha[14]*fUpwind_l[20]+fUpwind_l[14]*alpha[20]+alpha[9]*fUpwind_l[19]+fUpwind_l[9]*alpha[19]+alpha[10]*fUpwind_l[18]+fUpwind_l[10]*alpha[18]+alpha[11]*fUpwind_l[17]+fUpwind_l[11]*alpha[17]+alpha[0]*fUpwind_l[16]+fUpwind_l[0]*alpha[16]+alpha[1]*fUpwind_l[8]+fUpwind_l[1]*alpha[8]+alpha[2]*fUpwind_l[7]+fUpwind_l[2]*alpha[7]+alpha[3]*fUpwind_l[6]+fUpwind_l[3]*alpha[6]); 
  Ghat_l[17] += 0.1767766952966368*(alpha[14]*fUpwind_l[31]+alpha[21]*fUpwind_l[30]+alpha[22]*fUpwind_l[29]+alpha[5]*fUpwind_l[28]+fUpwind_l[25]*alpha[27]+alpha[3]*fUpwind_l[26]+fUpwind_l[3]*alpha[26]+alpha[12]*fUpwind_l[24]+alpha[13]*fUpwind_l[23]+fUpwind_l[15]*alpha[20]+alpha[7]*fUpwind_l[19]+fUpwind_l[7]*alpha[19]+alpha[8]*fUpwind_l[18]+fUpwind_l[8]*alpha[18]+alpha[0]*fUpwind_l[17]+fUpwind_l[0]*alpha[17]+alpha[11]*fUpwind_l[16]+fUpwind_l[11]*alpha[16]+alpha[1]*fUpwind_l[10]+fUpwind_l[1]*alpha[10]+alpha[2]*fUpwind_l[9]+fUpwind_l[2]*alpha[9]+alpha[4]*fUpwind_l[6]+fUpwind_l[4]*alpha[6]); 
  Ghat_l[18] += 0.1767766952966368*(alpha[13]*fUpwind_l[31]+alpha[20]*fUpwind_l[30]+alpha[5]*fUpwind_l[29]+alpha[22]*fUpwind_l[28]+fUpwind_l[24]*alpha[27]+alpha[2]*fUpwind_l[26]+fUpwind_l[2]*alpha[26]+alpha[12]*fUpwind_l[25]+alpha[14]*fUpwind_l[23]+fUpwind_l[15]*alpha[21]+alpha[6]*fUpwind_l[19]+fUpwind_l[6]*alpha[19]+alpha[0]*fUpwind_l[18]+fUpwind_l[0]*alpha[18]+alpha[8]*fUpwind_l[17]+fUpwind_l[8]*alpha[17]+alpha[10]*fUpwind_l[16]+fUpwind_l[10]*alpha[16]+alpha[1]*fUpwind_l[11]+fUpwind_l[1]*alpha[11]+alpha[3]*fUpwind_l[9]+fUpwind_l[3]*alpha[9]+alpha[4]*fUpwind_l[7]+fUpwind_l[4]*alpha[7]); 
  Ghat_l[19] += 0.1767766952966368*(alpha[12]*fUpwind_l[31]+alpha[5]*fUpwind_l[30]+alpha[20]*fUpwind_l[29]+alpha[21]*fUpwind_l[28]+fUpwind_l[23]*alpha[27]+alpha[1]*fUpwind_l[26]+fUpwind_l[1]*alpha[26]+alpha[13]*fUpwind_l[25]+alpha[14]*fUpwind_l[24]+fUpwind_l[15]*alpha[22]+alpha[0]*fUpwind_l[19]+fUpwind_l[0]*alpha[19]+alpha[6]*fUpwind_l[18]+fUpwind_l[6]*alpha[18]+alpha[7]*fUpwind_l[17]+fUpwind_l[7]*alpha[17]+alpha[9]*fUpwind_l[16]+fUpwind_l[9]*alpha[16]+alpha[2]*fUpwind_l[11]+fUpwind_l[2]*alpha[11]+alpha[3]*fUpwind_l[10]+fUpwind_l[3]*alpha[10]+alpha[4]*fUpwind_l[8]+fUpwind_l[4]*alpha[8]); 
  Ghat_l[20] += 0.1767766952966368*(alpha[11]*fUpwind_l[31]+alpha[18]*fUpwind_l[30]+alpha[19]*fUpwind_l[29]+alpha[4]*fUpwind_l[28]+alpha[3]*fUpwind_l[27]+fUpwind_l[3]*alpha[27]+fUpwind_l[25]*alpha[26]+alpha[9]*fUpwind_l[24]+alpha[10]*fUpwind_l[23]+alpha[7]*fUpwind_l[22]+fUpwind_l[7]*alpha[22]+alpha[8]*fUpwind_l[21]+fUpwind_l[8]*alpha[21]+alpha[0]*fUpwind_l[20]+fUpwind_l[0]*alpha[20]+fUpwind_l[15]*alpha[17]+alpha[14]*fUpwind_l[16]+fUpwind_l[14]*alpha[16]+alpha[1]*fUpwind_l[13]+fUpwind_l[1]*alpha[13]+alpha[2]*fUpwind_l[12]+fUpwind_l[2]*alpha[12]+alpha[5]*fUpwind_l[6]+fUpwind_l[5]*alpha[6]); 
  Ghat_l[21] += 0.1767766952966368*(alpha[10]*fUpwind_l[31]+alpha[17]*fUpwind_l[30]+alpha[4]*fUpwind_l[29]+alpha[19]*fUpwind_l[28]+alpha[2]*fUpwind_l[27]+fUpwind_l[2]*alpha[27]+fUpwind_l[24]*alpha[26]+alpha[9]*fUpwind_l[25]+alpha[11]*fUpwind_l[23]+alpha[6]*fUpwind_l[22]+fUpwind_l[6]*alpha[22]+alpha[0]*fUpwind_l[21]+fUpwind_l[0]*alpha[21]+alpha[8]*fUpwind_l[20]+fUpwind_l[8]*alpha[20]+fUpwind_l[15]*alpha[18]+alpha[13]*fUpwind_l[16]+fUpwind_l[13]*alpha[16]+alpha[1]*fUpwind_l[14]+fUpwind_l[1]*alpha[14]+alpha[3]*fUpwind_l[12]+fUpwind_l[3]*alpha[12]+alpha[5]*fUpwind_l[7]+fUpwind_l[5]*alpha[7]); 
  Ghat_l[22] += 0.1767766952966368*(alpha[9]*fUpwind_l[31]+alpha[4]*fUpwind_l[30]+alpha[17]*fUpwind_l[29]+alpha[18]*fUpwind_l[28]+alpha[1]*fUpwind_l[27]+fUpwind_l[1]*alpha[27]+fUpwind_l[23]*alpha[26]+alpha[10]*fUpwind_l[25]+alpha[11]*fUpwind_l[24]+alpha[0]*fUpwind_l[22]+fUpwind_l[0]*alpha[22]+alpha[6]*fUpwind_l[21]+fUpwind_l[6]*alpha[21]+alpha[7]*fUpwind_l[20]+fUpwind_l[7]*alpha[20]+fUpwind_l[15]*alpha[19]+alpha[12]*fUpwind_l[16]+fUpwind_l[12]*alpha[16]+alpha[2]*fUpwind_l[14]+fUpwind_l[2]*alpha[14]+alpha[3]*fUpwind_l[13]+fUpwind_l[3]*alpha[13]+alpha[5]*fUpwind_l[8]+fUpwind_l[5]*alpha[8]); 
  Ghat_l[23] += 0.1767766952966368*(alpha[8]*fUpwind_l[31]+alpha[16]*fUpwind_l[30]+alpha[3]*fUpwind_l[29]+alpha[2]*fUpwind_l[28]+alpha[19]*fUpwind_l[27]+fUpwind_l[19]*alpha[27]+alpha[22]*fUpwind_l[26]+fUpwind_l[22]*alpha[26]+alpha[7]*fUpwind_l[25]+alpha[6]*fUpwind_l[24]+alpha[0]*fUpwind_l[23]+alpha[11]*fUpwind_l[21]+fUpwind_l[11]*alpha[21]+alpha[10]*fUpwind_l[20]+fUpwind_l[10]*alpha[20]+alpha[14]*fUpwind_l[18]+fUpwind_l[14]*alpha[18]+alpha[13]*fUpwind_l[17]+fUpwind_l[13]*alpha[17]+alpha[1]*fUpwind_l[15]+alpha[4]*fUpwind_l[12]+fUpwind_l[4]*alpha[12]+alpha[5]*fUpwind_l[9]+fUpwind_l[5]*alpha[9]); 
  Ghat_l[24] += 0.1767766952966368*(alpha[7]*fUpwind_l[31]+alpha[3]*fUpwind_l[30]+alpha[16]*fUpwind_l[29]+alpha[1]*fUpwind_l[28]+alpha[18]*fUpwind_l[27]+fUpwind_l[18]*alpha[27]+alpha[21]*fUpwind_l[26]+fUpwind_l[21]*alpha[26]+alpha[8]*fUpwind_l[25]+alpha[0]*fUpwind_l[24]+alpha[6]*fUpwind_l[23]+alpha[11]*fUpwind_l[22]+fUpwind_l[11]*alpha[22]+alpha[9]*fUpwind_l[20]+fUpwind_l[9]*alpha[20]+alpha[14]*fUpwind_l[19]+fUpwind_l[14]*alpha[19]+alpha[12]*fUpwind_l[17]+fUpwind_l[12]*alpha[17]+alpha[2]*fUpwind_l[15]+alpha[4]*fUpwind_l[13]+fUpwind_l[4]*alpha[13]+alpha[5]*fUpwind_l[10]+fUpwind_l[5]*alpha[10]); 
  Ghat_l[25] += 0.1767766952966368*(alpha[6]*fUpwind_l[31]+alpha[2]*fUpwind_l[30]+alpha[1]*fUpwind_l[29]+alpha[16]*fUpwind_l[28]+alpha[17]*fUpwind_l[27]+fUpwind_l[17]*alpha[27]+alpha[20]*fUpwind_l[26]+fUpwind_l[20]*alpha[26]+alpha[0]*fUpwind_l[25]+alpha[8]*fUpwind_l[24]+alpha[7]*fUpwind_l[23]+alpha[10]*fUpwind_l[22]+fUpwind_l[10]*alpha[22]+alpha[9]*fUpwind_l[21]+fUpwind_l[9]*alpha[21]+alpha[13]*fUpwind_l[19]+fUpwind_l[13]*alpha[19]+alpha[12]*fUpwind_l[18]+fUpwind_l[12]*alpha[18]+alpha[3]*fUpwind_l[15]+alpha[4]*fUpwind_l[14]+fUpwind_l[4]*alpha[14]+alpha[5]*fUpwind_l[11]+fUpwind_l[5]*alpha[11]); 
  Ghat_l[26] += 0.1767766952966368*(alpha[5]*fUpwind_l[31]+alpha[12]*fUpwind_l[30]+alpha[13]*fUpwind_l[29]+alpha[14]*fUpwind_l[28]+fUpwind_l[15]*alpha[27]+alpha[0]*fUpwind_l[26]+fUpwind_l[0]*alpha[26]+alpha[20]*fUpwind_l[25]+alpha[21]*fUpwind_l[24]+alpha[22]*fUpwind_l[23]+alpha[1]*fUpwind_l[19]+fUpwind_l[1]*alpha[19]+alpha[2]*fUpwind_l[18]+fUpwind_l[2]*alpha[18]+alpha[3]*fUpwind_l[17]+fUpwind_l[3]*alpha[17]+alpha[4]*fUpwind_l[16]+fUpwind_l[4]*alpha[16]+alpha[6]*fUpwind_l[11]+fUpwind_l[6]*alpha[11]+alpha[7]*fUpwind_l[10]+fUpwind_l[7]*alpha[10]+alpha[8]*fUpwind_l[9]+fUpwind_l[8]*alpha[9]); 
  Ghat_l[27] += 0.1767766952966368*(alpha[4]*fUpwind_l[31]+alpha[9]*fUpwind_l[30]+alpha[10]*fUpwind_l[29]+alpha[11]*fUpwind_l[28]+alpha[0]*fUpwind_l[27]+fUpwind_l[0]*alpha[27]+fUpwind_l[15]*alpha[26]+alpha[17]*fUpwind_l[25]+alpha[18]*fUpwind_l[24]+alpha[19]*fUpwind_l[23]+alpha[1]*fUpwind_l[22]+fUpwind_l[1]*alpha[22]+alpha[2]*fUpwind_l[21]+fUpwind_l[2]*alpha[21]+alpha[3]*fUpwind_l[20]+fUpwind_l[3]*alpha[20]+alpha[5]*fUpwind_l[16]+fUpwind_l[5]*alpha[16]+alpha[6]*fUpwind_l[14]+fUpwind_l[6]*alpha[14]+alpha[7]*fUpwind_l[13]+fUpwind_l[7]*alpha[13]+alpha[8]*fUpwind_l[12]+fUpwind_l[8]*alpha[12]); 
  Ghat_l[28] += 0.1767766952966368*(alpha[3]*fUpwind_l[31]+alpha[7]*fUpwind_l[30]+alpha[8]*fUpwind_l[29]+alpha[0]*fUpwind_l[28]+alpha[11]*fUpwind_l[27]+fUpwind_l[11]*alpha[27]+alpha[14]*fUpwind_l[26]+fUpwind_l[14]*alpha[26]+alpha[16]*fUpwind_l[25]+alpha[1]*fUpwind_l[24]+alpha[2]*fUpwind_l[23]+alpha[18]*fUpwind_l[22]+fUpwind_l[18]*alpha[22]+alpha[19]*fUpwind_l[21]+fUpwind_l[19]*alpha[21]+alpha[4]*fUpwind_l[20]+fUpwind_l[4]*alpha[20]+alpha[5]*fUpwind_l[17]+fUpwind_l[5]*alpha[17]+alpha[6]*fUpwind_l[15]+alpha[9]*fUpwind_l[13]+fUpwind_l[9]*alpha[13]+alpha[10]*fUpwind_l[12]+fUpwind_l[10]*alpha[12]); 
  Ghat_l[29] += 0.1767766952966368*(alpha[2]*fUpwind_l[31]+alpha[6]*fUpwind_l[30]+alpha[0]*fUpwind_l[29]+alpha[8]*fUpwind_l[28]+alpha[10]*fUpwind_l[27]+fUpwind_l[10]*alpha[27]+alpha[13]*fUpwind_l[26]+fUpwind_l[13]*alpha[26]+alpha[1]*fUpwind_l[25]+alpha[16]*fUpwind_l[24]+alpha[3]*fUpwind_l[23]+alpha[17]*fUpwind_l[22]+fUpwind_l[17]*alpha[22]+alpha[4]*fUpwind_l[21]+fUpwind_l[4]*alpha[21]+alpha[19]*fUpwind_l[20]+fUpwind_l[19]*alpha[20]+alpha[5]*fUpwind_l[18]+fUpwind_l[5]*alpha[18]+alpha[7]*fUpwind_l[15]+alpha[9]*fUpwind_l[14]+fUpwind_l[9]*alpha[14]+alpha[11]*fUpwind_l[12]+fUpwind_l[11]*alpha[12]); 
  Ghat_l[30] += 0.1767766952966368*(alpha[1]*fUpwind_l[31]+alpha[0]*fUpwind_l[30]+alpha[6]*fUpwind_l[29]+alpha[7]*fUpwind_l[28]+alpha[9]*fUpwind_l[27]+fUpwind_l[9]*alpha[27]+alpha[12]*fUpwind_l[26]+fUpwind_l[12]*alpha[26]+alpha[2]*fUpwind_l[25]+alpha[3]*fUpwind_l[24]+alpha[16]*fUpwind_l[23]+alpha[4]*fUpwind_l[22]+fUpwind_l[4]*alpha[22]+alpha[17]*fUpwind_l[21]+fUpwind_l[17]*alpha[21]+alpha[18]*fUpwind_l[20]+fUpwind_l[18]*alpha[20]+alpha[5]*fUpwind_l[19]+fUpwind_l[5]*alpha[19]+alpha[8]*fUpwind_l[15]+alpha[10]*fUpwind_l[14]+fUpwind_l[10]*alpha[14]+alpha[11]*fUpwind_l[13]+fUpwind_l[11]*alpha[13]); 
  Ghat_l[31] += 0.1767766952966368*(alpha[0]*fUpwind_l[31]+alpha[1]*fUpwind_l[30]+alpha[2]*fUpwind_l[29]+alpha[3]*fUpwind_l[28]+alpha[4]*fUpwind_l[27]+fUpwind_l[4]*alpha[27]+alpha[5]*fUpwind_l[26]+fUpwind_l[5]*alpha[26]+alpha[6]*fUpwind_l[25]+alpha[7]*fUpwind_l[24]+alpha[8]*fUpwind_l[23]+alpha[9]*fUpwind_l[22]+fUpwind_l[9]*alpha[22]+alpha[10]*fUpwind_l[21]+fUpwind_l[10]*alpha[21]+alpha[11]*fUpwind_l[20]+fUpwind_l[11]*alpha[20]+alpha[12]*fUpwind_l[19]+fUpwind_l[12]*alpha[19]+alpha[13]*fUpwind_l[18]+fUpwind_l[13]*alpha[18]+alpha[14]*fUpwind_l[17]+fUpwind_l[14]*alpha[17]+fUpwind_l[15]*alpha[16]); 

  Ghat_r[0] += 0.1767766952966368*(alpha[27]*fUpwind_r[27]+alpha[26]*fUpwind_r[26]+alpha[22]*fUpwind_r[22]+alpha[21]*fUpwind_r[21]+alpha[20]*fUpwind_r[20]+alpha[19]*fUpwind_r[19]+alpha[18]*fUpwind_r[18]+alpha[17]*fUpwind_r[17]+alpha[16]*fUpwind_r[16]+alpha[14]*fUpwind_r[14]+alpha[13]*fUpwind_r[13]+alpha[12]*fUpwind_r[12]+alpha[11]*fUpwind_r[11]+alpha[10]*fUpwind_r[10]+alpha[9]*fUpwind_r[9]+alpha[8]*fUpwind_r[8]+alpha[7]*fUpwind_r[7]+alpha[6]*fUpwind_r[6]+alpha[5]*fUpwind_r[5]+alpha[4]*fUpwind_r[4]+alpha[3]*fUpwind_r[3]+alpha[2]*fUpwind_r[2]+alpha[1]*fUpwind_r[1]+alpha[0]*fUpwind_r[0]); 
  Ghat_r[1] += 0.1767766952966368*(alpha[22]*fUpwind_r[27]+fUpwind_r[22]*alpha[27]+alpha[19]*fUpwind_r[26]+fUpwind_r[19]*alpha[26]+alpha[14]*fUpwind_r[21]+fUpwind_r[14]*alpha[21]+alpha[13]*fUpwind_r[20]+fUpwind_r[13]*alpha[20]+alpha[11]*fUpwind_r[18]+fUpwind_r[11]*alpha[18]+alpha[10]*fUpwind_r[17]+fUpwind_r[10]*alpha[17]+alpha[8]*fUpwind_r[16]+fUpwind_r[8]*alpha[16]+alpha[5]*fUpwind_r[12]+fUpwind_r[5]*alpha[12]+alpha[4]*fUpwind_r[9]+fUpwind_r[4]*alpha[9]+alpha[3]*fUpwind_r[7]+fUpwind_r[3]*alpha[7]+alpha[2]*fUpwind_r[6]+fUpwind_r[2]*alpha[6]+alpha[0]*fUpwind_r[1]+fUpwind_r[0]*alpha[1]); 
  Ghat_r[2] += 0.1767766952966368*(alpha[21]*fUpwind_r[27]+fUpwind_r[21]*alpha[27]+alpha[18]*fUpwind_r[26]+fUpwind_r[18]*alpha[26]+alpha[14]*fUpwind_r[22]+fUpwind_r[14]*alpha[22]+alpha[12]*fUpwind_r[20]+fUpwind_r[12]*alpha[20]+alpha[11]*fUpwind_r[19]+fUpwind_r[11]*alpha[19]+alpha[9]*fUpwind_r[17]+fUpwind_r[9]*alpha[17]+alpha[7]*fUpwind_r[16]+fUpwind_r[7]*alpha[16]+alpha[5]*fUpwind_r[13]+fUpwind_r[5]*alpha[13]+alpha[4]*fUpwind_r[10]+fUpwind_r[4]*alpha[10]+alpha[3]*fUpwind_r[8]+fUpwind_r[3]*alpha[8]+alpha[1]*fUpwind_r[6]+fUpwind_r[1]*alpha[6]+alpha[0]*fUpwind_r[2]+fUpwind_r[0]*alpha[2]); 
  Ghat_r[3] += 0.1767766952966368*(alpha[20]*fUpwind_r[27]+fUpwind_r[20]*alpha[27]+alpha[17]*fUpwind_r[26]+fUpwind_r[17]*alpha[26]+alpha[13]*fUpwind_r[22]+fUpwind_r[13]*alpha[22]+alpha[12]*fUpwind_r[21]+fUpwind_r[12]*alpha[21]+alpha[10]*fUpwind_r[19]+fUpwind_r[10]*alpha[19]+alpha[9]*fUpwind_r[18]+fUpwind_r[9]*alpha[18]+alpha[6]*fUpwind_r[16]+fUpwind_r[6]*alpha[16]+alpha[5]*fUpwind_r[14]+fUpwind_r[5]*alpha[14]+alpha[4]*fUpwind_r[11]+fUpwind_r[4]*alpha[11]+alpha[2]*fUpwind_r[8]+fUpwind_r[2]*alpha[8]+alpha[1]*fUpwind_r[7]+fUpwind_r[1]*alpha[7]+alpha[0]*fUpwind_r[3]+fUpwind_r[0]*alpha[3]); 
  Ghat_r[4] += 0.1767766952966368*(alpha[27]*fUpwind_r[31]+alpha[22]*fUpwind_r[30]+alpha[21]*fUpwind_r[29]+alpha[20]*fUpwind_r[28]+alpha[16]*fUpwind_r[26]+fUpwind_r[16]*alpha[26]+alpha[14]*fUpwind_r[25]+alpha[13]*fUpwind_r[24]+alpha[12]*fUpwind_r[23]+alpha[8]*fUpwind_r[19]+fUpwind_r[8]*alpha[19]+alpha[7]*fUpwind_r[18]+fUpwind_r[7]*alpha[18]+alpha[6]*fUpwind_r[17]+fUpwind_r[6]*alpha[17]+alpha[5]*fUpwind_r[15]+alpha[3]*fUpwind_r[11]+fUpwind_r[3]*alpha[11]+alpha[2]*fUpwind_r[10]+fUpwind_r[2]*alpha[10]+alpha[1]*fUpwind_r[9]+fUpwind_r[1]*alpha[9]+alpha[0]*fUpwind_r[4]+fUpwind_r[0]*alpha[4]); 
  Ghat_r[5] += 0.1767766952966368*(alpha[26]*fUpwind_r[31]+alpha[19]*fUpwind_r[30]+alpha[18]*fUpwind_r[29]+alpha[17]*fUpwind_r[28]+alpha[16]*fUpwind_r[27]+fUpwind_r[16]*alpha[27]+alpha[11]*fUpwind_r[25]+alpha[10]*fUpwind_r[24]+alpha[9]*fUpwind_r[23]+alpha[8]*fUpwind_r[22]+fUpwind_r[8]*alpha[22]+alpha[7]*fUpwind_r[21]+fUpwind_r[7]*alpha[21]+alpha[6]*fUpwind_r[20]+fUpwind_r[6]*alpha[20]+alpha[4]*fUpwind_r[15]+alpha[3]*fUpwind_r[14]+fUpwind_r[3]*alpha[14]+alpha[2]*fUpwind_r[13]+fUpwind_r[2]*alpha[13]+alpha[1]*fUpwind_r[12]+fUpwind_r[1]*alpha[12]+alpha[0]*fUpwind_r[5]+fUpwind_r[0]*alpha[5]); 
  Ghat_r[6] += 0.1767766952966368*(alpha[14]*fUpwind_r[27]+fUpwind_r[14]*alpha[27]+alpha[11]*fUpwind_r[26]+fUpwind_r[11]*alpha[26]+alpha[21]*fUpwind_r[22]+fUpwind_r[21]*alpha[22]+alpha[5]*fUpwind_r[20]+fUpwind_r[5]*alpha[20]+alpha[18]*fUpwind_r[19]+fUpwind_r[18]*alpha[19]+alpha[4]*fUpwind_r[17]+fUpwind_r[4]*alpha[17]+alpha[3]*fUpwind_r[16]+fUpwind_r[3]*alpha[16]+alpha[12]*fUpwind_r[13]+fUpwind_r[12]*alpha[13]+alpha[9]*fUpwind_r[10]+fUpwind_r[9]*alpha[10]+alpha[7]*fUpwind_r[8]+fUpwind_r[7]*alpha[8]+alpha[0]*fUpwind_r[6]+fUpwind_r[0]*alpha[6]+alpha[1]*fUpwind_r[2]+fUpwind_r[1]*alpha[2]); 
  Ghat_r[7] += 0.1767766952966368*(alpha[13]*fUpwind_r[27]+fUpwind_r[13]*alpha[27]+alpha[10]*fUpwind_r[26]+fUpwind_r[10]*alpha[26]+alpha[20]*fUpwind_r[22]+fUpwind_r[20]*alpha[22]+alpha[5]*fUpwind_r[21]+fUpwind_r[5]*alpha[21]+alpha[17]*fUpwind_r[19]+fUpwind_r[17]*alpha[19]+alpha[4]*fUpwind_r[18]+fUpwind_r[4]*alpha[18]+alpha[2]*fUpwind_r[16]+fUpwind_r[2]*alpha[16]+alpha[12]*fUpwind_r[14]+fUpwind_r[12]*alpha[14]+alpha[9]*fUpwind_r[11]+fUpwind_r[9]*alpha[11]+alpha[6]*fUpwind_r[8]+fUpwind_r[6]*alpha[8]+alpha[0]*fUpwind_r[7]+fUpwind_r[0]*alpha[7]+alpha[1]*fUpwind_r[3]+fUpwind_r[1]*alpha[3]); 
  Ghat_r[8] += 0.1767766952966368*(alpha[12]*fUpwind_r[27]+fUpwind_r[12]*alpha[27]+alpha[9]*fUpwind_r[26]+fUpwind_r[9]*alpha[26]+alpha[5]*fUpwind_r[22]+fUpwind_r[5]*alpha[22]+alpha[20]*fUpwind_r[21]+fUpwind_r[20]*alpha[21]+alpha[4]*fUpwind_r[19]+fUpwind_r[4]*alpha[19]+alpha[17]*fUpwind_r[18]+fUpwind_r[17]*alpha[18]+alpha[1]*fUpwind_r[16]+fUpwind_r[1]*alpha[16]+alpha[13]*fUpwind_r[14]+fUpwind_r[13]*alpha[14]+alpha[10]*fUpwind_r[11]+fUpwind_r[10]*alpha[11]+alpha[0]*fUpwind_r[8]+fUpwind_r[0]*alpha[8]+alpha[6]*fUpwind_r[7]+fUpwind_r[6]*alpha[7]+alpha[2]*fUpwind_r[3]+fUpwind_r[2]*alpha[3]); 
  Ghat_r[9] += 0.1767766952966368*(alpha[22]*fUpwind_r[31]+alpha[27]*fUpwind_r[30]+alpha[14]*fUpwind_r[29]+alpha[13]*fUpwind_r[28]+alpha[8]*fUpwind_r[26]+fUpwind_r[8]*alpha[26]+alpha[21]*fUpwind_r[25]+alpha[20]*fUpwind_r[24]+alpha[5]*fUpwind_r[23]+alpha[16]*fUpwind_r[19]+fUpwind_r[16]*alpha[19]+alpha[3]*fUpwind_r[18]+fUpwind_r[3]*alpha[18]+alpha[2]*fUpwind_r[17]+fUpwind_r[2]*alpha[17]+alpha[12]*fUpwind_r[15]+alpha[7]*fUpwind_r[11]+fUpwind_r[7]*alpha[11]+alpha[6]*fUpwind_r[10]+fUpwind_r[6]*alpha[10]+alpha[0]*fUpwind_r[9]+fUpwind_r[0]*alpha[9]+alpha[1]*fUpwind_r[4]+fUpwind_r[1]*alpha[4]); 
  Ghat_r[10] += 0.1767766952966368*(alpha[21]*fUpwind_r[31]+alpha[14]*fUpwind_r[30]+alpha[27]*fUpwind_r[29]+alpha[12]*fUpwind_r[28]+alpha[7]*fUpwind_r[26]+fUpwind_r[7]*alpha[26]+alpha[22]*fUpwind_r[25]+alpha[5]*fUpwind_r[24]+alpha[20]*fUpwind_r[23]+alpha[3]*fUpwind_r[19]+fUpwind_r[3]*alpha[19]+alpha[16]*fUpwind_r[18]+fUpwind_r[16]*alpha[18]+alpha[1]*fUpwind_r[17]+fUpwind_r[1]*alpha[17]+alpha[13]*fUpwind_r[15]+alpha[8]*fUpwind_r[11]+fUpwind_r[8]*alpha[11]+alpha[0]*fUpwind_r[10]+fUpwind_r[0]*alpha[10]+alpha[6]*fUpwind_r[9]+fUpwind_r[6]*alpha[9]+alpha[2]*fUpwind_r[4]+fUpwind_r[2]*alpha[4]); 
  Ghat_r[11] += 0.1767766952966368*(alpha[20]*fUpwind_r[31]+alpha[13]*fUpwind_r[30]+alpha[12]*fUpwind_r[29]+alpha[27]*fUpwind_r[28]+alpha[6]*fUpwind_r[26]+fUpwind_r[6]*alpha[26]+alpha[5]*fUpwind_r[25]+alpha[22]*fUpwind_r[24]+alpha[21]*fUpwind_r[23]+alpha[2]*fUpwind_r[19]+fUpwind_r[2]*alpha[19]+alpha[1]*fUpwind_r[18]+fUpwind_r[1]*alpha[18]+alpha[16]*fUpwind_r[17]+fUpwind_r[16]*alpha[17]+alpha[14]*fUpwind_r[15]+alpha[0]*fUpwind_r[11]+fUpwind_r[0]*alpha[11]+alpha[8]*fUpwind_r[10]+fUpwind_r[8]*alpha[10]+alpha[7]*fUpwind_r[9]+fUpwind_r[7]*alpha[9]+alpha[3]*fUpwind_r[4]+fUpwind_r[3]*alpha[4]); 
  Ghat_r[12] += 0.1767766952966368*(alpha[19]*fUpwind_r[31]+alpha[26]*fUpwind_r[30]+alpha[11]*fUpwind_r[29]+alpha[10]*fUpwind_r[28]+alpha[8]*fUpwind_r[27]+fUpwind_r[8]*alpha[27]+alpha[18]*fUpwind_r[25]+alpha[17]*fUpwind_r[24]+alpha[4]*fUpwind_r[23]+alpha[16]*fUpwind_r[22]+fUpwind_r[16]*alpha[22]+alpha[3]*fUpwind_r[21]+fUpwind_r[3]*alpha[21]+alpha[2]*fUpwind_r[20]+fUpwind_r[2]*alpha[20]+alpha[9]*fUpwind_r[15]+alpha[7]*fUpwind_r[14]+fUpwind_r[7]*alpha[14]+alpha[6]*fUpwind_r[13]+fUpwind_r[6]*alpha[13]+alpha[0]*fUpwind_r[12]+fUpwind_r[0]*alpha[12]+alpha[1]*fUpwind_r[5]+fUpwind_r[1]*alpha[5]); 
  Ghat_r[13] += 0.1767766952966368*(alpha[18]*fUpwind_r[31]+alpha[11]*fUpwind_r[30]+alpha[26]*fUpwind_r[29]+alpha[9]*fUpwind_r[28]+alpha[7]*fUpwind_r[27]+fUpwind_r[7]*alpha[27]+alpha[19]*fUpwind_r[25]+alpha[4]*fUpwind_r[24]+alpha[17]*fUpwind_r[23]+alpha[3]*fUpwind_r[22]+fUpwind_r[3]*alpha[22]+alpha[16]*fUpwind_r[21]+fUpwind_r[16]*alpha[21]+alpha[1]*fUpwind_r[20]+fUpwind_r[1]*alpha[20]+alpha[10]*fUpwind_r[15]+alpha[8]*fUpwind_r[14]+fUpwind_r[8]*alpha[14]+alpha[0]*fUpwind_r[13]+fUpwind_r[0]*alpha[13]+alpha[6]*fUpwind_r[12]+fUpwind_r[6]*alpha[12]+alpha[2]*fUpwind_r[5]+fUpwind_r[2]*alpha[5]); 
  Ghat_r[14] += 0.1767766952966368*(alpha[17]*fUpwind_r[31]+alpha[10]*fUpwind_r[30]+alpha[9]*fUpwind_r[29]+alpha[26]*fUpwind_r[28]+alpha[6]*fUpwind_r[27]+fUpwind_r[6]*alpha[27]+alpha[4]*fUpwind_r[25]+alpha[19]*fUpwind_r[24]+alpha[18]*fUpwind_r[23]+alpha[2]*fUpwind_r[22]+fUpwind_r[2]*alpha[22]+alpha[1]*fUpwind_r[21]+fUpwind_r[1]*alpha[21]+alpha[16]*fUpwind_r[20]+fUpwind_r[16]*alpha[20]+alpha[11]*fUpwind_r[15]+alpha[0]*fUpwind_r[14]+fUpwind_r[0]*alpha[14]+alpha[8]*fUpwind_r[13]+fUpwind_r[8]*alpha[13]+alpha[7]*fUpwind_r[12]+fUpwind_r[7]*alpha[12]+alpha[3]*fUpwind_r[5]+fUpwind_r[3]*alpha[5]); 
  Ghat_r[15] += 0.1767766952966368*(alpha[16]*fUpwind_r[31]+alpha[8]*fUpwind_r[30]+alpha[7]*fUpwind_r[29]+alpha[6]*fUpwind_r[28]+alpha[26]*fUpwind_r[27]+fUpwind_r[26]*alpha[27]+alpha[3]*fUpwind_r[25]+alpha[2]*fUpwind_r[24]+alpha[1]*fUpwind_r[23]+alpha[19]*fUpwind_r[22]+fUpwind_r[19]*alpha[22]+alpha[18]*fUpwind_r[21]+fUpwind_r[18]*alpha[21]+alpha[17]*fUpwind_r[20]+fUpwind_r[17]*alpha[20]+alpha[0]*fUpwind_r[15]+alpha[11]*fUpwind_r[14]+fUpwind_r[11]*alpha[14]+alpha[10]*fUpwind_r[13]+fUpwind_r[10]*alpha[13]+alpha[9]*fUpwind_r[12]+fUpwind_r[9]*alpha[12]+alpha[4]*fUpwind_r[5]+fUpwind_r[4]*alpha[5]); 
  Ghat_r[16] += 0.1767766952966368*(alpha[5]*fUpwind_r[27]+fUpwind_r[5]*alpha[27]+alpha[4]*fUpwind_r[26]+fUpwind_r[4]*alpha[26]+alpha[12]*fUpwind_r[22]+fUpwind_r[12]*alpha[22]+alpha[13]*fUpwind_r[21]+fUpwind_r[13]*alpha[21]+alpha[14]*fUpwind_r[20]+fUpwind_r[14]*alpha[20]+alpha[9]*fUpwind_r[19]+fUpwind_r[9]*alpha[19]+alpha[10]*fUpwind_r[18]+fUpwind_r[10]*alpha[18]+alpha[11]*fUpwind_r[17]+fUpwind_r[11]*alpha[17]+alpha[0]*fUpwind_r[16]+fUpwind_r[0]*alpha[16]+alpha[1]*fUpwind_r[8]+fUpwind_r[1]*alpha[8]+alpha[2]*fUpwind_r[7]+fUpwind_r[2]*alpha[7]+alpha[3]*fUpwind_r[6]+fUpwind_r[3]*alpha[6]); 
  Ghat_r[17] += 0.1767766952966368*(alpha[14]*fUpwind_r[31]+alpha[21]*fUpwind_r[30]+alpha[22]*fUpwind_r[29]+alpha[5]*fUpwind_r[28]+fUpwind_r[25]*alpha[27]+alpha[3]*fUpwind_r[26]+fUpwind_r[3]*alpha[26]+alpha[12]*fUpwind_r[24]+alpha[13]*fUpwind_r[23]+fUpwind_r[15]*alpha[20]+alpha[7]*fUpwind_r[19]+fUpwind_r[7]*alpha[19]+alpha[8]*fUpwind_r[18]+fUpwind_r[8]*alpha[18]+alpha[0]*fUpwind_r[17]+fUpwind_r[0]*alpha[17]+alpha[11]*fUpwind_r[16]+fUpwind_r[11]*alpha[16]+alpha[1]*fUpwind_r[10]+fUpwind_r[1]*alpha[10]+alpha[2]*fUpwind_r[9]+fUpwind_r[2]*alpha[9]+alpha[4]*fUpwind_r[6]+fUpwind_r[4]*alpha[6]); 
  Ghat_r[18] += 0.1767766952966368*(alpha[13]*fUpwind_r[31]+alpha[20]*fUpwind_r[30]+alpha[5]*fUpwind_r[29]+alpha[22]*fUpwind_r[28]+fUpwind_r[24]*alpha[27]+alpha[2]*fUpwind_r[26]+fUpwind_r[2]*alpha[26]+alpha[12]*fUpwind_r[25]+alpha[14]*fUpwind_r[23]+fUpwind_r[15]*alpha[21]+alpha[6]*fUpwind_r[19]+fUpwind_r[6]*alpha[19]+alpha[0]*fUpwind_r[18]+fUpwind_r[0]*alpha[18]+alpha[8]*fUpwind_r[17]+fUpwind_r[8]*alpha[17]+alpha[10]*fUpwind_r[16]+fUpwind_r[10]*alpha[16]+alpha[1]*fUpwind_r[11]+fUpwind_r[1]*alpha[11]+alpha[3]*fUpwind_r[9]+fUpwind_r[3]*alpha[9]+alpha[4]*fUpwind_r[7]+fUpwind_r[4]*alpha[7]); 
  Ghat_r[19] += 0.1767766952966368*(alpha[12]*fUpwind_r[31]+alpha[5]*fUpwind_r[30]+alpha[20]*fUpwind_r[29]+alpha[21]*fUpwind_r[28]+fUpwind_r[23]*alpha[27]+alpha[1]*fUpwind_r[26]+fUpwind_r[1]*alpha[26]+alpha[13]*fUpwind_r[25]+alpha[14]*fUpwind_r[24]+fUpwind_r[15]*alpha[22]+alpha[0]*fUpwind_r[19]+fUpwind_r[0]*alpha[19]+alpha[6]*fUpwind_r[18]+fUpwind_r[6]*alpha[18]+alpha[7]*fUpwind_r[17]+fUpwind_r[7]*alpha[17]+alpha[9]*fUpwind_r[16]+fUpwind_r[9]*alpha[16]+alpha[2]*fUpwind_r[11]+fUpwind_r[2]*alpha[11]+alpha[3]*fUpwind_r[10]+fUpwind_r[3]*alpha[10]+alpha[4]*fUpwind_r[8]+fUpwind_r[4]*alpha[8]); 
  Ghat_r[20] += 0.1767766952966368*(alpha[11]*fUpwind_r[31]+alpha[18]*fUpwind_r[30]+alpha[19]*fUpwind_r[29]+alpha[4]*fUpwind_r[28]+alpha[3]*fUpwind_r[27]+fUpwind_r[3]*alpha[27]+fUpwind_r[25]*alpha[26]+alpha[9]*fUpwind_r[24]+alpha[10]*fUpwind_r[23]+alpha[7]*fUpwind_r[22]+fUpwind_r[7]*alpha[22]+alpha[8]*fUpwind_r[21]+fUpwind_r[8]*alpha[21]+alpha[0]*fUpwind_r[20]+fUpwind_r[0]*alpha[20]+fUpwind_r[15]*alpha[17]+alpha[14]*fUpwind_r[16]+fUpwind_r[14]*alpha[16]+alpha[1]*fUpwind_r[13]+fUpwind_r[1]*alpha[13]+alpha[2]*fUpwind_r[12]+fUpwind_r[2]*alpha[12]+alpha[5]*fUpwind_r[6]+fUpwind_r[5]*alpha[6]); 
  Ghat_r[21] += 0.1767766952966368*(alpha[10]*fUpwind_r[31]+alpha[17]*fUpwind_r[30]+alpha[4]*fUpwind_r[29]+alpha[19]*fUpwind_r[28]+alpha[2]*fUpwind_r[27]+fUpwind_r[2]*alpha[27]+fUpwind_r[24]*alpha[26]+alpha[9]*fUpwind_r[25]+alpha[11]*fUpwind_r[23]+alpha[6]*fUpwind_r[22]+fUpwind_r[6]*alpha[22]+alpha[0]*fUpwind_r[21]+fUpwind_r[0]*alpha[21]+alpha[8]*fUpwind_r[20]+fUpwind_r[8]*alpha[20]+fUpwind_r[15]*alpha[18]+alpha[13]*fUpwind_r[16]+fUpwind_r[13]*alpha[16]+alpha[1]*fUpwind_r[14]+fUpwind_r[1]*alpha[14]+alpha[3]*fUpwind_r[12]+fUpwind_r[3]*alpha[12]+alpha[5]*fUpwind_r[7]+fUpwind_r[5]*alpha[7]); 
  Ghat_r[22] += 0.1767766952966368*(alpha[9]*fUpwind_r[31]+alpha[4]*fUpwind_r[30]+alpha[17]*fUpwind_r[29]+alpha[18]*fUpwind_r[28]+alpha[1]*fUpwind_r[27]+fUpwind_r[1]*alpha[27]+fUpwind_r[23]*alpha[26]+alpha[10]*fUpwind_r[25]+alpha[11]*fUpwind_r[24]+alpha[0]*fUpwind_r[22]+fUpwind_r[0]*alpha[22]+alpha[6]*fUpwind_r[21]+fUpwind_r[6]*alpha[21]+alpha[7]*fUpwind_r[20]+fUpwind_r[7]*alpha[20]+fUpwind_r[15]*alpha[19]+alpha[12]*fUpwind_r[16]+fUpwind_r[12]*alpha[16]+alpha[2]*fUpwind_r[14]+fUpwind_r[2]*alpha[14]+alpha[3]*fUpwind_r[13]+fUpwind_r[3]*alpha[13]+alpha[5]*fUpwind_r[8]+fUpwind_r[5]*alpha[8]); 
  Ghat_r[23] += 0.1767766952966368*(alpha[8]*fUpwind_r[31]+alpha[16]*fUpwind_r[30]+alpha[3]*fUpwind_r[29]+alpha[2]*fUpwind_r[28]+alpha[19]*fUpwind_r[27]+fUpwind_r[19]*alpha[27]+alpha[22]*fUpwind_r[26]+fUpwind_r[22]*alpha[26]+alpha[7]*fUpwind_r[25]+alpha[6]*fUpwind_r[24]+alpha[0]*fUpwind_r[23]+alpha[11]*fUpwind_r[21]+fUpwind_r[11]*alpha[21]+alpha[10]*fUpwind_r[20]+fUpwind_r[10]*alpha[20]+alpha[14]*fUpwind_r[18]+fUpwind_r[14]*alpha[18]+alpha[13]*fUpwind_r[17]+fUpwind_r[13]*alpha[17]+alpha[1]*fUpwind_r[15]+alpha[4]*fUpwind_r[12]+fUpwind_r[4]*alpha[12]+alpha[5]*fUpwind_r[9]+fUpwind_r[5]*alpha[9]); 
  Ghat_r[24] += 0.1767766952966368*(alpha[7]*fUpwind_r[31]+alpha[3]*fUpwind_r[30]+alpha[16]*fUpwind_r[29]+alpha[1]*fUpwind_r[28]+alpha[18]*fUpwind_r[27]+fUpwind_r[18]*alpha[27]+alpha[21]*fUpwind_r[26]+fUpwind_r[21]*alpha[26]+alpha[8]*fUpwind_r[25]+alpha[0]*fUpwind_r[24]+alpha[6]*fUpwind_r[23]+alpha[11]*fUpwind_r[22]+fUpwind_r[11]*alpha[22]+alpha[9]*fUpwind_r[20]+fUpwind_r[9]*alpha[20]+alpha[14]*fUpwind_r[19]+fUpwind_r[14]*alpha[19]+alpha[12]*fUpwind_r[17]+fUpwind_r[12]*alpha[17]+alpha[2]*fUpwind_r[15]+alpha[4]*fUpwind_r[13]+fUpwind_r[4]*alpha[13]+alpha[5]*fUpwind_r[10]+fUpwind_r[5]*alpha[10]); 
  Ghat_r[25] += 0.1767766952966368*(alpha[6]*fUpwind_r[31]+alpha[2]*fUpwind_r[30]+alpha[1]*fUpwind_r[29]+alpha[16]*fUpwind_r[28]+alpha[17]*fUpwind_r[27]+fUpwind_r[17]*alpha[27]+alpha[20]*fUpwind_r[26]+fUpwind_r[20]*alpha[26]+alpha[0]*fUpwind_r[25]+alpha[8]*fUpwind_r[24]+alpha[7]*fUpwind_r[23]+alpha[10]*fUpwind_r[22]+fUpwind_r[10]*alpha[22]+alpha[9]*fUpwind_r[21]+fUpwind_r[9]*alpha[21]+alpha[13]*fUpwind_r[19]+fUpwind_r[13]*alpha[19]+alpha[12]*fUpwind_r[18]+fUpwind_r[12]*alpha[18]+alpha[3]*fUpwind_r[15]+alpha[4]*fUpwind_r[14]+fUpwind_r[4]*alpha[14]+alpha[5]*fUpwind_r[11]+fUpwind_r[5]*alpha[11]); 
  Ghat_r[26] += 0.1767766952966368*(alpha[5]*fUpwind_r[31]+alpha[12]*fUpwind_r[30]+alpha[13]*fUpwind_r[29]+alpha[14]*fUpwind_r[28]+fUpwind_r[15]*alpha[27]+alpha[0]*fUpwind_r[26]+fUpwind_r[0]*alpha[26]+alpha[20]*fUpwind_r[25]+alpha[21]*fUpwind_r[24]+alpha[22]*fUpwind_r[23]+alpha[1]*fUpwind_r[19]+fUpwind_r[1]*alpha[19]+alpha[2]*fUpwind_r[18]+fUpwind_r[2]*alpha[18]+alpha[3]*fUpwind_r[17]+fUpwind_r[3]*alpha[17]+alpha[4]*fUpwind_r[16]+fUpwind_r[4]*alpha[16]+alpha[6]*fUpwind_r[11]+fUpwind_r[6]*alpha[11]+alpha[7]*fUpwind_r[10]+fUpwind_r[7]*alpha[10]+alpha[8]*fUpwind_r[9]+fUpwind_r[8]*alpha[9]); 
  Ghat_r[27] += 0.1767766952966368*(alpha[4]*fUpwind_r[31]+alpha[9]*fUpwind_r[30]+alpha[10]*fUpwind_r[29]+alpha[11]*fUpwind_r[28]+alpha[0]*fUpwind_r[27]+fUpwind_r[0]*alpha[27]+fUpwind_r[15]*alpha[26]+alpha[17]*fUpwind_r[25]+alpha[18]*fUpwind_r[24]+alpha[19]*fUpwind_r[23]+alpha[1]*fUpwind_r[22]+fUpwind_r[1]*alpha[22]+alpha[2]*fUpwind_r[21]+fUpwind_r[2]*alpha[21]+alpha[3]*fUpwind_r[20]+fUpwind_r[3]*alpha[20]+alpha[5]*fUpwind_r[16]+fUpwind_r[5]*alpha[16]+alpha[6]*fUpwind_r[14]+fUpwind_r[6]*alpha[14]+alpha[7]*fUpwind_r[13]+fUpwind_r[7]*alpha[13]+alpha[8]*fUpwind_r[12]+fUpwind_r[8]*alpha[12]); 
  Ghat_r[28] += 0.1767766952966368*(alpha[3]*fUpwind_r[31]+alpha[7]*fUpwind_r[30]+alpha[8]*fUpwind_r[29]+alpha[0]*fUpwind_r[28]+alpha[11]*fUpwind_r[27]+fUpwind_r[11]*alpha[27]+alpha[14]*fUpwind_r[26]+fUpwind_r[14]*alpha[26]+alpha[16]*fUpwind_r[25]+alpha[1]*fUpwind_r[24]+alpha[2]*fUpwind_r[23]+alpha[18]*fUpwind_r[22]+fUpwind_r[18]*alpha[22]+alpha[19]*fUpwind_r[21]+fUpwind_r[19]*alpha[21]+alpha[4]*fUpwind_r[20]+fUpwind_r[4]*alpha[20]+alpha[5]*fUpwind_r[17]+fUpwind_r[5]*alpha[17]+alpha[6]*fUpwind_r[15]+alpha[9]*fUpwind_r[13]+fUpwind_r[9]*alpha[13]+alpha[10]*fUpwind_r[12]+fUpwind_r[10]*alpha[12]); 
  Ghat_r[29] += 0.1767766952966368*(alpha[2]*fUpwind_r[31]+alpha[6]*fUpwind_r[30]+alpha[0]*fUpwind_r[29]+alpha[8]*fUpwind_r[28]+alpha[10]*fUpwind_r[27]+fUpwind_r[10]*alpha[27]+alpha[13]*fUpwind_r[26]+fUpwind_r[13]*alpha[26]+alpha[1]*fUpwind_r[25]+alpha[16]*fUpwind_r[24]+alpha[3]*fUpwind_r[23]+alpha[17]*fUpwind_r[22]+fUpwind_r[17]*alpha[22]+alpha[4]*fUpwind_r[21]+fUpwind_r[4]*alpha[21]+alpha[19]*fUpwind_r[20]+fUpwind_r[19]*alpha[20]+alpha[5]*fUpwind_r[18]+fUpwind_r[5]*alpha[18]+alpha[7]*fUpwind_r[15]+alpha[9]*fUpwind_r[14]+fUpwind_r[9]*alpha[14]+alpha[11]*fUpwind_r[12]+fUpwind_r[11]*alpha[12]); 
  Ghat_r[30] += 0.1767766952966368*(alpha[1]*fUpwind_r[31]+alpha[0]*fUpwind_r[30]+alpha[6]*fUpwind_r[29]+alpha[7]*fUpwind_r[28]+alpha[9]*fUpwind_r[27]+fUpwind_r[9]*alpha[27]+alpha[12]*fUpwind_r[26]+fUpwind_r[12]*alpha[26]+alpha[2]*fUpwind_r[25]+alpha[3]*fUpwind_r[24]+alpha[16]*fUpwind_r[23]+alpha[4]*fUpwind_r[22]+fUpwind_r[4]*alpha[22]+alpha[17]*fUpwind_r[21]+fUpwind_r[17]*alpha[21]+alpha[18]*fUpwind_r[20]+fUpwind_r[18]*alpha[20]+alpha[5]*fUpwind_r[19]+fUpwind_r[5]*alpha[19]+alpha[8]*fUpwind_r[15]+alpha[10]*fUpwind_r[14]+fUpwind_r[10]*alpha[14]+alpha[11]*fUpwind_r[13]+fUpwind_r[11]*alpha[13]); 
  Ghat_r[31] += 0.1767766952966368*(alpha[0]*fUpwind_r[31]+alpha[1]*fUpwind_r[30]+alpha[2]*fUpwind_r[29]+alpha[3]*fUpwind_r[28]+alpha[4]*fUpwind_r[27]+fUpwind_r[4]*alpha[27]+alpha[5]*fUpwind_r[26]+fUpwind_r[5]*alpha[26]+alpha[6]*fUpwind_r[25]+alpha[7]*fUpwind_r[24]+alpha[8]*fUpwind_r[23]+alpha[9]*fUpwind_r[22]+fUpwind_r[9]*alpha[22]+alpha[10]*fUpwind_r[21]+fUpwind_r[10]*alpha[21]+alpha[11]*fUpwind_r[20]+fUpwind_r[11]*alpha[20]+alpha[12]*fUpwind_r[19]+fUpwind_r[12]*alpha[19]+alpha[13]*fUpwind_r[18]+fUpwind_r[13]*alpha[18]+alpha[14]*fUpwind_r[17]+fUpwind_r[14]*alpha[17]+fUpwind_r[15]*alpha[16]); 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv10; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv10; 
  out[2] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dv10; 
  out[3] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dv10; 
  out[4] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv10; 
  out[5] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*dv10; 
  out[6] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*dv10; 
  out[7] += (0.7071067811865475*Ghat_l[6]-0.7071067811865475*Ghat_r[6])*dv10; 
  out[8] += (0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*dv10; 
  out[9] += (0.7071067811865475*Ghat_l[8]-0.7071067811865475*Ghat_r[8])*dv10; 
  out[10] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv10; 
  out[11] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dv10; 
  out[12] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dv10; 
  out[13] += (0.7071067811865475*Ghat_l[9]-0.7071067811865475*Ghat_r[9])*dv10; 
  out[14] += (0.7071067811865475*Ghat_l[10]-0.7071067811865475*Ghat_r[10])*dv10; 
  out[15] += (0.7071067811865475*Ghat_l[11]-0.7071067811865475*Ghat_r[11])*dv10; 
  out[16] += -1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dv10; 
  out[17] += (0.7071067811865475*Ghat_l[12]-0.7071067811865475*Ghat_r[12])*dv10; 
  out[18] += (0.7071067811865475*Ghat_l[13]-0.7071067811865475*Ghat_r[13])*dv10; 
  out[19] += (0.7071067811865475*Ghat_l[14]-0.7071067811865475*Ghat_r[14])*dv10; 
  out[20] += -1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dv10; 
  out[21] += (0.7071067811865475*Ghat_l[15]-0.7071067811865475*Ghat_r[15])*dv10; 
  out[22] += (0.7071067811865475*Ghat_l[16]-0.7071067811865475*Ghat_r[16])*dv10; 
  out[23] += -1.224744871391589*(Ghat_r[6]+Ghat_l[6])*dv10; 
  out[24] += -1.224744871391589*(Ghat_r[7]+Ghat_l[7])*dv10; 
  out[25] += -1.224744871391589*(Ghat_r[8]+Ghat_l[8])*dv10; 
  out[26] += (0.7071067811865475*Ghat_l[17]-0.7071067811865475*Ghat_r[17])*dv10; 
  out[27] += (0.7071067811865475*Ghat_l[18]-0.7071067811865475*Ghat_r[18])*dv10; 
  out[28] += (0.7071067811865475*Ghat_l[19]-0.7071067811865475*Ghat_r[19])*dv10; 
  out[29] += -1.224744871391589*(Ghat_r[9]+Ghat_l[9])*dv10; 
  out[30] += -1.224744871391589*(Ghat_r[10]+Ghat_l[10])*dv10; 
  out[31] += -1.224744871391589*(Ghat_r[11]+Ghat_l[11])*dv10; 
  out[32] += (0.7071067811865475*Ghat_l[20]-0.7071067811865475*Ghat_r[20])*dv10; 
  out[33] += (0.7071067811865475*Ghat_l[21]-0.7071067811865475*Ghat_r[21])*dv10; 
  out[34] += (0.7071067811865475*Ghat_l[22]-0.7071067811865475*Ghat_r[22])*dv10; 
  out[35] += -1.224744871391589*(Ghat_r[12]+Ghat_l[12])*dv10; 
  out[36] += -1.224744871391589*(Ghat_r[13]+Ghat_l[13])*dv10; 
  out[37] += -1.224744871391589*(Ghat_r[14]+Ghat_l[14])*dv10; 
  out[38] += (0.7071067811865475*Ghat_l[23]-0.7071067811865475*Ghat_r[23])*dv10; 
  out[39] += (0.7071067811865475*Ghat_l[24]-0.7071067811865475*Ghat_r[24])*dv10; 
  out[40] += (0.7071067811865475*Ghat_l[25]-0.7071067811865475*Ghat_r[25])*dv10; 
  out[41] += -1.224744871391589*(Ghat_r[15]+Ghat_l[15])*dv10; 
  out[42] += -1.224744871391589*(Ghat_r[16]+Ghat_l[16])*dv10; 
  out[43] += (0.7071067811865475*Ghat_l[26]-0.7071067811865475*Ghat_r[26])*dv10; 
  out[44] += -1.224744871391589*(Ghat_r[17]+Ghat_l[17])*dv10; 
  out[45] += -1.224744871391589*(Ghat_r[18]+Ghat_l[18])*dv10; 
  out[46] += -1.224744871391589*(Ghat_r[19]+Ghat_l[19])*dv10; 
  out[47] += (0.7071067811865475*Ghat_l[27]-0.7071067811865475*Ghat_r[27])*dv10; 
  out[48] += -1.224744871391589*(Ghat_r[20]+Ghat_l[20])*dv10; 
  out[49] += -1.224744871391589*(Ghat_r[21]+Ghat_l[21])*dv10; 
  out[50] += -1.224744871391589*(Ghat_r[22]+Ghat_l[22])*dv10; 
  out[51] += (0.7071067811865475*Ghat_l[28]-0.7071067811865475*Ghat_r[28])*dv10; 
  out[52] += (0.7071067811865475*Ghat_l[29]-0.7071067811865475*Ghat_r[29])*dv10; 
  out[53] += (0.7071067811865475*Ghat_l[30]-0.7071067811865475*Ghat_r[30])*dv10; 
  out[54] += -1.224744871391589*(Ghat_r[23]+Ghat_l[23])*dv10; 
  out[55] += -1.224744871391589*(Ghat_r[24]+Ghat_l[24])*dv10; 
  out[56] += -1.224744871391589*(Ghat_r[25]+Ghat_l[25])*dv10; 
  out[57] += -1.224744871391589*(Ghat_r[26]+Ghat_l[26])*dv10; 
  out[58] += -1.224744871391589*(Ghat_r[27]+Ghat_l[27])*dv10; 
  out[59] += (0.7071067811865475*Ghat_l[31]-0.7071067811865475*Ghat_r[31])*dv10; 
  out[60] += -1.224744871391589*(Ghat_r[28]+Ghat_l[28])*dv10; 
  out[61] += -1.224744871391589*(Ghat_r[29]+Ghat_l[29])*dv10; 
  out[62] += -1.224744871391589*(Ghat_r[30]+Ghat_l[30])*dv10; 
  out[63] += -1.224744871391589*(Ghat_r[31]+Ghat_l[31])*dv10; 

} 
