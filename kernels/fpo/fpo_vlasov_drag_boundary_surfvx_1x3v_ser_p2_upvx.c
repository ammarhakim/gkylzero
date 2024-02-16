#include <gkyl_fpo_vlasov_kernels.h> 
#include <gkyl_basis_ser_4x_p2_surfx2_eval_quad.h> 
#include <gkyl_basis_ser_4x_p2_upwind_quad_to_modal.h> 

 
GKYL_CU_DH void fpo_vlasov_drag_boundary_surfvx_1x3v_ser_p2_upvx(const double *dxv, const double* drag_coeff_stencil[9], const double* f_stencil[9], double* out) { 
  // dxv[NDIM]: Cell spacing in each direction. 
  // drag_coeff_stencil[9]: 9-cell stencil of drag coefficient. 
  // f_stencil[9]: 9-cell stencil of distribution function. 
  // out: Incremented output. 


  double dv1 = 2.0/dxv[1]; 
 
  const double* aL = &drag_coeff_stencil[0][0]; 
  const double* aC = &drag_coeff_stencil[1][0]; 
  const double* fL = f_stencil[0]; 
  const double* fC = f_stencil[1]; 

  double fUpwindQuad[27] = {0.0};
  double fUpwind[20] = {0.0};
  double alphaDragSurf[20] = {0.0}; 
  double Ghat[20] = {0.0}; 

  alphaDragSurf[0] = 0.34587411908091625*aL[12]+0.34587411908091625*aC[12]+0.49755260400283263*aL[2]-0.49755260400283263*aC[2]+0.3535533905932737*aL[0]+0.3535533905932737*aC[0]; 
  alphaDragSurf[1] = 0.34587411908091625*aL[20]+0.34587411908091625*aC[20]+0.49755260400283263*aL[5]-0.49755260400283263*aC[5]+0.3535533905932737*aL[1]+0.3535533905932737*aC[1]; 
  alphaDragSurf[2] = 0.34587411908091625*aL[22]+0.34587411908091625*aC[22]+0.49755260400283263*aL[7]-0.49755260400283263*aC[7]+0.3535533905932737*aL[3]+0.3535533905932737*aC[3]; 
  alphaDragSurf[3] = 0.34587411908091625*aL[26]+0.34587411908091625*aC[26]+0.49755260400283263*aL[9]-0.49755260400283263*aC[9]+0.3535533905932737*aL[4]+0.3535533905932737*aC[4]; 
  alphaDragSurf[4] = 0.34587411908091625*aL[33]+0.34587411908091625*aC[33]+0.49755260400283263*aL[15]-0.49755260400283263*aC[15]+0.3535533905932737*aL[6]+0.3535533905932737*aC[6]; 
  alphaDragSurf[5] = 0.34587411908091625*aL[36]+0.34587411908091625*aC[36]+0.49755260400283263*aL[16]-0.49755260400283263*aC[16]+0.3535533905932737*aL[8]+0.3535533905932737*aC[8]; 
  alphaDragSurf[6] = 0.34587411908091625*aL[38]+0.34587411908091625*aC[38]+0.49755260400283263*aL[18]-0.49755260400283263*aC[18]+0.3535533905932737*aL[10]+0.3535533905932737*aC[10]; 
  alphaDragSurf[7] = 0.49755260400283263*aL[19]-0.49755260400283263*aC[19]+0.3535533905932737*aL[11]+0.3535533905932737*aC[11]; 
  alphaDragSurf[8] = 0.49755260400283263*aL[24]-0.49755260400283263*aC[24]+0.3535533905932737*aL[13]+0.3535533905932737*aC[13]; 
  alphaDragSurf[9] = 0.49755260400283263*aL[29]-0.49755260400283263*aC[29]+0.3535533905932737*aL[14]+0.3535533905932737*aC[14]; 
  alphaDragSurf[10] = 0.34587411908091625*aL[45]+0.34587411908091625*aC[45]+0.49755260400283263*aL[31]-0.49755260400283263*aC[31]+0.3535533905932737*aL[17]+0.3535533905932737*aC[17]; 
  alphaDragSurf[11] = 0.49755260400283263*aL[32]-0.49755260400283263*aC[32]+0.3535533905932737*aL[21]+0.3535533905932737*aC[21]; 
  alphaDragSurf[12] = 0.49755260400283263*aL[34]-0.49755260400283263*aC[34]+0.3535533905932737*aL[23]+0.3535533905932737*aC[23]; 
  alphaDragSurf[13] = 0.49755260400283263*aL[35]-0.49755260400283263*aC[35]+0.3535533905932737*aL[25]+0.3535533905932737*aC[25]; 
  alphaDragSurf[14] = 0.49755260400283263*aL[40]-0.49755260400283263*aC[40]+0.3535533905932737*aL[27]+0.3535533905932737*aC[27]; 
  alphaDragSurf[15] = 0.49755260400283263*aL[41]-0.49755260400283263*aC[41]+0.3535533905932737*aL[28]+0.3535533905932737*aC[28]; 
  alphaDragSurf[16] = 0.49755260400283263*aL[43]-0.49755260400283263*aC[43]+0.3535533905932737*aL[30]+0.3535533905932737*aC[30]; 
  alphaDragSurf[17] = 0.49755260400283263*aL[44]-0.49755260400283263*aC[44]+0.3535533905932737*aL[37]+0.3535533905932737*aC[37]; 
  alphaDragSurf[18] = 0.49755260400283263*aL[46]-0.49755260400283263*aC[46]+0.3535533905932737*aL[39]+0.3535533905932737*aC[39]; 
  alphaDragSurf[19] = 0.49755260400283263*aL[47]-0.49755260400283263*aC[47]+0.3535533905932737*aL[42]+0.3535533905932737*aC[42]; 

  if (0.5692099788303082*(alphaDragSurf[19]+alphaDragSurf[18]+alphaDragSurf[17])-0.42426406871192807*(alphaDragSurf[16]+alphaDragSurf[15])-0.42426406871192845*(alphaDragSurf[14]+alphaDragSurf[13])-0.42426406871192807*alphaDragSurf[12]-0.42426406871192845*alphaDragSurf[11]-0.853814968245462*alphaDragSurf[10]+0.3162277660168379*(alphaDragSurf[9]+alphaDragSurf[8]+alphaDragSurf[7])+0.6363961030678926*(alphaDragSurf[6]+alphaDragSurf[5]+alphaDragSurf[4])-0.4743416490252568*(alphaDragSurf[3]+alphaDragSurf[2]+alphaDragSurf[1])+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[0] = ser_4x_p2_surfx2_eval_quad_node_0_r(fL); 
  } else { 
    fUpwindQuad[0] = ser_4x_p2_surfx2_eval_quad_node_0_l(fC); 
  } 
  if (-(0.711512473537885*alphaDragSurf[19])+0.5303300858899104*(alphaDragSurf[16]+alphaDragSurf[15])-0.42426406871192807*alphaDragSurf[12]-0.42426406871192845*alphaDragSurf[11]-0.3952847075210473*alphaDragSurf[9]+0.3162277660168379*(alphaDragSurf[8]+alphaDragSurf[7])+0.6363961030678926*alphaDragSurf[4]-0.4743416490252568*(alphaDragSurf[2]+alphaDragSurf[1])+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[1] = ser_4x_p2_surfx2_eval_quad_node_1_r(fL); 
  } else { 
    fUpwindQuad[1] = ser_4x_p2_surfx2_eval_quad_node_1_l(fC); 
  } 
  if (0.5692099788303082*alphaDragSurf[19]-0.5692099788303082*(alphaDragSurf[18]+alphaDragSurf[17])-0.42426406871192807*(alphaDragSurf[16]+alphaDragSurf[15])+0.42426406871192845*(alphaDragSurf[14]+alphaDragSurf[13])-0.42426406871192807*alphaDragSurf[12]-0.42426406871192845*alphaDragSurf[11]+0.853814968245462*alphaDragSurf[10]+0.3162277660168379*(alphaDragSurf[9]+alphaDragSurf[8]+alphaDragSurf[7])-0.6363961030678926*(alphaDragSurf[6]+alphaDragSurf[5])+0.6363961030678926*alphaDragSurf[4]+0.4743416490252568*alphaDragSurf[3]-0.4743416490252568*(alphaDragSurf[2]+alphaDragSurf[1])+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[2] = ser_4x_p2_surfx2_eval_quad_node_2_r(fL); 
  } else { 
    fUpwindQuad[2] = ser_4x_p2_surfx2_eval_quad_node_2_l(fC); 
  } 
  if (-(0.711512473537885*alphaDragSurf[18])-0.42426406871192807*alphaDragSurf[15]+0.5303300858899104*alphaDragSurf[14]-0.42426406871192845*alphaDragSurf[13]+0.5303300858899104*alphaDragSurf[12]+0.3162277660168379*alphaDragSurf[9]-0.3952847075210473*alphaDragSurf[8]+0.3162277660168379*alphaDragSurf[7]+0.6363961030678926*alphaDragSurf[5]-0.4743416490252568*(alphaDragSurf[3]+alphaDragSurf[1])+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[3] = ser_4x_p2_surfx2_eval_quad_node_3_r(fL); 
  } else { 
    fUpwindQuad[3] = ser_4x_p2_surfx2_eval_quad_node_3_l(fC); 
  } 
  if (0.5303300858899104*(alphaDragSurf[15]+alphaDragSurf[12])-0.3952847075210473*(alphaDragSurf[9]+alphaDragSurf[8])+0.3162277660168379*alphaDragSurf[7]-0.4743416490252568*alphaDragSurf[1]+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[4] = ser_4x_p2_surfx2_eval_quad_node_4_r(fL); 
  } else { 
    fUpwindQuad[4] = ser_4x_p2_surfx2_eval_quad_node_4_l(fC); 
  } 
  if (0.711512473537885*alphaDragSurf[18]-0.42426406871192807*alphaDragSurf[15]-0.5303300858899104*alphaDragSurf[14]+0.42426406871192845*alphaDragSurf[13]+0.5303300858899104*alphaDragSurf[12]+0.3162277660168379*alphaDragSurf[9]-0.3952847075210473*alphaDragSurf[8]+0.3162277660168379*alphaDragSurf[7]-0.6363961030678926*alphaDragSurf[5]+0.4743416490252568*alphaDragSurf[3]-0.4743416490252568*alphaDragSurf[1]+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[5] = ser_4x_p2_surfx2_eval_quad_node_5_r(fL); 
  } else { 
    fUpwindQuad[5] = ser_4x_p2_surfx2_eval_quad_node_5_l(fC); 
  } 
  if (-(0.5692099788303082*alphaDragSurf[19])+0.5692099788303082*alphaDragSurf[18]-0.5692099788303082*alphaDragSurf[17]+0.42426406871192807*alphaDragSurf[16]-0.42426406871192807*alphaDragSurf[15]-0.42426406871192845*(alphaDragSurf[14]+alphaDragSurf[13])-0.42426406871192807*alphaDragSurf[12]+0.42426406871192845*alphaDragSurf[11]+0.853814968245462*alphaDragSurf[10]+0.3162277660168379*(alphaDragSurf[9]+alphaDragSurf[8]+alphaDragSurf[7])-0.6363961030678926*alphaDragSurf[6]+0.6363961030678926*alphaDragSurf[5]-0.6363961030678926*alphaDragSurf[4]-0.4743416490252568*alphaDragSurf[3]+0.4743416490252568*alphaDragSurf[2]-0.4743416490252568*alphaDragSurf[1]+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[6] = ser_4x_p2_surfx2_eval_quad_node_6_r(fL); 
  } else { 
    fUpwindQuad[6] = ser_4x_p2_surfx2_eval_quad_node_6_l(fC); 
  } 
  if (0.711512473537885*alphaDragSurf[19]-0.5303300858899104*alphaDragSurf[16]+0.5303300858899104*alphaDragSurf[15]-0.42426406871192807*alphaDragSurf[12]+0.42426406871192845*alphaDragSurf[11]-0.3952847075210473*alphaDragSurf[9]+0.3162277660168379*(alphaDragSurf[8]+alphaDragSurf[7])-0.6363961030678926*alphaDragSurf[4]+0.4743416490252568*alphaDragSurf[2]-0.4743416490252568*alphaDragSurf[1]+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[7] = ser_4x_p2_surfx2_eval_quad_node_7_r(fL); 
  } else { 
    fUpwindQuad[7] = ser_4x_p2_surfx2_eval_quad_node_7_l(fC); 
  } 
  if (-(0.5692099788303082*(alphaDragSurf[19]+alphaDragSurf[18]))+0.5692099788303082*alphaDragSurf[17]+0.42426406871192807*alphaDragSurf[16]-0.42426406871192807*alphaDragSurf[15]+0.42426406871192845*(alphaDragSurf[14]+alphaDragSurf[13])-0.42426406871192807*alphaDragSurf[12]+0.42426406871192845*alphaDragSurf[11]-0.853814968245462*alphaDragSurf[10]+0.3162277660168379*(alphaDragSurf[9]+alphaDragSurf[8]+alphaDragSurf[7])+0.6363961030678926*alphaDragSurf[6]-0.6363961030678926*(alphaDragSurf[5]+alphaDragSurf[4])+0.4743416490252568*(alphaDragSurf[3]+alphaDragSurf[2])-0.4743416490252568*alphaDragSurf[1]+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[8] = ser_4x_p2_surfx2_eval_quad_node_8_r(fL); 
  } else { 
    fUpwindQuad[8] = ser_4x_p2_surfx2_eval_quad_node_8_l(fC); 
  } 
  if (-(0.711512473537885*alphaDragSurf[17])-0.42426406871192807*alphaDragSurf[16]-0.42426406871192845*alphaDragSurf[14]+0.5303300858899104*(alphaDragSurf[13]+alphaDragSurf[11])+0.3162277660168379*(alphaDragSurf[9]+alphaDragSurf[8])-0.3952847075210473*alphaDragSurf[7]+0.6363961030678926*alphaDragSurf[6]-0.4743416490252568*(alphaDragSurf[3]+alphaDragSurf[2])+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[9] = ser_4x_p2_surfx2_eval_quad_node_9_r(fL); 
  } else { 
    fUpwindQuad[9] = ser_4x_p2_surfx2_eval_quad_node_9_l(fC); 
  } 
  if (0.5303300858899104*(alphaDragSurf[16]+alphaDragSurf[11])-0.3952847075210473*alphaDragSurf[9]+0.3162277660168379*alphaDragSurf[8]-0.3952847075210473*alphaDragSurf[7]-0.4743416490252568*alphaDragSurf[2]+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[10] = ser_4x_p2_surfx2_eval_quad_node_10_r(fL); 
  } else { 
    fUpwindQuad[10] = ser_4x_p2_surfx2_eval_quad_node_10_l(fC); 
  } 
  if (0.711512473537885*alphaDragSurf[17]-0.42426406871192807*alphaDragSurf[16]+0.42426406871192845*alphaDragSurf[14]-0.5303300858899104*alphaDragSurf[13]+0.5303300858899104*alphaDragSurf[11]+0.3162277660168379*(alphaDragSurf[9]+alphaDragSurf[8])-0.3952847075210473*alphaDragSurf[7]-0.6363961030678926*alphaDragSurf[6]+0.4743416490252568*alphaDragSurf[3]-0.4743416490252568*alphaDragSurf[2]+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[11] = ser_4x_p2_surfx2_eval_quad_node_11_r(fL); 
  } else { 
    fUpwindQuad[11] = ser_4x_p2_surfx2_eval_quad_node_11_l(fC); 
  } 
  if (0.5303300858899104*(alphaDragSurf[14]+alphaDragSurf[13])+0.3162277660168379*alphaDragSurf[9]-0.3952847075210473*(alphaDragSurf[8]+alphaDragSurf[7])-0.4743416490252568*alphaDragSurf[3]+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[12] = ser_4x_p2_surfx2_eval_quad_node_12_r(fL); 
  } else { 
    fUpwindQuad[12] = ser_4x_p2_surfx2_eval_quad_node_12_l(fC); 
  } 
  if (0.3535533905932737*alphaDragSurf[0]-0.3952847075210473*(alphaDragSurf[9]+alphaDragSurf[8]+alphaDragSurf[7]) > 0) { 
    fUpwindQuad[13] = ser_4x_p2_surfx2_eval_quad_node_13_r(fL); 
  } else { 
    fUpwindQuad[13] = ser_4x_p2_surfx2_eval_quad_node_13_l(fC); 
  } 
  if (-(0.5303300858899104*(alphaDragSurf[14]+alphaDragSurf[13]))+0.3162277660168379*alphaDragSurf[9]-0.3952847075210473*(alphaDragSurf[8]+alphaDragSurf[7])+0.4743416490252568*alphaDragSurf[3]+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[14] = ser_4x_p2_surfx2_eval_quad_node_14_r(fL); 
  } else { 
    fUpwindQuad[14] = ser_4x_p2_surfx2_eval_quad_node_14_l(fC); 
  } 
  if (0.711512473537885*alphaDragSurf[17]+0.42426406871192807*alphaDragSurf[16]-0.42426406871192845*alphaDragSurf[14]+0.5303300858899104*alphaDragSurf[13]-0.5303300858899104*alphaDragSurf[11]+0.3162277660168379*(alphaDragSurf[9]+alphaDragSurf[8])-0.3952847075210473*alphaDragSurf[7]-0.6363961030678926*alphaDragSurf[6]-0.4743416490252568*alphaDragSurf[3]+0.4743416490252568*alphaDragSurf[2]+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[15] = ser_4x_p2_surfx2_eval_quad_node_15_r(fL); 
  } else { 
    fUpwindQuad[15] = ser_4x_p2_surfx2_eval_quad_node_15_l(fC); 
  } 
  if (-(0.5303300858899104*(alphaDragSurf[16]+alphaDragSurf[11]))-0.3952847075210473*alphaDragSurf[9]+0.3162277660168379*alphaDragSurf[8]-0.3952847075210473*alphaDragSurf[7]+0.4743416490252568*alphaDragSurf[2]+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[16] = ser_4x_p2_surfx2_eval_quad_node_16_r(fL); 
  } else { 
    fUpwindQuad[16] = ser_4x_p2_surfx2_eval_quad_node_16_l(fC); 
  } 
  if (-(0.711512473537885*alphaDragSurf[17])+0.42426406871192807*alphaDragSurf[16]+0.42426406871192845*alphaDragSurf[14]-0.5303300858899104*(alphaDragSurf[13]+alphaDragSurf[11])+0.3162277660168379*(alphaDragSurf[9]+alphaDragSurf[8])-0.3952847075210473*alphaDragSurf[7]+0.6363961030678926*alphaDragSurf[6]+0.4743416490252568*(alphaDragSurf[3]+alphaDragSurf[2])+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[17] = ser_4x_p2_surfx2_eval_quad_node_17_r(fL); 
  } else { 
    fUpwindQuad[17] = ser_4x_p2_surfx2_eval_quad_node_17_l(fC); 
  } 
  if (-(0.5692099788303082*(alphaDragSurf[19]+alphaDragSurf[18]))+0.5692099788303082*alphaDragSurf[17]-0.42426406871192807*alphaDragSurf[16]+0.42426406871192807*alphaDragSurf[15]-0.42426406871192845*(alphaDragSurf[14]+alphaDragSurf[13])+0.42426406871192807*alphaDragSurf[12]-0.42426406871192845*alphaDragSurf[11]+0.853814968245462*alphaDragSurf[10]+0.3162277660168379*(alphaDragSurf[9]+alphaDragSurf[8]+alphaDragSurf[7])+0.6363961030678926*alphaDragSurf[6]-0.6363961030678926*(alphaDragSurf[5]+alphaDragSurf[4])-0.4743416490252568*(alphaDragSurf[3]+alphaDragSurf[2])+0.4743416490252568*alphaDragSurf[1]+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[18] = ser_4x_p2_surfx2_eval_quad_node_18_r(fL); 
  } else { 
    fUpwindQuad[18] = ser_4x_p2_surfx2_eval_quad_node_18_l(fC); 
  } 
  if (0.711512473537885*alphaDragSurf[19]+0.5303300858899104*alphaDragSurf[16]-0.5303300858899104*alphaDragSurf[15]+0.42426406871192807*alphaDragSurf[12]-0.42426406871192845*alphaDragSurf[11]-0.3952847075210473*alphaDragSurf[9]+0.3162277660168379*(alphaDragSurf[8]+alphaDragSurf[7])-0.6363961030678926*alphaDragSurf[4]-0.4743416490252568*alphaDragSurf[2]+0.4743416490252568*alphaDragSurf[1]+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[19] = ser_4x_p2_surfx2_eval_quad_node_19_r(fL); 
  } else { 
    fUpwindQuad[19] = ser_4x_p2_surfx2_eval_quad_node_19_l(fC); 
  } 
  if (-(0.5692099788303082*alphaDragSurf[19])+0.5692099788303082*alphaDragSurf[18]-0.5692099788303082*alphaDragSurf[17]-0.42426406871192807*alphaDragSurf[16]+0.42426406871192807*alphaDragSurf[15]+0.42426406871192845*(alphaDragSurf[14]+alphaDragSurf[13])+0.42426406871192807*alphaDragSurf[12]-0.42426406871192845*alphaDragSurf[11]-0.853814968245462*alphaDragSurf[10]+0.3162277660168379*(alphaDragSurf[9]+alphaDragSurf[8]+alphaDragSurf[7])-0.6363961030678926*alphaDragSurf[6]+0.6363961030678926*alphaDragSurf[5]-0.6363961030678926*alphaDragSurf[4]+0.4743416490252568*alphaDragSurf[3]-0.4743416490252568*alphaDragSurf[2]+0.4743416490252568*alphaDragSurf[1]+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[20] = ser_4x_p2_surfx2_eval_quad_node_20_r(fL); 
  } else { 
    fUpwindQuad[20] = ser_4x_p2_surfx2_eval_quad_node_20_l(fC); 
  } 
  if (0.711512473537885*alphaDragSurf[18]+0.42426406871192807*alphaDragSurf[15]+0.5303300858899104*alphaDragSurf[14]-0.42426406871192845*alphaDragSurf[13]-0.5303300858899104*alphaDragSurf[12]+0.3162277660168379*alphaDragSurf[9]-0.3952847075210473*alphaDragSurf[8]+0.3162277660168379*alphaDragSurf[7]-0.6363961030678926*alphaDragSurf[5]-0.4743416490252568*alphaDragSurf[3]+0.4743416490252568*alphaDragSurf[1]+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[21] = ser_4x_p2_surfx2_eval_quad_node_21_r(fL); 
  } else { 
    fUpwindQuad[21] = ser_4x_p2_surfx2_eval_quad_node_21_l(fC); 
  } 
  if (-(0.5303300858899104*(alphaDragSurf[15]+alphaDragSurf[12]))-0.3952847075210473*(alphaDragSurf[9]+alphaDragSurf[8])+0.3162277660168379*alphaDragSurf[7]+0.4743416490252568*alphaDragSurf[1]+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[22] = ser_4x_p2_surfx2_eval_quad_node_22_r(fL); 
  } else { 
    fUpwindQuad[22] = ser_4x_p2_surfx2_eval_quad_node_22_l(fC); 
  } 
  if (-(0.711512473537885*alphaDragSurf[18])+0.42426406871192807*alphaDragSurf[15]-0.5303300858899104*alphaDragSurf[14]+0.42426406871192845*alphaDragSurf[13]-0.5303300858899104*alphaDragSurf[12]+0.3162277660168379*alphaDragSurf[9]-0.3952847075210473*alphaDragSurf[8]+0.3162277660168379*alphaDragSurf[7]+0.6363961030678926*alphaDragSurf[5]+0.4743416490252568*(alphaDragSurf[3]+alphaDragSurf[1])+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[23] = ser_4x_p2_surfx2_eval_quad_node_23_r(fL); 
  } else { 
    fUpwindQuad[23] = ser_4x_p2_surfx2_eval_quad_node_23_l(fC); 
  } 
  if (0.5692099788303082*alphaDragSurf[19]-0.5692099788303082*(alphaDragSurf[18]+alphaDragSurf[17])+0.42426406871192807*(alphaDragSurf[16]+alphaDragSurf[15])-0.42426406871192845*(alphaDragSurf[14]+alphaDragSurf[13])+0.42426406871192807*alphaDragSurf[12]+0.42426406871192845*alphaDragSurf[11]-0.853814968245462*alphaDragSurf[10]+0.3162277660168379*(alphaDragSurf[9]+alphaDragSurf[8]+alphaDragSurf[7])-0.6363961030678926*(alphaDragSurf[6]+alphaDragSurf[5])+0.6363961030678926*alphaDragSurf[4]-0.4743416490252568*alphaDragSurf[3]+0.4743416490252568*(alphaDragSurf[2]+alphaDragSurf[1])+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[24] = ser_4x_p2_surfx2_eval_quad_node_24_r(fL); 
  } else { 
    fUpwindQuad[24] = ser_4x_p2_surfx2_eval_quad_node_24_l(fC); 
  } 
  if (-(0.711512473537885*alphaDragSurf[19])-0.5303300858899104*(alphaDragSurf[16]+alphaDragSurf[15])+0.42426406871192807*alphaDragSurf[12]+0.42426406871192845*alphaDragSurf[11]-0.3952847075210473*alphaDragSurf[9]+0.3162277660168379*(alphaDragSurf[8]+alphaDragSurf[7])+0.6363961030678926*alphaDragSurf[4]+0.4743416490252568*(alphaDragSurf[2]+alphaDragSurf[1])+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[25] = ser_4x_p2_surfx2_eval_quad_node_25_r(fL); 
  } else { 
    fUpwindQuad[25] = ser_4x_p2_surfx2_eval_quad_node_25_l(fC); 
  } 
  if (0.5692099788303082*(alphaDragSurf[19]+alphaDragSurf[18]+alphaDragSurf[17])+0.42426406871192807*(alphaDragSurf[16]+alphaDragSurf[15])+0.42426406871192845*(alphaDragSurf[14]+alphaDragSurf[13])+0.42426406871192807*alphaDragSurf[12]+0.42426406871192845*alphaDragSurf[11]+0.853814968245462*alphaDragSurf[10]+0.3162277660168379*(alphaDragSurf[9]+alphaDragSurf[8]+alphaDragSurf[7])+0.6363961030678926*(alphaDragSurf[6]+alphaDragSurf[5]+alphaDragSurf[4])+0.4743416490252568*(alphaDragSurf[3]+alphaDragSurf[2]+alphaDragSurf[1])+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[26] = ser_4x_p2_surfx2_eval_quad_node_26_r(fL); 
  } else { 
    fUpwindQuad[26] = ser_4x_p2_surfx2_eval_quad_node_26_l(fC); 
  } 
  ser_4x_p2_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.3535533905932737*alphaDragSurf[19]*fUpwind[19]+0.3535533905932737*alphaDragSurf[18]*fUpwind[18]+0.3535533905932737*alphaDragSurf[17]*fUpwind[17]+0.3535533905932737*alphaDragSurf[16]*fUpwind[16]+0.3535533905932737*alphaDragSurf[15]*fUpwind[15]+0.3535533905932737*alphaDragSurf[14]*fUpwind[14]+0.3535533905932737*alphaDragSurf[13]*fUpwind[13]+0.3535533905932737*alphaDragSurf[12]*fUpwind[12]+0.3535533905932737*alphaDragSurf[11]*fUpwind[11]+0.3535533905932737*alphaDragSurf[10]*fUpwind[10]+0.3535533905932737*alphaDragSurf[9]*fUpwind[9]+0.3535533905932737*alphaDragSurf[8]*fUpwind[8]+0.3535533905932737*alphaDragSurf[7]*fUpwind[7]+0.3535533905932737*alphaDragSurf[6]*fUpwind[6]+0.3535533905932737*alphaDragSurf[5]*fUpwind[5]+0.3535533905932737*alphaDragSurf[4]*fUpwind[4]+0.3535533905932737*alphaDragSurf[3]*fUpwind[3]+0.3535533905932737*alphaDragSurf[2]*fUpwind[2]+0.3535533905932737*alphaDragSurf[1]*fUpwind[1]+0.3535533905932737*alphaDragSurf[0]*fUpwind[0]; 
  Ghat[1] = 0.3535533905932737*alphaDragSurf[16]*fUpwind[19]+0.3535533905932737*fUpwind[16]*alphaDragSurf[19]+0.3535533905932737*alphaDragSurf[14]*fUpwind[18]+0.3535533905932737*fUpwind[14]*alphaDragSurf[18]+0.3162277660168379*alphaDragSurf[10]*fUpwind[17]+0.3162277660168379*fUpwind[10]*alphaDragSurf[17]+0.3535533905932737*alphaDragSurf[9]*fUpwind[15]+0.3535533905932737*fUpwind[9]*alphaDragSurf[15]+0.31622776601683794*alphaDragSurf[5]*fUpwind[13]+0.31622776601683794*fUpwind[5]*alphaDragSurf[13]+0.3535533905932737*alphaDragSurf[8]*fUpwind[12]+0.3535533905932737*fUpwind[8]*alphaDragSurf[12]+0.31622776601683794*alphaDragSurf[4]*fUpwind[11]+0.31622776601683794*fUpwind[4]*alphaDragSurf[11]+0.3535533905932737*alphaDragSurf[6]*fUpwind[10]+0.3535533905932737*fUpwind[6]*alphaDragSurf[10]+0.3162277660168379*alphaDragSurf[1]*fUpwind[7]+0.3162277660168379*fUpwind[1]*alphaDragSurf[7]+0.3535533905932737*alphaDragSurf[3]*fUpwind[5]+0.3535533905932737*fUpwind[3]*alphaDragSurf[5]+0.3535533905932737*alphaDragSurf[2]*fUpwind[4]+0.3535533905932737*fUpwind[2]*alphaDragSurf[4]+0.3535533905932737*alphaDragSurf[0]*fUpwind[1]+0.3535533905932737*fUpwind[0]*alphaDragSurf[1]; 
  Ghat[2] = 0.3535533905932737*alphaDragSurf[15]*fUpwind[19]+0.3535533905932737*fUpwind[15]*alphaDragSurf[19]+0.3162277660168379*alphaDragSurf[10]*fUpwind[18]+0.3162277660168379*fUpwind[10]*alphaDragSurf[18]+0.3535533905932737*alphaDragSurf[13]*fUpwind[17]+0.3535533905932737*fUpwind[13]*alphaDragSurf[17]+0.3535533905932737*alphaDragSurf[9]*fUpwind[16]+0.3535533905932737*fUpwind[9]*alphaDragSurf[16]+0.31622776601683794*alphaDragSurf[6]*fUpwind[14]+0.31622776601683794*fUpwind[6]*alphaDragSurf[14]+0.31622776601683794*alphaDragSurf[4]*fUpwind[12]+0.31622776601683794*fUpwind[4]*alphaDragSurf[12]+0.3535533905932737*alphaDragSurf[7]*fUpwind[11]+0.3535533905932737*fUpwind[7]*alphaDragSurf[11]+0.3535533905932737*alphaDragSurf[5]*fUpwind[10]+0.3535533905932737*fUpwind[5]*alphaDragSurf[10]+0.3162277660168379*alphaDragSurf[2]*fUpwind[8]+0.3162277660168379*fUpwind[2]*alphaDragSurf[8]+0.3535533905932737*alphaDragSurf[3]*fUpwind[6]+0.3535533905932737*fUpwind[3]*alphaDragSurf[6]+0.3535533905932737*alphaDragSurf[1]*fUpwind[4]+0.3535533905932737*fUpwind[1]*alphaDragSurf[4]+0.3535533905932737*alphaDragSurf[0]*fUpwind[2]+0.3535533905932737*fUpwind[0]*alphaDragSurf[2]; 
  Ghat[3] = 0.3162277660168379*alphaDragSurf[10]*fUpwind[19]+0.3162277660168379*fUpwind[10]*alphaDragSurf[19]+0.3535533905932737*alphaDragSurf[12]*fUpwind[18]+0.3535533905932737*fUpwind[12]*alphaDragSurf[18]+0.3535533905932737*alphaDragSurf[11]*fUpwind[17]+0.3535533905932737*fUpwind[11]*alphaDragSurf[17]+0.31622776601683794*alphaDragSurf[6]*fUpwind[16]+0.31622776601683794*fUpwind[6]*alphaDragSurf[16]+0.31622776601683794*alphaDragSurf[5]*fUpwind[15]+0.31622776601683794*fUpwind[5]*alphaDragSurf[15]+0.3535533905932737*alphaDragSurf[8]*fUpwind[14]+0.3535533905932737*fUpwind[8]*alphaDragSurf[14]+0.3535533905932737*alphaDragSurf[7]*fUpwind[13]+0.3535533905932737*fUpwind[7]*alphaDragSurf[13]+0.3535533905932737*alphaDragSurf[4]*fUpwind[10]+0.3535533905932737*fUpwind[4]*alphaDragSurf[10]+0.3162277660168379*alphaDragSurf[3]*fUpwind[9]+0.3162277660168379*fUpwind[3]*alphaDragSurf[9]+0.3535533905932737*alphaDragSurf[2]*fUpwind[6]+0.3535533905932737*fUpwind[2]*alphaDragSurf[6]+0.3535533905932737*alphaDragSurf[1]*fUpwind[5]+0.3535533905932737*fUpwind[1]*alphaDragSurf[5]+0.3535533905932737*alphaDragSurf[0]*fUpwind[3]+0.3535533905932737*fUpwind[0]*alphaDragSurf[3]; 
  Ghat[4] = 0.3535533905932737*alphaDragSurf[9]*fUpwind[19]+0.3535533905932737*fUpwind[9]*alphaDragSurf[19]+0.28284271247461906*alphaDragSurf[17]*fUpwind[18]+0.3162277660168379*alphaDragSurf[6]*fUpwind[18]+0.28284271247461906*fUpwind[17]*alphaDragSurf[18]+0.3162277660168379*fUpwind[6]*alphaDragSurf[18]+0.3162277660168379*alphaDragSurf[5]*fUpwind[17]+0.3162277660168379*fUpwind[5]*alphaDragSurf[17]+0.3535533905932737*alphaDragSurf[15]*fUpwind[16]+0.3535533905932737*fUpwind[15]*alphaDragSurf[16]+0.31622776601683794*alphaDragSurf[10]*fUpwind[14]+0.31622776601683794*fUpwind[10]*alphaDragSurf[14]+0.31622776601683794*alphaDragSurf[10]*fUpwind[13]+0.31622776601683794*fUpwind[10]*alphaDragSurf[13]+0.28284271247461906*alphaDragSurf[11]*fUpwind[12]+0.31622776601683794*alphaDragSurf[2]*fUpwind[12]+0.28284271247461906*fUpwind[11]*alphaDragSurf[12]+0.31622776601683794*fUpwind[2]*alphaDragSurf[12]+0.31622776601683794*alphaDragSurf[1]*fUpwind[11]+0.31622776601683794*fUpwind[1]*alphaDragSurf[11]+0.3535533905932737*alphaDragSurf[3]*fUpwind[10]+0.3535533905932737*fUpwind[3]*alphaDragSurf[10]+0.3162277660168379*alphaDragSurf[4]*fUpwind[8]+0.3162277660168379*fUpwind[4]*alphaDragSurf[8]+0.3162277660168379*alphaDragSurf[4]*fUpwind[7]+0.3162277660168379*fUpwind[4]*alphaDragSurf[7]+0.3535533905932737*alphaDragSurf[5]*fUpwind[6]+0.3535533905932737*fUpwind[5]*alphaDragSurf[6]+0.3535533905932737*alphaDragSurf[0]*fUpwind[4]+0.3535533905932737*fUpwind[0]*alphaDragSurf[4]+0.3535533905932737*alphaDragSurf[1]*fUpwind[2]+0.3535533905932737*fUpwind[1]*alphaDragSurf[2]; 
  Ghat[5] = 0.28284271247461906*alphaDragSurf[17]*fUpwind[19]+0.3162277660168379*alphaDragSurf[6]*fUpwind[19]+0.28284271247461906*fUpwind[17]*alphaDragSurf[19]+0.3162277660168379*fUpwind[6]*alphaDragSurf[19]+0.3535533905932737*alphaDragSurf[8]*fUpwind[18]+0.3535533905932737*fUpwind[8]*alphaDragSurf[18]+0.3162277660168379*alphaDragSurf[4]*fUpwind[17]+0.3162277660168379*fUpwind[4]*alphaDragSurf[17]+0.31622776601683794*alphaDragSurf[10]*fUpwind[16]+0.31622776601683794*fUpwind[10]*alphaDragSurf[16]+0.28284271247461906*alphaDragSurf[13]*fUpwind[15]+0.31622776601683794*alphaDragSurf[3]*fUpwind[15]+0.28284271247461906*fUpwind[13]*alphaDragSurf[15]+0.31622776601683794*fUpwind[3]*alphaDragSurf[15]+0.3535533905932737*alphaDragSurf[12]*fUpwind[14]+0.3535533905932737*fUpwind[12]*alphaDragSurf[14]+0.31622776601683794*alphaDragSurf[1]*fUpwind[13]+0.31622776601683794*fUpwind[1]*alphaDragSurf[13]+0.31622776601683794*alphaDragSurf[10]*fUpwind[11]+0.31622776601683794*fUpwind[10]*alphaDragSurf[11]+0.3535533905932737*alphaDragSurf[2]*fUpwind[10]+0.3535533905932737*fUpwind[2]*alphaDragSurf[10]+0.3162277660168379*alphaDragSurf[5]*fUpwind[9]+0.3162277660168379*fUpwind[5]*alphaDragSurf[9]+0.3162277660168379*alphaDragSurf[5]*fUpwind[7]+0.3162277660168379*fUpwind[5]*alphaDragSurf[7]+0.3535533905932737*alphaDragSurf[4]*fUpwind[6]+0.3535533905932737*fUpwind[4]*alphaDragSurf[6]+0.3535533905932737*alphaDragSurf[0]*fUpwind[5]+0.3535533905932737*fUpwind[0]*alphaDragSurf[5]+0.3535533905932737*alphaDragSurf[1]*fUpwind[3]+0.3535533905932737*fUpwind[1]*alphaDragSurf[3]; 
  Ghat[6] = 0.28284271247461906*alphaDragSurf[18]*fUpwind[19]+0.3162277660168379*alphaDragSurf[5]*fUpwind[19]+0.28284271247461906*fUpwind[18]*alphaDragSurf[19]+0.3162277660168379*fUpwind[5]*alphaDragSurf[19]+0.3162277660168379*alphaDragSurf[4]*fUpwind[18]+0.3162277660168379*fUpwind[4]*alphaDragSurf[18]+0.3535533905932737*alphaDragSurf[7]*fUpwind[17]+0.3535533905932737*fUpwind[7]*alphaDragSurf[17]+0.28284271247461906*alphaDragSurf[14]*fUpwind[16]+0.31622776601683794*alphaDragSurf[3]*fUpwind[16]+0.28284271247461906*fUpwind[14]*alphaDragSurf[16]+0.31622776601683794*fUpwind[3]*alphaDragSurf[16]+0.31622776601683794*alphaDragSurf[10]*fUpwind[15]+0.31622776601683794*fUpwind[10]*alphaDragSurf[15]+0.31622776601683794*alphaDragSurf[2]*fUpwind[14]+0.31622776601683794*fUpwind[2]*alphaDragSurf[14]+0.3535533905932737*alphaDragSurf[11]*fUpwind[13]+0.3535533905932737*fUpwind[11]*alphaDragSurf[13]+0.31622776601683794*alphaDragSurf[10]*fUpwind[12]+0.31622776601683794*fUpwind[10]*alphaDragSurf[12]+0.3535533905932737*alphaDragSurf[1]*fUpwind[10]+0.3535533905932737*fUpwind[1]*alphaDragSurf[10]+0.3162277660168379*alphaDragSurf[6]*fUpwind[9]+0.3162277660168379*fUpwind[6]*alphaDragSurf[9]+0.3162277660168379*alphaDragSurf[6]*fUpwind[8]+0.3162277660168379*fUpwind[6]*alphaDragSurf[8]+0.3535533905932737*alphaDragSurf[0]*fUpwind[6]+0.3535533905932737*fUpwind[0]*alphaDragSurf[6]+0.3535533905932737*alphaDragSurf[4]*fUpwind[5]+0.3535533905932737*fUpwind[4]*alphaDragSurf[5]+0.3535533905932737*alphaDragSurf[2]*fUpwind[3]+0.3535533905932737*fUpwind[2]*alphaDragSurf[3]; 
  Ghat[7] = 0.3162277660168379*alphaDragSurf[19]*fUpwind[19]+0.3162277660168379*alphaDragSurf[18]*fUpwind[18]+0.22587697572631277*alphaDragSurf[17]*fUpwind[17]+0.3535533905932737*alphaDragSurf[6]*fUpwind[17]+0.3535533905932737*fUpwind[6]*alphaDragSurf[17]+0.3162277660168379*alphaDragSurf[15]*fUpwind[15]+0.22587697572631277*alphaDragSurf[13]*fUpwind[13]+0.3535533905932737*alphaDragSurf[3]*fUpwind[13]+0.3535533905932737*fUpwind[3]*alphaDragSurf[13]+0.3162277660168379*alphaDragSurf[12]*fUpwind[12]+0.22587697572631277*alphaDragSurf[11]*fUpwind[11]+0.3535533905932737*alphaDragSurf[2]*fUpwind[11]+0.3535533905932737*fUpwind[2]*alphaDragSurf[11]+0.3162277660168379*alphaDragSurf[10]*fUpwind[10]+0.22587697572631277*alphaDragSurf[7]*fUpwind[7]+0.3535533905932737*alphaDragSurf[0]*fUpwind[7]+0.3535533905932737*fUpwind[0]*alphaDragSurf[7]+0.3162277660168379*alphaDragSurf[5]*fUpwind[5]+0.3162277660168379*alphaDragSurf[4]*fUpwind[4]+0.3162277660168379*alphaDragSurf[1]*fUpwind[1]; 
  Ghat[8] = 0.3162277660168379*alphaDragSurf[19]*fUpwind[19]+0.22587697572631277*alphaDragSurf[18]*fUpwind[18]+0.3535533905932737*alphaDragSurf[5]*fUpwind[18]+0.3535533905932737*fUpwind[5]*alphaDragSurf[18]+0.3162277660168379*alphaDragSurf[17]*fUpwind[17]+0.3162277660168379*alphaDragSurf[16]*fUpwind[16]+0.22587697572631277*alphaDragSurf[14]*fUpwind[14]+0.3535533905932737*alphaDragSurf[3]*fUpwind[14]+0.3535533905932737*fUpwind[3]*alphaDragSurf[14]+0.22587697572631277*alphaDragSurf[12]*fUpwind[12]+0.3535533905932737*alphaDragSurf[1]*fUpwind[12]+0.3535533905932737*fUpwind[1]*alphaDragSurf[12]+0.3162277660168379*alphaDragSurf[11]*fUpwind[11]+0.3162277660168379*alphaDragSurf[10]*fUpwind[10]+0.22587697572631277*alphaDragSurf[8]*fUpwind[8]+0.3535533905932737*alphaDragSurf[0]*fUpwind[8]+0.3535533905932737*fUpwind[0]*alphaDragSurf[8]+0.3162277660168379*alphaDragSurf[6]*fUpwind[6]+0.3162277660168379*alphaDragSurf[4]*fUpwind[4]+0.3162277660168379*alphaDragSurf[2]*fUpwind[2]; 
  Ghat[9] = 0.22587697572631277*alphaDragSurf[19]*fUpwind[19]+0.3535533905932737*alphaDragSurf[4]*fUpwind[19]+0.3535533905932737*fUpwind[4]*alphaDragSurf[19]+0.3162277660168379*alphaDragSurf[18]*fUpwind[18]+0.3162277660168379*alphaDragSurf[17]*fUpwind[17]+0.22587697572631277*alphaDragSurf[16]*fUpwind[16]+0.3535533905932737*alphaDragSurf[2]*fUpwind[16]+0.3535533905932737*fUpwind[2]*alphaDragSurf[16]+0.22587697572631277*alphaDragSurf[15]*fUpwind[15]+0.3535533905932737*alphaDragSurf[1]*fUpwind[15]+0.3535533905932737*fUpwind[1]*alphaDragSurf[15]+0.3162277660168379*alphaDragSurf[14]*fUpwind[14]+0.3162277660168379*alphaDragSurf[13]*fUpwind[13]+0.3162277660168379*alphaDragSurf[10]*fUpwind[10]+0.22587697572631277*alphaDragSurf[9]*fUpwind[9]+0.3535533905932737*alphaDragSurf[0]*fUpwind[9]+0.3535533905932737*fUpwind[0]*alphaDragSurf[9]+0.3162277660168379*alphaDragSurf[6]*fUpwind[6]+0.3162277660168379*alphaDragSurf[5]*fUpwind[5]+0.3162277660168379*alphaDragSurf[3]*fUpwind[3]; 
  Ghat[10] = 0.282842712474619*alphaDragSurf[14]*fUpwind[19]+0.282842712474619*alphaDragSurf[13]*fUpwind[19]+0.3162277660168379*alphaDragSurf[3]*fUpwind[19]+0.282842712474619*fUpwind[14]*alphaDragSurf[19]+0.282842712474619*fUpwind[13]*alphaDragSurf[19]+0.3162277660168379*fUpwind[3]*alphaDragSurf[19]+0.282842712474619*alphaDragSurf[16]*fUpwind[18]+0.282842712474619*alphaDragSurf[11]*fUpwind[18]+0.3162277660168379*alphaDragSurf[2]*fUpwind[18]+0.282842712474619*fUpwind[16]*alphaDragSurf[18]+0.282842712474619*fUpwind[11]*alphaDragSurf[18]+0.3162277660168379*fUpwind[2]*alphaDragSurf[18]+0.282842712474619*alphaDragSurf[15]*fUpwind[17]+0.282842712474619*alphaDragSurf[12]*fUpwind[17]+0.3162277660168379*alphaDragSurf[1]*fUpwind[17]+0.282842712474619*fUpwind[15]*alphaDragSurf[17]+0.282842712474619*fUpwind[12]*alphaDragSurf[17]+0.3162277660168379*fUpwind[1]*alphaDragSurf[17]+0.31622776601683794*alphaDragSurf[5]*fUpwind[16]+0.31622776601683794*fUpwind[5]*alphaDragSurf[16]+0.31622776601683794*alphaDragSurf[6]*fUpwind[15]+0.31622776601683794*fUpwind[6]*alphaDragSurf[15]+0.31622776601683794*alphaDragSurf[4]*fUpwind[14]+0.31622776601683794*fUpwind[4]*alphaDragSurf[14]+0.31622776601683794*alphaDragSurf[4]*fUpwind[13]+0.31622776601683794*fUpwind[4]*alphaDragSurf[13]+0.31622776601683794*alphaDragSurf[6]*fUpwind[12]+0.31622776601683794*fUpwind[6]*alphaDragSurf[12]+0.31622776601683794*alphaDragSurf[5]*fUpwind[11]+0.31622776601683794*fUpwind[5]*alphaDragSurf[11]+0.3162277660168379*alphaDragSurf[9]*fUpwind[10]+0.3162277660168379*alphaDragSurf[8]*fUpwind[10]+0.3162277660168379*alphaDragSurf[7]*fUpwind[10]+0.3535533905932737*alphaDragSurf[0]*fUpwind[10]+0.3162277660168379*fUpwind[9]*alphaDragSurf[10]+0.3162277660168379*fUpwind[8]*alphaDragSurf[10]+0.3162277660168379*fUpwind[7]*alphaDragSurf[10]+0.3535533905932737*fUpwind[0]*alphaDragSurf[10]+0.3535533905932737*alphaDragSurf[1]*fUpwind[6]+0.3535533905932737*fUpwind[1]*alphaDragSurf[6]+0.3535533905932737*alphaDragSurf[2]*fUpwind[5]+0.3535533905932737*fUpwind[2]*alphaDragSurf[5]+0.3535533905932737*alphaDragSurf[3]*fUpwind[4]+0.3535533905932737*fUpwind[3]*alphaDragSurf[4]; 
  Ghat[11] = 0.3162277660168379*alphaDragSurf[15]*fUpwind[19]+0.3162277660168379*fUpwind[15]*alphaDragSurf[19]+0.282842712474619*alphaDragSurf[10]*fUpwind[18]+0.282842712474619*fUpwind[10]*alphaDragSurf[18]+0.3162277660168379*alphaDragSurf[14]*fUpwind[17]+0.22587697572631277*alphaDragSurf[13]*fUpwind[17]+0.3535533905932737*alphaDragSurf[3]*fUpwind[17]+0.3162277660168379*fUpwind[14]*alphaDragSurf[17]+0.22587697572631277*fUpwind[13]*alphaDragSurf[17]+0.3535533905932737*fUpwind[3]*alphaDragSurf[17]+0.3535533905932737*alphaDragSurf[6]*fUpwind[13]+0.3535533905932737*fUpwind[6]*alphaDragSurf[13]+0.28284271247461906*alphaDragSurf[4]*fUpwind[12]+0.28284271247461906*fUpwind[4]*alphaDragSurf[12]+0.3162277660168379*alphaDragSurf[8]*fUpwind[11]+0.22587697572631277*alphaDragSurf[7]*fUpwind[11]+0.3535533905932737*alphaDragSurf[0]*fUpwind[11]+0.3162277660168379*fUpwind[8]*alphaDragSurf[11]+0.22587697572631277*fUpwind[7]*alphaDragSurf[11]+0.3535533905932737*fUpwind[0]*alphaDragSurf[11]+0.31622776601683794*alphaDragSurf[5]*fUpwind[10]+0.31622776601683794*fUpwind[5]*alphaDragSurf[10]+0.3535533905932737*alphaDragSurf[2]*fUpwind[7]+0.3535533905932737*fUpwind[2]*alphaDragSurf[7]+0.31622776601683794*alphaDragSurf[1]*fUpwind[4]+0.31622776601683794*fUpwind[1]*alphaDragSurf[4]; 
  Ghat[12] = 0.3162277660168379*alphaDragSurf[16]*fUpwind[19]+0.3162277660168379*fUpwind[16]*alphaDragSurf[19]+0.22587697572631277*alphaDragSurf[14]*fUpwind[18]+0.3162277660168379*alphaDragSurf[13]*fUpwind[18]+0.3535533905932737*alphaDragSurf[3]*fUpwind[18]+0.22587697572631277*fUpwind[14]*alphaDragSurf[18]+0.3162277660168379*fUpwind[13]*alphaDragSurf[18]+0.3535533905932737*fUpwind[3]*alphaDragSurf[18]+0.282842712474619*alphaDragSurf[10]*fUpwind[17]+0.282842712474619*fUpwind[10]*alphaDragSurf[17]+0.3535533905932737*alphaDragSurf[5]*fUpwind[14]+0.3535533905932737*fUpwind[5]*alphaDragSurf[14]+0.22587697572631277*alphaDragSurf[8]*fUpwind[12]+0.3162277660168379*alphaDragSurf[7]*fUpwind[12]+0.3535533905932737*alphaDragSurf[0]*fUpwind[12]+0.22587697572631277*fUpwind[8]*alphaDragSurf[12]+0.3162277660168379*fUpwind[7]*alphaDragSurf[12]+0.3535533905932737*fUpwind[0]*alphaDragSurf[12]+0.28284271247461906*alphaDragSurf[4]*fUpwind[11]+0.28284271247461906*fUpwind[4]*alphaDragSurf[11]+0.31622776601683794*alphaDragSurf[6]*fUpwind[10]+0.31622776601683794*fUpwind[6]*alphaDragSurf[10]+0.3535533905932737*alphaDragSurf[1]*fUpwind[8]+0.3535533905932737*fUpwind[1]*alphaDragSurf[8]+0.31622776601683794*alphaDragSurf[2]*fUpwind[4]+0.31622776601683794*fUpwind[2]*alphaDragSurf[4]; 
  Ghat[13] = 0.282842712474619*alphaDragSurf[10]*fUpwind[19]+0.282842712474619*fUpwind[10]*alphaDragSurf[19]+0.3162277660168379*alphaDragSurf[12]*fUpwind[18]+0.3162277660168379*fUpwind[12]*alphaDragSurf[18]+0.3162277660168379*alphaDragSurf[16]*fUpwind[17]+0.22587697572631277*alphaDragSurf[11]*fUpwind[17]+0.3535533905932737*alphaDragSurf[2]*fUpwind[17]+0.3162277660168379*fUpwind[16]*alphaDragSurf[17]+0.22587697572631277*fUpwind[11]*alphaDragSurf[17]+0.3535533905932737*fUpwind[2]*alphaDragSurf[17]+0.28284271247461906*alphaDragSurf[5]*fUpwind[15]+0.28284271247461906*fUpwind[5]*alphaDragSurf[15]+0.3162277660168379*alphaDragSurf[9]*fUpwind[13]+0.22587697572631277*alphaDragSurf[7]*fUpwind[13]+0.3535533905932737*alphaDragSurf[0]*fUpwind[13]+0.3162277660168379*fUpwind[9]*alphaDragSurf[13]+0.22587697572631277*fUpwind[7]*alphaDragSurf[13]+0.3535533905932737*fUpwind[0]*alphaDragSurf[13]+0.3535533905932737*alphaDragSurf[6]*fUpwind[11]+0.3535533905932737*fUpwind[6]*alphaDragSurf[11]+0.31622776601683794*alphaDragSurf[4]*fUpwind[10]+0.31622776601683794*fUpwind[4]*alphaDragSurf[10]+0.3535533905932737*alphaDragSurf[3]*fUpwind[7]+0.3535533905932737*fUpwind[3]*alphaDragSurf[7]+0.31622776601683794*alphaDragSurf[1]*fUpwind[5]+0.31622776601683794*fUpwind[1]*alphaDragSurf[5]; 
  Ghat[14] = 0.282842712474619*alphaDragSurf[10]*fUpwind[19]+0.282842712474619*fUpwind[10]*alphaDragSurf[19]+0.3162277660168379*alphaDragSurf[15]*fUpwind[18]+0.22587697572631277*alphaDragSurf[12]*fUpwind[18]+0.3535533905932737*alphaDragSurf[1]*fUpwind[18]+0.3162277660168379*fUpwind[15]*alphaDragSurf[18]+0.22587697572631277*fUpwind[12]*alphaDragSurf[18]+0.3535533905932737*fUpwind[1]*alphaDragSurf[18]+0.3162277660168379*alphaDragSurf[11]*fUpwind[17]+0.3162277660168379*fUpwind[11]*alphaDragSurf[17]+0.28284271247461906*alphaDragSurf[6]*fUpwind[16]+0.28284271247461906*fUpwind[6]*alphaDragSurf[16]+0.3162277660168379*alphaDragSurf[9]*fUpwind[14]+0.22587697572631277*alphaDragSurf[8]*fUpwind[14]+0.3535533905932737*alphaDragSurf[0]*fUpwind[14]+0.3162277660168379*fUpwind[9]*alphaDragSurf[14]+0.22587697572631277*fUpwind[8]*alphaDragSurf[14]+0.3535533905932737*fUpwind[0]*alphaDragSurf[14]+0.3535533905932737*alphaDragSurf[5]*fUpwind[12]+0.3535533905932737*fUpwind[5]*alphaDragSurf[12]+0.31622776601683794*alphaDragSurf[4]*fUpwind[10]+0.31622776601683794*fUpwind[4]*alphaDragSurf[10]+0.3535533905932737*alphaDragSurf[3]*fUpwind[8]+0.3535533905932737*fUpwind[3]*alphaDragSurf[8]+0.31622776601683794*alphaDragSurf[2]*fUpwind[6]+0.31622776601683794*fUpwind[2]*alphaDragSurf[6]; 
  Ghat[15] = 0.22587697572631277*alphaDragSurf[16]*fUpwind[19]+0.3162277660168379*alphaDragSurf[11]*fUpwind[19]+0.3535533905932737*alphaDragSurf[2]*fUpwind[19]+0.22587697572631277*fUpwind[16]*alphaDragSurf[19]+0.3162277660168379*fUpwind[11]*alphaDragSurf[19]+0.3535533905932737*fUpwind[2]*alphaDragSurf[19]+0.3162277660168379*alphaDragSurf[14]*fUpwind[18]+0.3162277660168379*fUpwind[14]*alphaDragSurf[18]+0.282842712474619*alphaDragSurf[10]*fUpwind[17]+0.282842712474619*fUpwind[10]*alphaDragSurf[17]+0.3535533905932737*alphaDragSurf[4]*fUpwind[16]+0.3535533905932737*fUpwind[4]*alphaDragSurf[16]+0.22587697572631277*alphaDragSurf[9]*fUpwind[15]+0.3162277660168379*alphaDragSurf[7]*fUpwind[15]+0.3535533905932737*alphaDragSurf[0]*fUpwind[15]+0.22587697572631277*fUpwind[9]*alphaDragSurf[15]+0.3162277660168379*fUpwind[7]*alphaDragSurf[15]+0.3535533905932737*fUpwind[0]*alphaDragSurf[15]+0.28284271247461906*alphaDragSurf[5]*fUpwind[13]+0.28284271247461906*fUpwind[5]*alphaDragSurf[13]+0.31622776601683794*alphaDragSurf[6]*fUpwind[10]+0.31622776601683794*fUpwind[6]*alphaDragSurf[10]+0.3535533905932737*alphaDragSurf[1]*fUpwind[9]+0.3535533905932737*fUpwind[1]*alphaDragSurf[9]+0.31622776601683794*alphaDragSurf[3]*fUpwind[5]+0.31622776601683794*fUpwind[3]*alphaDragSurf[5]; 
  Ghat[16] = 0.22587697572631277*alphaDragSurf[15]*fUpwind[19]+0.3162277660168379*alphaDragSurf[12]*fUpwind[19]+0.3535533905932737*alphaDragSurf[1]*fUpwind[19]+0.22587697572631277*fUpwind[15]*alphaDragSurf[19]+0.3162277660168379*fUpwind[12]*alphaDragSurf[19]+0.3535533905932737*fUpwind[1]*alphaDragSurf[19]+0.282842712474619*alphaDragSurf[10]*fUpwind[18]+0.282842712474619*fUpwind[10]*alphaDragSurf[18]+0.3162277660168379*alphaDragSurf[13]*fUpwind[17]+0.3162277660168379*fUpwind[13]*alphaDragSurf[17]+0.22587697572631277*alphaDragSurf[9]*fUpwind[16]+0.3162277660168379*alphaDragSurf[8]*fUpwind[16]+0.3535533905932737*alphaDragSurf[0]*fUpwind[16]+0.22587697572631277*fUpwind[9]*alphaDragSurf[16]+0.3162277660168379*fUpwind[8]*alphaDragSurf[16]+0.3535533905932737*fUpwind[0]*alphaDragSurf[16]+0.3535533905932737*alphaDragSurf[4]*fUpwind[15]+0.3535533905932737*fUpwind[4]*alphaDragSurf[15]+0.28284271247461906*alphaDragSurf[6]*fUpwind[14]+0.28284271247461906*fUpwind[6]*alphaDragSurf[14]+0.31622776601683794*alphaDragSurf[5]*fUpwind[10]+0.31622776601683794*fUpwind[5]*alphaDragSurf[10]+0.3535533905932737*alphaDragSurf[2]*fUpwind[9]+0.3535533905932737*fUpwind[2]*alphaDragSurf[9]+0.31622776601683794*alphaDragSurf[3]*fUpwind[6]+0.31622776601683794*fUpwind[3]*alphaDragSurf[6]; 
  Ghat[17] = 0.2529822128134704*alphaDragSurf[18]*fUpwind[19]+0.28284271247461906*alphaDragSurf[5]*fUpwind[19]+0.2529822128134704*fUpwind[18]*alphaDragSurf[19]+0.28284271247461906*fUpwind[5]*alphaDragSurf[19]+0.28284271247461906*alphaDragSurf[4]*fUpwind[18]+0.28284271247461906*fUpwind[4]*alphaDragSurf[18]+0.3162277660168379*alphaDragSurf[9]*fUpwind[17]+0.3162277660168379*alphaDragSurf[8]*fUpwind[17]+0.22587697572631277*alphaDragSurf[7]*fUpwind[17]+0.3535533905932737*alphaDragSurf[0]*fUpwind[17]+0.3162277660168379*fUpwind[9]*alphaDragSurf[17]+0.3162277660168379*fUpwind[8]*alphaDragSurf[17]+0.22587697572631277*fUpwind[7]*alphaDragSurf[17]+0.3535533905932737*fUpwind[0]*alphaDragSurf[17]+0.3162277660168379*alphaDragSurf[13]*fUpwind[16]+0.3162277660168379*fUpwind[13]*alphaDragSurf[16]+0.282842712474619*alphaDragSurf[10]*fUpwind[15]+0.282842712474619*fUpwind[10]*alphaDragSurf[15]+0.3162277660168379*alphaDragSurf[11]*fUpwind[14]+0.3162277660168379*fUpwind[11]*alphaDragSurf[14]+0.22587697572631277*alphaDragSurf[11]*fUpwind[13]+0.3535533905932737*alphaDragSurf[2]*fUpwind[13]+0.22587697572631277*fUpwind[11]*alphaDragSurf[13]+0.3535533905932737*fUpwind[2]*alphaDragSurf[13]+0.282842712474619*alphaDragSurf[10]*fUpwind[12]+0.282842712474619*fUpwind[10]*alphaDragSurf[12]+0.3535533905932737*alphaDragSurf[3]*fUpwind[11]+0.3535533905932737*fUpwind[3]*alphaDragSurf[11]+0.3162277660168379*alphaDragSurf[1]*fUpwind[10]+0.3162277660168379*fUpwind[1]*alphaDragSurf[10]+0.3535533905932737*alphaDragSurf[6]*fUpwind[7]+0.3535533905932737*fUpwind[6]*alphaDragSurf[7]+0.3162277660168379*alphaDragSurf[4]*fUpwind[5]+0.3162277660168379*fUpwind[4]*alphaDragSurf[5]; 
  Ghat[18] = 0.2529822128134704*alphaDragSurf[17]*fUpwind[19]+0.28284271247461906*alphaDragSurf[6]*fUpwind[19]+0.2529822128134704*fUpwind[17]*alphaDragSurf[19]+0.28284271247461906*fUpwind[6]*alphaDragSurf[19]+0.3162277660168379*alphaDragSurf[9]*fUpwind[18]+0.22587697572631277*alphaDragSurf[8]*fUpwind[18]+0.3162277660168379*alphaDragSurf[7]*fUpwind[18]+0.3535533905932737*alphaDragSurf[0]*fUpwind[18]+0.3162277660168379*fUpwind[9]*alphaDragSurf[18]+0.22587697572631277*fUpwind[8]*alphaDragSurf[18]+0.3162277660168379*fUpwind[7]*alphaDragSurf[18]+0.3535533905932737*fUpwind[0]*alphaDragSurf[18]+0.28284271247461906*alphaDragSurf[4]*fUpwind[17]+0.28284271247461906*fUpwind[4]*alphaDragSurf[17]+0.282842712474619*alphaDragSurf[10]*fUpwind[16]+0.282842712474619*fUpwind[10]*alphaDragSurf[16]+0.3162277660168379*alphaDragSurf[14]*fUpwind[15]+0.3162277660168379*fUpwind[14]*alphaDragSurf[15]+0.22587697572631277*alphaDragSurf[12]*fUpwind[14]+0.3535533905932737*alphaDragSurf[1]*fUpwind[14]+0.22587697572631277*fUpwind[12]*alphaDragSurf[14]+0.3535533905932737*fUpwind[1]*alphaDragSurf[14]+0.3162277660168379*alphaDragSurf[12]*fUpwind[13]+0.3162277660168379*fUpwind[12]*alphaDragSurf[13]+0.3535533905932737*alphaDragSurf[3]*fUpwind[12]+0.3535533905932737*fUpwind[3]*alphaDragSurf[12]+0.282842712474619*alphaDragSurf[10]*fUpwind[11]+0.282842712474619*fUpwind[10]*alphaDragSurf[11]+0.3162277660168379*alphaDragSurf[2]*fUpwind[10]+0.3162277660168379*fUpwind[2]*alphaDragSurf[10]+0.3535533905932737*alphaDragSurf[5]*fUpwind[8]+0.3535533905932737*fUpwind[5]*alphaDragSurf[8]+0.3162277660168379*alphaDragSurf[4]*fUpwind[6]+0.3162277660168379*fUpwind[4]*alphaDragSurf[6]; 
  Ghat[19] = 0.22587697572631277*alphaDragSurf[9]*fUpwind[19]+0.3162277660168379*alphaDragSurf[8]*fUpwind[19]+0.3162277660168379*alphaDragSurf[7]*fUpwind[19]+0.3535533905932737*alphaDragSurf[0]*fUpwind[19]+0.22587697572631277*fUpwind[9]*alphaDragSurf[19]+0.3162277660168379*fUpwind[8]*alphaDragSurf[19]+0.3162277660168379*fUpwind[7]*alphaDragSurf[19]+0.3535533905932737*fUpwind[0]*alphaDragSurf[19]+0.2529822128134704*alphaDragSurf[17]*fUpwind[18]+0.28284271247461906*alphaDragSurf[6]*fUpwind[18]+0.2529822128134704*fUpwind[17]*alphaDragSurf[18]+0.28284271247461906*fUpwind[6]*alphaDragSurf[18]+0.28284271247461906*alphaDragSurf[5]*fUpwind[17]+0.28284271247461906*fUpwind[5]*alphaDragSurf[17]+0.22587697572631277*alphaDragSurf[15]*fUpwind[16]+0.3162277660168379*alphaDragSurf[12]*fUpwind[16]+0.3535533905932737*alphaDragSurf[1]*fUpwind[16]+0.22587697572631277*fUpwind[15]*alphaDragSurf[16]+0.3162277660168379*fUpwind[12]*alphaDragSurf[16]+0.3535533905932737*fUpwind[1]*alphaDragSurf[16]+0.3162277660168379*alphaDragSurf[11]*fUpwind[15]+0.3535533905932737*alphaDragSurf[2]*fUpwind[15]+0.3162277660168379*fUpwind[11]*alphaDragSurf[15]+0.3535533905932737*fUpwind[2]*alphaDragSurf[15]+0.282842712474619*alphaDragSurf[10]*fUpwind[14]+0.282842712474619*fUpwind[10]*alphaDragSurf[14]+0.282842712474619*alphaDragSurf[10]*fUpwind[13]+0.282842712474619*fUpwind[10]*alphaDragSurf[13]+0.3162277660168379*alphaDragSurf[3]*fUpwind[10]+0.3162277660168379*fUpwind[3]*alphaDragSurf[10]+0.3535533905932737*alphaDragSurf[4]*fUpwind[9]+0.3535533905932737*fUpwind[4]*alphaDragSurf[9]+0.3162277660168379*alphaDragSurf[5]*fUpwind[6]+0.3162277660168379*fUpwind[5]*alphaDragSurf[6]; 

  out[0] += -(0.35355339059327373*Ghat[0]*dv1); 
  out[1] += -(0.35355339059327373*Ghat[1]*dv1); 
  out[2] += 0.6123724356957945*Ghat[0]*dv1; 
  out[3] += -(0.35355339059327373*Ghat[2]*dv1); 
  out[4] += -(0.35355339059327373*Ghat[3]*dv1); 
  out[5] += 0.6123724356957945*Ghat[1]*dv1; 
  out[6] += -(0.35355339059327373*Ghat[4]*dv1); 
  out[7] += 0.6123724356957945*Ghat[2]*dv1; 
  out[8] += -(0.35355339059327373*Ghat[5]*dv1); 
  out[9] += 0.6123724356957945*Ghat[3]*dv1; 
  out[10] += -(0.35355339059327373*Ghat[6]*dv1); 
  out[11] += -(0.35355339059327373*Ghat[7]*dv1); 
  out[12] += -(0.7905694150420948*Ghat[0]*dv1); 
  out[13] += -(0.35355339059327373*Ghat[8]*dv1); 
  out[14] += -(0.35355339059327373*Ghat[9]*dv1); 
  out[15] += 0.6123724356957945*Ghat[4]*dv1; 
  out[16] += 0.6123724356957945*Ghat[5]*dv1; 
  out[17] += -(0.35355339059327373*Ghat[10]*dv1); 
  out[18] += 0.6123724356957945*Ghat[6]*dv1; 
  out[19] += 0.6123724356957945*Ghat[7]*dv1; 
  out[20] += -(0.7905694150420949*Ghat[1]*dv1); 
  out[21] += -(0.35355339059327373*Ghat[11]*dv1); 
  out[22] += -(0.7905694150420949*Ghat[2]*dv1); 
  out[23] += -(0.35355339059327373*Ghat[12]*dv1); 
  out[24] += 0.6123724356957945*Ghat[8]*dv1; 
  out[25] += -(0.35355339059327373*Ghat[13]*dv1); 
  out[26] += -(0.7905694150420949*Ghat[3]*dv1); 
  out[27] += -(0.35355339059327373*Ghat[14]*dv1); 
  out[28] += -(0.35355339059327373*Ghat[15]*dv1); 
  out[29] += 0.6123724356957945*Ghat[9]*dv1; 
  out[30] += -(0.35355339059327373*Ghat[16]*dv1); 
  out[31] += 0.6123724356957945*Ghat[10]*dv1; 
  out[32] += 0.6123724356957945*Ghat[11]*dv1; 
  out[33] += -(0.7905694150420948*Ghat[4]*dv1); 
  out[34] += 0.6123724356957945*Ghat[12]*dv1; 
  out[35] += 0.6123724356957945*Ghat[13]*dv1; 
  out[36] += -(0.7905694150420948*Ghat[5]*dv1); 
  out[37] += -(0.35355339059327373*Ghat[17]*dv1); 
  out[38] += -(0.7905694150420948*Ghat[6]*dv1); 
  out[39] += -(0.35355339059327373*Ghat[18]*dv1); 
  out[40] += 0.6123724356957945*Ghat[14]*dv1; 
  out[41] += 0.6123724356957945*Ghat[15]*dv1; 
  out[42] += -(0.35355339059327373*Ghat[19]*dv1); 
  out[43] += 0.6123724356957945*Ghat[16]*dv1; 
  out[44] += 0.6123724356957945*Ghat[17]*dv1; 
  out[45] += -(0.7905694150420949*Ghat[10]*dv1); 
  out[46] += 0.6123724356957945*Ghat[18]*dv1; 
  out[47] += 0.6123724356957945*Ghat[19]*dv1; 
} 
