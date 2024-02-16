#include <gkyl_fpo_vlasov_kernels.h> 
#include <gkyl_basis_ser_4x_p2_surfx4_eval_quad.h> 
#include <gkyl_basis_ser_4x_p2_upwind_quad_to_modal.h> 

 
GKYL_CU_DH void fpo_vlasov_drag_boundary_surfvz_1x3v_ser_p2_lovz(const double *dxv, const double* drag_coeff_stencil[9], const double* f_stencil[9], double* out) { 
  // dxv[NDIM]: Cell spacing in each direction. 
  // drag_coeff_stencil[9]: 9-cell stencil of drag coefficient. 
  // f_stencil[9]: 9-cell stencil of distribution function. 
  // out: Incremented output. 


  double dv1 = 2.0/dxv[3]; 
 
  const double* aC = &drag_coeff_stencil[0][96]; 
  const double* aR = &drag_coeff_stencil[1][96]; 
  const double* fC = f_stencil[0]; 
  const double* fR = f_stencil[1]; 

  double fUpwindQuad[27] = {0.0};
  double fUpwind[20] = {0.0};
  double alphaDragSurf[20] = {0.0}; 
  double Ghat[20] = {0.0}; 

  alphaDragSurf[0] = 0.34587411908091625*aR[14]+0.34587411908091625*aC[14]-0.49755260400283263*aR[4]+0.49755260400283263*aC[4]+0.3535533905932737*aR[0]+0.3535533905932737*aC[0]; 
  alphaDragSurf[1] = 0.34587411908091625*aR[28]+0.34587411908091625*aC[28]-0.49755260400283263*aR[8]+0.49755260400283263*aC[8]+0.3535533905932737*aR[1]+0.3535533905932737*aC[1]; 
  alphaDragSurf[2] = 0.34587411908091625*aR[29]+0.34587411908091625*aC[29]-0.49755260400283263*aR[9]+0.49755260400283263*aC[9]+0.3535533905932737*aR[2]+0.3535533905932737*aC[2]; 
  alphaDragSurf[3] = 0.34587411908091625*aR[30]+0.34587411908091625*aC[30]-0.49755260400283263*aR[10]+0.49755260400283263*aC[10]+0.3535533905932737*aR[3]+0.3535533905932737*aC[3]; 
  alphaDragSurf[4] = 0.34587411908091625*aR[41]+0.34587411908091625*aC[41]-0.49755260400283263*aR[16]+0.49755260400283263*aC[16]+0.3535533905932737*aR[5]+0.3535533905932737*aC[5]; 
  alphaDragSurf[5] = 0.34587411908091625*aR[42]+0.34587411908091625*aC[42]-0.49755260400283263*aR[17]+0.49755260400283263*aC[17]+0.3535533905932737*aR[6]+0.3535533905932737*aC[6]; 
  alphaDragSurf[6] = 0.34587411908091625*aR[43]+0.34587411908091625*aC[43]-0.49755260400283263*aR[18]+0.49755260400283263*aC[18]+0.3535533905932737*aR[7]+0.3535533905932737*aC[7]; 
  alphaDragSurf[7] = -(0.49755260400283263*aR[25])+0.49755260400283263*aC[25]+0.3535533905932737*aR[11]+0.3535533905932737*aC[11]; 
  alphaDragSurf[8] = -(0.49755260400283263*aR[26])+0.49755260400283263*aC[26]+0.3535533905932737*aR[12]+0.3535533905932737*aC[12]; 
  alphaDragSurf[9] = -(0.49755260400283263*aR[27])+0.49755260400283263*aC[27]+0.3535533905932737*aR[13]+0.3535533905932737*aC[13]; 
  alphaDragSurf[10] = 0.34587411908091625*aR[47]+0.34587411908091625*aC[47]-0.49755260400283263*aR[31]+0.49755260400283263*aC[31]+0.3535533905932737*aR[15]+0.3535533905932737*aC[15]; 
  alphaDragSurf[11] = -(0.49755260400283263*aR[35])+0.49755260400283263*aC[35]+0.3535533905932737*aR[19]+0.3535533905932737*aC[19]; 
  alphaDragSurf[12] = -(0.49755260400283263*aR[36])+0.49755260400283263*aC[36]+0.3535533905932737*aR[20]+0.3535533905932737*aC[20]; 
  alphaDragSurf[13] = -(0.49755260400283263*aR[37])+0.49755260400283263*aC[37]+0.3535533905932737*aR[21]+0.3535533905932737*aC[21]; 
  alphaDragSurf[14] = -(0.49755260400283263*aR[38])+0.49755260400283263*aC[38]+0.3535533905932737*aR[22]+0.3535533905932737*aC[22]; 
  alphaDragSurf[15] = -(0.49755260400283263*aR[39])+0.49755260400283263*aC[39]+0.3535533905932737*aR[23]+0.3535533905932737*aC[23]; 
  alphaDragSurf[16] = -(0.49755260400283263*aR[40])+0.49755260400283263*aC[40]+0.3535533905932737*aR[24]+0.3535533905932737*aC[24]; 
  alphaDragSurf[17] = -(0.49755260400283263*aR[44])+0.49755260400283263*aC[44]+0.3535533905932737*aR[32]+0.3535533905932737*aC[32]; 
  alphaDragSurf[18] = -(0.49755260400283263*aR[45])+0.49755260400283263*aC[45]+0.3535533905932737*aR[33]+0.3535533905932737*aC[33]; 
  alphaDragSurf[19] = -(0.49755260400283263*aR[46])+0.49755260400283263*aC[46]+0.3535533905932737*aR[34]+0.3535533905932737*aC[34]; 

  if (0.5692099788303082*(alphaDragSurf[19]+alphaDragSurf[18]+alphaDragSurf[17])-0.42426406871192807*(alphaDragSurf[16]+alphaDragSurf[15])-0.42426406871192845*(alphaDragSurf[14]+alphaDragSurf[13])-0.42426406871192807*alphaDragSurf[12]-0.42426406871192845*alphaDragSurf[11]-0.853814968245462*alphaDragSurf[10]+0.3162277660168379*(alphaDragSurf[9]+alphaDragSurf[8]+alphaDragSurf[7])+0.6363961030678926*(alphaDragSurf[6]+alphaDragSurf[5]+alphaDragSurf[4])-0.4743416490252568*(alphaDragSurf[3]+alphaDragSurf[2]+alphaDragSurf[1])+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[0] = ser_4x_p2_surfx4_eval_quad_node_0_r(fC); 
  } else { 
    fUpwindQuad[0] = ser_4x_p2_surfx4_eval_quad_node_0_l(fR); 
  } 
  if (-(0.711512473537885*alphaDragSurf[19])+0.5303300858899104*(alphaDragSurf[16]+alphaDragSurf[15])-0.42426406871192807*alphaDragSurf[12]-0.42426406871192845*alphaDragSurf[11]-0.3952847075210473*alphaDragSurf[9]+0.3162277660168379*(alphaDragSurf[8]+alphaDragSurf[7])+0.6363961030678926*alphaDragSurf[4]-0.4743416490252568*(alphaDragSurf[2]+alphaDragSurf[1])+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[1] = ser_4x_p2_surfx4_eval_quad_node_1_r(fC); 
  } else { 
    fUpwindQuad[1] = ser_4x_p2_surfx4_eval_quad_node_1_l(fR); 
  } 
  if (0.5692099788303082*alphaDragSurf[19]-0.5692099788303082*(alphaDragSurf[18]+alphaDragSurf[17])-0.42426406871192807*(alphaDragSurf[16]+alphaDragSurf[15])+0.42426406871192845*(alphaDragSurf[14]+alphaDragSurf[13])-0.42426406871192807*alphaDragSurf[12]-0.42426406871192845*alphaDragSurf[11]+0.853814968245462*alphaDragSurf[10]+0.3162277660168379*(alphaDragSurf[9]+alphaDragSurf[8]+alphaDragSurf[7])-0.6363961030678926*(alphaDragSurf[6]+alphaDragSurf[5])+0.6363961030678926*alphaDragSurf[4]+0.4743416490252568*alphaDragSurf[3]-0.4743416490252568*(alphaDragSurf[2]+alphaDragSurf[1])+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[2] = ser_4x_p2_surfx4_eval_quad_node_2_r(fC); 
  } else { 
    fUpwindQuad[2] = ser_4x_p2_surfx4_eval_quad_node_2_l(fR); 
  } 
  if (-(0.711512473537885*alphaDragSurf[18])-0.42426406871192807*alphaDragSurf[15]+0.5303300858899104*alphaDragSurf[14]-0.42426406871192845*alphaDragSurf[13]+0.5303300858899104*alphaDragSurf[12]+0.3162277660168379*alphaDragSurf[9]-0.3952847075210473*alphaDragSurf[8]+0.3162277660168379*alphaDragSurf[7]+0.6363961030678926*alphaDragSurf[5]-0.4743416490252568*(alphaDragSurf[3]+alphaDragSurf[1])+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[3] = ser_4x_p2_surfx4_eval_quad_node_3_r(fC); 
  } else { 
    fUpwindQuad[3] = ser_4x_p2_surfx4_eval_quad_node_3_l(fR); 
  } 
  if (0.5303300858899104*(alphaDragSurf[15]+alphaDragSurf[12])-0.3952847075210473*(alphaDragSurf[9]+alphaDragSurf[8])+0.3162277660168379*alphaDragSurf[7]-0.4743416490252568*alphaDragSurf[1]+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[4] = ser_4x_p2_surfx4_eval_quad_node_4_r(fC); 
  } else { 
    fUpwindQuad[4] = ser_4x_p2_surfx4_eval_quad_node_4_l(fR); 
  } 
  if (0.711512473537885*alphaDragSurf[18]-0.42426406871192807*alphaDragSurf[15]-0.5303300858899104*alphaDragSurf[14]+0.42426406871192845*alphaDragSurf[13]+0.5303300858899104*alphaDragSurf[12]+0.3162277660168379*alphaDragSurf[9]-0.3952847075210473*alphaDragSurf[8]+0.3162277660168379*alphaDragSurf[7]-0.6363961030678926*alphaDragSurf[5]+0.4743416490252568*alphaDragSurf[3]-0.4743416490252568*alphaDragSurf[1]+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[5] = ser_4x_p2_surfx4_eval_quad_node_5_r(fC); 
  } else { 
    fUpwindQuad[5] = ser_4x_p2_surfx4_eval_quad_node_5_l(fR); 
  } 
  if (-(0.5692099788303082*alphaDragSurf[19])+0.5692099788303082*alphaDragSurf[18]-0.5692099788303082*alphaDragSurf[17]+0.42426406871192807*alphaDragSurf[16]-0.42426406871192807*alphaDragSurf[15]-0.42426406871192845*(alphaDragSurf[14]+alphaDragSurf[13])-0.42426406871192807*alphaDragSurf[12]+0.42426406871192845*alphaDragSurf[11]+0.853814968245462*alphaDragSurf[10]+0.3162277660168379*(alphaDragSurf[9]+alphaDragSurf[8]+alphaDragSurf[7])-0.6363961030678926*alphaDragSurf[6]+0.6363961030678926*alphaDragSurf[5]-0.6363961030678926*alphaDragSurf[4]-0.4743416490252568*alphaDragSurf[3]+0.4743416490252568*alphaDragSurf[2]-0.4743416490252568*alphaDragSurf[1]+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[6] = ser_4x_p2_surfx4_eval_quad_node_6_r(fC); 
  } else { 
    fUpwindQuad[6] = ser_4x_p2_surfx4_eval_quad_node_6_l(fR); 
  } 
  if (0.711512473537885*alphaDragSurf[19]-0.5303300858899104*alphaDragSurf[16]+0.5303300858899104*alphaDragSurf[15]-0.42426406871192807*alphaDragSurf[12]+0.42426406871192845*alphaDragSurf[11]-0.3952847075210473*alphaDragSurf[9]+0.3162277660168379*(alphaDragSurf[8]+alphaDragSurf[7])-0.6363961030678926*alphaDragSurf[4]+0.4743416490252568*alphaDragSurf[2]-0.4743416490252568*alphaDragSurf[1]+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[7] = ser_4x_p2_surfx4_eval_quad_node_7_r(fC); 
  } else { 
    fUpwindQuad[7] = ser_4x_p2_surfx4_eval_quad_node_7_l(fR); 
  } 
  if (-(0.5692099788303082*(alphaDragSurf[19]+alphaDragSurf[18]))+0.5692099788303082*alphaDragSurf[17]+0.42426406871192807*alphaDragSurf[16]-0.42426406871192807*alphaDragSurf[15]+0.42426406871192845*(alphaDragSurf[14]+alphaDragSurf[13])-0.42426406871192807*alphaDragSurf[12]+0.42426406871192845*alphaDragSurf[11]-0.853814968245462*alphaDragSurf[10]+0.3162277660168379*(alphaDragSurf[9]+alphaDragSurf[8]+alphaDragSurf[7])+0.6363961030678926*alphaDragSurf[6]-0.6363961030678926*(alphaDragSurf[5]+alphaDragSurf[4])+0.4743416490252568*(alphaDragSurf[3]+alphaDragSurf[2])-0.4743416490252568*alphaDragSurf[1]+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[8] = ser_4x_p2_surfx4_eval_quad_node_8_r(fC); 
  } else { 
    fUpwindQuad[8] = ser_4x_p2_surfx4_eval_quad_node_8_l(fR); 
  } 
  if (-(0.711512473537885*alphaDragSurf[17])-0.42426406871192807*alphaDragSurf[16]-0.42426406871192845*alphaDragSurf[14]+0.5303300858899104*(alphaDragSurf[13]+alphaDragSurf[11])+0.3162277660168379*(alphaDragSurf[9]+alphaDragSurf[8])-0.3952847075210473*alphaDragSurf[7]+0.6363961030678926*alphaDragSurf[6]-0.4743416490252568*(alphaDragSurf[3]+alphaDragSurf[2])+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[9] = ser_4x_p2_surfx4_eval_quad_node_9_r(fC); 
  } else { 
    fUpwindQuad[9] = ser_4x_p2_surfx4_eval_quad_node_9_l(fR); 
  } 
  if (0.5303300858899104*(alphaDragSurf[16]+alphaDragSurf[11])-0.3952847075210473*alphaDragSurf[9]+0.3162277660168379*alphaDragSurf[8]-0.3952847075210473*alphaDragSurf[7]-0.4743416490252568*alphaDragSurf[2]+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[10] = ser_4x_p2_surfx4_eval_quad_node_10_r(fC); 
  } else { 
    fUpwindQuad[10] = ser_4x_p2_surfx4_eval_quad_node_10_l(fR); 
  } 
  if (0.711512473537885*alphaDragSurf[17]-0.42426406871192807*alphaDragSurf[16]+0.42426406871192845*alphaDragSurf[14]-0.5303300858899104*alphaDragSurf[13]+0.5303300858899104*alphaDragSurf[11]+0.3162277660168379*(alphaDragSurf[9]+alphaDragSurf[8])-0.3952847075210473*alphaDragSurf[7]-0.6363961030678926*alphaDragSurf[6]+0.4743416490252568*alphaDragSurf[3]-0.4743416490252568*alphaDragSurf[2]+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[11] = ser_4x_p2_surfx4_eval_quad_node_11_r(fC); 
  } else { 
    fUpwindQuad[11] = ser_4x_p2_surfx4_eval_quad_node_11_l(fR); 
  } 
  if (0.5303300858899104*(alphaDragSurf[14]+alphaDragSurf[13])+0.3162277660168379*alphaDragSurf[9]-0.3952847075210473*(alphaDragSurf[8]+alphaDragSurf[7])-0.4743416490252568*alphaDragSurf[3]+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[12] = ser_4x_p2_surfx4_eval_quad_node_12_r(fC); 
  } else { 
    fUpwindQuad[12] = ser_4x_p2_surfx4_eval_quad_node_12_l(fR); 
  } 
  if (0.3535533905932737*alphaDragSurf[0]-0.3952847075210473*(alphaDragSurf[9]+alphaDragSurf[8]+alphaDragSurf[7]) > 0) { 
    fUpwindQuad[13] = ser_4x_p2_surfx4_eval_quad_node_13_r(fC); 
  } else { 
    fUpwindQuad[13] = ser_4x_p2_surfx4_eval_quad_node_13_l(fR); 
  } 
  if (-(0.5303300858899104*(alphaDragSurf[14]+alphaDragSurf[13]))+0.3162277660168379*alphaDragSurf[9]-0.3952847075210473*(alphaDragSurf[8]+alphaDragSurf[7])+0.4743416490252568*alphaDragSurf[3]+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[14] = ser_4x_p2_surfx4_eval_quad_node_14_r(fC); 
  } else { 
    fUpwindQuad[14] = ser_4x_p2_surfx4_eval_quad_node_14_l(fR); 
  } 
  if (0.711512473537885*alphaDragSurf[17]+0.42426406871192807*alphaDragSurf[16]-0.42426406871192845*alphaDragSurf[14]+0.5303300858899104*alphaDragSurf[13]-0.5303300858899104*alphaDragSurf[11]+0.3162277660168379*(alphaDragSurf[9]+alphaDragSurf[8])-0.3952847075210473*alphaDragSurf[7]-0.6363961030678926*alphaDragSurf[6]-0.4743416490252568*alphaDragSurf[3]+0.4743416490252568*alphaDragSurf[2]+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[15] = ser_4x_p2_surfx4_eval_quad_node_15_r(fC); 
  } else { 
    fUpwindQuad[15] = ser_4x_p2_surfx4_eval_quad_node_15_l(fR); 
  } 
  if (-(0.5303300858899104*(alphaDragSurf[16]+alphaDragSurf[11]))-0.3952847075210473*alphaDragSurf[9]+0.3162277660168379*alphaDragSurf[8]-0.3952847075210473*alphaDragSurf[7]+0.4743416490252568*alphaDragSurf[2]+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[16] = ser_4x_p2_surfx4_eval_quad_node_16_r(fC); 
  } else { 
    fUpwindQuad[16] = ser_4x_p2_surfx4_eval_quad_node_16_l(fR); 
  } 
  if (-(0.711512473537885*alphaDragSurf[17])+0.42426406871192807*alphaDragSurf[16]+0.42426406871192845*alphaDragSurf[14]-0.5303300858899104*(alphaDragSurf[13]+alphaDragSurf[11])+0.3162277660168379*(alphaDragSurf[9]+alphaDragSurf[8])-0.3952847075210473*alphaDragSurf[7]+0.6363961030678926*alphaDragSurf[6]+0.4743416490252568*(alphaDragSurf[3]+alphaDragSurf[2])+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[17] = ser_4x_p2_surfx4_eval_quad_node_17_r(fC); 
  } else { 
    fUpwindQuad[17] = ser_4x_p2_surfx4_eval_quad_node_17_l(fR); 
  } 
  if (-(0.5692099788303082*(alphaDragSurf[19]+alphaDragSurf[18]))+0.5692099788303082*alphaDragSurf[17]-0.42426406871192807*alphaDragSurf[16]+0.42426406871192807*alphaDragSurf[15]-0.42426406871192845*(alphaDragSurf[14]+alphaDragSurf[13])+0.42426406871192807*alphaDragSurf[12]-0.42426406871192845*alphaDragSurf[11]+0.853814968245462*alphaDragSurf[10]+0.3162277660168379*(alphaDragSurf[9]+alphaDragSurf[8]+alphaDragSurf[7])+0.6363961030678926*alphaDragSurf[6]-0.6363961030678926*(alphaDragSurf[5]+alphaDragSurf[4])-0.4743416490252568*(alphaDragSurf[3]+alphaDragSurf[2])+0.4743416490252568*alphaDragSurf[1]+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[18] = ser_4x_p2_surfx4_eval_quad_node_18_r(fC); 
  } else { 
    fUpwindQuad[18] = ser_4x_p2_surfx4_eval_quad_node_18_l(fR); 
  } 
  if (0.711512473537885*alphaDragSurf[19]+0.5303300858899104*alphaDragSurf[16]-0.5303300858899104*alphaDragSurf[15]+0.42426406871192807*alphaDragSurf[12]-0.42426406871192845*alphaDragSurf[11]-0.3952847075210473*alphaDragSurf[9]+0.3162277660168379*(alphaDragSurf[8]+alphaDragSurf[7])-0.6363961030678926*alphaDragSurf[4]-0.4743416490252568*alphaDragSurf[2]+0.4743416490252568*alphaDragSurf[1]+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[19] = ser_4x_p2_surfx4_eval_quad_node_19_r(fC); 
  } else { 
    fUpwindQuad[19] = ser_4x_p2_surfx4_eval_quad_node_19_l(fR); 
  } 
  if (-(0.5692099788303082*alphaDragSurf[19])+0.5692099788303082*alphaDragSurf[18]-0.5692099788303082*alphaDragSurf[17]-0.42426406871192807*alphaDragSurf[16]+0.42426406871192807*alphaDragSurf[15]+0.42426406871192845*(alphaDragSurf[14]+alphaDragSurf[13])+0.42426406871192807*alphaDragSurf[12]-0.42426406871192845*alphaDragSurf[11]-0.853814968245462*alphaDragSurf[10]+0.3162277660168379*(alphaDragSurf[9]+alphaDragSurf[8]+alphaDragSurf[7])-0.6363961030678926*alphaDragSurf[6]+0.6363961030678926*alphaDragSurf[5]-0.6363961030678926*alphaDragSurf[4]+0.4743416490252568*alphaDragSurf[3]-0.4743416490252568*alphaDragSurf[2]+0.4743416490252568*alphaDragSurf[1]+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[20] = ser_4x_p2_surfx4_eval_quad_node_20_r(fC); 
  } else { 
    fUpwindQuad[20] = ser_4x_p2_surfx4_eval_quad_node_20_l(fR); 
  } 
  if (0.711512473537885*alphaDragSurf[18]+0.42426406871192807*alphaDragSurf[15]+0.5303300858899104*alphaDragSurf[14]-0.42426406871192845*alphaDragSurf[13]-0.5303300858899104*alphaDragSurf[12]+0.3162277660168379*alphaDragSurf[9]-0.3952847075210473*alphaDragSurf[8]+0.3162277660168379*alphaDragSurf[7]-0.6363961030678926*alphaDragSurf[5]-0.4743416490252568*alphaDragSurf[3]+0.4743416490252568*alphaDragSurf[1]+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[21] = ser_4x_p2_surfx4_eval_quad_node_21_r(fC); 
  } else { 
    fUpwindQuad[21] = ser_4x_p2_surfx4_eval_quad_node_21_l(fR); 
  } 
  if (-(0.5303300858899104*(alphaDragSurf[15]+alphaDragSurf[12]))-0.3952847075210473*(alphaDragSurf[9]+alphaDragSurf[8])+0.3162277660168379*alphaDragSurf[7]+0.4743416490252568*alphaDragSurf[1]+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[22] = ser_4x_p2_surfx4_eval_quad_node_22_r(fC); 
  } else { 
    fUpwindQuad[22] = ser_4x_p2_surfx4_eval_quad_node_22_l(fR); 
  } 
  if (-(0.711512473537885*alphaDragSurf[18])+0.42426406871192807*alphaDragSurf[15]-0.5303300858899104*alphaDragSurf[14]+0.42426406871192845*alphaDragSurf[13]-0.5303300858899104*alphaDragSurf[12]+0.3162277660168379*alphaDragSurf[9]-0.3952847075210473*alphaDragSurf[8]+0.3162277660168379*alphaDragSurf[7]+0.6363961030678926*alphaDragSurf[5]+0.4743416490252568*(alphaDragSurf[3]+alphaDragSurf[1])+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[23] = ser_4x_p2_surfx4_eval_quad_node_23_r(fC); 
  } else { 
    fUpwindQuad[23] = ser_4x_p2_surfx4_eval_quad_node_23_l(fR); 
  } 
  if (0.5692099788303082*alphaDragSurf[19]-0.5692099788303082*(alphaDragSurf[18]+alphaDragSurf[17])+0.42426406871192807*(alphaDragSurf[16]+alphaDragSurf[15])-0.42426406871192845*(alphaDragSurf[14]+alphaDragSurf[13])+0.42426406871192807*alphaDragSurf[12]+0.42426406871192845*alphaDragSurf[11]-0.853814968245462*alphaDragSurf[10]+0.3162277660168379*(alphaDragSurf[9]+alphaDragSurf[8]+alphaDragSurf[7])-0.6363961030678926*(alphaDragSurf[6]+alphaDragSurf[5])+0.6363961030678926*alphaDragSurf[4]-0.4743416490252568*alphaDragSurf[3]+0.4743416490252568*(alphaDragSurf[2]+alphaDragSurf[1])+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[24] = ser_4x_p2_surfx4_eval_quad_node_24_r(fC); 
  } else { 
    fUpwindQuad[24] = ser_4x_p2_surfx4_eval_quad_node_24_l(fR); 
  } 
  if (-(0.711512473537885*alphaDragSurf[19])-0.5303300858899104*(alphaDragSurf[16]+alphaDragSurf[15])+0.42426406871192807*alphaDragSurf[12]+0.42426406871192845*alphaDragSurf[11]-0.3952847075210473*alphaDragSurf[9]+0.3162277660168379*(alphaDragSurf[8]+alphaDragSurf[7])+0.6363961030678926*alphaDragSurf[4]+0.4743416490252568*(alphaDragSurf[2]+alphaDragSurf[1])+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[25] = ser_4x_p2_surfx4_eval_quad_node_25_r(fC); 
  } else { 
    fUpwindQuad[25] = ser_4x_p2_surfx4_eval_quad_node_25_l(fR); 
  } 
  if (0.5692099788303082*(alphaDragSurf[19]+alphaDragSurf[18]+alphaDragSurf[17])+0.42426406871192807*(alphaDragSurf[16]+alphaDragSurf[15])+0.42426406871192845*(alphaDragSurf[14]+alphaDragSurf[13])+0.42426406871192807*alphaDragSurf[12]+0.42426406871192845*alphaDragSurf[11]+0.853814968245462*alphaDragSurf[10]+0.3162277660168379*(alphaDragSurf[9]+alphaDragSurf[8]+alphaDragSurf[7])+0.6363961030678926*(alphaDragSurf[6]+alphaDragSurf[5]+alphaDragSurf[4])+0.4743416490252568*(alphaDragSurf[3]+alphaDragSurf[2]+alphaDragSurf[1])+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[26] = ser_4x_p2_surfx4_eval_quad_node_26_r(fC); 
  } else { 
    fUpwindQuad[26] = ser_4x_p2_surfx4_eval_quad_node_26_l(fR); 
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

  out[0] += 0.35355339059327373*Ghat[0]*dv1; 
  out[1] += 0.35355339059327373*Ghat[1]*dv1; 
  out[2] += 0.35355339059327373*Ghat[2]*dv1; 
  out[3] += 0.35355339059327373*Ghat[3]*dv1; 
  out[4] += 0.6123724356957945*Ghat[0]*dv1; 
  out[5] += 0.35355339059327373*Ghat[4]*dv1; 
  out[6] += 0.35355339059327373*Ghat[5]*dv1; 
  out[7] += 0.35355339059327373*Ghat[6]*dv1; 
  out[8] += 0.6123724356957945*Ghat[1]*dv1; 
  out[9] += 0.6123724356957945*Ghat[2]*dv1; 
  out[10] += 0.6123724356957945*Ghat[3]*dv1; 
  out[11] += 0.35355339059327373*Ghat[7]*dv1; 
  out[12] += 0.35355339059327373*Ghat[8]*dv1; 
  out[13] += 0.35355339059327373*Ghat[9]*dv1; 
  out[14] += 0.7905694150420948*Ghat[0]*dv1; 
  out[15] += 0.35355339059327373*Ghat[10]*dv1; 
  out[16] += 0.6123724356957945*Ghat[4]*dv1; 
  out[17] += 0.6123724356957945*Ghat[5]*dv1; 
  out[18] += 0.6123724356957945*Ghat[6]*dv1; 
  out[19] += 0.35355339059327373*Ghat[11]*dv1; 
  out[20] += 0.35355339059327373*Ghat[12]*dv1; 
  out[21] += 0.35355339059327373*Ghat[13]*dv1; 
  out[22] += 0.35355339059327373*Ghat[14]*dv1; 
  out[23] += 0.35355339059327373*Ghat[15]*dv1; 
  out[24] += 0.35355339059327373*Ghat[16]*dv1; 
  out[25] += 0.6123724356957945*Ghat[7]*dv1; 
  out[26] += 0.6123724356957945*Ghat[8]*dv1; 
  out[27] += 0.6123724356957945*Ghat[9]*dv1; 
  out[28] += 0.7905694150420949*Ghat[1]*dv1; 
  out[29] += 0.7905694150420949*Ghat[2]*dv1; 
  out[30] += 0.7905694150420949*Ghat[3]*dv1; 
  out[31] += 0.6123724356957945*Ghat[10]*dv1; 
  out[32] += 0.35355339059327373*Ghat[17]*dv1; 
  out[33] += 0.35355339059327373*Ghat[18]*dv1; 
  out[34] += 0.35355339059327373*Ghat[19]*dv1; 
  out[35] += 0.6123724356957945*Ghat[11]*dv1; 
  out[36] += 0.6123724356957945*Ghat[12]*dv1; 
  out[37] += 0.6123724356957945*Ghat[13]*dv1; 
  out[38] += 0.6123724356957945*Ghat[14]*dv1; 
  out[39] += 0.6123724356957945*Ghat[15]*dv1; 
  out[40] += 0.6123724356957945*Ghat[16]*dv1; 
  out[41] += 0.7905694150420948*Ghat[4]*dv1; 
  out[42] += 0.7905694150420948*Ghat[5]*dv1; 
  out[43] += 0.7905694150420948*Ghat[6]*dv1; 
  out[44] += 0.6123724356957945*Ghat[17]*dv1; 
  out[45] += 0.6123724356957945*Ghat[18]*dv1; 
  out[46] += 0.6123724356957945*Ghat[19]*dv1; 
  out[47] += 0.7905694150420949*Ghat[10]*dv1; 
} 
