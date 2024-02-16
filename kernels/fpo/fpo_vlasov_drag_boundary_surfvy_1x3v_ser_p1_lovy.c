#include <gkyl_fpo_vlasov_kernels.h> 
#include <gkyl_basis_ser_4x_p1_surfx3_eval_quad.h> 
#include <gkyl_basis_ser_4x_p1_upwind_quad_to_modal.h> 

 
GKYL_CU_DH void fpo_vlasov_drag_boundary_surfvy_1x3v_ser_p1_lovy(const double *dxv, const double* drag_coeff_stencil[9], const double* f_stencil[9], double* out) { 
  // dxv[NDIM]: Cell spacing in each direction. 
  // drag_coeff_stencil[9]: 9-cell stencil of drag coefficient. 
  // f_stencil[9]: 9-cell stencil of distribution function. 
  // out: Incremented output. 


  double dv1 = 2.0/dxv[2]; 
 
  const double* aC = &drag_coeff_stencil[0][40]; 
  const double* aR = &drag_coeff_stencil[1][40]; 
  const double* fC = f_stencil[0]; 
  const double* fR = f_stencil[1]; 

  double fUpwindQuad[8] = {0.0};
  double fUpwind[8] = {0.0};
  double alphaDragSurf[8] = {0.0}; 
  double Ghat[8] = {0.0}; 

  alphaDragSurf[0] = -(0.408248290463863*aR[3])+0.408248290463863*aC[3]+0.3535533905932737*aR[0]+0.3535533905932737*aC[0]; 
  alphaDragSurf[1] = -(0.408248290463863*aR[6])+0.408248290463863*aC[6]+0.3535533905932737*aR[1]+0.3535533905932737*aC[1]; 
  alphaDragSurf[2] = -(0.408248290463863*aR[7])+0.408248290463863*aC[7]+0.3535533905932737*aR[2]+0.3535533905932737*aC[2]; 
  alphaDragSurf[3] = -(0.408248290463863*aR[10])+0.408248290463863*aC[10]+0.3535533905932737*aR[4]+0.3535533905932737*aC[4]; 
  alphaDragSurf[4] = -(0.408248290463863*aR[11])+0.408248290463863*aC[11]+0.3535533905932737*aR[5]+0.3535533905932737*aC[5]; 
  alphaDragSurf[5] = -(0.408248290463863*aR[13])+0.408248290463863*aC[13]+0.3535533905932737*aR[8]+0.3535533905932737*aC[8]; 
  alphaDragSurf[6] = -(0.408248290463863*aR[14])+0.408248290463863*aC[14]+0.3535533905932737*aR[9]+0.3535533905932737*aC[9]; 
  alphaDragSurf[7] = -(0.408248290463863*aR[15])+0.408248290463863*aC[15]+0.3535533905932737*aR[12]+0.3535533905932737*aC[12]; 

  if (-(0.3535533905932737*alphaDragSurf[7])+0.3535533905932737*(alphaDragSurf[6]+alphaDragSurf[5]+alphaDragSurf[4])-0.3535533905932737*(alphaDragSurf[3]+alphaDragSurf[2]+alphaDragSurf[1])+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[0] = ser_4x_p1_surfx3_eval_quad_node_0_r(fC); 
  } else { 
    fUpwindQuad[0] = ser_4x_p1_surfx3_eval_quad_node_0_l(fR); 
  } 
  if (0.3535533905932737*alphaDragSurf[7]-0.3535533905932737*(alphaDragSurf[6]+alphaDragSurf[5])+0.3535533905932737*(alphaDragSurf[4]+alphaDragSurf[3])-0.3535533905932737*(alphaDragSurf[2]+alphaDragSurf[1])+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[1] = ser_4x_p1_surfx3_eval_quad_node_1_r(fC); 
  } else { 
    fUpwindQuad[1] = ser_4x_p1_surfx3_eval_quad_node_1_l(fR); 
  } 
  if (0.3535533905932737*alphaDragSurf[7]-0.3535533905932737*alphaDragSurf[6]+0.3535533905932737*alphaDragSurf[5]-0.3535533905932737*(alphaDragSurf[4]+alphaDragSurf[3])+0.3535533905932737*alphaDragSurf[2]-0.3535533905932737*alphaDragSurf[1]+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[2] = ser_4x_p1_surfx3_eval_quad_node_2_r(fC); 
  } else { 
    fUpwindQuad[2] = ser_4x_p1_surfx3_eval_quad_node_2_l(fR); 
  } 
  if (-(0.3535533905932737*alphaDragSurf[7])+0.3535533905932737*alphaDragSurf[6]-0.3535533905932737*(alphaDragSurf[5]+alphaDragSurf[4])+0.3535533905932737*(alphaDragSurf[3]+alphaDragSurf[2])-0.3535533905932737*alphaDragSurf[1]+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[3] = ser_4x_p1_surfx3_eval_quad_node_3_r(fC); 
  } else { 
    fUpwindQuad[3] = ser_4x_p1_surfx3_eval_quad_node_3_l(fR); 
  } 
  if (0.3535533905932737*(alphaDragSurf[7]+alphaDragSurf[6])-0.3535533905932737*(alphaDragSurf[5]+alphaDragSurf[4]+alphaDragSurf[3]+alphaDragSurf[2])+0.3535533905932737*(alphaDragSurf[1]+alphaDragSurf[0]) > 0) { 
    fUpwindQuad[4] = ser_4x_p1_surfx3_eval_quad_node_4_r(fC); 
  } else { 
    fUpwindQuad[4] = ser_4x_p1_surfx3_eval_quad_node_4_l(fR); 
  } 
  if (-(0.3535533905932737*(alphaDragSurf[7]+alphaDragSurf[6]))+0.3535533905932737*alphaDragSurf[5]-0.3535533905932737*alphaDragSurf[4]+0.3535533905932737*alphaDragSurf[3]-0.3535533905932737*alphaDragSurf[2]+0.3535533905932737*(alphaDragSurf[1]+alphaDragSurf[0]) > 0) { 
    fUpwindQuad[5] = ser_4x_p1_surfx3_eval_quad_node_5_r(fC); 
  } else { 
    fUpwindQuad[5] = ser_4x_p1_surfx3_eval_quad_node_5_l(fR); 
  } 
  if (-(0.3535533905932737*(alphaDragSurf[7]+alphaDragSurf[6]+alphaDragSurf[5]))+0.3535533905932737*alphaDragSurf[4]-0.3535533905932737*alphaDragSurf[3]+0.3535533905932737*(alphaDragSurf[2]+alphaDragSurf[1]+alphaDragSurf[0]) > 0) { 
    fUpwindQuad[6] = ser_4x_p1_surfx3_eval_quad_node_6_r(fC); 
  } else { 
    fUpwindQuad[6] = ser_4x_p1_surfx3_eval_quad_node_6_l(fR); 
  } 
  if (0.3535533905932737*(alphaDragSurf[7]+alphaDragSurf[6]+alphaDragSurf[5]+alphaDragSurf[4]+alphaDragSurf[3]+alphaDragSurf[2]+alphaDragSurf[1]+alphaDragSurf[0]) > 0) { 
    fUpwindQuad[7] = ser_4x_p1_surfx3_eval_quad_node_7_r(fC); 
  } else { 
    fUpwindQuad[7] = ser_4x_p1_surfx3_eval_quad_node_7_l(fR); 
  } 
  ser_4x_p1_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.3535533905932737*alphaDragSurf[7]*fUpwind[7]+0.3535533905932737*alphaDragSurf[6]*fUpwind[6]+0.3535533905932737*alphaDragSurf[5]*fUpwind[5]+0.3535533905932737*alphaDragSurf[4]*fUpwind[4]+0.3535533905932737*alphaDragSurf[3]*fUpwind[3]+0.3535533905932737*alphaDragSurf[2]*fUpwind[2]+0.3535533905932737*alphaDragSurf[1]*fUpwind[1]+0.3535533905932737*alphaDragSurf[0]*fUpwind[0]; 
  Ghat[1] = 0.3535533905932737*alphaDragSurf[6]*fUpwind[7]+0.3535533905932737*fUpwind[6]*alphaDragSurf[7]+0.3535533905932737*alphaDragSurf[3]*fUpwind[5]+0.3535533905932737*fUpwind[3]*alphaDragSurf[5]+0.3535533905932737*alphaDragSurf[2]*fUpwind[4]+0.3535533905932737*fUpwind[2]*alphaDragSurf[4]+0.3535533905932737*alphaDragSurf[0]*fUpwind[1]+0.3535533905932737*fUpwind[0]*alphaDragSurf[1]; 
  Ghat[2] = 0.3535533905932737*alphaDragSurf[5]*fUpwind[7]+0.3535533905932737*fUpwind[5]*alphaDragSurf[7]+0.3535533905932737*alphaDragSurf[3]*fUpwind[6]+0.3535533905932737*fUpwind[3]*alphaDragSurf[6]+0.3535533905932737*alphaDragSurf[1]*fUpwind[4]+0.3535533905932737*fUpwind[1]*alphaDragSurf[4]+0.3535533905932737*alphaDragSurf[0]*fUpwind[2]+0.3535533905932737*fUpwind[0]*alphaDragSurf[2]; 
  Ghat[3] = 0.3535533905932737*alphaDragSurf[4]*fUpwind[7]+0.3535533905932737*fUpwind[4]*alphaDragSurf[7]+0.3535533905932737*alphaDragSurf[2]*fUpwind[6]+0.3535533905932737*fUpwind[2]*alphaDragSurf[6]+0.3535533905932737*alphaDragSurf[1]*fUpwind[5]+0.3535533905932737*fUpwind[1]*alphaDragSurf[5]+0.3535533905932737*alphaDragSurf[0]*fUpwind[3]+0.3535533905932737*fUpwind[0]*alphaDragSurf[3]; 
  Ghat[4] = 0.3535533905932737*alphaDragSurf[3]*fUpwind[7]+0.3535533905932737*fUpwind[3]*alphaDragSurf[7]+0.3535533905932737*alphaDragSurf[5]*fUpwind[6]+0.3535533905932737*fUpwind[5]*alphaDragSurf[6]+0.3535533905932737*alphaDragSurf[0]*fUpwind[4]+0.3535533905932737*fUpwind[0]*alphaDragSurf[4]+0.3535533905932737*alphaDragSurf[1]*fUpwind[2]+0.3535533905932737*fUpwind[1]*alphaDragSurf[2]; 
  Ghat[5] = 0.3535533905932737*alphaDragSurf[2]*fUpwind[7]+0.3535533905932737*fUpwind[2]*alphaDragSurf[7]+0.3535533905932737*alphaDragSurf[4]*fUpwind[6]+0.3535533905932737*fUpwind[4]*alphaDragSurf[6]+0.3535533905932737*alphaDragSurf[0]*fUpwind[5]+0.3535533905932737*fUpwind[0]*alphaDragSurf[5]+0.3535533905932737*alphaDragSurf[1]*fUpwind[3]+0.3535533905932737*fUpwind[1]*alphaDragSurf[3]; 
  Ghat[6] = 0.3535533905932737*alphaDragSurf[1]*fUpwind[7]+0.3535533905932737*fUpwind[1]*alphaDragSurf[7]+0.3535533905932737*alphaDragSurf[0]*fUpwind[6]+0.3535533905932737*fUpwind[0]*alphaDragSurf[6]+0.3535533905932737*alphaDragSurf[4]*fUpwind[5]+0.3535533905932737*fUpwind[4]*alphaDragSurf[5]+0.3535533905932737*alphaDragSurf[2]*fUpwind[3]+0.3535533905932737*fUpwind[2]*alphaDragSurf[3]; 
  Ghat[7] = 0.3535533905932737*alphaDragSurf[0]*fUpwind[7]+0.3535533905932737*fUpwind[0]*alphaDragSurf[7]+0.3535533905932737*alphaDragSurf[1]*fUpwind[6]+0.3535533905932737*fUpwind[1]*alphaDragSurf[6]+0.3535533905932737*alphaDragSurf[2]*fUpwind[5]+0.3535533905932737*fUpwind[2]*alphaDragSurf[5]+0.3535533905932737*alphaDragSurf[3]*fUpwind[4]+0.3535533905932737*fUpwind[3]*alphaDragSurf[4]; 

  out[0] += 0.35355339059327373*Ghat[0]*dv1; 
  out[1] += 0.35355339059327373*Ghat[1]*dv1; 
  out[2] += 0.35355339059327373*Ghat[2]*dv1; 
  out[3] += 0.6123724356957945*Ghat[0]*dv1; 
  out[4] += 0.35355339059327373*Ghat[3]*dv1; 
  out[5] += 0.35355339059327373*Ghat[4]*dv1; 
  out[6] += 0.6123724356957945*Ghat[1]*dv1; 
  out[7] += 0.6123724356957945*Ghat[2]*dv1; 
  out[8] += 0.35355339059327373*Ghat[5]*dv1; 
  out[9] += 0.35355339059327373*Ghat[6]*dv1; 
  out[10] += 0.6123724356957945*Ghat[3]*dv1; 
  out[11] += 0.6123724356957945*Ghat[4]*dv1; 
  out[12] += 0.35355339059327373*Ghat[7]*dv1; 
  out[13] += 0.6123724356957945*Ghat[5]*dv1; 
  out[14] += 0.6123724356957945*Ghat[6]*dv1; 
  out[15] += 0.6123724356957945*Ghat[7]*dv1; 
  out[24] += 0.7905694150420948*Ghat[0]*dv1; 
  out[25] += 0.7905694150420949*Ghat[1]*dv1; 
  out[26] += 0.7905694150420949*Ghat[2]*dv1; 
  out[27] += 0.7905694150420949*Ghat[3]*dv1; 
  out[28] += 0.7905694150420948*Ghat[4]*dv1; 
  out[29] += 0.7905694150420948*Ghat[5]*dv1; 
  out[30] += 0.7905694150420948*Ghat[6]*dv1; 
  out[31] += 0.7905694150420949*Ghat[7]*dv1; 
} 
