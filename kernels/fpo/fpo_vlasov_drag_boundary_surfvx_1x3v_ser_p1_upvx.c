#include <gkyl_fpo_vlasov_kernels.h> 
#include <gkyl_basis_ser_4x_p1_surfx2_eval_quad.h> 
#include <gkyl_basis_ser_4x_p1_upwind_quad_to_modal.h> 

 
GKYL_CU_DH void fpo_vlasov_drag_boundary_surfvx_1x3v_ser_p1_upvx(const double *dxv, const double* drag_coeff_stencil[9], const double* f_stencil[9], double* out) { 
  // dxv[NDIM]: Cell spacing in each direction. 
  // drag_coeff_stencil[9]: 9-cell stencil of drag coefficient. 
  // f_stencil[9]: 9-cell stencil of distribution function. 
  // out: Incremented output. 


  double dv1 = 2.0/dxv[1]; 
 
  const double* aL = &drag_coeff_stencil[0][0]; 
  const double* aC = &drag_coeff_stencil[1][0]; 
  const double* fL = f_stencil[0]; 
  const double* fC = f_stencil[1]; 

  double fUpwindQuad[8] = {0.0};
  double fUpwind[8] = {0.0};
  double alphaDragSurf[8] = {0.0}; 
  double Ghat[8] = {0.0}; 

  alphaDragSurf[0] = 0.408248290463863*aL[2]-0.408248290463863*aC[2]+0.3535533905932737*aL[0]+0.3535533905932737*aC[0]; 
  alphaDragSurf[1] = 0.408248290463863*aL[5]-0.408248290463863*aC[5]+0.3535533905932737*aL[1]+0.3535533905932737*aC[1]; 
  alphaDragSurf[2] = 0.408248290463863*aL[7]-0.408248290463863*aC[7]+0.3535533905932737*aL[3]+0.3535533905932737*aC[3]; 
  alphaDragSurf[3] = 0.408248290463863*aL[9]-0.408248290463863*aC[9]+0.3535533905932737*aL[4]+0.3535533905932737*aC[4]; 
  alphaDragSurf[4] = 0.408248290463863*aL[11]-0.408248290463863*aC[11]+0.3535533905932737*aL[6]+0.3535533905932737*aC[6]; 
  alphaDragSurf[5] = 0.408248290463863*aL[12]-0.408248290463863*aC[12]+0.3535533905932737*aL[8]+0.3535533905932737*aC[8]; 
  alphaDragSurf[6] = 0.408248290463863*aL[14]-0.408248290463863*aC[14]+0.3535533905932737*aL[10]+0.3535533905932737*aC[10]; 
  alphaDragSurf[7] = 0.408248290463863*aL[15]-0.408248290463863*aC[15]+0.3535533905932737*aL[13]+0.3535533905932737*aC[13]; 

  if (-(0.3535533905932737*alphaDragSurf[7])+0.3535533905932737*(alphaDragSurf[6]+alphaDragSurf[5]+alphaDragSurf[4])-0.3535533905932737*(alphaDragSurf[3]+alphaDragSurf[2]+alphaDragSurf[1])+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[0] = ser_4x_p1_surfx2_eval_quad_node_0_r(fL); 
  } else { 
    fUpwindQuad[0] = ser_4x_p1_surfx2_eval_quad_node_0_l(fC); 
  } 
  if (0.3535533905932737*alphaDragSurf[7]-0.3535533905932737*(alphaDragSurf[6]+alphaDragSurf[5])+0.3535533905932737*(alphaDragSurf[4]+alphaDragSurf[3])-0.3535533905932737*(alphaDragSurf[2]+alphaDragSurf[1])+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[1] = ser_4x_p1_surfx2_eval_quad_node_1_r(fL); 
  } else { 
    fUpwindQuad[1] = ser_4x_p1_surfx2_eval_quad_node_1_l(fC); 
  } 
  if (0.3535533905932737*alphaDragSurf[7]-0.3535533905932737*alphaDragSurf[6]+0.3535533905932737*alphaDragSurf[5]-0.3535533905932737*(alphaDragSurf[4]+alphaDragSurf[3])+0.3535533905932737*alphaDragSurf[2]-0.3535533905932737*alphaDragSurf[1]+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[2] = ser_4x_p1_surfx2_eval_quad_node_2_r(fL); 
  } else { 
    fUpwindQuad[2] = ser_4x_p1_surfx2_eval_quad_node_2_l(fC); 
  } 
  if (-(0.3535533905932737*alphaDragSurf[7])+0.3535533905932737*alphaDragSurf[6]-0.3535533905932737*(alphaDragSurf[5]+alphaDragSurf[4])+0.3535533905932737*(alphaDragSurf[3]+alphaDragSurf[2])-0.3535533905932737*alphaDragSurf[1]+0.3535533905932737*alphaDragSurf[0] > 0) { 
    fUpwindQuad[3] = ser_4x_p1_surfx2_eval_quad_node_3_r(fL); 
  } else { 
    fUpwindQuad[3] = ser_4x_p1_surfx2_eval_quad_node_3_l(fC); 
  } 
  if (0.3535533905932737*(alphaDragSurf[7]+alphaDragSurf[6])-0.3535533905932737*(alphaDragSurf[5]+alphaDragSurf[4]+alphaDragSurf[3]+alphaDragSurf[2])+0.3535533905932737*(alphaDragSurf[1]+alphaDragSurf[0]) > 0) { 
    fUpwindQuad[4] = ser_4x_p1_surfx2_eval_quad_node_4_r(fL); 
  } else { 
    fUpwindQuad[4] = ser_4x_p1_surfx2_eval_quad_node_4_l(fC); 
  } 
  if (-(0.3535533905932737*(alphaDragSurf[7]+alphaDragSurf[6]))+0.3535533905932737*alphaDragSurf[5]-0.3535533905932737*alphaDragSurf[4]+0.3535533905932737*alphaDragSurf[3]-0.3535533905932737*alphaDragSurf[2]+0.3535533905932737*(alphaDragSurf[1]+alphaDragSurf[0]) > 0) { 
    fUpwindQuad[5] = ser_4x_p1_surfx2_eval_quad_node_5_r(fL); 
  } else { 
    fUpwindQuad[5] = ser_4x_p1_surfx2_eval_quad_node_5_l(fC); 
  } 
  if (-(0.3535533905932737*(alphaDragSurf[7]+alphaDragSurf[6]+alphaDragSurf[5]))+0.3535533905932737*alphaDragSurf[4]-0.3535533905932737*alphaDragSurf[3]+0.3535533905932737*(alphaDragSurf[2]+alphaDragSurf[1]+alphaDragSurf[0]) > 0) { 
    fUpwindQuad[6] = ser_4x_p1_surfx2_eval_quad_node_6_r(fL); 
  } else { 
    fUpwindQuad[6] = ser_4x_p1_surfx2_eval_quad_node_6_l(fC); 
  } 
  if (0.3535533905932737*(alphaDragSurf[7]+alphaDragSurf[6]+alphaDragSurf[5]+alphaDragSurf[4]+alphaDragSurf[3]+alphaDragSurf[2]+alphaDragSurf[1]+alphaDragSurf[0]) > 0) { 
    fUpwindQuad[7] = ser_4x_p1_surfx2_eval_quad_node_7_r(fL); 
  } else { 
    fUpwindQuad[7] = ser_4x_p1_surfx2_eval_quad_node_7_l(fC); 
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
  out[11] += 0.6123724356957945*Ghat[4]*dv1; 
  out[12] += 0.6123724356957945*Ghat[5]*dv1; 
  out[13] += -(0.35355339059327373*Ghat[7]*dv1); 
  out[14] += 0.6123724356957945*Ghat[6]*dv1; 
  out[15] += 0.6123724356957945*Ghat[7]*dv1; 
  out[16] += -(0.7905694150420948*Ghat[0]*dv1); 
  out[17] += -(0.7905694150420949*Ghat[1]*dv1); 
  out[18] += -(0.7905694150420949*Ghat[2]*dv1); 
  out[19] += -(0.7905694150420949*Ghat[3]*dv1); 
  out[20] += -(0.7905694150420948*Ghat[4]*dv1); 
  out[21] += -(0.7905694150420948*Ghat[5]*dv1); 
  out[22] += -(0.7905694150420948*Ghat[6]*dv1); 
  out[23] += -(0.7905694150420949*Ghat[7]*dv1); 
} 