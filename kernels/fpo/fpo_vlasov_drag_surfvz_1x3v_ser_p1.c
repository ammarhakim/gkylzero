#include <gkyl_fpo_vlasov_kernels.h> 
#include <gkyl_basis_ser_4x_p1_surfx4_eval_quad.h> 
#include <gkyl_basis_ser_4x_p1_upwind_quad_to_modal.h> 

GKYL_CU_DH void fpo_vlasov_drag_surfvz_1x3v_ser_p1(const double* dxv, const double* drag_coeff_stencil[9], const double* f_stencil[9], double* GKYL_RESTRICT out) { 
  // dxv[NDIM]: Cell spacing in each direction. 
  // drag_coeff: Drag coefficient. 
  // f: Distribution function. 
  // out: Incremented output. 


  const double* fL = f_stencil[0]; 
  const double* fC = f_stencil[1]; 
  const double* fR = f_stencil[2]; 

  const double* aL = &drag_coeff_stencil[0][80]; 
  const double* aC = &drag_coeff_stencil[1][80]; 
  const double* aR = &drag_coeff_stencil[2][80]; 
  double dv1 = 2.0/dxv[3]; 

  double alphaDragSurf_lo[8] = {0.0}; 
  alphaDragSurf_lo[0] = 0.408248290463863*aL[4]-0.408248290463863*aC[4]+0.3535533905932737*aL[0]+0.3535533905932737*aC[0]; 
  alphaDragSurf_lo[1] = 0.408248290463863*aL[8]-0.408248290463863*aC[8]+0.3535533905932737*aL[1]+0.3535533905932737*aC[1]; 
  alphaDragSurf_lo[2] = 0.408248290463863*aL[9]-0.408248290463863*aC[9]+0.3535533905932737*aL[2]+0.3535533905932737*aC[2]; 
  alphaDragSurf_lo[3] = 0.408248290463863*aL[10]-0.408248290463863*aC[10]+0.3535533905932737*aL[3]+0.3535533905932737*aC[3]; 
  alphaDragSurf_lo[4] = 0.408248290463863*aL[12]-0.408248290463863*aC[12]+0.3535533905932737*aL[5]+0.3535533905932737*aC[5]; 
  alphaDragSurf_lo[5] = 0.408248290463863*aL[13]-0.408248290463863*aC[13]+0.3535533905932737*aL[6]+0.3535533905932737*aC[6]; 
  alphaDragSurf_lo[6] = 0.408248290463863*aL[14]-0.408248290463863*aC[14]+0.3535533905932737*aL[7]+0.3535533905932737*aC[7]; 
  alphaDragSurf_lo[7] = 0.408248290463863*aL[15]-0.408248290463863*aC[15]+0.3535533905932737*aL[11]+0.3535533905932737*aC[11]; 

  double alphaDragSurf_up[8] = {0.0}; 
  alphaDragSurf_up[0] = -(0.408248290463863*aR[4])+0.408248290463863*aC[4]+0.3535533905932737*aR[0]+0.3535533905932737*aC[0]; 
  alphaDragSurf_up[1] = -(0.408248290463863*aR[8])+0.408248290463863*aC[8]+0.3535533905932737*aR[1]+0.3535533905932737*aC[1]; 
  alphaDragSurf_up[2] = -(0.408248290463863*aR[9])+0.408248290463863*aC[9]+0.3535533905932737*aR[2]+0.3535533905932737*aC[2]; 
  alphaDragSurf_up[3] = -(0.408248290463863*aR[10])+0.408248290463863*aC[10]+0.3535533905932737*aR[3]+0.3535533905932737*aC[3]; 
  alphaDragSurf_up[4] = -(0.408248290463863*aR[12])+0.408248290463863*aC[12]+0.3535533905932737*aR[5]+0.3535533905932737*aC[5]; 
  alphaDragSurf_up[5] = -(0.408248290463863*aR[13])+0.408248290463863*aC[13]+0.3535533905932737*aR[6]+0.3535533905932737*aC[6]; 
  alphaDragSurf_up[6] = -(0.408248290463863*aR[14])+0.408248290463863*aC[14]+0.3535533905932737*aR[7]+0.3535533905932737*aC[7]; 
  alphaDragSurf_up[7] = -(0.408248290463863*aR[15])+0.408248290463863*aC[15]+0.3535533905932737*aR[11]+0.3535533905932737*aC[11]; 

  double fUpwindQuad_lo[8] = {0.0};
  double fUpwindQuad_up[8] = {0.0};
  double fUpwind_lo[8] = {0.0};
  double fUpwind_up[8] = {0.0};
  double Ghat_lo[8] = {0.0}; 
  double Ghat_up[8] = {0.0}; 

  if (-(0.3535533905932737*alphaDragSurf_lo[7])+0.3535533905932737*(alphaDragSurf_lo[6]+alphaDragSurf_lo[5]+alphaDragSurf_lo[4])-0.3535533905932737*(alphaDragSurf_lo[3]+alphaDragSurf_lo[2]+alphaDragSurf_lo[1])+0.3535533905932737*alphaDragSurf_lo[0] > 0) { 
    fUpwindQuad_lo[0] = ser_4x_p1_surfx4_eval_quad_node_0_r(fL); 
  } else { 
    fUpwindQuad_lo[0] = ser_4x_p1_surfx4_eval_quad_node_0_l(fC); 
  } 
  if (-(0.3535533905932737*alphaDragSurf_up[7])+0.3535533905932737*(alphaDragSurf_up[6]+alphaDragSurf_up[5]+alphaDragSurf_up[4])-0.3535533905932737*(alphaDragSurf_up[3]+alphaDragSurf_up[2]+alphaDragSurf_up[1])+0.3535533905932737*alphaDragSurf_up[0] > 0) { 
    fUpwindQuad_up[0] = ser_4x_p1_surfx4_eval_quad_node_0_r(fC); 
  } else { 
    fUpwindQuad_up[0] = ser_4x_p1_surfx4_eval_quad_node_0_l(fR); 
  } 
  if (0.3535533905932737*alphaDragSurf_lo[7]-0.3535533905932737*(alphaDragSurf_lo[6]+alphaDragSurf_lo[5])+0.3535533905932737*(alphaDragSurf_lo[4]+alphaDragSurf_lo[3])-0.3535533905932737*(alphaDragSurf_lo[2]+alphaDragSurf_lo[1])+0.3535533905932737*alphaDragSurf_lo[0] > 0) { 
    fUpwindQuad_lo[1] = ser_4x_p1_surfx4_eval_quad_node_1_r(fL); 
  } else { 
    fUpwindQuad_lo[1] = ser_4x_p1_surfx4_eval_quad_node_1_l(fC); 
  } 
  if (0.3535533905932737*alphaDragSurf_up[7]-0.3535533905932737*(alphaDragSurf_up[6]+alphaDragSurf_up[5])+0.3535533905932737*(alphaDragSurf_up[4]+alphaDragSurf_up[3])-0.3535533905932737*(alphaDragSurf_up[2]+alphaDragSurf_up[1])+0.3535533905932737*alphaDragSurf_up[0] > 0) { 
    fUpwindQuad_up[1] = ser_4x_p1_surfx4_eval_quad_node_1_r(fC); 
  } else { 
    fUpwindQuad_up[1] = ser_4x_p1_surfx4_eval_quad_node_1_l(fR); 
  } 
  if (0.3535533905932737*alphaDragSurf_lo[7]-0.3535533905932737*alphaDragSurf_lo[6]+0.3535533905932737*alphaDragSurf_lo[5]-0.3535533905932737*(alphaDragSurf_lo[4]+alphaDragSurf_lo[3])+0.3535533905932737*alphaDragSurf_lo[2]-0.3535533905932737*alphaDragSurf_lo[1]+0.3535533905932737*alphaDragSurf_lo[0] > 0) { 
    fUpwindQuad_lo[2] = ser_4x_p1_surfx4_eval_quad_node_2_r(fL); 
  } else { 
    fUpwindQuad_lo[2] = ser_4x_p1_surfx4_eval_quad_node_2_l(fC); 
  } 
  if (0.3535533905932737*alphaDragSurf_up[7]-0.3535533905932737*alphaDragSurf_up[6]+0.3535533905932737*alphaDragSurf_up[5]-0.3535533905932737*(alphaDragSurf_up[4]+alphaDragSurf_up[3])+0.3535533905932737*alphaDragSurf_up[2]-0.3535533905932737*alphaDragSurf_up[1]+0.3535533905932737*alphaDragSurf_up[0] > 0) { 
    fUpwindQuad_up[2] = ser_4x_p1_surfx4_eval_quad_node_2_r(fC); 
  } else { 
    fUpwindQuad_up[2] = ser_4x_p1_surfx4_eval_quad_node_2_l(fR); 
  } 
  if (-(0.3535533905932737*alphaDragSurf_lo[7])+0.3535533905932737*alphaDragSurf_lo[6]-0.3535533905932737*(alphaDragSurf_lo[5]+alphaDragSurf_lo[4])+0.3535533905932737*(alphaDragSurf_lo[3]+alphaDragSurf_lo[2])-0.3535533905932737*alphaDragSurf_lo[1]+0.3535533905932737*alphaDragSurf_lo[0] > 0) { 
    fUpwindQuad_lo[3] = ser_4x_p1_surfx4_eval_quad_node_3_r(fL); 
  } else { 
    fUpwindQuad_lo[3] = ser_4x_p1_surfx4_eval_quad_node_3_l(fC); 
  } 
  if (-(0.3535533905932737*alphaDragSurf_up[7])+0.3535533905932737*alphaDragSurf_up[6]-0.3535533905932737*(alphaDragSurf_up[5]+alphaDragSurf_up[4])+0.3535533905932737*(alphaDragSurf_up[3]+alphaDragSurf_up[2])-0.3535533905932737*alphaDragSurf_up[1]+0.3535533905932737*alphaDragSurf_up[0] > 0) { 
    fUpwindQuad_up[3] = ser_4x_p1_surfx4_eval_quad_node_3_r(fC); 
  } else { 
    fUpwindQuad_up[3] = ser_4x_p1_surfx4_eval_quad_node_3_l(fR); 
  } 
  if (0.3535533905932737*(alphaDragSurf_lo[7]+alphaDragSurf_lo[6])-0.3535533905932737*(alphaDragSurf_lo[5]+alphaDragSurf_lo[4]+alphaDragSurf_lo[3]+alphaDragSurf_lo[2])+0.3535533905932737*(alphaDragSurf_lo[1]+alphaDragSurf_lo[0]) > 0) { 
    fUpwindQuad_lo[4] = ser_4x_p1_surfx4_eval_quad_node_4_r(fL); 
  } else { 
    fUpwindQuad_lo[4] = ser_4x_p1_surfx4_eval_quad_node_4_l(fC); 
  } 
  if (0.3535533905932737*(alphaDragSurf_up[7]+alphaDragSurf_up[6])-0.3535533905932737*(alphaDragSurf_up[5]+alphaDragSurf_up[4]+alphaDragSurf_up[3]+alphaDragSurf_up[2])+0.3535533905932737*(alphaDragSurf_up[1]+alphaDragSurf_up[0]) > 0) { 
    fUpwindQuad_up[4] = ser_4x_p1_surfx4_eval_quad_node_4_r(fC); 
  } else { 
    fUpwindQuad_up[4] = ser_4x_p1_surfx4_eval_quad_node_4_l(fR); 
  } 
  if (-(0.3535533905932737*(alphaDragSurf_lo[7]+alphaDragSurf_lo[6]))+0.3535533905932737*alphaDragSurf_lo[5]-0.3535533905932737*alphaDragSurf_lo[4]+0.3535533905932737*alphaDragSurf_lo[3]-0.3535533905932737*alphaDragSurf_lo[2]+0.3535533905932737*(alphaDragSurf_lo[1]+alphaDragSurf_lo[0]) > 0) { 
    fUpwindQuad_lo[5] = ser_4x_p1_surfx4_eval_quad_node_5_r(fL); 
  } else { 
    fUpwindQuad_lo[5] = ser_4x_p1_surfx4_eval_quad_node_5_l(fC); 
  } 
  if (-(0.3535533905932737*(alphaDragSurf_up[7]+alphaDragSurf_up[6]))+0.3535533905932737*alphaDragSurf_up[5]-0.3535533905932737*alphaDragSurf_up[4]+0.3535533905932737*alphaDragSurf_up[3]-0.3535533905932737*alphaDragSurf_up[2]+0.3535533905932737*(alphaDragSurf_up[1]+alphaDragSurf_up[0]) > 0) { 
    fUpwindQuad_up[5] = ser_4x_p1_surfx4_eval_quad_node_5_r(fC); 
  } else { 
    fUpwindQuad_up[5] = ser_4x_p1_surfx4_eval_quad_node_5_l(fR); 
  } 
  if (-(0.3535533905932737*(alphaDragSurf_lo[7]+alphaDragSurf_lo[6]+alphaDragSurf_lo[5]))+0.3535533905932737*alphaDragSurf_lo[4]-0.3535533905932737*alphaDragSurf_lo[3]+0.3535533905932737*(alphaDragSurf_lo[2]+alphaDragSurf_lo[1]+alphaDragSurf_lo[0]) > 0) { 
    fUpwindQuad_lo[6] = ser_4x_p1_surfx4_eval_quad_node_6_r(fL); 
  } else { 
    fUpwindQuad_lo[6] = ser_4x_p1_surfx4_eval_quad_node_6_l(fC); 
  } 
  if (-(0.3535533905932737*(alphaDragSurf_up[7]+alphaDragSurf_up[6]+alphaDragSurf_up[5]))+0.3535533905932737*alphaDragSurf_up[4]-0.3535533905932737*alphaDragSurf_up[3]+0.3535533905932737*(alphaDragSurf_up[2]+alphaDragSurf_up[1]+alphaDragSurf_up[0]) > 0) { 
    fUpwindQuad_up[6] = ser_4x_p1_surfx4_eval_quad_node_6_r(fC); 
  } else { 
    fUpwindQuad_up[6] = ser_4x_p1_surfx4_eval_quad_node_6_l(fR); 
  } 
  if (0.3535533905932737*(alphaDragSurf_lo[7]+alphaDragSurf_lo[6]+alphaDragSurf_lo[5]+alphaDragSurf_lo[4]+alphaDragSurf_lo[3]+alphaDragSurf_lo[2]+alphaDragSurf_lo[1]+alphaDragSurf_lo[0]) > 0) { 
    fUpwindQuad_lo[7] = ser_4x_p1_surfx4_eval_quad_node_7_r(fL); 
  } else { 
    fUpwindQuad_lo[7] = ser_4x_p1_surfx4_eval_quad_node_7_l(fC); 
  } 
  if (0.3535533905932737*(alphaDragSurf_up[7]+alphaDragSurf_up[6]+alphaDragSurf_up[5]+alphaDragSurf_up[4]+alphaDragSurf_up[3]+alphaDragSurf_up[2]+alphaDragSurf_up[1]+alphaDragSurf_up[0]) > 0) { 
    fUpwindQuad_up[7] = ser_4x_p1_surfx4_eval_quad_node_7_r(fC); 
  } else { 
    fUpwindQuad_up[7] = ser_4x_p1_surfx4_eval_quad_node_7_l(fR); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_4x_p1_upwind_quad_to_modal(fUpwindQuad_lo, fUpwind_lo); 
  ser_4x_p1_upwind_quad_to_modal(fUpwindQuad_up, fUpwind_up); 

  Ghat_lo[0] = 0.3535533905932737*alphaDragSurf_lo[7]*fUpwind_lo[7]+0.3535533905932737*alphaDragSurf_lo[6]*fUpwind_lo[6]+0.3535533905932737*alphaDragSurf_lo[5]*fUpwind_lo[5]+0.3535533905932737*alphaDragSurf_lo[4]*fUpwind_lo[4]+0.3535533905932737*alphaDragSurf_lo[3]*fUpwind_lo[3]+0.3535533905932737*alphaDragSurf_lo[2]*fUpwind_lo[2]+0.3535533905932737*alphaDragSurf_lo[1]*fUpwind_lo[1]+0.3535533905932737*alphaDragSurf_lo[0]*fUpwind_lo[0]; 
  Ghat_lo[1] = 0.3535533905932737*alphaDragSurf_lo[6]*fUpwind_lo[7]+0.3535533905932737*fUpwind_lo[6]*alphaDragSurf_lo[7]+0.3535533905932737*alphaDragSurf_lo[3]*fUpwind_lo[5]+0.3535533905932737*fUpwind_lo[3]*alphaDragSurf_lo[5]+0.3535533905932737*alphaDragSurf_lo[2]*fUpwind_lo[4]+0.3535533905932737*fUpwind_lo[2]*alphaDragSurf_lo[4]+0.3535533905932737*alphaDragSurf_lo[0]*fUpwind_lo[1]+0.3535533905932737*fUpwind_lo[0]*alphaDragSurf_lo[1]; 
  Ghat_lo[2] = 0.3535533905932737*alphaDragSurf_lo[5]*fUpwind_lo[7]+0.3535533905932737*fUpwind_lo[5]*alphaDragSurf_lo[7]+0.3535533905932737*alphaDragSurf_lo[3]*fUpwind_lo[6]+0.3535533905932737*fUpwind_lo[3]*alphaDragSurf_lo[6]+0.3535533905932737*alphaDragSurf_lo[1]*fUpwind_lo[4]+0.3535533905932737*fUpwind_lo[1]*alphaDragSurf_lo[4]+0.3535533905932737*alphaDragSurf_lo[0]*fUpwind_lo[2]+0.3535533905932737*fUpwind_lo[0]*alphaDragSurf_lo[2]; 
  Ghat_lo[3] = 0.3535533905932737*alphaDragSurf_lo[4]*fUpwind_lo[7]+0.3535533905932737*fUpwind_lo[4]*alphaDragSurf_lo[7]+0.3535533905932737*alphaDragSurf_lo[2]*fUpwind_lo[6]+0.3535533905932737*fUpwind_lo[2]*alphaDragSurf_lo[6]+0.3535533905932737*alphaDragSurf_lo[1]*fUpwind_lo[5]+0.3535533905932737*fUpwind_lo[1]*alphaDragSurf_lo[5]+0.3535533905932737*alphaDragSurf_lo[0]*fUpwind_lo[3]+0.3535533905932737*fUpwind_lo[0]*alphaDragSurf_lo[3]; 
  Ghat_lo[4] = 0.3535533905932737*alphaDragSurf_lo[3]*fUpwind_lo[7]+0.3535533905932737*fUpwind_lo[3]*alphaDragSurf_lo[7]+0.3535533905932737*alphaDragSurf_lo[5]*fUpwind_lo[6]+0.3535533905932737*fUpwind_lo[5]*alphaDragSurf_lo[6]+0.3535533905932737*alphaDragSurf_lo[0]*fUpwind_lo[4]+0.3535533905932737*fUpwind_lo[0]*alphaDragSurf_lo[4]+0.3535533905932737*alphaDragSurf_lo[1]*fUpwind_lo[2]+0.3535533905932737*fUpwind_lo[1]*alphaDragSurf_lo[2]; 
  Ghat_lo[5] = 0.3535533905932737*alphaDragSurf_lo[2]*fUpwind_lo[7]+0.3535533905932737*fUpwind_lo[2]*alphaDragSurf_lo[7]+0.3535533905932737*alphaDragSurf_lo[4]*fUpwind_lo[6]+0.3535533905932737*fUpwind_lo[4]*alphaDragSurf_lo[6]+0.3535533905932737*alphaDragSurf_lo[0]*fUpwind_lo[5]+0.3535533905932737*fUpwind_lo[0]*alphaDragSurf_lo[5]+0.3535533905932737*alphaDragSurf_lo[1]*fUpwind_lo[3]+0.3535533905932737*fUpwind_lo[1]*alphaDragSurf_lo[3]; 
  Ghat_lo[6] = 0.3535533905932737*alphaDragSurf_lo[1]*fUpwind_lo[7]+0.3535533905932737*fUpwind_lo[1]*alphaDragSurf_lo[7]+0.3535533905932737*alphaDragSurf_lo[0]*fUpwind_lo[6]+0.3535533905932737*fUpwind_lo[0]*alphaDragSurf_lo[6]+0.3535533905932737*alphaDragSurf_lo[4]*fUpwind_lo[5]+0.3535533905932737*fUpwind_lo[4]*alphaDragSurf_lo[5]+0.3535533905932737*alphaDragSurf_lo[2]*fUpwind_lo[3]+0.3535533905932737*fUpwind_lo[2]*alphaDragSurf_lo[3]; 
  Ghat_lo[7] = 0.3535533905932737*alphaDragSurf_lo[0]*fUpwind_lo[7]+0.3535533905932737*fUpwind_lo[0]*alphaDragSurf_lo[7]+0.3535533905932737*alphaDragSurf_lo[1]*fUpwind_lo[6]+0.3535533905932737*fUpwind_lo[1]*alphaDragSurf_lo[6]+0.3535533905932737*alphaDragSurf_lo[2]*fUpwind_lo[5]+0.3535533905932737*fUpwind_lo[2]*alphaDragSurf_lo[5]+0.3535533905932737*alphaDragSurf_lo[3]*fUpwind_lo[4]+0.3535533905932737*fUpwind_lo[3]*alphaDragSurf_lo[4]; 

  Ghat_up[0] = 0.3535533905932737*alphaDragSurf_up[7]*fUpwind_up[7]+0.3535533905932737*alphaDragSurf_up[6]*fUpwind_up[6]+0.3535533905932737*alphaDragSurf_up[5]*fUpwind_up[5]+0.3535533905932737*alphaDragSurf_up[4]*fUpwind_up[4]+0.3535533905932737*alphaDragSurf_up[3]*fUpwind_up[3]+0.3535533905932737*alphaDragSurf_up[2]*fUpwind_up[2]+0.3535533905932737*alphaDragSurf_up[1]*fUpwind_up[1]+0.3535533905932737*alphaDragSurf_up[0]*fUpwind_up[0]; 
  Ghat_up[1] = 0.3535533905932737*alphaDragSurf_up[6]*fUpwind_up[7]+0.3535533905932737*fUpwind_up[6]*alphaDragSurf_up[7]+0.3535533905932737*alphaDragSurf_up[3]*fUpwind_up[5]+0.3535533905932737*fUpwind_up[3]*alphaDragSurf_up[5]+0.3535533905932737*alphaDragSurf_up[2]*fUpwind_up[4]+0.3535533905932737*fUpwind_up[2]*alphaDragSurf_up[4]+0.3535533905932737*alphaDragSurf_up[0]*fUpwind_up[1]+0.3535533905932737*fUpwind_up[0]*alphaDragSurf_up[1]; 
  Ghat_up[2] = 0.3535533905932737*alphaDragSurf_up[5]*fUpwind_up[7]+0.3535533905932737*fUpwind_up[5]*alphaDragSurf_up[7]+0.3535533905932737*alphaDragSurf_up[3]*fUpwind_up[6]+0.3535533905932737*fUpwind_up[3]*alphaDragSurf_up[6]+0.3535533905932737*alphaDragSurf_up[1]*fUpwind_up[4]+0.3535533905932737*fUpwind_up[1]*alphaDragSurf_up[4]+0.3535533905932737*alphaDragSurf_up[0]*fUpwind_up[2]+0.3535533905932737*fUpwind_up[0]*alphaDragSurf_up[2]; 
  Ghat_up[3] = 0.3535533905932737*alphaDragSurf_up[4]*fUpwind_up[7]+0.3535533905932737*fUpwind_up[4]*alphaDragSurf_up[7]+0.3535533905932737*alphaDragSurf_up[2]*fUpwind_up[6]+0.3535533905932737*fUpwind_up[2]*alphaDragSurf_up[6]+0.3535533905932737*alphaDragSurf_up[1]*fUpwind_up[5]+0.3535533905932737*fUpwind_up[1]*alphaDragSurf_up[5]+0.3535533905932737*alphaDragSurf_up[0]*fUpwind_up[3]+0.3535533905932737*fUpwind_up[0]*alphaDragSurf_up[3]; 
  Ghat_up[4] = 0.3535533905932737*alphaDragSurf_up[3]*fUpwind_up[7]+0.3535533905932737*fUpwind_up[3]*alphaDragSurf_up[7]+0.3535533905932737*alphaDragSurf_up[5]*fUpwind_up[6]+0.3535533905932737*fUpwind_up[5]*alphaDragSurf_up[6]+0.3535533905932737*alphaDragSurf_up[0]*fUpwind_up[4]+0.3535533905932737*fUpwind_up[0]*alphaDragSurf_up[4]+0.3535533905932737*alphaDragSurf_up[1]*fUpwind_up[2]+0.3535533905932737*fUpwind_up[1]*alphaDragSurf_up[2]; 
  Ghat_up[5] = 0.3535533905932737*alphaDragSurf_up[2]*fUpwind_up[7]+0.3535533905932737*fUpwind_up[2]*alphaDragSurf_up[7]+0.3535533905932737*alphaDragSurf_up[4]*fUpwind_up[6]+0.3535533905932737*fUpwind_up[4]*alphaDragSurf_up[6]+0.3535533905932737*alphaDragSurf_up[0]*fUpwind_up[5]+0.3535533905932737*fUpwind_up[0]*alphaDragSurf_up[5]+0.3535533905932737*alphaDragSurf_up[1]*fUpwind_up[3]+0.3535533905932737*fUpwind_up[1]*alphaDragSurf_up[3]; 
  Ghat_up[6] = 0.3535533905932737*alphaDragSurf_up[1]*fUpwind_up[7]+0.3535533905932737*fUpwind_up[1]*alphaDragSurf_up[7]+0.3535533905932737*alphaDragSurf_up[0]*fUpwind_up[6]+0.3535533905932737*fUpwind_up[0]*alphaDragSurf_up[6]+0.3535533905932737*alphaDragSurf_up[4]*fUpwind_up[5]+0.3535533905932737*fUpwind_up[4]*alphaDragSurf_up[5]+0.3535533905932737*alphaDragSurf_up[2]*fUpwind_up[3]+0.3535533905932737*fUpwind_up[2]*alphaDragSurf_up[3]; 
  Ghat_up[7] = 0.3535533905932737*alphaDragSurf_up[0]*fUpwind_up[7]+0.3535533905932737*fUpwind_up[0]*alphaDragSurf_up[7]+0.3535533905932737*alphaDragSurf_up[1]*fUpwind_up[6]+0.3535533905932737*fUpwind_up[1]*alphaDragSurf_up[6]+0.3535533905932737*alphaDragSurf_up[2]*fUpwind_up[5]+0.3535533905932737*fUpwind_up[2]*alphaDragSurf_up[5]+0.3535533905932737*alphaDragSurf_up[3]*fUpwind_up[4]+0.3535533905932737*fUpwind_up[3]*alphaDragSurf_up[4]; 

  out[0] += 0.35355339059327373*Ghat_up[0]*dv1-0.35355339059327373*Ghat_lo[0]*dv1; 
  out[1] += 0.35355339059327373*Ghat_up[1]*dv1-0.35355339059327373*Ghat_lo[1]*dv1; 
  out[2] += 0.35355339059327373*Ghat_up[2]*dv1-0.35355339059327373*Ghat_lo[2]*dv1; 
  out[3] += 0.35355339059327373*Ghat_up[3]*dv1-0.35355339059327373*Ghat_lo[3]*dv1; 
  out[4] += 0.6123724356957945*Ghat_up[0]*dv1+0.6123724356957945*Ghat_lo[0]*dv1; 
  out[5] += 0.35355339059327373*Ghat_up[4]*dv1-0.35355339059327373*Ghat_lo[4]*dv1; 
  out[6] += 0.35355339059327373*Ghat_up[5]*dv1-0.35355339059327373*Ghat_lo[5]*dv1; 
  out[7] += 0.35355339059327373*Ghat_up[6]*dv1-0.35355339059327373*Ghat_lo[6]*dv1; 
  out[8] += 0.6123724356957945*Ghat_up[1]*dv1+0.6123724356957945*Ghat_lo[1]*dv1; 
  out[9] += 0.6123724356957945*Ghat_up[2]*dv1+0.6123724356957945*Ghat_lo[2]*dv1; 
  out[10] += 0.6123724356957945*Ghat_up[3]*dv1+0.6123724356957945*Ghat_lo[3]*dv1; 
  out[11] += 0.35355339059327373*Ghat_up[7]*dv1-0.35355339059327373*Ghat_lo[7]*dv1; 
  out[12] += 0.6123724356957945*Ghat_up[4]*dv1+0.6123724356957945*Ghat_lo[4]*dv1; 
  out[13] += 0.6123724356957945*Ghat_up[5]*dv1+0.6123724356957945*Ghat_lo[5]*dv1; 
  out[14] += 0.6123724356957945*Ghat_up[6]*dv1+0.6123724356957945*Ghat_lo[6]*dv1; 
  out[15] += 0.6123724356957945*Ghat_up[7]*dv1+0.6123724356957945*Ghat_lo[7]*dv1; 
  out[32] += 0.7905694150420948*Ghat_up[0]*dv1-0.7905694150420948*Ghat_lo[0]*dv1; 
  out[33] += 0.7905694150420949*Ghat_up[1]*dv1-0.7905694150420949*Ghat_lo[1]*dv1; 
  out[34] += 0.7905694150420949*Ghat_up[2]*dv1-0.7905694150420949*Ghat_lo[2]*dv1; 
  out[35] += 0.7905694150420949*Ghat_up[3]*dv1-0.7905694150420949*Ghat_lo[3]*dv1; 
  out[36] += 0.7905694150420948*Ghat_up[4]*dv1-0.7905694150420948*Ghat_lo[4]*dv1; 
  out[37] += 0.7905694150420948*Ghat_up[5]*dv1-0.7905694150420948*Ghat_lo[5]*dv1; 
  out[38] += 0.7905694150420948*Ghat_up[6]*dv1-0.7905694150420948*Ghat_lo[6]*dv1; 
  out[39] += 0.7905694150420949*Ghat_up[7]*dv1-0.7905694150420949*Ghat_lo[7]*dv1; 

}
