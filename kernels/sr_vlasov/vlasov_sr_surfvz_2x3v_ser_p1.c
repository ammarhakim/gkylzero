#include <gkyl_vlasov_sr_kernels.h> 
#include <gkyl_basis_ser_5x_p1_surfx5_eval_quad.h> 
#include <gkyl_basis_ser_5x_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double vlasov_sr_surfvz_2x3v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // p_over_gamma:      p/gamma (velocity).
  // qmem:      q/m*EM fields.
  // fl/fc/fr:  Input Distribution function in left/center/right cells 
  // out:       Output distribution function in center cell 
  const double dv12 = 2/dxv[4]; 
  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv2 = dxv[3], wv2 = w[3]; 
  const double dv3 = dxv[4], wv3 = w[4]; 
  const double *E2 = &qmem[8]; 
  const double *p0_over_gamma = &p_over_gamma[0]; 
  const double *B0 = &qmem[12]; 
  const double *p1_over_gamma = &p_over_gamma[8]; 
  const double *B1 = &qmem[16]; 
  const double *p2_over_gamma = &p_over_gamma[16]; 
  const double *B2 = &qmem[20]; 

  double alpha_l[16] = {0.0}; 
  double alpha_r[16] = {0.0}; 

  alpha_l[0] = 1.224744871391589*B0[0]*p1_over_gamma[3]-1.224744871391589*B1[0]*p0_over_gamma[3]-0.7071067811865475*B0[0]*p1_over_gamma[0]+0.7071067811865475*B1[0]*p0_over_gamma[0]+2.0*E2[0]; 
  alpha_l[1] = 1.224744871391589*B0[1]*p1_over_gamma[3]-1.224744871391589*B1[1]*p0_over_gamma[3]+2.0*E2[1]+0.7071067811865475*p0_over_gamma[0]*B1[1]-0.7071067811865475*p1_over_gamma[0]*B0[1]; 
  alpha_l[2] = 1.224744871391589*B0[2]*p1_over_gamma[3]-1.224744871391589*B1[2]*p0_over_gamma[3]+2.0*E2[2]+0.7071067811865475*p0_over_gamma[0]*B1[2]-0.7071067811865475*p1_over_gamma[0]*B0[2]; 
  alpha_l[3] = 1.224744871391589*B0[0]*p1_over_gamma[5]-1.224744871391589*B1[0]*p0_over_gamma[5]-0.7071067811865475*B0[0]*p1_over_gamma[1]+0.7071067811865475*B1[0]*p0_over_gamma[1]; 
  alpha_l[4] = 1.224744871391589*B0[0]*p1_over_gamma[6]-1.224744871391589*B1[0]*p0_over_gamma[6]-0.7071067811865475*B0[0]*p1_over_gamma[2]+0.7071067811865475*B1[0]*p0_over_gamma[2]; 
  alpha_l[5] = 1.224744871391589*B0[3]*p1_over_gamma[3]-1.224744871391589*B1[3]*p0_over_gamma[3]+2.0*E2[3]+0.7071067811865475*p0_over_gamma[0]*B1[3]-0.7071067811865475*p1_over_gamma[0]*B0[3]; 
  alpha_l[6] = 1.224744871391589*B0[1]*p1_over_gamma[5]-1.224744871391589*B1[1]*p0_over_gamma[5]-0.7071067811865475*B0[1]*p1_over_gamma[1]+0.7071067811865475*B1[1]*p0_over_gamma[1]; 
  alpha_l[7] = 1.224744871391589*B0[2]*p1_over_gamma[5]-1.224744871391589*B1[2]*p0_over_gamma[5]+0.7071067811865475*p0_over_gamma[1]*B1[2]-0.7071067811865475*p1_over_gamma[1]*B0[2]; 
  alpha_l[8] = 1.224744871391589*B0[1]*p1_over_gamma[6]-1.224744871391589*B1[1]*p0_over_gamma[6]-0.7071067811865475*B0[1]*p1_over_gamma[2]+0.7071067811865475*B1[1]*p0_over_gamma[2]; 
  alpha_l[9] = 1.224744871391589*B0[2]*p1_over_gamma[6]-1.224744871391589*B1[2]*p0_over_gamma[6]-0.7071067811865475*B0[2]*p1_over_gamma[2]+0.7071067811865475*B1[2]*p0_over_gamma[2]; 
  alpha_l[10] = 1.224744871391589*B0[0]*p1_over_gamma[7]-1.224744871391589*B1[0]*p0_over_gamma[7]-0.7071067811865475*B0[0]*p1_over_gamma[4]+0.7071067811865475*B1[0]*p0_over_gamma[4]; 
  alpha_l[11] = 1.224744871391589*B0[3]*p1_over_gamma[5]-1.224744871391589*B1[3]*p0_over_gamma[5]+0.7071067811865475*p0_over_gamma[1]*B1[3]-0.7071067811865475*p1_over_gamma[1]*B0[3]; 
  alpha_l[12] = 1.224744871391589*B0[3]*p1_over_gamma[6]-1.224744871391589*B1[3]*p0_over_gamma[6]+0.7071067811865475*p0_over_gamma[2]*B1[3]-0.7071067811865475*p1_over_gamma[2]*B0[3]; 
  alpha_l[13] = 1.224744871391589*B0[1]*p1_over_gamma[7]-1.224744871391589*B1[1]*p0_over_gamma[7]-0.7071067811865475*B0[1]*p1_over_gamma[4]+0.7071067811865475*B1[1]*p0_over_gamma[4]; 
  alpha_l[14] = 1.224744871391589*B0[2]*p1_over_gamma[7]-1.224744871391589*B1[2]*p0_over_gamma[7]-0.7071067811865475*B0[2]*p1_over_gamma[4]+0.7071067811865475*B1[2]*p0_over_gamma[4]; 
  alpha_l[15] = 1.224744871391589*B0[3]*p1_over_gamma[7]-1.224744871391589*B1[3]*p0_over_gamma[7]-0.7071067811865475*B0[3]*p1_over_gamma[4]+0.7071067811865475*B1[3]*p0_over_gamma[4]; 

  alpha_r[0] = (-1.224744871391589*B0[0]*p1_over_gamma[3])+1.224744871391589*B1[0]*p0_over_gamma[3]-0.7071067811865475*B0[0]*p1_over_gamma[0]+0.7071067811865475*B1[0]*p0_over_gamma[0]+2.0*E2[0]; 
  alpha_r[1] = (-1.224744871391589*B0[1]*p1_over_gamma[3])+1.224744871391589*B1[1]*p0_over_gamma[3]+2.0*E2[1]+0.7071067811865475*p0_over_gamma[0]*B1[1]-0.7071067811865475*p1_over_gamma[0]*B0[1]; 
  alpha_r[2] = (-1.224744871391589*B0[2]*p1_over_gamma[3])+1.224744871391589*B1[2]*p0_over_gamma[3]+2.0*E2[2]+0.7071067811865475*p0_over_gamma[0]*B1[2]-0.7071067811865475*p1_over_gamma[0]*B0[2]; 
  alpha_r[3] = (-1.224744871391589*B0[0]*p1_over_gamma[5])+1.224744871391589*B1[0]*p0_over_gamma[5]-0.7071067811865475*B0[0]*p1_over_gamma[1]+0.7071067811865475*B1[0]*p0_over_gamma[1]; 
  alpha_r[4] = (-1.224744871391589*B0[0]*p1_over_gamma[6])+1.224744871391589*B1[0]*p0_over_gamma[6]-0.7071067811865475*B0[0]*p1_over_gamma[2]+0.7071067811865475*B1[0]*p0_over_gamma[2]; 
  alpha_r[5] = (-1.224744871391589*B0[3]*p1_over_gamma[3])+1.224744871391589*B1[3]*p0_over_gamma[3]+2.0*E2[3]+0.7071067811865475*p0_over_gamma[0]*B1[3]-0.7071067811865475*p1_over_gamma[0]*B0[3]; 
  alpha_r[6] = (-1.224744871391589*B0[1]*p1_over_gamma[5])+1.224744871391589*B1[1]*p0_over_gamma[5]-0.7071067811865475*B0[1]*p1_over_gamma[1]+0.7071067811865475*B1[1]*p0_over_gamma[1]; 
  alpha_r[7] = (-1.224744871391589*B0[2]*p1_over_gamma[5])+1.224744871391589*B1[2]*p0_over_gamma[5]+0.7071067811865475*p0_over_gamma[1]*B1[2]-0.7071067811865475*p1_over_gamma[1]*B0[2]; 
  alpha_r[8] = (-1.224744871391589*B0[1]*p1_over_gamma[6])+1.224744871391589*B1[1]*p0_over_gamma[6]-0.7071067811865475*B0[1]*p1_over_gamma[2]+0.7071067811865475*B1[1]*p0_over_gamma[2]; 
  alpha_r[9] = (-1.224744871391589*B0[2]*p1_over_gamma[6])+1.224744871391589*B1[2]*p0_over_gamma[6]-0.7071067811865475*B0[2]*p1_over_gamma[2]+0.7071067811865475*B1[2]*p0_over_gamma[2]; 
  alpha_r[10] = (-1.224744871391589*B0[0]*p1_over_gamma[7])+1.224744871391589*B1[0]*p0_over_gamma[7]-0.7071067811865475*B0[0]*p1_over_gamma[4]+0.7071067811865475*B1[0]*p0_over_gamma[4]; 
  alpha_r[11] = (-1.224744871391589*B0[3]*p1_over_gamma[5])+1.224744871391589*B1[3]*p0_over_gamma[5]+0.7071067811865475*p0_over_gamma[1]*B1[3]-0.7071067811865475*p1_over_gamma[1]*B0[3]; 
  alpha_r[12] = (-1.224744871391589*B0[3]*p1_over_gamma[6])+1.224744871391589*B1[3]*p0_over_gamma[6]+0.7071067811865475*p0_over_gamma[2]*B1[3]-0.7071067811865475*p1_over_gamma[2]*B0[3]; 
  alpha_r[13] = (-1.224744871391589*B0[1]*p1_over_gamma[7])+1.224744871391589*B1[1]*p0_over_gamma[7]-0.7071067811865475*B0[1]*p1_over_gamma[4]+0.7071067811865475*B1[1]*p0_over_gamma[4]; 
  alpha_r[14] = (-1.224744871391589*B0[2]*p1_over_gamma[7])+1.224744871391589*B1[2]*p0_over_gamma[7]-0.7071067811865475*B0[2]*p1_over_gamma[4]+0.7071067811865475*B1[2]*p0_over_gamma[4]; 
  alpha_r[15] = (-1.224744871391589*B0[3]*p1_over_gamma[7])+1.224744871391589*B1[3]*p0_over_gamma[7]-0.7071067811865475*B0[3]*p1_over_gamma[4]+0.7071067811865475*B1[3]*p0_over_gamma[4]; 

  double fUpwindQuad_l[16] = {0.0};
  double fUpwindQuad_r[16] = {0.0};
  double fUpwind_l[16] = {0.0};;
  double fUpwind_r[16] = {0.0};
  double Ghat_l[16] = {0.0}; 
  double Ghat_r[16] = {0.0}; 

  if (alpha_l[15]-alpha_l[14]-alpha_l[13]-alpha_l[12]-alpha_l[11]+alpha_l[10]+alpha_l[9]+alpha_l[8]+alpha_l[7]+alpha_l[6]+alpha_l[5]-alpha_l[4]-alpha_l[3]-alpha_l[2]-alpha_l[1]+alpha_l[0] > 0) { 
    fUpwindQuad_l[0] = ser_5x_p1_surfx5_eval_quad_node_0_r(fl); 
  } else { 
    fUpwindQuad_l[0] = ser_5x_p1_surfx5_eval_quad_node_0_l(fc); 
  } 
  if (alpha_r[15]-alpha_r[14]-alpha_r[13]-alpha_r[12]-alpha_r[11]+alpha_r[10]+alpha_r[9]+alpha_r[8]+alpha_r[7]+alpha_r[6]+alpha_r[5]-alpha_r[4]-alpha_r[3]-alpha_r[2]-alpha_r[1]+alpha_r[0] > 0) { 
    fUpwindQuad_r[0] = ser_5x_p1_surfx5_eval_quad_node_0_r(fc); 
  } else { 
    fUpwindQuad_r[0] = ser_5x_p1_surfx5_eval_quad_node_0_l(fr); 
  } 
  if ((-alpha_l[15])+alpha_l[14]+alpha_l[13]+alpha_l[12]-alpha_l[11]-alpha_l[10]-alpha_l[9]-alpha_l[8]+alpha_l[7]+alpha_l[6]+alpha_l[5]+alpha_l[4]-alpha_l[3]-alpha_l[2]-alpha_l[1]+alpha_l[0] > 0) { 
    fUpwindQuad_l[1] = ser_5x_p1_surfx5_eval_quad_node_1_r(fl); 
  } else { 
    fUpwindQuad_l[1] = ser_5x_p1_surfx5_eval_quad_node_1_l(fc); 
  } 
  if ((-alpha_r[15])+alpha_r[14]+alpha_r[13]+alpha_r[12]-alpha_r[11]-alpha_r[10]-alpha_r[9]-alpha_r[8]+alpha_r[7]+alpha_r[6]+alpha_r[5]+alpha_r[4]-alpha_r[3]-alpha_r[2]-alpha_r[1]+alpha_r[0] > 0) { 
    fUpwindQuad_r[1] = ser_5x_p1_surfx5_eval_quad_node_1_r(fc); 
  } else { 
    fUpwindQuad_r[1] = ser_5x_p1_surfx5_eval_quad_node_1_l(fr); 
  } 
  if ((-alpha_l[15])+alpha_l[14]+alpha_l[13]-alpha_l[12]+alpha_l[11]-alpha_l[10]+alpha_l[9]+alpha_l[8]-alpha_l[7]-alpha_l[6]+alpha_l[5]-alpha_l[4]+alpha_l[3]-alpha_l[2]-alpha_l[1]+alpha_l[0] > 0) { 
    fUpwindQuad_l[2] = ser_5x_p1_surfx5_eval_quad_node_2_r(fl); 
  } else { 
    fUpwindQuad_l[2] = ser_5x_p1_surfx5_eval_quad_node_2_l(fc); 
  } 
  if ((-alpha_r[15])+alpha_r[14]+alpha_r[13]-alpha_r[12]+alpha_r[11]-alpha_r[10]+alpha_r[9]+alpha_r[8]-alpha_r[7]-alpha_r[6]+alpha_r[5]-alpha_r[4]+alpha_r[3]-alpha_r[2]-alpha_r[1]+alpha_r[0] > 0) { 
    fUpwindQuad_r[2] = ser_5x_p1_surfx5_eval_quad_node_2_r(fc); 
  } else { 
    fUpwindQuad_r[2] = ser_5x_p1_surfx5_eval_quad_node_2_l(fr); 
  } 
  if (alpha_l[15]-alpha_l[14]-alpha_l[13]+alpha_l[12]+alpha_l[11]+alpha_l[10]-alpha_l[9]-alpha_l[8]-alpha_l[7]-alpha_l[6]+alpha_l[5]+alpha_l[4]+alpha_l[3]-alpha_l[2]-alpha_l[1]+alpha_l[0] > 0) { 
    fUpwindQuad_l[3] = ser_5x_p1_surfx5_eval_quad_node_3_r(fl); 
  } else { 
    fUpwindQuad_l[3] = ser_5x_p1_surfx5_eval_quad_node_3_l(fc); 
  } 
  if (alpha_r[15]-alpha_r[14]-alpha_r[13]+alpha_r[12]+alpha_r[11]+alpha_r[10]-alpha_r[9]-alpha_r[8]-alpha_r[7]-alpha_r[6]+alpha_r[5]+alpha_r[4]+alpha_r[3]-alpha_r[2]-alpha_r[1]+alpha_r[0] > 0) { 
    fUpwindQuad_r[3] = ser_5x_p1_surfx5_eval_quad_node_3_r(fc); 
  } else { 
    fUpwindQuad_r[3] = ser_5x_p1_surfx5_eval_quad_node_3_l(fr); 
  } 
  if ((-alpha_l[15])+alpha_l[14]-alpha_l[13]+alpha_l[12]+alpha_l[11]+alpha_l[10]-alpha_l[9]+alpha_l[8]-alpha_l[7]+alpha_l[6]-alpha_l[5]-alpha_l[4]-alpha_l[3]+alpha_l[2]-alpha_l[1]+alpha_l[0] > 0) { 
    fUpwindQuad_l[4] = ser_5x_p1_surfx5_eval_quad_node_4_r(fl); 
  } else { 
    fUpwindQuad_l[4] = ser_5x_p1_surfx5_eval_quad_node_4_l(fc); 
  } 
  if ((-alpha_r[15])+alpha_r[14]-alpha_r[13]+alpha_r[12]+alpha_r[11]+alpha_r[10]-alpha_r[9]+alpha_r[8]-alpha_r[7]+alpha_r[6]-alpha_r[5]-alpha_r[4]-alpha_r[3]+alpha_r[2]-alpha_r[1]+alpha_r[0] > 0) { 
    fUpwindQuad_r[4] = ser_5x_p1_surfx5_eval_quad_node_4_r(fc); 
  } else { 
    fUpwindQuad_r[4] = ser_5x_p1_surfx5_eval_quad_node_4_l(fr); 
  } 
  if (alpha_l[15]-alpha_l[14]+alpha_l[13]-alpha_l[12]+alpha_l[11]-alpha_l[10]+alpha_l[9]-alpha_l[8]-alpha_l[7]+alpha_l[6]-alpha_l[5]+alpha_l[4]-alpha_l[3]+alpha_l[2]-alpha_l[1]+alpha_l[0] > 0) { 
    fUpwindQuad_l[5] = ser_5x_p1_surfx5_eval_quad_node_5_r(fl); 
  } else { 
    fUpwindQuad_l[5] = ser_5x_p1_surfx5_eval_quad_node_5_l(fc); 
  } 
  if (alpha_r[15]-alpha_r[14]+alpha_r[13]-alpha_r[12]+alpha_r[11]-alpha_r[10]+alpha_r[9]-alpha_r[8]-alpha_r[7]+alpha_r[6]-alpha_r[5]+alpha_r[4]-alpha_r[3]+alpha_r[2]-alpha_r[1]+alpha_r[0] > 0) { 
    fUpwindQuad_r[5] = ser_5x_p1_surfx5_eval_quad_node_5_r(fc); 
  } else { 
    fUpwindQuad_r[5] = ser_5x_p1_surfx5_eval_quad_node_5_l(fr); 
  } 
  if (alpha_l[15]-alpha_l[14]+alpha_l[13]+alpha_l[12]-alpha_l[11]-alpha_l[10]-alpha_l[9]+alpha_l[8]+alpha_l[7]-alpha_l[6]-alpha_l[5]-alpha_l[4]+alpha_l[3]+alpha_l[2]-alpha_l[1]+alpha_l[0] > 0) { 
    fUpwindQuad_l[6] = ser_5x_p1_surfx5_eval_quad_node_6_r(fl); 
  } else { 
    fUpwindQuad_l[6] = ser_5x_p1_surfx5_eval_quad_node_6_l(fc); 
  } 
  if (alpha_r[15]-alpha_r[14]+alpha_r[13]+alpha_r[12]-alpha_r[11]-alpha_r[10]-alpha_r[9]+alpha_r[8]+alpha_r[7]-alpha_r[6]-alpha_r[5]-alpha_r[4]+alpha_r[3]+alpha_r[2]-alpha_r[1]+alpha_r[0] > 0) { 
    fUpwindQuad_r[6] = ser_5x_p1_surfx5_eval_quad_node_6_r(fc); 
  } else { 
    fUpwindQuad_r[6] = ser_5x_p1_surfx5_eval_quad_node_6_l(fr); 
  } 
  if ((-alpha_l[15])+alpha_l[14]-alpha_l[13]-alpha_l[12]-alpha_l[11]+alpha_l[10]+alpha_l[9]-alpha_l[8]+alpha_l[7]-alpha_l[6]-alpha_l[5]+alpha_l[4]+alpha_l[3]+alpha_l[2]-alpha_l[1]+alpha_l[0] > 0) { 
    fUpwindQuad_l[7] = ser_5x_p1_surfx5_eval_quad_node_7_r(fl); 
  } else { 
    fUpwindQuad_l[7] = ser_5x_p1_surfx5_eval_quad_node_7_l(fc); 
  } 
  if ((-alpha_r[15])+alpha_r[14]-alpha_r[13]-alpha_r[12]-alpha_r[11]+alpha_r[10]+alpha_r[9]-alpha_r[8]+alpha_r[7]-alpha_r[6]-alpha_r[5]+alpha_r[4]+alpha_r[3]+alpha_r[2]-alpha_r[1]+alpha_r[0] > 0) { 
    fUpwindQuad_r[7] = ser_5x_p1_surfx5_eval_quad_node_7_r(fc); 
  } else { 
    fUpwindQuad_r[7] = ser_5x_p1_surfx5_eval_quad_node_7_l(fr); 
  } 
  if ((-alpha_l[15])-alpha_l[14]+alpha_l[13]+alpha_l[12]+alpha_l[11]+alpha_l[10]+alpha_l[9]-alpha_l[8]+alpha_l[7]-alpha_l[6]-alpha_l[5]-alpha_l[4]-alpha_l[3]-alpha_l[2]+alpha_l[1]+alpha_l[0] > 0) { 
    fUpwindQuad_l[8] = ser_5x_p1_surfx5_eval_quad_node_8_r(fl); 
  } else { 
    fUpwindQuad_l[8] = ser_5x_p1_surfx5_eval_quad_node_8_l(fc); 
  } 
  if ((-alpha_r[15])-alpha_r[14]+alpha_r[13]+alpha_r[12]+alpha_r[11]+alpha_r[10]+alpha_r[9]-alpha_r[8]+alpha_r[7]-alpha_r[6]-alpha_r[5]-alpha_r[4]-alpha_r[3]-alpha_r[2]+alpha_r[1]+alpha_r[0] > 0) { 
    fUpwindQuad_r[8] = ser_5x_p1_surfx5_eval_quad_node_8_r(fc); 
  } else { 
    fUpwindQuad_r[8] = ser_5x_p1_surfx5_eval_quad_node_8_l(fr); 
  } 
  if (alpha_l[15]+alpha_l[14]-alpha_l[13]-alpha_l[12]+alpha_l[11]-alpha_l[10]-alpha_l[9]+alpha_l[8]+alpha_l[7]-alpha_l[6]-alpha_l[5]+alpha_l[4]-alpha_l[3]-alpha_l[2]+alpha_l[1]+alpha_l[0] > 0) { 
    fUpwindQuad_l[9] = ser_5x_p1_surfx5_eval_quad_node_9_r(fl); 
  } else { 
    fUpwindQuad_l[9] = ser_5x_p1_surfx5_eval_quad_node_9_l(fc); 
  } 
  if (alpha_r[15]+alpha_r[14]-alpha_r[13]-alpha_r[12]+alpha_r[11]-alpha_r[10]-alpha_r[9]+alpha_r[8]+alpha_r[7]-alpha_r[6]-alpha_r[5]+alpha_r[4]-alpha_r[3]-alpha_r[2]+alpha_r[1]+alpha_r[0] > 0) { 
    fUpwindQuad_r[9] = ser_5x_p1_surfx5_eval_quad_node_9_r(fc); 
  } else { 
    fUpwindQuad_r[9] = ser_5x_p1_surfx5_eval_quad_node_9_l(fr); 
  } 
  if (alpha_l[15]+alpha_l[14]-alpha_l[13]+alpha_l[12]-alpha_l[11]-alpha_l[10]+alpha_l[9]-alpha_l[8]-alpha_l[7]+alpha_l[6]-alpha_l[5]-alpha_l[4]+alpha_l[3]-alpha_l[2]+alpha_l[1]+alpha_l[0] > 0) { 
    fUpwindQuad_l[10] = ser_5x_p1_surfx5_eval_quad_node_10_r(fl); 
  } else { 
    fUpwindQuad_l[10] = ser_5x_p1_surfx5_eval_quad_node_10_l(fc); 
  } 
  if (alpha_r[15]+alpha_r[14]-alpha_r[13]+alpha_r[12]-alpha_r[11]-alpha_r[10]+alpha_r[9]-alpha_r[8]-alpha_r[7]+alpha_r[6]-alpha_r[5]-alpha_r[4]+alpha_r[3]-alpha_r[2]+alpha_r[1]+alpha_r[0] > 0) { 
    fUpwindQuad_r[10] = ser_5x_p1_surfx5_eval_quad_node_10_r(fc); 
  } else { 
    fUpwindQuad_r[10] = ser_5x_p1_surfx5_eval_quad_node_10_l(fr); 
  } 
  if ((-alpha_l[15])-alpha_l[14]+alpha_l[13]-alpha_l[12]-alpha_l[11]+alpha_l[10]-alpha_l[9]+alpha_l[8]-alpha_l[7]+alpha_l[6]-alpha_l[5]+alpha_l[4]+alpha_l[3]-alpha_l[2]+alpha_l[1]+alpha_l[0] > 0) { 
    fUpwindQuad_l[11] = ser_5x_p1_surfx5_eval_quad_node_11_r(fl); 
  } else { 
    fUpwindQuad_l[11] = ser_5x_p1_surfx5_eval_quad_node_11_l(fc); 
  } 
  if ((-alpha_r[15])-alpha_r[14]+alpha_r[13]-alpha_r[12]-alpha_r[11]+alpha_r[10]-alpha_r[9]+alpha_r[8]-alpha_r[7]+alpha_r[6]-alpha_r[5]+alpha_r[4]+alpha_r[3]-alpha_r[2]+alpha_r[1]+alpha_r[0] > 0) { 
    fUpwindQuad_r[11] = ser_5x_p1_surfx5_eval_quad_node_11_r(fc); 
  } else { 
    fUpwindQuad_r[11] = ser_5x_p1_surfx5_eval_quad_node_11_l(fr); 
  } 
  if (alpha_l[15]+alpha_l[14]+alpha_l[13]-alpha_l[12]-alpha_l[11]+alpha_l[10]-alpha_l[9]-alpha_l[8]-alpha_l[7]-alpha_l[6]+alpha_l[5]-alpha_l[4]-alpha_l[3]+alpha_l[2]+alpha_l[1]+alpha_l[0] > 0) { 
    fUpwindQuad_l[12] = ser_5x_p1_surfx5_eval_quad_node_12_r(fl); 
  } else { 
    fUpwindQuad_l[12] = ser_5x_p1_surfx5_eval_quad_node_12_l(fc); 
  } 
  if (alpha_r[15]+alpha_r[14]+alpha_r[13]-alpha_r[12]-alpha_r[11]+alpha_r[10]-alpha_r[9]-alpha_r[8]-alpha_r[7]-alpha_r[6]+alpha_r[5]-alpha_r[4]-alpha_r[3]+alpha_r[2]+alpha_r[1]+alpha_r[0] > 0) { 
    fUpwindQuad_r[12] = ser_5x_p1_surfx5_eval_quad_node_12_r(fc); 
  } else { 
    fUpwindQuad_r[12] = ser_5x_p1_surfx5_eval_quad_node_12_l(fr); 
  } 
  if ((-alpha_l[15])-alpha_l[14]-alpha_l[13]+alpha_l[12]-alpha_l[11]-alpha_l[10]+alpha_l[9]+alpha_l[8]-alpha_l[7]-alpha_l[6]+alpha_l[5]+alpha_l[4]-alpha_l[3]+alpha_l[2]+alpha_l[1]+alpha_l[0] > 0) { 
    fUpwindQuad_l[13] = ser_5x_p1_surfx5_eval_quad_node_13_r(fl); 
  } else { 
    fUpwindQuad_l[13] = ser_5x_p1_surfx5_eval_quad_node_13_l(fc); 
  } 
  if ((-alpha_r[15])-alpha_r[14]-alpha_r[13]+alpha_r[12]-alpha_r[11]-alpha_r[10]+alpha_r[9]+alpha_r[8]-alpha_r[7]-alpha_r[6]+alpha_r[5]+alpha_r[4]-alpha_r[3]+alpha_r[2]+alpha_r[1]+alpha_r[0] > 0) { 
    fUpwindQuad_r[13] = ser_5x_p1_surfx5_eval_quad_node_13_r(fc); 
  } else { 
    fUpwindQuad_r[13] = ser_5x_p1_surfx5_eval_quad_node_13_l(fr); 
  } 
  if ((-alpha_l[15])-alpha_l[14]-alpha_l[13]-alpha_l[12]+alpha_l[11]-alpha_l[10]-alpha_l[9]-alpha_l[8]+alpha_l[7]+alpha_l[6]+alpha_l[5]-alpha_l[4]+alpha_l[3]+alpha_l[2]+alpha_l[1]+alpha_l[0] > 0) { 
    fUpwindQuad_l[14] = ser_5x_p1_surfx5_eval_quad_node_14_r(fl); 
  } else { 
    fUpwindQuad_l[14] = ser_5x_p1_surfx5_eval_quad_node_14_l(fc); 
  } 
  if ((-alpha_r[15])-alpha_r[14]-alpha_r[13]-alpha_r[12]+alpha_r[11]-alpha_r[10]-alpha_r[9]-alpha_r[8]+alpha_r[7]+alpha_r[6]+alpha_r[5]-alpha_r[4]+alpha_r[3]+alpha_r[2]+alpha_r[1]+alpha_r[0] > 0) { 
    fUpwindQuad_r[14] = ser_5x_p1_surfx5_eval_quad_node_14_r(fc); 
  } else { 
    fUpwindQuad_r[14] = ser_5x_p1_surfx5_eval_quad_node_14_l(fr); 
  } 
  if (alpha_l[15]+alpha_l[14]+alpha_l[13]+alpha_l[12]+alpha_l[11]+alpha_l[10]+alpha_l[9]+alpha_l[8]+alpha_l[7]+alpha_l[6]+alpha_l[5]+alpha_l[4]+alpha_l[3]+alpha_l[2]+alpha_l[1]+alpha_l[0] > 0) { 
    fUpwindQuad_l[15] = ser_5x_p1_surfx5_eval_quad_node_15_r(fl); 
  } else { 
    fUpwindQuad_l[15] = ser_5x_p1_surfx5_eval_quad_node_15_l(fc); 
  } 
  if (alpha_r[15]+alpha_r[14]+alpha_r[13]+alpha_r[12]+alpha_r[11]+alpha_r[10]+alpha_r[9]+alpha_r[8]+alpha_r[7]+alpha_r[6]+alpha_r[5]+alpha_r[4]+alpha_r[3]+alpha_r[2]+alpha_r[1]+alpha_r[0] > 0) { 
    fUpwindQuad_r[15] = ser_5x_p1_surfx5_eval_quad_node_15_r(fc); 
  } else { 
    fUpwindQuad_r[15] = ser_5x_p1_surfx5_eval_quad_node_15_l(fr); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_5x_p1_upwind_quad_to_modal(fUpwindQuad_l, fUpwind_l); 
  ser_5x_p1_upwind_quad_to_modal(fUpwindQuad_r, fUpwind_r); 

  Ghat_l[0] = 0.25*(alpha_l[15]*fUpwind_l[15]+alpha_l[14]*fUpwind_l[14]+alpha_l[13]*fUpwind_l[13]+alpha_l[12]*fUpwind_l[12]+alpha_l[11]*fUpwind_l[11]+alpha_l[10]*fUpwind_l[10]+alpha_l[9]*fUpwind_l[9]+alpha_l[8]*fUpwind_l[8]+alpha_l[7]*fUpwind_l[7]+alpha_l[6]*fUpwind_l[6]+alpha_l[5]*fUpwind_l[5]+alpha_l[4]*fUpwind_l[4]+alpha_l[3]*fUpwind_l[3]+alpha_l[2]*fUpwind_l[2]+alpha_l[1]*fUpwind_l[1]+alpha_l[0]*fUpwind_l[0]); 
  Ghat_l[1] = 0.25*(alpha_l[14]*fUpwind_l[15]+fUpwind_l[14]*alpha_l[15]+alpha_l[10]*fUpwind_l[13]+fUpwind_l[10]*alpha_l[13]+alpha_l[9]*fUpwind_l[12]+fUpwind_l[9]*alpha_l[12]+alpha_l[7]*fUpwind_l[11]+fUpwind_l[7]*alpha_l[11]+alpha_l[4]*fUpwind_l[8]+fUpwind_l[4]*alpha_l[8]+alpha_l[3]*fUpwind_l[6]+fUpwind_l[3]*alpha_l[6]+alpha_l[2]*fUpwind_l[5]+fUpwind_l[2]*alpha_l[5]+alpha_l[0]*fUpwind_l[1]+fUpwind_l[0]*alpha_l[1]); 
  Ghat_l[2] = 0.25*(alpha_l[13]*fUpwind_l[15]+fUpwind_l[13]*alpha_l[15]+alpha_l[10]*fUpwind_l[14]+fUpwind_l[10]*alpha_l[14]+alpha_l[8]*fUpwind_l[12]+fUpwind_l[8]*alpha_l[12]+alpha_l[6]*fUpwind_l[11]+fUpwind_l[6]*alpha_l[11]+alpha_l[4]*fUpwind_l[9]+fUpwind_l[4]*alpha_l[9]+alpha_l[3]*fUpwind_l[7]+fUpwind_l[3]*alpha_l[7]+alpha_l[1]*fUpwind_l[5]+fUpwind_l[1]*alpha_l[5]+alpha_l[0]*fUpwind_l[2]+fUpwind_l[0]*alpha_l[2]); 
  Ghat_l[3] = 0.25*(alpha_l[12]*fUpwind_l[15]+fUpwind_l[12]*alpha_l[15]+alpha_l[9]*fUpwind_l[14]+fUpwind_l[9]*alpha_l[14]+alpha_l[8]*fUpwind_l[13]+fUpwind_l[8]*alpha_l[13]+alpha_l[5]*fUpwind_l[11]+fUpwind_l[5]*alpha_l[11]+alpha_l[4]*fUpwind_l[10]+fUpwind_l[4]*alpha_l[10]+alpha_l[2]*fUpwind_l[7]+fUpwind_l[2]*alpha_l[7]+alpha_l[1]*fUpwind_l[6]+fUpwind_l[1]*alpha_l[6]+alpha_l[0]*fUpwind_l[3]+fUpwind_l[0]*alpha_l[3]); 
  Ghat_l[4] = 0.25*(alpha_l[11]*fUpwind_l[15]+fUpwind_l[11]*alpha_l[15]+alpha_l[7]*fUpwind_l[14]+fUpwind_l[7]*alpha_l[14]+alpha_l[6]*fUpwind_l[13]+fUpwind_l[6]*alpha_l[13]+alpha_l[5]*fUpwind_l[12]+fUpwind_l[5]*alpha_l[12]+alpha_l[3]*fUpwind_l[10]+fUpwind_l[3]*alpha_l[10]+alpha_l[2]*fUpwind_l[9]+fUpwind_l[2]*alpha_l[9]+alpha_l[1]*fUpwind_l[8]+fUpwind_l[1]*alpha_l[8]+alpha_l[0]*fUpwind_l[4]+fUpwind_l[0]*alpha_l[4]); 
  Ghat_l[5] = 0.25*(alpha_l[10]*fUpwind_l[15]+fUpwind_l[10]*alpha_l[15]+alpha_l[13]*fUpwind_l[14]+fUpwind_l[13]*alpha_l[14]+alpha_l[4]*fUpwind_l[12]+fUpwind_l[4]*alpha_l[12]+alpha_l[3]*fUpwind_l[11]+fUpwind_l[3]*alpha_l[11]+alpha_l[8]*fUpwind_l[9]+fUpwind_l[8]*alpha_l[9]+alpha_l[6]*fUpwind_l[7]+fUpwind_l[6]*alpha_l[7]+alpha_l[0]*fUpwind_l[5]+fUpwind_l[0]*alpha_l[5]+alpha_l[1]*fUpwind_l[2]+fUpwind_l[1]*alpha_l[2]); 
  Ghat_l[6] = 0.25*(alpha_l[9]*fUpwind_l[15]+fUpwind_l[9]*alpha_l[15]+alpha_l[12]*fUpwind_l[14]+fUpwind_l[12]*alpha_l[14]+alpha_l[4]*fUpwind_l[13]+fUpwind_l[4]*alpha_l[13]+alpha_l[2]*fUpwind_l[11]+fUpwind_l[2]*alpha_l[11]+alpha_l[8]*fUpwind_l[10]+fUpwind_l[8]*alpha_l[10]+alpha_l[5]*fUpwind_l[7]+fUpwind_l[5]*alpha_l[7]+alpha_l[0]*fUpwind_l[6]+fUpwind_l[0]*alpha_l[6]+alpha_l[1]*fUpwind_l[3]+fUpwind_l[1]*alpha_l[3]); 
  Ghat_l[7] = 0.25*(alpha_l[8]*fUpwind_l[15]+fUpwind_l[8]*alpha_l[15]+alpha_l[4]*fUpwind_l[14]+fUpwind_l[4]*alpha_l[14]+alpha_l[12]*fUpwind_l[13]+fUpwind_l[12]*alpha_l[13]+alpha_l[1]*fUpwind_l[11]+fUpwind_l[1]*alpha_l[11]+alpha_l[9]*fUpwind_l[10]+fUpwind_l[9]*alpha_l[10]+alpha_l[0]*fUpwind_l[7]+fUpwind_l[0]*alpha_l[7]+alpha_l[5]*fUpwind_l[6]+fUpwind_l[5]*alpha_l[6]+alpha_l[2]*fUpwind_l[3]+fUpwind_l[2]*alpha_l[3]); 
  Ghat_l[8] = 0.25*(alpha_l[7]*fUpwind_l[15]+fUpwind_l[7]*alpha_l[15]+alpha_l[11]*fUpwind_l[14]+fUpwind_l[11]*alpha_l[14]+alpha_l[3]*fUpwind_l[13]+fUpwind_l[3]*alpha_l[13]+alpha_l[2]*fUpwind_l[12]+fUpwind_l[2]*alpha_l[12]+alpha_l[6]*fUpwind_l[10]+fUpwind_l[6]*alpha_l[10]+alpha_l[5]*fUpwind_l[9]+fUpwind_l[5]*alpha_l[9]+alpha_l[0]*fUpwind_l[8]+fUpwind_l[0]*alpha_l[8]+alpha_l[1]*fUpwind_l[4]+fUpwind_l[1]*alpha_l[4]); 
  Ghat_l[9] = 0.25*(alpha_l[6]*fUpwind_l[15]+fUpwind_l[6]*alpha_l[15]+alpha_l[3]*fUpwind_l[14]+fUpwind_l[3]*alpha_l[14]+alpha_l[11]*fUpwind_l[13]+fUpwind_l[11]*alpha_l[13]+alpha_l[1]*fUpwind_l[12]+fUpwind_l[1]*alpha_l[12]+alpha_l[7]*fUpwind_l[10]+fUpwind_l[7]*alpha_l[10]+alpha_l[0]*fUpwind_l[9]+fUpwind_l[0]*alpha_l[9]+alpha_l[5]*fUpwind_l[8]+fUpwind_l[5]*alpha_l[8]+alpha_l[2]*fUpwind_l[4]+fUpwind_l[2]*alpha_l[4]); 
  Ghat_l[10] = 0.25*(alpha_l[5]*fUpwind_l[15]+fUpwind_l[5]*alpha_l[15]+alpha_l[2]*fUpwind_l[14]+fUpwind_l[2]*alpha_l[14]+alpha_l[1]*fUpwind_l[13]+fUpwind_l[1]*alpha_l[13]+alpha_l[11]*fUpwind_l[12]+fUpwind_l[11]*alpha_l[12]+alpha_l[0]*fUpwind_l[10]+fUpwind_l[0]*alpha_l[10]+alpha_l[7]*fUpwind_l[9]+fUpwind_l[7]*alpha_l[9]+alpha_l[6]*fUpwind_l[8]+fUpwind_l[6]*alpha_l[8]+alpha_l[3]*fUpwind_l[4]+fUpwind_l[3]*alpha_l[4]); 
  Ghat_l[11] = 0.25*(alpha_l[4]*fUpwind_l[15]+fUpwind_l[4]*alpha_l[15]+alpha_l[8]*fUpwind_l[14]+fUpwind_l[8]*alpha_l[14]+alpha_l[9]*fUpwind_l[13]+fUpwind_l[9]*alpha_l[13]+alpha_l[10]*fUpwind_l[12]+fUpwind_l[10]*alpha_l[12]+alpha_l[0]*fUpwind_l[11]+fUpwind_l[0]*alpha_l[11]+alpha_l[1]*fUpwind_l[7]+fUpwind_l[1]*alpha_l[7]+alpha_l[2]*fUpwind_l[6]+fUpwind_l[2]*alpha_l[6]+alpha_l[3]*fUpwind_l[5]+fUpwind_l[3]*alpha_l[5]); 
  Ghat_l[12] = 0.25*(alpha_l[3]*fUpwind_l[15]+fUpwind_l[3]*alpha_l[15]+alpha_l[6]*fUpwind_l[14]+fUpwind_l[6]*alpha_l[14]+alpha_l[7]*fUpwind_l[13]+fUpwind_l[7]*alpha_l[13]+alpha_l[0]*fUpwind_l[12]+fUpwind_l[0]*alpha_l[12]+alpha_l[10]*fUpwind_l[11]+fUpwind_l[10]*alpha_l[11]+alpha_l[1]*fUpwind_l[9]+fUpwind_l[1]*alpha_l[9]+alpha_l[2]*fUpwind_l[8]+fUpwind_l[2]*alpha_l[8]+alpha_l[4]*fUpwind_l[5]+fUpwind_l[4]*alpha_l[5]); 
  Ghat_l[13] = 0.25*(alpha_l[2]*fUpwind_l[15]+fUpwind_l[2]*alpha_l[15]+alpha_l[5]*fUpwind_l[14]+fUpwind_l[5]*alpha_l[14]+alpha_l[0]*fUpwind_l[13]+fUpwind_l[0]*alpha_l[13]+alpha_l[7]*fUpwind_l[12]+fUpwind_l[7]*alpha_l[12]+alpha_l[9]*fUpwind_l[11]+fUpwind_l[9]*alpha_l[11]+alpha_l[1]*fUpwind_l[10]+fUpwind_l[1]*alpha_l[10]+alpha_l[3]*fUpwind_l[8]+fUpwind_l[3]*alpha_l[8]+alpha_l[4]*fUpwind_l[6]+fUpwind_l[4]*alpha_l[6]); 
  Ghat_l[14] = 0.25*(alpha_l[1]*fUpwind_l[15]+fUpwind_l[1]*alpha_l[15]+alpha_l[0]*fUpwind_l[14]+fUpwind_l[0]*alpha_l[14]+alpha_l[5]*fUpwind_l[13]+fUpwind_l[5]*alpha_l[13]+alpha_l[6]*fUpwind_l[12]+fUpwind_l[6]*alpha_l[12]+alpha_l[8]*fUpwind_l[11]+fUpwind_l[8]*alpha_l[11]+alpha_l[2]*fUpwind_l[10]+fUpwind_l[2]*alpha_l[10]+alpha_l[3]*fUpwind_l[9]+fUpwind_l[3]*alpha_l[9]+alpha_l[4]*fUpwind_l[7]+fUpwind_l[4]*alpha_l[7]); 
  Ghat_l[15] = 0.25*(alpha_l[0]*fUpwind_l[15]+fUpwind_l[0]*alpha_l[15]+alpha_l[1]*fUpwind_l[14]+fUpwind_l[1]*alpha_l[14]+alpha_l[2]*fUpwind_l[13]+fUpwind_l[2]*alpha_l[13]+alpha_l[3]*fUpwind_l[12]+fUpwind_l[3]*alpha_l[12]+alpha_l[4]*fUpwind_l[11]+fUpwind_l[4]*alpha_l[11]+alpha_l[5]*fUpwind_l[10]+fUpwind_l[5]*alpha_l[10]+alpha_l[6]*fUpwind_l[9]+fUpwind_l[6]*alpha_l[9]+alpha_l[7]*fUpwind_l[8]+fUpwind_l[7]*alpha_l[8]); 

  Ghat_r[0] = 0.25*(alpha_r[15]*fUpwind_r[15]+alpha_r[14]*fUpwind_r[14]+alpha_r[13]*fUpwind_r[13]+alpha_r[12]*fUpwind_r[12]+alpha_r[11]*fUpwind_r[11]+alpha_r[10]*fUpwind_r[10]+alpha_r[9]*fUpwind_r[9]+alpha_r[8]*fUpwind_r[8]+alpha_r[7]*fUpwind_r[7]+alpha_r[6]*fUpwind_r[6]+alpha_r[5]*fUpwind_r[5]+alpha_r[4]*fUpwind_r[4]+alpha_r[3]*fUpwind_r[3]+alpha_r[2]*fUpwind_r[2]+alpha_r[1]*fUpwind_r[1]+alpha_r[0]*fUpwind_r[0]); 
  Ghat_r[1] = 0.25*(alpha_r[14]*fUpwind_r[15]+fUpwind_r[14]*alpha_r[15]+alpha_r[10]*fUpwind_r[13]+fUpwind_r[10]*alpha_r[13]+alpha_r[9]*fUpwind_r[12]+fUpwind_r[9]*alpha_r[12]+alpha_r[7]*fUpwind_r[11]+fUpwind_r[7]*alpha_r[11]+alpha_r[4]*fUpwind_r[8]+fUpwind_r[4]*alpha_r[8]+alpha_r[3]*fUpwind_r[6]+fUpwind_r[3]*alpha_r[6]+alpha_r[2]*fUpwind_r[5]+fUpwind_r[2]*alpha_r[5]+alpha_r[0]*fUpwind_r[1]+fUpwind_r[0]*alpha_r[1]); 
  Ghat_r[2] = 0.25*(alpha_r[13]*fUpwind_r[15]+fUpwind_r[13]*alpha_r[15]+alpha_r[10]*fUpwind_r[14]+fUpwind_r[10]*alpha_r[14]+alpha_r[8]*fUpwind_r[12]+fUpwind_r[8]*alpha_r[12]+alpha_r[6]*fUpwind_r[11]+fUpwind_r[6]*alpha_r[11]+alpha_r[4]*fUpwind_r[9]+fUpwind_r[4]*alpha_r[9]+alpha_r[3]*fUpwind_r[7]+fUpwind_r[3]*alpha_r[7]+alpha_r[1]*fUpwind_r[5]+fUpwind_r[1]*alpha_r[5]+alpha_r[0]*fUpwind_r[2]+fUpwind_r[0]*alpha_r[2]); 
  Ghat_r[3] = 0.25*(alpha_r[12]*fUpwind_r[15]+fUpwind_r[12]*alpha_r[15]+alpha_r[9]*fUpwind_r[14]+fUpwind_r[9]*alpha_r[14]+alpha_r[8]*fUpwind_r[13]+fUpwind_r[8]*alpha_r[13]+alpha_r[5]*fUpwind_r[11]+fUpwind_r[5]*alpha_r[11]+alpha_r[4]*fUpwind_r[10]+fUpwind_r[4]*alpha_r[10]+alpha_r[2]*fUpwind_r[7]+fUpwind_r[2]*alpha_r[7]+alpha_r[1]*fUpwind_r[6]+fUpwind_r[1]*alpha_r[6]+alpha_r[0]*fUpwind_r[3]+fUpwind_r[0]*alpha_r[3]); 
  Ghat_r[4] = 0.25*(alpha_r[11]*fUpwind_r[15]+fUpwind_r[11]*alpha_r[15]+alpha_r[7]*fUpwind_r[14]+fUpwind_r[7]*alpha_r[14]+alpha_r[6]*fUpwind_r[13]+fUpwind_r[6]*alpha_r[13]+alpha_r[5]*fUpwind_r[12]+fUpwind_r[5]*alpha_r[12]+alpha_r[3]*fUpwind_r[10]+fUpwind_r[3]*alpha_r[10]+alpha_r[2]*fUpwind_r[9]+fUpwind_r[2]*alpha_r[9]+alpha_r[1]*fUpwind_r[8]+fUpwind_r[1]*alpha_r[8]+alpha_r[0]*fUpwind_r[4]+fUpwind_r[0]*alpha_r[4]); 
  Ghat_r[5] = 0.25*(alpha_r[10]*fUpwind_r[15]+fUpwind_r[10]*alpha_r[15]+alpha_r[13]*fUpwind_r[14]+fUpwind_r[13]*alpha_r[14]+alpha_r[4]*fUpwind_r[12]+fUpwind_r[4]*alpha_r[12]+alpha_r[3]*fUpwind_r[11]+fUpwind_r[3]*alpha_r[11]+alpha_r[8]*fUpwind_r[9]+fUpwind_r[8]*alpha_r[9]+alpha_r[6]*fUpwind_r[7]+fUpwind_r[6]*alpha_r[7]+alpha_r[0]*fUpwind_r[5]+fUpwind_r[0]*alpha_r[5]+alpha_r[1]*fUpwind_r[2]+fUpwind_r[1]*alpha_r[2]); 
  Ghat_r[6] = 0.25*(alpha_r[9]*fUpwind_r[15]+fUpwind_r[9]*alpha_r[15]+alpha_r[12]*fUpwind_r[14]+fUpwind_r[12]*alpha_r[14]+alpha_r[4]*fUpwind_r[13]+fUpwind_r[4]*alpha_r[13]+alpha_r[2]*fUpwind_r[11]+fUpwind_r[2]*alpha_r[11]+alpha_r[8]*fUpwind_r[10]+fUpwind_r[8]*alpha_r[10]+alpha_r[5]*fUpwind_r[7]+fUpwind_r[5]*alpha_r[7]+alpha_r[0]*fUpwind_r[6]+fUpwind_r[0]*alpha_r[6]+alpha_r[1]*fUpwind_r[3]+fUpwind_r[1]*alpha_r[3]); 
  Ghat_r[7] = 0.25*(alpha_r[8]*fUpwind_r[15]+fUpwind_r[8]*alpha_r[15]+alpha_r[4]*fUpwind_r[14]+fUpwind_r[4]*alpha_r[14]+alpha_r[12]*fUpwind_r[13]+fUpwind_r[12]*alpha_r[13]+alpha_r[1]*fUpwind_r[11]+fUpwind_r[1]*alpha_r[11]+alpha_r[9]*fUpwind_r[10]+fUpwind_r[9]*alpha_r[10]+alpha_r[0]*fUpwind_r[7]+fUpwind_r[0]*alpha_r[7]+alpha_r[5]*fUpwind_r[6]+fUpwind_r[5]*alpha_r[6]+alpha_r[2]*fUpwind_r[3]+fUpwind_r[2]*alpha_r[3]); 
  Ghat_r[8] = 0.25*(alpha_r[7]*fUpwind_r[15]+fUpwind_r[7]*alpha_r[15]+alpha_r[11]*fUpwind_r[14]+fUpwind_r[11]*alpha_r[14]+alpha_r[3]*fUpwind_r[13]+fUpwind_r[3]*alpha_r[13]+alpha_r[2]*fUpwind_r[12]+fUpwind_r[2]*alpha_r[12]+alpha_r[6]*fUpwind_r[10]+fUpwind_r[6]*alpha_r[10]+alpha_r[5]*fUpwind_r[9]+fUpwind_r[5]*alpha_r[9]+alpha_r[0]*fUpwind_r[8]+fUpwind_r[0]*alpha_r[8]+alpha_r[1]*fUpwind_r[4]+fUpwind_r[1]*alpha_r[4]); 
  Ghat_r[9] = 0.25*(alpha_r[6]*fUpwind_r[15]+fUpwind_r[6]*alpha_r[15]+alpha_r[3]*fUpwind_r[14]+fUpwind_r[3]*alpha_r[14]+alpha_r[11]*fUpwind_r[13]+fUpwind_r[11]*alpha_r[13]+alpha_r[1]*fUpwind_r[12]+fUpwind_r[1]*alpha_r[12]+alpha_r[7]*fUpwind_r[10]+fUpwind_r[7]*alpha_r[10]+alpha_r[0]*fUpwind_r[9]+fUpwind_r[0]*alpha_r[9]+alpha_r[5]*fUpwind_r[8]+fUpwind_r[5]*alpha_r[8]+alpha_r[2]*fUpwind_r[4]+fUpwind_r[2]*alpha_r[4]); 
  Ghat_r[10] = 0.25*(alpha_r[5]*fUpwind_r[15]+fUpwind_r[5]*alpha_r[15]+alpha_r[2]*fUpwind_r[14]+fUpwind_r[2]*alpha_r[14]+alpha_r[1]*fUpwind_r[13]+fUpwind_r[1]*alpha_r[13]+alpha_r[11]*fUpwind_r[12]+fUpwind_r[11]*alpha_r[12]+alpha_r[0]*fUpwind_r[10]+fUpwind_r[0]*alpha_r[10]+alpha_r[7]*fUpwind_r[9]+fUpwind_r[7]*alpha_r[9]+alpha_r[6]*fUpwind_r[8]+fUpwind_r[6]*alpha_r[8]+alpha_r[3]*fUpwind_r[4]+fUpwind_r[3]*alpha_r[4]); 
  Ghat_r[11] = 0.25*(alpha_r[4]*fUpwind_r[15]+fUpwind_r[4]*alpha_r[15]+alpha_r[8]*fUpwind_r[14]+fUpwind_r[8]*alpha_r[14]+alpha_r[9]*fUpwind_r[13]+fUpwind_r[9]*alpha_r[13]+alpha_r[10]*fUpwind_r[12]+fUpwind_r[10]*alpha_r[12]+alpha_r[0]*fUpwind_r[11]+fUpwind_r[0]*alpha_r[11]+alpha_r[1]*fUpwind_r[7]+fUpwind_r[1]*alpha_r[7]+alpha_r[2]*fUpwind_r[6]+fUpwind_r[2]*alpha_r[6]+alpha_r[3]*fUpwind_r[5]+fUpwind_r[3]*alpha_r[5]); 
  Ghat_r[12] = 0.25*(alpha_r[3]*fUpwind_r[15]+fUpwind_r[3]*alpha_r[15]+alpha_r[6]*fUpwind_r[14]+fUpwind_r[6]*alpha_r[14]+alpha_r[7]*fUpwind_r[13]+fUpwind_r[7]*alpha_r[13]+alpha_r[0]*fUpwind_r[12]+fUpwind_r[0]*alpha_r[12]+alpha_r[10]*fUpwind_r[11]+fUpwind_r[10]*alpha_r[11]+alpha_r[1]*fUpwind_r[9]+fUpwind_r[1]*alpha_r[9]+alpha_r[2]*fUpwind_r[8]+fUpwind_r[2]*alpha_r[8]+alpha_r[4]*fUpwind_r[5]+fUpwind_r[4]*alpha_r[5]); 
  Ghat_r[13] = 0.25*(alpha_r[2]*fUpwind_r[15]+fUpwind_r[2]*alpha_r[15]+alpha_r[5]*fUpwind_r[14]+fUpwind_r[5]*alpha_r[14]+alpha_r[0]*fUpwind_r[13]+fUpwind_r[0]*alpha_r[13]+alpha_r[7]*fUpwind_r[12]+fUpwind_r[7]*alpha_r[12]+alpha_r[9]*fUpwind_r[11]+fUpwind_r[9]*alpha_r[11]+alpha_r[1]*fUpwind_r[10]+fUpwind_r[1]*alpha_r[10]+alpha_r[3]*fUpwind_r[8]+fUpwind_r[3]*alpha_r[8]+alpha_r[4]*fUpwind_r[6]+fUpwind_r[4]*alpha_r[6]); 
  Ghat_r[14] = 0.25*(alpha_r[1]*fUpwind_r[15]+fUpwind_r[1]*alpha_r[15]+alpha_r[0]*fUpwind_r[14]+fUpwind_r[0]*alpha_r[14]+alpha_r[5]*fUpwind_r[13]+fUpwind_r[5]*alpha_r[13]+alpha_r[6]*fUpwind_r[12]+fUpwind_r[6]*alpha_r[12]+alpha_r[8]*fUpwind_r[11]+fUpwind_r[8]*alpha_r[11]+alpha_r[2]*fUpwind_r[10]+fUpwind_r[2]*alpha_r[10]+alpha_r[3]*fUpwind_r[9]+fUpwind_r[3]*alpha_r[9]+alpha_r[4]*fUpwind_r[7]+fUpwind_r[4]*alpha_r[7]); 
  Ghat_r[15] = 0.25*(alpha_r[0]*fUpwind_r[15]+fUpwind_r[0]*alpha_r[15]+alpha_r[1]*fUpwind_r[14]+fUpwind_r[1]*alpha_r[14]+alpha_r[2]*fUpwind_r[13]+fUpwind_r[2]*alpha_r[13]+alpha_r[3]*fUpwind_r[12]+fUpwind_r[3]*alpha_r[12]+alpha_r[4]*fUpwind_r[11]+fUpwind_r[4]*alpha_r[11]+alpha_r[5]*fUpwind_r[10]+fUpwind_r[5]*alpha_r[10]+alpha_r[6]*fUpwind_r[9]+fUpwind_r[6]*alpha_r[9]+alpha_r[7]*fUpwind_r[8]+fUpwind_r[7]*alpha_r[8]); 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv12; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv12; 
  out[2] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dv12; 
  out[3] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dv12; 
  out[4] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*dv12; 
  out[5] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv12; 
  out[6] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*dv12; 
  out[7] += (0.7071067811865475*Ghat_l[6]-0.7071067811865475*Ghat_r[6])*dv12; 
  out[8] += (0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*dv12; 
  out[9] += (0.7071067811865475*Ghat_l[8]-0.7071067811865475*Ghat_r[8])*dv12; 
  out[10] += (0.7071067811865475*Ghat_l[9]-0.7071067811865475*Ghat_r[9])*dv12; 
  out[11] += (0.7071067811865475*Ghat_l[10]-0.7071067811865475*Ghat_r[10])*dv12; 
  out[12] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv12; 
  out[13] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dv12; 
  out[14] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dv12; 
  out[15] += -1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dv12; 
  out[16] += (0.7071067811865475*Ghat_l[11]-0.7071067811865475*Ghat_r[11])*dv12; 
  out[17] += (0.7071067811865475*Ghat_l[12]-0.7071067811865475*Ghat_r[12])*dv12; 
  out[18] += (0.7071067811865475*Ghat_l[13]-0.7071067811865475*Ghat_r[13])*dv12; 
  out[19] += (0.7071067811865475*Ghat_l[14]-0.7071067811865475*Ghat_r[14])*dv12; 
  out[20] += -1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dv12; 
  out[21] += -1.224744871391589*(Ghat_r[6]+Ghat_l[6])*dv12; 
  out[22] += -1.224744871391589*(Ghat_r[7]+Ghat_l[7])*dv12; 
  out[23] += -1.224744871391589*(Ghat_r[8]+Ghat_l[8])*dv12; 
  out[24] += -1.224744871391589*(Ghat_r[9]+Ghat_l[9])*dv12; 
  out[25] += -1.224744871391589*(Ghat_r[10]+Ghat_l[10])*dv12; 
  out[26] += (0.7071067811865475*Ghat_l[15]-0.7071067811865475*Ghat_r[15])*dv12; 
  out[27] += -1.224744871391589*(Ghat_r[11]+Ghat_l[11])*dv12; 
  out[28] += -1.224744871391589*(Ghat_r[12]+Ghat_l[12])*dv12; 
  out[29] += -1.224744871391589*(Ghat_r[13]+Ghat_l[13])*dv12; 
  out[30] += -1.224744871391589*(Ghat_r[14]+Ghat_l[14])*dv12; 
  out[31] += -1.224744871391589*(Ghat_r[15]+Ghat_l[15])*dv12; 

  return 0.;

} 
