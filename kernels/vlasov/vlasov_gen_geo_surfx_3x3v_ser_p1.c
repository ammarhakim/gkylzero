#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_hyb_3x3v_p1_surfx1_eval_quad.h> 
#include <gkyl_basis_hyb_3x3v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH void vlasov_gen_geo_surfx_3x3v_ser_p1(const double *w, const double *dxv, const double *alpha_geo, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // alpha_geo:  General geometry alpha.
  // fl/fc/fr:  Input Distribution function in left/center/right cells.
  // out:       Incremented distribution function in center cell.

  const double dx10 = 2/dxv[0]; 
  const double dv = dxv[3], wv = w[3]; 
  double Ghat_r[80] = {0.0}; 
  double Ghat_l[80] = {0.0}; 

  const double *ag0 = &alpha_geo[0]; 

  double alphal[32] = {0.0}; 
  double alphar[32] = {0.0}; 

  alphal[0] = 0.7071067811865475*ag0[0]-1.224744871391589*ag0[1]; 
  alphal[1] = 0.7071067811865475*ag0[2]-1.224744871391589*ag0[7]; 
  alphal[2] = 0.7071067811865475*ag0[3]-1.224744871391589*ag0[8]; 
  alphal[3] = 0.7071067811865475*ag0[4]-1.224744871391589*ag0[10]; 
  alphal[4] = 0.7071067811865475*ag0[5]-1.224744871391589*ag0[13]; 
  alphal[5] = 0.7071067811865475*ag0[6]-1.224744871391589*ag0[17]; 
  alphal[6] = 0.7071067811865475*ag0[9]-1.224744871391589*ag0[22]; 
  alphal[7] = 0.7071067811865475*ag0[11]-1.224744871391589*ag0[23]; 
  alphal[8] = 0.7071067811865475*ag0[12]-1.224744871391589*ag0[24]; 
  alphal[9] = 0.7071067811865475*ag0[14]-1.224744871391589*ag0[26]; 
  alphal[10] = 0.7071067811865475*ag0[15]-1.224744871391589*ag0[27]; 
  alphal[12] = 0.7071067811865475*ag0[18]-1.224744871391589*ag0[32]; 
  alphal[13] = 0.7071067811865475*ag0[19]-1.224744871391589*ag0[33]; 
  alphal[16] = 0.7071067811865475*ag0[25]-1.224744871391589*ag0[42]; 
  alphal[17] = 0.7071067811865475*ag0[28]-1.224744871391589*ag0[43]; 
  alphal[20] = 0.7071067811865475*ag0[34]-1.224744871391589*ag0[47]; 
  alphar[0] = 1.224744871391589*ag0[1]+0.7071067811865475*ag0[0]; 
  alphar[1] = 1.224744871391589*ag0[7]+0.7071067811865475*ag0[2]; 
  alphar[2] = 1.224744871391589*ag0[8]+0.7071067811865475*ag0[3]; 
  alphar[3] = 1.224744871391589*ag0[10]+0.7071067811865475*ag0[4]; 
  alphar[4] = 1.224744871391589*ag0[13]+0.7071067811865475*ag0[5]; 
  alphar[5] = 1.224744871391589*ag0[17]+0.7071067811865475*ag0[6]; 
  alphar[6] = 1.224744871391589*ag0[22]+0.7071067811865475*ag0[9]; 
  alphar[7] = 1.224744871391589*ag0[23]+0.7071067811865475*ag0[11]; 
  alphar[8] = 1.224744871391589*ag0[24]+0.7071067811865475*ag0[12]; 
  alphar[9] = 1.224744871391589*ag0[26]+0.7071067811865475*ag0[14]; 
  alphar[10] = 1.224744871391589*ag0[27]+0.7071067811865475*ag0[15]; 
  alphar[12] = 1.224744871391589*ag0[32]+0.7071067811865475*ag0[18]; 
  alphar[13] = 1.224744871391589*ag0[33]+0.7071067811865475*ag0[19]; 
  alphar[16] = 1.224744871391589*ag0[42]+0.7071067811865475*ag0[25]; 
  alphar[17] = 1.224744871391589*ag0[43]+0.7071067811865475*ag0[28]; 
  alphar[20] = 1.224744871391589*ag0[47]+0.7071067811865475*ag0[34]; 

  double fUpwindQuad_l[108] = {0.0};
  double fUpwindQuad_r[108] = {0.0};
  double fUpwind_l[80] = {0.0};;
  double fUpwind_r[80] = {0.0};

  if ((-0.2371708245126284*(alphal[20]+alphal[17]+alphal[16]))+0.2371708245126284*(alphal[13]+alphal[12]+alphal[10]+alphal[9]+alphal[8]+alphal[7])+0.1767766952966368*alphal[6]-0.2371708245126284*(alphal[5]+alphal[4]+alphal[3])-0.1767766952966368*(alphal[2]+alphal[1])+0.1767766952966368*alphal[0] > 0) { 
    fUpwindQuad_l[0] = hyb_3x3v_p1_surfx1_eval_quad_node_0_r(fl); 
  } else { 
    fUpwindQuad_l[0] = hyb_3x3v_p1_surfx1_eval_quad_node_0_l(fc); 
  } 
  if ((-0.2371708245126284*(alphar[20]+alphar[17]+alphar[16]))+0.2371708245126284*(alphar[13]+alphar[12]+alphar[10]+alphar[9]+alphar[8]+alphar[7])+0.1767766952966368*alphar[6]-0.2371708245126284*(alphar[5]+alphar[4]+alphar[3])-0.1767766952966368*(alphar[2]+alphar[1])+0.1767766952966368*alphar[0] > 0) { 
    fUpwindQuad_r[0] = hyb_3x3v_p1_surfx1_eval_quad_node_0_r(fc); 
  } else { 
    fUpwindQuad_r[0] = hyb_3x3v_p1_surfx1_eval_quad_node_0_l(fr); 
  } 
  if ((-0.2371708245126284*(alphal[17]+alphal[16]))+0.2371708245126284*(alphal[10]+alphal[9]+alphal[8]+alphal[7])+0.1767766952966368*alphal[6]-0.2371708245126284*(alphal[4]+alphal[3])-0.1767766952966368*(alphal[2]+alphal[1])+0.1767766952966368*alphal[0] > 0) { 
    fUpwindQuad_l[1] = hyb_3x3v_p1_surfx1_eval_quad_node_1_r(fl); 
  } else { 
    fUpwindQuad_l[1] = hyb_3x3v_p1_surfx1_eval_quad_node_1_l(fc); 
  } 
  if ((-0.2371708245126284*(alphar[17]+alphar[16]))+0.2371708245126284*(alphar[10]+alphar[9]+alphar[8]+alphar[7])+0.1767766952966368*alphar[6]-0.2371708245126284*(alphar[4]+alphar[3])-0.1767766952966368*(alphar[2]+alphar[1])+0.1767766952966368*alphar[0] > 0) { 
    fUpwindQuad_r[1] = hyb_3x3v_p1_surfx1_eval_quad_node_1_r(fc); 
  } else { 
    fUpwindQuad_r[1] = hyb_3x3v_p1_surfx1_eval_quad_node_1_l(fr); 
  } 
  if (0.2371708245126284*alphal[20]-0.2371708245126284*(alphal[17]+alphal[16]+alphal[13]+alphal[12])+0.2371708245126284*(alphal[10]+alphal[9]+alphal[8]+alphal[7])+0.1767766952966368*alphal[6]+0.2371708245126284*alphal[5]-0.2371708245126284*(alphal[4]+alphal[3])-0.1767766952966368*(alphal[2]+alphal[1])+0.1767766952966368*alphal[0] > 0) { 
    fUpwindQuad_l[2] = hyb_3x3v_p1_surfx1_eval_quad_node_2_r(fl); 
  } else { 
    fUpwindQuad_l[2] = hyb_3x3v_p1_surfx1_eval_quad_node_2_l(fc); 
  } 
  if (0.2371708245126284*alphar[20]-0.2371708245126284*(alphar[17]+alphar[16]+alphar[13]+alphar[12])+0.2371708245126284*(alphar[10]+alphar[9]+alphar[8]+alphar[7])+0.1767766952966368*alphar[6]+0.2371708245126284*alphar[5]-0.2371708245126284*(alphar[4]+alphar[3])-0.1767766952966368*(alphar[2]+alphar[1])+0.1767766952966368*alphar[0] > 0) { 
    fUpwindQuad_r[2] = hyb_3x3v_p1_surfx1_eval_quad_node_2_r(fc); 
  } else { 
    fUpwindQuad_r[2] = hyb_3x3v_p1_surfx1_eval_quad_node_2_l(fr); 
  } 
  if ((-0.2371708245126284*(alphal[20]+alphal[16]))+0.2371708245126284*(alphal[13]+alphal[12]+alphal[8]+alphal[7])+0.1767766952966368*alphal[6]-0.2371708245126284*(alphal[5]+alphal[3])-0.1767766952966368*(alphal[2]+alphal[1])+0.1767766952966368*alphal[0] > 0) { 
    fUpwindQuad_l[3] = hyb_3x3v_p1_surfx1_eval_quad_node_3_r(fl); 
  } else { 
    fUpwindQuad_l[3] = hyb_3x3v_p1_surfx1_eval_quad_node_3_l(fc); 
  } 
  if ((-0.2371708245126284*(alphar[20]+alphar[16]))+0.2371708245126284*(alphar[13]+alphar[12]+alphar[8]+alphar[7])+0.1767766952966368*alphar[6]-0.2371708245126284*(alphar[5]+alphar[3])-0.1767766952966368*(alphar[2]+alphar[1])+0.1767766952966368*alphar[0] > 0) { 
    fUpwindQuad_r[3] = hyb_3x3v_p1_surfx1_eval_quad_node_3_r(fc); 
  } else { 
    fUpwindQuad_r[3] = hyb_3x3v_p1_surfx1_eval_quad_node_3_l(fr); 
  } 
  if ((-0.2371708245126284*alphal[16])+0.2371708245126284*(alphal[8]+alphal[7])+0.1767766952966368*alphal[6]-0.2371708245126284*alphal[3]-0.1767766952966368*(alphal[2]+alphal[1])+0.1767766952966368*alphal[0] > 0) { 
    fUpwindQuad_l[4] = hyb_3x3v_p1_surfx1_eval_quad_node_4_r(fl); 
  } else { 
    fUpwindQuad_l[4] = hyb_3x3v_p1_surfx1_eval_quad_node_4_l(fc); 
  } 
  if ((-0.2371708245126284*alphar[16])+0.2371708245126284*(alphar[8]+alphar[7])+0.1767766952966368*alphar[6]-0.2371708245126284*alphar[3]-0.1767766952966368*(alphar[2]+alphar[1])+0.1767766952966368*alphar[0] > 0) { 
    fUpwindQuad_r[4] = hyb_3x3v_p1_surfx1_eval_quad_node_4_r(fc); 
  } else { 
    fUpwindQuad_r[4] = hyb_3x3v_p1_surfx1_eval_quad_node_4_l(fr); 
  } 
  if (0.2371708245126284*alphal[20]-0.2371708245126284*(alphal[16]+alphal[13]+alphal[12])+0.2371708245126284*(alphal[8]+alphal[7])+0.1767766952966368*alphal[6]+0.2371708245126284*alphal[5]-0.2371708245126284*alphal[3]-0.1767766952966368*(alphal[2]+alphal[1])+0.1767766952966368*alphal[0] > 0) { 
    fUpwindQuad_l[5] = hyb_3x3v_p1_surfx1_eval_quad_node_5_r(fl); 
  } else { 
    fUpwindQuad_l[5] = hyb_3x3v_p1_surfx1_eval_quad_node_5_l(fc); 
  } 
  if (0.2371708245126284*alphar[20]-0.2371708245126284*(alphar[16]+alphar[13]+alphar[12])+0.2371708245126284*(alphar[8]+alphar[7])+0.1767766952966368*alphar[6]+0.2371708245126284*alphar[5]-0.2371708245126284*alphar[3]-0.1767766952966368*(alphar[2]+alphar[1])+0.1767766952966368*alphar[0] > 0) { 
    fUpwindQuad_r[5] = hyb_3x3v_p1_surfx1_eval_quad_node_5_r(fc); 
  } else { 
    fUpwindQuad_r[5] = hyb_3x3v_p1_surfx1_eval_quad_node_5_l(fr); 
  } 
  if ((-0.2371708245126284*alphal[20])+0.2371708245126284*alphal[17]-0.2371708245126284*alphal[16]+0.2371708245126284*(alphal[13]+alphal[12])-0.2371708245126284*(alphal[10]+alphal[9])+0.2371708245126284*(alphal[8]+alphal[7])+0.1767766952966368*alphal[6]-0.2371708245126284*alphal[5]+0.2371708245126284*alphal[4]-0.2371708245126284*alphal[3]-0.1767766952966368*(alphal[2]+alphal[1])+0.1767766952966368*alphal[0] > 0) { 
    fUpwindQuad_l[6] = hyb_3x3v_p1_surfx1_eval_quad_node_6_r(fl); 
  } else { 
    fUpwindQuad_l[6] = hyb_3x3v_p1_surfx1_eval_quad_node_6_l(fc); 
  } 
  if ((-0.2371708245126284*alphar[20])+0.2371708245126284*alphar[17]-0.2371708245126284*alphar[16]+0.2371708245126284*(alphar[13]+alphar[12])-0.2371708245126284*(alphar[10]+alphar[9])+0.2371708245126284*(alphar[8]+alphar[7])+0.1767766952966368*alphar[6]-0.2371708245126284*alphar[5]+0.2371708245126284*alphar[4]-0.2371708245126284*alphar[3]-0.1767766952966368*(alphar[2]+alphar[1])+0.1767766952966368*alphar[0] > 0) { 
    fUpwindQuad_r[6] = hyb_3x3v_p1_surfx1_eval_quad_node_6_r(fc); 
  } else { 
    fUpwindQuad_r[6] = hyb_3x3v_p1_surfx1_eval_quad_node_6_l(fr); 
  } 
  if (0.2371708245126284*alphal[17]-0.2371708245126284*(alphal[16]+alphal[10]+alphal[9])+0.2371708245126284*(alphal[8]+alphal[7])+0.1767766952966368*alphal[6]+0.2371708245126284*alphal[4]-0.2371708245126284*alphal[3]-0.1767766952966368*(alphal[2]+alphal[1])+0.1767766952966368*alphal[0] > 0) { 
    fUpwindQuad_l[7] = hyb_3x3v_p1_surfx1_eval_quad_node_7_r(fl); 
  } else { 
    fUpwindQuad_l[7] = hyb_3x3v_p1_surfx1_eval_quad_node_7_l(fc); 
  } 
  if (0.2371708245126284*alphar[17]-0.2371708245126284*(alphar[16]+alphar[10]+alphar[9])+0.2371708245126284*(alphar[8]+alphar[7])+0.1767766952966368*alphar[6]+0.2371708245126284*alphar[4]-0.2371708245126284*alphar[3]-0.1767766952966368*(alphar[2]+alphar[1])+0.1767766952966368*alphar[0] > 0) { 
    fUpwindQuad_r[7] = hyb_3x3v_p1_surfx1_eval_quad_node_7_r(fc); 
  } else { 
    fUpwindQuad_r[7] = hyb_3x3v_p1_surfx1_eval_quad_node_7_l(fr); 
  } 
  if (0.2371708245126284*(alphal[20]+alphal[17])-0.2371708245126284*(alphal[16]+alphal[13]+alphal[12]+alphal[10]+alphal[9])+0.2371708245126284*(alphal[8]+alphal[7])+0.1767766952966368*alphal[6]+0.2371708245126284*(alphal[5]+alphal[4])-0.2371708245126284*alphal[3]-0.1767766952966368*(alphal[2]+alphal[1])+0.1767766952966368*alphal[0] > 0) { 
    fUpwindQuad_l[8] = hyb_3x3v_p1_surfx1_eval_quad_node_8_r(fl); 
  } else { 
    fUpwindQuad_l[8] = hyb_3x3v_p1_surfx1_eval_quad_node_8_l(fc); 
  } 
  if (0.2371708245126284*(alphar[20]+alphar[17])-0.2371708245126284*(alphar[16]+alphar[13]+alphar[12]+alphar[10]+alphar[9])+0.2371708245126284*(alphar[8]+alphar[7])+0.1767766952966368*alphar[6]+0.2371708245126284*(alphar[5]+alphar[4])-0.2371708245126284*alphar[3]-0.1767766952966368*(alphar[2]+alphar[1])+0.1767766952966368*alphar[0] > 0) { 
    fUpwindQuad_r[8] = hyb_3x3v_p1_surfx1_eval_quad_node_8_r(fc); 
  } else { 
    fUpwindQuad_r[8] = hyb_3x3v_p1_surfx1_eval_quad_node_8_l(fr); 
  } 
  if ((-0.2371708245126284*(alphal[20]+alphal[17]))+0.2371708245126284*(alphal[13]+alphal[12]+alphal[10]+alphal[9])+0.1767766952966368*alphal[6]-0.2371708245126284*(alphal[5]+alphal[4])-0.1767766952966368*(alphal[2]+alphal[1])+0.1767766952966368*alphal[0] > 0) { 
    fUpwindQuad_l[9] = hyb_3x3v_p1_surfx1_eval_quad_node_9_r(fl); 
  } else { 
    fUpwindQuad_l[9] = hyb_3x3v_p1_surfx1_eval_quad_node_9_l(fc); 
  } 
  if ((-0.2371708245126284*(alphar[20]+alphar[17]))+0.2371708245126284*(alphar[13]+alphar[12]+alphar[10]+alphar[9])+0.1767766952966368*alphar[6]-0.2371708245126284*(alphar[5]+alphar[4])-0.1767766952966368*(alphar[2]+alphar[1])+0.1767766952966368*alphar[0] > 0) { 
    fUpwindQuad_r[9] = hyb_3x3v_p1_surfx1_eval_quad_node_9_r(fc); 
  } else { 
    fUpwindQuad_r[9] = hyb_3x3v_p1_surfx1_eval_quad_node_9_l(fr); 
  } 
  if ((-0.2371708245126284*alphal[17])+0.2371708245126284*(alphal[10]+alphal[9])+0.1767766952966368*alphal[6]-0.2371708245126284*alphal[4]-0.1767766952966368*(alphal[2]+alphal[1])+0.1767766952966368*alphal[0] > 0) { 
    fUpwindQuad_l[10] = hyb_3x3v_p1_surfx1_eval_quad_node_10_r(fl); 
  } else { 
    fUpwindQuad_l[10] = hyb_3x3v_p1_surfx1_eval_quad_node_10_l(fc); 
  } 
  if ((-0.2371708245126284*alphar[17])+0.2371708245126284*(alphar[10]+alphar[9])+0.1767766952966368*alphar[6]-0.2371708245126284*alphar[4]-0.1767766952966368*(alphar[2]+alphar[1])+0.1767766952966368*alphar[0] > 0) { 
    fUpwindQuad_r[10] = hyb_3x3v_p1_surfx1_eval_quad_node_10_r(fc); 
  } else { 
    fUpwindQuad_r[10] = hyb_3x3v_p1_surfx1_eval_quad_node_10_l(fr); 
  } 
  if (0.2371708245126284*alphal[20]-0.2371708245126284*(alphal[17]+alphal[13]+alphal[12])+0.2371708245126284*(alphal[10]+alphal[9])+0.1767766952966368*alphal[6]+0.2371708245126284*alphal[5]-0.2371708245126284*alphal[4]-0.1767766952966368*(alphal[2]+alphal[1])+0.1767766952966368*alphal[0] > 0) { 
    fUpwindQuad_l[11] = hyb_3x3v_p1_surfx1_eval_quad_node_11_r(fl); 
  } else { 
    fUpwindQuad_l[11] = hyb_3x3v_p1_surfx1_eval_quad_node_11_l(fc); 
  } 
  if (0.2371708245126284*alphar[20]-0.2371708245126284*(alphar[17]+alphar[13]+alphar[12])+0.2371708245126284*(alphar[10]+alphar[9])+0.1767766952966368*alphar[6]+0.2371708245126284*alphar[5]-0.2371708245126284*alphar[4]-0.1767766952966368*(alphar[2]+alphar[1])+0.1767766952966368*alphar[0] > 0) { 
    fUpwindQuad_r[11] = hyb_3x3v_p1_surfx1_eval_quad_node_11_r(fc); 
  } else { 
    fUpwindQuad_r[11] = hyb_3x3v_p1_surfx1_eval_quad_node_11_l(fr); 
  } 
  if ((-0.2371708245126284*alphal[20])+0.2371708245126284*(alphal[13]+alphal[12])+0.1767766952966368*alphal[6]-0.2371708245126284*alphal[5]-0.1767766952966368*(alphal[2]+alphal[1])+0.1767766952966368*alphal[0] > 0) { 
    fUpwindQuad_l[12] = hyb_3x3v_p1_surfx1_eval_quad_node_12_r(fl); 
  } else { 
    fUpwindQuad_l[12] = hyb_3x3v_p1_surfx1_eval_quad_node_12_l(fc); 
  } 
  if ((-0.2371708245126284*alphar[20])+0.2371708245126284*(alphar[13]+alphar[12])+0.1767766952966368*alphar[6]-0.2371708245126284*alphar[5]-0.1767766952966368*(alphar[2]+alphar[1])+0.1767766952966368*alphar[0] > 0) { 
    fUpwindQuad_r[12] = hyb_3x3v_p1_surfx1_eval_quad_node_12_r(fc); 
  } else { 
    fUpwindQuad_r[12] = hyb_3x3v_p1_surfx1_eval_quad_node_12_l(fr); 
  } 
  if (0.1767766952966368*alphal[6]-0.1767766952966368*(alphal[2]+alphal[1])+0.1767766952966368*alphal[0] > 0) { 
    fUpwindQuad_l[13] = hyb_3x3v_p1_surfx1_eval_quad_node_13_r(fl); 
  } else { 
    fUpwindQuad_l[13] = hyb_3x3v_p1_surfx1_eval_quad_node_13_l(fc); 
  } 
  if (0.1767766952966368*alphar[6]-0.1767766952966368*(alphar[2]+alphar[1])+0.1767766952966368*alphar[0] > 0) { 
    fUpwindQuad_r[13] = hyb_3x3v_p1_surfx1_eval_quad_node_13_r(fc); 
  } else { 
    fUpwindQuad_r[13] = hyb_3x3v_p1_surfx1_eval_quad_node_13_l(fr); 
  } 
  if (0.2371708245126284*alphal[20]-0.2371708245126284*(alphal[13]+alphal[12])+0.1767766952966368*alphal[6]+0.2371708245126284*alphal[5]-0.1767766952966368*(alphal[2]+alphal[1])+0.1767766952966368*alphal[0] > 0) { 
    fUpwindQuad_l[14] = hyb_3x3v_p1_surfx1_eval_quad_node_14_r(fl); 
  } else { 
    fUpwindQuad_l[14] = hyb_3x3v_p1_surfx1_eval_quad_node_14_l(fc); 
  } 
  if (0.2371708245126284*alphar[20]-0.2371708245126284*(alphar[13]+alphar[12])+0.1767766952966368*alphar[6]+0.2371708245126284*alphar[5]-0.1767766952966368*(alphar[2]+alphar[1])+0.1767766952966368*alphar[0] > 0) { 
    fUpwindQuad_r[14] = hyb_3x3v_p1_surfx1_eval_quad_node_14_r(fc); 
  } else { 
    fUpwindQuad_r[14] = hyb_3x3v_p1_surfx1_eval_quad_node_14_l(fr); 
  } 
  if ((-0.2371708245126284*alphal[20])+0.2371708245126284*(alphal[17]+alphal[13]+alphal[12])-0.2371708245126284*(alphal[10]+alphal[9])+0.1767766952966368*alphal[6]-0.2371708245126284*alphal[5]+0.2371708245126284*alphal[4]-0.1767766952966368*(alphal[2]+alphal[1])+0.1767766952966368*alphal[0] > 0) { 
    fUpwindQuad_l[15] = hyb_3x3v_p1_surfx1_eval_quad_node_15_r(fl); 
  } else { 
    fUpwindQuad_l[15] = hyb_3x3v_p1_surfx1_eval_quad_node_15_l(fc); 
  } 
  if ((-0.2371708245126284*alphar[20])+0.2371708245126284*(alphar[17]+alphar[13]+alphar[12])-0.2371708245126284*(alphar[10]+alphar[9])+0.1767766952966368*alphar[6]-0.2371708245126284*alphar[5]+0.2371708245126284*alphar[4]-0.1767766952966368*(alphar[2]+alphar[1])+0.1767766952966368*alphar[0] > 0) { 
    fUpwindQuad_r[15] = hyb_3x3v_p1_surfx1_eval_quad_node_15_r(fc); 
  } else { 
    fUpwindQuad_r[15] = hyb_3x3v_p1_surfx1_eval_quad_node_15_l(fr); 
  } 
  if (0.2371708245126284*alphal[17]-0.2371708245126284*(alphal[10]+alphal[9])+0.1767766952966368*alphal[6]+0.2371708245126284*alphal[4]-0.1767766952966368*(alphal[2]+alphal[1])+0.1767766952966368*alphal[0] > 0) { 
    fUpwindQuad_l[16] = hyb_3x3v_p1_surfx1_eval_quad_node_16_r(fl); 
  } else { 
    fUpwindQuad_l[16] = hyb_3x3v_p1_surfx1_eval_quad_node_16_l(fc); 
  } 
  if (0.2371708245126284*alphar[17]-0.2371708245126284*(alphar[10]+alphar[9])+0.1767766952966368*alphar[6]+0.2371708245126284*alphar[4]-0.1767766952966368*(alphar[2]+alphar[1])+0.1767766952966368*alphar[0] > 0) { 
    fUpwindQuad_r[16] = hyb_3x3v_p1_surfx1_eval_quad_node_16_r(fc); 
  } else { 
    fUpwindQuad_r[16] = hyb_3x3v_p1_surfx1_eval_quad_node_16_l(fr); 
  } 
  if (0.2371708245126284*(alphal[20]+alphal[17])-0.2371708245126284*(alphal[13]+alphal[12]+alphal[10]+alphal[9])+0.1767766952966368*alphal[6]+0.2371708245126284*(alphal[5]+alphal[4])-0.1767766952966368*(alphal[2]+alphal[1])+0.1767766952966368*alphal[0] > 0) { 
    fUpwindQuad_l[17] = hyb_3x3v_p1_surfx1_eval_quad_node_17_r(fl); 
  } else { 
    fUpwindQuad_l[17] = hyb_3x3v_p1_surfx1_eval_quad_node_17_l(fc); 
  } 
  if (0.2371708245126284*(alphar[20]+alphar[17])-0.2371708245126284*(alphar[13]+alphar[12]+alphar[10]+alphar[9])+0.1767766952966368*alphar[6]+0.2371708245126284*(alphar[5]+alphar[4])-0.1767766952966368*(alphar[2]+alphar[1])+0.1767766952966368*alphar[0] > 0) { 
    fUpwindQuad_r[17] = hyb_3x3v_p1_surfx1_eval_quad_node_17_r(fc); 
  } else { 
    fUpwindQuad_r[17] = hyb_3x3v_p1_surfx1_eval_quad_node_17_l(fr); 
  } 
  if ((-0.2371708245126284*(alphal[20]+alphal[17]))+0.2371708245126284*(alphal[16]+alphal[13]+alphal[12]+alphal[10]+alphal[9])-0.2371708245126284*(alphal[8]+alphal[7])+0.1767766952966368*alphal[6]-0.2371708245126284*(alphal[5]+alphal[4])+0.2371708245126284*alphal[3]-0.1767766952966368*(alphal[2]+alphal[1])+0.1767766952966368*alphal[0] > 0) { 
    fUpwindQuad_l[18] = hyb_3x3v_p1_surfx1_eval_quad_node_18_r(fl); 
  } else { 
    fUpwindQuad_l[18] = hyb_3x3v_p1_surfx1_eval_quad_node_18_l(fc); 
  } 
  if ((-0.2371708245126284*(alphar[20]+alphar[17]))+0.2371708245126284*(alphar[16]+alphar[13]+alphar[12]+alphar[10]+alphar[9])-0.2371708245126284*(alphar[8]+alphar[7])+0.1767766952966368*alphar[6]-0.2371708245126284*(alphar[5]+alphar[4])+0.2371708245126284*alphar[3]-0.1767766952966368*(alphar[2]+alphar[1])+0.1767766952966368*alphar[0] > 0) { 
    fUpwindQuad_r[18] = hyb_3x3v_p1_surfx1_eval_quad_node_18_r(fc); 
  } else { 
    fUpwindQuad_r[18] = hyb_3x3v_p1_surfx1_eval_quad_node_18_l(fr); 
  } 
  if ((-0.2371708245126284*alphal[17])+0.2371708245126284*(alphal[16]+alphal[10]+alphal[9])-0.2371708245126284*(alphal[8]+alphal[7])+0.1767766952966368*alphal[6]-0.2371708245126284*alphal[4]+0.2371708245126284*alphal[3]-0.1767766952966368*(alphal[2]+alphal[1])+0.1767766952966368*alphal[0] > 0) { 
    fUpwindQuad_l[19] = hyb_3x3v_p1_surfx1_eval_quad_node_19_r(fl); 
  } else { 
    fUpwindQuad_l[19] = hyb_3x3v_p1_surfx1_eval_quad_node_19_l(fc); 
  } 
  if ((-0.2371708245126284*alphar[17])+0.2371708245126284*(alphar[16]+alphar[10]+alphar[9])-0.2371708245126284*(alphar[8]+alphar[7])+0.1767766952966368*alphar[6]-0.2371708245126284*alphar[4]+0.2371708245126284*alphar[3]-0.1767766952966368*(alphar[2]+alphar[1])+0.1767766952966368*alphar[0] > 0) { 
    fUpwindQuad_r[19] = hyb_3x3v_p1_surfx1_eval_quad_node_19_r(fc); 
  } else { 
    fUpwindQuad_r[19] = hyb_3x3v_p1_surfx1_eval_quad_node_19_l(fr); 
  } 
  if (0.2371708245126284*alphal[20]-0.2371708245126284*alphal[17]+0.2371708245126284*alphal[16]-0.2371708245126284*(alphal[13]+alphal[12])+0.2371708245126284*(alphal[10]+alphal[9])-0.2371708245126284*(alphal[8]+alphal[7])+0.1767766952966368*alphal[6]+0.2371708245126284*alphal[5]-0.2371708245126284*alphal[4]+0.2371708245126284*alphal[3]-0.1767766952966368*(alphal[2]+alphal[1])+0.1767766952966368*alphal[0] > 0) { 
    fUpwindQuad_l[20] = hyb_3x3v_p1_surfx1_eval_quad_node_20_r(fl); 
  } else { 
    fUpwindQuad_l[20] = hyb_3x3v_p1_surfx1_eval_quad_node_20_l(fc); 
  } 
  if (0.2371708245126284*alphar[20]-0.2371708245126284*alphar[17]+0.2371708245126284*alphar[16]-0.2371708245126284*(alphar[13]+alphar[12])+0.2371708245126284*(alphar[10]+alphar[9])-0.2371708245126284*(alphar[8]+alphar[7])+0.1767766952966368*alphar[6]+0.2371708245126284*alphar[5]-0.2371708245126284*alphar[4]+0.2371708245126284*alphar[3]-0.1767766952966368*(alphar[2]+alphar[1])+0.1767766952966368*alphar[0] > 0) { 
    fUpwindQuad_r[20] = hyb_3x3v_p1_surfx1_eval_quad_node_20_r(fc); 
  } else { 
    fUpwindQuad_r[20] = hyb_3x3v_p1_surfx1_eval_quad_node_20_l(fr); 
  } 
  if ((-0.2371708245126284*alphal[20])+0.2371708245126284*(alphal[16]+alphal[13]+alphal[12])-0.2371708245126284*(alphal[8]+alphal[7])+0.1767766952966368*alphal[6]-0.2371708245126284*alphal[5]+0.2371708245126284*alphal[3]-0.1767766952966368*(alphal[2]+alphal[1])+0.1767766952966368*alphal[0] > 0) { 
    fUpwindQuad_l[21] = hyb_3x3v_p1_surfx1_eval_quad_node_21_r(fl); 
  } else { 
    fUpwindQuad_l[21] = hyb_3x3v_p1_surfx1_eval_quad_node_21_l(fc); 
  } 
  if ((-0.2371708245126284*alphar[20])+0.2371708245126284*(alphar[16]+alphar[13]+alphar[12])-0.2371708245126284*(alphar[8]+alphar[7])+0.1767766952966368*alphar[6]-0.2371708245126284*alphar[5]+0.2371708245126284*alphar[3]-0.1767766952966368*(alphar[2]+alphar[1])+0.1767766952966368*alphar[0] > 0) { 
    fUpwindQuad_r[21] = hyb_3x3v_p1_surfx1_eval_quad_node_21_r(fc); 
  } else { 
    fUpwindQuad_r[21] = hyb_3x3v_p1_surfx1_eval_quad_node_21_l(fr); 
  } 
  if (0.2371708245126284*alphal[16]-0.2371708245126284*(alphal[8]+alphal[7])+0.1767766952966368*alphal[6]+0.2371708245126284*alphal[3]-0.1767766952966368*(alphal[2]+alphal[1])+0.1767766952966368*alphal[0] > 0) { 
    fUpwindQuad_l[22] = hyb_3x3v_p1_surfx1_eval_quad_node_22_r(fl); 
  } else { 
    fUpwindQuad_l[22] = hyb_3x3v_p1_surfx1_eval_quad_node_22_l(fc); 
  } 
  if (0.2371708245126284*alphar[16]-0.2371708245126284*(alphar[8]+alphar[7])+0.1767766952966368*alphar[6]+0.2371708245126284*alphar[3]-0.1767766952966368*(alphar[2]+alphar[1])+0.1767766952966368*alphar[0] > 0) { 
    fUpwindQuad_r[22] = hyb_3x3v_p1_surfx1_eval_quad_node_22_r(fc); 
  } else { 
    fUpwindQuad_r[22] = hyb_3x3v_p1_surfx1_eval_quad_node_22_l(fr); 
  } 
  if (0.2371708245126284*(alphal[20]+alphal[16])-0.2371708245126284*(alphal[13]+alphal[12]+alphal[8]+alphal[7])+0.1767766952966368*alphal[6]+0.2371708245126284*(alphal[5]+alphal[3])-0.1767766952966368*(alphal[2]+alphal[1])+0.1767766952966368*alphal[0] > 0) { 
    fUpwindQuad_l[23] = hyb_3x3v_p1_surfx1_eval_quad_node_23_r(fl); 
  } else { 
    fUpwindQuad_l[23] = hyb_3x3v_p1_surfx1_eval_quad_node_23_l(fc); 
  } 
  if (0.2371708245126284*(alphar[20]+alphar[16])-0.2371708245126284*(alphar[13]+alphar[12]+alphar[8]+alphar[7])+0.1767766952966368*alphar[6]+0.2371708245126284*(alphar[5]+alphar[3])-0.1767766952966368*(alphar[2]+alphar[1])+0.1767766952966368*alphar[0] > 0) { 
    fUpwindQuad_r[23] = hyb_3x3v_p1_surfx1_eval_quad_node_23_r(fc); 
  } else { 
    fUpwindQuad_r[23] = hyb_3x3v_p1_surfx1_eval_quad_node_23_l(fr); 
  } 
  if ((-0.2371708245126284*alphal[20])+0.2371708245126284*(alphal[17]+alphal[16]+alphal[13]+alphal[12])-0.2371708245126284*(alphal[10]+alphal[9]+alphal[8]+alphal[7])+0.1767766952966368*alphal[6]-0.2371708245126284*alphal[5]+0.2371708245126284*(alphal[4]+alphal[3])-0.1767766952966368*(alphal[2]+alphal[1])+0.1767766952966368*alphal[0] > 0) { 
    fUpwindQuad_l[24] = hyb_3x3v_p1_surfx1_eval_quad_node_24_r(fl); 
  } else { 
    fUpwindQuad_l[24] = hyb_3x3v_p1_surfx1_eval_quad_node_24_l(fc); 
  } 
  if ((-0.2371708245126284*alphar[20])+0.2371708245126284*(alphar[17]+alphar[16]+alphar[13]+alphar[12])-0.2371708245126284*(alphar[10]+alphar[9]+alphar[8]+alphar[7])+0.1767766952966368*alphar[6]-0.2371708245126284*alphar[5]+0.2371708245126284*(alphar[4]+alphar[3])-0.1767766952966368*(alphar[2]+alphar[1])+0.1767766952966368*alphar[0] > 0) { 
    fUpwindQuad_r[24] = hyb_3x3v_p1_surfx1_eval_quad_node_24_r(fc); 
  } else { 
    fUpwindQuad_r[24] = hyb_3x3v_p1_surfx1_eval_quad_node_24_l(fr); 
  } 
  if (0.2371708245126284*(alphal[17]+alphal[16])-0.2371708245126284*(alphal[10]+alphal[9]+alphal[8]+alphal[7])+0.1767766952966368*alphal[6]+0.2371708245126284*(alphal[4]+alphal[3])-0.1767766952966368*(alphal[2]+alphal[1])+0.1767766952966368*alphal[0] > 0) { 
    fUpwindQuad_l[25] = hyb_3x3v_p1_surfx1_eval_quad_node_25_r(fl); 
  } else { 
    fUpwindQuad_l[25] = hyb_3x3v_p1_surfx1_eval_quad_node_25_l(fc); 
  } 
  if (0.2371708245126284*(alphar[17]+alphar[16])-0.2371708245126284*(alphar[10]+alphar[9]+alphar[8]+alphar[7])+0.1767766952966368*alphar[6]+0.2371708245126284*(alphar[4]+alphar[3])-0.1767766952966368*(alphar[2]+alphar[1])+0.1767766952966368*alphar[0] > 0) { 
    fUpwindQuad_r[25] = hyb_3x3v_p1_surfx1_eval_quad_node_25_r(fc); 
  } else { 
    fUpwindQuad_r[25] = hyb_3x3v_p1_surfx1_eval_quad_node_25_l(fr); 
  } 
  if (0.2371708245126284*(alphal[20]+alphal[17]+alphal[16])-0.2371708245126284*(alphal[13]+alphal[12]+alphal[10]+alphal[9]+alphal[8]+alphal[7])+0.1767766952966368*alphal[6]+0.2371708245126284*(alphal[5]+alphal[4]+alphal[3])-0.1767766952966368*(alphal[2]+alphal[1])+0.1767766952966368*alphal[0] > 0) { 
    fUpwindQuad_l[26] = hyb_3x3v_p1_surfx1_eval_quad_node_26_r(fl); 
  } else { 
    fUpwindQuad_l[26] = hyb_3x3v_p1_surfx1_eval_quad_node_26_l(fc); 
  } 
  if (0.2371708245126284*(alphar[20]+alphar[17]+alphar[16])-0.2371708245126284*(alphar[13]+alphar[12]+alphar[10]+alphar[9]+alphar[8]+alphar[7])+0.1767766952966368*alphar[6]+0.2371708245126284*(alphar[5]+alphar[4]+alphar[3])-0.1767766952966368*(alphar[2]+alphar[1])+0.1767766952966368*alphar[0] > 0) { 
    fUpwindQuad_r[26] = hyb_3x3v_p1_surfx1_eval_quad_node_26_r(fc); 
  } else { 
    fUpwindQuad_r[26] = hyb_3x3v_p1_surfx1_eval_quad_node_26_l(fr); 
  } 
  if (0.2371708245126284*(alphal[20]+alphal[17]+alphal[16])-0.2371708245126284*alphal[13]+0.2371708245126284*alphal[12]-0.2371708245126284*alphal[10]+0.2371708245126284*alphal[9]-0.2371708245126284*alphal[8]+0.2371708245126284*alphal[7]-0.1767766952966368*alphal[6]-0.2371708245126284*(alphal[5]+alphal[4]+alphal[3])+0.1767766952966368*alphal[2]-0.1767766952966368*alphal[1]+0.1767766952966368*alphal[0] > 0) { 
    fUpwindQuad_l[27] = hyb_3x3v_p1_surfx1_eval_quad_node_27_r(fl); 
  } else { 
    fUpwindQuad_l[27] = hyb_3x3v_p1_surfx1_eval_quad_node_27_l(fc); 
  } 
  if (0.2371708245126284*(alphar[20]+alphar[17]+alphar[16])-0.2371708245126284*alphar[13]+0.2371708245126284*alphar[12]-0.2371708245126284*alphar[10]+0.2371708245126284*alphar[9]-0.2371708245126284*alphar[8]+0.2371708245126284*alphar[7]-0.1767766952966368*alphar[6]-0.2371708245126284*(alphar[5]+alphar[4]+alphar[3])+0.1767766952966368*alphar[2]-0.1767766952966368*alphar[1]+0.1767766952966368*alphar[0] > 0) { 
    fUpwindQuad_r[27] = hyb_3x3v_p1_surfx1_eval_quad_node_27_r(fc); 
  } else { 
    fUpwindQuad_r[27] = hyb_3x3v_p1_surfx1_eval_quad_node_27_l(fr); 
  } 
  if (0.2371708245126284*(alphal[17]+alphal[16])-0.2371708245126284*alphal[10]+0.2371708245126284*alphal[9]-0.2371708245126284*alphal[8]+0.2371708245126284*alphal[7]-0.1767766952966368*alphal[6]-0.2371708245126284*(alphal[4]+alphal[3])+0.1767766952966368*alphal[2]-0.1767766952966368*alphal[1]+0.1767766952966368*alphal[0] > 0) { 
    fUpwindQuad_l[28] = hyb_3x3v_p1_surfx1_eval_quad_node_28_r(fl); 
  } else { 
    fUpwindQuad_l[28] = hyb_3x3v_p1_surfx1_eval_quad_node_28_l(fc); 
  } 
  if (0.2371708245126284*(alphar[17]+alphar[16])-0.2371708245126284*alphar[10]+0.2371708245126284*alphar[9]-0.2371708245126284*alphar[8]+0.2371708245126284*alphar[7]-0.1767766952966368*alphar[6]-0.2371708245126284*(alphar[4]+alphar[3])+0.1767766952966368*alphar[2]-0.1767766952966368*alphar[1]+0.1767766952966368*alphar[0] > 0) { 
    fUpwindQuad_r[28] = hyb_3x3v_p1_surfx1_eval_quad_node_28_r(fc); 
  } else { 
    fUpwindQuad_r[28] = hyb_3x3v_p1_surfx1_eval_quad_node_28_l(fr); 
  } 
  if ((-0.2371708245126284*alphal[20])+0.2371708245126284*(alphal[17]+alphal[16]+alphal[13])-0.2371708245126284*(alphal[12]+alphal[10])+0.2371708245126284*alphal[9]-0.2371708245126284*alphal[8]+0.2371708245126284*alphal[7]-0.1767766952966368*alphal[6]+0.2371708245126284*alphal[5]-0.2371708245126284*(alphal[4]+alphal[3])+0.1767766952966368*alphal[2]-0.1767766952966368*alphal[1]+0.1767766952966368*alphal[0] > 0) { 
    fUpwindQuad_l[29] = hyb_3x3v_p1_surfx1_eval_quad_node_29_r(fl); 
  } else { 
    fUpwindQuad_l[29] = hyb_3x3v_p1_surfx1_eval_quad_node_29_l(fc); 
  } 
  if ((-0.2371708245126284*alphar[20])+0.2371708245126284*(alphar[17]+alphar[16]+alphar[13])-0.2371708245126284*(alphar[12]+alphar[10])+0.2371708245126284*alphar[9]-0.2371708245126284*alphar[8]+0.2371708245126284*alphar[7]-0.1767766952966368*alphar[6]+0.2371708245126284*alphar[5]-0.2371708245126284*(alphar[4]+alphar[3])+0.1767766952966368*alphar[2]-0.1767766952966368*alphar[1]+0.1767766952966368*alphar[0] > 0) { 
    fUpwindQuad_r[29] = hyb_3x3v_p1_surfx1_eval_quad_node_29_r(fc); 
  } else { 
    fUpwindQuad_r[29] = hyb_3x3v_p1_surfx1_eval_quad_node_29_l(fr); 
  } 
  if (0.2371708245126284*(alphal[20]+alphal[16])-0.2371708245126284*alphal[13]+0.2371708245126284*alphal[12]-0.2371708245126284*alphal[8]+0.2371708245126284*alphal[7]-0.1767766952966368*alphal[6]-0.2371708245126284*(alphal[5]+alphal[3])+0.1767766952966368*alphal[2]-0.1767766952966368*alphal[1]+0.1767766952966368*alphal[0] > 0) { 
    fUpwindQuad_l[30] = hyb_3x3v_p1_surfx1_eval_quad_node_30_r(fl); 
  } else { 
    fUpwindQuad_l[30] = hyb_3x3v_p1_surfx1_eval_quad_node_30_l(fc); 
  } 
  if (0.2371708245126284*(alphar[20]+alphar[16])-0.2371708245126284*alphar[13]+0.2371708245126284*alphar[12]-0.2371708245126284*alphar[8]+0.2371708245126284*alphar[7]-0.1767766952966368*alphar[6]-0.2371708245126284*(alphar[5]+alphar[3])+0.1767766952966368*alphar[2]-0.1767766952966368*alphar[1]+0.1767766952966368*alphar[0] > 0) { 
    fUpwindQuad_r[30] = hyb_3x3v_p1_surfx1_eval_quad_node_30_r(fc); 
  } else { 
    fUpwindQuad_r[30] = hyb_3x3v_p1_surfx1_eval_quad_node_30_l(fr); 
  } 
  if (0.2371708245126284*alphal[16]-0.2371708245126284*alphal[8]+0.2371708245126284*alphal[7]-0.1767766952966368*alphal[6]-0.2371708245126284*alphal[3]+0.1767766952966368*alphal[2]-0.1767766952966368*alphal[1]+0.1767766952966368*alphal[0] > 0) { 
    fUpwindQuad_l[31] = hyb_3x3v_p1_surfx1_eval_quad_node_31_r(fl); 
  } else { 
    fUpwindQuad_l[31] = hyb_3x3v_p1_surfx1_eval_quad_node_31_l(fc); 
  } 
  if (0.2371708245126284*alphar[16]-0.2371708245126284*alphar[8]+0.2371708245126284*alphar[7]-0.1767766952966368*alphar[6]-0.2371708245126284*alphar[3]+0.1767766952966368*alphar[2]-0.1767766952966368*alphar[1]+0.1767766952966368*alphar[0] > 0) { 
    fUpwindQuad_r[31] = hyb_3x3v_p1_surfx1_eval_quad_node_31_r(fc); 
  } else { 
    fUpwindQuad_r[31] = hyb_3x3v_p1_surfx1_eval_quad_node_31_l(fr); 
  } 
  if ((-0.2371708245126284*alphal[20])+0.2371708245126284*(alphal[16]+alphal[13])-0.2371708245126284*(alphal[12]+alphal[8])+0.2371708245126284*alphal[7]-0.1767766952966368*alphal[6]+0.2371708245126284*alphal[5]-0.2371708245126284*alphal[3]+0.1767766952966368*alphal[2]-0.1767766952966368*alphal[1]+0.1767766952966368*alphal[0] > 0) { 
    fUpwindQuad_l[32] = hyb_3x3v_p1_surfx1_eval_quad_node_32_r(fl); 
  } else { 
    fUpwindQuad_l[32] = hyb_3x3v_p1_surfx1_eval_quad_node_32_l(fc); 
  } 
  if ((-0.2371708245126284*alphar[20])+0.2371708245126284*(alphar[16]+alphar[13])-0.2371708245126284*(alphar[12]+alphar[8])+0.2371708245126284*alphar[7]-0.1767766952966368*alphar[6]+0.2371708245126284*alphar[5]-0.2371708245126284*alphar[3]+0.1767766952966368*alphar[2]-0.1767766952966368*alphar[1]+0.1767766952966368*alphar[0] > 0) { 
    fUpwindQuad_r[32] = hyb_3x3v_p1_surfx1_eval_quad_node_32_r(fc); 
  } else { 
    fUpwindQuad_r[32] = hyb_3x3v_p1_surfx1_eval_quad_node_32_l(fr); 
  } 
  if (0.2371708245126284*alphal[20]-0.2371708245126284*alphal[17]+0.2371708245126284*alphal[16]-0.2371708245126284*alphal[13]+0.2371708245126284*(alphal[12]+alphal[10])-0.2371708245126284*(alphal[9]+alphal[8])+0.2371708245126284*alphal[7]-0.1767766952966368*alphal[6]-0.2371708245126284*alphal[5]+0.2371708245126284*alphal[4]-0.2371708245126284*alphal[3]+0.1767766952966368*alphal[2]-0.1767766952966368*alphal[1]+0.1767766952966368*alphal[0] > 0) { 
    fUpwindQuad_l[33] = hyb_3x3v_p1_surfx1_eval_quad_node_33_r(fl); 
  } else { 
    fUpwindQuad_l[33] = hyb_3x3v_p1_surfx1_eval_quad_node_33_l(fc); 
  } 
  if (0.2371708245126284*alphar[20]-0.2371708245126284*alphar[17]+0.2371708245126284*alphar[16]-0.2371708245126284*alphar[13]+0.2371708245126284*(alphar[12]+alphar[10])-0.2371708245126284*(alphar[9]+alphar[8])+0.2371708245126284*alphar[7]-0.1767766952966368*alphar[6]-0.2371708245126284*alphar[5]+0.2371708245126284*alphar[4]-0.2371708245126284*alphar[3]+0.1767766952966368*alphar[2]-0.1767766952966368*alphar[1]+0.1767766952966368*alphar[0] > 0) { 
    fUpwindQuad_r[33] = hyb_3x3v_p1_surfx1_eval_quad_node_33_r(fc); 
  } else { 
    fUpwindQuad_r[33] = hyb_3x3v_p1_surfx1_eval_quad_node_33_l(fr); 
  } 
  if ((-0.2371708245126284*alphal[17])+0.2371708245126284*(alphal[16]+alphal[10])-0.2371708245126284*(alphal[9]+alphal[8])+0.2371708245126284*alphal[7]-0.1767766952966368*alphal[6]+0.2371708245126284*alphal[4]-0.2371708245126284*alphal[3]+0.1767766952966368*alphal[2]-0.1767766952966368*alphal[1]+0.1767766952966368*alphal[0] > 0) { 
    fUpwindQuad_l[34] = hyb_3x3v_p1_surfx1_eval_quad_node_34_r(fl); 
  } else { 
    fUpwindQuad_l[34] = hyb_3x3v_p1_surfx1_eval_quad_node_34_l(fc); 
  } 
  if ((-0.2371708245126284*alphar[17])+0.2371708245126284*(alphar[16]+alphar[10])-0.2371708245126284*(alphar[9]+alphar[8])+0.2371708245126284*alphar[7]-0.1767766952966368*alphar[6]+0.2371708245126284*alphar[4]-0.2371708245126284*alphar[3]+0.1767766952966368*alphar[2]-0.1767766952966368*alphar[1]+0.1767766952966368*alphar[0] > 0) { 
    fUpwindQuad_r[34] = hyb_3x3v_p1_surfx1_eval_quad_node_34_r(fc); 
  } else { 
    fUpwindQuad_r[34] = hyb_3x3v_p1_surfx1_eval_quad_node_34_l(fr); 
  } 
  if ((-0.2371708245126284*(alphal[20]+alphal[17]))+0.2371708245126284*(alphal[16]+alphal[13])-0.2371708245126284*alphal[12]+0.2371708245126284*alphal[10]-0.2371708245126284*(alphal[9]+alphal[8])+0.2371708245126284*alphal[7]-0.1767766952966368*alphal[6]+0.2371708245126284*(alphal[5]+alphal[4])-0.2371708245126284*alphal[3]+0.1767766952966368*alphal[2]-0.1767766952966368*alphal[1]+0.1767766952966368*alphal[0] > 0) { 
    fUpwindQuad_l[35] = hyb_3x3v_p1_surfx1_eval_quad_node_35_r(fl); 
  } else { 
    fUpwindQuad_l[35] = hyb_3x3v_p1_surfx1_eval_quad_node_35_l(fc); 
  } 
  if ((-0.2371708245126284*(alphar[20]+alphar[17]))+0.2371708245126284*(alphar[16]+alphar[13])-0.2371708245126284*alphar[12]+0.2371708245126284*alphar[10]-0.2371708245126284*(alphar[9]+alphar[8])+0.2371708245126284*alphar[7]-0.1767766952966368*alphar[6]+0.2371708245126284*(alphar[5]+alphar[4])-0.2371708245126284*alphar[3]+0.1767766952966368*alphar[2]-0.1767766952966368*alphar[1]+0.1767766952966368*alphar[0] > 0) { 
    fUpwindQuad_r[35] = hyb_3x3v_p1_surfx1_eval_quad_node_35_r(fc); 
  } else { 
    fUpwindQuad_r[35] = hyb_3x3v_p1_surfx1_eval_quad_node_35_l(fr); 
  } 
  if (0.2371708245126284*(alphal[20]+alphal[17])-0.2371708245126284*alphal[13]+0.2371708245126284*alphal[12]-0.2371708245126284*alphal[10]+0.2371708245126284*alphal[9]-0.1767766952966368*alphal[6]-0.2371708245126284*(alphal[5]+alphal[4])+0.1767766952966368*alphal[2]-0.1767766952966368*alphal[1]+0.1767766952966368*alphal[0] > 0) { 
    fUpwindQuad_l[36] = hyb_3x3v_p1_surfx1_eval_quad_node_36_r(fl); 
  } else { 
    fUpwindQuad_l[36] = hyb_3x3v_p1_surfx1_eval_quad_node_36_l(fc); 
  } 
  if (0.2371708245126284*(alphar[20]+alphar[17])-0.2371708245126284*alphar[13]+0.2371708245126284*alphar[12]-0.2371708245126284*alphar[10]+0.2371708245126284*alphar[9]-0.1767766952966368*alphar[6]-0.2371708245126284*(alphar[5]+alphar[4])+0.1767766952966368*alphar[2]-0.1767766952966368*alphar[1]+0.1767766952966368*alphar[0] > 0) { 
    fUpwindQuad_r[36] = hyb_3x3v_p1_surfx1_eval_quad_node_36_r(fc); 
  } else { 
    fUpwindQuad_r[36] = hyb_3x3v_p1_surfx1_eval_quad_node_36_l(fr); 
  } 
  if (0.2371708245126284*alphal[17]-0.2371708245126284*alphal[10]+0.2371708245126284*alphal[9]-0.1767766952966368*alphal[6]-0.2371708245126284*alphal[4]+0.1767766952966368*alphal[2]-0.1767766952966368*alphal[1]+0.1767766952966368*alphal[0] > 0) { 
    fUpwindQuad_l[37] = hyb_3x3v_p1_surfx1_eval_quad_node_37_r(fl); 
  } else { 
    fUpwindQuad_l[37] = hyb_3x3v_p1_surfx1_eval_quad_node_37_l(fc); 
  } 
  if (0.2371708245126284*alphar[17]-0.2371708245126284*alphar[10]+0.2371708245126284*alphar[9]-0.1767766952966368*alphar[6]-0.2371708245126284*alphar[4]+0.1767766952966368*alphar[2]-0.1767766952966368*alphar[1]+0.1767766952966368*alphar[0] > 0) { 
    fUpwindQuad_r[37] = hyb_3x3v_p1_surfx1_eval_quad_node_37_r(fc); 
  } else { 
    fUpwindQuad_r[37] = hyb_3x3v_p1_surfx1_eval_quad_node_37_l(fr); 
  } 
  if ((-0.2371708245126284*alphal[20])+0.2371708245126284*(alphal[17]+alphal[13])-0.2371708245126284*(alphal[12]+alphal[10])+0.2371708245126284*alphal[9]-0.1767766952966368*alphal[6]+0.2371708245126284*alphal[5]-0.2371708245126284*alphal[4]+0.1767766952966368*alphal[2]-0.1767766952966368*alphal[1]+0.1767766952966368*alphal[0] > 0) { 
    fUpwindQuad_l[38] = hyb_3x3v_p1_surfx1_eval_quad_node_38_r(fl); 
  } else { 
    fUpwindQuad_l[38] = hyb_3x3v_p1_surfx1_eval_quad_node_38_l(fc); 
  } 
  if ((-0.2371708245126284*alphar[20])+0.2371708245126284*(alphar[17]+alphar[13])-0.2371708245126284*(alphar[12]+alphar[10])+0.2371708245126284*alphar[9]-0.1767766952966368*alphar[6]+0.2371708245126284*alphar[5]-0.2371708245126284*alphar[4]+0.1767766952966368*alphar[2]-0.1767766952966368*alphar[1]+0.1767766952966368*alphar[0] > 0) { 
    fUpwindQuad_r[38] = hyb_3x3v_p1_surfx1_eval_quad_node_38_r(fc); 
  } else { 
    fUpwindQuad_r[38] = hyb_3x3v_p1_surfx1_eval_quad_node_38_l(fr); 
  } 
  if (0.2371708245126284*alphal[20]-0.2371708245126284*alphal[13]+0.2371708245126284*alphal[12]-0.1767766952966368*alphal[6]-0.2371708245126284*alphal[5]+0.1767766952966368*alphal[2]-0.1767766952966368*alphal[1]+0.1767766952966368*alphal[0] > 0) { 
    fUpwindQuad_l[39] = hyb_3x3v_p1_surfx1_eval_quad_node_39_r(fl); 
  } else { 
    fUpwindQuad_l[39] = hyb_3x3v_p1_surfx1_eval_quad_node_39_l(fc); 
  } 
  if (0.2371708245126284*alphar[20]-0.2371708245126284*alphar[13]+0.2371708245126284*alphar[12]-0.1767766952966368*alphar[6]-0.2371708245126284*alphar[5]+0.1767766952966368*alphar[2]-0.1767766952966368*alphar[1]+0.1767766952966368*alphar[0] > 0) { 
    fUpwindQuad_r[39] = hyb_3x3v_p1_surfx1_eval_quad_node_39_r(fc); 
  } else { 
    fUpwindQuad_r[39] = hyb_3x3v_p1_surfx1_eval_quad_node_39_l(fr); 
  } 
  if ((-0.1767766952966368*alphal[6])+0.1767766952966368*alphal[2]-0.1767766952966368*alphal[1]+0.1767766952966368*alphal[0] > 0) { 
    fUpwindQuad_l[40] = hyb_3x3v_p1_surfx1_eval_quad_node_40_r(fl); 
  } else { 
    fUpwindQuad_l[40] = hyb_3x3v_p1_surfx1_eval_quad_node_40_l(fc); 
  } 
  if ((-0.1767766952966368*alphar[6])+0.1767766952966368*alphar[2]-0.1767766952966368*alphar[1]+0.1767766952966368*alphar[0] > 0) { 
    fUpwindQuad_r[40] = hyb_3x3v_p1_surfx1_eval_quad_node_40_r(fc); 
  } else { 
    fUpwindQuad_r[40] = hyb_3x3v_p1_surfx1_eval_quad_node_40_l(fr); 
  } 
  if ((-0.2371708245126284*alphal[20])+0.2371708245126284*alphal[13]-0.2371708245126284*alphal[12]-0.1767766952966368*alphal[6]+0.2371708245126284*alphal[5]+0.1767766952966368*alphal[2]-0.1767766952966368*alphal[1]+0.1767766952966368*alphal[0] > 0) { 
    fUpwindQuad_l[41] = hyb_3x3v_p1_surfx1_eval_quad_node_41_r(fl); 
  } else { 
    fUpwindQuad_l[41] = hyb_3x3v_p1_surfx1_eval_quad_node_41_l(fc); 
  } 
  if ((-0.2371708245126284*alphar[20])+0.2371708245126284*alphar[13]-0.2371708245126284*alphar[12]-0.1767766952966368*alphar[6]+0.2371708245126284*alphar[5]+0.1767766952966368*alphar[2]-0.1767766952966368*alphar[1]+0.1767766952966368*alphar[0] > 0) { 
    fUpwindQuad_r[41] = hyb_3x3v_p1_surfx1_eval_quad_node_41_r(fc); 
  } else { 
    fUpwindQuad_r[41] = hyb_3x3v_p1_surfx1_eval_quad_node_41_l(fr); 
  } 
  if (0.2371708245126284*alphal[20]-0.2371708245126284*(alphal[17]+alphal[13])+0.2371708245126284*(alphal[12]+alphal[10])-0.2371708245126284*alphal[9]-0.1767766952966368*alphal[6]-0.2371708245126284*alphal[5]+0.2371708245126284*alphal[4]+0.1767766952966368*alphal[2]-0.1767766952966368*alphal[1]+0.1767766952966368*alphal[0] > 0) { 
    fUpwindQuad_l[42] = hyb_3x3v_p1_surfx1_eval_quad_node_42_r(fl); 
  } else { 
    fUpwindQuad_l[42] = hyb_3x3v_p1_surfx1_eval_quad_node_42_l(fc); 
  } 
  if (0.2371708245126284*alphar[20]-0.2371708245126284*(alphar[17]+alphar[13])+0.2371708245126284*(alphar[12]+alphar[10])-0.2371708245126284*alphar[9]-0.1767766952966368*alphar[6]-0.2371708245126284*alphar[5]+0.2371708245126284*alphar[4]+0.1767766952966368*alphar[2]-0.1767766952966368*alphar[1]+0.1767766952966368*alphar[0] > 0) { 
    fUpwindQuad_r[42] = hyb_3x3v_p1_surfx1_eval_quad_node_42_r(fc); 
  } else { 
    fUpwindQuad_r[42] = hyb_3x3v_p1_surfx1_eval_quad_node_42_l(fr); 
  } 
  if ((-0.2371708245126284*alphal[17])+0.2371708245126284*alphal[10]-0.2371708245126284*alphal[9]-0.1767766952966368*alphal[6]+0.2371708245126284*alphal[4]+0.1767766952966368*alphal[2]-0.1767766952966368*alphal[1]+0.1767766952966368*alphal[0] > 0) { 
    fUpwindQuad_l[43] = hyb_3x3v_p1_surfx1_eval_quad_node_43_r(fl); 
  } else { 
    fUpwindQuad_l[43] = hyb_3x3v_p1_surfx1_eval_quad_node_43_l(fc); 
  } 
  if ((-0.2371708245126284*alphar[17])+0.2371708245126284*alphar[10]-0.2371708245126284*alphar[9]-0.1767766952966368*alphar[6]+0.2371708245126284*alphar[4]+0.1767766952966368*alphar[2]-0.1767766952966368*alphar[1]+0.1767766952966368*alphar[0] > 0) { 
    fUpwindQuad_r[43] = hyb_3x3v_p1_surfx1_eval_quad_node_43_r(fc); 
  } else { 
    fUpwindQuad_r[43] = hyb_3x3v_p1_surfx1_eval_quad_node_43_l(fr); 
  } 
  if ((-0.2371708245126284*(alphal[20]+alphal[17]))+0.2371708245126284*alphal[13]-0.2371708245126284*alphal[12]+0.2371708245126284*alphal[10]-0.2371708245126284*alphal[9]-0.1767766952966368*alphal[6]+0.2371708245126284*(alphal[5]+alphal[4])+0.1767766952966368*alphal[2]-0.1767766952966368*alphal[1]+0.1767766952966368*alphal[0] > 0) { 
    fUpwindQuad_l[44] = hyb_3x3v_p1_surfx1_eval_quad_node_44_r(fl); 
  } else { 
    fUpwindQuad_l[44] = hyb_3x3v_p1_surfx1_eval_quad_node_44_l(fc); 
  } 
  if ((-0.2371708245126284*(alphar[20]+alphar[17]))+0.2371708245126284*alphar[13]-0.2371708245126284*alphar[12]+0.2371708245126284*alphar[10]-0.2371708245126284*alphar[9]-0.1767766952966368*alphar[6]+0.2371708245126284*(alphar[5]+alphar[4])+0.1767766952966368*alphar[2]-0.1767766952966368*alphar[1]+0.1767766952966368*alphar[0] > 0) { 
    fUpwindQuad_r[44] = hyb_3x3v_p1_surfx1_eval_quad_node_44_r(fc); 
  } else { 
    fUpwindQuad_r[44] = hyb_3x3v_p1_surfx1_eval_quad_node_44_l(fr); 
  } 
  if (0.2371708245126284*(alphal[20]+alphal[17])-0.2371708245126284*(alphal[16]+alphal[13])+0.2371708245126284*alphal[12]-0.2371708245126284*alphal[10]+0.2371708245126284*(alphal[9]+alphal[8])-0.2371708245126284*alphal[7]-0.1767766952966368*alphal[6]-0.2371708245126284*(alphal[5]+alphal[4])+0.2371708245126284*alphal[3]+0.1767766952966368*alphal[2]-0.1767766952966368*alphal[1]+0.1767766952966368*alphal[0] > 0) { 
    fUpwindQuad_l[45] = hyb_3x3v_p1_surfx1_eval_quad_node_45_r(fl); 
  } else { 
    fUpwindQuad_l[45] = hyb_3x3v_p1_surfx1_eval_quad_node_45_l(fc); 
  } 
  if (0.2371708245126284*(alphar[20]+alphar[17])-0.2371708245126284*(alphar[16]+alphar[13])+0.2371708245126284*alphar[12]-0.2371708245126284*alphar[10]+0.2371708245126284*(alphar[9]+alphar[8])-0.2371708245126284*alphar[7]-0.1767766952966368*alphar[6]-0.2371708245126284*(alphar[5]+alphar[4])+0.2371708245126284*alphar[3]+0.1767766952966368*alphar[2]-0.1767766952966368*alphar[1]+0.1767766952966368*alphar[0] > 0) { 
    fUpwindQuad_r[45] = hyb_3x3v_p1_surfx1_eval_quad_node_45_r(fc); 
  } else { 
    fUpwindQuad_r[45] = hyb_3x3v_p1_surfx1_eval_quad_node_45_l(fr); 
  } 
  if (0.2371708245126284*alphal[17]-0.2371708245126284*(alphal[16]+alphal[10])+0.2371708245126284*(alphal[9]+alphal[8])-0.2371708245126284*alphal[7]-0.1767766952966368*alphal[6]-0.2371708245126284*alphal[4]+0.2371708245126284*alphal[3]+0.1767766952966368*alphal[2]-0.1767766952966368*alphal[1]+0.1767766952966368*alphal[0] > 0) { 
    fUpwindQuad_l[46] = hyb_3x3v_p1_surfx1_eval_quad_node_46_r(fl); 
  } else { 
    fUpwindQuad_l[46] = hyb_3x3v_p1_surfx1_eval_quad_node_46_l(fc); 
  } 
  if (0.2371708245126284*alphar[17]-0.2371708245126284*(alphar[16]+alphar[10])+0.2371708245126284*(alphar[9]+alphar[8])-0.2371708245126284*alphar[7]-0.1767766952966368*alphar[6]-0.2371708245126284*alphar[4]+0.2371708245126284*alphar[3]+0.1767766952966368*alphar[2]-0.1767766952966368*alphar[1]+0.1767766952966368*alphar[0] > 0) { 
    fUpwindQuad_r[46] = hyb_3x3v_p1_surfx1_eval_quad_node_46_r(fc); 
  } else { 
    fUpwindQuad_r[46] = hyb_3x3v_p1_surfx1_eval_quad_node_46_l(fr); 
  } 
  if ((-0.2371708245126284*alphal[20])+0.2371708245126284*alphal[17]-0.2371708245126284*alphal[16]+0.2371708245126284*alphal[13]-0.2371708245126284*(alphal[12]+alphal[10])+0.2371708245126284*(alphal[9]+alphal[8])-0.2371708245126284*alphal[7]-0.1767766952966368*alphal[6]+0.2371708245126284*alphal[5]-0.2371708245126284*alphal[4]+0.2371708245126284*alphal[3]+0.1767766952966368*alphal[2]-0.1767766952966368*alphal[1]+0.1767766952966368*alphal[0] > 0) { 
    fUpwindQuad_l[47] = hyb_3x3v_p1_surfx1_eval_quad_node_47_r(fl); 
  } else { 
    fUpwindQuad_l[47] = hyb_3x3v_p1_surfx1_eval_quad_node_47_l(fc); 
  } 
  if ((-0.2371708245126284*alphar[20])+0.2371708245126284*alphar[17]-0.2371708245126284*alphar[16]+0.2371708245126284*alphar[13]-0.2371708245126284*(alphar[12]+alphar[10])+0.2371708245126284*(alphar[9]+alphar[8])-0.2371708245126284*alphar[7]-0.1767766952966368*alphar[6]+0.2371708245126284*alphar[5]-0.2371708245126284*alphar[4]+0.2371708245126284*alphar[3]+0.1767766952966368*alphar[2]-0.1767766952966368*alphar[1]+0.1767766952966368*alphar[0] > 0) { 
    fUpwindQuad_r[47] = hyb_3x3v_p1_surfx1_eval_quad_node_47_r(fc); 
  } else { 
    fUpwindQuad_r[47] = hyb_3x3v_p1_surfx1_eval_quad_node_47_l(fr); 
  } 
  if (0.2371708245126284*alphal[20]-0.2371708245126284*(alphal[16]+alphal[13])+0.2371708245126284*(alphal[12]+alphal[8])-0.2371708245126284*alphal[7]-0.1767766952966368*alphal[6]-0.2371708245126284*alphal[5]+0.2371708245126284*alphal[3]+0.1767766952966368*alphal[2]-0.1767766952966368*alphal[1]+0.1767766952966368*alphal[0] > 0) { 
    fUpwindQuad_l[48] = hyb_3x3v_p1_surfx1_eval_quad_node_48_r(fl); 
  } else { 
    fUpwindQuad_l[48] = hyb_3x3v_p1_surfx1_eval_quad_node_48_l(fc); 
  } 
  if (0.2371708245126284*alphar[20]-0.2371708245126284*(alphar[16]+alphar[13])+0.2371708245126284*(alphar[12]+alphar[8])-0.2371708245126284*alphar[7]-0.1767766952966368*alphar[6]-0.2371708245126284*alphar[5]+0.2371708245126284*alphar[3]+0.1767766952966368*alphar[2]-0.1767766952966368*alphar[1]+0.1767766952966368*alphar[0] > 0) { 
    fUpwindQuad_r[48] = hyb_3x3v_p1_surfx1_eval_quad_node_48_r(fc); 
  } else { 
    fUpwindQuad_r[48] = hyb_3x3v_p1_surfx1_eval_quad_node_48_l(fr); 
  } 
  if ((-0.2371708245126284*alphal[16])+0.2371708245126284*alphal[8]-0.2371708245126284*alphal[7]-0.1767766952966368*alphal[6]+0.2371708245126284*alphal[3]+0.1767766952966368*alphal[2]-0.1767766952966368*alphal[1]+0.1767766952966368*alphal[0] > 0) { 
    fUpwindQuad_l[49] = hyb_3x3v_p1_surfx1_eval_quad_node_49_r(fl); 
  } else { 
    fUpwindQuad_l[49] = hyb_3x3v_p1_surfx1_eval_quad_node_49_l(fc); 
  } 
  if ((-0.2371708245126284*alphar[16])+0.2371708245126284*alphar[8]-0.2371708245126284*alphar[7]-0.1767766952966368*alphar[6]+0.2371708245126284*alphar[3]+0.1767766952966368*alphar[2]-0.1767766952966368*alphar[1]+0.1767766952966368*alphar[0] > 0) { 
    fUpwindQuad_r[49] = hyb_3x3v_p1_surfx1_eval_quad_node_49_r(fc); 
  } else { 
    fUpwindQuad_r[49] = hyb_3x3v_p1_surfx1_eval_quad_node_49_l(fr); 
  } 
  if ((-0.2371708245126284*(alphal[20]+alphal[16]))+0.2371708245126284*alphal[13]-0.2371708245126284*alphal[12]+0.2371708245126284*alphal[8]-0.2371708245126284*alphal[7]-0.1767766952966368*alphal[6]+0.2371708245126284*(alphal[5]+alphal[3])+0.1767766952966368*alphal[2]-0.1767766952966368*alphal[1]+0.1767766952966368*alphal[0] > 0) { 
    fUpwindQuad_l[50] = hyb_3x3v_p1_surfx1_eval_quad_node_50_r(fl); 
  } else { 
    fUpwindQuad_l[50] = hyb_3x3v_p1_surfx1_eval_quad_node_50_l(fc); 
  } 
  if ((-0.2371708245126284*(alphar[20]+alphar[16]))+0.2371708245126284*alphar[13]-0.2371708245126284*alphar[12]+0.2371708245126284*alphar[8]-0.2371708245126284*alphar[7]-0.1767766952966368*alphar[6]+0.2371708245126284*(alphar[5]+alphar[3])+0.1767766952966368*alphar[2]-0.1767766952966368*alphar[1]+0.1767766952966368*alphar[0] > 0) { 
    fUpwindQuad_r[50] = hyb_3x3v_p1_surfx1_eval_quad_node_50_r(fc); 
  } else { 
    fUpwindQuad_r[50] = hyb_3x3v_p1_surfx1_eval_quad_node_50_l(fr); 
  } 
  if (0.2371708245126284*alphal[20]-0.2371708245126284*(alphal[17]+alphal[16]+alphal[13])+0.2371708245126284*(alphal[12]+alphal[10])-0.2371708245126284*alphal[9]+0.2371708245126284*alphal[8]-0.2371708245126284*alphal[7]-0.1767766952966368*alphal[6]-0.2371708245126284*alphal[5]+0.2371708245126284*(alphal[4]+alphal[3])+0.1767766952966368*alphal[2]-0.1767766952966368*alphal[1]+0.1767766952966368*alphal[0] > 0) { 
    fUpwindQuad_l[51] = hyb_3x3v_p1_surfx1_eval_quad_node_51_r(fl); 
  } else { 
    fUpwindQuad_l[51] = hyb_3x3v_p1_surfx1_eval_quad_node_51_l(fc); 
  } 
  if (0.2371708245126284*alphar[20]-0.2371708245126284*(alphar[17]+alphar[16]+alphar[13])+0.2371708245126284*(alphar[12]+alphar[10])-0.2371708245126284*alphar[9]+0.2371708245126284*alphar[8]-0.2371708245126284*alphar[7]-0.1767766952966368*alphar[6]-0.2371708245126284*alphar[5]+0.2371708245126284*(alphar[4]+alphar[3])+0.1767766952966368*alphar[2]-0.1767766952966368*alphar[1]+0.1767766952966368*alphar[0] > 0) { 
    fUpwindQuad_r[51] = hyb_3x3v_p1_surfx1_eval_quad_node_51_r(fc); 
  } else { 
    fUpwindQuad_r[51] = hyb_3x3v_p1_surfx1_eval_quad_node_51_l(fr); 
  } 
  if ((-0.2371708245126284*(alphal[17]+alphal[16]))+0.2371708245126284*alphal[10]-0.2371708245126284*alphal[9]+0.2371708245126284*alphal[8]-0.2371708245126284*alphal[7]-0.1767766952966368*alphal[6]+0.2371708245126284*(alphal[4]+alphal[3])+0.1767766952966368*alphal[2]-0.1767766952966368*alphal[1]+0.1767766952966368*alphal[0] > 0) { 
    fUpwindQuad_l[52] = hyb_3x3v_p1_surfx1_eval_quad_node_52_r(fl); 
  } else { 
    fUpwindQuad_l[52] = hyb_3x3v_p1_surfx1_eval_quad_node_52_l(fc); 
  } 
  if ((-0.2371708245126284*(alphar[17]+alphar[16]))+0.2371708245126284*alphar[10]-0.2371708245126284*alphar[9]+0.2371708245126284*alphar[8]-0.2371708245126284*alphar[7]-0.1767766952966368*alphar[6]+0.2371708245126284*(alphar[4]+alphar[3])+0.1767766952966368*alphar[2]-0.1767766952966368*alphar[1]+0.1767766952966368*alphar[0] > 0) { 
    fUpwindQuad_r[52] = hyb_3x3v_p1_surfx1_eval_quad_node_52_r(fc); 
  } else { 
    fUpwindQuad_r[52] = hyb_3x3v_p1_surfx1_eval_quad_node_52_l(fr); 
  } 
  if ((-0.2371708245126284*(alphal[20]+alphal[17]+alphal[16]))+0.2371708245126284*alphal[13]-0.2371708245126284*alphal[12]+0.2371708245126284*alphal[10]-0.2371708245126284*alphal[9]+0.2371708245126284*alphal[8]-0.2371708245126284*alphal[7]-0.1767766952966368*alphal[6]+0.2371708245126284*(alphal[5]+alphal[4]+alphal[3])+0.1767766952966368*alphal[2]-0.1767766952966368*alphal[1]+0.1767766952966368*alphal[0] > 0) { 
    fUpwindQuad_l[53] = hyb_3x3v_p1_surfx1_eval_quad_node_53_r(fl); 
  } else { 
    fUpwindQuad_l[53] = hyb_3x3v_p1_surfx1_eval_quad_node_53_l(fc); 
  } 
  if ((-0.2371708245126284*(alphar[20]+alphar[17]+alphar[16]))+0.2371708245126284*alphar[13]-0.2371708245126284*alphar[12]+0.2371708245126284*alphar[10]-0.2371708245126284*alphar[9]+0.2371708245126284*alphar[8]-0.2371708245126284*alphar[7]-0.1767766952966368*alphar[6]+0.2371708245126284*(alphar[5]+alphar[4]+alphar[3])+0.1767766952966368*alphar[2]-0.1767766952966368*alphar[1]+0.1767766952966368*alphar[0] > 0) { 
    fUpwindQuad_r[53] = hyb_3x3v_p1_surfx1_eval_quad_node_53_r(fc); 
  } else { 
    fUpwindQuad_r[53] = hyb_3x3v_p1_surfx1_eval_quad_node_53_l(fr); 
  } 
  if (0.2371708245126284*(alphal[20]+alphal[17]+alphal[16]+alphal[13])-0.2371708245126284*alphal[12]+0.2371708245126284*alphal[10]-0.2371708245126284*alphal[9]+0.2371708245126284*alphal[8]-0.2371708245126284*alphal[7]-0.1767766952966368*alphal[6]-0.2371708245126284*(alphal[5]+alphal[4]+alphal[3])-0.1767766952966368*alphal[2]+0.1767766952966368*(alphal[1]+alphal[0]) > 0) { 
    fUpwindQuad_l[54] = hyb_3x3v_p1_surfx1_eval_quad_node_54_r(fl); 
  } else { 
    fUpwindQuad_l[54] = hyb_3x3v_p1_surfx1_eval_quad_node_54_l(fc); 
  } 
  if (0.2371708245126284*(alphar[20]+alphar[17]+alphar[16]+alphar[13])-0.2371708245126284*alphar[12]+0.2371708245126284*alphar[10]-0.2371708245126284*alphar[9]+0.2371708245126284*alphar[8]-0.2371708245126284*alphar[7]-0.1767766952966368*alphar[6]-0.2371708245126284*(alphar[5]+alphar[4]+alphar[3])-0.1767766952966368*alphar[2]+0.1767766952966368*(alphar[1]+alphar[0]) > 0) { 
    fUpwindQuad_r[54] = hyb_3x3v_p1_surfx1_eval_quad_node_54_r(fc); 
  } else { 
    fUpwindQuad_r[54] = hyb_3x3v_p1_surfx1_eval_quad_node_54_l(fr); 
  } 
  if (0.2371708245126284*(alphal[17]+alphal[16]+alphal[10])-0.2371708245126284*alphal[9]+0.2371708245126284*alphal[8]-0.2371708245126284*alphal[7]-0.1767766952966368*alphal[6]-0.2371708245126284*(alphal[4]+alphal[3])-0.1767766952966368*alphal[2]+0.1767766952966368*(alphal[1]+alphal[0]) > 0) { 
    fUpwindQuad_l[55] = hyb_3x3v_p1_surfx1_eval_quad_node_55_r(fl); 
  } else { 
    fUpwindQuad_l[55] = hyb_3x3v_p1_surfx1_eval_quad_node_55_l(fc); 
  } 
  if (0.2371708245126284*(alphar[17]+alphar[16]+alphar[10])-0.2371708245126284*alphar[9]+0.2371708245126284*alphar[8]-0.2371708245126284*alphar[7]-0.1767766952966368*alphar[6]-0.2371708245126284*(alphar[4]+alphar[3])-0.1767766952966368*alphar[2]+0.1767766952966368*(alphar[1]+alphar[0]) > 0) { 
    fUpwindQuad_r[55] = hyb_3x3v_p1_surfx1_eval_quad_node_55_r(fc); 
  } else { 
    fUpwindQuad_r[55] = hyb_3x3v_p1_surfx1_eval_quad_node_55_l(fr); 
  } 
  if ((-0.2371708245126284*alphal[20])+0.2371708245126284*(alphal[17]+alphal[16])-0.2371708245126284*alphal[13]+0.2371708245126284*(alphal[12]+alphal[10])-0.2371708245126284*alphal[9]+0.2371708245126284*alphal[8]-0.2371708245126284*alphal[7]-0.1767766952966368*alphal[6]+0.2371708245126284*alphal[5]-0.2371708245126284*(alphal[4]+alphal[3])-0.1767766952966368*alphal[2]+0.1767766952966368*(alphal[1]+alphal[0]) > 0) { 
    fUpwindQuad_l[56] = hyb_3x3v_p1_surfx1_eval_quad_node_56_r(fl); 
  } else { 
    fUpwindQuad_l[56] = hyb_3x3v_p1_surfx1_eval_quad_node_56_l(fc); 
  } 
  if ((-0.2371708245126284*alphar[20])+0.2371708245126284*(alphar[17]+alphar[16])-0.2371708245126284*alphar[13]+0.2371708245126284*(alphar[12]+alphar[10])-0.2371708245126284*alphar[9]+0.2371708245126284*alphar[8]-0.2371708245126284*alphar[7]-0.1767766952966368*alphar[6]+0.2371708245126284*alphar[5]-0.2371708245126284*(alphar[4]+alphar[3])-0.1767766952966368*alphar[2]+0.1767766952966368*(alphar[1]+alphar[0]) > 0) { 
    fUpwindQuad_r[56] = hyb_3x3v_p1_surfx1_eval_quad_node_56_r(fc); 
  } else { 
    fUpwindQuad_r[56] = hyb_3x3v_p1_surfx1_eval_quad_node_56_l(fr); 
  } 
  if (0.2371708245126284*(alphal[20]+alphal[16]+alphal[13])-0.2371708245126284*alphal[12]+0.2371708245126284*alphal[8]-0.2371708245126284*alphal[7]-0.1767766952966368*alphal[6]-0.2371708245126284*(alphal[5]+alphal[3])-0.1767766952966368*alphal[2]+0.1767766952966368*(alphal[1]+alphal[0]) > 0) { 
    fUpwindQuad_l[57] = hyb_3x3v_p1_surfx1_eval_quad_node_57_r(fl); 
  } else { 
    fUpwindQuad_l[57] = hyb_3x3v_p1_surfx1_eval_quad_node_57_l(fc); 
  } 
  if (0.2371708245126284*(alphar[20]+alphar[16]+alphar[13])-0.2371708245126284*alphar[12]+0.2371708245126284*alphar[8]-0.2371708245126284*alphar[7]-0.1767766952966368*alphar[6]-0.2371708245126284*(alphar[5]+alphar[3])-0.1767766952966368*alphar[2]+0.1767766952966368*(alphar[1]+alphar[0]) > 0) { 
    fUpwindQuad_r[57] = hyb_3x3v_p1_surfx1_eval_quad_node_57_r(fc); 
  } else { 
    fUpwindQuad_r[57] = hyb_3x3v_p1_surfx1_eval_quad_node_57_l(fr); 
  } 
  if (0.2371708245126284*(alphal[16]+alphal[8])-0.2371708245126284*alphal[7]-0.1767766952966368*alphal[6]-0.2371708245126284*alphal[3]-0.1767766952966368*alphal[2]+0.1767766952966368*(alphal[1]+alphal[0]) > 0) { 
    fUpwindQuad_l[58] = hyb_3x3v_p1_surfx1_eval_quad_node_58_r(fl); 
  } else { 
    fUpwindQuad_l[58] = hyb_3x3v_p1_surfx1_eval_quad_node_58_l(fc); 
  } 
  if (0.2371708245126284*(alphar[16]+alphar[8])-0.2371708245126284*alphar[7]-0.1767766952966368*alphar[6]-0.2371708245126284*alphar[3]-0.1767766952966368*alphar[2]+0.1767766952966368*(alphar[1]+alphar[0]) > 0) { 
    fUpwindQuad_r[58] = hyb_3x3v_p1_surfx1_eval_quad_node_58_r(fc); 
  } else { 
    fUpwindQuad_r[58] = hyb_3x3v_p1_surfx1_eval_quad_node_58_l(fr); 
  } 
  if ((-0.2371708245126284*alphal[20])+0.2371708245126284*alphal[16]-0.2371708245126284*alphal[13]+0.2371708245126284*(alphal[12]+alphal[8])-0.2371708245126284*alphal[7]-0.1767766952966368*alphal[6]+0.2371708245126284*alphal[5]-0.2371708245126284*alphal[3]-0.1767766952966368*alphal[2]+0.1767766952966368*(alphal[1]+alphal[0]) > 0) { 
    fUpwindQuad_l[59] = hyb_3x3v_p1_surfx1_eval_quad_node_59_r(fl); 
  } else { 
    fUpwindQuad_l[59] = hyb_3x3v_p1_surfx1_eval_quad_node_59_l(fc); 
  } 
  if ((-0.2371708245126284*alphar[20])+0.2371708245126284*alphar[16]-0.2371708245126284*alphar[13]+0.2371708245126284*(alphar[12]+alphar[8])-0.2371708245126284*alphar[7]-0.1767766952966368*alphar[6]+0.2371708245126284*alphar[5]-0.2371708245126284*alphar[3]-0.1767766952966368*alphar[2]+0.1767766952966368*(alphar[1]+alphar[0]) > 0) { 
    fUpwindQuad_r[59] = hyb_3x3v_p1_surfx1_eval_quad_node_59_r(fc); 
  } else { 
    fUpwindQuad_r[59] = hyb_3x3v_p1_surfx1_eval_quad_node_59_l(fr); 
  } 
  if (0.2371708245126284*alphal[20]-0.2371708245126284*alphal[17]+0.2371708245126284*(alphal[16]+alphal[13])-0.2371708245126284*(alphal[12]+alphal[10])+0.2371708245126284*(alphal[9]+alphal[8])-0.2371708245126284*alphal[7]-0.1767766952966368*alphal[6]-0.2371708245126284*alphal[5]+0.2371708245126284*alphal[4]-0.2371708245126284*alphal[3]-0.1767766952966368*alphal[2]+0.1767766952966368*(alphal[1]+alphal[0]) > 0) { 
    fUpwindQuad_l[60] = hyb_3x3v_p1_surfx1_eval_quad_node_60_r(fl); 
  } else { 
    fUpwindQuad_l[60] = hyb_3x3v_p1_surfx1_eval_quad_node_60_l(fc); 
  } 
  if (0.2371708245126284*alphar[20]-0.2371708245126284*alphar[17]+0.2371708245126284*(alphar[16]+alphar[13])-0.2371708245126284*(alphar[12]+alphar[10])+0.2371708245126284*(alphar[9]+alphar[8])-0.2371708245126284*alphar[7]-0.1767766952966368*alphar[6]-0.2371708245126284*alphar[5]+0.2371708245126284*alphar[4]-0.2371708245126284*alphar[3]-0.1767766952966368*alphar[2]+0.1767766952966368*(alphar[1]+alphar[0]) > 0) { 
    fUpwindQuad_r[60] = hyb_3x3v_p1_surfx1_eval_quad_node_60_r(fc); 
  } else { 
    fUpwindQuad_r[60] = hyb_3x3v_p1_surfx1_eval_quad_node_60_l(fr); 
  } 
  if ((-0.2371708245126284*alphal[17])+0.2371708245126284*alphal[16]-0.2371708245126284*alphal[10]+0.2371708245126284*(alphal[9]+alphal[8])-0.2371708245126284*alphal[7]-0.1767766952966368*alphal[6]+0.2371708245126284*alphal[4]-0.2371708245126284*alphal[3]-0.1767766952966368*alphal[2]+0.1767766952966368*(alphal[1]+alphal[0]) > 0) { 
    fUpwindQuad_l[61] = hyb_3x3v_p1_surfx1_eval_quad_node_61_r(fl); 
  } else { 
    fUpwindQuad_l[61] = hyb_3x3v_p1_surfx1_eval_quad_node_61_l(fc); 
  } 
  if ((-0.2371708245126284*alphar[17])+0.2371708245126284*alphar[16]-0.2371708245126284*alphar[10]+0.2371708245126284*(alphar[9]+alphar[8])-0.2371708245126284*alphar[7]-0.1767766952966368*alphar[6]+0.2371708245126284*alphar[4]-0.2371708245126284*alphar[3]-0.1767766952966368*alphar[2]+0.1767766952966368*(alphar[1]+alphar[0]) > 0) { 
    fUpwindQuad_r[61] = hyb_3x3v_p1_surfx1_eval_quad_node_61_r(fc); 
  } else { 
    fUpwindQuad_r[61] = hyb_3x3v_p1_surfx1_eval_quad_node_61_l(fr); 
  } 
  if ((-0.2371708245126284*(alphal[20]+alphal[17]))+0.2371708245126284*alphal[16]-0.2371708245126284*alphal[13]+0.2371708245126284*alphal[12]-0.2371708245126284*alphal[10]+0.2371708245126284*(alphal[9]+alphal[8])-0.2371708245126284*alphal[7]-0.1767766952966368*alphal[6]+0.2371708245126284*(alphal[5]+alphal[4])-0.2371708245126284*alphal[3]-0.1767766952966368*alphal[2]+0.1767766952966368*(alphal[1]+alphal[0]) > 0) { 
    fUpwindQuad_l[62] = hyb_3x3v_p1_surfx1_eval_quad_node_62_r(fl); 
  } else { 
    fUpwindQuad_l[62] = hyb_3x3v_p1_surfx1_eval_quad_node_62_l(fc); 
  } 
  if ((-0.2371708245126284*(alphar[20]+alphar[17]))+0.2371708245126284*alphar[16]-0.2371708245126284*alphar[13]+0.2371708245126284*alphar[12]-0.2371708245126284*alphar[10]+0.2371708245126284*(alphar[9]+alphar[8])-0.2371708245126284*alphar[7]-0.1767766952966368*alphar[6]+0.2371708245126284*(alphar[5]+alphar[4])-0.2371708245126284*alphar[3]-0.1767766952966368*alphar[2]+0.1767766952966368*(alphar[1]+alphar[0]) > 0) { 
    fUpwindQuad_r[62] = hyb_3x3v_p1_surfx1_eval_quad_node_62_r(fc); 
  } else { 
    fUpwindQuad_r[62] = hyb_3x3v_p1_surfx1_eval_quad_node_62_l(fr); 
  } 
  if (0.2371708245126284*(alphal[20]+alphal[17]+alphal[13])-0.2371708245126284*alphal[12]+0.2371708245126284*alphal[10]-0.2371708245126284*alphal[9]-0.1767766952966368*alphal[6]-0.2371708245126284*(alphal[5]+alphal[4])-0.1767766952966368*alphal[2]+0.1767766952966368*(alphal[1]+alphal[0]) > 0) { 
    fUpwindQuad_l[63] = hyb_3x3v_p1_surfx1_eval_quad_node_63_r(fl); 
  } else { 
    fUpwindQuad_l[63] = hyb_3x3v_p1_surfx1_eval_quad_node_63_l(fc); 
  } 
  if (0.2371708245126284*(alphar[20]+alphar[17]+alphar[13])-0.2371708245126284*alphar[12]+0.2371708245126284*alphar[10]-0.2371708245126284*alphar[9]-0.1767766952966368*alphar[6]-0.2371708245126284*(alphar[5]+alphar[4])-0.1767766952966368*alphar[2]+0.1767766952966368*(alphar[1]+alphar[0]) > 0) { 
    fUpwindQuad_r[63] = hyb_3x3v_p1_surfx1_eval_quad_node_63_r(fc); 
  } else { 
    fUpwindQuad_r[63] = hyb_3x3v_p1_surfx1_eval_quad_node_63_l(fr); 
  } 
  if (0.2371708245126284*(alphal[17]+alphal[10])-0.2371708245126284*alphal[9]-0.1767766952966368*alphal[6]-0.2371708245126284*alphal[4]-0.1767766952966368*alphal[2]+0.1767766952966368*(alphal[1]+alphal[0]) > 0) { 
    fUpwindQuad_l[64] = hyb_3x3v_p1_surfx1_eval_quad_node_64_r(fl); 
  } else { 
    fUpwindQuad_l[64] = hyb_3x3v_p1_surfx1_eval_quad_node_64_l(fc); 
  } 
  if (0.2371708245126284*(alphar[17]+alphar[10])-0.2371708245126284*alphar[9]-0.1767766952966368*alphar[6]-0.2371708245126284*alphar[4]-0.1767766952966368*alphar[2]+0.1767766952966368*(alphar[1]+alphar[0]) > 0) { 
    fUpwindQuad_r[64] = hyb_3x3v_p1_surfx1_eval_quad_node_64_r(fc); 
  } else { 
    fUpwindQuad_r[64] = hyb_3x3v_p1_surfx1_eval_quad_node_64_l(fr); 
  } 
  if ((-0.2371708245126284*alphal[20])+0.2371708245126284*alphal[17]-0.2371708245126284*alphal[13]+0.2371708245126284*(alphal[12]+alphal[10])-0.2371708245126284*alphal[9]-0.1767766952966368*alphal[6]+0.2371708245126284*alphal[5]-0.2371708245126284*alphal[4]-0.1767766952966368*alphal[2]+0.1767766952966368*(alphal[1]+alphal[0]) > 0) { 
    fUpwindQuad_l[65] = hyb_3x3v_p1_surfx1_eval_quad_node_65_r(fl); 
  } else { 
    fUpwindQuad_l[65] = hyb_3x3v_p1_surfx1_eval_quad_node_65_l(fc); 
  } 
  if ((-0.2371708245126284*alphar[20])+0.2371708245126284*alphar[17]-0.2371708245126284*alphar[13]+0.2371708245126284*(alphar[12]+alphar[10])-0.2371708245126284*alphar[9]-0.1767766952966368*alphar[6]+0.2371708245126284*alphar[5]-0.2371708245126284*alphar[4]-0.1767766952966368*alphar[2]+0.1767766952966368*(alphar[1]+alphar[0]) > 0) { 
    fUpwindQuad_r[65] = hyb_3x3v_p1_surfx1_eval_quad_node_65_r(fc); 
  } else { 
    fUpwindQuad_r[65] = hyb_3x3v_p1_surfx1_eval_quad_node_65_l(fr); 
  } 
  if (0.2371708245126284*(alphal[20]+alphal[13])-0.2371708245126284*alphal[12]-0.1767766952966368*alphal[6]-0.2371708245126284*alphal[5]-0.1767766952966368*alphal[2]+0.1767766952966368*(alphal[1]+alphal[0]) > 0) { 
    fUpwindQuad_l[66] = hyb_3x3v_p1_surfx1_eval_quad_node_66_r(fl); 
  } else { 
    fUpwindQuad_l[66] = hyb_3x3v_p1_surfx1_eval_quad_node_66_l(fc); 
  } 
  if (0.2371708245126284*(alphar[20]+alphar[13])-0.2371708245126284*alphar[12]-0.1767766952966368*alphar[6]-0.2371708245126284*alphar[5]-0.1767766952966368*alphar[2]+0.1767766952966368*(alphar[1]+alphar[0]) > 0) { 
    fUpwindQuad_r[66] = hyb_3x3v_p1_surfx1_eval_quad_node_66_r(fc); 
  } else { 
    fUpwindQuad_r[66] = hyb_3x3v_p1_surfx1_eval_quad_node_66_l(fr); 
  } 
  if (0.1767766952966368*(alphal[1]+alphal[0])-0.1767766952966368*(alphal[6]+alphal[2]) > 0) { 
    fUpwindQuad_l[67] = hyb_3x3v_p1_surfx1_eval_quad_node_67_r(fl); 
  } else { 
    fUpwindQuad_l[67] = hyb_3x3v_p1_surfx1_eval_quad_node_67_l(fc); 
  } 
  if (0.1767766952966368*(alphar[1]+alphar[0])-0.1767766952966368*(alphar[6]+alphar[2]) > 0) { 
    fUpwindQuad_r[67] = hyb_3x3v_p1_surfx1_eval_quad_node_67_r(fc); 
  } else { 
    fUpwindQuad_r[67] = hyb_3x3v_p1_surfx1_eval_quad_node_67_l(fr); 
  } 
  if ((-0.2371708245126284*(alphal[20]+alphal[13]))+0.2371708245126284*alphal[12]-0.1767766952966368*alphal[6]+0.2371708245126284*alphal[5]-0.1767766952966368*alphal[2]+0.1767766952966368*(alphal[1]+alphal[0]) > 0) { 
    fUpwindQuad_l[68] = hyb_3x3v_p1_surfx1_eval_quad_node_68_r(fl); 
  } else { 
    fUpwindQuad_l[68] = hyb_3x3v_p1_surfx1_eval_quad_node_68_l(fc); 
  } 
  if ((-0.2371708245126284*(alphar[20]+alphar[13]))+0.2371708245126284*alphar[12]-0.1767766952966368*alphar[6]+0.2371708245126284*alphar[5]-0.1767766952966368*alphar[2]+0.1767766952966368*(alphar[1]+alphar[0]) > 0) { 
    fUpwindQuad_r[68] = hyb_3x3v_p1_surfx1_eval_quad_node_68_r(fc); 
  } else { 
    fUpwindQuad_r[68] = hyb_3x3v_p1_surfx1_eval_quad_node_68_l(fr); 
  } 
  if (0.2371708245126284*alphal[20]-0.2371708245126284*alphal[17]+0.2371708245126284*alphal[13]-0.2371708245126284*(alphal[12]+alphal[10])+0.2371708245126284*alphal[9]-0.1767766952966368*alphal[6]-0.2371708245126284*alphal[5]+0.2371708245126284*alphal[4]-0.1767766952966368*alphal[2]+0.1767766952966368*(alphal[1]+alphal[0]) > 0) { 
    fUpwindQuad_l[69] = hyb_3x3v_p1_surfx1_eval_quad_node_69_r(fl); 
  } else { 
    fUpwindQuad_l[69] = hyb_3x3v_p1_surfx1_eval_quad_node_69_l(fc); 
  } 
  if (0.2371708245126284*alphar[20]-0.2371708245126284*alphar[17]+0.2371708245126284*alphar[13]-0.2371708245126284*(alphar[12]+alphar[10])+0.2371708245126284*alphar[9]-0.1767766952966368*alphar[6]-0.2371708245126284*alphar[5]+0.2371708245126284*alphar[4]-0.1767766952966368*alphar[2]+0.1767766952966368*(alphar[1]+alphar[0]) > 0) { 
    fUpwindQuad_r[69] = hyb_3x3v_p1_surfx1_eval_quad_node_69_r(fc); 
  } else { 
    fUpwindQuad_r[69] = hyb_3x3v_p1_surfx1_eval_quad_node_69_l(fr); 
  } 
  if ((-0.2371708245126284*(alphal[17]+alphal[10]))+0.2371708245126284*alphal[9]-0.1767766952966368*alphal[6]+0.2371708245126284*alphal[4]-0.1767766952966368*alphal[2]+0.1767766952966368*(alphal[1]+alphal[0]) > 0) { 
    fUpwindQuad_l[70] = hyb_3x3v_p1_surfx1_eval_quad_node_70_r(fl); 
  } else { 
    fUpwindQuad_l[70] = hyb_3x3v_p1_surfx1_eval_quad_node_70_l(fc); 
  } 
  if ((-0.2371708245126284*(alphar[17]+alphar[10]))+0.2371708245126284*alphar[9]-0.1767766952966368*alphar[6]+0.2371708245126284*alphar[4]-0.1767766952966368*alphar[2]+0.1767766952966368*(alphar[1]+alphar[0]) > 0) { 
    fUpwindQuad_r[70] = hyb_3x3v_p1_surfx1_eval_quad_node_70_r(fc); 
  } else { 
    fUpwindQuad_r[70] = hyb_3x3v_p1_surfx1_eval_quad_node_70_l(fr); 
  } 
  if ((-0.2371708245126284*(alphal[20]+alphal[17]+alphal[13]))+0.2371708245126284*alphal[12]-0.2371708245126284*alphal[10]+0.2371708245126284*alphal[9]-0.1767766952966368*alphal[6]+0.2371708245126284*(alphal[5]+alphal[4])-0.1767766952966368*alphal[2]+0.1767766952966368*(alphal[1]+alphal[0]) > 0) { 
    fUpwindQuad_l[71] = hyb_3x3v_p1_surfx1_eval_quad_node_71_r(fl); 
  } else { 
    fUpwindQuad_l[71] = hyb_3x3v_p1_surfx1_eval_quad_node_71_l(fc); 
  } 
  if ((-0.2371708245126284*(alphar[20]+alphar[17]+alphar[13]))+0.2371708245126284*alphar[12]-0.2371708245126284*alphar[10]+0.2371708245126284*alphar[9]-0.1767766952966368*alphar[6]+0.2371708245126284*(alphar[5]+alphar[4])-0.1767766952966368*alphar[2]+0.1767766952966368*(alphar[1]+alphar[0]) > 0) { 
    fUpwindQuad_r[71] = hyb_3x3v_p1_surfx1_eval_quad_node_71_r(fc); 
  } else { 
    fUpwindQuad_r[71] = hyb_3x3v_p1_surfx1_eval_quad_node_71_l(fr); 
  } 
  if (0.2371708245126284*(alphal[20]+alphal[17])-0.2371708245126284*alphal[16]+0.2371708245126284*alphal[13]-0.2371708245126284*alphal[12]+0.2371708245126284*alphal[10]-0.2371708245126284*(alphal[9]+alphal[8])+0.2371708245126284*alphal[7]-0.1767766952966368*alphal[6]-0.2371708245126284*(alphal[5]+alphal[4])+0.2371708245126284*alphal[3]-0.1767766952966368*alphal[2]+0.1767766952966368*(alphal[1]+alphal[0]) > 0) { 
    fUpwindQuad_l[72] = hyb_3x3v_p1_surfx1_eval_quad_node_72_r(fl); 
  } else { 
    fUpwindQuad_l[72] = hyb_3x3v_p1_surfx1_eval_quad_node_72_l(fc); 
  } 
  if (0.2371708245126284*(alphar[20]+alphar[17])-0.2371708245126284*alphar[16]+0.2371708245126284*alphar[13]-0.2371708245126284*alphar[12]+0.2371708245126284*alphar[10]-0.2371708245126284*(alphar[9]+alphar[8])+0.2371708245126284*alphar[7]-0.1767766952966368*alphar[6]-0.2371708245126284*(alphar[5]+alphar[4])+0.2371708245126284*alphar[3]-0.1767766952966368*alphar[2]+0.1767766952966368*(alphar[1]+alphar[0]) > 0) { 
    fUpwindQuad_r[72] = hyb_3x3v_p1_surfx1_eval_quad_node_72_r(fc); 
  } else { 
    fUpwindQuad_r[72] = hyb_3x3v_p1_surfx1_eval_quad_node_72_l(fr); 
  } 
  if (0.2371708245126284*alphal[17]-0.2371708245126284*alphal[16]+0.2371708245126284*alphal[10]-0.2371708245126284*(alphal[9]+alphal[8])+0.2371708245126284*alphal[7]-0.1767766952966368*alphal[6]-0.2371708245126284*alphal[4]+0.2371708245126284*alphal[3]-0.1767766952966368*alphal[2]+0.1767766952966368*(alphal[1]+alphal[0]) > 0) { 
    fUpwindQuad_l[73] = hyb_3x3v_p1_surfx1_eval_quad_node_73_r(fl); 
  } else { 
    fUpwindQuad_l[73] = hyb_3x3v_p1_surfx1_eval_quad_node_73_l(fc); 
  } 
  if (0.2371708245126284*alphar[17]-0.2371708245126284*alphar[16]+0.2371708245126284*alphar[10]-0.2371708245126284*(alphar[9]+alphar[8])+0.2371708245126284*alphar[7]-0.1767766952966368*alphar[6]-0.2371708245126284*alphar[4]+0.2371708245126284*alphar[3]-0.1767766952966368*alphar[2]+0.1767766952966368*(alphar[1]+alphar[0]) > 0) { 
    fUpwindQuad_r[73] = hyb_3x3v_p1_surfx1_eval_quad_node_73_r(fc); 
  } else { 
    fUpwindQuad_r[73] = hyb_3x3v_p1_surfx1_eval_quad_node_73_l(fr); 
  } 
  if ((-0.2371708245126284*alphal[20])+0.2371708245126284*alphal[17]-0.2371708245126284*(alphal[16]+alphal[13])+0.2371708245126284*(alphal[12]+alphal[10])-0.2371708245126284*(alphal[9]+alphal[8])+0.2371708245126284*alphal[7]-0.1767766952966368*alphal[6]+0.2371708245126284*alphal[5]-0.2371708245126284*alphal[4]+0.2371708245126284*alphal[3]-0.1767766952966368*alphal[2]+0.1767766952966368*(alphal[1]+alphal[0]) > 0) { 
    fUpwindQuad_l[74] = hyb_3x3v_p1_surfx1_eval_quad_node_74_r(fl); 
  } else { 
    fUpwindQuad_l[74] = hyb_3x3v_p1_surfx1_eval_quad_node_74_l(fc); 
  } 
  if ((-0.2371708245126284*alphar[20])+0.2371708245126284*alphar[17]-0.2371708245126284*(alphar[16]+alphar[13])+0.2371708245126284*(alphar[12]+alphar[10])-0.2371708245126284*(alphar[9]+alphar[8])+0.2371708245126284*alphar[7]-0.1767766952966368*alphar[6]+0.2371708245126284*alphar[5]-0.2371708245126284*alphar[4]+0.2371708245126284*alphar[3]-0.1767766952966368*alphar[2]+0.1767766952966368*(alphar[1]+alphar[0]) > 0) { 
    fUpwindQuad_r[74] = hyb_3x3v_p1_surfx1_eval_quad_node_74_r(fc); 
  } else { 
    fUpwindQuad_r[74] = hyb_3x3v_p1_surfx1_eval_quad_node_74_l(fr); 
  } 
  if (0.2371708245126284*alphal[20]-0.2371708245126284*alphal[16]+0.2371708245126284*alphal[13]-0.2371708245126284*(alphal[12]+alphal[8])+0.2371708245126284*alphal[7]-0.1767766952966368*alphal[6]-0.2371708245126284*alphal[5]+0.2371708245126284*alphal[3]-0.1767766952966368*alphal[2]+0.1767766952966368*(alphal[1]+alphal[0]) > 0) { 
    fUpwindQuad_l[75] = hyb_3x3v_p1_surfx1_eval_quad_node_75_r(fl); 
  } else { 
    fUpwindQuad_l[75] = hyb_3x3v_p1_surfx1_eval_quad_node_75_l(fc); 
  } 
  if (0.2371708245126284*alphar[20]-0.2371708245126284*alphar[16]+0.2371708245126284*alphar[13]-0.2371708245126284*(alphar[12]+alphar[8])+0.2371708245126284*alphar[7]-0.1767766952966368*alphar[6]-0.2371708245126284*alphar[5]+0.2371708245126284*alphar[3]-0.1767766952966368*alphar[2]+0.1767766952966368*(alphar[1]+alphar[0]) > 0) { 
    fUpwindQuad_r[75] = hyb_3x3v_p1_surfx1_eval_quad_node_75_r(fc); 
  } else { 
    fUpwindQuad_r[75] = hyb_3x3v_p1_surfx1_eval_quad_node_75_l(fr); 
  } 
  if ((-0.2371708245126284*(alphal[16]+alphal[8]))+0.2371708245126284*alphal[7]-0.1767766952966368*alphal[6]+0.2371708245126284*alphal[3]-0.1767766952966368*alphal[2]+0.1767766952966368*(alphal[1]+alphal[0]) > 0) { 
    fUpwindQuad_l[76] = hyb_3x3v_p1_surfx1_eval_quad_node_76_r(fl); 
  } else { 
    fUpwindQuad_l[76] = hyb_3x3v_p1_surfx1_eval_quad_node_76_l(fc); 
  } 
  if ((-0.2371708245126284*(alphar[16]+alphar[8]))+0.2371708245126284*alphar[7]-0.1767766952966368*alphar[6]+0.2371708245126284*alphar[3]-0.1767766952966368*alphar[2]+0.1767766952966368*(alphar[1]+alphar[0]) > 0) { 
    fUpwindQuad_r[76] = hyb_3x3v_p1_surfx1_eval_quad_node_76_r(fc); 
  } else { 
    fUpwindQuad_r[76] = hyb_3x3v_p1_surfx1_eval_quad_node_76_l(fr); 
  } 
  if ((-0.2371708245126284*(alphal[20]+alphal[16]+alphal[13]))+0.2371708245126284*alphal[12]-0.2371708245126284*alphal[8]+0.2371708245126284*alphal[7]-0.1767766952966368*alphal[6]+0.2371708245126284*(alphal[5]+alphal[3])-0.1767766952966368*alphal[2]+0.1767766952966368*(alphal[1]+alphal[0]) > 0) { 
    fUpwindQuad_l[77] = hyb_3x3v_p1_surfx1_eval_quad_node_77_r(fl); 
  } else { 
    fUpwindQuad_l[77] = hyb_3x3v_p1_surfx1_eval_quad_node_77_l(fc); 
  } 
  if ((-0.2371708245126284*(alphar[20]+alphar[16]+alphar[13]))+0.2371708245126284*alphar[12]-0.2371708245126284*alphar[8]+0.2371708245126284*alphar[7]-0.1767766952966368*alphar[6]+0.2371708245126284*(alphar[5]+alphar[3])-0.1767766952966368*alphar[2]+0.1767766952966368*(alphar[1]+alphar[0]) > 0) { 
    fUpwindQuad_r[77] = hyb_3x3v_p1_surfx1_eval_quad_node_77_r(fc); 
  } else { 
    fUpwindQuad_r[77] = hyb_3x3v_p1_surfx1_eval_quad_node_77_l(fr); 
  } 
  if (0.2371708245126284*alphal[20]-0.2371708245126284*(alphal[17]+alphal[16])+0.2371708245126284*alphal[13]-0.2371708245126284*(alphal[12]+alphal[10])+0.2371708245126284*alphal[9]-0.2371708245126284*alphal[8]+0.2371708245126284*alphal[7]-0.1767766952966368*alphal[6]-0.2371708245126284*alphal[5]+0.2371708245126284*(alphal[4]+alphal[3])-0.1767766952966368*alphal[2]+0.1767766952966368*(alphal[1]+alphal[0]) > 0) { 
    fUpwindQuad_l[78] = hyb_3x3v_p1_surfx1_eval_quad_node_78_r(fl); 
  } else { 
    fUpwindQuad_l[78] = hyb_3x3v_p1_surfx1_eval_quad_node_78_l(fc); 
  } 
  if (0.2371708245126284*alphar[20]-0.2371708245126284*(alphar[17]+alphar[16])+0.2371708245126284*alphar[13]-0.2371708245126284*(alphar[12]+alphar[10])+0.2371708245126284*alphar[9]-0.2371708245126284*alphar[8]+0.2371708245126284*alphar[7]-0.1767766952966368*alphar[6]-0.2371708245126284*alphar[5]+0.2371708245126284*(alphar[4]+alphar[3])-0.1767766952966368*alphar[2]+0.1767766952966368*(alphar[1]+alphar[0]) > 0) { 
    fUpwindQuad_r[78] = hyb_3x3v_p1_surfx1_eval_quad_node_78_r(fc); 
  } else { 
    fUpwindQuad_r[78] = hyb_3x3v_p1_surfx1_eval_quad_node_78_l(fr); 
  } 
  if ((-0.2371708245126284*(alphal[17]+alphal[16]+alphal[10]))+0.2371708245126284*alphal[9]-0.2371708245126284*alphal[8]+0.2371708245126284*alphal[7]-0.1767766952966368*alphal[6]+0.2371708245126284*(alphal[4]+alphal[3])-0.1767766952966368*alphal[2]+0.1767766952966368*(alphal[1]+alphal[0]) > 0) { 
    fUpwindQuad_l[79] = hyb_3x3v_p1_surfx1_eval_quad_node_79_r(fl); 
  } else { 
    fUpwindQuad_l[79] = hyb_3x3v_p1_surfx1_eval_quad_node_79_l(fc); 
  } 
  if ((-0.2371708245126284*(alphar[17]+alphar[16]+alphar[10]))+0.2371708245126284*alphar[9]-0.2371708245126284*alphar[8]+0.2371708245126284*alphar[7]-0.1767766952966368*alphar[6]+0.2371708245126284*(alphar[4]+alphar[3])-0.1767766952966368*alphar[2]+0.1767766952966368*(alphar[1]+alphar[0]) > 0) { 
    fUpwindQuad_r[79] = hyb_3x3v_p1_surfx1_eval_quad_node_79_r(fc); 
  } else { 
    fUpwindQuad_r[79] = hyb_3x3v_p1_surfx1_eval_quad_node_79_l(fr); 
  } 
  if ((-0.2371708245126284*(alphal[20]+alphal[17]+alphal[16]+alphal[13]))+0.2371708245126284*alphal[12]-0.2371708245126284*alphal[10]+0.2371708245126284*alphal[9]-0.2371708245126284*alphal[8]+0.2371708245126284*alphal[7]-0.1767766952966368*alphal[6]+0.2371708245126284*(alphal[5]+alphal[4]+alphal[3])-0.1767766952966368*alphal[2]+0.1767766952966368*(alphal[1]+alphal[0]) > 0) { 
    fUpwindQuad_l[80] = hyb_3x3v_p1_surfx1_eval_quad_node_80_r(fl); 
  } else { 
    fUpwindQuad_l[80] = hyb_3x3v_p1_surfx1_eval_quad_node_80_l(fc); 
  } 
  if ((-0.2371708245126284*(alphar[20]+alphar[17]+alphar[16]+alphar[13]))+0.2371708245126284*alphar[12]-0.2371708245126284*alphar[10]+0.2371708245126284*alphar[9]-0.2371708245126284*alphar[8]+0.2371708245126284*alphar[7]-0.1767766952966368*alphar[6]+0.2371708245126284*(alphar[5]+alphar[4]+alphar[3])-0.1767766952966368*alphar[2]+0.1767766952966368*(alphar[1]+alphar[0]) > 0) { 
    fUpwindQuad_r[80] = hyb_3x3v_p1_surfx1_eval_quad_node_80_r(fc); 
  } else { 
    fUpwindQuad_r[80] = hyb_3x3v_p1_surfx1_eval_quad_node_80_l(fr); 
  } 
  if ((-0.2371708245126284*(alphal[20]+alphal[17]+alphal[16]+alphal[13]+alphal[12]+alphal[10]+alphal[9]+alphal[8]+alphal[7]))+0.1767766952966368*alphal[6]-0.2371708245126284*(alphal[5]+alphal[4]+alphal[3])+0.1767766952966368*(alphal[2]+alphal[1]+alphal[0]) > 0) { 
    fUpwindQuad_l[81] = hyb_3x3v_p1_surfx1_eval_quad_node_81_r(fl); 
  } else { 
    fUpwindQuad_l[81] = hyb_3x3v_p1_surfx1_eval_quad_node_81_l(fc); 
  } 
  if ((-0.2371708245126284*(alphar[20]+alphar[17]+alphar[16]+alphar[13]+alphar[12]+alphar[10]+alphar[9]+alphar[8]+alphar[7]))+0.1767766952966368*alphar[6]-0.2371708245126284*(alphar[5]+alphar[4]+alphar[3])+0.1767766952966368*(alphar[2]+alphar[1]+alphar[0]) > 0) { 
    fUpwindQuad_r[81] = hyb_3x3v_p1_surfx1_eval_quad_node_81_r(fc); 
  } else { 
    fUpwindQuad_r[81] = hyb_3x3v_p1_surfx1_eval_quad_node_81_l(fr); 
  } 
  if ((-0.2371708245126284*(alphal[17]+alphal[16]+alphal[10]+alphal[9]+alphal[8]+alphal[7]))+0.1767766952966368*alphal[6]-0.2371708245126284*(alphal[4]+alphal[3])+0.1767766952966368*(alphal[2]+alphal[1]+alphal[0]) > 0) { 
    fUpwindQuad_l[82] = hyb_3x3v_p1_surfx1_eval_quad_node_82_r(fl); 
  } else { 
    fUpwindQuad_l[82] = hyb_3x3v_p1_surfx1_eval_quad_node_82_l(fc); 
  } 
  if ((-0.2371708245126284*(alphar[17]+alphar[16]+alphar[10]+alphar[9]+alphar[8]+alphar[7]))+0.1767766952966368*alphar[6]-0.2371708245126284*(alphar[4]+alphar[3])+0.1767766952966368*(alphar[2]+alphar[1]+alphar[0]) > 0) { 
    fUpwindQuad_r[82] = hyb_3x3v_p1_surfx1_eval_quad_node_82_r(fc); 
  } else { 
    fUpwindQuad_r[82] = hyb_3x3v_p1_surfx1_eval_quad_node_82_l(fr); 
  } 
  if (0.2371708245126284*alphal[20]-0.2371708245126284*(alphal[17]+alphal[16])+0.2371708245126284*(alphal[13]+alphal[12])-0.2371708245126284*(alphal[10]+alphal[9]+alphal[8]+alphal[7])+0.1767766952966368*alphal[6]+0.2371708245126284*alphal[5]-0.2371708245126284*(alphal[4]+alphal[3])+0.1767766952966368*(alphal[2]+alphal[1]+alphal[0]) > 0) { 
    fUpwindQuad_l[83] = hyb_3x3v_p1_surfx1_eval_quad_node_83_r(fl); 
  } else { 
    fUpwindQuad_l[83] = hyb_3x3v_p1_surfx1_eval_quad_node_83_l(fc); 
  } 
  if (0.2371708245126284*alphar[20]-0.2371708245126284*(alphar[17]+alphar[16])+0.2371708245126284*(alphar[13]+alphar[12])-0.2371708245126284*(alphar[10]+alphar[9]+alphar[8]+alphar[7])+0.1767766952966368*alphar[6]+0.2371708245126284*alphar[5]-0.2371708245126284*(alphar[4]+alphar[3])+0.1767766952966368*(alphar[2]+alphar[1]+alphar[0]) > 0) { 
    fUpwindQuad_r[83] = hyb_3x3v_p1_surfx1_eval_quad_node_83_r(fc); 
  } else { 
    fUpwindQuad_r[83] = hyb_3x3v_p1_surfx1_eval_quad_node_83_l(fr); 
  } 
  if ((-0.2371708245126284*(alphal[20]+alphal[16]+alphal[13]+alphal[12]+alphal[8]+alphal[7]))+0.1767766952966368*alphal[6]-0.2371708245126284*(alphal[5]+alphal[3])+0.1767766952966368*(alphal[2]+alphal[1]+alphal[0]) > 0) { 
    fUpwindQuad_l[84] = hyb_3x3v_p1_surfx1_eval_quad_node_84_r(fl); 
  } else { 
    fUpwindQuad_l[84] = hyb_3x3v_p1_surfx1_eval_quad_node_84_l(fc); 
  } 
  if ((-0.2371708245126284*(alphar[20]+alphar[16]+alphar[13]+alphar[12]+alphar[8]+alphar[7]))+0.1767766952966368*alphar[6]-0.2371708245126284*(alphar[5]+alphar[3])+0.1767766952966368*(alphar[2]+alphar[1]+alphar[0]) > 0) { 
    fUpwindQuad_r[84] = hyb_3x3v_p1_surfx1_eval_quad_node_84_r(fc); 
  } else { 
    fUpwindQuad_r[84] = hyb_3x3v_p1_surfx1_eval_quad_node_84_l(fr); 
  } 
  if ((-0.2371708245126284*(alphal[16]+alphal[8]+alphal[7]))+0.1767766952966368*alphal[6]-0.2371708245126284*alphal[3]+0.1767766952966368*(alphal[2]+alphal[1]+alphal[0]) > 0) { 
    fUpwindQuad_l[85] = hyb_3x3v_p1_surfx1_eval_quad_node_85_r(fl); 
  } else { 
    fUpwindQuad_l[85] = hyb_3x3v_p1_surfx1_eval_quad_node_85_l(fc); 
  } 
  if ((-0.2371708245126284*(alphar[16]+alphar[8]+alphar[7]))+0.1767766952966368*alphar[6]-0.2371708245126284*alphar[3]+0.1767766952966368*(alphar[2]+alphar[1]+alphar[0]) > 0) { 
    fUpwindQuad_r[85] = hyb_3x3v_p1_surfx1_eval_quad_node_85_r(fc); 
  } else { 
    fUpwindQuad_r[85] = hyb_3x3v_p1_surfx1_eval_quad_node_85_l(fr); 
  } 
  if (0.2371708245126284*alphal[20]-0.2371708245126284*alphal[16]+0.2371708245126284*(alphal[13]+alphal[12])-0.2371708245126284*(alphal[8]+alphal[7])+0.1767766952966368*alphal[6]+0.2371708245126284*alphal[5]-0.2371708245126284*alphal[3]+0.1767766952966368*(alphal[2]+alphal[1]+alphal[0]) > 0) { 
    fUpwindQuad_l[86] = hyb_3x3v_p1_surfx1_eval_quad_node_86_r(fl); 
  } else { 
    fUpwindQuad_l[86] = hyb_3x3v_p1_surfx1_eval_quad_node_86_l(fc); 
  } 
  if (0.2371708245126284*alphar[20]-0.2371708245126284*alphar[16]+0.2371708245126284*(alphar[13]+alphar[12])-0.2371708245126284*(alphar[8]+alphar[7])+0.1767766952966368*alphar[6]+0.2371708245126284*alphar[5]-0.2371708245126284*alphar[3]+0.1767766952966368*(alphar[2]+alphar[1]+alphar[0]) > 0) { 
    fUpwindQuad_r[86] = hyb_3x3v_p1_surfx1_eval_quad_node_86_r(fc); 
  } else { 
    fUpwindQuad_r[86] = hyb_3x3v_p1_surfx1_eval_quad_node_86_l(fr); 
  } 
  if ((-0.2371708245126284*alphal[20])+0.2371708245126284*alphal[17]-0.2371708245126284*(alphal[16]+alphal[13]+alphal[12])+0.2371708245126284*(alphal[10]+alphal[9])-0.2371708245126284*(alphal[8]+alphal[7])+0.1767766952966368*alphal[6]-0.2371708245126284*alphal[5]+0.2371708245126284*alphal[4]-0.2371708245126284*alphal[3]+0.1767766952966368*(alphal[2]+alphal[1]+alphal[0]) > 0) { 
    fUpwindQuad_l[87] = hyb_3x3v_p1_surfx1_eval_quad_node_87_r(fl); 
  } else { 
    fUpwindQuad_l[87] = hyb_3x3v_p1_surfx1_eval_quad_node_87_l(fc); 
  } 
  if ((-0.2371708245126284*alphar[20])+0.2371708245126284*alphar[17]-0.2371708245126284*(alphar[16]+alphar[13]+alphar[12])+0.2371708245126284*(alphar[10]+alphar[9])-0.2371708245126284*(alphar[8]+alphar[7])+0.1767766952966368*alphar[6]-0.2371708245126284*alphar[5]+0.2371708245126284*alphar[4]-0.2371708245126284*alphar[3]+0.1767766952966368*(alphar[2]+alphar[1]+alphar[0]) > 0) { 
    fUpwindQuad_r[87] = hyb_3x3v_p1_surfx1_eval_quad_node_87_r(fc); 
  } else { 
    fUpwindQuad_r[87] = hyb_3x3v_p1_surfx1_eval_quad_node_87_l(fr); 
  } 
  if (0.2371708245126284*alphal[17]-0.2371708245126284*alphal[16]+0.2371708245126284*(alphal[10]+alphal[9])-0.2371708245126284*(alphal[8]+alphal[7])+0.1767766952966368*alphal[6]+0.2371708245126284*alphal[4]-0.2371708245126284*alphal[3]+0.1767766952966368*(alphal[2]+alphal[1]+alphal[0]) > 0) { 
    fUpwindQuad_l[88] = hyb_3x3v_p1_surfx1_eval_quad_node_88_r(fl); 
  } else { 
    fUpwindQuad_l[88] = hyb_3x3v_p1_surfx1_eval_quad_node_88_l(fc); 
  } 
  if (0.2371708245126284*alphar[17]-0.2371708245126284*alphar[16]+0.2371708245126284*(alphar[10]+alphar[9])-0.2371708245126284*(alphar[8]+alphar[7])+0.1767766952966368*alphar[6]+0.2371708245126284*alphar[4]-0.2371708245126284*alphar[3]+0.1767766952966368*(alphar[2]+alphar[1]+alphar[0]) > 0) { 
    fUpwindQuad_r[88] = hyb_3x3v_p1_surfx1_eval_quad_node_88_r(fc); 
  } else { 
    fUpwindQuad_r[88] = hyb_3x3v_p1_surfx1_eval_quad_node_88_l(fr); 
  } 
  if (0.2371708245126284*(alphal[20]+alphal[17])-0.2371708245126284*alphal[16]+0.2371708245126284*(alphal[13]+alphal[12]+alphal[10]+alphal[9])-0.2371708245126284*(alphal[8]+alphal[7])+0.1767766952966368*alphal[6]+0.2371708245126284*(alphal[5]+alphal[4])-0.2371708245126284*alphal[3]+0.1767766952966368*(alphal[2]+alphal[1]+alphal[0]) > 0) { 
    fUpwindQuad_l[89] = hyb_3x3v_p1_surfx1_eval_quad_node_89_r(fl); 
  } else { 
    fUpwindQuad_l[89] = hyb_3x3v_p1_surfx1_eval_quad_node_89_l(fc); 
  } 
  if (0.2371708245126284*(alphar[20]+alphar[17])-0.2371708245126284*alphar[16]+0.2371708245126284*(alphar[13]+alphar[12]+alphar[10]+alphar[9])-0.2371708245126284*(alphar[8]+alphar[7])+0.1767766952966368*alphar[6]+0.2371708245126284*(alphar[5]+alphar[4])-0.2371708245126284*alphar[3]+0.1767766952966368*(alphar[2]+alphar[1]+alphar[0]) > 0) { 
    fUpwindQuad_r[89] = hyb_3x3v_p1_surfx1_eval_quad_node_89_r(fc); 
  } else { 
    fUpwindQuad_r[89] = hyb_3x3v_p1_surfx1_eval_quad_node_89_l(fr); 
  } 
  if ((-0.2371708245126284*(alphal[20]+alphal[17]+alphal[13]+alphal[12]+alphal[10]+alphal[9]))+0.1767766952966368*alphal[6]-0.2371708245126284*(alphal[5]+alphal[4])+0.1767766952966368*(alphal[2]+alphal[1]+alphal[0]) > 0) { 
    fUpwindQuad_l[90] = hyb_3x3v_p1_surfx1_eval_quad_node_90_r(fl); 
  } else { 
    fUpwindQuad_l[90] = hyb_3x3v_p1_surfx1_eval_quad_node_90_l(fc); 
  } 
  if ((-0.2371708245126284*(alphar[20]+alphar[17]+alphar[13]+alphar[12]+alphar[10]+alphar[9]))+0.1767766952966368*alphar[6]-0.2371708245126284*(alphar[5]+alphar[4])+0.1767766952966368*(alphar[2]+alphar[1]+alphar[0]) > 0) { 
    fUpwindQuad_r[90] = hyb_3x3v_p1_surfx1_eval_quad_node_90_r(fc); 
  } else { 
    fUpwindQuad_r[90] = hyb_3x3v_p1_surfx1_eval_quad_node_90_l(fr); 
  } 
  if ((-0.2371708245126284*(alphal[17]+alphal[10]+alphal[9]))+0.1767766952966368*alphal[6]-0.2371708245126284*alphal[4]+0.1767766952966368*(alphal[2]+alphal[1]+alphal[0]) > 0) { 
    fUpwindQuad_l[91] = hyb_3x3v_p1_surfx1_eval_quad_node_91_r(fl); 
  } else { 
    fUpwindQuad_l[91] = hyb_3x3v_p1_surfx1_eval_quad_node_91_l(fc); 
  } 
  if ((-0.2371708245126284*(alphar[17]+alphar[10]+alphar[9]))+0.1767766952966368*alphar[6]-0.2371708245126284*alphar[4]+0.1767766952966368*(alphar[2]+alphar[1]+alphar[0]) > 0) { 
    fUpwindQuad_r[91] = hyb_3x3v_p1_surfx1_eval_quad_node_91_r(fc); 
  } else { 
    fUpwindQuad_r[91] = hyb_3x3v_p1_surfx1_eval_quad_node_91_l(fr); 
  } 
  if (0.2371708245126284*alphal[20]-0.2371708245126284*alphal[17]+0.2371708245126284*(alphal[13]+alphal[12])-0.2371708245126284*(alphal[10]+alphal[9])+0.1767766952966368*alphal[6]+0.2371708245126284*alphal[5]-0.2371708245126284*alphal[4]+0.1767766952966368*(alphal[2]+alphal[1]+alphal[0]) > 0) { 
    fUpwindQuad_l[92] = hyb_3x3v_p1_surfx1_eval_quad_node_92_r(fl); 
  } else { 
    fUpwindQuad_l[92] = hyb_3x3v_p1_surfx1_eval_quad_node_92_l(fc); 
  } 
  if (0.2371708245126284*alphar[20]-0.2371708245126284*alphar[17]+0.2371708245126284*(alphar[13]+alphar[12])-0.2371708245126284*(alphar[10]+alphar[9])+0.1767766952966368*alphar[6]+0.2371708245126284*alphar[5]-0.2371708245126284*alphar[4]+0.1767766952966368*(alphar[2]+alphar[1]+alphar[0]) > 0) { 
    fUpwindQuad_r[92] = hyb_3x3v_p1_surfx1_eval_quad_node_92_r(fc); 
  } else { 
    fUpwindQuad_r[92] = hyb_3x3v_p1_surfx1_eval_quad_node_92_l(fr); 
  } 
  if ((-0.2371708245126284*(alphal[20]+alphal[13]+alphal[12]))+0.1767766952966368*alphal[6]-0.2371708245126284*alphal[5]+0.1767766952966368*(alphal[2]+alphal[1]+alphal[0]) > 0) { 
    fUpwindQuad_l[93] = hyb_3x3v_p1_surfx1_eval_quad_node_93_r(fl); 
  } else { 
    fUpwindQuad_l[93] = hyb_3x3v_p1_surfx1_eval_quad_node_93_l(fc); 
  } 
  if ((-0.2371708245126284*(alphar[20]+alphar[13]+alphar[12]))+0.1767766952966368*alphar[6]-0.2371708245126284*alphar[5]+0.1767766952966368*(alphar[2]+alphar[1]+alphar[0]) > 0) { 
    fUpwindQuad_r[93] = hyb_3x3v_p1_surfx1_eval_quad_node_93_r(fc); 
  } else { 
    fUpwindQuad_r[93] = hyb_3x3v_p1_surfx1_eval_quad_node_93_l(fr); 
  } 
  if (0.1767766952966368*(alphal[6]+alphal[2]+alphal[1]+alphal[0]) > 0) { 
    fUpwindQuad_l[94] = hyb_3x3v_p1_surfx1_eval_quad_node_94_r(fl); 
  } else { 
    fUpwindQuad_l[94] = hyb_3x3v_p1_surfx1_eval_quad_node_94_l(fc); 
  } 
  if (0.1767766952966368*(alphar[6]+alphar[2]+alphar[1]+alphar[0]) > 0) { 
    fUpwindQuad_r[94] = hyb_3x3v_p1_surfx1_eval_quad_node_94_r(fc); 
  } else { 
    fUpwindQuad_r[94] = hyb_3x3v_p1_surfx1_eval_quad_node_94_l(fr); 
  } 
  if (0.2371708245126284*(alphal[20]+alphal[13]+alphal[12])+0.1767766952966368*alphal[6]+0.2371708245126284*alphal[5]+0.1767766952966368*(alphal[2]+alphal[1]+alphal[0]) > 0) { 
    fUpwindQuad_l[95] = hyb_3x3v_p1_surfx1_eval_quad_node_95_r(fl); 
  } else { 
    fUpwindQuad_l[95] = hyb_3x3v_p1_surfx1_eval_quad_node_95_l(fc); 
  } 
  if (0.2371708245126284*(alphar[20]+alphar[13]+alphar[12])+0.1767766952966368*alphar[6]+0.2371708245126284*alphar[5]+0.1767766952966368*(alphar[2]+alphar[1]+alphar[0]) > 0) { 
    fUpwindQuad_r[95] = hyb_3x3v_p1_surfx1_eval_quad_node_95_r(fc); 
  } else { 
    fUpwindQuad_r[95] = hyb_3x3v_p1_surfx1_eval_quad_node_95_l(fr); 
  } 
  if ((-0.2371708245126284*alphal[20])+0.2371708245126284*alphal[17]-0.2371708245126284*(alphal[13]+alphal[12])+0.2371708245126284*(alphal[10]+alphal[9])+0.1767766952966368*alphal[6]-0.2371708245126284*alphal[5]+0.2371708245126284*alphal[4]+0.1767766952966368*(alphal[2]+alphal[1]+alphal[0]) > 0) { 
    fUpwindQuad_l[96] = hyb_3x3v_p1_surfx1_eval_quad_node_96_r(fl); 
  } else { 
    fUpwindQuad_l[96] = hyb_3x3v_p1_surfx1_eval_quad_node_96_l(fc); 
  } 
  if ((-0.2371708245126284*alphar[20])+0.2371708245126284*alphar[17]-0.2371708245126284*(alphar[13]+alphar[12])+0.2371708245126284*(alphar[10]+alphar[9])+0.1767766952966368*alphar[6]-0.2371708245126284*alphar[5]+0.2371708245126284*alphar[4]+0.1767766952966368*(alphar[2]+alphar[1]+alphar[0]) > 0) { 
    fUpwindQuad_r[96] = hyb_3x3v_p1_surfx1_eval_quad_node_96_r(fc); 
  } else { 
    fUpwindQuad_r[96] = hyb_3x3v_p1_surfx1_eval_quad_node_96_l(fr); 
  } 
  if (0.2371708245126284*(alphal[17]+alphal[10]+alphal[9])+0.1767766952966368*alphal[6]+0.2371708245126284*alphal[4]+0.1767766952966368*(alphal[2]+alphal[1]+alphal[0]) > 0) { 
    fUpwindQuad_l[97] = hyb_3x3v_p1_surfx1_eval_quad_node_97_r(fl); 
  } else { 
    fUpwindQuad_l[97] = hyb_3x3v_p1_surfx1_eval_quad_node_97_l(fc); 
  } 
  if (0.2371708245126284*(alphar[17]+alphar[10]+alphar[9])+0.1767766952966368*alphar[6]+0.2371708245126284*alphar[4]+0.1767766952966368*(alphar[2]+alphar[1]+alphar[0]) > 0) { 
    fUpwindQuad_r[97] = hyb_3x3v_p1_surfx1_eval_quad_node_97_r(fc); 
  } else { 
    fUpwindQuad_r[97] = hyb_3x3v_p1_surfx1_eval_quad_node_97_l(fr); 
  } 
  if (0.2371708245126284*(alphal[20]+alphal[17]+alphal[13]+alphal[12]+alphal[10]+alphal[9])+0.1767766952966368*alphal[6]+0.2371708245126284*(alphal[5]+alphal[4])+0.1767766952966368*(alphal[2]+alphal[1]+alphal[0]) > 0) { 
    fUpwindQuad_l[98] = hyb_3x3v_p1_surfx1_eval_quad_node_98_r(fl); 
  } else { 
    fUpwindQuad_l[98] = hyb_3x3v_p1_surfx1_eval_quad_node_98_l(fc); 
  } 
  if (0.2371708245126284*(alphar[20]+alphar[17]+alphar[13]+alphar[12]+alphar[10]+alphar[9])+0.1767766952966368*alphar[6]+0.2371708245126284*(alphar[5]+alphar[4])+0.1767766952966368*(alphar[2]+alphar[1]+alphar[0]) > 0) { 
    fUpwindQuad_r[98] = hyb_3x3v_p1_surfx1_eval_quad_node_98_r(fc); 
  } else { 
    fUpwindQuad_r[98] = hyb_3x3v_p1_surfx1_eval_quad_node_98_l(fr); 
  } 
  if ((-0.2371708245126284*(alphal[20]+alphal[17]))+0.2371708245126284*alphal[16]-0.2371708245126284*(alphal[13]+alphal[12]+alphal[10]+alphal[9])+0.2371708245126284*(alphal[8]+alphal[7])+0.1767766952966368*alphal[6]-0.2371708245126284*(alphal[5]+alphal[4])+0.2371708245126284*alphal[3]+0.1767766952966368*(alphal[2]+alphal[1]+alphal[0]) > 0) { 
    fUpwindQuad_l[99] = hyb_3x3v_p1_surfx1_eval_quad_node_99_r(fl); 
  } else { 
    fUpwindQuad_l[99] = hyb_3x3v_p1_surfx1_eval_quad_node_99_l(fc); 
  } 
  if ((-0.2371708245126284*(alphar[20]+alphar[17]))+0.2371708245126284*alphar[16]-0.2371708245126284*(alphar[13]+alphar[12]+alphar[10]+alphar[9])+0.2371708245126284*(alphar[8]+alphar[7])+0.1767766952966368*alphar[6]-0.2371708245126284*(alphar[5]+alphar[4])+0.2371708245126284*alphar[3]+0.1767766952966368*(alphar[2]+alphar[1]+alphar[0]) > 0) { 
    fUpwindQuad_r[99] = hyb_3x3v_p1_surfx1_eval_quad_node_99_r(fc); 
  } else { 
    fUpwindQuad_r[99] = hyb_3x3v_p1_surfx1_eval_quad_node_99_l(fr); 
  } 
  if ((-0.2371708245126284*alphal[17])+0.2371708245126284*alphal[16]-0.2371708245126284*(alphal[10]+alphal[9])+0.2371708245126284*(alphal[8]+alphal[7])+0.1767766952966368*alphal[6]-0.2371708245126284*alphal[4]+0.2371708245126284*alphal[3]+0.1767766952966368*(alphal[2]+alphal[1]+alphal[0]) > 0) { 
    fUpwindQuad_l[100] = hyb_3x3v_p1_surfx1_eval_quad_node_100_r(fl); 
  } else { 
    fUpwindQuad_l[100] = hyb_3x3v_p1_surfx1_eval_quad_node_100_l(fc); 
  } 
  if ((-0.2371708245126284*alphar[17])+0.2371708245126284*alphar[16]-0.2371708245126284*(alphar[10]+alphar[9])+0.2371708245126284*(alphar[8]+alphar[7])+0.1767766952966368*alphar[6]-0.2371708245126284*alphar[4]+0.2371708245126284*alphar[3]+0.1767766952966368*(alphar[2]+alphar[1]+alphar[0]) > 0) { 
    fUpwindQuad_r[100] = hyb_3x3v_p1_surfx1_eval_quad_node_100_r(fc); 
  } else { 
    fUpwindQuad_r[100] = hyb_3x3v_p1_surfx1_eval_quad_node_100_l(fr); 
  } 
  if (0.2371708245126284*alphal[20]-0.2371708245126284*alphal[17]+0.2371708245126284*(alphal[16]+alphal[13]+alphal[12])-0.2371708245126284*(alphal[10]+alphal[9])+0.2371708245126284*(alphal[8]+alphal[7])+0.1767766952966368*alphal[6]+0.2371708245126284*alphal[5]-0.2371708245126284*alphal[4]+0.2371708245126284*alphal[3]+0.1767766952966368*(alphal[2]+alphal[1]+alphal[0]) > 0) { 
    fUpwindQuad_l[101] = hyb_3x3v_p1_surfx1_eval_quad_node_101_r(fl); 
  } else { 
    fUpwindQuad_l[101] = hyb_3x3v_p1_surfx1_eval_quad_node_101_l(fc); 
  } 
  if (0.2371708245126284*alphar[20]-0.2371708245126284*alphar[17]+0.2371708245126284*(alphar[16]+alphar[13]+alphar[12])-0.2371708245126284*(alphar[10]+alphar[9])+0.2371708245126284*(alphar[8]+alphar[7])+0.1767766952966368*alphar[6]+0.2371708245126284*alphar[5]-0.2371708245126284*alphar[4]+0.2371708245126284*alphar[3]+0.1767766952966368*(alphar[2]+alphar[1]+alphar[0]) > 0) { 
    fUpwindQuad_r[101] = hyb_3x3v_p1_surfx1_eval_quad_node_101_r(fc); 
  } else { 
    fUpwindQuad_r[101] = hyb_3x3v_p1_surfx1_eval_quad_node_101_l(fr); 
  } 
  if ((-0.2371708245126284*alphal[20])+0.2371708245126284*alphal[16]-0.2371708245126284*(alphal[13]+alphal[12])+0.2371708245126284*(alphal[8]+alphal[7])+0.1767766952966368*alphal[6]-0.2371708245126284*alphal[5]+0.2371708245126284*alphal[3]+0.1767766952966368*(alphal[2]+alphal[1]+alphal[0]) > 0) { 
    fUpwindQuad_l[102] = hyb_3x3v_p1_surfx1_eval_quad_node_102_r(fl); 
  } else { 
    fUpwindQuad_l[102] = hyb_3x3v_p1_surfx1_eval_quad_node_102_l(fc); 
  } 
  if ((-0.2371708245126284*alphar[20])+0.2371708245126284*alphar[16]-0.2371708245126284*(alphar[13]+alphar[12])+0.2371708245126284*(alphar[8]+alphar[7])+0.1767766952966368*alphar[6]-0.2371708245126284*alphar[5]+0.2371708245126284*alphar[3]+0.1767766952966368*(alphar[2]+alphar[1]+alphar[0]) > 0) { 
    fUpwindQuad_r[102] = hyb_3x3v_p1_surfx1_eval_quad_node_102_r(fc); 
  } else { 
    fUpwindQuad_r[102] = hyb_3x3v_p1_surfx1_eval_quad_node_102_l(fr); 
  } 
  if (0.2371708245126284*(alphal[16]+alphal[8]+alphal[7])+0.1767766952966368*alphal[6]+0.2371708245126284*alphal[3]+0.1767766952966368*(alphal[2]+alphal[1]+alphal[0]) > 0) { 
    fUpwindQuad_l[103] = hyb_3x3v_p1_surfx1_eval_quad_node_103_r(fl); 
  } else { 
    fUpwindQuad_l[103] = hyb_3x3v_p1_surfx1_eval_quad_node_103_l(fc); 
  } 
  if (0.2371708245126284*(alphar[16]+alphar[8]+alphar[7])+0.1767766952966368*alphar[6]+0.2371708245126284*alphar[3]+0.1767766952966368*(alphar[2]+alphar[1]+alphar[0]) > 0) { 
    fUpwindQuad_r[103] = hyb_3x3v_p1_surfx1_eval_quad_node_103_r(fc); 
  } else { 
    fUpwindQuad_r[103] = hyb_3x3v_p1_surfx1_eval_quad_node_103_l(fr); 
  } 
  if (0.2371708245126284*(alphal[20]+alphal[16]+alphal[13]+alphal[12]+alphal[8]+alphal[7])+0.1767766952966368*alphal[6]+0.2371708245126284*(alphal[5]+alphal[3])+0.1767766952966368*(alphal[2]+alphal[1]+alphal[0]) > 0) { 
    fUpwindQuad_l[104] = hyb_3x3v_p1_surfx1_eval_quad_node_104_r(fl); 
  } else { 
    fUpwindQuad_l[104] = hyb_3x3v_p1_surfx1_eval_quad_node_104_l(fc); 
  } 
  if (0.2371708245126284*(alphar[20]+alphar[16]+alphar[13]+alphar[12]+alphar[8]+alphar[7])+0.1767766952966368*alphar[6]+0.2371708245126284*(alphar[5]+alphar[3])+0.1767766952966368*(alphar[2]+alphar[1]+alphar[0]) > 0) { 
    fUpwindQuad_r[104] = hyb_3x3v_p1_surfx1_eval_quad_node_104_r(fc); 
  } else { 
    fUpwindQuad_r[104] = hyb_3x3v_p1_surfx1_eval_quad_node_104_l(fr); 
  } 
  if ((-0.2371708245126284*alphal[20])+0.2371708245126284*(alphal[17]+alphal[16])-0.2371708245126284*(alphal[13]+alphal[12])+0.2371708245126284*(alphal[10]+alphal[9]+alphal[8]+alphal[7])+0.1767766952966368*alphal[6]-0.2371708245126284*alphal[5]+0.2371708245126284*(alphal[4]+alphal[3])+0.1767766952966368*(alphal[2]+alphal[1]+alphal[0]) > 0) { 
    fUpwindQuad_l[105] = hyb_3x3v_p1_surfx1_eval_quad_node_105_r(fl); 
  } else { 
    fUpwindQuad_l[105] = hyb_3x3v_p1_surfx1_eval_quad_node_105_l(fc); 
  } 
  if ((-0.2371708245126284*alphar[20])+0.2371708245126284*(alphar[17]+alphar[16])-0.2371708245126284*(alphar[13]+alphar[12])+0.2371708245126284*(alphar[10]+alphar[9]+alphar[8]+alphar[7])+0.1767766952966368*alphar[6]-0.2371708245126284*alphar[5]+0.2371708245126284*(alphar[4]+alphar[3])+0.1767766952966368*(alphar[2]+alphar[1]+alphar[0]) > 0) { 
    fUpwindQuad_r[105] = hyb_3x3v_p1_surfx1_eval_quad_node_105_r(fc); 
  } else { 
    fUpwindQuad_r[105] = hyb_3x3v_p1_surfx1_eval_quad_node_105_l(fr); 
  } 
  if (0.2371708245126284*(alphal[17]+alphal[16]+alphal[10]+alphal[9]+alphal[8]+alphal[7])+0.1767766952966368*alphal[6]+0.2371708245126284*(alphal[4]+alphal[3])+0.1767766952966368*(alphal[2]+alphal[1]+alphal[0]) > 0) { 
    fUpwindQuad_l[106] = hyb_3x3v_p1_surfx1_eval_quad_node_106_r(fl); 
  } else { 
    fUpwindQuad_l[106] = hyb_3x3v_p1_surfx1_eval_quad_node_106_l(fc); 
  } 
  if (0.2371708245126284*(alphar[17]+alphar[16]+alphar[10]+alphar[9]+alphar[8]+alphar[7])+0.1767766952966368*alphar[6]+0.2371708245126284*(alphar[4]+alphar[3])+0.1767766952966368*(alphar[2]+alphar[1]+alphar[0]) > 0) { 
    fUpwindQuad_r[106] = hyb_3x3v_p1_surfx1_eval_quad_node_106_r(fc); 
  } else { 
    fUpwindQuad_r[106] = hyb_3x3v_p1_surfx1_eval_quad_node_106_l(fr); 
  } 
  if (0.2371708245126284*(alphal[20]+alphal[17]+alphal[16]+alphal[13]+alphal[12]+alphal[10]+alphal[9]+alphal[8]+alphal[7])+0.1767766952966368*alphal[6]+0.2371708245126284*(alphal[5]+alphal[4]+alphal[3])+0.1767766952966368*(alphal[2]+alphal[1]+alphal[0]) > 0) { 
    fUpwindQuad_l[107] = hyb_3x3v_p1_surfx1_eval_quad_node_107_r(fl); 
  } else { 
    fUpwindQuad_l[107] = hyb_3x3v_p1_surfx1_eval_quad_node_107_l(fc); 
  } 
  if (0.2371708245126284*(alphar[20]+alphar[17]+alphar[16]+alphar[13]+alphar[12]+alphar[10]+alphar[9]+alphar[8]+alphar[7])+0.1767766952966368*alphar[6]+0.2371708245126284*(alphar[5]+alphar[4]+alphar[3])+0.1767766952966368*(alphar[2]+alphar[1]+alphar[0]) > 0) { 
    fUpwindQuad_r[107] = hyb_3x3v_p1_surfx1_eval_quad_node_107_r(fc); 
  } else { 
    fUpwindQuad_r[107] = hyb_3x3v_p1_surfx1_eval_quad_node_107_l(fr); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_3x3v_p1_xdir_upwind_quad_to_modal(fUpwindQuad_l, fUpwind_l); 
  hyb_3x3v_p1_xdir_upwind_quad_to_modal(fUpwindQuad_r, fUpwind_r); 

  Ghat_l[0] = 0.1767766952966368*alphal[20]*fUpwind_l[20]+0.1767766952966368*alphal[17]*fUpwind_l[17]+0.1767766952966368*alphal[16]*fUpwind_l[16]+0.1767766952966368*alphal[13]*fUpwind_l[13]+0.1767766952966368*alphal[12]*fUpwind_l[12]+0.1767766952966368*alphal[10]*fUpwind_l[10]+0.1767766952966368*alphal[9]*fUpwind_l[9]+0.1767766952966368*alphal[8]*fUpwind_l[8]+0.1767766952966368*alphal[7]*fUpwind_l[7]+0.1767766952966368*alphal[6]*fUpwind_l[6]+0.1767766952966368*alphal[5]*fUpwind_l[5]+0.1767766952966368*alphal[4]*fUpwind_l[4]+0.1767766952966368*alphal[3]*fUpwind_l[3]+0.1767766952966368*alphal[2]*fUpwind_l[2]+0.1767766952966368*alphal[1]*fUpwind_l[1]+0.1767766952966368*alphal[0]*fUpwind_l[0]; 
  Ghat_l[1] = 0.1767766952966368*alphal[13]*fUpwind_l[20]+0.1767766952966368*fUpwind_l[13]*alphal[20]+0.1767766952966368*alphal[10]*fUpwind_l[17]+0.1767766952966368*fUpwind_l[10]*alphal[17]+0.1767766952966368*alphal[8]*fUpwind_l[16]+0.1767766952966368*fUpwind_l[8]*alphal[16]+0.1767766952966368*alphal[5]*fUpwind_l[12]+0.1767766952966368*fUpwind_l[5]*alphal[12]+0.1767766952966368*alphal[4]*fUpwind_l[9]+0.1767766952966368*fUpwind_l[4]*alphal[9]+0.1767766952966368*alphal[3]*fUpwind_l[7]+0.1767766952966368*fUpwind_l[3]*alphal[7]+0.1767766952966368*alphal[2]*fUpwind_l[6]+0.1767766952966368*fUpwind_l[2]*alphal[6]+0.1767766952966368*alphal[0]*fUpwind_l[1]+0.1767766952966368*fUpwind_l[0]*alphal[1]; 
  Ghat_l[2] = 0.1767766952966368*alphal[12]*fUpwind_l[20]+0.1767766952966368*fUpwind_l[12]*alphal[20]+0.1767766952966368*alphal[9]*fUpwind_l[17]+0.1767766952966368*fUpwind_l[9]*alphal[17]+0.1767766952966368*alphal[7]*fUpwind_l[16]+0.1767766952966368*fUpwind_l[7]*alphal[16]+0.1767766952966368*alphal[5]*fUpwind_l[13]+0.1767766952966368*fUpwind_l[5]*alphal[13]+0.1767766952966368*alphal[4]*fUpwind_l[10]+0.1767766952966368*fUpwind_l[4]*alphal[10]+0.1767766952966368*alphal[3]*fUpwind_l[8]+0.1767766952966368*fUpwind_l[3]*alphal[8]+0.1767766952966368*alphal[1]*fUpwind_l[6]+0.1767766952966368*fUpwind_l[1]*alphal[6]+0.1767766952966368*alphal[0]*fUpwind_l[2]+0.1767766952966368*fUpwind_l[0]*alphal[2]; 
  Ghat_l[3] = 0.1581138830084189*alphal[16]*fUpwind_l[37]+0.1581138830084189*alphal[8]*fUpwind_l[34]+0.1581138830084189*alphal[7]*fUpwind_l[33]+0.1581138830084189*alphal[3]*fUpwind_l[32]+0.1767766952966368*alphal[20]*fUpwind_l[27]+0.1767766952966368*alphal[17]*fUpwind_l[26]+0.1767766952966368*alphal[13]*fUpwind_l[22]+0.1767766952966368*alphal[12]*fUpwind_l[21]+0.1767766952966368*alphal[10]*fUpwind_l[19]+0.1767766952966368*alphal[9]*fUpwind_l[18]+0.1767766952966368*alphal[6]*fUpwind_l[16]+0.1767766952966368*fUpwind_l[6]*alphal[16]+0.1767766952966368*alphal[5]*fUpwind_l[14]+0.1767766952966368*alphal[4]*fUpwind_l[11]+0.1767766952966368*alphal[2]*fUpwind_l[8]+0.1767766952966368*fUpwind_l[2]*alphal[8]+0.1767766952966368*alphal[1]*fUpwind_l[7]+0.1767766952966368*fUpwind_l[1]*alphal[7]+0.1767766952966368*alphal[0]*fUpwind_l[3]+0.1767766952966368*fUpwind_l[0]*alphal[3]; 
  Ghat_l[4] = 0.1581138830084189*alphal[17]*fUpwind_l[53]+0.1581138830084189*alphal[10]*fUpwind_l[50]+0.1581138830084189*alphal[9]*fUpwind_l[49]+0.1581138830084189*alphal[4]*fUpwind_l[48]+0.1767766952966368*alphal[20]*fUpwind_l[28]+0.1767766952966368*alphal[16]*fUpwind_l[26]+0.1767766952966368*alphal[13]*fUpwind_l[24]+0.1767766952966368*alphal[12]*fUpwind_l[23]+0.1767766952966368*alphal[8]*fUpwind_l[19]+0.1767766952966368*alphal[7]*fUpwind_l[18]+0.1767766952966368*alphal[6]*fUpwind_l[17]+0.1767766952966368*fUpwind_l[6]*alphal[17]+0.1767766952966368*alphal[5]*fUpwind_l[15]+0.1767766952966368*alphal[3]*fUpwind_l[11]+0.1767766952966368*alphal[2]*fUpwind_l[10]+0.1767766952966368*fUpwind_l[2]*alphal[10]+0.1767766952966368*alphal[1]*fUpwind_l[9]+0.1767766952966368*fUpwind_l[1]*alphal[9]+0.1767766952966368*alphal[0]*fUpwind_l[4]+0.1767766952966368*fUpwind_l[0]*alphal[4]; 
  Ghat_l[5] = 0.1581138830084189*alphal[20]*fUpwind_l[69]+0.1581138830084189*alphal[13]*fUpwind_l[66]+0.1581138830084189*alphal[12]*fUpwind_l[65]+0.1581138830084189*alphal[5]*fUpwind_l[64]+0.1767766952966368*alphal[17]*fUpwind_l[28]+0.1767766952966368*alphal[16]*fUpwind_l[27]+0.1767766952966368*alphal[10]*fUpwind_l[24]+0.1767766952966368*alphal[9]*fUpwind_l[23]+0.1767766952966368*alphal[8]*fUpwind_l[22]+0.1767766952966368*alphal[7]*fUpwind_l[21]+0.1767766952966368*alphal[6]*fUpwind_l[20]+0.1767766952966368*fUpwind_l[6]*alphal[20]+0.1767766952966368*alphal[4]*fUpwind_l[15]+0.1767766952966368*alphal[3]*fUpwind_l[14]+0.1767766952966368*alphal[2]*fUpwind_l[13]+0.1767766952966368*fUpwind_l[2]*alphal[13]+0.1767766952966368*alphal[1]*fUpwind_l[12]+0.1767766952966368*fUpwind_l[1]*alphal[12]+0.1767766952966368*alphal[0]*fUpwind_l[5]+0.1767766952966368*fUpwind_l[0]*alphal[5]; 
  Ghat_l[6] = 0.1767766952966368*alphal[5]*fUpwind_l[20]+0.1767766952966368*fUpwind_l[5]*alphal[20]+0.1767766952966368*alphal[4]*fUpwind_l[17]+0.1767766952966368*fUpwind_l[4]*alphal[17]+0.1767766952966368*alphal[3]*fUpwind_l[16]+0.1767766952966368*fUpwind_l[3]*alphal[16]+0.1767766952966368*alphal[12]*fUpwind_l[13]+0.1767766952966368*fUpwind_l[12]*alphal[13]+0.1767766952966368*alphal[9]*fUpwind_l[10]+0.1767766952966368*fUpwind_l[9]*alphal[10]+0.1767766952966368*alphal[7]*fUpwind_l[8]+0.1767766952966368*fUpwind_l[7]*alphal[8]+0.1767766952966368*alphal[0]*fUpwind_l[6]+0.1767766952966368*fUpwind_l[0]*alphal[6]+0.1767766952966368*alphal[1]*fUpwind_l[2]+0.1767766952966368*fUpwind_l[1]*alphal[2]; 
  Ghat_l[7] = 0.1581138830084189*alphal[8]*fUpwind_l[37]+0.1581138830084189*alphal[16]*fUpwind_l[34]+0.1581138830084189*alphal[3]*fUpwind_l[33]+0.1581138830084189*alphal[7]*fUpwind_l[32]+0.1767766952966368*alphal[13]*fUpwind_l[27]+0.1767766952966368*alphal[10]*fUpwind_l[26]+0.1767766952966368*alphal[20]*fUpwind_l[22]+0.1767766952966368*alphal[5]*fUpwind_l[21]+0.1767766952966368*alphal[17]*fUpwind_l[19]+0.1767766952966368*alphal[4]*fUpwind_l[18]+0.1767766952966368*alphal[2]*fUpwind_l[16]+0.1767766952966368*fUpwind_l[2]*alphal[16]+0.1767766952966368*alphal[12]*fUpwind_l[14]+0.1767766952966368*alphal[9]*fUpwind_l[11]+0.1767766952966368*alphal[6]*fUpwind_l[8]+0.1767766952966368*fUpwind_l[6]*alphal[8]+0.1767766952966368*alphal[0]*fUpwind_l[7]+0.1767766952966368*fUpwind_l[0]*alphal[7]+0.1767766952966368*alphal[1]*fUpwind_l[3]+0.1767766952966368*fUpwind_l[1]*alphal[3]; 
  Ghat_l[8] = 0.1581138830084189*alphal[7]*fUpwind_l[37]+0.1581138830084189*alphal[3]*fUpwind_l[34]+0.1581138830084189*alphal[16]*fUpwind_l[33]+0.1581138830084189*alphal[8]*fUpwind_l[32]+0.1767766952966368*alphal[12]*fUpwind_l[27]+0.1767766952966368*alphal[9]*fUpwind_l[26]+0.1767766952966368*alphal[5]*fUpwind_l[22]+0.1767766952966368*alphal[20]*fUpwind_l[21]+0.1767766952966368*alphal[4]*fUpwind_l[19]+0.1767766952966368*alphal[17]*fUpwind_l[18]+0.1767766952966368*alphal[1]*fUpwind_l[16]+0.1767766952966368*fUpwind_l[1]*alphal[16]+0.1767766952966368*alphal[13]*fUpwind_l[14]+0.1767766952966368*alphal[10]*fUpwind_l[11]+0.1767766952966368*alphal[0]*fUpwind_l[8]+0.1767766952966368*fUpwind_l[0]*alphal[8]+0.1767766952966368*alphal[6]*fUpwind_l[7]+0.1767766952966368*fUpwind_l[6]*alphal[7]+0.1767766952966368*alphal[2]*fUpwind_l[3]+0.1767766952966368*fUpwind_l[2]*alphal[3]; 
  Ghat_l[9] = 0.1581138830084189*alphal[10]*fUpwind_l[53]+0.1581138830084189*alphal[17]*fUpwind_l[50]+0.1581138830084189*alphal[4]*fUpwind_l[49]+0.1581138830084189*alphal[9]*fUpwind_l[48]+0.1767766952966368*alphal[13]*fUpwind_l[28]+0.1767766952966368*alphal[8]*fUpwind_l[26]+0.1767766952966368*alphal[20]*fUpwind_l[24]+0.1767766952966368*alphal[5]*fUpwind_l[23]+0.1767766952966368*alphal[16]*fUpwind_l[19]+0.1767766952966368*alphal[3]*fUpwind_l[18]+0.1767766952966368*alphal[2]*fUpwind_l[17]+0.1767766952966368*fUpwind_l[2]*alphal[17]+0.1767766952966368*alphal[12]*fUpwind_l[15]+0.1767766952966368*alphal[7]*fUpwind_l[11]+0.1767766952966368*alphal[6]*fUpwind_l[10]+0.1767766952966368*fUpwind_l[6]*alphal[10]+0.1767766952966368*alphal[0]*fUpwind_l[9]+0.1767766952966368*fUpwind_l[0]*alphal[9]+0.1767766952966368*alphal[1]*fUpwind_l[4]+0.1767766952966368*fUpwind_l[1]*alphal[4]; 
  Ghat_l[10] = 0.1581138830084189*alphal[9]*fUpwind_l[53]+0.1581138830084189*alphal[4]*fUpwind_l[50]+0.1581138830084189*alphal[17]*fUpwind_l[49]+0.1581138830084189*alphal[10]*fUpwind_l[48]+0.1767766952966368*alphal[12]*fUpwind_l[28]+0.1767766952966368*alphal[7]*fUpwind_l[26]+0.1767766952966368*alphal[5]*fUpwind_l[24]+0.1767766952966368*alphal[20]*fUpwind_l[23]+0.1767766952966368*alphal[3]*fUpwind_l[19]+0.1767766952966368*alphal[16]*fUpwind_l[18]+0.1767766952966368*alphal[1]*fUpwind_l[17]+0.1767766952966368*fUpwind_l[1]*alphal[17]+0.1767766952966368*alphal[13]*fUpwind_l[15]+0.1767766952966368*alphal[8]*fUpwind_l[11]+0.1767766952966368*alphal[0]*fUpwind_l[10]+0.1767766952966368*fUpwind_l[0]*alphal[10]+0.1767766952966368*alphal[6]*fUpwind_l[9]+0.1767766952966368*fUpwind_l[6]*alphal[9]+0.1767766952966368*alphal[2]*fUpwind_l[4]+0.1767766952966368*fUpwind_l[2]*alphal[4]; 
  Ghat_l[11] = 0.1581138830084189*alphal[17]*fUpwind_l[59]+0.1581138830084189*alphal[10]*fUpwind_l[55]+0.1581138830084189*alphal[9]*fUpwind_l[54]+0.1581138830084189*alphal[4]*fUpwind_l[51]+0.1581138830084189*alphal[16]*fUpwind_l[43]+0.1581138830084189*alphal[8]*fUpwind_l[39]+0.1581138830084189*alphal[7]*fUpwind_l[38]+0.1581138830084189*alphal[3]*fUpwind_l[35]+0.1767766952966368*alphal[20]*fUpwind_l[31]+0.1767766952966368*alphal[13]*fUpwind_l[30]+0.1767766952966368*alphal[12]*fUpwind_l[29]+0.1767766952966368*alphal[6]*fUpwind_l[26]+0.1767766952966368*alphal[5]*fUpwind_l[25]+0.1767766952966368*alphal[2]*fUpwind_l[19]+0.1767766952966368*alphal[1]*fUpwind_l[18]+0.1767766952966368*alphal[16]*fUpwind_l[17]+0.1767766952966368*fUpwind_l[16]*alphal[17]+0.1767766952966368*alphal[0]*fUpwind_l[11]+0.1767766952966368*alphal[8]*fUpwind_l[10]+0.1767766952966368*fUpwind_l[8]*alphal[10]+0.1767766952966368*alphal[7]*fUpwind_l[9]+0.1767766952966368*fUpwind_l[7]*alphal[9]+0.1767766952966368*alphal[3]*fUpwind_l[4]+0.1767766952966368*fUpwind_l[3]*alphal[4]; 
  Ghat_l[12] = 0.1581138830084189*alphal[13]*fUpwind_l[69]+0.1581138830084189*alphal[20]*fUpwind_l[66]+0.1581138830084189*alphal[5]*fUpwind_l[65]+0.1581138830084189*alphal[12]*fUpwind_l[64]+0.1767766952966368*alphal[10]*fUpwind_l[28]+0.1767766952966368*alphal[8]*fUpwind_l[27]+0.1767766952966368*alphal[17]*fUpwind_l[24]+0.1767766952966368*alphal[4]*fUpwind_l[23]+0.1767766952966368*alphal[16]*fUpwind_l[22]+0.1767766952966368*alphal[3]*fUpwind_l[21]+0.1767766952966368*alphal[2]*fUpwind_l[20]+0.1767766952966368*fUpwind_l[2]*alphal[20]+0.1767766952966368*alphal[9]*fUpwind_l[15]+0.1767766952966368*alphal[7]*fUpwind_l[14]+0.1767766952966368*alphal[6]*fUpwind_l[13]+0.1767766952966368*fUpwind_l[6]*alphal[13]+0.1767766952966368*alphal[0]*fUpwind_l[12]+0.1767766952966368*fUpwind_l[0]*alphal[12]+0.1767766952966368*alphal[1]*fUpwind_l[5]+0.1767766952966368*fUpwind_l[1]*alphal[5]; 
  Ghat_l[13] = 0.1581138830084189*alphal[12]*fUpwind_l[69]+0.1581138830084189*alphal[5]*fUpwind_l[66]+0.1581138830084189*alphal[20]*fUpwind_l[65]+0.1581138830084189*alphal[13]*fUpwind_l[64]+0.1767766952966368*alphal[9]*fUpwind_l[28]+0.1767766952966368*alphal[7]*fUpwind_l[27]+0.1767766952966368*alphal[4]*fUpwind_l[24]+0.1767766952966368*alphal[17]*fUpwind_l[23]+0.1767766952966368*alphal[3]*fUpwind_l[22]+0.1767766952966368*alphal[16]*fUpwind_l[21]+0.1767766952966368*alphal[1]*fUpwind_l[20]+0.1767766952966368*fUpwind_l[1]*alphal[20]+0.1767766952966368*alphal[10]*fUpwind_l[15]+0.1767766952966368*alphal[8]*fUpwind_l[14]+0.1767766952966368*alphal[0]*fUpwind_l[13]+0.1767766952966368*fUpwind_l[0]*alphal[13]+0.1767766952966368*alphal[6]*fUpwind_l[12]+0.1767766952966368*fUpwind_l[6]*alphal[12]+0.1767766952966368*alphal[2]*fUpwind_l[5]+0.1767766952966368*fUpwind_l[2]*alphal[5]; 
  Ghat_l[14] = 0.1581138830084189*alphal[20]*fUpwind_l[75]+0.1581138830084189*alphal[13]*fUpwind_l[71]+0.1581138830084189*alphal[12]*fUpwind_l[70]+0.1581138830084189*alphal[5]*fUpwind_l[67]+0.1581138830084189*alphal[16]*fUpwind_l[44]+0.1581138830084189*alphal[8]*fUpwind_l[41]+0.1581138830084189*alphal[7]*fUpwind_l[40]+0.1581138830084189*alphal[3]*fUpwind_l[36]+0.1767766952966368*alphal[17]*fUpwind_l[31]+0.1767766952966368*alphal[10]*fUpwind_l[30]+0.1767766952966368*alphal[9]*fUpwind_l[29]+0.1767766952966368*alphal[6]*fUpwind_l[27]+0.1767766952966368*alphal[4]*fUpwind_l[25]+0.1767766952966368*alphal[2]*fUpwind_l[22]+0.1767766952966368*alphal[1]*fUpwind_l[21]+0.1767766952966368*alphal[16]*fUpwind_l[20]+0.1767766952966368*fUpwind_l[16]*alphal[20]+0.1767766952966368*alphal[0]*fUpwind_l[14]+0.1767766952966368*alphal[8]*fUpwind_l[13]+0.1767766952966368*fUpwind_l[8]*alphal[13]+0.1767766952966368*alphal[7]*fUpwind_l[12]+0.1767766952966368*fUpwind_l[7]*alphal[12]+0.1767766952966368*alphal[3]*fUpwind_l[5]+0.1767766952966368*fUpwind_l[3]*alphal[5]; 
  Ghat_l[15] = 0.1581138830084189*alphal[20]*fUpwind_l[76]+0.1581138830084189*alphal[13]*fUpwind_l[73]+0.1581138830084189*alphal[12]*fUpwind_l[72]+0.1581138830084189*alphal[5]*fUpwind_l[68]+0.1581138830084189*alphal[17]*fUpwind_l[60]+0.1581138830084189*alphal[10]*fUpwind_l[57]+0.1581138830084189*alphal[9]*fUpwind_l[56]+0.1581138830084189*alphal[4]*fUpwind_l[52]+0.1767766952966368*alphal[16]*fUpwind_l[31]+0.1767766952966368*alphal[8]*fUpwind_l[30]+0.1767766952966368*alphal[7]*fUpwind_l[29]+0.1767766952966368*alphal[6]*fUpwind_l[28]+0.1767766952966368*alphal[3]*fUpwind_l[25]+0.1767766952966368*alphal[2]*fUpwind_l[24]+0.1767766952966368*alphal[1]*fUpwind_l[23]+0.1767766952966368*alphal[17]*fUpwind_l[20]+0.1767766952966368*fUpwind_l[17]*alphal[20]+0.1767766952966368*alphal[0]*fUpwind_l[15]+0.1767766952966368*alphal[10]*fUpwind_l[13]+0.1767766952966368*fUpwind_l[10]*alphal[13]+0.1767766952966368*alphal[9]*fUpwind_l[12]+0.1767766952966368*fUpwind_l[9]*alphal[12]+0.1767766952966368*alphal[4]*fUpwind_l[5]+0.1767766952966368*fUpwind_l[4]*alphal[5]; 
  Ghat_l[16] = 0.1581138830084189*alphal[3]*fUpwind_l[37]+0.1581138830084189*alphal[7]*fUpwind_l[34]+0.1581138830084189*alphal[8]*fUpwind_l[33]+0.1581138830084189*alphal[16]*fUpwind_l[32]+0.1767766952966368*alphal[5]*fUpwind_l[27]+0.1767766952966368*alphal[4]*fUpwind_l[26]+0.1767766952966368*alphal[12]*fUpwind_l[22]+0.1767766952966368*alphal[13]*fUpwind_l[21]+0.1767766952966368*fUpwind_l[14]*alphal[20]+0.1767766952966368*alphal[9]*fUpwind_l[19]+0.1767766952966368*alphal[10]*fUpwind_l[18]+0.1767766952966368*fUpwind_l[11]*alphal[17]+0.1767766952966368*alphal[0]*fUpwind_l[16]+0.1767766952966368*fUpwind_l[0]*alphal[16]+0.1767766952966368*alphal[1]*fUpwind_l[8]+0.1767766952966368*fUpwind_l[1]*alphal[8]+0.1767766952966368*alphal[2]*fUpwind_l[7]+0.1767766952966368*fUpwind_l[2]*alphal[7]+0.1767766952966368*alphal[3]*fUpwind_l[6]+0.1767766952966368*fUpwind_l[3]*alphal[6]; 
  Ghat_l[17] = 0.1581138830084189*alphal[4]*fUpwind_l[53]+0.1581138830084189*alphal[9]*fUpwind_l[50]+0.1581138830084189*alphal[10]*fUpwind_l[49]+0.1581138830084189*alphal[17]*fUpwind_l[48]+0.1767766952966368*alphal[5]*fUpwind_l[28]+0.1767766952966368*alphal[3]*fUpwind_l[26]+0.1767766952966368*alphal[12]*fUpwind_l[24]+0.1767766952966368*alphal[13]*fUpwind_l[23]+0.1767766952966368*fUpwind_l[15]*alphal[20]+0.1767766952966368*alphal[7]*fUpwind_l[19]+0.1767766952966368*alphal[8]*fUpwind_l[18]+0.1767766952966368*alphal[0]*fUpwind_l[17]+0.1767766952966368*fUpwind_l[0]*alphal[17]+0.1767766952966368*fUpwind_l[11]*alphal[16]+0.1767766952966368*alphal[1]*fUpwind_l[10]+0.1767766952966368*fUpwind_l[1]*alphal[10]+0.1767766952966368*alphal[2]*fUpwind_l[9]+0.1767766952966368*fUpwind_l[2]*alphal[9]+0.1767766952966368*alphal[4]*fUpwind_l[6]+0.1767766952966368*fUpwind_l[4]*alphal[6]; 
  Ghat_l[18] = 0.1581138830084189*alphal[10]*fUpwind_l[59]+0.1581138830084189*alphal[17]*fUpwind_l[55]+0.1581138830084189*alphal[4]*fUpwind_l[54]+0.1581138830084189*alphal[9]*fUpwind_l[51]+0.1581138830084189*alphal[8]*fUpwind_l[43]+0.1581138830084189*alphal[16]*fUpwind_l[39]+0.1581138830084189*alphal[3]*fUpwind_l[38]+0.1581138830084189*alphal[7]*fUpwind_l[35]+0.1767766952966368*alphal[13]*fUpwind_l[31]+0.1767766952966368*alphal[20]*fUpwind_l[30]+0.1767766952966368*alphal[5]*fUpwind_l[29]+0.1767766952966368*alphal[2]*fUpwind_l[26]+0.1767766952966368*alphal[12]*fUpwind_l[25]+0.1767766952966368*alphal[6]*fUpwind_l[19]+0.1767766952966368*alphal[0]*fUpwind_l[18]+0.1767766952966368*alphal[8]*fUpwind_l[17]+0.1767766952966368*fUpwind_l[8]*alphal[17]+0.1767766952966368*alphal[10]*fUpwind_l[16]+0.1767766952966368*fUpwind_l[10]*alphal[16]+0.1767766952966368*alphal[1]*fUpwind_l[11]+0.1767766952966368*alphal[3]*fUpwind_l[9]+0.1767766952966368*fUpwind_l[3]*alphal[9]+0.1767766952966368*alphal[4]*fUpwind_l[7]+0.1767766952966368*fUpwind_l[4]*alphal[7]; 
  Ghat_l[19] = 0.1581138830084189*alphal[9]*fUpwind_l[59]+0.1581138830084189*alphal[4]*fUpwind_l[55]+0.1581138830084189*alphal[17]*fUpwind_l[54]+0.1581138830084189*alphal[10]*fUpwind_l[51]+0.1581138830084189*alphal[7]*fUpwind_l[43]+0.1581138830084189*alphal[3]*fUpwind_l[39]+0.1581138830084189*alphal[16]*fUpwind_l[38]+0.1581138830084189*alphal[8]*fUpwind_l[35]+0.1767766952966368*alphal[12]*fUpwind_l[31]+0.1767766952966368*alphal[5]*fUpwind_l[30]+0.1767766952966368*alphal[20]*fUpwind_l[29]+0.1767766952966368*alphal[1]*fUpwind_l[26]+0.1767766952966368*alphal[13]*fUpwind_l[25]+0.1767766952966368*alphal[0]*fUpwind_l[19]+0.1767766952966368*alphal[6]*fUpwind_l[18]+0.1767766952966368*alphal[7]*fUpwind_l[17]+0.1767766952966368*fUpwind_l[7]*alphal[17]+0.1767766952966368*alphal[9]*fUpwind_l[16]+0.1767766952966368*fUpwind_l[9]*alphal[16]+0.1767766952966368*alphal[2]*fUpwind_l[11]+0.1767766952966368*alphal[3]*fUpwind_l[10]+0.1767766952966368*fUpwind_l[3]*alphal[10]+0.1767766952966368*alphal[4]*fUpwind_l[8]+0.1767766952966368*fUpwind_l[4]*alphal[8]; 
  Ghat_l[20] = 0.1581138830084189*alphal[5]*fUpwind_l[69]+0.1581138830084189*alphal[12]*fUpwind_l[66]+0.1581138830084189*alphal[13]*fUpwind_l[65]+0.1581138830084189*alphal[20]*fUpwind_l[64]+0.1767766952966368*alphal[4]*fUpwind_l[28]+0.1767766952966368*alphal[3]*fUpwind_l[27]+0.1767766952966368*alphal[9]*fUpwind_l[24]+0.1767766952966368*alphal[10]*fUpwind_l[23]+0.1767766952966368*alphal[7]*fUpwind_l[22]+0.1767766952966368*alphal[8]*fUpwind_l[21]+0.1767766952966368*alphal[0]*fUpwind_l[20]+0.1767766952966368*fUpwind_l[0]*alphal[20]+0.1767766952966368*fUpwind_l[15]*alphal[17]+0.1767766952966368*fUpwind_l[14]*alphal[16]+0.1767766952966368*alphal[1]*fUpwind_l[13]+0.1767766952966368*fUpwind_l[1]*alphal[13]+0.1767766952966368*alphal[2]*fUpwind_l[12]+0.1767766952966368*fUpwind_l[2]*alphal[12]+0.1767766952966368*alphal[5]*fUpwind_l[6]+0.1767766952966368*fUpwind_l[5]*alphal[6]; 
  Ghat_l[21] = 0.1581138830084189*alphal[13]*fUpwind_l[75]+0.1581138830084189*alphal[20]*fUpwind_l[71]+0.1581138830084189*alphal[5]*fUpwind_l[70]+0.1581138830084189*alphal[12]*fUpwind_l[67]+0.1581138830084189*alphal[8]*fUpwind_l[44]+0.1581138830084189*alphal[16]*fUpwind_l[41]+0.1581138830084189*alphal[3]*fUpwind_l[40]+0.1581138830084189*alphal[7]*fUpwind_l[36]+0.1767766952966368*alphal[10]*fUpwind_l[31]+0.1767766952966368*alphal[17]*fUpwind_l[30]+0.1767766952966368*alphal[4]*fUpwind_l[29]+0.1767766952966368*alphal[2]*fUpwind_l[27]+0.1767766952966368*alphal[9]*fUpwind_l[25]+0.1767766952966368*alphal[6]*fUpwind_l[22]+0.1767766952966368*alphal[0]*fUpwind_l[21]+0.1767766952966368*alphal[8]*fUpwind_l[20]+0.1767766952966368*fUpwind_l[8]*alphal[20]+0.1767766952966368*alphal[13]*fUpwind_l[16]+0.1767766952966368*fUpwind_l[13]*alphal[16]+0.1767766952966368*alphal[1]*fUpwind_l[14]+0.1767766952966368*alphal[3]*fUpwind_l[12]+0.1767766952966368*fUpwind_l[3]*alphal[12]+0.1767766952966368*alphal[5]*fUpwind_l[7]+0.1767766952966368*fUpwind_l[5]*alphal[7]; 
  Ghat_l[22] = 0.1581138830084189*alphal[12]*fUpwind_l[75]+0.1581138830084189*alphal[5]*fUpwind_l[71]+0.1581138830084189*alphal[20]*fUpwind_l[70]+0.1581138830084189*alphal[13]*fUpwind_l[67]+0.1581138830084189*alphal[7]*fUpwind_l[44]+0.1581138830084189*alphal[3]*fUpwind_l[41]+0.1581138830084189*alphal[16]*fUpwind_l[40]+0.1581138830084189*alphal[8]*fUpwind_l[36]+0.1767766952966368*alphal[9]*fUpwind_l[31]+0.1767766952966368*alphal[4]*fUpwind_l[30]+0.1767766952966368*alphal[17]*fUpwind_l[29]+0.1767766952966368*alphal[1]*fUpwind_l[27]+0.1767766952966368*alphal[10]*fUpwind_l[25]+0.1767766952966368*alphal[0]*fUpwind_l[22]+0.1767766952966368*alphal[6]*fUpwind_l[21]+0.1767766952966368*alphal[7]*fUpwind_l[20]+0.1767766952966368*fUpwind_l[7]*alphal[20]+0.1767766952966368*alphal[12]*fUpwind_l[16]+0.1767766952966368*fUpwind_l[12]*alphal[16]+0.1767766952966368*alphal[2]*fUpwind_l[14]+0.1767766952966368*alphal[3]*fUpwind_l[13]+0.1767766952966368*fUpwind_l[3]*alphal[13]+0.1767766952966368*alphal[5]*fUpwind_l[8]+0.1767766952966368*fUpwind_l[5]*alphal[8]; 
  Ghat_l[23] = 0.1581138830084189*alphal[13]*fUpwind_l[76]+0.1581138830084189*alphal[20]*fUpwind_l[73]+0.1581138830084189*alphal[5]*fUpwind_l[72]+0.1581138830084189*alphal[12]*fUpwind_l[68]+0.1581138830084189*alphal[10]*fUpwind_l[60]+0.1581138830084189*alphal[17]*fUpwind_l[57]+0.1581138830084189*alphal[4]*fUpwind_l[56]+0.1581138830084189*alphal[9]*fUpwind_l[52]+0.1767766952966368*alphal[8]*fUpwind_l[31]+0.1767766952966368*alphal[16]*fUpwind_l[30]+0.1767766952966368*alphal[3]*fUpwind_l[29]+0.1767766952966368*alphal[2]*fUpwind_l[28]+0.1767766952966368*alphal[7]*fUpwind_l[25]+0.1767766952966368*alphal[6]*fUpwind_l[24]+0.1767766952966368*alphal[0]*fUpwind_l[23]+0.1767766952966368*alphal[10]*fUpwind_l[20]+0.1767766952966368*fUpwind_l[10]*alphal[20]+0.1767766952966368*alphal[13]*fUpwind_l[17]+0.1767766952966368*fUpwind_l[13]*alphal[17]+0.1767766952966368*alphal[1]*fUpwind_l[15]+0.1767766952966368*alphal[4]*fUpwind_l[12]+0.1767766952966368*fUpwind_l[4]*alphal[12]+0.1767766952966368*alphal[5]*fUpwind_l[9]+0.1767766952966368*fUpwind_l[5]*alphal[9]; 
  Ghat_l[24] = 0.1581138830084189*alphal[12]*fUpwind_l[76]+0.1581138830084189*alphal[5]*fUpwind_l[73]+0.1581138830084189*alphal[20]*fUpwind_l[72]+0.1581138830084189*alphal[13]*fUpwind_l[68]+0.1581138830084189*alphal[9]*fUpwind_l[60]+0.1581138830084189*alphal[4]*fUpwind_l[57]+0.1581138830084189*alphal[17]*fUpwind_l[56]+0.1581138830084189*alphal[10]*fUpwind_l[52]+0.1767766952966368*alphal[7]*fUpwind_l[31]+0.1767766952966368*alphal[3]*fUpwind_l[30]+0.1767766952966368*alphal[16]*fUpwind_l[29]+0.1767766952966368*alphal[1]*fUpwind_l[28]+0.1767766952966368*alphal[8]*fUpwind_l[25]+0.1767766952966368*alphal[0]*fUpwind_l[24]+0.1767766952966368*alphal[6]*fUpwind_l[23]+0.1767766952966368*alphal[9]*fUpwind_l[20]+0.1767766952966368*fUpwind_l[9]*alphal[20]+0.1767766952966368*alphal[12]*fUpwind_l[17]+0.1767766952966368*fUpwind_l[12]*alphal[17]+0.1767766952966368*alphal[2]*fUpwind_l[15]+0.1767766952966368*alphal[4]*fUpwind_l[13]+0.1767766952966368*fUpwind_l[4]*alphal[13]+0.1767766952966368*alphal[5]*fUpwind_l[10]+0.1767766952966368*fUpwind_l[5]*alphal[10]; 
  Ghat_l[25] = 0.1581138830084189*alphal[20]*fUpwind_l[79]+0.1581138830084189*alphal[13]*fUpwind_l[78]+0.1581138830084189*alphal[12]*fUpwind_l[77]+0.1581138830084189*alphal[5]*fUpwind_l[74]+0.1581138830084189*alphal[17]*fUpwind_l[63]+0.1581138830084189*alphal[10]*fUpwind_l[62]+0.1581138830084189*alphal[9]*fUpwind_l[61]+0.1581138830084189*alphal[4]*fUpwind_l[58]+0.1581138830084189*alphal[16]*fUpwind_l[47]+0.1581138830084189*alphal[8]*fUpwind_l[46]+0.1581138830084189*alphal[7]*fUpwind_l[45]+0.1581138830084189*alphal[3]*fUpwind_l[42]+0.1767766952966368*alphal[6]*fUpwind_l[31]+0.1767766952966368*alphal[2]*fUpwind_l[30]+0.1767766952966368*alphal[1]*fUpwind_l[29]+0.1767766952966368*alphal[16]*fUpwind_l[28]+0.1767766952966368*alphal[17]*fUpwind_l[27]+0.1767766952966368*alphal[20]*fUpwind_l[26]+0.1767766952966368*alphal[0]*fUpwind_l[25]+0.1767766952966368*alphal[8]*fUpwind_l[24]+0.1767766952966368*alphal[7]*fUpwind_l[23]+0.1767766952966368*alphal[10]*fUpwind_l[22]+0.1767766952966368*alphal[9]*fUpwind_l[21]+0.1767766952966368*alphal[13]*fUpwind_l[19]+0.1767766952966368*alphal[12]*fUpwind_l[18]+0.1767766952966368*alphal[3]*fUpwind_l[15]+0.1767766952966368*alphal[4]*fUpwind_l[14]+0.1767766952966368*alphal[5]*fUpwind_l[11]; 
  Ghat_l[26] = 0.1581138830084189*alphal[4]*fUpwind_l[59]+0.1581138830084189*alphal[9]*fUpwind_l[55]+0.1581138830084189*alphal[10]*fUpwind_l[54]+0.1581138830084189*alphal[17]*fUpwind_l[51]+0.1581138830084189*alphal[3]*fUpwind_l[43]+0.1581138830084189*alphal[7]*fUpwind_l[39]+0.1581138830084189*alphal[8]*fUpwind_l[38]+0.1581138830084189*alphal[16]*fUpwind_l[35]+0.1767766952966368*alphal[5]*fUpwind_l[31]+0.1767766952966368*alphal[12]*fUpwind_l[30]+0.1767766952966368*alphal[13]*fUpwind_l[29]+0.1767766952966368*alphal[0]*fUpwind_l[26]+0.1767766952966368*alphal[20]*fUpwind_l[25]+0.1767766952966368*alphal[1]*fUpwind_l[19]+0.1767766952966368*alphal[2]*fUpwind_l[18]+0.1767766952966368*alphal[3]*fUpwind_l[17]+0.1767766952966368*fUpwind_l[3]*alphal[17]+0.1767766952966368*alphal[4]*fUpwind_l[16]+0.1767766952966368*fUpwind_l[4]*alphal[16]+0.1767766952966368*alphal[6]*fUpwind_l[11]+0.1767766952966368*alphal[7]*fUpwind_l[10]+0.1767766952966368*fUpwind_l[7]*alphal[10]+0.1767766952966368*alphal[8]*fUpwind_l[9]+0.1767766952966368*fUpwind_l[8]*alphal[9]; 
  Ghat_l[27] = 0.1581138830084189*alphal[5]*fUpwind_l[75]+0.1581138830084189*alphal[12]*fUpwind_l[71]+0.1581138830084189*alphal[13]*fUpwind_l[70]+0.1581138830084189*alphal[20]*fUpwind_l[67]+0.1581138830084189*alphal[3]*fUpwind_l[44]+0.1581138830084189*alphal[7]*fUpwind_l[41]+0.1581138830084189*alphal[8]*fUpwind_l[40]+0.1581138830084189*alphal[16]*fUpwind_l[36]+0.1767766952966368*alphal[4]*fUpwind_l[31]+0.1767766952966368*alphal[9]*fUpwind_l[30]+0.1767766952966368*alphal[10]*fUpwind_l[29]+0.1767766952966368*alphal[0]*fUpwind_l[27]+0.1767766952966368*alphal[17]*fUpwind_l[25]+0.1767766952966368*alphal[1]*fUpwind_l[22]+0.1767766952966368*alphal[2]*fUpwind_l[21]+0.1767766952966368*alphal[3]*fUpwind_l[20]+0.1767766952966368*fUpwind_l[3]*alphal[20]+0.1767766952966368*alphal[5]*fUpwind_l[16]+0.1767766952966368*fUpwind_l[5]*alphal[16]+0.1767766952966368*alphal[6]*fUpwind_l[14]+0.1767766952966368*alphal[7]*fUpwind_l[13]+0.1767766952966368*fUpwind_l[7]*alphal[13]+0.1767766952966368*alphal[8]*fUpwind_l[12]+0.1767766952966368*fUpwind_l[8]*alphal[12]; 
  Ghat_l[28] = 0.1581138830084189*alphal[5]*fUpwind_l[76]+0.1581138830084189*alphal[12]*fUpwind_l[73]+0.1581138830084189*alphal[13]*fUpwind_l[72]+0.1581138830084189*alphal[20]*fUpwind_l[68]+0.1581138830084189*alphal[4]*fUpwind_l[60]+0.1581138830084189*alphal[9]*fUpwind_l[57]+0.1581138830084189*alphal[10]*fUpwind_l[56]+0.1581138830084189*alphal[17]*fUpwind_l[52]+0.1767766952966368*alphal[3]*fUpwind_l[31]+0.1767766952966368*alphal[7]*fUpwind_l[30]+0.1767766952966368*alphal[8]*fUpwind_l[29]+0.1767766952966368*alphal[0]*fUpwind_l[28]+0.1767766952966368*alphal[16]*fUpwind_l[25]+0.1767766952966368*alphal[1]*fUpwind_l[24]+0.1767766952966368*alphal[2]*fUpwind_l[23]+0.1767766952966368*alphal[4]*fUpwind_l[20]+0.1767766952966368*fUpwind_l[4]*alphal[20]+0.1767766952966368*alphal[5]*fUpwind_l[17]+0.1767766952966368*fUpwind_l[5]*alphal[17]+0.1767766952966368*alphal[6]*fUpwind_l[15]+0.1767766952966368*alphal[9]*fUpwind_l[13]+0.1767766952966368*fUpwind_l[9]*alphal[13]+0.1767766952966368*alphal[10]*fUpwind_l[12]+0.1767766952966368*fUpwind_l[10]*alphal[12]; 
  Ghat_l[29] = 0.1581138830084189*alphal[13]*fUpwind_l[79]+0.1581138830084189*alphal[20]*fUpwind_l[78]+0.1581138830084189*alphal[5]*fUpwind_l[77]+0.1581138830084189*alphal[12]*fUpwind_l[74]+0.1581138830084189*alphal[10]*fUpwind_l[63]+0.1581138830084189*alphal[17]*fUpwind_l[62]+0.1581138830084189*alphal[4]*fUpwind_l[61]+0.1581138830084189*alphal[9]*fUpwind_l[58]+0.1581138830084189*alphal[8]*fUpwind_l[47]+0.1581138830084189*alphal[16]*fUpwind_l[46]+0.1581138830084189*alphal[3]*fUpwind_l[45]+0.1581138830084189*alphal[7]*fUpwind_l[42]+0.1767766952966368*alphal[2]*fUpwind_l[31]+0.1767766952966368*alphal[6]*fUpwind_l[30]+0.1767766952966368*alphal[0]*fUpwind_l[29]+0.1767766952966368*alphal[8]*fUpwind_l[28]+0.1767766952966368*alphal[10]*fUpwind_l[27]+0.1767766952966368*alphal[13]*fUpwind_l[26]+0.1767766952966368*alphal[1]*fUpwind_l[25]+0.1767766952966368*alphal[16]*fUpwind_l[24]+0.1767766952966368*alphal[3]*fUpwind_l[23]+0.1767766952966368*alphal[17]*fUpwind_l[22]+0.1767766952966368*alphal[4]*fUpwind_l[21]+0.1767766952966368*fUpwind_l[19]*alphal[20]+0.1767766952966368*alphal[5]*fUpwind_l[18]+0.1767766952966368*alphal[7]*fUpwind_l[15]+0.1767766952966368*alphal[9]*fUpwind_l[14]+0.1767766952966368*fUpwind_l[11]*alphal[12]; 
  Ghat_l[30] = 0.1581138830084189*alphal[12]*fUpwind_l[79]+0.1581138830084189*alphal[5]*fUpwind_l[78]+0.1581138830084189*alphal[20]*fUpwind_l[77]+0.1581138830084189*alphal[13]*fUpwind_l[74]+0.1581138830084189*alphal[9]*fUpwind_l[63]+0.1581138830084189*alphal[4]*fUpwind_l[62]+0.1581138830084189*alphal[17]*fUpwind_l[61]+0.1581138830084189*alphal[10]*fUpwind_l[58]+0.1581138830084189*alphal[7]*fUpwind_l[47]+0.1581138830084189*alphal[3]*fUpwind_l[46]+0.1581138830084189*alphal[16]*fUpwind_l[45]+0.1581138830084189*alphal[8]*fUpwind_l[42]+0.1767766952966368*alphal[1]*fUpwind_l[31]+0.1767766952966368*alphal[0]*fUpwind_l[30]+0.1767766952966368*alphal[6]*fUpwind_l[29]+0.1767766952966368*alphal[7]*fUpwind_l[28]+0.1767766952966368*alphal[9]*fUpwind_l[27]+0.1767766952966368*alphal[12]*fUpwind_l[26]+0.1767766952966368*alphal[2]*fUpwind_l[25]+0.1767766952966368*alphal[3]*fUpwind_l[24]+0.1767766952966368*alphal[16]*fUpwind_l[23]+0.1767766952966368*alphal[4]*fUpwind_l[22]+0.1767766952966368*alphal[17]*fUpwind_l[21]+0.1767766952966368*fUpwind_l[18]*alphal[20]+0.1767766952966368*alphal[5]*fUpwind_l[19]+0.1767766952966368*alphal[8]*fUpwind_l[15]+0.1767766952966368*alphal[10]*fUpwind_l[14]+0.1767766952966368*fUpwind_l[11]*alphal[13]; 
  Ghat_l[31] = 0.1581138830084189*alphal[5]*fUpwind_l[79]+0.1581138830084189*alphal[12]*fUpwind_l[78]+0.1581138830084189*alphal[13]*fUpwind_l[77]+0.1581138830084189*alphal[20]*fUpwind_l[74]+0.1581138830084189*alphal[4]*fUpwind_l[63]+0.1581138830084189*alphal[9]*fUpwind_l[62]+0.1581138830084189*alphal[10]*fUpwind_l[61]+0.1581138830084189*alphal[17]*fUpwind_l[58]+0.1581138830084189*alphal[3]*fUpwind_l[47]+0.1581138830084189*alphal[7]*fUpwind_l[46]+0.1581138830084189*alphal[8]*fUpwind_l[45]+0.1581138830084189*alphal[16]*fUpwind_l[42]+0.1767766952966368*alphal[0]*fUpwind_l[31]+0.1767766952966368*alphal[1]*fUpwind_l[30]+0.1767766952966368*alphal[2]*fUpwind_l[29]+0.1767766952966368*alphal[3]*fUpwind_l[28]+0.1767766952966368*alphal[4]*fUpwind_l[27]+0.1767766952966368*alphal[5]*fUpwind_l[26]+0.1767766952966368*alphal[6]*fUpwind_l[25]+0.1767766952966368*alphal[7]*fUpwind_l[24]+0.1767766952966368*alphal[8]*fUpwind_l[23]+0.1767766952966368*alphal[9]*fUpwind_l[22]+0.1767766952966368*alphal[10]*fUpwind_l[21]+0.1767766952966368*fUpwind_l[11]*alphal[20]+0.1767766952966368*alphal[12]*fUpwind_l[19]+0.1767766952966368*alphal[13]*fUpwind_l[18]+0.1767766952966368*fUpwind_l[14]*alphal[17]+0.1767766952966368*fUpwind_l[15]*alphal[16]; 
  Ghat_l[32] = 0.1767766952966368*alphal[20]*fUpwind_l[44]+0.1767766952966368*alphal[17]*fUpwind_l[43]+0.1767766952966368*alphal[13]*fUpwind_l[41]+0.1767766952966368*alphal[12]*fUpwind_l[40]+0.1767766952966368*alphal[10]*fUpwind_l[39]+0.1767766952966368*alphal[9]*fUpwind_l[38]+0.1767766952966368*alphal[6]*fUpwind_l[37]+0.1767766952966368*alphal[5]*fUpwind_l[36]+0.1767766952966368*alphal[4]*fUpwind_l[35]+0.1767766952966368*alphal[2]*fUpwind_l[34]+0.1767766952966368*alphal[1]*fUpwind_l[33]+0.1767766952966368*alphal[0]*fUpwind_l[32]+0.1581138830084189*alphal[16]*fUpwind_l[16]+0.1581138830084189*alphal[8]*fUpwind_l[8]+0.1581138830084189*alphal[7]*fUpwind_l[7]+0.1581138830084189*alphal[3]*fUpwind_l[3]; 
  Ghat_l[33] = 0.1767766952966368*alphal[13]*fUpwind_l[44]+0.1767766952966368*alphal[10]*fUpwind_l[43]+0.1767766952966368*alphal[20]*fUpwind_l[41]+0.1767766952966368*alphal[5]*fUpwind_l[40]+0.1767766952966368*alphal[17]*fUpwind_l[39]+0.1767766952966368*alphal[4]*fUpwind_l[38]+0.1767766952966368*alphal[2]*fUpwind_l[37]+0.1767766952966368*alphal[12]*fUpwind_l[36]+0.1767766952966368*alphal[9]*fUpwind_l[35]+0.1767766952966368*alphal[6]*fUpwind_l[34]+0.1767766952966368*alphal[0]*fUpwind_l[33]+0.1767766952966368*alphal[1]*fUpwind_l[32]+0.1581138830084189*alphal[8]*fUpwind_l[16]+0.1581138830084189*fUpwind_l[8]*alphal[16]+0.1581138830084189*alphal[3]*fUpwind_l[7]+0.1581138830084189*fUpwind_l[3]*alphal[7]; 
  Ghat_l[34] = 0.1767766952966368*alphal[12]*fUpwind_l[44]+0.1767766952966368*alphal[9]*fUpwind_l[43]+0.1767766952966368*alphal[5]*fUpwind_l[41]+0.1767766952966368*alphal[20]*fUpwind_l[40]+0.1767766952966368*alphal[4]*fUpwind_l[39]+0.1767766952966368*alphal[17]*fUpwind_l[38]+0.1767766952966368*alphal[1]*fUpwind_l[37]+0.1767766952966368*alphal[13]*fUpwind_l[36]+0.1767766952966368*alphal[10]*fUpwind_l[35]+0.1767766952966368*alphal[0]*fUpwind_l[34]+0.1767766952966368*alphal[6]*fUpwind_l[33]+0.1767766952966368*alphal[2]*fUpwind_l[32]+0.1581138830084189*alphal[7]*fUpwind_l[16]+0.1581138830084189*fUpwind_l[7]*alphal[16]+0.1581138830084189*alphal[3]*fUpwind_l[8]+0.1581138830084189*fUpwind_l[3]*alphal[8]; 
  Ghat_l[35] = 0.1767766952966368*alphal[20]*fUpwind_l[47]+0.1767766952966368*alphal[13]*fUpwind_l[46]+0.1767766952966368*alphal[12]*fUpwind_l[45]+0.1767766952966368*alphal[6]*fUpwind_l[43]+0.1767766952966368*alphal[5]*fUpwind_l[42]+0.1767766952966368*alphal[2]*fUpwind_l[39]+0.1767766952966368*alphal[1]*fUpwind_l[38]+0.1767766952966368*alphal[17]*fUpwind_l[37]+0.1767766952966368*alphal[0]*fUpwind_l[35]+0.1767766952966368*alphal[10]*fUpwind_l[34]+0.1767766952966368*alphal[9]*fUpwind_l[33]+0.1767766952966368*alphal[4]*fUpwind_l[32]+0.1581138830084189*alphal[16]*fUpwind_l[26]+0.1581138830084189*alphal[8]*fUpwind_l[19]+0.1581138830084189*alphal[7]*fUpwind_l[18]+0.1581138830084189*alphal[3]*fUpwind_l[11]; 
  Ghat_l[36] = 0.1767766952966368*alphal[17]*fUpwind_l[47]+0.1767766952966368*alphal[10]*fUpwind_l[46]+0.1767766952966368*alphal[9]*fUpwind_l[45]+0.1767766952966368*alphal[6]*fUpwind_l[44]+0.1767766952966368*alphal[4]*fUpwind_l[42]+0.1767766952966368*alphal[2]*fUpwind_l[41]+0.1767766952966368*alphal[1]*fUpwind_l[40]+0.1767766952966368*alphal[20]*fUpwind_l[37]+0.1767766952966368*alphal[0]*fUpwind_l[36]+0.1767766952966368*alphal[13]*fUpwind_l[34]+0.1767766952966368*alphal[12]*fUpwind_l[33]+0.1767766952966368*alphal[5]*fUpwind_l[32]+0.1581138830084189*alphal[16]*fUpwind_l[27]+0.1581138830084189*alphal[8]*fUpwind_l[22]+0.1581138830084189*alphal[7]*fUpwind_l[21]+0.1581138830084189*alphal[3]*fUpwind_l[14]; 
  Ghat_l[37] = 0.1767766952966368*alphal[5]*fUpwind_l[44]+0.1767766952966368*alphal[4]*fUpwind_l[43]+0.1767766952966368*alphal[12]*fUpwind_l[41]+0.1767766952966368*alphal[13]*fUpwind_l[40]+0.1767766952966368*alphal[9]*fUpwind_l[39]+0.1767766952966368*alphal[10]*fUpwind_l[38]+0.1767766952966368*alphal[0]*fUpwind_l[37]+0.1767766952966368*alphal[20]*fUpwind_l[36]+0.1767766952966368*alphal[17]*fUpwind_l[35]+0.1767766952966368*alphal[1]*fUpwind_l[34]+0.1767766952966368*alphal[2]*fUpwind_l[33]+0.1767766952966368*alphal[6]*fUpwind_l[32]+0.1581138830084189*alphal[3]*fUpwind_l[16]+0.1581138830084189*fUpwind_l[3]*alphal[16]+0.1581138830084189*alphal[7]*fUpwind_l[8]+0.1581138830084189*fUpwind_l[7]*alphal[8]; 
  Ghat_l[38] = 0.1767766952966368*alphal[13]*fUpwind_l[47]+0.1767766952966368*alphal[20]*fUpwind_l[46]+0.1767766952966368*alphal[5]*fUpwind_l[45]+0.1767766952966368*alphal[2]*fUpwind_l[43]+0.1767766952966368*alphal[12]*fUpwind_l[42]+0.1767766952966368*alphal[6]*fUpwind_l[39]+0.1767766952966368*alphal[0]*fUpwind_l[38]+0.1767766952966368*alphal[10]*fUpwind_l[37]+0.1767766952966368*alphal[1]*fUpwind_l[35]+0.1767766952966368*alphal[17]*fUpwind_l[34]+0.1767766952966368*alphal[4]*fUpwind_l[33]+0.1767766952966368*alphal[9]*fUpwind_l[32]+0.1581138830084189*alphal[8]*fUpwind_l[26]+0.1581138830084189*alphal[16]*fUpwind_l[19]+0.1581138830084189*alphal[3]*fUpwind_l[18]+0.1581138830084189*alphal[7]*fUpwind_l[11]; 
  Ghat_l[39] = 0.1767766952966368*alphal[12]*fUpwind_l[47]+0.1767766952966368*alphal[5]*fUpwind_l[46]+0.1767766952966368*alphal[20]*fUpwind_l[45]+0.1767766952966368*alphal[1]*fUpwind_l[43]+0.1767766952966368*alphal[13]*fUpwind_l[42]+0.1767766952966368*alphal[0]*fUpwind_l[39]+0.1767766952966368*alphal[6]*fUpwind_l[38]+0.1767766952966368*alphal[9]*fUpwind_l[37]+0.1767766952966368*alphal[2]*fUpwind_l[35]+0.1767766952966368*alphal[4]*fUpwind_l[34]+0.1767766952966368*alphal[17]*fUpwind_l[33]+0.1767766952966368*alphal[10]*fUpwind_l[32]+0.1581138830084189*alphal[7]*fUpwind_l[26]+0.1581138830084189*alphal[3]*fUpwind_l[19]+0.1581138830084189*alphal[16]*fUpwind_l[18]+0.1581138830084189*alphal[8]*fUpwind_l[11]; 
  Ghat_l[40] = 0.1767766952966368*alphal[10]*fUpwind_l[47]+0.1767766952966368*alphal[17]*fUpwind_l[46]+0.1767766952966368*alphal[4]*fUpwind_l[45]+0.1767766952966368*alphal[2]*fUpwind_l[44]+0.1767766952966368*alphal[9]*fUpwind_l[42]+0.1767766952966368*alphal[6]*fUpwind_l[41]+0.1767766952966368*alphal[0]*fUpwind_l[40]+0.1767766952966368*alphal[13]*fUpwind_l[37]+0.1767766952966368*alphal[1]*fUpwind_l[36]+0.1767766952966368*alphal[20]*fUpwind_l[34]+0.1767766952966368*alphal[5]*fUpwind_l[33]+0.1767766952966368*alphal[12]*fUpwind_l[32]+0.1581138830084189*alphal[8]*fUpwind_l[27]+0.1581138830084189*alphal[16]*fUpwind_l[22]+0.1581138830084189*alphal[3]*fUpwind_l[21]+0.1581138830084189*alphal[7]*fUpwind_l[14]; 
  Ghat_l[41] = 0.1767766952966368*alphal[9]*fUpwind_l[47]+0.1767766952966368*alphal[4]*fUpwind_l[46]+0.1767766952966368*alphal[17]*fUpwind_l[45]+0.1767766952966368*alphal[1]*fUpwind_l[44]+0.1767766952966368*alphal[10]*fUpwind_l[42]+0.1767766952966368*alphal[0]*fUpwind_l[41]+0.1767766952966368*alphal[6]*fUpwind_l[40]+0.1767766952966368*alphal[12]*fUpwind_l[37]+0.1767766952966368*alphal[2]*fUpwind_l[36]+0.1767766952966368*alphal[5]*fUpwind_l[34]+0.1767766952966368*alphal[20]*fUpwind_l[33]+0.1767766952966368*alphal[13]*fUpwind_l[32]+0.1581138830084189*alphal[7]*fUpwind_l[27]+0.1581138830084189*alphal[3]*fUpwind_l[22]+0.1581138830084189*alphal[16]*fUpwind_l[21]+0.1581138830084189*alphal[8]*fUpwind_l[14]; 
  Ghat_l[42] = 0.1767766952966368*alphal[6]*fUpwind_l[47]+0.1767766952966368*alphal[2]*fUpwind_l[46]+0.1767766952966368*alphal[1]*fUpwind_l[45]+0.1767766952966368*alphal[17]*fUpwind_l[44]+0.1767766952966368*alphal[20]*fUpwind_l[43]+0.1767766952966368*alphal[0]*fUpwind_l[42]+0.1767766952966368*alphal[10]*fUpwind_l[41]+0.1767766952966368*alphal[9]*fUpwind_l[40]+0.1767766952966368*alphal[13]*fUpwind_l[39]+0.1767766952966368*alphal[12]*fUpwind_l[38]+0.1767766952966368*alphal[4]*fUpwind_l[36]+0.1767766952966368*alphal[5]*fUpwind_l[35]+0.1581138830084189*alphal[16]*fUpwind_l[31]+0.1581138830084189*alphal[8]*fUpwind_l[30]+0.1581138830084189*alphal[7]*fUpwind_l[29]+0.1581138830084189*alphal[3]*fUpwind_l[25]; 
  Ghat_l[43] = 0.1767766952966368*alphal[5]*fUpwind_l[47]+0.1767766952966368*alphal[12]*fUpwind_l[46]+0.1767766952966368*alphal[13]*fUpwind_l[45]+0.1767766952966368*alphal[0]*fUpwind_l[43]+0.1767766952966368*alphal[20]*fUpwind_l[42]+0.1767766952966368*alphal[1]*fUpwind_l[39]+0.1767766952966368*alphal[2]*fUpwind_l[38]+0.1767766952966368*alphal[4]*fUpwind_l[37]+0.1767766952966368*alphal[6]*fUpwind_l[35]+0.1767766952966368*alphal[9]*fUpwind_l[34]+0.1767766952966368*alphal[10]*fUpwind_l[33]+0.1767766952966368*alphal[17]*fUpwind_l[32]+0.1581138830084189*alphal[3]*fUpwind_l[26]+0.1581138830084189*alphal[7]*fUpwind_l[19]+0.1581138830084189*alphal[8]*fUpwind_l[18]+0.1581138830084189*fUpwind_l[11]*alphal[16]; 
  Ghat_l[44] = 0.1767766952966368*alphal[4]*fUpwind_l[47]+0.1767766952966368*alphal[9]*fUpwind_l[46]+0.1767766952966368*alphal[10]*fUpwind_l[45]+0.1767766952966368*alphal[0]*fUpwind_l[44]+0.1767766952966368*alphal[17]*fUpwind_l[42]+0.1767766952966368*alphal[1]*fUpwind_l[41]+0.1767766952966368*alphal[2]*fUpwind_l[40]+0.1767766952966368*alphal[5]*fUpwind_l[37]+0.1767766952966368*alphal[6]*fUpwind_l[36]+0.1767766952966368*alphal[12]*fUpwind_l[34]+0.1767766952966368*alphal[13]*fUpwind_l[33]+0.1767766952966368*alphal[20]*fUpwind_l[32]+0.1581138830084189*alphal[3]*fUpwind_l[27]+0.1581138830084189*alphal[7]*fUpwind_l[22]+0.1581138830084189*alphal[8]*fUpwind_l[21]+0.1581138830084189*fUpwind_l[14]*alphal[16]; 
  Ghat_l[45] = 0.1767766952966368*alphal[2]*fUpwind_l[47]+0.1767766952966368*alphal[6]*fUpwind_l[46]+0.1767766952966368*alphal[0]*fUpwind_l[45]+0.1767766952966368*alphal[10]*fUpwind_l[44]+0.1767766952966368*alphal[13]*fUpwind_l[43]+0.1767766952966368*alphal[1]*fUpwind_l[42]+0.1767766952966368*alphal[17]*fUpwind_l[41]+0.1767766952966368*alphal[4]*fUpwind_l[40]+0.1767766952966368*alphal[20]*fUpwind_l[39]+0.1767766952966368*alphal[5]*fUpwind_l[38]+0.1767766952966368*alphal[9]*fUpwind_l[36]+0.1767766952966368*alphal[12]*fUpwind_l[35]+0.1581138830084189*alphal[8]*fUpwind_l[31]+0.1581138830084189*alphal[16]*fUpwind_l[30]+0.1581138830084189*alphal[3]*fUpwind_l[29]+0.1581138830084189*alphal[7]*fUpwind_l[25]; 
  Ghat_l[46] = 0.1767766952966368*alphal[1]*fUpwind_l[47]+0.1767766952966368*alphal[0]*fUpwind_l[46]+0.1767766952966368*alphal[6]*fUpwind_l[45]+0.1767766952966368*alphal[9]*fUpwind_l[44]+0.1767766952966368*alphal[12]*fUpwind_l[43]+0.1767766952966368*alphal[2]*fUpwind_l[42]+0.1767766952966368*alphal[4]*fUpwind_l[41]+0.1767766952966368*alphal[17]*fUpwind_l[40]+0.1767766952966368*alphal[5]*fUpwind_l[39]+0.1767766952966368*alphal[20]*fUpwind_l[38]+0.1767766952966368*alphal[10]*fUpwind_l[36]+0.1767766952966368*alphal[13]*fUpwind_l[35]+0.1581138830084189*alphal[7]*fUpwind_l[31]+0.1581138830084189*alphal[3]*fUpwind_l[30]+0.1581138830084189*alphal[16]*fUpwind_l[29]+0.1581138830084189*alphal[8]*fUpwind_l[25]; 
  Ghat_l[47] = 0.1767766952966368*alphal[0]*fUpwind_l[47]+0.1767766952966368*alphal[1]*fUpwind_l[46]+0.1767766952966368*alphal[2]*fUpwind_l[45]+0.1767766952966368*alphal[4]*fUpwind_l[44]+0.1767766952966368*alphal[5]*fUpwind_l[43]+0.1767766952966368*alphal[6]*fUpwind_l[42]+0.1767766952966368*alphal[9]*fUpwind_l[41]+0.1767766952966368*alphal[10]*fUpwind_l[40]+0.1767766952966368*alphal[12]*fUpwind_l[39]+0.1767766952966368*alphal[13]*fUpwind_l[38]+0.1767766952966368*alphal[17]*fUpwind_l[36]+0.1767766952966368*alphal[20]*fUpwind_l[35]+0.1581138830084189*alphal[3]*fUpwind_l[31]+0.1581138830084189*alphal[7]*fUpwind_l[30]+0.1581138830084189*alphal[8]*fUpwind_l[29]+0.1581138830084189*alphal[16]*fUpwind_l[25]; 
  Ghat_l[48] = 0.1767766952966368*alphal[20]*fUpwind_l[60]+0.1767766952966368*alphal[16]*fUpwind_l[59]+0.1767766952966368*alphal[13]*fUpwind_l[57]+0.1767766952966368*alphal[12]*fUpwind_l[56]+0.1767766952966368*alphal[8]*fUpwind_l[55]+0.1767766952966368*alphal[7]*fUpwind_l[54]+0.1767766952966368*alphal[6]*fUpwind_l[53]+0.1767766952966368*alphal[5]*fUpwind_l[52]+0.1767766952966368*alphal[3]*fUpwind_l[51]+0.1767766952966368*alphal[2]*fUpwind_l[50]+0.1767766952966368*alphal[1]*fUpwind_l[49]+0.1767766952966368*alphal[0]*fUpwind_l[48]+0.1581138830084189*alphal[17]*fUpwind_l[17]+0.1581138830084189*alphal[10]*fUpwind_l[10]+0.1581138830084189*alphal[9]*fUpwind_l[9]+0.1581138830084189*alphal[4]*fUpwind_l[4]; 
  Ghat_l[49] = 0.1767766952966368*alphal[13]*fUpwind_l[60]+0.1767766952966368*alphal[8]*fUpwind_l[59]+0.1767766952966368*alphal[20]*fUpwind_l[57]+0.1767766952966368*alphal[5]*fUpwind_l[56]+0.1767766952966368*alphal[16]*fUpwind_l[55]+0.1767766952966368*alphal[3]*fUpwind_l[54]+0.1767766952966368*alphal[2]*fUpwind_l[53]+0.1767766952966368*alphal[12]*fUpwind_l[52]+0.1767766952966368*alphal[7]*fUpwind_l[51]+0.1767766952966368*alphal[6]*fUpwind_l[50]+0.1767766952966368*alphal[0]*fUpwind_l[49]+0.1767766952966368*alphal[1]*fUpwind_l[48]+0.1581138830084189*alphal[10]*fUpwind_l[17]+0.1581138830084189*fUpwind_l[10]*alphal[17]+0.1581138830084189*alphal[4]*fUpwind_l[9]+0.1581138830084189*fUpwind_l[4]*alphal[9]; 
  Ghat_l[50] = 0.1767766952966368*alphal[12]*fUpwind_l[60]+0.1767766952966368*alphal[7]*fUpwind_l[59]+0.1767766952966368*alphal[5]*fUpwind_l[57]+0.1767766952966368*alphal[20]*fUpwind_l[56]+0.1767766952966368*alphal[3]*fUpwind_l[55]+0.1767766952966368*alphal[16]*fUpwind_l[54]+0.1767766952966368*alphal[1]*fUpwind_l[53]+0.1767766952966368*alphal[13]*fUpwind_l[52]+0.1767766952966368*alphal[8]*fUpwind_l[51]+0.1767766952966368*alphal[0]*fUpwind_l[50]+0.1767766952966368*alphal[6]*fUpwind_l[49]+0.1767766952966368*alphal[2]*fUpwind_l[48]+0.1581138830084189*alphal[9]*fUpwind_l[17]+0.1581138830084189*fUpwind_l[9]*alphal[17]+0.1581138830084189*alphal[4]*fUpwind_l[10]+0.1581138830084189*fUpwind_l[4]*alphal[10]; 
  Ghat_l[51] = 0.1767766952966368*alphal[20]*fUpwind_l[63]+0.1767766952966368*alphal[13]*fUpwind_l[62]+0.1767766952966368*alphal[12]*fUpwind_l[61]+0.1767766952966368*alphal[6]*fUpwind_l[59]+0.1767766952966368*alphal[5]*fUpwind_l[58]+0.1767766952966368*alphal[2]*fUpwind_l[55]+0.1767766952966368*alphal[1]*fUpwind_l[54]+0.1767766952966368*alphal[16]*fUpwind_l[53]+0.1767766952966368*alphal[0]*fUpwind_l[51]+0.1767766952966368*alphal[8]*fUpwind_l[50]+0.1767766952966368*alphal[7]*fUpwind_l[49]+0.1767766952966368*alphal[3]*fUpwind_l[48]+0.1581138830084189*alphal[17]*fUpwind_l[26]+0.1581138830084189*alphal[10]*fUpwind_l[19]+0.1581138830084189*alphal[9]*fUpwind_l[18]+0.1581138830084189*alphal[4]*fUpwind_l[11]; 
  Ghat_l[52] = 0.1767766952966368*alphal[16]*fUpwind_l[63]+0.1767766952966368*alphal[8]*fUpwind_l[62]+0.1767766952966368*alphal[7]*fUpwind_l[61]+0.1767766952966368*alphal[6]*fUpwind_l[60]+0.1767766952966368*alphal[3]*fUpwind_l[58]+0.1767766952966368*alphal[2]*fUpwind_l[57]+0.1767766952966368*alphal[1]*fUpwind_l[56]+0.1767766952966368*alphal[20]*fUpwind_l[53]+0.1767766952966368*alphal[0]*fUpwind_l[52]+0.1767766952966368*alphal[13]*fUpwind_l[50]+0.1767766952966368*alphal[12]*fUpwind_l[49]+0.1767766952966368*alphal[5]*fUpwind_l[48]+0.1581138830084189*alphal[17]*fUpwind_l[28]+0.1581138830084189*alphal[10]*fUpwind_l[24]+0.1581138830084189*alphal[9]*fUpwind_l[23]+0.1581138830084189*alphal[4]*fUpwind_l[15]; 
  Ghat_l[53] = 0.1767766952966368*alphal[5]*fUpwind_l[60]+0.1767766952966368*alphal[3]*fUpwind_l[59]+0.1767766952966368*alphal[12]*fUpwind_l[57]+0.1767766952966368*alphal[13]*fUpwind_l[56]+0.1767766952966368*alphal[7]*fUpwind_l[55]+0.1767766952966368*alphal[8]*fUpwind_l[54]+0.1767766952966368*alphal[0]*fUpwind_l[53]+0.1767766952966368*alphal[20]*fUpwind_l[52]+0.1767766952966368*alphal[16]*fUpwind_l[51]+0.1767766952966368*alphal[1]*fUpwind_l[50]+0.1767766952966368*alphal[2]*fUpwind_l[49]+0.1767766952966368*alphal[6]*fUpwind_l[48]+0.1581138830084189*alphal[4]*fUpwind_l[17]+0.1581138830084189*fUpwind_l[4]*alphal[17]+0.1581138830084189*alphal[9]*fUpwind_l[10]+0.1581138830084189*fUpwind_l[9]*alphal[10]; 
  Ghat_l[54] = 0.1767766952966368*alphal[13]*fUpwind_l[63]+0.1767766952966368*alphal[20]*fUpwind_l[62]+0.1767766952966368*alphal[5]*fUpwind_l[61]+0.1767766952966368*alphal[2]*fUpwind_l[59]+0.1767766952966368*alphal[12]*fUpwind_l[58]+0.1767766952966368*alphal[6]*fUpwind_l[55]+0.1767766952966368*alphal[0]*fUpwind_l[54]+0.1767766952966368*alphal[8]*fUpwind_l[53]+0.1767766952966368*alphal[1]*fUpwind_l[51]+0.1767766952966368*alphal[16]*fUpwind_l[50]+0.1767766952966368*alphal[3]*fUpwind_l[49]+0.1767766952966368*alphal[7]*fUpwind_l[48]+0.1581138830084189*alphal[10]*fUpwind_l[26]+0.1581138830084189*alphal[17]*fUpwind_l[19]+0.1581138830084189*alphal[4]*fUpwind_l[18]+0.1581138830084189*alphal[9]*fUpwind_l[11]; 
  Ghat_l[55] = 0.1767766952966368*alphal[12]*fUpwind_l[63]+0.1767766952966368*alphal[5]*fUpwind_l[62]+0.1767766952966368*alphal[20]*fUpwind_l[61]+0.1767766952966368*alphal[1]*fUpwind_l[59]+0.1767766952966368*alphal[13]*fUpwind_l[58]+0.1767766952966368*alphal[0]*fUpwind_l[55]+0.1767766952966368*alphal[6]*fUpwind_l[54]+0.1767766952966368*alphal[7]*fUpwind_l[53]+0.1767766952966368*alphal[2]*fUpwind_l[51]+0.1767766952966368*alphal[3]*fUpwind_l[50]+0.1767766952966368*alphal[16]*fUpwind_l[49]+0.1767766952966368*alphal[8]*fUpwind_l[48]+0.1581138830084189*alphal[9]*fUpwind_l[26]+0.1581138830084189*alphal[4]*fUpwind_l[19]+0.1581138830084189*alphal[17]*fUpwind_l[18]+0.1581138830084189*alphal[10]*fUpwind_l[11]; 
  Ghat_l[56] = 0.1767766952966368*alphal[8]*fUpwind_l[63]+0.1767766952966368*alphal[16]*fUpwind_l[62]+0.1767766952966368*alphal[3]*fUpwind_l[61]+0.1767766952966368*alphal[2]*fUpwind_l[60]+0.1767766952966368*alphal[7]*fUpwind_l[58]+0.1767766952966368*alphal[6]*fUpwind_l[57]+0.1767766952966368*alphal[0]*fUpwind_l[56]+0.1767766952966368*alphal[13]*fUpwind_l[53]+0.1767766952966368*alphal[1]*fUpwind_l[52]+0.1767766952966368*alphal[20]*fUpwind_l[50]+0.1767766952966368*alphal[5]*fUpwind_l[49]+0.1767766952966368*alphal[12]*fUpwind_l[48]+0.1581138830084189*alphal[10]*fUpwind_l[28]+0.1581138830084189*alphal[17]*fUpwind_l[24]+0.1581138830084189*alphal[4]*fUpwind_l[23]+0.1581138830084189*alphal[9]*fUpwind_l[15]; 
  Ghat_l[57] = 0.1767766952966368*alphal[7]*fUpwind_l[63]+0.1767766952966368*alphal[3]*fUpwind_l[62]+0.1767766952966368*alphal[16]*fUpwind_l[61]+0.1767766952966368*alphal[1]*fUpwind_l[60]+0.1767766952966368*alphal[8]*fUpwind_l[58]+0.1767766952966368*alphal[0]*fUpwind_l[57]+0.1767766952966368*alphal[6]*fUpwind_l[56]+0.1767766952966368*alphal[12]*fUpwind_l[53]+0.1767766952966368*alphal[2]*fUpwind_l[52]+0.1767766952966368*alphal[5]*fUpwind_l[50]+0.1767766952966368*alphal[20]*fUpwind_l[49]+0.1767766952966368*alphal[13]*fUpwind_l[48]+0.1581138830084189*alphal[9]*fUpwind_l[28]+0.1581138830084189*alphal[4]*fUpwind_l[24]+0.1581138830084189*alphal[17]*fUpwind_l[23]+0.1581138830084189*alphal[10]*fUpwind_l[15]; 
  Ghat_l[58] = 0.1767766952966368*alphal[6]*fUpwind_l[63]+0.1767766952966368*alphal[2]*fUpwind_l[62]+0.1767766952966368*alphal[1]*fUpwind_l[61]+0.1767766952966368*alphal[16]*fUpwind_l[60]+0.1767766952966368*alphal[20]*fUpwind_l[59]+0.1767766952966368*alphal[0]*fUpwind_l[58]+0.1767766952966368*alphal[8]*fUpwind_l[57]+0.1767766952966368*alphal[7]*fUpwind_l[56]+0.1767766952966368*alphal[13]*fUpwind_l[55]+0.1767766952966368*alphal[12]*fUpwind_l[54]+0.1767766952966368*alphal[3]*fUpwind_l[52]+0.1767766952966368*alphal[5]*fUpwind_l[51]+0.1581138830084189*alphal[17]*fUpwind_l[31]+0.1581138830084189*alphal[10]*fUpwind_l[30]+0.1581138830084189*alphal[9]*fUpwind_l[29]+0.1581138830084189*alphal[4]*fUpwind_l[25]; 
  Ghat_l[59] = 0.1767766952966368*alphal[5]*fUpwind_l[63]+0.1767766952966368*alphal[12]*fUpwind_l[62]+0.1767766952966368*alphal[13]*fUpwind_l[61]+0.1767766952966368*alphal[0]*fUpwind_l[59]+0.1767766952966368*alphal[20]*fUpwind_l[58]+0.1767766952966368*alphal[1]*fUpwind_l[55]+0.1767766952966368*alphal[2]*fUpwind_l[54]+0.1767766952966368*alphal[3]*fUpwind_l[53]+0.1767766952966368*alphal[6]*fUpwind_l[51]+0.1767766952966368*alphal[7]*fUpwind_l[50]+0.1767766952966368*alphal[8]*fUpwind_l[49]+0.1767766952966368*alphal[16]*fUpwind_l[48]+0.1581138830084189*alphal[4]*fUpwind_l[26]+0.1581138830084189*alphal[9]*fUpwind_l[19]+0.1581138830084189*alphal[10]*fUpwind_l[18]+0.1581138830084189*fUpwind_l[11]*alphal[17]; 
  Ghat_l[60] = 0.1767766952966368*alphal[3]*fUpwind_l[63]+0.1767766952966368*alphal[7]*fUpwind_l[62]+0.1767766952966368*alphal[8]*fUpwind_l[61]+0.1767766952966368*alphal[0]*fUpwind_l[60]+0.1767766952966368*alphal[16]*fUpwind_l[58]+0.1767766952966368*alphal[1]*fUpwind_l[57]+0.1767766952966368*alphal[2]*fUpwind_l[56]+0.1767766952966368*alphal[5]*fUpwind_l[53]+0.1767766952966368*alphal[6]*fUpwind_l[52]+0.1767766952966368*alphal[12]*fUpwind_l[50]+0.1767766952966368*alphal[13]*fUpwind_l[49]+0.1767766952966368*alphal[20]*fUpwind_l[48]+0.1581138830084189*alphal[4]*fUpwind_l[28]+0.1581138830084189*alphal[9]*fUpwind_l[24]+0.1581138830084189*alphal[10]*fUpwind_l[23]+0.1581138830084189*fUpwind_l[15]*alphal[17]; 
  Ghat_l[61] = 0.1767766952966368*alphal[2]*fUpwind_l[63]+0.1767766952966368*alphal[6]*fUpwind_l[62]+0.1767766952966368*alphal[0]*fUpwind_l[61]+0.1767766952966368*alphal[8]*fUpwind_l[60]+0.1767766952966368*alphal[13]*fUpwind_l[59]+0.1767766952966368*alphal[1]*fUpwind_l[58]+0.1767766952966368*alphal[16]*fUpwind_l[57]+0.1767766952966368*alphal[3]*fUpwind_l[56]+0.1767766952966368*alphal[20]*fUpwind_l[55]+0.1767766952966368*alphal[5]*fUpwind_l[54]+0.1767766952966368*alphal[7]*fUpwind_l[52]+0.1767766952966368*alphal[12]*fUpwind_l[51]+0.1581138830084189*alphal[10]*fUpwind_l[31]+0.1581138830084189*alphal[17]*fUpwind_l[30]+0.1581138830084189*alphal[4]*fUpwind_l[29]+0.1581138830084189*alphal[9]*fUpwind_l[25]; 
  Ghat_l[62] = 0.1767766952966368*alphal[1]*fUpwind_l[63]+0.1767766952966368*alphal[0]*fUpwind_l[62]+0.1767766952966368*alphal[6]*fUpwind_l[61]+0.1767766952966368*alphal[7]*fUpwind_l[60]+0.1767766952966368*alphal[12]*fUpwind_l[59]+0.1767766952966368*alphal[2]*fUpwind_l[58]+0.1767766952966368*alphal[3]*fUpwind_l[57]+0.1767766952966368*alphal[16]*fUpwind_l[56]+0.1767766952966368*alphal[5]*fUpwind_l[55]+0.1767766952966368*alphal[20]*fUpwind_l[54]+0.1767766952966368*alphal[8]*fUpwind_l[52]+0.1767766952966368*alphal[13]*fUpwind_l[51]+0.1581138830084189*alphal[9]*fUpwind_l[31]+0.1581138830084189*alphal[4]*fUpwind_l[30]+0.1581138830084189*alphal[17]*fUpwind_l[29]+0.1581138830084189*alphal[10]*fUpwind_l[25]; 
  Ghat_l[63] = 0.1767766952966368*alphal[0]*fUpwind_l[63]+0.1767766952966368*alphal[1]*fUpwind_l[62]+0.1767766952966368*alphal[2]*fUpwind_l[61]+0.1767766952966368*alphal[3]*fUpwind_l[60]+0.1767766952966368*alphal[5]*fUpwind_l[59]+0.1767766952966368*alphal[6]*fUpwind_l[58]+0.1767766952966368*alphal[7]*fUpwind_l[57]+0.1767766952966368*alphal[8]*fUpwind_l[56]+0.1767766952966368*alphal[12]*fUpwind_l[55]+0.1767766952966368*alphal[13]*fUpwind_l[54]+0.1767766952966368*alphal[16]*fUpwind_l[52]+0.1767766952966368*alphal[20]*fUpwind_l[51]+0.1581138830084189*alphal[4]*fUpwind_l[31]+0.1581138830084189*alphal[9]*fUpwind_l[30]+0.1581138830084189*alphal[10]*fUpwind_l[29]+0.1581138830084189*alphal[17]*fUpwind_l[25]; 
  Ghat_l[64] = 0.1767766952966368*alphal[17]*fUpwind_l[76]+0.1767766952966368*alphal[16]*fUpwind_l[75]+0.1767766952966368*alphal[10]*fUpwind_l[73]+0.1767766952966368*alphal[9]*fUpwind_l[72]+0.1767766952966368*alphal[8]*fUpwind_l[71]+0.1767766952966368*alphal[7]*fUpwind_l[70]+0.1767766952966368*alphal[6]*fUpwind_l[69]+0.1767766952966368*alphal[4]*fUpwind_l[68]+0.1767766952966368*alphal[3]*fUpwind_l[67]+0.1767766952966368*alphal[2]*fUpwind_l[66]+0.1767766952966368*alphal[1]*fUpwind_l[65]+0.1767766952966368*alphal[0]*fUpwind_l[64]+0.1581138830084189*alphal[20]*fUpwind_l[20]+0.1581138830084189*alphal[13]*fUpwind_l[13]+0.1581138830084189*alphal[12]*fUpwind_l[12]+0.1581138830084189*alphal[5]*fUpwind_l[5]; 
  Ghat_l[65] = 0.1767766952966368*alphal[10]*fUpwind_l[76]+0.1767766952966368*alphal[8]*fUpwind_l[75]+0.1767766952966368*alphal[17]*fUpwind_l[73]+0.1767766952966368*alphal[4]*fUpwind_l[72]+0.1767766952966368*alphal[16]*fUpwind_l[71]+0.1767766952966368*alphal[3]*fUpwind_l[70]+0.1767766952966368*alphal[2]*fUpwind_l[69]+0.1767766952966368*alphal[9]*fUpwind_l[68]+0.1767766952966368*alphal[7]*fUpwind_l[67]+0.1767766952966368*alphal[6]*fUpwind_l[66]+0.1767766952966368*alphal[0]*fUpwind_l[65]+0.1767766952966368*alphal[1]*fUpwind_l[64]+0.1581138830084189*alphal[13]*fUpwind_l[20]+0.1581138830084189*fUpwind_l[13]*alphal[20]+0.1581138830084189*alphal[5]*fUpwind_l[12]+0.1581138830084189*fUpwind_l[5]*alphal[12]; 
  Ghat_l[66] = 0.1767766952966368*alphal[9]*fUpwind_l[76]+0.1767766952966368*alphal[7]*fUpwind_l[75]+0.1767766952966368*alphal[4]*fUpwind_l[73]+0.1767766952966368*alphal[17]*fUpwind_l[72]+0.1767766952966368*alphal[3]*fUpwind_l[71]+0.1767766952966368*alphal[16]*fUpwind_l[70]+0.1767766952966368*alphal[1]*fUpwind_l[69]+0.1767766952966368*alphal[10]*fUpwind_l[68]+0.1767766952966368*alphal[8]*fUpwind_l[67]+0.1767766952966368*alphal[0]*fUpwind_l[66]+0.1767766952966368*alphal[6]*fUpwind_l[65]+0.1767766952966368*alphal[2]*fUpwind_l[64]+0.1581138830084189*alphal[12]*fUpwind_l[20]+0.1581138830084189*fUpwind_l[12]*alphal[20]+0.1581138830084189*alphal[5]*fUpwind_l[13]+0.1581138830084189*fUpwind_l[5]*alphal[13]; 
  Ghat_l[67] = 0.1767766952966368*alphal[17]*fUpwind_l[79]+0.1767766952966368*alphal[10]*fUpwind_l[78]+0.1767766952966368*alphal[9]*fUpwind_l[77]+0.1767766952966368*alphal[6]*fUpwind_l[75]+0.1767766952966368*alphal[4]*fUpwind_l[74]+0.1767766952966368*alphal[2]*fUpwind_l[71]+0.1767766952966368*alphal[1]*fUpwind_l[70]+0.1767766952966368*alphal[16]*fUpwind_l[69]+0.1767766952966368*alphal[0]*fUpwind_l[67]+0.1767766952966368*alphal[8]*fUpwind_l[66]+0.1767766952966368*alphal[7]*fUpwind_l[65]+0.1767766952966368*alphal[3]*fUpwind_l[64]+0.1581138830084189*alphal[20]*fUpwind_l[27]+0.1581138830084189*alphal[13]*fUpwind_l[22]+0.1581138830084189*alphal[12]*fUpwind_l[21]+0.1581138830084189*alphal[5]*fUpwind_l[14]; 
  Ghat_l[68] = 0.1767766952966368*alphal[16]*fUpwind_l[79]+0.1767766952966368*alphal[8]*fUpwind_l[78]+0.1767766952966368*alphal[7]*fUpwind_l[77]+0.1767766952966368*alphal[6]*fUpwind_l[76]+0.1767766952966368*alphal[3]*fUpwind_l[74]+0.1767766952966368*alphal[2]*fUpwind_l[73]+0.1767766952966368*alphal[1]*fUpwind_l[72]+0.1767766952966368*alphal[17]*fUpwind_l[69]+0.1767766952966368*alphal[0]*fUpwind_l[68]+0.1767766952966368*alphal[10]*fUpwind_l[66]+0.1767766952966368*alphal[9]*fUpwind_l[65]+0.1767766952966368*alphal[4]*fUpwind_l[64]+0.1581138830084189*alphal[20]*fUpwind_l[28]+0.1581138830084189*alphal[13]*fUpwind_l[24]+0.1581138830084189*alphal[12]*fUpwind_l[23]+0.1581138830084189*alphal[5]*fUpwind_l[15]; 
  Ghat_l[69] = 0.1767766952966368*alphal[4]*fUpwind_l[76]+0.1767766952966368*alphal[3]*fUpwind_l[75]+0.1767766952966368*alphal[9]*fUpwind_l[73]+0.1767766952966368*alphal[10]*fUpwind_l[72]+0.1767766952966368*alphal[7]*fUpwind_l[71]+0.1767766952966368*alphal[8]*fUpwind_l[70]+0.1767766952966368*alphal[0]*fUpwind_l[69]+0.1767766952966368*alphal[17]*fUpwind_l[68]+0.1767766952966368*alphal[16]*fUpwind_l[67]+0.1767766952966368*alphal[1]*fUpwind_l[66]+0.1767766952966368*alphal[2]*fUpwind_l[65]+0.1767766952966368*alphal[6]*fUpwind_l[64]+0.1581138830084189*alphal[5]*fUpwind_l[20]+0.1581138830084189*fUpwind_l[5]*alphal[20]+0.1581138830084189*alphal[12]*fUpwind_l[13]+0.1581138830084189*fUpwind_l[12]*alphal[13]; 
  Ghat_l[70] = 0.1767766952966368*alphal[10]*fUpwind_l[79]+0.1767766952966368*alphal[17]*fUpwind_l[78]+0.1767766952966368*alphal[4]*fUpwind_l[77]+0.1767766952966368*alphal[2]*fUpwind_l[75]+0.1767766952966368*alphal[9]*fUpwind_l[74]+0.1767766952966368*alphal[6]*fUpwind_l[71]+0.1767766952966368*alphal[0]*fUpwind_l[70]+0.1767766952966368*alphal[8]*fUpwind_l[69]+0.1767766952966368*alphal[1]*fUpwind_l[67]+0.1767766952966368*alphal[16]*fUpwind_l[66]+0.1767766952966368*alphal[3]*fUpwind_l[65]+0.1767766952966368*alphal[7]*fUpwind_l[64]+0.1581138830084189*alphal[13]*fUpwind_l[27]+0.1581138830084189*alphal[20]*fUpwind_l[22]+0.1581138830084189*alphal[5]*fUpwind_l[21]+0.1581138830084189*alphal[12]*fUpwind_l[14]; 
  Ghat_l[71] = 0.1767766952966368*alphal[9]*fUpwind_l[79]+0.1767766952966368*alphal[4]*fUpwind_l[78]+0.1767766952966368*alphal[17]*fUpwind_l[77]+0.1767766952966368*alphal[1]*fUpwind_l[75]+0.1767766952966368*alphal[10]*fUpwind_l[74]+0.1767766952966368*alphal[0]*fUpwind_l[71]+0.1767766952966368*alphal[6]*fUpwind_l[70]+0.1767766952966368*alphal[7]*fUpwind_l[69]+0.1767766952966368*alphal[2]*fUpwind_l[67]+0.1767766952966368*alphal[3]*fUpwind_l[66]+0.1767766952966368*alphal[16]*fUpwind_l[65]+0.1767766952966368*alphal[8]*fUpwind_l[64]+0.1581138830084189*alphal[12]*fUpwind_l[27]+0.1581138830084189*alphal[5]*fUpwind_l[22]+0.1581138830084189*alphal[20]*fUpwind_l[21]+0.1581138830084189*alphal[13]*fUpwind_l[14]; 
  Ghat_l[72] = 0.1767766952966368*alphal[8]*fUpwind_l[79]+0.1767766952966368*alphal[16]*fUpwind_l[78]+0.1767766952966368*alphal[3]*fUpwind_l[77]+0.1767766952966368*alphal[2]*fUpwind_l[76]+0.1767766952966368*alphal[7]*fUpwind_l[74]+0.1767766952966368*alphal[6]*fUpwind_l[73]+0.1767766952966368*alphal[0]*fUpwind_l[72]+0.1767766952966368*alphal[10]*fUpwind_l[69]+0.1767766952966368*alphal[1]*fUpwind_l[68]+0.1767766952966368*alphal[17]*fUpwind_l[66]+0.1767766952966368*alphal[4]*fUpwind_l[65]+0.1767766952966368*alphal[9]*fUpwind_l[64]+0.1581138830084189*alphal[13]*fUpwind_l[28]+0.1581138830084189*alphal[20]*fUpwind_l[24]+0.1581138830084189*alphal[5]*fUpwind_l[23]+0.1581138830084189*alphal[12]*fUpwind_l[15]; 
  Ghat_l[73] = 0.1767766952966368*alphal[7]*fUpwind_l[79]+0.1767766952966368*alphal[3]*fUpwind_l[78]+0.1767766952966368*alphal[16]*fUpwind_l[77]+0.1767766952966368*alphal[1]*fUpwind_l[76]+0.1767766952966368*alphal[8]*fUpwind_l[74]+0.1767766952966368*alphal[0]*fUpwind_l[73]+0.1767766952966368*alphal[6]*fUpwind_l[72]+0.1767766952966368*alphal[9]*fUpwind_l[69]+0.1767766952966368*alphal[2]*fUpwind_l[68]+0.1767766952966368*alphal[4]*fUpwind_l[66]+0.1767766952966368*alphal[17]*fUpwind_l[65]+0.1767766952966368*alphal[10]*fUpwind_l[64]+0.1581138830084189*alphal[12]*fUpwind_l[28]+0.1581138830084189*alphal[5]*fUpwind_l[24]+0.1581138830084189*alphal[20]*fUpwind_l[23]+0.1581138830084189*alphal[13]*fUpwind_l[15]; 
  Ghat_l[74] = 0.1767766952966368*alphal[6]*fUpwind_l[79]+0.1767766952966368*alphal[2]*fUpwind_l[78]+0.1767766952966368*alphal[1]*fUpwind_l[77]+0.1767766952966368*alphal[16]*fUpwind_l[76]+0.1767766952966368*alphal[17]*fUpwind_l[75]+0.1767766952966368*alphal[0]*fUpwind_l[74]+0.1767766952966368*alphal[8]*fUpwind_l[73]+0.1767766952966368*alphal[7]*fUpwind_l[72]+0.1767766952966368*alphal[10]*fUpwind_l[71]+0.1767766952966368*alphal[9]*fUpwind_l[70]+0.1767766952966368*alphal[3]*fUpwind_l[68]+0.1767766952966368*alphal[4]*fUpwind_l[67]+0.1581138830084189*alphal[20]*fUpwind_l[31]+0.1581138830084189*alphal[13]*fUpwind_l[30]+0.1581138830084189*alphal[12]*fUpwind_l[29]+0.1581138830084189*alphal[5]*fUpwind_l[25]; 
  Ghat_l[75] = 0.1767766952966368*alphal[4]*fUpwind_l[79]+0.1767766952966368*alphal[9]*fUpwind_l[78]+0.1767766952966368*alphal[10]*fUpwind_l[77]+0.1767766952966368*alphal[0]*fUpwind_l[75]+0.1767766952966368*alphal[17]*fUpwind_l[74]+0.1767766952966368*alphal[1]*fUpwind_l[71]+0.1767766952966368*alphal[2]*fUpwind_l[70]+0.1767766952966368*alphal[3]*fUpwind_l[69]+0.1767766952966368*alphal[6]*fUpwind_l[67]+0.1767766952966368*alphal[7]*fUpwind_l[66]+0.1767766952966368*alphal[8]*fUpwind_l[65]+0.1767766952966368*alphal[16]*fUpwind_l[64]+0.1581138830084189*alphal[5]*fUpwind_l[27]+0.1581138830084189*alphal[12]*fUpwind_l[22]+0.1581138830084189*alphal[13]*fUpwind_l[21]+0.1581138830084189*fUpwind_l[14]*alphal[20]; 
  Ghat_l[76] = 0.1767766952966368*alphal[3]*fUpwind_l[79]+0.1767766952966368*alphal[7]*fUpwind_l[78]+0.1767766952966368*alphal[8]*fUpwind_l[77]+0.1767766952966368*alphal[0]*fUpwind_l[76]+0.1767766952966368*alphal[16]*fUpwind_l[74]+0.1767766952966368*alphal[1]*fUpwind_l[73]+0.1767766952966368*alphal[2]*fUpwind_l[72]+0.1767766952966368*alphal[4]*fUpwind_l[69]+0.1767766952966368*alphal[6]*fUpwind_l[68]+0.1767766952966368*alphal[9]*fUpwind_l[66]+0.1767766952966368*alphal[10]*fUpwind_l[65]+0.1767766952966368*alphal[17]*fUpwind_l[64]+0.1581138830084189*alphal[5]*fUpwind_l[28]+0.1581138830084189*alphal[12]*fUpwind_l[24]+0.1581138830084189*alphal[13]*fUpwind_l[23]+0.1581138830084189*fUpwind_l[15]*alphal[20]; 
  Ghat_l[77] = 0.1767766952966368*alphal[2]*fUpwind_l[79]+0.1767766952966368*alphal[6]*fUpwind_l[78]+0.1767766952966368*alphal[0]*fUpwind_l[77]+0.1767766952966368*alphal[8]*fUpwind_l[76]+0.1767766952966368*alphal[10]*fUpwind_l[75]+0.1767766952966368*alphal[1]*fUpwind_l[74]+0.1767766952966368*alphal[16]*fUpwind_l[73]+0.1767766952966368*alphal[3]*fUpwind_l[72]+0.1767766952966368*alphal[17]*fUpwind_l[71]+0.1767766952966368*alphal[4]*fUpwind_l[70]+0.1767766952966368*alphal[7]*fUpwind_l[68]+0.1767766952966368*alphal[9]*fUpwind_l[67]+0.1581138830084189*alphal[13]*fUpwind_l[31]+0.1581138830084189*alphal[20]*fUpwind_l[30]+0.1581138830084189*alphal[5]*fUpwind_l[29]+0.1581138830084189*alphal[12]*fUpwind_l[25]; 
  Ghat_l[78] = 0.1767766952966368*alphal[1]*fUpwind_l[79]+0.1767766952966368*alphal[0]*fUpwind_l[78]+0.1767766952966368*alphal[6]*fUpwind_l[77]+0.1767766952966368*alphal[7]*fUpwind_l[76]+0.1767766952966368*alphal[9]*fUpwind_l[75]+0.1767766952966368*alphal[2]*fUpwind_l[74]+0.1767766952966368*alphal[3]*fUpwind_l[73]+0.1767766952966368*alphal[16]*fUpwind_l[72]+0.1767766952966368*alphal[4]*fUpwind_l[71]+0.1767766952966368*alphal[17]*fUpwind_l[70]+0.1767766952966368*alphal[8]*fUpwind_l[68]+0.1767766952966368*alphal[10]*fUpwind_l[67]+0.1581138830084189*alphal[12]*fUpwind_l[31]+0.1581138830084189*alphal[5]*fUpwind_l[30]+0.1581138830084189*alphal[20]*fUpwind_l[29]+0.1581138830084189*alphal[13]*fUpwind_l[25]; 
  Ghat_l[79] = 0.1767766952966368*alphal[0]*fUpwind_l[79]+0.1767766952966368*alphal[1]*fUpwind_l[78]+0.1767766952966368*alphal[2]*fUpwind_l[77]+0.1767766952966368*alphal[3]*fUpwind_l[76]+0.1767766952966368*alphal[4]*fUpwind_l[75]+0.1767766952966368*alphal[6]*fUpwind_l[74]+0.1767766952966368*alphal[7]*fUpwind_l[73]+0.1767766952966368*alphal[8]*fUpwind_l[72]+0.1767766952966368*alphal[9]*fUpwind_l[71]+0.1767766952966368*alphal[10]*fUpwind_l[70]+0.1767766952966368*alphal[16]*fUpwind_l[68]+0.1767766952966368*alphal[17]*fUpwind_l[67]+0.1581138830084189*alphal[5]*fUpwind_l[31]+0.1581138830084189*alphal[12]*fUpwind_l[30]+0.1581138830084189*alphal[13]*fUpwind_l[29]+0.1581138830084189*alphal[20]*fUpwind_l[25]; 

  Ghat_r[0] = 0.1767766952966368*alphar[20]*fUpwind_r[20]+0.1767766952966368*alphar[17]*fUpwind_r[17]+0.1767766952966368*alphar[16]*fUpwind_r[16]+0.1767766952966368*alphar[13]*fUpwind_r[13]+0.1767766952966368*alphar[12]*fUpwind_r[12]+0.1767766952966368*alphar[10]*fUpwind_r[10]+0.1767766952966368*alphar[9]*fUpwind_r[9]+0.1767766952966368*alphar[8]*fUpwind_r[8]+0.1767766952966368*alphar[7]*fUpwind_r[7]+0.1767766952966368*alphar[6]*fUpwind_r[6]+0.1767766952966368*alphar[5]*fUpwind_r[5]+0.1767766952966368*alphar[4]*fUpwind_r[4]+0.1767766952966368*alphar[3]*fUpwind_r[3]+0.1767766952966368*alphar[2]*fUpwind_r[2]+0.1767766952966368*alphar[1]*fUpwind_r[1]+0.1767766952966368*alphar[0]*fUpwind_r[0]; 
  Ghat_r[1] = 0.1767766952966368*alphar[13]*fUpwind_r[20]+0.1767766952966368*fUpwind_r[13]*alphar[20]+0.1767766952966368*alphar[10]*fUpwind_r[17]+0.1767766952966368*fUpwind_r[10]*alphar[17]+0.1767766952966368*alphar[8]*fUpwind_r[16]+0.1767766952966368*fUpwind_r[8]*alphar[16]+0.1767766952966368*alphar[5]*fUpwind_r[12]+0.1767766952966368*fUpwind_r[5]*alphar[12]+0.1767766952966368*alphar[4]*fUpwind_r[9]+0.1767766952966368*fUpwind_r[4]*alphar[9]+0.1767766952966368*alphar[3]*fUpwind_r[7]+0.1767766952966368*fUpwind_r[3]*alphar[7]+0.1767766952966368*alphar[2]*fUpwind_r[6]+0.1767766952966368*fUpwind_r[2]*alphar[6]+0.1767766952966368*alphar[0]*fUpwind_r[1]+0.1767766952966368*fUpwind_r[0]*alphar[1]; 
  Ghat_r[2] = 0.1767766952966368*alphar[12]*fUpwind_r[20]+0.1767766952966368*fUpwind_r[12]*alphar[20]+0.1767766952966368*alphar[9]*fUpwind_r[17]+0.1767766952966368*fUpwind_r[9]*alphar[17]+0.1767766952966368*alphar[7]*fUpwind_r[16]+0.1767766952966368*fUpwind_r[7]*alphar[16]+0.1767766952966368*alphar[5]*fUpwind_r[13]+0.1767766952966368*fUpwind_r[5]*alphar[13]+0.1767766952966368*alphar[4]*fUpwind_r[10]+0.1767766952966368*fUpwind_r[4]*alphar[10]+0.1767766952966368*alphar[3]*fUpwind_r[8]+0.1767766952966368*fUpwind_r[3]*alphar[8]+0.1767766952966368*alphar[1]*fUpwind_r[6]+0.1767766952966368*fUpwind_r[1]*alphar[6]+0.1767766952966368*alphar[0]*fUpwind_r[2]+0.1767766952966368*fUpwind_r[0]*alphar[2]; 
  Ghat_r[3] = 0.1581138830084189*alphar[16]*fUpwind_r[37]+0.1581138830084189*alphar[8]*fUpwind_r[34]+0.1581138830084189*alphar[7]*fUpwind_r[33]+0.1581138830084189*alphar[3]*fUpwind_r[32]+0.1767766952966368*alphar[20]*fUpwind_r[27]+0.1767766952966368*alphar[17]*fUpwind_r[26]+0.1767766952966368*alphar[13]*fUpwind_r[22]+0.1767766952966368*alphar[12]*fUpwind_r[21]+0.1767766952966368*alphar[10]*fUpwind_r[19]+0.1767766952966368*alphar[9]*fUpwind_r[18]+0.1767766952966368*alphar[6]*fUpwind_r[16]+0.1767766952966368*fUpwind_r[6]*alphar[16]+0.1767766952966368*alphar[5]*fUpwind_r[14]+0.1767766952966368*alphar[4]*fUpwind_r[11]+0.1767766952966368*alphar[2]*fUpwind_r[8]+0.1767766952966368*fUpwind_r[2]*alphar[8]+0.1767766952966368*alphar[1]*fUpwind_r[7]+0.1767766952966368*fUpwind_r[1]*alphar[7]+0.1767766952966368*alphar[0]*fUpwind_r[3]+0.1767766952966368*fUpwind_r[0]*alphar[3]; 
  Ghat_r[4] = 0.1581138830084189*alphar[17]*fUpwind_r[53]+0.1581138830084189*alphar[10]*fUpwind_r[50]+0.1581138830084189*alphar[9]*fUpwind_r[49]+0.1581138830084189*alphar[4]*fUpwind_r[48]+0.1767766952966368*alphar[20]*fUpwind_r[28]+0.1767766952966368*alphar[16]*fUpwind_r[26]+0.1767766952966368*alphar[13]*fUpwind_r[24]+0.1767766952966368*alphar[12]*fUpwind_r[23]+0.1767766952966368*alphar[8]*fUpwind_r[19]+0.1767766952966368*alphar[7]*fUpwind_r[18]+0.1767766952966368*alphar[6]*fUpwind_r[17]+0.1767766952966368*fUpwind_r[6]*alphar[17]+0.1767766952966368*alphar[5]*fUpwind_r[15]+0.1767766952966368*alphar[3]*fUpwind_r[11]+0.1767766952966368*alphar[2]*fUpwind_r[10]+0.1767766952966368*fUpwind_r[2]*alphar[10]+0.1767766952966368*alphar[1]*fUpwind_r[9]+0.1767766952966368*fUpwind_r[1]*alphar[9]+0.1767766952966368*alphar[0]*fUpwind_r[4]+0.1767766952966368*fUpwind_r[0]*alphar[4]; 
  Ghat_r[5] = 0.1581138830084189*alphar[20]*fUpwind_r[69]+0.1581138830084189*alphar[13]*fUpwind_r[66]+0.1581138830084189*alphar[12]*fUpwind_r[65]+0.1581138830084189*alphar[5]*fUpwind_r[64]+0.1767766952966368*alphar[17]*fUpwind_r[28]+0.1767766952966368*alphar[16]*fUpwind_r[27]+0.1767766952966368*alphar[10]*fUpwind_r[24]+0.1767766952966368*alphar[9]*fUpwind_r[23]+0.1767766952966368*alphar[8]*fUpwind_r[22]+0.1767766952966368*alphar[7]*fUpwind_r[21]+0.1767766952966368*alphar[6]*fUpwind_r[20]+0.1767766952966368*fUpwind_r[6]*alphar[20]+0.1767766952966368*alphar[4]*fUpwind_r[15]+0.1767766952966368*alphar[3]*fUpwind_r[14]+0.1767766952966368*alphar[2]*fUpwind_r[13]+0.1767766952966368*fUpwind_r[2]*alphar[13]+0.1767766952966368*alphar[1]*fUpwind_r[12]+0.1767766952966368*fUpwind_r[1]*alphar[12]+0.1767766952966368*alphar[0]*fUpwind_r[5]+0.1767766952966368*fUpwind_r[0]*alphar[5]; 
  Ghat_r[6] = 0.1767766952966368*alphar[5]*fUpwind_r[20]+0.1767766952966368*fUpwind_r[5]*alphar[20]+0.1767766952966368*alphar[4]*fUpwind_r[17]+0.1767766952966368*fUpwind_r[4]*alphar[17]+0.1767766952966368*alphar[3]*fUpwind_r[16]+0.1767766952966368*fUpwind_r[3]*alphar[16]+0.1767766952966368*alphar[12]*fUpwind_r[13]+0.1767766952966368*fUpwind_r[12]*alphar[13]+0.1767766952966368*alphar[9]*fUpwind_r[10]+0.1767766952966368*fUpwind_r[9]*alphar[10]+0.1767766952966368*alphar[7]*fUpwind_r[8]+0.1767766952966368*fUpwind_r[7]*alphar[8]+0.1767766952966368*alphar[0]*fUpwind_r[6]+0.1767766952966368*fUpwind_r[0]*alphar[6]+0.1767766952966368*alphar[1]*fUpwind_r[2]+0.1767766952966368*fUpwind_r[1]*alphar[2]; 
  Ghat_r[7] = 0.1581138830084189*alphar[8]*fUpwind_r[37]+0.1581138830084189*alphar[16]*fUpwind_r[34]+0.1581138830084189*alphar[3]*fUpwind_r[33]+0.1581138830084189*alphar[7]*fUpwind_r[32]+0.1767766952966368*alphar[13]*fUpwind_r[27]+0.1767766952966368*alphar[10]*fUpwind_r[26]+0.1767766952966368*alphar[20]*fUpwind_r[22]+0.1767766952966368*alphar[5]*fUpwind_r[21]+0.1767766952966368*alphar[17]*fUpwind_r[19]+0.1767766952966368*alphar[4]*fUpwind_r[18]+0.1767766952966368*alphar[2]*fUpwind_r[16]+0.1767766952966368*fUpwind_r[2]*alphar[16]+0.1767766952966368*alphar[12]*fUpwind_r[14]+0.1767766952966368*alphar[9]*fUpwind_r[11]+0.1767766952966368*alphar[6]*fUpwind_r[8]+0.1767766952966368*fUpwind_r[6]*alphar[8]+0.1767766952966368*alphar[0]*fUpwind_r[7]+0.1767766952966368*fUpwind_r[0]*alphar[7]+0.1767766952966368*alphar[1]*fUpwind_r[3]+0.1767766952966368*fUpwind_r[1]*alphar[3]; 
  Ghat_r[8] = 0.1581138830084189*alphar[7]*fUpwind_r[37]+0.1581138830084189*alphar[3]*fUpwind_r[34]+0.1581138830084189*alphar[16]*fUpwind_r[33]+0.1581138830084189*alphar[8]*fUpwind_r[32]+0.1767766952966368*alphar[12]*fUpwind_r[27]+0.1767766952966368*alphar[9]*fUpwind_r[26]+0.1767766952966368*alphar[5]*fUpwind_r[22]+0.1767766952966368*alphar[20]*fUpwind_r[21]+0.1767766952966368*alphar[4]*fUpwind_r[19]+0.1767766952966368*alphar[17]*fUpwind_r[18]+0.1767766952966368*alphar[1]*fUpwind_r[16]+0.1767766952966368*fUpwind_r[1]*alphar[16]+0.1767766952966368*alphar[13]*fUpwind_r[14]+0.1767766952966368*alphar[10]*fUpwind_r[11]+0.1767766952966368*alphar[0]*fUpwind_r[8]+0.1767766952966368*fUpwind_r[0]*alphar[8]+0.1767766952966368*alphar[6]*fUpwind_r[7]+0.1767766952966368*fUpwind_r[6]*alphar[7]+0.1767766952966368*alphar[2]*fUpwind_r[3]+0.1767766952966368*fUpwind_r[2]*alphar[3]; 
  Ghat_r[9] = 0.1581138830084189*alphar[10]*fUpwind_r[53]+0.1581138830084189*alphar[17]*fUpwind_r[50]+0.1581138830084189*alphar[4]*fUpwind_r[49]+0.1581138830084189*alphar[9]*fUpwind_r[48]+0.1767766952966368*alphar[13]*fUpwind_r[28]+0.1767766952966368*alphar[8]*fUpwind_r[26]+0.1767766952966368*alphar[20]*fUpwind_r[24]+0.1767766952966368*alphar[5]*fUpwind_r[23]+0.1767766952966368*alphar[16]*fUpwind_r[19]+0.1767766952966368*alphar[3]*fUpwind_r[18]+0.1767766952966368*alphar[2]*fUpwind_r[17]+0.1767766952966368*fUpwind_r[2]*alphar[17]+0.1767766952966368*alphar[12]*fUpwind_r[15]+0.1767766952966368*alphar[7]*fUpwind_r[11]+0.1767766952966368*alphar[6]*fUpwind_r[10]+0.1767766952966368*fUpwind_r[6]*alphar[10]+0.1767766952966368*alphar[0]*fUpwind_r[9]+0.1767766952966368*fUpwind_r[0]*alphar[9]+0.1767766952966368*alphar[1]*fUpwind_r[4]+0.1767766952966368*fUpwind_r[1]*alphar[4]; 
  Ghat_r[10] = 0.1581138830084189*alphar[9]*fUpwind_r[53]+0.1581138830084189*alphar[4]*fUpwind_r[50]+0.1581138830084189*alphar[17]*fUpwind_r[49]+0.1581138830084189*alphar[10]*fUpwind_r[48]+0.1767766952966368*alphar[12]*fUpwind_r[28]+0.1767766952966368*alphar[7]*fUpwind_r[26]+0.1767766952966368*alphar[5]*fUpwind_r[24]+0.1767766952966368*alphar[20]*fUpwind_r[23]+0.1767766952966368*alphar[3]*fUpwind_r[19]+0.1767766952966368*alphar[16]*fUpwind_r[18]+0.1767766952966368*alphar[1]*fUpwind_r[17]+0.1767766952966368*fUpwind_r[1]*alphar[17]+0.1767766952966368*alphar[13]*fUpwind_r[15]+0.1767766952966368*alphar[8]*fUpwind_r[11]+0.1767766952966368*alphar[0]*fUpwind_r[10]+0.1767766952966368*fUpwind_r[0]*alphar[10]+0.1767766952966368*alphar[6]*fUpwind_r[9]+0.1767766952966368*fUpwind_r[6]*alphar[9]+0.1767766952966368*alphar[2]*fUpwind_r[4]+0.1767766952966368*fUpwind_r[2]*alphar[4]; 
  Ghat_r[11] = 0.1581138830084189*alphar[17]*fUpwind_r[59]+0.1581138830084189*alphar[10]*fUpwind_r[55]+0.1581138830084189*alphar[9]*fUpwind_r[54]+0.1581138830084189*alphar[4]*fUpwind_r[51]+0.1581138830084189*alphar[16]*fUpwind_r[43]+0.1581138830084189*alphar[8]*fUpwind_r[39]+0.1581138830084189*alphar[7]*fUpwind_r[38]+0.1581138830084189*alphar[3]*fUpwind_r[35]+0.1767766952966368*alphar[20]*fUpwind_r[31]+0.1767766952966368*alphar[13]*fUpwind_r[30]+0.1767766952966368*alphar[12]*fUpwind_r[29]+0.1767766952966368*alphar[6]*fUpwind_r[26]+0.1767766952966368*alphar[5]*fUpwind_r[25]+0.1767766952966368*alphar[2]*fUpwind_r[19]+0.1767766952966368*alphar[1]*fUpwind_r[18]+0.1767766952966368*alphar[16]*fUpwind_r[17]+0.1767766952966368*fUpwind_r[16]*alphar[17]+0.1767766952966368*alphar[0]*fUpwind_r[11]+0.1767766952966368*alphar[8]*fUpwind_r[10]+0.1767766952966368*fUpwind_r[8]*alphar[10]+0.1767766952966368*alphar[7]*fUpwind_r[9]+0.1767766952966368*fUpwind_r[7]*alphar[9]+0.1767766952966368*alphar[3]*fUpwind_r[4]+0.1767766952966368*fUpwind_r[3]*alphar[4]; 
  Ghat_r[12] = 0.1581138830084189*alphar[13]*fUpwind_r[69]+0.1581138830084189*alphar[20]*fUpwind_r[66]+0.1581138830084189*alphar[5]*fUpwind_r[65]+0.1581138830084189*alphar[12]*fUpwind_r[64]+0.1767766952966368*alphar[10]*fUpwind_r[28]+0.1767766952966368*alphar[8]*fUpwind_r[27]+0.1767766952966368*alphar[17]*fUpwind_r[24]+0.1767766952966368*alphar[4]*fUpwind_r[23]+0.1767766952966368*alphar[16]*fUpwind_r[22]+0.1767766952966368*alphar[3]*fUpwind_r[21]+0.1767766952966368*alphar[2]*fUpwind_r[20]+0.1767766952966368*fUpwind_r[2]*alphar[20]+0.1767766952966368*alphar[9]*fUpwind_r[15]+0.1767766952966368*alphar[7]*fUpwind_r[14]+0.1767766952966368*alphar[6]*fUpwind_r[13]+0.1767766952966368*fUpwind_r[6]*alphar[13]+0.1767766952966368*alphar[0]*fUpwind_r[12]+0.1767766952966368*fUpwind_r[0]*alphar[12]+0.1767766952966368*alphar[1]*fUpwind_r[5]+0.1767766952966368*fUpwind_r[1]*alphar[5]; 
  Ghat_r[13] = 0.1581138830084189*alphar[12]*fUpwind_r[69]+0.1581138830084189*alphar[5]*fUpwind_r[66]+0.1581138830084189*alphar[20]*fUpwind_r[65]+0.1581138830084189*alphar[13]*fUpwind_r[64]+0.1767766952966368*alphar[9]*fUpwind_r[28]+0.1767766952966368*alphar[7]*fUpwind_r[27]+0.1767766952966368*alphar[4]*fUpwind_r[24]+0.1767766952966368*alphar[17]*fUpwind_r[23]+0.1767766952966368*alphar[3]*fUpwind_r[22]+0.1767766952966368*alphar[16]*fUpwind_r[21]+0.1767766952966368*alphar[1]*fUpwind_r[20]+0.1767766952966368*fUpwind_r[1]*alphar[20]+0.1767766952966368*alphar[10]*fUpwind_r[15]+0.1767766952966368*alphar[8]*fUpwind_r[14]+0.1767766952966368*alphar[0]*fUpwind_r[13]+0.1767766952966368*fUpwind_r[0]*alphar[13]+0.1767766952966368*alphar[6]*fUpwind_r[12]+0.1767766952966368*fUpwind_r[6]*alphar[12]+0.1767766952966368*alphar[2]*fUpwind_r[5]+0.1767766952966368*fUpwind_r[2]*alphar[5]; 
  Ghat_r[14] = 0.1581138830084189*alphar[20]*fUpwind_r[75]+0.1581138830084189*alphar[13]*fUpwind_r[71]+0.1581138830084189*alphar[12]*fUpwind_r[70]+0.1581138830084189*alphar[5]*fUpwind_r[67]+0.1581138830084189*alphar[16]*fUpwind_r[44]+0.1581138830084189*alphar[8]*fUpwind_r[41]+0.1581138830084189*alphar[7]*fUpwind_r[40]+0.1581138830084189*alphar[3]*fUpwind_r[36]+0.1767766952966368*alphar[17]*fUpwind_r[31]+0.1767766952966368*alphar[10]*fUpwind_r[30]+0.1767766952966368*alphar[9]*fUpwind_r[29]+0.1767766952966368*alphar[6]*fUpwind_r[27]+0.1767766952966368*alphar[4]*fUpwind_r[25]+0.1767766952966368*alphar[2]*fUpwind_r[22]+0.1767766952966368*alphar[1]*fUpwind_r[21]+0.1767766952966368*alphar[16]*fUpwind_r[20]+0.1767766952966368*fUpwind_r[16]*alphar[20]+0.1767766952966368*alphar[0]*fUpwind_r[14]+0.1767766952966368*alphar[8]*fUpwind_r[13]+0.1767766952966368*fUpwind_r[8]*alphar[13]+0.1767766952966368*alphar[7]*fUpwind_r[12]+0.1767766952966368*fUpwind_r[7]*alphar[12]+0.1767766952966368*alphar[3]*fUpwind_r[5]+0.1767766952966368*fUpwind_r[3]*alphar[5]; 
  Ghat_r[15] = 0.1581138830084189*alphar[20]*fUpwind_r[76]+0.1581138830084189*alphar[13]*fUpwind_r[73]+0.1581138830084189*alphar[12]*fUpwind_r[72]+0.1581138830084189*alphar[5]*fUpwind_r[68]+0.1581138830084189*alphar[17]*fUpwind_r[60]+0.1581138830084189*alphar[10]*fUpwind_r[57]+0.1581138830084189*alphar[9]*fUpwind_r[56]+0.1581138830084189*alphar[4]*fUpwind_r[52]+0.1767766952966368*alphar[16]*fUpwind_r[31]+0.1767766952966368*alphar[8]*fUpwind_r[30]+0.1767766952966368*alphar[7]*fUpwind_r[29]+0.1767766952966368*alphar[6]*fUpwind_r[28]+0.1767766952966368*alphar[3]*fUpwind_r[25]+0.1767766952966368*alphar[2]*fUpwind_r[24]+0.1767766952966368*alphar[1]*fUpwind_r[23]+0.1767766952966368*alphar[17]*fUpwind_r[20]+0.1767766952966368*fUpwind_r[17]*alphar[20]+0.1767766952966368*alphar[0]*fUpwind_r[15]+0.1767766952966368*alphar[10]*fUpwind_r[13]+0.1767766952966368*fUpwind_r[10]*alphar[13]+0.1767766952966368*alphar[9]*fUpwind_r[12]+0.1767766952966368*fUpwind_r[9]*alphar[12]+0.1767766952966368*alphar[4]*fUpwind_r[5]+0.1767766952966368*fUpwind_r[4]*alphar[5]; 
  Ghat_r[16] = 0.1581138830084189*alphar[3]*fUpwind_r[37]+0.1581138830084189*alphar[7]*fUpwind_r[34]+0.1581138830084189*alphar[8]*fUpwind_r[33]+0.1581138830084189*alphar[16]*fUpwind_r[32]+0.1767766952966368*alphar[5]*fUpwind_r[27]+0.1767766952966368*alphar[4]*fUpwind_r[26]+0.1767766952966368*alphar[12]*fUpwind_r[22]+0.1767766952966368*alphar[13]*fUpwind_r[21]+0.1767766952966368*fUpwind_r[14]*alphar[20]+0.1767766952966368*alphar[9]*fUpwind_r[19]+0.1767766952966368*alphar[10]*fUpwind_r[18]+0.1767766952966368*fUpwind_r[11]*alphar[17]+0.1767766952966368*alphar[0]*fUpwind_r[16]+0.1767766952966368*fUpwind_r[0]*alphar[16]+0.1767766952966368*alphar[1]*fUpwind_r[8]+0.1767766952966368*fUpwind_r[1]*alphar[8]+0.1767766952966368*alphar[2]*fUpwind_r[7]+0.1767766952966368*fUpwind_r[2]*alphar[7]+0.1767766952966368*alphar[3]*fUpwind_r[6]+0.1767766952966368*fUpwind_r[3]*alphar[6]; 
  Ghat_r[17] = 0.1581138830084189*alphar[4]*fUpwind_r[53]+0.1581138830084189*alphar[9]*fUpwind_r[50]+0.1581138830084189*alphar[10]*fUpwind_r[49]+0.1581138830084189*alphar[17]*fUpwind_r[48]+0.1767766952966368*alphar[5]*fUpwind_r[28]+0.1767766952966368*alphar[3]*fUpwind_r[26]+0.1767766952966368*alphar[12]*fUpwind_r[24]+0.1767766952966368*alphar[13]*fUpwind_r[23]+0.1767766952966368*fUpwind_r[15]*alphar[20]+0.1767766952966368*alphar[7]*fUpwind_r[19]+0.1767766952966368*alphar[8]*fUpwind_r[18]+0.1767766952966368*alphar[0]*fUpwind_r[17]+0.1767766952966368*fUpwind_r[0]*alphar[17]+0.1767766952966368*fUpwind_r[11]*alphar[16]+0.1767766952966368*alphar[1]*fUpwind_r[10]+0.1767766952966368*fUpwind_r[1]*alphar[10]+0.1767766952966368*alphar[2]*fUpwind_r[9]+0.1767766952966368*fUpwind_r[2]*alphar[9]+0.1767766952966368*alphar[4]*fUpwind_r[6]+0.1767766952966368*fUpwind_r[4]*alphar[6]; 
  Ghat_r[18] = 0.1581138830084189*alphar[10]*fUpwind_r[59]+0.1581138830084189*alphar[17]*fUpwind_r[55]+0.1581138830084189*alphar[4]*fUpwind_r[54]+0.1581138830084189*alphar[9]*fUpwind_r[51]+0.1581138830084189*alphar[8]*fUpwind_r[43]+0.1581138830084189*alphar[16]*fUpwind_r[39]+0.1581138830084189*alphar[3]*fUpwind_r[38]+0.1581138830084189*alphar[7]*fUpwind_r[35]+0.1767766952966368*alphar[13]*fUpwind_r[31]+0.1767766952966368*alphar[20]*fUpwind_r[30]+0.1767766952966368*alphar[5]*fUpwind_r[29]+0.1767766952966368*alphar[2]*fUpwind_r[26]+0.1767766952966368*alphar[12]*fUpwind_r[25]+0.1767766952966368*alphar[6]*fUpwind_r[19]+0.1767766952966368*alphar[0]*fUpwind_r[18]+0.1767766952966368*alphar[8]*fUpwind_r[17]+0.1767766952966368*fUpwind_r[8]*alphar[17]+0.1767766952966368*alphar[10]*fUpwind_r[16]+0.1767766952966368*fUpwind_r[10]*alphar[16]+0.1767766952966368*alphar[1]*fUpwind_r[11]+0.1767766952966368*alphar[3]*fUpwind_r[9]+0.1767766952966368*fUpwind_r[3]*alphar[9]+0.1767766952966368*alphar[4]*fUpwind_r[7]+0.1767766952966368*fUpwind_r[4]*alphar[7]; 
  Ghat_r[19] = 0.1581138830084189*alphar[9]*fUpwind_r[59]+0.1581138830084189*alphar[4]*fUpwind_r[55]+0.1581138830084189*alphar[17]*fUpwind_r[54]+0.1581138830084189*alphar[10]*fUpwind_r[51]+0.1581138830084189*alphar[7]*fUpwind_r[43]+0.1581138830084189*alphar[3]*fUpwind_r[39]+0.1581138830084189*alphar[16]*fUpwind_r[38]+0.1581138830084189*alphar[8]*fUpwind_r[35]+0.1767766952966368*alphar[12]*fUpwind_r[31]+0.1767766952966368*alphar[5]*fUpwind_r[30]+0.1767766952966368*alphar[20]*fUpwind_r[29]+0.1767766952966368*alphar[1]*fUpwind_r[26]+0.1767766952966368*alphar[13]*fUpwind_r[25]+0.1767766952966368*alphar[0]*fUpwind_r[19]+0.1767766952966368*alphar[6]*fUpwind_r[18]+0.1767766952966368*alphar[7]*fUpwind_r[17]+0.1767766952966368*fUpwind_r[7]*alphar[17]+0.1767766952966368*alphar[9]*fUpwind_r[16]+0.1767766952966368*fUpwind_r[9]*alphar[16]+0.1767766952966368*alphar[2]*fUpwind_r[11]+0.1767766952966368*alphar[3]*fUpwind_r[10]+0.1767766952966368*fUpwind_r[3]*alphar[10]+0.1767766952966368*alphar[4]*fUpwind_r[8]+0.1767766952966368*fUpwind_r[4]*alphar[8]; 
  Ghat_r[20] = 0.1581138830084189*alphar[5]*fUpwind_r[69]+0.1581138830084189*alphar[12]*fUpwind_r[66]+0.1581138830084189*alphar[13]*fUpwind_r[65]+0.1581138830084189*alphar[20]*fUpwind_r[64]+0.1767766952966368*alphar[4]*fUpwind_r[28]+0.1767766952966368*alphar[3]*fUpwind_r[27]+0.1767766952966368*alphar[9]*fUpwind_r[24]+0.1767766952966368*alphar[10]*fUpwind_r[23]+0.1767766952966368*alphar[7]*fUpwind_r[22]+0.1767766952966368*alphar[8]*fUpwind_r[21]+0.1767766952966368*alphar[0]*fUpwind_r[20]+0.1767766952966368*fUpwind_r[0]*alphar[20]+0.1767766952966368*fUpwind_r[15]*alphar[17]+0.1767766952966368*fUpwind_r[14]*alphar[16]+0.1767766952966368*alphar[1]*fUpwind_r[13]+0.1767766952966368*fUpwind_r[1]*alphar[13]+0.1767766952966368*alphar[2]*fUpwind_r[12]+0.1767766952966368*fUpwind_r[2]*alphar[12]+0.1767766952966368*alphar[5]*fUpwind_r[6]+0.1767766952966368*fUpwind_r[5]*alphar[6]; 
  Ghat_r[21] = 0.1581138830084189*alphar[13]*fUpwind_r[75]+0.1581138830084189*alphar[20]*fUpwind_r[71]+0.1581138830084189*alphar[5]*fUpwind_r[70]+0.1581138830084189*alphar[12]*fUpwind_r[67]+0.1581138830084189*alphar[8]*fUpwind_r[44]+0.1581138830084189*alphar[16]*fUpwind_r[41]+0.1581138830084189*alphar[3]*fUpwind_r[40]+0.1581138830084189*alphar[7]*fUpwind_r[36]+0.1767766952966368*alphar[10]*fUpwind_r[31]+0.1767766952966368*alphar[17]*fUpwind_r[30]+0.1767766952966368*alphar[4]*fUpwind_r[29]+0.1767766952966368*alphar[2]*fUpwind_r[27]+0.1767766952966368*alphar[9]*fUpwind_r[25]+0.1767766952966368*alphar[6]*fUpwind_r[22]+0.1767766952966368*alphar[0]*fUpwind_r[21]+0.1767766952966368*alphar[8]*fUpwind_r[20]+0.1767766952966368*fUpwind_r[8]*alphar[20]+0.1767766952966368*alphar[13]*fUpwind_r[16]+0.1767766952966368*fUpwind_r[13]*alphar[16]+0.1767766952966368*alphar[1]*fUpwind_r[14]+0.1767766952966368*alphar[3]*fUpwind_r[12]+0.1767766952966368*fUpwind_r[3]*alphar[12]+0.1767766952966368*alphar[5]*fUpwind_r[7]+0.1767766952966368*fUpwind_r[5]*alphar[7]; 
  Ghat_r[22] = 0.1581138830084189*alphar[12]*fUpwind_r[75]+0.1581138830084189*alphar[5]*fUpwind_r[71]+0.1581138830084189*alphar[20]*fUpwind_r[70]+0.1581138830084189*alphar[13]*fUpwind_r[67]+0.1581138830084189*alphar[7]*fUpwind_r[44]+0.1581138830084189*alphar[3]*fUpwind_r[41]+0.1581138830084189*alphar[16]*fUpwind_r[40]+0.1581138830084189*alphar[8]*fUpwind_r[36]+0.1767766952966368*alphar[9]*fUpwind_r[31]+0.1767766952966368*alphar[4]*fUpwind_r[30]+0.1767766952966368*alphar[17]*fUpwind_r[29]+0.1767766952966368*alphar[1]*fUpwind_r[27]+0.1767766952966368*alphar[10]*fUpwind_r[25]+0.1767766952966368*alphar[0]*fUpwind_r[22]+0.1767766952966368*alphar[6]*fUpwind_r[21]+0.1767766952966368*alphar[7]*fUpwind_r[20]+0.1767766952966368*fUpwind_r[7]*alphar[20]+0.1767766952966368*alphar[12]*fUpwind_r[16]+0.1767766952966368*fUpwind_r[12]*alphar[16]+0.1767766952966368*alphar[2]*fUpwind_r[14]+0.1767766952966368*alphar[3]*fUpwind_r[13]+0.1767766952966368*fUpwind_r[3]*alphar[13]+0.1767766952966368*alphar[5]*fUpwind_r[8]+0.1767766952966368*fUpwind_r[5]*alphar[8]; 
  Ghat_r[23] = 0.1581138830084189*alphar[13]*fUpwind_r[76]+0.1581138830084189*alphar[20]*fUpwind_r[73]+0.1581138830084189*alphar[5]*fUpwind_r[72]+0.1581138830084189*alphar[12]*fUpwind_r[68]+0.1581138830084189*alphar[10]*fUpwind_r[60]+0.1581138830084189*alphar[17]*fUpwind_r[57]+0.1581138830084189*alphar[4]*fUpwind_r[56]+0.1581138830084189*alphar[9]*fUpwind_r[52]+0.1767766952966368*alphar[8]*fUpwind_r[31]+0.1767766952966368*alphar[16]*fUpwind_r[30]+0.1767766952966368*alphar[3]*fUpwind_r[29]+0.1767766952966368*alphar[2]*fUpwind_r[28]+0.1767766952966368*alphar[7]*fUpwind_r[25]+0.1767766952966368*alphar[6]*fUpwind_r[24]+0.1767766952966368*alphar[0]*fUpwind_r[23]+0.1767766952966368*alphar[10]*fUpwind_r[20]+0.1767766952966368*fUpwind_r[10]*alphar[20]+0.1767766952966368*alphar[13]*fUpwind_r[17]+0.1767766952966368*fUpwind_r[13]*alphar[17]+0.1767766952966368*alphar[1]*fUpwind_r[15]+0.1767766952966368*alphar[4]*fUpwind_r[12]+0.1767766952966368*fUpwind_r[4]*alphar[12]+0.1767766952966368*alphar[5]*fUpwind_r[9]+0.1767766952966368*fUpwind_r[5]*alphar[9]; 
  Ghat_r[24] = 0.1581138830084189*alphar[12]*fUpwind_r[76]+0.1581138830084189*alphar[5]*fUpwind_r[73]+0.1581138830084189*alphar[20]*fUpwind_r[72]+0.1581138830084189*alphar[13]*fUpwind_r[68]+0.1581138830084189*alphar[9]*fUpwind_r[60]+0.1581138830084189*alphar[4]*fUpwind_r[57]+0.1581138830084189*alphar[17]*fUpwind_r[56]+0.1581138830084189*alphar[10]*fUpwind_r[52]+0.1767766952966368*alphar[7]*fUpwind_r[31]+0.1767766952966368*alphar[3]*fUpwind_r[30]+0.1767766952966368*alphar[16]*fUpwind_r[29]+0.1767766952966368*alphar[1]*fUpwind_r[28]+0.1767766952966368*alphar[8]*fUpwind_r[25]+0.1767766952966368*alphar[0]*fUpwind_r[24]+0.1767766952966368*alphar[6]*fUpwind_r[23]+0.1767766952966368*alphar[9]*fUpwind_r[20]+0.1767766952966368*fUpwind_r[9]*alphar[20]+0.1767766952966368*alphar[12]*fUpwind_r[17]+0.1767766952966368*fUpwind_r[12]*alphar[17]+0.1767766952966368*alphar[2]*fUpwind_r[15]+0.1767766952966368*alphar[4]*fUpwind_r[13]+0.1767766952966368*fUpwind_r[4]*alphar[13]+0.1767766952966368*alphar[5]*fUpwind_r[10]+0.1767766952966368*fUpwind_r[5]*alphar[10]; 
  Ghat_r[25] = 0.1581138830084189*alphar[20]*fUpwind_r[79]+0.1581138830084189*alphar[13]*fUpwind_r[78]+0.1581138830084189*alphar[12]*fUpwind_r[77]+0.1581138830084189*alphar[5]*fUpwind_r[74]+0.1581138830084189*alphar[17]*fUpwind_r[63]+0.1581138830084189*alphar[10]*fUpwind_r[62]+0.1581138830084189*alphar[9]*fUpwind_r[61]+0.1581138830084189*alphar[4]*fUpwind_r[58]+0.1581138830084189*alphar[16]*fUpwind_r[47]+0.1581138830084189*alphar[8]*fUpwind_r[46]+0.1581138830084189*alphar[7]*fUpwind_r[45]+0.1581138830084189*alphar[3]*fUpwind_r[42]+0.1767766952966368*alphar[6]*fUpwind_r[31]+0.1767766952966368*alphar[2]*fUpwind_r[30]+0.1767766952966368*alphar[1]*fUpwind_r[29]+0.1767766952966368*alphar[16]*fUpwind_r[28]+0.1767766952966368*alphar[17]*fUpwind_r[27]+0.1767766952966368*alphar[20]*fUpwind_r[26]+0.1767766952966368*alphar[0]*fUpwind_r[25]+0.1767766952966368*alphar[8]*fUpwind_r[24]+0.1767766952966368*alphar[7]*fUpwind_r[23]+0.1767766952966368*alphar[10]*fUpwind_r[22]+0.1767766952966368*alphar[9]*fUpwind_r[21]+0.1767766952966368*alphar[13]*fUpwind_r[19]+0.1767766952966368*alphar[12]*fUpwind_r[18]+0.1767766952966368*alphar[3]*fUpwind_r[15]+0.1767766952966368*alphar[4]*fUpwind_r[14]+0.1767766952966368*alphar[5]*fUpwind_r[11]; 
  Ghat_r[26] = 0.1581138830084189*alphar[4]*fUpwind_r[59]+0.1581138830084189*alphar[9]*fUpwind_r[55]+0.1581138830084189*alphar[10]*fUpwind_r[54]+0.1581138830084189*alphar[17]*fUpwind_r[51]+0.1581138830084189*alphar[3]*fUpwind_r[43]+0.1581138830084189*alphar[7]*fUpwind_r[39]+0.1581138830084189*alphar[8]*fUpwind_r[38]+0.1581138830084189*alphar[16]*fUpwind_r[35]+0.1767766952966368*alphar[5]*fUpwind_r[31]+0.1767766952966368*alphar[12]*fUpwind_r[30]+0.1767766952966368*alphar[13]*fUpwind_r[29]+0.1767766952966368*alphar[0]*fUpwind_r[26]+0.1767766952966368*alphar[20]*fUpwind_r[25]+0.1767766952966368*alphar[1]*fUpwind_r[19]+0.1767766952966368*alphar[2]*fUpwind_r[18]+0.1767766952966368*alphar[3]*fUpwind_r[17]+0.1767766952966368*fUpwind_r[3]*alphar[17]+0.1767766952966368*alphar[4]*fUpwind_r[16]+0.1767766952966368*fUpwind_r[4]*alphar[16]+0.1767766952966368*alphar[6]*fUpwind_r[11]+0.1767766952966368*alphar[7]*fUpwind_r[10]+0.1767766952966368*fUpwind_r[7]*alphar[10]+0.1767766952966368*alphar[8]*fUpwind_r[9]+0.1767766952966368*fUpwind_r[8]*alphar[9]; 
  Ghat_r[27] = 0.1581138830084189*alphar[5]*fUpwind_r[75]+0.1581138830084189*alphar[12]*fUpwind_r[71]+0.1581138830084189*alphar[13]*fUpwind_r[70]+0.1581138830084189*alphar[20]*fUpwind_r[67]+0.1581138830084189*alphar[3]*fUpwind_r[44]+0.1581138830084189*alphar[7]*fUpwind_r[41]+0.1581138830084189*alphar[8]*fUpwind_r[40]+0.1581138830084189*alphar[16]*fUpwind_r[36]+0.1767766952966368*alphar[4]*fUpwind_r[31]+0.1767766952966368*alphar[9]*fUpwind_r[30]+0.1767766952966368*alphar[10]*fUpwind_r[29]+0.1767766952966368*alphar[0]*fUpwind_r[27]+0.1767766952966368*alphar[17]*fUpwind_r[25]+0.1767766952966368*alphar[1]*fUpwind_r[22]+0.1767766952966368*alphar[2]*fUpwind_r[21]+0.1767766952966368*alphar[3]*fUpwind_r[20]+0.1767766952966368*fUpwind_r[3]*alphar[20]+0.1767766952966368*alphar[5]*fUpwind_r[16]+0.1767766952966368*fUpwind_r[5]*alphar[16]+0.1767766952966368*alphar[6]*fUpwind_r[14]+0.1767766952966368*alphar[7]*fUpwind_r[13]+0.1767766952966368*fUpwind_r[7]*alphar[13]+0.1767766952966368*alphar[8]*fUpwind_r[12]+0.1767766952966368*fUpwind_r[8]*alphar[12]; 
  Ghat_r[28] = 0.1581138830084189*alphar[5]*fUpwind_r[76]+0.1581138830084189*alphar[12]*fUpwind_r[73]+0.1581138830084189*alphar[13]*fUpwind_r[72]+0.1581138830084189*alphar[20]*fUpwind_r[68]+0.1581138830084189*alphar[4]*fUpwind_r[60]+0.1581138830084189*alphar[9]*fUpwind_r[57]+0.1581138830084189*alphar[10]*fUpwind_r[56]+0.1581138830084189*alphar[17]*fUpwind_r[52]+0.1767766952966368*alphar[3]*fUpwind_r[31]+0.1767766952966368*alphar[7]*fUpwind_r[30]+0.1767766952966368*alphar[8]*fUpwind_r[29]+0.1767766952966368*alphar[0]*fUpwind_r[28]+0.1767766952966368*alphar[16]*fUpwind_r[25]+0.1767766952966368*alphar[1]*fUpwind_r[24]+0.1767766952966368*alphar[2]*fUpwind_r[23]+0.1767766952966368*alphar[4]*fUpwind_r[20]+0.1767766952966368*fUpwind_r[4]*alphar[20]+0.1767766952966368*alphar[5]*fUpwind_r[17]+0.1767766952966368*fUpwind_r[5]*alphar[17]+0.1767766952966368*alphar[6]*fUpwind_r[15]+0.1767766952966368*alphar[9]*fUpwind_r[13]+0.1767766952966368*fUpwind_r[9]*alphar[13]+0.1767766952966368*alphar[10]*fUpwind_r[12]+0.1767766952966368*fUpwind_r[10]*alphar[12]; 
  Ghat_r[29] = 0.1581138830084189*alphar[13]*fUpwind_r[79]+0.1581138830084189*alphar[20]*fUpwind_r[78]+0.1581138830084189*alphar[5]*fUpwind_r[77]+0.1581138830084189*alphar[12]*fUpwind_r[74]+0.1581138830084189*alphar[10]*fUpwind_r[63]+0.1581138830084189*alphar[17]*fUpwind_r[62]+0.1581138830084189*alphar[4]*fUpwind_r[61]+0.1581138830084189*alphar[9]*fUpwind_r[58]+0.1581138830084189*alphar[8]*fUpwind_r[47]+0.1581138830084189*alphar[16]*fUpwind_r[46]+0.1581138830084189*alphar[3]*fUpwind_r[45]+0.1581138830084189*alphar[7]*fUpwind_r[42]+0.1767766952966368*alphar[2]*fUpwind_r[31]+0.1767766952966368*alphar[6]*fUpwind_r[30]+0.1767766952966368*alphar[0]*fUpwind_r[29]+0.1767766952966368*alphar[8]*fUpwind_r[28]+0.1767766952966368*alphar[10]*fUpwind_r[27]+0.1767766952966368*alphar[13]*fUpwind_r[26]+0.1767766952966368*alphar[1]*fUpwind_r[25]+0.1767766952966368*alphar[16]*fUpwind_r[24]+0.1767766952966368*alphar[3]*fUpwind_r[23]+0.1767766952966368*alphar[17]*fUpwind_r[22]+0.1767766952966368*alphar[4]*fUpwind_r[21]+0.1767766952966368*fUpwind_r[19]*alphar[20]+0.1767766952966368*alphar[5]*fUpwind_r[18]+0.1767766952966368*alphar[7]*fUpwind_r[15]+0.1767766952966368*alphar[9]*fUpwind_r[14]+0.1767766952966368*fUpwind_r[11]*alphar[12]; 
  Ghat_r[30] = 0.1581138830084189*alphar[12]*fUpwind_r[79]+0.1581138830084189*alphar[5]*fUpwind_r[78]+0.1581138830084189*alphar[20]*fUpwind_r[77]+0.1581138830084189*alphar[13]*fUpwind_r[74]+0.1581138830084189*alphar[9]*fUpwind_r[63]+0.1581138830084189*alphar[4]*fUpwind_r[62]+0.1581138830084189*alphar[17]*fUpwind_r[61]+0.1581138830084189*alphar[10]*fUpwind_r[58]+0.1581138830084189*alphar[7]*fUpwind_r[47]+0.1581138830084189*alphar[3]*fUpwind_r[46]+0.1581138830084189*alphar[16]*fUpwind_r[45]+0.1581138830084189*alphar[8]*fUpwind_r[42]+0.1767766952966368*alphar[1]*fUpwind_r[31]+0.1767766952966368*alphar[0]*fUpwind_r[30]+0.1767766952966368*alphar[6]*fUpwind_r[29]+0.1767766952966368*alphar[7]*fUpwind_r[28]+0.1767766952966368*alphar[9]*fUpwind_r[27]+0.1767766952966368*alphar[12]*fUpwind_r[26]+0.1767766952966368*alphar[2]*fUpwind_r[25]+0.1767766952966368*alphar[3]*fUpwind_r[24]+0.1767766952966368*alphar[16]*fUpwind_r[23]+0.1767766952966368*alphar[4]*fUpwind_r[22]+0.1767766952966368*alphar[17]*fUpwind_r[21]+0.1767766952966368*fUpwind_r[18]*alphar[20]+0.1767766952966368*alphar[5]*fUpwind_r[19]+0.1767766952966368*alphar[8]*fUpwind_r[15]+0.1767766952966368*alphar[10]*fUpwind_r[14]+0.1767766952966368*fUpwind_r[11]*alphar[13]; 
  Ghat_r[31] = 0.1581138830084189*alphar[5]*fUpwind_r[79]+0.1581138830084189*alphar[12]*fUpwind_r[78]+0.1581138830084189*alphar[13]*fUpwind_r[77]+0.1581138830084189*alphar[20]*fUpwind_r[74]+0.1581138830084189*alphar[4]*fUpwind_r[63]+0.1581138830084189*alphar[9]*fUpwind_r[62]+0.1581138830084189*alphar[10]*fUpwind_r[61]+0.1581138830084189*alphar[17]*fUpwind_r[58]+0.1581138830084189*alphar[3]*fUpwind_r[47]+0.1581138830084189*alphar[7]*fUpwind_r[46]+0.1581138830084189*alphar[8]*fUpwind_r[45]+0.1581138830084189*alphar[16]*fUpwind_r[42]+0.1767766952966368*alphar[0]*fUpwind_r[31]+0.1767766952966368*alphar[1]*fUpwind_r[30]+0.1767766952966368*alphar[2]*fUpwind_r[29]+0.1767766952966368*alphar[3]*fUpwind_r[28]+0.1767766952966368*alphar[4]*fUpwind_r[27]+0.1767766952966368*alphar[5]*fUpwind_r[26]+0.1767766952966368*alphar[6]*fUpwind_r[25]+0.1767766952966368*alphar[7]*fUpwind_r[24]+0.1767766952966368*alphar[8]*fUpwind_r[23]+0.1767766952966368*alphar[9]*fUpwind_r[22]+0.1767766952966368*alphar[10]*fUpwind_r[21]+0.1767766952966368*fUpwind_r[11]*alphar[20]+0.1767766952966368*alphar[12]*fUpwind_r[19]+0.1767766952966368*alphar[13]*fUpwind_r[18]+0.1767766952966368*fUpwind_r[14]*alphar[17]+0.1767766952966368*fUpwind_r[15]*alphar[16]; 
  Ghat_r[32] = 0.1767766952966368*alphar[20]*fUpwind_r[44]+0.1767766952966368*alphar[17]*fUpwind_r[43]+0.1767766952966368*alphar[13]*fUpwind_r[41]+0.1767766952966368*alphar[12]*fUpwind_r[40]+0.1767766952966368*alphar[10]*fUpwind_r[39]+0.1767766952966368*alphar[9]*fUpwind_r[38]+0.1767766952966368*alphar[6]*fUpwind_r[37]+0.1767766952966368*alphar[5]*fUpwind_r[36]+0.1767766952966368*alphar[4]*fUpwind_r[35]+0.1767766952966368*alphar[2]*fUpwind_r[34]+0.1767766952966368*alphar[1]*fUpwind_r[33]+0.1767766952966368*alphar[0]*fUpwind_r[32]+0.1581138830084189*alphar[16]*fUpwind_r[16]+0.1581138830084189*alphar[8]*fUpwind_r[8]+0.1581138830084189*alphar[7]*fUpwind_r[7]+0.1581138830084189*alphar[3]*fUpwind_r[3]; 
  Ghat_r[33] = 0.1767766952966368*alphar[13]*fUpwind_r[44]+0.1767766952966368*alphar[10]*fUpwind_r[43]+0.1767766952966368*alphar[20]*fUpwind_r[41]+0.1767766952966368*alphar[5]*fUpwind_r[40]+0.1767766952966368*alphar[17]*fUpwind_r[39]+0.1767766952966368*alphar[4]*fUpwind_r[38]+0.1767766952966368*alphar[2]*fUpwind_r[37]+0.1767766952966368*alphar[12]*fUpwind_r[36]+0.1767766952966368*alphar[9]*fUpwind_r[35]+0.1767766952966368*alphar[6]*fUpwind_r[34]+0.1767766952966368*alphar[0]*fUpwind_r[33]+0.1767766952966368*alphar[1]*fUpwind_r[32]+0.1581138830084189*alphar[8]*fUpwind_r[16]+0.1581138830084189*fUpwind_r[8]*alphar[16]+0.1581138830084189*alphar[3]*fUpwind_r[7]+0.1581138830084189*fUpwind_r[3]*alphar[7]; 
  Ghat_r[34] = 0.1767766952966368*alphar[12]*fUpwind_r[44]+0.1767766952966368*alphar[9]*fUpwind_r[43]+0.1767766952966368*alphar[5]*fUpwind_r[41]+0.1767766952966368*alphar[20]*fUpwind_r[40]+0.1767766952966368*alphar[4]*fUpwind_r[39]+0.1767766952966368*alphar[17]*fUpwind_r[38]+0.1767766952966368*alphar[1]*fUpwind_r[37]+0.1767766952966368*alphar[13]*fUpwind_r[36]+0.1767766952966368*alphar[10]*fUpwind_r[35]+0.1767766952966368*alphar[0]*fUpwind_r[34]+0.1767766952966368*alphar[6]*fUpwind_r[33]+0.1767766952966368*alphar[2]*fUpwind_r[32]+0.1581138830084189*alphar[7]*fUpwind_r[16]+0.1581138830084189*fUpwind_r[7]*alphar[16]+0.1581138830084189*alphar[3]*fUpwind_r[8]+0.1581138830084189*fUpwind_r[3]*alphar[8]; 
  Ghat_r[35] = 0.1767766952966368*alphar[20]*fUpwind_r[47]+0.1767766952966368*alphar[13]*fUpwind_r[46]+0.1767766952966368*alphar[12]*fUpwind_r[45]+0.1767766952966368*alphar[6]*fUpwind_r[43]+0.1767766952966368*alphar[5]*fUpwind_r[42]+0.1767766952966368*alphar[2]*fUpwind_r[39]+0.1767766952966368*alphar[1]*fUpwind_r[38]+0.1767766952966368*alphar[17]*fUpwind_r[37]+0.1767766952966368*alphar[0]*fUpwind_r[35]+0.1767766952966368*alphar[10]*fUpwind_r[34]+0.1767766952966368*alphar[9]*fUpwind_r[33]+0.1767766952966368*alphar[4]*fUpwind_r[32]+0.1581138830084189*alphar[16]*fUpwind_r[26]+0.1581138830084189*alphar[8]*fUpwind_r[19]+0.1581138830084189*alphar[7]*fUpwind_r[18]+0.1581138830084189*alphar[3]*fUpwind_r[11]; 
  Ghat_r[36] = 0.1767766952966368*alphar[17]*fUpwind_r[47]+0.1767766952966368*alphar[10]*fUpwind_r[46]+0.1767766952966368*alphar[9]*fUpwind_r[45]+0.1767766952966368*alphar[6]*fUpwind_r[44]+0.1767766952966368*alphar[4]*fUpwind_r[42]+0.1767766952966368*alphar[2]*fUpwind_r[41]+0.1767766952966368*alphar[1]*fUpwind_r[40]+0.1767766952966368*alphar[20]*fUpwind_r[37]+0.1767766952966368*alphar[0]*fUpwind_r[36]+0.1767766952966368*alphar[13]*fUpwind_r[34]+0.1767766952966368*alphar[12]*fUpwind_r[33]+0.1767766952966368*alphar[5]*fUpwind_r[32]+0.1581138830084189*alphar[16]*fUpwind_r[27]+0.1581138830084189*alphar[8]*fUpwind_r[22]+0.1581138830084189*alphar[7]*fUpwind_r[21]+0.1581138830084189*alphar[3]*fUpwind_r[14]; 
  Ghat_r[37] = 0.1767766952966368*alphar[5]*fUpwind_r[44]+0.1767766952966368*alphar[4]*fUpwind_r[43]+0.1767766952966368*alphar[12]*fUpwind_r[41]+0.1767766952966368*alphar[13]*fUpwind_r[40]+0.1767766952966368*alphar[9]*fUpwind_r[39]+0.1767766952966368*alphar[10]*fUpwind_r[38]+0.1767766952966368*alphar[0]*fUpwind_r[37]+0.1767766952966368*alphar[20]*fUpwind_r[36]+0.1767766952966368*alphar[17]*fUpwind_r[35]+0.1767766952966368*alphar[1]*fUpwind_r[34]+0.1767766952966368*alphar[2]*fUpwind_r[33]+0.1767766952966368*alphar[6]*fUpwind_r[32]+0.1581138830084189*alphar[3]*fUpwind_r[16]+0.1581138830084189*fUpwind_r[3]*alphar[16]+0.1581138830084189*alphar[7]*fUpwind_r[8]+0.1581138830084189*fUpwind_r[7]*alphar[8]; 
  Ghat_r[38] = 0.1767766952966368*alphar[13]*fUpwind_r[47]+0.1767766952966368*alphar[20]*fUpwind_r[46]+0.1767766952966368*alphar[5]*fUpwind_r[45]+0.1767766952966368*alphar[2]*fUpwind_r[43]+0.1767766952966368*alphar[12]*fUpwind_r[42]+0.1767766952966368*alphar[6]*fUpwind_r[39]+0.1767766952966368*alphar[0]*fUpwind_r[38]+0.1767766952966368*alphar[10]*fUpwind_r[37]+0.1767766952966368*alphar[1]*fUpwind_r[35]+0.1767766952966368*alphar[17]*fUpwind_r[34]+0.1767766952966368*alphar[4]*fUpwind_r[33]+0.1767766952966368*alphar[9]*fUpwind_r[32]+0.1581138830084189*alphar[8]*fUpwind_r[26]+0.1581138830084189*alphar[16]*fUpwind_r[19]+0.1581138830084189*alphar[3]*fUpwind_r[18]+0.1581138830084189*alphar[7]*fUpwind_r[11]; 
  Ghat_r[39] = 0.1767766952966368*alphar[12]*fUpwind_r[47]+0.1767766952966368*alphar[5]*fUpwind_r[46]+0.1767766952966368*alphar[20]*fUpwind_r[45]+0.1767766952966368*alphar[1]*fUpwind_r[43]+0.1767766952966368*alphar[13]*fUpwind_r[42]+0.1767766952966368*alphar[0]*fUpwind_r[39]+0.1767766952966368*alphar[6]*fUpwind_r[38]+0.1767766952966368*alphar[9]*fUpwind_r[37]+0.1767766952966368*alphar[2]*fUpwind_r[35]+0.1767766952966368*alphar[4]*fUpwind_r[34]+0.1767766952966368*alphar[17]*fUpwind_r[33]+0.1767766952966368*alphar[10]*fUpwind_r[32]+0.1581138830084189*alphar[7]*fUpwind_r[26]+0.1581138830084189*alphar[3]*fUpwind_r[19]+0.1581138830084189*alphar[16]*fUpwind_r[18]+0.1581138830084189*alphar[8]*fUpwind_r[11]; 
  Ghat_r[40] = 0.1767766952966368*alphar[10]*fUpwind_r[47]+0.1767766952966368*alphar[17]*fUpwind_r[46]+0.1767766952966368*alphar[4]*fUpwind_r[45]+0.1767766952966368*alphar[2]*fUpwind_r[44]+0.1767766952966368*alphar[9]*fUpwind_r[42]+0.1767766952966368*alphar[6]*fUpwind_r[41]+0.1767766952966368*alphar[0]*fUpwind_r[40]+0.1767766952966368*alphar[13]*fUpwind_r[37]+0.1767766952966368*alphar[1]*fUpwind_r[36]+0.1767766952966368*alphar[20]*fUpwind_r[34]+0.1767766952966368*alphar[5]*fUpwind_r[33]+0.1767766952966368*alphar[12]*fUpwind_r[32]+0.1581138830084189*alphar[8]*fUpwind_r[27]+0.1581138830084189*alphar[16]*fUpwind_r[22]+0.1581138830084189*alphar[3]*fUpwind_r[21]+0.1581138830084189*alphar[7]*fUpwind_r[14]; 
  Ghat_r[41] = 0.1767766952966368*alphar[9]*fUpwind_r[47]+0.1767766952966368*alphar[4]*fUpwind_r[46]+0.1767766952966368*alphar[17]*fUpwind_r[45]+0.1767766952966368*alphar[1]*fUpwind_r[44]+0.1767766952966368*alphar[10]*fUpwind_r[42]+0.1767766952966368*alphar[0]*fUpwind_r[41]+0.1767766952966368*alphar[6]*fUpwind_r[40]+0.1767766952966368*alphar[12]*fUpwind_r[37]+0.1767766952966368*alphar[2]*fUpwind_r[36]+0.1767766952966368*alphar[5]*fUpwind_r[34]+0.1767766952966368*alphar[20]*fUpwind_r[33]+0.1767766952966368*alphar[13]*fUpwind_r[32]+0.1581138830084189*alphar[7]*fUpwind_r[27]+0.1581138830084189*alphar[3]*fUpwind_r[22]+0.1581138830084189*alphar[16]*fUpwind_r[21]+0.1581138830084189*alphar[8]*fUpwind_r[14]; 
  Ghat_r[42] = 0.1767766952966368*alphar[6]*fUpwind_r[47]+0.1767766952966368*alphar[2]*fUpwind_r[46]+0.1767766952966368*alphar[1]*fUpwind_r[45]+0.1767766952966368*alphar[17]*fUpwind_r[44]+0.1767766952966368*alphar[20]*fUpwind_r[43]+0.1767766952966368*alphar[0]*fUpwind_r[42]+0.1767766952966368*alphar[10]*fUpwind_r[41]+0.1767766952966368*alphar[9]*fUpwind_r[40]+0.1767766952966368*alphar[13]*fUpwind_r[39]+0.1767766952966368*alphar[12]*fUpwind_r[38]+0.1767766952966368*alphar[4]*fUpwind_r[36]+0.1767766952966368*alphar[5]*fUpwind_r[35]+0.1581138830084189*alphar[16]*fUpwind_r[31]+0.1581138830084189*alphar[8]*fUpwind_r[30]+0.1581138830084189*alphar[7]*fUpwind_r[29]+0.1581138830084189*alphar[3]*fUpwind_r[25]; 
  Ghat_r[43] = 0.1767766952966368*alphar[5]*fUpwind_r[47]+0.1767766952966368*alphar[12]*fUpwind_r[46]+0.1767766952966368*alphar[13]*fUpwind_r[45]+0.1767766952966368*alphar[0]*fUpwind_r[43]+0.1767766952966368*alphar[20]*fUpwind_r[42]+0.1767766952966368*alphar[1]*fUpwind_r[39]+0.1767766952966368*alphar[2]*fUpwind_r[38]+0.1767766952966368*alphar[4]*fUpwind_r[37]+0.1767766952966368*alphar[6]*fUpwind_r[35]+0.1767766952966368*alphar[9]*fUpwind_r[34]+0.1767766952966368*alphar[10]*fUpwind_r[33]+0.1767766952966368*alphar[17]*fUpwind_r[32]+0.1581138830084189*alphar[3]*fUpwind_r[26]+0.1581138830084189*alphar[7]*fUpwind_r[19]+0.1581138830084189*alphar[8]*fUpwind_r[18]+0.1581138830084189*fUpwind_r[11]*alphar[16]; 
  Ghat_r[44] = 0.1767766952966368*alphar[4]*fUpwind_r[47]+0.1767766952966368*alphar[9]*fUpwind_r[46]+0.1767766952966368*alphar[10]*fUpwind_r[45]+0.1767766952966368*alphar[0]*fUpwind_r[44]+0.1767766952966368*alphar[17]*fUpwind_r[42]+0.1767766952966368*alphar[1]*fUpwind_r[41]+0.1767766952966368*alphar[2]*fUpwind_r[40]+0.1767766952966368*alphar[5]*fUpwind_r[37]+0.1767766952966368*alphar[6]*fUpwind_r[36]+0.1767766952966368*alphar[12]*fUpwind_r[34]+0.1767766952966368*alphar[13]*fUpwind_r[33]+0.1767766952966368*alphar[20]*fUpwind_r[32]+0.1581138830084189*alphar[3]*fUpwind_r[27]+0.1581138830084189*alphar[7]*fUpwind_r[22]+0.1581138830084189*alphar[8]*fUpwind_r[21]+0.1581138830084189*fUpwind_r[14]*alphar[16]; 
  Ghat_r[45] = 0.1767766952966368*alphar[2]*fUpwind_r[47]+0.1767766952966368*alphar[6]*fUpwind_r[46]+0.1767766952966368*alphar[0]*fUpwind_r[45]+0.1767766952966368*alphar[10]*fUpwind_r[44]+0.1767766952966368*alphar[13]*fUpwind_r[43]+0.1767766952966368*alphar[1]*fUpwind_r[42]+0.1767766952966368*alphar[17]*fUpwind_r[41]+0.1767766952966368*alphar[4]*fUpwind_r[40]+0.1767766952966368*alphar[20]*fUpwind_r[39]+0.1767766952966368*alphar[5]*fUpwind_r[38]+0.1767766952966368*alphar[9]*fUpwind_r[36]+0.1767766952966368*alphar[12]*fUpwind_r[35]+0.1581138830084189*alphar[8]*fUpwind_r[31]+0.1581138830084189*alphar[16]*fUpwind_r[30]+0.1581138830084189*alphar[3]*fUpwind_r[29]+0.1581138830084189*alphar[7]*fUpwind_r[25]; 
  Ghat_r[46] = 0.1767766952966368*alphar[1]*fUpwind_r[47]+0.1767766952966368*alphar[0]*fUpwind_r[46]+0.1767766952966368*alphar[6]*fUpwind_r[45]+0.1767766952966368*alphar[9]*fUpwind_r[44]+0.1767766952966368*alphar[12]*fUpwind_r[43]+0.1767766952966368*alphar[2]*fUpwind_r[42]+0.1767766952966368*alphar[4]*fUpwind_r[41]+0.1767766952966368*alphar[17]*fUpwind_r[40]+0.1767766952966368*alphar[5]*fUpwind_r[39]+0.1767766952966368*alphar[20]*fUpwind_r[38]+0.1767766952966368*alphar[10]*fUpwind_r[36]+0.1767766952966368*alphar[13]*fUpwind_r[35]+0.1581138830084189*alphar[7]*fUpwind_r[31]+0.1581138830084189*alphar[3]*fUpwind_r[30]+0.1581138830084189*alphar[16]*fUpwind_r[29]+0.1581138830084189*alphar[8]*fUpwind_r[25]; 
  Ghat_r[47] = 0.1767766952966368*alphar[0]*fUpwind_r[47]+0.1767766952966368*alphar[1]*fUpwind_r[46]+0.1767766952966368*alphar[2]*fUpwind_r[45]+0.1767766952966368*alphar[4]*fUpwind_r[44]+0.1767766952966368*alphar[5]*fUpwind_r[43]+0.1767766952966368*alphar[6]*fUpwind_r[42]+0.1767766952966368*alphar[9]*fUpwind_r[41]+0.1767766952966368*alphar[10]*fUpwind_r[40]+0.1767766952966368*alphar[12]*fUpwind_r[39]+0.1767766952966368*alphar[13]*fUpwind_r[38]+0.1767766952966368*alphar[17]*fUpwind_r[36]+0.1767766952966368*alphar[20]*fUpwind_r[35]+0.1581138830084189*alphar[3]*fUpwind_r[31]+0.1581138830084189*alphar[7]*fUpwind_r[30]+0.1581138830084189*alphar[8]*fUpwind_r[29]+0.1581138830084189*alphar[16]*fUpwind_r[25]; 
  Ghat_r[48] = 0.1767766952966368*alphar[20]*fUpwind_r[60]+0.1767766952966368*alphar[16]*fUpwind_r[59]+0.1767766952966368*alphar[13]*fUpwind_r[57]+0.1767766952966368*alphar[12]*fUpwind_r[56]+0.1767766952966368*alphar[8]*fUpwind_r[55]+0.1767766952966368*alphar[7]*fUpwind_r[54]+0.1767766952966368*alphar[6]*fUpwind_r[53]+0.1767766952966368*alphar[5]*fUpwind_r[52]+0.1767766952966368*alphar[3]*fUpwind_r[51]+0.1767766952966368*alphar[2]*fUpwind_r[50]+0.1767766952966368*alphar[1]*fUpwind_r[49]+0.1767766952966368*alphar[0]*fUpwind_r[48]+0.1581138830084189*alphar[17]*fUpwind_r[17]+0.1581138830084189*alphar[10]*fUpwind_r[10]+0.1581138830084189*alphar[9]*fUpwind_r[9]+0.1581138830084189*alphar[4]*fUpwind_r[4]; 
  Ghat_r[49] = 0.1767766952966368*alphar[13]*fUpwind_r[60]+0.1767766952966368*alphar[8]*fUpwind_r[59]+0.1767766952966368*alphar[20]*fUpwind_r[57]+0.1767766952966368*alphar[5]*fUpwind_r[56]+0.1767766952966368*alphar[16]*fUpwind_r[55]+0.1767766952966368*alphar[3]*fUpwind_r[54]+0.1767766952966368*alphar[2]*fUpwind_r[53]+0.1767766952966368*alphar[12]*fUpwind_r[52]+0.1767766952966368*alphar[7]*fUpwind_r[51]+0.1767766952966368*alphar[6]*fUpwind_r[50]+0.1767766952966368*alphar[0]*fUpwind_r[49]+0.1767766952966368*alphar[1]*fUpwind_r[48]+0.1581138830084189*alphar[10]*fUpwind_r[17]+0.1581138830084189*fUpwind_r[10]*alphar[17]+0.1581138830084189*alphar[4]*fUpwind_r[9]+0.1581138830084189*fUpwind_r[4]*alphar[9]; 
  Ghat_r[50] = 0.1767766952966368*alphar[12]*fUpwind_r[60]+0.1767766952966368*alphar[7]*fUpwind_r[59]+0.1767766952966368*alphar[5]*fUpwind_r[57]+0.1767766952966368*alphar[20]*fUpwind_r[56]+0.1767766952966368*alphar[3]*fUpwind_r[55]+0.1767766952966368*alphar[16]*fUpwind_r[54]+0.1767766952966368*alphar[1]*fUpwind_r[53]+0.1767766952966368*alphar[13]*fUpwind_r[52]+0.1767766952966368*alphar[8]*fUpwind_r[51]+0.1767766952966368*alphar[0]*fUpwind_r[50]+0.1767766952966368*alphar[6]*fUpwind_r[49]+0.1767766952966368*alphar[2]*fUpwind_r[48]+0.1581138830084189*alphar[9]*fUpwind_r[17]+0.1581138830084189*fUpwind_r[9]*alphar[17]+0.1581138830084189*alphar[4]*fUpwind_r[10]+0.1581138830084189*fUpwind_r[4]*alphar[10]; 
  Ghat_r[51] = 0.1767766952966368*alphar[20]*fUpwind_r[63]+0.1767766952966368*alphar[13]*fUpwind_r[62]+0.1767766952966368*alphar[12]*fUpwind_r[61]+0.1767766952966368*alphar[6]*fUpwind_r[59]+0.1767766952966368*alphar[5]*fUpwind_r[58]+0.1767766952966368*alphar[2]*fUpwind_r[55]+0.1767766952966368*alphar[1]*fUpwind_r[54]+0.1767766952966368*alphar[16]*fUpwind_r[53]+0.1767766952966368*alphar[0]*fUpwind_r[51]+0.1767766952966368*alphar[8]*fUpwind_r[50]+0.1767766952966368*alphar[7]*fUpwind_r[49]+0.1767766952966368*alphar[3]*fUpwind_r[48]+0.1581138830084189*alphar[17]*fUpwind_r[26]+0.1581138830084189*alphar[10]*fUpwind_r[19]+0.1581138830084189*alphar[9]*fUpwind_r[18]+0.1581138830084189*alphar[4]*fUpwind_r[11]; 
  Ghat_r[52] = 0.1767766952966368*alphar[16]*fUpwind_r[63]+0.1767766952966368*alphar[8]*fUpwind_r[62]+0.1767766952966368*alphar[7]*fUpwind_r[61]+0.1767766952966368*alphar[6]*fUpwind_r[60]+0.1767766952966368*alphar[3]*fUpwind_r[58]+0.1767766952966368*alphar[2]*fUpwind_r[57]+0.1767766952966368*alphar[1]*fUpwind_r[56]+0.1767766952966368*alphar[20]*fUpwind_r[53]+0.1767766952966368*alphar[0]*fUpwind_r[52]+0.1767766952966368*alphar[13]*fUpwind_r[50]+0.1767766952966368*alphar[12]*fUpwind_r[49]+0.1767766952966368*alphar[5]*fUpwind_r[48]+0.1581138830084189*alphar[17]*fUpwind_r[28]+0.1581138830084189*alphar[10]*fUpwind_r[24]+0.1581138830084189*alphar[9]*fUpwind_r[23]+0.1581138830084189*alphar[4]*fUpwind_r[15]; 
  Ghat_r[53] = 0.1767766952966368*alphar[5]*fUpwind_r[60]+0.1767766952966368*alphar[3]*fUpwind_r[59]+0.1767766952966368*alphar[12]*fUpwind_r[57]+0.1767766952966368*alphar[13]*fUpwind_r[56]+0.1767766952966368*alphar[7]*fUpwind_r[55]+0.1767766952966368*alphar[8]*fUpwind_r[54]+0.1767766952966368*alphar[0]*fUpwind_r[53]+0.1767766952966368*alphar[20]*fUpwind_r[52]+0.1767766952966368*alphar[16]*fUpwind_r[51]+0.1767766952966368*alphar[1]*fUpwind_r[50]+0.1767766952966368*alphar[2]*fUpwind_r[49]+0.1767766952966368*alphar[6]*fUpwind_r[48]+0.1581138830084189*alphar[4]*fUpwind_r[17]+0.1581138830084189*fUpwind_r[4]*alphar[17]+0.1581138830084189*alphar[9]*fUpwind_r[10]+0.1581138830084189*fUpwind_r[9]*alphar[10]; 
  Ghat_r[54] = 0.1767766952966368*alphar[13]*fUpwind_r[63]+0.1767766952966368*alphar[20]*fUpwind_r[62]+0.1767766952966368*alphar[5]*fUpwind_r[61]+0.1767766952966368*alphar[2]*fUpwind_r[59]+0.1767766952966368*alphar[12]*fUpwind_r[58]+0.1767766952966368*alphar[6]*fUpwind_r[55]+0.1767766952966368*alphar[0]*fUpwind_r[54]+0.1767766952966368*alphar[8]*fUpwind_r[53]+0.1767766952966368*alphar[1]*fUpwind_r[51]+0.1767766952966368*alphar[16]*fUpwind_r[50]+0.1767766952966368*alphar[3]*fUpwind_r[49]+0.1767766952966368*alphar[7]*fUpwind_r[48]+0.1581138830084189*alphar[10]*fUpwind_r[26]+0.1581138830084189*alphar[17]*fUpwind_r[19]+0.1581138830084189*alphar[4]*fUpwind_r[18]+0.1581138830084189*alphar[9]*fUpwind_r[11]; 
  Ghat_r[55] = 0.1767766952966368*alphar[12]*fUpwind_r[63]+0.1767766952966368*alphar[5]*fUpwind_r[62]+0.1767766952966368*alphar[20]*fUpwind_r[61]+0.1767766952966368*alphar[1]*fUpwind_r[59]+0.1767766952966368*alphar[13]*fUpwind_r[58]+0.1767766952966368*alphar[0]*fUpwind_r[55]+0.1767766952966368*alphar[6]*fUpwind_r[54]+0.1767766952966368*alphar[7]*fUpwind_r[53]+0.1767766952966368*alphar[2]*fUpwind_r[51]+0.1767766952966368*alphar[3]*fUpwind_r[50]+0.1767766952966368*alphar[16]*fUpwind_r[49]+0.1767766952966368*alphar[8]*fUpwind_r[48]+0.1581138830084189*alphar[9]*fUpwind_r[26]+0.1581138830084189*alphar[4]*fUpwind_r[19]+0.1581138830084189*alphar[17]*fUpwind_r[18]+0.1581138830084189*alphar[10]*fUpwind_r[11]; 
  Ghat_r[56] = 0.1767766952966368*alphar[8]*fUpwind_r[63]+0.1767766952966368*alphar[16]*fUpwind_r[62]+0.1767766952966368*alphar[3]*fUpwind_r[61]+0.1767766952966368*alphar[2]*fUpwind_r[60]+0.1767766952966368*alphar[7]*fUpwind_r[58]+0.1767766952966368*alphar[6]*fUpwind_r[57]+0.1767766952966368*alphar[0]*fUpwind_r[56]+0.1767766952966368*alphar[13]*fUpwind_r[53]+0.1767766952966368*alphar[1]*fUpwind_r[52]+0.1767766952966368*alphar[20]*fUpwind_r[50]+0.1767766952966368*alphar[5]*fUpwind_r[49]+0.1767766952966368*alphar[12]*fUpwind_r[48]+0.1581138830084189*alphar[10]*fUpwind_r[28]+0.1581138830084189*alphar[17]*fUpwind_r[24]+0.1581138830084189*alphar[4]*fUpwind_r[23]+0.1581138830084189*alphar[9]*fUpwind_r[15]; 
  Ghat_r[57] = 0.1767766952966368*alphar[7]*fUpwind_r[63]+0.1767766952966368*alphar[3]*fUpwind_r[62]+0.1767766952966368*alphar[16]*fUpwind_r[61]+0.1767766952966368*alphar[1]*fUpwind_r[60]+0.1767766952966368*alphar[8]*fUpwind_r[58]+0.1767766952966368*alphar[0]*fUpwind_r[57]+0.1767766952966368*alphar[6]*fUpwind_r[56]+0.1767766952966368*alphar[12]*fUpwind_r[53]+0.1767766952966368*alphar[2]*fUpwind_r[52]+0.1767766952966368*alphar[5]*fUpwind_r[50]+0.1767766952966368*alphar[20]*fUpwind_r[49]+0.1767766952966368*alphar[13]*fUpwind_r[48]+0.1581138830084189*alphar[9]*fUpwind_r[28]+0.1581138830084189*alphar[4]*fUpwind_r[24]+0.1581138830084189*alphar[17]*fUpwind_r[23]+0.1581138830084189*alphar[10]*fUpwind_r[15]; 
  Ghat_r[58] = 0.1767766952966368*alphar[6]*fUpwind_r[63]+0.1767766952966368*alphar[2]*fUpwind_r[62]+0.1767766952966368*alphar[1]*fUpwind_r[61]+0.1767766952966368*alphar[16]*fUpwind_r[60]+0.1767766952966368*alphar[20]*fUpwind_r[59]+0.1767766952966368*alphar[0]*fUpwind_r[58]+0.1767766952966368*alphar[8]*fUpwind_r[57]+0.1767766952966368*alphar[7]*fUpwind_r[56]+0.1767766952966368*alphar[13]*fUpwind_r[55]+0.1767766952966368*alphar[12]*fUpwind_r[54]+0.1767766952966368*alphar[3]*fUpwind_r[52]+0.1767766952966368*alphar[5]*fUpwind_r[51]+0.1581138830084189*alphar[17]*fUpwind_r[31]+0.1581138830084189*alphar[10]*fUpwind_r[30]+0.1581138830084189*alphar[9]*fUpwind_r[29]+0.1581138830084189*alphar[4]*fUpwind_r[25]; 
  Ghat_r[59] = 0.1767766952966368*alphar[5]*fUpwind_r[63]+0.1767766952966368*alphar[12]*fUpwind_r[62]+0.1767766952966368*alphar[13]*fUpwind_r[61]+0.1767766952966368*alphar[0]*fUpwind_r[59]+0.1767766952966368*alphar[20]*fUpwind_r[58]+0.1767766952966368*alphar[1]*fUpwind_r[55]+0.1767766952966368*alphar[2]*fUpwind_r[54]+0.1767766952966368*alphar[3]*fUpwind_r[53]+0.1767766952966368*alphar[6]*fUpwind_r[51]+0.1767766952966368*alphar[7]*fUpwind_r[50]+0.1767766952966368*alphar[8]*fUpwind_r[49]+0.1767766952966368*alphar[16]*fUpwind_r[48]+0.1581138830084189*alphar[4]*fUpwind_r[26]+0.1581138830084189*alphar[9]*fUpwind_r[19]+0.1581138830084189*alphar[10]*fUpwind_r[18]+0.1581138830084189*fUpwind_r[11]*alphar[17]; 
  Ghat_r[60] = 0.1767766952966368*alphar[3]*fUpwind_r[63]+0.1767766952966368*alphar[7]*fUpwind_r[62]+0.1767766952966368*alphar[8]*fUpwind_r[61]+0.1767766952966368*alphar[0]*fUpwind_r[60]+0.1767766952966368*alphar[16]*fUpwind_r[58]+0.1767766952966368*alphar[1]*fUpwind_r[57]+0.1767766952966368*alphar[2]*fUpwind_r[56]+0.1767766952966368*alphar[5]*fUpwind_r[53]+0.1767766952966368*alphar[6]*fUpwind_r[52]+0.1767766952966368*alphar[12]*fUpwind_r[50]+0.1767766952966368*alphar[13]*fUpwind_r[49]+0.1767766952966368*alphar[20]*fUpwind_r[48]+0.1581138830084189*alphar[4]*fUpwind_r[28]+0.1581138830084189*alphar[9]*fUpwind_r[24]+0.1581138830084189*alphar[10]*fUpwind_r[23]+0.1581138830084189*fUpwind_r[15]*alphar[17]; 
  Ghat_r[61] = 0.1767766952966368*alphar[2]*fUpwind_r[63]+0.1767766952966368*alphar[6]*fUpwind_r[62]+0.1767766952966368*alphar[0]*fUpwind_r[61]+0.1767766952966368*alphar[8]*fUpwind_r[60]+0.1767766952966368*alphar[13]*fUpwind_r[59]+0.1767766952966368*alphar[1]*fUpwind_r[58]+0.1767766952966368*alphar[16]*fUpwind_r[57]+0.1767766952966368*alphar[3]*fUpwind_r[56]+0.1767766952966368*alphar[20]*fUpwind_r[55]+0.1767766952966368*alphar[5]*fUpwind_r[54]+0.1767766952966368*alphar[7]*fUpwind_r[52]+0.1767766952966368*alphar[12]*fUpwind_r[51]+0.1581138830084189*alphar[10]*fUpwind_r[31]+0.1581138830084189*alphar[17]*fUpwind_r[30]+0.1581138830084189*alphar[4]*fUpwind_r[29]+0.1581138830084189*alphar[9]*fUpwind_r[25]; 
  Ghat_r[62] = 0.1767766952966368*alphar[1]*fUpwind_r[63]+0.1767766952966368*alphar[0]*fUpwind_r[62]+0.1767766952966368*alphar[6]*fUpwind_r[61]+0.1767766952966368*alphar[7]*fUpwind_r[60]+0.1767766952966368*alphar[12]*fUpwind_r[59]+0.1767766952966368*alphar[2]*fUpwind_r[58]+0.1767766952966368*alphar[3]*fUpwind_r[57]+0.1767766952966368*alphar[16]*fUpwind_r[56]+0.1767766952966368*alphar[5]*fUpwind_r[55]+0.1767766952966368*alphar[20]*fUpwind_r[54]+0.1767766952966368*alphar[8]*fUpwind_r[52]+0.1767766952966368*alphar[13]*fUpwind_r[51]+0.1581138830084189*alphar[9]*fUpwind_r[31]+0.1581138830084189*alphar[4]*fUpwind_r[30]+0.1581138830084189*alphar[17]*fUpwind_r[29]+0.1581138830084189*alphar[10]*fUpwind_r[25]; 
  Ghat_r[63] = 0.1767766952966368*alphar[0]*fUpwind_r[63]+0.1767766952966368*alphar[1]*fUpwind_r[62]+0.1767766952966368*alphar[2]*fUpwind_r[61]+0.1767766952966368*alphar[3]*fUpwind_r[60]+0.1767766952966368*alphar[5]*fUpwind_r[59]+0.1767766952966368*alphar[6]*fUpwind_r[58]+0.1767766952966368*alphar[7]*fUpwind_r[57]+0.1767766952966368*alphar[8]*fUpwind_r[56]+0.1767766952966368*alphar[12]*fUpwind_r[55]+0.1767766952966368*alphar[13]*fUpwind_r[54]+0.1767766952966368*alphar[16]*fUpwind_r[52]+0.1767766952966368*alphar[20]*fUpwind_r[51]+0.1581138830084189*alphar[4]*fUpwind_r[31]+0.1581138830084189*alphar[9]*fUpwind_r[30]+0.1581138830084189*alphar[10]*fUpwind_r[29]+0.1581138830084189*alphar[17]*fUpwind_r[25]; 
  Ghat_r[64] = 0.1767766952966368*alphar[17]*fUpwind_r[76]+0.1767766952966368*alphar[16]*fUpwind_r[75]+0.1767766952966368*alphar[10]*fUpwind_r[73]+0.1767766952966368*alphar[9]*fUpwind_r[72]+0.1767766952966368*alphar[8]*fUpwind_r[71]+0.1767766952966368*alphar[7]*fUpwind_r[70]+0.1767766952966368*alphar[6]*fUpwind_r[69]+0.1767766952966368*alphar[4]*fUpwind_r[68]+0.1767766952966368*alphar[3]*fUpwind_r[67]+0.1767766952966368*alphar[2]*fUpwind_r[66]+0.1767766952966368*alphar[1]*fUpwind_r[65]+0.1767766952966368*alphar[0]*fUpwind_r[64]+0.1581138830084189*alphar[20]*fUpwind_r[20]+0.1581138830084189*alphar[13]*fUpwind_r[13]+0.1581138830084189*alphar[12]*fUpwind_r[12]+0.1581138830084189*alphar[5]*fUpwind_r[5]; 
  Ghat_r[65] = 0.1767766952966368*alphar[10]*fUpwind_r[76]+0.1767766952966368*alphar[8]*fUpwind_r[75]+0.1767766952966368*alphar[17]*fUpwind_r[73]+0.1767766952966368*alphar[4]*fUpwind_r[72]+0.1767766952966368*alphar[16]*fUpwind_r[71]+0.1767766952966368*alphar[3]*fUpwind_r[70]+0.1767766952966368*alphar[2]*fUpwind_r[69]+0.1767766952966368*alphar[9]*fUpwind_r[68]+0.1767766952966368*alphar[7]*fUpwind_r[67]+0.1767766952966368*alphar[6]*fUpwind_r[66]+0.1767766952966368*alphar[0]*fUpwind_r[65]+0.1767766952966368*alphar[1]*fUpwind_r[64]+0.1581138830084189*alphar[13]*fUpwind_r[20]+0.1581138830084189*fUpwind_r[13]*alphar[20]+0.1581138830084189*alphar[5]*fUpwind_r[12]+0.1581138830084189*fUpwind_r[5]*alphar[12]; 
  Ghat_r[66] = 0.1767766952966368*alphar[9]*fUpwind_r[76]+0.1767766952966368*alphar[7]*fUpwind_r[75]+0.1767766952966368*alphar[4]*fUpwind_r[73]+0.1767766952966368*alphar[17]*fUpwind_r[72]+0.1767766952966368*alphar[3]*fUpwind_r[71]+0.1767766952966368*alphar[16]*fUpwind_r[70]+0.1767766952966368*alphar[1]*fUpwind_r[69]+0.1767766952966368*alphar[10]*fUpwind_r[68]+0.1767766952966368*alphar[8]*fUpwind_r[67]+0.1767766952966368*alphar[0]*fUpwind_r[66]+0.1767766952966368*alphar[6]*fUpwind_r[65]+0.1767766952966368*alphar[2]*fUpwind_r[64]+0.1581138830084189*alphar[12]*fUpwind_r[20]+0.1581138830084189*fUpwind_r[12]*alphar[20]+0.1581138830084189*alphar[5]*fUpwind_r[13]+0.1581138830084189*fUpwind_r[5]*alphar[13]; 
  Ghat_r[67] = 0.1767766952966368*alphar[17]*fUpwind_r[79]+0.1767766952966368*alphar[10]*fUpwind_r[78]+0.1767766952966368*alphar[9]*fUpwind_r[77]+0.1767766952966368*alphar[6]*fUpwind_r[75]+0.1767766952966368*alphar[4]*fUpwind_r[74]+0.1767766952966368*alphar[2]*fUpwind_r[71]+0.1767766952966368*alphar[1]*fUpwind_r[70]+0.1767766952966368*alphar[16]*fUpwind_r[69]+0.1767766952966368*alphar[0]*fUpwind_r[67]+0.1767766952966368*alphar[8]*fUpwind_r[66]+0.1767766952966368*alphar[7]*fUpwind_r[65]+0.1767766952966368*alphar[3]*fUpwind_r[64]+0.1581138830084189*alphar[20]*fUpwind_r[27]+0.1581138830084189*alphar[13]*fUpwind_r[22]+0.1581138830084189*alphar[12]*fUpwind_r[21]+0.1581138830084189*alphar[5]*fUpwind_r[14]; 
  Ghat_r[68] = 0.1767766952966368*alphar[16]*fUpwind_r[79]+0.1767766952966368*alphar[8]*fUpwind_r[78]+0.1767766952966368*alphar[7]*fUpwind_r[77]+0.1767766952966368*alphar[6]*fUpwind_r[76]+0.1767766952966368*alphar[3]*fUpwind_r[74]+0.1767766952966368*alphar[2]*fUpwind_r[73]+0.1767766952966368*alphar[1]*fUpwind_r[72]+0.1767766952966368*alphar[17]*fUpwind_r[69]+0.1767766952966368*alphar[0]*fUpwind_r[68]+0.1767766952966368*alphar[10]*fUpwind_r[66]+0.1767766952966368*alphar[9]*fUpwind_r[65]+0.1767766952966368*alphar[4]*fUpwind_r[64]+0.1581138830084189*alphar[20]*fUpwind_r[28]+0.1581138830084189*alphar[13]*fUpwind_r[24]+0.1581138830084189*alphar[12]*fUpwind_r[23]+0.1581138830084189*alphar[5]*fUpwind_r[15]; 
  Ghat_r[69] = 0.1767766952966368*alphar[4]*fUpwind_r[76]+0.1767766952966368*alphar[3]*fUpwind_r[75]+0.1767766952966368*alphar[9]*fUpwind_r[73]+0.1767766952966368*alphar[10]*fUpwind_r[72]+0.1767766952966368*alphar[7]*fUpwind_r[71]+0.1767766952966368*alphar[8]*fUpwind_r[70]+0.1767766952966368*alphar[0]*fUpwind_r[69]+0.1767766952966368*alphar[17]*fUpwind_r[68]+0.1767766952966368*alphar[16]*fUpwind_r[67]+0.1767766952966368*alphar[1]*fUpwind_r[66]+0.1767766952966368*alphar[2]*fUpwind_r[65]+0.1767766952966368*alphar[6]*fUpwind_r[64]+0.1581138830084189*alphar[5]*fUpwind_r[20]+0.1581138830084189*fUpwind_r[5]*alphar[20]+0.1581138830084189*alphar[12]*fUpwind_r[13]+0.1581138830084189*fUpwind_r[12]*alphar[13]; 
  Ghat_r[70] = 0.1767766952966368*alphar[10]*fUpwind_r[79]+0.1767766952966368*alphar[17]*fUpwind_r[78]+0.1767766952966368*alphar[4]*fUpwind_r[77]+0.1767766952966368*alphar[2]*fUpwind_r[75]+0.1767766952966368*alphar[9]*fUpwind_r[74]+0.1767766952966368*alphar[6]*fUpwind_r[71]+0.1767766952966368*alphar[0]*fUpwind_r[70]+0.1767766952966368*alphar[8]*fUpwind_r[69]+0.1767766952966368*alphar[1]*fUpwind_r[67]+0.1767766952966368*alphar[16]*fUpwind_r[66]+0.1767766952966368*alphar[3]*fUpwind_r[65]+0.1767766952966368*alphar[7]*fUpwind_r[64]+0.1581138830084189*alphar[13]*fUpwind_r[27]+0.1581138830084189*alphar[20]*fUpwind_r[22]+0.1581138830084189*alphar[5]*fUpwind_r[21]+0.1581138830084189*alphar[12]*fUpwind_r[14]; 
  Ghat_r[71] = 0.1767766952966368*alphar[9]*fUpwind_r[79]+0.1767766952966368*alphar[4]*fUpwind_r[78]+0.1767766952966368*alphar[17]*fUpwind_r[77]+0.1767766952966368*alphar[1]*fUpwind_r[75]+0.1767766952966368*alphar[10]*fUpwind_r[74]+0.1767766952966368*alphar[0]*fUpwind_r[71]+0.1767766952966368*alphar[6]*fUpwind_r[70]+0.1767766952966368*alphar[7]*fUpwind_r[69]+0.1767766952966368*alphar[2]*fUpwind_r[67]+0.1767766952966368*alphar[3]*fUpwind_r[66]+0.1767766952966368*alphar[16]*fUpwind_r[65]+0.1767766952966368*alphar[8]*fUpwind_r[64]+0.1581138830084189*alphar[12]*fUpwind_r[27]+0.1581138830084189*alphar[5]*fUpwind_r[22]+0.1581138830084189*alphar[20]*fUpwind_r[21]+0.1581138830084189*alphar[13]*fUpwind_r[14]; 
  Ghat_r[72] = 0.1767766952966368*alphar[8]*fUpwind_r[79]+0.1767766952966368*alphar[16]*fUpwind_r[78]+0.1767766952966368*alphar[3]*fUpwind_r[77]+0.1767766952966368*alphar[2]*fUpwind_r[76]+0.1767766952966368*alphar[7]*fUpwind_r[74]+0.1767766952966368*alphar[6]*fUpwind_r[73]+0.1767766952966368*alphar[0]*fUpwind_r[72]+0.1767766952966368*alphar[10]*fUpwind_r[69]+0.1767766952966368*alphar[1]*fUpwind_r[68]+0.1767766952966368*alphar[17]*fUpwind_r[66]+0.1767766952966368*alphar[4]*fUpwind_r[65]+0.1767766952966368*alphar[9]*fUpwind_r[64]+0.1581138830084189*alphar[13]*fUpwind_r[28]+0.1581138830084189*alphar[20]*fUpwind_r[24]+0.1581138830084189*alphar[5]*fUpwind_r[23]+0.1581138830084189*alphar[12]*fUpwind_r[15]; 
  Ghat_r[73] = 0.1767766952966368*alphar[7]*fUpwind_r[79]+0.1767766952966368*alphar[3]*fUpwind_r[78]+0.1767766952966368*alphar[16]*fUpwind_r[77]+0.1767766952966368*alphar[1]*fUpwind_r[76]+0.1767766952966368*alphar[8]*fUpwind_r[74]+0.1767766952966368*alphar[0]*fUpwind_r[73]+0.1767766952966368*alphar[6]*fUpwind_r[72]+0.1767766952966368*alphar[9]*fUpwind_r[69]+0.1767766952966368*alphar[2]*fUpwind_r[68]+0.1767766952966368*alphar[4]*fUpwind_r[66]+0.1767766952966368*alphar[17]*fUpwind_r[65]+0.1767766952966368*alphar[10]*fUpwind_r[64]+0.1581138830084189*alphar[12]*fUpwind_r[28]+0.1581138830084189*alphar[5]*fUpwind_r[24]+0.1581138830084189*alphar[20]*fUpwind_r[23]+0.1581138830084189*alphar[13]*fUpwind_r[15]; 
  Ghat_r[74] = 0.1767766952966368*alphar[6]*fUpwind_r[79]+0.1767766952966368*alphar[2]*fUpwind_r[78]+0.1767766952966368*alphar[1]*fUpwind_r[77]+0.1767766952966368*alphar[16]*fUpwind_r[76]+0.1767766952966368*alphar[17]*fUpwind_r[75]+0.1767766952966368*alphar[0]*fUpwind_r[74]+0.1767766952966368*alphar[8]*fUpwind_r[73]+0.1767766952966368*alphar[7]*fUpwind_r[72]+0.1767766952966368*alphar[10]*fUpwind_r[71]+0.1767766952966368*alphar[9]*fUpwind_r[70]+0.1767766952966368*alphar[3]*fUpwind_r[68]+0.1767766952966368*alphar[4]*fUpwind_r[67]+0.1581138830084189*alphar[20]*fUpwind_r[31]+0.1581138830084189*alphar[13]*fUpwind_r[30]+0.1581138830084189*alphar[12]*fUpwind_r[29]+0.1581138830084189*alphar[5]*fUpwind_r[25]; 
  Ghat_r[75] = 0.1767766952966368*alphar[4]*fUpwind_r[79]+0.1767766952966368*alphar[9]*fUpwind_r[78]+0.1767766952966368*alphar[10]*fUpwind_r[77]+0.1767766952966368*alphar[0]*fUpwind_r[75]+0.1767766952966368*alphar[17]*fUpwind_r[74]+0.1767766952966368*alphar[1]*fUpwind_r[71]+0.1767766952966368*alphar[2]*fUpwind_r[70]+0.1767766952966368*alphar[3]*fUpwind_r[69]+0.1767766952966368*alphar[6]*fUpwind_r[67]+0.1767766952966368*alphar[7]*fUpwind_r[66]+0.1767766952966368*alphar[8]*fUpwind_r[65]+0.1767766952966368*alphar[16]*fUpwind_r[64]+0.1581138830084189*alphar[5]*fUpwind_r[27]+0.1581138830084189*alphar[12]*fUpwind_r[22]+0.1581138830084189*alphar[13]*fUpwind_r[21]+0.1581138830084189*fUpwind_r[14]*alphar[20]; 
  Ghat_r[76] = 0.1767766952966368*alphar[3]*fUpwind_r[79]+0.1767766952966368*alphar[7]*fUpwind_r[78]+0.1767766952966368*alphar[8]*fUpwind_r[77]+0.1767766952966368*alphar[0]*fUpwind_r[76]+0.1767766952966368*alphar[16]*fUpwind_r[74]+0.1767766952966368*alphar[1]*fUpwind_r[73]+0.1767766952966368*alphar[2]*fUpwind_r[72]+0.1767766952966368*alphar[4]*fUpwind_r[69]+0.1767766952966368*alphar[6]*fUpwind_r[68]+0.1767766952966368*alphar[9]*fUpwind_r[66]+0.1767766952966368*alphar[10]*fUpwind_r[65]+0.1767766952966368*alphar[17]*fUpwind_r[64]+0.1581138830084189*alphar[5]*fUpwind_r[28]+0.1581138830084189*alphar[12]*fUpwind_r[24]+0.1581138830084189*alphar[13]*fUpwind_r[23]+0.1581138830084189*fUpwind_r[15]*alphar[20]; 
  Ghat_r[77] = 0.1767766952966368*alphar[2]*fUpwind_r[79]+0.1767766952966368*alphar[6]*fUpwind_r[78]+0.1767766952966368*alphar[0]*fUpwind_r[77]+0.1767766952966368*alphar[8]*fUpwind_r[76]+0.1767766952966368*alphar[10]*fUpwind_r[75]+0.1767766952966368*alphar[1]*fUpwind_r[74]+0.1767766952966368*alphar[16]*fUpwind_r[73]+0.1767766952966368*alphar[3]*fUpwind_r[72]+0.1767766952966368*alphar[17]*fUpwind_r[71]+0.1767766952966368*alphar[4]*fUpwind_r[70]+0.1767766952966368*alphar[7]*fUpwind_r[68]+0.1767766952966368*alphar[9]*fUpwind_r[67]+0.1581138830084189*alphar[13]*fUpwind_r[31]+0.1581138830084189*alphar[20]*fUpwind_r[30]+0.1581138830084189*alphar[5]*fUpwind_r[29]+0.1581138830084189*alphar[12]*fUpwind_r[25]; 
  Ghat_r[78] = 0.1767766952966368*alphar[1]*fUpwind_r[79]+0.1767766952966368*alphar[0]*fUpwind_r[78]+0.1767766952966368*alphar[6]*fUpwind_r[77]+0.1767766952966368*alphar[7]*fUpwind_r[76]+0.1767766952966368*alphar[9]*fUpwind_r[75]+0.1767766952966368*alphar[2]*fUpwind_r[74]+0.1767766952966368*alphar[3]*fUpwind_r[73]+0.1767766952966368*alphar[16]*fUpwind_r[72]+0.1767766952966368*alphar[4]*fUpwind_r[71]+0.1767766952966368*alphar[17]*fUpwind_r[70]+0.1767766952966368*alphar[8]*fUpwind_r[68]+0.1767766952966368*alphar[10]*fUpwind_r[67]+0.1581138830084189*alphar[12]*fUpwind_r[31]+0.1581138830084189*alphar[5]*fUpwind_r[30]+0.1581138830084189*alphar[20]*fUpwind_r[29]+0.1581138830084189*alphar[13]*fUpwind_r[25]; 
  Ghat_r[79] = 0.1767766952966368*alphar[0]*fUpwind_r[79]+0.1767766952966368*alphar[1]*fUpwind_r[78]+0.1767766952966368*alphar[2]*fUpwind_r[77]+0.1767766952966368*alphar[3]*fUpwind_r[76]+0.1767766952966368*alphar[4]*fUpwind_r[75]+0.1767766952966368*alphar[6]*fUpwind_r[74]+0.1767766952966368*alphar[7]*fUpwind_r[73]+0.1767766952966368*alphar[8]*fUpwind_r[72]+0.1767766952966368*alphar[9]*fUpwind_r[71]+0.1767766952966368*alphar[10]*fUpwind_r[70]+0.1767766952966368*alphar[16]*fUpwind_r[68]+0.1767766952966368*alphar[17]*fUpwind_r[67]+0.1581138830084189*alphar[5]*fUpwind_r[31]+0.1581138830084189*alphar[12]*fUpwind_r[30]+0.1581138830084189*alphar[13]*fUpwind_r[29]+0.1581138830084189*alphar[20]*fUpwind_r[25]; 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dx10; 
  out[1] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dx10; 
  out[2] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dx10; 
  out[3] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dx10; 
  out[4] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dx10; 
  out[5] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*dx10; 
  out[6] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*dx10; 
  out[7] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dx10; 
  out[8] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dx10; 
  out[9] += (0.7071067811865475*Ghat_l[6]-0.7071067811865475*Ghat_r[6])*dx10; 
  out[10] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dx10; 
  out[11] += (0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*dx10; 
  out[12] += (0.7071067811865475*Ghat_l[8]-0.7071067811865475*Ghat_r[8])*dx10; 
  out[13] += -1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dx10; 
  out[14] += (0.7071067811865475*Ghat_l[9]-0.7071067811865475*Ghat_r[9])*dx10; 
  out[15] += (0.7071067811865475*Ghat_l[10]-0.7071067811865475*Ghat_r[10])*dx10; 
  out[16] += (0.7071067811865475*Ghat_l[11]-0.7071067811865475*Ghat_r[11])*dx10; 
  out[17] += -1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dx10; 
  out[18] += (0.7071067811865475*Ghat_l[12]-0.7071067811865475*Ghat_r[12])*dx10; 
  out[19] += (0.7071067811865475*Ghat_l[13]-0.7071067811865475*Ghat_r[13])*dx10; 
  out[20] += (0.7071067811865475*Ghat_l[14]-0.7071067811865475*Ghat_r[14])*dx10; 
  out[21] += (0.7071067811865475*Ghat_l[15]-0.7071067811865475*Ghat_r[15])*dx10; 
  out[22] += -1.224744871391589*(Ghat_r[6]+Ghat_l[6])*dx10; 
  out[23] += -1.224744871391589*(Ghat_r[7]+Ghat_l[7])*dx10; 
  out[24] += -1.224744871391589*(Ghat_r[8]+Ghat_l[8])*dx10; 
  out[25] += (0.7071067811865475*Ghat_l[16]-0.7071067811865475*Ghat_r[16])*dx10; 
  out[26] += -1.224744871391589*(Ghat_r[9]+Ghat_l[9])*dx10; 
  out[27] += -1.224744871391589*(Ghat_r[10]+Ghat_l[10])*dx10; 
  out[28] += (0.7071067811865475*Ghat_l[17]-0.7071067811865475*Ghat_r[17])*dx10; 
  out[29] += -1.224744871391589*(Ghat_r[11]+Ghat_l[11])*dx10; 
  out[30] += (0.7071067811865475*Ghat_l[18]-0.7071067811865475*Ghat_r[18])*dx10; 
  out[31] += (0.7071067811865475*Ghat_l[19]-0.7071067811865475*Ghat_r[19])*dx10; 
  out[32] += -1.224744871391589*(Ghat_r[12]+Ghat_l[12])*dx10; 
  out[33] += -1.224744871391589*(Ghat_r[13]+Ghat_l[13])*dx10; 
  out[34] += (0.7071067811865475*Ghat_l[20]-0.7071067811865475*Ghat_r[20])*dx10; 
  out[35] += -1.224744871391589*(Ghat_r[14]+Ghat_l[14])*dx10; 
  out[36] += (0.7071067811865475*Ghat_l[21]-0.7071067811865475*Ghat_r[21])*dx10; 
  out[37] += (0.7071067811865475*Ghat_l[22]-0.7071067811865475*Ghat_r[22])*dx10; 
  out[38] += -1.224744871391589*(Ghat_r[15]+Ghat_l[15])*dx10; 
  out[39] += (0.7071067811865475*Ghat_l[23]-0.7071067811865475*Ghat_r[23])*dx10; 
  out[40] += (0.7071067811865475*Ghat_l[24]-0.7071067811865475*Ghat_r[24])*dx10; 
  out[41] += (0.7071067811865475*Ghat_l[25]-0.7071067811865475*Ghat_r[25])*dx10; 
  out[42] += -1.224744871391589*(Ghat_r[16]+Ghat_l[16])*dx10; 
  out[43] += -1.224744871391589*(Ghat_r[17]+Ghat_l[17])*dx10; 
  out[44] += -1.224744871391589*(Ghat_r[18]+Ghat_l[18])*dx10; 
  out[45] += -1.224744871391589*(Ghat_r[19]+Ghat_l[19])*dx10; 
  out[46] += (0.7071067811865475*Ghat_l[26]-0.7071067811865475*Ghat_r[26])*dx10; 
  out[47] += -1.224744871391589*(Ghat_r[20]+Ghat_l[20])*dx10; 
  out[48] += -1.224744871391589*(Ghat_r[21]+Ghat_l[21])*dx10; 
  out[49] += -1.224744871391589*(Ghat_r[22]+Ghat_l[22])*dx10; 
  out[50] += (0.7071067811865475*Ghat_l[27]-0.7071067811865475*Ghat_r[27])*dx10; 
  out[51] += -1.224744871391589*(Ghat_r[23]+Ghat_l[23])*dx10; 
  out[52] += -1.224744871391589*(Ghat_r[24]+Ghat_l[24])*dx10; 
  out[53] += (0.7071067811865475*Ghat_l[28]-0.7071067811865475*Ghat_r[28])*dx10; 
  out[54] += -1.224744871391589*(Ghat_r[25]+Ghat_l[25])*dx10; 
  out[55] += (0.7071067811865475*Ghat_l[29]-0.7071067811865475*Ghat_r[29])*dx10; 
  out[56] += (0.7071067811865475*Ghat_l[30]-0.7071067811865475*Ghat_r[30])*dx10; 
  out[57] += -1.224744871391589*(Ghat_r[26]+Ghat_l[26])*dx10; 
  out[58] += -1.224744871391589*(Ghat_r[27]+Ghat_l[27])*dx10; 
  out[59] += -1.224744871391589*(Ghat_r[28]+Ghat_l[28])*dx10; 
  out[60] += -1.224744871391589*(Ghat_r[29]+Ghat_l[29])*dx10; 
  out[61] += -1.224744871391589*(Ghat_r[30]+Ghat_l[30])*dx10; 
  out[62] += (0.7071067811865475*Ghat_l[31]-0.7071067811865475*Ghat_r[31])*dx10; 
  out[63] += -1.224744871391589*(Ghat_r[31]+Ghat_l[31])*dx10; 
  out[64] += (0.7071067811865475*Ghat_l[32]-0.7071067811865475*Ghat_r[32])*dx10; 
  out[65] += -1.224744871391589*(Ghat_r[32]+Ghat_l[32])*dx10; 
  out[66] += (0.7071067811865475*Ghat_l[33]-0.7071067811865475*Ghat_r[33])*dx10; 
  out[67] += (0.7071067811865475*Ghat_l[34]-0.7071067811865475*Ghat_r[34])*dx10; 
  out[68] += (0.7071067811865475*Ghat_l[35]-0.7071067811865475*Ghat_r[35])*dx10; 
  out[69] += (0.7071067811865475*Ghat_l[36]-0.7071067811865475*Ghat_r[36])*dx10; 
  out[70] += -1.224744871391589*(Ghat_r[33]+Ghat_l[33])*dx10; 
  out[71] += -1.224744871391589*(Ghat_r[34]+Ghat_l[34])*dx10; 
  out[72] += (0.7071067811865475*Ghat_l[37]-0.7071067811865475*Ghat_r[37])*dx10; 
  out[73] += -1.224744871391589*(Ghat_r[35]+Ghat_l[35])*dx10; 
  out[74] += (0.7071067811865475*Ghat_l[38]-0.7071067811865475*Ghat_r[38])*dx10; 
  out[75] += (0.7071067811865475*Ghat_l[39]-0.7071067811865475*Ghat_r[39])*dx10; 
  out[76] += -1.224744871391589*(Ghat_r[36]+Ghat_l[36])*dx10; 
  out[77] += (0.7071067811865475*Ghat_l[40]-0.7071067811865475*Ghat_r[40])*dx10; 
  out[78] += (0.7071067811865475*Ghat_l[41]-0.7071067811865475*Ghat_r[41])*dx10; 
  out[79] += (0.7071067811865475*Ghat_l[42]-0.7071067811865475*Ghat_r[42])*dx10; 
  out[80] += -1.224744871391589*(Ghat_r[37]+Ghat_l[37])*dx10; 
  out[81] += -1.224744871391589*(Ghat_r[38]+Ghat_l[38])*dx10; 
  out[82] += -1.224744871391589*(Ghat_r[39]+Ghat_l[39])*dx10; 
  out[83] += (0.7071067811865475*Ghat_l[43]-0.7071067811865475*Ghat_r[43])*dx10; 
  out[84] += -1.224744871391589*(Ghat_r[40]+Ghat_l[40])*dx10; 
  out[85] += -1.224744871391589*(Ghat_r[41]+Ghat_l[41])*dx10; 
  out[86] += (0.7071067811865475*Ghat_l[44]-0.7071067811865475*Ghat_r[44])*dx10; 
  out[87] += -1.224744871391589*(Ghat_r[42]+Ghat_l[42])*dx10; 
  out[88] += (0.7071067811865475*Ghat_l[45]-0.7071067811865475*Ghat_r[45])*dx10; 
  out[89] += (0.7071067811865475*Ghat_l[46]-0.7071067811865475*Ghat_r[46])*dx10; 
  out[90] += -1.224744871391589*(Ghat_r[43]+Ghat_l[43])*dx10; 
  out[91] += -1.224744871391589*(Ghat_r[44]+Ghat_l[44])*dx10; 
  out[92] += -1.224744871391589*(Ghat_r[45]+Ghat_l[45])*dx10; 
  out[93] += -1.224744871391589*(Ghat_r[46]+Ghat_l[46])*dx10; 
  out[94] += (0.7071067811865475*Ghat_l[47]-0.7071067811865475*Ghat_r[47])*dx10; 
  out[95] += -1.224744871391589*(Ghat_r[47]+Ghat_l[47])*dx10; 
  out[96] += (0.7071067811865475*Ghat_l[48]-0.7071067811865475*Ghat_r[48])*dx10; 
  out[97] += -1.224744871391589*(Ghat_r[48]+Ghat_l[48])*dx10; 
  out[98] += (0.7071067811865475*Ghat_l[49]-0.7071067811865475*Ghat_r[49])*dx10; 
  out[99] += (0.7071067811865475*Ghat_l[50]-0.7071067811865475*Ghat_r[50])*dx10; 
  out[100] += (0.7071067811865475*Ghat_l[51]-0.7071067811865475*Ghat_r[51])*dx10; 
  out[101] += (0.7071067811865475*Ghat_l[52]-0.7071067811865475*Ghat_r[52])*dx10; 
  out[102] += -1.224744871391589*(Ghat_r[49]+Ghat_l[49])*dx10; 
  out[103] += -1.224744871391589*(Ghat_r[50]+Ghat_l[50])*dx10; 
  out[104] += (0.7071067811865475*Ghat_l[53]-0.7071067811865475*Ghat_r[53])*dx10; 
  out[105] += -1.224744871391589*(Ghat_r[51]+Ghat_l[51])*dx10; 
  out[106] += (0.7071067811865475*Ghat_l[54]-0.7071067811865475*Ghat_r[54])*dx10; 
  out[107] += (0.7071067811865475*Ghat_l[55]-0.7071067811865475*Ghat_r[55])*dx10; 
  out[108] += -1.224744871391589*(Ghat_r[52]+Ghat_l[52])*dx10; 
  out[109] += (0.7071067811865475*Ghat_l[56]-0.7071067811865475*Ghat_r[56])*dx10; 
  out[110] += (0.7071067811865475*Ghat_l[57]-0.7071067811865475*Ghat_r[57])*dx10; 
  out[111] += (0.7071067811865475*Ghat_l[58]-0.7071067811865475*Ghat_r[58])*dx10; 
  out[112] += -1.224744871391589*(Ghat_r[53]+Ghat_l[53])*dx10; 
  out[113] += -1.224744871391589*(Ghat_r[54]+Ghat_l[54])*dx10; 
  out[114] += -1.224744871391589*(Ghat_r[55]+Ghat_l[55])*dx10; 
  out[115] += (0.7071067811865475*Ghat_l[59]-0.7071067811865475*Ghat_r[59])*dx10; 
  out[116] += -1.224744871391589*(Ghat_r[56]+Ghat_l[56])*dx10; 
  out[117] += -1.224744871391589*(Ghat_r[57]+Ghat_l[57])*dx10; 
  out[118] += (0.7071067811865475*Ghat_l[60]-0.7071067811865475*Ghat_r[60])*dx10; 
  out[119] += -1.224744871391589*(Ghat_r[58]+Ghat_l[58])*dx10; 
  out[120] += (0.7071067811865475*Ghat_l[61]-0.7071067811865475*Ghat_r[61])*dx10; 
  out[121] += (0.7071067811865475*Ghat_l[62]-0.7071067811865475*Ghat_r[62])*dx10; 
  out[122] += -1.224744871391589*(Ghat_r[59]+Ghat_l[59])*dx10; 
  out[123] += -1.224744871391589*(Ghat_r[60]+Ghat_l[60])*dx10; 
  out[124] += -1.224744871391589*(Ghat_r[61]+Ghat_l[61])*dx10; 
  out[125] += -1.224744871391589*(Ghat_r[62]+Ghat_l[62])*dx10; 
  out[126] += (0.7071067811865475*Ghat_l[63]-0.7071067811865475*Ghat_r[63])*dx10; 
  out[127] += -1.224744871391589*(Ghat_r[63]+Ghat_l[63])*dx10; 
  out[128] += (0.7071067811865475*Ghat_l[64]-0.7071067811865475*Ghat_r[64])*dx10; 
  out[129] += -1.224744871391589*(Ghat_r[64]+Ghat_l[64])*dx10; 
  out[130] += (0.7071067811865475*Ghat_l[65]-0.7071067811865475*Ghat_r[65])*dx10; 
  out[131] += (0.7071067811865475*Ghat_l[66]-0.7071067811865475*Ghat_r[66])*dx10; 
  out[132] += (0.7071067811865475*Ghat_l[67]-0.7071067811865475*Ghat_r[67])*dx10; 
  out[133] += (0.7071067811865475*Ghat_l[68]-0.7071067811865475*Ghat_r[68])*dx10; 
  out[134] += -1.224744871391589*(Ghat_r[65]+Ghat_l[65])*dx10; 
  out[135] += -1.224744871391589*(Ghat_r[66]+Ghat_l[66])*dx10; 
  out[136] += (0.7071067811865475*Ghat_l[69]-0.7071067811865475*Ghat_r[69])*dx10; 
  out[137] += -1.224744871391589*(Ghat_r[67]+Ghat_l[67])*dx10; 
  out[138] += (0.7071067811865475*Ghat_l[70]-0.7071067811865475*Ghat_r[70])*dx10; 
  out[139] += (0.7071067811865475*Ghat_l[71]-0.7071067811865475*Ghat_r[71])*dx10; 
  out[140] += -1.224744871391589*(Ghat_r[68]+Ghat_l[68])*dx10; 
  out[141] += (0.7071067811865475*Ghat_l[72]-0.7071067811865475*Ghat_r[72])*dx10; 
  out[142] += (0.7071067811865475*Ghat_l[73]-0.7071067811865475*Ghat_r[73])*dx10; 
  out[143] += (0.7071067811865475*Ghat_l[74]-0.7071067811865475*Ghat_r[74])*dx10; 
  out[144] += -1.224744871391589*(Ghat_r[69]+Ghat_l[69])*dx10; 
  out[145] += -1.224744871391589*(Ghat_r[70]+Ghat_l[70])*dx10; 
  out[146] += -1.224744871391589*(Ghat_r[71]+Ghat_l[71])*dx10; 
  out[147] += (0.7071067811865475*Ghat_l[75]-0.7071067811865475*Ghat_r[75])*dx10; 
  out[148] += -1.224744871391589*(Ghat_r[72]+Ghat_l[72])*dx10; 
  out[149] += -1.224744871391589*(Ghat_r[73]+Ghat_l[73])*dx10; 
  out[150] += (0.7071067811865475*Ghat_l[76]-0.7071067811865475*Ghat_r[76])*dx10; 
  out[151] += -1.224744871391589*(Ghat_r[74]+Ghat_l[74])*dx10; 
  out[152] += (0.7071067811865475*Ghat_l[77]-0.7071067811865475*Ghat_r[77])*dx10; 
  out[153] += (0.7071067811865475*Ghat_l[78]-0.7071067811865475*Ghat_r[78])*dx10; 
  out[154] += -1.224744871391589*(Ghat_r[75]+Ghat_l[75])*dx10; 
  out[155] += -1.224744871391589*(Ghat_r[76]+Ghat_l[76])*dx10; 
  out[156] += -1.224744871391589*(Ghat_r[77]+Ghat_l[77])*dx10; 
  out[157] += -1.224744871391589*(Ghat_r[78]+Ghat_l[78])*dx10; 
  out[158] += (0.7071067811865475*Ghat_l[79]-0.7071067811865475*Ghat_r[79])*dx10; 
  out[159] += -1.224744871391589*(Ghat_r[79]+Ghat_l[79])*dx10; 

} 
