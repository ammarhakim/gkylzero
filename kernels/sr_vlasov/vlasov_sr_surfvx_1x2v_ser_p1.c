#include <gkyl_vlasov_sr_kernels.h> 
#include <gkyl_basis_ser_3x_p1_surfx2_eval_quad.h> 
#include <gkyl_basis_ser_3x_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH void vlasov_sr_surfvx_1x2v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // p_over_gamma:      p/gamma (velocity).
  // qmem:      q/m*EM fields.
  // fl/fc/fr:  Input Distribution function in left/center/right cells 
  // out:       Output distribution function in center cell 
  const double dv10 = 2/dxv[1]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv2 = dxv[2], wv2 = w[2]; 
  const double *E0 = &qmem[0]; 
  const double *p0_over_gamma = &p_over_gamma[0]; 
  const double *p1_over_gamma = &p_over_gamma[4]; 
  const double *B2 = &qmem[10]; 

  double alpha_l[4] = {0.0}; 
  double alpha_r[4] = {0.0}; 

  alpha_l[0] = (-1.224744871391589*B2[0]*p1_over_gamma[1])+0.7071067811865475*B2[0]*p1_over_gamma[0]+1.414213562373095*E0[0]; 
  alpha_l[1] = (-1.224744871391589*B2[1]*p1_over_gamma[1])+1.414213562373095*E0[1]+0.7071067811865475*p1_over_gamma[0]*B2[1]; 
  alpha_l[2] = 0.7071067811865475*B2[0]*p1_over_gamma[2]-1.224744871391589*B2[0]*p1_over_gamma[3]; 
  alpha_l[3] = 0.7071067811865475*B2[1]*p1_over_gamma[2]-1.224744871391589*B2[1]*p1_over_gamma[3]; 

  alpha_r[0] = 1.224744871391589*B2[0]*p1_over_gamma[1]+0.7071067811865475*B2[0]*p1_over_gamma[0]+1.414213562373095*E0[0]; 
  alpha_r[1] = 1.224744871391589*B2[1]*p1_over_gamma[1]+1.414213562373095*E0[1]+0.7071067811865475*p1_over_gamma[0]*B2[1]; 
  alpha_r[2] = 1.224744871391589*B2[0]*p1_over_gamma[3]+0.7071067811865475*B2[0]*p1_over_gamma[2]; 
  alpha_r[3] = 1.224744871391589*B2[1]*p1_over_gamma[3]+0.7071067811865475*B2[1]*p1_over_gamma[2]; 

  double fUpwindQuad_l[4] = {0.0};
  double fUpwindQuad_r[4] = {0.0};
  double fUpwind_l[4] = {0.0};;
  double fUpwind_r[4] = {0.0};
  double Ghat_l[4] = {0.0}; 
  double Ghat_r[4] = {0.0}; 

  if (alpha_l[3]-alpha_l[2]-alpha_l[1]+alpha_l[0] > 0) { 
    fUpwindQuad_l[0] = ser_3x_p1_surfx2_eval_quad_node_0_r(fl); 
  } else { 
    fUpwindQuad_l[0] = ser_3x_p1_surfx2_eval_quad_node_0_l(fc); 
  } 
  if (alpha_r[3]-alpha_r[2]-alpha_r[1]+alpha_r[0] > 0) { 
    fUpwindQuad_r[0] = ser_3x_p1_surfx2_eval_quad_node_0_r(fc); 
  } else { 
    fUpwindQuad_r[0] = ser_3x_p1_surfx2_eval_quad_node_0_l(fr); 
  } 
  if ((-alpha_l[3])+alpha_l[2]-alpha_l[1]+alpha_l[0] > 0) { 
    fUpwindQuad_l[1] = ser_3x_p1_surfx2_eval_quad_node_1_r(fl); 
  } else { 
    fUpwindQuad_l[1] = ser_3x_p1_surfx2_eval_quad_node_1_l(fc); 
  } 
  if ((-alpha_r[3])+alpha_r[2]-alpha_r[1]+alpha_r[0] > 0) { 
    fUpwindQuad_r[1] = ser_3x_p1_surfx2_eval_quad_node_1_r(fc); 
  } else { 
    fUpwindQuad_r[1] = ser_3x_p1_surfx2_eval_quad_node_1_l(fr); 
  } 
  if ((-alpha_l[3])-alpha_l[2]+alpha_l[1]+alpha_l[0] > 0) { 
    fUpwindQuad_l[2] = ser_3x_p1_surfx2_eval_quad_node_2_r(fl); 
  } else { 
    fUpwindQuad_l[2] = ser_3x_p1_surfx2_eval_quad_node_2_l(fc); 
  } 
  if ((-alpha_r[3])-alpha_r[2]+alpha_r[1]+alpha_r[0] > 0) { 
    fUpwindQuad_r[2] = ser_3x_p1_surfx2_eval_quad_node_2_r(fc); 
  } else { 
    fUpwindQuad_r[2] = ser_3x_p1_surfx2_eval_quad_node_2_l(fr); 
  } 
  if (alpha_l[3]+alpha_l[2]+alpha_l[1]+alpha_l[0] > 0) { 
    fUpwindQuad_l[3] = ser_3x_p1_surfx2_eval_quad_node_3_r(fl); 
  } else { 
    fUpwindQuad_l[3] = ser_3x_p1_surfx2_eval_quad_node_3_l(fc); 
  } 
  if (alpha_r[3]+alpha_r[2]+alpha_r[1]+alpha_r[0] > 0) { 
    fUpwindQuad_r[3] = ser_3x_p1_surfx2_eval_quad_node_3_r(fc); 
  } else { 
    fUpwindQuad_r[3] = ser_3x_p1_surfx2_eval_quad_node_3_l(fr); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_3x_p1_upwind_quad_to_modal(fUpwindQuad_l, fUpwind_l); 
  ser_3x_p1_upwind_quad_to_modal(fUpwindQuad_r, fUpwind_r); 

  Ghat_l[0] = 0.5*(alpha_l[3]*fUpwind_l[3]+alpha_l[2]*fUpwind_l[2]+alpha_l[1]*fUpwind_l[1]+alpha_l[0]*fUpwind_l[0]); 
  Ghat_l[1] = 0.5*(alpha_l[2]*fUpwind_l[3]+fUpwind_l[2]*alpha_l[3]+alpha_l[0]*fUpwind_l[1]+fUpwind_l[0]*alpha_l[1]); 
  Ghat_l[2] = 0.5*(alpha_l[1]*fUpwind_l[3]+fUpwind_l[1]*alpha_l[3]+alpha_l[0]*fUpwind_l[2]+fUpwind_l[0]*alpha_l[2]); 
  Ghat_l[3] = 0.5*(alpha_l[0]*fUpwind_l[3]+fUpwind_l[0]*alpha_l[3]+alpha_l[1]*fUpwind_l[2]+fUpwind_l[1]*alpha_l[2]); 

  Ghat_r[0] = 0.5*(alpha_r[3]*fUpwind_r[3]+alpha_r[2]*fUpwind_r[2]+alpha_r[1]*fUpwind_r[1]+alpha_r[0]*fUpwind_r[0]); 
  Ghat_r[1] = 0.5*(alpha_r[2]*fUpwind_r[3]+fUpwind_r[2]*alpha_r[3]+alpha_r[0]*fUpwind_r[1]+fUpwind_r[0]*alpha_r[1]); 
  Ghat_r[2] = 0.5*(alpha_r[1]*fUpwind_r[3]+fUpwind_r[1]*alpha_r[3]+alpha_r[0]*fUpwind_r[2]+fUpwind_r[0]*alpha_r[2]); 
  Ghat_r[3] = 0.5*(alpha_r[0]*fUpwind_r[3]+fUpwind_r[0]*alpha_r[3]+alpha_r[1]*fUpwind_r[2]+fUpwind_r[1]*alpha_r[2]); 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv10; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv10; 
  out[2] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv10; 
  out[3] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dv10; 
  out[4] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv10; 
  out[5] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dv10; 
  out[6] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dv10; 
  out[7] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dv10; 

} 
