#include <gkyl_vlasov_sr_kernels.h> 
#include <gkyl_basis_tensor_3x_p2_surfx3_eval_quad.h> 
#include <gkyl_basis_tensor_3x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH double vlasov_sr_vmap_surfvy_1x2v_tensor_p2(const double *w, const double *dxv, const double *jacob_vel_inv, const double *gamma, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:             Cell-center coordinates.
  // dxv[NDIM]:           Cell spacing.
  // jacob_vel_inv[VDIM]: Inverse velocity space Jacobian in each direction.
  // gamma:               Particle Lorentz boost factor sqrt(1 + p^2).
  // qmem:                q/m*EM fields.
  // fl/fc/fr:            Input Distribution function in left/center/right cells 
  // out:                 Output distribution function in center cell 
  const double dv10 = 2.0/dxv[1]; 
  const double dv11 = 2.0/dxv[2]; 
  const double *E1 = &qmem[3]; 
  const double *jacob_vel_inv0 = &jacob_vel_inv[0]; 
  const double *jacob_vel_inv1 = &jacob_vel_inv[3]; 
  double p0_over_gamma_l[3] = {0.0}; 
  double p0_over_gamma_r[3] = {0.0}; 
  p0_over_gamma_l[0] = 4.330127018922193*jacob_vel_inv0[1]*gamma[8]*dv10+1.936491673103709*jacob_vel_inv0[0]*gamma[7]*dv10-3.354101966249684*jacob_vel_inv0[1]*gamma[6]*dv10+1.936491673103709*jacob_vel_inv0[1]*gamma[4]*dv10-1.5*jacob_vel_inv0[0]*gamma[3]*dv10+0.8660254037844386*jacob_vel_inv0[0]*gamma[1]*dv10; 
  p0_over_gamma_l[1] = 3.872983346207417*jacob_vel_inv0[2]*gamma[8]*dv10+4.330127018922193*jacob_vel_inv0[0]*gamma[8]*dv10+1.936491673103709*jacob_vel_inv0[1]*gamma[7]*dv10-3.0*jacob_vel_inv0[2]*gamma[6]*dv10-3.354101966249684*jacob_vel_inv0[0]*gamma[6]*dv10+1.732050807568877*jacob_vel_inv0[2]*gamma[4]*dv10+1.936491673103709*jacob_vel_inv0[0]*gamma[4]*dv10-1.5*jacob_vel_inv0[1]*gamma[3]*dv10+0.8660254037844386*jacob_vel_inv0[1]*gamma[1]*dv10; 
  p0_over_gamma_l[2] = 3.872983346207417*jacob_vel_inv0[1]*gamma[8]*dv10+1.936491673103709*jacob_vel_inv0[2]*gamma[7]*dv10-3.0*jacob_vel_inv0[1]*gamma[6]*dv10+1.732050807568877*jacob_vel_inv0[1]*gamma[4]*dv10-1.5*jacob_vel_inv0[2]*gamma[3]*dv10+0.8660254037844386*gamma[1]*jacob_vel_inv0[2]*dv10; 
  p0_over_gamma_r[0] = 4.330127018922193*jacob_vel_inv0[1]*gamma[8]*dv10+1.936491673103709*jacob_vel_inv0[0]*gamma[7]*dv10+3.354101966249684*jacob_vel_inv0[1]*gamma[6]*dv10+1.936491673103709*jacob_vel_inv0[1]*gamma[4]*dv10+1.5*jacob_vel_inv0[0]*gamma[3]*dv10+0.8660254037844386*jacob_vel_inv0[0]*gamma[1]*dv10; 
  p0_over_gamma_r[1] = 3.872983346207417*jacob_vel_inv0[2]*gamma[8]*dv10+4.330127018922193*jacob_vel_inv0[0]*gamma[8]*dv10+1.936491673103709*jacob_vel_inv0[1]*gamma[7]*dv10+3.0*jacob_vel_inv0[2]*gamma[6]*dv10+3.354101966249684*jacob_vel_inv0[0]*gamma[6]*dv10+1.732050807568877*jacob_vel_inv0[2]*gamma[4]*dv10+1.936491673103709*jacob_vel_inv0[0]*gamma[4]*dv10+1.5*jacob_vel_inv0[1]*gamma[3]*dv10+0.8660254037844386*jacob_vel_inv0[1]*gamma[1]*dv10; 
  p0_over_gamma_r[2] = 3.872983346207417*jacob_vel_inv0[1]*gamma[8]*dv10+1.936491673103709*jacob_vel_inv0[2]*gamma[7]*dv10+3.0*jacob_vel_inv0[1]*gamma[6]*dv10+1.732050807568877*jacob_vel_inv0[1]*gamma[4]*dv10+1.5*jacob_vel_inv0[2]*gamma[3]*dv10+0.8660254037844386*gamma[1]*jacob_vel_inv0[2]*dv10; 
  const double *B2 = &qmem[15]; 

  double alpha_l[9] = {0.0}; 
  double alpha_r[9] = {0.0}; 

  alpha_l[0] = 1.414213562373095*E1[0]-1.0*B2[0]*p0_over_gamma_l[0]; 
  alpha_l[1] = 1.414213562373095*E1[1]-1.0*p0_over_gamma_l[0]*B2[1]; 
  alpha_l[2] = -1.0*B2[0]*p0_over_gamma_l[1]; 
  alpha_l[3] = -1.0*B2[1]*p0_over_gamma_l[1]; 
  alpha_l[4] = 1.414213562373095*E1[2]-1.0*p0_over_gamma_l[0]*B2[2]; 
  alpha_l[5] = -1.0*B2[0]*p0_over_gamma_l[2]; 
  alpha_l[6] = -1.0*p0_over_gamma_l[1]*B2[2]; 
  alpha_l[7] = -1.0*B2[1]*p0_over_gamma_l[2]; 
  alpha_l[8] = -1.0*B2[2]*p0_over_gamma_l[2]; 

  alpha_r[0] = 1.414213562373095*E1[0]-1.0*B2[0]*p0_over_gamma_r[0]; 
  alpha_r[1] = 1.414213562373095*E1[1]-1.0*p0_over_gamma_r[0]*B2[1]; 
  alpha_r[2] = -1.0*B2[0]*p0_over_gamma_r[1]; 
  alpha_r[3] = -1.0*B2[1]*p0_over_gamma_r[1]; 
  alpha_r[4] = 1.414213562373095*E1[2]-1.0*p0_over_gamma_r[0]*B2[2]; 
  alpha_r[5] = -1.0*B2[0]*p0_over_gamma_r[2]; 
  alpha_r[6] = -1.0*p0_over_gamma_r[1]*B2[2]; 
  alpha_r[7] = -1.0*B2[1]*p0_over_gamma_r[2]; 
  alpha_r[8] = -1.0*B2[2]*p0_over_gamma_r[2]; 

  double fUpwindQuad_l[9] = {0.0};
  double fUpwindQuad_r[9] = {0.0};
  double fUpwind_l[9] = {0.0};;
  double fUpwind_r[9] = {0.0};
  double Ghat_l[9] = {0.0}; 
  double Ghat_r[9] = {0.0}; 

  if (0.4*alpha_l[8]-0.5999999999999995*alpha_l[7]-0.5999999999999999*alpha_l[6]+0.4472135954999579*(alpha_l[5]+alpha_l[4])+0.9*alpha_l[3]-0.6708203932499369*(alpha_l[2]+alpha_l[1])+0.5*alpha_l[0] > 0) { 
    fUpwindQuad_l[0] = tensor_3x_p2_surfx3_eval_quad_node_0_r(fl); 
  } else { 
    fUpwindQuad_l[0] = tensor_3x_p2_surfx3_eval_quad_node_0_l(fc); 
  } 
  if (0.4*alpha_r[8]-0.5999999999999995*alpha_r[7]-0.5999999999999999*alpha_r[6]+0.4472135954999579*(alpha_r[5]+alpha_r[4])+0.9*alpha_r[3]-0.6708203932499369*(alpha_r[2]+alpha_r[1])+0.5*alpha_r[0] > 0) { 
    fUpwindQuad_r[0] = tensor_3x_p2_surfx3_eval_quad_node_0_r(fc); 
  } else { 
    fUpwindQuad_r[0] = tensor_3x_p2_surfx3_eval_quad_node_0_l(fr); 
  } 
  if ((-0.5*alpha_l[8])+0.75*alpha_l[7]-0.5590169943749475*alpha_l[5]+0.4472135954999579*alpha_l[4]-0.6708203932499369*alpha_l[1]+0.5*alpha_l[0] > 0) { 
    fUpwindQuad_l[1] = tensor_3x_p2_surfx3_eval_quad_node_1_r(fl); 
  } else { 
    fUpwindQuad_l[1] = tensor_3x_p2_surfx3_eval_quad_node_1_l(fc); 
  } 
  if ((-0.5*alpha_r[8])+0.75*alpha_r[7]-0.5590169943749475*alpha_r[5]+0.4472135954999579*alpha_r[4]-0.6708203932499369*alpha_r[1]+0.5*alpha_r[0] > 0) { 
    fUpwindQuad_r[1] = tensor_3x_p2_surfx3_eval_quad_node_1_r(fc); 
  } else { 
    fUpwindQuad_r[1] = tensor_3x_p2_surfx3_eval_quad_node_1_l(fr); 
  } 
  if (0.4*alpha_l[8]-0.5999999999999995*alpha_l[7]+0.5999999999999999*alpha_l[6]+0.4472135954999579*(alpha_l[5]+alpha_l[4])-0.9*alpha_l[3]+0.6708203932499369*alpha_l[2]-0.6708203932499369*alpha_l[1]+0.5*alpha_l[0] > 0) { 
    fUpwindQuad_l[2] = tensor_3x_p2_surfx3_eval_quad_node_2_r(fl); 
  } else { 
    fUpwindQuad_l[2] = tensor_3x_p2_surfx3_eval_quad_node_2_l(fc); 
  } 
  if (0.4*alpha_r[8]-0.5999999999999995*alpha_r[7]+0.5999999999999999*alpha_r[6]+0.4472135954999579*(alpha_r[5]+alpha_r[4])-0.9*alpha_r[3]+0.6708203932499369*alpha_r[2]-0.6708203932499369*alpha_r[1]+0.5*alpha_r[0] > 0) { 
    fUpwindQuad_r[2] = tensor_3x_p2_surfx3_eval_quad_node_2_r(fc); 
  } else { 
    fUpwindQuad_r[2] = tensor_3x_p2_surfx3_eval_quad_node_2_l(fr); 
  } 
  if ((-0.5*alpha_l[8])+0.75*alpha_l[6]+0.4472135954999579*alpha_l[5]-0.5590169943749475*alpha_l[4]-0.6708203932499369*alpha_l[2]+0.5*alpha_l[0] > 0) { 
    fUpwindQuad_l[3] = tensor_3x_p2_surfx3_eval_quad_node_3_r(fl); 
  } else { 
    fUpwindQuad_l[3] = tensor_3x_p2_surfx3_eval_quad_node_3_l(fc); 
  } 
  if ((-0.5*alpha_r[8])+0.75*alpha_r[6]+0.4472135954999579*alpha_r[5]-0.5590169943749475*alpha_r[4]-0.6708203932499369*alpha_r[2]+0.5*alpha_r[0] > 0) { 
    fUpwindQuad_r[3] = tensor_3x_p2_surfx3_eval_quad_node_3_r(fc); 
  } else { 
    fUpwindQuad_r[3] = tensor_3x_p2_surfx3_eval_quad_node_3_l(fr); 
  } 
  if (0.625*alpha_l[8]-0.5590169943749475*(alpha_l[5]+alpha_l[4])+0.5*alpha_l[0] > 0) { 
    fUpwindQuad_l[4] = tensor_3x_p2_surfx3_eval_quad_node_4_r(fl); 
  } else { 
    fUpwindQuad_l[4] = tensor_3x_p2_surfx3_eval_quad_node_4_l(fc); 
  } 
  if (0.625*alpha_r[8]-0.5590169943749475*(alpha_r[5]+alpha_r[4])+0.5*alpha_r[0] > 0) { 
    fUpwindQuad_r[4] = tensor_3x_p2_surfx3_eval_quad_node_4_r(fc); 
  } else { 
    fUpwindQuad_r[4] = tensor_3x_p2_surfx3_eval_quad_node_4_l(fr); 
  } 
  if ((-0.5*alpha_l[8])-0.75*alpha_l[6]+0.4472135954999579*alpha_l[5]-0.5590169943749475*alpha_l[4]+0.6708203932499369*alpha_l[2]+0.5*alpha_l[0] > 0) { 
    fUpwindQuad_l[5] = tensor_3x_p2_surfx3_eval_quad_node_5_r(fl); 
  } else { 
    fUpwindQuad_l[5] = tensor_3x_p2_surfx3_eval_quad_node_5_l(fc); 
  } 
  if ((-0.5*alpha_r[8])-0.75*alpha_r[6]+0.4472135954999579*alpha_r[5]-0.5590169943749475*alpha_r[4]+0.6708203932499369*alpha_r[2]+0.5*alpha_r[0] > 0) { 
    fUpwindQuad_r[5] = tensor_3x_p2_surfx3_eval_quad_node_5_r(fc); 
  } else { 
    fUpwindQuad_r[5] = tensor_3x_p2_surfx3_eval_quad_node_5_l(fr); 
  } 
  if (0.4*alpha_l[8]+0.5999999999999995*alpha_l[7]-0.5999999999999999*alpha_l[6]+0.4472135954999579*(alpha_l[5]+alpha_l[4])-0.9*alpha_l[3]-0.6708203932499369*alpha_l[2]+0.6708203932499369*alpha_l[1]+0.5*alpha_l[0] > 0) { 
    fUpwindQuad_l[6] = tensor_3x_p2_surfx3_eval_quad_node_6_r(fl); 
  } else { 
    fUpwindQuad_l[6] = tensor_3x_p2_surfx3_eval_quad_node_6_l(fc); 
  } 
  if (0.4*alpha_r[8]+0.5999999999999995*alpha_r[7]-0.5999999999999999*alpha_r[6]+0.4472135954999579*(alpha_r[5]+alpha_r[4])-0.9*alpha_r[3]-0.6708203932499369*alpha_r[2]+0.6708203932499369*alpha_r[1]+0.5*alpha_r[0] > 0) { 
    fUpwindQuad_r[6] = tensor_3x_p2_surfx3_eval_quad_node_6_r(fc); 
  } else { 
    fUpwindQuad_r[6] = tensor_3x_p2_surfx3_eval_quad_node_6_l(fr); 
  } 
  if ((-0.5*alpha_l[8])-0.75*alpha_l[7]-0.5590169943749475*alpha_l[5]+0.4472135954999579*alpha_l[4]+0.6708203932499369*alpha_l[1]+0.5*alpha_l[0] > 0) { 
    fUpwindQuad_l[7] = tensor_3x_p2_surfx3_eval_quad_node_7_r(fl); 
  } else { 
    fUpwindQuad_l[7] = tensor_3x_p2_surfx3_eval_quad_node_7_l(fc); 
  } 
  if ((-0.5*alpha_r[8])-0.75*alpha_r[7]-0.5590169943749475*alpha_r[5]+0.4472135954999579*alpha_r[4]+0.6708203932499369*alpha_r[1]+0.5*alpha_r[0] > 0) { 
    fUpwindQuad_r[7] = tensor_3x_p2_surfx3_eval_quad_node_7_r(fc); 
  } else { 
    fUpwindQuad_r[7] = tensor_3x_p2_surfx3_eval_quad_node_7_l(fr); 
  } 
  if (0.4*alpha_l[8]+0.5999999999999995*alpha_l[7]+0.5999999999999999*alpha_l[6]+0.4472135954999579*(alpha_l[5]+alpha_l[4])+0.9*alpha_l[3]+0.6708203932499369*(alpha_l[2]+alpha_l[1])+0.5*alpha_l[0] > 0) { 
    fUpwindQuad_l[8] = tensor_3x_p2_surfx3_eval_quad_node_8_r(fl); 
  } else { 
    fUpwindQuad_l[8] = tensor_3x_p2_surfx3_eval_quad_node_8_l(fc); 
  } 
  if (0.4*alpha_r[8]+0.5999999999999995*alpha_r[7]+0.5999999999999999*alpha_r[6]+0.4472135954999579*(alpha_r[5]+alpha_r[4])+0.9*alpha_r[3]+0.6708203932499369*(alpha_r[2]+alpha_r[1])+0.5*alpha_r[0] > 0) { 
    fUpwindQuad_r[8] = tensor_3x_p2_surfx3_eval_quad_node_8_r(fc); 
  } else { 
    fUpwindQuad_r[8] = tensor_3x_p2_surfx3_eval_quad_node_8_l(fr); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  tensor_3x_p2_upwind_quad_to_modal(fUpwindQuad_l, fUpwind_l); 
  tensor_3x_p2_upwind_quad_to_modal(fUpwindQuad_r, fUpwind_r); 
  Ghat_l[0] = 0.5*alpha_l[8]*fUpwind_l[8]+0.5*alpha_l[7]*fUpwind_l[7]+0.5*alpha_l[6]*fUpwind_l[6]+0.5*alpha_l[5]*fUpwind_l[5]+0.5*alpha_l[4]*fUpwind_l[4]+0.5*alpha_l[3]*fUpwind_l[3]+0.5*alpha_l[2]*fUpwind_l[2]+0.5*alpha_l[1]*fUpwind_l[1]+0.5*alpha_l[0]*fUpwind_l[0]; 
  Ghat_l[1] = 0.447213595499958*alpha_l[7]*fUpwind_l[8]+0.447213595499958*fUpwind_l[7]*alpha_l[8]+0.5000000000000001*alpha_l[5]*fUpwind_l[7]+0.5000000000000001*fUpwind_l[5]*alpha_l[7]+0.447213595499958*alpha_l[3]*fUpwind_l[6]+0.447213595499958*fUpwind_l[3]*alpha_l[6]+0.4472135954999579*alpha_l[1]*fUpwind_l[4]+0.4472135954999579*fUpwind_l[1]*alpha_l[4]+0.5*alpha_l[2]*fUpwind_l[3]+0.5*fUpwind_l[2]*alpha_l[3]+0.5*alpha_l[0]*fUpwind_l[1]+0.5*fUpwind_l[0]*alpha_l[1]; 
  Ghat_l[2] = 0.447213595499958*alpha_l[6]*fUpwind_l[8]+0.447213595499958*fUpwind_l[6]*alpha_l[8]+0.447213595499958*alpha_l[3]*fUpwind_l[7]+0.447213595499958*fUpwind_l[3]*alpha_l[7]+0.5000000000000001*alpha_l[4]*fUpwind_l[6]+0.5000000000000001*fUpwind_l[4]*alpha_l[6]+0.4472135954999579*alpha_l[2]*fUpwind_l[5]+0.4472135954999579*fUpwind_l[2]*alpha_l[5]+0.5*alpha_l[1]*fUpwind_l[3]+0.5*fUpwind_l[1]*alpha_l[3]+0.5*alpha_l[0]*fUpwind_l[2]+0.5*fUpwind_l[0]*alpha_l[2]; 
  Ghat_l[3] = 0.4*alpha_l[3]*fUpwind_l[8]+0.4*fUpwind_l[3]*alpha_l[8]+0.4*alpha_l[6]*fUpwind_l[7]+0.447213595499958*alpha_l[2]*fUpwind_l[7]+0.4*fUpwind_l[6]*alpha_l[7]+0.447213595499958*fUpwind_l[2]*alpha_l[7]+0.447213595499958*alpha_l[1]*fUpwind_l[6]+0.447213595499958*fUpwind_l[1]*alpha_l[6]+0.4472135954999579*alpha_l[3]*fUpwind_l[5]+0.4472135954999579*fUpwind_l[3]*alpha_l[5]+0.4472135954999579*alpha_l[3]*fUpwind_l[4]+0.4472135954999579*fUpwind_l[3]*alpha_l[4]+0.5*alpha_l[0]*fUpwind_l[3]+0.5*fUpwind_l[0]*alpha_l[3]+0.5*alpha_l[1]*fUpwind_l[2]+0.5*fUpwind_l[1]*alpha_l[2]; 
  Ghat_l[4] = 0.31943828249997*alpha_l[8]*fUpwind_l[8]+0.5*alpha_l[5]*fUpwind_l[8]+0.5*fUpwind_l[5]*alpha_l[8]+0.4472135954999579*alpha_l[7]*fUpwind_l[7]+0.31943828249997*alpha_l[6]*fUpwind_l[6]+0.5000000000000001*alpha_l[2]*fUpwind_l[6]+0.5000000000000001*fUpwind_l[2]*alpha_l[6]+0.31943828249997*alpha_l[4]*fUpwind_l[4]+0.5*alpha_l[0]*fUpwind_l[4]+0.5*fUpwind_l[0]*alpha_l[4]+0.4472135954999579*alpha_l[3]*fUpwind_l[3]+0.4472135954999579*alpha_l[1]*fUpwind_l[1]; 
  Ghat_l[5] = 0.31943828249997*alpha_l[8]*fUpwind_l[8]+0.5*alpha_l[4]*fUpwind_l[8]+0.5*fUpwind_l[4]*alpha_l[8]+0.31943828249997*alpha_l[7]*fUpwind_l[7]+0.5000000000000001*alpha_l[1]*fUpwind_l[7]+0.5000000000000001*fUpwind_l[1]*alpha_l[7]+0.4472135954999579*alpha_l[6]*fUpwind_l[6]+0.31943828249997*alpha_l[5]*fUpwind_l[5]+0.5*alpha_l[0]*fUpwind_l[5]+0.5*fUpwind_l[0]*alpha_l[5]+0.4472135954999579*alpha_l[3]*fUpwind_l[3]+0.4472135954999579*alpha_l[2]*fUpwind_l[2]; 
  Ghat_l[6] = 0.2857142857142857*alpha_l[6]*fUpwind_l[8]+0.447213595499958*alpha_l[2]*fUpwind_l[8]+0.2857142857142857*fUpwind_l[6]*alpha_l[8]+0.447213595499958*fUpwind_l[2]*alpha_l[8]+0.4*alpha_l[3]*fUpwind_l[7]+0.4*fUpwind_l[3]*alpha_l[7]+0.4472135954999579*alpha_l[5]*fUpwind_l[6]+0.31943828249997*alpha_l[4]*fUpwind_l[6]+0.5*alpha_l[0]*fUpwind_l[6]+0.4472135954999579*fUpwind_l[5]*alpha_l[6]+0.31943828249997*fUpwind_l[4]*alpha_l[6]+0.5*fUpwind_l[0]*alpha_l[6]+0.5000000000000001*alpha_l[2]*fUpwind_l[4]+0.5000000000000001*fUpwind_l[2]*alpha_l[4]+0.447213595499958*alpha_l[1]*fUpwind_l[3]+0.447213595499958*fUpwind_l[1]*alpha_l[3]; 
  Ghat_l[7] = 0.2857142857142857*alpha_l[7]*fUpwind_l[8]+0.447213595499958*alpha_l[1]*fUpwind_l[8]+0.2857142857142857*fUpwind_l[7]*alpha_l[8]+0.447213595499958*fUpwind_l[1]*alpha_l[8]+0.31943828249997*alpha_l[5]*fUpwind_l[7]+0.4472135954999579*alpha_l[4]*fUpwind_l[7]+0.5*alpha_l[0]*fUpwind_l[7]+0.31943828249997*fUpwind_l[5]*alpha_l[7]+0.4472135954999579*fUpwind_l[4]*alpha_l[7]+0.5*fUpwind_l[0]*alpha_l[7]+0.4*alpha_l[3]*fUpwind_l[6]+0.4*fUpwind_l[3]*alpha_l[6]+0.5000000000000001*alpha_l[1]*fUpwind_l[5]+0.5000000000000001*fUpwind_l[1]*alpha_l[5]+0.447213595499958*alpha_l[2]*fUpwind_l[3]+0.447213595499958*fUpwind_l[2]*alpha_l[3]; 
  Ghat_l[8] = 0.2040816326530612*alpha_l[8]*fUpwind_l[8]+0.31943828249997*alpha_l[5]*fUpwind_l[8]+0.31943828249997*alpha_l[4]*fUpwind_l[8]+0.5*alpha_l[0]*fUpwind_l[8]+0.31943828249997*fUpwind_l[5]*alpha_l[8]+0.31943828249997*fUpwind_l[4]*alpha_l[8]+0.5*fUpwind_l[0]*alpha_l[8]+0.2857142857142857*alpha_l[7]*fUpwind_l[7]+0.447213595499958*alpha_l[1]*fUpwind_l[7]+0.447213595499958*fUpwind_l[1]*alpha_l[7]+0.2857142857142857*alpha_l[6]*fUpwind_l[6]+0.447213595499958*alpha_l[2]*fUpwind_l[6]+0.447213595499958*fUpwind_l[2]*alpha_l[6]+0.5*alpha_l[4]*fUpwind_l[5]+0.5*fUpwind_l[4]*alpha_l[5]+0.4*alpha_l[3]*fUpwind_l[3]; 

  Ghat_r[0] = 0.5*alpha_r[8]*fUpwind_r[8]+0.5*alpha_r[7]*fUpwind_r[7]+0.5*alpha_r[6]*fUpwind_r[6]+0.5*alpha_r[5]*fUpwind_r[5]+0.5*alpha_r[4]*fUpwind_r[4]+0.5*alpha_r[3]*fUpwind_r[3]+0.5*alpha_r[2]*fUpwind_r[2]+0.5*alpha_r[1]*fUpwind_r[1]+0.5*alpha_r[0]*fUpwind_r[0]; 
  Ghat_r[1] = 0.447213595499958*alpha_r[7]*fUpwind_r[8]+0.447213595499958*fUpwind_r[7]*alpha_r[8]+0.5000000000000001*alpha_r[5]*fUpwind_r[7]+0.5000000000000001*fUpwind_r[5]*alpha_r[7]+0.447213595499958*alpha_r[3]*fUpwind_r[6]+0.447213595499958*fUpwind_r[3]*alpha_r[6]+0.4472135954999579*alpha_r[1]*fUpwind_r[4]+0.4472135954999579*fUpwind_r[1]*alpha_r[4]+0.5*alpha_r[2]*fUpwind_r[3]+0.5*fUpwind_r[2]*alpha_r[3]+0.5*alpha_r[0]*fUpwind_r[1]+0.5*fUpwind_r[0]*alpha_r[1]; 
  Ghat_r[2] = 0.447213595499958*alpha_r[6]*fUpwind_r[8]+0.447213595499958*fUpwind_r[6]*alpha_r[8]+0.447213595499958*alpha_r[3]*fUpwind_r[7]+0.447213595499958*fUpwind_r[3]*alpha_r[7]+0.5000000000000001*alpha_r[4]*fUpwind_r[6]+0.5000000000000001*fUpwind_r[4]*alpha_r[6]+0.4472135954999579*alpha_r[2]*fUpwind_r[5]+0.4472135954999579*fUpwind_r[2]*alpha_r[5]+0.5*alpha_r[1]*fUpwind_r[3]+0.5*fUpwind_r[1]*alpha_r[3]+0.5*alpha_r[0]*fUpwind_r[2]+0.5*fUpwind_r[0]*alpha_r[2]; 
  Ghat_r[3] = 0.4*alpha_r[3]*fUpwind_r[8]+0.4*fUpwind_r[3]*alpha_r[8]+0.4*alpha_r[6]*fUpwind_r[7]+0.447213595499958*alpha_r[2]*fUpwind_r[7]+0.4*fUpwind_r[6]*alpha_r[7]+0.447213595499958*fUpwind_r[2]*alpha_r[7]+0.447213595499958*alpha_r[1]*fUpwind_r[6]+0.447213595499958*fUpwind_r[1]*alpha_r[6]+0.4472135954999579*alpha_r[3]*fUpwind_r[5]+0.4472135954999579*fUpwind_r[3]*alpha_r[5]+0.4472135954999579*alpha_r[3]*fUpwind_r[4]+0.4472135954999579*fUpwind_r[3]*alpha_r[4]+0.5*alpha_r[0]*fUpwind_r[3]+0.5*fUpwind_r[0]*alpha_r[3]+0.5*alpha_r[1]*fUpwind_r[2]+0.5*fUpwind_r[1]*alpha_r[2]; 
  Ghat_r[4] = 0.31943828249997*alpha_r[8]*fUpwind_r[8]+0.5*alpha_r[5]*fUpwind_r[8]+0.5*fUpwind_r[5]*alpha_r[8]+0.4472135954999579*alpha_r[7]*fUpwind_r[7]+0.31943828249997*alpha_r[6]*fUpwind_r[6]+0.5000000000000001*alpha_r[2]*fUpwind_r[6]+0.5000000000000001*fUpwind_r[2]*alpha_r[6]+0.31943828249997*alpha_r[4]*fUpwind_r[4]+0.5*alpha_r[0]*fUpwind_r[4]+0.5*fUpwind_r[0]*alpha_r[4]+0.4472135954999579*alpha_r[3]*fUpwind_r[3]+0.4472135954999579*alpha_r[1]*fUpwind_r[1]; 
  Ghat_r[5] = 0.31943828249997*alpha_r[8]*fUpwind_r[8]+0.5*alpha_r[4]*fUpwind_r[8]+0.5*fUpwind_r[4]*alpha_r[8]+0.31943828249997*alpha_r[7]*fUpwind_r[7]+0.5000000000000001*alpha_r[1]*fUpwind_r[7]+0.5000000000000001*fUpwind_r[1]*alpha_r[7]+0.4472135954999579*alpha_r[6]*fUpwind_r[6]+0.31943828249997*alpha_r[5]*fUpwind_r[5]+0.5*alpha_r[0]*fUpwind_r[5]+0.5*fUpwind_r[0]*alpha_r[5]+0.4472135954999579*alpha_r[3]*fUpwind_r[3]+0.4472135954999579*alpha_r[2]*fUpwind_r[2]; 
  Ghat_r[6] = 0.2857142857142857*alpha_r[6]*fUpwind_r[8]+0.447213595499958*alpha_r[2]*fUpwind_r[8]+0.2857142857142857*fUpwind_r[6]*alpha_r[8]+0.447213595499958*fUpwind_r[2]*alpha_r[8]+0.4*alpha_r[3]*fUpwind_r[7]+0.4*fUpwind_r[3]*alpha_r[7]+0.4472135954999579*alpha_r[5]*fUpwind_r[6]+0.31943828249997*alpha_r[4]*fUpwind_r[6]+0.5*alpha_r[0]*fUpwind_r[6]+0.4472135954999579*fUpwind_r[5]*alpha_r[6]+0.31943828249997*fUpwind_r[4]*alpha_r[6]+0.5*fUpwind_r[0]*alpha_r[6]+0.5000000000000001*alpha_r[2]*fUpwind_r[4]+0.5000000000000001*fUpwind_r[2]*alpha_r[4]+0.447213595499958*alpha_r[1]*fUpwind_r[3]+0.447213595499958*fUpwind_r[1]*alpha_r[3]; 
  Ghat_r[7] = 0.2857142857142857*alpha_r[7]*fUpwind_r[8]+0.447213595499958*alpha_r[1]*fUpwind_r[8]+0.2857142857142857*fUpwind_r[7]*alpha_r[8]+0.447213595499958*fUpwind_r[1]*alpha_r[8]+0.31943828249997*alpha_r[5]*fUpwind_r[7]+0.4472135954999579*alpha_r[4]*fUpwind_r[7]+0.5*alpha_r[0]*fUpwind_r[7]+0.31943828249997*fUpwind_r[5]*alpha_r[7]+0.4472135954999579*fUpwind_r[4]*alpha_r[7]+0.5*fUpwind_r[0]*alpha_r[7]+0.4*alpha_r[3]*fUpwind_r[6]+0.4*fUpwind_r[3]*alpha_r[6]+0.5000000000000001*alpha_r[1]*fUpwind_r[5]+0.5000000000000001*fUpwind_r[1]*alpha_r[5]+0.447213595499958*alpha_r[2]*fUpwind_r[3]+0.447213595499958*fUpwind_r[2]*alpha_r[3]; 
  Ghat_r[8] = 0.2040816326530612*alpha_r[8]*fUpwind_r[8]+0.31943828249997*alpha_r[5]*fUpwind_r[8]+0.31943828249997*alpha_r[4]*fUpwind_r[8]+0.5*alpha_r[0]*fUpwind_r[8]+0.31943828249997*fUpwind_r[5]*alpha_r[8]+0.31943828249997*fUpwind_r[4]*alpha_r[8]+0.5*fUpwind_r[0]*alpha_r[8]+0.2857142857142857*alpha_r[7]*fUpwind_r[7]+0.447213595499958*alpha_r[1]*fUpwind_r[7]+0.447213595499958*fUpwind_r[1]*alpha_r[7]+0.2857142857142857*alpha_r[6]*fUpwind_r[6]+0.447213595499958*alpha_r[2]*fUpwind_r[6]+0.447213595499958*fUpwind_r[2]*alpha_r[6]+0.5*alpha_r[4]*fUpwind_r[5]+0.5*fUpwind_r[4]*alpha_r[5]+0.4*alpha_r[3]*fUpwind_r[3]; 

  out[0] += ((1.118033988749895*Ghat_l[0]-1.118033988749895*Ghat_r[0])*jacob_vel_inv1[2]-0.8660254037844386*(Ghat_r[0]+Ghat_l[0])*jacob_vel_inv1[1]+(0.5*Ghat_l[0]-0.5*Ghat_r[0])*jacob_vel_inv1[0])*dv11; 
  out[1] += ((1.118033988749895*Ghat_l[1]-1.118033988749895*Ghat_r[1])*jacob_vel_inv1[2]-0.8660254037844386*(Ghat_r[1]+Ghat_l[1])*jacob_vel_inv1[1]+jacob_vel_inv1[0]*(0.5*Ghat_l[1]-0.5*Ghat_r[1]))*dv11; 
  out[2] += ((1.118033988749895*Ghat_l[2]-1.118033988749895*Ghat_r[2])*jacob_vel_inv1[2]+((-0.8660254037844386*jacob_vel_inv1[1])-0.5*jacob_vel_inv1[0])*Ghat_r[2]+(0.5*jacob_vel_inv1[0]-0.8660254037844386*jacob_vel_inv1[1])*Ghat_l[2])*dv11; 
  out[3] += ((-1.936491673103709*(Ghat_r[0]+Ghat_l[0])*jacob_vel_inv1[2])+(1.5*Ghat_l[0]-1.5*Ghat_r[0])*jacob_vel_inv1[1]-0.8660254037844386*(Ghat_r[0]+Ghat_l[0])*jacob_vel_inv1[0])*dv11; 
  out[4] += (((-1.118033988749895*jacob_vel_inv1[2])-0.8660254037844386*jacob_vel_inv1[1]-0.5*jacob_vel_inv1[0])*Ghat_r[3]+(1.118033988749895*jacob_vel_inv1[2]-0.8660254037844386*jacob_vel_inv1[1]+0.5*jacob_vel_inv1[0])*Ghat_l[3])*dv11; 
  out[5] += ((-1.936491673103709*(Ghat_r[1]+Ghat_l[1])*jacob_vel_inv1[2])+(1.5*Ghat_l[1]-1.5*Ghat_r[1])*jacob_vel_inv1[1]-0.8660254037844386*jacob_vel_inv1[0]*(Ghat_r[1]+Ghat_l[1]))*dv11; 
  out[6] += ((-1.936491673103709*(Ghat_r[2]+Ghat_l[2])*jacob_vel_inv1[2])+((-1.5*jacob_vel_inv1[1])-0.8660254037844386*jacob_vel_inv1[0])*Ghat_r[2]+(1.5*jacob_vel_inv1[1]-0.8660254037844386*jacob_vel_inv1[0])*Ghat_l[2])*dv11; 
  out[7] += (((-1.118033988749895*jacob_vel_inv1[2])-0.8660254037844386*jacob_vel_inv1[1]-0.5*jacob_vel_inv1[0])*Ghat_r[4]+(1.118033988749895*jacob_vel_inv1[2]-0.8660254037844386*jacob_vel_inv1[1]+0.5*jacob_vel_inv1[0])*Ghat_l[4])*dv11; 
  out[8] += (((-1.118033988749895*jacob_vel_inv1[2])-0.8660254037844386*jacob_vel_inv1[1]-0.5*jacob_vel_inv1[0])*Ghat_r[5]+(1.118033988749895*jacob_vel_inv1[2]-0.8660254037844386*jacob_vel_inv1[1]+0.5*jacob_vel_inv1[0])*Ghat_l[5])*dv11; 
  out[9] += ((2.5*Ghat_l[0]-2.5*Ghat_r[0])*jacob_vel_inv1[2]-1.936491673103709*(Ghat_r[0]+Ghat_l[0])*jacob_vel_inv1[1]+(1.118033988749895*Ghat_l[0]-1.118033988749895*Ghat_r[0])*jacob_vel_inv1[0])*dv11; 
  out[10] += (((-1.936491673103709*jacob_vel_inv1[2])-1.5*jacob_vel_inv1[1]-0.8660254037844386*jacob_vel_inv1[0])*Ghat_r[3]+((-1.936491673103709*jacob_vel_inv1[2])+1.5*jacob_vel_inv1[1]-0.8660254037844386*jacob_vel_inv1[0])*Ghat_l[3])*dv11; 
  out[11] += (((-1.118033988749895*jacob_vel_inv1[2])-0.8660254037844386*jacob_vel_inv1[1]-0.5*jacob_vel_inv1[0])*Ghat_r[6]+(1.118033988749895*jacob_vel_inv1[2]-0.8660254037844386*jacob_vel_inv1[1]+0.5*jacob_vel_inv1[0])*Ghat_l[6])*dv11; 
  out[12] += (((-1.118033988749895*jacob_vel_inv1[2])-0.8660254037844386*jacob_vel_inv1[1]-0.5*jacob_vel_inv1[0])*Ghat_r[7]+(1.118033988749895*jacob_vel_inv1[2]-0.8660254037844386*jacob_vel_inv1[1]+0.5*jacob_vel_inv1[0])*Ghat_l[7])*dv11; 
  out[13] += (((-1.936491673103709*jacob_vel_inv1[2])-1.5*jacob_vel_inv1[1]-0.8660254037844387*jacob_vel_inv1[0])*Ghat_r[4]+((-1.936491673103709*jacob_vel_inv1[2])+1.5*jacob_vel_inv1[1]-0.8660254037844387*jacob_vel_inv1[0])*Ghat_l[4])*dv11; 
  out[14] += (((-1.936491673103709*jacob_vel_inv1[2])-1.5*jacob_vel_inv1[1]-0.8660254037844387*jacob_vel_inv1[0])*Ghat_r[5]+((-1.936491673103709*jacob_vel_inv1[2])+1.5*jacob_vel_inv1[1]-0.8660254037844387*jacob_vel_inv1[0])*Ghat_l[5])*dv11; 
  out[15] += ((2.5*Ghat_l[1]-2.5*Ghat_r[1])*jacob_vel_inv1[2]-1.936491673103709*(Ghat_r[1]+Ghat_l[1])*jacob_vel_inv1[1]+jacob_vel_inv1[0]*(1.118033988749895*Ghat_l[1]-1.118033988749895*Ghat_r[1]))*dv11; 
  out[16] += ((2.5*Ghat_l[2]-2.5*Ghat_r[2])*jacob_vel_inv1[2]+((-1.936491673103709*jacob_vel_inv1[1])-1.118033988749895*jacob_vel_inv1[0])*Ghat_r[2]+(1.118033988749895*jacob_vel_inv1[0]-1.936491673103709*jacob_vel_inv1[1])*Ghat_l[2])*dv11; 
  out[17] += (((-1.936491673103709*jacob_vel_inv1[2])-1.5*jacob_vel_inv1[1]-0.8660254037844387*jacob_vel_inv1[0])*Ghat_r[6]+((-1.936491673103709*jacob_vel_inv1[2])+1.5*jacob_vel_inv1[1]-0.8660254037844387*jacob_vel_inv1[0])*Ghat_l[6])*dv11; 
  out[18] += (((-1.936491673103709*jacob_vel_inv1[2])-1.5*jacob_vel_inv1[1]-0.8660254037844387*jacob_vel_inv1[0])*Ghat_r[7]+((-1.936491673103709*jacob_vel_inv1[2])+1.5*jacob_vel_inv1[1]-0.8660254037844387*jacob_vel_inv1[0])*Ghat_l[7])*dv11; 
  out[19] += (((-2.5*jacob_vel_inv1[2])-1.936491673103709*jacob_vel_inv1[1]-1.118033988749895*jacob_vel_inv1[0])*Ghat_r[3]+(2.5*jacob_vel_inv1[2]-1.936491673103709*jacob_vel_inv1[1]+1.118033988749895*jacob_vel_inv1[0])*Ghat_l[3])*dv11; 
  out[20] += (((-1.118033988749895*jacob_vel_inv1[2])-0.8660254037844386*jacob_vel_inv1[1]-0.5*jacob_vel_inv1[0])*Ghat_r[8]+(1.118033988749895*jacob_vel_inv1[2]-0.8660254037844386*jacob_vel_inv1[1]+0.5*jacob_vel_inv1[0])*Ghat_l[8])*dv11; 
  out[21] += (((-2.5*jacob_vel_inv1[2])-1.936491673103709*jacob_vel_inv1[1]-1.118033988749895*jacob_vel_inv1[0])*Ghat_r[4]+(2.5*jacob_vel_inv1[2]-1.936491673103709*jacob_vel_inv1[1]+1.118033988749895*jacob_vel_inv1[0])*Ghat_l[4])*dv11; 
  out[22] += (((-2.5*jacob_vel_inv1[2])-1.936491673103709*jacob_vel_inv1[1]-1.118033988749895*jacob_vel_inv1[0])*Ghat_r[5]+(2.5*jacob_vel_inv1[2]-1.936491673103709*jacob_vel_inv1[1]+1.118033988749895*jacob_vel_inv1[0])*Ghat_l[5])*dv11; 
  out[23] += (((-1.936491673103709*jacob_vel_inv1[2])-1.5*jacob_vel_inv1[1]-0.8660254037844386*jacob_vel_inv1[0])*Ghat_r[8]+((-1.936491673103709*jacob_vel_inv1[2])+1.5*jacob_vel_inv1[1]-0.8660254037844386*jacob_vel_inv1[0])*Ghat_l[8])*dv11; 
  out[24] += (((-2.5*jacob_vel_inv1[2])-1.936491673103709*jacob_vel_inv1[1]-1.118033988749895*jacob_vel_inv1[0])*Ghat_r[6]+(2.5*jacob_vel_inv1[2]-1.936491673103709*jacob_vel_inv1[1]+1.118033988749895*jacob_vel_inv1[0])*Ghat_l[6])*dv11; 
  out[25] += (((-2.5*jacob_vel_inv1[2])-1.936491673103709*jacob_vel_inv1[1]-1.118033988749895*jacob_vel_inv1[0])*Ghat_r[7]+(2.5*jacob_vel_inv1[2]-1.936491673103709*jacob_vel_inv1[1]+1.118033988749895*jacob_vel_inv1[0])*Ghat_l[7])*dv11; 
  out[26] += (((-2.5*jacob_vel_inv1[2])-1.936491673103709*jacob_vel_inv1[1]-1.118033988749895*jacob_vel_inv1[0])*Ghat_r[8]+(2.5*jacob_vel_inv1[2]-1.936491673103709*jacob_vel_inv1[1]+1.118033988749895*jacob_vel_inv1[0])*Ghat_l[8])*dv11; 

  return 0.;

} 
