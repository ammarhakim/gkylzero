#include <gkyl_lbo_vlasov_kernels.h> 
#include <gkyl_basis_hyb_2x2v_p1_surfx4_eval_quad.h> 
#include <gkyl_basis_hyb_2x2v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double lbo_vlasov_vmap_drag_surfvy_2x2v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_vel_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[4]: cell-center coordinates. 
  // dxv[4]: cell spacing. 
  // vmap: Velocity-space nonuniform mapping in each dimension. 
  // jacob_vel_inv: Inverse of velocity space Jacobian in each dimension. 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[12]: sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 
  // fl/fc/fr: distribution function in cells 
  // out: incremented distribution function in cell 
  double rdv2 = 2.0/dxv[3]; 

  const double *v1 = &vmap[4]; 
  const double *jacob_vel_inv1 = &jacob_vel_inv[3]; 
  const double *sumNuUy = &nuPrimMomsSum[4]; 

  double alphaDrSurf_l[12] = {0.0}; 
  alphaDrSurf_l[0] = (-4.183300132670378*nuSum[0]*jacob_vel_inv1[2]*v1[3])+3.24037034920393*nuSum[0]*jacob_vel_inv1[1]*v1[3]-1.870828693386971*jacob_vel_inv1[0]*nuSum[0]*v1[3]+3.535533905932737*nuSum[0]*jacob_vel_inv1[2]*v1[2]-2.738612787525831*nuSum[0]*jacob_vel_inv1[1]*v1[2]+1.58113883008419*jacob_vel_inv1[0]*nuSum[0]*v1[2]-2.738612787525831*nuSum[0]*v1[1]*jacob_vel_inv1[2]+1.58113883008419*nuSum[0]*v1[0]*jacob_vel_inv1[2]-2.23606797749979*sumNuUy[0]*jacob_vel_inv1[2]+2.121320343559642*nuSum[0]*jacob_vel_inv1[1]*v1[1]-1.224744871391589*jacob_vel_inv1[0]*nuSum[0]*v1[1]-1.224744871391589*nuSum[0]*v1[0]*jacob_vel_inv1[1]+1.732050807568877*sumNuUy[0]*jacob_vel_inv1[1]+0.7071067811865475*jacob_vel_inv1[0]*nuSum[0]*v1[0]-1.0*jacob_vel_inv1[0]*sumNuUy[0]; 
  alphaDrSurf_l[1] = (-4.183300132670378*nuSum[1]*jacob_vel_inv1[2]*v1[3])+3.24037034920393*jacob_vel_inv1[1]*nuSum[1]*v1[3]-1.870828693386971*jacob_vel_inv1[0]*nuSum[1]*v1[3]+3.535533905932737*nuSum[1]*jacob_vel_inv1[2]*v1[2]-2.738612787525831*jacob_vel_inv1[1]*nuSum[1]*v1[2]+1.58113883008419*jacob_vel_inv1[0]*nuSum[1]*v1[2]-2.738612787525831*nuSum[1]*v1[1]*jacob_vel_inv1[2]-2.23606797749979*sumNuUy[1]*jacob_vel_inv1[2]+1.58113883008419*v1[0]*nuSum[1]*jacob_vel_inv1[2]+2.121320343559642*jacob_vel_inv1[1]*nuSum[1]*v1[1]-1.224744871391589*jacob_vel_inv1[0]*nuSum[1]*v1[1]+1.732050807568877*jacob_vel_inv1[1]*sumNuUy[1]-1.0*jacob_vel_inv1[0]*sumNuUy[1]-1.224744871391589*v1[0]*jacob_vel_inv1[1]*nuSum[1]+0.7071067811865475*jacob_vel_inv1[0]*v1[0]*nuSum[1]; 
  alphaDrSurf_l[2] = (-4.183300132670378*jacob_vel_inv1[2]*nuSum[2]*v1[3])+3.24037034920393*jacob_vel_inv1[1]*nuSum[2]*v1[3]-1.870828693386971*jacob_vel_inv1[0]*nuSum[2]*v1[3]+3.535533905932737*jacob_vel_inv1[2]*nuSum[2]*v1[2]-2.738612787525831*jacob_vel_inv1[1]*nuSum[2]*v1[2]+1.58113883008419*jacob_vel_inv1[0]*nuSum[2]*v1[2]-2.23606797749979*jacob_vel_inv1[2]*sumNuUy[2]+1.732050807568877*jacob_vel_inv1[1]*sumNuUy[2]-1.0*jacob_vel_inv1[0]*sumNuUy[2]-2.738612787525831*v1[1]*jacob_vel_inv1[2]*nuSum[2]+1.58113883008419*v1[0]*jacob_vel_inv1[2]*nuSum[2]+2.121320343559642*jacob_vel_inv1[1]*v1[1]*nuSum[2]-1.224744871391589*jacob_vel_inv1[0]*v1[1]*nuSum[2]-1.224744871391589*v1[0]*jacob_vel_inv1[1]*nuSum[2]+0.7071067811865475*jacob_vel_inv1[0]*v1[0]*nuSum[2]; 
  alphaDrSurf_l[4] = (-4.183300132670378*jacob_vel_inv1[2]*nuSum[3]*v1[3])+3.24037034920393*jacob_vel_inv1[1]*nuSum[3]*v1[3]-1.870828693386971*jacob_vel_inv1[0]*nuSum[3]*v1[3]-2.23606797749979*jacob_vel_inv1[2]*sumNuUy[3]+1.732050807568877*jacob_vel_inv1[1]*sumNuUy[3]-1.0*jacob_vel_inv1[0]*sumNuUy[3]+3.535533905932737*jacob_vel_inv1[2]*v1[2]*nuSum[3]-2.738612787525831*jacob_vel_inv1[1]*v1[2]*nuSum[3]+1.58113883008419*jacob_vel_inv1[0]*v1[2]*nuSum[3]-2.738612787525831*v1[1]*jacob_vel_inv1[2]*nuSum[3]+1.58113883008419*v1[0]*jacob_vel_inv1[2]*nuSum[3]+2.121320343559642*jacob_vel_inv1[1]*v1[1]*nuSum[3]-1.224744871391589*jacob_vel_inv1[0]*v1[1]*nuSum[3]-1.224744871391589*v1[0]*jacob_vel_inv1[1]*nuSum[3]+0.7071067811865475*jacob_vel_inv1[0]*v1[0]*nuSum[3]; 

  double alphaDrSurf_r[12] = {0.0}; 
  alphaDrSurf_r[0] = 4.183300132670378*nuSum[0]*jacob_vel_inv1[2]*v1[3]+3.24037034920393*nuSum[0]*jacob_vel_inv1[1]*v1[3]+1.870828693386971*jacob_vel_inv1[0]*nuSum[0]*v1[3]+3.535533905932737*nuSum[0]*jacob_vel_inv1[2]*v1[2]+2.738612787525831*nuSum[0]*jacob_vel_inv1[1]*v1[2]+1.58113883008419*jacob_vel_inv1[0]*nuSum[0]*v1[2]+2.738612787525831*nuSum[0]*v1[1]*jacob_vel_inv1[2]+1.58113883008419*nuSum[0]*v1[0]*jacob_vel_inv1[2]-2.23606797749979*sumNuUy[0]*jacob_vel_inv1[2]+2.121320343559642*nuSum[0]*jacob_vel_inv1[1]*v1[1]+1.224744871391589*jacob_vel_inv1[0]*nuSum[0]*v1[1]+1.224744871391589*nuSum[0]*v1[0]*jacob_vel_inv1[1]-1.732050807568877*sumNuUy[0]*jacob_vel_inv1[1]+0.7071067811865475*jacob_vel_inv1[0]*nuSum[0]*v1[0]-1.0*jacob_vel_inv1[0]*sumNuUy[0]; 
  alphaDrSurf_r[1] = 4.183300132670378*nuSum[1]*jacob_vel_inv1[2]*v1[3]+3.24037034920393*jacob_vel_inv1[1]*nuSum[1]*v1[3]+1.870828693386971*jacob_vel_inv1[0]*nuSum[1]*v1[3]+3.535533905932737*nuSum[1]*jacob_vel_inv1[2]*v1[2]+2.738612787525831*jacob_vel_inv1[1]*nuSum[1]*v1[2]+1.58113883008419*jacob_vel_inv1[0]*nuSum[1]*v1[2]+2.738612787525831*nuSum[1]*v1[1]*jacob_vel_inv1[2]-2.23606797749979*sumNuUy[1]*jacob_vel_inv1[2]+1.58113883008419*v1[0]*nuSum[1]*jacob_vel_inv1[2]+2.121320343559642*jacob_vel_inv1[1]*nuSum[1]*v1[1]+1.224744871391589*jacob_vel_inv1[0]*nuSum[1]*v1[1]-1.732050807568877*jacob_vel_inv1[1]*sumNuUy[1]-1.0*jacob_vel_inv1[0]*sumNuUy[1]+1.224744871391589*v1[0]*jacob_vel_inv1[1]*nuSum[1]+0.7071067811865475*jacob_vel_inv1[0]*v1[0]*nuSum[1]; 
  alphaDrSurf_r[2] = 4.183300132670378*jacob_vel_inv1[2]*nuSum[2]*v1[3]+3.24037034920393*jacob_vel_inv1[1]*nuSum[2]*v1[3]+1.870828693386971*jacob_vel_inv1[0]*nuSum[2]*v1[3]+3.535533905932737*jacob_vel_inv1[2]*nuSum[2]*v1[2]+2.738612787525831*jacob_vel_inv1[1]*nuSum[2]*v1[2]+1.58113883008419*jacob_vel_inv1[0]*nuSum[2]*v1[2]-2.23606797749979*jacob_vel_inv1[2]*sumNuUy[2]-1.732050807568877*jacob_vel_inv1[1]*sumNuUy[2]-1.0*jacob_vel_inv1[0]*sumNuUy[2]+2.738612787525831*v1[1]*jacob_vel_inv1[2]*nuSum[2]+1.58113883008419*v1[0]*jacob_vel_inv1[2]*nuSum[2]+2.121320343559642*jacob_vel_inv1[1]*v1[1]*nuSum[2]+1.224744871391589*jacob_vel_inv1[0]*v1[1]*nuSum[2]+1.224744871391589*v1[0]*jacob_vel_inv1[1]*nuSum[2]+0.7071067811865475*jacob_vel_inv1[0]*v1[0]*nuSum[2]; 
  alphaDrSurf_r[4] = 4.183300132670378*jacob_vel_inv1[2]*nuSum[3]*v1[3]+3.24037034920393*jacob_vel_inv1[1]*nuSum[3]*v1[3]+1.870828693386971*jacob_vel_inv1[0]*nuSum[3]*v1[3]-2.23606797749979*jacob_vel_inv1[2]*sumNuUy[3]-1.732050807568877*jacob_vel_inv1[1]*sumNuUy[3]-1.0*jacob_vel_inv1[0]*sumNuUy[3]+3.535533905932737*jacob_vel_inv1[2]*v1[2]*nuSum[3]+2.738612787525831*jacob_vel_inv1[1]*v1[2]*nuSum[3]+1.58113883008419*jacob_vel_inv1[0]*v1[2]*nuSum[3]+2.738612787525831*v1[1]*jacob_vel_inv1[2]*nuSum[3]+1.58113883008419*v1[0]*jacob_vel_inv1[2]*nuSum[3]+2.121320343559642*jacob_vel_inv1[1]*v1[1]*nuSum[3]+1.224744871391589*jacob_vel_inv1[0]*v1[1]*nuSum[3]+1.224744871391589*v1[0]*jacob_vel_inv1[1]*nuSum[3]+0.7071067811865475*jacob_vel_inv1[0]*v1[0]*nuSum[3]; 

  double fUpwindQuad_l[12] = {0.0};
  double fUpwindQuad_r[12] = {0.0};
  double fUpwind_l[12] = {0.0};
  double fUpwind_r[12] = {0.0};
  double Ghat_l[12] = {0.0}; 
  double Ghat_r[12] = {0.0}; 

  if (alphaDrSurf_l[4]-alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[0] = hyb_2x2v_p1_surfx4_eval_quad_node_0_r(fl); 
    fUpwindQuad_l[1] = hyb_2x2v_p1_surfx4_eval_quad_node_1_r(fl); 
    fUpwindQuad_l[2] = hyb_2x2v_p1_surfx4_eval_quad_node_2_r(fl); 
  } else { 
    fUpwindQuad_l[0] = hyb_2x2v_p1_surfx4_eval_quad_node_0_l(fc); 
    fUpwindQuad_l[1] = hyb_2x2v_p1_surfx4_eval_quad_node_1_l(fc); 
    fUpwindQuad_l[2] = hyb_2x2v_p1_surfx4_eval_quad_node_2_l(fc); 
  } 
  if (alphaDrSurf_r[4]-alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[0] = hyb_2x2v_p1_surfx4_eval_quad_node_0_r(fc); 
    fUpwindQuad_r[1] = hyb_2x2v_p1_surfx4_eval_quad_node_1_r(fc); 
    fUpwindQuad_r[2] = hyb_2x2v_p1_surfx4_eval_quad_node_2_r(fc); 
  } else { 
    fUpwindQuad_r[0] = hyb_2x2v_p1_surfx4_eval_quad_node_0_l(fr); 
    fUpwindQuad_r[1] = hyb_2x2v_p1_surfx4_eval_quad_node_1_l(fr); 
    fUpwindQuad_r[2] = hyb_2x2v_p1_surfx4_eval_quad_node_2_l(fr); 
  } 
  if (alphaDrSurf_l[4]-alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[3] = hyb_2x2v_p1_surfx4_eval_quad_node_3_r(fl); 
    fUpwindQuad_l[4] = hyb_2x2v_p1_surfx4_eval_quad_node_4_r(fl); 
    fUpwindQuad_l[5] = hyb_2x2v_p1_surfx4_eval_quad_node_5_r(fl); 
  } else { 
    fUpwindQuad_l[3] = hyb_2x2v_p1_surfx4_eval_quad_node_3_l(fc); 
    fUpwindQuad_l[4] = hyb_2x2v_p1_surfx4_eval_quad_node_4_l(fc); 
    fUpwindQuad_l[5] = hyb_2x2v_p1_surfx4_eval_quad_node_5_l(fc); 
  } 
  if (alphaDrSurf_r[4]-alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[3] = hyb_2x2v_p1_surfx4_eval_quad_node_3_r(fc); 
    fUpwindQuad_r[4] = hyb_2x2v_p1_surfx4_eval_quad_node_4_r(fc); 
    fUpwindQuad_r[5] = hyb_2x2v_p1_surfx4_eval_quad_node_5_r(fc); 
  } else { 
    fUpwindQuad_r[3] = hyb_2x2v_p1_surfx4_eval_quad_node_3_l(fr); 
    fUpwindQuad_r[4] = hyb_2x2v_p1_surfx4_eval_quad_node_4_l(fr); 
    fUpwindQuad_r[5] = hyb_2x2v_p1_surfx4_eval_quad_node_5_l(fr); 
  } 
  if (alphaDrSurf_l[4]-alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[6] = hyb_2x2v_p1_surfx4_eval_quad_node_6_r(fl); 
    fUpwindQuad_l[7] = hyb_2x2v_p1_surfx4_eval_quad_node_7_r(fl); 
    fUpwindQuad_l[8] = hyb_2x2v_p1_surfx4_eval_quad_node_8_r(fl); 
  } else { 
    fUpwindQuad_l[6] = hyb_2x2v_p1_surfx4_eval_quad_node_6_l(fc); 
    fUpwindQuad_l[7] = hyb_2x2v_p1_surfx4_eval_quad_node_7_l(fc); 
    fUpwindQuad_l[8] = hyb_2x2v_p1_surfx4_eval_quad_node_8_l(fc); 
  } 
  if (alphaDrSurf_r[4]-alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[6] = hyb_2x2v_p1_surfx4_eval_quad_node_6_r(fc); 
    fUpwindQuad_r[7] = hyb_2x2v_p1_surfx4_eval_quad_node_7_r(fc); 
    fUpwindQuad_r[8] = hyb_2x2v_p1_surfx4_eval_quad_node_8_r(fc); 
  } else { 
    fUpwindQuad_r[6] = hyb_2x2v_p1_surfx4_eval_quad_node_6_l(fr); 
    fUpwindQuad_r[7] = hyb_2x2v_p1_surfx4_eval_quad_node_7_l(fr); 
    fUpwindQuad_r[8] = hyb_2x2v_p1_surfx4_eval_quad_node_8_l(fr); 
  } 
  if ((-alphaDrSurf_l[4])+alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[9] = hyb_2x2v_p1_surfx4_eval_quad_node_9_r(fl); 
    fUpwindQuad_l[10] = hyb_2x2v_p1_surfx4_eval_quad_node_10_r(fl); 
    fUpwindQuad_l[11] = hyb_2x2v_p1_surfx4_eval_quad_node_11_r(fl); 
  } else { 
    fUpwindQuad_l[9] = hyb_2x2v_p1_surfx4_eval_quad_node_9_l(fc); 
    fUpwindQuad_l[10] = hyb_2x2v_p1_surfx4_eval_quad_node_10_l(fc); 
    fUpwindQuad_l[11] = hyb_2x2v_p1_surfx4_eval_quad_node_11_l(fc); 
  } 
  if ((-alphaDrSurf_r[4])+alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[9] = hyb_2x2v_p1_surfx4_eval_quad_node_9_r(fc); 
    fUpwindQuad_r[10] = hyb_2x2v_p1_surfx4_eval_quad_node_10_r(fc); 
    fUpwindQuad_r[11] = hyb_2x2v_p1_surfx4_eval_quad_node_11_r(fc); 
  } else { 
    fUpwindQuad_r[9] = hyb_2x2v_p1_surfx4_eval_quad_node_9_l(fr); 
    fUpwindQuad_r[10] = hyb_2x2v_p1_surfx4_eval_quad_node_10_l(fr); 
    fUpwindQuad_r[11] = hyb_2x2v_p1_surfx4_eval_quad_node_11_l(fr); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_2x2v_p1_vdir_upwind_quad_to_modal(fUpwindQuad_l, fUpwind_l); 
  hyb_2x2v_p1_vdir_upwind_quad_to_modal(fUpwindQuad_r, fUpwind_r); 

  Ghat_l[0] = 0.3535533905932737*(alphaDrSurf_l[4]*fUpwind_l[4]+alphaDrSurf_l[2]*fUpwind_l[2]+alphaDrSurf_l[1]*fUpwind_l[1]+alphaDrSurf_l[0]*fUpwind_l[0]); 
  Ghat_l[1] = 0.3535533905932737*(alphaDrSurf_l[2]*fUpwind_l[4]+fUpwind_l[2]*alphaDrSurf_l[4]+alphaDrSurf_l[0]*fUpwind_l[1]+fUpwind_l[0]*alphaDrSurf_l[1]); 
  Ghat_l[2] = 0.3535533905932737*(alphaDrSurf_l[1]*fUpwind_l[4]+fUpwind_l[1]*alphaDrSurf_l[4]+alphaDrSurf_l[0]*fUpwind_l[2]+fUpwind_l[0]*alphaDrSurf_l[2]); 
  Ghat_l[3] = 0.3535533905932737*(alphaDrSurf_l[4]*fUpwind_l[7]+alphaDrSurf_l[2]*fUpwind_l[6]+alphaDrSurf_l[1]*fUpwind_l[5]+alphaDrSurf_l[0]*fUpwind_l[3]); 
  Ghat_l[4] = 0.3535533905932737*(alphaDrSurf_l[0]*fUpwind_l[4]+fUpwind_l[0]*alphaDrSurf_l[4]+alphaDrSurf_l[1]*fUpwind_l[2]+fUpwind_l[1]*alphaDrSurf_l[2]); 
  Ghat_l[5] = 0.3535533905932737*(alphaDrSurf_l[2]*fUpwind_l[7]+alphaDrSurf_l[4]*fUpwind_l[6]+alphaDrSurf_l[0]*fUpwind_l[5]+alphaDrSurf_l[1]*fUpwind_l[3]); 
  Ghat_l[6] = 0.3535533905932737*(alphaDrSurf_l[1]*fUpwind_l[7]+alphaDrSurf_l[0]*fUpwind_l[6]+alphaDrSurf_l[4]*fUpwind_l[5]+alphaDrSurf_l[2]*fUpwind_l[3]); 
  Ghat_l[7] = 0.3535533905932737*(alphaDrSurf_l[0]*fUpwind_l[7]+alphaDrSurf_l[1]*fUpwind_l[6]+alphaDrSurf_l[2]*fUpwind_l[5]+fUpwind_l[3]*alphaDrSurf_l[4]); 
  Ghat_l[8] = 0.02357022603955158*(15.0*alphaDrSurf_l[4]*fUpwind_l[11]+15.0*(alphaDrSurf_l[2]*fUpwind_l[10]+alphaDrSurf_l[1]*fUpwind_l[9])+15.0*alphaDrSurf_l[0]*fUpwind_l[8]); 
  Ghat_l[9] = 0.02357022603955158*(15.0*alphaDrSurf_l[2]*fUpwind_l[11]+15.0*(alphaDrSurf_l[4]*fUpwind_l[10]+alphaDrSurf_l[0]*fUpwind_l[9])+15.0*alphaDrSurf_l[1]*fUpwind_l[8]); 
  Ghat_l[10] = 0.02357022603955158*(15.0*alphaDrSurf_l[1]*fUpwind_l[11]+15.0*(alphaDrSurf_l[0]*fUpwind_l[10]+alphaDrSurf_l[4]*fUpwind_l[9])+15.0*alphaDrSurf_l[2]*fUpwind_l[8]); 
  Ghat_l[11] = 0.02357022603955158*(15.0*alphaDrSurf_l[0]*fUpwind_l[11]+15.0*(alphaDrSurf_l[1]*fUpwind_l[10]+alphaDrSurf_l[2]*fUpwind_l[9])+15.0*alphaDrSurf_l[4]*fUpwind_l[8]); 

  Ghat_r[0] = 0.3535533905932737*(alphaDrSurf_r[4]*fUpwind_r[4]+alphaDrSurf_r[2]*fUpwind_r[2]+alphaDrSurf_r[1]*fUpwind_r[1]+alphaDrSurf_r[0]*fUpwind_r[0]); 
  Ghat_r[1] = 0.3535533905932737*(alphaDrSurf_r[2]*fUpwind_r[4]+fUpwind_r[2]*alphaDrSurf_r[4]+alphaDrSurf_r[0]*fUpwind_r[1]+fUpwind_r[0]*alphaDrSurf_r[1]); 
  Ghat_r[2] = 0.3535533905932737*(alphaDrSurf_r[1]*fUpwind_r[4]+fUpwind_r[1]*alphaDrSurf_r[4]+alphaDrSurf_r[0]*fUpwind_r[2]+fUpwind_r[0]*alphaDrSurf_r[2]); 
  Ghat_r[3] = 0.3535533905932737*(alphaDrSurf_r[4]*fUpwind_r[7]+alphaDrSurf_r[2]*fUpwind_r[6]+alphaDrSurf_r[1]*fUpwind_r[5]+alphaDrSurf_r[0]*fUpwind_r[3]); 
  Ghat_r[4] = 0.3535533905932737*(alphaDrSurf_r[0]*fUpwind_r[4]+fUpwind_r[0]*alphaDrSurf_r[4]+alphaDrSurf_r[1]*fUpwind_r[2]+fUpwind_r[1]*alphaDrSurf_r[2]); 
  Ghat_r[5] = 0.3535533905932737*(alphaDrSurf_r[2]*fUpwind_r[7]+alphaDrSurf_r[4]*fUpwind_r[6]+alphaDrSurf_r[0]*fUpwind_r[5]+alphaDrSurf_r[1]*fUpwind_r[3]); 
  Ghat_r[6] = 0.3535533905932737*(alphaDrSurf_r[1]*fUpwind_r[7]+alphaDrSurf_r[0]*fUpwind_r[6]+alphaDrSurf_r[4]*fUpwind_r[5]+alphaDrSurf_r[2]*fUpwind_r[3]); 
  Ghat_r[7] = 0.3535533905932737*(alphaDrSurf_r[0]*fUpwind_r[7]+alphaDrSurf_r[1]*fUpwind_r[6]+alphaDrSurf_r[2]*fUpwind_r[5]+fUpwind_r[3]*alphaDrSurf_r[4]); 
  Ghat_r[8] = 0.02357022603955158*(15.0*alphaDrSurf_r[4]*fUpwind_r[11]+15.0*(alphaDrSurf_r[2]*fUpwind_r[10]+alphaDrSurf_r[1]*fUpwind_r[9])+15.0*alphaDrSurf_r[0]*fUpwind_r[8]); 
  Ghat_r[9] = 0.02357022603955158*(15.0*alphaDrSurf_r[2]*fUpwind_r[11]+15.0*(alphaDrSurf_r[4]*fUpwind_r[10]+alphaDrSurf_r[0]*fUpwind_r[9])+15.0*alphaDrSurf_r[1]*fUpwind_r[8]); 
  Ghat_r[10] = 0.02357022603955158*(15.0*alphaDrSurf_r[1]*fUpwind_r[11]+15.0*(alphaDrSurf_r[0]*fUpwind_r[10]+alphaDrSurf_r[4]*fUpwind_r[9])+15.0*alphaDrSurf_r[2]*fUpwind_r[8]); 
  Ghat_r[11] = 0.02357022603955158*(15.0*alphaDrSurf_r[0]*fUpwind_r[11]+15.0*(alphaDrSurf_r[1]*fUpwind_r[10]+alphaDrSurf_r[2]*fUpwind_r[9])+15.0*alphaDrSurf_r[4]*fUpwind_r[8]); 

  out[0] += 0.7071067811865475*Ghat_r[0]*rdv2-0.7071067811865475*Ghat_l[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat_r[1]*rdv2-0.7071067811865475*Ghat_l[1]*rdv2; 
  out[2] += 0.7071067811865475*Ghat_r[2]*rdv2-0.7071067811865475*Ghat_l[2]*rdv2; 
  out[3] += 0.7071067811865475*Ghat_r[3]*rdv2-0.7071067811865475*Ghat_l[3]*rdv2; 
  out[4] += 1.224744871391589*Ghat_r[0]*rdv2+1.224744871391589*Ghat_l[0]*rdv2; 
  out[5] += 0.7071067811865475*Ghat_r[4]*rdv2-0.7071067811865475*Ghat_l[4]*rdv2; 
  out[6] += 0.7071067811865475*Ghat_r[5]*rdv2-0.7071067811865475*Ghat_l[5]*rdv2; 
  out[7] += 0.7071067811865475*Ghat_r[6]*rdv2-0.7071067811865475*Ghat_l[6]*rdv2; 
  out[8] += 1.224744871391589*Ghat_r[1]*rdv2+1.224744871391589*Ghat_l[1]*rdv2; 
  out[9] += 1.224744871391589*Ghat_r[2]*rdv2+1.224744871391589*Ghat_l[2]*rdv2; 
  out[10] += 1.224744871391589*Ghat_r[3]*rdv2+1.224744871391589*Ghat_l[3]*rdv2; 
  out[11] += 0.7071067811865475*Ghat_r[7]*rdv2-0.7071067811865475*Ghat_l[7]*rdv2; 
  out[12] += 1.224744871391589*Ghat_r[4]*rdv2+1.224744871391589*Ghat_l[4]*rdv2; 
  out[13] += 1.224744871391589*Ghat_r[5]*rdv2+1.224744871391589*Ghat_l[5]*rdv2; 
  out[14] += 1.224744871391589*Ghat_r[6]*rdv2+1.224744871391589*Ghat_l[6]*rdv2; 
  out[15] += 1.224744871391589*Ghat_r[7]*rdv2+1.224744871391589*Ghat_l[7]*rdv2; 
  out[16] += 0.7071067811865475*Ghat_r[8]*rdv2-0.7071067811865475*Ghat_l[8]*rdv2; 
  out[17] += 0.7071067811865475*Ghat_r[9]*rdv2-0.7071067811865475*Ghat_l[9]*rdv2; 
  out[18] += 0.7071067811865475*Ghat_r[10]*rdv2-0.7071067811865475*Ghat_l[10]*rdv2; 
  out[19] += 1.224744871391589*Ghat_r[8]*rdv2+1.224744871391589*Ghat_l[8]*rdv2; 
  out[20] += 0.7071067811865475*Ghat_r[11]*rdv2-0.7071067811865475*Ghat_l[11]*rdv2; 
  out[21] += 1.224744871391589*Ghat_r[9]*rdv2+1.224744871391589*Ghat_l[9]*rdv2; 
  out[22] += 1.224744871391589*Ghat_r[10]*rdv2+1.224744871391589*Ghat_l[10]*rdv2; 
  out[23] += 1.224744871391589*Ghat_r[11]*rdv2+1.224744871391589*Ghat_l[11]*rdv2; 
  out[24] += 1.58113883008419*Ghat_r[0]*rdv2-1.58113883008419*Ghat_l[0]*rdv2; 
  out[25] += 1.58113883008419*Ghat_r[1]*rdv2-1.58113883008419*Ghat_l[1]*rdv2; 
  out[26] += 1.58113883008419*Ghat_r[2]*rdv2-1.58113883008419*Ghat_l[2]*rdv2; 
  out[27] += 1.58113883008419*Ghat_r[3]*rdv2-1.58113883008419*Ghat_l[3]*rdv2; 
  out[28] += 1.58113883008419*Ghat_r[4]*rdv2-1.58113883008419*Ghat_l[4]*rdv2; 
  out[29] += 1.58113883008419*Ghat_r[5]*rdv2-1.58113883008419*Ghat_l[5]*rdv2; 
  out[30] += 1.58113883008419*Ghat_r[6]*rdv2-1.58113883008419*Ghat_l[6]*rdv2; 
  out[31] += 1.58113883008419*Ghat_r[7]*rdv2-1.58113883008419*Ghat_l[7]*rdv2; 

  return 0.;

} 
