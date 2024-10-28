#include <gkyl_lbo_vlasov_kernels.h> 
#include <gkyl_basis_hyb_1x2v_p1_surfx2_eval_quad.h> 
#include <gkyl_basis_hyb_1x2v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double lbo_vlasov_vmap_drag_surfvx_1x2v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_vel_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[3]: cell-center coordinates. 
  // dxv[3]: cell spacing. 
  // vmap: Velocity-space nonuniform mapping in each dimension. 
  // jacob_vel_inv: Inverse of velocity space Jacobian in each dimension. 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[6]: sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 
  // fl/fc/fr: distribution function in cells 
  // out: incremented distribution function in cell 
  double rdv2 = 2.0/dxv[1]; 

  const double *v0 = &vmap[0]; 
  const double *jacob_vel_inv0 = &jacob_vel_inv[0]; 
  const double *sumNuUx = &nuPrimMomsSum[0]; 

  double alphaDrSurf_l[6] = {0.0}; 
  alphaDrSurf_l[0] = (-4.183300132670378*nuSum[0]*jacob_vel_inv0[2]*v0[3])+3.24037034920393*nuSum[0]*jacob_vel_inv0[1]*v0[3]-1.870828693386971*jacob_vel_inv0[0]*nuSum[0]*v0[3]+3.535533905932737*nuSum[0]*jacob_vel_inv0[2]*v0[2]-2.738612787525831*nuSum[0]*jacob_vel_inv0[1]*v0[2]+1.58113883008419*jacob_vel_inv0[0]*nuSum[0]*v0[2]-2.738612787525831*nuSum[0]*v0[1]*jacob_vel_inv0[2]+1.58113883008419*nuSum[0]*v0[0]*jacob_vel_inv0[2]-2.23606797749979*sumNuUx[0]*jacob_vel_inv0[2]+2.121320343559642*nuSum[0]*jacob_vel_inv0[1]*v0[1]-1.224744871391589*jacob_vel_inv0[0]*nuSum[0]*v0[1]-1.224744871391589*nuSum[0]*v0[0]*jacob_vel_inv0[1]+1.732050807568877*sumNuUx[0]*jacob_vel_inv0[1]+0.7071067811865475*jacob_vel_inv0[0]*nuSum[0]*v0[0]-1.0*jacob_vel_inv0[0]*sumNuUx[0]; 
  alphaDrSurf_l[1] = (-4.183300132670378*nuSum[1]*jacob_vel_inv0[2]*v0[3])+3.24037034920393*jacob_vel_inv0[1]*nuSum[1]*v0[3]-1.870828693386971*jacob_vel_inv0[0]*nuSum[1]*v0[3]+3.535533905932737*nuSum[1]*jacob_vel_inv0[2]*v0[2]-2.738612787525831*jacob_vel_inv0[1]*nuSum[1]*v0[2]+1.58113883008419*jacob_vel_inv0[0]*nuSum[1]*v0[2]-2.738612787525831*nuSum[1]*v0[1]*jacob_vel_inv0[2]-2.23606797749979*sumNuUx[1]*jacob_vel_inv0[2]+1.58113883008419*v0[0]*nuSum[1]*jacob_vel_inv0[2]+2.121320343559642*jacob_vel_inv0[1]*nuSum[1]*v0[1]-1.224744871391589*jacob_vel_inv0[0]*nuSum[1]*v0[1]+1.732050807568877*jacob_vel_inv0[1]*sumNuUx[1]-1.0*jacob_vel_inv0[0]*sumNuUx[1]-1.224744871391589*v0[0]*jacob_vel_inv0[1]*nuSum[1]+0.7071067811865475*jacob_vel_inv0[0]*v0[0]*nuSum[1]; 

  double alphaDrSurf_r[6] = {0.0}; 
  alphaDrSurf_r[0] = 4.183300132670378*nuSum[0]*jacob_vel_inv0[2]*v0[3]+3.24037034920393*nuSum[0]*jacob_vel_inv0[1]*v0[3]+1.870828693386971*jacob_vel_inv0[0]*nuSum[0]*v0[3]+3.535533905932737*nuSum[0]*jacob_vel_inv0[2]*v0[2]+2.738612787525831*nuSum[0]*jacob_vel_inv0[1]*v0[2]+1.58113883008419*jacob_vel_inv0[0]*nuSum[0]*v0[2]+2.738612787525831*nuSum[0]*v0[1]*jacob_vel_inv0[2]+1.58113883008419*nuSum[0]*v0[0]*jacob_vel_inv0[2]-2.23606797749979*sumNuUx[0]*jacob_vel_inv0[2]+2.121320343559642*nuSum[0]*jacob_vel_inv0[1]*v0[1]+1.224744871391589*jacob_vel_inv0[0]*nuSum[0]*v0[1]+1.224744871391589*nuSum[0]*v0[0]*jacob_vel_inv0[1]-1.732050807568877*sumNuUx[0]*jacob_vel_inv0[1]+0.7071067811865475*jacob_vel_inv0[0]*nuSum[0]*v0[0]-1.0*jacob_vel_inv0[0]*sumNuUx[0]; 
  alphaDrSurf_r[1] = 4.183300132670378*nuSum[1]*jacob_vel_inv0[2]*v0[3]+3.24037034920393*jacob_vel_inv0[1]*nuSum[1]*v0[3]+1.870828693386971*jacob_vel_inv0[0]*nuSum[1]*v0[3]+3.535533905932737*nuSum[1]*jacob_vel_inv0[2]*v0[2]+2.738612787525831*jacob_vel_inv0[1]*nuSum[1]*v0[2]+1.58113883008419*jacob_vel_inv0[0]*nuSum[1]*v0[2]+2.738612787525831*nuSum[1]*v0[1]*jacob_vel_inv0[2]-2.23606797749979*sumNuUx[1]*jacob_vel_inv0[2]+1.58113883008419*v0[0]*nuSum[1]*jacob_vel_inv0[2]+2.121320343559642*jacob_vel_inv0[1]*nuSum[1]*v0[1]+1.224744871391589*jacob_vel_inv0[0]*nuSum[1]*v0[1]-1.732050807568877*jacob_vel_inv0[1]*sumNuUx[1]-1.0*jacob_vel_inv0[0]*sumNuUx[1]+1.224744871391589*v0[0]*jacob_vel_inv0[1]*nuSum[1]+0.7071067811865475*jacob_vel_inv0[0]*v0[0]*nuSum[1]; 

  double fUpwindQuad_l[6] = {0.0};
  double fUpwindQuad_r[6] = {0.0};
  double fUpwind_l[6] = {0.0};
  double fUpwind_r[6] = {0.0};
  double Ghat_l[6] = {0.0}; 
  double Ghat_r[6] = {0.0}; 

  if (alphaDrSurf_l[0]-alphaDrSurf_l[1] < 0) { 
    fUpwindQuad_l[0] = hyb_1x2v_p1_surfx2_eval_quad_node_0_r(fl); 
    fUpwindQuad_l[1] = hyb_1x2v_p1_surfx2_eval_quad_node_1_r(fl); 
    fUpwindQuad_l[2] = hyb_1x2v_p1_surfx2_eval_quad_node_2_r(fl); 
  } else { 
    fUpwindQuad_l[0] = hyb_1x2v_p1_surfx2_eval_quad_node_0_l(fc); 
    fUpwindQuad_l[1] = hyb_1x2v_p1_surfx2_eval_quad_node_1_l(fc); 
    fUpwindQuad_l[2] = hyb_1x2v_p1_surfx2_eval_quad_node_2_l(fc); 
  } 
  if (alphaDrSurf_r[0]-alphaDrSurf_r[1] < 0) { 
    fUpwindQuad_r[0] = hyb_1x2v_p1_surfx2_eval_quad_node_0_r(fc); 
    fUpwindQuad_r[1] = hyb_1x2v_p1_surfx2_eval_quad_node_1_r(fc); 
    fUpwindQuad_r[2] = hyb_1x2v_p1_surfx2_eval_quad_node_2_r(fc); 
  } else { 
    fUpwindQuad_r[0] = hyb_1x2v_p1_surfx2_eval_quad_node_0_l(fr); 
    fUpwindQuad_r[1] = hyb_1x2v_p1_surfx2_eval_quad_node_1_l(fr); 
    fUpwindQuad_r[2] = hyb_1x2v_p1_surfx2_eval_quad_node_2_l(fr); 
  } 
  if (alphaDrSurf_l[0]-alphaDrSurf_l[1] < 0) { 
    fUpwindQuad_l[3] = hyb_1x2v_p1_surfx2_eval_quad_node_3_r(fl); 
    fUpwindQuad_l[4] = hyb_1x2v_p1_surfx2_eval_quad_node_4_r(fl); 
    fUpwindQuad_l[5] = hyb_1x2v_p1_surfx2_eval_quad_node_5_r(fl); 
  } else { 
    fUpwindQuad_l[3] = hyb_1x2v_p1_surfx2_eval_quad_node_3_l(fc); 
    fUpwindQuad_l[4] = hyb_1x2v_p1_surfx2_eval_quad_node_4_l(fc); 
    fUpwindQuad_l[5] = hyb_1x2v_p1_surfx2_eval_quad_node_5_l(fc); 
  } 
  if (alphaDrSurf_r[0]-alphaDrSurf_r[1] < 0) { 
    fUpwindQuad_r[3] = hyb_1x2v_p1_surfx2_eval_quad_node_3_r(fc); 
    fUpwindQuad_r[4] = hyb_1x2v_p1_surfx2_eval_quad_node_4_r(fc); 
    fUpwindQuad_r[5] = hyb_1x2v_p1_surfx2_eval_quad_node_5_r(fc); 
  } else { 
    fUpwindQuad_r[3] = hyb_1x2v_p1_surfx2_eval_quad_node_3_l(fr); 
    fUpwindQuad_r[4] = hyb_1x2v_p1_surfx2_eval_quad_node_4_l(fr); 
    fUpwindQuad_r[5] = hyb_1x2v_p1_surfx2_eval_quad_node_5_l(fr); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_1x2v_p1_vdir_upwind_quad_to_modal(fUpwindQuad_l, fUpwind_l); 
  hyb_1x2v_p1_vdir_upwind_quad_to_modal(fUpwindQuad_r, fUpwind_r); 

  Ghat_l[0] = 0.5*(alphaDrSurf_l[1]*fUpwind_l[1]+alphaDrSurf_l[0]*fUpwind_l[0]); 
  Ghat_l[1] = 0.5*(alphaDrSurf_l[0]*fUpwind_l[1]+fUpwind_l[0]*alphaDrSurf_l[1]); 
  Ghat_l[2] = 0.5*(alphaDrSurf_l[1]*fUpwind_l[3]+alphaDrSurf_l[0]*fUpwind_l[2]); 
  Ghat_l[3] = 0.5*(alphaDrSurf_l[0]*fUpwind_l[3]+alphaDrSurf_l[1]*fUpwind_l[2]); 
  Ghat_l[4] = 0.03333333333333333*(15.0*alphaDrSurf_l[1]*fUpwind_l[5]+15.0*alphaDrSurf_l[0]*fUpwind_l[4]); 
  Ghat_l[5] = 0.03333333333333333*(15.0*alphaDrSurf_l[0]*fUpwind_l[5]+15.0*alphaDrSurf_l[1]*fUpwind_l[4]); 

  Ghat_r[0] = 0.5*(alphaDrSurf_r[1]*fUpwind_r[1]+alphaDrSurf_r[0]*fUpwind_r[0]); 
  Ghat_r[1] = 0.5*(alphaDrSurf_r[0]*fUpwind_r[1]+fUpwind_r[0]*alphaDrSurf_r[1]); 
  Ghat_r[2] = 0.5*(alphaDrSurf_r[1]*fUpwind_r[3]+alphaDrSurf_r[0]*fUpwind_r[2]); 
  Ghat_r[3] = 0.5*(alphaDrSurf_r[0]*fUpwind_r[3]+alphaDrSurf_r[1]*fUpwind_r[2]); 
  Ghat_r[4] = 0.03333333333333333*(15.0*alphaDrSurf_r[1]*fUpwind_r[5]+15.0*alphaDrSurf_r[0]*fUpwind_r[4]); 
  Ghat_r[5] = 0.03333333333333333*(15.0*alphaDrSurf_r[0]*fUpwind_r[5]+15.0*alphaDrSurf_r[1]*fUpwind_r[4]); 

  out[0] += 0.7071067811865475*Ghat_r[0]*rdv2-0.7071067811865475*Ghat_l[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat_r[1]*rdv2-0.7071067811865475*Ghat_l[1]*rdv2; 
  out[2] += 1.224744871391589*Ghat_r[0]*rdv2+1.224744871391589*Ghat_l[0]*rdv2; 
  out[3] += 0.7071067811865475*Ghat_r[2]*rdv2-0.7071067811865475*Ghat_l[2]*rdv2; 
  out[4] += 1.224744871391589*Ghat_r[1]*rdv2+1.224744871391589*Ghat_l[1]*rdv2; 
  out[5] += 0.7071067811865475*Ghat_r[3]*rdv2-0.7071067811865475*Ghat_l[3]*rdv2; 
  out[6] += 1.224744871391589*Ghat_r[2]*rdv2+1.224744871391589*Ghat_l[2]*rdv2; 
  out[7] += 1.224744871391589*Ghat_r[3]*rdv2+1.224744871391589*Ghat_l[3]*rdv2; 
  out[8] += 1.58113883008419*Ghat_r[0]*rdv2-1.58113883008419*Ghat_l[0]*rdv2; 
  out[9] += 1.58113883008419*Ghat_r[1]*rdv2-1.58113883008419*Ghat_l[1]*rdv2; 
  out[10] += 1.58113883008419*Ghat_r[2]*rdv2-1.58113883008419*Ghat_l[2]*rdv2; 
  out[11] += 1.58113883008419*Ghat_r[3]*rdv2-1.58113883008419*Ghat_l[3]*rdv2; 
  out[12] += 0.7071067811865475*Ghat_r[4]*rdv2-0.7071067811865475*Ghat_l[4]*rdv2; 
  out[13] += 0.7071067811865475*Ghat_r[5]*rdv2-0.7071067811865475*Ghat_l[5]*rdv2; 
  out[14] += 1.224744871391589*Ghat_r[4]*rdv2+1.224744871391589*Ghat_l[4]*rdv2; 
  out[15] += 1.224744871391589*Ghat_r[5]*rdv2+1.224744871391589*Ghat_l[5]*rdv2; 

  return 0.;

} 
