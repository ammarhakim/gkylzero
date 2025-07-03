#include <gkyl_lbo_vlasov_kernels.h> 
#include <gkyl_basis_ser_3x_p2_surfx3_eval_quad.h> 
#include <gkyl_basis_ser_3x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH double lbo_vlasov_vmap_drag_surfvx_2x1v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_vel_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[3]: cell-center coordinates. 
  // dxv[3]: cell spacing. 
  // vmap: Velocity-space nonuniform mapping in each dimension. 
  // jacob_vel_inv: Inverse of velocity space Jacobian in each dimension. 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[16]: sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 
  // fl/fc/fr: distribution function in cells 
  // out: incremented distribution function in cell 
  double rdv2 = 2.0/dxv[2]; 

  const double *v0 = &vmap[0]; 
  const double *jacob_vel_inv0 = &jacob_vel_inv[0]; 
  const double *sumNuUx = &nuPrimMomsSum[0]; 

  double alphaDrSurf_l[8] = {0.0}; 
  alphaDrSurf_l[0] = (-2.958039891549809*nuSum[0]*jacob_vel_inv0[2]*v0[3])+2.29128784747792*nuSum[0]*jacob_vel_inv0[1]*v0[3]-1.322875655532295*jacob_vel_inv0[0]*nuSum[0]*v0[3]+2.5*nuSum[0]*jacob_vel_inv0[2]*v0[2]-1.936491673103709*nuSum[0]*jacob_vel_inv0[1]*v0[2]+1.118033988749895*jacob_vel_inv0[0]*nuSum[0]*v0[2]-1.936491673103709*nuSum[0]*v0[1]*jacob_vel_inv0[2]+1.118033988749895*nuSum[0]*v0[0]*jacob_vel_inv0[2]-1.58113883008419*sumNuUx[0]*jacob_vel_inv0[2]+1.5*nuSum[0]*jacob_vel_inv0[1]*v0[1]-0.8660254037844386*jacob_vel_inv0[0]*nuSum[0]*v0[1]-0.8660254037844386*nuSum[0]*v0[0]*jacob_vel_inv0[1]+1.224744871391589*sumNuUx[0]*jacob_vel_inv0[1]+0.5*jacob_vel_inv0[0]*nuSum[0]*v0[0]-0.7071067811865475*jacob_vel_inv0[0]*sumNuUx[0]; 
  alphaDrSurf_l[1] = (-2.958039891549809*nuSum[1]*jacob_vel_inv0[2]*v0[3])+2.29128784747792*jacob_vel_inv0[1]*nuSum[1]*v0[3]-1.322875655532295*jacob_vel_inv0[0]*nuSum[1]*v0[3]+2.5*nuSum[1]*jacob_vel_inv0[2]*v0[2]-1.936491673103709*jacob_vel_inv0[1]*nuSum[1]*v0[2]+1.118033988749895*jacob_vel_inv0[0]*nuSum[1]*v0[2]-1.936491673103709*nuSum[1]*v0[1]*jacob_vel_inv0[2]-1.58113883008419*sumNuUx[1]*jacob_vel_inv0[2]+1.118033988749895*v0[0]*nuSum[1]*jacob_vel_inv0[2]+1.5*jacob_vel_inv0[1]*nuSum[1]*v0[1]-0.8660254037844386*jacob_vel_inv0[0]*nuSum[1]*v0[1]+1.224744871391589*jacob_vel_inv0[1]*sumNuUx[1]-0.7071067811865475*jacob_vel_inv0[0]*sumNuUx[1]-0.8660254037844386*v0[0]*jacob_vel_inv0[1]*nuSum[1]+0.5*jacob_vel_inv0[0]*v0[0]*nuSum[1]; 
  alphaDrSurf_l[2] = (-2.958039891549809*jacob_vel_inv0[2]*nuSum[2]*v0[3])+2.29128784747792*jacob_vel_inv0[1]*nuSum[2]*v0[3]-1.322875655532295*jacob_vel_inv0[0]*nuSum[2]*v0[3]+2.5*jacob_vel_inv0[2]*nuSum[2]*v0[2]-1.936491673103709*jacob_vel_inv0[1]*nuSum[2]*v0[2]+1.118033988749895*jacob_vel_inv0[0]*nuSum[2]*v0[2]-1.58113883008419*jacob_vel_inv0[2]*sumNuUx[2]+1.224744871391589*jacob_vel_inv0[1]*sumNuUx[2]-0.7071067811865475*jacob_vel_inv0[0]*sumNuUx[2]-1.936491673103709*v0[1]*jacob_vel_inv0[2]*nuSum[2]+1.118033988749895*v0[0]*jacob_vel_inv0[2]*nuSum[2]+1.5*jacob_vel_inv0[1]*v0[1]*nuSum[2]-0.8660254037844386*jacob_vel_inv0[0]*v0[1]*nuSum[2]-0.8660254037844386*v0[0]*jacob_vel_inv0[1]*nuSum[2]+0.5*jacob_vel_inv0[0]*v0[0]*nuSum[2]; 
  alphaDrSurf_l[3] = (-2.958039891549809*jacob_vel_inv0[2]*nuSum[3]*v0[3])+2.29128784747792*jacob_vel_inv0[1]*nuSum[3]*v0[3]-1.322875655532295*jacob_vel_inv0[0]*nuSum[3]*v0[3]-1.58113883008419*jacob_vel_inv0[2]*sumNuUx[3]+1.224744871391589*jacob_vel_inv0[1]*sumNuUx[3]-0.7071067811865475*jacob_vel_inv0[0]*sumNuUx[3]+2.5*jacob_vel_inv0[2]*v0[2]*nuSum[3]-1.936491673103709*jacob_vel_inv0[1]*v0[2]*nuSum[3]+1.118033988749895*jacob_vel_inv0[0]*v0[2]*nuSum[3]-1.936491673103709*v0[1]*jacob_vel_inv0[2]*nuSum[3]+1.118033988749895*v0[0]*jacob_vel_inv0[2]*nuSum[3]+1.5*jacob_vel_inv0[1]*v0[1]*nuSum[3]-0.8660254037844386*jacob_vel_inv0[0]*v0[1]*nuSum[3]-0.8660254037844386*v0[0]*jacob_vel_inv0[1]*nuSum[3]+0.5*jacob_vel_inv0[0]*v0[0]*nuSum[3]; 
  alphaDrSurf_l[4] = (-1.58113883008419*jacob_vel_inv0[2]*sumNuUx[4])+1.224744871391589*jacob_vel_inv0[1]*sumNuUx[4]-0.7071067811865475*jacob_vel_inv0[0]*sumNuUx[4]-2.958039891549809*jacob_vel_inv0[2]*v0[3]*nuSum[4]+2.29128784747792*jacob_vel_inv0[1]*v0[3]*nuSum[4]-1.322875655532295*jacob_vel_inv0[0]*v0[3]*nuSum[4]+2.5*jacob_vel_inv0[2]*v0[2]*nuSum[4]-1.936491673103709*jacob_vel_inv0[1]*v0[2]*nuSum[4]+1.118033988749895*jacob_vel_inv0[0]*v0[2]*nuSum[4]-1.936491673103709*v0[1]*jacob_vel_inv0[2]*nuSum[4]+1.118033988749895*v0[0]*jacob_vel_inv0[2]*nuSum[4]+1.5*jacob_vel_inv0[1]*v0[1]*nuSum[4]-0.8660254037844386*jacob_vel_inv0[0]*v0[1]*nuSum[4]-0.8660254037844386*v0[0]*jacob_vel_inv0[1]*nuSum[4]+0.5*jacob_vel_inv0[0]*v0[0]*nuSum[4]; 
  alphaDrSurf_l[5] = (-1.58113883008419*jacob_vel_inv0[2]*sumNuUx[5])+1.224744871391589*jacob_vel_inv0[1]*sumNuUx[5]-0.7071067811865475*jacob_vel_inv0[0]*sumNuUx[5]-2.958039891549809*jacob_vel_inv0[2]*v0[3]*nuSum[5]+2.29128784747792*jacob_vel_inv0[1]*v0[3]*nuSum[5]-1.322875655532295*jacob_vel_inv0[0]*v0[3]*nuSum[5]+2.5*jacob_vel_inv0[2]*v0[2]*nuSum[5]-1.936491673103709*jacob_vel_inv0[1]*v0[2]*nuSum[5]+1.118033988749895*jacob_vel_inv0[0]*v0[2]*nuSum[5]-1.936491673103709*v0[1]*jacob_vel_inv0[2]*nuSum[5]+1.118033988749895*v0[0]*jacob_vel_inv0[2]*nuSum[5]+1.5*jacob_vel_inv0[1]*v0[1]*nuSum[5]-0.8660254037844386*jacob_vel_inv0[0]*v0[1]*nuSum[5]-0.8660254037844386*v0[0]*jacob_vel_inv0[1]*nuSum[5]+0.5*jacob_vel_inv0[0]*v0[0]*nuSum[5]; 
  alphaDrSurf_l[6] = (-1.58113883008419*jacob_vel_inv0[2]*sumNuUx[6])+1.224744871391589*jacob_vel_inv0[1]*sumNuUx[6]-0.7071067811865475*jacob_vel_inv0[0]*sumNuUx[6]-2.958039891549809*jacob_vel_inv0[2]*v0[3]*nuSum[6]+2.29128784747792*jacob_vel_inv0[1]*v0[3]*nuSum[6]-1.322875655532295*jacob_vel_inv0[0]*v0[3]*nuSum[6]+2.5*jacob_vel_inv0[2]*v0[2]*nuSum[6]-1.936491673103709*jacob_vel_inv0[1]*v0[2]*nuSum[6]+1.118033988749895*jacob_vel_inv0[0]*v0[2]*nuSum[6]-1.936491673103709*v0[1]*jacob_vel_inv0[2]*nuSum[6]+1.118033988749895*v0[0]*jacob_vel_inv0[2]*nuSum[6]+1.5*jacob_vel_inv0[1]*v0[1]*nuSum[6]-0.8660254037844386*jacob_vel_inv0[0]*v0[1]*nuSum[6]-0.8660254037844386*v0[0]*jacob_vel_inv0[1]*nuSum[6]+0.5*jacob_vel_inv0[0]*v0[0]*nuSum[6]; 
  alphaDrSurf_l[7] = (-1.58113883008419*jacob_vel_inv0[2]*sumNuUx[7])+1.224744871391589*jacob_vel_inv0[1]*sumNuUx[7]-0.7071067811865475*jacob_vel_inv0[0]*sumNuUx[7]-2.958039891549809*jacob_vel_inv0[2]*v0[3]*nuSum[7]+2.29128784747792*jacob_vel_inv0[1]*v0[3]*nuSum[7]-1.322875655532295*jacob_vel_inv0[0]*v0[3]*nuSum[7]+2.5*jacob_vel_inv0[2]*v0[2]*nuSum[7]-1.936491673103709*jacob_vel_inv0[1]*v0[2]*nuSum[7]+1.118033988749895*jacob_vel_inv0[0]*v0[2]*nuSum[7]-1.936491673103709*v0[1]*jacob_vel_inv0[2]*nuSum[7]+1.118033988749895*v0[0]*jacob_vel_inv0[2]*nuSum[7]+1.5*jacob_vel_inv0[1]*v0[1]*nuSum[7]-0.8660254037844386*jacob_vel_inv0[0]*v0[1]*nuSum[7]-0.8660254037844386*v0[0]*jacob_vel_inv0[1]*nuSum[7]+0.5*jacob_vel_inv0[0]*v0[0]*nuSum[7]; 

  double alphaDrSurf_r[8] = {0.0}; 
  alphaDrSurf_r[0] = 2.958039891549809*nuSum[0]*jacob_vel_inv0[2]*v0[3]+2.29128784747792*nuSum[0]*jacob_vel_inv0[1]*v0[3]+1.322875655532295*jacob_vel_inv0[0]*nuSum[0]*v0[3]+2.5*nuSum[0]*jacob_vel_inv0[2]*v0[2]+1.936491673103709*nuSum[0]*jacob_vel_inv0[1]*v0[2]+1.118033988749895*jacob_vel_inv0[0]*nuSum[0]*v0[2]+1.936491673103709*nuSum[0]*v0[1]*jacob_vel_inv0[2]+1.118033988749895*nuSum[0]*v0[0]*jacob_vel_inv0[2]-1.58113883008419*sumNuUx[0]*jacob_vel_inv0[2]+1.5*nuSum[0]*jacob_vel_inv0[1]*v0[1]+0.8660254037844386*jacob_vel_inv0[0]*nuSum[0]*v0[1]+0.8660254037844386*nuSum[0]*v0[0]*jacob_vel_inv0[1]-1.224744871391589*sumNuUx[0]*jacob_vel_inv0[1]+0.5*jacob_vel_inv0[0]*nuSum[0]*v0[0]-0.7071067811865475*jacob_vel_inv0[0]*sumNuUx[0]; 
  alphaDrSurf_r[1] = 2.958039891549809*nuSum[1]*jacob_vel_inv0[2]*v0[3]+2.29128784747792*jacob_vel_inv0[1]*nuSum[1]*v0[3]+1.322875655532295*jacob_vel_inv0[0]*nuSum[1]*v0[3]+2.5*nuSum[1]*jacob_vel_inv0[2]*v0[2]+1.936491673103709*jacob_vel_inv0[1]*nuSum[1]*v0[2]+1.118033988749895*jacob_vel_inv0[0]*nuSum[1]*v0[2]+1.936491673103709*nuSum[1]*v0[1]*jacob_vel_inv0[2]-1.58113883008419*sumNuUx[1]*jacob_vel_inv0[2]+1.118033988749895*v0[0]*nuSum[1]*jacob_vel_inv0[2]+1.5*jacob_vel_inv0[1]*nuSum[1]*v0[1]+0.8660254037844386*jacob_vel_inv0[0]*nuSum[1]*v0[1]-1.224744871391589*jacob_vel_inv0[1]*sumNuUx[1]-0.7071067811865475*jacob_vel_inv0[0]*sumNuUx[1]+0.8660254037844386*v0[0]*jacob_vel_inv0[1]*nuSum[1]+0.5*jacob_vel_inv0[0]*v0[0]*nuSum[1]; 
  alphaDrSurf_r[2] = 2.958039891549809*jacob_vel_inv0[2]*nuSum[2]*v0[3]+2.29128784747792*jacob_vel_inv0[1]*nuSum[2]*v0[3]+1.322875655532295*jacob_vel_inv0[0]*nuSum[2]*v0[3]+2.5*jacob_vel_inv0[2]*nuSum[2]*v0[2]+1.936491673103709*jacob_vel_inv0[1]*nuSum[2]*v0[2]+1.118033988749895*jacob_vel_inv0[0]*nuSum[2]*v0[2]-1.58113883008419*jacob_vel_inv0[2]*sumNuUx[2]-1.224744871391589*jacob_vel_inv0[1]*sumNuUx[2]-0.7071067811865475*jacob_vel_inv0[0]*sumNuUx[2]+1.936491673103709*v0[1]*jacob_vel_inv0[2]*nuSum[2]+1.118033988749895*v0[0]*jacob_vel_inv0[2]*nuSum[2]+1.5*jacob_vel_inv0[1]*v0[1]*nuSum[2]+0.8660254037844386*jacob_vel_inv0[0]*v0[1]*nuSum[2]+0.8660254037844386*v0[0]*jacob_vel_inv0[1]*nuSum[2]+0.5*jacob_vel_inv0[0]*v0[0]*nuSum[2]; 
  alphaDrSurf_r[3] = 2.958039891549809*jacob_vel_inv0[2]*nuSum[3]*v0[3]+2.29128784747792*jacob_vel_inv0[1]*nuSum[3]*v0[3]+1.322875655532295*jacob_vel_inv0[0]*nuSum[3]*v0[3]-1.58113883008419*jacob_vel_inv0[2]*sumNuUx[3]-1.224744871391589*jacob_vel_inv0[1]*sumNuUx[3]-0.7071067811865475*jacob_vel_inv0[0]*sumNuUx[3]+2.5*jacob_vel_inv0[2]*v0[2]*nuSum[3]+1.936491673103709*jacob_vel_inv0[1]*v0[2]*nuSum[3]+1.118033988749895*jacob_vel_inv0[0]*v0[2]*nuSum[3]+1.936491673103709*v0[1]*jacob_vel_inv0[2]*nuSum[3]+1.118033988749895*v0[0]*jacob_vel_inv0[2]*nuSum[3]+1.5*jacob_vel_inv0[1]*v0[1]*nuSum[3]+0.8660254037844386*jacob_vel_inv0[0]*v0[1]*nuSum[3]+0.8660254037844386*v0[0]*jacob_vel_inv0[1]*nuSum[3]+0.5*jacob_vel_inv0[0]*v0[0]*nuSum[3]; 
  alphaDrSurf_r[4] = (-1.58113883008419*jacob_vel_inv0[2]*sumNuUx[4])-1.224744871391589*jacob_vel_inv0[1]*sumNuUx[4]-0.7071067811865475*jacob_vel_inv0[0]*sumNuUx[4]+2.958039891549809*jacob_vel_inv0[2]*v0[3]*nuSum[4]+2.29128784747792*jacob_vel_inv0[1]*v0[3]*nuSum[4]+1.322875655532295*jacob_vel_inv0[0]*v0[3]*nuSum[4]+2.5*jacob_vel_inv0[2]*v0[2]*nuSum[4]+1.936491673103709*jacob_vel_inv0[1]*v0[2]*nuSum[4]+1.118033988749895*jacob_vel_inv0[0]*v0[2]*nuSum[4]+1.936491673103709*v0[1]*jacob_vel_inv0[2]*nuSum[4]+1.118033988749895*v0[0]*jacob_vel_inv0[2]*nuSum[4]+1.5*jacob_vel_inv0[1]*v0[1]*nuSum[4]+0.8660254037844386*jacob_vel_inv0[0]*v0[1]*nuSum[4]+0.8660254037844386*v0[0]*jacob_vel_inv0[1]*nuSum[4]+0.5*jacob_vel_inv0[0]*v0[0]*nuSum[4]; 
  alphaDrSurf_r[5] = (-1.58113883008419*jacob_vel_inv0[2]*sumNuUx[5])-1.224744871391589*jacob_vel_inv0[1]*sumNuUx[5]-0.7071067811865475*jacob_vel_inv0[0]*sumNuUx[5]+2.958039891549809*jacob_vel_inv0[2]*v0[3]*nuSum[5]+2.29128784747792*jacob_vel_inv0[1]*v0[3]*nuSum[5]+1.322875655532295*jacob_vel_inv0[0]*v0[3]*nuSum[5]+2.5*jacob_vel_inv0[2]*v0[2]*nuSum[5]+1.936491673103709*jacob_vel_inv0[1]*v0[2]*nuSum[5]+1.118033988749895*jacob_vel_inv0[0]*v0[2]*nuSum[5]+1.936491673103709*v0[1]*jacob_vel_inv0[2]*nuSum[5]+1.118033988749895*v0[0]*jacob_vel_inv0[2]*nuSum[5]+1.5*jacob_vel_inv0[1]*v0[1]*nuSum[5]+0.8660254037844386*jacob_vel_inv0[0]*v0[1]*nuSum[5]+0.8660254037844386*v0[0]*jacob_vel_inv0[1]*nuSum[5]+0.5*jacob_vel_inv0[0]*v0[0]*nuSum[5]; 
  alphaDrSurf_r[6] = (-1.58113883008419*jacob_vel_inv0[2]*sumNuUx[6])-1.224744871391589*jacob_vel_inv0[1]*sumNuUx[6]-0.7071067811865475*jacob_vel_inv0[0]*sumNuUx[6]+2.958039891549809*jacob_vel_inv0[2]*v0[3]*nuSum[6]+2.29128784747792*jacob_vel_inv0[1]*v0[3]*nuSum[6]+1.322875655532295*jacob_vel_inv0[0]*v0[3]*nuSum[6]+2.5*jacob_vel_inv0[2]*v0[2]*nuSum[6]+1.936491673103709*jacob_vel_inv0[1]*v0[2]*nuSum[6]+1.118033988749895*jacob_vel_inv0[0]*v0[2]*nuSum[6]+1.936491673103709*v0[1]*jacob_vel_inv0[2]*nuSum[6]+1.118033988749895*v0[0]*jacob_vel_inv0[2]*nuSum[6]+1.5*jacob_vel_inv0[1]*v0[1]*nuSum[6]+0.8660254037844386*jacob_vel_inv0[0]*v0[1]*nuSum[6]+0.8660254037844386*v0[0]*jacob_vel_inv0[1]*nuSum[6]+0.5*jacob_vel_inv0[0]*v0[0]*nuSum[6]; 
  alphaDrSurf_r[7] = (-1.58113883008419*jacob_vel_inv0[2]*sumNuUx[7])-1.224744871391589*jacob_vel_inv0[1]*sumNuUx[7]-0.7071067811865475*jacob_vel_inv0[0]*sumNuUx[7]+2.958039891549809*jacob_vel_inv0[2]*v0[3]*nuSum[7]+2.29128784747792*jacob_vel_inv0[1]*v0[3]*nuSum[7]+1.322875655532295*jacob_vel_inv0[0]*v0[3]*nuSum[7]+2.5*jacob_vel_inv0[2]*v0[2]*nuSum[7]+1.936491673103709*jacob_vel_inv0[1]*v0[2]*nuSum[7]+1.118033988749895*jacob_vel_inv0[0]*v0[2]*nuSum[7]+1.936491673103709*v0[1]*jacob_vel_inv0[2]*nuSum[7]+1.118033988749895*v0[0]*jacob_vel_inv0[2]*nuSum[7]+1.5*jacob_vel_inv0[1]*v0[1]*nuSum[7]+0.8660254037844386*jacob_vel_inv0[0]*v0[1]*nuSum[7]+0.8660254037844386*v0[0]*jacob_vel_inv0[1]*nuSum[7]+0.5*jacob_vel_inv0[0]*v0[0]*nuSum[7]; 

  double fUpwindQuad_l[9] = {0.0};
  double fUpwindQuad_r[9] = {0.0};
  double fUpwind_l[8] = {0.0};
  double fUpwind_r[8] = {0.0};
  double Ghat_l[8] = {0.0}; 
  double Ghat_r[8] = {0.0}; 

  if ((-0.6*alphaDrSurf_l[7])-0.6*alphaDrSurf_l[6]+0.4472135954999572*alphaDrSurf_l[5]+0.4472135954999572*alphaDrSurf_l[4]+0.9*alphaDrSurf_l[3]-0.6708203932499357*alphaDrSurf_l[2]-0.6708203932499357*alphaDrSurf_l[1]+0.5*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[0] = ser_3x_p2_surfx3_eval_quad_node_0_r(fl); 
  } else { 
    fUpwindQuad_l[0] = ser_3x_p2_surfx3_eval_quad_node_0_l(fc); 
  } 
  if ((-0.6*alphaDrSurf_r[7])-0.6*alphaDrSurf_r[6]+0.4472135954999572*alphaDrSurf_r[5]+0.4472135954999572*alphaDrSurf_r[4]+0.9*alphaDrSurf_r[3]-0.6708203932499357*alphaDrSurf_r[2]-0.6708203932499357*alphaDrSurf_r[1]+0.5*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[0] = ser_3x_p2_surfx3_eval_quad_node_0_r(fc); 
  } else { 
    fUpwindQuad_r[0] = ser_3x_p2_surfx3_eval_quad_node_0_l(fr); 
  } 
  if (0.75*alphaDrSurf_l[7]-0.5590169943749465*alphaDrSurf_l[5]+0.4472135954999572*alphaDrSurf_l[4]-0.6708203932499357*alphaDrSurf_l[1]+0.5*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[1] = ser_3x_p2_surfx3_eval_quad_node_1_r(fl); 
  } else { 
    fUpwindQuad_l[1] = ser_3x_p2_surfx3_eval_quad_node_1_l(fc); 
  } 
  if (0.75*alphaDrSurf_r[7]-0.5590169943749465*alphaDrSurf_r[5]+0.4472135954999572*alphaDrSurf_r[4]-0.6708203932499357*alphaDrSurf_r[1]+0.5*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[1] = ser_3x_p2_surfx3_eval_quad_node_1_r(fc); 
  } else { 
    fUpwindQuad_r[1] = ser_3x_p2_surfx3_eval_quad_node_1_l(fr); 
  } 
  if ((-0.6*alphaDrSurf_l[7])+0.6*alphaDrSurf_l[6]+0.4472135954999572*alphaDrSurf_l[5]+0.4472135954999572*alphaDrSurf_l[4]-0.9*alphaDrSurf_l[3]+0.6708203932499357*alphaDrSurf_l[2]-0.6708203932499357*alphaDrSurf_l[1]+0.5*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[2] = ser_3x_p2_surfx3_eval_quad_node_2_r(fl); 
  } else { 
    fUpwindQuad_l[2] = ser_3x_p2_surfx3_eval_quad_node_2_l(fc); 
  } 
  if ((-0.6*alphaDrSurf_r[7])+0.6*alphaDrSurf_r[6]+0.4472135954999572*alphaDrSurf_r[5]+0.4472135954999572*alphaDrSurf_r[4]-0.9*alphaDrSurf_r[3]+0.6708203932499357*alphaDrSurf_r[2]-0.6708203932499357*alphaDrSurf_r[1]+0.5*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[2] = ser_3x_p2_surfx3_eval_quad_node_2_r(fc); 
  } else { 
    fUpwindQuad_r[2] = ser_3x_p2_surfx3_eval_quad_node_2_l(fr); 
  } 
  if (0.75*alphaDrSurf_l[6]+0.4472135954999572*alphaDrSurf_l[5]-0.5590169943749465*alphaDrSurf_l[4]-0.6708203932499357*alphaDrSurf_l[2]+0.5*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[3] = ser_3x_p2_surfx3_eval_quad_node_3_r(fl); 
  } else { 
    fUpwindQuad_l[3] = ser_3x_p2_surfx3_eval_quad_node_3_l(fc); 
  } 
  if (0.75*alphaDrSurf_r[6]+0.4472135954999572*alphaDrSurf_r[5]-0.5590169943749465*alphaDrSurf_r[4]-0.6708203932499357*alphaDrSurf_r[2]+0.5*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[3] = ser_3x_p2_surfx3_eval_quad_node_3_r(fc); 
  } else { 
    fUpwindQuad_r[3] = ser_3x_p2_surfx3_eval_quad_node_3_l(fr); 
  } 
  if ((-0.5590169943749465*alphaDrSurf_l[5])-0.5590169943749465*alphaDrSurf_l[4]+0.5*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[4] = ser_3x_p2_surfx3_eval_quad_node_4_r(fl); 
  } else { 
    fUpwindQuad_l[4] = ser_3x_p2_surfx3_eval_quad_node_4_l(fc); 
  } 
  if ((-0.5590169943749465*alphaDrSurf_r[5])-0.5590169943749465*alphaDrSurf_r[4]+0.5*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[4] = ser_3x_p2_surfx3_eval_quad_node_4_r(fc); 
  } else { 
    fUpwindQuad_r[4] = ser_3x_p2_surfx3_eval_quad_node_4_l(fr); 
  } 
  if ((-0.75*alphaDrSurf_l[6])+0.4472135954999572*alphaDrSurf_l[5]-0.5590169943749465*alphaDrSurf_l[4]+0.6708203932499357*alphaDrSurf_l[2]+0.5*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[5] = ser_3x_p2_surfx3_eval_quad_node_5_r(fl); 
  } else { 
    fUpwindQuad_l[5] = ser_3x_p2_surfx3_eval_quad_node_5_l(fc); 
  } 
  if ((-0.75*alphaDrSurf_r[6])+0.4472135954999572*alphaDrSurf_r[5]-0.5590169943749465*alphaDrSurf_r[4]+0.6708203932499357*alphaDrSurf_r[2]+0.5*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[5] = ser_3x_p2_surfx3_eval_quad_node_5_r(fc); 
  } else { 
    fUpwindQuad_r[5] = ser_3x_p2_surfx3_eval_quad_node_5_l(fr); 
  } 
  if (0.6*alphaDrSurf_l[7]-0.6*alphaDrSurf_l[6]+0.4472135954999572*alphaDrSurf_l[5]+0.4472135954999572*alphaDrSurf_l[4]-0.9*alphaDrSurf_l[3]-0.6708203932499357*alphaDrSurf_l[2]+0.6708203932499357*alphaDrSurf_l[1]+0.5*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[6] = ser_3x_p2_surfx3_eval_quad_node_6_r(fl); 
  } else { 
    fUpwindQuad_l[6] = ser_3x_p2_surfx3_eval_quad_node_6_l(fc); 
  } 
  if (0.6*alphaDrSurf_r[7]-0.6*alphaDrSurf_r[6]+0.4472135954999572*alphaDrSurf_r[5]+0.4472135954999572*alphaDrSurf_r[4]-0.9*alphaDrSurf_r[3]-0.6708203932499357*alphaDrSurf_r[2]+0.6708203932499357*alphaDrSurf_r[1]+0.5*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[6] = ser_3x_p2_surfx3_eval_quad_node_6_r(fc); 
  } else { 
    fUpwindQuad_r[6] = ser_3x_p2_surfx3_eval_quad_node_6_l(fr); 
  } 
  if ((-0.75*alphaDrSurf_l[7])-0.5590169943749465*alphaDrSurf_l[5]+0.4472135954999572*alphaDrSurf_l[4]+0.6708203932499357*alphaDrSurf_l[1]+0.5*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[7] = ser_3x_p2_surfx3_eval_quad_node_7_r(fl); 
  } else { 
    fUpwindQuad_l[7] = ser_3x_p2_surfx3_eval_quad_node_7_l(fc); 
  } 
  if ((-0.75*alphaDrSurf_r[7])-0.5590169943749465*alphaDrSurf_r[5]+0.4472135954999572*alphaDrSurf_r[4]+0.6708203932499357*alphaDrSurf_r[1]+0.5*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[7] = ser_3x_p2_surfx3_eval_quad_node_7_r(fc); 
  } else { 
    fUpwindQuad_r[7] = ser_3x_p2_surfx3_eval_quad_node_7_l(fr); 
  } 
  if (0.6*alphaDrSurf_l[7]+0.6*alphaDrSurf_l[6]+0.4472135954999572*alphaDrSurf_l[5]+0.4472135954999572*alphaDrSurf_l[4]+0.9*alphaDrSurf_l[3]+0.6708203932499357*alphaDrSurf_l[2]+0.6708203932499357*alphaDrSurf_l[1]+0.5*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[8] = ser_3x_p2_surfx3_eval_quad_node_8_r(fl); 
  } else { 
    fUpwindQuad_l[8] = ser_3x_p2_surfx3_eval_quad_node_8_l(fc); 
  } 
  if (0.6*alphaDrSurf_r[7]+0.6*alphaDrSurf_r[6]+0.4472135954999572*alphaDrSurf_r[5]+0.4472135954999572*alphaDrSurf_r[4]+0.9*alphaDrSurf_r[3]+0.6708203932499357*alphaDrSurf_r[2]+0.6708203932499357*alphaDrSurf_r[1]+0.5*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[8] = ser_3x_p2_surfx3_eval_quad_node_8_r(fc); 
  } else { 
    fUpwindQuad_r[8] = ser_3x_p2_surfx3_eval_quad_node_8_l(fr); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_3x_p2_upwind_quad_to_modal(fUpwindQuad_l, fUpwind_l); 
  ser_3x_p2_upwind_quad_to_modal(fUpwindQuad_r, fUpwind_r); 

  Ghat_l[0] = 0.5*alphaDrSurf_l[7]*fUpwind_l[7]+0.5*alphaDrSurf_l[6]*fUpwind_l[6]+0.5*alphaDrSurf_l[5]*fUpwind_l[5]+0.5*alphaDrSurf_l[4]*fUpwind_l[4]+0.5*alphaDrSurf_l[3]*fUpwind_l[3]+0.5*alphaDrSurf_l[2]*fUpwind_l[2]+0.5*alphaDrSurf_l[1]*fUpwind_l[1]+0.5*alphaDrSurf_l[0]*fUpwind_l[0]; 
  Ghat_l[1] = 0.5000000000000001*alphaDrSurf_l[5]*fUpwind_l[7]+0.5000000000000001*fUpwind_l[5]*alphaDrSurf_l[7]+0.447213595499958*alphaDrSurf_l[3]*fUpwind_l[6]+0.447213595499958*fUpwind_l[3]*alphaDrSurf_l[6]+0.4472135954999579*alphaDrSurf_l[1]*fUpwind_l[4]+0.4472135954999579*fUpwind_l[1]*alphaDrSurf_l[4]+0.5*alphaDrSurf_l[2]*fUpwind_l[3]+0.5*fUpwind_l[2]*alphaDrSurf_l[3]+0.5*alphaDrSurf_l[0]*fUpwind_l[1]+0.5*fUpwind_l[0]*alphaDrSurf_l[1]; 
  Ghat_l[2] = 0.447213595499958*alphaDrSurf_l[3]*fUpwind_l[7]+0.447213595499958*fUpwind_l[3]*alphaDrSurf_l[7]+0.5000000000000001*alphaDrSurf_l[4]*fUpwind_l[6]+0.5000000000000001*fUpwind_l[4]*alphaDrSurf_l[6]+0.4472135954999579*alphaDrSurf_l[2]*fUpwind_l[5]+0.4472135954999579*fUpwind_l[2]*alphaDrSurf_l[5]+0.5*alphaDrSurf_l[1]*fUpwind_l[3]+0.5*fUpwind_l[1]*alphaDrSurf_l[3]+0.5*alphaDrSurf_l[0]*fUpwind_l[2]+0.5*fUpwind_l[0]*alphaDrSurf_l[2]; 
  Ghat_l[3] = 0.4*alphaDrSurf_l[6]*fUpwind_l[7]+0.447213595499958*alphaDrSurf_l[2]*fUpwind_l[7]+0.4*fUpwind_l[6]*alphaDrSurf_l[7]+0.447213595499958*fUpwind_l[2]*alphaDrSurf_l[7]+0.447213595499958*alphaDrSurf_l[1]*fUpwind_l[6]+0.447213595499958*fUpwind_l[1]*alphaDrSurf_l[6]+0.4472135954999579*alphaDrSurf_l[3]*fUpwind_l[5]+0.4472135954999579*fUpwind_l[3]*alphaDrSurf_l[5]+0.4472135954999579*alphaDrSurf_l[3]*fUpwind_l[4]+0.4472135954999579*fUpwind_l[3]*alphaDrSurf_l[4]+0.5*alphaDrSurf_l[0]*fUpwind_l[3]+0.5*fUpwind_l[0]*alphaDrSurf_l[3]+0.5*alphaDrSurf_l[1]*fUpwind_l[2]+0.5*fUpwind_l[1]*alphaDrSurf_l[2]; 
  Ghat_l[4] = 0.4472135954999579*alphaDrSurf_l[7]*fUpwind_l[7]+0.31943828249997*alphaDrSurf_l[6]*fUpwind_l[6]+0.5000000000000001*alphaDrSurf_l[2]*fUpwind_l[6]+0.5000000000000001*fUpwind_l[2]*alphaDrSurf_l[6]+0.31943828249997*alphaDrSurf_l[4]*fUpwind_l[4]+0.5*alphaDrSurf_l[0]*fUpwind_l[4]+0.5*fUpwind_l[0]*alphaDrSurf_l[4]+0.4472135954999579*alphaDrSurf_l[3]*fUpwind_l[3]+0.4472135954999579*alphaDrSurf_l[1]*fUpwind_l[1]; 
  Ghat_l[5] = 0.31943828249997*alphaDrSurf_l[7]*fUpwind_l[7]+0.5000000000000001*alphaDrSurf_l[1]*fUpwind_l[7]+0.5000000000000001*fUpwind_l[1]*alphaDrSurf_l[7]+0.4472135954999579*alphaDrSurf_l[6]*fUpwind_l[6]+0.31943828249997*alphaDrSurf_l[5]*fUpwind_l[5]+0.5*alphaDrSurf_l[0]*fUpwind_l[5]+0.5*fUpwind_l[0]*alphaDrSurf_l[5]+0.4472135954999579*alphaDrSurf_l[3]*fUpwind_l[3]+0.4472135954999579*alphaDrSurf_l[2]*fUpwind_l[2]; 
  Ghat_l[6] = 0.4*alphaDrSurf_l[3]*fUpwind_l[7]+0.4*fUpwind_l[3]*alphaDrSurf_l[7]+0.4472135954999579*alphaDrSurf_l[5]*fUpwind_l[6]+0.31943828249997*alphaDrSurf_l[4]*fUpwind_l[6]+0.5*alphaDrSurf_l[0]*fUpwind_l[6]+0.4472135954999579*fUpwind_l[5]*alphaDrSurf_l[6]+0.31943828249997*fUpwind_l[4]*alphaDrSurf_l[6]+0.5*fUpwind_l[0]*alphaDrSurf_l[6]+0.5000000000000001*alphaDrSurf_l[2]*fUpwind_l[4]+0.5000000000000001*fUpwind_l[2]*alphaDrSurf_l[4]+0.447213595499958*alphaDrSurf_l[1]*fUpwind_l[3]+0.447213595499958*fUpwind_l[1]*alphaDrSurf_l[3]; 
  Ghat_l[7] = 0.31943828249997*alphaDrSurf_l[5]*fUpwind_l[7]+0.4472135954999579*alphaDrSurf_l[4]*fUpwind_l[7]+0.5*alphaDrSurf_l[0]*fUpwind_l[7]+0.31943828249997*fUpwind_l[5]*alphaDrSurf_l[7]+0.4472135954999579*fUpwind_l[4]*alphaDrSurf_l[7]+0.5*fUpwind_l[0]*alphaDrSurf_l[7]+0.4*alphaDrSurf_l[3]*fUpwind_l[6]+0.4*fUpwind_l[3]*alphaDrSurf_l[6]+0.5000000000000001*alphaDrSurf_l[1]*fUpwind_l[5]+0.5000000000000001*fUpwind_l[1]*alphaDrSurf_l[5]+0.447213595499958*alphaDrSurf_l[2]*fUpwind_l[3]+0.447213595499958*fUpwind_l[2]*alphaDrSurf_l[3]; 

  Ghat_r[0] = 0.5*alphaDrSurf_r[7]*fUpwind_r[7]+0.5*alphaDrSurf_r[6]*fUpwind_r[6]+0.5*alphaDrSurf_r[5]*fUpwind_r[5]+0.5*alphaDrSurf_r[4]*fUpwind_r[4]+0.5*alphaDrSurf_r[3]*fUpwind_r[3]+0.5*alphaDrSurf_r[2]*fUpwind_r[2]+0.5*alphaDrSurf_r[1]*fUpwind_r[1]+0.5*alphaDrSurf_r[0]*fUpwind_r[0]; 
  Ghat_r[1] = 0.5000000000000001*alphaDrSurf_r[5]*fUpwind_r[7]+0.5000000000000001*fUpwind_r[5]*alphaDrSurf_r[7]+0.447213595499958*alphaDrSurf_r[3]*fUpwind_r[6]+0.447213595499958*fUpwind_r[3]*alphaDrSurf_r[6]+0.4472135954999579*alphaDrSurf_r[1]*fUpwind_r[4]+0.4472135954999579*fUpwind_r[1]*alphaDrSurf_r[4]+0.5*alphaDrSurf_r[2]*fUpwind_r[3]+0.5*fUpwind_r[2]*alphaDrSurf_r[3]+0.5*alphaDrSurf_r[0]*fUpwind_r[1]+0.5*fUpwind_r[0]*alphaDrSurf_r[1]; 
  Ghat_r[2] = 0.447213595499958*alphaDrSurf_r[3]*fUpwind_r[7]+0.447213595499958*fUpwind_r[3]*alphaDrSurf_r[7]+0.5000000000000001*alphaDrSurf_r[4]*fUpwind_r[6]+0.5000000000000001*fUpwind_r[4]*alphaDrSurf_r[6]+0.4472135954999579*alphaDrSurf_r[2]*fUpwind_r[5]+0.4472135954999579*fUpwind_r[2]*alphaDrSurf_r[5]+0.5*alphaDrSurf_r[1]*fUpwind_r[3]+0.5*fUpwind_r[1]*alphaDrSurf_r[3]+0.5*alphaDrSurf_r[0]*fUpwind_r[2]+0.5*fUpwind_r[0]*alphaDrSurf_r[2]; 
  Ghat_r[3] = 0.4*alphaDrSurf_r[6]*fUpwind_r[7]+0.447213595499958*alphaDrSurf_r[2]*fUpwind_r[7]+0.4*fUpwind_r[6]*alphaDrSurf_r[7]+0.447213595499958*fUpwind_r[2]*alphaDrSurf_r[7]+0.447213595499958*alphaDrSurf_r[1]*fUpwind_r[6]+0.447213595499958*fUpwind_r[1]*alphaDrSurf_r[6]+0.4472135954999579*alphaDrSurf_r[3]*fUpwind_r[5]+0.4472135954999579*fUpwind_r[3]*alphaDrSurf_r[5]+0.4472135954999579*alphaDrSurf_r[3]*fUpwind_r[4]+0.4472135954999579*fUpwind_r[3]*alphaDrSurf_r[4]+0.5*alphaDrSurf_r[0]*fUpwind_r[3]+0.5*fUpwind_r[0]*alphaDrSurf_r[3]+0.5*alphaDrSurf_r[1]*fUpwind_r[2]+0.5*fUpwind_r[1]*alphaDrSurf_r[2]; 
  Ghat_r[4] = 0.4472135954999579*alphaDrSurf_r[7]*fUpwind_r[7]+0.31943828249997*alphaDrSurf_r[6]*fUpwind_r[6]+0.5000000000000001*alphaDrSurf_r[2]*fUpwind_r[6]+0.5000000000000001*fUpwind_r[2]*alphaDrSurf_r[6]+0.31943828249997*alphaDrSurf_r[4]*fUpwind_r[4]+0.5*alphaDrSurf_r[0]*fUpwind_r[4]+0.5*fUpwind_r[0]*alphaDrSurf_r[4]+0.4472135954999579*alphaDrSurf_r[3]*fUpwind_r[3]+0.4472135954999579*alphaDrSurf_r[1]*fUpwind_r[1]; 
  Ghat_r[5] = 0.31943828249997*alphaDrSurf_r[7]*fUpwind_r[7]+0.5000000000000001*alphaDrSurf_r[1]*fUpwind_r[7]+0.5000000000000001*fUpwind_r[1]*alphaDrSurf_r[7]+0.4472135954999579*alphaDrSurf_r[6]*fUpwind_r[6]+0.31943828249997*alphaDrSurf_r[5]*fUpwind_r[5]+0.5*alphaDrSurf_r[0]*fUpwind_r[5]+0.5*fUpwind_r[0]*alphaDrSurf_r[5]+0.4472135954999579*alphaDrSurf_r[3]*fUpwind_r[3]+0.4472135954999579*alphaDrSurf_r[2]*fUpwind_r[2]; 
  Ghat_r[6] = 0.4*alphaDrSurf_r[3]*fUpwind_r[7]+0.4*fUpwind_r[3]*alphaDrSurf_r[7]+0.4472135954999579*alphaDrSurf_r[5]*fUpwind_r[6]+0.31943828249997*alphaDrSurf_r[4]*fUpwind_r[6]+0.5*alphaDrSurf_r[0]*fUpwind_r[6]+0.4472135954999579*fUpwind_r[5]*alphaDrSurf_r[6]+0.31943828249997*fUpwind_r[4]*alphaDrSurf_r[6]+0.5*fUpwind_r[0]*alphaDrSurf_r[6]+0.5000000000000001*alphaDrSurf_r[2]*fUpwind_r[4]+0.5000000000000001*fUpwind_r[2]*alphaDrSurf_r[4]+0.447213595499958*alphaDrSurf_r[1]*fUpwind_r[3]+0.447213595499958*fUpwind_r[1]*alphaDrSurf_r[3]; 
  Ghat_r[7] = 0.31943828249997*alphaDrSurf_r[5]*fUpwind_r[7]+0.4472135954999579*alphaDrSurf_r[4]*fUpwind_r[7]+0.5*alphaDrSurf_r[0]*fUpwind_r[7]+0.31943828249997*fUpwind_r[5]*alphaDrSurf_r[7]+0.4472135954999579*fUpwind_r[4]*alphaDrSurf_r[7]+0.5*fUpwind_r[0]*alphaDrSurf_r[7]+0.4*alphaDrSurf_r[3]*fUpwind_r[6]+0.4*fUpwind_r[3]*alphaDrSurf_r[6]+0.5000000000000001*alphaDrSurf_r[1]*fUpwind_r[5]+0.5000000000000001*fUpwind_r[1]*alphaDrSurf_r[5]+0.447213595499958*alphaDrSurf_r[2]*fUpwind_r[3]+0.447213595499958*fUpwind_r[2]*alphaDrSurf_r[3]; 

  out[0] += 0.7071067811865475*Ghat_r[0]*rdv2-0.7071067811865475*Ghat_l[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat_r[1]*rdv2-0.7071067811865475*Ghat_l[1]*rdv2; 
  out[2] += 0.7071067811865475*Ghat_r[2]*rdv2-0.7071067811865475*Ghat_l[2]*rdv2; 
  out[3] += 1.224744871391589*Ghat_r[0]*rdv2+1.224744871391589*Ghat_l[0]*rdv2; 
  out[4] += 0.7071067811865475*Ghat_r[3]*rdv2-0.7071067811865475*Ghat_l[3]*rdv2; 
  out[5] += 1.224744871391589*Ghat_r[1]*rdv2+1.224744871391589*Ghat_l[1]*rdv2; 
  out[6] += 1.224744871391589*Ghat_r[2]*rdv2+1.224744871391589*Ghat_l[2]*rdv2; 
  out[7] += 0.7071067811865475*Ghat_r[4]*rdv2-0.7071067811865475*Ghat_l[4]*rdv2; 
  out[8] += 0.7071067811865475*Ghat_r[5]*rdv2-0.7071067811865475*Ghat_l[5]*rdv2; 
  out[9] += 1.58113883008419*Ghat_r[0]*rdv2-1.58113883008419*Ghat_l[0]*rdv2; 
  out[10] += 1.224744871391589*Ghat_r[3]*rdv2+1.224744871391589*Ghat_l[3]*rdv2; 
  out[11] += 0.7071067811865475*Ghat_r[6]*rdv2-0.7071067811865475*Ghat_l[6]*rdv2; 
  out[12] += 0.7071067811865475*Ghat_r[7]*rdv2-0.7071067811865475*Ghat_l[7]*rdv2; 
  out[13] += 1.224744871391589*Ghat_r[4]*rdv2+1.224744871391589*Ghat_l[4]*rdv2; 
  out[14] += 1.224744871391589*Ghat_r[5]*rdv2+1.224744871391589*Ghat_l[5]*rdv2; 
  out[15] += 1.58113883008419*Ghat_r[1]*rdv2-1.58113883008419*Ghat_l[1]*rdv2; 
  out[16] += 1.58113883008419*Ghat_r[2]*rdv2-1.58113883008419*Ghat_l[2]*rdv2; 
  out[17] += 1.224744871391589*Ghat_r[6]*rdv2+1.224744871391589*Ghat_l[6]*rdv2; 
  out[18] += 1.224744871391589*Ghat_r[7]*rdv2+1.224744871391589*Ghat_l[7]*rdv2; 
  out[19] += 1.58113883008419*Ghat_r[3]*rdv2-1.58113883008419*Ghat_l[3]*rdv2; 

  return 0.;

} 
