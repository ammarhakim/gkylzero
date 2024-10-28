#include <gkyl_lbo_vlasov_kernels.h> 
#include <gkyl_basis_hyb_1x1v_p1_surfx2_eval_quad.h> 
#include <gkyl_basis_hyb_1x1v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double lbo_vlasov_vmap_drag_boundary_surfvx_1x1v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_vel_inv, const double *nuSum, const double *nuPrimMomsSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[2]: Cell-center coordinates. 
  // dxv[2]: Cell spacing. 
  // vmap: Velocity-space nonuniform mapping in each dimension. 
  // jacob_vel_inv: Inverse of velocity space Jacobian in each dimension. 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[4]: sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 
  // fSkin/Edge: Distribution function in cells 
  // out: Incremented distribution function in cell 
  double rdv2 = 2.0/dxv[1]; 

  const double *v0 = &vmap[0]; 
  const double *jacob_vel_inv0 = &jacob_vel_inv[0]; 
  const double *sumNuUx = &nuPrimMomsSum[0]; 

  double alphaDrSurf[2] = {0.0}; 
  double fUpwindQuad[2] = {0.0};
  double fUpwind[2] = {0.0};;
  double Ghat[2] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = 0.5*(nuSum[0]*((5.916079783099617*jacob_vel_inv0[2]+4.58257569495584*jacob_vel_inv0[1]+2.645751311064591*jacob_vel_inv0[0])*v0[3]+(5.0*jacob_vel_inv0[2]+3.872983346207417*jacob_vel_inv0[1]+2.23606797749979*jacob_vel_inv0[0])*v0[2])+(nuSum[0]*(3.872983346207417*v0[1]+2.23606797749979*v0[0])-3.16227766016838*sumNuUx[0])*jacob_vel_inv0[2]+nuSum[0]*(3.0*jacob_vel_inv0[1]+1.732050807568877*jacob_vel_inv0[0])*v0[1]+(1.732050807568877*nuSum[0]*v0[0]-2.449489742783178*sumNuUx[0])*jacob_vel_inv0[1]+jacob_vel_inv0[0]*(nuSum[0]*v0[0]-1.414213562373095*sumNuUx[0])); 
  alphaDrSurf[1] = 0.5*(nuSum[1]*((5.916079783099617*jacob_vel_inv0[2]+4.58257569495584*jacob_vel_inv0[1]+2.645751311064591*jacob_vel_inv0[0])*v0[3]+(5.0*jacob_vel_inv0[2]+3.872983346207417*jacob_vel_inv0[1]+2.23606797749979*jacob_vel_inv0[0])*v0[2])+(3.872983346207417*nuSum[1]*v0[1]-3.16227766016838*sumNuUx[1]+2.23606797749979*v0[0]*nuSum[1])*jacob_vel_inv0[2]+(3.0*jacob_vel_inv0[1]+1.732050807568877*jacob_vel_inv0[0])*nuSum[1]*v0[1]+((-2.449489742783178*jacob_vel_inv0[1])-1.414213562373095*jacob_vel_inv0[0])*sumNuUx[1]+v0[0]*(1.732050807568877*jacob_vel_inv0[1]+jacob_vel_inv0[0])*nuSum[1]); 

  if (alphaDrSurf[0]-alphaDrSurf[1] < 0) { 
    fUpwindQuad[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_r(fSkin); 
  } else { 
    fUpwindQuad[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_l(fEdge); 
  } 
  if (alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_r(fSkin); 
  } else { 
    fUpwindQuad[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_l(fEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_1x1v_p1_vdir_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.7071067811865475*(alphaDrSurf[1]*fUpwind[1]+alphaDrSurf[0]*fUpwind[0]); 
  Ghat[1] = 0.7071067811865475*(alphaDrSurf[0]*fUpwind[1]+fUpwind[0]*alphaDrSurf[1]); 

  out[0] += 0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += 1.224744871391589*Ghat[0]*rdv2; 
  out[3] += 1.224744871391589*Ghat[1]*rdv2; 
  out[4] += 1.58113883008419*Ghat[0]*rdv2; 
  out[5] += 1.58113883008419*Ghat[1]*rdv2; 

  } else { 

  alphaDrSurf[0] = -0.5*(nuSum[0]*((5.916079783099617*jacob_vel_inv0[2]-4.58257569495584*jacob_vel_inv0[1]+2.645751311064591*jacob_vel_inv0[0])*v0[3]+((-5.0*jacob_vel_inv0[2])+3.872983346207417*jacob_vel_inv0[1]-2.23606797749979*jacob_vel_inv0[0])*v0[2])+(nuSum[0]*(3.872983346207417*v0[1]-2.23606797749979*v0[0])+3.16227766016838*sumNuUx[0])*jacob_vel_inv0[2]+nuSum[0]*(1.732050807568877*jacob_vel_inv0[0]-3.0*jacob_vel_inv0[1])*v0[1]+(1.732050807568877*nuSum[0]*v0[0]-2.449489742783178*sumNuUx[0])*jacob_vel_inv0[1]+jacob_vel_inv0[0]*(1.414213562373095*sumNuUx[0]-1.0*nuSum[0]*v0[0])); 
  alphaDrSurf[1] = -0.5*(nuSum[1]*((5.916079783099617*jacob_vel_inv0[2]-4.58257569495584*jacob_vel_inv0[1]+2.645751311064591*jacob_vel_inv0[0])*v0[3]+((-5.0*jacob_vel_inv0[2])+3.872983346207417*jacob_vel_inv0[1]-2.23606797749979*jacob_vel_inv0[0])*v0[2])+(3.872983346207417*nuSum[1]*v0[1]+3.16227766016838*sumNuUx[1]-2.23606797749979*v0[0]*nuSum[1])*jacob_vel_inv0[2]+(1.732050807568877*jacob_vel_inv0[0]-3.0*jacob_vel_inv0[1])*nuSum[1]*v0[1]+(1.414213562373095*jacob_vel_inv0[0]-2.449489742783178*jacob_vel_inv0[1])*sumNuUx[1]+v0[0]*(1.732050807568877*jacob_vel_inv0[1]-1.0*jacob_vel_inv0[0])*nuSum[1]); 

  if (alphaDrSurf[0]-alphaDrSurf[1] < 0) { 
    fUpwindQuad[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_r(fEdge); 
  } else { 
    fUpwindQuad[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_l(fSkin); 
  } 
  if (alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_r(fEdge); 
  } else { 
    fUpwindQuad[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_l(fSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_1x1v_p1_vdir_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.7071067811865475*(alphaDrSurf[1]*fUpwind[1]+alphaDrSurf[0]*fUpwind[0]); 
  Ghat[1] = 0.7071067811865475*(alphaDrSurf[0]*fUpwind[1]+fUpwind[0]*alphaDrSurf[1]); 

  out[0] += -0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += -0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += 1.224744871391589*Ghat[0]*rdv2; 
  out[3] += 1.224744871391589*Ghat[1]*rdv2; 
  out[4] += -1.58113883008419*Ghat[0]*rdv2; 
  out[5] += -1.58113883008419*Ghat[1]*rdv2; 

  } 

  return 0.;

} 
