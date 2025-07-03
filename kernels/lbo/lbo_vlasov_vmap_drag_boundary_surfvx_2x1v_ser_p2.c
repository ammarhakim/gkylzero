#include <gkyl_lbo_vlasov_kernels.h> 
#include <gkyl_basis_ser_3x_p2_surfx3_eval_quad.h> 
#include <gkyl_basis_ser_3x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH double lbo_vlasov_vmap_drag_boundary_surfvx_2x1v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_vel_inv, const double *nuSum, const double *nuPrimMomsSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[3]: Cell-center coordinates. 
  // dxv[3]: Cell spacing. 
  // vmap: Velocity-space nonuniform mapping in each dimension. 
  // jacob_vel_inv: Inverse of velocity space Jacobian in each dimension. 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[16]: sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 
  // fSkin/Edge: Distribution function in cells 
  // out: Incremented distribution function in cell 
  double rdv2 = 2.0/dxv[2]; 

  const double *v0 = &vmap[0]; 
  const double *jacob_vel_inv0 = &jacob_vel_inv[0]; 
  const double *sumNuUx = &nuPrimMomsSum[0]; 

  double alphaDrSurf[8] = {0.0}; 
  double fUpwindQuad[9] = {0.0};
  double fUpwind[8] = {0.0};;
  double Ghat[8] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = 0.3535533905932737*(nuSum[0]*((8.366600265340756*jacob_vel_inv0[2]+6.480740698407861*jacob_vel_inv0[1]+3.741657386773942*jacob_vel_inv0[0])*v0[3]+(7.071067811865476*jacob_vel_inv0[2]+5.477225575051662*jacob_vel_inv0[1]+3.16227766016838*jacob_vel_inv0[0])*v0[2])+(nuSum[0]*(5.477225575051662*v0[1]+3.16227766016838*v0[0])-4.47213595499958*sumNuUx[0])*jacob_vel_inv0[2]+nuSum[0]*(4.242640687119286*jacob_vel_inv0[1]+2.449489742783178*jacob_vel_inv0[0])*v0[1]+(2.449489742783178*nuSum[0]*v0[0]-3.464101615137754*sumNuUx[0])*jacob_vel_inv0[1]+jacob_vel_inv0[0]*(1.414213562373095*nuSum[0]*v0[0]-2.0*sumNuUx[0])); 
  alphaDrSurf[1] = 0.3535533905932737*(nuSum[1]*((8.366600265340756*jacob_vel_inv0[2]+6.480740698407861*jacob_vel_inv0[1]+3.741657386773942*jacob_vel_inv0[0])*v0[3]+(7.071067811865476*jacob_vel_inv0[2]+5.477225575051662*jacob_vel_inv0[1]+3.16227766016838*jacob_vel_inv0[0])*v0[2])+(5.477225575051662*nuSum[1]*v0[1]-4.47213595499958*sumNuUx[1]+3.16227766016838*v0[0]*nuSum[1])*jacob_vel_inv0[2]+(4.242640687119286*jacob_vel_inv0[1]+2.449489742783178*jacob_vel_inv0[0])*nuSum[1]*v0[1]+((-3.464101615137754*jacob_vel_inv0[1])-2.0*jacob_vel_inv0[0])*sumNuUx[1]+v0[0]*(2.449489742783178*jacob_vel_inv0[1]+1.414213562373095*jacob_vel_inv0[0])*nuSum[1]); 
  alphaDrSurf[2] = 0.3535533905932737*(nuSum[2]*((8.366600265340756*jacob_vel_inv0[2]+6.480740698407861*jacob_vel_inv0[1]+3.741657386773942*jacob_vel_inv0[0])*v0[3]+(7.071067811865476*jacob_vel_inv0[2]+5.477225575051662*jacob_vel_inv0[1]+3.16227766016838*jacob_vel_inv0[0])*v0[2])+((-4.47213595499958*jacob_vel_inv0[2])-3.464101615137754*jacob_vel_inv0[1]-2.0*jacob_vel_inv0[0])*sumNuUx[2]+((5.477225575051662*v0[1]+3.16227766016838*v0[0])*jacob_vel_inv0[2]+(4.242640687119286*jacob_vel_inv0[1]+2.449489742783178*jacob_vel_inv0[0])*v0[1]+v0[0]*(2.449489742783178*jacob_vel_inv0[1]+1.414213562373095*jacob_vel_inv0[0]))*nuSum[2]); 
  alphaDrSurf[3] = 0.3535533905932737*((8.366600265340756*jacob_vel_inv0[2]+6.480740698407861*jacob_vel_inv0[1]+3.741657386773942*jacob_vel_inv0[0])*nuSum[3]*v0[3]+((-4.47213595499958*jacob_vel_inv0[2])-3.464101615137754*jacob_vel_inv0[1]-2.0*jacob_vel_inv0[0])*sumNuUx[3]+((7.071067811865476*jacob_vel_inv0[2]+5.477225575051662*jacob_vel_inv0[1]+3.16227766016838*jacob_vel_inv0[0])*v0[2]+(5.477225575051662*v0[1]+3.16227766016838*v0[0])*jacob_vel_inv0[2]+(4.242640687119286*jacob_vel_inv0[1]+2.449489742783178*jacob_vel_inv0[0])*v0[1]+v0[0]*(2.449489742783178*jacob_vel_inv0[1]+1.414213562373095*jacob_vel_inv0[0]))*nuSum[3]); 
  alphaDrSurf[4] = -0.3535533905932737*((4.47213595499958*jacob_vel_inv0[2]+3.464101615137754*jacob_vel_inv0[1]+2.0*jacob_vel_inv0[0])*sumNuUx[4]+(((-8.366600265340756*jacob_vel_inv0[2])-6.480740698407861*jacob_vel_inv0[1]-3.741657386773942*jacob_vel_inv0[0])*v0[3]+((-7.071067811865476*jacob_vel_inv0[2])-5.477225575051662*jacob_vel_inv0[1]-3.16227766016838*jacob_vel_inv0[0])*v0[2]+((-5.477225575051662*v0[1])-3.16227766016838*v0[0])*jacob_vel_inv0[2]+((-4.242640687119286*jacob_vel_inv0[1])-2.449489742783178*jacob_vel_inv0[0])*v0[1]+v0[0]*((-2.449489742783178*jacob_vel_inv0[1])-1.414213562373095*jacob_vel_inv0[0]))*nuSum[4]); 
  alphaDrSurf[5] = -0.3535533905932737*((4.47213595499958*jacob_vel_inv0[2]+3.464101615137754*jacob_vel_inv0[1]+2.0*jacob_vel_inv0[0])*sumNuUx[5]+(((-8.366600265340756*jacob_vel_inv0[2])-6.480740698407861*jacob_vel_inv0[1]-3.741657386773942*jacob_vel_inv0[0])*v0[3]+((-7.071067811865476*jacob_vel_inv0[2])-5.477225575051662*jacob_vel_inv0[1]-3.16227766016838*jacob_vel_inv0[0])*v0[2]+((-5.477225575051662*v0[1])-3.16227766016838*v0[0])*jacob_vel_inv0[2]+((-4.242640687119286*jacob_vel_inv0[1])-2.449489742783178*jacob_vel_inv0[0])*v0[1]+v0[0]*((-2.449489742783178*jacob_vel_inv0[1])-1.414213562373095*jacob_vel_inv0[0]))*nuSum[5]); 
  alphaDrSurf[6] = -0.3535533905932737*((4.47213595499958*jacob_vel_inv0[2]+3.464101615137754*jacob_vel_inv0[1]+2.0*jacob_vel_inv0[0])*sumNuUx[6]+(((-8.366600265340756*jacob_vel_inv0[2])-6.480740698407861*jacob_vel_inv0[1]-3.741657386773942*jacob_vel_inv0[0])*v0[3]+((-7.071067811865476*jacob_vel_inv0[2])-5.477225575051662*jacob_vel_inv0[1]-3.16227766016838*jacob_vel_inv0[0])*v0[2]+((-5.477225575051662*v0[1])-3.16227766016838*v0[0])*jacob_vel_inv0[2]+((-4.242640687119286*jacob_vel_inv0[1])-2.449489742783178*jacob_vel_inv0[0])*v0[1]+v0[0]*((-2.449489742783178*jacob_vel_inv0[1])-1.414213562373095*jacob_vel_inv0[0]))*nuSum[6]); 
  alphaDrSurf[7] = -0.3535533905932737*((4.47213595499958*jacob_vel_inv0[2]+3.464101615137754*jacob_vel_inv0[1]+2.0*jacob_vel_inv0[0])*sumNuUx[7]+(((-8.366600265340756*jacob_vel_inv0[2])-6.480740698407861*jacob_vel_inv0[1]-3.741657386773942*jacob_vel_inv0[0])*v0[3]+((-7.071067811865476*jacob_vel_inv0[2])-5.477225575051662*jacob_vel_inv0[1]-3.16227766016838*jacob_vel_inv0[0])*v0[2]+((-5.477225575051662*v0[1])-3.16227766016838*v0[0])*jacob_vel_inv0[2]+((-4.242640687119286*jacob_vel_inv0[1])-2.449489742783178*jacob_vel_inv0[0])*v0[1]+v0[0]*((-2.449489742783178*jacob_vel_inv0[1])-1.414213562373095*jacob_vel_inv0[0]))*nuSum[7]); 

  if ((-0.6*alphaDrSurf[7])-0.6*alphaDrSurf[6]+0.4472135954999572*alphaDrSurf[5]+0.4472135954999572*alphaDrSurf[4]+0.9*alphaDrSurf[3]-0.6708203932499357*alphaDrSurf[2]-0.6708203932499357*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[0] = ser_3x_p2_surfx3_eval_quad_node_0_r(fSkin); 
  } else { 
    fUpwindQuad[0] = ser_3x_p2_surfx3_eval_quad_node_0_l(fEdge); 
  } 
  if (0.75*alphaDrSurf[7]-0.5590169943749465*alphaDrSurf[5]+0.4472135954999572*alphaDrSurf[4]-0.6708203932499357*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[1] = ser_3x_p2_surfx3_eval_quad_node_1_r(fSkin); 
  } else { 
    fUpwindQuad[1] = ser_3x_p2_surfx3_eval_quad_node_1_l(fEdge); 
  } 
  if ((-0.6*alphaDrSurf[7])+0.6*alphaDrSurf[6]+0.4472135954999572*alphaDrSurf[5]+0.4472135954999572*alphaDrSurf[4]-0.9*alphaDrSurf[3]+0.6708203932499357*alphaDrSurf[2]-0.6708203932499357*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[2] = ser_3x_p2_surfx3_eval_quad_node_2_r(fSkin); 
  } else { 
    fUpwindQuad[2] = ser_3x_p2_surfx3_eval_quad_node_2_l(fEdge); 
  } 
  if (0.75*alphaDrSurf[6]+0.4472135954999572*alphaDrSurf[5]-0.5590169943749465*alphaDrSurf[4]-0.6708203932499357*alphaDrSurf[2]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[3] = ser_3x_p2_surfx3_eval_quad_node_3_r(fSkin); 
  } else { 
    fUpwindQuad[3] = ser_3x_p2_surfx3_eval_quad_node_3_l(fEdge); 
  } 
  if ((-0.5590169943749465*alphaDrSurf[5])-0.5590169943749465*alphaDrSurf[4]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[4] = ser_3x_p2_surfx3_eval_quad_node_4_r(fSkin); 
  } else { 
    fUpwindQuad[4] = ser_3x_p2_surfx3_eval_quad_node_4_l(fEdge); 
  } 
  if ((-0.75*alphaDrSurf[6])+0.4472135954999572*alphaDrSurf[5]-0.5590169943749465*alphaDrSurf[4]+0.6708203932499357*alphaDrSurf[2]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[5] = ser_3x_p2_surfx3_eval_quad_node_5_r(fSkin); 
  } else { 
    fUpwindQuad[5] = ser_3x_p2_surfx3_eval_quad_node_5_l(fEdge); 
  } 
  if (0.6*alphaDrSurf[7]-0.6*alphaDrSurf[6]+0.4472135954999572*alphaDrSurf[5]+0.4472135954999572*alphaDrSurf[4]-0.9*alphaDrSurf[3]-0.6708203932499357*alphaDrSurf[2]+0.6708203932499357*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[6] = ser_3x_p2_surfx3_eval_quad_node_6_r(fSkin); 
  } else { 
    fUpwindQuad[6] = ser_3x_p2_surfx3_eval_quad_node_6_l(fEdge); 
  } 
  if ((-0.75*alphaDrSurf[7])-0.5590169943749465*alphaDrSurf[5]+0.4472135954999572*alphaDrSurf[4]+0.6708203932499357*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[7] = ser_3x_p2_surfx3_eval_quad_node_7_r(fSkin); 
  } else { 
    fUpwindQuad[7] = ser_3x_p2_surfx3_eval_quad_node_7_l(fEdge); 
  } 
  if (0.6*alphaDrSurf[7]+0.6*alphaDrSurf[6]+0.4472135954999572*alphaDrSurf[5]+0.4472135954999572*alphaDrSurf[4]+0.9*alphaDrSurf[3]+0.6708203932499357*alphaDrSurf[2]+0.6708203932499357*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[8] = ser_3x_p2_surfx3_eval_quad_node_8_r(fSkin); 
  } else { 
    fUpwindQuad[8] = ser_3x_p2_surfx3_eval_quad_node_8_l(fEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_3x_p2_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.5*alphaDrSurf[7]*fUpwind[7]+0.5*alphaDrSurf[6]*fUpwind[6]+0.5*alphaDrSurf[5]*fUpwind[5]+0.5*alphaDrSurf[4]*fUpwind[4]+0.5*alphaDrSurf[3]*fUpwind[3]+0.5*alphaDrSurf[2]*fUpwind[2]+0.5*alphaDrSurf[1]*fUpwind[1]+0.5*alphaDrSurf[0]*fUpwind[0]; 
  Ghat[1] = 0.5000000000000001*alphaDrSurf[5]*fUpwind[7]+0.5000000000000001*fUpwind[5]*alphaDrSurf[7]+0.447213595499958*alphaDrSurf[3]*fUpwind[6]+0.447213595499958*fUpwind[3]*alphaDrSurf[6]+0.4472135954999579*alphaDrSurf[1]*fUpwind[4]+0.4472135954999579*fUpwind[1]*alphaDrSurf[4]+0.5*alphaDrSurf[2]*fUpwind[3]+0.5*fUpwind[2]*alphaDrSurf[3]+0.5*alphaDrSurf[0]*fUpwind[1]+0.5*fUpwind[0]*alphaDrSurf[1]; 
  Ghat[2] = 0.447213595499958*alphaDrSurf[3]*fUpwind[7]+0.447213595499958*fUpwind[3]*alphaDrSurf[7]+0.5000000000000001*alphaDrSurf[4]*fUpwind[6]+0.5000000000000001*fUpwind[4]*alphaDrSurf[6]+0.4472135954999579*alphaDrSurf[2]*fUpwind[5]+0.4472135954999579*fUpwind[2]*alphaDrSurf[5]+0.5*alphaDrSurf[1]*fUpwind[3]+0.5*fUpwind[1]*alphaDrSurf[3]+0.5*alphaDrSurf[0]*fUpwind[2]+0.5*fUpwind[0]*alphaDrSurf[2]; 
  Ghat[3] = 0.4*alphaDrSurf[6]*fUpwind[7]+0.447213595499958*alphaDrSurf[2]*fUpwind[7]+0.4*fUpwind[6]*alphaDrSurf[7]+0.447213595499958*fUpwind[2]*alphaDrSurf[7]+0.447213595499958*alphaDrSurf[1]*fUpwind[6]+0.447213595499958*fUpwind[1]*alphaDrSurf[6]+0.4472135954999579*alphaDrSurf[3]*fUpwind[5]+0.4472135954999579*fUpwind[3]*alphaDrSurf[5]+0.4472135954999579*alphaDrSurf[3]*fUpwind[4]+0.4472135954999579*fUpwind[3]*alphaDrSurf[4]+0.5*alphaDrSurf[0]*fUpwind[3]+0.5*fUpwind[0]*alphaDrSurf[3]+0.5*alphaDrSurf[1]*fUpwind[2]+0.5*fUpwind[1]*alphaDrSurf[2]; 
  Ghat[4] = 0.4472135954999579*alphaDrSurf[7]*fUpwind[7]+0.31943828249997*alphaDrSurf[6]*fUpwind[6]+0.5000000000000001*alphaDrSurf[2]*fUpwind[6]+0.5000000000000001*fUpwind[2]*alphaDrSurf[6]+0.31943828249997*alphaDrSurf[4]*fUpwind[4]+0.5*alphaDrSurf[0]*fUpwind[4]+0.5*fUpwind[0]*alphaDrSurf[4]+0.4472135954999579*alphaDrSurf[3]*fUpwind[3]+0.4472135954999579*alphaDrSurf[1]*fUpwind[1]; 
  Ghat[5] = 0.31943828249997*alphaDrSurf[7]*fUpwind[7]+0.5000000000000001*alphaDrSurf[1]*fUpwind[7]+0.5000000000000001*fUpwind[1]*alphaDrSurf[7]+0.4472135954999579*alphaDrSurf[6]*fUpwind[6]+0.31943828249997*alphaDrSurf[5]*fUpwind[5]+0.5*alphaDrSurf[0]*fUpwind[5]+0.5*fUpwind[0]*alphaDrSurf[5]+0.4472135954999579*alphaDrSurf[3]*fUpwind[3]+0.4472135954999579*alphaDrSurf[2]*fUpwind[2]; 
  Ghat[6] = 0.4*alphaDrSurf[3]*fUpwind[7]+0.4*fUpwind[3]*alphaDrSurf[7]+0.4472135954999579*alphaDrSurf[5]*fUpwind[6]+0.31943828249997*alphaDrSurf[4]*fUpwind[6]+0.5*alphaDrSurf[0]*fUpwind[6]+0.4472135954999579*fUpwind[5]*alphaDrSurf[6]+0.31943828249997*fUpwind[4]*alphaDrSurf[6]+0.5*fUpwind[0]*alphaDrSurf[6]+0.5000000000000001*alphaDrSurf[2]*fUpwind[4]+0.5000000000000001*fUpwind[2]*alphaDrSurf[4]+0.447213595499958*alphaDrSurf[1]*fUpwind[3]+0.447213595499958*fUpwind[1]*alphaDrSurf[3]; 
  Ghat[7] = 0.31943828249997*alphaDrSurf[5]*fUpwind[7]+0.4472135954999579*alphaDrSurf[4]*fUpwind[7]+0.5*alphaDrSurf[0]*fUpwind[7]+0.31943828249997*fUpwind[5]*alphaDrSurf[7]+0.4472135954999579*fUpwind[4]*alphaDrSurf[7]+0.5*fUpwind[0]*alphaDrSurf[7]+0.4*alphaDrSurf[3]*fUpwind[6]+0.4*fUpwind[3]*alphaDrSurf[6]+0.5000000000000001*alphaDrSurf[1]*fUpwind[5]+0.5000000000000001*fUpwind[1]*alphaDrSurf[5]+0.447213595499958*alphaDrSurf[2]*fUpwind[3]+0.447213595499958*fUpwind[2]*alphaDrSurf[3]; 

  out[0] += 0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += 0.7071067811865475*Ghat[2]*rdv2; 
  out[3] += 1.224744871391589*Ghat[0]*rdv2; 
  out[4] += 0.7071067811865475*Ghat[3]*rdv2; 
  out[5] += 1.224744871391589*Ghat[1]*rdv2; 
  out[6] += 1.224744871391589*Ghat[2]*rdv2; 
  out[7] += 0.7071067811865475*Ghat[4]*rdv2; 
  out[8] += 0.7071067811865475*Ghat[5]*rdv2; 
  out[9] += 1.58113883008419*Ghat[0]*rdv2; 
  out[10] += 1.224744871391589*Ghat[3]*rdv2; 
  out[11] += 0.7071067811865475*Ghat[6]*rdv2; 
  out[12] += 0.7071067811865475*Ghat[7]*rdv2; 
  out[13] += 1.224744871391589*Ghat[4]*rdv2; 
  out[14] += 1.224744871391589*Ghat[5]*rdv2; 
  out[15] += 1.58113883008419*Ghat[1]*rdv2; 
  out[16] += 1.58113883008419*Ghat[2]*rdv2; 
  out[17] += 1.224744871391589*Ghat[6]*rdv2; 
  out[18] += 1.224744871391589*Ghat[7]*rdv2; 
  out[19] += 1.58113883008419*Ghat[3]*rdv2; 

  } else { 

  alphaDrSurf[0] = -0.3535533905932737*(nuSum[0]*((8.366600265340756*jacob_vel_inv0[2]-6.480740698407861*jacob_vel_inv0[1]+3.741657386773942*jacob_vel_inv0[0])*v0[3]+((-7.071067811865476*jacob_vel_inv0[2])+5.477225575051662*jacob_vel_inv0[1]-3.16227766016838*jacob_vel_inv0[0])*v0[2])+(nuSum[0]*(5.477225575051662*v0[1]-3.16227766016838*v0[0])+4.47213595499958*sumNuUx[0])*jacob_vel_inv0[2]+nuSum[0]*(2.449489742783178*jacob_vel_inv0[0]-4.242640687119286*jacob_vel_inv0[1])*v0[1]+(2.449489742783178*nuSum[0]*v0[0]-3.464101615137754*sumNuUx[0])*jacob_vel_inv0[1]+jacob_vel_inv0[0]*(2.0*sumNuUx[0]-1.414213562373095*nuSum[0]*v0[0])); 
  alphaDrSurf[1] = -0.3535533905932737*(nuSum[1]*((8.366600265340756*jacob_vel_inv0[2]-6.480740698407861*jacob_vel_inv0[1]+3.741657386773942*jacob_vel_inv0[0])*v0[3]+((-7.071067811865476*jacob_vel_inv0[2])+5.477225575051662*jacob_vel_inv0[1]-3.16227766016838*jacob_vel_inv0[0])*v0[2])+(5.477225575051662*nuSum[1]*v0[1]+4.47213595499958*sumNuUx[1]-3.16227766016838*v0[0]*nuSum[1])*jacob_vel_inv0[2]+(2.449489742783178*jacob_vel_inv0[0]-4.242640687119286*jacob_vel_inv0[1])*nuSum[1]*v0[1]+(2.0*jacob_vel_inv0[0]-3.464101615137754*jacob_vel_inv0[1])*sumNuUx[1]+v0[0]*(2.449489742783178*jacob_vel_inv0[1]-1.414213562373095*jacob_vel_inv0[0])*nuSum[1]); 
  alphaDrSurf[2] = -0.3535533905932737*(nuSum[2]*((8.366600265340756*jacob_vel_inv0[2]-6.480740698407861*jacob_vel_inv0[1]+3.741657386773942*jacob_vel_inv0[0])*v0[3]+((-7.071067811865476*jacob_vel_inv0[2])+5.477225575051662*jacob_vel_inv0[1]-3.16227766016838*jacob_vel_inv0[0])*v0[2])+(4.47213595499958*jacob_vel_inv0[2]-3.464101615137754*jacob_vel_inv0[1]+2.0*jacob_vel_inv0[0])*sumNuUx[2]+((5.477225575051662*v0[1]-3.16227766016838*v0[0])*jacob_vel_inv0[2]+(2.449489742783178*jacob_vel_inv0[0]-4.242640687119286*jacob_vel_inv0[1])*v0[1]+v0[0]*(2.449489742783178*jacob_vel_inv0[1]-1.414213562373095*jacob_vel_inv0[0]))*nuSum[2]); 
  alphaDrSurf[3] = -0.3535533905932737*((8.366600265340756*jacob_vel_inv0[2]-6.480740698407861*jacob_vel_inv0[1]+3.741657386773942*jacob_vel_inv0[0])*nuSum[3]*v0[3]+(4.47213595499958*jacob_vel_inv0[2]-3.464101615137754*jacob_vel_inv0[1]+2.0*jacob_vel_inv0[0])*sumNuUx[3]+(((-7.071067811865476*jacob_vel_inv0[2])+5.477225575051662*jacob_vel_inv0[1]-3.16227766016838*jacob_vel_inv0[0])*v0[2]+(5.477225575051662*v0[1]-3.16227766016838*v0[0])*jacob_vel_inv0[2]+(2.449489742783178*jacob_vel_inv0[0]-4.242640687119286*jacob_vel_inv0[1])*v0[1]+v0[0]*(2.449489742783178*jacob_vel_inv0[1]-1.414213562373095*jacob_vel_inv0[0]))*nuSum[3]); 
  alphaDrSurf[4] = -0.3535533905932737*((4.47213595499958*jacob_vel_inv0[2]-3.464101615137754*jacob_vel_inv0[1]+2.0*jacob_vel_inv0[0])*sumNuUx[4]+((8.366600265340756*jacob_vel_inv0[2]-6.480740698407861*jacob_vel_inv0[1]+3.741657386773942*jacob_vel_inv0[0])*v0[3]+((-7.071067811865476*jacob_vel_inv0[2])+5.477225575051662*jacob_vel_inv0[1]-3.16227766016838*jacob_vel_inv0[0])*v0[2]+(5.477225575051662*v0[1]-3.16227766016838*v0[0])*jacob_vel_inv0[2]+(2.449489742783178*jacob_vel_inv0[0]-4.242640687119286*jacob_vel_inv0[1])*v0[1]+v0[0]*(2.449489742783178*jacob_vel_inv0[1]-1.414213562373095*jacob_vel_inv0[0]))*nuSum[4]); 
  alphaDrSurf[5] = -0.3535533905932737*((4.47213595499958*jacob_vel_inv0[2]-3.464101615137754*jacob_vel_inv0[1]+2.0*jacob_vel_inv0[0])*sumNuUx[5]+((8.366600265340756*jacob_vel_inv0[2]-6.480740698407861*jacob_vel_inv0[1]+3.741657386773942*jacob_vel_inv0[0])*v0[3]+((-7.071067811865476*jacob_vel_inv0[2])+5.477225575051662*jacob_vel_inv0[1]-3.16227766016838*jacob_vel_inv0[0])*v0[2]+(5.477225575051662*v0[1]-3.16227766016838*v0[0])*jacob_vel_inv0[2]+(2.449489742783178*jacob_vel_inv0[0]-4.242640687119286*jacob_vel_inv0[1])*v0[1]+v0[0]*(2.449489742783178*jacob_vel_inv0[1]-1.414213562373095*jacob_vel_inv0[0]))*nuSum[5]); 
  alphaDrSurf[6] = -0.3535533905932737*((4.47213595499958*jacob_vel_inv0[2]-3.464101615137754*jacob_vel_inv0[1]+2.0*jacob_vel_inv0[0])*sumNuUx[6]+((8.366600265340756*jacob_vel_inv0[2]-6.480740698407861*jacob_vel_inv0[1]+3.741657386773942*jacob_vel_inv0[0])*v0[3]+((-7.071067811865476*jacob_vel_inv0[2])+5.477225575051662*jacob_vel_inv0[1]-3.16227766016838*jacob_vel_inv0[0])*v0[2]+(5.477225575051662*v0[1]-3.16227766016838*v0[0])*jacob_vel_inv0[2]+(2.449489742783178*jacob_vel_inv0[0]-4.242640687119286*jacob_vel_inv0[1])*v0[1]+v0[0]*(2.449489742783178*jacob_vel_inv0[1]-1.414213562373095*jacob_vel_inv0[0]))*nuSum[6]); 
  alphaDrSurf[7] = -0.3535533905932737*((4.47213595499958*jacob_vel_inv0[2]-3.464101615137754*jacob_vel_inv0[1]+2.0*jacob_vel_inv0[0])*sumNuUx[7]+((8.366600265340756*jacob_vel_inv0[2]-6.480740698407861*jacob_vel_inv0[1]+3.741657386773942*jacob_vel_inv0[0])*v0[3]+((-7.071067811865476*jacob_vel_inv0[2])+5.477225575051662*jacob_vel_inv0[1]-3.16227766016838*jacob_vel_inv0[0])*v0[2]+(5.477225575051662*v0[1]-3.16227766016838*v0[0])*jacob_vel_inv0[2]+(2.449489742783178*jacob_vel_inv0[0]-4.242640687119286*jacob_vel_inv0[1])*v0[1]+v0[0]*(2.449489742783178*jacob_vel_inv0[1]-1.414213562373095*jacob_vel_inv0[0]))*nuSum[7]); 

  if ((-0.6*alphaDrSurf[7])-0.6*alphaDrSurf[6]+0.4472135954999572*alphaDrSurf[5]+0.4472135954999572*alphaDrSurf[4]+0.9*alphaDrSurf[3]-0.6708203932499357*alphaDrSurf[2]-0.6708203932499357*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[0] = ser_3x_p2_surfx3_eval_quad_node_0_r(fEdge); 
  } else { 
    fUpwindQuad[0] = ser_3x_p2_surfx3_eval_quad_node_0_l(fSkin); 
  } 
  if (0.75*alphaDrSurf[7]-0.5590169943749465*alphaDrSurf[5]+0.4472135954999572*alphaDrSurf[4]-0.6708203932499357*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[1] = ser_3x_p2_surfx3_eval_quad_node_1_r(fEdge); 
  } else { 
    fUpwindQuad[1] = ser_3x_p2_surfx3_eval_quad_node_1_l(fSkin); 
  } 
  if ((-0.6*alphaDrSurf[7])+0.6*alphaDrSurf[6]+0.4472135954999572*alphaDrSurf[5]+0.4472135954999572*alphaDrSurf[4]-0.9*alphaDrSurf[3]+0.6708203932499357*alphaDrSurf[2]-0.6708203932499357*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[2] = ser_3x_p2_surfx3_eval_quad_node_2_r(fEdge); 
  } else { 
    fUpwindQuad[2] = ser_3x_p2_surfx3_eval_quad_node_2_l(fSkin); 
  } 
  if (0.75*alphaDrSurf[6]+0.4472135954999572*alphaDrSurf[5]-0.5590169943749465*alphaDrSurf[4]-0.6708203932499357*alphaDrSurf[2]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[3] = ser_3x_p2_surfx3_eval_quad_node_3_r(fEdge); 
  } else { 
    fUpwindQuad[3] = ser_3x_p2_surfx3_eval_quad_node_3_l(fSkin); 
  } 
  if ((-0.5590169943749465*alphaDrSurf[5])-0.5590169943749465*alphaDrSurf[4]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[4] = ser_3x_p2_surfx3_eval_quad_node_4_r(fEdge); 
  } else { 
    fUpwindQuad[4] = ser_3x_p2_surfx3_eval_quad_node_4_l(fSkin); 
  } 
  if ((-0.75*alphaDrSurf[6])+0.4472135954999572*alphaDrSurf[5]-0.5590169943749465*alphaDrSurf[4]+0.6708203932499357*alphaDrSurf[2]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[5] = ser_3x_p2_surfx3_eval_quad_node_5_r(fEdge); 
  } else { 
    fUpwindQuad[5] = ser_3x_p2_surfx3_eval_quad_node_5_l(fSkin); 
  } 
  if (0.6*alphaDrSurf[7]-0.6*alphaDrSurf[6]+0.4472135954999572*alphaDrSurf[5]+0.4472135954999572*alphaDrSurf[4]-0.9*alphaDrSurf[3]-0.6708203932499357*alphaDrSurf[2]+0.6708203932499357*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[6] = ser_3x_p2_surfx3_eval_quad_node_6_r(fEdge); 
  } else { 
    fUpwindQuad[6] = ser_3x_p2_surfx3_eval_quad_node_6_l(fSkin); 
  } 
  if ((-0.75*alphaDrSurf[7])-0.5590169943749465*alphaDrSurf[5]+0.4472135954999572*alphaDrSurf[4]+0.6708203932499357*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[7] = ser_3x_p2_surfx3_eval_quad_node_7_r(fEdge); 
  } else { 
    fUpwindQuad[7] = ser_3x_p2_surfx3_eval_quad_node_7_l(fSkin); 
  } 
  if (0.6*alphaDrSurf[7]+0.6*alphaDrSurf[6]+0.4472135954999572*alphaDrSurf[5]+0.4472135954999572*alphaDrSurf[4]+0.9*alphaDrSurf[3]+0.6708203932499357*alphaDrSurf[2]+0.6708203932499357*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[8] = ser_3x_p2_surfx3_eval_quad_node_8_r(fEdge); 
  } else { 
    fUpwindQuad[8] = ser_3x_p2_surfx3_eval_quad_node_8_l(fSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_3x_p2_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.5*alphaDrSurf[7]*fUpwind[7]+0.5*alphaDrSurf[6]*fUpwind[6]+0.5*alphaDrSurf[5]*fUpwind[5]+0.5*alphaDrSurf[4]*fUpwind[4]+0.5*alphaDrSurf[3]*fUpwind[3]+0.5*alphaDrSurf[2]*fUpwind[2]+0.5*alphaDrSurf[1]*fUpwind[1]+0.5*alphaDrSurf[0]*fUpwind[0]; 
  Ghat[1] = 0.5000000000000001*alphaDrSurf[5]*fUpwind[7]+0.5000000000000001*fUpwind[5]*alphaDrSurf[7]+0.447213595499958*alphaDrSurf[3]*fUpwind[6]+0.447213595499958*fUpwind[3]*alphaDrSurf[6]+0.4472135954999579*alphaDrSurf[1]*fUpwind[4]+0.4472135954999579*fUpwind[1]*alphaDrSurf[4]+0.5*alphaDrSurf[2]*fUpwind[3]+0.5*fUpwind[2]*alphaDrSurf[3]+0.5*alphaDrSurf[0]*fUpwind[1]+0.5*fUpwind[0]*alphaDrSurf[1]; 
  Ghat[2] = 0.447213595499958*alphaDrSurf[3]*fUpwind[7]+0.447213595499958*fUpwind[3]*alphaDrSurf[7]+0.5000000000000001*alphaDrSurf[4]*fUpwind[6]+0.5000000000000001*fUpwind[4]*alphaDrSurf[6]+0.4472135954999579*alphaDrSurf[2]*fUpwind[5]+0.4472135954999579*fUpwind[2]*alphaDrSurf[5]+0.5*alphaDrSurf[1]*fUpwind[3]+0.5*fUpwind[1]*alphaDrSurf[3]+0.5*alphaDrSurf[0]*fUpwind[2]+0.5*fUpwind[0]*alphaDrSurf[2]; 
  Ghat[3] = 0.4*alphaDrSurf[6]*fUpwind[7]+0.447213595499958*alphaDrSurf[2]*fUpwind[7]+0.4*fUpwind[6]*alphaDrSurf[7]+0.447213595499958*fUpwind[2]*alphaDrSurf[7]+0.447213595499958*alphaDrSurf[1]*fUpwind[6]+0.447213595499958*fUpwind[1]*alphaDrSurf[6]+0.4472135954999579*alphaDrSurf[3]*fUpwind[5]+0.4472135954999579*fUpwind[3]*alphaDrSurf[5]+0.4472135954999579*alphaDrSurf[3]*fUpwind[4]+0.4472135954999579*fUpwind[3]*alphaDrSurf[4]+0.5*alphaDrSurf[0]*fUpwind[3]+0.5*fUpwind[0]*alphaDrSurf[3]+0.5*alphaDrSurf[1]*fUpwind[2]+0.5*fUpwind[1]*alphaDrSurf[2]; 
  Ghat[4] = 0.4472135954999579*alphaDrSurf[7]*fUpwind[7]+0.31943828249997*alphaDrSurf[6]*fUpwind[6]+0.5000000000000001*alphaDrSurf[2]*fUpwind[6]+0.5000000000000001*fUpwind[2]*alphaDrSurf[6]+0.31943828249997*alphaDrSurf[4]*fUpwind[4]+0.5*alphaDrSurf[0]*fUpwind[4]+0.5*fUpwind[0]*alphaDrSurf[4]+0.4472135954999579*alphaDrSurf[3]*fUpwind[3]+0.4472135954999579*alphaDrSurf[1]*fUpwind[1]; 
  Ghat[5] = 0.31943828249997*alphaDrSurf[7]*fUpwind[7]+0.5000000000000001*alphaDrSurf[1]*fUpwind[7]+0.5000000000000001*fUpwind[1]*alphaDrSurf[7]+0.4472135954999579*alphaDrSurf[6]*fUpwind[6]+0.31943828249997*alphaDrSurf[5]*fUpwind[5]+0.5*alphaDrSurf[0]*fUpwind[5]+0.5*fUpwind[0]*alphaDrSurf[5]+0.4472135954999579*alphaDrSurf[3]*fUpwind[3]+0.4472135954999579*alphaDrSurf[2]*fUpwind[2]; 
  Ghat[6] = 0.4*alphaDrSurf[3]*fUpwind[7]+0.4*fUpwind[3]*alphaDrSurf[7]+0.4472135954999579*alphaDrSurf[5]*fUpwind[6]+0.31943828249997*alphaDrSurf[4]*fUpwind[6]+0.5*alphaDrSurf[0]*fUpwind[6]+0.4472135954999579*fUpwind[5]*alphaDrSurf[6]+0.31943828249997*fUpwind[4]*alphaDrSurf[6]+0.5*fUpwind[0]*alphaDrSurf[6]+0.5000000000000001*alphaDrSurf[2]*fUpwind[4]+0.5000000000000001*fUpwind[2]*alphaDrSurf[4]+0.447213595499958*alphaDrSurf[1]*fUpwind[3]+0.447213595499958*fUpwind[1]*alphaDrSurf[3]; 
  Ghat[7] = 0.31943828249997*alphaDrSurf[5]*fUpwind[7]+0.4472135954999579*alphaDrSurf[4]*fUpwind[7]+0.5*alphaDrSurf[0]*fUpwind[7]+0.31943828249997*fUpwind[5]*alphaDrSurf[7]+0.4472135954999579*fUpwind[4]*alphaDrSurf[7]+0.5*fUpwind[0]*alphaDrSurf[7]+0.4*alphaDrSurf[3]*fUpwind[6]+0.4*fUpwind[3]*alphaDrSurf[6]+0.5000000000000001*alphaDrSurf[1]*fUpwind[5]+0.5000000000000001*fUpwind[1]*alphaDrSurf[5]+0.447213595499958*alphaDrSurf[2]*fUpwind[3]+0.447213595499958*fUpwind[2]*alphaDrSurf[3]; 

  out[0] += -0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += -0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += -0.7071067811865475*Ghat[2]*rdv2; 
  out[3] += 1.224744871391589*Ghat[0]*rdv2; 
  out[4] += -0.7071067811865475*Ghat[3]*rdv2; 
  out[5] += 1.224744871391589*Ghat[1]*rdv2; 
  out[6] += 1.224744871391589*Ghat[2]*rdv2; 
  out[7] += -0.7071067811865475*Ghat[4]*rdv2; 
  out[8] += -0.7071067811865475*Ghat[5]*rdv2; 
  out[9] += -1.58113883008419*Ghat[0]*rdv2; 
  out[10] += 1.224744871391589*Ghat[3]*rdv2; 
  out[11] += -0.7071067811865475*Ghat[6]*rdv2; 
  out[12] += -0.7071067811865475*Ghat[7]*rdv2; 
  out[13] += 1.224744871391589*Ghat[4]*rdv2; 
  out[14] += 1.224744871391589*Ghat[5]*rdv2; 
  out[15] += -1.58113883008419*Ghat[1]*rdv2; 
  out[16] += -1.58113883008419*Ghat[2]*rdv2; 
  out[17] += 1.224744871391589*Ghat[6]*rdv2; 
  out[18] += 1.224744871391589*Ghat[7]*rdv2; 
  out[19] += -1.58113883008419*Ghat[3]*rdv2; 

  } 

  return 0.;

} 
