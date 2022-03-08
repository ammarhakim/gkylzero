#include <gkyl_lbo_vlasov_kernels.h> 
#include <gkyl_basis_ser_3x_p2_surfx2_quad.h> 
#include <gkyl_basis_ser_3x_p2_upwind.h> 
GKYL_CU_DH void lbo_vlasov_drag_surfvx_1x2v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[3]:         cell-center coordinates. 
  // dxv[3]:       cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[6]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[3]: sum of thermal speeds squared time their respective collisionalities. 
  // fl/fc/fr:      distribution function in cells 
  // out:           incremented distribution function in cell 
  double rdv2 = 2.0/dxv[1]; 

  const double *sumNuUx = &nuUSum[0]; 

  double alphaDrSurf_l[8] = {0.0}; 
  alphaDrSurf_l[0] = 1.414213562373095*nuSum[0]*w[1]-0.7071067811865475*nuSum[0]*dxv[1]-1.414213562373095*sumNuUx[0]; 
  alphaDrSurf_l[1] = 1.414213562373095*nuSum[1]*w[1]-1.414213562373095*sumNuUx[1]-0.7071067811865475*dxv[1]*nuSum[1]; 
  alphaDrSurf_l[4] = (-1.414213562373095*sumNuUx[2])+1.414213562373095*w[1]*nuSum[2]-0.7071067811865475*dxv[1]*nuSum[2]; 

  double alphaDrSurf_r[8] = {0.0}; 
  alphaDrSurf_r[0] = 1.414213562373095*nuSum[0]*w[1]+0.7071067811865475*nuSum[0]*dxv[1]-1.414213562373095*sumNuUx[0]; 
  alphaDrSurf_r[1] = 1.414213562373095*nuSum[1]*w[1]-1.414213562373095*sumNuUx[1]+0.7071067811865475*dxv[1]*nuSum[1]; 
  alphaDrSurf_r[4] = (-1.414213562373095*sumNuUx[2])+1.414213562373095*w[1]*nuSum[2]+0.7071067811865475*dxv[1]*nuSum[2]; 

  double fUpwindQuad_l[9] = {0.0};
  double fUpwindQuad_r[9] = {0.0};
  double fUpwind_l[8] = {0.0};
  double fUpwind_r[8] = {0.0};
  double Gdrag_l[8] = {0.0}; 
  double Gdrag_r[8] = {0.0}; 

  if (0.4472135954999579*alphaDrSurf_l[4]-0.6708203932499369*alphaDrSurf_l[1]+0.5*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[0] = ser_3x_p2_surfx2_quad_0_r(fl); 
    fUpwindQuad_l[3] = ser_3x_p2_surfx2_quad_3_r(fl); 
    fUpwindQuad_l[6] = ser_3x_p2_surfx2_quad_6_r(fl); 
  } else { 
    fUpwindQuad_l[0] = ser_3x_p2_surfx2_quad_0_l(fc); 
    fUpwindQuad_l[3] = ser_3x_p2_surfx2_quad_3_l(fc); 
    fUpwindQuad_l[6] = ser_3x_p2_surfx2_quad_6_l(fc); 
  } 
  if (0.4472135954999579*alphaDrSurf_r[4]-0.6708203932499369*alphaDrSurf_r[1]+0.5*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[0] = ser_3x_p2_surfx2_quad_0_r(fc); 
    fUpwindQuad_r[3] = ser_3x_p2_surfx2_quad_3_r(fc); 
    fUpwindQuad_r[6] = ser_3x_p2_surfx2_quad_6_r(fc); 
  } else { 
    fUpwindQuad_r[0] = ser_3x_p2_surfx2_quad_0_l(fr); 
    fUpwindQuad_r[3] = ser_3x_p2_surfx2_quad_3_l(fr); 
    fUpwindQuad_r[6] = ser_3x_p2_surfx2_quad_6_l(fr); 
  } 
  if (0.5*alphaDrSurf_l[0]-0.5590169943749475*alphaDrSurf_l[4] < 0) { 
    fUpwindQuad_l[1] = ser_3x_p2_surfx2_quad_1_r(fl); 
    fUpwindQuad_l[4] = ser_3x_p2_surfx2_quad_4_r(fl); 
    fUpwindQuad_l[7] = ser_3x_p2_surfx2_quad_7_r(fl); 
  } else { 
    fUpwindQuad_l[1] = ser_3x_p2_surfx2_quad_1_l(fc); 
    fUpwindQuad_l[4] = ser_3x_p2_surfx2_quad_4_l(fc); 
    fUpwindQuad_l[7] = ser_3x_p2_surfx2_quad_7_l(fc); 
  } 
  if (0.5*alphaDrSurf_r[0]-0.5590169943749475*alphaDrSurf_r[4] < 0) { 
    fUpwindQuad_r[1] = ser_3x_p2_surfx2_quad_1_r(fc); 
    fUpwindQuad_r[4] = ser_3x_p2_surfx2_quad_4_r(fc); 
    fUpwindQuad_r[7] = ser_3x_p2_surfx2_quad_7_r(fc); 
  } else { 
    fUpwindQuad_r[1] = ser_3x_p2_surfx2_quad_1_l(fr); 
    fUpwindQuad_r[4] = ser_3x_p2_surfx2_quad_4_l(fr); 
    fUpwindQuad_r[7] = ser_3x_p2_surfx2_quad_7_l(fr); 
  } 
  if (0.4472135954999579*alphaDrSurf_l[4]+0.6708203932499369*alphaDrSurf_l[1]+0.5*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[2] = ser_3x_p2_surfx2_quad_2_r(fl); 
    fUpwindQuad_l[5] = ser_3x_p2_surfx2_quad_5_r(fl); 
    fUpwindQuad_l[8] = ser_3x_p2_surfx2_quad_8_r(fl); 
  } else { 
    fUpwindQuad_l[2] = ser_3x_p2_surfx2_quad_2_l(fc); 
    fUpwindQuad_l[5] = ser_3x_p2_surfx2_quad_5_l(fc); 
    fUpwindQuad_l[8] = ser_3x_p2_surfx2_quad_8_l(fc); 
  } 
  if (0.4472135954999579*alphaDrSurf_r[4]+0.6708203932499369*alphaDrSurf_r[1]+0.5*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[2] = ser_3x_p2_surfx2_quad_2_r(fc); 
    fUpwindQuad_r[5] = ser_3x_p2_surfx2_quad_5_r(fc); 
    fUpwindQuad_r[8] = ser_3x_p2_surfx2_quad_8_r(fc); 
  } else { 
    fUpwindQuad_r[2] = ser_3x_p2_surfx2_quad_2_l(fr); 
    fUpwindQuad_r[5] = ser_3x_p2_surfx2_quad_5_l(fr); 
    fUpwindQuad_r[8] = ser_3x_p2_surfx2_quad_8_l(fr); 
  } 

  // Project nodal basis back onto modal basis. 
  ser_3x_p2_upwind(fUpwindQuad_l, fUpwind_l); 
  ser_3x_p2_upwind(fUpwindQuad_r, fUpwind_r); 

  Gdrag_l[0] = 0.5*alphaDrSurf_l[4]*fUpwind_l[4]+0.5*alphaDrSurf_l[1]*fUpwind_l[1]+0.5*alphaDrSurf_l[0]*fUpwind_l[0]; 
  Gdrag_l[1] = 0.4472135954999579*alphaDrSurf_l[1]*fUpwind_l[4]+0.4472135954999579*fUpwind_l[1]*alphaDrSurf_l[4]+0.5*alphaDrSurf_l[0]*fUpwind_l[1]+0.5*fUpwind_l[0]*alphaDrSurf_l[1]; 
  Gdrag_l[2] = 0.5000000000000001*alphaDrSurf_l[4]*fUpwind_l[6]+0.5*alphaDrSurf_l[1]*fUpwind_l[3]+0.5*alphaDrSurf_l[0]*fUpwind_l[2]; 
  Gdrag_l[3] = 0.447213595499958*alphaDrSurf_l[1]*fUpwind_l[6]+0.4472135954999579*fUpwind_l[3]*alphaDrSurf_l[4]+0.5*alphaDrSurf_l[0]*fUpwind_l[3]+0.5*alphaDrSurf_l[1]*fUpwind_l[2]; 
  Gdrag_l[4] = 0.31943828249997*alphaDrSurf_l[4]*fUpwind_l[4]+0.5*alphaDrSurf_l[0]*fUpwind_l[4]+0.5*fUpwind_l[0]*alphaDrSurf_l[4]+0.4472135954999579*alphaDrSurf_l[1]*fUpwind_l[1]; 
  Gdrag_l[5] = 0.5000000000000001*alphaDrSurf_l[1]*fUpwind_l[7]+0.5*alphaDrSurf_l[0]*fUpwind_l[5]; 
  Gdrag_l[6] = 0.31943828249997*alphaDrSurf_l[4]*fUpwind_l[6]+0.5*alphaDrSurf_l[0]*fUpwind_l[6]+0.5000000000000001*fUpwind_l[2]*alphaDrSurf_l[4]+0.447213595499958*alphaDrSurf_l[1]*fUpwind_l[3]; 
  Gdrag_l[7] = 0.4472135954999579*alphaDrSurf_l[4]*fUpwind_l[7]+0.5*alphaDrSurf_l[0]*fUpwind_l[7]+0.5000000000000001*alphaDrSurf_l[1]*fUpwind_l[5]; 

  Gdrag_r[0] = 0.5*alphaDrSurf_r[4]*fUpwind_r[4]+0.5*alphaDrSurf_r[1]*fUpwind_r[1]+0.5*alphaDrSurf_r[0]*fUpwind_r[0]; 
  Gdrag_r[1] = 0.4472135954999579*alphaDrSurf_r[1]*fUpwind_r[4]+0.4472135954999579*fUpwind_r[1]*alphaDrSurf_r[4]+0.5*alphaDrSurf_r[0]*fUpwind_r[1]+0.5*fUpwind_r[0]*alphaDrSurf_r[1]; 
  Gdrag_r[2] = 0.5000000000000001*alphaDrSurf_r[4]*fUpwind_r[6]+0.5*alphaDrSurf_r[1]*fUpwind_r[3]+0.5*alphaDrSurf_r[0]*fUpwind_r[2]; 
  Gdrag_r[3] = 0.447213595499958*alphaDrSurf_r[1]*fUpwind_r[6]+0.4472135954999579*fUpwind_r[3]*alphaDrSurf_r[4]+0.5*alphaDrSurf_r[0]*fUpwind_r[3]+0.5*alphaDrSurf_r[1]*fUpwind_r[2]; 
  Gdrag_r[4] = 0.31943828249997*alphaDrSurf_r[4]*fUpwind_r[4]+0.5*alphaDrSurf_r[0]*fUpwind_r[4]+0.5*fUpwind_r[0]*alphaDrSurf_r[4]+0.4472135954999579*alphaDrSurf_r[1]*fUpwind_r[1]; 
  Gdrag_r[5] = 0.5000000000000001*alphaDrSurf_r[1]*fUpwind_r[7]+0.5*alphaDrSurf_r[0]*fUpwind_r[5]; 
  Gdrag_r[6] = 0.31943828249997*alphaDrSurf_r[4]*fUpwind_r[6]+0.5*alphaDrSurf_r[0]*fUpwind_r[6]+0.5000000000000001*fUpwind_r[2]*alphaDrSurf_r[4]+0.447213595499958*alphaDrSurf_r[1]*fUpwind_r[3]; 
  Gdrag_r[7] = 0.4472135954999579*alphaDrSurf_r[4]*fUpwind_r[7]+0.5*alphaDrSurf_r[0]*fUpwind_r[7]+0.5000000000000001*alphaDrSurf_r[1]*fUpwind_r[5]; 

  out[0] += 0.7071067811865475*Gdrag_r[0]*rdv2-0.7071067811865475*Gdrag_l[0]*rdv2; 
  out[1] += 0.7071067811865475*Gdrag_r[1]*rdv2-0.7071067811865475*Gdrag_l[1]*rdv2; 
  out[2] += 1.224744871391589*Gdrag_r[0]*rdv2+1.224744871391589*Gdrag_l[0]*rdv2; 
  out[3] += 0.7071067811865475*Gdrag_r[2]*rdv2-0.7071067811865475*Gdrag_l[2]*rdv2; 
  out[4] += 1.224744871391589*Gdrag_r[1]*rdv2+1.224744871391589*Gdrag_l[1]*rdv2; 
  out[5] += 0.7071067811865475*Gdrag_r[3]*rdv2-0.7071067811865475*Gdrag_l[3]*rdv2; 
  out[6] += 1.224744871391589*Gdrag_r[2]*rdv2+1.224744871391589*Gdrag_l[2]*rdv2; 
  out[7] += 0.7071067811865475*Gdrag_r[4]*rdv2-0.7071067811865475*Gdrag_l[4]*rdv2; 
  out[8] += 1.58113883008419*Gdrag_r[0]*rdv2-1.58113883008419*Gdrag_l[0]*rdv2; 
  out[9] += 0.7071067811865475*Gdrag_r[5]*rdv2-0.7071067811865475*Gdrag_l[5]*rdv2; 
  out[10] += 1.224744871391589*Gdrag_r[3]*rdv2+1.224744871391589*Gdrag_l[3]*rdv2; 
  out[11] += 1.224744871391589*Gdrag_r[4]*rdv2+1.224744871391589*Gdrag_l[4]*rdv2; 
  out[12] += 1.58113883008419*Gdrag_r[1]*rdv2-1.58113883008419*Gdrag_l[1]*rdv2; 
  out[13] += 0.7071067811865475*Gdrag_r[6]*rdv2-0.7071067811865475*Gdrag_l[6]*rdv2; 
  out[14] += 1.58113883008419*Gdrag_r[2]*rdv2-1.58113883008419*Gdrag_l[2]*rdv2; 
  out[15] += 0.7071067811865475*Gdrag_r[7]*rdv2-0.7071067811865475*Gdrag_l[7]*rdv2; 
  out[16] += 1.224744871391589*Gdrag_r[5]*rdv2+1.224744871391589*Gdrag_l[5]*rdv2; 
  out[17] += 1.224744871391589*Gdrag_r[6]*rdv2+1.224744871391589*Gdrag_l[6]*rdv2; 
  out[18] += 1.58113883008419*Gdrag_r[3]*rdv2-1.58113883008419*Gdrag_l[3]*rdv2; 
  out[19] += 1.224744871391589*Gdrag_r[7]*rdv2+1.224744871391589*Gdrag_l[7]*rdv2; 
} 
