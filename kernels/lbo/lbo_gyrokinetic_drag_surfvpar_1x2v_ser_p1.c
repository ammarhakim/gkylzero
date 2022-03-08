#include <gkyl_lbo_gyrokinetic_kernels.h> 
#include <gkyl_basis_ser_3x_p1_surfx2_quad.h> 
#include <gkyl_basis_ser_3x_p1_upwind.h> 
GKYL_CU_DH void lbo_gyrokinetic_drag_surfvpar_1x2v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[3]:     cell-center coordinates. 
  // dxv[3]:   cell spacing. 
  // m_:        species mass.
  // bmag_inv:  1/(magnetic field magnitude). 
  // nuSum:     collisionalities added (self and cross species collisionalities). 
  // nuUSum[4]:sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum: sum of thermal speeds squared time their respective collisionalities. 
  // fl/fc/fr:  distribution function in cells 
  // out:       incremented distribution function in cell 
  double rdv2 = 2.0/dxv[1]; 

  double alphaDrSurf_l[4] = {0.0}; 
  alphaDrSurf_l[0] = 1.414213562373095*nuSum[0]*w[1]-0.7071067811865475*nuSum[0]*dxv[1]-1.414213562373095*nuUSum[0]; 
  alphaDrSurf_l[1] = 1.414213562373095*nuSum[1]*w[1]-1.414213562373095*nuUSum[1]-0.7071067811865475*dxv[1]*nuSum[1]; 

  double alphaDrSurf_r[4] = {0.0}; 
  alphaDrSurf_r[0] = 1.414213562373095*nuSum[0]*w[1]+0.7071067811865475*nuSum[0]*dxv[1]-1.414213562373095*nuUSum[0]; 
  alphaDrSurf_r[1] = 1.414213562373095*nuSum[1]*w[1]-1.414213562373095*nuUSum[1]+0.7071067811865475*dxv[1]*nuSum[1]; 

  double fUpwindQuad_l[4] = {0.0};
  double fUpwindQuad_r[4] = {0.0};
  double fUpwind_l[4] = {0.0};
  double fUpwind_r[4] = {0.0};
  double Gdrag_l[4] = {0.0}; 
  double Gdrag_r[4] = {0.0}; 

  if (alphaDrSurf_l[0]-alphaDrSurf_l[1] < 0) { 
    fUpwindQuad_l[0] = ser_3x_p1_surfx2_quad_0_r(fl); 
    fUpwindQuad_l[2] = ser_3x_p1_surfx2_quad_2_r(fl); 
  } else { 
    fUpwindQuad_l[0] = ser_3x_p1_surfx2_quad_0_l(fc); 
    fUpwindQuad_l[2] = ser_3x_p1_surfx2_quad_2_l(fc); 
  } 
  if (alphaDrSurf_r[0]-alphaDrSurf_r[1] < 0) { 
    fUpwindQuad_r[0] = ser_3x_p1_surfx2_quad_0_r(fc); 
    fUpwindQuad_r[2] = ser_3x_p1_surfx2_quad_2_r(fc); 
  } else { 
    fUpwindQuad_r[0] = ser_3x_p1_surfx2_quad_0_l(fr); 
    fUpwindQuad_r[2] = ser_3x_p1_surfx2_quad_2_l(fr); 
  } 
  if (alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[1] = ser_3x_p1_surfx2_quad_1_r(fl); 
    fUpwindQuad_l[3] = ser_3x_p1_surfx2_quad_3_r(fl); 
  } else { 
    fUpwindQuad_l[1] = ser_3x_p1_surfx2_quad_1_l(fc); 
    fUpwindQuad_l[3] = ser_3x_p1_surfx2_quad_3_l(fc); 
  } 
  if (alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[1] = ser_3x_p1_surfx2_quad_1_r(fc); 
    fUpwindQuad_r[3] = ser_3x_p1_surfx2_quad_3_r(fc); 
  } else { 
    fUpwindQuad_r[1] = ser_3x_p1_surfx2_quad_1_l(fr); 
    fUpwindQuad_r[3] = ser_3x_p1_surfx2_quad_3_l(fr); 
  } 

  // Project nodal basis back onto modal basis. 
  ser_3x_p1_upwind(fUpwindQuad_l, fUpwind_l); 
  ser_3x_p1_upwind(fUpwindQuad_r, fUpwind_r); 

  Gdrag_l[0] = 0.5*alphaDrSurf_l[1]*fUpwind_l[1]+0.5*alphaDrSurf_l[0]*fUpwind_l[0]; 
  Gdrag_l[1] = 0.5*alphaDrSurf_l[0]*fUpwind_l[1]+0.5*fUpwind_l[0]*alphaDrSurf_l[1]; 
  Gdrag_l[2] = 0.5*alphaDrSurf_l[1]*fUpwind_l[3]+0.5*alphaDrSurf_l[0]*fUpwind_l[2]; 
  Gdrag_l[3] = 0.5*alphaDrSurf_l[0]*fUpwind_l[3]+0.5*alphaDrSurf_l[1]*fUpwind_l[2]; 

  Gdrag_r[0] = 0.5*alphaDrSurf_r[1]*fUpwind_r[1]+0.5*alphaDrSurf_r[0]*fUpwind_r[0]; 
  Gdrag_r[1] = 0.5*alphaDrSurf_r[0]*fUpwind_r[1]+0.5*fUpwind_r[0]*alphaDrSurf_r[1]; 
  Gdrag_r[2] = 0.5*alphaDrSurf_r[1]*fUpwind_r[3]+0.5*alphaDrSurf_r[0]*fUpwind_r[2]; 
  Gdrag_r[3] = 0.5*alphaDrSurf_r[0]*fUpwind_r[3]+0.5*alphaDrSurf_r[1]*fUpwind_r[2]; 

  out[0] += 0.7071067811865475*Gdrag_r[0]*rdv2-0.7071067811865475*Gdrag_l[0]*rdv2; 
  out[1] += 0.7071067811865475*Gdrag_r[1]*rdv2-0.7071067811865475*Gdrag_l[1]*rdv2; 
  out[2] += 1.224744871391589*Gdrag_r[0]*rdv2+1.224744871391589*Gdrag_l[0]*rdv2; 
  out[3] += 0.7071067811865475*Gdrag_r[2]*rdv2-0.7071067811865475*Gdrag_l[2]*rdv2; 
  out[4] += 1.224744871391589*Gdrag_r[1]*rdv2+1.224744871391589*Gdrag_l[1]*rdv2; 
  out[5] += 0.7071067811865475*Gdrag_r[3]*rdv2-0.7071067811865475*Gdrag_l[3]*rdv2; 
  out[6] += 1.224744871391589*Gdrag_r[2]*rdv2+1.224744871391589*Gdrag_l[2]*rdv2; 
  out[7] += 1.224744871391589*Gdrag_r[3]*rdv2+1.224744871391589*Gdrag_l[3]*rdv2; 
} 
