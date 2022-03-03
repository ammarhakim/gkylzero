#include <gkyl_lbo_gyrokinetic_kernels.h> 
#include <gkyl_basis_ser_1x2v_p1_surfvx_quad.h> 
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
    fUpwindQuad_l[0] = ser_1x2v_p1_surfvx_quad_0(1, fl); 
    fUpwindQuad_l[2] = ser_1x2v_p1_surfvx_quad_2(1, fl); 
  } else { 
    fUpwindQuad_l[0] = ser_1x2v_p1_surfvx_quad_0(-1, fc); 
    fUpwindQuad_l[2] = ser_1x2v_p1_surfvx_quad_2(-1, fc); 
  } 
  if (alphaDrSurf_r[0]-alphaDrSurf_r[1] < 0) { 
    fUpwindQuad_r[0] = ser_1x2v_p1_surfvx_quad_0(1, fc); 
    fUpwindQuad_r[2] = ser_1x2v_p1_surfvx_quad_2(1, fc); 
  } else { 
    fUpwindQuad_r[0] = ser_1x2v_p1_surfvx_quad_0(-1, fr); 
    fUpwindQuad_r[2] = ser_1x2v_p1_surfvx_quad_2(-1, fr); 
  } 
  if (alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[1] = ser_1x2v_p1_surfvx_quad_1(1, fl); 
    fUpwindQuad_l[3] = ser_1x2v_p1_surfvx_quad_3(1, fl); 
  } else { 
    fUpwindQuad_l[1] = ser_1x2v_p1_surfvx_quad_1(-1, fc); 
    fUpwindQuad_l[3] = ser_1x2v_p1_surfvx_quad_3(-1, fc); 
  } 
  if (alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[1] = ser_1x2v_p1_surfvx_quad_1(1, fc); 
    fUpwindQuad_r[3] = ser_1x2v_p1_surfvx_quad_3(1, fc); 
  } else { 
    fUpwindQuad_r[1] = ser_1x2v_p1_surfvx_quad_1(-1, fr); 
    fUpwindQuad_r[3] = ser_1x2v_p1_surfvx_quad_3(-1, fr); 
  } 
  fUpwind_l[0] = 0.5*(fUpwindQuad_l[3]+fUpwindQuad_l[2]+fUpwindQuad_l[1]+fUpwindQuad_l[0]); 
  fUpwind_l[1] = 0.5*(fUpwindQuad_l[3]-1.0*fUpwindQuad_l[2]+fUpwindQuad_l[1]-1.0*fUpwindQuad_l[0]); 
  fUpwind_l[2] = 0.5*(fUpwindQuad_l[3]+fUpwindQuad_l[2]-1.0*(fUpwindQuad_l[1]+fUpwindQuad_l[0])); 
  fUpwind_l[3] = 0.5*(fUpwindQuad_l[3]-1.0*(fUpwindQuad_l[2]+fUpwindQuad_l[1])+fUpwindQuad_l[0]); 

  fUpwind_r[0] = 0.5*(fUpwindQuad_r[3]+fUpwindQuad_r[2]+fUpwindQuad_r[1]+fUpwindQuad_r[0]); 
  fUpwind_r[1] = 0.5*(fUpwindQuad_r[3]-1.0*fUpwindQuad_r[2]+fUpwindQuad_r[1]-1.0*fUpwindQuad_r[0]); 
  fUpwind_r[2] = 0.5*(fUpwindQuad_r[3]+fUpwindQuad_r[2]-1.0*(fUpwindQuad_r[1]+fUpwindQuad_r[0])); 
  fUpwind_r[3] = 0.5*(fUpwindQuad_r[3]-1.0*(fUpwindQuad_r[2]+fUpwindQuad_r[1])+fUpwindQuad_r[0]); 

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
