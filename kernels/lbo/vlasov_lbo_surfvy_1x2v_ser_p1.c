#include <gkyl_vlasov_lbo_kernels.h> 
#include <gkyl_basis_ser_1x2v_p1_surfvy_quad.h> 
GKYL_CU_DH void vlasov_lbo_surfvy_1x2v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[3]:         cell-center coordinates. 
  // dxv[3]:       cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[4]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[2]: sum of thermal speeds squared time their respective collisionalities. 
  // fl/fc/fr:      distribution function in cells 
  // out:           incremented distribution function in cell 
  double rdv2 = 2.0/dxv[2]; 
  double rdvSq4 = 4.0/(dxv[2]*dxv[2]); 

  const double *sumNuUy = &nuUSum[2]; 

  double alphaDrSurf_l[4] = {0.0}; 
  alphaDrSurf_l[0] = 1.414213562373095*nuSum[0]*w[2]-0.7071067811865475*nuSum[0]*dxv[2]-1.414213562373095*sumNuUy[0]; 
  alphaDrSurf_l[1] = 1.414213562373095*nuSum[1]*w[2]-0.7071067811865475*nuSum[1]*dxv[2]-1.414213562373095*sumNuUy[1]; 

  double alphaDrSurf_r[4] = {0.0}; 
  alphaDrSurf_r[0] = 1.414213562373095*nuSum[0]*w[2]+0.7071067811865475*nuSum[0]*dxv[2]-1.414213562373095*sumNuUy[0]; 
  alphaDrSurf_r[1] = 1.414213562373095*nuSum[1]*w[2]+0.7071067811865475*nuSum[1]*dxv[2]-1.414213562373095*sumNuUy[1]; 

  double fUpwindQuad_l[4] = {0.0};
  double fUpwindQuad_r[4] = {0.0};
  double fUpwind_l[4] = {0.0};;
  double fUpwind_r[4] = {0.0};
  double Ghat_l[4] = {0.0}; 
  double Ghat_r[4] = {0.0}; 
  double Gdiff_l[4] = {0.0}; 
  double Gdiff_r[4] = {0.0}; 
  double Gdiff2_l[4] = {0.0}; 
  double Gdiff2_r[4] = {0.0}; 

  if (alphaDrSurf_l[0]-alphaDrSurf_l[1] < 0) { 
    fUpwindQuad_l[0] = ser_1x2v_p1_surfvy_quad_0(1, fl); 
  } else { 
    fUpwindQuad_l[0] = ser_1x2v_p1_surfvy_quad_0(-1, fc); 
  } 
  if (alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[1] = ser_1x2v_p1_surfvy_quad_1(1, fl); 
  } else { 
    fUpwindQuad_l[1] = ser_1x2v_p1_surfvy_quad_1(-1, fc); 
  } 
  if (alphaDrSurf_l[0]-alphaDrSurf_l[1] < 0) { 
    fUpwindQuad_l[2] = ser_1x2v_p1_surfvy_quad_2(1, fl); 
  } else { 
    fUpwindQuad_l[2] = ser_1x2v_p1_surfvy_quad_2(-1, fc); 
  } 
  if (alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[3] = ser_1x2v_p1_surfvy_quad_3(1, fl); 
  } else { 
    fUpwindQuad_l[3] = ser_1x2v_p1_surfvy_quad_3(-1, fc); 
  } 
  if (alphaDrSurf_r[0]-alphaDrSurf_r[1] < 0) { 
    fUpwindQuad_r[0] = ser_1x2v_p1_surfvy_quad_0(1, fc); 
  } else { 
    fUpwindQuad_r[0] = ser_1x2v_p1_surfvy_quad_0(-1, fr); 
  } 
  if (alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[1] = ser_1x2v_p1_surfvy_quad_1(1, fc); 
  } else { 
    fUpwindQuad_r[1] = ser_1x2v_p1_surfvy_quad_1(-1, fr); 
  } 
  if (alphaDrSurf_r[0]-alphaDrSurf_r[1] < 0) { 
    fUpwindQuad_r[2] = ser_1x2v_p1_surfvy_quad_2(1, fc); 
  } else { 
    fUpwindQuad_r[2] = ser_1x2v_p1_surfvy_quad_2(-1, fr); 
  } 
  if (alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[3] = ser_1x2v_p1_surfvy_quad_3(1, fc); 
  } else { 
    fUpwindQuad_r[3] = ser_1x2v_p1_surfvy_quad_3(-1, fr); 
  } 
  fUpwind_l[0] = 0.5*(fUpwindQuad_l[3]+fUpwindQuad_l[2]+fUpwindQuad_l[1]+fUpwindQuad_l[0]); 
  fUpwind_l[1] = 0.5*(fUpwindQuad_l[3]-1.0*fUpwindQuad_l[2]+fUpwindQuad_l[1]-1.0*fUpwindQuad_l[0]); 
  fUpwind_l[2] = 0.5*(fUpwindQuad_l[3]+fUpwindQuad_l[2]-1.0*(fUpwindQuad_l[1]+fUpwindQuad_l[0])); 
  fUpwind_l[3] = 0.5*(fUpwindQuad_l[3]-1.0*(fUpwindQuad_l[2]+fUpwindQuad_l[1])+fUpwindQuad_l[0]); 

  fUpwind_r[0] = 0.5*(fUpwindQuad_r[3]+fUpwindQuad_r[2]+fUpwindQuad_r[1]+fUpwindQuad_r[0]); 
  fUpwind_r[1] = 0.5*(fUpwindQuad_r[3]-1.0*fUpwindQuad_r[2]+fUpwindQuad_r[1]-1.0*fUpwindQuad_r[0]); 
  fUpwind_r[2] = 0.5*(fUpwindQuad_r[3]+fUpwindQuad_r[2]-1.0*(fUpwindQuad_r[1]+fUpwindQuad_r[0])); 
  fUpwind_r[3] = 0.5*(fUpwindQuad_r[3]-1.0*(fUpwindQuad_r[2]+fUpwindQuad_r[1])+fUpwindQuad_r[0]); 

  Gdiff2_l[0] = 0.2886751345948129*nuVtSqSum[1]*fl[5]-0.2886751345948129*nuVtSqSum[1]*fc[5]+0.2886751345948129*nuVtSqSum[0]*fl[3]-0.2886751345948129*nuVtSqSum[0]*fc[3]+0.25*fl[1]*nuVtSqSum[1]+0.25*fc[1]*nuVtSqSum[1]+0.25*fl[0]*nuVtSqSum[0]+0.25*fc[0]*nuVtSqSum[0]; 
  Gdiff2_l[1] = 0.2886751345948129*nuVtSqSum[0]*fl[5]-0.2886751345948129*nuVtSqSum[0]*fc[5]+0.2886751345948129*nuVtSqSum[1]*fl[3]-0.2886751345948129*nuVtSqSum[1]*fc[3]+0.25*fl[0]*nuVtSqSum[1]+0.25*fc[0]*nuVtSqSum[1]+0.25*nuVtSqSum[0]*fl[1]+0.25*nuVtSqSum[0]*fc[1]; 
  Gdiff2_l[2] = 0.2886751345948129*nuVtSqSum[1]*fl[7]-0.2886751345948129*nuVtSqSum[1]*fc[7]+0.2886751345948129*nuVtSqSum[0]*fl[6]-0.2886751345948129*nuVtSqSum[0]*fc[6]+0.25*nuVtSqSum[1]*fl[4]+0.25*nuVtSqSum[1]*fc[4]+0.25*nuVtSqSum[0]*fl[2]+0.25*nuVtSqSum[0]*fc[2]; 
  Gdiff2_l[3] = 0.2886751345948129*nuVtSqSum[0]*fl[7]-0.2886751345948129*nuVtSqSum[0]*fc[7]+0.2886751345948129*nuVtSqSum[1]*fl[6]-0.2886751345948129*nuVtSqSum[1]*fc[6]+0.25*nuVtSqSum[0]*fl[4]+0.25*nuVtSqSum[0]*fc[4]+0.25*nuVtSqSum[1]*fl[2]+0.25*nuVtSqSum[1]*fc[2]; 

  Gdiff2_r[0] = (-0.2886751345948129*nuVtSqSum[1]*fr[5])+0.2886751345948129*nuVtSqSum[1]*fc[5]-0.2886751345948129*nuVtSqSum[0]*fr[3]+0.2886751345948129*nuVtSqSum[0]*fc[3]+0.25*fr[1]*nuVtSqSum[1]+0.25*fc[1]*nuVtSqSum[1]+0.25*fr[0]*nuVtSqSum[0]+0.25*fc[0]*nuVtSqSum[0]; 
  Gdiff2_r[1] = (-0.2886751345948129*nuVtSqSum[0]*fr[5])+0.2886751345948129*nuVtSqSum[0]*fc[5]-0.2886751345948129*nuVtSqSum[1]*fr[3]+0.2886751345948129*nuVtSqSum[1]*fc[3]+0.25*fr[0]*nuVtSqSum[1]+0.25*fc[0]*nuVtSqSum[1]+0.25*nuVtSqSum[0]*fr[1]+0.25*nuVtSqSum[0]*fc[1]; 
  Gdiff2_r[2] = (-0.2886751345948129*nuVtSqSum[1]*fr[7])+0.2886751345948129*nuVtSqSum[1]*fc[7]-0.2886751345948129*nuVtSqSum[0]*fr[6]+0.2886751345948129*nuVtSqSum[0]*fc[6]+0.25*nuVtSqSum[1]*fr[4]+0.25*nuVtSqSum[1]*fc[4]+0.25*nuVtSqSum[0]*fr[2]+0.25*nuVtSqSum[0]*fc[2]; 
  Gdiff2_r[3] = (-0.2886751345948129*nuVtSqSum[0]*fr[7])+0.2886751345948129*nuVtSqSum[0]*fc[7]-0.2886751345948129*nuVtSqSum[1]*fr[6]+0.2886751345948129*nuVtSqSum[1]*fc[6]+0.25*nuVtSqSum[0]*fr[4]+0.25*nuVtSqSum[0]*fc[4]+0.25*nuVtSqSum[1]*fr[2]+0.25*nuVtSqSum[1]*fc[2]; 

  Gdiff_l[0] = (-0.5412658773652741*nuVtSqSum[1]*fl[5])-0.5412658773652741*nuVtSqSum[1]*fc[5]-0.5412658773652741*nuVtSqSum[0]*fl[3]-0.5412658773652741*nuVtSqSum[0]*fc[3]-0.5625*fl[1]*nuVtSqSum[1]+0.5625*fc[1]*nuVtSqSum[1]-0.5625*fl[0]*nuVtSqSum[0]+0.5625*fc[0]*nuVtSqSum[0]; 
  Gdiff_l[1] = (-0.5412658773652741*nuVtSqSum[0]*fl[5])-0.5412658773652741*nuVtSqSum[0]*fc[5]-0.5412658773652741*nuVtSqSum[1]*fl[3]-0.5412658773652741*nuVtSqSum[1]*fc[3]-0.5625*fl[0]*nuVtSqSum[1]+0.5625*fc[0]*nuVtSqSum[1]-0.5625*nuVtSqSum[0]*fl[1]+0.5625*nuVtSqSum[0]*fc[1]; 
  Gdiff_l[2] = (-0.5412658773652741*nuVtSqSum[1]*fl[7])-0.5412658773652741*nuVtSqSum[1]*fc[7]-0.5412658773652741*nuVtSqSum[0]*fl[6]-0.5412658773652741*nuVtSqSum[0]*fc[6]-0.5625*nuVtSqSum[1]*fl[4]+0.5625*nuVtSqSum[1]*fc[4]-0.5625*nuVtSqSum[0]*fl[2]+0.5625*nuVtSqSum[0]*fc[2]; 
  Gdiff_l[3] = (-0.5412658773652741*nuVtSqSum[0]*fl[7])-0.5412658773652741*nuVtSqSum[0]*fc[7]-0.5412658773652741*nuVtSqSum[1]*fl[6]-0.5412658773652741*nuVtSqSum[1]*fc[6]-0.5625*nuVtSqSum[0]*fl[4]+0.5625*nuVtSqSum[0]*fc[4]-0.5625*nuVtSqSum[1]*fl[2]+0.5625*nuVtSqSum[1]*fc[2]; 

  Gdiff_r[0] = (-0.5412658773652741*nuVtSqSum[1]*fr[5])-0.5412658773652741*nuVtSqSum[1]*fc[5]-0.5412658773652741*nuVtSqSum[0]*fr[3]-0.5412658773652741*nuVtSqSum[0]*fc[3]+0.5625*fr[1]*nuVtSqSum[1]-0.5625*fc[1]*nuVtSqSum[1]+0.5625*fr[0]*nuVtSqSum[0]-0.5625*fc[0]*nuVtSqSum[0]; 
  Gdiff_r[1] = (-0.5412658773652741*nuVtSqSum[0]*fr[5])-0.5412658773652741*nuVtSqSum[0]*fc[5]-0.5412658773652741*nuVtSqSum[1]*fr[3]-0.5412658773652741*nuVtSqSum[1]*fc[3]+0.5625*fr[0]*nuVtSqSum[1]-0.5625*fc[0]*nuVtSqSum[1]+0.5625*nuVtSqSum[0]*fr[1]-0.5625*nuVtSqSum[0]*fc[1]; 
  Gdiff_r[2] = (-0.5412658773652741*nuVtSqSum[1]*fr[7])-0.5412658773652741*nuVtSqSum[1]*fc[7]-0.5412658773652741*nuVtSqSum[0]*fr[6]-0.5412658773652741*nuVtSqSum[0]*fc[6]+0.5625*nuVtSqSum[1]*fr[4]-0.5625*nuVtSqSum[1]*fc[4]+0.5625*nuVtSqSum[0]*fr[2]-0.5625*nuVtSqSum[0]*fc[2]; 
  Gdiff_r[3] = (-0.5412658773652741*nuVtSqSum[0]*fr[7])-0.5412658773652741*nuVtSqSum[0]*fc[7]-0.5412658773652741*nuVtSqSum[1]*fr[6]-0.5412658773652741*nuVtSqSum[1]*fc[6]+0.5625*nuVtSqSum[0]*fr[4]-0.5625*nuVtSqSum[0]*fc[4]+0.5625*nuVtSqSum[1]*fr[2]-0.5625*nuVtSqSum[1]*fc[2]; 

  Ghat_l[0] = Gdiff_l[0]*rdv2+0.5*alphaDrSurf_l[1]*fUpwind_l[1]+0.5*alphaDrSurf_l[0]*fUpwind_l[0]; 
  Ghat_l[1] = Gdiff_l[1]*rdv2+0.5*alphaDrSurf_l[0]*fUpwind_l[1]+0.5*fUpwind_l[0]*alphaDrSurf_l[1]; 
  Ghat_l[2] = Gdiff_l[2]*rdv2+0.5*alphaDrSurf_l[1]*fUpwind_l[3]+0.5*alphaDrSurf_l[0]*fUpwind_l[2]; 
  Ghat_l[3] = Gdiff_l[3]*rdv2+0.5*alphaDrSurf_l[0]*fUpwind_l[3]+0.5*alphaDrSurf_l[1]*fUpwind_l[2]; 

  Ghat_r[0] = Gdiff_r[0]*rdv2+0.5*alphaDrSurf_r[1]*fUpwind_r[1]+0.5*alphaDrSurf_r[0]*fUpwind_r[0]; 
  Ghat_r[1] = Gdiff_r[1]*rdv2+0.5*alphaDrSurf_r[0]*fUpwind_r[1]+0.5*fUpwind_r[0]*alphaDrSurf_r[1]; 
  Ghat_r[2] = Gdiff_r[2]*rdv2+0.5*alphaDrSurf_r[1]*fUpwind_r[3]+0.5*alphaDrSurf_r[0]*fUpwind_r[2]; 
  Ghat_r[3] = Gdiff_r[3]*rdv2+0.5*alphaDrSurf_r[0]*fUpwind_r[3]+0.5*alphaDrSurf_r[1]*fUpwind_r[2]; 

  out[0] += 0.7071067811865475*Ghat_l[0]*rdv2-0.7071067811865475*Ghat_r[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat_l[1]*rdv2-0.7071067811865475*Ghat_r[1]*rdv2; 
  out[2] += 0.7071067811865475*Ghat_l[2]*rdv2-0.7071067811865475*Ghat_r[2]*rdv2; 
  out[3] += 1.224744871391589*Gdiff2_r[0]*rdvSq4-1.224744871391589*Gdiff2_l[0]*rdvSq4-1.224744871391589*Ghat_r[0]*rdv2-1.224744871391589*Ghat_l[0]*rdv2; 
  out[4] += 0.7071067811865475*Ghat_l[3]*rdv2-0.7071067811865475*Ghat_r[3]*rdv2; 
  out[5] += 1.224744871391589*Gdiff2_r[1]*rdvSq4-1.224744871391589*Gdiff2_l[1]*rdvSq4-1.224744871391589*Ghat_r[1]*rdv2-1.224744871391589*Ghat_l[1]*rdv2; 
  out[6] += 1.224744871391589*Gdiff2_r[2]*rdvSq4-1.224744871391589*Gdiff2_l[2]*rdvSq4-1.224744871391589*Ghat_r[2]*rdv2-1.224744871391589*Ghat_l[2]*rdv2; 
  out[7] += 1.224744871391589*Gdiff2_r[3]*rdvSq4-1.224744871391589*Gdiff2_l[3]*rdvSq4-1.224744871391589*Ghat_r[3]*rdv2-1.224744871391589*Ghat_l[3]*rdv2; 
} 
