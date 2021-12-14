#include <gkyl_vlasov_lbo_kernels.h> 
#include <gkyl_basis_ser_2x2v_p1_surfvx_quad.h> 
GKYL_CU_DH void vlasov_lbo_surfvx_2x2v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[4]:         cell-center coordinates. 
  // dxv[4]:       cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[8]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[4]: sum of thermal speeds squared time their respective collisionalities. 
  // fl/fc/fr:      distribution function in cells 
  // out:           incremented distribution function in cell 
  double rdv2 = 2.0/dxv[2]; 
  double rdvSq4 = 4.0/(dxv[2]*dxv[2]); 

  const double *sumNuUx = &nuUSum[0]; 

  double alphaDrSurf_l[8] = {0.0}; 
  alphaDrSurf_l[0] = 1.414213562373095*nuSum[0]*w[2]-0.7071067811865475*nuSum[0]*dxv[2]-1.414213562373095*sumNuUx[0]; 
  alphaDrSurf_l[1] = 1.414213562373095*nuSum[1]*w[2]-0.7071067811865475*nuSum[1]*dxv[2]-1.414213562373095*sumNuUx[1]; 
  alphaDrSurf_l[2] = 1.414213562373095*nuSum[2]*w[2]-1.414213562373095*sumNuUx[2]-0.7071067811865475*dxv[2]*nuSum[2]; 
  alphaDrSurf_l[4] = (-1.414213562373095*sumNuUx[3])+1.414213562373095*w[2]*nuSum[3]-0.7071067811865475*dxv[2]*nuSum[3]; 

  double alphaDrSurf_r[8] = {0.0}; 
  alphaDrSurf_r[0] = 1.414213562373095*nuSum[0]*w[2]+0.7071067811865475*nuSum[0]*dxv[2]-1.414213562373095*sumNuUx[0]; 
  alphaDrSurf_r[1] = 1.414213562373095*nuSum[1]*w[2]+0.7071067811865475*nuSum[1]*dxv[2]-1.414213562373095*sumNuUx[1]; 
  alphaDrSurf_r[2] = 1.414213562373095*nuSum[2]*w[2]-1.414213562373095*sumNuUx[2]+0.7071067811865475*dxv[2]*nuSum[2]; 
  alphaDrSurf_r[4] = (-1.414213562373095*sumNuUx[3])+1.414213562373095*w[2]*nuSum[3]+0.7071067811865475*dxv[2]*nuSum[3]; 

  double fUpwindQuad_l[8] = {0.0};
  double fUpwindQuad_r[8] = {0.0};
  double fUpwind_l[8] = {0.0};;
  double fUpwind_r[8] = {0.0};
  double Gdrag_l[8] = {0.0}; 
  double Gdrag_r[8] = {0.0}; 

  if (alphaDrSurf_l[4]-alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[0] = ser_2x2v_p1_surfvx_quad_0(1, fl); 
  } else { 
    fUpwindQuad_l[0] = ser_2x2v_p1_surfvx_quad_0(-1, fc); 
  } 
  if ((-alphaDrSurf_l[4])-alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[1] = ser_2x2v_p1_surfvx_quad_1(1, fl); 
  } else { 
    fUpwindQuad_l[1] = ser_2x2v_p1_surfvx_quad_1(-1, fc); 
  } 
  if ((-alphaDrSurf_l[4])+alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[2] = ser_2x2v_p1_surfvx_quad_2(1, fl); 
  } else { 
    fUpwindQuad_l[2] = ser_2x2v_p1_surfvx_quad_2(-1, fc); 
  } 
  if (alphaDrSurf_l[4]+alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[3] = ser_2x2v_p1_surfvx_quad_3(1, fl); 
  } else { 
    fUpwindQuad_l[3] = ser_2x2v_p1_surfvx_quad_3(-1, fc); 
  } 
  if (alphaDrSurf_l[4]-alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[4] = ser_2x2v_p1_surfvx_quad_4(1, fl); 
  } else { 
    fUpwindQuad_l[4] = ser_2x2v_p1_surfvx_quad_4(-1, fc); 
  } 
  if ((-alphaDrSurf_l[4])-alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[5] = ser_2x2v_p1_surfvx_quad_5(1, fl); 
  } else { 
    fUpwindQuad_l[5] = ser_2x2v_p1_surfvx_quad_5(-1, fc); 
  } 
  if ((-alphaDrSurf_l[4])+alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[6] = ser_2x2v_p1_surfvx_quad_6(1, fl); 
  } else { 
    fUpwindQuad_l[6] = ser_2x2v_p1_surfvx_quad_6(-1, fc); 
  } 
  if (alphaDrSurf_l[4]+alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[7] = ser_2x2v_p1_surfvx_quad_7(1, fl); 
  } else { 
    fUpwindQuad_l[7] = ser_2x2v_p1_surfvx_quad_7(-1, fc); 
  } 
  if (alphaDrSurf_r[4]-alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[0] = ser_2x2v_p1_surfvx_quad_0(1, fc); 
  } else { 
    fUpwindQuad_r[0] = ser_2x2v_p1_surfvx_quad_0(-1, fr); 
  } 
  if ((-alphaDrSurf_r[4])-alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[1] = ser_2x2v_p1_surfvx_quad_1(1, fc); 
  } else { 
    fUpwindQuad_r[1] = ser_2x2v_p1_surfvx_quad_1(-1, fr); 
  } 
  if ((-alphaDrSurf_r[4])+alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[2] = ser_2x2v_p1_surfvx_quad_2(1, fc); 
  } else { 
    fUpwindQuad_r[2] = ser_2x2v_p1_surfvx_quad_2(-1, fr); 
  } 
  if (alphaDrSurf_r[4]+alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[3] = ser_2x2v_p1_surfvx_quad_3(1, fc); 
  } else { 
    fUpwindQuad_r[3] = ser_2x2v_p1_surfvx_quad_3(-1, fr); 
  } 
  if (alphaDrSurf_r[4]-alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[4] = ser_2x2v_p1_surfvx_quad_4(1, fc); 
  } else { 
    fUpwindQuad_r[4] = ser_2x2v_p1_surfvx_quad_4(-1, fr); 
  } 
  if ((-alphaDrSurf_r[4])-alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[5] = ser_2x2v_p1_surfvx_quad_5(1, fc); 
  } else { 
    fUpwindQuad_r[5] = ser_2x2v_p1_surfvx_quad_5(-1, fr); 
  } 
  if ((-alphaDrSurf_r[4])+alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[6] = ser_2x2v_p1_surfvx_quad_6(1, fc); 
  } else { 
    fUpwindQuad_r[6] = ser_2x2v_p1_surfvx_quad_6(-1, fr); 
  } 
  if (alphaDrSurf_r[4]+alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[7] = ser_2x2v_p1_surfvx_quad_7(1, fc); 
  } else { 
    fUpwindQuad_r[7] = ser_2x2v_p1_surfvx_quad_7(-1, fr); 
  } 
  fUpwind_l[0] = 0.3535533905932737*(fUpwindQuad_l[7]+fUpwindQuad_l[6]+fUpwindQuad_l[5]+fUpwindQuad_l[4]+fUpwindQuad_l[3]+fUpwindQuad_l[2]+fUpwindQuad_l[1]+fUpwindQuad_l[0]); 
  fUpwind_l[1] = 0.3535533905932737*(fUpwindQuad_l[7]-1.0*fUpwindQuad_l[6]+fUpwindQuad_l[5]-1.0*fUpwindQuad_l[4]+fUpwindQuad_l[3]-1.0*fUpwindQuad_l[2]+fUpwindQuad_l[1]-1.0*fUpwindQuad_l[0]); 
  fUpwind_l[2] = 0.3535533905932737*(fUpwindQuad_l[7]+fUpwindQuad_l[6]-1.0*(fUpwindQuad_l[5]+fUpwindQuad_l[4])+fUpwindQuad_l[3]+fUpwindQuad_l[2]-1.0*(fUpwindQuad_l[1]+fUpwindQuad_l[0])); 
  fUpwind_l[3] = 0.3535533905932737*(fUpwindQuad_l[7]+fUpwindQuad_l[6]+fUpwindQuad_l[5]+fUpwindQuad_l[4]-1.0*(fUpwindQuad_l[3]+fUpwindQuad_l[2]+fUpwindQuad_l[1]+fUpwindQuad_l[0])); 
  fUpwind_l[4] = 0.3535533905932737*(fUpwindQuad_l[7]-1.0*(fUpwindQuad_l[6]+fUpwindQuad_l[5])+fUpwindQuad_l[4]+fUpwindQuad_l[3]-1.0*(fUpwindQuad_l[2]+fUpwindQuad_l[1])+fUpwindQuad_l[0]); 
  fUpwind_l[5] = 0.3535533905932737*(fUpwindQuad_l[7]-1.0*fUpwindQuad_l[6]+fUpwindQuad_l[5]-1.0*(fUpwindQuad_l[4]+fUpwindQuad_l[3])+fUpwindQuad_l[2]-1.0*fUpwindQuad_l[1]+fUpwindQuad_l[0]); 
  fUpwind_l[6] = 0.3535533905932737*(fUpwindQuad_l[7]+fUpwindQuad_l[6]-1.0*(fUpwindQuad_l[5]+fUpwindQuad_l[4]+fUpwindQuad_l[3]+fUpwindQuad_l[2])+fUpwindQuad_l[1]+fUpwindQuad_l[0]); 
  fUpwind_l[7] = 0.3535533905932737*(fUpwindQuad_l[7]-1.0*(fUpwindQuad_l[6]+fUpwindQuad_l[5])+fUpwindQuad_l[4]-1.0*fUpwindQuad_l[3]+fUpwindQuad_l[2]+fUpwindQuad_l[1]-1.0*fUpwindQuad_l[0]); 

  fUpwind_r[0] = 0.3535533905932737*(fUpwindQuad_r[7]+fUpwindQuad_r[6]+fUpwindQuad_r[5]+fUpwindQuad_r[4]+fUpwindQuad_r[3]+fUpwindQuad_r[2]+fUpwindQuad_r[1]+fUpwindQuad_r[0]); 
  fUpwind_r[1] = 0.3535533905932737*(fUpwindQuad_r[7]-1.0*fUpwindQuad_r[6]+fUpwindQuad_r[5]-1.0*fUpwindQuad_r[4]+fUpwindQuad_r[3]-1.0*fUpwindQuad_r[2]+fUpwindQuad_r[1]-1.0*fUpwindQuad_r[0]); 
  fUpwind_r[2] = 0.3535533905932737*(fUpwindQuad_r[7]+fUpwindQuad_r[6]-1.0*(fUpwindQuad_r[5]+fUpwindQuad_r[4])+fUpwindQuad_r[3]+fUpwindQuad_r[2]-1.0*(fUpwindQuad_r[1]+fUpwindQuad_r[0])); 
  fUpwind_r[3] = 0.3535533905932737*(fUpwindQuad_r[7]+fUpwindQuad_r[6]+fUpwindQuad_r[5]+fUpwindQuad_r[4]-1.0*(fUpwindQuad_r[3]+fUpwindQuad_r[2]+fUpwindQuad_r[1]+fUpwindQuad_r[0])); 
  fUpwind_r[4] = 0.3535533905932737*(fUpwindQuad_r[7]-1.0*(fUpwindQuad_r[6]+fUpwindQuad_r[5])+fUpwindQuad_r[4]+fUpwindQuad_r[3]-1.0*(fUpwindQuad_r[2]+fUpwindQuad_r[1])+fUpwindQuad_r[0]); 
  fUpwind_r[5] = 0.3535533905932737*(fUpwindQuad_r[7]-1.0*fUpwindQuad_r[6]+fUpwindQuad_r[5]-1.0*(fUpwindQuad_r[4]+fUpwindQuad_r[3])+fUpwindQuad_r[2]-1.0*fUpwindQuad_r[1]+fUpwindQuad_r[0]); 
  fUpwind_r[6] = 0.3535533905932737*(fUpwindQuad_r[7]+fUpwindQuad_r[6]-1.0*(fUpwindQuad_r[5]+fUpwindQuad_r[4]+fUpwindQuad_r[3]+fUpwindQuad_r[2])+fUpwindQuad_r[1]+fUpwindQuad_r[0]); 
  fUpwind_r[7] = 0.3535533905932737*(fUpwindQuad_r[7]-1.0*(fUpwindQuad_r[6]+fUpwindQuad_r[5])+fUpwindQuad_r[4]-1.0*fUpwindQuad_r[3]+fUpwindQuad_r[2]+fUpwindQuad_r[1]-1.0*fUpwindQuad_r[0]); 

  Gdrag_l[0] = 0.3535533905932737*alphaDrSurf_l[4]*fUpwind_l[4]+0.3535533905932737*alphaDrSurf_l[2]*fUpwind_l[2]+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[1]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[0]; 
  Gdrag_l[1] = 0.3535533905932737*alphaDrSurf_l[2]*fUpwind_l[4]+0.3535533905932737*fUpwind_l[2]*alphaDrSurf_l[4]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[1]+0.3535533905932737*fUpwind_l[0]*alphaDrSurf_l[1]; 
  Gdrag_l[2] = 0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[4]+0.3535533905932737*fUpwind_l[1]*alphaDrSurf_l[4]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[2]+0.3535533905932737*fUpwind_l[0]*alphaDrSurf_l[2]; 
  Gdrag_l[3] = 0.3535533905932737*alphaDrSurf_l[4]*fUpwind_l[7]+0.3535533905932737*alphaDrSurf_l[2]*fUpwind_l[6]+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[5]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[3]; 
  Gdrag_l[4] = 0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[4]+0.3535533905932737*fUpwind_l[0]*alphaDrSurf_l[4]+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[2]+0.3535533905932737*fUpwind_l[1]*alphaDrSurf_l[2]; 
  Gdrag_l[5] = 0.3535533905932737*alphaDrSurf_l[2]*fUpwind_l[7]+0.3535533905932737*alphaDrSurf_l[4]*fUpwind_l[6]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[5]+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[3]; 
  Gdrag_l[6] = 0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[7]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[6]+0.3535533905932737*alphaDrSurf_l[4]*fUpwind_l[5]+0.3535533905932737*alphaDrSurf_l[2]*fUpwind_l[3]; 
  Gdrag_l[7] = 0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[7]+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[6]+0.3535533905932737*alphaDrSurf_l[2]*fUpwind_l[5]+0.3535533905932737*fUpwind_l[3]*alphaDrSurf_l[4]; 

  Gdrag_r[0] = 0.3535533905932737*alphaDrSurf_r[4]*fUpwind_r[4]+0.3535533905932737*alphaDrSurf_r[2]*fUpwind_r[2]+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[1]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[0]; 
  Gdrag_r[1] = 0.3535533905932737*alphaDrSurf_r[2]*fUpwind_r[4]+0.3535533905932737*fUpwind_r[2]*alphaDrSurf_r[4]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[1]+0.3535533905932737*fUpwind_r[0]*alphaDrSurf_r[1]; 
  Gdrag_r[2] = 0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[4]+0.3535533905932737*fUpwind_r[1]*alphaDrSurf_r[4]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[2]+0.3535533905932737*fUpwind_r[0]*alphaDrSurf_r[2]; 
  Gdrag_r[3] = 0.3535533905932737*alphaDrSurf_r[4]*fUpwind_r[7]+0.3535533905932737*alphaDrSurf_r[2]*fUpwind_r[6]+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[5]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[3]; 
  Gdrag_r[4] = 0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[4]+0.3535533905932737*fUpwind_r[0]*alphaDrSurf_r[4]+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[2]+0.3535533905932737*fUpwind_r[1]*alphaDrSurf_r[2]; 
  Gdrag_r[5] = 0.3535533905932737*alphaDrSurf_r[2]*fUpwind_r[7]+0.3535533905932737*alphaDrSurf_r[4]*fUpwind_r[6]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[5]+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[3]; 
  Gdrag_r[6] = 0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[7]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[6]+0.3535533905932737*alphaDrSurf_r[4]*fUpwind_r[5]+0.3535533905932737*alphaDrSurf_r[2]*fUpwind_r[3]; 
  Gdrag_r[7] = 0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[7]+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[6]+0.3535533905932737*alphaDrSurf_r[2]*fUpwind_r[5]+0.3535533905932737*fUpwind_r[3]*alphaDrSurf_r[4]; 

  out[0] += 0.7071067811865475*Gdrag_r[0]*rdv2-0.7071067811865475*Gdrag_l[0]*rdv2; 
  out[1] += 0.7071067811865475*Gdrag_r[1]*rdv2-0.7071067811865475*Gdrag_l[1]*rdv2; 
  out[2] += 0.7071067811865475*Gdrag_r[2]*rdv2-0.7071067811865475*Gdrag_l[2]*rdv2; 
  out[3] += 1.224744871391589*Gdrag_r[0]*rdv2+1.224744871391589*Gdrag_l[0]*rdv2; 
  out[4] += 0.7071067811865475*Gdrag_r[3]*rdv2-0.7071067811865475*Gdrag_l[3]*rdv2; 
  out[5] += 0.7071067811865475*Gdrag_r[4]*rdv2-0.7071067811865475*Gdrag_l[4]*rdv2; 
  out[6] += 1.224744871391589*Gdrag_r[1]*rdv2+1.224744871391589*Gdrag_l[1]*rdv2; 
  out[7] += 1.224744871391589*Gdrag_r[2]*rdv2+1.224744871391589*Gdrag_l[2]*rdv2; 
  out[8] += 0.7071067811865475*Gdrag_r[5]*rdv2-0.7071067811865475*Gdrag_l[5]*rdv2; 
  out[9] += 0.7071067811865475*Gdrag_r[6]*rdv2-0.7071067811865475*Gdrag_l[6]*rdv2; 
  out[10] += 1.224744871391589*Gdrag_r[3]*rdv2+1.224744871391589*Gdrag_l[3]*rdv2; 
  out[11] += 1.224744871391589*Gdrag_r[4]*rdv2+1.224744871391589*Gdrag_l[4]*rdv2; 
  out[12] += 0.7071067811865475*Gdrag_r[7]*rdv2-0.7071067811865475*Gdrag_l[7]*rdv2; 
  out[13] += 1.224744871391589*Gdrag_r[5]*rdv2+1.224744871391589*Gdrag_l[5]*rdv2; 
  out[14] += 1.224744871391589*Gdrag_r[6]*rdv2+1.224744871391589*Gdrag_l[6]*rdv2; 
  out[15] += 1.224744871391589*Gdrag_r[7]*rdv2+1.224744871391589*Gdrag_l[7]*rdv2; 
} 
