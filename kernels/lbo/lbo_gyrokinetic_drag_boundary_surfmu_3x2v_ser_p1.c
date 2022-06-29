#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH void lbo_gyrokinetic_drag_boundary_surfmu_3x2v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[5]:     cell-center coordinates. 
  // dxv[5]:   cell spacing. 
  // m_:        species mass.
  // bmag_inv:  1/(magnetic field magnitude). 
  // nuSum:     collisionalities added (self and cross species collisionalities). 
  // nuUSum[16]:sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum: sum of thermal speeds squared time their respective collisionalities. 
  // fSkin/Edge:    Distribution function in cells 
  // out:           Incremented distribution function in cell 
  double rdv2 = 2.0/dxv[4]; 

  double alphaDrSurf[16] = {0.0}; 
  double fUpwindQuad[16] = {0.0};
  double fUpwind[16] = {0.0};;
  double Ghat[16] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = nuSum[0]*(2.828427124746191*w[4]+1.414213562373095*dxv[4]); 
  alphaDrSurf[1] = nuSum[1]*(2.828427124746191*w[4]+1.414213562373095*dxv[4]); 
  alphaDrSurf[2] = nuSum[2]*(2.828427124746191*w[4]+1.414213562373095*dxv[4]); 
  alphaDrSurf[3] = nuSum[3]*(2.828427124746191*w[4]+1.414213562373095*dxv[4]); 
  alphaDrSurf[5] = nuSum[4]*(2.828427124746191*w[4]+1.414213562373095*dxv[4]); 
  alphaDrSurf[6] = (2.828427124746191*w[4]+1.414213562373095*dxv[4])*nuSum[5]; 
  alphaDrSurf[7] = (2.828427124746191*w[4]+1.414213562373095*dxv[4])*nuSum[6]; 
  alphaDrSurf[11] = (2.828427124746191*w[4]+1.414213562373095*dxv[4])*nuSum[7]; 

  Ghat[0] = -0.125*(2.449489742783178*(alphaDrSurf[11]*fEdge[27]+alphaDrSurf[7]*fEdge[22]+alphaDrSurf[6]*fEdge[21]+alphaDrSurf[5]*fEdge[20])-1.414213562373095*alphaDrSurf[11]*fEdge[16]+2.449489742783178*(alphaDrSurf[3]*fEdge[14]+alphaDrSurf[2]*fEdge[13]+alphaDrSurf[1]*fEdge[12])-1.414213562373095*(alphaDrSurf[7]*fEdge[8]+alphaDrSurf[6]*fEdge[7]+alphaDrSurf[5]*fEdge[6])+2.449489742783178*alphaDrSurf[0]*fEdge[5]-1.414213562373095*(alphaDrSurf[3]*fEdge[3]+alphaDrSurf[2]*fEdge[2]+alphaDrSurf[1]*fEdge[1]+alphaDrSurf[0]*fEdge[0])); 
  Ghat[1] = -0.125*(2.449489742783178*(alphaDrSurf[7]*fEdge[27]+alphaDrSurf[11]*fEdge[22]+alphaDrSurf[3]*fEdge[21]+alphaDrSurf[2]*fEdge[20])-1.414213562373095*alphaDrSurf[7]*fEdge[16]+2.449489742783178*(alphaDrSurf[6]*fEdge[14]+alphaDrSurf[5]*fEdge[13]+alphaDrSurf[0]*fEdge[12])-1.414213562373095*(fEdge[8]*alphaDrSurf[11]+alphaDrSurf[3]*fEdge[7]+alphaDrSurf[2]*fEdge[6]+fEdge[3]*alphaDrSurf[6])+2.449489742783178*alphaDrSurf[1]*fEdge[5]-1.414213562373095*(fEdge[2]*alphaDrSurf[5]+alphaDrSurf[0]*fEdge[1]+fEdge[0]*alphaDrSurf[1])); 
  Ghat[2] = -0.125*(2.449489742783178*(alphaDrSurf[6]*fEdge[27]+alphaDrSurf[3]*fEdge[22]+alphaDrSurf[11]*fEdge[21]+alphaDrSurf[1]*fEdge[20])-1.414213562373095*alphaDrSurf[6]*fEdge[16]+2.449489742783178*(alphaDrSurf[7]*fEdge[14]+alphaDrSurf[0]*fEdge[13]+alphaDrSurf[5]*fEdge[12])-1.414213562373095*(fEdge[7]*alphaDrSurf[11]+alphaDrSurf[3]*fEdge[8]+fEdge[3]*alphaDrSurf[7]+alphaDrSurf[1]*fEdge[6])+2.449489742783178*alphaDrSurf[2]*fEdge[5]-1.414213562373095*(fEdge[1]*alphaDrSurf[5]+alphaDrSurf[0]*fEdge[2]+fEdge[0]*alphaDrSurf[2])); 
  Ghat[3] = -0.125*(2.449489742783178*(alphaDrSurf[5]*fEdge[27]+alphaDrSurf[2]*fEdge[22]+alphaDrSurf[1]*fEdge[21]+alphaDrSurf[11]*fEdge[20])-1.414213562373095*alphaDrSurf[5]*fEdge[16]+2.449489742783178*(alphaDrSurf[0]*fEdge[14]+alphaDrSurf[7]*fEdge[13]+alphaDrSurf[6]*fEdge[12])-1.414213562373095*(fEdge[6]*alphaDrSurf[11]+alphaDrSurf[2]*fEdge[8]+alphaDrSurf[1]*fEdge[7]+fEdge[2]*alphaDrSurf[7]+fEdge[1]*alphaDrSurf[6])+2.449489742783178*alphaDrSurf[3]*fEdge[5]-1.414213562373095*(alphaDrSurf[0]*fEdge[3]+fEdge[0]*alphaDrSurf[3])); 
  Ghat[4] = -0.125*(2.449489742783178*(alphaDrSurf[11]*fEdge[31]+alphaDrSurf[7]*fEdge[30]+alphaDrSurf[6]*fEdge[29]+alphaDrSurf[5]*fEdge[28])-1.414213562373095*alphaDrSurf[11]*fEdge[26]+2.449489742783178*(alphaDrSurf[3]*fEdge[25]+alphaDrSurf[2]*fEdge[24]+alphaDrSurf[1]*fEdge[23])-1.414213562373095*(alphaDrSurf[7]*fEdge[19]+alphaDrSurf[6]*fEdge[18]+alphaDrSurf[5]*fEdge[17])+2.449489742783178*alphaDrSurf[0]*fEdge[15]-1.414213562373095*(alphaDrSurf[3]*fEdge[11]+alphaDrSurf[2]*fEdge[10]+alphaDrSurf[1]*fEdge[9]+alphaDrSurf[0]*fEdge[4])); 
  Ghat[5] = -0.125*(2.449489742783178*(alphaDrSurf[3]*fEdge[27]+alphaDrSurf[6]*fEdge[22]+alphaDrSurf[7]*fEdge[21]+alphaDrSurf[0]*fEdge[20])-1.414213562373095*alphaDrSurf[3]*fEdge[16]+2.449489742783178*(alphaDrSurf[11]*fEdge[14]+alphaDrSurf[1]*fEdge[13]+alphaDrSurf[2]*fEdge[12])-1.414213562373095*(fEdge[3]*alphaDrSurf[11]+alphaDrSurf[6]*fEdge[8]+alphaDrSurf[7]*fEdge[7]+alphaDrSurf[0]*fEdge[6])+2.449489742783178*alphaDrSurf[5]*fEdge[5]-1.414213562373095*(fEdge[0]*alphaDrSurf[5]+alphaDrSurf[1]*fEdge[2]+fEdge[1]*alphaDrSurf[2])); 
  Ghat[6] = -0.125*(2.449489742783178*(alphaDrSurf[2]*fEdge[27]+alphaDrSurf[5]*fEdge[22]+alphaDrSurf[0]*fEdge[21]+alphaDrSurf[7]*fEdge[20])-1.414213562373095*alphaDrSurf[2]*fEdge[16]+2.449489742783178*(alphaDrSurf[1]*fEdge[14]+alphaDrSurf[11]*fEdge[13]+alphaDrSurf[3]*fEdge[12])-1.414213562373095*(fEdge[2]*alphaDrSurf[11]+alphaDrSurf[5]*fEdge[8]+alphaDrSurf[0]*fEdge[7]+fEdge[6]*alphaDrSurf[7])+(2.449489742783178*fEdge[5]-1.414213562373095*fEdge[0])*alphaDrSurf[6]-1.414213562373095*(alphaDrSurf[1]*fEdge[3]+fEdge[1]*alphaDrSurf[3])); 
  Ghat[7] = -0.125*(2.449489742783178*(alphaDrSurf[1]*fEdge[27]+alphaDrSurf[0]*fEdge[22]+alphaDrSurf[5]*fEdge[21]+alphaDrSurf[6]*fEdge[20])-1.414213562373095*alphaDrSurf[1]*fEdge[16]+2.449489742783178*(alphaDrSurf[2]*fEdge[14]+alphaDrSurf[3]*fEdge[13]+alphaDrSurf[11]*fEdge[12])-1.414213562373095*(fEdge[1]*alphaDrSurf[11]+alphaDrSurf[0]*fEdge[8]+alphaDrSurf[5]*fEdge[7])+(2.449489742783178*fEdge[5]-1.414213562373095*fEdge[0])*alphaDrSurf[7]-1.414213562373095*(alphaDrSurf[6]*fEdge[6]+alphaDrSurf[2]*fEdge[3]+fEdge[2]*alphaDrSurf[3])); 
  Ghat[8] = -0.125*(2.449489742783178*(alphaDrSurf[7]*fEdge[31]+alphaDrSurf[11]*fEdge[30]+alphaDrSurf[3]*fEdge[29]+alphaDrSurf[2]*fEdge[28])-1.414213562373095*alphaDrSurf[7]*fEdge[26]+2.449489742783178*(alphaDrSurf[6]*fEdge[25]+alphaDrSurf[5]*fEdge[24]+alphaDrSurf[0]*fEdge[23])-1.414213562373095*(alphaDrSurf[11]*fEdge[19]+alphaDrSurf[3]*fEdge[18]+alphaDrSurf[2]*fEdge[17])+2.449489742783178*alphaDrSurf[1]*fEdge[15]-1.414213562373095*(alphaDrSurf[6]*fEdge[11]+alphaDrSurf[5]*fEdge[10]+alphaDrSurf[0]*fEdge[9]+alphaDrSurf[1]*fEdge[4])); 
  Ghat[9] = -0.125*(2.449489742783178*(alphaDrSurf[6]*fEdge[31]+alphaDrSurf[3]*fEdge[30]+alphaDrSurf[11]*fEdge[29]+alphaDrSurf[1]*fEdge[28])-1.414213562373095*alphaDrSurf[6]*fEdge[26]+2.449489742783178*(alphaDrSurf[7]*fEdge[25]+alphaDrSurf[0]*fEdge[24]+alphaDrSurf[5]*fEdge[23])-1.414213562373095*(alphaDrSurf[3]*fEdge[19]+alphaDrSurf[11]*fEdge[18]+alphaDrSurf[1]*fEdge[17])+2.449489742783178*alphaDrSurf[2]*fEdge[15]-1.414213562373095*(alphaDrSurf[7]*fEdge[11]+alphaDrSurf[0]*fEdge[10]+alphaDrSurf[5]*fEdge[9]+alphaDrSurf[2]*fEdge[4])); 
  Ghat[10] = -0.125*(2.449489742783178*(alphaDrSurf[5]*fEdge[31]+alphaDrSurf[2]*fEdge[30]+alphaDrSurf[1]*fEdge[29]+alphaDrSurf[11]*fEdge[28])-1.414213562373095*alphaDrSurf[5]*fEdge[26]+2.449489742783178*(alphaDrSurf[0]*fEdge[25]+alphaDrSurf[7]*fEdge[24]+alphaDrSurf[6]*fEdge[23])-1.414213562373095*(alphaDrSurf[2]*fEdge[19]+alphaDrSurf[1]*fEdge[18]+alphaDrSurf[11]*fEdge[17])+2.449489742783178*alphaDrSurf[3]*fEdge[15]-1.414213562373095*(alphaDrSurf[0]*fEdge[11]+alphaDrSurf[7]*fEdge[10]+alphaDrSurf[6]*fEdge[9]+alphaDrSurf[3]*fEdge[4])); 
  Ghat[11] = -0.125*(2.449489742783178*(alphaDrSurf[0]*fEdge[27]+alphaDrSurf[1]*fEdge[22]+alphaDrSurf[2]*fEdge[21]+alphaDrSurf[3]*fEdge[20])-1.414213562373095*alphaDrSurf[0]*fEdge[16]+2.449489742783178*(alphaDrSurf[5]*fEdge[14]+alphaDrSurf[6]*fEdge[13]+alphaDrSurf[7]*fEdge[12])+(2.449489742783178*fEdge[5]-1.414213562373095*fEdge[0])*alphaDrSurf[11]-1.414213562373095*(alphaDrSurf[1]*fEdge[8]+alphaDrSurf[2]*fEdge[7]+fEdge[1]*alphaDrSurf[7]+alphaDrSurf[3]*fEdge[6]+fEdge[2]*alphaDrSurf[6]+fEdge[3]*alphaDrSurf[5])); 
  Ghat[12] = -0.125*(2.449489742783178*(alphaDrSurf[3]*fEdge[31]+alphaDrSurf[6]*fEdge[30]+alphaDrSurf[7]*fEdge[29]+alphaDrSurf[0]*fEdge[28])-1.414213562373095*alphaDrSurf[3]*fEdge[26]+2.449489742783178*(alphaDrSurf[11]*fEdge[25]+alphaDrSurf[1]*fEdge[24]+alphaDrSurf[2]*fEdge[23])-1.414213562373095*(alphaDrSurf[6]*fEdge[19]+alphaDrSurf[7]*fEdge[18]+alphaDrSurf[0]*fEdge[17])+2.449489742783178*alphaDrSurf[5]*fEdge[15]-1.414213562373095*(alphaDrSurf[11]*fEdge[11]+alphaDrSurf[1]*fEdge[10]+alphaDrSurf[2]*fEdge[9]+fEdge[4]*alphaDrSurf[5])); 
  Ghat[13] = -0.125*(2.449489742783178*(alphaDrSurf[2]*fEdge[31]+alphaDrSurf[5]*fEdge[30]+alphaDrSurf[0]*fEdge[29]+alphaDrSurf[7]*fEdge[28])-1.414213562373095*alphaDrSurf[2]*fEdge[26]+2.449489742783178*(alphaDrSurf[1]*fEdge[25]+alphaDrSurf[11]*fEdge[24]+alphaDrSurf[3]*fEdge[23])-1.414213562373095*(alphaDrSurf[5]*fEdge[19]+alphaDrSurf[0]*fEdge[18]+alphaDrSurf[7]*fEdge[17])+2.449489742783178*alphaDrSurf[6]*fEdge[15]-1.414213562373095*(alphaDrSurf[1]*fEdge[11]+fEdge[10]*alphaDrSurf[11]+alphaDrSurf[3]*fEdge[9]+fEdge[4]*alphaDrSurf[6])); 
  Ghat[14] = -0.125*(2.449489742783178*(alphaDrSurf[1]*fEdge[31]+alphaDrSurf[0]*fEdge[30]+alphaDrSurf[5]*fEdge[29]+alphaDrSurf[6]*fEdge[28])-1.414213562373095*alphaDrSurf[1]*fEdge[26]+2.449489742783178*(alphaDrSurf[2]*fEdge[25]+alphaDrSurf[3]*fEdge[24]+alphaDrSurf[11]*fEdge[23])-1.414213562373095*(alphaDrSurf[0]*fEdge[19]+alphaDrSurf[5]*fEdge[18]+alphaDrSurf[6]*fEdge[17])+2.449489742783178*alphaDrSurf[7]*fEdge[15]-1.414213562373095*(alphaDrSurf[2]*fEdge[11]+fEdge[9]*alphaDrSurf[11]+alphaDrSurf[3]*fEdge[10]+fEdge[4]*alphaDrSurf[7])); 
  Ghat[15] = -0.125*(2.449489742783178*(alphaDrSurf[0]*fEdge[31]+alphaDrSurf[1]*fEdge[30]+alphaDrSurf[2]*fEdge[29]+alphaDrSurf[3]*fEdge[28])-1.414213562373095*alphaDrSurf[0]*fEdge[26]+2.449489742783178*(alphaDrSurf[5]*fEdge[25]+alphaDrSurf[6]*fEdge[24]+alphaDrSurf[7]*fEdge[23])-1.414213562373095*(alphaDrSurf[1]*fEdge[19]+alphaDrSurf[2]*fEdge[18]+alphaDrSurf[3]*fEdge[17])+2.449489742783178*alphaDrSurf[11]*fEdge[15]-1.414213562373095*(alphaDrSurf[5]*fEdge[11]+fEdge[4]*alphaDrSurf[11]+alphaDrSurf[6]*fEdge[10]+alphaDrSurf[7]*fEdge[9])); 

  out[0] += 0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += 0.7071067811865475*Ghat[2]*rdv2; 
  out[3] += 0.7071067811865475*Ghat[3]*rdv2; 
  out[4] += 0.7071067811865475*Ghat[4]*rdv2; 
  out[5] += 1.224744871391589*Ghat[0]*rdv2; 
  out[6] += 0.7071067811865475*Ghat[5]*rdv2; 
  out[7] += 0.7071067811865475*Ghat[6]*rdv2; 
  out[8] += 0.7071067811865475*Ghat[7]*rdv2; 
  out[9] += 0.7071067811865475*Ghat[8]*rdv2; 
  out[10] += 0.7071067811865475*Ghat[9]*rdv2; 
  out[11] += 0.7071067811865475*Ghat[10]*rdv2; 
  out[12] += 1.224744871391589*Ghat[1]*rdv2; 
  out[13] += 1.224744871391589*Ghat[2]*rdv2; 
  out[14] += 1.224744871391589*Ghat[3]*rdv2; 
  out[15] += 1.224744871391589*Ghat[4]*rdv2; 
  out[16] += 0.7071067811865475*Ghat[11]*rdv2; 
  out[17] += 0.7071067811865475*Ghat[12]*rdv2; 
  out[18] += 0.7071067811865475*Ghat[13]*rdv2; 
  out[19] += 0.7071067811865475*Ghat[14]*rdv2; 
  out[20] += 1.224744871391589*Ghat[5]*rdv2; 
  out[21] += 1.224744871391589*Ghat[6]*rdv2; 
  out[22] += 1.224744871391589*Ghat[7]*rdv2; 
  out[23] += 1.224744871391589*Ghat[8]*rdv2; 
  out[24] += 1.224744871391589*Ghat[9]*rdv2; 
  out[25] += 1.224744871391589*Ghat[10]*rdv2; 
  out[26] += 0.7071067811865475*Ghat[15]*rdv2; 
  out[27] += 1.224744871391589*Ghat[11]*rdv2; 
  out[28] += 1.224744871391589*Ghat[12]*rdv2; 
  out[29] += 1.224744871391589*Ghat[13]*rdv2; 
  out[30] += 1.224744871391589*Ghat[14]*rdv2; 
  out[31] += 1.224744871391589*Ghat[15]*rdv2; 

  } else { 

  alphaDrSurf[0] = nuSum[0]*(2.828427124746191*w[4]-1.414213562373095*dxv[4]); 
  alphaDrSurf[1] = nuSum[1]*(2.828427124746191*w[4]-1.414213562373095*dxv[4]); 
  alphaDrSurf[2] = nuSum[2]*(2.828427124746191*w[4]-1.414213562373095*dxv[4]); 
  alphaDrSurf[3] = nuSum[3]*(2.828427124746191*w[4]-1.414213562373095*dxv[4]); 
  alphaDrSurf[5] = nuSum[4]*(2.828427124746191*w[4]-1.414213562373095*dxv[4]); 
  alphaDrSurf[6] = (2.828427124746191*w[4]-1.414213562373095*dxv[4])*nuSum[5]; 
  alphaDrSurf[7] = (2.828427124746191*w[4]-1.414213562373095*dxv[4])*nuSum[6]; 
  alphaDrSurf[11] = (2.828427124746191*w[4]-1.414213562373095*dxv[4])*nuSum[7]; 

  Ghat[0] = -0.125*(2.449489742783178*(alphaDrSurf[11]*fSkin[27]+alphaDrSurf[7]*fSkin[22]+alphaDrSurf[6]*fSkin[21]+alphaDrSurf[5]*fSkin[20])-1.414213562373095*alphaDrSurf[11]*fSkin[16]+2.449489742783178*(alphaDrSurf[3]*fSkin[14]+alphaDrSurf[2]*fSkin[13]+alphaDrSurf[1]*fSkin[12])-1.414213562373095*(alphaDrSurf[7]*fSkin[8]+alphaDrSurf[6]*fSkin[7]+alphaDrSurf[5]*fSkin[6])+2.449489742783178*alphaDrSurf[0]*fSkin[5]-1.414213562373095*(alphaDrSurf[3]*fSkin[3]+alphaDrSurf[2]*fSkin[2]+alphaDrSurf[1]*fSkin[1]+alphaDrSurf[0]*fSkin[0])); 
  Ghat[1] = -0.125*(2.449489742783178*(alphaDrSurf[7]*fSkin[27]+alphaDrSurf[11]*fSkin[22]+alphaDrSurf[3]*fSkin[21]+alphaDrSurf[2]*fSkin[20])-1.414213562373095*alphaDrSurf[7]*fSkin[16]+2.449489742783178*(alphaDrSurf[6]*fSkin[14]+alphaDrSurf[5]*fSkin[13]+alphaDrSurf[0]*fSkin[12])-1.414213562373095*(fSkin[8]*alphaDrSurf[11]+alphaDrSurf[3]*fSkin[7]+alphaDrSurf[2]*fSkin[6]+fSkin[3]*alphaDrSurf[6])+2.449489742783178*alphaDrSurf[1]*fSkin[5]-1.414213562373095*(fSkin[2]*alphaDrSurf[5]+alphaDrSurf[0]*fSkin[1]+fSkin[0]*alphaDrSurf[1])); 
  Ghat[2] = -0.125*(2.449489742783178*(alphaDrSurf[6]*fSkin[27]+alphaDrSurf[3]*fSkin[22]+alphaDrSurf[11]*fSkin[21]+alphaDrSurf[1]*fSkin[20])-1.414213562373095*alphaDrSurf[6]*fSkin[16]+2.449489742783178*(alphaDrSurf[7]*fSkin[14]+alphaDrSurf[0]*fSkin[13]+alphaDrSurf[5]*fSkin[12])-1.414213562373095*(fSkin[7]*alphaDrSurf[11]+alphaDrSurf[3]*fSkin[8]+fSkin[3]*alphaDrSurf[7]+alphaDrSurf[1]*fSkin[6])+2.449489742783178*alphaDrSurf[2]*fSkin[5]-1.414213562373095*(fSkin[1]*alphaDrSurf[5]+alphaDrSurf[0]*fSkin[2]+fSkin[0]*alphaDrSurf[2])); 
  Ghat[3] = -0.125*(2.449489742783178*(alphaDrSurf[5]*fSkin[27]+alphaDrSurf[2]*fSkin[22]+alphaDrSurf[1]*fSkin[21]+alphaDrSurf[11]*fSkin[20])-1.414213562373095*alphaDrSurf[5]*fSkin[16]+2.449489742783178*(alphaDrSurf[0]*fSkin[14]+alphaDrSurf[7]*fSkin[13]+alphaDrSurf[6]*fSkin[12])-1.414213562373095*(fSkin[6]*alphaDrSurf[11]+alphaDrSurf[2]*fSkin[8]+alphaDrSurf[1]*fSkin[7]+fSkin[2]*alphaDrSurf[7]+fSkin[1]*alphaDrSurf[6])+2.449489742783178*alphaDrSurf[3]*fSkin[5]-1.414213562373095*(alphaDrSurf[0]*fSkin[3]+fSkin[0]*alphaDrSurf[3])); 
  Ghat[4] = -0.125*(2.449489742783178*(alphaDrSurf[11]*fSkin[31]+alphaDrSurf[7]*fSkin[30]+alphaDrSurf[6]*fSkin[29]+alphaDrSurf[5]*fSkin[28])-1.414213562373095*alphaDrSurf[11]*fSkin[26]+2.449489742783178*(alphaDrSurf[3]*fSkin[25]+alphaDrSurf[2]*fSkin[24]+alphaDrSurf[1]*fSkin[23])-1.414213562373095*(alphaDrSurf[7]*fSkin[19]+alphaDrSurf[6]*fSkin[18]+alphaDrSurf[5]*fSkin[17])+2.449489742783178*alphaDrSurf[0]*fSkin[15]-1.414213562373095*(alphaDrSurf[3]*fSkin[11]+alphaDrSurf[2]*fSkin[10]+alphaDrSurf[1]*fSkin[9]+alphaDrSurf[0]*fSkin[4])); 
  Ghat[5] = -0.125*(2.449489742783178*(alphaDrSurf[3]*fSkin[27]+alphaDrSurf[6]*fSkin[22]+alphaDrSurf[7]*fSkin[21]+alphaDrSurf[0]*fSkin[20])-1.414213562373095*alphaDrSurf[3]*fSkin[16]+2.449489742783178*(alphaDrSurf[11]*fSkin[14]+alphaDrSurf[1]*fSkin[13]+alphaDrSurf[2]*fSkin[12])-1.414213562373095*(fSkin[3]*alphaDrSurf[11]+alphaDrSurf[6]*fSkin[8]+alphaDrSurf[7]*fSkin[7]+alphaDrSurf[0]*fSkin[6])+2.449489742783178*alphaDrSurf[5]*fSkin[5]-1.414213562373095*(fSkin[0]*alphaDrSurf[5]+alphaDrSurf[1]*fSkin[2]+fSkin[1]*alphaDrSurf[2])); 
  Ghat[6] = -0.125*(2.449489742783178*(alphaDrSurf[2]*fSkin[27]+alphaDrSurf[5]*fSkin[22]+alphaDrSurf[0]*fSkin[21]+alphaDrSurf[7]*fSkin[20])-1.414213562373095*alphaDrSurf[2]*fSkin[16]+2.449489742783178*(alphaDrSurf[1]*fSkin[14]+alphaDrSurf[11]*fSkin[13]+alphaDrSurf[3]*fSkin[12])-1.414213562373095*(fSkin[2]*alphaDrSurf[11]+alphaDrSurf[5]*fSkin[8]+alphaDrSurf[0]*fSkin[7]+fSkin[6]*alphaDrSurf[7])+(2.449489742783178*fSkin[5]-1.414213562373095*fSkin[0])*alphaDrSurf[6]-1.414213562373095*(alphaDrSurf[1]*fSkin[3]+fSkin[1]*alphaDrSurf[3])); 
  Ghat[7] = -0.125*(2.449489742783178*(alphaDrSurf[1]*fSkin[27]+alphaDrSurf[0]*fSkin[22]+alphaDrSurf[5]*fSkin[21]+alphaDrSurf[6]*fSkin[20])-1.414213562373095*alphaDrSurf[1]*fSkin[16]+2.449489742783178*(alphaDrSurf[2]*fSkin[14]+alphaDrSurf[3]*fSkin[13]+alphaDrSurf[11]*fSkin[12])-1.414213562373095*(fSkin[1]*alphaDrSurf[11]+alphaDrSurf[0]*fSkin[8]+alphaDrSurf[5]*fSkin[7])+(2.449489742783178*fSkin[5]-1.414213562373095*fSkin[0])*alphaDrSurf[7]-1.414213562373095*(alphaDrSurf[6]*fSkin[6]+alphaDrSurf[2]*fSkin[3]+fSkin[2]*alphaDrSurf[3])); 
  Ghat[8] = -0.125*(2.449489742783178*(alphaDrSurf[7]*fSkin[31]+alphaDrSurf[11]*fSkin[30]+alphaDrSurf[3]*fSkin[29]+alphaDrSurf[2]*fSkin[28])-1.414213562373095*alphaDrSurf[7]*fSkin[26]+2.449489742783178*(alphaDrSurf[6]*fSkin[25]+alphaDrSurf[5]*fSkin[24]+alphaDrSurf[0]*fSkin[23])-1.414213562373095*(alphaDrSurf[11]*fSkin[19]+alphaDrSurf[3]*fSkin[18]+alphaDrSurf[2]*fSkin[17])+2.449489742783178*alphaDrSurf[1]*fSkin[15]-1.414213562373095*(alphaDrSurf[6]*fSkin[11]+alphaDrSurf[5]*fSkin[10]+alphaDrSurf[0]*fSkin[9]+alphaDrSurf[1]*fSkin[4])); 
  Ghat[9] = -0.125*(2.449489742783178*(alphaDrSurf[6]*fSkin[31]+alphaDrSurf[3]*fSkin[30]+alphaDrSurf[11]*fSkin[29]+alphaDrSurf[1]*fSkin[28])-1.414213562373095*alphaDrSurf[6]*fSkin[26]+2.449489742783178*(alphaDrSurf[7]*fSkin[25]+alphaDrSurf[0]*fSkin[24]+alphaDrSurf[5]*fSkin[23])-1.414213562373095*(alphaDrSurf[3]*fSkin[19]+alphaDrSurf[11]*fSkin[18]+alphaDrSurf[1]*fSkin[17])+2.449489742783178*alphaDrSurf[2]*fSkin[15]-1.414213562373095*(alphaDrSurf[7]*fSkin[11]+alphaDrSurf[0]*fSkin[10]+alphaDrSurf[5]*fSkin[9]+alphaDrSurf[2]*fSkin[4])); 
  Ghat[10] = -0.125*(2.449489742783178*(alphaDrSurf[5]*fSkin[31]+alphaDrSurf[2]*fSkin[30]+alphaDrSurf[1]*fSkin[29]+alphaDrSurf[11]*fSkin[28])-1.414213562373095*alphaDrSurf[5]*fSkin[26]+2.449489742783178*(alphaDrSurf[0]*fSkin[25]+alphaDrSurf[7]*fSkin[24]+alphaDrSurf[6]*fSkin[23])-1.414213562373095*(alphaDrSurf[2]*fSkin[19]+alphaDrSurf[1]*fSkin[18]+alphaDrSurf[11]*fSkin[17])+2.449489742783178*alphaDrSurf[3]*fSkin[15]-1.414213562373095*(alphaDrSurf[0]*fSkin[11]+alphaDrSurf[7]*fSkin[10]+alphaDrSurf[6]*fSkin[9]+alphaDrSurf[3]*fSkin[4])); 
  Ghat[11] = -0.125*(2.449489742783178*(alphaDrSurf[0]*fSkin[27]+alphaDrSurf[1]*fSkin[22]+alphaDrSurf[2]*fSkin[21]+alphaDrSurf[3]*fSkin[20])-1.414213562373095*alphaDrSurf[0]*fSkin[16]+2.449489742783178*(alphaDrSurf[5]*fSkin[14]+alphaDrSurf[6]*fSkin[13]+alphaDrSurf[7]*fSkin[12])+(2.449489742783178*fSkin[5]-1.414213562373095*fSkin[0])*alphaDrSurf[11]-1.414213562373095*(alphaDrSurf[1]*fSkin[8]+alphaDrSurf[2]*fSkin[7]+fSkin[1]*alphaDrSurf[7]+alphaDrSurf[3]*fSkin[6]+fSkin[2]*alphaDrSurf[6]+fSkin[3]*alphaDrSurf[5])); 
  Ghat[12] = -0.125*(2.449489742783178*(alphaDrSurf[3]*fSkin[31]+alphaDrSurf[6]*fSkin[30]+alphaDrSurf[7]*fSkin[29]+alphaDrSurf[0]*fSkin[28])-1.414213562373095*alphaDrSurf[3]*fSkin[26]+2.449489742783178*(alphaDrSurf[11]*fSkin[25]+alphaDrSurf[1]*fSkin[24]+alphaDrSurf[2]*fSkin[23])-1.414213562373095*(alphaDrSurf[6]*fSkin[19]+alphaDrSurf[7]*fSkin[18]+alphaDrSurf[0]*fSkin[17])+2.449489742783178*alphaDrSurf[5]*fSkin[15]-1.414213562373095*(alphaDrSurf[11]*fSkin[11]+alphaDrSurf[1]*fSkin[10]+alphaDrSurf[2]*fSkin[9]+fSkin[4]*alphaDrSurf[5])); 
  Ghat[13] = -0.125*(2.449489742783178*(alphaDrSurf[2]*fSkin[31]+alphaDrSurf[5]*fSkin[30]+alphaDrSurf[0]*fSkin[29]+alphaDrSurf[7]*fSkin[28])-1.414213562373095*alphaDrSurf[2]*fSkin[26]+2.449489742783178*(alphaDrSurf[1]*fSkin[25]+alphaDrSurf[11]*fSkin[24]+alphaDrSurf[3]*fSkin[23])-1.414213562373095*(alphaDrSurf[5]*fSkin[19]+alphaDrSurf[0]*fSkin[18]+alphaDrSurf[7]*fSkin[17])+2.449489742783178*alphaDrSurf[6]*fSkin[15]-1.414213562373095*(alphaDrSurf[1]*fSkin[11]+fSkin[10]*alphaDrSurf[11]+alphaDrSurf[3]*fSkin[9]+fSkin[4]*alphaDrSurf[6])); 
  Ghat[14] = -0.125*(2.449489742783178*(alphaDrSurf[1]*fSkin[31]+alphaDrSurf[0]*fSkin[30]+alphaDrSurf[5]*fSkin[29]+alphaDrSurf[6]*fSkin[28])-1.414213562373095*alphaDrSurf[1]*fSkin[26]+2.449489742783178*(alphaDrSurf[2]*fSkin[25]+alphaDrSurf[3]*fSkin[24]+alphaDrSurf[11]*fSkin[23])-1.414213562373095*(alphaDrSurf[0]*fSkin[19]+alphaDrSurf[5]*fSkin[18]+alphaDrSurf[6]*fSkin[17])+2.449489742783178*alphaDrSurf[7]*fSkin[15]-1.414213562373095*(alphaDrSurf[2]*fSkin[11]+fSkin[9]*alphaDrSurf[11]+alphaDrSurf[3]*fSkin[10]+fSkin[4]*alphaDrSurf[7])); 
  Ghat[15] = -0.125*(2.449489742783178*(alphaDrSurf[0]*fSkin[31]+alphaDrSurf[1]*fSkin[30]+alphaDrSurf[2]*fSkin[29]+alphaDrSurf[3]*fSkin[28])-1.414213562373095*alphaDrSurf[0]*fSkin[26]+2.449489742783178*(alphaDrSurf[5]*fSkin[25]+alphaDrSurf[6]*fSkin[24]+alphaDrSurf[7]*fSkin[23])-1.414213562373095*(alphaDrSurf[1]*fSkin[19]+alphaDrSurf[2]*fSkin[18]+alphaDrSurf[3]*fSkin[17])+2.449489742783178*alphaDrSurf[11]*fSkin[15]-1.414213562373095*(alphaDrSurf[5]*fSkin[11]+fSkin[4]*alphaDrSurf[11]+alphaDrSurf[6]*fSkin[10]+alphaDrSurf[7]*fSkin[9])); 

  out[0] += -0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += -0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += -0.7071067811865475*Ghat[2]*rdv2; 
  out[3] += -0.7071067811865475*Ghat[3]*rdv2; 
  out[4] += -0.7071067811865475*Ghat[4]*rdv2; 
  out[5] += 1.224744871391589*Ghat[0]*rdv2; 
  out[6] += -0.7071067811865475*Ghat[5]*rdv2; 
  out[7] += -0.7071067811865475*Ghat[6]*rdv2; 
  out[8] += -0.7071067811865475*Ghat[7]*rdv2; 
  out[9] += -0.7071067811865475*Ghat[8]*rdv2; 
  out[10] += -0.7071067811865475*Ghat[9]*rdv2; 
  out[11] += -0.7071067811865475*Ghat[10]*rdv2; 
  out[12] += 1.224744871391589*Ghat[1]*rdv2; 
  out[13] += 1.224744871391589*Ghat[2]*rdv2; 
  out[14] += 1.224744871391589*Ghat[3]*rdv2; 
  out[15] += 1.224744871391589*Ghat[4]*rdv2; 
  out[16] += -0.7071067811865475*Ghat[11]*rdv2; 
  out[17] += -0.7071067811865475*Ghat[12]*rdv2; 
  out[18] += -0.7071067811865475*Ghat[13]*rdv2; 
  out[19] += -0.7071067811865475*Ghat[14]*rdv2; 
  out[20] += 1.224744871391589*Ghat[5]*rdv2; 
  out[21] += 1.224744871391589*Ghat[6]*rdv2; 
  out[22] += 1.224744871391589*Ghat[7]*rdv2; 
  out[23] += 1.224744871391589*Ghat[8]*rdv2; 
  out[24] += 1.224744871391589*Ghat[9]*rdv2; 
  out[25] += 1.224744871391589*Ghat[10]*rdv2; 
  out[26] += -0.7071067811865475*Ghat[15]*rdv2; 
  out[27] += 1.224744871391589*Ghat[11]*rdv2; 
  out[28] += 1.224744871391589*Ghat[12]*rdv2; 
  out[29] += 1.224744871391589*Ghat[13]*rdv2; 
  out[30] += 1.224744871391589*Ghat[14]*rdv2; 
  out[31] += 1.224744871391589*Ghat[15]*rdv2; 

  } 
} 
