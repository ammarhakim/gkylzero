#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH void lbo_gyrokinetic_drag_boundary_surfmu_1x2v_tensor_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[3]:     cell-center coordinates. 
  // dxv[3]:   cell spacing. 
  // m_:        species mass.
  // bmag_inv:  1/(magnetic field magnitude). 
  // nuSum:     collisionalities added (self and cross species collisionalities). 
  // nuUSum[6]:sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum: sum of thermal speeds squared time their respective collisionalities. 
  // fSkin/Edge:    Distribution function in cells 
  // out:           Incremented distribution function in cell 
  double rdv2 = 2.0/dxv[2]; 

  double alphaDrSurf[9] = {0.0}; 
  double Ghat[9] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = nuSum[0]*(2.828427124746191*w[2]+1.414213562373095*dxv[2]); 
  alphaDrSurf[1] = nuSum[1]*(2.828427124746191*w[2]+1.414213562373095*dxv[2]); 
  alphaDrSurf[4] = nuSum[2]*(2.828427124746191*w[2]+1.414213562373095*dxv[2]); 

  Ghat[0] = 0.7905694150420947*alphaDrSurf[4]*fEdge[21]+0.7905694150420948*alphaDrSurf[1]*fEdge[15]-0.6123724356957944*alphaDrSurf[4]*fEdge[13]+0.7905694150420947*alphaDrSurf[0]*fEdge[9]+0.3535533905932737*alphaDrSurf[4]*fEdge[7]-0.6123724356957944*alphaDrSurf[1]*fEdge[5]-0.6123724356957944*alphaDrSurf[0]*fEdge[3]+0.3535533905932737*alphaDrSurf[1]*fEdge[1]+0.3535533905932737*alphaDrSurf[0]*fEdge[0]; 
  Ghat[1] = 0.7071067811865475*alphaDrSurf[1]*fEdge[21]+0.7071067811865475*alphaDrSurf[4]*fEdge[15]+0.7905694150420948*alphaDrSurf[0]*fEdge[15]-0.5477225575051661*alphaDrSurf[1]*fEdge[13]+0.7905694150420947*alphaDrSurf[1]*fEdge[9]+0.3162277660168379*alphaDrSurf[1]*fEdge[7]-0.5477225575051661*alphaDrSurf[4]*fEdge[5]-0.6123724356957944*alphaDrSurf[0]*fEdge[5]+0.3162277660168379*fEdge[1]*alphaDrSurf[4]-0.6123724356957944*alphaDrSurf[1]*fEdge[3]+0.3535533905932737*alphaDrSurf[0]*fEdge[1]+0.3535533905932737*fEdge[0]*alphaDrSurf[1]; 
  Ghat[2] = 0.7905694150420947*alphaDrSurf[4]*fEdge[24]+0.7905694150420947*alphaDrSurf[1]*fEdge[19]-0.6123724356957944*alphaDrSurf[4]*fEdge[17]+0.7905694150420948*alphaDrSurf[0]*fEdge[16]+0.3535533905932737*alphaDrSurf[4]*fEdge[11]-0.6123724356957944*alphaDrSurf[1]*fEdge[10]-0.6123724356957944*alphaDrSurf[0]*fEdge[6]+0.3535533905932737*alphaDrSurf[1]*fEdge[4]+0.3535533905932737*alphaDrSurf[0]*fEdge[2]; 
  Ghat[3] = 0.7071067811865475*alphaDrSurf[1]*fEdge[24]+0.7071067811865475*alphaDrSurf[4]*fEdge[19]+0.7905694150420947*alphaDrSurf[0]*fEdge[19]-0.5477225575051661*alphaDrSurf[1]*fEdge[17]+0.7905694150420948*alphaDrSurf[1]*fEdge[16]+0.3162277660168379*alphaDrSurf[1]*fEdge[11]-0.5477225575051661*alphaDrSurf[4]*fEdge[10]-0.6123724356957944*alphaDrSurf[0]*fEdge[10]-0.6123724356957944*alphaDrSurf[1]*fEdge[6]+0.3162277660168379*alphaDrSurf[4]*fEdge[4]+0.3535533905932737*alphaDrSurf[0]*fEdge[4]+0.3535533905932737*alphaDrSurf[1]*fEdge[2]; 
  Ghat[4] = 0.5050762722761053*alphaDrSurf[4]*fEdge[21]+0.7905694150420947*alphaDrSurf[0]*fEdge[21]+0.7071067811865475*alphaDrSurf[1]*fEdge[15]-0.3912303982179757*alphaDrSurf[4]*fEdge[13]-0.6123724356957944*alphaDrSurf[0]*fEdge[13]+0.7905694150420947*alphaDrSurf[4]*fEdge[9]+0.2258769757263128*alphaDrSurf[4]*fEdge[7]+0.3535533905932737*alphaDrSurf[0]*fEdge[7]-0.5477225575051661*alphaDrSurf[1]*fEdge[5]-0.6123724356957944*fEdge[3]*alphaDrSurf[4]+0.3535533905932737*fEdge[0]*alphaDrSurf[4]+0.3162277660168379*alphaDrSurf[1]*fEdge[1]; 
  Ghat[5] = 0.7905694150420947*alphaDrSurf[4]*fEdge[26]+0.7905694150420947*alphaDrSurf[1]*fEdge[25]-0.6123724356957944*alphaDrSurf[4]*fEdge[23]+0.7905694150420947*alphaDrSurf[0]*fEdge[22]+0.3535533905932737*alphaDrSurf[4]*fEdge[20]-0.6123724356957944*alphaDrSurf[1]*fEdge[18]-0.6123724356957944*alphaDrSurf[0]*fEdge[14]+0.3535533905932737*alphaDrSurf[1]*fEdge[12]+0.3535533905932737*alphaDrSurf[0]*fEdge[8]; 
  Ghat[6] = 0.5050762722761054*alphaDrSurf[4]*fEdge[24]+0.7905694150420948*alphaDrSurf[0]*fEdge[24]+0.7071067811865475*alphaDrSurf[1]*fEdge[19]-0.3912303982179757*alphaDrSurf[4]*fEdge[17]-0.6123724356957944*alphaDrSurf[0]*fEdge[17]+0.7905694150420947*alphaDrSurf[4]*fEdge[16]+0.2258769757263128*alphaDrSurf[4]*fEdge[11]+0.3535533905932737*alphaDrSurf[0]*fEdge[11]-0.5477225575051661*alphaDrSurf[1]*fEdge[10]-0.6123724356957944*alphaDrSurf[4]*fEdge[6]+0.3162277660168379*alphaDrSurf[1]*fEdge[4]+0.3535533905932737*fEdge[2]*alphaDrSurf[4]; 
  Ghat[7] = 0.7071067811865475*alphaDrSurf[1]*fEdge[26]+0.7071067811865475*alphaDrSurf[4]*fEdge[25]+0.7905694150420948*alphaDrSurf[0]*fEdge[25]-0.5477225575051661*alphaDrSurf[1]*fEdge[23]+0.7905694150420948*alphaDrSurf[1]*fEdge[22]+0.3162277660168379*alphaDrSurf[1]*fEdge[20]-0.5477225575051661*alphaDrSurf[4]*fEdge[18]-0.6123724356957944*alphaDrSurf[0]*fEdge[18]-0.6123724356957944*alphaDrSurf[1]*fEdge[14]+0.3162277660168379*alphaDrSurf[4]*fEdge[12]+0.3535533905932737*alphaDrSurf[0]*fEdge[12]+0.3535533905932737*alphaDrSurf[1]*fEdge[8]; 
  Ghat[8] = 0.5050762722761053*alphaDrSurf[4]*fEdge[26]+0.7905694150420947*alphaDrSurf[0]*fEdge[26]+0.7071067811865475*alphaDrSurf[1]*fEdge[25]-0.3912303982179757*alphaDrSurf[4]*fEdge[23]-0.6123724356957944*alphaDrSurf[0]*fEdge[23]+0.7905694150420947*alphaDrSurf[4]*fEdge[22]+0.2258769757263128*alphaDrSurf[4]*fEdge[20]+0.3535533905932737*alphaDrSurf[0]*fEdge[20]-0.5477225575051661*alphaDrSurf[1]*fEdge[18]-0.6123724356957944*alphaDrSurf[4]*fEdge[14]+0.3162277660168379*alphaDrSurf[1]*fEdge[12]+0.3535533905932737*alphaDrSurf[4]*fEdge[8]; 

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
  out[20] += 0.7071067811865475*Ghat[8]*rdv2; 
  out[21] += 1.58113883008419*Ghat[4]*rdv2; 
  out[22] += 1.58113883008419*Ghat[5]*rdv2; 
  out[23] += 1.224744871391589*Ghat[8]*rdv2; 
  out[24] += 1.58113883008419*Ghat[6]*rdv2; 
  out[25] += 1.58113883008419*Ghat[7]*rdv2; 
  out[26] += 1.58113883008419*Ghat[8]*rdv2; 

  } else { 

  alphaDrSurf[0] = nuSum[0]*(2.828427124746191*w[2]-1.414213562373095*dxv[2]); 
  alphaDrSurf[1] = nuSum[1]*(2.828427124746191*w[2]-1.414213562373095*dxv[2]); 
  alphaDrSurf[4] = nuSum[2]*(2.828427124746191*w[2]-1.414213562373095*dxv[2]); 

  Ghat[0] = 0.7905694150420947*alphaDrSurf[4]*fSkin[21]+0.7905694150420948*alphaDrSurf[1]*fSkin[15]-0.6123724356957944*alphaDrSurf[4]*fSkin[13]+0.7905694150420947*alphaDrSurf[0]*fSkin[9]+0.3535533905932737*alphaDrSurf[4]*fSkin[7]-0.6123724356957944*alphaDrSurf[1]*fSkin[5]-0.6123724356957944*alphaDrSurf[0]*fSkin[3]+0.3535533905932737*alphaDrSurf[1]*fSkin[1]+0.3535533905932737*alphaDrSurf[0]*fSkin[0]; 
  Ghat[1] = 0.7071067811865475*alphaDrSurf[1]*fSkin[21]+0.7071067811865475*alphaDrSurf[4]*fSkin[15]+0.7905694150420948*alphaDrSurf[0]*fSkin[15]-0.5477225575051661*alphaDrSurf[1]*fSkin[13]+0.7905694150420947*alphaDrSurf[1]*fSkin[9]+0.3162277660168379*alphaDrSurf[1]*fSkin[7]-0.5477225575051661*alphaDrSurf[4]*fSkin[5]-0.6123724356957944*alphaDrSurf[0]*fSkin[5]+0.3162277660168379*fSkin[1]*alphaDrSurf[4]-0.6123724356957944*alphaDrSurf[1]*fSkin[3]+0.3535533905932737*alphaDrSurf[0]*fSkin[1]+0.3535533905932737*fSkin[0]*alphaDrSurf[1]; 
  Ghat[2] = 0.7905694150420947*alphaDrSurf[4]*fSkin[24]+0.7905694150420947*alphaDrSurf[1]*fSkin[19]-0.6123724356957944*alphaDrSurf[4]*fSkin[17]+0.7905694150420948*alphaDrSurf[0]*fSkin[16]+0.3535533905932737*alphaDrSurf[4]*fSkin[11]-0.6123724356957944*alphaDrSurf[1]*fSkin[10]-0.6123724356957944*alphaDrSurf[0]*fSkin[6]+0.3535533905932737*alphaDrSurf[1]*fSkin[4]+0.3535533905932737*alphaDrSurf[0]*fSkin[2]; 
  Ghat[3] = 0.7071067811865475*alphaDrSurf[1]*fSkin[24]+0.7071067811865475*alphaDrSurf[4]*fSkin[19]+0.7905694150420947*alphaDrSurf[0]*fSkin[19]-0.5477225575051661*alphaDrSurf[1]*fSkin[17]+0.7905694150420948*alphaDrSurf[1]*fSkin[16]+0.3162277660168379*alphaDrSurf[1]*fSkin[11]-0.5477225575051661*alphaDrSurf[4]*fSkin[10]-0.6123724356957944*alphaDrSurf[0]*fSkin[10]-0.6123724356957944*alphaDrSurf[1]*fSkin[6]+0.3162277660168379*alphaDrSurf[4]*fSkin[4]+0.3535533905932737*alphaDrSurf[0]*fSkin[4]+0.3535533905932737*alphaDrSurf[1]*fSkin[2]; 
  Ghat[4] = 0.5050762722761053*alphaDrSurf[4]*fSkin[21]+0.7905694150420947*alphaDrSurf[0]*fSkin[21]+0.7071067811865475*alphaDrSurf[1]*fSkin[15]-0.3912303982179757*alphaDrSurf[4]*fSkin[13]-0.6123724356957944*alphaDrSurf[0]*fSkin[13]+0.7905694150420947*alphaDrSurf[4]*fSkin[9]+0.2258769757263128*alphaDrSurf[4]*fSkin[7]+0.3535533905932737*alphaDrSurf[0]*fSkin[7]-0.5477225575051661*alphaDrSurf[1]*fSkin[5]-0.6123724356957944*fSkin[3]*alphaDrSurf[4]+0.3535533905932737*fSkin[0]*alphaDrSurf[4]+0.3162277660168379*alphaDrSurf[1]*fSkin[1]; 
  Ghat[5] = 0.7905694150420947*alphaDrSurf[4]*fSkin[26]+0.7905694150420947*alphaDrSurf[1]*fSkin[25]-0.6123724356957944*alphaDrSurf[4]*fSkin[23]+0.7905694150420947*alphaDrSurf[0]*fSkin[22]+0.3535533905932737*alphaDrSurf[4]*fSkin[20]-0.6123724356957944*alphaDrSurf[1]*fSkin[18]-0.6123724356957944*alphaDrSurf[0]*fSkin[14]+0.3535533905932737*alphaDrSurf[1]*fSkin[12]+0.3535533905932737*alphaDrSurf[0]*fSkin[8]; 
  Ghat[6] = 0.5050762722761054*alphaDrSurf[4]*fSkin[24]+0.7905694150420948*alphaDrSurf[0]*fSkin[24]+0.7071067811865475*alphaDrSurf[1]*fSkin[19]-0.3912303982179757*alphaDrSurf[4]*fSkin[17]-0.6123724356957944*alphaDrSurf[0]*fSkin[17]+0.7905694150420947*alphaDrSurf[4]*fSkin[16]+0.2258769757263128*alphaDrSurf[4]*fSkin[11]+0.3535533905932737*alphaDrSurf[0]*fSkin[11]-0.5477225575051661*alphaDrSurf[1]*fSkin[10]-0.6123724356957944*alphaDrSurf[4]*fSkin[6]+0.3162277660168379*alphaDrSurf[1]*fSkin[4]+0.3535533905932737*fSkin[2]*alphaDrSurf[4]; 
  Ghat[7] = 0.7071067811865475*alphaDrSurf[1]*fSkin[26]+0.7071067811865475*alphaDrSurf[4]*fSkin[25]+0.7905694150420948*alphaDrSurf[0]*fSkin[25]-0.5477225575051661*alphaDrSurf[1]*fSkin[23]+0.7905694150420948*alphaDrSurf[1]*fSkin[22]+0.3162277660168379*alphaDrSurf[1]*fSkin[20]-0.5477225575051661*alphaDrSurf[4]*fSkin[18]-0.6123724356957944*alphaDrSurf[0]*fSkin[18]-0.6123724356957944*alphaDrSurf[1]*fSkin[14]+0.3162277660168379*alphaDrSurf[4]*fSkin[12]+0.3535533905932737*alphaDrSurf[0]*fSkin[12]+0.3535533905932737*alphaDrSurf[1]*fSkin[8]; 
  Ghat[8] = 0.5050762722761053*alphaDrSurf[4]*fSkin[26]+0.7905694150420947*alphaDrSurf[0]*fSkin[26]+0.7071067811865475*alphaDrSurf[1]*fSkin[25]-0.3912303982179757*alphaDrSurf[4]*fSkin[23]-0.6123724356957944*alphaDrSurf[0]*fSkin[23]+0.7905694150420947*alphaDrSurf[4]*fSkin[22]+0.2258769757263128*alphaDrSurf[4]*fSkin[20]+0.3535533905932737*alphaDrSurf[0]*fSkin[20]-0.5477225575051661*alphaDrSurf[1]*fSkin[18]-0.6123724356957944*alphaDrSurf[4]*fSkin[14]+0.3162277660168379*alphaDrSurf[1]*fSkin[12]+0.3535533905932737*alphaDrSurf[4]*fSkin[8]; 

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
  out[20] += -0.7071067811865475*Ghat[8]*rdv2; 
  out[21] += -1.58113883008419*Ghat[4]*rdv2; 
  out[22] += -1.58113883008419*Ghat[5]*rdv2; 
  out[23] += 1.224744871391589*Ghat[8]*rdv2; 
  out[24] += -1.58113883008419*Ghat[6]*rdv2; 
  out[25] += -1.58113883008419*Ghat[7]*rdv2; 
  out[26] += -1.58113883008419*Ghat[8]*rdv2; 

  } 
} 
