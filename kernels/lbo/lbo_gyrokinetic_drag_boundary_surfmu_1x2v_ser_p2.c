#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH void lbo_gyrokinetic_drag_boundary_surfmu_1x2v_ser_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
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

  double alphaDrSurf[8] = {0.0}; 
  double fUpwindQuad[9] = {0.0};
  double fUpwind[8] = {0.0};;
  double drag_incr[8] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = 0.5*(nuSum[0]*(2.828427124746191*w[1]+1.414213562373095*dxv[1])-2.828427124746191*nuUSum[0]); 
  alphaDrSurf[1] = 0.5*(2.828427124746191*nuSum[1]*w[1]-2.828427124746191*nuUSum[1]+1.414213562373095*dxv[1]*nuSum[1]); 
  alphaDrSurf[4] = -0.5*(2.828427124746191*nuUSum[2]+((-2.828427124746191*w[1])-1.414213562373095*dxv[1])*nuSum[2]); 

  drag_incr[0] = alphaDrSurf[1]*(0.7905694150420948*fEdge[15]-0.6123724356957944*fEdge[5]+0.3535533905932737*fEdge[1])-0.6123724356957944*alphaDrSurf[4]*fEdge[13]+alphaDrSurf[0]*(0.7905694150420947*fEdge[9]-0.6123724356957944*fEdge[3]+0.3535533905932737*fEdge[0])+0.3535533905932737*alphaDrSurf[4]*fEdge[7]; 
  drag_incr[1] = 0.7071067811865475*alphaDrSurf[4]*fEdge[15]+alphaDrSurf[0]*(0.7905694150420948*fEdge[15]-0.6123724356957944*fEdge[5]+0.3535533905932737*fEdge[1])+alphaDrSurf[1]*((-0.5477225575051661*fEdge[13])+0.7905694150420947*fEdge[9]+0.3162277660168379*fEdge[7]-0.6123724356957944*fEdge[3]+0.3535533905932737*fEdge[0])-0.5477225575051661*alphaDrSurf[4]*fEdge[5]+0.3162277660168379*fEdge[1]*alphaDrSurf[4]; 
  drag_incr[2] = alphaDrSurf[1]*(0.7905694150420947*fEdge[19]-0.6123724356957944*fEdge[10]+0.3535533905932737*fEdge[4])-0.6123724356957944*alphaDrSurf[4]*fEdge[17]+alphaDrSurf[0]*(0.7905694150420948*fEdge[16]-0.6123724356957944*fEdge[6]+0.3535533905932737*fEdge[2])+0.3535533905932737*alphaDrSurf[4]*fEdge[11]; 
  drag_incr[3] = 0.7071067811865475*alphaDrSurf[4]*fEdge[19]+alphaDrSurf[0]*(0.7905694150420947*fEdge[19]-0.6123724356957944*fEdge[10]+0.3535533905932737*fEdge[4])+alphaDrSurf[1]*((-0.5477225575051661*fEdge[17])+0.7905694150420948*fEdge[16]+0.3162277660168379*fEdge[11]-0.6123724356957944*fEdge[6]+0.3535533905932737*fEdge[2])-0.5477225575051661*alphaDrSurf[4]*fEdge[10]+0.3162277660168379*alphaDrSurf[4]*fEdge[4]; 
  drag_incr[4] = alphaDrSurf[1]*(0.7071067811865475*fEdge[15]-0.5477225575051661*fEdge[5]+0.3162277660168379*fEdge[1])-0.3912303982179757*alphaDrSurf[4]*fEdge[13]+alphaDrSurf[0]*(0.3535533905932737*fEdge[7]-0.6123724356957944*fEdge[13])+0.7905694150420947*alphaDrSurf[4]*fEdge[9]+0.2258769757263128*alphaDrSurf[4]*fEdge[7]-0.6123724356957944*fEdge[3]*alphaDrSurf[4]+0.3535533905932737*fEdge[0]*alphaDrSurf[4]; 
  drag_incr[5] = alphaDrSurf[1]*(0.3535533905932737*fEdge[12]-0.6123724356957944*fEdge[18])+alphaDrSurf[0]*(0.3535533905932737*fEdge[8]-0.6123724356957944*fEdge[14]); 
  drag_incr[6] = alphaDrSurf[1]*(0.7071067811865475*fEdge[19]-0.5477225575051661*fEdge[10]+0.3162277660168379*fEdge[4])-0.3912303982179757*alphaDrSurf[4]*fEdge[17]+alphaDrSurf[0]*(0.3535533905932737*fEdge[11]-0.6123724356957944*fEdge[17])+0.7905694150420947*alphaDrSurf[4]*fEdge[16]+0.2258769757263128*alphaDrSurf[4]*fEdge[11]-0.6123724356957944*alphaDrSurf[4]*fEdge[6]+0.3535533905932737*fEdge[2]*alphaDrSurf[4]; 
  drag_incr[7] = (-0.5477225575051661*alphaDrSurf[4]*fEdge[18])+alphaDrSurf[0]*(0.3535533905932737*fEdge[12]-0.6123724356957944*fEdge[18])+alphaDrSurf[1]*(0.3535533905932737*fEdge[8]-0.6123724356957944*fEdge[14])+0.3162277660168379*alphaDrSurf[4]*fEdge[12]; 

  out[0] += 0.7071067811865475*drag_incr[0]*rdv2; 
  out[1] += 0.7071067811865475*drag_incr[1]*rdv2; 
  out[2] += 0.7071067811865475*drag_incr[2]*rdv2; 
  out[3] += 1.224744871391589*drag_incr[0]*rdv2; 
  out[4] += 0.7071067811865475*drag_incr[3]*rdv2; 
  out[5] += 1.224744871391589*drag_incr[1]*rdv2; 
  out[6] += 1.224744871391589*drag_incr[2]*rdv2; 
  out[7] += 0.7071067811865475*drag_incr[4]*rdv2; 
  out[8] += 0.7071067811865475*drag_incr[5]*rdv2; 
  out[9] += 1.58113883008419*drag_incr[0]*rdv2; 
  out[10] += 1.224744871391589*drag_incr[3]*rdv2; 
  out[11] += 0.7071067811865475*drag_incr[6]*rdv2; 
  out[12] += 0.7071067811865475*drag_incr[7]*rdv2; 
  out[13] += 1.224744871391589*drag_incr[4]*rdv2; 
  out[14] += 1.224744871391589*drag_incr[5]*rdv2; 
  out[15] += 1.58113883008419*drag_incr[1]*rdv2; 
  out[16] += 1.58113883008419*drag_incr[2]*rdv2; 
  out[17] += 1.224744871391589*drag_incr[6]*rdv2; 
  out[18] += 1.224744871391589*drag_incr[7]*rdv2; 
  out[19] += 1.58113883008419*drag_incr[3]*rdv2; 

  } else { 

  alphaDrSurf[0] = 0.5*(nuSum[0]*(2.828427124746191*w[1]-1.414213562373095*dxv[1])-2.828427124746191*nuUSum[0]); 
  alphaDrSurf[1] = 0.5*(2.828427124746191*nuSum[1]*w[1]-2.828427124746191*nuUSum[1]-1.414213562373095*dxv[1]*nuSum[1]); 
  alphaDrSurf[4] = -0.5*(2.828427124746191*nuUSum[2]+(1.414213562373095*dxv[1]-2.828427124746191*w[1])*nuSum[2]); 

  drag_incr[0] = alphaDrSurf[1]*(0.7905694150420948*fSkin[15]-0.6123724356957944*fSkin[5]+0.3535533905932737*fSkin[1])-0.6123724356957944*alphaDrSurf[4]*fSkin[13]+alphaDrSurf[0]*(0.7905694150420947*fSkin[9]-0.6123724356957944*fSkin[3]+0.3535533905932737*fSkin[0])+0.3535533905932737*alphaDrSurf[4]*fSkin[7]; 
  drag_incr[1] = 0.7071067811865475*alphaDrSurf[4]*fSkin[15]+alphaDrSurf[0]*(0.7905694150420948*fSkin[15]-0.6123724356957944*fSkin[5]+0.3535533905932737*fSkin[1])+alphaDrSurf[1]*((-0.5477225575051661*fSkin[13])+0.7905694150420947*fSkin[9]+0.3162277660168379*fSkin[7]-0.6123724356957944*fSkin[3]+0.3535533905932737*fSkin[0])-0.5477225575051661*alphaDrSurf[4]*fSkin[5]+0.3162277660168379*fSkin[1]*alphaDrSurf[4]; 
  drag_incr[2] = alphaDrSurf[1]*(0.7905694150420947*fSkin[19]-0.6123724356957944*fSkin[10]+0.3535533905932737*fSkin[4])-0.6123724356957944*alphaDrSurf[4]*fSkin[17]+alphaDrSurf[0]*(0.7905694150420948*fSkin[16]-0.6123724356957944*fSkin[6]+0.3535533905932737*fSkin[2])+0.3535533905932737*alphaDrSurf[4]*fSkin[11]; 
  drag_incr[3] = 0.7071067811865475*alphaDrSurf[4]*fSkin[19]+alphaDrSurf[0]*(0.7905694150420947*fSkin[19]-0.6123724356957944*fSkin[10]+0.3535533905932737*fSkin[4])+alphaDrSurf[1]*((-0.5477225575051661*fSkin[17])+0.7905694150420948*fSkin[16]+0.3162277660168379*fSkin[11]-0.6123724356957944*fSkin[6]+0.3535533905932737*fSkin[2])-0.5477225575051661*alphaDrSurf[4]*fSkin[10]+0.3162277660168379*alphaDrSurf[4]*fSkin[4]; 
  drag_incr[4] = alphaDrSurf[1]*(0.7071067811865475*fSkin[15]-0.5477225575051661*fSkin[5]+0.3162277660168379*fSkin[1])-0.3912303982179757*alphaDrSurf[4]*fSkin[13]+alphaDrSurf[0]*(0.3535533905932737*fSkin[7]-0.6123724356957944*fSkin[13])+0.7905694150420947*alphaDrSurf[4]*fSkin[9]+0.2258769757263128*alphaDrSurf[4]*fSkin[7]-0.6123724356957944*fSkin[3]*alphaDrSurf[4]+0.3535533905932737*fSkin[0]*alphaDrSurf[4]; 
  drag_incr[5] = alphaDrSurf[1]*(0.3535533905932737*fSkin[12]-0.6123724356957944*fSkin[18])+alphaDrSurf[0]*(0.3535533905932737*fSkin[8]-0.6123724356957944*fSkin[14]); 
  drag_incr[6] = alphaDrSurf[1]*(0.7071067811865475*fSkin[19]-0.5477225575051661*fSkin[10]+0.3162277660168379*fSkin[4])-0.3912303982179757*alphaDrSurf[4]*fSkin[17]+alphaDrSurf[0]*(0.3535533905932737*fSkin[11]-0.6123724356957944*fSkin[17])+0.7905694150420947*alphaDrSurf[4]*fSkin[16]+0.2258769757263128*alphaDrSurf[4]*fSkin[11]-0.6123724356957944*alphaDrSurf[4]*fSkin[6]+0.3535533905932737*fSkin[2]*alphaDrSurf[4]; 
  drag_incr[7] = (-0.5477225575051661*alphaDrSurf[4]*fSkin[18])+alphaDrSurf[0]*(0.3535533905932737*fSkin[12]-0.6123724356957944*fSkin[18])+alphaDrSurf[1]*(0.3535533905932737*fSkin[8]-0.6123724356957944*fSkin[14])+0.3162277660168379*alphaDrSurf[4]*fSkin[12]; 

  out[0] += -0.7071067811865475*drag_incr[0]*rdv2; 
  out[1] += -0.7071067811865475*drag_incr[1]*rdv2; 
  out[2] += -0.7071067811865475*drag_incr[2]*rdv2; 
  out[3] += 1.224744871391589*drag_incr[0]*rdv2; 
  out[4] += -0.7071067811865475*drag_incr[3]*rdv2; 
  out[5] += 1.224744871391589*drag_incr[1]*rdv2; 
  out[6] += 1.224744871391589*drag_incr[2]*rdv2; 
  out[7] += -0.7071067811865475*drag_incr[4]*rdv2; 
  out[8] += -0.7071067811865475*drag_incr[5]*rdv2; 
  out[9] += -1.58113883008419*drag_incr[0]*rdv2; 
  out[10] += 1.224744871391589*drag_incr[3]*rdv2; 
  out[11] += -0.7071067811865475*drag_incr[6]*rdv2; 
  out[12] += -0.7071067811865475*drag_incr[7]*rdv2; 
  out[13] += 1.224744871391589*drag_incr[4]*rdv2; 
  out[14] += 1.224744871391589*drag_incr[5]*rdv2; 
  out[15] += -1.58113883008419*drag_incr[1]*rdv2; 
  out[16] += -1.58113883008419*drag_incr[2]*rdv2; 
  out[17] += 1.224744871391589*drag_incr[6]*rdv2; 
  out[18] += 1.224744871391589*drag_incr[7]*rdv2; 
  out[19] += -1.58113883008419*drag_incr[3]*rdv2; 

  } 
} 
