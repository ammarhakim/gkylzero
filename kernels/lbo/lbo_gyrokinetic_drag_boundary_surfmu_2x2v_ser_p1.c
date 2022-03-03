#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH void lbo_gyrokinetic_drag_boundary_surfmu_2x2v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[4]:     cell-center coordinates. 
  // dxv[4]:   cell spacing. 
  // m_:        species mass.
  // bmag_inv:  1/(magnetic field magnitude). 
  // nuSum:     collisionalities added (self and cross species collisionalities). 
  // nuUSum[8]:sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum: sum of thermal speeds squared time their respective collisionalities. 
  // fSkin/Edge:    Distribution function in cells 
  // out:           Incremented distribution function in cell 
  double rdv2 = 2.0/dxv[3]; 

  double alphaDrSurf[8] = {0.0}; 
  double fUpwindQuad[8] = {0.0};
  double fUpwind[8] = {0.0};;
  double drag_incr[8] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = 0.7071067811865475*(nuSum[0]*(2.0*w[2]+dxv[2])-2.0*nuUSum[0]); 
  alphaDrSurf[1] = 0.7071067811865475*(nuSum[1]*(2.0*w[2]+dxv[2])-2.0*nuUSum[1]); 
  alphaDrSurf[2] = 0.7071067811865475*(2.0*nuSum[2]*w[2]-2.0*nuUSum[2]+dxv[2]*nuSum[2]); 
  alphaDrSurf[4] = -0.7071067811865475*(2.0*nuUSum[3]+((-2.0*w[2])-1.0*dxv[2])*nuSum[3]); 

  drag_incr[0] = (-0.4330127018922193*alphaDrSurf[4]*fEdge[12])+alphaDrSurf[2]*(0.25*fEdge[2]-0.4330127018922193*fEdge[9])+alphaDrSurf[1]*(0.25*fEdge[1]-0.4330127018922193*fEdge[8])+0.25*alphaDrSurf[4]*fEdge[5]+alphaDrSurf[0]*(0.25*fEdge[0]-0.4330127018922193*fEdge[4]); 
  drag_incr[1] = alphaDrSurf[2]*(0.25*fEdge[5]-0.4330127018922193*fEdge[12])-0.4330127018922193*alphaDrSurf[4]*fEdge[9]+alphaDrSurf[0]*(0.25*fEdge[1]-0.4330127018922193*fEdge[8])+alphaDrSurf[1]*(0.25*fEdge[0]-0.4330127018922193*fEdge[4])+0.25*fEdge[2]*alphaDrSurf[4]; 
  drag_incr[2] = alphaDrSurf[1]*(0.25*fEdge[5]-0.4330127018922193*fEdge[12])+alphaDrSurf[0]*(0.25*fEdge[2]-0.4330127018922193*fEdge[9])-0.4330127018922193*alphaDrSurf[4]*fEdge[8]+alphaDrSurf[2]*(0.25*fEdge[0]-0.4330127018922193*fEdge[4])+0.25*fEdge[1]*alphaDrSurf[4]; 
  drag_incr[3] = (-0.4330127018922193*alphaDrSurf[4]*fEdge[15])+alphaDrSurf[2]*(0.25*fEdge[7]-0.4330127018922193*fEdge[14])+alphaDrSurf[1]*(0.25*fEdge[6]-0.4330127018922193*fEdge[13])+0.25*alphaDrSurf[4]*fEdge[11]+alphaDrSurf[0]*(0.25*fEdge[3]-0.4330127018922193*fEdge[10]); 
  drag_incr[4] = alphaDrSurf[0]*(0.25*fEdge[5]-0.4330127018922193*fEdge[12])+alphaDrSurf[1]*(0.25*fEdge[2]-0.4330127018922193*fEdge[9])+alphaDrSurf[2]*(0.25*fEdge[1]-0.4330127018922193*fEdge[8])-0.4330127018922193*alphaDrSurf[4]*fEdge[4]+0.25*fEdge[0]*alphaDrSurf[4]; 
  drag_incr[5] = alphaDrSurf[2]*(0.25*fEdge[11]-0.4330127018922193*fEdge[15])-0.4330127018922193*alphaDrSurf[4]*fEdge[14]+alphaDrSurf[0]*(0.25*fEdge[6]-0.4330127018922193*fEdge[13])+alphaDrSurf[1]*(0.25*fEdge[3]-0.4330127018922193*fEdge[10])+0.25*alphaDrSurf[4]*fEdge[7]; 
  drag_incr[6] = alphaDrSurf[1]*(0.25*fEdge[11]-0.4330127018922193*fEdge[15])+alphaDrSurf[0]*(0.25*fEdge[7]-0.4330127018922193*fEdge[14])-0.4330127018922193*alphaDrSurf[4]*fEdge[13]+alphaDrSurf[2]*(0.25*fEdge[3]-0.4330127018922193*fEdge[10])+0.25*alphaDrSurf[4]*fEdge[6]; 
  drag_incr[7] = alphaDrSurf[0]*(0.25*fEdge[11]-0.4330127018922193*fEdge[15])+alphaDrSurf[1]*(0.25*fEdge[7]-0.4330127018922193*fEdge[14])+alphaDrSurf[2]*(0.25*fEdge[6]-0.4330127018922193*fEdge[13])-0.4330127018922193*alphaDrSurf[4]*fEdge[10]+0.25*fEdge[3]*alphaDrSurf[4]; 

  out[0] += 0.7071067811865475*drag_incr[0]*rdv2; 
  out[1] += 0.7071067811865475*drag_incr[1]*rdv2; 
  out[2] += 0.7071067811865475*drag_incr[2]*rdv2; 
  out[3] += 0.7071067811865475*drag_incr[3]*rdv2; 
  out[4] += 1.224744871391589*drag_incr[0]*rdv2; 
  out[5] += 0.7071067811865475*drag_incr[4]*rdv2; 
  out[6] += 0.7071067811865475*drag_incr[5]*rdv2; 
  out[7] += 0.7071067811865475*drag_incr[6]*rdv2; 
  out[8] += 1.224744871391589*drag_incr[1]*rdv2; 
  out[9] += 1.224744871391589*drag_incr[2]*rdv2; 
  out[10] += 1.224744871391589*drag_incr[3]*rdv2; 
  out[11] += 0.7071067811865475*drag_incr[7]*rdv2; 
  out[12] += 1.224744871391589*drag_incr[4]*rdv2; 
  out[13] += 1.224744871391589*drag_incr[5]*rdv2; 
  out[14] += 1.224744871391589*drag_incr[6]*rdv2; 
  out[15] += 1.224744871391589*drag_incr[7]*rdv2; 

  } else { 

  alphaDrSurf[0] = 0.7071067811865475*(nuSum[0]*(2.0*w[2]-1.0*dxv[2])-2.0*nuUSum[0]); 
  alphaDrSurf[1] = 0.7071067811865475*(nuSum[1]*(2.0*w[2]-1.0*dxv[2])-2.0*nuUSum[1]); 
  alphaDrSurf[2] = 0.7071067811865475*(2.0*nuSum[2]*w[2]-2.0*nuUSum[2]-1.0*dxv[2]*nuSum[2]); 
  alphaDrSurf[4] = -0.7071067811865475*(2.0*nuUSum[3]+(dxv[2]-2.0*w[2])*nuSum[3]); 

  drag_incr[0] = (-0.4330127018922193*alphaDrSurf[4]*fSkin[12])+alphaDrSurf[2]*(0.25*fSkin[2]-0.4330127018922193*fSkin[9])+alphaDrSurf[1]*(0.25*fSkin[1]-0.4330127018922193*fSkin[8])+0.25*alphaDrSurf[4]*fSkin[5]+alphaDrSurf[0]*(0.25*fSkin[0]-0.4330127018922193*fSkin[4]); 
  drag_incr[1] = alphaDrSurf[2]*(0.25*fSkin[5]-0.4330127018922193*fSkin[12])-0.4330127018922193*alphaDrSurf[4]*fSkin[9]+alphaDrSurf[0]*(0.25*fSkin[1]-0.4330127018922193*fSkin[8])+alphaDrSurf[1]*(0.25*fSkin[0]-0.4330127018922193*fSkin[4])+0.25*fSkin[2]*alphaDrSurf[4]; 
  drag_incr[2] = alphaDrSurf[1]*(0.25*fSkin[5]-0.4330127018922193*fSkin[12])+alphaDrSurf[0]*(0.25*fSkin[2]-0.4330127018922193*fSkin[9])-0.4330127018922193*alphaDrSurf[4]*fSkin[8]+alphaDrSurf[2]*(0.25*fSkin[0]-0.4330127018922193*fSkin[4])+0.25*fSkin[1]*alphaDrSurf[4]; 
  drag_incr[3] = (-0.4330127018922193*alphaDrSurf[4]*fSkin[15])+alphaDrSurf[2]*(0.25*fSkin[7]-0.4330127018922193*fSkin[14])+alphaDrSurf[1]*(0.25*fSkin[6]-0.4330127018922193*fSkin[13])+0.25*alphaDrSurf[4]*fSkin[11]+alphaDrSurf[0]*(0.25*fSkin[3]-0.4330127018922193*fSkin[10]); 
  drag_incr[4] = alphaDrSurf[0]*(0.25*fSkin[5]-0.4330127018922193*fSkin[12])+alphaDrSurf[1]*(0.25*fSkin[2]-0.4330127018922193*fSkin[9])+alphaDrSurf[2]*(0.25*fSkin[1]-0.4330127018922193*fSkin[8])-0.4330127018922193*alphaDrSurf[4]*fSkin[4]+0.25*fSkin[0]*alphaDrSurf[4]; 
  drag_incr[5] = alphaDrSurf[2]*(0.25*fSkin[11]-0.4330127018922193*fSkin[15])-0.4330127018922193*alphaDrSurf[4]*fSkin[14]+alphaDrSurf[0]*(0.25*fSkin[6]-0.4330127018922193*fSkin[13])+alphaDrSurf[1]*(0.25*fSkin[3]-0.4330127018922193*fSkin[10])+0.25*alphaDrSurf[4]*fSkin[7]; 
  drag_incr[6] = alphaDrSurf[1]*(0.25*fSkin[11]-0.4330127018922193*fSkin[15])+alphaDrSurf[0]*(0.25*fSkin[7]-0.4330127018922193*fSkin[14])-0.4330127018922193*alphaDrSurf[4]*fSkin[13]+alphaDrSurf[2]*(0.25*fSkin[3]-0.4330127018922193*fSkin[10])+0.25*alphaDrSurf[4]*fSkin[6]; 
  drag_incr[7] = alphaDrSurf[0]*(0.25*fSkin[11]-0.4330127018922193*fSkin[15])+alphaDrSurf[1]*(0.25*fSkin[7]-0.4330127018922193*fSkin[14])+alphaDrSurf[2]*(0.25*fSkin[6]-0.4330127018922193*fSkin[13])-0.4330127018922193*alphaDrSurf[4]*fSkin[10]+0.25*fSkin[3]*alphaDrSurf[4]; 

  out[0] += -0.7071067811865475*drag_incr[0]*rdv2; 
  out[1] += -0.7071067811865475*drag_incr[1]*rdv2; 
  out[2] += -0.7071067811865475*drag_incr[2]*rdv2; 
  out[3] += -0.7071067811865475*drag_incr[3]*rdv2; 
  out[4] += 1.224744871391589*drag_incr[0]*rdv2; 
  out[5] += -0.7071067811865475*drag_incr[4]*rdv2; 
  out[6] += -0.7071067811865475*drag_incr[5]*rdv2; 
  out[7] += -0.7071067811865475*drag_incr[6]*rdv2; 
  out[8] += 1.224744871391589*drag_incr[1]*rdv2; 
  out[9] += 1.224744871391589*drag_incr[2]*rdv2; 
  out[10] += 1.224744871391589*drag_incr[3]*rdv2; 
  out[11] += -0.7071067811865475*drag_incr[7]*rdv2; 
  out[12] += 1.224744871391589*drag_incr[4]*rdv2; 
  out[13] += 1.224744871391589*drag_incr[5]*rdv2; 
  out[14] += 1.224744871391589*drag_incr[6]*rdv2; 
  out[15] += 1.224744871391589*drag_incr[7]*rdv2; 

  } 
} 
