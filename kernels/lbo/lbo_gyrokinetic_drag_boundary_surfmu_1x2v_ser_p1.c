#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH void lbo_gyrokinetic_drag_boundary_surfmu_1x2v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[3]:     cell-center coordinates. 
  // dxv[3]:   cell spacing. 
  // m_:        species mass.
  // bmag_inv:  1/(magnetic field magnitude). 
  // nuSum:     collisionalities added (self and cross species collisionalities). 
  // nuUSum[4]:sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum: sum of thermal speeds squared time their respective collisionalities. 
  // fSkin/Edge:    Distribution function in cells 
  // out:           Incremented distribution function in cell 
  double rdv2 = 2.0/dxv[2]; 

  double alphaDrSurf[4] = {0.0}; 
  double fUpwindQuad[4] = {0.0};
  double fUpwind[4] = {0.0};;
  double drag_incr[4] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = 0.5*(nuSum[0]*(2.828427124746191*w[1]+1.414213562373095*dxv[1])-2.828427124746191*nuUSum[0]); 
  alphaDrSurf[1] = 0.5*(2.828427124746191*nuSum[1]*w[1]-2.828427124746191*nuUSum[1]+1.414213562373095*dxv[1]*nuSum[1]); 

  drag_incr[0] = alphaDrSurf[1]*(0.3535533905932737*fEdge[1]-0.6123724356957944*fEdge[5])+alphaDrSurf[0]*(0.3535533905932737*fEdge[0]-0.6123724356957944*fEdge[3]); 
  drag_incr[1] = alphaDrSurf[0]*(0.3535533905932737*fEdge[1]-0.6123724356957944*fEdge[5])+alphaDrSurf[1]*(0.3535533905932737*fEdge[0]-0.6123724356957944*fEdge[3]); 
  drag_incr[2] = alphaDrSurf[1]*(0.3535533905932737*fEdge[4]-0.6123724356957944*fEdge[7])+alphaDrSurf[0]*(0.3535533905932737*fEdge[2]-0.6123724356957944*fEdge[6]); 
  drag_incr[3] = alphaDrSurf[0]*(0.3535533905932737*fEdge[4]-0.6123724356957944*fEdge[7])+alphaDrSurf[1]*(0.3535533905932737*fEdge[2]-0.6123724356957944*fEdge[6]); 

  out[0] += 0.7071067811865475*drag_incr[0]*rdv2; 
  out[1] += 0.7071067811865475*drag_incr[1]*rdv2; 
  out[2] += 0.7071067811865475*drag_incr[2]*rdv2; 
  out[3] += 1.224744871391589*drag_incr[0]*rdv2; 
  out[4] += 0.7071067811865475*drag_incr[3]*rdv2; 
  out[5] += 1.224744871391589*drag_incr[1]*rdv2; 
  out[6] += 1.224744871391589*drag_incr[2]*rdv2; 
  out[7] += 1.224744871391589*drag_incr[3]*rdv2; 

  } else { 

  alphaDrSurf[0] = 0.5*(nuSum[0]*(2.828427124746191*w[1]-1.414213562373095*dxv[1])-2.828427124746191*nuUSum[0]); 
  alphaDrSurf[1] = 0.5*(2.828427124746191*nuSum[1]*w[1]-2.828427124746191*nuUSum[1]-1.414213562373095*dxv[1]*nuSum[1]); 

  drag_incr[0] = alphaDrSurf[1]*(0.3535533905932737*fSkin[1]-0.6123724356957944*fSkin[5])+alphaDrSurf[0]*(0.3535533905932737*fSkin[0]-0.6123724356957944*fSkin[3]); 
  drag_incr[1] = alphaDrSurf[0]*(0.3535533905932737*fSkin[1]-0.6123724356957944*fSkin[5])+alphaDrSurf[1]*(0.3535533905932737*fSkin[0]-0.6123724356957944*fSkin[3]); 
  drag_incr[2] = alphaDrSurf[1]*(0.3535533905932737*fSkin[4]-0.6123724356957944*fSkin[7])+alphaDrSurf[0]*(0.3535533905932737*fSkin[2]-0.6123724356957944*fSkin[6]); 
  drag_incr[3] = alphaDrSurf[0]*(0.3535533905932737*fSkin[4]-0.6123724356957944*fSkin[7])+alphaDrSurf[1]*(0.3535533905932737*fSkin[2]-0.6123724356957944*fSkin[6]); 

  out[0] += -0.7071067811865475*drag_incr[0]*rdv2; 
  out[1] += -0.7071067811865475*drag_incr[1]*rdv2; 
  out[2] += -0.7071067811865475*drag_incr[2]*rdv2; 
  out[3] += 1.224744871391589*drag_incr[0]*rdv2; 
  out[4] += -0.7071067811865475*drag_incr[3]*rdv2; 
  out[5] += 1.224744871391589*drag_incr[1]*rdv2; 
  out[6] += 1.224744871391589*drag_incr[2]*rdv2; 
  out[7] += 1.224744871391589*drag_incr[3]*rdv2; 

  } 
} 
