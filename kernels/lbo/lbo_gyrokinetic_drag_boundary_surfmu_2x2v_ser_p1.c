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
  double Ghat[8] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = 0.7071067811865475*(nuSum[0]*(2.0*w[2]+dxv[2])-2.0*nuUSum[0]); 
  alphaDrSurf[1] = 0.7071067811865475*(nuSum[1]*(2.0*w[2]+dxv[2])-2.0*nuUSum[1]); 
  alphaDrSurf[2] = 0.7071067811865475*(2.0*nuSum[2]*w[2]-2.0*nuUSum[2]+dxv[2]*nuSum[2]); 
  alphaDrSurf[4] = -0.7071067811865475*(2.0*nuUSum[3]+((-2.0*w[2])-1.0*dxv[2])*nuSum[3]); 

  Ghat[0] = -0.25*(1.732050807568877*(alphaDrSurf[4]*fEdge[12]+alphaDrSurf[2]*fEdge[9]+alphaDrSurf[1]*fEdge[8])-1.0*alphaDrSurf[4]*fEdge[5]+1.732050807568877*alphaDrSurf[0]*fEdge[4]-1.0*(alphaDrSurf[2]*fEdge[2]+alphaDrSurf[1]*fEdge[1]+alphaDrSurf[0]*fEdge[0])); 
  Ghat[1] = -0.25*(1.732050807568877*(alphaDrSurf[2]*fEdge[12]+alphaDrSurf[4]*fEdge[9]+alphaDrSurf[0]*fEdge[8])-1.0*alphaDrSurf[2]*fEdge[5]+1.732050807568877*alphaDrSurf[1]*fEdge[4]-1.0*(fEdge[2]*alphaDrSurf[4]+alphaDrSurf[0]*fEdge[1]+fEdge[0]*alphaDrSurf[1])); 
  Ghat[2] = -0.25*(1.732050807568877*(alphaDrSurf[1]*fEdge[12]+alphaDrSurf[0]*fEdge[9]+alphaDrSurf[4]*fEdge[8])-1.0*alphaDrSurf[1]*fEdge[5]+1.732050807568877*alphaDrSurf[2]*fEdge[4]-1.0*(fEdge[1]*alphaDrSurf[4]+alphaDrSurf[0]*fEdge[2]+fEdge[0]*alphaDrSurf[2])); 
  Ghat[3] = -0.25*(1.732050807568877*(alphaDrSurf[4]*fEdge[15]+alphaDrSurf[2]*fEdge[14]+alphaDrSurf[1]*fEdge[13])-1.0*alphaDrSurf[4]*fEdge[11]+1.732050807568877*alphaDrSurf[0]*fEdge[10]-1.0*(alphaDrSurf[2]*fEdge[7]+alphaDrSurf[1]*fEdge[6]+alphaDrSurf[0]*fEdge[3])); 
  Ghat[4] = -0.25*(1.732050807568877*(alphaDrSurf[0]*fEdge[12]+alphaDrSurf[1]*fEdge[9]+alphaDrSurf[2]*fEdge[8])-1.0*alphaDrSurf[0]*fEdge[5]+1.732050807568877*alphaDrSurf[4]*fEdge[4]-1.0*(fEdge[0]*alphaDrSurf[4]+alphaDrSurf[1]*fEdge[2]+fEdge[1]*alphaDrSurf[2])); 
  Ghat[5] = -0.25*(1.732050807568877*(alphaDrSurf[2]*fEdge[15]+alphaDrSurf[4]*fEdge[14]+alphaDrSurf[0]*fEdge[13])-1.0*alphaDrSurf[2]*fEdge[11]+1.732050807568877*alphaDrSurf[1]*fEdge[10]-1.0*(alphaDrSurf[4]*fEdge[7]+alphaDrSurf[0]*fEdge[6]+alphaDrSurf[1]*fEdge[3])); 
  Ghat[6] = -0.25*(1.732050807568877*(alphaDrSurf[1]*fEdge[15]+alphaDrSurf[0]*fEdge[14]+alphaDrSurf[4]*fEdge[13])-1.0*alphaDrSurf[1]*fEdge[11]+1.732050807568877*alphaDrSurf[2]*fEdge[10]-1.0*(alphaDrSurf[0]*fEdge[7]+alphaDrSurf[4]*fEdge[6]+alphaDrSurf[2]*fEdge[3])); 
  Ghat[7] = -0.25*(1.732050807568877*(alphaDrSurf[0]*fEdge[15]+alphaDrSurf[1]*fEdge[14]+alphaDrSurf[2]*fEdge[13])-1.0*alphaDrSurf[0]*fEdge[11]+1.732050807568877*alphaDrSurf[4]*fEdge[10]-1.0*(alphaDrSurf[1]*fEdge[7]+alphaDrSurf[2]*fEdge[6]+fEdge[3]*alphaDrSurf[4])); 

  out[0] += 0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += 0.7071067811865475*Ghat[2]*rdv2; 
  out[3] += 0.7071067811865475*Ghat[3]*rdv2; 
  out[4] += 1.224744871391589*Ghat[0]*rdv2; 
  out[5] += 0.7071067811865475*Ghat[4]*rdv2; 
  out[6] += 0.7071067811865475*Ghat[5]*rdv2; 
  out[7] += 0.7071067811865475*Ghat[6]*rdv2; 
  out[8] += 1.224744871391589*Ghat[1]*rdv2; 
  out[9] += 1.224744871391589*Ghat[2]*rdv2; 
  out[10] += 1.224744871391589*Ghat[3]*rdv2; 
  out[11] += 0.7071067811865475*Ghat[7]*rdv2; 
  out[12] += 1.224744871391589*Ghat[4]*rdv2; 
  out[13] += 1.224744871391589*Ghat[5]*rdv2; 
  out[14] += 1.224744871391589*Ghat[6]*rdv2; 
  out[15] += 1.224744871391589*Ghat[7]*rdv2; 

  } else { 

  alphaDrSurf[0] = 0.7071067811865475*(nuSum[0]*(2.0*w[2]-1.0*dxv[2])-2.0*nuUSum[0]); 
  alphaDrSurf[1] = 0.7071067811865475*(nuSum[1]*(2.0*w[2]-1.0*dxv[2])-2.0*nuUSum[1]); 
  alphaDrSurf[2] = 0.7071067811865475*(2.0*nuSum[2]*w[2]-2.0*nuUSum[2]-1.0*dxv[2]*nuSum[2]); 
  alphaDrSurf[4] = -0.7071067811865475*(2.0*nuUSum[3]+(dxv[2]-2.0*w[2])*nuSum[3]); 

  Ghat[0] = -0.25*(1.732050807568877*(alphaDrSurf[4]*fSkin[12]+alphaDrSurf[2]*fSkin[9]+alphaDrSurf[1]*fSkin[8])-1.0*alphaDrSurf[4]*fSkin[5]+1.732050807568877*alphaDrSurf[0]*fSkin[4]-1.0*(alphaDrSurf[2]*fSkin[2]+alphaDrSurf[1]*fSkin[1]+alphaDrSurf[0]*fSkin[0])); 
  Ghat[1] = -0.25*(1.732050807568877*(alphaDrSurf[2]*fSkin[12]+alphaDrSurf[4]*fSkin[9]+alphaDrSurf[0]*fSkin[8])-1.0*alphaDrSurf[2]*fSkin[5]+1.732050807568877*alphaDrSurf[1]*fSkin[4]-1.0*(fSkin[2]*alphaDrSurf[4]+alphaDrSurf[0]*fSkin[1]+fSkin[0]*alphaDrSurf[1])); 
  Ghat[2] = -0.25*(1.732050807568877*(alphaDrSurf[1]*fSkin[12]+alphaDrSurf[0]*fSkin[9]+alphaDrSurf[4]*fSkin[8])-1.0*alphaDrSurf[1]*fSkin[5]+1.732050807568877*alphaDrSurf[2]*fSkin[4]-1.0*(fSkin[1]*alphaDrSurf[4]+alphaDrSurf[0]*fSkin[2]+fSkin[0]*alphaDrSurf[2])); 
  Ghat[3] = -0.25*(1.732050807568877*(alphaDrSurf[4]*fSkin[15]+alphaDrSurf[2]*fSkin[14]+alphaDrSurf[1]*fSkin[13])-1.0*alphaDrSurf[4]*fSkin[11]+1.732050807568877*alphaDrSurf[0]*fSkin[10]-1.0*(alphaDrSurf[2]*fSkin[7]+alphaDrSurf[1]*fSkin[6]+alphaDrSurf[0]*fSkin[3])); 
  Ghat[4] = -0.25*(1.732050807568877*(alphaDrSurf[0]*fSkin[12]+alphaDrSurf[1]*fSkin[9]+alphaDrSurf[2]*fSkin[8])-1.0*alphaDrSurf[0]*fSkin[5]+1.732050807568877*alphaDrSurf[4]*fSkin[4]-1.0*(fSkin[0]*alphaDrSurf[4]+alphaDrSurf[1]*fSkin[2]+fSkin[1]*alphaDrSurf[2])); 
  Ghat[5] = -0.25*(1.732050807568877*(alphaDrSurf[2]*fSkin[15]+alphaDrSurf[4]*fSkin[14]+alphaDrSurf[0]*fSkin[13])-1.0*alphaDrSurf[2]*fSkin[11]+1.732050807568877*alphaDrSurf[1]*fSkin[10]-1.0*(alphaDrSurf[4]*fSkin[7]+alphaDrSurf[0]*fSkin[6]+alphaDrSurf[1]*fSkin[3])); 
  Ghat[6] = -0.25*(1.732050807568877*(alphaDrSurf[1]*fSkin[15]+alphaDrSurf[0]*fSkin[14]+alphaDrSurf[4]*fSkin[13])-1.0*alphaDrSurf[1]*fSkin[11]+1.732050807568877*alphaDrSurf[2]*fSkin[10]-1.0*(alphaDrSurf[0]*fSkin[7]+alphaDrSurf[4]*fSkin[6]+alphaDrSurf[2]*fSkin[3])); 
  Ghat[7] = -0.25*(1.732050807568877*(alphaDrSurf[0]*fSkin[15]+alphaDrSurf[1]*fSkin[14]+alphaDrSurf[2]*fSkin[13])-1.0*alphaDrSurf[0]*fSkin[11]+1.732050807568877*alphaDrSurf[4]*fSkin[10]-1.0*(alphaDrSurf[1]*fSkin[7]+alphaDrSurf[2]*fSkin[6]+fSkin[3]*alphaDrSurf[4])); 

  out[0] += -0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += -0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += -0.7071067811865475*Ghat[2]*rdv2; 
  out[3] += -0.7071067811865475*Ghat[3]*rdv2; 
  out[4] += 1.224744871391589*Ghat[0]*rdv2; 
  out[5] += -0.7071067811865475*Ghat[4]*rdv2; 
  out[6] += -0.7071067811865475*Ghat[5]*rdv2; 
  out[7] += -0.7071067811865475*Ghat[6]*rdv2; 
  out[8] += 1.224744871391589*Ghat[1]*rdv2; 
  out[9] += 1.224744871391589*Ghat[2]*rdv2; 
  out[10] += 1.224744871391589*Ghat[3]*rdv2; 
  out[11] += -0.7071067811865475*Ghat[7]*rdv2; 
  out[12] += 1.224744871391589*Ghat[4]*rdv2; 
  out[13] += 1.224744871391589*Ghat[5]*rdv2; 
  out[14] += 1.224744871391589*Ghat[6]*rdv2; 
  out[15] += 1.224744871391589*Ghat[7]*rdv2; 

  } 
} 
