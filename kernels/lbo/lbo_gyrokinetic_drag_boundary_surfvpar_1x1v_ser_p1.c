#include <gkyl_lbo_gyrokinetic_kernels.h> 
#include <gkyl_basis_ser_2x_p1_surfx2_quad.h> 
#include <gkyl_basis_ser_2x_p1_upwind.h> 
GKYL_CU_DH void lbo_gyrokinetic_drag_boundary_surfvpar_1x1v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[2]:     cell-center coordinates. 
  // dxv[2]:   cell spacing. 
  // m_:        species mass.
  // bmag_inv:  1/(magnetic field magnitude). 
  // nuSum:     collisionalities added (self and cross species collisionalities). 
  // nuUSum[2]:sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum: sum of thermal speeds squared time their respective collisionalities. 
  // fSkin/Edge:    Distribution function in cells 
  // out:           Incremented distribution function in cell 
  double rdv2 = 2.0/dxv[1]; 


  double alphaDrSurf[2] = {0.0}; 
  double fUpwindQuad[2] = {0.0};
  double fUpwind[2] = {0.0};;
  double drag_incr[2] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = 0.5*(nuSum[0]*(2.0*w[1]+dxv[1])-2.0*nuUSum[0]); 
  alphaDrSurf[1] = 0.5*(2.0*nuSum[1]*w[1]-2.0*nuUSum[1]+dxv[1]*nuSum[1]); 

  if (alphaDrSurf[0]-alphaDrSurf[1] < 0) { 
    fUpwindQuad[0] = ser_2x_p1_surfx2_quad_0_r(fSkin); 
  } else { 
    fUpwindQuad[0] = ser_2x_p1_surfx2_quad_0_l(fEdge); 
  } 
  if (alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[1] = ser_2x_p1_surfx2_quad_1_r(fSkin); 
  } else { 
    fUpwindQuad[1] = ser_2x_p1_surfx2_quad_1_l(fEdge); 
  } 

  // Project nodal basis back onto modal basis. 
  ser_2x_p1_upwind(fUpwindQuad, fUpwind); 

  drag_incr[0] = 0.7071067811865475*alphaDrSurf[1]*fUpwind[1]+0.7071067811865475*alphaDrSurf[0]*fUpwind[0]; 
  drag_incr[1] = 0.7071067811865475*alphaDrSurf[0]*fUpwind[1]+0.7071067811865475*fUpwind[0]*alphaDrSurf[1]; 

  out[0] += 0.7071067811865475*drag_incr[0]*rdv2; 
  out[1] += 0.7071067811865475*drag_incr[1]*rdv2; 
  out[2] += 1.224744871391589*drag_incr[0]*rdv2; 
  out[3] += 1.224744871391589*drag_incr[1]*rdv2; 

  } else { 

  alphaDrSurf[0] = 0.5*(nuSum[0]*(2.0*w[1]-1.0*dxv[1])-2.0*nuUSum[0]); 
  alphaDrSurf[1] = 0.5*(2.0*nuSum[1]*w[1]-2.0*nuUSum[1]-1.0*dxv[1]*nuSum[1]); 

  if (alphaDrSurf[0]-alphaDrSurf[1] < 0) { 
    fUpwindQuad[0] = ser_2x_p1_surfx2_quad_0_r(fEdge); 
  } else { 
    fUpwindQuad[0] = ser_2x_p1_surfx2_quad_0_l(fSkin); 
  } 
  if (alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[1] = ser_2x_p1_surfx2_quad_1_r(fEdge); 
  } else { 
    fUpwindQuad[1] = ser_2x_p1_surfx2_quad_1_l(fSkin); 
  } 

  // Project nodal basis back onto modal basis. 
  ser_2x_p1_upwind(fUpwindQuad, fUpwind); 

  drag_incr[0] = 0.7071067811865475*alphaDrSurf[1]*fUpwind[1]+0.7071067811865475*alphaDrSurf[0]*fUpwind[0]; 
  drag_incr[1] = 0.7071067811865475*alphaDrSurf[0]*fUpwind[1]+0.7071067811865475*fUpwind[0]*alphaDrSurf[1]; 

  out[0] += -0.7071067811865475*drag_incr[0]*rdv2; 
  out[1] += -0.7071067811865475*drag_incr[1]*rdv2; 
  out[2] += 1.224744871391589*drag_incr[0]*rdv2; 
  out[3] += 1.224744871391589*drag_incr[1]*rdv2; 

  } 
} 
