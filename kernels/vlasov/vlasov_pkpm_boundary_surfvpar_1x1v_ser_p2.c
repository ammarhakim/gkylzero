#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_ser_2x_p2_surfx2_eval_quad.h> 
#include <gkyl_basis_ser_2x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH void vlasov_pkpm_boundary_surfvpar_1x1v_ser_p2(const double *w, const double *dxv, 
     const double *u_i, const double *p_ij, const double *bvar,
     const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // u_i:      bulk flow velocity (ux, uy, uz).
  // p_ij:     pressure tensor (P_xx, P_xy, P_xz, P_yy, P_yz, P_zz).
  // bvar:      magnetic field unit vector (nine components; first three components, b_i, other six components, b_i b_j.) 
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:       Incremented distribution function in center cell.
  const double dv1par = 2.0/dxv[1]; 
  const double dvpar = dxv[1], wvpar = w[1]; 
  const double *ux = &u_i[0]; 
  const double *uy = &u_i[3]; 
  const double *uz = &u_i[6]; 

  const double *Pxx = &p_ij[0]; 
  const double *Pxy = &p_ij[3]; 
  const double *Pxz = &p_ij[6]; 
  const double *Pyy = &p_ij[9]; 
  const double *Pyz = &p_ij[12]; 
  const double *Pzz = &p_ij[15]; 

  const double *bx = &bvar[0]; 
  const double *by = &bvar[3]; 
  const double *bz = &bvar[6]; 
  const double *bxbx = &bvar[9]; 
  const double *bxby = &bvar[12]; 
  const double *bxbz = &bvar[15]; 
  const double *byby = &bvar[18]; 
  const double *bybz = &bvar[21]; 
  const double *bzbz = &bvar[24]; 
  double alphaSurf[3] = {0.0}; 
  double fUpwindQuad[3] = {0.0};
  double fUpwind[3] = {0.0};;
  double Ghat[3] = {0.0}; 

  if (edge == -1) { 

  alphaSurf[0] = 0.3535533905932737*((7.745966692414834*(bxbz[1]*uz[2]+bxby[1]*uy[2]+bxbx[1]*ux[2])+3.464101615137754*(bxbz[0]*uz[1]+bxby[0]*uy[1]+bxbx[0]*ux[1]))*wvpar+(3.872983346207417*(bxbz[1]*uz[2]+bxby[1]*uy[2]+bxbx[1]*ux[2])+1.732050807568877*(bxbz[0]*uz[1]+bxby[0]*uy[1]+bxbx[0]*ux[1]))*dvpar+7.745966692414834*(bz[1]*Pxz[2]+by[1]*Pxy[2]+bx[1]*Pxx[2])+3.464101615137754*(bz[0]*Pxz[1]+by[0]*Pxy[1]+bx[0]*Pxx[1])); 
  alphaSurf[1] = 0.3535533905932737*(((6.928203230275509*bxbz[2]+7.745966692414834*bxbz[0])*uz[2]+(6.928203230275509*bxby[2]+7.745966692414834*bxby[0])*uy[2]+(6.928203230275509*bxbx[2]+7.745966692414834*bxbx[0])*ux[2]+3.464101615137754*(bxbz[1]*uz[1]+bxby[1]*uy[1]+bxbx[1]*ux[1]))*wvpar+((3.464101615137754*bxbz[2]+3.872983346207417*bxbz[0])*uz[2]+(3.464101615137754*bxby[2]+3.872983346207417*bxby[0])*uy[2]+(3.464101615137754*bxbx[2]+3.872983346207417*bxbx[0])*ux[2]+1.732050807568877*(bxbz[1]*uz[1]+bxby[1]*uy[1]+bxbx[1]*ux[1]))*dvpar+6.928203230275509*(Pxz[2]*bz[2]+Pxy[2]*by[2]+Pxx[2]*bx[2])+7.745966692414834*(bz[0]*Pxz[2]+by[0]*Pxy[2]+bx[0]*Pxx[2])+3.464101615137754*(Pxz[1]*bz[1]+Pxy[1]*by[1]+Pxx[1]*bx[1])); 
  alphaSurf[2] = 0.3535533905932737*((6.928203230275509*(bxbz[1]*uz[2]+bxby[1]*uy[2]+bxbx[1]*ux[2])+3.464101615137754*(uz[1]*bxbz[2]+uy[1]*bxby[2]+ux[1]*bxbx[2]))*wvpar+(3.464101615137754*(bxbz[1]*uz[2]+bxby[1]*uy[2]+bxbx[1]*ux[2])+1.732050807568877*(uz[1]*bxbz[2]+uy[1]*bxby[2]+ux[1]*bxbx[2]))*dvpar+3.464101615137754*(Pxz[1]*bz[2]+Pxy[1]*by[2]+Pxx[1]*bx[2])+6.928203230275509*(bz[1]*Pxz[2]+by[1]*Pxy[2]+bx[1]*Pxx[2])); 

  if (0.6324555320336759*alphaSurf[2]-0.9486832980505137*alphaSurf[1]+0.7071067811865475*alphaSurf[0] > 0) { 
    fUpwindQuad[0] = ser_2x_p2_surfx2_eval_quad_node_0_r(fSkin); 
  } else { 
    fUpwindQuad[0] = ser_2x_p2_surfx2_eval_quad_node_0_l(fEdge); 
  } 
  if (0.7071067811865475*alphaSurf[0]-0.7905694150420947*alphaSurf[2] > 0) { 
    fUpwindQuad[1] = ser_2x_p2_surfx2_eval_quad_node_1_r(fSkin); 
  } else { 
    fUpwindQuad[1] = ser_2x_p2_surfx2_eval_quad_node_1_l(fEdge); 
  } 
  if (0.6324555320336759*alphaSurf[2]+0.9486832980505137*alphaSurf[1]+0.7071067811865475*alphaSurf[0] > 0) { 
    fUpwindQuad[2] = ser_2x_p2_surfx2_eval_quad_node_2_r(fSkin); 
  } else { 
    fUpwindQuad[2] = ser_2x_p2_surfx2_eval_quad_node_2_l(fEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_2x_p2_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.7071067811865475*alphaSurf[2]*fUpwind[2]+0.7071067811865475*alphaSurf[1]*fUpwind[1]+0.7071067811865475*alphaSurf[0]*fUpwind[0]; 
  Ghat[1] = 0.6324555320336759*alphaSurf[1]*fUpwind[2]+0.6324555320336759*fUpwind[1]*alphaSurf[2]+0.7071067811865475*alphaSurf[0]*fUpwind[1]+0.7071067811865475*fUpwind[0]*alphaSurf[1]; 
  Ghat[2] = 0.4517539514526256*alphaSurf[2]*fUpwind[2]+0.7071067811865475*alphaSurf[0]*fUpwind[2]+0.7071067811865475*fUpwind[0]*alphaSurf[2]+0.6324555320336759*alphaSurf[1]*fUpwind[1]; 

  out[0] += -0.7071067811865475*Ghat[0]*dv1par; 
  out[1] += -0.7071067811865475*Ghat[1]*dv1par; 
  out[2] += -1.224744871391589*Ghat[0]*dv1par; 
  out[3] += -1.224744871391589*Ghat[1]*dv1par; 
  out[4] += -0.7071067811865475*Ghat[2]*dv1par; 
  out[5] += -1.58113883008419*Ghat[0]*dv1par; 
  out[6] += -1.224744871391589*Ghat[2]*dv1par; 
  out[7] += -1.58113883008419*Ghat[1]*dv1par; 

  } else { 

  alphaSurf[0] = 0.3535533905932737*((7.745966692414834*(bxbz[1]*uz[2]+bxby[1]*uy[2]+bxbx[1]*ux[2])+3.464101615137754*(bxbz[0]*uz[1]+bxby[0]*uy[1]+bxbx[0]*ux[1]))*wvpar+((-3.872983346207417*(bxbz[1]*uz[2]+bxby[1]*uy[2]+bxbx[1]*ux[2]))-1.732050807568877*(bxbz[0]*uz[1]+bxby[0]*uy[1]+bxbx[0]*ux[1]))*dvpar+7.745966692414834*(bz[1]*Pxz[2]+by[1]*Pxy[2]+bx[1]*Pxx[2])+3.464101615137754*(bz[0]*Pxz[1]+by[0]*Pxy[1]+bx[0]*Pxx[1])); 
  alphaSurf[1] = 0.3535533905932737*(((6.928203230275509*bxbz[2]+7.745966692414834*bxbz[0])*uz[2]+(6.928203230275509*bxby[2]+7.745966692414834*bxby[0])*uy[2]+(6.928203230275509*bxbx[2]+7.745966692414834*bxbx[0])*ux[2]+3.464101615137754*(bxbz[1]*uz[1]+bxby[1]*uy[1]+bxbx[1]*ux[1]))*wvpar+(((-3.464101615137754*bxbz[2])-3.872983346207417*bxbz[0])*uz[2]+((-3.464101615137754*bxby[2])-3.872983346207417*bxby[0])*uy[2]+((-3.464101615137754*bxbx[2])-3.872983346207417*bxbx[0])*ux[2]-1.732050807568877*(bxbz[1]*uz[1]+bxby[1]*uy[1]+bxbx[1]*ux[1]))*dvpar+6.928203230275509*(Pxz[2]*bz[2]+Pxy[2]*by[2]+Pxx[2]*bx[2])+7.745966692414834*(bz[0]*Pxz[2]+by[0]*Pxy[2]+bx[0]*Pxx[2])+3.464101615137754*(Pxz[1]*bz[1]+Pxy[1]*by[1]+Pxx[1]*bx[1])); 
  alphaSurf[2] = 0.3535533905932737*((6.928203230275509*(bxbz[1]*uz[2]+bxby[1]*uy[2]+bxbx[1]*ux[2])+3.464101615137754*(uz[1]*bxbz[2]+uy[1]*bxby[2]+ux[1]*bxbx[2]))*wvpar+((-3.464101615137754*(bxbz[1]*uz[2]+bxby[1]*uy[2]+bxbx[1]*ux[2]))-1.732050807568877*(uz[1]*bxbz[2]+uy[1]*bxby[2]+ux[1]*bxbx[2]))*dvpar+3.464101615137754*(Pxz[1]*bz[2]+Pxy[1]*by[2]+Pxx[1]*bx[2])+6.928203230275509*(bz[1]*Pxz[2]+by[1]*Pxy[2]+bx[1]*Pxx[2])); 

  if (0.6324555320336759*alphaSurf[2]-0.9486832980505137*alphaSurf[1]+0.7071067811865475*alphaSurf[0] > 0) { 
    fUpwindQuad[0] = ser_2x_p2_surfx2_eval_quad_node_0_r(fEdge); 
  } else { 
    fUpwindQuad[0] = ser_2x_p2_surfx2_eval_quad_node_0_l(fSkin); 
  } 
  if (0.7071067811865475*alphaSurf[0]-0.7905694150420947*alphaSurf[2] > 0) { 
    fUpwindQuad[1] = ser_2x_p2_surfx2_eval_quad_node_1_r(fEdge); 
  } else { 
    fUpwindQuad[1] = ser_2x_p2_surfx2_eval_quad_node_1_l(fSkin); 
  } 
  if (0.6324555320336759*alphaSurf[2]+0.9486832980505137*alphaSurf[1]+0.7071067811865475*alphaSurf[0] > 0) { 
    fUpwindQuad[2] = ser_2x_p2_surfx2_eval_quad_node_2_r(fEdge); 
  } else { 
    fUpwindQuad[2] = ser_2x_p2_surfx2_eval_quad_node_2_l(fSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_2x_p2_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.7071067811865475*alphaSurf[2]*fUpwind[2]+0.7071067811865475*alphaSurf[1]*fUpwind[1]+0.7071067811865475*alphaSurf[0]*fUpwind[0]; 
  Ghat[1] = 0.6324555320336759*alphaSurf[1]*fUpwind[2]+0.6324555320336759*fUpwind[1]*alphaSurf[2]+0.7071067811865475*alphaSurf[0]*fUpwind[1]+0.7071067811865475*fUpwind[0]*alphaSurf[1]; 
  Ghat[2] = 0.4517539514526256*alphaSurf[2]*fUpwind[2]+0.7071067811865475*alphaSurf[0]*fUpwind[2]+0.7071067811865475*fUpwind[0]*alphaSurf[2]+0.6324555320336759*alphaSurf[1]*fUpwind[1]; 

  out[0] += 0.7071067811865475*Ghat[0]*dv1par; 
  out[1] += 0.7071067811865475*Ghat[1]*dv1par; 
  out[2] += -1.224744871391589*Ghat[0]*dv1par; 
  out[3] += -1.224744871391589*Ghat[1]*dv1par; 
  out[4] += 0.7071067811865475*Ghat[2]*dv1par; 
  out[5] += 1.58113883008419*Ghat[0]*dv1par; 
  out[6] += -1.224744871391589*Ghat[2]*dv1par; 
  out[7] += 1.58113883008419*Ghat[1]*dv1par; 

  } 
} 
