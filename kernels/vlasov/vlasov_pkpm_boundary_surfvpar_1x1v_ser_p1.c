#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_hyb_1x1v_p1_surfx2_eval_quad.h> 
#include <gkyl_basis_hyb_1x1v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH void vlasov_pkpm_boundary_surfvpar_1x1v_ser_p1(const double *w, const double *dxv, 
     const double *u_i, const double *p_ij, const double *bvar, const double *rho_inv_b, 
     const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // u_i:      bulk flow velocity (ux, uy, uz).
  // p_ij:     pressure tensor (P_xx, P_xy, P_xz, P_yy, P_yz, P_zz).
  // bvar:      magnetic field unit vector (nine components; first three components, b_i, other six components, b_i b_j.) 
  // rho_inv_b: b_i/rho (for pressure force 1/rho * b . div(P)).
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:       Incremented distribution function in center cell.
  const double dx0 = 2.0/dxv[0]; 
  const double dv1par = 2.0/dxv[1]; 
  const double dvpar = dxv[1], wvpar = w[1]; 
  const double *ux = &u_i[0]; 
  const double *uy = &u_i[2]; 
  const double *uz = &u_i[4]; 

  const double *Pxx = &p_ij[0]; 
  const double *Pxy = &p_ij[2]; 
  const double *Pxz = &p_ij[4]; 
  const double *Pyy = &p_ij[6]; 
  const double *Pyz = &p_ij[8]; 
  const double *Pzz = &p_ij[10]; 

  const double *bx = &bvar[0]; 
  const double *by = &bvar[2]; 
  const double *bz = &bvar[4]; 
  const double *bxbx = &bvar[6]; 
  const double *bxby = &bvar[8]; 
  const double *bxbz = &bvar[10]; 
  const double *byby = &bvar[12]; 
  const double *bybz = &bvar[14]; 
  const double *bzbz = &bvar[16]; 

  const double *rho_inv_bx = &rho_inv_b[0]; 
  const double *rho_inv_by = &rho_inv_b[2]; 
  const double *rho_inv_bz = &rho_inv_b[4]; 

  double alphaSurf[2] = {0.0}; 
  double fUpwindQuad[2] = {0.0};
  double fUpwind[2] = {0.0};;
  double Ghat[2] = {0.0}; 

  if (edge == -1) { 

  alphaSurf[0] = ((-1.224744871391589*bxbz[0]*uz[1])-1.224744871391589*bxby[0]*uy[1]-1.224744871391589*bxbx[0]*ux[1])*dx0*wvpar+((-0.6123724356957944*bxbz[0]*uz[1])-0.6123724356957944*bxby[0]*uy[1]-0.6123724356957944*bxbx[0]*ux[1])*dvpar*dx0+(1.224744871391589*rho_inv_bz[0]*Pxz[1]+1.224744871391589*rho_inv_by[0]*Pxy[1]+1.224744871391589*rho_inv_bx[0]*Pxx[1])*dx0; 
  alphaSurf[1] = ((-1.224744871391589*bxbz[1]*uz[1])-1.224744871391589*bxby[1]*uy[1]-1.224744871391589*bxbx[1]*ux[1])*dx0*wvpar+((-0.6123724356957944*bxbz[1]*uz[1])-0.6123724356957944*bxby[1]*uy[1]-0.6123724356957944*bxbx[1]*ux[1])*dvpar*dx0+(1.224744871391589*Pxz[1]*rho_inv_bz[1]+1.224744871391589*Pxy[1]*rho_inv_by[1]+1.224744871391589*Pxx[1]*rho_inv_bx[1])*dx0; 

  if (0.7071067811865475*alphaSurf[0]-0.7071067811865475*alphaSurf[1] > 0) { 
    fUpwindQuad[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_r(fSkin); 
  } else { 
    fUpwindQuad[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_l(fEdge); 
  } 
  if (0.7071067811865475*(alphaSurf[1]+alphaSurf[0]) > 0) { 
    fUpwindQuad[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_r(fSkin); 
  } else { 
    fUpwindQuad[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_l(fEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_1x1v_p1_vdir_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.7071067811865475*alphaSurf[1]*fUpwind[1]+0.7071067811865475*alphaSurf[0]*fUpwind[0]; 
  Ghat[1] = 0.7071067811865475*alphaSurf[0]*fUpwind[1]+0.7071067811865475*fUpwind[0]*alphaSurf[1]; 

  out[0] += -0.7071067811865475*Ghat[0]*dv1par; 
  out[1] += -0.7071067811865475*Ghat[1]*dv1par; 
  out[2] += -1.224744871391589*Ghat[0]*dv1par; 
  out[3] += -1.224744871391589*Ghat[1]*dv1par; 
  out[4] += -1.58113883008419*Ghat[0]*dv1par; 
  out[5] += -1.58113883008419*Ghat[1]*dv1par; 

  } else { 

  alphaSurf[0] = ((-1.224744871391589*bxbz[0]*uz[1])-1.224744871391589*bxby[0]*uy[1]-1.224744871391589*bxbx[0]*ux[1])*dx0*wvpar+(0.6123724356957944*bxbz[0]*uz[1]+0.6123724356957944*bxby[0]*uy[1]+0.6123724356957944*bxbx[0]*ux[1])*dvpar*dx0+(1.224744871391589*rho_inv_bz[0]*Pxz[1]+1.224744871391589*rho_inv_by[0]*Pxy[1]+1.224744871391589*rho_inv_bx[0]*Pxx[1])*dx0; 
  alphaSurf[1] = ((-1.224744871391589*bxbz[1]*uz[1])-1.224744871391589*bxby[1]*uy[1]-1.224744871391589*bxbx[1]*ux[1])*dx0*wvpar+(0.6123724356957944*bxbz[1]*uz[1]+0.6123724356957944*bxby[1]*uy[1]+0.6123724356957944*bxbx[1]*ux[1])*dvpar*dx0+(1.224744871391589*Pxz[1]*rho_inv_bz[1]+1.224744871391589*Pxy[1]*rho_inv_by[1]+1.224744871391589*Pxx[1]*rho_inv_bx[1])*dx0; 

  if (0.7071067811865475*alphaSurf[0]-0.7071067811865475*alphaSurf[1] > 0) { 
    fUpwindQuad[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_r(fEdge); 
  } else { 
    fUpwindQuad[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_l(fSkin); 
  } 
  if (0.7071067811865475*(alphaSurf[1]+alphaSurf[0]) > 0) { 
    fUpwindQuad[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_r(fEdge); 
  } else { 
    fUpwindQuad[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_l(fSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_1x1v_p1_vdir_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.7071067811865475*alphaSurf[1]*fUpwind[1]+0.7071067811865475*alphaSurf[0]*fUpwind[0]; 
  Ghat[1] = 0.7071067811865475*alphaSurf[0]*fUpwind[1]+0.7071067811865475*fUpwind[0]*alphaSurf[1]; 

  out[0] += 0.7071067811865475*Ghat[0]*dv1par; 
  out[1] += 0.7071067811865475*Ghat[1]*dv1par; 
  out[2] += -1.224744871391589*Ghat[0]*dv1par; 
  out[3] += -1.224744871391589*Ghat[1]*dv1par; 
  out[4] += 1.58113883008419*Ghat[0]*dv1par; 
  out[5] += 1.58113883008419*Ghat[1]*dv1par; 

  } 
} 
