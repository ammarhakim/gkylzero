#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_hyb_1x1v_p1_surfx2_eval_quad.h> 
#include <gkyl_basis_hyb_1x1v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH void vlasov_pkpm_surfvpar_1x1v_ser_p1(const double *w, const double *dxv, 
     const double *u_i, const double *p_ij, const double *bvar, const double *rho_inv_b, 
     const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // u_i:      bulk flow velocity (ux, uy, uz).
  // p_ij:     pressure tensor (P_xx, P_xy, P_xz, P_yy, P_yz, P_zz).
  // bvar:      magnetic field unit vector (nine components; first three components, b_i, other six components, b_i b_j.) 
  // rho_inv_b: b_i/rho (for pressure force 1/rho * b . div(P)).
  // fl/fc/fr:  Input Distribution function in left/center/right cells.
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

  double alphaSurf_l[2] = {0.0}; 
  alphaSurf_l[0] = ((-1.224744871391589*bxbz[0]*uz[1])-1.224744871391589*bxby[0]*uy[1]-1.224744871391589*bxbx[0]*ux[1])*dx0*wvpar+(0.6123724356957944*bxbz[0]*uz[1]+0.6123724356957944*bxby[0]*uy[1]+0.6123724356957944*bxbx[0]*ux[1])*dvpar*dx0+(1.224744871391589*rho_inv_bz[0]*Pxz[1]+1.224744871391589*rho_inv_by[0]*Pxy[1]+1.224744871391589*rho_inv_bx[0]*Pxx[1])*dx0; 
  alphaSurf_l[1] = ((-1.224744871391589*bxbz[1]*uz[1])-1.224744871391589*bxby[1]*uy[1]-1.224744871391589*bxbx[1]*ux[1])*dx0*wvpar+(0.6123724356957944*bxbz[1]*uz[1]+0.6123724356957944*bxby[1]*uy[1]+0.6123724356957944*bxbx[1]*ux[1])*dvpar*dx0+(1.224744871391589*Pxz[1]*rho_inv_bz[1]+1.224744871391589*Pxy[1]*rho_inv_by[1]+1.224744871391589*Pxx[1]*rho_inv_bx[1])*dx0; 

  double alphaSurf_r[2] = {0.0}; 
  alphaSurf_r[0] = ((-1.224744871391589*bxbz[0]*uz[1])-1.224744871391589*bxby[0]*uy[1]-1.224744871391589*bxbx[0]*ux[1])*dx0*wvpar+((-0.6123724356957944*bxbz[0]*uz[1])-0.6123724356957944*bxby[0]*uy[1]-0.6123724356957944*bxbx[0]*ux[1])*dvpar*dx0+(1.224744871391589*rho_inv_bz[0]*Pxz[1]+1.224744871391589*rho_inv_by[0]*Pxy[1]+1.224744871391589*rho_inv_bx[0]*Pxx[1])*dx0; 
  alphaSurf_r[1] = ((-1.224744871391589*bxbz[1]*uz[1])-1.224744871391589*bxby[1]*uy[1]-1.224744871391589*bxbx[1]*ux[1])*dx0*wvpar+((-0.6123724356957944*bxbz[1]*uz[1])-0.6123724356957944*bxby[1]*uy[1]-0.6123724356957944*bxbx[1]*ux[1])*dvpar*dx0+(1.224744871391589*Pxz[1]*rho_inv_bz[1]+1.224744871391589*Pxy[1]*rho_inv_by[1]+1.224744871391589*Pxx[1]*rho_inv_bx[1])*dx0; 

  double fUpwindQuad_l[2] = {0.0};
  double fUpwindQuad_r[2] = {0.0};
  double fUpwind_l[2] = {0.0};
  double fUpwind_r[2] = {0.0};
  double Ghat_l[2] = {0.0}; 
  double Ghat_r[2] = {0.0}; 

  if (0.7071067811865475*alphaSurf_l[0]-0.7071067811865475*alphaSurf_l[1] > 0) { 
    fUpwindQuad_l[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_r(fl); 
  } else { 
    fUpwindQuad_l[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_l(fc); 
  } 
  if (0.7071067811865475*alphaSurf_r[0]-0.7071067811865475*alphaSurf_r[1] > 0) { 
    fUpwindQuad_r[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_r(fc); 
  } else { 
    fUpwindQuad_r[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_l(fr); 
  } 
  if (0.7071067811865475*(alphaSurf_l[1]+alphaSurf_l[0]) > 0) { 
    fUpwindQuad_l[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_r(fl); 
  } else { 
    fUpwindQuad_l[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_l(fc); 
  } 
  if (0.7071067811865475*(alphaSurf_r[1]+alphaSurf_r[0]) > 0) { 
    fUpwindQuad_r[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_r(fc); 
  } else { 
    fUpwindQuad_r[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_l(fr); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_1x1v_p1_vdir_upwind_quad_to_modal(fUpwindQuad_l, fUpwind_l); 
  hyb_1x1v_p1_vdir_upwind_quad_to_modal(fUpwindQuad_r, fUpwind_r); 

  Ghat_l[0] = 0.7071067811865475*alphaSurf_l[1]*fUpwind_l[1]+0.7071067811865475*alphaSurf_l[0]*fUpwind_l[0]; 
  Ghat_l[1] = 0.7071067811865475*alphaSurf_l[0]*fUpwind_l[1]+0.7071067811865475*fUpwind_l[0]*alphaSurf_l[1]; 

  Ghat_r[0] = 0.7071067811865475*alphaSurf_r[1]*fUpwind_r[1]+0.7071067811865475*alphaSurf_r[0]*fUpwind_r[0]; 
  Ghat_r[1] = 0.7071067811865475*alphaSurf_r[0]*fUpwind_r[1]+0.7071067811865475*fUpwind_r[0]*alphaSurf_r[1]; 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv1par; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv1par; 
  out[2] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv1par; 
  out[3] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv1par; 
  out[4] += (1.58113883008419*Ghat_l[0]-1.58113883008419*Ghat_r[0])*dv1par; 
  out[5] += (1.58113883008419*Ghat_l[1]-1.58113883008419*Ghat_r[1])*dv1par; 

} 
