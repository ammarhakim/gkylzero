#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_hyb_1x1v_p1_surfx2_eval_quad.h> 
#include <gkyl_basis_hyb_1x1v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH void vlasov_pkpm_surfvpar_1x1v_ser_p1(const double *w, const double *dxv, 
     const double *u_i, const double *div_p, const double *bvar, const double *rho_inv_b, 
     const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // u_i:      bulk flow velocity (ux, uy, uz).
  // div_p:     divergence of the pressure tensor.
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

  const double *div_p_x = &div_p[0]; 
  const double *div_p_y = &div_p[2]; 
  const double *div_p_z = &div_p[4]; 

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
  alphaSurf_l[0] = ((-1.224744871391589*bxbz[0]*uz[1])-1.224744871391589*bxby[0]*uy[1]-1.224744871391589*bxbx[0]*ux[1])*dx0*wvpar+(0.6123724356957944*bxbz[0]*uz[1]+0.6123724356957944*bxby[0]*uy[1]+0.6123724356957944*bxbx[0]*ux[1])*dvpar*dx0+0.7071067811865475*div_p_z[1]*rho_inv_bz[1]+0.7071067811865475*div_p_y[1]*rho_inv_by[1]+0.7071067811865475*div_p_x[1]*rho_inv_bx[1]+0.7071067811865475*div_p_z[0]*rho_inv_bz[0]+0.7071067811865475*div_p_y[0]*rho_inv_by[0]+0.7071067811865475*div_p_x[0]*rho_inv_bx[0]; 
  alphaSurf_l[1] = ((-1.224744871391589*bxbz[1]*uz[1])-1.224744871391589*bxby[1]*uy[1]-1.224744871391589*bxbx[1]*ux[1])*dx0*wvpar+(0.6123724356957944*bxbz[1]*uz[1]+0.6123724356957944*bxby[1]*uy[1]+0.6123724356957944*bxbx[1]*ux[1])*dvpar*dx0+0.7071067811865475*div_p_z[0]*rho_inv_bz[1]+0.7071067811865475*div_p_y[0]*rho_inv_by[1]+0.7071067811865475*div_p_x[0]*rho_inv_bx[1]+0.7071067811865475*rho_inv_bz[0]*div_p_z[1]+0.7071067811865475*rho_inv_by[0]*div_p_y[1]+0.7071067811865475*rho_inv_bx[0]*div_p_x[1]; 

  double alphaSurf_r[2] = {0.0}; 
  alphaSurf_r[0] = ((-1.224744871391589*bxbz[0]*uz[1])-1.224744871391589*bxby[0]*uy[1]-1.224744871391589*bxbx[0]*ux[1])*dx0*wvpar+((-0.6123724356957944*bxbz[0]*uz[1])-0.6123724356957944*bxby[0]*uy[1]-0.6123724356957944*bxbx[0]*ux[1])*dvpar*dx0+0.7071067811865475*div_p_z[1]*rho_inv_bz[1]+0.7071067811865475*div_p_y[1]*rho_inv_by[1]+0.7071067811865475*div_p_x[1]*rho_inv_bx[1]+0.7071067811865475*div_p_z[0]*rho_inv_bz[0]+0.7071067811865475*div_p_y[0]*rho_inv_by[0]+0.7071067811865475*div_p_x[0]*rho_inv_bx[0]; 
  alphaSurf_r[1] = ((-1.224744871391589*bxbz[1]*uz[1])-1.224744871391589*bxby[1]*uy[1]-1.224744871391589*bxbx[1]*ux[1])*dx0*wvpar+((-0.6123724356957944*bxbz[1]*uz[1])-0.6123724356957944*bxby[1]*uy[1]-0.6123724356957944*bxbx[1]*ux[1])*dvpar*dx0+0.7071067811865475*div_p_z[0]*rho_inv_bz[1]+0.7071067811865475*div_p_y[0]*rho_inv_by[1]+0.7071067811865475*div_p_x[0]*rho_inv_bx[1]+0.7071067811865475*rho_inv_bz[0]*div_p_z[1]+0.7071067811865475*rho_inv_by[0]*div_p_y[1]+0.7071067811865475*rho_inv_bx[0]*div_p_x[1]; 

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
