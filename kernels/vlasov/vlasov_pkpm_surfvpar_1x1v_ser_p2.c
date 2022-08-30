#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_ser_2x_p2_surfx2_eval_quad.h> 
#include <gkyl_basis_ser_2x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH void vlasov_pkpm_surfvpar_1x1v_ser_p2(const double *w, const double *dxv, 
     const double *u_i, const double *p_ij, const double *bvar,
     const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // u_i:      bulk flow velocity (ux, uy, uz).
  // p_ij:     pressure tensor (P_xx, P_xy, P_xz, P_yy, P_yz, P_zz).
  // bvar:      magnetic field unit vector (nine components; first three components, b_i, other six components, b_i b_j.) 
  // fl/fc/fr:  Input Distribution function in left/center/right cells.
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
  double alphaSurf_l[3] = {0.0}; 
  alphaSurf_l[0] = (-2.738612787525831*bxbz[1]*uz[2]*wvpar)-2.738612787525831*bxby[1]*uy[2]*wvpar-2.738612787525831*bxbx[1]*ux[2]*wvpar-1.224744871391589*bxbz[0]*uz[1]*wvpar-1.224744871391589*bxby[0]*uy[1]*wvpar-1.224744871391589*bxbx[0]*ux[1]*wvpar+1.369306393762915*bxbz[1]*uz[2]*dvpar+1.369306393762915*bxby[1]*uy[2]*dvpar+1.369306393762915*bxbx[1]*ux[2]*dvpar+0.6123724356957944*bxbz[0]*uz[1]*dvpar+0.6123724356957944*bxby[0]*uy[1]*dvpar+0.6123724356957944*bxbx[0]*ux[1]*dvpar+2.738612787525831*bz[1]*Pxz[2]+2.738612787525831*by[1]*Pxy[2]+2.738612787525831*bx[1]*Pxx[2]+1.224744871391589*bz[0]*Pxz[1]+1.224744871391589*by[0]*Pxy[1]+1.224744871391589*bx[0]*Pxx[1]; 
  alphaSurf_l[1] = (-2.449489742783178*bxbz[2]*uz[2]*wvpar)-2.738612787525831*bxbz[0]*uz[2]*wvpar-2.449489742783178*bxby[2]*uy[2]*wvpar-2.738612787525831*bxby[0]*uy[2]*wvpar-2.449489742783178*bxbx[2]*ux[2]*wvpar-2.738612787525831*bxbx[0]*ux[2]*wvpar-1.224744871391589*bxbz[1]*uz[1]*wvpar-1.224744871391589*bxby[1]*uy[1]*wvpar-1.224744871391589*bxbx[1]*ux[1]*wvpar+1.224744871391589*bxbz[2]*uz[2]*dvpar+1.369306393762915*bxbz[0]*uz[2]*dvpar+1.224744871391589*bxby[2]*uy[2]*dvpar+1.369306393762915*bxby[0]*uy[2]*dvpar+1.224744871391589*bxbx[2]*ux[2]*dvpar+1.369306393762915*bxbx[0]*ux[2]*dvpar+0.6123724356957944*bxbz[1]*uz[1]*dvpar+0.6123724356957944*bxby[1]*uy[1]*dvpar+0.6123724356957944*bxbx[1]*ux[1]*dvpar+2.449489742783178*Pxz[2]*bz[2]+2.449489742783178*Pxy[2]*by[2]+2.449489742783178*Pxx[2]*bx[2]+2.738612787525831*bz[0]*Pxz[2]+2.738612787525831*by[0]*Pxy[2]+2.738612787525831*bx[0]*Pxx[2]+1.224744871391589*Pxz[1]*bz[1]+1.224744871391589*Pxy[1]*by[1]+1.224744871391589*Pxx[1]*bx[1]; 
  alphaSurf_l[2] = (-2.449489742783178*bxbz[1]*uz[2]*wvpar)-2.449489742783178*bxby[1]*uy[2]*wvpar-2.449489742783178*bxbx[1]*ux[2]*wvpar-1.224744871391589*uz[1]*bxbz[2]*wvpar-1.224744871391589*uy[1]*bxby[2]*wvpar-1.224744871391589*ux[1]*bxbx[2]*wvpar+1.224744871391589*bxbz[1]*uz[2]*dvpar+1.224744871391589*bxby[1]*uy[2]*dvpar+1.224744871391589*bxbx[1]*ux[2]*dvpar+0.6123724356957944*uz[1]*bxbz[2]*dvpar+0.6123724356957944*uy[1]*bxby[2]*dvpar+0.6123724356957944*ux[1]*bxbx[2]*dvpar+1.224744871391589*Pxz[1]*bz[2]+1.224744871391589*Pxy[1]*by[2]+1.224744871391589*Pxx[1]*bx[2]+2.449489742783178*bz[1]*Pxz[2]+2.449489742783178*by[1]*Pxy[2]+2.449489742783178*bx[1]*Pxx[2]; 

  double alphaSurf_r[3] = {0.0}; 
  alphaSurf_r[0] = (-2.738612787525831*bxbz[1]*uz[2]*wvpar)-2.738612787525831*bxby[1]*uy[2]*wvpar-2.738612787525831*bxbx[1]*ux[2]*wvpar-1.224744871391589*bxbz[0]*uz[1]*wvpar-1.224744871391589*bxby[0]*uy[1]*wvpar-1.224744871391589*bxbx[0]*ux[1]*wvpar-1.369306393762915*bxbz[1]*uz[2]*dvpar-1.369306393762915*bxby[1]*uy[2]*dvpar-1.369306393762915*bxbx[1]*ux[2]*dvpar-0.6123724356957944*bxbz[0]*uz[1]*dvpar-0.6123724356957944*bxby[0]*uy[1]*dvpar-0.6123724356957944*bxbx[0]*ux[1]*dvpar+2.738612787525831*bz[1]*Pxz[2]+2.738612787525831*by[1]*Pxy[2]+2.738612787525831*bx[1]*Pxx[2]+1.224744871391589*bz[0]*Pxz[1]+1.224744871391589*by[0]*Pxy[1]+1.224744871391589*bx[0]*Pxx[1]; 
  alphaSurf_r[1] = (-2.449489742783178*bxbz[2]*uz[2]*wvpar)-2.738612787525831*bxbz[0]*uz[2]*wvpar-2.449489742783178*bxby[2]*uy[2]*wvpar-2.738612787525831*bxby[0]*uy[2]*wvpar-2.449489742783178*bxbx[2]*ux[2]*wvpar-2.738612787525831*bxbx[0]*ux[2]*wvpar-1.224744871391589*bxbz[1]*uz[1]*wvpar-1.224744871391589*bxby[1]*uy[1]*wvpar-1.224744871391589*bxbx[1]*ux[1]*wvpar-1.224744871391589*bxbz[2]*uz[2]*dvpar-1.369306393762915*bxbz[0]*uz[2]*dvpar-1.224744871391589*bxby[2]*uy[2]*dvpar-1.369306393762915*bxby[0]*uy[2]*dvpar-1.224744871391589*bxbx[2]*ux[2]*dvpar-1.369306393762915*bxbx[0]*ux[2]*dvpar-0.6123724356957944*bxbz[1]*uz[1]*dvpar-0.6123724356957944*bxby[1]*uy[1]*dvpar-0.6123724356957944*bxbx[1]*ux[1]*dvpar+2.449489742783178*Pxz[2]*bz[2]+2.449489742783178*Pxy[2]*by[2]+2.449489742783178*Pxx[2]*bx[2]+2.738612787525831*bz[0]*Pxz[2]+2.738612787525831*by[0]*Pxy[2]+2.738612787525831*bx[0]*Pxx[2]+1.224744871391589*Pxz[1]*bz[1]+1.224744871391589*Pxy[1]*by[1]+1.224744871391589*Pxx[1]*bx[1]; 
  alphaSurf_r[2] = (-2.449489742783178*bxbz[1]*uz[2]*wvpar)-2.449489742783178*bxby[1]*uy[2]*wvpar-2.449489742783178*bxbx[1]*ux[2]*wvpar-1.224744871391589*uz[1]*bxbz[2]*wvpar-1.224744871391589*uy[1]*bxby[2]*wvpar-1.224744871391589*ux[1]*bxbx[2]*wvpar-1.224744871391589*bxbz[1]*uz[2]*dvpar-1.224744871391589*bxby[1]*uy[2]*dvpar-1.224744871391589*bxbx[1]*ux[2]*dvpar-0.6123724356957944*uz[1]*bxbz[2]*dvpar-0.6123724356957944*uy[1]*bxby[2]*dvpar-0.6123724356957944*ux[1]*bxbx[2]*dvpar+1.224744871391589*Pxz[1]*bz[2]+1.224744871391589*Pxy[1]*by[2]+1.224744871391589*Pxx[1]*bx[2]+2.449489742783178*bz[1]*Pxz[2]+2.449489742783178*by[1]*Pxy[2]+2.449489742783178*bx[1]*Pxx[2]; 

  double fUpwindQuad_l[3] = {0.0};
  double fUpwindQuad_r[3] = {0.0};
  double fUpwind_l[3] = {0.0};;
  double fUpwind_r[3] = {0.0};
  double Ghat_l[3] = {0.0}; 
  double Ghat_r[3] = {0.0}; 

  if (0.6324555320336759*alphaSurf_l[2]-0.9486832980505137*alphaSurf_l[1]+0.7071067811865475*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[0] = ser_2x_p2_surfx2_eval_quad_node_0_r(fl); 
  } else { 
    fUpwindQuad_l[0] = ser_2x_p2_surfx2_eval_quad_node_0_l(fc); 
  } 
  if (0.6324555320336759*alphaSurf_r[2]-0.9486832980505137*alphaSurf_r[1]+0.7071067811865475*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[0] = ser_2x_p2_surfx2_eval_quad_node_0_r(fc); 
  } else { 
    fUpwindQuad_r[0] = ser_2x_p2_surfx2_eval_quad_node_0_l(fr); 
  } 
  if (0.7071067811865475*alphaSurf_l[0]-0.7905694150420947*alphaSurf_l[2] > 0) { 
    fUpwindQuad_l[1] = ser_2x_p2_surfx2_eval_quad_node_1_r(fl); 
  } else { 
    fUpwindQuad_l[1] = ser_2x_p2_surfx2_eval_quad_node_1_l(fc); 
  } 
  if (0.7071067811865475*alphaSurf_r[0]-0.7905694150420947*alphaSurf_r[2] > 0) { 
    fUpwindQuad_r[1] = ser_2x_p2_surfx2_eval_quad_node_1_r(fc); 
  } else { 
    fUpwindQuad_r[1] = ser_2x_p2_surfx2_eval_quad_node_1_l(fr); 
  } 
  if (0.6324555320336759*alphaSurf_l[2]+0.9486832980505137*alphaSurf_l[1]+0.7071067811865475*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[2] = ser_2x_p2_surfx2_eval_quad_node_2_r(fl); 
  } else { 
    fUpwindQuad_l[2] = ser_2x_p2_surfx2_eval_quad_node_2_l(fc); 
  } 
  if (0.6324555320336759*alphaSurf_r[2]+0.9486832980505137*alphaSurf_r[1]+0.7071067811865475*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[2] = ser_2x_p2_surfx2_eval_quad_node_2_r(fc); 
  } else { 
    fUpwindQuad_r[2] = ser_2x_p2_surfx2_eval_quad_node_2_l(fr); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_2x_p2_upwind_quad_to_modal(fUpwindQuad_l, fUpwind_l); 
  ser_2x_p2_upwind_quad_to_modal(fUpwindQuad_r, fUpwind_r); 

  Ghat_l[0] = 0.7071067811865475*alphaSurf_l[2]*fUpwind_l[2]+0.7071067811865475*alphaSurf_l[1]*fUpwind_l[1]+0.7071067811865475*alphaSurf_l[0]*fUpwind_l[0]; 
  Ghat_l[1] = 0.6324555320336759*alphaSurf_l[1]*fUpwind_l[2]+0.6324555320336759*fUpwind_l[1]*alphaSurf_l[2]+0.7071067811865475*alphaSurf_l[0]*fUpwind_l[1]+0.7071067811865475*fUpwind_l[0]*alphaSurf_l[1]; 
  Ghat_l[2] = 0.4517539514526256*alphaSurf_l[2]*fUpwind_l[2]+0.7071067811865475*alphaSurf_l[0]*fUpwind_l[2]+0.7071067811865475*fUpwind_l[0]*alphaSurf_l[2]+0.6324555320336759*alphaSurf_l[1]*fUpwind_l[1]; 

  Ghat_r[0] = 0.7071067811865475*alphaSurf_r[2]*fUpwind_r[2]+0.7071067811865475*alphaSurf_r[1]*fUpwind_r[1]+0.7071067811865475*alphaSurf_r[0]*fUpwind_r[0]; 
  Ghat_r[1] = 0.6324555320336759*alphaSurf_r[1]*fUpwind_r[2]+0.6324555320336759*fUpwind_r[1]*alphaSurf_r[2]+0.7071067811865475*alphaSurf_r[0]*fUpwind_r[1]+0.7071067811865475*fUpwind_r[0]*alphaSurf_r[1]; 
  Ghat_r[2] = 0.4517539514526256*alphaSurf_r[2]*fUpwind_r[2]+0.7071067811865475*alphaSurf_r[0]*fUpwind_r[2]+0.7071067811865475*fUpwind_r[0]*alphaSurf_r[2]+0.6324555320336759*alphaSurf_r[1]*fUpwind_r[1]; 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv1par; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv1par; 
  out[2] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv1par; 
  out[3] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv1par; 
  out[4] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dv1par; 
  out[5] += (1.58113883008419*Ghat_l[0]-1.58113883008419*Ghat_r[0])*dv1par; 
  out[6] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dv1par; 
  out[7] += (1.58113883008419*Ghat_l[1]-1.58113883008419*Ghat_r[1])*dv1par; 

} 
