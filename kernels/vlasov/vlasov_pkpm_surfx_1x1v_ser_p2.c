#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_ser_2x_p2_surfx1_eval_quad.h> 
#include <gkyl_basis_ser_2x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH void vlasov_pkpm_surfx_1x1v_ser_p2(const double *w, const double *dxv, 
     const double *u_i, const double *bvar, 
     const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // u_i: Input bulk velocity (ux,uy,uz) in cell being updated (ASSUMED TO BE CONTINUOUS).
  // bvar: Input magnetic field unit vector in cell being updated (ASSUMED TO BE CONTINUOUS).
  // fl/fc/fr:  Input Distribution function in left/center/right cells.
  // out:       Incremented distribution function in center cell.
  const double dx1 = 2.0/dxv[0]; 
  const double dvpar = dxv[1], wvpar = w[1]; 
  const double *uc = &u_i[0]; 
  const double *bc = &bvar[0]; 
  double alphaSurf_l[3] = {0.0}; 
  alphaSurf_l[0] = 2.23606797749979*bc[2]*wvpar-1.732050807568877*bc[1]*wvpar+bc[0]*wvpar+2.23606797749979*uc[2]-1.732050807568877*uc[1]+uc[0]; 
  alphaSurf_l[1] = 0.6454972243679029*bc[2]*dvpar-0.5*bc[1]*dvpar+0.2886751345948129*bc[0]*dvpar; 

  double alphaSurf_r[3] = {0.0}; 
  alphaSurf_r[0] = 2.23606797749979*bc[2]*wvpar+1.732050807568877*bc[1]*wvpar+bc[0]*wvpar+2.23606797749979*uc[2]+1.732050807568877*uc[1]+uc[0]; 
  alphaSurf_r[1] = 0.6454972243679029*bc[2]*dvpar+0.5*bc[1]*dvpar+0.2886751345948129*bc[0]*dvpar; 

  double fUpwindQuad_l[3] = {0.0};
  double fUpwindQuad_r[3] = {0.0};
  double fUpwind_l[3] = {0.0};
  double fUpwind_r[3] = {0.0};
  double Ghat_l[3] = {0.0}; 
  double Ghat_r[3] = {0.0}; 

  if (0.7071067811865475*alphaSurf_l[0]-0.9486832980505137*alphaSurf_l[1] > 0) { 
    fUpwindQuad_l[0] = ser_2x_p2_surfx1_eval_quad_node_0_r(fl); 
  } else { 
    fUpwindQuad_l[0] = ser_2x_p2_surfx1_eval_quad_node_0_l(fc); 
  } 
  if (0.7071067811865475*alphaSurf_r[0]-0.9486832980505137*alphaSurf_r[1] > 0) { 
    fUpwindQuad_r[0] = ser_2x_p2_surfx1_eval_quad_node_0_r(fc); 
  } else { 
    fUpwindQuad_r[0] = ser_2x_p2_surfx1_eval_quad_node_0_l(fr); 
  } 
  if (0.7071067811865475*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[1] = ser_2x_p2_surfx1_eval_quad_node_1_r(fl); 
  } else { 
    fUpwindQuad_l[1] = ser_2x_p2_surfx1_eval_quad_node_1_l(fc); 
  } 
  if (0.7071067811865475*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[1] = ser_2x_p2_surfx1_eval_quad_node_1_r(fc); 
  } else { 
    fUpwindQuad_r[1] = ser_2x_p2_surfx1_eval_quad_node_1_l(fr); 
  } 
  if (0.9486832980505137*alphaSurf_l[1]+0.7071067811865475*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[2] = ser_2x_p2_surfx1_eval_quad_node_2_r(fl); 
  } else { 
    fUpwindQuad_l[2] = ser_2x_p2_surfx1_eval_quad_node_2_l(fc); 
  } 
  if (0.9486832980505137*alphaSurf_r[1]+0.7071067811865475*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[2] = ser_2x_p2_surfx1_eval_quad_node_2_r(fc); 
  } else { 
    fUpwindQuad_r[2] = ser_2x_p2_surfx1_eval_quad_node_2_l(fr); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_2x_p2_upwind_quad_to_modal(fUpwindQuad_l, fUpwind_l); 
  ser_2x_p2_upwind_quad_to_modal(fUpwindQuad_r, fUpwind_r); 

  Ghat_l[0] = 0.7071067811865475*alphaSurf_l[1]*fUpwind_l[1]+0.7071067811865475*alphaSurf_l[0]*fUpwind_l[0]; 
  Ghat_l[1] = 0.6324555320336759*alphaSurf_l[1]*fUpwind_l[2]+0.7071067811865475*alphaSurf_l[0]*fUpwind_l[1]+0.7071067811865475*fUpwind_l[0]*alphaSurf_l[1]; 
  Ghat_l[2] = 0.7071067811865475*alphaSurf_l[0]*fUpwind_l[2]+0.6324555320336759*alphaSurf_l[1]*fUpwind_l[1]; 

  Ghat_r[0] = 0.7071067811865475*alphaSurf_r[1]*fUpwind_r[1]+0.7071067811865475*alphaSurf_r[0]*fUpwind_r[0]; 
  Ghat_r[1] = 0.6324555320336759*alphaSurf_r[1]*fUpwind_r[2]+0.7071067811865475*alphaSurf_r[0]*fUpwind_r[1]+0.7071067811865475*fUpwind_r[0]*alphaSurf_r[1]; 
  Ghat_r[2] = 0.7071067811865475*alphaSurf_r[0]*fUpwind_r[2]+0.6324555320336759*alphaSurf_r[1]*fUpwind_r[1]; 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dx1; 
  out[1] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dx1; 
  out[2] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dx1; 
  out[3] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dx1; 
  out[4] += (1.58113883008419*Ghat_l[0]-1.58113883008419*Ghat_r[0])*dx1; 
  out[5] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dx1; 
  out[6] += (1.58113883008419*Ghat_l[1]-1.58113883008419*Ghat_r[1])*dx1; 
  out[7] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dx1; 

} 
