#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_hyb_3x1v_p1_surfx1_eval_quad.h> 
#include <gkyl_basis_hyb_3x1v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH void vlasov_pkpm_surfx_3x1v_ser_p1(const double *w, const double *dxv, 
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
  const double dvpar = dxv[3], wvpar = w[3]; 
  const double *uc = &u_i[0]; 
  const double *bc = &bvar[0]; 
  double alphaSurf_l[12] = {0.0}; 
  alphaSurf_l[0] = (-1.732050807568877*bc[1]*wvpar)+bc[0]*wvpar-1.732050807568877*uc[1]+uc[0]; 
  alphaSurf_l[1] = (-1.732050807568877*bc[4]*wvpar)+bc[2]*wvpar-1.732050807568877*uc[4]+uc[2]; 
  alphaSurf_l[2] = (-1.732050807568877*bc[5]*wvpar)+bc[3]*wvpar-1.732050807568877*uc[5]+uc[3]; 
  alphaSurf_l[3] = 0.2886751345948129*bc[0]*dvpar-0.5*bc[1]*dvpar; 
  alphaSurf_l[4] = (-1.732050807568877*bc[7]*wvpar)+bc[6]*wvpar-1.732050807568877*uc[7]+uc[6]; 
  alphaSurf_l[5] = 0.2886751345948129*bc[2]*dvpar-0.5*bc[4]*dvpar; 
  alphaSurf_l[6] = 0.2886751345948129*bc[3]*dvpar-0.5*bc[5]*dvpar; 
  alphaSurf_l[7] = 0.2886751345948129*bc[6]*dvpar-0.5*bc[7]*dvpar; 

  double alphaSurf_r[12] = {0.0}; 
  alphaSurf_r[0] = 1.732050807568877*bc[1]*wvpar+bc[0]*wvpar+1.732050807568877*uc[1]+uc[0]; 
  alphaSurf_r[1] = 1.732050807568877*bc[4]*wvpar+bc[2]*wvpar+1.732050807568877*uc[4]+uc[2]; 
  alphaSurf_r[2] = 1.732050807568877*bc[5]*wvpar+bc[3]*wvpar+1.732050807568877*uc[5]+uc[3]; 
  alphaSurf_r[3] = 0.5*bc[1]*dvpar+0.2886751345948129*bc[0]*dvpar; 
  alphaSurf_r[4] = 1.732050807568877*bc[7]*wvpar+bc[6]*wvpar+1.732050807568877*uc[7]+uc[6]; 
  alphaSurf_r[5] = 0.5*bc[4]*dvpar+0.2886751345948129*bc[2]*dvpar; 
  alphaSurf_r[6] = 0.5*bc[5]*dvpar+0.2886751345948129*bc[3]*dvpar; 
  alphaSurf_r[7] = 0.5*bc[7]*dvpar+0.2886751345948129*bc[6]*dvpar; 

  double fUpwindQuad_l[12] = {0.0};
  double fUpwindQuad_r[12] = {0.0};
  double fUpwind_l[12] = {0.0};
  double fUpwind_r[12] = {0.0};
  double Ghat_l[12] = {0.0}; 
  double Ghat_r[12] = {0.0}; 

  if ((-0.4743416490252568*alphaSurf_l[7])+0.4743416490252568*(alphaSurf_l[6]+alphaSurf_l[5])+0.3535533905932737*alphaSurf_l[4]-0.4743416490252568*alphaSurf_l[3]-0.3535533905932737*(alphaSurf_l[2]+alphaSurf_l[1])+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[0] = hyb_3x1v_p1_surfx1_eval_quad_node_0_r(fl); 
  } else { 
    fUpwindQuad_l[0] = hyb_3x1v_p1_surfx1_eval_quad_node_0_l(fc); 
  } 
  if ((-0.4743416490252568*alphaSurf_r[7])+0.4743416490252568*(alphaSurf_r[6]+alphaSurf_r[5])+0.3535533905932737*alphaSurf_r[4]-0.4743416490252568*alphaSurf_r[3]-0.3535533905932737*(alphaSurf_r[2]+alphaSurf_r[1])+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[0] = hyb_3x1v_p1_surfx1_eval_quad_node_0_r(fc); 
  } else { 
    fUpwindQuad_r[0] = hyb_3x1v_p1_surfx1_eval_quad_node_0_l(fr); 
  } 
  if (0.3535533905932737*alphaSurf_l[4]-0.3535533905932737*(alphaSurf_l[2]+alphaSurf_l[1])+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[1] = hyb_3x1v_p1_surfx1_eval_quad_node_1_r(fl); 
  } else { 
    fUpwindQuad_l[1] = hyb_3x1v_p1_surfx1_eval_quad_node_1_l(fc); 
  } 
  if (0.3535533905932737*alphaSurf_r[4]-0.3535533905932737*(alphaSurf_r[2]+alphaSurf_r[1])+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[1] = hyb_3x1v_p1_surfx1_eval_quad_node_1_r(fc); 
  } else { 
    fUpwindQuad_r[1] = hyb_3x1v_p1_surfx1_eval_quad_node_1_l(fr); 
  } 
  if (0.4743416490252568*alphaSurf_l[7]-0.4743416490252568*(alphaSurf_l[6]+alphaSurf_l[5])+0.3535533905932737*alphaSurf_l[4]+0.4743416490252568*alphaSurf_l[3]-0.3535533905932737*(alphaSurf_l[2]+alphaSurf_l[1])+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[2] = hyb_3x1v_p1_surfx1_eval_quad_node_2_r(fl); 
  } else { 
    fUpwindQuad_l[2] = hyb_3x1v_p1_surfx1_eval_quad_node_2_l(fc); 
  } 
  if (0.4743416490252568*alphaSurf_r[7]-0.4743416490252568*(alphaSurf_r[6]+alphaSurf_r[5])+0.3535533905932737*alphaSurf_r[4]+0.4743416490252568*alphaSurf_r[3]-0.3535533905932737*(alphaSurf_r[2]+alphaSurf_r[1])+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[2] = hyb_3x1v_p1_surfx1_eval_quad_node_2_r(fc); 
  } else { 
    fUpwindQuad_r[2] = hyb_3x1v_p1_surfx1_eval_quad_node_2_l(fr); 
  } 
  if (0.4743416490252568*alphaSurf_l[7]-0.4743416490252568*alphaSurf_l[6]+0.4743416490252568*alphaSurf_l[5]-0.3535533905932737*alphaSurf_l[4]-0.4743416490252568*alphaSurf_l[3]+0.3535533905932737*alphaSurf_l[2]-0.3535533905932737*alphaSurf_l[1]+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[3] = hyb_3x1v_p1_surfx1_eval_quad_node_3_r(fl); 
  } else { 
    fUpwindQuad_l[3] = hyb_3x1v_p1_surfx1_eval_quad_node_3_l(fc); 
  } 
  if (0.4743416490252568*alphaSurf_r[7]-0.4743416490252568*alphaSurf_r[6]+0.4743416490252568*alphaSurf_r[5]-0.3535533905932737*alphaSurf_r[4]-0.4743416490252568*alphaSurf_r[3]+0.3535533905932737*alphaSurf_r[2]-0.3535533905932737*alphaSurf_r[1]+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[3] = hyb_3x1v_p1_surfx1_eval_quad_node_3_r(fc); 
  } else { 
    fUpwindQuad_r[3] = hyb_3x1v_p1_surfx1_eval_quad_node_3_l(fr); 
  } 
  if ((-0.3535533905932737*alphaSurf_l[4])+0.3535533905932737*alphaSurf_l[2]-0.3535533905932737*alphaSurf_l[1]+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[4] = hyb_3x1v_p1_surfx1_eval_quad_node_4_r(fl); 
  } else { 
    fUpwindQuad_l[4] = hyb_3x1v_p1_surfx1_eval_quad_node_4_l(fc); 
  } 
  if ((-0.3535533905932737*alphaSurf_r[4])+0.3535533905932737*alphaSurf_r[2]-0.3535533905932737*alphaSurf_r[1]+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[4] = hyb_3x1v_p1_surfx1_eval_quad_node_4_r(fc); 
  } else { 
    fUpwindQuad_r[4] = hyb_3x1v_p1_surfx1_eval_quad_node_4_l(fr); 
  } 
  if ((-0.4743416490252568*alphaSurf_l[7])+0.4743416490252568*alphaSurf_l[6]-0.4743416490252568*alphaSurf_l[5]-0.3535533905932737*alphaSurf_l[4]+0.4743416490252568*alphaSurf_l[3]+0.3535533905932737*alphaSurf_l[2]-0.3535533905932737*alphaSurf_l[1]+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[5] = hyb_3x1v_p1_surfx1_eval_quad_node_5_r(fl); 
  } else { 
    fUpwindQuad_l[5] = hyb_3x1v_p1_surfx1_eval_quad_node_5_l(fc); 
  } 
  if ((-0.4743416490252568*alphaSurf_r[7])+0.4743416490252568*alphaSurf_r[6]-0.4743416490252568*alphaSurf_r[5]-0.3535533905932737*alphaSurf_r[4]+0.4743416490252568*alphaSurf_r[3]+0.3535533905932737*alphaSurf_r[2]-0.3535533905932737*alphaSurf_r[1]+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[5] = hyb_3x1v_p1_surfx1_eval_quad_node_5_r(fc); 
  } else { 
    fUpwindQuad_r[5] = hyb_3x1v_p1_surfx1_eval_quad_node_5_l(fr); 
  } 
  if (0.4743416490252568*(alphaSurf_l[7]+alphaSurf_l[6])-0.4743416490252568*alphaSurf_l[5]-0.3535533905932737*alphaSurf_l[4]-0.4743416490252568*alphaSurf_l[3]-0.3535533905932737*alphaSurf_l[2]+0.3535533905932737*(alphaSurf_l[1]+alphaSurf_l[0]) > 0) { 
    fUpwindQuad_l[6] = hyb_3x1v_p1_surfx1_eval_quad_node_6_r(fl); 
  } else { 
    fUpwindQuad_l[6] = hyb_3x1v_p1_surfx1_eval_quad_node_6_l(fc); 
  } 
  if (0.4743416490252568*(alphaSurf_r[7]+alphaSurf_r[6])-0.4743416490252568*alphaSurf_r[5]-0.3535533905932737*alphaSurf_r[4]-0.4743416490252568*alphaSurf_r[3]-0.3535533905932737*alphaSurf_r[2]+0.3535533905932737*(alphaSurf_r[1]+alphaSurf_r[0]) > 0) { 
    fUpwindQuad_r[6] = hyb_3x1v_p1_surfx1_eval_quad_node_6_r(fc); 
  } else { 
    fUpwindQuad_r[6] = hyb_3x1v_p1_surfx1_eval_quad_node_6_l(fr); 
  } 
  if (0.3535533905932737*(alphaSurf_l[1]+alphaSurf_l[0])-0.3535533905932737*(alphaSurf_l[4]+alphaSurf_l[2]) > 0) { 
    fUpwindQuad_l[7] = hyb_3x1v_p1_surfx1_eval_quad_node_7_r(fl); 
  } else { 
    fUpwindQuad_l[7] = hyb_3x1v_p1_surfx1_eval_quad_node_7_l(fc); 
  } 
  if (0.3535533905932737*(alphaSurf_r[1]+alphaSurf_r[0])-0.3535533905932737*(alphaSurf_r[4]+alphaSurf_r[2]) > 0) { 
    fUpwindQuad_r[7] = hyb_3x1v_p1_surfx1_eval_quad_node_7_r(fc); 
  } else { 
    fUpwindQuad_r[7] = hyb_3x1v_p1_surfx1_eval_quad_node_7_l(fr); 
  } 
  if ((-0.4743416490252568*(alphaSurf_l[7]+alphaSurf_l[6]))+0.4743416490252568*alphaSurf_l[5]-0.3535533905932737*alphaSurf_l[4]+0.4743416490252568*alphaSurf_l[3]-0.3535533905932737*alphaSurf_l[2]+0.3535533905932737*(alphaSurf_l[1]+alphaSurf_l[0]) > 0) { 
    fUpwindQuad_l[8] = hyb_3x1v_p1_surfx1_eval_quad_node_8_r(fl); 
  } else { 
    fUpwindQuad_l[8] = hyb_3x1v_p1_surfx1_eval_quad_node_8_l(fc); 
  } 
  if ((-0.4743416490252568*(alphaSurf_r[7]+alphaSurf_r[6]))+0.4743416490252568*alphaSurf_r[5]-0.3535533905932737*alphaSurf_r[4]+0.4743416490252568*alphaSurf_r[3]-0.3535533905932737*alphaSurf_r[2]+0.3535533905932737*(alphaSurf_r[1]+alphaSurf_r[0]) > 0) { 
    fUpwindQuad_r[8] = hyb_3x1v_p1_surfx1_eval_quad_node_8_r(fc); 
  } else { 
    fUpwindQuad_r[8] = hyb_3x1v_p1_surfx1_eval_quad_node_8_l(fr); 
  } 
  if ((-0.4743416490252568*(alphaSurf_l[7]+alphaSurf_l[6]+alphaSurf_l[5]))+0.3535533905932737*alphaSurf_l[4]-0.4743416490252568*alphaSurf_l[3]+0.3535533905932737*(alphaSurf_l[2]+alphaSurf_l[1]+alphaSurf_l[0]) > 0) { 
    fUpwindQuad_l[9] = hyb_3x1v_p1_surfx1_eval_quad_node_9_r(fl); 
  } else { 
    fUpwindQuad_l[9] = hyb_3x1v_p1_surfx1_eval_quad_node_9_l(fc); 
  } 
  if ((-0.4743416490252568*(alphaSurf_r[7]+alphaSurf_r[6]+alphaSurf_r[5]))+0.3535533905932737*alphaSurf_r[4]-0.4743416490252568*alphaSurf_r[3]+0.3535533905932737*(alphaSurf_r[2]+alphaSurf_r[1]+alphaSurf_r[0]) > 0) { 
    fUpwindQuad_r[9] = hyb_3x1v_p1_surfx1_eval_quad_node_9_r(fc); 
  } else { 
    fUpwindQuad_r[9] = hyb_3x1v_p1_surfx1_eval_quad_node_9_l(fr); 
  } 
  if (0.3535533905932737*(alphaSurf_l[4]+alphaSurf_l[2]+alphaSurf_l[1]+alphaSurf_l[0]) > 0) { 
    fUpwindQuad_l[10] = hyb_3x1v_p1_surfx1_eval_quad_node_10_r(fl); 
  } else { 
    fUpwindQuad_l[10] = hyb_3x1v_p1_surfx1_eval_quad_node_10_l(fc); 
  } 
  if (0.3535533905932737*(alphaSurf_r[4]+alphaSurf_r[2]+alphaSurf_r[1]+alphaSurf_r[0]) > 0) { 
    fUpwindQuad_r[10] = hyb_3x1v_p1_surfx1_eval_quad_node_10_r(fc); 
  } else { 
    fUpwindQuad_r[10] = hyb_3x1v_p1_surfx1_eval_quad_node_10_l(fr); 
  } 
  if (0.4743416490252568*(alphaSurf_l[7]+alphaSurf_l[6]+alphaSurf_l[5])+0.3535533905932737*alphaSurf_l[4]+0.4743416490252568*alphaSurf_l[3]+0.3535533905932737*(alphaSurf_l[2]+alphaSurf_l[1]+alphaSurf_l[0]) > 0) { 
    fUpwindQuad_l[11] = hyb_3x1v_p1_surfx1_eval_quad_node_11_r(fl); 
  } else { 
    fUpwindQuad_l[11] = hyb_3x1v_p1_surfx1_eval_quad_node_11_l(fc); 
  } 
  if (0.4743416490252568*(alphaSurf_r[7]+alphaSurf_r[6]+alphaSurf_r[5])+0.3535533905932737*alphaSurf_r[4]+0.4743416490252568*alphaSurf_r[3]+0.3535533905932737*(alphaSurf_r[2]+alphaSurf_r[1]+alphaSurf_r[0]) > 0) { 
    fUpwindQuad_r[11] = hyb_3x1v_p1_surfx1_eval_quad_node_11_r(fc); 
  } else { 
    fUpwindQuad_r[11] = hyb_3x1v_p1_surfx1_eval_quad_node_11_l(fr); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_3x1v_p1_xdir_upwind_quad_to_modal(fUpwindQuad_l, fUpwind_l); 
  hyb_3x1v_p1_xdir_upwind_quad_to_modal(fUpwindQuad_r, fUpwind_r); 

  Ghat_l[0] = 0.3535533905932737*alphaSurf_l[7]*fUpwind_l[7]+0.3535533905932737*alphaSurf_l[6]*fUpwind_l[6]+0.3535533905932737*alphaSurf_l[5]*fUpwind_l[5]+0.3535533905932737*alphaSurf_l[4]*fUpwind_l[4]+0.3535533905932737*alphaSurf_l[3]*fUpwind_l[3]+0.3535533905932737*alphaSurf_l[2]*fUpwind_l[2]+0.3535533905932737*alphaSurf_l[1]*fUpwind_l[1]+0.3535533905932737*alphaSurf_l[0]*fUpwind_l[0]; 
  Ghat_l[1] = 0.3535533905932737*alphaSurf_l[6]*fUpwind_l[7]+0.3535533905932737*fUpwind_l[6]*alphaSurf_l[7]+0.3535533905932737*alphaSurf_l[3]*fUpwind_l[5]+0.3535533905932737*fUpwind_l[3]*alphaSurf_l[5]+0.3535533905932737*alphaSurf_l[2]*fUpwind_l[4]+0.3535533905932737*fUpwind_l[2]*alphaSurf_l[4]+0.3535533905932737*alphaSurf_l[0]*fUpwind_l[1]+0.3535533905932737*fUpwind_l[0]*alphaSurf_l[1]; 
  Ghat_l[2] = 0.3535533905932737*alphaSurf_l[5]*fUpwind_l[7]+0.3535533905932737*fUpwind_l[5]*alphaSurf_l[7]+0.3535533905932737*alphaSurf_l[3]*fUpwind_l[6]+0.3535533905932737*fUpwind_l[3]*alphaSurf_l[6]+0.3535533905932737*alphaSurf_l[1]*fUpwind_l[4]+0.3535533905932737*fUpwind_l[1]*alphaSurf_l[4]+0.3535533905932737*alphaSurf_l[0]*fUpwind_l[2]+0.3535533905932737*fUpwind_l[0]*alphaSurf_l[2]; 
  Ghat_l[3] = 0.3162277660168379*alphaSurf_l[7]*fUpwind_l[11]+0.3162277660168379*alphaSurf_l[6]*fUpwind_l[10]+0.3162277660168379*alphaSurf_l[5]*fUpwind_l[9]+0.3162277660168379*alphaSurf_l[3]*fUpwind_l[8]+0.3535533905932737*alphaSurf_l[4]*fUpwind_l[7]+0.3535533905932737*fUpwind_l[4]*alphaSurf_l[7]+0.3535533905932737*alphaSurf_l[2]*fUpwind_l[6]+0.3535533905932737*fUpwind_l[2]*alphaSurf_l[6]+0.3535533905932737*alphaSurf_l[1]*fUpwind_l[5]+0.3535533905932737*fUpwind_l[1]*alphaSurf_l[5]+0.3535533905932737*alphaSurf_l[0]*fUpwind_l[3]+0.3535533905932737*fUpwind_l[0]*alphaSurf_l[3]; 
  Ghat_l[4] = 0.3535533905932737*alphaSurf_l[3]*fUpwind_l[7]+0.3535533905932737*fUpwind_l[3]*alphaSurf_l[7]+0.3535533905932737*alphaSurf_l[5]*fUpwind_l[6]+0.3535533905932737*fUpwind_l[5]*alphaSurf_l[6]+0.3535533905932737*alphaSurf_l[0]*fUpwind_l[4]+0.3535533905932737*fUpwind_l[0]*alphaSurf_l[4]+0.3535533905932737*alphaSurf_l[1]*fUpwind_l[2]+0.3535533905932737*fUpwind_l[1]*alphaSurf_l[2]; 
  Ghat_l[5] = 0.3162277660168379*alphaSurf_l[6]*fUpwind_l[11]+0.3162277660168379*alphaSurf_l[7]*fUpwind_l[10]+0.3162277660168379*alphaSurf_l[3]*fUpwind_l[9]+0.3162277660168379*alphaSurf_l[5]*fUpwind_l[8]+0.3535533905932737*alphaSurf_l[2]*fUpwind_l[7]+0.3535533905932737*fUpwind_l[2]*alphaSurf_l[7]+0.3535533905932737*alphaSurf_l[4]*fUpwind_l[6]+0.3535533905932737*fUpwind_l[4]*alphaSurf_l[6]+0.3535533905932737*alphaSurf_l[0]*fUpwind_l[5]+0.3535533905932737*fUpwind_l[0]*alphaSurf_l[5]+0.3535533905932737*alphaSurf_l[1]*fUpwind_l[3]+0.3535533905932737*fUpwind_l[1]*alphaSurf_l[3]; 
  Ghat_l[6] = 0.3162277660168379*alphaSurf_l[5]*fUpwind_l[11]+0.3162277660168379*alphaSurf_l[3]*fUpwind_l[10]+0.3162277660168379*alphaSurf_l[7]*fUpwind_l[9]+0.3162277660168379*alphaSurf_l[6]*fUpwind_l[8]+0.3535533905932737*alphaSurf_l[1]*fUpwind_l[7]+0.3535533905932737*fUpwind_l[1]*alphaSurf_l[7]+0.3535533905932737*alphaSurf_l[0]*fUpwind_l[6]+0.3535533905932737*fUpwind_l[0]*alphaSurf_l[6]+0.3535533905932737*alphaSurf_l[4]*fUpwind_l[5]+0.3535533905932737*fUpwind_l[4]*alphaSurf_l[5]+0.3535533905932737*alphaSurf_l[2]*fUpwind_l[3]+0.3535533905932737*fUpwind_l[2]*alphaSurf_l[3]; 
  Ghat_l[7] = 0.3162277660168379*alphaSurf_l[3]*fUpwind_l[11]+0.3162277660168379*alphaSurf_l[5]*fUpwind_l[10]+0.3162277660168379*alphaSurf_l[6]*fUpwind_l[9]+0.3162277660168379*alphaSurf_l[7]*fUpwind_l[8]+0.3535533905932737*alphaSurf_l[0]*fUpwind_l[7]+0.3535533905932737*fUpwind_l[0]*alphaSurf_l[7]+0.3535533905932737*alphaSurf_l[1]*fUpwind_l[6]+0.3535533905932737*fUpwind_l[1]*alphaSurf_l[6]+0.3535533905932737*alphaSurf_l[2]*fUpwind_l[5]+0.3535533905932737*fUpwind_l[2]*alphaSurf_l[5]+0.3535533905932737*alphaSurf_l[3]*fUpwind_l[4]+0.3535533905932737*fUpwind_l[3]*alphaSurf_l[4]; 
  Ghat_l[8] = 0.3535533905932737*alphaSurf_l[4]*fUpwind_l[11]+0.3535533905932737*alphaSurf_l[2]*fUpwind_l[10]+0.3535533905932737*alphaSurf_l[1]*fUpwind_l[9]+0.3535533905932737*alphaSurf_l[0]*fUpwind_l[8]+0.3162277660168379*alphaSurf_l[7]*fUpwind_l[7]+0.3162277660168379*alphaSurf_l[6]*fUpwind_l[6]+0.3162277660168379*alphaSurf_l[5]*fUpwind_l[5]+0.3162277660168379*alphaSurf_l[3]*fUpwind_l[3]; 
  Ghat_l[9] = 0.3535533905932737*alphaSurf_l[2]*fUpwind_l[11]+0.3535533905932737*alphaSurf_l[4]*fUpwind_l[10]+0.3535533905932737*alphaSurf_l[0]*fUpwind_l[9]+0.3535533905932737*alphaSurf_l[1]*fUpwind_l[8]+0.3162277660168379*alphaSurf_l[6]*fUpwind_l[7]+0.3162277660168379*fUpwind_l[6]*alphaSurf_l[7]+0.3162277660168379*alphaSurf_l[3]*fUpwind_l[5]+0.3162277660168379*fUpwind_l[3]*alphaSurf_l[5]; 
  Ghat_l[10] = 0.3535533905932737*alphaSurf_l[1]*fUpwind_l[11]+0.3535533905932737*alphaSurf_l[0]*fUpwind_l[10]+0.3535533905932737*alphaSurf_l[4]*fUpwind_l[9]+0.3535533905932737*alphaSurf_l[2]*fUpwind_l[8]+0.3162277660168379*alphaSurf_l[5]*fUpwind_l[7]+0.3162277660168379*fUpwind_l[5]*alphaSurf_l[7]+0.3162277660168379*alphaSurf_l[3]*fUpwind_l[6]+0.3162277660168379*fUpwind_l[3]*alphaSurf_l[6]; 
  Ghat_l[11] = 0.3535533905932737*alphaSurf_l[0]*fUpwind_l[11]+0.3535533905932737*alphaSurf_l[1]*fUpwind_l[10]+0.3535533905932737*alphaSurf_l[2]*fUpwind_l[9]+0.3535533905932737*alphaSurf_l[4]*fUpwind_l[8]+0.3162277660168379*alphaSurf_l[3]*fUpwind_l[7]+0.3162277660168379*fUpwind_l[3]*alphaSurf_l[7]+0.3162277660168379*alphaSurf_l[5]*fUpwind_l[6]+0.3162277660168379*fUpwind_l[5]*alphaSurf_l[6]; 

  Ghat_r[0] = 0.3535533905932737*alphaSurf_r[7]*fUpwind_r[7]+0.3535533905932737*alphaSurf_r[6]*fUpwind_r[6]+0.3535533905932737*alphaSurf_r[5]*fUpwind_r[5]+0.3535533905932737*alphaSurf_r[4]*fUpwind_r[4]+0.3535533905932737*alphaSurf_r[3]*fUpwind_r[3]+0.3535533905932737*alphaSurf_r[2]*fUpwind_r[2]+0.3535533905932737*alphaSurf_r[1]*fUpwind_r[1]+0.3535533905932737*alphaSurf_r[0]*fUpwind_r[0]; 
  Ghat_r[1] = 0.3535533905932737*alphaSurf_r[6]*fUpwind_r[7]+0.3535533905932737*fUpwind_r[6]*alphaSurf_r[7]+0.3535533905932737*alphaSurf_r[3]*fUpwind_r[5]+0.3535533905932737*fUpwind_r[3]*alphaSurf_r[5]+0.3535533905932737*alphaSurf_r[2]*fUpwind_r[4]+0.3535533905932737*fUpwind_r[2]*alphaSurf_r[4]+0.3535533905932737*alphaSurf_r[0]*fUpwind_r[1]+0.3535533905932737*fUpwind_r[0]*alphaSurf_r[1]; 
  Ghat_r[2] = 0.3535533905932737*alphaSurf_r[5]*fUpwind_r[7]+0.3535533905932737*fUpwind_r[5]*alphaSurf_r[7]+0.3535533905932737*alphaSurf_r[3]*fUpwind_r[6]+0.3535533905932737*fUpwind_r[3]*alphaSurf_r[6]+0.3535533905932737*alphaSurf_r[1]*fUpwind_r[4]+0.3535533905932737*fUpwind_r[1]*alphaSurf_r[4]+0.3535533905932737*alphaSurf_r[0]*fUpwind_r[2]+0.3535533905932737*fUpwind_r[0]*alphaSurf_r[2]; 
  Ghat_r[3] = 0.3162277660168379*alphaSurf_r[7]*fUpwind_r[11]+0.3162277660168379*alphaSurf_r[6]*fUpwind_r[10]+0.3162277660168379*alphaSurf_r[5]*fUpwind_r[9]+0.3162277660168379*alphaSurf_r[3]*fUpwind_r[8]+0.3535533905932737*alphaSurf_r[4]*fUpwind_r[7]+0.3535533905932737*fUpwind_r[4]*alphaSurf_r[7]+0.3535533905932737*alphaSurf_r[2]*fUpwind_r[6]+0.3535533905932737*fUpwind_r[2]*alphaSurf_r[6]+0.3535533905932737*alphaSurf_r[1]*fUpwind_r[5]+0.3535533905932737*fUpwind_r[1]*alphaSurf_r[5]+0.3535533905932737*alphaSurf_r[0]*fUpwind_r[3]+0.3535533905932737*fUpwind_r[0]*alphaSurf_r[3]; 
  Ghat_r[4] = 0.3535533905932737*alphaSurf_r[3]*fUpwind_r[7]+0.3535533905932737*fUpwind_r[3]*alphaSurf_r[7]+0.3535533905932737*alphaSurf_r[5]*fUpwind_r[6]+0.3535533905932737*fUpwind_r[5]*alphaSurf_r[6]+0.3535533905932737*alphaSurf_r[0]*fUpwind_r[4]+0.3535533905932737*fUpwind_r[0]*alphaSurf_r[4]+0.3535533905932737*alphaSurf_r[1]*fUpwind_r[2]+0.3535533905932737*fUpwind_r[1]*alphaSurf_r[2]; 
  Ghat_r[5] = 0.3162277660168379*alphaSurf_r[6]*fUpwind_r[11]+0.3162277660168379*alphaSurf_r[7]*fUpwind_r[10]+0.3162277660168379*alphaSurf_r[3]*fUpwind_r[9]+0.3162277660168379*alphaSurf_r[5]*fUpwind_r[8]+0.3535533905932737*alphaSurf_r[2]*fUpwind_r[7]+0.3535533905932737*fUpwind_r[2]*alphaSurf_r[7]+0.3535533905932737*alphaSurf_r[4]*fUpwind_r[6]+0.3535533905932737*fUpwind_r[4]*alphaSurf_r[6]+0.3535533905932737*alphaSurf_r[0]*fUpwind_r[5]+0.3535533905932737*fUpwind_r[0]*alphaSurf_r[5]+0.3535533905932737*alphaSurf_r[1]*fUpwind_r[3]+0.3535533905932737*fUpwind_r[1]*alphaSurf_r[3]; 
  Ghat_r[6] = 0.3162277660168379*alphaSurf_r[5]*fUpwind_r[11]+0.3162277660168379*alphaSurf_r[3]*fUpwind_r[10]+0.3162277660168379*alphaSurf_r[7]*fUpwind_r[9]+0.3162277660168379*alphaSurf_r[6]*fUpwind_r[8]+0.3535533905932737*alphaSurf_r[1]*fUpwind_r[7]+0.3535533905932737*fUpwind_r[1]*alphaSurf_r[7]+0.3535533905932737*alphaSurf_r[0]*fUpwind_r[6]+0.3535533905932737*fUpwind_r[0]*alphaSurf_r[6]+0.3535533905932737*alphaSurf_r[4]*fUpwind_r[5]+0.3535533905932737*fUpwind_r[4]*alphaSurf_r[5]+0.3535533905932737*alphaSurf_r[2]*fUpwind_r[3]+0.3535533905932737*fUpwind_r[2]*alphaSurf_r[3]; 
  Ghat_r[7] = 0.3162277660168379*alphaSurf_r[3]*fUpwind_r[11]+0.3162277660168379*alphaSurf_r[5]*fUpwind_r[10]+0.3162277660168379*alphaSurf_r[6]*fUpwind_r[9]+0.3162277660168379*alphaSurf_r[7]*fUpwind_r[8]+0.3535533905932737*alphaSurf_r[0]*fUpwind_r[7]+0.3535533905932737*fUpwind_r[0]*alphaSurf_r[7]+0.3535533905932737*alphaSurf_r[1]*fUpwind_r[6]+0.3535533905932737*fUpwind_r[1]*alphaSurf_r[6]+0.3535533905932737*alphaSurf_r[2]*fUpwind_r[5]+0.3535533905932737*fUpwind_r[2]*alphaSurf_r[5]+0.3535533905932737*alphaSurf_r[3]*fUpwind_r[4]+0.3535533905932737*fUpwind_r[3]*alphaSurf_r[4]; 
  Ghat_r[8] = 0.3535533905932737*alphaSurf_r[4]*fUpwind_r[11]+0.3535533905932737*alphaSurf_r[2]*fUpwind_r[10]+0.3535533905932737*alphaSurf_r[1]*fUpwind_r[9]+0.3535533905932737*alphaSurf_r[0]*fUpwind_r[8]+0.3162277660168379*alphaSurf_r[7]*fUpwind_r[7]+0.3162277660168379*alphaSurf_r[6]*fUpwind_r[6]+0.3162277660168379*alphaSurf_r[5]*fUpwind_r[5]+0.3162277660168379*alphaSurf_r[3]*fUpwind_r[3]; 
  Ghat_r[9] = 0.3535533905932737*alphaSurf_r[2]*fUpwind_r[11]+0.3535533905932737*alphaSurf_r[4]*fUpwind_r[10]+0.3535533905932737*alphaSurf_r[0]*fUpwind_r[9]+0.3535533905932737*alphaSurf_r[1]*fUpwind_r[8]+0.3162277660168379*alphaSurf_r[6]*fUpwind_r[7]+0.3162277660168379*fUpwind_r[6]*alphaSurf_r[7]+0.3162277660168379*alphaSurf_r[3]*fUpwind_r[5]+0.3162277660168379*fUpwind_r[3]*alphaSurf_r[5]; 
  Ghat_r[10] = 0.3535533905932737*alphaSurf_r[1]*fUpwind_r[11]+0.3535533905932737*alphaSurf_r[0]*fUpwind_r[10]+0.3535533905932737*alphaSurf_r[4]*fUpwind_r[9]+0.3535533905932737*alphaSurf_r[2]*fUpwind_r[8]+0.3162277660168379*alphaSurf_r[5]*fUpwind_r[7]+0.3162277660168379*fUpwind_r[5]*alphaSurf_r[7]+0.3162277660168379*alphaSurf_r[3]*fUpwind_r[6]+0.3162277660168379*fUpwind_r[3]*alphaSurf_r[6]; 
  Ghat_r[11] = 0.3535533905932737*alphaSurf_r[0]*fUpwind_r[11]+0.3535533905932737*alphaSurf_r[1]*fUpwind_r[10]+0.3535533905932737*alphaSurf_r[2]*fUpwind_r[9]+0.3535533905932737*alphaSurf_r[4]*fUpwind_r[8]+0.3162277660168379*alphaSurf_r[3]*fUpwind_r[7]+0.3162277660168379*fUpwind_r[3]*alphaSurf_r[7]+0.3162277660168379*alphaSurf_r[5]*fUpwind_r[6]+0.3162277660168379*fUpwind_r[5]*alphaSurf_r[6]; 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dx1; 
  out[1] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dx1; 
  out[2] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dx1; 
  out[3] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dx1; 
  out[4] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dx1; 
  out[5] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dx1; 
  out[6] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dx1; 
  out[7] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*dx1; 
  out[8] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dx1; 
  out[9] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*dx1; 
  out[10] += (0.7071067811865475*Ghat_l[6]-0.7071067811865475*Ghat_r[6])*dx1; 
  out[11] += -1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dx1; 
  out[12] += -1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dx1; 
  out[13] += -1.224744871391589*(Ghat_r[6]+Ghat_l[6])*dx1; 
  out[14] += (0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*dx1; 
  out[15] += -1.224744871391589*(Ghat_r[7]+Ghat_l[7])*dx1; 
  out[16] += (0.7071067811865475*Ghat_l[8]-0.7071067811865475*Ghat_r[8])*dx1; 
  out[17] += -1.224744871391589*(Ghat_r[8]+Ghat_l[8])*dx1; 
  out[18] += (0.7071067811865475*Ghat_l[9]-0.7071067811865475*Ghat_r[9])*dx1; 
  out[19] += (0.7071067811865475*Ghat_l[10]-0.7071067811865475*Ghat_r[10])*dx1; 
  out[20] += -1.224744871391589*(Ghat_r[9]+Ghat_l[9])*dx1; 
  out[21] += -1.224744871391589*(Ghat_r[10]+Ghat_l[10])*dx1; 
  out[22] += (0.7071067811865475*Ghat_l[11]-0.7071067811865475*Ghat_r[11])*dx1; 
  out[23] += -1.224744871391589*(Ghat_r[11]+Ghat_l[11])*dx1; 

} 
