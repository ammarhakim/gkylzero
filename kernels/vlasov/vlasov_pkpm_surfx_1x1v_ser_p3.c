#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_ser_2x_p3_surfx1_eval_quad.h> 
#include <gkyl_basis_ser_2x_p3_upwind_quad_to_modal.h> 
GKYL_CU_DH void vlasov_pkpm_surfx_1x1v_ser_p3(const double *w, const double *dxv, 
     const double *bvarl, const double *bvarc, const double *bvarr, 
     const double *u_il, const double *u_ic, const double *u_ir, 
     const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // bvarl/bvarc/bvarr:  Input magnetic field unit vector in left/center/right cells.
  // u_il/u_ic/u_ir:  Input bulk velocity (ux,uy,uz) in left/center/right cells.
  // fl/fc/fr:  Input Distribution function in left/center/right cells.
  // out:       Incremented distribution function in center cell.
  const double dx1 = 2.0/dxv[0]; 
  const double dvpar = dxv[1], wvpar = w[1]; 
  const double *ul = &u_il[0]; 
  const double *uc = &u_ic[0]; 
  const double *ur = &u_ir[0]; 
  const double *bl = &bvarl[0]; 
  const double *bc = &bvarc[0]; 
  const double *br = &bvarr[0]; 
  double alpha_l[12] = {0.0}; 
  double alpha_c[12] = {0.0}; 
  double alpha_r[12] = {0.0}; 
  alpha_l[0] = 1.414213562373095*bl[0]*wvpar+1.414213562373095*ul[0]; 
  alpha_l[1] = 1.414213562373095*bl[1]*wvpar+1.414213562373095*ul[1]; 
  alpha_l[2] = 0.408248290463863*bl[0]*dvpar; 
  alpha_l[3] = 0.408248290463863*bl[1]*dvpar; 
  alpha_l[4] = 1.414213562373095*bl[2]*wvpar+1.414213562373095*ul[2]; 
  alpha_l[6] = 0.408248290463863*bl[2]*dvpar; 
  alpha_l[8] = 1.414213562373095*bl[3]*wvpar+1.414213562373095*ul[3]; 
  alpha_l[10] = 0.4082482904638629*bl[3]*dvpar; 

  alpha_c[0] = 1.414213562373095*bc[0]*wvpar+1.414213562373095*uc[0]; 
  alpha_c[1] = 1.414213562373095*bc[1]*wvpar+1.414213562373095*uc[1]; 
  alpha_c[2] = 0.408248290463863*bc[0]*dvpar; 
  alpha_c[3] = 0.408248290463863*bc[1]*dvpar; 
  alpha_c[4] = 1.414213562373095*bc[2]*wvpar+1.414213562373095*uc[2]; 
  alpha_c[6] = 0.408248290463863*bc[2]*dvpar; 
  alpha_c[8] = 1.414213562373095*bc[3]*wvpar+1.414213562373095*uc[3]; 
  alpha_c[10] = 0.4082482904638629*bc[3]*dvpar; 

  alpha_r[0] = 1.414213562373095*br[0]*wvpar+1.414213562373095*ur[0]; 
  alpha_r[1] = 1.414213562373095*br[1]*wvpar+1.414213562373095*ur[1]; 
  alpha_r[2] = 0.408248290463863*br[0]*dvpar; 
  alpha_r[3] = 0.408248290463863*br[1]*dvpar; 
  alpha_r[4] = 1.414213562373095*br[2]*wvpar+1.414213562373095*ur[2]; 
  alpha_r[6] = 0.408248290463863*br[2]*dvpar; 
  alpha_r[8] = 1.414213562373095*br[3]*wvpar+1.414213562373095*ur[3]; 
  alpha_r[10] = 0.4082482904638629*br[3]*dvpar; 

  double alphaSurf_l[4] = {0.0}; 
  alphaSurf_l[0] = 0.2672612419124243*alpha_l[8]-0.2672612419124243*alpha_c[8]+0.4941058844013091*alpha_l[4]+0.4941058844013091*alpha_c[4]+0.5358258812338199*alpha_l[1]-0.5358258812338199*alpha_c[1]+0.3535533905932737*alpha_l[0]+0.3535533905932737*alpha_c[0]; 
  alphaSurf_l[1] = 0.2672612419124243*alpha_l[10]-0.2672612419124243*alpha_c[10]+0.4941058844013091*alpha_l[6]+0.4941058844013091*alpha_c[6]+0.5358258812338199*alpha_l[3]-0.5358258812338199*alpha_c[3]+0.3535533905932737*alpha_l[2]+0.3535533905932737*alpha_c[2]; 

  double alphaSurf_r[4] = {0.0}; 
  alphaSurf_r[0] = (-0.2672612419124243*alpha_r[8])+0.2672612419124243*alpha_c[8]+0.4941058844013091*alpha_r[4]+0.4941058844013091*alpha_c[4]-0.5358258812338199*alpha_r[1]+0.5358258812338199*alpha_c[1]+0.3535533905932737*alpha_r[0]+0.3535533905932737*alpha_c[0]; 
  alphaSurf_r[1] = (-0.2672612419124243*alpha_r[10])+0.2672612419124243*alpha_c[10]+0.4941058844013091*alpha_r[6]+0.4941058844013091*alpha_c[6]-0.5358258812338199*alpha_r[3]+0.5358258812338199*alpha_c[3]+0.3535533905932737*alpha_r[2]+0.3535533905932737*alpha_c[2]; 

  double fUpwindQuad_l[4] = {0.0};
  double fUpwindQuad_r[4] = {0.0};
  double fUpwind_l[4] = {0.0};
  double fUpwind_r[4] = {0.0};
  double Ghat_l[4] = {0.0}; 
  double Ghat_r[4] = {0.0}; 

  if (0.7071067811865475*alphaSurf_l[0]-1.054672281193885*alphaSurf_l[1] > 0) { 
    fUpwindQuad_l[0] = ser_2x_p3_surfx1_eval_quad_node_0_r(fl); 
  } else { 
    fUpwindQuad_l[0] = ser_2x_p3_surfx1_eval_quad_node_0_l(fc); 
  } 
  if (0.7071067811865475*alphaSurf_r[0]-1.054672281193885*alphaSurf_r[1] > 0) { 
    fUpwindQuad_r[0] = ser_2x_p3_surfx1_eval_quad_node_0_r(fc); 
  } else { 
    fUpwindQuad_r[0] = ser_2x_p3_surfx1_eval_quad_node_0_l(fr); 
  } 
  if (0.7071067811865475*alphaSurf_l[0]-0.4163900395009129*alphaSurf_l[1] > 0) { 
    fUpwindQuad_l[1] = ser_2x_p3_surfx1_eval_quad_node_1_r(fl); 
  } else { 
    fUpwindQuad_l[1] = ser_2x_p3_surfx1_eval_quad_node_1_l(fc); 
  } 
  if (0.7071067811865475*alphaSurf_r[0]-0.4163900395009129*alphaSurf_r[1] > 0) { 
    fUpwindQuad_r[1] = ser_2x_p3_surfx1_eval_quad_node_1_r(fc); 
  } else { 
    fUpwindQuad_r[1] = ser_2x_p3_surfx1_eval_quad_node_1_l(fr); 
  } 
  if (0.4163900395009129*alphaSurf_l[1]+0.7071067811865475*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[2] = ser_2x_p3_surfx1_eval_quad_node_2_r(fl); 
  } else { 
    fUpwindQuad_l[2] = ser_2x_p3_surfx1_eval_quad_node_2_l(fc); 
  } 
  if (0.4163900395009129*alphaSurf_r[1]+0.7071067811865475*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[2] = ser_2x_p3_surfx1_eval_quad_node_2_r(fc); 
  } else { 
    fUpwindQuad_r[2] = ser_2x_p3_surfx1_eval_quad_node_2_l(fr); 
  } 
  if (1.054672281193885*alphaSurf_l[1]+0.7071067811865475*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[3] = ser_2x_p3_surfx1_eval_quad_node_3_r(fl); 
  } else { 
    fUpwindQuad_l[3] = ser_2x_p3_surfx1_eval_quad_node_3_l(fc); 
  } 
  if (1.054672281193885*alphaSurf_r[1]+0.7071067811865475*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[3] = ser_2x_p3_surfx1_eval_quad_node_3_r(fc); 
  } else { 
    fUpwindQuad_r[3] = ser_2x_p3_surfx1_eval_quad_node_3_l(fr); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_2x_p3_upwind_quad_to_modal(fUpwindQuad_l, fUpwind_l); 
  ser_2x_p3_upwind_quad_to_modal(fUpwindQuad_r, fUpwind_r); 

  Ghat_l[0] = 0.7071067811865475*alphaSurf_l[1]*fUpwind_l[1]+0.7071067811865475*alphaSurf_l[0]*fUpwind_l[0]; 
  Ghat_l[1] = 0.6324555320336759*alphaSurf_l[1]*fUpwind_l[2]+0.7071067811865475*alphaSurf_l[0]*fUpwind_l[1]+0.7071067811865475*fUpwind_l[0]*alphaSurf_l[1]; 
  Ghat_l[2] = 0.6210590034081186*alphaSurf_l[1]*fUpwind_l[3]+0.7071067811865475*alphaSurf_l[0]*fUpwind_l[2]+0.6324555320336759*alphaSurf_l[1]*fUpwind_l[1]; 
  Ghat_l[3] = 0.7071067811865475*alphaSurf_l[0]*fUpwind_l[3]+0.6210590034081186*alphaSurf_l[1]*fUpwind_l[2]; 

  Ghat_r[0] = 0.7071067811865475*alphaSurf_r[1]*fUpwind_r[1]+0.7071067811865475*alphaSurf_r[0]*fUpwind_r[0]; 
  Ghat_r[1] = 0.6324555320336759*alphaSurf_r[1]*fUpwind_r[2]+0.7071067811865475*alphaSurf_r[0]*fUpwind_r[1]+0.7071067811865475*fUpwind_r[0]*alphaSurf_r[1]; 
  Ghat_r[2] = 0.6210590034081186*alphaSurf_r[1]*fUpwind_r[3]+0.7071067811865475*alphaSurf_r[0]*fUpwind_r[2]+0.6324555320336759*alphaSurf_r[1]*fUpwind_r[1]; 
  Ghat_r[3] = 0.7071067811865475*alphaSurf_r[0]*fUpwind_r[3]+0.6210590034081186*alphaSurf_r[1]*fUpwind_r[2]; 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dx1; 
  out[1] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dx1; 
  out[2] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dx1; 
  out[3] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dx1; 
  out[4] += (1.58113883008419*Ghat_l[0]-1.58113883008419*Ghat_r[0])*dx1; 
  out[5] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dx1; 
  out[6] += (1.58113883008419*Ghat_l[1]-1.58113883008419*Ghat_r[1])*dx1; 
  out[7] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dx1; 
  out[8] += -1.870828693386971*(Ghat_r[0]+Ghat_l[0])*dx1; 
  out[9] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dx1; 
  out[10] += -1.870828693386971*(Ghat_r[1]+Ghat_l[1])*dx1; 
  out[11] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dx1; 

} 
