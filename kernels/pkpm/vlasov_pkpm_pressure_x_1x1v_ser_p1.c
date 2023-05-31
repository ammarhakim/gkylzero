#include <gkyl_vlasov_pkpm_kernels.h> 
#include <gkyl_basis_hyb_1x1v_p1_surfx1_eval_quad.h> 
#include <gkyl_basis_hyb_1x1v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH void vlasov_pkpm_pressure_x_1x1v_ser_p1(const double *w, const double *dxv, 
     const double *bvarl, const double *bvarc, const double *bvarr, 
     const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:                 Cell-center coordinates.
  // dxv[NDIM]:               Cell spacing.
  // bvarl/bvarc/bvarr:       Input magnetic field unit vector in left/center/right cells.
  // fl/fc/fr:                Input Distribution function [F_0, T_perp G = T_perp (F_1 - F_0)] in left/center/right cells.
  // out:                     Incremented distribution function in center cell.
  const double dx1 = 2.0/dxv[0]; 
  const double dvpar = dxv[1], wvpar = w[1]; 
  const double volFact = dxv[1]/2.0; 
  const double *bl = &bvarl[0]; 
  const double *bc = &bvarc[0]; 
  const double *br = &bvarr[0]; 
  const double *F_0l = &fl[0]; 
  const double *F_0c = &fc[0]; 
  const double *F_0r = &fr[0]; 
  double *out_pressure = &out[0]; 
  double alpha_l[6] = {0.0}; 
  double alpha_c[6] = {0.0}; 
  double alpha_r[6] = {0.0}; 
  alpha_l[0] = 1.414213562373095*bl[0]*wvpar; 
  alpha_l[1] = 1.414213562373095*bl[1]*wvpar; 
  alpha_l[2] = 0.408248290463863*bl[0]*dvpar; 
  alpha_l[3] = 0.408248290463863*bl[1]*dvpar; 

  alpha_c[0] = 1.414213562373095*bc[0]*wvpar; 
  alpha_c[1] = 1.414213562373095*bc[1]*wvpar; 
  alpha_c[2] = 0.408248290463863*bc[0]*dvpar; 
  alpha_c[3] = 0.408248290463863*bc[1]*dvpar; 

  alpha_r[0] = 1.414213562373095*br[0]*wvpar; 
  alpha_r[1] = 1.414213562373095*br[1]*wvpar; 
  alpha_r[2] = 0.408248290463863*br[0]*dvpar; 
  alpha_r[3] = 0.408248290463863*br[1]*dvpar; 

  double alphaSurf_l[3] = {0.0}; 
  alphaSurf_l[0] = 0.408248290463863*alpha_l[1]-0.408248290463863*alpha_c[1]+0.3535533905932737*alpha_l[0]+0.3535533905932737*alpha_c[0]; 
  alphaSurf_l[1] = 0.408248290463863*alpha_l[3]-0.408248290463863*alpha_c[3]+0.3535533905932737*alpha_l[2]+0.3535533905932737*alpha_c[2]; 

  double alphaSurf_r[3] = {0.0}; 
  alphaSurf_r[0] = (-0.408248290463863*alpha_r[1])+0.408248290463863*alpha_c[1]+0.3535533905932737*alpha_r[0]+0.3535533905932737*alpha_c[0]; 
  alphaSurf_r[1] = (-0.408248290463863*alpha_r[3])+0.408248290463863*alpha_c[3]+0.3535533905932737*alpha_r[2]+0.3535533905932737*alpha_c[2]; 

  double F_0_UpwindQuad_l[3] = {0.0};
  double F_0_UpwindQuad_r[3] = {0.0};
  double F_0_Upwind_l[3] = {0.0};
  double F_0_Upwind_r[3] = {0.0};
  double Ghat_F_0_l[3] = {0.0}; 
  double Ghat_F_0_r[3] = {0.0}; 
  if (0.7071067811865475*alphaSurf_l[0]-0.9486832980505137*alphaSurf_l[1] > 0) { 
    F_0_UpwindQuad_l[0] = hyb_1x1v_p1_surfx1_eval_quad_node_0_r(F_0l); 
  } else { 
    F_0_UpwindQuad_l[0] = hyb_1x1v_p1_surfx1_eval_quad_node_0_l(F_0c); 
  } 
  if (0.7071067811865475*alphaSurf_r[0]-0.9486832980505137*alphaSurf_r[1] > 0) { 
    F_0_UpwindQuad_r[0] = hyb_1x1v_p1_surfx1_eval_quad_node_0_r(F_0c); 
  } else { 
    F_0_UpwindQuad_r[0] = hyb_1x1v_p1_surfx1_eval_quad_node_0_l(F_0r); 
  } 
  if (0.7071067811865475*alphaSurf_l[0] > 0) { 
    F_0_UpwindQuad_l[1] = hyb_1x1v_p1_surfx1_eval_quad_node_1_r(F_0l); 
  } else { 
    F_0_UpwindQuad_l[1] = hyb_1x1v_p1_surfx1_eval_quad_node_1_l(F_0c); 
  } 
  if (0.7071067811865475*alphaSurf_r[0] > 0) { 
    F_0_UpwindQuad_r[1] = hyb_1x1v_p1_surfx1_eval_quad_node_1_r(F_0c); 
  } else { 
    F_0_UpwindQuad_r[1] = hyb_1x1v_p1_surfx1_eval_quad_node_1_l(F_0r); 
  } 
  if (0.9486832980505137*alphaSurf_l[1]+0.7071067811865475*alphaSurf_l[0] > 0) { 
    F_0_UpwindQuad_l[2] = hyb_1x1v_p1_surfx1_eval_quad_node_2_r(F_0l); 
  } else { 
    F_0_UpwindQuad_l[2] = hyb_1x1v_p1_surfx1_eval_quad_node_2_l(F_0c); 
  } 
  if (0.9486832980505137*alphaSurf_r[1]+0.7071067811865475*alphaSurf_r[0] > 0) { 
    F_0_UpwindQuad_r[2] = hyb_1x1v_p1_surfx1_eval_quad_node_2_r(F_0c); 
  } else { 
    F_0_UpwindQuad_r[2] = hyb_1x1v_p1_surfx1_eval_quad_node_2_l(F_0r); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_1x1v_p1_xdir_upwind_quad_to_modal(F_0_UpwindQuad_l, F_0_Upwind_l); 
  hyb_1x1v_p1_xdir_upwind_quad_to_modal(F_0_UpwindQuad_r, F_0_Upwind_r); 

  Ghat_F_0_l[0] = 0.7071067811865475*F_0_Upwind_l[1]*alphaSurf_l[1]+0.7071067811865475*F_0_Upwind_l[0]*alphaSurf_l[0]; 
  Ghat_F_0_l[1] = 0.6324555320336759*alphaSurf_l[1]*F_0_Upwind_l[2]+0.7071067811865475*F_0_Upwind_l[0]*alphaSurf_l[1]+0.7071067811865475*alphaSurf_l[0]*F_0_Upwind_l[1]; 
  Ghat_F_0_l[2] = 0.7071067811865475*alphaSurf_l[0]*F_0_Upwind_l[2]+0.6324555320336759*F_0_Upwind_l[1]*alphaSurf_l[1]; 

  Ghat_F_0_r[0] = 0.7071067811865475*F_0_Upwind_r[1]*alphaSurf_r[1]+0.7071067811865475*F_0_Upwind_r[0]*alphaSurf_r[0]; 
  Ghat_F_0_r[1] = 0.6324555320336759*alphaSurf_r[1]*F_0_Upwind_r[2]+0.7071067811865475*F_0_Upwind_r[0]*alphaSurf_r[1]+0.7071067811865475*alphaSurf_r[0]*F_0_Upwind_r[1]; 
  Ghat_F_0_r[2] = 0.7071067811865475*alphaSurf_r[0]*F_0_Upwind_r[2]+0.6324555320336759*F_0_Upwind_r[1]*alphaSurf_r[1]; 

  out_pressure[0] += dx1*volFact*((Ghat_F_0_r[0]-1.0*Ghat_F_0_l[0])*wvpar+(0.2886751345948129*Ghat_F_0_r[1]-0.2886751345948129*Ghat_F_0_l[1])*dvpar); 
  out_pressure[1] += dx1*volFact*((1.732050807568877*(Ghat_F_0_r[0]+Ghat_F_0_l[0])-1.224744871391589*(F_0c[3]*alpha_c[3]+F_0c[2]*alpha_c[2]+F_0c[1]*alpha_c[1]+F_0c[0]*alpha_c[0]))*wvpar+((-0.3162277660168379*alpha_c[3]*F_0c[5])-0.3162277660168379*alpha_c[2]*F_0c[4]-0.3535533905932737*(F_0c[1]*alpha_c[3]+alpha_c[1]*F_0c[3]+F_0c[0]*alpha_c[2]+alpha_c[0]*F_0c[2])+0.5*(Ghat_F_0_r[1]+Ghat_F_0_l[1]))*dvpar); 

} 
