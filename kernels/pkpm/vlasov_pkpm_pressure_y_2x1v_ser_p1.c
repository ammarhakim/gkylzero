#include <gkyl_vlasov_pkpm_kernels.h> 
#include <gkyl_basis_hyb_2x1v_p1_surfx2_eval_quad.h> 
#include <gkyl_basis_hyb_2x1v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH void vlasov_pkpm_pressure_y_2x1v_ser_p1(const double *w, const double *dxv, 
     const double *bvarl, const double *bvarc, const double *bvarr, 
     const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:                 Cell-center coordinates.
  // dxv[NDIM]:               Cell spacing.
  // bvarl/bvarc/bvarr:       Input magnetic field unit vector in left/center/right cells.
  // fl/fc/fr:                Input Distribution function [F_0, T_perp G = T_perp (F_1 - F_0)] in left/center/right cells.
  // out:                     Incremented distribution function in center cell.
  const double dx1 = 2.0/dxv[1]; 
  const double dvpar = dxv[2], wvpar = w[2]; 
  const double volFact = dxv[2]/2.0; 
  const double *bl = &bvarl[4]; 
  const double *bc = &bvarc[4]; 
  const double *br = &bvarr[4]; 
  const double *F_0l = &fl[0]; 
  const double *F_0c = &fc[0]; 
  const double *F_0r = &fr[0]; 
  double *out_pressure = &out[4]; 
  double alpha_l[12] = {0.0}; 
  double alpha_c[12] = {0.0}; 
  double alpha_r[12] = {0.0}; 
  alpha_l[0] = 1.414213562373095*bl[0]*wvpar; 
  alpha_l[1] = 1.414213562373095*bl[1]*wvpar; 
  alpha_l[2] = 1.414213562373095*bl[2]*wvpar; 
  alpha_l[3] = 0.408248290463863*bl[0]*dvpar; 
  alpha_l[4] = 1.414213562373095*bl[3]*wvpar; 
  alpha_l[5] = 0.408248290463863*bl[1]*dvpar; 
  alpha_l[6] = 0.408248290463863*bl[2]*dvpar; 
  alpha_l[7] = 0.408248290463863*bl[3]*dvpar; 

  alpha_c[0] = 1.414213562373095*bc[0]*wvpar; 
  alpha_c[1] = 1.414213562373095*bc[1]*wvpar; 
  alpha_c[2] = 1.414213562373095*bc[2]*wvpar; 
  alpha_c[3] = 0.408248290463863*bc[0]*dvpar; 
  alpha_c[4] = 1.414213562373095*bc[3]*wvpar; 
  alpha_c[5] = 0.408248290463863*bc[1]*dvpar; 
  alpha_c[6] = 0.408248290463863*bc[2]*dvpar; 
  alpha_c[7] = 0.408248290463863*bc[3]*dvpar; 

  alpha_r[0] = 1.414213562373095*br[0]*wvpar; 
  alpha_r[1] = 1.414213562373095*br[1]*wvpar; 
  alpha_r[2] = 1.414213562373095*br[2]*wvpar; 
  alpha_r[3] = 0.408248290463863*br[0]*dvpar; 
  alpha_r[4] = 1.414213562373095*br[3]*wvpar; 
  alpha_r[5] = 0.408248290463863*br[1]*dvpar; 
  alpha_r[6] = 0.408248290463863*br[2]*dvpar; 
  alpha_r[7] = 0.408248290463863*br[3]*dvpar; 

  double alphaSurf_l[6] = {0.0}; 
  alphaSurf_l[0] = 0.408248290463863*alpha_l[2]-0.408248290463863*alpha_c[2]+0.3535533905932737*alpha_l[0]+0.3535533905932737*alpha_c[0]; 
  alphaSurf_l[1] = 0.408248290463863*alpha_l[4]-0.408248290463863*alpha_c[4]+0.3535533905932737*alpha_l[1]+0.3535533905932737*alpha_c[1]; 
  alphaSurf_l[2] = 0.408248290463863*alpha_l[6]-0.408248290463863*alpha_c[6]+0.3535533905932737*alpha_l[3]+0.3535533905932737*alpha_c[3]; 
  alphaSurf_l[3] = 0.408248290463863*alpha_l[7]-0.408248290463863*alpha_c[7]+0.3535533905932737*alpha_l[5]+0.3535533905932737*alpha_c[5]; 

  double alphaSurf_r[6] = {0.0}; 
  alphaSurf_r[0] = (-0.408248290463863*alpha_r[2])+0.408248290463863*alpha_c[2]+0.3535533905932737*alpha_r[0]+0.3535533905932737*alpha_c[0]; 
  alphaSurf_r[1] = (-0.408248290463863*alpha_r[4])+0.408248290463863*alpha_c[4]+0.3535533905932737*alpha_r[1]+0.3535533905932737*alpha_c[1]; 
  alphaSurf_r[2] = (-0.408248290463863*alpha_r[6])+0.408248290463863*alpha_c[6]+0.3535533905932737*alpha_r[3]+0.3535533905932737*alpha_c[3]; 
  alphaSurf_r[3] = (-0.408248290463863*alpha_r[7])+0.408248290463863*alpha_c[7]+0.3535533905932737*alpha_r[5]+0.3535533905932737*alpha_c[5]; 

  double F_0_UpwindQuad_l[6] = {0.0};
  double F_0_UpwindQuad_r[6] = {0.0};
  double F_0_Upwind_l[6] = {0.0};
  double F_0_Upwind_r[6] = {0.0};
  double Ghat_F_0_l[6] = {0.0}; 
  double Ghat_F_0_r[6] = {0.0}; 
  if (0.6708203932499369*alphaSurf_l[3]-0.6708203932499369*alphaSurf_l[2]-0.5*alphaSurf_l[1]+0.5*alphaSurf_l[0] > 0) { 
    F_0_UpwindQuad_l[0] = hyb_2x1v_p1_surfx2_eval_quad_node_0_r(F_0l); 
  } else { 
    F_0_UpwindQuad_l[0] = hyb_2x1v_p1_surfx2_eval_quad_node_0_l(F_0c); 
  } 
  if (0.6708203932499369*alphaSurf_r[3]-0.6708203932499369*alphaSurf_r[2]-0.5*alphaSurf_r[1]+0.5*alphaSurf_r[0] > 0) { 
    F_0_UpwindQuad_r[0] = hyb_2x1v_p1_surfx2_eval_quad_node_0_r(F_0c); 
  } else { 
    F_0_UpwindQuad_r[0] = hyb_2x1v_p1_surfx2_eval_quad_node_0_l(F_0r); 
  } 
  if (0.5*alphaSurf_l[0]-0.5*alphaSurf_l[1] > 0) { 
    F_0_UpwindQuad_l[1] = hyb_2x1v_p1_surfx2_eval_quad_node_1_r(F_0l); 
  } else { 
    F_0_UpwindQuad_l[1] = hyb_2x1v_p1_surfx2_eval_quad_node_1_l(F_0c); 
  } 
  if (0.5*alphaSurf_r[0]-0.5*alphaSurf_r[1] > 0) { 
    F_0_UpwindQuad_r[1] = hyb_2x1v_p1_surfx2_eval_quad_node_1_r(F_0c); 
  } else { 
    F_0_UpwindQuad_r[1] = hyb_2x1v_p1_surfx2_eval_quad_node_1_l(F_0r); 
  } 
  if ((-0.6708203932499369*alphaSurf_l[3])+0.6708203932499369*alphaSurf_l[2]-0.5*alphaSurf_l[1]+0.5*alphaSurf_l[0] > 0) { 
    F_0_UpwindQuad_l[2] = hyb_2x1v_p1_surfx2_eval_quad_node_2_r(F_0l); 
  } else { 
    F_0_UpwindQuad_l[2] = hyb_2x1v_p1_surfx2_eval_quad_node_2_l(F_0c); 
  } 
  if ((-0.6708203932499369*alphaSurf_r[3])+0.6708203932499369*alphaSurf_r[2]-0.5*alphaSurf_r[1]+0.5*alphaSurf_r[0] > 0) { 
    F_0_UpwindQuad_r[2] = hyb_2x1v_p1_surfx2_eval_quad_node_2_r(F_0c); 
  } else { 
    F_0_UpwindQuad_r[2] = hyb_2x1v_p1_surfx2_eval_quad_node_2_l(F_0r); 
  } 
  if (0.5*(alphaSurf_l[1]+alphaSurf_l[0])-0.6708203932499369*(alphaSurf_l[3]+alphaSurf_l[2]) > 0) { 
    F_0_UpwindQuad_l[3] = hyb_2x1v_p1_surfx2_eval_quad_node_3_r(F_0l); 
  } else { 
    F_0_UpwindQuad_l[3] = hyb_2x1v_p1_surfx2_eval_quad_node_3_l(F_0c); 
  } 
  if (0.5*(alphaSurf_r[1]+alphaSurf_r[0])-0.6708203932499369*(alphaSurf_r[3]+alphaSurf_r[2]) > 0) { 
    F_0_UpwindQuad_r[3] = hyb_2x1v_p1_surfx2_eval_quad_node_3_r(F_0c); 
  } else { 
    F_0_UpwindQuad_r[3] = hyb_2x1v_p1_surfx2_eval_quad_node_3_l(F_0r); 
  } 
  if (0.5*(alphaSurf_l[1]+alphaSurf_l[0]) > 0) { 
    F_0_UpwindQuad_l[4] = hyb_2x1v_p1_surfx2_eval_quad_node_4_r(F_0l); 
  } else { 
    F_0_UpwindQuad_l[4] = hyb_2x1v_p1_surfx2_eval_quad_node_4_l(F_0c); 
  } 
  if (0.5*(alphaSurf_r[1]+alphaSurf_r[0]) > 0) { 
    F_0_UpwindQuad_r[4] = hyb_2x1v_p1_surfx2_eval_quad_node_4_r(F_0c); 
  } else { 
    F_0_UpwindQuad_r[4] = hyb_2x1v_p1_surfx2_eval_quad_node_4_l(F_0r); 
  } 
  if (0.6708203932499369*(alphaSurf_l[3]+alphaSurf_l[2])+0.5*(alphaSurf_l[1]+alphaSurf_l[0]) > 0) { 
    F_0_UpwindQuad_l[5] = hyb_2x1v_p1_surfx2_eval_quad_node_5_r(F_0l); 
  } else { 
    F_0_UpwindQuad_l[5] = hyb_2x1v_p1_surfx2_eval_quad_node_5_l(F_0c); 
  } 
  if (0.6708203932499369*(alphaSurf_r[3]+alphaSurf_r[2])+0.5*(alphaSurf_r[1]+alphaSurf_r[0]) > 0) { 
    F_0_UpwindQuad_r[5] = hyb_2x1v_p1_surfx2_eval_quad_node_5_r(F_0c); 
  } else { 
    F_0_UpwindQuad_r[5] = hyb_2x1v_p1_surfx2_eval_quad_node_5_l(F_0r); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_2x1v_p1_xdir_upwind_quad_to_modal(F_0_UpwindQuad_l, F_0_Upwind_l); 
  hyb_2x1v_p1_xdir_upwind_quad_to_modal(F_0_UpwindQuad_r, F_0_Upwind_r); 

  Ghat_F_0_l[0] = 0.5*F_0_Upwind_l[3]*alphaSurf_l[3]+0.5*F_0_Upwind_l[2]*alphaSurf_l[2]+0.5*F_0_Upwind_l[1]*alphaSurf_l[1]+0.5*F_0_Upwind_l[0]*alphaSurf_l[0]; 
  Ghat_F_0_l[1] = 0.5*F_0_Upwind_l[2]*alphaSurf_l[3]+0.5*alphaSurf_l[2]*F_0_Upwind_l[3]+0.5*F_0_Upwind_l[0]*alphaSurf_l[1]+0.5*alphaSurf_l[0]*F_0_Upwind_l[1]; 
  Ghat_F_0_l[2] = 0.447213595499958*alphaSurf_l[3]*F_0_Upwind_l[5]+0.4472135954999579*alphaSurf_l[2]*F_0_Upwind_l[4]+0.5*F_0_Upwind_l[1]*alphaSurf_l[3]+0.5*alphaSurf_l[1]*F_0_Upwind_l[3]+0.5*F_0_Upwind_l[0]*alphaSurf_l[2]+0.5*alphaSurf_l[0]*F_0_Upwind_l[2]; 
  Ghat_F_0_l[3] = 0.447213595499958*alphaSurf_l[2]*F_0_Upwind_l[5]+0.4472135954999579*alphaSurf_l[3]*F_0_Upwind_l[4]+0.5*F_0_Upwind_l[0]*alphaSurf_l[3]+0.5*alphaSurf_l[0]*F_0_Upwind_l[3]+0.5*F_0_Upwind_l[1]*alphaSurf_l[2]+0.5*alphaSurf_l[1]*F_0_Upwind_l[2]; 
  Ghat_F_0_l[4] = 0.5000000000000001*alphaSurf_l[1]*F_0_Upwind_l[5]+0.5*alphaSurf_l[0]*F_0_Upwind_l[4]+0.4472135954999579*F_0_Upwind_l[3]*alphaSurf_l[3]+0.4472135954999579*F_0_Upwind_l[2]*alphaSurf_l[2]; 
  Ghat_F_0_l[5] = 0.5*alphaSurf_l[0]*F_0_Upwind_l[5]+0.5000000000000001*alphaSurf_l[1]*F_0_Upwind_l[4]+0.447213595499958*F_0_Upwind_l[2]*alphaSurf_l[3]+0.447213595499958*alphaSurf_l[2]*F_0_Upwind_l[3]; 

  Ghat_F_0_r[0] = 0.5*F_0_Upwind_r[3]*alphaSurf_r[3]+0.5*F_0_Upwind_r[2]*alphaSurf_r[2]+0.5*F_0_Upwind_r[1]*alphaSurf_r[1]+0.5*F_0_Upwind_r[0]*alphaSurf_r[0]; 
  Ghat_F_0_r[1] = 0.5*F_0_Upwind_r[2]*alphaSurf_r[3]+0.5*alphaSurf_r[2]*F_0_Upwind_r[3]+0.5*F_0_Upwind_r[0]*alphaSurf_r[1]+0.5*alphaSurf_r[0]*F_0_Upwind_r[1]; 
  Ghat_F_0_r[2] = 0.447213595499958*alphaSurf_r[3]*F_0_Upwind_r[5]+0.4472135954999579*alphaSurf_r[2]*F_0_Upwind_r[4]+0.5*F_0_Upwind_r[1]*alphaSurf_r[3]+0.5*alphaSurf_r[1]*F_0_Upwind_r[3]+0.5*F_0_Upwind_r[0]*alphaSurf_r[2]+0.5*alphaSurf_r[0]*F_0_Upwind_r[2]; 
  Ghat_F_0_r[3] = 0.447213595499958*alphaSurf_r[2]*F_0_Upwind_r[5]+0.4472135954999579*alphaSurf_r[3]*F_0_Upwind_r[4]+0.5*F_0_Upwind_r[0]*alphaSurf_r[3]+0.5*alphaSurf_r[0]*F_0_Upwind_r[3]+0.5*F_0_Upwind_r[1]*alphaSurf_r[2]+0.5*alphaSurf_r[1]*F_0_Upwind_r[2]; 
  Ghat_F_0_r[4] = 0.5000000000000001*alphaSurf_r[1]*F_0_Upwind_r[5]+0.5*alphaSurf_r[0]*F_0_Upwind_r[4]+0.4472135954999579*F_0_Upwind_r[3]*alphaSurf_r[3]+0.4472135954999579*F_0_Upwind_r[2]*alphaSurf_r[2]; 
  Ghat_F_0_r[5] = 0.5*alphaSurf_r[0]*F_0_Upwind_r[5]+0.5000000000000001*alphaSurf_r[1]*F_0_Upwind_r[4]+0.447213595499958*F_0_Upwind_r[2]*alphaSurf_r[3]+0.447213595499958*alphaSurf_r[2]*F_0_Upwind_r[3]; 

  out_pressure[0] += dx1*volFact*((Ghat_F_0_r[0]-1.0*Ghat_F_0_l[0])*wvpar+(0.2886751345948129*Ghat_F_0_r[2]-0.2886751345948129*Ghat_F_0_l[2])*dvpar); 
  out_pressure[1] += dx1*volFact*((Ghat_F_0_r[1]-1.0*Ghat_F_0_l[1])*wvpar+(0.2886751345948129*Ghat_F_0_r[3]-0.2886751345948129*Ghat_F_0_l[3])*dvpar); 
  out_pressure[2] += dx1*volFact*((1.732050807568877*(Ghat_F_0_r[0]+Ghat_F_0_l[0])-0.8660254037844386*(F_0c[7]*alpha_c[7]+F_0c[6]*alpha_c[6]+F_0c[5]*alpha_c[5]+F_0c[4]*alpha_c[4]+F_0c[3]*alpha_c[3]+F_0c[2]*alpha_c[2]+F_0c[1]*alpha_c[1]+F_0c[0]*alpha_c[0]))*wvpar+((-0.223606797749979*alpha_c[7]*F_0c[11])-0.223606797749979*(alpha_c[6]*F_0c[10]+alpha_c[5]*F_0c[9])-0.223606797749979*alpha_c[3]*F_0c[8]-0.25*(F_0c[4]*alpha_c[7]+alpha_c[4]*F_0c[7]+F_0c[2]*alpha_c[6]+alpha_c[2]*F_0c[6]+F_0c[1]*alpha_c[5]+alpha_c[1]*F_0c[5]+F_0c[0]*alpha_c[3]+alpha_c[0]*F_0c[3])+0.5*(Ghat_F_0_r[2]+Ghat_F_0_l[2]))*dvpar); 
  out_pressure[3] += dx1*volFact*(((-0.8660254037844386*(F_0c[6]*alpha_c[7]+alpha_c[6]*F_0c[7]+F_0c[3]*alpha_c[5]+alpha_c[3]*F_0c[5]+F_0c[2]*alpha_c[4]+alpha_c[2]*F_0c[4]+F_0c[0]*alpha_c[1]))+1.732050807568877*(Ghat_F_0_r[1]+Ghat_F_0_l[1])-0.8660254037844386*alpha_c[0]*F_0c[1])*wvpar+((-0.223606797749979*alpha_c[6]*F_0c[11])-0.223606797749979*(alpha_c[7]*F_0c[10]+alpha_c[3]*F_0c[9])-0.223606797749979*alpha_c[5]*F_0c[8]-0.25*(F_0c[2]*alpha_c[7]+alpha_c[2]*F_0c[7]+F_0c[4]*alpha_c[6]+alpha_c[4]*F_0c[6]+F_0c[0]*alpha_c[5]+alpha_c[0]*F_0c[5]+F_0c[1]*alpha_c[3])+0.5*(Ghat_F_0_r[3]+Ghat_F_0_l[3])-0.25*alpha_c[1]*F_0c[3])*dvpar); 

} 
