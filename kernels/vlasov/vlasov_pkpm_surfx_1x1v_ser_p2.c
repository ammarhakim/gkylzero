#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_ser_2x_p2_surfx1_eval_quad.h> 
#include <gkyl_basis_ser_2x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH void vlasov_pkpm_surfx_1x1v_ser_p2(const double *w, const double *dxv, 
     const double *u_il, const double *u_ic, const double *u_ir, const double *bvarl, const double *bvarc, const double *bvarr, 
     const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // u_il/u_ic/u_ir:  Input bulk velocity (ux,uy,uz) in left/center/right cells.
  // bvarl/bvarc/bvarr:  Input magnetic field unit vector in left/center/right cells.
  // fl/fc/fr:  Input Distribution function in left/center/right cells.
  // out:       Incremented distribution function in center cell.
  const double dx1 = 2.0/dxv[0]; 
  const double dvpar = dxv[1], wvpar = w[1]; 
  const double *ul_0 = &u_il[0]; 
  const double *uc_0 = &u_ic[0]; 
  const double *ur_0 = &u_ir[0]; 
  const double *bl_0 = &bvarl[0]; 
  const double *bc_0 = &bvarc[0]; 
  const double *br_0 = &bvarr[0]; 
  double alpha_lr[3] = {0.0}; 
  double alpha_cl[3] = {0.0}; 
  double alpha_cr[3] = {0.0}; 
  double alpha_rl[3] = {0.0}; 

  alpha_lr[0] = 2.23606797749979*bl_0[2]*wvpar+1.732050807568877*bl_0[1]*wvpar+bl_0[0]*wvpar+2.23606797749979*ul_0[2]+1.732050807568877*ul_0[1]+ul_0[0]; 
  alpha_lr[1] = 0.6454972243679029*bl_0[2]*dvpar+0.5*bl_0[1]*dvpar+0.2886751345948129*bl_0[0]*dvpar; 

  alpha_cl[0] = 2.23606797749979*bc_0[2]*wvpar-1.732050807568877*bc_0[1]*wvpar+bc_0[0]*wvpar+2.23606797749979*uc_0[2]-1.732050807568877*uc_0[1]+uc_0[0]; 
  alpha_cl[1] = 0.6454972243679029*bc_0[2]*dvpar-0.5*bc_0[1]*dvpar+0.2886751345948129*bc_0[0]*dvpar; 

  alpha_cr[0] = 2.23606797749979*bc_0[2]*wvpar+1.732050807568877*bc_0[1]*wvpar+bc_0[0]*wvpar+2.23606797749979*uc_0[2]+1.732050807568877*uc_0[1]+uc_0[0]; 
  alpha_cr[1] = 0.6454972243679029*bc_0[2]*dvpar+0.5*bc_0[1]*dvpar+0.2886751345948129*bc_0[0]*dvpar; 

  alpha_rl[0] = 2.23606797749979*br_0[2]*wvpar-1.732050807568877*br_0[1]*wvpar+br_0[0]*wvpar+2.23606797749979*ur_0[2]-1.732050807568877*ur_0[1]+ur_0[0]; 
  alpha_rl[1] = 0.6454972243679029*br_0[2]*dvpar-0.5*br_0[1]*dvpar+0.2886751345948129*br_0[0]*dvpar; 

  double alphaQuad_l[3] = {0.0};
  double alphaQuad_r[3] = {0.0};
  double alphaMax_l[3] = {0.0};;
  double alphaMax_r[3] = {0.0};
  double Ghat_l[3] = {0.0}; 
  double Ghat_r[3] = {0.0}; 
  double alpha_l_r = 0.0; 
  double alpha_c_l = 0.0; 
  double alpha_c_r = 0.0; 
  double alpha_r_l = 0.0; 

  alpha_l_r = ser_2x_p2_surfx1_eval_quad_node_0_r(alpha_lr); 
  alpha_c_l = ser_2x_p2_surfx1_eval_quad_node_0_l(alpha_cl); 
  alpha_c_r = ser_2x_p2_surfx1_eval_quad_node_0_r(alpha_cr); 
  alpha_r_l = ser_2x_p2_surfx1_eval_quad_node_0_l(alpha_rl); 
  alphaQuad_l[0] = fmax(fabs(alpha_l_r), fabs(alpha_c_l)); 
  alphaQuad_r[0] = fmax(fabs(alpha_c_r), fabs(alpha_r_l)); 
  alpha_l_r = ser_2x_p2_surfx1_eval_quad_node_1_r(alpha_lr); 
  alpha_c_l = ser_2x_p2_surfx1_eval_quad_node_1_l(alpha_cl); 
  alpha_c_r = ser_2x_p2_surfx1_eval_quad_node_1_r(alpha_cr); 
  alpha_r_l = ser_2x_p2_surfx1_eval_quad_node_1_l(alpha_rl); 
  alphaQuad_l[1] = fmax(fabs(alpha_l_r), fabs(alpha_c_l)); 
  alphaQuad_r[1] = fmax(fabs(alpha_c_r), fabs(alpha_r_l)); 
  alpha_l_r = ser_2x_p2_surfx1_eval_quad_node_2_r(alpha_lr); 
  alpha_c_l = ser_2x_p2_surfx1_eval_quad_node_2_l(alpha_cl); 
  alpha_c_r = ser_2x_p2_surfx1_eval_quad_node_2_r(alpha_cr); 
  alpha_r_l = ser_2x_p2_surfx1_eval_quad_node_2_l(alpha_rl); 
  alphaQuad_l[2] = fmax(fabs(alpha_l_r), fabs(alpha_c_l)); 
  alphaQuad_r[2] = fmax(fabs(alpha_c_r), fabs(alpha_r_l)); 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_2x_p2_upwind_quad_to_modal(alphaQuad_l, alphaMax_l); 
  ser_2x_p2_upwind_quad_to_modal(alphaQuad_r, alphaMax_r); 
  Ghat_l[0] = 0.4330127018922194*alphaMax_l[2]*fl[7]+0.4330127018922194*alphaMax_l[2]*fc[7]+0.5590169943749476*alpha_lr[1]*fl[6]+0.5590169943749476*alphaMax_l[1]*fl[6]+0.5590169943749476*alpha_cl[1]*fc[6]-0.5590169943749476*alphaMax_l[1]*fc[6]+0.25*alphaMax_l[2]*fl[5]-0.25*alphaMax_l[2]*fc[5]+0.5590169943749475*alpha_lr[0]*fl[4]+0.5590169943749475*alphaMax_l[0]*fl[4]+0.5590169943749475*alpha_cl[0]*fc[4]-0.5590169943749475*alphaMax_l[0]*fc[4]+0.4330127018922193*alpha_lr[1]*fl[3]+0.4330127018922193*alphaMax_l[1]*fl[3]-0.4330127018922193*alpha_cl[1]*fc[3]+0.4330127018922193*alphaMax_l[1]*fc[3]+0.25*alpha_lr[1]*fl[2]+0.25*alphaMax_l[1]*fl[2]+0.25*alpha_cl[1]*fc[2]-0.25*alphaMax_l[1]*fc[2]+0.4330127018922193*alpha_lr[0]*fl[1]+0.4330127018922193*alphaMax_l[0]*fl[1]-0.4330127018922193*alpha_cl[0]*fc[1]+0.4330127018922193*alphaMax_l[0]*fc[1]+0.25*alpha_lr[0]*fl[0]+0.25*alphaMax_l[0]*fl[0]+0.25*alpha_cl[0]*fc[0]-0.25*alphaMax_l[0]*fc[0]; 
  Ghat_l[1] = 0.3872983346207417*alpha_lr[1]*fl[7]+0.3872983346207417*alphaMax_l[1]*fl[7]-0.3872983346207417*alpha_cl[1]*fc[7]+0.3872983346207417*alphaMax_l[1]*fc[7]+0.5000000000000001*alphaMax_l[2]*fl[6]+0.5590169943749476*alpha_lr[0]*fl[6]+0.5590169943749476*alphaMax_l[0]*fl[6]-0.5000000000000001*alphaMax_l[2]*fc[6]+0.5590169943749476*alpha_cl[0]*fc[6]-0.5590169943749476*alphaMax_l[0]*fc[6]+0.223606797749979*alpha_lr[1]*fl[5]+0.223606797749979*alphaMax_l[1]*fl[5]+0.223606797749979*alpha_cl[1]*fc[5]-0.223606797749979*alphaMax_l[1]*fc[5]+0.5590169943749475*alpha_lr[1]*fl[4]+0.5590169943749475*alphaMax_l[1]*fl[4]+0.5590169943749475*alpha_cl[1]*fc[4]-0.5590169943749475*alphaMax_l[1]*fc[4]+0.3872983346207416*alphaMax_l[2]*fl[3]+0.4330127018922193*alpha_lr[0]*fl[3]+0.4330127018922193*alphaMax_l[0]*fl[3]+0.3872983346207416*alphaMax_l[2]*fc[3]-0.4330127018922193*alpha_cl[0]*fc[3]+0.4330127018922193*alphaMax_l[0]*fc[3]+0.223606797749979*alphaMax_l[2]*fl[2]+0.25*alpha_lr[0]*fl[2]+0.25*alphaMax_l[0]*fl[2]-0.223606797749979*alphaMax_l[2]*fc[2]+0.25*alpha_cl[0]*fc[2]-0.25*alphaMax_l[0]*fc[2]+0.4330127018922193*alpha_lr[1]*fl[1]+0.4330127018922193*alphaMax_l[1]*fl[1]-0.4330127018922193*alpha_cl[1]*fc[1]+0.4330127018922193*alphaMax_l[1]*fc[1]+0.25*fl[0]*alpha_lr[1]+0.25*fc[0]*alpha_cl[1]+0.25*fl[0]*alphaMax_l[1]-0.25*fc[0]*alphaMax_l[1]; 
  Ghat_l[2] = 0.276641667586244*alphaMax_l[2]*fl[7]+0.4330127018922194*alpha_lr[0]*fl[7]+0.4330127018922194*alphaMax_l[0]*fl[7]+0.276641667586244*alphaMax_l[2]*fc[7]-0.4330127018922194*alpha_cl[0]*fc[7]+0.4330127018922194*alphaMax_l[0]*fc[7]+0.5000000000000001*alpha_lr[1]*fl[6]+0.5000000000000001*alphaMax_l[1]*fl[6]+0.5000000000000001*alpha_cl[1]*fc[6]-0.5000000000000001*alphaMax_l[1]*fc[6]+0.159719141249985*alphaMax_l[2]*fl[5]+0.25*alpha_lr[0]*fl[5]+0.25*alphaMax_l[0]*fl[5]-0.159719141249985*alphaMax_l[2]*fc[5]+0.25*alpha_cl[0]*fc[5]-0.25*alphaMax_l[0]*fc[5]+0.5590169943749475*alphaMax_l[2]*fl[4]-0.5590169943749475*alphaMax_l[2]*fc[4]+0.3872983346207416*alpha_lr[1]*fl[3]+0.3872983346207416*alphaMax_l[1]*fl[3]-0.3872983346207416*alpha_cl[1]*fc[3]+0.3872983346207416*alphaMax_l[1]*fc[3]+0.223606797749979*alpha_lr[1]*fl[2]+0.223606797749979*alphaMax_l[1]*fl[2]+0.223606797749979*alpha_cl[1]*fc[2]-0.223606797749979*alphaMax_l[1]*fc[2]+0.4330127018922193*fl[1]*alphaMax_l[2]+0.4330127018922193*fc[1]*alphaMax_l[2]+0.25*fl[0]*alphaMax_l[2]-0.25*fc[0]*alphaMax_l[2]; 

  Ghat_r[0] = 0.4330127018922194*alphaMax_r[2]*fr[7]+0.4330127018922194*alphaMax_r[2]*fc[7]+0.5590169943749476*alpha_rl[1]*fr[6]-0.5590169943749476*alphaMax_r[1]*fr[6]+0.5590169943749476*alpha_cr[1]*fc[6]+0.5590169943749476*alphaMax_r[1]*fc[6]-0.25*alphaMax_r[2]*fr[5]+0.25*alphaMax_r[2]*fc[5]+0.5590169943749475*alpha_rl[0]*fr[4]-0.5590169943749475*alphaMax_r[0]*fr[4]+0.5590169943749475*alpha_cr[0]*fc[4]+0.5590169943749475*alphaMax_r[0]*fc[4]-0.4330127018922193*alpha_rl[1]*fr[3]+0.4330127018922193*alphaMax_r[1]*fr[3]+0.4330127018922193*alpha_cr[1]*fc[3]+0.4330127018922193*alphaMax_r[1]*fc[3]+0.25*alpha_rl[1]*fr[2]-0.25*alphaMax_r[1]*fr[2]+0.25*alpha_cr[1]*fc[2]+0.25*alphaMax_r[1]*fc[2]-0.4330127018922193*alpha_rl[0]*fr[1]+0.4330127018922193*alphaMax_r[0]*fr[1]+0.4330127018922193*alpha_cr[0]*fc[1]+0.4330127018922193*alphaMax_r[0]*fc[1]+0.25*alpha_rl[0]*fr[0]-0.25*alphaMax_r[0]*fr[0]+0.25*alpha_cr[0]*fc[0]+0.25*alphaMax_r[0]*fc[0]; 
  Ghat_r[1] = (-0.3872983346207417*alpha_rl[1]*fr[7])+0.3872983346207417*alphaMax_r[1]*fr[7]+0.3872983346207417*alpha_cr[1]*fc[7]+0.3872983346207417*alphaMax_r[1]*fc[7]-0.5000000000000001*alphaMax_r[2]*fr[6]+0.5590169943749476*alpha_rl[0]*fr[6]-0.5590169943749476*alphaMax_r[0]*fr[6]+0.5000000000000001*alphaMax_r[2]*fc[6]+0.5590169943749476*alpha_cr[0]*fc[6]+0.5590169943749476*alphaMax_r[0]*fc[6]+0.223606797749979*alpha_rl[1]*fr[5]-0.223606797749979*alphaMax_r[1]*fr[5]+0.223606797749979*alpha_cr[1]*fc[5]+0.223606797749979*alphaMax_r[1]*fc[5]+0.5590169943749475*alpha_rl[1]*fr[4]-0.5590169943749475*alphaMax_r[1]*fr[4]+0.5590169943749475*alpha_cr[1]*fc[4]+0.5590169943749475*alphaMax_r[1]*fc[4]+0.3872983346207416*alphaMax_r[2]*fr[3]-0.4330127018922193*alpha_rl[0]*fr[3]+0.4330127018922193*alphaMax_r[0]*fr[3]+0.3872983346207416*alphaMax_r[2]*fc[3]+0.4330127018922193*alpha_cr[0]*fc[3]+0.4330127018922193*alphaMax_r[0]*fc[3]-0.223606797749979*alphaMax_r[2]*fr[2]+0.25*alpha_rl[0]*fr[2]-0.25*alphaMax_r[0]*fr[2]+0.223606797749979*alphaMax_r[2]*fc[2]+0.25*alpha_cr[0]*fc[2]+0.25*alphaMax_r[0]*fc[2]-0.4330127018922193*alpha_rl[1]*fr[1]+0.4330127018922193*alphaMax_r[1]*fr[1]+0.4330127018922193*alpha_cr[1]*fc[1]+0.4330127018922193*alphaMax_r[1]*fc[1]+0.25*fr[0]*alpha_rl[1]+0.25*fc[0]*alpha_cr[1]-0.25*fr[0]*alphaMax_r[1]+0.25*fc[0]*alphaMax_r[1]; 
  Ghat_r[2] = 0.276641667586244*alphaMax_r[2]*fr[7]-0.4330127018922194*alpha_rl[0]*fr[7]+0.4330127018922194*alphaMax_r[0]*fr[7]+0.276641667586244*alphaMax_r[2]*fc[7]+0.4330127018922194*alpha_cr[0]*fc[7]+0.4330127018922194*alphaMax_r[0]*fc[7]+0.5000000000000001*alpha_rl[1]*fr[6]-0.5000000000000001*alphaMax_r[1]*fr[6]+0.5000000000000001*alpha_cr[1]*fc[6]+0.5000000000000001*alphaMax_r[1]*fc[6]-0.159719141249985*alphaMax_r[2]*fr[5]+0.25*alpha_rl[0]*fr[5]-0.25*alphaMax_r[0]*fr[5]+0.159719141249985*alphaMax_r[2]*fc[5]+0.25*alpha_cr[0]*fc[5]+0.25*alphaMax_r[0]*fc[5]-0.5590169943749475*alphaMax_r[2]*fr[4]+0.5590169943749475*alphaMax_r[2]*fc[4]-0.3872983346207416*alpha_rl[1]*fr[3]+0.3872983346207416*alphaMax_r[1]*fr[3]+0.3872983346207416*alpha_cr[1]*fc[3]+0.3872983346207416*alphaMax_r[1]*fc[3]+0.223606797749979*alpha_rl[1]*fr[2]-0.223606797749979*alphaMax_r[1]*fr[2]+0.223606797749979*alpha_cr[1]*fc[2]+0.223606797749979*alphaMax_r[1]*fc[2]+0.4330127018922193*fr[1]*alphaMax_r[2]+0.4330127018922193*fc[1]*alphaMax_r[2]-0.25*fr[0]*alphaMax_r[2]+0.25*fc[0]*alphaMax_r[2]; 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dx1; 
  out[1] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dx1; 
  out[2] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dx1; 
  out[3] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dx1; 
  out[4] += (1.58113883008419*Ghat_l[0]-1.58113883008419*Ghat_r[0])*dx1; 
  out[5] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dx1; 
  out[6] += (1.58113883008419*Ghat_l[1]-1.58113883008419*Ghat_r[1])*dx1; 
  out[7] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dx1; 

} 
