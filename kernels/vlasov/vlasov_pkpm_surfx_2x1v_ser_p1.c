#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_hyb_2x1v_p1_surfx1_eval_quad.h> 
#include <gkyl_basis_hyb_2x1v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH void vlasov_pkpm_surfx_2x1v_ser_p1(const double *w, const double *dxv, 
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
  const double dvpar = dxv[2], wvpar = w[2]; 
  const double *ul = &u_il[0]; 
  const double *uc = &u_ic[0]; 
  const double *ur = &u_ir[0]; 
  const double *bl = &bvarl[0]; 
  const double *bc = &bvarc[0]; 
  const double *br = &bvarr[0]; 
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
  alphaSurf_l[0] = 0.408248290463863*alpha_l[1]-0.408248290463863*alpha_c[1]+0.3535533905932737*alpha_l[0]+0.3535533905932737*alpha_c[0]; 
  alphaSurf_l[1] = 0.408248290463863*alpha_l[4]-0.408248290463863*alpha_c[4]+0.3535533905932737*alpha_l[2]+0.3535533905932737*alpha_c[2]; 
  alphaSurf_l[2] = 0.408248290463863*alpha_l[5]-0.408248290463863*alpha_c[5]+0.3535533905932737*alpha_l[3]+0.3535533905932737*alpha_c[3]; 
  alphaSurf_l[3] = 0.408248290463863*alpha_l[7]-0.408248290463863*alpha_c[7]+0.3535533905932737*alpha_l[6]+0.3535533905932737*alpha_c[6]; 

  double alphaSurf_r[6] = {0.0}; 
  alphaSurf_r[0] = (-0.408248290463863*alpha_r[1])+0.408248290463863*alpha_c[1]+0.3535533905932737*alpha_r[0]+0.3535533905932737*alpha_c[0]; 
  alphaSurf_r[1] = (-0.408248290463863*alpha_r[4])+0.408248290463863*alpha_c[4]+0.3535533905932737*alpha_r[2]+0.3535533905932737*alpha_c[2]; 
  alphaSurf_r[2] = (-0.408248290463863*alpha_r[5])+0.408248290463863*alpha_c[5]+0.3535533905932737*alpha_r[3]+0.3535533905932737*alpha_c[3]; 
  alphaSurf_r[3] = (-0.408248290463863*alpha_r[7])+0.408248290463863*alpha_c[7]+0.3535533905932737*alpha_r[6]+0.3535533905932737*alpha_c[6]; 

  double fUpwindQuad_l[6] = {0.0};
  double fUpwindQuad_r[6] = {0.0};
  double fUpwind_l[6] = {0.0};
  double fUpwind_r[6] = {0.0};
  double Ghat_l[6] = {0.0}; 
  double Ghat_r[6] = {0.0}; 

  if (0.6708203932499369*alphaSurf_l[3]-0.6708203932499369*alphaSurf_l[2]-0.5*alphaSurf_l[1]+0.5*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[0] = hyb_2x1v_p1_surfx1_eval_quad_node_0_r(fl); 
  } else { 
    fUpwindQuad_l[0] = hyb_2x1v_p1_surfx1_eval_quad_node_0_l(fc); 
  } 
  if (0.6708203932499369*alphaSurf_r[3]-0.6708203932499369*alphaSurf_r[2]-0.5*alphaSurf_r[1]+0.5*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[0] = hyb_2x1v_p1_surfx1_eval_quad_node_0_r(fc); 
  } else { 
    fUpwindQuad_r[0] = hyb_2x1v_p1_surfx1_eval_quad_node_0_l(fr); 
  } 
  if (0.5*alphaSurf_l[0]-0.5*alphaSurf_l[1] > 0) { 
    fUpwindQuad_l[1] = hyb_2x1v_p1_surfx1_eval_quad_node_1_r(fl); 
  } else { 
    fUpwindQuad_l[1] = hyb_2x1v_p1_surfx1_eval_quad_node_1_l(fc); 
  } 
  if (0.5*alphaSurf_r[0]-0.5*alphaSurf_r[1] > 0) { 
    fUpwindQuad_r[1] = hyb_2x1v_p1_surfx1_eval_quad_node_1_r(fc); 
  } else { 
    fUpwindQuad_r[1] = hyb_2x1v_p1_surfx1_eval_quad_node_1_l(fr); 
  } 
  if ((-0.6708203932499369*alphaSurf_l[3])+0.6708203932499369*alphaSurf_l[2]-0.5*alphaSurf_l[1]+0.5*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[2] = hyb_2x1v_p1_surfx1_eval_quad_node_2_r(fl); 
  } else { 
    fUpwindQuad_l[2] = hyb_2x1v_p1_surfx1_eval_quad_node_2_l(fc); 
  } 
  if ((-0.6708203932499369*alphaSurf_r[3])+0.6708203932499369*alphaSurf_r[2]-0.5*alphaSurf_r[1]+0.5*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[2] = hyb_2x1v_p1_surfx1_eval_quad_node_2_r(fc); 
  } else { 
    fUpwindQuad_r[2] = hyb_2x1v_p1_surfx1_eval_quad_node_2_l(fr); 
  } 
  if (0.5*(alphaSurf_l[1]+alphaSurf_l[0])-0.6708203932499369*(alphaSurf_l[3]+alphaSurf_l[2]) > 0) { 
    fUpwindQuad_l[3] = hyb_2x1v_p1_surfx1_eval_quad_node_3_r(fl); 
  } else { 
    fUpwindQuad_l[3] = hyb_2x1v_p1_surfx1_eval_quad_node_3_l(fc); 
  } 
  if (0.5*(alphaSurf_r[1]+alphaSurf_r[0])-0.6708203932499369*(alphaSurf_r[3]+alphaSurf_r[2]) > 0) { 
    fUpwindQuad_r[3] = hyb_2x1v_p1_surfx1_eval_quad_node_3_r(fc); 
  } else { 
    fUpwindQuad_r[3] = hyb_2x1v_p1_surfx1_eval_quad_node_3_l(fr); 
  } 
  if (0.5*(alphaSurf_l[1]+alphaSurf_l[0]) > 0) { 
    fUpwindQuad_l[4] = hyb_2x1v_p1_surfx1_eval_quad_node_4_r(fl); 
  } else { 
    fUpwindQuad_l[4] = hyb_2x1v_p1_surfx1_eval_quad_node_4_l(fc); 
  } 
  if (0.5*(alphaSurf_r[1]+alphaSurf_r[0]) > 0) { 
    fUpwindQuad_r[4] = hyb_2x1v_p1_surfx1_eval_quad_node_4_r(fc); 
  } else { 
    fUpwindQuad_r[4] = hyb_2x1v_p1_surfx1_eval_quad_node_4_l(fr); 
  } 
  if (0.6708203932499369*(alphaSurf_l[3]+alphaSurf_l[2])+0.5*(alphaSurf_l[1]+alphaSurf_l[0]) > 0) { 
    fUpwindQuad_l[5] = hyb_2x1v_p1_surfx1_eval_quad_node_5_r(fl); 
  } else { 
    fUpwindQuad_l[5] = hyb_2x1v_p1_surfx1_eval_quad_node_5_l(fc); 
  } 
  if (0.6708203932499369*(alphaSurf_r[3]+alphaSurf_r[2])+0.5*(alphaSurf_r[1]+alphaSurf_r[0]) > 0) { 
    fUpwindQuad_r[5] = hyb_2x1v_p1_surfx1_eval_quad_node_5_r(fc); 
  } else { 
    fUpwindQuad_r[5] = hyb_2x1v_p1_surfx1_eval_quad_node_5_l(fr); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_2x1v_p1_xdir_upwind_quad_to_modal(fUpwindQuad_l, fUpwind_l); 
  hyb_2x1v_p1_xdir_upwind_quad_to_modal(fUpwindQuad_r, fUpwind_r); 

  Ghat_l[0] = 0.5*alphaSurf_l[3]*fUpwind_l[3]+0.5*alphaSurf_l[2]*fUpwind_l[2]+0.5*alphaSurf_l[1]*fUpwind_l[1]+0.5*alphaSurf_l[0]*fUpwind_l[0]; 
  Ghat_l[1] = 0.5*alphaSurf_l[2]*fUpwind_l[3]+0.5*fUpwind_l[2]*alphaSurf_l[3]+0.5*alphaSurf_l[0]*fUpwind_l[1]+0.5*fUpwind_l[0]*alphaSurf_l[1]; 
  Ghat_l[2] = 0.447213595499958*alphaSurf_l[3]*fUpwind_l[5]+0.4472135954999579*alphaSurf_l[2]*fUpwind_l[4]+0.5*alphaSurf_l[1]*fUpwind_l[3]+0.5*fUpwind_l[1]*alphaSurf_l[3]+0.5*alphaSurf_l[0]*fUpwind_l[2]+0.5*fUpwind_l[0]*alphaSurf_l[2]; 
  Ghat_l[3] = 0.447213595499958*alphaSurf_l[2]*fUpwind_l[5]+0.4472135954999579*alphaSurf_l[3]*fUpwind_l[4]+0.5*alphaSurf_l[0]*fUpwind_l[3]+0.5*fUpwind_l[0]*alphaSurf_l[3]+0.5*alphaSurf_l[1]*fUpwind_l[2]+0.5*fUpwind_l[1]*alphaSurf_l[2]; 
  Ghat_l[4] = 0.5000000000000001*alphaSurf_l[1]*fUpwind_l[5]+0.5*alphaSurf_l[0]*fUpwind_l[4]+0.4472135954999579*alphaSurf_l[3]*fUpwind_l[3]+0.4472135954999579*alphaSurf_l[2]*fUpwind_l[2]; 
  Ghat_l[5] = 0.5*alphaSurf_l[0]*fUpwind_l[5]+0.5000000000000001*alphaSurf_l[1]*fUpwind_l[4]+0.447213595499958*alphaSurf_l[2]*fUpwind_l[3]+0.447213595499958*fUpwind_l[2]*alphaSurf_l[3]; 

  Ghat_r[0] = 0.5*alphaSurf_r[3]*fUpwind_r[3]+0.5*alphaSurf_r[2]*fUpwind_r[2]+0.5*alphaSurf_r[1]*fUpwind_r[1]+0.5*alphaSurf_r[0]*fUpwind_r[0]; 
  Ghat_r[1] = 0.5*alphaSurf_r[2]*fUpwind_r[3]+0.5*fUpwind_r[2]*alphaSurf_r[3]+0.5*alphaSurf_r[0]*fUpwind_r[1]+0.5*fUpwind_r[0]*alphaSurf_r[1]; 
  Ghat_r[2] = 0.447213595499958*alphaSurf_r[3]*fUpwind_r[5]+0.4472135954999579*alphaSurf_r[2]*fUpwind_r[4]+0.5*alphaSurf_r[1]*fUpwind_r[3]+0.5*fUpwind_r[1]*alphaSurf_r[3]+0.5*alphaSurf_r[0]*fUpwind_r[2]+0.5*fUpwind_r[0]*alphaSurf_r[2]; 
  Ghat_r[3] = 0.447213595499958*alphaSurf_r[2]*fUpwind_r[5]+0.4472135954999579*alphaSurf_r[3]*fUpwind_r[4]+0.5*alphaSurf_r[0]*fUpwind_r[3]+0.5*fUpwind_r[0]*alphaSurf_r[3]+0.5*alphaSurf_r[1]*fUpwind_r[2]+0.5*fUpwind_r[1]*alphaSurf_r[2]; 
  Ghat_r[4] = 0.5000000000000001*alphaSurf_r[1]*fUpwind_r[5]+0.5*alphaSurf_r[0]*fUpwind_r[4]+0.4472135954999579*alphaSurf_r[3]*fUpwind_r[3]+0.4472135954999579*alphaSurf_r[2]*fUpwind_r[2]; 
  Ghat_r[5] = 0.5*alphaSurf_r[0]*fUpwind_r[5]+0.5000000000000001*alphaSurf_r[1]*fUpwind_r[4]+0.447213595499958*alphaSurf_r[2]*fUpwind_r[3]+0.447213595499958*fUpwind_r[2]*alphaSurf_r[3]; 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dx1; 
  out[1] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dx1; 
  out[2] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dx1; 
  out[3] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dx1; 
  out[4] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dx1; 
  out[5] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dx1; 
  out[6] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dx1; 
  out[7] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dx1; 
  out[8] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*dx1; 
  out[9] += -1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dx1; 
  out[10] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*dx1; 
  out[11] += -1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dx1; 

  double alpha_u_l[12] = {0.0}; 
  double alpha_u_c[12] = {0.0}; 
  double alpha_u_r[12] = {0.0}; 
  alpha_u_l[0] = 1.414213562373095*ul[0]; 
  alpha_u_l[1] = 1.414213562373095*ul[1]; 
  alpha_u_l[2] = 1.414213562373095*ul[2]; 
  alpha_u_l[4] = 1.414213562373095*ul[3]; 

  alpha_u_c[0] = 1.414213562373095*uc[0]; 
  alpha_u_c[1] = 1.414213562373095*uc[1]; 
  alpha_u_c[2] = 1.414213562373095*uc[2]; 
  alpha_u_c[4] = 1.414213562373095*uc[3]; 

  alpha_u_r[0] = 1.414213562373095*ur[0]; 
  alpha_u_r[1] = 1.414213562373095*ur[1]; 
  alpha_u_r[2] = 1.414213562373095*ur[2]; 
  alpha_u_r[4] = 1.414213562373095*ur[3]; 

  double alphaSurf_lr[6] = {0.0}; 
  double alphaSurf_cl[6] = {0.0}; 
  double alphaSurf_cr[6] = {0.0}; 
  double alphaSurf_rl[6] = {0.0}; 
  alphaSurf_lr[0] = 1.732050807568877*ul[1]+ul[0]; 
  alphaSurf_lr[1] = 1.732050807568877*ul[3]+ul[2]; 

  alphaSurf_cl[0] = uc[0]-1.732050807568877*uc[1]; 
  alphaSurf_cl[1] = uc[2]-1.732050807568877*uc[3]; 

  alphaSurf_cr[0] = 1.732050807568877*uc[1]+uc[0]; 
  alphaSurf_cr[1] = 1.732050807568877*uc[3]+uc[2]; 

  alphaSurf_rl[0] = ur[0]-1.732050807568877*ur[1]; 
  alphaSurf_rl[1] = ur[2]-1.732050807568877*ur[3]; 

  double fSurf_lr[6] = {0.0}; 
  double fSurf_cl[6] = {0.0}; 
  double fSurf_cr[6] = {0.0}; 
  double fSurf_rl[6] = {0.0}; 
  fSurf_lr[0] = 1.224744871391589*fl[1]+0.7071067811865475*fl[0]; 
  fSurf_lr[1] = 1.224744871391589*fl[4]+0.7071067811865475*fl[2]; 
  fSurf_lr[2] = 1.224744871391589*fl[5]+0.7071067811865475*fl[3]; 
  fSurf_lr[3] = 1.224744871391589*fl[7]+0.7071067811865475*fl[6]; 
  fSurf_lr[4] = 1.224744871391589*fl[9]+0.7071067811865475*fl[8]; 
  fSurf_lr[5] = 1.224744871391589*fl[11]+0.7071067811865475*fl[10]; 

  fSurf_cl[0] = 0.7071067811865475*fc[0]-1.224744871391589*fc[1]; 
  fSurf_cl[1] = 0.7071067811865475*fc[2]-1.224744871391589*fc[4]; 
  fSurf_cl[2] = 0.7071067811865475*fc[3]-1.224744871391589*fc[5]; 
  fSurf_cl[3] = 0.7071067811865475*fc[6]-1.224744871391589*fc[7]; 
  fSurf_cl[4] = 0.7071067811865475*fc[8]-1.224744871391589*fc[9]; 
  fSurf_cl[5] = 0.7071067811865475*fc[10]-1.224744871391589*fc[11]; 

  fSurf_cr[0] = 1.224744871391589*fc[1]+0.7071067811865475*fc[0]; 
  fSurf_cr[1] = 1.224744871391589*fc[4]+0.7071067811865475*fc[2]; 
  fSurf_cr[2] = 1.224744871391589*fc[5]+0.7071067811865475*fc[3]; 
  fSurf_cr[3] = 1.224744871391589*fc[7]+0.7071067811865475*fc[6]; 
  fSurf_cr[4] = 1.224744871391589*fc[9]+0.7071067811865475*fc[8]; 
  fSurf_cr[5] = 1.224744871391589*fc[11]+0.7071067811865475*fc[10]; 

  fSurf_rl[0] = 0.7071067811865475*fr[0]-1.224744871391589*fr[1]; 
  fSurf_rl[1] = 0.7071067811865475*fr[2]-1.224744871391589*fr[4]; 
  fSurf_rl[2] = 0.7071067811865475*fr[3]-1.224744871391589*fr[5]; 
  fSurf_rl[3] = 0.7071067811865475*fr[6]-1.224744871391589*fr[7]; 
  fSurf_rl[4] = 0.7071067811865475*fr[8]-1.224744871391589*fr[9]; 
  fSurf_rl[5] = 0.7071067811865475*fr[10]-1.224744871391589*fr[11]; 

  double jump_f_quad_l[6] = {0.0};
  double jump_f_quad_r[6] = {0.0};
  double jump_f_modal_l[6] = {0.0};
  double jump_f_modal_r[6] = {0.0};

  double Ghat_u_l[6] = {0.0}; 
  double Ghat_u_r[6] = {0.0}; 

  double alpha_l_r = 0.0; 
  double alpha_c_l = 0.0; 
  double alpha_c_r = 0.0; 
  double alpha_r_l = 0.0; 
  double alphaQuad_l = 0.0; 
  double alphaQuad_r = 0.0; 
  double f_l_r = 0.0; 
  double f_c_l = 0.0; 
  double f_c_r = 0.0; 
  double f_r_l = 0.0; 
  alpha_l_r = hyb_2x1v_p1_surfx1_eval_quad_node_0_r(alpha_u_l); 
  alpha_c_l = hyb_2x1v_p1_surfx1_eval_quad_node_0_l(alpha_u_c); 
  alpha_c_r = hyb_2x1v_p1_surfx1_eval_quad_node_0_r(alpha_u_c); 
  alpha_r_l = hyb_2x1v_p1_surfx1_eval_quad_node_0_l(alpha_u_r); 
  alphaQuad_l = fmax(fabs(alpha_l_r), fabs(alpha_c_l)); 
  alphaQuad_r = fmax(fabs(alpha_c_r), fabs(alpha_r_l)); 
  f_l_r = hyb_2x1v_p1_surfx1_eval_quad_node_0_r(fl); 
  f_c_l = hyb_2x1v_p1_surfx1_eval_quad_node_0_l(fc); 
  f_c_r = hyb_2x1v_p1_surfx1_eval_quad_node_0_r(fc); 
  f_r_l = hyb_2x1v_p1_surfx1_eval_quad_node_0_l(fr); 
  jump_f_quad_l[0] = 0.5*alphaQuad_l*(f_c_l - f_l_r); 
  jump_f_quad_r[0] = 0.5*alphaQuad_r*(f_r_l - f_c_r); 

  alpha_l_r = hyb_2x1v_p1_surfx1_eval_quad_node_1_r(alpha_u_l); 
  alpha_c_l = hyb_2x1v_p1_surfx1_eval_quad_node_1_l(alpha_u_c); 
  alpha_c_r = hyb_2x1v_p1_surfx1_eval_quad_node_1_r(alpha_u_c); 
  alpha_r_l = hyb_2x1v_p1_surfx1_eval_quad_node_1_l(alpha_u_r); 
  alphaQuad_l = fmax(fabs(alpha_l_r), fabs(alpha_c_l)); 
  alphaQuad_r = fmax(fabs(alpha_c_r), fabs(alpha_r_l)); 
  f_l_r = hyb_2x1v_p1_surfx1_eval_quad_node_1_r(fl); 
  f_c_l = hyb_2x1v_p1_surfx1_eval_quad_node_1_l(fc); 
  f_c_r = hyb_2x1v_p1_surfx1_eval_quad_node_1_r(fc); 
  f_r_l = hyb_2x1v_p1_surfx1_eval_quad_node_1_l(fr); 
  jump_f_quad_l[1] = 0.5*alphaQuad_l*(f_c_l - f_l_r); 
  jump_f_quad_r[1] = 0.5*alphaQuad_r*(f_r_l - f_c_r); 

  alpha_l_r = hyb_2x1v_p1_surfx1_eval_quad_node_2_r(alpha_u_l); 
  alpha_c_l = hyb_2x1v_p1_surfx1_eval_quad_node_2_l(alpha_u_c); 
  alpha_c_r = hyb_2x1v_p1_surfx1_eval_quad_node_2_r(alpha_u_c); 
  alpha_r_l = hyb_2x1v_p1_surfx1_eval_quad_node_2_l(alpha_u_r); 
  alphaQuad_l = fmax(fabs(alpha_l_r), fabs(alpha_c_l)); 
  alphaQuad_r = fmax(fabs(alpha_c_r), fabs(alpha_r_l)); 
  f_l_r = hyb_2x1v_p1_surfx1_eval_quad_node_2_r(fl); 
  f_c_l = hyb_2x1v_p1_surfx1_eval_quad_node_2_l(fc); 
  f_c_r = hyb_2x1v_p1_surfx1_eval_quad_node_2_r(fc); 
  f_r_l = hyb_2x1v_p1_surfx1_eval_quad_node_2_l(fr); 
  jump_f_quad_l[2] = 0.5*alphaQuad_l*(f_c_l - f_l_r); 
  jump_f_quad_r[2] = 0.5*alphaQuad_r*(f_r_l - f_c_r); 

  alpha_l_r = hyb_2x1v_p1_surfx1_eval_quad_node_3_r(alpha_u_l); 
  alpha_c_l = hyb_2x1v_p1_surfx1_eval_quad_node_3_l(alpha_u_c); 
  alpha_c_r = hyb_2x1v_p1_surfx1_eval_quad_node_3_r(alpha_u_c); 
  alpha_r_l = hyb_2x1v_p1_surfx1_eval_quad_node_3_l(alpha_u_r); 
  alphaQuad_l = fmax(fabs(alpha_l_r), fabs(alpha_c_l)); 
  alphaQuad_r = fmax(fabs(alpha_c_r), fabs(alpha_r_l)); 
  f_l_r = hyb_2x1v_p1_surfx1_eval_quad_node_3_r(fl); 
  f_c_l = hyb_2x1v_p1_surfx1_eval_quad_node_3_l(fc); 
  f_c_r = hyb_2x1v_p1_surfx1_eval_quad_node_3_r(fc); 
  f_r_l = hyb_2x1v_p1_surfx1_eval_quad_node_3_l(fr); 
  jump_f_quad_l[3] = 0.5*alphaQuad_l*(f_c_l - f_l_r); 
  jump_f_quad_r[3] = 0.5*alphaQuad_r*(f_r_l - f_c_r); 

  alpha_l_r = hyb_2x1v_p1_surfx1_eval_quad_node_4_r(alpha_u_l); 
  alpha_c_l = hyb_2x1v_p1_surfx1_eval_quad_node_4_l(alpha_u_c); 
  alpha_c_r = hyb_2x1v_p1_surfx1_eval_quad_node_4_r(alpha_u_c); 
  alpha_r_l = hyb_2x1v_p1_surfx1_eval_quad_node_4_l(alpha_u_r); 
  alphaQuad_l = fmax(fabs(alpha_l_r), fabs(alpha_c_l)); 
  alphaQuad_r = fmax(fabs(alpha_c_r), fabs(alpha_r_l)); 
  f_l_r = hyb_2x1v_p1_surfx1_eval_quad_node_4_r(fl); 
  f_c_l = hyb_2x1v_p1_surfx1_eval_quad_node_4_l(fc); 
  f_c_r = hyb_2x1v_p1_surfx1_eval_quad_node_4_r(fc); 
  f_r_l = hyb_2x1v_p1_surfx1_eval_quad_node_4_l(fr); 
  jump_f_quad_l[4] = 0.5*alphaQuad_l*(f_c_l - f_l_r); 
  jump_f_quad_r[4] = 0.5*alphaQuad_r*(f_r_l - f_c_r); 

  alpha_l_r = hyb_2x1v_p1_surfx1_eval_quad_node_5_r(alpha_u_l); 
  alpha_c_l = hyb_2x1v_p1_surfx1_eval_quad_node_5_l(alpha_u_c); 
  alpha_c_r = hyb_2x1v_p1_surfx1_eval_quad_node_5_r(alpha_u_c); 
  alpha_r_l = hyb_2x1v_p1_surfx1_eval_quad_node_5_l(alpha_u_r); 
  alphaQuad_l = fmax(fabs(alpha_l_r), fabs(alpha_c_l)); 
  alphaQuad_r = fmax(fabs(alpha_c_r), fabs(alpha_r_l)); 
  f_l_r = hyb_2x1v_p1_surfx1_eval_quad_node_5_r(fl); 
  f_c_l = hyb_2x1v_p1_surfx1_eval_quad_node_5_l(fc); 
  f_c_r = hyb_2x1v_p1_surfx1_eval_quad_node_5_r(fc); 
  f_r_l = hyb_2x1v_p1_surfx1_eval_quad_node_5_l(fr); 
  jump_f_quad_l[5] = 0.5*alphaQuad_l*(f_c_l - f_l_r); 
  jump_f_quad_r[5] = 0.5*alphaQuad_r*(f_r_l - f_c_r); 


  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_2x1v_p1_xdir_upwind_quad_to_modal(jump_f_quad_l, jump_f_modal_l); 
  hyb_2x1v_p1_xdir_upwind_quad_to_modal(jump_f_quad_r, jump_f_modal_r); 
  Ghat_u_l[0] = 0.25*alphaSurf_lr[1]*fSurf_lr[1]+0.25*alphaSurf_cl[1]*fSurf_cl[1]-1.0*jump_f_modal_l[0]+0.25*alphaSurf_lr[0]*fSurf_lr[0]+0.25*alphaSurf_cl[0]*fSurf_cl[0]; 
  Ghat_u_l[1] = (-1.0*jump_f_modal_l[1])+0.25*alphaSurf_lr[0]*fSurf_lr[1]+0.25*alphaSurf_cl[0]*fSurf_cl[1]+0.25*fSurf_lr[0]*alphaSurf_lr[1]+0.25*fSurf_cl[0]*alphaSurf_cl[1]; 
  Ghat_u_l[2] = 0.25*alphaSurf_lr[1]*fSurf_lr[3]+0.25*alphaSurf_cl[1]*fSurf_cl[3]-1.0*jump_f_modal_l[2]+0.25*alphaSurf_lr[0]*fSurf_lr[2]+0.25*alphaSurf_cl[0]*fSurf_cl[2]; 
  Ghat_u_l[3] = (-1.0*jump_f_modal_l[3])+0.25*alphaSurf_lr[0]*fSurf_lr[3]+0.25*alphaSurf_cl[0]*fSurf_cl[3]+0.25*alphaSurf_lr[1]*fSurf_lr[2]+0.25*alphaSurf_cl[1]*fSurf_cl[2]; 
  Ghat_u_l[4] = 0.2500000000000001*alphaSurf_lr[1]*fSurf_lr[5]+0.2500000000000001*alphaSurf_cl[1]*fSurf_cl[5]-1.0*jump_f_modal_l[4]+0.25*alphaSurf_lr[0]*fSurf_lr[4]+0.25*alphaSurf_cl[0]*fSurf_cl[4]; 
  Ghat_u_l[5] = (-1.0*jump_f_modal_l[5])+0.25*alphaSurf_lr[0]*fSurf_lr[5]+0.25*alphaSurf_cl[0]*fSurf_cl[5]+0.2500000000000001*alphaSurf_lr[1]*fSurf_lr[4]+0.2500000000000001*alphaSurf_cl[1]*fSurf_cl[4]; 

  Ghat_u_r[0] = 0.25*alphaSurf_rl[1]*fSurf_rl[1]+0.25*alphaSurf_cr[1]*fSurf_cr[1]-1.0*jump_f_modal_r[0]+0.25*alphaSurf_rl[0]*fSurf_rl[0]+0.25*alphaSurf_cr[0]*fSurf_cr[0]; 
  Ghat_u_r[1] = (-1.0*jump_f_modal_r[1])+0.25*alphaSurf_rl[0]*fSurf_rl[1]+0.25*alphaSurf_cr[0]*fSurf_cr[1]+0.25*fSurf_rl[0]*alphaSurf_rl[1]+0.25*fSurf_cr[0]*alphaSurf_cr[1]; 
  Ghat_u_r[2] = 0.25*alphaSurf_rl[1]*fSurf_rl[3]+0.25*alphaSurf_cr[1]*fSurf_cr[3]-1.0*jump_f_modal_r[2]+0.25*alphaSurf_rl[0]*fSurf_rl[2]+0.25*alphaSurf_cr[0]*fSurf_cr[2]; 
  Ghat_u_r[3] = (-1.0*jump_f_modal_r[3])+0.25*alphaSurf_rl[0]*fSurf_rl[3]+0.25*alphaSurf_cr[0]*fSurf_cr[3]+0.25*alphaSurf_rl[1]*fSurf_rl[2]+0.25*alphaSurf_cr[1]*fSurf_cr[2]; 
  Ghat_u_r[4] = 0.2500000000000001*alphaSurf_rl[1]*fSurf_rl[5]+0.2500000000000001*alphaSurf_cr[1]*fSurf_cr[5]-1.0*jump_f_modal_r[4]+0.25*alphaSurf_rl[0]*fSurf_rl[4]+0.25*alphaSurf_cr[0]*fSurf_cr[4]; 
  Ghat_u_r[5] = (-1.0*jump_f_modal_r[5])+0.25*alphaSurf_rl[0]*fSurf_rl[5]+0.25*alphaSurf_cr[0]*fSurf_cr[5]+0.2500000000000001*alphaSurf_rl[1]*fSurf_rl[4]+0.2500000000000001*alphaSurf_cr[1]*fSurf_cr[4]; 

  out[0] += (0.7071067811865475*Ghat_u_l[0]-0.7071067811865475*Ghat_u_r[0])*dx1; 
  out[1] += -1.224744871391589*(Ghat_u_r[0]+Ghat_u_l[0])*dx1; 
  out[2] += (0.7071067811865475*Ghat_u_l[1]-0.7071067811865475*Ghat_u_r[1])*dx1; 
  out[3] += (0.7071067811865475*Ghat_u_l[2]-0.7071067811865475*Ghat_u_r[2])*dx1; 
  out[4] += -1.224744871391589*(Ghat_u_r[1]+Ghat_u_l[1])*dx1; 
  out[5] += -1.224744871391589*(Ghat_u_r[2]+Ghat_u_l[2])*dx1; 
  out[6] += (0.7071067811865475*Ghat_u_l[3]-0.7071067811865475*Ghat_u_r[3])*dx1; 
  out[7] += -1.224744871391589*(Ghat_u_r[3]+Ghat_u_l[3])*dx1; 
  out[8] += (0.7071067811865475*Ghat_u_l[4]-0.7071067811865475*Ghat_u_r[4])*dx1; 
  out[9] += -1.224744871391589*(Ghat_u_r[4]+Ghat_u_l[4])*dx1; 
  out[10] += (0.7071067811865475*Ghat_u_l[5]-0.7071067811865475*Ghat_u_r[5])*dx1; 
  out[11] += -1.224744871391589*(Ghat_u_r[5]+Ghat_u_l[5])*dx1; 


} 
