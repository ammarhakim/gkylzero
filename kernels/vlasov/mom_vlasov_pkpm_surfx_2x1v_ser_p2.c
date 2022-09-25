#include <gkyl_mom_vlasov_kernels.h> 
#include <gkyl_basis_ser_3x_p2_surfx1_eval_quad.h> 
#include <gkyl_basis_ser_3x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH void mom_vlasov_pkpm_surfx_2x1v_ser_p2(const double *w, const double *dxv, const int *idx, double mass, const double *u_i, const double *bvar, const double *fl, const double *fr, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[0]/2.0; 
  const double wvpar = w[0], dvpar = dxv[0]; 
  const double wvpar_sq = wvpar*wvpar, dvpar_sq = dvpar*dvpar; 
  const double wvpar_cu = wvpar*wvpar*wvpar, dvpar_cu = dvpar*dvpar*dvpar; 

  const double *ul = &u_i[0]; 
  const double *bl = &bvar[0]; 
  double *rho_flux = &out[0]; 
  double *heat_flux = &out[6]; 
  double u[8] = {0.0}; 
  double q[8] = {0.0}; 

  u[0] = 2.23606797749979*ul[4]+1.732050807568877*ul[1]+ul[0]; 
  u[1] = 2.23606797749979*ul[6]+1.732050807568877*ul[3]+ul[2]; 
  u[4] = 1.732050807568877*ul[7]+ul[5]; 

  q[0] = 1.118033988749895*bl[4]*wvpar_cu+0.8660254037844386*bl[1]*wvpar_cu+0.5*bl[0]*wvpar_cu+0.2795084971874737*bl[4]*dvpar_sq*wvpar+0.2165063509461096*bl[1]*dvpar_sq*wvpar+0.125*bl[0]*dvpar_sq*wvpar; 
  q[1] = 1.118033988749895*bl[6]*wvpar_cu+0.8660254037844386*bl[3]*wvpar_cu+0.5*bl[2]*wvpar_cu+0.2795084971874738*bl[6]*dvpar_sq*wvpar+0.2165063509461096*bl[3]*dvpar_sq*wvpar+0.125*bl[2]*dvpar_sq*wvpar; 
  q[2] = 0.9682458365518543*bl[4]*dvpar*wvpar_sq+0.75*bl[1]*dvpar*wvpar_sq+0.4330127018922193*bl[0]*dvpar*wvpar_sq+0.04841229182759271*bl[4]*dvpar_cu+0.0375*bl[1]*dvpar_cu+0.02165063509461097*bl[0]*dvpar_cu; 
  q[3] = 0.9682458365518543*bl[6]*dvpar*wvpar_sq+0.75*bl[3]*dvpar*wvpar_sq+0.4330127018922193*bl[2]*dvpar*wvpar_sq+0.04841229182759271*bl[6]*dvpar_cu+0.0375*bl[3]*dvpar_cu+0.02165063509461097*bl[2]*dvpar_cu; 
  q[4] = 0.8660254037844387*bl[7]*wvpar_cu+0.5*bl[5]*wvpar_cu+0.2165063509461097*bl[7]*dvpar_sq*wvpar+0.125*bl[5]*dvpar_sq*wvpar; 
  q[5] = 0.25*bl[4]*dvpar_sq*wvpar+0.1936491673103708*bl[1]*dvpar_sq*wvpar+0.1118033988749895*bl[0]*dvpar_sq*wvpar; 
  q[6] = 0.75*bl[7]*dvpar*wvpar_sq+0.4330127018922194*bl[5]*dvpar*wvpar_sq+0.0375*bl[7]*dvpar_cu+0.02165063509461096*bl[5]*dvpar_cu; 
  q[7] = 0.25*bl[6]*dvpar_sq*wvpar+0.1936491673103709*bl[3]*dvpar_sq*wvpar+0.1118033988749895*bl[2]*dvpar_sq*wvpar; 

  double fUpwindQuad_u[9] = {0.0};
  double fUpwindQuad_q[9] = {0.0};
  double fUpwind_u[8] = {0.0};;
  double fUpwind_q[8] = {0.0};
  double Ghat_u[8] = {0.0}; 
  double Ghat_q[8] = {0.0}; 

  if (0.4472135954999579*u[4]-0.6708203932499369*u[1]+0.5*u[0] > 0) { 
    fUpwindQuad_u[0] = ser_3x_p2_surfx1_eval_quad_node_0_r(fl); 
  } else { 
    fUpwindQuad_u[0] = ser_3x_p2_surfx1_eval_quad_node_0_l(fr); 
  } 
  if ((-0.5999999999999995*q[7])-0.5999999999999999*q[6]+0.4472135954999579*(q[5]+q[4])+0.9*q[3]-0.6708203932499369*(q[2]+q[1])+0.5*q[0] > 0) { 
    fUpwindQuad_q[0] = ser_3x_p2_surfx1_eval_quad_node_0_r(fl); 
  } else { 
    fUpwindQuad_q[0] = ser_3x_p2_surfx1_eval_quad_node_0_l(fr); 
  } 
  if (0.4472135954999579*u[4]-0.6708203932499369*u[1]+0.5*u[0] > 0) { 
    fUpwindQuad_u[1] = ser_3x_p2_surfx1_eval_quad_node_1_r(fl); 
  } else { 
    fUpwindQuad_u[1] = ser_3x_p2_surfx1_eval_quad_node_1_l(fr); 
  } 
  if (0.75*q[7]-0.5590169943749475*q[5]+0.4472135954999579*q[4]-0.6708203932499369*q[1]+0.5*q[0] > 0) { 
    fUpwindQuad_q[1] = ser_3x_p2_surfx1_eval_quad_node_1_r(fl); 
  } else { 
    fUpwindQuad_q[1] = ser_3x_p2_surfx1_eval_quad_node_1_l(fr); 
  } 
  if (0.4472135954999579*u[4]-0.6708203932499369*u[1]+0.5*u[0] > 0) { 
    fUpwindQuad_u[2] = ser_3x_p2_surfx1_eval_quad_node_2_r(fl); 
  } else { 
    fUpwindQuad_u[2] = ser_3x_p2_surfx1_eval_quad_node_2_l(fr); 
  } 
  if ((-0.5999999999999995*q[7])+0.5999999999999999*q[6]+0.4472135954999579*(q[5]+q[4])-0.9*q[3]+0.6708203932499369*q[2]-0.6708203932499369*q[1]+0.5*q[0] > 0) { 
    fUpwindQuad_q[2] = ser_3x_p2_surfx1_eval_quad_node_2_r(fl); 
  } else { 
    fUpwindQuad_q[2] = ser_3x_p2_surfx1_eval_quad_node_2_l(fr); 
  } 
  if (0.5*u[0]-0.5590169943749475*u[4] > 0) { 
    fUpwindQuad_u[3] = ser_3x_p2_surfx1_eval_quad_node_3_r(fl); 
  } else { 
    fUpwindQuad_u[3] = ser_3x_p2_surfx1_eval_quad_node_3_l(fr); 
  } 
  if (0.75*q[6]+0.4472135954999579*q[5]-0.5590169943749475*q[4]-0.6708203932499369*q[2]+0.5*q[0] > 0) { 
    fUpwindQuad_q[3] = ser_3x_p2_surfx1_eval_quad_node_3_r(fl); 
  } else { 
    fUpwindQuad_q[3] = ser_3x_p2_surfx1_eval_quad_node_3_l(fr); 
  } 
  if (0.5*u[0]-0.5590169943749475*u[4] > 0) { 
    fUpwindQuad_u[4] = ser_3x_p2_surfx1_eval_quad_node_4_r(fl); 
  } else { 
    fUpwindQuad_u[4] = ser_3x_p2_surfx1_eval_quad_node_4_l(fr); 
  } 
  if (0.5*q[0]-0.5590169943749475*(q[5]+q[4]) > 0) { 
    fUpwindQuad_q[4] = ser_3x_p2_surfx1_eval_quad_node_4_r(fl); 
  } else { 
    fUpwindQuad_q[4] = ser_3x_p2_surfx1_eval_quad_node_4_l(fr); 
  } 
  if (0.5*u[0]-0.5590169943749475*u[4] > 0) { 
    fUpwindQuad_u[5] = ser_3x_p2_surfx1_eval_quad_node_5_r(fl); 
  } else { 
    fUpwindQuad_u[5] = ser_3x_p2_surfx1_eval_quad_node_5_l(fr); 
  } 
  if ((-0.75*q[6])+0.4472135954999579*q[5]-0.5590169943749475*q[4]+0.6708203932499369*q[2]+0.5*q[0] > 0) { 
    fUpwindQuad_q[5] = ser_3x_p2_surfx1_eval_quad_node_5_r(fl); 
  } else { 
    fUpwindQuad_q[5] = ser_3x_p2_surfx1_eval_quad_node_5_l(fr); 
  } 
  if (0.4472135954999579*u[4]+0.6708203932499369*u[1]+0.5*u[0] > 0) { 
    fUpwindQuad_u[6] = ser_3x_p2_surfx1_eval_quad_node_6_r(fl); 
  } else { 
    fUpwindQuad_u[6] = ser_3x_p2_surfx1_eval_quad_node_6_l(fr); 
  } 
  if (0.5999999999999995*q[7]-0.5999999999999999*q[6]+0.4472135954999579*(q[5]+q[4])-0.9*q[3]-0.6708203932499369*q[2]+0.6708203932499369*q[1]+0.5*q[0] > 0) { 
    fUpwindQuad_q[6] = ser_3x_p2_surfx1_eval_quad_node_6_r(fl); 
  } else { 
    fUpwindQuad_q[6] = ser_3x_p2_surfx1_eval_quad_node_6_l(fr); 
  } 
  if (0.4472135954999579*u[4]+0.6708203932499369*u[1]+0.5*u[0] > 0) { 
    fUpwindQuad_u[7] = ser_3x_p2_surfx1_eval_quad_node_7_r(fl); 
  } else { 
    fUpwindQuad_u[7] = ser_3x_p2_surfx1_eval_quad_node_7_l(fr); 
  } 
  if ((-0.75*q[7])-0.5590169943749475*q[5]+0.4472135954999579*q[4]+0.6708203932499369*q[1]+0.5*q[0] > 0) { 
    fUpwindQuad_q[7] = ser_3x_p2_surfx1_eval_quad_node_7_r(fl); 
  } else { 
    fUpwindQuad_q[7] = ser_3x_p2_surfx1_eval_quad_node_7_l(fr); 
  } 
  if (0.4472135954999579*u[4]+0.6708203932499369*u[1]+0.5*u[0] > 0) { 
    fUpwindQuad_u[8] = ser_3x_p2_surfx1_eval_quad_node_8_r(fl); 
  } else { 
    fUpwindQuad_u[8] = ser_3x_p2_surfx1_eval_quad_node_8_l(fr); 
  } 
  if (0.5999999999999995*q[7]+0.5999999999999999*q[6]+0.4472135954999579*(q[5]+q[4])+0.9*q[3]+0.6708203932499369*(q[2]+q[1])+0.5*q[0] > 0) { 
    fUpwindQuad_q[8] = ser_3x_p2_surfx1_eval_quad_node_8_r(fl); 
  } else { 
    fUpwindQuad_q[8] = ser_3x_p2_surfx1_eval_quad_node_8_l(fr); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_3x_p2_upwind_quad_to_modal(fUpwindQuad_u, fUpwind_u); 
  ser_3x_p2_upwind_quad_to_modal(fUpwindQuad_q, fUpwind_q); 

  rho_flux[0] += (0.7071067811865475*fUpwind_u[4]*u[4]+0.7071067811865475*fUpwind_u[1]*u[1]+0.7071067811865475*fUpwind_u[0]*u[0])*mass*volFact; 
  rho_flux[1] += (0.6324555320336759*fUpwind_u[1]*u[4]+0.6324555320336759*u[1]*fUpwind_u[4]+0.7071067811865475*fUpwind_u[0]*u[1]+0.7071067811865475*u[0]*fUpwind_u[1])*mass*volFact; 
  rho_flux[2] += (0.4517539514526256*fUpwind_u[4]*u[4]+0.7071067811865475*fUpwind_u[0]*u[4]+0.7071067811865475*u[0]*fUpwind_u[4]+0.6324555320336759*fUpwind_u[1]*u[1])*mass*volFact; 
  heat_flux[0] += (0.7071067811865475*fUpwind_q[7]*q[7]+0.7071067811865475*fUpwind_q[6]*q[6]+0.7071067811865475*fUpwind_q[5]*q[5]+0.7071067811865475*fUpwind_q[4]*q[4]+0.7071067811865475*fUpwind_q[3]*q[3]+0.7071067811865475*fUpwind_q[2]*q[2]+0.7071067811865475*fUpwind_q[1]*q[1]+0.7071067811865475*fUpwind_q[0]*q[0])*mass*volFact; 
  heat_flux[1] += (0.7071067811865475*fUpwind_q[5]*q[7]+0.7071067811865475*q[5]*fUpwind_q[7]+0.632455532033676*fUpwind_q[3]*q[6]+0.632455532033676*q[3]*fUpwind_q[6]+0.6324555320336759*fUpwind_q[1]*q[4]+0.6324555320336759*q[1]*fUpwind_q[4]+0.7071067811865475*fUpwind_q[2]*q[3]+0.7071067811865475*q[2]*fUpwind_q[3]+0.7071067811865475*fUpwind_q[0]*q[1]+0.7071067811865475*q[0]*fUpwind_q[1])*mass*volFact; 
  heat_flux[2] += (0.6324555320336759*fUpwind_q[7]*q[7]+0.4517539514526256*fUpwind_q[6]*q[6]+0.7071067811865475*fUpwind_q[2]*q[6]+0.7071067811865475*q[2]*fUpwind_q[6]+0.4517539514526256*fUpwind_q[4]*q[4]+0.7071067811865475*fUpwind_q[0]*q[4]+0.7071067811865475*q[0]*fUpwind_q[4]+0.6324555320336759*fUpwind_q[3]*q[3]+0.6324555320336759*fUpwind_q[1]*q[1])*mass*volFact; 
} 
