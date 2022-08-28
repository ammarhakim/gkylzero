#include <gkyl_mom_vlasov_kernels.h> 
#include <gkyl_basis_ser_2x_p2_surfx1_eval_quad.h> 
#include <gkyl_basis_ser_2x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH void mom_vlasov_pkpm_surfx_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, double mass, const double *u_il, const double *u_ir, const double *bvarl, const double *bvarr, const double *fl, const double *fr, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wvpar = w[1], dvpar = dxv[1]; 
  const double wvpar_sq = wvpar*wvpar, dvpar_sq = dvpar*dvpar; 
  const double wvpar_cu = wvpar*wvpar*wvpar, dvpar_cu = dvpar*dvpar*dvpar; 

  const double *ul_0 = &u_il[0]; 
  const double *ur_0 = &u_ir[0]; 
  const double *bl_0 = &bvarl[0]; 
  const double *br_0 = &bvarr[0]; 
  double u_l[3] = {0.0}; 
  double u_r[3] = {0.0}; 
  double q_l[3] = {0.0}; 
  double q_r[3] = {0.0}; 

  u_l[0] = 2.23606797749979*ul_0[2]+1.732050807568877*ul_0[1]+ul_0[0]; 

  u_r[0] = 2.23606797749979*ur_0[2]-1.732050807568877*ur_0[1]+ur_0[0]; 

  q_l[0] = 1.118033988749895*bl_0[2]*wvpar_cu+0.8660254037844386*bl_0[1]*wvpar_cu+0.5*bl_0[0]*wvpar_cu+0.2795084971874737*bl_0[2]*dvpar_sq*wvpar+0.2165063509461096*bl_0[1]*dvpar_sq*wvpar+0.125*bl_0[0]*dvpar_sq*wvpar; 
  q_l[1] = 0.9682458365518543*bl_0[2]*dvpar*wvpar_sq+0.75*bl_0[1]*dvpar*wvpar_sq+0.4330127018922193*bl_0[0]*dvpar*wvpar_sq+0.04841229182759271*bl_0[2]*dvpar_cu+0.0375*bl_0[1]*dvpar_cu+0.02165063509461097*bl_0[0]*dvpar_cu; 
  q_l[2] = 0.25*bl_0[2]*dvpar_sq*wvpar+0.1936491673103708*bl_0[1]*dvpar_sq*wvpar+0.1118033988749895*bl_0[0]*dvpar_sq*wvpar; 

  q_r[0] = 1.118033988749895*br_0[2]*wvpar_cu-0.8660254037844386*br_0[1]*wvpar_cu+0.5*br_0[0]*wvpar_cu+0.2795084971874737*br_0[2]*dvpar_sq*wvpar-0.2165063509461096*br_0[1]*dvpar_sq*wvpar+0.125*br_0[0]*dvpar_sq*wvpar; 
  q_r[1] = 0.9682458365518543*br_0[2]*dvpar*wvpar_sq-0.75*br_0[1]*dvpar*wvpar_sq+0.4330127018922193*br_0[0]*dvpar*wvpar_sq+0.04841229182759271*br_0[2]*dvpar_cu-0.0375*br_0[1]*dvpar_cu+0.02165063509461097*br_0[0]*dvpar_cu; 
  q_r[2] = 0.25*br_0[2]*dvpar_sq*wvpar-0.1936491673103708*br_0[1]*dvpar_sq*wvpar+0.1118033988749895*br_0[0]*dvpar_sq*wvpar; 

  double uQuad[3] = {0.0};
  double qQuad[3] = {0.0};
  double uMax[3] = {0.0};;
  double qMax[3] = {0.0};
  double ev_u_l = 0.0; 
  double ev_u_r = 0.0; 
  double ev_q_l = 0.0; 
  double ev_q_r = 0.0; 

  ev_u_l = ser_2x_p2_surfx1_eval_quad_node_0_r(u_l); 
  ev_u_r = ser_2x_p2_surfx1_eval_quad_node_0_l(u_r); 
  ev_q_l = ser_2x_p2_surfx1_eval_quad_node_0_r(q_l); 
  ev_q_r = ser_2x_p2_surfx1_eval_quad_node_0_l(q_r); 
  uQuad[0] = fmax(fabs(ev_u_l), fabs(ev_u_r)); 
  qQuad[0] = fmax(fabs(ev_q_l), fabs(ev_q_r)); 
  ev_u_l = ser_2x_p2_surfx1_eval_quad_node_1_r(u_l); 
  ev_u_r = ser_2x_p2_surfx1_eval_quad_node_1_l(u_r); 
  ev_q_l = ser_2x_p2_surfx1_eval_quad_node_1_r(q_l); 
  ev_q_r = ser_2x_p2_surfx1_eval_quad_node_1_l(q_r); 
  uQuad[1] = fmax(fabs(ev_u_l), fabs(ev_u_r)); 
  qQuad[1] = fmax(fabs(ev_q_l), fabs(ev_q_r)); 
  ev_u_l = ser_2x_p2_surfx1_eval_quad_node_2_r(u_l); 
  ev_u_r = ser_2x_p2_surfx1_eval_quad_node_2_l(u_r); 
  ev_q_l = ser_2x_p2_surfx1_eval_quad_node_2_r(q_l); 
  ev_q_r = ser_2x_p2_surfx1_eval_quad_node_2_l(q_r); 
  uQuad[2] = fmax(fabs(ev_u_l), fabs(ev_u_r)); 
  qQuad[2] = fmax(fabs(ev_q_l), fabs(ev_q_r)); 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_2x_p2_upwind_quad_to_modal(uQuad, uMax); 
  ser_2x_p2_upwind_quad_to_modal(qQuad, qMax); 
  out[0] += (0.6123724356957944*uMax[2]*fr[7]+0.6123724356957944*uMax[2]*fl[7]-0.7905694150420948*uMax[1]*fr[6]+0.7905694150420948*uMax[1]*fl[6]-0.3535533905932737*uMax[2]*fr[5]+0.3535533905932737*uMax[2]*fl[5]+0.7905694150420947*u_r[0]*fr[4]-0.7905694150420947*uMax[0]*fr[4]+0.7905694150420947*u_l[0]*fl[4]+0.7905694150420947*uMax[0]*fl[4]+0.6123724356957944*uMax[1]*fr[3]+0.6123724356957944*uMax[1]*fl[3]-0.3535533905932737*uMax[1]*fr[2]+0.3535533905932737*uMax[1]*fl[2]-0.6123724356957944*u_r[0]*fr[1]+0.6123724356957944*uMax[0]*fr[1]+0.6123724356957944*u_l[0]*fl[1]+0.6123724356957944*uMax[0]*fl[1]+0.3535533905932737*fr[0]*u_r[0]+0.3535533905932737*fl[0]*u_l[0]-0.3535533905932737*fr[0]*uMax[0]+0.3535533905932737*fl[0]*uMax[0])*mass*volFact; 
  out[1] += ((-0.6123724356957944*q_r[2]*fr[7])+0.6123724356957944*qMax[2]*fr[7]+0.6123724356957944*q_l[2]*fl[7]+0.6123724356957944*qMax[2]*fl[7]+0.7905694150420948*q_r[1]*fr[6]-0.7905694150420948*qMax[1]*fr[6]+0.7905694150420948*q_l[1]*fl[6]+0.7905694150420948*qMax[1]*fl[6]+0.3535533905932737*q_r[2]*fr[5]-0.3535533905932737*qMax[2]*fr[5]+0.3535533905932737*q_l[2]*fl[5]+0.3535533905932737*qMax[2]*fl[5]+0.7905694150420947*q_r[0]*fr[4]-0.7905694150420947*qMax[0]*fr[4]+0.7905694150420947*q_l[0]*fl[4]+0.7905694150420947*qMax[0]*fl[4]-0.6123724356957944*q_r[1]*fr[3]+0.6123724356957944*qMax[1]*fr[3]+0.6123724356957944*q_l[1]*fl[3]+0.6123724356957944*qMax[1]*fl[3]+0.3535533905932737*q_r[1]*fr[2]-0.3535533905932737*qMax[1]*fr[2]+0.3535533905932737*q_l[1]*fl[2]+0.3535533905932737*qMax[1]*fl[2]-0.6123724356957944*q_r[0]*fr[1]+0.6123724356957944*qMax[0]*fr[1]+0.6123724356957944*q_l[0]*fl[1]+0.6123724356957944*qMax[0]*fl[1]+0.3535533905932737*fr[0]*q_r[0]+0.3535533905932737*fl[0]*q_l[0]-0.3535533905932737*fr[0]*qMax[0]+0.3535533905932737*fl[0]*qMax[0])*mass*volFact; 
} 
