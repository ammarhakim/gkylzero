#include <gkyl_mom_vlasov_kernels.h> 
#include <gkyl_basis_ser_2x_p3_surfx1_eval_quad.h> 
#include <gkyl_basis_ser_2x_p3_upwind_quad_to_modal.h> 
GKYL_CU_DH void mom_vlasov_pkpm_surfx_1x1v_ser_p3(const double *w, const double *dxv, const int *idx, double mass, const double *u_i, const double *bvar, const double *fl, const double *fr, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[0]/2.0; 
  const double wvpar = w[0], dvpar = dxv[0]; 
  const double wvpar_sq = wvpar*wvpar, dvpar_sq = dvpar*dvpar; 
  const double wvpar_cu = wvpar*wvpar*wvpar, dvpar_cu = dvpar*dvpar*dvpar; 

  const double *ul = &u_i[0]; 
  const double *bl = &bvar[0]; 
  double *rho_flux = &out[0]; 
  double *heat_flux = &out[1]; 
  double u[4] = {0.0}; 
  double q[4] = {0.0}; 

  u[0] = 2.645751311064591*ul[3]+2.23606797749979*ul[2]+1.732050807568877*ul[1]+ul[0]; 

  q[0] = 1.322875655532295*bl[3]*wvpar_cu+1.118033988749895*bl[2]*wvpar_cu+0.8660254037844386*bl[1]*wvpar_cu+0.5*bl[0]*wvpar_cu+0.3307189138830738*bl[3]*dvpar_sq*wvpar+0.2795084971874737*bl[2]*dvpar_sq*wvpar+0.2165063509461096*bl[1]*dvpar_sq*wvpar+0.125*bl[0]*dvpar_sq*wvpar; 
  q[1] = 1.14564392373896*bl[3]*dvpar*wvpar_sq+0.9682458365518543*bl[2]*dvpar*wvpar_sq+0.75*bl[1]*dvpar*wvpar_sq+0.4330127018922193*bl[0]*dvpar*wvpar_sq+0.05728219618694801*bl[3]*dvpar_cu+0.04841229182759271*bl[2]*dvpar_cu+0.0375*bl[1]*dvpar_cu+0.02165063509461097*bl[0]*dvpar_cu; 
  q[2] = 0.2958039891549808*bl[3]*dvpar_sq*wvpar+0.25*bl[2]*dvpar_sq*wvpar+0.1936491673103708*bl[1]*dvpar_sq*wvpar+0.1118033988749895*bl[0]*dvpar_sq*wvpar; 
  q[3] = 0.025*bl[3]*dvpar_cu+0.02112885636821291*bl[2]*dvpar_cu+0.01636634176769943*bl[1]*dvpar_cu+0.00944911182523068*bl[0]*dvpar_cu; 

  double fUpwindQuad_u[4] = {0.0};
  double fUpwindQuad_q[4] = {0.0};
  double fUpwind_u[4] = {0.0};;
  double fUpwind_q[4] = {0.0};
  double Ghat_u[4] = {0.0}; 
  double Ghat_q[4] = {0.0}; 

  if (0.7071067811865475*u[0] > 0) { 
    fUpwindQuad_u[0] = ser_2x_p3_surfx1_eval_quad_node_0_r(fl); 
  } else { 
    fUpwindQuad_u[0] = ser_2x_p3_surfx1_eval_quad_node_0_l(fr); 
  } 
  if ((-0.5701294036773671*q[3])+0.9681844646844029*q[2]-1.054672281193885*q[1]+0.7071067811865475*q[0] > 0) { 
    fUpwindQuad_q[0] = ser_2x_p3_surfx1_eval_quad_node_0_r(fl); 
  } else { 
    fUpwindQuad_q[0] = ser_2x_p3_surfx1_eval_quad_node_0_l(fr); 
  } 
  if (0.7071067811865475*u[0] > 0) { 
    fUpwindQuad_u[1] = ser_2x_p3_surfx1_eval_quad_node_1_r(fl); 
  } else { 
    fUpwindQuad_u[1] = ser_2x_p3_surfx1_eval_quad_node_1_l(fr); 
  } 
  if (0.7702725556588816*q[3]-0.5164305132317774*q[2]-0.4163900395009129*q[1]+0.7071067811865475*q[0] > 0) { 
    fUpwindQuad_q[1] = ser_2x_p3_surfx1_eval_quad_node_1_r(fl); 
  } else { 
    fUpwindQuad_q[1] = ser_2x_p3_surfx1_eval_quad_node_1_l(fr); 
  } 
  if (0.7071067811865475*u[0] > 0) { 
    fUpwindQuad_u[2] = ser_2x_p3_surfx1_eval_quad_node_2_r(fl); 
  } else { 
    fUpwindQuad_u[2] = ser_2x_p3_surfx1_eval_quad_node_2_l(fr); 
  } 
  if ((-0.7702725556588816*q[3])-0.5164305132317774*q[2]+0.4163900395009129*q[1]+0.7071067811865475*q[0] > 0) { 
    fUpwindQuad_q[2] = ser_2x_p3_surfx1_eval_quad_node_2_r(fl); 
  } else { 
    fUpwindQuad_q[2] = ser_2x_p3_surfx1_eval_quad_node_2_l(fr); 
  } 
  if (0.7071067811865475*u[0] > 0) { 
    fUpwindQuad_u[3] = ser_2x_p3_surfx1_eval_quad_node_3_r(fl); 
  } else { 
    fUpwindQuad_u[3] = ser_2x_p3_surfx1_eval_quad_node_3_l(fr); 
  } 
  if (0.5701294036773671*q[3]+0.9681844646844029*q[2]+1.054672281193885*q[1]+0.7071067811865475*q[0] > 0) { 
    fUpwindQuad_q[3] = ser_2x_p3_surfx1_eval_quad_node_3_r(fl); 
  } else { 
    fUpwindQuad_q[3] = ser_2x_p3_surfx1_eval_quad_node_3_l(fr); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_2x_p3_upwind_quad_to_modal(fUpwindQuad_u, fUpwind_u); 
  ser_2x_p3_upwind_quad_to_modal(fUpwindQuad_q, fUpwind_q); 

  rho_flux[0] += fUpwind_u[0]*u[0]*mass*volFact; 
  heat_flux[0] += (fUpwind_q[3]*q[3]+fUpwind_q[2]*q[2]+fUpwind_q[1]*q[1]+fUpwind_q[0]*q[0])*mass*volFact; 
} 
