#include <gkyl_isoeuler_kernels.h> 
#include <gkyl_basis_ser_2x_p1_upwind_quad_to_modal.h> 
#include <gkyl_basis_ser_2x_p1_surfx1_eval_quad.h> 
GKYL_CU_DH void isoeuler_surfx_2x_ser_p1(const double *w, const double *dxv, double vth, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out)  
{ 
  // w: Cell-center coordinates. dxv[NDIM]: Cell spacing. statevec(l/c/r): [rho, rho ux, rho uy, rho uz] in (left/center/right) cell, uvar(l/c/r): [ux, uv, uz]  in (left/center/right) cell
  // out: output in center cell 
  const double dx1 = 2.0/dxv[0]; 
  const double *rhol = &statevecl[0]; 
  const double *rhou0l = &statevecl[4]; 
  const double *rhou1l = &statevecl[8]; 
  const double *rhou2l = &statevecl[12]; 
  const double *rhoc = &statevecc[0]; 
  const double *rhou0c = &statevecc[4]; 
  const double *rhou1c = &statevecc[8]; 
  const double *rhou2c = &statevecc[12]; 
  const double *rhor = &statevecr[0]; 
  const double *rhou0r = &statevecr[4]; 
  const double *rhou1r = &statevecr[8]; 
  const double *rhou2r = &statevecr[12]; 
  const double *uvar0l = &uvarl[0]; 
  const double *uvar1l = &uvarl[4]; 
  const double *uvar2l = &uvarl[8]; 
  const double *uvar0c = &uvarc[0]; 
  const double *uvar1c = &uvarc[4]; 
  const double *uvar2c = &uvarc[8]; 
  const double *uvar0r = &uvarr[0]; 
  const double *uvar1r = &uvarr[4]; 
  const double *uvar2r = &uvarr[8]; 
  double *outrho = &out[0]; 
  double *outrhoux = &out[4]; 
  double *outrhouy = &out[8]; 
  double *outrhouz = &out[12]; 
  double incr[4]; 

  double vthsq = vth*vth; 
  double umax_l[2] = {0.0};
  double umax_r[2] = {0.0};
  double uquad_l[2] = {0.0};
  double uquad_r[2] = {0.0};
  double urec_0_l[2] = {0.0}; 
  double urec_0_r[2] = {0.0}; 
  double urec_1_l[2] = {0.0}; 
  double urec_1_r[2] = {0.0}; 
  double urec_2_l[2] = {0.0}; 
  double urec_2_r[2] = {0.0}; 
  double rhourec_0_l[2] = {0.0}; 
  double rhourec_0_r[2] = {0.0}; 
  double rhourec_1_l[2] = {0.0}; 
  double rhourec_1_r[2] = {0.0}; 
  double rhourec_2_l[2] = {0.0}; 
  double rhourec_2_r[2] = {0.0}; 
  double rhorec_l[2] = {0.0}; 
  double rhorec_r[2] = {0.0}; 
  double ghat_rho_l[2] = {0.0}; 
  double ghat_rho_r[2] = {0.0}; 
  double u_l_r = 0.0; 
  double u_c_l = 0.0; 
  double u_c_r = 0.0; 
  double u_r_l = 0.0; 

  u_l_r = ser_2x_p1_surfx1_eval_quad_node_0_r(uvar0l); 
  u_c_l = ser_2x_p1_surfx1_eval_quad_node_0_l(uvar0c); 
  u_c_r = ser_2x_p1_surfx1_eval_quad_node_0_r(uvar0c); 
  u_r_l = ser_2x_p1_surfx1_eval_quad_node_0_l(uvar0r); 
  uquad_l[0] = fmax(fabs(u_l_r), fabs(u_c_l)); 
  uquad_r[0] = fmax(fabs(u_c_r), fabs(u_r_l)); 
  u_l_r = ser_2x_p1_surfx1_eval_quad_node_1_r(uvar0l); 
  u_c_l = ser_2x_p1_surfx1_eval_quad_node_1_l(uvar0c); 
  u_c_r = ser_2x_p1_surfx1_eval_quad_node_1_r(uvar0c); 
  u_r_l = ser_2x_p1_surfx1_eval_quad_node_1_l(uvar0r); 
  uquad_l[1] = fmax(fabs(u_l_r), fabs(u_c_l)); 
  uquad_r[1] = fmax(fabs(u_c_r), fabs(u_r_l)); 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_2x_p1_upwind_quad_to_modal(uquad_l, umax_l); 
  ser_2x_p1_upwind_quad_to_modal(uquad_r, umax_r); 

  urec_0_l[0] = 0.408248290463863*uvar0l[1]-0.408248290463863*uvar0c[1]+0.3535533905932737*uvar0l[0]+0.3535533905932737*uvar0c[0]; 
  urec_0_l[1] = 0.408248290463863*uvar0l[3]-0.408248290463863*uvar0c[3]+0.3535533905932737*uvar0l[2]+0.3535533905932737*uvar0c[2]; 

  urec_1_l[0] = 0.408248290463863*uvar1l[1]-0.408248290463863*uvar1c[1]+0.3535533905932737*uvar1l[0]+0.3535533905932737*uvar1c[0]; 
  urec_1_l[1] = 0.408248290463863*uvar1l[3]-0.408248290463863*uvar1c[3]+0.3535533905932737*uvar1l[2]+0.3535533905932737*uvar1c[2]; 

  urec_2_l[0] = 0.408248290463863*uvar2l[1]-0.408248290463863*uvar2c[1]+0.3535533905932737*uvar2l[0]+0.3535533905932737*uvar2c[0]; 
  urec_2_l[1] = 0.408248290463863*uvar2l[3]-0.408248290463863*uvar2c[3]+0.3535533905932737*uvar2l[2]+0.3535533905932737*uvar2c[2]; 

  rhourec_0_l[0] = 0.408248290463863*rhou0l[1]-0.408248290463863*rhou0c[1]+0.3535533905932737*rhou0l[0]+0.3535533905932737*rhou0c[0]; 
  rhourec_0_l[1] = 0.408248290463863*rhou0l[3]-0.408248290463863*rhou0c[3]+0.3535533905932737*rhou0l[2]+0.3535533905932737*rhou0c[2]; 

  rhourec_1_l[0] = 0.408248290463863*rhou1l[1]-0.408248290463863*rhou1c[1]+0.3535533905932737*rhou1l[0]+0.3535533905932737*rhou1c[0]; 
  rhourec_1_l[1] = 0.408248290463863*rhou1l[3]-0.408248290463863*rhou1c[3]+0.3535533905932737*rhou1l[2]+0.3535533905932737*rhou1c[2]; 

  rhourec_2_l[0] = 0.408248290463863*rhou2l[1]-0.408248290463863*rhou2c[1]+0.3535533905932737*rhou2l[0]+0.3535533905932737*rhou2c[0]; 
  rhourec_2_l[1] = 0.408248290463863*rhou2l[3]-0.408248290463863*rhou2c[3]+0.3535533905932737*rhou2l[2]+0.3535533905932737*rhou2c[2]; 

  rhorec_l[0] = 0.408248290463863*rhol[1]-0.408248290463863*rhoc[1]+0.3535533905932737*rhol[0]+0.3535533905932737*rhoc[0]; 
  rhorec_l[1] = 0.408248290463863*rhol[3]-0.408248290463863*rhoc[3]+0.3535533905932737*rhol[2]+0.3535533905932737*rhoc[2]; 

  urec_0_r[0] = (-0.408248290463863*uvar0r[1])+0.408248290463863*uvar0c[1]+0.3535533905932737*uvar0r[0]+0.3535533905932737*uvar0c[0]; 
  urec_0_r[1] = (-0.408248290463863*uvar0r[3])+0.408248290463863*uvar0c[3]+0.3535533905932737*uvar0r[2]+0.3535533905932737*uvar0c[2]; 

  urec_1_r[0] = (-0.408248290463863*uvar1r[1])+0.408248290463863*uvar1c[1]+0.3535533905932737*uvar1r[0]+0.3535533905932737*uvar1c[0]; 
  urec_1_r[1] = (-0.408248290463863*uvar1r[3])+0.408248290463863*uvar1c[3]+0.3535533905932737*uvar1r[2]+0.3535533905932737*uvar1c[2]; 

  urec_2_r[0] = (-0.408248290463863*uvar2l[1])+0.408248290463863*uvar2c[1]+0.3535533905932737*uvar2l[0]+0.3535533905932737*uvar2c[0]; 
  urec_2_r[1] = (-0.408248290463863*uvar2l[3])+0.408248290463863*uvar2c[3]+0.3535533905932737*uvar2l[2]+0.3535533905932737*uvar2c[2]; 

  rhourec_0_r[0] = (-0.408248290463863*rhou0r[1])+0.408248290463863*rhou0c[1]+0.3535533905932737*rhou0r[0]+0.3535533905932737*rhou0c[0]; 
  rhourec_0_r[1] = (-0.408248290463863*rhou0r[3])+0.408248290463863*rhou0c[3]+0.3535533905932737*rhou0r[2]+0.3535533905932737*rhou0c[2]; 

  rhourec_1_r[0] = (-0.408248290463863*rhou1r[1])+0.408248290463863*rhou1c[1]+0.3535533905932737*rhou1r[0]+0.3535533905932737*rhou1c[0]; 
  rhourec_1_r[1] = (-0.408248290463863*rhou1r[3])+0.408248290463863*rhou1c[3]+0.3535533905932737*rhou1r[2]+0.3535533905932737*rhou1c[2]; 

  rhourec_2_r[0] = (-0.408248290463863*rhou2r[1])+0.408248290463863*rhou2c[1]+0.3535533905932737*rhou2r[0]+0.3535533905932737*rhou2c[0]; 
  rhourec_2_r[1] = (-0.408248290463863*rhou2r[3])+0.408248290463863*rhou2c[3]+0.3535533905932737*rhou2r[2]+0.3535533905932737*rhou2c[2]; 

  rhorec_r[0] = (-0.408248290463863*rhor[1])+0.408248290463863*rhoc[1]+0.3535533905932737*rhor[0]+0.3535533905932737*rhoc[0]; 
  rhorec_r[1] = (-0.408248290463863*rhor[3])+0.408248290463863*rhoc[3]+0.3535533905932737*rhor[2]+0.3535533905932737*rhoc[2]; 

  ghat_rho_l[0] = 0.4330127018922193*umax_l[1]*rhol[3]+0.4330127018922193*umax_l[1]*rhoc[3]+0.25*umax_l[1]*rhol[2]-0.25*umax_l[1]*rhoc[2]+0.4330127018922193*umax_l[0]*rhol[1]+0.4330127018922193*umax_l[0]*rhoc[1]+0.25*rhol[0]*umax_l[0]-0.25*rhoc[0]*umax_l[0]+rhourec_0_l[0]; 
  ghat_rho_l[1] = 0.4330127018922193*umax_l[0]*rhol[3]+0.4330127018922193*umax_l[0]*rhoc[3]+0.25*umax_l[0]*rhol[2]-0.25*umax_l[0]*rhoc[2]+0.4330127018922193*rhol[1]*umax_l[1]+0.4330127018922193*rhoc[1]*umax_l[1]+0.25*rhol[0]*umax_l[1]-0.25*rhoc[0]*umax_l[1]+rhourec_0_l[1]; 

  ghat_rho_r[0] = 0.4330127018922193*umax_r[1]*rhor[3]+0.4330127018922193*umax_r[1]*rhoc[3]-0.25*umax_r[1]*rhor[2]+0.25*umax_r[1]*rhoc[2]+0.4330127018922193*umax_r[0]*rhor[1]+0.4330127018922193*umax_r[0]*rhoc[1]-0.25*rhor[0]*umax_r[0]+0.25*rhoc[0]*umax_r[0]+rhourec_0_r[0]; 
  ghat_rho_r[1] = 0.4330127018922193*umax_r[0]*rhor[3]+0.4330127018922193*umax_r[0]*rhoc[3]-0.25*umax_r[0]*rhor[2]+0.25*umax_r[0]*rhoc[2]+0.4330127018922193*rhor[1]*umax_r[1]+0.4330127018922193*rhoc[1]*umax_r[1]-0.25*rhor[0]*umax_r[1]+0.25*rhoc[0]*umax_r[1]+rhourec_0_r[1]; 

  outrho[0] += 0.7071067811865475*ghat_rho_l[0]*dx1-0.7071067811865475*ghat_rho_r[0]*dx1; 
  outrho[1] += (-1.224744871391589*ghat_rho_r[0]*dx1)-1.224744871391589*ghat_rho_l[0]*dx1; 
  outrho[2] += 0.7071067811865475*ghat_rho_l[1]*dx1-0.7071067811865475*ghat_rho_r[1]*dx1; 
  outrho[3] += (-1.224744871391589*ghat_rho_r[1]*dx1)-1.224744871391589*ghat_rho_l[1]*dx1; 

  outrhoux[0] += (-0.7071067811865475*rhorec_r[0]*dx1*vthsq)+0.7071067811865475*rhorec_l[0]*dx1*vthsq-0.3061862178478971*umax_r[1]*rhou0r[3]*dx1+0.3061862178478971*umax_l[1]*rhou0l[3]*dx1-0.3061862178478971*umax_r[1]*rhou0c[3]*dx1+0.3061862178478971*umax_l[1]*rhou0c[3]*dx1+0.1767766952966368*umax_r[1]*rhou0r[2]*dx1+0.1767766952966368*umax_l[1]*rhou0l[2]*dx1-0.1767766952966368*umax_r[1]*rhou0c[2]*dx1-0.1767766952966368*umax_l[1]*rhou0c[2]*dx1-0.5*ghat_rho_r[1]*urec_0_r[1]*dx1+0.5*ghat_rho_l[1]*urec_0_l[1]*dx1-0.3061862178478971*umax_r[0]*rhou0r[1]*dx1+0.3061862178478971*umax_l[0]*rhou0l[1]*dx1-0.3061862178478971*umax_r[0]*rhou0c[1]*dx1+0.3061862178478971*umax_l[0]*rhou0c[1]*dx1-0.5*ghat_rho_r[0]*urec_0_r[0]*dx1+0.5*ghat_rho_l[0]*urec_0_l[0]*dx1+0.1767766952966368*rhou0r[0]*umax_r[0]*dx1-0.1767766952966368*rhou0c[0]*umax_r[0]*dx1+0.1767766952966368*rhou0l[0]*umax_l[0]*dx1-0.1767766952966368*rhou0c[0]*umax_l[0]*dx1; 
  outrhoux[1] += (-1.224744871391589*rhorec_r[0]*dx1*vthsq)-1.224744871391589*rhorec_l[0]*dx1*vthsq-0.5303300858899105*umax_r[1]*rhou0r[3]*dx1-0.5303300858899105*umax_l[1]*rhou0l[3]*dx1-0.5303300858899105*umax_r[1]*rhou0c[3]*dx1-0.5303300858899105*umax_l[1]*rhou0c[3]*dx1+0.3061862178478971*umax_r[1]*rhou0r[2]*dx1-0.3061862178478971*umax_l[1]*rhou0l[2]*dx1-0.3061862178478971*umax_r[1]*rhou0c[2]*dx1+0.3061862178478971*umax_l[1]*rhou0c[2]*dx1-0.8660254037844386*ghat_rho_r[1]*urec_0_r[1]*dx1-0.8660254037844386*ghat_rho_l[1]*urec_0_l[1]*dx1-0.5303300858899105*umax_r[0]*rhou0r[1]*dx1-0.5303300858899105*umax_l[0]*rhou0l[1]*dx1-0.5303300858899105*umax_r[0]*rhou0c[1]*dx1-0.5303300858899105*umax_l[0]*rhou0c[1]*dx1-0.8660254037844386*ghat_rho_r[0]*urec_0_r[0]*dx1-0.8660254037844386*ghat_rho_l[0]*urec_0_l[0]*dx1+0.3061862178478971*rhou0r[0]*umax_r[0]*dx1-0.3061862178478971*rhou0c[0]*umax_r[0]*dx1-0.3061862178478971*rhou0l[0]*umax_l[0]*dx1+0.3061862178478971*rhou0c[0]*umax_l[0]*dx1; 
  outrhoux[2] += (-0.7071067811865475*rhorec_r[1]*dx1*vthsq)+0.7071067811865475*rhorec_l[1]*dx1*vthsq-0.3061862178478971*umax_r[0]*rhou0r[3]*dx1+0.3061862178478971*umax_l[0]*rhou0l[3]*dx1-0.3061862178478971*umax_r[0]*rhou0c[3]*dx1+0.3061862178478971*umax_l[0]*rhou0c[3]*dx1+0.1767766952966368*umax_r[0]*rhou0r[2]*dx1+0.1767766952966368*umax_l[0]*rhou0l[2]*dx1-0.1767766952966368*umax_r[0]*rhou0c[2]*dx1-0.1767766952966368*umax_l[0]*rhou0c[2]*dx1-0.5*ghat_rho_r[0]*urec_0_r[1]*dx1+0.5*ghat_rho_l[0]*urec_0_l[1]*dx1-0.3061862178478971*rhou0r[1]*umax_r[1]*dx1-0.3061862178478971*rhou0c[1]*umax_r[1]*dx1+0.1767766952966368*rhou0r[0]*umax_r[1]*dx1-0.1767766952966368*rhou0c[0]*umax_r[1]*dx1+0.3061862178478971*rhou0l[1]*umax_l[1]*dx1+0.3061862178478971*rhou0c[1]*umax_l[1]*dx1+0.1767766952966368*rhou0l[0]*umax_l[1]*dx1-0.1767766952966368*rhou0c[0]*umax_l[1]*dx1-0.5*urec_0_r[0]*ghat_rho_r[1]*dx1+0.5*urec_0_l[0]*ghat_rho_l[1]*dx1; 
  outrhoux[3] += (-1.224744871391589*rhorec_r[1]*dx1*vthsq)-1.224744871391589*rhorec_l[1]*dx1*vthsq-0.5303300858899105*umax_r[0]*rhou0r[3]*dx1-0.5303300858899105*umax_l[0]*rhou0l[3]*dx1-0.5303300858899105*umax_r[0]*rhou0c[3]*dx1-0.5303300858899105*umax_l[0]*rhou0c[3]*dx1+0.3061862178478971*umax_r[0]*rhou0r[2]*dx1-0.3061862178478971*umax_l[0]*rhou0l[2]*dx1-0.3061862178478971*umax_r[0]*rhou0c[2]*dx1+0.3061862178478971*umax_l[0]*rhou0c[2]*dx1-0.8660254037844386*ghat_rho_r[0]*urec_0_r[1]*dx1-0.8660254037844386*ghat_rho_l[0]*urec_0_l[1]*dx1-0.5303300858899105*rhou0r[1]*umax_r[1]*dx1-0.5303300858899105*rhou0c[1]*umax_r[1]*dx1+0.3061862178478971*rhou0r[0]*umax_r[1]*dx1-0.3061862178478971*rhou0c[0]*umax_r[1]*dx1-0.5303300858899105*rhou0l[1]*umax_l[1]*dx1-0.5303300858899105*rhou0c[1]*umax_l[1]*dx1-0.3061862178478971*rhou0l[0]*umax_l[1]*dx1+0.3061862178478971*rhou0c[0]*umax_l[1]*dx1-0.8660254037844386*urec_0_r[0]*ghat_rho_r[1]*dx1-0.8660254037844386*urec_0_l[0]*ghat_rho_l[1]*dx1; 

  outrhouy[0] += (-0.3061862178478971*umax_r[1]*rhou1r[3]*dx1)+0.3061862178478971*umax_l[1]*rhou1l[3]*dx1-0.3061862178478971*umax_r[1]*rhou1c[3]*dx1+0.3061862178478971*umax_l[1]*rhou1c[3]*dx1+0.1767766952966368*umax_r[1]*rhou1r[2]*dx1+0.1767766952966368*umax_l[1]*rhou1l[2]*dx1-0.1767766952966368*umax_r[1]*rhou1c[2]*dx1-0.1767766952966368*umax_l[1]*rhou1c[2]*dx1-0.5*ghat_rho_r[1]*urec_1_r[1]*dx1+0.5*ghat_rho_l[1]*urec_1_l[1]*dx1-0.3061862178478971*umax_r[0]*rhou1r[1]*dx1+0.3061862178478971*umax_l[0]*rhou1l[1]*dx1-0.3061862178478971*umax_r[0]*rhou1c[1]*dx1+0.3061862178478971*umax_l[0]*rhou1c[1]*dx1-0.5*ghat_rho_r[0]*urec_1_r[0]*dx1+0.5*ghat_rho_l[0]*urec_1_l[0]*dx1+0.1767766952966368*rhou1r[0]*umax_r[0]*dx1-0.1767766952966368*rhou1c[0]*umax_r[0]*dx1+0.1767766952966368*rhou1l[0]*umax_l[0]*dx1-0.1767766952966368*rhou1c[0]*umax_l[0]*dx1; 
  outrhouy[1] += (-0.5303300858899105*umax_r[1]*rhou1r[3]*dx1)-0.5303300858899105*umax_l[1]*rhou1l[3]*dx1-0.5303300858899105*umax_r[1]*rhou1c[3]*dx1-0.5303300858899105*umax_l[1]*rhou1c[3]*dx1+0.3061862178478971*umax_r[1]*rhou1r[2]*dx1-0.3061862178478971*umax_l[1]*rhou1l[2]*dx1-0.3061862178478971*umax_r[1]*rhou1c[2]*dx1+0.3061862178478971*umax_l[1]*rhou1c[2]*dx1-0.8660254037844386*ghat_rho_r[1]*urec_1_r[1]*dx1-0.8660254037844386*ghat_rho_l[1]*urec_1_l[1]*dx1-0.5303300858899105*umax_r[0]*rhou1r[1]*dx1-0.5303300858899105*umax_l[0]*rhou1l[1]*dx1-0.5303300858899105*umax_r[0]*rhou1c[1]*dx1-0.5303300858899105*umax_l[0]*rhou1c[1]*dx1-0.8660254037844386*ghat_rho_r[0]*urec_1_r[0]*dx1-0.8660254037844386*ghat_rho_l[0]*urec_1_l[0]*dx1+0.3061862178478971*rhou1r[0]*umax_r[0]*dx1-0.3061862178478971*rhou1c[0]*umax_r[0]*dx1-0.3061862178478971*rhou1l[0]*umax_l[0]*dx1+0.3061862178478971*rhou1c[0]*umax_l[0]*dx1; 
  outrhouy[2] += (-0.3061862178478971*umax_r[0]*rhou1r[3]*dx1)+0.3061862178478971*umax_l[0]*rhou1l[3]*dx1-0.3061862178478971*umax_r[0]*rhou1c[3]*dx1+0.3061862178478971*umax_l[0]*rhou1c[3]*dx1+0.1767766952966368*umax_r[0]*rhou1r[2]*dx1+0.1767766952966368*umax_l[0]*rhou1l[2]*dx1-0.1767766952966368*umax_r[0]*rhou1c[2]*dx1-0.1767766952966368*umax_l[0]*rhou1c[2]*dx1-0.5*ghat_rho_r[0]*urec_1_r[1]*dx1+0.5*ghat_rho_l[0]*urec_1_l[1]*dx1-0.3061862178478971*rhou1r[1]*umax_r[1]*dx1-0.3061862178478971*rhou1c[1]*umax_r[1]*dx1+0.1767766952966368*rhou1r[0]*umax_r[1]*dx1-0.1767766952966368*rhou1c[0]*umax_r[1]*dx1+0.3061862178478971*rhou1l[1]*umax_l[1]*dx1+0.3061862178478971*rhou1c[1]*umax_l[1]*dx1+0.1767766952966368*rhou1l[0]*umax_l[1]*dx1-0.1767766952966368*rhou1c[0]*umax_l[1]*dx1-0.5*urec_1_r[0]*ghat_rho_r[1]*dx1+0.5*urec_1_l[0]*ghat_rho_l[1]*dx1; 
  outrhouy[3] += (-0.5303300858899105*umax_r[0]*rhou1r[3]*dx1)-0.5303300858899105*umax_l[0]*rhou1l[3]*dx1-0.5303300858899105*umax_r[0]*rhou1c[3]*dx1-0.5303300858899105*umax_l[0]*rhou1c[3]*dx1+0.3061862178478971*umax_r[0]*rhou1r[2]*dx1-0.3061862178478971*umax_l[0]*rhou1l[2]*dx1-0.3061862178478971*umax_r[0]*rhou1c[2]*dx1+0.3061862178478971*umax_l[0]*rhou1c[2]*dx1-0.8660254037844386*ghat_rho_r[0]*urec_1_r[1]*dx1-0.8660254037844386*ghat_rho_l[0]*urec_1_l[1]*dx1-0.5303300858899105*rhou1r[1]*umax_r[1]*dx1-0.5303300858899105*rhou1c[1]*umax_r[1]*dx1+0.3061862178478971*rhou1r[0]*umax_r[1]*dx1-0.3061862178478971*rhou1c[0]*umax_r[1]*dx1-0.5303300858899105*rhou1l[1]*umax_l[1]*dx1-0.5303300858899105*rhou1c[1]*umax_l[1]*dx1-0.3061862178478971*rhou1l[0]*umax_l[1]*dx1+0.3061862178478971*rhou1c[0]*umax_l[1]*dx1-0.8660254037844386*urec_1_r[0]*ghat_rho_r[1]*dx1-0.8660254037844386*urec_1_l[0]*ghat_rho_l[1]*dx1; 

  outrhouz[0] += (-0.3061862178478971*umax_r[1]*rhou2r[3]*dx1)+0.3061862178478971*umax_l[1]*rhou2l[3]*dx1-0.3061862178478971*umax_r[1]*rhou2c[3]*dx1+0.3061862178478971*umax_l[1]*rhou2c[3]*dx1+0.1767766952966368*umax_r[1]*rhou2r[2]*dx1+0.1767766952966368*umax_l[1]*rhou2l[2]*dx1-0.1767766952966368*umax_r[1]*rhou2c[2]*dx1-0.1767766952966368*umax_l[1]*rhou2c[2]*dx1-0.5*ghat_rho_r[1]*urec_2_r[1]*dx1+0.5*ghat_rho_l[1]*urec_2_l[1]*dx1-0.3061862178478971*umax_r[0]*rhou2r[1]*dx1+0.3061862178478971*umax_l[0]*rhou2l[1]*dx1-0.3061862178478971*umax_r[0]*rhou2c[1]*dx1+0.3061862178478971*umax_l[0]*rhou2c[1]*dx1-0.5*ghat_rho_r[0]*urec_2_r[0]*dx1+0.5*ghat_rho_l[0]*urec_2_l[0]*dx1+0.1767766952966368*rhou2r[0]*umax_r[0]*dx1-0.1767766952966368*rhou2c[0]*umax_r[0]*dx1+0.1767766952966368*rhou2l[0]*umax_l[0]*dx1-0.1767766952966368*rhou2c[0]*umax_l[0]*dx1; 
  outrhouz[1] += (-0.5303300858899105*umax_r[1]*rhou2r[3]*dx1)-0.5303300858899105*umax_l[1]*rhou2l[3]*dx1-0.5303300858899105*umax_r[1]*rhou2c[3]*dx1-0.5303300858899105*umax_l[1]*rhou2c[3]*dx1+0.3061862178478971*umax_r[1]*rhou2r[2]*dx1-0.3061862178478971*umax_l[1]*rhou2l[2]*dx1-0.3061862178478971*umax_r[1]*rhou2c[2]*dx1+0.3061862178478971*umax_l[1]*rhou2c[2]*dx1-0.8660254037844386*ghat_rho_r[1]*urec_2_r[1]*dx1-0.8660254037844386*ghat_rho_l[1]*urec_2_l[1]*dx1-0.5303300858899105*umax_r[0]*rhou2r[1]*dx1-0.5303300858899105*umax_l[0]*rhou2l[1]*dx1-0.5303300858899105*umax_r[0]*rhou2c[1]*dx1-0.5303300858899105*umax_l[0]*rhou2c[1]*dx1-0.8660254037844386*ghat_rho_r[0]*urec_2_r[0]*dx1-0.8660254037844386*ghat_rho_l[0]*urec_2_l[0]*dx1+0.3061862178478971*rhou2r[0]*umax_r[0]*dx1-0.3061862178478971*rhou2c[0]*umax_r[0]*dx1-0.3061862178478971*rhou2l[0]*umax_l[0]*dx1+0.3061862178478971*rhou2c[0]*umax_l[0]*dx1; 
  outrhouz[2] += (-0.3061862178478971*umax_r[0]*rhou2r[3]*dx1)+0.3061862178478971*umax_l[0]*rhou2l[3]*dx1-0.3061862178478971*umax_r[0]*rhou2c[3]*dx1+0.3061862178478971*umax_l[0]*rhou2c[3]*dx1+0.1767766952966368*umax_r[0]*rhou2r[2]*dx1+0.1767766952966368*umax_l[0]*rhou2l[2]*dx1-0.1767766952966368*umax_r[0]*rhou2c[2]*dx1-0.1767766952966368*umax_l[0]*rhou2c[2]*dx1-0.5*ghat_rho_r[0]*urec_2_r[1]*dx1+0.5*ghat_rho_l[0]*urec_2_l[1]*dx1-0.3061862178478971*rhou2r[1]*umax_r[1]*dx1-0.3061862178478971*rhou2c[1]*umax_r[1]*dx1+0.1767766952966368*rhou2r[0]*umax_r[1]*dx1-0.1767766952966368*rhou2c[0]*umax_r[1]*dx1+0.3061862178478971*rhou2l[1]*umax_l[1]*dx1+0.3061862178478971*rhou2c[1]*umax_l[1]*dx1+0.1767766952966368*rhou2l[0]*umax_l[1]*dx1-0.1767766952966368*rhou2c[0]*umax_l[1]*dx1-0.5*urec_2_r[0]*ghat_rho_r[1]*dx1+0.5*urec_2_l[0]*ghat_rho_l[1]*dx1; 
  outrhouz[3] += (-0.5303300858899105*umax_r[0]*rhou2r[3]*dx1)-0.5303300858899105*umax_l[0]*rhou2l[3]*dx1-0.5303300858899105*umax_r[0]*rhou2c[3]*dx1-0.5303300858899105*umax_l[0]*rhou2c[3]*dx1+0.3061862178478971*umax_r[0]*rhou2r[2]*dx1-0.3061862178478971*umax_l[0]*rhou2l[2]*dx1-0.3061862178478971*umax_r[0]*rhou2c[2]*dx1+0.3061862178478971*umax_l[0]*rhou2c[2]*dx1-0.8660254037844386*ghat_rho_r[0]*urec_2_r[1]*dx1-0.8660254037844386*ghat_rho_l[0]*urec_2_l[1]*dx1-0.5303300858899105*rhou2r[1]*umax_r[1]*dx1-0.5303300858899105*rhou2c[1]*umax_r[1]*dx1+0.3061862178478971*rhou2r[0]*umax_r[1]*dx1-0.3061862178478971*rhou2c[0]*umax_r[1]*dx1-0.5303300858899105*rhou2l[1]*umax_l[1]*dx1-0.5303300858899105*rhou2c[1]*umax_l[1]*dx1-0.3061862178478971*rhou2l[0]*umax_l[1]*dx1+0.3061862178478971*rhou2c[0]*umax_l[1]*dx1-0.8660254037844386*urec_2_r[0]*ghat_rho_r[1]*dx1-0.8660254037844386*urec_2_l[0]*ghat_rho_l[1]*dx1; 

} 
#include <gkyl_basis_ser_2x_p1_surfx2_eval_quad.h> 
GKYL_CU_DH void isoeuler_surfy_2x_ser_p1(const double *w, const double *dxv, double vth, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out)  
{ 
  // w: Cell-center coordinates. dxv[NDIM]: Cell spacing. statevec(l/c/r): [rho, rho ux, rho uy, rho uz] in (left/center/right) cell, uvar(l/c/r): [ux, uv, uz]  in (left/center/right) cell
  // out: output in center cell 
  const double dx1 = 2.0/dxv[1]; 
  const double *rhol = &statevecl[0]; 
  const double *rhou0l = &statevecl[4]; 
  const double *rhou1l = &statevecl[8]; 
  const double *rhou2l = &statevecl[12]; 
  const double *rhoc = &statevecc[0]; 
  const double *rhou0c = &statevecc[4]; 
  const double *rhou1c = &statevecc[8]; 
  const double *rhou2c = &statevecc[12]; 
  const double *rhor = &statevecr[0]; 
  const double *rhou0r = &statevecr[4]; 
  const double *rhou1r = &statevecr[8]; 
  const double *rhou2r = &statevecr[12]; 
  const double *uvar0l = &uvarl[0]; 
  const double *uvar1l = &uvarl[4]; 
  const double *uvar2l = &uvarl[8]; 
  const double *uvar0c = &uvarc[0]; 
  const double *uvar1c = &uvarc[4]; 
  const double *uvar2c = &uvarc[8]; 
  const double *uvar0r = &uvarr[0]; 
  const double *uvar1r = &uvarr[4]; 
  const double *uvar2r = &uvarr[8]; 
  double *outrho = &out[0]; 
  double *outrhoux = &out[4]; 
  double *outrhouy = &out[8]; 
  double *outrhouz = &out[12]; 
  double incr[4]; 

  double vthsq = vth*vth; 
  double umax_l[2] = {0.0};
  double umax_r[2] = {0.0};
  double uquad_l[2] = {0.0};
  double uquad_r[2] = {0.0};
  double urec_0_l[2] = {0.0}; 
  double urec_0_r[2] = {0.0}; 
  double urec_1_l[2] = {0.0}; 
  double urec_1_r[2] = {0.0}; 
  double urec_2_l[2] = {0.0}; 
  double urec_2_r[2] = {0.0}; 
  double rhourec_0_l[2] = {0.0}; 
  double rhourec_0_r[2] = {0.0}; 
  double rhourec_1_l[2] = {0.0}; 
  double rhourec_1_r[2] = {0.0}; 
  double rhourec_2_l[2] = {0.0}; 
  double rhourec_2_r[2] = {0.0}; 
  double rhorec_l[2] = {0.0}; 
  double rhorec_r[2] = {0.0}; 
  double ghat_rho_l[2] = {0.0}; 
  double ghat_rho_r[2] = {0.0}; 
  double u_l_r = 0.0; 
  double u_c_l = 0.0; 
  double u_c_r = 0.0; 
  double u_r_l = 0.0; 

  u_l_r = ser_2x_p1_surfx2_eval_quad_node_0_r(uvar1l); 
  u_c_l = ser_2x_p1_surfx2_eval_quad_node_0_l(uvar1c); 
  u_c_r = ser_2x_p1_surfx2_eval_quad_node_0_r(uvar1c); 
  u_r_l = ser_2x_p1_surfx2_eval_quad_node_0_l(uvar1r); 
  uquad_l[0] = fmax(fabs(u_l_r), fabs(u_c_l)); 
  uquad_r[0] = fmax(fabs(u_c_r), fabs(u_r_l)); 
  u_l_r = ser_2x_p1_surfx2_eval_quad_node_1_r(uvar1l); 
  u_c_l = ser_2x_p1_surfx2_eval_quad_node_1_l(uvar1c); 
  u_c_r = ser_2x_p1_surfx2_eval_quad_node_1_r(uvar1c); 
  u_r_l = ser_2x_p1_surfx2_eval_quad_node_1_l(uvar1r); 
  uquad_l[1] = fmax(fabs(u_l_r), fabs(u_c_l)); 
  uquad_r[1] = fmax(fabs(u_c_r), fabs(u_r_l)); 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_2x_p1_upwind_quad_to_modal(uquad_l, umax_l); 
  ser_2x_p1_upwind_quad_to_modal(uquad_r, umax_r); 

  urec_0_l[0] = 0.408248290463863*uvar0l[2]-0.408248290463863*uvar0c[2]+0.3535533905932737*uvar0l[0]+0.3535533905932737*uvar0c[0]; 
  urec_0_l[1] = 0.408248290463863*uvar0l[3]-0.408248290463863*uvar0c[3]+0.3535533905932737*uvar0l[1]+0.3535533905932737*uvar0c[1]; 

  urec_1_l[0] = 0.408248290463863*uvar1l[2]-0.408248290463863*uvar1c[2]+0.3535533905932737*uvar1l[0]+0.3535533905932737*uvar1c[0]; 
  urec_1_l[1] = 0.408248290463863*uvar1l[3]-0.408248290463863*uvar1c[3]+0.3535533905932737*uvar1l[1]+0.3535533905932737*uvar1c[1]; 

  urec_2_l[0] = 0.408248290463863*uvar2l[2]-0.408248290463863*uvar2c[2]+0.3535533905932737*uvar2l[0]+0.3535533905932737*uvar2c[0]; 
  urec_2_l[1] = 0.408248290463863*uvar2l[3]-0.408248290463863*uvar2c[3]+0.3535533905932737*uvar2l[1]+0.3535533905932737*uvar2c[1]; 

  rhourec_0_l[0] = 0.408248290463863*rhou0l[2]-0.408248290463863*rhou0c[2]+0.3535533905932737*rhou0l[0]+0.3535533905932737*rhou0c[0]; 
  rhourec_0_l[1] = 0.408248290463863*rhou0l[3]-0.408248290463863*rhou0c[3]+0.3535533905932737*rhou0l[1]+0.3535533905932737*rhou0c[1]; 

  rhourec_1_l[0] = 0.408248290463863*rhou1l[2]-0.408248290463863*rhou1c[2]+0.3535533905932737*rhou1l[0]+0.3535533905932737*rhou1c[0]; 
  rhourec_1_l[1] = 0.408248290463863*rhou1l[3]-0.408248290463863*rhou1c[3]+0.3535533905932737*rhou1l[1]+0.3535533905932737*rhou1c[1]; 

  rhourec_2_l[0] = 0.408248290463863*rhou2l[2]-0.408248290463863*rhou2c[2]+0.3535533905932737*rhou2l[0]+0.3535533905932737*rhou2c[0]; 
  rhourec_2_l[1] = 0.408248290463863*rhou2l[3]-0.408248290463863*rhou2c[3]+0.3535533905932737*rhou2l[1]+0.3535533905932737*rhou2c[1]; 

  rhorec_l[0] = 0.408248290463863*rhol[2]-0.408248290463863*rhoc[2]+0.3535533905932737*rhol[0]+0.3535533905932737*rhoc[0]; 
  rhorec_l[1] = 0.408248290463863*rhol[3]-0.408248290463863*rhoc[3]+0.3535533905932737*rhol[1]+0.3535533905932737*rhoc[1]; 

  urec_0_r[0] = (-0.408248290463863*uvar0r[2])+0.408248290463863*uvar0c[2]+0.3535533905932737*uvar0r[0]+0.3535533905932737*uvar0c[0]; 
  urec_0_r[1] = (-0.408248290463863*uvar0r[3])+0.408248290463863*uvar0c[3]+0.3535533905932737*uvar0r[1]+0.3535533905932737*uvar0c[1]; 

  urec_1_r[0] = (-0.408248290463863*uvar1r[2])+0.408248290463863*uvar1c[2]+0.3535533905932737*uvar1r[0]+0.3535533905932737*uvar1c[0]; 
  urec_1_r[1] = (-0.408248290463863*uvar1r[3])+0.408248290463863*uvar1c[3]+0.3535533905932737*uvar1r[1]+0.3535533905932737*uvar1c[1]; 

  urec_2_r[0] = (-0.408248290463863*uvar2l[2])+0.408248290463863*uvar2c[2]+0.3535533905932737*uvar2l[0]+0.3535533905932737*uvar2c[0]; 
  urec_2_r[1] = (-0.408248290463863*uvar2l[3])+0.408248290463863*uvar2c[3]+0.3535533905932737*uvar2l[1]+0.3535533905932737*uvar2c[1]; 

  rhourec_0_r[0] = (-0.408248290463863*rhou0r[2])+0.408248290463863*rhou0c[2]+0.3535533905932737*rhou0r[0]+0.3535533905932737*rhou0c[0]; 
  rhourec_0_r[1] = (-0.408248290463863*rhou0r[3])+0.408248290463863*rhou0c[3]+0.3535533905932737*rhou0r[1]+0.3535533905932737*rhou0c[1]; 

  rhourec_1_r[0] = (-0.408248290463863*rhou1r[2])+0.408248290463863*rhou1c[2]+0.3535533905932737*rhou1r[0]+0.3535533905932737*rhou1c[0]; 
  rhourec_1_r[1] = (-0.408248290463863*rhou1r[3])+0.408248290463863*rhou1c[3]+0.3535533905932737*rhou1r[1]+0.3535533905932737*rhou1c[1]; 

  rhourec_2_r[0] = (-0.408248290463863*rhou2r[2])+0.408248290463863*rhou2c[2]+0.3535533905932737*rhou2r[0]+0.3535533905932737*rhou2c[0]; 
  rhourec_2_r[1] = (-0.408248290463863*rhou2r[3])+0.408248290463863*rhou2c[3]+0.3535533905932737*rhou2r[1]+0.3535533905932737*rhou2c[1]; 

  rhorec_r[0] = (-0.408248290463863*rhor[2])+0.408248290463863*rhoc[2]+0.3535533905932737*rhor[0]+0.3535533905932737*rhoc[0]; 
  rhorec_r[1] = (-0.408248290463863*rhor[3])+0.408248290463863*rhoc[3]+0.3535533905932737*rhor[1]+0.3535533905932737*rhoc[1]; 

  ghat_rho_l[0] = 0.4330127018922193*umax_l[1]*rhol[3]+0.4330127018922193*umax_l[1]*rhoc[3]+0.4330127018922193*umax_l[0]*rhol[2]+0.4330127018922193*umax_l[0]*rhoc[2]+0.25*rhol[1]*umax_l[1]-0.25*rhoc[1]*umax_l[1]+0.25*rhol[0]*umax_l[0]-0.25*rhoc[0]*umax_l[0]+rhourec_1_l[0]; 
  ghat_rho_l[1] = 0.4330127018922193*umax_l[0]*rhol[3]+0.4330127018922193*umax_l[0]*rhoc[3]+0.4330127018922193*umax_l[1]*rhol[2]+0.4330127018922193*umax_l[1]*rhoc[2]+0.25*rhol[0]*umax_l[1]-0.25*rhoc[0]*umax_l[1]+rhourec_1_l[1]+0.25*umax_l[0]*rhol[1]-0.25*umax_l[0]*rhoc[1]; 

  ghat_rho_r[0] = 0.4330127018922193*umax_r[1]*rhor[3]+0.4330127018922193*umax_r[1]*rhoc[3]+0.4330127018922193*umax_r[0]*rhor[2]+0.4330127018922193*umax_r[0]*rhoc[2]-0.25*rhor[1]*umax_r[1]+0.25*rhoc[1]*umax_r[1]-0.25*rhor[0]*umax_r[0]+0.25*rhoc[0]*umax_r[0]+rhourec_1_r[0]; 
  ghat_rho_r[1] = 0.4330127018922193*umax_r[0]*rhor[3]+0.4330127018922193*umax_r[0]*rhoc[3]+0.4330127018922193*umax_r[1]*rhor[2]+0.4330127018922193*umax_r[1]*rhoc[2]-0.25*rhor[0]*umax_r[1]+0.25*rhoc[0]*umax_r[1]+rhourec_1_r[1]-0.25*umax_r[0]*rhor[1]+0.25*umax_r[0]*rhoc[1]; 

  outrho[0] += 0.7071067811865475*ghat_rho_l[0]*dx1-0.7071067811865475*ghat_rho_r[0]*dx1; 
  outrho[1] += 0.7071067811865475*ghat_rho_l[1]*dx1-0.7071067811865475*ghat_rho_r[1]*dx1; 
  outrho[2] += (-1.224744871391589*ghat_rho_r[0]*dx1)-1.224744871391589*ghat_rho_l[0]*dx1; 
  outrho[3] += (-1.224744871391589*ghat_rho_r[1]*dx1)-1.224744871391589*ghat_rho_l[1]*dx1; 

  outrhoux[0] += (-0.3061862178478971*umax_r[1]*rhou0r[3]*dx1)+0.3061862178478971*umax_l[1]*rhou0l[3]*dx1-0.3061862178478971*umax_r[1]*rhou0c[3]*dx1+0.3061862178478971*umax_l[1]*rhou0c[3]*dx1-0.3061862178478971*umax_r[0]*rhou0r[2]*dx1+0.3061862178478971*umax_l[0]*rhou0l[2]*dx1-0.3061862178478971*umax_r[0]*rhou0c[2]*dx1+0.3061862178478971*umax_l[0]*rhou0c[2]*dx1-0.5*ghat_rho_r[1]*urec_0_r[1]*dx1+0.5*ghat_rho_l[1]*urec_0_l[1]*dx1+0.1767766952966368*rhou0r[1]*umax_r[1]*dx1-0.1767766952966368*rhou0c[1]*umax_r[1]*dx1+0.1767766952966368*rhou0l[1]*umax_l[1]*dx1-0.1767766952966368*rhou0c[1]*umax_l[1]*dx1-0.5*ghat_rho_r[0]*urec_0_r[0]*dx1+0.5*ghat_rho_l[0]*urec_0_l[0]*dx1+0.1767766952966368*rhou0r[0]*umax_r[0]*dx1-0.1767766952966368*rhou0c[0]*umax_r[0]*dx1+0.1767766952966368*rhou0l[0]*umax_l[0]*dx1-0.1767766952966368*rhou0c[0]*umax_l[0]*dx1; 
  outrhoux[1] += (-0.3061862178478971*umax_r[0]*rhou0r[3]*dx1)+0.3061862178478971*umax_l[0]*rhou0l[3]*dx1-0.3061862178478971*umax_r[0]*rhou0c[3]*dx1+0.3061862178478971*umax_l[0]*rhou0c[3]*dx1-0.3061862178478971*umax_r[1]*rhou0r[2]*dx1+0.3061862178478971*umax_l[1]*rhou0l[2]*dx1-0.3061862178478971*umax_r[1]*rhou0c[2]*dx1+0.3061862178478971*umax_l[1]*rhou0c[2]*dx1-0.5*ghat_rho_r[0]*urec_0_r[1]*dx1+0.5*ghat_rho_l[0]*urec_0_l[1]*dx1+0.1767766952966368*rhou0r[0]*umax_r[1]*dx1-0.1767766952966368*rhou0c[0]*umax_r[1]*dx1+0.1767766952966368*rhou0l[0]*umax_l[1]*dx1-0.1767766952966368*rhou0c[0]*umax_l[1]*dx1+0.1767766952966368*umax_r[0]*rhou0r[1]*dx1+0.1767766952966368*umax_l[0]*rhou0l[1]*dx1-0.1767766952966368*umax_r[0]*rhou0c[1]*dx1-0.1767766952966368*umax_l[0]*rhou0c[1]*dx1-0.5*urec_0_r[0]*ghat_rho_r[1]*dx1+0.5*urec_0_l[0]*ghat_rho_l[1]*dx1; 
  outrhoux[2] += (-0.5303300858899105*umax_r[1]*rhou0r[3]*dx1)-0.5303300858899105*umax_l[1]*rhou0l[3]*dx1-0.5303300858899105*umax_r[1]*rhou0c[3]*dx1-0.5303300858899105*umax_l[1]*rhou0c[3]*dx1-0.5303300858899105*umax_r[0]*rhou0r[2]*dx1-0.5303300858899105*umax_l[0]*rhou0l[2]*dx1-0.5303300858899105*umax_r[0]*rhou0c[2]*dx1-0.5303300858899105*umax_l[0]*rhou0c[2]*dx1-0.8660254037844386*ghat_rho_r[1]*urec_0_r[1]*dx1-0.8660254037844386*ghat_rho_l[1]*urec_0_l[1]*dx1+0.3061862178478971*rhou0r[1]*umax_r[1]*dx1-0.3061862178478971*rhou0c[1]*umax_r[1]*dx1-0.3061862178478971*rhou0l[1]*umax_l[1]*dx1+0.3061862178478971*rhou0c[1]*umax_l[1]*dx1-0.8660254037844386*ghat_rho_r[0]*urec_0_r[0]*dx1-0.8660254037844386*ghat_rho_l[0]*urec_0_l[0]*dx1+0.3061862178478971*rhou0r[0]*umax_r[0]*dx1-0.3061862178478971*rhou0c[0]*umax_r[0]*dx1-0.3061862178478971*rhou0l[0]*umax_l[0]*dx1+0.3061862178478971*rhou0c[0]*umax_l[0]*dx1; 
  outrhoux[3] += (-0.5303300858899105*umax_r[0]*rhou0r[3]*dx1)-0.5303300858899105*umax_l[0]*rhou0l[3]*dx1-0.5303300858899105*umax_r[0]*rhou0c[3]*dx1-0.5303300858899105*umax_l[0]*rhou0c[3]*dx1-0.5303300858899105*umax_r[1]*rhou0r[2]*dx1-0.5303300858899105*umax_l[1]*rhou0l[2]*dx1-0.5303300858899105*umax_r[1]*rhou0c[2]*dx1-0.5303300858899105*umax_l[1]*rhou0c[2]*dx1-0.8660254037844386*ghat_rho_r[0]*urec_0_r[1]*dx1-0.8660254037844386*ghat_rho_l[0]*urec_0_l[1]*dx1+0.3061862178478971*rhou0r[0]*umax_r[1]*dx1-0.3061862178478971*rhou0c[0]*umax_r[1]*dx1-0.3061862178478971*rhou0l[0]*umax_l[1]*dx1+0.3061862178478971*rhou0c[0]*umax_l[1]*dx1+0.3061862178478971*umax_r[0]*rhou0r[1]*dx1-0.3061862178478971*umax_l[0]*rhou0l[1]*dx1-0.3061862178478971*umax_r[0]*rhou0c[1]*dx1+0.3061862178478971*umax_l[0]*rhou0c[1]*dx1-0.8660254037844386*urec_0_r[0]*ghat_rho_r[1]*dx1-0.8660254037844386*urec_0_l[0]*ghat_rho_l[1]*dx1; 

  outrhouy[0] += (-0.7071067811865475*rhorec_r[0]*dx1*vthsq)+0.7071067811865475*rhorec_l[0]*dx1*vthsq-0.3061862178478971*umax_r[1]*rhou1r[3]*dx1+0.3061862178478971*umax_l[1]*rhou1l[3]*dx1-0.3061862178478971*umax_r[1]*rhou1c[3]*dx1+0.3061862178478971*umax_l[1]*rhou1c[3]*dx1-0.3061862178478971*umax_r[0]*rhou1r[2]*dx1+0.3061862178478971*umax_l[0]*rhou1l[2]*dx1-0.3061862178478971*umax_r[0]*rhou1c[2]*dx1+0.3061862178478971*umax_l[0]*rhou1c[2]*dx1-0.5*ghat_rho_r[1]*urec_1_r[1]*dx1+0.5*ghat_rho_l[1]*urec_1_l[1]*dx1+0.1767766952966368*rhou1r[1]*umax_r[1]*dx1-0.1767766952966368*rhou1c[1]*umax_r[1]*dx1+0.1767766952966368*rhou1l[1]*umax_l[1]*dx1-0.1767766952966368*rhou1c[1]*umax_l[1]*dx1-0.5*ghat_rho_r[0]*urec_1_r[0]*dx1+0.5*ghat_rho_l[0]*urec_1_l[0]*dx1+0.1767766952966368*rhou1r[0]*umax_r[0]*dx1-0.1767766952966368*rhou1c[0]*umax_r[0]*dx1+0.1767766952966368*rhou1l[0]*umax_l[0]*dx1-0.1767766952966368*rhou1c[0]*umax_l[0]*dx1; 
  outrhouy[1] += (-0.7071067811865475*rhorec_r[1]*dx1*vthsq)+0.7071067811865475*rhorec_l[1]*dx1*vthsq-0.3061862178478971*umax_r[0]*rhou1r[3]*dx1+0.3061862178478971*umax_l[0]*rhou1l[3]*dx1-0.3061862178478971*umax_r[0]*rhou1c[3]*dx1+0.3061862178478971*umax_l[0]*rhou1c[3]*dx1-0.3061862178478971*umax_r[1]*rhou1r[2]*dx1+0.3061862178478971*umax_l[1]*rhou1l[2]*dx1-0.3061862178478971*umax_r[1]*rhou1c[2]*dx1+0.3061862178478971*umax_l[1]*rhou1c[2]*dx1-0.5*ghat_rho_r[0]*urec_1_r[1]*dx1+0.5*ghat_rho_l[0]*urec_1_l[1]*dx1+0.1767766952966368*rhou1r[0]*umax_r[1]*dx1-0.1767766952966368*rhou1c[0]*umax_r[1]*dx1+0.1767766952966368*rhou1l[0]*umax_l[1]*dx1-0.1767766952966368*rhou1c[0]*umax_l[1]*dx1+0.1767766952966368*umax_r[0]*rhou1r[1]*dx1+0.1767766952966368*umax_l[0]*rhou1l[1]*dx1-0.1767766952966368*umax_r[0]*rhou1c[1]*dx1-0.1767766952966368*umax_l[0]*rhou1c[1]*dx1-0.5*urec_1_r[0]*ghat_rho_r[1]*dx1+0.5*urec_1_l[0]*ghat_rho_l[1]*dx1; 
  outrhouy[2] += (-1.224744871391589*rhorec_r[0]*dx1*vthsq)-1.224744871391589*rhorec_l[0]*dx1*vthsq-0.5303300858899105*umax_r[1]*rhou1r[3]*dx1-0.5303300858899105*umax_l[1]*rhou1l[3]*dx1-0.5303300858899105*umax_r[1]*rhou1c[3]*dx1-0.5303300858899105*umax_l[1]*rhou1c[3]*dx1-0.5303300858899105*umax_r[0]*rhou1r[2]*dx1-0.5303300858899105*umax_l[0]*rhou1l[2]*dx1-0.5303300858899105*umax_r[0]*rhou1c[2]*dx1-0.5303300858899105*umax_l[0]*rhou1c[2]*dx1-0.8660254037844386*ghat_rho_r[1]*urec_1_r[1]*dx1-0.8660254037844386*ghat_rho_l[1]*urec_1_l[1]*dx1+0.3061862178478971*rhou1r[1]*umax_r[1]*dx1-0.3061862178478971*rhou1c[1]*umax_r[1]*dx1-0.3061862178478971*rhou1l[1]*umax_l[1]*dx1+0.3061862178478971*rhou1c[1]*umax_l[1]*dx1-0.8660254037844386*ghat_rho_r[0]*urec_1_r[0]*dx1-0.8660254037844386*ghat_rho_l[0]*urec_1_l[0]*dx1+0.3061862178478971*rhou1r[0]*umax_r[0]*dx1-0.3061862178478971*rhou1c[0]*umax_r[0]*dx1-0.3061862178478971*rhou1l[0]*umax_l[0]*dx1+0.3061862178478971*rhou1c[0]*umax_l[0]*dx1; 
  outrhouy[3] += (-1.224744871391589*rhorec_r[1]*dx1*vthsq)-1.224744871391589*rhorec_l[1]*dx1*vthsq-0.5303300858899105*umax_r[0]*rhou1r[3]*dx1-0.5303300858899105*umax_l[0]*rhou1l[3]*dx1-0.5303300858899105*umax_r[0]*rhou1c[3]*dx1-0.5303300858899105*umax_l[0]*rhou1c[3]*dx1-0.5303300858899105*umax_r[1]*rhou1r[2]*dx1-0.5303300858899105*umax_l[1]*rhou1l[2]*dx1-0.5303300858899105*umax_r[1]*rhou1c[2]*dx1-0.5303300858899105*umax_l[1]*rhou1c[2]*dx1-0.8660254037844386*ghat_rho_r[0]*urec_1_r[1]*dx1-0.8660254037844386*ghat_rho_l[0]*urec_1_l[1]*dx1+0.3061862178478971*rhou1r[0]*umax_r[1]*dx1-0.3061862178478971*rhou1c[0]*umax_r[1]*dx1-0.3061862178478971*rhou1l[0]*umax_l[1]*dx1+0.3061862178478971*rhou1c[0]*umax_l[1]*dx1+0.3061862178478971*umax_r[0]*rhou1r[1]*dx1-0.3061862178478971*umax_l[0]*rhou1l[1]*dx1-0.3061862178478971*umax_r[0]*rhou1c[1]*dx1+0.3061862178478971*umax_l[0]*rhou1c[1]*dx1-0.8660254037844386*urec_1_r[0]*ghat_rho_r[1]*dx1-0.8660254037844386*urec_1_l[0]*ghat_rho_l[1]*dx1; 

  outrhouz[0] += (-0.3061862178478971*umax_r[1]*rhou2r[3]*dx1)+0.3061862178478971*umax_l[1]*rhou2l[3]*dx1-0.3061862178478971*umax_r[1]*rhou2c[3]*dx1+0.3061862178478971*umax_l[1]*rhou2c[3]*dx1-0.3061862178478971*umax_r[0]*rhou2r[2]*dx1+0.3061862178478971*umax_l[0]*rhou2l[2]*dx1-0.3061862178478971*umax_r[0]*rhou2c[2]*dx1+0.3061862178478971*umax_l[0]*rhou2c[2]*dx1-0.5*ghat_rho_r[1]*urec_2_r[1]*dx1+0.5*ghat_rho_l[1]*urec_2_l[1]*dx1+0.1767766952966368*rhou2r[1]*umax_r[1]*dx1-0.1767766952966368*rhou2c[1]*umax_r[1]*dx1+0.1767766952966368*rhou2l[1]*umax_l[1]*dx1-0.1767766952966368*rhou2c[1]*umax_l[1]*dx1-0.5*ghat_rho_r[0]*urec_2_r[0]*dx1+0.5*ghat_rho_l[0]*urec_2_l[0]*dx1+0.1767766952966368*rhou2r[0]*umax_r[0]*dx1-0.1767766952966368*rhou2c[0]*umax_r[0]*dx1+0.1767766952966368*rhou2l[0]*umax_l[0]*dx1-0.1767766952966368*rhou2c[0]*umax_l[0]*dx1; 
  outrhouz[1] += (-0.3061862178478971*umax_r[0]*rhou2r[3]*dx1)+0.3061862178478971*umax_l[0]*rhou2l[3]*dx1-0.3061862178478971*umax_r[0]*rhou2c[3]*dx1+0.3061862178478971*umax_l[0]*rhou2c[3]*dx1-0.3061862178478971*umax_r[1]*rhou2r[2]*dx1+0.3061862178478971*umax_l[1]*rhou2l[2]*dx1-0.3061862178478971*umax_r[1]*rhou2c[2]*dx1+0.3061862178478971*umax_l[1]*rhou2c[2]*dx1-0.5*ghat_rho_r[0]*urec_2_r[1]*dx1+0.5*ghat_rho_l[0]*urec_2_l[1]*dx1+0.1767766952966368*rhou2r[0]*umax_r[1]*dx1-0.1767766952966368*rhou2c[0]*umax_r[1]*dx1+0.1767766952966368*rhou2l[0]*umax_l[1]*dx1-0.1767766952966368*rhou2c[0]*umax_l[1]*dx1+0.1767766952966368*umax_r[0]*rhou2r[1]*dx1+0.1767766952966368*umax_l[0]*rhou2l[1]*dx1-0.1767766952966368*umax_r[0]*rhou2c[1]*dx1-0.1767766952966368*umax_l[0]*rhou2c[1]*dx1-0.5*urec_2_r[0]*ghat_rho_r[1]*dx1+0.5*urec_2_l[0]*ghat_rho_l[1]*dx1; 
  outrhouz[2] += (-0.5303300858899105*umax_r[1]*rhou2r[3]*dx1)-0.5303300858899105*umax_l[1]*rhou2l[3]*dx1-0.5303300858899105*umax_r[1]*rhou2c[3]*dx1-0.5303300858899105*umax_l[1]*rhou2c[3]*dx1-0.5303300858899105*umax_r[0]*rhou2r[2]*dx1-0.5303300858899105*umax_l[0]*rhou2l[2]*dx1-0.5303300858899105*umax_r[0]*rhou2c[2]*dx1-0.5303300858899105*umax_l[0]*rhou2c[2]*dx1-0.8660254037844386*ghat_rho_r[1]*urec_2_r[1]*dx1-0.8660254037844386*ghat_rho_l[1]*urec_2_l[1]*dx1+0.3061862178478971*rhou2r[1]*umax_r[1]*dx1-0.3061862178478971*rhou2c[1]*umax_r[1]*dx1-0.3061862178478971*rhou2l[1]*umax_l[1]*dx1+0.3061862178478971*rhou2c[1]*umax_l[1]*dx1-0.8660254037844386*ghat_rho_r[0]*urec_2_r[0]*dx1-0.8660254037844386*ghat_rho_l[0]*urec_2_l[0]*dx1+0.3061862178478971*rhou2r[0]*umax_r[0]*dx1-0.3061862178478971*rhou2c[0]*umax_r[0]*dx1-0.3061862178478971*rhou2l[0]*umax_l[0]*dx1+0.3061862178478971*rhou2c[0]*umax_l[0]*dx1; 
  outrhouz[3] += (-0.5303300858899105*umax_r[0]*rhou2r[3]*dx1)-0.5303300858899105*umax_l[0]*rhou2l[3]*dx1-0.5303300858899105*umax_r[0]*rhou2c[3]*dx1-0.5303300858899105*umax_l[0]*rhou2c[3]*dx1-0.5303300858899105*umax_r[1]*rhou2r[2]*dx1-0.5303300858899105*umax_l[1]*rhou2l[2]*dx1-0.5303300858899105*umax_r[1]*rhou2c[2]*dx1-0.5303300858899105*umax_l[1]*rhou2c[2]*dx1-0.8660254037844386*ghat_rho_r[0]*urec_2_r[1]*dx1-0.8660254037844386*ghat_rho_l[0]*urec_2_l[1]*dx1+0.3061862178478971*rhou2r[0]*umax_r[1]*dx1-0.3061862178478971*rhou2c[0]*umax_r[1]*dx1-0.3061862178478971*rhou2l[0]*umax_l[1]*dx1+0.3061862178478971*rhou2c[0]*umax_l[1]*dx1+0.3061862178478971*umax_r[0]*rhou2r[1]*dx1-0.3061862178478971*umax_l[0]*rhou2l[1]*dx1-0.3061862178478971*umax_r[0]*rhou2c[1]*dx1+0.3061862178478971*umax_l[0]*rhou2c[1]*dx1-0.8660254037844386*urec_2_r[0]*ghat_rho_r[1]*dx1-0.8660254037844386*urec_2_l[0]*ghat_rho_l[1]*dx1; 

} 
