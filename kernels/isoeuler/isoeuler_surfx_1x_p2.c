#include <gkyl_isoeuler_kernels.h> 
#include <gkyl_basis_ser_1x_p2_surfx1_eval_quad.h> 
GKYL_CU_DH void isoeuler_surfx_1x_ser_p2(const double *w, const double *dxv, double vth, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out)  
{ 
  // w: Cell-center coordinates. dxv[NDIM]: Cell spacing. statevec(l/c/r): [rho, rho ux, rho uy, rho uz] in (left/center/right) cell, uvar(l/c/r): [ux, uv, uz]  in (left/center/right) cell
  // out: output in center cell 
  const double dx1 = 2.0/dxv[0]; 
  const double *rhol = &statevecl[0]; 
  const double *rhou0l = &statevecl[3]; 
  const double *rhou1l = &statevecl[6]; 
  const double *rhou2l = &statevecl[9]; 
  const double *rhoc = &statevecc[0]; 
  const double *rhou0c = &statevecc[3]; 
  const double *rhou1c = &statevecc[6]; 
  const double *rhou2c = &statevecc[9]; 
  const double *rhor = &statevecr[0]; 
  const double *rhou0r = &statevecr[3]; 
  const double *rhou1r = &statevecr[6]; 
  const double *rhou2r = &statevecr[9]; 
  const double *uvar0l = &uvarl[0]; 
  const double *uvar1l = &uvarl[3]; 
  const double *uvar2l = &uvarl[6]; 
  const double *uvar0c = &uvarc[0]; 
  const double *uvar1c = &uvarc[3]; 
  const double *uvar2c = &uvarc[6]; 
  const double *uvar0r = &uvarr[0]; 
  const double *uvar1r = &uvarr[3]; 
  const double *uvar2r = &uvarr[6]; 
  double *outrho = &out[0]; 
  double *outrhoux = &out[3]; 
  double *outrhouy = &out[6]; 
  double *outrhouz = &out[9]; 
  double incr[3]; 

  double vthsq = vth*vth; 
  double u_l_r = ser_1x_p2_surfx1_eval_quad_node_0_r(uvar0l); 
  double u_c_l = ser_1x_p2_surfx1_eval_quad_node_0_l(uvar0c); 
  double u_c_r = ser_1x_p2_surfx1_eval_quad_node_0_r(uvar0c); 
  double u_r_l = ser_1x_p2_surfx1_eval_quad_node_0_l(uvar0r); 
  double u_max_l = fmax(fabs(u_l_r), fabs(u_c_l)); 
  double u_max_r = fmax(fabs(u_c_r), fabs(u_r_l)); 

  double urec_0_l = 0.3458741190809163*uvar0l[2]+0.3458741190809163*uvar0c[2]+0.4975526040028326*uvar0l[1]-0.4975526040028326*uvar0c[1]+0.3535533905932737*uvar0l[0]+0.3535533905932737*uvar0c[0]; 
  double urec_0_r = 0.3458741190809163*uvar0r[2]+0.3458741190809163*uvar0c[2]-0.4975526040028326*uvar0r[1]+0.4975526040028326*uvar0c[1]+0.3535533905932737*uvar0r[0]+0.3535533905932737*uvar0c[0]; 
  double urec_1_l = 0.3458741190809163*uvar1l[2]+0.3458741190809163*uvar1c[2]+0.4975526040028326*uvar1l[1]-0.4975526040028326*uvar1c[1]+0.3535533905932737*uvar1l[0]+0.3535533905932737*uvar1c[0]; 
  double urec_1_r = 0.3458741190809163*uvar1r[2]+0.3458741190809163*uvar1c[2]-0.4975526040028326*uvar1r[1]+0.4975526040028326*uvar1c[1]+0.3535533905932737*uvar1r[0]+0.3535533905932737*uvar1c[0]; 
  double urec_2_l = 0.3458741190809163*uvar2l[2]+0.3458741190809163*uvar2c[2]+0.4975526040028326*uvar2l[1]-0.4975526040028326*uvar2c[1]+0.3535533905932737*uvar2l[0]+0.3535533905932737*uvar2c[0]; 
  double urec_2_r = 0.3458741190809163*uvar2l[2]+0.3458741190809163*uvar2c[2]-0.4975526040028326*uvar2l[1]+0.4975526040028326*uvar2c[1]+0.3535533905932737*uvar2l[0]+0.3535533905932737*uvar2c[0]; 
  double rhourec_0_l = 0.3458741190809163*rhou0l[2]+0.3458741190809163*rhou0c[2]+0.4975526040028326*rhou0l[1]-0.4975526040028326*rhou0c[1]+0.3535533905932737*rhou0l[0]+0.3535533905932737*rhou0c[0]; 
  double rhourec_0_r = 0.3458741190809163*rhou0r[2]+0.3458741190809163*rhou0c[2]-0.4975526040028326*rhou0r[1]+0.4975526040028326*rhou0c[1]+0.3535533905932737*rhou0r[0]+0.3535533905932737*rhou0c[0]; 
  double rhourec_1_l = 0.3458741190809163*rhou1l[2]+0.3458741190809163*rhou1c[2]+0.4975526040028326*rhou1l[1]-0.4975526040028326*rhou1c[1]+0.3535533905932737*rhou1l[0]+0.3535533905932737*rhou1c[0]; 
  double rhourec_1_r = 0.3458741190809163*rhou1r[2]+0.3458741190809163*rhou1c[2]-0.4975526040028326*rhou1r[1]+0.4975526040028326*rhou1c[1]+0.3535533905932737*rhou1r[0]+0.3535533905932737*rhou1c[0]; 
  double rhourec_2_l = 0.3458741190809163*rhou2l[2]+0.3458741190809163*rhou2c[2]+0.4975526040028326*rhou2l[1]-0.4975526040028326*rhou2c[1]+0.3535533905932737*rhou2l[0]+0.3535533905932737*rhou2c[0]; 
  double rhourec_2_r = 0.3458741190809163*rhou2r[2]+0.3458741190809163*rhou2c[2]-0.4975526040028326*rhou2r[1]+0.4975526040028326*rhou2c[1]+0.3535533905932737*rhou2r[0]+0.3535533905932737*rhou2c[0]; 
  double rhorec_l = 0.3458741190809163*rhol[2]+0.3458741190809163*rhoc[2]+0.4975526040028326*rhol[1]-0.4975526040028326*rhoc[1]+0.3535533905932737*rhol[0]+0.3535533905932737*rhoc[0]; 
  double rhorec_r = 0.3458741190809163*rhor[2]+0.3458741190809163*rhoc[2]-0.4975526040028326*rhor[1]+0.4975526040028326*rhoc[1]+0.3535533905932737*rhor[0]+0.3535533905932737*rhoc[0]; 

  double ghat_rho_l = 0.7905694150420948*rhol[2]*u_max_l-0.7905694150420948*rhoc[2]*u_max_l+0.6123724356957945*rhol[1]*u_max_l+0.6123724356957945*rhoc[1]*u_max_l+0.3535533905932737*rhol[0]*u_max_l-0.3535533905932737*rhoc[0]*u_max_l+rhourec_0_l; 
  double ghat_rho_r = (-0.7905694150420948*rhor[2]*u_max_r)+0.7905694150420948*rhoc[2]*u_max_r+0.6123724356957945*rhor[1]*u_max_r+0.6123724356957945*rhoc[1]*u_max_r-0.3535533905932737*rhor[0]*u_max_r+0.3535533905932737*rhoc[0]*u_max_r+rhourec_0_r; 
  outrho[0] += 0.7071067811865475*dx1*ghat_rho_l-0.7071067811865475*dx1*ghat_rho_r; 
  outrho[1] += (-1.224744871391589*dx1*ghat_rho_r)-1.224744871391589*dx1*ghat_rho_l; 
  outrho[2] += 1.58113883008419*dx1*ghat_rho_l-1.58113883008419*dx1*ghat_rho_r; 

  outrhoux[0] += (-0.7071067811865475*dx1*rhorec_r*vthsq)+0.7071067811865475*dx1*rhorec_l*vthsq-0.7071067811865475*dx1*ghat_rho_r*urec_0_r+0.7071067811865475*dx1*ghat_rho_l*urec_0_l+0.5590169943749475*rhou0r[2]*dx1*u_max_r-0.5590169943749475*rhou0c[2]*dx1*u_max_r-0.4330127018922193*rhou0r[1]*dx1*u_max_r-0.4330127018922193*rhou0c[1]*dx1*u_max_r+0.25*rhou0r[0]*dx1*u_max_r-0.25*rhou0c[0]*dx1*u_max_r+0.5590169943749475*rhou0l[2]*dx1*u_max_l-0.5590169943749475*rhou0c[2]*dx1*u_max_l+0.4330127018922193*rhou0l[1]*dx1*u_max_l+0.4330127018922193*rhou0c[1]*dx1*u_max_l+0.25*rhou0l[0]*dx1*u_max_l-0.25*rhou0c[0]*dx1*u_max_l; 
  outrhoux[1] += (-1.224744871391589*dx1*rhorec_r*vthsq)-1.224744871391589*dx1*rhorec_l*vthsq-1.224744871391589*dx1*ghat_rho_r*urec_0_r-1.224744871391589*dx1*ghat_rho_l*urec_0_l+0.9682458365518543*rhou0r[2]*dx1*u_max_r-0.9682458365518543*rhou0c[2]*dx1*u_max_r-0.75*rhou0r[1]*dx1*u_max_r-0.75*rhou0c[1]*dx1*u_max_r+0.4330127018922193*rhou0r[0]*dx1*u_max_r-0.4330127018922193*rhou0c[0]*dx1*u_max_r-0.9682458365518543*rhou0l[2]*dx1*u_max_l+0.9682458365518543*rhou0c[2]*dx1*u_max_l-0.75*rhou0l[1]*dx1*u_max_l-0.75*rhou0c[1]*dx1*u_max_l-0.4330127018922193*rhou0l[0]*dx1*u_max_l+0.4330127018922193*rhou0c[0]*dx1*u_max_l; 
  outrhoux[2] += (-1.58113883008419*dx1*rhorec_r*vthsq)+1.58113883008419*dx1*rhorec_l*vthsq-1.58113883008419*dx1*ghat_rho_r*urec_0_r+1.58113883008419*dx1*ghat_rho_l*urec_0_l+1.25*rhou0r[2]*dx1*u_max_r-1.25*rhou0c[2]*dx1*u_max_r-0.9682458365518543*rhou0r[1]*dx1*u_max_r-0.9682458365518543*rhou0c[1]*dx1*u_max_r+0.5590169943749475*rhou0r[0]*dx1*u_max_r-0.5590169943749475*rhou0c[0]*dx1*u_max_r+1.25*rhou0l[2]*dx1*u_max_l-1.25*rhou0c[2]*dx1*u_max_l+0.9682458365518543*rhou0l[1]*dx1*u_max_l+0.9682458365518543*rhou0c[1]*dx1*u_max_l+0.5590169943749475*rhou0l[0]*dx1*u_max_l-0.5590169943749475*rhou0c[0]*dx1*u_max_l; 

  outrhouy[0] += (-0.7071067811865475*dx1*ghat_rho_r*urec_1_r)+0.7071067811865475*dx1*ghat_rho_l*urec_1_l+0.5590169943749475*rhou1r[2]*dx1*u_max_r-0.5590169943749475*rhou1c[2]*dx1*u_max_r-0.4330127018922193*rhou1r[1]*dx1*u_max_r-0.4330127018922193*rhou1c[1]*dx1*u_max_r+0.25*rhou1r[0]*dx1*u_max_r-0.25*rhou1c[0]*dx1*u_max_r+0.5590169943749475*rhou1l[2]*dx1*u_max_l-0.5590169943749475*rhou1c[2]*dx1*u_max_l+0.4330127018922193*rhou1l[1]*dx1*u_max_l+0.4330127018922193*rhou1c[1]*dx1*u_max_l+0.25*rhou1l[0]*dx1*u_max_l-0.25*rhou1c[0]*dx1*u_max_l; 
  outrhouy[1] += (-1.224744871391589*dx1*ghat_rho_r*urec_1_r)-1.224744871391589*dx1*ghat_rho_l*urec_1_l+0.9682458365518543*rhou1r[2]*dx1*u_max_r-0.9682458365518543*rhou1c[2]*dx1*u_max_r-0.75*rhou1r[1]*dx1*u_max_r-0.75*rhou1c[1]*dx1*u_max_r+0.4330127018922193*rhou1r[0]*dx1*u_max_r-0.4330127018922193*rhou1c[0]*dx1*u_max_r-0.9682458365518543*rhou1l[2]*dx1*u_max_l+0.9682458365518543*rhou1c[2]*dx1*u_max_l-0.75*rhou1l[1]*dx1*u_max_l-0.75*rhou1c[1]*dx1*u_max_l-0.4330127018922193*rhou1l[0]*dx1*u_max_l+0.4330127018922193*rhou1c[0]*dx1*u_max_l; 
  outrhouy[2] += (-1.58113883008419*dx1*ghat_rho_r*urec_1_r)+1.58113883008419*dx1*ghat_rho_l*urec_1_l+1.25*rhou1r[2]*dx1*u_max_r-1.25*rhou1c[2]*dx1*u_max_r-0.9682458365518543*rhou1r[1]*dx1*u_max_r-0.9682458365518543*rhou1c[1]*dx1*u_max_r+0.5590169943749475*rhou1r[0]*dx1*u_max_r-0.5590169943749475*rhou1c[0]*dx1*u_max_r+1.25*rhou1l[2]*dx1*u_max_l-1.25*rhou1c[2]*dx1*u_max_l+0.9682458365518543*rhou1l[1]*dx1*u_max_l+0.9682458365518543*rhou1c[1]*dx1*u_max_l+0.5590169943749475*rhou1l[0]*dx1*u_max_l-0.5590169943749475*rhou1c[0]*dx1*u_max_l; 

  outrhouz[0] += (-0.7071067811865475*dx1*ghat_rho_r*urec_2_r)+0.7071067811865475*dx1*ghat_rho_l*urec_2_l+0.5590169943749475*rhou2r[2]*dx1*u_max_r-0.5590169943749475*rhou2c[2]*dx1*u_max_r-0.4330127018922193*rhou2r[1]*dx1*u_max_r-0.4330127018922193*rhou2c[1]*dx1*u_max_r+0.25*rhou2r[0]*dx1*u_max_r-0.25*rhou2c[0]*dx1*u_max_r+0.5590169943749475*rhou2l[2]*dx1*u_max_l-0.5590169943749475*rhou2c[2]*dx1*u_max_l+0.4330127018922193*rhou2l[1]*dx1*u_max_l+0.4330127018922193*rhou2c[1]*dx1*u_max_l+0.25*rhou2l[0]*dx1*u_max_l-0.25*rhou2c[0]*dx1*u_max_l; 
  outrhouz[1] += (-1.224744871391589*dx1*ghat_rho_r*urec_2_r)-1.224744871391589*dx1*ghat_rho_l*urec_2_l+0.9682458365518543*rhou2r[2]*dx1*u_max_r-0.9682458365518543*rhou2c[2]*dx1*u_max_r-0.75*rhou2r[1]*dx1*u_max_r-0.75*rhou2c[1]*dx1*u_max_r+0.4330127018922193*rhou2r[0]*dx1*u_max_r-0.4330127018922193*rhou2c[0]*dx1*u_max_r-0.9682458365518543*rhou2l[2]*dx1*u_max_l+0.9682458365518543*rhou2c[2]*dx1*u_max_l-0.75*rhou2l[1]*dx1*u_max_l-0.75*rhou2c[1]*dx1*u_max_l-0.4330127018922193*rhou2l[0]*dx1*u_max_l+0.4330127018922193*rhou2c[0]*dx1*u_max_l; 
  outrhouz[2] += (-1.58113883008419*dx1*ghat_rho_r*urec_2_r)+1.58113883008419*dx1*ghat_rho_l*urec_2_l+1.25*rhou2r[2]*dx1*u_max_r-1.25*rhou2c[2]*dx1*u_max_r-0.9682458365518543*rhou2r[1]*dx1*u_max_r-0.9682458365518543*rhou2c[1]*dx1*u_max_r+0.5590169943749475*rhou2r[0]*dx1*u_max_r-0.5590169943749475*rhou2c[0]*dx1*u_max_r+1.25*rhou2l[2]*dx1*u_max_l-1.25*rhou2c[2]*dx1*u_max_l+0.9682458365518543*rhou2l[1]*dx1*u_max_l+0.9682458365518543*rhou2c[1]*dx1*u_max_l+0.5590169943749475*rhou2l[0]*dx1*u_max_l-0.5590169943749475*rhou2c[0]*dx1*u_max_l; 

} 
