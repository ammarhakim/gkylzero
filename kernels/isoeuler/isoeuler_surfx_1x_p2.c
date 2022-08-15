#include <gkyl_isoeuler_kernels.h> 
#include <gkyl_basis_ser_1x_p2_surfx1_eval_quad.h> 
GKYL_CU_DH void isoeuler_surfx_1x_ser_p2(const double *w, const double *dxv, double vth, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out)  
{ 
  // w: Cell-center coordinates. dxv[NDIM]: Cell spacing. statevec(l/c/r): [rho, rho ux, rho uy, rho uz] in (left/center/right) cell, uvar(l/c/r): [ux, uv, uz]  in (left/center/right) cell
  // out: output in left/right cells 
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
  //FluxRho; 
  double ghat_rho_l = 0.7905694150420948*rhol[2]*u_max_l-0.7905694150420948*rhoc[2]*u_max_l+0.6123724356957945*rhol[1]*u_max_l+0.6123724356957945*rhoc[1]*u_max_l+0.3535533905932737*rhol[0]*u_max_l-0.3535533905932737*rhoc[0]*u_max_l+0.7905694150420948*rhou0l[2]+0.7905694150420948*rhou0c[2]+0.6123724356957945*rhou0l[1]-0.6123724356957945*rhou0c[1]+0.3535533905932737*rhou0l[0]+0.3535533905932737*rhou0c[0]; 
  double ghat_rho_r = (-0.7905694150420948*rhor[2]*u_max_r)+0.7905694150420948*rhoc[2]*u_max_r+0.6123724356957945*rhor[1]*u_max_r+0.6123724356957945*rhoc[1]*u_max_r-0.3535533905932737*rhor[0]*u_max_r+0.3535533905932737*rhoc[0]*u_max_r+0.7905694150420948*rhou0r[2]+0.7905694150420948*rhou0c[2]-0.6123724356957945*rhou0r[1]+0.6123724356957945*rhou0c[1]+0.3535533905932737*rhou0r[0]+0.3535533905932737*rhou0c[0]; 
  outrho[0] += 0.7071067811865475*dx1*ghat_rho_l-0.7071067811865475*dx1*ghat_rho_r; 
  outrho[1] += (-1.224744871391589*dx1*ghat_rho_r)-1.224744871391589*dx1*ghat_rho_l; 
  outrho[2] += 1.58113883008419*dx1*ghat_rho_l-1.58113883008419*dx1*ghat_rho_r; 

  outrhoux[0] += (-0.5590169943749475*rhor[2]*dx1*vthsq)+0.5590169943749475*rhol[2]*dx1*vthsq+0.4330127018922193*rhor[1]*dx1*vthsq+0.4330127018922193*rhol[1]*dx1*vthsq-0.8660254037844386*rhoc[1]*dx1*vthsq-0.25*rhor[0]*dx1*vthsq+0.25*rhol[0]*dx1*vthsq-0.5590169943749475*uvar0r[2]*dx1*ghat_rho_l+0.5590169943749475*uvar0l[2]*dx1*ghat_rho_l+0.4330127018922193*uvar0r[1]*dx1*ghat_rho_l+0.4330127018922193*uvar0l[1]*dx1*ghat_rho_l-0.8660254037844386*uvar0c[1]*dx1*ghat_rho_l-0.25*uvar0r[0]*dx1*ghat_rho_l+0.25*uvar0l[0]*dx1*ghat_rho_l; 
  outrhoux[1] += (-0.9682458365518543*rhor[2]*dx1*vthsq)-0.9682458365518543*rhol[2]*dx1*vthsq-1.936491673103709*rhoc[2]*dx1*vthsq+0.75*rhor[1]*dx1*vthsq-0.75*rhol[1]*dx1*vthsq-0.4330127018922193*rhor[0]*dx1*vthsq-0.4330127018922193*rhol[0]*dx1*vthsq-0.8660254037844386*rhoc[0]*dx1*vthsq-0.9682458365518543*uvar0r[2]*dx1*ghat_rho_l-0.9682458365518543*uvar0l[2]*dx1*ghat_rho_l-1.936491673103709*uvar0c[2]*dx1*ghat_rho_l+0.75*uvar0r[1]*dx1*ghat_rho_l-0.75*uvar0l[1]*dx1*ghat_rho_l-0.4330127018922193*uvar0r[0]*dx1*ghat_rho_l-0.4330127018922193*uvar0l[0]*dx1*ghat_rho_l-0.8660254037844386*uvar0c[0]*dx1*ghat_rho_l; 
  outrhoux[2] += (-1.25*rhor[2]*dx1*vthsq)+1.25*rhol[2]*dx1*vthsq+0.9682458365518543*rhor[1]*dx1*vthsq+0.9682458365518543*rhol[1]*dx1*vthsq-1.936491673103709*rhoc[1]*dx1*vthsq-0.5590169943749475*rhor[0]*dx1*vthsq+0.5590169943749475*rhol[0]*dx1*vthsq-1.25*uvar0r[2]*dx1*ghat_rho_l+1.25*uvar0l[2]*dx1*ghat_rho_l+0.9682458365518543*uvar0r[1]*dx1*ghat_rho_l+0.9682458365518543*uvar0l[1]*dx1*ghat_rho_l-1.936491673103709*uvar0c[1]*dx1*ghat_rho_l-0.5590169943749475*uvar0r[0]*dx1*ghat_rho_l+0.5590169943749475*uvar0l[0]*dx1*ghat_rho_l; 

  outrhouy[0] += (-0.5590169943749475*uvar1r[2]*dx1*ghat_rho_l)+0.5590169943749475*uvar1l[2]*dx1*ghat_rho_l+0.4330127018922193*uvar1r[1]*dx1*ghat_rho_l+0.4330127018922193*uvar1l[1]*dx1*ghat_rho_l-0.8660254037844386*uvar1c[1]*dx1*ghat_rho_l-0.25*uvar1r[0]*dx1*ghat_rho_l+0.25*uvar1l[0]*dx1*ghat_rho_l; 
  outrhouy[1] += (-0.9682458365518543*uvar1r[2]*dx1*ghat_rho_l)-0.9682458365518543*uvar1l[2]*dx1*ghat_rho_l-1.936491673103709*uvar1c[2]*dx1*ghat_rho_l+0.75*uvar1r[1]*dx1*ghat_rho_l-0.75*uvar1l[1]*dx1*ghat_rho_l-0.4330127018922193*uvar1r[0]*dx1*ghat_rho_l-0.4330127018922193*uvar1l[0]*dx1*ghat_rho_l-0.8660254037844386*uvar1c[0]*dx1*ghat_rho_l; 
  outrhouy[2] += (-1.25*uvar1r[2]*dx1*ghat_rho_l)+1.25*uvar1l[2]*dx1*ghat_rho_l+0.9682458365518543*uvar1r[1]*dx1*ghat_rho_l+0.9682458365518543*uvar1l[1]*dx1*ghat_rho_l-1.936491673103709*uvar1c[1]*dx1*ghat_rho_l-0.5590169943749475*uvar1r[0]*dx1*ghat_rho_l+0.5590169943749475*uvar1l[0]*dx1*ghat_rho_l; 

  outrhouz[0] += (-0.5590169943749475*uvar2r[2]*dx1*ghat_rho_l)+0.5590169943749475*uvar2l[2]*dx1*ghat_rho_l+0.4330127018922193*uvar2r[1]*dx1*ghat_rho_l+0.4330127018922193*uvar2l[1]*dx1*ghat_rho_l-0.8660254037844386*uvar2c[1]*dx1*ghat_rho_l-0.25*uvar2r[0]*dx1*ghat_rho_l+0.25*uvar2l[0]*dx1*ghat_rho_l; 
  outrhouz[1] += (-0.9682458365518543*uvar2r[2]*dx1*ghat_rho_l)-0.9682458365518543*uvar2l[2]*dx1*ghat_rho_l-1.936491673103709*uvar2c[2]*dx1*ghat_rho_l+0.75*uvar2r[1]*dx1*ghat_rho_l-0.75*uvar2l[1]*dx1*ghat_rho_l-0.4330127018922193*uvar2r[0]*dx1*ghat_rho_l-0.4330127018922193*uvar2l[0]*dx1*ghat_rho_l-0.8660254037844386*uvar2c[0]*dx1*ghat_rho_l; 
  outrhouz[2] += (-1.25*uvar2r[2]*dx1*ghat_rho_l)+1.25*uvar2l[2]*dx1*ghat_rho_l+0.9682458365518543*uvar2r[1]*dx1*ghat_rho_l+0.9682458365518543*uvar2l[1]*dx1*ghat_rho_l-1.936491673103709*uvar2c[1]*dx1*ghat_rho_l-0.5590169943749475*uvar2r[0]*dx1*ghat_rho_l+0.5590169943749475*uvar2l[0]*dx1*ghat_rho_l; 

} 
