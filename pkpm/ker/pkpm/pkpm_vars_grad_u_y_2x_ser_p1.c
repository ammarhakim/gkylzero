#include <gkyl_pkpm_vars_kernels.h> 
GKYL_CU_DH void pkpm_vars_grad_u_y_2x_ser_p1(const double *dxv, const double *jacob_pos, const double *u_l, const double *u_c, const double *u_r, double* GKYL_RESTRICT grad_u) 
{ 
  // dxv[NDIM]: Cell spacing (computational). 
  // jacob_pos: Position-map Jacobian per configuration direction (cdim blocks of 2 equal entries). 
  // u_l/u_c/u_r: Flow velocity [ux, uy, uz] in the left/current/right cell in direction y. 
  // grad_u:    Output weak gradient of u, block of direction y: [d_y ux, d_y uy, d_y uz]. 
  const double dx1 = 2.0/(dxv[1]*jacob_pos[2]); 
  const double *ux_l = &u_l[0]; 
  const double *ux_c = &u_c[0]; 
  const double *ux_r = &u_r[0]; 
  double *grad_ux = &grad_u[12]; 
  const double *uy_l = &u_l[4]; 
  const double *uy_c = &u_c[4]; 
  const double *uy_r = &u_r[4]; 
  double *grad_uy = &grad_u[16]; 
  const double *uz_l = &u_l[8]; 
  const double *uz_c = &u_c[8]; 
  const double *uz_r = &u_r[8]; 
  double *grad_uz = &grad_u[20]; 
 
  grad_ux[0] = -(0.4330127018922193*ux_r[2]*dx1)-0.4330127018922193*ux_l[2]*dx1+0.8660254037844386*ux_c[2]*dx1+0.25*ux_r[0]*dx1-0.25*ux_l[0]*dx1; 
  grad_ux[1] = -(0.4330127018922193*ux_r[3]*dx1)-0.4330127018922193*ux_l[3]*dx1+0.8660254037844386*ux_c[3]*dx1+0.25*ux_r[1]*dx1-0.25*ux_l[1]*dx1; 
  grad_ux[2] = -(0.75*ux_r[2]*dx1)+0.75*ux_l[2]*dx1+0.4330127018922193*ux_r[0]*dx1+0.4330127018922193*ux_l[0]*dx1-0.8660254037844385*ux_c[0]*dx1; 
  grad_ux[3] = -(0.75*ux_r[3]*dx1)+0.75*ux_l[3]*dx1+0.4330127018922193*ux_r[1]*dx1+0.4330127018922193*ux_l[1]*dx1-0.8660254037844385*ux_c[1]*dx1; 
  grad_uy[0] = -(0.4330127018922193*uy_r[2]*dx1)-0.4330127018922193*uy_l[2]*dx1+0.8660254037844386*uy_c[2]*dx1+0.25*uy_r[0]*dx1-0.25*uy_l[0]*dx1; 
  grad_uy[1] = -(0.4330127018922193*uy_r[3]*dx1)-0.4330127018922193*uy_l[3]*dx1+0.8660254037844386*uy_c[3]*dx1+0.25*uy_r[1]*dx1-0.25*uy_l[1]*dx1; 
  grad_uy[2] = -(0.75*uy_r[2]*dx1)+0.75*uy_l[2]*dx1+0.4330127018922193*uy_r[0]*dx1+0.4330127018922193*uy_l[0]*dx1-0.8660254037844385*uy_c[0]*dx1; 
  grad_uy[3] = -(0.75*uy_r[3]*dx1)+0.75*uy_l[3]*dx1+0.4330127018922193*uy_r[1]*dx1+0.4330127018922193*uy_l[1]*dx1-0.8660254037844385*uy_c[1]*dx1; 
  grad_uz[0] = -(0.4330127018922193*uz_r[2]*dx1)-0.4330127018922193*uz_l[2]*dx1+0.8660254037844386*uz_c[2]*dx1+0.25*uz_r[0]*dx1-0.25*uz_l[0]*dx1; 
  grad_uz[1] = -(0.4330127018922193*uz_r[3]*dx1)-0.4330127018922193*uz_l[3]*dx1+0.8660254037844386*uz_c[3]*dx1+0.25*uz_r[1]*dx1-0.25*uz_l[1]*dx1; 
  grad_uz[2] = -(0.75*uz_r[2]*dx1)+0.75*uz_l[2]*dx1+0.4330127018922193*uz_r[0]*dx1+0.4330127018922193*uz_l[0]*dx1-0.8660254037844385*uz_c[0]*dx1; 
  grad_uz[3] = -(0.75*uz_r[3]*dx1)+0.75*uz_l[3]*dx1+0.4330127018922193*uz_r[1]*dx1+0.4330127018922193*uz_l[1]*dx1-0.8660254037844385*uz_c[1]*dx1; 
} 
