#include <gkyl_pkpm_vars_kernels.h> 
GKYL_CU_DH void pkpm_vars_grad_u_z_3x_ser_p1(const double *dxv, const double *u_l, const double *u_c, const double *u_r, double* GKYL_RESTRICT grad_u) 
{ 
  // dxv[NDIM]: Cell spacing. 
  // u_l/u_c/u_r: Flow velocity [ux, uy, uz] in the left/current/right cell in direction z. 
  // grad_u:    Output weak gradient of u, block of direction z: [d_z ux, d_z uy, d_z uz]. 
  const double dx1 = 2.0/dxv[2]; 
  const double *ux_l = &u_l[0]; 
  const double *ux_c = &u_c[0]; 
  const double *ux_r = &u_r[0]; 
  double *grad_ux = &grad_u[48]; 
  const double *uy_l = &u_l[8]; 
  const double *uy_c = &u_c[8]; 
  const double *uy_r = &u_r[8]; 
  double *grad_uy = &grad_u[56]; 
  const double *uz_l = &u_l[16]; 
  const double *uz_c = &u_c[16]; 
  const double *uz_r = &u_r[16]; 
  double *grad_uz = &grad_u[64]; 
 
  grad_ux[0] = -(0.4330127018922193*ux_r[3]*dx1)-0.4330127018922193*ux_l[3]*dx1+0.8660254037844386*ux_c[3]*dx1+0.25*ux_r[0]*dx1-0.25*ux_l[0]*dx1; 
  grad_ux[1] = -(0.4330127018922193*ux_r[5]*dx1)-0.4330127018922193*ux_l[5]*dx1+0.8660254037844386*ux_c[5]*dx1+0.25*ux_r[1]*dx1-0.25*ux_l[1]*dx1; 
  grad_ux[2] = -(0.4330127018922193*ux_r[6]*dx1)-0.4330127018922193*ux_l[6]*dx1+0.8660254037844386*ux_c[6]*dx1+0.25*ux_r[2]*dx1-0.25*ux_l[2]*dx1; 
  grad_ux[3] = -(0.75*ux_r[3]*dx1)+0.75*ux_l[3]*dx1+0.4330127018922193*ux_r[0]*dx1+0.4330127018922193*ux_l[0]*dx1-0.8660254037844385*ux_c[0]*dx1; 
  grad_ux[4] = -(0.4330127018922193*ux_r[7]*dx1)-0.4330127018922193*ux_l[7]*dx1+0.8660254037844386*ux_c[7]*dx1+0.25*ux_r[4]*dx1-0.25*ux_l[4]*dx1; 
  grad_ux[5] = -(0.75*ux_r[5]*dx1)+0.75*ux_l[5]*dx1+0.4330127018922193*ux_r[1]*dx1+0.4330127018922193*ux_l[1]*dx1-0.8660254037844385*ux_c[1]*dx1; 
  grad_ux[6] = -(0.75*ux_r[6]*dx1)+0.75*ux_l[6]*dx1+0.4330127018922193*ux_r[2]*dx1+0.4330127018922193*ux_l[2]*dx1-0.8660254037844385*ux_c[2]*dx1; 
  grad_ux[7] = -(0.75*ux_r[7]*dx1)+0.75*ux_l[7]*dx1+0.4330127018922193*ux_r[4]*dx1+0.4330127018922193*ux_l[4]*dx1-0.8660254037844385*ux_c[4]*dx1; 
  grad_uy[0] = -(0.4330127018922193*uy_r[3]*dx1)-0.4330127018922193*uy_l[3]*dx1+0.8660254037844386*uy_c[3]*dx1+0.25*uy_r[0]*dx1-0.25*uy_l[0]*dx1; 
  grad_uy[1] = -(0.4330127018922193*uy_r[5]*dx1)-0.4330127018922193*uy_l[5]*dx1+0.8660254037844386*uy_c[5]*dx1+0.25*uy_r[1]*dx1-0.25*uy_l[1]*dx1; 
  grad_uy[2] = -(0.4330127018922193*uy_r[6]*dx1)-0.4330127018922193*uy_l[6]*dx1+0.8660254037844386*uy_c[6]*dx1+0.25*uy_r[2]*dx1-0.25*uy_l[2]*dx1; 
  grad_uy[3] = -(0.75*uy_r[3]*dx1)+0.75*uy_l[3]*dx1+0.4330127018922193*uy_r[0]*dx1+0.4330127018922193*uy_l[0]*dx1-0.8660254037844385*uy_c[0]*dx1; 
  grad_uy[4] = -(0.4330127018922193*uy_r[7]*dx1)-0.4330127018922193*uy_l[7]*dx1+0.8660254037844386*uy_c[7]*dx1+0.25*uy_r[4]*dx1-0.25*uy_l[4]*dx1; 
  grad_uy[5] = -(0.75*uy_r[5]*dx1)+0.75*uy_l[5]*dx1+0.4330127018922193*uy_r[1]*dx1+0.4330127018922193*uy_l[1]*dx1-0.8660254037844385*uy_c[1]*dx1; 
  grad_uy[6] = -(0.75*uy_r[6]*dx1)+0.75*uy_l[6]*dx1+0.4330127018922193*uy_r[2]*dx1+0.4330127018922193*uy_l[2]*dx1-0.8660254037844385*uy_c[2]*dx1; 
  grad_uy[7] = -(0.75*uy_r[7]*dx1)+0.75*uy_l[7]*dx1+0.4330127018922193*uy_r[4]*dx1+0.4330127018922193*uy_l[4]*dx1-0.8660254037844385*uy_c[4]*dx1; 
  grad_uz[0] = -(0.4330127018922193*uz_r[3]*dx1)-0.4330127018922193*uz_l[3]*dx1+0.8660254037844386*uz_c[3]*dx1+0.25*uz_r[0]*dx1-0.25*uz_l[0]*dx1; 
  grad_uz[1] = -(0.4330127018922193*uz_r[5]*dx1)-0.4330127018922193*uz_l[5]*dx1+0.8660254037844386*uz_c[5]*dx1+0.25*uz_r[1]*dx1-0.25*uz_l[1]*dx1; 
  grad_uz[2] = -(0.4330127018922193*uz_r[6]*dx1)-0.4330127018922193*uz_l[6]*dx1+0.8660254037844386*uz_c[6]*dx1+0.25*uz_r[2]*dx1-0.25*uz_l[2]*dx1; 
  grad_uz[3] = -(0.75*uz_r[3]*dx1)+0.75*uz_l[3]*dx1+0.4330127018922193*uz_r[0]*dx1+0.4330127018922193*uz_l[0]*dx1-0.8660254037844385*uz_c[0]*dx1; 
  grad_uz[4] = -(0.4330127018922193*uz_r[7]*dx1)-0.4330127018922193*uz_l[7]*dx1+0.8660254037844386*uz_c[7]*dx1+0.25*uz_r[4]*dx1-0.25*uz_l[4]*dx1; 
  grad_uz[5] = -(0.75*uz_r[5]*dx1)+0.75*uz_l[5]*dx1+0.4330127018922193*uz_r[1]*dx1+0.4330127018922193*uz_l[1]*dx1-0.8660254037844385*uz_c[1]*dx1; 
  grad_uz[6] = -(0.75*uz_r[6]*dx1)+0.75*uz_l[6]*dx1+0.4330127018922193*uz_r[2]*dx1+0.4330127018922193*uz_l[2]*dx1-0.8660254037844385*uz_c[2]*dx1; 
  grad_uz[7] = -(0.75*uz_r[7]*dx1)+0.75*uz_l[7]*dx1+0.4330127018922193*uz_r[4]*dx1+0.4330127018922193*uz_l[4]*dx1-0.8660254037844385*uz_c[4]*dx1; 
} 
