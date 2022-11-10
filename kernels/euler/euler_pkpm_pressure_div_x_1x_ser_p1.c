#include <gkyl_euler_kernels.h> 
GKYL_CU_DH void euler_pkpm_pressure_div_x_1x_ser_p1(const double *dxv, const double *p_ijl, const double *p_ijc, const double *p_ijr, double* GKYL_RESTRICT div_p) 
{ 
  // dxv[NDIM]: Cell spacing.
  // p_ijl/p_ijc/p_ijr: Pressure tensor in left/center/right cells.
  // div_p: Volume expansion of div(p).

  const double dx1 = 2.0/dxv[0]; 

  const double *Pxx_l = &p_ijl[0]; 
  const double *Pxy_l = &p_ijl[2]; 
  const double *Pxz_l = &p_ijl[4]; 
  const double *Pyy_l = &p_ijl[6]; 
  const double *Pyz_l = &p_ijl[8]; 
  const double *Pzz_l = &p_ijl[10]; 

  const double *Pxx_c = &p_ijc[0]; 
  const double *Pxy_c = &p_ijc[2]; 
  const double *Pxz_c = &p_ijc[4]; 
  const double *Pyy_c = &p_ijc[6]; 
  const double *Pyz_c = &p_ijc[8]; 
  const double *Pzz_c = &p_ijc[10]; 

  const double *Pxx_r = &p_ijr[0]; 
  const double *Pxy_r = &p_ijr[2]; 
  const double *Pxz_r = &p_ijr[4]; 
  const double *Pyy_r = &p_ijr[6]; 
  const double *Pyz_r = &p_ijr[8]; 
  const double *Pzz_r = &p_ijr[10]; 

  double *div_p_x = &div_p[0]; 
  double *div_p_y = &div_p[2]; 
  double *div_p_z = &div_p[4]; 
  div_p_x[0] = (-0.2886751345948129*Pxx_r[1]*dx1)-0.2886751345948129*Pxx_l[1]*dx1+0.5773502691896258*Pxx_c[1]*dx1+0.25*Pxx_r[0]*dx1-0.25*Pxx_l[0]*dx1; 
  div_p_x[1] = (-0.5*Pxx_r[1]*dx1)+0.5*Pxx_l[1]*dx1+0.4330127018922193*Pxx_r[0]*dx1+0.4330127018922193*Pxx_l[0]*dx1-0.8660254037844386*Pxx_c[0]*dx1; 

  div_p_y[0] = (-0.2886751345948129*Pxy_r[1]*dx1)-0.2886751345948129*Pxy_l[1]*dx1+0.5773502691896258*Pxy_c[1]*dx1+0.25*Pxy_r[0]*dx1-0.25*Pxy_l[0]*dx1; 
  div_p_y[1] = (-0.5*Pxy_r[1]*dx1)+0.5*Pxy_l[1]*dx1+0.4330127018922193*Pxy_r[0]*dx1+0.4330127018922193*Pxy_l[0]*dx1-0.8660254037844386*Pxy_c[0]*dx1; 

  div_p_z[0] = (-0.2886751345948129*Pxz_r[1]*dx1)-0.2886751345948129*Pxz_l[1]*dx1+0.5773502691896258*Pxz_c[1]*dx1+0.25*Pxz_r[0]*dx1-0.25*Pxz_l[0]*dx1; 
  div_p_z[1] = (-0.5*Pxz_r[1]*dx1)+0.5*Pxz_l[1]*dx1+0.4330127018922193*Pxz_r[0]*dx1+0.4330127018922193*Pxz_l[0]*dx1-0.8660254037844386*Pxz_c[0]*dx1; 

} 
