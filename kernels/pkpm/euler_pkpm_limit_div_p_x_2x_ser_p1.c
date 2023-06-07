#include <gkyl_euler_pkpm_kernels.h> 
GKYL_CU_DH void euler_pkpm_limit_div_p_x_2x_ser_p1(const double *dxv, 
  const double *u_il, const double *u_ic, const double *u_ir, 
  const double *p_ijl, const double *p_ijc, const double *p_ijr, 
  double* div_p_cell_avg) 
{ 
  // dxv[NDIM]:             Cell spacing.
  // u_il/c/r:              Input bulk velocity (ux,uy,uz) in left/center/right cells.
  // p_ijl/c/r:             Input pressure tensor in left/center/right cells.
  // div_p_cell_avg:        Output d/dx_i p_ij cell average for use in limiter.

  const double dx1 = 2.0/dxv[0]; 
  const double *ux_l = &u_il[0]; 
  const double *uy_l = &u_il[4]; 
  const double *uz_l = &u_il[8]; 

  const double *ux_c = &u_ic[0]; 
  const double *uy_c = &u_ic[4]; 
  const double *uz_c = &u_ic[8]; 

  const double *ux_r = &u_ir[0]; 
  const double *uy_r = &u_ir[4]; 
  const double *uz_r = &u_ir[8]; 

  const double *Pxx_l = &p_ijl[0]; 
  const double *Pxy_l = &p_ijl[4]; 
  const double *Pxz_l = &p_ijl[8]; 
  const double *Pyy_l = &p_ijl[12]; 
  const double *Pyz_l = &p_ijl[16]; 
  const double *Pzz_l = &p_ijl[20]; 

  const double *Pxx_c = &p_ijc[0]; 
  const double *Pxy_c = &p_ijc[4]; 
  const double *Pxz_c = &p_ijc[8]; 
  const double *Pyy_c = &p_ijc[12]; 
  const double *Pyz_c = &p_ijc[16]; 
  const double *Pzz_c = &p_ijc[20]; 

  const double *Pxx_r = &p_ijr[0]; 
  const double *Pxy_r = &p_ijr[4]; 
  const double *Pxz_r = &p_ijr[8]; 
  const double *Pyy_r = &p_ijr[12]; 
  const double *Pyz_r = &p_ijr[16]; 
  const double *Pzz_r = &p_ijr[20]; 

  double *div_p_x = &div_p_cell_avg[0]; 
  double *div_p_y = &div_p_cell_avg[1]; 
  double *div_p_z = &div_p_cell_avg[2]; 

  double *grad_u_x = &div_p_cell_avg[9]; 
  double *grad_u_y = &div_p_cell_avg[10]; 
  double *grad_u_z = &div_p_cell_avg[11]; 

  div_p_x[0] = ((-0.2886751345948129*Pxx_r[1])-0.2886751345948129*Pxx_l[1]+0.5773502691896258*Pxx_c[1]+0.25*Pxx_r[0]-0.25*Pxx_l[0])*dx1; 
  div_p_y[0] = ((-0.2886751345948129*Pxy_r[1])-0.2886751345948129*Pxy_l[1]+0.5773502691896258*Pxy_c[1]+0.25*Pxy_r[0]-0.25*Pxy_l[0])*dx1; 
  div_p_z[0] = ((-0.2886751345948129*Pxz_r[1])-0.2886751345948129*Pxz_l[1]+0.5773502691896258*Pxz_c[1]+0.25*Pxz_r[0]-0.25*Pxz_l[0])*dx1; 
  grad_u_x[0] = ((-0.2886751345948129*ux_r[1])-0.2886751345948129*ux_l[1]+0.5773502691896258*ux_c[1]+0.25*ux_r[0]-0.25*ux_l[0])*dx1; 
  grad_u_y[0] = ((-0.2886751345948129*uy_r[1])-0.2886751345948129*uy_l[1]+0.5773502691896258*uy_c[1]+0.25*uy_r[0]-0.25*uy_l[0])*dx1; 
  grad_u_z[0] = ((-0.2886751345948129*uz_r[1])-0.2886751345948129*uz_l[1]+0.5773502691896258*uz_c[1]+0.25*uz_r[0]-0.25*uz_l[0])*dx1; 
} 
