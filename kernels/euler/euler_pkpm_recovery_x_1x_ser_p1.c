#include <gkyl_euler_kernels.h> 
GKYL_CU_DH void euler_pkpm_recovery_x_1x_ser_p1(const double *dxv, const double *bvarl, const double *bvarc, const double *bvarr, const double *u_il, const double *u_ic, const double *u_ir, const double *p_ijl, const double *p_ijc, const double *p_ijr, double* div_b, double* bb_grad_u, double* div_p) 
{ 
  // dxv[NDIM]: Cell spacing.
  // Al/Ac/Ar:  Inpute vector in left/center/right cells.
  // out:       Increment to volume expansion of div(A) in one direction.

  const double dx1 = 2.0/dxv[0]; 

  const double *b_l = &bvarl[0]; 
  const double *b_c = &bvarc[0]; 
  const double *b_r = &bvarr[0]; 

  const double *ux_l = &u_il[0]; 
  const double *uy_l = &u_il[2]; 
  const double *uz_l = &u_il[4]; 

  const double *ux_c = &u_ic[0]; 
  const double *uy_c = &u_ic[2]; 
  const double *uz_c = &u_ic[4]; 

  const double *ux_r = &u_ir[0]; 
  const double *uy_r = &u_ir[2]; 
  const double *uz_r = &u_ir[4]; 

  const double *bxbx = &bvarc[6]; 
  const double *bxby = &bvarc[8]; 
  const double *bxbz = &bvarc[10]; 
  const double *byby = &bvarc[12]; 
  const double *bybz = &bvarc[14]; 
  const double *bzbz = &bvarc[16]; 

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
  double grad_u_x[2] = {0.0}; 
  double grad_u_y[2] = {0.0}; 
  double grad_u_z[2] = {0.0}; 
  grad_u_x[0] = (-0.2886751345948129*ux_r[1])-0.2886751345948129*ux_l[1]+0.5773502691896258*ux_c[1]+0.25*ux_r[0]-0.25*ux_l[0]; 
  grad_u_x[1] = (-0.5*ux_r[1])+0.5*ux_l[1]+0.4330127018922193*ux_r[0]+0.4330127018922193*ux_l[0]-0.8660254037844386*ux_c[0]; 

  grad_u_y[0] = (-0.2886751345948129*uy_r[1])-0.2886751345948129*uy_l[1]+0.5773502691896258*uy_c[1]+0.25*uy_r[0]-0.25*uy_l[0]; 
  grad_u_y[1] = (-0.5*uy_r[1])+0.5*uy_l[1]+0.4330127018922193*uy_r[0]+0.4330127018922193*uy_l[0]-0.8660254037844386*uy_c[0]; 

  grad_u_z[0] = (-0.2886751345948129*uz_r[1])-0.2886751345948129*uz_l[1]+0.5773502691896258*uz_c[1]+0.25*uz_r[0]-0.25*uz_l[0]; 
  grad_u_z[1] = (-0.5*uz_r[1])+0.5*uz_l[1]+0.4330127018922193*uz_r[0]+0.4330127018922193*uz_l[0]-0.8660254037844386*uz_c[0]; 

  div_b[0] += ((-0.2886751345948129*(b_r[1]+b_l[1]))+0.5773502691896258*b_c[1]+0.25*b_r[0]-0.25*b_l[0])*dx1; 
  div_b[1] += ((-0.5*b_r[1])+0.5*b_l[1]+0.4330127018922193*(b_r[0]+b_l[0])-0.8660254037844386*b_c[0])*dx1; 

  bb_grad_u[0] += 0.7071067811865475*(bxbz[1]*grad_u_z[1]+bxby[1]*grad_u_y[1]+bxbx[1]*grad_u_x[1]+bxbz[0]*grad_u_z[0]+bxby[0]*grad_u_y[0]+bxbx[0]*grad_u_x[0])*dx1; 
  bb_grad_u[1] += 0.7071067811865475*(bxbz[0]*grad_u_z[1]+bxby[0]*grad_u_y[1]+bxbx[0]*grad_u_x[1]+grad_u_z[0]*bxbz[1]+grad_u_y[0]*bxby[1]+grad_u_x[0]*bxbx[1])*dx1; 

  div_p_x[0] += ((-0.2886751345948129*(Pxx_r[1]+Pxx_l[1]))+0.5773502691896258*Pxx_c[1]+0.25*Pxx_r[0]-0.25*Pxx_l[0])*dx1; 
  div_p_x[1] += ((-0.5*Pxx_r[1])+0.5*Pxx_l[1]+0.4330127018922193*(Pxx_r[0]+Pxx_l[0])-0.8660254037844386*Pxx_c[0])*dx1; 

  div_p_y[0] += ((-0.2886751345948129*(Pxy_r[1]+Pxy_l[1]))+0.5773502691896258*Pxy_c[1]+0.25*Pxy_r[0]-0.25*Pxy_l[0])*dx1; 
  div_p_y[1] += ((-0.5*Pxy_r[1])+0.5*Pxy_l[1]+0.4330127018922193*(Pxy_r[0]+Pxy_l[0])-0.8660254037844386*Pxy_c[0])*dx1; 

  div_p_z[0] += ((-0.2886751345948129*(Pxz_r[1]+Pxz_l[1]))+0.5773502691896258*Pxz_c[1]+0.25*Pxz_r[0]-0.25*Pxz_l[0])*dx1; 
  div_p_z[1] += ((-0.5*Pxz_r[1])+0.5*Pxz_l[1]+0.4330127018922193*(Pxz_r[0]+Pxz_l[0])-0.8660254037844386*Pxz_c[0])*dx1; 

} 
