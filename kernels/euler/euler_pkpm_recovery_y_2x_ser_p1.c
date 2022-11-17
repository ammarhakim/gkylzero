#include <gkyl_euler_kernels.h> 
GKYL_CU_DH void euler_pkpm_recovery_y_2x_ser_p1(const double *dxv, const double *bvarl, const double *bvarc, const double *bvarr, const double *u_il, const double *u_ic, const double *u_ir, const double *p_ijl, const double *p_ijc, const double *p_ijr, double* div_b, double* bb_grad_u, double* div_p) 
{ 
  // dxv[NDIM]: Cell spacing.
  // Al/Ac/Ar:  Inpute vector in left/center/right cells.
  // out:       Increment to volume expansion of div(A) in one direction.

  const double dx1 = 2.0/dxv[1]; 

  const double *b_l = &bvarl[4]; 
  const double *b_c = &bvarc[4]; 
  const double *b_r = &bvarr[4]; 

  const double *ux_l = &u_il[0]; 
  const double *uy_l = &u_il[4]; 
  const double *uz_l = &u_il[8]; 

  const double *ux_c = &u_ic[0]; 
  const double *uy_c = &u_ic[4]; 
  const double *uz_c = &u_ic[8]; 

  const double *ux_r = &u_ir[0]; 
  const double *uy_r = &u_ir[4]; 
  const double *uz_r = &u_ir[8]; 

  const double *bxbx = &bvarc[12]; 
  const double *bxby = &bvarc[16]; 
  const double *bxbz = &bvarc[20]; 
  const double *byby = &bvarc[24]; 
  const double *bybz = &bvarc[28]; 
  const double *bzbz = &bvarc[32]; 

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

  double *div_p_x = &div_p[0]; 
  double *div_p_y = &div_p[4]; 
  double *div_p_z = &div_p[8]; 
  double grad_u_x[4] = {0.0}; 
  double grad_u_y[4] = {0.0}; 
  double grad_u_z[4] = {0.0}; 
  grad_u_x[0] = (-0.2886751345948129*ux_r[2])-0.2886751345948129*ux_l[2]+0.5773502691896258*ux_c[2]+0.25*ux_r[0]-0.25*ux_l[0]; 
  grad_u_x[1] = (-0.2886751345948129*ux_r[3])-0.2886751345948129*ux_l[3]+0.5773502691896258*ux_c[3]+0.25*ux_r[1]-0.25*ux_l[1]; 
  grad_u_x[2] = (-0.5*ux_r[2])+0.5*ux_l[2]+0.4330127018922193*ux_r[0]+0.4330127018922193*ux_l[0]-0.8660254037844386*ux_c[0]; 
  grad_u_x[3] = (-0.5*ux_r[3])+0.5*ux_l[3]+0.4330127018922193*ux_r[1]+0.4330127018922193*ux_l[1]-0.8660254037844386*ux_c[1]; 

  grad_u_y[0] = (-0.2886751345948129*uy_r[2])-0.2886751345948129*uy_l[2]+0.5773502691896258*uy_c[2]+0.25*uy_r[0]-0.25*uy_l[0]; 
  grad_u_y[1] = (-0.2886751345948129*uy_r[3])-0.2886751345948129*uy_l[3]+0.5773502691896258*uy_c[3]+0.25*uy_r[1]-0.25*uy_l[1]; 
  grad_u_y[2] = (-0.5*uy_r[2])+0.5*uy_l[2]+0.4330127018922193*uy_r[0]+0.4330127018922193*uy_l[0]-0.8660254037844386*uy_c[0]; 
  grad_u_y[3] = (-0.5*uy_r[3])+0.5*uy_l[3]+0.4330127018922193*uy_r[1]+0.4330127018922193*uy_l[1]-0.8660254037844386*uy_c[1]; 

  grad_u_z[0] = (-0.2886751345948129*uz_r[2])-0.2886751345948129*uz_l[2]+0.5773502691896258*uz_c[2]+0.25*uz_r[0]-0.25*uz_l[0]; 
  grad_u_z[1] = (-0.2886751345948129*uz_r[3])-0.2886751345948129*uz_l[3]+0.5773502691896258*uz_c[3]+0.25*uz_r[1]-0.25*uz_l[1]; 
  grad_u_z[2] = (-0.5*uz_r[2])+0.5*uz_l[2]+0.4330127018922193*uz_r[0]+0.4330127018922193*uz_l[0]-0.8660254037844386*uz_c[0]; 
  grad_u_z[3] = (-0.5*uz_r[3])+0.5*uz_l[3]+0.4330127018922193*uz_r[1]+0.4330127018922193*uz_l[1]-0.8660254037844386*uz_c[1]; 

  div_b[0] += ((-0.2886751345948129*(b_r[2]+b_l[2]))+0.5773502691896258*b_c[2]+0.25*b_r[0]-0.25*b_l[0])*dx1; 
  div_b[1] += ((-0.2886751345948129*(b_r[3]+b_l[3]))+0.5773502691896258*b_c[3]+0.25*b_r[1]-0.25*b_l[1])*dx1; 
  div_b[2] += ((-0.5*b_r[2])+0.5*b_l[2]+0.4330127018922193*(b_r[0]+b_l[0])-0.8660254037844386*b_c[0])*dx1; 
  div_b[3] += ((-0.5*b_r[3])+0.5*b_l[3]+0.4330127018922193*(b_r[1]+b_l[1])-0.8660254037844386*b_c[1])*dx1; 

  bb_grad_u[0] += 0.5*(bybz[3]*grad_u_z[3]+byby[3]*grad_u_y[3]+bxby[3]*grad_u_x[3]+bybz[2]*grad_u_z[2]+byby[2]*grad_u_y[2]+bxby[2]*grad_u_x[2]+bybz[1]*grad_u_z[1]+byby[1]*grad_u_y[1]+bxby[1]*grad_u_x[1]+bybz[0]*grad_u_z[0]+byby[0]*grad_u_y[0]+bxby[0]*grad_u_x[0])*dx1; 
  bb_grad_u[1] += 0.5*(bybz[2]*grad_u_z[3]+byby[2]*grad_u_y[3]+bxby[2]*grad_u_x[3]+grad_u_z[2]*bybz[3]+grad_u_y[2]*byby[3]+grad_u_x[2]*bxby[3]+bybz[0]*grad_u_z[1]+byby[0]*grad_u_y[1]+bxby[0]*grad_u_x[1]+grad_u_z[0]*bybz[1]+grad_u_y[0]*byby[1]+grad_u_x[0]*bxby[1])*dx1; 
  bb_grad_u[2] += 0.5*(bybz[1]*grad_u_z[3]+byby[1]*grad_u_y[3]+bxby[1]*grad_u_x[3]+grad_u_z[1]*bybz[3]+grad_u_y[1]*byby[3]+grad_u_x[1]*bxby[3]+bybz[0]*grad_u_z[2]+byby[0]*grad_u_y[2]+bxby[0]*grad_u_x[2]+grad_u_z[0]*bybz[2]+grad_u_y[0]*byby[2]+grad_u_x[0]*bxby[2])*dx1; 
  bb_grad_u[3] += 0.5*(bybz[0]*grad_u_z[3]+byby[0]*grad_u_y[3]+bxby[0]*grad_u_x[3]+grad_u_z[0]*bybz[3]+grad_u_y[0]*byby[3]+grad_u_x[0]*bxby[3]+bybz[1]*grad_u_z[2]+byby[1]*grad_u_y[2]+bxby[1]*grad_u_x[2]+grad_u_z[1]*bybz[2]+grad_u_y[1]*byby[2]+grad_u_x[1]*bxby[2])*dx1; 

  div_p_x[0] += ((-0.2886751345948129*(Pxy_r[2]+Pxy_l[2]))+0.5773502691896258*Pxy_c[2]+0.25*Pxy_r[0]-0.25*Pxy_l[0])*dx1; 
  div_p_x[1] += ((-0.2886751345948129*(Pxy_r[3]+Pxy_l[3]))+0.5773502691896258*Pxy_c[3]+0.25*Pxy_r[1]-0.25*Pxy_l[1])*dx1; 
  div_p_x[2] += ((-0.5*Pxy_r[2])+0.5*Pxy_l[2]+0.4330127018922193*(Pxy_r[0]+Pxy_l[0])-0.8660254037844386*Pxy_c[0])*dx1; 
  div_p_x[3] += ((-0.5*Pxy_r[3])+0.5*Pxy_l[3]+0.4330127018922193*(Pxy_r[1]+Pxy_l[1])-0.8660254037844386*Pxy_c[1])*dx1; 

  div_p_y[0] += ((-0.2886751345948129*(Pyy_r[2]+Pyy_l[2]))+0.5773502691896258*Pyy_c[2]+0.25*Pyy_r[0]-0.25*Pyy_l[0])*dx1; 
  div_p_y[1] += ((-0.2886751345948129*(Pyy_r[3]+Pyy_l[3]))+0.5773502691896258*Pyy_c[3]+0.25*Pyy_r[1]-0.25*Pyy_l[1])*dx1; 
  div_p_y[2] += ((-0.5*Pyy_r[2])+0.5*Pyy_l[2]+0.4330127018922193*(Pyy_r[0]+Pyy_l[0])-0.8660254037844386*Pyy_c[0])*dx1; 
  div_p_y[3] += ((-0.5*Pyy_r[3])+0.5*Pyy_l[3]+0.4330127018922193*(Pyy_r[1]+Pyy_l[1])-0.8660254037844386*Pyy_c[1])*dx1; 

  div_p_z[0] += ((-0.2886751345948129*(Pyz_r[2]+Pyz_l[2]))+0.5773502691896258*Pyz_c[2]+0.25*Pyz_r[0]-0.25*Pyz_l[0])*dx1; 
  div_p_z[1] += ((-0.2886751345948129*(Pyz_r[3]+Pyz_l[3]))+0.5773502691896258*Pyz_c[3]+0.25*Pyz_r[1]-0.25*Pyz_l[1])*dx1; 
  div_p_z[2] += ((-0.5*Pyz_r[2])+0.5*Pyz_l[2]+0.4330127018922193*(Pyz_r[0]+Pyz_l[0])-0.8660254037844386*Pyz_c[0])*dx1; 
  div_p_z[3] += ((-0.5*Pyz_r[3])+0.5*Pyz_l[3]+0.4330127018922193*(Pyz_r[1]+Pyz_l[1])-0.8660254037844386*Pyz_c[1])*dx1; 

} 
