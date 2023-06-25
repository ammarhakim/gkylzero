#include <gkyl_euler_pkpm_kernels.h> 
GKYL_CU_DH void euler_pkpm_recovery_y_3x_ser_p1(const double *dxv, 
  const double *bvarl, const double *bvarc, const double *bvarr, 
  const double *u_il, const double *u_ic, const double *u_ir, 
  const double *p_ijl, const double *p_ijc, const double *p_ijr, 
  const double *pkpm_div_ppar, const double *rho_inv, const double *T_perp_over_m, const double *nu, 
  double* div_p, double* pkpm_accel_vars) 
{ 
  // dxv[NDIM]:             Cell spacing.
  // bvarl/c/r:             Input magnetic field unit vector in left/center/right cells.
  // u_il/c/r:              Input bulk velocity (ux,uy,uz) in left/center/right cells.
  // p_ijl/c/r:             Input pressure tensor in left/center/right cells.
  // pkpm_div_ppar:         Input div(p_parallel b_hat) in center cell.
  // rho_inv:               Input 1/rho in center cell.
  // T_perp_over_m:         Input p_perp/rho = T_perp/m in center cell.
  // nu:                    Input collisionality in center cell.
  // div_p:                 Increment to volume expansion of div(p) in one direction; includes hyper-diffusion for momentum.
  // pkpm_accel_vars:       Increment to volume expansion of pkpm acceleration variables.

  const double dx1 = 2.0/dxv[1]; 
  const double *b_l = &bvarl[8]; 
  const double *b_c = &bvarc[8]; 
  const double *b_r = &bvarr[8]; 

  const double *pkpm_div_ppar_c = &pkpm_div_ppar[8]; 

  const double *ux_l = &u_il[0]; 
  const double *uy_l = &u_il[8]; 
  const double *uz_l = &u_il[16]; 

  const double *ux_c = &u_ic[0]; 
  const double *uy_c = &u_ic[8]; 
  const double *uz_c = &u_ic[16]; 

  const double *ux_r = &u_ir[0]; 
  const double *uy_r = &u_ir[8]; 
  const double *uz_r = &u_ir[16]; 

  const double *bxbx = &bvarc[24]; 
  const double *bxby = &bvarc[32]; 
  const double *bxbz = &bvarc[40]; 
  const double *byby = &bvarc[48]; 
  const double *bybz = &bvarc[56]; 
  const double *bzbz = &bvarc[64]; 

  const double *Pxx_l = &p_ijl[0]; 
  const double *Pxy_l = &p_ijl[8]; 
  const double *Pxz_l = &p_ijl[16]; 
  const double *Pyy_l = &p_ijl[24]; 
  const double *Pyz_l = &p_ijl[32]; 
  const double *Pzz_l = &p_ijl[40]; 

  const double *Pxx_c = &p_ijc[0]; 
  const double *Pxy_c = &p_ijc[8]; 
  const double *Pxz_c = &p_ijc[16]; 
  const double *Pyy_c = &p_ijc[24]; 
  const double *Pyz_c = &p_ijc[32]; 
  const double *Pzz_c = &p_ijc[40]; 

  const double *Pxx_r = &p_ijr[0]; 
  const double *Pxy_r = &p_ijr[8]; 
  const double *Pxz_r = &p_ijr[16]; 
  const double *Pyy_r = &p_ijr[24]; 
  const double *Pyz_r = &p_ijr[32]; 
  const double *Pzz_r = &p_ijr[40]; 

  double *div_p_x = &div_p[0]; 
  double *div_p_y = &div_p[8]; 
  double *div_p_z = &div_p[16]; 

  double *div_b = &pkpm_accel_vars[0]; 
  double *bb_grad_u = &pkpm_accel_vars[8]; 
  double *p_force = &pkpm_accel_vars[16]; 
  double *p_perp_source = &pkpm_accel_vars[24]; 
  double *p_perp_div_b = &pkpm_accel_vars[32]; 

  double grad_u_x[8] = {0.0}; 
  double grad_u_y[8] = {0.0}; 
  double grad_u_z[8] = {0.0}; 
  grad_u_x[0] = (-0.4330127018922193*ux_r[2]*dx1)-0.4330127018922193*ux_l[2]*dx1+0.8660254037844386*ux_c[2]*dx1+0.25*ux_r[0]*dx1-0.25*ux_l[0]*dx1; 
  grad_u_x[1] = (-0.4330127018922193*ux_r[4]*dx1)-0.4330127018922193*ux_l[4]*dx1+0.8660254037844386*ux_c[4]*dx1+0.25*ux_r[1]*dx1-0.25*ux_l[1]*dx1; 
  grad_u_x[2] = (-0.75*ux_r[2]*dx1)+0.75*ux_l[2]*dx1+0.4330127018922193*ux_r[0]*dx1+0.4330127018922193*ux_l[0]*dx1-0.8660254037844385*ux_c[0]*dx1; 
  grad_u_x[3] = (-0.4330127018922193*ux_r[6]*dx1)-0.4330127018922193*ux_l[6]*dx1+0.8660254037844386*ux_c[6]*dx1+0.25*ux_r[3]*dx1-0.25*ux_l[3]*dx1; 
  grad_u_x[4] = (-0.75*ux_r[4]*dx1)+0.75*ux_l[4]*dx1+0.4330127018922193*ux_r[1]*dx1+0.4330127018922193*ux_l[1]*dx1-0.8660254037844385*ux_c[1]*dx1; 
  grad_u_x[5] = (-0.4330127018922193*ux_r[7]*dx1)-0.4330127018922193*ux_l[7]*dx1+0.8660254037844386*ux_c[7]*dx1+0.25*ux_r[5]*dx1-0.25*ux_l[5]*dx1; 
  grad_u_x[6] = (-0.75*ux_r[6]*dx1)+0.75*ux_l[6]*dx1+0.4330127018922193*ux_r[3]*dx1+0.4330127018922193*ux_l[3]*dx1-0.8660254037844385*ux_c[3]*dx1; 
  grad_u_x[7] = (-0.75*ux_r[7]*dx1)+0.75*ux_l[7]*dx1+0.4330127018922193*ux_r[5]*dx1+0.4330127018922193*ux_l[5]*dx1-0.8660254037844385*ux_c[5]*dx1; 

  grad_u_y[0] = (-0.4330127018922193*uy_r[2]*dx1)-0.4330127018922193*uy_l[2]*dx1+0.8660254037844386*uy_c[2]*dx1+0.25*uy_r[0]*dx1-0.25*uy_l[0]*dx1; 
  grad_u_y[1] = (-0.4330127018922193*uy_r[4]*dx1)-0.4330127018922193*uy_l[4]*dx1+0.8660254037844386*uy_c[4]*dx1+0.25*uy_r[1]*dx1-0.25*uy_l[1]*dx1; 
  grad_u_y[2] = (-0.75*uy_r[2]*dx1)+0.75*uy_l[2]*dx1+0.4330127018922193*uy_r[0]*dx1+0.4330127018922193*uy_l[0]*dx1-0.8660254037844385*uy_c[0]*dx1; 
  grad_u_y[3] = (-0.4330127018922193*uy_r[6]*dx1)-0.4330127018922193*uy_l[6]*dx1+0.8660254037844386*uy_c[6]*dx1+0.25*uy_r[3]*dx1-0.25*uy_l[3]*dx1; 
  grad_u_y[4] = (-0.75*uy_r[4]*dx1)+0.75*uy_l[4]*dx1+0.4330127018922193*uy_r[1]*dx1+0.4330127018922193*uy_l[1]*dx1-0.8660254037844385*uy_c[1]*dx1; 
  grad_u_y[5] = (-0.4330127018922193*uy_r[7]*dx1)-0.4330127018922193*uy_l[7]*dx1+0.8660254037844386*uy_c[7]*dx1+0.25*uy_r[5]*dx1-0.25*uy_l[5]*dx1; 
  grad_u_y[6] = (-0.75*uy_r[6]*dx1)+0.75*uy_l[6]*dx1+0.4330127018922193*uy_r[3]*dx1+0.4330127018922193*uy_l[3]*dx1-0.8660254037844385*uy_c[3]*dx1; 
  grad_u_y[7] = (-0.75*uy_r[7]*dx1)+0.75*uy_l[7]*dx1+0.4330127018922193*uy_r[5]*dx1+0.4330127018922193*uy_l[5]*dx1-0.8660254037844385*uy_c[5]*dx1; 

  grad_u_z[0] = (-0.4330127018922193*uz_r[2]*dx1)-0.4330127018922193*uz_l[2]*dx1+0.8660254037844386*uz_c[2]*dx1+0.25*uz_r[0]*dx1-0.25*uz_l[0]*dx1; 
  grad_u_z[1] = (-0.4330127018922193*uz_r[4]*dx1)-0.4330127018922193*uz_l[4]*dx1+0.8660254037844386*uz_c[4]*dx1+0.25*uz_r[1]*dx1-0.25*uz_l[1]*dx1; 
  grad_u_z[2] = (-0.75*uz_r[2]*dx1)+0.75*uz_l[2]*dx1+0.4330127018922193*uz_r[0]*dx1+0.4330127018922193*uz_l[0]*dx1-0.8660254037844385*uz_c[0]*dx1; 
  grad_u_z[3] = (-0.4330127018922193*uz_r[6]*dx1)-0.4330127018922193*uz_l[6]*dx1+0.8660254037844386*uz_c[6]*dx1+0.25*uz_r[3]*dx1-0.25*uz_l[3]*dx1; 
  grad_u_z[4] = (-0.75*uz_r[4]*dx1)+0.75*uz_l[4]*dx1+0.4330127018922193*uz_r[1]*dx1+0.4330127018922193*uz_l[1]*dx1-0.8660254037844385*uz_c[1]*dx1; 
  grad_u_z[5] = (-0.4330127018922193*uz_r[7]*dx1)-0.4330127018922193*uz_l[7]*dx1+0.8660254037844386*uz_c[7]*dx1+0.25*uz_r[5]*dx1-0.25*uz_l[5]*dx1; 
  grad_u_z[6] = (-0.75*uz_r[6]*dx1)+0.75*uz_l[6]*dx1+0.4330127018922193*uz_r[3]*dx1+0.4330127018922193*uz_l[3]*dx1-0.8660254037844385*uz_c[3]*dx1; 
  grad_u_z[7] = (-0.75*uz_r[7]*dx1)+0.75*uz_l[7]*dx1+0.4330127018922193*uz_r[5]*dx1+0.4330127018922193*uz_l[5]*dx1-0.8660254037844385*uz_c[5]*dx1; 

  double div_b_comp[8] = {0.0}; 
  double div_p_x_comp[8] = {0.0}; 
  double div_p_y_comp[8] = {0.0}; 
  double div_p_z_comp[8] = {0.0}; 
  double bb_grad_u_comp[8] = {0.0}; 
  div_b_comp[0] = (-0.2886751345948129*b_r[2]*dx1)-0.2886751345948129*b_l[2]*dx1+0.5773502691896258*b_c[2]*dx1+0.25*b_r[0]*dx1-0.25*b_l[0]*dx1; 
  div_b_comp[1] = (-0.2886751345948129*b_r[4]*dx1)-0.2886751345948129*b_l[4]*dx1+0.5773502691896258*b_c[4]*dx1+0.25*b_r[1]*dx1-0.25*b_l[1]*dx1; 
  div_b_comp[2] = (-0.5*b_r[2]*dx1)+0.5*b_l[2]*dx1+0.4330127018922193*b_r[0]*dx1+0.4330127018922193*b_l[0]*dx1-0.8660254037844386*b_c[0]*dx1; 
  div_b_comp[3] = (-0.2886751345948129*b_r[6]*dx1)-0.2886751345948129*b_l[6]*dx1+0.5773502691896258*b_c[6]*dx1+0.25*b_r[3]*dx1-0.25*b_l[3]*dx1; 
  div_b_comp[4] = (-0.5*b_r[4]*dx1)+0.5*b_l[4]*dx1+0.4330127018922193*b_r[1]*dx1+0.4330127018922193*b_l[1]*dx1-0.8660254037844386*b_c[1]*dx1; 
  div_b_comp[5] = (-0.2886751345948129*b_r[7]*dx1)-0.2886751345948129*b_l[7]*dx1+0.5773502691896258*b_c[7]*dx1+0.25*b_r[5]*dx1-0.25*b_l[5]*dx1; 
  div_b_comp[6] = (-0.5*b_r[6]*dx1)+0.5*b_l[6]*dx1+0.4330127018922193*b_r[3]*dx1+0.4330127018922193*b_l[3]*dx1-0.8660254037844386*b_c[3]*dx1; 
  div_b_comp[7] = (-0.5*b_r[7]*dx1)+0.5*b_l[7]*dx1+0.4330127018922193*b_r[5]*dx1+0.4330127018922193*b_l[5]*dx1-0.8660254037844386*b_c[5]*dx1; 

  div_b[0] += div_b_comp[0]; 
  div_b[1] += div_b_comp[1]; 
  div_b[2] += div_b_comp[2]; 
  div_b[3] += div_b_comp[3]; 
  div_b[4] += div_b_comp[4]; 
  div_b[5] += div_b_comp[5]; 
  div_b[6] += div_b_comp[6]; 
  div_b[7] += div_b_comp[7]; 

  div_p_x_comp[0] += ((-0.4330127018922193*(Pxy_r[2]+Pxy_l[2]))+0.8660254037844386*Pxy_c[2]+0.25*Pxy_r[0]-0.25*Pxy_l[0])*dx1; 
  div_p_x_comp[1] += ((-0.4330127018922193*(Pxy_r[4]+Pxy_l[4]))+0.8660254037844386*Pxy_c[4]+0.25*Pxy_r[1]-0.25*Pxy_l[1])*dx1; 
  div_p_x_comp[2] += ((-0.75*Pxy_r[2])+0.75*Pxy_l[2]+0.4330127018922193*(Pxy_r[0]+Pxy_l[0])-0.8660254037844385*Pxy_c[0])*dx1; 
  div_p_x_comp[3] += ((-0.4330127018922193*(Pxy_r[6]+Pxy_l[6]))+0.8660254037844386*Pxy_c[6]+0.25*Pxy_r[3]-0.25*Pxy_l[3])*dx1; 
  div_p_x_comp[4] += ((-0.75*Pxy_r[4])+0.75*Pxy_l[4]+0.4330127018922193*(Pxy_r[1]+Pxy_l[1])-0.8660254037844385*Pxy_c[1])*dx1; 
  div_p_x_comp[5] += ((-0.4330127018922193*(Pxy_r[7]+Pxy_l[7]))+0.8660254037844386*Pxy_c[7]+0.25*Pxy_r[5]-0.25*Pxy_l[5])*dx1; 
  div_p_x_comp[6] += ((-0.75*Pxy_r[6])+0.75*Pxy_l[6]+0.4330127018922193*(Pxy_r[3]+Pxy_l[3])-0.8660254037844385*Pxy_c[3])*dx1; 
  div_p_x_comp[7] += ((-0.75*Pxy_r[7])+0.75*Pxy_l[7]+0.4330127018922193*(Pxy_r[5]+Pxy_l[5])-0.8660254037844385*Pxy_c[5])*dx1; 

  div_p_y_comp[0] += ((-0.4330127018922193*(Pyy_r[2]+Pyy_l[2]))+0.8660254037844386*Pyy_c[2]+0.25*Pyy_r[0]-0.25*Pyy_l[0])*dx1; 
  div_p_y_comp[1] += ((-0.4330127018922193*(Pyy_r[4]+Pyy_l[4]))+0.8660254037844386*Pyy_c[4]+0.25*Pyy_r[1]-0.25*Pyy_l[1])*dx1; 
  div_p_y_comp[2] += ((-0.75*Pyy_r[2])+0.75*Pyy_l[2]+0.4330127018922193*(Pyy_r[0]+Pyy_l[0])-0.8660254037844385*Pyy_c[0])*dx1; 
  div_p_y_comp[3] += ((-0.4330127018922193*(Pyy_r[6]+Pyy_l[6]))+0.8660254037844386*Pyy_c[6]+0.25*Pyy_r[3]-0.25*Pyy_l[3])*dx1; 
  div_p_y_comp[4] += ((-0.75*Pyy_r[4])+0.75*Pyy_l[4]+0.4330127018922193*(Pyy_r[1]+Pyy_l[1])-0.8660254037844385*Pyy_c[1])*dx1; 
  div_p_y_comp[5] += ((-0.4330127018922193*(Pyy_r[7]+Pyy_l[7]))+0.8660254037844386*Pyy_c[7]+0.25*Pyy_r[5]-0.25*Pyy_l[5])*dx1; 
  div_p_y_comp[6] += ((-0.75*Pyy_r[6])+0.75*Pyy_l[6]+0.4330127018922193*(Pyy_r[3]+Pyy_l[3])-0.8660254037844385*Pyy_c[3])*dx1; 
  div_p_y_comp[7] += ((-0.75*Pyy_r[7])+0.75*Pyy_l[7]+0.4330127018922193*(Pyy_r[5]+Pyy_l[5])-0.8660254037844385*Pyy_c[5])*dx1; 

  div_p_z_comp[0] += ((-0.4330127018922193*(Pyz_r[2]+Pyz_l[2]))+0.8660254037844386*Pyz_c[2]+0.25*Pyz_r[0]-0.25*Pyz_l[0])*dx1; 
  div_p_z_comp[1] += ((-0.4330127018922193*(Pyz_r[4]+Pyz_l[4]))+0.8660254037844386*Pyz_c[4]+0.25*Pyz_r[1]-0.25*Pyz_l[1])*dx1; 
  div_p_z_comp[2] += ((-0.75*Pyz_r[2])+0.75*Pyz_l[2]+0.4330127018922193*(Pyz_r[0]+Pyz_l[0])-0.8660254037844385*Pyz_c[0])*dx1; 
  div_p_z_comp[3] += ((-0.4330127018922193*(Pyz_r[6]+Pyz_l[6]))+0.8660254037844386*Pyz_c[6]+0.25*Pyz_r[3]-0.25*Pyz_l[3])*dx1; 
  div_p_z_comp[4] += ((-0.75*Pyz_r[4])+0.75*Pyz_l[4]+0.4330127018922193*(Pyz_r[1]+Pyz_l[1])-0.8660254037844385*Pyz_c[1])*dx1; 
  div_p_z_comp[5] += ((-0.4330127018922193*(Pyz_r[7]+Pyz_l[7]))+0.8660254037844386*Pyz_c[7]+0.25*Pyz_r[5]-0.25*Pyz_l[5])*dx1; 
  div_p_z_comp[6] += ((-0.75*Pyz_r[6])+0.75*Pyz_l[6]+0.4330127018922193*(Pyz_r[3]+Pyz_l[3])-0.8660254037844385*Pyz_c[3])*dx1; 
  div_p_z_comp[7] += ((-0.75*Pyz_r[7])+0.75*Pyz_l[7]+0.4330127018922193*(Pyz_r[5]+Pyz_l[5])-0.8660254037844385*Pyz_c[5])*dx1; 

  bb_grad_u_comp[0] = 0.3535533905932737*bybz[7]*grad_u_z[7]+0.3535533905932737*byby[7]*grad_u_y[7]+0.3535533905932737*bxby[7]*grad_u_x[7]+0.3535533905932737*bybz[6]*grad_u_z[6]+0.3535533905932737*byby[6]*grad_u_y[6]+0.3535533905932737*bxby[6]*grad_u_x[6]+0.3535533905932737*bybz[5]*grad_u_z[5]+0.3535533905932737*byby[5]*grad_u_y[5]+0.3535533905932737*bxby[5]*grad_u_x[5]+0.3535533905932737*bybz[4]*grad_u_z[4]+0.3535533905932737*byby[4]*grad_u_y[4]+0.3535533905932737*bxby[4]*grad_u_x[4]+0.3535533905932737*bybz[3]*grad_u_z[3]+0.3535533905932737*byby[3]*grad_u_y[3]+0.3535533905932737*bxby[3]*grad_u_x[3]+0.3535533905932737*bybz[2]*grad_u_z[2]+0.3535533905932737*byby[2]*grad_u_y[2]+0.3535533905932737*bxby[2]*grad_u_x[2]+0.3535533905932737*bybz[1]*grad_u_z[1]+0.3535533905932737*byby[1]*grad_u_y[1]+0.3535533905932737*bxby[1]*grad_u_x[1]+0.3535533905932737*bybz[0]*grad_u_z[0]+0.3535533905932737*byby[0]*grad_u_y[0]+0.3535533905932737*bxby[0]*grad_u_x[0]; 
  bb_grad_u_comp[1] = 0.3535533905932737*bybz[6]*grad_u_z[7]+0.3535533905932737*byby[6]*grad_u_y[7]+0.3535533905932737*bxby[6]*grad_u_x[7]+0.3535533905932737*grad_u_z[6]*bybz[7]+0.3535533905932737*grad_u_y[6]*byby[7]+0.3535533905932737*grad_u_x[6]*bxby[7]+0.3535533905932737*bybz[3]*grad_u_z[5]+0.3535533905932737*byby[3]*grad_u_y[5]+0.3535533905932737*bxby[3]*grad_u_x[5]+0.3535533905932737*grad_u_z[3]*bybz[5]+0.3535533905932737*grad_u_y[3]*byby[5]+0.3535533905932737*grad_u_x[3]*bxby[5]+0.3535533905932737*bybz[2]*grad_u_z[4]+0.3535533905932737*byby[2]*grad_u_y[4]+0.3535533905932737*bxby[2]*grad_u_x[4]+0.3535533905932737*grad_u_z[2]*bybz[4]+0.3535533905932737*grad_u_y[2]*byby[4]+0.3535533905932737*grad_u_x[2]*bxby[4]+0.3535533905932737*bybz[0]*grad_u_z[1]+0.3535533905932737*byby[0]*grad_u_y[1]+0.3535533905932737*bxby[0]*grad_u_x[1]+0.3535533905932737*grad_u_z[0]*bybz[1]+0.3535533905932737*grad_u_y[0]*byby[1]+0.3535533905932737*grad_u_x[0]*bxby[1]; 
  bb_grad_u_comp[2] = 0.3535533905932737*bybz[5]*grad_u_z[7]+0.3535533905932737*byby[5]*grad_u_y[7]+0.3535533905932737*bxby[5]*grad_u_x[7]+0.3535533905932737*grad_u_z[5]*bybz[7]+0.3535533905932737*grad_u_y[5]*byby[7]+0.3535533905932737*grad_u_x[5]*bxby[7]+0.3535533905932737*bybz[3]*grad_u_z[6]+0.3535533905932737*byby[3]*grad_u_y[6]+0.3535533905932737*bxby[3]*grad_u_x[6]+0.3535533905932737*grad_u_z[3]*bybz[6]+0.3535533905932737*grad_u_y[3]*byby[6]+0.3535533905932737*grad_u_x[3]*bxby[6]+0.3535533905932737*bybz[1]*grad_u_z[4]+0.3535533905932737*byby[1]*grad_u_y[4]+0.3535533905932737*bxby[1]*grad_u_x[4]+0.3535533905932737*grad_u_z[1]*bybz[4]+0.3535533905932737*grad_u_y[1]*byby[4]+0.3535533905932737*grad_u_x[1]*bxby[4]+0.3535533905932737*bybz[0]*grad_u_z[2]+0.3535533905932737*byby[0]*grad_u_y[2]+0.3535533905932737*bxby[0]*grad_u_x[2]+0.3535533905932737*grad_u_z[0]*bybz[2]+0.3535533905932737*grad_u_y[0]*byby[2]+0.3535533905932737*grad_u_x[0]*bxby[2]; 
  bb_grad_u_comp[3] = 0.3535533905932737*bybz[4]*grad_u_z[7]+0.3535533905932737*byby[4]*grad_u_y[7]+0.3535533905932737*bxby[4]*grad_u_x[7]+0.3535533905932737*grad_u_z[4]*bybz[7]+0.3535533905932737*grad_u_y[4]*byby[7]+0.3535533905932737*grad_u_x[4]*bxby[7]+0.3535533905932737*bybz[2]*grad_u_z[6]+0.3535533905932737*byby[2]*grad_u_y[6]+0.3535533905932737*bxby[2]*grad_u_x[6]+0.3535533905932737*grad_u_z[2]*bybz[6]+0.3535533905932737*grad_u_y[2]*byby[6]+0.3535533905932737*grad_u_x[2]*bxby[6]+0.3535533905932737*bybz[1]*grad_u_z[5]+0.3535533905932737*byby[1]*grad_u_y[5]+0.3535533905932737*bxby[1]*grad_u_x[5]+0.3535533905932737*grad_u_z[1]*bybz[5]+0.3535533905932737*grad_u_y[1]*byby[5]+0.3535533905932737*grad_u_x[1]*bxby[5]+0.3535533905932737*bybz[0]*grad_u_z[3]+0.3535533905932737*byby[0]*grad_u_y[3]+0.3535533905932737*bxby[0]*grad_u_x[3]+0.3535533905932737*grad_u_z[0]*bybz[3]+0.3535533905932737*grad_u_y[0]*byby[3]+0.3535533905932737*grad_u_x[0]*bxby[3]; 
  bb_grad_u_comp[4] = 0.3535533905932737*bybz[3]*grad_u_z[7]+0.3535533905932737*byby[3]*grad_u_y[7]+0.3535533905932737*bxby[3]*grad_u_x[7]+0.3535533905932737*grad_u_z[3]*bybz[7]+0.3535533905932737*grad_u_y[3]*byby[7]+0.3535533905932737*grad_u_x[3]*bxby[7]+0.3535533905932737*bybz[5]*grad_u_z[6]+0.3535533905932737*byby[5]*grad_u_y[6]+0.3535533905932737*bxby[5]*grad_u_x[6]+0.3535533905932737*grad_u_z[5]*bybz[6]+0.3535533905932737*grad_u_y[5]*byby[6]+0.3535533905932737*grad_u_x[5]*bxby[6]+0.3535533905932737*bybz[0]*grad_u_z[4]+0.3535533905932737*byby[0]*grad_u_y[4]+0.3535533905932737*bxby[0]*grad_u_x[4]+0.3535533905932737*grad_u_z[0]*bybz[4]+0.3535533905932737*grad_u_y[0]*byby[4]+0.3535533905932737*grad_u_x[0]*bxby[4]+0.3535533905932737*bybz[1]*grad_u_z[2]+0.3535533905932737*byby[1]*grad_u_y[2]+0.3535533905932737*bxby[1]*grad_u_x[2]+0.3535533905932737*grad_u_z[1]*bybz[2]+0.3535533905932737*grad_u_y[1]*byby[2]+0.3535533905932737*grad_u_x[1]*bxby[2]; 
  bb_grad_u_comp[5] = 0.3535533905932737*bybz[2]*grad_u_z[7]+0.3535533905932737*byby[2]*grad_u_y[7]+0.3535533905932737*bxby[2]*grad_u_x[7]+0.3535533905932737*grad_u_z[2]*bybz[7]+0.3535533905932737*grad_u_y[2]*byby[7]+0.3535533905932737*grad_u_x[2]*bxby[7]+0.3535533905932737*bybz[4]*grad_u_z[6]+0.3535533905932737*byby[4]*grad_u_y[6]+0.3535533905932737*bxby[4]*grad_u_x[6]+0.3535533905932737*grad_u_z[4]*bybz[6]+0.3535533905932737*grad_u_y[4]*byby[6]+0.3535533905932737*grad_u_x[4]*bxby[6]+0.3535533905932737*bybz[0]*grad_u_z[5]+0.3535533905932737*byby[0]*grad_u_y[5]+0.3535533905932737*bxby[0]*grad_u_x[5]+0.3535533905932737*grad_u_z[0]*bybz[5]+0.3535533905932737*grad_u_y[0]*byby[5]+0.3535533905932737*grad_u_x[0]*bxby[5]+0.3535533905932737*bybz[1]*grad_u_z[3]+0.3535533905932737*byby[1]*grad_u_y[3]+0.3535533905932737*bxby[1]*grad_u_x[3]+0.3535533905932737*grad_u_z[1]*bybz[3]+0.3535533905932737*grad_u_y[1]*byby[3]+0.3535533905932737*grad_u_x[1]*bxby[3]; 
  bb_grad_u_comp[6] = 0.3535533905932737*bybz[1]*grad_u_z[7]+0.3535533905932737*byby[1]*grad_u_y[7]+0.3535533905932737*bxby[1]*grad_u_x[7]+0.3535533905932737*grad_u_z[1]*bybz[7]+0.3535533905932737*grad_u_y[1]*byby[7]+0.3535533905932737*grad_u_x[1]*bxby[7]+0.3535533905932737*bybz[0]*grad_u_z[6]+0.3535533905932737*byby[0]*grad_u_y[6]+0.3535533905932737*bxby[0]*grad_u_x[6]+0.3535533905932737*grad_u_z[0]*bybz[6]+0.3535533905932737*grad_u_y[0]*byby[6]+0.3535533905932737*grad_u_x[0]*bxby[6]+0.3535533905932737*bybz[4]*grad_u_z[5]+0.3535533905932737*byby[4]*grad_u_y[5]+0.3535533905932737*bxby[4]*grad_u_x[5]+0.3535533905932737*grad_u_z[4]*bybz[5]+0.3535533905932737*grad_u_y[4]*byby[5]+0.3535533905932737*grad_u_x[4]*bxby[5]+0.3535533905932737*bybz[2]*grad_u_z[3]+0.3535533905932737*byby[2]*grad_u_y[3]+0.3535533905932737*bxby[2]*grad_u_x[3]+0.3535533905932737*grad_u_z[2]*bybz[3]+0.3535533905932737*grad_u_y[2]*byby[3]+0.3535533905932737*grad_u_x[2]*bxby[3]; 
  bb_grad_u_comp[7] = 0.3535533905932737*bybz[0]*grad_u_z[7]+0.3535533905932737*byby[0]*grad_u_y[7]+0.3535533905932737*bxby[0]*grad_u_x[7]+0.3535533905932737*grad_u_z[0]*bybz[7]+0.3535533905932737*grad_u_y[0]*byby[7]+0.3535533905932737*grad_u_x[0]*bxby[7]+0.3535533905932737*bybz[1]*grad_u_z[6]+0.3535533905932737*byby[1]*grad_u_y[6]+0.3535533905932737*bxby[1]*grad_u_x[6]+0.3535533905932737*grad_u_z[1]*bybz[6]+0.3535533905932737*grad_u_y[1]*byby[6]+0.3535533905932737*grad_u_x[1]*bxby[6]+0.3535533905932737*bybz[2]*grad_u_z[5]+0.3535533905932737*byby[2]*grad_u_y[5]+0.3535533905932737*bxby[2]*grad_u_x[5]+0.3535533905932737*grad_u_z[2]*bybz[5]+0.3535533905932737*grad_u_y[2]*byby[5]+0.3535533905932737*grad_u_x[2]*bxby[5]+0.3535533905932737*bybz[3]*grad_u_z[4]+0.3535533905932737*byby[3]*grad_u_y[4]+0.3535533905932737*bxby[3]*grad_u_x[4]+0.3535533905932737*grad_u_z[3]*bybz[4]+0.3535533905932737*grad_u_y[3]*byby[4]+0.3535533905932737*grad_u_x[3]*bxby[4]; 

  bb_grad_u[0] += bb_grad_u_comp[0]; 
  bb_grad_u[1] += bb_grad_u_comp[1]; 
  bb_grad_u[2] += bb_grad_u_comp[2]; 
  bb_grad_u[3] += bb_grad_u_comp[3]; 
  bb_grad_u[4] += bb_grad_u_comp[4]; 
  bb_grad_u[5] += bb_grad_u_comp[5]; 
  bb_grad_u[6] += bb_grad_u_comp[6]; 
  bb_grad_u[7] += bb_grad_u_comp[7]; 

  div_p_x[0] += div_p_x_comp[0]; 
  div_p_x[1] += div_p_x_comp[1]; 
  div_p_x[2] += div_p_x_comp[2]; 
  div_p_x[3] += div_p_x_comp[3]; 
  div_p_x[4] += div_p_x_comp[4]; 
  div_p_x[5] += div_p_x_comp[5]; 
  div_p_x[6] += div_p_x_comp[6]; 
  div_p_x[7] += div_p_x_comp[7]; 

  div_p_y[0] += div_p_y_comp[0]; 
  div_p_y[1] += div_p_y_comp[1]; 
  div_p_y[2] += div_p_y_comp[2]; 
  div_p_y[3] += div_p_y_comp[3]; 
  div_p_y[4] += div_p_y_comp[4]; 
  div_p_y[5] += div_p_y_comp[5]; 
  div_p_y[6] += div_p_y_comp[6]; 
  div_p_y[7] += div_p_y_comp[7]; 

  div_p_z[0] += div_p_z_comp[0]; 
  div_p_z[1] += div_p_z_comp[1]; 
  div_p_z[2] += div_p_z_comp[2]; 
  div_p_z[3] += div_p_z_comp[3]; 
  div_p_z[4] += div_p_z_comp[4]; 
  div_p_z[5] += div_p_z_comp[5]; 
  div_p_z[6] += div_p_z_comp[6]; 
  div_p_z[7] += div_p_z_comp[7]; 

  p_force[0] += 0.3535533905932737*pkpm_div_ppar_c[7]*rho_inv[7]-0.3535533905932737*T_perp_over_m[7]*div_b_comp[7]+0.3535533905932737*pkpm_div_ppar_c[6]*rho_inv[6]-0.3535533905932737*T_perp_over_m[6]*div_b_comp[6]+0.3535533905932737*pkpm_div_ppar_c[5]*rho_inv[5]-0.3535533905932737*T_perp_over_m[5]*div_b_comp[5]+0.3535533905932737*pkpm_div_ppar_c[4]*rho_inv[4]-0.3535533905932737*T_perp_over_m[4]*div_b_comp[4]+0.3535533905932737*pkpm_div_ppar_c[3]*rho_inv[3]-0.3535533905932737*T_perp_over_m[3]*div_b_comp[3]+0.3535533905932737*pkpm_div_ppar_c[2]*rho_inv[2]-0.3535533905932737*T_perp_over_m[2]*div_b_comp[2]+0.3535533905932737*pkpm_div_ppar_c[1]*rho_inv[1]-0.3535533905932737*T_perp_over_m[1]*div_b_comp[1]+0.3535533905932737*pkpm_div_ppar_c[0]*rho_inv[0]-0.3535533905932737*T_perp_over_m[0]*div_b_comp[0]; 
  p_force[1] += 0.3535533905932737*(pkpm_div_ppar_c[6]*rho_inv[7]+rho_inv[6]*pkpm_div_ppar_c[7])-0.3535533905932737*(T_perp_over_m[6]*div_b_comp[7]+div_b_comp[6]*T_perp_over_m[7])+0.3535533905932737*(pkpm_div_ppar_c[3]*rho_inv[5]+rho_inv[3]*pkpm_div_ppar_c[5])-0.3535533905932737*(T_perp_over_m[3]*div_b_comp[5]+div_b_comp[3]*T_perp_over_m[5])+0.3535533905932737*(pkpm_div_ppar_c[2]*rho_inv[4]+rho_inv[2]*pkpm_div_ppar_c[4])-0.3535533905932737*(T_perp_over_m[2]*div_b_comp[4]+div_b_comp[2]*T_perp_over_m[4])+0.3535533905932737*(pkpm_div_ppar_c[0]*rho_inv[1]+rho_inv[0]*pkpm_div_ppar_c[1])-0.3535533905932737*(T_perp_over_m[0]*div_b_comp[1]+div_b_comp[0]*T_perp_over_m[1]); 
  p_force[2] += 0.3535533905932737*(pkpm_div_ppar_c[5]*rho_inv[7]+rho_inv[5]*pkpm_div_ppar_c[7])-0.3535533905932737*(T_perp_over_m[5]*div_b_comp[7]+div_b_comp[5]*T_perp_over_m[7])+0.3535533905932737*(pkpm_div_ppar_c[3]*rho_inv[6]+rho_inv[3]*pkpm_div_ppar_c[6])-0.3535533905932737*(T_perp_over_m[3]*div_b_comp[6]+div_b_comp[3]*T_perp_over_m[6])+0.3535533905932737*(pkpm_div_ppar_c[1]*rho_inv[4]+rho_inv[1]*pkpm_div_ppar_c[4])-0.3535533905932737*(T_perp_over_m[1]*div_b_comp[4]+div_b_comp[1]*T_perp_over_m[4])+0.3535533905932737*(pkpm_div_ppar_c[0]*rho_inv[2]+rho_inv[0]*pkpm_div_ppar_c[2])-0.3535533905932737*(T_perp_over_m[0]*div_b_comp[2]+div_b_comp[0]*T_perp_over_m[2]); 
  p_force[3] += 0.3535533905932737*(pkpm_div_ppar_c[4]*rho_inv[7]+rho_inv[4]*pkpm_div_ppar_c[7])-0.3535533905932737*(T_perp_over_m[4]*div_b_comp[7]+div_b_comp[4]*T_perp_over_m[7])+0.3535533905932737*(pkpm_div_ppar_c[2]*rho_inv[6]+rho_inv[2]*pkpm_div_ppar_c[6])-0.3535533905932737*(T_perp_over_m[2]*div_b_comp[6]+div_b_comp[2]*T_perp_over_m[6])+0.3535533905932737*(pkpm_div_ppar_c[1]*rho_inv[5]+rho_inv[1]*pkpm_div_ppar_c[5])-0.3535533905932737*(T_perp_over_m[1]*div_b_comp[5]+div_b_comp[1]*T_perp_over_m[5])+0.3535533905932737*(pkpm_div_ppar_c[0]*rho_inv[3]+rho_inv[0]*pkpm_div_ppar_c[3])-0.3535533905932737*(T_perp_over_m[0]*div_b_comp[3]+div_b_comp[0]*T_perp_over_m[3]); 
  p_force[4] += 0.3535533905932737*(pkpm_div_ppar_c[3]*rho_inv[7]+rho_inv[3]*pkpm_div_ppar_c[7])-0.3535533905932737*(T_perp_over_m[3]*div_b_comp[7]+div_b_comp[3]*T_perp_over_m[7])+0.3535533905932737*(pkpm_div_ppar_c[5]*rho_inv[6]+rho_inv[5]*pkpm_div_ppar_c[6])-0.3535533905932737*(T_perp_over_m[5]*div_b_comp[6]+div_b_comp[5]*T_perp_over_m[6])+0.3535533905932737*(pkpm_div_ppar_c[0]*rho_inv[4]+rho_inv[0]*pkpm_div_ppar_c[4])-0.3535533905932737*(T_perp_over_m[0]*div_b_comp[4]+div_b_comp[0]*T_perp_over_m[4])+0.3535533905932737*(pkpm_div_ppar_c[1]*rho_inv[2]+rho_inv[1]*pkpm_div_ppar_c[2])-0.3535533905932737*(T_perp_over_m[1]*div_b_comp[2]+div_b_comp[1]*T_perp_over_m[2]); 
  p_force[5] += 0.3535533905932737*(pkpm_div_ppar_c[2]*rho_inv[7]+rho_inv[2]*pkpm_div_ppar_c[7])-0.3535533905932737*(T_perp_over_m[2]*div_b_comp[7]+div_b_comp[2]*T_perp_over_m[7])+0.3535533905932737*(pkpm_div_ppar_c[4]*rho_inv[6]+rho_inv[4]*pkpm_div_ppar_c[6])-0.3535533905932737*(T_perp_over_m[4]*div_b_comp[6]+div_b_comp[4]*T_perp_over_m[6])+0.3535533905932737*(pkpm_div_ppar_c[0]*rho_inv[5]+rho_inv[0]*pkpm_div_ppar_c[5])-0.3535533905932737*(T_perp_over_m[0]*div_b_comp[5]+div_b_comp[0]*T_perp_over_m[5])+0.3535533905932737*(pkpm_div_ppar_c[1]*rho_inv[3]+rho_inv[1]*pkpm_div_ppar_c[3])-0.3535533905932737*(T_perp_over_m[1]*div_b_comp[3]+div_b_comp[1]*T_perp_over_m[3]); 
  p_force[6] += 0.3535533905932737*(pkpm_div_ppar_c[1]*rho_inv[7]+rho_inv[1]*pkpm_div_ppar_c[7])-0.3535533905932737*(T_perp_over_m[1]*div_b_comp[7]+div_b_comp[1]*T_perp_over_m[7])+0.3535533905932737*(pkpm_div_ppar_c[0]*rho_inv[6]+rho_inv[0]*pkpm_div_ppar_c[6])-0.3535533905932737*(T_perp_over_m[0]*div_b_comp[6]+div_b_comp[0]*T_perp_over_m[6])+0.3535533905932737*(pkpm_div_ppar_c[4]*rho_inv[5]+rho_inv[4]*pkpm_div_ppar_c[5])-0.3535533905932737*(T_perp_over_m[4]*div_b_comp[5]+div_b_comp[4]*T_perp_over_m[5])+0.3535533905932737*(pkpm_div_ppar_c[2]*rho_inv[3]+rho_inv[2]*pkpm_div_ppar_c[3])-0.3535533905932737*(T_perp_over_m[2]*div_b_comp[3]+div_b_comp[2]*T_perp_over_m[3]); 
  p_force[7] += 0.3535533905932737*(pkpm_div_ppar_c[0]*rho_inv[7]+rho_inv[0]*pkpm_div_ppar_c[7])-0.3535533905932737*(T_perp_over_m[0]*div_b_comp[7]+div_b_comp[0]*T_perp_over_m[7])+0.3535533905932737*(pkpm_div_ppar_c[1]*rho_inv[6]+rho_inv[1]*pkpm_div_ppar_c[6])-0.3535533905932737*(T_perp_over_m[1]*div_b_comp[6]+div_b_comp[1]*T_perp_over_m[6])+0.3535533905932737*(pkpm_div_ppar_c[2]*rho_inv[5]+rho_inv[2]*pkpm_div_ppar_c[5])-0.3535533905932737*(T_perp_over_m[2]*div_b_comp[5]+div_b_comp[2]*T_perp_over_m[5])+0.3535533905932737*(pkpm_div_ppar_c[3]*rho_inv[4]+rho_inv[3]*pkpm_div_ppar_c[4])-0.3535533905932737*(T_perp_over_m[3]*div_b_comp[4]+div_b_comp[3]*T_perp_over_m[4]); 

  p_perp_source[0] += (-0.6666666666666666*nu[0])-1.0*grad_u_y[0]+bb_grad_u_comp[0]; 
  p_perp_source[1] += (-0.6666666666666666*nu[1])-1.0*grad_u_y[1]+bb_grad_u_comp[1]; 
  p_perp_source[2] += (-0.6666666666666666*nu[2])-1.0*grad_u_y[2]+bb_grad_u_comp[2]; 
  p_perp_source[3] += (-0.6666666666666666*nu[3])-1.0*grad_u_y[3]+bb_grad_u_comp[3]; 
  p_perp_source[4] += (-0.6666666666666666*nu[4])-1.0*grad_u_y[4]+bb_grad_u_comp[4]; 
  p_perp_source[5] += (-0.6666666666666666*nu[5])-1.0*grad_u_y[5]+bb_grad_u_comp[5]; 
  p_perp_source[6] += (-0.6666666666666666*nu[6])-1.0*grad_u_y[6]+bb_grad_u_comp[6]; 
  p_perp_source[7] += (-0.6666666666666666*nu[7])-1.0*grad_u_y[7]+bb_grad_u_comp[7]; 

  p_perp_div_b[0] += 0.3535533905932737*(T_perp_over_m[7]*div_b_comp[7]+T_perp_over_m[6]*div_b_comp[6]+T_perp_over_m[5]*div_b_comp[5]+T_perp_over_m[4]*div_b_comp[4]+T_perp_over_m[3]*div_b_comp[3]+T_perp_over_m[2]*div_b_comp[2]+T_perp_over_m[1]*div_b_comp[1]+T_perp_over_m[0]*div_b_comp[0]); 
  p_perp_div_b[1] += 0.3535533905932737*(T_perp_over_m[6]*div_b_comp[7]+div_b_comp[6]*T_perp_over_m[7]+T_perp_over_m[3]*div_b_comp[5]+div_b_comp[3]*T_perp_over_m[5]+T_perp_over_m[2]*div_b_comp[4]+div_b_comp[2]*T_perp_over_m[4]+T_perp_over_m[0]*div_b_comp[1]+div_b_comp[0]*T_perp_over_m[1]); 
  p_perp_div_b[2] += 0.3535533905932737*(T_perp_over_m[5]*div_b_comp[7]+div_b_comp[5]*T_perp_over_m[7]+T_perp_over_m[3]*div_b_comp[6]+div_b_comp[3]*T_perp_over_m[6]+T_perp_over_m[1]*div_b_comp[4]+div_b_comp[1]*T_perp_over_m[4]+T_perp_over_m[0]*div_b_comp[2]+div_b_comp[0]*T_perp_over_m[2]); 
  p_perp_div_b[3] += 0.3535533905932737*(T_perp_over_m[4]*div_b_comp[7]+div_b_comp[4]*T_perp_over_m[7]+T_perp_over_m[2]*div_b_comp[6]+div_b_comp[2]*T_perp_over_m[6]+T_perp_over_m[1]*div_b_comp[5]+div_b_comp[1]*T_perp_over_m[5]+T_perp_over_m[0]*div_b_comp[3]+div_b_comp[0]*T_perp_over_m[3]); 
  p_perp_div_b[4] += 0.3535533905932737*(T_perp_over_m[3]*div_b_comp[7]+div_b_comp[3]*T_perp_over_m[7]+T_perp_over_m[5]*div_b_comp[6]+div_b_comp[5]*T_perp_over_m[6]+T_perp_over_m[0]*div_b_comp[4]+div_b_comp[0]*T_perp_over_m[4]+T_perp_over_m[1]*div_b_comp[2]+div_b_comp[1]*T_perp_over_m[2]); 
  p_perp_div_b[5] += 0.3535533905932737*(T_perp_over_m[2]*div_b_comp[7]+div_b_comp[2]*T_perp_over_m[7]+T_perp_over_m[4]*div_b_comp[6]+div_b_comp[4]*T_perp_over_m[6]+T_perp_over_m[0]*div_b_comp[5]+div_b_comp[0]*T_perp_over_m[5]+T_perp_over_m[1]*div_b_comp[3]+div_b_comp[1]*T_perp_over_m[3]); 
  p_perp_div_b[6] += 0.3535533905932737*(T_perp_over_m[1]*div_b_comp[7]+div_b_comp[1]*T_perp_over_m[7]+T_perp_over_m[0]*div_b_comp[6]+div_b_comp[0]*T_perp_over_m[6]+T_perp_over_m[4]*div_b_comp[5]+div_b_comp[4]*T_perp_over_m[5]+T_perp_over_m[2]*div_b_comp[3]+div_b_comp[2]*T_perp_over_m[3]); 
  p_perp_div_b[7] += 0.3535533905932737*(T_perp_over_m[0]*div_b_comp[7]+div_b_comp[0]*T_perp_over_m[7]+T_perp_over_m[1]*div_b_comp[6]+div_b_comp[1]*T_perp_over_m[6]+T_perp_over_m[2]*div_b_comp[5]+div_b_comp[2]*T_perp_over_m[5]+T_perp_over_m[3]*div_b_comp[4]+div_b_comp[3]*T_perp_over_m[4]); 

} 
