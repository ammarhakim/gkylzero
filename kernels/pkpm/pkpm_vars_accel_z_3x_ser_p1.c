#include <gkyl_euler_pkpm_kernels.h> 
GKYL_CU_DH void pkpm_vars_accel_z_3x_ser_p1(const double *dxv, 
  const double *bvarl, const double *bvarc, const double *bvarr, 
  const double *priml, const double *primc, const double *primr, 
  const double *nu, double* GKYL_RESTRICT pkpm_accel) 
{ 
  // dxv[NDIM]:         Cell spacing.
  // bvarl/c/r:         Input magnetic field unit vector in left/center/right cells.
  // priml/primc/primr: Input [ux, uy, uz, 3*Txx/m, 3*Tyy/m, 3*Tzz/m, 1/rho div(p_par b), T_perp/m, m/T_perp] in left/center/right cells.
  // nu:                Input collisionality in center cell.
  // pkpm_accel:        Volume expansion of pkpm acceleration variables.

  const double dx1 = 2.0/dxv[2]; 
  const double *b_l = &bvarl[16]; 
  const double *b_c = &bvarc[16]; 
  const double *b_r = &bvarr[16]; 

  const double *ux_l = &priml[0]; 
  const double *uy_l = &priml[8]; 
  const double *uz_l = &priml[16]; 

  const double *ux_c = &primc[0]; 
  const double *uy_c = &primc[8]; 
  const double *uz_c = &primc[16]; 

  const double *ux_r = &primr[0]; 
  const double *uy_r = &primr[8]; 
  const double *uz_r = &primr[16]; 

  const double *bxbx = &bvarc[24]; 
  const double *bxby = &bvarc[32]; 
  const double *bxbz = &bvarc[40]; 
  const double *byby = &bvarc[48]; 
  const double *bybz = &bvarc[56]; 
  const double *bzbz = &bvarc[64]; 

  const double *pkpm_div_ppar = &primc[48]; 
  const double *T_perp_over_m = &primc[56]; 

  double *div_b = &pkpm_accel[0]; 
  double *bb_grad_u = &pkpm_accel[8]; 
  double *p_force = &pkpm_accel[16]; 
  double *p_perp_source = &pkpm_accel[24]; 
  double *p_perp_div_b = &pkpm_accel[32]; 

  double grad_u_x[8] = {0.0}; 
  double grad_u_y[8] = {0.0}; 
  double grad_u_z[8] = {0.0}; 
  grad_u_x[0] = (-0.4330127018922193*ux_r[3]*dx1)-0.4330127018922193*ux_l[3]*dx1+0.8660254037844386*ux_c[3]*dx1+0.25*ux_r[0]*dx1-0.25*ux_l[0]*dx1; 
  grad_u_x[1] = (-0.4330127018922193*ux_r[5]*dx1)-0.4330127018922193*ux_l[5]*dx1+0.8660254037844386*ux_c[5]*dx1+0.25*ux_r[1]*dx1-0.25*ux_l[1]*dx1; 
  grad_u_x[2] = (-0.4330127018922193*ux_r[6]*dx1)-0.4330127018922193*ux_l[6]*dx1+0.8660254037844386*ux_c[6]*dx1+0.25*ux_r[2]*dx1-0.25*ux_l[2]*dx1; 
  grad_u_x[3] = (-0.75*ux_r[3]*dx1)+0.75*ux_l[3]*dx1+0.4330127018922193*ux_r[0]*dx1+0.4330127018922193*ux_l[0]*dx1-0.8660254037844385*ux_c[0]*dx1; 
  grad_u_x[4] = (-0.4330127018922193*ux_r[7]*dx1)-0.4330127018922193*ux_l[7]*dx1+0.8660254037844386*ux_c[7]*dx1+0.25*ux_r[4]*dx1-0.25*ux_l[4]*dx1; 
  grad_u_x[5] = (-0.75*ux_r[5]*dx1)+0.75*ux_l[5]*dx1+0.4330127018922193*ux_r[1]*dx1+0.4330127018922193*ux_l[1]*dx1-0.8660254037844385*ux_c[1]*dx1; 
  grad_u_x[6] = (-0.75*ux_r[6]*dx1)+0.75*ux_l[6]*dx1+0.4330127018922193*ux_r[2]*dx1+0.4330127018922193*ux_l[2]*dx1-0.8660254037844385*ux_c[2]*dx1; 
  grad_u_x[7] = (-0.75*ux_r[7]*dx1)+0.75*ux_l[7]*dx1+0.4330127018922193*ux_r[4]*dx1+0.4330127018922193*ux_l[4]*dx1-0.8660254037844385*ux_c[4]*dx1; 

  grad_u_y[0] = (-0.4330127018922193*uy_r[3]*dx1)-0.4330127018922193*uy_l[3]*dx1+0.8660254037844386*uy_c[3]*dx1+0.25*uy_r[0]*dx1-0.25*uy_l[0]*dx1; 
  grad_u_y[1] = (-0.4330127018922193*uy_r[5]*dx1)-0.4330127018922193*uy_l[5]*dx1+0.8660254037844386*uy_c[5]*dx1+0.25*uy_r[1]*dx1-0.25*uy_l[1]*dx1; 
  grad_u_y[2] = (-0.4330127018922193*uy_r[6]*dx1)-0.4330127018922193*uy_l[6]*dx1+0.8660254037844386*uy_c[6]*dx1+0.25*uy_r[2]*dx1-0.25*uy_l[2]*dx1; 
  grad_u_y[3] = (-0.75*uy_r[3]*dx1)+0.75*uy_l[3]*dx1+0.4330127018922193*uy_r[0]*dx1+0.4330127018922193*uy_l[0]*dx1-0.8660254037844385*uy_c[0]*dx1; 
  grad_u_y[4] = (-0.4330127018922193*uy_r[7]*dx1)-0.4330127018922193*uy_l[7]*dx1+0.8660254037844386*uy_c[7]*dx1+0.25*uy_r[4]*dx1-0.25*uy_l[4]*dx1; 
  grad_u_y[5] = (-0.75*uy_r[5]*dx1)+0.75*uy_l[5]*dx1+0.4330127018922193*uy_r[1]*dx1+0.4330127018922193*uy_l[1]*dx1-0.8660254037844385*uy_c[1]*dx1; 
  grad_u_y[6] = (-0.75*uy_r[6]*dx1)+0.75*uy_l[6]*dx1+0.4330127018922193*uy_r[2]*dx1+0.4330127018922193*uy_l[2]*dx1-0.8660254037844385*uy_c[2]*dx1; 
  grad_u_y[7] = (-0.75*uy_r[7]*dx1)+0.75*uy_l[7]*dx1+0.4330127018922193*uy_r[4]*dx1+0.4330127018922193*uy_l[4]*dx1-0.8660254037844385*uy_c[4]*dx1; 

  grad_u_z[0] = (-0.4330127018922193*uz_r[3]*dx1)-0.4330127018922193*uz_l[3]*dx1+0.8660254037844386*uz_c[3]*dx1+0.25*uz_r[0]*dx1-0.25*uz_l[0]*dx1; 
  grad_u_z[1] = (-0.4330127018922193*uz_r[5]*dx1)-0.4330127018922193*uz_l[5]*dx1+0.8660254037844386*uz_c[5]*dx1+0.25*uz_r[1]*dx1-0.25*uz_l[1]*dx1; 
  grad_u_z[2] = (-0.4330127018922193*uz_r[6]*dx1)-0.4330127018922193*uz_l[6]*dx1+0.8660254037844386*uz_c[6]*dx1+0.25*uz_r[2]*dx1-0.25*uz_l[2]*dx1; 
  grad_u_z[3] = (-0.75*uz_r[3]*dx1)+0.75*uz_l[3]*dx1+0.4330127018922193*uz_r[0]*dx1+0.4330127018922193*uz_l[0]*dx1-0.8660254037844385*uz_c[0]*dx1; 
  grad_u_z[4] = (-0.4330127018922193*uz_r[7]*dx1)-0.4330127018922193*uz_l[7]*dx1+0.8660254037844386*uz_c[7]*dx1+0.25*uz_r[4]*dx1-0.25*uz_l[4]*dx1; 
  grad_u_z[5] = (-0.75*uz_r[5]*dx1)+0.75*uz_l[5]*dx1+0.4330127018922193*uz_r[1]*dx1+0.4330127018922193*uz_l[1]*dx1-0.8660254037844385*uz_c[1]*dx1; 
  grad_u_z[6] = (-0.75*uz_r[6]*dx1)+0.75*uz_l[6]*dx1+0.4330127018922193*uz_r[2]*dx1+0.4330127018922193*uz_l[2]*dx1-0.8660254037844385*uz_c[2]*dx1; 
  grad_u_z[7] = (-0.75*uz_r[7]*dx1)+0.75*uz_l[7]*dx1+0.4330127018922193*uz_r[4]*dx1+0.4330127018922193*uz_l[4]*dx1-0.8660254037844385*uz_c[4]*dx1; 

  double div_b_comp[8] = {0.0}; 
  double div_p_x_comp[8] = {0.0}; 
  double div_p_y_comp[8] = {0.0}; 
  double div_p_z_comp[8] = {0.0}; 
  double bb_grad_u_comp[8] = {0.0}; 
  div_b_comp[0] = (-0.2886751345948129*b_r[3]*dx1)-0.2886751345948129*b_l[3]*dx1+0.5773502691896258*b_c[3]*dx1+0.25*b_r[0]*dx1-0.25*b_l[0]*dx1; 
  div_b_comp[1] = (-0.2886751345948129*b_r[5]*dx1)-0.2886751345948129*b_l[5]*dx1+0.5773502691896258*b_c[5]*dx1+0.25*b_r[1]*dx1-0.25*b_l[1]*dx1; 
  div_b_comp[2] = (-0.2886751345948129*b_r[6]*dx1)-0.2886751345948129*b_l[6]*dx1+0.5773502691896258*b_c[6]*dx1+0.25*b_r[2]*dx1-0.25*b_l[2]*dx1; 
  div_b_comp[3] = (-0.5*b_r[3]*dx1)+0.5*b_l[3]*dx1+0.4330127018922193*b_r[0]*dx1+0.4330127018922193*b_l[0]*dx1-0.8660254037844386*b_c[0]*dx1; 
  div_b_comp[4] = (-0.2886751345948129*b_r[7]*dx1)-0.2886751345948129*b_l[7]*dx1+0.5773502691896258*b_c[7]*dx1+0.25*b_r[4]*dx1-0.25*b_l[4]*dx1; 
  div_b_comp[5] = (-0.5*b_r[5]*dx1)+0.5*b_l[5]*dx1+0.4330127018922193*b_r[1]*dx1+0.4330127018922193*b_l[1]*dx1-0.8660254037844386*b_c[1]*dx1; 
  div_b_comp[6] = (-0.5*b_r[6]*dx1)+0.5*b_l[6]*dx1+0.4330127018922193*b_r[2]*dx1+0.4330127018922193*b_l[2]*dx1-0.8660254037844386*b_c[2]*dx1; 
  div_b_comp[7] = (-0.5*b_r[7]*dx1)+0.5*b_l[7]*dx1+0.4330127018922193*b_r[4]*dx1+0.4330127018922193*b_l[4]*dx1-0.8660254037844386*b_c[4]*dx1; 

  div_b[0] += div_b_comp[0]; 
  div_b[1] += div_b_comp[1]; 
  div_b[2] += div_b_comp[2]; 
  div_b[3] += div_b_comp[3]; 
  div_b[4] += div_b_comp[4]; 
  div_b[5] += div_b_comp[5]; 
  div_b[6] += div_b_comp[6]; 
  div_b[7] += div_b_comp[7]; 

  bb_grad_u_comp[0] = 0.3535533905932737*bzbz[7]*grad_u_z[7]+0.3535533905932737*bybz[7]*grad_u_y[7]+0.3535533905932737*bxbz[7]*grad_u_x[7]+0.3535533905932737*bzbz[6]*grad_u_z[6]+0.3535533905932737*bybz[6]*grad_u_y[6]+0.3535533905932737*bxbz[6]*grad_u_x[6]+0.3535533905932737*bzbz[5]*grad_u_z[5]+0.3535533905932737*bybz[5]*grad_u_y[5]+0.3535533905932737*bxbz[5]*grad_u_x[5]+0.3535533905932737*bzbz[4]*grad_u_z[4]+0.3535533905932737*bybz[4]*grad_u_y[4]+0.3535533905932737*bxbz[4]*grad_u_x[4]+0.3535533905932737*bzbz[3]*grad_u_z[3]+0.3535533905932737*bybz[3]*grad_u_y[3]+0.3535533905932737*bxbz[3]*grad_u_x[3]+0.3535533905932737*bzbz[2]*grad_u_z[2]+0.3535533905932737*bybz[2]*grad_u_y[2]+0.3535533905932737*bxbz[2]*grad_u_x[2]+0.3535533905932737*bzbz[1]*grad_u_z[1]+0.3535533905932737*bybz[1]*grad_u_y[1]+0.3535533905932737*bxbz[1]*grad_u_x[1]+0.3535533905932737*bzbz[0]*grad_u_z[0]+0.3535533905932737*bybz[0]*grad_u_y[0]+0.3535533905932737*bxbz[0]*grad_u_x[0]; 
  bb_grad_u_comp[1] = 0.3535533905932737*bzbz[6]*grad_u_z[7]+0.3535533905932737*bybz[6]*grad_u_y[7]+0.3535533905932737*bxbz[6]*grad_u_x[7]+0.3535533905932737*grad_u_z[6]*bzbz[7]+0.3535533905932737*grad_u_y[6]*bybz[7]+0.3535533905932737*grad_u_x[6]*bxbz[7]+0.3535533905932737*bzbz[3]*grad_u_z[5]+0.3535533905932737*bybz[3]*grad_u_y[5]+0.3535533905932737*bxbz[3]*grad_u_x[5]+0.3535533905932737*grad_u_z[3]*bzbz[5]+0.3535533905932737*grad_u_y[3]*bybz[5]+0.3535533905932737*grad_u_x[3]*bxbz[5]+0.3535533905932737*bzbz[2]*grad_u_z[4]+0.3535533905932737*bybz[2]*grad_u_y[4]+0.3535533905932737*bxbz[2]*grad_u_x[4]+0.3535533905932737*grad_u_z[2]*bzbz[4]+0.3535533905932737*grad_u_y[2]*bybz[4]+0.3535533905932737*grad_u_x[2]*bxbz[4]+0.3535533905932737*bzbz[0]*grad_u_z[1]+0.3535533905932737*bybz[0]*grad_u_y[1]+0.3535533905932737*bxbz[0]*grad_u_x[1]+0.3535533905932737*grad_u_z[0]*bzbz[1]+0.3535533905932737*grad_u_y[0]*bybz[1]+0.3535533905932737*grad_u_x[0]*bxbz[1]; 
  bb_grad_u_comp[2] = 0.3535533905932737*bzbz[5]*grad_u_z[7]+0.3535533905932737*bybz[5]*grad_u_y[7]+0.3535533905932737*bxbz[5]*grad_u_x[7]+0.3535533905932737*grad_u_z[5]*bzbz[7]+0.3535533905932737*grad_u_y[5]*bybz[7]+0.3535533905932737*grad_u_x[5]*bxbz[7]+0.3535533905932737*bzbz[3]*grad_u_z[6]+0.3535533905932737*bybz[3]*grad_u_y[6]+0.3535533905932737*bxbz[3]*grad_u_x[6]+0.3535533905932737*grad_u_z[3]*bzbz[6]+0.3535533905932737*grad_u_y[3]*bybz[6]+0.3535533905932737*grad_u_x[3]*bxbz[6]+0.3535533905932737*bzbz[1]*grad_u_z[4]+0.3535533905932737*bybz[1]*grad_u_y[4]+0.3535533905932737*bxbz[1]*grad_u_x[4]+0.3535533905932737*grad_u_z[1]*bzbz[4]+0.3535533905932737*grad_u_y[1]*bybz[4]+0.3535533905932737*grad_u_x[1]*bxbz[4]+0.3535533905932737*bzbz[0]*grad_u_z[2]+0.3535533905932737*bybz[0]*grad_u_y[2]+0.3535533905932737*bxbz[0]*grad_u_x[2]+0.3535533905932737*grad_u_z[0]*bzbz[2]+0.3535533905932737*grad_u_y[0]*bybz[2]+0.3535533905932737*grad_u_x[0]*bxbz[2]; 
  bb_grad_u_comp[3] = 0.3535533905932737*bzbz[4]*grad_u_z[7]+0.3535533905932737*bybz[4]*grad_u_y[7]+0.3535533905932737*bxbz[4]*grad_u_x[7]+0.3535533905932737*grad_u_z[4]*bzbz[7]+0.3535533905932737*grad_u_y[4]*bybz[7]+0.3535533905932737*grad_u_x[4]*bxbz[7]+0.3535533905932737*bzbz[2]*grad_u_z[6]+0.3535533905932737*bybz[2]*grad_u_y[6]+0.3535533905932737*bxbz[2]*grad_u_x[6]+0.3535533905932737*grad_u_z[2]*bzbz[6]+0.3535533905932737*grad_u_y[2]*bybz[6]+0.3535533905932737*grad_u_x[2]*bxbz[6]+0.3535533905932737*bzbz[1]*grad_u_z[5]+0.3535533905932737*bybz[1]*grad_u_y[5]+0.3535533905932737*bxbz[1]*grad_u_x[5]+0.3535533905932737*grad_u_z[1]*bzbz[5]+0.3535533905932737*grad_u_y[1]*bybz[5]+0.3535533905932737*grad_u_x[1]*bxbz[5]+0.3535533905932737*bzbz[0]*grad_u_z[3]+0.3535533905932737*bybz[0]*grad_u_y[3]+0.3535533905932737*bxbz[0]*grad_u_x[3]+0.3535533905932737*grad_u_z[0]*bzbz[3]+0.3535533905932737*grad_u_y[0]*bybz[3]+0.3535533905932737*grad_u_x[0]*bxbz[3]; 
  bb_grad_u_comp[4] = 0.3535533905932737*bzbz[3]*grad_u_z[7]+0.3535533905932737*bybz[3]*grad_u_y[7]+0.3535533905932737*bxbz[3]*grad_u_x[7]+0.3535533905932737*grad_u_z[3]*bzbz[7]+0.3535533905932737*grad_u_y[3]*bybz[7]+0.3535533905932737*grad_u_x[3]*bxbz[7]+0.3535533905932737*bzbz[5]*grad_u_z[6]+0.3535533905932737*bybz[5]*grad_u_y[6]+0.3535533905932737*bxbz[5]*grad_u_x[6]+0.3535533905932737*grad_u_z[5]*bzbz[6]+0.3535533905932737*grad_u_y[5]*bybz[6]+0.3535533905932737*grad_u_x[5]*bxbz[6]+0.3535533905932737*bzbz[0]*grad_u_z[4]+0.3535533905932737*bybz[0]*grad_u_y[4]+0.3535533905932737*bxbz[0]*grad_u_x[4]+0.3535533905932737*grad_u_z[0]*bzbz[4]+0.3535533905932737*grad_u_y[0]*bybz[4]+0.3535533905932737*grad_u_x[0]*bxbz[4]+0.3535533905932737*bzbz[1]*grad_u_z[2]+0.3535533905932737*bybz[1]*grad_u_y[2]+0.3535533905932737*bxbz[1]*grad_u_x[2]+0.3535533905932737*grad_u_z[1]*bzbz[2]+0.3535533905932737*grad_u_y[1]*bybz[2]+0.3535533905932737*grad_u_x[1]*bxbz[2]; 
  bb_grad_u_comp[5] = 0.3535533905932737*bzbz[2]*grad_u_z[7]+0.3535533905932737*bybz[2]*grad_u_y[7]+0.3535533905932737*bxbz[2]*grad_u_x[7]+0.3535533905932737*grad_u_z[2]*bzbz[7]+0.3535533905932737*grad_u_y[2]*bybz[7]+0.3535533905932737*grad_u_x[2]*bxbz[7]+0.3535533905932737*bzbz[4]*grad_u_z[6]+0.3535533905932737*bybz[4]*grad_u_y[6]+0.3535533905932737*bxbz[4]*grad_u_x[6]+0.3535533905932737*grad_u_z[4]*bzbz[6]+0.3535533905932737*grad_u_y[4]*bybz[6]+0.3535533905932737*grad_u_x[4]*bxbz[6]+0.3535533905932737*bzbz[0]*grad_u_z[5]+0.3535533905932737*bybz[0]*grad_u_y[5]+0.3535533905932737*bxbz[0]*grad_u_x[5]+0.3535533905932737*grad_u_z[0]*bzbz[5]+0.3535533905932737*grad_u_y[0]*bybz[5]+0.3535533905932737*grad_u_x[0]*bxbz[5]+0.3535533905932737*bzbz[1]*grad_u_z[3]+0.3535533905932737*bybz[1]*grad_u_y[3]+0.3535533905932737*bxbz[1]*grad_u_x[3]+0.3535533905932737*grad_u_z[1]*bzbz[3]+0.3535533905932737*grad_u_y[1]*bybz[3]+0.3535533905932737*grad_u_x[1]*bxbz[3]; 
  bb_grad_u_comp[6] = 0.3535533905932737*bzbz[1]*grad_u_z[7]+0.3535533905932737*bybz[1]*grad_u_y[7]+0.3535533905932737*bxbz[1]*grad_u_x[7]+0.3535533905932737*grad_u_z[1]*bzbz[7]+0.3535533905932737*grad_u_y[1]*bybz[7]+0.3535533905932737*grad_u_x[1]*bxbz[7]+0.3535533905932737*bzbz[0]*grad_u_z[6]+0.3535533905932737*bybz[0]*grad_u_y[6]+0.3535533905932737*bxbz[0]*grad_u_x[6]+0.3535533905932737*grad_u_z[0]*bzbz[6]+0.3535533905932737*grad_u_y[0]*bybz[6]+0.3535533905932737*grad_u_x[0]*bxbz[6]+0.3535533905932737*bzbz[4]*grad_u_z[5]+0.3535533905932737*bybz[4]*grad_u_y[5]+0.3535533905932737*bxbz[4]*grad_u_x[5]+0.3535533905932737*grad_u_z[4]*bzbz[5]+0.3535533905932737*grad_u_y[4]*bybz[5]+0.3535533905932737*grad_u_x[4]*bxbz[5]+0.3535533905932737*bzbz[2]*grad_u_z[3]+0.3535533905932737*bybz[2]*grad_u_y[3]+0.3535533905932737*bxbz[2]*grad_u_x[3]+0.3535533905932737*grad_u_z[2]*bzbz[3]+0.3535533905932737*grad_u_y[2]*bybz[3]+0.3535533905932737*grad_u_x[2]*bxbz[3]; 
  bb_grad_u_comp[7] = 0.3535533905932737*bzbz[0]*grad_u_z[7]+0.3535533905932737*bybz[0]*grad_u_y[7]+0.3535533905932737*bxbz[0]*grad_u_x[7]+0.3535533905932737*grad_u_z[0]*bzbz[7]+0.3535533905932737*grad_u_y[0]*bybz[7]+0.3535533905932737*grad_u_x[0]*bxbz[7]+0.3535533905932737*bzbz[1]*grad_u_z[6]+0.3535533905932737*bybz[1]*grad_u_y[6]+0.3535533905932737*bxbz[1]*grad_u_x[6]+0.3535533905932737*grad_u_z[1]*bzbz[6]+0.3535533905932737*grad_u_y[1]*bybz[6]+0.3535533905932737*grad_u_x[1]*bxbz[6]+0.3535533905932737*bzbz[2]*grad_u_z[5]+0.3535533905932737*bybz[2]*grad_u_y[5]+0.3535533905932737*bxbz[2]*grad_u_x[5]+0.3535533905932737*grad_u_z[2]*bzbz[5]+0.3535533905932737*grad_u_y[2]*bybz[5]+0.3535533905932737*grad_u_x[2]*bxbz[5]+0.3535533905932737*bzbz[3]*grad_u_z[4]+0.3535533905932737*bybz[3]*grad_u_y[4]+0.3535533905932737*bxbz[3]*grad_u_x[4]+0.3535533905932737*grad_u_z[3]*bzbz[4]+0.3535533905932737*grad_u_y[3]*bybz[4]+0.3535533905932737*grad_u_x[3]*bxbz[4]; 

  bb_grad_u[0] += bb_grad_u_comp[0]; 
  bb_grad_u[1] += bb_grad_u_comp[1]; 
  bb_grad_u[2] += bb_grad_u_comp[2]; 
  bb_grad_u[3] += bb_grad_u_comp[3]; 
  bb_grad_u[4] += bb_grad_u_comp[4]; 
  bb_grad_u[5] += bb_grad_u_comp[5]; 
  bb_grad_u[6] += bb_grad_u_comp[6]; 
  bb_grad_u[7] += bb_grad_u_comp[7]; 

  p_force[0] += (-0.3535533905932737*(T_perp_over_m[7]*div_b_comp[7]+T_perp_over_m[6]*div_b_comp[6]+T_perp_over_m[5]*div_b_comp[5]+T_perp_over_m[4]*div_b_comp[4]+T_perp_over_m[3]*div_b_comp[3]+T_perp_over_m[2]*div_b_comp[2]+T_perp_over_m[1]*div_b_comp[1]))+0.3333333333333333*pkpm_div_ppar[0]-0.3535533905932737*T_perp_over_m[0]*div_b_comp[0]; 
  p_force[1] += (-0.3535533905932737*(T_perp_over_m[6]*div_b_comp[7]+div_b_comp[6]*T_perp_over_m[7]+T_perp_over_m[3]*div_b_comp[5]+div_b_comp[3]*T_perp_over_m[5]+T_perp_over_m[2]*div_b_comp[4]+div_b_comp[2]*T_perp_over_m[4]))+0.3333333333333333*pkpm_div_ppar[1]-0.3535533905932737*(T_perp_over_m[0]*div_b_comp[1]+div_b_comp[0]*T_perp_over_m[1]); 
  p_force[2] += (-0.3535533905932737*(T_perp_over_m[5]*div_b_comp[7]+div_b_comp[5]*T_perp_over_m[7]+T_perp_over_m[3]*div_b_comp[6]+div_b_comp[3]*T_perp_over_m[6]+T_perp_over_m[1]*div_b_comp[4]+div_b_comp[1]*T_perp_over_m[4]))+0.3333333333333333*pkpm_div_ppar[2]-0.3535533905932737*(T_perp_over_m[0]*div_b_comp[2]+div_b_comp[0]*T_perp_over_m[2]); 
  p_force[3] += (-0.3535533905932737*(T_perp_over_m[4]*div_b_comp[7]+div_b_comp[4]*T_perp_over_m[7]+T_perp_over_m[2]*div_b_comp[6]+div_b_comp[2]*T_perp_over_m[6]+T_perp_over_m[1]*div_b_comp[5]+div_b_comp[1]*T_perp_over_m[5]))+0.3333333333333333*pkpm_div_ppar[3]-0.3535533905932737*(T_perp_over_m[0]*div_b_comp[3]+div_b_comp[0]*T_perp_over_m[3]); 
  p_force[4] += (-0.3535533905932737*(T_perp_over_m[3]*div_b_comp[7]+div_b_comp[3]*T_perp_over_m[7]+T_perp_over_m[5]*div_b_comp[6]+div_b_comp[5]*T_perp_over_m[6]))+0.3333333333333333*pkpm_div_ppar[4]-0.3535533905932737*(T_perp_over_m[0]*div_b_comp[4]+div_b_comp[0]*T_perp_over_m[4]+T_perp_over_m[1]*div_b_comp[2]+div_b_comp[1]*T_perp_over_m[2]); 
  p_force[5] += (-0.3535533905932737*(T_perp_over_m[2]*div_b_comp[7]+div_b_comp[2]*T_perp_over_m[7]+T_perp_over_m[4]*div_b_comp[6]+div_b_comp[4]*T_perp_over_m[6]))+0.3333333333333333*pkpm_div_ppar[5]-0.3535533905932737*(T_perp_over_m[0]*div_b_comp[5]+div_b_comp[0]*T_perp_over_m[5]+T_perp_over_m[1]*div_b_comp[3]+div_b_comp[1]*T_perp_over_m[3]); 
  p_force[6] += (-0.3535533905932737*(T_perp_over_m[1]*div_b_comp[7]+div_b_comp[1]*T_perp_over_m[7]))+0.3333333333333333*pkpm_div_ppar[6]-0.3535533905932737*(T_perp_over_m[0]*div_b_comp[6]+div_b_comp[0]*T_perp_over_m[6]+T_perp_over_m[4]*div_b_comp[5]+div_b_comp[4]*T_perp_over_m[5]+T_perp_over_m[2]*div_b_comp[3]+div_b_comp[2]*T_perp_over_m[3]); 
  p_force[7] += 0.3333333333333333*pkpm_div_ppar[7]-0.3535533905932737*(T_perp_over_m[0]*div_b_comp[7]+div_b_comp[0]*T_perp_over_m[7]+T_perp_over_m[1]*div_b_comp[6]+div_b_comp[1]*T_perp_over_m[6]+T_perp_over_m[2]*div_b_comp[5]+div_b_comp[2]*T_perp_over_m[5]+T_perp_over_m[3]*div_b_comp[4]+div_b_comp[3]*T_perp_over_m[4]); 

  p_perp_source[0] += (-0.6666666666666666*nu[0])-1.0*grad_u_z[0]+bb_grad_u_comp[0]; 
  p_perp_source[1] += (-0.6666666666666666*nu[1])-1.0*grad_u_z[1]+bb_grad_u_comp[1]; 
  p_perp_source[2] += (-0.6666666666666666*nu[2])-1.0*grad_u_z[2]+bb_grad_u_comp[2]; 
  p_perp_source[3] += (-0.6666666666666666*nu[3])-1.0*grad_u_z[3]+bb_grad_u_comp[3]; 
  p_perp_source[4] += (-0.6666666666666666*nu[4])-1.0*grad_u_z[4]+bb_grad_u_comp[4]; 
  p_perp_source[5] += (-0.6666666666666666*nu[5])-1.0*grad_u_z[5]+bb_grad_u_comp[5]; 
  p_perp_source[6] += (-0.6666666666666666*nu[6])-1.0*grad_u_z[6]+bb_grad_u_comp[6]; 
  p_perp_source[7] += (-0.6666666666666666*nu[7])-1.0*grad_u_z[7]+bb_grad_u_comp[7]; 

  p_perp_div_b[0] += 0.3535533905932737*(T_perp_over_m[7]*div_b_comp[7]+T_perp_over_m[6]*div_b_comp[6]+T_perp_over_m[5]*div_b_comp[5]+T_perp_over_m[4]*div_b_comp[4]+T_perp_over_m[3]*div_b_comp[3]+T_perp_over_m[2]*div_b_comp[2]+T_perp_over_m[1]*div_b_comp[1]+T_perp_over_m[0]*div_b_comp[0]); 
  p_perp_div_b[1] += 0.3535533905932737*(T_perp_over_m[6]*div_b_comp[7]+div_b_comp[6]*T_perp_over_m[7]+T_perp_over_m[3]*div_b_comp[5]+div_b_comp[3]*T_perp_over_m[5]+T_perp_over_m[2]*div_b_comp[4]+div_b_comp[2]*T_perp_over_m[4]+T_perp_over_m[0]*div_b_comp[1]+div_b_comp[0]*T_perp_over_m[1]); 
  p_perp_div_b[2] += 0.3535533905932737*(T_perp_over_m[5]*div_b_comp[7]+div_b_comp[5]*T_perp_over_m[7]+T_perp_over_m[3]*div_b_comp[6]+div_b_comp[3]*T_perp_over_m[6]+T_perp_over_m[1]*div_b_comp[4]+div_b_comp[1]*T_perp_over_m[4]+T_perp_over_m[0]*div_b_comp[2]+div_b_comp[0]*T_perp_over_m[2]); 
  p_perp_div_b[3] += 0.3535533905932737*(T_perp_over_m[4]*div_b_comp[7]+div_b_comp[4]*T_perp_over_m[7]+T_perp_over_m[2]*div_b_comp[6]+div_b_comp[2]*T_perp_over_m[6]+T_perp_over_m[1]*div_b_comp[5]+div_b_comp[1]*T_perp_over_m[5]+T_perp_over_m[0]*div_b_comp[3]+div_b_comp[0]*T_perp_over_m[3]); 
  p_perp_div_b[4] += 0.3535533905932737*(T_perp_over_m[3]*div_b_comp[7]+div_b_comp[3]*T_perp_over_m[7]+T_perp_over_m[5]*div_b_comp[6]+div_b_comp[5]*T_perp_over_m[6]+T_perp_over_m[0]*div_b_comp[4]+div_b_comp[0]*T_perp_over_m[4]+T_perp_over_m[1]*div_b_comp[2]+div_b_comp[1]*T_perp_over_m[2]); 
  p_perp_div_b[5] += 0.3535533905932737*(T_perp_over_m[2]*div_b_comp[7]+div_b_comp[2]*T_perp_over_m[7]+T_perp_over_m[4]*div_b_comp[6]+div_b_comp[4]*T_perp_over_m[6]+T_perp_over_m[0]*div_b_comp[5]+div_b_comp[0]*T_perp_over_m[5]+T_perp_over_m[1]*div_b_comp[3]+div_b_comp[1]*T_perp_over_m[3]); 
  p_perp_div_b[6] += 0.3535533905932737*(T_perp_over_m[1]*div_b_comp[7]+div_b_comp[1]*T_perp_over_m[7]+T_perp_over_m[0]*div_b_comp[6]+div_b_comp[0]*T_perp_over_m[6]+T_perp_over_m[4]*div_b_comp[5]+div_b_comp[4]*T_perp_over_m[5]+T_perp_over_m[2]*div_b_comp[3]+div_b_comp[2]*T_perp_over_m[3]); 
  p_perp_div_b[7] += 0.3535533905932737*(T_perp_over_m[0]*div_b_comp[7]+div_b_comp[0]*T_perp_over_m[7]+T_perp_over_m[1]*div_b_comp[6]+div_b_comp[1]*T_perp_over_m[6]+T_perp_over_m[2]*div_b_comp[5]+div_b_comp[2]*T_perp_over_m[5]+T_perp_over_m[3]*div_b_comp[4]+div_b_comp[3]*T_perp_over_m[4]); 

} 
