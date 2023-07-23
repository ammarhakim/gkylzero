#include <gkyl_euler_pkpm_kernels.h> 
GKYL_CU_DH void pkpm_vars_accel_x_2x_ser_p1(const double *dxv, 
  const double *bvarl, const double *bvarc, const double *bvarr, 
  const double *priml, const double *primc, const double *primr, 
  const double *nu, double* GKYL_RESTRICT pkpm_accel) 
{ 
  // dxv[NDIM]:         Cell spacing.
  // bvarl/c/r:         Input magnetic field unit vector in left/center/right cells.
  // priml/primc/primr: Input [ux, uy, uz, 3*Txx/m, 3*Tyy/m, 3*Tzz/m, 1/rho div(p_par b), T_perp/m, m/T_perp] in left/center/right cells.
  // nu:                Input collisionality in center cell.
  // pkpm_accel:        Volume expansion of pkpm acceleration variables.

  const double dx1 = 2.0/dxv[0]; 
  const double *b_l = &bvarl[0]; 
  const double *b_c = &bvarc[0]; 
  const double *b_r = &bvarr[0]; 

  const double *ux_l = &priml[0]; 
  const double *uy_l = &priml[4]; 
  const double *uz_l = &priml[8]; 

  const double *ux_c = &primc[0]; 
  const double *uy_c = &primc[4]; 
  const double *uz_c = &primc[8]; 

  const double *ux_r = &primr[0]; 
  const double *uy_r = &primr[4]; 
  const double *uz_r = &primr[8]; 

  const double *bxbx = &bvarc[12]; 
  const double *bxby = &bvarc[16]; 
  const double *bxbz = &bvarc[20]; 
  const double *byby = &bvarc[24]; 
  const double *bybz = &bvarc[28]; 
  const double *bzbz = &bvarc[32]; 

  const double *pkpm_div_ppar = &primc[24]; 
  const double *T_perp_over_m = &primc[28]; 

  double *div_b = &pkpm_accel[0]; 
  double *bb_grad_u = &pkpm_accel[4]; 
  double *p_force = &pkpm_accel[8]; 
  double *p_perp_source = &pkpm_accel[12]; 
  double *p_perp_div_b = &pkpm_accel[16]; 

  double grad_u_x[4] = {0.0}; 
  double grad_u_y[4] = {0.0}; 
  double grad_u_z[4] = {0.0}; 
  grad_u_x[0] = (-0.4330127018922193*ux_r[1]*dx1)-0.4330127018922193*ux_l[1]*dx1+0.8660254037844386*ux_c[1]*dx1+0.25*ux_r[0]*dx1-0.25*ux_l[0]*dx1; 
  grad_u_x[1] = (-0.75*ux_r[1]*dx1)+0.75*ux_l[1]*dx1+0.4330127018922193*ux_r[0]*dx1+0.4330127018922193*ux_l[0]*dx1-0.8660254037844385*ux_c[0]*dx1; 
  grad_u_x[2] = (-0.4330127018922193*ux_r[3]*dx1)-0.4330127018922193*ux_l[3]*dx1+0.8660254037844386*ux_c[3]*dx1+0.25*ux_r[2]*dx1-0.25*ux_l[2]*dx1; 
  grad_u_x[3] = (-0.75*ux_r[3]*dx1)+0.75*ux_l[3]*dx1+0.4330127018922193*ux_r[2]*dx1+0.4330127018922193*ux_l[2]*dx1-0.8660254037844385*ux_c[2]*dx1; 

  grad_u_y[0] = (-0.4330127018922193*uy_r[1]*dx1)-0.4330127018922193*uy_l[1]*dx1+0.8660254037844386*uy_c[1]*dx1+0.25*uy_r[0]*dx1-0.25*uy_l[0]*dx1; 
  grad_u_y[1] = (-0.75*uy_r[1]*dx1)+0.75*uy_l[1]*dx1+0.4330127018922193*uy_r[0]*dx1+0.4330127018922193*uy_l[0]*dx1-0.8660254037844385*uy_c[0]*dx1; 
  grad_u_y[2] = (-0.4330127018922193*uy_r[3]*dx1)-0.4330127018922193*uy_l[3]*dx1+0.8660254037844386*uy_c[3]*dx1+0.25*uy_r[2]*dx1-0.25*uy_l[2]*dx1; 
  grad_u_y[3] = (-0.75*uy_r[3]*dx1)+0.75*uy_l[3]*dx1+0.4330127018922193*uy_r[2]*dx1+0.4330127018922193*uy_l[2]*dx1-0.8660254037844385*uy_c[2]*dx1; 

  grad_u_z[0] = (-0.4330127018922193*uz_r[1]*dx1)-0.4330127018922193*uz_l[1]*dx1+0.8660254037844386*uz_c[1]*dx1+0.25*uz_r[0]*dx1-0.25*uz_l[0]*dx1; 
  grad_u_z[1] = (-0.75*uz_r[1]*dx1)+0.75*uz_l[1]*dx1+0.4330127018922193*uz_r[0]*dx1+0.4330127018922193*uz_l[0]*dx1-0.8660254037844385*uz_c[0]*dx1; 
  grad_u_z[2] = (-0.4330127018922193*uz_r[3]*dx1)-0.4330127018922193*uz_l[3]*dx1+0.8660254037844386*uz_c[3]*dx1+0.25*uz_r[2]*dx1-0.25*uz_l[2]*dx1; 
  grad_u_z[3] = (-0.75*uz_r[3]*dx1)+0.75*uz_l[3]*dx1+0.4330127018922193*uz_r[2]*dx1+0.4330127018922193*uz_l[2]*dx1-0.8660254037844385*uz_c[2]*dx1; 

  double div_b_comp[4] = {0.0}; 
  double div_p_x_comp[4] = {0.0}; 
  double div_p_y_comp[4] = {0.0}; 
  double div_p_z_comp[4] = {0.0}; 
  double bb_grad_u_comp[4] = {0.0}; 
  div_b_comp[0] = (-0.2886751345948129*b_r[1]*dx1)-0.2886751345948129*b_l[1]*dx1+0.5773502691896258*b_c[1]*dx1+0.25*b_r[0]*dx1-0.25*b_l[0]*dx1; 
  div_b_comp[1] = (-0.5*b_r[1]*dx1)+0.5*b_l[1]*dx1+0.4330127018922193*b_r[0]*dx1+0.4330127018922193*b_l[0]*dx1-0.8660254037844386*b_c[0]*dx1; 
  div_b_comp[2] = (-0.2886751345948129*b_r[3]*dx1)-0.2886751345948129*b_l[3]*dx1+0.5773502691896258*b_c[3]*dx1+0.25*b_r[2]*dx1-0.25*b_l[2]*dx1; 
  div_b_comp[3] = (-0.5*b_r[3]*dx1)+0.5*b_l[3]*dx1+0.4330127018922193*b_r[2]*dx1+0.4330127018922193*b_l[2]*dx1-0.8660254037844386*b_c[2]*dx1; 

  div_b[0] += div_b_comp[0]; 
  div_b[1] += div_b_comp[1]; 
  div_b[2] += div_b_comp[2]; 
  div_b[3] += div_b_comp[3]; 

  bb_grad_u_comp[0] = 0.5*bxbz[3]*grad_u_z[3]+0.5*bxby[3]*grad_u_y[3]+0.5*bxbx[3]*grad_u_x[3]+0.5*bxbz[2]*grad_u_z[2]+0.5*bxby[2]*grad_u_y[2]+0.5*bxbx[2]*grad_u_x[2]+0.5*bxbz[1]*grad_u_z[1]+0.5*bxby[1]*grad_u_y[1]+0.5*bxbx[1]*grad_u_x[1]+0.5*bxbz[0]*grad_u_z[0]+0.5*bxby[0]*grad_u_y[0]+0.5*bxbx[0]*grad_u_x[0]; 
  bb_grad_u_comp[1] = 0.5*bxbz[2]*grad_u_z[3]+0.5*bxby[2]*grad_u_y[3]+0.5*bxbx[2]*grad_u_x[3]+0.5*grad_u_z[2]*bxbz[3]+0.5*grad_u_y[2]*bxby[3]+0.5*grad_u_x[2]*bxbx[3]+0.5*bxbz[0]*grad_u_z[1]+0.5*bxby[0]*grad_u_y[1]+0.5*bxbx[0]*grad_u_x[1]+0.5*grad_u_z[0]*bxbz[1]+0.5*grad_u_y[0]*bxby[1]+0.5*grad_u_x[0]*bxbx[1]; 
  bb_grad_u_comp[2] = 0.5*bxbz[1]*grad_u_z[3]+0.5*bxby[1]*grad_u_y[3]+0.5*bxbx[1]*grad_u_x[3]+0.5*grad_u_z[1]*bxbz[3]+0.5*grad_u_y[1]*bxby[3]+0.5*grad_u_x[1]*bxbx[3]+0.5*bxbz[0]*grad_u_z[2]+0.5*bxby[0]*grad_u_y[2]+0.5*bxbx[0]*grad_u_x[2]+0.5*grad_u_z[0]*bxbz[2]+0.5*grad_u_y[0]*bxby[2]+0.5*grad_u_x[0]*bxbx[2]; 
  bb_grad_u_comp[3] = 0.5*bxbz[0]*grad_u_z[3]+0.5*bxby[0]*grad_u_y[3]+0.5*bxbx[0]*grad_u_x[3]+0.5*grad_u_z[0]*bxbz[3]+0.5*grad_u_y[0]*bxby[3]+0.5*grad_u_x[0]*bxbx[3]+0.5*bxbz[1]*grad_u_z[2]+0.5*bxby[1]*grad_u_y[2]+0.5*bxbx[1]*grad_u_x[2]+0.5*grad_u_z[1]*bxbz[2]+0.5*grad_u_y[1]*bxby[2]+0.5*grad_u_x[1]*bxbx[2]; 

  bb_grad_u[0] += bb_grad_u_comp[0]; 
  bb_grad_u[1] += bb_grad_u_comp[1]; 
  bb_grad_u[2] += bb_grad_u_comp[2]; 
  bb_grad_u[3] += bb_grad_u_comp[3]; 

  p_force[0] += (-0.5*(T_perp_over_m[3]*div_b_comp[3]+T_perp_over_m[2]*div_b_comp[2]+T_perp_over_m[1]*div_b_comp[1]))+0.5*pkpm_div_ppar[0]-0.5*T_perp_over_m[0]*div_b_comp[0]; 
  p_force[1] += (-0.5*(T_perp_over_m[2]*div_b_comp[3]+div_b_comp[2]*T_perp_over_m[3]))+0.5*pkpm_div_ppar[1]-0.5*(T_perp_over_m[0]*div_b_comp[1]+div_b_comp[0]*T_perp_over_m[1]); 
  p_force[2] += (-0.5*(T_perp_over_m[1]*div_b_comp[3]+div_b_comp[1]*T_perp_over_m[3]))+0.5*pkpm_div_ppar[2]-0.5*(T_perp_over_m[0]*div_b_comp[2]+div_b_comp[0]*T_perp_over_m[2]); 
  p_force[3] += 0.5*pkpm_div_ppar[3]-0.5*(T_perp_over_m[0]*div_b_comp[3]+div_b_comp[0]*T_perp_over_m[3]+T_perp_over_m[1]*div_b_comp[2]+div_b_comp[1]*T_perp_over_m[2]); 

  p_perp_source[0] += bb_grad_u_comp[0]-1.0*(nu[0]+grad_u_x[0]); 
  p_perp_source[1] += bb_grad_u_comp[1]-1.0*(nu[1]+grad_u_x[1]); 
  p_perp_source[2] += bb_grad_u_comp[2]-1.0*(nu[2]+grad_u_x[2]); 
  p_perp_source[3] += bb_grad_u_comp[3]-1.0*(nu[3]+grad_u_x[3]); 

  p_perp_div_b[0] += 0.5*(T_perp_over_m[3]*div_b_comp[3]+T_perp_over_m[2]*div_b_comp[2]+T_perp_over_m[1]*div_b_comp[1]+T_perp_over_m[0]*div_b_comp[0]); 
  p_perp_div_b[1] += 0.5*(T_perp_over_m[2]*div_b_comp[3]+div_b_comp[2]*T_perp_over_m[3]+T_perp_over_m[0]*div_b_comp[1]+div_b_comp[0]*T_perp_over_m[1]); 
  p_perp_div_b[2] += 0.5*(T_perp_over_m[1]*div_b_comp[3]+div_b_comp[1]*T_perp_over_m[3]+T_perp_over_m[0]*div_b_comp[2]+div_b_comp[0]*T_perp_over_m[2]); 
  p_perp_div_b[3] += 0.5*(T_perp_over_m[0]*div_b_comp[3]+div_b_comp[0]*T_perp_over_m[3]+T_perp_over_m[1]*div_b_comp[2]+div_b_comp[1]*T_perp_over_m[2]); 

} 
