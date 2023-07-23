#include <gkyl_euler_pkpm_kernels.h> 
GKYL_CU_DH void pkpm_vars_accel_x_1x_ser_p1(const double *dxv, 
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
  const double *uy_l = &priml[2]; 
  const double *uz_l = &priml[4]; 

  const double *ux_c = &primc[0]; 
  const double *uy_c = &primc[2]; 
  const double *uz_c = &primc[4]; 

  const double *ux_r = &primr[0]; 
  const double *uy_r = &primr[2]; 
  const double *uz_r = &primr[4]; 

  const double *bxbx = &bvarc[6]; 
  const double *bxby = &bvarc[8]; 
  const double *bxbz = &bvarc[10]; 
  const double *byby = &bvarc[12]; 
  const double *bybz = &bvarc[14]; 
  const double *bzbz = &bvarc[16]; 

  const double *pkpm_div_ppar = &primc[12]; 
  const double *T_perp_over_m = &primc[14]; 

  double *div_b = &pkpm_accel[0]; 
  double *bb_grad_u = &pkpm_accel[2]; 
  double *p_force = &pkpm_accel[4]; 
  double *p_perp_source = &pkpm_accel[6]; 
  double *p_perp_div_b = &pkpm_accel[8]; 

  double grad_u_x[2] = {0.0}; 
  double grad_u_y[2] = {0.0}; 
  double grad_u_z[2] = {0.0}; 
  grad_u_x[0] = (-0.4330127018922193*ux_r[1]*dx1)-0.4330127018922193*ux_l[1]*dx1+0.8660254037844386*ux_c[1]*dx1+0.25*ux_r[0]*dx1-0.25*ux_l[0]*dx1; 
  grad_u_x[1] = (-0.75*ux_r[1]*dx1)+0.75*ux_l[1]*dx1+0.4330127018922193*ux_r[0]*dx1+0.4330127018922193*ux_l[0]*dx1-0.8660254037844385*ux_c[0]*dx1; 

  grad_u_y[0] = (-0.4330127018922193*uy_r[1]*dx1)-0.4330127018922193*uy_l[1]*dx1+0.8660254037844386*uy_c[1]*dx1+0.25*uy_r[0]*dx1-0.25*uy_l[0]*dx1; 
  grad_u_y[1] = (-0.75*uy_r[1]*dx1)+0.75*uy_l[1]*dx1+0.4330127018922193*uy_r[0]*dx1+0.4330127018922193*uy_l[0]*dx1-0.8660254037844385*uy_c[0]*dx1; 

  grad_u_z[0] = (-0.4330127018922193*uz_r[1]*dx1)-0.4330127018922193*uz_l[1]*dx1+0.8660254037844386*uz_c[1]*dx1+0.25*uz_r[0]*dx1-0.25*uz_l[0]*dx1; 
  grad_u_z[1] = (-0.75*uz_r[1]*dx1)+0.75*uz_l[1]*dx1+0.4330127018922193*uz_r[0]*dx1+0.4330127018922193*uz_l[0]*dx1-0.8660254037844385*uz_c[0]*dx1; 

  double div_b_comp[2] = {0.0}; 
  double div_p_x_comp[2] = {0.0}; 
  double div_p_y_comp[2] = {0.0}; 
  double div_p_z_comp[2] = {0.0}; 
  double bb_grad_u_comp[2] = {0.0}; 
  div_b_comp[0] = (-0.2886751345948129*b_r[1]*dx1)-0.2886751345948129*b_l[1]*dx1+0.5773502691896258*b_c[1]*dx1+0.25*b_r[0]*dx1-0.25*b_l[0]*dx1; 
  div_b_comp[1] = (-0.5*b_r[1]*dx1)+0.5*b_l[1]*dx1+0.4330127018922193*b_r[0]*dx1+0.4330127018922193*b_l[0]*dx1-0.8660254037844386*b_c[0]*dx1; 

  div_b[0] += div_b_comp[0]; 
  div_b[1] += div_b_comp[1]; 

  bb_grad_u_comp[0] = 0.7071067811865475*bxbz[1]*grad_u_z[1]+0.7071067811865475*bxby[1]*grad_u_y[1]+0.7071067811865475*bxbx[1]*grad_u_x[1]+0.7071067811865475*bxbz[0]*grad_u_z[0]+0.7071067811865475*bxby[0]*grad_u_y[0]+0.7071067811865475*bxbx[0]*grad_u_x[0]; 
  bb_grad_u_comp[1] = 0.7071067811865475*bxbz[0]*grad_u_z[1]+0.7071067811865475*bxby[0]*grad_u_y[1]+0.7071067811865475*bxbx[0]*grad_u_x[1]+0.7071067811865475*grad_u_z[0]*bxbz[1]+0.7071067811865475*grad_u_y[0]*bxby[1]+0.7071067811865475*grad_u_x[0]*bxbx[1]; 

  bb_grad_u[0] += bb_grad_u_comp[0]; 
  bb_grad_u[1] += bb_grad_u_comp[1]; 

  p_force[0] += (-0.7071067811865475*T_perp_over_m[1]*div_b_comp[1])+pkpm_div_ppar[0]-0.7071067811865475*T_perp_over_m[0]*div_b_comp[0]; 
  p_force[1] += pkpm_div_ppar[1]-0.7071067811865475*(T_perp_over_m[0]*div_b_comp[1]+div_b_comp[0]*T_perp_over_m[1]); 

  p_perp_source[0] += (-2.0*nu[0])-1.0*grad_u_x[0]+bb_grad_u_comp[0]; 
  p_perp_source[1] += (-2.0*nu[1])-1.0*grad_u_x[1]+bb_grad_u_comp[1]; 

  p_perp_div_b[0] += 0.7071067811865475*(T_perp_over_m[1]*div_b_comp[1]+T_perp_over_m[0]*div_b_comp[0]); 
  p_perp_div_b[1] += 0.7071067811865475*(T_perp_over_m[0]*div_b_comp[1]+div_b_comp[0]*T_perp_over_m[1]); 

} 
