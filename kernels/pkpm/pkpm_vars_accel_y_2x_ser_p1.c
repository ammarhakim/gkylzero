#include <gkyl_euler_pkpm_kernels.h> 
#include <gkyl_basis_ser_2x_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH void pkpm_vars_accel_y_2x_ser_p1(const double *dxv, 
  const double *bvar_l, const double *bvar_c, const double *bvar_r, 
  const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
  const double *prim_c, const double *nu_c, 
  double* GKYL_RESTRICT pkpm_lax, double* GKYL_RESTRICT pkpm_accel) 
{ 
  // dxv[NDIM]:      Cell spacing.
  // bvar_l/c/r:      Input magnetic field unit vector in left/center/right cells.
  // prim_surf_l/c/r: Input surface primitive variables [u_i, 3*T_ii/m] in left/center/right cells in each direction.
  // primc:          Input volume expansion of primitive variables [ux, uy, uz, 1/rho div(p_par b), T_perp/m, m/T_perp] in center cell.
  // nuc:            Input volume expansion of collisionality in center cell.
  // pkpm_lax:       Surface expansion of pkpm Lax penalization: lambda_i = |u_i| + sqrt(3.0*T_ii/m).
  // pkpm_accel:     Volume expansion of pkpm acceleration variables.

  const double dx1 = 2.0/dxv[1]; 
  const double *b_l = &bvar_l[4]; 
  const double *b_c = &bvar_c[4]; 
  const double *b_r = &bvar_r[4]; 

  const double *ux_c = &prim_c[0]; 
  const double *uy_c = &prim_c[4]; 
  const double *uz_c = &prim_c[8]; 

  const double *bxbx = &bvar_c[12]; 
  const double *bxby = &bvar_c[16]; 
  const double *bxbz = &bvar_c[20]; 
  const double *byby = &bvar_c[24]; 
  const double *bybz = &bvar_c[28]; 
  const double *bzbz = &bvar_c[32]; 

  const double *ux_surf_lr = &prim_surf_l[18]; 
  const double *uy_surf_lr = &prim_surf_l[22]; 
  const double *uz_surf_lr = &prim_surf_l[26]; 

  const double *ux_surf_cl = &prim_surf_c[16]; 
  const double *uy_surf_cl = &prim_surf_c[20]; 
  const double *uz_surf_cl = &prim_surf_c[24]; 

  const double *ux_surf_cr = &prim_surf_c[18]; 
  const double *uy_surf_cr = &prim_surf_c[22]; 
  const double *uz_surf_cr = &prim_surf_c[26]; 

  const double *ux_surf_rl = &prim_surf_r[16]; 
  const double *uy_surf_rl = &prim_surf_r[20]; 
  const double *uz_surf_rl = &prim_surf_r[24]; 

  const double *Tii_surf_lr = &prim_surf_l[30]; 
  const double *Tii_surf_cl = &prim_surf_c[28]; 
  const double *Tii_surf_cr = &prim_surf_c[30]; 
  const double *Tii_surf_rl = &prim_surf_r[28]; 

  const double *pkpm_div_ppar = &prim_c[12]; 
  const double *T_perp_over_m = &prim_c[16]; 

  double *pkpm_lax_l = &pkpm_lax[4]; 
  double *pkpm_lax_r = &pkpm_lax[6]; 

  double *div_b = &pkpm_accel[0]; 
  double *bb_grad_u = &pkpm_accel[4]; 
  double *p_force = &pkpm_accel[8]; 
  double *p_perp_source = &pkpm_accel[12]; 
  double *p_perp_div_b = &pkpm_accel[16]; 

  double ul_r = 0.0; 
  double uc_l = 0.0; 
  double uc_r = 0.0; 
  double ur_l = 0.0; 
  double uQuad_l = 0.0; 
  double uQuad_r = 0.0; 
  double Tiil_r = 0.0; 
  double Tiic_l = 0.0; 
  double Tiic_r = 0.0; 
  double Tiir_l = 0.0; 
  double TiiQuad_l = 0.0; 
  double TiiQuad_r = 0.0; 
  double pkpm_lax_quad_l[2] = {0.0}; 
  double pkpm_lax_quad_r[2] = {0.0}; 

  ul_r = 0.7071067811865475*uy_surf_lr[0]-0.7071067811865475*uy_surf_lr[1]; 
  uc_l = 0.7071067811865475*uy_surf_cl[0]-0.7071067811865475*uy_surf_cl[1]; 
  uc_r = 0.7071067811865475*uy_surf_cr[0]-0.7071067811865475*uy_surf_cr[1]; 
  ur_l = 0.7071067811865475*uy_surf_rl[0]-0.7071067811865475*uy_surf_rl[1]; 
  uQuad_l = fmax(fabs(ul_r), fabs(uc_l)); 
  uQuad_r = fmax(fabs(uc_r), fabs(ur_l)); 
  Tiil_r = 0.7071067811865475*Tii_surf_lr[0]-0.7071067811865475*Tii_surf_lr[1]; 
  Tiic_l = 0.7071067811865475*Tii_surf_cl[0]-0.7071067811865475*Tii_surf_cl[1]; 
  Tiic_r = 0.7071067811865475*Tii_surf_cr[0]-0.7071067811865475*Tii_surf_cr[1]; 
  Tiir_l = 0.7071067811865475*Tii_surf_rl[0]-0.7071067811865475*Tii_surf_rl[1]; 
  TiiQuad_l = fmax(sqrt(fabs(Tiil_r)), sqrt(fabs(Tiic_l))); 
  TiiQuad_r = fmax(sqrt(fabs(Tiic_r)), sqrt(fabs(Tiir_l))); 
  pkpm_lax_quad_l[0] = uQuad_l + TiiQuad_l; 
  pkpm_lax_quad_r[0] = uQuad_r + TiiQuad_r; 

  ul_r = 0.7071067811865475*uy_surf_lr[1]+0.7071067811865475*uy_surf_lr[0]; 
  uc_l = 0.7071067811865475*uy_surf_cl[1]+0.7071067811865475*uy_surf_cl[0]; 
  uc_r = 0.7071067811865475*uy_surf_cr[1]+0.7071067811865475*uy_surf_cr[0]; 
  ur_l = 0.7071067811865475*uy_surf_rl[1]+0.7071067811865475*uy_surf_rl[0]; 
  uQuad_l = fmax(fabs(ul_r), fabs(uc_l)); 
  uQuad_r = fmax(fabs(uc_r), fabs(ur_l)); 
  Tiil_r = 0.7071067811865475*Tii_surf_lr[1]+0.7071067811865475*Tii_surf_lr[0]; 
  Tiic_l = 0.7071067811865475*Tii_surf_cl[1]+0.7071067811865475*Tii_surf_cl[0]; 
  Tiic_r = 0.7071067811865475*Tii_surf_cr[1]+0.7071067811865475*Tii_surf_cr[0]; 
  Tiir_l = 0.7071067811865475*Tii_surf_rl[1]+0.7071067811865475*Tii_surf_rl[0]; 
  TiiQuad_l = fmax(sqrt(fabs(Tiil_r)), sqrt(fabs(Tiic_l))); 
  TiiQuad_r = fmax(sqrt(fabs(Tiic_r)), sqrt(fabs(Tiir_l))); 
  pkpm_lax_quad_l[1] = uQuad_l + TiiQuad_l; 
  pkpm_lax_quad_r[1] = uQuad_r + TiiQuad_r; 

  ser_2x_p1_upwind_quad_to_modal(pkpm_lax_quad_l, pkpm_lax_l); 
  ser_2x_p1_upwind_quad_to_modal(pkpm_lax_quad_r, pkpm_lax_r); 

  double grad_u_x[4] = {0.0}; 
  double grad_u_y[4] = {0.0}; 
  double grad_u_z[4] = {0.0}; 
  grad_u_x[0] = (0.3535533905932737*ux_surf_rl[0]-0.3535533905932737*ux_surf_lr[0]+0.3535533905932737*ux_surf_cr[0]-0.3535533905932737*ux_surf_cl[0])*dx1; 
  grad_u_x[1] = (0.3535533905932737*ux_surf_rl[1]-0.3535533905932737*ux_surf_lr[1]+0.3535533905932737*ux_surf_cr[1]-0.3535533905932737*ux_surf_cl[1])*dx1; 
  grad_u_x[2] = (0.6123724356957944*(ux_surf_rl[0]+ux_surf_lr[0]+ux_surf_cr[0]+ux_surf_cl[0])-1.732050807568877*ux_c[0])*dx1; 
  grad_u_x[3] = (0.6123724356957944*(ux_surf_rl[1]+ux_surf_lr[1]+ux_surf_cr[1]+ux_surf_cl[1])-1.732050807568877*ux_c[1])*dx1; 

  grad_u_y[0] = (0.3535533905932737*uy_surf_rl[0]-0.3535533905932737*uy_surf_lr[0]+0.3535533905932737*uy_surf_cr[0]-0.3535533905932737*uy_surf_cl[0])*dx1; 
  grad_u_y[1] = (0.3535533905932737*uy_surf_rl[1]-0.3535533905932737*uy_surf_lr[1]+0.3535533905932737*uy_surf_cr[1]-0.3535533905932737*uy_surf_cl[1])*dx1; 
  grad_u_y[2] = (0.6123724356957944*(uy_surf_rl[0]+uy_surf_lr[0]+uy_surf_cr[0]+uy_surf_cl[0])-1.732050807568877*uy_c[0])*dx1; 
  grad_u_y[3] = (0.6123724356957944*(uy_surf_rl[1]+uy_surf_lr[1]+uy_surf_cr[1]+uy_surf_cl[1])-1.732050807568877*uy_c[1])*dx1; 

  grad_u_z[0] = (0.3535533905932737*uz_surf_rl[0]-0.3535533905932737*uz_surf_lr[0]+0.3535533905932737*uz_surf_cr[0]-0.3535533905932737*uz_surf_cl[0])*dx1; 
  grad_u_z[1] = (0.3535533905932737*uz_surf_rl[1]-0.3535533905932737*uz_surf_lr[1]+0.3535533905932737*uz_surf_cr[1]-0.3535533905932737*uz_surf_cl[1])*dx1; 
  grad_u_z[2] = (0.6123724356957944*(uz_surf_rl[0]+uz_surf_lr[0]+uz_surf_cr[0]+uz_surf_cl[0])-1.732050807568877*uz_c[0])*dx1; 
  grad_u_z[3] = (0.6123724356957944*(uz_surf_rl[1]+uz_surf_lr[1]+uz_surf_cr[1]+uz_surf_cl[1])-1.732050807568877*uz_c[1])*dx1; 

  double div_b_comp[4] = {0.0}; 
  double bb_grad_u_comp[4] = {0.0}; 
  div_b_comp[0] = ((-0.2886751345948129*(b_r[2]+b_l[2]))+0.5773502691896258*b_c[2]+0.25*b_r[0]-0.25*b_l[0])*dx1; 
  div_b_comp[1] = ((-0.2886751345948129*(b_r[3]+b_l[3]))+0.5773502691896258*b_c[3]+0.25*b_r[1]-0.25*b_l[1])*dx1; 
  div_b_comp[2] = ((-0.5*b_r[2])+0.5*b_l[2]+0.4330127018922193*(b_r[0]+b_l[0])-0.8660254037844386*b_c[0])*dx1; 
  div_b_comp[3] = ((-0.5*b_r[3])+0.5*b_l[3]+0.4330127018922193*(b_r[1]+b_l[1])-0.8660254037844386*b_c[1])*dx1; 

  div_b[0] += div_b_comp[0]; 
  div_b[1] += div_b_comp[1]; 
  div_b[2] += div_b_comp[2]; 
  div_b[3] += div_b_comp[3]; 

  bb_grad_u_comp[0] = 0.5*bybz[3]*grad_u_z[3]+0.5*byby[3]*grad_u_y[3]+0.5*bxby[3]*grad_u_x[3]+0.5*bybz[2]*grad_u_z[2]+0.5*byby[2]*grad_u_y[2]+0.5*bxby[2]*grad_u_x[2]+0.5*bybz[1]*grad_u_z[1]+0.5*byby[1]*grad_u_y[1]+0.5*bxby[1]*grad_u_x[1]+0.5*bybz[0]*grad_u_z[0]+0.5*byby[0]*grad_u_y[0]+0.5*bxby[0]*grad_u_x[0]; 
  bb_grad_u_comp[1] = 0.5*bybz[2]*grad_u_z[3]+0.5*byby[2]*grad_u_y[3]+0.5*bxby[2]*grad_u_x[3]+0.5*grad_u_z[2]*bybz[3]+0.5*grad_u_y[2]*byby[3]+0.5*grad_u_x[2]*bxby[3]+0.5*bybz[0]*grad_u_z[1]+0.5*byby[0]*grad_u_y[1]+0.5*bxby[0]*grad_u_x[1]+0.5*grad_u_z[0]*bybz[1]+0.5*grad_u_y[0]*byby[1]+0.5*grad_u_x[0]*bxby[1]; 
  bb_grad_u_comp[2] = 0.5*bybz[1]*grad_u_z[3]+0.5*byby[1]*grad_u_y[3]+0.5*bxby[1]*grad_u_x[3]+0.5*grad_u_z[1]*bybz[3]+0.5*grad_u_y[1]*byby[3]+0.5*grad_u_x[1]*bxby[3]+0.5*bybz[0]*grad_u_z[2]+0.5*byby[0]*grad_u_y[2]+0.5*bxby[0]*grad_u_x[2]+0.5*grad_u_z[0]*bybz[2]+0.5*grad_u_y[0]*byby[2]+0.5*grad_u_x[0]*bxby[2]; 
  bb_grad_u_comp[3] = 0.5*bybz[0]*grad_u_z[3]+0.5*byby[0]*grad_u_y[3]+0.5*bxby[0]*grad_u_x[3]+0.5*grad_u_z[0]*bybz[3]+0.5*grad_u_y[0]*byby[3]+0.5*grad_u_x[0]*bxby[3]+0.5*bybz[1]*grad_u_z[2]+0.5*byby[1]*grad_u_y[2]+0.5*bxby[1]*grad_u_x[2]+0.5*grad_u_z[1]*bybz[2]+0.5*grad_u_y[1]*byby[2]+0.5*grad_u_x[1]*bxby[2]; 

  bb_grad_u[0] += bb_grad_u_comp[0]; 
  bb_grad_u[1] += bb_grad_u_comp[1]; 
  bb_grad_u[2] += bb_grad_u_comp[2]; 
  bb_grad_u[3] += bb_grad_u_comp[3]; 

  p_force[0] += (-0.5*(T_perp_over_m[3]*div_b_comp[3]+T_perp_over_m[2]*div_b_comp[2]+T_perp_over_m[1]*div_b_comp[1]))+0.5*pkpm_div_ppar[0]-0.5*T_perp_over_m[0]*div_b_comp[0]; 
  p_force[1] += (-0.5*(T_perp_over_m[2]*div_b_comp[3]+div_b_comp[2]*T_perp_over_m[3]))+0.5*pkpm_div_ppar[1]-0.5*(T_perp_over_m[0]*div_b_comp[1]+div_b_comp[0]*T_perp_over_m[1]); 
  p_force[2] += (-0.5*(T_perp_over_m[1]*div_b_comp[3]+div_b_comp[1]*T_perp_over_m[3]))+0.5*pkpm_div_ppar[2]-0.5*(T_perp_over_m[0]*div_b_comp[2]+div_b_comp[0]*T_perp_over_m[2]); 
  p_force[3] += 0.5*pkpm_div_ppar[3]-0.5*(T_perp_over_m[0]*div_b_comp[3]+div_b_comp[0]*T_perp_over_m[3]+T_perp_over_m[1]*div_b_comp[2]+div_b_comp[1]*T_perp_over_m[2]); 

  p_perp_source[0] += bb_grad_u_comp[0]-1.0*(nu_c[0]+grad_u_y[0]); 
  p_perp_source[1] += bb_grad_u_comp[1]-1.0*(nu_c[1]+grad_u_y[1]); 
  p_perp_source[2] += bb_grad_u_comp[2]-1.0*(nu_c[2]+grad_u_y[2]); 
  p_perp_source[3] += bb_grad_u_comp[3]-1.0*(nu_c[3]+grad_u_y[3]); 

  p_perp_div_b[0] += 0.5*(T_perp_over_m[3]*div_b_comp[3]+T_perp_over_m[2]*div_b_comp[2]+T_perp_over_m[1]*div_b_comp[1]+T_perp_over_m[0]*div_b_comp[0]); 
  p_perp_div_b[1] += 0.5*(T_perp_over_m[2]*div_b_comp[3]+div_b_comp[2]*T_perp_over_m[3]+T_perp_over_m[0]*div_b_comp[1]+div_b_comp[0]*T_perp_over_m[1]); 
  p_perp_div_b[2] += 0.5*(T_perp_over_m[1]*div_b_comp[3]+div_b_comp[1]*T_perp_over_m[3]+T_perp_over_m[0]*div_b_comp[2]+div_b_comp[0]*T_perp_over_m[2]); 
  p_perp_div_b[3] += 0.5*(T_perp_over_m[0]*div_b_comp[3]+div_b_comp[0]*T_perp_over_m[3]+T_perp_over_m[1]*div_b_comp[2]+div_b_comp[1]*T_perp_over_m[2]); 

} 
