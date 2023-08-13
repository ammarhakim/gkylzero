#include <gkyl_euler_pkpm_kernels.h> 
#include <gkyl_basis_ser_3x_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH void pkpm_vars_accel_z_3x_ser_p1(const double *dxv, 
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

  const double dx1 = 2.0/dxv[2]; 
  const double *b_l = &bvar_l[16]; 
  const double *b_c = &bvar_c[16]; 
  const double *b_r = &bvar_r[16]; 

  const double *ux_c = &prim_c[0]; 
  const double *uy_c = &prim_c[8]; 
  const double *uz_c = &prim_c[16]; 

  const double *bxbx = &bvar_c[24]; 
  const double *bxby = &bvar_c[32]; 
  const double *bxbz = &bvar_c[40]; 
  const double *byby = &bvar_c[48]; 
  const double *bybz = &bvar_c[56]; 
  const double *bzbz = &bvar_c[64]; 

  const double *ux_surf_lr = &prim_surf_l[68]; 
  const double *uy_surf_lr = &prim_surf_l[76]; 
  const double *uz_surf_lr = &prim_surf_l[84]; 

  const double *ux_surf_cl = &prim_surf_c[64]; 
  const double *uy_surf_cl = &prim_surf_c[72]; 
  const double *uz_surf_cl = &prim_surf_c[80]; 

  const double *ux_surf_cr = &prim_surf_c[68]; 
  const double *uy_surf_cr = &prim_surf_c[76]; 
  const double *uz_surf_cr = &prim_surf_c[84]; 

  const double *ux_surf_rl = &prim_surf_r[64]; 
  const double *uy_surf_rl = &prim_surf_r[72]; 
  const double *uz_surf_rl = &prim_surf_r[80]; 

  const double *Tii_surf_lr = &prim_surf_l[92]; 
  const double *Tii_surf_cl = &prim_surf_c[88]; 
  const double *Tii_surf_cr = &prim_surf_c[92]; 
  const double *Tii_surf_rl = &prim_surf_r[88]; 

  const double *pkpm_div_ppar = &prim_c[24]; 
  const double *T_perp_over_m = &prim_c[32]; 

  double *pkpm_lax_l = &pkpm_lax[16]; 
  double *pkpm_lax_r = &pkpm_lax[20]; 

  double *div_b = &pkpm_accel[0]; 
  double *bb_grad_u = &pkpm_accel[8]; 
  double *p_force = &pkpm_accel[16]; 
  double *p_perp_source = &pkpm_accel[24]; 
  double *p_perp_div_b = &pkpm_accel[32]; 

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
  double pkpm_lax_quad_l[4] = {0.0}; 
  double pkpm_lax_quad_r[4] = {0.0}; 

  ul_r = 0.5*uz_surf_lr[3]-0.5*uz_surf_lr[2]-0.5*uz_surf_lr[1]+0.5*uz_surf_lr[0]; 
  uc_l = 0.5*uz_surf_cl[3]-0.5*uz_surf_cl[2]-0.5*uz_surf_cl[1]+0.5*uz_surf_cl[0]; 
  uc_r = 0.5*uz_surf_cr[3]-0.5*uz_surf_cr[2]-0.5*uz_surf_cr[1]+0.5*uz_surf_cr[0]; 
  ur_l = 0.5*uz_surf_rl[3]-0.5*uz_surf_rl[2]-0.5*uz_surf_rl[1]+0.5*uz_surf_rl[0]; 
  uQuad_l = fmax(fabs(ul_r), fabs(uc_l)); 
  uQuad_r = fmax(fabs(uc_r), fabs(ur_l)); 
  Tiil_r = 0.5*Tii_surf_lr[3]-0.5*Tii_surf_lr[2]-0.5*Tii_surf_lr[1]+0.5*Tii_surf_lr[0]; 
  Tiic_l = 0.5*Tii_surf_cl[3]-0.5*Tii_surf_cl[2]-0.5*Tii_surf_cl[1]+0.5*Tii_surf_cl[0]; 
  Tiic_r = 0.5*Tii_surf_cr[3]-0.5*Tii_surf_cr[2]-0.5*Tii_surf_cr[1]+0.5*Tii_surf_cr[0]; 
  Tiir_l = 0.5*Tii_surf_rl[3]-0.5*Tii_surf_rl[2]-0.5*Tii_surf_rl[1]+0.5*Tii_surf_rl[0]; 
  TiiQuad_l = fmax(sqrt(fabs(Tiil_r)), sqrt(fabs(Tiic_l))); 
  TiiQuad_r = fmax(sqrt(fabs(Tiic_r)), sqrt(fabs(Tiir_l))); 
  pkpm_lax_quad_l[0] = uQuad_l + TiiQuad_l; 
  pkpm_lax_quad_r[0] = uQuad_r + TiiQuad_r; 

  ul_r = (-0.5*uz_surf_lr[3])+0.5*uz_surf_lr[2]-0.5*uz_surf_lr[1]+0.5*uz_surf_lr[0]; 
  uc_l = (-0.5*uz_surf_cl[3])+0.5*uz_surf_cl[2]-0.5*uz_surf_cl[1]+0.5*uz_surf_cl[0]; 
  uc_r = (-0.5*uz_surf_cr[3])+0.5*uz_surf_cr[2]-0.5*uz_surf_cr[1]+0.5*uz_surf_cr[0]; 
  ur_l = (-0.5*uz_surf_rl[3])+0.5*uz_surf_rl[2]-0.5*uz_surf_rl[1]+0.5*uz_surf_rl[0]; 
  uQuad_l = fmax(fabs(ul_r), fabs(uc_l)); 
  uQuad_r = fmax(fabs(uc_r), fabs(ur_l)); 
  Tiil_r = (-0.5*Tii_surf_lr[3])+0.5*Tii_surf_lr[2]-0.5*Tii_surf_lr[1]+0.5*Tii_surf_lr[0]; 
  Tiic_l = (-0.5*Tii_surf_cl[3])+0.5*Tii_surf_cl[2]-0.5*Tii_surf_cl[1]+0.5*Tii_surf_cl[0]; 
  Tiic_r = (-0.5*Tii_surf_cr[3])+0.5*Tii_surf_cr[2]-0.5*Tii_surf_cr[1]+0.5*Tii_surf_cr[0]; 
  Tiir_l = (-0.5*Tii_surf_rl[3])+0.5*Tii_surf_rl[2]-0.5*Tii_surf_rl[1]+0.5*Tii_surf_rl[0]; 
  TiiQuad_l = fmax(sqrt(fabs(Tiil_r)), sqrt(fabs(Tiic_l))); 
  TiiQuad_r = fmax(sqrt(fabs(Tiic_r)), sqrt(fabs(Tiir_l))); 
  pkpm_lax_quad_l[1] = uQuad_l + TiiQuad_l; 
  pkpm_lax_quad_r[1] = uQuad_r + TiiQuad_r; 

  ul_r = (-0.5*uz_surf_lr[3])-0.5*uz_surf_lr[2]+0.5*uz_surf_lr[1]+0.5*uz_surf_lr[0]; 
  uc_l = (-0.5*uz_surf_cl[3])-0.5*uz_surf_cl[2]+0.5*uz_surf_cl[1]+0.5*uz_surf_cl[0]; 
  uc_r = (-0.5*uz_surf_cr[3])-0.5*uz_surf_cr[2]+0.5*uz_surf_cr[1]+0.5*uz_surf_cr[0]; 
  ur_l = (-0.5*uz_surf_rl[3])-0.5*uz_surf_rl[2]+0.5*uz_surf_rl[1]+0.5*uz_surf_rl[0]; 
  uQuad_l = fmax(fabs(ul_r), fabs(uc_l)); 
  uQuad_r = fmax(fabs(uc_r), fabs(ur_l)); 
  Tiil_r = (-0.5*Tii_surf_lr[3])-0.5*Tii_surf_lr[2]+0.5*Tii_surf_lr[1]+0.5*Tii_surf_lr[0]; 
  Tiic_l = (-0.5*Tii_surf_cl[3])-0.5*Tii_surf_cl[2]+0.5*Tii_surf_cl[1]+0.5*Tii_surf_cl[0]; 
  Tiic_r = (-0.5*Tii_surf_cr[3])-0.5*Tii_surf_cr[2]+0.5*Tii_surf_cr[1]+0.5*Tii_surf_cr[0]; 
  Tiir_l = (-0.5*Tii_surf_rl[3])-0.5*Tii_surf_rl[2]+0.5*Tii_surf_rl[1]+0.5*Tii_surf_rl[0]; 
  TiiQuad_l = fmax(sqrt(fabs(Tiil_r)), sqrt(fabs(Tiic_l))); 
  TiiQuad_r = fmax(sqrt(fabs(Tiic_r)), sqrt(fabs(Tiir_l))); 
  pkpm_lax_quad_l[2] = uQuad_l + TiiQuad_l; 
  pkpm_lax_quad_r[2] = uQuad_r + TiiQuad_r; 

  ul_r = 0.5*uz_surf_lr[3]+0.5*uz_surf_lr[2]+0.5*uz_surf_lr[1]+0.5*uz_surf_lr[0]; 
  uc_l = 0.5*uz_surf_cl[3]+0.5*uz_surf_cl[2]+0.5*uz_surf_cl[1]+0.5*uz_surf_cl[0]; 
  uc_r = 0.5*uz_surf_cr[3]+0.5*uz_surf_cr[2]+0.5*uz_surf_cr[1]+0.5*uz_surf_cr[0]; 
  ur_l = 0.5*uz_surf_rl[3]+0.5*uz_surf_rl[2]+0.5*uz_surf_rl[1]+0.5*uz_surf_rl[0]; 
  uQuad_l = fmax(fabs(ul_r), fabs(uc_l)); 
  uQuad_r = fmax(fabs(uc_r), fabs(ur_l)); 
  Tiil_r = 0.5*Tii_surf_lr[3]+0.5*Tii_surf_lr[2]+0.5*Tii_surf_lr[1]+0.5*Tii_surf_lr[0]; 
  Tiic_l = 0.5*Tii_surf_cl[3]+0.5*Tii_surf_cl[2]+0.5*Tii_surf_cl[1]+0.5*Tii_surf_cl[0]; 
  Tiic_r = 0.5*Tii_surf_cr[3]+0.5*Tii_surf_cr[2]+0.5*Tii_surf_cr[1]+0.5*Tii_surf_cr[0]; 
  Tiir_l = 0.5*Tii_surf_rl[3]+0.5*Tii_surf_rl[2]+0.5*Tii_surf_rl[1]+0.5*Tii_surf_rl[0]; 
  TiiQuad_l = fmax(sqrt(fabs(Tiil_r)), sqrt(fabs(Tiic_l))); 
  TiiQuad_r = fmax(sqrt(fabs(Tiic_r)), sqrt(fabs(Tiir_l))); 
  pkpm_lax_quad_l[3] = uQuad_l + TiiQuad_l; 
  pkpm_lax_quad_r[3] = uQuad_r + TiiQuad_r; 

  ser_3x_p1_upwind_quad_to_modal(pkpm_lax_quad_l, pkpm_lax_l); 
  ser_3x_p1_upwind_quad_to_modal(pkpm_lax_quad_r, pkpm_lax_r); 

  double grad_u_x[8] = {0.0}; 
  double grad_u_y[8] = {0.0}; 
  double grad_u_z[8] = {0.0}; 
  grad_u_x[0] = (0.3535533905932737*ux_surf_rl[0]-0.3535533905932737*ux_surf_lr[0]+0.3535533905932737*ux_surf_cr[0]-0.3535533905932737*ux_surf_cl[0])*dx1; 
  grad_u_x[1] = (0.3535533905932737*ux_surf_rl[1]-0.3535533905932737*ux_surf_lr[1]+0.3535533905932737*ux_surf_cr[1]-0.3535533905932737*ux_surf_cl[1])*dx1; 
  grad_u_x[2] = (0.3535533905932737*ux_surf_rl[2]-0.3535533905932737*ux_surf_lr[2]+0.3535533905932737*ux_surf_cr[2]-0.3535533905932737*ux_surf_cl[2])*dx1; 
  grad_u_x[3] = (0.6123724356957944*(ux_surf_rl[0]+ux_surf_lr[0]+ux_surf_cr[0]+ux_surf_cl[0])-1.732050807568877*ux_c[0])*dx1; 
  grad_u_x[4] = (0.3535533905932737*ux_surf_rl[3]-0.3535533905932737*ux_surf_lr[3]+0.3535533905932737*ux_surf_cr[3]-0.3535533905932737*ux_surf_cl[3])*dx1; 
  grad_u_x[5] = (0.6123724356957944*(ux_surf_rl[1]+ux_surf_lr[1]+ux_surf_cr[1]+ux_surf_cl[1])-1.732050807568877*ux_c[1])*dx1; 
  grad_u_x[6] = (0.6123724356957944*(ux_surf_rl[2]+ux_surf_lr[2]+ux_surf_cr[2]+ux_surf_cl[2])-1.732050807568877*ux_c[2])*dx1; 
  grad_u_x[7] = (0.6123724356957944*(ux_surf_rl[3]+ux_surf_lr[3]+ux_surf_cr[3]+ux_surf_cl[3])-1.732050807568877*ux_c[4])*dx1; 

  grad_u_y[0] = (0.3535533905932737*uy_surf_rl[0]-0.3535533905932737*uy_surf_lr[0]+0.3535533905932737*uy_surf_cr[0]-0.3535533905932737*uy_surf_cl[0])*dx1; 
  grad_u_y[1] = (0.3535533905932737*uy_surf_rl[1]-0.3535533905932737*uy_surf_lr[1]+0.3535533905932737*uy_surf_cr[1]-0.3535533905932737*uy_surf_cl[1])*dx1; 
  grad_u_y[2] = (0.3535533905932737*uy_surf_rl[2]-0.3535533905932737*uy_surf_lr[2]+0.3535533905932737*uy_surf_cr[2]-0.3535533905932737*uy_surf_cl[2])*dx1; 
  grad_u_y[3] = (0.6123724356957944*(uy_surf_rl[0]+uy_surf_lr[0]+uy_surf_cr[0]+uy_surf_cl[0])-1.732050807568877*uy_c[0])*dx1; 
  grad_u_y[4] = (0.3535533905932737*uy_surf_rl[3]-0.3535533905932737*uy_surf_lr[3]+0.3535533905932737*uy_surf_cr[3]-0.3535533905932737*uy_surf_cl[3])*dx1; 
  grad_u_y[5] = (0.6123724356957944*(uy_surf_rl[1]+uy_surf_lr[1]+uy_surf_cr[1]+uy_surf_cl[1])-1.732050807568877*uy_c[1])*dx1; 
  grad_u_y[6] = (0.6123724356957944*(uy_surf_rl[2]+uy_surf_lr[2]+uy_surf_cr[2]+uy_surf_cl[2])-1.732050807568877*uy_c[2])*dx1; 
  grad_u_y[7] = (0.6123724356957944*(uy_surf_rl[3]+uy_surf_lr[3]+uy_surf_cr[3]+uy_surf_cl[3])-1.732050807568877*uy_c[4])*dx1; 

  grad_u_z[0] = (0.3535533905932737*uz_surf_rl[0]-0.3535533905932737*uz_surf_lr[0]+0.3535533905932737*uz_surf_cr[0]-0.3535533905932737*uz_surf_cl[0])*dx1; 
  grad_u_z[1] = (0.3535533905932737*uz_surf_rl[1]-0.3535533905932737*uz_surf_lr[1]+0.3535533905932737*uz_surf_cr[1]-0.3535533905932737*uz_surf_cl[1])*dx1; 
  grad_u_z[2] = (0.3535533905932737*uz_surf_rl[2]-0.3535533905932737*uz_surf_lr[2]+0.3535533905932737*uz_surf_cr[2]-0.3535533905932737*uz_surf_cl[2])*dx1; 
  grad_u_z[3] = (0.6123724356957944*(uz_surf_rl[0]+uz_surf_lr[0]+uz_surf_cr[0]+uz_surf_cl[0])-1.732050807568877*uz_c[0])*dx1; 
  grad_u_z[4] = (0.3535533905932737*uz_surf_rl[3]-0.3535533905932737*uz_surf_lr[3]+0.3535533905932737*uz_surf_cr[3]-0.3535533905932737*uz_surf_cl[3])*dx1; 
  grad_u_z[5] = (0.6123724356957944*(uz_surf_rl[1]+uz_surf_lr[1]+uz_surf_cr[1]+uz_surf_cl[1])-1.732050807568877*uz_c[1])*dx1; 
  grad_u_z[6] = (0.6123724356957944*(uz_surf_rl[2]+uz_surf_lr[2]+uz_surf_cr[2]+uz_surf_cl[2])-1.732050807568877*uz_c[2])*dx1; 
  grad_u_z[7] = (0.6123724356957944*(uz_surf_rl[3]+uz_surf_lr[3]+uz_surf_cr[3]+uz_surf_cl[3])-1.732050807568877*uz_c[4])*dx1; 

  double div_b_comp[8] = {0.0}; 
  double bb_grad_u_comp[8] = {0.0}; 
  div_b_comp[0] = ((-0.2886751345948129*(b_r[3]+b_l[3]))+0.5773502691896258*b_c[3]+0.25*b_r[0]-0.25*b_l[0])*dx1; 
  div_b_comp[1] = ((-0.2886751345948129*(b_r[5]+b_l[5]))+0.5773502691896258*b_c[5]+0.25*b_r[1]-0.25*b_l[1])*dx1; 
  div_b_comp[2] = ((-0.2886751345948129*(b_r[6]+b_l[6]))+0.5773502691896258*b_c[6]+0.25*b_r[2]-0.25*b_l[2])*dx1; 
  div_b_comp[3] = ((-0.5*b_r[3])+0.5*b_l[3]+0.4330127018922193*(b_r[0]+b_l[0])-0.8660254037844386*b_c[0])*dx1; 
  div_b_comp[4] = ((-0.2886751345948129*(b_r[7]+b_l[7]))+0.5773502691896258*b_c[7]+0.25*b_r[4]-0.25*b_l[4])*dx1; 
  div_b_comp[5] = ((-0.5*b_r[5])+0.5*b_l[5]+0.4330127018922193*(b_r[1]+b_l[1])-0.8660254037844386*b_c[1])*dx1; 
  div_b_comp[6] = ((-0.5*b_r[6])+0.5*b_l[6]+0.4330127018922193*(b_r[2]+b_l[2])-0.8660254037844386*b_c[2])*dx1; 
  div_b_comp[7] = ((-0.5*b_r[7])+0.5*b_l[7]+0.4330127018922193*(b_r[4]+b_l[4])-0.8660254037844386*b_c[4])*dx1; 

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

  p_perp_source[0] += (-0.6666666666666666*nu_c[0])-1.0*grad_u_z[0]+bb_grad_u_comp[0]; 
  p_perp_source[1] += (-0.6666666666666666*nu_c[1])-1.0*grad_u_z[1]+bb_grad_u_comp[1]; 
  p_perp_source[2] += (-0.6666666666666666*nu_c[2])-1.0*grad_u_z[2]+bb_grad_u_comp[2]; 
  p_perp_source[3] += (-0.6666666666666666*nu_c[3])-1.0*grad_u_z[3]+bb_grad_u_comp[3]; 
  p_perp_source[4] += (-0.6666666666666666*nu_c[4])-1.0*grad_u_z[4]+bb_grad_u_comp[4]; 
  p_perp_source[5] += (-0.6666666666666666*nu_c[5])-1.0*grad_u_z[5]+bb_grad_u_comp[5]; 
  p_perp_source[6] += (-0.6666666666666666*nu_c[6])-1.0*grad_u_z[6]+bb_grad_u_comp[6]; 
  p_perp_source[7] += (-0.6666666666666666*nu_c[7])-1.0*grad_u_z[7]+bb_grad_u_comp[7]; 

  p_perp_div_b[0] += 0.3535533905932737*(T_perp_over_m[7]*div_b_comp[7]+T_perp_over_m[6]*div_b_comp[6]+T_perp_over_m[5]*div_b_comp[5]+T_perp_over_m[4]*div_b_comp[4]+T_perp_over_m[3]*div_b_comp[3]+T_perp_over_m[2]*div_b_comp[2]+T_perp_over_m[1]*div_b_comp[1]+T_perp_over_m[0]*div_b_comp[0]); 
  p_perp_div_b[1] += 0.3535533905932737*(T_perp_over_m[6]*div_b_comp[7]+div_b_comp[6]*T_perp_over_m[7]+T_perp_over_m[3]*div_b_comp[5]+div_b_comp[3]*T_perp_over_m[5]+T_perp_over_m[2]*div_b_comp[4]+div_b_comp[2]*T_perp_over_m[4]+T_perp_over_m[0]*div_b_comp[1]+div_b_comp[0]*T_perp_over_m[1]); 
  p_perp_div_b[2] += 0.3535533905932737*(T_perp_over_m[5]*div_b_comp[7]+div_b_comp[5]*T_perp_over_m[7]+T_perp_over_m[3]*div_b_comp[6]+div_b_comp[3]*T_perp_over_m[6]+T_perp_over_m[1]*div_b_comp[4]+div_b_comp[1]*T_perp_over_m[4]+T_perp_over_m[0]*div_b_comp[2]+div_b_comp[0]*T_perp_over_m[2]); 
  p_perp_div_b[3] += 0.3535533905932737*(T_perp_over_m[4]*div_b_comp[7]+div_b_comp[4]*T_perp_over_m[7]+T_perp_over_m[2]*div_b_comp[6]+div_b_comp[2]*T_perp_over_m[6]+T_perp_over_m[1]*div_b_comp[5]+div_b_comp[1]*T_perp_over_m[5]+T_perp_over_m[0]*div_b_comp[3]+div_b_comp[0]*T_perp_over_m[3]); 
  p_perp_div_b[4] += 0.3535533905932737*(T_perp_over_m[3]*div_b_comp[7]+div_b_comp[3]*T_perp_over_m[7]+T_perp_over_m[5]*div_b_comp[6]+div_b_comp[5]*T_perp_over_m[6]+T_perp_over_m[0]*div_b_comp[4]+div_b_comp[0]*T_perp_over_m[4]+T_perp_over_m[1]*div_b_comp[2]+div_b_comp[1]*T_perp_over_m[2]); 
  p_perp_div_b[5] += 0.3535533905932737*(T_perp_over_m[2]*div_b_comp[7]+div_b_comp[2]*T_perp_over_m[7]+T_perp_over_m[4]*div_b_comp[6]+div_b_comp[4]*T_perp_over_m[6]+T_perp_over_m[0]*div_b_comp[5]+div_b_comp[0]*T_perp_over_m[5]+T_perp_over_m[1]*div_b_comp[3]+div_b_comp[1]*T_perp_over_m[3]); 
  p_perp_div_b[6] += 0.3535533905932737*(T_perp_over_m[1]*div_b_comp[7]+div_b_comp[1]*T_perp_over_m[7]+T_perp_over_m[0]*div_b_comp[6]+div_b_comp[0]*T_perp_over_m[6]+T_perp_over_m[4]*div_b_comp[5]+div_b_comp[4]*T_perp_over_m[5]+T_perp_over_m[2]*div_b_comp[3]+div_b_comp[2]*T_perp_over_m[3]); 
  p_perp_div_b[7] += 0.3535533905932737*(T_perp_over_m[0]*div_b_comp[7]+div_b_comp[0]*T_perp_over_m[7]+T_perp_over_m[1]*div_b_comp[6]+div_b_comp[1]*T_perp_over_m[6]+T_perp_over_m[2]*div_b_comp[5]+div_b_comp[2]*T_perp_over_m[5]+T_perp_over_m[3]*div_b_comp[4]+div_b_comp[3]*T_perp_over_m[4]); 

} 
