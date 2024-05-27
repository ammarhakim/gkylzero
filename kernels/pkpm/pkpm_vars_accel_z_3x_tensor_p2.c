#include <gkyl_euler_pkpm_kernels.h> 
#include <gkyl_basis_tensor_3x_p2_surfx3_eval_quad.h> 
#include <gkyl_basis_tensor_3x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH void pkpm_vars_accel_z_3x_tensor_p2(const double *dxv, 
  const double *u_surf_l, const double *u_surf_c, const double *u_surf_r, 
  const double *prim_l, const double *prim_c, const double *prim_r, 
  const double *pkpm_u_c, const double *bvar_c, const double *nu_c, 
  double* GKYL_RESTRICT pkpm_lax, double* GKYL_RESTRICT pkpm_accel) 
{ 
  // dxv[NDIM]:       Cell spacing.
  // u_surf_l/c/r: Input surface flow velocity expansion in left/center/right cells in each direction.
  //               [ux_xl, ux_xr, uy_xl, uy_xr, uz_xl, uz_xr, 
  //                ux_yl, ux_yr, uy_yl, uy_yr, uz_yl, uz_yr, 
  //                ux_zl, ux_zr, uy_zl, uy_zr, uz_zl, uz_zr]  
  // prim_l/c/r:   Input volume expansion of primitive variables [1/rho div(p_par b), T_perp/m, m/T_perp, 3*Txx/m, 3*Tyy/m, 3*Tzz/m] in left/center/right cells.
  // pkpm_u_c:     Input volume expansion of flow velocity in center cell.
  // bvar_c:       Input volume expansion of magnetic field unit vector and tensor in center cell.
  // nu_c:         Input volume expansion of collisionality in center cell.
  // pkpm_lax:     Surface expansion of pkpm Lax penalization: lambda_i = |u_i| + sqrt(3.0*T_ii/m).
  // pkpm_accel:   Volume expansion of pkpm acceleration variables.

  const double dx1 = 2.0/dxv[2]; 
  const double *ux_c = &pkpm_u_c[0]; 
  const double *uy_c = &pkpm_u_c[8]; 
  const double *uz_c = &pkpm_u_c[16]; 

  const double *bxbx = &bvar_c[81]; 
  const double *bxby = &bvar_c[108]; 
  const double *bxbz = &bvar_c[135]; 
  const double *byby = &bvar_c[162]; 
  const double *bybz = &bvar_c[189]; 
  const double *bzbz = &bvar_c[216]; 

  const double *ux_surf_lr = &u_surf_l[52]; 
  const double *uy_surf_lr = &u_surf_l[60]; 
  const double *uz_surf_lr = &u_surf_l[68]; 

  const double *ux_surf_cl = &u_surf_c[48]; 
  const double *uy_surf_cl = &u_surf_c[56]; 
  const double *uz_surf_cl = &u_surf_c[64]; 

  const double *ux_surf_cr = &u_surf_c[52]; 
  const double *uy_surf_cr = &u_surf_c[60]; 
  const double *uz_surf_cr = &u_surf_c[68]; 

  const double *ux_surf_rl = &u_surf_r[48]; 
  const double *uy_surf_rl = &u_surf_r[56]; 
  const double *uz_surf_rl = &u_surf_r[64]; 

  const double *Tii_l = &prim_l[135]; 
  const double *Tii_c = &prim_c[135]; 
  const double *Tii_r = &prim_r[135]; 

  double *pkpm_lax_l = &pkpm_lax[36]; 
  double *pkpm_lax_r = &pkpm_lax[45]; 

  double *bb_grad_u = &pkpm_accel[27]; 
  double *p_perp_source = &pkpm_accel[81]; 

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
  double pkpm_lax_quad_l[9] = {0.0}; 
  double pkpm_lax_quad_r[9] = {0.0}; 

  ul_r = 0.9*uz_surf_lr[3]-0.6708203932499369*uz_surf_lr[2]-0.6708203932499369*uz_surf_lr[1]+0.5*uz_surf_lr[0]; 
  uc_l = 0.9*uz_surf_cl[3]-0.6708203932499369*uz_surf_cl[2]-0.6708203932499369*uz_surf_cl[1]+0.5*uz_surf_cl[0]; 
  uc_r = 0.9*uz_surf_cr[3]-0.6708203932499369*uz_surf_cr[2]-0.6708203932499369*uz_surf_cr[1]+0.5*uz_surf_cr[0]; 
  ur_l = 0.9*uz_surf_rl[3]-0.6708203932499369*uz_surf_rl[2]-0.6708203932499369*uz_surf_rl[1]+0.5*uz_surf_rl[0]; 
  uQuad_l = fmax(fabs(ul_r), fabs(uc_l)); 
  uQuad_r = fmax(fabs(uc_r), fabs(ur_l)); 
  Tiil_r = tensor_3x_p2_surfx3_eval_quad_node_0_r(Tii_l); 
  Tiic_l = tensor_3x_p2_surfx3_eval_quad_node_0_l(Tii_c); 
  Tiic_r = tensor_3x_p2_surfx3_eval_quad_node_0_r(Tii_c); 
  Tiir_l = tensor_3x_p2_surfx3_eval_quad_node_0_l(Tii_r); 
  TiiQuad_l = fmax(sqrt(fabs(Tiil_r)), sqrt(fabs(Tiic_l))); 
  TiiQuad_r = fmax(sqrt(fabs(Tiic_r)), sqrt(fabs(Tiir_l))); 
  pkpm_lax_quad_l[0] = uQuad_l + TiiQuad_l; 
  pkpm_lax_quad_r[0] = uQuad_r + TiiQuad_r; 

  ul_r = 0.5*uz_surf_lr[0]-0.6708203932499369*uz_surf_lr[1]; 
  uc_l = 0.5*uz_surf_cl[0]-0.6708203932499369*uz_surf_cl[1]; 
  uc_r = 0.5*uz_surf_cr[0]-0.6708203932499369*uz_surf_cr[1]; 
  ur_l = 0.5*uz_surf_rl[0]-0.6708203932499369*uz_surf_rl[1]; 
  uQuad_l = fmax(fabs(ul_r), fabs(uc_l)); 
  uQuad_r = fmax(fabs(uc_r), fabs(ur_l)); 
  Tiil_r = tensor_3x_p2_surfx3_eval_quad_node_1_r(Tii_l); 
  Tiic_l = tensor_3x_p2_surfx3_eval_quad_node_1_l(Tii_c); 
  Tiic_r = tensor_3x_p2_surfx3_eval_quad_node_1_r(Tii_c); 
  Tiir_l = tensor_3x_p2_surfx3_eval_quad_node_1_l(Tii_r); 
  TiiQuad_l = fmax(sqrt(fabs(Tiil_r)), sqrt(fabs(Tiic_l))); 
  TiiQuad_r = fmax(sqrt(fabs(Tiic_r)), sqrt(fabs(Tiir_l))); 
  pkpm_lax_quad_l[1] = uQuad_l + TiiQuad_l; 
  pkpm_lax_quad_r[1] = uQuad_r + TiiQuad_r; 

  ul_r = (-0.9*uz_surf_lr[3])+0.6708203932499369*uz_surf_lr[2]-0.6708203932499369*uz_surf_lr[1]+0.5*uz_surf_lr[0]; 
  uc_l = (-0.9*uz_surf_cl[3])+0.6708203932499369*uz_surf_cl[2]-0.6708203932499369*uz_surf_cl[1]+0.5*uz_surf_cl[0]; 
  uc_r = (-0.9*uz_surf_cr[3])+0.6708203932499369*uz_surf_cr[2]-0.6708203932499369*uz_surf_cr[1]+0.5*uz_surf_cr[0]; 
  ur_l = (-0.9*uz_surf_rl[3])+0.6708203932499369*uz_surf_rl[2]-0.6708203932499369*uz_surf_rl[1]+0.5*uz_surf_rl[0]; 
  uQuad_l = fmax(fabs(ul_r), fabs(uc_l)); 
  uQuad_r = fmax(fabs(uc_r), fabs(ur_l)); 
  Tiil_r = tensor_3x_p2_surfx3_eval_quad_node_2_r(Tii_l); 
  Tiic_l = tensor_3x_p2_surfx3_eval_quad_node_2_l(Tii_c); 
  Tiic_r = tensor_3x_p2_surfx3_eval_quad_node_2_r(Tii_c); 
  Tiir_l = tensor_3x_p2_surfx3_eval_quad_node_2_l(Tii_r); 
  TiiQuad_l = fmax(sqrt(fabs(Tiil_r)), sqrt(fabs(Tiic_l))); 
  TiiQuad_r = fmax(sqrt(fabs(Tiic_r)), sqrt(fabs(Tiir_l))); 
  pkpm_lax_quad_l[2] = uQuad_l + TiiQuad_l; 
  pkpm_lax_quad_r[2] = uQuad_r + TiiQuad_r; 

  ul_r = 0.5*uz_surf_lr[0]-0.6708203932499369*uz_surf_lr[2]; 
  uc_l = 0.5*uz_surf_cl[0]-0.6708203932499369*uz_surf_cl[2]; 
  uc_r = 0.5*uz_surf_cr[0]-0.6708203932499369*uz_surf_cr[2]; 
  ur_l = 0.5*uz_surf_rl[0]-0.6708203932499369*uz_surf_rl[2]; 
  uQuad_l = fmax(fabs(ul_r), fabs(uc_l)); 
  uQuad_r = fmax(fabs(uc_r), fabs(ur_l)); 
  Tiil_r = tensor_3x_p2_surfx3_eval_quad_node_3_r(Tii_l); 
  Tiic_l = tensor_3x_p2_surfx3_eval_quad_node_3_l(Tii_c); 
  Tiic_r = tensor_3x_p2_surfx3_eval_quad_node_3_r(Tii_c); 
  Tiir_l = tensor_3x_p2_surfx3_eval_quad_node_3_l(Tii_r); 
  TiiQuad_l = fmax(sqrt(fabs(Tiil_r)), sqrt(fabs(Tiic_l))); 
  TiiQuad_r = fmax(sqrt(fabs(Tiic_r)), sqrt(fabs(Tiir_l))); 
  pkpm_lax_quad_l[3] = uQuad_l + TiiQuad_l; 
  pkpm_lax_quad_r[3] = uQuad_r + TiiQuad_r; 

  ul_r = 0.5*uz_surf_lr[0]; 
  uc_l = 0.5*uz_surf_cl[0]; 
  uc_r = 0.5*uz_surf_cr[0]; 
  ur_l = 0.5*uz_surf_rl[0]; 
  uQuad_l = fmax(fabs(ul_r), fabs(uc_l)); 
  uQuad_r = fmax(fabs(uc_r), fabs(ur_l)); 
  Tiil_r = tensor_3x_p2_surfx3_eval_quad_node_4_r(Tii_l); 
  Tiic_l = tensor_3x_p2_surfx3_eval_quad_node_4_l(Tii_c); 
  Tiic_r = tensor_3x_p2_surfx3_eval_quad_node_4_r(Tii_c); 
  Tiir_l = tensor_3x_p2_surfx3_eval_quad_node_4_l(Tii_r); 
  TiiQuad_l = fmax(sqrt(fabs(Tiil_r)), sqrt(fabs(Tiic_l))); 
  TiiQuad_r = fmax(sqrt(fabs(Tiic_r)), sqrt(fabs(Tiir_l))); 
  pkpm_lax_quad_l[4] = uQuad_l + TiiQuad_l; 
  pkpm_lax_quad_r[4] = uQuad_r + TiiQuad_r; 

  ul_r = 0.6708203932499369*uz_surf_lr[2]+0.5*uz_surf_lr[0]; 
  uc_l = 0.6708203932499369*uz_surf_cl[2]+0.5*uz_surf_cl[0]; 
  uc_r = 0.6708203932499369*uz_surf_cr[2]+0.5*uz_surf_cr[0]; 
  ur_l = 0.6708203932499369*uz_surf_rl[2]+0.5*uz_surf_rl[0]; 
  uQuad_l = fmax(fabs(ul_r), fabs(uc_l)); 
  uQuad_r = fmax(fabs(uc_r), fabs(ur_l)); 
  Tiil_r = tensor_3x_p2_surfx3_eval_quad_node_5_r(Tii_l); 
  Tiic_l = tensor_3x_p2_surfx3_eval_quad_node_5_l(Tii_c); 
  Tiic_r = tensor_3x_p2_surfx3_eval_quad_node_5_r(Tii_c); 
  Tiir_l = tensor_3x_p2_surfx3_eval_quad_node_5_l(Tii_r); 
  TiiQuad_l = fmax(sqrt(fabs(Tiil_r)), sqrt(fabs(Tiic_l))); 
  TiiQuad_r = fmax(sqrt(fabs(Tiic_r)), sqrt(fabs(Tiir_l))); 
  pkpm_lax_quad_l[5] = uQuad_l + TiiQuad_l; 
  pkpm_lax_quad_r[5] = uQuad_r + TiiQuad_r; 

  ul_r = (-0.9*uz_surf_lr[3])-0.6708203932499369*uz_surf_lr[2]+0.6708203932499369*uz_surf_lr[1]+0.5*uz_surf_lr[0]; 
  uc_l = (-0.9*uz_surf_cl[3])-0.6708203932499369*uz_surf_cl[2]+0.6708203932499369*uz_surf_cl[1]+0.5*uz_surf_cl[0]; 
  uc_r = (-0.9*uz_surf_cr[3])-0.6708203932499369*uz_surf_cr[2]+0.6708203932499369*uz_surf_cr[1]+0.5*uz_surf_cr[0]; 
  ur_l = (-0.9*uz_surf_rl[3])-0.6708203932499369*uz_surf_rl[2]+0.6708203932499369*uz_surf_rl[1]+0.5*uz_surf_rl[0]; 
  uQuad_l = fmax(fabs(ul_r), fabs(uc_l)); 
  uQuad_r = fmax(fabs(uc_r), fabs(ur_l)); 
  Tiil_r = tensor_3x_p2_surfx3_eval_quad_node_6_r(Tii_l); 
  Tiic_l = tensor_3x_p2_surfx3_eval_quad_node_6_l(Tii_c); 
  Tiic_r = tensor_3x_p2_surfx3_eval_quad_node_6_r(Tii_c); 
  Tiir_l = tensor_3x_p2_surfx3_eval_quad_node_6_l(Tii_r); 
  TiiQuad_l = fmax(sqrt(fabs(Tiil_r)), sqrt(fabs(Tiic_l))); 
  TiiQuad_r = fmax(sqrt(fabs(Tiic_r)), sqrt(fabs(Tiir_l))); 
  pkpm_lax_quad_l[6] = uQuad_l + TiiQuad_l; 
  pkpm_lax_quad_r[6] = uQuad_r + TiiQuad_r; 

  ul_r = 0.6708203932499369*uz_surf_lr[1]+0.5*uz_surf_lr[0]; 
  uc_l = 0.6708203932499369*uz_surf_cl[1]+0.5*uz_surf_cl[0]; 
  uc_r = 0.6708203932499369*uz_surf_cr[1]+0.5*uz_surf_cr[0]; 
  ur_l = 0.6708203932499369*uz_surf_rl[1]+0.5*uz_surf_rl[0]; 
  uQuad_l = fmax(fabs(ul_r), fabs(uc_l)); 
  uQuad_r = fmax(fabs(uc_r), fabs(ur_l)); 
  Tiil_r = tensor_3x_p2_surfx3_eval_quad_node_7_r(Tii_l); 
  Tiic_l = tensor_3x_p2_surfx3_eval_quad_node_7_l(Tii_c); 
  Tiic_r = tensor_3x_p2_surfx3_eval_quad_node_7_r(Tii_c); 
  Tiir_l = tensor_3x_p2_surfx3_eval_quad_node_7_l(Tii_r); 
  TiiQuad_l = fmax(sqrt(fabs(Tiil_r)), sqrt(fabs(Tiic_l))); 
  TiiQuad_r = fmax(sqrt(fabs(Tiic_r)), sqrt(fabs(Tiir_l))); 
  pkpm_lax_quad_l[7] = uQuad_l + TiiQuad_l; 
  pkpm_lax_quad_r[7] = uQuad_r + TiiQuad_r; 

  ul_r = 0.9*uz_surf_lr[3]+0.6708203932499369*uz_surf_lr[2]+0.6708203932499369*uz_surf_lr[1]+0.5*uz_surf_lr[0]; 
  uc_l = 0.9*uz_surf_cl[3]+0.6708203932499369*uz_surf_cl[2]+0.6708203932499369*uz_surf_cl[1]+0.5*uz_surf_cl[0]; 
  uc_r = 0.9*uz_surf_cr[3]+0.6708203932499369*uz_surf_cr[2]+0.6708203932499369*uz_surf_cr[1]+0.5*uz_surf_cr[0]; 
  ur_l = 0.9*uz_surf_rl[3]+0.6708203932499369*uz_surf_rl[2]+0.6708203932499369*uz_surf_rl[1]+0.5*uz_surf_rl[0]; 
  uQuad_l = fmax(fabs(ul_r), fabs(uc_l)); 
  uQuad_r = fmax(fabs(uc_r), fabs(ur_l)); 
  Tiil_r = tensor_3x_p2_surfx3_eval_quad_node_8_r(Tii_l); 
  Tiic_l = tensor_3x_p2_surfx3_eval_quad_node_8_l(Tii_c); 
  Tiic_r = tensor_3x_p2_surfx3_eval_quad_node_8_r(Tii_c); 
  Tiir_l = tensor_3x_p2_surfx3_eval_quad_node_8_l(Tii_r); 
  TiiQuad_l = fmax(sqrt(fabs(Tiil_r)), sqrt(fabs(Tiic_l))); 
  TiiQuad_r = fmax(sqrt(fabs(Tiic_r)), sqrt(fabs(Tiir_l))); 
  pkpm_lax_quad_l[8] = uQuad_l + TiiQuad_l; 
  pkpm_lax_quad_r[8] = uQuad_r + TiiQuad_r; 

  tensor_3x_p2_upwind_quad_to_modal(pkpm_lax_quad_l, pkpm_lax_l); 
  tensor_3x_p2_upwind_quad_to_modal(pkpm_lax_quad_r, pkpm_lax_r); 

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

  double bb_grad_u_comp[27] = {0.0}; 
  bb_grad_u_comp[0] = 0.3535533905932737*grad_u_z[7]*bzbz[10]+0.3535533905932737*grad_u_y[7]*bybz[10]+0.3535533905932737*grad_u_x[7]*bxbz[10]+0.3535533905932737*bzbz[6]*grad_u_z[6]+0.3535533905932737*bybz[6]*grad_u_y[6]+0.3535533905932737*bxbz[6]*grad_u_x[6]+0.3535533905932737*bzbz[5]*grad_u_z[5]+0.3535533905932737*bybz[5]*grad_u_y[5]+0.3535533905932737*bxbz[5]*grad_u_x[5]+0.3535533905932737*bzbz[4]*grad_u_z[4]+0.3535533905932737*bybz[4]*grad_u_y[4]+0.3535533905932737*bxbz[4]*grad_u_x[4]+0.3535533905932737*bzbz[3]*grad_u_z[3]+0.3535533905932737*bybz[3]*grad_u_y[3]+0.3535533905932737*bxbz[3]*grad_u_x[3]+0.3535533905932737*bzbz[2]*grad_u_z[2]+0.3535533905932737*bybz[2]*grad_u_y[2]+0.3535533905932737*bxbz[2]*grad_u_x[2]+0.3535533905932737*bzbz[1]*grad_u_z[1]+0.3535533905932737*bybz[1]*grad_u_y[1]+0.3535533905932737*bxbz[1]*grad_u_x[1]+0.3535533905932737*bzbz[0]*grad_u_z[0]+0.3535533905932737*bybz[0]*grad_u_y[0]+0.3535533905932737*bxbz[0]*grad_u_x[0]; 
  bb_grad_u_comp[1] = 0.3162277660168379*grad_u_z[7]*bzbz[17]+0.3162277660168379*grad_u_y[7]*bybz[17]+0.3162277660168379*grad_u_x[7]*bxbz[17]+0.3162277660168379*grad_u_z[5]*bzbz[13]+0.3162277660168379*grad_u_y[5]*bybz[13]+0.3162277660168379*grad_u_x[5]*bxbz[13]+0.3162277660168379*grad_u_z[4]*bzbz[11]+0.3162277660168379*grad_u_y[4]*bybz[11]+0.3162277660168379*grad_u_x[4]*bxbz[11]+0.3535533905932737*grad_u_z[6]*bzbz[10]+0.3535533905932737*grad_u_y[6]*bybz[10]+0.3535533905932737*grad_u_x[6]*bxbz[10]+0.3535533905932737*bzbz[6]*grad_u_z[7]+0.3535533905932737*bybz[6]*grad_u_y[7]+0.3535533905932737*bxbz[6]*grad_u_x[7]+0.3162277660168379*grad_u_z[1]*bzbz[7]+0.3162277660168379*grad_u_y[1]*bybz[7]+0.3162277660168379*grad_u_x[1]*bxbz[7]+0.3535533905932737*bzbz[3]*grad_u_z[5]+0.3535533905932737*bybz[3]*grad_u_y[5]+0.3535533905932737*bxbz[3]*grad_u_x[5]+0.3535533905932737*grad_u_z[3]*bzbz[5]+0.3535533905932737*grad_u_y[3]*bybz[5]+0.3535533905932737*grad_u_x[3]*bxbz[5]+0.3535533905932737*bzbz[2]*grad_u_z[4]+0.3535533905932737*bybz[2]*grad_u_y[4]+0.3535533905932737*bxbz[2]*grad_u_x[4]+0.3535533905932737*grad_u_z[2]*bzbz[4]+0.3535533905932737*grad_u_y[2]*bybz[4]+0.3535533905932737*grad_u_x[2]*bxbz[4]+0.3535533905932737*bzbz[0]*grad_u_z[1]+0.3535533905932737*bybz[0]*grad_u_y[1]+0.3535533905932737*bxbz[0]*grad_u_x[1]+0.3535533905932737*grad_u_z[0]*bzbz[1]+0.3535533905932737*grad_u_y[0]*bybz[1]+0.3535533905932737*grad_u_x[0]*bxbz[1]; 
  bb_grad_u_comp[2] = 0.3162277660168379*grad_u_z[7]*bzbz[18]+0.3162277660168379*grad_u_y[7]*bybz[18]+0.3162277660168379*grad_u_x[7]*bxbz[18]+0.3162277660168379*grad_u_z[6]*bzbz[14]+0.3162277660168379*grad_u_y[6]*bybz[14]+0.3162277660168379*grad_u_x[6]*bxbz[14]+0.3162277660168379*grad_u_z[4]*bzbz[12]+0.3162277660168379*grad_u_y[4]*bybz[12]+0.3162277660168379*grad_u_x[4]*bxbz[12]+0.3535533905932737*grad_u_z[5]*bzbz[10]+0.3535533905932737*grad_u_y[5]*bybz[10]+0.3535533905932737*grad_u_x[5]*bxbz[10]+0.3162277660168379*grad_u_z[2]*bzbz[8]+0.3162277660168379*grad_u_y[2]*bybz[8]+0.3162277660168379*grad_u_x[2]*bxbz[8]+0.3535533905932737*bzbz[5]*grad_u_z[7]+0.3535533905932737*bybz[5]*grad_u_y[7]+0.3535533905932737*bxbz[5]*grad_u_x[7]+0.3535533905932737*bzbz[3]*grad_u_z[6]+0.3535533905932737*bybz[3]*grad_u_y[6]+0.3535533905932737*bxbz[3]*grad_u_x[6]+0.3535533905932737*grad_u_z[3]*bzbz[6]+0.3535533905932737*grad_u_y[3]*bybz[6]+0.3535533905932737*grad_u_x[3]*bxbz[6]+0.3535533905932737*bzbz[1]*grad_u_z[4]+0.3535533905932737*bybz[1]*grad_u_y[4]+0.3535533905932737*bxbz[1]*grad_u_x[4]+0.3535533905932737*grad_u_z[1]*bzbz[4]+0.3535533905932737*grad_u_y[1]*bybz[4]+0.3535533905932737*grad_u_x[1]*bxbz[4]+0.3535533905932737*bzbz[0]*grad_u_z[2]+0.3535533905932737*bybz[0]*grad_u_y[2]+0.3535533905932737*bxbz[0]*grad_u_x[2]+0.3535533905932737*grad_u_z[0]*bzbz[2]+0.3535533905932737*grad_u_y[0]*bybz[2]+0.3535533905932737*grad_u_x[0]*bxbz[2]; 
  bb_grad_u_comp[3] = 0.3162277660168379*grad_u_z[7]*bzbz[19]+0.3162277660168379*grad_u_y[7]*bybz[19]+0.3162277660168379*grad_u_x[7]*bxbz[19]+0.3162277660168379*grad_u_z[6]*bzbz[16]+0.3162277660168379*grad_u_y[6]*bybz[16]+0.3162277660168379*grad_u_x[6]*bxbz[16]+0.3162277660168379*grad_u_z[5]*bzbz[15]+0.3162277660168379*grad_u_y[5]*bybz[15]+0.3162277660168379*grad_u_x[5]*bxbz[15]+0.3535533905932737*grad_u_z[4]*bzbz[10]+0.3535533905932737*grad_u_y[4]*bybz[10]+0.3535533905932737*grad_u_x[4]*bxbz[10]+0.3162277660168379*grad_u_z[3]*bzbz[9]+0.3162277660168379*grad_u_y[3]*bybz[9]+0.3162277660168379*grad_u_x[3]*bxbz[9]+0.3535533905932737*bzbz[4]*grad_u_z[7]+0.3535533905932737*bybz[4]*grad_u_y[7]+0.3535533905932737*bxbz[4]*grad_u_x[7]+0.3535533905932737*bzbz[2]*grad_u_z[6]+0.3535533905932737*bybz[2]*grad_u_y[6]+0.3535533905932737*bxbz[2]*grad_u_x[6]+0.3535533905932737*grad_u_z[2]*bzbz[6]+0.3535533905932737*grad_u_y[2]*bybz[6]+0.3535533905932737*grad_u_x[2]*bxbz[6]+0.3535533905932737*bzbz[1]*grad_u_z[5]+0.3535533905932737*bybz[1]*grad_u_y[5]+0.3535533905932737*bxbz[1]*grad_u_x[5]+0.3535533905932737*grad_u_z[1]*bzbz[5]+0.3535533905932737*grad_u_y[1]*bybz[5]+0.3535533905932737*grad_u_x[1]*bxbz[5]+0.3535533905932737*bzbz[0]*grad_u_z[3]+0.3535533905932737*bybz[0]*grad_u_y[3]+0.3535533905932737*bxbz[0]*grad_u_x[3]+0.3535533905932737*grad_u_z[0]*bzbz[3]+0.3535533905932737*grad_u_y[0]*bybz[3]+0.3535533905932737*grad_u_x[0]*bxbz[3]; 
  bb_grad_u_comp[4] = 0.2828427124746191*grad_u_z[7]*bzbz[23]+0.2828427124746191*grad_u_y[7]*bybz[23]+0.2828427124746191*grad_u_x[7]*bxbz[23]+0.2828427124746191*grad_u_z[4]*bzbz[20]+0.2828427124746191*grad_u_y[4]*bybz[20]+0.2828427124746191*grad_u_x[4]*bxbz[20]+0.3162277660168379*grad_u_z[6]*bzbz[18]+0.3162277660168379*grad_u_y[6]*bybz[18]+0.3162277660168379*grad_u_x[6]*bxbz[18]+0.3162277660168379*grad_u_z[5]*bzbz[17]+0.3162277660168379*grad_u_y[5]*bybz[17]+0.3162277660168379*grad_u_x[5]*bxbz[17]+0.3162277660168379*grad_u_z[7]*bzbz[14]+0.3162277660168379*grad_u_y[7]*bybz[14]+0.3162277660168379*grad_u_x[7]*bxbz[14]+0.3162277660168379*grad_u_z[7]*bzbz[13]+0.3162277660168379*grad_u_y[7]*bybz[13]+0.3162277660168379*grad_u_x[7]*bxbz[13]+0.3162277660168379*grad_u_z[2]*bzbz[12]+0.3162277660168379*grad_u_y[2]*bybz[12]+0.3162277660168379*grad_u_x[2]*bxbz[12]+0.3162277660168379*grad_u_z[1]*bzbz[11]+0.3162277660168379*grad_u_y[1]*bybz[11]+0.3162277660168379*grad_u_x[1]*bxbz[11]+0.3535533905932737*grad_u_z[3]*bzbz[10]+0.3535533905932737*grad_u_y[3]*bybz[10]+0.3535533905932737*grad_u_x[3]*bxbz[10]+0.3162277660168379*grad_u_z[4]*bzbz[8]+0.3162277660168379*grad_u_y[4]*bybz[8]+0.3162277660168379*grad_u_x[4]*bxbz[8]+0.3535533905932737*bzbz[3]*grad_u_z[7]+0.3535533905932737*bybz[3]*grad_u_y[7]+0.3535533905932737*bxbz[3]*grad_u_x[7]+0.3162277660168379*grad_u_z[4]*bzbz[7]+0.3162277660168379*grad_u_y[4]*bybz[7]+0.3162277660168379*grad_u_x[4]*bxbz[7]+0.3535533905932737*bzbz[5]*grad_u_z[6]+0.3535533905932737*bybz[5]*grad_u_y[6]+0.3535533905932737*bxbz[5]*grad_u_x[6]+0.3535533905932737*grad_u_z[5]*bzbz[6]+0.3535533905932737*grad_u_y[5]*bybz[6]+0.3535533905932737*grad_u_x[5]*bxbz[6]+0.3535533905932737*bzbz[0]*grad_u_z[4]+0.3535533905932737*bybz[0]*grad_u_y[4]+0.3535533905932737*bxbz[0]*grad_u_x[4]+0.3535533905932737*grad_u_z[0]*bzbz[4]+0.3535533905932737*grad_u_y[0]*bybz[4]+0.3535533905932737*grad_u_x[0]*bxbz[4]+0.3535533905932737*bzbz[1]*grad_u_z[2]+0.3535533905932737*bybz[1]*grad_u_y[2]+0.3535533905932737*bxbz[1]*grad_u_x[2]+0.3535533905932737*grad_u_z[1]*bzbz[2]+0.3535533905932737*grad_u_y[1]*bybz[2]+0.3535533905932737*grad_u_x[1]*bxbz[2]; 
  bb_grad_u_comp[5] = 0.2828427124746191*grad_u_z[7]*bzbz[24]+0.2828427124746191*grad_u_y[7]*bybz[24]+0.2828427124746191*grad_u_x[7]*bxbz[24]+0.2828427124746191*grad_u_z[5]*bzbz[21]+0.2828427124746191*grad_u_y[5]*bybz[21]+0.2828427124746191*grad_u_x[5]*bxbz[21]+0.3162277660168379*grad_u_z[6]*bzbz[19]+0.3162277660168379*grad_u_y[6]*bybz[19]+0.3162277660168379*grad_u_x[6]*bxbz[19]+0.3162277660168379*grad_u_z[4]*bzbz[17]+0.3162277660168379*grad_u_y[4]*bybz[17]+0.3162277660168379*grad_u_x[4]*bxbz[17]+0.3162277660168379*grad_u_z[7]*bzbz[16]+0.3162277660168379*grad_u_y[7]*bybz[16]+0.3162277660168379*grad_u_x[7]*bxbz[16]+0.3162277660168379*grad_u_z[3]*bzbz[15]+0.3162277660168379*grad_u_y[3]*bybz[15]+0.3162277660168379*grad_u_x[3]*bxbz[15]+0.3162277660168379*grad_u_z[1]*bzbz[13]+0.3162277660168379*grad_u_y[1]*bybz[13]+0.3162277660168379*grad_u_x[1]*bxbz[13]+0.3162277660168379*grad_u_z[7]*bzbz[11]+0.3162277660168379*grad_u_y[7]*bybz[11]+0.3162277660168379*grad_u_x[7]*bxbz[11]+0.3535533905932737*grad_u_z[2]*bzbz[10]+0.3535533905932737*grad_u_y[2]*bybz[10]+0.3535533905932737*grad_u_x[2]*bxbz[10]+0.3162277660168379*grad_u_z[5]*bzbz[9]+0.3162277660168379*grad_u_y[5]*bybz[9]+0.3162277660168379*grad_u_x[5]*bxbz[9]+0.3535533905932737*bzbz[2]*grad_u_z[7]+0.3535533905932737*bybz[2]*grad_u_y[7]+0.3535533905932737*bxbz[2]*grad_u_x[7]+0.3162277660168379*grad_u_z[5]*bzbz[7]+0.3162277660168379*grad_u_y[5]*bybz[7]+0.3162277660168379*grad_u_x[5]*bxbz[7]+0.3535533905932737*bzbz[4]*grad_u_z[6]+0.3535533905932737*bybz[4]*grad_u_y[6]+0.3535533905932737*bxbz[4]*grad_u_x[6]+0.3535533905932737*grad_u_z[4]*bzbz[6]+0.3535533905932737*grad_u_y[4]*bybz[6]+0.3535533905932737*grad_u_x[4]*bxbz[6]+0.3535533905932737*bzbz[0]*grad_u_z[5]+0.3535533905932737*bybz[0]*grad_u_y[5]+0.3535533905932737*bxbz[0]*grad_u_x[5]+0.3535533905932737*grad_u_z[0]*bzbz[5]+0.3535533905932737*grad_u_y[0]*bybz[5]+0.3535533905932737*grad_u_x[0]*bxbz[5]+0.3535533905932737*bzbz[1]*grad_u_z[3]+0.3535533905932737*bybz[1]*grad_u_y[3]+0.3535533905932737*bxbz[1]*grad_u_x[3]+0.3535533905932737*grad_u_z[1]*bzbz[3]+0.3535533905932737*grad_u_y[1]*bybz[3]+0.3535533905932737*grad_u_x[1]*bxbz[3]; 
  bb_grad_u_comp[6] = 0.2828427124746191*grad_u_z[7]*bzbz[25]+0.2828427124746191*grad_u_y[7]*bybz[25]+0.2828427124746191*grad_u_x[7]*bxbz[25]+0.2828427124746191*grad_u_z[6]*bzbz[22]+0.2828427124746191*grad_u_y[6]*bybz[22]+0.2828427124746191*grad_u_x[6]*bxbz[22]+0.3162277660168379*grad_u_z[5]*bzbz[19]+0.3162277660168379*grad_u_y[5]*bybz[19]+0.3162277660168379*grad_u_x[5]*bxbz[19]+0.3162277660168379*grad_u_z[4]*bzbz[18]+0.3162277660168379*grad_u_y[4]*bybz[18]+0.3162277660168379*grad_u_x[4]*bxbz[18]+0.3162277660168379*grad_u_z[3]*bzbz[16]+0.3162277660168379*grad_u_y[3]*bybz[16]+0.3162277660168379*grad_u_x[3]*bxbz[16]+0.3162277660168379*grad_u_z[7]*bzbz[15]+0.3162277660168379*grad_u_y[7]*bybz[15]+0.3162277660168379*grad_u_x[7]*bxbz[15]+0.3162277660168379*grad_u_z[2]*bzbz[14]+0.3162277660168379*grad_u_y[2]*bybz[14]+0.3162277660168379*grad_u_x[2]*bxbz[14]+0.3162277660168379*grad_u_z[7]*bzbz[12]+0.3162277660168379*grad_u_y[7]*bybz[12]+0.3162277660168379*grad_u_x[7]*bxbz[12]+0.3535533905932737*grad_u_z[1]*bzbz[10]+0.3535533905932737*grad_u_y[1]*bybz[10]+0.3535533905932737*grad_u_x[1]*bxbz[10]+0.3162277660168379*grad_u_z[6]*bzbz[9]+0.3162277660168379*grad_u_y[6]*bybz[9]+0.3162277660168379*grad_u_x[6]*bxbz[9]+0.3162277660168379*grad_u_z[6]*bzbz[8]+0.3162277660168379*grad_u_y[6]*bybz[8]+0.3162277660168379*grad_u_x[6]*bxbz[8]+0.3535533905932737*bzbz[1]*grad_u_z[7]+0.3535533905932737*bybz[1]*grad_u_y[7]+0.3535533905932737*bxbz[1]*grad_u_x[7]+0.3535533905932737*bzbz[0]*grad_u_z[6]+0.3535533905932737*bybz[0]*grad_u_y[6]+0.3535533905932737*bxbz[0]*grad_u_x[6]+0.3535533905932737*grad_u_z[0]*bzbz[6]+0.3535533905932737*grad_u_y[0]*bybz[6]+0.3535533905932737*grad_u_x[0]*bxbz[6]+0.3535533905932737*bzbz[4]*grad_u_z[5]+0.3535533905932737*bybz[4]*grad_u_y[5]+0.3535533905932737*bxbz[4]*grad_u_x[5]+0.3535533905932737*grad_u_z[4]*bzbz[5]+0.3535533905932737*grad_u_y[4]*bybz[5]+0.3535533905932737*grad_u_x[4]*bxbz[5]+0.3535533905932737*bzbz[2]*grad_u_z[3]+0.3535533905932737*bybz[2]*grad_u_y[3]+0.3535533905932737*bxbz[2]*grad_u_x[3]+0.3535533905932737*grad_u_z[2]*bzbz[3]+0.3535533905932737*grad_u_y[2]*bybz[3]+0.3535533905932737*grad_u_x[2]*bxbz[3]; 
  bb_grad_u_comp[7] = 0.3535533905932737*grad_u_z[6]*bzbz[17]+0.3535533905932737*grad_u_y[6]*bybz[17]+0.3535533905932737*grad_u_x[6]*bxbz[17]+0.3535533905932737*grad_u_z[3]*bzbz[13]+0.3535533905932737*grad_u_y[3]*bybz[13]+0.3535533905932737*grad_u_x[3]*bxbz[13]+0.3535533905932737*grad_u_z[2]*bzbz[11]+0.3535533905932737*grad_u_y[2]*bybz[11]+0.3535533905932737*grad_u_x[2]*bxbz[11]+0.3162277660168379*grad_u_z[7]*bzbz[10]+0.3162277660168379*grad_u_y[7]*bybz[10]+0.3162277660168379*grad_u_x[7]*bxbz[10]+0.3535533905932737*grad_u_z[0]*bzbz[7]+0.3535533905932737*grad_u_y[0]*bybz[7]+0.3535533905932737*grad_u_x[0]*bxbz[7]+0.3162277660168379*bzbz[5]*grad_u_z[5]+0.3162277660168379*bybz[5]*grad_u_y[5]+0.3162277660168379*bxbz[5]*grad_u_x[5]+0.3162277660168379*bzbz[4]*grad_u_z[4]+0.3162277660168379*bybz[4]*grad_u_y[4]+0.3162277660168379*bxbz[4]*grad_u_x[4]+0.3162277660168379*bzbz[1]*grad_u_z[1]+0.3162277660168379*bybz[1]*grad_u_y[1]+0.3162277660168379*bxbz[1]*grad_u_x[1]; 
  bb_grad_u_comp[8] = 0.3535533905932737*grad_u_z[5]*bzbz[18]+0.3535533905932737*grad_u_y[5]*bybz[18]+0.3535533905932737*grad_u_x[5]*bxbz[18]+0.3535533905932737*grad_u_z[3]*bzbz[14]+0.3535533905932737*grad_u_y[3]*bybz[14]+0.3535533905932737*grad_u_x[3]*bxbz[14]+0.3535533905932737*grad_u_z[1]*bzbz[12]+0.3535533905932737*grad_u_y[1]*bybz[12]+0.3535533905932737*grad_u_x[1]*bxbz[12]+0.3162277660168379*grad_u_z[7]*bzbz[10]+0.3162277660168379*grad_u_y[7]*bybz[10]+0.3162277660168379*grad_u_x[7]*bxbz[10]+0.3535533905932737*grad_u_z[0]*bzbz[8]+0.3535533905932737*grad_u_y[0]*bybz[8]+0.3535533905932737*grad_u_x[0]*bxbz[8]+0.3162277660168379*bzbz[6]*grad_u_z[6]+0.3162277660168379*bybz[6]*grad_u_y[6]+0.3162277660168379*bxbz[6]*grad_u_x[6]+0.3162277660168379*bzbz[4]*grad_u_z[4]+0.3162277660168379*bybz[4]*grad_u_y[4]+0.3162277660168379*bxbz[4]*grad_u_x[4]+0.3162277660168379*bzbz[2]*grad_u_z[2]+0.3162277660168379*bybz[2]*grad_u_y[2]+0.3162277660168379*bxbz[2]*grad_u_x[2]; 
  bb_grad_u_comp[9] = 0.3535533905932737*grad_u_z[4]*bzbz[19]+0.3535533905932737*grad_u_y[4]*bybz[19]+0.3535533905932737*grad_u_x[4]*bxbz[19]+0.3535533905932737*grad_u_z[2]*bzbz[16]+0.3535533905932737*grad_u_y[2]*bybz[16]+0.3535533905932737*grad_u_x[2]*bxbz[16]+0.3535533905932737*grad_u_z[1]*bzbz[15]+0.3535533905932737*grad_u_y[1]*bybz[15]+0.3535533905932737*grad_u_x[1]*bxbz[15]+0.3162277660168379*grad_u_z[7]*bzbz[10]+0.3162277660168379*grad_u_y[7]*bybz[10]+0.3162277660168379*grad_u_x[7]*bxbz[10]+0.3535533905932737*grad_u_z[0]*bzbz[9]+0.3535533905932737*grad_u_y[0]*bybz[9]+0.3535533905932737*grad_u_x[0]*bxbz[9]+0.3162277660168379*bzbz[6]*grad_u_z[6]+0.3162277660168379*bybz[6]*grad_u_y[6]+0.3162277660168379*bxbz[6]*grad_u_x[6]+0.3162277660168379*bzbz[5]*grad_u_z[5]+0.3162277660168379*bybz[5]*grad_u_y[5]+0.3162277660168379*bxbz[5]*grad_u_x[5]+0.3162277660168379*bzbz[3]*grad_u_z[3]+0.3162277660168379*bybz[3]*grad_u_y[3]+0.3162277660168379*bxbz[3]*grad_u_x[3]; 
  bb_grad_u_comp[10] = 0.2529822128134704*grad_u_z[7]*bzbz[26]+0.2529822128134704*grad_u_y[7]*bybz[26]+0.2529822128134704*grad_u_x[7]*bxbz[26]+0.2828427124746191*grad_u_z[6]*bzbz[25]+0.2828427124746191*grad_u_y[6]*bybz[25]+0.2828427124746191*grad_u_x[6]*bxbz[25]+0.2828427124746191*grad_u_z[5]*bzbz[24]+0.2828427124746191*grad_u_y[5]*bybz[24]+0.2828427124746191*grad_u_x[5]*bxbz[24]+0.2828427124746191*grad_u_z[4]*bzbz[23]+0.2828427124746191*grad_u_y[4]*bybz[23]+0.2828427124746191*grad_u_x[4]*bxbz[23]+0.2828427124746191*grad_u_z[7]*bzbz[22]+0.2828427124746191*grad_u_y[7]*bybz[22]+0.2828427124746191*grad_u_x[7]*bxbz[22]+0.2828427124746191*grad_u_z[7]*bzbz[21]+0.2828427124746191*grad_u_y[7]*bybz[21]+0.2828427124746191*grad_u_x[7]*bxbz[21]+0.2828427124746191*grad_u_z[7]*bzbz[20]+0.2828427124746191*grad_u_y[7]*bybz[20]+0.2828427124746191*grad_u_x[7]*bxbz[20]+0.3162277660168379*grad_u_z[3]*bzbz[19]+0.3162277660168379*grad_u_y[3]*bybz[19]+0.3162277660168379*grad_u_x[3]*bxbz[19]+0.3162277660168379*grad_u_z[2]*bzbz[18]+0.3162277660168379*grad_u_y[2]*bybz[18]+0.3162277660168379*grad_u_x[2]*bxbz[18]+0.3162277660168379*grad_u_z[1]*bzbz[17]+0.3162277660168379*grad_u_y[1]*bybz[17]+0.3162277660168379*grad_u_x[1]*bxbz[17]+0.3162277660168379*grad_u_z[5]*bzbz[16]+0.3162277660168379*grad_u_y[5]*bybz[16]+0.3162277660168379*grad_u_x[5]*bxbz[16]+0.3162277660168379*grad_u_z[6]*bzbz[15]+0.3162277660168379*grad_u_y[6]*bybz[15]+0.3162277660168379*grad_u_x[6]*bxbz[15]+0.3162277660168379*grad_u_z[4]*bzbz[14]+0.3162277660168379*grad_u_y[4]*bybz[14]+0.3162277660168379*grad_u_x[4]*bxbz[14]+0.3162277660168379*grad_u_z[4]*bzbz[13]+0.3162277660168379*grad_u_y[4]*bybz[13]+0.3162277660168379*grad_u_x[4]*bxbz[13]+0.3162277660168379*grad_u_z[6]*bzbz[12]+0.3162277660168379*grad_u_y[6]*bybz[12]+0.3162277660168379*grad_u_x[6]*bxbz[12]+0.3162277660168379*grad_u_z[5]*bzbz[11]+0.3162277660168379*grad_u_y[5]*bybz[11]+0.3162277660168379*grad_u_x[5]*bxbz[11]+0.3535533905932737*grad_u_z[0]*bzbz[10]+0.3535533905932737*grad_u_y[0]*bybz[10]+0.3535533905932737*grad_u_x[0]*bxbz[10]+0.3162277660168379*grad_u_z[7]*bzbz[9]+0.3162277660168379*grad_u_y[7]*bybz[9]+0.3162277660168379*grad_u_x[7]*bxbz[9]+0.3162277660168379*grad_u_z[7]*bzbz[8]+0.3162277660168379*grad_u_y[7]*bybz[8]+0.3162277660168379*grad_u_x[7]*bxbz[8]+0.3162277660168379*bzbz[7]*grad_u_z[7]+0.3535533905932737*bzbz[0]*grad_u_z[7]+0.3162277660168379*bybz[7]*grad_u_y[7]+0.3535533905932737*bybz[0]*grad_u_y[7]+0.3162277660168379*bxbz[7]*grad_u_x[7]+0.3535533905932737*bxbz[0]*grad_u_x[7]+0.3535533905932737*bzbz[1]*grad_u_z[6]+0.3535533905932737*bybz[1]*grad_u_y[6]+0.3535533905932737*bxbz[1]*grad_u_x[6]+0.3535533905932737*grad_u_z[1]*bzbz[6]+0.3535533905932737*grad_u_y[1]*bybz[6]+0.3535533905932737*grad_u_x[1]*bxbz[6]+0.3535533905932737*bzbz[2]*grad_u_z[5]+0.3535533905932737*bybz[2]*grad_u_y[5]+0.3535533905932737*bxbz[2]*grad_u_x[5]+0.3535533905932737*grad_u_z[2]*bzbz[5]+0.3535533905932737*grad_u_y[2]*bybz[5]+0.3535533905932737*grad_u_x[2]*bxbz[5]+0.3535533905932737*bzbz[3]*grad_u_z[4]+0.3535533905932737*bybz[3]*grad_u_y[4]+0.3535533905932737*bxbz[3]*grad_u_x[4]+0.3535533905932737*grad_u_z[3]*bzbz[4]+0.3535533905932737*grad_u_y[3]*bybz[4]+0.3535533905932737*grad_u_x[3]*bxbz[4]; 
  bb_grad_u_comp[11] = 0.3162277660168379*grad_u_z[6]*bzbz[23]+0.3162277660168379*grad_u_y[6]*bybz[23]+0.3162277660168379*grad_u_x[6]*bxbz[23]+0.3162277660168379*grad_u_z[2]*bzbz[20]+0.3162277660168379*grad_u_y[2]*bybz[20]+0.3162277660168379*grad_u_x[2]*bxbz[20]+0.282842712474619*grad_u_z[7]*bzbz[18]+0.282842712474619*grad_u_y[7]*bybz[18]+0.282842712474619*grad_u_x[7]*bxbz[18]+0.3535533905932737*grad_u_z[3]*bzbz[17]+0.3535533905932737*grad_u_y[3]*bybz[17]+0.3535533905932737*grad_u_x[3]*bxbz[17]+0.3535533905932737*grad_u_z[6]*bzbz[13]+0.3535533905932737*grad_u_y[6]*bybz[13]+0.3535533905932737*grad_u_x[6]*bxbz[13]+0.2828427124746191*grad_u_z[4]*bzbz[12]+0.2828427124746191*grad_u_y[4]*bybz[12]+0.2828427124746191*grad_u_x[4]*bxbz[12]+0.3535533905932737*grad_u_z[0]*bzbz[11]+0.3535533905932737*grad_u_y[0]*bybz[11]+0.3535533905932737*grad_u_x[0]*bxbz[11]+0.3162277660168379*grad_u_z[5]*bzbz[10]+0.3162277660168379*grad_u_y[5]*bybz[10]+0.3162277660168379*grad_u_x[5]*bxbz[10]+0.3162277660168379*bzbz[5]*grad_u_z[7]+0.3162277660168379*bybz[5]*grad_u_y[7]+0.3162277660168379*bxbz[5]*grad_u_x[7]+0.3535533905932737*grad_u_z[2]*bzbz[7]+0.3535533905932737*grad_u_y[2]*bybz[7]+0.3535533905932737*grad_u_x[2]*bxbz[7]+0.3162277660168379*bzbz[1]*grad_u_z[4]+0.3162277660168379*bybz[1]*grad_u_y[4]+0.3162277660168379*bxbz[1]*grad_u_x[4]+0.3162277660168379*grad_u_z[1]*bzbz[4]+0.3162277660168379*grad_u_y[1]*bybz[4]+0.3162277660168379*grad_u_x[1]*bxbz[4]; 
  bb_grad_u_comp[12] = 0.3162277660168379*grad_u_z[5]*bzbz[23]+0.3162277660168379*grad_u_y[5]*bybz[23]+0.3162277660168379*grad_u_x[5]*bxbz[23]+0.3162277660168379*grad_u_z[1]*bzbz[20]+0.3162277660168379*grad_u_y[1]*bybz[20]+0.3162277660168379*grad_u_x[1]*bxbz[20]+0.3535533905932737*grad_u_z[3]*bzbz[18]+0.3535533905932737*grad_u_y[3]*bybz[18]+0.3535533905932737*grad_u_x[3]*bxbz[18]+0.282842712474619*grad_u_z[7]*bzbz[17]+0.282842712474619*grad_u_y[7]*bybz[17]+0.282842712474619*grad_u_x[7]*bxbz[17]+0.3535533905932737*grad_u_z[5]*bzbz[14]+0.3535533905932737*grad_u_y[5]*bybz[14]+0.3535533905932737*grad_u_x[5]*bxbz[14]+0.3535533905932737*grad_u_z[0]*bzbz[12]+0.3535533905932737*grad_u_y[0]*bybz[12]+0.3535533905932737*grad_u_x[0]*bxbz[12]+0.2828427124746191*grad_u_z[4]*bzbz[11]+0.2828427124746191*grad_u_y[4]*bybz[11]+0.2828427124746191*grad_u_x[4]*bxbz[11]+0.3162277660168379*grad_u_z[6]*bzbz[10]+0.3162277660168379*grad_u_y[6]*bybz[10]+0.3162277660168379*grad_u_x[6]*bxbz[10]+0.3535533905932737*grad_u_z[1]*bzbz[8]+0.3535533905932737*grad_u_y[1]*bybz[8]+0.3535533905932737*grad_u_x[1]*bxbz[8]+0.3162277660168379*bzbz[6]*grad_u_z[7]+0.3162277660168379*bybz[6]*grad_u_y[7]+0.3162277660168379*bxbz[6]*grad_u_x[7]+0.3162277660168379*bzbz[2]*grad_u_z[4]+0.3162277660168379*bybz[2]*grad_u_y[4]+0.3162277660168379*bxbz[2]*grad_u_x[4]+0.3162277660168379*grad_u_z[2]*bzbz[4]+0.3162277660168379*grad_u_y[2]*bybz[4]+0.3162277660168379*grad_u_x[2]*bxbz[4]; 
  bb_grad_u_comp[13] = 0.3162277660168379*grad_u_z[6]*bzbz[24]+0.3162277660168379*grad_u_y[6]*bybz[24]+0.3162277660168379*grad_u_x[6]*bxbz[24]+0.3162277660168379*grad_u_z[3]*bzbz[21]+0.3162277660168379*grad_u_y[3]*bybz[21]+0.3162277660168379*grad_u_x[3]*bxbz[21]+0.282842712474619*grad_u_z[7]*bzbz[19]+0.282842712474619*grad_u_y[7]*bybz[19]+0.282842712474619*grad_u_x[7]*bxbz[19]+0.3535533905932737*grad_u_z[2]*bzbz[17]+0.3535533905932737*grad_u_y[2]*bybz[17]+0.3535533905932737*grad_u_x[2]*bxbz[17]+0.2828427124746191*grad_u_z[5]*bzbz[15]+0.2828427124746191*grad_u_y[5]*bybz[15]+0.2828427124746191*grad_u_x[5]*bxbz[15]+0.3535533905932737*grad_u_z[0]*bzbz[13]+0.3535533905932737*grad_u_y[0]*bybz[13]+0.3535533905932737*grad_u_x[0]*bxbz[13]+0.3535533905932737*grad_u_z[6]*bzbz[11]+0.3535533905932737*grad_u_y[6]*bybz[11]+0.3535533905932737*grad_u_x[6]*bxbz[11]+0.3162277660168379*grad_u_z[4]*bzbz[10]+0.3162277660168379*grad_u_y[4]*bybz[10]+0.3162277660168379*grad_u_x[4]*bxbz[10]+0.3162277660168379*bzbz[4]*grad_u_z[7]+0.3162277660168379*bybz[4]*grad_u_y[7]+0.3162277660168379*bxbz[4]*grad_u_x[7]+0.3535533905932737*grad_u_z[3]*bzbz[7]+0.3535533905932737*grad_u_y[3]*bybz[7]+0.3535533905932737*grad_u_x[3]*bxbz[7]+0.3162277660168379*bzbz[1]*grad_u_z[5]+0.3162277660168379*bybz[1]*grad_u_y[5]+0.3162277660168379*bxbz[1]*grad_u_x[5]+0.3162277660168379*grad_u_z[1]*bzbz[5]+0.3162277660168379*grad_u_y[1]*bybz[5]+0.3162277660168379*grad_u_x[1]*bxbz[5]; 
  bb_grad_u_comp[14] = 0.3162277660168379*grad_u_z[5]*bzbz[25]+0.3162277660168379*grad_u_y[5]*bybz[25]+0.3162277660168379*grad_u_x[5]*bxbz[25]+0.3162277660168379*grad_u_z[3]*bzbz[22]+0.3162277660168379*grad_u_y[3]*bybz[22]+0.3162277660168379*grad_u_x[3]*bxbz[22]+0.282842712474619*grad_u_z[7]*bzbz[19]+0.282842712474619*grad_u_y[7]*bybz[19]+0.282842712474619*grad_u_x[7]*bxbz[19]+0.3535533905932737*grad_u_z[1]*bzbz[18]+0.3535533905932737*grad_u_y[1]*bybz[18]+0.3535533905932737*grad_u_x[1]*bxbz[18]+0.2828427124746191*grad_u_z[6]*bzbz[16]+0.2828427124746191*grad_u_y[6]*bybz[16]+0.2828427124746191*grad_u_x[6]*bxbz[16]+0.3535533905932737*grad_u_z[0]*bzbz[14]+0.3535533905932737*grad_u_y[0]*bybz[14]+0.3535533905932737*grad_u_x[0]*bxbz[14]+0.3535533905932737*grad_u_z[5]*bzbz[12]+0.3535533905932737*grad_u_y[5]*bybz[12]+0.3535533905932737*grad_u_x[5]*bxbz[12]+0.3162277660168379*grad_u_z[4]*bzbz[10]+0.3162277660168379*grad_u_y[4]*bybz[10]+0.3162277660168379*grad_u_x[4]*bxbz[10]+0.3535533905932737*grad_u_z[3]*bzbz[8]+0.3535533905932737*grad_u_y[3]*bybz[8]+0.3535533905932737*grad_u_x[3]*bxbz[8]+0.3162277660168379*bzbz[4]*grad_u_z[7]+0.3162277660168379*bybz[4]*grad_u_y[7]+0.3162277660168379*bxbz[4]*grad_u_x[7]+0.3162277660168379*bzbz[2]*grad_u_z[6]+0.3162277660168379*bybz[2]*grad_u_y[6]+0.3162277660168379*bxbz[2]*grad_u_x[6]+0.3162277660168379*grad_u_z[2]*bzbz[6]+0.3162277660168379*grad_u_y[2]*bybz[6]+0.3162277660168379*grad_u_x[2]*bxbz[6]; 
  bb_grad_u_comp[15] = 0.3162277660168379*grad_u_z[4]*bzbz[24]+0.3162277660168379*grad_u_y[4]*bybz[24]+0.3162277660168379*grad_u_x[4]*bxbz[24]+0.3162277660168379*grad_u_z[1]*bzbz[21]+0.3162277660168379*grad_u_y[1]*bybz[21]+0.3162277660168379*grad_u_x[1]*bxbz[21]+0.3535533905932737*grad_u_z[2]*bzbz[19]+0.3535533905932737*grad_u_y[2]*bybz[19]+0.3535533905932737*grad_u_x[2]*bxbz[19]+0.282842712474619*grad_u_z[7]*bzbz[17]+0.282842712474619*grad_u_y[7]*bybz[17]+0.282842712474619*grad_u_x[7]*bxbz[17]+0.3535533905932737*grad_u_z[4]*bzbz[16]+0.3535533905932737*grad_u_y[4]*bybz[16]+0.3535533905932737*grad_u_x[4]*bxbz[16]+0.3535533905932737*grad_u_z[0]*bzbz[15]+0.3535533905932737*grad_u_y[0]*bybz[15]+0.3535533905932737*grad_u_x[0]*bxbz[15]+0.2828427124746191*grad_u_z[5]*bzbz[13]+0.2828427124746191*grad_u_y[5]*bybz[13]+0.2828427124746191*grad_u_x[5]*bxbz[13]+0.3162277660168379*grad_u_z[6]*bzbz[10]+0.3162277660168379*grad_u_y[6]*bybz[10]+0.3162277660168379*grad_u_x[6]*bxbz[10]+0.3535533905932737*grad_u_z[1]*bzbz[9]+0.3535533905932737*grad_u_y[1]*bybz[9]+0.3535533905932737*grad_u_x[1]*bxbz[9]+0.3162277660168379*bzbz[6]*grad_u_z[7]+0.3162277660168379*bybz[6]*grad_u_y[7]+0.3162277660168379*bxbz[6]*grad_u_x[7]+0.3162277660168379*bzbz[3]*grad_u_z[5]+0.3162277660168379*bybz[3]*grad_u_y[5]+0.3162277660168379*bxbz[3]*grad_u_x[5]+0.3162277660168379*grad_u_z[3]*bzbz[5]+0.3162277660168379*grad_u_y[3]*bybz[5]+0.3162277660168379*grad_u_x[3]*bxbz[5]; 
  bb_grad_u_comp[16] = 0.3162277660168379*grad_u_z[4]*bzbz[25]+0.3162277660168379*grad_u_y[4]*bybz[25]+0.3162277660168379*grad_u_x[4]*bxbz[25]+0.3162277660168379*grad_u_z[2]*bzbz[22]+0.3162277660168379*grad_u_y[2]*bybz[22]+0.3162277660168379*grad_u_x[2]*bxbz[22]+0.3535533905932737*grad_u_z[1]*bzbz[19]+0.3535533905932737*grad_u_y[1]*bybz[19]+0.3535533905932737*grad_u_x[1]*bxbz[19]+0.282842712474619*grad_u_z[7]*bzbz[18]+0.282842712474619*grad_u_y[7]*bybz[18]+0.282842712474619*grad_u_x[7]*bxbz[18]+0.3535533905932737*grad_u_z[0]*bzbz[16]+0.3535533905932737*grad_u_y[0]*bybz[16]+0.3535533905932737*grad_u_x[0]*bxbz[16]+0.3535533905932737*grad_u_z[4]*bzbz[15]+0.3535533905932737*grad_u_y[4]*bybz[15]+0.3535533905932737*grad_u_x[4]*bxbz[15]+0.2828427124746191*grad_u_z[6]*bzbz[14]+0.2828427124746191*grad_u_y[6]*bybz[14]+0.2828427124746191*grad_u_x[6]*bxbz[14]+0.3162277660168379*grad_u_z[5]*bzbz[10]+0.3162277660168379*grad_u_y[5]*bybz[10]+0.3162277660168379*grad_u_x[5]*bxbz[10]+0.3535533905932737*grad_u_z[2]*bzbz[9]+0.3535533905932737*grad_u_y[2]*bybz[9]+0.3535533905932737*grad_u_x[2]*bxbz[9]+0.3162277660168379*bzbz[5]*grad_u_z[7]+0.3162277660168379*bybz[5]*grad_u_y[7]+0.3162277660168379*bxbz[5]*grad_u_x[7]+0.3162277660168379*bzbz[3]*grad_u_z[6]+0.3162277660168379*bybz[3]*grad_u_y[6]+0.3162277660168379*bxbz[3]*grad_u_x[6]+0.3162277660168379*grad_u_z[3]*bzbz[6]+0.3162277660168379*grad_u_y[3]*bybz[6]+0.3162277660168379*grad_u_x[3]*bxbz[6]; 
  bb_grad_u_comp[17] = 0.2828427124746191*grad_u_z[6]*bzbz[26]+0.2828427124746191*grad_u_y[6]*bybz[26]+0.2828427124746191*grad_u_x[6]*bxbz[26]+0.2529822128134704*grad_u_z[7]*bzbz[25]+0.2529822128134704*grad_u_y[7]*bybz[25]+0.2529822128134704*grad_u_x[7]*bxbz[25]+0.3162277660168379*grad_u_z[3]*bzbz[24]+0.3162277660168379*grad_u_y[3]*bybz[24]+0.3162277660168379*grad_u_x[3]*bxbz[24]+0.3162277660168379*grad_u_z[2]*bzbz[23]+0.3162277660168379*grad_u_y[2]*bybz[23]+0.3162277660168379*grad_u_x[2]*bxbz[23]+0.3162277660168379*grad_u_z[6]*bzbz[21]+0.3162277660168379*grad_u_y[6]*bybz[21]+0.3162277660168379*grad_u_x[6]*bxbz[21]+0.3162277660168379*grad_u_z[6]*bzbz[20]+0.3162277660168379*grad_u_y[6]*bybz[20]+0.3162277660168379*grad_u_x[6]*bxbz[20]+0.2828427124746191*grad_u_z[5]*bzbz[19]+0.2828427124746191*grad_u_y[5]*bybz[19]+0.2828427124746191*grad_u_x[5]*bxbz[19]+0.2828427124746191*grad_u_z[4]*bzbz[18]+0.2828427124746191*grad_u_y[4]*bybz[18]+0.2828427124746191*grad_u_x[4]*bxbz[18]+0.3535533905932737*grad_u_z[0]*bzbz[17]+0.3535533905932737*grad_u_y[0]*bybz[17]+0.3535533905932737*grad_u_x[0]*bxbz[17]+0.282842712474619*grad_u_z[7]*bzbz[15]+0.282842712474619*grad_u_y[7]*bybz[15]+0.282842712474619*grad_u_x[7]*bxbz[15]+0.3535533905932737*grad_u_z[2]*bzbz[13]+0.3535533905932737*grad_u_y[2]*bybz[13]+0.3535533905932737*grad_u_x[2]*bxbz[13]+0.282842712474619*grad_u_z[7]*bzbz[12]+0.282842712474619*grad_u_y[7]*bybz[12]+0.282842712474619*grad_u_x[7]*bxbz[12]+0.3535533905932737*grad_u_z[3]*bzbz[11]+0.3535533905932737*grad_u_y[3]*bybz[11]+0.3535533905932737*grad_u_x[3]*bxbz[11]+0.3162277660168379*grad_u_z[1]*bzbz[10]+0.3162277660168379*grad_u_y[1]*bybz[10]+0.3162277660168379*grad_u_x[1]*bxbz[10]+0.3162277660168379*bzbz[1]*grad_u_z[7]+0.3162277660168379*bybz[1]*grad_u_y[7]+0.3162277660168379*bxbz[1]*grad_u_x[7]+0.3535533905932737*grad_u_z[6]*bzbz[7]+0.3535533905932737*grad_u_y[6]*bybz[7]+0.3535533905932737*grad_u_x[6]*bxbz[7]+0.3162277660168379*bzbz[4]*grad_u_z[5]+0.3162277660168379*bybz[4]*grad_u_y[5]+0.3162277660168379*bxbz[4]*grad_u_x[5]+0.3162277660168379*grad_u_z[4]*bzbz[5]+0.3162277660168379*grad_u_y[4]*bybz[5]+0.3162277660168379*grad_u_x[4]*bxbz[5]; 
  bb_grad_u_comp[18] = 0.2828427124746191*grad_u_z[5]*bzbz[26]+0.2828427124746191*grad_u_y[5]*bybz[26]+0.2828427124746191*grad_u_x[5]*bxbz[26]+0.3162277660168379*grad_u_z[3]*bzbz[25]+0.3162277660168379*grad_u_y[3]*bybz[25]+0.3162277660168379*grad_u_x[3]*bxbz[25]+0.2529822128134704*grad_u_z[7]*bzbz[24]+0.2529822128134704*grad_u_y[7]*bybz[24]+0.2529822128134704*grad_u_x[7]*bxbz[24]+0.3162277660168379*grad_u_z[1]*bzbz[23]+0.3162277660168379*grad_u_y[1]*bybz[23]+0.3162277660168379*grad_u_x[1]*bxbz[23]+0.3162277660168379*grad_u_z[5]*bzbz[22]+0.3162277660168379*grad_u_y[5]*bybz[22]+0.3162277660168379*grad_u_x[5]*bxbz[22]+0.3162277660168379*grad_u_z[5]*bzbz[20]+0.3162277660168379*grad_u_y[5]*bybz[20]+0.3162277660168379*grad_u_x[5]*bxbz[20]+0.2828427124746191*grad_u_z[6]*bzbz[19]+0.2828427124746191*grad_u_y[6]*bybz[19]+0.2828427124746191*grad_u_x[6]*bxbz[19]+0.3535533905932737*grad_u_z[0]*bzbz[18]+0.3535533905932737*grad_u_y[0]*bybz[18]+0.3535533905932737*grad_u_x[0]*bxbz[18]+0.2828427124746191*grad_u_z[4]*bzbz[17]+0.2828427124746191*grad_u_y[4]*bybz[17]+0.2828427124746191*grad_u_x[4]*bxbz[17]+0.282842712474619*grad_u_z[7]*bzbz[16]+0.282842712474619*grad_u_y[7]*bybz[16]+0.282842712474619*grad_u_x[7]*bxbz[16]+0.3535533905932737*grad_u_z[1]*bzbz[14]+0.3535533905932737*grad_u_y[1]*bybz[14]+0.3535533905932737*grad_u_x[1]*bxbz[14]+0.3535533905932737*grad_u_z[3]*bzbz[12]+0.3535533905932737*grad_u_y[3]*bybz[12]+0.3535533905932737*grad_u_x[3]*bxbz[12]+0.282842712474619*grad_u_z[7]*bzbz[11]+0.282842712474619*grad_u_y[7]*bybz[11]+0.282842712474619*grad_u_x[7]*bxbz[11]+0.3162277660168379*grad_u_z[2]*bzbz[10]+0.3162277660168379*grad_u_y[2]*bybz[10]+0.3162277660168379*grad_u_x[2]*bxbz[10]+0.3535533905932737*grad_u_z[5]*bzbz[8]+0.3535533905932737*grad_u_y[5]*bybz[8]+0.3535533905932737*grad_u_x[5]*bxbz[8]+0.3162277660168379*bzbz[2]*grad_u_z[7]+0.3162277660168379*bybz[2]*grad_u_y[7]+0.3162277660168379*bxbz[2]*grad_u_x[7]+0.3162277660168379*bzbz[4]*grad_u_z[6]+0.3162277660168379*bybz[4]*grad_u_y[6]+0.3162277660168379*bxbz[4]*grad_u_x[6]+0.3162277660168379*grad_u_z[4]*bzbz[6]+0.3162277660168379*grad_u_y[4]*bybz[6]+0.3162277660168379*grad_u_x[4]*bxbz[6]; 
  bb_grad_u_comp[19] = 0.2828427124746191*grad_u_z[4]*bzbz[26]+0.2828427124746191*grad_u_y[4]*bybz[26]+0.2828427124746191*grad_u_x[4]*bxbz[26]+0.3162277660168379*grad_u_z[2]*bzbz[25]+0.3162277660168379*grad_u_y[2]*bybz[25]+0.3162277660168379*grad_u_x[2]*bxbz[25]+0.3162277660168379*grad_u_z[1]*bzbz[24]+0.3162277660168379*grad_u_y[1]*bybz[24]+0.3162277660168379*grad_u_x[1]*bxbz[24]+0.2529822128134704*grad_u_z[7]*bzbz[23]+0.2529822128134704*grad_u_y[7]*bybz[23]+0.2529822128134704*grad_u_x[7]*bxbz[23]+0.3162277660168379*grad_u_z[4]*bzbz[22]+0.3162277660168379*grad_u_y[4]*bybz[22]+0.3162277660168379*grad_u_x[4]*bxbz[22]+0.3162277660168379*grad_u_z[4]*bzbz[21]+0.3162277660168379*grad_u_y[4]*bybz[21]+0.3162277660168379*grad_u_x[4]*bxbz[21]+0.3535533905932737*grad_u_z[0]*bzbz[19]+0.3535533905932737*grad_u_y[0]*bybz[19]+0.3535533905932737*grad_u_x[0]*bxbz[19]+0.2828427124746191*grad_u_z[6]*bzbz[18]+0.2828427124746191*grad_u_y[6]*bybz[18]+0.2828427124746191*grad_u_x[6]*bxbz[18]+0.2828427124746191*grad_u_z[5]*bzbz[17]+0.2828427124746191*grad_u_y[5]*bybz[17]+0.2828427124746191*grad_u_x[5]*bxbz[17]+0.3535533905932737*grad_u_z[1]*bzbz[16]+0.3535533905932737*grad_u_y[1]*bybz[16]+0.3535533905932737*grad_u_x[1]*bxbz[16]+0.3535533905932737*grad_u_z[2]*bzbz[15]+0.3535533905932737*grad_u_y[2]*bybz[15]+0.3535533905932737*grad_u_x[2]*bxbz[15]+0.282842712474619*grad_u_z[7]*bzbz[14]+0.282842712474619*grad_u_y[7]*bybz[14]+0.282842712474619*grad_u_x[7]*bxbz[14]+0.282842712474619*grad_u_z[7]*bzbz[13]+0.282842712474619*grad_u_y[7]*bybz[13]+0.282842712474619*grad_u_x[7]*bxbz[13]+0.3162277660168379*grad_u_z[3]*bzbz[10]+0.3162277660168379*grad_u_y[3]*bybz[10]+0.3162277660168379*grad_u_x[3]*bxbz[10]+0.3535533905932737*grad_u_z[4]*bzbz[9]+0.3535533905932737*grad_u_y[4]*bybz[9]+0.3535533905932737*grad_u_x[4]*bxbz[9]+0.3162277660168379*bzbz[3]*grad_u_z[7]+0.3162277660168379*bybz[3]*grad_u_y[7]+0.3162277660168379*bxbz[3]*grad_u_x[7]+0.3162277660168379*bzbz[5]*grad_u_z[6]+0.3162277660168379*bybz[5]*grad_u_y[6]+0.3162277660168379*bxbz[5]*grad_u_x[6]+0.3162277660168379*grad_u_z[5]*bzbz[6]+0.3162277660168379*grad_u_y[5]*bybz[6]+0.3162277660168379*grad_u_x[5]*bxbz[6]; 
  bb_grad_u_comp[20] = 0.3535533905932737*grad_u_z[3]*bzbz[23]+0.3535533905932737*grad_u_y[3]*bybz[23]+0.3535533905932737*grad_u_x[3]*bxbz[23]+0.3535533905932737*grad_u_z[0]*bzbz[20]+0.3535533905932737*grad_u_y[0]*bybz[20]+0.3535533905932737*grad_u_x[0]*bxbz[20]+0.3162277660168379*grad_u_z[5]*bzbz[18]+0.3162277660168379*grad_u_y[5]*bybz[18]+0.3162277660168379*grad_u_x[5]*bxbz[18]+0.3162277660168379*grad_u_z[6]*bzbz[17]+0.3162277660168379*grad_u_y[6]*bybz[17]+0.3162277660168379*grad_u_x[6]*bxbz[17]+0.3162277660168379*grad_u_z[1]*bzbz[12]+0.3162277660168379*grad_u_y[1]*bybz[12]+0.3162277660168379*grad_u_x[1]*bxbz[12]+0.3162277660168379*grad_u_z[2]*bzbz[11]+0.3162277660168379*grad_u_y[2]*bybz[11]+0.3162277660168379*grad_u_x[2]*bxbz[11]+0.2828427124746191*grad_u_z[7]*bzbz[10]+0.2828427124746191*grad_u_y[7]*bybz[10]+0.2828427124746191*grad_u_x[7]*bxbz[10]+0.2828427124746191*bzbz[4]*grad_u_z[4]+0.2828427124746191*bybz[4]*grad_u_y[4]+0.2828427124746191*bxbz[4]*grad_u_x[4]; 
  bb_grad_u_comp[21] = 0.3535533905932737*grad_u_z[2]*bzbz[24]+0.3535533905932737*grad_u_y[2]*bybz[24]+0.3535533905932737*grad_u_x[2]*bxbz[24]+0.3535533905932737*grad_u_z[0]*bzbz[21]+0.3535533905932737*grad_u_y[0]*bybz[21]+0.3535533905932737*grad_u_x[0]*bxbz[21]+0.3162277660168379*grad_u_z[4]*bzbz[19]+0.3162277660168379*grad_u_y[4]*bybz[19]+0.3162277660168379*grad_u_x[4]*bxbz[19]+0.3162277660168379*grad_u_z[6]*bzbz[17]+0.3162277660168379*grad_u_y[6]*bybz[17]+0.3162277660168379*grad_u_x[6]*bxbz[17]+0.3162277660168379*grad_u_z[1]*bzbz[15]+0.3162277660168379*grad_u_y[1]*bybz[15]+0.3162277660168379*grad_u_x[1]*bxbz[15]+0.3162277660168379*grad_u_z[3]*bzbz[13]+0.3162277660168379*grad_u_y[3]*bybz[13]+0.3162277660168379*grad_u_x[3]*bxbz[13]+0.2828427124746191*grad_u_z[7]*bzbz[10]+0.2828427124746191*grad_u_y[7]*bybz[10]+0.2828427124746191*grad_u_x[7]*bxbz[10]+0.2828427124746191*bzbz[5]*grad_u_z[5]+0.2828427124746191*bybz[5]*grad_u_y[5]+0.2828427124746191*bxbz[5]*grad_u_x[5]; 
  bb_grad_u_comp[22] = 0.3535533905932737*grad_u_z[1]*bzbz[25]+0.3535533905932737*grad_u_y[1]*bybz[25]+0.3535533905932737*grad_u_x[1]*bxbz[25]+0.3535533905932737*grad_u_z[0]*bzbz[22]+0.3535533905932737*grad_u_y[0]*bybz[22]+0.3535533905932737*grad_u_x[0]*bxbz[22]+0.3162277660168379*grad_u_z[4]*bzbz[19]+0.3162277660168379*grad_u_y[4]*bybz[19]+0.3162277660168379*grad_u_x[4]*bxbz[19]+0.3162277660168379*grad_u_z[5]*bzbz[18]+0.3162277660168379*grad_u_y[5]*bybz[18]+0.3162277660168379*grad_u_x[5]*bxbz[18]+0.3162277660168379*grad_u_z[2]*bzbz[16]+0.3162277660168379*grad_u_y[2]*bybz[16]+0.3162277660168379*grad_u_x[2]*bxbz[16]+0.3162277660168379*grad_u_z[3]*bzbz[14]+0.3162277660168379*grad_u_y[3]*bybz[14]+0.3162277660168379*grad_u_x[3]*bxbz[14]+0.2828427124746191*grad_u_z[7]*bzbz[10]+0.2828427124746191*grad_u_y[7]*bybz[10]+0.2828427124746191*grad_u_x[7]*bxbz[10]+0.2828427124746191*bzbz[6]*grad_u_z[6]+0.2828427124746191*bybz[6]*grad_u_y[6]+0.2828427124746191*bxbz[6]*grad_u_x[6]; 
  bb_grad_u_comp[23] = 0.3162277660168379*grad_u_z[3]*bzbz[26]+0.3162277660168379*grad_u_y[3]*bybz[26]+0.3162277660168379*grad_u_x[3]*bxbz[26]+0.2828427124746191*grad_u_z[5]*bzbz[25]+0.2828427124746191*grad_u_y[5]*bybz[25]+0.2828427124746191*grad_u_x[5]*bxbz[25]+0.2828427124746191*grad_u_z[6]*bzbz[24]+0.2828427124746191*grad_u_y[6]*bybz[24]+0.2828427124746191*grad_u_x[6]*bxbz[24]+0.3535533905932737*grad_u_z[0]*bzbz[23]+0.3535533905932737*grad_u_y[0]*bybz[23]+0.3535533905932737*grad_u_x[0]*bxbz[23]+0.3535533905932737*grad_u_z[3]*bzbz[20]+0.3535533905932737*grad_u_y[3]*bybz[20]+0.3535533905932737*grad_u_x[3]*bxbz[20]+0.2529822128134704*grad_u_z[7]*bzbz[19]+0.2529822128134704*grad_u_y[7]*bybz[19]+0.2529822128134704*grad_u_x[7]*bxbz[19]+0.3162277660168379*grad_u_z[1]*bzbz[18]+0.3162277660168379*grad_u_y[1]*bybz[18]+0.3162277660168379*grad_u_x[1]*bxbz[18]+0.3162277660168379*grad_u_z[2]*bzbz[17]+0.3162277660168379*grad_u_y[2]*bybz[17]+0.3162277660168379*grad_u_x[2]*bxbz[17]+0.3162277660168379*grad_u_z[5]*bzbz[12]+0.3162277660168379*grad_u_y[5]*bybz[12]+0.3162277660168379*grad_u_x[5]*bxbz[12]+0.3162277660168379*grad_u_z[6]*bzbz[11]+0.3162277660168379*grad_u_y[6]*bybz[11]+0.3162277660168379*grad_u_x[6]*bxbz[11]+0.2828427124746191*grad_u_z[4]*bzbz[10]+0.2828427124746191*grad_u_y[4]*bybz[10]+0.2828427124746191*grad_u_x[4]*bxbz[10]+0.2828427124746191*bzbz[4]*grad_u_z[7]+0.2828427124746191*bybz[4]*grad_u_y[7]+0.2828427124746191*bxbz[4]*grad_u_x[7]; 
  bb_grad_u_comp[24] = 0.3162277660168379*grad_u_z[2]*bzbz[26]+0.3162277660168379*grad_u_y[2]*bybz[26]+0.3162277660168379*grad_u_x[2]*bxbz[26]+0.2828427124746191*grad_u_z[4]*bzbz[25]+0.2828427124746191*grad_u_y[4]*bybz[25]+0.2828427124746191*grad_u_x[4]*bxbz[25]+0.3535533905932737*grad_u_z[0]*bzbz[24]+0.3535533905932737*grad_u_y[0]*bybz[24]+0.3535533905932737*grad_u_x[0]*bxbz[24]+0.2828427124746191*grad_u_z[6]*bzbz[23]+0.2828427124746191*grad_u_y[6]*bybz[23]+0.2828427124746191*grad_u_x[6]*bxbz[23]+0.3535533905932737*grad_u_z[2]*bzbz[21]+0.3535533905932737*grad_u_y[2]*bybz[21]+0.3535533905932737*grad_u_x[2]*bxbz[21]+0.3162277660168379*grad_u_z[1]*bzbz[19]+0.3162277660168379*grad_u_y[1]*bybz[19]+0.3162277660168379*grad_u_x[1]*bxbz[19]+0.2529822128134704*grad_u_z[7]*bzbz[18]+0.2529822128134704*grad_u_y[7]*bybz[18]+0.2529822128134704*grad_u_x[7]*bxbz[18]+0.3162277660168379*grad_u_z[3]*bzbz[17]+0.3162277660168379*grad_u_y[3]*bybz[17]+0.3162277660168379*grad_u_x[3]*bxbz[17]+0.3162277660168379*grad_u_z[4]*bzbz[15]+0.3162277660168379*grad_u_y[4]*bybz[15]+0.3162277660168379*grad_u_x[4]*bxbz[15]+0.3162277660168379*grad_u_z[6]*bzbz[13]+0.3162277660168379*grad_u_y[6]*bybz[13]+0.3162277660168379*grad_u_x[6]*bxbz[13]+0.2828427124746191*grad_u_z[5]*bzbz[10]+0.2828427124746191*grad_u_y[5]*bybz[10]+0.2828427124746191*grad_u_x[5]*bxbz[10]+0.2828427124746191*bzbz[5]*grad_u_z[7]+0.2828427124746191*bybz[5]*grad_u_y[7]+0.2828427124746191*bxbz[5]*grad_u_x[7]; 
  bb_grad_u_comp[25] = 0.3162277660168379*grad_u_z[1]*bzbz[26]+0.3162277660168379*grad_u_y[1]*bybz[26]+0.3162277660168379*grad_u_x[1]*bxbz[26]+0.3535533905932737*grad_u_z[0]*bzbz[25]+0.3535533905932737*grad_u_y[0]*bybz[25]+0.3535533905932737*grad_u_x[0]*bxbz[25]+0.2828427124746191*grad_u_z[4]*bzbz[24]+0.2828427124746191*grad_u_y[4]*bybz[24]+0.2828427124746191*grad_u_x[4]*bxbz[24]+0.2828427124746191*grad_u_z[5]*bzbz[23]+0.2828427124746191*grad_u_y[5]*bybz[23]+0.2828427124746191*grad_u_x[5]*bxbz[23]+0.3535533905932737*grad_u_z[1]*bzbz[22]+0.3535533905932737*grad_u_y[1]*bybz[22]+0.3535533905932737*grad_u_x[1]*bxbz[22]+0.3162277660168379*grad_u_z[2]*bzbz[19]+0.3162277660168379*grad_u_y[2]*bybz[19]+0.3162277660168379*grad_u_x[2]*bxbz[19]+0.3162277660168379*grad_u_z[3]*bzbz[18]+0.3162277660168379*grad_u_y[3]*bybz[18]+0.3162277660168379*grad_u_x[3]*bxbz[18]+0.2529822128134704*grad_u_z[7]*bzbz[17]+0.2529822128134704*grad_u_y[7]*bybz[17]+0.2529822128134704*grad_u_x[7]*bxbz[17]+0.3162277660168379*grad_u_z[4]*bzbz[16]+0.3162277660168379*grad_u_y[4]*bybz[16]+0.3162277660168379*grad_u_x[4]*bxbz[16]+0.3162277660168379*grad_u_z[5]*bzbz[14]+0.3162277660168379*grad_u_y[5]*bybz[14]+0.3162277660168379*grad_u_x[5]*bxbz[14]+0.2828427124746191*grad_u_z[6]*bzbz[10]+0.2828427124746191*grad_u_y[6]*bybz[10]+0.2828427124746191*grad_u_x[6]*bxbz[10]+0.2828427124746191*bzbz[6]*grad_u_z[7]+0.2828427124746191*bybz[6]*grad_u_y[7]+0.2828427124746191*bxbz[6]*grad_u_x[7]; 
  bb_grad_u_comp[26] = 0.3535533905932737*grad_u_z[0]*bzbz[26]+0.3535533905932737*grad_u_y[0]*bybz[26]+0.3535533905932737*grad_u_x[0]*bxbz[26]+0.3162277660168379*grad_u_z[1]*bzbz[25]+0.3162277660168379*grad_u_y[1]*bybz[25]+0.3162277660168379*grad_u_x[1]*bxbz[25]+0.3162277660168379*grad_u_z[2]*bzbz[24]+0.3162277660168379*grad_u_y[2]*bybz[24]+0.3162277660168379*grad_u_x[2]*bxbz[24]+0.3162277660168379*grad_u_z[3]*bzbz[23]+0.3162277660168379*grad_u_y[3]*bybz[23]+0.3162277660168379*grad_u_x[3]*bxbz[23]+0.2828427124746191*grad_u_z[4]*bzbz[19]+0.2828427124746191*grad_u_y[4]*bybz[19]+0.2828427124746191*grad_u_x[4]*bxbz[19]+0.2828427124746191*grad_u_z[5]*bzbz[18]+0.2828427124746191*grad_u_y[5]*bybz[18]+0.2828427124746191*grad_u_x[5]*bxbz[18]+0.2828427124746191*grad_u_z[6]*bzbz[17]+0.2828427124746191*grad_u_y[6]*bybz[17]+0.2828427124746191*grad_u_x[6]*bxbz[17]+0.2529822128134704*grad_u_z[7]*bzbz[10]+0.2529822128134704*grad_u_y[7]*bybz[10]+0.2529822128134704*grad_u_x[7]*bxbz[10]; 

  bb_grad_u[0] += bb_grad_u_comp[0]; 
  bb_grad_u[1] += bb_grad_u_comp[1]; 
  bb_grad_u[2] += bb_grad_u_comp[2]; 
  bb_grad_u[3] += bb_grad_u_comp[3]; 
  bb_grad_u[4] += bb_grad_u_comp[4]; 
  bb_grad_u[5] += bb_grad_u_comp[5]; 
  bb_grad_u[6] += bb_grad_u_comp[6]; 
  bb_grad_u[7] += bb_grad_u_comp[7]; 
  bb_grad_u[8] += bb_grad_u_comp[8]; 
  bb_grad_u[9] += bb_grad_u_comp[9]; 
  bb_grad_u[10] += bb_grad_u_comp[10]; 
  bb_grad_u[11] += bb_grad_u_comp[11]; 
  bb_grad_u[12] += bb_grad_u_comp[12]; 
  bb_grad_u[13] += bb_grad_u_comp[13]; 
  bb_grad_u[14] += bb_grad_u_comp[14]; 
  bb_grad_u[15] += bb_grad_u_comp[15]; 
  bb_grad_u[16] += bb_grad_u_comp[16]; 
  bb_grad_u[17] += bb_grad_u_comp[17]; 
  bb_grad_u[18] += bb_grad_u_comp[18]; 
  bb_grad_u[19] += bb_grad_u_comp[19]; 
  bb_grad_u[20] += bb_grad_u_comp[20]; 
  bb_grad_u[21] += bb_grad_u_comp[21]; 
  bb_grad_u[22] += bb_grad_u_comp[22]; 
  bb_grad_u[23] += bb_grad_u_comp[23]; 
  bb_grad_u[24] += bb_grad_u_comp[24]; 
  bb_grad_u[25] += bb_grad_u_comp[25]; 
  bb_grad_u[26] += bb_grad_u_comp[26]; 

  p_perp_source[0] += (-0.6666666666666666*nu_c[0])-1.0*grad_u_z[0]+bb_grad_u_comp[0]; 
  p_perp_source[1] += (-0.6666666666666666*nu_c[1])-1.0*grad_u_z[1]+bb_grad_u_comp[1]; 
  p_perp_source[2] += (-0.6666666666666666*nu_c[2])-1.0*grad_u_z[2]+bb_grad_u_comp[2]; 
  p_perp_source[3] += (-0.6666666666666666*nu_c[3])-1.0*grad_u_z[3]+bb_grad_u_comp[3]; 
  p_perp_source[4] += (-0.6666666666666666*nu_c[4])-1.0*grad_u_z[4]+bb_grad_u_comp[4]; 
  p_perp_source[5] += (-0.6666666666666666*nu_c[5])-1.0*grad_u_z[5]+bb_grad_u_comp[5]; 
  p_perp_source[6] += (-0.6666666666666666*nu_c[6])-1.0*grad_u_z[6]+bb_grad_u_comp[6]; 
  p_perp_source[7] += bb_grad_u_comp[7]-0.6666666666666666*nu_c[7]; 
  p_perp_source[8] += bb_grad_u_comp[8]-0.6666666666666666*nu_c[8]; 
  p_perp_source[9] += bb_grad_u_comp[9]-0.6666666666666666*nu_c[9]; 
  p_perp_source[10] += (-0.6666666666666666*nu_c[10])+bb_grad_u_comp[10]-1.0*grad_u_z[7]; 
  p_perp_source[11] += bb_grad_u_comp[11]-0.6666666666666666*nu_c[11]; 
  p_perp_source[12] += bb_grad_u_comp[12]-0.6666666666666666*nu_c[12]; 
  p_perp_source[13] += bb_grad_u_comp[13]-0.6666666666666666*nu_c[13]; 
  p_perp_source[14] += bb_grad_u_comp[14]-0.6666666666666666*nu_c[14]; 
  p_perp_source[15] += bb_grad_u_comp[15]-0.6666666666666666*nu_c[15]; 
  p_perp_source[16] += bb_grad_u_comp[16]-0.6666666666666666*nu_c[16]; 
  p_perp_source[17] += bb_grad_u_comp[17]-0.6666666666666666*nu_c[17]; 
  p_perp_source[18] += bb_grad_u_comp[18]-0.6666666666666666*nu_c[18]; 
  p_perp_source[19] += bb_grad_u_comp[19]-0.6666666666666666*nu_c[19]; 
  p_perp_source[20] += bb_grad_u_comp[20]-0.6666666666666666*nu_c[20]; 
  p_perp_source[21] += bb_grad_u_comp[21]-0.6666666666666666*nu_c[21]; 
  p_perp_source[22] += bb_grad_u_comp[22]-0.6666666666666666*nu_c[22]; 
  p_perp_source[23] += bb_grad_u_comp[23]-0.6666666666666666*nu_c[23]; 
  p_perp_source[24] += bb_grad_u_comp[24]-0.6666666666666666*nu_c[24]; 
  p_perp_source[25] += bb_grad_u_comp[25]-0.6666666666666666*nu_c[25]; 
  p_perp_source[26] += bb_grad_u_comp[26]-0.6666666666666666*nu_c[26]; 

} 
