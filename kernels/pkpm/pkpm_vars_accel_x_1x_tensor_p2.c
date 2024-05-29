#include <gkyl_euler_pkpm_kernels.h> 
#include <gkyl_basis_tensor_1x_p2_surfx1_eval_quad.h> 
GKYL_CU_DH void pkpm_vars_accel_x_1x_tensor_p2(const double *dxv, 
  const double *u_surf_l, const double *u_surf_c, const double *u_surf_r, 
  const double *prim_l, const double *prim_c, const double *prim_r, 
  const double *pkpm_u_c, const double *bb_c, const double *nu_c, 
  double* GKYL_RESTRICT pkpm_lax, double* GKYL_RESTRICT pkpm_accel) 
{ 
  // dxv[NDIM]:       Cell spacing.
  // u_surf_l/c/r: Input surface flow velocity expansion in left/center/right cells in each direction.
  //               [ux_xl, ux_xr, uy_xl, uy_xr, uz_xl, uz_xr, 
  //                ux_yl, ux_yr, uy_yl, uy_yr, uz_yl, uz_yr, 
  //                ux_zl, ux_zr, uy_zl, uy_zr, uz_zl, uz_zr]  
  // prim_l/c/r:   Input volume expansion of primitive variables [1/rho div(p_par b), T_perp/m, m/T_perp, 3*Txx/m, 3*Tyy/m, 3*Tzz/m] in left/center/right cells.
  // pkpm_u_c:     Input volume expansion of flow velocity in center cell.
  // bb_c:         Input volume expansion of magnetic field unit tensor in center cell.
  // nu_c:         Input volume expansion of collisionality in center cell.
  // pkpm_lax:     Surface expansion of pkpm Lax penalization: lambda_i = |u_i| + sqrt(3.0*T_ii/m).
  // pkpm_accel:   Volume expansion of pkpm acceleration variables.

  const double dx1 = 2.0/dxv[0]; 
  const double *ux_c = &pkpm_u_c[0]; 
  const double *uy_c = &pkpm_u_c[2]; 
  const double *uz_c = &pkpm_u_c[4]; 

  const double *bxbx = &bb_c[0]; 
  const double *bxby = &bb_c[3]; 
  const double *bxbz = &bb_c[6]; 
  const double *byby = &bb_c[9]; 
  const double *bybz = &bb_c[12]; 
  const double *bzbz = &bb_c[15]; 

  const double *ux_surf_lr = &u_surf_l[1]; 
  const double *uy_surf_lr = &u_surf_l[3]; 
  const double *uz_surf_lr = &u_surf_l[5]; 

  const double *ux_surf_cl = &u_surf_c[0]; 
  const double *uy_surf_cl = &u_surf_c[2]; 
  const double *uz_surf_cl = &u_surf_c[4]; 

  const double *ux_surf_cr = &u_surf_c[1]; 
  const double *uy_surf_cr = &u_surf_c[3]; 
  const double *uz_surf_cr = &u_surf_c[5]; 

  const double *ux_surf_rl = &u_surf_r[0]; 
  const double *uy_surf_rl = &u_surf_r[2]; 
  const double *uz_surf_rl = &u_surf_r[4]; 

  const double *Tii_l = &prim_l[9]; 
  const double *Tii_c = &prim_c[9]; 
  const double *Tii_r = &prim_r[9]; 

  double *pkpm_lax_l = &pkpm_lax[0]; 
  double *pkpm_lax_r = &pkpm_lax[1]; 

  double *bb_grad_u = &pkpm_accel[3]; 
  double *p_perp_source = &pkpm_accel[9]; 

  double max_u_l = fmax(fabs(ux_surf_lr[0]), fabs(ux_surf_cl[0])); 
  double max_u_r = fmax(fabs(ux_surf_cr[0]), fabs(ux_surf_rl[0])); 
  double Tiil_r = tensor_1x_p2_surfx1_eval_quad_node_0_r(Tii_l); 
  double Tiic_l = tensor_1x_p2_surfx1_eval_quad_node_0_l(Tii_c); 
  double Tiic_r = tensor_1x_p2_surfx1_eval_quad_node_0_r(Tii_c); 
  double Tiir_l = tensor_1x_p2_surfx1_eval_quad_node_0_l(Tii_r); 
  double max_vth_l = fmax(sqrt(fabs(Tiil_r)), sqrt(fabs(Tiic_l))); 
  double max_vth_r = fmax(sqrt(fabs(Tiic_r)), sqrt(fabs(Tiir_l))); 
  pkpm_lax_l[0] = max_u_l + max_vth_l; 
  pkpm_lax_r[0] = max_u_r + max_vth_r; 

  double grad_u_x[3] = {0.0}; 
  double grad_u_y[3] = {0.0}; 
  double grad_u_z[3] = {0.0}; 
  grad_u_x[0] = (0.3535533905932737*ux_surf_rl[0]-0.3535533905932737*ux_surf_lr[0]+0.3535533905932737*ux_surf_cr[0]-0.3535533905932737*ux_surf_cl[0])*dx1; 
  grad_u_x[1] = (0.6123724356957944*(ux_surf_rl[0]+ux_surf_lr[0]+ux_surf_cr[0]+ux_surf_cl[0])-1.732050807568877*ux_c[0])*dx1; 
  grad_u_x[2] = ((-3.872983346207417*ux_c[1])+0.7905694150420947*ux_surf_rl[0]-0.7905694150420947*ux_surf_lr[0]+0.7905694150420947*ux_surf_cr[0]-0.7905694150420947*ux_surf_cl[0])*dx1; 

  grad_u_y[0] = (0.3535533905932737*uy_surf_rl[0]-0.3535533905932737*uy_surf_lr[0]+0.3535533905932737*uy_surf_cr[0]-0.3535533905932737*uy_surf_cl[0])*dx1; 
  grad_u_y[1] = (0.6123724356957944*(uy_surf_rl[0]+uy_surf_lr[0]+uy_surf_cr[0]+uy_surf_cl[0])-1.732050807568877*uy_c[0])*dx1; 
  grad_u_y[2] = ((-3.872983346207417*uy_c[1])+0.7905694150420947*uy_surf_rl[0]-0.7905694150420947*uy_surf_lr[0]+0.7905694150420947*uy_surf_cr[0]-0.7905694150420947*uy_surf_cl[0])*dx1; 

  grad_u_z[0] = (0.3535533905932737*uz_surf_rl[0]-0.3535533905932737*uz_surf_lr[0]+0.3535533905932737*uz_surf_cr[0]-0.3535533905932737*uz_surf_cl[0])*dx1; 
  grad_u_z[1] = (0.6123724356957944*(uz_surf_rl[0]+uz_surf_lr[0]+uz_surf_cr[0]+uz_surf_cl[0])-1.732050807568877*uz_c[0])*dx1; 
  grad_u_z[2] = ((-3.872983346207417*uz_c[1])+0.7905694150420947*uz_surf_rl[0]-0.7905694150420947*uz_surf_lr[0]+0.7905694150420947*uz_surf_cr[0]-0.7905694150420947*uz_surf_cl[0])*dx1; 

  double bb_grad_u_comp[3] = {0.0}; 
  bb_grad_u_comp[0] = 0.7071067811865475*bxbz[2]*grad_u_z[2]+0.7071067811865475*bxby[2]*grad_u_y[2]+0.7071067811865475*bxbx[2]*grad_u_x[2]+0.7071067811865475*bxbz[1]*grad_u_z[1]+0.7071067811865475*bxby[1]*grad_u_y[1]+0.7071067811865475*bxbx[1]*grad_u_x[1]+0.7071067811865475*bxbz[0]*grad_u_z[0]+0.7071067811865475*bxby[0]*grad_u_y[0]+0.7071067811865475*bxbx[0]*grad_u_x[0]; 
  bb_grad_u_comp[1] = 0.6324555320336759*bxbz[1]*grad_u_z[2]+0.6324555320336759*bxby[1]*grad_u_y[2]+0.6324555320336759*bxbx[1]*grad_u_x[2]+0.6324555320336759*grad_u_z[1]*bxbz[2]+0.6324555320336759*grad_u_y[1]*bxby[2]+0.6324555320336759*grad_u_x[1]*bxbx[2]+0.7071067811865475*bxbz[0]*grad_u_z[1]+0.7071067811865475*bxby[0]*grad_u_y[1]+0.7071067811865475*bxbx[0]*grad_u_x[1]+0.7071067811865475*grad_u_z[0]*bxbz[1]+0.7071067811865475*grad_u_y[0]*bxby[1]+0.7071067811865475*grad_u_x[0]*bxbx[1]; 
  bb_grad_u_comp[2] = 0.4517539514526256*bxbz[2]*grad_u_z[2]+0.7071067811865475*bxbz[0]*grad_u_z[2]+0.4517539514526256*bxby[2]*grad_u_y[2]+0.7071067811865475*bxby[0]*grad_u_y[2]+0.4517539514526256*bxbx[2]*grad_u_x[2]+0.7071067811865475*bxbx[0]*grad_u_x[2]+0.7071067811865475*grad_u_z[0]*bxbz[2]+0.7071067811865475*grad_u_y[0]*bxby[2]+0.7071067811865475*grad_u_x[0]*bxbx[2]+0.6324555320336759*bxbz[1]*grad_u_z[1]+0.6324555320336759*bxby[1]*grad_u_y[1]+0.6324555320336759*bxbx[1]*grad_u_x[1]; 

  bb_grad_u[0] += bb_grad_u_comp[0]; 
  bb_grad_u[1] += bb_grad_u_comp[1]; 
  bb_grad_u[2] += bb_grad_u_comp[2]; 

  p_perp_source[0] += (-2.0*nu_c[0])-1.0*grad_u_x[0]+bb_grad_u_comp[0]; 
  p_perp_source[1] += (-2.0*nu_c[1])-1.0*grad_u_x[1]+bb_grad_u_comp[1]; 
  p_perp_source[2] += (-2.0*nu_c[2])-1.0*grad_u_x[2]+bb_grad_u_comp[2]; 

} 
