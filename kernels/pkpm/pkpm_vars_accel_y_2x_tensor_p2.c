#include <gkyl_euler_pkpm_kernels.h> 
#include <gkyl_basis_tensor_2x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH void pkpm_vars_accel_y_2x_tensor_p2(const double *dxv, 
  const double *u_surf_l, const double *u_surf_c, const double *u_surf_r, 
  const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
  const double *pkpm_u_c, const double *bb_c, const double *nu_c, 
  double* GKYL_RESTRICT pkpm_lax, double* GKYL_RESTRICT pkpm_accel) 
{ 
  // dxv[NDIM]:       Cell spacing.
  // u_surf_l/c/r: Input surface flow velocity expansion in left/center/right cells in each direction.
  //               [ux_xl, ux_xr, uy_xl, uy_xr, uz_xl, uz_xr, 
  //                ux_yl, ux_yr, uy_yl, uy_yr, uz_yl, uz_yr, 
  //                ux_zl, ux_zr, uy_zl, uy_zr, uz_zl, uz_zr]  
  // prim_surf_l/c/r: Input surface primitive variables [3*T_ii/m] in left/center/right cells in each direction.
  // pkpm_u_c:     Input volume expansion of flow velocity in center cell.
  // bb_c:         Input volume expansion of magnetic field unit tensor in center cell.
  // nu_c:         Input volume expansion of collisionality in center cell.
  // pkpm_lax:     Surface expansion of pkpm Lax penalization: lambda_i = |u_i| + sqrt(3.0*T_ii/m).
  // pkpm_accel:   Volume expansion of pkpm acceleration variables.

  const double dx1 = 2.0/dxv[1]; 
  const double *ux_c = &pkpm_u_c[0]; 
  const double *uy_c = &pkpm_u_c[4]; 
  const double *uz_c = &pkpm_u_c[8]; 

  const double *bxbx = &bb_c[0]; 
  const double *bxby = &bb_c[9]; 
  const double *bxbz = &bb_c[18]; 
  const double *byby = &bb_c[27]; 
  const double *bybz = &bb_c[36]; 
  const double *bzbz = &bb_c[45]; 

  const double *ux_surf_lr = &u_surf_l[14]; 
  const double *uy_surf_lr = &u_surf_l[18]; 
  const double *uz_surf_lr = &u_surf_l[22]; 

  const double *ux_surf_cl = &u_surf_c[12]; 
  const double *uy_surf_cl = &u_surf_c[16]; 
  const double *uz_surf_cl = &u_surf_c[20]; 

  const double *ux_surf_cr = &u_surf_c[14]; 
  const double *uy_surf_cr = &u_surf_c[18]; 
  const double *uz_surf_cr = &u_surf_c[22]; 

  const double *ux_surf_rl = &u_surf_r[12]; 
  const double *uy_surf_rl = &u_surf_r[16]; 
  const double *uz_surf_rl = &u_surf_r[20]; 

  const double *Tii_surf_lr = &prim_surf_l[9]; 
  const double *Tii_surf_cl = &prim_surf_c[6]; 
  const double *Tii_surf_cr = &prim_surf_c[9]; 
  const double *Tii_surf_rl = &prim_surf_r[6]; 

  double *pkpm_lax_l = &pkpm_lax[6]; 
  double *pkpm_lax_r = &pkpm_lax[9]; 

  double *bb_grad_u = &pkpm_accel[9]; 
  double *p_perp_source = &pkpm_accel[27]; 

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
  double pkpm_lax_quad_l[3] = {0.0}; 
  double pkpm_lax_quad_r[3] = {0.0}; 

  ul_r = 0.7071067811865475*uy_surf_lr[0]-0.9486832980505137*uy_surf_lr[1]; 
  uc_l = 0.7071067811865475*uy_surf_cl[0]-0.9486832980505137*uy_surf_cl[1]; 
  uc_r = 0.7071067811865475*uy_surf_cr[0]-0.9486832980505137*uy_surf_cr[1]; 
  ur_l = 0.7071067811865475*uy_surf_rl[0]-0.9486832980505137*uy_surf_rl[1]; 
  uQuad_l = fmax(fabs(ul_r), fabs(uc_l)); 
  uQuad_r = fmax(fabs(uc_r), fabs(ur_l)); 
  Tiil_r = 0.6324555320336759*Tii_surf_lr[2]-0.9486832980505137*Tii_surf_lr[1]+0.7071067811865475*Tii_surf_lr[0]; 
  Tiic_l = 0.6324555320336759*Tii_surf_cl[2]-0.9486832980505137*Tii_surf_cl[1]+0.7071067811865475*Tii_surf_cl[0]; 
  Tiic_r = 0.6324555320336759*Tii_surf_cr[2]-0.9486832980505137*Tii_surf_cr[1]+0.7071067811865475*Tii_surf_cr[0]; 
  Tiir_l = 0.6324555320336759*Tii_surf_rl[2]-0.9486832980505137*Tii_surf_rl[1]+0.7071067811865475*Tii_surf_rl[0]; 
  TiiQuad_l = fmax(sqrt(fabs(Tiil_r)), sqrt(fabs(Tiic_l))); 
  TiiQuad_r = fmax(sqrt(fabs(Tiic_r)), sqrt(fabs(Tiir_l))); 
  pkpm_lax_quad_l[0] = uQuad_l + TiiQuad_l; 
  pkpm_lax_quad_r[0] = uQuad_r + TiiQuad_r; 

  ul_r = 0.7071067811865475*uy_surf_lr[0]; 
  uc_l = 0.7071067811865475*uy_surf_cl[0]; 
  uc_r = 0.7071067811865475*uy_surf_cr[0]; 
  ur_l = 0.7071067811865475*uy_surf_rl[0]; 
  uQuad_l = fmax(fabs(ul_r), fabs(uc_l)); 
  uQuad_r = fmax(fabs(uc_r), fabs(ur_l)); 
  Tiil_r = 0.7071067811865475*Tii_surf_lr[0]-0.7905694150420947*Tii_surf_lr[2]; 
  Tiic_l = 0.7071067811865475*Tii_surf_cl[0]-0.7905694150420947*Tii_surf_cl[2]; 
  Tiic_r = 0.7071067811865475*Tii_surf_cr[0]-0.7905694150420947*Tii_surf_cr[2]; 
  Tiir_l = 0.7071067811865475*Tii_surf_rl[0]-0.7905694150420947*Tii_surf_rl[2]; 
  TiiQuad_l = fmax(sqrt(fabs(Tiil_r)), sqrt(fabs(Tiic_l))); 
  TiiQuad_r = fmax(sqrt(fabs(Tiic_r)), sqrt(fabs(Tiir_l))); 
  pkpm_lax_quad_l[1] = uQuad_l + TiiQuad_l; 
  pkpm_lax_quad_r[1] = uQuad_r + TiiQuad_r; 

  ul_r = 0.9486832980505137*uy_surf_lr[1]+0.7071067811865475*uy_surf_lr[0]; 
  uc_l = 0.9486832980505137*uy_surf_cl[1]+0.7071067811865475*uy_surf_cl[0]; 
  uc_r = 0.9486832980505137*uy_surf_cr[1]+0.7071067811865475*uy_surf_cr[0]; 
  ur_l = 0.9486832980505137*uy_surf_rl[1]+0.7071067811865475*uy_surf_rl[0]; 
  uQuad_l = fmax(fabs(ul_r), fabs(uc_l)); 
  uQuad_r = fmax(fabs(uc_r), fabs(ur_l)); 
  Tiil_r = 0.6324555320336759*Tii_surf_lr[2]+0.9486832980505137*Tii_surf_lr[1]+0.7071067811865475*Tii_surf_lr[0]; 
  Tiic_l = 0.6324555320336759*Tii_surf_cl[2]+0.9486832980505137*Tii_surf_cl[1]+0.7071067811865475*Tii_surf_cl[0]; 
  Tiic_r = 0.6324555320336759*Tii_surf_cr[2]+0.9486832980505137*Tii_surf_cr[1]+0.7071067811865475*Tii_surf_cr[0]; 
  Tiir_l = 0.6324555320336759*Tii_surf_rl[2]+0.9486832980505137*Tii_surf_rl[1]+0.7071067811865475*Tii_surf_rl[0]; 
  TiiQuad_l = fmax(sqrt(fabs(Tiil_r)), sqrt(fabs(Tiic_l))); 
  TiiQuad_r = fmax(sqrt(fabs(Tiic_r)), sqrt(fabs(Tiir_l))); 
  pkpm_lax_quad_l[2] = uQuad_l + TiiQuad_l; 
  pkpm_lax_quad_r[2] = uQuad_r + TiiQuad_r; 

  tensor_2x_p2_upwind_quad_to_modal(pkpm_lax_quad_l, pkpm_lax_l); 
  tensor_2x_p2_upwind_quad_to_modal(pkpm_lax_quad_r, pkpm_lax_r); 

  double grad_u_x[9] = {0.0}; 
  double grad_u_y[9] = {0.0}; 
  double grad_u_z[9] = {0.0}; 
  grad_u_x[0] = (0.3535533905932737*ux_surf_rl[0]-0.3535533905932737*ux_surf_lr[0]+0.3535533905932737*ux_surf_cr[0]-0.3535533905932737*ux_surf_cl[0])*dx1; 
  grad_u_x[1] = (0.3535533905932737*ux_surf_rl[1]-0.3535533905932737*ux_surf_lr[1]+0.3535533905932737*ux_surf_cr[1]-0.3535533905932737*ux_surf_cl[1])*dx1; 
  grad_u_x[2] = (0.6123724356957944*(ux_surf_rl[0]+ux_surf_lr[0]+ux_surf_cr[0]+ux_surf_cl[0])-1.732050807568877*ux_c[0])*dx1; 
  grad_u_x[3] = (0.6123724356957944*(ux_surf_rl[1]+ux_surf_lr[1]+ux_surf_cr[1]+ux_surf_cl[1])-1.732050807568877*ux_c[1])*dx1; 
  grad_u_x[5] = ((-3.872983346207417*ux_c[2])+0.7905694150420947*ux_surf_rl[0]-0.7905694150420947*ux_surf_lr[0]+0.7905694150420947*ux_surf_cr[0]-0.7905694150420947*ux_surf_cl[0])*dx1; 
  grad_u_x[7] = ((-3.872983346207417*ux_c[3])+0.7905694150420948*ux_surf_rl[1]-0.7905694150420948*ux_surf_lr[1]+0.7905694150420948*ux_surf_cr[1]-0.7905694150420948*ux_surf_cl[1])*dx1; 

  grad_u_y[0] = (0.3535533905932737*uy_surf_rl[0]-0.3535533905932737*uy_surf_lr[0]+0.3535533905932737*uy_surf_cr[0]-0.3535533905932737*uy_surf_cl[0])*dx1; 
  grad_u_y[1] = (0.3535533905932737*uy_surf_rl[1]-0.3535533905932737*uy_surf_lr[1]+0.3535533905932737*uy_surf_cr[1]-0.3535533905932737*uy_surf_cl[1])*dx1; 
  grad_u_y[2] = (0.6123724356957944*(uy_surf_rl[0]+uy_surf_lr[0]+uy_surf_cr[0]+uy_surf_cl[0])-1.732050807568877*uy_c[0])*dx1; 
  grad_u_y[3] = (0.6123724356957944*(uy_surf_rl[1]+uy_surf_lr[1]+uy_surf_cr[1]+uy_surf_cl[1])-1.732050807568877*uy_c[1])*dx1; 
  grad_u_y[5] = ((-3.872983346207417*uy_c[2])+0.7905694150420947*uy_surf_rl[0]-0.7905694150420947*uy_surf_lr[0]+0.7905694150420947*uy_surf_cr[0]-0.7905694150420947*uy_surf_cl[0])*dx1; 
  grad_u_y[7] = ((-3.872983346207417*uy_c[3])+0.7905694150420948*uy_surf_rl[1]-0.7905694150420948*uy_surf_lr[1]+0.7905694150420948*uy_surf_cr[1]-0.7905694150420948*uy_surf_cl[1])*dx1; 

  grad_u_z[0] = (0.3535533905932737*uz_surf_rl[0]-0.3535533905932737*uz_surf_lr[0]+0.3535533905932737*uz_surf_cr[0]-0.3535533905932737*uz_surf_cl[0])*dx1; 
  grad_u_z[1] = (0.3535533905932737*uz_surf_rl[1]-0.3535533905932737*uz_surf_lr[1]+0.3535533905932737*uz_surf_cr[1]-0.3535533905932737*uz_surf_cl[1])*dx1; 
  grad_u_z[2] = (0.6123724356957944*(uz_surf_rl[0]+uz_surf_lr[0]+uz_surf_cr[0]+uz_surf_cl[0])-1.732050807568877*uz_c[0])*dx1; 
  grad_u_z[3] = (0.6123724356957944*(uz_surf_rl[1]+uz_surf_lr[1]+uz_surf_cr[1]+uz_surf_cl[1])-1.732050807568877*uz_c[1])*dx1; 
  grad_u_z[5] = ((-3.872983346207417*uz_c[2])+0.7905694150420947*uz_surf_rl[0]-0.7905694150420947*uz_surf_lr[0]+0.7905694150420947*uz_surf_cr[0]-0.7905694150420947*uz_surf_cl[0])*dx1; 
  grad_u_z[7] = ((-3.872983346207417*uz_c[3])+0.7905694150420948*uz_surf_rl[1]-0.7905694150420948*uz_surf_lr[1]+0.7905694150420948*uz_surf_cr[1]-0.7905694150420948*uz_surf_cl[1])*dx1; 

  double bb_grad_u_comp[9] = {0.0}; 
  bb_grad_u_comp[0] = 0.5*bybz[7]*grad_u_z[7]+0.5*byby[7]*grad_u_y[7]+0.5*bxby[7]*grad_u_x[7]+0.5*bybz[5]*grad_u_z[5]+0.5*byby[5]*grad_u_y[5]+0.5*bxby[5]*grad_u_x[5]+0.5*bybz[3]*grad_u_z[3]+0.5*byby[3]*grad_u_y[3]+0.5*bxby[3]*grad_u_x[3]+0.5*bybz[2]*grad_u_z[2]+0.5*byby[2]*grad_u_y[2]+0.5*bxby[2]*grad_u_x[2]+0.5*bybz[1]*grad_u_z[1]+0.5*byby[1]*grad_u_y[1]+0.5*bxby[1]*grad_u_x[1]+0.5*bybz[0]*grad_u_z[0]+0.5*byby[0]*grad_u_y[0]+0.5*bxby[0]*grad_u_x[0]; 
  bb_grad_u_comp[1] = 0.447213595499958*grad_u_z[7]*bybz[8]+0.447213595499958*grad_u_y[7]*byby[8]+0.447213595499958*grad_u_x[7]*bxby[8]+0.5000000000000001*bybz[5]*grad_u_z[7]+0.5000000000000001*byby[5]*grad_u_y[7]+0.5000000000000001*bxby[5]*grad_u_x[7]+0.5000000000000001*grad_u_z[5]*bybz[7]+0.5000000000000001*grad_u_y[5]*byby[7]+0.5000000000000001*grad_u_x[5]*bxby[7]+0.447213595499958*grad_u_z[3]*bybz[6]+0.447213595499958*grad_u_y[3]*byby[6]+0.447213595499958*grad_u_x[3]*bxby[6]+0.4472135954999579*grad_u_z[1]*bybz[4]+0.4472135954999579*grad_u_y[1]*byby[4]+0.4472135954999579*grad_u_x[1]*bxby[4]+0.5*bybz[2]*grad_u_z[3]+0.5*byby[2]*grad_u_y[3]+0.5*bxby[2]*grad_u_x[3]+0.5*grad_u_z[2]*bybz[3]+0.5*grad_u_y[2]*byby[3]+0.5*grad_u_x[2]*bxby[3]+0.5*bybz[0]*grad_u_z[1]+0.5*byby[0]*grad_u_y[1]+0.5*bxby[0]*grad_u_x[1]+0.5*grad_u_z[0]*bybz[1]+0.5*grad_u_y[0]*byby[1]+0.5*grad_u_x[0]*bxby[1]; 
  bb_grad_u_comp[2] = 0.447213595499958*bybz[3]*grad_u_z[7]+0.447213595499958*byby[3]*grad_u_y[7]+0.447213595499958*bxby[3]*grad_u_x[7]+0.447213595499958*grad_u_z[3]*bybz[7]+0.447213595499958*grad_u_y[3]*byby[7]+0.447213595499958*grad_u_x[3]*bxby[7]+0.4472135954999579*bybz[2]*grad_u_z[5]+0.4472135954999579*byby[2]*grad_u_y[5]+0.4472135954999579*bxby[2]*grad_u_x[5]+0.4472135954999579*grad_u_z[2]*bybz[5]+0.4472135954999579*grad_u_y[2]*byby[5]+0.4472135954999579*grad_u_x[2]*bxby[5]+0.5*bybz[1]*grad_u_z[3]+0.5*byby[1]*grad_u_y[3]+0.5*bxby[1]*grad_u_x[3]+0.5*grad_u_z[1]*bybz[3]+0.5*grad_u_y[1]*byby[3]+0.5*grad_u_x[1]*bxby[3]+0.5*bybz[0]*grad_u_z[2]+0.5*byby[0]*grad_u_y[2]+0.5*bxby[0]*grad_u_x[2]+0.5*grad_u_z[0]*bybz[2]+0.5*grad_u_y[0]*byby[2]+0.5*grad_u_x[0]*bxby[2]; 
  bb_grad_u_comp[3] = 0.4*grad_u_z[3]*bybz[8]+0.4*grad_u_y[3]*byby[8]+0.4*grad_u_x[3]*bxby[8]+0.4*bybz[6]*grad_u_z[7]+0.447213595499958*bybz[2]*grad_u_z[7]+0.4*byby[6]*grad_u_y[7]+0.447213595499958*byby[2]*grad_u_y[7]+0.4*bxby[6]*grad_u_x[7]+0.447213595499958*bxby[2]*grad_u_x[7]+0.447213595499958*grad_u_z[2]*bybz[7]+0.447213595499958*grad_u_y[2]*byby[7]+0.447213595499958*grad_u_x[2]*bxby[7]+0.447213595499958*grad_u_z[1]*bybz[6]+0.447213595499958*grad_u_y[1]*byby[6]+0.447213595499958*grad_u_x[1]*bxby[6]+0.4472135954999579*bybz[3]*grad_u_z[5]+0.4472135954999579*byby[3]*grad_u_y[5]+0.4472135954999579*bxby[3]*grad_u_x[5]+0.4472135954999579*grad_u_z[3]*bybz[5]+0.4472135954999579*grad_u_y[3]*byby[5]+0.4472135954999579*grad_u_x[3]*bxby[5]+0.4472135954999579*grad_u_z[3]*bybz[4]+0.4472135954999579*grad_u_y[3]*byby[4]+0.4472135954999579*grad_u_x[3]*bxby[4]+0.5*bybz[0]*grad_u_z[3]+0.5*byby[0]*grad_u_y[3]+0.5*bxby[0]*grad_u_x[3]+0.5*grad_u_z[0]*bybz[3]+0.5*grad_u_y[0]*byby[3]+0.5*grad_u_x[0]*bxby[3]+0.5*bybz[1]*grad_u_z[2]+0.5*byby[1]*grad_u_y[2]+0.5*bxby[1]*grad_u_x[2]+0.5*grad_u_z[1]*bybz[2]+0.5*grad_u_y[1]*byby[2]+0.5*grad_u_x[1]*bxby[2]; 
  bb_grad_u_comp[4] = 0.5*grad_u_z[5]*bybz[8]+0.5*grad_u_y[5]*byby[8]+0.5*grad_u_x[5]*bxby[8]+0.4472135954999579*bybz[7]*grad_u_z[7]+0.4472135954999579*byby[7]*grad_u_y[7]+0.4472135954999579*bxby[7]*grad_u_x[7]+0.5000000000000001*grad_u_z[2]*bybz[6]+0.5000000000000001*grad_u_y[2]*byby[6]+0.5000000000000001*grad_u_x[2]*bxby[6]+0.5*grad_u_z[0]*bybz[4]+0.5*grad_u_y[0]*byby[4]+0.5*grad_u_x[0]*bxby[4]+0.4472135954999579*bybz[3]*grad_u_z[3]+0.4472135954999579*byby[3]*grad_u_y[3]+0.4472135954999579*bxby[3]*grad_u_x[3]+0.4472135954999579*bybz[1]*grad_u_z[1]+0.4472135954999579*byby[1]*grad_u_y[1]+0.4472135954999579*bxby[1]*grad_u_x[1]; 
  bb_grad_u_comp[5] = 0.31943828249997*bybz[7]*grad_u_z[7]+0.5000000000000001*bybz[1]*grad_u_z[7]+0.31943828249997*byby[7]*grad_u_y[7]+0.5000000000000001*byby[1]*grad_u_y[7]+0.31943828249997*bxby[7]*grad_u_x[7]+0.5000000000000001*bxby[1]*grad_u_x[7]+0.5000000000000001*grad_u_z[1]*bybz[7]+0.5000000000000001*grad_u_y[1]*byby[7]+0.5000000000000001*grad_u_x[1]*bxby[7]+0.31943828249997*bybz[5]*grad_u_z[5]+0.5*bybz[0]*grad_u_z[5]+0.31943828249997*byby[5]*grad_u_y[5]+0.5*byby[0]*grad_u_y[5]+0.31943828249997*bxby[5]*grad_u_x[5]+0.5*bxby[0]*grad_u_x[5]+0.5*grad_u_z[0]*bybz[5]+0.5*grad_u_y[0]*byby[5]+0.5*grad_u_x[0]*bxby[5]+0.4472135954999579*bybz[3]*grad_u_z[3]+0.4472135954999579*byby[3]*grad_u_y[3]+0.4472135954999579*bxby[3]*grad_u_x[3]+0.4472135954999579*bybz[2]*grad_u_z[2]+0.4472135954999579*byby[2]*grad_u_y[2]+0.4472135954999579*bxby[2]*grad_u_x[2]; 
  bb_grad_u_comp[6] = 0.447213595499958*grad_u_z[2]*bybz[8]+0.447213595499958*grad_u_y[2]*byby[8]+0.447213595499958*grad_u_x[2]*bxby[8]+0.4*bybz[3]*grad_u_z[7]+0.4*byby[3]*grad_u_y[7]+0.4*bxby[3]*grad_u_x[7]+0.4*grad_u_z[3]*bybz[7]+0.4*grad_u_y[3]*byby[7]+0.4*grad_u_x[3]*bxby[7]+0.4472135954999579*grad_u_z[5]*bybz[6]+0.5*grad_u_z[0]*bybz[6]+0.4472135954999579*grad_u_y[5]*byby[6]+0.5*grad_u_y[0]*byby[6]+0.4472135954999579*grad_u_x[5]*bxby[6]+0.5*grad_u_x[0]*bxby[6]+0.5000000000000001*grad_u_z[2]*bybz[4]+0.5000000000000001*grad_u_y[2]*byby[4]+0.5000000000000001*grad_u_x[2]*bxby[4]+0.447213595499958*bybz[1]*grad_u_z[3]+0.447213595499958*byby[1]*grad_u_y[3]+0.447213595499958*bxby[1]*grad_u_x[3]+0.447213595499958*grad_u_z[1]*bybz[3]+0.447213595499958*grad_u_y[1]*byby[3]+0.447213595499958*grad_u_x[1]*bxby[3]; 
  bb_grad_u_comp[7] = 0.2857142857142857*grad_u_z[7]*bybz[8]+0.447213595499958*grad_u_z[1]*bybz[8]+0.2857142857142857*grad_u_y[7]*byby[8]+0.447213595499958*grad_u_y[1]*byby[8]+0.2857142857142857*grad_u_x[7]*bxby[8]+0.447213595499958*grad_u_x[1]*bxby[8]+0.31943828249997*bybz[5]*grad_u_z[7]+0.4472135954999579*bybz[4]*grad_u_z[7]+0.5*bybz[0]*grad_u_z[7]+0.31943828249997*byby[5]*grad_u_y[7]+0.4472135954999579*byby[4]*grad_u_y[7]+0.5*byby[0]*grad_u_y[7]+0.31943828249997*bxby[5]*grad_u_x[7]+0.4472135954999579*bxby[4]*grad_u_x[7]+0.5*bxby[0]*grad_u_x[7]+0.31943828249997*grad_u_z[5]*bybz[7]+0.5*grad_u_z[0]*bybz[7]+0.31943828249997*grad_u_y[5]*byby[7]+0.5*grad_u_y[0]*byby[7]+0.31943828249997*grad_u_x[5]*bxby[7]+0.5*grad_u_x[0]*bxby[7]+0.4*grad_u_z[3]*bybz[6]+0.4*grad_u_y[3]*byby[6]+0.4*grad_u_x[3]*bxby[6]+0.5000000000000001*bybz[1]*grad_u_z[5]+0.5000000000000001*byby[1]*grad_u_y[5]+0.5000000000000001*bxby[1]*grad_u_x[5]+0.5000000000000001*grad_u_z[1]*bybz[5]+0.5000000000000001*grad_u_y[1]*byby[5]+0.5000000000000001*grad_u_x[1]*bxby[5]+0.447213595499958*bybz[2]*grad_u_z[3]+0.447213595499958*byby[2]*grad_u_y[3]+0.447213595499958*bxby[2]*grad_u_x[3]+0.447213595499958*grad_u_z[2]*bybz[3]+0.447213595499958*grad_u_y[2]*byby[3]+0.447213595499958*grad_u_x[2]*bxby[3]; 
  bb_grad_u_comp[8] = 0.31943828249997*grad_u_z[5]*bybz[8]+0.5*grad_u_z[0]*bybz[8]+0.31943828249997*grad_u_y[5]*byby[8]+0.5*grad_u_y[0]*byby[8]+0.31943828249997*grad_u_x[5]*bxby[8]+0.5*grad_u_x[0]*bxby[8]+0.2857142857142857*bybz[7]*grad_u_z[7]+0.447213595499958*bybz[1]*grad_u_z[7]+0.2857142857142857*byby[7]*grad_u_y[7]+0.447213595499958*byby[1]*grad_u_y[7]+0.2857142857142857*bxby[7]*grad_u_x[7]+0.447213595499958*bxby[1]*grad_u_x[7]+0.447213595499958*grad_u_z[1]*bybz[7]+0.447213595499958*grad_u_y[1]*byby[7]+0.447213595499958*grad_u_x[1]*bxby[7]+0.447213595499958*grad_u_z[2]*bybz[6]+0.447213595499958*grad_u_y[2]*byby[6]+0.447213595499958*grad_u_x[2]*bxby[6]+0.5*bybz[4]*grad_u_z[5]+0.5*byby[4]*grad_u_y[5]+0.5*bxby[4]*grad_u_x[5]+0.4*bybz[3]*grad_u_z[3]+0.4*byby[3]*grad_u_y[3]+0.4*bxby[3]*grad_u_x[3]; 

  bb_grad_u[0] += bb_grad_u_comp[0]; 
  bb_grad_u[1] += bb_grad_u_comp[1]; 
  bb_grad_u[2] += bb_grad_u_comp[2]; 
  bb_grad_u[3] += bb_grad_u_comp[3]; 
  bb_grad_u[4] += bb_grad_u_comp[4]; 
  bb_grad_u[5] += bb_grad_u_comp[5]; 
  bb_grad_u[6] += bb_grad_u_comp[6]; 
  bb_grad_u[7] += bb_grad_u_comp[7]; 
  bb_grad_u[8] += bb_grad_u_comp[8]; 

  p_perp_source[0] += bb_grad_u_comp[0]-1.0*(nu_c[0]+grad_u_y[0]); 
  p_perp_source[1] += bb_grad_u_comp[1]-1.0*(nu_c[1]+grad_u_y[1]); 
  p_perp_source[2] += bb_grad_u_comp[2]-1.0*(nu_c[2]+grad_u_y[2]); 
  p_perp_source[3] += bb_grad_u_comp[3]-1.0*(nu_c[3]+grad_u_y[3]); 
  p_perp_source[4] += bb_grad_u_comp[4]-1.0*nu_c[4]; 
  p_perp_source[5] += bb_grad_u_comp[5]-1.0*(nu_c[5]+grad_u_y[5]); 
  p_perp_source[6] += bb_grad_u_comp[6]-1.0*nu_c[6]; 
  p_perp_source[7] += bb_grad_u_comp[7]-1.0*(nu_c[7]+grad_u_y[7]); 
  p_perp_source[8] += bb_grad_u_comp[8]-1.0*nu_c[8]; 

} 
