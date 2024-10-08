#include <gkyl_maxwell_kernels.h> 
#include <gkyl_basis_ser_2x_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH void em_div_b_x_2x_ser_p1(const double *dxv, 
  const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
  const double *bvar_c, double* GKYL_RESTRICT max_b, double* GKYL_RESTRICT div_b) 
{ 
  // dxv[NDIM]:       Cell spacing.
  // bvar_surf_l/c/r: Input surface magnetic field unit vector in left/center/right cells in each direction. 
  // bvar_c:          Input volume expansion of magnetic field unit vector in center cell. 
  // max_b:           Output surface expansion of max |b| for Lax penalization of streaming: lambda_i = |b_i|. 
  // div_b:           Output volume expansion of div(b).

  const double dx1 = 2.0/dxv[0]; 
  const double *b_c = &bvar_c[0]; 
  const double *b_surf_lr = &bvar_surf_l[2]; 
  const double *b_surf_cl = &bvar_surf_c[0]; 
  const double *b_surf_cr = &bvar_surf_c[2]; 
  const double *b_surf_rl = &bvar_surf_r[0]; 

  double *max_b_l = &max_b[0]; 
  double *max_b_r = &max_b[2]; 

  double bl_r = 0.0; 
  double bc_l = 0.0; 
  double bc_r = 0.0; 
  double br_l = 0.0; 
  double max_b_quad_l[2] = {0.0}; 
  double max_b_quad_r[2] = {0.0}; 

  bl_r = 0.7071067811865475*b_surf_lr[0]-0.7071067811865475*b_surf_lr[1]; 
  bc_l = 0.7071067811865475*b_surf_cl[0]-0.7071067811865475*b_surf_cl[1]; 
  bc_r = 0.7071067811865475*b_surf_cr[0]-0.7071067811865475*b_surf_cr[1]; 
  br_l = 0.7071067811865475*b_surf_rl[0]-0.7071067811865475*b_surf_rl[1]; 
  max_b_quad_l[0] = fmax(fabs(bl_r), fabs(bc_l)); 
  max_b_quad_r[0] = fmax(fabs(bc_r), fabs(br_l)); 

  bl_r = 0.7071067811865475*b_surf_lr[1]+0.7071067811865475*b_surf_lr[0]; 
  bc_l = 0.7071067811865475*b_surf_cl[1]+0.7071067811865475*b_surf_cl[0]; 
  bc_r = 0.7071067811865475*b_surf_cr[1]+0.7071067811865475*b_surf_cr[0]; 
  br_l = 0.7071067811865475*b_surf_rl[1]+0.7071067811865475*b_surf_rl[0]; 
  max_b_quad_l[1] = fmax(fabs(bl_r), fabs(bc_l)); 
  max_b_quad_r[1] = fmax(fabs(bc_r), fabs(br_l)); 

  ser_2x_p1_upwind_quad_to_modal(max_b_quad_l, max_b_l); 
  ser_2x_p1_upwind_quad_to_modal(max_b_quad_r, max_b_r); 

  div_b[0] += (0.3535533905932737*b_surf_rl[0]-0.3535533905932737*b_surf_lr[0]+0.3535533905932737*b_surf_cr[0]-0.3535533905932737*b_surf_cl[0])*dx1; 
  div_b[1] += (0.6123724356957944*(b_surf_rl[0]+b_surf_lr[0]+b_surf_cr[0]+b_surf_cl[0])-1.732050807568877*b_c[0])*dx1; 
  div_b[2] += (0.3535533905932737*b_surf_rl[1]-0.3535533905932737*b_surf_lr[1]+0.3535533905932737*b_surf_cr[1]-0.3535533905932737*b_surf_cl[1])*dx1; 
  div_b[3] += (0.6123724356957944*(b_surf_rl[1]+b_surf_lr[1]+b_surf_cr[1]+b_surf_cl[1])-1.732050807568877*b_c[2])*dx1; 

} 
