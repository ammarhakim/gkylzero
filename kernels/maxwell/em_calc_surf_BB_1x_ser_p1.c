#include <gkyl_maxwell_kernels.h> 
GKYL_CU_DH void em_calc_surf_BB_1x_ser_p1(const double *em, double* GKYL_RESTRICT out_surf) 
{ 
  // em:             Input electromagnetic fields. 
  // out_surf:       Output [BxBx_xl, BxBx_xr, ByBy_xl, ByBy_xr, BzBz_xl, BzBz_xr, 
  //                         BxBx_yl, BxBx_yr, ByBy_yl, ByBy_yr, BzBz_yl, BzBz_yr, 
  //                         BxBx_zl, BxBx_zr, ByBy_zl, ByBy_zr, BzBz_zl, BzBz_zr]. 
 
  double *Bx_sq_xl = &out_surf[0]; 
  double *Bx_sq_xr = &out_surf[1]; 
  double *By_sq_xl = &out_surf[2]; 
  double *By_sq_xr = &out_surf[3]; 
  double *Bz_sq_xl = &out_surf[4]; 
  double *Bz_sq_xr = &out_surf[5]; 
  const double *B_x = &em[6]; 
  const double *B_y = &em[8]; 
  const double *B_z = &em[10]; 
 
  double B_x_l = 0.7071067811865475*B_x[0]-1.224744871391589*B_x[1]; 
  double B_x_r = 1.224744871391589*B_x[1]+0.7071067811865475*B_x[0]; 
  double B_y_l = 0.7071067811865475*B_y[0]-1.224744871391589*B_y[1]; 
  double B_y_r = 1.224744871391589*B_y[1]+0.7071067811865475*B_y[0]; 
  double B_z_l = 0.7071067811865475*B_z[0]-1.224744871391589*B_z[1]; 
  double B_z_r = 1.224744871391589*B_z[1]+0.7071067811865475*B_z[0]; 
  Bx_sq_xl[0] = B_x_l*B_x_l; 
  By_sq_xl[0] = B_y_l*B_y_l; 
  Bz_sq_xl[0] = B_z_l*B_z_l; 
  Bx_sq_xr[0] = B_x_r*B_x_r; 
  By_sq_xr[0] = B_y_r*B_y_r; 
  Bz_sq_xr[0] = B_z_r*B_z_r; 
 
} 
 
