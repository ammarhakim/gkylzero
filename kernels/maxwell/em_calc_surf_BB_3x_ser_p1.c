#include <gkyl_maxwell_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
GKYL_CU_DH void em_calc_surf_BB_3x_ser_p1(const double *em, double* GKYL_RESTRICT out_surf) 
{ 
  // em:             Input electromagnetic fields. 
  // out_surf:       Output [BxBx_xl, BxBx_xr, ByBy_xl, ByBy_xr, BzBz_xl, BzBz_xr, 
  //                         BxBx_yl, BxBx_yr, ByBy_yl, ByBy_yr, BzBz_yl, BzBz_yr, 
  //                         BxBx_zl, BxBx_zr, ByBy_zl, ByBy_zr, BzBz_zl, BzBz_zr]. 
 
  double *Bx_sq_xl = &out_surf[0]; 
  double *Bx_sq_xr = &out_surf[4]; 
  double *By_sq_xl = &out_surf[8]; 
  double *By_sq_xr = &out_surf[12]; 
  double *Bz_sq_xl = &out_surf[16]; 
  double *Bz_sq_xr = &out_surf[20]; 
 
  double *Bx_sq_yl = &out_surf[24]; 
  double *Bx_sq_yr = &out_surf[28]; 
  double *By_sq_yl = &out_surf[32]; 
  double *By_sq_yr = &out_surf[36]; 
  double *Bz_sq_yl = &out_surf[40]; 
  double *Bz_sq_yr = &out_surf[44]; 
 
  double *Bx_sq_zl = &out_surf[48]; 
  double *Bx_sq_zr = &out_surf[52]; 
  double *By_sq_zl = &out_surf[56]; 
  double *By_sq_zr = &out_surf[60]; 
  double *Bz_sq_zl = &out_surf[64]; 
  double *Bz_sq_zr = &out_surf[68]; 
 
  const double *B_x = &em[24]; 
  const double *B_y = &em[32]; 
  const double *B_z = &em[40]; 
 
  double B_x_xl[4] = {0.0}; 
  double B_x_xr[4] = {0.0}; 
  double B_y_xl[4] = {0.0}; 
  double B_y_xr[4] = {0.0}; 
  double B_z_xl[4] = {0.0}; 
  double B_z_xr[4] = {0.0}; 
 
  B_x_xl[0] = 0.7071067811865475*B_x[0]-1.224744871391589*B_x[1]; 
  B_x_xl[1] = 0.7071067811865475*B_x[2]-1.224744871391589*B_x[4]; 
  B_x_xl[2] = 0.7071067811865475*B_x[3]-1.224744871391589*B_x[5]; 
  B_x_xl[3] = 0.7071067811865475*B_x[6]-1.224744871391589*B_x[7]; 
  B_x_xr[0] = 1.224744871391589*B_x[1]+0.7071067811865475*B_x[0]; 
  B_x_xr[1] = 1.224744871391589*B_x[4]+0.7071067811865475*B_x[2]; 
  B_x_xr[2] = 1.224744871391589*B_x[5]+0.7071067811865475*B_x[3]; 
  B_x_xr[3] = 1.224744871391589*B_x[7]+0.7071067811865475*B_x[6]; 
  B_y_xl[0] = 0.7071067811865475*B_y[0]-1.224744871391589*B_y[1]; 
  B_y_xl[1] = 0.7071067811865475*B_y[2]-1.224744871391589*B_y[4]; 
  B_y_xl[2] = 0.7071067811865475*B_y[3]-1.224744871391589*B_y[5]; 
  B_y_xl[3] = 0.7071067811865475*B_y[6]-1.224744871391589*B_y[7]; 
  B_y_xr[0] = 1.224744871391589*B_y[1]+0.7071067811865475*B_y[0]; 
  B_y_xr[1] = 1.224744871391589*B_y[4]+0.7071067811865475*B_y[2]; 
  B_y_xr[2] = 1.224744871391589*B_y[5]+0.7071067811865475*B_y[3]; 
  B_y_xr[3] = 1.224744871391589*B_y[7]+0.7071067811865475*B_y[6]; 
  B_z_xl[0] = 0.7071067811865475*B_z[0]-1.224744871391589*B_z[1]; 
  B_z_xl[1] = 0.7071067811865475*B_z[2]-1.224744871391589*B_z[4]; 
  B_z_xl[2] = 0.7071067811865475*B_z[3]-1.224744871391589*B_z[5]; 
  B_z_xl[3] = 0.7071067811865475*B_z[6]-1.224744871391589*B_z[7]; 
  B_z_xr[0] = 1.224744871391589*B_z[1]+0.7071067811865475*B_z[0]; 
  B_z_xr[1] = 1.224744871391589*B_z[4]+0.7071067811865475*B_z[2]; 
  B_z_xr[2] = 1.224744871391589*B_z[5]+0.7071067811865475*B_z[3]; 
  B_z_xr[3] = 1.224744871391589*B_z[7]+0.7071067811865475*B_z[6]; 
  // Calculate B_i B_i on the left x interface. 
  binop_mul_2d_ser_p1(B_x_xl, B_x_xl, Bx_sq_xl); 
  binop_mul_2d_ser_p1(B_y_xl, B_y_xl, By_sq_xl); 
  binop_mul_2d_ser_p1(B_z_xl, B_z_xl, Bz_sq_xl); 
 
  // Calculate B_i B_i on the right x interface. 
  binop_mul_2d_ser_p1(B_x_xr, B_x_xr, Bx_sq_xr); 
  binop_mul_2d_ser_p1(B_y_xr, B_y_xr, By_sq_xr); 
  binop_mul_2d_ser_p1(B_z_xr, B_z_xr, Bz_sq_xr); 
 
  double B_x_yl[4] = {0.0}; 
  double B_x_yr[4] = {0.0}; 
  double B_y_yl[4] = {0.0}; 
  double B_y_yr[4] = {0.0}; 
  double B_z_yl[4] = {0.0}; 
  double B_z_yr[4] = {0.0}; 
 
  B_x_yl[0] = 0.7071067811865475*B_x[0]-1.224744871391589*B_x[2]; 
  B_x_yl[1] = 0.7071067811865475*B_x[1]-1.224744871391589*B_x[4]; 
  B_x_yl[2] = 0.7071067811865475*B_x[3]-1.224744871391589*B_x[6]; 
  B_x_yl[3] = 0.7071067811865475*B_x[5]-1.224744871391589*B_x[7]; 
  B_x_yr[0] = 1.224744871391589*B_x[2]+0.7071067811865475*B_x[0]; 
  B_x_yr[1] = 1.224744871391589*B_x[4]+0.7071067811865475*B_x[1]; 
  B_x_yr[2] = 1.224744871391589*B_x[6]+0.7071067811865475*B_x[3]; 
  B_x_yr[3] = 1.224744871391589*B_x[7]+0.7071067811865475*B_x[5]; 
  B_y_yl[0] = 0.7071067811865475*B_y[0]-1.224744871391589*B_y[2]; 
  B_y_yl[1] = 0.7071067811865475*B_y[1]-1.224744871391589*B_y[4]; 
  B_y_yl[2] = 0.7071067811865475*B_y[3]-1.224744871391589*B_y[6]; 
  B_y_yl[3] = 0.7071067811865475*B_y[5]-1.224744871391589*B_y[7]; 
  B_y_yr[0] = 1.224744871391589*B_y[2]+0.7071067811865475*B_y[0]; 
  B_y_yr[1] = 1.224744871391589*B_y[4]+0.7071067811865475*B_y[1]; 
  B_y_yr[2] = 1.224744871391589*B_y[6]+0.7071067811865475*B_y[3]; 
  B_y_yr[3] = 1.224744871391589*B_y[7]+0.7071067811865475*B_y[5]; 
  B_z_yl[0] = 0.7071067811865475*B_z[0]-1.224744871391589*B_z[2]; 
  B_z_yl[1] = 0.7071067811865475*B_z[1]-1.224744871391589*B_z[4]; 
  B_z_yl[2] = 0.7071067811865475*B_z[3]-1.224744871391589*B_z[6]; 
  B_z_yl[3] = 0.7071067811865475*B_z[5]-1.224744871391589*B_z[7]; 
  B_z_yr[0] = 1.224744871391589*B_z[2]+0.7071067811865475*B_z[0]; 
  B_z_yr[1] = 1.224744871391589*B_z[4]+0.7071067811865475*B_z[1]; 
  B_z_yr[2] = 1.224744871391589*B_z[6]+0.7071067811865475*B_z[3]; 
  B_z_yr[3] = 1.224744871391589*B_z[7]+0.7071067811865475*B_z[5]; 
  // Calculate B_i B_i on the left y interface. 
  binop_mul_2d_ser_p1(B_x_yl, B_x_yl, Bx_sq_yl); 
  binop_mul_2d_ser_p1(B_y_yl, B_y_yl, By_sq_yl); 
  binop_mul_2d_ser_p1(B_z_yl, B_z_yl, Bz_sq_yl); 
 
  // Calculate B_i B_i on the right y interface. 
  binop_mul_2d_ser_p1(B_x_yr, B_x_yr, Bx_sq_yr); 
  binop_mul_2d_ser_p1(B_y_yr, B_y_yr, By_sq_yr); 
  binop_mul_2d_ser_p1(B_z_yr, B_z_yr, Bz_sq_yr); 
 
  double B_x_zl[4] = {0.0}; 
  double B_x_zr[4] = {0.0}; 
  double B_y_zl[4] = {0.0}; 
  double B_y_zr[4] = {0.0}; 
  double B_z_zl[4] = {0.0}; 
  double B_z_zr[4] = {0.0}; 
 
  B_x_zl[0] = 0.7071067811865475*B_x[0]-1.224744871391589*B_x[3]; 
  B_x_zl[1] = 0.7071067811865475*B_x[1]-1.224744871391589*B_x[5]; 
  B_x_zl[2] = 0.7071067811865475*B_x[2]-1.224744871391589*B_x[6]; 
  B_x_zl[3] = 0.7071067811865475*B_x[4]-1.224744871391589*B_x[7]; 
  B_x_zr[0] = 1.224744871391589*B_x[3]+0.7071067811865475*B_x[0]; 
  B_x_zr[1] = 1.224744871391589*B_x[5]+0.7071067811865475*B_x[1]; 
  B_x_zr[2] = 1.224744871391589*B_x[6]+0.7071067811865475*B_x[2]; 
  B_x_zr[3] = 1.224744871391589*B_x[7]+0.7071067811865475*B_x[4]; 
  B_y_zl[0] = 0.7071067811865475*B_y[0]-1.224744871391589*B_y[3]; 
  B_y_zl[1] = 0.7071067811865475*B_y[1]-1.224744871391589*B_y[5]; 
  B_y_zl[2] = 0.7071067811865475*B_y[2]-1.224744871391589*B_y[6]; 
  B_y_zl[3] = 0.7071067811865475*B_y[4]-1.224744871391589*B_y[7]; 
  B_y_zr[0] = 1.224744871391589*B_y[3]+0.7071067811865475*B_y[0]; 
  B_y_zr[1] = 1.224744871391589*B_y[5]+0.7071067811865475*B_y[1]; 
  B_y_zr[2] = 1.224744871391589*B_y[6]+0.7071067811865475*B_y[2]; 
  B_y_zr[3] = 1.224744871391589*B_y[7]+0.7071067811865475*B_y[4]; 
  B_z_zl[0] = 0.7071067811865475*B_z[0]-1.224744871391589*B_z[3]; 
  B_z_zl[1] = 0.7071067811865475*B_z[1]-1.224744871391589*B_z[5]; 
  B_z_zl[2] = 0.7071067811865475*B_z[2]-1.224744871391589*B_z[6]; 
  B_z_zl[3] = 0.7071067811865475*B_z[4]-1.224744871391589*B_z[7]; 
  B_z_zr[0] = 1.224744871391589*B_z[3]+0.7071067811865475*B_z[0]; 
  B_z_zr[1] = 1.224744871391589*B_z[5]+0.7071067811865475*B_z[1]; 
  B_z_zr[2] = 1.224744871391589*B_z[6]+0.7071067811865475*B_z[2]; 
  B_z_zr[3] = 1.224744871391589*B_z[7]+0.7071067811865475*B_z[4]; 
  // Calculate B_i B_i on the left z interface. 
  binop_mul_2d_ser_p1(B_x_zl, B_x_zl, Bx_sq_zl); 
  binop_mul_2d_ser_p1(B_y_zl, B_y_zl, By_sq_zl); 
  binop_mul_2d_ser_p1(B_z_zl, B_z_zl, Bz_sq_zl); 
 
  // Calculate B_i B_i on the right z interface. 
  binop_mul_2d_ser_p1(B_x_zr, B_x_zr, Bx_sq_zr); 
  binop_mul_2d_ser_p1(B_y_zr, B_y_zr, By_sq_zr); 
  binop_mul_2d_ser_p1(B_z_zr, B_z_zr, Bz_sq_zr); 
 
} 
 
