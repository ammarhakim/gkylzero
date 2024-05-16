#include <gkyl_maxwell_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
GKYL_CU_DH void em_calc_surf_BB_2x_tensor_p2(const double *em, double* GKYL_RESTRICT out_surf) 
{ 
  // em:             Input electromagnetic fields. 
  // out_surf:       Output [BxBx_xl, BxBx_xr, ByBy_xl, ByBy_xr, BzBz_xl, BzBz_xr, 
  //                         BxBx_yl, BxBx_yr, ByBy_yl, ByBy_yr, BzBz_yl, BzBz_yr, 
  //                         BxBx_zl, BxBx_zr, ByBy_zl, ByBy_zr, BzBz_zl, BzBz_zr]. 
 
  double *Bx_sq_xl = &out_surf[0]; 
  double *Bx_sq_xr = &out_surf[3]; 
  double *By_sq_xl = &out_surf[6]; 
  double *By_sq_xr = &out_surf[9]; 
  double *Bz_sq_xl = &out_surf[12]; 
  double *Bz_sq_xr = &out_surf[15]; 
 
  double *Bx_sq_yl = &out_surf[18]; 
  double *Bx_sq_yr = &out_surf[21]; 
  double *By_sq_yl = &out_surf[24]; 
  double *By_sq_yr = &out_surf[27]; 
  double *Bz_sq_yl = &out_surf[30]; 
  double *Bz_sq_yr = &out_surf[33]; 
 
  const double *B_x = &em[27]; 
  const double *B_y = &em[36]; 
  const double *B_z = &em[45]; 
 
  double B_x_xl[3] = {0.0}; 
  double B_x_xr[3] = {0.0}; 
  double B_y_xl[3] = {0.0}; 
  double B_y_xr[3] = {0.0}; 
  double B_z_xl[3] = {0.0}; 
  double B_z_xr[3] = {0.0}; 
 
  B_x_xl[0] = 1.58113883008419*B_x[4]-1.224744871391589*B_x[1]+0.7071067811865475*B_x[0]; 
  B_x_xl[1] = 1.58113883008419*B_x[6]-1.224744871391589*B_x[3]+0.7071067811865475*B_x[2]; 
  B_x_xl[2] = 1.58113883008419*B_x[8]-1.224744871391589*B_x[7]+0.7071067811865475*B_x[5]; 
  B_x_xr[0] = 1.58113883008419*B_x[4]+1.224744871391589*B_x[1]+0.7071067811865475*B_x[0]; 
  B_x_xr[1] = 1.58113883008419*B_x[6]+1.224744871391589*B_x[3]+0.7071067811865475*B_x[2]; 
  B_x_xr[2] = 1.58113883008419*B_x[8]+1.224744871391589*B_x[7]+0.7071067811865475*B_x[5]; 
  B_y_xl[0] = 1.58113883008419*B_y[4]-1.224744871391589*B_y[1]+0.7071067811865475*B_y[0]; 
  B_y_xl[1] = 1.58113883008419*B_y[6]-1.224744871391589*B_y[3]+0.7071067811865475*B_y[2]; 
  B_y_xl[2] = 1.58113883008419*B_y[8]-1.224744871391589*B_y[7]+0.7071067811865475*B_y[5]; 
  B_y_xr[0] = 1.58113883008419*B_y[4]+1.224744871391589*B_y[1]+0.7071067811865475*B_y[0]; 
  B_y_xr[1] = 1.58113883008419*B_y[6]+1.224744871391589*B_y[3]+0.7071067811865475*B_y[2]; 
  B_y_xr[2] = 1.58113883008419*B_y[8]+1.224744871391589*B_y[7]+0.7071067811865475*B_y[5]; 
  B_z_xl[0] = 1.58113883008419*B_z[4]-1.224744871391589*B_z[1]+0.7071067811865475*B_z[0]; 
  B_z_xl[1] = 1.58113883008419*B_z[6]-1.224744871391589*B_z[3]+0.7071067811865475*B_z[2]; 
  B_z_xl[2] = 1.58113883008419*B_z[8]-1.224744871391589*B_z[7]+0.7071067811865475*B_z[5]; 
  B_z_xr[0] = 1.58113883008419*B_z[4]+1.224744871391589*B_z[1]+0.7071067811865475*B_z[0]; 
  B_z_xr[1] = 1.58113883008419*B_z[6]+1.224744871391589*B_z[3]+0.7071067811865475*B_z[2]; 
  B_z_xr[2] = 1.58113883008419*B_z[8]+1.224744871391589*B_z[7]+0.7071067811865475*B_z[5]; 
  // Calculate B_i B_i on the left x interface. 
  binop_mul_1d_ser_p2(B_x_xl, B_x_xl, Bx_sq_xl); 
  binop_mul_1d_ser_p2(B_y_xl, B_y_xl, By_sq_xl); 
  binop_mul_1d_ser_p2(B_z_xl, B_z_xl, Bz_sq_xl); 
 
  // Calculate B_i B_i on the right x interface. 
  binop_mul_1d_ser_p2(B_x_xr, B_x_xr, Bx_sq_xr); 
  binop_mul_1d_ser_p2(B_y_xr, B_y_xr, By_sq_xr); 
  binop_mul_1d_ser_p2(B_z_xr, B_z_xr, Bz_sq_xr); 
 
  double B_x_yl[3] = {0.0}; 
  double B_x_yr[3] = {0.0}; 
  double B_y_yl[3] = {0.0}; 
  double B_y_yr[3] = {0.0}; 
  double B_z_yl[3] = {0.0}; 
  double B_z_yr[3] = {0.0}; 
 
  B_x_yl[0] = 1.58113883008419*B_x[5]-1.224744871391589*B_x[2]+0.7071067811865475*B_x[0]; 
  B_x_yl[1] = 1.58113883008419*B_x[7]-1.224744871391589*B_x[3]+0.7071067811865475*B_x[1]; 
  B_x_yl[2] = 1.58113883008419*B_x[8]-1.224744871391589*B_x[6]+0.7071067811865475*B_x[4]; 
  B_x_yr[0] = 1.58113883008419*B_x[5]+1.224744871391589*B_x[2]+0.7071067811865475*B_x[0]; 
  B_x_yr[1] = 1.58113883008419*B_x[7]+1.224744871391589*B_x[3]+0.7071067811865475*B_x[1]; 
  B_x_yr[2] = 1.58113883008419*B_x[8]+1.224744871391589*B_x[6]+0.7071067811865475*B_x[4]; 
  B_y_yl[0] = 1.58113883008419*B_y[5]-1.224744871391589*B_y[2]+0.7071067811865475*B_y[0]; 
  B_y_yl[1] = 1.58113883008419*B_y[7]-1.224744871391589*B_y[3]+0.7071067811865475*B_y[1]; 
  B_y_yl[2] = 1.58113883008419*B_y[8]-1.224744871391589*B_y[6]+0.7071067811865475*B_y[4]; 
  B_y_yr[0] = 1.58113883008419*B_y[5]+1.224744871391589*B_y[2]+0.7071067811865475*B_y[0]; 
  B_y_yr[1] = 1.58113883008419*B_y[7]+1.224744871391589*B_y[3]+0.7071067811865475*B_y[1]; 
  B_y_yr[2] = 1.58113883008419*B_y[8]+1.224744871391589*B_y[6]+0.7071067811865475*B_y[4]; 
  B_z_yl[0] = 1.58113883008419*B_z[5]-1.224744871391589*B_z[2]+0.7071067811865475*B_z[0]; 
  B_z_yl[1] = 1.58113883008419*B_z[7]-1.224744871391589*B_z[3]+0.7071067811865475*B_z[1]; 
  B_z_yl[2] = 1.58113883008419*B_z[8]-1.224744871391589*B_z[6]+0.7071067811865475*B_z[4]; 
  B_z_yr[0] = 1.58113883008419*B_z[5]+1.224744871391589*B_z[2]+0.7071067811865475*B_z[0]; 
  B_z_yr[1] = 1.58113883008419*B_z[7]+1.224744871391589*B_z[3]+0.7071067811865475*B_z[1]; 
  B_z_yr[2] = 1.58113883008419*B_z[8]+1.224744871391589*B_z[6]+0.7071067811865475*B_z[4]; 
  // Calculate B_i B_i on the left y interface. 
  binop_mul_1d_ser_p2(B_x_yl, B_x_yl, Bx_sq_yl); 
  binop_mul_1d_ser_p2(B_y_yl, B_y_yl, By_sq_yl); 
  binop_mul_1d_ser_p2(B_z_yl, B_z_yl, Bz_sq_yl); 
 
  // Calculate B_i B_i on the right y interface. 
  binop_mul_1d_ser_p2(B_x_yr, B_x_yr, Bx_sq_yr); 
  binop_mul_1d_ser_p2(B_y_yr, B_y_yr, By_sq_yr); 
  binop_mul_1d_ser_p2(B_z_yr, B_z_yr, Bz_sq_yr); 
 
} 
 
