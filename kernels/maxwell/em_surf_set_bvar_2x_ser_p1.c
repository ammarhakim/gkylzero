#include <gkyl_mat.h> 
#include <gkyl_maxwell_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
#include <gkyl_basis_ser_1x_p1_inv.h> 
GKYL_CU_DH void em_surf_set_bvar_2x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *BB_surf, int* cell_avg_magB2_surf) 
{ 
  // count:   integer to indicate which matrix being fetched. 
  // A:       preallocated LHS matrix. 
  // rhs:     preallocated RHS vector. 
  // BB_surf: Surface B_i B_j [BxBx_xl, BxBx_xr, ByBy_xl, ByBy_xr, BzBz_xl, BzBz_xr, 
  //                           BxBx_yl, BxBx_yr, ByBy_yl, ByBy_yr, BzBz_yl, BzBz_yr,  
  //                           BxBx_zl, BxBx_zr, ByBy_zl, ByBy_zr, BzBz_zl, BzBz_zr]. 
  // cell_avg_magB2_surf:      Output flag for cell average if 1/|B|^2 at a surface only used cell averages. 

  struct gkyl_mat rhs_bxbx_xl = gkyl_nmat_get(rhs, count); 
  struct gkyl_mat rhs_bxbx_xr = gkyl_nmat_get(rhs, count+1); 
  gkyl_mat_clear(&rhs_bxbx_xl, 0.0); 
  gkyl_mat_clear(&rhs_bxbx_xr, 0.0); 
  const double *Bx_sq_xl = &BB_surf[0]; 
  const double *Bx_sq_xr = &BB_surf[2]; 
  const double *By_sq_xl = &BB_surf[4]; 
  const double *By_sq_xr = &BB_surf[6]; 
  const double *Bz_sq_xl = &BB_surf[8]; 
  const double *Bz_sq_xr = &BB_surf[10]; 
  int *cell_avg_magB2_xl = &cell_avg_magB2_surf[0]; 
  int *cell_avg_magB2_xr = &cell_avg_magB2_surf[1]; 
 
  struct gkyl_mat rhs_byby_yl = gkyl_nmat_get(rhs, count+2); 
  struct gkyl_mat rhs_byby_yr = gkyl_nmat_get(rhs, count+3); 
  gkyl_mat_clear(&rhs_byby_yl, 0.0); 
  gkyl_mat_clear(&rhs_byby_yr, 0.0); 
  const double *Bx_sq_yl = &BB_surf[12]; 
  const double *Bx_sq_yr = &BB_surf[14]; 
  const double *By_sq_yl = &BB_surf[16]; 
  const double *By_sq_yr = &BB_surf[18]; 
  const double *Bz_sq_yl = &BB_surf[20]; 
  const double *Bz_sq_yr = &BB_surf[22]; 
  int *cell_avg_magB2_yl = &cell_avg_magB2_surf[2]; 
  int *cell_avg_magB2_yr = &cell_avg_magB2_surf[3]; 
 
  double magB2_xl[2] = {0.0}; 
  double magB2_xr[2] = {0.0}; 
  double magB2_yl[2] = {0.0}; 
  double magB2_yr[2] = {0.0}; 
  magB2_xl[0] = Bx_sq_xl[0] + By_sq_xl[0] + Bz_sq_xl[0]; 
  magB2_xr[0] = Bx_sq_xr[0] + By_sq_xr[0] + Bz_sq_xr[0]; 
  magB2_yl[0] = Bx_sq_yl[0] + By_sq_yl[0] + Bz_sq_yl[0]; 
  magB2_yr[0] = Bx_sq_yr[0] + By_sq_yr[0] + Bz_sq_yr[0]; 
  magB2_xl[1] = Bx_sq_xl[1] + By_sq_xl[1] + Bz_sq_xl[1]; 
  magB2_xr[1] = Bx_sq_xr[1] + By_sq_xr[1] + Bz_sq_xr[1]; 
  magB2_yl[1] = Bx_sq_yl[1] + By_sq_yl[1] + Bz_sq_yl[1]; 
  magB2_yr[1] = Bx_sq_yr[1] + By_sq_yr[1] + Bz_sq_yr[1]; 
  // If |B|^2 < 0 at control points along a surface, only use cell average to get 1/|B|^2. 
  // Each surface is checked independently. 
  int cell_avg_xl = 0;
  double bxbx_xl[2] = {0.0}; 
 
  int cell_avg_xr = 0;
  double bxbx_xr[2] = {0.0}; 
 
  int cell_avg_yl = 0;
  double byby_yl[2] = {0.0}; 
 
  int cell_avg_yr = 0;
  double byby_yr[2] = {0.0}; 
 
  if (0.7071067811865475*magB2_xl[0]-1.224744871391589*magB2_xl[1] < 0.0) cell_avg_xl = 1; 
  if (0.7071067811865475*magB2_xr[0]-1.224744871391589*magB2_xr[1] < 0.0) cell_avg_xr = 1; 
  if (0.7071067811865475*magB2_yl[0]-1.224744871391589*magB2_yl[1] < 0.0) cell_avg_yl = 1; 
  if (0.7071067811865475*magB2_yr[0]-1.224744871391589*magB2_yr[1] < 0.0) cell_avg_yr = 1; 
  if (1.224744871391589*magB2_xl[1]+0.7071067811865475*magB2_xl[0] < 0.0) cell_avg_xl = 1; 
  if (1.224744871391589*magB2_xr[1]+0.7071067811865475*magB2_xr[0] < 0.0) cell_avg_xr = 1; 
  if (1.224744871391589*magB2_yl[1]+0.7071067811865475*magB2_yl[0] < 0.0) cell_avg_yl = 1; 
  if (1.224744871391589*magB2_yr[1]+0.7071067811865475*magB2_yr[0] < 0.0) cell_avg_yr = 1; 
 
  cell_avg_magB2_xl[0] = cell_avg_xl; 
  cell_avg_magB2_xr[0] = cell_avg_xr; 
  cell_avg_magB2_yl[0] = cell_avg_yl; 
  cell_avg_magB2_yr[0] = cell_avg_yr; 
 
  double magB2_inv_xl[2] = {0.0}; 

  if (cell_avg_xl) { 
  magB2_inv_xl[0] = 2.0/magB2_xl[0]; 
  } else { 
  ser_1x_p1_inv(magB2_xl, magB2_inv_xl); 
  } 
  binop_mul_1d_ser_p1(magB2_inv_xl, Bx_sq_xl, bxbx_xl); 
 
  double magB2_inv_xr[2] = {0.0}; 

  if (cell_avg_xr) { 
  magB2_inv_xr[0] = 2.0/magB2_xr[0]; 
  } else { 
  ser_1x_p1_inv(magB2_xr, magB2_inv_xr); 
  } 
  binop_mul_1d_ser_p1(magB2_inv_xr, Bx_sq_xr, bxbx_xr); 
 
  double magB2_inv_yl[2] = {0.0}; 

  if (cell_avg_yl) { 
  magB2_inv_yl[0] = 2.0/magB2_yl[0]; 
  } else { 
  ser_1x_p1_inv(magB2_yl, magB2_inv_yl); 
  } 
  binop_mul_1d_ser_p1(magB2_inv_yl, By_sq_yl, byby_yl); 
 
  double magB2_inv_yr[2] = {0.0}; 

  if (cell_avg_yr) { 
  magB2_inv_yr[0] = 2.0/magB2_yr[0]; 
  } else { 
  ser_1x_p1_inv(magB2_yr, magB2_inv_yr); 
  } 
  binop_mul_1d_ser_p1(magB2_inv_yr, By_sq_yr, byby_yr); 
 
  gkyl_mat_set(&rhs_bxbx_xl,0,0,bxbx_xl[0]); 
  gkyl_mat_set(&rhs_bxbx_xr,0,0,bxbx_xr[0]); 
 
  gkyl_mat_set(&rhs_byby_yl,0,0,byby_yl[0]); 
  gkyl_mat_set(&rhs_byby_yr,0,0,byby_yr[0]); 
 
  gkyl_mat_set(&rhs_bxbx_xl,1,0,bxbx_xl[1]); 
  gkyl_mat_set(&rhs_bxbx_xr,1,0,bxbx_xr[1]); 
 
  gkyl_mat_set(&rhs_byby_yl,1,0,byby_yl[1]); 
  gkyl_mat_set(&rhs_byby_yr,1,0,byby_yr[1]); 
 
} 
