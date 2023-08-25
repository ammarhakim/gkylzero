#include <gkyl_mat.h> 
#include <gkyl_maxwell_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
GKYL_CU_DH void em_surf_set_bvar_2x_tensor_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *BB_surf, int* cell_avg_magB2_surf) 
{ 
  // count:   integer to indicate which matrix being fetched. 
  // A:       preallocated LHS matrix. 
  // rhs:     preallocated RHS vector. 
  // BB_surf: Surface B_i B_j [BxBx_xl, BxBx_xr, ByBy_xl, ByBy_xr, BzBz_xl, BzBz_xr, BxBy_xl, BxBy_xr, BxBz_xl, BxBz_xr, 
  //                           BxBx_yl, BxBx_yr, ByBy_yl, ByBy_yr, BzBz_yl, BzBz_yr, BxBy_yl, BxBy_yr, ByBz_yl, ByBz_yr, 
  //                           BxBx_zl, BxBx_zr, ByBy_zl, ByBy_zr, BzBz_zl, BzBz_zr, BxBz_zl, BxBz_zr, ByBz_zl, ByBz_zr]. 
  // cell_avg_magB2_surf:      Output flag for cell average if 1/|B|^2 at a surface only used cell averages. 

  struct gkyl_mat A_bxbx_xl = gkyl_nmat_get(A, count); 
  struct gkyl_mat A_bxbx_xr = gkyl_nmat_get(A, count+1); 
  struct gkyl_mat A_bxby_xl = gkyl_nmat_get(A, count+2); 
  struct gkyl_mat A_bxby_xr = gkyl_nmat_get(A, count+3); 
  struct gkyl_mat A_bxbz_xl = gkyl_nmat_get(A, count+4); 
  struct gkyl_mat A_bxbz_xr = gkyl_nmat_get(A, count+5); 
  struct gkyl_mat rhs_bxbx_xl = gkyl_nmat_get(rhs, count); 
  struct gkyl_mat rhs_bxbx_xr = gkyl_nmat_get(rhs, count+1); 
  struct gkyl_mat rhs_bxby_xl = gkyl_nmat_get(rhs, count+2); 
  struct gkyl_mat rhs_bxby_xr = gkyl_nmat_get(rhs, count+3); 
  struct gkyl_mat rhs_bxbz_xl = gkyl_nmat_get(rhs, count+4); 
  struct gkyl_mat rhs_bxbz_xr = gkyl_nmat_get(rhs, count+5); 
  gkyl_mat_clear(&A_bxbx_xl, 0.0); gkyl_mat_clear(&rhs_bxbx_xl, 0.0); 
  gkyl_mat_clear(&A_bxbx_xr, 0.0); gkyl_mat_clear(&rhs_bxbx_xr, 0.0); 
  gkyl_mat_clear(&A_bxby_xl, 0.0); gkyl_mat_clear(&rhs_bxby_xl, 0.0); 
  gkyl_mat_clear(&A_bxby_xr, 0.0); gkyl_mat_clear(&rhs_bxby_xr, 0.0); 
  gkyl_mat_clear(&A_bxbz_xl, 0.0); gkyl_mat_clear(&rhs_bxbz_xl, 0.0); 
  gkyl_mat_clear(&A_bxbz_xr, 0.0); gkyl_mat_clear(&rhs_bxbz_xr, 0.0); 
  const double *Bx_sq_xl = &BB_surf[0]; 
  const double *Bx_sq_xr = &BB_surf[3]; 
  const double *By_sq_xl = &BB_surf[6]; 
  const double *By_sq_xr = &BB_surf[9]; 
  const double *Bz_sq_xl = &BB_surf[12]; 
  const double *Bz_sq_xr = &BB_surf[15]; 
  const double *B_x_B_y_xl = &BB_surf[18]; 
  const double *B_x_B_y_xr = &BB_surf[21]; 
  const double *B_x_B_z_xl = &BB_surf[24]; 
  const double *B_x_B_z_xr = &BB_surf[27]; 
  int *cell_avg_magB2_xl = &cell_avg_magB2_surf[0]; 
  int *cell_avg_magB2_xr = &cell_avg_magB2_surf[1]; 
 
  struct gkyl_mat A_byby_yl = gkyl_nmat_get(A, count+6); 
  struct gkyl_mat A_byby_yr = gkyl_nmat_get(A, count+7); 
  struct gkyl_mat A_bxby_yl = gkyl_nmat_get(A, count+8); 
  struct gkyl_mat A_bxby_yr = gkyl_nmat_get(A, count+9); 
  struct gkyl_mat A_bybz_yl = gkyl_nmat_get(A, count+10); 
  struct gkyl_mat A_bybz_yr = gkyl_nmat_get(A, count+11); 
  struct gkyl_mat rhs_byby_yl = gkyl_nmat_get(rhs, count+6); 
  struct gkyl_mat rhs_byby_yr = gkyl_nmat_get(rhs, count+7); 
  struct gkyl_mat rhs_bxby_yl = gkyl_nmat_get(rhs, count+8); 
  struct gkyl_mat rhs_bxby_yr = gkyl_nmat_get(rhs, count+9); 
  struct gkyl_mat rhs_bybz_yl = gkyl_nmat_get(rhs, count+10); 
  struct gkyl_mat rhs_bybz_yr = gkyl_nmat_get(rhs, count+11); 
  gkyl_mat_clear(&A_byby_yl, 0.0); gkyl_mat_clear(&rhs_byby_yl, 0.0); 
  gkyl_mat_clear(&A_byby_yr, 0.0); gkyl_mat_clear(&rhs_byby_yr, 0.0); 
  gkyl_mat_clear(&A_bxby_yl, 0.0); gkyl_mat_clear(&rhs_bxby_yl, 0.0); 
  gkyl_mat_clear(&A_bxby_yr, 0.0); gkyl_mat_clear(&rhs_bxby_yr, 0.0); 
  gkyl_mat_clear(&A_bybz_yl, 0.0); gkyl_mat_clear(&rhs_bybz_yl, 0.0); 
  gkyl_mat_clear(&A_bybz_yr, 0.0); gkyl_mat_clear(&rhs_bybz_yr, 0.0); 
  const double *Bx_sq_yl = &BB_surf[30]; 
  const double *Bx_sq_yr = &BB_surf[33]; 
  const double *By_sq_yl = &BB_surf[36]; 
  const double *By_sq_yr = &BB_surf[39]; 
  const double *Bz_sq_yl = &BB_surf[42]; 
  const double *Bz_sq_yr = &BB_surf[45]; 
  const double *B_x_B_y_yl = &BB_surf[48]; 
  const double *B_x_B_y_yr = &BB_surf[51]; 
  const double *B_y_B_z_yl = &BB_surf[54]; 
  const double *B_y_B_z_yr = &BB_surf[57]; 
  int *cell_avg_magB2_yl = &cell_avg_magB2_surf[2]; 
  int *cell_avg_magB2_yr = &cell_avg_magB2_surf[3]; 
 
  double magB2_xl[3] = {0.0}; 
  double magB2_xr[3] = {0.0}; 
  double magB2_yl[3] = {0.0}; 
  double magB2_yr[3] = {0.0}; 
  magB2_xl[0] = Bx_sq_xl[0] + By_sq_xl[0] + Bz_sq_xl[0]; 
  magB2_xr[0] = Bx_sq_xr[0] + By_sq_xr[0] + Bz_sq_xr[0]; 
  magB2_yl[0] = Bx_sq_yl[0] + By_sq_yl[0] + Bz_sq_yl[0]; 
  magB2_yr[0] = Bx_sq_yr[0] + By_sq_yr[0] + Bz_sq_yr[0]; 
  magB2_xl[1] = Bx_sq_xl[1] + By_sq_xl[1] + Bz_sq_xl[1]; 
  magB2_xr[1] = Bx_sq_xr[1] + By_sq_xr[1] + Bz_sq_xr[1]; 
  magB2_yl[1] = Bx_sq_yl[1] + By_sq_yl[1] + Bz_sq_yl[1]; 
  magB2_yr[1] = Bx_sq_yr[1] + By_sq_yr[1] + Bz_sq_yr[1]; 
  magB2_xl[2] = Bx_sq_xl[2] + By_sq_xl[2] + Bz_sq_xl[2]; 
  magB2_xr[2] = Bx_sq_xr[2] + By_sq_xr[2] + Bz_sq_xr[2]; 
  magB2_yl[2] = Bx_sq_yl[2] + By_sq_yl[2] + Bz_sq_yl[2]; 
  magB2_yr[2] = Bx_sq_yr[2] + By_sq_yr[2] + Bz_sq_yr[2]; 
  // If |B|^2 < 0 at control points along a surface, only use cell average to get 1/|B|^2. 
  // Each surface is checked independently. 
  int cell_avg_xl = 0;
  int cell_avg_xr = 0;
  int cell_avg_yl = 0;
  int cell_avg_yr = 0;
 
  if (1.58113883008419*magB2_xl[2]-1.224744871391589*magB2_xl[1]+0.7071067811865475*magB2_xl[0] < 0.0) cell_avg_xl = 1; 
  if (1.58113883008419*magB2_xr[2]-1.224744871391589*magB2_xr[1]+0.7071067811865475*magB2_xr[0] < 0.0) cell_avg_xr = 1; 
  if (1.58113883008419*magB2_yl[2]-1.224744871391589*magB2_yl[1]+0.7071067811865475*magB2_yl[0] < 0.0) cell_avg_yl = 1; 
  if (1.58113883008419*magB2_yr[2]-1.224744871391589*magB2_yr[1]+0.7071067811865475*magB2_yr[0] < 0.0) cell_avg_yr = 1; 
  if (0.7071067811865475*magB2_xl[0]-0.7905694150420947*magB2_xl[2] < 0.0) cell_avg_xl = 1; 
  if (0.7071067811865475*magB2_xr[0]-0.7905694150420947*magB2_xr[2] < 0.0) cell_avg_xr = 1; 
  if (0.7071067811865475*magB2_yl[0]-0.7905694150420947*magB2_yl[2] < 0.0) cell_avg_yl = 1; 
  if (0.7071067811865475*magB2_yr[0]-0.7905694150420947*magB2_yr[2] < 0.0) cell_avg_yr = 1; 
  if (1.58113883008419*magB2_xl[2]+1.224744871391589*magB2_xl[1]+0.7071067811865475*magB2_xl[0] < 0.0) cell_avg_xl = 1; 
  if (1.58113883008419*magB2_xr[2]+1.224744871391589*magB2_xr[1]+0.7071067811865475*magB2_xr[0] < 0.0) cell_avg_xr = 1; 
  if (1.58113883008419*magB2_yl[2]+1.224744871391589*magB2_yl[1]+0.7071067811865475*magB2_yl[0] < 0.0) cell_avg_yl = 1; 
  if (1.58113883008419*magB2_yr[2]+1.224744871391589*magB2_yr[1]+0.7071067811865475*magB2_yr[0] < 0.0) cell_avg_yr = 1; 
 
  cell_avg_magB2_xl[0] = cell_avg_xl; 
  cell_avg_magB2_xr[0] = cell_avg_xr; 
  cell_avg_magB2_yl[0] = cell_avg_yl; 
  cell_avg_magB2_yr[0] = cell_avg_yr; 
 
  if (cell_avg_xl) { 
  magB2_xl[1] = 0.0; 
  magB2_xl[2] = 0.0; 
  } 
 
  if (cell_avg_xr) { 
  magB2_xr[1] = 0.0; 
  magB2_xr[2] = 0.0; 
  } 
 
  if (cell_avg_yl) { 
  magB2_yl[1] = 0.0; 
  magB2_yl[2] = 0.0; 
  } 
 
  if (cell_avg_yr) { 
  magB2_yr[1] = 0.0; 
  magB2_yr[2] = 0.0; 
  } 
 
  gkyl_mat_set(&rhs_bxbx_xl,0,0,Bx_sq_xl[0]); 
  gkyl_mat_set(&rhs_bxbx_xr,0,0,Bx_sq_xr[0]); 
  gkyl_mat_set(&rhs_bxby_xl,0,0,B_x_B_y_xl[0]); 
  gkyl_mat_set(&rhs_bxby_xr,0,0,B_x_B_y_xr[0]); 
  gkyl_mat_set(&rhs_bxbz_xl,0,0,B_x_B_z_xl[0]); 
  gkyl_mat_set(&rhs_bxbz_xr,0,0,B_x_B_z_xr[0]); 
 
  gkyl_mat_set(&rhs_byby_yl,0,0,By_sq_yl[0]); 
  gkyl_mat_set(&rhs_byby_yr,0,0,By_sq_yr[0]); 
  gkyl_mat_set(&rhs_bxby_yl,0,0,B_x_B_y_yl[0]); 
  gkyl_mat_set(&rhs_bxby_yr,0,0,B_x_B_y_yr[0]); 
  gkyl_mat_set(&rhs_bybz_yl,0,0,B_y_B_z_yl[0]); 
  gkyl_mat_set(&rhs_bybz_yr,0,0,B_y_B_z_yr[0]); 
 
  gkyl_mat_set(&rhs_bxbx_xl,1,0,Bx_sq_xl[1]); 
  gkyl_mat_set(&rhs_bxbx_xr,1,0,Bx_sq_xr[1]); 
  gkyl_mat_set(&rhs_bxby_xl,1,0,B_x_B_y_xl[1]); 
  gkyl_mat_set(&rhs_bxby_xr,1,0,B_x_B_y_xr[1]); 
  gkyl_mat_set(&rhs_bxbz_xl,1,0,B_x_B_z_xl[1]); 
  gkyl_mat_set(&rhs_bxbz_xr,1,0,B_x_B_z_xr[1]); 
 
  gkyl_mat_set(&rhs_byby_yl,1,0,By_sq_yl[1]); 
  gkyl_mat_set(&rhs_byby_yr,1,0,By_sq_yr[1]); 
  gkyl_mat_set(&rhs_bxby_yl,1,0,B_x_B_y_yl[1]); 
  gkyl_mat_set(&rhs_bxby_yr,1,0,B_x_B_y_yr[1]); 
  gkyl_mat_set(&rhs_bybz_yl,1,0,B_y_B_z_yl[1]); 
  gkyl_mat_set(&rhs_bybz_yr,1,0,B_y_B_z_yr[1]); 
 
  gkyl_mat_set(&rhs_bxbx_xl,2,0,Bx_sq_xl[2]); 
  gkyl_mat_set(&rhs_bxbx_xr,2,0,Bx_sq_xr[2]); 
  gkyl_mat_set(&rhs_bxby_xl,2,0,B_x_B_y_xl[2]); 
  gkyl_mat_set(&rhs_bxby_xr,2,0,B_x_B_y_xr[2]); 
  gkyl_mat_set(&rhs_bxbz_xl,2,0,B_x_B_z_xl[2]); 
  gkyl_mat_set(&rhs_bxbz_xr,2,0,B_x_B_z_xr[2]); 
 
  gkyl_mat_set(&rhs_byby_yl,2,0,By_sq_yl[2]); 
  gkyl_mat_set(&rhs_byby_yr,2,0,By_sq_yr[2]); 
  gkyl_mat_set(&rhs_bxby_yl,2,0,B_x_B_y_yl[2]); 
  gkyl_mat_set(&rhs_bxby_yr,2,0,B_x_B_y_yr[2]); 
  gkyl_mat_set(&rhs_bybz_yl,2,0,B_y_B_z_yl[2]); 
  gkyl_mat_set(&rhs_bybz_yr,2,0,B_y_B_z_yr[2]); 
 
  double temp_magB2_xl = 0.0; 
  double temp_magB2_xr = 0.0; 
  double temp_magB2_yl = 0.0; 
  double temp_magB2_yr = 0.0; 
  temp_magB2_xl = 0.7071067811865475*magB2_xl[0]; 
  gkyl_mat_set(&A_bxbx_xl,0,0,temp_magB2_xl); 
  gkyl_mat_set(&A_bxby_xl,0,0,temp_magB2_xl); 
  gkyl_mat_set(&A_bxbz_xl,0,0,temp_magB2_xl); 
 
  temp_magB2_xr = 0.7071067811865475*magB2_xr[0]; 
  gkyl_mat_set(&A_bxbx_xr,0,0,temp_magB2_xr); 
  gkyl_mat_set(&A_bxby_xr,0,0,temp_magB2_xr); 
  gkyl_mat_set(&A_bxbz_xr,0,0,temp_magB2_xr); 
 
  temp_magB2_yl = 0.7071067811865475*magB2_yl[0]; 
  gkyl_mat_set(&A_byby_yl,0,0,temp_magB2_yl); 
  gkyl_mat_set(&A_bxby_yl,0,0,temp_magB2_yl); 
  gkyl_mat_set(&A_bybz_yl,0,0,temp_magB2_yl); 
 
  temp_magB2_yr = 0.7071067811865475*magB2_yr[0]; 
  gkyl_mat_set(&A_byby_yr,0,0,temp_magB2_yr); 
  gkyl_mat_set(&A_bxby_yr,0,0,temp_magB2_yr); 
  gkyl_mat_set(&A_bybz_yr,0,0,temp_magB2_yr); 
 
  temp_magB2_xl = 0.7071067811865475*magB2_xl[1]; 
  gkyl_mat_set(&A_bxbx_xl,0,1,temp_magB2_xl); 
  gkyl_mat_set(&A_bxby_xl,0,1,temp_magB2_xl); 
  gkyl_mat_set(&A_bxbz_xl,0,1,temp_magB2_xl); 
 
  temp_magB2_xr = 0.7071067811865475*magB2_xr[1]; 
  gkyl_mat_set(&A_bxbx_xr,0,1,temp_magB2_xr); 
  gkyl_mat_set(&A_bxby_xr,0,1,temp_magB2_xr); 
  gkyl_mat_set(&A_bxbz_xr,0,1,temp_magB2_xr); 
 
  temp_magB2_yl = 0.7071067811865475*magB2_yl[1]; 
  gkyl_mat_set(&A_byby_yl,0,1,temp_magB2_yl); 
  gkyl_mat_set(&A_bxby_yl,0,1,temp_magB2_yl); 
  gkyl_mat_set(&A_bybz_yl,0,1,temp_magB2_yl); 
 
  temp_magB2_yr = 0.7071067811865475*magB2_yr[1]; 
  gkyl_mat_set(&A_byby_yr,0,1,temp_magB2_yr); 
  gkyl_mat_set(&A_bxby_yr,0,1,temp_magB2_yr); 
  gkyl_mat_set(&A_bybz_yr,0,1,temp_magB2_yr); 
 
  temp_magB2_xl = 0.7071067811865475*magB2_xl[2]; 
  gkyl_mat_set(&A_bxbx_xl,0,2,temp_magB2_xl); 
  gkyl_mat_set(&A_bxby_xl,0,2,temp_magB2_xl); 
  gkyl_mat_set(&A_bxbz_xl,0,2,temp_magB2_xl); 
 
  temp_magB2_xr = 0.7071067811865475*magB2_xr[2]; 
  gkyl_mat_set(&A_bxbx_xr,0,2,temp_magB2_xr); 
  gkyl_mat_set(&A_bxby_xr,0,2,temp_magB2_xr); 
  gkyl_mat_set(&A_bxbz_xr,0,2,temp_magB2_xr); 
 
  temp_magB2_yl = 0.7071067811865475*magB2_yl[2]; 
  gkyl_mat_set(&A_byby_yl,0,2,temp_magB2_yl); 
  gkyl_mat_set(&A_bxby_yl,0,2,temp_magB2_yl); 
  gkyl_mat_set(&A_bybz_yl,0,2,temp_magB2_yl); 
 
  temp_magB2_yr = 0.7071067811865475*magB2_yr[2]; 
  gkyl_mat_set(&A_byby_yr,0,2,temp_magB2_yr); 
  gkyl_mat_set(&A_bxby_yr,0,2,temp_magB2_yr); 
  gkyl_mat_set(&A_bybz_yr,0,2,temp_magB2_yr); 
 
  temp_magB2_xl = 0.7071067811865475*magB2_xl[1]; 
  gkyl_mat_set(&A_bxbx_xl,1,0,temp_magB2_xl); 
  gkyl_mat_set(&A_bxby_xl,1,0,temp_magB2_xl); 
  gkyl_mat_set(&A_bxbz_xl,1,0,temp_magB2_xl); 
 
  temp_magB2_xr = 0.7071067811865475*magB2_xr[1]; 
  gkyl_mat_set(&A_bxbx_xr,1,0,temp_magB2_xr); 
  gkyl_mat_set(&A_bxby_xr,1,0,temp_magB2_xr); 
  gkyl_mat_set(&A_bxbz_xr,1,0,temp_magB2_xr); 
 
  temp_magB2_yl = 0.7071067811865475*magB2_yl[1]; 
  gkyl_mat_set(&A_byby_yl,1,0,temp_magB2_yl); 
  gkyl_mat_set(&A_bxby_yl,1,0,temp_magB2_yl); 
  gkyl_mat_set(&A_bybz_yl,1,0,temp_magB2_yl); 
 
  temp_magB2_yr = 0.7071067811865475*magB2_yr[1]; 
  gkyl_mat_set(&A_byby_yr,1,0,temp_magB2_yr); 
  gkyl_mat_set(&A_bxby_yr,1,0,temp_magB2_yr); 
  gkyl_mat_set(&A_bybz_yr,1,0,temp_magB2_yr); 
 
  temp_magB2_xl = 0.6324555320336759*magB2_xl[2]+0.7071067811865475*magB2_xl[0]; 
  gkyl_mat_set(&A_bxbx_xl,1,1,temp_magB2_xl); 
  gkyl_mat_set(&A_bxby_xl,1,1,temp_magB2_xl); 
  gkyl_mat_set(&A_bxbz_xl,1,1,temp_magB2_xl); 
 
  temp_magB2_xr = 0.6324555320336759*magB2_xr[2]+0.7071067811865475*magB2_xr[0]; 
  gkyl_mat_set(&A_bxbx_xr,1,1,temp_magB2_xr); 
  gkyl_mat_set(&A_bxby_xr,1,1,temp_magB2_xr); 
  gkyl_mat_set(&A_bxbz_xr,1,1,temp_magB2_xr); 
 
  temp_magB2_yl = 0.6324555320336759*magB2_yl[2]+0.7071067811865475*magB2_yl[0]; 
  gkyl_mat_set(&A_byby_yl,1,1,temp_magB2_yl); 
  gkyl_mat_set(&A_bxby_yl,1,1,temp_magB2_yl); 
  gkyl_mat_set(&A_bybz_yl,1,1,temp_magB2_yl); 
 
  temp_magB2_yr = 0.6324555320336759*magB2_yr[2]+0.7071067811865475*magB2_yr[0]; 
  gkyl_mat_set(&A_byby_yr,1,1,temp_magB2_yr); 
  gkyl_mat_set(&A_bxby_yr,1,1,temp_magB2_yr); 
  gkyl_mat_set(&A_bybz_yr,1,1,temp_magB2_yr); 
 
  temp_magB2_xl = 0.6324555320336759*magB2_xl[1]; 
  gkyl_mat_set(&A_bxbx_xl,1,2,temp_magB2_xl); 
  gkyl_mat_set(&A_bxby_xl,1,2,temp_magB2_xl); 
  gkyl_mat_set(&A_bxbz_xl,1,2,temp_magB2_xl); 
 
  temp_magB2_xr = 0.6324555320336759*magB2_xr[1]; 
  gkyl_mat_set(&A_bxbx_xr,1,2,temp_magB2_xr); 
  gkyl_mat_set(&A_bxby_xr,1,2,temp_magB2_xr); 
  gkyl_mat_set(&A_bxbz_xr,1,2,temp_magB2_xr); 
 
  temp_magB2_yl = 0.6324555320336759*magB2_yl[1]; 
  gkyl_mat_set(&A_byby_yl,1,2,temp_magB2_yl); 
  gkyl_mat_set(&A_bxby_yl,1,2,temp_magB2_yl); 
  gkyl_mat_set(&A_bybz_yl,1,2,temp_magB2_yl); 
 
  temp_magB2_yr = 0.6324555320336759*magB2_yr[1]; 
  gkyl_mat_set(&A_byby_yr,1,2,temp_magB2_yr); 
  gkyl_mat_set(&A_bxby_yr,1,2,temp_magB2_yr); 
  gkyl_mat_set(&A_bybz_yr,1,2,temp_magB2_yr); 
 
  temp_magB2_xl = 0.7071067811865475*magB2_xl[2]; 
  gkyl_mat_set(&A_bxbx_xl,2,0,temp_magB2_xl); 
  gkyl_mat_set(&A_bxby_xl,2,0,temp_magB2_xl); 
  gkyl_mat_set(&A_bxbz_xl,2,0,temp_magB2_xl); 
 
  temp_magB2_xr = 0.7071067811865475*magB2_xr[2]; 
  gkyl_mat_set(&A_bxbx_xr,2,0,temp_magB2_xr); 
  gkyl_mat_set(&A_bxby_xr,2,0,temp_magB2_xr); 
  gkyl_mat_set(&A_bxbz_xr,2,0,temp_magB2_xr); 
 
  temp_magB2_yl = 0.7071067811865475*magB2_yl[2]; 
  gkyl_mat_set(&A_byby_yl,2,0,temp_magB2_yl); 
  gkyl_mat_set(&A_bxby_yl,2,0,temp_magB2_yl); 
  gkyl_mat_set(&A_bybz_yl,2,0,temp_magB2_yl); 
 
  temp_magB2_yr = 0.7071067811865475*magB2_yr[2]; 
  gkyl_mat_set(&A_byby_yr,2,0,temp_magB2_yr); 
  gkyl_mat_set(&A_bxby_yr,2,0,temp_magB2_yr); 
  gkyl_mat_set(&A_bybz_yr,2,0,temp_magB2_yr); 
 
  temp_magB2_xl = 0.6324555320336759*magB2_xl[1]; 
  gkyl_mat_set(&A_bxbx_xl,2,1,temp_magB2_xl); 
  gkyl_mat_set(&A_bxby_xl,2,1,temp_magB2_xl); 
  gkyl_mat_set(&A_bxbz_xl,2,1,temp_magB2_xl); 
 
  temp_magB2_xr = 0.6324555320336759*magB2_xr[1]; 
  gkyl_mat_set(&A_bxbx_xr,2,1,temp_magB2_xr); 
  gkyl_mat_set(&A_bxby_xr,2,1,temp_magB2_xr); 
  gkyl_mat_set(&A_bxbz_xr,2,1,temp_magB2_xr); 
 
  temp_magB2_yl = 0.6324555320336759*magB2_yl[1]; 
  gkyl_mat_set(&A_byby_yl,2,1,temp_magB2_yl); 
  gkyl_mat_set(&A_bxby_yl,2,1,temp_magB2_yl); 
  gkyl_mat_set(&A_bybz_yl,2,1,temp_magB2_yl); 
 
  temp_magB2_yr = 0.6324555320336759*magB2_yr[1]; 
  gkyl_mat_set(&A_byby_yr,2,1,temp_magB2_yr); 
  gkyl_mat_set(&A_bxby_yr,2,1,temp_magB2_yr); 
  gkyl_mat_set(&A_bybz_yr,2,1,temp_magB2_yr); 
 
  temp_magB2_xl = 0.4517539514526256*magB2_xl[2]+0.7071067811865475*magB2_xl[0]; 
  gkyl_mat_set(&A_bxbx_xl,2,2,temp_magB2_xl); 
  gkyl_mat_set(&A_bxby_xl,2,2,temp_magB2_xl); 
  gkyl_mat_set(&A_bxbz_xl,2,2,temp_magB2_xl); 
 
  temp_magB2_xr = 0.4517539514526256*magB2_xr[2]+0.7071067811865475*magB2_xr[0]; 
  gkyl_mat_set(&A_bxbx_xr,2,2,temp_magB2_xr); 
  gkyl_mat_set(&A_bxby_xr,2,2,temp_magB2_xr); 
  gkyl_mat_set(&A_bxbz_xr,2,2,temp_magB2_xr); 
 
  temp_magB2_yl = 0.4517539514526256*magB2_yl[2]+0.7071067811865475*magB2_yl[0]; 
  gkyl_mat_set(&A_byby_yl,2,2,temp_magB2_yl); 
  gkyl_mat_set(&A_bxby_yl,2,2,temp_magB2_yl); 
  gkyl_mat_set(&A_bybz_yl,2,2,temp_magB2_yl); 
 
  temp_magB2_yr = 0.4517539514526256*magB2_yr[2]+0.7071067811865475*magB2_yr[0]; 
  gkyl_mat_set(&A_byby_yr,2,2,temp_magB2_yr); 
  gkyl_mat_set(&A_bxby_yr,2,2,temp_magB2_yr); 
  gkyl_mat_set(&A_bybz_yr,2,2,temp_magB2_yr); 
 
} 
