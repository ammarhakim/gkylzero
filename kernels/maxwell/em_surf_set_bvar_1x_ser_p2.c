#include <gkyl_mat.h> 
#include <gkyl_maxwell_kernels.h> 
GKYL_CU_DH void em_surf_set_bvar_1x_ser_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *BB_surf, int* cell_avg_magB2_surf) 
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
  // Clear rhs for each component of surface magnetic field unit tensor and vector being solved for 
  gkyl_mat_clear(&rhs_bxbx_xl, 0.0); 
  gkyl_mat_clear(&rhs_bxbx_xr, 0.0); 
  const double *Bx_sq_xl = &BB_surf[0]; 
  const double *Bx_sq_xr = &BB_surf[1]; 
  const double *By_sq_xl = &BB_surf[2]; 
  const double *By_sq_xr = &BB_surf[3]; 
  const double *Bz_sq_xl = &BB_surf[4]; 
  const double *Bz_sq_xr = &BB_surf[5]; 
  int *cell_avg_magB2_xl = &cell_avg_magB2_surf[0]; 
  int *cell_avg_magB2_xr = &cell_avg_magB2_surf[1]; 
 
  double magB2_xl = Bx_sq_xl[0] + By_sq_xl[0] + Bz_sq_xl[0]; 
  double magB2_xr = Bx_sq_xr[0] + By_sq_xr[0] + Bz_sq_xr[0]; 
 
  cell_avg_magB2_xl[0] = 0; 
  cell_avg_magB2_xr[0] = 0; 
 
  double bxbx_xl = Bx_sq_xl[0]/magB2_xl; 
  double bxbx_xr = Bx_sq_xr[0]/magB2_xr; 
  gkyl_mat_set(&rhs_bxbx_xl,0,0,bxbx_xl); 
  gkyl_mat_set(&rhs_bxbx_xr,0,0,bxbx_xr); 
} 
