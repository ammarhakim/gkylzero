#include <gkyl_mat.h> 
#include <gkyl_maxwell_kernels.h> 
GKYL_CU_DH void em_surf_copy_bvar_1x_ser_p1(int count, struct gkyl_nmat *x, const double *em, int* cell_avg_magB2_surf, double* GKYL_RESTRICT bvar_surf) 
{ 
  // count:               Integer to indicate which matrix being fetched. 
  // x:                   Input solution vector. 
  // em:                  Input electromagnetic fields. 
  // cell_avg_magB2_surf: Output flag for cell average if 1/|B|^2 at a surface only used cell averages. 
  // bvar_surf:           Output magnetic field unit tensor and unit vector at surfaces. 
  //                      [bx_xl, bx_xr, by_yl, by_yr, bz_zl, bz_zr] 
 
  double bxbx_xl[1] = {0.0}; 
  double bxbx_xr[1] = {0.0}; 
  struct gkyl_mat x_bxbx_xl = gkyl_nmat_get(x, count); 
  struct gkyl_mat x_bxbx_xr = gkyl_nmat_get(x, count+1); 
 
  double *bx_xl = &bvar_surf[0]; 
  double *bx_xr = &bvar_surf[1]; 
  int *cell_avg_magB2_xl = &cell_avg_magB2_surf[0]; 
  int *cell_avg_magB2_xr = &cell_avg_magB2_surf[1]; 
 
  bxbx_xl[0] = gkyl_mat_get(&x_bxbx_xl,0,0); 
  bxbx_xr[0] = gkyl_mat_get(&x_bxbx_xr,0,0); 
 
  const double *B_x = &em[6]; 
  const double *B_y = &em[8]; 
  const double *B_z = &em[10]; 
 
  double B_x_xl = 0.7071067811865475*B_x[0]-1.224744871391589*B_x[1]; 
  double B_x_xr = 1.224744871391589*B_x[1]+0.7071067811865475*B_x[0]; 
  if (B_x_xl < 0.0) 
    bx_xl[0] = -sqrt(bxbx_xl[0]); 
  else 
    bx_xl[0] = sqrt(bxbx_xl[0]); 
 
  if (B_x_xr < 0.0) 
    bx_xr[0] = -sqrt(bxbx_xr[0]); 
  else 
    bx_xr[0] = sqrt(bxbx_xr[0]); 
} 
 
