#include <gkyl_mat.h> 
#include <gkyl_maxwell_kernels.h> 
#include <gkyl_basis_ser_1x_p2_sqrt_with_sign.h> 
GKYL_CU_DH void em_surf_copy_bvar_2x_tensor_p2(int count, struct gkyl_nmat *x, const double *em, int* cell_avg_magB2_surf, double* GKYL_RESTRICT bvar_surf) 
{ 
  // count:               Integer to indicate which matrix being fetched. 
  // x:                   Input solution vector. 
  // em:                  Input electromagnetic fields. 
  // cell_avg_magB2_surf: Output flag for cell average if 1/|B|^2 at a surface only used cell averages. 
  // bvar_surf:           Output magnetic field unit tensor and unit vector at surfaces. 
  //                      [bx_xl, bx_xr, by_yl, by_yr, bz_zl, bz_zr] 
 
  double bxbx_xl[3] = {0.0}; 
  double bxbx_xr[3] = {0.0}; 
  struct gkyl_mat x_bxbx_xl = gkyl_nmat_get(x, count); 
  struct gkyl_mat x_bxbx_xr = gkyl_nmat_get(x, count+1); 
 
  double *bx_xl = &bvar_surf[0]; 
  double *bx_xr = &bvar_surf[3]; 
  int *cell_avg_magB2_xl = &cell_avg_magB2_surf[0]; 
  int *cell_avg_magB2_xr = &cell_avg_magB2_surf[1]; 
 
  double byby_yl[3] = {0.0}; 
  double byby_yr[3] = {0.0}; 
  struct gkyl_mat x_byby_yl = gkyl_nmat_get(x, count+2); 
  struct gkyl_mat x_byby_yr = gkyl_nmat_get(x, count+3); 
 
  double *by_yl = &bvar_surf[6]; 
  double *by_yr = &bvar_surf[9]; 
  int *cell_avg_magB2_yl = &cell_avg_magB2_surf[2]; 
  int *cell_avg_magB2_yr = &cell_avg_magB2_surf[3]; 
 
  bxbx_xl[0] = gkyl_mat_get(&x_bxbx_xl,0,0); 
  bxbx_xr[0] = gkyl_mat_get(&x_bxbx_xr,0,0); 
 
  byby_yl[0] = gkyl_mat_get(&x_byby_yl,0,0); 
  byby_yr[0] = gkyl_mat_get(&x_byby_yr,0,0); 
 
  bxbx_xl[1] = gkyl_mat_get(&x_bxbx_xl,1,0); 
  bxbx_xr[1] = gkyl_mat_get(&x_bxbx_xr,1,0); 
 
  byby_yl[1] = gkyl_mat_get(&x_byby_yl,1,0); 
  byby_yr[1] = gkyl_mat_get(&x_byby_yr,1,0); 
 
  bxbx_xl[2] = gkyl_mat_get(&x_bxbx_xl,2,0); 
  bxbx_xr[2] = gkyl_mat_get(&x_bxbx_xr,2,0); 
 
  byby_yl[2] = gkyl_mat_get(&x_byby_yl,2,0); 
  byby_yr[2] = gkyl_mat_get(&x_byby_yr,2,0); 
 
  const double *B_x = &em[27]; 
  const double *B_y = &em[36]; 
  const double *B_z = &em[45]; 
 
  int cell_avg_xl = 0;
  int cell_avg_xr = 0;
  if (0.6324555320336759*bxbx_xl[2]-0.9486832980505137*bxbx_xl[1]+0.7071067811865475*bxbx_xl[0] < 0.0) cell_avg_xl = 1; 
  if (0.6324555320336759*bxbx_xr[2]-0.9486832980505137*bxbx_xr[1]+0.7071067811865475*bxbx_xr[0] < 0.0) cell_avg_xr = 1; 
  if (0.7071067811865475*bxbx_xl[0]-0.7905694150420947*bxbx_xl[2] < 0.0) cell_avg_xl = 1; 
  if (0.7071067811865475*bxbx_xr[0]-0.7905694150420947*bxbx_xr[2] < 0.0) cell_avg_xr = 1; 
  if (0.6324555320336759*bxbx_xl[2]+0.9486832980505137*bxbx_xl[1]+0.7071067811865475*bxbx_xl[0] < 0.0) cell_avg_xl = 1; 
  if (0.6324555320336759*bxbx_xr[2]+0.9486832980505137*bxbx_xr[1]+0.7071067811865475*bxbx_xr[0] < 0.0) cell_avg_xr = 1; 
 
  if (cell_avg_xl || cell_avg_magB2_xl[0]) { 
    bxbx_xl[1] = 0.0; 
    bxbx_xl[2] = 0.0; 
    // If bxbx < 0.0 at the lower x surface quadrature points, 
    // set cell_avg_magB2_xl to be true in case it was not true before. 
    cell_avg_magB2_xl[0] = 1; 
  } 
 
  if (cell_avg_xr || cell_avg_magB2_xr[0]) { 
    bxbx_xr[1] = 0.0; 
    bxbx_xr[2] = 0.0; 
    // If bxbx < 0.0 at the upper x surface quadrature points, 
    // set cell_avg_magB2_xr to be true in case it was not true before. 
    cell_avg_magB2_xr[0] = 1; 
  } 
 
  int cell_avg_yl = 0;
  int cell_avg_yr = 0;
  if (0.6324555320336759*byby_yl[2]-0.9486832980505137*byby_yl[1]+0.7071067811865475*byby_yl[0] < 0.0) cell_avg_yl = 1; 
  if (0.6324555320336759*byby_yr[2]-0.9486832980505137*byby_yr[1]+0.7071067811865475*byby_yr[0] < 0.0) cell_avg_yr = 1; 
  if (0.7071067811865475*byby_yl[0]-0.7905694150420947*byby_yl[2] < 0.0) cell_avg_yl = 1; 
  if (0.7071067811865475*byby_yr[0]-0.7905694150420947*byby_yr[2] < 0.0) cell_avg_yr = 1; 
  if (0.6324555320336759*byby_yl[2]+0.9486832980505137*byby_yl[1]+0.7071067811865475*byby_yl[0] < 0.0) cell_avg_yl = 1; 
  if (0.6324555320336759*byby_yr[2]+0.9486832980505137*byby_yr[1]+0.7071067811865475*byby_yr[0] < 0.0) cell_avg_yr = 1; 
 
  if (cell_avg_yl || cell_avg_magB2_yl[0]) { 
    byby_yl[1] = 0.0; 
    byby_yl[2] = 0.0; 
    // If byby < 0.0 at the lower y surface quadrature points, 
    // set cell_avg_magB2_yl to be true in case it was not true before. 
    cell_avg_magB2_yl[0] = 1; 
  } 
 
  if (cell_avg_yr || cell_avg_magB2_yr[0]) { 
    byby_yr[1] = 0.0; 
    byby_yr[2] = 0.0; 
    // If byby < 0.0 at the upper y surface quadrature points, 
    // set cell_avg_magB2_yr to be true in case it was not true before. 
    cell_avg_magB2_yr[0] = 1; 
  } 
 
  // Calculate b_i = B_i/|B| by taking square root of B_i^2/|B|^2 at quadrature points. 
  // Uses the sign of B_i at quadrature points *on the surface* to get the correct sign of b_i *on the surface*. 
  // Note: positivity check already happened, so only uses cell average if needed to avoid imaginary values of b_hat. 
  double B_x_xl[3] = {0.0}; 
  double B_x_xr[3] = {0.0}; 
 
  B_x_xl[0] = 1.58113883008419*B_x[4]-1.224744871391589*B_x[1]+0.7071067811865475*B_x[0]; 
  B_x_xl[1] = 1.58113883008419*B_x[6]-1.224744871391589*B_x[3]+0.7071067811865475*B_x[2]; 
  B_x_xl[2] = 1.58113883008419*B_x[8]-1.224744871391589*B_x[7]+0.7071067811865475*B_x[5]; 
  B_x_xr[0] = 1.58113883008419*B_x[4]+1.224744871391589*B_x[1]+0.7071067811865475*B_x[0]; 
  B_x_xr[1] = 1.58113883008419*B_x[6]+1.224744871391589*B_x[3]+0.7071067811865475*B_x[2]; 
  B_x_xr[2] = 1.58113883008419*B_x[8]+1.224744871391589*B_x[7]+0.7071067811865475*B_x[5]; 
  ser_1x_p2_sqrt_with_sign(B_x_xl, bxbx_xl, bx_xl); 
  ser_1x_p2_sqrt_with_sign(B_x_xr, bxbx_xr, bx_xr); 
 
  double B_y_yl[3] = {0.0}; 
  double B_y_yr[3] = {0.0}; 
 
  B_y_yl[0] = 1.58113883008419*B_y[5]-1.224744871391589*B_y[2]+0.7071067811865475*B_y[0]; 
  B_y_yl[1] = 1.58113883008419*B_y[7]-1.224744871391589*B_y[3]+0.7071067811865475*B_y[1]; 
  B_y_yl[2] = 1.58113883008419*B_y[8]-1.224744871391589*B_y[6]+0.7071067811865475*B_y[4]; 
  B_y_yr[0] = 1.58113883008419*B_y[5]+1.224744871391589*B_y[2]+0.7071067811865475*B_y[0]; 
  B_y_yr[1] = 1.58113883008419*B_y[7]+1.224744871391589*B_y[3]+0.7071067811865475*B_y[1]; 
  B_y_yr[2] = 1.58113883008419*B_y[8]+1.224744871391589*B_y[6]+0.7071067811865475*B_y[4]; 
  ser_1x_p2_sqrt_with_sign(B_y_yl, byby_yl, by_yl); 
  ser_1x_p2_sqrt_with_sign(B_y_yr, byby_yr, by_yr); 
 
} 
 
