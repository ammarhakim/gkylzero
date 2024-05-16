#include <gkyl_mat.h> 
#include <gkyl_maxwell_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
GKYL_CU_DH void em_surf_set_bvar_3x_tensor_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *BB_surf, int* cell_avg_magB2_surf) 
{ 
  // count:   integer to indicate which matrix being fetched. 
  // A:       preallocated LHS matrix. 
  // rhs:     preallocated RHS vector. 
  // BB_surf: Surface B_i B_j [BxBx_xl, BxBx_xr, ByBy_xl, ByBy_xr, BzBz_xl, BzBz_xr,  
  //                           BxBx_yl, BxBx_yr, ByBy_yl, ByBy_yr, BzBz_yl, BzBz_yr,  
  //                           BxBx_zl, BxBx_zr, ByBy_zl, ByBy_zr, BzBz_zl, BzBz_zr]. 
  // cell_avg_magB2_surf:      Output flag for cell average if 1/|B|^2 at a surface only used cell averages. 

  struct gkyl_mat A_bxbx_xl = gkyl_nmat_get(A, count); 
  struct gkyl_mat A_bxbx_xr = gkyl_nmat_get(A, count+1); 
  struct gkyl_mat rhs_bxbx_xl = gkyl_nmat_get(rhs, count); 
  struct gkyl_mat rhs_bxbx_xr = gkyl_nmat_get(rhs, count+1); 
  gkyl_mat_clear(&A_bxbx_xl, 0.0); gkyl_mat_clear(&rhs_bxbx_xl, 0.0); 
  gkyl_mat_clear(&A_bxbx_xr, 0.0); gkyl_mat_clear(&rhs_bxbx_xr, 0.0); 
  const double *Bx_sq_xl = &BB_surf[0]; 
  const double *Bx_sq_xr = &BB_surf[9]; 
  const double *By_sq_xl = &BB_surf[18]; 
  const double *By_sq_xr = &BB_surf[27]; 
  const double *Bz_sq_xl = &BB_surf[36]; 
  const double *Bz_sq_xr = &BB_surf[45]; 
  int *cell_avg_magB2_xl = &cell_avg_magB2_surf[0]; 
  int *cell_avg_magB2_xr = &cell_avg_magB2_surf[1]; 
 
  struct gkyl_mat A_byby_yl = gkyl_nmat_get(A, count+2); 
  struct gkyl_mat A_byby_yr = gkyl_nmat_get(A, count+3); 
  struct gkyl_mat rhs_byby_yl = gkyl_nmat_get(rhs, count+2); 
  struct gkyl_mat rhs_byby_yr = gkyl_nmat_get(rhs, count+3); 
  gkyl_mat_clear(&A_byby_yl, 0.0); gkyl_mat_clear(&rhs_byby_yl, 0.0); 
  gkyl_mat_clear(&A_byby_yr, 0.0); gkyl_mat_clear(&rhs_byby_yr, 0.0); 
  const double *Bx_sq_yl = &BB_surf[54]; 
  const double *Bx_sq_yr = &BB_surf[63]; 
  const double *By_sq_yl = &BB_surf[72]; 
  const double *By_sq_yr = &BB_surf[81]; 
  const double *Bz_sq_yl = &BB_surf[90]; 
  const double *Bz_sq_yr = &BB_surf[99]; 
  int *cell_avg_magB2_yl = &cell_avg_magB2_surf[2]; 
  int *cell_avg_magB2_yr = &cell_avg_magB2_surf[3]; 
 
  struct gkyl_mat A_bzbz_zl = gkyl_nmat_get(A, count+4); 
  struct gkyl_mat A_bzbz_zr = gkyl_nmat_get(A, count+5); 
  struct gkyl_mat rhs_bzbz_zl = gkyl_nmat_get(rhs, count+4); 
  struct gkyl_mat rhs_bzbz_zr = gkyl_nmat_get(rhs, count+5); 
  gkyl_mat_clear(&A_bzbz_zl, 0.0); gkyl_mat_clear(&rhs_bzbz_zl, 0.0); 
  gkyl_mat_clear(&A_bzbz_zr, 0.0); gkyl_mat_clear(&rhs_bzbz_zr, 0.0); 
  const double *Bx_sq_zl = &BB_surf[108]; 
  const double *Bx_sq_zr = &BB_surf[117]; 
  const double *By_sq_zl = &BB_surf[126]; 
  const double *By_sq_zr = &BB_surf[135]; 
  const double *Bz_sq_zl = &BB_surf[144]; 
  const double *Bz_sq_zr = &BB_surf[153]; 
  int *cell_avg_magB2_zl = &cell_avg_magB2_surf[4]; 
  int *cell_avg_magB2_zr = &cell_avg_magB2_surf[5]; 
 
  double magB2_xl[9] = {0.0}; 
  double magB2_xr[9] = {0.0}; 
  double magB2_yl[9] = {0.0}; 
  double magB2_yr[9] = {0.0}; 
  double magB2_zl[9] = {0.0}; 
  double magB2_zr[9] = {0.0}; 
  magB2_xl[0] = Bx_sq_xl[0] + By_sq_xl[0] + Bz_sq_xl[0]; 
  magB2_xr[0] = Bx_sq_xr[0] + By_sq_xr[0] + Bz_sq_xr[0]; 
  magB2_yl[0] = Bx_sq_yl[0] + By_sq_yl[0] + Bz_sq_yl[0]; 
  magB2_yr[0] = Bx_sq_yr[0] + By_sq_yr[0] + Bz_sq_yr[0]; 
  magB2_zl[0] = Bx_sq_zl[0] + By_sq_zl[0] + Bz_sq_zl[0]; 
  magB2_zr[0] = Bx_sq_zr[0] + By_sq_zr[0] + Bz_sq_zr[0]; 
  magB2_xl[1] = Bx_sq_xl[1] + By_sq_xl[1] + Bz_sq_xl[1]; 
  magB2_xr[1] = Bx_sq_xr[1] + By_sq_xr[1] + Bz_sq_xr[1]; 
  magB2_yl[1] = Bx_sq_yl[1] + By_sq_yl[1] + Bz_sq_yl[1]; 
  magB2_yr[1] = Bx_sq_yr[1] + By_sq_yr[1] + Bz_sq_yr[1]; 
  magB2_zl[1] = Bx_sq_zl[1] + By_sq_zl[1] + Bz_sq_zl[1]; 
  magB2_zr[1] = Bx_sq_zr[1] + By_sq_zr[1] + Bz_sq_zr[1]; 
  magB2_xl[2] = Bx_sq_xl[2] + By_sq_xl[2] + Bz_sq_xl[2]; 
  magB2_xr[2] = Bx_sq_xr[2] + By_sq_xr[2] + Bz_sq_xr[2]; 
  magB2_yl[2] = Bx_sq_yl[2] + By_sq_yl[2] + Bz_sq_yl[2]; 
  magB2_yr[2] = Bx_sq_yr[2] + By_sq_yr[2] + Bz_sq_yr[2]; 
  magB2_zl[2] = Bx_sq_zl[2] + By_sq_zl[2] + Bz_sq_zl[2]; 
  magB2_zr[2] = Bx_sq_zr[2] + By_sq_zr[2] + Bz_sq_zr[2]; 
  magB2_xl[3] = Bx_sq_xl[3] + By_sq_xl[3] + Bz_sq_xl[3]; 
  magB2_xr[3] = Bx_sq_xr[3] + By_sq_xr[3] + Bz_sq_xr[3]; 
  magB2_yl[3] = Bx_sq_yl[3] + By_sq_yl[3] + Bz_sq_yl[3]; 
  magB2_yr[3] = Bx_sq_yr[3] + By_sq_yr[3] + Bz_sq_yr[3]; 
  magB2_zl[3] = Bx_sq_zl[3] + By_sq_zl[3] + Bz_sq_zl[3]; 
  magB2_zr[3] = Bx_sq_zr[3] + By_sq_zr[3] + Bz_sq_zr[3]; 
  magB2_xl[4] = Bx_sq_xl[4] + By_sq_xl[4] + Bz_sq_xl[4]; 
  magB2_xr[4] = Bx_sq_xr[4] + By_sq_xr[4] + Bz_sq_xr[4]; 
  magB2_yl[4] = Bx_sq_yl[4] + By_sq_yl[4] + Bz_sq_yl[4]; 
  magB2_yr[4] = Bx_sq_yr[4] + By_sq_yr[4] + Bz_sq_yr[4]; 
  magB2_zl[4] = Bx_sq_zl[4] + By_sq_zl[4] + Bz_sq_zl[4]; 
  magB2_zr[4] = Bx_sq_zr[4] + By_sq_zr[4] + Bz_sq_zr[4]; 
  magB2_xl[5] = Bx_sq_xl[5] + By_sq_xl[5] + Bz_sq_xl[5]; 
  magB2_xr[5] = Bx_sq_xr[5] + By_sq_xr[5] + Bz_sq_xr[5]; 
  magB2_yl[5] = Bx_sq_yl[5] + By_sq_yl[5] + Bz_sq_yl[5]; 
  magB2_yr[5] = Bx_sq_yr[5] + By_sq_yr[5] + Bz_sq_yr[5]; 
  magB2_zl[5] = Bx_sq_zl[5] + By_sq_zl[5] + Bz_sq_zl[5]; 
  magB2_zr[5] = Bx_sq_zr[5] + By_sq_zr[5] + Bz_sq_zr[5]; 
  magB2_xl[6] = Bx_sq_xl[6] + By_sq_xl[6] + Bz_sq_xl[6]; 
  magB2_xr[6] = Bx_sq_xr[6] + By_sq_xr[6] + Bz_sq_xr[6]; 
  magB2_yl[6] = Bx_sq_yl[6] + By_sq_yl[6] + Bz_sq_yl[6]; 
  magB2_yr[6] = Bx_sq_yr[6] + By_sq_yr[6] + Bz_sq_yr[6]; 
  magB2_zl[6] = Bx_sq_zl[6] + By_sq_zl[6] + Bz_sq_zl[6]; 
  magB2_zr[6] = Bx_sq_zr[6] + By_sq_zr[6] + Bz_sq_zr[6]; 
  magB2_xl[7] = Bx_sq_xl[7] + By_sq_xl[7] + Bz_sq_xl[7]; 
  magB2_xr[7] = Bx_sq_xr[7] + By_sq_xr[7] + Bz_sq_xr[7]; 
  magB2_yl[7] = Bx_sq_yl[7] + By_sq_yl[7] + Bz_sq_yl[7]; 
  magB2_yr[7] = Bx_sq_yr[7] + By_sq_yr[7] + Bz_sq_yr[7]; 
  magB2_zl[7] = Bx_sq_zl[7] + By_sq_zl[7] + Bz_sq_zl[7]; 
  magB2_zr[7] = Bx_sq_zr[7] + By_sq_zr[7] + Bz_sq_zr[7]; 
  magB2_xl[8] = Bx_sq_xl[8] + By_sq_xl[8] + Bz_sq_xl[8]; 
  magB2_xr[8] = Bx_sq_xr[8] + By_sq_xr[8] + Bz_sq_xr[8]; 
  magB2_yl[8] = Bx_sq_yl[8] + By_sq_yl[8] + Bz_sq_yl[8]; 
  magB2_yr[8] = Bx_sq_yr[8] + By_sq_yr[8] + Bz_sq_yr[8]; 
  magB2_zl[8] = Bx_sq_zl[8] + By_sq_zl[8] + Bz_sq_zl[8]; 
  magB2_zr[8] = Bx_sq_zr[8] + By_sq_zr[8] + Bz_sq_zr[8]; 
  // If |B|^2 < 0 at control points along a surface, only use cell average to get 1/|B|^2. 
  // Each surface is checked independently. 
  int cell_avg_xl = 0;
  int cell_avg_xr = 0;
  int cell_avg_yl = 0;
  int cell_avg_yr = 0;
  int cell_avg_zl = 0;
  int cell_avg_zr = 0;
 
  if (2.5*magB2_xl[8]-1.936491673103709*magB2_xl[7]-1.936491673103709*magB2_xl[6]+1.118033988749895*magB2_xl[5]+1.118033988749895*magB2_xl[4]+1.5*magB2_xl[3]-0.8660254037844386*magB2_xl[2]-0.8660254037844386*magB2_xl[1]+0.5*magB2_xl[0] < 0.0) cell_avg_xl = 1; 
  if (2.5*magB2_xr[8]-1.936491673103709*magB2_xr[7]-1.936491673103709*magB2_xr[6]+1.118033988749895*magB2_xr[5]+1.118033988749895*magB2_xr[4]+1.5*magB2_xr[3]-0.8660254037844386*magB2_xr[2]-0.8660254037844386*magB2_xr[1]+0.5*magB2_xr[0] < 0.0) cell_avg_xr = 1; 
  if (2.5*magB2_yl[8]-1.936491673103709*magB2_yl[7]-1.936491673103709*magB2_yl[6]+1.118033988749895*magB2_yl[5]+1.118033988749895*magB2_yl[4]+1.5*magB2_yl[3]-0.8660254037844386*magB2_yl[2]-0.8660254037844386*magB2_yl[1]+0.5*magB2_yl[0] < 0.0) cell_avg_yl = 1; 
  if (2.5*magB2_yr[8]-1.936491673103709*magB2_yr[7]-1.936491673103709*magB2_yr[6]+1.118033988749895*magB2_yr[5]+1.118033988749895*magB2_yr[4]+1.5*magB2_yr[3]-0.8660254037844386*magB2_yr[2]-0.8660254037844386*magB2_yr[1]+0.5*magB2_yr[0] < 0.0) cell_avg_yr = 1; 
  if (2.5*magB2_zl[8]-1.936491673103709*magB2_zl[7]-1.936491673103709*magB2_zl[6]+1.118033988749895*magB2_zl[5]+1.118033988749895*magB2_zl[4]+1.5*magB2_zl[3]-0.8660254037844386*magB2_zl[2]-0.8660254037844386*magB2_zl[1]+0.5*magB2_zl[0] < 0.0) cell_avg_zl = 1; 
  if (2.5*magB2_zr[8]-1.936491673103709*magB2_zr[7]-1.936491673103709*magB2_zr[6]+1.118033988749895*magB2_zr[5]+1.118033988749895*magB2_zr[4]+1.5*magB2_zr[3]-0.8660254037844386*magB2_zr[2]-0.8660254037844386*magB2_zr[1]+0.5*magB2_zr[0] < 0.0) cell_avg_zr = 1; 
  if ((-1.25*magB2_xl[8])+0.9682458365518543*magB2_xl[6]+1.118033988749895*magB2_xl[5]-0.5590169943749475*magB2_xl[4]-0.8660254037844386*magB2_xl[2]+0.5*magB2_xl[0] < 0.0) cell_avg_xl = 1; 
  if ((-1.25*magB2_xr[8])+0.9682458365518543*magB2_xr[6]+1.118033988749895*magB2_xr[5]-0.5590169943749475*magB2_xr[4]-0.8660254037844386*magB2_xr[2]+0.5*magB2_xr[0] < 0.0) cell_avg_xr = 1; 
  if ((-1.25*magB2_yl[8])+0.9682458365518543*magB2_yl[6]+1.118033988749895*magB2_yl[5]-0.5590169943749475*magB2_yl[4]-0.8660254037844386*magB2_yl[2]+0.5*magB2_yl[0] < 0.0) cell_avg_yl = 1; 
  if ((-1.25*magB2_yr[8])+0.9682458365518543*magB2_yr[6]+1.118033988749895*magB2_yr[5]-0.5590169943749475*magB2_yr[4]-0.8660254037844386*magB2_yr[2]+0.5*magB2_yr[0] < 0.0) cell_avg_yr = 1; 
  if ((-1.25*magB2_zl[8])+0.9682458365518543*magB2_zl[6]+1.118033988749895*magB2_zl[5]-0.5590169943749475*magB2_zl[4]-0.8660254037844386*magB2_zl[2]+0.5*magB2_zl[0] < 0.0) cell_avg_zl = 1; 
  if ((-1.25*magB2_zr[8])+0.9682458365518543*magB2_zr[6]+1.118033988749895*magB2_zr[5]-0.5590169943749475*magB2_zr[4]-0.8660254037844386*magB2_zr[2]+0.5*magB2_zr[0] < 0.0) cell_avg_zr = 1; 
  if (2.5*magB2_xl[8]+1.936491673103709*magB2_xl[7]-1.936491673103709*magB2_xl[6]+1.118033988749895*magB2_xl[5]+1.118033988749895*magB2_xl[4]-1.5*magB2_xl[3]-0.8660254037844386*magB2_xl[2]+0.8660254037844386*magB2_xl[1]+0.5*magB2_xl[0] < 0.0) cell_avg_xl = 1; 
  if (2.5*magB2_xr[8]+1.936491673103709*magB2_xr[7]-1.936491673103709*magB2_xr[6]+1.118033988749895*magB2_xr[5]+1.118033988749895*magB2_xr[4]-1.5*magB2_xr[3]-0.8660254037844386*magB2_xr[2]+0.8660254037844386*magB2_xr[1]+0.5*magB2_xr[0] < 0.0) cell_avg_xr = 1; 
  if (2.5*magB2_yl[8]+1.936491673103709*magB2_yl[7]-1.936491673103709*magB2_yl[6]+1.118033988749895*magB2_yl[5]+1.118033988749895*magB2_yl[4]-1.5*magB2_yl[3]-0.8660254037844386*magB2_yl[2]+0.8660254037844386*magB2_yl[1]+0.5*magB2_yl[0] < 0.0) cell_avg_yl = 1; 
  if (2.5*magB2_yr[8]+1.936491673103709*magB2_yr[7]-1.936491673103709*magB2_yr[6]+1.118033988749895*magB2_yr[5]+1.118033988749895*magB2_yr[4]-1.5*magB2_yr[3]-0.8660254037844386*magB2_yr[2]+0.8660254037844386*magB2_yr[1]+0.5*magB2_yr[0] < 0.0) cell_avg_yr = 1; 
  if (2.5*magB2_zl[8]+1.936491673103709*magB2_zl[7]-1.936491673103709*magB2_zl[6]+1.118033988749895*magB2_zl[5]+1.118033988749895*magB2_zl[4]-1.5*magB2_zl[3]-0.8660254037844386*magB2_zl[2]+0.8660254037844386*magB2_zl[1]+0.5*magB2_zl[0] < 0.0) cell_avg_zl = 1; 
  if (2.5*magB2_zr[8]+1.936491673103709*magB2_zr[7]-1.936491673103709*magB2_zr[6]+1.118033988749895*magB2_zr[5]+1.118033988749895*magB2_zr[4]-1.5*magB2_zr[3]-0.8660254037844386*magB2_zr[2]+0.8660254037844386*magB2_zr[1]+0.5*magB2_zr[0] < 0.0) cell_avg_zr = 1; 
  if ((-1.25*magB2_xl[8])+0.9682458365518543*magB2_xl[7]-0.5590169943749475*magB2_xl[5]+1.118033988749895*magB2_xl[4]-0.8660254037844386*magB2_xl[1]+0.5*magB2_xl[0] < 0.0) cell_avg_xl = 1; 
  if ((-1.25*magB2_xr[8])+0.9682458365518543*magB2_xr[7]-0.5590169943749475*magB2_xr[5]+1.118033988749895*magB2_xr[4]-0.8660254037844386*magB2_xr[1]+0.5*magB2_xr[0] < 0.0) cell_avg_xr = 1; 
  if ((-1.25*magB2_yl[8])+0.9682458365518543*magB2_yl[7]-0.5590169943749475*magB2_yl[5]+1.118033988749895*magB2_yl[4]-0.8660254037844386*magB2_yl[1]+0.5*magB2_yl[0] < 0.0) cell_avg_yl = 1; 
  if ((-1.25*magB2_yr[8])+0.9682458365518543*magB2_yr[7]-0.5590169943749475*magB2_yr[5]+1.118033988749895*magB2_yr[4]-0.8660254037844386*magB2_yr[1]+0.5*magB2_yr[0] < 0.0) cell_avg_yr = 1; 
  if ((-1.25*magB2_zl[8])+0.9682458365518543*magB2_zl[7]-0.5590169943749475*magB2_zl[5]+1.118033988749895*magB2_zl[4]-0.8660254037844386*magB2_zl[1]+0.5*magB2_zl[0] < 0.0) cell_avg_zl = 1; 
  if ((-1.25*magB2_zr[8])+0.9682458365518543*magB2_zr[7]-0.5590169943749475*magB2_zr[5]+1.118033988749895*magB2_zr[4]-0.8660254037844386*magB2_zr[1]+0.5*magB2_zr[0] < 0.0) cell_avg_zr = 1; 
  if (0.625*magB2_xl[8]-0.5590169943749475*magB2_xl[5]-0.5590169943749475*magB2_xl[4]+0.5*magB2_xl[0] < 0.0) cell_avg_xl = 1; 
  if (0.625*magB2_xr[8]-0.5590169943749475*magB2_xr[5]-0.5590169943749475*magB2_xr[4]+0.5*magB2_xr[0] < 0.0) cell_avg_xr = 1; 
  if (0.625*magB2_yl[8]-0.5590169943749475*magB2_yl[5]-0.5590169943749475*magB2_yl[4]+0.5*magB2_yl[0] < 0.0) cell_avg_yl = 1; 
  if (0.625*magB2_yr[8]-0.5590169943749475*magB2_yr[5]-0.5590169943749475*magB2_yr[4]+0.5*magB2_yr[0] < 0.0) cell_avg_yr = 1; 
  if (0.625*magB2_zl[8]-0.5590169943749475*magB2_zl[5]-0.5590169943749475*magB2_zl[4]+0.5*magB2_zl[0] < 0.0) cell_avg_zl = 1; 
  if (0.625*magB2_zr[8]-0.5590169943749475*magB2_zr[5]-0.5590169943749475*magB2_zr[4]+0.5*magB2_zr[0] < 0.0) cell_avg_zr = 1; 
  if ((-1.25*magB2_xl[8])-0.9682458365518543*magB2_xl[7]-0.5590169943749475*magB2_xl[5]+1.118033988749895*magB2_xl[4]+0.8660254037844386*magB2_xl[1]+0.5*magB2_xl[0] < 0.0) cell_avg_xl = 1; 
  if ((-1.25*magB2_xr[8])-0.9682458365518543*magB2_xr[7]-0.5590169943749475*magB2_xr[5]+1.118033988749895*magB2_xr[4]+0.8660254037844386*magB2_xr[1]+0.5*magB2_xr[0] < 0.0) cell_avg_xr = 1; 
  if ((-1.25*magB2_yl[8])-0.9682458365518543*magB2_yl[7]-0.5590169943749475*magB2_yl[5]+1.118033988749895*magB2_yl[4]+0.8660254037844386*magB2_yl[1]+0.5*magB2_yl[0] < 0.0) cell_avg_yl = 1; 
  if ((-1.25*magB2_yr[8])-0.9682458365518543*magB2_yr[7]-0.5590169943749475*magB2_yr[5]+1.118033988749895*magB2_yr[4]+0.8660254037844386*magB2_yr[1]+0.5*magB2_yr[0] < 0.0) cell_avg_yr = 1; 
  if ((-1.25*magB2_zl[8])-0.9682458365518543*magB2_zl[7]-0.5590169943749475*magB2_zl[5]+1.118033988749895*magB2_zl[4]+0.8660254037844386*magB2_zl[1]+0.5*magB2_zl[0] < 0.0) cell_avg_zl = 1; 
  if ((-1.25*magB2_zr[8])-0.9682458365518543*magB2_zr[7]-0.5590169943749475*magB2_zr[5]+1.118033988749895*magB2_zr[4]+0.8660254037844386*magB2_zr[1]+0.5*magB2_zr[0] < 0.0) cell_avg_zr = 1; 
  if (2.5*magB2_xl[8]-1.936491673103709*magB2_xl[7]+1.936491673103709*magB2_xl[6]+1.118033988749895*magB2_xl[5]+1.118033988749895*magB2_xl[4]-1.5*magB2_xl[3]+0.8660254037844386*magB2_xl[2]-0.8660254037844386*magB2_xl[1]+0.5*magB2_xl[0] < 0.0) cell_avg_xl = 1; 
  if (2.5*magB2_xr[8]-1.936491673103709*magB2_xr[7]+1.936491673103709*magB2_xr[6]+1.118033988749895*magB2_xr[5]+1.118033988749895*magB2_xr[4]-1.5*magB2_xr[3]+0.8660254037844386*magB2_xr[2]-0.8660254037844386*magB2_xr[1]+0.5*magB2_xr[0] < 0.0) cell_avg_xr = 1; 
  if (2.5*magB2_yl[8]-1.936491673103709*magB2_yl[7]+1.936491673103709*magB2_yl[6]+1.118033988749895*magB2_yl[5]+1.118033988749895*magB2_yl[4]-1.5*magB2_yl[3]+0.8660254037844386*magB2_yl[2]-0.8660254037844386*magB2_yl[1]+0.5*magB2_yl[0] < 0.0) cell_avg_yl = 1; 
  if (2.5*magB2_yr[8]-1.936491673103709*magB2_yr[7]+1.936491673103709*magB2_yr[6]+1.118033988749895*magB2_yr[5]+1.118033988749895*magB2_yr[4]-1.5*magB2_yr[3]+0.8660254037844386*magB2_yr[2]-0.8660254037844386*magB2_yr[1]+0.5*magB2_yr[0] < 0.0) cell_avg_yr = 1; 
  if (2.5*magB2_zl[8]-1.936491673103709*magB2_zl[7]+1.936491673103709*magB2_zl[6]+1.118033988749895*magB2_zl[5]+1.118033988749895*magB2_zl[4]-1.5*magB2_zl[3]+0.8660254037844386*magB2_zl[2]-0.8660254037844386*magB2_zl[1]+0.5*magB2_zl[0] < 0.0) cell_avg_zl = 1; 
  if (2.5*magB2_zr[8]-1.936491673103709*magB2_zr[7]+1.936491673103709*magB2_zr[6]+1.118033988749895*magB2_zr[5]+1.118033988749895*magB2_zr[4]-1.5*magB2_zr[3]+0.8660254037844386*magB2_zr[2]-0.8660254037844386*magB2_zr[1]+0.5*magB2_zr[0] < 0.0) cell_avg_zr = 1; 
  if ((-1.25*magB2_xl[8])-0.9682458365518543*magB2_xl[6]+1.118033988749895*magB2_xl[5]-0.5590169943749475*magB2_xl[4]+0.8660254037844386*magB2_xl[2]+0.5*magB2_xl[0] < 0.0) cell_avg_xl = 1; 
  if ((-1.25*magB2_xr[8])-0.9682458365518543*magB2_xr[6]+1.118033988749895*magB2_xr[5]-0.5590169943749475*magB2_xr[4]+0.8660254037844386*magB2_xr[2]+0.5*magB2_xr[0] < 0.0) cell_avg_xr = 1; 
  if ((-1.25*magB2_yl[8])-0.9682458365518543*magB2_yl[6]+1.118033988749895*magB2_yl[5]-0.5590169943749475*magB2_yl[4]+0.8660254037844386*magB2_yl[2]+0.5*magB2_yl[0] < 0.0) cell_avg_yl = 1; 
  if ((-1.25*magB2_yr[8])-0.9682458365518543*magB2_yr[6]+1.118033988749895*magB2_yr[5]-0.5590169943749475*magB2_yr[4]+0.8660254037844386*magB2_yr[2]+0.5*magB2_yr[0] < 0.0) cell_avg_yr = 1; 
  if ((-1.25*magB2_zl[8])-0.9682458365518543*magB2_zl[6]+1.118033988749895*magB2_zl[5]-0.5590169943749475*magB2_zl[4]+0.8660254037844386*magB2_zl[2]+0.5*magB2_zl[0] < 0.0) cell_avg_zl = 1; 
  if ((-1.25*magB2_zr[8])-0.9682458365518543*magB2_zr[6]+1.118033988749895*magB2_zr[5]-0.5590169943749475*magB2_zr[4]+0.8660254037844386*magB2_zr[2]+0.5*magB2_zr[0] < 0.0) cell_avg_zr = 1; 
  if (2.5*magB2_xl[8]+1.936491673103709*magB2_xl[7]+1.936491673103709*magB2_xl[6]+1.118033988749895*magB2_xl[5]+1.118033988749895*magB2_xl[4]+1.5*magB2_xl[3]+0.8660254037844386*magB2_xl[2]+0.8660254037844386*magB2_xl[1]+0.5*magB2_xl[0] < 0.0) cell_avg_xl = 1; 
  if (2.5*magB2_xr[8]+1.936491673103709*magB2_xr[7]+1.936491673103709*magB2_xr[6]+1.118033988749895*magB2_xr[5]+1.118033988749895*magB2_xr[4]+1.5*magB2_xr[3]+0.8660254037844386*magB2_xr[2]+0.8660254037844386*magB2_xr[1]+0.5*magB2_xr[0] < 0.0) cell_avg_xr = 1; 
  if (2.5*magB2_yl[8]+1.936491673103709*magB2_yl[7]+1.936491673103709*magB2_yl[6]+1.118033988749895*magB2_yl[5]+1.118033988749895*magB2_yl[4]+1.5*magB2_yl[3]+0.8660254037844386*magB2_yl[2]+0.8660254037844386*magB2_yl[1]+0.5*magB2_yl[0] < 0.0) cell_avg_yl = 1; 
  if (2.5*magB2_yr[8]+1.936491673103709*magB2_yr[7]+1.936491673103709*magB2_yr[6]+1.118033988749895*magB2_yr[5]+1.118033988749895*magB2_yr[4]+1.5*magB2_yr[3]+0.8660254037844386*magB2_yr[2]+0.8660254037844386*magB2_yr[1]+0.5*magB2_yr[0] < 0.0) cell_avg_yr = 1; 
  if (2.5*magB2_zl[8]+1.936491673103709*magB2_zl[7]+1.936491673103709*magB2_zl[6]+1.118033988749895*magB2_zl[5]+1.118033988749895*magB2_zl[4]+1.5*magB2_zl[3]+0.8660254037844386*magB2_zl[2]+0.8660254037844386*magB2_zl[1]+0.5*magB2_zl[0] < 0.0) cell_avg_zl = 1; 
  if (2.5*magB2_zr[8]+1.936491673103709*magB2_zr[7]+1.936491673103709*magB2_zr[6]+1.118033988749895*magB2_zr[5]+1.118033988749895*magB2_zr[4]+1.5*magB2_zr[3]+0.8660254037844386*magB2_zr[2]+0.8660254037844386*magB2_zr[1]+0.5*magB2_zr[0] < 0.0) cell_avg_zr = 1; 
 
  cell_avg_magB2_xl[0] = cell_avg_xl; 
  cell_avg_magB2_xr[0] = cell_avg_xr; 
  cell_avg_magB2_yl[0] = cell_avg_yl; 
  cell_avg_magB2_yr[0] = cell_avg_yr; 
  cell_avg_magB2_zl[0] = cell_avg_zl; 
  cell_avg_magB2_zr[0] = cell_avg_zr; 
 
  if (cell_avg_xl) { 
  magB2_xl[1] = 0.0; 
  magB2_xl[2] = 0.0; 
  magB2_xl[3] = 0.0; 
  magB2_xl[4] = 0.0; 
  magB2_xl[5] = 0.0; 
  magB2_xl[6] = 0.0; 
  magB2_xl[7] = 0.0; 
  magB2_xl[8] = 0.0; 
  } 
 
  if (cell_avg_xr) { 
  magB2_xr[1] = 0.0; 
  magB2_xr[2] = 0.0; 
  magB2_xr[3] = 0.0; 
  magB2_xr[4] = 0.0; 
  magB2_xr[5] = 0.0; 
  magB2_xr[6] = 0.0; 
  magB2_xr[7] = 0.0; 
  magB2_xr[8] = 0.0; 
  } 
 
  if (cell_avg_yl) { 
  magB2_yl[1] = 0.0; 
  magB2_yl[2] = 0.0; 
  magB2_yl[3] = 0.0; 
  magB2_yl[4] = 0.0; 
  magB2_yl[5] = 0.0; 
  magB2_yl[6] = 0.0; 
  magB2_yl[7] = 0.0; 
  magB2_yl[8] = 0.0; 
  } 
 
  if (cell_avg_yr) { 
  magB2_yr[1] = 0.0; 
  magB2_yr[2] = 0.0; 
  magB2_yr[3] = 0.0; 
  magB2_yr[4] = 0.0; 
  magB2_yr[5] = 0.0; 
  magB2_yr[6] = 0.0; 
  magB2_yr[7] = 0.0; 
  magB2_yr[8] = 0.0; 
  } 
 
  if (cell_avg_zl) { 
  magB2_zl[1] = 0.0; 
  magB2_zl[2] = 0.0; 
  magB2_zl[3] = 0.0; 
  magB2_zl[4] = 0.0; 
  magB2_zl[5] = 0.0; 
  magB2_zl[6] = 0.0; 
  magB2_zl[7] = 0.0; 
  magB2_zl[8] = 0.0; 
  } 
 
  if (cell_avg_zr) { 
  magB2_zr[1] = 0.0; 
  magB2_zr[2] = 0.0; 
  magB2_zr[3] = 0.0; 
  magB2_zr[4] = 0.0; 
  magB2_zr[5] = 0.0; 
  magB2_zr[6] = 0.0; 
  magB2_zr[7] = 0.0; 
  magB2_zr[8] = 0.0; 
  } 
 
  gkyl_mat_set(&rhs_bxbx_xl,0,0,Bx_sq_xl[0]); 
  gkyl_mat_set(&rhs_bxbx_xr,0,0,Bx_sq_xr[0]); 
 
  gkyl_mat_set(&rhs_byby_yl,0,0,By_sq_yl[0]); 
  gkyl_mat_set(&rhs_byby_yr,0,0,By_sq_yr[0]); 
 
  gkyl_mat_set(&rhs_bzbz_zl,0,0,Bz_sq_zl[0]); 
  gkyl_mat_set(&rhs_bzbz_zr,0,0,Bz_sq_zr[0]); 
 
  gkyl_mat_set(&rhs_bxbx_xl,1,0,Bx_sq_xl[1]); 
  gkyl_mat_set(&rhs_bxbx_xr,1,0,Bx_sq_xr[1]); 
 
  gkyl_mat_set(&rhs_byby_yl,1,0,By_sq_yl[1]); 
  gkyl_mat_set(&rhs_byby_yr,1,0,By_sq_yr[1]); 
 
  gkyl_mat_set(&rhs_bzbz_zl,1,0,Bz_sq_zl[1]); 
  gkyl_mat_set(&rhs_bzbz_zr,1,0,Bz_sq_zr[1]); 
 
  gkyl_mat_set(&rhs_bxbx_xl,2,0,Bx_sq_xl[2]); 
  gkyl_mat_set(&rhs_bxbx_xr,2,0,Bx_sq_xr[2]); 
 
  gkyl_mat_set(&rhs_byby_yl,2,0,By_sq_yl[2]); 
  gkyl_mat_set(&rhs_byby_yr,2,0,By_sq_yr[2]); 
 
  gkyl_mat_set(&rhs_bzbz_zl,2,0,Bz_sq_zl[2]); 
  gkyl_mat_set(&rhs_bzbz_zr,2,0,Bz_sq_zr[2]); 
 
  gkyl_mat_set(&rhs_bxbx_xl,3,0,Bx_sq_xl[3]); 
  gkyl_mat_set(&rhs_bxbx_xr,3,0,Bx_sq_xr[3]); 
 
  gkyl_mat_set(&rhs_byby_yl,3,0,By_sq_yl[3]); 
  gkyl_mat_set(&rhs_byby_yr,3,0,By_sq_yr[3]); 
 
  gkyl_mat_set(&rhs_bzbz_zl,3,0,Bz_sq_zl[3]); 
  gkyl_mat_set(&rhs_bzbz_zr,3,0,Bz_sq_zr[3]); 
 
  gkyl_mat_set(&rhs_bxbx_xl,4,0,Bx_sq_xl[4]); 
  gkyl_mat_set(&rhs_bxbx_xr,4,0,Bx_sq_xr[4]); 
 
  gkyl_mat_set(&rhs_byby_yl,4,0,By_sq_yl[4]); 
  gkyl_mat_set(&rhs_byby_yr,4,0,By_sq_yr[4]); 
 
  gkyl_mat_set(&rhs_bzbz_zl,4,0,Bz_sq_zl[4]); 
  gkyl_mat_set(&rhs_bzbz_zr,4,0,Bz_sq_zr[4]); 
 
  gkyl_mat_set(&rhs_bxbx_xl,5,0,Bx_sq_xl[5]); 
  gkyl_mat_set(&rhs_bxbx_xr,5,0,Bx_sq_xr[5]); 
 
  gkyl_mat_set(&rhs_byby_yl,5,0,By_sq_yl[5]); 
  gkyl_mat_set(&rhs_byby_yr,5,0,By_sq_yr[5]); 
 
  gkyl_mat_set(&rhs_bzbz_zl,5,0,Bz_sq_zl[5]); 
  gkyl_mat_set(&rhs_bzbz_zr,5,0,Bz_sq_zr[5]); 
 
  gkyl_mat_set(&rhs_bxbx_xl,6,0,Bx_sq_xl[6]); 
  gkyl_mat_set(&rhs_bxbx_xr,6,0,Bx_sq_xr[6]); 
 
  gkyl_mat_set(&rhs_byby_yl,6,0,By_sq_yl[6]); 
  gkyl_mat_set(&rhs_byby_yr,6,0,By_sq_yr[6]); 
 
  gkyl_mat_set(&rhs_bzbz_zl,6,0,Bz_sq_zl[6]); 
  gkyl_mat_set(&rhs_bzbz_zr,6,0,Bz_sq_zr[6]); 
 
  gkyl_mat_set(&rhs_bxbx_xl,7,0,Bx_sq_xl[7]); 
  gkyl_mat_set(&rhs_bxbx_xr,7,0,Bx_sq_xr[7]); 
 
  gkyl_mat_set(&rhs_byby_yl,7,0,By_sq_yl[7]); 
  gkyl_mat_set(&rhs_byby_yr,7,0,By_sq_yr[7]); 
 
  gkyl_mat_set(&rhs_bzbz_zl,7,0,Bz_sq_zl[7]); 
  gkyl_mat_set(&rhs_bzbz_zr,7,0,Bz_sq_zr[7]); 
 
  gkyl_mat_set(&rhs_bxbx_xl,8,0,Bx_sq_xl[8]); 
  gkyl_mat_set(&rhs_bxbx_xr,8,0,Bx_sq_xr[8]); 
 
  gkyl_mat_set(&rhs_byby_yl,8,0,By_sq_yl[8]); 
  gkyl_mat_set(&rhs_byby_yr,8,0,By_sq_yr[8]); 
 
  gkyl_mat_set(&rhs_bzbz_zl,8,0,Bz_sq_zl[8]); 
  gkyl_mat_set(&rhs_bzbz_zr,8,0,Bz_sq_zr[8]); 
 
  double temp_magB2_xl = 0.0; 
  double temp_magB2_xr = 0.0; 
  double temp_magB2_yl = 0.0; 
  double temp_magB2_yr = 0.0; 
  double temp_magB2_zl = 0.0; 
  double temp_magB2_zr = 0.0; 
  temp_magB2_xl = 0.5*magB2_xl[0]; 
  gkyl_mat_set(&A_bxbx_xl,0,0,temp_magB2_xl); 
 
  temp_magB2_xr = 0.5*magB2_xr[0]; 
  gkyl_mat_set(&A_bxbx_xr,0,0,temp_magB2_xr); 
 
  temp_magB2_yl = 0.5*magB2_yl[0]; 
  gkyl_mat_set(&A_byby_yl,0,0,temp_magB2_yl); 
 
  temp_magB2_yr = 0.5*magB2_yr[0]; 
  gkyl_mat_set(&A_byby_yr,0,0,temp_magB2_yr); 
 
  temp_magB2_zl = 0.5*magB2_zl[0]; 
  gkyl_mat_set(&A_bzbz_zl,0,0,temp_magB2_zl); 
 
  temp_magB2_zr = 0.5*magB2_zr[0]; 
  gkyl_mat_set(&A_bzbz_zr,0,0,temp_magB2_zr); 
 
  temp_magB2_xl = 0.5*magB2_xl[1]; 
  gkyl_mat_set(&A_bxbx_xl,0,1,temp_magB2_xl); 
 
  temp_magB2_xr = 0.5*magB2_xr[1]; 
  gkyl_mat_set(&A_bxbx_xr,0,1,temp_magB2_xr); 
 
  temp_magB2_yl = 0.5*magB2_yl[1]; 
  gkyl_mat_set(&A_byby_yl,0,1,temp_magB2_yl); 
 
  temp_magB2_yr = 0.5*magB2_yr[1]; 
  gkyl_mat_set(&A_byby_yr,0,1,temp_magB2_yr); 
 
  temp_magB2_zl = 0.5*magB2_zl[1]; 
  gkyl_mat_set(&A_bzbz_zl,0,1,temp_magB2_zl); 
 
  temp_magB2_zr = 0.5*magB2_zr[1]; 
  gkyl_mat_set(&A_bzbz_zr,0,1,temp_magB2_zr); 
 
  temp_magB2_xl = 0.5*magB2_xl[2]; 
  gkyl_mat_set(&A_bxbx_xl,0,2,temp_magB2_xl); 
 
  temp_magB2_xr = 0.5*magB2_xr[2]; 
  gkyl_mat_set(&A_bxbx_xr,0,2,temp_magB2_xr); 
 
  temp_magB2_yl = 0.5*magB2_yl[2]; 
  gkyl_mat_set(&A_byby_yl,0,2,temp_magB2_yl); 
 
  temp_magB2_yr = 0.5*magB2_yr[2]; 
  gkyl_mat_set(&A_byby_yr,0,2,temp_magB2_yr); 
 
  temp_magB2_zl = 0.5*magB2_zl[2]; 
  gkyl_mat_set(&A_bzbz_zl,0,2,temp_magB2_zl); 
 
  temp_magB2_zr = 0.5*magB2_zr[2]; 
  gkyl_mat_set(&A_bzbz_zr,0,2,temp_magB2_zr); 
 
  temp_magB2_xl = 0.5*magB2_xl[3]; 
  gkyl_mat_set(&A_bxbx_xl,0,3,temp_magB2_xl); 
 
  temp_magB2_xr = 0.5*magB2_xr[3]; 
  gkyl_mat_set(&A_bxbx_xr,0,3,temp_magB2_xr); 
 
  temp_magB2_yl = 0.5*magB2_yl[3]; 
  gkyl_mat_set(&A_byby_yl,0,3,temp_magB2_yl); 
 
  temp_magB2_yr = 0.5*magB2_yr[3]; 
  gkyl_mat_set(&A_byby_yr,0,3,temp_magB2_yr); 
 
  temp_magB2_zl = 0.5*magB2_zl[3]; 
  gkyl_mat_set(&A_bzbz_zl,0,3,temp_magB2_zl); 
 
  temp_magB2_zr = 0.5*magB2_zr[3]; 
  gkyl_mat_set(&A_bzbz_zr,0,3,temp_magB2_zr); 
 
  temp_magB2_xl = 0.5*magB2_xl[4]; 
  gkyl_mat_set(&A_bxbx_xl,0,4,temp_magB2_xl); 
 
  temp_magB2_xr = 0.5*magB2_xr[4]; 
  gkyl_mat_set(&A_bxbx_xr,0,4,temp_magB2_xr); 
 
  temp_magB2_yl = 0.5*magB2_yl[4]; 
  gkyl_mat_set(&A_byby_yl,0,4,temp_magB2_yl); 
 
  temp_magB2_yr = 0.5*magB2_yr[4]; 
  gkyl_mat_set(&A_byby_yr,0,4,temp_magB2_yr); 
 
  temp_magB2_zl = 0.5*magB2_zl[4]; 
  gkyl_mat_set(&A_bzbz_zl,0,4,temp_magB2_zl); 
 
  temp_magB2_zr = 0.5*magB2_zr[4]; 
  gkyl_mat_set(&A_bzbz_zr,0,4,temp_magB2_zr); 
 
  temp_magB2_xl = 0.5*magB2_xl[5]; 
  gkyl_mat_set(&A_bxbx_xl,0,5,temp_magB2_xl); 
 
  temp_magB2_xr = 0.5*magB2_xr[5]; 
  gkyl_mat_set(&A_bxbx_xr,0,5,temp_magB2_xr); 
 
  temp_magB2_yl = 0.5*magB2_yl[5]; 
  gkyl_mat_set(&A_byby_yl,0,5,temp_magB2_yl); 
 
  temp_magB2_yr = 0.5*magB2_yr[5]; 
  gkyl_mat_set(&A_byby_yr,0,5,temp_magB2_yr); 
 
  temp_magB2_zl = 0.5*magB2_zl[5]; 
  gkyl_mat_set(&A_bzbz_zl,0,5,temp_magB2_zl); 
 
  temp_magB2_zr = 0.5*magB2_zr[5]; 
  gkyl_mat_set(&A_bzbz_zr,0,5,temp_magB2_zr); 
 
  temp_magB2_xl = 0.5*magB2_xl[6]; 
  gkyl_mat_set(&A_bxbx_xl,0,6,temp_magB2_xl); 
 
  temp_magB2_xr = 0.5*magB2_xr[6]; 
  gkyl_mat_set(&A_bxbx_xr,0,6,temp_magB2_xr); 
 
  temp_magB2_yl = 0.5*magB2_yl[6]; 
  gkyl_mat_set(&A_byby_yl,0,6,temp_magB2_yl); 
 
  temp_magB2_yr = 0.5*magB2_yr[6]; 
  gkyl_mat_set(&A_byby_yr,0,6,temp_magB2_yr); 
 
  temp_magB2_zl = 0.5*magB2_zl[6]; 
  gkyl_mat_set(&A_bzbz_zl,0,6,temp_magB2_zl); 
 
  temp_magB2_zr = 0.5*magB2_zr[6]; 
  gkyl_mat_set(&A_bzbz_zr,0,6,temp_magB2_zr); 
 
  temp_magB2_xl = 0.5*magB2_xl[7]; 
  gkyl_mat_set(&A_bxbx_xl,0,7,temp_magB2_xl); 
 
  temp_magB2_xr = 0.5*magB2_xr[7]; 
  gkyl_mat_set(&A_bxbx_xr,0,7,temp_magB2_xr); 
 
  temp_magB2_yl = 0.5*magB2_yl[7]; 
  gkyl_mat_set(&A_byby_yl,0,7,temp_magB2_yl); 
 
  temp_magB2_yr = 0.5*magB2_yr[7]; 
  gkyl_mat_set(&A_byby_yr,0,7,temp_magB2_yr); 
 
  temp_magB2_zl = 0.5*magB2_zl[7]; 
  gkyl_mat_set(&A_bzbz_zl,0,7,temp_magB2_zl); 
 
  temp_magB2_zr = 0.5*magB2_zr[7]; 
  gkyl_mat_set(&A_bzbz_zr,0,7,temp_magB2_zr); 
 
  temp_magB2_xl = 0.5*magB2_xl[8]; 
  gkyl_mat_set(&A_bxbx_xl,0,8,temp_magB2_xl); 
 
  temp_magB2_xr = 0.5*magB2_xr[8]; 
  gkyl_mat_set(&A_bxbx_xr,0,8,temp_magB2_xr); 
 
  temp_magB2_yl = 0.5*magB2_yl[8]; 
  gkyl_mat_set(&A_byby_yl,0,8,temp_magB2_yl); 
 
  temp_magB2_yr = 0.5*magB2_yr[8]; 
  gkyl_mat_set(&A_byby_yr,0,8,temp_magB2_yr); 
 
  temp_magB2_zl = 0.5*magB2_zl[8]; 
  gkyl_mat_set(&A_bzbz_zl,0,8,temp_magB2_zl); 
 
  temp_magB2_zr = 0.5*magB2_zr[8]; 
  gkyl_mat_set(&A_bzbz_zr,0,8,temp_magB2_zr); 
 
  temp_magB2_xl = 0.5*magB2_xl[1]; 
  gkyl_mat_set(&A_bxbx_xl,1,0,temp_magB2_xl); 
 
  temp_magB2_xr = 0.5*magB2_xr[1]; 
  gkyl_mat_set(&A_bxbx_xr,1,0,temp_magB2_xr); 
 
  temp_magB2_yl = 0.5*magB2_yl[1]; 
  gkyl_mat_set(&A_byby_yl,1,0,temp_magB2_yl); 
 
  temp_magB2_yr = 0.5*magB2_yr[1]; 
  gkyl_mat_set(&A_byby_yr,1,0,temp_magB2_yr); 
 
  temp_magB2_zl = 0.5*magB2_zl[1]; 
  gkyl_mat_set(&A_bzbz_zl,1,0,temp_magB2_zl); 
 
  temp_magB2_zr = 0.5*magB2_zr[1]; 
  gkyl_mat_set(&A_bzbz_zr,1,0,temp_magB2_zr); 
 
  temp_magB2_xl = 0.4472135954999579*magB2_xl[4]+0.5*magB2_xl[0]; 
  gkyl_mat_set(&A_bxbx_xl,1,1,temp_magB2_xl); 
 
  temp_magB2_xr = 0.4472135954999579*magB2_xr[4]+0.5*magB2_xr[0]; 
  gkyl_mat_set(&A_bxbx_xr,1,1,temp_magB2_xr); 
 
  temp_magB2_yl = 0.4472135954999579*magB2_yl[4]+0.5*magB2_yl[0]; 
  gkyl_mat_set(&A_byby_yl,1,1,temp_magB2_yl); 
 
  temp_magB2_yr = 0.4472135954999579*magB2_yr[4]+0.5*magB2_yr[0]; 
  gkyl_mat_set(&A_byby_yr,1,1,temp_magB2_yr); 
 
  temp_magB2_zl = 0.4472135954999579*magB2_zl[4]+0.5*magB2_zl[0]; 
  gkyl_mat_set(&A_bzbz_zl,1,1,temp_magB2_zl); 
 
  temp_magB2_zr = 0.4472135954999579*magB2_zr[4]+0.5*magB2_zr[0]; 
  gkyl_mat_set(&A_bzbz_zr,1,1,temp_magB2_zr); 
 
  temp_magB2_xl = 0.5*magB2_xl[3]; 
  gkyl_mat_set(&A_bxbx_xl,1,2,temp_magB2_xl); 
 
  temp_magB2_xr = 0.5*magB2_xr[3]; 
  gkyl_mat_set(&A_bxbx_xr,1,2,temp_magB2_xr); 
 
  temp_magB2_yl = 0.5*magB2_yl[3]; 
  gkyl_mat_set(&A_byby_yl,1,2,temp_magB2_yl); 
 
  temp_magB2_yr = 0.5*magB2_yr[3]; 
  gkyl_mat_set(&A_byby_yr,1,2,temp_magB2_yr); 
 
  temp_magB2_zl = 0.5*magB2_zl[3]; 
  gkyl_mat_set(&A_bzbz_zl,1,2,temp_magB2_zl); 
 
  temp_magB2_zr = 0.5*magB2_zr[3]; 
  gkyl_mat_set(&A_bzbz_zr,1,2,temp_magB2_zr); 
 
  temp_magB2_xl = 0.447213595499958*magB2_xl[6]+0.5*magB2_xl[2]; 
  gkyl_mat_set(&A_bxbx_xl,1,3,temp_magB2_xl); 
 
  temp_magB2_xr = 0.447213595499958*magB2_xr[6]+0.5*magB2_xr[2]; 
  gkyl_mat_set(&A_bxbx_xr,1,3,temp_magB2_xr); 
 
  temp_magB2_yl = 0.447213595499958*magB2_yl[6]+0.5*magB2_yl[2]; 
  gkyl_mat_set(&A_byby_yl,1,3,temp_magB2_yl); 
 
  temp_magB2_yr = 0.447213595499958*magB2_yr[6]+0.5*magB2_yr[2]; 
  gkyl_mat_set(&A_byby_yr,1,3,temp_magB2_yr); 
 
  temp_magB2_zl = 0.447213595499958*magB2_zl[6]+0.5*magB2_zl[2]; 
  gkyl_mat_set(&A_bzbz_zl,1,3,temp_magB2_zl); 
 
  temp_magB2_zr = 0.447213595499958*magB2_zr[6]+0.5*magB2_zr[2]; 
  gkyl_mat_set(&A_bzbz_zr,1,3,temp_magB2_zr); 
 
  temp_magB2_xl = 0.4472135954999579*magB2_xl[1]; 
  gkyl_mat_set(&A_bxbx_xl,1,4,temp_magB2_xl); 
 
  temp_magB2_xr = 0.4472135954999579*magB2_xr[1]; 
  gkyl_mat_set(&A_bxbx_xr,1,4,temp_magB2_xr); 
 
  temp_magB2_yl = 0.4472135954999579*magB2_yl[1]; 
  gkyl_mat_set(&A_byby_yl,1,4,temp_magB2_yl); 
 
  temp_magB2_yr = 0.4472135954999579*magB2_yr[1]; 
  gkyl_mat_set(&A_byby_yr,1,4,temp_magB2_yr); 
 
  temp_magB2_zl = 0.4472135954999579*magB2_zl[1]; 
  gkyl_mat_set(&A_bzbz_zl,1,4,temp_magB2_zl); 
 
  temp_magB2_zr = 0.4472135954999579*magB2_zr[1]; 
  gkyl_mat_set(&A_bzbz_zr,1,4,temp_magB2_zr); 
 
  temp_magB2_xl = 0.5000000000000001*magB2_xl[7]; 
  gkyl_mat_set(&A_bxbx_xl,1,5,temp_magB2_xl); 
 
  temp_magB2_xr = 0.5000000000000001*magB2_xr[7]; 
  gkyl_mat_set(&A_bxbx_xr,1,5,temp_magB2_xr); 
 
  temp_magB2_yl = 0.5000000000000001*magB2_yl[7]; 
  gkyl_mat_set(&A_byby_yl,1,5,temp_magB2_yl); 
 
  temp_magB2_yr = 0.5000000000000001*magB2_yr[7]; 
  gkyl_mat_set(&A_byby_yr,1,5,temp_magB2_yr); 
 
  temp_magB2_zl = 0.5000000000000001*magB2_zl[7]; 
  gkyl_mat_set(&A_bzbz_zl,1,5,temp_magB2_zl); 
 
  temp_magB2_zr = 0.5000000000000001*magB2_zr[7]; 
  gkyl_mat_set(&A_bzbz_zr,1,5,temp_magB2_zr); 
 
  temp_magB2_xl = 0.447213595499958*magB2_xl[3]; 
  gkyl_mat_set(&A_bxbx_xl,1,6,temp_magB2_xl); 
 
  temp_magB2_xr = 0.447213595499958*magB2_xr[3]; 
  gkyl_mat_set(&A_bxbx_xr,1,6,temp_magB2_xr); 
 
  temp_magB2_yl = 0.447213595499958*magB2_yl[3]; 
  gkyl_mat_set(&A_byby_yl,1,6,temp_magB2_yl); 
 
  temp_magB2_yr = 0.447213595499958*magB2_yr[3]; 
  gkyl_mat_set(&A_byby_yr,1,6,temp_magB2_yr); 
 
  temp_magB2_zl = 0.447213595499958*magB2_zl[3]; 
  gkyl_mat_set(&A_bzbz_zl,1,6,temp_magB2_zl); 
 
  temp_magB2_zr = 0.447213595499958*magB2_zr[3]; 
  gkyl_mat_set(&A_bzbz_zr,1,6,temp_magB2_zr); 
 
  temp_magB2_xl = 0.447213595499958*magB2_xl[8]+0.5000000000000001*magB2_xl[5]; 
  gkyl_mat_set(&A_bxbx_xl,1,7,temp_magB2_xl); 
 
  temp_magB2_xr = 0.447213595499958*magB2_xr[8]+0.5000000000000001*magB2_xr[5]; 
  gkyl_mat_set(&A_bxbx_xr,1,7,temp_magB2_xr); 
 
  temp_magB2_yl = 0.447213595499958*magB2_yl[8]+0.5000000000000001*magB2_yl[5]; 
  gkyl_mat_set(&A_byby_yl,1,7,temp_magB2_yl); 
 
  temp_magB2_yr = 0.447213595499958*magB2_yr[8]+0.5000000000000001*magB2_yr[5]; 
  gkyl_mat_set(&A_byby_yr,1,7,temp_magB2_yr); 
 
  temp_magB2_zl = 0.447213595499958*magB2_zl[8]+0.5000000000000001*magB2_zl[5]; 
  gkyl_mat_set(&A_bzbz_zl,1,7,temp_magB2_zl); 
 
  temp_magB2_zr = 0.447213595499958*magB2_zr[8]+0.5000000000000001*magB2_zr[5]; 
  gkyl_mat_set(&A_bzbz_zr,1,7,temp_magB2_zr); 
 
  temp_magB2_xl = 0.447213595499958*magB2_xl[7]; 
  gkyl_mat_set(&A_bxbx_xl,1,8,temp_magB2_xl); 
 
  temp_magB2_xr = 0.447213595499958*magB2_xr[7]; 
  gkyl_mat_set(&A_bxbx_xr,1,8,temp_magB2_xr); 
 
  temp_magB2_yl = 0.447213595499958*magB2_yl[7]; 
  gkyl_mat_set(&A_byby_yl,1,8,temp_magB2_yl); 
 
  temp_magB2_yr = 0.447213595499958*magB2_yr[7]; 
  gkyl_mat_set(&A_byby_yr,1,8,temp_magB2_yr); 
 
  temp_magB2_zl = 0.447213595499958*magB2_zl[7]; 
  gkyl_mat_set(&A_bzbz_zl,1,8,temp_magB2_zl); 
 
  temp_magB2_zr = 0.447213595499958*magB2_zr[7]; 
  gkyl_mat_set(&A_bzbz_zr,1,8,temp_magB2_zr); 
 
  temp_magB2_xl = 0.5*magB2_xl[2]; 
  gkyl_mat_set(&A_bxbx_xl,2,0,temp_magB2_xl); 
 
  temp_magB2_xr = 0.5*magB2_xr[2]; 
  gkyl_mat_set(&A_bxbx_xr,2,0,temp_magB2_xr); 
 
  temp_magB2_yl = 0.5*magB2_yl[2]; 
  gkyl_mat_set(&A_byby_yl,2,0,temp_magB2_yl); 
 
  temp_magB2_yr = 0.5*magB2_yr[2]; 
  gkyl_mat_set(&A_byby_yr,2,0,temp_magB2_yr); 
 
  temp_magB2_zl = 0.5*magB2_zl[2]; 
  gkyl_mat_set(&A_bzbz_zl,2,0,temp_magB2_zl); 
 
  temp_magB2_zr = 0.5*magB2_zr[2]; 
  gkyl_mat_set(&A_bzbz_zr,2,0,temp_magB2_zr); 
 
  temp_magB2_xl = 0.5*magB2_xl[3]; 
  gkyl_mat_set(&A_bxbx_xl,2,1,temp_magB2_xl); 
 
  temp_magB2_xr = 0.5*magB2_xr[3]; 
  gkyl_mat_set(&A_bxbx_xr,2,1,temp_magB2_xr); 
 
  temp_magB2_yl = 0.5*magB2_yl[3]; 
  gkyl_mat_set(&A_byby_yl,2,1,temp_magB2_yl); 
 
  temp_magB2_yr = 0.5*magB2_yr[3]; 
  gkyl_mat_set(&A_byby_yr,2,1,temp_magB2_yr); 
 
  temp_magB2_zl = 0.5*magB2_zl[3]; 
  gkyl_mat_set(&A_bzbz_zl,2,1,temp_magB2_zl); 
 
  temp_magB2_zr = 0.5*magB2_zr[3]; 
  gkyl_mat_set(&A_bzbz_zr,2,1,temp_magB2_zr); 
 
  temp_magB2_xl = 0.4472135954999579*magB2_xl[5]+0.5*magB2_xl[0]; 
  gkyl_mat_set(&A_bxbx_xl,2,2,temp_magB2_xl); 
 
  temp_magB2_xr = 0.4472135954999579*magB2_xr[5]+0.5*magB2_xr[0]; 
  gkyl_mat_set(&A_bxbx_xr,2,2,temp_magB2_xr); 
 
  temp_magB2_yl = 0.4472135954999579*magB2_yl[5]+0.5*magB2_yl[0]; 
  gkyl_mat_set(&A_byby_yl,2,2,temp_magB2_yl); 
 
  temp_magB2_yr = 0.4472135954999579*magB2_yr[5]+0.5*magB2_yr[0]; 
  gkyl_mat_set(&A_byby_yr,2,2,temp_magB2_yr); 
 
  temp_magB2_zl = 0.4472135954999579*magB2_zl[5]+0.5*magB2_zl[0]; 
  gkyl_mat_set(&A_bzbz_zl,2,2,temp_magB2_zl); 
 
  temp_magB2_zr = 0.4472135954999579*magB2_zr[5]+0.5*magB2_zr[0]; 
  gkyl_mat_set(&A_bzbz_zr,2,2,temp_magB2_zr); 
 
  temp_magB2_xl = 0.447213595499958*magB2_xl[7]+0.5*magB2_xl[1]; 
  gkyl_mat_set(&A_bxbx_xl,2,3,temp_magB2_xl); 
 
  temp_magB2_xr = 0.447213595499958*magB2_xr[7]+0.5*magB2_xr[1]; 
  gkyl_mat_set(&A_bxbx_xr,2,3,temp_magB2_xr); 
 
  temp_magB2_yl = 0.447213595499958*magB2_yl[7]+0.5*magB2_yl[1]; 
  gkyl_mat_set(&A_byby_yl,2,3,temp_magB2_yl); 
 
  temp_magB2_yr = 0.447213595499958*magB2_yr[7]+0.5*magB2_yr[1]; 
  gkyl_mat_set(&A_byby_yr,2,3,temp_magB2_yr); 
 
  temp_magB2_zl = 0.447213595499958*magB2_zl[7]+0.5*magB2_zl[1]; 
  gkyl_mat_set(&A_bzbz_zl,2,3,temp_magB2_zl); 
 
  temp_magB2_zr = 0.447213595499958*magB2_zr[7]+0.5*magB2_zr[1]; 
  gkyl_mat_set(&A_bzbz_zr,2,3,temp_magB2_zr); 
 
  temp_magB2_xl = 0.5000000000000001*magB2_xl[6]; 
  gkyl_mat_set(&A_bxbx_xl,2,4,temp_magB2_xl); 
 
  temp_magB2_xr = 0.5000000000000001*magB2_xr[6]; 
  gkyl_mat_set(&A_bxbx_xr,2,4,temp_magB2_xr); 
 
  temp_magB2_yl = 0.5000000000000001*magB2_yl[6]; 
  gkyl_mat_set(&A_byby_yl,2,4,temp_magB2_yl); 
 
  temp_magB2_yr = 0.5000000000000001*magB2_yr[6]; 
  gkyl_mat_set(&A_byby_yr,2,4,temp_magB2_yr); 
 
  temp_magB2_zl = 0.5000000000000001*magB2_zl[6]; 
  gkyl_mat_set(&A_bzbz_zl,2,4,temp_magB2_zl); 
 
  temp_magB2_zr = 0.5000000000000001*magB2_zr[6]; 
  gkyl_mat_set(&A_bzbz_zr,2,4,temp_magB2_zr); 
 
  temp_magB2_xl = 0.4472135954999579*magB2_xl[2]; 
  gkyl_mat_set(&A_bxbx_xl,2,5,temp_magB2_xl); 
 
  temp_magB2_xr = 0.4472135954999579*magB2_xr[2]; 
  gkyl_mat_set(&A_bxbx_xr,2,5,temp_magB2_xr); 
 
  temp_magB2_yl = 0.4472135954999579*magB2_yl[2]; 
  gkyl_mat_set(&A_byby_yl,2,5,temp_magB2_yl); 
 
  temp_magB2_yr = 0.4472135954999579*magB2_yr[2]; 
  gkyl_mat_set(&A_byby_yr,2,5,temp_magB2_yr); 
 
  temp_magB2_zl = 0.4472135954999579*magB2_zl[2]; 
  gkyl_mat_set(&A_bzbz_zl,2,5,temp_magB2_zl); 
 
  temp_magB2_zr = 0.4472135954999579*magB2_zr[2]; 
  gkyl_mat_set(&A_bzbz_zr,2,5,temp_magB2_zr); 
 
  temp_magB2_xl = 0.447213595499958*magB2_xl[8]+0.5000000000000001*magB2_xl[4]; 
  gkyl_mat_set(&A_bxbx_xl,2,6,temp_magB2_xl); 
 
  temp_magB2_xr = 0.447213595499958*magB2_xr[8]+0.5000000000000001*magB2_xr[4]; 
  gkyl_mat_set(&A_bxbx_xr,2,6,temp_magB2_xr); 
 
  temp_magB2_yl = 0.447213595499958*magB2_yl[8]+0.5000000000000001*magB2_yl[4]; 
  gkyl_mat_set(&A_byby_yl,2,6,temp_magB2_yl); 
 
  temp_magB2_yr = 0.447213595499958*magB2_yr[8]+0.5000000000000001*magB2_yr[4]; 
  gkyl_mat_set(&A_byby_yr,2,6,temp_magB2_yr); 
 
  temp_magB2_zl = 0.447213595499958*magB2_zl[8]+0.5000000000000001*magB2_zl[4]; 
  gkyl_mat_set(&A_bzbz_zl,2,6,temp_magB2_zl); 
 
  temp_magB2_zr = 0.447213595499958*magB2_zr[8]+0.5000000000000001*magB2_zr[4]; 
  gkyl_mat_set(&A_bzbz_zr,2,6,temp_magB2_zr); 
 
  temp_magB2_xl = 0.447213595499958*magB2_xl[3]; 
  gkyl_mat_set(&A_bxbx_xl,2,7,temp_magB2_xl); 
 
  temp_magB2_xr = 0.447213595499958*magB2_xr[3]; 
  gkyl_mat_set(&A_bxbx_xr,2,7,temp_magB2_xr); 
 
  temp_magB2_yl = 0.447213595499958*magB2_yl[3]; 
  gkyl_mat_set(&A_byby_yl,2,7,temp_magB2_yl); 
 
  temp_magB2_yr = 0.447213595499958*magB2_yr[3]; 
  gkyl_mat_set(&A_byby_yr,2,7,temp_magB2_yr); 
 
  temp_magB2_zl = 0.447213595499958*magB2_zl[3]; 
  gkyl_mat_set(&A_bzbz_zl,2,7,temp_magB2_zl); 
 
  temp_magB2_zr = 0.447213595499958*magB2_zr[3]; 
  gkyl_mat_set(&A_bzbz_zr,2,7,temp_magB2_zr); 
 
  temp_magB2_xl = 0.447213595499958*magB2_xl[6]; 
  gkyl_mat_set(&A_bxbx_xl,2,8,temp_magB2_xl); 
 
  temp_magB2_xr = 0.447213595499958*magB2_xr[6]; 
  gkyl_mat_set(&A_bxbx_xr,2,8,temp_magB2_xr); 
 
  temp_magB2_yl = 0.447213595499958*magB2_yl[6]; 
  gkyl_mat_set(&A_byby_yl,2,8,temp_magB2_yl); 
 
  temp_magB2_yr = 0.447213595499958*magB2_yr[6]; 
  gkyl_mat_set(&A_byby_yr,2,8,temp_magB2_yr); 
 
  temp_magB2_zl = 0.447213595499958*magB2_zl[6]; 
  gkyl_mat_set(&A_bzbz_zl,2,8,temp_magB2_zl); 
 
  temp_magB2_zr = 0.447213595499958*magB2_zr[6]; 
  gkyl_mat_set(&A_bzbz_zr,2,8,temp_magB2_zr); 
 
  temp_magB2_xl = 0.5*magB2_xl[3]; 
  gkyl_mat_set(&A_bxbx_xl,3,0,temp_magB2_xl); 
 
  temp_magB2_xr = 0.5*magB2_xr[3]; 
  gkyl_mat_set(&A_bxbx_xr,3,0,temp_magB2_xr); 
 
  temp_magB2_yl = 0.5*magB2_yl[3]; 
  gkyl_mat_set(&A_byby_yl,3,0,temp_magB2_yl); 
 
  temp_magB2_yr = 0.5*magB2_yr[3]; 
  gkyl_mat_set(&A_byby_yr,3,0,temp_magB2_yr); 
 
  temp_magB2_zl = 0.5*magB2_zl[3]; 
  gkyl_mat_set(&A_bzbz_zl,3,0,temp_magB2_zl); 
 
  temp_magB2_zr = 0.5*magB2_zr[3]; 
  gkyl_mat_set(&A_bzbz_zr,3,0,temp_magB2_zr); 
 
  temp_magB2_xl = 0.447213595499958*magB2_xl[6]+0.5*magB2_xl[2]; 
  gkyl_mat_set(&A_bxbx_xl,3,1,temp_magB2_xl); 
 
  temp_magB2_xr = 0.447213595499958*magB2_xr[6]+0.5*magB2_xr[2]; 
  gkyl_mat_set(&A_bxbx_xr,3,1,temp_magB2_xr); 
 
  temp_magB2_yl = 0.447213595499958*magB2_yl[6]+0.5*magB2_yl[2]; 
  gkyl_mat_set(&A_byby_yl,3,1,temp_magB2_yl); 
 
  temp_magB2_yr = 0.447213595499958*magB2_yr[6]+0.5*magB2_yr[2]; 
  gkyl_mat_set(&A_byby_yr,3,1,temp_magB2_yr); 
 
  temp_magB2_zl = 0.447213595499958*magB2_zl[6]+0.5*magB2_zl[2]; 
  gkyl_mat_set(&A_bzbz_zl,3,1,temp_magB2_zl); 
 
  temp_magB2_zr = 0.447213595499958*magB2_zr[6]+0.5*magB2_zr[2]; 
  gkyl_mat_set(&A_bzbz_zr,3,1,temp_magB2_zr); 
 
  temp_magB2_xl = 0.447213595499958*magB2_xl[7]+0.5*magB2_xl[1]; 
  gkyl_mat_set(&A_bxbx_xl,3,2,temp_magB2_xl); 
 
  temp_magB2_xr = 0.447213595499958*magB2_xr[7]+0.5*magB2_xr[1]; 
  gkyl_mat_set(&A_bxbx_xr,3,2,temp_magB2_xr); 
 
  temp_magB2_yl = 0.447213595499958*magB2_yl[7]+0.5*magB2_yl[1]; 
  gkyl_mat_set(&A_byby_yl,3,2,temp_magB2_yl); 
 
  temp_magB2_yr = 0.447213595499958*magB2_yr[7]+0.5*magB2_yr[1]; 
  gkyl_mat_set(&A_byby_yr,3,2,temp_magB2_yr); 
 
  temp_magB2_zl = 0.447213595499958*magB2_zl[7]+0.5*magB2_zl[1]; 
  gkyl_mat_set(&A_bzbz_zl,3,2,temp_magB2_zl); 
 
  temp_magB2_zr = 0.447213595499958*magB2_zr[7]+0.5*magB2_zr[1]; 
  gkyl_mat_set(&A_bzbz_zr,3,2,temp_magB2_zr); 
 
  temp_magB2_xl = 0.4*magB2_xl[8]+0.4472135954999579*magB2_xl[5]+0.4472135954999579*magB2_xl[4]+0.5*magB2_xl[0]; 
  gkyl_mat_set(&A_bxbx_xl,3,3,temp_magB2_xl); 
 
  temp_magB2_xr = 0.4*magB2_xr[8]+0.4472135954999579*magB2_xr[5]+0.4472135954999579*magB2_xr[4]+0.5*magB2_xr[0]; 
  gkyl_mat_set(&A_bxbx_xr,3,3,temp_magB2_xr); 
 
  temp_magB2_yl = 0.4*magB2_yl[8]+0.4472135954999579*magB2_yl[5]+0.4472135954999579*magB2_yl[4]+0.5*magB2_yl[0]; 
  gkyl_mat_set(&A_byby_yl,3,3,temp_magB2_yl); 
 
  temp_magB2_yr = 0.4*magB2_yr[8]+0.4472135954999579*magB2_yr[5]+0.4472135954999579*magB2_yr[4]+0.5*magB2_yr[0]; 
  gkyl_mat_set(&A_byby_yr,3,3,temp_magB2_yr); 
 
  temp_magB2_zl = 0.4*magB2_zl[8]+0.4472135954999579*magB2_zl[5]+0.4472135954999579*magB2_zl[4]+0.5*magB2_zl[0]; 
  gkyl_mat_set(&A_bzbz_zl,3,3,temp_magB2_zl); 
 
  temp_magB2_zr = 0.4*magB2_zr[8]+0.4472135954999579*magB2_zr[5]+0.4472135954999579*magB2_zr[4]+0.5*magB2_zr[0]; 
  gkyl_mat_set(&A_bzbz_zr,3,3,temp_magB2_zr); 
 
  temp_magB2_xl = 0.4472135954999579*magB2_xl[3]; 
  gkyl_mat_set(&A_bxbx_xl,3,4,temp_magB2_xl); 
 
  temp_magB2_xr = 0.4472135954999579*magB2_xr[3]; 
  gkyl_mat_set(&A_bxbx_xr,3,4,temp_magB2_xr); 
 
  temp_magB2_yl = 0.4472135954999579*magB2_yl[3]; 
  gkyl_mat_set(&A_byby_yl,3,4,temp_magB2_yl); 
 
  temp_magB2_yr = 0.4472135954999579*magB2_yr[3]; 
  gkyl_mat_set(&A_byby_yr,3,4,temp_magB2_yr); 
 
  temp_magB2_zl = 0.4472135954999579*magB2_zl[3]; 
  gkyl_mat_set(&A_bzbz_zl,3,4,temp_magB2_zl); 
 
  temp_magB2_zr = 0.4472135954999579*magB2_zr[3]; 
  gkyl_mat_set(&A_bzbz_zr,3,4,temp_magB2_zr); 
 
  temp_magB2_xl = 0.4472135954999579*magB2_xl[3]; 
  gkyl_mat_set(&A_bxbx_xl,3,5,temp_magB2_xl); 
 
  temp_magB2_xr = 0.4472135954999579*magB2_xr[3]; 
  gkyl_mat_set(&A_bxbx_xr,3,5,temp_magB2_xr); 
 
  temp_magB2_yl = 0.4472135954999579*magB2_yl[3]; 
  gkyl_mat_set(&A_byby_yl,3,5,temp_magB2_yl); 
 
  temp_magB2_yr = 0.4472135954999579*magB2_yr[3]; 
  gkyl_mat_set(&A_byby_yr,3,5,temp_magB2_yr); 
 
  temp_magB2_zl = 0.4472135954999579*magB2_zl[3]; 
  gkyl_mat_set(&A_bzbz_zl,3,5,temp_magB2_zl); 
 
  temp_magB2_zr = 0.4472135954999579*magB2_zr[3]; 
  gkyl_mat_set(&A_bzbz_zr,3,5,temp_magB2_zr); 
 
  temp_magB2_xl = 0.4*magB2_xl[7]+0.447213595499958*magB2_xl[1]; 
  gkyl_mat_set(&A_bxbx_xl,3,6,temp_magB2_xl); 
 
  temp_magB2_xr = 0.4*magB2_xr[7]+0.447213595499958*magB2_xr[1]; 
  gkyl_mat_set(&A_bxbx_xr,3,6,temp_magB2_xr); 
 
  temp_magB2_yl = 0.4*magB2_yl[7]+0.447213595499958*magB2_yl[1]; 
  gkyl_mat_set(&A_byby_yl,3,6,temp_magB2_yl); 
 
  temp_magB2_yr = 0.4*magB2_yr[7]+0.447213595499958*magB2_yr[1]; 
  gkyl_mat_set(&A_byby_yr,3,6,temp_magB2_yr); 
 
  temp_magB2_zl = 0.4*magB2_zl[7]+0.447213595499958*magB2_zl[1]; 
  gkyl_mat_set(&A_bzbz_zl,3,6,temp_magB2_zl); 
 
  temp_magB2_zr = 0.4*magB2_zr[7]+0.447213595499958*magB2_zr[1]; 
  gkyl_mat_set(&A_bzbz_zr,3,6,temp_magB2_zr); 
 
  temp_magB2_xl = 0.4*magB2_xl[6]+0.447213595499958*magB2_xl[2]; 
  gkyl_mat_set(&A_bxbx_xl,3,7,temp_magB2_xl); 
 
  temp_magB2_xr = 0.4*magB2_xr[6]+0.447213595499958*magB2_xr[2]; 
  gkyl_mat_set(&A_bxbx_xr,3,7,temp_magB2_xr); 
 
  temp_magB2_yl = 0.4*magB2_yl[6]+0.447213595499958*magB2_yl[2]; 
  gkyl_mat_set(&A_byby_yl,3,7,temp_magB2_yl); 
 
  temp_magB2_yr = 0.4*magB2_yr[6]+0.447213595499958*magB2_yr[2]; 
  gkyl_mat_set(&A_byby_yr,3,7,temp_magB2_yr); 
 
  temp_magB2_zl = 0.4*magB2_zl[6]+0.447213595499958*magB2_zl[2]; 
  gkyl_mat_set(&A_bzbz_zl,3,7,temp_magB2_zl); 
 
  temp_magB2_zr = 0.4*magB2_zr[6]+0.447213595499958*magB2_zr[2]; 
  gkyl_mat_set(&A_bzbz_zr,3,7,temp_magB2_zr); 
 
  temp_magB2_xl = 0.4*magB2_xl[3]; 
  gkyl_mat_set(&A_bxbx_xl,3,8,temp_magB2_xl); 
 
  temp_magB2_xr = 0.4*magB2_xr[3]; 
  gkyl_mat_set(&A_bxbx_xr,3,8,temp_magB2_xr); 
 
  temp_magB2_yl = 0.4*magB2_yl[3]; 
  gkyl_mat_set(&A_byby_yl,3,8,temp_magB2_yl); 
 
  temp_magB2_yr = 0.4*magB2_yr[3]; 
  gkyl_mat_set(&A_byby_yr,3,8,temp_magB2_yr); 
 
  temp_magB2_zl = 0.4*magB2_zl[3]; 
  gkyl_mat_set(&A_bzbz_zl,3,8,temp_magB2_zl); 
 
  temp_magB2_zr = 0.4*magB2_zr[3]; 
  gkyl_mat_set(&A_bzbz_zr,3,8,temp_magB2_zr); 
 
  temp_magB2_xl = 0.5*magB2_xl[4]; 
  gkyl_mat_set(&A_bxbx_xl,4,0,temp_magB2_xl); 
 
  temp_magB2_xr = 0.5*magB2_xr[4]; 
  gkyl_mat_set(&A_bxbx_xr,4,0,temp_magB2_xr); 
 
  temp_magB2_yl = 0.5*magB2_yl[4]; 
  gkyl_mat_set(&A_byby_yl,4,0,temp_magB2_yl); 
 
  temp_magB2_yr = 0.5*magB2_yr[4]; 
  gkyl_mat_set(&A_byby_yr,4,0,temp_magB2_yr); 
 
  temp_magB2_zl = 0.5*magB2_zl[4]; 
  gkyl_mat_set(&A_bzbz_zl,4,0,temp_magB2_zl); 
 
  temp_magB2_zr = 0.5*magB2_zr[4]; 
  gkyl_mat_set(&A_bzbz_zr,4,0,temp_magB2_zr); 
 
  temp_magB2_xl = 0.4472135954999579*magB2_xl[1]; 
  gkyl_mat_set(&A_bxbx_xl,4,1,temp_magB2_xl); 
 
  temp_magB2_xr = 0.4472135954999579*magB2_xr[1]; 
  gkyl_mat_set(&A_bxbx_xr,4,1,temp_magB2_xr); 
 
  temp_magB2_yl = 0.4472135954999579*magB2_yl[1]; 
  gkyl_mat_set(&A_byby_yl,4,1,temp_magB2_yl); 
 
  temp_magB2_yr = 0.4472135954999579*magB2_yr[1]; 
  gkyl_mat_set(&A_byby_yr,4,1,temp_magB2_yr); 
 
  temp_magB2_zl = 0.4472135954999579*magB2_zl[1]; 
  gkyl_mat_set(&A_bzbz_zl,4,1,temp_magB2_zl); 
 
  temp_magB2_zr = 0.4472135954999579*magB2_zr[1]; 
  gkyl_mat_set(&A_bzbz_zr,4,1,temp_magB2_zr); 
 
  temp_magB2_xl = 0.5000000000000001*magB2_xl[6]; 
  gkyl_mat_set(&A_bxbx_xl,4,2,temp_magB2_xl); 
 
  temp_magB2_xr = 0.5000000000000001*magB2_xr[6]; 
  gkyl_mat_set(&A_bxbx_xr,4,2,temp_magB2_xr); 
 
  temp_magB2_yl = 0.5000000000000001*magB2_yl[6]; 
  gkyl_mat_set(&A_byby_yl,4,2,temp_magB2_yl); 
 
  temp_magB2_yr = 0.5000000000000001*magB2_yr[6]; 
  gkyl_mat_set(&A_byby_yr,4,2,temp_magB2_yr); 
 
  temp_magB2_zl = 0.5000000000000001*magB2_zl[6]; 
  gkyl_mat_set(&A_bzbz_zl,4,2,temp_magB2_zl); 
 
  temp_magB2_zr = 0.5000000000000001*magB2_zr[6]; 
  gkyl_mat_set(&A_bzbz_zr,4,2,temp_magB2_zr); 
 
  temp_magB2_xl = 0.4472135954999579*magB2_xl[3]; 
  gkyl_mat_set(&A_bxbx_xl,4,3,temp_magB2_xl); 
 
  temp_magB2_xr = 0.4472135954999579*magB2_xr[3]; 
  gkyl_mat_set(&A_bxbx_xr,4,3,temp_magB2_xr); 
 
  temp_magB2_yl = 0.4472135954999579*magB2_yl[3]; 
  gkyl_mat_set(&A_byby_yl,4,3,temp_magB2_yl); 
 
  temp_magB2_yr = 0.4472135954999579*magB2_yr[3]; 
  gkyl_mat_set(&A_byby_yr,4,3,temp_magB2_yr); 
 
  temp_magB2_zl = 0.4472135954999579*magB2_zl[3]; 
  gkyl_mat_set(&A_bzbz_zl,4,3,temp_magB2_zl); 
 
  temp_magB2_zr = 0.4472135954999579*magB2_zr[3]; 
  gkyl_mat_set(&A_bzbz_zr,4,3,temp_magB2_zr); 
 
  temp_magB2_xl = 0.31943828249997*magB2_xl[4]+0.5*magB2_xl[0]; 
  gkyl_mat_set(&A_bxbx_xl,4,4,temp_magB2_xl); 
 
  temp_magB2_xr = 0.31943828249997*magB2_xr[4]+0.5*magB2_xr[0]; 
  gkyl_mat_set(&A_bxbx_xr,4,4,temp_magB2_xr); 
 
  temp_magB2_yl = 0.31943828249997*magB2_yl[4]+0.5*magB2_yl[0]; 
  gkyl_mat_set(&A_byby_yl,4,4,temp_magB2_yl); 
 
  temp_magB2_yr = 0.31943828249997*magB2_yr[4]+0.5*magB2_yr[0]; 
  gkyl_mat_set(&A_byby_yr,4,4,temp_magB2_yr); 
 
  temp_magB2_zl = 0.31943828249997*magB2_zl[4]+0.5*magB2_zl[0]; 
  gkyl_mat_set(&A_bzbz_zl,4,4,temp_magB2_zl); 
 
  temp_magB2_zr = 0.31943828249997*magB2_zr[4]+0.5*magB2_zr[0]; 
  gkyl_mat_set(&A_bzbz_zr,4,4,temp_magB2_zr); 
 
  temp_magB2_xl = 0.5*magB2_xl[8]; 
  gkyl_mat_set(&A_bxbx_xl,4,5,temp_magB2_xl); 
 
  temp_magB2_xr = 0.5*magB2_xr[8]; 
  gkyl_mat_set(&A_bxbx_xr,4,5,temp_magB2_xr); 
 
  temp_magB2_yl = 0.5*magB2_yl[8]; 
  gkyl_mat_set(&A_byby_yl,4,5,temp_magB2_yl); 
 
  temp_magB2_yr = 0.5*magB2_yr[8]; 
  gkyl_mat_set(&A_byby_yr,4,5,temp_magB2_yr); 
 
  temp_magB2_zl = 0.5*magB2_zl[8]; 
  gkyl_mat_set(&A_bzbz_zl,4,5,temp_magB2_zl); 
 
  temp_magB2_zr = 0.5*magB2_zr[8]; 
  gkyl_mat_set(&A_bzbz_zr,4,5,temp_magB2_zr); 
 
  temp_magB2_xl = 0.31943828249997*magB2_xl[6]+0.5000000000000001*magB2_xl[2]; 
  gkyl_mat_set(&A_bxbx_xl,4,6,temp_magB2_xl); 
 
  temp_magB2_xr = 0.31943828249997*magB2_xr[6]+0.5000000000000001*magB2_xr[2]; 
  gkyl_mat_set(&A_bxbx_xr,4,6,temp_magB2_xr); 
 
  temp_magB2_yl = 0.31943828249997*magB2_yl[6]+0.5000000000000001*magB2_yl[2]; 
  gkyl_mat_set(&A_byby_yl,4,6,temp_magB2_yl); 
 
  temp_magB2_yr = 0.31943828249997*magB2_yr[6]+0.5000000000000001*magB2_yr[2]; 
  gkyl_mat_set(&A_byby_yr,4,6,temp_magB2_yr); 
 
  temp_magB2_zl = 0.31943828249997*magB2_zl[6]+0.5000000000000001*magB2_zl[2]; 
  gkyl_mat_set(&A_bzbz_zl,4,6,temp_magB2_zl); 
 
  temp_magB2_zr = 0.31943828249997*magB2_zr[6]+0.5000000000000001*magB2_zr[2]; 
  gkyl_mat_set(&A_bzbz_zr,4,6,temp_magB2_zr); 
 
  temp_magB2_xl = 0.4472135954999579*magB2_xl[7]; 
  gkyl_mat_set(&A_bxbx_xl,4,7,temp_magB2_xl); 
 
  temp_magB2_xr = 0.4472135954999579*magB2_xr[7]; 
  gkyl_mat_set(&A_bxbx_xr,4,7,temp_magB2_xr); 
 
  temp_magB2_yl = 0.4472135954999579*magB2_yl[7]; 
  gkyl_mat_set(&A_byby_yl,4,7,temp_magB2_yl); 
 
  temp_magB2_yr = 0.4472135954999579*magB2_yr[7]; 
  gkyl_mat_set(&A_byby_yr,4,7,temp_magB2_yr); 
 
  temp_magB2_zl = 0.4472135954999579*magB2_zl[7]; 
  gkyl_mat_set(&A_bzbz_zl,4,7,temp_magB2_zl); 
 
  temp_magB2_zr = 0.4472135954999579*magB2_zr[7]; 
  gkyl_mat_set(&A_bzbz_zr,4,7,temp_magB2_zr); 
 
  temp_magB2_xl = 0.31943828249997*magB2_xl[8]+0.5*magB2_xl[5]; 
  gkyl_mat_set(&A_bxbx_xl,4,8,temp_magB2_xl); 
 
  temp_magB2_xr = 0.31943828249997*magB2_xr[8]+0.5*magB2_xr[5]; 
  gkyl_mat_set(&A_bxbx_xr,4,8,temp_magB2_xr); 
 
  temp_magB2_yl = 0.31943828249997*magB2_yl[8]+0.5*magB2_yl[5]; 
  gkyl_mat_set(&A_byby_yl,4,8,temp_magB2_yl); 
 
  temp_magB2_yr = 0.31943828249997*magB2_yr[8]+0.5*magB2_yr[5]; 
  gkyl_mat_set(&A_byby_yr,4,8,temp_magB2_yr); 
 
  temp_magB2_zl = 0.31943828249997*magB2_zl[8]+0.5*magB2_zl[5]; 
  gkyl_mat_set(&A_bzbz_zl,4,8,temp_magB2_zl); 
 
  temp_magB2_zr = 0.31943828249997*magB2_zr[8]+0.5*magB2_zr[5]; 
  gkyl_mat_set(&A_bzbz_zr,4,8,temp_magB2_zr); 
 
  temp_magB2_xl = 0.5*magB2_xl[5]; 
  gkyl_mat_set(&A_bxbx_xl,5,0,temp_magB2_xl); 
 
  temp_magB2_xr = 0.5*magB2_xr[5]; 
  gkyl_mat_set(&A_bxbx_xr,5,0,temp_magB2_xr); 
 
  temp_magB2_yl = 0.5*magB2_yl[5]; 
  gkyl_mat_set(&A_byby_yl,5,0,temp_magB2_yl); 
 
  temp_magB2_yr = 0.5*magB2_yr[5]; 
  gkyl_mat_set(&A_byby_yr,5,0,temp_magB2_yr); 
 
  temp_magB2_zl = 0.5*magB2_zl[5]; 
  gkyl_mat_set(&A_bzbz_zl,5,0,temp_magB2_zl); 
 
  temp_magB2_zr = 0.5*magB2_zr[5]; 
  gkyl_mat_set(&A_bzbz_zr,5,0,temp_magB2_zr); 
 
  temp_magB2_xl = 0.5000000000000001*magB2_xl[7]; 
  gkyl_mat_set(&A_bxbx_xl,5,1,temp_magB2_xl); 
 
  temp_magB2_xr = 0.5000000000000001*magB2_xr[7]; 
  gkyl_mat_set(&A_bxbx_xr,5,1,temp_magB2_xr); 
 
  temp_magB2_yl = 0.5000000000000001*magB2_yl[7]; 
  gkyl_mat_set(&A_byby_yl,5,1,temp_magB2_yl); 
 
  temp_magB2_yr = 0.5000000000000001*magB2_yr[7]; 
  gkyl_mat_set(&A_byby_yr,5,1,temp_magB2_yr); 
 
  temp_magB2_zl = 0.5000000000000001*magB2_zl[7]; 
  gkyl_mat_set(&A_bzbz_zl,5,1,temp_magB2_zl); 
 
  temp_magB2_zr = 0.5000000000000001*magB2_zr[7]; 
  gkyl_mat_set(&A_bzbz_zr,5,1,temp_magB2_zr); 
 
  temp_magB2_xl = 0.4472135954999579*magB2_xl[2]; 
  gkyl_mat_set(&A_bxbx_xl,5,2,temp_magB2_xl); 
 
  temp_magB2_xr = 0.4472135954999579*magB2_xr[2]; 
  gkyl_mat_set(&A_bxbx_xr,5,2,temp_magB2_xr); 
 
  temp_magB2_yl = 0.4472135954999579*magB2_yl[2]; 
  gkyl_mat_set(&A_byby_yl,5,2,temp_magB2_yl); 
 
  temp_magB2_yr = 0.4472135954999579*magB2_yr[2]; 
  gkyl_mat_set(&A_byby_yr,5,2,temp_magB2_yr); 
 
  temp_magB2_zl = 0.4472135954999579*magB2_zl[2]; 
  gkyl_mat_set(&A_bzbz_zl,5,2,temp_magB2_zl); 
 
  temp_magB2_zr = 0.4472135954999579*magB2_zr[2]; 
  gkyl_mat_set(&A_bzbz_zr,5,2,temp_magB2_zr); 
 
  temp_magB2_xl = 0.4472135954999579*magB2_xl[3]; 
  gkyl_mat_set(&A_bxbx_xl,5,3,temp_magB2_xl); 
 
  temp_magB2_xr = 0.4472135954999579*magB2_xr[3]; 
  gkyl_mat_set(&A_bxbx_xr,5,3,temp_magB2_xr); 
 
  temp_magB2_yl = 0.4472135954999579*magB2_yl[3]; 
  gkyl_mat_set(&A_byby_yl,5,3,temp_magB2_yl); 
 
  temp_magB2_yr = 0.4472135954999579*magB2_yr[3]; 
  gkyl_mat_set(&A_byby_yr,5,3,temp_magB2_yr); 
 
  temp_magB2_zl = 0.4472135954999579*magB2_zl[3]; 
  gkyl_mat_set(&A_bzbz_zl,5,3,temp_magB2_zl); 
 
  temp_magB2_zr = 0.4472135954999579*magB2_zr[3]; 
  gkyl_mat_set(&A_bzbz_zr,5,3,temp_magB2_zr); 
 
  temp_magB2_xl = 0.5*magB2_xl[8]; 
  gkyl_mat_set(&A_bxbx_xl,5,4,temp_magB2_xl); 
 
  temp_magB2_xr = 0.5*magB2_xr[8]; 
  gkyl_mat_set(&A_bxbx_xr,5,4,temp_magB2_xr); 
 
  temp_magB2_yl = 0.5*magB2_yl[8]; 
  gkyl_mat_set(&A_byby_yl,5,4,temp_magB2_yl); 
 
  temp_magB2_yr = 0.5*magB2_yr[8]; 
  gkyl_mat_set(&A_byby_yr,5,4,temp_magB2_yr); 
 
  temp_magB2_zl = 0.5*magB2_zl[8]; 
  gkyl_mat_set(&A_bzbz_zl,5,4,temp_magB2_zl); 
 
  temp_magB2_zr = 0.5*magB2_zr[8]; 
  gkyl_mat_set(&A_bzbz_zr,5,4,temp_magB2_zr); 
 
  temp_magB2_xl = 0.31943828249997*magB2_xl[5]+0.5*magB2_xl[0]; 
  gkyl_mat_set(&A_bxbx_xl,5,5,temp_magB2_xl); 
 
  temp_magB2_xr = 0.31943828249997*magB2_xr[5]+0.5*magB2_xr[0]; 
  gkyl_mat_set(&A_bxbx_xr,5,5,temp_magB2_xr); 
 
  temp_magB2_yl = 0.31943828249997*magB2_yl[5]+0.5*magB2_yl[0]; 
  gkyl_mat_set(&A_byby_yl,5,5,temp_magB2_yl); 
 
  temp_magB2_yr = 0.31943828249997*magB2_yr[5]+0.5*magB2_yr[0]; 
  gkyl_mat_set(&A_byby_yr,5,5,temp_magB2_yr); 
 
  temp_magB2_zl = 0.31943828249997*magB2_zl[5]+0.5*magB2_zl[0]; 
  gkyl_mat_set(&A_bzbz_zl,5,5,temp_magB2_zl); 
 
  temp_magB2_zr = 0.31943828249997*magB2_zr[5]+0.5*magB2_zr[0]; 
  gkyl_mat_set(&A_bzbz_zr,5,5,temp_magB2_zr); 
 
  temp_magB2_xl = 0.4472135954999579*magB2_xl[6]; 
  gkyl_mat_set(&A_bxbx_xl,5,6,temp_magB2_xl); 
 
  temp_magB2_xr = 0.4472135954999579*magB2_xr[6]; 
  gkyl_mat_set(&A_bxbx_xr,5,6,temp_magB2_xr); 
 
  temp_magB2_yl = 0.4472135954999579*magB2_yl[6]; 
  gkyl_mat_set(&A_byby_yl,5,6,temp_magB2_yl); 
 
  temp_magB2_yr = 0.4472135954999579*magB2_yr[6]; 
  gkyl_mat_set(&A_byby_yr,5,6,temp_magB2_yr); 
 
  temp_magB2_zl = 0.4472135954999579*magB2_zl[6]; 
  gkyl_mat_set(&A_bzbz_zl,5,6,temp_magB2_zl); 
 
  temp_magB2_zr = 0.4472135954999579*magB2_zr[6]; 
  gkyl_mat_set(&A_bzbz_zr,5,6,temp_magB2_zr); 
 
  temp_magB2_xl = 0.31943828249997*magB2_xl[7]+0.5000000000000001*magB2_xl[1]; 
  gkyl_mat_set(&A_bxbx_xl,5,7,temp_magB2_xl); 
 
  temp_magB2_xr = 0.31943828249997*magB2_xr[7]+0.5000000000000001*magB2_xr[1]; 
  gkyl_mat_set(&A_bxbx_xr,5,7,temp_magB2_xr); 
 
  temp_magB2_yl = 0.31943828249997*magB2_yl[7]+0.5000000000000001*magB2_yl[1]; 
  gkyl_mat_set(&A_byby_yl,5,7,temp_magB2_yl); 
 
  temp_magB2_yr = 0.31943828249997*magB2_yr[7]+0.5000000000000001*magB2_yr[1]; 
  gkyl_mat_set(&A_byby_yr,5,7,temp_magB2_yr); 
 
  temp_magB2_zl = 0.31943828249997*magB2_zl[7]+0.5000000000000001*magB2_zl[1]; 
  gkyl_mat_set(&A_bzbz_zl,5,7,temp_magB2_zl); 
 
  temp_magB2_zr = 0.31943828249997*magB2_zr[7]+0.5000000000000001*magB2_zr[1]; 
  gkyl_mat_set(&A_bzbz_zr,5,7,temp_magB2_zr); 
 
  temp_magB2_xl = 0.31943828249997*magB2_xl[8]+0.5*magB2_xl[4]; 
  gkyl_mat_set(&A_bxbx_xl,5,8,temp_magB2_xl); 
 
  temp_magB2_xr = 0.31943828249997*magB2_xr[8]+0.5*magB2_xr[4]; 
  gkyl_mat_set(&A_bxbx_xr,5,8,temp_magB2_xr); 
 
  temp_magB2_yl = 0.31943828249997*magB2_yl[8]+0.5*magB2_yl[4]; 
  gkyl_mat_set(&A_byby_yl,5,8,temp_magB2_yl); 
 
  temp_magB2_yr = 0.31943828249997*magB2_yr[8]+0.5*magB2_yr[4]; 
  gkyl_mat_set(&A_byby_yr,5,8,temp_magB2_yr); 
 
  temp_magB2_zl = 0.31943828249997*magB2_zl[8]+0.5*magB2_zl[4]; 
  gkyl_mat_set(&A_bzbz_zl,5,8,temp_magB2_zl); 
 
  temp_magB2_zr = 0.31943828249997*magB2_zr[8]+0.5*magB2_zr[4]; 
  gkyl_mat_set(&A_bzbz_zr,5,8,temp_magB2_zr); 
 
  temp_magB2_xl = 0.5*magB2_xl[6]; 
  gkyl_mat_set(&A_bxbx_xl,6,0,temp_magB2_xl); 
 
  temp_magB2_xr = 0.5*magB2_xr[6]; 
  gkyl_mat_set(&A_bxbx_xr,6,0,temp_magB2_xr); 
 
  temp_magB2_yl = 0.5*magB2_yl[6]; 
  gkyl_mat_set(&A_byby_yl,6,0,temp_magB2_yl); 
 
  temp_magB2_yr = 0.5*magB2_yr[6]; 
  gkyl_mat_set(&A_byby_yr,6,0,temp_magB2_yr); 
 
  temp_magB2_zl = 0.5*magB2_zl[6]; 
  gkyl_mat_set(&A_bzbz_zl,6,0,temp_magB2_zl); 
 
  temp_magB2_zr = 0.5*magB2_zr[6]; 
  gkyl_mat_set(&A_bzbz_zr,6,0,temp_magB2_zr); 
 
  temp_magB2_xl = 0.447213595499958*magB2_xl[3]; 
  gkyl_mat_set(&A_bxbx_xl,6,1,temp_magB2_xl); 
 
  temp_magB2_xr = 0.447213595499958*magB2_xr[3]; 
  gkyl_mat_set(&A_bxbx_xr,6,1,temp_magB2_xr); 
 
  temp_magB2_yl = 0.447213595499958*magB2_yl[3]; 
  gkyl_mat_set(&A_byby_yl,6,1,temp_magB2_yl); 
 
  temp_magB2_yr = 0.447213595499958*magB2_yr[3]; 
  gkyl_mat_set(&A_byby_yr,6,1,temp_magB2_yr); 
 
  temp_magB2_zl = 0.447213595499958*magB2_zl[3]; 
  gkyl_mat_set(&A_bzbz_zl,6,1,temp_magB2_zl); 
 
  temp_magB2_zr = 0.447213595499958*magB2_zr[3]; 
  gkyl_mat_set(&A_bzbz_zr,6,1,temp_magB2_zr); 
 
  temp_magB2_xl = 0.447213595499958*magB2_xl[8]+0.5000000000000001*magB2_xl[4]; 
  gkyl_mat_set(&A_bxbx_xl,6,2,temp_magB2_xl); 
 
  temp_magB2_xr = 0.447213595499958*magB2_xr[8]+0.5000000000000001*magB2_xr[4]; 
  gkyl_mat_set(&A_bxbx_xr,6,2,temp_magB2_xr); 
 
  temp_magB2_yl = 0.447213595499958*magB2_yl[8]+0.5000000000000001*magB2_yl[4]; 
  gkyl_mat_set(&A_byby_yl,6,2,temp_magB2_yl); 
 
  temp_magB2_yr = 0.447213595499958*magB2_yr[8]+0.5000000000000001*magB2_yr[4]; 
  gkyl_mat_set(&A_byby_yr,6,2,temp_magB2_yr); 
 
  temp_magB2_zl = 0.447213595499958*magB2_zl[8]+0.5000000000000001*magB2_zl[4]; 
  gkyl_mat_set(&A_bzbz_zl,6,2,temp_magB2_zl); 
 
  temp_magB2_zr = 0.447213595499958*magB2_zr[8]+0.5000000000000001*magB2_zr[4]; 
  gkyl_mat_set(&A_bzbz_zr,6,2,temp_magB2_zr); 
 
  temp_magB2_xl = 0.4*magB2_xl[7]+0.447213595499958*magB2_xl[1]; 
  gkyl_mat_set(&A_bxbx_xl,6,3,temp_magB2_xl); 
 
  temp_magB2_xr = 0.4*magB2_xr[7]+0.447213595499958*magB2_xr[1]; 
  gkyl_mat_set(&A_bxbx_xr,6,3,temp_magB2_xr); 
 
  temp_magB2_yl = 0.4*magB2_yl[7]+0.447213595499958*magB2_yl[1]; 
  gkyl_mat_set(&A_byby_yl,6,3,temp_magB2_yl); 
 
  temp_magB2_yr = 0.4*magB2_yr[7]+0.447213595499958*magB2_yr[1]; 
  gkyl_mat_set(&A_byby_yr,6,3,temp_magB2_yr); 
 
  temp_magB2_zl = 0.4*magB2_zl[7]+0.447213595499958*magB2_zl[1]; 
  gkyl_mat_set(&A_bzbz_zl,6,3,temp_magB2_zl); 
 
  temp_magB2_zr = 0.4*magB2_zr[7]+0.447213595499958*magB2_zr[1]; 
  gkyl_mat_set(&A_bzbz_zr,6,3,temp_magB2_zr); 
 
  temp_magB2_xl = 0.31943828249997*magB2_xl[6]+0.5000000000000001*magB2_xl[2]; 
  gkyl_mat_set(&A_bxbx_xl,6,4,temp_magB2_xl); 
 
  temp_magB2_xr = 0.31943828249997*magB2_xr[6]+0.5000000000000001*magB2_xr[2]; 
  gkyl_mat_set(&A_bxbx_xr,6,4,temp_magB2_xr); 
 
  temp_magB2_yl = 0.31943828249997*magB2_yl[6]+0.5000000000000001*magB2_yl[2]; 
  gkyl_mat_set(&A_byby_yl,6,4,temp_magB2_yl); 
 
  temp_magB2_yr = 0.31943828249997*magB2_yr[6]+0.5000000000000001*magB2_yr[2]; 
  gkyl_mat_set(&A_byby_yr,6,4,temp_magB2_yr); 
 
  temp_magB2_zl = 0.31943828249997*magB2_zl[6]+0.5000000000000001*magB2_zl[2]; 
  gkyl_mat_set(&A_bzbz_zl,6,4,temp_magB2_zl); 
 
  temp_magB2_zr = 0.31943828249997*magB2_zr[6]+0.5000000000000001*magB2_zr[2]; 
  gkyl_mat_set(&A_bzbz_zr,6,4,temp_magB2_zr); 
 
  temp_magB2_xl = 0.4472135954999579*magB2_xl[6]; 
  gkyl_mat_set(&A_bxbx_xl,6,5,temp_magB2_xl); 
 
  temp_magB2_xr = 0.4472135954999579*magB2_xr[6]; 
  gkyl_mat_set(&A_bxbx_xr,6,5,temp_magB2_xr); 
 
  temp_magB2_yl = 0.4472135954999579*magB2_yl[6]; 
  gkyl_mat_set(&A_byby_yl,6,5,temp_magB2_yl); 
 
  temp_magB2_yr = 0.4472135954999579*magB2_yr[6]; 
  gkyl_mat_set(&A_byby_yr,6,5,temp_magB2_yr); 
 
  temp_magB2_zl = 0.4472135954999579*magB2_zl[6]; 
  gkyl_mat_set(&A_bzbz_zl,6,5,temp_magB2_zl); 
 
  temp_magB2_zr = 0.4472135954999579*magB2_zr[6]; 
  gkyl_mat_set(&A_bzbz_zr,6,5,temp_magB2_zr); 
 
  temp_magB2_xl = 0.2857142857142857*magB2_xl[8]+0.4472135954999579*magB2_xl[5]+0.31943828249997*magB2_xl[4]+0.5*magB2_xl[0]; 
  gkyl_mat_set(&A_bxbx_xl,6,6,temp_magB2_xl); 
 
  temp_magB2_xr = 0.2857142857142857*magB2_xr[8]+0.4472135954999579*magB2_xr[5]+0.31943828249997*magB2_xr[4]+0.5*magB2_xr[0]; 
  gkyl_mat_set(&A_bxbx_xr,6,6,temp_magB2_xr); 
 
  temp_magB2_yl = 0.2857142857142857*magB2_yl[8]+0.4472135954999579*magB2_yl[5]+0.31943828249997*magB2_yl[4]+0.5*magB2_yl[0]; 
  gkyl_mat_set(&A_byby_yl,6,6,temp_magB2_yl); 
 
  temp_magB2_yr = 0.2857142857142857*magB2_yr[8]+0.4472135954999579*magB2_yr[5]+0.31943828249997*magB2_yr[4]+0.5*magB2_yr[0]; 
  gkyl_mat_set(&A_byby_yr,6,6,temp_magB2_yr); 
 
  temp_magB2_zl = 0.2857142857142857*magB2_zl[8]+0.4472135954999579*magB2_zl[5]+0.31943828249997*magB2_zl[4]+0.5*magB2_zl[0]; 
  gkyl_mat_set(&A_bzbz_zl,6,6,temp_magB2_zl); 
 
  temp_magB2_zr = 0.2857142857142857*magB2_zr[8]+0.4472135954999579*magB2_zr[5]+0.31943828249997*magB2_zr[4]+0.5*magB2_zr[0]; 
  gkyl_mat_set(&A_bzbz_zr,6,6,temp_magB2_zr); 
 
  temp_magB2_xl = 0.4*magB2_xl[3]; 
  gkyl_mat_set(&A_bxbx_xl,6,7,temp_magB2_xl); 
 
  temp_magB2_xr = 0.4*magB2_xr[3]; 
  gkyl_mat_set(&A_bxbx_xr,6,7,temp_magB2_xr); 
 
  temp_magB2_yl = 0.4*magB2_yl[3]; 
  gkyl_mat_set(&A_byby_yl,6,7,temp_magB2_yl); 
 
  temp_magB2_yr = 0.4*magB2_yr[3]; 
  gkyl_mat_set(&A_byby_yr,6,7,temp_magB2_yr); 
 
  temp_magB2_zl = 0.4*magB2_zl[3]; 
  gkyl_mat_set(&A_bzbz_zl,6,7,temp_magB2_zl); 
 
  temp_magB2_zr = 0.4*magB2_zr[3]; 
  gkyl_mat_set(&A_bzbz_zr,6,7,temp_magB2_zr); 
 
  temp_magB2_xl = 0.2857142857142857*magB2_xl[6]+0.447213595499958*magB2_xl[2]; 
  gkyl_mat_set(&A_bxbx_xl,6,8,temp_magB2_xl); 
 
  temp_magB2_xr = 0.2857142857142857*magB2_xr[6]+0.447213595499958*magB2_xr[2]; 
  gkyl_mat_set(&A_bxbx_xr,6,8,temp_magB2_xr); 
 
  temp_magB2_yl = 0.2857142857142857*magB2_yl[6]+0.447213595499958*magB2_yl[2]; 
  gkyl_mat_set(&A_byby_yl,6,8,temp_magB2_yl); 
 
  temp_magB2_yr = 0.2857142857142857*magB2_yr[6]+0.447213595499958*magB2_yr[2]; 
  gkyl_mat_set(&A_byby_yr,6,8,temp_magB2_yr); 
 
  temp_magB2_zl = 0.2857142857142857*magB2_zl[6]+0.447213595499958*magB2_zl[2]; 
  gkyl_mat_set(&A_bzbz_zl,6,8,temp_magB2_zl); 
 
  temp_magB2_zr = 0.2857142857142857*magB2_zr[6]+0.447213595499958*magB2_zr[2]; 
  gkyl_mat_set(&A_bzbz_zr,6,8,temp_magB2_zr); 
 
  temp_magB2_xl = 0.5*magB2_xl[7]; 
  gkyl_mat_set(&A_bxbx_xl,7,0,temp_magB2_xl); 
 
  temp_magB2_xr = 0.5*magB2_xr[7]; 
  gkyl_mat_set(&A_bxbx_xr,7,0,temp_magB2_xr); 
 
  temp_magB2_yl = 0.5*magB2_yl[7]; 
  gkyl_mat_set(&A_byby_yl,7,0,temp_magB2_yl); 
 
  temp_magB2_yr = 0.5*magB2_yr[7]; 
  gkyl_mat_set(&A_byby_yr,7,0,temp_magB2_yr); 
 
  temp_magB2_zl = 0.5*magB2_zl[7]; 
  gkyl_mat_set(&A_bzbz_zl,7,0,temp_magB2_zl); 
 
  temp_magB2_zr = 0.5*magB2_zr[7]; 
  gkyl_mat_set(&A_bzbz_zr,7,0,temp_magB2_zr); 
 
  temp_magB2_xl = 0.447213595499958*magB2_xl[8]+0.5000000000000001*magB2_xl[5]; 
  gkyl_mat_set(&A_bxbx_xl,7,1,temp_magB2_xl); 
 
  temp_magB2_xr = 0.447213595499958*magB2_xr[8]+0.5000000000000001*magB2_xr[5]; 
  gkyl_mat_set(&A_bxbx_xr,7,1,temp_magB2_xr); 
 
  temp_magB2_yl = 0.447213595499958*magB2_yl[8]+0.5000000000000001*magB2_yl[5]; 
  gkyl_mat_set(&A_byby_yl,7,1,temp_magB2_yl); 
 
  temp_magB2_yr = 0.447213595499958*magB2_yr[8]+0.5000000000000001*magB2_yr[5]; 
  gkyl_mat_set(&A_byby_yr,7,1,temp_magB2_yr); 
 
  temp_magB2_zl = 0.447213595499958*magB2_zl[8]+0.5000000000000001*magB2_zl[5]; 
  gkyl_mat_set(&A_bzbz_zl,7,1,temp_magB2_zl); 
 
  temp_magB2_zr = 0.447213595499958*magB2_zr[8]+0.5000000000000001*magB2_zr[5]; 
  gkyl_mat_set(&A_bzbz_zr,7,1,temp_magB2_zr); 
 
  temp_magB2_xl = 0.447213595499958*magB2_xl[3]; 
  gkyl_mat_set(&A_bxbx_xl,7,2,temp_magB2_xl); 
 
  temp_magB2_xr = 0.447213595499958*magB2_xr[3]; 
  gkyl_mat_set(&A_bxbx_xr,7,2,temp_magB2_xr); 
 
  temp_magB2_yl = 0.447213595499958*magB2_yl[3]; 
  gkyl_mat_set(&A_byby_yl,7,2,temp_magB2_yl); 
 
  temp_magB2_yr = 0.447213595499958*magB2_yr[3]; 
  gkyl_mat_set(&A_byby_yr,7,2,temp_magB2_yr); 
 
  temp_magB2_zl = 0.447213595499958*magB2_zl[3]; 
  gkyl_mat_set(&A_bzbz_zl,7,2,temp_magB2_zl); 
 
  temp_magB2_zr = 0.447213595499958*magB2_zr[3]; 
  gkyl_mat_set(&A_bzbz_zr,7,2,temp_magB2_zr); 
 
  temp_magB2_xl = 0.4*magB2_xl[6]+0.447213595499958*magB2_xl[2]; 
  gkyl_mat_set(&A_bxbx_xl,7,3,temp_magB2_xl); 
 
  temp_magB2_xr = 0.4*magB2_xr[6]+0.447213595499958*magB2_xr[2]; 
  gkyl_mat_set(&A_bxbx_xr,7,3,temp_magB2_xr); 
 
  temp_magB2_yl = 0.4*magB2_yl[6]+0.447213595499958*magB2_yl[2]; 
  gkyl_mat_set(&A_byby_yl,7,3,temp_magB2_yl); 
 
  temp_magB2_yr = 0.4*magB2_yr[6]+0.447213595499958*magB2_yr[2]; 
  gkyl_mat_set(&A_byby_yr,7,3,temp_magB2_yr); 
 
  temp_magB2_zl = 0.4*magB2_zl[6]+0.447213595499958*magB2_zl[2]; 
  gkyl_mat_set(&A_bzbz_zl,7,3,temp_magB2_zl); 
 
  temp_magB2_zr = 0.4*magB2_zr[6]+0.447213595499958*magB2_zr[2]; 
  gkyl_mat_set(&A_bzbz_zr,7,3,temp_magB2_zr); 
 
  temp_magB2_xl = 0.4472135954999579*magB2_xl[7]; 
  gkyl_mat_set(&A_bxbx_xl,7,4,temp_magB2_xl); 
 
  temp_magB2_xr = 0.4472135954999579*magB2_xr[7]; 
  gkyl_mat_set(&A_bxbx_xr,7,4,temp_magB2_xr); 
 
  temp_magB2_yl = 0.4472135954999579*magB2_yl[7]; 
  gkyl_mat_set(&A_byby_yl,7,4,temp_magB2_yl); 
 
  temp_magB2_yr = 0.4472135954999579*magB2_yr[7]; 
  gkyl_mat_set(&A_byby_yr,7,4,temp_magB2_yr); 
 
  temp_magB2_zl = 0.4472135954999579*magB2_zl[7]; 
  gkyl_mat_set(&A_bzbz_zl,7,4,temp_magB2_zl); 
 
  temp_magB2_zr = 0.4472135954999579*magB2_zr[7]; 
  gkyl_mat_set(&A_bzbz_zr,7,4,temp_magB2_zr); 
 
  temp_magB2_xl = 0.31943828249997*magB2_xl[7]+0.5000000000000001*magB2_xl[1]; 
  gkyl_mat_set(&A_bxbx_xl,7,5,temp_magB2_xl); 
 
  temp_magB2_xr = 0.31943828249997*magB2_xr[7]+0.5000000000000001*magB2_xr[1]; 
  gkyl_mat_set(&A_bxbx_xr,7,5,temp_magB2_xr); 
 
  temp_magB2_yl = 0.31943828249997*magB2_yl[7]+0.5000000000000001*magB2_yl[1]; 
  gkyl_mat_set(&A_byby_yl,7,5,temp_magB2_yl); 
 
  temp_magB2_yr = 0.31943828249997*magB2_yr[7]+0.5000000000000001*magB2_yr[1]; 
  gkyl_mat_set(&A_byby_yr,7,5,temp_magB2_yr); 
 
  temp_magB2_zl = 0.31943828249997*magB2_zl[7]+0.5000000000000001*magB2_zl[1]; 
  gkyl_mat_set(&A_bzbz_zl,7,5,temp_magB2_zl); 
 
  temp_magB2_zr = 0.31943828249997*magB2_zr[7]+0.5000000000000001*magB2_zr[1]; 
  gkyl_mat_set(&A_bzbz_zr,7,5,temp_magB2_zr); 
 
  temp_magB2_xl = 0.4*magB2_xl[3]; 
  gkyl_mat_set(&A_bxbx_xl,7,6,temp_magB2_xl); 
 
  temp_magB2_xr = 0.4*magB2_xr[3]; 
  gkyl_mat_set(&A_bxbx_xr,7,6,temp_magB2_xr); 
 
  temp_magB2_yl = 0.4*magB2_yl[3]; 
  gkyl_mat_set(&A_byby_yl,7,6,temp_magB2_yl); 
 
  temp_magB2_yr = 0.4*magB2_yr[3]; 
  gkyl_mat_set(&A_byby_yr,7,6,temp_magB2_yr); 
 
  temp_magB2_zl = 0.4*magB2_zl[3]; 
  gkyl_mat_set(&A_bzbz_zl,7,6,temp_magB2_zl); 
 
  temp_magB2_zr = 0.4*magB2_zr[3]; 
  gkyl_mat_set(&A_bzbz_zr,7,6,temp_magB2_zr); 
 
  temp_magB2_xl = 0.2857142857142857*magB2_xl[8]+0.31943828249997*magB2_xl[5]+0.4472135954999579*magB2_xl[4]+0.5*magB2_xl[0]; 
  gkyl_mat_set(&A_bxbx_xl,7,7,temp_magB2_xl); 
 
  temp_magB2_xr = 0.2857142857142857*magB2_xr[8]+0.31943828249997*magB2_xr[5]+0.4472135954999579*magB2_xr[4]+0.5*magB2_xr[0]; 
  gkyl_mat_set(&A_bxbx_xr,7,7,temp_magB2_xr); 
 
  temp_magB2_yl = 0.2857142857142857*magB2_yl[8]+0.31943828249997*magB2_yl[5]+0.4472135954999579*magB2_yl[4]+0.5*magB2_yl[0]; 
  gkyl_mat_set(&A_byby_yl,7,7,temp_magB2_yl); 
 
  temp_magB2_yr = 0.2857142857142857*magB2_yr[8]+0.31943828249997*magB2_yr[5]+0.4472135954999579*magB2_yr[4]+0.5*magB2_yr[0]; 
  gkyl_mat_set(&A_byby_yr,7,7,temp_magB2_yr); 
 
  temp_magB2_zl = 0.2857142857142857*magB2_zl[8]+0.31943828249997*magB2_zl[5]+0.4472135954999579*magB2_zl[4]+0.5*magB2_zl[0]; 
  gkyl_mat_set(&A_bzbz_zl,7,7,temp_magB2_zl); 
 
  temp_magB2_zr = 0.2857142857142857*magB2_zr[8]+0.31943828249997*magB2_zr[5]+0.4472135954999579*magB2_zr[4]+0.5*magB2_zr[0]; 
  gkyl_mat_set(&A_bzbz_zr,7,7,temp_magB2_zr); 
 
  temp_magB2_xl = 0.2857142857142857*magB2_xl[7]+0.447213595499958*magB2_xl[1]; 
  gkyl_mat_set(&A_bxbx_xl,7,8,temp_magB2_xl); 
 
  temp_magB2_xr = 0.2857142857142857*magB2_xr[7]+0.447213595499958*magB2_xr[1]; 
  gkyl_mat_set(&A_bxbx_xr,7,8,temp_magB2_xr); 
 
  temp_magB2_yl = 0.2857142857142857*magB2_yl[7]+0.447213595499958*magB2_yl[1]; 
  gkyl_mat_set(&A_byby_yl,7,8,temp_magB2_yl); 
 
  temp_magB2_yr = 0.2857142857142857*magB2_yr[7]+0.447213595499958*magB2_yr[1]; 
  gkyl_mat_set(&A_byby_yr,7,8,temp_magB2_yr); 
 
  temp_magB2_zl = 0.2857142857142857*magB2_zl[7]+0.447213595499958*magB2_zl[1]; 
  gkyl_mat_set(&A_bzbz_zl,7,8,temp_magB2_zl); 
 
  temp_magB2_zr = 0.2857142857142857*magB2_zr[7]+0.447213595499958*magB2_zr[1]; 
  gkyl_mat_set(&A_bzbz_zr,7,8,temp_magB2_zr); 
 
  temp_magB2_xl = 0.5*magB2_xl[8]; 
  gkyl_mat_set(&A_bxbx_xl,8,0,temp_magB2_xl); 
 
  temp_magB2_xr = 0.5*magB2_xr[8]; 
  gkyl_mat_set(&A_bxbx_xr,8,0,temp_magB2_xr); 
 
  temp_magB2_yl = 0.5*magB2_yl[8]; 
  gkyl_mat_set(&A_byby_yl,8,0,temp_magB2_yl); 
 
  temp_magB2_yr = 0.5*magB2_yr[8]; 
  gkyl_mat_set(&A_byby_yr,8,0,temp_magB2_yr); 
 
  temp_magB2_zl = 0.5*magB2_zl[8]; 
  gkyl_mat_set(&A_bzbz_zl,8,0,temp_magB2_zl); 
 
  temp_magB2_zr = 0.5*magB2_zr[8]; 
  gkyl_mat_set(&A_bzbz_zr,8,0,temp_magB2_zr); 
 
  temp_magB2_xl = 0.447213595499958*magB2_xl[7]; 
  gkyl_mat_set(&A_bxbx_xl,8,1,temp_magB2_xl); 
 
  temp_magB2_xr = 0.447213595499958*magB2_xr[7]; 
  gkyl_mat_set(&A_bxbx_xr,8,1,temp_magB2_xr); 
 
  temp_magB2_yl = 0.447213595499958*magB2_yl[7]; 
  gkyl_mat_set(&A_byby_yl,8,1,temp_magB2_yl); 
 
  temp_magB2_yr = 0.447213595499958*magB2_yr[7]; 
  gkyl_mat_set(&A_byby_yr,8,1,temp_magB2_yr); 
 
  temp_magB2_zl = 0.447213595499958*magB2_zl[7]; 
  gkyl_mat_set(&A_bzbz_zl,8,1,temp_magB2_zl); 
 
  temp_magB2_zr = 0.447213595499958*magB2_zr[7]; 
  gkyl_mat_set(&A_bzbz_zr,8,1,temp_magB2_zr); 
 
  temp_magB2_xl = 0.447213595499958*magB2_xl[6]; 
  gkyl_mat_set(&A_bxbx_xl,8,2,temp_magB2_xl); 
 
  temp_magB2_xr = 0.447213595499958*magB2_xr[6]; 
  gkyl_mat_set(&A_bxbx_xr,8,2,temp_magB2_xr); 
 
  temp_magB2_yl = 0.447213595499958*magB2_yl[6]; 
  gkyl_mat_set(&A_byby_yl,8,2,temp_magB2_yl); 
 
  temp_magB2_yr = 0.447213595499958*magB2_yr[6]; 
  gkyl_mat_set(&A_byby_yr,8,2,temp_magB2_yr); 
 
  temp_magB2_zl = 0.447213595499958*magB2_zl[6]; 
  gkyl_mat_set(&A_bzbz_zl,8,2,temp_magB2_zl); 
 
  temp_magB2_zr = 0.447213595499958*magB2_zr[6]; 
  gkyl_mat_set(&A_bzbz_zr,8,2,temp_magB2_zr); 
 
  temp_magB2_xl = 0.4*magB2_xl[3]; 
  gkyl_mat_set(&A_bxbx_xl,8,3,temp_magB2_xl); 
 
  temp_magB2_xr = 0.4*magB2_xr[3]; 
  gkyl_mat_set(&A_bxbx_xr,8,3,temp_magB2_xr); 
 
  temp_magB2_yl = 0.4*magB2_yl[3]; 
  gkyl_mat_set(&A_byby_yl,8,3,temp_magB2_yl); 
 
  temp_magB2_yr = 0.4*magB2_yr[3]; 
  gkyl_mat_set(&A_byby_yr,8,3,temp_magB2_yr); 
 
  temp_magB2_zl = 0.4*magB2_zl[3]; 
  gkyl_mat_set(&A_bzbz_zl,8,3,temp_magB2_zl); 
 
  temp_magB2_zr = 0.4*magB2_zr[3]; 
  gkyl_mat_set(&A_bzbz_zr,8,3,temp_magB2_zr); 
 
  temp_magB2_xl = 0.31943828249997*magB2_xl[8]+0.5*magB2_xl[5]; 
  gkyl_mat_set(&A_bxbx_xl,8,4,temp_magB2_xl); 
 
  temp_magB2_xr = 0.31943828249997*magB2_xr[8]+0.5*magB2_xr[5]; 
  gkyl_mat_set(&A_bxbx_xr,8,4,temp_magB2_xr); 
 
  temp_magB2_yl = 0.31943828249997*magB2_yl[8]+0.5*magB2_yl[5]; 
  gkyl_mat_set(&A_byby_yl,8,4,temp_magB2_yl); 
 
  temp_magB2_yr = 0.31943828249997*magB2_yr[8]+0.5*magB2_yr[5]; 
  gkyl_mat_set(&A_byby_yr,8,4,temp_magB2_yr); 
 
  temp_magB2_zl = 0.31943828249997*magB2_zl[8]+0.5*magB2_zl[5]; 
  gkyl_mat_set(&A_bzbz_zl,8,4,temp_magB2_zl); 
 
  temp_magB2_zr = 0.31943828249997*magB2_zr[8]+0.5*magB2_zr[5]; 
  gkyl_mat_set(&A_bzbz_zr,8,4,temp_magB2_zr); 
 
  temp_magB2_xl = 0.31943828249997*magB2_xl[8]+0.5*magB2_xl[4]; 
  gkyl_mat_set(&A_bxbx_xl,8,5,temp_magB2_xl); 
 
  temp_magB2_xr = 0.31943828249997*magB2_xr[8]+0.5*magB2_xr[4]; 
  gkyl_mat_set(&A_bxbx_xr,8,5,temp_magB2_xr); 
 
  temp_magB2_yl = 0.31943828249997*magB2_yl[8]+0.5*magB2_yl[4]; 
  gkyl_mat_set(&A_byby_yl,8,5,temp_magB2_yl); 
 
  temp_magB2_yr = 0.31943828249997*magB2_yr[8]+0.5*magB2_yr[4]; 
  gkyl_mat_set(&A_byby_yr,8,5,temp_magB2_yr); 
 
  temp_magB2_zl = 0.31943828249997*magB2_zl[8]+0.5*magB2_zl[4]; 
  gkyl_mat_set(&A_bzbz_zl,8,5,temp_magB2_zl); 
 
  temp_magB2_zr = 0.31943828249997*magB2_zr[8]+0.5*magB2_zr[4]; 
  gkyl_mat_set(&A_bzbz_zr,8,5,temp_magB2_zr); 
 
  temp_magB2_xl = 0.2857142857142857*magB2_xl[6]+0.447213595499958*magB2_xl[2]; 
  gkyl_mat_set(&A_bxbx_xl,8,6,temp_magB2_xl); 
 
  temp_magB2_xr = 0.2857142857142857*magB2_xr[6]+0.447213595499958*magB2_xr[2]; 
  gkyl_mat_set(&A_bxbx_xr,8,6,temp_magB2_xr); 
 
  temp_magB2_yl = 0.2857142857142857*magB2_yl[6]+0.447213595499958*magB2_yl[2]; 
  gkyl_mat_set(&A_byby_yl,8,6,temp_magB2_yl); 
 
  temp_magB2_yr = 0.2857142857142857*magB2_yr[6]+0.447213595499958*magB2_yr[2]; 
  gkyl_mat_set(&A_byby_yr,8,6,temp_magB2_yr); 
 
  temp_magB2_zl = 0.2857142857142857*magB2_zl[6]+0.447213595499958*magB2_zl[2]; 
  gkyl_mat_set(&A_bzbz_zl,8,6,temp_magB2_zl); 
 
  temp_magB2_zr = 0.2857142857142857*magB2_zr[6]+0.447213595499958*magB2_zr[2]; 
  gkyl_mat_set(&A_bzbz_zr,8,6,temp_magB2_zr); 
 
  temp_magB2_xl = 0.2857142857142857*magB2_xl[7]+0.447213595499958*magB2_xl[1]; 
  gkyl_mat_set(&A_bxbx_xl,8,7,temp_magB2_xl); 
 
  temp_magB2_xr = 0.2857142857142857*magB2_xr[7]+0.447213595499958*magB2_xr[1]; 
  gkyl_mat_set(&A_bxbx_xr,8,7,temp_magB2_xr); 
 
  temp_magB2_yl = 0.2857142857142857*magB2_yl[7]+0.447213595499958*magB2_yl[1]; 
  gkyl_mat_set(&A_byby_yl,8,7,temp_magB2_yl); 
 
  temp_magB2_yr = 0.2857142857142857*magB2_yr[7]+0.447213595499958*magB2_yr[1]; 
  gkyl_mat_set(&A_byby_yr,8,7,temp_magB2_yr); 
 
  temp_magB2_zl = 0.2857142857142857*magB2_zl[7]+0.447213595499958*magB2_zl[1]; 
  gkyl_mat_set(&A_bzbz_zl,8,7,temp_magB2_zl); 
 
  temp_magB2_zr = 0.2857142857142857*magB2_zr[7]+0.447213595499958*magB2_zr[1]; 
  gkyl_mat_set(&A_bzbz_zr,8,7,temp_magB2_zr); 
 
  temp_magB2_xl = 0.2040816326530612*magB2_xl[8]+0.31943828249997*magB2_xl[5]+0.31943828249997*magB2_xl[4]+0.5*magB2_xl[0]; 
  gkyl_mat_set(&A_bxbx_xl,8,8,temp_magB2_xl); 
 
  temp_magB2_xr = 0.2040816326530612*magB2_xr[8]+0.31943828249997*magB2_xr[5]+0.31943828249997*magB2_xr[4]+0.5*magB2_xr[0]; 
  gkyl_mat_set(&A_bxbx_xr,8,8,temp_magB2_xr); 
 
  temp_magB2_yl = 0.2040816326530612*magB2_yl[8]+0.31943828249997*magB2_yl[5]+0.31943828249997*magB2_yl[4]+0.5*magB2_yl[0]; 
  gkyl_mat_set(&A_byby_yl,8,8,temp_magB2_yl); 
 
  temp_magB2_yr = 0.2040816326530612*magB2_yr[8]+0.31943828249997*magB2_yr[5]+0.31943828249997*magB2_yr[4]+0.5*magB2_yr[0]; 
  gkyl_mat_set(&A_byby_yr,8,8,temp_magB2_yr); 
 
  temp_magB2_zl = 0.2040816326530612*magB2_zl[8]+0.31943828249997*magB2_zl[5]+0.31943828249997*magB2_zl[4]+0.5*magB2_zl[0]; 
  gkyl_mat_set(&A_bzbz_zl,8,8,temp_magB2_zl); 
 
  temp_magB2_zr = 0.2040816326530612*magB2_zr[8]+0.31943828249997*magB2_zr[5]+0.31943828249997*magB2_zr[4]+0.5*magB2_zr[0]; 
  gkyl_mat_set(&A_bzbz_zr,8,8,temp_magB2_zr); 
 
} 
