#include <gkyl_calc_metric.h>
#include <gkyl_rect_grid.h>
#include <gkyl_calc_metric_kernels.h>
#include <assert.h>
#include <math.h>

typedef void (*metric_kernel)(const double **xyz, double *gFld);

typedef struct { metric_kernel kernels[27]; } metric_kernel_list;  // For use in kernel tables.

GKYL_CU_DH
static const metric_kernel ser_metric_kernel_list[] = {
gij_inx_iny_inz_3x_Ser_p1, gij_lox_iny_inz_3x_Ser_p1, gij_upx_iny_inz_3x_Ser_p1, gij_inx_loy_inz_3x_Ser_p1, gij_inx_upy_inz_3x_Ser_p1, gij_lox_loy_inz_3x_Ser_p1, gij_lox_upy_inz_3x_Ser_p1, gij_upx_loy_inz_3x_Ser_p1, gij_upx_upy_inz_3x_Ser_p1, gij_inx_iny_loz_3x_Ser_p1, gij_inx_iny_upz_3x_Ser_p1, gij_lox_iny_loz_3x_Ser_p1, gij_lox_iny_upz_3x_Ser_p1, gij_upx_iny_loz_3x_Ser_p1, gij_upx_iny_upz_3x_Ser_p1, gij_inx_loy_loz_3x_Ser_p1, gij_inx_loy_upz_3x_Ser_p1, gij_inx_upy_loz_3x_Ser_p1, gij_inx_upy_upz_3x_Ser_p1, gij_lox_loy_loz_3x_Ser_p1, gij_lox_loy_upz_3x_Ser_p1, gij_lox_upy_loz_3x_Ser_p1, gij_lox_upy_upz_3x_Ser_p1, gij_upx_loy_loz_3x_Ser_p1, gij_upx_loy_upz_3x_Ser_p1, gij_upx_upy_loz_3x_Ser_p1, gij_upx_upy_upz_3x_Ser_p1
};



struct gkyl_calc_metric {
  unsigned cdim; // Configuration-space dimension.
  unsigned cnum_basis; // Number of conf-space basis functions.
  unsigned poly_order; // Polynomial order of the basis.
  struct gkyl_rect_grid* grid;
  bool use_gpu;
  metric_kernel kernel;
  int *num_cells;
};



GKYL_CU_DH
static inline int idx_to_inloup_ker(const int dim, const int *num_cells, const int *idx) {
  // Return the index of the kernel (in the array of kernels) needed given the grid index.
  // This function is for kernels that differentiate between lower, interior
  // and upper cells.
  int iout = 0;
  for (int d=0; d<dim; d++) {
    if (idx[d] == 1) {
      iout = 2*iout+(int)(pow(3,d)+0.5);
    } else if (idx[d] == num_cells[d]) {
      iout = 2*iout+(int)(pow(3,d)+0.5)+1;
    }
  }
  return iout;
}


GKYL_CU_DH
//static metric_kernel_list
//metric_choose_kernel(int dim, int basis_type, int poly_order, int lidx)
static metric_kernel
metric_choose_kernel(int lidx)
{
      return ser_metric_kernel_list[lidx];
  //switch (basis_type) {
  //  case GKYL_BASIS_MODAL_SERENDIPITY:
  //    //return ser_metric_kernel_list[dim].kernels[poly_order];
  //    return ser_metric_kernel_list[lidx];
  //  default:
  //    assert(false);
  //    break;
  //}
}


