#include <gkyl_calc_metric.h>
#include <gkyl_rect_grid.h>
#include <gkyl_calc_metric_kernels.h>
#include <assert.h>
#include <math.h>

typedef void (*metric_kernel)(const double **xyz, double *gFld);

typedef struct { metric_kernel kernels[27]; } metric_kernel_list;  // For use in kernel tables.

GKYL_CU_DH
static const metric_kernel ser_metric_kernel_list[] = {
  //gij_3x_Ser_p1,
  //gij_lo1_3x_Ser_p1, gij_up1_3x_Ser_p1,
  //gij_lo2_3x_Ser_p1, gij_up2_3x_Ser_p1,
  //gij_lo3_3x_Ser_p1, gij_up3_3x_Ser_p1,

  ////7
  //gij_lo1_lo2_3x_Ser_p1, gij_lo1_up2_3x_Ser_p1, gij_up1_lo2_3x_Ser_p1, gij_up1_up2_3x_Ser_p1, 
  //gij_lo1_lo3_3x_Ser_p1, gij_lo1_up3_3x_Ser_p1, gij_up1_lo3_3x_Ser_p1, gij_up1_up3_3x_Ser_p1, 
  //gij_lo2_lo3_3x_Ser_p1, gij_lo2_up3_3x_Ser_p1, gij_up2_lo3_3x_Ser_p1, gij_up2_up3_3x_Ser_p1, 

  ////19
  //gij_lo1_lo2_lo3_3x_Ser_p1, gij_lo1_lo2_up3_3x_Ser_p1, gij_lo1_up2_lo3_3x_Ser_p1, gij_lo1_up2_up3_3x_Ser_p1, gij_up1_lo2_lo3_3x_Ser_p1, gij_up1_lo2_up3_3x_Ser_p1, gij_up1_up2_lo3_3x_Ser_p1, gij_up1_up2_up3_3x_Ser_p1

  //// reorder
  gij_3x_Ser_p1,
  gij_lo1_3x_Ser_p1, gij_up1_3x_Ser_p1,
  gij_lo2_3x_Ser_p1, gij_up2_3x_Ser_p1,
  gij_lo1_lo2_3x_Ser_p1, gij_lo1_up2_3x_Ser_p1,

  //7
  gij_up1_lo2_3x_Ser_p1, gij_up1_up2_3x_Ser_p1,

  //9
  gij_lo3_3x_Ser_p1, gij_up3_3x_Ser_p1,
  
  //11  
  gij_lo1_lo3_3x_Ser_p1, gij_lo1_up3_3x_Ser_p1,

  // 13
  gij_up1_lo3_3x_Ser_p1, gij_up1_up3_3x_Ser_p1, 

  //15
  gij_lo2_lo3_3x_Ser_p1, gij_lo2_up3_3x_Ser_p1,

  //17
  gij_up2_lo3_3x_Ser_p1, gij_up2_up3_3x_Ser_p1, 

  //19
  gij_lo1_lo2_lo3_3x_Ser_p1, gij_lo1_lo2_up3_3x_Ser_p1, gij_lo1_up2_lo3_3x_Ser_p1,

  //22
  gij_lo1_up2_up3_3x_Ser_p1, gij_up1_lo2_lo3_3x_Ser_p1,
  //24
  gij_up1_lo2_up3_3x_Ser_p1, gij_up1_up2_lo3_3x_Ser_p1, gij_up1_up2_up3_3x_Ser_p1
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


