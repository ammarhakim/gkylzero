#include <gkyl_calc_metric.h>
#include <gkyl_rect_grid.h>
#include <gkyl_calc_metric_kernels.h>
#include <assert.h>

typedef void (*metric_kernel)(const double **xyz, double *gFld);

typedef struct { metric_kernel kernels[3]; } metric_kernel_list;  // For use in kernel tables.

GKYL_CU_DH
static const metric_kernel_list ser_metric_kernel_list[] = {
  { NULL, NULL, NULL }, // 0x No 0D basis functions
  { NULL, NULL, NULL}, // 1x Not tested yet
  { NULL, gij_2x_Ser_p1, gij_2x_Ser_p2},
  { NULL, gij_3x_Ser_p1, gij_3x_Ser_p2}
};

struct gkyl_calc_metric {
  unsigned cdim; // Configuration-space dimension.
  unsigned cnum_basis; // Number of conf-space basis functions.
  unsigned poly_order; // Polynomial order of the basis.
  struct gkyl_rect_grid* grid;
  bool use_gpu;
  metric_kernel kernel;
};

GKYL_CU_DH
static metric_kernel
metric_choose_kernel(int dim, int basis_type, int poly_order)
{
  switch (basis_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_metric_kernel_list[dim].kernels[poly_order];
    default:
      assert(false);
      break;
  }
}




