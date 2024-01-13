#include <gkyl_deflate_zsurf.h>
#include <gkyl_rect_grid.h>
#include <assert.h>

#include <gkyl_deflate_zsurf_kernels.h>

typedef void (*deflate_zsurf_kernel)(const double *fld, double* deflated_fld);


typedef struct { deflate_zsurf_kernel kernels[3]; } deflate_zsurf_kernel_list;


GKYL_CU_DH
static const deflate_zsurf_kernel_list ser_deflate_zsurf_kernel_list[] = {
  {NULL, deflate_zsurf_lo_2x_ser_p1_remy, deflate_zsurf_lo_2x_ser_p2_remy},
  {NULL, deflate_zsurf_up_2x_ser_p1_remy, deflate_zsurf_up_2x_ser_p2_remy},
};

struct gkyl_deflate_zsurf{
  const struct gkyl_basis *basis;
  const struct gkyl_basis *deflated_basis;
  const struct gkyl_rect_grid* grid;
  const struct gkyl_rect_grid* deflated_grid;
  deflate_zsurf_kernel kernel;
  int *rem_dirs;
  bool use_gpu;
};

GKYL_CU_DH
static deflate_zsurf_kernel
deflate_zsurf_choose_kernel(const int basis_type, const int edge, const int poly_order)
{
  switch (basis_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_deflate_zsurf_kernel_list[edge].kernels[poly_order];

    default:
      assert(false);
      break;
  }
}





