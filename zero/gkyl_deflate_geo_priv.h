#include <gkyl_deflate_geo.h>
#include <gkyl_rect_grid.h>
#include <assert.h>

#include <gkyl_deflate_geo_kernels.h>

typedef void (*deflate_geo_kernel)(const double *fld, double* deflated_fld);


typedef struct { deflate_geo_kernel kernels[3]; } deflate_geo_kernel_list;  // For use in kernel tables.


GKYL_CU_DH
static const deflate_geo_kernel_list ser_deflate_geo_kernel_list[] = {
  { NULL, NULL, NULL }, // 0x No 0D basis functions
  { NULL, deflate_geo_1x_Ser_p1, deflate_geo_Ser_1x_p2}, // 1x Not tested yet
  { NULL, deflate_geo_2x_Ser_p1, deflate_geo_Ser_2x_p2}, // 1x Not tested yet
  { NULL, NULL, NULL} // don't need for 3x
};

struct gkyl_deflate_geo{
  const struct gkyl_basis *basis;
  const struct gkyl_basis *deflated_basis;
  const struct gkyl_rect_grid* grid;
  const struct gkyl_rect_grid* deflated_grid;
  bool use_gpu;
};

GKYL_CU_DH
static deflate_geo_kernel
deflate_geo_choose_kernel(int dim, int basis_type, int poly_order)
{
  switch (basis_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_deflate_geo_kernel_list[dim].kernels[poly_order];
    default:
      assert(false);
      break;
  }
}





