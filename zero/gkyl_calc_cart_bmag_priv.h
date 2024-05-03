#include <gkyl_calc_cart_bmag.h>
#include <gkyl_rect_grid.h>
#include <gkyl_calc_cart_bmag_kernels.h>
#include <assert.h>

typedef void (*cart_bmag_kernel)( const double **psibyr, const double *psibyr2, double *bmag_out, double *br_out, double *bz_out, double scale_factorR, double scale_factorZ);

typedef struct { cart_bmag_kernel kernels[3]; } cart_bmag_kernel_list;  // For use in kernel tables.

GKYL_CU_DH
static const cart_bmag_kernel_list ser_cart_bmag_kernel_list[] = {
  { NULL, NULL, NULL }, // 0x No 0D basis functions
  { NULL, NULL, NULL}, // 1x Not tested yet
  { NULL, cart_bmag_2x_Ser_p1, cart_bmag_2x_Ser_p2}, //Only 2x makes sense
  { NULL, NULL, NULL}
};

GKYL_CU_DH
static const bmag_kernel_list tensor_cart_bmag_kernel_list[] = {
  { NULL, NULL, NULL }, // 0x No 0D basis functions
  { NULL, NULL, NULL}, // 1x Not tested yet
  { NULL, cart_bmag_2x_Tensor_p1, cart_bmag_2x_Tensor_p2}, //Only 2x makes sense
  { NULL, NULL, NULL}
};

struct gkyl_calc_cart_bmag {
  const struct gkyl_basis* basis; //physical RZ basis
  const struct gkyl_rect_grid* grid; // physical RZ grid
  struct gkyl_range local; // physical rz range
  struct gkyl_range local_ext; // extended physical rz range
  bool use_gpu;
  cart_bmag_kernel kernel;
  struct gkyl_array *bmag_rz;
  struct gkyl_array *br_rz;
  struct gkyl_array *bz_rz;
};

GKYL_CU_DH
static cart_bmag_kernel
cart_bmag_choose_kernel(int dim, int basis_type, int poly_order)
{
  switch (basis_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_cart_bmag_kernel_list[dim].kernels[poly_order];
    case GKYL_BASIS_MODAL_TENSOR:
      return ser_cart_bmag_kernel_list[dim].kernels[poly_order];
    default:
      assert(false);
      break;
  }
  return 0;
}




