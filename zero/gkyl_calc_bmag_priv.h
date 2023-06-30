#include <gkyl_calc_bmag.h>
#include <gkyl_rect_grid.h>
#include <gkyl_calc_bmag_kernels.h>
#include <assert.h>

typedef void (*bmag_kernel)(const double **xyz, double *gFld);

typedef struct { bmag_kernel kernels[3]; } bmag_kernel_list;  // For use in kernel tables.

GKYL_CU_DH
static const bmag_kernel_list ser_bmag_kernel_list[] = {
  { NULL, NULL, NULL }, // 0x No 0D basis functions
  { NULL, NULL, NULL}, // 1x Not tested yet
  { NULL, bmag_2x_Ser_p1, gij_2x_Ser_p2}, //Only 2x makes sense
  { NULL, NULL, NULL}
};

struct gkyl_calc_bmag {
  unsigned cdim; // Configuration-space dimension.
  unsigned cnum_basis; // Number of conf-space basis functions.
  unsigned poly_order; // Polynomial order of the basis.
  struct gkyl_rect_grid* grid;
  bool use_gpu;
  bmag_kernel kernel;
  void *bmag_ctx;
  void (*mapc2p)(gkgeom_app *app, const double *xn, double *ret);
  evalf_t psi;
  evalf_t psibyr;
  evalf_t psibyr2;
};

GKYL_CU_DH
static bmag_kernel
bmag_choose_kernel(int dim, int basis_type, int poly_order)
{
  switch (basis_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_bmag_kernel_list[dim].kernels[poly_order];
    default:
      assert(false);
      break;
  }
}




