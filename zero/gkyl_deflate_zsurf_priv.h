// Private header: not for direct use
#pragma once

#include <math.h>
#include <gkyl_deflate_zsurf_kernels.h>
#include <assert.h>

typedef void (*deflate_zsurf_kernel)(const double *fld, double* deflated_fld);

typedef struct { deflate_zsurf_kernel kernels[3]; } deflate_zsurf_kernel_list;

GKYL_CU_D
static const deflate_zsurf_kernel_list ser_deflate_zsurf_kernel_list[] = {
  {NULL, deflate_zsurf_lo_2x_ser_p1_remy, deflate_zsurf_lo_2x_ser_p2_remy},
  {NULL, deflate_zsurf_up_2x_ser_p1_remy, deflate_zsurf_up_2x_ser_p2_remy},
};

struct gkyl_deflate_zsurf {
  int num_basis; // Number of basis functions in full basis
  int num_deflated_basis; // Number of basis functions in deflated basis
  deflate_zsurf_kernel kernel;
  uint32_t flags;
  struct gkyl_deflate_zsurf *on_dev; // pointer to itself or device data
};

GKYL_CU_D
static deflate_zsurf_kernel
deflate_zsurf_choose_kernel(enum gkyl_basis_type basis_type, int edge, int poly_order)
{
  switch (basis_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_deflate_zsurf_kernel_list[edge].kernels[poly_order];

    default:
      assert(false);
      break;
  }
}




