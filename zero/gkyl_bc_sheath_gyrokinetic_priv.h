#pragma once

// Private header for bc_sheath_gyrokinetic updater, not for direct use in user code.

#include <gkyl_bc_sheath_gyrokinetic.h>
#include <gkyl_bc_sheath_gyrokinetic_kernels.h>
#include <assert.h>

// Function pointer type for sheath reflection kernels.
typedef void (*sheath_reflectedf_t)(const double wv, const double dv,
  const double vlowerSq, const double vupperSq, const double q2Dm,
  const double *phi, const double *phiWall, const double *f, double *fRefl);

typedef struct { sheath_reflectedf_t kernels[3]; } sheath_reflectedf_kern_list;  // For use in kernel tables.
typedef struct { sheath_reflectedf_kern_list list[4]; } edged_sheath_reflectedf_kern_list;

// Serendipity  kernels.
GKYL_CU_D
static const edged_sheath_reflectedf_kern_list ser_sheath_reflect_list[] = {
  { .list={
           { bc_sheath_gyrokinetic_reflectedf_lower_1x1v_ser_p1, bc_sheath_gyrokinetic_reflectedf_lower_1x1v_ser_p2 },
           { bc_sheath_gyrokinetic_reflectedf_lower_1x2v_ser_p1, bc_sheath_gyrokinetic_reflectedf_lower_1x2v_ser_p2 },
           { NULL, NULL},
           { bc_sheath_gyrokinetic_reflectedf_lower_3x2v_ser_p1, bc_sheath_gyrokinetic_reflectedf_lower_3x2v_ser_p2 },
          },
  },
  { .list={
           { bc_sheath_gyrokinetic_reflectedf_upper_1x1v_ser_p1, bc_sheath_gyrokinetic_reflectedf_upper_1x1v_ser_p2 },
           { bc_sheath_gyrokinetic_reflectedf_upper_1x2v_ser_p1, bc_sheath_gyrokinetic_reflectedf_upper_1x2v_ser_p2 },
           { NULL, NULL},
           { bc_sheath_gyrokinetic_reflectedf_upper_3x2v_ser_p1, bc_sheath_gyrokinetic_reflectedf_upper_3x2v_ser_p2 },
          },
  },
};

// Serendipity  kernels.
GKYL_CU_D
static const edged_sheath_reflectedf_kern_list tensor_sheath_reflect_list[] = {
  { .list={
           { NULL, bc_sheath_gyrokinetic_reflectedf_lower_1x1v_tensor_p2 },
           { NULL, bc_sheath_gyrokinetic_reflectedf_lower_1x2v_tensor_p2 },
           { NULL, NULL},
           { NULL, NULL},
          },
  },
  { .list={
           { NULL, bc_sheath_gyrokinetic_reflectedf_upper_1x1v_tensor_p2 },
           { NULL, bc_sheath_gyrokinetic_reflectedf_upper_1x2v_tensor_p2 },
           { NULL, NULL},
           { NULL, NULL},
          },
  },
};

// Primary struct in this updater.
struct gkyl_bc_sheath_gyrokinetic {
  int dir, cdim;
  enum gkyl_edge_loc edge;
  struct gkyl_range skin_r, ghost_r;
  const struct gkyl_basis *basis;
  bool use_gpu;
  double q2Dm; // charge-to-mass ratio times 2.
  sheath_reflectedf_t ker_reflectedf;  // reflectedf kernel.
  const struct gkyl_rect_grid *grid;
  struct gkyl_range conf_r;
};

GKYL_CU_D
static sheath_reflectedf_t
bc_gksheath_choose_reflectedf_kernel(const int dim, const int basis_type, const int poly_order, enum gkyl_edge_loc edge)
{
  switch (basis_type) {
    case GKYL_BASIS_MODAL_GKHYBRID:
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_sheath_reflect_list[edge].list[dim-2].kernels[poly_order-1];
    case GKYL_BASIS_MODAL_TENSOR:
      return tensor_sheath_reflect_list[edge].list[dim-2].kernels[poly_order-1];
    default:
      assert(false);
      break;
  }
}

GKYL_CU_D
static void
bc_gksheath_reflect(int dir, const struct gkyl_basis *basis, int cdim, double *out, const double *inp)
{
  basis->flip_odd_sign(dir, inp, out);
  basis->flip_odd_sign(cdim, out, out); // cdim is the vpar direction.
}

