#pragma once

// Private header for skin_surf_from_ghost updater, not for direct use in user code.

#include <gkyl_skin_surf_from_ghost.h>
#include <gkyl_skin_surf_from_ghost_kernels.h>
#include <assert.h>

// Function pointer type for the kernel setting the skin value equal to the ghost value at surface
typedef void (*skin_surf_from_ghost_t)(const double *fghost, double *fskin);

typedef struct { skin_surf_from_ghost_t kernels[2]; } skin_surf_from_ghost_kern_list;  // For use in kernel tables.
typedef struct { skin_surf_from_ghost_kern_list list[3]; } edged_skin_surf_from_ghost_kern_list;

// Serendipity  kernels.
GKYL_CU_D
static const edged_skin_surf_from_ghost_kern_list ser_skin_surf_from_ghost_list[] = {
  { .list={
           { NULL, NULL }, // For 1D, to be implemented
           { skin_surf_from_ghost_lower_2x_ser_p1, NULL },
           { skin_surf_from_ghost_lower_3x_ser_p1, NULL },
          },
  },
  { .list={
           { NULL, NULL }, // For 1D, to be implemented
           { skin_surf_from_ghost_upper_2x_ser_p1, NULL },
           { skin_surf_from_ghost_upper_3x_ser_p1, NULL },
          },
  },
};
 
struct gkyl_skin_surf_from_ghost_kernels {
  skin_surf_from_ghost_t ghost_to_skin;
};

// Primary struct in this updater.
struct gkyl_skin_surf_from_ghost {
  int dir, cdim;
  enum gkyl_edge_loc edge;
  bool use_gpu;
  struct gkyl_skin_surf_from_ghost_kernels *kernels;  // skin surface from ghost kernel.
  struct gkyl_skin_surf_from_ghost_kernels *kernels_cu;  // device copy.
  const struct gkyl_range *skin_r, *ghost_r;
};

#ifdef GKYL_HAVE_CUDA

void skin_surf_from_ghost_choose_kernel_cu(const struct gkyl_basis basis, enum gkyl_edge_loc edge, 
  struct gkyl_skin_surf_from_ghost_kernels *kers);

void skin_surf_from_ghost_advance_cu(const struct gkyl_skin_surf_from_ghost *up, struct gkyl_array *field);
#endif

GKYL_CU_D
static void skin_surf_from_ghost_choose_kernel(const struct gkyl_basis basis, enum gkyl_edge_loc edge, 
  bool use_gpu, struct gkyl_skin_surf_from_ghost_kernels *kernels)
{
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    skin_surf_from_ghost_choose_kernel_cu(basis, edge, kernels);
    return;
  }
#endif

  int dim = basis.ndim;
  enum gkyl_basis_type basis_type = basis.b_type;
  int poly_order = basis.poly_order;
  switch (basis_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      kernels->ghost_to_skin = ser_skin_surf_from_ghost_list[edge].list[dim-1].kernels[poly_order-1];
      break;
    default:
      assert(false);
      break;
  }
}
