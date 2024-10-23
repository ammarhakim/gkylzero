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
           { NULL, NULL }, // For 1D
           { skin_surf_from_ghost_lower_2x_ser_p1, NULL },
           { skin_surf_from_ghost_lower_3x_ser_p1, NULL },
          },
  },
  { .list={
           { NULL, NULL }, // For 1D
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
  const struct gkyl_basis *basis;
  bool use_gpu;
  struct gkyl_skin_surf_from_ghost_kernels *kernels;  // skin surface from ghost kernel.
  struct gkyl_skin_surf_from_ghost_kernels *kernels_cu;  // device copy.
  const struct gkyl_range *skin_r, *ghost_r;
};

void
gkyl_skin_surf_from_ghost_choose_kernel_cu(const struct gkyl_basis *basis, enum gkyl_edge_loc edge, struct gkyl_skin_surf_from_ghost_kernels *kers);

GKYL_CU_D
static skin_surf_from_ghost_t
skin_surf_from_ghost_choose_kernel(const struct gkyl_basis *basis, enum gkyl_edge_loc edge)
{
  int dim = basis->ndim;
  enum gkyl_basis_type basis_type = basis->b_type;
  int poly_order = basis->poly_order;
  switch (basis_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_skin_surf_from_ghost_list[edge].list[dim-1].kernels[poly_order-1];
    default:
      assert(false);
      break;
  }
  return 0;
}

#ifdef GKYL_HAVE_CUDA

/**
 * CUDA device function to copy the values from the ghost region to the adjacent skin region (boundary cells).
 *
 * @param up Pointer to the boundary condition updater.
 * @param field Array representing the field values to update (currently works only in configuration space).
 */
void gkyl_skin_surf_from_ghost_advance_cu(const struct gkyl_skin_surf_from_ghost *up, struct gkyl_array *field);

#endif
