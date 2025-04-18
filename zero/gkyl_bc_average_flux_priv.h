#pragma once

// Private header for bc_average_flux updater, not for direct use in user code.

#include <gkyl_bc_average_flux.h>
// #include <gkyl_bc_average_flux_kernels.h>
#include <assert.h>

// // Function pointer type for average flux kernels.
typedef void (*average_fluxf_t)(const double *vmap, const double q2Dm,
  const double *phi, const double *phiWall, const double *f, double *fRefl);

typedef struct { average_fluxf_t kernels[3]; } average_fluxf_kern_list;  // For use in kernel tables.
typedef struct { average_fluxf_kern_list list[4]; } edged_average_fluxf_kern_list;

// Serendipity  kernels.
GKYL_CU_D
static const edged_average_fluxf_kern_list ser_average_flux_list[] = {
  { .list={
           { NULL, NULL },
           { NULL, NULL },
           { NULL, NULL },
           { NULL, NULL },
          },
  },
  { .list={
           { NULL, NULL },
           { NULL, NULL },
           { NULL, NULL },
           { NULL, NULL },
          },
  },
};

struct gkyl_bc_average_flux_kernels {
  // average_fluxf_t reflectedf;  // reflectedf kernel.
};

// Primary struct in this updater.
struct gkyl_bc_average_flux {
  int dir; // Direction perpendicular to the sheath boundary.
  int cdim; // Conf-space dimensionality.
  enum gkyl_edge_loc edge; // Lower or upper boundary.
  const struct gkyl_basis *basis; // Phase-space basis.
  bool use_gpu; // Whether to run on GPU.
  struct gkyl_bc_average_flux_kernels *kernels;  // reflectedf kernel.
  struct gkyl_bc_average_flux_kernels *kernels_cu;  // device copy.
  const struct gkyl_range *skin_r, *ghost_r; // Skin and ghost ranges.
  const struct gkyl_velocity_map *vel_map; // Velocity space mapping.
};

void
gkyl_bc_gksheath_choose_reflectedf_kernel_cu(const struct gkyl_basis *basis, enum gkyl_edge_loc edge, struct gkyl_bc_average_flux_kernels *kers);

GKYL_CU_D
static average_fluxf_t
bc_gksheath_choose_reflectedf_kernel(const struct gkyl_basis *basis, enum gkyl_edge_loc edge)
{
  int dim = basis->ndim;
  enum gkyl_basis_type basis_type = basis->b_type;
  int poly_order = basis->poly_order;
  switch (basis_type) {
    case GKYL_BASIS_MODAL_GKHYBRID:
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_average_flux_list[edge].list[dim-2].kernels[poly_order-1];
    default:
      assert(false);
      break;
  }
  return 0;
}

#ifdef GKYL_HAVE_CUDA

/**
 * CUDA device function to apply the sheath BC.

 * @param up BC updater.
 * @param phi Electrostatic potential.
 * @param phi_wall Wall potential.
 * @param distf Distribution function array to apply BC to.
 * @param conf_r Configuration space range (to index phi).
 */
void gkyl_bc_average_flux_advance_cu(const struct gkyl_bc_average_flux *up, const struct gkyl_array *phi,
  const struct gkyl_array *phi_wall, struct gkyl_array *distf, const struct gkyl_range *conf_r);

#endif
