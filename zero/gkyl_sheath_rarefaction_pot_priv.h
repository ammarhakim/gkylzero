#pragma once

// Private header for sheath_rarefaction_pot updater, not for direct use in user code.

#include <gkyl_sheath_rarefaction_pot.h>
#include <gkyl_sheath_rarefaction_pot_kernels.h>
#include <assert.h>

// Function pointer type for kernels.
typedef void (*rarefaction_phimod_t)(double elem_q, double mElc, const double *momsElc, const double *m2parElc, double mIon, const double *momsIon, const double *m2parIon, const double *phiWall, double *phi);

typedef struct { rarefaction_phimod_t kernels[3]; } rarefaction_phimod_kern_list;  // For use in kernel tables.
typedef struct { rarefaction_phimod_kern_list list[4]; } edged_rarefaction_phimod_kern_list;

// Serendipity  kernels.
GKYL_CU_D
static const edged_rarefaction_phimod_kern_list ser_sheath_rarepot_list[] = {
  { .list={
           { sheath_rarefaction_phi_mod_lower_1x_ser_p1, sheath_rarefaction_phi_mod_lower_1x_ser_p2 },
           { NULL, NULL},
//           { sheath_rarefaction_phi_mod_lower_3x_ser_p1, sheath_rarefaction_phi_mod_lower_3x_ser_p2 },
           { NULL, NULL},
          },
  },
  { .list={
           { sheath_rarefaction_phi_mod_upper_1x_ser_p1, sheath_rarefaction_phi_mod_upper_1x_ser_p2 },
           { NULL, NULL},
//           { sheath_rarefaction_phi_mod_upper_3x_ser_p1, sheath_rarefaction_phi_mod_upper_3x_ser_p2 },
           { NULL, NULL},
          },
  },
};

// Tensor kernels.
GKYL_CU_D
static const edged_rarefaction_phimod_kern_list tensor_sheath_rarepot_list[] = {
  { .list={
           { NULL, sheath_rarefaction_phi_mod_lower_1x_tensor_p2 },
           { NULL, NULL},
           { NULL, NULL},
          },
  },
  { .list={
           { NULL, sheath_rarefaction_phi_mod_upper_1x_tensor_p2 },
           { NULL, NULL},
           { NULL, NULL},
          },
  },
};

struct gkyl_sheath_rarefaction_pot_kernels {
  rarefaction_phimod_t phimod;
};

// Primary struct in this updater.
struct gkyl_sheath_rarefaction_pot {
  enum gkyl_edge_loc edge;
  struct gkyl_range skin_r;
  const struct gkyl_basis *basis;
  bool use_gpu;
  double elem_charge; // elementary charge.
  double mass_e, mass_i; // electron and ion mass.
  struct gkyl_sheath_rarefaction_pot_kernels *kernels;
  struct gkyl_sheath_rarefaction_pot_kernels *kernels_cu;  // device copy.
  const struct gkyl_rect_grid *grid;
};

void
gkyl_sheath_rarepot_choose_phimod_kernel_cu(const struct gkyl_basis *basis, enum gkyl_edge_loc edge, struct gkyl_sheath_rarefaction_pot_kernels *kers);

GKYL_CU_D
static rarefaction_phimod_t
sheath_rarepot_choose_phimod_kernel(const struct gkyl_basis *basis, enum gkyl_edge_loc edge)
{
  int dim = basis->ndim;
  enum gkyl_basis_type basis_type = basis->b_type;
  int poly_order = basis->poly_order;
  switch (basis_type) {
    case GKYL_BASIS_MODAL_GKHYBRID:
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_sheath_rarepot_list[edge].list[dim-2].kernels[poly_order-1];
    case GKYL_BASIS_MODAL_TENSOR:
      return tensor_sheath_rarepot_list[edge].list[dim-2].kernels[poly_order-1];
    default:
      assert(false);
      break;
  }
}

#ifdef GKYL_HAVE_CUDA

/**
 * CUDA device function to modify the potential at the boundary.

 * @param up updater object.
 * @param moms_e first threee moments of the electron distribution.
 * @param m2par_e v_par^2 moment of the electron distribution.
 * @param moms_i first threee moments of the ion distribution.
 * @param m2par_i v_par^2 moment of the ion distribution.
 * @param phi_wall Wall potential.
 * @param phi Electrostatic potential.
 */
void gkyl_sheath_rarefaction_pot_advance_cu(const struct gkyl_sheath_rarefaction_pot *up,
  const struct gkyl_array *moms_e, const struct gkyl_array *m2par_e,
  const struct gkyl_array *moms_i, const struct gkyl_array *m2par_i,
  const struct gkyl_array *phi_wall, struct gkyl_array *phi);

#endif
