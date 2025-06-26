#pragma once

// Private header for sheath_rarefaction_pot updater, not for direct use in user code.

#include <gkyl_sheath_rarefaction_pot.h>
#include <gkyl_sheath_rarefaction_pot_kernels.h>
#include <assert.h>

// Function pointer type for kernels.
typedef void (*rarefaction_phimod_t)(double elem_q, double mElc, const double *momsElc,
  double mIon, const double *momsIon, const double *phiWall, double *phi);

typedef struct { rarefaction_phimod_t kernels[2]; } rarefaction_phimod_kern_list;  // For use in kernel tables.
typedef struct { rarefaction_phimod_kern_list list[3]; } edged_rarefaction_phimod_kern_list;

// Serendipity  kernels.
GKYL_CU_D
static const edged_rarefaction_phimod_kern_list ser_sheath_rarepot_list[] = {
  { .list={
           { sheath_rarefaction_phi_mod_lower_1x_ser_p1, sheath_rarefaction_phi_mod_lower_1x_ser_p2 },
           { sheath_rarefaction_phi_mod_lower_2x_ser_p1, sheath_rarefaction_phi_mod_lower_2x_ser_p2 },
           { sheath_rarefaction_phi_mod_lower_3x_ser_p1, sheath_rarefaction_phi_mod_lower_3x_ser_p2 },
          },
  },
  { .list={
           { sheath_rarefaction_phi_mod_upper_1x_ser_p1, sheath_rarefaction_phi_mod_upper_1x_ser_p2 },
           { sheath_rarefaction_phi_mod_upper_2x_ser_p1, sheath_rarefaction_phi_mod_upper_2x_ser_p2 },
           { sheath_rarefaction_phi_mod_upper_3x_ser_p1, sheath_rarefaction_phi_mod_upper_3x_ser_p2 },
          },
  },
};

// Tensor kernels.
GKYL_CU_D
static const edged_rarefaction_phimod_kern_list tensor_sheath_rarepot_list[] = {
  { .list={
           { NULL, NULL },
           { NULL, NULL },
           { NULL, NULL },
          },
  },
  { .list={
           { NULL, NULL },
           { NULL, NULL },
           { NULL, NULL },
          },
  },
};

struct gkyl_sheath_rarefaction_pot_kernels {
  rarefaction_phimod_t phimod;
};

// Primary struct in this updater.
struct gkyl_sheath_rarefaction_pot {
  bool use_gpu;
  double elem_charge; // elementary charge.
  double mass_e, mass_i; // electron and ion mass.
  struct gkyl_sheath_rarefaction_pot_kernels *kernels;
};

#ifdef GKYL_HAVE_CUDA
void
gkyl_sheath_rarepot_choose_phimod_kernel_cu(const struct gkyl_basis *basis, enum gkyl_edge_loc edge,
  struct gkyl_sheath_rarefaction_pot_kernels *kers);
#endif

GKYL_CU_D
static void
sheath_rarepot_choose_phimod_kernel(const struct gkyl_basis *basis, enum gkyl_edge_loc edge,
  struct gkyl_sheath_rarefaction_pot_kernels *kers, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    gkyl_sheath_rarepot_choose_phimod_kernel_cu(basis, edge, kers);
    return;
  }
#endif

  int dim = basis->ndim;
  enum gkyl_basis_type basis_type = basis->b_type;
  int poly_order = basis->poly_order;
  switch (basis_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      kers->phimod = ser_sheath_rarepot_list[edge].list[dim-1].kernels[poly_order-1];
      break;
//    case GKYL_BASIS_MODAL_TENSOR:
//      kers->phimod = tensor_sheath_rarepot_list[edge].list[dim-1].kernels[poly_order-1];
//      break;
    default:
      assert(false);
      break;
  }

  assert(kers->phimod);
}

#ifdef GKYL_HAVE_CUDA

/**
 * CUDA device function to modify the potential at the boundary.

 * @param up updater object.
 * @param skin_range Range of the skin where potential is modified.
 * @param surf_range Surface range the wall potential lives on.
 * @param moms_e first four moments (m0, m1, m2par, m2perp) of the electron distribution.
 * @param moms_i first four moments (m0, m1, m2par, m2perp) of the ion distribution.
 * @param phi_wall Wall potential.
 * @param phi Electrostatic potential.
 */
void gkyl_sheath_rarefaction_pot_advance_cu(const struct gkyl_sheath_rarefaction_pot *up,
  const struct gkyl_range *skin_range, const struct gkyl_range *surf_range,
  const struct gkyl_array *moms_e, const struct gkyl_array *moms_i,
  const struct gkyl_array *phi_wall, struct gkyl_array *phi);

#endif
