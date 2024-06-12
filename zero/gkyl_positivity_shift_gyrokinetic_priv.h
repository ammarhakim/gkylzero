#pragma once

// Private header for positivity_shift_gyrokinetic updater, not for direct use in user code.

#include <gkyl_positivity_shift_gyrokinetic.h>
#include <gkyl_positivity_shift_gyrokinetic_kernels.h>
#include <gkyl_mom_gyrokinetic_kernels.h>
#include <gkyl_mom_gyrokinetic_priv.h>
#include <gkyl_util.h>
#include <assert.h>

// Function pointer type for sheath reflection kernels.
typedef bool (*shift_t)(double ffloor, double *distf, double *Deltaf);
typedef void (*intmom_t)(const double *dxv, const double *vmap,
  double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out);

typedef struct { shift_t kernels[3]; } pos_shift_gk_kern_list_shift;  // For use in kernel tables.
typedef struct { intmom_t kernels[3]; } pos_shift_gk_kern_list_intmom;  // For use in kernel tables.

// Serendipity  kernels.
GKYL_CU_D
static const pos_shift_gk_kern_list_shift pos_shift_gk_kern_list_shift_ser[] = {
  { positivity_shift_gyrokinetic_1x1v_ser_p1, NULL, NULL },
  { positivity_shift_gyrokinetic_1x2v_ser_p1, NULL, NULL },
  { positivity_shift_gyrokinetic_2x2v_ser_p1, NULL, NULL },
  { positivity_shift_gyrokinetic_3x2v_ser_p1, NULL, NULL },
};

GKYL_CU_D
static const pos_shift_gk_kern_list_intmom pos_shift_gk_kern_list_intmom_ser[] = {
  { gyrokinetic_int_mom_1x1v_ser_p1, NULL, NULL },
  { gyrokinetic_int_mom_1x2v_ser_p1, NULL, NULL },
  { gyrokinetic_int_mom_2x2v_ser_p1, NULL, NULL },
  { gyrokinetic_int_mom_3x2v_ser_p1, NULL, NULL },
};

struct gkyl_positivity_shift_gyrokinetic_kernels {
  shift_t shift;  // Kernel that reflects f and computes Deltaf.
  intmom_t int_mom;  // Kernel that computes integrated moments.
};

// Primary struct in this updater.
struct gkyl_positivity_shift_gyrokinetic {
  int num_basis;  // Number of phase-space basis monomials.
  struct gkyl_rect_grid grid;  // Phase-space grid.
  double mass;  // Species mass.
  double *ffloor;  // Minimum f to shift distribution to when it's <0.
  double ffloor_fac;  // ffloor = max(f)*ffloor_fac.
  double cellav_fac; // Factor multiplying 0th DG coefficient to give cellav.
  const struct gk_geometry *gk_geom; // Pointer to geometry object.
  const struct gkyl_velocity_map *vel_map; // Pointer to velocity mapping object.
  bool use_gpu;
  struct gkyl_positivity_shift_gyrokinetic_kernels *kernels;
};

#ifdef GKYL_HAVE_CUDA
// Declaration of cuda device functions.

void pos_shift_gk_choose_shift_kernel_cu(struct gkyl_positivity_shift_gyrokinetic_kernels *kernels,
  struct gkyl_basis pbasis);

void gkyl_positivity_shift_gyrokinetic_advance_cu(gkyl_positivity_shift_gyrokinetic* up,
  const struct gkyl_range *phase_rng, const struct gkyl_range *conf_rng,
  struct gkyl_array *GKYL_RESTRICT distf, struct gkyl_array *GKYL_RESTRICT mom);
#endif

GKYL_CU_D
static void pos_shift_gk_choose_shift_kernel(struct gkyl_positivity_shift_gyrokinetic_kernels *kernels,
  struct gkyl_basis pbasis, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    pos_shift_gk_choose_shift_kernel_cu(kernels, pbasis);
    return;
  }
#endif

  enum gkyl_basis_type basis_type = pbasis.b_type;
  int pdim = pbasis.ndim;
  int poly_order = pbasis.poly_order;

  switch (basis_type) {
    case GKYL_BASIS_MODAL_GKHYBRID:
    case GKYL_BASIS_MODAL_SERENDIPITY:
      kernels->shift = pos_shift_gk_kern_list_shift_ser[pdim-2].kernels[poly_order-1];
      kernels->int_mom = pos_shift_gk_kern_list_intmom_ser[pdim-2].kernels[poly_order-1];
      break;
    default:
      assert(false);
      break;
  }
}
