#pragma once

#include <gkyl_gk_fluid_moments.h>
#include <gkyl_dg_bin_ops.h>

struct gkyl_gk_fluid_moments {
  struct gkyl_basis basis;
  int num_mom, vdim, component;
  long num_cells; // Batch size of the weak division solves.
  double mass, charge, cell_factor, factor;
  struct gkyl_array *bmag, *upar, *ppar, *pperp, *tmp, *result;
  gkyl_dg_bin_op_mem *mem;
  void (*calc)(const struct gkyl_gk_fluid_moments *, const struct gkyl_range *, const struct gkyl_array *, const struct gkyl_array *);
  void (*finish)(const struct gkyl_gk_fluid_moments *, const struct gkyl_range *, struct gkyl_array *);
  void (*perp)(const struct gkyl_gk_fluid_moments *, const struct gkyl_range *, const struct gkyl_array *);
  void (*project_perp)(const struct gkyl_gk_fluid_moments *, const struct gkyl_range *, const struct gkyl_array *, const struct gkyl_array *);
};

#ifdef GKYL_HAVE_CUDA
/** Extract constant coefficients as cell integrals on the GPU.
 * @param up Fluid updater with DG results in its scratch array.
 * @param range Configuration subrange.
 * @param out Output cell integrals.
 */
void gkyl_gk_fluid_moments_integrate_cu(
  const struct gkyl_gk_fluid_moments *up, const struct gkyl_range *range, struct gkyl_array *out
);
#endif
