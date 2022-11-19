#include <gkyl_bgk_collisions.h>
#include <gkyl_dg_bin_ops_priv.h>

struct gkyl_bgk_collisions {
  int cdim; // Configuration-space dimension.
  int pdim; // Phase-space dimension.
  int cnum_basis; // Number of conf-space basis functions.
  int pnum_basis; // Number of phase-space basis functions.

  mul_op_t mul_op; // Conf*phase bin_op multiplication kernel.

  bool use_gpu;
};

void
gkyl_bgk_collisions_advance_cu(const gkyl_bgk_collisions *up,
  const struct gkyl_range *crange, const struct gkyl_range *prange,
  const struct gkyl_array *nu, const struct gkyl_array *nufM, const struct gkyl_array *fin,
  struct gkyl_array *out, struct gkyl_array *cflfreq);
