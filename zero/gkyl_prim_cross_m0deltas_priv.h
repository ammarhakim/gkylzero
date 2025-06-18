#include <gkyl_prim_cross_m0deltas.h>
#include <gkyl_dg_bin_ops_priv.h>

struct gkyl_prim_cross_m0deltas {
  struct gkyl_dg_bin_op_mem *mem;
  double betap1;
  bool use_gpu;
};

void
gkyl_prim_cross_m0deltas_advance_cu(gkyl_prim_cross_m0deltas *up, struct gkyl_basis basis,
  double massself, const struct gkyl_array* m0self, const struct gkyl_array* nuself,
  double massother, const struct gkyl_array* m0other, const struct gkyl_array* nuother,
  const struct gkyl_array* prem0s, const struct gkyl_range *range, struct gkyl_array* out);
