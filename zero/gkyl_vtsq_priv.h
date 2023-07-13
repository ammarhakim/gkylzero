#include <gkyl_vtsq.h>
#include <gkyl_dg_bin_ops_priv.h>

struct gkyl_vtsq {
  struct gkyl_dg_bin_op_mem *mem;
  int m1comps, vdim_phys;
  bool use_gpu;
};

void
gkyl_vtsq_advance_cu(struct gkyl_vtsq *up, struct gkyl_basis basis,
  const struct gkyl_array* moms, const struct gkyl_range *range,
  struct gkyl_array* out);
