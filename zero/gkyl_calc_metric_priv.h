#include <gkyl_calc_metric.h>
#include <gkyl_rect_grid.h>
#include <assert.h>
#include <math.h>
GKYL_CU_DH

struct gkyl_calc_metric {
  unsigned cdim; // Configuration-space dimension.
  unsigned cnum_basis; // Number of conf-space basis functions.
  unsigned poly_order; // Polynomial order of the basis.
  const struct gkyl_rect_grid* grid;
  bool use_gpu;
  const int *num_cells;
  const struct gkyl_basis* cbasis;
};
