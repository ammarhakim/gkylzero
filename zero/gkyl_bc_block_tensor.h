#include <gkyl_rect_grid.h>
#include <gkyl_range.h>
#include <gkyl_basis.h>

typedef struct bc_block_tensor  bc_block_tensor;


struct bc_block_tensor*
gkyl_bc_block_tensor_new(const struct gkyl_rect_grid* grid, const struct gkyl_range *range, const struct gkyl_range* range_ext, 
  const struct gkyl_basis* basis, bool use_gpu);


void calc_tensor(struct bc_block_tensor *up, int dir, int edge1, int edge2, const double *ej, const double *e_i, double *tj_i);
