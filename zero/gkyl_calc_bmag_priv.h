#include <gkyl_calc_bmag.h>
#include <gkyl_rect_grid.h>
#include <assert.h>

struct gkyl_calc_bmag {
  struct gkyl_basis* cbasis; //comp basis
  struct gkyl_basis* pbasis; //physical RZ basis
  struct gkyl_rect_grid* cgrid; // computational grid
  struct gkyl_rect_grid* pgrid; // physical RZ grid
  bool use_gpu;
  struct gkyl_bmag_ctx* bmag_ctx;
};

