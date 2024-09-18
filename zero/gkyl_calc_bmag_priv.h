#include <gkyl_calc_bmag.h>
#include <gkyl_rect_grid.h>
#include <assert.h>

struct bmag_ctx{
   const struct gkyl_rect_grid* grid;
   const struct gkyl_rect_grid* cgrid;
   const struct gkyl_range* range;
   const struct gkyl_range* crange;
   const struct gkyl_range* crange_global;
   const struct gkyl_basis* basis;
   const struct gkyl_basis* cbasis;
   const struct gkyl_array* bmagdg;
   const struct gkyl_array* mapc2p;
};

struct gkyl_calc_bmag {
  const struct gkyl_basis* cbasis; //comp basis
  const struct gkyl_basis* pbasis; //physical RZ basis
  const struct gkyl_rect_grid* cgrid; // computational grid
  const struct gkyl_rect_grid* pgrid; // physical RZ grid
  bool use_gpu;
  bmag_ctx* bmag_ctx;
};

