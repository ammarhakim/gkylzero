#include <gkyl_util.h>
typedef struct gkyl_gk_geometry gkyl_gk_geometry;

gkyl_gk_geometry* gkyl_gk_geometry_new(const struct gkyl_rect_grid* grid, const struct gkyl_range *range, const struct gkyl_range* range_ext, const struct gkyl_basis* basis, evalf_t mapc2p_func, evalf_t bmag_func);

void gkyl_gk_geometry_advance(struct gkyl_gk_geometry* up);

void gkyl_gk_geometry_release(struct gkyl_gk_geometry* up);
