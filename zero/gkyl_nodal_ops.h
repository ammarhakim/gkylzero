#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_basis.h>

void gkyl_nodal_ops_n2m(const struct gkyl_basis* cbasis, const struct gkyl_rect_grid *grid, const struct gkyl_range *nrange, const struct gkyl_range *update_range, int num_comp, const struct gkyl_array* nodal_fld, struct gkyl_array *modal_fld);


void gkyl_nodal_ops_m2n(const struct gkyl_basis* cbasis, const struct gkyl_rect_grid *grid, const struct gkyl_range *nrange, const struct gkyl_range *update_range, int num_comp, struct gkyl_array* nodal_fld, const struct gkyl_array *modal_fld);
