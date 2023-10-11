#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_rect_grid.h>
#include <gkyl_deflate_geo.h>
#include <gkyl_deflate_geo_priv.h>

#include <gkyl_array_ops_priv.h>


gkyl_deflate_geo* gkyl_deflate_geo_new(const struct gkyl_basis *cbasis,
  const struct gkyl_rect_grid *grid, const struct gkyl_rect_grid *deflated_grid, bool use_gpu){

  gkyl_deflate_geo *up = gkyl_malloc(sizeof(gkyl_deflate_geo));
  up->grid = grid;
  up->deflated_grid = deflated_grid;
  up->basis = cbasis;

}


void gkyl_deflate_geo_advance(const gkyl_deflate_geo *up, const struct gkyl_range *crange,
    const struct gkyl_array *field, struct gkyl_array *deflated_field){
  // Inflated grid will always be 3 long in other directions
  // So deflate reange to the cell at index 2 in the ignored directions (z or z and y)
  // Then loop through and call the kernel on the inflated field and deflated field
}

void gkyl_deflate_geo_release(gkyl_deflate_geo* up){
  gkyl_free(up);

}


