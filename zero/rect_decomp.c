#include <gkyl_rect_decomp.h>

void
gkyl_create_grid_ranges(const struct gkyl_rect_grid *grid,
  const int *nghost, struct gkyl_range *ext_range, struct gkyl_range *range)
{
  int lower_ext[GKYL_MAX_DIM], upper_ext[GKYL_MAX_DIM];
  int lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  
  for (int i=0; i<grid->ndim; ++i) {
    lower_ext[i] = 1-nghost[i];
    upper_ext[i] = grid->cells[i]+nghost[i];

    lower[i] = 1;
    upper[i] = grid->cells[i];
  }
  gkyl_range_init(ext_range, grid->ndim, lower_ext, upper_ext);
  gkyl_sub_range_init(range, ext_range, lower, upper);
}

void gkyl_create_ranges(const struct gkyl_range *inrange,
  const int *nghost, struct gkyl_range *ext_range, struct gkyl_range *range)
{
  int lower_ext[GKYL_MAX_DIM], upper_ext[GKYL_MAX_DIM];
  int lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  
  for (int i=0; i<inrange->ndim; ++i) {
    lower_ext[i] = inrange->lower[i]-nghost[i];
    upper_ext[i] = inrange->upper[i]+nghost[i];

    lower[i] = inrange->lower[i];
    upper[i] = inrange->upper[i];
  }
  gkyl_range_init(ext_range, inrange->ndim, lower_ext, upper_ext);
  gkyl_sub_range_init(range, ext_range, lower, upper);  
}
