#pragma once

#include <gkyl_mom_type.h>
#include <gkyl_rect_grid.h>

struct gkyl_mom_calc_bcorr {
  struct gkyl_rect_grid grid; // grid object
  int ndim; // number of dimensions
  int num_basis; // number of basis functions
  int num_up_dirs; // number of update directions
  int update_dirs[GKYL_MAX_DIM]; // directions to update
  int space;
  const struct gkyl_mom_type *momt; // moment type object

  uint32_t flags;
  struct gkyl_mom_calc_bcorr *on_dev;
};
