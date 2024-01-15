#pragma once
#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_basis.h>
#include <gkyl_fem_poisson_bctype.h>


// Updater type
struct gkyl_line_fem_poisson {
  struct gkyl_rect_grid grid;
  struct gkyl_basis basis;
  struct gkyl_range local;
  struct gkyl_range local_ext;
  struct gkyl_array *epsilon;
  struct gkyl_poisson_bc poisson_bc;
  bool use_gpu;
};
