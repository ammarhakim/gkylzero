#pragma once
#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_basis.h>
#include <gkyl_fem_poisson_bctype.h>

struct deflated_fem_data {
  struct gkyl_array *deflated_field;
  struct gkyl_fem_poisson *fem_poisson;
};

// Updater type
struct gkyl_line_fem_poisson {
  struct gkyl_rect_grid grid;
  struct gkyl_basis basis;
  struct gkyl_range local;
  struct gkyl_range local_ext;
  struct gkyl_array *epsilon;
  struct gkyl_poisson_bc poisson_bc;
  struct deflated_fem_data *d_fem_data;
  bool use_gpu;
};
