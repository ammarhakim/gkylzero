#pragma once

#include <gkyl_rect_grid.h>

struct gkyl_mom_pkpm_surf_calc {
  struct gkyl_rect_grid vel_grid;
  double mass;

  uint32_t flags;
  struct gkyl_mom_pkpm_surf_calc *on_dev; // pointer to itself or device data
};
