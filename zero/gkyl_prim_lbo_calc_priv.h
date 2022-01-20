#pragma once

#include <gkyl_prim_lbo.h>
#include <gkyl_rect_grid.h>

struct gkyl_prim_lbo_calc {
  struct gkyl_rect_grid grid;
  const struct gkyl_prim_lbo *prim;

  uint32_t flags;
  struct gkyl_prim_lbo_calc *on_dev; // pointer to itself or device data
};
