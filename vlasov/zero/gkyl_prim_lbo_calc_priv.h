#pragma once

#include <gkyl_prim_lbo_type.h>
#include <gkyl_rect_grid.h>

struct gkyl_prim_lbo_calc {
  struct gkyl_rect_grid grid;
  const struct gkyl_prim_lbo_type *prim;

  bool is_first; // flag to indicate first call to update
  struct gkyl_nmat *As, *xs; // matrices for LHS and RHS
  gkyl_nmat_mem *mem; // memory for use in batched linear solve

  uint32_t flags;
  struct gkyl_prim_lbo_calc *on_dev; // pointer to itself or device data
};
