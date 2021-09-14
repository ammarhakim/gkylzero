// Private header, not for direct use as public API
#pragma once

#include <gkyl_alloc.h>
#include <gkyl_evalf_def.h>
#include <gkyl_rect_grid.h>
#include <gkyl_wave_geom.h>

/**
 * Create a new wave geometry object. 
 *
 * @param grid Grid on which geometry lives
 * @param range Range on which geometry should be constructed
 * @param mapc2p Mapping from computational to physical space
 * @param ctx Context for use in mapping
 */
struct gkyl_wave_geom* gkyl_wave_geom_new(const struct gkyl_rect_grid *grid,
  struct gkyl_range *range, evalf_t mapc2p, void *ctx);

/**
 * Release geometry object.
 *
 * @param wg Wave geometry object to release.
 */
void gkyl_wave_geom_release(struct gkyl_wave_geom *wg);
