#pragma once

#include <gkyl_gk_geometry.h>
/**
 * Create a new geometry object using tokamak input (efit)
 *
 * @param geometry_inp geometry input struct containing grid, range, and other geo info
 */
struct gk_geometry* gkyl_gk_geometry_tok_new(struct gkyl_gk_geometry_inp *geometry_inp);



