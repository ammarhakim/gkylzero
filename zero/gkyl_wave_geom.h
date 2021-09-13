#pragma once

// Object type of geometry
typedef struct gkyl_wave_geom gkyl_wave_geom;

/**
 * Set index into the geometry object. This method must be called
 * before any geometry information is actuall accessed.
 */
void gkyl_wave_geom_set_idx(const gkyl_wave_geom *wg, const int *idx);
