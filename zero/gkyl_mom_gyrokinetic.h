#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_range.h>
#include <gkyl_mom_type.h>

/**
 * Create new gyrokinetic moment type object. Valid 'mom' strings are "M0",
 * "M1", "M2", "M2par", "M2perp", "M3par", "M3perp", "ThreeMoments"
 *
 * @param cbasis Configuration-space basis-functions
 * @param pbasis Phase-space basis-functions
 * @param conf_range Configuration-space range
 * @param vel_range Velocity space range for use in indexing velocity mappping.
 * @param mass Mass of species
 * @param vmap Velocity space mapping.
 * @param gk_geom Geometry object
 * @param mom Name of moment to compute.
 * @param use_gpu bool to determine if on GPU
 */
struct gkyl_mom_type* 
gkyl_mom_gyrokinetic_new(const struct gkyl_basis* cbasis, 
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range,
  const struct gkyl_range* vel_range, double mass, const struct gkyl_array *vmap,
  const struct gk_geometry *gk_geom, const char *mom, bool use_gpu);

/**
 * Create new integrated gyrokinetic moment type object. Lab-frame
 * integrated moments (M0, M1 (M1par), M2par, and M2perp (for vdim>1)) are computed.
 *
 * @param cbasis Configuration-space basis-functions.
 * @param pbasis Phase-space basis-functions.
 * @param conf_range Configuration-space range.
 * @param vel_range Velocity-space range.
 * @param mass Mass of species.
 * @param vmap Velocity space mapping.
 * @param gk_geom Geometry object.
 * @param use_gpu bool to determine if on GPU.
 */
struct gkyl_mom_type* 
gkyl_int_mom_gyrokinetic_new(const struct gkyl_basis* cbasis, 
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range, 
  const struct gkyl_range* vel_range, double mass,
  const struct gkyl_array* vmap, const struct gk_geometry *gk_geom, bool use_gpu);
