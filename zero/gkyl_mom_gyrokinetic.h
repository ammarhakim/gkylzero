#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_mom_type.h>

/**
 * Create new Gyrokinetic moment type object. Valid 'mom' strings are "GkM0",
 * "GkM1", "GkM2", "GkM2par", "GkM2perp", "GkM3par", "GkM3perp", "ThreeMoments"
 *
 * @param cbasis Configuration-space basis-functions
 * @param pbasis Phase-space basis-functions
 * @param conf_range Configuration-space range
 * @param mass Mass of species
 * @param mom Name of moment to compute.
 */
struct gkyl_mom_type* gkyl_mom_gyrokinetic_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* conf_range, double mass, const char *mom);

/**
 * Create new Gyrokinetic moment type object on NV-GPU: see new() method
 * above for documentation.
 */
struct gkyl_mom_type* gkyl_mom_gyrokinetic_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* conf_range, double mass, const char *mom);

/**
 * Set magnitude of the magnetic field, bmag, needed in computing moments.
 * 
 * @param momt Moment type pointer
 * @param bmag Pointer to magnitude of magnetic field
 */
void gkyl_gyrokinetic_set_bmag(const struct gkyl_mom_type *momt, const struct gkyl_array *bmag);

#ifdef GKYL_HAVE_CUDA
/**
 * CUDA device functions to set magnitude of the magnetic field, bmag, needed in computing moments.
 * 
 * @param momt Moment type pointer
 * @param bmag Pointer to magnitude of magnetic field
 */
void gkyl_gyrokinetic_set_bmag_cu(const struct gkyl_mom_type *momt, const struct gkyl_array *bmag);

#endif
