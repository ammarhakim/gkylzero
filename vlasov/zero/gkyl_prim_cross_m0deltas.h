#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>

// Object type
typedef struct gkyl_prim_cross_m0deltas gkyl_prim_cross_m0deltas;

/**
 * Create a new updater that computes the m0_s*delta_s*(beta+1) prefactor
 * in the calculation of cross-primitive moments for LBO and BGK
 * collisions. That is:
 *   m0_s*delta_s*(beta+1) = m0_s*2*m_r*m0_r*nu_rs*(beta+1)/(m_s*m0_s*nu_sr+m_r*m0_r*nu_rs)
 *
 * @param basis Basis object (configuration space).
 * @param range Range in which we'll compute m0_s*delta_s.
 * @param betap1 beta+1 parameter in Greene's formulism.
 * @param use_gpu boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */
gkyl_prim_cross_m0deltas* gkyl_prim_cross_m0deltas_new(
  const struct gkyl_basis *basis, const struct gkyl_range *range,
  double betap1, bool use_gpu);

/**
 * Compute m0_s*delta_s.
 *
 * @param up Struct defining this updater..
 * @param basis Basis object (configuration space).
 * @param massself Mass of this species, m_s.
 * @param m0self Number density of this species, m0_s.
 * @param nuself Cross collision frequency of this species, nu_sr.
 * @param massother Mass of the other species, m_r.
 * @param m0other Number density of the other species, m0_r.
 * @param nuother Cross collision frequency of the other species, nu_rs.
 * @param prem0s Prefactor, M_0s for LBO, 1 for BGK.
 * @param range Range in which we'll compute m0_s*delta_s.
 * @param out Output array.
 * @return New updater pointer.
 */
void gkyl_prim_cross_m0deltas_advance(gkyl_prim_cross_m0deltas *up, struct gkyl_basis basis,
  double massself, const struct gkyl_array* m0self, const struct gkyl_array* nuself,
  double massother, const struct gkyl_array* m0other, const struct gkyl_array* nuother,
  const struct gkyl_array* prem0s, const struct gkyl_range *range, struct gkyl_array* out);

/**
 * Delete updater.
 *
 * @param pob Updater to delete.
 */
void gkyl_prim_cross_m0deltas_release(gkyl_prim_cross_m0deltas* up);
