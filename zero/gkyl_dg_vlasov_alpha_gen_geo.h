#pragma once

#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_rect_grid.h>

/**
 * Compute alpha_gen_geo = v^i components for vlasov neutral
 * streaming term. 
 *
 * @param basis Basis functions used in expansions
 * @param conf_rng Update range (Configuration space)
 * @param phase_rng Phase-space range
 * @param grid Phase-space grid for cell center and width
 * @param tv_comp Tangent vector components dX/dx, dX/dy, dX/dz, ...
 * @param gxx Metric coefficient g^xx
 * @param gxy Metric coefficient g^xy
 * @param gxz Metric coefficient g^xz
 * @param gyy Metric coefficient g^yy
 * @param gyz Metric coefficient g^yz
 * @param gzz Metric coefficient g^zz
 * @param alpha_geo Output contravariant velocity components v^i
 */
void gkyl_dg_alpha_gen_geo(struct gkyl_basis basis, const struct gkyl_range *conf_rng,
  const struct gkyl_range *phase_rng, const struct gkyl_array *tv_comp,
  const struct gkyl_rect_grid *grid,
  const struct gkyl_array *gxx, const struct gkyl_array *gxy,
  const struct gkyl_array *gxz, const struct gkyl_array *gyy,
  const struct gkyl_array *gyz, const struct gkyl_array *gzz,
  struct gkyl_array *alpha_geo);
