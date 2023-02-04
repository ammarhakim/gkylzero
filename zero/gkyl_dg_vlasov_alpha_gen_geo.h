#pragma once

#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_rect_grid.h>

/**
 * Compute alpha_gen_geo = v^i components for vlasov neutral
 * streaming term. 
 *
 * @param conf_basis Basis functions used in conf space expansions
 * @param phase_basis Basis functions used in phase space expansions
 * @param conf_rng Update range (Configuration space)
 * @param phase_rng Phase-space range
 * @param grid Phase-space grid for cell center and width
 * @param tv_comp Tangent vector components dX/dx, dX/dy, dX/dz, ...
 * @param gij Metric coefficients g^ij
 * @param alpha_geo Output contravariant velocity components v^i
 */
void
gkyl_dg_alpha_gen_geo(struct gkyl_basis *conf_basis,
  struct gkyl_basis *phase_basis,
  const struct gkyl_range *conf_rng, const struct gkyl_range *phase_rng,
  const struct gkyl_rect_grid *grid,const struct gkyl_array* tv_comp,
  const struct gkyl_array* gij, struct gkyl_array* alpha_geo);
