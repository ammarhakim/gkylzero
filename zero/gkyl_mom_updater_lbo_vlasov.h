#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_mom_updater_lbo_vlasov gkyl_mom_updater_lbo_vlasov;

gkyl_mom_updater_lbo_vlasov*
gkyl_mom_updater_lbo_vlasov_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *cbasis,
  const struct gkyl_basis *pbasis, const double *v_bounds, bool collides_with_fluid, bool use_gpu);

void
gkyl_mom_updater_lbo_vlasov_advance(gkyl_mom_updater_lbo_vlasov *up, struct gkyl_basis cbasis,
  struct gkyl_range *phase_rng, struct gkyl_range *conf_rng, 
  const struct gkyl_array *fin, const struct gkyl_array *moms, struct gkyl_array *boundary_corrections, 
  struct gkyl_array *uout, struct gkyl_array *vtSqout);

void
gkyl_mom_updater_lbo_vlasov_advance_cu(gkyl_mom_updater_lbo_vlasov *up, struct gkyl_basis cbasis,
  struct gkyl_range *phase_rng, struct gkyl_range *conf_rng, 
  const struct gkyl_array *fin, const struct gkyl_array *moms, struct gkyl_array *boundary_corrections, 
  struct gkyl_array *uout, struct gkyl_array *vtSqout);

void
gkyl_mom_updater_lbo_vlasov_release(gkyl_mom_updater_lbo_vlasov* up);

