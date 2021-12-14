#include <assert.h>
#include <math.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_dg_lbo_updater.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_util.h>

void
gkyl_dg_lbo_updater_advance(gkyl_dg_lbo_updater *lbo, struct gkyl_range update_rng,
  const struct gkyl_array *fIn, struct gkyl_array *cflrate, struct gkyl_array *rhs)
{
  gkyl_hyper_dg_advance(lbo->diff, update_rng, fIn, cflrate, rhs);
  gkyl_hyper_dg_advance(lbo->drag, update_rng, fIn, cflrate, rhs);
}

gkyl_dg_lbo_updater*
gkyl_dg_lbo_updater_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *basis,
  const int vdim, const struct gkyl_dg_eqn *drag, const struct gkyl_dg_eqn *diff)
{
  gkyl_dg_lbo_updater *lbo = gkyl_malloc(sizeof(gkyl_dg_lbo_updater));

  int num_up_dirs = vdim;
  int up_dirs[GKYL_MAX_DIM], zero_flux_flags[GKYL_MAX_DIM];
  for (int d=0; d<vdim; ++d) {
    up_dirs[d] = d + basis->ndim - vdim;
    zero_flux_flags[d] = 1;
  }
  
  lbo->diff = gkyl_hyper_dg_new(grid, basis, diff, num_up_dirs, up_dirs, zero_flux_flags, 1);
  lbo->drag = gkyl_hyper_dg_new(grid, basis, drag, num_up_dirs, up_dirs, zero_flux_flags, 1);
  return lbo;
}

void
gkyl_dg_lbo_updater_release(gkyl_dg_lbo_updater* lbo)
{
  gkyl_hyper_dg_release(lbo->drag);
  gkyl_hyper_dg_release(lbo->diff);
  gkyl_free(lbo);
}
