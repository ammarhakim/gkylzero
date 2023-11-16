#include <assert.h>
#include <math.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_dg_rad_gyrokinetic_drag.h>
#include <gkyl_dg_updater_rad_gyrokinetic.h>
#include <gkyl_dg_updater_collisions_priv.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_util.h>

// Done first pass
struct gkyl_dg_updater_collisions*
gkyl_dg_updater_rad_gyrokinetic_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *cbasis,
				    const struct gkyl_basis *pbasis, const struct gkyl_range *conf_range, const struct gkyl_array *bmag, const struct gkyl_array *fit_params, bool use_gpu)
{
  printf("Before allocation\n");
  struct gkyl_dg_updater_collisions *up = gkyl_malloc(sizeof(gkyl_dg_updater_collisions));
  printf("Before coll creation\n");
  up->coll_drag = gkyl_dg_rad_gyrokinetic_drag_new(cbasis, pbasis, conf_range, grid, bmag, fit_params, use_gpu);
  printf("After coll creation\n");
  int cdim = cbasis->ndim, pdim = pbasis->ndim;
  int vdim = pdim-cdim;
  int num_up_dirs = vdim;
  int up_dirs[GKYL_MAX_DIM] = { 0 };
  for (int d=0; d<vdim; ++d)
    up_dirs[d] = d + pbasis->ndim - vdim;

  int zero_flux_flags[GKYL_MAX_DIM] = { 0 };
  for (int d=cdim; d<pdim; ++d)
    zero_flux_flags[d] = 1;
  printf("Before hyper dg new\n");
  up->drag = gkyl_hyper_dg_new(grid, pbasis, up->coll_drag, num_up_dirs, up_dirs, zero_flux_flags, 1, use_gpu);
  
  return up;
}

void
gkyl_dg_updater_rad_gyrokinetic_advance(struct gkyl_dg_updater_collisions *rad,
  const struct gkyl_range *update_rng,
  const struct gkyl_array *nI,
  const struct gkyl_array* GKYL_RESTRICT fIn,
  struct gkyl_array* GKYL_RESTRICT cflrate, struct gkyl_array* GKYL_RESTRICT rhs)
{
  printf("In updater rad advance\n");
  // Set arrays needed
  gkyl_rad_gyrokinetic_drag_set_auxfields(rad->coll_drag,
    (struct gkyl_dg_rad_gyrokinetic_drag_auxfields) { .nI = nI });
  printf("Aux fields set(updater-advance)\n");
  gkyl_hyper_dg_advance(rad->drag, update_rng, fIn, cflrate, rhs);
  printf("After hyper advance (updater-advance)\n");
}

struct gkyl_dg_updater_rad_gyrokinetic_tm
gkyl_dg_updater_rad_gyrokinetic_get_tm(const struct gkyl_dg_updater_collisions *coll)
{
  return (struct gkyl_dg_updater_rad_gyrokinetic_tm) {
    .drag_tm = coll->drag_tm
  };
}

void
gkyl_dg_updater_rad_gyrokinetic_release(struct gkyl_dg_updater_collisions* coll)
{
  gkyl_dg_eqn_release(coll->coll_drag);
  gkyl_hyper_dg_release(coll->drag);
  gkyl_free(coll);
}

#ifdef GKYL_HAVE_CUDA

void
gkyl_dg_updater_rad_gyrokinetic_advance_cu(struct gkyl_dg_updater_collisions *rad,
  const struct gkyl_range *update_rng,
  const struct gkyl_array *nI, 
  const struct gkyl_array* GKYL_RESTRICT fIn, struct gkyl_array* GKYL_RESTRICT cflrate,
  struct gkyl_array* GKYL_RESTRICT rhs)
{
  // Set arrays needed
  gkyl_rad_gyrokinetic_drag_set_auxfields(rad->coll_drag,
    (struct gkyl_dg_rad_gyrokinetic_drag_auxfields) { .nI = nI});

  gkyl_hyper_dg_advance_cu(rad->drag, update_rng, fIn, cflrate, rhs);
}

#endif

#ifndef GKYL_HAVE_CUDA

void
gkyl_dg_updater_rad_gyrokinetic_advance_cu(struct gkyl_dg_updater_collisions *rad,
  const struct gkyl_range *update_rng,
					   const struct gkyl_array *nI,// const struct gkyl_array *vnu, 
					   //const struct gkyl_array *vsqnu, 
  const struct gkyl_array *fIn, struct gkyl_array *cflrate, struct gkyl_array *rhs)
{
  assert(false);
  }

#endif
