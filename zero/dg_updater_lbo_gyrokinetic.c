#include <assert.h>
#include <math.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_dg_lbo_gyrokinetic_drag.h>
#include <gkyl_dg_lbo_gyrokinetic_diff.h>
#include <gkyl_dg_updater_lbo_gyrokinetic.h>
#include <gkyl_dg_updater_collisions_priv.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_util.h>

struct gkyl_dg_updater_collisions*
gkyl_dg_updater_lbo_gyrokinetic_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *cbasis,
  const struct gkyl_basis *pbasis, const struct gkyl_range *conf_range, double mass, bool use_gpu)
{
  struct gkyl_dg_updater_collisions *up = gkyl_malloc(sizeof(gkyl_dg_updater_collisions));

  up->coll_drag = gkyl_dg_lbo_gyrokinetic_drag_new(cbasis, pbasis, conf_range, mass, use_gpu);
  up->coll_diff = gkyl_dg_lbo_gyrokinetic_diff_new(cbasis, pbasis, conf_range, mass, use_gpu);

  int cdim = cbasis->ndim, pdim = pbasis->ndim;
  int vdim = pdim-cdim;
  int num_up_dirs = vdim;
  int up_dirs[GKYL_MAX_DIM] = { 0 };
  for (int d=0; d<vdim; ++d)
    up_dirs[d] = d + pbasis->ndim - vdim;

  int zero_flux_flags[GKYL_MAX_DIM] = { 0 };
  for (int d=cdim; d<pdim; ++d)
    zero_flux_flags[d] = 1;
  
  up->diff = gkyl_hyper_dg_new(grid, pbasis, up->coll_diff, num_up_dirs, up_dirs, zero_flux_flags, 1, use_gpu);
  up->drag = gkyl_hyper_dg_new(grid, pbasis, up->coll_drag, num_up_dirs, up_dirs, zero_flux_flags, 1, use_gpu);
  
  return up;
}

void
gkyl_dg_updater_lbo_gyrokinetic_advance(struct gkyl_dg_updater_collisions *lbo,
  const struct gkyl_range *update_rng,
  const struct gkyl_array *bmag_inv,
  const struct gkyl_array *nu_sum, const struct gkyl_array *nu_u, const struct gkyl_array *nu_vthsq,
  const struct gkyl_array* GKYL_RESTRICT fIn,
  struct gkyl_array* GKYL_RESTRICT cflrate, struct gkyl_array* GKYL_RESTRICT rhs)
{
  // Set arrays needed
  gkyl_lbo_gyrokinetic_drag_set_auxfields(lbo->coll_drag,
    (struct gkyl_dg_lbo_gyrokinetic_drag_auxfields) { .bmag_inv = bmag_inv, .nuSum = nu_sum, .nuUSum = nu_u, .nuVtSqSum = nu_vthsq });
  gkyl_lbo_gyrokinetic_diff_set_auxfields(lbo->coll_diff,
    (struct gkyl_dg_lbo_gyrokinetic_diff_auxfields) { .bmag_inv = bmag_inv, .nuSum = nu_sum, .nuUSum = nu_u, .nuVtSqSum = nu_vthsq });
  
  gkyl_hyper_dg_advance(lbo->diff, update_rng, fIn, cflrate, rhs);
  gkyl_hyper_dg_advance(lbo->drag, update_rng, fIn, cflrate, rhs);
}

struct gkyl_dg_updater_lbo_gyrokinetic_tm
gkyl_dg_updater_lbo_gyrokinetic_get_tm(const struct gkyl_dg_updater_collisions *coll)
{
  return (struct gkyl_dg_updater_lbo_gyrokinetic_tm) {
    .diff_tm = coll->diff_tm,
    .drag_tm = coll->drag_tm
  };
}

void
gkyl_dg_updater_lbo_gyrokinetic_release(struct gkyl_dg_updater_collisions* coll)
{
  gkyl_dg_eqn_release(coll->coll_diff);
  gkyl_dg_eqn_release(coll->coll_drag);
  gkyl_hyper_dg_release(coll->drag);
  gkyl_hyper_dg_release(coll->diff);
  gkyl_free(coll);
}

#ifdef GKYL_HAVE_CUDA

void
gkyl_dg_updater_lbo_gyrokinetic_advance_cu(struct gkyl_dg_updater_collisions *lbo,
  const struct gkyl_range *update_rng,
  const struct gkyl_array *bmag_inv,
  const struct gkyl_array *nu_sum, const struct gkyl_array *nu_u, const struct gkyl_array *nu_vthsq,
  const struct gkyl_array* GKYL_RESTRICT fIn, struct gkyl_array* GKYL_RESTRICT cflrate,
  struct gkyl_array* GKYL_RESTRICT rhs)
{
  // Set arrays needed
  gkyl_lbo_gyrokinetic_drag_set_auxfields(lbo->coll_drag,
    (struct gkyl_dg_lbo_gyrokinetic_drag_auxfields) { .bmag_inv = bmag_inv, .nuSum = nu_sum, .nuUSum = nu_u, .nuVtSqSum = nu_vthsq });
  gkyl_lbo_gyrokinetic_diff_set_auxfields(lbo->coll_diff,
    (struct gkyl_dg_lbo_gyrokinetic_diff_auxfields) { .bmag_inv = bmag_inv, .nuSum = nu_sum, .nuUSum = nu_u, .nuVtSqSum = nu_vthsq });

  gkyl_hyper_dg_advance_cu(lbo->diff, update_rng, fIn, cflrate, rhs);
  gkyl_hyper_dg_advance_cu(lbo->drag, update_rng, fIn, cflrate, rhs);
}

#endif

#ifndef GKYL_HAVE_CUDA

void
gkyl_dg_updater_lbo_gyrokinetic_advance_cu(struct gkyl_dg_updater_collisions *lbo,
  const struct gkyl_range *update_rng,
  const struct gkyl_array *bmag_inv,
  const struct gkyl_array *nu_sum, const struct gkyl_array *nu_u, const struct gkyl_array *nu_vthsq, 
  const struct gkyl_array *fIn, struct gkyl_array *cflrate, struct gkyl_array *rhs)
{
  assert(false);
}

#endif
