#include <assert.h>
#include <math.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_dg_lbo_vlasov_diff.h>
#include <gkyl_dg_lbo_vlasov_drag.h>
#include <gkyl_dg_updater_lbo_vlasov.h>
#include <gkyl_dg_updater_lbo_vlasov_priv.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_util.h>

gkyl_dg_updater_lbo_vlasov*
gkyl_dg_updater_lbo_vlasov_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *cbasis,
  const struct gkyl_basis *pbasis, const struct gkyl_range *conf_range, bool use_gpu)
{
  gkyl_dg_updater_lbo_vlasov *up = gkyl_malloc(sizeof(gkyl_dg_updater_lbo_vlasov));

  up->coll_drag = gkyl_dg_lbo_vlasov_drag_new(cbasis, pbasis, conf_range, grid, use_gpu);
  up->coll_diff = gkyl_dg_lbo_vlasov_diff_new(cbasis, pbasis, conf_range, grid, use_gpu);

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

  up->diff_tm = 0.0; up->drag_tm = 0.0;
  
  return up;
}

void
gkyl_dg_updater_lbo_vlasov_advance(gkyl_dg_updater_lbo_vlasov *lbo,
  const struct gkyl_range *update_rng,
  const struct gkyl_array *nu_sum, const struct gkyl_array *nu_prim_moms,
  const struct gkyl_array* GKYL_RESTRICT fIn,
  struct gkyl_array* GKYL_RESTRICT cflrate, struct gkyl_array* GKYL_RESTRICT rhs)
{
  // Set arrays needed
  gkyl_lbo_vlasov_drag_set_auxfields(lbo->coll_drag,
    (struct gkyl_dg_lbo_vlasov_drag_auxfields) { .nuSum = nu_sum, .nuPrimMomsSum = nu_prim_moms });
  gkyl_lbo_vlasov_diff_set_auxfields(lbo->coll_diff,
    (struct gkyl_dg_lbo_vlasov_diff_auxfields) { .nuSum = nu_sum, .nuPrimMomsSum = nu_prim_moms });

  struct timespec wst = gkyl_wall_clock();
  gkyl_hyper_dg_advance(lbo->diff, update_rng, fIn, cflrate, rhs);
  lbo->diff_tm += gkyl_time_diff_now_sec(wst);

  wst = gkyl_wall_clock();
  gkyl_hyper_dg_advance(lbo->drag, update_rng, fIn, cflrate, rhs);
  lbo->drag_tm += gkyl_time_diff_now_sec(wst);
}

struct gkyl_dg_updater_lbo_vlasov_tm
gkyl_dg_updater_lbo_vlasov_get_tm(const gkyl_dg_updater_lbo_vlasov *lbo)
{
  return (struct gkyl_dg_updater_lbo_vlasov_tm) {
    .diff_tm = lbo->diff_tm,
    .drag_tm = lbo->drag_tm
  };
}

void
gkyl_dg_updater_lbo_vlasov_release(gkyl_dg_updater_lbo_vlasov* lbo)
{
  gkyl_dg_eqn_release(lbo->coll_diff);
  gkyl_dg_eqn_release(lbo->coll_drag);
  gkyl_hyper_dg_release(lbo->drag);
  gkyl_hyper_dg_release(lbo->diff);
  gkyl_free(lbo);
}

#ifdef GKYL_HAVE_CUDA

void
gkyl_dg_updater_lbo_vlasov_advance_cu(gkyl_dg_updater_lbo_vlasov *lbo,
  const struct gkyl_range *update_rng,
  const struct gkyl_array *nu_sum, const struct gkyl_array *nu_prim_moms,
  const struct gkyl_array* GKYL_RESTRICT fIn, struct gkyl_array* GKYL_RESTRICT cflrate,
  struct gkyl_array* GKYL_RESTRICT rhs)
{
  // Set arrays needed
  gkyl_lbo_vlasov_drag_set_auxfields(lbo->coll_drag,
    (struct gkyl_dg_lbo_vlasov_drag_auxfields) { .nuSum = nu_sum, .nuPrimMomsSum = nu_prim_moms });
  gkyl_lbo_vlasov_diff_set_auxfields(lbo->coll_diff,
    (struct gkyl_dg_lbo_vlasov_diff_auxfields) { .nuSum = nu_sum, .nuPrimMomsSum = nu_prim_moms });

  struct timespec wst = gkyl_wall_clock();
  gkyl_hyper_dg_advance_cu(lbo->diff, update_rng, fIn, cflrate, rhs);
  lbo->diff_tm += gkyl_time_diff_now_sec(wst);

  wst = gkyl_wall_clock();
  gkyl_hyper_dg_advance_cu(lbo->drag, update_rng, fIn, cflrate, rhs);
  lbo->drag_tm += gkyl_time_diff_now_sec(wst);
}

#endif

#ifndef GKYL_HAVE_CUDA

void
gkyl_dg_updater_lbo_vlasov_advance_cu(gkyl_dg_updater_lbo_vlasov *lbo,
  const struct gkyl_range *update_rng,
  const struct gkyl_array *nu_sum, const struct gkyl_array *nu_prim_moms, 
  const struct gkyl_array *fIn, struct gkyl_array *cflrate, struct gkyl_array *rhs)
{
  assert(false);
}

#endif
