#include <assert.h>
#include <math.h>
#include <time.h>
#include <gkyl_alloc.h>
#include <gkyl_dg_updater_gyrokinetic.h>
#include <gkyl_dg_updater_bflux_gyrokinetic.h>
#include <gkyl_dg_updater_bflux_gyrokinetic_priv.h>
#include <gkyl_util.h>

struct gkyl_dg_updater_bflux_gyrokinetic*
gkyl_dg_updater_bflux_gyrokinetic_new(const struct gkyl_rect_grid *grid, 
  int cdim, const gkyl_dg_updater_gyrokinetic *gyrokinetic, bool use_gpu)
{
  struct gkyl_dg_updater_bflux_gyrokinetic *up = gkyl_malloc(sizeof(struct gkyl_dg_updater_bflux_gyrokinetic));
  up->slvr = gkyl_ghost_surf_calc_new(grid, gkyl_dg_updater_gyrokinetic_acquire_eqn(gyrokinetic), cdim, use_gpu);
  up->bflux_tm = 0.0;
  
  return up;
}

void
gkyl_dg_updater_bflux_gyrokinetic_advance(struct gkyl_dg_updater_bflux_gyrokinetic *up,
  const struct gkyl_range *update_rng,
  const struct gkyl_array* GKYL_RESTRICT fIn, struct gkyl_array* GKYL_RESTRICT rhs)
{
  struct timespec wst = gkyl_wall_clock();
  gkyl_ghost_surf_calc_advance(up->slvr, update_rng, fIn, rhs);
  up->bflux_tm += gkyl_time_diff_now_sec(wst);
}

struct gkyl_dg_updater_bflux_gyrokinetic_tm
OAgkyl_dg_updater_bflux_gyrokinetic_get_tm(const struct gkyl_dg_updater_bflux_gyrokinetic *up)
{
  return (struct gkyl_dg_updater_bflux_gyrokinetic_tm) {
    .bflux_tm = up->bflux_tm,
  };
}

void
gkyl_dg_updater_bflux_gyrokinetic_release(struct gkyl_dg_updater_bflux_gyrokinetic* up)
{
  gkyl_ghost_surf_calc_release(up->slvr);
  gkyl_free(up);
}

#ifdef GKYL_HAVE_CUDA

void
gkyl_dg_updater_bflux_gyrokinetic_advance_cu(struct gkyl_dg_updater_bflux_gyrokinetic *up,
  const struct gkyl_range *update_rng,
  const struct gkyl_array* GKYL_RESTRICT fIn, struct gkyl_array* GKYL_RESTRICT rhs)
{
  struct timespec wst = gkyl_wall_clock();
  gkyl_ghost_surf_calc_advance_cu(up->slvr, update_rng, fIn, rhs);
  up->bflux_tm += gkyl_time_diff_now_sec(wst);
}

#endif

#ifndef GKYL_HAVE_CUDA

void
gkyl_dg_updater_bflux_gyrokinetic_advance_cu(struct gkyl_dg_updater_bflux_gyrokinetic *up,
  const struct gkyl_range *update_rng,
  const struct gkyl_array* GKYL_RESTRICT fIn, struct gkyl_array* GKYL_RESTRICT rhs)
{
  assert(false);
}

#endif
