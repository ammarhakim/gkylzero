#include <assert.h>
#include <math.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_dg_updater_vlasov.h>
#include <gkyl_dg_updater_bflux_vlasov.h>
#include <gkyl_dg_updater_bflux_vlasov_priv.h>
#include <gkyl_util.h>

gkyl_dg_updater_bflux_vlasov*
gkyl_dg_updater_bflux_vlasov_new(const struct gkyl_rect_grid *grid, 
  int cdim, const gkyl_dg_updater_vlasov *vlasov, bool use_gpu)
{
  printf("Gkeyll was here!\n");
  gkyl_dg_updater_bflux_vlasov *up = gkyl_malloc(sizeof(gkyl_dg_updater_bflux_vlasov));
  up->slvr = gkyl_ghost_surf_calc_new(grid, gkyl_dg_updater_vlasov_acquire_eqn(vlasov), cdim, use_gpu);
  up->bflux_tm = 0.0;
  
  return up;
}

void
gkyl_dg_updater_bflux_vlasov_advance(gkyl_dg_updater_bflux_vlasov *bflux,
  const struct gkyl_range *update_rng,
  const struct gkyl_array* GKYL_RESTRICT fIn, struct gkyl_array* GKYL_RESTRICT rhs)
{
  printf("Advance the Gkeyll!\n");
  struct timespec wst = gkyl_wall_clock();
  gkyl_ghost_surf_calc_advance(bflux->slvr, update_rng, fIn, rhs);
  bflux->bflux_tm += gkyl_time_diff_now_sec(wst);
}

struct gkyl_dg_updater_bflux_vlasov_tm
OAgkyl_dg_updater_bflux_vlasov_get_tm(const gkyl_dg_updater_bflux_vlasov *bflux)
{
  return (struct gkyl_dg_updater_bflux_vlasov_tm) {
    .bflux_tm = bflux->bflux_tm,
  };
}

void
gkyl_dg_updater_bflux_vlasov_release(gkyl_dg_updater_bflux_vlasov* bflux)
{
  printf("Here, bitch!\n");
  gkyl_ghost_surf_calc_release(bflux->slvr);
  gkyl_free(bflux);
}

#ifdef GKYL_HAVE_CUDA

void
gkyl_dg_updater_bflux_vlasov_advance_cu(gkyl_dg_updater_bflux_vlasov *bflux,
  const struct gkyl_range *update_rng,
  const struct gkyl_array* GKYL_RESTRICT fIn, struct gkyl_array* GKYL_RESTRICT rhs)
{
  struct timespec wst = gkyl_wall_clock();
  gkyl_ghost_surf_calc_advance_cu(bflux->slvr, update_rng, fIn, rhs);
  bflux->bflux_tm += gkyl_time_diff_now_sec(wst);
}

#endif

#ifndef GKYL_HAVE_CUDA

void
gkyl_dg_updater_bflux_vlasov_advance_cu(gkyl_dg_updater_bflux_vlasov *bflux,
  const struct gkyl_range *update_rng,
  const struct gkyl_array* GKYL_RESTRICT fIn, struct gkyl_array* GKYL_RESTRICT rhs)
{
  assert(false);
}

#endif
