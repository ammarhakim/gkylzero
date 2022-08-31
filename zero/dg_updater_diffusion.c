#include <assert.h>
#include <math.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_dg_diffusion.h>
#include <gkyl_dg_diffusion_euler_iso.h>
#include <gkyl_dg_gen_diffusion.h>
#include <gkyl_dg_updater_diffusion.h>
#include <gkyl_dg_updater_diffusion_priv.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_util.h>

struct gkyl_dg_eqn*
gkyl_dg_updater_diffusion_acquire_eqn(const gkyl_dg_updater_diffusion* diffusion)
{
  return gkyl_dg_eqn_acquire(diffusion->eqn_diffusion);
}

gkyl_dg_updater_diffusion*
gkyl_dg_updater_diffusion_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *cbasis, const struct gkyl_range *conf_range,
  enum gkyl_diffusion_id diffusion_id, bool use_gpu)
{
  gkyl_dg_updater_diffusion *up = gkyl_malloc(sizeof(gkyl_dg_updater_diffusion));

  if (diffusion_id == GKYL_ISO_DIFFUSION)
    up->eqn_diffusion = gkyl_dg_diffusion_new(cbasis, conf_range, use_gpu);
  else if (diffusion_id == GKYL_ANISO_DIFFUSION)
    up->eqn_diffusion = gkyl_dg_gen_diffusion_new(cbasis, conf_range, use_gpu);
  else if (diffusion_id == GKYL_EULER_ISO_DIFFUSION)
    up->eqn_diffusion = gkyl_dg_diffusion_e_i_new(cbasis, conf_range, use_gpu);

  int cdim = cbasis->ndim;
  int up_dirs[GKYL_MAX_DIM], zero_flux_flags[GKYL_MAX_DIM];
  for (int d=0; d<cdim; ++d) {
    up_dirs[d] = d;
    zero_flux_flags[d] = 0;
  }
  int num_up_dirs = cdim;

  up->up_diffusion = gkyl_hyper_dg_new(grid, cbasis, up->eqn_diffusion, num_up_dirs, up_dirs, zero_flux_flags, 1, use_gpu);

  up->diffusion_tm = 0.0;

  return up;
}

void
gkyl_dg_updater_diffusion_advance(gkyl_dg_updater_diffusion *diffusion,
  enum gkyl_diffusion_id diffusion_id, const struct gkyl_range *update_rng,
  const struct gkyl_array *D,
  const struct gkyl_array *u,
  const struct gkyl_array* GKYL_RESTRICT fIn,
  struct gkyl_array* GKYL_RESTRICT cflrate, struct gkyl_array* GKYL_RESTRICT rhs)
{
  struct timespec wst = gkyl_wall_clock();
  // Set arrays needed and call the specific advance method required
  if (diffusion_id == GKYL_ISO_DIFFUSION) {
    gkyl_diffusion_set_auxfields(diffusion->eqn_diffusion,
      (struct gkyl_dg_diffusion_auxfields) { .D = D });
    gkyl_hyper_dg_advance(diffusion->up_diffusion, update_rng, fIn, cflrate, rhs);
  }
  else if (diffusion_id == GKYL_ANISO_DIFFUSION) {
    gkyl_gen_diffusion_set_auxfields(diffusion->eqn_diffusion,
      (struct gkyl_dg_gen_diffusion_auxfields) { .Dij = D });
    gkyl_hyper_dg_gen_stencil_advance(diffusion->up_diffusion, update_rng, fIn, cflrate, rhs);
  }
  else if (diffusion_id == GKYL_EULER_ISO_DIFFUSION) {
    gkyl_diffusion_e_i_set_auxfields(diffusion->eqn_diffusion,
      (struct gkyl_dg_diffusion_e_i_auxfields) { .D = D, .u_i = u });
    gkyl_hyper_dg_advance(diffusion->up_diffusion, update_rng, fIn, cflrate, rhs);
  }
  diffusion->diffusion_tm += gkyl_time_diff_now_sec(wst);
}

struct gkyl_dg_updater_diffusion_tm
gkyl_dg_updater_diffusion_get_tm(const gkyl_dg_updater_diffusion *diffusion)
{
  return (struct gkyl_dg_updater_diffusion_tm) {
    .diffusion_tm = diffusion->diffusion_tm,
  };
}

void
gkyl_dg_updater_diffusion_release(gkyl_dg_updater_diffusion* diffusion)
{
  gkyl_dg_eqn_release(diffusion->eqn_diffusion);
  gkyl_hyper_dg_release(diffusion->up_diffusion);
  gkyl_free(diffusion);
}

#ifdef GKYL_HAVE_CUDA

void
gkyl_dg_updater_diffusion_advance_cu(gkyl_dg_updater_diffusion *diffusion,
  enum gkyl_diffusion_id diffusion_id, const struct gkyl_range *update_rng,
  const struct gkyl_array *D,
  const struct gkyl_array *u,
  const struct gkyl_array* GKYL_RESTRICT fIn,
  struct gkyl_array* GKYL_RESTRICT cflrate, struct gkyl_array* GKYL_RESTRICT rhs)
{
  struct timespec wst = gkyl_wall_clock();
  // Set arrays needed and call the specific advance method required
  if (diffusion_id == GKYL_ISO_DIFFUSION) {
    gkyl_diffusion_set_auxfields(diffusion->eqn_diffusion,
      (struct gkyl_dg_diffusion_auxfields) { .D = D });
    gkyl_hyper_dg_advance_cu(diffusion->up_diffusion, update_rng, fIn, cflrate, rhs);
  }
  else if (diffusion_id == GKYL_ANISO_DIFFUSION) {
    gkyl_gen_diffusion_set_auxfields(diffusion->eqn_diffusion,
      (struct gkyl_dg_gen_diffusion_auxfields) { .Dij = D });
    gkyl_hyper_dg_gen_stencil_advance_cu(diffusion->up_diffusion, update_rng, fIn, cflrate, rhs);
  }
  else if (diffusion_id == GKYL_EULER_ISO_DIFFUSION) {
    gkyl_diffusion_e_i_set_auxfields(diffusion->eqn_diffusion,
      (struct gkyl_dg_diffusion_e_i_auxfields) { .D = D, .u_i = u });//TODO: use consistent name for 'u'/'u_i'
    gkyl_hyper_dg_advance(diffusion->up_diffusion, update_rng, fIn, cflrate, rhs);
  }
  diffusion->diffusion_tm += gkyl_time_diff_now_sec(wst);
}

#endif

#ifndef GKYL_HAVE_CUDA

void
gkyl_dg_updater_diffusion_advance_cu(gkyl_dg_updater_diffusion *diffusion,
  enum gkyl_diffusion_id diffusion_id, const struct gkyl_range *update_rng,
  const struct gkyl_array *D,
  const struct gkyl_array *u,
  const struct gkyl_array* GKYL_RESTRICT fIn,
  struct gkyl_array* GKYL_RESTRICT cflrate, struct gkyl_array* GKYL_RESTRICT rhs)
{
  assert(false);
}

#endif
