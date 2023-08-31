#include <assert.h>
#include <math.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_dg_diffusion.h>
#include <gkyl_dg_diffusion_gen.h>
#include <gkyl_dg_updater_diffusion.h>
#include <gkyl_dg_updater_diffusion_priv.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_util.h>

struct gkyl_dg_eqn*
gkyl_dg_updater_diffusion_acquire_eqn(const struct gkyl_dg_updater_diffusion *up)
{
  return gkyl_dg_eqn_acquire(up->dgeqn);
}

struct gkyl_dg_updater_diffusion*
gkyl_dg_updater_diffusion_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *basis, const struct gkyl_basis *cbasis,
  enum gkyl_diffusion_id diffusion_id, bool *diff_in_dir, int diff_order,
  const struct gkyl_range *conf_range, bool use_gpu)
{
  struct gkyl_dg_updater_diffusion *up = gkyl_malloc(sizeof(struct gkyl_dg_updater_diffusion));

  up->diffid = diffusion_id;

  if (diffusion_id == GKYL_DIFFUSION_GEN)
    up->dgeqn = gkyl_dg_diffusion_gen_new(cbasis, conf_range, use_gpu);
  else 
    up->dgeqn = gkyl_dg_diffusion_new(basis, cbasis, diffusion_id, diff_in_dir,
                                      diff_order, conf_range, use_gpu);

  int cdim = cbasis->ndim;
  int num_up_dirs = cdim;
  int up_dirs[GKYL_MAX_DIM], zero_flux_flags[GKYL_MAX_DIM];
  for (int d=0; d<cdim; ++d) {
    if (diff_in_dir) {
      if (diff_in_dir[d]) up_dirs[d] = d;
    } else {
      up_dirs[d] = d;
    }
    zero_flux_flags[d] = 0;
  }

  up->hyperdg = gkyl_hyper_dg_new(grid, cbasis, up->dgeqn, num_up_dirs, up_dirs, zero_flux_flags, 1, use_gpu);

  up->diffusion_tm = 0.0;

  return up;
}

void
gkyl_dg_updater_diffusion_advance(struct gkyl_dg_updater_diffusion *up,
  const struct gkyl_range *update_rng, const struct gkyl_array *coeff,
  const struct gkyl_array* GKYL_RESTRICT fIn, struct gkyl_array* GKYL_RESTRICT cflrate,
  struct gkyl_array* GKYL_RESTRICT rhs)
{
  struct timespec wst = gkyl_wall_clock();
  // Set arrays needed and call the specific advance method required
  if (up->diffid == GKYL_DIFFUSION_GEN) {
    gkyl_diffusion_gen_set_auxfields(up->dgeqn, (struct gkyl_dg_diffusion_gen_auxfields) { .Dij = coeff });
#ifdef GKYL_HAVE_CUDA
//    if (up->use_gpu)
//      // hyper_dg_gen_stencil NOT YET IMPLEMENTED ON DEVICE
//      gkyl_hyper_dg_gen_stencil_advance_cu(up->hyperdg, update_rng, fIn, cflrate, rhs);
//    else
//      gkyl_hyper_dg_gen_stencil_advance(up->hyperdg, update_rng, fIn, cflrate, rhs);
    assert(false);
#else
    gkyl_hyper_dg_gen_stencil_advance(up->hyperdg, update_rng, fIn, cflrate, rhs);
#endif
  } else {
    gkyl_diffusion_set_auxfields(up->dgeqn, (struct gkyl_dg_diffusion_auxfields) { .D = coeff });
#ifdef GKYL_HAVE_CUDA
    if (up->use_gpu)
      gkyl_hyper_dg_advance_cu(up->hyperdg, update_rng, fIn, cflrate, rhs);
    else
      gkyl_hyper_dg_advance(up->hyperdg, update_rng, fIn, cflrate, rhs);
#else
    gkyl_hyper_dg_advance(up->hyperdg, update_rng, fIn, cflrate, rhs);
#endif
  }
  up->diffusion_tm += gkyl_time_diff_now_sec(wst);
}

struct gkyl_dg_updater_diffusion_tm
gkyl_dg_updater_diffusion_get_tm(const struct gkyl_dg_updater_diffusion *up)
{
  return (struct gkyl_dg_updater_diffusion_tm) {
    .diffusion_tm = up->diffusion_tm,
  };
}

void
gkyl_dg_updater_diffusion_release(struct gkyl_dg_updater_diffusion *up)
{
  gkyl_dg_eqn_release(up->dgeqn);
  gkyl_hyper_dg_release(up->hyperdg);
  gkyl_free(up);
}
