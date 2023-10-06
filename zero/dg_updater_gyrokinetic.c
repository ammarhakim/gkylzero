#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_dg_gyrokinetic.h>
#include <gkyl_dg_updater_gyrokinetic.h>
#include <gkyl_dg_updater_gyrokinetic_priv.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_util.h>

struct gkyl_dg_eqn*
gkyl_dg_updater_gyrokinetic_acquire_eqn(const gkyl_dg_updater_gyrokinetic* gyrokinetic)
{
  return gkyl_dg_eqn_acquire(gyrokinetic->dgeqn);
}

struct gkyl_dg_updater_gyrokinetic*
gkyl_dg_updater_gyrokinetic_new(const struct gkyl_rect_grid *grid, 
  const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis, 
  const struct gkyl_range *conf_range, const bool *is_zero_flux_dir,
  enum gkyl_gkeqn_id eqn_id, double charge, double mass, bool use_gpu)
{
  struct gkyl_dg_updater_gyrokinetic *up = gkyl_malloc(sizeof(struct gkyl_dg_updater_gyrokinetic));

  up->use_gpu = use_gpu;

  up->eqn_id = eqn_id;
  if (eqn_id == GKYL_GK_DEFAULT)
    up->dgeqn = gkyl_dg_gyrokinetic_new(cbasis, pbasis, conf_range, charge, mass, use_gpu);
//  else if (eqn_id == GKYL_GK_SOFTMIRROR)
//    up->dgeqn = gkyl_dg_gyrokinetic_softmirror_new(cbasis, pbasis, conf_range, charge, mass, use_gpu);

  int cdim = cbasis->ndim, pdim = pbasis->ndim;
  int vdim = pdim-cdim;
  int up_dirs[GKYL_MAX_DIM] = {0};
  int num_up_dirs = cdim+1;
  for (int d=0; d<num_up_dirs; ++d) up_dirs[d] = d;

  int zero_flux_flags[GKYL_MAX_DIM] = {0};
  for (int d=0; d<cdim; ++d)
    zero_flux_flags[d] = is_zero_flux_dir[d]? 1 : 0;
  for (int d=cdim; d<pdim; ++d)
    zero_flux_flags[d] = 1; // zero-flux BCs in vel-space

  up->hyperdg = gkyl_hyper_dg_new(grid, pbasis, up->dgeqn,
    num_up_dirs, up_dirs, zero_flux_flags, 1, use_gpu);

  return up;
}

void
gkyl_dg_updater_gyrokinetic_advance(struct gkyl_dg_updater_gyrokinetic *up,
  const struct gkyl_range *update_rng,
  const struct gkyl_array *bmag, const struct gkyl_array *jacobtot_inv,
  const struct gkyl_array *cmag, const struct gkyl_array *b_i,
  const struct gkyl_array *phi, const struct gkyl_array *apar,
  const struct gkyl_array *apardot, const struct gkyl_array* GKYL_RESTRICT fIn,
  struct gkyl_array* GKYL_RESTRICT cflrate, struct gkyl_array* GKYL_RESTRICT rhs)
{
  // Set arrays needed.
  if (up->eqn_id == GKYL_GK_DEFAULT)
    gkyl_gyrokinetic_set_auxfields(up->dgeqn, 
      (struct gkyl_dg_gyrokinetic_auxfields) { .bmag = bmag, .jacobtot_inv = jacobtot_inv,
        .cmag = cmag, .b_i = b_i, .phi = phi, .apar = apar, .apardot = apardot });

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    gkyl_hyper_dg_advance_cu(up->hyperdg, update_rng, fIn, cflrate, rhs);
  else
    gkyl_hyper_dg_advance(up->hyperdg, update_rng, fIn, cflrate, rhs);
#else
  gkyl_hyper_dg_advance(up->hyperdg, update_rng, fIn, cflrate, rhs);
#endif
}

void
gkyl_dg_updater_gyrokinetic_release(struct gkyl_dg_updater_gyrokinetic* up)
{
  gkyl_dg_eqn_release(up->dgeqn);
  gkyl_hyper_dg_release(up->hyperdg);
  gkyl_free(up);
}
