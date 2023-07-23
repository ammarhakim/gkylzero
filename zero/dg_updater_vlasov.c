#include <assert.h>
#include <math.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_dg_vlasov.h>
#include <gkyl_dg_vlasov_pkpm.h>
#include <gkyl_dg_vlasov_sr.h>
#include <gkyl_dg_updater_vlasov.h>
#include <gkyl_dg_updater_vlasov_priv.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_util.h>

struct gkyl_dg_eqn*
gkyl_dg_updater_vlasov_acquire_eqn(const gkyl_dg_updater_vlasov* vlasov)
{
  return gkyl_dg_eqn_acquire(vlasov->eqn_vlasov);
}

gkyl_dg_updater_vlasov*
gkyl_dg_updater_vlasov_new(const struct gkyl_rect_grid *grid, 
  const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis, 
  const struct gkyl_range *conf_range, const struct gkyl_range *vel_range, const struct gkyl_range *phase_range,
  enum gkyl_model_id model_id, enum gkyl_field_id field_id, bool use_gpu)
{
  gkyl_dg_updater_vlasov *up = gkyl_malloc(sizeof(gkyl_dg_updater_vlasov));
  up->model_id = model_id;
  up->field_id = field_id;
  if (up->model_id == GKYL_MODEL_SR)
    up->eqn_vlasov = gkyl_dg_vlasov_sr_new(cbasis, pbasis, conf_range, vel_range, up->field_id, use_gpu);
  else if (up->model_id == GKYL_MODEL_PKPM)
    up->eqn_vlasov = gkyl_dg_vlasov_pkpm_new(cbasis, pbasis, conf_range, phase_range, use_gpu);
  else
    up->eqn_vlasov = gkyl_dg_vlasov_new(cbasis, pbasis, conf_range, phase_range, up->model_id, up->field_id, use_gpu);

  int cdim = cbasis->ndim, pdim = pbasis->ndim;
  int vdim = pdim-cdim;
  int up_dirs[GKYL_MAX_DIM], zero_flux_flags[GKYL_MAX_DIM];
  for (int d=0; d<cdim; ++d) {
    up_dirs[d] = d;
    zero_flux_flags[d] = 0;
  }
  int num_up_dirs = cdim;
  // update velocity space only when field is present (or pkpm model, which always has force update)
  if (field_id != GKYL_FIELD_NULL || up->model_id == GKYL_MODEL_PKPM) {
    for (int d=cdim; d<pdim; ++d) {
      up_dirs[d] = d;
      zero_flux_flags[d] = 1; // zero-flux BCs in vel-space
    }
    num_up_dirs = pdim;
  }
  up->up_vlasov = gkyl_hyper_dg_new(grid, pbasis, up->eqn_vlasov, num_up_dirs, up_dirs, zero_flux_flags, 1, use_gpu);

  up->vlasov_tm = 0.0;
  
  return up;
}

void
gkyl_dg_updater_vlasov_advance(gkyl_dg_updater_vlasov *vlasov,
  const struct gkyl_range *update_rng,
  const struct gkyl_array *aux1, const struct gkyl_array *aux2, 
  const struct gkyl_array *aux3, const struct gkyl_array *aux4, 
  const struct gkyl_array *aux5,  
  const struct gkyl_array* GKYL_RESTRICT fIn,
  struct gkyl_array* GKYL_RESTRICT cflrate, struct gkyl_array* GKYL_RESTRICT rhs)
{
  // Set arrays needed
  // Assumes a particular order of the arrays
  // TO DO: More intelligent way to do these aux field sets? (JJ: 09/08/22)
  if (vlasov->model_id == GKYL_MODEL_SR) {
    gkyl_vlasov_sr_set_auxfields(vlasov->eqn_vlasov, 
      (struct gkyl_dg_vlasov_sr_auxfields) { .qmem = aux1, .p_over_gamma = aux2 });
  }
  else if (vlasov->model_id == GKYL_MODEL_PKPM) {
    gkyl_vlasov_pkpm_set_auxfields(vlasov->eqn_vlasov, 
      (struct gkyl_dg_vlasov_pkpm_auxfields) { 
        .bvar = aux1, .pkpm_prim = aux2, 
        .pkpm_accel_vars = aux3, .g_dist_source = aux4});
  }
  else {
    gkyl_vlasov_set_auxfields(vlasov->eqn_vlasov, 
      (struct gkyl_dg_vlasov_auxfields) { 
        .field = aux1, .ext_field = aux3, 
        .cot_vec = aux4, .alpha_geo = aux5 }); 
  }
  struct timespec wst = gkyl_wall_clock();
  gkyl_hyper_dg_advance(vlasov->up_vlasov, update_rng, fIn, cflrate, rhs);
  vlasov->vlasov_tm += gkyl_time_diff_now_sec(wst);
}

struct gkyl_dg_updater_vlasov_tm
gkyl_dg_updater_vlasov_get_tm(const gkyl_dg_updater_vlasov *vlasov)
{
  return (struct gkyl_dg_updater_vlasov_tm) {
    .vlasov_tm = vlasov->vlasov_tm,
  };
}

void
gkyl_dg_updater_vlasov_release(gkyl_dg_updater_vlasov* vlasov)
{
  gkyl_dg_eqn_release(vlasov->eqn_vlasov);
  gkyl_hyper_dg_release(vlasov->up_vlasov);
  gkyl_free(vlasov);
}

#ifdef GKYL_HAVE_CUDA

void
gkyl_dg_updater_vlasov_advance_cu(gkyl_dg_updater_vlasov *vlasov,
  const struct gkyl_range *update_rng,
  const struct gkyl_array *aux1, const struct gkyl_array *aux2, 
  const struct gkyl_array *aux3, const struct gkyl_array *aux4, 
  const struct gkyl_array *aux5,  
  const struct gkyl_array* GKYL_RESTRICT fIn,
  struct gkyl_array* GKYL_RESTRICT cflrate, struct gkyl_array* GKYL_RESTRICT rhs)
{
  // Set arrays needed
  // Assumes a particular order of the arrays
  // TO DO: More intelligent way to do these aux field sets? (JJ: 09/08/22)
  if (vlasov->model_id == GKYL_MODEL_SR) {
    gkyl_vlasov_sr_set_auxfields(vlasov->eqn_vlasov, 
      (struct gkyl_dg_vlasov_sr_auxfields) { .qmem = aux1, .p_over_gamma = aux2 });
  }
  else if (vlasov->model_id == GKYL_MODEL_PKPM) {
    gkyl_vlasov_pkpm_set_auxfields(vlasov->eqn_vlasov, 
      (struct gkyl_dg_vlasov_pkpm_auxfields) { 
        .bvar = aux1, .pkpm_prim = aux2, 
        .pkpm_accel_vars = aux3, .g_dist_source = aux4});
  }
  else {
    gkyl_vlasov_set_auxfields(vlasov->eqn_vlasov, 
      (struct gkyl_dg_vlasov_auxfields) { 
        .field = aux1, .ext_field = aux3, 
        .cot_vec = aux4, .alpha_geo = aux5 }); 
  }

  struct timespec wst = gkyl_wall_clock();
  gkyl_hyper_dg_advance_cu(vlasov->up_vlasov, update_rng, fIn, cflrate, rhs);
  vlasov->vlasov_tm += gkyl_time_diff_now_sec(wst);
}

#endif

#ifndef GKYL_HAVE_CUDA

void
gkyl_dg_updater_vlasov_advance_cu(gkyl_dg_updater_vlasov *vlasov,
  const struct gkyl_range *update_rng,
  const struct gkyl_array *aux1, const struct gkyl_array *aux2, 
  const struct gkyl_array *aux3, const struct gkyl_array *aux4, 
  const struct gkyl_array *aux5,  
  const struct gkyl_array* GKYL_RESTRICT fIn,
  struct gkyl_array* GKYL_RESTRICT cflrate, struct gkyl_array* GKYL_RESTRICT rhs)
{
  assert(false);
}

#endif
