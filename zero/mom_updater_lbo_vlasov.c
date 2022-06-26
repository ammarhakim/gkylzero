#include <assert.h>
#include <math.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_mom_bcorr_lbo_vlasov.h>
#include <gkyl_mom_calc_bcorr.h>
#include <gkyl_prim_lbo_calc.h>
#include <gkyl_prim_lbo_cross_calc.h>
#include <gkyl_prim_lbo_type.h>
#include <gkyl_prim_lbo_vlasov.h>
#include <gkyl_prim_lbo_vlasov_with_fluid.h>
#include <gkyl_mom_updater_lbo_vlasov.h>
#include <gkyl_mom_updater_lbo_vlasov_priv.h>
#include <gkyl_util.h>

gkyl_mom_updater_lbo_vlasov*
gkyl_mom_updater_lbo_vlasov_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *cbasis,
  const struct gkyl_basis *pbasis, const struct gkyl_range *conf_rng, const double *v_bounds, bool collides_with_fluid, bool use_gpu)
{
  gkyl_mom_updater_lbo_vlasov *up = gkyl_malloc(sizeof(gkyl_mom_updater_lbo_vlasov));

  if(use_gpu) {
    // edge of velocity space corrections to momentum and energy 
    up->bcorr_type = gkyl_mom_bcorr_lbo_vlasov_cu_dev_new(cbasis, pbasis, v_bounds);
    up->bcorr_calc = gkyl_mom_calc_bcorr_cu_dev_new(grid, up->bcorr_type);

    // primitive moment calculators
    if (collides_with_fluid)
      up->coll_prim = gkyl_prim_lbo_vlasov_with_fluid_cu_dev_new(cbasis, pbasis, conf_rng);
    else
      up->coll_prim = gkyl_prim_lbo_vlasov_cu_dev_new(cbasis, pbasis);
    up->coll_pcalc = gkyl_prim_lbo_calc_cu_dev_new(grid, up->coll_prim);
    up->cross_calc = gkyl_prim_lbo_cross_calc_cu_dev_new(grid, up->coll_prim);
  } else {
    // edge of velocity space corrections to momentum and energy 
    up->bcorr_type = gkyl_mom_bcorr_lbo_vlasov_new(cbasis, pbasis, v_bounds);
    up->bcorr_calc = gkyl_mom_calc_bcorr_new(grid, up->bcorr_type);

    // primitive moment calculators
    if (collides_with_fluid)
      up->coll_prim = gkyl_prim_lbo_vlasov_with_fluid_new(cbasis, pbasis, conf_rng);
    else
      up->coll_prim = gkyl_prim_lbo_vlasov_new(cbasis, pbasis);
    up->coll_pcalc = gkyl_prim_lbo_calc_new(grid, up->coll_prim);
    up->cross_calc = gkyl_prim_lbo_cross_calc_new(grid, up->coll_prim);
  }

  return up;
}

void
gkyl_mom_updater_lbo_vlasov_advance(gkyl_mom_updater_lbo_vlasov *up, struct gkyl_basis cbasis,
  const struct gkyl_range *phase_rng, const struct gkyl_range *conf_rng, 
  const struct gkyl_array *fin, const struct gkyl_array *moms, struct gkyl_array *boundary_corrections, 
  struct gkyl_array *uout, struct gkyl_array *vtSqout)
{
  // construct boundary corrections
  gkyl_mom_calc_bcorr_advance(up->bcorr_calc,
    phase_rng, conf_rng, fin, boundary_corrections);
  
  // construct primitive moments
  gkyl_prim_lbo_calc_advance(up->coll_pcalc, cbasis, *conf_rng, 
    moms, boundary_corrections,
    uout, vtSqout);
}

void
gkyl_mom_updater_lbo_cross_vlasov_advance(gkyl_mom_updater_lbo_vlasov *up, struct gkyl_basis cbasis,
  const struct gkyl_range *conf_rng, const struct gkyl_array *greene,
  double self_m, const struct gkyl_array *self_u, const struct gkyl_array *self_vtsq,
  double cross_m, const struct gkyl_array *cross_u, const struct gkyl_array *cross_vtsq, 
  const struct gkyl_array *moms, const struct gkyl_array *boundary_corrections, 
  struct gkyl_array *u_out, struct gkyl_array *vtsq_out)
{
  gkyl_prim_lbo_cross_calc_advance(up->cross_calc,
    cbasis, *conf_rng, greene, 
    self_m, self_u, self_vtsq, 
    cross_m , cross_u, cross_vtsq,
    moms, boundary_corrections,
    u_out, vtsq_out);
}

void
gkyl_mom_updater_lbo_vlasov_release(gkyl_mom_updater_lbo_vlasov* up)
{
  gkyl_mom_type_release(up->bcorr_type);
  gkyl_mom_calc_bcorr_release(up->bcorr_calc);
  gkyl_prim_lbo_type_release(up->coll_prim);
  gkyl_prim_lbo_calc_release(up->coll_pcalc);
  gkyl_free(up);
}

#ifdef GKYL_HAVE_CUDA

void
gkyl_mom_updater_lbo_vlasov_advance_cu(gkyl_mom_updater_lbo_vlasov *up, struct gkyl_basis cbasis,
  const struct gkyl_range *phase_rng, const struct gkyl_range *conf_rng, 
  const struct gkyl_array *fin, const struct gkyl_array *moms, struct gkyl_array *boundary_corrections, 
  struct gkyl_array *uout, struct gkyl_array *vtSqout)
{
  // construct boundary corrections
  gkyl_mom_calc_bcorr_advance_cu(up->bcorr_calc,
    phase_rng, conf_rng, fin, boundary_corrections);
  
  // construct primitive moments
  gkyl_prim_lbo_calc_advance_cu(up->coll_pcalc, cbasis, *conf_rng, 
    moms, boundary_corrections,
    uout, vtSqout);
}

void
gkyl_mom_updater_lbo_cross_vlasov_advance_cu(gkyl_mom_updater_lbo_vlasov *up, struct gkyl_basis cbasis,
  const struct gkyl_range *conf_rng, const struct gkyl_array *greene,
  double self_m, const struct gkyl_array *self_u, const struct gkyl_array *self_vtsq,
  double cross_m, const struct gkyl_array *cross_u, const struct gkyl_array *cross_vtsq, 
  const struct gkyl_array *moms, const struct gkyl_array *boundary_corrections, 
  struct gkyl_array *u_out, struct gkyl_array *vtsq_out)
{
  gkyl_prim_lbo_cross_calc_advance_cu(up->cross_calc,
    cbasis, *conf_rng, greene, 
    self_m, self_u, self_vtsq, 
    cross_m , cross_u, cross_vtsq,
    moms, boundary_corrections,
    u_out, vtsq_out);
}

#endif

#ifndef GKYL_HAVE_CUDA

void
gkyl_mom_updater_lbo_vlasov_advance_cu(gkyl_mom_updater_lbo_vlasov *up, struct gkyl_basis cbasis,
  const struct gkyl_range *phase_rng, const struct gkyl_range *conf_rng, 
  const struct gkyl_array *fin, const struct gkyl_array *moms, struct gkyl_array *boundary_corrections, 
  struct gkyl_array *uout, struct gkyl_array *vtSqout)
{
  assert(false);
}

void
gkyl_mom_updater_lbo_cross_vlasov_advance_cu(gkyl_mom_updater_lbo_vlasov *up, struct gkyl_basis cbasis,
  const struct gkyl_range *conf_rng, const struct gkyl_array *greene,
  double self_m, const struct gkyl_array *self_u, const struct gkyl_array *self_vtsq,
  double cross_m, const struct gkyl_array *cross_u, const struct gkyl_array *cross_vtsq, 
  const struct gkyl_array *moms, const struct gkyl_array *boundary_corrections, 
  struct gkyl_array *u_out, struct gkyl_array *vtsq_out)
{
  assert(false);
}

#endif
