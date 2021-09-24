#include <assert.h>
#include <gkyl_vlasov_priv.h>

void 
vm_species_lbo_init(struct gkyl_vlasov_app *app, struct vm_species *s, struct lbo_collisions *lbo)
{
  // TO DO: Expose nu_u and nu_vthsq arrays above species object
  //        for cross-species collisions. Just testing for now JJ 09/24/21
  int cdim = app->cdim, vdim = app->vdim;
  int num_up_dirs = vdim;
  int up_dirs[GKYL_MAX_DIM], zero_flux_flags[GKYL_MAX_DIM];
  for (int d=0; d<vdim; ++d) {
    up_dirs[d] = d+cdim;
    zero_flux_flags[d] = 1;
  }

  lbo->nu_sum = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  lbo->nu_u = mkarr(app->use_gpu, vdim*app->confBasis.num_basis, app->local_ext.volume);
  lbo->nu_vthsq = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume); 
  // create collision equation object and solver
  if (app->use_gpu)
  {
    lbo->coll_eqn = gkyl_dg_vlasov_lbo_cu_dev_new(&app->confBasis, &app->basis, &app->local);
    lbo->coll_slvr = gkyl_hyper_dg_cu_dev_new(&s->grid, &app->basis, lbo->coll_eqn, num_up_dirs, up_dirs, zero_flux_flags, 1);
  }
  else
  {
    lbo->coll_eqn = gkyl_dg_vlasov_lbo_new(&app->confBasis, &app->basis, &app->local);
    lbo->coll_slvr = gkyl_hyper_dg_new(&s->grid, &app->basis, lbo->coll_eqn, num_up_dirs, up_dirs, zero_flux_flags, 1);
  }
}

void 
vm_species_lbo_release(const struct gkyl_vlasov_app *app, const struct lbo_collisions *lbo)
{
  gkyl_array_release(lbo->nu_sum);
  gkyl_array_release(lbo->nu_u);
  gkyl_array_release(lbo->nu_vthsq);
  if (app->use_gpu) {
  }
  else {
    gkyl_dg_eqn_release(lbo->coll_eqn);
    gkyl_hyper_dg_release(lbo->coll_slvr);
  }
}
