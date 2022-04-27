#include <assert.h>
#include <gkyl_vlasov_priv.h>

// list of valid moment names
static const char *const valid_moment_names[] = { "M0", "M1i", "M2ij", "M2", "M3i", "M3ijk", "FiveMoments" };

// check if name of moment is valid or not
static bool
is_moment_name_valid(const char *nm)
{
  int n = sizeof(valid_moment_names)/sizeof(valid_moment_names[0]);
  for (int i=0; i<n; ++i)
    if (strcmp(valid_moment_names[i], nm) == 0)
      return 1;
  return 0;
}

// initialize species moment object
void
vm_species_moment_init(struct gkyl_vlasov_app *app, struct vm_species *s,
  struct vm_species_moment *sm, const char *nm)
{
  assert(is_moment_name_valid(nm));

  sm->use_gpu = app->use_gpu;
  
  if (app->use_gpu) {
    struct gkyl_mom_type *mtype;
    if (s->field_id == GKYL_FIELD_SR_E_B) {
      mtype = gkyl_mom_vlasov_sr_cu_dev_new(&app->confBasis, &app->basis, &s->local_vel, nm);
      gkyl_mom_vlasov_sr_set_auxfields(mtype, 
        (struct gkyl_mom_vlasov_sr_auxfields) { .p_over_gamma = s->p_over_gamma });
    }
    else {
      mtype = gkyl_mom_vlasov_cu_dev_new(&app->confBasis, &app->basis, nm);
    }
    sm->mcalc = gkyl_mom_calc_cu_dev_new(&s->grid, mtype);

    sm->marr = mkarr(app->use_gpu, mtype->num_mom*app->confBasis.num_basis,
      app->local_ext.volume);

    sm->marr_host = mkarr(false, mtype->num_mom*app->confBasis.num_basis,
      app->local_ext.volume);

    gkyl_mom_type_release(mtype);
  }
  else {
    struct gkyl_mom_type *mtype;
    if (s->field_id == GKYL_FIELD_SR_E_B) {
      mtype = gkyl_mom_vlasov_sr_new(&app->confBasis, &app->basis, &s->local_vel, nm);
      gkyl_mom_vlasov_sr_set_auxfields(mtype, 
        (struct gkyl_mom_vlasov_sr_auxfields) { .p_over_gamma = s->p_over_gamma });
    }
    else {
      mtype = gkyl_mom_vlasov_new(&app->confBasis, &app->basis, nm);
    }
    sm->mcalc = gkyl_mom_calc_new(&s->grid, mtype);

    sm->marr = mkarr(app->use_gpu, mtype->num_mom*app->confBasis.num_basis,
      app->local_ext.volume);
    sm->marr_host = sm->marr;

    gkyl_mom_type_release(mtype);
  }
}

void
vm_species_moment_calc(const struct vm_species_moment *sm,
  const struct gkyl_range phase_rng, const struct gkyl_range conf_rng,
  const struct gkyl_array *fin)
{
  if (sm->use_gpu)
    gkyl_mom_calc_advance_cu(sm->mcalc, &phase_rng, &conf_rng, fin, sm->marr);
  else
    gkyl_mom_calc_advance(sm->mcalc, &phase_rng, &conf_rng, fin, sm->marr);
}

// release memory for moment data object
void
vm_species_moment_release(const struct gkyl_vlasov_app *app, const struct vm_species_moment *sm)
{
  if (app->use_gpu)
    gkyl_array_release(sm->marr_host);

  gkyl_mom_calc_release(sm->mcalc);
  gkyl_array_release(sm->marr);
}
