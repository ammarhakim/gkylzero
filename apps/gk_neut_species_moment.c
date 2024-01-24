#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

// initialize neutral species moment object
void
gk_neut_species_moment_init(struct gkyl_gyrokinetic_oneb *oneb, struct gk_neut_species *s,
  struct gk_neut_species_moment *sm, const char *nm)
{
  assert(is_neut_moment_name_valid(nm));

  bool is_integrated = strcmp(nm, "Integrated") == 0;

  sm->mcalc = gkyl_dg_updater_moment_new(&s->grid, &oneb->confBasis, 
    &oneb->neut_basis, &oneb->local, &s->local_vel, s->model_id, 0,
    nm, is_integrated, s->info.mass, oneb->use_gpu);    

  int num_mom = gkyl_dg_updater_moment_num_mom(sm->mcalc);

  if (is_integrated) {
    sm->marr = mkarr(oneb->use_gpu, num_mom, oneb->local_ext.volume);
    sm->marr_host = sm->marr;
    if (oneb->use_gpu)
      sm->marr_host = mkarr(false, num_mom, oneb->local_ext.volume); 
    // Bin Op memory for rescaling moment by inverse of Jacobian
    if (oneb->use_gpu)
      sm->mem_geo = gkyl_dg_bin_op_mem_cu_dev_new(oneb->local.volume, num_mom);  
    else   
      sm->mem_geo = gkyl_dg_bin_op_mem_new(oneb->local.volume, num_mom);  
  }
  else {
    sm->marr = mkarr(oneb->use_gpu, num_mom*oneb->confBasis.num_basis, oneb->local_ext.volume);
    sm->marr_host = sm->marr;
    if (oneb->use_gpu)
      sm->marr_host = mkarr(false, num_mom*oneb->confBasis.num_basis, oneb->local_ext.volume);
    // Bin Op memory for rescaling moment by inverse of Jacobian
    if (oneb->use_gpu)
      sm->mem_geo = gkyl_dg_bin_op_mem_cu_dev_new(oneb->local.volume, num_mom*oneb->confBasis.num_basis);
    else
      sm->mem_geo = gkyl_dg_bin_op_mem_new(oneb->local.volume, num_mom*oneb->confBasis.num_basis);
  }
}

void
gk_neut_species_moment_calc(const struct gk_neut_species_moment *sm,
  const struct gkyl_range phase_rng, const struct gkyl_range conf_rng,
  const struct gkyl_array *fin)
{
  gkyl_dg_updater_moment_advance(sm->mcalc, &phase_rng, &conf_rng, fin, sm->marr);
}

// release memory for moment data object
void
gk_neut_species_moment_release(const struct gkyl_gyrokinetic_oneb *app, const struct gk_neut_species_moment *sm)
{
  if (app->use_gpu)
    gkyl_array_release(sm->marr_host);

  gkyl_dg_updater_moment_release(sm->mcalc);
  gkyl_array_release(sm->marr);

  gkyl_dg_bin_op_mem_release(sm->mem_geo);
}
