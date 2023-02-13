#include <assert.h>
#include <gkyl_vlasov_priv.h>

void 
vm_species_bflux_init(struct gkyl_vlasov_app *app, struct vm_species *s, struct vm_boundary_fluxes *bflux)
{ 
  // allocate solver
  bflux->flux_slvr = gkyl_ghost_surf_calc_new(&s->grid, s->eqn_vlasov, app->cdim, app->use_gpu);

  // initialize moment solver
  for (int i=0; i<2*app->cdim; ++i) {
    vm_species_moment_init(app, s, &bflux->integ_moms[i], "Integrated");
  }
}

// computes rhs of the boundary flux
void
vm_species_bflux_rhs(gkyl_vlasov_app *app, const struct vm_species *species,
  struct vm_boundary_fluxes *bflux, const struct gkyl_array *fin,
  struct gkyl_array *rhs)
{
  // zero ghost cells before calculation to ensure there's no residual data
  for (int j=0; j<app->cdim; ++j) {
    gkyl_array_clear_range(rhs, 0.0, species->skin_ghost.lower_ghost[j]);
    gkyl_array_clear_range(rhs, 0.0, species->skin_ghost.upper_ghost[j]);
  }
  // ghost cells of the rhs array are filled with the bflux
  // This is overwritten by the boundary conditions and is not being stored,
  // it is only currently used to calculate moments for other applications
  if (app->use_gpu) {
    gkyl_ghost_surf_calc_advance_cu(bflux->flux_slvr, &species->local_ext, fin, rhs);
  } else {
    gkyl_ghost_surf_calc_advance(bflux->flux_slvr, &species->local_ext, fin, rhs);
  }

  // only calculating integrated moments for use in the bflux source for now,
  // others can be added if applications require
  for (int j=0; j<app->cdim; ++j) {
    vm_species_moment_calc(&bflux->integ_moms[2*j], species->skin_ghost.lower_ghost[j], app->skin_ghost.lower_ghost[j], rhs);
    vm_species_moment_calc(&bflux->integ_moms[2*j+1], species->skin_ghost.upper_ghost[j], app->skin_ghost.upper_ghost[j], rhs);
  }
}

void
vm_species_bflux_release(const struct gkyl_vlasov_app *app, const struct vm_boundary_fluxes *bflux)
{
  gkyl_ghost_surf_calc_release(bflux->flux_slvr);
  for (int i=0; i<2*app->cdim; ++i) {
    vm_species_moment_release(app, &bflux->integ_moms[i]);
  }
}
