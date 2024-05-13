#include <assert.h>
#include <gkyl_vlasov_poisson_priv.h>

void 
vp_species_bflux_init(struct gkyl_vlasov_poisson_app *app, struct vp_species *s, struct vp_boundary_fluxes *bflux)
{ 
  // allocate solver
  bflux->flux_slvr = gkyl_ghost_surf_calc_new(&s->grid, s->eqn_vlasov, app->cdim, app->use_gpu);

  // initialize moment solver
  for (int i=0; i<2*app->cdim; ++i) {
    vp_species_moment_init(app, s, &bflux->integ_moms[i], "Integrated");
  }
}

// computes rhs of the boundary flux
void
vp_species_bflux_rhs(gkyl_vlasov_poisson_app *app, const struct vp_species *species,
  struct vp_boundary_fluxes *bflux, const struct gkyl_array *fin,
  struct gkyl_array *rhs)
{
  // zero ghost cells before calculation to ensure there's no residual data
  for (int j=0; j<app->cdim; ++j) {
    gkyl_array_clear_range(rhs, 0.0, &(species->lower_ghost[j]));
    gkyl_array_clear_range(rhs, 0.0, &(species->upper_ghost[j]));
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
    vp_species_moment_calc(&bflux->integ_moms[2*j], species->lower_ghost[j], app->lower_ghost[j], rhs);
    vp_species_moment_calc(&bflux->integ_moms[2*j+1], species->upper_ghost[j], app->upper_ghost[j], rhs);
  }
}

void
vp_species_bflux_release(const struct gkyl_vlasov_poisson_app *app, const struct vp_boundary_fluxes *bflux)
{
  gkyl_ghost_surf_calc_release(bflux->flux_slvr);
  for (int i=0; i<2*app->cdim; ++i) {
    vp_species_moment_release(app, &bflux->integ_moms[i]);
  }
}
