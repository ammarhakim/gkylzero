#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

void 
gk_species_bflux_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s, struct gk_boundary_fluxes *bflux)
{ 
  // Allocate solver.
  bflux->flux_slvr = gkyl_ghost_surf_calc_new(&s->grid, s->eqn_gyrokinetic, app->cdim, app->use_gpu);

  // Initialize moment solver.
  for (int i=0; i<2*app->cdim; ++i) {
    gk_species_moment_init(app, s, &bflux->gammai[i], "M0");
  }
}

// Computes rhs of the boundary flux.
void
gk_species_bflux_rhs(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_boundary_fluxes *bflux, const struct gkyl_array *fin,
  struct gkyl_array *rhs)
{
  // Zero ghost cells before calculation to ensure there's no residual data.
  for (int j=0; j<app->cdim; ++j) {
    gkyl_array_clear_range(rhs, 0.0, &species->lower_ghost[j]);
    gkyl_array_clear_range(rhs, 0.0, &species->upper_ghost[j]);
  }
  // Ghost cells of the rhs array are filled with the bflux
  // This is overwritten by the boundary conditions and is not being stored,
  // it is only currently used to calculate moments for other applications.
  if (app->use_gpu) {
    gkyl_ghost_surf_calc_advance_cu(bflux->flux_slvr, &species->local_ext, fin, rhs);
  } else {
    gkyl_ghost_surf_calc_advance(bflux->flux_slvr, &species->local_ext, fin, rhs);
  }

  // Only calculating density for use in ambipotential solve.
  for (int j=0; j<app->cdim; ++j) {
    gk_species_moment_calc(&bflux->gammai[2*j], species->lower_ghost[j], app->lower_ghost[j], rhs);
    gk_species_moment_calc(&bflux->gammai[2*j+1], species->upper_ghost[j], app->upper_ghost[j], rhs);
  }
}

void
gk_species_bflux_release(const struct gkyl_gyrokinetic_app *app, const struct gk_boundary_fluxes *bflux)
{
  gkyl_ghost_surf_calc_release(bflux->flux_slvr);
  for (int i=0; i<2*app->cdim; ++i) {
    gk_species_moment_release(app, &bflux->gammai[i]);
  }
}
