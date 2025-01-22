#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

void 
gk_species_bflux_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s, struct gk_boundary_fluxes *bflux)
{ 
  for (int d=0; d<app->cdim; ++d) {
    for (int e=0; e<2; ++e) {
      // Allocate solver.
      struct gkyl_range *skin_r = e==0? &s->lower_skin[d] : &s->upper_skin[d];
      struct gkyl_range *ghost_r = e==0? &s->lower_ghost[d] : &s->upper_ghost[d];
      bflux->flux_slvr[2*d+e] = gkyl_boundary_flux_new(d, e, &s->grid, skin_r, ghost_r, s->eqn_gyrokinetic, app->use_gpu);

      // Initialize moment solver.
      gk_species_moment_init(app, s, &bflux->gammai[2*d+e], "M0");
    }
  }
}

// Computes rhs of the boundary flux.
void
gk_species_bflux_rhs(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_boundary_fluxes *bflux, const struct gkyl_array *fin,
  struct gkyl_array *rhs)
{
  // Zero ghost cells before calculation to ensure there's no residual data.
  for (int d=0; d<app->cdim; ++d) {
    gkyl_array_clear_range(rhs, 0.0, &species->lower_ghost[d]);
    gkyl_array_clear_range(rhs, 0.0, &species->upper_ghost[d]);
  }
  // Only calculating density for use in ambipotential solve.
  for (int d=0; d<app->cdim; ++d) {
    for (int e=0; e<2; ++e) {
      // Ghost cells of the rhs array are filled with the bflux
      // This is overwritten by the boundary conditions and is not being stored,
      // it is only currently used to calculate moments for other applications.
      gkyl_boundary_flux_advance(bflux->flux_slvr[2*d+e], fin, rhs);
    }

    gk_species_moment_calc(&bflux->gammai[2*d+0], species->lower_ghost[d], app->lower_ghost[d], rhs);
    gk_species_moment_calc(&bflux->gammai[2*d+1], species->upper_ghost[d], app->upper_ghost[d], rhs);
  }
}

void
gk_species_bflux_release(const struct gkyl_gyrokinetic_app *app, const struct gk_boundary_fluxes *bflux)
{
  for (int d=0; d<app->cdim; ++d) {
    for (int e=0; e<2; ++e) {
      gkyl_boundary_flux_release(bflux->flux_slvr[2*d+e]);
      gk_species_moment_release(app, &bflux->gammai[2*d+e]);
    }
  }
}
