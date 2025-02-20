#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

void 
gk_species_bflux_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s, struct gk_boundary_fluxes *bflux)
{
  int cdim = app->cdim;
  int ndim = app->cdim + app->vdim;
  int cells[GKYL_MAX_DIM];
  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  bflux->ghost_flux_slvr = gkyl_ghost_surf_calc_new(&s->grid, s->eqn_gyrokinetic, app->cdim, app->use_gpu);
  
  for (int d=0; d<app->cdim; ++d) {
    for (int i=0; i<ndim; ++i) {
      cells[i] = s->grid.cells[i]; // reset cell values
      lower[i] = s->grid.lower[i];
      upper[i] = s->grid.upper[i];
    }
    cells[d] = 1;
    for (int e=0; e<2; ++e) {
      // Allocate solver.
      struct gkyl_range *skin_r = e==0? &s->lower_skin[d] : &s->upper_skin[d];
      struct gkyl_range *ghost_r = e==0? &s->lower_ghost[d] : &s->upper_ghost[d];
      bflux->flux_slvr[2*d+e] = gkyl_boundary_flux_new(d, e, &s->grid, skin_r, ghost_r, s->eqn_gyrokinetic, app->use_gpu);

      // Initialize moment solver.
      gk_species_moment_init(app, s, &bflux->gammai[2*d+e], "M0");

      // Allocate arrays and ranges for recycling
      lower[d] = e==0? s->grid.lower[d] - s->grid.dx[d] : s->grid.upper[d];
      upper[d] = e==0? s->grid.lower[d] : s->grid.upper[d] + s->grid.dx[d];

      bflux->flux_arr[2*d+e] = mkarr(app->use_gpu, s->basis.num_basis, ghost_r->volume);

      gkyl_range_init(&bflux->flux_r[2*d+e], ndim, ghost_r->lower, ghost_r->upper);
      gkyl_range_init(&bflux->conf_r[2*d+e], cdim, ghost_r->lower, ghost_r->upper);

      gkyl_rect_grid_init(&bflux->boundary_grid[2*d+e], ndim, lower, upper, cells);
      gkyl_rect_grid_init(&bflux->conf_boundary_grid[2*d+e], cdim, lower, upper, cells);
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
    /* if (app->use_gpu) { */
    /*   gkyl_ghost_surf_calc_advance_cu(bflux->ghost_flux_slvr, &species->local_ext, fin, rhs); */
    /* } else { */
    /*   gkyl_ghost_surf_calc_advance(bflux->ghost_flux_slvr, &species->local_ext, fin, rhs); */
    /* } */

    gkyl_array_copy_range_to_range(bflux->flux_arr[2*d+0], rhs, &bflux->flux_r[2*d+0],
      &species->lower_ghost[d]);
    gkyl_array_copy_range_to_range(bflux->flux_arr[2*d+1], rhs, &bflux->flux_r[2*d+1],
      &species->upper_ghost[d]);
    
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
      gkyl_array_release(bflux->flux_arr[2*d+e]);
    }
  }
}
