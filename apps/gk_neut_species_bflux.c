#include <assert.h>
#include <gkyl_mom_canonical_pb.h>
#include <gkyl_dg_updater_moment.h>
#include <gkyl_gyrokinetic_priv.h>

void 
gk_neut_species_bflux_init(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s, struct gk_boundary_fluxes *bflux)
{ 
  // allocate solver
  bflux->ghost_flux_slvr = gkyl_ghost_surf_calc_new(&s->grid, s->eqn_vlasov, app->cdim, app->use_gpu);
  int cdim = app->cdim;
  int ndim = app->cdim + app->vdim + 1;
  int cells[GKYL_MAX_DIM], ghost[GKYL_MAX_DIM];
  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  for (int d=0; d<ndim; ++d) {
    cells[d] = s->grid.cells[d];
    lower[d] = s->grid.lower[d];
    upper[d] = s->grid.upper[d];
    ghost[d] = 0;
  }

  // initialize moment solver
  for (int i=0; i<app->cdim; ++i) {
    for (int d=0; d<cdim; ++d) {
      cells[d] = s->grid.cells[d];
    }
    cells[i] = 1; 
    
    bflux->flux_arr[2*i] = mkarr(app->use_gpu, s->basis.num_basis, s->lower_ghost[i].volume);
    bflux->flux_arr[2*i+1] = mkarr(app->use_gpu, s->basis.num_basis, s->upper_ghost[i].volume);

    gkyl_range_init(&bflux->flux_r[2*i], ndim, s->lower_ghost[i].lower, s->lower_ghost[i].upper);
    gkyl_range_init(&bflux->flux_r[2*i+1], ndim, s->upper_ghost[i].lower, s->upper_ghost[i].upper);

    gkyl_range_init(&bflux->conf_r[2*i], app->cdim, s->lower_ghost[i].lower,
      s->lower_ghost[i].upper);
    gkyl_range_init(&bflux->conf_r[2*i+1], app->cdim, s->upper_ghost[i].lower,
      s->upper_ghost[i].upper);

    upper[i] = s->grid.lower[i] + s->grid.dx[i];

    gkyl_rect_grid_init(&bflux->boundary_grid[2*i], ndim, lower, upper, cells);
    gkyl_rect_grid_init(&bflux->conf_boundary_grid[2*i], cdim, lower, upper, cells);

    upper[i] = s->grid.upper[i];
    lower[i] = s->grid.upper[i] - s->grid.dx[i];

    gkyl_rect_grid_init(&bflux->boundary_grid[2*i+1], ndim, lower, upper, cells);
    gkyl_rect_grid_init(&bflux->conf_boundary_grid[2*i+1], cdim, lower, upper, cells);
  }
}

// computes rhs of the boundary flux
void
gk_neut_species_bflux_rhs(gkyl_gyrokinetic_app *app, const struct gk_neut_species *species,
  struct gk_boundary_fluxes *bflux, const struct gkyl_array *fin,
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
    gkyl_ghost_surf_calc_advance_cu(bflux->ghost_flux_slvr, &species->local_ext, fin, rhs);
  } else {
    gkyl_ghost_surf_calc_advance(bflux->ghost_flux_slvr, &species->local_ext, fin, rhs);
  }

  // only calculating integrated moments for use in the bflux source for now,
  // others can be added if applications require
  for (int j=0; j<app->cdim; ++j) {
    gkyl_array_copy_range_to_range(bflux->flux_arr[2*j], rhs, &bflux->flux_r[2*j],
      &species->lower_ghost[j]);
    gkyl_array_copy_range_to_range(bflux->flux_arr[2*j+1], rhs, &bflux->flux_r[2*j+1],
      &species->upper_ghost[j]);
  }
}

void
gk_neut_species_bflux_release(const struct gkyl_gyrokinetic_app *app, const struct gk_boundary_fluxes *bflux)
{
  gkyl_ghost_surf_calc_release(bflux->ghost_flux_slvr);
  for (int i=0; i<2*app->cdim; ++i) {
    gkyl_array_release(bflux->flux_arr[i]);
  }
}
