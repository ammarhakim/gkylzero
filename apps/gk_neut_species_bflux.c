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
  int cells[GKYL_MAX_DIM];
  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];

  // initialize moment solver
  for (int d=0; d<app->cdim; ++d) {
    for (int i=0; i<ndim; ++i) {
      cells[i] = s->grid.cells[i];
      lower[i] = s->grid.lower[i];
      upper[i] = s->grid.upper[i];
    }
    cells[d] = 1;
    for (int e=0; e<2; ++e) {
      // Allocate arrays and ranges for recycling
      struct gkyl_range *skin_r = e==0? &s->lower_skin[d] : &s->upper_skin[d];
      struct gkyl_range *ghost_r = e==0? &s->lower_ghost[d] : &s->upper_ghost[d];
      
      lower[d] = e==0? s->grid.lower[d] : (s->grid.upper[d] - s->grid.dx[d]);
      upper[d] = e==0? (s->grid.lower[d] + s->grid.dx[d]) : s->grid.upper[d];

      bflux->flux_arr[2*d+e] = mkarr(app->use_gpu, s->basis.num_basis, ghost_r->volume);
      gkyl_range_init(&bflux->flux_r[2*d+e], ndim, ghost_r->lower, ghost_r->upper);
      gkyl_range_init(&bflux->conf_r[2*d+e], cdim, ghost_r->lower, ghost_r->upper);

      gkyl_rect_grid_init(&bflux->boundary_grid[2*d+e], ndim, lower, upper, cells);
      gkyl_rect_grid_init(&bflux->conf_boundary_grid[2*d+e], cdim, lower, upper, cells);
    }
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
