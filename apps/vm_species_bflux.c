#include <assert.h>
#include <gkyl_vlasov_priv.h>

void 
vm_species_bflux_init(struct gkyl_vlasov_app *app, struct vm_species *s, struct vm_boundary_fluxes *bflux)
{
  app->bflux_species[app->num_bflux_species] = s;
  app->num_bflux_species = app->num_bflux_species + 1;
  bflux->flux_slvr = gkyl_ghost_surf_calc_new(&s->grid, s->eqn_vlasov);
}

// computes rhs of the boundary flux
double
vm_species_bflux_rhs(gkyl_vlasov_app *app, const struct vm_species *species,
  struct vm_boundary_fluxes *bflux, const struct gkyl_array *fin, struct gkyl_array *rhs,
  struct gkyl_array *bfluxin[], struct gkyl_array *bfluxout[], int species_idx)
{
  gkyl_ghost_surf_calc_advance(bflux->flux_slvr, &species->local_ext, &app->local_ext, fin, rhs);
  int edge_idx_lower, edge_idx_upper;
  for (int j=0; j<app->cdim; ++j) {
    edge_idx_lower = species_idx + (2*j)*app->num_bflux_species;
    edge_idx_upper = species_idx + (2*j+1)*app->num_bflux_species;
    gkyl_array_copy_range(bfluxin[edge_idx_lower], fin, species->skin_ghost.lower_ghost[j]);
    gkyl_array_copy_range(bfluxin[edge_idx_upper], fin, species->skin_ghost.upper_ghost[j]);
    gkyl_array_copy_range(bfluxout[edge_idx_lower], rhs, species->skin_ghost.lower_ghost[j]);
    gkyl_array_copy_range(bfluxout[edge_idx_upper], rhs, species->skin_ghost.upper_ghost[j]);
  }
  return 0;
}

/* void  */
/* vm_species_bflux_release(const struct gkyl_vlasov_app *app, struct vm_boundary_fluxes *bflux) */
/* { */
  
/* } */
