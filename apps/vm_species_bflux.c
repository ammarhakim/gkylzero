#include <assert.h>
#include <gkyl_vlasov_priv.h>

void 
vm_species_bflux_init(struct gkyl_vlasov_app *app, struct vm_species *s, struct vm_boundary_fluxes *bflux)
{
  // allocate buffer for applying periodic BCs
  long buff_sz = 0;
  // compute buffer size needed
  for (int d=0; d<app->cdim; ++d) {
    long vol = s->skin_ghost.lower_skin[d].volume;
    buff_sz = buff_sz > vol ? buff_sz : vol;
  }
  bflux->ghost_buffer = mkarr(app->use_gpu, app->basis.num_basis, buff_sz);
  
  app->bflux_species[app->num_bflux_species] = s;
  app->num_bflux_species = app->num_bflux_species + 1;
  bflux->flux_slvr = gkyl_ghost_surf_calc_new(&s->grid, s->eqn_vlasov);

  for (int i=0; i<2*app->cdim; ++i) {
    bflux->bf[i] = mkarr(app->use_gpu, app->basis.num_basis, buff_sz);
    bflux->bf1[i] = mkarr(app->use_gpu, app->basis.num_basis, buff_sz);
    bflux->bfnew[i] = mkarr(app->use_gpu, app->basis.num_basis, buff_sz);
    
    vm_species_moment_init(app, s, &bflux->moms[i], "FiveMoments");
    vm_species_moment_init(app, s, &bflux->m0[i], "M0");
    vm_species_moment_init(app, s, &bflux->m1i[i], "M1i");
    vm_species_moment_init(app, s, &bflux->integ_moms[i], "Integrated");
  }
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

    gkyl_array_copy_to_buffer(bflux->ghost_buffer->data, fin, species->skin_ghost.lower_ghost[j]);
    gkyl_array_copy_from_buffer(bfluxin[edge_idx_lower], bflux->ghost_buffer->data, species->skin_ghost.lower_ghost[j]);

    gkyl_array_copy_to_buffer(bflux->ghost_buffer->data, fin, species->skin_ghost.upper_ghost[j]);
    gkyl_array_copy_from_buffer(bfluxin[edge_idx_upper], bflux->ghost_buffer->data, species->skin_ghost.upper_ghost[j]);
    
    gkyl_array_copy_to_buffer(bflux->ghost_buffer->data, rhs, species->skin_ghost.lower_ghost[j]);
    gkyl_array_copy_from_buffer(bfluxout[edge_idx_lower], bflux->ghost_buffer->data, species->skin_ghost.lower_ghost[j]);
    
    gkyl_array_copy_to_buffer(bflux->ghost_buffer->data, rhs, species->skin_ghost.upper_ghost[j]);
    gkyl_array_copy_from_buffer(bfluxout[edge_idx_upper], bflux->ghost_buffer->data, species->skin_ghost.upper_ghost[j]);
    
    vm_species_moment_calc(&bflux->moms[2*j], species->skin_ghost.lower_ghost[j], app->skin_ghost.lower_ghost[j], rhs);
    vm_species_moment_calc(&bflux->moms[2*j+1], species->skin_ghost.upper_ghost[j], app->skin_ghost.upper_ghost[j], rhs);

    vm_species_moment_calc(&bflux->m0[2*j], species->skin_ghost.lower_ghost[j], app->skin_ghost.lower_ghost[j], rhs);
    vm_species_moment_calc(&bflux->m0[2*j+1], species->skin_ghost.upper_ghost[j], app->skin_ghost.upper_ghost[j], rhs);

    vm_species_moment_calc(&bflux->m1i[2*j], species->skin_ghost.lower_ghost[j], app->skin_ghost.lower_ghost[j], rhs);
    vm_species_moment_calc(&bflux->m1i[2*j+1], species->skin_ghost.upper_ghost[j], app->skin_ghost.upper_ghost[j], rhs);

    vm_species_moment_calc(&bflux->integ_moms[2*j], species->skin_ghost.lower_ghost[j], app->skin_ghost.lower_ghost[j], rhs);
    vm_species_moment_calc(&bflux->integ_moms[2*j+1], species->skin_ghost.upper_ghost[j], app->skin_ghost.upper_ghost[j], rhs);
  }
  return 0;
}

/* void  */
/* vm_species_bflux_release(const struct gkyl_vlasov_app *app, struct vm_boundary_fluxes *bflux) */
/* { */
  
/* } */
