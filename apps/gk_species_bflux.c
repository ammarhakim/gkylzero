#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

void 
gk_species_bflux_init(struct gkyl_gyrokinetic_app *app, struct gk_species *gk_s, struct gk_boundary_fluxes *bflux, bool is_diagnostic)
{ 
  for (int d=0; d<app->cdim; ++d) {
    for (int e=0; e<2; ++e) {
      // Allocate solver.
      struct gkyl_range *skin_r = e==0? &gk_s->lower_skin[d] : &gk_s->upper_skin[d];
      struct gkyl_range *ghost_r = e==0? &gk_s->lower_ghost[d] : &gk_s->upper_ghost[d];
      bflux->flux_slvr[2*d+e] = gkyl_boundary_flux_new(d, e, &gk_s->grid,
        skin_r, ghost_r, gk_s->eqn_gyrokinetic, is_diagnostic, app->use_gpu);

      // Initialize moment solver.
      gk_species_moment_init(app, gk_s, &bflux->gammai[2*d+e], "M0", false);
    }
  }

  bflux->is_diagnostic = false;
  if (is_diagnostic) {
    bflux->is_diagnostic = true;
    bflux->f = gkyl_malloc((2*app->cdim)*sizeof(struct gkyl_array *));
    bflux->f1 = gkyl_malloc((2*app->cdim)*sizeof(struct gkyl_array *));
    bflux->fnew = gkyl_malloc((2*app->cdim)*sizeof(struct gkyl_array *));
    for (int d=0; d<app->cdim; ++d) {
      for (int e=0; e<2; ++e) {
        bflux->f[2*d+e] = mkarr(app->use_gpu, app->basis.num_basis, gk_s->local_ext.volume);
        bflux->f1[2*d+e] = mkarr(app->use_gpu, app->basis.num_basis, gk_s->local_ext.volume);
        bflux->fnew[2*d+e] = mkarr(app->use_gpu, app->basis.num_basis, gk_s->local_ext.volume);
      }
    }

    gk_species_moment_init(app, gk_s, &bflux->integ_moms,
      gk_s->info.integrated_hamiltonian_moments? "HamiltonianMoments" : "FourMoments", true);
    for (int d=0; d<app->cdim; ++d) {
      for (int e=0; e<2; ++e)
        bflux->intmom[2*d+e] = gkyl_dynvec_new(GKYL_DOUBLE, bflux->integ_moms.num_mom);
    }
    bflux->is_first_intmom_write_call = true;
  }
}

void
gk_species_bflux_rhs(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_boundary_fluxes *bflux, const struct gkyl_array *fin,
  struct gkyl_array *rhs)
{
  // Computes rhs of the boundary flux.

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
gk_species_bflux_rhs_diag(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_boundary_fluxes *bflux, const struct gkyl_array *fin,
  struct gkyl_array **bflux_out)
{
  // Computes rhs of the boundary flux for diagnostics.

  // Only calculating density for use in ambipotential solve.
  for (int d=0; d<app->cdim; ++d) {
    for (int e=0; e<2; ++e) {
      // Zero ghost cells before calculation to ensure there's no residual data.
      gkyl_array_clear_range(bflux_out[2*d+e], 0.0, e==0? &species->lower_ghost[d] : &species->upper_ghost[d]);

      // Ghost cells of the rhs array are filled with the bflux
      // This is overwritten by the boundary conditions and is not being stored,
      // it is only currently used to calculate moments for other applications.
      gkyl_boundary_flux_advance(bflux->flux_slvr[2*d+e], fin, bflux_out[2*d+e]);
    }
  }
}

void
gk_species_bflux_clear(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **bflux_in, double val)
{
  if (bflux->is_diagnostic) {
    for (int d=0; d<app->cdim; ++d) {
      for (int e=0; e<2; ++e) {
        gkyl_array_clear(bflux_in[2*d+e], val);
      }
    }
  }
}

void
gk_species_bflux_scale(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **bflux_in, double val)
{
  if (bflux->is_diagnostic) {
    for (int d=0; d<app->cdim; ++d) {
      for (int e=0; e<2; ++e) {
        gkyl_array_scale(bflux_in[2*d+e], val);
      }
    }
  }
}

void
gk_species_bflux_forward_euler(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **bflux_out, double dt, const struct gkyl_array **bflux_in)
{
  if (bflux->is_diagnostic) {
    for (int d=0; d<app->cdim; ++d) {
      for (int e=0; e<2; ++e) {
        gkyl_array_accumulate(gkyl_array_scale(bflux_out[2*d+e], dt), 1.0, bflux_in[2*d+e]);
      }
    }
  }
}

void
gk_species_bflux_combine(gkyl_gyrokinetic_app *app, struct gk_species *gk_s, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **bflux_out, double fac1, struct gkyl_array **bflux_in1, double fac2, struct gkyl_array **bflux_in2)
{
  if (bflux->is_diagnostic) {
    for (int d=0; d<app->cdim; ++d) {
      for (int e=0; e<2; ++e) {
        array_combine(bflux_out[2*d+e], fac1, bflux_in1[2*d+e], fac2, bflux_in2[2*d+e], &gk_s->local_ext);
      }
    }
  }
}

void
gk_species_bflux_copy(gkyl_gyrokinetic_app *app, struct gk_species *gk_s, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **bflux_out, struct gkyl_array **bflux_in)
{
  if (bflux->is_diagnostic) {
    for (int d=0; d<app->cdim; ++d) {
      for (int e=0; e<2; ++e) {
        gkyl_array_copy_range(bflux_out[2*d+e], bflux_in[2*d+e], &gk_s->local_ext);
      }
    }
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

  if (bflux->is_diagnostic) {
    for (int d=0; d<app->cdim; ++d) {
      for (int e=0; e<2; ++e) {
        gkyl_array_release(bflux->f[2*d+e]);
        gkyl_array_release(bflux->f1[2*d+e]);
        gkyl_array_release(bflux->fnew[2*d+e]);
      }
    }
    gkyl_free(bflux->f);
    gkyl_free(bflux->f1);
    gkyl_free(bflux->fnew);

    gk_species_moment_release(app, &bflux->integ_moms); 

    for (int d=0; d<app->cdim; ++d) {
      for (int e=0; e<2; ++e) {
        gkyl_dynvec_release(bflux->intmom[2*d+e]);
      }
    }
  }
}
