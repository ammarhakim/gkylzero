#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

void 
gk_species_bflux_init(struct gkyl_gyrokinetic_app *app, struct gk_species *gk_s, struct gk_boundary_fluxes *bflux, bool is_diagnostic)
{ 
  bflux->is_diagnostic = is_diagnostic;

  for (int d=0; d<app->cdim; ++d) {
    for (int e=0; e<2; ++e) {
      // Allocate solver.
      struct gkyl_range *skin_r = e==0? &gk_s->lower_skin[d] : &gk_s->upper_skin[d];
      struct gkyl_range *ghost_r = e==0? &gk_s->lower_ghost[d] : &gk_s->upper_ghost[d];
      bflux->flux_slvr[2*d+e] = gkyl_boundary_flux_new(d, e, &gk_s->grid,
        skin_r, ghost_r, gk_s->eqn_gyrokinetic, is_diagnostic, app->use_gpu);

      if (!is_diagnostic) {
        // Initialize moment solver.
        gk_species_moment_init(app, gk_s, &bflux->gammai[2*d+e], "M0", false);
      }
    }
  }

  if (is_diagnostic) {
    gk_species_moment_init(app, gk_s, &bflux->moms_op,
      gk_s->info.integrated_hamiltonian_moments? "HamiltonianMoments" : "FourMoments", false);

    bflux->f = gkyl_malloc((2*app->cdim)*sizeof(struct gkyl_array *));
    bflux->f1 = gkyl_malloc((2*app->cdim)*sizeof(struct gkyl_array *));
    bflux->fnew = gkyl_malloc((2*app->cdim)*sizeof(struct gkyl_array *));
    for (int d=0; d<app->cdim; ++d) {
      for (int e=0; e<2; ++e) {
        bflux->f[2*d+e] = mkarr(app->use_gpu, bflux->moms_op.num_mom*app->confBasis.num_basis, app->local_ext.volume);
        bflux->f1[2*d+e] = mkarr(app->use_gpu, bflux->moms_op.num_mom*app->confBasis.num_basis, app->local_ext.volume);
        bflux->fnew[2*d+e] = mkarr(app->use_gpu, bflux->moms_op.num_mom*app->confBasis.num_basis, app->local_ext.volume);
      }
    }

    bflux->integ_op = gkyl_array_integrate_new(&app->grid, &app->confBasis,
      bflux->moms_op.num_mom, GKYL_ARRAY_INTEGRATE_OP_NONE, app->use_gpu);
    if (app->use_gpu) {
      bflux->int_moms_local = gkyl_cu_malloc(bflux->moms_op.num_mom*sizeof(double));
      bflux->int_moms_global = gkyl_cu_malloc(bflux->moms_op.num_mom*sizeof(double));
    }
    else {
      bflux->int_moms_local = gkyl_malloc(bflux->moms_op.num_mom*sizeof(double));
      bflux->int_moms_global = gkyl_malloc(bflux->moms_op.num_mom*sizeof(double));
    }
    for (int d=0; d<app->cdim; ++d) {
      for (int e=0; e<2; ++e)
        bflux->intmom[2*d+e] = gkyl_dynvec_new(GKYL_DOUBLE, bflux->moms_op.num_mom);
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
  struct gk_boundary_fluxes *bflux, const struct gkyl_array *fin, struct gkyl_array *rhs,
  struct gkyl_array **bflux_moms)
{
  // Computes rhs of the boundary flux for diagnostics.

  // Only calculating density for use in ambipotential solve.
  for (int d=0; d<app->cdim; ++d) {
    for (int e=0; e<2; ++e) {
      if (species->info.integrated_hamiltonian_moments) {
        // Apply BC to phi so it is defined in the ghost cell.
        // Fill the ghost with the skin evaluated at the boundary.
        gkyl_bc_basic_advance(app->bc_op[2*d+e], app->bc_buffer, app->field->phi_smooth);
      }
      gkyl_array_scale_range(app->field->phi_smooth, 0.5, e==0? &app->lower_ghost[d] : &app->upper_ghost[d]);

      // Ghost cells of the rhs array are filled with the bflux
      // This is overwritten by the boundary conditions and is not being stored,
      // it is only currently used to calculate moments for other applications.
      gkyl_boundary_flux_advance(bflux->flux_slvr[2*d+e], fin, rhs);

      // Compute moments of boundary fluxes.
      gk_species_moment_calc(&bflux->moms_op, e==0? species->lower_ghost[d] : species->upper_ghost[d],
        e==0? app->lower_ghost[d] : app->upper_ghost[d], rhs);

      gkyl_array_copy_range(bflux_moms[2*d+e], bflux->moms_op.marr, e==0? &app->lower_ghost[d] : &app->upper_ghost[d]);
    }
  }
}

void
gk_species_bflux_clear(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **fin, double val)
{
  if (bflux->is_diagnostic) {
    for (int d=0; d<app->cdim; ++d) {
      for (int e=0; e<2; ++e) {
        gkyl_array_clear(fin[2*d+e], val);
      }
    }
  }
}

void
gk_species_bflux_scale(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **fin, double val)
{
  if (bflux->is_diagnostic) {
    for (int d=0; d<app->cdim; ++d) {
      for (int e=0; e<2; ++e) {
        gkyl_array_scale(fin[2*d+e], val);
      }
    }
  }
}

void
gk_species_bflux_forward_euler(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **fout, double dt, const struct gkyl_array **fin)
{
  if (bflux->is_diagnostic) {
    for (int d=0; d<app->cdim; ++d) {
      for (int e=0; e<2; ++e) {
        gkyl_array_accumulate(gkyl_array_scale(fout[2*d+e], dt), 1.0, fin[2*d+e]);
      }
    }
  }
}

void
gk_species_bflux_combine_range(gkyl_gyrokinetic_app *app, struct gk_species *gk_s, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **fout, double fac1, struct gkyl_array **fin1, double fac2, struct gkyl_array **fin2, struct gkyl_range *range)
{
  if (bflux->is_diagnostic) {
    for (int d=0; d<app->cdim; ++d) {
      for (int e=0; e<2; ++e) {
        array_combine(fout[2*d+e], fac1, fin1[2*d+e], fac2, fin2[2*d+e], range);
      }
    }
  }
}

void
gk_species_bflux_copy_range(gkyl_gyrokinetic_app *app, struct gk_species *gk_s, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **fout, struct gkyl_array **fin, struct gkyl_range *range)
{
  if (bflux->is_diagnostic) {
    for (int d=0; d<app->cdim; ++d) {
      for (int e=0; e<2; ++e) {
        gkyl_array_copy_range(fout[2*d+e], fin[2*d+e], range);
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
    }
  }

  if (bflux->is_diagnostic) {
    gk_species_moment_release(app, &bflux->moms_op); 
    gkyl_array_integrate_release(bflux->integ_op);

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

    if (app->use_gpu) {
      gkyl_cu_free(bflux->int_moms_local);
      gkyl_cu_free(bflux->int_moms_global);
    }
    else {
      gkyl_free(bflux->int_moms_local);
      gkyl_free(bflux->int_moms_global);
    }
    for (int d=0; d<app->cdim; ++d) {
      for (int e=0; e<2; ++e) {
        gkyl_dynvec_release(bflux->intmom[2*d+e]);
      }
    }
  }
  else {
    for (int d=0; d<app->cdim; ++d) {
      for (int e=0; e<2; ++e) {
        gk_species_moment_release(app, &bflux->gammai[2*d+e]);
      }
    }
  }
}
