#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

static void
gk_species_bflux_clear_dynamic(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **fin, double val)
{
  for (int d=0; d<app->cdim; ++d) {
    for (int e=0; e<2; ++e) {
      gkyl_array_clear(fin[2*d+e], val);
    }
  }
}

static void
gk_species_bflux_clear_static(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **fin, double val)
{
}

void
gk_species_bflux_clear(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **fin, double val)
{
  bflux->bflux_clear_func(app, bflux, fin, val);
}

static void
gk_species_bflux_scale_dynamic(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **fin, double val)
{
  for (int d=0; d<app->cdim; ++d) {
    for (int e=0; e<2; ++e) {
      gkyl_array_scale(fin[2*d+e], val);
    }
  }
}

static void
gk_species_bflux_scale_static(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **fin, double val)
{
}

void
gk_species_bflux_scale(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **fin, double val)
{
  bflux->bflux_scale_func(app, bflux, fin, val);
}

static void
gk_species_bflux_step_f_dynamic(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **fout, double dt, const struct gkyl_array **fin)
{
  for (int d=0; d<app->cdim; ++d) {
    for (int e=0; e<2; ++e) {
      gkyl_array_accumulate(gkyl_array_scale(fout[2*d+e], dt), 1.0, fin[2*d+e]);
    }
  }
}

static void
gk_species_bflux_step_f_static(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **fout, double dt, const struct gkyl_array **fin)
{
}

void
gk_species_bflux_step_f(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **fout, double dt, const struct gkyl_array **fin)
{
  bflux->bflux_step_f_func(app, bflux, fout, dt, fin);
}

static void
gk_species_bflux_combine_range_dynamic(gkyl_gyrokinetic_app *app, struct gk_species *gk_s, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **fout, double fac1, struct gkyl_array **fin1, double fac2, struct gkyl_array **fin2, struct gkyl_range *range)
{
  for (int d=0; d<app->cdim; ++d) {
    for (int e=0; e<2; ++e) {
      gkyl_array_accumulate_range(gkyl_array_set_range(fout[2*d+e], fac1, fin1[2*d+e], range), fac2, fin2[2*d+e], range);
    }
  }
}

static void
gk_species_bflux_combine_range_static(gkyl_gyrokinetic_app *app, struct gk_species *gk_s, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **fout, double fac1, struct gkyl_array **fin1, double fac2, struct gkyl_array **fin2, struct gkyl_range *range)
{
}

void
gk_species_bflux_combine_range(gkyl_gyrokinetic_app *app, struct gk_species *gk_s, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **fout, double fac1, struct gkyl_array **fin1, double fac2, struct gkyl_array **fin2, struct gkyl_range *range)
{
  bflux->bflux_combine_range_func(app, gk_s, bflux, fout, fac1, fin1, fac2, fin2, range);
}

static void
gk_species_bflux_copy_range_dynamic(gkyl_gyrokinetic_app *app, struct gk_species *gk_s, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **fout, struct gkyl_array **fin, struct gkyl_range *range)
{
  for (int d=0; d<app->cdim; ++d) {
    for (int e=0; e<2; ++e) {
      gkyl_array_copy_range(fout[2*d+e], fin[2*d+e], range);
    }
  }
}

static void
gk_species_bflux_copy_range_static(gkyl_gyrokinetic_app *app, struct gk_species *gk_s, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **fout, struct gkyl_array **fin, struct gkyl_range *range)
{
}

void
gk_species_bflux_copy_range(gkyl_gyrokinetic_app *app, struct gk_species *gk_s, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **fout, struct gkyl_array **fin, struct gkyl_range *range)
{
  bflux->bflux_copy_range_func(app, gk_s, bflux, fout, fin, range);
}

static void
gk_species_bflux_rhs_solver(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_boundary_fluxes *bflux, const struct gkyl_array *fin, struct gkyl_array *rhs,
  struct gkyl_array **bflux_moms)
{
  // Computes rhs of the boundary flux to be used by another solver.

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

static void
gk_species_bflux_rhs_diag(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_boundary_fluxes *bflux, const struct gkyl_array *fin, struct gkyl_array *rhs,
  struct gkyl_array **bflux_moms)
{
  // Computes rhs of the boundary flux for diagnostics.
  for (int d=0; d<app->cdim; ++d) {
    gkyl_array_clear_range(rhs, 0.0, &species->lower_ghost[d]);
    gkyl_array_clear_range(rhs, 0.0, &species->upper_ghost[d]);
  }

  // Only calculating density for use in ambipotential solve.
  for (int d=0; d<app->cdim; ++d) {
    for (int e=0; e<2; ++e) {
      if (species->info.integrated_hamiltonian_moments) {
        // Apply BC to phi so it is defined in the ghost cell.
        // Fill the ghost with the skin evaluated at the boundary.
        gkyl_bc_basic_advance(app->gfss_bc_op[2*d+e], app->bc_buffer, app->field->phi_smooth);

        // Scale phi by 0.5 in the ghost cell so that the Hamiltonian moment
        // becomes the 0.5*m*v^2+0.5*q*phi moment. Only this way does the energy
        // balance come out reasonable.
        gkyl_array_scale_range(app->field->phi_smooth, 0.5, e==0? &app->lower_ghost[d] : &app->upper_ghost[d]);
      }

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

static void
gk_species_bflux_rhs_static(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_boundary_fluxes *bflux, const struct gkyl_array *fin, struct gkyl_array *rhs,
  struct gkyl_array **bflux_moms)
{
}

static void
gk_species_bflux_rhs_solver_init(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_boundary_fluxes *bflux, const struct gkyl_array *fin, struct gkyl_array *rhs,
  struct gkyl_array **bflux_moms)
{
  gk_species_bflux_rhs_solver(app, species, bflux, fin, rhs, bflux_moms);

  // Only need to compute the initial boundary fluxes for static species; don't
  // need to compute them in time.
  if ( (species->info.is_static) || 
       !(app->field->update_field && app->field->gkfield_id == GKYL_GK_FIELD_BOLTZMANN) )
    bflux->bflux_rhs_func = gk_species_bflux_rhs_static; 
}

void
gk_species_bflux_rhs(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_boundary_fluxes *bflux, const struct gkyl_array *fin, struct gkyl_array *rhs,
  struct gkyl_array **bflux_moms)
{
  bflux->bflux_rhs_func(app, species, bflux, fin, rhs, bflux_moms);
}

void
gk_species_bflux_calc_boundary_flux_integrated_mom_dynamic(gkyl_gyrokinetic_app* app,
  const struct gk_species *gk_s, struct gk_boundary_fluxes *bflux, double tm)
{
  int vdim = app->vdim;
  int num_mom = bflux->moms_op.num_mom; 
  double avals_global[num_mom];

  for (int d=0; d<app->cdim; ++d) {
    for (int e=0; e<2; ++e) {
      // Integrated moment of the boundary flux.
      gkyl_array_integrate_advance(bflux->integ_op, bflux->f[2*d+e], 1., 0,
        e==0? &app->lower_ghost[d] : &app->upper_ghost[d], 0, bflux->int_moms_local);

      gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, num_mom, 
        bflux->int_moms_local, bflux->int_moms_global);
      if (app->use_gpu) {
        gkyl_cu_memcpy(avals_global, bflux->int_moms_global, sizeof(double[num_mom]), GKYL_CU_MEMCPY_D2H);
      }
      else {
        memcpy(avals_global, bflux->int_moms_global, sizeof(double[num_mom]));
      }

      gkyl_dynvec_append(bflux->intmom[2*d+e], tm, avals_global);
    }
  } 
}

static void
gk_species_bflux_calc_boundary_flux_integrated_mom_time_integrated_dynamic(gkyl_gyrokinetic_app* app,
  const struct gk_species *gk_s, struct gk_boundary_fluxes *bflux, double tm)
{
  int vdim = app->vdim;
  int num_mom = bflux->moms_op.num_mom; 
  double avals_global[num_mom];

  for (int d=0; d<app->cdim; ++d) {
    for (int e=0; e<2; ++e) {
      // Integrated moment of the boundary flux.
      gkyl_array_integrate_advance(bflux->integ_op, bflux->f[2*d+e], 1., 0,
        e==0? &app->lower_ghost[d] : &app->upper_ghost[d], 0, bflux->int_moms_local);

      gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, num_mom, 
        bflux->int_moms_local, bflux->int_moms_global);
      if (app->use_gpu) {
        gkyl_cu_memcpy(avals_global, bflux->int_moms_global, sizeof(double[num_mom]), GKYL_CU_MEMCPY_D2H);
      }
      else {
        memcpy(avals_global, bflux->int_moms_global, sizeof(double[num_mom]));
      }

      for (int k=0; k<num_mom; k++) {
        bflux->intmom_cumm_buff[(2*d+e)*num_mom+k] += avals_global[k];
      }
    }
  } 
}

static void
gk_species_bflux_calc_boundary_flux_integrated_mom_time_integrated_static(gkyl_gyrokinetic_app* app,
  const struct gk_species *gk_s, struct gk_boundary_fluxes *bflux, double tm)
{
}

void
gk_species_bflux_calc_boundary_flux_integrated_mom_time_integrated(gkyl_gyrokinetic_app* app,
  const struct gk_species *gk_s, struct gk_boundary_fluxes *bflux, double tm)
{
  bflux->calc_bflux_int_mom_time_integrate_func(app, gk_s, bflux, tm);
}

static void
gk_species_bflux_append_boundary_flux_integrated_mom(gkyl_gyrokinetic_app* app,
  const struct gk_species *gk_s, struct gk_boundary_fluxes *bflux, double tm)
{
  int num_mom = bflux->moms_op.num_mom; 

  for (int d=0; d<app->cdim; ++d) {
    for (int e=0; e<2; ++e) {
      // Append the time integrated moment of the boundary flux.
      gkyl_dynvec_append(bflux->intmom[2*d+e], tm, &bflux->intmom_cumm_buff[(2*d+e)*num_mom]);
    }
  } 
}

static void
gk_species_bflux_calc_boundary_flux_integrated_mom_static(gkyl_gyrokinetic_app* app,
  const struct gk_species *gk_s, struct gk_boundary_fluxes *bflux, double tm)
{
}

void
gk_species_bflux_calc_boundary_flux_integrated_mom(gkyl_gyrokinetic_app* app,
  const struct gk_species *gk_s, struct gk_boundary_fluxes *bflux, double tm)
{
  bflux->calc_boundary_flux_integrated_mom_func(app, gk_s, bflux, tm);
}
  
static void
gk_species_bflux_write_boundary_flux_integrated_mom_dynamic(gkyl_gyrokinetic_app *app,
  const struct gk_species *gks, struct gk_boundary_fluxes *bflux)
{
  int rank;
  gkyl_comm_get_rank(app->comm, &rank);

  const char *vars[] = {"x","y","z"};
  const char *edge[] = {"lower","upper"};

  if (rank == 0) {
    for (int d=0; d<app->cdim; ++d) {
      for (int e=0; e<2; ++e) {
        // Write integrated moments of the boundary fluxes.
        const char *fmt = "%s-%s_bflux_%s%s_%s.gkyl";

        int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, vars[d], edge[e], "integrated_moms");
        char fileNm[sz+1]; // ensures no buffer overflow
        snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, vars[d], edge[e], "integrated_moms");

        struct timespec wtm = gkyl_wall_clock();
        if (bflux->is_first_intmom_write_call) {
          gkyl_dynvec_write(bflux->intmom[2*d+e], fileNm);
        }
        else {
          gkyl_dynvec_awrite(bflux->intmom[2*d+e], fileNm);
        }

        app->stat.diag_io_tm += gkyl_time_diff_now_sec(wtm);
        app->stat.n_diag_io += 1;

        gkyl_dynvec_clear(bflux->intmom[2*d+e]);
      }
    }
    if (bflux->is_first_intmom_write_call)
      bflux->is_first_intmom_write_call = false;
  }
}

static void
gk_species_bflux_write_boundary_flux_integrated_mom_static(gkyl_gyrokinetic_app *app,
  const struct gk_species *gks, struct gk_boundary_fluxes *bflux)
{
}

void
gk_species_bflux_write_boundary_flux_integrated_mom(gkyl_gyrokinetic_app *app,
  const struct gk_species *gks, struct gk_boundary_fluxes *bflux)
{
  bflux->write_boundary_flux_integrated_mom_func(app, gks, bflux);
}

void 
gk_species_bflux_init(struct gkyl_gyrokinetic_app *app, struct gk_species *gk_s, struct gk_boundary_fluxes *bflux, bool is_diagnostic)
{ 
  // Set methods for time-stepping boundary fluxes.
  bflux->bflux_clear_func = gk_species_bflux_clear_static;
  bflux->bflux_scale_func = gk_species_bflux_scale_static;
  bflux->bflux_step_f_func = gk_species_bflux_step_f_static;
  bflux->bflux_combine_range_func = gk_species_bflux_combine_range_static;
  bflux->bflux_copy_range_func = gk_species_bflux_copy_range_static;
  bflux->bflux_rhs_func = gk_species_bflux_rhs_static; 
  bflux->calc_boundary_flux_integrated_mom_func = gk_species_bflux_calc_boundary_flux_integrated_mom_static;
  bflux->write_boundary_flux_integrated_mom_func = gk_species_bflux_write_boundary_flux_integrated_mom_static;
  bflux->calc_bflux_int_mom_time_integrate_func = gk_species_bflux_calc_boundary_flux_integrated_mom_time_integrated_static;
  bflux->allocated_diagnostic = false;
  bflux->allocated_solver = false;

  if ( is_diagnostic &&
       (gk_s->info.boundary_flux_diagnostics.num_diag_moments > 0 ||
        gk_s->info.boundary_flux_diagnostics.num_integrated_diag_moments > 0) ) {
    bflux->allocated_diagnostic = true;

    bflux->bflux_clear_func = gk_species_bflux_clear_dynamic;
    bflux->bflux_scale_func = gk_species_bflux_scale_dynamic;
    bflux->bflux_step_f_func = gk_species_bflux_step_f_dynamic;
    bflux->bflux_combine_range_func = gk_species_bflux_combine_range_dynamic;
    bflux->bflux_copy_range_func = gk_species_bflux_copy_range_dynamic;
    bflux->bflux_rhs_func = gk_species_bflux_rhs_diag; 
    if (gk_s->info.boundary_flux_diagnostics.time_integrated) {
      bflux->calc_boundary_flux_integrated_mom_func = gk_species_bflux_append_boundary_flux_integrated_mom;
      bflux->calc_bflux_int_mom_time_integrate_func = gk_species_bflux_calc_boundary_flux_integrated_mom_time_integrated_dynamic;
    }
    else
      bflux->calc_boundary_flux_integrated_mom_func = gk_species_bflux_calc_boundary_flux_integrated_mom_dynamic;
    bflux->write_boundary_flux_integrated_mom_func = gk_species_bflux_write_boundary_flux_integrated_mom_dynamic;

    // Object computing moments of the boundary flux.
    assert(gk_s->info.boundary_flux_diagnostics.num_diag_moments == 0); // Moments NYI.
    assert(gk_s->info.boundary_flux_diagnostics.num_integrated_diag_moments == 1); // 1 int moment allowed now.
    gk_species_moment_init(app, gk_s, &bflux->moms_op,
      gk_s->info.boundary_flux_diagnostics.integrated_diag_moments[0], false);
  
    int num_mom = bflux->moms_op.num_mom;
  
    bflux->f = gkyl_malloc((2*app->cdim)*sizeof(struct gkyl_array *));
    bflux->f1 = gkyl_malloc((2*app->cdim)*sizeof(struct gkyl_array *));
    bflux->fnew = gkyl_malloc((2*app->cdim)*sizeof(struct gkyl_array *));
    for (int d=0; d<app->cdim; ++d) {
      for (int e=0; e<2; ++e) {
        // Allocate arrays storing moments of the boundary flux.
        bflux->f[2*d+e] = mkarr(app->use_gpu, num_mom*app->confBasis.num_basis, app->local_ext.volume);
        bflux->f1[2*d+e] = mkarr(app->use_gpu, num_mom*app->confBasis.num_basis, app->local_ext.volume);
        bflux->fnew[2*d+e] = mkarr(app->use_gpu, num_mom*app->confBasis.num_basis, app->local_ext.volume);
  
        // Allocate solver.
        struct gkyl_range *skin_r = e==0? &gk_s->lower_skin[d] : &gk_s->upper_skin[d];
        struct gkyl_range *ghost_r = e==0? &gk_s->lower_ghost[d] : &gk_s->upper_ghost[d];
        bflux->flux_slvr[2*d+e] = gkyl_boundary_flux_new(d, e, &gk_s->grid,
          skin_r, ghost_r, gk_s->eqn_gyrokinetic, true, app->use_gpu);
      }
    }
  
    // Updater to compute the volume integral of the boundary flux moments.
    bflux->integ_op = gkyl_array_integrate_new(&app->grid, &app->confBasis,
      num_mom, GKYL_ARRAY_INTEGRATE_OP_NONE, app->use_gpu);
    if (app->use_gpu) {
      bflux->int_moms_local = gkyl_cu_malloc(num_mom*sizeof(double));
      bflux->int_moms_global = gkyl_cu_malloc(num_mom*sizeof(double));
    }
    else {
      bflux->int_moms_local = gkyl_malloc(num_mom*sizeof(double));
      bflux->int_moms_global = gkyl_malloc(num_mom*sizeof(double));
    }

    for (int d=0; d<app->cdim; ++d) {
      for (int e=0; e<2; ++e)
        bflux->intmom[2*d+e] = gkyl_dynvec_new(GKYL_DOUBLE, num_mom);
    }

    // Cummulative integrated moments of boundary fluxes.
    bflux->intmom_cumm_buff = gkyl_malloc(2*app->cdim*num_mom*sizeof(double));
    for (int d=0; d<app->cdim; ++d) {
      for (int e=0; e<2; ++e) {
        for (int k=0; k<num_mom; k++)
          bflux->intmom_cumm_buff[(2*d+e)*num_mom+k] = 0.0;
      }
    }
  
    bflux->is_first_intmom_write_call = true;
  }

  if (!is_diagnostic) {
    bflux->allocated_solver = true;

    if (app->field->update_field && app->field->gkfield_id == GKYL_GK_FIELD_BOLTZMANN)
      bflux->bflux_rhs_func = gk_species_bflux_rhs_solver; 
    else
      bflux->bflux_rhs_func = gk_species_bflux_rhs_solver_init; 

    for (int d=0; d<app->cdim; ++d) {
      for (int e=0; e<2; ++e) {
        // Allocate solver.
        struct gkyl_range *skin_r = e==0? &gk_s->lower_skin[d] : &gk_s->upper_skin[d];
        struct gkyl_range *ghost_r = e==0? &gk_s->lower_ghost[d] : &gk_s->upper_ghost[d];
        bflux->flux_slvr[2*d+e] = gkyl_boundary_flux_new(d, e, &gk_s->grid,
          skin_r, ghost_r, gk_s->eqn_gyrokinetic, false, app->use_gpu);

        // M0 moment calculator to get particle flux.
        gk_species_moment_init(app, gk_s, &bflux->gammai[2*d+e], "M0", false);
      }
    }
  }
}

void
gk_species_bflux_release(const struct gkyl_gyrokinetic_app *app, const struct gk_boundary_fluxes *bflux)
{
  if (bflux->allocated_diagnostic) {
    gk_species_moment_release(app, &bflux->moms_op); 
    gkyl_array_integrate_release(bflux->integ_op);
  
    for (int d=0; d<app->cdim; ++d) {
      for (int e=0; e<2; ++e) {
        gkyl_array_release(bflux->f[2*d+e]);
        gkyl_array_release(bflux->f1[2*d+e]);
        gkyl_array_release(bflux->fnew[2*d+e]);

        gkyl_boundary_flux_release(bflux->flux_slvr[2*d+e]);
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
    gkyl_free(bflux->intmom_cumm_buff);
  }

  if (bflux->allocated_solver) {
    for (int d=0; d<app->cdim; ++d) {
      for (int e=0; e<2; ++e) {
        gkyl_boundary_flux_release(bflux->flux_slvr[2*d+e]);
        gk_species_moment_release(app, &bflux->gammai[2*d+e]);
      }
    }
  }
}
