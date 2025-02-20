#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

static void
gk_species_bflux_clear_dynamic(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **fin, double val)
{
  for (int b=0; b<bflux->num_boundaries; ++b)
    gkyl_array_clear_range(fin[b], val, bflux->boundaries_conf_ghost[b]);
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
  for (int b=0; b<bflux->num_boundaries; ++b)
    gkyl_array_scale_range(fin[b], val, bflux->boundaries_conf_ghost[b]);
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
  for (int b=0; b<bflux->num_boundaries; ++b)
    gkyl_array_accumulate_range(gkyl_array_scale_range(fout[b], dt, bflux->boundaries_conf_ghost[b]),
      1.0, fin[b], bflux->boundaries_conf_ghost[b]);
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
gk_species_bflux_combine_dynamic(gkyl_gyrokinetic_app *app, struct gk_species *gk_s, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **fout, double fac1, struct gkyl_array **fin1, double fac2, struct gkyl_array **fin2)
{
  for (int b=0; b<bflux->num_boundaries; ++b)
    gkyl_array_accumulate_range(gkyl_array_set_range(fout[b], fac1, fin1[b], bflux->boundaries_conf_ghost[b]),
      fac2, fin2[b], bflux->boundaries_conf_ghost[b]);
}

static void
gk_species_bflux_combine_static(gkyl_gyrokinetic_app *app, struct gk_species *gk_s, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **fout, double fac1, struct gkyl_array **fin1, double fac2, struct gkyl_array **fin2)
{
}

void
gk_species_bflux_combine(gkyl_gyrokinetic_app *app, struct gk_species *gk_s, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **fout, double fac1, struct gkyl_array **fin1, double fac2, struct gkyl_array **fin2)
{
  bflux->bflux_combine_func(app, gk_s, bflux, fout, fac1, fin1, fac2, fin2);
}

static void
gk_species_bflux_copy_dynamic(gkyl_gyrokinetic_app *app, struct gk_species *gk_s, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **fout, struct gkyl_array **fin)
{
  for (int b=0; b<bflux->num_boundaries; ++b)
    gkyl_array_copy_range(fout[b], fin[b], bflux->boundaries_conf_ghost[b]);
}

static void
gk_species_bflux_copy_static(gkyl_gyrokinetic_app *app, struct gk_species *gk_s, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **fout, struct gkyl_array **fin)
{
}

void
gk_species_bflux_copy(gkyl_gyrokinetic_app *app, struct gk_species *gk_s, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **fout, struct gkyl_array **fin)
{
  bflux->bflux_copy_func(app, gk_s, bflux, fout, fin);
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

  // Only calculating density for use in ambipotential solve.
  for (int b=0; b<bflux->num_boundaries; ++b) {
    const struct gkyl_range *phase_ghost_r, *conf_ghost_r;
    int dir = bflux->boundaries_dir[b];
    if (bflux->boundaries_edge[b]==GKYL_LOWER_EDGE) {
      phase_ghost_r = &species->lower_ghost[dir];
      conf_ghost_r = &app->lower_ghost[dir];
    }
    else {
      phase_ghost_r = &species->upper_ghost[dir];
      conf_ghost_r = &app->upper_ghost[dir];
    }

    gkyl_array_clear_range(rhs, 0.0, conf_ghost_r);

    if (species->info.integrated_hamiltonian_moments) {
      // Apply BC to phi so it is defined in the ghost cell.
      // Fill the ghost with the skin evaluated at the boundary.
      gkyl_bc_basic_advance(bflux->gfss_bc_op[b], bflux->bc_buffer, app->field->phi_smooth);

      // Scale phi by 0.5 in the ghost cell so that the Hamiltonian moment
      // becomes the 0.5*m*v^2+0.5*q*phi moment. Only this way does the energy
      // balance come out reasonable.
      gkyl_array_scale_range(app->field->phi_smooth, 0.5, conf_ghost_r);
    }

    // Ghost cells of the rhs array are filled with the bflux
    // This is overwritten by the boundary conditions and is not being stored,
    // it is only currently used to calculate moments for other applications.
    gkyl_boundary_flux_advance(bflux->flux_slvr[b], fin, rhs);

    // Compute moments of boundary fluxes.
    gk_species_moment_calc(&bflux->moms_op, *phase_ghost_r, *conf_ghost_r, rhs);

    gkyl_array_copy_range(bflux_moms[b], bflux->moms_op.marr, conf_ghost_r);
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
gk_species_bflux_calc_integrated_mom_dynamic(gkyl_gyrokinetic_app* app,
  const struct gk_species *gk_s, struct gk_boundary_fluxes *bflux, double tm)
{
  if (bflux->is_not_first_restart_write_call) {
    int vdim = app->vdim;
    int num_mom = bflux->moms_op.num_mom; 
    double avals_global[num_mom];

    for (int b=0; b<bflux->num_boundaries; ++b) {
      // Integrated moment of the boundary flux.
      int dir = bflux->boundaries_dir[b];
      gkyl_array_integrate_advance(bflux->integ_op, bflux->f[b], 1., 0,
        bflux->boundaries_edge[b]==GKYL_LOWER_EDGE? &app->lower_ghost[dir] : &app->upper_ghost[dir],
        0, bflux->int_moms_local);

      gkyl_comm_allreduce(app->comm_plane[dir], GKYL_DOUBLE, GKYL_SUM, num_mom, 
        bflux->int_moms_local, bflux->int_moms_global);
      if (app->use_gpu) {
        gkyl_cu_memcpy(avals_global, bflux->int_moms_global, sizeof(double[num_mom]), GKYL_CU_MEMCPY_D2H);
      }
      else {
        memcpy(avals_global, bflux->int_moms_global, sizeof(double[num_mom]));
      }

      gkyl_dynvec_append(bflux->intmom[b], tm, avals_global);
      
    } 
  } 
}

static void
gk_species_bflux_calc_voltime_integrated_mom_dynamic(gkyl_gyrokinetic_app* app,
  const struct gk_species *gk_s, struct gk_boundary_fluxes *bflux, double tm)
{
  if (bflux->is_not_first_restart_write_call) {
    int vdim = app->vdim;
    int num_mom = bflux->moms_op.num_mom; 
    double avals_global[num_mom];

    for (int b=0; b<bflux->num_boundaries; ++b) {
      // Integrated moment of the boundary flux.
      int dir = bflux->boundaries_dir[b];
      gkyl_array_integrate_advance(bflux->integ_op, bflux->f[b], 1., 0,
        bflux->boundaries_edge[b]==GKYL_LOWER_EDGE? &app->lower_ghost[dir] : &app->upper_ghost[dir], 0, bflux->int_moms_local);

      gkyl_comm_allreduce(app->comm_plane[dir], GKYL_DOUBLE, GKYL_SUM, num_mom, 
        bflux->int_moms_local, bflux->int_moms_global);
      if (app->use_gpu) {
        gkyl_cu_memcpy(avals_global, bflux->int_moms_global, sizeof(double[num_mom]), GKYL_CU_MEMCPY_D2H);
      }
      else {
        memcpy(avals_global, bflux->int_moms_global, sizeof(double[num_mom]));
      }

      for (int k=0; k<num_mom; k++) {
        bflux->intmom_cumm_buff[b*num_mom+k] += avals_global[k];
      }
    }
  }
}

static void
gk_species_bflux_calc_voltime_integrated_mom_static(gkyl_gyrokinetic_app* app,
  const struct gk_species *gk_s, struct gk_boundary_fluxes *bflux, double tm)
{
}

void
gk_species_bflux_calc_voltime_integrated_mom(gkyl_gyrokinetic_app* app,
  const struct gk_species *gk_s, struct gk_boundary_fluxes *bflux, double tm)
{
  bflux->bflux_calc_voltime_int_mom_func(app, gk_s, bflux, tm);
}

static void
gk_species_bflux_append_integrated_mom(gkyl_gyrokinetic_app* app,
  const struct gk_species *gk_s, struct gk_boundary_fluxes *bflux, double tm)
{
  int num_mom = bflux->moms_op.num_mom; 

  // Append the time integrated moment of the boundary flux.
  for (int b=0; b<bflux->num_boundaries; ++b)
    gkyl_dynvec_append(bflux->intmom[b], tm, &bflux->intmom_cumm_buff[b*num_mom]);
}

static void
gk_species_bflux_calc_integrated_mom_static(gkyl_gyrokinetic_app* app,
  const struct gk_species *gk_s, struct gk_boundary_fluxes *bflux, double tm)
{
}

void
gk_species_bflux_calc_integrated_mom(gkyl_gyrokinetic_app* app,
  const struct gk_species *gk_s, struct gk_boundary_fluxes *bflux, double tm)
{
  bflux->bflux_calc_integrated_mom_func(app, gk_s, bflux, tm);
}
  
static void
gk_species_bflux_write_integrated_mom_dynamic(gkyl_gyrokinetic_app *app,
  const struct gk_species *gks, struct gk_boundary_fluxes *bflux)
{
  int rank, comm_size;
  gkyl_comm_get_rank(app->comm, &rank);
  gkyl_comm_get_size(app->comm, &comm_size);

  const char *vars[] = {"x","y","z"};
  const char *edge[] = {"lower","upper"};

  for (int b=0; b<bflux->num_boundaries; ++b) {
    int dir = bflux->boundaries_dir[b];
    int edi = bflux->boundaries_edge[b]==GKYL_LOWER_EDGE? 0 : 1;

    if ((edi == 0 && rank == 0) || (edi == 1 && rank == comm_size-1)) {
      // Write integrated moments of the boundary fluxes.
      const char *fmt = "%s-%s_bflux_%s%s_%s.gkyl";

      int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, vars[dir], edge[edi], "integrated_moms");
      char fileNm[sz+1]; // ensures no buffer overflow
      snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, vars[dir], edge[edi], "integrated_moms");

      if (!bflux->is_not_first_restart_write_call) {
        // We didn't calculate int_mom at restart, so read the value from the previous sim
        // and append it. This is the best solution given the flow in present input files.
        bool res = gkyl_dynvec_read(bflux->intmom[b], fileNm);
        int num_mom = gkyl_dynvec_ncomp(bflux->intmom[b]);

        double vals_prev[num_mom];
        gkyl_dynvec_getlast(bflux->intmom[b], vals_prev);
        double tm = gkyl_dynvec_getlast_tm(bflux->intmom[b]);
        gkyl_dynvec_clear(bflux->intmom[b]);

        if (gks->info.boundary_flux_diagnostics.time_integrated) {
          for (int k=0; k<num_mom; k++)
            bflux->intmom_cumm_buff[b*num_mom+k] += vals_prev[k];
        }

        gkyl_dynvec_append(bflux->intmom[b], tm, vals_prev);
      }

      struct timespec wtm = gkyl_wall_clock();
      if (bflux->is_first_intmom_write_call[b]) {
        gkyl_dynvec_write(bflux->intmom[b], fileNm);
      }
      else {
        gkyl_dynvec_awrite(bflux->intmom[b], fileNm);
      }

      app->stat.diag_io_tm += gkyl_time_diff_now_sec(wtm);
      app->stat.n_diag_io += 1;

      gkyl_dynvec_clear(bflux->intmom[b]);
    
      if (bflux->is_first_intmom_write_call[b])
        bflux->is_first_intmom_write_call[b] = false;
    }
  }

  if (!bflux->is_not_first_restart_write_call)
    bflux->is_not_first_restart_write_call = true;
}

static void
gk_species_bflux_write_integrated_mom_static(gkyl_gyrokinetic_app *app,
  const struct gk_species *gks, struct gk_boundary_fluxes *bflux)
{
}

void
gk_species_bflux_write_integrated_mom(gkyl_gyrokinetic_app *app,
  const struct gk_species *gks, struct gk_boundary_fluxes *bflux)
{
  bflux->bflux_write_integrated_mom_func(app, gks, bflux);
}

void 
gk_species_bflux_init(struct gkyl_gyrokinetic_app *app, struct gk_species *gk_s, struct gk_boundary_fluxes *bflux, bool is_diagnostic)
{ 
  // Set methods for time-stepping boundary fluxes.
  bflux->bflux_clear_func = gk_species_bflux_clear_static;
  bflux->bflux_scale_func = gk_species_bflux_scale_static;
  bflux->bflux_step_f_func = gk_species_bflux_step_f_static;
  bflux->bflux_combine_func = gk_species_bflux_combine_static;
  bflux->bflux_copy_func = gk_species_bflux_copy_static;
  bflux->bflux_rhs_func = gk_species_bflux_rhs_static; 
  bflux->bflux_calc_integrated_mom_func = gk_species_bflux_calc_integrated_mom_static;
  bflux->bflux_write_integrated_mom_func = gk_species_bflux_write_integrated_mom_static;
  bflux->bflux_calc_voltime_int_mom_func = gk_species_bflux_calc_voltime_integrated_mom_static;
  bflux->allocated_diagnostic = false;
  bflux->allocated_solver = false;

  if ( is_diagnostic &&
       (gk_s->info.boundary_flux_diagnostics.num_diag_moments > 0 ||
        gk_s->info.boundary_flux_diagnostics.num_integrated_diag_moments > 0) ) {
    bflux->allocated_diagnostic = true;

    bflux->bflux_clear_func = gk_species_bflux_clear_dynamic;
    bflux->bflux_scale_func = gk_species_bflux_scale_dynamic;
    bflux->bflux_step_f_func = gk_species_bflux_step_f_dynamic;
    bflux->bflux_combine_func = gk_species_bflux_combine_dynamic;
    bflux->bflux_copy_func = gk_species_bflux_copy_dynamic;
    bflux->bflux_rhs_func = gk_species_bflux_rhs_diag; 
    if (gk_s->info.boundary_flux_diagnostics.time_integrated) {
      bflux->bflux_calc_integrated_mom_func = gk_species_bflux_append_integrated_mom;
      bflux->bflux_calc_voltime_int_mom_func = gk_species_bflux_calc_voltime_integrated_mom_dynamic;
    }
    else
      bflux->bflux_calc_integrated_mom_func = gk_species_bflux_calc_integrated_mom_dynamic;
    bflux->bflux_write_integrated_mom_func = gk_species_bflux_write_integrated_mom_dynamic;

    // Object computing moments of the boundary flux.
    assert(gk_s->info.boundary_flux_diagnostics.num_diag_moments == 0); // Moments NYI.
    assert(gk_s->info.boundary_flux_diagnostics.num_integrated_diag_moments == 1); // 1 int moment allowed now.
    gk_species_moment_init(app, gk_s, &bflux->moms_op,
      gk_s->info.boundary_flux_diagnostics.integrated_diag_moments[0], false);
  
    int num_mom = bflux->moms_op.num_mom;

    // Identify the non-periodic, non-zero-flux boundaries to compute
    // diagnostics at.
    bflux->num_boundaries = 0;
    for (int d=0; d<app->cdim; ++d) {
      if (gk_s->bc_is_np[d]) {
        for (int e=0; e<2; ++e) {
          if (gk_s->lower_bc[d].type != GKYL_SPECIES_ZERO_FLUX || gk_s->upper_bc[d].type != GKYL_SPECIES_ZERO_FLUX) {
            bflux->boundaries_dir[bflux->num_boundaries] = d;
            bflux->boundaries_edge[bflux->num_boundaries] = e==0? GKYL_LOWER_EDGE : GKYL_UPPER_EDGE;
            bflux->boundaries_conf_ghost[bflux->num_boundaries] = e==0? &app->lower_ghost[d] : &app->global_upper_ghost[d];
            bflux->num_boundaries++;
          }
        }
      }
    }
  
    // Ghost from skin surf (gfss): Fill the global ghost cells the field evaluated at the boundary.
    // Needed for bmag and phi to be nonzero in the ghost cell so we can take
    // the Hamiltonian moment of the boundaryb fluxes.
    long buff_sz = 0;
    for (int b=0; b<bflux->num_boundaries; ++b) {
      int dir = bflux->boundaries_dir[b];
      struct gkyl_range *skin_r = bflux->boundaries_edge[b]==GKYL_LOWER_EDGE? &app->global_lower_skin[dir]
                                                                            : &app->global_upper_skin[dir];
      struct gkyl_range *ghost_r = bflux->boundaries_edge[b]==GKYL_LOWER_EDGE? &app->global_lower_ghost[dir]
                                                                             : &app->global_upper_ghost[dir];

      bflux->gfss_bc_op[b] = gkyl_bc_basic_new(bflux->boundaries_dir[b], bflux->boundaries_edge[b], GKYL_BC_CONF_BOUNDARY_VALUE,
        &app->confBasis, skin_r, ghost_r, 1, app->cdim, app->use_gpu);
      
      long vol = skin_r->volume;
      buff_sz = buff_sz > vol ? buff_sz : vol;
    }
    bflux->bc_buffer = mkarr(app->use_gpu, app->confBasis.num_basis, buff_sz);

    // Fill ghost cell of bmag.
    for (int b=0; b<bflux->num_boundaries; ++b)
      gkyl_bc_basic_advance(bflux->gfss_bc_op[b], bflux->bc_buffer, app->gk_geom->bmag);

    bflux->f = gkyl_malloc(bflux->num_boundaries*sizeof(struct gkyl_array *));
    bflux->f1 = gkyl_malloc(bflux->num_boundaries*sizeof(struct gkyl_array *));
    bflux->fnew = gkyl_malloc(bflux->num_boundaries*sizeof(struct gkyl_array *));
    for (int b=0; b<bflux->num_boundaries; ++b) {
      // Allocate arrays storing moments of the boundary flux.
      bflux->f[b] = mkarr(app->use_gpu, num_mom*app->confBasis.num_basis, app->local_ext.volume);
      bflux->f1[b] = mkarr(app->use_gpu, num_mom*app->confBasis.num_basis, app->local_ext.volume);
      bflux->fnew[b] = mkarr(app->use_gpu, num_mom*app->confBasis.num_basis, app->local_ext.volume);
  
      // Allocate solver.
      int dir = bflux->boundaries_dir[b];
      struct gkyl_range *skin_r = bflux->boundaries_edge[b]==GKYL_LOWER_EDGE? &gk_s->lower_skin[dir] : &gk_s->upper_skin[dir];
      struct gkyl_range *ghost_r = bflux->boundaries_edge[b]==GKYL_LOWER_EDGE? &gk_s->lower_ghost[dir] : &gk_s->upper_ghost[dir];
      bflux->flux_slvr[b] = gkyl_boundary_flux_new(dir, bflux->boundaries_edge[b], &gk_s->grid,
        skin_r, ghost_r, gk_s->eqn_gyrokinetic, true, app->use_gpu);
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

    for (int b=0; b<bflux->num_boundaries; ++b) {
      bflux->intmom[b] = gkyl_dynvec_new(GKYL_DOUBLE, num_mom);
      bflux->is_first_intmom_write_call[b] = true;
    }
    bflux->is_not_first_restart_write_call = true;

    // Cummulative integrated moments of boundary fluxes.
    bflux->intmom_cumm_buff = gkyl_malloc(bflux->num_boundaries*num_mom*sizeof(double));
    for (int b=0; b<bflux->num_boundaries; ++b) {
      for (int k=0; k<num_mom; k++)
        bflux->intmom_cumm_buff[b*num_mom+k] = 0.0;
    }
  
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
    gkyl_array_release(bflux->bc_buffer);
    for (int b=0; b<bflux->num_boundaries; ++b)
      gkyl_bc_basic_release(bflux->gfss_bc_op[b]);

    gk_species_moment_release(app, &bflux->moms_op); 
    gkyl_array_integrate_release(bflux->integ_op);
  
    for (int b=0; b<bflux->num_boundaries; ++b) {
      gkyl_array_release(bflux->f[b]);
      gkyl_array_release(bflux->f1[b]);
      gkyl_array_release(bflux->fnew[b]);

      gkyl_boundary_flux_release(bflux->flux_slvr[b]);
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
    for (int b=0; b<bflux->num_boundaries; ++b) {
      gkyl_dynvec_release(bflux->intmom[b]);
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
