#include <assert.h>
#include <float.h>

#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_array_rio_priv.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_util.h>
#include <gkyl_gyrokinetic_priv.h>

static void
gk_fluid_neut_species_prim_var(gkyl_gyrokinetic_app *app, struct gk_fluid_neut_species *fns,
  const struct gkyl_array *moms)
{
  struct timespec tm = gkyl_wall_clock();

//  // Compute flow velocity in both the volume and on surfaces
//  gkyl_dg_calc_fluid_vars_advance(f->calc_fluid_vars_ext,
//    fluid, f->cell_avg_prim, f->u, f->u_surf); 
//
//  // Compute scalar pressure in the volume and at needed surfaces
//  gkyl_dg_calc_fluid_vars_pressure(f->calc_fluid_vars, 
//    &app->local_ext, fluid, f->u, f->p, f->p_surf);

//  app->stat.fluid_species_vars_tm += gkyl_time_diff_now_sec(tm);  
}

static void
gk_fluid_neut_species_write_static(gkyl_gyrokinetic_app *app, struct gk_fluid_neut_species *fns, 
  double tm, int frame)
{
  struct gkyl_msgpack_data *mt = gk_array_meta_new( (struct gyrokinetic_output_meta) {
      .frame = frame,
      .stime = tm,
      .poly_order = app->poly_order,
      .basis_type = app->basis.id
    }
  );

  const char *fmt = "%s-%s_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, fns->info.name, frame);
  char fileNm[sz+1]; // Ensures no buffer overflow.
  snprintf(fileNm, sizeof fileNm, fmt, app->name, fns->info.name, frame);

  // Copy data from device to host.
  if (app->use_gpu)
    gkyl_array_copy(fns->f_host, fns->f);

  gkyl_comm_array_write(app->comm, &app->grid, &app->local, 
    mt, fns->f_host, fileNm); 

  // Write primitive variables (udrift, pressure).
  fns->prim_var_func(app, fns, fns->f);
  const char *fmt_prim = "%s-%s_prim_var_%d.gkyl";
  int sz_prim = gkyl_calc_strlen(fmt_prim, app->name, fns->info.name, frame);
  char fileNm_prim[sz_prim+1]; // Ensures no buffer overflow.
  snprintf(fileNm_prim, sizeof fileNm_prim, fmt_prim, app->name, fns->info.name, frame);

  // Copy data from device to host.
  if (app->use_gpu)
    gkyl_array_copy(fns->prim_var_host, fns->prim_var);

  gkyl_comm_array_write(app->comm, &app->grid, &app->local, 
    mt, fns->prim_var_host, fileNm_prim); 

  gk_array_meta_release(mt);       
}

static void
gk_fluid_neut_species_release_static(const gkyl_gyrokinetic_app *app, struct gk_fluid_neut_species *fns)
{
  // Release solver and fluid arrays for update.
  gkyl_array_release(fns->f);
  gkyl_array_release(fns->f1);
  gkyl_array_release(fns->fnew);

  gkyl_array_release(fns->f_host);

  gkyl_array_release(fns->cflrate);
  if (app->use_gpu)
    gkyl_cu_free(fns->omega_cfl);
  else
    gkyl_free(fns->omega_cfl);

  gkyl_array_release(fns->prim_var);
  gkyl_array_release(fns->prim_var_host);

  gkyl_array_release(fns->integ_mom);
  gkyl_dynvec_release(fns->integ_diag);
  if (app->use_gpu)
    gkyl_cu_free(fns->red_integ_diag);
  else
    gkyl_free(fns->red_integ_diag);

}

static void
gk_fluid_neut_species_calc_integrated_mom_enabled(struct gkyl_gyrokinetic_app *app,
  struct gk_fluid_neut_species *fns, double tm)
{
  // Integrated moments: rho, rho*ux, rho*uy, rho*uy, flowE=0.5*rho*u^2, thermalE=p/(gas_gamma-1).
  fns->prim_var_func(app, fns, fns->f);

  gkyl_array_clear(fns->integ_mom, 0.0);
//  gkyl_dg_calc_fluid_integrated_vars(fns->calc_fluid_vars, &app->local, 
//    fns->fluid, fns->u, fns->p, fns->integ_mom);

  double avals_mom[fns->num_integ_mom], avals_mom_global[fns->num_integ_mom];
  gkyl_array_scale_range(fns->integ_mom, app->grid.cellVolume, &app->local);
  gkyl_array_reduce_range(fns->red_integ_diag, fns->integ_mom, GKYL_SUM, &app->local);
  if (app->use_gpu)
    gkyl_cu_memcpy(avals_mom, fns->red_integ_diag, sizeof(double[fns->num_integ_mom]), GKYL_CU_MEMCPY_D2H);
  else 
    memcpy(avals_mom, fns->red_integ_diag, sizeof(double[fns->num_integ_mom]));

  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_SUM, fns->num_integ_mom, avals_mom, avals_mom_global);
  gkyl_dynvec_append(fns->integ_diag, tm, avals_mom_global);
}

void
gk_fluid_neut_species_write_integrated_mom_enabled(gkyl_gyrokinetic_app *app,
  struct gk_fluid_neut_species *fns)
{
  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  if (rank == 0) {
    // Write out integrated diagnostic moments.
    const char *fmt = "%s-%s_%s.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, fns->info.name, "integrated_moms");
    char fileNm[sz+1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, app->name, fns->info.name, "integrated_moms");

    if (fns->is_first_integ_write_call) {
      gkyl_dynvec_write(fns->integ_diag, fileNm);
      fns->is_first_integ_write_call = false;
    }
    else {
      gkyl_dynvec_awrite(fns->integ_diag, fileNm);
    }
  }
  gkyl_dynvec_clear(fns->integ_diag);
}

// Initialize fluid species object.
void
gk_fluid_neut_species_init(struct gkyl_gk *gk, struct gkyl_gyrokinetic_app *app, struct gk_fluid_neut_species *fns)
{
  int cdim = app->cdim;

  // Number of moments depends on eqn_type.
  fns->num_moments = 5; // rho, rho*ux, rho*uy, rho*uy, totE

  // Moment arrays.
  fns->f = mkarr(app->use_gpu, fns->num_moments*app->basis.num_basis, app->local_ext.volume);
  fns->f1 = mkarr(app->use_gpu, fns->num_moments*app->basis.num_basis, app->local_ext.volume);
  fns->fnew = mkarr(app->use_gpu, fns->num_moments*app->basis.num_basis, app->local_ext.volume);

  fns->f_host = app->use_gpu? mkarr(false, fns->f->ncomp, fns->f->size)
                            : gkyl_array_acquire(fns->f);

  // CFL rate as a function of space, and (max) reduced rate.
  fns->cflrate = mkarr(app->use_gpu, 1, app->local_ext.volume);
  if (app->use_gpu) {
    fns->omega_cfl = gkyl_cu_malloc(sizeof(double));
  }
  else {
    fns->omega_cfl = gkyl_malloc(sizeof(double));
  }

  // Primitive variables (udrift, pressure).
  fns->num_prim_var = 4;
  fns->prim_var = mkarr(app->use_gpu, fns->num_prim_var*app->basis.num_basis, app->local_ext.volume);
  fns->prim_var_host = app->use_gpu? mkarr(false, fns->prim_var->ncomp, fns->prim_var->size)
                                   : gkyl_array_acquire(fns->prim_var);

  // Integrated moments: rho, rho*ux, rho*uy, rho*uy, flowE=0.5*rho*u^2, thermalE=p/(gas_gamma-1).
  fns->num_integ_mom = 6;
  fns->integ_mom = mkarr(app->use_gpu, fns->num_integ_mom, app->local_ext.volume);
  if (app->use_gpu) {
    fns->red_integ_diag = gkyl_cu_malloc(sizeof(double[fns->num_integ_mom]));
  }
  else {
    fns->red_integ_diag = gkyl_malloc(sizeof(double[fns->num_integ_mom]));
  }
  // Dynamic-vector for integrated moments.
  fns->integ_diag = gkyl_dynvec_new(GKYL_DOUBLE, fns->num_integ_mom);
  fns->is_first_integ_write_call = true;

  // Function pointers chosen at runtime.
  fns->prim_var_func = gk_fluid_neut_species_prim_var; 
  fns->write_func = gk_fluid_neut_species_write_static; 
  fns->calc_integrated_mom_func = gk_fluid_neut_species_calc_integrated_mom_enabled;  
  fns->write_integrated_mom_func = gk_fluid_neut_species_write_integrated_mom_enabled;  
  fns->release_func = gk_fluid_neut_species_release_static; 
}

void
gk_fluid_neut_species_write(gkyl_gyrokinetic_app *app, struct gk_fluid_neut_species *fns, 
  double tm, int frame)
{
  fns->write_func(app, fns, tm, frame);
}

void
gk_fluid_neut_species_calc_integrated_mom(struct gkyl_gyrokinetic_app *app,
  struct gk_fluid_neut_species *fns, double tm)
{
  fns->calc_integrated_mom_func(app, fns, tm);
}

void
gk_fluid_neut_species_write_integrated_mom(struct gkyl_gyrokinetic_app *app,
  struct gk_fluid_neut_species *fns)
{
  fns->write_integrated_mom_func(app, fns);
}

// Release resources for fluid species.
void
gk_fluid_neut_species_release(const gkyl_gyrokinetic_app* app, struct gk_fluid_neut_species *fns)
{
  fns->release_func(app, fns);
}
