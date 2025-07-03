#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_util.h>

static void
proj_on_basis_c2p_phase_func(const double *xcomp, double *xphys, void *ctx)
{
  struct gk_proj_on_basis_c2p_func_ctx *c2p_ctx = ctx;
  int cdim = c2p_ctx->cdim; // Assumes update range is a phase range.
  gkyl_position_map_eval_mc2nu(c2p_ctx->pos_map, xcomp, xphys);
  gkyl_velocity_map_eval_c2p(c2p_ctx->vel_map, &xcomp[cdim], &xphys[cdim]);
}

static void
proj_on_basis_c2p_position_func(const double *xcomp, double *xphys, void *ctx)
{
  struct gk_proj_on_basis_c2p_func_ctx *c2p_ctx = ctx;
  gkyl_position_map_eval_mc2nu(c2p_ctx->pos_map, xcomp, xphys);
}

struct func_gaussian_ctx {
  bool is_dir_periodic[GKYL_MAX_CDIM]; // Periodicity in configuration space.
  double box_size[GKYL_MAX_CDIM]; // Size of the box in each direction
  double gaussian_mean[GKYL_MAX_CDIM]; // Center in configuration space.
  double gaussian_std_dev[GKYL_MAX_CDIM]; // Sigma in configuration space, function is constant if sigma is 0.
  double f_floor; // Floor value of the distribution.
};

static void 
func_gaussian(double t, const double* xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct func_gaussian_ctx *inp = ctx;
  double envelope = 1.0;
  for (int dir = 0; dir < GKYL_MAX_CDIM; ++dir) {
    double dx = xn[dir] - inp->gaussian_mean[dir];
    double L = inp->box_size[dir];
    if (inp->is_dir_periodic[dir]) { 
      // Periodic wrapping
      dx = fmod(dx + L/2.0, L);
      if (dx < 0) dx += L;
      dx -= L/2.0;
    }
    if (inp->gaussian_std_dev[dir] > 0.0)
      envelope *= exp(-dx*dx/(2.0*inp->gaussian_std_dev[dir]*inp->gaussian_std_dev[dir]));
  }
  fout[0] = envelope + inp->f_floor;
}

struct proj_read_field_inp {
  struct gkyl_gyrokinetic_app *app; // Gyrokinetic app.
  struct gk_proj *proj; // Projection app.
  // Function (and context) defining the field in space.
  void (*mom_func)(double t, const double *xn, double *fout, void *ctx);
  void *mom_func_ctx;
  struct gkyl_proj_on_basis *proj_op; // Projection operator.
  struct gkyl_gyrokinetic_ic_import *ic_import; // Info for reading from file.
  struct gkyl_array *field_host; // Array on the host to store field.
};

static void
proj_read_field(double tm, struct proj_read_field_inp *inp)
{
  // Project a field or read it from file.
  struct gkyl_gyrokinetic_app *app = inp->app;

  if (inp->mom_func && inp->mom_func_ctx) {
    gkyl_proj_on_basis_advance(inp->proj_op, tm, &app->local, inp->field_host);
  }
  else if (inp->ic_import->type) {
    struct gkyl_app_restart_status rstat = gyrokinetic_header_from_file(inp->ic_import->file_name);

    if (rstat.io_status == GKYL_ARRAY_RIO_SUCCESS) {
      struct gkyl_app_restart_status rstat;
      rstat.io_status = gkyl_comm_array_read(app->comm, &app->grid, &app->local, inp->field_host, inp->ic_import->file_name);
    }
    else
      assert(false);
  }
  else
    assert(false);
}

static void
gk_species_projection_calc_proj_func(gkyl_gyrokinetic_app *app, struct gk_species *s, 
  struct gk_proj *proj, struct gkyl_array *f, double tm) 
{
  if (app->use_gpu) {
    gkyl_proj_on_basis_advance(proj->proj_func, tm, &s->local, proj->proj_host);
    gkyl_array_copy(f, proj->proj_host);
  }
  else {
    gkyl_proj_on_basis_advance(proj->proj_func, tm, &s->local, f);
  }

  // Multiply by the gyrocenter coord jacobian (bmag).
  gkyl_dg_mul_conf_phase_op_range(&app->basis, &s->basis, f, 
      app->gk_geom->bmag, f, &app->local, &s->local); 

  // Multiply by the velocity-space jacobian. 
  gkyl_array_scale_by_cell(f, s->vel_map->jacobvel);
}

static void
gk_species_projection_calc_max_prim(gkyl_gyrokinetic_app *app, struct gk_species *s, 
  struct gk_proj *proj, struct gkyl_array *f, double tm)
{
  proj_read_field(tm, &(struct proj_read_field_inp){
      .app = app,
      .mom_func = proj->info.density,
      .mom_func_ctx = proj->info.ctx_density,
      .proj_op = proj->proj_dens,
      .ic_import = &proj->info.density_from_file,
      .field_host = proj->dens,
    }
  );
  proj_read_field(tm, &(struct proj_read_field_inp){
      .app = app,
      .mom_func = proj->info.upar,
      .mom_func_ctx = proj->info.ctx_upar,
      .proj_op = proj->proj_upar,
      .ic_import = &proj->info.upar_from_file,
      .field_host = proj->upar,
    }
  );
  proj_read_field(tm, &(struct proj_read_field_inp){
      .app = app,
      .mom_func = proj->info.temp,
      .mom_func_ctx = proj->info.ctx_temp,
      .proj_op = proj->proj_temp,
      .ic_import = &proj->info.temp_from_file,
      .field_host = proj->vtsq,
    }
  );

  gkyl_array_scale(proj->vtsq, 1.0/s->info.mass);

  // Maxwellian projection expects primitive moments as one array (on GPU).
  gkyl_array_set_offset(proj->prim_moms_host, 1.0, proj->dens, 0*app->basis.num_basis);
  gkyl_array_set_offset(proj->prim_moms_host, 1.0, proj->upar, 1*app->basis.num_basis);
  gkyl_array_set_offset(proj->prim_moms_host, 1.0, proj->vtsq, 2*app->basis.num_basis);
  gkyl_array_copy(proj->prim_moms, proj->prim_moms_host);

  gkyl_gk_maxwellian_proj_on_basis_advance(proj->proj_max,
    &s->local, &app->local, proj->prim_moms, false, f);
}

static void
gk_species_projection_calc_bimax(gkyl_gyrokinetic_app *app, struct gk_species *s, 
  struct gk_proj *proj, struct gkyl_array *f, double tm)
{
  proj_read_field(tm, &(struct proj_read_field_inp){
      .app = app,
      .mom_func = proj->info.density,
      .mom_func_ctx = proj->info.ctx_density,
      .proj_op = proj->proj_dens,
      .ic_import = &proj->info.density_from_file,
      .field_host = proj->dens,
    }
  );
  proj_read_field(tm, &(struct proj_read_field_inp){
      .app = app,
      .mom_func = proj->info.upar,
      .mom_func_ctx = proj->info.ctx_upar,
      .proj_op = proj->proj_upar,
      .ic_import = &proj->info.upar_from_file,
      .field_host = proj->upar,
    }
  );
  proj_read_field(tm, &(struct proj_read_field_inp){
      .app = app,
      .mom_func = proj->info.temppar,
      .mom_func_ctx = proj->info.ctx_temppar,
      .proj_op = proj->proj_temppar,
      .ic_import = &proj->info.temppar_from_file,
      .field_host = proj->vtsqpar,
    }
  );
  proj_read_field(tm, &(struct proj_read_field_inp){
      .app = app,
      .mom_func = proj->info.tempperp,
      .mom_func_ctx = proj->info.ctx_tempperp,
      .proj_op = proj->proj_tempperp,
      .ic_import = &proj->info.tempperp_from_file,
      .field_host = proj->vtsqperp,
    }
  );

  gkyl_array_scale(proj->vtsqpar, 1.0/s->info.mass);
  gkyl_array_scale(proj->vtsqperp, 1.0/s->info.mass);

  // BiMaxwellian projection expects primitive moments as one array (on GPU).
  gkyl_array_set_offset(proj->prim_moms_host, 1.0, proj->dens, 0*app->basis.num_basis);
  gkyl_array_set_offset(proj->prim_moms_host, 1.0, proj->upar, 1*app->basis.num_basis);
  gkyl_array_set_offset(proj->prim_moms_host, 1.0, proj->vtsqpar , 2*app->basis.num_basis);
  gkyl_array_set_offset(proj->prim_moms_host, 1.0, proj->vtsqperp , 3*app->basis.num_basis);  
  gkyl_array_copy(proj->prim_moms, proj->prim_moms_host);

  gkyl_gk_maxwellian_proj_on_basis_advance(proj->proj_max,
    &s->local, &app->local, proj->prim_moms, false, f);
}

static void
gk_species_projection_calc_max_gauss(gkyl_gyrokinetic_app *app, struct gk_species *s, 
  struct gk_proj *proj, struct gkyl_array *f, double tm)
{
  gk_species_lte_from_moms(app, s, &s->lte, proj->prim_moms);
  gkyl_array_copy(f, s->lte.f_lte);
}

static void
gk_species_projection_calc_none(gkyl_gyrokinetic_app *app, struct gk_species *s, 
  struct gk_proj *proj, struct gkyl_array *f, double tm)
{
}

static void
gk_species_projection_correct_all_moms(gkyl_gyrokinetic_app *app, struct gk_species *s, 
  struct gk_proj *proj, struct gkyl_array *f, double tm)
{
  struct gkyl_gk_maxwellian_correct_status status_corr;
  status_corr = gkyl_gk_maxwellian_correct_all_moments(proj->corr_max, 
    f, proj->prim_moms, &s->local, &app->local);
}

static void
gk_species_projection_correct_all_moms_none(gkyl_gyrokinetic_app *app, struct gk_species *s, 
  struct gk_proj *proj, struct gkyl_array *f, double tm)
{
}

static struct gkyl_proj_on_basis *
alloc_mom_proj_op(struct proj_read_field_inp *inp)
{
  // Allocate operator to project a moment (if needed).
  struct gkyl_gyrokinetic_app *app = inp->app;

  struct gkyl_proj_on_basis *proj_out = 0;
  if (inp->mom_func && inp->mom_func_ctx) {
    proj_out = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
        .grid = &app->grid,
        .basis = &app->basis,
        .qtype = GKYL_GAUSS_QUAD,
        .num_quad = app->basis.poly_order+1,
        .num_ret_vals = 1,
        .eval = inp->mom_func,
        .ctx = inp->mom_func_ctx,
        .c2p_func = proj_on_basis_c2p_position_func,
        .c2p_func_ctx = &inp->proj->proj_on_basis_c2p_ctx,
      }
    );
  }

  return proj_out;
}

static void
init_maxwellian_bimaxwellian(struct gkyl_gyrokinetic_app *app, struct gk_species *s, 
  struct gkyl_gyrokinetic_projection inp, struct gk_proj *proj)
{
  proj->dens = mkarr(false, app->basis.num_basis, app->local_ext.volume);
  proj->proj_dens = alloc_mom_proj_op(&(struct proj_read_field_inp) {
      .app = app,
      .proj = proj,
      .mom_func = inp.density,
      .mom_func_ctx = inp.ctx_density,
    }
  );

  proj->upar = mkarr(false, app->basis.num_basis, app->local_ext.volume);
  proj->proj_upar = alloc_mom_proj_op(&(struct proj_read_field_inp) {
      .app = app,
      .proj = proj,
      .mom_func = inp.upar,
      .mom_func_ctx = inp.ctx_upar,
    }
  );

  bool bimaxwellian = false;
  if (proj->info.proj_id == GKYL_PROJ_MAXWELLIAN_PRIM) {
    proj->vtsq = mkarr(false, app->basis.num_basis, app->local_ext.volume);
    proj->proj_temp = alloc_mom_proj_op(&(struct proj_read_field_inp) {
        .app = app,
        .proj = proj,
        .mom_func = inp.temp,
        .mom_func_ctx = inp.ctx_temp,
      }
    );

    proj->prim_moms_host = mkarr(false, 3*app->basis.num_basis, app->local_ext.volume);
    proj->prim_moms = mkarr(app->use_gpu, 3*app->basis.num_basis, app->local_ext.volume);
  }
  else {
    bimaxwellian = true;
    proj->vtsqpar = mkarr(false, app->basis.num_basis, app->local_ext.volume);
    proj->proj_temppar = alloc_mom_proj_op(&(struct proj_read_field_inp) {
        .app = app,
        .proj = proj,
        .mom_func = inp.temppar,
        .mom_func_ctx = inp.ctx_temppar,
      }
    );

    proj->vtsqperp = mkarr(false, app->basis.num_basis, app->local_ext.volume);
    proj->proj_tempperp = alloc_mom_proj_op(&(struct proj_read_field_inp) {
        .app = app,
        .proj = proj,
        .mom_func = inp.tempperp,
        .mom_func_ctx = inp.ctx_tempperp,
      }
    );

    proj->prim_moms_host = mkarr(false, 4*app->basis.num_basis, app->local_ext.volume);
    proj->prim_moms = mkarr(app->use_gpu, 4*app->basis.num_basis, app->local_ext.volume);      
  }

  // Maxwellian (or BiMaxwellian) projection updater.
  struct gkyl_gk_maxwellian_proj_on_basis_inp inp_proj = {
    .phase_grid = &s->grid,
    .conf_basis = &app->basis,
    .phase_basis = &s->basis,
    .conf_range =  &app->local,
    .conf_range_ext = &app->local_ext,
    .vel_range = &s->local_vel, 
    .gk_geom = app->gk_geom,
    .vel_map = s->vel_map,
    .quad_type = inp.quad_type,
    .mass = s->info.mass,
    .bimaxwellian = bimaxwellian, 
    .divide_jacobgeo = false, // Jacobian multiplication handled in advance.
    .use_gpu = app->use_gpu,
  };
  proj->proj_max = gkyl_gk_maxwellian_proj_on_basis_inew( &inp_proj );

  proj->correct_all_moms = false; 
  if (inp.correct_all_moms) {
    proj->correct_all_moms = true;

    int max_iter = inp.max_iter > 0 ? inp.max_iter : 100;
    double iter_eps = inp.iter_eps > 0 ? inp.iter_eps  : 1e-12;
    bool use_last_converged = inp.use_last_converged; 

    // Maxwellian correction updater
    struct gkyl_gk_maxwellian_correct_inp inp_corr = {
      .phase_grid = &s->grid,
      .conf_basis = &app->basis,
      .phase_basis = &s->basis,
      .conf_range =  &app->local,
      .conf_range_ext = &app->local_ext,
      .vel_range = &s->local_vel, 
      .gk_geom = app->gk_geom,
      .vel_map = s->vel_map,
      .quad_type = inp.quad_type,
      .mass = s->info.mass,
      .bimaxwellian = bimaxwellian, 
      .divide_jacobgeo = false, // final Jacobian multiplication will be handled in advance
      .max_iter = max_iter,
      .eps = iter_eps,
      .use_last_converged = use_last_converged, 
      .use_gpu = app->use_gpu,
    };
    proj->corr_max = gkyl_gk_maxwellian_correct_inew( &inp_corr );
  }

}

static void
init_maxwellian_gaussian(struct gkyl_gyrokinetic_app *app, struct gk_species *s, 
  struct gkyl_gyrokinetic_projection inp, struct gk_proj *proj)
{
  // Fill the box_size attribute of the projection (used for periodicity).
  struct func_gaussian_ctx fg_ctx;
  fg_ctx.f_floor = inp.f_floor;
  for (int dir = 0; dir < app->cdim; ++dir) {
    fg_ctx.gaussian_mean[dir] = inp.gaussian_mean[dir];
    fg_ctx.gaussian_std_dev[dir] = inp.gaussian_std_dev[dir];
    fg_ctx.box_size[dir] = app->grid.upper[dir] - app->grid.lower[dir];
  }
  // Set periodicity for last dim if we are in IWL, and all other directions defined by the user.
  for (int dir = 0; dir < GKYL_MAX_CDIM; ++dir)
    fg_ctx.is_dir_periodic[dir] = false;
  fg_ctx.is_dir_periodic[app->cdim-1] = app->field->info.gkfield_id == GKYL_GK_FIELD_ES_IWL;
  for (int i=0; i < app->num_periodic_dir; ++i)
    fg_ctx.is_dir_periodic[app->periodic_dirs[i]] = true;

  struct gkyl_array *shape_ho = mkarr(false, app->basis.num_basis, app->local_ext.volume);
  struct gkyl_proj_on_basis *proj_gaussian = gkyl_proj_on_basis_new(&app->grid, &app->basis, app->poly_order + 1, 1, func_gaussian, &fg_ctx);
  proj->gaussian_profile = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);

  gkyl_proj_on_basis_advance(proj_gaussian, 0, &app->local, shape_ho);
  gkyl_array_copy(proj->gaussian_profile, shape_ho);

  gkyl_proj_on_basis_release(proj_gaussian);
  gkyl_array_release(shape_ho);

  // Build the integrant Jacobian * s(x), to normalize the shape function, and integrate it.
  struct gkyl_array *integrant = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  double *integral = app->use_gpu? gkyl_cu_malloc(sizeof(double)) : gkyl_malloc(sizeof(double));
  double *red_integral = app->use_gpu? gkyl_cu_malloc(sizeof(double)) : gkyl_malloc(sizeof(double));
  struct gkyl_array_integrate *int_op = gkyl_array_integrate_new(&app->grid, &app->basis, 1, 
    GKYL_ARRAY_INTEGRATE_OP_NONE, app->use_gpu);
  double red_integral_ho[1];

  gkyl_dg_mul_op_range(app->basis, 0, integrant, 0, app->gk_geom->jacobgeo, 0, proj->gaussian_profile, &app->local);
  gkyl_array_integrate_advance(int_op, integrant, 1.0, NULL, &app->local, NULL, integral);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, integral, red_integral);
  
  if (app->use_gpu) {
    gkyl_cu_memcpy(red_integral_ho, red_integral, sizeof(double), GKYL_CU_MEMCPY_D2H);
    gkyl_cu_free(integral);
    gkyl_cu_free(red_integral);
  } else {
    memcpy(red_integral_ho, red_integral, sizeof(double));
    gkyl_free(integral);
    gkyl_free(red_integral);
  }
  gkyl_array_release(integrant);
  gkyl_array_integrate_release(int_op);

  // Scale the shape configuration function
  gkyl_array_scale(proj->gaussian_profile, 1.0/red_integral_ho[0]);
  
  // We can now build the moments of the projection.
  proj->prim_moms = mkarr(app->use_gpu, 4*app->basis.num_basis, app->local_ext.volume);      

  // Density
  gkyl_array_set_offset(proj->prim_moms, inp.total_num_particles, proj->gaussian_profile, 0*app->basis.num_basis);

  // Parallel velocity
  gkyl_array_set_offset(proj->prim_moms, 0.0, proj->gaussian_profile, 1*app->basis.num_basis);
  
  // Temperature
  assert(inp.temp_max > 0);
  double temp = inp.total_num_particles == 0 ? inp.temp_max/2.0 : 2./3. * inp.total_kin_energy/inp.total_num_particles;
  temp = temp > inp.temp_max ? inp.temp_max : temp; // saturate to max temperature.
  gkyl_array_shiftc(proj->prim_moms, temp/s->info.mass, 2*app->basis.num_basis);
}

void 
gk_species_projection_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s, 
  struct gkyl_gyrokinetic_projection inp, struct gk_proj *proj)
{
  proj->info = inp;
  // Context for c2p function passed to proj_on_basis.
  proj->proj_on_basis_c2p_ctx.cdim = app->cdim;
  proj->proj_on_basis_c2p_ctx.vdim = s->local_vel.ndim;
  proj->proj_on_basis_c2p_ctx.vel_map = s->vel_map;
  proj->proj_on_basis_c2p_ctx.pos_map = app->position_map;

  if (proj->info.proj_id == GKYL_PROJ_FUNC) {
    proj->proj_func = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
        .grid = &s->grid,
        .basis = &s->basis,
        .qtype = GKYL_GAUSS_QUAD,
        .num_quad = app->basis.poly_order+1,
        .num_ret_vals = 1,
        .eval = inp.func,
        .ctx = inp.ctx_func,
        .c2p_func = proj_on_basis_c2p_phase_func,
        .c2p_func_ctx = &proj->proj_on_basis_c2p_ctx,
      }
    );
    if (app->use_gpu) {
      proj->proj_host = mkarr(false, s->basis.num_basis, s->local_ext.volume);
    }
    proj->projection_calc = gk_species_projection_calc_proj_func;
  }
  else if (proj->info.proj_id == GKYL_PROJ_MAXWELLIAN_PRIM || proj->info.proj_id == GKYL_PROJ_BIMAXWELLIAN) {
    init_maxwellian_bimaxwellian(app, s, inp, proj);
    proj->projection_calc = proj->info.proj_id == GKYL_PROJ_MAXWELLIAN_PRIM ?
      gk_species_projection_calc_max_prim : gk_species_projection_calc_bimax;
  } else if (proj->info.proj_id == GKYL_PROJ_MAXWELLIAN_GAUSSIAN) {
    init_maxwellian_gaussian(app, s, inp, proj);
    proj->projection_calc = gk_species_projection_calc_max_gauss;
  }

  proj->moms_correct = proj->correct_all_moms ? 
    gk_species_projection_correct_all_moms : gk_species_projection_correct_all_moms_none;
}

void
gk_species_projection_calc(gkyl_gyrokinetic_app *app, struct gk_species *s, 
  struct gk_proj *proj, struct gkyl_array *f, double tm)
{
  proj->projection_calc(app, s, proj, f, tm);
  proj->moms_correct(app, s, proj, f, tm);  
  // Multiply by the configuration space jacobian.
  gkyl_dg_mul_conf_phase_op_range(&app->basis, &s->basis, f, 
    app->gk_geom->jacobgeo, f, &app->local, &s->local);
}

void
gk_species_projection_release(const struct gkyl_gyrokinetic_app *app, const struct gk_proj *proj)
{
  if (proj->info.proj_id == GKYL_PROJ_FUNC) {
    gkyl_proj_on_basis_release(proj->proj_func);
    if (app->use_gpu) {
      gkyl_array_release(proj->proj_host);
    }
  }
  else if (proj->info.proj_id == GKYL_PROJ_MAXWELLIAN_PRIM || proj->info.proj_id == GKYL_PROJ_BIMAXWELLIAN) { 
    gkyl_array_release(proj->dens);
    gkyl_array_release(proj->upar); 
    gkyl_array_release(proj->prim_moms_host);
    gkyl_array_release(proj->prim_moms);
    if (proj->info.density && proj->info.ctx_density) 
      gkyl_proj_on_basis_release(proj->proj_dens);

    if (proj->info.upar && proj->info.ctx_upar) 
      gkyl_proj_on_basis_release(proj->proj_upar);

    gkyl_gk_maxwellian_proj_on_basis_release(proj->proj_max);

    if (proj->info.proj_id == GKYL_PROJ_MAXWELLIAN_PRIM) {
      gkyl_array_release(proj->vtsq);
      if (proj->info.temp && proj->info.ctx_temp) 
        gkyl_proj_on_basis_release(proj->proj_temp);
    }
    else if (proj->info.proj_id == GKYL_PROJ_BIMAXWELLIAN) {
      gkyl_array_release(proj->vtsqpar);
      gkyl_array_release(proj->vtsqperp);
      if (proj->info.temppar && proj->info.ctx_temppar)
        gkyl_proj_on_basis_release(proj->proj_temppar);

      if (proj->info.tempperp && proj->info.ctx_tempperp)
        gkyl_proj_on_basis_release(proj->proj_tempperp);
    }

    if (proj->correct_all_moms) {
      gkyl_gk_maxwellian_correct_release(proj->corr_max);    
    }
  } 
  else if (proj->info.proj_id == GKYL_PROJ_MAXWELLIAN_GAUSSIAN) {
    gkyl_array_release(proj->gaussian_profile);
    gkyl_array_release(proj->prim_moms);
  }
}
