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

// Default shaping function for Maxwellian Gaussian projection (to be normalized later).
void func_gaussian(double t, const double* xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gkyl_gyrokinetic_projection *inp = ctx;
  double envelope = 1.0;
  for (int dir = 0; dir < GKYL_MAX_CDIM; ++dir) {
    double t = xn[dir] - inp->center_gauss[dir];
    if (inp->periodic[dir])
      t = fmod(fabs(xn[dir]) - inp->center_gauss[dir], inp->box_size[dir]);
    if (inp->sigma_gauss[dir] > 0.0)
      envelope *= exp(-(pow(t,2))/(2.*pow(inp->sigma_gauss[dir],2)));
  }
  fout[0] = (envelope + inp->floor);
}

void func_one(double t, const double* xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 1.0;
}

void 
gk_species_projection_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s, 
  struct gkyl_gyrokinetic_projection inp, struct gk_proj *proj)
{
  proj->proj_id = inp.proj_id;
  // Context for c2p function passed to proj_on_basis.
  proj->proj_on_basis_c2p_ctx.cdim = app->cdim;
  proj->proj_on_basis_c2p_ctx.vdim = s->local_vel.ndim;
  proj->proj_on_basis_c2p_ctx.vel_map = s->vel_map;
  proj->proj_on_basis_c2p_ctx.pos_map = app->position_map;
  if (proj->proj_id == GKYL_PROJ_FUNC) {
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
  }
  else if (proj->proj_id == GKYL_PROJ_MAXWELLIAN_PRIM || proj->proj_id == GKYL_PROJ_BIMAXWELLIAN) {
    proj->dens = mkarr(false, app->basis.num_basis, app->local_ext.volume);
    proj->upar = mkarr(false, app->basis.num_basis, app->local_ext.volume);
    proj->proj_dens = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
        .grid = &app->grid,
        .basis = &app->basis,
        .qtype = GKYL_GAUSS_QUAD,
        .num_quad = app->basis.poly_order+1,
        .num_ret_vals = 1,
        .eval = inp.density,
        .ctx = inp.ctx_density,
        .c2p_func = proj_on_basis_c2p_position_func,
        .c2p_func_ctx = &proj->proj_on_basis_c2p_ctx,
      }
    );
    proj->proj_upar = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
        .grid = &app->grid,
        .basis = &app->basis,
        .qtype = GKYL_GAUSS_QUAD,
        .num_quad = app->basis.poly_order+1,
        .num_ret_vals = 1,
        .eval = inp.upar,
        .ctx = inp.ctx_upar,
        .c2p_func = proj_on_basis_c2p_position_func,
        .c2p_func_ctx = &proj->proj_on_basis_c2p_ctx,
      }
    );

    bool bimaxwellian = false;
    if (proj->proj_id == GKYL_PROJ_MAXWELLIAN_PRIM) {
      proj->vtsq = mkarr(false, app->basis.num_basis, app->local_ext.volume);
      proj->proj_temp = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
          .grid = &app->grid,
          .basis = &app->basis,
          .qtype = GKYL_GAUSS_QUAD,
          .num_quad = app->basis.poly_order+1,
          .num_ret_vals = 1,
          .eval = inp.temp,
          .ctx = inp.ctx_temp,
          .c2p_func = proj_on_basis_c2p_position_func,
          .c2p_func_ctx = &proj->proj_on_basis_c2p_ctx,
        }
      );

      proj->prim_moms_host = mkarr(false, 3*app->basis.num_basis, app->local_ext.volume);
      proj->prim_moms = mkarr(app->use_gpu, 3*app->basis.num_basis, app->local_ext.volume);
    }
    else {
      bimaxwellian = true;
      proj->vtsqpar = mkarr(false, app->basis.num_basis, app->local_ext.volume);
      proj->vtsqperp = mkarr(false, app->basis.num_basis, app->local_ext.volume);
      proj->proj_temppar = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
          .grid = &app->grid,
          .basis = &app->basis,
          .qtype = GKYL_GAUSS_QUAD,
          .num_quad = app->basis.poly_order+1,
          .num_ret_vals = 1,
          .eval = inp.temppar,
          .ctx = inp.ctx_temppar,
          .c2p_func = proj_on_basis_c2p_position_func,
          .c2p_func_ctx = &proj->proj_on_basis_c2p_ctx,
        }
      );
      proj->proj_tempperp = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
          .grid = &app->grid,
          .basis = &app->basis,
          .qtype = GKYL_GAUSS_QUAD,
          .num_quad = app->basis.poly_order+1,
          .num_ret_vals = 1,
          .eval = inp.tempperp,
          .ctx = inp.ctx_tempperp,
          .c2p_func = proj_on_basis_c2p_position_func,
          .c2p_func_ctx = &proj->proj_on_basis_c2p_ctx,
        }
      );

      proj->prim_moms_host = mkarr(false, 4*app->basis.num_basis, app->local_ext.volume);
      proj->prim_moms = mkarr(app->use_gpu, 4*app->basis.num_basis, app->local_ext.volume);      
    }

    // Maxwellian (or bi-Maxwellian) projection updater.
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
      .divide_jacobgeo = false, // final Jacobian multiplication will be handled in advance
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
  } else if (proj->proj_id == GKYL_PROJ_MAXWELLIAN_GAUSSIAN) {

    // Fill the box_size attribute of the projection (used for periodicity).
    for (int dir = 0; dir < app->cdim; ++dir)
      inp.box_size[dir] = app->grid.upper[dir] - app->grid.lower[dir];

    // First project the shape function, s(x), onto the DG basis and send it to device.
    proj->proj_shape = gkyl_proj_on_basis_new(&app->grid, &app->basis, app->poly_order + 1, 1, func_gaussian, &inp);
      
    struct gkyl_array *shape_ho = mkarr(false, app->basis.num_basis, app->local_ext.volume);
    gkyl_proj_on_basis_advance(proj->proj_shape, 0, &app->local, shape_ho);

    proj->shape_conf = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    gkyl_array_copy(proj->shape_conf, shape_ho);

    gkyl_array_release(shape_ho);

    // Build the integrant Jxyz * s(x), to normalize the shape function, and integrate it.
    struct gkyl_array *integrant = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    gkyl_dg_mul_op_range(app->basis, 0, integrant, 0, app->gk_geom->jacobgeo, 0, proj->shape_conf, &app->local);

    struct gkyl_array_integrate *int_op = gkyl_array_integrate_new(&app->grid, &app->basis, 1, 
      GKYL_ARRAY_INTEGRATE_OP_NONE, app->use_gpu);
    double *integral = app->use_gpu? gkyl_cu_malloc(sizeof(double)) : gkyl_malloc(sizeof(double));
    gkyl_array_integrate_advance(int_op, integrant, 1.0, NULL, &app->local, NULL, integral);

    double *red_integral = app->use_gpu? gkyl_cu_malloc(sizeof(double)) : gkyl_malloc(sizeof(double));
    gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, integral, red_integral);
    
    double red_integral_ho[1];
    if (app->use_gpu)
      gkyl_cu_memcpy(red_integral_ho, red_integral, sizeof(double), GKYL_CU_MEMCPY_D2H);
    else
      memcpy(red_integral_ho, red_integral, sizeof(double));

    // Scale the shape configuration function
    gkyl_array_scale(proj->shape_conf, 1.0/red_integral_ho[0]);
    
    // We can now build the moments of the projection.
    proj->dens = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    proj->upar = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    proj->vtsq = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    proj->prim_moms = mkarr(app->use_gpu, 4*app->basis.num_basis, app->local_ext.volume);      

    // n(x) =  #particle * shape(x)
    gkyl_array_set(proj->dens, inp.particle, proj->shape_conf);

    //  u(x) = 0
    gkyl_array_set(proj->upar, 0.0, proj->shape_conf);
    
    // T(x) = const
    // Compute the new temperature (no meaning if no particle).
    assert(inp.temp_max > 0);
    double temp = inp.particle == 0 ? inp.temp_max/2.0 : 2./3. * inp.energy/inp.particle;
    temp = temp > inp.temp_max ? inp.temp_max : temp; // saturate to max temperature.

    // Define a constant function = 1
    proj->proj_one = gkyl_proj_on_basis_new(&app->grid, &app->basis, app->poly_order + 1, 1, func_one, &inp);
    struct gkyl_array *one_ho = mkarr(false, app->basis.num_basis, app->local_ext.volume);
    gkyl_proj_on_basis_advance(proj->proj_one, 0, &app->local, one_ho);
    proj->one_conf = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    gkyl_array_copy(proj->one_conf, one_ho);
    gkyl_array_release(one_ho);

    gkyl_array_set(proj->vtsq, temp/s->info.mass, proj->one_conf);

    // release 
    if (app->use_gpu) {
      gkyl_cu_free(integral);
      gkyl_cu_free(red_integral);
    } else {
      gkyl_free(integral);
      gkyl_free(red_integral);
    }
    gkyl_array_release(integrant);
    gkyl_array_integrate_release(int_op);
  }
}

void
gk_species_projection_calc(gkyl_gyrokinetic_app *app, struct gk_species *s, 
  struct gk_proj *proj, struct gkyl_array *f, double tm)
{
  if (proj->proj_id == GKYL_PROJ_FUNC) {
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
  else if (proj->proj_id == GKYL_PROJ_MAXWELLIAN_PRIM) { 
    gkyl_proj_on_basis_advance(proj->proj_dens, tm, &app->local, proj->dens); 
    gkyl_proj_on_basis_advance(proj->proj_upar, tm, &app->local, proj->upar);
    gkyl_proj_on_basis_advance(proj->proj_temp, tm, &app->local, proj->vtsq);
    gkyl_array_scale(proj->vtsq, 1.0/s->info.mass);

    // proj_maxwellian expects the primitive moments as a single array.
    gkyl_array_set_offset(proj->prim_moms_host, 1.0, proj->dens, 0*app->basis.num_basis);
    gkyl_array_set_offset(proj->prim_moms_host, 1.0, proj->upar, 1*app->basis.num_basis);
    gkyl_array_set_offset(proj->prim_moms_host, 1.0, proj->vtsq, 2*app->basis.num_basis);

    // Copy the contents into the array we will use (potentially on GPUs).
    gkyl_array_copy(proj->prim_moms, proj->prim_moms_host);
    gkyl_gk_maxwellian_proj_on_basis_advance(proj->proj_max,
      &s->local, &app->local, proj->prim_moms, false, f);
  }
  else if (proj->proj_id == GKYL_PROJ_BIMAXWELLIAN) {
    gkyl_proj_on_basis_advance(proj->proj_dens, tm, &app->local, proj->dens); 
    gkyl_proj_on_basis_advance(proj->proj_upar, tm, &app->local, proj->upar);
    gkyl_proj_on_basis_advance(proj->proj_temppar, tm, &app->local, proj->vtsqpar);
    gkyl_proj_on_basis_advance(proj->proj_tempperp, tm, &app->local, proj->vtsqperp);
    gkyl_array_scale(proj->vtsqpar, 1.0/s->info.mass);
    gkyl_array_scale(proj->vtsqperp, 1.0/s->info.mass);

    // proj_bimaxwellian expects the primitive moments as a single array.
    gkyl_array_set_offset(proj->prim_moms_host, 1.0, proj->dens, 0*app->basis.num_basis);
    gkyl_array_set_offset(proj->prim_moms_host, 1.0, proj->upar, 1*app->basis.num_basis);
    gkyl_array_set_offset(proj->prim_moms_host, 1.0, proj->vtsqpar , 2*app->basis.num_basis);
    gkyl_array_set_offset(proj->prim_moms_host, 1.0, proj->vtsqperp , 3*app->basis.num_basis);  

    // Copy the contents into the array we will use (potentially on GPUs).
    gkyl_array_copy(proj->prim_moms, proj->prim_moms_host);
    gkyl_gk_maxwellian_proj_on_basis_advance(proj->proj_max,
      &s->local, &app->local, proj->prim_moms, false, f);
  } 
  else if (proj->proj_id == GKYL_PROJ_MAXWELLIAN_GAUSSIAN) {
    // LTE projection expects the primitive moments as a single array.
    gkyl_array_set_offset(proj->prim_moms, 1.0, proj->dens, 0*app->basis.num_basis);
    gkyl_array_set_offset(proj->prim_moms, 1.0, proj->upar, 1*app->basis.num_basis);
    gkyl_array_set_offset(proj->prim_moms, 1.0, proj->vtsq, 2*app->basis.num_basis);

    // Project the LTE distribution function to obtain f_lte.
    // Projection routine also corrects the density of the projected distribution function.
    gk_species_lte_from_moms(app, s, &s->lte, proj->prim_moms);
    gkyl_array_copy(f, s->lte.f_lte);
  }
  
  if (proj->proj_id == GKYL_PROJ_MAXWELLIAN_PRIM || proj->proj_id == GKYL_PROJ_BIMAXWELLIAN) {
    // Correct all the moments of the projected Maxwellian (or bi-Maxwellian) distribution function.
    if (proj->correct_all_moms) {
      struct gkyl_gk_maxwellian_correct_status status_corr;
      status_corr = gkyl_gk_maxwellian_correct_all_moments(proj->corr_max, 
        f, proj->prim_moms, &s->local, &app->local);
    } 
  }
  // Multiply by the configuration space jacobian.
  gkyl_dg_mul_conf_phase_op_range(&app->basis, &s->basis, f, 
    app->gk_geom->jacobgeo, f, &app->local, &s->local);
}

void
gk_species_projection_release(const struct gkyl_gyrokinetic_app *app, const struct gk_proj *proj)
{
  if (proj->proj_id == GKYL_PROJ_FUNC) {
    gkyl_proj_on_basis_release(proj->proj_func);
    if (app->use_gpu) {
      gkyl_array_release(proj->proj_host);
    }
  }
  else if (proj->proj_id == GKYL_PROJ_MAXWELLIAN_PRIM || proj->proj_id == GKYL_PROJ_BIMAXWELLIAN) { 
    gkyl_array_release(proj->dens);
    gkyl_array_release(proj->upar); 
    gkyl_array_release(proj->prim_moms_host);
    gkyl_array_release(proj->prim_moms);

    gkyl_proj_on_basis_release(proj->proj_dens);
    gkyl_proj_on_basis_release(proj->proj_upar);
    if (proj->proj_id == GKYL_PROJ_MAXWELLIAN_PRIM) {
      gkyl_array_release(proj->vtsq);
      gkyl_proj_on_basis_release(proj->proj_temp);
    }
    else if (proj->proj_id == GKYL_PROJ_BIMAXWELLIAN) {
      gkyl_array_release(proj->vtsqpar);
      gkyl_array_release(proj->vtsqperp);
      gkyl_proj_on_basis_release(proj->proj_temppar);
      gkyl_proj_on_basis_release(proj->proj_tempperp);
    }
    if (proj->proj_id == GKYL_PROJ_MAXWELLIAN_GAUSSIAN) {
      gkyl_array_release(proj->dens);
      gkyl_array_release(proj->upar); 
      gkyl_array_release(proj->vtsq);
      gkyl_array_release(proj->prim_moms);
      gkyl_array_release(proj->shape_conf);
      gkyl_array_release(proj->one_conf);

      gkyl_proj_on_basis_release(proj->proj_shape);
      gkyl_proj_on_basis_release(proj->proj_one);
    } else {
      gkyl_gk_maxwellian_proj_on_basis_release(proj->proj_max);
    }
    if (proj->correct_all_moms) {
      gkyl_gk_maxwellian_correct_release(proj->corr_max);    
    }
  }
}
