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
  }
}

void
gk_species_projection_calc(gkyl_gyrokinetic_app *app, const struct gk_species *s, 
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
    gkyl_gk_maxwellian_proj_on_basis_release(proj->proj_max);
    if (proj->correct_all_moms) {
      gkyl_gk_maxwellian_correct_release(proj->corr_max);    
    }
  }
}
