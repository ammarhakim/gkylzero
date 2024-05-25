#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

void 
gk_species_projection_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s, 
  struct gkyl_gyrokinetic_projection inp, struct gk_proj *proj)
{
  proj->proj_id = inp.proj_id;
  if (proj->proj_id == GKYL_PROJ_FUNC) {
    proj->proj_func = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
        .grid = &s->grid,
        .basis = &app->basis,
        .qtype = GKYL_GAUSS_QUAD,
        .num_quad = app->basis.poly_order+1,
        .num_ret_vals = 1,
        .eval = inp.func,
        .ctx = inp.ctx_func,
        .vel_map = s->vel_map,
      }
    );
    if (app->use_gpu) {
      proj->proj_host = mkarr(false, app->basis.num_basis, s->local_ext.volume);
    }
  }
  else if (proj->proj_id == GKYL_PROJ_MAXWELLIAN_PRIM) {
    proj->dens = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
    proj->upar = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
    proj->vtsq = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
    proj->prim_moms_host = mkarr(false, 3*app->confBasis.num_basis, app->local_ext.volume);
    proj->prim_moms = mkarr(app->use_gpu, 3*app->confBasis.num_basis, app->local_ext.volume);

    proj->proj_dens = gkyl_proj_on_basis_new(&app->grid, &app->confBasis,
      app->basis.poly_order+1, 1, inp.density, inp.ctx_density);
    proj->proj_upar = gkyl_proj_on_basis_new(&app->grid, &app->confBasis,
      app->basis.poly_order+1, 1, inp.upar, inp.ctx_upar);
    proj->proj_temp = gkyl_proj_on_basis_new(&app->grid, &app->confBasis,
      app->basis.poly_order+1, 1, inp.temp, inp.ctx_temp);

    // Maxwellian correction updater
    struct gkyl_gyrokinetic_maxwellian_correct_inp inp_corr = {
      .phase_grid = &s->grid,
      .conf_basis = &app->confBasis,
      .phase_basis = &app->basis,
      .conf_range =  &app->local,
      .conf_range_ext = &app->local_ext,
      .gk_geom = app->gk_geom,
      .vel_map = s->vel_map,
      .mass = s->info.mass,
      .divide_jacobgeo = false, // final Jacobian multiplication will be handled in advance
      .use_last_converged = false, // do not use the unconverged moments if the scheme fails to converge
      .use_gpu = app->use_gpu,
      .max_iter = 50,
      .eps = 1e-10,
    };
    proj->corr_max = gkyl_gyrokinetic_maxwellian_correct_inew( &inp_corr );

    proj->correct_all_moms = false;
    if (inp.correct_all_moms) {
      proj->correct_all_moms = true;
    }

    // Maxwellian projection updater.
    proj->proj_max_prim = gkyl_proj_maxwellian_on_basis_new(&s->grid,
      &app->confBasis, &app->basis, app->basis.poly_order+1, s->vel_map, app->use_gpu);
  }
  else if (proj->proj_id == GKYL_PROJ_BIMAXWELLIAN) {
    proj->dens = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
    proj->upar = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
    proj->vtsqpar = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
    proj->vtsqperp = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
    proj->prim_moms_host = mkarr(false, 4*app->confBasis.num_basis, app->local_ext.volume);
    proj->prim_moms = mkarr(app->use_gpu, 4*app->confBasis.num_basis, app->local_ext.volume);

    proj->proj_dens = gkyl_proj_on_basis_new(&app->grid, &app->confBasis,
      app->basis.poly_order+1, 1, inp.density, inp.ctx_density);
    proj->proj_upar = gkyl_proj_on_basis_new(&app->grid, &app->confBasis,
      app->basis.poly_order+1, 1, inp.upar, inp.ctx_upar);
    proj->proj_temppar = gkyl_proj_on_basis_new(&app->grid, &app->confBasis,
      app->basis.poly_order+1, 1, inp.temppar, inp.ctx_temppar);
    proj->proj_tempperp = gkyl_proj_on_basis_new(&app->grid, &app->confBasis,
      app->basis.poly_order+1, 1, inp.tempperp, inp.ctx_tempperp);

    // Maxwellian correction updater
    struct gkyl_gyrokinetic_maxwellian_correct_inp inp_corr = {
      .phase_grid = &s->grid,
      .conf_basis = &app->confBasis,
      .phase_basis = &app->basis,
      .conf_range =  &app->local,
      .conf_range_ext = &app->local_ext,
      .gk_geom = app->gk_geom,
      .vel_map = s->vel_map,
      .divide_jacobgeo = false, // final Jacobian multiplication will be handled in advance
      .use_last_converged = false, // do not use the unconverged moments if the scheme fails to converge
      .correct_bimaxwellian = true, // correct the Bimaxwellian moments
      .mass = s->info.mass,
      .use_gpu = app->use_gpu,
      .max_iter = 50,
      .eps = 1e-10,
    };
    proj->corr_max = gkyl_gyrokinetic_maxwellian_correct_inew( &inp_corr );

    proj->correct_all_moms = false;
    if (inp.correct_all_moms) {
      proj->correct_all_moms = true;
    }

    // Maxwellian projection updater.
    proj->proj_bimax = gkyl_proj_bimaxwellian_on_basis_new(&s->grid,
      &app->confBasis, &app->basis, app->basis.poly_order+1, s->vel_map, app->use_gpu);
  }
}

void
gk_species_projection_calc(gkyl_gyrokinetic_app *app, const struct gk_species *s, 
  struct gk_proj *proj, struct gkyl_array *f, double tm)
{
  if (proj->proj_id == GKYL_PROJ_FUNC) {
    if (app->use_gpu) {
      gkyl_proj_on_basis_advance(proj->proj_func, tm, &s->local_ext, proj->proj_host);
      gkyl_array_copy(f, proj->proj_host);
    }
    else {
      gkyl_proj_on_basis_advance(proj->proj_func, tm, &s->local_ext, f);
    }

    // Multiply by the gyrocenter coord jacobian (bmag).
    gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->basis, f, 
        app->gk_geom->bmag, f, &app->local, &s->local);      
  }
  else if (proj->proj_id == GKYL_PROJ_MAXWELLIAN_PRIM) { 
    gkyl_proj_on_basis_advance(proj->proj_dens, tm, &app->local_ext, proj->dens); 
    gkyl_proj_on_basis_advance(proj->proj_upar, tm, &app->local_ext, proj->upar);
    gkyl_proj_on_basis_advance(proj->proj_temp, tm, &app->local_ext, proj->vtsq);
    gkyl_array_scale(proj->vtsq, 1.0/s->info.mass);

    // proj_maxwellian expects the primitive moments as a single array.
    gkyl_array_set_offset(proj->prim_moms_host, 1.0, proj->dens, 0*app->confBasis.num_basis);
    gkyl_array_set_offset(proj->prim_moms_host, 1.0, proj->upar, 1*app->confBasis.num_basis);
    gkyl_array_set_offset(proj->prim_moms_host, 1.0, proj->vtsq, 2*app->confBasis.num_basis);

    // Copy the contents into the array we will use (potentially on GPUs).
    gkyl_array_copy(proj->prim_moms, proj->prim_moms_host);
    gkyl_proj_gkmaxwellian_on_basis_prim_mom(proj->proj_max_prim, &s->local_ext, &app->local_ext, 
      proj->prim_moms, app->gk_geom->bmag, app->gk_geom->bmag, s->info.mass, f);
  }
  else if (proj->proj_id == GKYL_PROJ_BIMAXWELLIAN) {
    gkyl_proj_on_basis_advance(proj->proj_dens, tm, &app->local_ext, proj->dens); 
    gkyl_proj_on_basis_advance(proj->proj_upar, tm, &app->local_ext, proj->upar);
    gkyl_proj_on_basis_advance(proj->proj_temppar, tm, &app->local_ext, proj->vtsqpar);
    gkyl_proj_on_basis_advance(proj->proj_tempperp, tm, &app->local_ext, proj->vtsqperp);
    gkyl_array_scale(proj->vtsqpar, 1.0/s->info.mass);
    gkyl_array_scale(proj->vtsqperp, 1.0/s->info.mass);

    // proj_bimaxwellian expects the primitive moments as a single array.
    gkyl_array_set_offset(proj->prim_moms_host, 1.0, proj->dens, 0*app->confBasis.num_basis);
    gkyl_array_set_offset(proj->prim_moms_host, 1.0, proj->upar, 1*app->confBasis.num_basis);
    gkyl_array_set_offset(proj->prim_moms_host, 1.0, proj->vtsqpar , 2*app->confBasis.num_basis);
    gkyl_array_set_offset(proj->prim_moms_host, 1.0, proj->vtsqperp , 3*app->confBasis.num_basis);  

    // Copy the contents into the array we will use (potentially on GPUs).
    gkyl_array_copy(proj->prim_moms, proj->prim_moms_host);
    gkyl_proj_bimaxwellian_on_basis_gyrokinetic_prim_mom(proj->proj_bimax, &s->local_ext, &app->local_ext, 
      proj->prim_moms, app->gk_geom->bmag, app->gk_geom->bmag, s->info.mass, f);
  }

  // Multiply by the velocity space jacobian.
  gkyl_array_scale_by_cell(f, s->vel_map->jacobvel);

  if (proj->proj_id == GKYL_PROJ_MAXWELLIAN_PRIM || proj->proj_id == GKYL_PROJ_BIMAXWELLIAN) {
    // Now compute and scale the density to the desired density function 
    // based on input density from Maxwellian projection and potentially correct
    // u_par and T/m as well using gyrokinetic_maxwellian_correct updater.
    gkyl_gyrokinetic_maxwellian_correct_density_moment(proj->corr_max, 
      f, proj->prim_moms, &s->local, &app->local);
    if (proj->correct_all_moms) {
      struct gkyl_gyrokinetic_maxwellian_correct_status status_corr;
      if (proj->proj_id == GKYL_PROJ_MAXWELLIAN_PRIM) {
        status_corr = gkyl_gyrokinetic_maxwellian_correct_all_moments(proj->corr_max, 
          f, proj->prim_moms, &s->local, &app->local);    
      }  
      else if (proj->proj_id == GKYL_PROJ_BIMAXWELLIAN) {
        status_corr = gkyl_gyrokinetic_bimaxwellian_correct_all_moments(proj->corr_max, 
          f, proj->prim_moms, &s->local, &app->local);    
      }  
    }
  }
  // Multiply by the configuration space jacobian.
  gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->basis, f, 
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
      gkyl_proj_maxwellian_on_basis_release(proj->proj_max_prim);
      gkyl_gyrokinetic_maxwellian_correct_release(proj->corr_max);
    }
    else if (proj->proj_id == GKYL_PROJ_BIMAXWELLIAN) {
      gkyl_array_release(proj->vtsqpar);
      gkyl_array_release(proj->vtsqperp);
      gkyl_proj_on_basis_release(proj->proj_temppar);
      gkyl_proj_on_basis_release(proj->proj_tempperp);
      gkyl_proj_bimaxwellian_on_basis_release(proj->proj_bimax);
      gkyl_gyrokinetic_maxwellian_correct_release(proj->corr_max);
    }
  }
}
