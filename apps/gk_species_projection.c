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
    proj->prim_moms = mkarr(false, 2*app->confBasis.num_basis, app->local_ext.volume);
    // for correcting the density
    proj->dens_mod = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    if (app->use_gpu) {
      proj->prim_moms_dev = mkarr(app->use_gpu, 2*app->confBasis.num_basis, app->local_ext.volume);
      proj->dens_dev = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
      proj->mem = gkyl_dg_bin_op_mem_cu_dev_new(app->local.volume, app->confBasis.num_basis);
    }
    else {
      proj->mem = gkyl_dg_bin_op_mem_new(app->local.volume, app->confBasis.num_basis);
    }

    proj->proj_dens = gkyl_proj_on_basis_new(&app->grid, &app->confBasis,
      app->basis.poly_order+1, 1, inp.density, inp.ctx_density);
    proj->proj_upar = gkyl_proj_on_basis_new(&app->grid, &app->confBasis,
      app->basis.poly_order+1, 1, inp.upar, inp.ctx_upar);
    proj->proj_temp = gkyl_proj_on_basis_new(&app->grid, &app->confBasis,
      app->basis.poly_order+1, 1, inp.temp, inp.ctx_temp);

    struct gkyl_proj_maxwellian_on_basis_inp proj_inp = {
      .grid = &s->grid,
      .conf_basis = &app->confBasis,
      .phase_basis = &app->basis,
      .num_quad = app->basis.poly_order+1,
      .vel_map = s->vel_map,
      .use_gpu = app->use_gpu,
    };
    proj->proj_max_prim = gkyl_proj_maxwellian_on_basis_inew( &proj_inp );
  }
  else if (proj->proj_id == GKYL_PROJ_MAXWELLIAN_LAB) {
    proj->lab_moms = mkarr(app->use_gpu, 3*app->confBasis.num_basis, app->local_ext.volume);
    proj->lab_moms_host = proj->lab_moms;
    if (app->use_gpu) {
      proj->lab_moms_host = mkarr(false, 3*app->confBasis.num_basis, app->local_ext.volume);
    }

    proj->proj_lab_moms = gkyl_proj_on_basis_new(&app->grid, &app->confBasis,
      app->basis.poly_order+1, 3, inp.lab_moms, inp.ctx_lab_moms);
    // Maxwellian correction updater
    struct gkyl_correct_maxwellian_gyrokinetic_inp inp = {
      .phase_grid = &s->grid,
      .conf_grid = &app->grid,
      .phase_basis = &app->basis,
      .conf_basis = &app->confBasis,

      .conf_local = &app->local,
      .conf_local_ext = &app->local_ext,
      .mass = s->info.mass, 
      .gk_geom = app->gk_geom,
      .max_iter = 50, 
      .eps_err = 1.0e-14, 
      .use_gpu = app->use_gpu
    };
    proj->corr_max_lab = gkyl_correct_maxwellian_gyrokinetic_new(&inp);  

    struct gkyl_proj_maxwellian_on_basis_inp proj_inp = {
      .grid = &s->grid,
      .conf_basis = &app->confBasis,
      .phase_basis = &app->basis,
      .num_quad = app->basis.poly_order+1,
      .vel_map = s->vel_map,
      .use_gpu = app->use_gpu,
    };
    proj->proj_max_lab = gkyl_proj_maxwellian_on_basis_inew( &proj_inp );
  }
  else if (proj->proj_id == GKYL_PROJ_BIMAXWELLIAN) {
    proj->dens = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
    proj->upar = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
    proj->vtsqpar = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
    proj->vtsqperp = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
    proj->prim_moms = mkarr(false, 4*app->confBasis.num_basis, app->local_ext.volume);
    // for correcting the density
    proj->dens_mod = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    if (app->use_gpu) {
      proj->prim_moms_dev = mkarr(app->use_gpu, 4*app->confBasis.num_basis, app->local_ext.volume);
      proj->dens_dev = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
      proj->mem = gkyl_dg_bin_op_mem_cu_dev_new(app->local.volume, app->confBasis.num_basis);
    }
    else {
      proj->mem = gkyl_dg_bin_op_mem_new(app->local.volume, app->confBasis.num_basis);
    }

    proj->proj_dens = gkyl_proj_on_basis_new(&app->grid, &app->confBasis,
      app->basis.poly_order+1, 1, inp.density, inp.ctx_density);
    proj->proj_upar = gkyl_proj_on_basis_new(&app->grid, &app->confBasis,
      app->basis.poly_order+1, 1, inp.upar, inp.ctx_upar);
    proj->proj_temppar = gkyl_proj_on_basis_new(&app->grid, &app->confBasis,
      app->basis.poly_order+1, 1, inp.temppar, inp.ctx_temppar);
    proj->proj_tempperp = gkyl_proj_on_basis_new(&app->grid, &app->confBasis,
      app->basis.poly_order+1, 1, inp.tempperp, inp.ctx_tempperp);

    struct gkyl_proj_bimaxwellian_on_basis_inp proj_inp = {
      .grid = &s->grid,
      .conf_basis = &app->confBasis,
      .phase_basis = &app->basis,
      .num_quad = app->basis.poly_order+1,
      .vel_map = s->vel_map,
      .use_gpu = app->use_gpu,
    };
    proj->proj_bimax = gkyl_proj_bimaxwellian_on_basis_inew( &proj_inp );
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

    // Multiply by the velocity space jacobian.
    gkyl_array_scale_by_cell(f, s->vel_map->jacobvel);
  }
  else if (proj->proj_id == GKYL_PROJ_MAXWELLIAN_PRIM) { 
    gkyl_proj_on_basis_advance(proj->proj_dens, tm, &app->local_ext, proj->dens); 
    gkyl_proj_on_basis_advance(proj->proj_upar, tm, &app->local_ext, proj->upar);
    gkyl_proj_on_basis_advance(proj->proj_temp, tm, &app->local_ext, proj->vtsq);
    gkyl_array_scale(proj->vtsq, 1/s->info.mass);

    // proj_maxwellian expects the primitive moments as a single array.
    gkyl_array_set_offset(proj->prim_moms, 1.0, proj->upar, 0*app->confBasis.num_basis);
    gkyl_array_set_offset(proj->prim_moms, 1.0, proj->vtsq, 1*app->confBasis.num_basis);

    if (app->use_gpu) {
      gkyl_array_copy(proj->prim_moms_dev, proj->prim_moms);
      gkyl_array_copy(proj->dens_dev, proj->dens);
      gkyl_proj_gkmaxwellian_on_basis_prim_mom(proj->proj_max_prim, &s->local_ext, &app->local_ext, 
        proj->dens_dev, proj->prim_moms_dev,
        app->gk_geom->bmag, app->gk_geom->bmag, s->info.mass, f);
    }
    else {
      gkyl_proj_gkmaxwellian_on_basis_prim_mom(proj->proj_max_prim, &s->local_ext, &app->local_ext, 
        proj->dens, proj->prim_moms, app->gk_geom->bmag, app->gk_geom->bmag, s->info.mass, f);
    }
  }
  else if (proj->proj_id == GKYL_PROJ_MAXWELLIAN_LAB) { 
    gkyl_proj_on_basis_advance(proj->proj_lab_moms, tm, &app->local, proj->lab_moms_host); 
    if (app->use_gpu) {
      gkyl_array_copy(proj->lab_moms, proj->lab_moms_host);
    }

    gkyl_proj_gkmaxwellian_on_basis_lab_mom(proj->proj_max_lab, &s->local, &app->local, proj->lab_moms,
      app->gk_geom->bmag, app->gk_geom->bmag, s->info.mass, f);

    // Multiply by the velocity space jacobian. This has to happen before
    // trying to correct the moments otherwise moments will be wrong.
    gkyl_array_scale_by_cell(f, s->vel_map->jacobvel);

    gkyl_correct_maxwellian_gyrokinetic_advance(proj->corr_max_lab, f, proj->lab_moms, &app->local, &s->local);
  }
  else if (proj->proj_id == GKYL_PROJ_BIMAXWELLIAN) {
    gkyl_proj_on_basis_advance(proj->proj_dens, tm, &app->local_ext, proj->dens); 
    gkyl_proj_on_basis_advance(proj->proj_upar, tm, &app->local_ext, proj->upar);
    gkyl_proj_on_basis_advance(proj->proj_temppar, tm, &app->local_ext, proj->vtsqpar);
    gkyl_proj_on_basis_advance(proj->proj_tempperp, tm, &app->local_ext, proj->vtsqperp);
    gkyl_array_scale(proj->vtsqpar, 1.0/s->info.mass);
    gkyl_array_scale(proj->vtsqperp, 1.0/s->info.mass);

    // proj_bimaxwellian expects the primitive moments as a single array.
    gkyl_array_set_offset(proj->prim_moms, 1.0, proj->dens, 0*app->confBasis.num_basis);
    gkyl_array_set_offset(proj->prim_moms, 1.0, proj->upar, 1*app->confBasis.num_basis);
    gkyl_array_set_offset(proj->prim_moms, 1.0, proj->vtsqpar , 2*app->confBasis.num_basis);
    gkyl_array_set_offset(proj->prim_moms, 1.0, proj->vtsqperp , 3*app->confBasis.num_basis);  

    if (app->use_gpu) {
      gkyl_array_copy(proj->prim_moms_dev, proj->prim_moms);
      gkyl_array_copy(proj->dens_dev, proj->dens);
      gkyl_proj_bimaxwellian_on_basis_gyrokinetic_prim_mom(proj->proj_bimax, &s->local_ext, &app->local_ext, 
        proj->prim_moms_dev, app->gk_geom->bmag, app->gk_geom->bmag, s->info.mass, f);
    }
    else {
      gkyl_proj_bimaxwellian_on_basis_gyrokinetic_prim_mom(proj->proj_bimax, &s->local_ext, &app->local_ext, 
        proj->prim_moms, app->gk_geom->bmag, app->gk_geom->bmag, s->info.mass, f);
    }
  }

  if (proj->proj_id == GKYL_PROJ_MAXWELLIAN_PRIM || proj->proj_id == GKYL_PROJ_BIMAXWELLIAN) {
    // Multiply by the velocity space jacobian. This has to happen before
    // trying to correct the moments otherwise moments will be wrong.
    gkyl_array_scale_by_cell(f, s->vel_map->jacobvel);

    // Now compute and scale the density to the desired density function 
    // based on input density from Maxwellian projection.
    gk_species_moment_calc(&s->m0, s->local_ext, app->local_ext, f); 
    if (app->use_gpu) {
      gkyl_dg_div_op_range(proj->mem, app->confBasis, 
        0, proj->dens_mod, 0, proj->dens_dev, 0, s->m0.marr, &app->local);
    }
    else {
      gkyl_dg_div_op_range(proj->mem, app->confBasis, 
        0, proj->dens_mod, 0, proj->dens, 0, s->m0.marr, &app->local);
    }

    gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->basis, f, 
      proj->dens_mod, f, &app->local_ext, &s->local_ext);
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
    gkyl_array_release(proj->prim_moms);
    gkyl_array_release(proj->dens_mod); 
    if (app->use_gpu) {
      gkyl_array_release(proj->dens_dev);
      gkyl_array_release(proj->prim_moms_dev);      
    }
    gkyl_proj_on_basis_release(proj->proj_dens);
    gkyl_proj_on_basis_release(proj->proj_upar);
    gkyl_dg_bin_op_mem_release(proj->mem);
    if (proj->proj_id == GKYL_PROJ_MAXWELLIAN_PRIM) {
      gkyl_array_release(proj->vtsq);
      gkyl_proj_on_basis_release(proj->proj_temp);
      gkyl_proj_maxwellian_on_basis_release(proj->proj_max_prim);
    }
    else if (proj->proj_id == GKYL_PROJ_BIMAXWELLIAN) {
      gkyl_array_release(proj->vtsqpar);
      gkyl_array_release(proj->vtsqperp);
      gkyl_proj_on_basis_release(proj->proj_temppar);
      gkyl_proj_on_basis_release(proj->proj_tempperp);
      gkyl_proj_bimaxwellian_on_basis_release(proj->proj_bimax);
    }
  }
  else if (proj->proj_id == GKYL_PROJ_MAXWELLIAN_LAB) {
    gkyl_array_release(proj->lab_moms);
    if (app->use_gpu) {
      gkyl_array_release(proj->lab_moms_host);
    }

    gkyl_proj_on_basis_release(proj->proj_lab_moms);
    gkyl_correct_maxwellian_gyrokinetic_release(proj->corr_max_lab);
    gkyl_proj_maxwellian_on_basis_release(proj->proj_max_lab);
  }
}
