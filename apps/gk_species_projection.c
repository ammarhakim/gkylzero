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
      }
    );
    if (app->use_gpu)
      proj->proj_host = mkarr(false, app->basis.num_basis, s->local_ext.volume);
  }
  else if (proj->proj_id == GKYL_PROJ_MAXWELLIAN) {
    proj->m0 = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
    proj->upar = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
    proj->vtsq = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
    proj->prim_moms = mkarr(false, 2*app->confBasis.num_basis, app->local_ext.volume);
    // for correcting the density
    proj->m0mod = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    if (app->use_gpu) {
      proj->prim_moms_dev = mkarr(app->use_gpu, 2*app->confBasis.num_basis, app->local_ext.volume);
      proj->m0_dev = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
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
      .vel_range = &s->local_vel,
      .vmap_basis = s->vmap_basis,
      .vmap = s->vmap,
      .use_gpu = app->use_gpu,
    };
    proj->proj_max = gkyl_proj_maxwellian_on_basis_inew( &proj_inp );
  }
  else if (proj->proj_id == GKYL_PROJ_BIMAXWELLIAN) {
    proj->m0 = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
    proj->upar = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
    proj->vtsqpar = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
    proj->vtsqperp = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
    proj->prim_moms = mkarr(false, 4*app->confBasis.num_basis, app->local_ext.volume);
    // for correcting the density
    proj->m0mod = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    if (app->use_gpu) {
      proj->prim_moms_dev = mkarr(app->use_gpu, 4*app->confBasis.num_basis, app->local_ext.volume);
      proj->m0_dev = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
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

    proj->proj_bimax = gkyl_proj_bimaxwellian_on_basis_new(&s->grid,
      &app->confBasis, &app->basis, app->basis.poly_order+1, app->use_gpu);
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
  }
  else if (proj->proj_id == GKYL_PROJ_MAXWELLIAN) { 
    gkyl_proj_on_basis_advance(proj->proj_dens, tm, &app->local_ext, proj->m0); 
    gkyl_proj_on_basis_advance(proj->proj_upar, tm, &app->local_ext, proj->upar);
    gkyl_proj_on_basis_advance(proj->proj_temp, tm, &app->local_ext, proj->vtsq);
    gkyl_array_scale(proj->vtsq, 1/s->info.mass);

    // proj_maxwellian expects the primitive moments as a single array.
    gkyl_array_set_offset(proj->prim_moms, 1.0, proj->upar, 0*app->confBasis.num_basis);
    gkyl_array_set_offset(proj->prim_moms, 1.0, proj->vtsq  , 1*app->confBasis.num_basis);

    if (app->use_gpu) {
      gkyl_array_copy(proj->prim_moms_dev, proj->prim_moms);
      gkyl_array_copy(proj->m0_dev, proj->m0);
      gkyl_proj_gkmaxwellian_on_basis_prim_mom(proj->proj_max, &s->local_ext, &app->local_ext, 
        proj->m0_dev, proj->prim_moms_dev,
        app->gk_geom->bmag, app->gk_geom->bmag, s->info.mass, f);
    }
    else {
      gkyl_proj_gkmaxwellian_on_basis_prim_mom(proj->proj_max, &s->local_ext, &app->local_ext, 
        proj->m0, proj->prim_moms,
        app->gk_geom->bmag, app->gk_geom->bmag, s->info.mass, f);
    }  
  }
  else if (proj->proj_id == GKYL_PROJ_BIMAXWELLIAN) {
    gkyl_proj_on_basis_advance(proj->proj_dens, tm, &app->local_ext, proj->m0); 
    gkyl_proj_on_basis_advance(proj->proj_upar, tm, &app->local_ext, proj->upar);
    gkyl_proj_on_basis_advance(proj->proj_temppar, tm, &app->local_ext, proj->vtsqpar);
    gkyl_proj_on_basis_advance(proj->proj_tempperp, tm, &app->local_ext, proj->vtsqperp);
    gkyl_array_scale(proj->vtsqpar, 1.0/s->info.mass);
    gkyl_array_scale(proj->vtsqperp, 1.0/s->info.mass);

    // proj_bimaxwellian expects the primitive moments as a single array.
    gkyl_array_set_offset(proj->prim_moms, 1.0, proj->m0, 0*app->confBasis.num_basis);
    gkyl_array_set_offset(proj->prim_moms, 1.0, proj->upar, 1*app->confBasis.num_basis);
    gkyl_array_set_offset(proj->prim_moms, 1.0, proj->vtsqpar , 2*app->confBasis.num_basis);
    gkyl_array_set_offset(proj->prim_moms, 1.0, proj->vtsqperp , 3*app->confBasis.num_basis);  

    if (app->use_gpu) {
      gkyl_array_copy(proj->prim_moms_dev, proj->prim_moms);
      gkyl_array_copy(proj->m0_dev, proj->m0);
      gkyl_proj_bimaxwellian_on_basis_gyrokinetic_prim_mom(proj->proj_bimax, &s->local_ext, &app->local_ext, 
        proj->prim_moms_dev, app->gk_geom->bmag, app->gk_geom->bmag, s->info.mass, f);
    }
    else {
      gkyl_proj_bimaxwellian_on_basis_gyrokinetic_prim_mom(proj->proj_bimax, &s->local_ext, &app->local_ext, 
        proj->prim_moms, app->gk_geom->bmag, app->gk_geom->bmag, s->info.mass, f);
    } 
  }

  // Multiply by the velocity space jacobian. This has to happen before
  // trying to correct the moments otherwise moments will be wrong.
  gkyl_array_scale_by_cell(f, s->jacobvel);

  // Now compute and scale the density to the desired density function 
  // based on input density from Maxwellian projection. Also multiplies the
  // final distribution function by the Jacobian since we evolve J*f
  if (proj->proj_id == GKYL_PROJ_MAXWELLIAN || proj->proj_id == GKYL_PROJ_BIMAXWELLIAN) {
    gk_species_moment_calc(&s->m0, s->local_ext, app->local_ext, f); 
    if (app->use_gpu)
      gkyl_dg_div_op_range(proj->mem, app->confBasis, 
        0, proj->m0mod, 0, proj->m0_dev, 0, s->m0.marr, &app->local);
    else
      gkyl_dg_div_op_range(proj->mem, app->confBasis, 
        0, proj->m0mod, 0, proj->m0, 0, s->m0.marr, &app->local);

    gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->basis, f, 
        proj->m0mod, f, &app->local_ext, &s->local_ext);
  }

  // Multiply by the configuration space jacobian.
  gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->basis, f, 
    app->gk_geom->jacobgeo, f, &app->local_ext, &s->local_ext);  
}

void
gk_species_projection_release(const struct gkyl_gyrokinetic_app *app, const struct gk_proj *proj)
{
  if (proj->proj_id == GKYL_PROJ_FUNC) {
    gkyl_proj_on_basis_release(proj->proj_func);
    if (app->use_gpu)
      gkyl_array_release(proj->proj_host);
  }
  else if (proj->proj_id == GKYL_PROJ_MAXWELLIAN || proj->proj_id == GKYL_PROJ_BIMAXWELLIAN) { 
    gkyl_array_release(proj->m0);
    gkyl_array_release(proj->upar); 
    gkyl_array_release(proj->prim_moms);
    gkyl_array_release(proj->m0mod); 
    if (app->use_gpu) {
      gkyl_array_release(proj->m0_dev);
      gkyl_array_release(proj->prim_moms_dev);      
    }
    gkyl_proj_on_basis_release(proj->proj_dens);
    gkyl_proj_on_basis_release(proj->proj_upar);
    gkyl_dg_bin_op_mem_release(proj->mem);
    if (proj->proj_id == GKYL_PROJ_MAXWELLIAN) {
      gkyl_array_release(proj->vtsq);
      gkyl_proj_on_basis_release(proj->proj_temp);
      gkyl_proj_maxwellian_on_basis_release(proj->proj_max);
    }
    else if (proj->proj_id == GKYL_PROJ_BIMAXWELLIAN) {
      gkyl_array_release(proj->vtsqpar);
      gkyl_array_release(proj->vtsqperp);
      gkyl_proj_on_basis_release(proj->proj_temppar);
      gkyl_proj_on_basis_release(proj->proj_tempperp);
      gkyl_proj_bimaxwellian_on_basis_release(proj->proj_bimax);
    }
  } 
}
