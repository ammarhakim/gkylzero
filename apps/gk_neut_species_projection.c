#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

void 
gk_neut_species_projection_init(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s, 
  struct gkyl_gyrokinetic_projection inp, struct gk_proj *proj)
{
  proj->proj_id = inp.proj_id;
  if (proj->proj_id == GKYL_PROJ_FUNC) {
    proj->proj_func = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
        .grid = &s->grid,
        .basis = &app->neut_basis,
        .qtype = GKYL_GAUSS_QUAD,
        .num_quad = app->neut_basis.poly_order+1,
        .num_ret_vals = 1,
        .eval = inp.func,
        .ctx = inp.ctx_func,
      }
    );
    if (app->use_gpu) {
      proj->proj_host = mkarr(false, app->neut_basis.num_basis, s->local_ext.volume);
    }
  }
  else if (proj->proj_id == GKYL_PROJ_MAXWELLIAN_PRIM) {
    int vdim = app->vdim+1; // neutral species are 3v otherwise
    proj->dens = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
    proj->udrift = mkarr(false, vdim*app->confBasis.num_basis, app->local_ext.volume);
    proj->vtsq = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
    proj->prim_moms = mkarr(false, (1+vdim)*app->confBasis.num_basis, app->local_ext.volume);
    // for correcting the density
    proj->dens_mod = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    if (app->use_gpu) {
      proj->prim_moms_dev = mkarr(app->use_gpu, (1+vdim)*app->confBasis.num_basis, app->local_ext.volume);
      proj->dens_dev = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
      proj->mem = gkyl_dg_bin_op_mem_cu_dev_new(app->local.volume, app->confBasis.num_basis);
    }
    else {
      proj->mem = gkyl_dg_bin_op_mem_new(app->local.volume, app->confBasis.num_basis);
    }

    proj->proj_dens = gkyl_proj_on_basis_new(&app->grid, &app->confBasis,
      app->neut_basis.poly_order+1, 1, inp.density, inp.ctx_density);
    proj->proj_udrift = gkyl_proj_on_basis_new(&app->grid, &app->confBasis,
      app->neut_basis.poly_order+1, vdim, inp.udrift, inp.ctx_udrift);
    proj->proj_temp = gkyl_proj_on_basis_new(&app->grid, &app->confBasis,
      app->neut_basis.poly_order+1, 1, inp.temp, inp.ctx_temp);

    proj->proj_max_prim = gkyl_proj_maxwellian_on_basis_new(&s->grid,
      &app->confBasis, &app->neut_basis, app->neut_basis.poly_order+1, s->vel_map, app->use_gpu);
  }
}

void
gk_neut_species_projection_calc(gkyl_gyrokinetic_app *app, const struct gk_neut_species *s, 
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
  else if (proj->proj_id == GKYL_PROJ_MAXWELLIAN_PRIM) {
    int vdim = app->vdim+1;
    gkyl_proj_on_basis_advance(proj->proj_dens, tm, &app->local_ext, proj->dens); 
    gkyl_proj_on_basis_advance(proj->proj_udrift, tm, &app->local_ext, proj->udrift);
    gkyl_proj_on_basis_advance(proj->proj_temp, tm, &app->local_ext, proj->vtsq);
    gkyl_array_scale(proj->vtsq, 1/s->info.mass);

    // proj_maxwellian expects the primitive moments as a single array.
    gkyl_array_set_offset(proj->prim_moms, 1.0, proj->udrift, 0*app->confBasis.num_basis);
    gkyl_array_set_offset(proj->prim_moms, 1.0, proj->vtsq, vdim*app->confBasis.num_basis);

    if (app->use_gpu) {
      gkyl_array_copy(proj->prim_moms_dev, proj->prim_moms);
      gkyl_array_copy(proj->dens_dev, proj->dens);
      gkyl_proj_maxwellian_on_basis_prim_mom(proj->proj_max_prim, &s->local_ext, &app->local_ext, 
        proj->dens_dev, proj->prim_moms_dev, f);
    }
    else {
      gkyl_proj_maxwellian_on_basis_prim_mom(proj->proj_max_prim, &s->local_ext, &app->local_ext, 
        proj->dens, proj->prim_moms, f);
    }

    // Now compute and scale the density to the desired density function 
    // based on input density from Maxwellian projection.
    gk_neut_species_moment_calc(&s->m0, s->local_ext, app->local_ext, f); 
    if (app->use_gpu) {
      gkyl_dg_div_op_range(proj->mem, app->confBasis, 
        0, proj->dens_mod, 0, proj->dens_dev, 0, s->m0.marr, &app->local);
    }
    else {
      gkyl_dg_div_op_range(proj->mem, app->confBasis, 
        0, proj->dens_mod, 0, proj->dens, 0, s->m0.marr, &app->local);
    }

    gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->neut_basis, f, 
      proj->dens_mod, f, &app->local_ext, &s->local_ext);
  }

  // Multiply by the configuration space jacobian.
  gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->neut_basis, f, 
    app->gk_geom->jacobgeo, f, &app->local_ext, &s->local_ext);  
}

void
gk_neut_species_projection_release(const struct gkyl_gyrokinetic_app *app, const struct gk_proj *proj)
{
  if (proj->proj_id == GKYL_PROJ_FUNC) {
    gkyl_proj_on_basis_release(proj->proj_func);
    if (app->use_gpu) {
      gkyl_array_release(proj->proj_host);
    }
  }
  else if (proj->proj_id == GKYL_PROJ_MAXWELLIAN_PRIM) { 
    gkyl_array_release(proj->dens);
    gkyl_array_release(proj->udrift); 
    gkyl_array_release(proj->vtsq);
    gkyl_array_release(proj->prim_moms);
    gkyl_array_release(proj->dens_mod); 
    if (app->use_gpu) {
      gkyl_array_release(proj->dens_dev);
      gkyl_array_release(proj->prim_moms_dev);      
    }
    gkyl_proj_on_basis_release(proj->proj_dens);
    gkyl_proj_on_basis_release(proj->proj_udrift);
    gkyl_proj_on_basis_release(proj->proj_temp);
    gkyl_proj_maxwellian_on_basis_release(proj->proj_max_prim);
    gkyl_dg_bin_op_mem_release(proj->mem);
  } 
}
