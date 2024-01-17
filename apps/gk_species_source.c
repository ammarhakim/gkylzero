#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

void 
gk_species_source_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s, struct gk_source *src)
{
  // we need to ensure source has same shape as distribution function
  src->source = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume);
  src->source_id = s->source_id;
  src->source_host = src->source;
  if (app->use_gpu)
    src->source_host = mkarr(false, app->basis.num_basis, s->local_ext.volume);

  src->write_source = s->info.source.write_source; // optional flag to write out source

  if (src->source_id == GKYL_FUNC_SOURCE) {
    src->source_proj = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
        .grid = &s->grid,
        .basis = &app->basis,
        .qtype = GKYL_GAUSS_QUAD,
        .num_quad = app->basis.poly_order+1,
        .num_ret_vals = 1,
        .eval = s->info.source.profile,
        .ctx = s->info.source.ctx_profile
      }
    );
  }
}

void
gk_species_source_calc(gkyl_gyrokinetic_app *app, struct gk_species *species, double tm)
{
  if (species->source_id == GKYL_MAXWELLIAN_SOURCE) {
    int poly_order = app->poly_order;
    // Project n, upar, and vt^2 based on input functions
    struct gkyl_array *m0 = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
    struct gkyl_array *upar = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
    struct gkyl_array *vtsq = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
    gkyl_proj_on_basis *proj_m0 = gkyl_proj_on_basis_new(&app->grid, &app->confBasis,
      poly_order+1, 1, species->info.source.density_profile, species->info.source.ctx_density);
    gkyl_proj_on_basis *proj_upar = gkyl_proj_on_basis_new(&app->grid, &app->confBasis,
      poly_order+1, 1, species->info.source.upar_profile, species->info.source.ctx_upar);
    gkyl_proj_on_basis *proj_vtsq = gkyl_proj_on_basis_new(&app->grid, &app->confBasis,
      poly_order+1, 1, species->info.source.temp_profile, species->info.source.ctx_temp);

    gkyl_proj_on_basis_advance(proj_m0, 0.0, &app->local_ext, m0); 
    gkyl_proj_on_basis_advance(proj_upar, 0.0, &app->local_ext, upar);
    gkyl_proj_on_basis_advance(proj_vtsq, 0.0, &app->local_ext, vtsq);
    gkyl_array_scale(vtsq, 1/species->info.mass);

    // proj_maxwellian expects the primitive moments as a single array.
    struct gkyl_array *prim_moms = mkarr(false, 2*app->confBasis.num_basis, app->local_ext.volume);
    gkyl_array_set_offset(prim_moms, 1.0, upar, 0*app->confBasis.num_basis);
    gkyl_array_set_offset(prim_moms, 1.0, vtsq  , 1*app->confBasis.num_basis);

    // Initialize Maxwellian projection object
    gkyl_proj_maxwellian_on_basis *proj_max = gkyl_proj_maxwellian_on_basis_new(&species->grid,
        &app->confBasis, &app->basis, poly_order+1, app->use_gpu);

    // If on GPUs, need to copy n, upar, and vt^2 onto device
    struct gkyl_array *prim_moms_dev, *m0_dev;
    if (app->use_gpu) {
      prim_moms_dev = mkarr(app->use_gpu, 2*app->confBasis.num_basis, app->local_ext.volume);
      m0_dev = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);

      gkyl_array_copy(prim_moms_dev, prim_moms);
      gkyl_array_copy(m0_dev, m0);
      gkyl_proj_gkmaxwellian_on_basis_prim_mom(proj_max, &species->local_ext, &app->local_ext, m0_dev, prim_moms_dev,
          app->gk_geom->bmag, app->gk_geom->bmag, species->info.mass, species->src.source);
    }
    else {
      gkyl_proj_gkmaxwellian_on_basis_prim_mom(proj_max, &species->local_ext, &app->local_ext, m0, prim_moms,
          app->gk_geom->bmag, app->gk_geom->bmag, species->info.mass, species->src.source);
    }
    // Now compute and scale the density to the desired density function based on input density from Maxwellian projection
    gk_species_moment_calc(&species->m0, species->local_ext, app->local_ext, species->src.source); 

    // Rescale projected density to desired input density function
    struct gkyl_array *m0mod = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    struct gkyl_dg_bin_op_mem *mem;
    if (app->use_gpu) {
      mem = gkyl_dg_bin_op_mem_cu_dev_new(app->local.volume, app->confBasis.num_basis);
      gkyl_dg_div_op_range(mem, app->confBasis, 0, m0mod, 0, m0_dev, 0, species->m0.marr, &app->local);
    }
    else {
      mem = gkyl_dg_bin_op_mem_new(app->local.volume, app->confBasis.num_basis);
      gkyl_dg_div_op_range(mem, app->confBasis, 0, m0mod, 0, m0, 0, species->m0.marr, &app->local);
    }
    gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->basis, species->src.source, 
        m0mod, species->src.source, &app->local_ext, &species->local_ext);

    // multiply final distribution function by Jacobian
    gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->basis, species->src.source, 
        app->gk_geom->jacobgeo, species->src.source, &app->local_ext, &species->local_ext);

    // Free temporary variables and projection objects
    gkyl_array_release(m0);
    gkyl_array_release(upar); 
    gkyl_array_release(vtsq);
    gkyl_array_release(prim_moms);
    if (app->use_gpu) {
      gkyl_array_release(m0_dev);
      gkyl_array_release(prim_moms_dev);      
    }
    gkyl_proj_on_basis_release(proj_m0);
    gkyl_proj_on_basis_release(proj_upar);
    gkyl_proj_on_basis_release(proj_vtsq);
    gkyl_array_release(m0mod); 
    gkyl_dg_bin_op_mem_release(mem);
    gkyl_proj_maxwellian_on_basis_release(proj_max);
  }
  else {
    gkyl_proj_on_basis_advance(species->src.source_proj, tm, &species->local_ext, species->src.source_host);
    if (app->use_gpu) // note: source_host is same as source when not on GPUs
      gkyl_array_copy(species->src.source, species->src.source_host);
  }
}

// computes rhs of the source
void
gk_species_source_rhs(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_source *src, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  gkyl_array_accumulate(rhs, 1.0, src->source);
}

void
gk_species_source_release(const struct gkyl_gyrokinetic_app *app, const struct gk_source *src)
{
  gkyl_array_release(src->source);
  if (app->use_gpu) 
    gkyl_array_release(src->source_host);
  if (src->source_id == GKYL_FUNC_SOURCE)   
    gkyl_proj_on_basis_release(src->source_proj);
}
