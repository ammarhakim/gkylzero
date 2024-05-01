#include <gkyl_alloc.h>
#include <gkyl_app.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_bc_basic.h>
#include <gkyl_dynvec.h>
#include <gkyl_elem_type.h>
#include <gkyl_eqn_type.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_eval_on_nodes.h>

#include <assert.h>
#include <time.h>

struct mapc2p_vel_nomap_ctx {
  int vdim;
};

void mapc2p_vel_nomap(double t, const double *vc, double* GKYL_RESTRICT vp, void *ctx)
{
  struct mapc2p_vel_nomap_ctx *nomap_ctx = ctx;
  int vdim = nomap_ctx->vdim;
  for (int d=0; d<vdim; d++) vp[d] = vc[d];
}

void
gk_species_init_vmap(struct gkyl_gk *gk, struct gkyl_gyrokinetic_app *app, struct gk_species *s)
{
  // Project velocity mappings, compute their derivative
  // and the corresponding Jacobian.
  int vmap_poly_order = 1;
  int vmapSq_poly_order = 2;

  int cdim = app->cdim, vdim = app->vdim;

  // Create the velocity space bases.
  if (app->use_gpu) {
    s->vmap_basis = gkyl_cart_modal_serendip_cu_dev_new(1, vmap_poly_order);
  } else {
    s->vmap_basis = gkyl_cart_modal_serendip_new(1, vmap_poly_order);
  }
  struct gkyl_basis vmap_basis_ho, vmap_basis2d_ho, vmapSq_basis_ho;
  gkyl_cart_modal_serendip(&vmap_basis_ho, 1, vmap_poly_order);
  gkyl_cart_modal_serendip(&vmap_basis2d_ho, vdim, vmap_poly_order);
  gkyl_cart_modal_serendip(&vmapSq_basis_ho, 1, vmapSq_poly_order);

  s->vmap = mkarr(app->use_gpu, vdim*vmap_basis_ho.num_basis, s->local_ext_vel.volume);
  s->vmapSq = mkarr(app->use_gpu, vdim*vmapSq_basis_ho.num_basis, s->local_ext_vel.volume);
  s->vmap_prime = mkarr(app->use_gpu, vdim, s->local_ext_vel.volume);
  s->jacobvel = mkarr(app->use_gpu, 1, s->local_ext.volume);

  // Project the velocity mapping (onto a 2D basis for vdim=2).
  struct gkyl_array *vmap2d = mkarr(app->use_gpu, vdim*vmap_basis2d_ho.num_basis, s->vmap->size);
  struct gkyl_array *vmap2d_ho = mkarr(false, vmap2d->ncomp, vmap2d->size);
  struct gkyl_array *vmapc1 = mkarr(app->use_gpu, 1, s->vmap->size);

  gkyl_eval_on_nodes *evup;
  if (s->info.mapc2p.is_mapped) {
    evup = gkyl_eval_on_nodes_new(&s->grid_vel, &vmap_basis2d_ho,
      vdim, s->info.mapc2p.mapping, s->info.mapc2p.ctx); 
  } else {
    struct mapc2p_vel_nomap_ctx nomap_ctx = { .vdim = vdim };
    evup = gkyl_eval_on_nodes_new(&s->grid_vel, &vmap_basis2d_ho,
      vdim, mapc2p_vel_nomap, &nomap_ctx); 
  }

  gkyl_eval_on_nodes_advance(evup, 0., &s->local_vel, vmap2d_ho);
  gkyl_array_copy(vmap2d, vmap2d_ho);

  // Extract the 1D basis expansions from the 2D basis expansion.
  for (int d=0; d<vdim; d++) {
    int numb_1d = vmap_basis_ho.num_basis;
    int numb_2d = vmap_basis2d_ho.num_basis;
    gkyl_array_set_offset(vmapc1, 1./sqrt(vdim), vmap2d, d*numb_2d);
    gkyl_array_set_offset(s->vmap, 1., vmapc1, d*numb_1d);
    gkyl_array_set_offset(vmapc1, 1./sqrt(vdim), vmap2d, d*numb_2d+d+1);
    gkyl_array_set_offset(s->vmap, 1., vmapc1, d*numb_1d+1);
  }

  gkyl_eval_on_nodes_release(evup);
  gkyl_array_release(vmapc1);
  gkyl_array_release(vmap2d);
  gkyl_array_release(vmap2d_ho);

  struct gkyl_array *vmap1d = mkarr(app->use_gpu, vmap_basis_ho.num_basis, s->vmap->size);
  struct gkyl_array *vmap1d_p2 = mkarr(app->use_gpu, vmapSq_basis_ho.num_basis, s->vmap->size);
  struct gkyl_array *vmap_prime1d = mkarr(app->use_gpu, 1, s->vmap_prime->size);
  gkyl_array_clear(vmap1d_p2, 0.);
  for (int d=0; d<vdim; d++) {
    gkyl_array_set_offset(vmap1d, 1., s->vmap, d*vmap_basis_ho.num_basis);

    // Compute the square mapping via weak multiplication.
    gkyl_array_set_offset(vmap1d_p2, 1., vmap1d, 0);
    gkyl_dg_mul_op(vmapSq_basis_ho, d, s->vmapSq, 0, vmap1d_p2, 0, vmap1d_p2);

    // Compute the derivative of the mapping.
    gkyl_array_set_offset(vmap_prime1d, sqrt(6.)/s->grid_vel.dx[d], vmap1d, 1);
    gkyl_array_set_offset(s->vmap_prime, 1., vmap_prime1d, d);
  }
  gkyl_array_release(vmap1d);
  gkyl_array_release(vmap_prime1d);
  gkyl_array_release(vmap1d_p2);

  // Compute the velocity space Jacobian (cell-wise constant).
  struct gkyl_array *jacobvel_ho = mkarr(false, s->jacobvel->ncomp, s->jacobvel->size);
  
  struct gkyl_array *vmap_prime_ho = mkarr(false, s->vmap_prime->ncomp, s->vmap_prime->size);
  gkyl_array_copy(vmap_prime_ho, s->vmap_prime);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &s->local);
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&s->local, iter.idx);
    double *jacv_d = gkyl_array_fetch(jacobvel_ho, linidx);
    jacv_d[0] = 1.;

    int vidx[vdim];
    for (int d=0; d<vdim; d++) vidx[d] = iter.idx[cdim+d];
    long vlinidx = gkyl_range_idx(&s->local_vel, vidx);
    double *vprime_d = gkyl_array_fetch(vmap_prime_ho, vlinidx);
    for (int d=0; d<vdim; d++)
      jacv_d[0] *= fabs(vprime_d[d]);
  }
  gkyl_array_copy(s->jacobvel, jacobvel_ho);

  struct gkyl_array *vmap_ho = mkarr(false, s->vmap->ncomp, s->vmap->size);
  gkyl_array_copy(vmap_ho, s->vmap);

  int rank, sz;
  int err = gkyl_comm_get_rank(s->comm, &rank);
  if (rank == 0) {
    // Write out the velocity mapping.
    const char *fmt0 = "%s-%s_mapc2p_vel.gkyl";
    sz = gkyl_calc_strlen(fmt0, app->name, s->info.name);
    char fileNm0[sz+1]; // ensures no buffer overflow
    snprintf(fileNm0, sizeof fileNm0, fmt0, app->name, s->info.name);
    gkyl_grid_sub_array_write(&s->grid_vel, &s->local_vel, NULL, vmap_ho, fileNm0);
  }

  // Write out the velocity space Jacobian.
  const char *fmt1 = "%s-%s_jacobvel.gkyl";
  sz = gkyl_calc_strlen(fmt1, app->name, s->info.name);
  char fileNm1[sz+1]; // ensures no buffer overflow
  snprintf(fileNm1, sizeof fileNm1, fmt1, app->name, s->info.name);
  gkyl_comm_array_write(s->comm, &s->grid, &s->local, NULL, jacobvel_ho, fileNm1);

  gkyl_array_release(vmap_ho);
  gkyl_array_release(vmap_prime_ho);
  gkyl_array_release(jacobvel_ho);

}

void
gk_species_release_vmap(const struct gkyl_gyrokinetic_app *app, const struct gk_species *s)
{
  // Free memory allocated for v-space mappings.
  gkyl_array_release(s->vmap);
  gkyl_array_release(s->vmap_prime);
  gkyl_array_release(s->vmapSq);
  if (app->use_gpu)
    gkyl_cart_modal_basis_release_cu(s->vmap_basis);
  else
    gkyl_cart_modal_basis_release(s->vmap_basis);
  gkyl_array_release(s->jacobvel);
}

// initialize species object
void
gk_species_init(struct gkyl_gk *gk_app_inp, struct gkyl_gyrokinetic_app *app, struct gk_species *gks)
{
  int cdim = app->cdim, vdim = app->vdim;
  int pdim = cdim+vdim;

  int cells[GKYL_MAX_DIM], ghost[GKYL_MAX_DIM];
  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];

  int cells_vel[GKYL_MAX_DIM], ghost_vel[GKYL_MAX_DIM];
  double lower_vel[GKYL_MAX_DIM], upper_vel[GKYL_MAX_DIM];

  for (int d=0; d<cdim; ++d) {
    cells[d] = gk_app_inp->cells[d];
    lower[d] = gk_app_inp->lower[d];
    upper[d] = gk_app_inp->upper[d];
    ghost[d] = 1;
  }
  for (int d=0; d<vdim; ++d) {
    // full phase space grid
    cells[cdim+d] = gks->info.cells[d];
    lower[cdim+d] = gks->info.lower[d];
    upper[cdim+d] = gks->info.upper[d];
    ghost[cdim+d] = 0; // no ghost-cells in velocity space

    // only velocity space
    cells_vel[d] = gks->info.cells[d];
    lower_vel[d] = gks->info.lower[d];
    upper_vel[d] = gks->info.upper[d];
    ghost_vel[d] = 0; // no ghost-cells in velocity space
  }
  // full phase space grid
  gkyl_rect_grid_init(&gks->grid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&gks->grid, ghost, &gks->global_ext, &gks->global);
  
  // velocity space grid
  gkyl_rect_grid_init(&gks->grid_vel, vdim, lower_vel, upper_vel, cells_vel);
  gkyl_create_grid_ranges(&gks->grid_vel, ghost_vel, &gks->local_ext_vel, &gks->local_vel);

  // phase-space communicator
  gks->comm = gkyl_comm_extend_comm(app->comm, &gks->local_vel);

  // create local and local_ext from app local range
  struct gkyl_range local;
  // local = conf-local X local_vel
  gkyl_range_ten_prod(&local, &app->local, &gks->local_vel);
  gkyl_create_ranges(&local, ghost, &gks->local_ext, &gks->local);

  // Velocity space mapping.
  gk_species_init_vmap(gk_app_inp, app, gks);

  // determine field-type 
  gks->gkfield_id = app->field->gkfield_id;
  if (gks->info.no_by) {
    gks->gkmodel_id = GKYL_GK_MODEL_NO_BY;
  }
  else {
    gks->gkmodel_id = GKYL_GK_MODEL_GEN_GEO;
  }

  // allocate distribution function arrays
  gks->f = mkarr(app->use_gpu, app->basis.num_basis, gks->local_ext.volume);
  gks->f1 = mkarr(app->use_gpu, app->basis.num_basis, gks->local_ext.volume);
  gks->fnew = mkarr(app->use_gpu, app->basis.num_basis, gks->local_ext.volume);

  gks->f_host = gks->f;
  if (app->use_gpu)
    gks->f_host = mkarr(false, app->basis.num_basis, gks->local_ext.volume);

  // allocate cflrate (scalar array)
  gks->cflrate = mkarr(app->use_gpu, 1, gks->local_ext.volume);

  if (app->use_gpu)
    gks->omega_cfl = gkyl_cu_malloc(sizeof(double));
  else 
    gks->omega_cfl = gkyl_malloc(sizeof(double));

  // Need to figure out size of alpha_surf and sgn_alpha_surf by finding size of surface basis set 
  struct gkyl_basis surf_basis, surf_quad_basis;
  if (app->poly_order > 1) {
    gkyl_cart_modal_serendip(&surf_basis, pdim-1, app->poly_order);
    gkyl_cart_modal_tensor(&surf_quad_basis, pdim-1, app->poly_order);
  }
  else {
    if (vdim>1) {
      gkyl_cart_modal_gkhybrid(&surf_basis, cdim-1, vdim); // p=2 in vparallel
      gkyl_cart_modal_gkhybrid(&surf_quad_basis, cdim-1, vdim); 
    }
    else {
      gkyl_cart_modal_serendip(&surf_basis, pdim-1, 2); // p=2 in vparallel
      gkyl_cart_modal_tensor(&surf_quad_basis, pdim-1, 2); 
    }
  }
  int alpha_surf_sz = (cdim+1)*surf_basis.num_basis;
  int sgn_alpha_surf_sz = (cdim+1)*surf_quad_basis.num_basis; // sign(alpha) is store at quadrature points

  // allocate arrays to store fields: 
  // 1. alpha_surf (surface phase space flux)
  // 2. sgn_alpha_surf (sign(alpha_surf) at quadrature points)
  // 3. const_sgn_alpha (boolean for if sign(alpha_surf) is a constant, either +1 or -1)
  gks->alpha_surf = mkarr(app->use_gpu, alpha_surf_sz, gks->local_ext.volume);
  gks->sgn_alpha_surf = mkarr(app->use_gpu, sgn_alpha_surf_sz, gks->local_ext.volume);
  gks->const_sgn_alpha = mk_int_arr(app->use_gpu, (cdim+1), gks->local_ext.volume);
  // 4. EM fields: phi and (if EM GK) Apar and d/dt Apar  
  gks->phi = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  gks->apar = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  gks->apardot = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);    

  gks->calc_gk_vars = gkyl_dg_calc_gyrokinetic_vars_new(&gks->grid, &app->confBasis, &app->basis, 
    gks->info.charge, gks->info.mass, gks->gkmodel_id, app->gk_geom, gks->vmap, gks->vmapSq, app->use_gpu);

  // by default, we do not have zero-flux boundary conditions in any direction
  bool is_zero_flux[2*GKYL_MAX_DIM] = {false};

  // Determine which directions are not periodic. If any BCs are zero-flux,
  // need to set it in is_zero_flux.
  // Keep a copy of num_periodic_dir and periodic_dirs in species so we can
  // modify it in GK_IWL BCs without modifying the app's.
  gks->num_periodic_dir = app->num_periodic_dir;
  for (int d=0; d<gks->num_periodic_dir; ++d)
    gks->periodic_dirs[d] = app->periodic_dirs[d];

  for (int d=0; d<app->cdim; ++d) gks->bc_is_np[d] = true;
  for (int d=0; d<gks->num_periodic_dir; ++d)
    gks->bc_is_np[gks->periodic_dirs[d]] = false;

  for (int dir=0; dir<app->cdim; ++dir) {
    gks->lower_bc[dir].type = gks->upper_bc[dir].type = GKYL_SPECIES_COPY;
    if (gks->bc_is_np[dir]) {
      const struct gkyl_gyrokinetic_bcs *bc;
      if (dir == 0)
        bc = &gks->info.bcx;
      else if (dir == 1)
        bc = &gks->info.bcy;
      else
        bc = &gks->info.bcz;

      gks->lower_bc[dir] = bc->lower;
      gks->upper_bc[dir] = bc->upper;
      if (gks->lower_bc[dir].type == GKYL_SPECIES_ZERO_FLUX) {
        is_zero_flux[dir] = true;
      }
      if (gks->upper_bc[dir].type == GKYL_SPECIES_ZERO_FLUX) {
        is_zero_flux[dir+pdim] = true;
      }
      if (gks->lower_bc[dir].type == GKYL_SPECIES_GK_IWL || gks->upper_bc[dir].type == GKYL_SPECIES_GK_IWL) {
        // Make the parallel direction periodic so that we sync the core before
        // applying sheath BCs in the SOL.
        gks->periodic_dirs[gks->num_periodic_dir] = app->cdim-1; // The last direction is the parallel one.
        gks->num_periodic_dir += 1;
        // Check that the LCFS is the same on both BCs and that it's on a cell boundary within our grid.
        double xLCFS = gks->lower_bc[dir].aux_parameter;
        assert(fabs(xLCFS-gks->upper_bc[dir].aux_parameter) < 1e-14);
        // Assume the split happens at a cell boundary and within the domain.
        assert((app->grid.lower[0]<xLCFS) && (xLCFS<app->grid.upper[0]));
        double needint = (xLCFS-app->grid.lower[0])/app->grid.dx[0];
        assert(floor(fabs(needint-floor(needint))) < 1.);
      }
    }
  }

  struct gkyl_dg_gyrokinetic_auxfields aux_inp = { .alpha_surf = gks->alpha_surf, 
    .sgn_alpha_surf = gks->sgn_alpha_surf, .const_sgn_alpha = gks->const_sgn_alpha, 
    .phi = gks->phi, .apar = gks->apar, .apardot = gks->apardot };
  // create solver
  gks->slvr = gkyl_dg_updater_gyrokinetic_new(&gks->grid, &app->confBasis, &app->basis, 
    &app->local, &gks->local_vel, &gks->local, is_zero_flux, 
    gks->info.charge, gks->info.mass, gks->gkmodel_id, 
    app->gk_geom, gks->vmap, gks->vmapSq, gks->vmap_prime, &aux_inp, app->use_gpu);

  // acquire equation object
  gks->eqn_gyrokinetic = gkyl_dg_updater_gyrokinetic_acquire_eqn(gks->slvr);

  // allocate data for density (for use in charge density accumulation and weak division for upar)
  gk_species_moment_init(app, gks, &gks->m0, "M0");
  // allocate data for integrated moments
  gk_species_moment_init(app, gks, &gks->integ_moms, "Integrated");

  // allocate data for diagnostic moments
  int ndm = gks->info.num_diag_moments;
  gks->moms = gkyl_malloc(sizeof(struct gk_species_moment[ndm]));
  for (int m=0; m<ndm; ++m)
    gk_species_moment_init(app, gks, &gks->moms[m], gks->info.diag_moments[m]);

  if (app->use_gpu) {
    gks->red_integ_diag = gkyl_cu_malloc(sizeof(double[vdim+2]));
    gks->red_integ_diag_global = gkyl_cu_malloc(sizeof(double[vdim+2]));
  } else {
    gks->red_integ_diag = gkyl_malloc(sizeof(double[vdim+2]));
    gks->red_integ_diag_global = gkyl_malloc(sizeof(double[vdim+2]));
  }
  // allocate dynamic-vector to store all-reduced integrated moments 
  gks->integ_diag = gkyl_dynvec_new(GKYL_DOUBLE, vdim+2);
  gks->is_first_integ_write_call = true;

  // initialize projection routine for initial conditions
  gk_species_projection_init(app, gks, gks->info.projection, &gks->proj_init);

  // set species source id
  gks->source_id = gks->info.source.source_id;
  
  // determine collision type to use in gyrokinetic update and initialize structs
  gks->collision_id = gks->info.collisions.collision_id;
  gks->lbo = (struct gk_lbo_collisions) { };
  gks->bgk = (struct gk_bgk_collisions) { };
  if (gks->collision_id == GKYL_LBO_COLLISIONS) 
    gk_species_lbo_init(app, gks, &gks->lbo);
  else if (gks->collision_id == GKYL_BGK_COLLISIONS)
    gk_species_bgk_init(app, gks, &gks->bgk);

  gks->has_reactions = false;
  gks->has_neutral_reactions = false;
  if (gks->info.react.num_react) {
    gks->has_reactions = true;
    gk_species_react_init(app, gks, gks->info.react, &gks->react, true);
  }
  if (gks->info.react_neut.num_react) {
    gks->has_neutral_reactions = true;
    gk_species_react_init(app, gks, gks->info.react_neut, &gks->react_neut, false);
  }

  // determine radiation type to use in gyrokinetic update
  gks->radiation_id = gks->info.radiation.radiation_id;

  // initialize boundary fluxes for diagnostics and, if present,
  // ambipolar potential solve
  gk_species_bflux_init(app, gks, &gks->bflux); 

  // initialize diffusion if present
  gks->has_diffusion = false;  
  if (gks->info.diffusion.num_diff_dir) {
    gks->has_diffusion = true;
    int diffusion_order = gks->info.diffusion.order ? gks->info.diffusion.order : 2;

    int szD = cdim*app->confBasis.num_basis;
    gks->diffD = mkarr(app->use_gpu, szD, app->local_ext.volume);
    bool diff_dir[GKYL_MAX_CDIM] = {false};

    int num_diff_dir = gks->info.diffusion.num_diff_dir ? gks->info.diffusion.num_diff_dir : app->cdim;
    // Assuming diffusion along x only for now.
    assert(num_diff_dir == 1);
    assert(gks->info.diffusion.diff_dirs[0] == 0);
    for (int d=0; d<num_diff_dir; ++d) {
      int dir = gks->info.diffusion.diff_dirs[d]; 
      diff_dir[dir] = 1; 
      gkyl_array_shiftc(gks->diffD, gks->info.diffusion.D[d]*pow(sqrt(2),app->cdim), dir);
    }
    // Multiply diffD by g^xx*jacobgeo.
    gkyl_dg_mul_op(app->confBasis, 0, gks->diffD, 0, app->gk_geom->gxxj, 0, gks->diffD);

    gks->diff_slvr = gkyl_dg_updater_diffusion_gyrokinetic_new(&gks->grid, &app->basis, &app->confBasis, 
      false, diff_dir, diffusion_order, &app->local, is_zero_flux, app->use_gpu);
  }

  // create ranges and allocate buffers for applying periodic and non-periodic BCs
  long buff_sz = 0;
  // compute buffer size needed
  for (int dir=0; dir<cdim; ++dir) {
    // Create local lower skin and ghost ranges for distribution function
    gkyl_skin_ghost_ranges(&gks->lower_skin[dir], &gks->lower_ghost[dir], dir, GKYL_LOWER_EDGE, &gks->local_ext, ghost);
    // Create local upper skin and ghost ranges for distribution function
    gkyl_skin_ghost_ranges(&gks->upper_skin[dir], &gks->upper_ghost[dir], dir, GKYL_UPPER_EDGE, &gks->local_ext, ghost);

    long vol = GKYL_MAX2(gks->lower_skin[dir].volume, gks->upper_skin[dir].volume);
    buff_sz = buff_sz > vol ? buff_sz : vol;
  }
  gks->bc_buffer = mkarr(app->use_gpu, app->basis.num_basis, buff_sz);
  // buffer arrays for fixed function boundary conditions on distribution function
  gks->bc_buffer_lo_fixed = mkarr(app->use_gpu, app->basis.num_basis, buff_sz);
  gks->bc_buffer_up_fixed = mkarr(app->use_gpu, app->basis.num_basis, buff_sz);

  for (int d=0; d<cdim; ++d) {
    // Copy BCs by default.
    enum gkyl_bc_basic_type bctype = GKYL_BC_COPY;

    // Lower BC.
    if (gks->lower_bc[d].type == GKYL_SPECIES_GK_SHEATH) {
      gks->bc_sheath_lo = gkyl_bc_sheath_gyrokinetic_new(d, GKYL_LOWER_EDGE, app->basis_on_dev.basis, 
        &gks->lower_skin[d], &gks->lower_ghost[d], &gks->local_vel, gks->vmap,
        cdim, 2.0*(gks->info.charge/gks->info.mass), app->use_gpu);
    }
    else if (gks->lower_bc[d].type == GKYL_SPECIES_GK_IWL) {
      double xLCFS = gks->lower_bc[d].aux_parameter;
      // Index of the cell that abuts the xLCFS from below.
      int idxLCFS_m = (xLCFS-1e-8 - app->grid.lower[0])/app->grid.dx[0]+1;
      gkyl_range_shorten_from_below(&gks->lower_skin_par_sol, &gks->lower_skin[d], 0, app->grid.cells[0]-idxLCFS_m+1);
      gkyl_range_shorten_from_below(&gks->lower_ghost_par_sol, &gks->lower_ghost[d], 0, app->grid.cells[0]-idxLCFS_m+1);

      gks->bc_sheath_lo = gkyl_bc_sheath_gyrokinetic_new(d, GKYL_LOWER_EDGE, app->basis_on_dev.basis, 
        &gks->lower_skin_par_sol, &gks->lower_ghost_par_sol, &gks->local_vel, gks->vmap,
        cdim, 2.0*(gks->info.charge/gks->info.mass), app->use_gpu);
    }
    else { 
      struct gkyl_gyrokinetic_projection proj_inp_bc_lo;

      if (gks->lower_bc[d].type == GKYL_SPECIES_COPY) 
        bctype = GKYL_BC_COPY;
      else if (gks->lower_bc[d].type == GKYL_SPECIES_ABSORB) 
        bctype = GKYL_BC_ABSORB;
      else if (gks->lower_bc[d].type == GKYL_SPECIES_REFLECT) 
        bctype = GKYL_BC_REFLECT;
      else if (gks->lower_bc[d].type == GKYL_SPECIES_INITIAL_SKIN) {
        bctype = GKYL_BC_FIXED_FUNC;
        proj_inp_bc_lo = gks->info.projection;
      }
      else if (gks->lower_bc[d].type == GKYL_SPECIES_FIXED_FUNC) {
        bctype = GKYL_BC_FIXED_FUNC;
        proj_inp_bc_lo = gks->lower_bc[d].projection;
      }

      gks->bc_lo[d] = gkyl_bc_basic_new(d, GKYL_LOWER_EDGE, bctype, app->basis_on_dev.basis,
        &gks->lower_skin[d], &gks->lower_ghost[d], gks->f->ncomp, app->cdim, app->use_gpu);

      if (gks->lower_bc[d].type == GKYL_SPECIES_INITIAL_SKIN || gks->lower_bc[d].type == GKYL_SPECIES_FIXED_FUNC) {
        // Fill the buffer used for BCs.
        struct gk_proj gk_proj_bc_lo;
        gk_species_projection_init(app, gks, proj_inp_bc_lo, &gk_proj_bc_lo);
        gk_species_projection_calc(app, gks, &gk_proj_bc_lo, gks->f1, 0.0); // Temporarily use f1.
        gkyl_bc_basic_buffer_fixed_func(gks->bc_lo[d], gks->bc_buffer_lo_fixed, gks->f1);
        gkyl_array_clear(gks->f1, 0.0);
        gk_species_projection_release(app, &gk_proj_bc_lo);
      }
    }

    // Upper BC.
    if (gks->upper_bc[d].type == GKYL_SPECIES_GK_SHEATH) {
      gks->bc_sheath_up = gkyl_bc_sheath_gyrokinetic_new(d, GKYL_UPPER_EDGE, app->basis_on_dev.basis, 
        &gks->upper_skin[d], &gks->upper_ghost[d], &gks->local_vel, gks->vmap,
        cdim, 2.0*(gks->info.charge/gks->info.mass), app->use_gpu);
    }
    else if (gks->lower_bc[d].type == GKYL_SPECIES_GK_IWL) {
      double xLCFS = gks->upper_bc[d].aux_parameter;
      // Index of the cell that abuts the xLCFS from below.
      int idxLCFS_m = (xLCFS-1e-8 - app->grid.lower[0])/app->grid.dx[0]+1;
      gkyl_range_shorten_from_below(&gks->upper_skin_par_sol, &gks->upper_skin[d], 0, app->grid.cells[0]-idxLCFS_m+1);
      gkyl_range_shorten_from_below(&gks->upper_ghost_par_sol, &gks->upper_ghost[d], 0, app->grid.cells[0]-idxLCFS_m+1);

      gks->bc_sheath_up = gkyl_bc_sheath_gyrokinetic_new(d, GKYL_UPPER_EDGE, app->basis_on_dev.basis, 
        &gks->upper_skin_par_sol, &gks->upper_ghost_par_sol, &gks->local_vel, gks->vmap,
        cdim, 2.0*(gks->info.charge/gks->info.mass), app->use_gpu);
    }
    else {
      struct gkyl_gyrokinetic_projection proj_inp_bc_up;

      if (gks->upper_bc[d].type == GKYL_SPECIES_COPY) 
        bctype = GKYL_BC_COPY;
      else if (gks->upper_bc[d].type == GKYL_SPECIES_ABSORB) 
        bctype = GKYL_BC_ABSORB;
      else if (gks->upper_bc[d].type == GKYL_SPECIES_REFLECT) 
        bctype = GKYL_BC_REFLECT;
      else if (gks->upper_bc[d].type == GKYL_SPECIES_INITIAL_SKIN) {
        bctype = GKYL_BC_FIXED_FUNC;
        proj_inp_bc_up = gks->info.projection;
      }
      else if (gks->upper_bc[d].type == GKYL_SPECIES_FIXED_FUNC) {
        bctype = GKYL_BC_FIXED_FUNC;
        proj_inp_bc_up = gks->upper_bc[d].projection;
      }

      gks->bc_up[d] = gkyl_bc_basic_new(d, GKYL_UPPER_EDGE, bctype, app->basis_on_dev.basis,
        &gks->upper_skin[d], &gks->upper_ghost[d], gks->f->ncomp, app->cdim, app->use_gpu);

      if (gks->upper_bc[d].type == GKYL_SPECIES_INITIAL_SKIN || gks->upper_bc[d].type == GKYL_SPECIES_FIXED_FUNC) {
        // Fill the buffer used for BCs.
        struct gk_proj gk_proj_bc_up;
        struct gkyl_gyrokinetic_projection proj_inp_bc_up;
        gk_species_projection_init(app, gks, proj_inp_bc_up, &gk_proj_bc_up);
        gk_species_projection_calc(app, gks, &gk_proj_bc_up, gks->f1, 0.0); // Temporarily use f1.
        gkyl_bc_basic_buffer_fixed_func(gks->bc_up[d], gks->bc_buffer_up_fixed, gks->f1);
        gkyl_array_clear(gks->f1, 0.0);
        gk_species_projection_release(app, &gk_proj_bc_up);
      }
    }
  }

}

void
gk_species_apply_ic(gkyl_gyrokinetic_app *app, struct gk_species *gks, double t0)
{
  gk_species_projection_calc(app, gks, &gks->proj_init, gks->f, t0);

  // We are pre-computing source for now as it is time-independent.
  if (gks->source_id)
    gk_species_source_calc(app, gks, &gks->src, t0);
}

// Compute the RHS for species update, returning maximum stable
// time-step.
double
gk_species_rhs(gkyl_gyrokinetic_app *app, struct gk_species *species,
  const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  gkyl_array_set(species->phi, 1.0, app->field->phi_smooth);

  gkyl_array_clear(species->cflrate, 0.0);
  gkyl_array_clear(rhs, 0.0);

  // Compute the surface expansion of the phase space flux
  // Note: Each cell stores the *lower* surface expansions of the 
  // phase space flux, so local_ext range needed to index the output
  // values of alpha_surf even though we only loop over local ranges
  // to avoid evaluating quantities such as geometry in ghost cells
  // where they are not defined.
  gkyl_dg_calc_gyrokinetic_vars_alpha_surf(species->calc_gk_vars, 
    &app->local, &species->local_vel, &species->local, &species->local_ext, 
    species->phi, species->alpha_surf, species->sgn_alpha_surf, species->const_sgn_alpha);

  gkyl_dg_updater_gyrokinetic_advance(species->slvr, &species->local, 
    fin, species->cflrate, rhs);

  if (species->collision_id == GKYL_LBO_COLLISIONS)
    gk_species_lbo_rhs(app, species, &species->lbo, fin, rhs);
  else if (species->collision_id == GKYL_BGK_COLLISIONS)
    gk_species_bgk_rhs(app, species, &species->bgk, fin, rhs);
  
  if (species->has_diffusion)
    gkyl_dg_updater_diffusion_gyrokinetic_advance(species->diff_slvr, &species->local, 
      species->diffD, app->gk_geom->jacobgeo_inv, fin, species->cflrate, rhs);

  if (species->radiation_id == GKYL_GK_RADIATION)
    gk_species_radiation_rhs(app, species, &species->rad, fin, rhs);

  if (species->has_reactions)
    gk_species_react_rhs(app, species, &species->react, fin, rhs);
  if (species->has_neutral_reactions)
    gk_species_react_rhs(app, species, &species->react_neut, fin, rhs);
  
  app->stat.nspecies_omega_cfl +=1;
  struct timespec tm = gkyl_wall_clock();
  gkyl_array_reduce_range(species->omega_cfl, species->cflrate, GKYL_MAX, &species->local);

  double omega_cfl_ho[1];
  if (app->use_gpu)
    gkyl_cu_memcpy(omega_cfl_ho, species->omega_cfl, sizeof(double), GKYL_CU_MEMCPY_D2H);
  else
    omega_cfl_ho[0] = species->omega_cfl[0];

  app->stat.species_omega_cfl_tm += gkyl_time_diff_now_sec(tm);
  
  return app->cfl/omega_cfl_ho[0];
}

// Apply boundary conditions to the distribution function.
void
gk_species_apply_bc(gkyl_gyrokinetic_app *app, const struct gk_species *species, struct gkyl_array *f)
{
  struct timespec wst = gkyl_wall_clock();
  
  int num_periodic_dir = species->num_periodic_dir, cdim = app->cdim;
  gkyl_comm_array_per_sync(species->comm, &species->local, &species->local_ext,
    num_periodic_dir, species->periodic_dirs, f); 
  
  for (int d=0; d<cdim; ++d) {
    if (species->bc_is_np[d]) {

      switch (species->lower_bc[d].type) {
        case GKYL_SPECIES_GK_SHEATH:
        case GKYL_SPECIES_GK_IWL:
          gkyl_bc_sheath_gyrokinetic_advance(species->bc_sheath_lo, app->field->phi_smooth, 
            app->field->phi_wall_lo, f, &app->local);
          break;
        case GKYL_SPECIES_COPY:
        case GKYL_SPECIES_REFLECT:
        case GKYL_SPECIES_ABSORB:
          gkyl_bc_basic_advance(species->bc_lo[d], species->bc_buffer, f);
          break;
        case GKYL_SPECIES_INITIAL_SKIN:
        case GKYL_SPECIES_FIXED_FUNC:
          gkyl_bc_basic_advance(species->bc_lo[d], species->bc_buffer_lo_fixed, f);
          break;
        case GKYL_SPECIES_ZERO_FLUX:
          break; // do nothing, BCs already applied in hyper_dg loop by not updating flux
          break;
        default:
          break;
      }

      switch (species->upper_bc[d].type) {
        case GKYL_SPECIES_GK_SHEATH:
        case GKYL_SPECIES_GK_IWL:
          gkyl_bc_sheath_gyrokinetic_advance(species->bc_sheath_up, app->field->phi_smooth, 
            app->field->phi_wall_up, f, &app->local);
          break;
        case GKYL_SPECIES_COPY:
        case GKYL_SPECIES_REFLECT:
        case GKYL_SPECIES_ABSORB:
          gkyl_bc_basic_advance(species->bc_up[d], species->bc_buffer, f);
          break;
        case GKYL_SPECIES_INITIAL_SKIN:
        case GKYL_SPECIES_FIXED_FUNC:
          gkyl_bc_basic_advance(species->bc_up[d], species->bc_buffer_up_fixed, f);
          break;
        case GKYL_SPECIES_ZERO_FLUX:
          break; // do nothing, BCs already applied in hyper_dg loop by not updating flux
        default:
          break;
      }      
    }
  }

  gkyl_comm_array_sync(species->comm, &species->local, &species->local_ext, f);

  app->stat.species_bc_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_species_coll_tm(gkyl_gyrokinetic_app *app)
{
  for (int i=0; i<app->num_species; ++i) {
    if (app->species[i].collision_id == GKYL_LBO_COLLISIONS) {
      struct gkyl_dg_updater_lbo_gyrokinetic_tm tm =
        gkyl_dg_updater_lbo_gyrokinetic_get_tm(app->species[i].lbo.coll_slvr);
      app->stat.species_lbo_coll_diff_tm[i] = tm.diff_tm;
      app->stat.species_lbo_coll_drag_tm[i] = tm.drag_tm;
    }
  }
}

void
gk_species_tm(gkyl_gyrokinetic_app *app)
{
  app->stat.species_rhs_tm = 0.0;
  for (int i=0; i<app->num_species; ++i) {
    struct gkyl_dg_updater_gyrokinetic_tm tm =
      gkyl_dg_updater_gyrokinetic_get_tm(app->species[i].slvr);
    app->stat.species_rhs_tm += tm.gyrokinetic_tm;
  }
}

// release resources for species
void
gk_species_release(const gkyl_gyrokinetic_app* app, const struct gk_species *s)
{
  // release various arrays and species objects
  gkyl_array_release(s->f);
  gkyl_array_release(s->f1);
  gkyl_array_release(s->fnew);
  gkyl_array_release(s->cflrate);
  gkyl_array_release(s->bc_buffer);
  gkyl_array_release(s->bc_buffer_lo_fixed);
  gkyl_array_release(s->bc_buffer_up_fixed);

  gk_species_projection_release(app, &s->proj_init);

  gkyl_comm_release(s->comm);

  if (app->use_gpu)
    gkyl_array_release(s->f_host);

  gkyl_array_release(s->alpha_surf);
  gkyl_array_release(s->sgn_alpha_surf);
  gkyl_array_release(s->const_sgn_alpha);
  gkyl_array_release(s->phi);
  gkyl_array_release(s->apar);
  gkyl_array_release(s->apardot);

  gkyl_dg_calc_gyrokinetic_vars_release(s->calc_gk_vars);

  // release equation object and solver
  gkyl_dg_eqn_release(s->eqn_gyrokinetic);
  gkyl_dg_updater_gyrokinetic_release(s->slvr);

  // release moment data
  gk_species_moment_release(app, &s->m0);
  for (int i=0; i<s->info.num_diag_moments; ++i)
    gk_species_moment_release(app, &s->moms[i]);
  gkyl_free(s->moms);
  gk_species_moment_release(app, &s->integ_moms); 
  gkyl_dynvec_release(s->integ_diag);

  if (s->source_id) {
    gk_species_source_release(app, &s->src);
  }

  if (s->collision_id == GKYL_LBO_COLLISIONS)
    gk_species_lbo_release(app, &s->lbo);
  else if (s->collision_id == GKYL_BGK_COLLISIONS)
    gk_species_bgk_release(app, &s->bgk);

  if (s->has_reactions)
    gk_species_react_release(app, &s->react);
  if (s->has_neutral_reactions)
    gk_species_react_release(app, &s->react_neut);  

  if (s->radiation_id == GKYL_GK_RADIATION) 
    gk_species_radiation_release(app, &s->rad);

  gk_species_bflux_release(app, &s->bflux);

  // Copy BCs are allocated by default. Need to free.
  for (int d=0; d<app->cdim; ++d) {
    if ((s->lower_bc[d].type == GKYL_SPECIES_GK_SHEATH) ||
        (s->lower_bc[d].type == GKYL_SPECIES_GK_IWL))
      gkyl_bc_sheath_gyrokinetic_release(s->bc_sheath_lo);
    else 
      gkyl_bc_basic_release(s->bc_lo[d]);
    
    if ((s->upper_bc[d].type == GKYL_SPECIES_GK_SHEATH) ||
        (s->upper_bc[d].type == GKYL_SPECIES_GK_IWL))
      gkyl_bc_sheath_gyrokinetic_release(s->bc_sheath_up);
    else 
      gkyl_bc_basic_release(s->bc_up[d]);
  }
  
  if (app->use_gpu) {
    gkyl_cu_free(s->omega_cfl);
    gkyl_cu_free(s->red_integ_diag);
    gkyl_cu_free(s->red_integ_diag_global);
  }
  else {
    gkyl_free(s->red_integ_diag);
    gkyl_free(s->red_integ_diag_global);
    gkyl_free(s->omega_cfl);
  }

  gk_species_release_vmap(app, s);
}
