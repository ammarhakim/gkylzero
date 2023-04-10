#include <gkyl_bc_see.h>
#include <gkyl_bc_see_priv.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <assert.h>
#include <math.h>

// CHANGE NEEDED - Instead of passing gain coefficients as a precalculated double array, ideally we should calculate them in initialization. If possible, I would like to merge this object with bc_emission, though attempts to do so so far have proven messy. Should the user import the SEY function, or pick it from a list and supply fitting parameters? 

void
gkyl_bc_see_calc_normalization(double *tot_flux, void *ctx)
{
  struct bc_see_ctx *mc = (struct bc_see_ctx*) ctx;
  
  double effective_gamma = mc->numerator/mc->denominator;
  mc->k = 6.0*effective_gamma*tot_flux[0]*mc->phi*mc->phi*mc->mass/fabs(mc->charge);
}

// CHANGE NEEDED - function assumes 1X1V
static void 
chung_everhart(double t, const double *xn, double *out, void *ctx)
{
  struct bc_see_ctx *mc = (struct bc_see_ctx*) ctx;
  double k = mc->k;
  double phi = mc->phi;
  double E = 0.5*mc->mass*xn[1]*xn[1]/fabs(mc->charge);
  out[0] = k*E/((E + phi)*(E + phi)*(E + phi)*(E + phi));
}

struct gkyl_bc_see*
gkyl_bc_see_new(struct gkyl_rect_grid *grid, int dir, enum gkyl_edge_loc edge, const struct gkyl_range *conf_range_ext, const struct gkyl_range *local_range_ext,
		const int *num_ghosts, const struct gkyl_basis *cbasis, const struct gkyl_basis *basis,
		int num_comp, int cdim, int vdim, const double *bc_param, double *gain, double *elastic, bool use_gpu)
{
  // Allocate space for new updater.
  struct gkyl_bc_see *up = gkyl_malloc(sizeof(struct gkyl_bc_see));

  int vdir = dir + cdim;

  // CHANGE NEEDED - number of positive and negative cell currently assumes a symmetric domain
  int num_pos[GKYL_MAX_DIM];
  int num_neg[GKYL_MAX_DIM];
  for (int d=0; d<vdim; ++d) {
    // full phase space grid
    num_pos[cdim+d] = grid->cells[cdim+d]/2;
    num_neg[cdim+d] = grid->cells[cdim+d]/2;
  }

  up->dir = dir;
  up->cdim = cdim;
  up->edge = edge;
  up->use_gpu = use_gpu;
  up->cbasis = cbasis;
  up->grid = grid;

  // Create the skin/ghost ranges.
  gkyl_skin_ghost_ranges(&up->skin_r, &up->ghost_r, dir, edge,
    local_range_ext, num_ghosts);
  gkyl_skin_ghost_ranges(&up->cskin_r, &up->cghost_r, dir, edge,
    conf_range_ext, num_ghosts);
  gkyl_directional_ranges(&up->pos_r, &up->neg_r, dir, vdir, edge,
    local_range_ext, num_ghosts, num_pos, num_neg);

  struct bc_see_ctx *ctx = gkyl_malloc(sizeof(*ctx));
  up->func = bc_see;
  ctx->basis = basis;
  ctx->dir = dir;
  ctx->cdim = cdim;
  ctx->ncomp = num_comp;
  ctx->gain = gain;
  ctx->elastic = elastic;
  ctx->mass = bc_param[7];
  ctx->charge = bc_param[8];
  ctx->phi = bc_param[0];
  ctx->edge = edge;

  up->boundary_proj = gkyl_proj_on_basis_new(grid, basis,
    basis->poly_order+1, 1, chung_everhart, ctx);

  up->mcalc  = gkyl_dg_updater_moment_new(grid, cbasis, 
    basis, NULL, NULL, GKYL_MODEL_DEFAULT, "M1i", false, 0.0, use_gpu);
  int num_mom = gkyl_dg_updater_moment_num_mom(up->mcalc);
  up->flux = gkyl_array_new(GKYL_DOUBLE, num_mom*cbasis->num_basis, conf_range_ext->volume);
  up->f_store = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, local_range_ext->volume);
  up->ctx = ctx;
  up->ctx_on_dev = up->ctx;

  return up;
}

/* Modeled after gkyl_array_copy_to_buffer_fn */
void
gkyl_bc_see_advance(const struct gkyl_bc_see *up, struct gkyl_array *buff_arr, struct gkyl_array *f_arr)
{
  if (up->edge == GKYL_LOWER_EDGE) {
    gkyl_dg_updater_moment_advance(up->mcalc, &up->neg_r, &up->cskin_r, NULL, NULL, NULL, NULL, NULL, NULL, f_arr, up->flux);
  } else {
    gkyl_dg_updater_moment_advance(up->mcalc, &up->pos_r, &up->cskin_r, NULL, NULL, NULL, NULL, NULL, NULL, f_arr, up->flux);
  }
  struct gkyl_range_iter iter;
  // CHANGE NEEDED - these should all be stored in ctx so they don't have to be passed to the function
  double numerator[1] = {0.0};
  double denominator[1] = {0.0};
  double tot_flux[1] = {0.0};
  double xc[GKYL_MAX_DIM];
  double z[up->cbasis->num_basis];
  if (up->edge == GKYL_LOWER_EDGE) {
    z[up->dir] = -1.0;
  } else {
    z[up->dir] = 1.0;
  }
  const double *in_flux = gkyl_array_cfetch(up->flux, up->cskin_r.upper[0]);
  tot_flux[0] = up->cbasis->eval_expand(z, in_flux);
  double k[1] = {0.0};
  
  gkyl_range_iter_init(&iter, &up->skin_r);

  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&up->skin_r, iter.idx);

    const double *inp = gkyl_array_cfetch(f_arr, loc);
    const double *in_flux = gkyl_array_cfetch(up->flux, up->cskin_r.upper[0]);
    double *out = flat_fetch(buff_arr->data, f_arr->esznc*count);
    
    gkyl_rect_grid_cell_center(up->grid, iter.idx, xc);
    
    up->func(out, inp, up->ctx, iter.idx, xc);
    count += 1;
  }
  gkyl_bc_see_calc_normalization(tot_flux, up->ctx);

  struct bc_see_ctx *mc = (struct bc_see_ctx*) up->ctx;

  mc->k = k[0];

  // CHANGE NEEDED - order should be the opposite, or it causes problems when there are no elastic particles
  gkyl_array_copy_from_buffer(f_arr, buff_arr->data, up->ghost_r);

  // CHANGE NEEDED - projection should only be done one, and simply scaled with each advance
  gkyl_proj_on_basis_advance(up->boundary_proj, 0.0, &up->ghost_r, up->f_store);
  gkyl_array_accumulate_range(f_arr, 1.0, up->f_store, up->ghost_r);

  // gkyl_proj_on_basis_advance(up->boundary_proj, 0.0, &up->ghost_r, f_arr);
}

void gkyl_bc_see_release(struct gkyl_bc_see *up)
{
  // Release updater memory.
  gkyl_free(up);
}
