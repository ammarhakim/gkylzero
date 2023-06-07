#include <gkyl_bc_emission_spectrum.h>
#include <gkyl_bc_emission_spectrum_priv.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <assert.h>
#include <math.h>

// CHANGE NEEDED - Instead of passing gain coefficients as a precalculated double array,
// ideally we should calculate them in initialization.

struct gkyl_bc_emission_spectrum*
gkyl_bc_emission_spectrum_new(int dir, enum gkyl_edge_loc edge,
  const struct gkyl_range *ghost_r,
  enum gkyl_bc_emission_spectrum_type bctype, const struct gkyl_basis *cbasis, const struct gkyl_basis *basis,
  int cdim, int vdim, bool use_gpu)
{
  // Allocate space for new updater.
  struct gkyl_bc_emission_spectrum *up = gkyl_malloc(sizeof(struct gkyl_bc_emission_spectrum));

  int vdir = dir + cdim;

  up->dir = dir;
  up->cdim = cdim;
  up->edge = edge;
  up->use_gpu = use_gpu;
  up->cbasis = cbasis;
  up->func = bc_weighted_gamma;
  up->ghost_r = *ghost_r;

  switch (bctype) {
    case GKYL_BC_CHUNG_EVERHART:
      up->norm = chung_everhart_norm;
      break;

    case GKYL_BC_GAUSSIAN:
      up->norm = gaussian_norm;
      break;
      
    default:
      assert(false);
      break;
  }

  struct bc_emission_spectrum_ctx *ctx = gkyl_malloc(sizeof(*ctx));
  ctx->basis = basis;
  ctx->dir = dir;
  ctx->cdim = cdim;
  ctx->vdim = vdim;
  ctx->edge = edge;
  
  up->ctx = ctx;
  up->ctx_on_dev = up->ctx;

  return up;
}

double
gkyl_bc_emission_spectrum_advance_cross(const struct gkyl_bc_emission_spectrum *up,
  struct gkyl_array *f_self, struct gkyl_array *f_other, struct gkyl_array *f_proj,
  double *bc_param, double flux, struct gkyl_rect_grid *grid, double *gain, const struct gkyl_range *other_r)
{
  double weight[2] = {0.0, 0.0};

  double xc[GKYL_MAX_DIM];

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, other_r);

  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(other_r, iter.idx);

    const double *inp = gkyl_array_cfetch(f_other, loc);
    
    gkyl_rect_grid_cell_center(grid, iter.idx, xc);

    up->func(inp, up->ctx, iter.idx, xc, gain, weight);
    count += 1;
  }
  double effective_gamma = weight[0]/weight[1];
  double k = up->norm(flux, bc_param, effective_gamma);
  return k;
}

void gkyl_bc_emission_spectrum_release(struct gkyl_bc_emission_spectrum *up)
{
  // Release updater memory.
  gkyl_free(up);
}
