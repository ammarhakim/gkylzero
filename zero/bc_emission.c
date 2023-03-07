#include <gkyl_bc_emission.h>
#include <gkyl_bc_emission_priv.h>
#include <gkyl_alloc.h>
#include <assert.h>

struct gkyl_bc_emission*
gkyl_bc_emission_new(int dir, enum gkyl_edge_loc edge, const struct gkyl_range *local_range_ext,
  const int *num_ghosts, enum gkyl_bc_emission_type bctype, const struct gkyl_basis *basis,
  int num_comp, int cdim, const double *bc_param, const struct gkyl_array *bc_field, bool use_gpu)
{
  // Allocate space for new updater.
  struct gkyl_bc_emission *up = gkyl_malloc(sizeof(struct gkyl_bc_emission));

  up->dir = dir;
  up->cdim = cdim;
  up->edge = edge;
  up->bctype = bctype;
  up->use_gpu = use_gpu;

  // Create the skin/ghost ranges.
  gkyl_skin_ghost_ranges(&up->skin_r, &up->ghost_r, dir, edge,
    local_range_ext, num_ghosts);

  switch (bctype) {
    case GKYL_BC_CONSTANT_GAIN:
      struct bc_gain_ctx *ctx = gkyl_malloc(sizeof(*ctx));
      ctx->gain = bc_param[0];
      up->func = bc_emission_gain;
      ctx->basis = basis;
      ctx->dir = dir;
      ctx->cdim = cdim;
      ctx->ncomp = num_comp;
  
      up->ctx = ctx;
      up->ctx_on_dev = up->ctx;
      break;
      
    default:
      assert(false);
      break;
  }
  return up;
}

/* Modeled after gkyl_array_flip_copy_to_buffer_fn */
void
gkyl_bc_emission_advance(const struct gkyl_bc_emission *up, struct gkyl_array *buff_arr, struct gkyl_array *f_arr)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_bc_emission_advance_cu(up, buff_arr, f_arr);
    return;
  }
#endif
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->skin_r);

  int dir = up->dir + up->cdim;

  int fidx[GKYL_MAX_DIM]; // flipped index
  struct gkyl_range buff_range;
  gkyl_range_init(&buff_range, up->skin_r.ndim, up->skin_r.lower, up->skin_r.upper);

  int uplo = up->skin_r.upper[dir] + up->skin_r.lower[dir];

  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&up->skin_r, iter.idx);
    
    gkyl_copy_int_arr(up->skin_r.ndim, iter.idx, fidx);
    fidx[dir] = uplo - iter.idx[dir];
    
    long count = gkyl_range_idx(&buff_range, fidx);

    const double *inp = gkyl_array_cfetch(f_arr, loc);
    double *out = flat_fetch(buff_arr->data, f_arr->esznc*count);
    up->func(out, inp, up->ctx, iter.idx, fidx);
  }
  gkyl_array_copy_from_buffer(f_arr, buff_arr->data, up->ghost_r);
}

void gkyl_bc_emission_release(struct gkyl_bc_emission *up)
{
  // Release updater memory.
  gkyl_free(up);
}
