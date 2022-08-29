#include <gkyl_bc_external.h>
#include <gkyl_bc_basic_priv.h>
#include <gkyl_alloc.h>
#include <assert.h>

void
gkyl_bc_external_furman_pivi_advance(void *data, const struct gkyl_array *arr,
  int dir, struct gkyl_range range, struct gkyl_array_copy_func *cf)
{
  // gkyl_array_clear_range(out, 0.0, &range);
  
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);

  int fidx[GKYL_MAX_DIM]; // flipped index
  struct gkyl_range buff_range;
  gkyl_range_init(&buff_range, range.ndim, range.lower, range.upper);

  int uplo = range.upper[dir]+range.lower[dir];
  struct dg_bc_ctx *mc = (struct dg_bc_ctx*) cf->ctx;

  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);
    const double *inp = gkyl_array_cfetch(arr, loc);
    mc->in_idx = iter.idx[dir];

    struct gkyl_range_iter out_iter;
    gkyl_range_iter_init(&out_iter, &range);
    
    while (gkyl_range_iter_next(&out_iter)) {

      gkyl_copy_int_arr(range.ndim, out_iter.idx, fidx);
      fidx[dir] = uplo - out_iter.idx[dir];
    
      long count = gkyl_range_idx(&buff_range, fidx);
      mc->out_idx = out_iter.idx[dir];

      double *out = flat_fetch(data, arr->esznc*count);
      cf->func(arr->ncomp, out, inp, cf->ctx);
    }
  }  
}
