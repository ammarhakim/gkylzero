#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_rect_grid.h>
#include <gkyl_deflate_zsurf.h>
#include <gkyl_deflate_zsurf_priv.h>



struct gkyl_deflate_zsurf* gkyl_deflate_zsurf_new(const struct gkyl_basis *cbasis,const struct gkyl_basis *deflated_cbasis, const struct gkyl_rect_grid *grid, const struct gkyl_rect_grid *deflated_grid, int edge, bool use_gpu){

  gkyl_deflate_zsurf *up = gkyl_malloc(sizeof(gkyl_deflate_zsurf));
  up->grid = grid;
  up->deflated_grid = deflated_grid;
  up->basis = cbasis;
  up->deflated_basis = deflated_cbasis;
  up->kernel = deflate_zsurf_choose_kernel(cbasis->b_type, edge, cbasis->poly_order); // edge = 0,1 = lo, up

  return up;
}


void gkyl_deflate_zsurf_advance(const gkyl_deflate_zsurf *up, int zidx, const struct gkyl_range *range, const struct gkyl_range* deflated_range, const struct gkyl_array *field, struct gkyl_array *deflated_field, int ncomp){

  //we loop along

  struct gkyl_range update_range;
  int do_idx[2];
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, deflated_range);

  while(gkyl_range_iter_next(&iter)){
      do_idx[0] = iter.idx[0];
      do_idx[1] = zidx;

      long loc = gkyl_range_idx(range, do_idx);
      const double *fld = gkyl_array_cfetch(field, loc);

      long loc_deflated = gkyl_range_idx(deflated_range, iter.idx);
      double *fld_deflated = gkyl_array_fetch(deflated_field, loc_deflated);
      for(int c = 0; c<ncomp; c++){
        up->kernel(&fld[c*up->basis->num_basis], &fld_deflated[c*up->deflated_basis->num_basis]);
      }
 }

}

void gkyl_deflate_zsurf_release(gkyl_deflate_zsurf* up){
  gkyl_free(up);
}


