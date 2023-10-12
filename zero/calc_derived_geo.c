#include <gkyl_calc_derived_geo.h>
#include <gkyl_calc_derived_geo_priv.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_range.h>

#include <gkyl_array_ops_priv.h>

gkyl_calc_derived_geo*
gkyl_calc_derived_geo_new(const struct gkyl_basis *cbasis, const struct gkyl_rect_grid *grid, bool use_gpu)
{
  gkyl_calc_derived_geo *up = gkyl_malloc(sizeof(gkyl_calc_derived_geo));
  up->cdim = cbasis->ndim;
  up->cnum_basis = cbasis->num_basis;
  up->poly_order = cbasis->poly_order;
  up->grid = grid;
  up->use_gpu = use_gpu;
  up->kernel = derived_geo_choose_kernel(up->cdim, cbasis->b_type, up->poly_order);
  up->adjustment_kernel = adjust_bmag_choose_kernel(up->cdim, cbasis->b_type, up->poly_order);
  return up;
}


void
gkyl_calc_derived_geo_advance(const gkyl_calc_derived_geo *up, const struct gkyl_range *crange, struct gkyl_array *gFld, struct gkyl_array *bmagFld, struct gkyl_array *jFld, struct gkyl_array *jinvFld, struct gkyl_array *grFld, struct gkyl_array *biFld, struct gkyl_array *cmagFld, struct gkyl_array *jtotFld, struct gkyl_array *jtotinvFld, struct gkyl_array *bmaginvFld, struct gkyl_array *bmaginvsqFld, struct gkyl_array *gxxJFld,  struct gkyl_array *gxyJFld, struct gkyl_array *gyyJFld)
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, crange);
  while (gkyl_range_iter_next(&iter)) {
    //printf("iter.idx = %d,%d,%d\n", iter.idx[0],iter.idx[1],iter.idx[2]);
    long loc = gkyl_range_idx(crange, iter.idx);
    const double *gij = gkyl_array_cfetch(gFld, loc);
    const double *bmag_i = gkyl_array_cfetch(bmagFld, loc);
    double *j_i = gkyl_array_fetch(jFld, loc);
    double *jinv_i = gkyl_array_fetch(jinvFld, loc);
    double *grij = gkyl_array_fetch(grFld, loc);
    double *bi_i = gkyl_array_fetch(biFld, loc);
    double *cmag_i = gkyl_array_fetch(cmagFld, loc);
    double *jtot_i = gkyl_array_fetch(jtotFld, loc);
    double *jtotinv_i = gkyl_array_fetch(jtotinvFld, loc);
    double *bmaginv_i = gkyl_array_fetch(bmaginvFld, loc);
    double *bmaginvsq_i = gkyl_array_fetch(bmaginvsqFld, loc);
    double *gxxJ_i= gkyl_array_fetch(gxxJFld, loc);
    double *gxyJ_i= gkyl_array_fetch(gxyJFld, loc);
    double *gyyJ_i= gkyl_array_fetch(gyyJFld, loc);

    up->kernel(gij, bmag_i, j_i, jinv_i, grij, bi_i, cmag_i, jtot_i, jtotinv_i, bmaginv_i, bmaginvsq_i, gxxJ_i, gxyJ_i, gyyJ_i);
  }

  struct gkyl_range_iter cmag_iter;
  gkyl_range_iter_init(&iter, crange);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(crange, iter.idx);
    const double *bmag_i = gkyl_array_cfetch(bmagFld, loc); // bmag we want to change
    const double *j_i = gkyl_array_cfetch(jFld, loc);
    double *gij_i = gkyl_array_fetch(gFld, loc);
    double *gzz_i = &gij_i[(2*up->cdim - 1)*up->cnum_basis];
    const double *cmag_i = gkyl_array_cfetch(cmagFld, loc);

    cmag_iter.idx[0] = iter.idx[0];
    for(int i = 1; i <up->cdim; i++)
      cmag_iter.idx[i] = (crange->upper[i] - crange->lower[i])/2 + 1; // for using middle node
    long cmag_loc = gkyl_range_idx(crange, cmag_iter.idx);

    const double *cmag_ref_i = gkyl_array_cfetch(cmagFld, cmag_loc); // cmag we want for yz plane
    // now call a kernel that takes j, gzz, and cmag as inputs and calculates bmag
    up->adjustment_kernel(cmag_i, cmag_ref_i, gzz_i, j_i, bmag_i, gij_i);

  }

  // recalculate the rest now
  gkyl_range_iter_init(&iter, crange);
  while (gkyl_range_iter_next(&iter)) {
    //printf("iter.idx = %d,%d,%d\n", iter.idx[0],iter.idx[1],iter.idx[2]);
    long loc = gkyl_range_idx(crange, iter.idx);
    const double *gij = gkyl_array_cfetch(gFld, loc);
    const double *bmag_i = gkyl_array_cfetch(bmagFld, loc);
    double *j_i = gkyl_array_fetch(jFld, loc);
    double *jinv_i = gkyl_array_fetch(jinvFld, loc);
    double *grij = gkyl_array_fetch(grFld, loc);
    double *bi_i = gkyl_array_fetch(biFld, loc);
    double *cmag_i = gkyl_array_fetch(cmagFld, loc);

    double *jtot_i = gkyl_array_fetch(jtotFld, loc);
    double *jtotinv_i = gkyl_array_fetch(jtotinvFld, loc);
    double *bmaginv_i = gkyl_array_fetch(bmaginvFld, loc);
    double *bmaginvsq_i = gkyl_array_fetch(bmaginvsqFld, loc);
    double *gxxJ_i= gkyl_array_fetch(gxxJFld, loc);
    double *gxyJ_i= gkyl_array_fetch(gxyJFld, loc);
    double *gyyJ_i= gkyl_array_fetch(gyyJFld, loc);

    up->kernel(gij, bmag_i, j_i, jinv_i, grij, bi_i, cmag_i, jtot_i, jtotinv_i, bmaginv_i, bmaginvsq_i,gxxJ_i, gxyJ_i, gyyJ_i);
    //printf("gij[1] = %g\n",gij[1]);
    //printf("j[1] = %g\n",j_i[0]);
  }

}

void
gkyl_calc_derived_geo_release(gkyl_calc_derived_geo* up)
{
  gkyl_free(up);
}
