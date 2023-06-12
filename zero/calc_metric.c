#include <gkyl_calc_metric.h>
#include <gkyl_calc_metric_priv.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_range.h>

#include <gkyl_array_ops_priv.h>

gkyl_calc_metric*
gkyl_calc_metric_new(const struct gkyl_basis *cbasis, struct gkyl_rect_grid grid, bool use_gpu)
{
  gkyl_calc_metric *up = gkyl_malloc(sizeof(gkyl_calc_metric));
  up->cdim = cbasis->ndim;
  up->cnum_basis = cbasis->num_basis;
  up->poly_order = cbasis->poly_order;
  up->grid = grid;
  up->use_gpu = use_gpu;
  up->kernel = metric_choose_kernel(up->cdim, cbasis->b_type, up->poly_order);
  return up;
}

void
gkyl_calc_metric_advance(const gkyl_calc_metric *up, const struct gkyl_range *crange, struct gkyl_array *XYZ, struct gkyl_array *gFld)
{
  const double **xyz = gkyl_malloc(7*sizeof(double*));
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, crange);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(crange, iter.idx);
    double *gij = gkyl_array_fetch(gFld, loc);
    xyz[0] = gkyl_array_cfetch(XYZ,loc);
    int count = 1;
    int idx_temp[3] = {iter.idx[0], iter.idx[1], iter.idx[2]};
    for(int i = 0; i<3; i++){
      idx_temp[0] = iter.idx[0];
      idx_temp[1] = iter.idx[1];
      idx_temp[2] = iter.idx[2];
      for(int j = -1; j<3; j+=2){
        idx_temp[i] = iter.idx[i] + j;
        loc = gkyl_range_idx(crange, idx_temp);
        xyz[count] = gkyl_array_cfetch(XYZ, loc);
        count = count+1;
      }
    }
    up->kernel(xyz,gij);
  }

  double scale_factor[6];
  int count = 0;
  for(int i=0; i<3; i++){
    for(int j=i; j<3; j++){
      scale_factor[count] = 4.0/(up->grid.dx[i]*up->grid.dx[j]);
      count = count+1;
    }
  }

  gkyl_range_iter_init(&iter, crange);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(crange, iter.idx);
    double *gij = gkyl_array_fetch(gFld, loc);
    for(int i=0; i<6; i++){
      double *gcomp = &gij[i*(up->cnum_basis)];
      for(int j=0; j < (up->cnum_basis); j++){
        gcomp[j] = gcomp[j]*scale_factor[i];
      }
    }
  }
}

void
gkyl_calc_metric_release(gkyl_calc_metric* up)
{
  gkyl_free(up);
}
