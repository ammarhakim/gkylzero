#include <gkyl_calc_metric.h>
#include <gkyl_calc_metric_priv.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_range.h>

#include <gkyl_array_ops_priv.h>

gkyl_calc_metric*
gkyl_calc_metric_new(const struct gkyl_basis *cbasis, struct gkyl_rect_grid *grid, int *bcs, bool use_gpu)
{
  gkyl_calc_metric *up = gkyl_malloc(sizeof(gkyl_calc_metric));
  up->cdim = cbasis->ndim;
  up->cnum_basis = cbasis->num_basis;
  up->poly_order = cbasis->poly_order;
  up->grid = grid;
  up->use_gpu = use_gpu;
  up->num_cells = up->grid->cells;
  up->bcs = bcs;
  up->kernels = metric_choose_kernel(up->cdim, up->poly_order, up->bcs);

  return up;
}

void
gkyl_calc_metric_advance(const gkyl_calc_metric *up, const struct gkyl_range *crange, struct gkyl_array *XYZ, struct gkyl_array *gFld)
{
  const double **xyz = gkyl_malloc((1+2*up->cdim)*sizeof(double*));
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, crange);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(crange, iter.idx);
    double *gij = gkyl_array_fetch(gFld, loc);
    xyz[0] = gkyl_array_cfetch(XYZ,loc);
    int count = 1;
    int idx_temp[up->cdim];
    for(int l = 0; l<up->cdim; l++){idx_temp[l] = iter.idx[l]; }
    for(int i = 0; i<up->cdim; i++){
      for(int l = 0; l<up->cdim; l++){idx_temp[l] = iter.idx[l]; }
      for(int j = -1; j<3; j+=2){
        idx_temp[i] = iter.idx[i] + j;
        loc = gkyl_range_idx(crange, idx_temp);
        xyz[count] = gkyl_array_cfetch(XYZ, loc);
        count = count+1;
      }
    }
    int linker_idx = idx_to_inloup_ker(up->cdim, up->num_cells, iter.idx);
    //linker_idx = 0; // to always use two sided
    up->kernels.kernels[linker_idx](xyz,gij);
  }

  double scale_factor[up->cdim * (up->cdim+1)/2];
  int count = 0;
  for(int i=0; i<up->cdim; i++){
    for(int j=i; j<up->cdim; j++){
      scale_factor[count] = 4.0/(up->grid->dx[i]*up->grid->dx[j]);
      count = count+1;
    }
  }

  gkyl_range_iter_init(&iter, crange);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(crange, iter.idx);
    double *gij = gkyl_array_fetch(gFld, loc);
    for(int i=0; i<up->cdim * (up->cdim+1)/2; i++){
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
