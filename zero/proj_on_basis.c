#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>

struct gkyl_proj_on_basis {
  struct gkyl_rect_grid grid;
  int num_quad; // number of quadrature points to use in each direction
  int num_ret_vals; // number of values returned by eval function
  evalf_t eval; // function to project
  void *ctx; // evaluation context

  int numBasis; // number of basis functions
  int tot_quad; // total number of quadrature points
  struct gkyl_array *ordinates; // ordinates for quadrature
  struct gkyl_array *weights; // weights for quadrature
  struct gkyl_array *basis_at_ords; // basis functions at ordinates
};

gkyl_proj_on_basis*
gkyl_proj_on_basis_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *basis,
  int num_quad, int num_ret_vals, evalf_t eval, void *ctx)
{
  gkyl_proj_on_basis *up = gkyl_malloc(sizeof(gkyl_proj_on_basis));

  up->grid = *grid;
  up->num_quad = num_quad;
  up->num_ret_vals = num_ret_vals;
  up->eval = eval;
  up->ctx = ctx;
  up->numBasis = basis->numBasis;

  double ordinates1[num_quad], weights1[num_quad];

  if (num_quad <= gkyl_gauss_max) {
    // use pre-computed values if possible (these are more accurate
    // than computing them on the fly)
    memcpy(ordinates1, gkyl_gauss_ordinates[num_quad], sizeof(double[num_quad]));
    memcpy(weights1, gkyl_gauss_weights[num_quad], sizeof(double[num_quad]));
  }
  else {
    gkyl_gauleg(-1, 1, ordinates1, weights1, num_quad);
  }

  // create range to loop over quadrature points
  int qshape[GKYL_MAX_DIM];
  for (int i=0; i<grid->ndim; ++i) qshape[i] = num_quad;
  struct gkyl_range qrange;
  gkyl_range_init_from_shape(&qrange, grid->ndim, qshape);

  int tot_quad = up->tot_quad = qrange.volume;

  // create ordinates and weights for multi-D quadrature  
  up->ordinates = gkyl_array_new(GKYL_DOUBLE, grid->ndim, tot_quad);
  up->weights = gkyl_array_new(GKYL_DOUBLE, 1, tot_quad);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &qrange);

  while (gkyl_range_iter_next(&iter)) {
    long node = gkyl_range_idx(&qrange, iter.idx);
    
    // set ordinates
    double *ord = gkyl_array_fetch(up->ordinates, node);
    for (int i=0; i<grid->ndim; ++i)
      ord[i] = ordinates1[iter.idx[i]-qrange.lower[i]];
    
    // set weights
    double *wgt = gkyl_array_fetch(up->weights, node);
    wgt[0] = 1.0;
    for (int i=0; i<qrange.ndim; ++i)
      wgt[0] *= weights1[iter.idx[i]-qrange.lower[i]];
  }

  // pre-compute basis functions at ordinates
  up->basis_at_ords = gkyl_array_new(GKYL_DOUBLE, basis->numBasis, tot_quad);
  for (int n=0; n<tot_quad; ++n)
    basis->eval(gkyl_array_fetch(up->ordinates, n),
      gkyl_array_fetch(up->basis_at_ords, n));

  return up;
}

static inline void
comp_to_phys(int ndim, const double *eta,
  const double * restrict dx, const double * restrict xc,
  double* restrict xout)
{
  for (int d=0; d<ndim; ++d) xout[d] = 0.5*dx[d]*eta[d]+xc[d];
}

static void
proj_on_basis(const gkyl_proj_on_basis *up, const struct gkyl_array *fun_at_ords, double* f)
{
  int numBasis = up->numBasis;
  int tot_quad = up->tot_quad;
  int num_ret_vals = up->num_ret_vals;

  const double *restrict weights = up->weights->data;
  const double *restrict basis_at_ords = up->basis_at_ords->data;
  const double *restrict func_at_ords = fun_at_ords->data;

  // arrangement of f is as:
  // c0[0], c0[1], ... c1[0], c1[1], ....
  // where c0, c1, ... are components of f (num_ret_vals)
  int offset = 0;
  for (int n=0; n<num_ret_vals; ++n) {
    for (int k=0; k<numBasis; ++k) f[offset+k] = 0.0;

    for (int imu=0; imu<tot_quad; ++imu) {
      double tmp = weights[imu]*func_at_ords[n+num_ret_vals*imu];
      for(int k=0; k<numBasis; ++k)
        f[offset+k] += tmp*basis_at_ords[k+numBasis*imu];
    }
    offset += numBasis;
  }
}

void
gkyl_proj_on_basis_advance(const gkyl_proj_on_basis *up,
  double tm, const struct gkyl_range *update_range, struct gkyl_array *arr)
{
  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM];

  int num_ret_vals = up->num_ret_vals;
  int tot_quad = up->tot_quad;
  struct gkyl_array *fun_at_ords = gkyl_array_new(GKYL_DOUBLE, num_ret_vals, tot_quad);
  
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_range);
  
  while (gkyl_range_iter_next(&iter)) {
    gkyl_rect_grid_cell_center(&up->grid, iter.idx, xc);

    for (int i=0; i<tot_quad; ++i) {
      comp_to_phys(up->grid.ndim, gkyl_array_cfetch(up->ordinates, i),
        up->grid.dx, xc, xmu);
      up->eval(tm, xmu, gkyl_array_fetch(fun_at_ords, i), up->ctx);
    }

    long lidx = gkyl_range_idx(update_range, iter.idx);
    proj_on_basis(up, fun_at_ords, gkyl_array_fetch(arr, lidx));
  }

  gkyl_array_release(fun_at_ords);
}

void
gkyl_proj_on_basis_release(gkyl_proj_on_basis* up)
{
  gkyl_array_release(up->ordinates);
  gkyl_array_release(up->weights);
  gkyl_array_release(up->basis_at_ords);
  gkyl_free(up);
}
