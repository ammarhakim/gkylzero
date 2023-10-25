#include <gkyl_array_integrate.h>
#include <gkyl_array_integrate_priv.h>
#include <gkyl_alloc.h>
#include <assert.h>

struct gkyl_array_integrate*
gkyl_array_integrate_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *basis,
  int num_comp, enum gkyl_array_integrate_op op, bool use_gpu)
{
  // Allocate space for new updater.
  struct gkyl_array_integrate *up = gkyl_malloc(sizeof(struct gkyl_array_integrate));

  up->num_basis = basis->num_basis;
  up->num_comp = num_comp;
  up->use_gpu = use_gpu;

  int ndim = basis->ndim;
  up->vol = 1.0;
  if (basis->poly_order > 0) {
    for (unsigned d=0; d<ndim; ++d)
      up->vol *= op == GKYL_ARRAY_INTEGRATE_OP_SQ? grid->dx[d]/2.0 : grid->dx[d]/sqrt(2.0);
  } else {
    // for polyOrder = 0 no normalization is applied.
    for (unsigned d=0; d<ndim; ++d)
      up->vol *= grid->dx[d];
  }

  // Choose the kernel that performs the desired operation within the integral.
  up->kernels = gkyl_malloc(sizeof(struct gkyl_array_integrate_kernels));
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    up->kernels_cu = gkyl_cu_malloc(sizeof(struct gkyl_array_integrate_kernels));
    gkyl_array_integrate_choose_kernel_cu(op, up->kernels_cu);
  } else {
    gkyl_array_integrate_choose_kernel(op, up->kernels);
    assert(up->kernels->intker);
    up->kernels_cu = up->kernels;
  }
#else
  gkyl_array_integrate_choose_kernel(op, up->kernels);
  assert(up->kernels->intker);
  up->kernels_cu = up->kernels;
#endif

  return up;
}

void gkyl_array_integrate_advance(gkyl_array_integrate *up, const struct gkyl_array *arr,
  double weight, const struct gkyl_range *range, double *out)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_array_integrate_advance_cu(up, arr, weight, out);
    return;
  }
#endif

  for (int k=0; k<up->num_comp; k++) out[k] = 0;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {

    long linidx = gkyl_range_idx(range, iter.idx);
    const double *arr_d = (const double*) gkyl_array_cfetch(arr, linidx);

    up->kernels->intker(up->vol, up->num_comp, up->num_basis, arr_d, out);
  }

  for (int k=0; k<up->num_comp; k++) out[k] *= weight;
}

void gkyl_array_integrate_release(gkyl_array_integrate *up)
{
  // Release memory associated with this updater.
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    gkyl_cu_free(up->kernels_cu);
#endif
  gkyl_free(up->kernels);
  gkyl_free(up);
}
