#include <gkyl_translate_dim_gyrokinetic.h>
#include <gkyl_translate_dim_gyrokinetic_priv.h>
#include <gkyl_alloc.h>

struct gkyl_translate_dim_gyrokinetic*
gkyl_translate_dim_gyrokinetic_new(int cdim_do, struct gkyl_basis pbasis_do,
  int cdim_tar, struct gkyl_basis pbasis_tar, bool use_gpu)
{
  // Allocate space for new updater.
  struct gkyl_translate_dim_gyrokinetic *up = gkyl_malloc(sizeof(*up));

  up->use_gpu = use_gpu;
  up->cdim_do = cdim_do;
  up->cdim_tar = cdim_tar;
  up->vdim_do = pbasis_do.ndim - cdim_do;
  up->vdim_tar = pbasis_tar.ndim - cdim_tar;

  // Perform some basic checks:
  assert(cdim_tar > 1 );
  assert(up->vdim_do == up->vdim_tar);
  assert(pbasis_do.b_type == pbasis_tar.b_type);
  assert(pbasis_do.poly_order == pbasis_tar.poly_order);

  if (!use_gpu)
    up->kernels = gkyl_malloc(sizeof(struct gkyl_translate_dim_gyrokinetic_kernels));
#ifdef GKYL_HAVE_CUDA
  if (use_gpu)
    up->kernels = gkyl_cu_malloc(sizeof(struct gkyl_translate_dim_gyrokinetic_kernels));
#endif

  // Choose kernels that translates the DG coefficients.
  trans_dim_gk_choose_kernel(up->kernels, cdim_do, pbasis_do, cdim_tar, pbasis_tar, use_gpu);

  return up;
}

void
gkyl_translate_dim_gyrokinetic_advance(gkyl_translate_dim_gyrokinetic* up,
  const struct gkyl_range *phase_rng_do, const struct gkyl_range *phase_rng_tar,
  const struct gkyl_array *GKYL_RESTRICT fdo, struct gkyl_array *GKYL_RESTRICT ftar)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_translate_dim_gyrokinetic_advance_cu(up, phase_rng_do, phase_rng_tar, fdo, ftar);
    return;
  }
#endif

  // Perform some basic checks:
  for (int d=0; d<up->cdim_do; d++)
    assert(phase_rng_do->lower[d] == phase_rng_tar->lower[d]);
  for (int d=0; d<up->vdim_do; d++)
    assert(phase_rng_do->lower[up->cdim_do+d] == phase_rng_tar->lower[up->cdim_tar+d]);

  int pidx_do[GKYL_MAX_DIM] = {-1};

  struct gkyl_range_iter phase_iter;
  gkyl_range_iter_init(&phase_iter, phase_rng_tar);
  while (gkyl_range_iter_next(&phase_iter)) {

    // Translate the target idx to the donor idx:
    for (int d=0; d<up->cdim_do; d++) pidx_do[d] = phase_iter.idx[d]; 
    for (int d=0; d<up->vdim_do; d++) pidx_do[up->cdim_do+d] = phase_iter.idx[up->cdim_tar+d]; 

    long plinidx_tar = gkyl_range_idx(phase_rng_tar, phase_iter.idx);
    long plinidx_do = gkyl_range_idx(phase_rng_do, pidx_do);

    const double *fdo_c = gkyl_array_cfetch(fdo, plinidx_do);
    double *ftar_c = gkyl_array_fetch(ftar, plinidx_tar);

    up->kernels->translate(fdo_c, ftar_c);
  }
}

void
gkyl_translate_dim_gyrokinetic_release(gkyl_translate_dim_gyrokinetic* up)
{
  // Release memory associated with this updater.
  if (!up->use_gpu)
    gkyl_free(up->kernels);
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    gkyl_cu_free(up->kernels);
#endif
  gkyl_free(up);
}
