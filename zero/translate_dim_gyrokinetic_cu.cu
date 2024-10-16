/* -*- c++ -*- */

extern "C" {
#include <gkyl_translate_dim_gyrokinetic.h>
#include <gkyl_translate_dim_gyrokinetic_priv.h>
#include <gkyl_array_ops.h>
#include <float.h>
}

// CUDA kernel to set device pointers to kernels.
__global__ static void
gkyl_trans_dim_gk_set_cu_ker_ptrs(struct gkyl_translate_dim_gyrokinetic_kernels *kernels,
  int cdim_do, struct gkyl_basis pbasis_do, int cdim_tar, struct gkyl_basis pbasis_tar)
{
  enum gkyl_basis_type basis_type = pbasis_tar.b_type;
  int poly_order = pbasis_tar.poly_order;

  switch (basis_type) {
    case GKYL_BASIS_MODAL_GKHYBRID:
    case GKYL_BASIS_MODAL_SERENDIPITY:
      kernels->translate = trans_dim_gk_kern_list_shift_ser[cdim_do+cdim_tar-3].kernels[poly_order-1];
      break;
    default:
      assert(false);
  }
};

void
trans_dim_gk_choose_shift_kernel_cu(struct gkyl_translate_dim_gyrokinetic_kernels *kernels,
  int cdim_do, struct gkyl_basis pbasis_do, int cdim_tar, struct gkyl_basis pbasis_tar)
{
  gkyl_trans_dim_gk_set_cu_ker_ptrs<<<1,1>>>(kernels, cdim_do, pbasis_do, cdim_tar, pbasis_tar);
}

__global__ static void
gkyl_translate_dim_gyrokinetic_advance_cu_ker(int cdim_do, int cdim_tar, int vdim_do, int vdim_tar,
  struct gkyl_translate_dim_gyrokinetic_kernels *kernels,
  const struct gkyl_range phase_rng_do, const struct gkyl_range phase_rng_tar,
  const struct gkyl_array *GKYL_RESTRICT fdo, struct gkyl_array *GKYL_RESTRICT ftar)
{
  int pidx_do[GKYL_MAX_DIM] = {-1};
  int pidx_tar[GKYL_MAX_DIM] = {-1};

  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < phase_rng_tar.volume; tid += blockDim.x*gridDim.x) {
    gkyl_sub_range_inv_idx(&phase_rng_tar, tid, pidx_tar);

    // Translate the target idx to the donor idx:
    for (int d=0; d<cdim_do-1; d++) pidx_do[d] = pidx_tar[d]; 
    pidx_do[cdim_do-1] = pidx_tar[cdim_tar-1]; 
    for (int d=0; d<vdim_do; d++) pidx_do[cdim_do+d] = pidx_tar[cdim_tar+d]; 

    long plinidx_do = gkyl_range_idx(&phase_rng_do, pidx_do);
    long plinidx_tar = gkyl_range_idx(&phase_rng_tar, pidx_tar);

    const double *fdo_c = (const double *) gkyl_array_cfetch(fdo, plinidx_do);
    double *ftar_c = (double *) gkyl_array_fetch(ftar, plinidx_tar);

    kernels->translate(fdo_c, ftar_c);

  }
}

void
gkyl_translate_dim_gyrokinetic_advance_cu(gkyl_translate_dim_gyrokinetic* up,
  const struct gkyl_range *phase_rng_do, const struct gkyl_range *phase_rng_tar,
  const struct gkyl_array *GKYL_RESTRICT fdo, struct gkyl_array *GKYL_RESTRICT ftar)
{
  int nblocks = phase_rng_tar->nblocks, nthreads = phase_rng_tar->nthreads;

  gkyl_translate_dim_gyrokinetic_advance_cu_ker<<<nblocks, nthreads>>>
    (up->cdim_do, up->cdim_tar, up->vdim_do, up->vdim_tar, up->kernels,
     *phase_rng_do, *phase_rng_tar, fdo->on_dev, ftar->on_dev);
}
