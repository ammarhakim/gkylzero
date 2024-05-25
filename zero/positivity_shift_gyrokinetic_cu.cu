/* -*- c++ -*- */

extern "C" {
#include <gkyl_positivity_shift_gyrokinetic.h>
#include <gkyl_positivity_shift_gyrokinetic_priv.h>
#include <gkyl_array_ops.h>
#include <float.h>
}

// CUDA kernel to set device pointers to kernels.
__global__ static void
gkyl_pos_shift_gk_set_cu_ker_ptrs(struct gkyl_positivity_shift_gyrokinetic_kernels *kernels,
  struct gkyl_basis pbasis)
{
  int pdim = pbasis.ndim;
  enum gkyl_basis_type b_type = pbasis.b_type;
  int poly_order = pbasis.poly_order;

  switch (b_type) {
    case GKYL_BASIS_MODAL_GKHYBRID:
    case GKYL_BASIS_MODAL_SERENDIPITY:
      kernels->shift = pos_shift_gk_kern_list_shift_ser[pdim-2].kernels[poly_order-1];
      kernels->int_mom = pos_shift_gk_kern_list_intmom_ser[pdim-2].kernels[poly_order-1];
      break;
    default:
      assert(false);
  }
};

void
pos_shift_gk_choose_shift_kernel_cu(struct gkyl_positivity_shift_gyrokinetic_kernels *kernels,
  struct gkyl_basis pbasis)
{
  gkyl_pos_shift_gk_set_cu_ker_ptrs<<<1,1>>>(kernels, pbasis);
}

// Function borrowed from array_reduce_cu.cu.
__device__ static __forceinline__ double
pos_shift_atomicMax_double(double *address, double val)
{
  unsigned long long int ret = __double_as_longlong(*address);
  while(val > __longlong_as_double(ret))
  {
    unsigned long long int old = ret;
    if((ret = atomicCAS((unsigned long long int*)address, old, __double_as_longlong(val))) == old)
      break;
  }
  return __longlong_as_double(ret);
}

__global__ static void
gkyl_positivity_shift_gyrokinetic_advance_cu_ker(
  struct gkyl_positivity_shift_gyrokinetic_kernels *kers,
  const struct gkyl_rect_grid grid,
  const struct gkyl_range phase_range, const struct gkyl_range conf_range,
  double *ffloor, double ffloor_fac, double cellav_fac, double mass, const struct gkyl_array* GKYL_RESTRICT bmag, 
  struct gkyl_array* GKYL_RESTRICT distf, struct gkyl_array* GKYL_RESTRICT mom)
{
  double xc[GKYL_MAX_DIM];
  int pidx[GKYL_MAX_DIM];
  double distf_max = -DBL_MAX;

  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < phase_range.volume; tid += blockDim.x*gridDim.x) {
    gkyl_sub_range_inv_idx(&phase_range, tid, pidx);

    long plinidx = gkyl_range_idx(&phase_range, pidx);
    double *fptr = (double*) gkyl_array_fetch(distf, plinidx);

    const int max_num_basis = 80; // hard-coded to 4 * max confBasis.num_basis (3x p=2 Ser) for now.

    // Shift f if needed, and compute Deltaf.
    double Deltaf[max_num_basis];
    bool shiftedf = kers->shift(ffloor[0], fptr, Deltaf);

    distf_max = fmax(distf_max, fptr[0]);

    if (shiftedf) {
      gkyl_rect_grid_cell_center(&grid, pidx, xc);

      double momLocal[max_num_basis] = {0.0};

      long clinidx = gkyl_range_idx(&conf_range, pidx);
      double *mptr = (double*) gkyl_array_fetch(mom, clinidx);
      const double *bmagptr = (const double *) gkyl_array_cfetch(bmag, clinidx);

      kers->int_mom(xc, grid.dx, pidx, mass, bmagptr, Deltaf, momLocal);

      for (unsigned int k = 0; k < mom->ncomp; ++k) {
         if (tid < phase_range.volume)
           atomicAdd(&mptr[k], momLocal[k]);
      }
    }
  }

  pos_shift_atomicMax_double(ffloor, ffloor_fac * distf_max * cellav_fac);
}

void
gkyl_positivity_shift_gyrokinetic_advance_cu(gkyl_positivity_shift_gyrokinetic* up,
  const struct gkyl_range *phase_rng, const struct gkyl_range *conf_rng,
  struct gkyl_array *GKYL_RESTRICT distf, struct gkyl_array *GKYL_RESTRICT mom)
{
  int nblocks = phase_rng->nblocks, nthreads = phase_rng->nthreads;
  gkyl_array_clear_range(mom, 0.0, conf_rng);

  gkyl_positivity_shift_gyrokinetic_advance_cu_ker<<<nblocks, nthreads>>>
    (up->kernels, up->grid, *phase_rng, *conf_rng, up->ffloor, up->ffloor_fac, up->cellav_fac,
     up->mass, up->gk_geom->bmag->on_dev, distf->on_dev, mom->on_dev);
}
