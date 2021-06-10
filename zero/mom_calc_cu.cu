extern "C" {
#include <gkyl_mom_calc.h>
#include <gkyl_util.h>
}

__global__ void gkyl_mom_calc_advance_cu_ker(const gkyl_mom_calc* mcalc,
  const struct gkyl_range phase_range, const struct gkyl_range conf_range,
  const struct gkyl_array* GKYL_RESTRICT fin, struct gkyl_array* GKYL_RESTRICT mout){

  double xc[GKYL_MAX_DIM];
  int pidx[GKYL_MAX_DIM], cidx[GKYL_MAX_CDIM];

  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < phase_range.volume;
      tid += blockDim.x*gridDim.x)
  {
    gkyl_sub_range_inv_idx(&phase_range, tid, pidx);
    gkyl_rect_grid_cell_center(&mcalc->grid, pidx, xc);

    long lincP = gkyl_range_idx(&phase_range, pidx);
    const double* fptr = (const double*) gkyl_array_cfetch(fin, lincP);
    double momLocal[96]; // hard-coded to max confBasis.numBasis (3x p=3 Ser) for now.
    for (unsigned int k=0; k<96; ++k){
      momLocal[k] = 0.0;
    }

    // reduce local f to local mom
    mcalc->momt->kernel(xc, mcalc->grid.dx, pidx, fptr, &momLocal[0]);

    // get conf-space linear index.
    for (unsigned int k = 0; k < conf_range.ndim; k++) {
      cidx[k] = pidx[k];
    }
    long lincC = gkyl_range_idx(&conf_range, cidx);

    double* mptr = (double*) gkyl_array_fetch(mout, lincC);

    for (unsigned int k = 0; k < mout->ncomp; ++k) {
       if (tid < phase_range.volume) atomicAdd(&mptr[k], momLocal[k]);
    };
  };

}

void
gkyl_mom_calc_advance_cu(const gkyl_mom_calc* mcalc,
  const struct gkyl_range phase_range, const struct gkyl_range conf_range,
  const struct gkyl_array* GKYL_RESTRICT fin, struct gkyl_array* GKYL_RESTRICT mout)
{

  int nblocks = phase_range.nblocks, nthreads = phase_range.nthreads;

  gkyl_mom_calc_advance_cu_ker<<<nblocks, nthreads>>>(mcalc, phase_range, conf_range, fin, mout);

}
