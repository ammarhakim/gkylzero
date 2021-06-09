extern "C" {
#include <gkyl_mom_calc.h>
#include <gkyl_util.h>
}

__global__ void gkyl_mom_calc_advance_cu_ker(const gkyl_mom_calc* mcalc,
  const struct gkyl_range phase_range, const struct gkyl_range conf_range,
  const struct gkyl_array* GKYL_RESTRICT fin, struct gkyl_array* GKYL_RESTRICT mout){

  unsigned long linc = blockIdx.x * blockDim.x + threadIdx.x;

  double xc[GKYL_MAX_DIM];
  int pidx[GKYL_MAX_DIM];
  gkyl_sub_range_inv_idx(&phase_range, linc, pidx);
  gkyl_rect_grid_cell_center(&mcalc->grid, pidx, xc);

  long linc_r = gkyl_range_idx(&phase_range, pidx);
  const double* fptr = (const double*) gkyl_array_cfetch(fin, linc_r);
  double momLocal[32]; // hard-coded to max confBasis.numBasis (3x p=3) for now.

  // reduce local f to local mom
  mcalc->momt->kernel(xc, mcalc->grid.dx, pidx, fptr, &momLocal[0]);

  // get conf-space linear index.
  int cidx[GKYL_MAX_DIM];
  for (unsigned int k = 0; k < 3; k++) {
    cidx[k] = pidx[k];
  }
  long lincConf = gkyl_range_idx(&conf_range, cidx);

  double* mptr = (double*) gkyl_array_fetch(mout, lincConf);
  for (int k = 0; k < mout->ncomp; ++k) {
     if (linc < phase_range.volume) atomicAdd(&mptr[k], momLocal[k]);
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
