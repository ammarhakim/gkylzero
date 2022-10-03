/* -*- c++ -*- */

extern "C" {
#include <gkyl_spitzer_coll_freq.h>
#include <gkyl_spitzer_coll_freq_priv.h>
#include <gkyl_const.h>
#include <gkyl_range.h>
}

__global__ static void
gkyl_spitzer_coll_freq_advance_normnu_cu_ker(int num_quad, 
  const struct gkyl_range range, const struct gkyl_array* GKYL_RESTRICT basis_at_ords,
  const struct gkyl_array* GKYL_RESTRICT weights, const struct gkyl_array* GKYL_RESTRICT vtSqSelf,
  const struct gkyl_array* GKYL_RESTRICT m0Other, const struct gkyl_array* GKYL_RESTRICT vtSqOther,
  double normNu, struct gkyl_array* GKYL_RESTRICT nuOut)
{
  int num_basis = basis_at_ords->ncomp;
  int tot_quad = basis_at_ords->size;

  int idx[3];

  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < range.volume; tid += blockDim.x*gridDim.x) {

    gkyl_sub_range_inv_idx(&range, tid, idx);
    long linidx = gkyl_range_idx(&range, idx);

    const double *vtSqSelf_d = (const double *) gkyl_array_cfetch(vtSqSelf, linidx);
    const double *m0Other_d = (const double *) gkyl_array_cfetch(m0Other, linidx);
    const double *vtSqOther_d = (const double *) gkyl_array_cfetch(vtSqOther, linidx);

    double *nuOut_d = (double *) gkyl_array_fetch(nuOut, linidx);
    for (int k=0; k<num_basis; ++k) nuOut_d[k] = 0.0;

    // Compute expansion coefficients of nuOut using quadrature.
    const double *w_d = (const double *) weights->data;
    const double *bo_d = (const double *) basis_at_ords->data;

    for (int n=0; n<tot_quad; ++n) {

      const double *b_ord = (const double *) gkyl_array_cfetch(basis_at_ords, n);

      // Evaluate densities and thermal speeds (squared) at quad point.
      double vtSqSelf_q=0., m0Other_q=0., vtSqOther_q=0.;
      for (int k=0; k<num_basis; ++k) {
        vtSqSelf_q += vtSqSelf_d[k]*b_ord[k];
        m0Other_q += m0Other_d[k]*b_ord[k];
        vtSqOther_q += vtSqOther_d[k]*b_ord[k];
      }

      double nu_o = normNu*m0Other_q/pow(sqrt(vtSqSelf_q+vtSqOther_q),3);

      double tmp = w_d[n]*nu_o;
      for (int k=0; k<num_basis; ++k)
        nuOut_d[k] += tmp*bo_d[k+num_basis*n];
    }

  }
}

void
gkyl_spitzer_coll_freq_advance_normnu_cu(const gkyl_spitzer_coll_freq *up,
  const struct gkyl_range *range, const struct gkyl_array *vtSqSelf,
  const struct gkyl_array *m0Other, const struct gkyl_array *vtSqOther,
  double normNu, struct gkyl_array *nuOut)
{
  int nblocks = range->nblocks, nthreads = range->nthreads;
  gkyl_spitzer_coll_freq_advance_normnu_cu_ker<<<nblocks, nthreads>>>
    (up->num_quad, *range, up->basis_at_ords->on_dev,
     up->weights->on_dev, vtSqSelf->on_dev,
     m0Other->on_dev, vtSqOther->on_dev, normNu, nuOut->on_dev);
}
