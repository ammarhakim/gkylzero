/* -*- c++ -*- */

extern "C" {
#include <gkyl_spitzer_coll_freq.h>
#include <gkyl_spitzer_coll_freq_priv.h>
#include <gkyl_const.h>
#include <gkyl_range.h>
}

GKYL_CU_D static inline void
calc_nu_cu(const struct gkyl_array* GKYL_RESTRICT basis_at_ords, const struct gkyl_array* GKYL_RESTRICT weights,
  const struct gkyl_array* GKYL_RESTRICT vtSqSelf_d, const struct gkyl_array* GKYL_RESTRICT m0Other_d,
  const struct gkyl_array* GKYL_RESTRICT vtSqOther_d, double normNu, long linidx,
  struct gkyl_array* GKYL_RESTRICT nuOut_d)
{
  // Perform the multiplication of normNu*n_r/(v_ts^2+v_tr^2)^(3/2) via
  // quadrature in one cell.
  int num_basis = basis_at_ords->ncomp;
  int tot_quad = basis_at_ords->size;

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

    double nu_o = ((m0Other_q<0.) || (vtSqSelf_q<0.) || (vtSqOther_q<0.)) ?
      1.e-40 : normNu*m0Other_q/pow(sqrt(vtSqSelf_q+vtSqOther_q),3);

    double tmp = w_d[n]*nu_o;
    for (int k=0; k<num_basis; ++k)
      nuOut_d[k] += tmp*bo_d[k+num_basis*n];
  }
}


__global__ static void
gkyl_spitzer_coll_freq_advance_normnu_cu_ker(const struct gkyl_range range,
  const struct gkyl_array* GKYL_RESTRICT basis_at_ords, const struct gkyl_array* GKYL_RESTRICT weights,
  const struct gkyl_array* GKYL_RESTRICT vtSqSelf, const struct gkyl_array* GKYL_RESTRICT m0Other,
  const struct gkyl_array* GKYL_RESTRICT vtSqOther, double normNu, struct gkyl_array* GKYL_RESTRICT nuOut)
{
  int idx[3];
  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < range.volume; tid += blockDim.x*gridDim.x) {

    gkyl_sub_range_inv_idx(&range, tid, idx);
    long linidx = gkyl_range_idx(&range, idx);

    const double *vtSqSelf_d = (const double *) gkyl_array_cfetch(vtSqSelf, linidx);
    const double *m0Other_d = (const double *) gkyl_array_cfetch(m0Other, linidx);
    const double *vtSqOther_d = (const double *) gkyl_array_cfetch(vtSqOther, linidx);

    calc_nu_cu(basis_at_ords, weights, vtSqSelf_d, m0Other_d, vtSqOther_d, normNu, linidx, nuOut_d);
  }
}

__global__ static void
gkyl_spitzer_coll_freq_advance_cu_ker(const struct gkyl_range range,
  const struct gkyl_array* GKYL_RESTRICT basis_at_ords, const struct gkyl_array* GKYL_RESTRICT weights,
  double nufraceps0_fac, double cellav_fac, double r4pieps0_fac, double hbar_fac, double eps0,
  const struct gkyl_array* GKYL_RESTRICT bmag,
  double qSelf, double mSelf, const struct gkyl_array* GKYL_RESTRICT m0Self, const struct gkyl_array* GKYL_RESTRICT vtSqSelf,
  double qOther, double mOther, const struct gkyl_array* GKYL_RESTRICT m0Other, const struct gkyl_array* GKYL_RESTRICT vtSqOther,
  struct gkyl_array* GKYL_RESTRICT nuOut)
{
  double mReduced = 1./(1./mSelf+1./mOther);
  double timeConstFac = nufraceps0_fac*pow(qSelf*qOther,2)/(mSelf*mReduced);

  int idx[3];
  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < range.volume; tid += blockDim.x*gridDim.x) {

    gkyl_sub_range_inv_idx(&range, tid, idx);
    long linidx = gkyl_range_idx(&range, idx);

    const double *bmag_d = gkyl_array_cfetch(bmag, linidx);
    const double *m0Self_d = gkyl_array_cfetch(m0Self, linidx);
    const double *vtSqSelf_d = (const double *) gkyl_array_cfetch(vtSqSelf, linidx);
    const double *m0Other_d = (const double *) gkyl_array_cfetch(m0Other, linidx);
    const double *vtSqOther_d = (const double *) gkyl_array_cfetch(vtSqOther, linidx);

    // Compute the Coulomb logarithm using cell-average values.
    double bmagAv = bmag_d[0]*cellav_fac;
    double m0SelfAv = m0Self_d[0]*cellav_fac;
    double vtSqSelfAv = vtSqSelf_d[0]*cellav_fac;
    double m0OtherAv = m0Other_d[0]*cellav_fac;
    double vtSqOtherAv = vtSqOther_d[0]*cellav_fac;

    double omegaSqSumSelf  = m0SelfAv*pow(qSelf,2)/(eps0*mSelf)+pow(qSelf*bmagAv/mSelf,2);
    double omegaSqSumOther = m0OtherAv*pow(qOther,2)/(eps0*mOther)+pow(qOther*bmagAv/mOther,2);

    double rmaxSumSelf = omegaSqSumSelf/(vtSqSelfAv+3.*vtSqSelfAv)+omegaSqSumOther/(vtSqOtherAv+3.*vtSqSelfAv);
    double rmaxSumOther = omegaSqSumSelf/(vtSqSelfAv+3.*vtSqOtherAv)+omegaSqSumOther/(vtSqOtherAv+3.*vtSqOtherAv);

    double rmaxSelf = 1./sqrt(rmaxSumSelf);
    double rmaxOther = 1./sqrt(rmaxSumOther);

    double uRelSq = 3.*(vtSqOtherAv+vtSqSelfAv);

    double rMin = GKYL_MAX(fabs(qSelf*qOther)*r4pieps0_fac/(mReduced*uRelSq), hbar_fac/(mReduced*sqrt(uRelSq)));

    double logLambda = 0.5*(0.5*log(1.+pow(rmaxSelf/rMin,2))+0.5*log(1.+pow(rmaxOther/rMin,2)));

    // Normalized nu (nu missing density and temperature factors).
    double normNu = timeConstFac*logLambda;

    calc_nu_cu(basis_at_ords, weights, vtSqSelf_d, m0Other_d, vtSqOther_d, normNu, linidx, nuOut_d);
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
    (*range, up->basis_at_ords->on_dev,
     up->weights->on_dev, vtSqSelf->on_dev,
     m0Other->on_dev, vtSqOther->on_dev, normNu, nuOut->on_dev);
}

void
gkyl_spitzer_coll_freq_advance_cu(const gkyl_spitzer_coll_freq *up,
  const struct gkyl_range *range, const struct gkyl_array *bmag,
  double qSelf, double mSelf, const struct gkyl_array *m0Self, const struct gkyl_array *vtSqSelf,
  double qOther, double mOther, const struct gkyl_array *m0Other, const struct gkyl_array *vtSqOther,
  struct gkyl_array *nuOut);
{
  int nblocks = range->nblocks, nthreads = range->nthreads;
  gkyl_spitzer_coll_freq_advance_cu_ker<<<nblocks, nthreads>>>
    (*range, up->basis_at_ords->on_dev, up->weights->on_dev,
     up->nufraceps0_fac, up->cellav_fac, up->r4pieps0_fac, up->hbar_fac, up->eps0,
     bmag->on_dev,
     qSelf, mSelf, m0Self->on_dev, vtSqSelf->on_dev,
     qOther, mOther, m0Other->on_dev, vtSqOther->on_dev, nuOut->on_dev);
}
