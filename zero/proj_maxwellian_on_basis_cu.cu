/* -*- c++ -*- */

extern "C" {
#include <gkyl_proj_maxwellian_on_basis.h>
#include <gkyl_proj_maxwellian_on_basis_priv.h>
#include <gkyl_const.h>
#include <gkyl_range.h>
}

__global__ static void
gkyl_proj_maxwellian_on_basis_lab_mom_cu_ker(const struct gkyl_rect_grid grid,
  const struct gkyl_range phase_r, const struct gkyl_range conf_r,
  const struct gkyl_array* GKYL_RESTRICT conf_basis_at_ords, 
  const struct gkyl_array* GKYL_RESTRICT phase_basis_at_ords, 
  const struct gkyl_array* GKYL_RESTRICT phase_ordinates, 
  const struct gkyl_array* GKYL_RESTRICT phase_weights, const int *p2c_qidx,
  const struct gkyl_array* GKYL_RESTRICT M0, const struct gkyl_array* GKYL_RESTRICT M1i,
  const struct gkyl_array* GKYL_RESTRICT M2, struct gkyl_array* GKYL_RESTRICT fmax)
{
  int pdim = phase_r.ndim, cdim = conf_r.ndim;
  int vdim = pdim-cdim;

  int num_conf_basis = conf_basis_at_ords->ncomp;
  int num_phase_basis = fmax->ncomp;
  int tot_conf_quad = conf_basis_at_ords->size;
  int tot_phase_quad = phase_basis_at_ords->size;

  // double exp_amp[tot_conf_quad], vel[tot_conf_quad][vdim], vtsq[tot_conf_quad];
  // MF 2022/08/09: hard-coded to 3x, vdim=3, p=2 for now.
  double exp_amp[27], vel[27][3], vtsq[27];

  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM];
  int pidx[GKYL_MAX_DIM], cidx[GKYL_MAX_CDIM];

  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < phase_r.volume; tid += blockDim.x*gridDim.x) {
    gkyl_sub_range_inv_idx(&phase_r, tid, pidx);

    // get conf-space linear index.
    for (unsigned int k = 0; k < conf_r.ndim; k++) cidx[k] = pidx[k];
    long lincC = gkyl_range_idx(&conf_r, cidx);

    const double *M0_d = (const double *) gkyl_array_cfetch(M0, lincC);
    const double *M1i_d = (const double *) gkyl_array_cfetch(M1i, lincC);
    const double *M2_d = (const double *) gkyl_array_cfetch(M2, lincC);

    // compute primitive moments at quadrature nodes
    for (int n=0; n<tot_conf_quad; ++n) {
      const double *b_ord = (const double *) gkyl_array_cfetch(conf_basis_at_ords, n);

      // number density
      double den = 0.0;
      for (int k=0; k<num_conf_basis; ++k) den += M0_d[k]*b_ord[k];

      // drift velocity vector
      for (int d=0; d<vdim; ++d) {
        double M1i_n = 0.0;
        for (int k=0; k<num_conf_basis; ++k) M1i_n += M1i_d[num_conf_basis*d+k]*b_ord[k];
        vel[n][d] = M1i_n/den;
      }

      // thermal speed squared
      double M2_n = 0.0;
      for (int k=0; k<num_conf_basis; ++k) M2_n += M2_d[k]*b_ord[k];
      double vsq = 0.0; // velocity^2
      for (int d=0; d<vdim; ++d) vsq += vel[n][d]*vel[n][d];
      vtsq[n] = (M2_n - den*vsq)/(den*vdim);

      // Amplitude of the exponential.
      exp_amp[n] = den/sqrt(pow(2.0*GKYL_PI*vtsq[n], vdim));
    }

    gkyl_rect_grid_cell_center(&grid, pidx, xc);

    long lidx = gkyl_range_idx(&phase_r, pidx);
    double *fm = (double *) gkyl_array_fetch(fmax, lidx);

    for (int k=0; k<num_phase_basis; ++k) fm[k] = 0.0;

    // compute expansion coefficients of Maxwellian on basis
    // The following is modeled after proj_on_basis in the private header.
    const double *phase_w = (const double *) phase_weights->data;
    const double *phaseb_o = (const double *) phase_basis_at_ords->data;
  
    // compute Maxwellian at phase-space quadrature nodes
    for (int n=0; n<tot_phase_quad; ++n) {

      int cqidx = p2c_qidx[n];

      comp_to_phys(pdim, (const double *) gkyl_array_cfetch(phase_ordinates, n),
        grid.dx, xc, &xmu[0]);

      double efact = 0.0;
      for (int d=0; d<vdim; ++d)
        efact += pow(xmu[cdim+d]-vel[cqidx][d],2);

      double fmax_o = exp_amp[cqidx]*exp(-efact/(2.0*vtsq[cqidx]));

      double tmp = phase_w[n]*fmax_o;
      for (int k=0; k<num_phase_basis; ++k)
        fm[k] += tmp*phaseb_o[k+num_phase_basis*n];
    }

  }
}

__global__ static void
gkyl_proj_maxwellian_on_basis_prim_mom_cu_ker(const struct gkyl_rect_grid grid,
  const struct gkyl_range phase_r, const struct gkyl_range conf_r,
  const struct gkyl_array* GKYL_RESTRICT conf_basis_at_ords, 
  const struct gkyl_array* GKYL_RESTRICT phase_basis_at_ords, 
  const struct gkyl_array* GKYL_RESTRICT phase_ordinates, 
  const struct gkyl_array* GKYL_RESTRICT phase_weights, const int *p2c_qidx,
  const struct gkyl_array* GKYL_RESTRICT m0, const struct gkyl_array* GKYL_RESTRICT udrift,
  const struct gkyl_array* GKYL_RESTRICT vtsq, struct gkyl_array* GKYL_RESTRICT fmax)
{
  int pdim = phase_r.ndim, cdim = conf_r.ndim;
  int vdim = pdim-cdim;

  int num_conf_basis = conf_basis_at_ords->ncomp;
  int num_phase_basis = fmax->ncomp;
  int tot_conf_quad = conf_basis_at_ords->size;
  int tot_phase_quad = phase_basis_at_ords->size;

  // double expamp_o[tot_conf_quad], udrift_o[tot_conf_quad][vdim], vtsq_o[tot_conf_quad];
  // MF 2022/08/09: hard-coded to 3x, vdim=3, p=2 for now.
  double expamp_o[27], udrift_o[27][3], vtsq_o[27];

  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM];
  int pidx[GKYL_MAX_DIM], cidx[GKYL_MAX_CDIM];

  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < phase_r.volume; tid += blockDim.x*gridDim.x) {
    gkyl_sub_range_inv_idx(&phase_r, tid, pidx);

    // get conf-space linear index.
    for (unsigned int k = 0; k < conf_r.ndim; k++) cidx[k] = pidx[k];
    long lincC = gkyl_range_idx(&conf_r, cidx);

    const double *m0_d = (const double *) gkyl_array_cfetch(m0, lincC);
    const double *udrift_d = (const double *) gkyl_array_cfetch(udrift, lincC);
    const double *vtsq_d = (const double *) gkyl_array_cfetch(vtsq, lincC);

    // compute primitive moments at quadrature nodes
    for (int n=0; n<tot_conf_quad; ++n) {
      const double *b_ord = (const double *) gkyl_array_cfetch(conf_basis_at_ords, n);

      double m0_o = 0.0;  // number density
      for (int d=0; d<vdim; ++d) udrift_o[n][d] = 0.0;
      vtsq_o[n] = 0.0;
      for (int k=0; k<num_conf_basis; ++k) {
        m0_o += m0_d[k]*b_ord[k];

        for (int d=0; d<vdim; ++d)
          udrift_o[n][d] += udrift_d[num_conf_basis*d+k]*b_ord[k];

        vtsq_o[n] += vtsq_d[k]*b_ord[k];
      }
      // Amplitude of the exponential.
      expamp_o[n] = m0_o/sqrt(pow(2.0*GKYL_PI*vtsq_o[n], vdim));
    }

    gkyl_rect_grid_cell_center(&grid, pidx, xc);

    long lidx = gkyl_range_idx(&phase_r, pidx);
    double *fm = (double *) gkyl_array_fetch(fmax, lidx);

    for (int k=0; k<num_phase_basis; ++k) fm[k] = 0.0;

    // compute expansion coefficients of Maxwellian on basis
    // The following is modeled after proj_on_basis in the private header.
    const double *phase_w = (const double *) phase_weights->data;
    const double *phaseb_o = (const double *) phase_basis_at_ords->data;
  
    // compute Maxwellian at phase-space quadrature nodes
    for (int n=0; n<tot_phase_quad; ++n) {

      int cqidx = p2c_qidx[n];

      comp_to_phys(pdim, (const double *) gkyl_array_cfetch(phase_ordinates, n),
        grid.dx, xc, &xmu[0]);

      double efact = 0.0;
      for (int d=0; d<vdim; ++d)
        efact += pow(xmu[cdim+d]-udrift_o[cqidx][d],2);

      double fmax_o = expamp_o[cqidx]*exp(-efact/(2.0*vtsq_o[cqidx]));

      double tmp = phase_w[n]*fmax_o;
      for (int k=0; k<num_phase_basis; ++k)
        fm[k] += tmp*phaseb_o[k+num_phase_basis*n];
    }

  }
}

void
gkyl_proj_maxwellian_on_basis_lab_mom_cu(const gkyl_proj_maxwellian_on_basis *up,
  const struct gkyl_range *phase_r, const struct gkyl_range *conf_r,
  const struct gkyl_array *M0, const struct gkyl_array *M1i, const struct gkyl_array *M2,
  struct gkyl_array *fmax)
{
  int nblocks = phase_r->nblocks, nthreads = phase_r->nthreads;
  gkyl_proj_maxwellian_on_basis_lab_mom_cu_ker<<<nblocks, nthreads>>>
    (up->grid, *phase_r, *conf_r, up->conf_basis_at_ords->on_dev, up->basis_at_ords->on_dev,
     up->ordinates->on_dev, up->weights->on_dev, up->p2c_qidx,
     M0->on_dev, M1i->on_dev, M2->on_dev, fmax->on_dev);
}

void
gkyl_proj_maxwellian_on_basis_prim_mom_cu(const gkyl_proj_maxwellian_on_basis *up,
  const struct gkyl_range *phase_r, const struct gkyl_range *conf_r,
  const struct gkyl_array *m0, const struct gkyl_array *udrift, const struct gkyl_array *vtsq,
  struct gkyl_array *fmax)
{
  int nblocks = phase_r->nblocks, nthreads = phase_r->nthreads;
  gkyl_proj_maxwellian_on_basis_prim_mom_cu_ker<<<nblocks, nthreads>>>
    (up->grid, *phase_r, *conf_r, up->conf_basis_at_ords->on_dev, up->basis_at_ords->on_dev,
     up->ordinates->on_dev, up->weights->on_dev, up->p2c_qidx,
     m0->on_dev, udrift->on_dev, vtsq->on_dev, fmax->on_dev);
}

void
gkyl_proj_gkmaxwellian_on_basis_lab_mom_cu(const gkyl_proj_maxwellian_on_basis *up,
  const struct gkyl_range *phase_r, const struct gkyl_range *conf_r,
  const struct gkyl_array *m0, const struct gkyl_array *m1, const struct gkyl_array *m2,
  const struct gkyl_array *bmag, const struct gkyl_array *jacob_tot, double mass,
  struct gkyl_array *fmax)
{
  int nblocks = phase_r->nblocks, nthreads = phase_r->nthreads;
//  gkyl_proj_gkmaxwellian_on_basis_lab_mom_cu_ker<<<nblocks, nthreads>>>
//    (up->grid, *phase_r, *conf_r, up->conf_basis_at_ords->on_dev, up->basis_at_ords->on_dev,
//     up->ordinates->on_dev, up->weights->on_dev, up->p2c_qidx,
//     m0->on_dev, m1->on_dev, m2->on_dev, bmag->on_dev, jacob_tot->on_dev, mass, fmax->on_dev);
}

void
gkyl_proj_gkmaxwellian_on_basis_prim_mom_cu(const gkyl_proj_maxwellian_on_basis *up,
  const struct gkyl_range *phase_r, const struct gkyl_range *conf_r,
  const struct gkyl_array *m0, const struct gkyl_array *upar, const struct gkyl_array *vtsq,
  const struct gkyl_array *bmag, const struct gkyl_array *jacob_tot, double mass,
  struct gkyl_array *fmax)
{
  int nblocks = phase_r->nblocks, nthreads = phase_r->nthreads;
//  gkyl_proj_gkmaxwellian_on_basis_prim_mom_cu_ker<<<nblocks, nthreads>>>
//    (up->grid, *phase_r, *conf_r, up->conf_basis_at_ords->on_dev, up->basis_at_ords->on_dev,
//     up->ordinates->on_dev, up->weights->on_dev, up->p2c_qidx,
//     m0->on_dev, upar->on_dev, vtsq->on_dev, bmag->on_dev, jacob_tot->on_dev, mass, fmax->on_dev);
}
