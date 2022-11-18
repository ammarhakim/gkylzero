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

__global__ static void
gkyl_proj_gkmaxwellian_on_basis_lab_mom_cu_ker(const struct gkyl_rect_grid grid,
  const struct gkyl_range phase_r, const struct gkyl_range conf_r,
  const struct gkyl_array* GKYL_RESTRICT conf_basis_at_ords, 
  const struct gkyl_array* GKYL_RESTRICT phase_basis_at_ords, 
  const struct gkyl_array* GKYL_RESTRICT phase_ordinates, 
  const struct gkyl_array* GKYL_RESTRICT phase_weights, const int *p2c_qidx,
  const struct gkyl_array* GKYL_RESTRICT m0, const struct gkyl_array* GKYL_RESTRICT m1,
  const struct gkyl_array* GKYL_RESTRICT m2, const struct gkyl_array* GKYL_RESTRICT bmag,
  const struct gkyl_array* GKYL_RESTRICT jacob_tot, double mass,
  struct gkyl_array* GKYL_RESTRICT fmax)
{
  double fJacB_floor = 1.e-40;
  int pdim = phase_r.ndim, cdim = conf_r.ndim;
  int vdim = pdim-cdim;
  int vdim_phys = vdim==1 ? 1 : 3;

  int num_conf_basis = conf_basis_at_ords->ncomp;
  int num_phase_basis = fmax->ncomp;
  int tot_conf_quad = conf_basis_at_ords->size;
  int tot_phase_quad = phase_basis_at_ords->size;

  // double exp_amp[tot_conf_quad], vel[tot_conf_quad][vdim], vtsq[tot_conf_quad];
  // MF 2022/08/09: hard-coded to 3x, vdim=3, p=2 for now.
  double exp_amp[27], upar[27], vtsq[27], bfield[27];

  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM] = {0.};
  int pidx[GKYL_MAX_DIM], cidx[GKYL_MAX_CDIM];

  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < phase_r.volume; tid += blockDim.x*gridDim.x) {
    gkyl_sub_range_inv_idx(&phase_r, tid, pidx);

    // get conf-space linear index.
    for (unsigned int k = 0; k < conf_r.ndim; k++) cidx[k] = pidx[k];
    long lincC = gkyl_range_idx(&conf_r, cidx);

    const double *m0_d = (const double *) gkyl_array_cfetch(m0, lincC);
    const double *m1_d = (const double *) gkyl_array_cfetch(m1, lincC);
    const double *m2_d = (const double *) gkyl_array_cfetch(m2, lincC);
    const double *bmag_d = (const double *) gkyl_array_cfetch(bmag, lincC);
    const double *jactot_d = (const double *) gkyl_array_cfetch(jacob_tot, lincC);

    // compute primitive moments at quadrature nodes
    for (int n=0; n<tot_conf_quad; ++n) {
      const double *b_ord = (const double *) gkyl_array_cfetch(conf_basis_at_ords, n);
    
      double den = 0.;   // number density
      double m1_n = 0.;  // momentum density.
      double m2_n = 0.;  // kinetic energy density.
      double jac_n = 0.; // total jacobian (conf * guiding center jacobian).
      bfield[n] = 0.; // magnetic field amplitude.
      for (int k=0; k<num_conf_basis; ++k) {
        den += m0_d[k]*b_ord[k];
        m1_n += m1_d[k]*b_ord[k];
        m2_n += m2_d[k]*b_ord[k];
        jac_n += jactot_d[k]*b_ord[k];
        bfield[n] += bmag_d[k]*b_ord[k];
      }
      // parallel drift speed.
      upar[n] = m1_n/den;
      // thermal speed squared.
      vtsq[n] = (m2_n - den*upar[n]*upar[n])/(den*vdim_phys);

      // Amplitude of the exponential.
      if ((den > 0.) && (vtsq[n]>0.))
        exp_amp[n] = jac_n*den/sqrt(pow(2.0*GKYL_PI*vtsq[n], vdim_phys));
      else
        exp_amp[n] = 0.;
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
      // vpar term.
      efact += pow(xmu[cdim]-upar[cqidx],2);
      // mu term (only for 2v, vdim_phys=3).
      efact += (vdim_phys-1)*xmu[cdim+1]*bfield[cqidx]/mass;

      double fmax_o = (fJacB_floor+exp_amp[cqidx])*exp(-efact/(2.0*vtsq[cqidx]));

      double tmp = phase_w[n]*fmax_o;
      for (int k=0; k<num_phase_basis; ++k)
        fm[k] += tmp*phaseb_o[k+num_phase_basis*n];
    }
  }
}

__global__ static void
gkyl_proj_gkmaxwellian_on_basis_prim_mom_cu_ker(const struct gkyl_rect_grid grid,
  const struct gkyl_range phase_r, const struct gkyl_range conf_r,
  const struct gkyl_array* GKYL_RESTRICT conf_basis_at_ords, 
  const struct gkyl_array* GKYL_RESTRICT phase_basis_at_ords, 
  const struct gkyl_array* GKYL_RESTRICT phase_ordinates, 
  const struct gkyl_array* GKYL_RESTRICT phase_weights, const int *p2c_qidx,
  const struct gkyl_array* GKYL_RESTRICT m0, const struct gkyl_array* GKYL_RESTRICT upar,
  const struct gkyl_array* GKYL_RESTRICT vtsq, const struct gkyl_array* GKYL_RESTRICT bmag,
  const struct gkyl_array* GKYL_RESTRICT jacob_tot, double mass,
  struct gkyl_array* GKYL_RESTRICT fmax)
{
  double fJacB_floor = 1.e-40;
  int pdim = phase_r.ndim, cdim = conf_r.ndim;
  int vdim = pdim-cdim;
  int vdim_phys = vdim==1 ? 1 : 3;

  int num_conf_basis = conf_basis_at_ords->ncomp;
  int num_phase_basis = fmax->ncomp;
  int tot_conf_quad = conf_basis_at_ords->size;
  int tot_phase_quad = phase_basis_at_ords->size;

  // double expamp_o[tot_conf_quad], upar_o[tot_conf_quad][vdim], vtsq_o[tot_conf_quad];
  // MF 2022/08/09: hard-coded to 3x, vdim=3, p=2 for now.
  double expamp_o[27], upar_o[27], vtsq_o[27], bmag_o[27];

  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM] = {0.};
  int pidx[GKYL_MAX_DIM], cidx[GKYL_MAX_CDIM];

  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < phase_r.volume; tid += blockDim.x*gridDim.x) {
    gkyl_sub_range_inv_idx(&phase_r, tid, pidx);

    // get conf-space linear index.
    for (unsigned int k = 0; k < conf_r.ndim; k++) cidx[k] = pidx[k];
    long lincC = gkyl_range_idx(&conf_r, cidx);

    const double *m0_d = (const double *) gkyl_array_cfetch(m0, lincC);
    const double *upar_d = (const double *) gkyl_array_cfetch(upar, lincC);
    const double *vtsq_d = (const double *) gkyl_array_cfetch(vtsq, lincC);
    const double *bmag_d = (const double *) gkyl_array_cfetch(bmag, lincC);
    const double *jactot_d = (const double *) gkyl_array_cfetch(jacob_tot, lincC);

    // compute primitive moments at quadrature nodes
    for (int n=0; n<tot_conf_quad; ++n) {
      const double *b_ord = (const double *) gkyl_array_cfetch(conf_basis_at_ords, n);

      double m0_o = 0.0;  // number density
      double jac_o = 0.0;  // total jacobian (conf * guiding center jacobian).
      upar_o[n] = 0.0;
      vtsq_o[n] = 0.0;
      bmag_o[n] = 0.0;
      for (int k=0; k<num_conf_basis; ++k) {
        m0_o += m0_d[k]*b_ord[k];
        jac_o += jactot_d[k]*b_ord[k];
        bmag_o[n] += bmag_d[k]*b_ord[k];
        upar_o[n] += upar_d[k]*b_ord[k];
        vtsq_o[n] += vtsq_d[k]*b_ord[k];
      }
      // Amplitude of the exponential.
      if ((m0_o > 0.) && (vtsq_o[n]>0.))
        expamp_o[n] = jac_o*m0_o/sqrt(pow(2.0*GKYL_PI*vtsq_o[n], vdim_phys));
      else
        expamp_o[n] = 0.;
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
      // vpar term.
      efact += pow(xmu[cdim]-upar_o[cqidx],2);
      // mu term (only for 2v, vdim_phys=3).
      efact += (vdim_phys-1)*xmu[cdim+1]*bmag_o[cqidx]/mass;

      double fmax_o = (fJacB_floor+expamp_o[cqidx])*exp(-efact/(2.0*vtsq_o[cqidx]));

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
  gkyl_proj_gkmaxwellian_on_basis_lab_mom_cu_ker<<<nblocks, nthreads>>>
    (up->grid, *phase_r, *conf_r, up->conf_basis_at_ords->on_dev, up->basis_at_ords->on_dev,
     up->ordinates->on_dev, up->weights->on_dev, up->p2c_qidx,
     m0->on_dev, m1->on_dev, m2->on_dev, bmag->on_dev, jacob_tot->on_dev, mass, fmax->on_dev);
}

void
gkyl_proj_gkmaxwellian_on_basis_prim_mom_cu(const gkyl_proj_maxwellian_on_basis *up,
  const struct gkyl_range *phase_r, const struct gkyl_range *conf_r,
  const struct gkyl_array *m0, const struct gkyl_array *upar, const struct gkyl_array *vtsq,
  const struct gkyl_array *bmag, const struct gkyl_array *jacob_tot, double mass,
  struct gkyl_array *fmax)
{
  int nblocks = phase_r->nblocks, nthreads = phase_r->nthreads;
  gkyl_proj_gkmaxwellian_on_basis_prim_mom_cu_ker<<<nblocks, nthreads>>>
    (up->grid, *phase_r, *conf_r, up->conf_basis_at_ords->on_dev, up->basis_at_ords->on_dev,
     up->ordinates->on_dev, up->weights->on_dev, up->p2c_qidx,
     m0->on_dev, upar->on_dev, vtsq->on_dev, bmag->on_dev, jacob_tot->on_dev, mass, fmax->on_dev);
}
