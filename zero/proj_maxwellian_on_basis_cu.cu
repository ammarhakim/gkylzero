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
  const struct gkyl_array* GKYL_RESTRICT moms, struct gkyl_array* GKYL_RESTRICT fmax)
{
  double f_floor = 1.e-40;
  int pdim = phase_r.ndim, cdim = conf_r.ndim;
  int vdim = pdim-cdim;

  int num_conf_basis = conf_basis_at_ords->ncomp;
  int num_phase_basis = fmax->ncomp;
  int tot_conf_quad = conf_basis_at_ords->size;
  int tot_phase_quad = phase_basis_at_ords->size;

  // double exp_amp[tot_conf_quad], udrift[tot_conf_quad][vdim], vtsq[tot_conf_quad];
  // MF 2022/08/09: hard-coded to 3x, vdim=3, p=2 for now.
  double exp_amp[27], udrift[27][3], vtsq[27];

  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM];
  int pidx[GKYL_MAX_DIM], cidx[GKYL_MAX_CDIM];

  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < phase_r.volume; tid += blockDim.x*gridDim.x) {
    gkyl_sub_range_inv_idx(&phase_r, tid, pidx);

    // get conf-space linear index.
    for (unsigned int k = 0; k < conf_r.ndim; k++) cidx[k] = pidx[k];
    long lincC = gkyl_range_idx(&conf_r, cidx);

    const double *moms_d = (const double *) gkyl_array_cfetch(moms, lincC);
    const double *m0_d  = moms_d;
    const double *m1i_d = &moms_d[num_conf_basis];
    const double *m2_d  = &moms_d[(vdim+1)*num_conf_basis];

    // compute primitive moments at quadrature nodes
    for (int n=0; n<tot_conf_quad; ++n) {
      const double *b_ord = (const double *) gkyl_array_cfetch(conf_basis_at_ords, n);

      // number density
      double den = 0.0;
      for (int k=0; k<num_conf_basis; ++k) den += m0_d[k]*b_ord[k];

      // drift velocity vector
      for (int d=0; d<vdim; ++d) {
        double m1i_n = 0.0;
        for (int k=0; k<num_conf_basis; ++k) m1i_n += m1i_d[num_conf_basis*d+k]*b_ord[k];
        udrift[n][d] = m1i_n/den;
      }

      // thermal speed squared
      double m2_n = 0.0;
      for (int k=0; k<num_conf_basis; ++k) m2_n += m2_d[k]*b_ord[k];
      double usq = 0.0; // velocity^2
      for (int d=0; d<vdim; ++d) usq += udrift[n][d]*udrift[n][d];
      vtsq[n] = (m2_n - den*usq)/(den*vdim);

      // Amplitude of the exponential.
      if ((den > 0.) && (vtsq[n]>0.))
         exp_amp[n] = den/sqrt(pow(2.0*GKYL_PI*vtsq[n], vdim));
      else
         exp_amp[n] = 0.0;
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
        efact += pow(xmu[cdim+d]-udrift[cqidx][d],2);

      double fmax_o = vtsq[cqidx] > 0.0 ? f_floor + exp_amp[cqidx]*exp(-efact/(2.0*vtsq[cqidx])) : f_floor;

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
  const struct gkyl_array* GKYL_RESTRICT moms, const struct gkyl_array* GKYL_RESTRICT prim_moms,
  struct gkyl_array* GKYL_RESTRICT fmax)
{
  double f_floor = 1.e-40;
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

    const double *m0_d = (const double *) gkyl_array_cfetch(moms, lincC);
    const double *prim_moms_d = (const double *) gkyl_array_cfetch(prim_moms, lincC);
    const double *udrift_d = prim_moms_d;
    const double *vtsq_d   = &prim_moms_d[vdim*num_conf_basis];

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
      if ((m0_o > 0.) && (vtsq_o[n]>0.))
         expamp_o[n] = m0_o/sqrt(pow(2.0*GKYL_PI*vtsq_o[n], vdim));
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
      for (int d=0; d<vdim; ++d)
        efact += pow(xmu[cdim+d]-udrift_o[cqidx][d],2);

      double fmax_o = vtsq_o[cqidx] > 0.0 ? f_floor + expamp_o[cqidx]*exp(-efact/(2.0*vtsq_o[cqidx])) : f_floor;

      double tmp = phase_w[n]*fmax_o;
      for (int k=0; k<num_phase_basis; ++k)
        fm[k] += tmp*phaseb_o[k+num_phase_basis*n];
    }

  }
}

__global__ static void
gkyl_proj_gkmaxwellian_on_basis_lab_mom_cu_ker(const struct gkyl_rect_grid grid,
  const struct gkyl_range phase_r, const struct gkyl_range conf_r, struct gkyl_range vel_r,
  const struct gkyl_array* GKYL_RESTRICT conf_basis_at_ords, 
  const struct gkyl_array* GKYL_RESTRICT phase_basis_at_ords, 
  const struct gkyl_array* GKYL_RESTRICT phase_ordinates, 
  const struct gkyl_array* GKYL_RESTRICT phase_weights, const int *p2c_qidx,
  const struct gkyl_array* GKYL_RESTRICT moms, const struct gkyl_array* GKYL_RESTRICT bmag,
  const struct gkyl_array* GKYL_RESTRICT jacob_tot,
  struct gkyl_array* GKYL_RESTRICT vmap, struct gkyl_basis* GKYL_RESTRICT vmap_basis,
  double mass, struct gkyl_array* GKYL_RESTRICT fmax)
{
  double fJacB_floor = 1.e-40;
  int pdim = phase_r.ndim, cdim = conf_r.ndim;
  int vdim = pdim-cdim;
  int vdim_phys = vdim==1 ? 1 : 3;

  int num_conf_basis = conf_basis_at_ords->ncomp;
  int num_phase_basis = fmax->ncomp;
  int tot_conf_quad = conf_basis_at_ords->size;
  int tot_phase_quad = phase_basis_at_ords->size;

  // double exp_amp[tot_conf_quad], udrift[tot_conf_quad][vdim], vtsq[tot_conf_quad];
  // MF 2022/08/09: hard-coded to 3x, vdim=3, p=2 for now.
  double exp_amp[27], upar[27], vtsq[27], bfield[27];

  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM] = {0.};
  int pidx[GKYL_MAX_DIM], cidx[GKYL_MAX_CDIM], vidx[2];

  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < phase_r.volume; tid += blockDim.x*gridDim.x) {
    gkyl_sub_range_inv_idx(&phase_r, tid, pidx);

    // get conf-space linear index.
    for (unsigned int k = 0; k < conf_r.ndim; k++) cidx[k] = pidx[k];
    long lincC = gkyl_range_idx(&conf_r, cidx);

    const double *moms_d = (const double *) gkyl_array_cfetch(moms, lincC);
    const double *m0_d = moms_d;
    const double *m1_d = &moms_d[num_conf_basis];
    const double *m2_d = &moms_d[2*num_conf_basis];
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
  
    for (unsigned int d=cdim; d<pdim; d++) vidx[d-cdim] = pidx[d];
    long vlinidx = gkyl_range_idx(&vel_r, vidx);
    const double *vmap_d = (const double *) gkyl_array_cfetch(vmap, vlinidx);

    // compute Maxwellian at phase-space quadrature nodes
    for (int n=0; n<tot_phase_quad; ++n) {

      int cqidx = p2c_qidx[n];
      const double *xcomp_d = (const double *) gkyl_array_cfetch(phase_ordinates, n);

      // Convert comp velocity coordinate to phys velocity coord.
      double xcomp[1];
      for (int vd=0; vd<vdim; vd++) {
        xcomp[0] = xcomp_d[cdim+vd];
        xmu[cdim+vd] = vmap_basis->eval_expand(xcomp, vmap_d+vd*vmap_basis->num_basis);
      }

      double efact = 0.0;
      // vpar term.
      efact += pow(xmu[cdim]-upar[cqidx],2);
      // mu term (only for 2v, vdim_phys=3).
      efact += (vdim_phys-1)*xmu[cdim+1]*bfield[cqidx]/mass;

      double fmax_o = vtsq[cqidx] > 0.0 ?  fJacB_floor+exp_amp[cqidx]*exp(-efact/(2.0*vtsq[cqidx])): fJacB_floor;

      double tmp = phase_w[n]*fmax_o;
      for (int k=0; k<num_phase_basis; ++k)
        fm[k] += tmp*phaseb_o[k+num_phase_basis*n];
    }
  }
}

__global__ static void
gkyl_proj_gkmaxwellian_on_basis_prim_mom_cu_ker(const struct gkyl_rect_grid grid,
  const struct gkyl_range phase_r, const struct gkyl_range conf_r, struct gkyl_range vel_r,
  const struct gkyl_array* GKYL_RESTRICT conf_basis_at_ords, 
  const struct gkyl_array* GKYL_RESTRICT phase_basis_at_ords, 
  const struct gkyl_array* GKYL_RESTRICT phase_ordinates, 
  const struct gkyl_array* GKYL_RESTRICT phase_weights, const int *p2c_qidx,
  const struct gkyl_array* GKYL_RESTRICT prim_moms,
  const struct gkyl_array* GKYL_RESTRICT bmag, const struct gkyl_array* GKYL_RESTRICT jacob_tot,
  struct gkyl_array* GKYL_RESTRICT vmap, struct gkyl_basis* GKYL_RESTRICT vmap_basis,
  double mass, struct gkyl_array* GKYL_RESTRICT fmax)
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
  int pidx[GKYL_MAX_DIM], cidx[GKYL_MAX_CDIM], vidx[2];

  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < phase_r.volume; tid += blockDim.x*gridDim.x) {
    gkyl_sub_range_inv_idx(&phase_r, tid, pidx);

    // get conf-space linear index.
    for (unsigned int k = 0; k < cdim; k++) cidx[k] = pidx[k];
    long lincC = gkyl_range_idx(&conf_r, cidx);

    const double *prim_moms_d = (const double *) gkyl_array_cfetch(prim_moms, lincC);
    const double *m0_d = prim_moms_d;
    const double *upar_d = &prim_moms_d[num_conf_basis];
    const double *vtsq_d = &prim_moms_d[2*num_conf_basis];
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
  
    for (unsigned int d=cdim; d<pdim; d++) vidx[d-cdim] = pidx[d];
    long vlinidx = gkyl_range_idx(&vel_r, vidx);
    const double *vmap_d = (const double *) gkyl_array_cfetch(vmap, vlinidx);

    // compute Maxwellian at phase-space quadrature nodes
    for (int n=0; n<tot_phase_quad; ++n) {

      int cqidx = p2c_qidx[n];
      const double *xcomp_d = (const double *) gkyl_array_cfetch(phase_ordinates, n);

      // Convert comp velocity coordinate to phys velocity coord.
      double xcomp[1];
      for (int vd=0; vd<vdim; vd++) {
        xcomp[0] = xcomp_d[cdim+vd];
        xmu[cdim+vd] = vmap_basis->eval_expand(xcomp, vmap_d+vd*vmap_basis->num_basis);
      }

      double efact = 0.0;
      // vpar term.
      efact += pow(xmu[cdim]-upar_o[cqidx],2);
      // mu term (only for 2v, vdim_phys=3).
      efact += (vdim_phys-1)*xmu[cdim+1]*bmag_o[cqidx]/mass;

      double fmax_o = vtsq_o[cqidx] > 0.0 ? fJacB_floor+expamp_o[cqidx]*exp(-efact/(2.0*vtsq_o[cqidx])) : fJacB_floor;

      double tmp = phase_w[n]*fmax_o;
      for (int k=0; k<num_phase_basis; ++k)
        fm[k] += tmp*phaseb_o[k+num_phase_basis*n];
    }
  }
}

void
gkyl_proj_maxwellian_on_basis_lab_mom_cu(const gkyl_proj_maxwellian_on_basis *up,
  const struct gkyl_range *phase_r, const struct gkyl_range *conf_r,
  const struct gkyl_array *moms, struct gkyl_array *fmax)
{
  int nblocks = phase_r->nblocks, nthreads = phase_r->nthreads;
  gkyl_proj_maxwellian_on_basis_lab_mom_cu_ker<<<nblocks, nthreads>>>
    (up->grid, *phase_r, *conf_r, up->conf_basis_at_ords->on_dev, up->basis_at_ords->on_dev,
     up->ordinates->on_dev, up->weights->on_dev, up->p2c_qidx,
     moms->on_dev, fmax->on_dev);
}

void
gkyl_proj_maxwellian_on_basis_prim_mom_cu(const gkyl_proj_maxwellian_on_basis *up,
  const struct gkyl_range *phase_r, const struct gkyl_range *conf_r,
  const struct gkyl_array *moms, const struct gkyl_array *prim_moms,
  struct gkyl_array *fmax)
{
  int nblocks = phase_r->nblocks, nthreads = phase_r->nthreads;
  gkyl_proj_maxwellian_on_basis_prim_mom_cu_ker<<<nblocks, nthreads>>>
    (up->grid, *phase_r, *conf_r, up->conf_basis_at_ords->on_dev, up->basis_at_ords->on_dev,
     up->ordinates->on_dev, up->weights->on_dev, up->p2c_qidx,
     moms->on_dev, prim_moms->on_dev, fmax->on_dev);
}

void
gkyl_proj_gkmaxwellian_on_basis_lab_mom_cu(const gkyl_proj_maxwellian_on_basis *up,
  const struct gkyl_range *phase_r, const struct gkyl_range *conf_r,
  const struct gkyl_array *moms, const struct gkyl_array *bmag,
  const struct gkyl_array *jacob_tot, double mass, struct gkyl_array *fmax)
{
  const struct gkyl_velocity_map *gvm = up->vel_map;

  int nblocks = phase_r->nblocks, nthreads = phase_r->nthreads;
  gkyl_proj_gkmaxwellian_on_basis_lab_mom_cu_ker<<<nblocks, nthreads>>>
    (up->grid, *phase_r, *conf_r, gvm->local_ext_vel, up->conf_basis_at_ords->on_dev,
     up->basis_at_ords->on_dev, up->ordinates->on_dev, up->weights->on_dev, up->p2c_qidx,
     moms->on_dev, bmag->on_dev, jacob_tot->on_dev,
     gvm->vmap->on_dev, gvm->vmap_basis, mass, fmax->on_dev);
}

void
gkyl_proj_gkmaxwellian_on_basis_prim_mom_cu(const gkyl_proj_maxwellian_on_basis *up,
  const struct gkyl_range *phase_r, const struct gkyl_range *conf_r,
  const struct gkyl_array *prim_moms,
  const struct gkyl_array *bmag, const struct gkyl_array *jacob_tot, double mass,
  struct gkyl_array *fmax)
{
  const struct gkyl_velocity_map *gvm = up->vel_map;

  int nblocks = phase_r->nblocks, nthreads = phase_r->nthreads;
  gkyl_proj_gkmaxwellian_on_basis_prim_mom_cu_ker<<<nblocks, nthreads>>>
    (up->grid, *phase_r, *conf_r, gvm->local_ext_vel, up->conf_basis_at_ords->on_dev,
     up->basis_at_ords->on_dev, up->ordinates->on_dev, up->weights->on_dev, up->p2c_qidx,
     prim_moms->on_dev, bmag->on_dev, jacob_tot->on_dev,
     gvm->vmap->on_dev, gvm->vmap_basis, mass, fmax->on_dev);
}
