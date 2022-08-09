/* -*- c++ -*- */

extern "C" {
#include <gkyl_proj_maxwellian_on_basis.h>
#include <gkyl_proj_maxwellian_on_basis_priv.h>
}

__global__ static void
gkyl_proj_maxwellian_on_basis_lab_mom_cu_ker(int num_quad, const struct gkyl_rect_grid *grid,
  const struct gkyl_range phase_r, const struct gkyl_range conf_r,
  const struct gkyl_array* GKYL_RESTRICT conf_basis_at_ords, 
  const struct gkyl_array* GKYL_RESTRICT phase_basis_at_ords, 
  const struct gkyl_array* GKYL_RESTRICT phase_ordinates, 
  const struct gkyl_array* GKYL_RESTRICT phase_weights, 
  const struct gkyl_array* GKYL_RESTRICT fun_at_ords, 
  const struct gkyl_array* GKYL_RESTRICT M0, struct gkyl_array* GKYL_RESTRICT M1i
  const struct gkyl_array* GKYL_RESTRICT M2, const gkyl_mom_calc* fmax)
{
  int pdim = phase_r.ndim, cdim = conf_r.ndim;
  int vdim = pdim-cdim;

  int num_conf_basis = conf_basis_at_ords->ncomp;
  int num_phase_basis = fmax->ncomp;
  int tot_conf_quad = conf_basis_at_ords->size;
  int tot_phase_quad = phase_basis_at_ords->size;

  double num[tot_conf_quad], vel[tot_conf_quad][vdim], vtsq[tot_conf_quad];

  // create range to loop over config-space and phase-space quadrature points
  int qshape[GKYL_MAX_DIM];
  for (int i=0; i<cdim; ++i) qshape[i] = num_quad;
  struct gkyl_range conf_qrange;
  gkyl_range_init_from_shape(&conf_qrange, cdim, qshape);

  for (int i=0; i<pdim; ++i) qshape[i] = num_quad;
  struct gkyl_range phase_qrange;
  gkyl_range_init_from_shape(&phase_qrange, pdim, qshape);

  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM];
  int pidx[GKYL_MAX_DIM], cidx[GKYL_MAX_CDIM];

  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < phase_r.volume; tid += blockDim.x*gridDim.x) {
    gkyl_sub_range_inv_idx(&phase_r, tid, pidx);

    // get conf-space linear index.
    for (unsigned int k = 0; k < conf_r.ndim; k++) cidx[k] = pidx[k];
    long lincC = gkyl_range_idx(&conf_r, cidx);

    const double *M0_d = gkyl_array_cfetch(M0, lincC);
    const double *M1i_d = gkyl_array_cfetch(M1i, lincC);
    const double *M2_d = gkyl_array_cfetch(M2, lincC);

    // compute primitive moments at quadrature nodes
    for (int n=0; n<tot_conf_quad; ++n) {
      const double *b_ord = gkyl_array_cfetch(conf_basis_at_ords, n);

      // number density
      num[n] = 0.0;
      for (int k=0; k<num_conf_basis; ++k) num[n] += M0_d[k]*b_ord[k];

      // drift velocity vector
      for (int d=0; d<vdim; ++d) {
        double M1i_n = 0.0;
        for (int k=0; k<num_conf_basis; ++k)
          M1i_n += M1i_d[num_conf_basis*d+k]*b_ord[k];
        vel[n][d] = M1i_n/num[n];
      }

      // thermal speed squared
      double M2_n = 0.0;
      for (int k=0; k<num_conf_basis; ++k)
        M2_n += M2_d[k]*b_ord[k];

      double vsq = 0.0; // velocity^2
      for (int d=0; d<vdim; ++d) vsq += vel[n][d]*vel[n][d];
      vtsq[n] = (M2_n - num[n]*vsq)/(num[n]*vdim);
    }

    gkyl_rect_grid_cell_center(grid, pidx, xc);

    // compute Maxwellian at phase-space quadrature nodes
    struct gkyl_range_iter qiter;
    gkyl_range_iter_init(&qiter, &phase_qrange);
    while (gkyl_range_iter_next(&qiter)) {

      long cqidx = gkyl_range_idx(&conf_qrange, qiter.idx);
      double nvtsq_q = num[cqidx]/pow(2*GKYL_PI*vtsq[cqidx], vdim/2.0);

      long pqidx = gkyl_range_idx(&phase_qrange, qiter.idx);

      comp_to_phys(pdim, gkyl_array_cfetch(phase_ordinates, pqidx),
        up->grid.dx, xc, xmu);

      double efact = 0.0;
      for (int d=0; d<vdim; ++d)
        efact += (vel[cqidx][d]-xmu[cdim+d])*(vel[cqidx][d]-xmu[cdim+d]);

      double *fq = gkyl_array_fetch(fun_at_ords, pqidx);
      fq[0] = nvtsq_q*exp(-efact/(2*vtsq[cqidx]));
    }

    // compute expansion coefficients of Maxwellian on basis
    // The following is modeled after proj_on_basis in the private header.
    long lidx = gkyl_range_idx(&phase_r, pidx);
    double *fm = gkyl_array_fetch(fmax, lidx);
  
    for (int k=0; k<num_phase_basis; ++k) fm[k] = 0.0;
  
    for (int imu=0; imu<tot_phase_quad; ++imu) {
      double tmp = phase_weights->data[imu]*fun_at_ords->data[imu];
      for (int k=0; k<num_phase_basis; ++k)
        fm[k] += tmp*phase_basis_at_ords->data[k+num_phase_basis*imu];
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
  gkyl_proj_maxwellian_on_basis_advance_cu_ker<<<nblocks, nthreads>>>
    (up->num_quad, &up->grid, *phase_r, *conf_r, up->conf_basis_at_ords, up->basis_at_ords, up->ordinates, 
     up->weights, up->fun_at_ords, M0->on_dev, M1i->on_dev, M2->on_dev, fmax->on_dev);
}
