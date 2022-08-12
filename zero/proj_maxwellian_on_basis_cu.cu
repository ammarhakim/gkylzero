/* -*- c++ -*- */

extern "C" {
#include <gkyl_proj_maxwellian_on_basis.h>
#include <gkyl_proj_maxwellian_on_basis_priv.h>
#include <gkyl_const.h>
#include <gkyl_range.h>
}

__global__ static void
gkyl_proj_maxwellian_on_basis_lab_mom_cu_ker(int num_quad, const struct gkyl_rect_grid grid,
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

  // double den[tot_conf_quad], vel[tot_conf_quad][vdim], vtsq[tot_conf_quad];
  // MF 2022/08/09: hard-coded to 3x, vdim=3, p=2 for now.
  double den[27], vel[27][3], vtsq[27];

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
      den[n] = 0.0;
      for (int k=0; k<num_conf_basis; ++k) den[n] += M0_d[k]*b_ord[k];

      // drift velocity vector
      for (int d=0; d<vdim; ++d) {
        double M1i_n = 0.0;
        for (int k=0; k<num_conf_basis; ++k) M1i_n += M1i_d[num_conf_basis*d+k]*b_ord[k];
        vel[n][d] = M1i_n/den[n];
      }

      // thermal speed squared
      double M2_n = 0.0;
      for (int k=0; k<num_conf_basis; ++k) M2_n += M2_d[k]*b_ord[k];
      double vsq = 0.0; // velocity^2
      for (int d=0; d<vdim; ++d) vsq += vel[n][d]*vel[n][d];
      vtsq[n] = (M2_n - den[n]*vsq)/(den[n]*vdim);
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
      double nvtsq_q = den[cqidx]/pow(2.0*GKYL_PI*vtsq[cqidx], vdim/2.0);

      comp_to_phys(pdim, (const double *) gkyl_array_cfetch(phase_ordinates, n),
        grid.dx, xc, &xmu[0]);

      double efact = 0.0;
      for (int d=0; d<vdim; ++d)
        efact += pow(xmu[cdim+d]-vel[cqidx][d],2);

      double fmax_o = nvtsq_q*exp(-efact/(2.0*vtsq[cqidx]));

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
    (up->num_quad, up->grid, *phase_r, *conf_r, up->conf_basis_at_ords->on_dev, up->basis_at_ords->on_dev,
     up->ordinates->on_dev, up->weights->on_dev, up->p2c_qidx,
     M0->on_dev, M1i->on_dev, M2->on_dev, fmax->on_dev);
}
