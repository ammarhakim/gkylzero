#include "gkyl_const.h"
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_proj_maxwellian_pots_on_basis.h>
#include <gkyl_proj_maxwellian_pots_on_basis_priv.h>

// Create range to loop over quadrature points
static inline struct gkyl_range get_qrange(int cdim, int dim, int num_quad) {
  int qshape[GKYL_MAX_DIM];
  for (int i=0; i<dim; ++i) qshape[i] = num_quad;
  for (int i=cdim; i<dim; ++i) qshape[i] = num_quad;
  struct gkyl_range qrange;
  gkyl_range_init_from_shape(&qrange, dim, qshape);
  return qrange;
}

static int init_quad_values(int cdim, const struct gkyl_basis *basis, int num_quad, struct gkyl_array **ordinates,
  struct gkyl_array **weights, struct gkyl_array **basis_at_ords)
{
  int ndim = basis->ndim;
  int vdim = ndim-cdim;

  double ordinates1[num_quad], weights1[num_quad];

  // Use pre-computed values of ords and weights if possible
  if (num_quad <= gkyl_gauss_max) {
    memcpy(ordinates1, gkyl_gauss_ordinates[num_quad], sizeof(double[num_quad]));
    memcpy(weights1, gkyl_gauss_weights[num_quad], sizeof(double[num_quad]));
  } else {
    gkyl_gauleg(-1, 1, ordinates1, weights1, num_quad);
  }

  struct gkyl_range qrange = get_qrange(cdim, ndim, num_quad);

  int tot_quad = qrange.volume;

  // create ordinates and weights for multi-dimensional quadrature
  struct gkyl_array *ordinates_ho = gkyl_array_new(GKYL_DOUBLE, ndim, tot_quad);
  struct gkyl_array *weights_ho = gkyl_array_new(GKYL_DOUBLE, 1, tot_quad);

  // Not supporting CUDA yet, will need to check for use_gpu and use correct array initializer
  *ordinates = gkyl_array_new(GKYL_DOUBLE, ndim, tot_quad);
  *weights = gkyl_array_new(GKYL_DOUBLE, 1, tot_quad);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &qrange);

  while (gkyl_range_iter_next(&iter)) {
    int node = gkyl_range_idx(&qrange, iter.idx);

    // set ordinates
    double *ord = gkyl_array_fetch(ordinates_ho, node);
    
    for (int i=0; i<ndim; ++i)
      ord[i] = ordinates1[iter.idx[i]-qrange.lower[i]];

    // set weights
    double *wgt = gkyl_array_fetch(weights_ho, node);
    wgt[0] = 1.0;
    for (int i=0; i<ndim; ++i)
      wgt[0] *= weights1[iter.idx[i]-qrange.lower[i]];
  }

  // pre-compute basis functions at ordinates
  struct gkyl_array *basis_at_ords_ho = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, tot_quad);
  // CUDA support to be added
  *basis_at_ords = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, tot_quad);
  for (int n=0; n<tot_quad; ++n)
    basis->eval(gkyl_array_fetch(ordinates_ho, n), gkyl_array_fetch(basis_at_ords_ho, n));

  // copy host array to device array
  gkyl_array_copy(*ordinates, ordinates_ho);
  gkyl_array_copy(*weights, weights_ho);
  gkyl_array_copy(*basis_at_ords, basis_at_ords_ho);

  gkyl_array_release(ordinates_ho);
  gkyl_array_release(weights_ho);
  gkyl_array_release(basis_at_ords_ho);

  return tot_quad;
}

static void
proj_on_basis(const gkyl_proj_maxwellian_pots_on_basis *up, const struct gkyl_array *fun_at_ords, double* f)
{
  int num_basis = up->num_phase_basis;
  int tot_quad = up->tot_quad;

  const double* GKYL_RESTRICT weights = up->weights->data;
  const double* GKYL_RESTRICT basis_at_ords = up->basis_at_ords->data;
  const double* GKYL_RESTRICT func_at_ords = fun_at_ords->data;

  for (int k=0; k<num_basis; ++k) f[k] = 0.0;
  
  for (int imu=0; imu<tot_quad; ++imu) {
    double tmp = weights[imu]*func_at_ords[imu];
    for (int k=0; k<num_basis; ++k)
      f[k] += tmp*basis_at_ords[k+num_basis*imu];
  }
}

gkyl_proj_maxwellian_pots_on_basis* gkyl_proj_maxwellian_pots_on_basis_new(const struct gkyl_rect_grid *grid, 
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, 
  int num_quad) 
{
  gkyl_proj_maxwellian_pots_on_basis *up = gkyl_malloc(sizeof(gkyl_proj_maxwellian_pots_on_basis));
  up->grid = *grid;
  up->cdim = conf_basis->ndim;
  up->pdim = phase_basis->ndim;

  up->num_conf_basis = conf_basis->num_basis;
  up->num_phase_basis = phase_basis->num_basis;

  // Initialize data neede for configuration space quadrature
  up->tot_conf_quad = init_quad_values(up->cdim, conf_basis, num_quad,
    &up->conf_ordinates, &up->conf_weights, &up->conf_basis_at_ords);

  // Initialize data neede for phase space quadrature
  up->tot_quad = init_quad_values(up->pdim, phase_basis, num_quad,
    &up->ordinates, &up->weights, &up->basis_at_ords);

  up->fpo_h_at_ords = gkyl_array_new(GKYL_DOUBLE, 1, up->tot_quad);
  up->fpo_g_at_ords = gkyl_array_new(GKYL_DOUBLE, 1, up->tot_quad);

  up->conf_qrange = get_qrange(up->cdim, up->cdim, num_quad);
  up->phase_qrange = get_qrange(up->cdim, up->pdim, num_quad);
  return up;
}

void gkyl_proj_maxwellian_pots_on_basis_lab_mom(const gkyl_proj_maxwellian_pots_on_basis *up,
    const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
    const struct gkyl_array *moms, const struct gkyl_array *gamma, const double mass, struct gkyl_array *max_pots) 
{
  // Calculate Maxwellian potentials using primitive moments
  int cdim = up->cdim, pdim = up->pdim;
  int vdim = pdim-cdim;
  int tot_quad = up->tot_quad;
  int num_phase_basis = up->num_phase_basis;

  int tot_conf_quad = up->tot_conf_quad;
  int num_conf_basis = up->num_conf_basis;

  struct gkyl_range vel_range;
  struct gkyl_range_iter conf_iter, vel_iter;

  int pidx[GKYL_MAX_DIM], rem_dir[GKYL_MAX_DIM] = { 0 };
  for (int d=0; d<conf_range->ndim; ++d) rem_dir[d] = 1;

  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM];
  double den[tot_conf_quad], udrift[tot_conf_quad][vdim], temp[tot_conf_quad], gamma_ev[tot_conf_quad];
  double sqrt2temp_over_m[tot_conf_quad], fpo_g_prefact[tot_conf_quad];

  // Outer loop over configuration space cells
  gkyl_range_iter_init(&conf_iter, conf_range);
  while (gkyl_range_iter_next(&conf_iter)) {
    long midx = gkyl_range_idx(conf_range, conf_iter.idx);

    const double *moms_d = gkyl_array_cfetch(moms, midx);
    const double *m1i_d = &moms_d[num_conf_basis];
    const double *m2_d = &moms_d[num_conf_basis*(vdim+1)];
    const double *gamma_d = gkyl_array_cfetch(gamma, midx);

    // Compute primitive moments at quadrature nodes
    for (int n=0; n<tot_conf_quad; ++n) {
      const double *b_conf_ord = gkyl_array_cfetch(up->conf_basis_at_ords, n);

      // Number density and gamma evaluation
      double den_n = 0.0, gamma_evn = 0.0;
      for (int k=0; k<num_conf_basis; ++k) {
        den_n += moms_d[k]*b_conf_ord[k];
        gamma_evn += gamma_d[k]*b_conf_ord[k];
      }
      den[n] = den_n;
      gamma_ev[n] = gamma_evn;

      // Bulk velocity vector
      for (int d=0; d<vdim; ++d) {
        double m1i_n = 0.0;
        for (int k=0; k<num_conf_basis; ++k) 
          m1i_n += m1i_d[num_conf_basis*d+k]*b_conf_ord[k];
        udrift[n][d] = m1i_n/den[n];
      }

      // vth^2 and drift velocity to get temperature
      double m2_n = 0.0;
      for (int k=0; k<num_conf_basis; ++k)
         m2_n += m2_d[k]*b_conf_ord[k];

      double usq = 0.0;
      for (int d=0; d<vdim; ++d)
        usq += udrift[n][d]*udrift[n][d];

      double vtsq = (m2_n - den[n]*usq)/(den[n]*vdim);
      temp[n] = vtsq*mass;

      sqrt2temp_over_m[n] = sqrt(2.0*vtsq);
      fpo_g_prefact[n] = gamma_evn*den_n*sqrt2temp_over_m[n];
    }
    
    // Inner loop over velocity space
    gkyl_range_deflate(&vel_range, phase_range, rem_dir, conf_iter.idx);
    gkyl_range_iter_no_split_init(&vel_iter, &vel_range);
    while (gkyl_range_iter_next(&vel_iter)) {
      copy_idx_arrays(conf_range->ndim, phase_range->ndim, conf_iter.idx, vel_iter.idx, pidx);
      gkyl_rect_grid_cell_center(&up->grid, pidx, xc);
      
      struct gkyl_range_iter qiter;
      gkyl_range_iter_init(&qiter, &up->phase_qrange);
      // Compute potentials at phase space quadrature nodes
      while (gkyl_range_iter_next(&qiter)) {
        int cqidx = gkyl_range_idx(&up->conf_qrange, qiter.idx);
        int pqidx = gkyl_range_idx(&up->phase_qrange, qiter.idx);

        comp_to_phys(pdim, gkyl_array_cfetch(up->ordinates, pqidx), up->grid.dx, xc, xmu);

        double rel_speedsq_q = 0.0;
        for (int d=0; d<vdim; ++d)
          rel_speedsq_q += pow(udrift[cqidx][d]-xmu[cdim+d],2);

        // Retrieve values of necessary quantities at quadrature node
        double rel_speed_q = sqrt(rel_speedsq_q);
        double den_q = den[cqidx];
        double temp_q = temp[cqidx];
        double gamma_q = gamma_ev[cqidx];
        double sqrt2temp_over_m_q = sqrt2temp_over_m[cqidx];

        // Compute H and G at quadrature node
        double *fpo_h_q = gkyl_array_fetch(up->fpo_h_at_ords, pqidx);
        double *fpo_g_q = gkyl_array_fetch(up->fpo_g_at_ords, pqidx);
        
        fpo_h_q[0] = gamma_q*den_q/rel_speed_q * erf(rel_speed_q/sqrt2temp_over_m_q);
        fpo_g_q[0] = fpo_g_prefact[cqidx]*(1.0/(sqrt(GKYL_PI))*exp(-mass*rel_speedsq_q/(2.0*temp_q)) + 
          erf(rel_speed_q/sqrt2temp_over_m_q)*(sqrt2temp_over_m_q/rel_speed_q + rel_speed_q/sqrt2temp_over_m_q));
      }

      // Project potentials onto basis
      long lidx = gkyl_range_idx(&vel_range, vel_iter.idx);
      double *fpo_h_d = gkyl_array_fetch(max_pots, lidx);
      double *fpo_g_d = &fpo_h_d[num_phase_basis];
      proj_on_basis(up, up->fpo_h_at_ords, fpo_h_d);
      proj_on_basis(up, up->fpo_g_at_ords, fpo_g_d);
    }
  }
}


