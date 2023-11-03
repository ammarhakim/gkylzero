#include "gkyl_eval_on_nodes.h"
#include <math.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_const.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_proj_maxwellian_pots_on_basis.h>
#include <gkyl_proj_maxwellian_pots_on_basis_priv.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

static inline double eval_fpo_h(double gamma_q, double den_q, 
  double rel_speed_q, double sqrt2temp_over_m_q) 
{
 return gamma_q*den_q/rel_speed_q * erf(rel_speed_q/sqrt2temp_over_m_q);
}

static inline double eval_fpo_g(double gamma_q, double den_q, double rel_speed_q, 
  double sqrt2temp_over_m_q) 
{
  double rel_speedsq_q = pow(rel_speed_q, 2);
  return  gamma_q*den_q*sqrt2temp_over_m_q*
    (1.0/(sqrt(GKYL_PI))*exp(-rel_speedsq_q/pow(sqrt2temp_over_m_q, 2)) + 
    erf(rel_speed_q/sqrt2temp_over_m_q)*(sqrt2temp_over_m_q/(2.0*rel_speed_q) + 
    rel_speed_q/sqrt2temp_over_m_q));
}

static inline double eval_fpo_dhdv(double gamma_q, double den_q, 
  double rel_vel_in_dir_q, double sqrt2temp_over_m_q, double rel_speed_q) 
{
  double rel_speedsq_q = pow(rel_speed_q, 2);
  return gamma_q*den_q*rel_vel_in_dir_q * (
    2.0*exp(-rel_speedsq_q/pow(sqrt2temp_over_m_q,2))/(sqrt(GKYL_PI)*sqrt2temp_over_m_q*rel_speedsq_q) -
    erf(rel_speed_q/sqrt2temp_over_m_q)/pow(rel_speedsq_q, 1.5)
  );
}

static inline double eval_fpo_dgdv(double gamma_q, double den_q, double rel_vel_in_dir_q,
  double sqrt2temp_over_m_q, double rel_speed_q) {
  double rel_speedsq_q = pow(rel_speed_q,2);
  return gamma_q*den_q*rel_vel_in_dir_q * (
    exp(-rel_speedsq_q/pow(sqrt2temp_over_m_q,2))*sqrt2temp_over_m_q/(sqrt(GKYL_PI)*rel_speedsq_q) -
    erf(rel_speed_q/sqrt2temp_over_m_q)*(pow(sqrt2temp_over_m_q,2)/(2.0*pow(rel_speed_q,3)) - 1.0/rel_speed_q)
  );
}

static inline double eval_fpo_d2gdv2(double gamma_q, double den_q, 
  double rel_vel_in_dir_q, double sqrt2temp_over_m_q, double rel_speed_q) {
  double rel_speedsq_q = pow(rel_speed_q,2);
  return gamma_q*den_q*(exp(-rel_speedsq_q/pow(sqrt2temp_over_m_q,2))/sqrt(GKYL_PI)*(sqrt2temp_over_m_q/rel_speedsq_q - 
    pow(rel_vel_in_dir_q,2)*2.0*sqrt2temp_over_m_q*(1.0/pow(rel_speedsq_q, 2) +
    1.0/(pow(sqrt2temp_over_m_q,2)*rel_speedsq_q) + (pow(sqrt2temp_over_m_q,2) - 
    2.0*rel_speedsq_q)/(2.0*pow(sqrt2temp_over_m_q,2)*pow(rel_speedsq_q, 2)))) +
    erf(rel_speed_q/sqrt2temp_over_m_q)*(pow(rel_vel_in_dir_q,2)*
    (3.0*pow(sqrt2temp_over_m_q,2)-2.0*rel_speedsq_q)/(2.0*pow(rel_speed_q, 5)) - 
    pow(sqrt2temp_over_m_q, 2)/(2.0*pow(rel_speed_q, 3)) + 1.0/rel_speed_q));
}

static inline double eval_fpo_d2gdv2_cross(double gamma_q, double den_q,
  double rel_vel_in_dir1_q, double rel_vel_in_dir2_q, double rel_speed_q,
  double sqrt2temp_over_m_q) {
  double rel_speedsq_q = pow(rel_speed_q,2);
  return gamma_q*den_q*rel_vel_in_dir1_q*rel_vel_in_dir2_q/pow(rel_speedsq_q,2)*(
    erf(rel_speed_q/sqrt2temp_over_m_q)*(3.0*pow(sqrt2temp_over_m_q,2)-2.0*rel_speedsq_q)/(2.0*rel_speed_q) -
    3.0*exp(-rel_speedsq_q/pow(sqrt2temp_over_m_q,2))*sqrt2temp_over_m_q/sqrt(GKYL_PI)
  );
}

// create range to loop over quadrature points.
static inline struct gkyl_range get_qrange(int cdim, int dim, int num_quad, int num_quad_v, bool *is_vdim_p2) {
  int qshape[GKYL_MAX_DIM];
  for (int i=0; i<cdim; ++i) qshape[i] = num_quad;
  for (int i=cdim; i<dim; ++i) qshape[i] = is_vdim_p2[i-cdim]? num_quad_v : num_quad;
  struct gkyl_range qrange;
  gkyl_range_init_from_shape(&qrange, dim, qshape);
  return qrange;
}

// Sets ordinates, weights and basis functions at ords.
// Returns the total number of quadrature nodes
static int
init_quad_values(int cdim, const struct gkyl_basis *basis, int num_quad, struct gkyl_array **ordinates,
  struct gkyl_array **weights, struct gkyl_array **basis_at_ords, bool use_gpu)
{
  int ndim = basis->ndim;
  int vdim = ndim-cdim;
  int num_quad_v = num_quad;
  // hybrid basis have p=2 in velocity space.
  bool is_vdim_p2[] = {false, false, false};  // 3 is the max vdim.
  if (basis->b_type == GKYL_BASIS_MODAL_HYBRID) {
    num_quad_v = num_quad+1;
    for (int d=0; d<vdim; d++) is_vdim_p2[d] = true;
  }

  double ordinates1[num_quad], weights1[num_quad];
  double ordinates1_v[num_quad_v], weights1_v[num_quad_v];
  if (num_quad <= gkyl_gauss_max) {
    // use pre-computed values if possible (these are more accurate
    // than computing them on the fly)
    memcpy(ordinates1, gkyl_gauss_ordinates[num_quad], sizeof(double[num_quad]));
    memcpy(weights1, gkyl_gauss_weights[num_quad], sizeof(double[num_quad]));
  } else {
    gkyl_gauleg(-1, 1, ordinates1, weights1, num_quad);
  }
  if (num_quad_v <= gkyl_gauss_max) {
    memcpy(ordinates1_v, gkyl_gauss_ordinates[num_quad_v], sizeof(double[num_quad_v]));
    memcpy(weights1_v, gkyl_gauss_weights[num_quad_v], sizeof(double[num_quad_v]));
  } else {
    gkyl_gauleg(-1, 1, ordinates1_v, weights1_v, num_quad_v);
  }

  struct gkyl_range qrange = get_qrange(cdim, ndim, num_quad, num_quad_v, is_vdim_p2);

  int tot_quad = qrange.volume;

  // create ordinates and weights for multi-D quadrature
  struct gkyl_array *ordinates_ho = gkyl_array_new(GKYL_DOUBLE, ndim, tot_quad);
  struct gkyl_array *weights_ho = gkyl_array_new(GKYL_DOUBLE, 1, tot_quad);
  if (use_gpu) {
    *ordinates = gkyl_array_cu_dev_new(GKYL_DOUBLE, ndim, tot_quad);
    *weights = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, tot_quad);
  } else {
    *ordinates = gkyl_array_new(GKYL_DOUBLE, ndim, tot_quad);
    *weights = gkyl_array_new(GKYL_DOUBLE, 1, tot_quad);
  }

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &qrange);

  while (gkyl_range_iter_next(&iter)) {
    int node = gkyl_range_idx(&qrange, iter.idx);
    
    // set ordinates
    double *ord = gkyl_array_fetch(ordinates_ho, node);
    for (int i=0; i<cdim; ++i)
      ord[i] = ordinates1[iter.idx[i]-qrange.lower[i]];
    for (int i=cdim; i<ndim; ++i)
      ord[i] = is_vdim_p2[i-cdim]? ordinates1_v[iter.idx[i]-qrange.lower[i]] :
                                   ordinates1[iter.idx[i]-qrange.lower[i]];
    
    // set weights
    double *wgt = gkyl_array_fetch(weights_ho, node);
    wgt[0] = 1.0;
    for (int i=0; i<cdim; ++i)
      wgt[0] *= weights1[iter.idx[i]-qrange.lower[i]];
    for (int i=cdim; i<ndim; ++i)
      wgt[0] *= is_vdim_p2[i-cdim]? weights1_v[iter.idx[i]-qrange.lower[i]] :
                                    weights1[iter.idx[i]-qrange.lower[i]];
  }

  // pre-compute basis functions at ordinates
  struct gkyl_array *basis_at_ords_ho = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, tot_quad);
  if (use_gpu) 
    *basis_at_ords = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis->num_basis, tot_quad);
  else
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

static void
proj_on_surf_basis(const gkyl_proj_maxwellian_pots_on_basis *up, int offset, const struct gkyl_array *fun_at_ords, double* f)
{
  int num_surf_basis = up->num_surf_basis;
  int tot_surf_quad = up->tot_surf_quad;

  const double* GKYL_RESTRICT weights = up->surf_weights->data;
  const double* GKYL_RESTRICT surf_basis_at_ords = up->surf_basis_at_ords->data;
  const double* GKYL_RESTRICT func_at_ords = fun_at_ords->data;

  for (int k=0; k<num_surf_basis; ++k) f[offset*num_surf_basis+k] = 0.0;
  
  for (int imu=0; imu<tot_surf_quad; ++imu) {
    double tmp = weights[imu]*func_at_ords[imu];
    for (int k=0; k<num_surf_basis; ++k)
      f[offset*num_surf_basis+k] += tmp*surf_basis_at_ords[num_surf_basis*imu+k];
  }
}

static void
eval_on_nodes_nod2mod(int num_ret_vals, const struct gkyl_basis *basis, const struct gkyl_array *fun_at_nodes, double *f) {
  const double *fao = gkyl_array_cfetch(fun_at_nodes, 0);

  int num_basis = basis->num_basis;
  double fnodal[num_basis];
  for (int i=0; i<num_ret_vals; ++i) {
    for (int k=0; k<num_basis; ++k) {
      fnodal[k] = fao[num_ret_vals*k+i];
    }

    basis->nodal_to_modal(fnodal, &f[num_basis*i]);
  }
}

gkyl_proj_maxwellian_pots_on_basis* gkyl_proj_maxwellian_pots_on_basis_new(const struct gkyl_rect_grid *grid, 
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, const struct gkyl_basis *surf_basis,
  int num_quad) 
{
  gkyl_proj_maxwellian_pots_on_basis *up = gkyl_malloc(sizeof(gkyl_proj_maxwellian_pots_on_basis));
  up->grid = *grid;
  up->cdim = conf_basis->ndim;
  up->pdim = phase_basis->ndim;
  up->num_quad = num_quad;

  up->surf_basis = surf_basis;

  up->num_conf_basis = conf_basis->num_basis;
  up->num_phase_basis = phase_basis->num_basis;
  up->num_surf_basis = surf_basis->num_basis;

  up->surf_nodes = gkyl_array_new(GKYL_DOUBLE, grid->ndim-1, surf_basis->num_basis);
  surf_basis->node_list(gkyl_array_fetch(up->surf_nodes, 0));

  bool use_gpu = false;

  // Initialize data needed for configuration space quadrature
  up->tot_conf_quad = init_quad_values(up->cdim, conf_basis, num_quad,
    &up->conf_ordinates, &up->conf_weights, &up->conf_basis_at_ords, use_gpu);

  // Initialize data needed for phase space quadrature
  up->tot_quad = init_quad_values(up->cdim, phase_basis, num_quad,
    &up->ordinates, &up->weights, &up->basis_at_ords, use_gpu);

  // Initialize quadrature for surface expansion
  up->tot_surf_quad = init_quad_values(up->cdim, surf_basis, num_quad, 
    &up->surf_ordinates, &up->surf_weights, &up->surf_basis_at_ords, use_gpu);

  up->fpo_h_at_ords = gkyl_array_new(GKYL_DOUBLE, 1, up->tot_quad);
  up->fpo_g_at_ords = gkyl_array_new(GKYL_DOUBLE, 1, up->tot_quad);

  up->fpo_h_at_surf_ords = gkyl_array_new(GKYL_DOUBLE, 1, up->tot_surf_quad);
  up->fpo_g_at_surf_ords = gkyl_array_new(GKYL_DOUBLE, 1, up->tot_surf_quad);
  up->fpo_dhdv_at_surf_ords = gkyl_array_new(GKYL_DOUBLE, 1, up->tot_surf_quad);
  up->fpo_dgdv_at_surf_ords = gkyl_array_new(GKYL_DOUBLE, 1, up->tot_surf_quad);
  up->fpo_d2gdv2_at_surf_ords = gkyl_array_new(GKYL_DOUBLE, 1, up->tot_surf_quad);

  // Hybrid basis support: uses p=2 in velocity space
  int vdim = up->pdim-up->cdim;
  int num_quad_v = num_quad;
  bool is_vdim_p2[] = {false, false, false};  // 3 is the max vdim.
  if ((phase_basis->b_type == GKYL_BASIS_MODAL_HYBRID) ||
      (phase_basis->b_type == GKYL_BASIS_MODAL_GKHYBRID)) {
    num_quad_v = num_quad+1;
    is_vdim_p2[0] = true;  // for gkhybrid.
    if (phase_basis->b_type == GKYL_BASIS_MODAL_HYBRID)
      for (int d=0; d<vdim; d++) is_vdim_p2[d] = true;
  }

  up->conf_qrange = get_qrange(up->cdim, up->cdim, num_quad, num_quad_v, is_vdim_p2);
  up->phase_qrange = get_qrange(up->cdim, up->pdim, num_quad, num_quad_v, is_vdim_p2);
  up->surf_qrange = get_qrange(up->cdim, up->pdim-1, num_quad, num_quad_v, is_vdim_p2);
  return up;
}

void gkyl_proj_maxwellian_pots_on_basis_lab_mom(const gkyl_proj_maxwellian_pots_on_basis *up,
    const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
    const struct gkyl_array *moms, const struct gkyl_array *gamma, const double mass, 
    struct gkyl_array *fpo_h, struct gkyl_array *fpo_g,
    struct gkyl_array *fpo_h_surf, struct gkyl_array *fpo_g_surf,
    struct gkyl_array *fpo_dhdv_surf, struct gkyl_array *fpo_dgdv_surf, struct gkyl_array *fpo_d2gdv2_surf)
{
  // Calculate Maxwellian potentials using primitive moments
  int cdim = up->cdim, pdim = up->pdim;
  int vdim = pdim-cdim;
  int tot_quad = up->tot_quad;
  int num_phase_basis = up->num_phase_basis;

  int tot_conf_quad = up->tot_conf_quad;
  int num_conf_basis = up->num_conf_basis;
  int num_surf_basis = up->num_surf_basis;

  struct gkyl_range vel_range;
  struct gkyl_range_iter conf_iter, vel_iter;

  int pidx[GKYL_MAX_DIM], rem_dir[GKYL_MAX_DIM] = { 0 };
  for (int d=0; d<conf_range->ndim; ++d) rem_dir[d] = 1;

  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM];
  double den[tot_conf_quad], udrift[tot_conf_quad][vdim], temp[tot_conf_quad], gamma_ev[tot_conf_quad];
  double sqrt2temp_over_m[tot_conf_quad]; 

  // Loop over configuration space cells for quad integration of moments
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
      den[n] = 0.0;
      for (int k=0; k<num_conf_basis; ++k) {
        den[n] += moms_d[k]*b_conf_ord[k];
      }

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
    }
  
      
 
    gkyl_range_deflate(&vel_range, phase_range, rem_dir, conf_iter.idx);
    gkyl_range_iter_no_split_init(&vel_iter, &vel_range);
    while (gkyl_range_iter_next(&vel_iter)) {
      copy_idx_arrays(conf_range->ndim, phase_range->ndim, conf_iter.idx, vel_iter.idx, pidx);
      gkyl_rect_grid_cell_center(&up->grid, pidx, xc);
      
      struct gkyl_range_iter qiter;
      gkyl_range_iter_init(&qiter, &up->phase_qrange);
  
      long lidx = gkyl_range_idx(&vel_range, vel_iter.idx);
  
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
        double gamma_q = gamma_d[0];
        double sqrt2temp_over_m_q = sqrt2temp_over_m[cqidx];
  
        // Compute H and G at quadrature node
        double *fpo_h_q = gkyl_array_fetch(up->fpo_h_at_ords, pqidx);
        double *fpo_g_q = gkyl_array_fetch(up->fpo_g_at_ords, pqidx);
  
        fpo_h_q[0] = eval_fpo_h(gamma_q, den_q, rel_speed_q, sqrt2temp_over_m_q); 
  
        fpo_g_q[0] = eval_fpo_g(gamma_q, den_q, rel_speed_q, sqrt2temp_over_m_q);
      }
  
      struct gkyl_array *fpo_g_at_nodes = gkyl_array_new(GKYL_DOUBLE, 1, num_surf_basis);
  
      // Check if we're at a velocity space boundary in each velocity direction
      for (int d=0; d<vdim; ++d) {
        int dir = d + cdim;
  
        if (pidx[dir] == vel_range.lower[d] || pidx[dir] == vel_range.upper[d]) { 
          // Velocity value at boundary for surface expansion
          int vmax = pidx[dir] == vel_range.lower[d] ? up->grid.lower[dir] : up->grid.upper[dir]; 
  
          for (int i=0; i<num_surf_basis; ++i) {
            comp_to_phys(pdim, gkyl_array_cfetch(up->surf_nodes, i), up->grid.dx, xc, xmu);
  
            xmu[dir] = vmax;
            double rel_speedsq_n = 0.0;
            for (int d=0; d<vdim; ++d)
              rel_speedsq_n += pow(xmu[cdim+d]-udrift[0][d], 2);
  
            double rel_vel_in_dir_n = xmu[dir]-udrift[0][d];
            double rel_speed_n = sqrt(rel_speedsq_n);
            double den_n = den[0];
            double temp_n = temp[0];
            double gamma_n = gamma_d[0];
            double sqrt2temp_over_m_n = sqrt2temp_over_m[0];
  
            double* fpo_g_at_nodes_n = gkyl_array_fetch(fpo_g_at_nodes, i);
  
            fpo_g_at_nodes_n[0] = eval_fpo_g(gamma_n, den_n, rel_speed_n, sqrt2temp_over_m_n);
          }
          
          double *fpo_dgdv_surf_n = gkyl_array_fetch(fpo_dgdv_surf, lidx);
          eval_on_nodes_nod2mod(1, up->surf_basis, fpo_g_at_nodes, &fpo_dgdv_surf_n[d*num_surf_basis]);
  
          struct gkyl_range_iter surf_qiter;
          gkyl_range_iter_init(&surf_qiter, &up->surf_qrange);
  
          // Iterate over surface quadrature points to calculate:
          // H, G, dH/dv, dG/dv
          while (gkyl_range_iter_next(&surf_qiter)) {
            int surf_qidx = gkyl_range_idx(&up->surf_qrange, surf_qiter.idx);
            int surf_cqidx = gkyl_range_idx(&up->conf_qrange, surf_qiter.idx);
  
            // Have to map pdim-1 surface quadrature index to pdim phase quadrature index
            // to get correct phase space variables.
            int edge_idx = pidx[dir] == vel_range.lower[d] ? 0 : up->num_quad-1;
            int phase_idx[GKYL_MAX_DIM];
            edge_idx_to_phase_idx(pdim, dir, surf_qiter.idx, edge_idx, phase_idx);
            int phase_lidx = gkyl_range_idx(&up->phase_qrange, phase_idx);
            comp_to_phys(pdim, gkyl_array_cfetch(up->ordinates, phase_lidx), up->grid.dx, xc, xmu);
  
            xmu[dir] = vmax;
            double rel_speedsq_q = 0.0;
            for (int d=0; d<vdim; ++d)
              rel_speedsq_q += pow(xmu[cdim+d]-udrift[surf_cqidx][d],2);
  
            // Retrieve values of necessary quantities at quadrature node
            double rel_vel_in_dir_q = xmu[dir]-udrift[surf_cqidx][d];
            double rel_speed_q = sqrt(rel_speedsq_q);
            double den_q = den[surf_cqidx];
            double temp_q = temp[surf_cqidx];
            double gamma_q = gamma_d[0];
            double sqrt2temp_over_m_q = sqrt2temp_over_m[surf_cqidx];
  
            double *fpo_h_at_surf_ords_q = gkyl_array_fetch(up->fpo_h_at_surf_ords, surf_qidx);
            double *fpo_g_at_surf_ords_q = gkyl_array_fetch(up->fpo_g_at_surf_ords, surf_qidx);
            double *fpo_dhdv_at_surf_ords_q = gkyl_array_fetch(up->fpo_dhdv_at_surf_ords, surf_qidx);
            double *fpo_dgdv_at_surf_ords_q = gkyl_array_fetch(up->fpo_dgdv_at_surf_ords, surf_qidx);
            double *fpo_d2gdv2_at_surf_ords_q = gkyl_array_fetch(up->fpo_d2gdv2_at_surf_ords, surf_qidx);
           
            fpo_h_at_surf_ords_q[0] = eval_fpo_h(gamma_q, den_q, rel_speed_q, sqrt2temp_over_m_q); 
  
            fpo_g_at_surf_ords_q[0] = eval_fpo_g(gamma_q, den_q, rel_speed_q, sqrt2temp_over_m_q);
  
            fpo_dhdv_at_surf_ords_q[0] = eval_fpo_dhdv(gamma_q, den_q, rel_vel_in_dir_q, 
              sqrt2temp_over_m_q, rel_speed_q); 
  
            fpo_dgdv_at_surf_ords_q[0] = eval_fpo_dgdv(gamma_q, den_q, rel_vel_in_dir_q, 
              sqrt2temp_over_m_q, rel_speed_q);
          }
  
          // Iterate over second velocity space direction to calculate d2G/dv2
          for (int d2=0; d2<vdim; ++d2) {
            int dir2 = d2 + cdim;
            // Iterating over velocity space direction pairs within a phase space cell
            // relevant constants are defined at quadrature nodes?
            // int num_surf_basis = up->num_surf_basis;
            // struct gkyl_array *fpo_g_at_nodes = gkyl_array_new(GKYL_DOUBLE, 1, num_surf_basis);
  
            gkyl_range_iter_init(&surf_qiter, &up->surf_qrange);
            while (gkyl_range_iter_next(&surf_qiter)) {
              int surf_qidx = gkyl_range_idx(&up->surf_qrange, surf_qiter.idx);
              int surf_cqidx = gkyl_range_idx(&up->conf_qrange, surf_qiter.idx);
  
              // Have to map pdim-1 surface quadrature index to pdim phase quadrature index
              // to get correct phase space variables.
              int edge_idx = pidx[dir] == vel_range.lower[d] ? 0 : up->num_quad-1;
              int phase_idx[GKYL_MAX_DIM];
              edge_idx_to_phase_idx(pdim, dir, surf_qiter.idx, edge_idx, phase_idx);
              int phase_lidx = gkyl_range_idx(&up->phase_qrange, phase_idx);
              comp_to_phys(pdim, gkyl_array_cfetch(up->ordinates, phase_lidx), up->grid.dx, xc, xmu);
  
              xmu[dir] = vmax;
              double rel_speedsq_q = 0.0;
              for (int d=0; d<vdim; ++d)
                rel_speedsq_q += pow(xmu[cdim+d]-udrift[surf_cqidx][d],2);
  
              // Retrieve values of necessary quantities at quadrature node
              double rel_vel_in_dir1_q = xmu[dir]-udrift[surf_cqidx][d];
              double rel_vel_in_dir2_q = xmu[dir2]-udrift[surf_cqidx][d2];
              double rel_speed_q = sqrt(rel_speedsq_q);
              double den_q = den[surf_cqidx];
              double temp_q = temp[surf_cqidx];
              double gamma_q = gamma_d[0];
              double sqrt2temp_over_m_q = sqrt2temp_over_m[surf_cqidx];
  
              double *fpo_d2gdv2_at_surf_ords_q = gkyl_array_fetch(up->fpo_d2gdv2_at_surf_ords, surf_qidx);
  
              if (d == d2) {
                fpo_d2gdv2_at_surf_ords_q[0] = eval_fpo_d2gdv2(gamma_q, den_q, 
                  rel_vel_in_dir1_q, 
                  sqrt2temp_over_m_q, rel_speed_q);
              }
              else {
                fpo_d2gdv2_at_surf_ords_q[0] = eval_fpo_d2gdv2_cross(gamma_q, den_q,
                  rel_vel_in_dir1_q, rel_vel_in_dir2_q, rel_speed_q, sqrt2temp_over_m_q); 
              }
            }
            proj_on_surf_basis(up, d*vdim+d2, up->fpo_d2gdv2_at_surf_ords, gkyl_array_fetch(fpo_d2gdv2_surf, lidx));
          }
  
          // Project surface expansions onto surface basis
          proj_on_surf_basis(up, d, up->fpo_h_at_surf_ords, gkyl_array_fetch(fpo_h_surf, lidx));
          proj_on_surf_basis(up, d, up->fpo_g_at_surf_ords, gkyl_array_fetch(fpo_g_surf, lidx));
          proj_on_surf_basis(up, d, up->fpo_dhdv_at_surf_ords, gkyl_array_fetch(fpo_dhdv_surf, lidx));
          // proj_on_surf_basis(up, d, up->fpo_dgdv_at_surf_ords, gkyl_array_fetch(fpo_dgdv_surf, lidx))    
        }
      } 
      // Project potentials onto basis
      proj_on_basis(up, up->fpo_h_at_ords, gkyl_array_fetch(fpo_h, lidx));
      proj_on_basis(up, up->fpo_g_at_ords, gkyl_array_fetch(fpo_g, lidx));
      gkyl_array_release(fpo_g_at_nodes);
    }
  }
}

void gkyl_proj_maxwellian_pots_on_basis_release(gkyl_proj_maxwellian_pots_on_basis *up)
{
  gkyl_array_release(up->ordinates);
  gkyl_array_release(up->weights);
  gkyl_array_release(up->basis_at_ords);

  gkyl_array_release(up->conf_ordinates);
  gkyl_array_release(up->conf_weights);
  gkyl_array_release(up->conf_basis_at_ords);

  gkyl_array_release(up->surf_ordinates);
  gkyl_array_release(up->surf_weights);
  gkyl_array_release(up->surf_basis_at_ords); 

  gkyl_array_release(up->fpo_h_at_ords);
  gkyl_array_release(up->fpo_g_at_ords);

  gkyl_array_release(up->fpo_h_at_surf_ords);
  gkyl_array_release(up->fpo_g_at_surf_ords);
  gkyl_array_release(up->fpo_dhdv_at_surf_ords);
  gkyl_array_release(up->fpo_dgdv_at_surf_ords);
  gkyl_array_release(up->fpo_d2gdv2_at_surf_ords);

  gkyl_free(up);
}
