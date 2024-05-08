#include <string.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_const.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_vlasov_lte_proj_on_basis.h>
#include <gkyl_vlasov_lte_proj_on_basis_priv.h>
#include <gkyl_range.h>
#include <assert.h>

// create range to loop over quadrature points.
static inline struct gkyl_range get_qrange(int cdim, int dim, int num_quad, int num_quad_v, bool *is_vdim_p2) {
  int qshape[GKYL_MAX_DIM];
  for (int i=0; i<cdim; ++i) qshape[i] = num_quad;
  for (int i=cdim; i<dim; ++i) qshape[i] = is_vdim_p2[i-cdim] ? num_quad_v : num_quad;
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
    for (int d=0; d<vdim; d++) {
      is_vdim_p2[d] = true;
    }
  }

  double ordinates1[num_quad], weights1[num_quad];
  double ordinates1_v[num_quad_v], weights1_v[num_quad_v];
  if (num_quad <= gkyl_gauss_max) {
    // use pre-computed values if possible (these are more accurate
    // than computing them on the fly)
    memcpy(ordinates1, gkyl_gauss_ordinates[num_quad], sizeof(double[num_quad]));
    memcpy(weights1, gkyl_gauss_weights[num_quad], sizeof(double[num_quad]));
  } 
  else {
    gkyl_gauleg(-1, 1, ordinates1, weights1, num_quad);
  }
  if (num_quad_v <= gkyl_gauss_max) {
    memcpy(ordinates1_v, gkyl_gauss_ordinates[num_quad_v], sizeof(double[num_quad_v]));
    memcpy(weights1_v, gkyl_gauss_weights[num_quad_v], sizeof(double[num_quad_v]));
  } 
  else {
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
  } 
  else {
    *ordinates = gkyl_array_new(GKYL_DOUBLE, ndim, tot_quad);
    *weights = gkyl_array_new(GKYL_DOUBLE, 1, tot_quad);
  }

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &qrange);

  while (gkyl_range_iter_next(&iter)) {
    int node = gkyl_range_idx(&qrange, iter.idx);
    
    // set ordinates
    double *ord = gkyl_array_fetch(ordinates_ho, node);
    for (int i=0; i<cdim; ++i) {
      ord[i] = ordinates1[iter.idx[i]-qrange.lower[i]];
    }
    for (int i=cdim; i<ndim; ++i) {
      ord[i] = is_vdim_p2[i-cdim] ? 
        ordinates1_v[iter.idx[i]-qrange.lower[i]] : ordinates1[iter.idx[i]-qrange.lower[i]];
    }
    
    // set weights
    double *wgt = gkyl_array_fetch(weights_ho, node);
    wgt[0] = 1.0;
    for (int i=0; i<cdim; ++i) {
      wgt[0] *= weights1[iter.idx[i]-qrange.lower[i]];
    }
    for (int i=cdim; i<ndim; ++i) {
      wgt[0] *= is_vdim_p2[i-cdim] ? 
        weights1_v[iter.idx[i]-qrange.lower[i]] : weights1[iter.idx[i]-qrange.lower[i]];
    }
  }

  // pre-compute basis functions at ordinates
  struct gkyl_array *basis_at_ords_ho = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, tot_quad);
  if (use_gpu) {
    *basis_at_ords = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis->num_basis, tot_quad);
  }
  else {
    *basis_at_ords = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, tot_quad);
  }

  for (int n=0; n<tot_quad; ++n) {
    basis->eval(gkyl_array_fetch(ordinates_ho, n), gkyl_array_fetch(basis_at_ords_ho, n));
  }

  // copy host array to device array
  gkyl_array_copy(*ordinates, ordinates_ho);
  gkyl_array_copy(*weights, weights_ho);
  gkyl_array_copy(*basis_at_ords, basis_at_ords_ho);

  gkyl_array_release(ordinates_ho);
  gkyl_array_release(weights_ho);
  gkyl_array_release(basis_at_ords_ho);

  return tot_quad;
}

struct gkyl_vlasov_lte_proj_on_basis* 
gkyl_vlasov_lte_proj_on_basis_inew(const struct gkyl_vlasov_lte_proj_on_basis_inp *inp)
{
  gkyl_vlasov_lte_proj_on_basis *up = gkyl_malloc(sizeof(*up));

  up->phase_grid = *inp->phase_grid;
  up->conf_basis = *inp->conf_basis;
  up->phase_basis = *inp->phase_basis;
  up->h_ij_inv = inp->h_ij_inv;

  up->cdim = up->conf_basis.ndim;
  up->pdim = up->phase_basis.ndim;
  int vdim = up->pdim - up->cdim;
  up->num_conf_basis = up->conf_basis.num_basis;
  up->num_phase_basis = up->phase_basis.num_basis;
  up->use_gpu = inp->use_gpu;

  up->is_relativistic = false;
  if (inp->model_id == GKYL_MODEL_SR) {
    up->is_relativistic = true;
  }
  up->is_canonical_pb = false;
  if (inp->model_id == GKYL_MODEL_CANONICAL_PB) {
    up->is_canonical_pb = true;
  }

  // JJ 2024/03/23: device kernel has arrays hard-coded to 3x, vdim=3, p=2 for now.
  if (up->use_gpu) {
    assert(up->cdim<3 && up->conf_basis.poly_order<3);
  }

  int num_quad = up->conf_basis.poly_order+1;
  // initialize data needed for conf-space quadrature 
  up->tot_conf_quad = init_quad_values(up->cdim, &up->conf_basis, num_quad,
    &up->conf_ordinates, &up->conf_weights, &up->conf_basis_at_ords, up->use_gpu);

  // initialize data needed for phase-space quadrature 
  up->tot_quad = init_quad_values(up->cdim, &up->phase_basis, num_quad,
    &up->ordinates, &up->weights, &up->basis_at_ords, up->use_gpu);

  up->fun_at_ords = gkyl_array_new(GKYL_DOUBLE, 1, up->tot_quad); // Only used in CPU implementation.

  // To avoid creating iterators over ranges in device kernel, we'll
  // create a map between phase-space and conf-space ordinates.
  int num_quad_v = num_quad;  // Hybrid basis have p=2 in velocity space.
  // hybrid basis have p=2 in velocity space.
  bool is_vdim_p2[] = {false, false, false};  // 3 is the max vdim.
  if (up->phase_basis.b_type == GKYL_BASIS_MODAL_HYBRID) {
    num_quad_v = num_quad+1;
    for (int d=0; d<vdim; d++) {
      is_vdim_p2[d] = true;
    }
  }
  up->conf_qrange = get_qrange(up->cdim, up->cdim, num_quad, num_quad_v, is_vdim_p2);
  up->phase_qrange = get_qrange(up->cdim, up->pdim, num_quad, num_quad_v, is_vdim_p2);

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    // Allocate device copies of arrays needed for quadrature.

    int p2c_qidx_ho[up->phase_qrange.volume];
    up->p2c_qidx = (int*) gkyl_cu_malloc(sizeof(int)*up->phase_qrange.volume);

    int pidx[GKYL_MAX_DIM];
    for (int n=0; n<up->tot_quad; ++n) {
      gkyl_range_inv_idx(&up->phase_qrange, n, pidx);
      int cqidx = gkyl_range_idx(&up->conf_qrange, pidx);
      p2c_qidx_ho[n] = cqidx;
    }
    gkyl_cu_memcpy(up->p2c_qidx, p2c_qidx_ho, sizeof(int)*up->phase_qrange.volume, GKYL_CU_MEMCPY_H2D);
  }
#endif

  // Store a LTE moment calculation updater to compute and correct the density
  struct gkyl_vlasov_lte_moments_inp inp_mom = {
    .phase_grid = inp->phase_grid,
    .conf_basis = inp->conf_basis,
    .phase_basis = inp->phase_basis,
    .conf_range =  inp->conf_range,
    .conf_range_ext = inp->conf_range_ext,
    .vel_range = inp->vel_range,
    .p_over_gamma = inp->p_over_gamma,
    .gamma = inp->gamma,
    .gamma_inv = inp->gamma_inv,
    .h_ij_inv = inp->h_ij_inv,
    .model_id = inp->model_id,
    .mass = inp->mass,
    .use_gpu = inp->use_gpu,
  };
  up->moments_up = gkyl_vlasov_lte_moments_inew( &inp_mom );

  long conf_local_ncells = inp->conf_range->volume;
  long conf_local_ext_ncells = inp->conf_range_ext->volume;

  // Number density ratio: num_ratio = n_target/n0 and bin_op memory to compute ratio
  // Used for fixing the density with simple rescaling
  if (up->use_gpu) {
    up->num_ratio = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->conf_basis.num_basis, conf_local_ext_ncells);
    up->mem = gkyl_dg_bin_op_mem_cu_dev_new(conf_local_ncells, up->conf_basis.num_basis);
  }
  else {
    up->num_ratio = gkyl_array_new(GKYL_DOUBLE, up->conf_basis.num_basis, conf_local_ext_ncells);
    up->mem = gkyl_dg_bin_op_mem_new(conf_local_ncells, up->conf_basis.num_basis);
  }

  return up;
}

static void
proj_on_basis(const gkyl_vlasov_lte_proj_on_basis *up, const struct gkyl_array *fun_at_ords, double* f)
{
  int num_basis = up->num_phase_basis;
  int tot_quad = up->tot_quad;

  const double* GKYL_RESTRICT weights = up->weights->data;
  const double* GKYL_RESTRICT basis_at_ords = up->basis_at_ords->data;
  const double* GKYL_RESTRICT func_at_ords = fun_at_ords->data;

  for (int k=0; k<num_basis; ++k) f[k] = 0.0;
  
  for (int imu=0; imu<tot_quad; ++imu) {
    double tmp = weights[imu]*func_at_ords[imu];
    for (int k=0; k<num_basis; ++k) {
      f[k] += tmp*basis_at_ords[k+num_basis*imu];
    }
  }
}

void
gkyl_vlasov_lte_proj_on_basis_advance(gkyl_vlasov_lte_proj_on_basis *up,
  const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
  const struct gkyl_array *moms_lte, struct gkyl_array *f_lte)
{

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    return gkyl_vlasov_lte_proj_on_basis_advance_cu(up, phase_range, conf_range, moms_lte, f_lte);
#endif

  double f_floor = 1.e-40;  
  int cdim = up->cdim, pdim = up->pdim;
  int vdim = pdim-cdim;
  int tot_quad = up->tot_quad;
  int num_phase_basis = up->num_phase_basis;  

  int tot_conf_quad = up->tot_conf_quad;
  int num_conf_basis = up->num_conf_basis;

  struct gkyl_range vel_rng;
  struct gkyl_range_iter conf_iter, vel_iter;

  int pidx[GKYL_MAX_DIM], rem_dir[GKYL_MAX_DIM] = { 0 };
  for (int d=0; d<conf_range->ndim; ++d) rem_dir[d] = 1;

  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM];
  double n_quad[tot_conf_quad], V_drift_quad[tot_conf_quad][vdim], T_over_m_quad[tot_conf_quad];
  double h_ij_inv_quad[tot_conf_quad][vdim*(vdim + 1)/2];
  double V_drift_quad_cell_avg[tot_conf_quad][vdim];
  double expamp_quad[tot_conf_quad];

  // outer loop over configuration space cells; for each
  // config-space cell inner loop walks over velocity space
  gkyl_range_iter_init(&conf_iter, conf_range);
  while (gkyl_range_iter_next(&conf_iter)) {
    long midx = gkyl_range_idx(conf_range, conf_iter.idx);

    const double *moms_lte_d = gkyl_array_cfetch(moms_lte, midx);
    const double *n_d = moms_lte_d;
    const double *V_drift_d = &moms_lte_d[num_conf_basis];
    const double *T_over_m_d = &moms_lte_d[num_conf_basis*(vdim+1)];
    const double *h_ij_inv_d;
    if (up->is_canonical_pb) {
      h_ij_inv_d = gkyl_array_cfetch(up->h_ij_inv, midx);
    }

    // Sum over basis for given LTE moments (n, V_drift, T/m) in the stationary frame
    for (int n=0; n<tot_conf_quad; ++n) {
      const double *b_ord = gkyl_array_cfetch(up->conf_basis_at_ords, n);

      // Zero out quadrature values
      n_quad[n] = 0.0;
      for (int d=0; d<vdim; ++d) {
        V_drift_quad[n][d] = 0.0;
        // Store the cell average of V_drift to use if V_drift^2 > c^2 at quadrature points
        V_drift_quad_cell_avg[n][d] = V_drift_d[num_conf_basis*d]*b_ord[0];
      }
      T_over_m_quad[n] = 0.0;

      // Compute the configuration-space quadrature
      for (int k=0; k<num_conf_basis; ++k) {
        n_quad[n] += n_d[k]*b_ord[k];
        for (int d=0; d<vdim; ++d) {
          V_drift_quad[n][d] += V_drift_d[num_conf_basis*d+k]*b_ord[k];
        }
        T_over_m_quad[n] += T_over_m_d[k]*b_ord[k];
      }

      if (up->is_canonical_pb) {
        for (int k=0; k<num_conf_basis; ++k) {
          for (int j=0; j<vdim*(vdim+1)/2; ++j) {
            h_ij_inv_quad[n][j] = 0;
          }
        }
        for (int k=0; k<num_conf_basis; ++k) {
          for (int j=0; j<vdim*(vdim+1)/2; ++j) {
            h_ij_inv_quad[n][j] += h_ij_inv_d[num_conf_basis*j+k]*b_ord[k];
          }
        }
      }

      // Amplitude of the exponential.
      if ((n_quad[n] > 0.0) && (T_over_m_quad[n] > 0.0)) {
        if (up->is_relativistic) {
          expamp_quad[n] = n_quad[n]*(1.0/(4.0*GKYL_PI*T_over_m_quad[n]))*(sqrt(2*T_over_m_quad[n]/GKYL_PI));;
        }
        else {
          expamp_quad[n] = n_quad[n]/sqrt(pow(2.0*GKYL_PI*T_over_m_quad[n], vdim));
        }
      }
      else {
        expamp_quad[n] = 0.0;
      }      
    }

    // inner loop over velocity space
    gkyl_range_deflate(&vel_rng, phase_range, rem_dir, conf_iter.idx);
    gkyl_range_iter_no_split_init(&vel_iter, &vel_rng);
    while (gkyl_range_iter_next(&vel_iter)) {
      
      copy_idx_arrays(conf_range->ndim, phase_range->ndim, conf_iter.idx, vel_iter.idx, pidx);
      gkyl_rect_grid_cell_center(&up->phase_grid, pidx, xc);

      struct gkyl_range_iter qiter;
      // compute Maxwellian at phase-space quadrature nodes
      gkyl_range_iter_init(&qiter, &up->phase_qrange);
      while (gkyl_range_iter_next(&qiter)) {

        int cqidx = gkyl_range_idx(&up->conf_qrange, qiter.idx);
        int pqidx = gkyl_range_idx(&up->phase_qrange, qiter.idx);

        comp_to_phys(pdim, gkyl_array_cfetch(up->ordinates, pqidx),
          up->phase_grid.dx, xc, xmu);

        double *fq = gkyl_array_fetch(up->fun_at_ords, pqidx);
        fq[0] = f_floor;
        if (T_over_m_quad[cqidx] > 0.0) {
          if (up->is_relativistic) {
            double uu = 0.0;
            double vu = 0.0;
            double vv = 0.0;
            for (int d=0; d<vdim; ++d) {
              vv += (V_drift_quad[cqidx][d]*V_drift_quad[cqidx][d]);
              vu += (V_drift_quad[cqidx][d]*xmu[cdim+d]);
              uu += (xmu[cdim+d]*xmu[cdim+d]);
            }
            double gamma_shifted = 0.0;
            if (vv > 1.0) {
              // Check if V_drift^2 > c^2 (where c = 1.0) at quadrature points 
              // If it is, switch to just using the cell average of V_drift for
              // computing the Lorentz boost factor
              double V_drift_sq_avg = 0.0;
              for (int d=0; d<vdim; ++d) { 
                V_drift_sq_avg += (V_drift_quad_cell_avg[cqidx][d]*V_drift_quad_cell_avg[cqidx][d]);
              }
              gamma_shifted = 1.0/sqrt(1.0-V_drift_sq_avg);
            } 
            else {
              gamma_shifted = 1.0/sqrt(1.0-vv);
            }

            fq[0] += expamp_quad[cqidx]*exp( (1.0/T_over_m_quad[cqidx]) 
              - (gamma_shifted/T_over_m_quad[cqidx])*(sqrt(1+uu) - vu) );
          }
          // Assumes a (particle) hamiltonian in canocial form: g = 1/2g^{ij}w_i_w_j
          else if (up->is_canonical_pb) {
            double efact = 0.0;
            for (int d0=0; d0<vdim; ++d0) {
              for (int d1=d0; d1<vdim; ++d1) {
                int sym_tensor_index = (d0*(2*vdim - d0 + 1))/2 + (d1-d0);
                // Grab the spatial metric component, the ctx includes geometry that isn't 
                // part of the canonical set of variables, like R on the surf of a sphere
                // q_can includes the canonical variables list
                double h_ij_inv_loc = h_ij_inv_quad[cqidx][sym_tensor_index]; 
                // For off-diagnol components, we need to count these twice, due to symmetry
                int sym_fact = 2;
                if (d0 == d1){
                  sym_fact = 1;
                }
                efact += sym_fact*h_ij_inv_loc*(xmu[cdim+d0]-V_drift_quad[cqidx][d0])*(xmu[cdim+d1]-V_drift_quad[cqidx][d1]);
              }
            }
            // Accuracy of the prefactor doesn't really matter since it will 
            // be fixed by the correct routine
            fq[0] += expamp_quad[cqidx]*exp(-efact/(2.0*T_over_m_quad[cqidx]));
          }
          else {
            double efact = 0.0;        
            for (int d=0; d<vdim; ++d) {
              efact += (xmu[cdim+d]-V_drift_quad[cqidx][d])*(xmu[cdim+d]-V_drift_quad[cqidx][d]);
            }
            fq[0] += expamp_quad[cqidx]*exp(-efact/(2.0*T_over_m_quad[cqidx]));
          }
        }
      }
      // compute expansion coefficients of Maxwell-Juttner on basis
      long lidx = gkyl_range_idx(&vel_rng, vel_iter.idx);
      proj_on_basis(up, up->fun_at_ords, gkyl_array_fetch(f_lte, lidx));
    }
  }
  // Correct the density of the projected LTE distribution function through rescaling.
  // This correction is needed especially for the relativistic LTE, whose pre-factor
  // we construct through an expansion of the Bessel functions to avoid finite 
  // precision effects in such a way that we can recover arbitrary temperature 
  // relativistic LTE distributions by rescaling the distribution to the desired density.  
  gkyl_vlasov_lte_density_moment_advance(up->moments_up, phase_range, conf_range, f_lte, up->num_ratio);

  // compute number density ratio: num_ratio = n/n0
  // 0th component of moms_target is the target density
  gkyl_dg_div_op_range(up->mem, up->conf_basis, 0, up->num_ratio,
    0, moms_lte, 0, up->num_ratio, conf_range);

  // rescale distribution function
  gkyl_dg_mul_conf_phase_op_range(&up->conf_basis, &up->phase_basis,
    f_lte, up->num_ratio, f_lte, conf_range, phase_range);
}

void
gkyl_vlasov_lte_proj_on_basis_release(gkyl_vlasov_lte_proj_on_basis* up)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    gkyl_cu_free(up->p2c_qidx);
#endif
  gkyl_array_release(up->ordinates);
  gkyl_array_release(up->weights);
  gkyl_array_release(up->basis_at_ords);
  gkyl_array_release(up->conf_ordinates);
  gkyl_array_release(up->conf_weights);
  gkyl_array_release(up->conf_basis_at_ords);
  gkyl_array_release(up->fun_at_ords);

  gkyl_vlasov_lte_moments_release(up->moments_up);
  gkyl_array_release(up->num_ratio);
  gkyl_dg_bin_op_mem_release(up->mem);

  gkyl_free(up);
}
