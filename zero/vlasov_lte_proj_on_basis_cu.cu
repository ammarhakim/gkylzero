/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_const.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_vlasov_lte_proj_on_basis.h>
#include <gkyl_vlasov_lte_proj_on_basis_priv.h>
#include <gkyl_range.h>
}

__global__ static void
gkyl_vlasov_lte_proj_on_basis_advance_cu_ker(const struct gkyl_rect_grid phase_grid,
  const struct gkyl_range phase_range, const struct gkyl_range conf_range,
  const struct gkyl_array* GKYL_RESTRICT conf_basis_at_ords, 
  const struct gkyl_array* GKYL_RESTRICT phase_basis_at_ords, 
  const struct gkyl_array* GKYL_RESTRICT phase_ordinates, 
  const struct gkyl_array* GKYL_RESTRICT phase_weights, const int *p2c_qidx, bool is_relativistic, 
  bool is_canonical_pb, const struct gkyl_array* GKYL_RESTRICT h_ij_inv,  const struct gkyl_array* GKYL_RESTRICT det_h,
  const struct gkyl_array* GKYL_RESTRICT moms_lte, struct gkyl_array* GKYL_RESTRICT f_lte)
{
  double f_floor = 1.e-40;
  int pdim = phase_range.ndim, cdim = conf_range.ndim;
  int vdim = pdim-cdim;

  int num_conf_basis = conf_basis_at_ords->ncomp;
  int num_phase_basis = f_lte->ncomp;
  int tot_conf_quad = conf_basis_at_ords->size;
  int tot_phase_quad = phase_basis_at_ords->size;

  // JJ 2024/03/23: hard-coded to 3x, vdim=3, p=2 for now.
  // This hard-coding avoids issues with GPUs and dynamic memory allocation.
  double n_quad[27], V_drift_quad[27][3], T_over_m_quad[27];
  double V_drift_quad_cell_avg[27][3];
  double expamp_quad[27];
  double h_ij_inv_quad[27][6];
  double det_h_quad[27];

  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM];
  int pidx[GKYL_MAX_DIM], cidx[GKYL_MAX_CDIM];

  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < phase_range.volume; tid += blockDim.x*gridDim.x) {
    gkyl_sub_range_inv_idx(&phase_range, tid, pidx);

    // get configuration-space linear index.
    for (unsigned int k = 0; k < conf_range.ndim; k++) {
      cidx[k] = pidx[k];
    }
    long lincC = gkyl_range_idx(&conf_range, cidx);

    const double *moms_lte_d = (const double*) gkyl_array_cfetch(moms_lte, lincC);
    const double *n_d = moms_lte_d;
    const double *V_drift_d = &moms_lte_d[num_conf_basis];
    const double *T_over_m_d = &moms_lte_d[num_conf_basis*(vdim+1)];
    const double *h_ij_inv_d;
    const double *det_h_d;
    if (is_canonical_pb) {
      h_ij_inv_d = (const double*) gkyl_array_cfetch(h_ij_inv, lincC);
      det_h_d = (const double*) gkyl_array_cfetch(det_h, lincC);
    }

    // Sum over basis for given LTE moments (n, V_drift, T/m) in the stationary frame
    for (int n=0; n<tot_conf_quad; ++n) {
      const double *b_ord = (const double*) gkyl_array_cfetch(conf_basis_at_ords, n);

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

      if (is_canonical_pb) {
        det_h_quad[n] = 0.0;
        for (int k=0; k<num_conf_basis; ++k) {
          for (int j=0; j<vdim*(vdim+1)/2; ++j) {
            h_ij_inv_quad[n][j] = 0;
          }
        }
        for (int k=0; k<num_conf_basis; ++k) {
          det_h_quad[n] += det_h_d[k]*b_ord[k];
          for (int j=0; j<vdim*(vdim+1)/2; ++j) {
            h_ij_inv_quad[n][j] += h_ij_inv_d[num_conf_basis*j+k]*b_ord[k];
          }
        }
      }

      // Amplitude of the exponential.
      if ((n_quad[n] > 0.0) && (T_over_m_quad[n] > 0.0)) {
        if (is_relativistic) {
          expamp_quad[n] = n_quad[n]*(1.0/(4.0*GKYL_PI*T_over_m_quad[n]))*(sqrt(2*T_over_m_quad[n]/GKYL_PI));;
        }
        else if (is_canonical_pb) { 
          expamp_quad[n] = (1/det_h_quad[n])*n_quad[n]/sqrt(pow(2.0*GKYL_PI*T_over_m_quad[n], vdim));
        }
        else {
          expamp_quad[n] = n_quad[n]/sqrt(pow(2.0*GKYL_PI*T_over_m_quad[n], vdim));
        }
      }
      else {
        expamp_quad[n] = 0.0;
      }      
    }

    gkyl_rect_grid_cell_center(&phase_grid, pidx, xc);

    long lidx = gkyl_range_idx(&phase_range, pidx);
    double *f_lte_d = (double*) gkyl_array_fetch(f_lte, lidx);

    for (int k=0; k<num_phase_basis; ++k) {
      f_lte_d[k] = 0.0;
    }

    // compute expansion coefficients of LTE distribution function on basis
    // The following is modeled after proj_on_basis in the private header.
    const double *phase_w = (const double*) phase_weights->data;
    const double *phaseb_o = (const double*) phase_basis_at_ords->data;
  
    // compute Maxwellian at phase-space quadrature nodes
    for (int n=0; n<tot_phase_quad; ++n) {

      int cqidx = p2c_qidx[n];

      comp_to_phys(pdim, (const double*) gkyl_array_cfetch(phase_ordinates, n),
        phase_grid.dx, xc, &xmu[0]);

      double fq = f_floor;
      if (T_over_m_quad[cqidx] > 0.0) {
        if (is_relativistic) {
          double uu = 0.0;
          double vu = 0.0;
          double vv = 0.0;
          for (int d=0; d<vdim; ++d){
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

          fq += expamp_quad[cqidx]*exp( (1.0/T_over_m_quad[cqidx]) 
            - (gamma_shifted/T_over_m_quad[cqidx])*(sqrt(1+uu) - vu) );
        }
          // Assumes a (particle) hamiltonian in canocial form: g = 1/2g^{ij}w_i_w_j
          else if (is_canonical_pb) {
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
            fq += expamp_quad[cqidx]*exp(-efact/(2.0*T_over_m_quad[cqidx]));
        }
        else {
          double efact = 0.0;        
          for (int d=0; d<vdim; ++d) {
            efact += (xmu[cdim+d]-V_drift_quad[cqidx][d])*(xmu[cdim+d]-V_drift_quad[cqidx][d]);
          }
          fq += expamp_quad[cqidx]*exp(-efact/(2.0*T_over_m_quad[cqidx]));
        }
      }

      double tmp = phase_w[n]*fq;
      for (int k=0; k<num_phase_basis; ++k) {
        f_lte_d[k] += tmp*phaseb_o[k+num_phase_basis*n];
      }
    }
  }
}

void
gkyl_vlasov_lte_proj_on_basis_advance_cu(gkyl_vlasov_lte_proj_on_basis *up,
  const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
  const struct gkyl_array *moms_lte, struct gkyl_array *f_lte)
{
  int nblocks = phase_range->nblocks, nthreads = phase_range->nthreads;
  gkyl_vlasov_lte_proj_on_basis_advance_cu_ker<<<nblocks, nthreads>>>
    (up->phase_grid, *phase_range, *conf_range, up->conf_basis_at_ords->on_dev, up->basis_at_ords->on_dev,
     up->ordinates->on_dev, up->weights->on_dev, up->p2c_qidx,
     up->is_relativistic, up->is_canonical_pb, up->h_ij_inv->on_dev,  
     up->det_h->on_dev, moms_lte->on_dev, f_lte->on_dev);

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
