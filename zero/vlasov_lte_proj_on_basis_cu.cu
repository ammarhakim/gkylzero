/* -*- c++ -*- */
#include <cuda_runtime.h>
#include <cublas_v2.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_const.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_vlasov_lte_proj_on_basis.h>
#include <gkyl_vlasov_lte_proj_on_basis_priv.h>
#include <gkyl_range.h>

#include <gkyl_mat.h>
#include <gkyl_mat_priv.h>
}

static void
gkyl_parallelize_components_kernel_launch_dims(dim3* dimGrid, dim3* dimBlock, gkyl_range range, int ncomp)
{
  // Create a 2D thread grid so we launch ncomp*range.volume number of threads and can parallelize over components too
  dimBlock->y = GKYL_MIN2(ncomp, GKYL_DEFAULT_NUM_THREADS);
  dimGrid->y = gkyl_int_div_up(ncomp, dimBlock->y);
  dimBlock->x = GKYL_DEFAULT_NUM_THREADS/ncomp;
  dimGrid->x = gkyl_int_div_up(range.volume, dimBlock->x);
}

__global__ static void
gkyl_vlasov_lte_proj_on_basis_advance_cu_ker(const struct gkyl_rect_grid phase_grid,
  const struct gkyl_range phase_range, const struct gkyl_range conf_range,
  const struct gkyl_array* GKYL_RESTRICT conf_basis_at_ords, 
  struct gkyl_array* GKYL_RESTRICT n_quad_d, 
  struct gkyl_array* GKYL_RESTRICT V_drift_quad_d, 
  struct gkyl_array* GKYL_RESTRICT T_over_m_quad_d, 
  struct gkyl_array* GKYL_RESTRICT V_drift_quad_cell_avg_d, 
  struct gkyl_array* GKYL_RESTRICT h_ij_inv_quad_d, 
  struct gkyl_array* GKYL_RESTRICT det_h_quad_d, const int *p2c_qidx, bool is_relativistic, 
  bool is_canonical_pb, const struct gkyl_array* GKYL_RESTRICT h_ij_inv,  const struct gkyl_array* GKYL_RESTRICT det_h,
  const struct gkyl_array* GKYL_RESTRICT moms_lte, const struct gkyl_basis* conf_basis_on_dev,
   struct gkyl_array* GKYL_RESTRICT f_lte_at_nodes, int num_conf_basis, int tot_conf_quad)
{
  int pdim = phase_range.ndim, cdim = conf_range.ndim;
  int vdim = pdim-cdim;

  int cidx[GKYL_MAX_CDIM];

  // 2D thread grid
  // linc2 = c where c is the component index (from 0 to tot_conf_quad)
  long linc2 = threadIdx.y + blockIdx.y*blockDim.y;
  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < conf_range.volume; tid += blockDim.x*gridDim.x) {
    gkyl_sub_range_inv_idx(&conf_range, tid, cidx);

    long lincC = gkyl_range_idx(&conf_range, cidx);

    double *n_quad = (double*) gkyl_array_fetch(n_quad_d, lincC);
    double *V_drift_quad = (double*) gkyl_array_fetch(V_drift_quad_d, lincC);
    double *T_over_m_quad = (double*) gkyl_array_fetch(T_over_m_quad_d, lincC);
    double *V_drift_quad_cell_avg = (double*) gkyl_array_fetch(V_drift_quad_cell_avg_d, lincC);
    double *h_ij_inv_quad = (double*) gkyl_array_fetch(h_ij_inv_quad_d, lincC);
    double *det_h_quad = (double*) gkyl_array_fetch(det_h_quad_d, lincC);

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
    int n = linc2;
    const double *b_ord = (const double*) gkyl_array_cfetch(conf_basis_at_ords, n);

    // Compute the configuration-space quadrature
    conf_basis_on_dev->modal_to_quad_nodal(n_d, n_quad, n);
    for (int d=0; d<vdim; ++d) {
      conf_basis_on_dev->modal_to_quad_nodal(&V_drift_d[num_conf_basis*d], &V_drift_quad[tot_conf_quad*d], n);
      V_drift_quad_cell_avg[tot_conf_quad*d + n] = V_drift_d[num_conf_basis*d]*b_ord[0];
    }
    // Update T at nodal points
    conf_basis_on_dev->modal_to_quad_nodal(T_over_m_d, T_over_m_quad, n);

    if (is_canonical_pb) {
      conf_basis_on_dev->modal_to_quad_nodal(det_h_d, det_h_quad, n);
      for (int j=0; j<vdim*(vdim+1)/2; ++j) {
        conf_basis_on_dev->modal_to_quad_nodal(&h_ij_inv_d[num_conf_basis*j], &h_ij_inv_quad[tot_conf_quad*j], n);
      }
    }     
  }
}

__global__ static void
gkyl_vlasov_lte_proj_on_basis_accumulate_cu_ker(const struct gkyl_rect_grid phase_grid,
  const struct gkyl_range phase_range, const struct gkyl_range conf_range,
  const struct gkyl_array* GKYL_RESTRICT phase_ordinates, 
  struct gkyl_array* GKYL_RESTRICT n_quad_d, 
  struct gkyl_array* GKYL_RESTRICT V_drift_quad_d, 
  struct gkyl_array* GKYL_RESTRICT T_over_m_quad_d, 
  struct gkyl_array* GKYL_RESTRICT V_drift_quad_cell_avg_d, 
  struct gkyl_array* GKYL_RESTRICT h_ij_inv_quad_d, 
  struct gkyl_array* GKYL_RESTRICT det_h_quad_d,  const int *p2c_qidx, bool is_relativistic, 
  bool is_canonical_pb, const struct gkyl_array* GKYL_RESTRICT h_ij_inv,  const struct gkyl_array* GKYL_RESTRICT det_h,
  const struct gkyl_array* GKYL_RESTRICT moms_lte, struct gkyl_array* GKYL_RESTRICT f_lte_at_nodes, int num_conf_basis, int tot_conf_quad)
{
  double f_floor = 1.e-40;
  int pdim = phase_range.ndim, cdim = conf_range.ndim;
  int vdim = pdim-cdim;

  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM];
  int pidx[GKYL_MAX_DIM], cidx[GKYL_MAX_CDIM];

  // 2D thread grid
  // linc2 = c where c is the component index (from 0 to tot_phase_quad)
  long linc2 = threadIdx.y + blockIdx.y*blockDim.y;
  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < phase_range.volume; tid += blockDim.x*gridDim.x) {
    gkyl_sub_range_inv_idx(&phase_range, tid, pidx);

    // get configuration-space linear index.
    for (unsigned int k = 0; k < conf_range.ndim; k++) {
      cidx[k] = pidx[k];
    }
    long lincC = gkyl_range_idx(&conf_range, cidx);

    double *n_quad = (double*) gkyl_array_fetch(n_quad_d, lincC);
    double *V_drift_quad = (double*) gkyl_array_fetch(V_drift_quad_d, lincC);
    double *T_over_m_quad = (double*) gkyl_array_fetch(T_over_m_quad_d, lincC);
    double *V_drift_quad_cell_avg = (double*) gkyl_array_fetch(V_drift_quad_cell_avg_d, lincC);
    double *h_ij_inv_quad = (double*) gkyl_array_fetch(h_ij_inv_quad_d, lincC);
    double *det_h_quad = (double*) gkyl_array_fetch(det_h_quad_d, lincC);

    gkyl_rect_grid_cell_center(&phase_grid, pidx, xc);

    long lidx = gkyl_range_idx(&phase_range, pidx);

    // Select for a phase space index fq
    double *fq = (double*) gkyl_array_fetch(f_lte_at_nodes, lidx);
  
    int n = linc2;

    int cqidx = p2c_qidx[n];

    comp_to_phys(pdim, (const double*) gkyl_array_cfetch(phase_ordinates, n),
      phase_grid.dx, xc, &xmu[0]);

    if (T_over_m_quad[cqidx] > 0.0) {
      if (is_relativistic) {
        double uu = 0.0;
        double vu = 0.0;
        double vv = 0.0;
        for (int d=0; d<vdim; ++d){
          vv += (V_drift_quad[tot_conf_quad*d + cqidx]*V_drift_quad[tot_conf_quad*d + cqidx]);
          vu += (V_drift_quad[tot_conf_quad*d + cqidx]*xmu[cdim+d]);
          uu += (xmu[cdim+d]*xmu[cdim+d]);
        }
        double gamma_shifted = 0.0;
        if (vv > 1.0) {
          // Check if V_drift^2 > c^2 (where c = 1.0) at quadrature points 
          // If it is, switch to just using the cell average of V_drift for
          // computing the Lorentz boost factor
          double V_drift_sq_avg = 0.0;
          for (int d=0; d<vdim; ++d) { 
            V_drift_sq_avg += (V_drift_quad_cell_avg[tot_conf_quad*d + cqidx]*V_drift_quad_cell_avg[tot_conf_quad*d + cqidx]);
          }
          gamma_shifted = 1.0/sqrt(1.0-V_drift_sq_avg);
        } 
        else {
          gamma_shifted = 1.0/sqrt(1.0-vv);
        }

        fq[n] = f_floor + exp( (1.0/T_over_m_quad[cqidx]) 
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
              double h_ij_inv_loc = h_ij_inv_quad[tot_conf_quad*sym_tensor_index + cqidx]; 
              // For off-diagnol components, we need to count these twice, due to symmetry
              int sym_fact = 2;
              if (d0 == d1){
                sym_fact = 1;
              }
              efact += sym_fact*h_ij_inv_loc*(xmu[cdim+d0]-V_drift_quad[tot_conf_quad*d0 + cqidx])*(xmu[cdim+d1]-V_drift_quad[tot_conf_quad*d1 + cqidx]);
            }
          }
          // Accuracy of the prefactor doesn't really matter since it will 
          // be fixed by the correct routine
          fq[n] = f_floor + exp(-efact/(2.0*T_over_m_quad[cqidx]));
      }
      else {
        double efact = 0.0;        
        for (int d=0; d<vdim; ++d) {
          efact += (xmu[cdim+d]-V_drift_quad[tot_conf_quad*d + cqidx])*(xmu[cdim+d]-V_drift_quad[tot_conf_quad*d + cqidx]);
        }
        fq[n] = f_floor + exp(-efact/(2.0*T_over_m_quad[cqidx]));
      }
    }
  }
}

void
gkyl_vlasov_lte_proj_on_basis_advance_cu(gkyl_vlasov_lte_proj_on_basis *up,
  const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
  const struct gkyl_array *moms_lte, struct gkyl_array *f_lte)
{
  
  dim3 dimGrid, dimBlock;
  int tot_phase_quad = up->basis_at_ords->size;
  int tot_conf_quad = up->conf_basis_at_ords->size;
  int num_conf_basis = up->num_conf_basis;
  //gkyl_parallelize_components_kernel_launch_dims(&dimGrid, &dimBlock, *phase_range, tot_phase_quad);
  gkyl_parallelize_components_kernel_launch_dims(&dimGrid, &dimBlock, *conf_range, tot_conf_quad);
  gkyl_vlasov_lte_proj_on_basis_advance_cu_ker<<<dimGrid, dimBlock>>>
    (up->phase_grid, *phase_range, *conf_range, up->conf_basis_at_ords->on_dev, 
     up->n_quad->on_dev, 
     up->V_drift_quad->on_dev, 
     up->T_over_m_quad->on_dev, 
     up->V_drift_quad_cell_avg->on_dev, 
     up->h_ij_inv_quad->on_dev, 
     up->det_h_quad->on_dev, up->p2c_qidx,
     up->is_relativistic, up->is_canonical_pb, 
     up->is_canonical_pb ? up->h_ij_inv->on_dev : 0, 
     up->is_canonical_pb ? up->det_h->on_dev : 0, 
     moms_lte->on_dev, up->conf_basis_on_dev, up->f_lte_at_nodes->on_dev, num_conf_basis, tot_conf_quad);

  //int nblocks = phase_range->nblocks, nthreads = phase_range->nthreads;
  gkyl_parallelize_components_kernel_launch_dims(&dimGrid, &dimBlock, *phase_range, tot_phase_quad);
  gkyl_vlasov_lte_proj_on_basis_accumulate_cu_ker<<<dimGrid, dimBlock>>>
    (up->phase_grid, *phase_range, *conf_range,
     up->ordinates->on_dev,
     up->n_quad->on_dev, 
     up->V_drift_quad->on_dev, 
     up->T_over_m_quad->on_dev, 
     up->V_drift_quad_cell_avg->on_dev, 
     up->h_ij_inv_quad->on_dev, 
     up->det_h_quad->on_dev, up->p2c_qidx,
     up->is_relativistic, up->is_canonical_pb, 
     up->is_canonical_pb ? up->h_ij_inv->on_dev : 0, 
     up->is_canonical_pb ? up->det_h->on_dev : 0, 
     moms_lte->on_dev, up->f_lte_at_nodes->on_dev, num_conf_basis, tot_conf_quad);

  // Call cublas to do the nodal to modal conversion
  cu_mat_mm_array(up->cuh, up->phase_nodal_to_modal_mem, up->f_lte_at_nodes, f_lte);

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
