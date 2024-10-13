#include <cuda_runtime.h>
#include <cublas_v2.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_fpo_proj_maxwellian_pots_on_basis.h>
#include <gkyl_fpo_proj_maxwellian_pots_on_basis_priv.h>
#include <gkyl_basis.h>
#include <gkyl_const.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

#include <gkyl_rect_grid.h>
#include <gkyl_array_rio.h>
}

// Enum to choose potentials to compute
enum POT_TYPE { FPO_H, FPO_G };

static void
gkyl_parallelize_components_kernel_launch_dims(dim3* dimGrid, dim3* dimBlock, gkyl_range range, int ncomp)
{
  // Create a 2D thread grid so we launch ncomp*range.volume number of threads
  // so we can parallelize over components too
  dimBlock->y = ncomp; // ncomp *must* be less than 256
  dimGrid->y = 1;
  dimBlock->x = GKYL_DEFAULT_NUM_THREADS/ncomp;
  dimGrid->x = gkyl_int_div_up(range.volume, dimBlock->x);
}

__global__ static void
gkyl_proj_maxwellian_pots_on_basis_conf_quad_ker(struct gkyl_range conf_range,
  const struct gkyl_array *conf_basis_at_ords, const struct gkyl_array *prim_moms,
  struct gkyl_array *prim_moms_conf_quad)
{
  int vdim = 3;
  int num_conf_basis = conf_basis_at_ords->ncomp;
  int tot_conf_quad = conf_basis_at_ords->size;

  int cidx[GKYL_MAX_CDIM];

  // 2D thread grid
  // linq goes from 0 to tot_conf_quad
  long linq = threadIdx.y + blockIdx.y*blockDim.y;
  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < conf_range.volume; tid += blockDim.x*gridDim.x) {
    gkyl_sub_range_inv_idx(&conf_range, tid, cidx);

    long linc = gkyl_range_idx(&conf_range, cidx);

    const double *prim_moms_d = (const double *)gkyl_array_cfetch(prim_moms, linc);
    double *prim_moms_conf_quad_d = (double *)gkyl_array_fetch(prim_moms_conf_quad, linc);

    // Sum over basis of primitive moments (n, u_drift, T/m) at configuration space points
    const double *conf_basis_at_ord_q = (const double *) gkyl_array_cfetch(conf_basis_at_ords, linq);
    for (int k=0; k<num_conf_basis; ++k) {
      for (int d=0; d<vdim+2; ++d) {
        prim_moms_conf_quad_d[d*tot_conf_quad+linq] += prim_moms_d[num_conf_basis*d+k]*conf_basis_at_ord_q[k];
      }
    }
  }
}

__global__ static void
gkyl_proj_maxwellian_pots_on_basis_conf_nodes_ker(struct gkyl_range conf_range,
  const struct gkyl_array *conf_basis_at_nodes, const struct gkyl_array *prim_moms,
  struct gkyl_array *prim_moms_conf_nodes)
{
  int vdim = 3;
  int num_conf_basis = conf_basis_at_nodes->ncomp;

  int cidx[GKYL_MAX_CDIM];

  // 2D thread grid
  // linq goes from 0 to num_conf_basis
  long linn = threadIdx.y + blockIdx.y*blockDim.y;
  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < conf_range.volume; tid += blockDim.x*gridDim.x) {
    gkyl_sub_range_inv_idx(&conf_range, tid, cidx);

    long linc = gkyl_range_idx(&conf_range, cidx);

    const double *prim_moms_d = (const double *)gkyl_array_cfetch(prim_moms, linc);
    double *prim_moms_conf_nodes_d = (double *)gkyl_array_fetch(prim_moms_conf_nodes, linc);

    // Sum over basis of primitive moments (n, u_drift, T/m) at configuration space points
    const double *conf_basis_at_node_n = (const double *)gkyl_array_cfetch(conf_basis_at_nodes, linn);
    for (int k=0; k<num_conf_basis; ++k) {
      for (int d=0; d<vdim+2; ++d) {
        prim_moms_conf_nodes_d[d*num_conf_basis+linn] += prim_moms_d[num_conf_basis*d+k]*conf_basis_at_node_n[k];
      }
    }
  }
}

__global__ static void
gkyl_proj_maxwellian_pots_on_basis_phase_quad_ker(POT_TYPE pots_to_calc, struct gkyl_rect_grid phase_grid,
  struct gkyl_range phase_range, struct gkyl_range conf_range, const struct gkyl_array *ordinates,
  const struct gkyl_array *conf_basis_at_ords,
  const struct gkyl_array *prim_moms_conf_quad, int *p2c_qidx, struct gkyl_array *pot_phase_quad)
{
  int pdim = phase_range.ndim, cdim = conf_range.ndim;
  int vdim = 3;
  int tot_conf_quad = conf_basis_at_ords->size;

  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM];
  int pidx[GKYL_MAX_DIM], cidx[GKYL_MAX_CDIM];

  // 2D thread grid
  // linq goes from 0 to tot_phase_quad
  long linq = threadIdx.y + blockIdx.y*blockDim.y;
  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < phase_range.volume; tid += blockDim.x*gridDim.x) {
    // Phase and configuration space indices
    gkyl_sub_range_inv_idx(&phase_range, tid, pidx);
    for (int i=0; i<cdim; ++i) cidx[i] = pidx[i];

    gkyl_rect_grid_cell_center(&phase_grid, pidx, xc);
    long linc = gkyl_range_idx(&conf_range, cidx);
    long linp = gkyl_range_idx(&phase_range, pidx);

    // Fetch primitive moments evaluated at quadrature points
    const double *prim_moms_conf_quad_d = (const double *)gkyl_array_cfetch(prim_moms_conf_quad, linc);
    const double *den_d = &prim_moms_conf_quad_d[0];
    const double *u_drift_d = &prim_moms_conf_quad_d[tot_conf_quad];
    const double *vtsq_d = &prim_moms_conf_quad_d[tot_conf_quad*(vdim+1)];

    double *pot_q = (double *)gkyl_array_fetch(pot_phase_quad, linp);

    comp_to_phys(pdim, (const double*)gkyl_array_cfetch(ordinates, linq),
      phase_grid.dx, xc, &xmu[0]);
    int cqidx = p2c_qidx[linq];

    double rel_speedsq_d = 0.0;
    for (int d=0; d<vdim; ++d) {
      rel_speedsq_d += pow(xmu[cdim+d]-u_drift_d[tot_conf_quad*d+cqidx],2);
    }
    double rel_speed_d = sqrt(rel_speedsq_d);

    if (pots_to_calc == FPO_H)
      pot_q[linq] += eval_fpo_h(den_d[cqidx], rel_speed_d, vtsq_d[cqidx]);
    else if (pots_to_calc == FPO_G)
      pot_q[linq] += eval_fpo_g(den_d[cqidx], rel_speed_d, vtsq_d[cqidx]);
  }
}

__global__ static void
gkyl_proj_maxwellian_pots_on_basis_surf_quad_ker(POT_TYPE pots_to_calc, struct gkyl_rect_grid phase_grid,
  struct gkyl_range phase_range, struct gkyl_range conf_range, int tot_conf_quad,
  const struct gkyl_array *surf_basis_at_ords, const struct gkyl_array *surf_ordinates,
  const struct gkyl_array *prim_moms_conf_quad, int *surf2c_qidx,
  struct gkyl_array *pot_surf_quad, struct gkyl_array *pot_deriv_surf_quad)
{
  int pdim = phase_range.ndim, cdim = conf_range.ndim;
  int vdim = 3;
  int tot_surf_quad = surf_basis_at_ords->size;

  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM];
  int pidx[GKYL_MAX_DIM], cidx[GKYL_MAX_CDIM];

  // 2D thread grid
  // linq goes from 0 to tot_surf_quad
  long linq = threadIdx.y + blockIdx.y*blockDim.y;
  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < phase_range.volume; tid += blockDim.x*gridDim.x) {
    // Phase and configuration space indices
    gkyl_sub_range_inv_idx(&phase_range, tid, pidx);
    for (int i=0; i<cdim; ++i) cidx[i] = pidx[i];

    gkyl_rect_grid_cell_center(&phase_grid, pidx, xc);
    long linc = gkyl_range_idx(&conf_range, cidx);
    long linp = gkyl_range_idx(&phase_range, pidx);

    // Fetch primitive moments evaluated at quadrature points
    const double *prim_moms_conf_quad_d = (const double *)gkyl_array_cfetch(prim_moms_conf_quad, linc);
    const double *den_d = &prim_moms_conf_quad_d[0];
    const double *u_drift_d = &prim_moms_conf_quad_d[tot_conf_quad];
    const double *vtsq_d = &prim_moms_conf_quad_d[tot_conf_quad*(vdim+1)];

    double *pot_q = (double *)gkyl_array_fetch(pot_surf_quad, linp);
    double *pot_deriv_q = (double *)gkyl_array_fetch(pot_deriv_surf_quad, linp);

    // Iterate over velocity directions for H, dH/dv, G, d2G/dv2 at each dir1 domain boundary
    // Derivatives dH/dv and d2G/dv2 are normal to surface
    for (int d1=0; d1<vdim; ++d1) {
      int dir1 = d1 + cdim;
      bool is_edge_in_dir1 = 0;
      is_edge_in_dir1 = is_edge_in_dir1 || (pidx[dir1] == phase_range.lower[dir1]);
      is_edge_in_dir1 = is_edge_in_dir1 || (pidx[dir1] == phase_range.upper[dir1]);

      const double* surf_ord = (const double*)gkyl_array_cfetch(surf_ordinates, linq);
      int cqidx = surf2c_qidx[linq];

      if (is_edge_in_dir1) {
        // Velocity value at dir1 boundary
        double vbound = pidx[dir1] == phase_range.lower[dir1] ?
          phase_grid.lower[dir1] : phase_grid.upper[dir1];

        surf_comp_to_phys(dir1, pdim, surf_ord, phase_grid.dx, xc, xmu);
        xmu[dir1] = vbound;

        double rel_speedsq_d = 0.0;
        double rel_vel_in_dir1_d = xmu[dir1] - u_drift_d[tot_conf_quad*d1+cqidx];
        for (int d=0; d<vdim; ++d) {
          rel_speedsq_d += pow(xmu[cdim+d]-u_drift_d[tot_conf_quad*d+cqidx],2);
        }
        double rel_speed_d = sqrt(rel_speedsq_d);

        if (pots_to_calc == FPO_H) {
          pot_q[d1*tot_surf_quad + linq] += eval_fpo_h(den_d[cqidx], rel_speed_d, vtsq_d[cqidx]);
          pot_deriv_q[d1*tot_surf_quad + linq] += eval_fpo_dhdv(den_d[cqidx], rel_vel_in_dir1_d, vtsq_d[cqidx], rel_speed_d);
        }
        else if (pots_to_calc == FPO_G) {
          pot_q[d1*tot_surf_quad + linq] += eval_fpo_g(den_d[cqidx], rel_speed_d, vtsq_d[cqidx]);
          pot_deriv_q[d1*tot_surf_quad + linq] += eval_fpo_d2gdv2(den_d[cqidx], rel_vel_in_dir1_d, vtsq_d[cqidx], rel_speed_d);
        }
      }
    }
  }
}

__global__ static void
gkyl_proj_dgdv_on_basis_surf_nodal_ker(struct gkyl_rect_grid phase_grid,
  struct gkyl_range phase_range, struct gkyl_range conf_range, const struct gkyl_array* surf_nodes,
  const struct gkyl_basis surf_basis, const struct gkyl_basis conf_basis,
  const struct gkyl_array *prim_moms_conf_nodes, int *surf2c_nidx,
  struct gkyl_array *fpo_dgdv_at_surf_nodes)
{
  int pdim = phase_range.ndim, cdim = conf_range.ndim;
  int vdim = 3;
  int num_surf_basis = surf_basis.num_basis;
  int num_conf_basis = conf_basis.num_basis;

  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM], xmu_conf[GKYL_MAX_DIM];
  int pidx[GKYL_MAX_DIM], cidx[GKYL_MAX_CDIM];

  // 2D thread grid
  // linn goes from 0 to num_surf_basis
  long linn = threadIdx.y + blockIdx.y*blockDim.y;
  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < phase_range.volume; tid += blockDim.x*gridDim.x) {
    // Phase and configuration space indices
    gkyl_sub_range_inv_idx(&phase_range, tid, pidx);
    for (int i=0; i<cdim; ++i) cidx[i] = pidx[i];

    gkyl_rect_grid_cell_center(&phase_grid, pidx, xc);
    long linc = gkyl_range_idx(&conf_range, cidx);
    long linp = gkyl_range_idx(&phase_range, pidx);

    const double *node = (const double *)gkyl_array_cfetch(surf_nodes, linn);
    comp_to_phys(cdim, node, phase_grid.dx, xc, &xmu_conf[0]);
    int cnidx = surf2c_nidx[linn];

    // Fetch primitive moments evaluated at nodes
    const double *prim_moms_conf_nodes_d = (const double *)gkyl_array_cfetch(prim_moms_conf_nodes, linc);
    const double *den_d = &prim_moms_conf_nodes_d[0];
    const double *u_drift_d = &prim_moms_conf_nodes_d[num_conf_basis];
    const double *vtsq_d = &prim_moms_conf_nodes_d[num_conf_basis*(vdim+1)];

    double *fpo_dgdv_at_surf_nodes_n = (double *)gkyl_array_fetch(fpo_dgdv_at_surf_nodes, linp);

    for (int d1=0; d1<vdim; ++d1) {
      int dir1 = d1 + cdim;
      bool is_edge_in_dir1 = 0;
      is_edge_in_dir1 = is_edge_in_dir1 || (pidx[dir1] == phase_range.lower[dir1]);
      is_edge_in_dir1 = is_edge_in_dir1 || (pidx[dir1] == phase_range.upper[dir1]);

      if (is_edge_in_dir1) {
        for (int d2=0; d2<vdim; ++d2) {
          if (d1 == d2) continue;
          int dir2 = d2+cdim;

          // Velocity value at dir boundary
          double vbound = pidx[dir1] == phase_range.lower[dir1] ?
            phase_grid.lower[dir1] : phase_grid.upper[dir1];

          surf_comp_to_phys(dir1, pdim, node, phase_grid.dx, xc, xmu);
          xmu[dir1] = vbound;

          double rel_speedsq_d = 0.0;
          double rel_vel_in_dir2_d = xmu[dir2] - u_drift_d[num_conf_basis*d2+cnidx];
          for (int d=0; d<vdim; ++d) {
            rel_speedsq_d += pow(xmu[cdim+d]-u_drift_d[num_conf_basis*d+cnidx],2);
          }

          int dgdv_index = (d1 < d2) ? (vdim-1)*d1+d2-1 : (vdim-1)*d1+d2;
          fpo_dgdv_at_surf_nodes_n[dgdv_index*num_surf_basis+linn] += eval_fpo_dgdv(den_d[cnidx],
            rel_vel_in_dir2_d, vtsq_d[cnidx], rel_speedsq_d);
        }
      }
    }
  }
}

__global__ static void
fpo_dgdv_nod2mod_ker(const struct gkyl_range phase_range, struct gkyl_basis *basis,
  const struct gkyl_array *fpo_dgdv_at_surf_nodes, struct gkyl_array *fpo_dgdv_surf)
{
  int pdim = phase_range.ndim, vdim = 3;
  int cdim = pdim - vdim;
  int pidx[GKYL_MAX_DIM];

  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < phase_range.volume; tid += blockDim.x*gridDim.x) {
    gkyl_sub_range_inv_idx(&phase_range, tid, pidx);
    long linp = gkyl_range_idx(&phase_range, pidx);

    const double *nodal_d = (const double*)gkyl_array_cfetch(fpo_dgdv_at_surf_nodes, linp);
    double *fpo_dgdv_surf_d = (double *)gkyl_array_fetch(fpo_dgdv_surf, linp);

    int num_basis = basis->num_basis;
    for (int d1=0; d1<vdim; ++d1) {
      int dir1 = d1 + cdim;
      bool is_edge_in_dir1 = 0;
      is_edge_in_dir1 = is_edge_in_dir1 || (pidx[dir1] == phase_range.lower[dir1]);
      is_edge_in_dir1 = is_edge_in_dir1 || (pidx[dir1] == phase_range.upper[dir1]);

      if (is_edge_in_dir1) {
        for (int d2=0; d2<vdim; ++d2) {
          if (d1 == d2) continue;

          int dgdv_index = (d1 < d2) ? (vdim-1)*d1+d2-1 : (vdim-1)*d1+d2;
          basis->nodal_to_modal(&nodal_d[dgdv_index*num_basis], &fpo_dgdv_surf_d[dgdv_index*num_basis]);
        }
      }
    }
  }
}

void
gkyl_proj_maxwellian_pots_on_basis_advance_cu(const gkyl_proj_maxwellian_pots_on_basis *up,
  const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
  const struct gkyl_array* prim_moms,
  struct gkyl_array *fpo_h, struct gkyl_array *fpo_g,
  struct gkyl_array *fpo_h_surf, struct gkyl_array *fpo_g_surf,
  struct gkyl_array *fpo_dhdv_surf, struct gkyl_array *fpo_dgdv_surf,
  struct gkyl_array *fpo_d2gdv2_surf)
{
  // Zero out input arrays
  gkyl_array_clear(fpo_h, 0.0);
  gkyl_array_clear(fpo_g, 0.0);
  gkyl_array_clear(fpo_h_surf, 0.0);
  gkyl_array_clear(fpo_g_surf, 0.0);
  gkyl_array_clear(fpo_dhdv_surf, 0.0);
  gkyl_array_clear(fpo_dgdv_surf, 0.0);
  gkyl_array_clear(fpo_d2gdv2_surf, 0.0);

  // Evaluate primitive moments at configuration space quadrature points and nodes
  dim3 dimGrid_conf, dimBlock_conf;
  int tot_conf_quad = up->tot_conf_quad;
  gkyl_array_clear(up->prim_moms_conf_quad, 0.0);
  gkyl_parallelize_components_kernel_launch_dims(&dimGrid_conf, &dimBlock_conf, *conf_range, tot_conf_quad);
  gkyl_proj_maxwellian_pots_on_basis_conf_quad_ker<<<dimGrid_conf, dimBlock_conf>>>(*conf_range,
    up->conf_basis_at_ords->on_dev, prim_moms->on_dev, up->prim_moms_conf_quad->on_dev);

  gkyl_array_clear(up->prim_moms_conf_nodes, 0.0);
  gkyl_parallelize_components_kernel_launch_dims(&dimGrid_conf, &dimBlock_conf, *conf_range, up->num_conf_basis);
  gkyl_proj_maxwellian_pots_on_basis_conf_nodes_ker<<<dimGrid_conf, dimBlock_conf>>>(*conf_range,
    up->conf_basis_at_nodes->on_dev, prim_moms->on_dev, up->prim_moms_conf_nodes->on_dev);

  // Phase space quadrature quantities. First compute and project H, then the same for G.
  // Computed individually like this to avoid having to allocate more arrays
  dim3 dimGrid, dimBlock;
  int tot_quad = up->tot_quad;
  gkyl_array_clear(up->pot_phase_quad, 0.0);
  gkyl_parallelize_components_kernel_launch_dims(&dimGrid, &dimBlock, *phase_range, tot_quad);
  gkyl_proj_maxwellian_pots_on_basis_phase_quad_ker<<<dimGrid, dimBlock>>>(FPO_H, up->grid,
    *phase_range, *conf_range, up->ordinates->on_dev, up->conf_basis_at_ords->on_dev,
    up->prim_moms_conf_quad->on_dev, up->p2c_qidx, up->pot_phase_quad->on_dev);

  // Matrix multiplication for nodal to modal conversion
  gkyl_mat_mm_array(up->phase_quad_nodal_to_modal_mem, up->pot_phase_quad, fpo_h);

  // Compute G
  gkyl_array_clear(up->pot_phase_quad, 0.0);
  gkyl_proj_maxwellian_pots_on_basis_phase_quad_ker<<<dimGrid, dimBlock>>>(FPO_G, up->grid,
    *phase_range, *conf_range, up->ordinates->on_dev, up->conf_basis_at_ords->on_dev,
    up->prim_moms_conf_quad->on_dev, up->p2c_qidx, up->pot_phase_quad->on_dev);
  gkyl_mat_mm_array(up->phase_quad_nodal_to_modal_mem, up->pot_phase_quad, fpo_g);

  // Phase space surface quadrature quantities.
  // First call is to compute H, dH/dv surface expansions
  dim3 dimGrid_surf, dimBlock_surf;
  int tot_surf_quad = up->tot_surf_quad;
  gkyl_array_clear(up->pot_surf_quad, 0.0);
  gkyl_array_clear(up->pot_deriv_surf_quad, 0.0);
  gkyl_parallelize_components_kernel_launch_dims(&dimGrid_surf, &dimBlock_surf, *phase_range, tot_surf_quad);
  gkyl_proj_maxwellian_pots_on_basis_surf_quad_ker<<<dimGrid_surf, dimBlock_surf>>>(FPO_H,
    up->grid, *phase_range, *conf_range, up->tot_conf_quad, up->surf_basis_at_ords->on_dev,
    up->surf_ordinates->on_dev, up->prim_moms_conf_quad->on_dev, up->surf2c_qidx,
    up->pot_surf_quad->on_dev, up->pot_deriv_surf_quad->on_dev);

  gkyl_mat_mm_array(up->surf_quad_nodal_to_modal_mem, up->pot_surf_quad, up->sol_pot_surf_modal);
  gkyl_array_clear(fpo_h_surf, 0.0);
  gkyl_array_accumulate(fpo_h_surf, 1.0, up->sol_pot_surf_modal);
  gkyl_mat_mm_array(up->surf_quad_nodal_to_modal_mem, up->pot_deriv_surf_quad, up->sol_pot_surf_modal);
  gkyl_array_accumulate(fpo_dhdv_surf, 1.0, up->sol_pot_surf_modal);

  // Second call to phase space surface quadrature is for G, d2G/dv2 surface expansions
  gkyl_array_clear(up->pot_surf_quad, 0.0);
  gkyl_array_clear(up->pot_deriv_surf_quad, 0.0);
  gkyl_parallelize_components_kernel_launch_dims(&dimGrid_surf, &dimBlock_surf, *phase_range, tot_surf_quad);
  gkyl_proj_maxwellian_pots_on_basis_surf_quad_ker<<<dimGrid_surf, dimBlock_surf>>>(FPO_G,
    up->grid, *phase_range, *conf_range, up->tot_conf_quad, up->surf_basis_at_ords->on_dev,
    up->surf_ordinates->on_dev, up->prim_moms_conf_quad->on_dev, up->surf2c_qidx,
    up->pot_surf_quad->on_dev, up->pot_deriv_surf_quad->on_dev);

  gkyl_mat_mm_array(up->surf_quad_nodal_to_modal_mem, up->pot_surf_quad, up->sol_pot_surf_modal);
  gkyl_array_accumulate(fpo_g_surf, 1.0, up->sol_pot_surf_modal);
  gkyl_mat_mm_array(up->surf_quad_nodal_to_modal_mem, up->pot_deriv_surf_quad, up->sol_pot_surf_modal);
  gkyl_array_accumulate(fpo_d2gdv2_surf, 1.0, up->sol_pot_surf_modal);

  // Phase space nodal evaluation of dG/dv
  gkyl_array_clear(up->fpo_dgdv_at_surf_nodes, 0.0);
  dim3 dimGrid_surf_nodal, dimBlock_surf_nodal;
  gkyl_parallelize_components_kernel_launch_dims(&dimGrid_surf_nodal, &dimBlock_surf_nodal, *phase_range, up->num_surf_basis);
  gkyl_proj_dgdv_on_basis_surf_nodal_ker<<<dimGrid_surf_nodal, dimBlock_surf_nodal>>>(up->grid,
    *phase_range, *conf_range, up->surf_nodes->on_dev, up->surf_basis, up->conf_basis,
    up->prim_moms_conf_nodes->on_dev, up->surf2c_nidx,
    up->fpo_dgdv_at_surf_nodes->on_dev);

  fpo_dgdv_nod2mod_ker<<<dimGrid_surf_nodal, dimBlock_surf_nodal>>>(*phase_range,
    up->surf_basis_dev, up->fpo_dgdv_at_surf_nodes->on_dev, fpo_dgdv_surf->on_dev);
}

