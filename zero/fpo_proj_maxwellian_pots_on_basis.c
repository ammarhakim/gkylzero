#include <math.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_const.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_fpo_proj_maxwellian_pots_on_basis.h>
#include <gkyl_fpo_proj_maxwellian_pots_on_basis_priv.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_mat.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

static inline void
copy_idx_arrays(int cdim, int pdim, const int *cidx, const int *vidx, int *out)
{
  for (int i=0; i<cdim; ++i)
    out[i] = cidx[i];
  for (int i=cdim; i<pdim; ++i)
    out[i] = vidx[i-cdim];
}

struct gkyl_proj_maxwellian_pots_on_basis* 
gkyl_proj_maxwellian_pots_on_basis_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range,
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, int num_quad, bool use_gpu) 
{
  gkyl_proj_maxwellian_pots_on_basis *up = gkyl_malloc(sizeof(gkyl_proj_maxwellian_pots_on_basis));
  up->grid = *grid;
  up->cdim = conf_basis->ndim;
  up->pdim = phase_basis->ndim;
  up->num_quad = num_quad;

  up->conf_basis = *conf_basis;

  up->use_gpu = use_gpu;

  if (phase_basis->poly_order == 1) {  
    gkyl_cart_modal_hybrid(&up->surf_basis, up->cdim, up->pdim-up->cdim-1);
  } else {
    gkyl_cart_modal_serendip(&up->surf_basis, up->pdim-1, phase_basis->poly_order);
  }

  up->num_conf_basis = conf_basis->num_basis;
  up->num_phase_basis = phase_basis->num_basis;
  up->num_surf_basis = up->surf_basis.num_basis;
  
  // Quadrature points for phase space expansion
  up->tot_quad = init_quad_values(up->cdim, phase_basis, num_quad,
    &up->ordinates, &up->weights, &up->basis_at_ords, false);

  // Initialize quadrature for surface expansion
  up->tot_surf_quad = init_quad_values(up->cdim, &up->surf_basis, num_quad,
    &up->surf_ordinates, &up->surf_weights, &up->surf_basis_at_ords, false);

  // Hybrid basis support: uses p=2 in velocity space
  int vdim = up->pdim-up->cdim;
  int num_quad_v = num_quad;
  bool is_vdim_p2[] = {false, false, false};  // 3 is the max vdim.
  if (phase_basis->b_type == GKYL_BASIS_MODAL_HYBRID) {
    num_quad_v = num_quad+1;
    // Maxwellian potentials are always 3V
    is_vdim_p2[0] = true, is_vdim_p2[1] = true, is_vdim_p2[2] = true;
  }

  up->phase_qrange = get_qrange(up->cdim, up->pdim, num_quad, num_quad_v, is_vdim_p2);
  up->surf_qrange = get_qrange(up->cdim, up->pdim-1, num_quad, num_quad_v, is_vdim_p2);

  // Nodes for nodal surface expansions
  struct gkyl_array *surf_nodes_ho = gkyl_array_new(GKYL_DOUBLE, grid->ndim-1, up->surf_basis.num_basis);
  up->surf_basis.node_list(gkyl_array_fetch(surf_nodes_ho, 0));

  if (!up->use_gpu) {
    up->fpo_h_at_ords = gkyl_array_new(GKYL_DOUBLE, 1, up->tot_quad);
    up->fpo_g_at_ords = gkyl_array_new(GKYL_DOUBLE, 1, up->tot_quad);
    up->fpo_dhdv_at_ords = gkyl_array_new(GKYL_DOUBLE, 1, up->tot_quad);
    up->fpo_d2gdv2_at_ords = gkyl_array_new(GKYL_DOUBLE, 1, up->tot_quad);
  
    up->fpo_h_at_surf_ords = gkyl_array_new(GKYL_DOUBLE, 1, up->tot_surf_quad);
    up->fpo_g_at_surf_ords = gkyl_array_new(GKYL_DOUBLE, 1, up->tot_surf_quad);
    up->fpo_dhdv_at_surf_ords = gkyl_array_new(GKYL_DOUBLE, 1, up->tot_surf_quad);
    up->fpo_d2gdv2_at_surf_ords = gkyl_array_new(GKYL_DOUBLE, 1, up->tot_surf_quad);

    up->surf_nodes = surf_nodes_ho;
    up->fpo_dgdv_at_surf_nodes = gkyl_array_new(GKYL_DOUBLE, 1, up->num_surf_basis);
  }

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    // Evaluate conf basis at nodes
    struct gkyl_array *conf_nodes_ho = gkyl_array_new(GKYL_DOUBLE,
      up->cdim, up->surf_basis.num_basis);
    up->conf_basis.node_list((double *)gkyl_array_fetch(conf_nodes_ho, 0));

    struct gkyl_array *conf_basis_at_nodes_ho = gkyl_array_new(GKYL_DOUBLE, 
      up->num_conf_basis, up->num_conf_basis);
    for (int k=0; k<up->num_conf_basis; ++k) {
      up->conf_basis.eval((double *)gkyl_array_fetch(conf_nodes_ho,k),
        (double *)gkyl_array_fetch(conf_basis_at_nodes_ho,k));
    }
    up->conf_basis_at_nodes = gkyl_array_cu_dev_new(GKYL_DOUBLE, 
      up->num_conf_basis, up->num_conf_basis);
    gkyl_array_copy(up->conf_basis_at_nodes, conf_basis_at_nodes_ho);
    gkyl_array_release(conf_basis_at_nodes_ho);
    gkyl_array_release(conf_nodes_ho);

    up->surf_nodes = gkyl_array_cu_dev_new(GKYL_DOUBLE, grid->ndim-1, up->surf_basis.num_basis);
    gkyl_array_copy(up->surf_nodes, surf_nodes_ho);
    gkyl_array_release(surf_nodes_ho);

    // Allocate the memory for computing the specific phase and surface nodal to modal calculations
    struct gkyl_mat_mm_array_mem *phase_quad_nodal_to_modal_mem_ho;
    phase_quad_nodal_to_modal_mem_ho = gkyl_mat_mm_array_mem_new(up->num_phase_basis, up->tot_quad, 1.0, 0.0, 
      GKYL_NO_TRANS, GKYL_NO_TRANS, false);

    struct gkyl_mat_mm_array_mem *surf_quad_nodal_to_modal_mem_ho;
    surf_quad_nodal_to_modal_mem_ho = gkyl_mat_mm_array_mem_new(vdim*up->num_surf_basis, vdim*up->tot_surf_quad, 1.0, 0.0, 
      GKYL_NO_TRANS, GKYL_NO_TRANS, false);

    // Compute the matrix A for the phase and surface nodal to modal memory
    const double *phase_w = (const double*) up->weights->data;
    const double *phaseb_o = (const double*) up->basis_at_ords->data;
    for (int n=0; n<up->tot_quad; ++n){
      for (int k=0; k<up->num_phase_basis; ++k){
        gkyl_mat_set(phase_quad_nodal_to_modal_mem_ho->A, k, n, phase_w[n]*phaseb_o[k+up->num_phase_basis*n]);
      }
    }

    // Block matrix for vdim component quad nodal to modal matrix multiplication
    const double *surf_w = (const double*) up->surf_weights->data;
    const double *surfb_o = (const double*) up->surf_basis_at_ords->data;
    for (int n=0; n<vdim*up->tot_surf_quad; ++n){
      for (int k=0; k<vdim*up->num_surf_basis; ++k){
        bool block = !((n-n%up->tot_surf_quad)/up->tot_surf_quad - (k-k%up->num_surf_basis)/up->num_surf_basis);
        if (block) {
          gkyl_mat_set(surf_quad_nodal_to_modal_mem_ho->A, k, n, 
            surf_w[n%up->tot_surf_quad]*surfb_o[(k%up->num_surf_basis)+up->num_surf_basis*(n%up->tot_surf_quad)]);
        }
      }
    }

    // Copy matrix memory to device
    up->phase_quad_nodal_to_modal_mem = gkyl_mat_mm_array_mem_new(up->num_phase_basis, up->tot_quad, 1.0, 0.0, 
      GKYL_NO_TRANS, GKYL_NO_TRANS, up->use_gpu);
    gkyl_mat_copy(up->phase_quad_nodal_to_modal_mem->A, phase_quad_nodal_to_modal_mem_ho->A);
    gkyl_mat_mm_array_mem_release(phase_quad_nodal_to_modal_mem_ho);

    up->surf_quad_nodal_to_modal_mem = gkyl_mat_mm_array_mem_new(vdim*up->num_surf_basis, vdim*up->tot_surf_quad, 1.0, 0.0, 
      GKYL_NO_TRANS, GKYL_NO_TRANS, up->use_gpu);
    gkyl_mat_copy(up->surf_quad_nodal_to_modal_mem->A, surf_quad_nodal_to_modal_mem_ho->A);
    gkyl_mat_mm_array_mem_release(surf_quad_nodal_to_modal_mem_ho);

    // Initialize quadrature data on device
    // Quadrature points for configuration space expansion
    up->tot_conf_quad = init_quad_values(up->cdim, &up->conf_basis, num_quad,
      &up->conf_ordinates, &up->conf_weights, &up->conf_basis_at_ords, up->use_gpu);

    // Quadrature points for phase space expansion
    up->tot_quad = init_quad_values(up->cdim, phase_basis, num_quad,
      &up->ordinates, &up->weights, &up->basis_at_ords, up->use_gpu);
  
    // Initialize quadrature for surface expansion
    up->tot_surf_quad = init_quad_values(up->cdim, &up->surf_basis, num_quad,
      &up->surf_ordinates, &up->surf_weights, &up->surf_basis_at_ords, up->use_gpu);

    // Allocate arrays to store quantities at quadrature points
    up->prim_moms_conf_quad = gkyl_array_cu_dev_new(GKYL_DOUBLE,
      (vdim+2)*up->tot_conf_quad, conf_range->volume);
    up->pot_phase_quad = gkyl_array_cu_dev_new(GKYL_DOUBLE,
      up->tot_quad, phase_range->volume);

    // Primitive moments evaluated at nodes
    up->prim_moms_conf_nodes = gkyl_array_cu_dev_new(GKYL_DOUBLE,
      (vdim+2)*up->num_conf_basis, conf_range->volume);

    // H and G surface quadrature
    up->pot_surf_quad = gkyl_array_cu_dev_new(GKYL_DOUBLE,
      vdim*up->tot_surf_quad, phase_range->volume);
    // dH/dv and d2G/dv2 surface quadrature 
    up->pot_deriv_surf_quad = gkyl_array_cu_dev_new(GKYL_DOUBLE,
      vdim*up->tot_surf_quad, phase_range->volume); 
    // Nodal evaluation of dG/dv at surfaces in transverse directions
    up->fpo_dgdv_at_surf_nodes = gkyl_array_cu_dev_new(GKYL_DOUBLE,
      2*vdim*up->num_surf_basis, phase_range->volume); 

    // Arrays to store computed potentials for copying
    up->sol_pot_surf_modal = gkyl_array_cu_dev_new(GKYL_DOUBLE,
      vdim*up->num_surf_basis, phase_range->volume); 

    // Mappings between phase and configuration space quadrature points 
    int p2c_qidx_ho[up->phase_qrange.volume];
    int surf2c_qidx_ho[up->surf_qrange.volume];
    up->p2c_qidx = (int *) gkyl_cu_malloc(sizeof(int)*up->phase_qrange.volume);
    up->surf2c_qidx = (int *) gkyl_cu_malloc(sizeof(int)*up->surf_qrange.volume);

    struct gkyl_range conf_qrange = get_qrange(up->cdim, up->cdim, num_quad, num_quad_v, is_vdim_p2);

    int pidx[GKYL_MAX_DIM];
    for (int n=0; n<up->tot_quad; ++n) {
      gkyl_range_inv_idx(&up->phase_qrange, n, pidx);
      int cqidx = gkyl_range_idx(&conf_qrange, pidx);
      p2c_qidx_ho[n] = cqidx;
    }
    gkyl_cu_memcpy(up->p2c_qidx, p2c_qidx_ho, sizeof(int)*up->phase_qrange.volume, GKYL_CU_MEMCPY_H2D);

    for (int n=0; n<up->tot_surf_quad; ++n) {
      gkyl_range_inv_idx(&up->surf_qrange, n, pidx);
      int cqidx = gkyl_range_idx(&conf_qrange, pidx);
      surf2c_qidx_ho[n] = cqidx;
    }
    gkyl_cu_memcpy(up->surf2c_qidx, surf2c_qidx_ho, sizeof(int)*up->surf_qrange.volume, GKYL_CU_MEMCPY_H2D);

    // Mapping between phase and configuration space nodes
    int surf2c_nidx_ho[up->num_conf_basis];
    up->surf2c_nidx = (int *) gkyl_cu_malloc(sizeof(int)*up->num_surf_basis);

    struct gkyl_range surf_nrange = get_nrange(up->pdim-1, up->num_surf_basis);
    struct gkyl_range conf_nrange = get_nrange(up->cdim, up->num_conf_basis);
    for(int k=0; k<up->num_surf_basis; ++k) {
      gkyl_range_inv_idx(&surf_nrange, k, pidx);
      int cnidx = gkyl_range_idx(&conf_nrange, pidx);
      surf2c_nidx_ho[k] = cnidx;
    }
    gkyl_cu_memcpy(up->surf2c_nidx, surf2c_nidx_ho, sizeof(int)*up->num_surf_basis, GKYL_CU_MEMCPY_H2D);

    up->surf_basis_dev = gkyl_cart_modal_serendip_cu_dev_new(up->pdim-1, phase_basis->poly_order);
  }
#endif

  return up;
}

void 
gkyl_proj_maxwellian_pots_on_basis_advance(const gkyl_proj_maxwellian_pots_on_basis *up,
  const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
  const struct gkyl_array* prim_moms,
  struct gkyl_array *fpo_h, struct gkyl_array *fpo_g,
  struct gkyl_array *fpo_h_surf, struct gkyl_array *fpo_g_surf,
  struct gkyl_array *fpo_dhdv_surf, struct gkyl_array *fpo_dgdv_surf,
  struct gkyl_array *fpo_d2gdv2_surf)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    return gkyl_proj_maxwellian_pots_on_basis_advance_cu(up, phase_range, conf_range, prim_moms,
      fpo_h, fpo_g, fpo_h_surf, fpo_g_surf, fpo_dhdv_surf, fpo_dgdv_surf, fpo_d2gdv2_surf);
  }
#endif

  // Calculate Maxwellian potentials using primitive moments
  int cdim = up->cdim, pdim = up->pdim;
  int vdim = pdim-cdim;
  int tot_quad = up->tot_quad;
  int num_phase_basis = up->num_phase_basis;
  int num_conf_basis = up->num_conf_basis;
  int num_surf_basis = up->num_surf_basis;

  struct gkyl_range_iter phase_iter, qiter, surf_qiter; 

  int pidx[GKYL_MAX_DIM], rem_dir[GKYL_MAX_DIM] = { 0 };
  for (int d=0; d<cdim; ++d) rem_dir[d] = 1;

  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM];

  // Loop over configuration space cells
  gkyl_range_iter_init(&phase_iter, phase_range);
  while (gkyl_range_iter_next(&phase_iter)) {
    long linc = gkyl_range_idx(conf_range, phase_iter.idx);
    long linp = gkyl_range_idx(phase_range, phase_iter.idx);
    gkyl_copy_int_arr(pdim, phase_iter.idx, pidx);
    gkyl_rect_grid_cell_center(&up->grid, pidx, xc);

    const double *prim_moms_d = gkyl_array_cfetch(prim_moms, linc);
    const double *m0_d = &prim_moms_d[0];
    const double *u_drift_d = &prim_moms_d[num_conf_basis];
    const double *vtsq_d = &prim_moms_d[num_conf_basis*(vdim+1)];

    // Iterate over quadrature points and compute compute potentials H and G
    gkyl_range_iter_init(&qiter, &up->phase_qrange);
    while (gkyl_range_iter_next(&qiter)) {
      long linq = gkyl_range_idx(&up->phase_qrange, qiter.idx);
      const double* quad_node = gkyl_array_cfetch(up->ordinates, linq);

      comp_to_phys(pdim, quad_node, up->grid.dx, xc, xmu);

      // Evaluate moments at quadrature node
      double den_q = up->conf_basis.eval_expand(quad_node, m0_d);
      double vtsq_q = up->conf_basis.eval_expand(quad_node, vtsq_d);

      double rel_speedsq_q = 0.0;
      for (int d=0; d<vdim; ++d) {
        double udrift_q = up->conf_basis.eval_expand(quad_node, &u_drift_d[d*num_conf_basis]);
        rel_speedsq_q += pow(xmu[cdim+d]-udrift_q,2);
      }
      double rel_speed_q = sqrt(rel_speedsq_q);

      double* fpo_h_q = gkyl_array_fetch(up->fpo_h_at_ords, linq);
      double* fpo_g_q = gkyl_array_fetch(up->fpo_g_at_ords, linq);

      fpo_h_q[0] = eval_fpo_h(den_q, rel_speed_q, vtsq_q);
      fpo_g_q[0] = eval_fpo_g(den_q, rel_speed_q, vtsq_q);
    }

    // Project potentials onto basis
    proj_on_basis(up, up->fpo_h_at_ords, gkyl_array_fetch(fpo_h, linp));
    proj_on_basis(up, up->fpo_g_at_ords, gkyl_array_fetch(fpo_g, linp));

    // Iterate over velocity directions for surface expansions at domain boundaries
    for (int d1=0; d1<vdim; ++d1) {
      int dir1 = d1 + cdim;
      bool is_edge_in_dir1 = 0;
      is_edge_in_dir1 = is_edge_in_dir1 || (pidx[dir1] == phase_range->lower[dir1]);
      is_edge_in_dir1 = is_edge_in_dir1 || (pidx[dir1] == phase_range->upper[dir1]);

      if (is_edge_in_dir1) {
        // Compute H, dH/dv, G, d2G/dv2 at dir1 boundaries. 
        // Velocity value at dir1 boundary
        double vbound = pidx[dir1] == phase_range->lower[dir1] ? 
          up->grid.lower[dir1] : up->grid.upper[dir1];

        // Calculate surface expansions using Gauss-Legendre quadrature,
        // similar to proj_on_basis routine.
        gkyl_range_iter_init(&surf_qiter, &up->surf_qrange);
        while (gkyl_range_iter_next(&surf_qiter)) {
          int surf_qidx = gkyl_range_idx(&up->surf_qrange, surf_qiter.idx);
          const double* surf_ord = gkyl_array_cfetch(up->surf_ordinates, surf_qidx);
          surf_comp_to_phys(dir1, pdim, surf_ord, up->grid.dx, xc, xmu);
          xmu[dir1] = vbound;

          // Have to map pdim-1 surface quadrature index to pdim phase quadrature index
          // to get correct phase space variables.
          int edge_idx = pidx[dir1] == phase_range->lower[dir1] ? 0 : up->num_quad-1;
          int phase_idx[GKYL_MAX_DIM];
          edge_idx_to_phase_idx(pdim, dir1, surf_qiter.idx, edge_idx, phase_idx);
          int phase_qidx = gkyl_range_idx(&up->phase_qrange, phase_idx);
          comp_to_phys(pdim, gkyl_array_cfetch(up->ordinates, phase_qidx), up->grid.dx, xc, xmu);
          xmu[dir1] = vbound;

          double den_q = up->conf_basis.eval_expand(surf_ord, m0_d);
          double vtsq_q = up->conf_basis.eval_expand(surf_ord, vtsq_d);

          double rel_speedsq_q = 0.0;
          double rel_vel_in_dir_q = 0.0;
          for (int d=0; d<vdim; ++d) {
            double udrift_q = up->conf_basis.eval_expand(surf_ord, &u_drift_d[d*num_conf_basis]);
              rel_speedsq_q += pow(xmu[cdim+d]-udrift_q, 2);
              if (d == d1)
                rel_vel_in_dir_q += xmu[dir1]-udrift_q;
          }
          double rel_speed_q = sqrt(rel_speedsq_q);

          double* fpo_h_at_surf_ords_q = gkyl_array_fetch(up->fpo_h_at_surf_ords, surf_qidx);
          double* fpo_g_at_surf_ords_q = gkyl_array_fetch(up->fpo_g_at_surf_ords, surf_qidx);
          double* fpo_dhdv_at_surf_ords_q = gkyl_array_fetch(up->fpo_dhdv_at_surf_ords, surf_qidx);
          double* fpo_d2gdv2_at_surf_ords_q = gkyl_array_fetch(up->fpo_d2gdv2_at_surf_ords, surf_qidx);

          fpo_h_at_surf_ords_q[0] = eval_fpo_h(den_q, rel_speed_q, vtsq_q);
          fpo_g_at_surf_ords_q[0] = eval_fpo_g(den_q, rel_speed_q, vtsq_q);
          fpo_dhdv_at_surf_ords_q[0] = eval_fpo_dhdv(den_q, rel_vel_in_dir_q, vtsq_q, rel_speed_q);
          fpo_d2gdv2_at_surf_ords_q[0] = eval_fpo_d2gdv2(den_q, rel_vel_in_dir_q, vtsq_q, rel_speed_q);
        }

        proj_on_surf_basis(up, d1, up->fpo_h_at_surf_ords, gkyl_array_fetch(fpo_h_surf, linp));
        proj_on_surf_basis(up, d1, up->fpo_g_at_surf_ords, gkyl_array_fetch(fpo_g_surf, linp));
        proj_on_surf_basis(up, d1, up->fpo_dhdv_at_surf_ords, gkyl_array_fetch(fpo_dhdv_surf, linp));
        proj_on_surf_basis(up, d1, up->fpo_d2gdv2_at_surf_ords, gkyl_array_fetch(fpo_d2gdv2_surf, linp));

        // Iterate over transverse directions at dir1 boundary for dG/dv
        for (int d2=0; d2<vdim; ++d2) {
          if (d1 == d2) continue;
          int dir2 = d2+cdim;

          // Calculate dG/dv surface expansions with Gauss-Lobatto quadrature for continuity
          // i.e. similar to eval_on_nodes
          for (int i=0; i<num_surf_basis; ++i) {
            const double* surf_node = gkyl_array_cfetch(up->surf_nodes, i);
            surf_comp_to_phys(dir1, pdim, surf_node, up->grid.dx, xc, xmu);
            xmu[dir1] = vbound;

            double den_n = up->conf_basis.eval_expand(surf_node, m0_d);
            double vtsq_n = up->conf_basis.eval_expand(surf_node, vtsq_d);

            double rel_speedsq_n = 0.0;
            double rel_vel_in_dir2_n = 0.0;
            for (int d=0; d<vdim; ++d) {
              double udrift_n = up->conf_basis.eval_expand(surf_node, &u_drift_d[d*num_conf_basis]);
              rel_speedsq_n += pow(xmu[cdim+d]-udrift_n,2);
              if (d == d2)
                rel_vel_in_dir2_n += xmu[dir2]-udrift_n;
            }
            double rel_speed_n = sqrt(rel_speedsq_n);

            double *fpo_dgdv_surf_n = gkyl_array_fetch(up->fpo_dgdv_at_surf_nodes, i);
            fpo_dgdv_surf_n[0] = eval_fpo_dgdv(den_n, rel_vel_in_dir2_n,
              vtsq_n, rel_speedsq_n);
          } 

          double *fpo_dgdv_surf_d = gkyl_array_fetch(fpo_dgdv_surf, linp);

          int num_ret_vals = 1;
          int dgdv_index = (d1 < d2) ? d2-1 : d2;
          dgdv_index += (vdim-1)*d1;
          nod2mod(num_ret_vals, &up->surf_basis, up->fpo_dgdv_at_surf_nodes, 
            &fpo_dgdv_surf_d[dgdv_index*num_surf_basis]);

          gkyl_array_clear(up->fpo_dgdv_at_surf_nodes, 0.0);
        }
      }
    } 
  }
}

void 
gkyl_proj_maxwellian_pots_deriv_on_basis_advance(const gkyl_proj_maxwellian_pots_on_basis *up,
  const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
  const struct gkyl_array *prim_moms,
  struct gkyl_array *fpo_dhdv, struct gkyl_array *fpo_d2gdv2)
{
  // Calculate Maxwellian potentials using primitive moments
  int cdim = up->cdim, pdim = up->pdim;
  int vdim = pdim-cdim;
  int tot_quad = up->tot_quad;
  int num_phase_basis = up->num_phase_basis;
  int num_conf_basis = up->num_conf_basis;
  int num_surf_basis = up->num_surf_basis;

  struct gkyl_range vel_range;
  struct gkyl_range_iter conf_iter, vel_iter;

  int pidx[GKYL_MAX_DIM], rem_dir[GKYL_MAX_DIM] = { 0 };
  for (int d=0; d<conf_range->ndim; ++d) rem_dir[d] = 1;

  // struct gkyl_array *fpo_dgdv_at_surf_nodes = gkyl_array_new(GKYL_DOUBLE, 1, num_surf_basis);

  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM];

  // Loop over configuration space cells for quad integration of moments
  gkyl_range_iter_init(&conf_iter, conf_range);
  while (gkyl_range_iter_next(&conf_iter)) {
    long midx = gkyl_range_idx(conf_range, conf_iter.idx);

    const double *prim_moms_d = gkyl_array_cfetch(prim_moms, midx);
    const double *m0_d = &prim_moms_d[0];
    const double *u_drift_d = &prim_moms_d[num_conf_basis];
    const double *vtsq_d = &prim_moms_d[num_conf_basis*(vdim+1)];

    // Inner loop over velocity space
    gkyl_range_deflate(&vel_range, phase_range, rem_dir, conf_iter.idx);
    gkyl_range_iter_no_split_init(&vel_iter, &vel_range);
    while (gkyl_range_iter_next(&vel_iter)) {
      copy_idx_arrays(conf_range->ndim, phase_range->ndim, conf_iter.idx, vel_iter.idx, pidx);
      gkyl_rect_grid_cell_center(&up->grid, pidx, xc);
 
      long lidx = gkyl_range_idx(&vel_range, vel_iter.idx);

      // Project potential onto basis
      proj_on_basis(up, up->fpo_dhdv_at_ords, gkyl_array_fetch(fpo_dhdv, lidx));

      // Iterate over directions to calculate derivatives of H and G 
      for (int d1=0; d1<vdim; ++d1) {
        for (int d2=0; d2<vdim; ++d2) {
          struct gkyl_range_iter qiter;
          gkyl_range_iter_init(&qiter, &up->phase_qrange);
    
          // Iterate over quadrature points and compute compute potentials
          while (gkyl_range_iter_next(&qiter)) {
            int pqidx = gkyl_range_idx(&up->phase_qrange, qiter.idx);
            const double* quad_node = gkyl_array_cfetch(up->ordinates, pqidx);
    
            comp_to_phys(pdim, quad_node, up->grid.dx, xc, xmu);
    
            // Evaluate moments at quadrature node
            double den_q = up->conf_basis.eval_expand(quad_node, m0_d);
            double u_drift_q = up->conf_basis.eval_expand(quad_node, u_drift_d);
            double vtsq_q = up->conf_basis.eval_expand(quad_node, vtsq_d);
    
            double rel_speedsq_q = 0.0;
            double rel_vel_in_dir1_q = 0.0;
            double rel_vel_in_dir2_q = 0.0;
            for (int d=0; d<vdim; ++d) {
              double udrift_q = up->conf_basis.eval_expand(quad_node, &u_drift_d[d*num_conf_basis]);
              rel_speedsq_q += pow(xmu[cdim+d]-udrift_q,2);
              if (d == d1)
                rel_vel_in_dir1_q += xmu[d1+cdim]-udrift_q;
              if (d == d2)
                rel_vel_in_dir2_q += xmu[d2+cdim]-udrift_q;
            }
            double rel_speed_q = sqrt(rel_speedsq_q);
    
            double* fpo_dhdv_q = gkyl_array_fetch(up->fpo_dhdv_at_ords, pqidx);
            double* fpo_d2gdv2_q = gkyl_array_fetch(up->fpo_d2gdv2_at_ords, pqidx);

            if (d1 == d2) {
              fpo_dhdv_q[0] = eval_fpo_dhdv(den_q, rel_vel_in_dir1_q, vtsq_q, rel_speed_q);

              fpo_d2gdv2_q[0] = eval_fpo_d2gdv2(den_q,
                rel_vel_in_dir1_q, vtsq_q, rel_speed_q);
            }
            else {
              fpo_d2gdv2_q[0] = eval_fpo_d2gdv2_cross(den_q,
                rel_vel_in_dir1_q, rel_vel_in_dir2_q, rel_speed_q, vtsq_q); 
            }
          }
          double *fpo_d2gdv2_d = gkyl_array_fetch(fpo_d2gdv2, lidx);
          proj_on_basis(up, up->fpo_d2gdv2_at_ords, &fpo_d2gdv2_d[(d1*vdim+d2)*up->num_phase_basis]);

          if (d1 == d2) {
            double *fpo_dhdv_d = gkyl_array_fetch(fpo_dhdv, lidx);
            proj_on_basis(up, up->fpo_dhdv_at_ords, &fpo_dhdv_d[d1*up->num_phase_basis]);
          }
        }
      }
    }
  }
}

void 
gkyl_proj_maxwellian_pots_on_basis_release(gkyl_proj_maxwellian_pots_on_basis *up)
{
  gkyl_array_release(up->ordinates);
  gkyl_array_release(up->weights);
  gkyl_array_release(up->basis_at_ords);

  gkyl_array_release(up->surf_ordinates);
  gkyl_array_release(up->surf_weights);
  gkyl_array_release(up->surf_basis_at_ords);

  gkyl_array_release(up->surf_nodes);

  if (up->use_gpu) {
    gkyl_array_release(up->conf_basis_at_nodes);

    gkyl_array_release(up->prim_moms_conf_quad);
    gkyl_array_release(up->prim_moms_conf_nodes);
  
    gkyl_array_release(up->pot_phase_quad);
    gkyl_array_release(up->pot_surf_quad);
    gkyl_array_release(up->pot_deriv_surf_quad);
    gkyl_array_release(up->sol_pot_surf_modal);
  
    gkyl_mat_mm_array_mem_release(up->phase_quad_nodal_to_modal_mem);
    gkyl_mat_mm_array_mem_release(up->surf_quad_nodal_to_modal_mem);

    gkyl_cart_modal_basis_release_cu(up->surf_basis_dev);

    gkyl_array_release(up->conf_ordinates);
    gkyl_array_release(up->conf_weights);
    gkyl_array_release(up->conf_basis_at_ords);
  }
  else {
    gkyl_array_release(up->fpo_h_at_ords);
    gkyl_array_release(up->fpo_g_at_ords);
    gkyl_array_release(up->fpo_dhdv_at_ords);
    gkyl_array_release(up->fpo_d2gdv2_at_ords);

    gkyl_array_release(up->fpo_h_at_surf_ords);
    gkyl_array_release(up->fpo_g_at_surf_ords);
    gkyl_array_release(up->fpo_dhdv_at_surf_ords);
    gkyl_array_release(up->fpo_d2gdv2_at_surf_ords);
    gkyl_array_release(up->fpo_dgdv_at_surf_nodes);
  }

  gkyl_free(up);
}
