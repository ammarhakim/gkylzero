extern "C" {
  #include <gkyl_fpo_proj_maxwellian_pots_on_basis.h>
  #include <gkyl_fpo_proj_maxwellian_pots_on_basis_priv.h>
  #include <gkyl_const.h>
  #include <gkyl_range.h>
}

__global__ static void
gkyl_proj_maxwellian_pots_on_basis_advance_cu_ker(const struct gkyl_rect_grid grid,
  const struct gkyl_range phase_range, const struct gkyl_range conf_range,
  const struct gkyl_basis conf_basis,  const struct gkyl_basis phase_basis,
  const struct gkyl_array* GKYL_RESTRICT phase_basis_at_ords, 
  const struct gkyl_array* GKYL_RESTRICT phase_ordinates, 
  const struct gkyl_array* GKYL_RESTRICT phase_weights,
  const struct gkyl_array* GKYL_RESTRICT surf_basis_at_ords, 
  const struct gkyl_array* GKYL_RESTRICT surf_ordinates, 
  const struct gkyl_array* GKYL_RESTRICT surf_weights,
  const struct gkyl_array* GKYL_RESTRICT surf_nodes,
  const int *p2c_qidx,
  const struct gkyl_array* GKYL_RESTRICT prim_moms,
  struct gkyl_array* GKYL_RESTRICT fpo_h,
  struct gkyl_array* GKYL_RESTRICT fpo_g,
  struct gkyl_array* GKYL_RESTRICT fpo_h_surf,
  struct gkyl_array* GKYL_RESTRICT fpo_g_surf,
  struct gkyl_array* GKYL_RESTRICT fpo_dhdv_surf,
  struct gkyl_array* GKYL_RESTRICT fpo_dgdv_surf,
  struct gkyl_array* GKYL_RESTRICT fpo_d2gdv2_surf)
{
  int cdim = conf_range.ndim;
  int pdim = phase_range.ndim;
  int vdim = pdim - cdim;

  int num_conf_basis = prim_moms->ncomp;
  int num_phase_basis = fpo_h->ncomp;
  int num_surf_basis = (fpo_h_surf->ncomp)/vdim;
  int tot_phase_quad = phase_basis_at_ords->size;
  int tot_surf_quad = surf_basis_at_ords->size;

  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM];

  int idxc[GKYL_MAX_DIM], idxp[GKYL_MAX_DIM];

  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < phase_range.volume; tid += blockDim.x*gridDim.x) {
    gkyl_sub_range_inv_idx(&phase_range, tid, idxp);
    for (int i=0; i<cdim; ++i) idxc[i] = idxp[i];

    // Configuration and phase space linear indices
    long linc = gkyl_range_idx(&conf_range, idxc);
    long linp = gkyl_range_idx(&phase_range, idxp);

    // Index into moment arrays
    const double *prim_moms_d = (const double*)gkyl_array_cfetch(prim_moms, linc);
    const double *m0_d = &prim_moms_d[0];
    const double *u_drift_d = &prim_moms_d[num_conf_basis];
    const double *vtsq_d = &prim_moms_d[num_conf_basis*(vdim+1)];

    // Index into output potential arrays
    double *fpo_h_d = (double*)gkyl_array_fetch(fpo_h, linp);
    double *fpo_g_d = (double*)gkyl_array_fetch(fpo_g, linp);
    double *fpo_h_surf_d = (double*)gkyl_array_fetch(fpo_h_surf, linp);
    double *fpo_g_surf_d = (double*)gkyl_array_fetch(fpo_g_surf, linp);
    double *fpo_dhdv_surf_d = (double*)gkyl_array_fetch(fpo_dhdv_surf, linp);
    double *fpo_dgdv_surf_d = (double*)gkyl_array_fetch(fpo_dgdv_surf, linp);
    double *fpo_d2gdv2_surf_d = (double*)gkyl_array_fetch(fpo_d2gdv2_surf, linp);

    // Clear components of potentials
    for (int k=0; k<num_phase_basis; ++k) {
      fpo_h_d[k] = 0.0;
      fpo_g_d[k] = 0.0;
    }

    // The following is modeled after proj_on_basis in the private header.
    const double *phase_weights_arr = (const double *) phase_weights->data;
    const double *phase_basis_at_ords_arr = (const double *) phase_basis_at_ords->data;
    const double *surf_weights_arr = (const double *) surf_weights->data;
    const double *surf_basis_at_ords_arr = (const double *) surf_basis_at_ords->data;

    // Iterate over phase space quadrature points to compute potentials H, G
    for (int n=0; n<tot_phase_quad; ++n) {
      const double* quad_ord = (const double*)gkyl_array_cfetch(phase_ordinates, n);

      comp_to_phys(pdim, quad_ord, grid.dx, xc, xmu);

      // Evaluate moments at quadrature node
      double den_q = conf_basis.eval_expand(quad_ord, m0_d);
      double vtsq_q = conf_basis.eval_expand(quad_ord, vtsq_d);

      double rel_speedsq_q = 0.0;

      for (int d=0; d<vdim; ++d) {
        double udrift_q = conf_basis.eval_expand(quad_ord, &u_drift_d[d*num_conf_basis]);
        rel_speedsq_q += pow(xmu[cdim+d]-udrift_q,2);
      }
      double rel_speed_q = sqrt(rel_speedsq_q);

      double fpo_h_q = eval_fpo_h(den_q, rel_speed_q, vtsq_q);
      double fpo_g_q = eval_fpo_g(den_q, rel_speed_q, vtsq_q);

      // Accumulate quadrature node contribution to expansion coefficients
      double tmp_h = phase_weights_arr[n]*fpo_h_q;
      double tmp_g = phase_weights_arr[n]*fpo_g_q;
      for (int k=0; k<num_phase_basis; ++k) {
        fpo_h_d[k] += tmp_h*phase_basis_at_ords_arr[k+num_phase_basis*n];
        fpo_g_d[k] += tmp_g*phase_basis_at_ords_arr[k+num_phase_basis*n];
      }
    }

  // Iterate over velocity directions for surface expansions at domain boundaries
  for (int d1=0; d1<vdim; ++d1) {
      int dir1 = d1 + cdim;
      bool is_edge_in_dir1 = 0;
      is_edge_in_dir1 = is_edge_in_dir1 || (idxp[dir1] == phase_range.lower[dir1]);
      is_edge_in_dir1 = is_edge_in_dir1 || (idxp[dir1] == phase_range.upper[dir1]);

      // Compute H, dH/dv, G, d2G/dv2 at dir1 boundaries. 
      if (is_edge_in_dir1) {
        // Clear components of potentials
        for (int k=0; k<num_surf_basis; ++k) {
          fpo_h_surf_d[num_surf_basis*d1 + k] = 0.0;
          fpo_g_surf_d[num_surf_basis*d1 + k] = 0.0;
          fpo_dhdv_surf_d[num_surf_basis*d1 + k] = 0.0;
          fpo_d2gdv2_surf_d[num_surf_basis*d1 + k] = 0.0;
        }

        // Velocity value at dir1 boundary
        double vbound = idxp[dir1] == phase_range.lower[dir1] ? 
          grid.lower[dir1] : grid.upper[dir1];

        for (int n=0; n<tot_surf_quad; ++n) {
          const double* surf_ord = (const double*)gkyl_array_cfetch(surf_ordinates, n);

          surf_comp_to_phys(d1+cdim, pdim, surf_ord, grid.dx, xc, xmu);
          xmu[d1] = vbound;

          // Evaluate moments at quadrature node
          double den_q = conf_basis.eval_expand(surf_ord, m0_d);
          double vtsq_q = conf_basis.eval_expand(surf_ord, vtsq_d);

          double rel_speedsq_q = 0.0;
          double rel_vel_in_dir_q = 0.0;
          for (int d=0; d<vdim; ++d) {
            double udrift_q = conf_basis.eval_expand(surf_ord, &u_drift_d[d*num_conf_basis]);
            rel_speedsq_q += pow(xmu[cdim+d]-udrift_q,2);
            if (d == d1)
              rel_vel_in_dir_q += xmu[dir1]-udrift_q;
          }
          double rel_speed_q = sqrt(rel_speedsq_q);

          double fpo_h_q = eval_fpo_h(den_q, rel_speed_q, vtsq_q);
          double fpo_g_q = eval_fpo_g(den_q, rel_speed_q, vtsq_q);
          double fpo_dhdv_q = eval_fpo_dhdv(den_q, rel_vel_in_dir_q, vtsq_q, rel_speed_q);
          double fpo_d2gdv2_q = eval_fpo_d2gdv2(den_q, rel_vel_in_dir_q, vtsq_q, rel_speed_q);

          // Accumulate quadrature node contribution to expansion coefficients
          double tmp_h_surf = surf_weights_arr[n]*fpo_h_q;
          double tmp_g_surf = surf_weights_arr[n]*fpo_g_q;
          double tmp_dhdv_surf = surf_weights_arr[n]*fpo_dhdv_q;
          double tmp_d2gdv2_surf = surf_weights_arr[n]*fpo_d2gdv2_q;
          for (int k=0; k<num_surf_basis; ++k) {
            fpo_h_surf_d[d1*num_surf_basis+k] += tmp_h_surf*surf_basis_at_ords_arr[k+num_surf_basis*n];
            fpo_g_surf_d[d1*num_surf_basis+k] += tmp_g_surf*surf_basis_at_ords_arr[k+num_surf_basis*n];
            fpo_dhdv_surf_d[d1*num_surf_basis+k] += tmp_dhdv_surf*surf_basis_at_ords_arr[k+num_surf_basis*n];
            fpo_d2gdv2_surf_d[d1*num_surf_basis+k] += tmp_d2gdv2_surf*surf_basis_at_ords_arr[k+num_surf_basis*n];
          }

          // Iterate over transverse directions at dir1 boundary for dG/dv
          for (int d2=0; d2<vdim; ++d2) {
            if (d1 == d2) continue;
            int dir2 = d2 + cdim;

            // Calculate dG/dv surface expansions with Gauss-Lobatto quadrature for continuity
            // i.e. similar to eval_on_nodes
            double tmp_dgdv_surf_nodal[256];
            for (int i=0; i<num_surf_basis; ++i) {
              const double *surf_node = (const double*)gkyl_array_cfetch(surf_nodes, i); 

              surf_comp_to_phys(dir1, pdim, surf_node, grid.dx, xc, xmu);
              xmu[dir1] = vbound;

              // Evaluate moments at node
              double den_n = conf_basis.eval_expand(surf_node, m0_d);
              double vtsq_n = conf_basis.eval_expand(surf_node, vtsq_d);

              double rel_speedsq_n = 0.0;
              double rel_vel_in_dir2_n = 0.0;
              for (int d=0; d<vdim; ++d) {
                double udrift_n = conf_basis.eval_expand(surf_node, &u_drift_d[d*num_conf_basis]);
                rel_speedsq_n += pow(xmu[cdim+d]-udrift_n,2);
                if (d == d2)
                  rel_vel_in_dir2_n += xmu[dir2]-udrift_n;
              }
              double rel_speed_n = sqrt(rel_speedsq_n);

              tmp_dgdv_surf_nodal[i] = eval_fpo_dgdv(den_n, rel_vel_in_dir2_n, vtsq_n, rel_speedsq_n);
            }
            
            // Convert nodal evaluation to modal
            int dgdv_index = (d1 < d2) ? d2-1 : d2;
            dgdv_index += (vdim-1)*d1;
            phase_basis.nodal_to_modal(tmp_dgdv_surf_nodal, &fpo_dgdv_surf_d[dgdv_index*num_surf_basis]);
          }
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
  int nblocks = phase_range->nblocks, nthreads = phase_range->nthreads;

  gkyl_proj_maxwellian_pots_on_basis_advance_cu_ker<<<nblocks, nthreads>>>
    (up->grid, *phase_range, *conf_range, *(up->conf_basis), *(up->phase_basis), 
     up->basis_at_ords->on_dev, up->ordinates->on_dev, up->weights->on_dev,
     up->surf_basis_at_ords->on_dev, up->surf_ordinates->on_dev, up->surf_weights->on_dev,
     up->surf_nodes->on_dev, up->p2c_qidx, prim_moms->on_dev, fpo_h->on_dev, fpo_g->on_dev,
     fpo_h_surf->on_dev, fpo_g_surf->on_dev, fpo_dhdv_surf->on_dev, fpo_dhdv_surf->on_dev,
     fpo_d2gdv2_surf->on_dev); 
}
