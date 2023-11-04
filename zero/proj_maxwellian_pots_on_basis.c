#include "gkyl_eval_on_nodes.h"
#include <math.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_const.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_proj_maxwellian_pots_on_basis.h>
#include <gkyl_proj_maxwellian_pots_on_basis_priv.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

static void
nod2mod(int num_ret_vals, const struct gkyl_basis *basis, const struct gkyl_array *fun_at_nodes, double *f) {
  const double *fao = gkyl_array_cfetch(fun_at_nodes, 0);

  int num_basis = basis->num_basis;
  double fnodal[num_basis];
  for (int i=0; i<num_ret_vals; ++i) {
    for (int k=0; k<num_basis; ++k) {
      fnodal[k] = fao[k];
    }

    basis->nodal_to_modal(fnodal, &f[num_basis*i]);
  }
}

gkyl_proj_maxwellian_pots_on_basis* gkyl_proj_maxwellian_pots_on_basis_new(const struct gkyl_rect_grid *grid, 
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, const struct gkyl_basis *surf_basis) 
{
  gkyl_proj_maxwellian_pots_on_basis *up = gkyl_malloc(sizeof(gkyl_proj_maxwellian_pots_on_basis));
  up->grid = *grid;
  up->cdim = conf_basis->ndim;
  up->pdim = phase_basis->ndim;

  up->phase_basis = phase_basis;
  up->conf_basis = conf_basis;
  up->surf_basis = surf_basis;

  up->num_conf_basis = conf_basis->num_basis;
  up->num_phase_basis = phase_basis->num_basis;
  up->num_surf_basis = surf_basis->num_basis;

  up->surf_nodes = gkyl_array_new(GKYL_DOUBLE, grid->ndim-1, surf_basis->num_basis);
  surf_basis->node_list(gkyl_array_fetch(up->surf_nodes, 0));

  up->nodes = gkyl_array_new(GKYL_DOUBLE, grid->ndim, phase_basis->num_basis);
  phase_basis->node_list(gkyl_array_fetch(up->nodes, 0));
  return up;
}

void gkyl_proj_maxwellian_pots_on_basis_lab_mom(const gkyl_proj_maxwellian_pots_on_basis *up,
    const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
    const struct gkyl_array* m0, const struct gkyl_array* u_drift, const struct gkyl_array* vtsq,
    const struct gkyl_array *gamma, const double mass, 
    struct gkyl_array *fpo_h, struct gkyl_array *fpo_g,
    struct gkyl_array *fpo_h_surf, struct gkyl_array *fpo_g_surf,
    struct gkyl_array *fpo_dhdv_surf, struct gkyl_array *fpo_d2gdv2_surf)
{
  // Calculate Maxwellian potentials using primitive moments
  int cdim = up->cdim, pdim = up->pdim;
  int vdim = pdim-cdim;
  int num_phase_basis = up->num_phase_basis;
  int num_conf_basis = up->num_conf_basis;
  int num_surf_basis = up->num_surf_basis;

  struct gkyl_range vel_range;
  struct gkyl_range_iter conf_iter, vel_iter;

  int pidx[GKYL_MAX_DIM], rem_dir[GKYL_MAX_DIM] = { 0 };
  for (int d=0; d<conf_range->ndim; ++d) rem_dir[d] = 1;

  struct gkyl_array *fpo_g_at_nodes = gkyl_array_new(GKYL_DOUBLE, 1, num_phase_basis);
  struct gkyl_array *fpo_h_at_nodes = gkyl_array_new(GKYL_DOUBLE, 1, num_phase_basis);
  struct gkyl_array *fpo_g_at_surf_nodes = gkyl_array_new(GKYL_DOUBLE, 1, num_surf_basis);
  struct gkyl_array *fpo_h_at_surf_nodes = gkyl_array_new(GKYL_DOUBLE, 1, num_surf_basis);
  struct gkyl_array *fpo_dhdv_at_surf_nodes = gkyl_array_new(GKYL_DOUBLE, 1, num_surf_basis);
  struct gkyl_array *fpo_d2gdv2_at_surf_nodes = gkyl_array_new(GKYL_DOUBLE, 1, num_surf_basis);

  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM];

  // Loop over configuration space cells for quad integration of moments
  gkyl_range_iter_init(&conf_iter, conf_range);
  while (gkyl_range_iter_next(&conf_iter)) {
    long midx = gkyl_range_idx(conf_range, conf_iter.idx);

    const double *m0_d = gkyl_array_cfetch(m0, midx);
    const double *u_drift_d = gkyl_array_cfetch(u_drift, midx);
    const double *vtsq_d = gkyl_array_cfetch(vtsq, midx);
    const double *gamma_d = gkyl_array_cfetch(gamma, midx);

    // Inner loop over velocity space
    // Should be able to replace this with a phase space loop
    // Since we don't have the quadrature anymore
    gkyl_range_deflate(&vel_range, phase_range, rem_dir, conf_iter.idx);
    gkyl_range_iter_no_split_init(&vel_iter, &vel_range);
    while (gkyl_range_iter_next(&vel_iter)) {
      copy_idx_arrays(conf_range->ndim, phase_range->ndim, conf_iter.idx, vel_iter.idx, pidx);
      gkyl_rect_grid_cell_center(&up->grid, pidx, xc);
 
      long lidx = gkyl_range_idx(&vel_range, vel_iter.idx);

      for (int i=0; i<num_phase_basis; ++i) {
        const double* node = gkyl_array_cfetch(up->nodes,i);
        comp_to_phys(pdim, node, up->grid.dx, xc, xmu);

        double den_n = up->conf_basis->eval_expand(node, m0_d);
        double u_drift_n = up->conf_basis->eval_expand(node, u_drift_d);
        double vtsq_n = up->conf_basis->eval_expand(node, vtsq_d);
        
        double rel_speedsq_n = 0.0;
        for (int d=0; d<vdim; ++d) {
          double udrift_n = up->conf_basis->eval_expand(node, &u_drift_d[d*num_conf_basis]);
          rel_speedsq_n += pow(xmu[cdim+d]-udrift_n,2);
        }
        double rel_speed_n = sqrt(rel_speedsq_n);

        double* fpo_h_at_nodes_n = gkyl_array_fetch(fpo_h_at_nodes, i);
        double* fpo_g_at_nodes_n = gkyl_array_fetch(fpo_g_at_nodes, i);
        fpo_h_at_nodes_n[0] = eval_fpo_h(gamma_d[0], den_n, rel_speed_n, vtsq_n);
        fpo_g_at_nodes_n[0] = eval_fpo_g(gamma_d[0], den_n, rel_speed_n, vtsq_n);
      }

      nod2mod(1, up->phase_basis, fpo_h_at_nodes, gkyl_array_fetch(fpo_h, lidx));
      nod2mod(1, up->phase_basis, fpo_g_at_nodes, gkyl_array_fetch(fpo_g, lidx));
      gkyl_array_clear(fpo_h_at_nodes, 0.0);
      gkyl_array_clear(fpo_g_at_nodes, 0.0); 
   
      // Check if we're at a velocity space boundary in each velocity direction
      for (int d1=0; d1<vdim; ++d1) {
        int dir1 = d1 + cdim;
  
        if (pidx[dir1] == vel_range.lower[d1] || pidx[dir1] == vel_range.upper[d1]) { 
          // Velocity value at boundary for surface expansion
          int vmax = pidx[dir1] == vel_range.lower[d1] ? up->grid.lower[dir1] : up->grid.upper[dir1]; 
  
          for (int i=0; i<num_surf_basis; ++i) {
            const double* surf_node = gkyl_array_cfetch(up->surf_nodes, i);
            surf_comp_to_phys(dir1, pdim, surf_node, up->grid.dx, xc, xmu);
            xmu[dir1] = vmax;

            double node[pdim];
            for (int i=0; i<pdim; ++i) {
              if (i < dir1)
                node[i] = surf_node[i];
              else if (i == dir1)
                node[i] = vmax > 0 ? 1 : -1;
              else if (i > dir1)
                node[i] = surf_node[i-1];
            }

            double den_n = up->conf_basis->eval_expand(node, m0_d);
            double u_drift_n = up->conf_basis->eval_expand(node, u_drift_d);
            double vtsq_n = up->conf_basis->eval_expand(node, vtsq_d);

            double rel_speedsq_n = 0.0;
            double rel_vel_in_dir_n = 0.0;
            for (int d=0; d<vdim; ++d) {
              double udrift_n = up->conf_basis->eval_expand(node, &u_drift_d[d*num_conf_basis]);
              rel_speedsq_n += pow(xmu[cdim+d]-udrift_n,2);
              if (d == d1)
                  rel_vel_in_dir_n += xmu[dir1]-udrift_n;
            }
            double rel_speed_n = sqrt(rel_speedsq_n);
  
            double* fpo_h_at_surf_nodes_n = gkyl_array_fetch(fpo_h_at_surf_nodes, i);
            double* fpo_g_at_surf_nodes_n = gkyl_array_fetch(fpo_g_at_surf_nodes, i);
            double* fpo_dhdv_at_surf_nodes_n = gkyl_array_fetch(fpo_dhdv_at_surf_nodes, i);
            fpo_h_at_surf_nodes_n[0] = eval_fpo_h(gamma_d[0], den_n, rel_speed_n, vtsq_n);
            fpo_g_at_surf_nodes_n[0] = eval_fpo_g(gamma_d[0], den_n, rel_speed_n, vtsq_n);
            fpo_dhdv_at_surf_nodes_n[0] = eval_fpo_dhdv(gamma_d[0], den_n, rel_vel_in_dir_n, vtsq_n, rel_speed_n);
          }
         
          double *fpo_h_surf_n = gkyl_array_fetch(fpo_h_surf, lidx);
          double *fpo_g_surf_n = gkyl_array_fetch(fpo_g_surf, lidx);
          double *fpo_dhdv_surf_n = gkyl_array_fetch(fpo_dhdv_surf, lidx);
          nod2mod(1, up->surf_basis, fpo_h_at_surf_nodes, &fpo_h_surf_n[d1*num_surf_basis]);
          nod2mod(1, up->surf_basis, fpo_g_at_surf_nodes, &fpo_g_surf_n[d1*num_surf_basis]);
          nod2mod(1, up->surf_basis, fpo_dhdv_at_surf_nodes, &fpo_dhdv_surf_n[d1*num_surf_basis]);

          gkyl_array_clear(fpo_h_at_surf_nodes, 0.0);
          gkyl_array_clear(fpo_g_at_surf_nodes, 0.0); 
          gkyl_array_clear(fpo_dhdv_at_surf_nodes, 0.0); 
  
          // Iterate over second velocity space direction to calculate d2G/dv2
          for (int d2=0; d2<vdim; ++d2) {
            int dir2 = d2 + cdim; 
            for (int i=0; i<num_surf_basis; ++i) {
              const double* surf_node = gkyl_array_cfetch(up->surf_nodes, i);
              surf_comp_to_phys(dir1, pdim, surf_node, up->grid.dx, xc, xmu);
              xmu[dir1] = vmax;

              double node[pdim];
              for (int i=0; i<pdim; ++i) {
                if (i < dir1)
                  node[i] = surf_node[i];
                else if (i == dir1)
                  node[i] = vmax > 0 ? 1 : -1;
                else if (i > dir1)
                  node[i] = surf_node[i-1];
              }

              double den_n = up->conf_basis->eval_expand(node, m0_d);
              double u_drift_n = up->conf_basis->eval_expand(node, u_drift_d);
              double vtsq_n = up->conf_basis->eval_expand(node, vtsq_d);
  
              double rel_speedsq_n = 0.0;
              double rel_vel_in_dir1_n = 0.0;
              double rel_vel_in_dir2_n = 0.0;
              for (int d=0; d<vdim; ++d) {
                double udrift_n = up->conf_basis->eval_expand(node, &u_drift_d[d*num_conf_basis]);
                rel_speedsq_n += pow(xmu[cdim+d]-udrift_n,2);
                if (d == d1)
                  rel_vel_in_dir1_n += xmu[cdim+d]-udrift_n;
                if (d == d2)
                  rel_vel_in_dir2_n += xmu[cdim+d]-udrift_n;
              }
              double rel_speed_n = sqrt(rel_speedsq_n);

              double *fpo_d2gdv2_surf_n = gkyl_array_fetch(fpo_d2gdv2_at_surf_nodes, i);

              if (d1 == d2) {
                fpo_d2gdv2_surf_n[0] = eval_fpo_d2gdv2(gamma_d[0], den_n, 
                  rel_vel_in_dir1_n, 
                  vtsq_n, rel_speed_n);
              }
              else {
                fpo_d2gdv2_surf_n[0] = eval_fpo_d2gdv2_cross(gamma_d[0], den_n,
                  rel_vel_in_dir1_n, rel_vel_in_dir2_n, rel_speed_n, vtsq_n); 
              } 
            } 

            double *fpo_d2gdv2_surf_n = gkyl_array_fetch(fpo_d2gdv2_surf, lidx);
            nod2mod(1, up->surf_basis, fpo_d2gdv2_at_surf_nodes, &fpo_d2gdv2_surf_n[(d1*vdim+d2)*num_surf_basis]);
            gkyl_array_clear(fpo_d2gdv2_at_surf_nodes, 0.0);
          } 
        }
      } 
    }
  }
  gkyl_array_release(fpo_g_at_nodes);
  gkyl_array_release(fpo_h_at_nodes);
  gkyl_array_release(fpo_g_at_surf_nodes);
  gkyl_array_release(fpo_h_at_surf_nodes);
  gkyl_array_release(fpo_dhdv_at_surf_nodes);
  gkyl_array_release(fpo_d2gdv2_at_surf_nodes);
}

void gkyl_proj_maxwellian_pots_on_basis_release(gkyl_proj_maxwellian_pots_on_basis *up)
{
  gkyl_array_release(up->nodes);
  gkyl_array_release(up->surf_nodes);
  gkyl_free(up);
}
