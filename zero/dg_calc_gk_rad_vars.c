#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_bin_ops_priv.h>
#include <gkyl_dg_calc_gk_rad_vars.h>
#include <gkyl_dg_calc_gk_rad_vars_priv.h>
#include <gkyl_util.h>

gkyl_dg_calc_gk_rad_vars* gkyl_dg_calc_gk_rad_vars_new(const struct gkyl_rect_grid *phase_grid, 
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, 
  double charge, double mass, const struct gk_geometry *gk_geom, 
  double a, double alpha, double beta, double gamma, double v0) 
{
  gkyl_dg_calc_gk_rad_vars *up = gkyl_malloc(sizeof(gkyl_dg_calc_gk_rad_vars));
  up->phase_grid = *phase_grid;
  up->cdim = conf_basis->ndim;
  up->pdim = phase_basis->ndim;

  up->phase_basis = phase_basis;
  up->conf_basis = conf_basis;

  up->num_conf_basis = conf_basis->num_basis;
  up->num_phase_basis = phase_basis->num_basis;

  up->nodes = gkyl_array_new(GKYL_DOUBLE, phase_grid->ndim, phase_basis->num_basis);
  phase_basis->node_list(gkyl_array_fetch(up->nodes, 0));

  up->charge = charge;
  up->mass = mass;
  up->gk_geom = gkyl_gk_geometry_acquire(gk_geom);

  // Fitting parameters for a given collision type
  up->a = a;
  up->alpha = alpha;
  up->beta = beta;
  up->gamma = gamma;
  up->v0 = v0;

  return up;
}

void gkyl_dg_calc_gk_rad_vars_advance(const struct gkyl_dg_calc_gk_rad_vars *up,
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range, 
  struct gkyl_array* vnu, struct gkyl_array* vsqnu)
{
  // Calculate radiation drag coefficients
  int cdim = up->cdim, pdim = up->pdim;
  int vdim = pdim-cdim;
  int num_phase_basis = up->num_phase_basis;
  int num_conf_basis = up->num_conf_basis;

  struct gkyl_array *vnu_at_nodes = gkyl_array_new(GKYL_DOUBLE, 1, num_phase_basis);
  struct gkyl_array *vsqnu_at_nodes = gkyl_array_new(GKYL_DOUBLE, 1, num_phase_basis);

  int idx[GKYL_MAX_DIM];
  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM];
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, phase_range);

  while (gkyl_range_iter_next(&iter)) {
    gkyl_copy_int_arr(pdim, iter.idx, idx);
    long loc_conf = gkyl_range_idx(conf_range, idx);
    long loc_phase = gkyl_range_idx(phase_range, idx);
    gkyl_rect_grid_cell_center(&up->phase_grid, idx, xc);

    const double *bmag_d = gkyl_array_cfetch(up->gk_geom->bmag, loc_conf);

    for (int i=0; i<num_phase_basis; ++i) {
      const double* node = gkyl_array_cfetch(up->nodes,i);
      comp_to_phys(pdim, node, up->phase_grid.dx, xc, xmu);

      double bmag_n = up->conf_basis->eval_expand(node, bmag_d);
      double vpar = xmu[cdim];
      double mu = 0.0;
      if (vdim > 1)
        mu = xmu[cdim+1];
      double* vnu_at_nodes_n = gkyl_array_fetch(vnu_at_nodes, i);
      double* vsqnu_at_nodes_n = gkyl_array_fetch(vsqnu_at_nodes, i);
      double nu;
      vnu_at_nodes_n[0] = eval_vnu(up->charge, up->mass, 
        up->a, up->alpha, up->beta, up->gamma, up->v0,
				   vpar, mu, bmag_n, &nu);
      vsqnu_at_nodes_n[0] = eval_vsqnu(mu, nu);
    }

    nod2mod(1, up->phase_basis, vnu_at_nodes, gkyl_array_fetch(vnu, loc_phase));
    nod2mod(1, up->phase_basis, vsqnu_at_nodes, gkyl_array_fetch(vsqnu, loc_phase));
    gkyl_array_clear(vnu_at_nodes, 0.0);
    gkyl_array_clear(vsqnu_at_nodes, 0.0); 
  }  
}

void gkyl_dg_calc_gk_rad_vars_release(gkyl_dg_calc_gk_rad_vars *up)
{
  gkyl_array_release(up->nodes);
  gkyl_gk_geometry_release(up->gk_geom);
  gkyl_free(up);
}
