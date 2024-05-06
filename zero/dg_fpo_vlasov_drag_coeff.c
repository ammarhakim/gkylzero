#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_fpo_vlasov_drag_coeff.h>
#include <gkyl_dg_fpo_vlasov_drag_coeff_priv.h>
#include <gkyl_util.h>


void gkyl_calc_fpo_drag_coeff_recovery(const struct gkyl_rect_grid *grid, 
  struct gkyl_basis pbasis, const struct gkyl_range *range, const struct gkyl_range *conf_range,
  const struct gkyl_array* gamma, const struct gkyl_array* fpo_h, 
  const struct gkyl_array* fpo_dhdv_surf, struct gkyl_array* fpo_drag_coeff,
  struct gkyl_array* fpo_drag_coeff_surf, struct gkyl_array* sgn_drag_coeff_surf,
  struct gkyl_array* const_sgn_drag_coeff_surf)
{
  int pdim = pbasis.ndim;
  int vdim = 3;
  int cdim = pdim - vdim;

  int poly_order = pbasis.poly_order; 
 
  fpo_drag_coeff_t drag_coeff_recovery_stencil[3][3];
 
  // Fetch kernels in each direction
  for (int d=0; d<vdim; ++d) {
    for (int idx=0; idx<3; ++idx) {
      drag_coeff_recovery_stencil[d][idx] = 
       choose_ser_fpo_drag_coeff_recovery_kern(d, cdim, poly_order, idx);
    }
  }
 
  // Indices in each direction
  int idxl[GKYL_MAX_DIM], idxc[GKYL_MAX_DIM], idxr[GKYL_MAX_DIM], conf_idxc[GKYL_MAX_DIM];
  int idx_edge[GKYL_MAX_DIM], idx_skin[GKYL_MAX_DIM];
  int edge;
 
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);

  while (gkyl_range_iter_next(&iter)) {
    gkyl_copy_int_arr(pdim, iter.idx, idxc);
    gkyl_copy_int_arr(cdim, iter.idx, conf_idxc);

    long linc = gkyl_range_idx(range, idxc);
    long conf_linc = gkyl_range_idx(conf_range, conf_idxc);
 
    const double *fpo_dhdv_surf_c = gkyl_array_cfetch(fpo_dhdv_surf, linc);
    double *fpo_drag_coeff_c = gkyl_array_fetch(fpo_drag_coeff, linc);
    double *fpo_drag_coeff_surf_c = gkyl_array_fetch(fpo_drag_coeff_surf, linc);
    double *sgn_drag_coeff_surf_c = gkyl_array_fetch(sgn_drag_coeff_surf, linc);
    int *const_sgn_drag_coeff_surf_c = gkyl_array_fetch(const_sgn_drag_coeff_surf, linc);

    const double *gamma_c = gkyl_array_cfetch(gamma, conf_linc);

    // Iterate through velocity space directions
    for (int d=0; d<vdim; ++d) {
      int dir = d + cdim;

      // Always a 1D, 3-cell stencil.
      long offsets[3] = {0};
      int update_dir[] = {dir};

      bool is_edge_upper[1], is_edge_lower[1];
      is_edge_lower[0] = idxc[dir] == range->lower[dir]; 
      is_edge_upper[0] = idxc[dir] == range->upper[dir];

      // Create offsets from center cell to stencil and index into kernel list.
      create_offsets(1, is_edge_lower, is_edge_upper, update_dir, range, offsets);
      int keri = idx_to_inloup_ker(1, idxc, update_dir, range->upper);

      const double* fpo_h_stencil[3];
      int idx[3][GKYL_MAX_DIM];
      int in_grid = 1;
      for (int i=0; i<3; ++i) {
        gkyl_range_inv_idx(range, linc+offsets[i], idx[i]);
        if (!(idx[i][dir] < range->lower[dir] || idx[i][dir] > range->upper[dir])) {
          fpo_h_stencil[i] = gkyl_array_cfetch(fpo_h, linc+offsets[i]);
        }
      }

      const_sgn_drag_coeff_surf_c[d] = drag_coeff_recovery_stencil[d][keri](
        grid->dx, gamma_c, fpo_h_stencil, 
        fpo_dhdv_surf_c, fpo_drag_coeff_c, fpo_drag_coeff_surf_c,
        sgn_drag_coeff_surf_c);
    }
  }
}
