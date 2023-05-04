#include <gkyl_range.h>
#include <gkyl_util.h>
#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_fpo_vlasov_diff_coeff.h>
#include <gkyl_dg_fpo_vlasov_diff_coeff_priv.h>

void gkyl_calc_fpo_diff_coeff_recovery(const struct gkyl_rect_grid *grid, 
  struct gkyl_basis pbasis, const struct gkyl_range *range, 
  const struct gkyl_array *fpo_g, struct gkyl_array *fpo_diff_coeff)
{
  int pdim = pbasis.ndim;
  int vdim = 3;
  int cdim = pdim - vdim; 

  int poly_order = pbasis.poly_order;

  fpo_diff_coeff_t diff_coeff_recovery[3];
  fpo_diff_coeff_surf_t diff_coeff_recovery_surf[3];

  // Fetch kernels in each direction
  for (int d=0; d<vdim; ++d) {
    diff_coeff_recovery[d] = choose_ser_fpo_diff_coeff_recovery_kern(d, cdim, poly_order);
    diff_coeff_recovery_surf[d] = choose_ser_fpo_diff_coeff_recovery_surf_kern(d, cdim, poly_order);
  }

  // Indices in each direction
  int idxl[GKYL_MAX_DIM], idxc[GKYL_MAX_DIM], idxr[GKYL_MAX_DIM];
  int idx_edge[GKYL_MAX_DIM], idx_skin[GKYL_MAX_DIM];
  int edge;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);

  while (gkyl_range_iter_next(&iter)) {
    gkyl_copy_int_arr(pdim, iter.idx, idxc);
    long linc = gkyl_range_idx(range, idxc);

    const double *fpo_g_c = gkyl_array_cfetch(fpo_g, linc);
    double *fpo_diff_coeff_d = gkyl_array_fetch(fpo_diff_coeff, linc);

    // Iterate over velocity space directions
    for (int d=0; d<vdim; ++d) {
      int dir = d + cdim;

      if (idxc[dir] == range->lower[dir] || idxc[dir] == range->upper[dir]) {
        gkyl_copy_int_arr(pdim, iter.idx, idx_edge);
        gkyl_copy_int_arr(pdim, iter.idx, idx_skin);
        edge = (idxc[dir] == range->lower[dir]) ? -1 : 1;
        idx_skin[dir] = idx_edge[dir] - edge;

        long lin_skin = gkyl_range_idx(range, idx_skin);
        long lin_edge = gkyl_range_idx(range, idx_edge);

        const double *fpo_g_skin = gkyl_array_cfetch(fpo_g, lin_skin);
        const double *fpo_g_edge = gkyl_array_cfetch(fpo_g, lin_edge);

        // Surface update
        diff_coeff_recovery_surf[d](edge, grid->dx, fpo_g_skin, fpo_g_edge, fpo_diff_coeff_d);
      }
      else {
        // Volume update
        gkyl_copy_int_arr(pdim, iter.idx, idxl);
        gkyl_copy_int_arr(pdim, iter.idx, idxr);

        idxl[dir] = idxl[dir]-1; idxr[dir] = idxr[dir]+1;

        long linl = gkyl_range_idx(range, idxl);
        long linr = gkyl_range_idx(range, idxr);

        const double *fpo_g_l = gkyl_array_cfetch(fpo_g, linl);
        const double *fpo_g_r = gkyl_array_cfetch(fpo_g, linr);

        diff_coeff_recovery[d](grid->dx, fpo_g_l, fpo_g_c, fpo_g_r, fpo_diff_coeff_d);
      }
    }
  }
}
