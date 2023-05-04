#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_fpo_vlasov_drag_coeff.h>
#include <gkyl_dg_fpo_vlasov_drag_coeff_priv.h>
#include <gkyl_util.h>


void gkyl_calc_fpo_drag_coeff_recovery(const struct gkyl_rect_grid *grid, 
  struct gkyl_basis pbasis, const struct gkyl_range *range, 
  const struct gkyl_array* fpo_h, struct gkyl_array* fpo_drag_coeff)
{
  int pdim = pbasis.ndim;
  int vdim = 3;
  int cdim = pdim - vdim;

  int poly_order = pbasis.poly_order; 
 
  fpo_drag_coeff_t drag_coeff_recovery[3];
  fpo_drag_coeff_surf_t drag_coeff_recovery_surf[3];
 
  // Fetch kernels in each direction
  for (int d=0; d<vdim; ++d) {
    drag_coeff_recovery[d] = choose_ser_fpo_drag_coeff_recovery_kern(d, cdim, poly_order);
    drag_coeff_recovery_surf[d] = choose_ser_fpo_drag_coeff_recovery_surf_kern(d, cdim, poly_order);
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
 
    const double *fpo_h_c = gkyl_array_cfetch(fpo_h, linc);
    double *fpo_drag_coeff_d = gkyl_array_fetch(fpo_drag_coeff, linc);

    // Iterate through velocity space directions
    for (int d=0; d<vdim; ++d) {
      int dir = d + cdim;
      // Are we at a domain boundary?
      if (idxc[dir] == range->lower[dir] || idxc[dir] == range->upper[dir]) {
        gkyl_copy_int_arr(pdim, iter.idx, idx_edge);
        gkyl_copy_int_arr(pdim, iter.idx, idx_skin);
        edge = (idxc[dir] == range->lower[dir]) ? -1 : 1;
        idx_skin[dir] = idx_edge[dir] - edge;
        
        long lin_skin = gkyl_range_idx(range, idx_skin);
        long lin_edge = gkyl_range_idx(range, idx_edge);

        const double *fpo_h_skin = gkyl_array_cfetch(fpo_h, lin_skin);
        const double *fpo_h_edge = gkyl_array_cfetch(fpo_h, lin_edge);

        drag_coeff_recovery_surf[d](edge, grid->dx, fpo_h_skin, fpo_h_edge, fpo_drag_coeff_d);
      } 
      else {
        // Volume update
        gkyl_copy_int_arr(pdim, iter.idx, idxl);
        gkyl_copy_int_arr(pdim, iter.idx, idxr);
  
        idxl[dir] = idxl[dir]-1; idxr[dir] = idxr[dir]+1;
  
        long linl = gkyl_range_idx(range, idxl);
        long linr = gkyl_range_idx(range, idxr);
  
        const double *fpo_h_l = gkyl_array_cfetch(fpo_h, linl);
        const double *fpo_h_r = gkyl_array_cfetch(fpo_h, linr);

        drag_coeff_recovery[d](grid->dx, fpo_h_l, fpo_h_c, fpo_h_r, fpo_drag_coeff_d);
      }
    }
   }
}
