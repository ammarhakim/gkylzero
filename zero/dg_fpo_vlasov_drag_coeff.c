#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_fpo_vlasov_drag_coeff.h>
#include <gkyl_dg_fpo_vlasov_drag_coeff_priv.h>
#include <gkyl_util.h>


void gkyl_calc_fpo_drag_coeff_recovery(bool use_gpu, const struct gkyl_rect_grid *grid, 
  struct gkyl_basis pbasis, const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
  const struct gkyl_array* offsets,
  const struct gkyl_array* gamma, const struct gkyl_array* fpo_h, 
  const struct gkyl_array* fpo_dhdv_surf, struct gkyl_array* fpo_drag_coeff,
  struct gkyl_array* fpo_drag_coeff_surf)
{
#ifdef GKYL_HAVE_CUDA
  if (use_gpu)
    return gkyl_calc_fpo_drag_coeff_recovery_cu(grid, pbasis, phase_range, conf_range, offsets, gamma, fpo_h, fpo_dhdv_surf, fpo_drag_coeff, fpo_drag_coeff_surf);
#endif

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
 
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, phase_range);

  while (gkyl_range_iter_next(&iter)) {
    gkyl_copy_int_arr(pdim, iter.idx, idxc);
    gkyl_copy_int_arr(cdim, iter.idx, conf_idxc);

    long linc = gkyl_range_idx(phase_range, idxc);
    long conf_linc = gkyl_range_idx(conf_range, conf_idxc);
 
    const double *fpo_dhdv_surf_c = gkyl_array_cfetch(fpo_dhdv_surf, linc);
    double *fpo_drag_coeff_c = gkyl_array_fetch(fpo_drag_coeff, linc);
    double *fpo_drag_coeff_surf_c = gkyl_array_fetch(fpo_drag_coeff_surf, linc);

    const double *gamma_c = gkyl_array_cfetch(gamma, conf_linc);

    // Iterate through velocity space directions
    for (int d=0; d<vdim; ++d) {
      int dir = d + cdim;

      // Always a 1D, 3-cell stencil.
      long offsets[3] = {0};
      int update_dir[] = {dir};

      // Create offsets from center cell to stencil and index into kernel list.
      create_offsets(1, update_dir, phase_range, idxc, offsets);
      int keri = idx_to_inloup_ker(1, idxc, update_dir, phase_range->upper);

      const double* fpo_h_stencil[3];
      int idx[3][GKYL_MAX_DIM];
      int in_grid = 1;
      for (int i=0; i<3; ++i) {
        gkyl_range_inv_idx(phase_range, linc+offsets[i], idx[i]);
        if (!(idx[i][dir] < phase_range->lower[dir] || idx[i][dir] > phase_range->upper[dir])) {
          fpo_h_stencil[i] = gkyl_array_cfetch(fpo_h, linc+offsets[i]);
        }
      }

      drag_coeff_recovery_stencil[d][keri](
        grid->dx, gamma_c, fpo_h_stencil, 
        fpo_dhdv_surf_c, fpo_drag_coeff_c, fpo_drag_coeff_surf_c);
    }
  }
}

GKYL_CU_DH
void gkyl_calc_fpo_sgn_drag_coeff(struct gkyl_basis pbasis, 
  const struct gkyl_range *phase_range, struct gkyl_array* fpo_drag_coeff_surf, 
  struct gkyl_array* sgn_drag_coeff_surf, struct gkyl_array* const_sgn_drag_coeff_surf) 
{
  int pdim = pbasis.ndim;
  int vdim = 3;
  int cdim = pdim - vdim;

  int poly_order = pbasis.poly_order; 

  fpo_sgn_drag_coeff_t sgn_drag_coeff_stencil[3][3];
 
  // Fetch kernels in each direction
  for (int d=0; d<vdim; ++d) {
    for (int idx=0; idx<3; ++idx) {
      sgn_drag_coeff_stencil[d][idx] = 
       choose_ser_fpo_sgn_drag_coeff_recovery_kern(d, cdim, poly_order, idx);
    }
  }
  
  int idxp[GKYL_MAX_DIM];

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, phase_range);

  while (gkyl_range_iter_next(&iter)) {
    gkyl_copy_int_arr(pdim, iter.idx, idxp);

    long linp = gkyl_range_idx(phase_range, idxp);

    const double *fpo_drag_coeff_surf_c = gkyl_array_cfetch(fpo_drag_coeff_surf, linp);
    double *sgn_drag_coeff_surf_c = gkyl_array_fetch(sgn_drag_coeff_surf, linp);
    int *const_sgn_drag_coeff_surf_c = gkyl_array_fetch(const_sgn_drag_coeff_surf, linp);

    // Iterate through velocity space directions
    for (int d=0; d<vdim; ++d) {
      int dir = d + cdim;
      int update_dir[] = {dir};

      int keri = idx_to_inloup_ker(1, idxp, update_dir, phase_range->upper);

      const_sgn_drag_coeff_surf_c[d] = sgn_drag_coeff_stencil[d][keri](
        fpo_drag_coeff_surf_c, sgn_drag_coeff_surf_c);
    }
  }
}
