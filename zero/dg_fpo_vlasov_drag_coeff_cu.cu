extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_fpo_vlasov_drag_coeff.h>
#include <gkyl_dg_fpo_vlasov_drag_coeff_priv.h>
#include <gkyl_util.h>
}


__global__ void
gkyl_calc_fpo_drag_coeff_recovery_cu_kernel(const struct gkyl_rect_grid grid,
  struct gkyl_basis pbasis, const struct gkyl_range phase_range, 
  const struct gkyl_range conf_range, const struct gkyl_array *offsets,
  const struct gkyl_array *gamma, const struct gkyl_array* fpo_h,
  const struct gkyl_array* fpo_dhdv_surf, struct gkyl_array* fpo_drag_coeff,
  struct gkyl_array* fpo_drag_coeff_surf)
{
  int pdim = phase_range.ndim;
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

  int idxc[GKYL_MAX_DIM];

  for (unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < phase_range.volume;
      tid += gridDim.x*blockDim.x)
  {
    gkyl_sub_range_inv_idx(&phase_range, tid, idxc);  
    long linp = gkyl_range_idx(&phase_range, idxc);

    const double *fpo_dhdv_surf_c = (const double*)gkyl_array_cfetch(fpo_dhdv_surf, linp);
    double *fpo_drag_coeff_c = (double *)gkyl_array_fetch(fpo_drag_coeff, linp);
    double *fpo_drag_coeff_surf_c = (double *)gkyl_array_fetch(fpo_drag_coeff_surf, linp);

    const long *offsets_d = (long *)gkyl_array_cfetch(offsets, linp);

    const double *gamma_c = (const double *)gkyl_array_cfetch(gamma, linp);

    for (int d=0; d<vdim; ++d) {
      int dir = d + cdim;

      // Always a 1D, 3-cell stencil.
      int update_dir[] = {dir};

      // Create offsets from center cell to stencil and index into kernel list.
      int keri = idx_to_inloup_ker(1, idxc, update_dir, phase_range.upper);

      const double* fpo_h_stencil[3];
      int idx[3][GKYL_MAX_DIM];

      for (int i=0; i<3; ++i) {
        gkyl_range_inv_idx(&phase_range, linp+offsets_d[3*d + i], idx[i]);
        if (!(idx[i][dir] < phase_range.lower[dir] || idx[i][dir] > phase_range.upper[dir])) {
          fpo_h_stencil[i] = (const double *)gkyl_array_cfetch(fpo_h, linp+offsets_d[3*d+i]);
        }
      }

      drag_coeff_recovery_stencil[d][keri](
        grid.dx, gamma_c, fpo_h_stencil, 
        fpo_dhdv_surf_c, fpo_drag_coeff_c, fpo_drag_coeff_surf_c);
    }
  }
}

// Host-side wrapper
void gkyl_calc_fpo_drag_coeff_recovery_cu(const struct gkyl_rect_grid *grid, 
  struct gkyl_basis pbasis, const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
  const struct gkyl_array* offsets,
  const struct gkyl_array* gamma, const struct gkyl_array* fpo_h, 
  const struct gkyl_array* fpo_dhdv_surf, struct gkyl_array* fpo_drag_coeff,
  struct gkyl_array* fpo_drag_coeff_surf)
{
  int nblocks = phase_range->nblocks;
  int nthreads = phase_range->nthreads;

  gkyl_calc_fpo_drag_coeff_recovery_cu_kernel<<<nblocks, nthreads>>>(*grid, pbasis, *phase_range,
    *conf_range, offsets->on_dev, gamma->on_dev, fpo_h->on_dev, fpo_dhdv_surf->on_dev, 
    fpo_drag_coeff->on_dev, fpo_drag_coeff_surf->on_dev);
}
