#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_mom_pkpm_surf_calc.h>
#include <gkyl_mom_pkpm_surf_calc_priv.h>
#include <gkyl_mom_vlasov_kernels.h>

#include <assert.h>

struct gkyl_mom_pkpm_surf_calc*
gkyl_mom_pkpm_surf_calc_new(const struct gkyl_rect_grid *vel_grid, double mass)
{
  gkyl_mom_pkpm_surf_calc *up = gkyl_malloc(sizeof(gkyl_mom_pkpm_surf_calc));
  up->vel_grid = *vel_grid;
  up->mass = mass;

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->on_dev = up; // self-reference on host

  return up;
}

// pkpm solver is only 1V in velocity space
static inline void
copy_idx_arrays(int cdim, const int *cidx, const int *vidx, int *out)
{
  for (int i=0; i<cdim; ++i)
    out[i] = cidx[i];

  out[cdim] = vidx[0];
}

void
gkyl_mom_pkpm_surf_calc_advance(const struct gkyl_mom_pkpm_surf_calc* calc,
  const struct gkyl_range *conf_edge_rng, const struct gkyl_range *vel_rng, 
  const struct gkyl_range *conf_cell_rng, const struct gkyl_range *phase_cell_rng,
  const struct gkyl_array *GKYL_RESTRICT uin, const struct gkyl_array *GKYL_RESTRICT bin, 
  const struct gkyl_array *GKYL_RESTRICT fin, struct gkyl_array *GKYL_RESTRICT mout)
{
  double xc[GKYL_MAX_DIM]; 
  struct gkyl_range_iter conf_edge_iter, vel_iter;
  
  int cidxl[GKYL_MAX_DIM], cidx[GKYL_MAX_DIM]; // indices for configuration space inputs to left and right of interface
  int pidxl[GKYL_MAX_DIM], pidx[GKYL_MAX_DIM]; // indices for distribution function to left and right of interface

  gkyl_array_clear_range(mout, 0.0, *conf_edge_rng);

  // the outer loop is over configuration space *edges*.
  // For each config-space edge the inner loop walks over the velocity space
  // computing the contribution to the moment
  gkyl_range_iter_init(&conf_edge_iter, conf_edge_rng);
  gkyl_range_iter_init(&vel_iter, vel_rng);
  while (gkyl_range_iter_next(&conf_edge_iter)) {
    long midx = gkyl_range_idx(conf_edge_rng, conf_edge_iter.idx);

    while (gkyl_range_iter_next(&vel_iter)) {
      // only need the cell center in velocity space
      gkyl_rect_grid_cell_center(&calc->vel_grid, vel_iter.idx, xc);

      for (int d=0; d<conf_edge_rng->ndim; ++d) {
        copy_idx_arrays(conf_edge_rng->ndim, conf_edge_iter.idx, vel_iter.idx, pidxl);
        copy_idx_arrays(conf_edge_rng->ndim, conf_edge_iter.idx, vel_iter.idx, pidx);
        gkyl_copy_int_arr(conf_edge_rng->ndim, conf_edge_iter.idx, cidxl);
        gkyl_copy_int_arr(conf_edge_rng->ndim, conf_edge_iter.idx, cidx);

        // need both sides of the configuration space interface
        pidxl[d] = pidxl[d] - 1;
        cidxl[d] = cidxl[d] - 1;

        // fetch linear index of quantities to the left and right of the configuration space interface
        long ccell_idxl = gkyl_range_idx(conf_cell_rng, cidxl);
        long ccell_idx = gkyl_range_idx(conf_cell_rng, cidx);
        long pcell_idxl = gkyl_range_idx(phase_cell_rng, pidxl);
        long pcell_idx = gkyl_range_idx(phase_cell_rng, pidx);

        mom_vlasov_pkpm_surfx_1x1v_ser_p2(xc, calc->vel_grid.dx, pidx, calc->mass,
          gkyl_array_cfetch(uin, ccell_idxl), gkyl_array_cfetch(uin, ccell_idx), 
          gkyl_array_cfetch(bin, ccell_idxl), gkyl_array_cfetch(bin, ccell_idx), 
          gkyl_array_cfetch(fin, pcell_idxl), gkyl_array_cfetch(fin, pcell_idx), 
          gkyl_array_fetch(mout, midx)
        );
      }
    }
  }
}

void gkyl_mom_pkpm_surf_calc_release(gkyl_mom_pkpm_surf_calc* up)
{
  if (GKYL_IS_CU_ALLOC(up->flags))
    gkyl_cu_free(up->on_dev);
  gkyl_free(up);
}

#ifndef GKYL_HAVE_CUDA

void
gkyl_mom_pkpm_surf_calc_advance_cu(const struct gkyl_mom_pkpm_surf_calc* mcalc,
  const struct gkyl_range *conf_edge_rng, const struct gkyl_range *vel_rng, 
  const struct gkyl_range *conf_cell_rng, const struct gkyl_range *phase_cell_rng,
  const struct gkyl_array *GKYL_RESTRICT uin, const struct gkyl_array *GKYL_RESTRICT bin, 
  const struct gkyl_array *GKYL_RESTRICT fin, struct gkyl_array *GKYL_RESTRICT mout)
{
  assert(false);
}

struct gkyl_mom_pkpm_surf_calc*
gkyl_mom_pkpm_surf_calc_cu_dev_new(const struct gkyl_rect_grid *vel_grid, double mass)
{
  assert(false);
}

#endif
