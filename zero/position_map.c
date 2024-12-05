#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_comm_io.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_position_map.h>
#include <gkyl_position_map_priv.h>

#include <float.h>
#include <assert.h>

// Remove with the print statements at the bottom
#include <gkyl_util.h>


void
gkyl_position_map_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_position_map *gpm = container_of(ref, struct gkyl_position_map, ref_count);
  gkyl_array_release(gpm->pmap);
  gkyl_free(gpm);
}

struct gkyl_position_map*
gkyl_position_map_new(struct gkyl_nonuniform_position_map_info mapc2p_in, struct gkyl_rect_grid grid,
  struct gkyl_range local, struct gkyl_range local_ext, struct gkyl_range global, struct gkyl_basis basis)
{
  struct gkyl_position_map *gpm = gkyl_malloc(sizeof(*gpm));
  gpm->is_identity = (mapc2p_in.mapping == 0 && (fabs(mapc2p_in.numerical_mapping_fraction) < 1e-14));
  gpm->grid = grid;
  gpm->local = local;
  gpm->local_ext = local_ext;
  gpm->global = global;
  gpm->basis = basis;
  int cdim = grid.ndim; 
  gpm->pmap = mkarr(false, cdim*gpm->basis.num_basis, gpm->local_ext.volume);
  gpm->flags = 0;
  GKYL_CLEAR_CU_ALLOC(gpm->flags);
  gpm->ref_count = gkyl_ref_count_init(gkyl_position_map_free);

  struct gkyl_position_map *gpm_out = gpm;
  return gpm_out;
}

void
gkyl_position_map_set(struct gkyl_position_map* gpm, struct gkyl_array* pmap)
{
  // Should be a copy, but there are issues I will look into later
  gkyl_array_release(gpm->pmap);
  gpm->pmap = gkyl_array_acquire(pmap);
}

// How do I unit test this function?
void gkyl_position_map_eval_c2p(const struct gkyl_position_map* gpm, const double *x_comp, double *x_fa)
{
  if (gpm->is_identity) {
    for (int i=0; i<gpm->grid.ndim; i++) {
      x_fa[i] = x_comp[i];
    };
  } else
  {
    int cidx[GKYL_MAX_CDIM];
    for(int i = 0; i < gpm->grid.ndim; i++){
      int idxtemp = gpm->global.lower[i] + (int) floor((x_comp[i] - (gpm->grid.lower[i]) )/gpm->grid.dx[i]);
      idxtemp = GKYL_MAX2(GKYL_MIN2(idxtemp, gpm->local.upper[i]), gpm->local.lower[i]);
      cidx[i] = idxtemp;
    }
    long lidx = gkyl_range_idx(&gpm->local, cidx);
    const double *pmap_coeffs = gkyl_array_cfetch(gpm->pmap, lidx);
    double cxc[gpm->grid.ndim];
    double x_log[gpm->grid.ndim];
    gkyl_rect_grid_cell_center(&gpm->grid, cidx, cxc);
    for(int i = 0; i < gpm->grid.ndim; i++){
      x_log[i] = (x_comp[i]-cxc[i])/(gpm->grid.dx[i]*0.5);
    }
    double xyz_fa[3];
    for(int i = 0; i < 3; i++){
      xyz_fa[i] = gpm->basis.eval_expand(x_log, &pmap_coeffs[i*gpm->basis.num_basis]);
    }
    for (int i=0; i<gpm->grid.ndim; i++) {
      x_fa[i] = xyz_fa[i];
    }
    x_fa[gpm->grid.ndim-1] = xyz_fa[2];
  }
}


struct gkyl_position_map*
gkyl_position_map_acquire(const struct gkyl_position_map* gpm)
{
  gkyl_ref_count_inc(&gpm->ref_count);
  return (struct gkyl_position_map*) gpm;
}

void
gkyl_position_map_release(const struct gkyl_position_map *gpm)
{
  gkyl_ref_count_dec(&gpm->ref_count);
}
