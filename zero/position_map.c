#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_calc_bmag.h>
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
gkyl_position_map_identity(double t, const double *xn, double *fout, void *ctx)
{
  fout[0] = xn[0];
}

struct gkyl_position_map*
gkyl_position_map_new(struct gkyl_position_map_inp pmap_info, struct gkyl_rect_grid grid,
  struct gkyl_range local, struct gkyl_range local_ext, struct gkyl_range global, struct gkyl_range global_ext,
  struct gkyl_basis basis)
{
  struct gkyl_position_map *gpm = gkyl_malloc(sizeof(*gpm));
  gpm->id = pmap_info.id;

  gpm->bmag_ctx = gkyl_malloc(sizeof(struct gkyl_bmag_ctx));
  gpm->constB_ctx = gkyl_malloc(sizeof(struct gkyl_position_map_const_B_ctx));

  switch (pmap_info.id)
  {
    case GKYL_PMAP_FUNC:
      if (pmap_info.map_x == 0)
      { gpm->map_x = gkyl_position_map_identity;
        gpm->map_x_ctx = NULL;
      } else
      { gpm->map_x = pmap_info.map_x;
        gpm->map_x_ctx = pmap_info.ctx_x;
      }

      if (pmap_info.map_y == 0)
      { gpm->map_y = gkyl_position_map_identity;
        gpm->map_y_ctx = NULL;
      } else
      { gpm->map_y = pmap_info.map_y;
        gpm->map_y_ctx = pmap_info.ctx_y;
      }

      if (pmap_info.map_z == 0)
      { gpm->map_z = gkyl_position_map_identity;
        gpm->map_z_ctx = NULL;
      } else
      { gpm->map_z = pmap_info.map_z;
        gpm->map_z_ctx = pmap_info.ctx_z;
      }

    case GKYL_PMAP_UNIFORM_B_POLYNOMIAL:
      gpm->map_x = gkyl_position_map_identity;
      gpm->map_x_ctx = NULL;
      if (pmap_info.map_x == 0)
      { gpm->constB_ctx->map_x_backup = gkyl_position_map_identity;
        gpm->constB_ctx->map_x_ctx_backup = NULL;
      } else
      { gpm->constB_ctx->map_x_backup = pmap_info.map_x;
        gpm->constB_ctx->map_x_ctx_backup = pmap_info.ctx_x;
      }

      gpm->map_y = gkyl_position_map_identity;
      gpm->map_y_ctx = NULL;
      if (pmap_info.map_y == 0)
      { gpm->constB_ctx->map_y_backup = gkyl_position_map_identity;
        gpm->constB_ctx->map_y_ctx_backup = NULL;
      } else
      { gpm->constB_ctx->map_y_backup = pmap_info.map_y;
        gpm->constB_ctx->map_y_ctx_backup = pmap_info.ctx_y;
      }
      gpm->map_z = gkyl_position_map_identity;
      gpm->map_z_ctx = gpm->constB_ctx;
      gpm->constB_ctx->map_strength = pmap_info.map_strength;

    case GKYL_PMAP_UNIFORM_B_NUMERIC:
      gpm->map_x = gkyl_position_map_identity;
      gpm->map_x_ctx = NULL;
      if (pmap_info.map_x == 0)
      { gpm->constB_ctx->map_x_backup = gkyl_position_map_identity;
        gpm->constB_ctx->map_x_ctx_backup = NULL;
      } else
      { gpm->constB_ctx->map_x_backup = pmap_info.map_x;
        gpm->constB_ctx->map_x_ctx_backup = pmap_info.ctx_x;
      }

      gpm->map_y = gkyl_position_map_identity;
      gpm->map_y_ctx = NULL;
      if (pmap_info.map_y == 0)
      { gpm->constB_ctx->map_y_backup = gkyl_position_map_identity;
        gpm->constB_ctx->map_y_ctx_backup = NULL;
      } else
      { gpm->constB_ctx->map_y_backup = pmap_info.map_y;
        gpm->constB_ctx->map_y_ctx_backup = pmap_info.ctx_y;
      }
      gpm->map_z = gkyl_position_map_identity;
      gpm->map_z_ctx = gpm->constB_ctx;
      gpm->constB_ctx->map_strength = pmap_info.map_strength;
  }

  gpm->grid = grid;
  gpm->local = local;
  gpm->local_ext = local_ext;
  gpm->global = global;
  gpm->global_ext = global_ext;
  gpm->basis = basis;
  gpm->cdim = grid.ndim; 
  gpm->mc2nu = mkarr(false, gpm->cdim*gpm->basis.num_basis, gpm->local_ext.volume);
  gpm->flags = 0;
  GKYL_CLEAR_CU_ALLOC(gpm->flags);
  gpm->ref_count = gkyl_ref_count_init(gkyl_position_map_free);

  struct gkyl_position_map *gpm_out = gpm;
  return gpm_out;
}

void
gkyl_position_map_set(struct gkyl_position_map* gpm, struct gkyl_array* mc2nu)
{
  // Should be a copy, but there are issues I will look into later
  gkyl_array_release(gpm->mc2nu);
  gpm->mc2nu = gkyl_array_acquire(mc2nu);
}

// How do I unit test this function?
void gkyl_position_map_eval_mc2nu(const struct gkyl_position_map* gpm, const double *x_comp, double *x_fa)
{
  int cidx[GKYL_MAX_CDIM];
  for(int i = 0; i < gpm->grid.ndim; i++){
    int idxtemp = gpm->global.lower[i] + (int) floor((x_comp[i] - (gpm->grid.lower[i]) )/gpm->grid.dx[i]);
    idxtemp = GKYL_MAX2(GKYL_MIN2(idxtemp, gpm->local.upper[i]), gpm->local.lower[i]);
    cidx[i] = idxtemp;
  }
  long lidx = gkyl_range_idx(&gpm->local, cidx);
  const double *pmap_coeffs = gkyl_array_cfetch(gpm->mc2nu, lidx);
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

void
gkyl_position_map_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_position_map *gpm = container_of(ref, struct gkyl_position_map, ref_count);
  gkyl_array_release(gpm->mc2nu);
  gkyl_free(gpm->bmag_ctx);
  gkyl_free(gpm->constB_ctx);
  gkyl_free(gpm);
}

void
gkyl_position_map_optimize(struct gkyl_position_map* gpm)
{
  if (gpm->id == GKYL_PMAP_UNIFORM_B_POLYNOMIAL && gpm->bmag_ctx->bmag != 0)
  {
    gpm->map_x     = gpm->constB_ctx->map_x_backup;
    gpm->map_x_ctx = gpm->constB_ctx->map_x_ctx_backup;
    gpm->map_y     = gpm->constB_ctx->map_y_backup;
    gpm->map_y_ctx = gpm->constB_ctx->map_y_ctx_backup;
    gpm->map_z     = position_map_constB_z_polynomial;
    gpm->map_z_ctx = gpm->constB_ctx;

    gpm->bmag_ctx->crange_global = &gpm->global;
    gpm->bmag_ctx->cbasis = &gpm->basis;
    gpm->bmag_ctx->cgrid = &gpm->grid;

    gpm->constB_ctx->psi = (gpm->constB_ctx->psi_min + gpm->constB_ctx->psi_max) / 2;
    gpm->constB_ctx->alpha = (gpm->constB_ctx->alpha_min + gpm->constB_ctx->alpha_max) / 2;

    calculate_mirror_throat_location_polynomial(gpm->constB_ctx, gpm->bmag_ctx);
    calculate_optimal_mapping_polynomial(gpm->constB_ctx, gpm->bmag_ctx);
  }
  else if (gpm->id == GKYL_PMAP_UNIFORM_B_NUMERIC && gpm->bmag_ctx->bmag != 0)
  {
    gpm->map_x     = gpm->constB_ctx->map_x_backup;
    gpm->map_x_ctx = gpm->constB_ctx->map_x_ctx_backup;
    gpm->map_y     = gpm->constB_ctx->map_y_backup;
    gpm->map_y_ctx = gpm->constB_ctx->map_y_ctx_backup;
    gpm->map_z     = position_map_constB_z_numeric;
    gpm->map_z_ctx = gpm;

    gpm->bmag_ctx->crange_global = &gpm->global;
    gpm->bmag_ctx->cbasis = &gpm->basis;
    gpm->bmag_ctx->cgrid = &gpm->grid;

    gpm->constB_ctx->psi = (gpm->constB_ctx->psi_min + gpm->constB_ctx->psi_max) / 2;
    gpm->constB_ctx->alpha = (gpm->constB_ctx->alpha_min + gpm->constB_ctx->alpha_max) / 2;

    find_B_field_extrema(gpm);
  }
}

