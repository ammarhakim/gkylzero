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
  struct gkyl_range local, struct gkyl_range local_ext, struct gkyl_range global, struct gkyl_range global_ext, struct gkyl_basis basis)
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

    case GKYL_PMAP_UNIFORM_B:
      gpm->map_x = gkyl_position_map_identity;
      gpm->map_x_ctx = NULL;
      if (pmap_info.map_x == 0)
      { gpm->map_x_backup = gkyl_position_map_identity;
        gpm->map_x_ctx_backup = NULL;
      } else
      { gpm->map_x_backup = pmap_info.map_x;
        gpm->map_x_ctx_backup = pmap_info.ctx_x;
      }

      gpm->map_y = gkyl_position_map_identity;
      gpm->map_y_ctx = NULL;
      if (pmap_info.map_y == 0)
      { gpm->map_y_backup = gkyl_position_map_identity;
        gpm->map_y_ctx_backup = NULL;
      } else
      { gpm->map_y_backup = pmap_info.map_y;
        gpm->map_y_ctx_backup = pmap_info.ctx_y;
      }
      gpm->map_z = gkyl_position_map_identity;
      gpm->map_z_ctx = gpm->constB_ctx;
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
  gkyl_free(gpm);
}


////////////////////////// Utility functions for constant B mapping //////////////////////////

void
gkyl_position_map_gather_bmag_global(struct gkyl_position_map* gpm, struct gkyl_array* bmag)
{
  struct gkyl_array *bmag_global = gkyl_array_new(GKYL_DOUBLE, gpm->basis.num_basis, gpm->global_ext.volume);
  if (gpm->use_gpu) { 
    // If on GPU, allgather is only a GPU operation, so we must copy these arrays to GPU, then back to CPU
    struct gkyl_array *bmag_global_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, gpm->basis.num_basis, gpm->global_ext.volume);
    struct gkyl_array *bmag_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, gpm->basis.num_basis, gpm->local_ext.volume);
    gkyl_array_copy(bmag_dev, bmag);
    gkyl_comm_array_allgather(gpm->comm, &gpm->local, &gpm->global, bmag_dev, bmag_global_dev);
    gkyl_array_copy(bmag_global, bmag_global_dev);
    gkyl_array_release(bmag_global_dev);
    gkyl_array_release(bmag_dev);
  }
  else {
    gkyl_comm_array_allgather(gpm->comm, &gpm->local, &gpm->global, bmag, bmag_global);
  }
  gpm->bmag_global = bmag_global;
}

void
gkyl_position_map_optimize(struct gkyl_position_map* gpm)
{
  if (gpm->id == GKYL_PMAP_UNIFORM_B && gpm->bmag_global != 0)
  {
    gpm->map_x = gpm->map_x_backup;
    gpm->map_x_ctx = gpm->map_x_ctx_backup;
    gpm->map_y = gpm->map_y_backup;
    gpm->map_y_ctx = gpm->map_y_ctx_backup;
    gpm->map_z = gkyl_position_map_constB_z;
    gpm->map_z_ctx = gpm->constB_ctx;

    gpm->bmag_ctx->crange_global = &gpm->global;
    gpm->bmag_ctx->cbasis = &gpm->basis;
    gpm->bmag_ctx->cgrid = &gpm->grid;
    gpm->bmag_ctx->bmag = gkyl_array_acquire(gpm->bmag_global);

    gpm->constB_ctx->psi = (gpm->constB_ctx->psi_min + gpm->constB_ctx->psi_max) / 2;
    gpm->constB_ctx->alpha = (gpm->constB_ctx->alpha_min + gpm->constB_ctx->alpha_max) / 2;
    gpm->constB_ctx->map_strength = gpm->map_strength;

    calculate_mirror_throat_location(gpm->constB_ctx, gpm->bmag_ctx);
    calculate_optimal_mapping(gpm->constB_ctx, gpm->bmag_ctx);
  }
}

void
calculate_mirror_throat_location(struct gkyl_position_map_const_B_ctx *constB_ctx, struct gkyl_bmag_ctx *bmag_ctx)
{
  // Assumes symmetry along theta centered at 0 and two local maxima in Bmag that are symmetric

  enum { X_IDX, Y_IDX, Z_IDX }; // arrangement of cartesian coordinates
  double psi = constB_ctx->psi;
  double alpha = constB_ctx->alpha;
  double *xp = malloc(3*sizeof(double));
  xp[X_IDX] = psi;
  xp[Y_IDX] = alpha;
  int itterations = 10;
  double interval_left = 0.0;
  double interval_right = constB_ctx->theta_max;
  int points_per_level = 20;
  double maximum_Bmag = 0.0;
  double maximum_Bmag_location = 0.0;
  double *fout = malloc(3*sizeof(double));
  for (int j = 0; j < itterations; j++)
  {
    double dz = (interval_right - interval_left) / points_per_level;
    maximum_Bmag = 0.0;
    maximum_Bmag_location = 0.0;
    for (int i = 0; i < points_per_level; i++)
    {
      double z = interval_left + i * dz;
      xp[Z_IDX] = z;
      gkyl_calc_bmag_global(0.0, xp, fout, &bmag_ctx);
      double Bmag = fout[0];
      if (Bmag > maximum_Bmag)
      {
        maximum_Bmag = Bmag;
        maximum_Bmag_location = z;
      }
    }
    interval_left = maximum_Bmag_location - dz;
    interval_right = maximum_Bmag_location + dz;
  }
  constB_ctx->theta_throat = maximum_Bmag_location;
  constB_ctx->Bmag_throat = maximum_Bmag;
  free(xp);
  free(fout);
}

void
calculate_optimal_mapping(struct gkyl_position_map_const_B_ctx *constB_ctx, struct gkyl_bmag_ctx *bmag_ctx)
{
  // Could be refined further by doing midpoint root finding, like the mirror throat finding does
  // Expander region
  enum { X_IDX, Y_IDX, Z_IDX }; // arrangement of cartesian coordinates
  double psi = constB_ctx->psi;
  double alpha = constB_ctx->alpha;
  double *xp = malloc(3*sizeof(double));
  double *fout = malloc(3*sizeof(double));
  xp[X_IDX] = psi;
  xp[Y_IDX] = alpha;
  constB_ctx->map_order_center = 1;
  double scan_cells = 50;
  double scan_left = constB_ctx->theta_throat;
  double scan_right = constB_ctx->theta_max;
  double scan_dxi = (scan_right - scan_left) / scan_cells;
  int expander_order = 1;
  double max_dB_dCell_prior = 99999999.99;
  double max_dB_dCell;
  double max_dB_dCell_order1 = 0.0;
  while (1)
  {
    max_dB_dCell = 0.0;
    constB_ctx->map_order_expander = expander_order;
    for (int iz = 0; iz < scan_cells; iz++)
    {
      double left_xi = scan_left + iz * scan_dxi;
      double right_xi = scan_left + (iz + 1) * scan_dxi;
      double psi = constB_ctx->psi;
      double alpha = constB_ctx->alpha;

      double left_theta[1], right_theta[1];
      gkyl_position_map_constB_z(0.0, &left_xi, left_theta, constB_ctx);
      gkyl_position_map_constB_z(0.0, &right_xi, right_theta, constB_ctx);

      xp[Z_IDX] = left_theta[0];
      gkyl_calc_bmag_global(0.0, xp, fout, bmag_ctx);
      double Bmag_left = fout[0];
      xp[Z_IDX] = right_theta[0];
      gkyl_calc_bmag_global(0.0, xp, fout, bmag_ctx);
      double Bmag_right = fout[0];
      double dB_dCell = (Bmag_right - Bmag_left);
      if (fabs(dB_dCell) > max_dB_dCell)
      {
        max_dB_dCell = fabs(dB_dCell);
      }
    }
    double improvement = max_dB_dCell_prior - max_dB_dCell;
    if (improvement > 1e-3)
    {
      expander_order++;
      max_dB_dCell_prior = max_dB_dCell;
    }
    else if (improvement < 0)
    {
      expander_order--;
      constB_ctx->map_order_expander = expander_order;
      break;
    }
    else
    {
      break;
    }
  }
  double max_dB_dCell_expander = max_dB_dCell;
  //Center region
  scan_left = 0.0;
  scan_right = constB_ctx->theta_throat;
  scan_dxi = (scan_right - scan_left) / scan_cells;
  int center_order = 1;
  max_dB_dCell_prior = 99999999.99;
  while (1)
  {
    max_dB_dCell = 0.0;
    constB_ctx->map_order_center = center_order;
    for (int iz = 0; iz < scan_cells; iz++)
    {
      double left_xi = scan_left + iz * scan_dxi;
      double right_xi = scan_left + (iz + 1) * scan_dxi;

      double left_theta[1], right_theta[1];
      gkyl_position_map_constB_z(0.0, &left_xi, left_theta, constB_ctx);
      gkyl_position_map_constB_z(0.0, &right_xi, right_theta, constB_ctx);

      xp[Z_IDX] = left_theta[0];
      gkyl_calc_bmag_global(0.0, xp, fout, bmag_ctx);
      double Bmag_left = fout[0];
      xp[Z_IDX] = right_theta[0];
      gkyl_calc_bmag_global(0.0, xp, fout, bmag_ctx);
      double Bmag_right = fout[0];

      double dB_dCell = (Bmag_right - Bmag_left);
      if (fabs(dB_dCell) > max_dB_dCell)
      {
        max_dB_dCell = fabs(dB_dCell);
      }
    }
    double improvement = max_dB_dCell_prior - max_dB_dCell;
    if (improvement > 1e-3 & max_dB_dCell > max_dB_dCell_expander)
    {
      center_order++;
      max_dB_dCell_prior = max_dB_dCell;
    }
    else if (improvement < 0)
    {
      center_order--;
      constB_ctx->map_order_center = center_order;
      break;
    }
    else
    {
      break;
    }
  }
  free(xp);
  free(fout);
}

void
gkyl_position_map_constB_z(double t, const double *xn, double *fout, void *ctx)
{
  // Converts our uniform coordinate along field line length to a non-uniform coordinate
  // according to the polynomial mapping of arbitrary order.
  // Notation: I switch from theta to z here. Both are computational coordinate.
  struct gkyl_position_map_const_B_ctx *app = ctx;
  int n_ex = app->map_order_expander;//app->mapping_order_expander;
  int n_ct = app->map_order_center;//app->mapping_order_center;
  double z_min = app->theta_min;
  double z_max = app->theta_max;
  double z_m = app->theta_throat;
  double frac = app->map_strength; // 1 is full mapping, 0 is no mapping
  double uniform_coordinate = xn[0];
  double nonuniform_coordinate, left, right;
  int n;
  if (uniform_coordinate >= z_min && uniform_coordinate <= z_max)
  {
    if (uniform_coordinate <= -z_m)
    {
      left = -z_m;
      right = z_min;
      n = n_ex;
    }
    else if (uniform_coordinate <= 0.0)
    {
      left = -z_m;
      right = 0.0;
      n = n_ct;
    }
    else if (uniform_coordinate <= z_m)
    {
      left = z_m;
      right = 0.0;
      n = n_ct;
    }
    else
    {
      left = z_m;
      right = z_max;
      n = n_ex;
    }
    nonuniform_coordinate = (pow(right - left, 1 - n) * pow(uniform_coordinate - left, n) + left) * frac + uniform_coordinate * (1 - frac);
  }
  else
  {
    nonuniform_coordinate = uniform_coordinate;
  }
  fout[0] = nonuniform_coordinate;
}
