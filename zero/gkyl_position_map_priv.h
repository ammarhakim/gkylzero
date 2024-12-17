#include <gkyl_array.h>

struct opt_Theta_ctx
{
  struct gkyl_position_map *gpm;
  struct gkyl_bmag_ctx *bmag_ctx;
  double dB_target;
  double dB_global_lower;
  double B_lower_region;
};

// Allocate double array (filled with zeros).
static struct gkyl_array*
mkarr(bool on_gpu, long nc, long size)
{
  struct gkyl_array* a;
  a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

/**
 * Function that actually frees memory associated with this
 * object when the number of references has decreased to zero.
 *
 * @param ref Reference counter for this object.
 */
void gkyl_position_map_free(const struct gkyl_ref_count *ref);

// Utility functions for polynomial based B mapping

void
calculate_mirror_throat_location_polynomial(struct gkyl_position_map_const_B_ctx *constB_ctx, struct gkyl_bmag_ctx *bmag_ctx)
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
      gkyl_calc_bmag_global(0.0, xp, fout, bmag_ctx);
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
calculate_optimal_mapping_polynomial(struct gkyl_position_map_const_B_ctx *constB_ctx, struct gkyl_bmag_ctx *bmag_ctx)
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
      position_map_constB_z_polynomial(0.0, &left_xi, left_theta, constB_ctx);
      position_map_constB_z_polynomial(0.0, &right_xi, right_theta, constB_ctx);

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
      position_map_constB_z_polynomial(0.0, &left_xi, left_theta, constB_ctx);
      position_map_constB_z_polynomial(0.0, &right_xi, right_theta, constB_ctx);

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
position_map_constB_z_polynomial(double t, const double *xn, double *fout, void *ctx)
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


// Utility functions for numeric B mapping


double
calc_bmag_global_derivative(double xn, void *ctx)
{
  struct gkyl_position_map *gpm = ctx;
  struct gkyl_bmag_ctx *bmag_ctx = gpm->bmag_ctx;
  double h = 1e-6;
  double xh[3];
  double fout[3];
  xh[0] = gpm->constB_ctx->psi;
  xh[1] = gpm->constB_ctx->alpha;
  xh[2] = xn - h;
  gkyl_calc_bmag_global(0.0, xh, fout, bmag_ctx);
  double Bmag_plus = fout[0];
  xh[2] = xn - 2*h;
  gkyl_calc_bmag_global(0.0, xh, fout, bmag_ctx);
  double Bmag_minus = fout[0];
  return (Bmag_plus - Bmag_minus) / (h);
}

void
find_B_field_extrema(struct gkyl_position_map *gpm)
{
  // Assumes we are P1 in z, which means maxima and minima can only be in the center or edge of cells
  struct gkyl_position_map_const_B_ctx *constB_ctx = gpm->constB_ctx;
  struct gkyl_bmag_ctx *bmag_ctx = gpm->bmag_ctx;
  enum { X_IDX, Y_IDX, Z_IDX }; // arrangement of cartesian coordinates
  double psi = constB_ctx->psi;
  double alpha = constB_ctx->alpha;
  double xp[3];
  xp[X_IDX] = psi;
  xp[Y_IDX] = alpha;
  int npts = 2 * constB_ctx->N_theta_boundaries;
  double theta_lo = constB_ctx->theta_min;
  double theta_hi = constB_ctx->theta_max;
  double theta_dxi = (theta_hi - theta_lo) / npts;
  double bmag_vals[npts];
  double dbmag_vals[npts];

  int extrema = 0;
  double theta_extrema[16];
  double bmag_extrema[16];

  for (int i = 0; i <= npts; i++){
    double theta = theta_lo + i * theta_dxi;
    xp[Z_IDX] = theta;
    gkyl_calc_bmag_global(0.0, xp, &bmag_vals[i], bmag_ctx);
    dbmag_vals[i] = calc_bmag_global_derivative(theta, gpm);
    if (i == 0 || i == npts)
    {
      theta_extrema[extrema] = theta;
      bmag_extrema[extrema] = bmag_vals[i];
      extrema++;
    }

    // Minima
    if (dbmag_vals[i] > 0 && dbmag_vals[i-1] < 0){
      if (bmag_vals[i] < bmag_vals[i-1])
      {
        theta_extrema[extrema] = theta;
        bmag_extrema[extrema] = bmag_vals[i];
        extrema++;
      }
      else
      {
        theta_extrema[extrema] = theta - theta_dxi;
        bmag_extrema[extrema] = bmag_vals[i-1];
        extrema++;
      }
    }

    // Maxima
    if (dbmag_vals[i] < 0 && dbmag_vals[i-1] > 0){
      if (bmag_vals[i] > bmag_vals[i-1])
      {
        theta_extrema[extrema] = theta;
        bmag_extrema[extrema] = bmag_vals[i];
        extrema++;
      }
      else
      {
        theta_extrema[extrema] = theta - theta_dxi;
        bmag_extrema[extrema] = bmag_vals[i-1];
        extrema++;
      }
    }
  }
  gpm->constB_ctx->num_extrema = extrema;
  for (int i = 0; i < extrema; i++)
  {
    gpm->constB_ctx->theta_extrema[i] = theta_extrema[i];
    gpm->constB_ctx->bmag_extrema[i] = bmag_extrema[i];
  }
}

double
position_map_numeric_optimization_function(double theta, void *ctx)
{
  struct opt_Theta_ctx *ridders_ctx = ctx;
  struct gkyl_position_map *gpm = ridders_ctx->gpm;
  struct gkyl_bmag_ctx *bmag_ctx = ridders_ctx->bmag_ctx;

  double psi = gpm->constB_ctx->psi;
  double alpha = gpm->constB_ctx->alpha;
  double xp[3];
  xp[0] = psi;
  xp[1] = alpha;
  xp[2] = theta;
  double Bmag;
  gkyl_calc_bmag_global(0.0, xp, &Bmag, bmag_ctx);

  double dB_target = ridders_ctx->dB_target;
  double dB_global_lower = ridders_ctx->dB_global_lower;
  double B_lower_region = ridders_ctx->B_lower_region;

  double result = fabs(Bmag - B_lower_region) + dB_global_lower - dB_target;

  return result;
}


void
position_map_constB_z_numeric(double t, const double *xn, double *fout, void *ctx)
{

  struct gkyl_position_map *gpm = ctx;
  int num_boundaries = gpm->constB_ctx->N_theta_boundaries;
  double *theta_extrema = gpm->constB_ctx->theta_extrema;
  int num_extrema = gpm->constB_ctx->num_extrema;
  double theta_hi = gpm->constB_ctx->theta_max;
  double theta_lo = gpm->constB_ctx->theta_min;
  double theta_range = theta_hi - theta_lo;
  double theta_dxi = theta_range / num_boundaries;
  double theta = xn[0];
  int it = (theta - theta_lo) / theta_dxi;

  if (it ==0 || it == num_boundaries-1)
  {
    fout[0] = theta;
    return;
  }

  // Determine which region theta is in
  // Regions start at 0 and count up to num_extrema-1
  int region = 0;
  for (int i = 1; i <= num_extrema-2; i++)
  {
    if (theta >= theta_extrema[i])
    {
      region = i;
    }
    else
    {
      break;
    }
  }


  double B_tot; // Total change in magnetic field
  for (int i = 0; i <= num_extrema-1; i++)
  {
    B_tot += fabs(gpm->constB_ctx->bmag_extrema[i+1] - gpm->constB_ctx->bmag_extrema[i]);
  }
  double dB_cell = B_tot / (num_boundaries-1);
  double dB_target = dB_cell * it;

  double dB_global_lower = 0.0;
  for (int i = 0; i < region; i++)
  {
    dB_global_lower += fabs(gpm->constB_ctx->bmag_extrema[i+1] - gpm->constB_ctx->bmag_extrema[i]);
  }
  double B_lower_region = gpm->constB_ctx->bmag_extrema[region];

  struct opt_Theta_ctx ridders_ctx = {
    .gpm = gpm,
    .bmag_ctx = gpm->bmag_ctx,
    .dB_target = dB_target,
    .dB_global_lower = dB_global_lower,
    .B_lower_region = B_lower_region,
  };

  double interval_lower = theta_extrema[region];
  double interval_upper = theta_extrema[region+1];
  double interval_lower_eval = position_map_numeric_optimization_function(interval_lower, &ridders_ctx);
  double interval_upper_eval = position_map_numeric_optimization_function(interval_upper, &ridders_ctx);

  if (interval_lower_eval * interval_upper_eval > 0)
  {
    fout[0] = theta;
    return;
  }

  struct gkyl_qr_res res = gkyl_ridders(position_map_numeric_optimization_function, &ridders_ctx,
    interval_lower, interval_upper, interval_lower_eval, interval_upper_eval, 10, 1e-6);
  double Theta = res.res;
  fout[0] = Theta*gpm->constB_ctx->map_strength + theta*(1-gpm->constB_ctx->map_strength);
}
