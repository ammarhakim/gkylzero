#include <gkyl_array.h>
#include <gkyl_calc_bmag.h>

// Context for numeric root finding B mapping
struct opt_Theta_ctx
{
  struct gkyl_position_map *gpm;
  struct gkyl_bmag_ctx *bmag_ctx;
  double dB_target; // How much B should change in 1 cell
  double dB_global_lower; // How much dB came before this region
  double B_lower_region; // The B value at the lower end of the region
};

/**
 * Function that actually frees memory associated with this
 * object when the number of references has decreased to zero.
 *
 * @param ref Reference counter for this object.
 */
static void gkyl_position_map_free(const struct gkyl_ref_count *ref);

// Utility functions for polynomial based B mapping

/**
 * Calculates the location of the throat of the mirror in the theta direction
 * and the value of Bmag at that location.
 * 
 * @param constB_ctx Context for the constant B mapping
 * @param bmag_ctx Context for the magnetic field calculation
 */
static void
calculate_mirror_throat_location_polynomial(struct gkyl_position_map_const_B_ctx *constB_ctx, struct gkyl_bmag_ctx *bmag_ctx)
{
  // Parameters to use for the midpoint rule root finding algorithm to find the throat of the mirror
  int itterations = 10;
  int points_per_level = 20;

  // Assumes symmetry along theta, centered at 0, and two local maxima in Bmag that are symmetric
  enum { X_IDX, Y_IDX, Z_IDX }; // arrangement of cartesian coordinates
  double psi = constB_ctx->psi;
  double alpha = constB_ctx->alpha;
  double xp[3];
  xp[X_IDX] = psi;
  xp[Y_IDX] = alpha;

  double interval_left = 0.0;
  double interval_right = constB_ctx->theta_max;
  double maximum_Bmag = 0.0;
  double maximum_Bmag_location = 0.0;
  double fout[3];
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
}

/**
 * Converts our uniform coordinate along field line length to a non-uniform coordinate
 * according to the polynomial mapping of arbitrary order.
 * Notation: We switch from theta to z here. Both are the third computational coordinate.
 * 
 * @param t Time
 * @param xn Uniform coordinate
 * @param fout Non-uniform coordinate
 * @param ctx position_map_constB_ctx context for the constant B mapping
 */
static void
position_map_constB_z_polynomial(double t, const double *xn, double *fout, void *ctx)
{
  struct gkyl_position_map_const_B_ctx *app = ctx;
  int n_ex = app->map_order_expander;
  int n_ct = app->map_order_center;
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

/**
 * Calculates the optimal orders for the polynomail mapping of the B field.
 * The optimal orders are the orders that minimize the maximum dB/dTheta in that region.
 * 
 * @param constB_ctx Context for the constant B mapping
 * @param bmag_ctx Context for the magnetic field calculation
 */
static void
calculate_optimal_mapping_polynomial(struct gkyl_position_map_const_B_ctx *constB_ctx, struct gkyl_bmag_ctx *bmag_ctx)
{
  // Could be refined further by doing midpoint root finding for maximum dB/dz
  // Expander region
  enum { X_IDX, Y_IDX, Z_IDX }; // arrangement of cartesian coordinates
  double psi = constB_ctx->psi;
  double alpha = constB_ctx->alpha;
  double xp[3];
  double fout[3];
  xp[X_IDX] = psi;
  xp[Y_IDX] = alpha;
  constB_ctx->map_order_center = 1;
  double scan_cells = 10;
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
    if (improvement > 1e-3)
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
}


// Utility functions for numeric root finding B mapping

/**
 * Calculates dB/dTheta numerically at a given xn value. Calculates
 * the derivative to the left of the point.
 * 
 * @param theta The theta value to calculate the derivative at
 * @param ctx The context for the position map
 */
static double
calc_bmag_global_derivative(double theta, void *ctx)
{
  struct gkyl_position_map *gpm = ctx;
  struct gkyl_bmag_ctx *bmag_ctx = gpm->bmag_ctx;
  double dtheta_cell = (gpm->constB_ctx->theta_max - gpm->constB_ctx->theta_min)/gpm->constB_ctx->N_theta_boundaries;
  double h = 1e-2 * dtheta_cell;
  double xh[3];
  double fout[3];
  xh[0] = gpm->constB_ctx->psi;
  xh[1] = gpm->constB_ctx->alpha;
  xh[2] = theta - h;
  gkyl_calc_bmag_global(0.0, xh, fout, bmag_ctx);
  double Bmag_plus = fout[0];
  xh[2] = theta - 2*h;
  gkyl_calc_bmag_global(0.0, xh, fout, bmag_ctx);
  double Bmag_minus = fout[0];
  return (Bmag_plus - Bmag_minus) / (h);
}

/**
 * Finds the local min and max of the B field along the field line
 * specified by the input psi and alpha values.
 * 
 * @param gpm The position map object
 */
static void
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
  double *bmag_vals = gkyl_malloc(sizeof(double) * (npts + 1));
  double *dbmag_vals = gkyl_malloc(sizeof(double) * (npts + 1));

  int extrema = 1; // Offset by 1 for the first point
  double *theta_extrema = gkyl_malloc(sizeof(double) * (npts + 1));
  double *bmag_extrema = gkyl_malloc(sizeof(double) * (npts + 1));

  for (int i = 0; i <= npts; i++){
    double theta = theta_lo + i * theta_dxi;
    xp[Z_IDX] = theta;
    gkyl_calc_bmag_global(0.0, xp, &bmag_vals[i], bmag_ctx);
    dbmag_vals[i] = calc_bmag_global_derivative(theta, gpm);
    if (i==0) continue;

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
  
  // Set final extrema after the loop. MR April 22 2025
  theta_extrema[0] = constB_ctx->theta_min;
  xp[Z_IDX] = constB_ctx->theta_min;
  gkyl_calc_bmag_global(0.0, xp, &bmag_extrema[0], bmag_ctx);

  theta_extrema[extrema] = constB_ctx->theta_max;
  xp[Z_IDX] = constB_ctx->theta_max;
  gkyl_calc_bmag_global(0.0, xp, &bmag_extrema[extrema], bmag_ctx);
  extrema++;

  gpm->constB_ctx->num_extrema = extrema;
  for (int i = 0; i < extrema; i++)
  {
    gpm->constB_ctx->theta_extrema[i] = theta_extrema[i];
    gpm->constB_ctx->bmag_extrema[i] = bmag_extrema[i];
  }

  // Identify 1 for maxima, 0 for minima

  // Left edge
  if (bmag_extrema[0] > bmag_extrema[1])
  {    gpm->constB_ctx->min_or_max[0] = 1;  } // Maximum
  else if (bmag_extrema[0] < bmag_extrema[1])
  {    gpm->constB_ctx->min_or_max[0] = 0;  } // Minimum
  else
  {    printf("Error: Extrema is not an extrema. Position_map optimization failed\n");  }

  // Middle points
  for (int i = 1; i < extrema - 1; i++)
  {
    if (bmag_extrema[i] > bmag_extrema[i-1] && bmag_extrema[i] > bmag_extrema[i+1])
    {      gpm->constB_ctx->min_or_max[i] = 1;    } // Maximum
    else if (bmag_extrema[i] < bmag_extrema[i-1] && bmag_extrema[i] < bmag_extrema[i+1])
    {      gpm->constB_ctx->min_or_max[i] = 0;    } // Minimum
    else
    {      printf("Error: Extrema is not an extrema. Position_map optimization failed\n");  }
  }

  // Right edge
  if (bmag_extrema[extrema-1] > bmag_extrema[extrema-2])
  {    gpm->constB_ctx->min_or_max[extrema-1] = 1; } // Maximum
  else if (bmag_extrema[extrema-1] < bmag_extrema[extrema-2])
  {    gpm->constB_ctx->min_or_max[extrema-1] = 0; } // Minimum
  else  
  {    printf("Error: Extrema is not an extrema. Position_map optimization failed\n");  }

  // Free mallocs
  gkyl_free(bmag_vals);
  gkyl_free(dbmag_vals);
  gkyl_free(theta_extrema);
  gkyl_free(bmag_extrema);
}

/**
 * Refines the extrema found in the B field along the field line
 * specified by the input psi and alpha values.
 * 
 * @param gpm The position map object
 */
static void
refine_B_field_extrema(struct gkyl_position_map *gpm)
{
  int num_points_per_level = 10; // Number of points to evaluate per level for midpoint rule
  int num_iterations = 22; // Number of iterations to refine the extrema with midpoint rule

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

  for (int i = 1; i < gpm->constB_ctx->num_extrema - 1; i++)
  {
    double theta = gpm->constB_ctx->theta_extrema[i];
    xp[Z_IDX] = theta;
    double bmag_cent, bmag_left, bmag_right;
    gkyl_calc_bmag_global(0.0, xp, &bmag_cent, bmag_ctx);
    xp[Z_IDX] = theta - theta_dxi;
    gkyl_calc_bmag_global(0.0, xp, &bmag_left, bmag_ctx);
    xp[Z_IDX] = theta + theta_dxi;
    gkyl_calc_bmag_global(0.0, xp, &bmag_right, bmag_ctx);

    double interval_left = theta - theta_dxi;
    double interval_right = theta + theta_dxi;
    double extrema_Bmag;
    double extrema_Bmag_location;
    double bmag_out;
    bool is_maximum;
    if (bmag_cent > bmag_left && bmag_cent > bmag_right)
    { is_maximum = true; } // Local maxima
    else if (bmag_cent < bmag_left && bmag_cent < bmag_right)
    { is_maximum = false; } // Local minima
    else
    { printf("Error: Extrema is not an extrema. Position_map optimization failed\n");
      break;
    }

    // Midpoint rule refinement
    for (int j = 0; j < num_iterations; j++)
    {
      double dz = (interval_right - interval_left) / num_points_per_level;

      if (is_maximum)
      { extrema_Bmag = 0.0; }
      else
      { extrema_Bmag = 99999999999999999.; }

      extrema_Bmag_location = 0.0;
      for (int k = 0; k <= num_points_per_level; k++)
      {
        double z = interval_left + k * dz;
        xp[Z_IDX] = z;
        gkyl_calc_bmag_global(0.0, xp, &bmag_out, bmag_ctx);
        if (is_maximum && bmag_out > extrema_Bmag)
        {
          extrema_Bmag = bmag_out;
          extrema_Bmag_location = z;
        }
        else if (!is_maximum && bmag_out < extrema_Bmag)
        {
          extrema_Bmag = bmag_out;
          extrema_Bmag_location = z;
        }
      }
      interval_left = extrema_Bmag_location - dz;
      interval_right = extrema_Bmag_location + dz;
    }
    gpm->constB_ctx->theta_extrema[i] = extrema_Bmag_location;
    gpm->constB_ctx->bmag_extrema[i] = extrema_Bmag;
  }

  // Find the change in B over each cell
  double B_total_change = 0.0; // Total change in magnetic field
  for (int i = 1; i < gpm->constB_ctx->num_extrema; i++)
  {
    B_total_change += fabs(gpm->constB_ctx->bmag_extrema[i] - gpm->constB_ctx->bmag_extrema[i-1]);
  }
  gpm->constB_ctx->dB_cell = B_total_change / (gpm->constB_ctx->N_theta_boundaries);
}

/**
 * Function used for root finding to determine the optimal theta value
 * for the numeric constant dB mapping.
 * 
 * @param theta The theta value to evaluate
 * @param ctx The context for the root finder. Type opt_Theta_ctx
 */
static double
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

/**
 * Maps the uniform computational coordinate to a non-uniform coordinate
 * according to the numeric constant B mapping.
 * 
 * @param t Time
 * @param xn Uniform coordinate
 * @param fout Non-uniform coordinate
 * @param ctx The context for the position map
 */
static void
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
  double dB_cell = gpm->constB_ctx->dB_cell;
  double it = (theta - theta_lo) / theta_dxi;

  // Set strict floor and ceiling limits for theta
  // This is to prevent the root finding algorithm from going out of bounds
  // Not fout[0] = theta because of the finite differences and can lead to jumps
  if (it <= 0 || it >= num_boundaries)
  {
    fout[0] = (it <= 0) ? theta_lo : theta_hi;
    return;
  }

  // Determine which region theta is in
  // Regions start at 0 and count up to num_extrema-1
  // Initial guess is not accurate because the theta_extrema are not Theta_extrema
  // We use itteration to further refine this, but it's a good initial guess
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

  double dB_target, dB_global_lower, B_lower_region;
  double interval_lower, interval_upper, interval_lower_eval, interval_upper_eval;
  struct opt_Theta_ctx ridders_ctx = {
    .gpm = gpm,
    .bmag_ctx = gpm->bmag_ctx,
  };
  dB_target = dB_cell * it;

  bool outside_region = true; // Asuume that we identified the region incorrectly
  while (outside_region)
  {
    dB_global_lower = 0.0;
    for (int i = 0; i < region; i++)
    {
      dB_global_lower += fabs(gpm->constB_ctx->bmag_extrema[i+1] - gpm->constB_ctx->bmag_extrema[i]);
    }
    B_lower_region = gpm->constB_ctx->bmag_extrema[region];

    ridders_ctx.dB_target = dB_target;
    ridders_ctx.dB_global_lower = dB_global_lower;
    ridders_ctx.B_lower_region = B_lower_region;

    interval_lower = theta_extrema[region];
    interval_upper = theta_extrema[region+1];
    interval_lower_eval = position_map_numeric_optimization_function(interval_lower, &ridders_ctx);
    interval_upper_eval = position_map_numeric_optimization_function(interval_upper, &ridders_ctx);

    if (interval_lower_eval * interval_upper_eval < 0) {
      // If the interval changes sign, then there is a zero in between. We can find the root and are in the correct region
      outside_region = false;
    }
    else {
      // It means we are in the wrong region
      if (interval_lower_eval > 0.0 && interval_upper_eval > 0.0) {
        // If the bounds on the interval are both positive, we should move down a region to make it pass through zero
        region--;
        if (region < 0) {
          // If we can't move down any regions and leave the simulation domain, we are likely on the lower limit of the domain and should just return the input theta
          fout[0] = theta_lo;
          return;
        }
      }
      else if (interval_lower_eval < 0.0 && interval_upper_eval < 0.0) {
        // If the bounds on the interval are both negative, we should move up a region to make it pass through zero
        region++;
        if (region > num_extrema-2) {
          // If we can't move up any regions and leave the simulation domain, we are likely on the upper limit of the domain and should just return the input theta
          fout[0] = theta_hi;
          return;
        }
      }
    }
  }

  struct gkyl_qr_res res = gkyl_ridders(position_map_numeric_optimization_function, &ridders_ctx,
    interval_lower, interval_upper, interval_lower_eval, interval_upper_eval, 10, 1e-6);
  double Theta = res.res;
  fout[0] = Theta*gpm->constB_ctx->map_strength + theta*(1-gpm->constB_ctx->map_strength); 

  bool enable_limits_min_B = gpm->constB_ctx->enable_maximum_slope_limits_at_min_B;
  bool enable_limits_max_B = gpm->constB_ctx->enable_maximum_slope_limits_at_max_B;

  if (enable_limits_min_B || enable_limits_max_B)
  {
    // Set a minimum cell size on the edges
    // Assume that at inflection points, Theta = theta. This should be true
    double Theta_left  = interval_lower;
    double Theta_right = interval_upper;
    double theta_middle = 0.5 * (interval_lower + interval_upper);

    bool left_is_maximum = gpm->constB_ctx->min_or_max[region];
    bool right_is_maximum = gpm->constB_ctx->min_or_max[region+1];

    if (theta > theta_middle && left_is_maximum)
    {
      enable_limits_max_B = false;
    }
    if (theta < theta_middle && right_is_maximum)
    {
      enable_limits_max_B = false;
    }

    double max_slope_min_B = gpm->constB_ctx->maximum_slope_at_min_B;
    double max_slope_max_B = gpm->constB_ctx->maximum_slope_at_max_B;

    double right_straight_line_value, left_straight_line_value;
    if (left_is_maximum){
      left_straight_line_value = max_slope_max_B * theta + (1-max_slope_max_B) * Theta_left;
    }
    else {
      left_straight_line_value = max_slope_min_B * theta + (1-max_slope_min_B) * Theta_left;
    }

    if (right_is_maximum){
      right_straight_line_value = max_slope_max_B * theta + (1-max_slope_max_B) * Theta_right;
    }
    else {
      right_straight_line_value = max_slope_min_B * theta + (1-max_slope_min_B) * Theta_right;
    }

    if ( fout[0] < right_straight_line_value && 
      ((right_is_maximum && enable_limits_max_B) || 
      ((!right_is_maximum) && enable_limits_min_B))) 
    {
      fout[0] = right_straight_line_value;
    }

    if (fout[0] > left_straight_line_value && 
      ((left_is_maximum && enable_limits_max_B) ||
      ((!left_is_maximum) && enable_limits_min_B))) 
    {
      fout[0] = left_straight_line_value;
    }
  }
}
