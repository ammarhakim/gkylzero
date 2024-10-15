#include <gkyl_mirror_geo_priv.h>
#include <gkyl_calc_bmag.h>

void
mirror_find_endpoints(struct gkyl_mirror_geo_grid_inp* inp, struct gkyl_mirror_geo *geo, struct arc_length_ctx* arc_ctx, struct plate_ctx* pctx, double psi_curr, double alpha_curr, double* zmin, double* zmax, double* arc_memo){
  enum { PH_IDX, AL_IDX, TH_IDX }; // arrangement of computational coordinates
  enum { X_IDX, Y_IDX, Z_IDX }; // arrangement of cartesian coordinates

  double rclose = inp->rclose;

  double arcL, arcL_curr, arcL_lo;
  double arcL_l, arcL_r;
  double phi_r, phi_l;

  if (geo->plate_spec){ // if we dont have a fixed zmin and zmax
    double rzplate[2];
    pctx->psi_curr = psi_curr;
    pctx->lower=false;
    double a = 0;
    double b = 1;
    double fa = mirror_plate_psi_func(a, pctx);
    double fb = mirror_plate_psi_func(b, pctx);
    struct gkyl_qr_res res = gkyl_ridders(mirror_plate_psi_func, pctx,
      a, b, fa, fb, geo->root_param.max_iter, 1e-10);
    double smax = res.res;
    geo->plate_func_upper(smax, rzplate);
    *zmax = rzplate[1];

    pctx->lower=true;
    a = 0;
    b = 1;
    fa = mirror_plate_psi_func(a, pctx);
    fb = mirror_plate_psi_func(b, pctx);
    res = gkyl_ridders(mirror_plate_psi_func, pctx,
      a, b, fa, fb, geo->root_param.max_iter, 1e-10);
    double smin = res.res;
    geo->plate_func_lower(smin, rzplate);
    *zmin = rzplate[1];
  }

  arcL = integrate_psi_contour_memo(geo, psi_curr, *zmin, *zmax, rclose, true, true, arc_memo);
  arc_ctx->arcL_tot = arcL;
}



void
mirror_set_ridders(struct gkyl_mirror_geo_grid_inp* inp, struct arc_length_ctx* arc_ctx, double psi_curr, double arcL, double arcL_curr, double zmin, double zmax, double* rclose, double *ridders_min, double* ridders_max){

  *ridders_min = -arcL_curr;
  *ridders_max = arcL-arcL_curr;
  arc_ctx->zmin = zmin;
  arc_ctx->zmax = zmax;

  arc_ctx->psi = psi_curr;
  arc_ctx->rclose = *rclose;
  arc_ctx->arcL = arcL_curr;
}



double map_theta_to_z(const double uniform_coordinate, void *arc_ctx)
{
  // Converts our uniform coordinate along field line length to a non-uniform coordinate
  // according to the polynomial mapping of arbitrary order.
  // Notation: I switch from theta to z here. Both are computational coordinate.
  struct arc_length_ctx *arc_app = arc_ctx;
  int n_ex = arc_app->mapping_order_expander;//app->mapping_order_expander;
  int n_ct = arc_app->mapping_order_center;//app->mapping_order_center;
  double z_min = arc_app->theta_min * 0.95;
  double z_max = arc_app->theta_max * 0.95;
  double z_m = arc_app->theta_throat;
  double frac = arc_app->mapping_frac; // 1 is full mapping, 0 is no mapping
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
  return nonuniform_coordinate;
}

void
calculate_mirror_throat_location(void *arc_ctx, void *bmag_ctx_inp)
{
  enum { X_IDX, Y_IDX, Z_IDX }; // arrangement of cartesian coordinates
  struct arc_length_ctx *app = arc_ctx;
  double psi = app->psi;
  double alpha = app->alpha;
  double *xp = malloc(3*sizeof(double));
  xp[X_IDX] = psi;
  xp[Y_IDX] = alpha;
  int itterations = 10;
  double interval_left = 0.0;
  double interval_right = app->theta_max;
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
      gkyl_calc_bmag_global(0.0, xp, fout, bmag_ctx_inp);
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
  app->theta_throat = maximum_Bmag_location;
  app->Bmag_throat = maximum_Bmag;
  free(xp);
  free(fout);
}

void
calculate_optimal_mapping(void *arc_ctx, void *bmag_ctx_inp)
{
  // Could be refined further by doing midpoint root finding, like the mirror throat finding does
  // Expander region
  enum { X_IDX, Y_IDX, Z_IDX }; // arrangement of cartesian coordinates
  struct arc_length_ctx *arc_app = arc_ctx;
  double psi = arc_app->psi;
  double alpha = arc_app->alpha;
  double *xp = malloc(3*sizeof(double));
  double *fout = malloc(3*sizeof(double));
  xp[X_IDX] = psi;
  xp[Y_IDX] = alpha;
  arc_app->mapping_order_center = 1;
  double scan_cells = 50;
  double scan_left = arc_app->theta_throat;
  double scan_right = arc_app->theta_max;
  double scan_dxi = (scan_right - scan_left) / scan_cells;
  int expander_order = 1;
  double max_dB_dCell_prior = 99999999.99;
  double max_dB_dCell;
  double max_dB_dCell_order1 = 0.0;
  while (1)
  {
    max_dB_dCell = 0.0;
    arc_app->mapping_order_expander = expander_order;
    for (int iz = 0; iz < scan_cells; iz++)
    {
      double left_xi = scan_left + iz * scan_dxi;
      double right_xi = scan_left + (iz + 1) * scan_dxi;
      double psi = arc_app->psi;
      double alpha = arc_app->alpha;
      double left_theta = map_theta_to_z(left_xi, arc_ctx);
      double right_theta = map_theta_to_z(right_xi, arc_ctx);
      xp[Z_IDX] = left_theta;
      gkyl_calc_bmag_global(0.0, xp, fout, bmag_ctx_inp);
      double Bmag_left = fout[0];
      xp[Z_IDX] = right_theta;
      gkyl_calc_bmag_global(0.0, xp, fout, bmag_ctx_inp);
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
      arc_app->mapping_order_expander = expander_order;
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
  scan_right = arc_app->theta_throat;
  scan_dxi = (scan_right - scan_left) / scan_cells;
  int center_order = 1;
  max_dB_dCell_prior = 99999999.99;
  while (1)
  {
    max_dB_dCell = 0.0;
    arc_app->mapping_order_center = center_order;
    for (int iz = 0; iz < scan_cells; iz++)
    {
      double left_xi = scan_left + iz * scan_dxi;
      double right_xi = scan_left + (iz + 1) * scan_dxi;
      double left_theta = map_theta_to_z(left_xi, arc_ctx);
      double right_theta = map_theta_to_z(right_xi, arc_ctx);
      xp[Z_IDX] = left_theta;
      gkyl_calc_bmag_global(0.0, xp, fout, bmag_ctx_inp);
      double Bmag_left = fout[0];
      xp[Z_IDX] = right_theta;
      gkyl_calc_bmag_global(0.0, xp, fout, bmag_ctx_inp);
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
      arc_app->mapping_order_center = center_order;
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