#include <gkyl_mirror_geo_priv.h>

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

double map_theta_to_z(const double uniform_coordinate, const double z_min, const double z_max, const double z_m)
{
  int n_ex = 2;//app->mapping_order_expander;
  int n_ct = 2;//app->mapping_order_center;
  int n;
  double frac = 1;//app->mapping_frac; // 1 is full mapping, 0 is no mapping
  double nonuniform_coordinate, left, right;
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