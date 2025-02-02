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
