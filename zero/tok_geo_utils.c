#include <gkyl_tok_geo_priv.h>

void
find_endpoints(struct gkyl_tok_geo_grid_inp* inp, struct gkyl_tok_geo *geo, struct arc_length_ctx* arc_ctx, struct plate_ctx* pctx, double psi_curr, double alpha_curr, double* zmin, double* zmax, double* zmin_left, double* zmin_right, double* arc_memo, double* arc_memo_left, double* arc_memo_right){
  enum { PH_IDX, AL_IDX, TH_IDX }; // arrangement of computational coordinates
  enum { X_IDX, Y_IDX, Z_IDX }; // arrangement of cartesian coordinates

  double rclose = inp->rclose;
  double rright = inp->rright;
  double rleft = inp->rleft;

  double arcL, arcL_curr, arcL_lo;
  double arcL_l, arcL_r;
  double phi_r, phi_l;

  if(inp->ftype == GKYL_CORE){
    //Find the turning points
    double zlo, zup, zlo_last;
    zlo = 0.01;
    zup=*zmax;
    zlo_last = zlo;
    double R[4], dR[4];
    while(true){
      int nlo = R_psiZ(geo, psi_curr, zlo, 4, R, dR);
      if(nlo==2){
        if(fabs(zlo-zup)<1e-12){
          *zmax = zlo;
          break;
        }
        zlo_last = zlo;
        zlo = (zlo+zup)/2;
      }
      if(nlo==0){
        zup = zlo;
        zlo = zlo_last;
      }
    }
    for(int i =0; i<4; i++){
      R[i] = 0.0;
      dR[i] = 0.0;
    }
    int nup = 0;
    double zup_last;
    zup = -0.01;
    zlo=*zmin;
    zup_last = zup;
    while(true){
      int nup = R_psiZ(geo, psi_curr, zup, 4, R, dR);
      if(nup==2){
        if(fabs(zlo-zup)<1e-12){
          *zmin = zup;
          break;
        }
        zup_last = zup;
        zup = (zlo+zup)/2;
      }
      if(nup==0){
        zlo = zup;
        zup = zup_last;
      }
    }
    // Done finding turning points
    arcL_r = integrate_psi_contour_memo(geo, psi_curr, *zmin, *zmax, rright,
      true, true, arc_memo_right);
    arc_ctx->arcL_right = arcL_r;
    arc_ctx->right = false;
    arcL_l = integrate_psi_contour_memo(geo, psi_curr, *zmin, *zmax, rleft,
      true, true, arc_memo_left);
    arcL = arcL_l + arcL_r;
    arc_ctx->arcL_tot = arcL;

    arc_ctx->right = true;
    arc_ctx->phi_right = 0.0;
    arc_ctx->rclose = rright;
    arc_ctx->psi = psi_curr;
    phi_r = phi_func(alpha_curr, *zmax, arc_ctx);
    arc_ctx->phi_right = phi_r - alpha_curr; // otherwise alpha will get added on twice
  }
  else if(inp->ftype == GKYL_PF_LO){
    //Find the  upper turning point
    double zlo, zup, zlo_last;
    zlo = fmax(inp->zmin_left, inp->zmin_right);
    zup=*zmax;
    zlo_last = zlo;
    double R[4], dR[4];
    while(true){
      int nlo = R_psiZ(geo, psi_curr, zlo, 4, R, dR);
      if(nlo>=2){
        if(fabs(zlo-zup)<1e-12){
          *zmax = zlo;
          break;
        }
        zlo_last = zlo;
        zlo = (zlo+zup)/2;
      }
      if(nlo<2){
        zup = zlo;
        zlo = zlo_last;
      }
    }
    // Done finding turning point
    arcL_r = integrate_psi_contour_memo(geo, psi_curr, inp->zmin_right, *zmax, rright,
      true, true, arc_memo);
    arc_ctx->arcL_right = arcL_r;
    arc_ctx->right = false;
    arcL_l = integrate_psi_contour_memo(geo, psi_curr, inp->zmin_left, *zmax, rleft,
      true, true, arc_memo);
    arcL = arcL_l + arcL_r;
    arc_ctx->arcL_tot = arcL;

    arc_ctx->right = true;
    arc_ctx->phi_right = 0.0;
    arc_ctx->rclose = rright;
    arc_ctx->psi = psi_curr;
    arc_ctx->zmin = inp->zmin_right;
    phi_r = phi_func(alpha_curr, *zmax, arc_ctx);
    arc_ctx->phi_right = phi_r - alpha_curr; // otherwise alpha will get added on twice
  }
  else if(inp->ftype == GKYL_PF_UP){
    //Find the lower turning point
    double zlo, zup, zlo_last;
    double zup_last;
    zup = fmin(inp->zmax_left, inp->zmax_right);
    zlo=*zmin;
    zup_last = zup;
    double R[4], dR[4];
    while(true){
      int nup = R_psiZ(geo, psi_curr, zup, 4, R, dR);
      if(nup>=2){
        if(fabs(zlo-zup)<1e-12){
          *zmin = zup;
          break;
        }
        zup_last = zup;
        zup = (zlo+zup)/2;
      }
      if(nup<2){
        zlo = zup;
        zup = zup_last;
      }
    }
    // Done finding turning point
    arcL_r = integrate_psi_contour_memo(geo, psi_curr, *zmin, inp->zmax_right, rright,
      true, true, arc_memo);
    arc_ctx->right = false;
    arcL_l = integrate_psi_contour_memo(geo, psi_curr, *zmin, inp->zmax_left, rleft,
      true, true, arc_memo);
    arc_ctx->arcL_left= arcL_l;
    arcL = arcL_l + arcL_r;
    arc_ctx->arcL_tot = arcL;

    arc_ctx->right = false;
    arc_ctx->phi_left = 0.0;
    arc_ctx->rclose = rleft;
    arc_ctx->psi = psi_curr;
    arc_ctx->zmax = inp->zmax_left;
    phi_l = phi_func(alpha_curr, *zmin, arc_ctx);
    arc_ctx->phi_left = phi_l - alpha_curr; // otherwise alpha will get added on twice
  }
  else if(inp->ftype==GKYL_SOL_DN_OUT){
    if (geo->plate_spec){ // if we dont have a fixed zmin and zmax
      double rzplate[2];
      pctx->psi_curr = psi_curr;
      pctx->lower=false;
      double a = 0;
      double b = 1;
      double fa = plate_psi_func(a, pctx);
      double fb = plate_psi_func(b, pctx);
      struct gkyl_qr_res res = gkyl_ridders(plate_psi_func, pctx,
        a, b, fa, fb, geo->root_param.max_iter, 1e-10);
      double smax = res.res;
      geo->plate_func_upper(smax, rzplate);
      *zmax = rzplate[1];

      pctx->lower=true;
      a = 0;
      b = 1;
      fa = plate_psi_func(a, pctx);
      fb = plate_psi_func(b, pctx);
      res = gkyl_ridders(plate_psi_func, pctx,
        a, b, fa, fb, geo->root_param.max_iter, 1e-10);
      double smin = res.res;
      geo->plate_func_lower(smin, rzplate);
      *zmin = rzplate[1];
    }

    arc_ctx->phi_right = 0.0;
    arcL = integrate_psi_contour_memo(geo, psi_curr, *zmin, *zmax, rclose, true, true, arc_memo);
    arc_ctx->arcL_tot = arcL;
  }
  else if(inp->ftype==GKYL_SOL_DN_IN){
    arc_ctx->phi_right = 0.0;
    arcL = integrate_psi_contour_memo(geo, psi_curr, *zmin, *zmax, rclose, true, true, arc_memo);
    arc_ctx->arcL_tot = arcL;
  }
  else if(inp->ftype == GKYL_SOL_SN_LO){
    //Find the  upper turning point
    double zlo, zup, zlo_last;
    zlo = fmax(inp->zmin_left, inp->zmin_right);
    zup = *zmax;
    zlo_last = zlo;
    double R[4], dR[4];
    while(true){
      int nlo = R_psiZ(geo, psi_curr, zlo, 4, R, dR);
      if(nlo>=2){
        if(fabs(zlo-zup)<1e-12){
          *zmax = zlo;
          break;
        }
        zlo_last = zlo;
        zlo = (zlo+zup)/2;
      }
      if(nlo<2){
        zup = zlo;
        zlo = zlo_last;
      }
    }
    if (geo->plate_spec){ // if we dont have a fixed zmin, set based on plate func
      double rzplate[2];
      pctx->psi_curr = psi_curr;
      pctx->lower=false;
      double a = 0;
      double b = 1;
      double fa = plate_psi_func(a, pctx);
      double fb = plate_psi_func(b, pctx);
      struct gkyl_qr_res res = gkyl_ridders(plate_psi_func, pctx,
        a, b, fa, fb, geo->root_param.max_iter, 1e-10);
      double smax = res.res;
      geo->plate_func_upper(smax, rzplate);
      *zmin_left = rzplate[1];

      pctx->lower=true;
      a = 0;
      b = 1;
      fa = plate_psi_func(a, pctx);
      fb = plate_psi_func(b, pctx);
      res = gkyl_ridders(plate_psi_func, pctx,
        a, b, fa, fb, geo->root_param.max_iter, 1e-10);
      double smin = res.res;
      geo->plate_func_lower(smin, rzplate);
      *zmin_right = rzplate[1];
    }

    // Done finding turning point
    arcL_r = integrate_psi_contour_memo(geo, psi_curr, *zmin_right, *zmax, rright,
      true, true, arc_memo);
    arc_ctx->arcL_right = arcL_r;
    arc_ctx->right = false;
    arcL_l = integrate_psi_contour_memo(geo, psi_curr, *zmin_left, *zmax, rleft,
      true, true, arc_memo);
    arcL = arcL_l + arcL_r;
    arc_ctx->arcL_tot = arcL;

    arc_ctx->right = true;
    arc_ctx->phi_right = 0.0;
    arc_ctx->rclose = rright;
    arc_ctx->psi = psi_curr;
    arc_ctx->zmin = inp->zmin_right;
    phi_r = phi_func(alpha_curr, *zmax, arc_ctx);
    arc_ctx->phi_right = phi_r - alpha_curr; // otherwise alpha will get added on twice
  }


}



void
set_ridders(struct gkyl_tok_geo_grid_inp* inp, struct arc_length_ctx* arc_ctx, double psi_curr, double arcL, double arcL_curr, double zmin, double zmax, double zmin_left, double zmin_right, double rright, double rleft, double* rclose, double *ridders_min, double* ridders_max){


  if(inp->ftype==GKYL_CORE){
    if(arcL_curr <= arc_ctx->arcL_right){
      *rclose = rright;
      arc_ctx->right = true;
      *ridders_min = -arcL_curr;
      *ridders_max = arcL-arcL_curr;
      arc_ctx->zmin = zmin;
      arc_ctx->zmax = zmax;
    }
    else{
      *rclose = rleft;
      arc_ctx->right = false;
      *ridders_min = arcL - arcL_curr;
      *ridders_max = -arcL_curr + arc_ctx->arcL_right;
      arc_ctx->zmin = zmin;
      arc_ctx->zmax = zmax;
    }
  }
  if(inp->ftype==GKYL_PF_LO){
    if(arcL_curr <= arc_ctx->arcL_right){
      *rclose = rright;
      arc_ctx->right = true;
      *ridders_min = -arcL_curr;
      *ridders_max = arcL-arcL_curr;
      arc_ctx->zmin = inp->zmin_right;
      arc_ctx->zmax = zmax;
    }
    else{
      *rclose = rleft;
      arc_ctx->right = false;
      *ridders_min = arcL - arcL_curr;
      *ridders_max = -arcL_curr + arc_ctx->arcL_right;
      arc_ctx->zmin = inp->zmin_left;
      arc_ctx->zmax = zmax;
    }
  }
  if(inp->ftype==GKYL_PF_UP){
    if(arcL_curr > arc_ctx->arcL_left){
      *rclose = rright;
      arc_ctx->right = true;
      *ridders_min = arc_ctx->arcL_left - arcL_curr;
      *ridders_max = arcL - arcL_curr;
      arc_ctx->zmin = zmin;
      arc_ctx->zmax = inp->zmax_right;
    }
    else{
      *rclose = rleft;
      arc_ctx->right = false;
      *ridders_min = arc_ctx->arcL_left - arcL_curr;
      *ridders_max = -arcL_curr;
      arc_ctx->zmin = zmin;
      arc_ctx->zmax = inp->zmax_left;
    }
  }
  if(arc_ctx->ftype==GKYL_SOL_DN_OUT){
    *ridders_min = -arcL_curr;
    *ridders_max = arcL-arcL_curr;
    arc_ctx->right = false;
    arc_ctx->zmin = zmin;
    arc_ctx->zmax = zmax;
  }
  if(arc_ctx->ftype==GKYL_SOL_DN_IN){
    *ridders_min = arcL-arcL_curr;
    *ridders_max = -arcL_curr;
    arc_ctx->right = false;
    arc_ctx->zmin = zmin;
    arc_ctx->zmax = zmax;
  }
  if(arc_ctx->ftype==GKYL_SOL_SN_LO){
    if(arcL_curr <= arc_ctx->arcL_right){
      *rclose = rright;
      arc_ctx->right = true;
      *ridders_min = -arcL_curr;
      *ridders_max = arcL-arcL_curr;
      arc_ctx->zmin = zmin_right;
      arc_ctx->zmax = zmax;
    }
    else{
      *rclose = rleft;
      arc_ctx->right = false;
      *ridders_min = arcL - arcL_curr;
      *ridders_max = -arcL_curr + arc_ctx->arcL_right;
      arc_ctx->zmin = zmin_left;
      arc_ctx->zmax = zmax;
    }
  }

  arc_ctx->psi = psi_curr;
  arc_ctx->rclose = *rclose;
  arc_ctx->arcL = arcL_curr;
}
