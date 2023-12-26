#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_rio.h>
#include <gkyl_basis.h>
#include <gkyl_tok_geo.h>
#include <gkyl_math.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_nodal_ops.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_efit.h>
#include <gkyl_tok_geo_priv.h>

#include <math.h>
#include <string.h>



// Function context to pass to root finder
struct arc_length_ctx {
  const gkyl_tok_geo *geo;
  double *arc_memo;
  double *arc_memo_left;
  double *arc_memo_right;
  double psi, rclose, zmin, arcL;
  double rleft, rright, zmax;
  double arcL_right; // this is for when we need to switch sides
  double arcL_left; // this is for when we need to switch sides
  double phi_right; // this is for when we need to switch sides
  double phi_left; // this is for when we need to switch sides
  double phi_bot; // For new way of trying to do core
  bool right;
  double zmaxis;
  enum gkyl_tok_geo_type ftype; // type of geometry
};
struct plate_ctx{
  const struct gkyl_tok_geo* geo;
  double psi_curr;
  bool lower;
};

static inline double
plate_psi_func(double s, void *ctx){
  // uses a pointer to the plate function to get R(s), Z(s)
  // Then calculates psi(R, Z)
  // will be used by ridders later
  
  struct plate_ctx *gc = ctx;
  double RZ[2];
  if(gc->lower==true)
    gc->geo->plate_func_lower(s, RZ);
  else
    gc->geo->plate_func_upper(s, RZ);

  double R = RZ[0];
  double Z = RZ[1];

  // Now find the cell where this R and Z is
  int rzidx[2];
  rzidx[0] = fmin(gc->geo->rzlocal.lower[0] + (int) floor((R - gc->geo->rzgrid.lower[0])/gc->geo->rzgrid.dx[0]), gc->geo->rzlocal.upper[0]);
  rzidx[1] = fmin(gc->geo->rzlocal.lower[1] + (int) floor((Z - gc->geo->rzgrid.lower[1])/gc->geo->rzgrid.dx[1]), gc->geo->rzlocal.upper[1]);
  long loc = gkyl_range_idx(&gc->geo->rzlocal, rzidx);
  const double *coeffs = gkyl_array_cfetch(gc->geo->psiRZ,loc);

  double xc[2];
  gkyl_rect_grid_cell_center(&gc->geo->rzgrid, rzidx, xc);
  double xy[2];
  xy[0] = (R-xc[0])/(gc->geo->rzgrid.dx[0]*0.5);
  xy[1] = (Z-xc[1])/(gc->geo->rzgrid.dx[1]*0.5);
  double psi = gc->geo->rzbasis.eval_expand(xy, coeffs);
  return psi - gc->psi_curr;
}



// Function to pass to root-finder to find Z location for given arc-length
static inline double
arc_length_func(double Z, void *ctx)
{
  struct arc_length_ctx *actx = ctx;
  double *arc_memo;
  double psi = actx->psi, rclose = actx->rclose, zmin = actx->zmin, arcL = actx->arcL;
  double zmax = actx->zmax;
  double ival = 0.0;

  if(actx->ftype==GKYL_SOL_DN_OUT){
    double *arc_memo = actx->arc_memo;
    ival = integrate_psi_contour_memo(actx->geo, psi, zmin, Z, rclose, true, false, arc_memo) - arcL;
  }
  if(actx->ftype==GKYL_SOL_DN_IN){
    ival = integrate_psi_contour_memo(actx->geo, psi, Z, zmax, rclose, true, false, arc_memo) - arcL;
  }
  else if(actx->ftype==GKYL_SOL_SN_LO){
    if(actx->right==true){
      ival = integrate_psi_contour_memo(actx->geo, psi, zmin, Z, rclose, false, false, arc_memo) - arcL;
    }
    else{
      ival = integrate_psi_contour_memo(actx->geo, psi, Z, zmax, rclose, false, false, arc_memo)  - arcL + actx->arcL_right;
    }
  }
  else if(actx->ftype==GKYL_CORE){
    if(actx->right==true){
      double *arc_memo = actx->arc_memo_right;
      ival = integrate_psi_contour_memo(actx->geo, psi, zmin, Z, rclose, true, false, arc_memo) - arcL;
    }
    else{
      double *arc_memo = actx->arc_memo_left;
      ival = integrate_psi_contour_memo(actx->geo, psi, Z, zmax, rclose, true, false, arc_memo)  - arcL + actx->arcL_right;
    }
  }
  else if(actx->ftype==GKYL_PF_LO){
    if(actx->right==true){
      ival = integrate_psi_contour_memo(actx->geo, psi, zmin, Z, rclose, false, false, arc_memo) - arcL;
    }
    else{
      ival = integrate_psi_contour_memo(actx->geo, psi, Z, zmax, rclose, false, false, arc_memo)  - arcL + actx->arcL_right;
    }
  }
  else if(actx->ftype==GKYL_PF_UP){ // not yet implemented
    if(actx->right==true){
      ival = integrate_psi_contour_memo(actx->geo, psi, zmin, Z, rclose, false, false, arc_memo) + actx->arcL_left - arcL;
    }
    else{
      ival = integrate_psi_contour_memo(actx->geo, psi, Z, zmax, rclose, false, false, arc_memo) - arcL;
    }
  }

  return ival;
}

// Function to calculate phi given alpha
static inline double
phi_func(double alpha_curr, double Z, void *ctx)
{
  struct arc_length_ctx *actx = ctx;
  double *arc_memo = actx->arc_memo;
  double psi = actx->psi, rclose = actx->rclose, zmin = actx->zmin, arcL = actx->arcL;

  double ival = 0;
  double phi_ref = 0.0;
  if(actx->ftype==GKYL_SOL_DN_OUT){ // Using convention from Noah Mandell's thesis Eq 5.104 phi = alpha at midplane
    if (Z<actx->zmaxis)
      ival = -integrate_phi_along_psi_contour_memo(actx->geo, psi, Z, actx->zmaxis, rclose, false, false, arc_memo);
    else
      ival = integrate_phi_along_psi_contour_memo(actx->geo, psi, actx->zmaxis, Z, rclose, false, false, arc_memo);
  }
  if(actx->ftype==GKYL_SOL_DN_IN){ // alpha = phi at inboard midplane
    if (Z<actx->zmaxis)
      ival = integrate_phi_along_psi_contour_memo(actx->geo, psi, Z, actx->zmaxis, rclose, false, false, arc_memo);
    else
      ival = -integrate_phi_along_psi_contour_memo(actx->geo, psi, actx->zmaxis, Z, rclose, false, false, arc_memo);
  }
  else if(actx->ftype==GKYL_SOL_SN_LO){ // alpha = phi at outboard midplane
    if (actx->right==true){
      if (Z<actx->zmaxis)
        ival = -integrate_phi_along_psi_contour_memo(actx->geo, psi, Z, actx->zmaxis, rclose, false, false, arc_memo);
      else
        ival = integrate_phi_along_psi_contour_memo(actx->geo, psi, actx->zmaxis, Z, rclose, false, false, arc_memo);
    }
    else{
      ival = integrate_phi_along_psi_contour_memo(actx->geo, psi, Z, actx->zmax, rclose, false, false, arc_memo);
      phi_ref = actx->phi_right;
    }
  }
  else if (actx->ftype==GKYL_CORE){ 
    if(actx->right==true){ // phi = alpha at outboard midplane
      if(Z<actx->zmaxis)
        ival = -integrate_phi_along_psi_contour_memo(actx->geo, psi, Z, actx->zmaxis, rclose, false, false, arc_memo);
      else
        ival = integrate_phi_along_psi_contour_memo(actx->geo, psi, actx->zmaxis, Z, rclose, false, false, arc_memo);
    }
    else{// alpha = phi at inboard midplane
      if (Z<actx->zmaxis)
        ival = integrate_phi_along_psi_contour_memo(actx->geo, psi, Z, actx->zmaxis, rclose, false, false, arc_memo);
      else
        ival = -integrate_phi_along_psi_contour_memo(actx->geo, psi, actx->zmaxis, Z, rclose, false, false, arc_memo);
    }
  }
  else if(actx->ftype==GKYL_PF_LO){
    if(actx->right==true){
      ival = integrate_phi_along_psi_contour_memo(actx->geo, psi, actx->zmin, Z, rclose, false, false, arc_memo);
    }
    else{
      ival = -integrate_phi_along_psi_contour_memo(actx->geo, psi, Z, actx->zmax, rclose, false, false, arc_memo);// + actx->phi_right;
      phi_ref = actx->phi_right;
    }
  }
  else if(actx->ftype==GKYL_PF_UP){
    if(actx->right==true){
      ival = integrate_phi_along_psi_contour_memo(actx->geo, psi, actx->zmin, Z, rclose, false, false, arc_memo);
      phi_ref = actx->phi_left;
    }
    else{
      ival = integrate_phi_along_psi_contour_memo(actx->geo, psi, Z, actx->zmax, rclose, false, false, arc_memo);// + actx->phi_right;
    }
  }
  // Now multiply by fpol
  double R[4] = {0};
  double dR[4] = {0};
  int nr = R_psiZ(actx->geo, psi, Z, 4, R, dR);
  double r_curr = nr == 1 ? R[0] : choose_closest(rclose, R, R, nr);
  double psi_fpol = psi;
  if (psi_fpol < actx->geo->psisep) // F = F(psi_sep) in the SOL. Convention of psi increases inward
    psi_fpol = actx->geo->psisep;
  int idx = fmin(actx->geo->frange.lower[0] + (int) floor((psi_fpol - actx->geo->fgrid.lower[0])/actx->geo->fgrid.dx[0]), actx->geo->frange.upper[0]);
  long loc = gkyl_range_idx(&actx->geo->frange, &idx);
  const double *coeffs = gkyl_array_cfetch(actx->geo->fpoldg,loc);
  double fxc;
  gkyl_rect_grid_cell_center(&actx->geo->fgrid, &idx, &fxc);
  double fx = (psi_fpol-fxc)/(actx->geo->fgrid.dx[0]*0.5);
  double fpol = actx->geo->fbasis.eval_expand(&fx, coeffs);
  ival = ival*fpol;

  while(ival < -M_PI){
    ival +=2*M_PI;
  }
  while(ival > M_PI){
    ival -=2*M_PI;
  }
  return alpha_curr + ival + phi_ref;
}



gkyl_tok_geo*
gkyl_tok_geo_new(const struct gkyl_tok_geo_inp *inp)
{
  struct gkyl_tok_geo *geo = gkyl_malloc(sizeof(*geo));

  geo->efit = gkyl_efit_new(inp->filepath, inp->rzpoly_order, inp->fluxpoly_order, false);

  geo->plate_spec = inp->plate_spec;
  geo->plate_func_lower = inp->plate_func_lower;
  geo->plate_func_upper = inp->plate_func_upper;

  geo->rzbasis= *geo->efit->rzbasis;
  geo->rzgrid = *geo->efit->rzgrid;
  geo->psiRZ = gkyl_array_acquire(geo->efit->psizr);
  geo->psibyrRZ = gkyl_array_acquire(geo->efit->psibyrzr);
  geo->psibyr2RZ = gkyl_array_acquire(geo->efit->psibyr2zr);

  geo->num_rzbasis = geo->efit->rzbasis->num_basis;
  geo->rzlocal = *geo->efit->rzlocal;
  geo->rzlocal_ext = *geo->efit->rzlocal_ext;
  geo->fgrid = *geo->efit->fluxgrid;
  geo->fbasis = *geo->efit->fluxbasis;
  geo->frange = *geo->efit->fluxlocal;
  geo->frange_ext = *geo->efit->fluxlocal_ext;
  geo->fpoldg= gkyl_array_acquire(geo->efit->fpolflux);
  geo->qdg= gkyl_array_acquire(geo->efit->qflux);
  geo->psisep = geo->efit->sibry;
  geo->zmaxis = geo->efit->zmaxis;

  geo->root_param.eps =
    inp->root_param.eps > 0 ? inp->root_param.eps : 1e-10;
  geo->root_param.max_iter =
    inp->root_param.max_iter > 0 ? inp->root_param.max_iter : 100;

  geo->quad_param.max_level =
    inp->quad_param.max_levels > 0 ? inp->quad_param.max_levels : 10;
  geo->quad_param.eps =
    inp->quad_param.eps > 0 ? inp->quad_param.eps : 1e-10;

  if (geo->efit->rzbasis->poly_order == 1)
    geo->calc_roots = calc_RdR_p1;
  else if (geo->efit->rzbasis->poly_order == 2)
    geo->calc_roots = calc_RdR_p2;

  geo->stat = (struct gkyl_tok_geo_stat) { };

  
  return geo;
}

double
gkyl_tok_geo_integrate_psi_contour(const gkyl_tok_geo *geo, double psi,
  double zmin, double zmax, double rclose)
{
  return integrate_psi_contour_memo(geo, psi, zmin, zmax, rclose,
    false, false, 0);
}

int
gkyl_tok_geo_R_psiZ(const gkyl_tok_geo *geo, double psi, double Z, int nmaxroots,
  double *R, double *dR)
{
  return R_psiZ(geo, psi, Z, nmaxroots, R, dR);
}

// write out nodal coordinates 
static void
write_nodal_coordinates(const char *nm, struct gkyl_range *nrange,
  struct gkyl_array *nodes)
{
  double lower[3] = { 0.0, 0.0, 0.0 };
  double upper[3] = { 1.0, 1.0, 1.0 };
  int cells[3];
  for (int i=0; i<nrange->ndim; ++i)
    cells[i] = gkyl_range_shape(nrange, i);
  
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 3, lower, upper, cells);

  gkyl_grid_sub_array_write(&grid, nrange, nodes, nm);
}


void gkyl_tok_geo_calc(struct gk_geometry* up, struct gkyl_range *nrange, double dzc[3], 
  evalf_t mapc2p_func, void* mapc2p_ctx, evalf_t bmag_func, void *bmag_ctx, 
  struct gkyl_array *mc2p_nodal_fd, struct gkyl_array *mc2p_nodal, struct gkyl_array *mc2p)
{

  struct gkyl_tok_geo *geo = mapc2p_ctx;
  struct gkyl_tok_geo_geo_inp *inp = bmag_ctx;

  enum { PH_IDX, AL_IDX, TH_IDX }; // arrangement of computational coordinates
  enum { X_IDX, Y_IDX, Z_IDX }; // arrangement of cartesian coordinates
  
  double dtheta = inp->cgrid.dx[TH_IDX],
    dpsi = inp->cgrid.dx[PH_IDX],
    dalpha = inp->cgrid.dx[AL_IDX];
  
  double theta_lo = inp->cgrid.lower[TH_IDX],
    phi_lo = inp->cgrid.lower[PH_IDX],
    alpha_lo = inp->cgrid.lower[AL_IDX];

  double dx_fact = up->basis.poly_order == 1 ? 1 : 0.5;
  dtheta *= dx_fact; dpsi *= dx_fact; dalpha *= dx_fact;

  // used for finite differences 
  double delta_alpha = dalpha*1e-4;
  double delta_psi = dpsi*1e-6;
  double delta_theta = dtheta*1e-4;
  dzc[0] = delta_psi;
  dzc[1] = delta_alpha;
  dzc[2] = delta_theta;
  int modifiers[5] = {0, -1, 1, -2, 2};

  double rclose = inp->rclose;
  double rright = inp->rright;
  double rleft = inp->rleft;

  int nzcells = geo->rzgrid.cells[1];
  double *arc_memo = gkyl_malloc(sizeof(double[nzcells]));
  double *arc_memo_left = gkyl_malloc(sizeof(double[nzcells]));
  double *arc_memo_right = gkyl_malloc(sizeof(double[nzcells]));

  struct arc_length_ctx arc_ctx = {
    .geo = geo,
    .arc_memo = arc_memo,
    .arc_memo_right = arc_memo_right,
    .arc_memo_left = arc_memo_left,
    .ftype = inp->ftype,
    .zmaxis = geo->zmaxis
  };
  struct plate_ctx pctx = {
    .geo = geo
  };

  int cidx[3] = { 0 };
  for(int ia=nrange->lower[AL_IDX]; ia<=nrange->upper[AL_IDX]; ++ia){
    cidx[AL_IDX] = ia;
    for(int ia_delta = 0; ia_delta < 5; ia_delta++){ // should be <5
      if(ia == nrange->lower[AL_IDX]){
        if(ia_delta == 1 || ia_delta == 3)
          continue; // want to use one sided stencils at edge
      }
      else if(ia == nrange->upper[AL_IDX]){
          if(ia_delta == 2 || ia_delta == 4)
            continue; // want to use one sided stencils at edge
      }
      else{ //interior
        if( ia_delta == 3 || ia_delta == 4)
          continue; //dont do two away
      }

      double alpha_curr = alpha_lo + ia*dalpha + modifiers[ia_delta]*delta_alpha;

      for (int ip=nrange->lower[PH_IDX]; ip<=nrange->upper[PH_IDX]; ++ip) {
        int ip_delta_max = 5;// should be 5
        if(ia_delta != 0)
          ip_delta_max = 1;
        for(int ip_delta = 0; ip_delta < ip_delta_max; ip_delta++){
          if(ip == nrange->lower[PH_IDX]){
            if(ip_delta == 1 || ip_delta == 3)
              continue; // want to use one sided stencils at edge
          }
          else if(ip == nrange->upper[PH_IDX]){
            if(ip_delta == 2 || ip_delta == 4)
              continue; // want to use one sided stencils at edge
          }
          else{ // interior 
            if( ip_delta == 3 || ip_delta == 4)
              continue; //dont do two away
          }

          double psi_curr = phi_lo + ip*dpsi + modifiers[ip_delta]*delta_psi;
          double zmin = inp->zmin, zmax = inp->zmax;

          double arcL, darcL, arcL_curr, arcL_lo;
          double arcL_l, arcL_r;
          double phi_r, phi_l;
          if(inp->ftype == GKYL_CORE){
            //Find the turning points
            double zlo, zup, zlo_last;
            zlo = 0.01;
            zup=zmax;
            zlo_last = zlo;
            double R[4], dR[4];
            while(true){
              int nlo = R_psiZ(geo, psi_curr, zlo, 4, R, dR);
              if(nlo==2){
                if(fabs(zlo-zup)<1e-12){
                  zmax = zlo;
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
            zlo=zmin;
            zup_last = zup;
            while(true){
              int nup = R_psiZ(geo, psi_curr, zup, 4, R, dR);
              if(nup==2){
                if(fabs(zlo-zup)<1e-12){
                  zmin = zup;
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
            arcL_r = integrate_psi_contour_memo(geo, psi_curr, zmin, zmax, rright,
              true, true, arc_memo_right);
            arc_ctx.arcL_right = arcL_r;
            arc_ctx.right = false;
            arcL_l = integrate_psi_contour_memo(geo, psi_curr, zmin, zmax, rleft,
              true, true, arc_memo_left);
            arcL = arcL_l + arcL_r;
            darcL = arcL/(up->basis.poly_order*inp->cgrid.cells[TH_IDX]) * (inp->cgrid.upper[TH_IDX] - inp->cgrid.lower[TH_IDX])/2/M_PI;

            arc_ctx.right = true;
            arc_ctx.phi_right = 0.0;
            arc_ctx.rclose = rright;
            arc_ctx.psi = psi_curr;
            phi_r = phi_func(alpha_curr, zmax, &arc_ctx);
            arc_ctx.phi_right = phi_r - alpha_curr; // otherwise alpha will get added on twice
          }
          else if(inp->ftype == GKYL_PF_LO){
            //Find the  upper turning point
            double zlo, zup, zlo_last;
            zlo = fmax(inp->zmin_left, inp->zmin_right);
            zup=zmax;
            zlo_last = zlo;
            double R[4], dR[4];
            while(true){
              int nlo = R_psiZ(geo, psi_curr, zlo, 4, R, dR);
              if(nlo>=2){
                if(fabs(zlo-zup)<1e-12){
                  zmax = zlo;
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
            arcL_r = integrate_psi_contour_memo(geo, psi_curr, inp->zmin_right, zmax, rright,
              true, true, arc_memo);
            arc_ctx.arcL_right = arcL_r;
            arc_ctx.right = false;
            arcL_l = integrate_psi_contour_memo(geo, psi_curr, inp->zmin_left, zmax, rleft,
              true, true, arc_memo);
            arcL = arcL_l + arcL_r;
            darcL = arcL/(up->basis.poly_order*inp->cgrid.cells[TH_IDX]) * (inp->cgrid.upper[TH_IDX] - inp->cgrid.lower[TH_IDX])/2/M_PI;

            arc_ctx.right = true;
            arc_ctx.phi_right = 0.0;
            arc_ctx.rclose = rright;
            arc_ctx.psi = psi_curr;
            arc_ctx.zmin = inp->zmin_right;
            phi_r = phi_func(alpha_curr, zmax, &arc_ctx);
            arc_ctx.phi_right = phi_r - alpha_curr; // otherwise alpha will get added on twice
          }
          else if(inp->ftype == GKYL_PF_UP){
            //Find the lower turning point
            double zlo, zup, zlo_last;
            double zup_last;
            zup = fmin(inp->zmax_left, inp->zmax_right);
            zlo=zmin;
            zup_last = zup;
            double R[4], dR[4];
            while(true){
              int nup = R_psiZ(geo, psi_curr, zup, 4, R, dR);
              if(nup>=2){
                if(fabs(zlo-zup)<1e-12){
                  zmin = zup;
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
            arcL_r = integrate_psi_contour_memo(geo, psi_curr, zmin, inp->zmax_right, rright,
              true, true, arc_memo);
            arc_ctx.right = false;
            arcL_l = integrate_psi_contour_memo(geo, psi_curr, zmin, inp->zmax_left, rleft,
              true, true, arc_memo);
            arc_ctx.arcL_left= arcL_l;
            arcL = arcL_l + arcL_r;
            darcL = arcL/(up->basis.poly_order*inp->cgrid.cells[TH_IDX]) * (inp->cgrid.upper[TH_IDX] - inp->cgrid.lower[TH_IDX])/2/M_PI;

            arc_ctx.right = false;
            arc_ctx.phi_left = 0.0;
            arc_ctx.rclose = rleft;
            arc_ctx.psi = psi_curr;
            arc_ctx.zmax = inp->zmax_left;
            phi_l = phi_func(alpha_curr, zmin, &arc_ctx);
            arc_ctx.phi_left = phi_l - alpha_curr; // otherwise alpha will get added on twice
          }
          else if(inp->ftype==GKYL_SOL_DN_OUT){
            if (geo->plate_spec){ // if we dont have a fixed zmin and zmax
              double rzplate[2];
              pctx.psi_curr = psi_curr;
              pctx.lower=false;
              double a = 0;
              double b = 1;
              double fa = plate_psi_func(a, &pctx);
              double fb = plate_psi_func(b, &pctx);
              struct gkyl_qr_res res = gkyl_ridders(plate_psi_func, &pctx,
                a, b, fa, fb, geo->root_param.max_iter, 1e-10);
              double smax = res.res;
              geo->plate_func_upper(smax, rzplate);
              zmax = rzplate[1];

              pctx.lower=true;
              a = 0;
              b = 1;
              fa = plate_psi_func(a, &pctx);
              fb = plate_psi_func(b, &pctx);
              res = gkyl_ridders(plate_psi_func, &pctx,
                a, b, fa, fb, geo->root_param.max_iter, 1e-10);
              double smin = res.res;
              geo->plate_func_lower(smin, rzplate);
              zmin = rzplate[1];
            }

            arc_ctx.phi_right = 0.0;
            arcL = integrate_psi_contour_memo(geo, psi_curr, zmin, zmax, rclose, true, true, arc_memo);
            darcL = arcL/(up->basis.poly_order*inp->cgrid.cells[TH_IDX]) * (inp->cgrid.upper[TH_IDX] - inp->cgrid.lower[TH_IDX])/2/M_PI;
          }
          else if(inp->ftype==GKYL_SOL_DN_IN){
            arc_ctx.phi_right = 0.0;
            arcL = integrate_psi_contour_memo(geo, psi_curr, zmin, zmax, rclose, true, true, arc_memo);
            darcL = arcL/(up->basis.poly_order*inp->cgrid.cells[TH_IDX]) * (inp->cgrid.upper[TH_IDX] - inp->cgrid.lower[TH_IDX])/2/M_PI;
          }
          else if(inp->ftype == GKYL_SOL_SN_LO){
            //Find the  upper turning point
            double zlo, zup, zlo_last;
            zlo = fmax(inp->zmin_left, inp->zmin_right);
            zup=zmax;
            zlo_last = zlo;
            double R[4], dR[4];
            while(true){
              int nlo = R_psiZ(geo, psi_curr, zlo, 4, R, dR);
              if(nlo>=2){
                if(fabs(zlo-zup)<1e-12){
                  zmax = zlo;
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
            arcL_r = integrate_psi_contour_memo(geo, psi_curr, inp->zmin_right, zmax, rright,
              true, true, arc_memo);
            arc_ctx.arcL_right = arcL_r;
            arc_ctx.right = false;
            arcL_l = integrate_psi_contour_memo(geo, psi_curr, inp->zmin_left, zmax, rleft,
              true, true, arc_memo);
            arcL = arcL_l + arcL_r;
            darcL = arcL/(up->basis.poly_order*inp->cgrid.cells[TH_IDX]) * (inp->cgrid.upper[TH_IDX] - inp->cgrid.lower[TH_IDX])/2/M_PI;

            arc_ctx.right = true;
            arc_ctx.phi_right = 0.0;
            arc_ctx.rclose = rright;
            arc_ctx.psi = psi_curr;
            arc_ctx.zmin = inp->zmin_right;
            phi_r = phi_func(alpha_curr, zmax, &arc_ctx);
            arc_ctx.phi_right = phi_r - alpha_curr; // otherwise alpha will get added on twice
          }

          // at the beginning of each theta loop we need to reset things
          cidx[PH_IDX] = ip;
          arcL_curr = 0.0;
          arcL_lo = (theta_lo + M_PI)/2/M_PI*arcL;
          double ridders_min, ridders_max;
          // set node coordinates
          for (int it=nrange->lower[TH_IDX]; it<=nrange->upper[TH_IDX]; ++it) {
            int it_delta_max = 5; // should be 5
            if(ia_delta != 0 || ip_delta != 0 )
              it_delta_max = 1;
            for(int it_delta = 0; it_delta < it_delta_max; it_delta++){
              if(it == nrange->lower[TH_IDX]){
                if(it_delta == 1 || it_delta == 3)
                  continue; // want to use one sided stencils at edge
              }
              else if(it == nrange->upper[TH_IDX]){
                if(it_delta == 2 || it_delta == 4)
                  continue; // want to use one sided stencils at edge
              }
              else{
                if( it_delta == 3 || it_delta == 4)
                  continue; //dont do two away
              }
              arcL_curr = arcL_lo + it*darcL + modifiers[it_delta]*delta_theta*(arcL/2/M_PI);
              double theta_curr = arcL_curr*(2*M_PI/arcL) - M_PI ; 
              if(inp->ftype==GKYL_CORE){
                if(arcL_curr <= arcL_r){
                  rclose = rright;
                  arc_ctx.right = true;
                  ridders_min = -arcL_curr;
                  ridders_max = arcL-arcL_curr;
                  arc_ctx.zmin = zmin;
                  arc_ctx.zmax = zmax;
                }
                else{
                  rclose = rleft;
                  arc_ctx.right = false;
                  ridders_min = arcL - arcL_curr;
                  ridders_max = -arcL_curr + arc_ctx.arcL_right;
                  arc_ctx.zmin = zmin;
                  arc_ctx.zmax = zmax;
                }
              }
              if(inp->ftype==GKYL_PF_LO){
                if(arcL_curr <= arcL_r){
                  rclose = rright;
                  arc_ctx.right = true;
                  ridders_min = -arcL_curr;
                  ridders_max = arcL-arcL_curr;
                  arc_ctx.zmin = inp->zmin_right;
                  arc_ctx.zmax = zmax;
                }
                else{
                  rclose = rleft;
                  arc_ctx.right = false;
                  ridders_min = arcL - arcL_curr;
                  ridders_max = -arcL_curr + arc_ctx.arcL_right;
                  arc_ctx.zmin = inp->zmin_left;
                  arc_ctx.zmax = zmax;
                }
              }
              if(inp->ftype==GKYL_PF_UP){
                if(arcL_curr > arcL_l){
                  rclose = rright;
                  arc_ctx.right = true;
                  ridders_min = arc_ctx.arcL_left - arcL_curr;
                  ridders_max = arcL - arcL_curr;
                  arc_ctx.zmin = zmin;
                  arc_ctx.zmax = inp->zmax_right;
                }
                else{
                  rclose = rleft;
                  arc_ctx.right = false;
                  ridders_min = arc_ctx.arcL_left - arcL_curr;
                  ridders_max = -arcL_curr;
                  arc_ctx.zmin = zmin;
                  arc_ctx.zmax = inp->zmax_left;
                }
              }
              if(arc_ctx.ftype==GKYL_SOL_DN_OUT){
                ridders_min = -arcL_curr;
                ridders_max = arcL-arcL_curr;
                arc_ctx.right = false;
                arc_ctx.zmin = zmin;
                arc_ctx.zmax = zmax;
              }
              if(arc_ctx.ftype==GKYL_SOL_DN_IN){
                ridders_min = arcL-arcL_curr;
                ridders_max = -arcL_curr;
                arc_ctx.right = false;
                arc_ctx.zmin = zmin;
                arc_ctx.zmax = zmax;
              }
              if(arc_ctx.ftype==GKYL_SOL_SN_LO){
                if(arcL_curr <= arcL_r){
                  rclose = rright;
                  arc_ctx.right = true;
                  ridders_min = -arcL_curr;
                  ridders_max = arcL-arcL_curr;
                  arc_ctx.zmin = inp->zmin_right;
                  arc_ctx.zmax = zmax;
                }
                else{
                  rclose = rleft;
                  arc_ctx.right = false;
                  ridders_min = arcL - arcL_curr;
                  ridders_max = -arcL_curr + arc_ctx.arcL_right;
                  arc_ctx.zmin = inp->zmin_left;
                  arc_ctx.zmax = zmax;
                }
              }

              arc_ctx.psi = psi_curr;
              arc_ctx.rclose = rclose;
              arc_ctx.arcL = arcL_curr;
              struct gkyl_qr_res res = gkyl_ridders(arc_length_func, &arc_ctx,
                arc_ctx.zmin, arc_ctx.zmax, ridders_min, ridders_max,
                geo->root_param.max_iter, 1e-10);
              double z_curr = res.res;
              ((gkyl_tok_geo *)geo)->stat.nroot_cont_calls += res.nevals;
              double R[4] = { 0 }, dR[4] = { 0 };
              int nr = R_psiZ(geo, psi_curr, z_curr, 4, R, dR);
              double r_curr = choose_closest(rclose, R, R, nr);

              cidx[TH_IDX] = it;
              int lidx = 0;
              if (ip_delta != 0)
                lidx = 3 + 3*(ip_delta-1);
              if (ia_delta != 0)
                lidx = 15 + 3*(ia_delta-1);
              if (it_delta != 0)
                lidx = 27 + 3*(it_delta-1);

              if(ia_delta==0 && ip_delta==0 && it_delta==0)
                lidx = 0;

              double phi_curr = phi_func(alpha_curr, z_curr, &arc_ctx);
              double *mc2p_fd_n = gkyl_array_fetch(mc2p_nodal_fd, gkyl_range_idx(nrange, cidx));
              double *mc2p_n = gkyl_array_fetch(mc2p_nodal, gkyl_range_idx(nrange, cidx));
              mc2p_fd_n[lidx+X_IDX] = r_curr*cos(phi_curr);
              mc2p_fd_n[lidx+Y_IDX] = r_curr*sin(phi_curr);
              mc2p_fd_n[lidx+Z_IDX] = z_curr;

              if(ip_delta==0 && ia_delta==0 && it_delta==0){
                mc2p_n[X_IDX] = r_curr*cos(phi_curr);
                mc2p_n[Y_IDX] = r_curr*sin(phi_curr);
                mc2p_n[Z_IDX] = z_curr;
              }
            }
          }
        }
      }
    }
  }
  gkyl_nodal_ops_n2m(&inp->cbasis, &inp->cgrid, nrange, &up->range, 3, mc2p_nodal, mc2p);

  char str1[50] = "xyz";
  char str2[50] = "allxyz";
  if (inp->write_node_coord_array){
    write_nodal_coordinates(strcat(str1, inp->node_file_nm), nrange, mc2p_nodal);
    write_nodal_coordinates(strcat(str2, inp->node_file_nm), nrange, mc2p_nodal_fd);
  }
  gkyl_free(arc_memo);
  gkyl_free(arc_memo_left);
  gkyl_free(arc_memo_right);
}


struct gkyl_tok_geo_stat
gkyl_tok_geo_get_stat(const gkyl_tok_geo *geo)
{
  return geo->stat;
}

void
gkyl_tok_geo_release(gkyl_tok_geo *geo)
{
  gkyl_efit_release(geo->efit);
  gkyl_free(geo);
}
