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

#include <math.h>
#include <string.h>


// some helper functions
double
choose_closest(double ref, double* R, double* out, int nr)
{
  //return fabs(R[0]-ref) < fabs(R[1]-ref) ? out[0] : out[1];
  int imin = 0;
  double min = fabs(R[0]-ref);
  for(int i = 1; i< nr; i++){
    if( fabs(R[i] - ref) < min){
      imin = i;
      min = fabs(R[i] - ref);
    }
  }
  return out[imin];
}

static inline double SQ(double x) { return x*x; }

static inline int
get_idx(int dir, double x, const struct gkyl_rect_grid *grid, const struct gkyl_range *range)
{
  double xlower = grid->lower[dir], dx = grid->dx[dir];
  int idx = range->lower[dir] + (int) floor((x-xlower)/dx);
  return idx <= range->upper[dir] ? idx : range->upper[dir];
}

// struct for solutions to roots
struct RdRdZ_sol {
  int nsol;
  double R[2], dRdZ[2];
};

// Compute roots R(psi,Z) and dR/dZ(psi,Z) in a p=1 DG cell
static inline struct RdRdZ_sol
calc_RdR_p1(const double *psi, double psi0, double Z, double xc[2], double dx[2])
{
  struct RdRdZ_sol sol = { .nsol = 0 };

  double y = (Z-xc[1])/(dx[1]*0.5);
  
  double rnorm = (-(1.732050807568877*psi[2]*y)/(3.0*psi[3]*y+1.732050807568877*psi[1]))+(2.0*psi0)/(3.0*psi[3]*y+1.732050807568877*psi[1])-(1.0*psi[0])/(3.0*psi[3]*y+1.732050807568877*psi[1]) ;

  if ((-1<=rnorm) && (rnorm < 1)) {
    double drdznorm = -(3.0*(2.0*psi[3]*psi0-1.0*psi[0]*psi[3]+psi[1]*psi[2]))/SQ(3.0*psi[3]*y+1.732050807568877*psi[1]) ;
    
    sol.nsol = 1;
    sol.R[0] = rnorm*dx[0]*0.5 + xc[0];
    sol.dRdZ[0] = drdznorm*dx[0]/dx[1];
  }
  return sol;
}

// Compute roots R(psi,Z) and dR/dZ(psi,Z) in a p=2 DG cell
static inline struct RdRdZ_sol
calc_RdR_p2(const double *psi, double psi0, double Z, double xc[2], double dx[2])
{
  struct RdRdZ_sol sol = { .nsol = 0 };
  double y = (Z-xc[1])/(dx[1]*0.5);

  double aq = 2.904737509655563*psi[6]*y+1.677050983124842*psi[4]; 
  double bq = 2.904737509655563*psi[7]*SQ(y)+1.5*psi[3]*y-0.9682458365518543*psi[7]+0.8660254037844386*psi[1]; 
  double cq = 1.677050983124842*psi[5]*SQ(y)-0.9682458365518543*psi[6]*y+0.8660254037844386*psi[2]*y-1.0*psi0-0.5590169943749475*psi[5]-0.5590169943749475*psi[4]+0.5*psi[0]; 
  double delta2 = bq*bq - 4*aq*cq;

  if (delta2 > 0) {
    double r1, r2;
    double delta = sqrt(delta2);
    // compute both roots
    if (bq>=0) {
      r1 = (-bq-delta)/(2*aq);
      r2 = 2*cq/(-bq-delta);
    }
    else {
      r1 = 2*cq/(-bq+delta);
      r2 = (-bq+delta)/(2*aq);
    }

    int sidx = 0;
    if ((-1<=r1) && (r1 < 1)) {
      sol.nsol += 1;
      sol.R[sidx] = r1*dx[0]*0.5 + xc[0];

      double x = r1;
      double C = 5.809475019311126*psi[7]*x*y+3.354101966249685*psi[5]*y+2.904737509655563*psi[6]*SQ(x)+1.5*psi[3]*x-0.9682458365518543*psi[6]+0.8660254037844386*psi[2]; 
      double A = 2.904737509655563*psi[7]*SQ(y)+5.809475019311126*psi[6]*x*y+1.5*psi[3]*y+3.354101966249685*psi[4]*x-0.9682458365518543*psi[7]+0.8660254037844386*psi[1];
      sol.dRdZ[sidx] = -C/A*dx[0]/dx[1];
      
      sidx += 1;
    }
    if ((-1<=r2) && (r2 < 1)) {
      sol.nsol += 1;
      sol.R[sidx] = r2*dx[0]*0.5 + xc[0];

      double x = r2;
      double C = 5.809475019311126*psi[7]*x*y+3.354101966249685*psi[5]*y+2.904737509655563*psi[6]*SQ(x)+1.5*psi[3]*x-0.9682458365518543*psi[6]+0.8660254037844386*psi[2]; 
      double A = 2.904737509655563*psi[7]*SQ(y)+5.809475019311126*psi[6]*x*y+1.5*psi[3]*y+3.354101966249685*psi[4]*x-0.9682458365518543*psi[7]+0.8660254037844386*psi[1];
      sol.dRdZ[sidx] = -C/A*dx[0]/dx[1];
      
      sidx += 1;
    }
  }
  return sol;
}

// Compute R(psi,Z) given a psi and Z. Can return multiple solutions
// or no solutions. The number of roots found is returned and are
// copied in the array R and dR. The calling function must ensure that
// these arrays are big enough to hold all roots required
static int
R_psiZ(const gkyl_tok_geo *geo, double psi, double Z, int nmaxroots,
  double *R, double *dR)
{
  int zcell = get_idx(1, Z, geo->rzgrid, geo->rzlocal);

  int sidx = 0;
  int idx[2] = { 0, zcell };
  double dx[2] = { geo->rzgrid->dx[0], geo->rzgrid->dx[1] };
  
  struct gkyl_range rangeR;
  gkyl_range_deflate(&rangeR, geo->rzlocal, (int[]) { 0, 1 }, (int[]) { 0, zcell });

  struct gkyl_range_iter riter;
  gkyl_range_iter_init(&riter, &rangeR);
  
  // loop over all R cells to find psi crossing
  while (gkyl_range_iter_next(&riter) && sidx<=nmaxroots) {
    long loc = gkyl_range_idx(&rangeR, riter.idx);
    const double *psih = gkyl_array_cfetch(geo->psiRZ, loc);

    double xc[2];
    idx[0] = riter.idx[0];
    gkyl_rect_grid_cell_center(geo->rzgrid, idx, xc);

    struct RdRdZ_sol sol = geo->calc_roots(psih, psi, Z, xc, dx);
    
    if (sol.nsol > 0)
      for (int s=0; s<sol.nsol; ++s) {
        R[sidx] = sol.R[s];
        dR[sidx] = sol.dRdZ[s];
        sidx += 1;
      }
  }
  return sidx;
}


// Function context to pass to coutour integration function
struct contour_ctx {
  const gkyl_tok_geo *geo;
  double psi, last_R;
  long ncall;
};

// Function to pass to numerical quadrature to integrate along a contour
static inline double
contour_func(double Z, void *ctx)
{
  struct contour_ctx *c = ctx;
  c->ncall += 1;
  double R[4] = { 0 }, dR[4] = { 0 };
  
  int nr = R_psiZ(c->geo, c->psi, Z, 4, R, dR);
  double dRdZ = nr == 1 ? dR[0] : choose_closest(c->last_R, R, dR,nr);
  
  return nr>0 ? sqrt(1+dRdZ*dRdZ) : 0.0;
}

static inline double
phi_contour_func(double Z, void *ctx)
{
  struct contour_ctx *c = ctx;
  c->ncall += 1;
  double R[4] = { 0 }, dR[4] = { 0 };
  
  int nr = R_psiZ(c->geo, c->psi, Z, 4, R, dR);
  double dRdZ = nr == 1 ? dR[0] : choose_closest(c->last_R, R, dR, nr);
  double r_curr = nr == 1 ? R[0] : choose_closest(c->last_R, R, R, nr);

  struct gkyl_range_iter iter;
  iter.idx[0] = fmin(c->geo->rzlocal->lower[0] + (int) floor((r_curr - c->geo->rzgrid->lower[0])/c->geo->rzgrid->dx[0]), c->geo->rzlocal->upper[0]);
  iter.idx[1] = fmin(c->geo->rzlocal->lower[1] + (int) floor((Z - c->geo->rzgrid->lower[1])/c->geo->rzgrid->dx[1]), c->geo->rzlocal->upper[1]);
  long loc = gkyl_range_idx((c->geo->rzlocal), iter.idx);
  const double *psih = gkyl_array_cfetch(c->geo->psiRZ, loc);

  double xc[2];
  gkyl_rect_grid_cell_center((c->geo->rzgrid), iter.idx, xc);
  double x = (r_curr-xc[0])/(c->geo->rzgrid->dx[0]*0.5);
  double y = (Z-xc[1])/(c->geo->rzgrid->dx[1]*0.5);

  // if psi is polyorder 2 we can get grad psi
  // in cylindrical coords it is grad psi = dpsi/dR Rhat + dpsi/dZ zhat
  double dpsidx = 2.904737509655563*psih[7]*(y*y-0.3333333333333333)+5.809475019311126*psih[6]*x*y+1.5*psih[3]*y+3.354101966249684*psih[4]*x+0.8660254037844386*psih[1]; 
  double dpsidy =	5.809475019311126*psih[7]*x*y+3.354101966249684*psih[5]*y+2.904737509655563*psih[6]*(x*x-0.3333333333333333)+1.5*psih[3]*x+0.8660254037844386*psih[2];
  dpsidx = dpsidx*2.0/c->geo->rzgrid->dx[0];
  dpsidy = dpsidy*2.0/c->geo->rzgrid->dx[1];
  double grad_psi_mag = sqrt(dpsidx*dpsidx + dpsidy*dpsidy);
  double result  = (1/r_curr/sqrt(dpsidx*dpsidx + dpsidy*dpsidy)) *sqrt(1+dRdZ*dRdZ) ;
  return nr>0 ? result : 0.0;
}

// Integrates along a specified contour, optionally using a "memory"
// of previously computed values, or storing computed values in
// memory. The function basically breaks up the integral into a loop
// over z-cells. This needs to be done as the DG representation is,
// well, discontinuous, and adaptive quadrature struggles with such
// functions.
static double
integrate_psi_contour_memo(const gkyl_tok_geo *geo, double psi,
  double zmin, double zmax, double rclose,
  bool use_memo, bool fill_memo, double *memo)
{
  struct contour_ctx ctx = {
    .geo = geo,
    .psi = psi,
    .ncall = 0,
    .last_R = rclose
  };

  int nlevels = geo->quad_param.max_level;
  double eps = geo->quad_param.eps;
  
  double dz = geo->rzgrid->dx[1];
  double zlo = geo->rzgrid->lower[1];
  int izlo = geo->rzlocal->lower[1], izup = geo->rzlocal->upper[1];
  
  int ilo = get_idx(1, zmin, geo->rzgrid, geo->rzlocal);
  int iup = get_idx(1, zmax, geo->rzgrid, geo->rzlocal);

  double res = 0.0;
  for (int i=ilo; i<=iup; ++i) {
    double z1 = gkyl_median(zmin, zlo+(i-izlo)*dz, zlo+(i-izlo+1)*dz);
    double z2 = gkyl_median(zmax, zlo+(i-izlo)*dz, zlo+(i-izlo+1)*dz);
    
    if (z1 < z2) {
      if (use_memo) {
        if (fill_memo) {
          struct gkyl_qr_res res_local =
            gkyl_dbl_exp(contour_func, &ctx, z1, z2, nlevels, eps);
          memo[i-izlo] = res_local.res;
          res += res_local.res;
        }
        else {
          if (z2-z1 == dz) {
            res += memo[i-izlo];
          }
          else {
            struct gkyl_qr_res res_local =
              gkyl_dbl_exp(contour_func, &ctx, z1, z2, nlevels, eps);
            res += res_local.res;
          }
        }
      }
      else {
        struct gkyl_qr_res res_local =
          gkyl_dbl_exp(contour_func, &ctx, z1, z2, nlevels, eps);
        res += res_local.res;
      }
    }
  }

  ((gkyl_tok_geo *)geo)->stat.nquad_cont_calls += ctx.ncall;
  return res;
}

static double
integrate_phi_along_psi_contour_memo(const gkyl_tok_geo *geo, double psi,
  double zmin, double zmax, double rclose,
  bool use_memo, bool fill_memo, double *memo)
{
  struct contour_ctx ctx = {
    .geo = geo,
    .psi = psi,
    .ncall = 0,
    .last_R = rclose
  };

  int nlevels = geo->quad_param.max_level;
  double eps = geo->quad_param.eps;
  
  double dz = geo->rzgrid->dx[1];
  double zlo = geo->rzgrid->lower[1];
  int izlo = geo->rzlocal->lower[1], izup = geo->rzlocal->upper[1];
  
  int ilo = get_idx(1, zmin, geo->rzgrid, geo->rzlocal);
  int iup = get_idx(1, zmax, geo->rzgrid, geo->rzlocal);

  double res = 0.0;
  for (int i=ilo; i<=iup; ++i) {
    double z1 = gkyl_median(zmin, zlo+(i-izlo)*dz, zlo+(i-izlo+1)*dz);
    double z2 = gkyl_median(zmax, zlo+(i-izlo)*dz, zlo+(i-izlo+1)*dz);
    
    if (z1 < z2) {
      if (use_memo) {
        if (fill_memo) {
          struct gkyl_qr_res res_local =
            gkyl_dbl_exp(phi_contour_func, &ctx, z1, z2, nlevels, eps);
          memo[i-izlo] = res_local.res;
          res += res_local.res;
        }
        else {
          if (z2-z1 == dz) {
            res += memo[i-izlo];
          }
          else {
            struct gkyl_qr_res res_local =
              gkyl_dbl_exp(phi_contour_func, &ctx, z1, z2, nlevels, eps);
            res += res_local.res;
          }
        }
      }
      else {
        struct gkyl_qr_res res_local =
          gkyl_dbl_exp(phi_contour_func, &ctx, z1, z2, nlevels, eps);
        res += res_local.res;
      }
    }
  }

  ((gkyl_tok_geo *)geo)->stat.nquad_cont_calls += ctx.ncall;
  return res;
}

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
plate_psi_func(double R, void *ctx){
  // uses a pointer to the plate function to get Z(R)
  // Then calculates psi(R, Z(R))
  // will be used by ridders later
  
  struct plate_ctx *gc = ctx;
  // First find R(z)
  double Z;
  if(gc->lower==true)
    Z = gc->geo->plate_func_lower(R);
  else
    Z = gc->geo->plate_func_upper(R);

  // Now find the cell where this R and Z is
  int rzidx[2];
  rzidx[0] = fmin(gc->geo->rzlocal->lower[0] + (int) floor((R - gc->geo->rzgrid->lower[0])/gc->geo->rzgrid->dx[0]), gc->geo->rzlocal->upper[0]);
  rzidx[1] = fmin(gc->geo->rzlocal->lower[1] + (int) floor((Z - gc->geo->rzgrid->lower[1])/gc->geo->rzgrid->dx[1]), gc->geo->rzlocal->upper[1]);
  long loc = gkyl_range_idx(gc->geo->rzlocal, rzidx);
  const double *coeffs = gkyl_array_cfetch(gc->geo->psiRZ,loc);

  double xc[2];
  gkyl_rect_grid_cell_center(gc->geo->rzgrid, rzidx, xc);
  double xy[2];
  xy[0] = (R-xc[0])/(gc->geo->rzgrid->dx[0]*0.5);
  xy[1] = (Z-xc[1])/(gc->geo->rzgrid->dx[1]*0.5);
  double psi = gc->geo->rzbasis->eval_expand(xy, coeffs);
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
  int idx = fmin(actx->geo->frange->lower[0] + (int) floor((psi_fpol - actx->geo->fgrid->lower[0])/actx->geo->fgrid->dx[0]), actx->geo->frange->upper[0]);
  long loc = gkyl_range_idx(actx->geo->frange, &idx);
  const double *coeffs = gkyl_array_cfetch(actx->geo->fpoldg,loc);
  double fxc;
  gkyl_rect_grid_cell_center(actx->geo->fgrid, &idx, &fxc);
  double fx = (psi_fpol-fxc)/(actx->geo->fgrid->dx[0]*0.5);
  double fpol = actx->geo->fbasis->eval_expand(&fx, coeffs);
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

  struct gkyl_efit* efit = gkyl_efit_new(inp->filepath, inp->rzpoly_order, inp->fluxpoly_order, false);

  geo->plate_spec = inp->plate_spec;
  geo->plate_func_lower = inp->plate_func_lower;
  geo->plate_func_upper = inp->plate_func_upper;
  geo->plate_lower_Rl = inp->plate_lower_Rl;
  geo->plate_lower_Rr = inp->plate_lower_Rr;
  geo->plate_upper_Rl = inp->plate_upper_Rl;
  geo->plate_upper_Rr = inp->plate_upper_Rr;

  geo->rzbasis= efit->rzbasis;
  geo->rzgrid = efit->rzgrid;
  geo->psiRZ = gkyl_array_acquire(efit->psizr);
  geo->psibyrRZ = gkyl_array_acquire(efit->psibyrzr);
  geo->psibyr2RZ = gkyl_array_acquire(efit->psibyr2zr);

  geo->num_rzbasis = efit->rzbasis->num_basis;
  geo->rzlocal = efit->rzlocal;
  geo->rzlocal_ext = efit->rzlocal_ext;
  geo->fgrid = efit->fluxgrid;
  geo->fbasis = efit->fluxbasis;
  geo->frange = efit->fluxlocal;
  geo->frange_ext = efit->fluxlocal_ext;
  geo->fpoldg= gkyl_array_acquire(efit->fpolflux);
  geo->qdg= gkyl_array_acquire(efit->qflux);
  geo->psisep = efit->sibry;
  geo->zmaxis = efit->zmaxis;

  geo->root_param.eps =
    inp->root_param.eps > 0 ? inp->root_param.eps : 1e-10;
  geo->root_param.max_iter =
    inp->root_param.max_iter > 0 ? inp->root_param.max_iter : 100;

  geo->quad_param.max_level =
    inp->quad_param.max_levels > 0 ? inp->quad_param.max_levels : 10;
  geo->quad_param.eps =
    inp->quad_param.eps > 0 ? inp->quad_param.eps : 1e-10;

  if (efit->rzbasis->poly_order == 1)
    geo->calc_roots = calc_RdR_p1;
  else if (efit->rzbasis->poly_order == 2)
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


void gkyl_tok_geo_advance(struct gk_geometry* up, struct gkyl_range *nrange, double dzc[3], 
  evalf_t mapc2p_func, void* mapc2p_ctx, evalf_t bmag_func, void *bmag_ctx, 
  struct gkyl_array *mc2p_nodal_fd, struct gkyl_array *mc2p_nodal, struct gkyl_array *mc2p, 
  struct gkyl_array *bmag, struct gkyl_array *g_ij, 
  struct gkyl_array *jacobgeo, struct gkyl_array *jacobgeo_inv, 
  struct gkyl_array *gij, struct gkyl_array *b_i, struct gkyl_array *cmag, struct gkyl_array *jacobtot, 
  struct gkyl_array *jacobtot_inv, struct gkyl_array *bmag_inv, struct gkyl_array *bmag_inv_sq, 
  struct gkyl_array *gxxj, struct gkyl_array *gxyj, struct gkyl_array *gyyj) 
{

  struct gkyl_tok_geo *geo = mapc2p_ctx;
  struct gkyl_tok_geo_geo_inp *inp = bmag_ctx;

  enum { PH_IDX, AL_IDX, TH_IDX }; // arrangement of computational coordinates
  enum { X_IDX, Y_IDX, Z_IDX }; // arrangement of cartesian coordinates
  
  double dtheta = inp->cgrid->dx[TH_IDX],
    dpsi = inp->cgrid->dx[PH_IDX],
    dalpha = inp->cgrid->dx[AL_IDX];
  
  double theta_lo = inp->cgrid->lower[TH_IDX],
    phi_lo = inp->cgrid->lower[PH_IDX],
    alpha_lo = inp->cgrid->lower[AL_IDX];

  double dx_fact = up->basis->poly_order == 1 ? 1 : 0.5;
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

  int nzcells = geo->rzgrid->cells[1];
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
            darcL = arcL/(up->basis->poly_order*inp->cgrid->cells[TH_IDX]) * (inp->cgrid->upper[TH_IDX] - inp->cgrid->lower[TH_IDX])/2/M_PI;

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
            darcL = arcL/(up->basis->poly_order*inp->cgrid->cells[TH_IDX]) * (inp->cgrid->upper[TH_IDX] - inp->cgrid->lower[TH_IDX])/2/M_PI;

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
            darcL = arcL/(up->basis->poly_order*inp->cgrid->cells[TH_IDX]) * (inp->cgrid->upper[TH_IDX] - inp->cgrid->lower[TH_IDX])/2/M_PI;

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
              pctx.psi_curr = psi_curr;
              pctx.lower=false;
              double a = geo->plate_upper_Rl;
              double b = geo->plate_upper_Rr;
              double fa = plate_psi_func(a, &pctx);
              double fb = plate_psi_func(b, &pctx);
              struct gkyl_qr_res res = gkyl_ridders(plate_psi_func, &pctx,
                a, b, fa, fb, geo->root_param.max_iter, 1e-10);
              double rmax = res.res;
              zmax = geo->plate_func_upper(rmax);

              pctx.lower=true;
              a = geo->plate_lower_Rl;
              b = geo->plate_lower_Rr;
              fa = plate_psi_func(a, &pctx);
              fb = plate_psi_func(b, &pctx);
              res = gkyl_ridders(plate_psi_func, &pctx,
                a, b, fa, fb, geo->root_param.max_iter, 1e-10);
              double rmin = res.res;
              zmin = geo->plate_func_lower(rmin);
            }

            arc_ctx.phi_right = 0.0;
            arcL = integrate_psi_contour_memo(geo, psi_curr, zmin, zmax, rclose, true, true, arc_memo);
            darcL = arcL/(up->basis->poly_order*inp->cgrid->cells[TH_IDX]) * (inp->cgrid->upper[TH_IDX] - inp->cgrid->lower[TH_IDX])/2/M_PI;
          }
          else if(inp->ftype==GKYL_SOL_DN_IN){
            arc_ctx.phi_right = 0.0;
            arcL = integrate_psi_contour_memo(geo, psi_curr, zmin, zmax, rclose, true, true, arc_memo);
            darcL = arcL/(up->basis->poly_order*inp->cgrid->cells[TH_IDX]) * (inp->cgrid->upper[TH_IDX] - inp->cgrid->lower[TH_IDX])/2/M_PI;
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
            darcL = arcL/(up->basis->poly_order*inp->cgrid->cells[TH_IDX]) * (inp->cgrid->upper[TH_IDX] - inp->cgrid->lower[TH_IDX])/2/M_PI;

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
  gkyl_nodal_ops_n2m(inp->cbasis, inp->cgrid, nrange, up->range, 3, mc2p_nodal, mc2p);

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
  gkyl_free(geo);
}
