#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_rio.h>
#include <gkyl_basis.h>
#include <gkyl_geo_gyrokinetic.h>
#include <gkyl_math.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_nodal_ops.h>

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
R_psiZ(const gkyl_geo_gyrokinetic *geo, double psi, double Z, int nmaxroots,
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
  const gkyl_geo_gyrokinetic *geo;
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
integrate_psi_contour_memo(const gkyl_geo_gyrokinetic *geo, double psi,
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

  ((gkyl_geo_gyrokinetic *)geo)->stat.nquad_cont_calls += ctx.ncall;
  return res;
}

static double
integrate_phi_along_psi_contour_memo(const gkyl_geo_gyrokinetic *geo, double psi,
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

  ((gkyl_geo_gyrokinetic *)geo)->stat.nquad_cont_calls += ctx.ncall;
  return res;
}

// Function context to pass to root finder
struct arc_length_ctx {
  const gkyl_geo_gyrokinetic *geo;
  double *arc_memo;
  double psi, rclose, zmin, arcL;
  double rleft, rright, zmax;
  double arcL_right; // this is for when we need to switch sides
  double arcL_left; // this is for when we need to switch sides
  double phi_right; // this is for when we need to switch sides
  double phi_left; // this is for when we need to switch sides
  bool right;
  double zmaxis;
  enum gkyl_geo_gyrokinetic_type ftype; // type of geometry
};


// Function to pass to root-finder to find Z location for given arc-length
static inline double
arc_length_func(double Z, void *ctx)
{
  struct arc_length_ctx *actx = ctx;
  double *arc_memo = actx->arc_memo;
  double psi = actx->psi, rclose = actx->rclose, zmin = actx->zmin, arcL = actx->arcL;
  double zmax = actx->zmax;
  double ival = 0.0;

  if(actx->ftype==GKYL_SOL_DN_OUT){
    ival = integrate_psi_contour_memo(actx->geo, psi, zmin, Z, rclose, true, true, arc_memo) - arcL;
  }
  if(actx->ftype==GKYL_SOL_DN_IN){
    ival = integrate_psi_contour_memo(actx->geo, psi, Z, zmax, rclose, true, true, arc_memo) - arcL;
  }
  else if(actx->ftype==GKYL_CORE){
    if(actx->right==true){
      ival = integrate_psi_contour_memo(actx->geo, psi, zmin, Z, rclose, false, false, arc_memo) - arcL;
    }
    else{
      ival = integrate_psi_contour_memo(actx->geo, psi, Z, zmax, rclose, false, false, arc_memo)  - arcL + actx->arcL_right;
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

  // Using convention from Noah Mandell's thesis Eq 5.104 phi = alpha at midplane
  double ival = 0;
  double phi_ref = 0.0;
  if(actx->ftype==GKYL_SOL_DN_OUT){
    if(Z<actx->zmaxis)
      ival = -integrate_phi_along_psi_contour_memo(actx->geo, psi, Z, actx->zmaxis, rclose, false, false, arc_memo);
    else
      ival = integrate_phi_along_psi_contour_memo(actx->geo, psi, actx->zmaxis, Z, rclose, false, false, arc_memo);
  }
  if(actx->ftype==GKYL_SOL_DN_IN){
    ival = integrate_phi_along_psi_contour_memo(actx->geo, psi, Z, actx->zmax, rclose, false, false, arc_memo);
  }
  else if(actx->ftype==GKYL_CORE){
    if(actx->right==true){
      if(Z<actx->zmaxis)
        ival = -integrate_phi_along_psi_contour_memo(actx->geo, psi, Z, actx->zmaxis, rclose, false, false, arc_memo);
      else
        ival = integrate_phi_along_psi_contour_memo(actx->geo, psi, actx->zmaxis, Z, rclose, false, false, arc_memo);
    }
    else{
      ival = integrate_phi_along_psi_contour_memo(actx->geo, psi, Z, actx->zmax, rclose, false, false, arc_memo);// + actx->phi_right;
      phi_ref = actx->phi_right;
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
      printf("Z, zmax, ival = %g %g %g\n", Z, actx->zmax, ival);
    }
  }
  // Now multiply by RBphi
  double R[4] = {0};
  double dR[4] = {0};
  int nr = R_psiZ(actx->geo, psi, Z, 4, R, dR);
  double r_curr = nr == 1 ? R[0] : choose_closest(rclose, R, R, nr);
  double Bphi = actx->geo->B0*actx->geo->R0/r_curr;
  ival = ival*r_curr*Bphi;

  // now keep in range 2pi
  while(ival < -M_PI){
    ival +=2*M_PI;
  }
  while(ival > M_PI){
    ival -=2*M_PI;
  }
  return alpha_curr + ival + phi_ref;
}



gkyl_geo_gyrokinetic*
gkyl_geo_gyrokinetic_new(const struct gkyl_geo_gyrokinetic_inp *inp)
{
  struct gkyl_geo_gyrokinetic *geo = gkyl_malloc(sizeof(*geo));

  geo->B0 = inp->B0;
  geo->R0 = inp->R0;

  geo->rzgrid = inp->rzgrid;
  geo->psiRZ = gkyl_array_acquire(inp->psiRZ);
  geo->num_rzbasis = inp->rzbasis->num_basis;
  //memcpy(&geo->rzlocal, inp->rzlocal, sizeof(struct gkyl_range));
  geo->rzlocal = inp->rzlocal;

  geo->root_param.eps =
    inp->root_param.eps > 0 ? inp->root_param.eps : 1e-10;
  geo->root_param.max_iter =
    inp->root_param.max_iter > 0 ? inp->root_param.max_iter : 100;

  geo->quad_param.max_level =
    inp->quad_param.max_levels > 0 ? inp->quad_param.max_levels : 10;
  geo->quad_param.eps =
    inp->quad_param.eps > 0 ? inp->quad_param.eps : 1e-10;

  if (inp->rzbasis->poly_order == 1)
    geo->calc_roots = calc_RdR_p1;
  else if (inp->rzbasis->poly_order == 2)
    geo->calc_roots = calc_RdR_p2;

  geo->stat = (struct gkyl_geo_gyrokinetic_stat) { };
  
  return geo;
}

double
gkyl_geo_gyrokinetic_integrate_psi_contour(const gkyl_geo_gyrokinetic *geo, double psi,
  double zmin, double zmax, double rclose)
{
  return integrate_psi_contour_memo(geo, psi, zmin, zmax, rclose,
    false, false, 0);
}

int
gkyl_geo_gyrokinetic_R_psiZ(const gkyl_geo_gyrokinetic *geo, double psi, double Z, int nmaxroots,
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


void
gkyl_geo_gyrokinetic_calcgeom(gkyl_geo_gyrokinetic *geo,
  const struct gkyl_geo_gyrokinetic_geo_inp *inp, struct gkyl_array *mapc2p, struct gkyl_range *conversion_range)
{
  int poly_order = inp->cbasis->poly_order;
  int nodes[3] = { 1, 1, 1 };
  if (poly_order == 1){
    for (int d=0; d<inp->cgrid->ndim; ++d)
      nodes[d] = inp->cgrid->cells[d] + 1;
    if(inp->bcs[1]==1)
      nodes[1] += 2; // specically alpha gets ghosts if nonperiodic in y
  }
                   
  if (poly_order == 2){
    for (int d=0; d<inp->cgrid->ndim; ++d)
      nodes[d] = 2*(inp->cgrid->cells[d]) + 1;
    nodes[1] += 4; // specically alpha gets ghosts
  }

  for(int d=0; d<inp->cgrid->ndim; d++){
    printf("d[%d] = %d\n", d, nodes[d]);
  }


  geo->nrange = gkyl_malloc(sizeof(struct gkyl_range));
  gkyl_range_init_from_shape(geo->nrange, inp->cgrid->ndim, nodes);
  printf("nrange lower = %d %d %d\n", geo->nrange->lower[0], geo->nrange->lower[1], geo->nrange->lower[2]);
  printf("nrange upper = %d %d %d\n", geo->nrange->upper[0], geo->nrange->upper[1], geo->nrange->upper[2]);
  printf("cgrid ndim  = %d\n", inp->cgrid->ndim);
  printf("Checking the range volumes\n nrange.volume = %ld\n conversion_range.volume = %ld\n", geo->nrange->volume, conversion_range->volume);
  geo->mc2p_nodal_fd = gkyl_array_new(GKYL_DOUBLE, inp->cgrid->ndim*13, geo->nrange->volume);
  struct gkyl_array *mc2p = gkyl_array_new(GKYL_DOUBLE, inp->cgrid->ndim, geo->nrange->volume);

  enum { PH_IDX, AL_IDX, TH_IDX }; // arrangement of computational coordinates
  enum { X_IDX, Y_IDX, Z_IDX }; // arrangement of cartesian coordinates
  
  double dtheta = inp->cgrid->dx[TH_IDX],
    dpsi = inp->cgrid->dx[PH_IDX],
    dalpha = inp->cgrid->dx[AL_IDX];
  
  double theta_lo = inp->cgrid->lower[TH_IDX],
    phi_lo = inp->cgrid->lower[PH_IDX],
    alpha_lo = inp->cgrid->lower[AL_IDX];
  if(inp->bcs[1]==1)
    alpha_lo -= dalpha;

  double dx_fact = poly_order == 1 ? 1 : 0.5;
  dtheta *= dx_fact; dpsi *= dx_fact; dalpha *= dx_fact;

  // used for finite differences 
  double delta_alpha = dalpha*1e-4;
  double delta_psi = dpsi*1e-8;
  double delta_theta = dtheta*1e-4;
  geo->dzc = gkyl_malloc(3*sizeof(double));
  geo->dzc[0] = delta_psi;
  geo->dzc[1] = delta_alpha;
  geo->dzc[2] = delta_theta;
  int modifiers[5] = {0, -1, 1, -2, 2};

  double rclose = inp->rclose;
  double rright = inp->rright;
  double rleft = inp->rleft;
  printf("rleft, rright = %g, %g\n", rleft, rright);

  int nzcells = geo->rzgrid->cells[1];
  double *arc_memo = gkyl_malloc(sizeof(double[nzcells]));

  struct arc_length_ctx arc_ctx = {
    .geo = geo,
    .arc_memo = arc_memo,
    .ftype = inp->ftype,
    .zmaxis = inp->zmaxis
  };

  int cidx[3] = { 0 };
  for(int ia=geo->nrange->lower[AL_IDX]; ia<=geo->nrange->upper[AL_IDX]; ++ia){
    cidx[AL_IDX] = ia;
    for(int ia_delta = 0; ia_delta < 5; ia_delta++){ // should be <5
      if(ia == geo->nrange->lower[AL_IDX]){
        if(ia_delta == 1 || ia_delta == 3)
          continue; // want to use one sided stencils at edge
      }
      else if(ia == geo->nrange->upper[AL_IDX]){
          if(ia_delta == 2 || ia_delta == 4)
            continue; // want to use one sided stencils at edge
      }
      else{ //interior
        if( ia_delta == 3 || ia_delta == 4)
          continue; //dont do two away
      }

      double alpha_curr = alpha_lo + ia*dalpha + modifiers[ia_delta]*delta_alpha;
      printf("alpha_curr = %g\n", alpha_curr);

      for (int ip=geo->nrange->lower[PH_IDX]; ip<=geo->nrange->upper[PH_IDX]; ++ip) {
        int ip_delta_max = 5;// should be 5
        if(ia_delta != 0)
          ip_delta_max = 1;
        for(int ip_delta = 0; ip_delta < ip_delta_max; ip_delta++){
          if(ip == geo->nrange->lower[PH_IDX]){
            if(ip_delta == 1 || ip_delta == 3)
              continue; // want to use one sided stencils at edge
          }
          else if(ip == geo->nrange->upper[PH_IDX]){
            if(ip_delta == 2 || ip_delta == 4)
              continue; // want to use one sided stencils at edge
          }
          else{ // interior 
            if( ip_delta == 3 || ip_delta == 4)
              continue; //dont do two away
          }

          double zmin = inp->zmin, zmax = inp->zmax;
          double psi_curr = phi_lo + ip*dpsi + modifiers[ip_delta]*delta_psi;
          printf("psi_curr = %g\n", psi_curr);

          double arcL, darcL, arcL_curr, arcL_lo;
          double arcL_l, arcL_r;
          double phi_r, phi_l;

          if(inp->ftype == GKYL_CORE){
            //Find the turning point
            double zlo, zup, zlo_last;
            zlo = 0.01;
            zup=zmax;
            zlo_last = zlo;
            double R[4], dR[4];
            while(true){
              int nlo = R_psiZ(geo, psi_curr, zlo, 4, R, dR);
              //printf("zlo, zup nlo = %g, %g, %d\n",zlo,zup, nlo);
              if(nlo==2){
                if(fabs(zlo-zup)<1e-12){
                  printf("terminating, nlo = %d\n", nlo); 
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
            // lower one
            int nup = 0;
            double zup_last;
            zup = -0.01;
            zlo=zmin;
            zup_last = zup;
            while(true){
              int nup = R_psiZ(geo, psi_curr, zup, 4, R, dR);
              //printf("zlo, zup nlo = %g, %g, %d\n",zlo,zup, nlo);
              if(nup==2){
                if(fabs(zlo-zup)<1e-12){
                  printf("terminating, nup = %d\n", nup); 
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
            //end lower one
            printf("found zmin, zmax = %g, %g\n", zmin, zmax);
            printf("psi_curr, Z = %g, %1.16f\n", psi_curr, zlo);
            // Done finding turning point
            arcL_r = integrate_psi_contour_memo(geo, psi_curr, zmin, zmax, rright,
              true, true, arc_memo);
            arc_ctx.arcL_right = arcL_r;
            arc_ctx.right = false;
            arcL_l = integrate_psi_contour_memo(geo, psi_curr, zmin, zmax, rleft,
              true, true, arc_memo);
            arcL = arcL_l + arcL_r;
            printf("arcL total, L, R = %g, %g, %g\n", arcL, arcL_l, arcL_r);
            darcL = arcL/(poly_order*inp->cgrid->cells[TH_IDX]) * (inp->cgrid->upper[TH_IDX] - inp->cgrid->lower[TH_IDX])/2/M_PI;

            arc_ctx.right = true;
            arc_ctx.phi_right = 0.0;
            arc_ctx.rclose = rright;
            arc_ctx.psi = psi_curr;
            phi_r = phi_func(alpha_curr, zmax, &arc_ctx);
            arc_ctx.phi_right = phi_r - alpha_curr; // otherwise alpha will get added on twice
            printf("phi right = %g\n", arc_ctx.phi_right);
          }
          else if(inp->ftype == GKYL_PF_LO){
            //Find the  upper turning point
            double zlo, zup, zlo_last;
            zlo = zmin;
            zup=zmax;
            zlo_last = zlo;
            double R[4], dR[4];
            while(true){
              int nlo = R_psiZ(geo, psi_curr, zlo, 4, R, dR);
              if(nlo>=2){
                if(fabs(zlo-zup)<1e-12){
                  printf("terminating, nlo = %d\n", nlo); 
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
            printf("found zmin, zmax = %g, %g\n", zmin, zmax);
            printf("psi_curr, Z = %g, %1.16f\n", psi_curr, zlo);
            // Done finding turning point
            arcL_r = integrate_psi_contour_memo(geo, psi_curr, zmin, zmax, rright,
              true, true, arc_memo);
            arc_ctx.arcL_right = arcL_r;
            arc_ctx.right = false;
            arcL_l = integrate_psi_contour_memo(geo, psi_curr, zmin, zmax, rleft,
              true, true, arc_memo);
            arcL = arcL_l + arcL_r;
            printf("arcL total, L, R = %g, %g, %g\n", arcL, arcL_l, arcL_r);
            darcL = arcL/(poly_order*inp->cgrid->cells[TH_IDX]) * (inp->cgrid->upper[TH_IDX] - inp->cgrid->lower[TH_IDX])/2/M_PI;

            arc_ctx.right = true;
            arc_ctx.phi_right = 0.0;
            arc_ctx.rclose = rright;
            arc_ctx.psi = psi_curr;
            arc_ctx.zmin = zmin;
            phi_r = phi_func(alpha_curr, zmax, &arc_ctx);
            arc_ctx.phi_right = phi_r - alpha_curr; // otherwise alpha will get added on twice
            printf("phi right = %g\n", arc_ctx.phi_right);
          }
          else if(inp->ftype == GKYL_PF_UP){
            //Find the lower turning point
            double zlo, zup, zlo_last;
            double zup_last;
            zup = zmax;
            zlo=zmin;
            zup_last = zup;
            double R[4], dR[4];
            while(true){
              int nup = R_psiZ(geo, psi_curr, zup, 4, R, dR);
              //printf("zlo, zup nlo = %g, %g, %d\n",zlo,zup, nlo);
              if(nup>=2){
                if(fabs(zlo-zup)<1e-12){
                  printf("terminating, nup = %d\n", nup); 
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
            //end lower one
            printf("found zmin, zmax = %g, %g\n", zmin, zmax);
            printf("psi_curr, Z = %g, %1.16f\n", psi_curr, zlo);
            // Done finding turning point
            arcL_r = integrate_psi_contour_memo(geo, psi_curr, zmin, zmax, rright,
              true, true, arc_memo);
            arc_ctx.right = false;
            arcL_l = integrate_psi_contour_memo(geo, psi_curr, zmin, zmax, rleft,
              true, true, arc_memo);
            arc_ctx.arcL_left= arcL_l;
            arcL = arcL_l + arcL_r;
            printf("arcL total, L, R = %g, %g, %g\n", arcL, arcL_l, arcL_r);
            darcL = arcL/(poly_order*inp->cgrid->cells[TH_IDX]) * (inp->cgrid->upper[TH_IDX] - inp->cgrid->lower[TH_IDX])/2/M_PI;

            arc_ctx.right = false;
            arc_ctx.phi_left = 0.0;
            arc_ctx.rclose = rleft;
            arc_ctx.psi = psi_curr;
            arc_ctx.zmax = zmax;
            phi_l = phi_func(alpha_curr, zmin, &arc_ctx);
            arc_ctx.phi_left = phi_l - alpha_curr; // otherwise alpha will get added on twice
            printf("phi left= %g\n", arc_ctx.phi_left);
          }
          else if(inp->ftype==GKYL_SOL_DN_OUT){
            arc_ctx.phi_right = 0.0;
            arcL = integrate_psi_contour_memo(geo, psi_curr, zmin, zmax, rclose, true, true, arc_memo);
            darcL = arcL/(poly_order*inp->cgrid->cells[TH_IDX]) * (inp->cgrid->upper[TH_IDX] - inp->cgrid->lower[TH_IDX])/2/M_PI;
            printf("SOL_DN_OUT, arcL = %g\n", arcL);
          }
          else if(inp->ftype==GKYL_SOL_DN_IN){
            arc_ctx.phi_right = 0.0;
            arcL = integrate_psi_contour_memo(geo, psi_curr, zmin, zmax, rclose, true, true, arc_memo);
            darcL = arcL/(poly_order*inp->cgrid->cells[TH_IDX]) * (inp->cgrid->upper[TH_IDX] - inp->cgrid->lower[TH_IDX])/2/M_PI;
            printf("SOL_DN_IN, arcL = %g\n", arcL);
          }

          // at the beginning of each theta loop we need to reset things
          cidx[PH_IDX] = ip;
          arcL_curr = 0.0;
          arcL_lo = (theta_lo + M_PI)/2/M_PI*arcL;
          double ridders_min, ridders_max;
          // set node coordinates
          for (int it=geo->nrange->lower[TH_IDX]; it<=geo->nrange->upper[TH_IDX]; ++it) {
            int it_delta_max = 5; // should be 5
            if(ia_delta != 0 || ip_delta != 0 )
              it_delta_max = 1;
            for(int it_delta = 0; it_delta < it_delta_max; it_delta++){
              if(it == geo->nrange->lower[TH_IDX]){
                if(it_delta == 1 || it_delta == 3)
                  continue; // want to use one sided stencils at edge
              }
              else if(it == geo->nrange->upper[TH_IDX]){
                if(it_delta == 2 || it_delta == 4)
                  continue; // want to use one sided stencils at edge
              }
              else{
                if( it_delta == 3 || it_delta == 4)
                  continue; //dont do two away
              }
              arcL_curr = arcL_lo + it*darcL + modifiers[it_delta]*delta_theta*(arcL/2/M_PI);
              double theta_curr = arcL_curr*(2*M_PI/arcL) - M_PI ; // this is wrong need total arcL factor. Edit: 8/23 AS Not sure about this comment, shold have put a date in original. Seems to work fine.
              //printf("  it, theta_curr, arcL_curr = %d, %g %g\n", it, theta_curr, arcL_curr);
              //printf("  left and right arcL = %g %g\n", arcL_l, arcL_r);
              //printf("  rclose = %g\n", rclose);
              if(inp->ftype==GKYL_CORE){
                if(arcL_curr <= arcL_r){
                  rclose = rright;
                  arc_ctx.right = true;
                  ridders_min = -arcL_curr;
                  ridders_max = arcL-arcL_curr;
                }
                else{
                  rclose = rleft;
                  arc_ctx.right = false;
                  ridders_min = arcL - arcL_curr;
                  ridders_max = -arcL_curr + arc_ctx.arcL_right;
                }
              }
              if(inp->ftype==GKYL_PF_LO){
                if(arcL_curr <= arcL_r){
                  rclose = rright;
                  arc_ctx.right = true;
                  ridders_min = -arcL_curr;
                  ridders_max = arcL-arcL_curr;
                }
                else{
                  rclose = rleft;
                  arc_ctx.right = false;
                  ridders_min = arcL - arcL_curr;
                  ridders_max = -arcL_curr + arc_ctx.arcL_right;
                }
              }
              if(inp->ftype==GKYL_PF_UP){
                if(arcL_curr > arcL_l){
                  rclose = rright;
                  arc_ctx.right = true;
                  ridders_min = arc_ctx.arcL_left - arcL_curr;
                  ridders_max = arcL - arcL_curr;
                }
                else{
                  rclose = rleft;
                  arc_ctx.right = false;
                  ridders_min = arc_ctx.arcL_left - arcL_curr;
                  ridders_max = -arcL_curr;
                }
              }
              if(arc_ctx.ftype==GKYL_SOL_DN_OUT){
                ridders_min = -arcL_curr;
                ridders_max = arcL-arcL_curr;
                arc_ctx.right = false;
              }
              if(arc_ctx.ftype==GKYL_SOL_DN_IN){
                ridders_min = arcL-arcL_curr;
                ridders_max = -arcL_curr;
                arc_ctx.right = false;
              }

              arc_ctx.psi = psi_curr;
              arc_ctx.rclose = rclose;
              arc_ctx.zmin = zmin;
              arc_ctx.zmax = zmax;
              arc_ctx.arcL = arcL_curr;
              struct gkyl_qr_res res = gkyl_ridders(arc_length_func, &arc_ctx,
                zmin, zmax, ridders_min, ridders_max,
                geo->root_param.max_iter, 1e-10);
              double z_curr = res.res;
              ((gkyl_geo_gyrokinetic *)geo)->stat.nroot_cont_calls += res.nevals;
              double R[4] = { 0 }, dR[4] = { 0 };
              int nr = R_psiZ(geo, psi_curr, z_curr, 4, R, dR);
              double r_curr = choose_closest(rclose, R, R, nr);
              //printf("rcurr, zcurr = %g, %g\n\n", r_curr, z_curr);

              cidx[TH_IDX] = it;
              int lidx = 0;
              if (ip_delta != 0)
                lidx = 3 + 3*(ip_delta-1);
              if (ia_delta != 0)
                lidx = 15 + 3*(ia_delta-1);
              if (it_delta != 0)
                lidx = 27 + 3*(it_delta-1);

              double phi_curr = phi_func(alpha_curr, z_curr, &arc_ctx);
              //printf("ip,ia,it = %d, %d, %d\n", ip, ia, it);
              //printf("deltas : ip,ia,it = %d, %d, %d\n", ip_delta, ia_delta, it_delta);
              //printf("PHICURR = %g\n", phi_curr);
              // convert to x,y,z
              double *mc2p_fd_n = gkyl_array_fetch(geo->mc2p_nodal_fd, gkyl_range_idx(geo->nrange, cidx));
              double *mc2p_n = gkyl_array_fetch(mc2p, gkyl_range_idx(geo->nrange, cidx));
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
  gkyl_nodal_ops_n2m(inp->cbasis, inp->cgrid, geo->nrange, conversion_range, 3, mc2p, mapc2p, inp->bcs);

  char str1[50] = "xyz";
  char str2[50] = "allxyz";
  char str3[50] = "all";
  if (inp->write_node_coord_array){
    write_nodal_coordinates(strcat(str1, inp->node_file_nm), geo->nrange, mc2p);
    write_nodal_coordinates(strcat(str2, inp->node_file_nm), geo->nrange, geo->mc2p_nodal_fd);
    write_nodal_coordinates(strcat(str3, inp->node_file_nm), geo->nrange, mc2p);
  }
  gkyl_free(arc_memo);
  gkyl_array_release(mc2p);  
}


struct gkyl_geo_gyrokinetic_stat
gkyl_geo_gyrokinetic_get_stat(const gkyl_geo_gyrokinetic *geo)
{
  return geo->stat;
}

void
gkyl_geo_gyrokinetic_release(gkyl_geo_gyrokinetic *geo)
{
  gkyl_array_release(geo->psiRZ);
  gkyl_array_release(geo->mc2p_nodal_fd);
  gkyl_free(geo);
}


struct gkyl_range* gkyl_geo_gyrokinetic_get_nrange(gkyl_geo_gyrokinetic* geo){
  return geo->nrange;
}

struct gkyl_array* gkyl_geo_gyrokinetic_get_mc2p_nodal_fd(gkyl_geo_gyrokinetic* geo){
  return geo->mc2p_nodal_fd;
}

double* gkyl_geo_gyrokinetic_get_dzc(gkyl_geo_gyrokinetic* geo){
  return geo->dzc;
}
