#include <gkyl_tok_geo.h>

// Function context to pass to root finder
struct arc_length_ctx {
  const struct gkyl_tok_geo *geo;
  double *arc_memo;
  double *arc_memo_left;
  double *arc_memo_right;
  double psi, rclose, zmin, arcL;
  double rleft, rright, zmax;
  double zmin_left, zmin_right; // for single null full SOL only
  double arcL_right; // this is for when we need to switch sides
  double arcL_left; // this is for when we need to switch sides
  double arcL_tot; // total arc length
  double phi_right; // this is for when we need to switch sides
  double phi_left; // this is for when we need to switch sides
  double phi_bot; // For new way of trying to do core
  bool right;
  double zmaxis;
  enum gkyl_tok_geo_type ftype; // type of geometry
};


// Context to pass to endpoint finder
struct plate_ctx{
  const struct gkyl_tok_geo* geo;
  double psi_curr;
  bool lower;
};

// some helper functions
static double
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

// Compute roots R(psi,Z) and dR/dZ(psi,Z) in a p=2 DG cell with tensor basis
static inline struct RdRdZ_sol
calc_RdR_p2_tensor(const double *psi, double psi0, double Z, double xc[2], double dx[2])
{
  struct RdRdZ_sol sol = { .nsol = 0 };
  double y = (Z-xc[1])/(dx[1]*0.5);

  double aq = 0.125*(45.0*psi[8]*SQ(y)+23.2379000772445*psi[6]*y-15.0*psi[8]+13.41640786499874*psi[4]);
  double bq = 0.125*(23.2379000772445*psi[7]*SQ(y)+12.0*psi[3]*y-7.745966692414834*psi[7]+6.928203230275509*psi[1]) ;
  double cq = 0.125*((13.41640786499874*psi[5]-15.0*psi[8])*SQ(y)+(6.928203230275509*psi[2]-7.745966692414834*psi[6])*y+5.0*psi[8]- 4.47213595499958*psi[5]-4.47213595499958*psi[4]+4.0*psi[0] ) - psi0;

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
      double C = 0.125*(SQ(x)*(90.0*psi[8]*y+23.2379000772445*psi[6])+x*(46.47580015448901*psi[7]*y+12.0*psi[3])+2* (13.41640786499874*psi[5]-15.0*psi[8])*y-7.745966692414834*psi[6]+6.928203230275509*psi[2]) ;
      double A = 0.125*(2*x*(45.0*psi[8]*SQ(y)+23.2379000772445*psi[6]*y-15.0*psi[8]+13.41640786499874*psi[4])+23.2379000772445*psi[7]*SQ(y)+12.0*psi[3]*y-7.745966692414834*psi[7]+6.928203230275509*psi[1]); 
      sol.dRdZ[sidx] = -C/A*dx[0]/dx[1];
      
      sidx += 1;
    }
    if ((-1<=r2) && (r2 < 1)) {
      sol.nsol += 1;
      sol.R[sidx] = r2*dx[0]*0.5 + xc[0];

      double x = r2;
      double C = 0.125*(SQ(x)*(90.0*psi[8]*y+23.2379000772445*psi[6])+x*(46.47580015448901*psi[7]*y+12.0*psi[3])+2* (13.41640786499874*psi[5]-15.0*psi[8])*y-7.745966692414834*psi[6]+6.928203230275509*psi[2]) ;
      double A = 0.125*(2*x*(45.0*psi[8]*SQ(y)+23.2379000772445*psi[6]*y-15.0*psi[8]+13.41640786499874*psi[4])+23.2379000772445*psi[7]*SQ(y)+12.0*psi[3]*y-7.745966692414834*psi[7]+6.928203230275509*psi[1]); 
      sol.dRdZ[sidx] = -C/A*dx[0]/dx[1];
      
      sidx += 1;
    }
  }
  return sol;
}

// Compute roots R(psi,Z) and dR/dZ(psi,Z) in a p=2 DG cell with tensor basis
// Use more accurate roots from numerical recipes in C 2007 section 5.6
static inline struct RdRdZ_sol
calc_RdR_p2_tensor_nrc(const double *psi, double psi0, double Z, double xc[2], double dx[2])
{
  struct RdRdZ_sol sol = { .nsol = 0 };
  double y = (Z-xc[1])/(dx[1]*0.5);

  double aq = 0.125*(45.0*psi[8]*SQ(y)+23.2379000772445*psi[6]*y-15.0*psi[8]+13.41640786499874*psi[4]);
  double bq = 0.125*(23.2379000772445*psi[7]*SQ(y)+12.0*psi[3]*y-7.745966692414834*psi[7]+6.928203230275509*psi[1]) ;
  double cq = 0.125*((13.41640786499874*psi[5]-15.0*psi[8])*SQ(y)+(6.928203230275509*psi[2]-7.745966692414834*psi[6])*y+5.0*psi[8]- 4.47213595499958*psi[5]-4.47213595499958*psi[4]+4.0*psi[0] ) - psi0;

  double delta2 = bq*bq - 4*aq*cq;

  if (delta2 > 0) {
    double r1, r2;
    double delta = sqrt(delta2);
    //// compute both roots
    double qq = -0.5*(bq + (bq/fabs(bq)) * delta);
    r1 = qq/aq;
    r2 = cq/qq;

    int sidx = 0;
    if ((-1<=r1) && (r1 < 1)) {
      sol.nsol += 1;
      sol.R[sidx] = r1*dx[0]*0.5 + xc[0];

      double x = r1;
      double C = 0.125*(SQ(x)*(90.0*psi[8]*y+23.2379000772445*psi[6])+x*(46.47580015448901*psi[7]*y+12.0*psi[3])+2* (13.41640786499874*psi[5]-15.0*psi[8])*y-7.745966692414834*psi[6]+6.928203230275509*psi[2]) ;
      double A = 0.125*(2*x*(45.0*psi[8]*SQ(y)+23.2379000772445*psi[6]*y-15.0*psi[8]+13.41640786499874*psi[4])+23.2379000772445*psi[7]*SQ(y)+12.0*psi[3]*y-7.745966692414834*psi[7]+6.928203230275509*psi[1]); 
      sol.dRdZ[sidx] = -C/A*dx[0]/dx[1];
      
      sidx += 1;
    }
    if ((-1<=r2) && (r2 < 1)) {
      sol.nsol += 1;
      sol.R[sidx] = r2*dx[0]*0.5 + xc[0];

      double x = r2;
      double C = 0.125*(SQ(x)*(90.0*psi[8]*y+23.2379000772445*psi[6])+x*(46.47580015448901*psi[7]*y+12.0*psi[3])+2* (13.41640786499874*psi[5]-15.0*psi[8])*y-7.745966692414834*psi[6]+6.928203230275509*psi[2]) ;
      double A = 0.125*(2*x*(45.0*psi[8]*SQ(y)+23.2379000772445*psi[6]*y-15.0*psi[8]+13.41640786499874*psi[4])+23.2379000772445*psi[7]*SQ(y)+12.0*psi[3]*y-7.745966692414834*psi[7]+6.928203230275509*psi[1]); 
      sol.dRdZ[sidx] = -C/A*dx[0]/dx[1];
      
      sidx += 1;
    }
  }
  return sol;
}

// Compute roots R(psi,Z) and dR/dZ(psi,Z) in a p=2 DG cell with tensor basis if delta2 is negative byt very small
// Use more accurate roots from numerical recipes in C 2007 section 5.6
static inline struct RdRdZ_sol
calc_RdR_p2_tensor_nrc_none(const double *psi, double psi0, double Z, double xc[2], double dx[2])
{
  struct RdRdZ_sol sol = { .nsol = 0 };
  double y = (Z-xc[1])/(dx[1]*0.5);

  double aq = 0.125*(45.0*psi[8]*SQ(y)+23.2379000772445*psi[6]*y-15.0*psi[8]+13.41640786499874*psi[4]);
  double bq = 0.125*(23.2379000772445*psi[7]*SQ(y)+12.0*psi[3]*y-7.745966692414834*psi[7]+6.928203230275509*psi[1]) ;
  double cq = 0.125*((13.41640786499874*psi[5]-15.0*psi[8])*SQ(y)+(6.928203230275509*psi[2]-7.745966692414834*psi[6])*y+5.0*psi[8]- 4.47213595499958*psi[5]-4.47213595499958*psi[4]+4.0*psi[0] ) - psi0;

  double delta2 = bq*bq - 4*aq*cq;
  if(delta2 > 0)
    return sol;

  if (fabs(delta2) < 1e-8) {
    double r1, r2;
    double delta = 0.0;
    //// compute both roots
    double qq = -0.5*bq;
    r1 = qq/aq;
    r2 = cq/qq;
    // The two roots should really both be equal to -sqrt(c/a) but numerically they are quite different
    // It seems that only the expression for r2 = cq/aq is robust when both c and a are very small

    int sidx = 0;
    //if ((-1<=r1) && (r1 < 1)) {
    //  sol.nsol += 1;
    //  sol.R[sidx] = r1*dx[0]*0.5 + xc[0];

    //  double x = r1;
    //  double C = 0.125*(SQ(x)*(90.0*psi[8]*y+23.2379000772445*psi[6])+x*(46.47580015448901*psi[7]*y+12.0*psi[3])+2* (13.41640786499874*psi[5]-15.0*psi[8])*y-7.745966692414834*psi[6]+6.928203230275509*psi[2]) ;
    //  double A = 0.125*(2*x*(45.0*psi[8]*SQ(y)+23.2379000772445*psi[6]*y-15.0*psi[8]+13.41640786499874*psi[4])+23.2379000772445*psi[7]*SQ(y)+12.0*psi[3]*y-7.745966692414834*psi[7]+6.928203230275509*psi[1]); 
    //  sol.dRdZ[sidx] = -C/A*dx[0]/dx[1];
    //  
    //  sidx += 1;
    //}
    if ((-1<=r2) && (r2 < 1)) {
      sol.nsol += 1;
      sol.R[sidx] = r2*dx[0]*0.5 + xc[0];

      double x = r2;
      double C = 0.125*(SQ(x)*(90.0*psi[8]*y+23.2379000772445*psi[6])+x*(46.47580015448901*psi[7]*y+12.0*psi[3])+2* (13.41640786499874*psi[5]-15.0*psi[8])*y-7.745966692414834*psi[6]+6.928203230275509*psi[2]) ;
      double A = 0.125*(2*x*(45.0*psi[8]*SQ(y)+23.2379000772445*psi[6]*y-15.0*psi[8]+13.41640786499874*psi[4])+23.2379000772445*psi[7]*SQ(y)+12.0*psi[3]*y-7.745966692414834*psi[7]+6.928203230275509*psi[1]); 
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
R_psiZ(const struct gkyl_tok_geo *geo, double psi, double Z, int nmaxroots,
  double *R, double *dR)
{
  int zcell = get_idx(1, Z, &geo->rzgrid, &geo->rzlocal);

  int sidx = 0;
  int idx[2] = { 0, zcell };
  double dx[2] = { geo->rzgrid.dx[0], geo->rzgrid.dx[1] };
  
  struct gkyl_range rangeR;
  gkyl_range_deflate(&rangeR, &geo->rzlocal, (int[]) { 0, 1 }, (int[]) { 0, zcell });

  struct gkyl_range_iter riter;
  gkyl_range_iter_init(&riter, &rangeR);
  
  // loop over all R cells to find psi crossing
  while (gkyl_range_iter_next(&riter) && sidx<=nmaxroots) {
    long loc = gkyl_range_idx(&rangeR, riter.idx);
    const double *psih = gkyl_array_cfetch(geo->psiRZ, loc);

    double xc[2];
    idx[0] = riter.idx[0];
    gkyl_rect_grid_cell_center(&geo->rzgrid, idx, xc);

    struct RdRdZ_sol sol = geo->calc_roots(psih, psi, Z, xc, dx);
    
    if (sol.nsol > 0)
      for (int s=0; s<sol.nsol; ++s) {
        if( (sol.R[s] > geo->rmin) && (sol.R[s] < geo->rmax) ) {
          R[sidx] = sol.R[s];
          dR[sidx] = sol.dRdZ[s];
          sidx += 1;
        }
      }
  }

  // Try again if we didn't find any
  //if (sidx==0 && geo->tol_no_roots) {
  if (sidx==0 && geo->tol_no_roots) {
    gkyl_range_iter_init(&riter, &rangeR);
    while (gkyl_range_iter_next(&riter) && sidx<=nmaxroots) {
      long loc = gkyl_range_idx(&rangeR, riter.idx);
      const double *psih = gkyl_array_cfetch(geo->psiRZ, loc);

      double xc[2];
      idx[0] = riter.idx[0];
      gkyl_rect_grid_cell_center(&geo->rzgrid, idx, xc);

      struct RdRdZ_sol sol = calc_RdR_p2_tensor_nrc_none(psih, psi, Z, xc, dx);
      
      if (sol.nsol > 0)
        for (int s=0; s<sol.nsol; ++s) {
          if( (sol.R[s] > geo->rmin) && (sol.R[s] < geo->rmax) ) {
            R[sidx] = sol.R[s];
            dR[sidx] = sol.dRdZ[s];
            sidx += 1;
          }
        }
    }
  }

  return sidx;
}


// Function context to pass to coutour integration function
struct contour_ctx {
  const struct gkyl_tok_geo *geo;
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
  iter.idx[0] = fmin(c->geo->rzlocal.lower[0] + (int) floor((r_curr - c->geo->rzgrid.lower[0])/c->geo->rzgrid.dx[0]), c->geo->rzlocal.upper[0]);
  iter.idx[1] = fmin(c->geo->rzlocal.lower[1] + (int) floor((Z - c->geo->rzgrid.lower[1])/c->geo->rzgrid.dx[1]), c->geo->rzlocal.upper[1]);
  long loc = gkyl_range_idx((&c->geo->rzlocal), iter.idx);
  const double *psih = gkyl_array_cfetch(c->geo->psiRZ, loc);

  double xc[2];
  gkyl_rect_grid_cell_center((&c->geo->rzgrid), iter.idx, xc);
  double x = (r_curr-xc[0])/(c->geo->rzgrid.dx[0]*0.5);
  double y = (Z-xc[1])/(c->geo->rzgrid.dx[1]*0.5);

  // if psi is polyorder 2 we can get grad psi
  // in cylindrical coords it is grad psi = dpsi/dR Rhat + dpsi/dZ zhat
  double dpsidx = 2.904737509655563*psih[7]*(y*y-0.3333333333333333)+5.809475019311126*psih[6]*x*y+1.5*psih[3]*y+3.354101966249684*psih[4]*x+0.8660254037844386*psih[1]; 
  double dpsidy =	5.809475019311126*psih[7]*x*y+3.354101966249684*psih[5]*y+2.904737509655563*psih[6]*(x*x-0.3333333333333333)+1.5*psih[3]*x+0.8660254037844386*psih[2];
  dpsidx = dpsidx*2.0/c->geo->rzgrid.dx[0];
  dpsidy = dpsidy*2.0/c->geo->rzgrid.dx[1];
  double grad_psi_mag = sqrt(dpsidx*dpsidx + dpsidy*dpsidy);
  double result  = (1/r_curr/grad_psi_mag) *sqrt(1+dRdZ*dRdZ) ;
  return nr>0 ? result : 0.0;
}

static inline double
dphidtheta_integrand(double Z, void *ctx)
{
  struct contour_ctx *c = ctx;
  c->ncall += 1;
  double R[4] = { 0 }, dR[4] = { 0 };
  
  int nr = R_psiZ(c->geo, c->psi, Z, 4, R, dR);
  double dRdZ = nr == 1 ? dR[0] : choose_closest(c->last_R, R, dR, nr);
  double r_curr = nr == 1 ? R[0] : choose_closest(c->last_R, R, R, nr);

  struct gkyl_range_iter iter;
  iter.idx[0] = fmin(c->geo->rzlocal.lower[0] + (int) floor((r_curr - c->geo->rzgrid.lower[0])/c->geo->rzgrid.dx[0]), c->geo->rzlocal.upper[0]);
  iter.idx[1] = fmin(c->geo->rzlocal.lower[1] + (int) floor((Z - c->geo->rzgrid.lower[1])/c->geo->rzgrid.dx[1]), c->geo->rzlocal.upper[1]);
  long loc = gkyl_range_idx((&c->geo->rzlocal), iter.idx);
  const double *psih = gkyl_array_cfetch(c->geo->psiRZ, loc);

  double xc[2];
  gkyl_rect_grid_cell_center((&c->geo->rzgrid), iter.idx, xc);
  double x = (r_curr-xc[0])/(c->geo->rzgrid.dx[0]*0.5);
  double y = (Z-xc[1])/(c->geo->rzgrid.dx[1]*0.5);

  // if psi is polyorder 2 we can get grad psi
  // in cylindrical coords it is grad psi = dpsi/dR Rhat + dpsi/dZ zhat
  double dpsidx = 2.904737509655563*psih[7]*(y*y-0.3333333333333333)+5.809475019311126*psih[6]*x*y+1.5*psih[3]*y+3.354101966249684*psih[4]*x+0.8660254037844386*psih[1]; 
  double dpsidy =	5.809475019311126*psih[7]*x*y+3.354101966249684*psih[5]*y+2.904737509655563*psih[6]*(x*x-0.3333333333333333)+1.5*psih[3]*x+0.8660254037844386*psih[2];
  dpsidx = dpsidx*2.0/c->geo->rzgrid.dx[0];
  dpsidy = dpsidy*2.0/c->geo->rzgrid.dx[1];
  double grad_psi_mag = sqrt(dpsidx*dpsidx + dpsidy*dpsidy);
  double result  = (1/r_curr/grad_psi_mag);
  return nr>0 ? result : 0.0;
}

// Integrates along a specified contour, optionally using a "memory"
// of previously computed values, or storing computed values in
// memory. The function basically breaks up the integral into a loop
// over z-cells. This needs to be done as the DG representation is,
// well, discontinuous, and adaptive quadrature struggles with such
// functions.
static double
integrate_psi_contour_memo(const struct gkyl_tok_geo *geo, double psi,
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
  
  double dz = geo->rzgrid.dx[1];
  double zlo = geo->rzgrid.lower[1];
  int izlo = geo->rzlocal.lower[1], izup = geo->rzlocal.upper[1];
  
  int ilo = get_idx(1, zmin, &geo->rzgrid, &geo->rzlocal);
  int iup = get_idx(1, zmax, &geo->rzgrid, &geo->rzlocal);

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

  ((struct gkyl_tok_geo *)geo)->stat.nquad_cont_calls += ctx.ncall;
  return res;
}

static double
integrate_phi_along_psi_contour_memo(const struct gkyl_tok_geo *geo, double psi,
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
  
  double dz = geo->rzgrid.dx[1];
  double zlo = geo->rzgrid.lower[1];
  int izlo = geo->rzlocal.lower[1], izup = geo->rzlocal.upper[1];
  
  int ilo = get_idx(1, zmin, &geo->rzgrid, &geo->rzlocal);
  int iup = get_idx(1, zmax, &geo->rzgrid, &geo->rzlocal);

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

  ((struct gkyl_tok_geo *)geo)->stat.nquad_cont_calls += ctx.ncall;
  return res;
}


double phi_func(double alpha_curr, double Z, void *ctx);
double tok_plate_psi_func(double s, void *ctx);

/*
 * Used to set zmin and zmax and attributes of arc_ctx before looping over arc length
*/
void tok_find_endpoints(struct gkyl_tok_geo_grid_inp* inp, struct gkyl_tok_geo *geo, struct arc_length_ctx* arc_ctx, struct plate_ctx* pctx, double psi_curr, double alpha_curr, double* arc_memo, double* arc_memo_left, double* arc_memo_right);


/*
 * Used to set arc_ctx attributes before using ridders to find z
*/
void tok_set_ridders(struct gkyl_tok_geo_grid_inp* inp, struct arc_length_ctx* arc_ctx, double psi_curr, double arcL_curr, double* rclose, double *ridders_min, double* ridders_max);

