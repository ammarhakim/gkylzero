#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_wv_cold_sr_fluid.h>

// Files for polynomial solve
#include <gkyl_math.h>

#define NUX 1
#define NUY 2
#define NUZ 3

#define SQ(x) ((x)*(x))

struct wv_cold_sr_fluid {
  struct gkyl_wv_eqn eqn; // base object
  double clight;
};

// relativistic case
struct cold_sr {
  double rhol, rhor;
  double vl[3], vr[3];
};

void
cold_sr_fluid_flux(const double q[4], double *flux)
{
  // Vx = NUx/sqrt(N^2 + NU^2/c^2) 
  const double c = 1.0; //c = 299792458.0; 
  //const double c = 299792458.0; 
  printf("TEMP SPEED IN COLD_SR_FLUX\n");
  double Vx = q[NUX]/sqrt(q[0]*q[0] + (q[NUX]*q[NUX] + q[NUY]*q[NUY] + q[NUZ]*q[NUZ])/(c*c)); 
  flux[0] = q[0]*Vx; // N*Vx
  flux[NUX] = q[NUX]*Vx; // N*Ux*Vx
  flux[NUY] = q[NUY]*Vx; // N*Uy*Vx
  flux[NUZ] = q[NUZ]*Vx; // N*Uz*Vx
}

static inline void
cold_sr_fluid_cons_to_diag(const struct gkyl_wv_eqn *eqn,
  const double *qin, double *diag)
{
  // density and moment is copied as-is
  for (int i=0; i<4; ++i) diag[i] = qin[i];
  // Ke-sr density (gamma-1)
  const double c = 299792458.0;
  double ke = sqrt(qin[0]*qin[0] + (qin[1]*qin[1] + qin[2]*qin[2] + qin[3]*qin[3])/(c*c))/qin[0] - 1.0; 
  diag[4] = ke;
}

static void
cold_sr_fluid_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_wv_eqn *base = container_of(ref, struct gkyl_wv_eqn, ref_count);
  struct wv_cold_sr_fluid *cold_sr_fluid = container_of(base, struct wv_cold_sr_fluid, eqn);
  gkyl_free(cold_sr_fluid);
}

static inline void
cons_to_riem_sr(const struct gkyl_wv_eqn *eqn,
  const double *qstate, const double *qin, double *wout)
{
  // TODO: this should use proper L matrix
  for (int i=0; i<4; ++i)
    wout[i] = qin[i];
}
static inline void
riem_to_cons_sr(const struct gkyl_wv_eqn *eqn,
  const double *qstate, const double *win, double *qout)
{
  // TODO: this should use proper L matrix
  for (int i=0; i<4; ++i)
    qout[i] = win[i];
}

static inline void
rot_to_local(const double *tau1, const double *tau2, const double *norm,
  const double *GKYL_RESTRICT qglobal, double *GKYL_RESTRICT qlocal)
{
  qlocal[0] = qglobal[0];
  qlocal[1] = qglobal[1]*norm[0] + qglobal[2]*norm[1] + qglobal[3]*norm[2];
  qlocal[2] = qglobal[1]*tau1[0] + qglobal[2]*tau1[1] + qglobal[3]*tau1[2];
  qlocal[3] = qglobal[1]*tau2[0] + qglobal[2]*tau2[1] + qglobal[3]*tau2[2];
}

static inline void
rot_to_global(const double *tau1, const double *tau2, const double *norm,
  const double *GKYL_RESTRICT qlocal, double *GKYL_RESTRICT qglobal)
{
  qglobal[0] = qlocal[0];
  qglobal[1] = qlocal[1]*norm[0] + qlocal[2]*tau1[0] + qlocal[3]*tau2[0];
  qglobal[2] = qlocal[1]*norm[1] + qlocal[2]*tau1[1] + qlocal[3]*tau2[1];
  qglobal[3] = qlocal[1]*norm[2] + qlocal[2]*tau1[2] + qlocal[3]*tau2[2];
}

static inline void
calc_ufromv(const double v[3], double *u)
{
  double v2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
  double gamma1 = 1/sqrt(1-v2);
  u[0] = v[0]*gamma1;
  u[1] = v[1]*gamma1;
  u[2] = v[2]*gamma1;
}

// returns 0 if it passes, 1 if fails
static bool 
entropy_check(const double UL[3], const double UR[3], double u[3])
{
  bool cond = 1;
  double tol = 1.e-15;
  for(int i=0; i<3; ++i){
    cond = cond && (fmin(UR[i],UL[i]) - tol) <= u[i] && u[i] <= (fmax(UR[i],UL[i]) + tol);
  }
  return !cond;
}

// returns 0 if it passes, 1 if fails, strict (no wiggle for overshoot)
static bool 
entropy_check_strict(const double UL[3], const double UR[3], double u[3])
{
  bool cond = 1;
  for(int i=0; i<3; ++i){
    cond = cond && (fmin(UR[i],UL[i])) <= u[i] && u[i] <= (fmax(UR[i],UL[i]));
  }
  return !cond;
}


static void
calc_u_from_result(double *uhat, double v_hat, void *ctx)
{

  // compute the coefficients
  struct cold_sr *cnr = ctx;
  double rhol = cnr->rhol, rhor = cnr->rhor;
  double vl[3], vr[3];

  for (int i=0; i<3; ++i) {
    vl[i] = cnr->vl[i];
    vr[i] = cnr->vr[i];
  }

  double ul[3], ur[3];
  calc_ufromv(vl, ul);
  calc_ufromv(vr, ur);
  double a, drho, A[3], B[3];

  for (int i=0; i<3; ++i) {
    A[i] = -rhol*ul[i]*vl[0] + rhor*ur[i]*vr[0];
    B[i] = rhol*ul[i] - rhor*ur[i];
  }

  drho = rhor-rhol;
  a = rhol*vl[0] - rhor*vr[0];

  // Compute the resulting U:
  for (int i=0; i<3; ++i) {
    uhat[i] = -(B[i]*v_hat+A[i])/(drho*v_hat+a);
  }
}

static inline void
sr_poly_coeffs(void *ctx, double *coeff)
{

  // compute the coefficients
  struct cold_sr *cnr = ctx;
  double rhol = cnr->rhol, rhor = cnr->rhor;
  double vl[3], vr[3];

  for (int i=0; i<3; ++i) {
    vl[i] = cnr->vl[i];
    vr[i] = cnr->vr[i];
  }

  double ul[3], ur[3];
  calc_ufromv(vl, ul);
  calc_ufromv(vr, ur);
  double a, drho, A[3], B[3];

  for (int i=0; i<3; ++i) {
    A[i] = -rhol*ul[i]*vl[0] + rhor*ur[i]*vr[0];
    B[i] = rhol*ul[i] - rhor*ur[i];
  }

  drho = rhor-rhol;
  a = rhol*vl[0] - rhor*vr[0];

  double p4coeff = SQ(drho)+SQ(B[2])+SQ(B[1])+SQ(B[0]); 
  double p3coeff = 2.0*a*drho+2.0*A[2]*B[2]+2.0*A[1]*B[1]+2.0*A[0]*B[0]; 
  double p2coeff = SQ(a)+SQ(A[2])+SQ(A[1])-SQ(B[0])+SQ(A[0]); 
  double p1coeff = -2.0*A[0]*B[0]; 
  double p0coeff = -SQ(A[0]); 

  // Find the four roots to the polynomial Vx equation:
  // f(x)=x^{4}+ax^{3}+bx^{2}+cx+d,
  // via durand-kerner: https://en.wikipedia.org/wiki/Durand–Kerner_method
  coeff[3] = p3coeff/p4coeff;
  coeff[2] = p2coeff/p4coeff;
  coeff[1] = p1coeff/p4coeff;
  coeff[0] = p0coeff/p4coeff;

}


static double 
root_select(struct gkyl_root_intervals *roots, double ul[3], double ur[3], int *status,
void *ctx)
{

  // status and number of roots
  status[0] = 0; //Passed 
  int num_roots = roots->nroots;
  double chosen_velocity = 1.1;
  int roots_satisfying_cond = 0;

  // Stop on problem:
  if (num_roots <= 0 || num_roots > 4){
    //printf("Issue in root select: num_roots is invalid %d\n",num_roots);
    status[0] = 1;

  } else if (roots->status){
    //printf("Root_select failed with status 1 %d\n",num_roots);
    status[0] = 2;

  } else {
    // Iterate over each root
    for (int i=0; i<num_roots; ++i) {

      // Start with zero root
      double root = 0.0;
      
      // Select the root
      if (roots->status_ridders[i] == 0){
        root = roots->real_roots_ridders[i];  
        status[0] = 0;
      // else replace with average of right and left bounds after refinement 
      //} else if (roots->status_refinement[i] == 0) {
      //  root = roots->root_bound_lower[i]
      //  + (roots->root_bound_upper[i]-roots->root_bound_lower[i])/2.0;
      //  status[0] = 1; 
      } else {
        status[0] = 3;
        //printf("Failure to isolate a specifc root\n");
      }

      // See if the solution is the right one
      double u_root[3];
      bool entropy_test_passed;
      if (status[0] == 0){

        // Reconstruct U_root 
        calc_u_from_result(u_root,root,ctx);

        // Test if U_root satisfied UL_i < U_root_i < UR_i
        bool entropy_test_passed = entropy_check(ul,ur,u_root);

        // If the entropy test is passed, then this is the right root
        if (entropy_test_passed == 0){
          chosen_velocity = root;
          roots_satisfying_cond = roots_satisfying_cond + 1;
        } else if (roots_satisfying_cond == 0 && entropy_test_passed == 1) {
          status[0] = 4; 
        }
      }

      // Roots satisfying conditions
      if (roots_satisfying_cond < 1){
        status[0] = 5; 
        //printf("Invalid number of roots satisfying entropy condition: %d\n",roots_satisfying_cond);
      } else if (roots_satisfying_cond > 1) { 

        // isolate the right root if multiple satisfy the condition
        roots_satisfying_cond = 0;
        for (int i=0; i<num_roots; ++i) {
          if (roots->status_ridders[i] == 0){
            root = roots->real_roots_ridders[i];  
            double u_root[3];
            calc_u_from_result(u_root,root,ctx);
            bool entropy_test_strict_passed = entropy_check_strict(ul,ur,u_root);
            if (entropy_test_strict_passed == 0){
              chosen_velocity = root;
              roots_satisfying_cond = roots_satisfying_cond + 1;
            } else if (roots_satisfying_cond == 0 && entropy_test_strict_passed == 1) {
              status[0] = 6; 
            }
          }
        }
      }
    } // end iterating over roots
  } // end if

  // Return the root which best matches
  return chosen_velocity;
}

double 
compute_sr_roe_averaged_velocity_cold_limit(const double ql[4], const double qr[4], const double c, double *u, double *gamma_avg) 
{

  // As long as one density is positive on both sides:
  if (ql[0] > 0.0 || qr[0] > 0.0){

    // output 
    double v_hat;

    // isolate variables (right/left), normalize, u & v by c
    double rhor = qr[0];
    double urx = qr[1]/(c*qr[0]);
    double ury = qr[2]/(c*qr[0]);
    double urz = qr[3]/(c*qr[0]);
    double rhol = ql[0];
    double ulx = ql[1]/(c*ql[0]);
    double uly = ql[2]/(c*ql[0]);
    double ulz = ql[3]/(c*ql[0]);

    // compute the constants:
    double gammal = sqrt(1.0 + (ulx*ulx + uly*uly + ulz*ulz));
    double gammar = sqrt(1.0 + (urx*urx + ury*ury + urz*urz));
    double vlx = ulx/gammal;
    double vrx = urx/gammar;
    double vly = uly/gammal;
    double vry = ury/gammar;
    double vlz = ulz/gammal;
    double vrz = urz/gammar;

    //compute the primative vars
    double rhol_prim = rhol/gammal;
    double rhor_prim = rhor/gammar;

    // Compute the primative-parameterization state vector w
    // these are the averages of the left and right states
    double k = sqrt(rhol_prim) + sqrt(rhor_prim);
    double w0 = sqrt(rhol_prim)*gammal + sqrt(rhor_prim)*gammar;
    double w1 = sqrt(rhol_prim)*gammal*vlx + sqrt(rhor_prim)*gammar*vrx;
    double w2 = sqrt(rhol_prim)*gammal*vly + sqrt(rhor_prim)*gammar*vry;
    double w3 = sqrt(rhol_prim)*gammal*vlz + sqrt(rhor_prim)*gammar*vrz;

    // Compute F^0 which is in terms of k, w state (These are our new conserved variables)
    // q_avg = [N, NUx, NUy, Nuz]
    double q_avg[4];
    q_avg[0] = k*w0;
    q_avg[1] = w0*w1;
    q_avg[2] = w0*w2;
    q_avg[3] = w0*w3;

    // Recover u and v
    u[0] = q_avg[1]/q_avg[0];
    u[1] = q_avg[2]/q_avg[0];
    u[2] = q_avg[3]/q_avg[0];


    // convert back from u/c -> u
    v_hat = u[0]/sqrt(1.0 + u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
    double other_v1 = w1/w0;
    double other_v2 = 2*w0*w1/( k*k + w0*w0 + w1*w1 + w2*w2 + w3*w3 );
    double other_v3 = u[0]/sqrt(1.0 + u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);;
    printf("Choice, Vx = w1/w0: %1.16e or cons-reconstruction: %1.16e, 2*w0*w1/...: %1.16e\n",other_v1, other_v3, other_v2);
    u[0] = c*u[0];
    u[1] = c*u[1];
    u[2] = c*u[2];

    // return the velocity
    return v_hat*c;

  } else {
    return 0.0;
  }
}

double
compute_sr_roe_averaged_velocity(const double ql[4], const double qr[4], const double c, double *u) 
{

  // output 
  double v_hat;

  // Isolate variables (right/left), normalize, u & v by c
  double rhoR = qr[0];
  double uRx = qr[1]/(c*qr[0]);
  double uRy = qr[2]/(c*qr[0]);
  double uRz = qr[3]/(c*qr[0]);
  double rhoL = ql[0];
  double uLx = ql[1]/(c*ql[0]);
  double uLy = ql[2]/(c*ql[0]);
  double uLz = ql[3]/(c*ql[0]);

  // Compute the constants:
  double gammaL = sqrt(1.0 + (uLx*uLx + uLy*uLy + uLz*uLz));
  double gammaR = sqrt(1.0 + (uRx*uRx + uRy*uRy + uRz*uRz));
  double vLx = uLx/gammaL;
  double vRx = uRx/gammaR;
  double vLy = uLy/gammaL;
  double vRy = uRy/gammaR;
  double vLz = uLz/gammaL;
  double vRz = uRz/gammaR;

  // Normalize density, verify positivity/nonzero:
  if (ql[0] > 0.0 && qr[0] > 0.0){

    // Normalize density by dividing by rhoL
    rhoL = 1.0;
    rhoR = qr[0]/ql[0];

    struct cold_sr csr = {
      .rhol = rhoL,
      .rhor = rhoR,
      .vl = { vLx, vLy, vLz },
      .vr = { vRx, vRy, vRz }
    };

    // Compute a default answer, the v-hat-avg, this is a worst case-everything else failed
    double u_avg[3] = {(uLx+uRx)/2.0,(uLy+uRy)/2.0,(uLz+uRz)/2.0};
    double v_hat_avg = u_avg[0]/sqrt(1.0 + u_avg[0]*u_avg[0] + u_avg[1]*u_avg[1] + u_avg[2]*u_avg[2] );

    // Check the two states are distinguishable
    double UL[3] = {uLx,uLy,uLz};
    double UR[3] = {uRx,uRy,uRz};
    double drho = fabs(rhoR-rhoL);
    double du = fabs(uLx-uRx)+fabs(uLy-uRy)+fabs(uLz-uRz);
    if ((drho == 0) && (du == 0)) {
      // Identical states
      v_hat = vLx;

    } else if ((drho + du > 1e-15) && (du > 1e-15)) {

      // compute the coefficients of the polynomial
      double coeff[4];
      sr_poly_coeffs(&csr, coeff);

      // compute Vx_hat
      struct gkyl_root_intervals root_intervals; 
      double domain[2] = {-1.1, 1.1};
      double tol = 1e-16;

      // compute root inverals, refine them via bisection, compute the result via riddrrs
      root_intervals = gkyl_calc_quartic_root_intervals( coeff, domain, tol);
      //gkyl_refine_root_intervals_bisection(&root_intervals, tol);
      gkyl_root_isolation_from_intervals_via_ridders(&root_intervals, tol);

      // Select the proper root
      int status; 
      v_hat = root_select(&root_intervals,UL,UR,&status,&csr);

      // Compute U from Vx
      calc_u_from_result(u,v_hat,&csr);
      bool pass_entropy_check = entropy_check(UL, UR, u);
      if (pass_entropy_check == 0 && status == 0){
      } else {
        printf("Entropy test fails with by entropy_test: %d, status: %d (0 is pass, 1-n is fail) \n",pass_entropy_check, status);
        //printf("v_computed = %1.16e; uxhat = %1.16e; uyhat = %1.16e; uzhat = %1.16e;\n",v_hat,u[0],u[1],u[2]);
        //printf("rhol = %1.16e; ulx = %1.16e; uly = %1.16e; ulz = %1.16e;\n",rhoL,uLx,uLy,uLz);
        //printf("rhor = %1.16e; urx = %1.16e; ury = %1.16e; urz = %1.16e;\n",rhoR,uRx,uRy,uRz);
        //printf("Poly Coeff: [");
        for (int i=0; i<4; ++i) //printf("%1.16e,",coeff[i]);
        //printf("1.0]\n");

        // Check v_hat_avg
        //printf("\nReplacing with v_hat_avg\n");
        calc_u_from_result(u,v_hat_avg,&csr);
        pass_entropy_check = entropy_check(UL, UR, u);
        ////printf("v_avg = %1.16e; uxhat = %1.16e; uyhat = %1.16e; uzhat = %1.16e;\n",v_hat_avg,u[0],u[1],u[2]);
        return c*v_hat_avg;
     }

    } else {
      ////printf("\nReplacing with v_hat_avg\n");
      return c*v_hat_avg;
    }

  } else if (ql[0] <= 0.0 && qr[0] > 0.0) {
    v_hat = vRx;
  } else if (qr[0] <= 0.0 && ql[0] > 0.0) {
    v_hat = vLx;
  } else { // qr[0] and ql[0] are zero/negative
    v_hat = 0.0;
  }

  // Return to units of c
  u[0] = c*u[0];
  u[1] = c*u[1];
  u[2] = c*u[2];

  // Renormalize to c and return
  return c*v_hat;
}


// Waves and speeds using Roe averaging
// corrective terms
static double
wave_roe_sr(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *delta, const double *ql, const double *qr, double *waves, double *s)
{
  double f[4];
  printf("TEMP Set c = 1 in wave_roe_sr\n");
  const double c = 1.0; //299792458.0;
  double vl = ql[NUX]/sqrt(ql[0]*ql[0] + (ql[NUX]*ql[NUX] + ql[NUY]*ql[NUY] + ql[NUZ]*ql[NUZ])/(c*c));
  double vr = qr[NUX]/sqrt(qr[0]*qr[0] + (qr[NUX]*qr[NUX] + qr[NUY]*qr[NUY] + qr[NUZ]*qr[NUZ])/(c*c)); 

  double *wv = 0;

  // corrective term
  if ((vl < 0.0) && (0.0 < vr)) { // vacuum intermediate state will be formed
    cold_sr_fluid_flux(ql, f);
    wv = &waves[0];
    for(int m=0; m<4; ++m) wv[m] = -f[m];
    s[0] = vl;

    cold_sr_fluid_flux(qr, f);
    wv = &waves[4];
    for(int m=0; m<4; ++m) wv[m] = f[m];
    s[1] = vr;
  }
  else {
    // no vacuum state
    double rl = ql[0];
    double rr = qr[0];
    double u_diag[3];
    double gamma_avg;

    // compute Roe averaged speed
    //double vel = compute_sr_roe_averaged_velocity(ql,qr,c,u_diag);
    ////printf("vL: %1.16f, vs %1.16f, vR: %1.16f\n",vl,vel,vr);

    // Compute Roe-averaged speed using the cold limit of the conservation laws
    double vel = compute_sr_roe_averaged_velocity_cold_limit(ql,qr,c,u_diag,&gamma_avg);
            
    if(vel<0) {
      wv = &waves[0];
      for(int m=0; m<4; ++m)
        wv[m] = delta[m];

      wv = &waves[4];
      for(int m=0; m<4; ++m)
        wv[m] = 0.0;
    }
    else {
      wv = &waves[0];
      for(int m=0; m<4; ++m)
        wv[m] = 0.0;

      wv = &waves[4];
      for(int m=0; m<4; ++m)
        wv[m] = delta[m];
    }
    s[0] = vel;
    s[1] = vel;
  }

  return fmax(fabs(s[0]), fabs(s[1]));
}

static void
qfluct_roe(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *ql, const double *qr, const double *waves, const double *s,
  double *amdq, double *apdq)
{
  int meqn = 4, mwaves = 2;
  
  for (int m=0; m<meqn; ++m) {
    amdq[m] = 0.0; apdq[m] = 0.0;

    for (int mw=0; mw<mwaves; ++mw) {
      const double *wv = &waves[mw*meqn];
      
      if (s[mw] < 0.0)
        amdq[m] += s[mw]*wv[m];
      else
        apdq[m] += s[mw]*wv[m];
    }
  }
}

// First order solve
static void
ffluct_roe(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *ql, const double *qr, const double *waves, const double *s,
  double *amdq, double *apdq)
{
  int meqn = 4, mwaves = 2;
  
  for (int m=0; m<meqn; ++m) {
    amdq[m] = 0.0; apdq[m] = 0.0;

    for (int mw=0; mw<mwaves; ++mw) {
      const double *wv = &waves[mw*meqn];
      
      if (s[mw] < 0.0) {
        amdq[m] += wv[m];
      }
      else if (s[mw] > 0.0) {
        apdq[m] += wv[m];
      }
      else {
        amdq[m] += 0.5*wv[m];
        apdq[m] += 0.5*wv[m];
      }
    }
  }
}

static double
flux_jump_sr(const struct gkyl_wv_eqn *eqn, const double *ql, const double *qr, double *flux_jump_sr)
{
  double fr[4], fl[4];
  cold_sr_fluid_flux(ql, fl);
  cold_sr_fluid_flux(qr, fr);

  for (int m=0; m<4; ++m) flux_jump_sr[m] = fr[m]-fl[m];

 // Vn = NUn/sqrt(N^2 + NU^2/c^2)
  const double c = 299792458.0;
  double amaxl = ql[NUX]/sqrt(ql[0]*ql[0] + (ql[NUX]*ql[NUX] + ql[NUY]*ql[NUY] + ql[NUZ]*ql[NUZ])/(c*c));
  double amaxr = qr[NUX]/sqrt(qr[0]*qr[0] + (qr[NUX]*qr[NUX] + qr[NUY]*qr[NUY] + qr[NUZ]*qr[NUZ])/(c*c)); 

  return fmax(amaxl, amaxr);
}

static bool
check_inv(const struct gkyl_wv_eqn *eqn, const double *q)
{
  return q[0] > 0.0;
}

static double
max_speed_sr(const struct gkyl_wv_eqn *eqn, const double *q)
{
  const struct wv_cold_sr_fluid *cold_sr_fluid = container_of(eqn, struct wv_cold_sr_fluid, eqn);
  const double c = 299792458.0;
  return fabs(q[NUX]/sqrt(q[0]*q[0] + (q[NUX]*q[NUX] + q[NUY]*q[NUY] + q[NUZ]*q[NUZ])/(c*c)));
}

struct gkyl_wv_eqn*
gkyl_wv_cold_sr_fluid_new(void)
{
  struct wv_cold_sr_fluid *cold_sr_fluid = gkyl_malloc(sizeof(struct wv_cold_sr_fluid));

  cold_sr_fluid->eqn.type = GKYL_EQN_COLDFLUID_SR;
  cold_sr_fluid->eqn.num_equations = 4;
  cold_sr_fluid->eqn.num_waves = 2;
  cold_sr_fluid->eqn.num_diag = 5; // KE is final component
  
  cold_sr_fluid->eqn.waves_func = wave_roe_sr;
  cold_sr_fluid->eqn.qfluct_func = qfluct_roe;
  cold_sr_fluid->eqn.ffluct_func = ffluct_roe;
  cold_sr_fluid->eqn.flux_jump = flux_jump_sr;
  
  cold_sr_fluid->eqn.check_inv_func = check_inv;
  cold_sr_fluid->eqn.max_speed_func = max_speed_sr;

  cold_sr_fluid->eqn.rotate_to_local_func = rot_to_local;
  cold_sr_fluid->eqn.rotate_to_global_func = rot_to_global;

  cold_sr_fluid->eqn.cons_to_riem = cons_to_riem_sr;
  cold_sr_fluid->eqn.riem_to_cons = riem_to_cons_sr;

  cold_sr_fluid->eqn.cons_to_diag = cold_sr_fluid_cons_to_diag;

  cold_sr_fluid->eqn.ref_count = gkyl_ref_count_init(cold_sr_fluid_free);

  return &cold_sr_fluid->eqn;
}
