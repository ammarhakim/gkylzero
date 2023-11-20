#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_wv_cold_sr_fluid.h>

// Files for roe-averged velocity 
#include <gkyl_math.h>

#define NUX 1
#define NUY 2
#define NUZ 3

#define SQ(x) ((x)*(x))

struct wv_cold_sr_fluid {
  struct gkyl_wv_eqn eqn; // base object
};

// relativistic case
struct cold_sr {
  double rhol, rhor;
  double vl[3], vr[3];
};

void
cold_sr_fluid_flux(const double q[5], double flux[5])
{
  // Vx = NUx/(N^2 + NU^2/c^2) 
  const double c = 299792458.0;
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
calc_ufromv(const double v[3], double u[3])
{
  double v2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
  double gamma1 = 1/sqrt(1-v2);
  u[0] = v[0]*gamma1;
  u[1] = v[1]*gamma1;
  u[2] = v[2]*gamma1;
}

// returns 1 if it passes, 0 if fails
static bool 
entropy_check(const double UL[3], const double UR[4], double u[3], const double c)
{
  bool cond = 1;
  double tol = 1e-15;
  for(int i=0; i<3; ++i){
    cond = cond && (fmin(UR[i],UL[i]) - tol) <= u[i] && u[i] <= (fmax(UR[i],UL[i]) + tol);
  }
  return cond;
}


static void
calc_u_from_result(double uhat[3], double v_hat, void *ctx)
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

static void
sr_poly_coeffs(void *ctx, double coeff[4])
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
  double p3coeff = 2*a*drho+2*A[2]*B[2]+2*A[1]*B[1]+2*A[0]*B[0]; 
  double p2coeff = SQ(a)+SQ(A[2])+SQ(A[1])-SQ(B[0])+SQ(A[0]); 
  double p1coeff = -2*A[0]*B[0]; 
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
root_select(struct gkyl_lo_poly_roots *roots, double vrx, double vlx, int status)
{

  // status 
  status = 1;

  // Check that each root has a small complex component
  double real_comps[4], imag_comps[4], tol[4];
  for (int i=0; i<4; ++i) {
    real_comps[i] = roots->rpart[i];
    imag_comps[i] = roots->impart[i];
    tol[i] = roots->err[i];
  }

  // Find the closest real root
  double v_mean = 0.5*(vrx + vlx);;
  int real_root_only[4] = { 0, 0, 0, 0 };
  double v_diff_real[4];
  int i_choosen_root = 5;
  double min_val = 1.0;
  for (int i=0; i<4; ++i) {
    if (fabs(2.0*tol[i]) > fabs(imag_comps[i])) {

      // Then the root is considered purely real
      real_root_only[i] = 1;

      // compute the difference of the real part from mean vl
      v_diff_real[i] = fabs( real_comps[i] - v_mean );

      if(min_val > v_diff_real[i]) {
        i_choosen_root = i;
        min_val = v_diff_real[i];
      }
    }
  }

  // If choosen i_choosen_root = 5 then throw an error:
  if (i_choosen_root == 5){
    printf("THERE WAS AN ERROR ISOLATING THE REAL ROOT\n");
    status = 0;
  }

  // Return the root which best matches
  return real_comps[i_choosen_root];

}


static double
compute_sr_roe_averaged_velocity(const double ql[4], const double qr[4], const double c, double u[3]) 
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

    // compute the coefficients of the polynomial
    double coeff[4];
    sr_poly_coeffs(&csr, coeff);

    // compute Vx_hat
    enum gkyl_lo_poly_order order = GKYL_LO_POLY_4;
    struct gkyl_lo_poly_roots roots = gkyl_calc_lo_poly_roots(order, coeff);

    // Select the proper root
    int status = 0;
    double v_hat = root_select(&roots,vRx,vLx,status);

    // Compute U from Vx
    double u[3];
    calc_u_from_result(u,v_hat,&csr);

    // Compute a default answer, the v-hat-avg
    double u_avg[3] = {(uLx+uRx)/2.0,(uLy+uRy)/2.0,(uLz+uRz)/2.0};
    double v_hat_avg = u_avg[0]/sqrt(1.0 + u_avg[0]*u_avg[0] + u_avg[1]*u_avg[1] + u_avg[2]*u_avg[2] );

    // Do the entropy check
    double UL[3] = {uLx,uLy,uLz};
    double UR[3] = {uRx,uRy,uRz};
    double drho = fabs(rhoR-rhoL);
    double du = fabs(uLx-uRx)+fabs(uLy-uRy)+fabs(uLz-uRz);
    if ((drho + du > 1e-13) && (du > 1e-15)) {
      bool pass_entropy_check = entropy_check(UL, UR, u, c);
      if (pass_entropy_check && status){
      } else {
        printf("Entropy test fails with by entropy_test: %d, status: %d (1 is pass, 0 is fail) \n",pass_entropy_check, status);
        printf("v_computed = %1.16e; uxhat = %1.16e; uyhat = %1.16e; uzhat = %1.16e;\n",v_hat,u[0],u[1],u[2]);
        printf("rhol = %1.16e; ulx = %1.16e; uly = %1.16e; ulz = %1.16e;\n",rhoL,uLx,uLy,uLz);
        printf("rhor = %1.16e; urx = %1.16e; ury = %1.16e; urz = %1.16e;\n",rhoR,uRx,uRy,uRz);
        return c*v_hat_avg;
      }
    } else {
      return c*v_hat_avg;
    }

  } else if (ql[0] <= 0.0 && qr[0] > 0.0) {
    v_hat = vRx;
  } else if (qr[0] <= 0.0 && ql[0] > 0.0) {
    v_hat = vLx;
  } else { // qr[0] and ql[0] are zero/negative
    v_hat = 0.0;
  }

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
  const double c = 299792458.0;
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
    // compute Roe averaged speed
    double vel = compute_sr_roe_averaged_velocity(ql,qr,c,u_diag);
    //printf("vL: %1.16f, vs %1.16f, vR: %1.16f\n",vl,vel,vr);
            
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

 // Vn = NUn/(N^2 + NU^2/c^2)
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
