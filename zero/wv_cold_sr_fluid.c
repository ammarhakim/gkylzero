#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_wv_cold_sr_fluid.h>

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

static inline void
cold_sr_fluid_flux(const double q[4], double *flux)
{
  // Vx = NUx/sqrt(N^2 + NU^2/c^2) 
  ////printf("TEMP SPEED IN COLD_SR_FLUX\n");
  //const double c = 1.0; //c = 299792458.0; 
  const double c = 299792458.0; 
  //if ((q[0]*q[0] + (q[NUX]*q[NUX] + q[NUY]*q[NUY] + q[NUZ]*q[NUZ])/(c*c)) < 0) //printf("invalid sqrt");
  //for (int i=0; i<4; ++i) if isnan(q[i]) //printf("q[%d] is nan\n",i);
  double Vx = q[NUX]/sqrt(q[0]*q[0] + (q[NUX]*q[NUX] + q[NUY]*q[NUY] + q[NUZ]*q[NUZ])/(c*c)); 
  flux[0] =  q[0]*Vx; // N*Vx
  flux[NUX] =  q[NUX]*Vx; // N*Ux*Vx
  flux[NUY] =  q[NUY]*Vx; // N*Uy*Vx
  flux[NUZ] =  q[NUZ]*Vx; // N*Uz*Vx
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



// Waves and speeds using Roe averaging
// corrective terms
static double
wave_roe_sr(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *delta, const double *ql, const double *qr, double *waves, double *s)
{
  ////printf("TEMP Set c = 1 in wave_roe_sr\n");
  //const double c = 1.0;
  const double c = 299792458.0;

  double *wv = 0;


  // As long as one density is positive on both sides:
  if (ql[0] > 0.0 || qr[0] > 0.0){

    // isolate left and right states
    double rhor = qr[0];
    if (qr[0] < 0.0) {
      rhor = 0.0;
      //printf("rhor fix!\n");
    }
    double urx = qr[1]/(rhor);
    double ury = qr[2]/(rhor);
    double urz = qr[3]/(rhor);
    if (qr[0] < 0.0) {
      //urx = 0.0;
      //ury = 0.0;
      //urz = 0.0;
    }
    double rhol = ql[0];
    // Fix density if negative:
    if (ql[0] < 0.0) {
      rhol = 0.0;
      //printf("rhol fix!\n");
    }
    double ulx = ql[1]/(rhol);
    double uly = ql[2]/(rhol);
    double ulz = ql[3]/(rhol);
    if (ql[0] < 0.0) {
      //ulx = 0.0;
      //uly = 0.0;
      //ulz = 0.0;
    }



    // compute the constants:
    double gammal = sqrt(1.0 + (ulx*ulx + uly*uly + ulz*ulz)/(c*c));
    double gammar = sqrt(1.0 + (urx*urx + ury*ury + urz*urz)/(c*c));
    double vlx = ulx/gammal;
    double vrx = urx/gammar;
    double vly = uly/gammal;
    double vry = ury/gammar;
    double vlz = ulz/gammal;
    double vrz = urz/gammar;

    // Primative rho
    double rhol_prim = rhol/gammal;
    double rhor_prim = rhor/gammar;

    // TEMP: Norm of density prim (not working!)
    //double max_prim = fmin(rhol_prim,rhor_prim);
    //rhol_prim = rhol_prim/max_prim;
    //rhor_prim = rhor_prim/max_prim;

    // Compute the primative-parameterization state vector w
    // these are the averages of the left and right states
    //double k = (sqrt(rhol_prim) + sqrt(rhor_prim))/(c);
    //double w0 = sqrt(rhol_prim)*gammal + sqrt(rhor_prim)*gammar;
    //double w1 = sqrt(rhol_prim)*gammal*vlx/c + sqrt(rhor_prim)*gammar*vrx/c;
    //double w2 = sqrt(rhol_prim)*gammal*vly/c + sqrt(rhor_prim)*gammar*vry/c;
    //double w3 = sqrt(rhol_prim)*gammal*vlz/c + sqrt(rhor_prim)*gammar*vrz/c;

    // Can /1.0e12 each term as a sort of normalization, but doesn't change the answer

    double k = (sqrt(rhol_prim) + sqrt(rhor_prim))/(2.0*c);
    double w0 = (sqrt(rhol_prim)*gammal + sqrt(rhor_prim)*gammar)/(2.0);
    double w1 = (sqrt(rhol_prim)*gammal*vlx/c + sqrt(rhor_prim)*gammar*vrx/c)/(2.0);
    double w2 = (sqrt(rhol_prim)*gammal*vly/c + sqrt(rhor_prim)*gammar*vry/c)/(2.0);
    double w3 = (sqrt(rhol_prim)*gammal*vlz/c + sqrt(rhor_prim)*gammar*vrz/c)/(2.0);

    // Assign the jump in the state vector, d = (d0,d1,d2,d3) 
    double d0 = qr[0] - ql[0]; 
    double d1 = qr[1] - ql[1]; 
    double d2 = qr[2] - ql[2]; 
    double d3 = qr[3] - ql[3]; 

      

    //if ( (ql[0] != qr[0])  ){ // 
    if  ((w1 != 0) || (w2 != 0) || (w3 != 0)){
      //printf("\n*****\n");
      //printf("(STATE JUMP): k: %1.16e, w[0]: %1.16e, w[1]: %1.16e, w[2]: %1.16e, w[3]: %1.16e\n",k,w0,w1,w2,w3);
      //printf("DENOM1: %1.16e, DENOM2: %1.16e\n",(w0 * (c * c * k * k - w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3)), ((c * c * k * k - w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3) * (c * c * k * k + w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3)));
      //printf("NUM Wv0[0]: %1.16e, NUM Wv0[1]: %1.16e\n",(c * (d1 * c * c * k * k * k - d0 * c * c * k * k * w1 + d1 * k * w0 * w0 - d1 * k * w1 * w1 - 2 * d2 * k * w1 * w2 - 2 * d3 * k * w1 * w3 + d1 * k * w2 * w2 + d1 * k * w3 * w3 - d0 * w0 * w0 * w1 + d0 * w1 * w1 * w1 + d0 * w1 * w2 * w2 + d0 * w1 * w3 * w3)),(2 * c * w1 * (d1 * c * c * k * k - d0 * w1 * c * c * k + d1 * w2 * w2 - d2 * w1 * w2 + d1 * w3 * w3 - d3 * w1 * w3)));
      //printf("NUM Wv1[0]: %1.16e, NUM Wv1[1]: %1.16e\n",-(2 * c *  k * w0 * (d1 * c * c * k * k - 2 * d0 * c * c * k * w1 + d1 * w0 * w0 - d1 * w1 * w1 - 2 * d2 * w1 * w2 - 2 * d3 * w1 * w3 + d1 * w2 * w2 + d1 * w3 * w3)),-(2 * c * w0 * w1 * (d1 * c * c * k * k - 2 * d0 * c * c * k * w1 + d1 * w0 * w0 - d1 * w1 * w1 - 2 * d2 * w1 * w2 - 2 * d3 * w1 * w3 + d1 * w2 * w2 + d1 * w3 * w3)));
      //printf("Wv0[0]: %1.16e, Wv0[1]: %1.16e\n",(c * (d1 * c * c * k * k * k - d0 * c * c * k * k * w1 + d1 * k * w0 * w0 - d1 * k * w1 * w1 - 2 * d2 * k * w1 * w2 - 2 * d3 * k * w1 * w3 + d1 * k * w2 * w2 + d1 * k * w3 * w3 - d0 * w0 * w0 * w1 + d0 * w1 * w1 * w1 + d0 * w1 * w2 * w2 + d0 * w1 * w3 * w3)) / (w0 * (c * c * k * k - w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3)),(2 * c * w1 * (d1 * c * c * k * k - d0 * w1 * c * c * k + d1 * w2 * w2 - d2 * w1 * w2 + d1 * w3 * w3 - d3 * w1 * w3)) / (w0 * (c * c * k * k - w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3)));
      //printf("Wv1[0]: %1.16e, Wv1[1]: %1.16e\n",-(2 * c *  k * w0 * (d1 * c * c * k * k - 2 * d0 * c * c * k * w1 + d1 * w0 * w0 - d1 * w1 * w1 - 2 * d2 * w1 * w2 - 2 * d3 * w1 * w3 + d1 * w2 * w2 + d1 * w3 * w3)) / ((c * c * k * k - w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3) * (c * c * k * k + w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3)),-(2 * c * w0 * w1 * (d1 * c * c * k * k - 2 * d0 * c * c * k * w1 + d1 * w0 * w0 - d1 * w1 * w1 - 2 * d2 * w1 * w2 - 2 * d3 * w1 * w3 + d1 * w2 * w2 + d1 * w3 * w3)) / ((c * c * k * k - w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3) * (c * c * k * k + w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3)));
      //printf("(c * c * k * k): %1.16e, - (w0 * w0): %1.16e, w1 * w1: %1.16e, (w2 * w2): %1.16e, (w3 * w3): %1.16e\n",(c * c * k * k), - (w0 * w0),w1 * w1, w2 * w2, w3 * w3);
      //printf("ql[0]: %1.16e, ql[1]: %1.16e, ql[2]: %1.16e, qr[0]: %1.16e, qr[1]: %1.16e, qr[2]: %1.16e\n",ql[0],ql[1],ql[2],qr[0],qr[1],qr[2]);
      //printf("*****\n");
    }

    // Wave 1: eigenvalue is w1/w0 repeated, three waves are lumped into one
    // waves = Vx*[N, NUx, NUy, Nuz]
    wv = &waves[0];

    // Diagnostics
    //double wv_val[4] = { 0.0, 0.0, 0.0, 0.0}; 
    //double wv_num[4] = { 0.0, 0.0, 0.0, 0.0}; 
    //double wv_denom;
    //wv_denom = (w0 * (c * c * k * k - w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3));
    //wv_num[0] = (c * (d1 * c * c * k * k * k - d0 * c * c * k * k * w1 + d1 * k * w0 * w0 - d1 * k * w1 * w1 - 2 * d2 * k * w1 * w2 - 2 * d3 * k * w1 * w3 + d1 * k * w2 * w2 + d1 * k * w3 * w3 - d0 * w0 * w0 * w1 + d0 * w1 * w1 * w1 + d0 * w1 * w2 * w2 + d0 * w1 * w3 * w3)) ;
    //wv_num[1] = (2 * c * w1 * (d1 * c * c * k * k - d0 * w1 * c * c * k + d1 * w2 * w2 - d2 * w1 * w2 + d1 * w3 * w3 - d3 * w1 * w3)) / (w0 * (c * c * k * k - w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3));
    //wv_num[2] = (c * (d2 * c * c * k * k * w1 + d1 * c * c * k * k * w2 - 2 * d0 * c * c * k * w1 * w2 - d2 * w0 * w0 * w1 + d1 * w0 * w0 * w2 + d2 * w1 * w1 * w1 - d1 * w1 * w1 * w2 - d2 * w1 * w2 * w2 - 2 * d3 * w1 * w2 * w3 + d2 * w1 * w3 * w3 + d1 * w2 * w2 * w2 + d1 * w2 * w3 * w3));
    //wv_num[3] = (c * (d3 * c * c * k * k * w1 + d1 * c * c * k * k * w3 - 2 * d0 * c * c * k * w1 * w3 - d3 * w0 * w0 * w1 + d1 * w0 * w0 * w3 + d3 * w1 * w1 * w1 - d1 * w1 * w1 * w3 + d3 * w1 * w2 * w2 - 2 * d2 * w1 * w2 * w3 - d3 * w1 * w3 * w3 + d1 * w2 * w2 * w3 + d1 * w3 * w3 * w3));
    //wv_val[0] = wv_num[0]/wv_denom;
    //wv_val[1] = wv_num[1]/wv_denom;
    //wv_val[2] = wv_num[2]/wv_denom;
    //wv_val[3] = wv_num[3]/wv_denom;


    wv[0] =  (c * (d1 * c * c * k * k * k - d0 * c * c * k * k * w1 + d1 * k * w0 * w0 - d1 * k * w1 * w1 - 2 * d2 * k * w1 * w2 - 2 * d3 * k * w1 * w3 + d1 * k * w2 * w2 + d1 * k * w3 * w3 - d0 * w0 * w0 * w1 + d0 * w1 * w1 * w1 + d0 * w1 * w2 * w2 + d0 * w1 * w3 * w3)) / (w0 * (c * c * k * k - w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3)); 
    wv[1] =  (2 * c * w1 * (d1 * c * c * k * k - d0 * w1 * c * c * k + d1 * w2 * w2 - d2 * w1 * w2 + d1 * w3 * w3 - d3 * w1 * w3)) / (w0 * (c * c * k * k - w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3));
    wv[2] =  (c * (d2 * c * c * k * k * w1 + d1 * c * c * k * k * w2 - 2 * d0 * c * c * k * w1 * w2 - d2 * w0 * w0 * w1 + d1 * w0 * w0 * w2 + d2 * w1 * w1 * w1 - d1 * w1 * w1 * w2 - d2 * w1 * w2 * w2 - 2 * d3 * w1 * w2 * w3 + d2 * w1 * w3 * w3 + d1 * w2 * w2 * w2 + d1 * w2 * w3 * w3)) / (w0 * (c * c * k * k - w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3));
    wv[3] =  (c * (d3 * c * c * k * k * w1 + d1 * c * c * k * k * w3 - 2 * d0 * c * c * k * w1 * w3 - d3 * w0 * w0 * w1 + d1 * w0 * w0 * w3 + d3 * w1 * w1 * w1 - d1 * w1 * w1 * w3 + d3 * w1 * w2 * w2 - 2 * d2 * w1 * w2 * w3 - d3 * w1 * w3 * w3 + d1 * w2 * w2 * w3 + d1 * w3 * w3 * w3)) / (w0 * (c * c * k * k - w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3));
    for (int i=0; i<4; ++i) if (isnan(wv[i])) {
      //if ( (ql[0] != qr[0])  )
        //printf("wv[%d]: %1.16e, wv_num/denom: %1.16e/%1.16e\n",i,wv[i],wv_num[i],wv_denom);
      wv[i] = 0;
      //if ( (ql[0] != qr[0])  )
        //printf("S1 RESET WV[%d] is nan\n",i);
      if ((w1 != 0) || (w2 != 0) || (w3 != 0)){
        //printf("S1 RESET WV[%d] is nan\n",i);
        //printf("ql[0]: %1.16e, ql[1]: %1.16e, ql[2]: %1.16e, qr[0]: %1.16e, qr[1]: %1.16e, qr[2]: %1.16e\n",ql[0],ql[1],ql[2],qr[0],qr[1],qr[2]);
      }
    }
    for (int i=0; i<4; ++i) if (!(isfinite(wv[i]))) {
      //if ( (ql[0] != qr[0])  )
        //printf("wv[%d]: %1.16e, wv_num/denom: %1.16e/%1.16e\n",i,wv[i],wv_num[i],wv_denom);
      wv[i] = 0;
      //if ( (ql[0] != qr[0])  )
        //printf("S1 RESET WV[%d] is inf\n",i);
      if ((w1 != 0) || (w2 != 0) || (w3 != 0)){
        //printf("S1 RESET WV[%d] is inf\n",i);
        //printf("ql[0]: %1.16e, ql[1]: %1.16e, ql[2]: %1.16e, qr[0]: %1.16e, qr[1]: %1.16e, qr[2]: %1.16e\n",ql[0],ql[1],ql[2],qr[0],qr[1],qr[2]);
      }
    }
    s[0] = (c*w1)/w0;
    if (isnan(s[0])) {
      s[0] = 0;
      //printf("RESET s[0] is nan\n");
    }
    if (!(isfinite(s[0]))) {
      s[0] = 0;
      //printf("RESET s[0] is inf\n");
    }

    if ((w1 != 0) || (w2 != 0) || (w3 != 0)){
    //if ( (ql[0] != qr[0])  ){
      //printf("(SET 1): wv[0]: %1.16e, wv[1]: %1.16e, wv[2]: %1.16e, wv[3]: %1.16e, s[0]: %1.16e\n",wv[0],wv[1],wv[2],wv[3],s[0]);
      //printf("(SET 1): vlx: %1.16e, s1: %1.16e, vrx: %1.16e\n",vlx,s[0],vrx);
    }

    // Wave 2: eigenvalue is 2*w0*w1/( k*k + w0*w0 + w1*w1 + w2*w2 + w3*w3 );
    wv = &waves[4];
    wv[0] =  -(2 * c *  k * w0 * (d1 * c * c * k * k - 2 * d0 * c * c * k * w1 + d1 * w0 * w0 - d1 * w1 * w1 - 2 * d2 * w1 * w2 - 2 * d3 * w1 * w3 + d1 * w2 * w2 + d1 * w3 * w3)) / ((c * c * k * k - w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3) * (c * c * k * k + w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3));
    wv[1] =  -(2 * c * w0 * w1 * (d1 * c * c * k * k - 2 * d0 * c * c * k * w1 + d1 * w0 * w0 - d1 * w1 * w1 - 2 * d2 * w1 * w2 - 2 * d3 * w1 * w3 + d1 * w2 * w2 + d1 * w3 * w3)) / ((c * c * k * k - w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3) * (c * c * k * k + w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3));
    wv[2] =  -(2 * c * w0 * w2 * (d1 * c * c * k * k - 2 * d0 * c * c * k * w1 + d1 * w0 * w0 - d1 * w1 * w1 - 2 * d2 * w1 * w2 - 2 * d3 * w1 * w3 + d1 * w2 * w2 + d1 * w3 * w3)) / ((c * c * k * k - w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3) * (c * c * k * k + w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3));
    wv[3] =  -(2 * c * w0 * w3 * (d1 * c * c * k * k - 2 * d0 * c * c * k * w1 + d1 * w0 * w0 - d1 * w1 * w1 - 2 * d2 * w1 * w2 - 2 * d3 * w1 * w3 + d1 * w2 * w2 + d1 * w3 * w3)) / ((c * c * k * k - w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3) * (c * c * k * k + w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3));
    for (int i=0; i<4; ++i) if (isnan(wv[i])){
      wv[i] = 0;
      //if ( (ql[0] != qr[0])  )
        //printf("S2 RESET WV[%d] is nan\n",i);
      if ((w1 != 0) || (w2 != 0) || (w3 != 0)){
        //printf("S2 RESET WV[%d] is nan\n",i);
        //printf("ql[0]: %1.16e, ql[1]: %1.16e, ql[2]: %1.16e, qr[0]: %1.16e, qr[1]: %1.16e, qr[2]: %1.16e\n",ql[0],ql[1],ql[2],qr[0],qr[1],qr[2]); 
      }
    }
    for (int i=0; i<4; ++i) if (!(isfinite(wv[i]))) {
      wv[i] = 0;
      //if ( (ql[0] != qr[0])  )
        //printf("S2 RESET WV[%d] is inf\n",i);
      if ((w1 != 0) || (w2 != 0) || (w3 != 0)){
        //printf("S2 RESET WV[%d] is inf\n",i);
        //printf("ql[0]: %1.16e, ql[1]: %1.16e, ql[2]: %1.16e, qr[0]: %1.16e, qr[1]: %1.16e, qr[2]: %1.16e\n",ql[0],ql[1],ql[2],qr[0],qr[1],qr[2]);
      }
    }
    s[1] =  2.0*c*w0*w1/( c*c*k*k + w0*w0 + w1*w1 + w2*w2 + w3*w3 ); // (c*w1)/w0;
    if (isnan(s[1])) {
      s[1] = 0;
      //printf("RESET s[1] is nan\n");
    }
    if (!(isfinite(s[1]))) {
      s[1] = 0;
      //printf("RESET s[1] is inf\n");
    }




    //double other_v3 = u[0]/sqrt(1.0 + u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);;
    //if (w1 != 0){
      //printf("Eigenvalues, Vx = w1/w0: %1.16e, 2*w0*w1/...: %1.16e\n",s[0], s[1]);
      //printf("vxl: %1.16e, vxr: %1.16e\n\n",vlx,vrx);
    //}


    if ((w1 != 0) || (w2 != 0) || (w3 != 0)){
    //if ( (ql[0] != qr[0])  ){
      //wv = &waves[0];
      //double fr[4], fl[4];
      //cold_sr_fluid_flux(ql, fl);
      //cold_sr_fluid_flux(qr, fr);
      //printf("(SET 2): wv[0]: %1.16e, wv[1]: %1.16e, wv[2]: %1.16e, wv[3]: %1.16e, s[1]: %1.16e\n",wv[0],wv[1],wv[2],wv[3],s[1]);
      //printf("(SET 2): vlx: %1.16e, s2: %1.16e, vrx: %1.16e\n",vlx,s[1],vrx);
      //printf("ql[0]: %1.16e, ql[1]: %1.16e, ql[2]: %1.16e, qr[0]: %1.16e, qr[1]: %1.16e, qr[2]: %1.16e\n",ql[0],ql[1],ql[2],qr[0],qr[1],qr[2]); 
      //for (int i=0; i<4; ++i) printf("(FLUX) FL[%d] %1.16e, F_tot[%d] %1.16e, FR[%d] %1.16e, || F0[%d] %1.16e, F1[%d] %1.16e,\n",i,fl[i],i,wv[i]+ wv[4+i],i,fr[i],i,wv[i],i,wv[4+i]);
    }


  } else {
    wv = &waves[0];
    wv[0] = 0.0;
    wv[1] = 0.0;
    wv[2] = 0.0;
    wv[3] = 0.0;
    s[0] = 0.0;
    wv = &waves[4];
    wv[0] = 0.0;
    wv[1] = 0.0;
    wv[2] = 0.0;
    wv[3] = 0.0;
    s[1] = 0.0;
  }

  // Print if the waves are nan
  //for (int i=0; i<8; ++i) if (isnan(waves[i])) //printf("waves[%d] is nan\n",i);
  //for (int i=0; i<2; ++i) if (isnan(s[i])) //printf("s[%d] is nan\n",i);


  return fmax(fabs(s[0]), fabs(s[1]));
}

static void
qfluct_roe(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *ql, const double *qr, const double *waves, const double *s,
  double *amdq, double *apdq)
{
  //printf("Q-Waves will not work with this system (L,R Eigenvectors are not unique)\n");
}

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

  //for (int i=0; i<2; ++i) if (isnan(s[i])) //printf("s[%d] is nan ffluc_roe\n",i);
  //for (int i=0; i<8; ++i) if (isnan(waves[i])) //printf("wv[%d] is nan ffluc_roe\n",i);
  //for (int i=0; i<4; ++i) if (isnan(amdq[i])) //printf("amdq[%d] is nan ffluc_roe\n",i);
  //for (int i=0; i<4; ++i) if (isnan(amdq[i])) //printf("waves1[%d] is: %1.16e, waves2[%d] is: %1.16e, \n",i,waves[i],i+4,waves[i*meqn]);
  //for (int i=0; i<4; ++i) if (isnan(apdq[i])) //printf("apdq[%d] is nan ffluc_roe\n",i);
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
  //if ((ql[0]*ql[0] + (ql[NUX]*ql[NUX] + ql[NUY]*ql[NUY] + ql[NUZ]*ql[NUZ])/(c*c)) < 0) //printf("invalid sqrt");
  //if ((qr[0]*qr[0] + (qr[NUX]*qr[NUX] + qr[NUY]*qr[NUY] + qr[NUZ]*qr[NUZ])/(c*c)) < 0) //printf("invalid sqrt");
  double amaxl =  ql[NUX]/sqrt(ql[0]*ql[0] + (ql[NUX]*ql[NUX] + ql[NUY]*ql[NUY] + ql[NUZ]*ql[NUZ])/(c*c));
  double amaxr =  qr[NUX]/sqrt(qr[0]*qr[0] + (qr[NUX]*qr[NUX] + qr[NUY]*qr[NUY] + qr[NUZ]*qr[NUZ])/(c*c)); 

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
  //if ((q[0]*q[0] + (q[NUX]*q[NUX] + q[NUY]*q[NUY] + q[NUZ]*q[NUZ])/(c*c)) < 0) //printf("invalid sqrt");
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
