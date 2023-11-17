#pragma once

#include <complex.h>
#include <gkyl_wv_eqn.h>
#include <gkyl_math.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


static double
sqr(double x)
{
  return x*x;
}


/**
 * Create a new cold-relativistic-fluid equation object.
 * 
 * @return Pointer to cold-relativistic-fluid equation object.
 */
struct gkyl_wv_eqn* gkyl_wv_cold_sr_fluid_new(void);


// relativistic case
struct cold_sr {
  double rhol, rhor;
  double vl[3], vr[3];
};

static inline void
calc_ufromv(const double v[3], double u[3])
{
  double v2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
  double gamma1 = 1/sqrt(1-v2);
  u[0] = v[0]*gamma1;
  u[1] = v[1]*gamma1;
  u[2] = v[2]*gamma1;
}

double tolerance_dk_4th_order(double complex p, double complex q, double complex r, double complex s,
  double complex p_old, double complex q_old, double complex r_old, double complex s_old) {

  // Calculate the absolute differences for each element
  double diff_p = cabs(p - p_old);
  double diff_q = cabs(q - q_old);
  double diff_r = cabs(r - r_old);
  double diff_s = cabs(s - s_old);

  // Sum up the absolute differences
  double sum_of_differences = diff_p + diff_q + diff_r + diff_s;

  return sum_of_differences;
}

static double complex inline
f_poly(double complex x, double a, double b, double c, double d)
{
  return cpow(x,4)+a*cpow(x,3)+b*cpow(x,2)+c*x+d;
}

void
durand_kerner_method(void *ctx, double complex roots[4])
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

  double p4coeff = sqr(drho)+sqr(B[2])+sqr(B[1])+sqr(B[0]); 
  double p3coeff = 2*a*drho+2*A[2]*B[2]+2*A[1]*B[1]+2*A[0]*B[0]; 
  double p2coeff = sqr(a)+sqr(A[2])+sqr(A[1])-sqr(B[0])+sqr(A[0]); 
  double p1coeff = -2*A[0]*B[0]; 
  double p0coeff = -sqr(A[0]); 

  // Find the four roots to the polynomial Vx equation:
  // f(x)=x^{4}+ax^{3}+bx^{2}+cx+d,
  // via durand-kerner: https://en.wikipedia.org/wiki/Durand–Kerner_method
  double a_poly = p3coeff/p4coeff;
  double b_poly = p2coeff/p4coeff;
  double c_poly = p1coeff/p4coeff;
  double d_poly = p0coeff/p4coeff;

  // Compute the intial conditions for solutions
  double complex p, q, r, s;
  double complex p_old = 1.0;
  double complex q_old = I;
  double complex r_old = -1.0;
  double complex s_old = -I;

  // iterate until the desired precision is reached
  int iter = 0;
  double tol = 1.0;
  while(tol > 1e-10 && iter < 10000){

    // Update n-1 -> n for (p, q, r, s)
    p = p_old - f_poly(p_old,a_poly,b_poly,c_poly,d_poly)/( (p_old - q_old)*(p_old - r_old)*(p_old - s_old) );
    q = q_old - f_poly(q_old,a_poly,b_poly,c_poly,d_poly)/( (q_old - p)*(q_old - r_old)*(q_old - s_old) );
    r = r_old - f_poly(r_old,a_poly,b_poly,c_poly,d_poly)/( (r_old - p)*(r_old - q)*(r_old - s_old) );
    s = s_old - f_poly(s_old,a_poly,b_poly,c_poly,d_poly)/( (s_old - p)*(s_old - q)*(s_old - r) );

    // Compute tolerance 
    tol = tolerance_dk_4th_order(p,q,r,s,p_old,q_old,r_old,s_old);

    // Update using old values
    p_old = p;
    q_old = q;
    r_old = r;
    s_old = s;

    // update the iterator
    iter = iter + 1;
  }

  // Save the roots to the outputs
  roots[0] = p;
  roots[1] = q;
  roots[2] = r;
  roots[3] = s;

}

// returns 1 if it passes, 0 if fails
static bool 
entropy_check(const double UL[3], const double UR[4], double u[3], const double c)
{
  bool cond = 1;
  double tol = 1e-15;
  for(int i=0; i<3; ++i){
      //if (fabs(ql[0]) > 1.e-15 && fabs(qr[0]) > 1.e-15){
        cond = cond && (fmin(UR[i],UL[i]) - 1e-15) <= u[i] && u[i] <= (fmax(UR[i],UL[i]) + 1e-15);
      //}
      //cond = cond && (fabs(NU_cons_func( ql, qr, i, u, c)) <= 1.e-10);
  }
  return cond;
}


static double
cold_sr_func(double vhat, void *ctx)
{
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

  double denom = (drho*vhat+a);
  //printf("drho = %lg, vhat = %lg, a = %lg, denom = %lg\n", drho, vhat, a, denom);
  double uhat[3], uhat2 = 0.0;
  for (int i=0; i<3; ++i) {
    uhat[i] = -(B[i]*vhat+A[i])/(drho*vhat+a);
    uhat2 += uhat[i]*uhat[i];
  }
  //printf("vhat=%lg uhat0=%lg uhat2=%lg\n", vhat, uhat[0], uhat2);

  return vhat - uhat[0]/sqrt(1+uhat2);
}

static double
cold_sr_vhat(int status, struct cold_sr *csr)
{
  double rhol = csr->rhol, rhor = csr->rhor;
  double vl[3], vr[3];

  for (int i=0; i<3; ++i) {
    vl[i] = csr->vl[i];
    vr[i] = csr->vr[i];
  }

  double drho = fabs(rhor-rhol);
  double dv = fabs(vl[0]-vr[0])+fabs(vl[1]-vr[1])+fabs(vl[2]-vr[2]);
  if (drho + dv < 1e-14)
    return 0.5*(vl[0]+vr[0]);

  // Originally between +/- c:
  //double fl = cold_sr_func(-1.0+1e-10, csr);
  //double fr = cold_sr_func(1.0-1e-10, csr);
  //struct gkyl_qr_res res = gkyl_ridders(cold_sr_func, csr, -1.0, 1.0, fl, fr, 100, 1e-10);
  double fl = cold_sr_func(vl[0], csr);
  double fr = cold_sr_func(vr[0], csr);
  struct gkyl_qr_res res = gkyl_ridders(cold_sr_func, csr, vl[0], vr[0], fl, fr, 100, 1e-16);

  if(res.status){
    if (res.res > 0.0){
      double fl = cold_sr_func(0.0, csr);
      double fr = cold_sr_func(0.99, csr);
      struct gkyl_qr_res res = gkyl_ridders(cold_sr_func, csr, 0.0, 0.99, fl, fr, 100, 1e-16);
    } else if (res.res < 0.0) {
      double fl = cold_sr_func(0.0, csr);
      double fr = cold_sr_func(-0.99, csr);
      struct gkyl_qr_res res = gkyl_ridders(cold_sr_func, csr, 0.0, -0.99, fl, fr, 100, 1e-16);
    }
  }
  status = res.status;

  return res.res;
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

/**
 * Compute the roe-averaged velocity for the sr-cold fluid equations, solved via
 * Ridders' method.
 *
 * @param ql Conserved variables (left)
 * @param qr Conserved variables (right)
 * @param c speed of light
 * @param u roe averaged momentum (output)
 */
double
compute_sr_roe_averaged_velocity_via_ridders(const double ql[4], const double qr[4], const double c, double u[3]) 
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

    // Return to sepped of light units
    int status = 0;
    v_hat = cold_sr_vhat(status,&csr);

    // Compute U from Vx
    double u[3];
    calc_u_from_result(u,v_hat,&csr);

    // Do the entropy check
    double UL[3] = {uLx,uLy,uLz};
    double UR[3] = {uRx,uRy,uRz};
    double drho = fabs(rhoR-rhoL);
    double du = fabs(uLx-uRx)+fabs(uLy-uRy)+fabs(uLz-uRz);
    if ((drho + du > 1e-13) && (du > 1e-15)) {
      bool pass_entropy_check = entropy_check(UL, UR, u, c);
      if (pass_entropy_check && status == 0){
      } else {
        printf("Entropy test fails with by, status: %d entropy_test: %d \n", status,pass_entropy_check);
        printf("v_computed = %1.16e; uxhat = %1.16e; uyhat = %1.16e; uzhat = %1.16e;\n",v_hat,u[0],u[1],u[2]);
        printf("rhol = %1.16e; ulx = %1.16e; uly = %1.16e; ulz = %1.16e;\n",rhoL,uLx,uLy,uLz);
        printf("rhor = %1.16e; urx = %1.16e; ury = %1.16e; urz = %1.16e;\n",rhoR,uRx,uRy,uRz);

        // give the momentum averaged velocitity
        double u_avg[3] = {(uLx+uRx)/2.0,(uLy+uRy)/2.0,(uLz+uRz)/2.0};
        double v_hat_avg = u_avg[0]/sqrt(1.0 + u_avg[0]*u_avg[0] + u_avg[1]*u_avg[1] + u_avg[2]*u_avg[2] );
        return c*v_hat_avg;
      }
    } else {
      // give the momentum averaged velocitity
      double u_avg[3] = {(uLx+uRx)/2.0,(uLy+uRy)/2.0,(uLz+uRz)/2.0};
      double v_hat_avg = u_avg[0]/sqrt(1.0 + u_avg[0]*u_avg[0] + u_avg[1]*u_avg[1] + u_avg[2]*u_avg[2] );
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