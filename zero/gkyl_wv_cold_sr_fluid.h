#pragma once

#include <gkyl_wv_eqn.h>
#include <gkyl_math.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>




/**
 * Create a new cold-relativistic-fluid equation object.
 * 
 * @return Pointer to cold-relativistic-fluid equation object.
 */
struct gkyl_wv_eqn* gkyl_wv_cold_sr_fluid_new(void);



/**
 * Compute the jacobian and its inverse, and multiply it by F(U) (= 0) for the
 * newton iteration of a system of equations 
 * 
 * @param u inital momentum guess
 * @param rr density (right)
 * @param rl density (left)
 * @param ur momentum (right)
 * @param ul momentum (left)
 * @param j_inv_fx final computation of the inv. jacobian and F(U) (output)
 * @param c speed of light
 * 
 */
static inline void
j_inverse_times_fx(const double u[3], const double rr, const double rl, 
const double ur[3], const double ul[3], double j_inv_fx[3], const double c)
{

  // save temporary variables
  double fx[3];
  double gamma = sqrt(1.0 + (u[0]*u[0] + u[1]*u[1] + u[2]*u[2])/c*c);
  double gammal = sqrt(1.0 + (ul[0]*ul[0] + ul[1]*ul[1] + ul[2]*ul[2])/c*c);
  double gammar = sqrt(1.0 + (ur[0]*ur[0] + ur[1]*ur[1] + ur[2]*ur[2])/c*c);
  double vlx = ul[0]/gammal;
  double vrx = ur[0]/gammar;

  // compute the elements of the 3x3 Jacobian
  double j00 = -(u[0]*u[0]*u[0]*rl - u[0]*u[0]*u[0]*rr - c*c*rl*ul[0] + c*c*rr*ur[0] + 2.0*u[0]*u[1]*u[1]*rl -
    2.0*u[0]*u[1]*u[1]*rr + 2.0*u[0]*u[2]*u[2]*rl - 2.0*u[0]*u[2]*u[2]*rr + 2.0*u[0]*c*c*rl - 2.0*u[0]*c*c*rr -
    u[1]*u[1]*rl*ul[0] - u[2]*u[2]*rl*ul[0] + u[1]*u[1]*rr*ur[0] + u[2]*u[2]*rr*ur[0] - c*c*gamma*gamma*gamma*rl*vlx +
    c*c*gamma*gamma*gamma*rr*vrx)/(c*c*gamma*gamma*gamma);
  double j01 = (u[0]*u[1]*(u[0]*rl - u[0]*rr - rl*ul[0] + rr*ur[0]))/(c*c*gamma*gamma*gamma);
  double j02 = (u[0]*u[2]*(u[0]*rl - u[0]*rr - rl*ul[0] + rr*ur[0]))/(c*c*gamma*gamma*gamma);

  double j10 = -((u[1]*u[1] + u[2]*u[2] + c*c)*(u[1]*rl - u[1]*rr - rl*ul[1] + rr*ur[1]))/(c*c*gamma*gamma*gamma);
  double j11 = -(u[0]*u[0]*u[0]*rl - u[0]*u[0]*u[0]*rr + u[0]*u[2]*u[2]*rl - u[0]*u[2]*u[2]*rr + u[0]*c*c*rl - 
    u[0]*c*c*rr - c*c*gamma*gamma*gamma*rl*vlx + c*c*gamma*gamma*gamma*rr*vrx + u[0]*u[1]*rl*ul[1] - u[0]*u[1]*rr*ur[1])/(c*c*gamma*gamma*gamma);
  double j12 = (u[0]*u[2]*(u[1]*rl - u[1]*rr - rl*ul[1] + rr*ur[1]))/(c*c*gamma*gamma*gamma);

  double j20 = -((u[1]*u[1] + u[2]*u[2] + c*c)*(u[2]*rl - u[2]*rr - rl*ul[2] + rr*ur[2]))/(c*c*gamma*gamma*gamma);
  double j21 = (u[0]*u[1]*(u[2]*rl - u[2]*rr - rl*ul[2] + rr*ur[2]))/(c*c*gamma*gamma*gamma);
  double j22 = -(u[0]*u[0]*u[0]*rl - u[0]*u[0]*u[0]*rr + u[0]*u[1]*u[1]*rl - u[0]*u[1]*u[1]*rr + u[0]*c*c*rl -
   u[0]*c*c*rr - c*c*gamma*gamma*gamma*rl*vlx + c*c*gamma*gamma*gamma*rr*vrx + u[0]*u[2]*rl*ul[2] - u[0]*u[2]*rr*ur[2])/(c*c*gamma*gamma*gamma);

  // Compute the inverse
  double one_div_det_j = 1.0/(j00*(j11*j22 - j12*j21) - j01*(j10*j22 - j12*j20) + j02*(j10*j21-j11*j20));

  // compute the elements of the 3x3 Jacobian-inverse
  double j_inv00 = one_div_det_j*(j11*j22 - j12*j21);
  double j_inv01 = one_div_det_j*(j02*j21 - j01*j22);
  double j_inv02 = one_div_det_j*(j01*j12 - j02*j11);

  double j_inv10 = one_div_det_j*(j12*j20 - j10*j22);
  double j_inv11 = one_div_det_j*(j00*j22 - j02*j20);
  double j_inv12 = one_div_det_j*(j02*j10 - j00*j12);

  double j_inv20 = one_div_det_j*(j10*j21 - j11*j20);
  double j_inv21 = one_div_det_j*(j01*j20 - j00*j21);
  double j_inv22 = one_div_det_j*(j00*j11 - j01*j10);

  // Compute fx
  double vsx = u[0]/gamma;
  fx[0] = u[0]*( (rr-rl)*vsx + rl*vlx - rr*vrx ) + rl*ul[0]*(vsx - vlx) + rr*ur[0]*(vrx - vsx);
  fx[1] = u[1]*( (rr-rl)*vsx + rl*vlx - rr*vrx ) + rl*ul[1]*(vsx - vlx) + rr*ur[1]*(vrx - vsx);
  fx[2] = u[2]*( (rr-rl)*vsx + rl*vlx - rr*vrx ) + rl*ul[2]*(vsx - vlx) + rr*ur[2]*(vrx - vsx);

  // save j_inv_fx
  j_inv_fx[0] = j_inv00*fx[0] + j_inv01*fx[1] + j_inv02*fx[2];
  j_inv_fx[1] = j_inv10*fx[0] + j_inv11*fx[1] + j_inv12*fx[2];
  j_inv_fx[2] = j_inv20*fx[0] + j_inv21*fx[1] + j_inv22*fx[2];

  //if (isnan(j_inv_fx[0]) || !isfinite(j_inv_fx[0])) {
  //  printf("Nan/Inf (PRE)");
  //}

  //for(int m=0; m<3; ++m) printf("DIAG (PRE): %1.16e ",j_inv_fx[m]);
  //printf("\n");

}


/**
 * Compute the primitive variables (N,NU) -> (N,U)
 * 
 * @param q conserved quantity vector 
 * @param r density (output)
 * @param u momentum (output)
 */
static inline void
isolate_prims(const double q[4], double r, double u[3])
{
  r = q[0];
  u[0] = q[1]/q[0];
  u[1] = q[2]/q[0];
  u[2] = q[3]/q[0];
}


/**
 * Compute the roe-averaged velocity for the sr-cold fluid equations, solved via
 * newton method.
 * 
 * This requries 3 velocity dimensions, and assumes we are calculating the first index 
 * interface (for instance the direction along Vx for Ux, Uy, Uz)
 *
 * @param ql Conserved variables (left)
 * @param qr Conserved variables (right)
 * @param c Speed of light 
 */
static double
compute_sr_roe_avereged_velocity(const double ql[4], const double qr[4], const double c) 
{
  double error = 1.; 
  double tol = 1.e-12;
  size_t iter = 0;
  int vdim = 3;
  double u[3], u0[3];
  double rr, rl;
  double ul[3], ur[3];
  double j_inv_fx[3];
  double v;

  // Isolate N (r) and U (u) from q
  rl = ql[0];
  ul[0] = ql[1]/ql[0];
  ul[1] = ql[2]/ql[0];
  ul[2] = ql[3]/ql[0];
  rr = qr[0];
  ur[0] = qr[1]/qr[0];
  ur[1] = qr[2]/qr[0];
  ur[2] = qr[3]/qr[0];

  // Normalize the density
  rr = qr[0]/ql[0];
  rl = 1.0;

  // check rr and rl are not zero - stop for both zero or either negative densities
  if ( ((rr != 0.0) && (rl != 0.0)) && ((rr >= 0.0) && (rl >= 0.0)) ) {

    // intial guess (non-relativistic rho-averaged velocity)
    for(int m=0; m<3; ++m)
      u0[m] = (sqrt(rl)*ul[m] + sqrt(rr)*ur[m])/(sqrt(rl)+sqrt(rr)); 
    
    while (error > tol || iter < 10) {
      j_inverse_times_fx(u0,rr,rl,ur,ul,j_inv_fx,c);
      //for(int m=0; m<3; ++m) printf("DIAG: %1.16e ",j_inv_fx[m]);
      //printf("\n");
      error = 0.0;
      for(int m=0; m<3; ++m){
        u0[m] = u0[m] - j_inv_fx[m];
        error = fmax(fabs(j_inv_fx[m]),error);
      }
      iter = iter + 1;
      v = u0[0]/sqrt(1.0 + (u0[0]*u0[0] + u0[1]*u0[1] + u0[2]*u0[2]));
      if (iter > 15) error = 0.0; // break the loop
      
      //printf("---> Iteration %ld (%g) \n", iter, ps2);
      //printf(" %lg %lg %lg %lg %lg\n", q[0], q[1], q[2], q[3], q[4]);
    }
    //printf("Iterations %ld. Error %lg \n", iter, fabs(ps2-ps));

  } else {
    // If rho is zero then it doesn't matter what the velocities are so we set them to zero
    v = 0.0;
  }

  // Return the velocity
  return v;
}

// Params: //(Uy,a,Ay,dp,By,Ax,Bx,Az,Bz,c)
static inline double
case_1(double Ux, void *ctx)
{ 
  double *constants = (double*)ctx;
  double a = constants[0];
  double dp = constants[2];
  double Ax = constants[4];
  double Bx = constants[5];
  double c = constants[8];

  // Compute values with old u0
  //double Vx = (-Ux*a - Ax)/(Ux*dp + Bx);

  // Eval of the function
  return pow((Ux*a + Ax),2)*(1.0 + pow(Ux/c,2)) - pow(Ux,2)*pow(Ux*dp + Bx,2);//Vx -  Ux / sqrt(1.0 + (Ux*Ux)/(c*c));
}

// Params: //(Uz,a,Ay,dp,By,Ax,Bx,Az,Bz,c)
static inline double
case_2(double Uz, void *ctx)
{ 
  double *constants = (double*)ctx;
  double a = constants[0];
  double dp = constants[2];
  double Ax = constants[4];
  double Bx = constants[5];
  double Az = constants[6];
  double Bz = constants[7];
  double c = constants[8];

  // Compute values with old u0
  double Vx = (-Uz*a - Az)/(Uz*dp + Bz);
  double Ux = (Ax + Bx*Vx)/(-dp*Vx - a);

  // Eval of the function
  return Vx -  Ux / sqrt(1.0 + (Ux*Ux + Uz*Uz)/(c*c));
}

// Params: //(Uy,a,Ay,dp,By,Ax,Bx,Az,Bz,c)
static inline double
case_3(double Uy, void *ctx)
{ 
  double *constants = (double*)ctx;
  double a = constants[0];
  double Ay = constants[1];
  double dp = constants[2];
  double By = constants[3];
  double Ax = constants[4];
  double Bx = constants[5];
  double c = constants[8];

  // Compute values with old u0
  double Vx = (-Uy*a - Ay)/(Uy*dp + By);
  double Ux = (Ax + Bx*Vx)/(-dp*Vx - a);

  // Eval of the function
  return Vx -  Ux / sqrt(1.0 + (Ux*Ux + Uy*Uy)/(c*c));
}

// Params: //(Uy,a,Ay,dp,By,Ax,Bx,Az,Bz,c)
static inline double
case_4(double Uy, void *ctx)
{ 
  double *constants = (double*)ctx;
  double a = constants[0];
  double Ay = constants[1];
  double dp = constants[2];
  double By = constants[3];
  double Ax = constants[4];
  double Bx = constants[5];
  double Az = constants[6];
  double Bz = constants[7];
  double c = constants[8];

  // Compute values with old u0
  double Vx = (-Uy*a - Ay)/(Uy*dp + By);
  double Ux = (Ax + Bx*Vx)/(-dp*Vx - a);
  double Uz = (Az + Bz*Vx)/(-dp*Vx - a);

  // Eval of the function
  return Vx -  Ux / sqrt(1.0 + (Ux*Ux + Uy*Uy + Uz*Uz)/(c*c));
}


/**
 * Compute the roe-averaged velocity for the sr-cold fluid equations, solved via
 * Ridders' method.
 * 
 * This requries 3 velocity dimensions, and 
 * A. Assumes we are calculating the first index 
 * interface (for instance the direction along Vx for Ux, Uy, Uz)
 * B. Assumes density has already been checked and ql[0] > 0, and qr[0] > 0
 * C. If all velocities are zero, the calculation returns zero
 * D. Non gaurds against div by zero in case_#(...) functions
 * E. c is defined as one.
 * 
 * v_result can be obtained by simple relations because of the contraints in the NL
 * solve
 *
 * @param ql Conserved variables (left)
 * @param qr Conserved variables (right)
 * @param c speed of light
 */
static double
compute_sr_roe_averaged_velocity_via_ridders(const double ql[4], const double qr[4], const double c) 
{

  // speed of light
  double v_result;

  // Isolate varaibles (right/left)
  double rhoR = qr[0];
  double uRx = qr[1]/qr[0];
  double uRy = qr[2]/qr[0];
  double uRz = qr[3]/qr[0];
  double rhoL = ql[0];
  double uLx = ql[1]/ql[0];
  double uLy = ql[2]/ql[0];
  double uLz = ql[3]/ql[0];

  // if left = right then return the left values
  if (qr[0] == ql[0] && qr[1] == ql[1] && qr[2] == ql[2] && qr[3] == ql[3]) {
    v_result = uLx/sqrt(1.0 + (uLx*uLx + uLy*uLy + uLz*uLz));

  // If density is zero on either side
  } else if (qr[0] == 0.0 && ql[0] == 0.0) {
    v_result = 0.0;
  } else if (qr[0] == 0.0) {
    v_result = uLx/sqrt(1.0 + (uLx*uLx + uLy*uLy + uLz*uLz));
  } else if (ql[0] == 0.0) {
    v_result = uRx/sqrt(1.0 + (uRx*uRx + uRy*uRy + uRz*uRz));

  } else {
    // Compute the constants:
    double vLx = uLx/sqrt(1.0 + (uLx*uLx + uLy*uLy + uLz*uLz)/(c*c));
    double vRx = uRx/sqrt(1.0 + (uRx*uRx + uRy*uRy + uRz*uRz)/(c*c));

    // Not division out by rhoL
    rhoR = rhoR/rhoL;
    rhoL = 1.0;
    double dp = rhoR - rhoL;
    double a = rhoL*vLx - rhoR*vRx;
    double Ax = - rhoL*vLx*uLx + rhoR*vRx*uRx;
    double Ay = - rhoL*vLx*uLy + rhoR*vRx*uRy;
    double Az = - rhoL*vLx*uLz + rhoR*vRx*uRz;
    double Bx = rhoL*uLx - rhoR*uRx;
    double By = rhoL*uLy - rhoR*uRy;
    double Bz = rhoL*uLz - rhoR*uRz;

    // Division out by rhoL (doesn't effect the equations)
    double consts[] = {a,Ay,dp,By,Ax,Bx,Az,Bz,c};
    double x1, x2;
    double v_result_default;
  

    //Case 1: Ux /= 0, Uy == 0, Uz == 0, 
    if (uLz == 0.0 && uRz == 0.0 && uLy == 0.0 && uRy == 0.0 && (dp != 0.0 || Bx != 0.0 )) {

      // Compute Ridders, case 1
      //printf("Case 1: "); 
      (uLx > uRx) ? (x2 = uLx, x1 = uRx) : (x1 = uLx, x2 = uRx);
      double f1 = case_1(x1, consts), f2 = case_1(x2, consts);
      struct gkyl_qr_res res = gkyl_ridders(case_1, consts, x1, x2, f1, f2, 100, 1e-12);
      double Ux = res.res;
      v_result = Ux/sqrt(1.0 + Ux*Ux/(c*c));

      // Compute a default result if all other methods fail
      Ux = 0.5*(uLx + uRx); 
      v_result_default = Ux/sqrt(1.0 + Ux*Ux/(c*c));

    // Case 2: Ux /= 0, Uy == 0, Uz /= 0,
    } else if (uLy == 0.0 && uRy == 0.0 && (dp != 0.0 || Bx != 0.0 || Bz != 0.0)){

      // Compute Ridders, case 2
      //printf("Case 2: ");
      (uLz > uRz) ? (x2 = uLz, x1 = uRz) : (x1 = uLz, x2 = uRz);
      double f1 = case_2(x1, consts), f2 = case_2(x2, consts);
      struct gkyl_qr_res res = gkyl_ridders(case_2, consts, x1, x2, f1, f2, 100, 1e-12);
      double Uz = res.res;
      v_result = (-Uz*a - Az)/(Uz*dp + Bz);

    // Case 3: Ux /= 0, Uy /= 0, Uz == 0,
    } else if (uLz == 0.0 && uRz == 0.0 && (dp != 0.0 || Bx != 0.0 || By != 0.0)){

      // Compute Ridders, case 3
      //printf("Case 3: ");
      (uLy > uRy) ? (x2 = uLy, x1 = uRy) : (x1 = uLy, x2 = uRy);
      double f1 = case_3(x1, consts), f2 = case_3(x2, consts);
      struct gkyl_qr_res res = gkyl_ridders(case_3, consts, x1, x2, f1, f2, 100, 1e-12);
      double Uy = res.res;
      v_result = (-Uy*a - Ay)/(Uy*dp + By);

    //Case 4: Ux /= 0, Uy /= 0, Uz /= 0,
    } else if (((uRx != 0.0 ) || (uLx != 0.0 )) && ((uRy != 0.0 ) || (uLy != 0.0 )) && ((uRz != 0.0 ) || (uLz != 0.0 )) && (dp != 0.0 || Bx != 0.0 || By != 0.0 || Bz != 0.0)) {

      // Compute Ridders, case 4
      //printf("Case 4: ");
      (uLy > uRy) ? (x2 = uLy, x1 = uRy) : (x1 = uLy, x2 = uRy);
      double f1 = case_4(x1, consts), f2 = case_4(x2, consts);
      struct gkyl_qr_res res = gkyl_ridders(case_4, consts, x1, x2, f1, f2, 100, 1e-12);
      double Uy = res.res;
      v_result = (-Uy*a - Ay)/(Uy*dp + By);

    } else {

      //Return 0 if the cases all fail
      v_result = 0.0;

    } // end cases


    // if the Ridders fails to converge, run a newton method
    // if both methods fail, then take a central averge
    if (isnan(v_result) || !isfinite(v_result)) {
      double old_result = v_result;
      v_result = compute_sr_roe_avereged_velocity(ql, qr, c);
      //printf("v_result: %1.16f fails to pass, retrying Newton: %1.16f \n",old_result,v_result);

      // If everything fails, its usualy because the tolerance is too low
      if (isnan(v_result) || !isfinite(v_result)) {
        //printf("(***) Defaulting to avg: pl %1.16e, pr %1.16e, ulx %1.16e, urx %1.16e\n",rhoL,rhoR,uLx,uRx);
        v_result = v_result_default;
      }

      //printf("The value is NaN.\n");
      //double Ux_res = (sqrt(rhoL)*uLx + sqrt(rhoR)*uRx)/(sqrt(rhoL)+sqrt(rhoR));
      //v_result = Ux/sqrt(1.0 + Ux*Ux/(c*c));
    }
  }

  return v_result;
}