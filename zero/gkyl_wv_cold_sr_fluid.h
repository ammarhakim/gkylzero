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
  double Vx = (-Ux*a - Ax)/(Ux*dp + Bx);

  // Eval of the function
  return Vx -  Ux / sqrt(1.0 + (Ux*Ux)/(c*c));
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
 */
static double
compute_sr_roe_averaged_velocity_via_ridders(const double ql[4], const double qr[4]) 
{

  // speed of light
  double c = 1.0;

  // Isolate varaibles (right/left)
  double rhoR = qr[0];
  double uRx = qr[1]/qr[0];
  double uRy = qr[2]/qr[0];
  double uRz = qr[3]/qr[0];
  double rhoL = ql[0];
  double uLx = ql[1]/ql[0];
  double uLy = ql[2]/ql[0];
  double uLz = ql[3]/ql[0];

  // Compute the constants:
  double vLx = uLx/sqrt(1.0 + (uLx*uLx + uLy*uLy + uLz*uLz)/(c*c));
  double vRx = uRx/sqrt(1.0 + (uRx*uRx + uRy*uRy + uRz*uRz)/(c*c));
  double dp = rhoR - rhoL;
  double a = rhoL*vLx - rhoR*vRx;
  double Ax = - rhoL*vLx*uLx + rhoR*vRx*uRx;
  double Ay = - rhoL*vLx*uLy + rhoR*vRx*uRy;
  double Az = - rhoL*vLx*uLz + rhoR*vRx*uRz;
  double Bx = rhoL*uLx - rhoR*uRx;
  double By = rhoL*uLy - rhoR*uRy;
  double Bz = rhoL*uLz - rhoR*uRz;
  double consts[] = {a,Ay,dp,By,Ax,Bx,Az,Bz,c};
  double x1, x2;
  double v_result;

  //Case 1: Ux /= 0, Uy == 0, Uz == 0, 
  if (uLz == 0.0 && uRz == 0.0 && uLy == 0.0 && uRy == 0.0) {

    // Compute Ridders, case 1
    (uLx > uRx) ? (x2 = uLx, x1 = uRx) : (x1 = uLx, x2 = uRx);
    double f1 = case_1(x1, consts), f2 = case_1(x2, consts);
    struct gkyl_qr_res res = gkyl_ridders(case_1, consts, x1, x2, f1, f2, 100, 1e-12);
    double Ux = res.res;
    v_result = (-Ux*a - Ax)/(Ux*dp + Bx);

  // Case 2: Ux /= 0, Uy == 0, Uz /= 0,
  } else if (uLy == 0.0 && uRy == 0.0){

    // Compute Ridders, case 2
    (uLz > uRz) ? (x2 = uLz, x1 = uRz) : (x1 = uLz, x2 = uRz);
    double f1 = case_2(x1, consts), f2 = case_2(x2, consts);
    struct gkyl_qr_res res = gkyl_ridders(case_2, consts, x1, x2, f1, f2, 100, 1e-12);
    double Uz = res.res;
    v_result = (-Uz*a - Az)/(Uz*dp + Bz);

  // Case 3: Ux /= 0, Uy /= 0, Uz == 0,
  } else if (uLz == 0.0 && uRz == 0.0){

    // Compute Ridders, case 3
    (uLy > uRy) ? (x2 = uLy, x1 = uRy) : (x1 = uLy, x2 = uRy);
    double f1 = case_3(x1, consts), f2 = case_3(x2, consts);
    struct gkyl_qr_res res = gkyl_ridders(case_3, consts, x1, x2, f1, f2, 100, 1e-12);
    double Uy = res.res;
    v_result = (-Uy*a - Ay)/(Uy*dp + By);

  //Case 4: Ux /= 0, Uy /= 0, Uz /= 0,
  } else if (((uRx != 0.0 ) || (uLx != 0.0 )) && ((uRy != 0.0 ) || (uLy != 0.0 )) && ((uRz != 0.0 ) || (uLz != 0.0 ))) {

    // Compute Ridders, case 4
    (uLy > uRy) ? (x2 = uLy, x1 = uRy) : (x1 = uLy, x2 = uRy);
    double f1 = case_4(x1, consts), f2 = case_4(x2, consts);
    struct gkyl_qr_res res = gkyl_ridders(case_4, consts, x1, x2, f1, f2, 100, 1e-12);
    double Uy = res.res;
    v_result = (-Uy*a - Ay)/(Uy*dp + By);

  } else {

    //Return 0 if the cases all fail
    v_result = 0.0;

  } // end cases

  return v_result;
}