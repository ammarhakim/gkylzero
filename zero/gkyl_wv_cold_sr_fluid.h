#pragma once

#include <gkyl_wv_sr_eqn.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>




/**
 * Create a new cold-relativistic-fluid equation object.
 * 
 * @return Pointer to cold-relativistic-fluid equation object.
 */
struct gkyl_wv_sr_eqn* gkyl_wv_cold_sr_fluid_new(void);




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
 */
static inline void
j_inverse_times_fx(const double u[3], const double rr, const double rl, 
const double ur[3], const double ul[3], double j_inv_fx[3])
{

  // save temporary variables
  double fx[3];
  double c = 1.0; // assumed, so it can be changed later without redoing the eqs.
  double gamma = sqrt(1.0 + (u[0]^2 + u[1]^2 + u[2]^2)/c^2);

  // compute the elements of the 3x3 Jacobian
  double j00 = -(u[0]^3*rl - u[0]^3*rr - c^2*rl*ul[0] + c^2*rr*ur[0] + 2*u[0]*u[1]^2*rl -
    2*u[0]*u[1]^2*rr + 2*u[0]*u[2]^2*rl - 2*u[0]*u[2]^2*rr + 2*u[0]*c^2*rl - 2*u[0]*c^2*rr -
    u[1]^2*rl*ul[0] - u[2]^2*rl*ul[0] + u[1]^2*rr*ur[0] + u[2]^2*rr*ur[0] - c^2*gamma^3*rl*vlx +
    c^2*gamma^3*rr*vrx)/(c^2*gamma^3);
  double j01 = (u[0]*u[1]*(u[0]*rl - u[0]*rr - rl*ul[0] + rr*ur[0]))/(c^2*gamma^3);
  double j02 = (u[0]*u[2]*(u[0]*rl - u[0]*rr - rl*ul[0] + rr*ur[0]))/(c^2*gamma^3);

  double j10 = -((u[1]^2 + u[2]^2 + c^2)*(u[1]*rl - u[1]*rr - rl*ul[1] + rr*ur[1]))/(c^2*gamma^3);
  double j11 = -(u[0]^3*rl - u[0]^3*rr + u[0]*u[2]^2*rl - u[0]*u[2]^2*rr + u[0]*c^2*rl - 
    u[0]*c^2*rr - c^2*gamma^3*rl*vlx + c^2*gamma^3*rr*vrx + u[0]*u[1]*rl*ul[1] - u[0]*u[1]*rr*ur[1])/(c^2*gamma^3);
  double j12 = (u[0]*u[2]*(u[1]*rl - u[1]*rr - rl*ul[1] + rr*ur[1]))/(c^2*gamma^3);

  double j20 = -((u[1]^2 + u[2]^2 + c^2)*(u[2]*rl - u[2]*rr - rl*ul[2] + rr*ur[2]))/(c^2*gamma^3);
  double j21 = (u[0]*u[1]*(u[2]*rl - u[2]*rr - rl*ul[2] + rr*ur[2]))/(c^2*gamma^3);
  double j22 = -(u[0]^3*rl - u[0]^3*rr + u[0]*u[1]^2*rl - u[0]*u[1]^2*rr + u[0]*c^2*rl -
   u[0]*c^2*rr - c^2*gamma^3*rl*vlx + c^2*gamma^3*rr*vrx + u[0]*u[2]*rl*ul[2] - u[0]*u[2]*rr*ur[2])/(c^2*gamma^3);

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
  double vlx = ul[0]/sqrt(1.0 + (ul[0]^2 + ul[1]^2 + ul[2]^2)/c^2);
  double vrx = ur[0]/sqrt(1.0 + (ur[0]^2 + ur[1]^2 + ur[2]^2)/c^2);
  f1 = u[0]*( (rr-rl)*vsx + rl*vlx - rr*vrx ) + rl*ul[0]*(vsx - vlx) + rr*ur[0]*(vrx - vsx);
  f2 = u[1]*( (rr-rl)*vsx + rl*vlx - rr*vrx ) + rl*ul[1]*(vsx - vlx) + rr*ur[1]*(vrx - vsx);
  f3 = u[2]*( (rr-rl)*vsx + rl*vlx - rr*vrx ) + rl*ul[2]*(vsx - vlx) + rr*ur[2]*(vrx - vsx);

  // save j_inv_fx
  j_inv_fx[0] = j_inv00*fx[0] + j_inv01*fx[1] + j_inv02*fx[2];
  j_inv_fx[1] = j_inv10*fx[0] + j_inv11*fx[1] + j_inv12*fx[2];
  j_inv_fx[2] = j_inv20*fx[0] + j_inv21*fx[1] + j_inv22*fx[2];

}


/**
 * Compute the primitive variables (N,NU) -> (N,U)
 * 
 * @param q conserved quantity vector 
 * @param r desnity (output)
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
 */
static double
compute_sr_roe_avereged_velocity(const double ql[4], const double qr[4]) 
{
  double error = 1.; 
  double tol = 1.e-12;
  size_t iter = 0;
  int vdim = 3;
  double u[3], u0[3];
  double rr, rl;
  double ul[3], ur[3];
  double j_inv_fx[3];

  // Isolate N (r) and U (u) from q
  isolate_prims(ql,rl,ul); // TODO: make this function
  isolate_prims(qr,rr,ur);

  // check rr and rl are not zero - stop for both zero or either negative densities
  if ( ((rr != 0.0) && (rl != 0.0)) && ((rr >= 0.0) && (rl >= 0.0)) ) {

    // intial guess (non-relativistic rho-averaged velocity)
    for(int m=0; m<3; ++m)
      u0[m] = (sqrt(rl)*ul[m] + sqrt(rr)*ur[m])/(sqrt(rl)+sqrt(rr)); 
    
    while (error > tol) {

      j_inverse_times_fx(u0,rr,rl,ur,ul,j_inv_fx);
      error = 0.0;
      for(int m=0; m<3; ++m){
        u0[m] = u0[m] - j_inv_fx[m];
        error = max(abs(j_inv_fx[m]),error);
      }
      iter = iter + 1;
      v = u0[0]/sqrt(1.0 + (u0[0]^2 + u0[1]^2 + u0[2]^2));
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