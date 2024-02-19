
#include <acutest.h>
#include <wv_cold_sr_fluid.c>
#include <math.h>

void
calcq(const double c, const double pv[4], double *q)
{
  double rho = pv[0], u = pv[1], v = pv[2], w = pv[3];
  double gamma = 1 / sqrt(1 + (- u*u - v*v - w*w)/(c*c));
  double rhoh = rho;
  q[0] = gamma*rho;
  q[1] = gamma*gamma*rhoh*u;
  q[2] = gamma*gamma*rhoh*v;
  q[3] = gamma*gamma*rhoh*w;
  
}

void
calcu(double q[4], double *u)
{
  
  u[0] = q[0];
  u[1] = q[1]/q[0];
  u[2] = q[2]/q[0];
  u[3] = q[3]/q[0];
  
}

void
compute_flux_jacobian(double A1[4][4], double *u, double c)
{
  // grab elements
  double Ux = u[0];
  double Uy = u[1];
  double Uz = u[2];
  double gamma = sqrt(1.0 + (Ux*Ux + Uy*Uy + Uz*Uz)/(c*c));
  double a = c*c*gamma*gamma*gamma;
  double Ux2 = Ux * Ux;
  double Uy2 = Uy * Uy;
  double Uz2 = Uz * Uz;
  double c2 = c*c;

  A1[0][0] = (Ux * Uz2 + Ux * Uy2 + pow(Ux, 3.0)) / a;
  A1[0][1] = (c2 + Uz2 + Uy2) / a;
  A1[0][2] = - (Ux * Uy) / a;
  A1[0][3] = - (Ux * Uz) / a;

  A1[1][0] = - Ux2 / pow(gamma, 3.0);
  A1[1][1] = (2.0 * Ux * c2 + 2.0 * Ux * Uz2 + 2.0 * Ux * Uy2 + pow(Ux, 3.0)) / a;
  A1[1][2] = - (Ux2 * Uy) / a;
  A1[1][3] = - (Ux2 * Uz) / a;

  A1[2][0] = - (Ux * Uy) / pow(gamma, 3.0);
  A1[2][1] = (Uy * c2 + Uy * Uz2 + pow(Uy, 3.0)) / a;
  A1[2][2] = (Ux * c2 + Ux * Uz2 + pow(Ux, 3.0)) / a;
  A1[2][3] = - (Ux * Uy * Uz) / a;

  A1[3][0] = - (Ux * Uz) / pow(gamma, 3.0);
  A1[3][1] = (Uz * c2 + Uz2 * Uz + Uy2 * Uz) / a; //Had an error!
  A1[3][2] = - (Ux * Uy * Uz) / a;
  A1[3][3] = (Ux * c2 + Ux * Uy2 + pow(Ux, 3.0)) / a;
}

void
compute_approximated_jacobian_c_0(double A1[4][4], double ql[4], double qr[4], double c)
{
  // isolate left and right states
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

  // Primative rho
  double rhol_prim = rhol/gammal;
  double rhor_prim = rhor/gammar;

  // Compute the primative-parameterization state vector w
  // these are the averages of the left and right states
  double k = sqrt(rhol_prim) + sqrt(rhor_prim);
  double w0 = sqrt(rhol_prim)*gammal + sqrt(rhor_prim)*gammar;
  double w1 = sqrt(rhol_prim)*gammal*vlx + sqrt(rhor_prim)*gammar*vrx;
  double w2 = sqrt(rhol_prim)*gammal*vly + sqrt(rhor_prim)*gammar*vry;
  double w3 = sqrt(rhol_prim)*gammal*vlz + sqrt(rhor_prim)*gammar*vrz;

  A1[0][0] = (w1 * (-k * k + w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3)) / (w0 * (k * k + w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3));
  A1[0][1] = (k * (k * k + w0 * w0 - w1 * w1 + w2 * w2 + w3 * w3)) / (w0 * (k * k + w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3));
  A1[0][2] = -(2 * k * w1 * w2) / (w0 * (k * k + w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3));
  A1[0][3] = -(2 * k * w1 * w3) / (w0 * (k * k + w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3));

  A1[1][0] = -(2 * k * w1 * w1) / (w0 * (k * k + w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3));
  A1[1][1] = (2 * w1 * (k * k + w0 * w0 + w2 * w2 + w3 * w3)) / (w0 * (k * k + w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3));
  A1[1][2] = -(2 * w1 * w1 * w2) / (w0 * (k * k + w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3));
  A1[1][3] = -(2 * w1 * w1 * w3) / (w0 * (k * k + w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3));

  A1[2][0] = -(2 * k * w1 * w2) / (w0 * (k * k + w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3));
  A1[2][1] = (w2 * (k * k + w0 * w0 - w1 * w1 + w2 * w2 + w3 * w3)) / (w0 * (k * k + w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3));
  A1[2][2] = (w1 * (k * k + w0 * w0 + w1 * w1 - w2 * w2 + w3 * w3)) / (w0 * (k * k + w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3));
  A1[2][3] = -(2 * w1 * w2 * w3) / (w0 * (k * k + w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3));

  A1[3][0] = -(2 * k * w1 * w3) / (w0 * (k * k + w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3));
  A1[3][1] = (w3 * (k * k + w0 * w0 - w1 * w1 + w2 * w2 + w3 * w3)) / (w0 * (k * k + w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3));
  A1[3][2] = -(2 * w1 * w2 * w3) / (w0 * (k * k + w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3));
  A1[3][3] = (w1 * (k * k + w0 * w0 + w1 * w1 + w2 * w2 - w3 * w3)) / (w0 * (k * k + w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3));
}


void
compute_approximated_jacobian(double A[4][4], double ql[4], double qr[4], double c)
{
  // isolate left and right states
  double rhor = qr[0];
  double urx = qr[1]/(qr[0]);
  double ury = qr[2]/(qr[0]);
  double urz = qr[3]/(qr[0]);
  double rhol = ql[0];
  double ulx = ql[1]/(ql[0]);
  double uly = ql[2]/(ql[0]);
  double ulz = ql[3]/(ql[0]);

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

  // Compute the primative-parameterization state vector w
  // these are the averages of the left and right states
  double k = (sqrt(rhol_prim) + sqrt(rhor_prim))/(c);
  double w0 = sqrt(rhol_prim)*gammal + sqrt(rhor_prim)*gammar;
  double w1 = sqrt(rhol_prim)*gammal*vlx/c + sqrt(rhor_prim)*gammar*vrx/c;
  double w2 = sqrt(rhol_prim)*gammal*vly/c + sqrt(rhor_prim)*gammar*vry/c;
  double w3 = sqrt(rhol_prim)*gammal*vlz/c + sqrt(rhor_prim)*gammar*vrz/c;

  A[0][0] = (c * w1 * (-c * c * k * k + w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3)) / (w0 * (c * c * k * k + w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3));
  A[0][1] = (c * k * (c * c * k * k + w0 * w0 - w1 * w1 + w2 * w2 + w3 * w3)) / (w0 * (c * c * k * k + w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3));
  A[0][2] = -(2 * c * k * w1 * w2) / (w0 * (c * c * k * k + w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3));
  A[0][3] = -(2 * c * k * w1 * w3) / (w0 * (c * c * k * k + w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3));

  A[1][0] = -(2 * c * c * c * k * w1 * w1) / (w0 * (c * c * k * k + w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3));
  A[1][1] = (2 * c * w1 * (c * c * k * k + w0 * w0 + w2 * w2 + w3 * w3)) / (w0 * (c * c * k * k + w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3));
  A[1][2] = -(2 * c * w1 * w1 * w2) / (w0 * (c * c * k * k + w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3));
  A[1][3] = -(2 * c * w1 * w1 * w3) / (w0 * (c * c * k * k + w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3));

  A[2][0] = -(2 * c * c * c * k * w1 * w2) / (w0 * (c * c * k * k + w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3));
  A[2][1] = (c * w2 * (c * c * k * k + w0 * w0 - w1 * w1 + w2 * w2 + w3 * w3)) / (w0 * (c * c * k * k + w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3));
  A[2][2] = (c * w1 * (c * c * k * k + w0 * w0 + w1 * w1 - w2 * w2 + w3 * w3)) / (w0 * (c * c * k * k + w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3));
  A[2][3] = -(2 * c * w1 * w2 * w3) / (w0 * (c * c * k * k + w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3));

  A[3][0] = -(2 * c * c * c * k * w1 * w3) / (w0 * (c * c * k * k + w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3));
  A[3][1] = (c * w3 * (c * c * k * k + w0 * w0 - w1 * w1 + w2 * w2 + w3 * w3)) / (w0 * (c * c * k * k + w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3));
  A[3][2] = -(2 * c * w1 * w2 * w3) / (w0 * (c * c * k * k + w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3));
  A[3][3] = (c * w1 * (c * c * k * k + w0 * w0 + w1 * w1 + w2 * w2 - w3 * w3)) / (w0 * (c * c * k * k + w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3));

}


void
compute_prim_approximated_jacobian(double A[4][4], double ql[4], double qr[4], double c)
{
  // isolate left and right states
  double rhor = qr[0];
  double urx = qr[1]/(qr[0]);
  double ury = qr[2]/(qr[0]);
  double urz = qr[3]/(qr[0]);
  double rhol = ql[0];
  double ulx = ql[1]/(ql[0]);
  double uly = ql[2]/(ql[0]);
  double ulz = ql[3]/(ql[0]);

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

  // Compute the primative-parameterization state vector w
  // these are the averages of the left and right states
  double k = (sqrt(rhol_prim) + sqrt(rhor_prim))/(c);
  double w0 = sqrt(rhol_prim)*gammal + sqrt(rhor_prim)*gammar;
  double w1 = sqrt(rhol_prim)*gammal*vlx/c + sqrt(rhor_prim)*gammar*vrx/c;
  double w2 = sqrt(rhol_prim)*gammal*vly/c + sqrt(rhor_prim)*gammar*vry/c;
  double w3 = sqrt(rhol_prim)*gammal*vlz/c + sqrt(rhor_prim)*gammar*vrz/c;

  // Compute A-tilde 
  double denominator = c*c*k*k + w0*w0 + w1*w1 + w2*w2 + w3*w3;
  double kstar = sqrt(rhol_prim)*sqrt(rhor_prim)/(c*c);

  // First row
  A[0][0] = (2*c*w0*w1) / denominator;
  A[0][1] = (c*c*kstar*(c*c*k*k + w0*w0 - w1*w1 + w2*w2 + w3*w3)) / denominator;
  A[0][2] = -(2*c*c*kstar*w1*w2) / denominator;
  A[0][3] = -(2*c*c*kstar*w1*w3) / denominator;

  // Second row
  A[1][0] = (2*c*w0*w1*w1) / (k*denominator);
  A[1][1] = (2*c*c*kstar*w1*(c*c*k*k + w0*w0 + w2*w2 + w3*w3)) / (k*denominator);
  A[1][2] = -(2*c*c*kstar*w1*w1*w2) / (k*denominator);
  A[1][3] = -(2*c*c*kstar*w1*w1*w3) / (k*denominator);

  // Third row
  A[2][0] = (2*c*w0*w1*w2) / (k*denominator);
  A[2][1] = (c*c*kstar*w2*(c*c*k*k + w0*w0 - w1*w1 + w2*w2 + w3*w3)) / (k*denominator);
  A[2][2] = (c*c*kstar*w1*(c*c*k*k + w0*w0 + w1*w1 - w2*w2 + w3*w3)) / (k*denominator);
  A[2][3] = -(2*c*c*kstar*w1*w2*w3) / (k*denominator);

  // Fourth row
  A[3][0] = (2*c*w0*w1*w3) / (k*denominator);
  A[3][1] = (c*c*kstar*w3*(c*c*k*k + w0*w0 - w1*w1 + w2*w2 + w3*w3)) / (k*denominator);
  A[3][2] = -(2*c*c*kstar*w1*w2*w3) / (k*denominator);
  A[3][3] = (c*c*kstar*w1*(c*c*k*k + w0*w0 + w1*w1 + w2*w2 - w3*w3)) / (k*denominator);

}


void 
multiplyMatrixVector(double A[4][4], double b[4], double *result) {
  for (int i = 0; i < 4; ++i) {
    result[i] = 0.0;
    for (int j = 0; j < 4; ++j) {
      result[i] += A[i][j] * b[j];
    }
  }
}


void 
test_flux_jump()
{

  // TEST: A(u-hat)(qr - ql) = F(qr) - F(ql) 

  // Averages are computed via:
  // Martí, José María, and Ewald Müller. "Grid-based methods 
  // in relativistic hydrodynamics and magnetohydrodynamics." 
  // Living reviews in computational astrophysics 1.1 (2015): 1-182.

  // Give constants for the test:
  double c = 1.0;
  c = 299792458.0;

  double vl[4] = { 1.0, 0.0999*c, 0.02*c, 0.03*c };
  double vr[4] = { 0.1, 0.7*c, 0.2*c, 0.3*c  };
  double ql[4], qr[4], ul[4], ur[4];
  calcq(c, vl, ql); calcq(c, vr, qr);
  calcu(ql, ul); calcu(qr, ur);
  double u_hat[3];

  // Compute the flux jacobian:
  double A_flux_jacobian[4][4];
  double A_tilde[4][4];
  double A_tilde_prim[4][4];
  compute_approximated_jacobian(A_tilde, ql, qr, c);

  // Compute fluxes
  double fl[4], fr[4];
  cold_sr_fluid_flux(ql,fl);
  cold_sr_fluid_flux(qr,fr);

  // Compute f(ql) - f(qr) and (qr - ql)
  double dfq[4], dq[4], du[4];
  for(int i=0; i<4; ++i) {
    dfq[i] =  fr[i] - fl[i];
    dq[i] =  qr[i] - ql[i];
    du[i] =  ur[i] - ul[i];
  }

  double Adq[4], Adu[4];

  // compute A(u_hat)(qr - ql) - ( f(ql) - f(qr) ) ~ 0 (retruns error in comp)
  double error[4], amdq[4], apdq[4];
  multiplyMatrixVector(A_tilde, dq, Adq);
  multiplyMatrixVector(A_tilde, ql, amdq);
  multiplyMatrixVector(A_tilde, qr, apdq);
  for(int i=0; i<4; ++i){
    error[i] = Adq[i] - dfq[i];
    amdq[i] = -amdq[i];
  } 

  // Condition (ii) as ul -> ur -> u, then A(ul,ur)_tilde -> A(u)
  // take ul := ur as the state u
  double ur_cons[3];
  double diff_A_tilde_A_cond2[4][4];
  ur_cons[0] = qr[1]/qr[0];
  ur_cons[1] = qr[2]/qr[0];
  ur_cons[2] = qr[3]/qr[0];
  compute_flux_jacobian(A_flux_jacobian, ur_cons, c);
  compute_approximated_jacobian(A_tilde, qr, qr, c);
   for(int i=0; i<4; ++i)
     for(int j=0; j<4; ++j)
      diff_A_tilde_A_cond2[i][j] = A_flux_jacobian[i][j] - A_tilde[i][j];

  for(int i=0; i<4; ++i)
     for(int j=0; j<4; ++j)
      TEST_CHECK( gkyl_compare((diff_A_tilde_A_cond2[i][j]/A_tilde[i][j]), 0, 1e-12) );


  // (PRIMATIVE) Condition (ii) as ul -> ur -> u, then A(ul,ur)_tilde -> A(u)
  // take ul := ur as the state u
  compute_approximated_jacobian(A_tilde, ql, qr, c);
  compute_prim_approximated_jacobian(A_tilde_prim, ql, qr, c);
  multiplyMatrixVector(A_tilde_prim, du, Adu);
  multiplyMatrixVector(A_tilde, dq, Adq);

   for(int i=0; i<4; ++i)
    TEST_CHECK( gkyl_compare((Adq[i] - Adu[i])/Adu[i], 0, 1e-12) );
}


void
test_sr_cold_fluid()
{

    // Give constants for the test:
  //const double c = 1.0;

  struct gkyl_wv_eqn *sr_cold_fluid = gkyl_wv_cold_sr_fluid_new();
  const double c = 299792458.0;
  const double n0 = 2.0e24;
  double vl[4] = { 1.0*n0, 0.0999*c, 0.02*c, 0.03*c };
  double vr[4] = { 0.1*n0, 0.7*c, 0.2*c, 0.3*c  };
  double ql[4], qr[4];
  double ql_local[4], qr_local[4];
  calcq(c, vl, ql); calcq(c, vr, qr);

  double norm[3][3] = {
    { 1.0, 0.0, 0.0 },
    { 0.0, -1.0, 0.0 },
    { 0.0, 0.0, 1.0 }
  };

  double tau1[3][3] = {
    { 0.0, 1.0, 0.0 },
    { 1.0, 0.0, 0.0 },
    { 1.0, 0.0, 0.0 }
  };

  double tau2[3][3] = {
    { 0.0, 0.0, 1.0 },
    { 0.0, 0.0, 1.0 },
    { 0.0, 1.0, 0.0 }
  };  

  for (int d=0; d<3; ++d) {
    double speeds[2], waves[2*4], waves_local[2*4];
    // rotate to local tangent-normal frame
    gkyl_wv_eqn_rotate_to_local(sr_cold_fluid, tau1[d], tau2[d], norm[d], ql, ql_local);
    gkyl_wv_eqn_rotate_to_local(sr_cold_fluid, tau1[d], tau2[d], norm[d], qr, qr_local);

    // check if sum of left/right going fluctuations sum to jump in flux
    double fl_local[4], fr_local[4];
    cold_sr_fluid_flux(ql_local, fl_local);
    cold_sr_fluid_flux(qr_local, fr_local);

    double delta[4];
    for (int i=0; i<4; ++i) delta[i] = fr_local[i]-fl_local[i];
    
    gkyl_wv_eqn_waves(sr_cold_fluid, GKYL_WV_HIGH_ORDER_FLUX, delta, ql_local, qr_local, waves_local, speeds);

    // rotate waves back to global frame
    for (int mw=0; mw<2; ++mw)
      gkyl_wv_eqn_rotate_to_global(sr_cold_fluid, tau1[d], tau2[d], norm[d], &waves_local[mw*4], &waves[mw*4]);

    double apdq[4], amdq[4];
    gkyl_wv_eqn_ffluct(sr_cold_fluid, GKYL_WV_HIGH_ORDER_FLUX, ql, qr, waves, speeds, amdq, apdq);
  
    double fl[4], fr[4];
    gkyl_wv_eqn_rotate_to_global(sr_cold_fluid, tau1[d], tau2[d], norm[d], fl_local, fl);
    gkyl_wv_eqn_rotate_to_global(sr_cold_fluid, tau1[d], tau2[d], norm[d], fr_local, fr);
    
    for (int i=0; i<4; ++i){
      //printf("fr[%d]: %1.4e, fl[%d]: %1.4e, amdq[%d]: %1.4e, apdq[%d]: %1.4e ",i,fr[i],i,fl[i],i,amdq[i],i,apdq[i]);
      //printf("  df[%d]: %1.4e, dadq[%d]: %1.4e,     total_diff: %1.4e\n",i,fr[i]-fl[i],i,amdq[i]+apdq[i], fr[i]-fl[i]-(amdq[i]+apdq[i]));
      TEST_CHECK( gkyl_compare((fr[i]-fl[i]), (amdq[i]+apdq[i]), 1e-8) ); 
    }
  }

  gkyl_wv_eqn_release(sr_cold_fluid);
}


void test_flux_jump_relation() { test_flux_jump(); }
void test_sr_cold_fluids_internal() { test_sr_cold_fluid(); }

TEST_LIST = {
  {"test_cold_sr_fluids_flux_jump", test_flux_jump_relation},
  {"test_cold_sr_fluids_ffluct", test_sr_cold_fluids_internal},
  {NULL, NULL},
};