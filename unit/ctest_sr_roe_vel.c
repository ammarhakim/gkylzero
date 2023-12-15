
#include <acutest.h>
#include <wv_cold_sr_fluid.c>
#include <math.h>

void
calcq(const double c, const double pv[4], double q[4])
{
  double rho = pv[0], u = pv[1], v = pv[2], w = pv[3];
  double gamma = 1 / sqrt(1 - u*u - v*v - w*w);
  double rhoh = rho;
  printf("v/c: %1.4e\n", sqrt(u*u + v*v + w*w) );
  
  q[0] = gamma*rho;
  q[1] = gamma*gamma*rhoh*u;
  q[2] = gamma*gamma*rhoh*v;
  q[3] = gamma*gamma*rhoh*w;
  
}

void
cold_sr_fluid_flux_local_copy(const double q[4], double *flux)
{
  // Vx = NUx/sqrt(N^2 + NU^2/c^2) 
  const double c = 1.0; 
  double Vx = q[NUX]/sqrt(q[0]*q[0] + (q[NUX]*q[NUX] + q[NUY]*q[NUY] + q[NUZ]*q[NUZ])/(c*c)); 
  flux[0] = q[0]*Vx; // N*Vx
  flux[NUX] = q[NUX]*Vx; // N*Ux*Vx
  flux[NUY] = q[NUY]*Vx; // N*Uy*Vx
  flux[NUZ] = q[NUZ]*Vx; // N*Uz*Vx
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
compute_approximated_jacobian(double A1[4][4], double ql[4], double qr[4], double c)
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

  // CASE 1: 

  // Give constants for the test:
  const double c = 1.0;

  double vl[4] = { 1.0, 0.0999, 0.02, 0.03 };
  double vr[4] = { 0.1, 0.7, 0.2, 0.3  };
  double ql[4], qr[4];
  calcq(c, vl, ql); calcq(c, vr, qr);

  //const double ql[4] = { 1.0, 15.0*c, 0.01*c, -5.0*c };
  //const double qr[4] = { 0.1, -0.5*c,  0.01*c, 25.0*c };
  double u_hat[3];
  double gamma_hat = 0.0;


  // Call the roe-averaged velocity routine
  double v_avg = compute_sr_roe_averaged_velocity_cold_limit(ql, qr, c, u_hat, &gamma_hat);
  // double v_avg = compute_sr_roe_averaged_velocity(ql, qr, c, u_hat);


  // Compute the flux jacobian:
  double A_flux_jacobian[4][4];
  double A_tilde[4][4];
  compute_flux_jacobian(A_flux_jacobian, u_hat, c);
  compute_approximated_jacobian(A_tilde, ql, qr, c);

  // Compute fluxes
  double fl[4], fr[4];
  cold_sr_fluid_flux(ql,fl);
  cold_sr_fluid_flux(qr,fr);

  // Compute f(ql) - f(qr) and (qr - ql)
  double dfq[4], dq[4];
  for(int i=0; i<4; ++i) {
    dfq[i] =  fr[i] - fl[i];
    dq[i] =  qr[i] - ql[i];
  }

  // Compute A(u_hat)(qr - ql)
  double Adq[4];
  multiplyMatrixVector(A_flux_jacobian, dq, Adq);

  // compute A(u_hat)(qr - ql) - ( f(ql) - f(qr) ) ~ 0 (retruns error in comp)
  double error[4], amdq[4], apdq[4];
  multiplyMatrixVector(A_flux_jacobian, ql, amdq);
  multiplyMatrixVector(A_flux_jacobian, qr, apdq);
  printf("\n");
  for(int i=0; i<4; ++i){
    error[i] = Adq[i] - dfq[i];
    amdq[i] = -amdq[i];
    printf("(Flux_Jacobian Method) fr[%d]: %1.4e, fl[%d]: %1.4e, amdq[%d]: %1.4e, apdq[%d]: %1.4e ",i,fr[i],i,fl[i],i,amdq[i],i,apdq[i]);
    printf("  df[%d]: %1.4e, dadq[%d]: %1.4e,     total_diff: %1.4e\n",i,fr[i]-fl[i],i,amdq[i]+apdq[i], fr[i]-fl[i]-(amdq[i]+apdq[i]));
  } 
  multiplyMatrixVector(A_tilde, dq, Adq);
  multiplyMatrixVector(A_tilde, ql, amdq);
  multiplyMatrixVector(A_tilde, qr, apdq);
  printf("\n");
  for(int i=0; i<4; ++i){
    error[i] = Adq[i] - dfq[i];
    amdq[i] = -amdq[i];
    printf("(Approximated A-tilde) fr[%d]: %1.4e, fl[%d]: %1.4e, amdq[%d]: %1.4e, apdq[%d]: %1.4e ",i,fr[i],i,fl[i],i,amdq[i],i,apdq[i]);
    printf("  df[%d]: %1.4e, dadq[%d]: %1.4e,     total_diff: %1.4e\n",i,fr[i]-fl[i],i,amdq[i]+apdq[i], fr[i]-fl[i]-(amdq[i]+apdq[i]));
  } 

    // Condition (ii) as ul -> ur -> u, then A(ul,ur)_tilde -> A(u)
  // take ul := ur as the state u
  double ur[3];
  double diff_A_tilde_A_cond2[4][4];
  ur[0] = qr[1]/qr[0];
  ur[1] = qr[2]/qr[0];
  ur[2] = qr[3]/qr[0];
  compute_flux_jacobian(A_flux_jacobian, ur, c);
  compute_approximated_jacobian(A_tilde, qr, qr, c);
  printf("------------\n\n");
  printf("A: \n");
   for(int i=0; i<4; ++i){
    printf("[");
     for(int j=0; j<4; ++j){
      printf("%1.4e " ,A_flux_jacobian[i][j]);
      diff_A_tilde_A_cond2[i][j] = A_flux_jacobian[i][j] - A_tilde[i][j];
     }
     printf("]\n");
   }
   printf("\n");
   printf("A_tilde: \n");
   for(int i=0; i<4; ++i){
    printf("[");
     for(int j=0; j<4; ++j){
      printf("%1.4e " ,A_tilde[i][j]);
     }
     printf("]\n");
   }
   printf("\n");
      printf("error: A-A_tilde: \n");
   for(int i=0; i<4; ++i){
    printf("[");
     for(int j=0; j<4; ++j){
      printf("%1.4e " ,diff_A_tilde_A_cond2[i][j]);
     }
     printf("]\n");
   }
   printf("\n");

}


void
test_sr_cold_fluid()
{

    // Give constants for the test:
  const double c = 1.0;

  struct gkyl_wv_eqn *sr_cold_fluid = gkyl_wv_cold_sr_fluid_new();
  double vl[4] = { 1.0, 0.0999, 0.02, 0.03 };
  double vr[4] = { 0.1, 0.7, 0.2, 0.3  };
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
    cold_sr_fluid_flux_local_copy(ql_local, fl_local);
    cold_sr_fluid_flux_local_copy(qr_local, fr_local);

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
      printf("fr[%d]: %1.4e, fl[%d]: %1.4e, amdq[%d]: %1.4e, apdq[%d]: %1.4e ",i,fr[i],i,fl[i],i,amdq[i],i,apdq[i]);
      printf("  df[%d]: %1.4e, dadq[%d]: %1.4e,     total_diff: %1.4e\n",i,fr[i]-fl[i],i,amdq[i]+apdq[i], fr[i]-fl[i]-(amdq[i]+apdq[i]));
      TEST_CHECK( gkyl_compare(fr[i]-fl[i], amdq[i]+apdq[i], 1e-8) ); 
    }
  }

  gkyl_wv_eqn_release(sr_cold_fluid);
}


void test_flux_jump_relation() { test_flux_jump(); }
void test_sr_cold_fluids_internal() { test_sr_cold_fluid(); }

TEST_LIST = {
  {"test_flux_jump", test_flux_jump_relation},
  {"test_cold_sr_fluids", test_sr_cold_fluids_internal},
  {NULL, NULL},
};