#include <acutest.h>
#include <gkyl_prim_sr_euler.h>
#include <gkyl_wv_sr_euler.h>

void
calcq(double gas_gamma, const double pv[5], double q[5])
{
  double rho = pv[0], u = pv[1], v = pv[2], w = pv[3], pr = pv[4];
  double gamma = 1 / sqrt(1 - u*u - v*v - w*w);
  double rhoh = gas_gamma * pr / (gas_gamma - 1)  + rho;
  
  q[0] = gamma*rho;
  q[1] = gamma*gamma*rhoh - pr;
  q[2] = gamma*gamma*rhoh*u;
  q[3] = gamma*gamma*rhoh*v;
  q[4] = gamma*gamma*rhoh*w;
  
}

void
test_sr_euler_prim1()
{
  double gas_gamma = 1.333;
  struct gkyl_wv_eqn *sr_euler = gkyl_wv_sr_euler_new(gas_gamma);

  TEST_CHECK( sr_euler->num_equations == 5 );
  TEST_CHECK( sr_euler->num_waves == 3 );

  double rho = 1.0, u = 0.9999, v = 0.5, w = 0.01, pr = 1.5;
  double q[5], pv2[5], pv[5] = { rho, u, v, w, pr };
  calcq(gas_gamma, pv, q);
  
  gkyl_sr_euler_prim_vars(gas_gamma, q, pv2);
  
  TEST_CHECK ( gkyl_compare(rho, pv2[0], 1e-15) );
  TEST_CHECK ( gkyl_compare(pr, pv2[1], 1e-15) );
  TEST_CHECK ( gkyl_compare(u, pv2[2], 1e-15) );
  TEST_CHECK ( gkyl_compare(v, pv2[3], 1e-15) );
  TEST_CHECK ( gkyl_compare(w, pv2[4], 1e-15) );

  double flux[5];
  double gamma = 1 / sqrt(1 - u*u - v*v - w*w);
  double rhoh = gas_gamma * pr / (gas_gamma - 1)  + rho;

  gkyl_sr_euler_flux(0, gas_gamma, q, flux);
  TEST_CHECK( gkyl_compare(flux[0], gamma*rho*u, 1e-15) );
  TEST_CHECK( gkyl_compare(flux[1], gamma*gamma*rhoh*u, 1e-15) );
  TEST_CHECK( gkyl_compare(flux[2], gamma*gamma*rhoh*u*u + pr ,1e-15));
  TEST_CHECK( gkyl_compare(flux[3], gamma*gamma*rhoh*u*v ,1e-15) );
  TEST_CHECK( gkyl_compare(flux[4], gamma*gamma*rhoh*u*w, 1e-15) );

  gkyl_sr_euler_flux(1, gas_gamma, q, flux);
  TEST_CHECK( gkyl_compare(flux[0], gamma*rho*v, 1e-15) );
  TEST_CHECK( gkyl_compare(flux[1], gamma*gamma*rhoh*v, 1e-15) );
  TEST_CHECK( gkyl_compare(flux[2], gamma*gamma*rhoh*v*u ,1e-15));
  TEST_CHECK( gkyl_compare(flux[3], gamma*gamma*rhoh*v*v + pr ,1e-15) );
  TEST_CHECK( gkyl_compare(flux[4], gamma*gamma*rhoh*v*w, 1e-15) );

  gkyl_sr_euler_flux(2, gas_gamma, q, flux);
  TEST_CHECK( gkyl_compare(flux[0], gamma*rho*w, 1e-15) );
  TEST_CHECK( gkyl_compare(flux[1], gamma*gamma*rhoh*w, 1e-15) );
  TEST_CHECK( gkyl_compare(flux[2], gamma*gamma*rhoh*w*u ,1e-15));
  TEST_CHECK( gkyl_compare(flux[3], gamma*gamma*rhoh*w*v ,1e-15) );
  TEST_CHECK( gkyl_compare(flux[4], gamma*gamma*rhoh*w*w + pr, 1e-15) );
    
  gkyl_wv_eqn_release(sr_euler);

  
}

void
test_sr_euler_waves()
{
  double gas_gamma = 1.333;
  struct gkyl_wv_eqn *sr_euler = gkyl_wv_sr_euler_new(gas_gamma);

  double vl[5] = { 1.0, 0.09, 0.02, 0.03, 1.5};
  double vr[5] = { 0.1, .9, 0.2, 0.3, 0.15};

  double ql[5], qr[5];
  double ql_local[5], qr_local[5];
  calcq(gas_gamma, vl, ql); calcq(gas_gamma, vr, qr);

  for (int d=0; d<3; ++d) {
    double speeds[3], waves[3*5], waves_local[3*5];
    // rotate to local tangent-normal frame
    gkyl_wv_eqn_rotate_to_local(sr_euler, d, 0, 0, 0, ql, ql_local);
    gkyl_wv_eqn_rotate_to_local(sr_euler, d, 0, 0, 0, qr, qr_local);

    double delta[5];
    for (int i=0; i<5; ++i) delta[i] = qr_local[i]-ql_local[i];
    
    gkyl_wv_eqn_waves(sr_euler, d, delta, ql_local, qr_local, waves_local, speeds);

    // rotate waves back to global frame
    for (int mw=0; mw<3; ++mw)
      gkyl_wv_eqn_rotate_to_global(sr_euler, d, 0, 0, 0, &waves_local[mw*5], &waves[mw*5]);
    
    double apdq[5], amdq[5];
    gkyl_wv_eqn_qfluct(sr_euler, d, ql, qr, waves, speeds, amdq, apdq);
    
    // check if sum of left/right going fluctuations sum to jump in flux
    double fl[5], fr[5];
    gkyl_sr_euler_flux(d, gas_gamma, ql, fl);
    gkyl_sr_euler_flux(d, gas_gamma, qr, fr);
    
    for (int i=0; i<5; ++i)
      TEST_CHECK( gkyl_compare(fr[i]-fl[i], amdq[i]+apdq[i], 1e-14) );
  }
    
  gkyl_wv_eqn_release(sr_euler);
}

void
test_sr_euler_waves2()
{
  double gas_gamma = 1.3333;
  struct gkyl_wv_eqn *sr_euler = gkyl_wv_sr_euler_new(gas_gamma);

  double vl[5] = { 1.0, 0.0999, 0.02, 0.03, 1500.};
  double vr[5] = { 0.1, 0.999, 0.2, 0.3, 0.015};
  
  double ql[5], qr[5];
  double ql_local[5], qr_local[5];
  calcq(gas_gamma, vl, ql); calcq(gas_gamma, vr, qr);

  for (int d=0; d<3; ++d) {
    double speeds[3], waves[3*5], waves_local[3*5];
    // rotate to local tangent-normal frame
    gkyl_wv_eqn_rotate_to_local(sr_euler, d, 0, 0, 0, ql, ql_local);
    gkyl_wv_eqn_rotate_to_local(sr_euler, d, 0, 0, 0, qr, qr_local);

    double delta[5];
    for (int i=0; i<5; ++i) delta[i] = qr_local[i]-ql_local[i];
    
    gkyl_wv_eqn_waves(sr_euler, d, delta, ql_local, qr_local, waves_local, speeds);

    // rotate waves back to global frame
    for (int mw=0; mw<3; ++mw)
      gkyl_wv_eqn_rotate_to_global(sr_euler, d, 0, 0, 0, &waves_local[mw*5], &waves[mw*5]);

    double apdq[5], amdq[5];
    gkyl_wv_eqn_qfluct(sr_euler, d, ql, qr, waves, speeds, amdq, apdq);
    
    // check if sum of left/right going fluctuations sum to jump in flux
    double fl[5], fr[5];
    gkyl_sr_euler_flux(d, gas_gamma, ql, fl);
    gkyl_sr_euler_flux(d, gas_gamma, qr, fr);
    
    for (int i=0; i<5; ++i)
      TEST_CHECK( gkyl_compare(fr[i]-fl[i], amdq[i]+apdq[i], 1e-14) );
  }

  gkyl_wv_eqn_release(sr_euler);
}


TEST_LIST = {
  { "euler_sr_prim1", test_sr_euler_prim1 },
  { "test_sr_euler_waves", test_sr_euler_waves},
  { "test_sr_euler_waves2", test_sr_euler_waves2},
  { NULL, NULL },
};
