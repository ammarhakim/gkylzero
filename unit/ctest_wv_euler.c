#include <acutest.h>
#include <gkyl_euler_prim.h>
#include <gkyl_wv_euler.h>

void
calcq(double gas_gamma, const double pv[5], double q[5])
{
  double rho = pv[0], u = pv[1], v = pv[2], w = pv[3], pr = pv[4];
  q[0] = rho;
  q[1] = rho*u; q[2] = rho*v; q[3] = rho*w;
  q[4] = pr/(gas_gamma-1) + 0.5*rho*(u*u+v*v+w*w);
}

void
test_euler_basic()
{
  double gas_gamma = 1.4;
  struct gkyl_wv_eqn *euler = gkyl_wv_euler_new(gas_gamma);

  TEST_CHECK( euler->num_equations == 5 );
  TEST_CHECK( euler->num_waves == 3 );

  double rho = 1.0, u = 0.1, v = 0.2, w = 0.3, pr = 1.5;
  double q[5], pv[5] = { rho, u, v, w, pr };
  calcq(gas_gamma, pv, q);

  TEST_CHECK ( pr == gkyl_euler_pressure(gas_gamma, q) );

  double E = q[4], flux[5];

  gkyl_euler_flux(0, gas_gamma, q, flux);
  TEST_CHECK( flux[0] == rho*u );
  TEST_CHECK( flux[1] == rho*u*u + pr );
  TEST_CHECK( flux[2] == rho*u*v );
  TEST_CHECK( flux[3] == rho*u*w );
  TEST_CHECK( flux[4] == (E+pr)*u );

  gkyl_euler_flux(1, gas_gamma, q, flux);
  TEST_CHECK( flux[0] == rho*v );
  TEST_CHECK( flux[1] == rho*v*u );
  TEST_CHECK( flux[2] == rho*v*v + pr );
  TEST_CHECK( flux[3] == rho*v*w );
  TEST_CHECK( flux[4] == (E+pr)*v );

  gkyl_euler_flux(2, gas_gamma, q, flux);
  TEST_CHECK( flux[0] == rho*w );
  TEST_CHECK( flux[1] == rho*w*u );
  TEST_CHECK( flux[2] == rho*w*v );
  TEST_CHECK( flux[3] == rho*w*w + pr );
  TEST_CHECK( flux[4] == (E+pr)*w );
  
  gkyl_wv_eqn_release(euler);
}

void
test_euler_waves()
{
  double gas_gamma = 1.4;
  struct gkyl_wv_eqn *euler = gkyl_wv_euler_new(gas_gamma);

  double vl[5] = { 1.0, 0.1, 0.2, 0.3, 1.5};
  double vr[5] = { 0.1, 1.0, 2.0, 3.0, 0.15};

  double ql[5], qr[5];
  calcq(gas_gamma, vl, ql); calcq(gas_gamma, vr, qr);

  double delta[5];
  for (int i=0; i<5; ++i) delta[i] = qr[i]-ql[i];

  for (int d=0; d<3; ++d) {
    double speeds[3], waves[3*5];
    gkyl_wv_eqn_waves(euler, d, delta, ql, qr, waves, speeds);

    double apdq[5], amdq[5];
    gkyl_wv_eqn_qfluct(euler, d, ql, qr, waves, speeds, amdq, apdq);
    
    // check if sum of left/right going fluctuations sum to jump in flux
    double fl[5], fr[5];
    gkyl_euler_flux(d, gas_gamma, ql, fl);
    gkyl_euler_flux(d, gas_gamma, qr, fr);
    
    for (int i=0; i<5; ++i)
      TEST_CHECK( gkyl_compare(fr[i]-fl[i], amdq[i]+apdq[i], 1e-14) );
  }
    
  gkyl_wv_eqn_release(euler);
}

void
test_euler_waves_2()
{
  double gas_gamma = 1.4;
  struct gkyl_wv_eqn *euler = gkyl_wv_euler_new(gas_gamma);

  double vl[5] = { 1.0, 0.1, 0.2, 0.3, 1.5};
  double vr[5] = { 0.01, 1.0, 2.0, 3.0, 15.0};

  double ql[5], qr[5];
  calcq(gas_gamma, vl, ql); calcq(gas_gamma, vr, qr);

  double delta[5];
  for (int i=0; i<5; ++i) delta[i] = qr[i]-ql[i];

  for (int d=0; d<3; ++d) {
    double speeds[3], waves[3*5];
    gkyl_wv_eqn_waves(euler, d, delta, ql, qr, waves, speeds);

    double apdq[5], amdq[5];
    gkyl_wv_eqn_qfluct(euler, d, ql, qr, waves, speeds, amdq, apdq);
    
    // check if sum of left/right going fluctuations sum to jump in flux
    double fl[5], fr[5];
    gkyl_euler_flux(d, gas_gamma, ql, fl);
    gkyl_euler_flux(d, gas_gamma, qr, fr);
    
    for (int i=0; i<5; ++i)
      TEST_CHECK( gkyl_compare(fr[i]-fl[i], amdq[i]+apdq[i], 1e-14) );
  }
  
  gkyl_wv_eqn_release(euler);
}

TEST_LIST = {
  { "euler_basic", test_euler_basic },
  { "euler_waves", test_euler_waves },
  { "euler_waves_2", test_euler_waves_2 },
  { NULL, NULL },
};
