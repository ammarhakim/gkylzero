#include <acutest.h>
#include <gkyl_iso_euler_prim.h>
#include <gkyl_wv_iso_euler.h>

void
calcq(const double pv[4], double q[4])
{
  double rho = pv[0], u = pv[1], v = pv[2], w = pv[3];
  q[0] = rho;
  q[1] = rho*u; q[2] = rho*v; q[3] = rho*w;
}

void
test_iso_euler_basic()
{
  double vt = 1.0;
  struct gkyl_wv_eqn *iso_euler = gkyl_wv_iso_euler_new(vt);

  TEST_CHECK( iso_euler->num_equations == 4 );
  TEST_CHECK( iso_euler->num_waves == 3 );

  double rho = 1.0, u = 0.1, v = 0.2, w = 0.3;
  double q[4], pv[4] = { rho, u, v, w };
  calcq(pv, q);

  double flux[4];

  gkyl_iso_euler_flux(0, vt, q, flux);
  TEST_CHECK( flux[0] == rho*u );
  TEST_CHECK( flux[1] == rho*(u*u + vt*vt) );
  TEST_CHECK( flux[2] == rho*u*v );
  TEST_CHECK( flux[3] == rho*u*w );

  gkyl_iso_euler_flux(1, vt, q, flux);
  TEST_CHECK( flux[0] == rho*v );
  TEST_CHECK( flux[1] == rho*v*u );
  TEST_CHECK( flux[2] == rho*(v*v + vt*vt) );
  TEST_CHECK( flux[3] == rho*v*w );

  gkyl_iso_euler_flux(2, vt, q, flux);
  TEST_CHECK( flux[0] == rho*w );
  TEST_CHECK( flux[1] == rho*w*u );
  TEST_CHECK( flux[2] == rho*w*v );
  TEST_CHECK( flux[3] == rho*(w*w + vt*vt) );
  
  gkyl_wv_eqn_release(iso_euler);
}

void
test_iso_euler_waves()
{
  double vt = 1.0;
  struct gkyl_wv_eqn *iso_euler = gkyl_wv_iso_euler_new(vt);

  double vl[4] = { 1.0, 0.1, 0.2, 0.3};
  double vr[4] = { 0.1, 1.0, 2.0, 3.0};

  double ql[4], qr[4];
  calcq(vl, ql); calcq(vr, qr);

  double delta[4];
  for (int i=0; i<4; ++i) delta[i] = qr[i]-ql[i];

  for (int d=0; d<3; ++d) {
    double speeds[3], waves[3*4];
    gkyl_wv_eqn_waves(iso_euler, d, delta, ql, qr, waves, speeds);

    double apdq[4], amdq[4];
    gkyl_wv_eqn_qfluct(iso_euler, d, ql, qr, waves, speeds, amdq, apdq);
    
    // check if sum of left/right going fluctuations sum to jump in flux
    double fl[4], fr[4];
    gkyl_iso_euler_flux(d, vt, ql, fl);
    gkyl_iso_euler_flux(d, vt, qr, fr);
    
    for (int i=0; i<4; ++i)
      TEST_CHECK( gkyl_compare(fr[i]-fl[i], amdq[i]+apdq[i], 1e-14) );
  }
    
  gkyl_wv_eqn_release(iso_euler);
}

void
test_iso_euler_waves_2()
{
  double vt = 10.0;
  struct gkyl_wv_eqn *iso_euler = gkyl_wv_iso_euler_new(vt);

  double vl[4] = { 1.0, 0.1, 0.2, 0.3};
  double vr[4] = { 0.01, 1.0, 2.0, 3.0};

  double ql[4], qr[4];
  calcq(vl, ql); calcq(vr, qr);

  double delta[4];
  for (int i=0; i<4; ++i) delta[i] = qr[i]-ql[i];

  for (int d=0; d<3; ++d) {
    double speeds[3], waves[3*4];
    gkyl_wv_eqn_waves(iso_euler, d, delta, ql, qr, waves, speeds);

    double apdq[4], amdq[4];
    gkyl_wv_eqn_qfluct(iso_euler, d, ql, qr, waves, speeds, amdq, apdq);
    
    // check if sum of left/right going fluctuations sum to jump in flux
    double fl[4], fr[4];
    gkyl_iso_euler_flux(d, vt, ql, fl);
    gkyl_iso_euler_flux(d, vt, qr, fr);
    
    for (int i=0; i<4; ++i)
      TEST_CHECK( gkyl_compare(fr[i]-fl[i], amdq[i]+apdq[i], 1e-14) );
  }
  
  gkyl_wv_eqn_release(iso_euler);
}

TEST_LIST = {
  { "iso_euler_basic", test_iso_euler_basic },
  { "iso_euler_waves", test_iso_euler_waves },
  { "iso_euler_waves_2", test_iso_euler_waves_2 },
  { NULL, NULL },
};
