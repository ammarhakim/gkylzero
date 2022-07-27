#include <acutest.h>
#include <gkyl_moment_prim_euler.h>
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

  //double rho = 1.0, u = 0.1, v = 0.2, w = 0.3, pr = 1.5;
  double rho = 1.0, u = 0.1, v = 0.2, w = 0.3, pr = 0.0;
  double q[5], pv[5] = { rho, u, v, w, pr };

  calcq(gas_gamma, pv, q);
  double E = q[4];

  double fluxes[3][5] = {
    { rho*u, rho*u*u+pr, rho*u*v, rho*u*w, (E+pr)*u },
    { rho*v, rho*u*v, rho*v*v+pr, rho*v*w, (E+pr)*v },
    { rho*w, rho*u*w, rho*v*w, rho*w*w, (E+pr)*w },
  };

  double norm[3][3] = {
    { 1.0, 0.0, 0.0 },
    { 0.0, 1.0, 0.0 },
    { 0.0, 0.0, 1.0 }
  };

  double tau1[3][3] = {
    { 0.0, 1.0, 0.0 },
    { 1.0, 0.0, 0.0 },
    { 1.0, 0.0, 0.0 }
  };

  double tau2[3][3] = {
    { 0.0, 0.0, 1.0 },
    { 0.0, 0.0, -1.0 },
    { 0.0, 1.0, 0.0 }
  };

  TEST_CHECK ( pr == gkyl_euler_pressure(gas_gamma, q) );

  double q_local[5], flux_local[5], flux[5];
  for (int d=1; d<2; ++d) {
    euler->rotate_to_local_func(tau1[d], tau2[d], norm[d], q, q_local);
    gkyl_euler_flux(gas_gamma, q_local, flux_local);
    euler->rotate_to_global_func(tau1[d], tau2[d], norm[d], flux_local, flux);
    
    for (int m=0; m<5; ++m)
      TEST_CHECK( gkyl_compare(flux[m], fluxes[d][m], 1e-15) );
  }

  double q_l[5], q_g[5];
  for (int d=1; d<3; ++d) {
    gkyl_wv_eqn_rotate_to_local(euler, tau1[d], tau2[d], norm[d], q, q_l);
    gkyl_wv_eqn_rotate_to_global(euler, tau1[d], tau2[d], norm[d], q_l, q_g);

    for (int m=0; m<5; ++m) TEST_CHECK( q[m] == q_g[m] );
  }
  
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
  double ql_local[5], qr_local[5];
  calcq(gas_gamma, vl, ql); calcq(gas_gamma, vr, qr);

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
    double speeds[3], waves[3*5], waves_local[3*5];
    // rotate to local tangent-normal frame
    gkyl_wv_eqn_rotate_to_local(euler, tau1[d], tau2[d], norm[d], ql, ql_local);
    gkyl_wv_eqn_rotate_to_local(euler, tau1[d], tau2[d], norm[d], qr, qr_local);

    double delta[5];
    for (int i=0; i<5; ++i) delta[i] = qr_local[i]-ql_local[i];
    
    gkyl_wv_eqn_waves(euler, GKYL_WV_HIGH_ORDER_FLUX, delta, ql_local, qr_local, waves_local, speeds);

    // rotate waves back to global frame
    for (int mw=0; mw<3; ++mw)
      gkyl_wv_eqn_rotate_to_global(euler, tau1[d], tau2[d], norm[d], &waves_local[mw*5], &waves[mw*5]);
    
    double apdq[5], amdq[5];
    gkyl_wv_eqn_qfluct(euler, GKYL_WV_HIGH_ORDER_FLUX, ql, qr, waves, speeds, amdq, apdq);
    
    // check if sum of left/right going fluctuations sum to jump in flux
    double fl_local[5], fr_local[5];
    gkyl_euler_flux(gas_gamma, ql_local, fl_local);
    gkyl_euler_flux(gas_gamma, qr_local, fr_local);

    double fl[5], fr[5];
    gkyl_wv_eqn_rotate_to_global(euler, tau1[d], tau2[d], norm[d], fl_local, fl);
    gkyl_wv_eqn_rotate_to_global(euler, tau1[d], tau2[d], norm[d], fr_local, fr);
    
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
  double ql_local[5], qr_local[5];
  calcq(gas_gamma, vl, ql); calcq(gas_gamma, vr, qr);

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
    double speeds[3], waves[3*5], waves_local[3*5];
    // rotate to local tangent-normal frame
    gkyl_wv_eqn_rotate_to_local(euler, tau1[d], tau2[d], norm[d], ql, ql_local);
    gkyl_wv_eqn_rotate_to_local(euler, tau1[d], tau2[d], norm[d], qr, qr_local);

    double delta[5];
    for (int i=0; i<5; ++i) delta[i] = qr_local[i]-ql_local[i];
    
    gkyl_wv_eqn_waves(euler, GKYL_WV_HIGH_ORDER_FLUX, delta, ql_local, qr_local, waves_local, speeds);

    // rotate waves back to global frame
    for (int mw=0; mw<3; ++mw)
      gkyl_wv_eqn_rotate_to_global(euler, tau1[d], tau2[d], norm[d], &waves_local[mw*5], &waves[mw*5]);

    double apdq[5], amdq[5];
    gkyl_wv_eqn_qfluct(euler, GKYL_WV_HIGH_ORDER_FLUX, ql, qr, waves, speeds, amdq, apdq);
    
    // check if sum of left/right going fluctuations sum to jump in flux
    double fl_local[5], fr_local[5];
    gkyl_euler_flux(gas_gamma, ql_local, fl_local);
    gkyl_euler_flux(gas_gamma, qr_local, fr_local);

    double fl[5], fr[5];
    gkyl_wv_eqn_rotate_to_global(euler, tau1[d], tau2[d], norm[d], fl_local, fl);
    gkyl_wv_eqn_rotate_to_global(euler, tau1[d], tau2[d], norm[d], fr_local, fr);
    
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
