#include <acutest.h>
#include <gkyl_wv_iso_euler.h>
#include <gkyl_wv_iso_euler_priv.h>

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
  double vt = 2.0;
  struct gkyl_wv_eqn *iso_euler = gkyl_wv_iso_euler_new(vt);

  TEST_CHECK( iso_euler->num_equations == 4 );
  TEST_CHECK( iso_euler->num_waves == 3 );

  double rho = 1.0, u = 0.1, v = 0.2, w = 0.3;
  double q[4], pv[4] = { rho, u, v, w };
  calcq(pv, q);

  double fluxes[3][4] = {
    { rho*u, rho*(u*u+vt*vt), rho*u*v, rho*u*w },
    { rho*v, rho*u*v, rho*(v*v+vt*vt), rho*v*w },
    { rho*w, rho*u*w, rho*v*w, rho*(w*w+vt*vt) },
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

  double q_local[4], flux_local[4], flux[4];
  for (int d=0; d<3; ++d) {
    iso_euler->rotate_to_local_func(tau1[d], tau2[d], norm[d], q, q_local);
    gkyl_iso_euler_flux(vt, q_local, flux_local);
    iso_euler->rotate_to_global_func(tau1[d], tau2[d], norm[d], flux_local, flux);
    
    for (int m=0; m<4; ++m)
      TEST_CHECK( gkyl_compare(flux[m], fluxes[d][m], 1e-14) );
  }  

  double q_l[4], q_g[4];
  for (int d=0; d<3; ++d) {
    iso_euler->rotate_to_local_func(tau1[d], tau2[d], norm[d], q, q_l);
    iso_euler->rotate_to_global_func(tau1[d], tau2[d], norm[d], q_l, q_g);

    for (int m=0; m<4; ++m) TEST_CHECK( q[m] == q_g[m] );
  }

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
  double ql_local[4], qr_local[4];
  calcq(vl, ql); calcq(vr, qr);

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

  for (int d=0; d<3; ++d) {
    double speeds[3], waves[3*4], waves_local[3*4];
    // rotate to local tangent-normal frame
    iso_euler->rotate_to_local_func(tau1[d], tau2[d], norm[d], ql, ql_local);
    iso_euler->rotate_to_local_func(tau1[d], tau2[d], norm[d], qr, qr_local);

    double delta[4];
    for (int i=0; i<4; ++i) delta[i] = qr_local[i]-ql_local[i];
    
    iso_euler->waves_func(iso_euler, delta, ql_local, qr_local, waves_local, speeds);

    // rotate waves back to global frame
    for (int mw=0; mw<3; ++mw)
      iso_euler->rotate_to_global_func(tau1[d], tau2[d], norm[d], &waves_local[mw*4], &waves[mw*4]);
    
    double apdq[4], amdq[4];
    iso_euler->qfluct_func(iso_euler, ql, qr, waves, speeds, amdq, apdq);
    
    // check if sum of left/right going fluctuations sum to jump in flux
    double fl_local[4], fr_local[4];
    gkyl_iso_euler_flux(vt, ql_local, fl_local);
    gkyl_iso_euler_flux(vt, qr_local, fr_local);
    
    double fl[4], fr[4];
    iso_euler->rotate_to_global_func(tau1[d], tau2[d], norm[d], fl_local, fl);
    iso_euler->rotate_to_global_func(tau1[d], tau2[d], norm[d], fr_local, fr);
    
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
  double ql_local[4], qr_local[4];
  calcq(vl, ql); calcq(vr, qr);

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

  for (int d=0; d<3; ++d) {
    double speeds[3], waves[3*4], waves_local[3*4];
    // rotate to local tangent-normal frame
    iso_euler->rotate_to_local_func(tau1[d], tau2[d], norm[d], ql, ql_local);
    iso_euler->rotate_to_local_func(tau1[d], tau2[d], norm[d], qr, qr_local);

    double delta[4];
    for (int i=0; i<4; ++i) delta[i] = qr_local[i]-ql_local[i];
    
    iso_euler->waves_func(iso_euler, delta, ql_local, qr_local, waves_local, speeds);

    // rotate waves back to global frame
    for (int mw=0; mw<3; ++mw)
      iso_euler->rotate_to_global_func(tau1[d], tau2[d], norm[d], &waves_local[mw*4], &waves[mw*4]);

    double apdq[4], amdq[4];
    iso_euler->qfluct_func(iso_euler, ql, qr, waves, speeds, amdq, apdq);
    
    // check if sum of left/right going fluctuations sum to jump in flux
    double fl_local[4], fr_local[4];
    gkyl_iso_euler_flux(vt, ql_local, fl_local);
    gkyl_iso_euler_flux(vt, qr_local, fr_local);
    
    double fl[4], fr[4];
    iso_euler->rotate_to_global_func(tau1[d], tau2[d], norm[d], fl_local, fl);
    iso_euler->rotate_to_global_func(tau1[d], tau2[d], norm[d], fr_local, fr);
    
    for (int i=0; i<4; ++i)
      TEST_CHECK( gkyl_compare(fr[i]-fl[i], amdq[i]+apdq[i], 1e-14) );
  }
  
  gkyl_wv_eqn_release(iso_euler);
}

#ifdef GKYL_HAVE_CUDA

int cu_wv_iso_euler_test(const struct gkyl_wv_eqn *eqn);

void
test_cu_wv_iso_euler()
{
  double vt = 2.0;
  struct gkyl_wv_eqn *eqn = gkyl_wv_iso_euler_cu_dev_new(vt);

  // this is not possible from user code and should NOT be done. This
  // is for testing only
  struct wv_iso_euler *iso_euler = container_of(eqn, struct wv_iso_euler, eqn);

  TEST_CHECK( iso_euler->vt == 2.0 ); 
  
  // call CUDA test
  int nfail = cu_wv_iso_euler_test(eqn->on_dev);

  TEST_CHECK( nfail == 0 );

  gkyl_wv_eqn_release(eqn);
}

#endif

TEST_LIST = {
  { "iso_euler_basic", test_iso_euler_basic },
  { "iso_euler_waves", test_iso_euler_waves },
  { "iso_euler_waves_2", test_iso_euler_waves_2 },
#ifdef GKYL_HAVE_CUDA
  { "cu_wv_iso_euler", test_cu_wv_iso_euler },
#endif  
  { NULL, NULL },
};
