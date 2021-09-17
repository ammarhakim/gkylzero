#include <acutest.h>
#include <gkyl_prim_mhd.h>
#include <gkyl_wv_mhd.h>

void
calcq(double gas_gamma, const double pv[8], double q[8])
{
  double rho = pv[0], u = pv[1], v = pv[2], w = pv[3], pr = pv[4];
  q[0] = rho;
  q[1] = rho*u; q[2] = rho*v; q[3] = rho*w;
  q[5] = pv[5]; q[6] = pv[6]; q[7] = pv[7];  // B field
  double pb = 0.5*(pv[5]*pv[5]+pv[6]*pv[6]+pv[7]*pv[7]); // magnetic pressure
  q[4] = pr/(gas_gamma-1) + 0.5*rho*(u*u+v*v+w*w) + pb;
}

/* Check flux function implementation */
void
test_mhd_basic()
{
  double gas_gamma = 1.4;
  int has_eight_waves = 0;
  struct gkyl_wv_eqn *mhd = gkyl_wv_mhd_new(gas_gamma, has_eight_waves);

  TEST_CHECK( mhd->num_equations == 8 );
  TEST_CHECK( mhd->num_waves == 7 );

  double rho = 1.0, u = 0.1, v = 0.2, w = 0.3, pr = 1.5;
  double bx = 0.11, by = 0.22, bz = 0.33;
  double q[8], pv[8] = { rho, u, v, w, pr, bx, by, bz };
  calcq(gas_gamma, pv, q);
  double pb = 0.5*(bx*bx+by*by+bz*bz);
  double u_dot_b = u*bx+v*by+w*bz;
  double E = q[4];

  TEST_CHECK ( pr == gkyl_mhd_pressure(gas_gamma, q) );

  double flux[8];
 
  gkyl_mhd_flux(0, gas_gamma, q, flux);
  TEST_CHECK( flux[0] == rho*u );
  TEST_CHECK( gkyl_compare(flux[1], rho*u*u - bx*bx + pr + pb, 1e-15) );
  TEST_CHECK( flux[2] == rho*u*v - bx*by );
  TEST_CHECK( flux[3] == rho*u*w - bx*bz );
  TEST_CHECK( flux[4] == (E+pr+pb)*u-bx*u_dot_b);
  TEST_CHECK( flux[5] == 0.0 );
  TEST_CHECK( flux[6] == u*by - v*bx );
  TEST_CHECK( flux[7] == u*bz - w*bx );

  gkyl_mhd_flux(1, gas_gamma, q, flux);
  TEST_CHECK( flux[0] == rho*v );
  TEST_CHECK( flux[1] == rho*v*u - bx*by );
  TEST_CHECK( gkyl_compare(flux[2], rho*v*v - by*by + pr +pb, 1e-15) );
  TEST_CHECK( flux[3] == rho*v*w - by*bz);
  TEST_CHECK( flux[4] == (E+pr+pb)*v-by*u_dot_b );
  TEST_CHECK( flux[5] == v*bx - u*by );
  TEST_CHECK( flux[6] == 0.0 );
  TEST_CHECK( flux[7] == v*bz - w*by );

  gkyl_mhd_flux(2, gas_gamma, q, flux);
  TEST_CHECK( flux[0] == rho*w );
  TEST_CHECK( flux[1] == rho*w*u - bx*bz );
  TEST_CHECK( flux[2] == rho*w*v - by*bz);
  TEST_CHECK( gkyl_compare(flux[3], rho*w*w - bz*bz + pr + pb, 1e-15) );
  TEST_CHECK( flux[4] == (E+pr+pb)*w-bz*u_dot_b );
  TEST_CHECK( flux[5] == w*bx - u*bz );
  TEST_CHECK( flux[6] == w*by - v*bz );
  TEST_CHECK( flux[6] == 0.0 );
  
  gkyl_wv_eqn_release(mhd);
}

/* check if sum of left/right going fluctuations sum to jump in flux */
void
test_mhd_waves()
{
  double gas_gamma = 1.4;
  int has_eight_waves = 0;
  struct gkyl_wv_eqn *mhd = gkyl_wv_mhd_new(gas_gamma, has_eight_waves);

  double ql[8], qr[8];
  double delta[8];

  // Bx must be the same since presently the wv_mhd does not have the divB wave
  double vl[8] = { 1.0,  0.1,  0.2,  0.3,  1.5, 0.4, 0.4,  0.3};
  double vr[8] = { 1.1, 0.13, 0.25, 0.34, 1.54, 0.4, 0.44, 0.34};

  calcq(gas_gamma, vl, ql); calcq(gas_gamma, vr, qr);

  for (int i=0; i<8; ++i) delta[i] = qr[i]-ql[i];

  for (int d=0; d<1; ++d) {
    double speeds[7], waves[7*8];
    gkyl_wv_eqn_waves(mhd, d, delta, ql, qr, waves, speeds);

    double apdq[8], amdq[8];
    gkyl_wv_eqn_qfluct(mhd, d, ql, qr, waves, speeds, amdq, apdq);
    
    double fl[8], fr[8];
    gkyl_mhd_flux(d, gas_gamma, ql, fl);
    gkyl_mhd_flux(d, gas_gamma, qr, fr);

    for (int i=0; i<8; ++i)
      TEST_CHECK( gkyl_compare(fr[i]-fl[i], amdq[i]+apdq[i], 1e-8) );
  }
    
  gkyl_wv_eqn_release(mhd);
}

TEST_LIST = {
  { "mhd_basic", test_mhd_basic },
  { "mhd_waves", test_mhd_waves },
  { NULL, NULL },
};
