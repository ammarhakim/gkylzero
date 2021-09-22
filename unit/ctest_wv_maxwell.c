#include <acutest.h>
#include <gkyl_prim_maxwell.h>
#include <gkyl_wv_maxwell.h>

static const int dir_shuffle[][6] = {
  {0, 1, 2, 3, 4, 5},
  {1, 2, 0, 4, 5, 3},
  {2, 0, 1, 5, 3, 4}
};

// Make indexing cleaner with the dir_shuffle
#define EX d[0]
#define EY d[1]
#define EZ d[2]
#define BX d[3]
#define BY d[4]
#define BZ d[5]

void
test_maxwell_basic()
{
  // speed of light in SI units so tests are non-trivial
  double c = 299792458.0;
  double c2 = c*c;
  double e_fact = 1.0;
  double b_fact = 1.0;
  struct gkyl_wv_eqn *maxwell = gkyl_wv_maxwell_new(c, e_fact, b_fact);

  TEST_CHECK( maxwell->num_equations == 8 );
  TEST_CHECK( maxwell->num_waves == 6 );

  double Ex = 1.0, Ey = 0.1, Ez = 0.2;
  double Bx = 1.0, By = 0.1, Bz = 0.2;
  double phi = 0.01, psi = 0.01;
  double q[8] = { Ex, Ey, Ez, Bx, By, Bz, phi, psi };

  double flux[8];

  gkyl_maxwell_flux(0, c, e_fact, b_fact, q, flux);
  const int *d = dir_shuffle[0];
  TEST_CHECK( flux[EX] == e_fact*c2*q[6] );
  TEST_CHECK( flux[EY] == c2*q[BZ] );
  TEST_CHECK( flux[EZ] == -c2*q[BY] );
  TEST_CHECK( flux[BX] == b_fact*q[7] );
  TEST_CHECK( flux[BY] == -q[EZ] );
  TEST_CHECK( flux[BZ] == q[EY] );
  TEST_CHECK( flux[6] == e_fact*q[EX] );
  TEST_CHECK( flux[7] == b_fact*c2*q[BX] );

  gkyl_maxwell_flux(1, c, e_fact, b_fact, q, flux);
  d = dir_shuffle[1];
  TEST_CHECK( flux[EX] == e_fact*c2*q[6] );
  TEST_CHECK( flux[EY] == c2*q[BZ] );
  TEST_CHECK( flux[EZ] == -c2*q[BY] );
  TEST_CHECK( flux[BX] == b_fact*q[7] );
  TEST_CHECK( flux[BY] == -q[EZ] );
  TEST_CHECK( flux[BZ] == q[EY] );
  TEST_CHECK( flux[6] == e_fact*q[EX] );
  TEST_CHECK( flux[7] == b_fact*c2*q[BX] );

  gkyl_maxwell_flux(2, c, e_fact, b_fact, q, flux);
  d = dir_shuffle[2];
  TEST_CHECK( flux[EX] == e_fact*c2*q[6] );
  TEST_CHECK( flux[EY] == c2*q[BZ] );
  TEST_CHECK( flux[EZ] == -c2*q[BY] );
  TEST_CHECK( flux[BX] == b_fact*q[7] );
  TEST_CHECK( flux[BY] == -q[EZ] );
  TEST_CHECK( flux[BZ] == q[EY] );
  TEST_CHECK( flux[6] == e_fact*q[EX] );
  TEST_CHECK( flux[7] == b_fact*c2*q[BX] );
  
  gkyl_wv_eqn_release(maxwell);
}

void
test_maxwell_waves()
{
  // speed of light in SI units so tests are non-trivial
  double c = 299792458.0;
  double c2 = c*c;
  double e_fact = 1.0;
  double b_fact = 1.0;
  struct gkyl_wv_eqn *maxwell = gkyl_wv_maxwell_new(c, e_fact, b_fact);

  double ql[8] = { 0.0, 1.0, 0.0, 1.0, -0.75, 0.0, 0.0, 0.0};
  double qr[8] = { 0.0, -1.0, 0.0, 1.0, 0.75, 0.0, 0.0, 0.0};
  double ql_local[8] = { 0.0, 1.0, 0.0, 1.0, -0.75, 0.0, 0.0, 0.0};
  double qr_local[8] = { 0.0, -1.0, 0.0, 1.0, 0.75, 0.0, 0.0, 0.0};

  for (int d=0; d<3; ++d) {
    double speeds[6], waves[6*8], waves_local[6*8];
    // rotate to local tangent-normal frame
    gkyl_wv_eqn_rotate_to_local(maxwell, d, 0, 0, 0, ql, ql_local);
    gkyl_wv_eqn_rotate_to_local(maxwell, d, 0, 0, 0, qr, qr_local);

    double delta[8];
    for (int i=0; i<8; ++i) delta[i] = qr_local[i]-ql_local[i];
    
    gkyl_wv_eqn_waves(maxwell, d, delta, ql_local, qr_local, waves_local, speeds);

    // rotate waves back to global frame
    for (int mw=0; mw<6; ++mw)
      gkyl_wv_eqn_rotate_to_global(maxwell, d, 0, 0, 0, &waves_local[mw*8], &waves[mw*8]);

    double apdq[8], amdq[8];
    gkyl_wv_eqn_qfluct(maxwell, d, ql, qr, waves, speeds, amdq, apdq);
    
    // check if sum of left/right going fluctuations sum to jump in flux
    double fl[8], fr[8];
    gkyl_maxwell_flux(d, c, e_fact, b_fact, ql, fl);
    gkyl_maxwell_flux(d, c, e_fact, b_fact, qr, fr);
    
    for (int i=0; i<8; ++i)
      TEST_CHECK( gkyl_compare(fr[i]-fl[i], amdq[i]+apdq[i], 1e-14) );
  }
    
  gkyl_wv_eqn_release(maxwell);
}

TEST_LIST = {
  { "maxwell_basic", test_maxwell_basic },
  { "maxwell_waves", test_maxwell_waves },
  { NULL, NULL },
};
