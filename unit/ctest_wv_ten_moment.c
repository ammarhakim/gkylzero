#include <acutest.h>
#include <gkyl_prim_ten_moment.h>
#include <gkyl_wv_ten_moment.h>

static const int dir_u_shuffle[][3] = {
  {1, 2, 3},
  {2, 3, 1},
  {3, 1, 2}
};

static const int dir_p_shuffle[][6] = {
  {4, 5, 6, 7, 8, 9},
  {7, 8, 5, 9, 6, 4},
  {9, 6, 8, 4, 5, 7}
};

// Make indexing cleaner with the dir_shuffle
#define RHOU d[0]
#define RHOV d[1]
#define RHOW d[2]

#define PXX dp[0]
#define PXY dp[1]
#define PXZ dp[2]
#define PYY dp[3]
#define PYZ dp[4]
#define PZZ dp[5]

void
calcq(const double pv[10], double q[10])
{
  double rho = pv[0], u = pv[1], v = pv[2], w = pv[3];
  double pxx = pv[4], pxy = pv[5], pxz = pv[6], pyy = pv[7], pyz = pv[8], pzz = pv[9];

  q[0] = rho;
  q[1] = rho*u; q[2] = rho*v; q[3] = rho*w;
  q[4] = pxx + rho*u*u; q[5] = pxy + rho*u*v; q[6] = pxz + rho*u*w;
  q[7] = pyy + rho*v*v; q[8] = pyz + rho*v*w; q[9] = pzz + rho*w*w;
}

void
test_ten_moment_basic()
{
  struct gkyl_wv_eqn *ten_moment = gkyl_wv_ten_moment_new();

  TEST_CHECK( ten_moment->num_equations == 10 );
  TEST_CHECK( ten_moment->num_waves == 5 );

  double rho = 1.0, u = 0.1, v = 0.2, w = 0.3;
  double pxx = 0.5, pxy = 0.1, pxz = 0.2, pyy = 1.0, pyz = 0.3, pzz = 1.5;
  double q[10], pv[10] = { rho, u, v, w, pxx, pxy, pxz, pyy, pyz, pzz };
  calcq(pv, q);

  // new array for re-computed primitive variables
  double var[10];
  gkyl_ten_moment_primitive(q, var);
  // check implementation of primitive variable calculation
  TEST_CHECK ( var[0] == pv[0] );
  TEST_CHECK ( var[1] == pv[1] );
  TEST_CHECK ( var[2] == pv[2] );
  TEST_CHECK ( var[3] == pv[3] );
  TEST_CHECK ( var[4] == pv[4] );
  TEST_CHECK ( var[5] == pv[5] );
  TEST_CHECK ( var[6] == pv[6] );
  TEST_CHECK ( var[7] == pv[7] );
  TEST_CHECK ( var[8] == pv[8] );
  TEST_CHECK ( var[9] == pv[9] );

  double flux[10];

  gkyl_ten_moment_flux(0, q, flux);
  const int *d = dir_u_shuffle[0];
  const int *dp = dir_p_shuffle[0];
  TEST_CHECK( flux[0] == q[RHOU] );
  TEST_CHECK( flux[RHOU] == q[PXX] );
  TEST_CHECK( flux[RHOV] == q[PXY] );
  TEST_CHECK( flux[RHOW] == q[PXZ] );
  TEST_CHECK( flux[PXX] == var[0]*var[RHOU]*var[RHOU]*var[RHOU] + 3*var[RHOU]*var[PXX] );
  TEST_CHECK( flux[PXY] == var[0]*var[RHOU]*var[RHOU]*var[RHOV] + 2*var[RHOU]*var[PXY] + var[RHOV]*var[PXX] );
  TEST_CHECK( flux[PXZ] == var[0]*var[RHOU]*var[RHOU]*var[RHOW] + 2*var[RHOU]*var[PXZ] + var[RHOW]*var[PXX] );
  TEST_CHECK( flux[PYY] == var[0]*var[RHOU]*var[RHOV]*var[RHOV] + 2*var[RHOV]*var[PXY] + var[RHOU]*var[PYY] );
  TEST_CHECK( flux[PYZ] == var[0]*var[RHOU]*var[RHOV]*var[RHOW] + var[RHOU]*var[PYZ] + var[RHOV]*var[PXZ] + var[RHOW]*var[PXY] );
  TEST_CHECK( flux[PZZ] == var[0]*var[RHOU]*var[RHOW]*var[RHOW] + 2*var[RHOW]*var[PXZ] + var[RHOU]*var[PZZ] );

  gkyl_ten_moment_flux(1, q, flux);
  d = dir_u_shuffle[1];
  dp = dir_p_shuffle[1];
  TEST_CHECK( flux[0] == q[RHOU] );
  TEST_CHECK( flux[RHOU] == q[PXX] );
  TEST_CHECK( flux[RHOV] == q[PXY] );
  TEST_CHECK( flux[RHOW] == q[PXZ] );
  TEST_CHECK( flux[PXX] == var[0]*var[RHOU]*var[RHOU]*var[RHOU] + 3*var[RHOU]*var[PXX] );
  TEST_CHECK( flux[PXY] == var[0]*var[RHOU]*var[RHOU]*var[RHOV] + 2*var[RHOU]*var[PXY] + var[RHOV]*var[PXX] );
  TEST_CHECK( flux[PXZ] == var[0]*var[RHOU]*var[RHOU]*var[RHOW] + 2*var[RHOU]*var[PXZ] + var[RHOW]*var[PXX] );
  TEST_CHECK( flux[PYY] == var[0]*var[RHOU]*var[RHOV]*var[RHOV] + 2*var[RHOV]*var[PXY] + var[RHOU]*var[PYY] );
  TEST_CHECK( flux[PYZ] == var[0]*var[RHOU]*var[RHOV]*var[RHOW] + var[RHOU]*var[PYZ] + var[RHOV]*var[PXZ] + var[RHOW]*var[PXY] );
  TEST_CHECK( flux[PZZ] == var[0]*var[RHOU]*var[RHOW]*var[RHOW] + 2*var[RHOW]*var[PXZ] + var[RHOU]*var[PZZ] );

  gkyl_ten_moment_flux(2, q, flux);
  d = dir_u_shuffle[2];
  dp = dir_p_shuffle[2];
  TEST_CHECK( flux[0] == q[RHOU] );
  TEST_CHECK( flux[RHOU] == q[PXX] );
  TEST_CHECK( flux[RHOV] == q[PXY] );
  TEST_CHECK( flux[RHOW] == q[PXZ] );
  TEST_CHECK( flux[PXX] == var[0]*var[RHOU]*var[RHOU]*var[RHOU] + 3*var[RHOU]*var[PXX] );
  TEST_CHECK( flux[PXY] == var[0]*var[RHOU]*var[RHOU]*var[RHOV] + 2*var[RHOU]*var[PXY] + var[RHOV]*var[PXX] );
  TEST_CHECK( flux[PXZ] == var[0]*var[RHOU]*var[RHOU]*var[RHOW] + 2*var[RHOU]*var[PXZ] + var[RHOW]*var[PXX] );
  TEST_CHECK( flux[PYY] == var[0]*var[RHOU]*var[RHOV]*var[RHOV] + 2*var[RHOV]*var[PXY] + var[RHOU]*var[PYY] );
  TEST_CHECK( flux[PYZ] == var[0]*var[RHOU]*var[RHOV]*var[RHOW] + var[RHOU]*var[PYZ] + var[RHOV]*var[PXZ] + var[RHOW]*var[PXY] );
  TEST_CHECK( flux[PZZ] == var[0]*var[RHOU]*var[RHOW]*var[RHOW] + 2*var[RHOW]*var[PXZ] + var[RHOU]*var[PZZ] );
  
  gkyl_wv_eqn_release(ten_moment);
}

void
test_ten_moment_waves()
{
  struct gkyl_wv_eqn *ten_moment = gkyl_wv_ten_moment_new();

  double vl[10] = { 1.0, 0.1, 0.2, 0.3, 0.5, 0.0, 0.0, 1.0, 0.0, 1.5};
  double vr[10] = { 0.1, 1.0, 2.0, 3.0, 0.1, 0.0, 0.0, 0.2, 0.0, 0.3};

  double ql[10], qr[10];
  calcq(vl, ql); calcq(vr, qr);

  double delta[10];
  for (int i=0; i<10; ++i) delta[i] = qr[i]-ql[i];

  for (int d=0; d<3; ++d) {
    double speeds[5], waves[5*10];
    gkyl_wv_eqn_waves(ten_moment, d, delta, ql, qr, waves, speeds);

    double apdq[10], amdq[10];
    gkyl_wv_eqn_qfluct(ten_moment, d, ql, qr, waves, speeds, amdq, apdq);
    
    // check if sum of left/right going fluctuations sum to jump in flux
    double fl[10], fr[10];
    gkyl_ten_moment_flux(d, ql, fl);
    gkyl_ten_moment_flux(d, qr, fr);
    
    for (int i=0; i<10; ++i)
      TEST_CHECK( gkyl_compare(fr[i]-fl[i], amdq[i]+apdq[i], 1e-14) );
  }
    
  gkyl_wv_eqn_release(ten_moment);
}

TEST_LIST = {
  { "ten_moment_basic", test_ten_moment_basic },
  { "ten_moment_waves", test_ten_moment_waves },
  { NULL, NULL },
};
