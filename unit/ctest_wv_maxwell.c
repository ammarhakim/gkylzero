#include <acutest.h>
#include <gkyl_wv_maxwell.h>
#include <gkyl_wv_maxwell_priv.h>

// Make indexing cleaner with the dir_shuffle
#define EX 0
#define EY 1
#define EZ 2
#define BX 3
#define BY 4
#define BZ 5

void
test_maxwell_basic()
{
  // speed of light in SI units so tests are non-trivial
  double c = 299792458.0;
  double c2 = c*c;
  double e_fact = 2.0;
  double b_fact = 2.5;
  struct gkyl_wv_eqn *maxwell = gkyl_wv_maxwell_new(c, e_fact, b_fact);

  TEST_CHECK( maxwell->num_equations == 8 );
  TEST_CHECK( maxwell->num_waves == 6 );

  double Ex = 1.0, Ey = 0.1, Ez = 0.2;
  double Bx = 10.0, By = 10.1, Bz = 10.2;
  double phi = 0.01, psi = 0.02;
  double q[8] = { Ex, Ey, Ez, Bx, By, Bz, phi, psi };

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

  double fluxes[3][8] = {
    {
      e_fact*c2*q[6],
      c2*q[BZ],
      -c2*q[BY],
      
      b_fact*q[7],
      -q[EZ],
      q[EY],
      
      e_fact*q[EX],
      b_fact*c2*q[BX]
    },
    {
      -c2*q[BZ],
      e_fact*c2*q[6],
      c2*q[BX],
      
      q[EZ],
      b_fact*q[7],
      -q[EX],
      
      e_fact*q[EY],
      b_fact*c2*q[BY]
    },
    {
      c2*q[BY],
      -c2*q[BX],      
      e_fact*c2*q[6],
      
      -q[EY],
      q[EX],
      b_fact*q[7],
      
      e_fact*q[EZ],
      b_fact*c2*q[BZ]
    },

  };


  double q_local[8], flux_local[8], flux[8];

  for (int d=0; d<3; ++d) {
    maxwell->rotate_to_local_func(tau1[d], tau2[d], norm[d], q, q_local);
    gkyl_maxwell_flux(c, e_fact, b_fact, q_local, flux_local);
    maxwell->rotate_to_global_func(tau1[d], tau2[d], norm[d], flux_local, flux);

    for (int m=0; m<8; ++m)
      TEST_CHECK( gkyl_compare(flux[m], fluxes[d][m], 1e-15) );
  }

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
    double speeds[6], waves[6*8], waves_local[6*8];
    double ql_local[8], qr_local[8];
    
    // rotate to local tangent-normal frame
    maxwell->rotate_to_local_func(tau1[d], tau2[d], norm[d], ql, ql_local);
    maxwell->rotate_to_local_func(tau1[d], tau2[d], norm[d], qr, qr_local);

    double delta[8];
    for (int i=0; i<8; ++i) delta[i] = qr_local[i]-ql_local[i];
    
    maxwell->waves_func(maxwell, delta, ql_local, qr_local, waves_local, speeds);

    // rotate waves back to global frame
    for (int mw=0; mw<6; ++mw)
      maxwell->rotate_to_global_func(tau1[d], tau2[d], norm[d], &waves_local[mw*8], &waves[mw*8]);

    double apdq[8], amdq[8];
    maxwell->qfluct_func(maxwell, ql, qr, waves, speeds, amdq, apdq);
    
    // check if sum of left/right going fluctuations sum to jump in flux
    double fl_local[8], fr_local[8];
    gkyl_maxwell_flux(c, e_fact, b_fact, ql_local, fl_local);
    gkyl_maxwell_flux(c, e_fact, b_fact, qr_local, fr_local);
    
    double fl[8], fr[8];
    maxwell->rotate_to_global_func(tau1[d], tau2[d], norm[d], fl_local, fl);
    maxwell->rotate_to_global_func(tau1[d], tau2[d], norm[d], fr_local, fr);
    
    for (int i=0; i<8; ++i)
      TEST_CHECK( gkyl_compare(fr[i]-fl[i], amdq[i]+apdq[i], 1e-14) );
  }
    
  gkyl_wv_eqn_release(maxwell);
}

#ifdef GKYL_HAVE_CUDA

int cu_maxwell_test(const struct gkyl_wv_eqn *eqn);

void
test_cu_wv_maxwell()
{
  // speed of light in SI units so tests are non-trivial
  double c = 299792458.0;
  double c2 = c*c;
  double e_fact = 2.0;
  double b_fact = 2.5;
  struct gkyl_wv_eqn *eqn = gkyl_wv_maxwell_cu_dev_new(c, e_fact, b_fact);

  // this is not possible from user code and should NOT be done. This
  // is for testing only
  struct wv_maxwell *maxwell = container_of(eqn, struct wv_maxwell, eqn);

  TEST_CHECK( maxwell->c == 299792458.0 );
  TEST_CHECK( maxwell->e_fact == 0.5 );
  TEST_CHECK( maxwell->b_fact == 0.25 );  
  
  // call CUDA test
  int nfail = cu_maxwell_test(eqn->on_dev);

  TEST_CHECK( nfail == 0 );

  gkyl_wv_eqn_release(eqn);
}

#endif

TEST_LIST = {
  { "maxwell_basic", test_maxwell_basic },
  { "maxwell_waves", test_maxwell_waves },
#ifdef GKYL_HAVE_CUDA
  { "cu_wv_maxwell", test_cu_wv_maxwell },
#endif  
  { NULL, NULL },
};
