/* -*- c -*- */

#include <stdio.h>
#include <gkylzero.h>
#include <gkyl_wv_mhd.h>
#include <gkyl_wv_mhd_priv.h>

extern "C" {
    int cu_wv_mhd_test(const struct gkyl_wv_eqn *eqn);
}

// FIXME: duplicate
GKYL_CU_D static void
calcq(double gas_gamma, const double pv[8], double q[8])
{
  double rho = pv[0], u = pv[1], v = pv[2], w = pv[3], pr = pv[4];
  q[0] = rho;
  q[1] = rho*u; q[2] = rho*v; q[3] = rho*w;
  q[5] = pv[5]; q[6] = pv[6]; q[7] = pv[7];  // B field
  double pb = 0.5*(pv[5]*pv[5]+pv[6]*pv[6]+pv[7]*pv[7]); // magnetic pressure
  q[4] = pr/(gas_gamma-1) + 0.5*rho*(u*u+v*v+w*w) + pb;
}

__global__
void ker_cu_wv_mhd_test(const struct gkyl_wv_eqn *eqn, int *nfail)
{
  *nfail = 0;

  double gas_gamma = 5.0/3.0;
  GKYL_CU_CHECK( eqn->num_equations == 8, nfail );
  GKYL_CU_CHECK( eqn->num_waves == 7, nfail );

  // DO NOT DO THIS IN PRODUCTION! ONLY FOR TESTING
  struct wv_mhd *mhd = container_of(eqn, struct wv_mhd, eqn);

  GKYL_CU_CHECK( mhd->gas_gamma == gas_gamma, nfail );

  double ql[8], qr[8];
  double ql_local[8], qr_local[8];

  // Bx must be the same since presently the wv_mhd does not have the divB wave
  double vl[8] = { 1.0,  0.1,  0.2,  0.3,  1.5, 0.4, 0.4,  0.3};
  double vr[8] = { 1.1, 0.13, 0.25, 0.34, 1.54, 0.4, 0.44, 0.34};

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

  for (int d=0; d<1; ++d) {
    double speeds[7], waves[7*8], waves_local[7*8];
    // rotate to local tangent-normal frame
    eqn->rotate_to_local_func(tau1[d], tau2[d], norm[d], ql, ql_local);
    eqn->rotate_to_local_func(tau1[d], tau2[d], norm[d], qr, qr_local);

    double delta[8];
    for (int i=0; i<8; ++i) delta[i] = qr_local[i]-ql_local[i];
    
    eqn->waves_func(eqn, delta, ql_local, qr_local, waves_local, speeds);

    // rotate waves back to global frame
    for (int mw=0; mw<7; ++mw)
      eqn->rotate_to_global_func(tau1[d], tau2[d], norm[d], &waves_local[mw*8], &waves[mw*8]);
    
    double apdq[8], amdq[8];
    eqn->qfluct_func(eqn, ql, qr, waves, speeds, amdq, apdq);
    
    // check if sum of left/right going fluctuations sum to jump in flux
    double fl_local[8], fr_local[8];
    gkyl_mhd_flux(gas_gamma, ql_local, fl_local);
    gkyl_mhd_flux(gas_gamma, qr_local, fr_local);

    double fl[8], fr[8];
    eqn->rotate_to_global_func(tau1[d], tau2[d], norm[d], fl_local, fl);
    eqn->rotate_to_global_func(tau1[d], tau2[d], norm[d], fr_local, fr);
    
    for (int i=0; i<8; ++i)
      GKYL_CU_CHECK( fabs((fr[i]-fl[i]) - (amdq[i]+apdq[i])) < 1e-13, nfail );
  }
}

int cu_wv_mhd_test(const struct gkyl_wv_eqn *eqn)
{
  int *nfail_dev = (int *) gkyl_cu_malloc(sizeof(int));
  ker_cu_wv_mhd_test<<<1,1>>>(eqn, nfail_dev);

  int nfail;
  gkyl_cu_memcpy(&nfail, nfail_dev, sizeof(int), GKYL_CU_MEMCPY_D2H);
  gkyl_cu_free(nfail_dev);

  return nfail;
}

