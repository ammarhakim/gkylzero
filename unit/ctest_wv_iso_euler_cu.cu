/* -*- c -*- */

#include <stdio.h>
#include <gkylzero.h>
#include <gkyl_wv_iso_euler.h>
#include <gkyl_wv_iso_euler_priv.h>

extern "C" {
    int cu_wv_iso_euler_test(const struct gkyl_wv_eqn *eqn);
}

__global__
void ker_cu_wv_iso_euler_test(const struct gkyl_wv_eqn *eqn, int *nfail)
{
  *nfail = 0;

  GKYL_CU_CHECK( eqn->num_equations == 4, nfail );
  GKYL_CU_CHECK( eqn->num_waves == 3, nfail );

  // DO NOT DO THIS IN PRODUCTION! ONLY FOR TESTING
  struct wv_iso_euler *iso_euler = container_of(eqn, struct wv_iso_euler, eqn);

  GKYL_CU_CHECK( iso_euler->vt == 2.0, nfail );
  double vt = iso_euler->vt;
  
  double rho = 1.0, u = 0.1, v = 0.2, w = 0.3;
  double q[4];

  q[0] = rho;
  q[1] = rho*u; q[2] = rho*v; q[3] = rho*w;

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
    eqn->rotate_to_local_func(tau1[d], tau2[d], norm[d], q, q_local);
    gkyl_iso_euler_flux(vt, q_local, flux_local);
    eqn->rotate_to_global_func(tau1[d], tau2[d], norm[d], flux_local, flux);
    
    for (int m=0; m<4; ++m)
      GKYL_CU_CHECK( flux[m] == fluxes[d][m], nfail );
  }

  double q_l[4], q_g[4];
  for (int d=0; d<3; ++d) {
    eqn->rotate_to_local_func(tau1[d], tau2[d], norm[d], q, q_l);
    eqn->rotate_to_global_func(tau1[d], tau2[d], norm[d], q_l, q_g);

    for (int m=0; m<4; ++m) GKYL_CU_CHECK( q[m] == q_g[m] , nfail );
  }
}

int cu_wv_iso_euler_test(const struct gkyl_wv_eqn *eqn)
{
  int *nfail_dev = (int *) gkyl_cu_malloc(sizeof(int));
  ker_cu_wv_iso_euler_test<<<1,1>>>(eqn, nfail_dev);

  int nfail;
  gkyl_cu_memcpy(&nfail, nfail_dev, sizeof(int), GKYL_CU_MEMCPY_D2H);
  gkyl_cu_free(nfail_dev);

  return nfail;
}
