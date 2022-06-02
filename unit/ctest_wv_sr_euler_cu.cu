/* -*- c -*- */

#include <stdio.h>
#include <gkylzero.h>
#include <gkyl_wv_sr_euler.h>
#include <gkyl_wv_sr_euler_priv.h>

extern "C" {
    int cu_wv_sr_euler_test(const struct gkyl_wv_eqn *eqn);
}

__global__
void ker_cu_wv_sr_euler_test(const struct gkyl_wv_eqn *eqn, int *nfail)
{
  *nfail = 0;

  GKYL_CU_CHECK( eqn->num_equations == 5, nfail );
  GKYL_CU_CHECK( eqn->num_waves == 3, nfail );

  // DO NOT DO THIS IN PRODUCTION! ONLY FOR TESTING
  struct wv_sr_euler *sr_euler = container_of(eqn, struct wv_sr_euler, eqn);

  GKYL_CU_CHECK( sr_euler->gas_gamma == 1.333, nfail );

  double rho = 1.0, u = 0.9999, v = 0.5, w = 0.01, pr = 1.5;
  double q[5], pv2[5];

  double gamma = 1 / sqrt(1 - u*u - v*v - w*w);
  double rhoh = gas_gamma * pr / (gas_gamma - 1)  + rho;
  
  q[0] = gamma*rho;
  q[1] = gamma*gamma*rhoh - pr;
  q[2] = gamma*gamma*rhoh*u;
  q[3] = gamma*gamma*rhoh*v;
  q[4] = gamma*gamma*rhoh*w;

  gkyl_sr_euler_prim_vars(gas_gamma, q, pv2);
  
  GKYL_CU_CHECK( rho == pv2[0], nfail );
  GKYL_CU_CHECK( pr == pv2[1], nfail );
  GKYL_CU_CHECK( u == pv2[2], nfail );
  GKYL_CU_CHECK( v == pv2[3], nfail );
  GKYL_CU_CHECK( w == pv2[4], nfail );

  double gamma = 1 / sqrt(1 - u*u - v*v - w*w);
  double rhoh = gas_gamma * pr / (gas_gamma - 1)  + rho;

  double fluxes[3][5] = {
    { gamma*rho*u, gamma*gamma*rhoh*u, gamma*gamma*rhoh*u*u + pr, gamma*gamma*rhoh*u*v, gamma*gamma*rhoh*u*w },
    { gamma*rho*v, gamma*gamma*rhoh*v, gamma*gamma*rhoh*v*u, gamma*gamma*rhoh*v*v + pr, gamma*gamma*rhoh*v*w },
    { gamma*rho*w, gamma*gamma*rhoh*w, gamma*gamma*rhoh*w*u, gamma*gamma*rhoh*w*v, gamma*gamma*rhoh*w*w + pr },
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

  double q_local[5], flux_local[5], flux[5];
  for (int d=0; d<3; ++d) {
    eqn->rotate_to_local_func(tau1[d], tau2[d], norm[d], q, q_local);
    gkyl_sr_euler_flux(gas_gamma, q_local, flux_local);
    eqn->rotate_to_global_func(tau1[d], tau2[d], norm[d], flux_local, flux);
    
    for (int m=0; m<5; ++m)
      GKYL_CU_CHECK( flux[m] == fluxes[d][m], nfail );
  }

  double q_l[5], q_g[5];
  for (int d=0; d<3; ++d) {
    eqn->rotate_to_local_func(tau1[d], tau2[d], norm[d], q, q_l);
    eqn->rotate_to_global_func(tau1[d], tau2[d], norm[d], q_l, q_g);

    for (int m=0; m<5; ++m) GKYL_CU_CHECK( q[m] == q_g[m] , nfail );
  }
}

int cu_wv_sr_euler_test(const struct gkyl_wv_eqn *eqn)
{
  int *nfail_dev = (int *) gkyl_cu_malloc(sizeof(int));
  ker_cu_wv_sr_euler_test<<<1,1>>>(eqn, nfail_dev);

  int nfail;
  gkyl_cu_memcpy(&nfail, nfail_dev, sizeof(int), GKYL_CU_MEMCPY_D2H);
  gkyl_cu_free(nfail_dev);

  return nfail;
}
