/* -*- c -*- */

#include <stdio.h>
#include <gkylzero.h>
#include <gkyl_wv_ten_moment.h>
#include <gkyl_wv_ten_moment_priv.h>

extern "C" {
    int cu_wv_ten_moment_test(const struct gkyl_wv_eqn *eqn);
}

__global__
void ker_cu_wv_ten_moment_test(const struct gkyl_wv_eqn *eqn, int *nfail)
{
  *nfail = 0;

  GKYL_CU_CHECK( eqn->num_equations == 10, nfail );
  GKYL_CU_CHECK( eqn->num_waves == 5, nfail );

  // DO NOT DO THIS IN PRODUCTION! ONLY FOR TESTING
  struct wv_ten_moment *ten_moment = container_of(eqn, struct wv_ten_moment, eqn);

  GKYL_CU_CHECK( ten_moment->k0 == 1.0, nfail );

  double rho = 1.0, u = 0.1, v = 0.2, w = 0.3;
  double pxx = 0.5, pxy = 0.1, pxz = 0.2, pyy = 1.0, pyz = 0.3, pzz = 1.5;
  double q[10], pv[10] = { rho, u, v, w, pxx, pxy, pxz, pyy, pyz, pzz };

  q[0] = rho;
  q[1] = rho*u; q[2] = rho*v; q[3] = rho*w;
  q[4] = pxx + rho*u*u; q[5] = pxy + rho*u*v; q[6] = pxz + rho*u*w;
  q[7] = pyy + rho*v*v; q[8] = pyz + rho*v*w; q[9] = pzz + rho*w*w;

  // new array for re-computed primitive variables
  double var[10];
  gkyl_ten_moment_primitive(q, var);
  // check implementation of primitive variable calculation
  GKYL_CU_CHECK ( var[0] == pv[0], nfail );
  GKYL_CU_CHECK ( var[1] == pv[1], nfail );
  GKYL_CU_CHECK ( var[2] == pv[2], nfail );
  GKYL_CU_CHECK ( var[3] == pv[3], nfail );
  GKYL_CU_CHECK ( var[4] == pv[4], nfail );
  GKYL_CU_CHECK ( var[5] == pv[5], nfail );
  GKYL_CU_CHECK ( var[6] == pv[6], nfail );
  GKYL_CU_CHECK ( var[7] == pv[7], nfail );
  GKYL_CU_CHECK ( var[8] == pv[8], nfail );
  GKYL_CU_CHECK ( var[9] == pv[9], nfail );

  double fluxes[3][10] = {
    { rho*u, rho*u*u + pxx, rho*u*v + pxy, rho*u*w + pxz, 
      rho*u*u*u + 3*u*pxx, rho*u*u*v + 2*u*pxy + v*pxx, rho*u*u*w + 2*u*pxz + w*pxx,
      rho*u*v*v + 2*v*pxy + u*pyy, rho*u*v*w + u*pyz + v*pxz + w*pxy, rho*u*w*w + 2*w*pxz + u*pzz },
    { rho*v, rho*u*v + pxy, rho*v*v + pyy, rho*v*w + pyz, 
      rho*v*u*u + 2*u*pxy + v*pxx, rho*u*v*v + 2*v*pxy + u*pyy, rho*u*v*w + u*pyz + v*pxz + w*pxy,
      rho*v*v*v + 3*v*pyy, rho*v*v*w + 2*v*pyz + w*pyy, rho*v*w*w + 2*w*pyz + v*pzz},
    { rho*w, rho*u*w + pxz, rho*v*w + pyz, rho*w*w + pzz, 
      rho*u*u*w + 2*u*pxz + w*pxx, rho*u*v*w + u*pyz + v*pxz + w*pxy, rho*u*w*w + 2*w*pxz + u*pzz,
      rho*v*v*w + 2*v*pyz + w*pyy, rho*v*w*w + 2*w*pyz + v*pzz, rho*w*w*w + 3*w*pzz },
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

  double q_local[10], flux_local[10], flux[10];
  // finite precision error for d = 0
  // need better GKYL_CU_CHECK method, for now start at d=1
  for (int d=1; d<3; ++d) {
    eqn->rotate_to_local_func(tau1[d], tau2[d], norm[d], q, q_local);
    gkyl_ten_moment_flux(q_local, flux_local);
    eqn->rotate_to_global_func(tau1[d], tau2[d], norm[d], flux_local, flux);
    
    for (int m=0; m<10; ++m)
      GKYL_CU_CHECK( flux[m] == fluxes[d][m], nfail );
  }

  double q_l[10], q_g[10];
  for (int d=0; d<3; ++d) {
    eqn->rotate_to_local_func(tau1[d], tau2[d], norm[d], q, q_l);
    eqn->rotate_to_global_func(tau1[d], tau2[d], norm[d], q_l, q_g);

    for (int m=0; m<10; ++m) GKYL_CU_CHECK( q[m] == q_g[m] , nfail );
  }
}

int cu_wv_ten_moment_test(const struct gkyl_wv_eqn *eqn)
{
  int *nfail_dev = (int *) gkyl_cu_malloc(sizeof(int));
  ker_cu_wv_ten_moment_test<<<1,1>>>(eqn, nfail_dev);

  int nfail;
  gkyl_cu_memcpy(&nfail, nfail_dev, sizeof(int), GKYL_CU_MEMCPY_D2H);
  gkyl_cu_free(nfail_dev);

  return nfail;
}
