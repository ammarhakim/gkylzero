/* -*- c -*- */

#include <stdio.h>
#include <gkylzero.h>
#include <gkyl_wv_maxwell.h>
#include <gkyl_wv_maxwell_priv.h>

// Make indexing cleaner with the dir_shuffle
#define EX 0
#define EY 1
#define EZ 2
#define BX 3
#define BY 4
#define BZ 5

extern "C" {
    int cu_maxwell_test(const struct gkyl_wv_eqn *eqn);
}

__global__
void ker_cu_maxwell_test(const struct gkyl_wv_eqn *eqn, int *nfail)
{
  *nfail = 0;

  GKYL_CU_CHECK( eqn->num_equations == 8, nfail );
  GKYL_CU_CHECK( eqn->num_waves == 6 );

  // DO NOT DO THIS IN PRODUCTION! ONLY FOR TESTING
  struct wv_maxwell *maxwell = container_of(eqn, struct wv_maxwell, eqn);

  GKYL_CU_CHECK( maxwell->c == 299792458.0, nfail );
  GKYL_CU_CHECK( maxwell->e_fact == 2.0, nfail );
  GKYL_CU_CHECK( maxwell->b_fact == 2.5, nfail );

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
    eqn->rotate_to_local_func(tau1[d], tau2[d], norm[d], q, q_local);
    gkyl_maxwell_flux(c, e_fact, b_fact, q_local, flux_local);
    eqn->rotate_to_global_func(tau1[d], tau2[d], norm[d], flux_local, flux);

    for (int m=0; m<8; ++m)
      GKYL_CU_CHECK( gkyl_compare(flux[m], fluxes[d][m], 1e-15) , nfail );
  }
}

int cu_maxwell_test(const struct gkyl_wv_eqn *eqn)
{
  int *nfail_dev = (int *) gkyl_cu_malloc(sizeof(int));
  ker_cu_maxwell_test<<<1,1>>>(eqn, nfail_dev);

  int nfail;
  gkyl_cu_memcpy(&nfail, nfail_dev, sizeof(int), GKYL_CU_MEMCPY_D2H);
  gkyl_cu_free(nfail_dev);

  return nfail;
}

