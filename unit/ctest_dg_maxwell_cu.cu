/* -*- c -*- */

#include <gkylzero.h>
#include <stdio.h>

extern "C" {
    int cu_maxwell_test(const struct gkyl_dg_eqn *eqn);
}

__global__
void ker_cu_maxwell_test(const struct gkyl_dg_eqn *eqn, int *nfail)
{
  *nfail = 0;

  GKYL_CU_CHECK( eqn->num_equations == 8, nfail );
}

int cu_maxwell_test(const struct gkyl_dg_eqn *arr)
{
  int *nfail_dev = (int *) gkyl_cu_malloc(sizeof(int));
  ker_cu_maxwell_test<<<1,1>>>(arr, nfail_dev);

  int nfail;
  gkyl_cu_memcpy(&nfail, nfail_dev, sizeof(int), GKYL_CU_MEMCPY_D2H);
  gkyl_cu_free(nfail_dev);

  return nfail;
}

