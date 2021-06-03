/* -*- c -*- */

#include <stdio.h>
#include <gkylzero.h>
#include <gkyl_dg_maxwell_priv.h>

extern "C" {
    int cu_maxwell_test(const struct gkyl_dg_eqn *eqn);
}

__global__
void ker_cu_maxwell_test(const struct gkyl_dg_eqn *eqn, int *nfail)
{
  *nfail = 0;

  GKYL_CU_CHECK( eqn->num_equations == 8, nfail );
  struct dg_maxwell *maxwell = container_of(eqn, struct dg_maxwell, eqn);

  GKYL_CU_CHECK( maxwell->maxwell_data.c == 1.0, nfail );
  GKYL_CU_CHECK( maxwell->maxwell_data.chi == 0.5, nfail );
  GKYL_CU_CHECK( maxwell->maxwell_data.gamma == 0.25, nfail );
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

