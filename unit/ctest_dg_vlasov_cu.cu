/* -*- c -*- */

#include <stdio.h>
#include <gkylzero.h>

extern "C" {
#include <gkyl_dg_vlasov_priv.h>    
int cu_vlasov_test(const struct gkyl_dg_eqn *eqn);
}

__global__
void ker_cu_vlasov_test(const struct gkyl_dg_eqn *eqn, int *nfail)
{
  *nfail = 0;

  GKYL_CU_CHECK( eqn->num_equations == 1, nfail );

  // DO NOT DO THIS IN PRODUCTION! ONLY FOR TESTING
  struct dg_vlasov *vlasov = container_of(eqn, struct dg_vlasov, eqn);

  GKYL_CU_CHECK( vlasov->cdim == 1, nfail );
  GKYL_CU_CHECK( vlasov->pdim == 2, nfail );
  GKYL_CU_CHECK( vlasov->conf_range.volume == 100, nfail );
  
}

int cu_vlasov_test(const struct gkyl_dg_eqn *eqn)
{
  int *nfail_dev = (int *) gkyl_cu_malloc(sizeof(int));
  ker_cu_vlasov_test<<<1,1>>>(eqn, nfail_dev);

  int nfail;
  gkyl_cu_memcpy(&nfail, nfail_dev, sizeof(int), GKYL_CU_MEMCPY_D2H);
  gkyl_cu_free(nfail_dev);

  return nfail;
}

