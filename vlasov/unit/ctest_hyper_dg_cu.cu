#include <gkylzero.h>

#define TEST_NO_MAIN
#include <acutest.h>

extern "C" {
#include <gkyl_hyper_dg_priv.h>
#include <gkyl_dg_vlasov_priv.h>
    
int hyper_dg_kernel_test(const gkyl_hyper_dg *slvr);
}

__global__
void ker_cu_hyper_dg_kernel_test(const gkyl_hyper_dg *slvr, int *nfail)
{
  *nfail = 0;
  
  GKYL_CU_CHECK( slvr->ndim == 5, nfail );
  GKYL_CU_CHECK( slvr->num_basis == 32, nfail );
  GKYL_CU_CHECK( slvr->num_up_dirs == 5, nfail );

  GKYL_CU_CHECK( slvr->grid.ndim == 5, nfail );

  GKYL_CU_CHECK( slvr->equation->num_equations == 1, nfail );

  // DO NOT DO THIS IN PRODUCTION! ONLY FOR TESTING
  struct dg_vlasov *vlasov = container_of(slvr->equation, struct dg_vlasov, eqn);

  GKYL_CU_CHECK( vlasov->cdim == 2, nfail );
  GKYL_CU_CHECK( vlasov->pdim == 5, nfail );
  GKYL_CU_CHECK( vlasov->conf_range.volume == 8*8, nfail );
}

int
hyper_dg_kernel_test(const gkyl_hyper_dg *slvr)
{
  int *nfail_dev = (int *) gkyl_cu_malloc(sizeof(int));
  ker_cu_hyper_dg_kernel_test<<<1,1>>>(slvr, nfail_dev);

  int nfail;
  gkyl_cu_memcpy(&nfail, nfail_dev, sizeof(int), GKYL_CU_MEMCPY_D2H);
  gkyl_cu_free(nfail_dev);

  return nfail;  
}

