/* -*- c -*- */

#include <stdio.h>

extern "C" {
#include <gkyl_mom_type.h>
#include <gkyl_alloc.h>
#include <gkyl_util.h>
  int cu_mom_vlasov_test(const struct gkyl_mom_type *mom);
}

__global__
void ker_cu_mom_vlasov_test(const struct gkyl_mom_type *m2ij, int *nfail)
{
  *nfail = 0;

  GKYL_CU_CHECK( m2ij->cdim == 1, nfail );
  GKYL_CU_CHECK( m2ij->pdim == 4, nfail );
  GKYL_CU_CHECK( m2ij->poly_order == 2, nfail );
  GKYL_CU_CHECK( m2ij->num_config == 3, nfail );
  GKYL_CU_CHECK( m2ij->num_phase == 48, nfail );
  GKYL_CU_CHECK( m2ij->num_mom == 6, nfail );
}

int cu_mom_vlasov_test(const struct gkyl_mom_type *mom)
{
  int *nfail_dev = (int *) gkyl_cu_malloc(sizeof(int));
  ker_cu_mom_vlasov_test<<<1,1>>>(mom, nfail_dev);

  int nfail;
  gkyl_cu_memcpy(&nfail, nfail_dev, sizeof(int), GKYL_CU_MEMCPY_D2H);
  gkyl_cu_free(nfail_dev);

  return nfail;
}

