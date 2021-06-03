#include <acutest.h>

#include <gkyl_alloc.h>
#include <gkyl_basis.h>
#include <gkyl_dg_vlasov.h>
#include <gkyl_dg_vlasov_priv.h>

void
test_dg_vlasov()
{
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, 1, 1);
  gkyl_cart_modal_serendip(&pbasis, 2, 1);

  struct gkyl_range crange;
  gkyl_range_init_from_shape(&crange, 1, (int[]) { 100 } );

  struct gkyl_dg_eqn* eqn = gkyl_dg_vlasov_new(&cbasis, &pbasis, &crange);

  TEST_CHECK( eqn->num_equations == 1 );

  // this is not possible from user code and should NOT be done. This
  // is for testing only
  struct dg_vlasov *vlasov = container_of(eqn, struct dg_vlasov, eqn);

  TEST_CHECK( vlasov->cdim == 1 );
  TEST_CHECK( vlasov->pdim == 2 );
  TEST_CHECK( vlasov->conf_range.volume == 100 );

  gkyl_dg_eqn_release(eqn);
}

#ifdef GKYL_HAVE_CUDA

int cu_vlasov_test(const struct gkyl_dg_eqn *eqn);

void
test_cu_dg_vlasov()
{
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, 1, 1);
  gkyl_cart_modal_serendip(&pbasis, 2, 1);

  struct gkyl_range crange;
  gkyl_range_init_from_shape(&crange, 1, (int[]) { 100 } );

  struct gkyl_range *crange_cu = gkyl_range_clone_on_cu_dev(&crange);
  struct gkyl_dg_eqn* eqn = gkyl_dg_vlasov_cu_dev_new(&cbasis, &pbasis, crange_cu);

  int nfail = cu_vlasov_test(eqn);

  TEST_CHECK( nfail == 0 );

  gkyl_cu_free(crange_cu);
  gkyl_cu_free(eqn);
}

#endif

TEST_LIST = {
  { "dg_vlasov", test_dg_vlasov },
#ifdef GKYL_HAVE_CUDA
  { "cu_dg_vlasov", test_cu_dg_vlasov },
#endif  
  { NULL, NULL },
};
