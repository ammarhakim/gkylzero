#include <acutest.h>

#include <gkyl_alloc.h>
#include <gkyl_basis.h>
#include <gkyl_vlasov_mom.h>

void
test_vlasov_mom()
{
  int polyOrder = 2;
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, 1, polyOrder); // 1X
  gkyl_cart_modal_serendip(&pbasis, 4, polyOrder); // 1X3V

  struct gkyl_mom_type *m2ij = gkyl_vlasov_mom_new(&cbasis, &pbasis, "M2ij");

  TEST_CHECK( m2ij->cdim == 1 );
  TEST_CHECK( m2ij->pdim == 4 );
  TEST_CHECK( m2ij->polyOrder == 2 );
  TEST_CHECK( m2ij->num_config == cbasis.numBasis );
  TEST_CHECK( m2ij->num_phase == pbasis.numBasis );
  TEST_CHECK( m2ij->num_mom == 6 );
  
  struct gkyl_mom_type *m3ijk = gkyl_vlasov_mom_new(&cbasis, &pbasis, "M3ijk");
  TEST_CHECK( m3ijk->num_mom == 10 );

  gkyl_mom_type_release(m2ij);
  gkyl_mom_type_release(m3ijk);
}

#ifdef GKYL_HAVE_CUDA

int cu_vlasov_mom_test(const struct gkyl_mom_type *mom);

void
test_cu_vlasov_mom()
{
  int polyOrder = 2;
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, 1, polyOrder); // 1X
  gkyl_cart_modal_serendip(&pbasis, 4, polyOrder); // 1X3V

  struct gkyl_mom_type *m2ij = gkyl_vlasov_mom_cu_dev_new(&cbasis, &pbasis, "M2ij");

  int nfail = cu_vlasov_mom_test(m2ij);
  TEST_CHECK( nfail == 0 );

  gkyl_cu_free(m2ij);
}

#endif

TEST_LIST = {
  { "vlasov_mom", test_vlasov_mom },
#ifdef GKYL_HAVE_CUDA
  { "cu_vlasov_mom", test_cu_vlasov_mom },
#endif  
  { NULL, NULL },
};
