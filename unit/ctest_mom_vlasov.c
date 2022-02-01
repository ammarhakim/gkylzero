#include <acutest.h>

#include <gkyl_alloc.h>
#include <gkyl_basis.h>
#include <gkyl_mom_vlasov.h>
#include <gkyl_mom_vlasov_priv.h>

void
test_mom_vlasov()
{
  int poly_order = 2;
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, 1, poly_order); // 1X
  gkyl_cart_modal_serendip(&pbasis, 4, poly_order); // 1X3V

  struct gkyl_mom_type *m2ij = gkyl_mom_vlasov_new(&cbasis, &pbasis, "M2ij");

  TEST_CHECK( m2ij->cdim == 1 );
  TEST_CHECK( m2ij->pdim == 4 );
  TEST_CHECK( m2ij->poly_order == 2 );
  TEST_CHECK( m2ij->num_config == cbasis.num_basis );
  TEST_CHECK( m2ij->num_phase == pbasis.num_basis );
  TEST_CHECK( m2ij->num_mom == 6 );
  
  struct gkyl_mom_type *m3ijk = gkyl_mom_vlasov_new(&cbasis, &pbasis, "M3ijk");
  TEST_CHECK( m3ijk->num_mom == 10 );

  gkyl_mom_type_release(m2ij);
  gkyl_mom_type_release(m3ijk);
}

#ifdef GKYL_HAVE_CUDA

int cu_mom_vlasov_test(const struct gkyl_mom_type *mom);

void
test_cu_mom_vlasov()
{
  int poly_order = 2;
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, 1, poly_order); // 1X
  gkyl_cart_modal_serendip(&pbasis, 4, poly_order); // 1X3V

  struct gkyl_mom_type *m2ij = gkyl_mom_vlasov_cu_dev_new(&cbasis, &pbasis, "M2ij");

  int nfail = cu_mom_vlasov_test(m2ij->on_dev);
  TEST_CHECK( nfail == 0 );

  gkyl_cu_free(m2ij);
}

#endif

TEST_LIST = {
  { "mom_vlasov", test_mom_vlasov },
#ifdef GKYL_HAVE_CUDA
  { "cu_mom_vlasov", test_cu_mom_vlasov },
#endif  
  { NULL, NULL },
};
