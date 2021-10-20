#include <acutest.h>

#include <gkyl_alloc.h>
#include <gkyl_basis.h>
#include <gkyl_vlasov_mom.h>
#include <gkyl_vlasov_mom_priv.h>

void
test_vlasov_mom()
{
  int poly_order = 2;
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, 1, poly_order); // 1X
  gkyl_cart_modal_serendip(&pbasis, 4, poly_order); // 1X3V

  struct gkyl_mom_type *m2ij = gkyl_vlasov_mom_new(&cbasis, &pbasis, "M2ij");
  struct gkyl_mom_type *m3ijk = gkyl_vlasov_mom_new(&cbasis, &pbasis, "M3ijk");
  // this is not possible from user code and should NOT be done. This
  // is for testing only
  struct vlasov_mom_type *vlasov_mom_m2ij = container_of(m2ij, struct vlasov_mom_type, momt);
  struct vlasov_mom_type *vlasov_mom_m3ijk = container_of(m3ijk, struct vlasov_mom_type, momt);
  
  TEST_CHECK( vlasov_mom_m2ij->cdim == 1 );
  TEST_CHECK( vlasov_mom_m2ij->pdim == 4 );
  TEST_CHECK( vlasov_mom_m2ij->poly_order == 2 );
  TEST_CHECK( vlasov_mom_m2ij->num_config == cbasis.num_basis );
  TEST_CHECK( vlasov_mom_m2ij->num_phase == pbasis.num_basis );
  TEST_CHECK( vlasov_mom_m2ij->momt.num_mom == 6 );
  
  TEST_CHECK( vlasov_mom_m3ijk->momt.num_mom == 10 );

  gkyl_mom_type_release(m2ij);
  gkyl_mom_type_release(m3ijk);
}

#ifdef GKYL_HAVE_CUDA

int cu_vlasov_mom_test(const struct gkyl_mom_type *mom);

void
test_cu_vlasov_mom()
{
  int poly_order = 2;
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, 1, poly_order); // 1X
  gkyl_cart_modal_serendip(&pbasis, 4, poly_order); // 1X3V

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
