#include <acutest.h>

#include <gkyl_alloc.h>
#include <gkyl_basis.h>
#include <gkyl_dg_gyrokinetic.h>
#include <gkyl_dg_gyrokinetic_priv.h>

void
test_dg_gyrokinetic()
{
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, 1, 1);
  gkyl_cart_modal_serendip(&pbasis, 2, 1);

  struct gkyl_range crange;
  gkyl_range_init_from_shape(&crange, 1, (int[]) { 100 } );

  double charge = 1.;
  double mass = 1.;

  struct gkyl_dg_eqn* eqn = gkyl_dg_gyrokinetic_new(&cbasis, &pbasis, &crange, charge, mass, false);

  TEST_CHECK( eqn->num_equations == 1 );

  // this is not possible from user code and should NOT be done. This
  // is for testing only
  struct dg_gyrokinetic *gyrokinetic = container_of(eqn, struct dg_gyrokinetic, eqn);

  TEST_CHECK( gyrokinetic->cdim == 1 );
  TEST_CHECK( gyrokinetic->pdim == 2 );
  TEST_CHECK( gyrokinetic->conf_range.volume == 100 );

  gkyl_dg_eqn_release(eqn);
}

#ifdef GKYL_HAVE_CUDA

/* int cu_gyrokinetic_test(const struct gkyl_dg_eqn *eqn); */

/* void */
/* test_cu_dg_gyrokinetic() */
/* { */
/*   struct gkyl_basis cbasis, pbasis; */
/*   gkyl_cart_modal_serendip(&cbasis, 1, 1); */
/*   gkyl_cart_modal_serendip(&pbasis, 2, 1); */

/*   struct gkyl_range crange; */
/*   gkyl_range_init_from_shape(&crange, 1, (int[]) { 100 } ); */

/*   struct gkyl_dg_eqn* eqn = gkyl_dg_gyrokinetic_cu_dev_new(&cbasis, &pbasis, &crange); */

/*   // this is not possible from user code and should NOT be done. This */
/*   // is for testing only */
/*   struct dg_gyrokinetic *gyrokinetic = container_of(eqn, struct dg_gyrokinetic, eqn); */

/*   TEST_CHECK( gyrokinetic->cdim == 1 ); */
/*   TEST_CHECK( gyrokinetic->pdim == 2 ); */
/*   TEST_CHECK( gyrokinetic->conf_range.volume == 100 ); */

/*   int nfail = cu_gyrokinetic_test(eqn->on_dev); */

/*   TEST_CHECK( nfail == 0 ); */

/*   gkyl_dg_eqn_release(eqn); */
/* } */

#endif

TEST_LIST = {
  { "dg_gyrokinetic", test_dg_gyrokinetic },
#ifdef GKYL_HAVE_CUDA
/*  { "cu_dg_gyrokinetic", test_cu_dg_gyrokinetic }, */
#endif  
  { NULL, NULL },
};
