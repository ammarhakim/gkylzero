#include <acutest.h>
#include <gkyl_basis.h>
#include <gkyl_dg_maxwell.h>

void
test_dg_max()
{
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 1, 1);

  struct gkyl_dg_eqn* eqn = gkyl_dg_maxwell_new(&basis, 1.0, 0.5, 0.25);

  TEST_CHECK( eqn->num_equations == 8 );

  gkyl_dg_eqn_release(eqn);
}

#ifdef GKYL_HAVE_CUDA

void
test_cu_dg_max()
{
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 1, 1);

  struct gkyl_dg_eqn* eqn = gkyl_dg_maxwell_cu_dev_new(&basis, 1.0, 0.5, 0.25);
  
}

#endif

TEST_LIST = {
  { "dg_max", test_dg_max },
#ifdef GKYL_HAVE_CUDA
  { "cu_dg_max", test_cu_dg_max },
#endif  
  { NULL, NULL },
};
