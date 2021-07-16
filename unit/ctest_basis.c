#include <math.h>

#include <acutest.h>
#include <gkyl_basis.h>
#include <gkyl_util.h>

void test_ser_1d()
{
  struct gkyl_basis basis1;
  gkyl_cart_modal_serendip(&basis1, 1, 1);

  TEST_CHECK( basis1.ndim == 1 );
  TEST_CHECK( basis1.polyOrder == 1 );
  TEST_CHECK( basis1.numBasis == 2 );
  TEST_CHECK( strcmp(basis1.id, "serendipity") == 0 );

  double z[1], b[basis1.numBasis];

  z[0] = 0.0; basis1.eval(z, b);

  TEST_CHECK( gkyl_compare(b[0], 1/sqrt(2.0), 1e-15) );
  TEST_CHECK( b[1] == 0.0 );

  z[0] = 0.5; basis1.eval(z, b);
  
  TEST_CHECK( gkyl_compare(b[0], 1/sqrt(2.0), 1e-15) );
  TEST_CHECK( b[1] == sqrt(3.0/2.0)*0.5 );
}

TEST_LIST = {
  { "ser_1d", test_ser_1d },
  { NULL, NULL },
};
