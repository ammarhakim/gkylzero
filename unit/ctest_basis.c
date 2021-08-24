#include <math.h>

#include <acutest.h>
#include <gkyl_basis.h>
#include <gkyl_util.h>

void
test_ser_1d()
{
  struct gkyl_basis basis1;
  gkyl_cart_modal_serendip(&basis1, 1, 1);

  TEST_CHECK( basis1.ndim == 1 );
  TEST_CHECK( basis1.poly_order == 1 );
  TEST_CHECK( basis1.num_basis == 2 );
  TEST_CHECK( strcmp(basis1.id, "serendipity") == 0 );

  double z[1], b[basis1.num_basis];

  z[0] = 0.0; basis1.eval(z, b);

  TEST_CHECK( gkyl_compare(b[0], 1/sqrt(2.0), 1e-15) );
  TEST_CHECK( b[1] == 0.0 );

  z[0] = 0.5; basis1.eval(z, b);
  
  TEST_CHECK( gkyl_compare(b[0], 1/sqrt(2.0), 1e-15) );
  TEST_CHECK( b[1] == sqrt(3.0/2.0)*0.5 );
}

void
test_ser_2d()
{
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 2, 2);

  TEST_CHECK( basis.ndim == 2 );
  TEST_CHECK( basis.poly_order == 2 );
  TEST_CHECK( basis.num_basis == 8 );
  TEST_CHECK( strcmp(basis.id, "serendipity") == 0 );

  double z[basis.ndim], b[basis.num_basis];

  z[1] = 0.0; z[2] = 0.0;
  basis.eval(z, b);

  TEST_CHECK( gkyl_compare(0.5, b[0], 1e-15) );
  TEST_CHECK( gkyl_compare(0.0, b[1], 1e-15) );
  TEST_CHECK( gkyl_compare(0.0, b[2], 1e-15) );
  TEST_CHECK( gkyl_compare(0.0, b[3], 1e-15) );
  TEST_CHECK( gkyl_compare(-sqrt(5.0)/4, b[4], 1e-15) );
  TEST_CHECK( gkyl_compare(-sqrt(5.0)/4, b[5], 1e-15) );
  TEST_CHECK( gkyl_compare(0.0, b[6], 1e-15) );
  TEST_CHECK( gkyl_compare(0.0, b[7], 1e-15) );

  double fin[basis.num_basis], fout[basis.num_basis];
  for (int i=0; i<basis.num_basis; ++i) {
    fin[i] = 1.0;
    fout[i] = 0.0;
  }

  basis.flip_sign(0, fin, fout);
  TEST_CHECK( fin[0] == fout[0] );
  TEST_CHECK( -fin[1] == fout[1] );
  TEST_CHECK( fin[2] == fout[2] );  
  TEST_CHECK( -fin[3] == fout[3] );
  TEST_CHECK( fin[4] == fout[4] );
  TEST_CHECK( fin[5] == fout[5] );
  TEST_CHECK( fin[6] == fout[6] );
  TEST_CHECK( -fin[7] == fout[7] );

  basis.flip_sign(1, fin, fout);
  TEST_CHECK( fin[0] == fout[0] );
  TEST_CHECK( fin[1] == fout[1] );
  TEST_CHECK( -fin[2] == fout[2] );  
  TEST_CHECK( -fin[3] == fout[3] );
  TEST_CHECK( fin[4] == fout[4] );
  TEST_CHECK( fin[5] == fout[5] );
  TEST_CHECK( -fin[6] == fout[6] );
  TEST_CHECK( fin[7] == fout[7] );
  
}

void
test_ten_2d()
{
  struct gkyl_basis basis;
  gkyl_cart_modal_tensor(&basis, 2, 2);

  TEST_CHECK( basis.ndim == 2 );
  TEST_CHECK( basis.poly_order == 2 );
  TEST_CHECK( basis.num_basis == 9 );
  TEST_CHECK( strcmp(basis.id, "tensor") == 0 );

  double z[basis.ndim], b[basis.num_basis];

  z[1] = 0.0; z[2] = 0.0;
  basis.eval(z, b);

  TEST_CHECK( gkyl_compare(0.5, b[0], 1e-15) );
  TEST_CHECK( gkyl_compare(0.0, b[1], 1e-15) );
  TEST_CHECK( gkyl_compare(0.0, b[2], 1e-15) );
  TEST_CHECK( gkyl_compare(0.0, b[3], 1e-15) );
  TEST_CHECK( gkyl_compare(-sqrt(5.0)/4, b[4], 1e-15) );
  TEST_CHECK( gkyl_compare(-sqrt(5.0)/4, b[5], 1e-15) );
  TEST_CHECK( gkyl_compare(0.0, b[6], 1e-15) );
  TEST_CHECK( gkyl_compare(0.0, b[7], 1e-15) );
  TEST_CHECK( gkyl_compare(5.0/8.0, b[8], 1e-15) );
}

TEST_LIST = {
  { "ser_1d", test_ser_1d },
  { "ser_2d", test_ser_2d },
  { "ten_2d", test_ten_2d },
  { NULL, NULL },
};
