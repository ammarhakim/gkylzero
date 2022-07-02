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
  TEST_CHECK( basis1.b_type == GKYL_BASIS_MODAL_SERENDIPITY );

  double z[basis1.num_basis], b[basis1.num_basis];

  z[0] = 0.0;
  basis1.eval(z, b);

  TEST_CHECK( gkyl_compare(b[0], 1/sqrt(2.0), 1e-15) );
  TEST_CHECK( b[1] == 0.0 );

  z[0] = 0.5; basis1.eval(z, b);
  
  TEST_CHECK( gkyl_compare(b[0], 1/sqrt(2.0), 1e-15) );
  TEST_CHECK( b[1] == sqrt(3.0/2.0)*0.5 );

  double nodes[basis1.ndim*basis1.num_basis];
  basis1.node_list(nodes);

  TEST_CHECK( nodes[0] == -1 );
  TEST_CHECK( nodes[1] == 1 );

  double f[basis1.num_basis];
  for (int i=0; i<basis1.num_basis; ++i) f[i] = 1.0;

  z[0] = 0.5;
  TEST_CHECK ( gkyl_compare(1.319479216882342, basis1.eval_expand(z, f), 1e-15) );
  TEST_CHECK ( gkyl_compare(1.224744871391589, basis1.eval_grad_expand(0, z, f), 1e-15) );  

  z[0] = 0.25;
  TEST_CHECK ( gkyl_compare(1.013292999034445, basis1.eval_expand(z, f), 1e-15) );
  TEST_CHECK ( gkyl_compare(1.224744871391589, basis1.eval_grad_expand(0, z, f), 1e-15) );
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
  TEST_CHECK( basis.b_type == GKYL_BASIS_MODAL_SERENDIPITY );

  double z[basis.ndim], b[basis.num_basis];

  z[0] = 0.0; z[1] = 0.0;
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

  basis.flip_odd_sign(0, fin, fout);
  TEST_CHECK( fin[0] == fout[0] );
  TEST_CHECK( -fin[1] == fout[1] );
  TEST_CHECK( fin[2] == fout[2] );  
  TEST_CHECK( -fin[3] == fout[3] );
  TEST_CHECK( fin[4] == fout[4] );
  TEST_CHECK( fin[5] == fout[5] );
  TEST_CHECK( fin[6] == fout[6] );
  TEST_CHECK( -fin[7] == fout[7] );

  basis.flip_odd_sign(1, fin, fout);
  TEST_CHECK( fin[0] == fout[0] );
  TEST_CHECK( fin[1] == fout[1] );
  TEST_CHECK( -fin[2] == fout[2] );  
  TEST_CHECK( -fin[3] == fout[3] );
  TEST_CHECK( fin[4] == fout[4] );
  TEST_CHECK( fin[5] == fout[5] );
  TEST_CHECK( -fin[6] == fout[6] );
  TEST_CHECK( fin[7] == fout[7] );

  basis.flip_even_sign(0, fin, fout);
  TEST_CHECK( -fin[0] == fout[0] );
  TEST_CHECK( fin[1] == fout[1] );
  TEST_CHECK( -fin[2] == fout[2] );  
  TEST_CHECK( fin[3] == fout[3] );
  TEST_CHECK( -fin[4] == fout[4] );
  TEST_CHECK( -fin[5] == fout[5] );
  TEST_CHECK( -fin[6] == fout[6] );
  TEST_CHECK( fin[7] == fout[7] );

  basis.flip_even_sign(1, fin, fout);
  TEST_CHECK( -fin[0] == fout[0] );
  TEST_CHECK( -fin[1] == fout[1] );
  TEST_CHECK( fin[2] == fout[2] );  
  TEST_CHECK( fin[3] == fout[3] );
  TEST_CHECK( -fin[4] == fout[4] );
  TEST_CHECK( -fin[5] == fout[5] );
  TEST_CHECK( fin[6] == fout[6] );
  TEST_CHECK( -fin[7] == fout[7] );  

  double nodes[basis.ndim*basis.num_basis];
  basis.node_list(nodes);

  TEST_CHECK( nodes[0] == -1 );
  TEST_CHECK( nodes[1] == -1 );

  TEST_CHECK( nodes[2] == 0 );
  TEST_CHECK( nodes[3] == -1 );

  TEST_CHECK( nodes[4] == 1 );
  TEST_CHECK( nodes[5] == -1 );

  double f[basis.num_basis];
  for (int i=0; i<basis.num_basis; ++i) f[i] = 1.0;

  z[0] = 0.5; z[1] = 0.5;
  TEST_CHECK ( gkyl_compare(1.219455447459001, basis.eval_expand(z, f), 1e-15) );
  TEST_CHECK ( gkyl_compare(4.503383682599098, basis.eval_grad_expand(0, z, f), 1e-15) );
  TEST_CHECK ( gkyl_compare(4.503383682599098, basis.eval_grad_expand(1, z, f), 1e-15) );

  z[0] = 0.25; z[1] = -0.75;
  TEST_CHECK ( gkyl_compare(0.4723022336170485, basis.eval_expand(z, f), 1e-15) );
  TEST_CHECK ( gkyl_compare(0.1559433418554234, basis.eval_grad_expand(0, z, f), 1e-15) );
  TEST_CHECK ( gkyl_compare(-3.150527379222043, basis.eval_grad_expand(1, z, f), 1e-15) );
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
  TEST_CHECK( basis.b_type == GKYL_BASIS_MODAL_TENSOR );

  double z[basis.ndim], b[basis.num_basis];

  z[0] = 0.0; z[1] = 0.0;
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

  double f[basis.num_basis];
  for (int i=0; i<basis.num_basis; ++i) f[i] = 1.0;

  z[0] = 0.5; z[1] = 0.5;
  TEST_CHECK ( gkyl_compare(1.258517947459001, basis.eval_expand(z, f), 1e-15) );
  TEST_CHECK ( gkyl_compare(4.034633682599098, basis.eval_grad_expand(0, z, f), 1e-15) );
  TEST_CHECK ( gkyl_compare(4.034633682599098, basis.eval_grad_expand(1, z, f), 1e-15) );  

  z[0] = 0.25; z[1] = -0.75;
  TEST_CHECK ( gkyl_compare(0.1231811398670484, basis.eval_expand(z, f), 1e-15) );
  TEST_CHECK ( gkyl_compare(0.8004745918554234, basis.eval_grad_expand(0, z, f), 1e-15) );
  TEST_CHECK ( gkyl_compare(-0.8653711292220426, basis.eval_grad_expand(1, z, f), 1e-15) );
}

void
test_gk_hyb()
{
  struct gkyl_basis basis;
  gkyl_cart_modal_gk_hybrid(&basis, 2);

  TEST_CHECK( basis.ndim == 2 );
  TEST_CHECK( basis.poly_order == 1 );
  TEST_CHECK( basis.num_basis == 6 );
  TEST_CHECK( strcmp(basis.id, "gk_hybrid") == 0 );
  TEST_CHECK( basis.b_type == GKYL_BASIS_MODAL_GK_HYBRID );

  double z[basis.ndim], b[basis.num_basis];

  z[0] = 0.0; z[1] = 0.0;
  basis.eval(z, b);

  TEST_CHECK( gkyl_compare(0.5, b[0], 1e-15) );
  TEST_CHECK( gkyl_compare(0.0, b[1], 1e-15) );
  TEST_CHECK( gkyl_compare(0.0, b[2], 1e-15) );
  TEST_CHECK( gkyl_compare(0.0, b[3], 1e-15) );
  TEST_CHECK( gkyl_compare(-sqrt(5.0)/4, b[4], 1e-15) );
  TEST_CHECK( gkyl_compare(0.0, b[6], 1e-15) );

  double fin[basis.num_basis], fout[basis.num_basis];
  for (int i=0; i<basis.num_basis; ++i) {
    fin[i] = 1.0;
    fout[i] = 0.0;
  }

  basis.flip_odd_sign(0, fin, fout);
  TEST_CHECK( fin[0] == fout[0] );
  TEST_CHECK( -fin[1] == fout[1] );
  TEST_CHECK( fin[2] == fout[2] );  
  TEST_CHECK( -fin[3] == fout[3] );
  TEST_CHECK( fin[4] == fout[4] );
  TEST_CHECK( -fin[5] == fout[5] );

}

#ifdef GKYL_HAVE_CUDA

int dev_cu_ser_2d(struct gkyl_basis *basis);

void
test_cu_ser_2d()
{
  struct gkyl_basis *basis = gkyl_cu_malloc(sizeof(struct gkyl_basis));
  gkyl_cart_modal_serendip_cu_dev(basis, 2, 2);

  int nfail = dev_cu_ser_2d(basis);

  TEST_CHECK(nfail == 0);

  gkyl_cu_free(basis);
}

#endif

TEST_LIST = {
  { "ser_1d", test_ser_1d },
  { "ser_2d", test_ser_2d },
  { "ten_2d", test_ten_2d },
  { "gk_hyb", test_gk_hyb },
#ifdef GKYL_HAVE_CUDA
  { "cu_ser_2d", test_cu_ser_2d },
#endif    
  { NULL, NULL },
};
