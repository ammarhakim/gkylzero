#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_mat.h>

void
test_mat_base()
{
  struct gkyl_mat *m = gkyl_mat_new(10, 20, 0.25);

  TEST_CHECK( 10 == m->nr );
  TEST_CHECK( 20 == m->nc );

  for (size_t j=0; j<m->nc; ++j)
    for (size_t i=0; i<m->nr; ++i)
      TEST_CHECK ( 0.25 == gkyl_mat_get(m, i, j) );

  gkyl_mat_clear(m, 0.1);

  for (size_t j=0; j<m->nc; ++j)
    for (size_t i=0; i<m->nr; ++i)
      TEST_CHECK ( 0.1 == gkyl_mat_get(m, i, j) );

  size_t count = 0;
  for (size_t j=0; j<m->nc; ++j)
    for (size_t i=0; i<m->nr; ++i)
      gkyl_mat_set(m, i, j, count++);

  count = 0;
  for (size_t j=0; j<m->nc; ++j) {
    const double* col = gkyl_mat_get_ccol(m, j);
    for (size_t i=0; i<m->nr; ++i)
      TEST_CHECK( col[i] == count++ );
  }

  count = 0;
  for (size_t i=0; i<m->nr*m->nc; ++i)
    TEST_CHECK( m->data[i] == i );

  gkyl_mat_diag(m, 1.0);
  
  for (size_t j=0; j<m->nc; ++j)
    for (size_t i=0; i<m->nr; ++i)
      TEST_CHECK( gkyl_mat_get(m, i, j) == ( i==j ? 1.0 : 0.0 ) );

  struct gkyl_mat *m2 = gkyl_mat_clone(m);

  TEST_CHECK( m2->nr == m->nr );
  TEST_CHECK( m2->nc == m->nc );

  for (size_t j=0; j<m->nc; ++j)
    for (size_t i=0; i<m->nr; ++i)
      TEST_CHECK( gkyl_mat_get(m, i, j) == gkyl_mat_get(m2, i, j) );

  gkyl_mat_release(m);
  gkyl_mat_release(m2);
}

void
test_mat_mm_op()
{
  struct gkyl_mat *A = gkyl_mat_new(2, 3, 0.0);
  struct gkyl_mat *B = gkyl_mat_new(3, 2, 0.0);
  struct gkyl_mat *C = gkyl_mat_new(2, 2, 0.0);
  struct gkyl_mat *D = gkyl_mat_new(3, 3, 0.0);

  double val = 1.0;

  // A : matrix( [1,2,3], [4,5,6] );
  for (int i=0; i<A->nr; ++i)  
    for (int j=0; j<A->nc; ++j) {
      gkyl_mat_set(A,i,j,val);
      val += 1.0;
    }

  // B : matrix( [7,8], [9,10], [11,12] );
  for (int i=0; i<B->nr; ++i)  
    for (int j=0; j<B->nc; ++j) {
      gkyl_mat_set(B,i,j,val);
      val += 1.0;
    }

  // C = 0.5*A*B + 0.0*C
  gkyl_mat_mm(0.5, 0.0, GKYL_NO_TRANS, A, GKYL_NO_TRANS, B, C);

  // C : matrix( [29.0, 32.0], [69.5, 77.0] )
  TEST_CHECK( gkyl_mat_get(C, 0, 0) == 29.0 );
  TEST_CHECK( gkyl_mat_get(C, 0, 1) == 32.0 );
  TEST_CHECK( gkyl_mat_get(C, 1, 0) == 69.5 );
  TEST_CHECK( gkyl_mat_get(C, 1, 1) == 77.0 );

  // D = 0.5*A'*B'
  gkyl_mat_mm(0.5, 0.0, GKYL_TRANS, A, GKYL_TRANS, B, D);

  // D : matrix( [ 19.5  24.5  29.5 ], [ 27.0  34.0  41.0 ], [ 34.5  43.5  52.5 ] )
  TEST_CHECK( gkyl_mat_get(D, 0, 0) == 19.5 );
  TEST_CHECK( gkyl_mat_get(D, 0, 1) == 24.5 );
  TEST_CHECK( gkyl_mat_get(D, 0, 2) == 29.5 );

  TEST_CHECK( gkyl_mat_get(D, 1, 0) == 27.0 );
  TEST_CHECK( gkyl_mat_get(D, 1, 1) == 34.0 );
  TEST_CHECK( gkyl_mat_get(D, 1, 2) == 41.0 );

  TEST_CHECK( gkyl_mat_get(D, 2, 0) == 34.5 );
  TEST_CHECK( gkyl_mat_get(D, 2, 1) == 43.5 );
  TEST_CHECK( gkyl_mat_get(D, 2, 2) == 52.5 );
  
  gkyl_mat_release(A);
  gkyl_mat_release(B);
  gkyl_mat_release(C);
  gkyl_mat_release(D);
}

void
test_mat_linsolve()
{
  struct gkyl_mat *A = gkyl_mat_new(3, 3, 0.0);
  struct gkyl_mat *x = gkyl_mat_new(3, 1, 0.0);

  double val = 1.0;

  // A : matrix( [1,2,3], [4,5,6], [7,8,10] );
  for (int i=0; i<A->nr; ++i)  
    for (int j=0; j<A->nc; ++j) {
      gkyl_mat_set(A,i,j,val);
      val += 1.0;
    }
  gkyl_mat_set(A,2,2,10.0); // ensures determinant is not zero

  //gkyl_mat_show("A", stdout, A);

  // x = matrix( [1, 1, 1] );
  gkyl_mat_clear(x, 1.0);

  // memory for pivot vector
  void *ipiv = gkyl_malloc(sizeof(long[A->nr]));

  // solve linear system: sol : matrix( [-1, 1, 0] )
  bool status = gkyl_mat_linsolve_lu(A, x, ipiv);

  TEST_CHECK( status );

  //gkyl_mat_show("A", stdout, A);
  //gkyl_mat_show("x", stdout, x);

  TEST_CHECK( gkyl_compare(gkyl_mat_get(x,0,0), -1.0, 1e-15) );
  TEST_CHECK( gkyl_compare(gkyl_mat_get(x,1,0), 1.0, 1e-15) );
  TEST_CHECK( gkyl_compare(gkyl_mat_get(x,2,0), 0.0, 1e-15) );
  
  gkyl_mat_release(A);
  gkyl_mat_release(x);
  gkyl_free(ipiv);
}

TEST_LIST = {
  { "mat_base", test_mat_base },
  { "mat_mm_op", test_mat_mm_op },
  { "mat_linsolve", test_mat_linsolve },
  { NULL, NULL },
};
