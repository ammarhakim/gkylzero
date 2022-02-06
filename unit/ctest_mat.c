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
  gkyl_mem_buff ipiv = gkyl_mem_buff_new(sizeof(long[3]));

  // solve linear system: sol : matrix( [-1, 1, 0] )
  bool status = gkyl_mat_linsolve_lu(A, x, gkyl_mem_buff_data(ipiv));

  TEST_CHECK( status );

  //gkyl_mat_show("A", stdout, A);
  //gkyl_mat_show("x", stdout, x);

  TEST_CHECK( gkyl_compare(gkyl_mat_get(x,0,0), -1.0, 1e-15) );
  TEST_CHECK( gkyl_compare(gkyl_mat_get(x,1,0), 1.0, 1e-15) );
  TEST_CHECK( gkyl_compare(gkyl_mat_get(x,2,0), 0.0, 1e-15) );
  
  gkyl_mat_release(A);
  gkyl_mat_release(x);
  gkyl_mem_buff_release(ipiv);
}

void
test_nmat_base()
{
  // 5 matrices with shape 10x20
  struct gkyl_nmat *nmat = gkyl_nmat_new(5, 10, 20);

  TEST_CHECK( false == gkyl_nmat_is_cu_dev(nmat) );

  TEST_CHECK( 5 == nmat->num );
  TEST_CHECK( 10 == nmat->nr );
  TEST_CHECK( 20 == nmat->nc );

  struct gkyl_mat m = gkyl_nmat_get(nmat, 0);
  TEST_CHECK( 10 == m.nr );
  TEST_CHECK( 20 == m.nc );

  for (size_t n=0; n<nmat->num; ++n) {
    struct gkyl_mat m = gkyl_nmat_get(nmat, n);
    
    for (size_t j=0; j<nmat->nc; ++j)
      for (size_t i=0; i<nmat->nr; ++i)
        gkyl_mat_set(&m, i, j, n*0.5);
  }

  for (size_t n=0; n<nmat->num; ++n) {
    struct gkyl_mat m = gkyl_nmat_get(nmat, n);
    
    for (size_t j=0; j<nmat->nc; ++j)
      for (size_t i=0; i<nmat->nr; ++i)
        TEST_CHECK ( n*0.5 == gkyl_mat_get(&m, i, j) );
  }  

  for (size_t n=0; n<nmat->num; ++n)
    for (size_t i=0; i<nmat->nr*nmat->nc; ++i)
      TEST_CHECK( nmat->mptr[n][i] == n*0.5 );

  // copy matrix
  struct gkyl_nmat *ncpy = gkyl_nmat_new(nmat->num, nmat->nr, nmat->nc);
  gkyl_nmat_copy(ncpy, nmat);

  for (size_t n=0; n<ncpy->num; ++n) {
    struct gkyl_mat m = gkyl_nmat_get(ncpy, n);
    
    for (size_t j=0; j<ncpy->nc; ++j)
      for (size_t i=0; i<ncpy->nr; ++i)
        TEST_CHECK ( n*0.5 == gkyl_mat_get(&m, i, j) );
  }  

  for (size_t n=0; n<ncpy->num; ++n)
    for (size_t i=0; i<ncpy->nr*ncpy->nc; ++i)
      TEST_CHECK( ncpy->mptr[n][i] == n*0.5 );  

  gkyl_nmat_release(nmat);
  gkyl_nmat_release(ncpy);
}

void
test_nmat_linsolve_(bool pre_alloc)
{
  struct gkyl_nmat *As = gkyl_nmat_new(5, 3, 3);
  struct gkyl_nmat *xs = gkyl_nmat_new(5, 3, 1);

  for (int n=0; n<As->num; ++n) {
    struct gkyl_mat A = gkyl_nmat_get(As, n);
    
    double val = 1.0;
    // A : matrix( [1,2,3], [4,5,6], [7,8,10] );
    for (int i=0; i<A.nr; ++i)  
      for (int j=0; j<A.nc; ++j) {
        gkyl_mat_set(&A,i,j,val);
        val += 1.0;
      }
    gkyl_mat_set(&A,2,2,10.0); // ensures determinant is not zero
  }

  // x = matrix( [1, 1, 1] );
  for (size_t n=0; n<As->num; ++n) {
    struct gkyl_mat x = gkyl_nmat_get(xs, n);
    gkyl_mat_clear(&x, 1.0);
  }

  // solve all linear systems: sol : matrix( [-1, 1, 0] )
  bool status = false;

  if (pre_alloc) {
    gkyl_nmat_mem *mem = gkyl_nmat_linsolve_lu_alloc(As->num, As->nr);
    status = gkyl_nmat_linsolve_lu_pa(mem, As, xs);
    gkyl_nmat_linsolve_lu_release(mem);
  }
  else {
    status = gkyl_nmat_linsolve_lu(As, xs);
  }

  TEST_CHECK( status );

  //gkyl_mat_show("A", stdout, A);
  //gkyl_mat_show("x", stdout, x);

  for (int n=0; n<xs->num; ++n) {
    struct gkyl_mat x = gkyl_nmat_get(xs, n);

    TEST_CHECK( gkyl_compare(gkyl_mat_get(&x,0,0), -1.0, 1e-15) );
    TEST_CHECK( gkyl_compare(gkyl_mat_get(&x,1,0), 1.0, 1e-15) );
    TEST_CHECK( gkyl_compare(gkyl_mat_get(&x,2,0), 0.0, 1e-15) );
  }
  
  gkyl_nmat_release(As);
  gkyl_nmat_release(xs);
}

void test_nmat_linsolve() { test_nmat_linsolve_(false); }
void test_nmat_linsolve_pa() { test_nmat_linsolve_(true); }

#ifdef GKYL_HAVE_CUDA

void
test_cu_nmat_base()
{
  // 5 matrices with shape 10x20
  struct gkyl_nmat *nmat = gkyl_nmat_cu_dev_new(5, 10, 20);

  TEST_CHECK( 5 == nmat->num );
  TEST_CHECK( 10 == nmat->nr );
  TEST_CHECK( 20 == nmat->nc );

  TEST_CHECK( gkyl_nmat_is_cu_dev(nmat) == true );

  // create host-side matrix
  struct gkyl_nmat *h1 = gkyl_nmat_new(5, 10, 20);
  for (size_t n=0; n<h1->num; ++n) {
    struct gkyl_mat m = gkyl_nmat_get(h1, n);
    
    for (size_t j=0; j<h1->nc; ++j)
      for (size_t i=0; i<h1->nr; ++i)
        gkyl_mat_set(&m, i, j, n*0.5);
  }

  // copy to device
  gkyl_nmat_copy(nmat, h1);

  // copy it back to host
  struct gkyl_nmat *h2 = gkyl_nmat_new(5, 10, 20);
  gkyl_nmat_copy(h2, nmat);

  // check
  for (size_t n=0; n<h1->num; ++n) {
    struct gkyl_mat m1 = gkyl_nmat_get(h1, n);
    struct gkyl_mat m2 = gkyl_nmat_get(h2, n);
    
    for (size_t j=0; j<h1->nc; ++j)
      for (size_t i=0; i<h1->nr; ++i)
        TEST_CHECK( gkyl_mat_get(&m1, i, j) == gkyl_mat_get(&m2, i, j) );
  }

  gkyl_nmat_release(nmat);
  gkyl_nmat_release(h1);
  gkyl_nmat_release(h2);
}

void
test_cu_nmat_linsolve_(bool pre_alloc)
{
  struct gkyl_nmat *As = gkyl_nmat_new(5, 3, 3);
  struct gkyl_nmat *xs = gkyl_nmat_new(5, 3, 1);

  for (int n=0; n<As->num; ++n) {
    struct gkyl_mat A = gkyl_nmat_get(As, n);
    
    double val = 1.0;
    // A : matrix( [1,2,3], [4,5,6], [7,8,10] );
    for (int i=0; i<A.nr; ++i)  
      for (int j=0; j<A.nc; ++j) {
        gkyl_mat_set(&A,i,j,val);
        val += 1.0;
      }
    gkyl_mat_set(&A,2,2,10.0); // ensures determinant is not zero
  }

  // x = matrix( [1, 1, 1] );
  for (size_t n=0; n<As->num; ++n) {
    struct gkyl_mat x = gkyl_nmat_get(xs, n);
    gkyl_mat_clear(&x, 1.0);
  }

  // copy data to GPU matrices
  struct gkyl_nmat *As_d = gkyl_nmat_cu_dev_new(5, 3, 3);
  struct gkyl_nmat *xs_d = gkyl_nmat_cu_dev_new(5, 3, 1);

  gkyl_nmat_copy(As_d, As);
  gkyl_nmat_copy(xs_d, xs);

  // solve all linear systems: sol : matrix( [-1, 1, 0] )
  bool status = false;

  if (pre_alloc) {
    gkyl_nmat_mem *mem = gkyl_nmat_linsolve_lu_alloc_cu_dev(As->num, As->nr);
    status = gkyl_nmat_linsolve_lu_pa(mem, As_d, xs_d);
    gkyl_nmat_linsolve_lu_release(mem);
  }
  else {
    status = gkyl_nmat_linsolve_lu(As_d, xs_d);
  }

  TEST_CHECK( status );

  // copy solution back
  gkyl_nmat_copy(As, As_d);
  gkyl_nmat_copy(xs, xs_d);

  for (int n=0; n<xs->num; ++n) {
    struct gkyl_mat x = gkyl_nmat_get(xs, n);

    TEST_CHECK( gkyl_compare(gkyl_mat_get(&x,0,0), -1.0, 1e-15) );
    TEST_CHECK( gkyl_compare(gkyl_mat_get(&x,1,0), 1.0, 1e-15) );
    TEST_CHECK( gkyl_compare(gkyl_mat_get(&x,2,0), 0.0, 1e-15) );
  }
  
  gkyl_nmat_release(As);
  gkyl_nmat_release(xs);

  gkyl_nmat_release(As_d);
  gkyl_nmat_release(xs_d);
}

void test_cu_nmat_linsolve() { test_cu_nmat_linsolve_(false); }
void test_cu_nmat_linsolve_pa() { test_cu_nmat_linsolve_(true); }

#endif

TEST_LIST = {
  { "mat_base", test_mat_base },
  { "mat_mm_op", test_mat_mm_op },
  { "mat_linsolve", test_mat_linsolve },
  { "nmat_base", test_nmat_base },
  { "nmat_linsolve", test_nmat_linsolve },
  { "nmat_linsolve_pa", test_nmat_linsolve_pa },
#ifdef GKYL_HAVE_CUDA
  { "cu_nmat_base", test_cu_nmat_base },
  { "cu_nmat_linsolve", test_cu_nmat_linsolve },
#endif
  { NULL, NULL },
};
