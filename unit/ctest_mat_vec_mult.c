#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_mat.h>
//#include <gkyl_rect_grid.h>
//#include <gkyl_rect_decomp.h>
//#include <gkyl_basis.h>
//#include <gkyl_array.h>
//#include <gkyl_array_rio.h>
//#include <gkyl_range.h>


void test_mat_mm()
{
  struct gkyl_mat *A = gkyl_mat_new(4, 4, 1);
  struct gkyl_mat *B = gkyl_mat_new(4, 1, 2);
  struct gkyl_mat *C = gkyl_mat_new(4, 1, 0);
  // C = 1.0*A*B + 0.0*C
  gkyl_mat_mm(1.0, 0.0, GKYL_NO_TRANS, A, GKYL_NO_TRANS, B, C);

  for (size_t j=0; j<C->nc; ++j)
    for (size_t i=0; i<C->nr; ++i)
      TEST_CHECK ( 8 == gkyl_mat_get(C, i, j) );
}

void test_mat_mv()
{

  struct gkyl_mat *A = gkyl_mat_new(4, 4, 1);
  struct gkyl_mat *x = gkyl_mat_new(4, 1, 2);
  struct gkyl_mat *y = gkyl_mat_new(4, 1, 1);
  // y = 1.0*A*x + 0.0*C
  gkyl_mat_mv(1.0, 0.0, GKYL_NO_TRANS, A, x, y);
  for (size_t j=0; j<y->nc; ++j)
    for (size_t i=0; i<y->nr; ++i)
      TEST_CHECK ( 8 == gkyl_mat_get(y, i, j) );
}

void test_nmat_mv()
{

  // n_do matrices with shape 4x4
  int n_do = 3;
  struct gkyl_nmat *nmat_A = gkyl_nmat_new(n_do, 4, 4);
  struct gkyl_nmat *nmat_x = gkyl_nmat_new(n_do, 4, 1);
  struct gkyl_nmat *nmat_y = gkyl_nmat_new(n_do, 4, 1);

  //fill each A matrix with its number
  for (size_t n=0; n<nmat_A->num; ++n) {
    struct gkyl_mat m = gkyl_nmat_get(nmat_A, n);
    for (size_t j=0; j<nmat_A->nc; ++j)
      for (size_t i=0; i<nmat_A->nr; ++i)
        gkyl_mat_set(&m, i, j, n);
  }

  //fill each x vecotr with ones
  for (size_t n=0; n<nmat_x->num; ++n) {
    struct gkyl_mat m = gkyl_nmat_get(nmat_x, n);
    for (size_t j=0; j<nmat_x->nc; ++j)
      for (size_t i=0; i<nmat_x->nr; ++i)
        gkyl_mat_set(&m, i, j, 1);
  }

  //fill each y matrix with 0
  for (size_t n=0; n<nmat_y->num; ++n) {
    struct gkyl_mat m = gkyl_nmat_get(nmat_y, n);
    for (size_t j=0; j<nmat_y->nc; ++j)
      for (size_t i=0; i<nmat_y->nr; ++i)
        gkyl_mat_set(&m, i, j, 0);
  }

  enum gkyl_mat_trans transa = GKYL_NO_TRANS;
  double alpha = 1.0;
  double beta = 0.0;

  gkyl_nmat_mv(alpha, beta, transa, nmat_A, nmat_x, nmat_y);


  //check the expected result
  for (size_t n=0; n<nmat_y->num; ++n) {
    //printf("n = %ld\n", n);
    struct gkyl_mat y = gkyl_nmat_get(nmat_y,n);
    double expected = 4*n;
    //printf("expected = %g\n", expected);
    for (size_t j=0; j<y.nc; ++j)
      for (size_t i=0; i<y.nr; ++i)
        //printf("actual =%g\n",gkyl_mat_get(&y, i, j));
        TEST_CHECK ( expected == gkyl_mat_get(&y, i, j) );
  }


}

void test_nmat_mm()
{

  int n_do = 3;
  struct gkyl_nmat *nmat_A = gkyl_nmat_new(n_do, 4, 3);
  struct gkyl_nmat *nmat_x = gkyl_nmat_new(n_do, 3, 2);
  struct gkyl_nmat *nmat_y = gkyl_nmat_new(n_do, 4, 2);

  //fill each A matrix with its number
  for (size_t n=0; n<nmat_A->num; ++n) {
    struct gkyl_mat m = gkyl_nmat_get(nmat_A, n);
    for (size_t j=0; j<nmat_A->nc; ++j)
      for (size_t i=0; i<nmat_A->nr; ++i)
        gkyl_mat_set(&m, i, j, n);
  }

  //fill each x vector with ones
  for (size_t n=0; n<nmat_x->num; ++n) {
    struct gkyl_mat m = gkyl_nmat_get(nmat_x, n);
    for (size_t j=0; j<nmat_x->nc; ++j)
      for (size_t i=0; i<nmat_x->nr; ++i)
        gkyl_mat_set(&m, i, j, 1);
  }

  //fill each y matrix with 0
  for (size_t n=0; n<nmat_y->num; ++n) {
    struct gkyl_mat m = gkyl_nmat_get(nmat_y, n);
    for (size_t j=0; j<nmat_y->nc; ++j)
      for (size_t i=0; i<nmat_y->nr; ++i)
        gkyl_mat_set(&m, i, j, 0);
  }

  enum gkyl_mat_trans transa = GKYL_NO_TRANS;
  enum gkyl_mat_trans transb = GKYL_NO_TRANS;
  double alpha = 1.0;
  double beta = 0.0;

  gkyl_nmat_mm(alpha, beta, transa, nmat_A, transb, nmat_x, nmat_y);


  //check the expected result
  for (size_t n=0; n<nmat_y->num; ++n) {
    //printf("n = %ld\n", n);
    struct gkyl_mat y = gkyl_nmat_get(nmat_y,n);
    double expected = 3*n;
    //printf("expected = %g\n", expected);
    for (size_t j=0; j<y.nc; ++j)
      for (size_t i=0; i<y.nr; ++i)
        //printf("actual =%g\n",gkyl_mat_get(&y, i, j));
        TEST_CHECK ( expected == gkyl_mat_get(&y, i, j) );
  }


}




#ifdef GKYL_HAVE_CUDA
void test_cu_nmat_mv()
{

  // n_do matrices with shape 4x4
  int n_do = 3;
  struct gkyl_nmat *nmat_A = gkyl_nmat_new(n_do, 4, 4);
  struct gkyl_nmat *nmat_x = gkyl_nmat_new(n_do, 4, 1);
  struct gkyl_nmat *nmat_y = gkyl_nmat_new(n_do, 4, 1);

  struct gkyl_nmat *nmat_Acu = gkyl_nmat_cu_dev_new(n_do, 4, 4);
  struct gkyl_nmat *nmat_xcu = gkyl_nmat_cu_dev_new(n_do, 4, 1);
  struct gkyl_nmat *nmat_ycu = gkyl_nmat_cu_dev_new(n_do, 4, 1);

  //fill each A matrix with its number
  for (size_t n=0; n<nmat_A->num; ++n) {
    struct gkyl_mat m = gkyl_nmat_get(nmat_A, n);
    for (size_t j=0; j<nmat_A->nc; ++j)
      for (size_t i=0; i<nmat_A->nr; ++i)
        gkyl_mat_set(&m, i, j, n);
  }

  //fill each x vecotr with ones
  for (size_t n=0; n<nmat_x->num; ++n) {
    struct gkyl_mat m = gkyl_nmat_get(nmat_x, n);
    for (size_t j=0; j<nmat_x->nc; ++j)
      for (size_t i=0; i<nmat_x->nr; ++i)
        gkyl_mat_set(&m, i, j, 1);
  }

  //fill each y matrix with 0
  for (size_t n=0; n<nmat_y->num; ++n) {
    struct gkyl_mat m = gkyl_nmat_get(nmat_y, n);
    for (size_t j=0; j<nmat_y->nc; ++j)
      for (size_t i=0; i<nmat_y->nr; ++i)
        gkyl_mat_set(&m, i, j, 0);
  }

  enum gkyl_mat_trans transa = GKYL_NO_TRANS;
  enum gkyl_mat_trans transb = GKYL_NO_TRANS;
  double alpha = 1.0;
  double beta = 0.0;


  //copy to device
  gkyl_nmat_copy(nmat_Acu, nmat_A);
  gkyl_nmat_copy(nmat_xcu, nmat_x);
  gkyl_nmat_copy(nmat_ycu, nmat_y);

  gkyl_nmat_mm(alpha, beta, transa, nmat_Acu, transb, nmat_xcu, nmat_ycu);

  //copy to host
  gkyl_nmat_copy(nmat_A, nmat_Acu);
  gkyl_nmat_copy(nmat_x, nmat_xcu);
  gkyl_nmat_copy(nmat_y, nmat_ycu);


  //check the expected result
  for (size_t n=0; n<nmat_y->num; ++n) {
    //printf("n = %ld\n", n);
    struct gkyl_mat y = gkyl_nmat_get(nmat_y,n);
    double expected = 4*n;
    //printf("expected = %g\n", expected);
    for (size_t j=0; j<y.nc; ++j){
      for (size_t i=0; i<y.nr; ++i){
        //printf("actual =%g\n",gkyl_mat_get(&y, i, j));
        double actual = gkyl_mat_get(&y, i, j);
        TEST_CHECK ( expected == actual );
      }
    }
  }
  gkyl_nmat_release(nmat_A);
  gkyl_nmat_release(nmat_x);
  gkyl_nmat_release(nmat_y);
  gkyl_nmat_release(nmat_Acu);
  gkyl_nmat_release(nmat_xcu);
  gkyl_nmat_release(nmat_ycu);


}

void test_cu_nmat_mm()
{

  int n_do = 3;
  struct gkyl_nmat *nmat_A = gkyl_nmat_new(n_do, 4, 3);
  struct gkyl_nmat *nmat_x = gkyl_nmat_new(n_do, 3, 2);
  struct gkyl_nmat *nmat_y = gkyl_nmat_new(n_do, 4, 2);

  struct gkyl_nmat *nmat_Acu = gkyl_nmat_cu_dev_new(n_do, 4, 3);
  struct gkyl_nmat *nmat_xcu = gkyl_nmat_cu_dev_new(n_do, 3, 2);
  struct gkyl_nmat *nmat_ycu = gkyl_nmat_cu_dev_new(n_do, 4, 2);

  //fill each A matrix with its number
  for (size_t n=0; n<nmat_A->num; ++n) {
    struct gkyl_mat m = gkyl_nmat_get(nmat_A, n);
    for (size_t j=0; j<nmat_A->nc; ++j)
      for (size_t i=0; i<nmat_A->nr; ++i)
        gkyl_mat_set(&m, i, j, n);
  }

  //fill each x vecotr with ones
  for (size_t n=0; n<nmat_x->num; ++n) {
    struct gkyl_mat m = gkyl_nmat_get(nmat_x, n);
    for (size_t j=0; j<nmat_x->nc; ++j)
      for (size_t i=0; i<nmat_x->nr; ++i)
        gkyl_mat_set(&m, i, j, 1);
  }

  //fill each y matrix with 0
  for (size_t n=0; n<nmat_y->num; ++n) {
    struct gkyl_mat m = gkyl_nmat_get(nmat_y, n);
    for (size_t j=0; j<nmat_y->nc; ++j)
      for (size_t i=0; i<nmat_y->nr; ++i)
        gkyl_mat_set(&m, i, j, 0);
  }

  enum gkyl_mat_trans transa = GKYL_NO_TRANS;
  double alpha = 1.0;
  double beta = 0.0;


  //copy to device
  gkyl_nmat_copy(nmat_Acu, nmat_A);
  gkyl_nmat_copy(nmat_xcu, nmat_x);
  gkyl_nmat_copy(nmat_ycu, nmat_y);

  gkyl_nmat_mv(alpha, beta, transa, nmat_Acu, nmat_xcu, nmat_ycu);

  //copy to host
  gkyl_nmat_copy(nmat_A, nmat_Acu);
  gkyl_nmat_copy(nmat_x, nmat_xcu);
  gkyl_nmat_copy(nmat_y, nmat_ycu);


  //check the expected result
  for (size_t n=0; n<nmat_y->num; ++n) {
    //printf("n = %ld\n", n);
    struct gkyl_mat y = gkyl_nmat_get(nmat_y,n);
    double expected = 3*n;
    //printf("expected = %g\n", expected);
    for (size_t j=0; j<y.nc; ++j){
      for (size_t i=0; i<y.nr; ++i){
        //printf("actual =%g\n",gkyl_mat_get(&y, i, j));
        double actual = gkyl_mat_get(&y, i, j);
        TEST_CHECK ( expected == actual );
      }
    }
  }
  gkyl_nmat_release(nmat_A);
  gkyl_nmat_release(nmat_x);
  gkyl_nmat_release(nmat_y);
  gkyl_nmat_release(nmat_Acu);
  gkyl_nmat_release(nmat_xcu);
  gkyl_nmat_release(nmat_ycu);


}

#endif




TEST_LIST = {
  { "mv", test_mat_mv},
  { "mm", test_mat_mm},
  { "nmat_mv", test_nmat_mv},
  { "nmat_mm", test_nmat_mm},
  #ifdef GKYL_HAVE_CUDA
  { "cu_nmat_mv", test_cu_nmat_mv},
  { "cu_nmat_mm", test_cu_nmat_mm},
  #endif
  { NULL, NULL },
};
