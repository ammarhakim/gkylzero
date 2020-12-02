#include <acutest.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_util.h>

void test_array_base()
{
  size_t shape[2] = {10, 20};
  struct gkyl_array *arr = gkyl_array_new(2, sizeof(double), shape);

  TEST_CHECK( arr->rank == 2 );
  TEST_CHECK( arr->shape[0] == 10 );
  TEST_CHECK( arr->shape[1] == 20 );
  TEST_CHECK( arr->elemSz == sizeof(double) );
  TEST_CHECK( arr->size == 20*10 );
  TEST_CHECK( arr->ref_count.count == 1 );

  double *arrData  = arr->data;
  for (unsigned i=0; i<arr->size; ++i)
    arrData[i] = (i+0.5)*0.1;

  // clone array
  struct gkyl_array *brr = gkyl_array_clone(arr);

  TEST_CHECK( brr->rank == 2 );
  TEST_CHECK( brr->shape[0] == 10 );
  TEST_CHECK( brr->shape[1] == 20 );
  TEST_CHECK( brr->elemSz == sizeof(double) );
  TEST_CHECK( brr->size == 20*10 );
  TEST_CHECK( brr->ref_count.count == 1 );

  double *brrData  = brr->data;
  for (unsigned i=0; i<brr->size; ++i)
    TEST_CHECK( brrData[i] == arrData[i] );

  // reset values in brr
  for (unsigned i=0; i<brr->size; ++i)
    brrData[i] = (i-0.5)*0.5;

  gkyl_array_copy(arr, brr);

  for (unsigned i=0; i<arr->size; ++i)
    TEST_CHECK( arrData[i] == brrData[i] );

  // aquire pointer
  struct gkyl_array *crr = gkyl_array_aquire(arr);

  TEST_CHECK( crr->ref_count.count == 2 );
  TEST_CHECK( arr->ref_count.count == 2 );

  struct gkyl_array *drr = gkyl_array_aquire(crr);

  TEST_CHECK( drr->ref_count.count == 3 );
  TEST_CHECK( crr->ref_count.count == 3 );  
  TEST_CHECK( arr->ref_count.count == 3 );  
  
  gkyl_array_release(crr);
  TEST_CHECK( arr->ref_count.count == 2 );
  gkyl_array_release(drr);
  TEST_CHECK( arr->ref_count.count == 1 );
  
  gkyl_array_release(arr);
  gkyl_array_release(brr);
}

void test_array_reshape()
{
  size_t shape[1] = {20};
  struct gkyl_array *arr = gkyl_array_new(1, sizeof(double), shape);

  double *arrData  = arr->data;
  for (unsigned i=0; i<arr->size; ++i)
    arrData[i] = (i+0.5)*0.1;

  size_t shape2d[2] = { 4, 5 };
  gkyl_array_reshape(arr, 2, shape2d);

  TEST_CHECK( arr->rank == 2 );
  TEST_CHECK( arr->shape[0] == shape2d[0] );
  TEST_CHECK( arr->shape[1] == shape2d[1] );
  TEST_CHECK( arr->elemSz == sizeof(double) );
  TEST_CHECK( arr->size == 20 );
  TEST_CHECK( arr->ref_count.count == 1 );

  size_t shape1d[1] = { 45 };
  gkyl_array_reshape(arr, 1, shape1d);

  TEST_CHECK( arr->rank == 1 );
  TEST_CHECK( arr->shape[0] == shape1d[0] );
  TEST_CHECK( arr->elemSz == sizeof(double) );
  TEST_CHECK( arr->size == 45 );
  TEST_CHECK( arr->ref_count.count == 1 );

  arrData  = arr->data; // refetch data pointer
  
  for (unsigned i=0; i<20; ++i)
    TEST_CHECK( arrData[i] == (i+0.5)*0.1 );
  for (unsigned i=20; i<40; ++i)
    TEST_CHECK( arrData[i] == (i-20+0.5)*0.1 );
  for (unsigned i=40; i<45; ++i)
    TEST_CHECK( arrData[i] == (i-40+0.5)*0.1 );

  gkyl_array_release(arr);
}

void test_array_reshape_2()
{
  size_t shape[1] = {20};
  struct gkyl_array *arr = gkyl_array_new(1, sizeof(double), shape);

  double *arrData  = arr->data;
  for (unsigned i=0; i<arr->size; ++i)
    arrData[i] = (i+0.5)*0.1;

  size_t shape2d[2] = { 2, 2 };
  gkyl_array_reshape(arr, 2, shape2d);

  TEST_CHECK( arr->rank == 2 );
  TEST_CHECK( arr->shape[0] == shape2d[0] );
  TEST_CHECK( arr->shape[1] == shape2d[1] );
  TEST_CHECK( arr->elemSz == sizeof(double) );
  TEST_CHECK( arr->size == 4 );
  TEST_CHECK( arr->ref_count.count == 1 );

  arrData  = arr->data; // refetch data pointer
  
  for (unsigned i=0; i<4; ++i)
    TEST_CHECK( arrData[i] == (i+0.5)*0.1 );

  gkyl_array_release(arr);
}

void test_array_reshape_3()
{
  size_t shape[1] = {20};
  struct gkyl_array *arr = gkyl_array_new(1, sizeof(double), shape);

  double *arrData  = arr->data;
  for (unsigned i=0; i<arr->size; ++i)
    arrData[i] = (i+0.5)*0.1;

  size_t shape2d[2] = { 2, 2 };
  // create a new reshaped array from arr (don't modify arr)
  struct gkyl_array *brr = gkyl_array_reshape(gkyl_array_clone(arr),
    2, shape2d);

  TEST_CHECK( arr->rank == 1 );
  TEST_CHECK( arr->shape[0] == shape[0] );

  TEST_CHECK( brr->rank == 2 );
  TEST_CHECK( brr->shape[0] == shape2d[0] );
  TEST_CHECK( brr->shape[1] == shape2d[1] );
  TEST_CHECK( brr->elemSz == sizeof(double) );
  TEST_CHECK( brr->size == 4 );
  TEST_CHECK( brr->ref_count.count == 1 );

  double *brrData  = brr->data; // refetch data pointer
  
  for (unsigned i=0; i<4; ++i)
    TEST_CHECK( brrData[i] == (i+0.5)*0.1 );

  gkyl_array_release(arr);
  gkyl_array_release(brr);  
}

void test_array_fetch()
{
  size_t shape[1] = {20};
  struct gkyl_array *arr = gkyl_array_new(1, sizeof(double), shape);

  double *arrData  = arr->data;
  for (unsigned i=0; i<arr->size; ++i)
    arrData[i] = (i+0.5)*0.1;

  double *arrDataLh = gkyl_array_fetch(arr, 0);
  TEST_CHECK( arrDataLh[0] == 0.05 );

  double *arrDataUh = gkyl_array_fetch(arr, 10);
  TEST_CHECK( arrDataUh[0] == (10+0.5)*0.1);
 
  gkyl_array_release(arr);
}

void test_array_clear_double()
{
  size_t shape[] = {10};
  struct gkyl_array *a1 = gkyl_array_new(1, sizeof(double), shape);

  gkyl_array_clear(a1, 0.5);
  double *a1_d  = a1->data;  

  for (unsigned i=0; i<a1->size; ++i)
    TEST_CHECK( gkyl_compare(a1_d[i], 0.5, 1e-14) );

  gkyl_array_release(a1);
}

void test_array_clear_float()
{
  size_t shape[] = {10};
  struct gkyl_array *a1 = gkyl_array_new(1, sizeof(float), shape);

  gkyl_array_clear(a1, 0.5f);

  float *a1_d  = a1->data;  
  for (unsigned i=0; i<a1->size; ++i)
    TEST_CHECK( gkyl_compare(a1_d[i], 0.5, 1e-10) );

  gkyl_array_release(a1);
}

void test_array_accumulate_double()
{
  size_t shape[] = {10};
  struct gkyl_array *a1 = gkyl_array_new(1, sizeof(double), shape);
  struct gkyl_array *a2 = gkyl_array_new(1, sizeof(double), shape);

  double *a1_d  = a1->data, *a2_d = a2->data;
  for (unsigned i=0; i<a1->size; ++i) {
    a1_d[i] = i*1.0;
    a2_d[i] = i*0.1;
  }

  gkyl_array_accumulate(a1, 0.5, a2);

  for (unsigned i=0; i<a1->size; ++i)
    TEST_CHECK( gkyl_compare(a1_d[i], i*1.0+0.5*i*0.1, 1e-14) );

  gkyl_array_release(a1);
  gkyl_array_release(a2);
}

void test_array_accumulate_float()
{
  size_t shape[] = {10};
  struct gkyl_array *a1 = gkyl_array_new(1, sizeof(float), shape);
  struct gkyl_array *a2 = gkyl_array_new(1, sizeof(float), shape);

  float *a1_d  = a1->data, *a2_d = a2->data;
  for (unsigned i=0; i<a1->size; ++i) {
    a1_d[i] = i*1.0;
    a2_d[i] = i*0.1;
  }

  gkyl_array_accumulate(a1, 0.5f, a2);

  for (unsigned i=0; i<a1->size; ++i)
    TEST_CHECK( gkyl_compare(a1_d[i], i*1.0+0.5*i*0.1, 1e-10) );

  gkyl_array_release(a1);
  gkyl_array_release(a2);
}

void test_array_set_double()
{
  size_t shape[] = {10};
  struct gkyl_array *a1 = gkyl_array_new(1, sizeof(double), shape);
  struct gkyl_array *a2 = gkyl_array_new(1, sizeof(double), shape);

  double *a1_d  = a1->data, *a2_d = a2->data;
  for (unsigned i=0; i<a1->size; ++i) {
    a1_d[i] = i*1.0;
    a2_d[i] = i*0.1;
  }

  gkyl_array_set(a1, 0.5, a2);

  for (unsigned i=0; i<a1->size; ++i)
    TEST_CHECK( gkyl_compare(a1_d[i], 0.5*i*0.1, 1e-14) );

  gkyl_array_release(a1);
  gkyl_array_release(a2);
}

void test_array_set_float()
{
  size_t shape[] = {10};
  struct gkyl_array *a1 = gkyl_array_new(1, sizeof(float), shape);
  struct gkyl_array *a2 = gkyl_array_new(1, sizeof(float), shape);

  float *a1_d  = a1->data, *a2_d = a2->data;
  for (unsigned i=0; i<a1->size; ++i) {
    a1_d[i] = i*1.0;
    a2_d[i] = i*0.1;
  }

  gkyl_array_set(a1, 0.5f, a2);

  for (unsigned i=0; i<a1->size; ++i)
    TEST_CHECK( gkyl_compare(a1_d[i], 0.5*i*0.1, 1e-10) );

  gkyl_array_release(a1);
  gkyl_array_release(a2);
}

void test_array_scale_double()
{
  size_t shape[] = {10};
  struct gkyl_array *a1 = gkyl_array_new(1, sizeof(double), shape);

  double *a1_d  = a1->data;
  for (unsigned i=0; i<a1->size; ++i) {
    a1_d[i] = i*1.0;
  }

  gkyl_array_scale(a1, 0.25);

  for (unsigned i=0; i<a1->size; ++i)
    TEST_CHECK( gkyl_compare(a1_d[i], i*0.25, 1e-14) );

  gkyl_array_release(a1);
}

void test_array_opcombine()
{
  size_t shape[] = {10};
  struct gkyl_array *a1 = gkyl_array_new(1, sizeof(double), shape);
  struct gkyl_array *a2 = gkyl_array_new(1, sizeof(double), shape);

  double *a1_d  = a1->data, *a2_d = a2->data;
  for (unsigned i=0; i<a1->size; ++i) {
    a1_d[i] = i*1.0;
    a2_d[i] = i*0.1;
  }

  // a1 <- 0.25*(a1 + 0.5*a2)
  gkyl_array_scale(gkyl_array_set(a1, 0.5, a2), 0.25);

  for (unsigned i=0; i<a1->size; ++i)
    TEST_CHECK( gkyl_compare(a1_d[i], 0.25*0.5*i*0.1, 1e-14) );

  gkyl_array_release(a1);
  gkyl_array_release(a2);
}

void test_array_ops()
{
  size_t shape[] = {10};
  struct gkyl_array *a1 = gkyl_array_new(1, sizeof(double), shape);
  struct gkyl_array *a2 = gkyl_array_new(1, sizeof(double), shape);

  double *a1_d  = a1->data, *a2_d = a2->data;
  for (unsigned i=0; i<a1->size; ++i)
    a2_d[i] = i*0.1;

  gkyl_array_uniop("square", 0.0, a1, 1.0, a2);

  for (unsigned i=0; i<a1->size; ++i)
    TEST_CHECK( gkyl_compare(a1_d[i], (i*0.1)*(i*0.1), 1e-14) );

  gkyl_array_release(a1);
  gkyl_array_release(a2);  
}

TEST_LIST = {
  { "array_base", test_array_base },
  { "array_reshape", test_array_reshape },
  { "array_reshape_2", test_array_reshape_2 },
  { "array_reshape_3", test_array_reshape_3 },
  { "array_fetch", test_array_fetch },
  { "array_clear_double", test_array_clear_double },
  { "array_clear_float", test_array_clear_float },
  { "array_accumulate_double", test_array_accumulate_double },
  { "array_accumulate_float", test_array_accumulate_float },
  { "array_set_double", test_array_set_double },
  { "array_set_float", test_array_set_float },
  { "array_scale_double", test_array_scale_double },
  { "array_opcombine", test_array_opcombine },
  { "array_ops", test_array_ops },
  { NULL, NULL },
};
