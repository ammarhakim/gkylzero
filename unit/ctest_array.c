#include <acutest.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_util.h>

void test_array_base()
{
  struct gkyl_array *arr = gkyl_array_new(sizeof(double), 200);

  TEST_CHECK( arr->elemSz == sizeof(double) );
  TEST_CHECK( arr->size == 20*10 );
  TEST_CHECK( arr->ref_count.count == 1 );

  double *arrData  = arr->data;
  for (unsigned i=0; i<arr->size; ++i)
    arrData[i] = (i+0.5)*0.1;

  // clone array
  struct gkyl_array *brr = gkyl_array_clone(arr);

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

void test_array_fetch()
{
  struct gkyl_array *arr = gkyl_array_new(sizeof(double), 20);

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
  struct gkyl_array *a1 = gkyl_array_new(sizeof(double), 10);

  gkyl_array_clear(a1, 0.5);
  double *a1_d  = a1->data;  

  for (unsigned i=0; i<a1->size; ++i)
    TEST_CHECK( gkyl_compare(a1_d[i], 0.5, 1e-14) );

  gkyl_array_release(a1);
}

void test_array_clear_float()
{
  struct gkyl_array *a1 = gkyl_array_new(sizeof(float), 10);

  gkyl_array_clear(a1, 0.5f);

  float *a1_d  = a1->data;  
  for (unsigned i=0; i<a1->size; ++i)
    TEST_CHECK( gkyl_compare(a1_d[i], 0.5, 1e-10) );

  gkyl_array_release(a1);
}

void test_array_accumulate_double()
{
  struct gkyl_array *a1 = gkyl_array_new(sizeof(double), 10);
  struct gkyl_array *a2 = gkyl_array_new(sizeof(double), 10);

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
  struct gkyl_array *a1 = gkyl_array_new(sizeof(float), 10);
  struct gkyl_array *a2 = gkyl_array_new(sizeof(float), 10);

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
  struct gkyl_array *a1 = gkyl_array_new(sizeof(double), 10);
  struct gkyl_array *a2 = gkyl_array_new(sizeof(double), 10);

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
  struct gkyl_array *a1 = gkyl_array_new(sizeof(float), 10);
  struct gkyl_array *a2 = gkyl_array_new(sizeof(float), 10);

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
  struct gkyl_array *a1 = gkyl_array_new(sizeof(double), 10);

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
  struct gkyl_array *a1 = gkyl_array_new(sizeof(double), 10);  
  struct gkyl_array *a2 = gkyl_array_new(sizeof(double), 10);

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
  struct gkyl_array *a1 = gkyl_array_new(sizeof(double), 10);
  struct gkyl_array *a2 = gkyl_array_new(sizeof(double), 10);

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
