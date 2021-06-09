#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_util.h>

void test_array_base()
{
  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, 1, 200);

  TEST_CHECK( arr->type = GKYL_DOUBLE );
  TEST_CHECK( arr->elemsz == sizeof(double) );
  TEST_CHECK( arr->ncomp == 1 );
  TEST_CHECK( arr->size == 20*10 );
  TEST_CHECK( arr->ref_count.count == 1 );

  TEST_CHECK( gkyl_array_is_cu_dev(arr) == false );

  double *arrData  = arr->data;
  for (unsigned i=0; i<arr->size; ++i)
    arrData[i] = (i+0.5)*0.1;

  // clone array
  struct gkyl_array *brr = gkyl_array_clone(arr);

  TEST_CHECK( brr->elemsz == sizeof(double) );
  TEST_CHECK( arr->ncomp == 1 );  
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
  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, 1, 20);

  double *arrData  = arr->data;
  for (unsigned i=0; i<arr->size; ++i)
    arrData[i] = (i+0.5)*0.1;

  double *arrDataLh = gkyl_array_fetch(arr, 0);
  TEST_CHECK( arrDataLh[0] == 0.05 );

  double *arrDataUh = gkyl_array_fetch(arr, 10);
  TEST_CHECK( arrDataUh[0] == (10+0.5)*0.1);
 
  gkyl_array_release(arr);
}

void test_array_clear()
{
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 1, 10);

  gkyl_array_clear(a1, 0.5);
  double *a1_d  = a1->data;  

  for (unsigned i=0; i<a1->size; ++i)
    TEST_CHECK( gkyl_compare(a1_d[i], 0.5, 1e-14) );

  gkyl_array_release(a1);
}

void test_array_clear_range()
{
  int shape[] = {10, 20};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 1, range.volume);
  gkyl_array_clear_range(a1, 0.5, range);

  double *a1_d = a1->data;
  for (unsigned i=0; i<a1->size; ++i)
    TEST_CHECK( gkyl_compare(a1_d[i], 0.5, 1e-14) );

  gkyl_array_release(a1);
}

void test_array_accumulate()
{
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 1, 10);
  struct gkyl_array *a2 = gkyl_array_new(GKYL_DOUBLE, 1, 10);

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

void test_array_accumulate_range()
{
  int shape[] = {10, 20};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 8, range.volume);
  struct gkyl_array *a2 = gkyl_array_new(GKYL_DOUBLE, 3, range.volume);

  // test a1 = a1 + 0.5*a2
  gkyl_array_clear(a1, 0.5);
  gkyl_array_clear(a2, 1.5);

  gkyl_array_accumulate_range(a1, 0.5, a2, range);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);
    double *a1d = gkyl_array_fetch(a1, loc);
    for (int i=0; i<3; ++i)
      TEST_CHECK( a1d[i] == 0.5 + 0.5*1.5 );
    for (int i=3; i<8; ++i)
      TEST_CHECK( a1d[i] == 0.5);
  }

  // test a2 = a2 + 0.5*a
  gkyl_array_clear(a1, 0.5);
  gkyl_array_clear(a2, 1.5);

  gkyl_array_accumulate_range(a2, 0.5, a1, range);

  gkyl_range_iter_init(&iter, &range);
  
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);
    double *a2d = gkyl_array_fetch(a2, loc);
    for (int i=0; i<3; ++i)
      TEST_CHECK( a2d[i] == 1.5 + 0.5*0.5 );
  }  

  gkyl_array_release(a1);
  gkyl_array_release(a2);
}

void test_array_combine()
{
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 1, 10);
  struct gkyl_array *a2 = gkyl_array_new(GKYL_DOUBLE, 1, 10);
  struct gkyl_array *b = gkyl_array_new(GKYL_DOUBLE, 1, 10);

  double *a1_d  = a1->data, *a2_d = a2->data, *b_d = b->data;
  for (unsigned i=0; i<a1->size; ++i) {
    a1_d[i] = i*1.0;
    a2_d[i] = i*0.1;
    b_d[i] = 10.5;
  }

  // b = 0.5*a1 + 2.5*a2
  gkyl_array_accumulate(gkyl_array_set(b, 0.5, a1), 2.5, a2);

  for (unsigned i=0; i<a1->size; ++i)
    TEST_CHECK( gkyl_compare(b_d[i], 0.5*i*1.0+2.5*i*0.1, 1e-14) );

  gkyl_array_release(a1);
  gkyl_array_release(a2);
  gkyl_array_release(b);
}

void test_array_set()
{
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 1, 10);
  struct gkyl_array *a2 = gkyl_array_new(GKYL_DOUBLE, 1, 10);

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

void test_array_set_range()
{
  int shape[] = {10, 20};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 8, range.volume);
  struct gkyl_array *a2 = gkyl_array_new(GKYL_DOUBLE, 3, range.volume);

  // test a1 = 0.5*a2
  gkyl_array_clear_range(a1, 0.5, range);
  gkyl_array_clear_range(a2, 1.5, range);

  gkyl_array_set_range(a1, 0.5, a2, range);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);
    double *a1d = gkyl_array_fetch(a1, loc);
    for (int i=0; i<3; ++i)
      TEST_CHECK( a1d[i] == 0.5*1.5 );
    for (int i=3; i<8; ++i)
      TEST_CHECK( a1d[i] == 0.5);
  }

  // test a2 = 0.5*a1
  gkyl_array_clear(a1, 0.5);
  gkyl_array_clear(a2, 1.5);

  gkyl_array_set_range(a2, 0.5, a1, range);

  gkyl_range_iter_init(&iter, &range);
  
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);
    double *a2d = gkyl_array_fetch(a2, loc);
    for (int i=0; i<3; ++i)
      TEST_CHECK( a2d[i] == 0.5*0.5 );
  }  

  gkyl_array_release(a1);
  gkyl_array_release(a2);
}

void test_array_scale()
{
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 1, 10);

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
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 1, 10);  
  struct gkyl_array *a2 = gkyl_array_new(GKYL_DOUBLE, 1, 10);

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

void test_array_ops_comp() // more than 1 "component" in array
{
  int nc = 5; // number of "components"
  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, nc, 10);

  for (unsigned i=0; i<arr->size; ++i) {
    double *d = gkyl_array_fetch(arr, i);
    for (int k=0; k<nc; ++k)
      d[k] = i*1.0;
  }

  gkyl_array_clear(arr, 1.5f);

  for (unsigned i=0; i<arr->size; ++i) {
    const double *d = gkyl_array_fetch(arr, i);
    for (int k=0; k<nc; ++k)
      TEST_CHECK( gkyl_compare(d[k], 1.5f, 1e-10) );
  }

  gkyl_array_release(arr);
}

void test_array_copy()
{
  int shape[] = {10, 20};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, 1, range.volume);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    double *d = gkyl_array_fetch(arr, gkyl_range_idx(&range, iter.idx));
    d[0] = iter.idx[0] + 10.5*iter.idx[1];
  }

  int lower[] = {1, 1}, upper[] = {5, 10};
  struct gkyl_range sub_range;
  gkyl_sub_range_init(&sub_range, &range, lower, upper);

  double *buff = gkyl_malloc(sizeof(double)*sub_range.volume);
  gkyl_array_copy_to_buffer(buff, arr, sub_range);

  long count = 0;
  gkyl_range_iter_init(&iter, &sub_range);
  while (gkyl_range_iter_next(&iter))
    TEST_CHECK( buff[count++] == iter.idx[0] + 10.5*iter.idx[1] );

  gkyl_array_clear(arr, 0.0);
  // copy back from buffer
  gkyl_array_copy_from_buffer(arr, buff, sub_range);

  gkyl_range_iter_init(&iter, &sub_range);
  while (gkyl_range_iter_next(&iter)) {
    double *d = gkyl_array_fetch(arr, gkyl_range_idx(&sub_range, iter.idx));
    TEST_CHECK( d[0]  == iter.idx[0] + 10.5*iter.idx[1] );
  }

  gkyl_array_release(arr);
  gkyl_free(buff);
}

void test_array_copy_split()
{
  int shape[] = {10, 20};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, 1, range.volume);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    double *d = gkyl_array_fetch(arr, gkyl_range_idx(&range, iter.idx));
    d[0] = iter.idx[0] + 10.5*iter.idx[1];
  }

  int lower[] = {1, 1}, upper[] = {5, 10};
  struct gkyl_range sub_range;
  gkyl_sub_range_init(&sub_range, &range, lower, upper);

  double *buff = gkyl_malloc(sizeof(double[sub_range.volume]));

  int num_split = 3;
  
  long loc = 0;
  for (int tid=0; tid<num_split; ++tid) {
    gkyl_range_set_split(&sub_range, num_split, tid);
    gkyl_array_copy_to_buffer(buff+loc, arr, sub_range);
    
    loc += gkyl_range_split_len(&sub_range);
  }

  long count = 0;
  gkyl_range_set_split(&sub_range, 1, 0);
  gkyl_range_iter_init(&iter, &sub_range);
  while (gkyl_range_iter_next(&iter))
    TEST_CHECK( buff[count++] == iter.idx[0] + 10.5*iter.idx[1] );

  gkyl_array_clear(arr, 0.0);

  loc = 0;
  for (int tid=0; tid<num_split; ++tid) {
    // copy back from buffer
    gkyl_range_set_split(&sub_range, num_split, tid);
    gkyl_array_copy_from_buffer(arr, buff+loc, sub_range);
    loc += gkyl_range_split_len(&sub_range);
  }

  gkyl_range_set_split(&sub_range, 1, 0); // reset

  gkyl_range_iter_init(&iter, &sub_range);
  while (gkyl_range_iter_next(&iter)) {
    double *d = gkyl_array_fetch(arr, gkyl_range_idx(&sub_range, iter.idx));
    TEST_CHECK( d[0]  == iter.idx[0] + 10.5*iter.idx[1] );
  }

  gkyl_array_release(arr);
  gkyl_free(buff);
}

void test_non_numeric()
{
  struct euler { double rho, u, E; };
  struct gkyl_array *arr = gkyl_array_new(GKYL_USER, sizeof(struct euler), 10);

  for (unsigned i=0; i<arr->size; ++i) {
    struct euler *e = gkyl_array_fetch(arr, i);
    e->rho = 1.0; e->u = 0.0; e->E = 100.5;
  }
  
  struct gkyl_array *brr = gkyl_array_new(GKYL_USER, sizeof(struct euler), 10);
  gkyl_array_copy(brr, arr);

  for (unsigned i=0; i<arr->size; ++i) {
    struct euler *e = gkyl_array_fetch(brr, i);
    TEST_CHECK( e->rho == 1.0 );
    TEST_CHECK( e->u == 0.0 );
    TEST_CHECK( e->E == 100.5 );
  }  

  gkyl_array_release(arr);
  gkyl_array_release(brr);
}

void test_reduce()
{
  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, 3, 200);

  for (size_t i=0; i<arr->size; ++i) {
    double *d = gkyl_array_fetch(arr, i);
    for (size_t c=0; c<3; ++c)
      d[c] = 0.5*i + 0.1*c;
  }
  
  double amin, amax;
  gkyl_array_reduce(arr, GKYL_MIN, &amin);
  gkyl_array_reduce(arr, GKYL_MAX, &amax);

  TEST_CHECK( amin == 0.0 );
  TEST_CHECK( amax == 0.5*199 + 0.1*2 );

  gkyl_array_release(arr);
}

void test_reduce_range()
{
  int shape[] = {10, 20};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, 3, range.volume);

  int count = -1000;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);
    double *d = gkyl_array_fetch(arr, loc);
    d[0] = count+0.5;
    d[1] = count+1.5;
    d[2] = count+2.5;

    count += 1;
  }

  double amin[3], amax[3];
  gkyl_array_reduce_range(amin, arr, GKYL_MIN, range);

  TEST_CHECK( amin[0] == -999.5 );
  TEST_CHECK( amin[1] == -998.5 );
  TEST_CHECK( amin[2] == -997.5 );

  gkyl_array_reduce_range(amax, arr, GKYL_MAX, range);
  
  TEST_CHECK( amax[0] == -800.5 );
  TEST_CHECK( amax[1] == -799.5 );
  TEST_CHECK( amax[2] == -798.5 );

  gkyl_array_release(arr);
}

// CUDA specific tests
#ifdef GKYL_HAVE_CUDA

/* Function signatures of kernel calls */
int cu_array_test_and_flip_sign( struct gkyl_array *arr);

void test_cu_array_base()
{
  struct gkyl_array *arr_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, 200);

  TEST_CHECK( arr_cu->type = GKYL_DOUBLE );
  TEST_CHECK( arr_cu->elemsz == sizeof(double) );
  TEST_CHECK( arr_cu->ncomp == 1 );
  TEST_CHECK( arr_cu->size == 20*10 );
  TEST_CHECK( arr_cu->ref_count.count == 1 );

  TEST_CHECK( gkyl_array_is_cu_dev(arr_cu) == true );

  // create host array and initialize it
  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, 1, 200);

  double *arrData  = arr->data;
  for (unsigned i=0; i<arr->size; ++i)
    arrData[i] = (i+0.5)*0.1;  

  // copy to device
  gkyl_array_copy(arr_cu, arr);

  // reset host array
  for (unsigned i=0; i<arr->size; ++i)
    arrData[i] = 0.0;

  // copy from device and check if things are ok
  gkyl_array_copy(arr, arr_cu);

  for (unsigned i=0; i<arr->size; ++i)
    TEST_CHECK( arrData[i] == (i+0.5)*0.1 );

  gkyl_array_release(arr);
  gkyl_array_release(arr_cu);
}

void test_cu_array_dev_kernel()
{
  // create a host array struct containing device data
  struct gkyl_array *arr_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, 20);

  TEST_CHECK( arr_cu->type = GKYL_DOUBLE );
  TEST_CHECK( arr_cu->elemsz == sizeof(double) );
  TEST_CHECK( arr_cu->ncomp == 1 );
  TEST_CHECK( arr_cu->size == 20 );
  TEST_CHECK( arr_cu->ref_count.count == 1 );
  TEST_CHECK( gkyl_array_is_cu_dev(arr_cu) == true );

  // create host array and initialize it
  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, 1, 20);

  double *arrData  = arr->data;
  for (unsigned i=0; i<arr->size; ++i)
    arrData[i] = (i+0.5)*0.1;  

  // copy arr data to device data in arr_cu
  gkyl_array_copy(arr_cu, arr);

  // create new struct arr_cu_cl that is a clone of arr_cu
  struct gkyl_array *arr_cu_cl = gkyl_array_clone(arr_cu);

  // check arr_cu on device and flip sign
  int nfail = cu_array_test_and_flip_sign(arr_cu->on_device);
  TEST_CHECK( nfail == 0 );

  // restore arr_cu by copying from arr_cu_cl, and test again
  gkyl_array_copy(arr_cu, arr_cu_cl);
  nfail = cu_array_test_and_flip_sign(arr_cu->on_device);
  TEST_CHECK( nfail == 0 );

  // check arr_cu_cl on device and flip sign
  nfail = cu_array_test_and_flip_sign(arr_cu_cl->on_device);
  TEST_CHECK( nfail == 0 );

  // copy arr_cu back to host and check
  gkyl_array_copy(arr, arr_cu);
  for (unsigned i=0; i<arr->size; ++i)
    TEST_CHECK( arrData[i] == -(i+0.5)*0.1);  

  // copy arr_cu_cl back to host and check
  // zero out arr first (no cheating)
  gkyl_array_clear(arr, 0.);
  gkyl_array_copy(arr, arr_cu_cl);
  for (unsigned i=0; i<arr->size; ++i)
    TEST_CHECK( arrData[i] == -(i+0.5)*0.1);  

  // release all data
  gkyl_array_release(arr);
  gkyl_array_release(arr_cu);  
  gkyl_array_release(arr_cu_cl);  
}

void test_cu_array_clear()
{
  // create host array 
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 1, 10);
  // make device copy
  struct gkyl_array *a1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, 10);
  gkyl_array_copy(a1_cu, a1);

  // clear array on device
  gkyl_array_clear_cu(a1_cu, 0.5);
  // copy from device and check if things are ok
  gkyl_array_copy(a1, a1_cu);
  double *a1_d  = a1->data;  
  for (unsigned i=0; i<a1->size; ++i)
    TEST_CHECK( gkyl_compare(a1_d[i], 0.5, 1e-14) );

  gkyl_array_release(a1);
  gkyl_array_release(a1_cu);
}

void test_cu_array_clear_range()
{
  int shape[] = {10, 20};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);

  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 1, range.volume);
  // make device copy of array
  struct gkyl_array *a1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, range.volume);
  gkyl_array_copy(a1_cu, a1);
  gkyl_array_clear_range_cu(a1_cu, 0.5, range);

  // copy from device and check if things are ok
  gkyl_array_copy(a1, a1_cu);
  double *a1_d  = a1->data;  
  for (unsigned i=0; i<a1->size; ++i)
    TEST_CHECK( gkyl_compare(a1_d[i], 0.5, 1e-14) );

  gkyl_array_release(a1);
  gkyl_array_release(a1_cu);
}

void test_cu_array_accumulate()
{
  // create host arrays 
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 1, 10);
  struct gkyl_array *a2 = gkyl_array_new(GKYL_DOUBLE, 1, 10);
  // make device copies
  struct gkyl_array *a1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, 10);
  struct gkyl_array *a2_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, 10);

  // initialize data
  double *a1_d  = a1->data, *a2_d = a2->data;
  for (unsigned i=0; i<a1->size; ++i) {
    a1_d[i] = i*1.0;
    a2_d[i] = i*0.1;
  }

  // copy initialized arrays to device
  gkyl_array_copy(a1_cu, a1);
  gkyl_array_copy(a2_cu, a2);

  gkyl_array_accumulate_cu(a1_cu, 0.5, a2_cu);

  // copy from device and check if things are ok
  gkyl_array_copy(a1, a1_cu);
  for (unsigned i=0; i<a1->size; ++i)
    TEST_CHECK( gkyl_compare(a1_d[i], i*1.0+0.5*i*0.1, 1e-14) );

  gkyl_array_release(a1);
  gkyl_array_release(a2);
  gkyl_array_release(a1_cu);
  gkyl_array_release(a2_cu);
}

void test_cu_array_accumulate_range()
{
  int shape[] = {10, 20};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 8, range.volume);
  struct gkyl_array *a2 = gkyl_array_new(GKYL_DOUBLE, 3, range.volume);

  // make device copies of arrays
  struct gkyl_array *a1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 8, range.volume);
  struct gkyl_array *a2_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3, range.volume);

  // copy host arrays to device
  gkyl_array_copy(a1_cu, a1);
  gkyl_array_copy(a2_cu, a2);

  // test a1 = a1 + 0.5*a2
  gkyl_array_clear_cu(a1_cu, 0.5);
  gkyl_array_clear_cu(a2_cu, 1.5);

  gkyl_array_accumulate_range_cu(a1_cu, 0.5, a2_cu, range);

 // copy from device and check if things are ok
  gkyl_array_copy(a1, a1_cu);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);
    double *a1d = gkyl_array_fetch(a1, loc);
    for (int i=0; i<3; ++i)
      TEST_CHECK( a1d[i] == 0.5 + 0.5*1.5 );
    for (int i=3; i<8; ++i)
      TEST_CHECK( a1d[i] == 0.5);
  }

  // test a2 = a2 + 0.5*a
  gkyl_array_clear_cu(a1_cu, 0.5);
  gkyl_array_clear_cu(a2_cu, 1.5);

  gkyl_array_accumulate_range_cu(a2_cu, 0.5, a1_cu, range);

 // copy from device and check if things are ok
  gkyl_array_copy(a2, a2_cu);
  gkyl_range_iter_init(&iter, &range);

  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);
    double *a2d = gkyl_array_fetch(a2, loc);
    for (int i=0; i<3; ++i)
      TEST_CHECK( a2d[i] == 1.5 + 0.5*0.5 );
  }  

  gkyl_array_release(a1);
  gkyl_array_release(a2);
  gkyl_array_release(a1_cu);
  gkyl_array_release(a2_cu);
}

void test_cu_array_combine()
{
  // create host arrays 
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 1, 10);
  struct gkyl_array *a2 = gkyl_array_new(GKYL_DOUBLE, 1, 10);
  struct gkyl_array *b = gkyl_array_new(GKYL_DOUBLE, 1, 10);

  // make device copies
  struct gkyl_array *a1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, 10);
  struct gkyl_array *a2_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, 10);
  struct gkyl_array *b_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, 10);

  double *a1_d  = a1->data, *a2_d = a2->data, *b_d = b->data;
  for (unsigned i=0; i<a1->size; ++i) {
    a1_d[i] = i*1.0;
    a2_d[i] = i*0.1;
    b_d[i] = 10.5;
  }

  // copy initialized arrays to device
  gkyl_array_copy(a1_cu, a1);
  gkyl_array_copy(a2_cu, a2);
  gkyl_array_copy(b_cu, b);

  // b = 0.5*a1 + 2.5*a2
  gkyl_array_set_cu(b_cu, 0.5, a1_cu);
  gkyl_array_accumulate_cu(b_cu, 2.5, a2_cu);

  // copy from device and check if things are ok
  gkyl_array_copy(b, b_cu);
  for (unsigned i=0; i<a1->size; ++i)
    TEST_CHECK( gkyl_compare(b_d[i], 0.5*i*1.0+2.5*i*0.1, 1e-14) );

  gkyl_array_release(a1);
  gkyl_array_release(a2);
  gkyl_array_release(b);
  gkyl_array_release(a1_cu);
  gkyl_array_release(a2_cu);
  gkyl_array_release(b_cu);
}

void test_cu_array_set()
{
  // create host arrays 
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 1, 10);
  struct gkyl_array *a2 = gkyl_array_new(GKYL_DOUBLE, 1, 10);
  // make device copies
  struct gkyl_array *a1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, 10);
  struct gkyl_array *a2_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, 10);

  // initialize data
  double *a1_d  = a1->data, *a2_d = a2->data;
  for (unsigned i=0; i<a1->size; ++i) {
    a1_d[i] = i*1.0;
    a2_d[i] = i*0.1;
  }

  // copy initialized arrays to device
  gkyl_array_copy(a1_cu, a1);
  gkyl_array_copy(a2_cu, a2);

  gkyl_array_set_cu(a1_cu, 0.5, a2_cu);

  // copy from device and check if things are ok
  gkyl_array_copy(a1, a1_cu);
  for (unsigned i=0; i<a1->size; ++i)
    TEST_CHECK( gkyl_compare(a1_d[i], 0.5*i*0.1, 1e-14) );

  gkyl_array_release(a1);
  gkyl_array_release(a2);
  gkyl_array_release(a1_cu);
  gkyl_array_release(a2_cu);
}

void test_cu_array_set_range()
{
  int shape[] = {10, 20};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 8, range.volume);
  struct gkyl_array *a2 = gkyl_array_new(GKYL_DOUBLE, 3, range.volume);

  // make device copies of arrays
  struct gkyl_array *a1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 8, range.volume);
  struct gkyl_array *a2_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3, range.volume);

  // copy host arrays to device
  gkyl_array_copy(a1_cu, a1);
  gkyl_array_copy(a2_cu, a2);

  // test a1 = 0.5*a2
  gkyl_array_clear_cu(a1_cu, 0.5);
  gkyl_array_clear_cu(a2_cu, 1.5);

  gkyl_array_set_range_cu(a1_cu, 0.5, a2_cu, range);

 // copy from device and check if things are ok
  gkyl_array_copy(a1, a1_cu);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);
    double *a1d = gkyl_array_fetch(a1, loc);
    for (int i=0; i<3; ++i)
      TEST_CHECK( a1d[i] == 0.5*1.5 );
    for (int i=3; i<8; ++i)
      TEST_CHECK( a1d[i] == 0.5);
  }

  // test a2 = 0.5*a1
  gkyl_array_clear_cu(a1_cu, 0.5);
  gkyl_array_clear_cu(a2_cu, 1.5);

  gkyl_array_set_range_cu(a2_cu, 0.5, a1_cu, range);

  // copy from device and check if things are ok
  gkyl_array_copy(a2, a2_cu);
  gkyl_range_iter_init(&iter, &range);
  
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);
    double *a2d = gkyl_array_fetch(a2, loc);
    for (int i=0; i<3; ++i)
      TEST_CHECK( a2d[i] == 0.5*0.5 );
  }  

  gkyl_array_release(a1);
  gkyl_array_release(a2);
  gkyl_array_release(a1_cu);
  gkyl_array_release(a2_cu);
}

void test_cu_array_scale()
{
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 1, 10);
  // make device copies
  struct gkyl_array *a1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, 10);

  // initialize data
  double *a1_d  = a1->data;
  for (unsigned i=0; i<a1->size; ++i) {
    a1_d[i] = i*1.0;
  }

  // copy host arrays to device
  gkyl_array_copy(a1_cu, a1);

  gkyl_array_scale_cu(a1_cu, 0.25);

 // copy from device and check if things are ok
  gkyl_array_copy(a1, a1_cu);
  for (unsigned i=0; i<a1->size; ++i)
    TEST_CHECK( gkyl_compare(a1_d[i], i*0.25, 1e-14) );

  gkyl_array_release(a1);
  gkyl_array_release(a1_cu);
}

void test_cu_array_copy()
{
  int shape[] = {10, 20};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, 1, range.volume);
  // make device copies of arrays
  struct gkyl_array *arr_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, range.volume);

  // initialize array that will be copied
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    double *d = gkyl_array_fetch(arr, gkyl_range_idx(&range, iter.idx));
    d[0] = iter.idx[0] + 10.5*iter.idx[1];
  }

  // copy host array to device
  gkyl_array_copy(arr_cu, arr);

  int lower[] = {1, 1}, upper[] = {5, 10};
  struct gkyl_range sub_range;
  gkyl_sub_range_init(&sub_range, &range, lower, upper);

  double *buff_cu = gkyl_cu_malloc(sizeof(double)*sub_range.volume);
  gkyl_array_copy_to_buffer_cu(buff_cu, arr_cu, sub_range);

  gkyl_array_clear_cu(arr, 0.0);
  // copy back from buffer
  gkyl_array_copy_from_buffer_cu(arr_cu, buff_cu, sub_range);

  // copy from device and check if things are ok
  gkyl_array_copy(arr, arr_cu);
  gkyl_range_iter_init(&iter, &sub_range);
  while (gkyl_range_iter_next(&iter)) {
    double *d = gkyl_array_fetch(arr, gkyl_range_idx(&sub_range, iter.idx));
    TEST_CHECK( d[0]  == iter.idx[0] + 10.5*iter.idx[1] );
  }

  gkyl_array_release(arr);
  gkyl_array_release(arr_cu);
  gkyl_cu_free(buff_cu);
}

#endif

TEST_LIST = {
  { "array_base", test_array_base },
  { "array_fetch", test_array_fetch },
  { "array_clear", test_array_clear },
  { "array_clear_range", test_array_clear_range },
  { "array_accumulate", test_array_accumulate },
  { "array_accumulate_range", test_array_accumulate_range },
  { "array_combine", test_array_combine },
  { "array_set", test_array_set },
  { "array_set_range", test_array_set_range },
  { "array_scale", test_array_scale },
  { "array_opcombine", test_array_opcombine },
  { "array_ops_comp", test_array_ops_comp },
  { "array_copy", test_array_copy },
  { "array_copy_split", test_array_copy_split },
  { "non_numeric", test_non_numeric },
  { "reduce", test_reduce },
  { "reduce_range", test_reduce_range },
#ifdef GKYL_HAVE_CUDA
  { "cu_array_base", test_cu_array_base },
  { "cu_array_clear", test_cu_array_clear},
  { "cu_array_clear_range", test_cu_array_clear_range},
  { "cu_array_accumulate", test_cu_array_accumulate},
  { "cu_array_accumulate_range", test_cu_array_accumulate_range},
  { "cu_array_combine", test_cu_array_combine},
  { "cu_array_set", test_cu_array_set },
  { "cu_array_set_range", test_cu_array_set_range },
  { "cu_array_scale", test_cu_array_scale },
  { "cu_array_copy", test_cu_array_copy },
  { "cu_array_dev_kernel", test_cu_array_dev_kernel },
#endif
  { NULL, NULL },
};
