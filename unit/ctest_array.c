#include <acutest.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_util.h>
#include <gkyl_alloc.h>

void test_array_base()
{
  struct gkyl_array *arr = gkyl_array_new(sizeof(double), 200);

  TEST_CHECK( arr->elemsz == sizeof(double) );
  TEST_CHECK( arr->size == 20*10 );
  TEST_CHECK( arr->ref_count.count == 1 );

  double *arrData  = arr->data;
  for (unsigned i=0; i<arr->size; ++i)
    arrData[i] = (i+0.5)*0.1;

  // clone array
  struct gkyl_array *brr = gkyl_array_clone(arr);

  TEST_CHECK( brr->elemsz == sizeof(double) );
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

void test_array_clear()
{
  struct gkyl_array *a1 = gkyl_array_new(sizeof(double), 10);

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
  
  struct gkyl_array *a1 = gkyl_array_new(sizeof(double[1]), range.volume);
  gkyl_array_clear_range(a1, 0.5, &range);

  double *a1_d = a1->data;
  for (unsigned i=0; i<a1->size; ++i)
    TEST_CHECK( gkyl_compare(a1_d[i], 0.5, 1e-14) );

  gkyl_array_release(a1);
}

void test_array_accumulate()
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

void test_array_accumulate_range()
{
  int shape[] = {10, 20};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  struct gkyl_array *a1 = gkyl_array_new(sizeof(double[8]), range.volume);
  struct gkyl_array *a2 = gkyl_array_new(sizeof(double[3]), range.volume);

  // test a1 = a1 + 0.5*a2
  gkyl_array_clear(a1, 0.5);
  gkyl_array_clear(a2, 1.5);

  gkyl_array_accumulate_range(a1, 0.5, a2, &range);

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

  gkyl_array_accumulate_range(a2, 0.5, a1, &range);

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
  struct gkyl_array *a1 = gkyl_array_new(sizeof(double), 10);
  struct gkyl_array *a2 = gkyl_array_new(sizeof(double), 10);
  struct gkyl_array *b = gkyl_array_new(sizeof(double), 10);

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

void test_array_set_range()
{
  int shape[] = {10, 20};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  struct gkyl_array *a1 = gkyl_array_new(sizeof(double[8]), range.volume);
  struct gkyl_array *a2 = gkyl_array_new(sizeof(double[3]), range.volume);

  // test a1 = 0.5*a2
  gkyl_array_clear_range(a1, 0.5, &range);
  gkyl_array_clear_range(a2, 1.5, &range);

  gkyl_array_set_range(a1, 0.5, a2, &range);

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

  // test a2 = a2 + 0.5*a
  gkyl_array_clear(a1, 0.5);
  gkyl_array_clear(a2, 1.5);

  gkyl_array_accumulate_range(a2, 0.5, a1, &range);

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

void test_array_scale()
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

void test_array_ops_comp() // more than 1 "component" in array
{
  int nc = 5; // number of "components"
  struct gkyl_array *arr = gkyl_array_new(sizeof(double)*nc, 10);

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
  
  struct gkyl_array *arr = gkyl_array_new(sizeof(double), range.volume);

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
  gkyl_array_copy_to_buffer(buff, arr, &sub_range);

  long count = 0;
  gkyl_range_iter_init(&iter, &sub_range);
  while (gkyl_range_iter_next(&iter))
    TEST_CHECK( buff[count++] == iter.idx[0] + 10.5*iter.idx[1] );

  gkyl_array_clear(arr, 0.0);
  // copy back from buffer
  gkyl_array_copy_from_buffer(arr, buff, &sub_range);

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
  
  struct gkyl_array *arr = gkyl_array_new(sizeof(double), range.volume);

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
    gkyl_array_copy_to_buffer(buff+loc, arr, &sub_range);
    
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
    gkyl_array_copy_from_buffer(arr, buff+loc, &sub_range);
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
  struct gkyl_array *arr = gkyl_array_new(sizeof(struct euler), 10);

  for (unsigned i=0; i<arr->size; ++i) {
    struct euler *e = gkyl_array_fetch(arr, i);
    e->rho = 1.0; e->u = 0.0; e->E = 100.5;
  }
  
  struct gkyl_array *brr = gkyl_array_new(sizeof(struct euler), 10);
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
  struct gkyl_array *arr = gkyl_array_new(sizeof(double[3]), 200);

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
  
  struct gkyl_array *arr = gkyl_array_new(sizeof(double[3]), range.volume);

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
  gkyl_array_reduce_range(amin, arr, GKYL_MIN, &range);

  TEST_CHECK( amin[0] == -999.5 );
  TEST_CHECK( amin[1] == -998.5 );
  TEST_CHECK( amin[2] == -997.5 );

  gkyl_array_reduce_range(amax, arr, GKYL_MAX, &range);
  
  TEST_CHECK( amax[0] == -800.5 );
  TEST_CHECK( amax[1] == -799.5 );
  TEST_CHECK( amax[2] == -798.5 );

  gkyl_array_release(arr);
}

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
  { NULL, NULL },
};
