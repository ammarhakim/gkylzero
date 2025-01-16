#include <acutest.h>
#include <mpack.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_array_rio_format_desc.h>
#include <gkyl_array_rio_priv.h>
#include <gkyl_elem_type_priv.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_util.h>

void test_array_0()
{
  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, 1, 200);
  gkyl_array_release(arr);
}

void test_array_base()
{
  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, 1, 200);

  TEST_CHECK( gkyl_array_is_using_buffer(arr) == false );

  TEST_CHECK( arr->type = GKYL_DOUBLE );
  TEST_CHECK( arr->elemsz == sizeof(double) );
  TEST_CHECK( arr->ncomp == 1 );
  TEST_CHECK( arr->size == 20*10 );
  TEST_CHECK( arr->ref_count.count == 1 );

  TEST_CHECK( arr->on_dev == arr );

  TEST_CHECK( gkyl_array_is_cu_dev(arr) == false );

  double *arrData  = arr->data;
  for (unsigned i=0; i<arr->size; ++i){
    TEST_CHECK( arrData[i] == 0. );
    arrData[i] = (i+0.5)*0.1;
  }

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

  // acquire pointer
  struct gkyl_array *crr = gkyl_array_acquire(arr);

  TEST_CHECK( crr->ref_count.count == 2 );
  TEST_CHECK( arr->ref_count.count == 2 );

  struct gkyl_array *drr = gkyl_array_acquire(crr);

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
  gkyl_array_clear_range(a1, 0.5, &range);

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

void test_array_accumulate_offset()
{
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 2, 10);
  struct gkyl_array *a2 = gkyl_array_new(GKYL_DOUBLE, 3*a1->ncomp, 10);

  double *a1_d  = a1->data, *a2_d = a2->data;

  // test a1 = 0.1*a2[a1->ncomp]
  for (unsigned i=0; i<a1->size; ++i)
    for (unsigned j=0; j<a1->ncomp; ++j)
      a1_d[i*a1->ncomp+j] = i*1.0+j;

  for (unsigned i=0; i<a2->size; ++i)
    for (unsigned j=0; j<a2->ncomp/a1->ncomp; ++j)
      for (unsigned k=0; k<a1->ncomp; ++k)
        a2_d[i*a2->ncomp+j*a1->ncomp+k] = i*0.1+k;

  gkyl_array_accumulate_offset(a1, 0.5, a2, 1*a1->ncomp);

  for (unsigned i=0; i<a1->size; ++i)
    for (unsigned j=0; j<a1->ncomp; ++j)
      TEST_CHECK( gkyl_compare(a1_d[i*a1->ncomp+j], i*1.0+j+0.5*(i*0.1+j), 1e-14) );

  // test a2[a1->ncomp] = 0.1*a1
  for (unsigned i=0; i<a1->size; ++i)
    for (unsigned j=0; j<a1->ncomp; ++j)
      a1_d[i*a1->ncomp+j] = i*1.0+j;

  gkyl_array_accumulate_offset(a2, 0.5, a1, 1*a1->ncomp);

  for (unsigned i=0; i<a1->size; ++i) {
    for (unsigned j=0; j<a1->ncomp; ++j) {
      TEST_CHECK( gkyl_compare(a2_d[i*a2->ncomp+0*a1->ncomp+j], i*0.1+j, 1e-14) );
      TEST_CHECK( gkyl_compare(a2_d[i*a2->ncomp+1*a1->ncomp+j], i*0.1+j+0.5*(i*1.0+j), 1e-14) );
      TEST_CHECK( gkyl_compare(a2_d[i*a2->ncomp+2*a1->ncomp+j], i*0.1+j, 1e-14) );
    }
  }

  gkyl_array_release(a1);
  gkyl_array_release(a2);
}

void test_array_accumulate_offset_range()
{
  int shape[] = {10, 20};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 2, range.volume);
  struct gkyl_array *a2 = gkyl_array_new(GKYL_DOUBLE, 3*a1->ncomp, range.volume);

  // test a1 = a1+0.5*a2[a1->ncomp]
  gkyl_array_clear(a1, 0.5);
  gkyl_array_clear(a2, 1.5);

  gkyl_array_accumulate_offset_range(a1, 0.5, a2, 1*a1->ncomp, &range);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);
    double *a1d = gkyl_array_fetch(a1, loc);
    for (int i=0; i<a1->ncomp; ++i)
      TEST_CHECK( a1d[i] == 0.5+0.5*1.5 );
  }

  // test a2[a1->ncomp] = a2[a1->ncomp]+0.5*a1
  gkyl_array_clear(a1, 0.5);
  gkyl_array_clear(a2, 1.5);

  gkyl_array_accumulate_offset_range(a2, 0.5, a1, 1*a1->ncomp, &range);

  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);
    double *a2d = gkyl_array_fetch(a2, loc);
    for (int i=0; i<a1->ncomp; ++i) {
      TEST_CHECK( a2d[i+0*a1->ncomp] == 1.5 );
      TEST_CHECK( a2d[i+1*a1->ncomp] == 1.5+0.5*0.5 );
      TEST_CHECK( a2d[i+2*a1->ncomp] == 1.5 );
    }
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

  // test a2 = 0.5*a1
  gkyl_array_clear(a1, 0.5);
  gkyl_array_clear(a2, 1.5);

  gkyl_array_set_range(a2, 0.5, a1, &range);

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

void test_array_set_offset()
{
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 2, 10);
  struct gkyl_array *a2 = gkyl_array_new(GKYL_DOUBLE, 3*a1->ncomp, a1->size);

  double *a1_d = a1->data, *a2_d = a2->data;

  // Assign a component of the vector to the scalar.
  for (unsigned i=0; i<a1->size; ++i)
    for (unsigned j=0; j<a1->ncomp; ++j)
      a1_d[i*a1->ncomp+j] = i*0.2+2*j;

  for (unsigned i=0; i<a2->size; ++i)
    for (unsigned j=0; j<a2->ncomp/a1->ncomp; ++j)
      for (unsigned k=0; k<a1->ncomp; ++k)
        a2_d[i*a2->ncomp+j*a1->ncomp+k] = i*0.1+k;

  gkyl_array_set_offset(a1, 0.5, a2, 1*a1->ncomp);

  for (unsigned i=0; i<a1->size; ++i)
    for (unsigned j=0; j<a1->ncomp; ++j)
      TEST_CHECK( gkyl_compare(a1_d[i*a1->ncomp+j], 0.5*(i*0.1+j), 1e-14) );

  // Assign the scalar to a component of the vector.
  for (unsigned i=0; i<a1->size; ++i)
    for (unsigned j=0; j<a1->ncomp; ++j)
      a1_d[i*a1->ncomp+j] = i*0.2+2*j;

  gkyl_array_set_offset(a2, 2., a1, 1*a1->ncomp);

  for (unsigned i=0; i<a1->size; ++i) {
    for (unsigned j=0; j<a1->ncomp; ++j) {
      TEST_CHECK( gkyl_compare(a2_d[i*a2->ncomp+0*a1->ncomp+j], i*0.1+j, 1e-14) );
      TEST_CHECK( gkyl_compare(a2_d[i*a2->ncomp+1*a1->ncomp+j], 2.*(i*0.2+2*j), 1e-14) );
      TEST_CHECK( gkyl_compare(a2_d[i*a2->ncomp+2*a1->ncomp+j], i*0.1+j, 1e-14) );
    }
  }

  gkyl_array_release(a1);
  gkyl_array_release(a2);
}

void test_array_set_offset_range()
{
  int shape[] = {10, 20};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 2, range.volume);
  struct gkyl_array *a2 = gkyl_array_new(GKYL_DOUBLE, 3*a1->ncomp, range.volume);

  // test a1 = 0.1*a2[a1->ncomp]
  gkyl_array_clear_range(a1, 0.5, &range);
  gkyl_array_clear_range(a2, 1.5, &range);

  gkyl_array_set_offset_range(a1, 0.1, a2, 1*a1->ncomp, &range);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);
    double *a1d = gkyl_array_fetch(a1, loc);
    for (int i=0; i<a1->ncomp; ++i)
      TEST_CHECK( a1d[i] == 0.1*1.5 );
  }

  // test a2[a1->ncomp] = 0.1*a1
  gkyl_array_clear(a1, 0.5);
  gkyl_array_clear(a2, 1.5);

  gkyl_array_set_offset_range(a2, 0.1, a1, 1*a1->ncomp, &range);

  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);
    double *a2d = gkyl_array_fetch(a2, loc);
    for (int i=0; i<a1->ncomp; ++i) {
      TEST_CHECK( a2d[i+0*a1->ncomp] == 1.5 );
      TEST_CHECK( a2d[i+1*a1->ncomp] == 0.1*0.5 );
      TEST_CHECK( a2d[i+2*a1->ncomp] == 1.5 );
    }
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

void test_array_scale_by_cell()
{
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 3, 10);
  struct gkyl_array *s = gkyl_array_new(GKYL_DOUBLE, 1, 10);

  double *a1_d  = a1->data;
  for (unsigned i=0; i<a1->size; ++i) {
    a1_d[i] = i*1.0;
  }
  double *s_d  = s->data;
  for (unsigned i=0; i<s->size; ++i) {
    s_d[i] = i*1.0;
  }

  gkyl_array_scale_by_cell(a1, s);

  for (unsigned i=0; i<a1->size; ++i) {
    int fact = (i/a1->ncomp);
    TEST_CHECK( gkyl_compare(a1_d[i], i*1.0*fact, 1e-14) );
  }

  gkyl_array_release(a1);
  gkyl_array_release(s);
}

void test_array_shiftc()
{
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 3, 10);

  double s = -0.5;
  double *a1_d = a1->data;
  for (unsigned i=0; i<a1->size; ++i) {
    for (size_t k=0; k<a1->ncomp; ++k) a1_d[i*a1->ncomp+k] = i*2.0+k;
  }

  gkyl_array_shiftc(a1, s, 0);

  TEST_CHECK( gkyl_compare(a1_d[0], 0*1.0+0+s, 1e-14) );
  for (unsigned i=0; i<a1->size; ++i) {
    for (size_t k=1; k<a1->ncomp; ++k) 
      TEST_CHECK( gkyl_compare(a1_d[i*a1->ncomp+k], i*2.0+k, 1e-14) );
  }

  gkyl_array_release(a1);

  // Repeat the test but shifting another coefficient as well.
  int shiftks[] = {0, 2};
  int nks = sizeof(shiftks)/sizeof(shiftks[0]);

  struct gkyl_array *a2 = gkyl_array_new(GKYL_DOUBLE, 4, 8);
  double *a2_d = a2->data;
  for (unsigned i=0; i<a2->size; ++i) {
    for (size_t k=0; k<a2->ncomp; ++k) a2_d[i*a2->ncomp+k] = i*2.0+k;
  }

  for (size_t l=0; l<nks; l++)
    gkyl_array_shiftc(a2, s, shiftks[l]);

  for (unsigned i=0; i<a2->size; ++i) {
    for (size_t k=0; k<a2->ncomp; ++k) {
      bool isshifted = false;
      for (size_t l=0; l<nks; l++) {
        if (shiftks[l]==k) {isshifted = true; break;}
      }
      if (isshifted)
        TEST_CHECK( gkyl_compare(a2_d[i*a2->ncomp+k], i*2.0+k+s, 1e-14) );
      else
        TEST_CHECK( gkyl_compare(a2_d[i*a2->ncomp+k], i*2.0+k, 1e-14) );
    }
  }

  gkyl_array_release(a2);
}

void test_array_shiftc_range(bool on_gpu)
{
  int lower[] = {1}, upper[] = {10};
  struct gkyl_range range;
  gkyl_range_init(&range, 1, lower, upper);

  struct gkyl_array *a1_ho = gkyl_array_new(GKYL_DOUBLE, 3, range.volume);
  struct gkyl_array *a1 = a1_ho;
  if (on_gpu) a1 = gkyl_array_cu_dev_new(GKYL_DOUBLE, a1_ho->ncomp, a1_ho->size);

  double s = -0.5;
  double *a1_ho_d = a1_ho->data;
  for (unsigned i=0; i<a1_ho->size; ++i)
    for (size_t k=0; k<a1_ho->ncomp; ++k) a1_ho_d[i*a1_ho->ncomp+k] = i*2.0+k;

  if (on_gpu) gkyl_array_copy(a1, a1_ho);

  int lowerSub[] = {2}, upperSub[] = {6};
  struct gkyl_range subrange;
  gkyl_sub_range_init(&subrange, &range, lowerSub, upperSub);

  gkyl_array_shiftc_range(a1, s, 0, &subrange);

  if (on_gpu) gkyl_array_copy(a1_ho, a1);

  for (size_t k=1; k<a1_ho->ncomp; ++k) {
    for (unsigned i=0; i<1; ++i)
      TEST_CHECK( gkyl_compare(a1_ho_d[i*a1_ho->ncomp+k], i*2.0+k, 1e-14) );
    for (unsigned i=6; i<10; ++i)
      TEST_CHECK( gkyl_compare(a1_ho_d[i*a1_ho->ncomp+k], i*2.0+k, 1e-14) );
  }
  for (unsigned i=1; i<6; ++i)
    TEST_CHECK( gkyl_compare(a1_ho_d[i*a1_ho->ncomp+0], i*2.0+0+s, 1e-14) );

  gkyl_array_release(a1_ho);
  if (on_gpu) gkyl_array_release(a1);

  // Repeat the test but shifting another coefficient as well.
  int shiftks[] = {0, 2};
  int nks = sizeof(shiftks)/sizeof(shiftks[0]);

  lower[0] = 1;  upper[0] = 8;
  struct gkyl_range range2;
  gkyl_range_init(&range2, 1, lower, upper);

  struct gkyl_array *a2_ho = gkyl_array_new(GKYL_DOUBLE, 4, range2.volume);
  struct gkyl_array *a2 = a2_ho;
  if (on_gpu) a2 = gkyl_array_cu_dev_new(GKYL_DOUBLE, a2_ho->ncomp, a2_ho->size);

  double *a2_ho_d = a2_ho->data;
  for (unsigned i=0; i<a2_ho->size; ++i) {
    for (size_t k=0; k<a2_ho->ncomp; ++k) a2_ho_d[i*a2_ho->ncomp+k] = i*2.0+k;
  }

  if (on_gpu) gkyl_array_copy(a2, a2_ho);

  for (size_t l=0; l<nks; l++)
    gkyl_array_shiftc_range(a2, s, shiftks[l], &subrange);

  if (on_gpu) gkyl_array_copy(a2_ho, a2);

  for (size_t k=0; k<a2_ho->ncomp; ++k) {
    for (unsigned i=0; i<1; ++i)
      TEST_CHECK( gkyl_compare(a2_ho_d[i*a2_ho->ncomp+k], i*2.0+k, 1e-14) );
    for (unsigned i=6; i<8; ++i)
      TEST_CHECK( gkyl_compare(a2_ho_d[i*a2_ho->ncomp+k], i*2.0+k, 1e-14) );
  }
  for (unsigned i=1; i<6; ++i) {
    for (size_t k=0; k<a2_ho->ncomp; ++k) {
      bool isshifted = false;
      for (size_t l=0; l<nks; l++) {
        if (shiftks[l]==k) {isshifted = true; break;}
      }
      if (isshifted)
        TEST_CHECK( gkyl_compare(a2_ho_d[i*a2_ho->ncomp+k], i*2.0+k+s, 1e-14) );
      else
        TEST_CHECK( gkyl_compare(a2_ho_d[i*a2_ho->ncomp+k], i*2.0+k, 1e-14) );
    }
  }

  gkyl_array_release(a2);
  if (on_gpu) gkyl_array_release(a2);
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

void test_array_copy_buffer()
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

// function for use in the buffer_fn method
GKYL_CU_DH static void
buffer_fn(size_t nc, double *out, const double *inp, void *ctx)
{
  for (size_t i=0; i<nc; ++i)
    out[i] = 2*inp[i];
}
    
void test_array_copy_buffer_fn()
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
  gkyl_array_copy_to_buffer_fn(buff, arr, &sub_range,
    &(struct gkyl_array_copy_func) { .func = buffer_fn, .ctx = 0 }
  );

  long count = 0;
  gkyl_range_iter_init(&iter, &sub_range);
  while (gkyl_range_iter_next(&iter))
    TEST_CHECK( buff[count++] == 2*(iter.idx[0] + 10.5*iter.idx[1]) );

  gkyl_array_clear(arr, 0.0);
  // copy back from buffer
  gkyl_array_copy_from_buffer(arr, buff, &sub_range);

  gkyl_range_iter_init(&iter, &sub_range);
  while (gkyl_range_iter_next(&iter)) {
    double *d = gkyl_array_fetch(arr, gkyl_range_idx(&sub_range, iter.idx));
    TEST_CHECK( d[0]  == 2*(iter.idx[0] + 10.5*iter.idx[1]) );
  }

  gkyl_array_release(arr);
  gkyl_free(buff);
}

void test_array_flip_copy_buffer_fn()
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

  gkyl_array_flip_copy_to_buffer_fn(buff, arr, 0, &sub_range,
    &(struct gkyl_array_copy_func) { .func = buffer_fn, .ctx = 0 }
  );
  long count = 0;
  gkyl_range_iter_init(&iter, &sub_range);
  while (gkyl_range_iter_next(&iter))
    TEST_CHECK( buff[count++] == 2*((5+1)-iter.idx[0] + 10.5*iter.idx[1]) );

  gkyl_array_flip_copy_to_buffer_fn(buff, arr, 1, &sub_range,     
    &(struct gkyl_array_copy_func) { .func = buffer_fn, .ctx = 0 }
  );
  count = 0;
  gkyl_range_iter_init(&iter, &sub_range);
  while (gkyl_range_iter_next(&iter))
    TEST_CHECK( buff[count++] == 2*(iter.idx[0] + 10.5*((10+1)-iter.idx[1])) );

  gkyl_array_clear(arr, 0.0);
  // copy back from buffer
  gkyl_array_copy_from_buffer(arr, buff, &sub_range);

  gkyl_range_iter_init(&iter, &sub_range);
  while (gkyl_range_iter_next(&iter)) {
    double *d = gkyl_array_fetch(arr, gkyl_range_idx(&sub_range, iter.idx));
    //printf("%lg %lg\n", d[0], 2*(iter.idx[0] + 10.5*((10+1)-iter.idx[1])) );
    TEST_CHECK( d[0]  == 2*(iter.idx[0] + 10.5*((10+1)-iter.idx[1])) );
  }

  gkyl_array_release(arr);
  gkyl_free(buff);
}

void test_array_copy_range()
{
  int shape[] = {10, 20};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);

  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 10, range.volume);
  struct gkyl_array *a2 = gkyl_array_new(GKYL_DOUBLE, 10, range.volume);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    double *d = gkyl_array_fetch(a1, gkyl_range_idx(&range, iter.idx));
    d[0] = iter.idx[0] + 10.5*iter.idx[1];
  }

  // copy the contents of a1 into a2 over the specified range
  gkyl_array_copy_range(a2, a1, &range);
  
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    double *d = gkyl_array_fetch(a2, gkyl_range_idx(&range, iter.idx));
    TEST_CHECK( d[0]  == iter.idx[0] + 10.5*iter.idx[1] );
  }

  // clear array for second test
  gkyl_array_clear(a2, 0.0);
  // initialize left sub-range
  int lower_l[] = {range.lower[0], range.lower[1]}, upper_l[] = {range.lower[0], shape[1]/2-1};
  struct gkyl_range sub_range_l;
  gkyl_sub_range_init(&sub_range_l, &range, lower_l, upper_l);
  
  // initialize right sub-range
  int lower_r[] = {range.upper[0], shape[1]/2}, upper_r[] = {range.upper[0], range.upper[1]};
  struct gkyl_range sub_range_r;
  gkyl_sub_range_init(&sub_range_r, &range, lower_r, upper_r);

  gkyl_array_copy_range_to_range(a2, a1, &sub_range_r, &sub_range_l);

  gkyl_range_iter_init(&iter, &sub_range_r);
  while (gkyl_range_iter_next(&iter)) {
    int idx_l[GKYL_MAX_DIM];
    idx_l[0] = range.lower[0], idx_l[1] = iter.idx[1]-shape[1]/2;
    double *d = gkyl_array_fetch(a2, gkyl_range_idx(&sub_range_r, iter.idx));
    TEST_CHECK( d[0]  == idx_l[0] + 10.5*idx_l[1] );
    TEST_MSG("Expected: %.13e in cell (%d,%d)", iter.idx[0] + 10.5*iter.idx[1], iter.idx[0], iter.idx[1]);
    TEST_MSG("Produced: %.13e", d[0]);
  }

  gkyl_array_release(a1);
  gkyl_array_release(a2);
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
    struct gkyl_range sr = gkyl_range_split(&sub_range, num_split, tid);
    gkyl_array_copy_to_buffer(buff+loc, arr, &sr);
    
    loc += gkyl_range_split_len(&sr);
  }

  long count = 0;
  gkyl_range_iter_init(&iter, &sub_range);
  while (gkyl_range_iter_next(&iter))
    TEST_CHECK( buff[count++] == iter.idx[0] + 10.5*iter.idx[1] );

  gkyl_array_clear(arr, 0.0);

  loc = 0;
  for (int tid=0; tid<num_split; ++tid) {
    // copy back from buffer
    struct gkyl_range sr = gkyl_range_split(&sub_range, num_split, tid);
    gkyl_array_copy_from_buffer(arr, buff+loc, &sr);
    loc += gkyl_range_split_len(&sr);
  }

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
  int ncomp = 3, ncells = 200;
  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, ncomp, ncells);

  for (size_t i=0; i<arr->size; ++i) {
    double *d = gkyl_array_fetch(arr, i);
    for (size_t c=0; c<ncomp; ++c)
      d[c] = 0.5*i + 0.1*c;
  }
  
  double amin[ncomp], amax[ncomp];
  gkyl_array_reduce(&amin[0], arr, GKYL_MIN);
  gkyl_array_reduce(&amax[0], arr, GKYL_MAX);

  for (size_t c=0; c<ncomp; ++c) {
    TEST_CHECK( amin[c] == 0.1*c );
    TEST_CHECK( amax[c] == 0.5*(ncells-1) + 0.1*c );
  }

  gkyl_array_release(arr);
}

void test_reduce_dg()
{
  int poly_order = 1;
  int ncomp = 3;
  double lower[] = {-M_PI}, upper[] = {M_PI};
  int cells[] = {20};

  int ndim = sizeof(lower)/sizeof(lower[0]);

  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  int ghost[ndim];
  for (int d=0; d<ndim; d++) ghost[d] = 1;
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, ncomp*basis.num_basis, local_ext.volume);

  // Load 1D Gauss-Legendre nodes.
  int num_quad = poly_order+1;
  double ordinates1[num_quad];
  memcpy(ordinates1, gkyl_gauss_ordinates[num_quad], sizeof(double[num_quad]));

  // Create range to loop over nodes.
  int qshape[GKYL_MAX_DIM];
  for (int i=0; i<ndim; ++i) qshape[i] = num_quad;
  struct gkyl_range qrange;
  gkyl_range_init_from_shape(&qrange, ndim, qshape);

  // Create nodes array.
  int num_nodes = qrange.volume;
  struct gkyl_array *nodes = gkyl_array_new(GKYL_DOUBLE, ndim, num_nodes);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &qrange);
  while (gkyl_range_iter_next(&iter)) {
    long linc = gkyl_range_idx(&qrange, iter.idx);
    double *nod = gkyl_array_fetch(nodes, linc);
    for (int i=0; i<ndim; ++i)
      nod[i] = ordinates1[iter.idx[i]-qrange.lower[i]];
  }

  // Polulate arr with a function evaluated at nodes (transformed to modal).
  double arr_max[ncomp], arr_min[ncomp], arr_sum[ncomp];
  for (size_t i=0; i<ncomp; ++i) {
    arr_max[i] = -1e20;
    arr_min[i] = 1e20;
    arr_sum[i] = 0.0;
  }

  for (size_t i=0; i<arr->size; ++i) {

    double *arr_c = gkyl_array_fetch(arr, i);

    int idx[ndim];
    gkyl_range_inv_idx(&local_ext, i, idx);

    double xc[ndim];
    gkyl_rect_grid_cell_center(&grid, idx, xc);

    for (int ci=0; ci<ncomp; ci++) {
      double arr_nodal[num_nodes];
      for (size_t k=0; k<num_nodes; k++) {
        const double *nod = gkyl_array_cfetch(nodes, k);
        double x[ndim];
        for (int d=0; d<ndim; d++) x[d] = xc[d] + 0.5*grid.dx[0]*nod[d];
  
        arr_nodal[k] = (ci+1);
        for (int d=0; d<ndim; d++) arr_nodal[k] *= sin(((d+1)*2.0*M_PI/(upper[d]-lower[d]))*x[d]);
  
        arr_max[ci] = GKYL_MAX2(arr_max[ci], arr_nodal[k]);
        arr_min[ci] = GKYL_MIN2(arr_min[ci], arr_nodal[k]);
        arr_sum[ci] += arr_nodal[k];
      }
  
      for (size_t k=0; k<basis.num_basis; k++) {
        basis.quad_nodal_to_modal(arr_nodal, &arr_c[ci*basis.num_basis], k);
      }
    }
  }
  
  double amin[ncomp], amax[ncomp], asum[ncomp];
  for (int ci=0; ci<ncomp; ci++) {
    gkyl_array_reducec_dg(&amin[ci], arr, ci, GKYL_MIN, &basis);
    gkyl_array_reducec_dg(&amax[ci], arr, ci, GKYL_MAX, &basis);
    gkyl_array_reducec_dg(&asum[ci], arr, ci, GKYL_SUM, &basis);
  }

  for (int c=0; c<ncomp; ++c) {
    TEST_CHECK( gkyl_compare(amin[c], arr_min[c], 1e-14) );
    TEST_MSG( "%d MIN: Expected: %g | Got: %g\n", c, arr_min[c], amin[c] );
    TEST_CHECK( gkyl_compare(amax[c], arr_max[c], 1e-14) );
    TEST_MSG( "%d MAX: Expected: %g | Got: %g\n", c, arr_max[c], amax[c] );
    TEST_CHECK( gkyl_compare(asum[c], arr_sum[c], 1e-14) );
    TEST_MSG( "%d MAX: Expected: %g | Got: %g\n", c, arr_sum[c], asum[c] );
  }

  gkyl_array_release(arr);
  gkyl_array_release(nodes);
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

void test_sum_reduce_range()
{
  int shape[] = {10, 20};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, 3, range.volume);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);
    double *d = gkyl_array_fetch(arr, loc);
    d[0] = 0.5;
    d[1] = 1.5;
    d[2] = 2.5;
  }

  double asum[3];
  gkyl_array_reduce_range(asum, arr, GKYL_SUM, &range);

  TEST_CHECK( asum[0] == 0.5*range.volume );
  TEST_CHECK( asum[1] == 1.5*range.volume );
  TEST_CHECK( asum[2] == 2.5*range.volume );

  gkyl_array_release(arr);
}

void test_reduce_dg_range()
{
  int poly_order = 1;
  int ncomp = 3;
  double lower[] = {-M_PI}, upper[] = {M_PI};
  int cells[] = {20};

  int ndim = sizeof(lower)/sizeof(lower[0]);

  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  int ghost[ndim];
  for (int d=0; d<ndim; d++) ghost[d] = 1;
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, ncomp*basis.num_basis, local_ext.volume);

  // Load 1D Gauss-Legendre nodes.
  int num_quad = poly_order+1;
  double ordinates1[num_quad];
  memcpy(ordinates1, gkyl_gauss_ordinates[num_quad], sizeof(double[num_quad]));

  // Create range to loop over nodes.
  int qshape[GKYL_MAX_DIM];
  for (int i=0; i<ndim; ++i) qshape[i] = num_quad;
  struct gkyl_range qrange;
  gkyl_range_init_from_shape(&qrange, ndim, qshape);

  // Create nodes array.
  int num_nodes = qrange.volume;
  struct gkyl_array *nodes = gkyl_array_new(GKYL_DOUBLE, ndim, num_nodes);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &qrange);
  while (gkyl_range_iter_next(&iter)) {
    long linc = gkyl_range_idx(&qrange, iter.idx);
    double *nod = gkyl_array_fetch(nodes, linc);
    for (int i=0; i<ndim; ++i)
      nod[i] = ordinates1[iter.idx[i]-qrange.lower[i]];
  }

  // Polulate arr with a function evaluated at nodes (transformed to modal).
  double arr_max[ncomp], arr_min[ncomp], arr_sum[ncomp];
  for (size_t i=0; i<ncomp; ++i) {
    arr_max[i] = -1e20;
    arr_min[i] = 1e20;
    arr_sum[i] = 0.0;
  }

  for (size_t i=0; i<arr->size; ++i) {

    double *arr_c = gkyl_array_fetch(arr, i);

    int idx[ndim];
    gkyl_range_inv_idx(&local_ext, i, idx);

    double xc[ndim];
    gkyl_rect_grid_cell_center(&grid, idx, xc);

    for (int ci=0; ci<ncomp; ci++) {
      double arr_nodal[num_nodes];
      for (size_t k=0; k<num_nodes; k++) {
        const double *nod = gkyl_array_cfetch(nodes, k);
        double x[ndim];
        for (int d=0; d<ndim; d++) x[d] = xc[d] + 0.5*grid.dx[0]*nod[d];
  
        arr_nodal[k] = (ci+1);
        for (int d=0; d<ndim; d++) arr_nodal[k] *= sin(((d+1)*2.0*M_PI/(upper[d]-lower[d]))*x[d]);
  
        arr_max[ci] = GKYL_MAX2(arr_max[ci], arr_nodal[k]);
        arr_min[ci] = GKYL_MIN2(arr_min[ci], arr_nodal[k]);
        arr_sum[ci] += arr_nodal[k];
      }
  
      for (size_t k=0; k<basis.num_basis; k++) {
        basis.quad_nodal_to_modal(arr_nodal, &arr_c[ci*basis.num_basis], k);
      }
    }
  }
  
  double amin[ncomp], amax[ncomp], asum[ncomp];
  for (int ci=0; ci<ncomp; ci++) {
    gkyl_array_reducec_dg_range(&amin[ci], arr, ci, GKYL_MIN, &basis, &local);
    gkyl_array_reducec_dg_range(&amax[ci], arr, ci, GKYL_MAX, &basis, &local);
    gkyl_array_reducec_dg_range(&asum[ci], arr, ci, GKYL_SUM, &basis, &local);
  }

  for (int c=0; c<ncomp; ++c) {
    TEST_CHECK( gkyl_compare(amin[c], arr_min[c], 1e-14) );
    TEST_MSG( "%d MIN: Expected: %g | Got: %g\n", c, arr_min[c], amin[c] );
    TEST_CHECK( gkyl_compare(amax[c], arr_max[c], 1e-14) );
    TEST_MSG( "%d MAX: Expected: %g | Got: %g\n", c, arr_max[c], amax[c] );
    TEST_CHECK( gkyl_compare(asum[c], arr_sum[c], 1e-14) );
    TEST_MSG( "%d MAX: Expected: %g | Got: %g\n", c, arr_sum[c], asum[c] );
  }

  gkyl_array_release(arr);
  gkyl_array_release(nodes);
}


void
test_grid_sub_array_read_1()
{
  double lower[] = {1.0, 1.0}, upper[] = {2.5, 5.0};
  int cells[] = {20, 60};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);

  int nghost[] = { 1, 2 };
  struct gkyl_range range, ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &ext_range, &range);

  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, 2, ext_range.volume);
  gkyl_array_clear(arr, 0.0);

  // set some values in array
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);

    double *d = gkyl_array_fetch(arr, loc);
    for (int k=0; k<2; ++k)
      d[k] = (10.5*iter.idx[0] + 220.5*iter.idx[1])*(k+0.5);
  }

  gkyl_grid_sub_array_write(&grid, &range, 0, arr, "ctest_grid_sub_array_1.gkyl");

  struct gkyl_array *arr2 = gkyl_array_new(GKYL_DOUBLE, 2, ext_range.volume);
  gkyl_array_clear(arr2, 0.0);

  struct gkyl_rect_grid grid2;

  // read just header
  FILE *fp;
  with_file(fp, "ctest_grid_sub_array_1.gkyl", "r") {
    struct gkyl_array_header_info hdr;
    
    int status = gkyl_grid_sub_array_header_read_fp(&grid2, &hdr, fp);
    TEST_CHECK( status == 0 );

    TEST_CHECK( hdr.file_type == 1);
    TEST_CHECK( hdr.etype == GKYL_DOUBLE);

    long tot_cells = 1L;
    TEST_CHECK( grid.ndim == grid2.ndim );
    for (int d=0; d<grid.ndim; ++d) {
      TEST_CHECK( grid.lower[d] == grid2.lower[d] );
      TEST_CHECK( grid.upper[d] == grid2.upper[d] );
      TEST_CHECK( grid.cells[d] == grid2.cells[d] );
      TEST_CHECK( grid.dx[d] == grid2.dx[d] );

      tot_cells *= grid.cells[d];
    }
    TEST_CHECK( grid.cellVolume == grid2.cellVolume );

    TEST_CHECK( hdr.esznc = arr->esznc );
    TEST_CHECK( hdr.tot_cells == tot_cells );

    TEST_CHECK( 0 == hdr.meta_size );
    TEST_CHECK( 1 == hdr.nrange );
  }

  int file_type = gkyl_get_gkyl_file_type("ctest_grid_sub_array_1.gkyl");
  TEST_CHECK( 1 == file_type );
  
  // read back the grid and the array
  int err =
    gkyl_grid_sub_array_read(&grid2, &range, arr2, "ctest_grid_sub_array_1.gkyl");

  TEST_CHECK( err < 1 );
  
  if (err < 1) {

    TEST_CHECK( grid.ndim == grid2.ndim );
    for (int d=0; d<grid.ndim; ++d) {
      TEST_CHECK( grid.lower[d] == grid2.lower[d] );
      TEST_CHECK( grid.upper[d] == grid2.upper[d] );
      TEST_CHECK( grid.cells[d] == grid2.cells[d] );
      TEST_CHECK( grid.dx[d] == grid2.dx[d] );
    }
    TEST_CHECK( grid.cellVolume == grid2.cellVolume );
    
    gkyl_range_iter_init(&iter, &range);
    while (gkyl_range_iter_next(&iter)) {
      long loc = gkyl_range_idx(&range, iter.idx);
      
      const double *rhs = gkyl_array_cfetch(arr, loc);
      const double *lhs = gkyl_array_cfetch(arr2, loc);
      for (int k=0; k<2; ++k)
        TEST_CHECK( lhs[k] == rhs[k] );
    }
  }
  
  gkyl_array_release(arr);
  gkyl_array_release(arr2);
}

void
test_grid_sub_array_read_2()
{
  double lower[] = {1.0, 1.0}, upper[] = {2.5, 5.0};
  int cells[] = {20, 60};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);

  int nghost[] = { 1, 2 };
  struct gkyl_range range, ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &ext_range, &range);

  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, 2, ext_range.volume);
  gkyl_array_clear(arr, 0.0);

  // set some values in array
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);

    double *d = gkyl_array_fetch(arr, loc);
    for (int k=0; k<2; ++k)
      d[k] = (10.5*iter.idx[0] + 220.5*iter.idx[1])*(k+0.5);
  }

  gkyl_grid_sub_array_write(&grid, &range, 0, arr, "ctest_grid_sub_array_2.gkyl");

  struct gkyl_rect_grid grid2;

  struct gkyl_range srange;
  gkyl_range_init(&srange, grid.ndim, (int[]) { 5, 5 }, (int[]) { 10, 15 });

  struct gkyl_array *arr2 = gkyl_array_new(GKYL_DOUBLE, 2, srange.volume);
  gkyl_array_clear(arr2, 0.0);

  // read back the grid and the array
  int err =
    gkyl_grid_sub_array_read(&grid2, &srange, arr2, "ctest_grid_sub_array_2.gkyl");

  TEST_CHECK( err < 1 );
  
  if (err < 1) {

    TEST_CHECK( grid.ndim == grid2.ndim );
    for (int d=0; d<grid.ndim; ++d) {
      TEST_CHECK( grid.lower[d] == grid2.lower[d] );
      TEST_CHECK( grid.upper[d] == grid2.upper[d] );
      TEST_CHECK( grid.cells[d] == grid2.cells[d] );
      TEST_CHECK( grid.dx[d] == grid2.dx[d] );
    }
    TEST_CHECK( grid.cellVolume == grid2.cellVolume );
    
    gkyl_range_iter_init(&iter, &srange);
    while (gkyl_range_iter_next(&iter)) {
      const double *rhs = gkyl_array_cfetch(arr, gkyl_range_idx(&range, iter.idx));
      const double *lhs = gkyl_array_cfetch(arr2, gkyl_range_idx(&srange, iter.idx));
      for (int k=0; k<2; ++k)
        TEST_CHECK( lhs[k] == rhs[k] );
    }
  }
  
  gkyl_array_release(arr);
  gkyl_array_release(arr2);
}

void
test_grid_array_new_from_file_1()
{
  double lower[] = {1.0, 1.0}, upper[] = {2.5, 5.0};
  int cells[] = {20, 60};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);

  int nghost[] = { 1, 2 };
  struct gkyl_range range, ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &ext_range, &range);

  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, 2, ext_range.volume);
  gkyl_array_clear(arr, 0.0);

  // set some values in array
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);

    double *d = gkyl_array_fetch(arr, loc);
    for (int k=0; k<2; ++k)
      d[k] = (10.5*iter.idx[0] + 220.5*iter.idx[1])*(k+0.5);
  }

  gkyl_grid_sub_array_write(&grid, &range, 0, arr, "ctest_grid_array_new_from_file_1.gkyl");

  struct gkyl_rect_grid grid2;
  struct gkyl_array *arr2 = gkyl_grid_array_new_from_file(&grid2, "ctest_grid_array_new_from_file_1.gkyl");

  TEST_CHECK( arr2->type == GKYL_DOUBLE );

  TEST_CHECK( grid.ndim == grid2.ndim );
  for (int d=0; d<grid.ndim; ++d) {
    TEST_CHECK( grid.lower[d] == grid2.lower[d] );
    TEST_CHECK( grid.upper[d] == grid2.upper[d] );
    TEST_CHECK( grid.cells[d] == grid2.cells[d] );
    TEST_CHECK( grid.dx[d] == grid2.dx[d] );
  }
  TEST_CHECK( grid.cellVolume == grid2.cellVolume );  

  // NOTE: array read from file is not the same shape as "arr". This
  // is because the new_from_file does not read ghost cells.
  long lhs_loc = 0;
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);

    const double *rhs = gkyl_array_cfetch(arr, loc);
    const double *lhs = gkyl_array_cfetch(arr2, lhs_loc++);
    for (int k=0; k<2; ++k)
      TEST_CHECK( lhs[k] == rhs[k] );
  }
  
  gkyl_array_release(arr);
  gkyl_array_release(arr2);
}

void
test_grid_array_read_p1(void)
{
  // read just header
  struct gkyl_rect_grid grid;  
  struct gkyl_array_header_info hdr;
  FILE *fp = 0;  
  with_file(fp, "data/unit/ser-euler_riem_2d_hllc-euler_1.gkyl", "r") {
    
    int status = gkyl_grid_sub_array_header_read_fp(&grid, &hdr, fp);
    TEST_CHECK( status == 0 );

    TEST_CHECK( hdr.file_type == gkyl_file_type_int[GKYL_FIELD_DATA_FILE]);
    TEST_CHECK( hdr.etype == GKYL_DOUBLE);

    TEST_CHECK( hdr.esznc = 5*sizeof(double));
    TEST_CHECK( 50*50 == hdr.tot_cells );

    TEST_CHECK( 50 == grid.cells[0] );
    TEST_CHECK( 50 == grid.cells[1] );
    
    TEST_CHECK( 0.0 == grid.lower[0] );
    TEST_CHECK( 0.0 == grid.lower[1] );

    TEST_CHECK( 1.0 == grid.upper[0] );
    TEST_CHECK( 1.0 == grid.upper[1] );

    TEST_CHECK( 1 == hdr.nrange );

    TEST_CHECK( hdr.meta_size > 0 );

    mpack_tree_t tree;
    mpack_tree_init_data(&tree, hdr.meta, hdr.meta_size);
    mpack_tree_parse(&tree);

    mpack_node_t root = mpack_tree_root(&tree);
    TEST_CHECK(mpack_node_type(root) == mpack_type_map);

    mpack_node_t tm_node = mpack_node_map_cstr(root, "time");
    TEST_CHECK( mpack_node_double(tm_node) >  0.80675 );

    mpack_node_t fr_node = mpack_node_map_cstr(root, "frame");
    TEST_CHECK( mpack_node_i64(fr_node) == 1 );

    status = mpack_tree_destroy(&tree);
    TEST_CHECK( mpack_ok == status );

    gkyl_free(hdr.meta);
  }

  size_t nc = hdr.esznc/gkyl_elem_type_size[hdr.etype];

  int nghost[] = { 1, 2 };
  struct gkyl_range range, ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &ext_range, &range);

  // read serial data for comparison
  struct gkyl_rect_grid s_grid;
  struct gkyl_array *s_arr = gkyl_array_new(hdr.etype, nc, ext_range.volume);
  int s_status = gkyl_grid_sub_array_read(&s_grid, &range, s_arr,
    "data/unit/ser-euler_riem_2d_hllc-euler_1.gkyl");

  // read parallel data (whole domain)
  do {
    struct gkyl_rect_grid p_grid;  
    struct gkyl_array *p_arr = gkyl_array_new(hdr.etype,  nc, ext_range.volume);
    int p_status = gkyl_grid_sub_array_read(&p_grid, &range, p_arr,
      "data/unit/euler_riem_2d_hllc-euler_1.gkyl");

    TEST_CHECK( 0 == p_status );
    
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &range);
    while (gkyl_range_iter_next(&iter)) {
      long loc = gkyl_range_idx(&range, iter.idx);
      const double *s_dat = gkyl_array_fetch(s_arr, loc);
      const double *p_dat = gkyl_array_fetch(p_arr, loc);
      
      for (int c=0; c<nc; ++c)
        TEST_CHECK( gkyl_compare_double(s_dat[c], p_dat[c], 1e-15) );
    }
    gkyl_array_release(p_arr);
  } while (0);

  // read parallel data (partial domain)
  do {
    struct gkyl_range prange;
    gkyl_range_init(&prange, 2, (int[]) { 10, 10 }, (int[]) { 30, 40 });
    
    struct gkyl_rect_grid p_grid;
    struct gkyl_array *p_arr = gkyl_array_new(hdr.etype,  nc, prange.volume);
    int p_status = gkyl_grid_sub_array_read(&p_grid, &prange, p_arr,
      "data/unit/euler_riem_2d_hllc-euler_1.gkyl");

    TEST_CHECK( 0 == p_status );
    
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &prange);
    while (gkyl_range_iter_next(&iter)) {

      const double *s_dat = gkyl_array_fetch(s_arr, gkyl_range_idx(&range, iter.idx));
      const double *p_dat = gkyl_array_fetch(p_arr, gkyl_range_idx(&prange, iter.idx));
      
      for (int c=0; c<nc; ++c)
        TEST_CHECK( gkyl_compare_double(s_dat[c], p_dat[c], 1e-15) );
    }
    gkyl_array_release(p_arr);
  } while (0);

  // read parallel data (partial domain)
  do {
    struct gkyl_range prange;
    gkyl_range_init(&prange, 2, (int[]) { 4, 5 }, (int[]) { 10, 10 });
    
    struct gkyl_rect_grid p_grid;
    struct gkyl_array *p_arr = gkyl_array_new(hdr.etype,  nc, prange.volume);
    int p_status = gkyl_grid_sub_array_read(&p_grid, &prange, p_arr,
      "data/unit/euler_riem_2d_hllc-euler_1.gkyl");

    TEST_CHECK( 0 == p_status );
    
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &prange);
    while (gkyl_range_iter_next(&iter)) {

      const double *s_dat = gkyl_array_fetch(s_arr, gkyl_range_idx(&range, iter.idx));
      const double *p_dat = gkyl_array_fetch(p_arr, gkyl_range_idx(&prange, iter.idx));
      
      for (int c=0; c<nc; ++c)
        TEST_CHECK( gkyl_compare_double(s_dat[c], p_dat[c], 1e-15) );
    }
    gkyl_array_release(p_arr);
  } while (0);

  gkyl_array_release(s_arr);
}

static void
test_array_from_buff(void)
{
  double *buff = gkyl_malloc(sizeof(double[400]));
  
  struct gkyl_array *arr = gkyl_array_new_from_buff(GKYL_DOUBLE, 1, 200, buff);

  TEST_CHECK( gkyl_array_is_using_buffer(arr) == true );

  TEST_CHECK( arr->type = GKYL_DOUBLE );
  TEST_CHECK( arr->elemsz == sizeof(double) );
  TEST_CHECK( arr->ncomp == 1 );
  TEST_CHECK( arr->size == 20*10 );
  TEST_CHECK( arr->ref_count.count == 1 );

  TEST_CHECK( arr->on_dev == arr );

  TEST_CHECK( gkyl_array_is_cu_dev(arr) == false );

  gkyl_array_clear(arr, 0.0);
  
  double *arrData  = arr->data;
  for (unsigned i=0; i<arr->size; ++i){
    TEST_CHECK( arrData[i] == 0. );
    arrData[i] = (i+0.5)*0.1;
  }

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

  // acquire pointer
  struct gkyl_array *crr = gkyl_array_acquire(arr);

  TEST_CHECK( crr->ref_count.count == 2 );
  TEST_CHECK( arr->ref_count.count == 2 );

  struct gkyl_array *drr = gkyl_array_acquire(crr);

  TEST_CHECK( drr->ref_count.count == 3 );
  TEST_CHECK( crr->ref_count.count == 3 );  
  TEST_CHECK( arr->ref_count.count == 3 );
  
  gkyl_array_release(crr);
  TEST_CHECK( arr->ref_count.count == 2 );
  gkyl_array_release(drr);
  TEST_CHECK( arr->ref_count.count == 1 );

  gkyl_free(buff);
  
  gkyl_array_release(arr);
  gkyl_array_release(brr);
}

// Cuda specific tests
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

  gkyl_array_copy(arr, arr_cu);

  double *arrData = arr->data;
  for (unsigned i=0; i<arr->size; ++i) {
    TEST_CHECK( arrData[i] == 0. );
    arrData[i] = (i+0.5)*0.1;  
  }

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
  int nfail = cu_array_test_and_flip_sign(arr_cu->on_dev);
  TEST_CHECK( nfail == 0 );

  // restore arr_cu by copying from arr_cu_cl, and test again
  gkyl_array_copy(arr_cu, arr_cu_cl);
  nfail = cu_array_test_and_flip_sign(arr_cu->on_dev);
  TEST_CHECK( nfail == 0 );

  // check arr_cu_cl on device and flip sign
  nfail = cu_array_test_and_flip_sign(arr_cu_cl->on_dev);
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
  gkyl_array_clear(a1_cu, 0.5);
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
  gkyl_array_clear_range(a1_cu, 0.5, &range);

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

  gkyl_array_accumulate(a1_cu, 0.5, a2_cu);

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
  int shape[] = {20, 10};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 8, range.volume);
  struct gkyl_array *a2 = gkyl_array_new(GKYL_DOUBLE, 3, range.volume);

  // make device copies of arrays
  struct gkyl_array *a1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 8, range.volume);
  struct gkyl_array *a2_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3, range.volume);

  // initialize data
  double *a1_d  = a1->data, *a2_d = a2->data;
  for (unsigned i=0; i<a1->size; ++i) {
    for(unsigned c=0; c<a1->ncomp; ++c) 
      a1_d[c+a1->ncomp*i] = i*1.0 + .01*c;
    for(unsigned c=0; c<a2->ncomp; ++c) 
      a2_d[c+a2->ncomp*i] = i*0.1 + .01*c;
  }

  // copy initialized arrays to device
  gkyl_array_copy(a1_cu, a1);
  gkyl_array_copy(a2_cu, a2);

  // a1 = a1 + 0.5*a2
  gkyl_array_accumulate_range(a1_cu, 0.5, a2_cu, &range);

 // copy from device and check if things are ok
  gkyl_array_copy(a1, a1_cu);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);
    double *a1d = gkyl_array_fetch(a1, loc);
    for (int c=0; c<a2->ncomp; ++c) {
      TEST_CHECK( gkyl_compare( a1d[c], (.01+.5*.01)*c + loc*1.0 + .5*loc*.1, 1e-14 ) );
    }
    for (int c=a2->ncomp; c<a1->ncomp; ++c)
      TEST_CHECK( a1d[c] == (.01)*c + loc*1.0 );
  }

  // test a2 = a2 + 0.5*a
  gkyl_array_clear(a1_cu, 0.5);
  gkyl_array_clear(a2_cu, 1.5);

  gkyl_array_accumulate_range(a2_cu, 0.5, a1_cu, &range);

 // copy from device and check if things are ok
  gkyl_array_copy(a2, a2_cu);
  gkyl_range_iter_init(&iter, &range);

  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);
    double *a2d = gkyl_array_fetch(a2, loc);
    for (int i=0; i<a2->ncomp; ++i)
      TEST_CHECK( a2d[i] == 1.5 + 0.5*0.5 );
  }  

  gkyl_array_release(a1);
  gkyl_array_release(a2);
  gkyl_array_release(a1_cu);
  gkyl_array_release(a2_cu);
}

void test_cu_array_accumulate_offset()
{
  // create host arrays 
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 2, 10);
  struct gkyl_array *a2 = gkyl_array_new(GKYL_DOUBLE, 3*a1->ncomp, a1->size);
  // make device copies
  struct gkyl_array *a1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, a1->ncomp, a1->size);
  struct gkyl_array *a2_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, a2->ncomp, a2->size);

  // initialize data
  double *a1_d  = a1->data, *a2_d = a2->data;
  for (unsigned i=0; i<a1->size; ++i)
    for (unsigned j=0; j<a1->ncomp; ++j)
      a1_d[i*a1->ncomp+j] = i*1.0+j;

  for (unsigned i=0; i<a2->size; ++i)
    for (unsigned j=0; j<a2->ncomp/a1->ncomp; ++j)
      for (unsigned k=0; k<a1->ncomp; ++k)
        a2_d[i*a2->ncomp+j*a1->ncomp+k] = i*0.1+k;

  // copy initialized arrays to device
  gkyl_array_copy(a1_cu, a1);
  gkyl_array_copy(a2_cu, a2);

  // test a1 = 0.1*a2[a1->ncomp]
  gkyl_array_accumulate_offset(a1_cu, 0.5, a2_cu, 1*a1->ncomp);

  gkyl_array_copy(a1, a1_cu);
  for (unsigned i=0; i<a1->size; ++i)
    for (unsigned j=0; j<a1->ncomp; ++j)
      TEST_CHECK( gkyl_compare(a1_d[i*a1->ncomp+j], i*1.0+j+0.5*(i*0.1+j), 1e-14) );

  // test a2[a1->ncomp] = 0.1*a1
  for (unsigned i=0; i<a1->size; ++i)
    for (unsigned j=0; j<a1->ncomp; ++j)
      a1_d[i*a1->ncomp+j] = i*1.0+j;
  gkyl_array_copy(a1_cu, a1);

  gkyl_array_accumulate_offset(a2_cu, 0.5, a1_cu, 1*a1->ncomp);

  gkyl_array_copy(a2, a2_cu);
  for (unsigned i=0; i<a1->size; ++i) {
    for (unsigned j=0; j<a1->ncomp; ++j) {
      TEST_CHECK( gkyl_compare(a2_d[i*a2->ncomp+0*a1->ncomp+j], i*0.1+j, 1e-14) );
      TEST_CHECK( gkyl_compare(a2_d[i*a2->ncomp+1*a1->ncomp+j], i*0.1+j+0.5*(i*1.0+j), 1e-14) );
      TEST_CHECK( gkyl_compare(a2_d[i*a2->ncomp+2*a1->ncomp+j], i*0.1+j, 1e-14) );
    }
  }

  gkyl_array_release(a1);
  gkyl_array_release(a2);
  gkyl_array_release(a1_cu);
  gkyl_array_release(a2_cu);
}

void test_cu_array_accumulate_offset_range()
{
  int shape[] = {20, 10};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 2, range.volume);
  struct gkyl_array *a2 = gkyl_array_new(GKYL_DOUBLE, 3*a1->ncomp, range.volume);

  // make device copies of arrays
  struct gkyl_array *a1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, a1->ncomp, range.volume);
  struct gkyl_array *a2_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, a2->ncomp, range.volume);

  // initialize data
  double *a1_d  = a1->data, *a2_d = a2->data;
  for (unsigned i=0; i<a1->size; ++i) {
    for(unsigned c=0; c<a1->ncomp; ++c) 
      a1_d[i*a1->ncomp+c] = i*1.0 + .01*c;
    for(unsigned j=0; j<a2->ncomp/a1->ncomp; ++j)
      for(unsigned c=0; c<a1->ncomp; ++c) 
        a2_d[i*a2->ncomp+j*a1->ncomp+c] = i*0.1 + .01*c;
  }

  // copy initialized arrays to device
  gkyl_array_copy(a1_cu, a1);
  gkyl_array_copy(a2_cu, a2);

  // a1 = a1 + 0.5*a2[1*a1->ncomp]
  gkyl_array_accumulate_offset_range(a1_cu, 0.5, a2_cu, 1*a1->ncomp, &range);

  gkyl_array_copy(a1, a1_cu);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);
    double *a1d = gkyl_array_fetch(a1, loc);
    for (int c=0; c<a1->ncomp; ++c)
      TEST_CHECK( gkyl_compare( a1d[c], loc*1.0+.01*c + 0.5*(loc*0.1+0.01*c), 1e-14 ) );
  }

  // test a2[1*a1->ncomp] = a2[1*a1->ncomp] + 0.5*a1
  gkyl_array_clear(a1_cu, 0.5);
  gkyl_array_clear(a2_cu, 1.5);

  gkyl_array_accumulate_offset_range(a2_cu, 0.5, a1_cu, 1*a1->ncomp, &range);

  gkyl_array_copy(a2, a2_cu);
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);
    double *a2d = gkyl_array_fetch(a2, loc);
    for (int i=0; i<a1->ncomp; ++i) {
      TEST_CHECK( a2d[i+0*a1->ncomp] == 1.5 );
      TEST_CHECK( a2d[i+1*a1->ncomp] == 1.5+0.5*0.5 );
      TEST_CHECK( a2d[i+2*a1->ncomp] == 1.5 );
    }
  }  

  gkyl_array_release(a1);
  gkyl_array_release(a2);
  gkyl_array_release(a1_cu);
  gkyl_array_release(a2_cu);
}

void test_cu_array_accumulate_range_4d()
{
  int lower[] = { 1, 1, 1, 1 };
  int upper[] = { 46, 46, 32, 32};
  struct gkyl_range range;
  gkyl_range_init(&range, 4, lower, upper);
  
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 8, range.volume);
  struct gkyl_array *a2 = gkyl_array_new(GKYL_DOUBLE, 3, range.volume);

  // make device copies of arrays
  struct gkyl_array *a1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 8, range.volume);
  struct gkyl_array *a2_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3, range.volume);

  // initialize data
  gkyl_array_clear(a1, 0.5);
  gkyl_array_clear(a2, 1.5);
  
  // copy initialized arrays to device
  gkyl_array_copy(a1_cu, a1);
  gkyl_array_copy(a2_cu, a2);

  int slower[] = { 2, 2, 1, 1 };
  int supper[] = { 45, 45, 32, 32 };
  struct gkyl_range sub_range;
  gkyl_sub_range_init(&sub_range, &range, slower, supper);

  // a1 = a1 + 0.5*a2 (only first 3 components of a1 are modified)
  gkyl_array_accumulate_range(a1_cu, 0.5, a2_cu, &sub_range);

 // copy from device and check if things are ok
  gkyl_array_clear(a1, 0.0); 
  gkyl_array_copy(a1, a1_cu);
  
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &sub_range);
  
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);
    double *a1d = gkyl_array_fetch(a1, loc);
    for (int c=0; c<3; ++c)
      TEST_CHECK( a1d[c] == 0.5 + 0.5*1.5 );
    for (int c=3; c<8; ++c)
      TEST_CHECK( a1d[c] == 0.5 );
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
  gkyl_array_set(b_cu, 0.5, a1_cu);
  gkyl_array_accumulate(b_cu, 2.5, a2_cu);

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

  gkyl_array_set(a1_cu, 0.5, a2_cu);

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
  gkyl_array_clear(a1_cu, 0.5);
  gkyl_array_clear(a2_cu, 1.5);

  gkyl_array_set_range(a1_cu, 0.5, a2_cu, &range);

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
  gkyl_array_clear(a1_cu, 0.5);
  gkyl_array_clear(a2_cu, 1.5);

  gkyl_array_set_range(a2_cu, 0.5, a1_cu, &range);

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

void test_cu_array_set_offset()
{
  // create host arrays 
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 2, 10);
  struct gkyl_array *a2 = gkyl_array_new(GKYL_DOUBLE, 3*a1->ncomp, a1->size);
  // make device copies
  struct gkyl_array *a1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, a1->ncomp, a1->size);
  struct gkyl_array *a2_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, a2->ncomp, a2->size);

  // initialize data
  double *a1_d  = a1->data, *a2_d = a2->data;
  for (unsigned i=0; i<a1->size; ++i)
    for (unsigned j=0; j<a1->ncomp; ++j)
      a1_d[i*a1->ncomp+j] = i*1.0+j;

  for (unsigned i=0; i<a2->size; ++i)
    for (unsigned j=0; j<a2->ncomp/a1->ncomp; ++j)
      for (unsigned k=0; k<a1->ncomp; ++k)
        a2_d[i*a2->ncomp+j*a1->ncomp+k] = i*0.1+k;

  // copy initialized arrays to device
  gkyl_array_copy(a1_cu, a1);
  gkyl_array_copy(a2_cu, a2);

  // test a1 = 0.5*a2[a1->ncomp]
  gkyl_array_set_offset(a1_cu, 0.5, a2_cu, 1*a1->ncomp);

  gkyl_array_copy(a1, a1_cu);
  for (unsigned i=0; i<a1->size; ++i)
    for (unsigned j=0; j<a1->ncomp; ++j)
      TEST_CHECK( gkyl_compare(a1_d[i*a1->ncomp+j], 0.5*(i*0.1+j), 1e-14) );

  // test a2[a1->ncomp] = 0.1*a1
  for (unsigned i=0; i<a1->size; ++i)
    for (unsigned j=0; j<a1->ncomp; ++j)
      a1_d[i*a1->ncomp+j] = i*1.0+j;
  gkyl_array_copy(a1_cu, a1);

  gkyl_array_set_offset(a2_cu, 0.5, a1_cu, 1*a1->ncomp);

  gkyl_array_copy(a2, a2_cu);
  for (unsigned i=0; i<a1->size; ++i) {
    for (unsigned j=0; j<a1->ncomp; ++j) {
      TEST_CHECK( gkyl_compare(a2_d[i*a2->ncomp+0*a1->ncomp+j], i*0.1+j, 1e-14) );
      TEST_CHECK( gkyl_compare(a2_d[i*a2->ncomp+1*a1->ncomp+j], 0.5*(i*1.0+j), 1e-14) );
      TEST_CHECK( gkyl_compare(a2_d[i*a2->ncomp+2*a1->ncomp+j], i*0.1+j, 1e-14) );
    }
  }

  gkyl_array_release(a1);
  gkyl_array_release(a2);
  gkyl_array_release(a1_cu);
  gkyl_array_release(a2_cu);
}

void test_cu_array_set_offset_range()
{
  int shape[] = {10, 20};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 2, range.volume);
  struct gkyl_array *a2 = gkyl_array_new(GKYL_DOUBLE, 3*a1->ncomp, range.volume);

  // make device copies of arrays
  struct gkyl_array *a1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, a1->ncomp, range.volume);
  struct gkyl_array *a2_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, a2->ncomp, range.volume);

  // copy host arrays to device
  gkyl_array_copy(a1_cu, a1);
  gkyl_array_copy(a2_cu, a2);

  // test a1 = 0.5*a2[a1->ncomp]
  gkyl_array_clear(a1_cu, 0.5);
  gkyl_array_clear(a2_cu, 1.5);

  gkyl_array_set_offset_range(a1_cu, 0.5, a2_cu, 1*a1->ncomp, &range);

  gkyl_array_copy(a1, a1_cu);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);
    double *a1d = gkyl_array_fetch(a1, loc);
    for (int i=0; i<a1->ncomp; ++i)
      TEST_CHECK( a1d[i] == 0.5*1.5 );
  }

  // test a2[a1->ncomp] = 0.5*a1
  gkyl_array_clear(a1_cu, 0.5);
  gkyl_array_clear(a2_cu, 1.5);

  gkyl_array_set_offset_range(a2_cu, 0.5, a1_cu, 1*a1->ncomp, &range);

  gkyl_array_copy(a2, a2_cu);
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);
    double *a2d = gkyl_array_fetch(a2, loc);
    for (int i=0; i<a1->ncomp; ++i) {
      TEST_CHECK( a2d[i+0*a1->ncomp] == 1.5 );
      TEST_CHECK( a2d[i+1*a1->ncomp] == 0.5*0.5 );
      TEST_CHECK( a2d[i+2*a1->ncomp] == 1.5 );
    }
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

  gkyl_array_scale(a1_cu, 0.25);

 // copy from device and check if things are ok
  gkyl_array_copy(a1, a1_cu);
  for (unsigned i=0; i<a1->size; ++i)
    TEST_CHECK( gkyl_compare(a1_d[i], i*0.25, 1e-14) );

  gkyl_array_release(a1);
  gkyl_array_release(a1_cu);
}

void test_cu_array_scale_by_cell()
{
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 3, 10);
  // make device copies
  struct gkyl_array *a1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3, 10);

  struct gkyl_array *s = gkyl_array_new(GKYL_DOUBLE, 1, 10);
  // make device copies
  struct gkyl_array *s_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, 10);

  // initialize data
  double *a1_d  = a1->data;
  for (unsigned i=0; i<a1->size; ++i) {
    a1_d[i] = i*1.0;
  }
  double *s_d  = s->data;
  for (unsigned i=0; i<s->size; ++i) {
    s_d[i] = i*1.0;
  }

  // copy host arrays to device
  gkyl_array_copy(a1_cu, a1);
  gkyl_array_copy(s_cu, s);

  gkyl_array_scale_by_cell(a1_cu, s_cu);

 // copy from device and check if things are ok
  gkyl_array_copy(a1, a1_cu);
  for (unsigned i=0; i<a1->size; ++i) {
    int fact = (i/a1->ncomp);    
    TEST_CHECK( gkyl_compare(a1_d[i], i*1.0*fact, 1e-14) );
  }

  gkyl_array_release(a1);
  gkyl_array_release(a1_cu);
  gkyl_array_release(s);
  gkyl_array_release(s_cu);
}

void test_cu_array_shiftc()
{
  double s = -0.5;

  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 3, 10);
  // make device copies
  struct gkyl_array *a1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3, 10);

  // initialize data
  double *a1_d  = a1->data;
  for (unsigned i=0; i<a1->size; ++i) {
    for (size_t k=0; k<a1->ncomp; ++k) a1_d[i*a1->ncomp+k] = i*2.0+k;
  }

  // copy host arrays to device
  gkyl_array_copy(a1_cu, a1);

  gkyl_array_shiftc(a1_cu, s, 0);

  // copy from device and check if things are ok
  gkyl_array_copy(a1, a1_cu);
  TEST_CHECK( gkyl_compare(a1_d[0], 0*1.0+0-0.5, 1e-14) );
  for (unsigned i=0; i<a1->size; ++i) {
    for (size_t k=1; k<a1->ncomp; ++k)
      TEST_CHECK( gkyl_compare(a1_d[i*a1->ncomp+k], i*2.0+k, 1e-14) );
  }

  gkyl_array_release(a1);
  gkyl_array_release(a1_cu);

  // Repeat the test but shifting another coefficient as well.
  int shiftks[] = {0, 2};
  int nks = sizeof(shiftks)/sizeof(shiftks[0]);

  struct gkyl_array *a2 = gkyl_array_new(GKYL_DOUBLE, 4, 8);
  struct gkyl_array *a2_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 4, 8);
  double *a2_d = a2->data;
  for (unsigned i=0; i<a2->size; ++i) {
    for (size_t k=0; k<a2->ncomp; ++k) a2_d[i*a2->ncomp+k] = i*2.0+k;
  }
  gkyl_array_copy(a2_cu, a2);

  for (size_t l=0; l<nks; l++)
    gkyl_array_shiftc(a2_cu, s, shiftks[l]);

  gkyl_array_copy(a2, a2_cu);
  for (unsigned i=0; i<a2->size; ++i) {
    for (size_t k=0; k<a2->ncomp; ++k) {
      bool isshifted = false;
      for (size_t l=0; l<nks; l++) {
        if (shiftks[l]==k) {isshifted = true; break;}
      }
      if (isshifted)
        TEST_CHECK( gkyl_compare(a2_d[i*a2->ncomp+k], i*2.0+k+s, 1e-14) );
      else
        TEST_CHECK( gkyl_compare(a2_d[i*a2->ncomp+k], i*2.0+k, 1e-14) );
    }
  }

  gkyl_array_release(a2);
  gkyl_array_release(a2_cu);
}

void test_cu_array_copy_buffer()
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
  gkyl_array_copy_to_buffer(buff_cu, arr_cu, &sub_range);

  gkyl_array_clear(arr, 0.0);
  // copy back from buffer
  gkyl_array_copy_from_buffer(arr_cu, buff_cu, &sub_range);

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

// declare so we can use below
void set_array_copy_fn(struct gkyl_array_copy_func *fn);

void test_cu_array_copy_buffer_fn()
{
  int shape[] = {10, 20};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, 1, range.volume);

  // make device copy of array
  struct gkyl_array *arr_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, range.volume);

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

  // create function pointer on device
  struct gkyl_array_copy_func *fn = gkyl_cu_malloc(sizeof(*fn));
  set_array_copy_fn(fn);
  
  double *buff_cu = gkyl_cu_malloc(sizeof(double)*sub_range.volume);
  gkyl_array_copy_to_buffer_fn_cu(buff_cu, arr_cu, &sub_range, fn);
  // copy back from buffer
  gkyl_array_copy_from_buffer(arr_cu, buff_cu, &sub_range);

  // copy from device and check if things are ok
  gkyl_array_copy(arr, arr_cu);
  gkyl_range_iter_init(&iter, &sub_range);
  while (gkyl_range_iter_next(&iter)) {
    double *d = gkyl_array_fetch(arr, gkyl_range_idx(&sub_range, iter.idx));
    TEST_CHECK( d[0]  == 2*(iter.idx[0] + 10.5*iter.idx[1]) );
  }

  gkyl_array_release(arr);
  gkyl_array_release(arr_cu);
  gkyl_cu_free(buff_cu);
  gkyl_cu_free(fn);
}

void test_cu_array_flip_copy_buffer_fn()
{
  int shape[] = {10, 20};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);
  
  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, 1, range.volume);

  // make device copy of array
  struct gkyl_array *arr_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, range.volume);

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

  // create function pointer on device
  struct gkyl_array_copy_func *fn = gkyl_cu_malloc(sizeof(*fn));
  set_array_copy_fn(fn);
  
  double *buff_cu = gkyl_cu_malloc(sizeof(double)*sub_range.volume);

  // test flip copy on first direction of 2D array
  gkyl_array_flip_copy_to_buffer_fn_cu(buff_cu, arr_cu, 0, &sub_range, fn);

  gkyl_array_clear(arr_cu, 0.0);
  // copy back from buffer
  gkyl_array_copy_from_buffer(arr_cu, buff_cu, &sub_range);

  gkyl_array_clear(arr, 0.0);  
  // copy from device and check if things are ok
  gkyl_array_copy(arr, arr_cu);
  gkyl_range_iter_init(&iter, &sub_range);
  while (gkyl_range_iter_next(&iter)) {
    double *d = gkyl_array_fetch(arr, gkyl_range_idx(&sub_range, iter.idx));
    TEST_CHECK( d[0]  == 2*((5+1)-iter.idx[0] + 10.5*iter.idx[1]) );
  }

  // re-initialize the array
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    double *d = gkyl_array_fetch(arr, gkyl_range_idx(&range, iter.idx));
    d[0] = iter.idx[0] + 10.5*iter.idx[1];
  }
  // copy host array to device
  gkyl_array_copy(arr_cu, arr);
  
  // test flip copy on second direction of 2D array
  gkyl_array_flip_copy_to_buffer_fn(buff_cu, arr_cu, 1, &sub_range, fn);

  gkyl_array_clear(arr_cu, 0.0);
  // copy back from buffer
  gkyl_array_copy_from_buffer(arr_cu, buff_cu, &sub_range);

  gkyl_array_clear(arr, 0.0);
  // copy from device and check if things are ok
  gkyl_array_copy(arr, arr_cu);
  gkyl_range_iter_init(&iter, &sub_range);
  while (gkyl_range_iter_next(&iter)) {
    double *d = gkyl_array_fetch(arr, gkyl_range_idx(&sub_range, iter.idx));
    //printf("%lg %lg\n", d[0], 2*(iter.idx[0] + 10.5*((10+1)-iter.idx[1])) );
    TEST_CHECK( d[0]  == 2*(iter.idx[0] + 10.5*((10+1)-iter.idx[1])) );
  }

  gkyl_array_release(arr);
  gkyl_array_release(arr_cu);
  gkyl_cu_free(buff_cu);
  gkyl_cu_free(fn);
}

void test_cu_array_copy_range()
{
  int shape[] = {10, 20};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);

  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, 10, range.volume);
  struct gkyl_array *a2 = gkyl_array_new(GKYL_DOUBLE, 10, range.volume);

  // make device copies of arrays
  struct gkyl_array *a1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 10, range.volume);
  struct gkyl_array *a2_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 10, range.volume);  
  
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    double *d = gkyl_array_fetch(a1, gkyl_range_idx(&range, iter.idx));
    d[0] = iter.idx[0] + 10.5*iter.idx[1];
  }

  // copy host array to device
  gkyl_array_copy(a1_cu, a1);
  gkyl_array_copy(a2_cu, a2);

  // copy the contents of a1 into a2 over the specified range
  gkyl_array_copy_range_cu(a2_cu, a1_cu, &range);

  /// copy back to host to check contents
  gkyl_array_copy(a2, a2_cu);
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    double *d = gkyl_array_fetch(a2, gkyl_range_idx(&range, iter.idx));
    TEST_CHECK( d[0]  == iter.idx[0] + 10.5*iter.idx[1] );
  }

  // clear array for second test
  gkyl_array_clear_cu(a2_cu, 0.0);
  // initialize left sub-range
  int lower_l[] = {range.lower[0], range.lower[1]}, upper_l[] = {range.lower[0], shape[1]/2-1};
  struct gkyl_range sub_range_l;
  gkyl_sub_range_init(&sub_range_l, &range, lower_l, upper_l);
  
  // initialize right sub-range
  int lower_r[] = {range.upper[0], shape[1]/2}, upper_r[] = {range.upper[0], range.upper[1]};
  struct gkyl_range sub_range_r;
  gkyl_sub_range_init(&sub_range_r, &range, lower_r, upper_r);

  gkyl_array_copy_range_to_range_cu(a2_cu, a1_cu, &sub_range_r, &sub_range_l);

  // copy back to host to check contents
  gkyl_array_copy(a2, a2_cu);
  gkyl_range_iter_init(&iter, &sub_range_r);
  while (gkyl_range_iter_next(&iter)) {
    int idx_l[GKYL_MAX_DIM];
    idx_l[0] = range.lower[0], idx_l[1] = iter.idx[1]-shape[1]/2;
    double *d = gkyl_array_fetch(a2, gkyl_range_idx(&sub_range_r, iter.idx));
    TEST_CHECK( d[0]  == idx_l[0] + 10.5*idx_l[1] );
    TEST_MSG("Expected: %.13e in cell (%d,%d)", iter.idx[0] + 10.5*iter.idx[1], iter.idx[0], iter.idx[1]);
    TEST_MSG("Produced: %.13e", d[0]);
  }

  gkyl_array_release(a1);
  gkyl_array_release(a2);
  gkyl_array_release(a1_cu);
  gkyl_array_release(a2_cu);
}

#endif

void test_array_shiftc_range_ho() {
  test_array_shiftc_range(false);
}

void test_array_shiftc_range_dev() {
  test_array_shiftc_range(true);
}

TEST_LIST = {
  { "array_0", test_array_0 },  
  { "array_base", test_array_base },
  { "array_fetch", test_array_fetch },
  { "array_clear", test_array_clear },
  { "array_clear_range", test_array_clear_range },
  { "array_accumulate", test_array_accumulate },
  { "array_accumulate_range", test_array_accumulate_range },
  { "array_accumulate_offset", test_array_accumulate_offset },
  { "array_accumulate_offset_range", test_array_accumulate_offset_range },
  { "array_combine", test_array_combine },
  { "array_set", test_array_set },
  { "array_set_range", test_array_set_range },
  { "array_set_offset", test_array_set_offset },
  { "array_set_offset_range", test_array_set_offset_range },
  { "array_scale", test_array_scale },
  { "array_scale_by_cell", test_array_scale_by_cell },
  { "array_shiftc", test_array_shiftc },
  { "array_shiftc_range", test_array_shiftc_range_ho },
  { "array_opcombine", test_array_opcombine },
  { "array_ops_comp", test_array_ops_comp },
  { "array_copy_buffer", test_array_copy_buffer },
  { "array_copy_buffer_fn", test_array_copy_buffer_fn },
  { "array_flip_copy_buffer_fn", test_array_flip_copy_buffer_fn },
  { "array_copy_range", test_array_copy_range},
  { "array_copy_split", test_array_copy_split },
  { "non_numeric", test_non_numeric },
  { "reduce", test_reduce },
  { "reduce_dg", test_reduce_dg },
  { "reduce_range", test_reduce_range },
  { "sum_reduce_range", test_sum_reduce_range },
  { "reduce_dg_range", test_reduce_dg_range },
  { "grid_sub_array_read_1", test_grid_sub_array_read_1 },
  { "grid_sub_array_read_2", test_grid_sub_array_read_2 },  
  { "grid_array_new_from_file_1", test_grid_array_new_from_file_1 },
  { "grid_array_read_1", test_grid_array_read_p1 },
  { "array_from_buff", test_array_from_buff },  
#ifdef GKYL_HAVE_CUDA
  { "cu_array_base", test_cu_array_base },
  { "cu_array_clear", test_cu_array_clear},
  { "cu_array_clear_range", test_cu_array_clear_range},
  { "cu_array_accumulate", test_cu_array_accumulate},
  { "cu_array_accumulate_range", test_cu_array_accumulate_range},
  { "cu_array_accumulate_offset", test_cu_array_accumulate_offset},
  { "cu_array_accumulate_offset_range", test_cu_array_accumulate_offset_range},
  { "cu_array_accumulate_range_4d", test_cu_array_accumulate_range_4d  },
  { "cu_array_combine", test_cu_array_combine},
  { "cu_array_set", test_cu_array_set },
  { "cu_array_set_range", test_cu_array_set_range },
  { "cu_array_set_offset", test_cu_array_set_offset },
  { "cu_array_set_offset_range", test_cu_array_set_offset_range },
  { "cu_array_scale", test_cu_array_scale },
  { "cu_array_scale_by_cell", test_cu_array_scale_by_cell },
  { "cu_array_shiftc", test_cu_array_shiftc },
  { "cu_array_shiftc_range", test_array_shiftc_range_dev },
  { "cu_array_copy_buffer", test_cu_array_copy_buffer },
  { "cu_array_copy_buffer_fn", test_cu_array_copy_buffer_fn },
  { "cu_array_flip_copy_buffer_fn", test_cu_array_flip_copy_buffer_fn },
  { "cu_array_copy_range", test_cu_array_copy_range },
  { "cu_array_dev_kernel", test_cu_array_dev_kernel },
#endif
  { NULL, NULL },
};
