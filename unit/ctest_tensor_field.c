#include <acutest.h>
#include <mpack.h>

#include <gkyl_alloc.h>
#include <gkyl_tensor_field.h>
#include <gkyl_elem_type_priv.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>


void test_tensor_field()
{

  // Tensor field size
  int rank = 2; 
  int ndim = 3;
  int size = 10;

  // All covaraint indices
  enum gkyl_tensor_index_loc iloc[GKYL_MAX_DIM];
  for (int i=0; i<rank; ++i)  
    iloc[i] = GKYL_TENSOR_INDEX_LOWER;

  // Build and release the tensor
  struct gkyl_tensor_field *tfld = gkyl_tensor_field_new(rank,ndim,size,iloc);
  gkyl_tensor_field_release(tfld);
}

void test_tensor_field_base()
{
  int rank = 2; 
  int ndim = 3;
  int size = 10;

  enum gkyl_tensor_index_loc iloc[GKYL_MAX_DIM];
  for (int i=0; i<rank; ++i)  
    iloc[i] = GKYL_TENSOR_INDEX_LOWER;

  struct gkyl_tensor_field *tfld = gkyl_tensor_field_new(rank,ndim,size,iloc);

  TEST_CHECK( tfld->rank == 2 );
  TEST_CHECK( tfld->ndim == 3 );
  TEST_CHECK( tfld->size == 10 );
  TEST_CHECK( tfld->ref_count.count == 1 );

  TEST_CHECK( tfld->tdata->on_dev == tfld->tdata );

  TEST_CHECK( gkyl_array_is_cu_dev(tfld->tdata) == false );

  TEST_CHECK( tfld->tdata->size == size );
  TEST_CHECK( tfld->tdata->ncomp == pow(ndim,rank) );

  double *tfldData  = tfld->tdata->data;
  // Iterate over the array (with is size*(ndim)^rank)
  for (unsigned i=0; i<tfld->tdata->size; ++i){
    TEST_CHECK( tfldData[i] == 0. );
    tfldData[i] = (i+0.5)*0.1;
  }

  // acquire pointer
  struct gkyl_tensor_field *crr = gkyl_tensor_field_acquire(tfld);

  TEST_CHECK( crr->ref_count.count == 2 );
  TEST_CHECK( tfld->ref_count.count == 2 );

  struct gkyl_tensor_field *drr = gkyl_tensor_field_acquire(crr);

  TEST_CHECK( drr->ref_count.count == 3 );
  TEST_CHECK( crr->ref_count.count == 3 );  
  TEST_CHECK( tfld->ref_count.count == 3 );  
  
  gkyl_tensor_field_release(crr);
  TEST_CHECK( tfld->ref_count.count == 2 );
  gkyl_tensor_field_release(drr);
  TEST_CHECK( tfld->ref_count.count == 1 );
  
  gkyl_tensor_field_release(tfld);
}

void test_tensor_field_fetch()
{
  int rank = 2; 
  int ndim = 3;
  int size = 10;

  enum gkyl_tensor_index_loc iloc[GKYL_MAX_DIM];
  for (int i=0; i<rank; ++i)  
    iloc[i] = GKYL_TENSOR_INDEX_LOWER;

  struct gkyl_tensor_field *tfld = gkyl_tensor_field_new(rank,ndim,size,iloc);

  double *tfldData  = tfld->tdata->data;
  for (unsigned i=0; i<tfld->size; ++i){
    double *tensor = gkyl_array_fetch(tfld->tdata, i);
    for (unsigned j=0; j<tfld->tdata->ncomp; ++j) {
     tensor[j] = i + j;
    }
  }

  // Tensor fetch method
  double *tfldDataLh = gkyl_tensor_field_fetch(tfld, 0);
  TEST_CHECK( tfldDataLh[0] == (0.0 + 0.0) );

  double *tfldDataUh = gkyl_tensor_field_fetch(tfld, 8);
  TEST_CHECK( tfldDataUh[3] == (8.0 + 3.0) );

  // Tensor element fetch method
  int indxLh[GKYL_MAX_DIM] = {0.0, 0.0}; 
  indxLh[0] = 0; indxLh[1] = 2;
  double tfldDataLhElem = gkyl_tensor_field_elem_fetch(tfld, 0, indxLh);
  TEST_CHECK( tfldDataLhElem == (0.0 + 2.0) );

  int indxUh[GKYL_MAX_DIM] = {0.0, 0.0}; 
  indxUh[0] = 2; indxUh[1] = 0;
  double tfldDataUhElem = gkyl_tensor_field_elem_fetch(tfld, 5, indxUh);
  TEST_CHECK( tfldDataUhElem == (5.0 + 6.0) );
 
  gkyl_tensor_field_release(tfld);
}


TEST_LIST = {
  { "test_tensor_field", test_tensor_field },  
  { "test_tensor_field_base", test_tensor_field_base },
  { "test_tensor_field_fetch", test_tensor_field_fetch },
  { NULL, NULL },
};