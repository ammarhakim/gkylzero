#include <acutest.h>
#include <mpack.h>

#include <gkyl_alloc.h>
#include <gkyl_tensor_field.h>
#include <gkyl_elem_type_priv.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>


void test_tensor_field_0()
{

  // Tensor field size
  int rank = 2; 
  int ndim = 3;
  int size = 10;

  // All covaraint indices
  int iloc[GKYL_MAX_DIM] = {0.0};

  // Build and release the tensor
  struct gkyl_tensor_field *tfld = gkyl_tensor_field_new(rank,ndim,size,iloc);
  gkyl_tensor_field_release(tfld);
}


TEST_LIST = {
  { "tensor_field_0", test_tensor_field_0 },  
  { NULL, NULL },
};