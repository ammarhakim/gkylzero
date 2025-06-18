#include <gkyl_tensor_field.h>
#include <gkyl_alloc.h>
  
static void
tensor_field_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_tensor_field *tfld = container_of(ref, struct gkyl_tensor_field, ref_count);

  // release the tensorfield data
  gkyl_array_release(tfld->tdata);

  gkyl_free(tfld);
}

struct gkyl_tensor_field *
gkyl_tensor_field_new(size_t ndim, size_t rank, size_t size, int *iloc)
{
  struct gkyl_tensor_field *tfld = gkyl_malloc(sizeof *tfld);

  tfld->ndim = ndim;
  tfld->size = size;

  size_t ncomp = 1;
  int shape[rank];
  for (int i=0; i<rank; ++i) {
    ncomp *= ndim;
    shape[i] = ndim;
  }
  
  tfld->tdata = gkyl_array_new(GKYL_DOUBLE, ncomp, size);
  gkyl_range_init_from_shape(&tfld->trange, ndim, shape);

  for (int i=0; i<GKYL_MAX_DIM; ++i) {
    tfld->iloc[i] = iloc[i]; // covariant index 0, contravariant index 1
  }

  tfld->ref_count = gkyl_ref_count_init(tensor_field_free);
  
  return tfld;
}

void
gkyl_tensor_field_release(const struct gkyl_tensor_field* ten)
{
  if (ten) 
    gkyl_ref_count_dec(&ten->ref_count);
}