#pragma once

#include <gkyl_array.h>
#include <gkyl_range.h>

enum gkyl_tensor_index_loc {
  GKYL_TENSOR_INDEX_LOWER,
  GKYL_TENSOR_INDEX_UPPER, 
};
 
struct gkyl_tensor_field {
  size_t ndim, rank; // dimension and rank of tensor
  size_t size; // number of indices

  struct gkyl_range trange; // range object for indexing tensor
  struct gkyl_array *tdata; // tensor data

  enum gkyl_tensor_index_loc iloc[GKYL_MAX_DIM]; // covariant index 0, contravariant index 1
  // note: maximum number of indices is fixed

  struct gkyl_ref_count ref_count;  
};

/**
 * Create new tensor field. Delete using gkyl_tensor_field_release method.
 * 
 * @param rank Rank of the tensor field
 * @param ndim Number of dimensions
 * @param size Number of indices 
 * @param iloc Integer index array of size GKYL_MAX_DIM which is either
 *             0 or 1, for covaraint or contravaraint elements respectively
 * @return Pointer to newly allocated tensor field.
 */
struct gkyl_tensor_field *gkyl_tensor_field_new(size_t rank, size_t ndim, size_t size, int *iloc);

/**
 * Fetches a pointer to the tensor stored at the index 'loc'.
 * 
 * @param ten Tensor field
 * @param loc Tensor to fetch
 * @return Tensor at loc
 */
GKYL_CU_DH
static inline double *
gkyl_tensor_field_fetch(struct gkyl_tensor_field *ten, long loc)
{
  return gkyl_array_fetch(ten->tdata, loc);
}

/** Same as above, except fetches a constant pointer */
GKYL_CU_DH
static inline const double *
gkyl_tensor_field_cfetch(const struct gkyl_tensor_field *ten, long loc)
{
  return gkyl_array_cfetch(ten->tdata, loc);
}

/**
 * Fetches the range index to the tensor stored at the index idx.
 * 
 * @param ten Tensor field
 * @param idx Element of the tensor
 * @return Index to the array storing the tensor field data
 */
static inline long
gkyl_tensor_field_idx(const struct gkyl_tensor_field *ten, int idx[GKYL_MAX_DIM])
{
  return gkyl_range_idx(&ten->trange, idx);
}  

/**
 * Fetches the element to the tensor stored at the index 'loc' and index idx.
 * 
 * @param ten Tensor field
 * @param loc Tensor to fetch
 * @param idx Element to fetch
 * @return Tensor at loc and element at idx
 */
static inline double
gkyl_tensor_field_elem_fetch(const struct gkyl_tensor_field *ten, long loc, int idx[GKYL_MAX_DIM])
{
  const double *val = gkyl_tensor_field_cfetch(ten, loc);
  return val[ gkyl_range_idx(&ten->trange, idx) ];
}

/**
 * Sets the element to the tensor stored at the index 'loc' and index idx.
 * 
 * @param ten Tensor field
 * @param loc Tensor to set
 * @param idx Element to set
 */
static inline void
gkyl_tensor_field_elem_set(struct gkyl_tensor_field *ten, long loc, int idx[GKYL_MAX_DIM], double ev)
{
  double *val = gkyl_tensor_field_fetch(ten, loc);
  val[ gkyl_range_idx(&ten->trange, idx) ] = ev;
}

/**
 * Release pointer to tensor field
 *
 * @param tensor_field Tensor field to release.
 */
void gkyl_tensor_field_release(const struct gkyl_tensor_field* tensor_field);