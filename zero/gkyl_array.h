#pragma once

#include <gkyl_ref_count.h>
#include <gkyl_util.h>

/**
 * Array object. This is an untyped, reference counted array object.
 */
struct gkyl_array {
    int rank; // rank of array (dimensions)
    size_t elemSz; // size of element stored in array
    size_t shape[GKYL_MAX_DIM+2]; // shape of array
    size_t size; // total number of elements
    struct gkyl_ref_count ref_count; // reference count
    void *data; // pointer to data
};

/**
 * Create new array. Delete using gkyl_array_release method.
 * 
 * @param rank Rank (dimension) of array to create
 * @param elemSz Size of objects (elements) stored in array.
 * @param shape shape[d] gives number of elements in direction d
 * @return Pointer to newly allocated array.
 */
struct gkyl_array* gkyl_array_new(int rank, size_t elemSz, const size_t *shape);

/**
 * Reshape array: array is reshaped in place. Note this method changes
 * (reduces or increases) the size of array in case new size is not
 * same as original size. Data is replicated if the new size is
 * bigger.
 * 
 * @param arr Array to reshape.
 * @param rank Dimension of reshaped array
 * @param shape shape[d] is number of elements in direction d.
 * @return Pointer to reshaped array (same as 'arr')
 */
struct gkyl_array* gkyl_array_reshape(struct gkyl_array* arr,
  int rank, const size_t *shape);

/**
 * Copy into array: pointer to dest array is returned. 'dest' and
 * 'src' must not point to same data. Arrays may be of different
 * shapes and sizes.
 * 
 * @param dest Destination for copy.
 * @param src Srouce to copy from.
 */
struct gkyl_array* gkyl_array_copy(struct gkyl_array* dest,
  const struct gkyl_array* src);

/**
 * Clone array: pointer to newly created array is returned.
 * 
 * @param arr Array to clone
 * @return Pointer to clone
 */
struct gkyl_array* gkyl_array_clone(const struct gkyl_array* arr);


/**
 * Fetches a pointer to the element stored at the linear index
 * 'loc'. The calling method will usually need to cast the returned
 * pointer to the appropriate type before using it.
 *
 * @param arr Array to fetch from
 * @param loc Element to fetch
 * @return Element at location 'loc'
 */
static inline void*
gkyl_array_fetch(const struct gkyl_array* arr, long loc) {
  return ((char*) arr->data) + loc*arr->elemSz;
}

/**
 * Aquire pointer to array. The pointer must be released using
 * gkyl_array_release method.
 *
 * @param arr Array to which a pointer is needed
 * @return Pointer to aquired array
 */
struct gkyl_array* gkyl_array_aquire(const struct gkyl_array* arr);

/**
 * Release pointer to array
 *
 * @param arr Array to release.
 */
void gkyl_array_release(const struct gkyl_array* arr);
