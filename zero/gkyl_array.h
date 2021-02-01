#pragma once

#include <gkyl_ref_count.h>
#include <gkyl_util.h>

#include <stdint.h>

// Type of element stored in array
enum gkyl_elem_type { GKYL_INT, GKYL_FLOAT, GKYL_DOUBLE, GKYL_USER };

/**
 * Array object. This is an untype, undimensioned, reference counted
 * array object. All additional structure is provided else where,
 * mainly by the range object.
 */
struct gkyl_array {
    enum gkyl_elem_type type; // type stored in array
    size_t elemsz, ncomp; // size of elements, number of 'components'
    size_t size; // number of indices
    
    uint32_t flags;
    size_t esznc;
    void *data;
    struct gkyl_ref_count ref_count;
};

/**
 * Create new array. Delete using gkyl_array_release method.
 * 
 * @param type Type of data in array
 * @param ncomp Number of components at each index
 * @param size Number indices 
 * @return Pointer to newly allocated array.
 */
struct gkyl_array* gkyl_array_new(enum gkyl_elem_type type, size_t ncomp, size_t size);

/**
 * Copy into array: pointer to dest array is returned. 'dest' and
 * 'src' must not point to same data.
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
 * Fetches a pointer to the element stored at the index 'loc'.
 *
 * @param arr Array to fetch from
 * @param loc Element to fetch
 * @return Element at location 'loc'
 */
static inline void*
gkyl_array_fetch(struct gkyl_array* arr, long loc)
{
  return ((char*) arr->data) + loc*arr->esznc;
}

/** Same as above, except fetches a constant pointer */
static inline const void*
gkyl_array_cfetch(const struct gkyl_array* arr, long loc)
{
  return ((const char*) arr->data) + loc*arr->esznc;
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
