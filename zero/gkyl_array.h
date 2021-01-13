#pragma once

#include <gkyl_ref_count.h>
#include <gkyl_util.h>

/**
 * Array object. This is an untyped, undimensioned, reference counted
 * array object. All additional structure is provided else where,
 * mainly by the range object.
 */
struct gkyl_array {
    size_t elemsz; // size of elements
    size_t size; // number of elements of size elemsz
    void *data; // pointer to data
    struct gkyl_ref_count ref_count; // reference count
};

/**
 * Create new array. Delete using gkyl_array_release method.
 * 
 * @param elemsz Size of objects (elements) stored in array.
 * @param size Number elements to store in array
 * @return Pointer to newly allocated array.
 */
struct gkyl_array* gkyl_array_new(size_t elemsz, size_t size);

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
gkyl_array_fetch(const struct gkyl_array* arr, long loc)
{
  return ((char*) arr->data) + loc*arr->elemsz;
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
