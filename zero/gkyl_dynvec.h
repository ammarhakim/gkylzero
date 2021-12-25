#pragma once

#include <gkyl_elem_type.h>

#include <stdbool.h>
#include <stddef.h>

/** Dynamic vector to store time-dependent diagnostics */
typedef struct gkyl_dynvec_tag* gkyl_dynvec;

/**
 * Create a new new dynvec. Delete using gkyl_dynvec_release method.
 * 
 * @param type Type of data in vector
 * @param ncomp Number of components
 * @return Newly allocated vector.
 */
gkyl_dynvec gkyl_dynvec_new(enum gkyl_elem_type type, size_t ncomp);

/**
 * Append data to vector. You must ensure the data has the proper type
 * and ncomp number of elements.
 * 
 * @param vec Vector to append to
 * @param tm Time-stamp of data
 * @param data to append
 */
void gkyl_dynvec_append(gkyl_dynvec vec, double tm, const void *data);

/**
 * Get data at @a idx location. You must ensure the data has the
 * proper type and ncomp number of elements.
 * 
 * @param vec Vector
 * @param idx Index of data to fetch
 * @param data On return, data is copied in this buffer
 * @return If idx is not in range, false is returned
 */
bool gkyl_dynvec_get(const gkyl_dynvec vec, size_t idx, void *data);

/**
 * Get last appended data. You must ensure the data has the proper
 * type and ncomp number of elements.
 * 
 * @param vec Vector
 * @param data On return, data is copied in this buffer
 * @return If vector is empty, return false.
 */
bool gkyl_dynvec_getlast(const gkyl_dynvec vec, void *data);

/**
 * Get time-stamp for data at index idx.
 * 
 * @param vec Vector
 * @return Time-stamp of last data at idx
 */
double gkyl_dynvec_get_tm(const gkyl_dynvec vec, size_t idx);

/**
 * Get last appended data time-stamp.
 * 
 * @param vec Vector
 * @return Time-stamp of last appened data
 */
double gkyl_dynvec_getlast_tm(const gkyl_dynvec vec);

/**
 * Get size of dynvec.
 * 
 * @param vec Vector 
 * @return Number of elements in vector
 */
size_t gkyl_dynvec_size(const gkyl_dynvec vec);

/**
 * Clear contents of the vector
 *
 * @param vec Vector to clear
 */
void gkyl_dynvec_clear(gkyl_dynvec vec);

/**
 * Clear contents of the vector, but keep the last @a num inserted
 * elements, shifting them to the start of the vector. This is
 * typically useful when the vector has been written to file and the
 * data needs to be flushed, but the vector still is in use.
 *
 * @param vec Vector to clear
 * @param num Number of final elements to keep.
 */
void gkyl_dynvec_clear_all_but(gkyl_dynvec vec, size_t num);

/**
 * Acquire a reference to the dynvec. Delete using gkyl_dynvec_release method.
 *
 * @param vec Vector to acquire reference from
 * @return Dynamic vector
 */
gkyl_dynvec gkyl_dynvec_aquire(const gkyl_dynvec vec);

/**
 * Release dynvec.
 *
 * @param vec Vector to release
 */
void gkyl_dynvec_release(gkyl_dynvec vec);

