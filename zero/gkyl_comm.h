#pragma once

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_elem_type.h>
#include <gkyl_ref_count.h>

// The return value of the functions is an error code. Success is
// denoted by 0 and failure by other values.

// Forward declaration
struct gkyl_comm;

// Get local "rank"
typedef int (*get_rank_t)(struct gkyl_comm *comm, int *rank);

// "Reduce" all elements of @a type in array @a data and store output in @a out
typedef int (*all_reduce_t)(struct gkyl_comm *comm, enum gkyl_elem_type type,
  enum gkyl_array_op op, int nelem, const void *inp, void *out);

// "Synchronize" @a array across the regions or blocks.
typedef int (*gkyl_array_sync_t)(struct gkyl_comm *comm, struct gkyl_array *array);

// Structure holding data and function pointers to communicate various
// Gkeyll objects across multi-region or multi-block domains
struct gkyl_comm {

  get_rank_t get_rank; // get local rank function
  all_reduce_t all_reduce; // all reduce function
  gkyl_array_sync_t gkyl_array_sync; // sync array

  struct gkyl_ref_count ref_count; // reference count
};

/**
 * Acquire pointer to communicator
 *
 * @param comm Communicator to to get acquire
 * @return Acquired comm obj pointer
 */
struct gkyl_comm* gkyl_comm_acquire(const struct gkyl_comm *comm);

/**
 * Release communicator memory.
 *
 * @param comm Communicator to release
 */
void gkyl_comm_release(const struct gkyl_comm *comm);
