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

// Get number of ranks
typedef int (*get_size_t)(struct gkyl_comm *comm, int *sz);

// "Reduce" all elements of @a type in array @a data and store output in @a out
typedef int (*all_reduce_t)(struct gkyl_comm *comm, enum gkyl_elem_type type,
  enum gkyl_array_op op, int nelem, const void *inp, void *out);

// "Synchronize" @a array across the regions or blocks.
typedef int (*gkyl_array_sync_t)(struct gkyl_comm *comm, int ndim, const int *nghost,
  struct gkyl_array *array);

// Barrier
typedef int (*barrier_t)(struct gkyl_comm *comm);

// Structure holding data and function pointers to communicate various
// Gkeyll objects across multi-region or multi-block domains
struct gkyl_comm {

  get_rank_t get_rank; // get local rank function
  get_size_t get_size; // get number of ranks
  all_reduce_t all_reduce; // all reduce function
  gkyl_array_sync_t gkyl_array_sync; // sync array
  barrier_t barrier; // barrier

  struct gkyl_ref_count ref_count; // reference count
};

/**
 * Get rank of communicator.
 *
 * @param comm Communicator
 * @param rank On output, the rank
 * @return error code: 0 for success
 */
static int
gkyl_comm_get_rank(struct gkyl_comm *comm, int *rank)
{
  return comm->get_rank(comm, rank);
}

/**
 * Get number of ranks in communcator
 *
 * @param comm Communicator
 * @param rank On output, the rank
 * @return error code: 0 for success
 */
static int
gkyl_comm_get_size(struct gkyl_comm *comm, int *sz)
{
  return comm->get_size(comm, sz);
}

/**
 * All reduce values across domains.
 *
 * @param comm Communicator
 * @param type Data-type of element
 * @param op Operator to use in reduction
 * @param nelem Number of elemets in inp and out
 * @param inp Local values on domain
 * @param out Reduced values
 * @return error code: 0 for success
 */
static int
gkyl_comm_all_reduce(struct gkyl_comm *comm, enum gkyl_elem_type type,
  enum gkyl_array_op op, int nelem, const void *inp, void *out)
{
  return comm->all_reduce(comm, type, op, nelem, inp, out);
}

/**
 * Synchronize array across domain.
 *
 * @param comm Communicator
 * @param array Array to synchronize
 * @return error code: 0 for success
 */
static int
gkyl_comm_gkyl_array_sync(struct gkyl_comm *comm, int ndim, const int *nghost,
  struct gkyl_array *array)
{
  return comm->gkyl_array_sync(comm, ndim, nghost, array);
}

/**
 * Barrier across domains
 *
 * @param comm Communcator
 * @return error code: 0 for success
 */
static int
gkyl_comm_barrier(struct gkyl_comm *comm)
{
  return comm->barrier(comm);
}

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
