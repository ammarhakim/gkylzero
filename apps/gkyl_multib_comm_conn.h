#pragma once

#include <gkyl_ref_count.h>

// Information for block send/recv
struct gkyl_comm_conn {
  int rank; // send/recv rank
};

// List of send/recv for a given rank
struct gkyl_multib_comm_conn {
  int num_comm_conn; // number of send/recv
  struct gkyl_comm_conn *comm_conn; // communication connections


  struct gkyl_ref_count ref_count;  
};

/**
 * Create new communication connection object. The list @a comm_conn
 * is managed by the caller.
 *
 * @param num Number of communication connections
 * @param comm_conn List of individual communication connections
 * @return New communication connection object
 */
struct gkyl_multib_comm_conn *gkyl_multib_comm_conn_new(int num,
  const struct gkyl_comm_conn *comm_conn);

/**
 * Release comm-conn object.
 *
 * @param cconn Object to release
 */
void gkyl_multib_comm_conn_release(const struct gkyl_multib_comm_conn *cconn);
