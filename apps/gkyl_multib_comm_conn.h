#pragma once

#include <gkyl_block_topo.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_ref_count.h>

// Information for block send/recv
struct gkyl_comm_conn {
  int block_id; // send/recv block ID
  int rank; // send/recv range ID in block
  struct gkyl_range range; // send/recv range
};

// List of send/recv for a given rank
struct gkyl_multib_comm_conn {
  int num_comm_conn; // number of send/recv
  struct gkyl_comm_conn *comm_conn; // communication connections (size num_comm_conn)
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
 * Construct the send communication connections for a rank from its
 * local block rank and topology connection.
 *
 * @param block_id ID of block
 * @param block_rank Local rank in block
 * @param block_conn Topological connections for block
 * @param nghost Number of ghost cells in direction d is nghost[d]
 * @param decomp List of decomposition objects for each block
 * @return New communication connection object for sends
 */
struct gkyl_multib_comm_conn *gkyl_multib_comm_conn_new_send(
  int block_id, int block_rank, const int *nghost,
  const struct gkyl_block_connections *block_conn, struct gkyl_rect_decomp **decomp);

/**
 * Construct the received communication connections for a rank from
 * its local block rank and topology connection.
 *
 * @param block_id ID of block
 * @param block_rank Local rank in block
 * @param block_conn Topological connections for block
 * @param nghost Number of ghost cells in direction d is nghost[d]
 * @param decomp List of decomposition objects for each block
 * @return New communication connection object receives
 */
struct gkyl_multib_comm_conn *gkyl_multib_comm_conn_new_recv(
  int block_id, int block_rank, const int *nghost,
  const struct gkyl_block_connections *block_conn, struct gkyl_rect_decomp **decomp);

/**
 * Release comm-conn object.
 *
 * @param cconn Object to release
 */
void gkyl_multib_comm_conn_release(const struct gkyl_multib_comm_conn *cconn);
