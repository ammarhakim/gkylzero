#include <gkyl_gyrokinetic_multib.h>
#include <gkyl_gyrokinetic_multib_priv.h>
#include <gkyl_multib_comm_conn.h>


// Identifiers for connection type
enum gkyl_conn_id {
  GKYL_CONN_NEIGHBOR = 0, // Adjacent blocks
  GKYL_CONN_ALL = 1, // Blocks connected along one direction
  GKYL_CONN_CORNER = 2, // Blocks connected by a corner
};

/** 
 * Given a block topology, connection type, block id, and direction, 
 * get the number of connected blocks
 *
 * @param block_topo block topology object
 * @param bidx block index
 * @param dir direction in which to find neighbors or connected blocks
 * @param conn_id type of connection : GKYL_CONN_NEIGHBOR, _ALL, or _CORNER
 * return number of connected blocks
 */
int gkyl_multib_conn_get_num_connected(struct gkyl_block_topo *block_topo, int bidx, int dir, enum gkyl_conn_id conn_id);


/** 
 * Given a block topology, connection type, block id, and direction, 
 * get the number of connected blocks and a list of their ids
 *
 * @param block_topo block topology object
 * @param bidx block index
 * @param dir direction in which to find neighbors or connected blocks
 * @param conn_id type of connection : GKYL_CONN_NEIGHBOR, _ALL, or _CORNER
 * @param block_list on output, list of connected block ids
 * return number of connected blocks
 */
int gkyl_multib_conn_get_connection(struct gkyl_block_topo *block_topo, int bidx, int dir, enum gkyl_conn_id conn_id, int *block_list);

