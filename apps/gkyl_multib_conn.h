#include <gkyl_gyrokinetic_multib.h>
#include <gkyl_gyrokinetic_multib_priv.h>
#include <gkyl_multib_comm_conn.h>


// Identifiers for connection type
enum gkyl_conn_id {
  GKYL_CONN_NEIGHBOR = 0, // Adjacent blocks
  GKYL_CONN_ALL = 1, // Blocks connected along one direction
  GKYL_CONN_CORNER = 2, // Blocks connected by a corner
};


void gkyl_multib_conn_get_connection(gkyl_gyrokinetic_multib_app *mbapp, int bidx, int dir, enum gkyl_conn_id conn_id);
