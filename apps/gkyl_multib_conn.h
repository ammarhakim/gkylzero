#include <gkyl_gyrokinetic_multib.h>
#include <gkyl_gyrokinetic_multib_priv.h>
#include <gkyl_multib_comm_conn.h>


// Identifiers for connection type
enum gkyl_conn_id {
  GKYL_CONN_ZNEIGHBOR = 0, // Adjacent blocks
  GKYL_CONN_XNEIGHBOR = 1, // Adjacent blocks
  GKYL_CONN_Z = 2, // Bocks connected along Z
  GKYL_CONN_X = 3, // Blocks connected along x
  GKYL_CONN_CORNER = 4, // Blocks connected by a corner
};


void get_connection(gkyl_gyrokinetic_multib_app *mbapp, int bidx, enum gkyl_conn_id conn_id);
