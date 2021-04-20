#include <gkyl_alloc.h>
#include <gkyl_block_topo.h>

struct gkyl_block_topo*
gkyl_block_topo_new(int nblocks)
{
  struct gkyl_block_topo *btopo = gkyl_malloc(sizeof(struct gkyl_block_topo));
  btopo->num_blocks = nblocks;
  btopo->conn = gkyl_malloc(sizeof(struct gkyl_block_connections[nblocks]));

  return btopo;
}

void
gkyl_block_topo_free(struct gkyl_block_topo* btopo)
{
  gkyl_free(btopo->conn);
  gkyl_free(btopo); btopo = NULL;
}
