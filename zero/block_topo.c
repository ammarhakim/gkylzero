#include <gkyl_alloc.h>
#include <gkyl_block_topo.h>

// for use in consistency checking
static const enum gkyl_oriented_edge complimentary_edges[] = {
  [0] = 0, // can't happen for fully-specified edges
  [GKYL_LOWER_POSITIVE] = GKYL_UPPER_POSITIVE,
  [GKYL_LOWER_NEGATIVE] = GKYL_UPPER_NEGATIVE,
  [GKYL_UPPER_POSITIVE] = GKYL_LOWER_POSITIVE,
  [GKYL_UPPER_NEGATIVE] = GKYL_LOWER_NEGATIVE,
  [GKYL_PHYSICAL] = GKYL_PHYSICAL, 
};

static void
block_topo_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_block_topo *btopo = container_of(ref, struct gkyl_block_topo, ref_count);
  gkyl_free(btopo->conn);
  gkyl_free(btopo);
}

struct gkyl_block_topo*
gkyl_block_topo_new(int ndim, int nblocks)
{
  struct gkyl_block_topo *btopo = gkyl_malloc(sizeof(struct gkyl_block_topo));
  btopo->ndim = ndim;
  btopo->num_blocks = nblocks;
  btopo->conn = gkyl_calloc(sizeof(struct gkyl_block_connections), nblocks);

  btopo->ref_count = gkyl_ref_count_init(block_topo_free);

  return btopo;
}

int
gkyl_block_topo_check_consistency(const struct gkyl_block_topo *btopo)
{
  for (int i=0; i<btopo->num_blocks; ++i) {
    for (int d=0; d<btopo->ndim; ++d) {
      
      const struct gkyl_target_edge *te = btopo->conn[i].connections[d];
      
      for (int e=0; e<2; ++e) { // 0: lower, 1: upper
        if (te[e].edge < 1) // unspecified edges are defaulted to 0
          return 0;

        // check consistency
        if (te[e].edge != GKYL_PHYSICAL) {
          if (te[e].bid < 0 || te[e].bid >= btopo->num_blocks) // improperly numbered block
            return 0;

          // fetch edge which should point back to ith-block
          const struct gkyl_target_edge *te_back =
            &btopo->conn[te[e].bid].connections[d][(e+1)%2];

          // check if the edge belongs to block 'i'
          if (te_back->bid != i)
            return 0;

          // check if edge orientation is complimentary
          if (te_back->edge != complimentary_edges[te[e].edge])
            return 0;
        }
        
      }
    }
  }

  return 1;
}

struct gkyl_block_topo *
gkyl_block_topo_acquire(const struct gkyl_block_topo* btopo)
{
  gkyl_ref_count_inc(&btopo->ref_count);
  return (struct gkyl_block_topo*) btopo;
}     

void
gkyl_block_topo_release(struct gkyl_block_topo* btopo)
{
  gkyl_ref_count_dec(&btopo->ref_count);
}
