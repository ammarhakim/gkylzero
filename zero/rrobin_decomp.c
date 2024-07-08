#include <gkyl_alloc.h>
#include <gkyl_rrobin_decomp.h>

// internal struct to represent the round-robin algo data
struct rrobin_decomp {
  struct gkyl_rrobin_decomp rrobin;
};

const struct gkyl_rrobin_decomp*
gkyl_rrobin_decomp_new(int total_ranks, int nblocks, const int *branks)
{
  struct rrobin_decomp *rr = gkyl_malloc(sizeof(*rr));
  rr->rrobin.total_ranks = total_ranks;  
  rr->rrobin.nblocks = nblocks;

  return &rr->rrobin;
}

void
gkyl_rrobin_decomp_release(const struct gkyl_rrobin_decomp *rr)
{
  struct rrobin_decomp *rrd = container_of(rr, struct rrobin_decomp, rrobin);
  gkyl_free(rrd);
}
