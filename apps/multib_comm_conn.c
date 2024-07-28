#include <gkyl_alloc.h>
#include <gkyl_multib_comm_conn.h>
#include <gkyl_util.h>

#include <string.h>

struct multib_comm_conn {
  struct gkyl_multib_comm_conn mcc;
};

static void
multib_comm_conn_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_multib_comm_conn *mcc = container_of(ref, struct gkyl_multib_comm_conn, ref_count);
  struct multib_comm_conn *cconn = container_of(mcc, struct multib_comm_conn, mcc);

  if (cconn->mcc.num_comm_conn > 0)
    gkyl_free(cconn->mcc.comm_conn);
  gkyl_free(cconn);
}

struct gkyl_multib_comm_conn *
gkyl_multib_comm_conn_new(int num, const struct gkyl_comm_conn *comm_conn)
{
  struct multib_comm_conn *cconn = gkyl_malloc(sizeof *cconn);
  cconn->mcc.num_comm_conn = num;

  cconn->mcc.comm_conn = 0;
  if (num > 0)
    cconn->mcc.comm_conn = gkyl_malloc(sizeof(struct gkyl_comm_conn[num]));

  for (int i=0; i<num; ++i)
    memcpy(&cconn->mcc.comm_conn[i], &comm_conn[i], sizeof(struct gkyl_comm_conn));

  cconn->mcc.ref_count = gkyl_ref_count_init(multib_comm_conn_free);

  return &cconn->mcc;
}

struct gkyl_multib_comm_conn *
gkyl_multib_comm_conn_new_send(
  int block_id, int block_rank, const int *nghost,
  const struct gkyl_block_connections *block_conn, struct gkyl_rect_decomp **decomp)
{
  int ndim = decomp[0]->ndim;
  // determine maximum number of ranks we can send data to
  int max_send_ranks = 0;
  for (int d=0; d<ndim; ++d) {
    for (int e=0; e<2; ++e) {
    
      if (block_conn->connections[d][e].edge != GKYL_PHYSICAL) {
        int nranks = decomp[block_conn->connections[d][e].bid]->ndecomp;
        max_send_ranks += nranks;
      }
    }
  }

  int comm_conn_idx = 0;
  struct gkyl_comm_conn *comm_conn
    = gkyl_malloc(sizeof(struct gkyl_comm_conn[max_send_ranks]));

  const struct gkyl_range *src_pr = &decomp[block_id]->parent_range;
  const struct gkyl_range *src_br = &decomp[block_id]->ranges[block_rank];

  for (int dir=0; dir<ndim; ++dir) {
    for (int e=0; e<2; ++e) {
      int tar_bid = block_conn->connections[dir][e].bid;
      int tar_dir = block_conn->connections[dir][e].dir;

      if (block_conn->connections[dir][e].edge != GKYL_PHYSICAL) {

        const struct gkyl_rect_decomp *tar_decomp = decomp[tar_bid];
        // to find intersection with the neighbor we first need to get
        // the neigbor block into the index-space to source block: we
        // do this by reseting the lower indices of the neigbor's
        // parent
        const struct gkyl_range *tar_pr = &tar_decomp->parent_range;

        int new_lower[GKYL_MAX_DIM];
        for (int d=0; d<ndim; ++d)
          new_lower[d] = src_pr->lower[d];

        if (GKYL_UPPER_EDGE == e)
          new_lower[dir] += gkyl_range_shape(src_pr, dir);
        else
          new_lower[dir] -= gkyl_range_shape(tar_pr, dir);
        
        struct gkyl_range reset_tar_pr;
        gkyl_range_reset_lower(&reset_tar_pr, tar_pr, new_lower);

        int delta[GKYL_MAX_DIM] = { 0 };
        for (int d=0; d<ndim; ++d)
          delta[d] = reset_tar_pr.lower[d] - tar_pr->lower[d];

        struct gkyl_range_dir_edge dir_edge = gkyl_range_edge_match(&reset_tar_pr, src_br);
        if (dir_edge.eloc != GKYL_NO_EDGE) {
          // source block range touches the parent of the neigbor
          // block: find the appropriate sub-blocks it intersects

          for (int nn=0; nn<tar_decomp->ndecomp; ++nn) {
            struct gkyl_range srange;
            gkyl_range_shift(&srange, &tar_decomp->ranges[nn], delta);

            int elo[GKYL_MAX_DIM] = { 0 }, eup[GKYL_MAX_DIM] = { 0 };
            elo[dir] = eup[dir] = 1; // only extend in 1 direction

            struct gkyl_range srange_ext;
            gkyl_range_extend(&srange_ext, &srange, elo, eup);

            struct gkyl_range irng;
            int is_inter = gkyl_range_intersect(&irng, &srange_ext, src_br);
            if (is_inter) {
              comm_conn[comm_conn_idx].rank = nn;
              comm_conn[comm_conn_idx].block_id = tar_bid;
              memcpy(&comm_conn[comm_conn_idx].range, &irng, sizeof(struct gkyl_range));

              comm_conn_idx += 1;
            }
          }

        }
      }
    }
  }

  struct gkyl_multib_comm_conn *mbcc =
    gkyl_multib_comm_conn_new(comm_conn_idx, comm_conn);

  gkyl_free(comm_conn);
  
  return mbcc;
}

void
gkyl_multib_comm_conn_release(const struct gkyl_multib_comm_conn *cconn)
{
  if (cconn)
    gkyl_ref_count_dec(&cconn->ref_count);
}
