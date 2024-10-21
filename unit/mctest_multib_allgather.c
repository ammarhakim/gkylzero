#include <acutest.h>

#include <gkyl_block_geom.h>
#include <gkyl_multib_comm_conn.h>
#include <gkyl_rrobin_decomp.h>

#ifdef GKYL_HAVE_MPI
#include <gkyl_mpi_comm.h>
#endif

#ifdef GKYL_HAVE_CUDA
#ifdef GKYL_HAVE_NCCL
#include <gkyl_nccl_comm.h>
#endif
#endif


#include <gkyl_block_geom.h>

static struct gkyl_block_geom *
create_L_domain(const int *cuts)
{
    // 2D with 3 blocks
  struct gkyl_block_geom *bgeom = gkyl_block_geom_new(2, 3);
  
  /* Block layout

     +------+
     |0     |
     |      |
     +------+-----+
     |1     |2    |
     |      |     |
     +------+-----+
    
  */

  // block 0
  gkyl_block_geom_set_block(bgeom, 0, &(struct gkyl_block_geom_info) {
      .lower = { 0, 0 },
      .upper = { 1, 1 },
      .cells = { 300, 300 },
      .cuts = { cuts[0], cuts[1] },
      
      .connections[0] = { // x-direction connections
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL }, // physical boundary
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL }  // physical boundary
      },
      .connections[1] = { // y-direction connections
        { .bid = 1, .dir = 1, .edge = GKYL_UPPER_POSITIVE },
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL } // physical boundary
      }
    }
  );
  
  // block 1
  gkyl_block_geom_set_block(bgeom, 1, &(struct gkyl_block_geom_info) {
      .lower = { 0, 0 },
      .upper = { 1, 1 },
      .cells = { 300, 300 },
      .cuts = { cuts[0], cuts[1] },
      
      .connections[0] = { // x-direction connections
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL }, // physical boundary
        { .bid = 2, .dir = 0, .edge = GKYL_LOWER_POSITIVE }
      },
      .connections[1] = { // y-direction connections
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL }, // physical boundary
        { .bid = 0, .dir = 1, .edge = GKYL_LOWER_POSITIVE }
      }
    }
  );

  // block 2
  gkyl_block_geom_set_block(bgeom, 2, &(struct gkyl_block_geom_info) {
      .lower = { 0, 0 },
      .upper = { 1, 1 },
      .cells = { 300, 300 },
      .cuts = { cuts[0], cuts[1] },
      
      .connections[0] = { // x-direction connections
        { .bid = 1, .dir = 0, .edge = GKYL_UPPER_POSITIVE },
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL } // physical boundary
      },
      .connections[1] = { // y-direction connections
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL }, // physical boundary
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL } // physical boundary
      }
    }
  );

  return bgeom;
}

// Allocate array (filled with zeros).
static struct gkyl_array*
mkarr(bool on_gpu, long nc, long size)
{
  struct gkyl_array* a;
  if (on_gpu)
    a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
  else
    a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

static void
test_L_domain_send_connections_dir0(void)
{
  struct gkyl_block_geom *geom = create_L_domain((int[]) { 1, 1 } );
  struct gkyl_block_topo *topo = gkyl_block_geom_topo(geom);

  int num_blocks = topo->num_blocks;
  int num_cuts[num_blocks];
  int nghost[] = { 1, 1 };

  // Setup for a gather along x
  int block_list[3][2] = {{0},{1,2},{1,2}};
  int dir = 0;
  int nconnected[3] = {1,2,2};
  
  // construct decomp objects
  struct gkyl_rect_decomp **decomp =
    gkyl_malloc(sizeof(struct gkyl_rect_decomp*[num_blocks]));
  for (int i=0; i<num_blocks; ++i) {
    const struct gkyl_block_geom_info *ginfo = gkyl_block_geom_get_block(geom, i);

    num_cuts[i] = 1;
    for (int d=0; d<topo->ndim; ++d)
      num_cuts[i] *= ginfo->cuts[d];
    
    struct gkyl_range range;
    gkyl_create_global_range(2, ginfo->cells, &range);
    decomp[i] = gkyl_rect_decomp_new_from_cuts(2, ginfo->cuts, &range);
  }

  // for testing
  int num_send_neigh[] = { 0, 1, 1 };

  // for testing (these hard-coded values depend on how the algorithm
  // is implemented)  
  struct gkyl_comm_conn conn_0[] = { };

  struct gkyl_comm_conn conn_1[] = {
    { .block_id = 2, .rank = 0 },    
  };
  gkyl_range_init(&conn_1[0].range, 2, (int[]) { 1, 1 }, (int[]) { 300, 300 });
  
  struct gkyl_comm_conn conn_2[] = {
    { .block_id = 1, .rank = 0 },
  };
  gkyl_range_init(&conn_2[0].range, 2, (int[]) { 301, 1 }, (int[]) { 600, 300 });
  
  struct gkyl_comm_conn *block_conn[] = { conn_0, conn_1, conn_2 };
  
  for (int bid=0; bid<num_blocks; ++bid) {

    for (int brank=0; brank<num_cuts[bid]; ++brank) {
      struct gkyl_multib_comm_conn *mbcc = gkyl_multib_comm_conn_new_send_from_connections(bid, brank, nconnected[bid], block_list[bid], dir, decomp);

      TEST_CHECK( num_send_neigh[bid] == mbcc->num_comm_conn );
      if (mbcc->num_comm_conn > 0) {
        for (int ns=0; ns<mbcc->num_comm_conn; ++ns) {
          TEST_CHECK( block_conn[bid][ns].block_id == mbcc->comm_conn[ns].block_id);
          TEST_CHECK( block_conn[bid][ns].rank == mbcc->comm_conn[ns].rank);
          TEST_CHECK( gkyl_range_compare(&block_conn[bid][ns].range, &mbcc->comm_conn[ns].range) );
        }
      }

      gkyl_multib_comm_conn_release(mbcc);
    }
  }

  for (int i=0; i<num_blocks; ++i)
    gkyl_rect_decomp_release(decomp[i]);
  gkyl_free(decomp);
  
  gkyl_block_topo_release(topo);
  gkyl_block_geom_release(geom);
}

static void
test_L_domain_recv_connections_dir0(void)
{
  struct gkyl_block_geom *geom = create_L_domain((int[]) { 1, 1 } );
  struct gkyl_block_topo *topo = gkyl_block_geom_topo(geom);

  int num_blocks = topo->num_blocks;
  int num_cuts[num_blocks];
  int nghost[] = { 1, 1 };

  // Setup for a gather along x
  int block_list[3][2] = {{0},{1,2},{1,2}};
  int dir = 0;
  int nconnected[3] = {1,2,2};
  
  // construct decomp objects
  struct gkyl_rect_decomp **decomp =
    gkyl_malloc(sizeof(struct gkyl_rect_decomp*[num_blocks]));
  for (int i=0; i<num_blocks; ++i) {
    const struct gkyl_block_geom_info *ginfo = gkyl_block_geom_get_block(geom, i);

    num_cuts[i] = 1;
    for (int d=0; d<topo->ndim; ++d)
      num_cuts[i] *= ginfo->cuts[d];
    
    struct gkyl_range range;
    gkyl_create_global_range(2, ginfo->cells, &range);
    decomp[i] = gkyl_rect_decomp_new_from_cuts(2, ginfo->cuts, &range);
  }

  // for testing
  int num_send_neigh[] = { 0, 1, 1 };

  // for testing (these hard-coded values depend on how the algorithm
  // is implemented)  
  struct gkyl_comm_conn conn_0[] = { };

  struct gkyl_comm_conn conn_1[] = {
    { .block_id = 2, .rank = 0 },    
  };
  gkyl_range_init(&conn_1[0].range, 2, (int[]) { 301, 1 }, (int[]) { 600, 300 });
  
  struct gkyl_comm_conn conn_2[] = {
    { .block_id = 1, .rank = 0 },
  };
  gkyl_range_init(&conn_2[0].range, 2, (int[]) { 1, 1 }, (int[]) { 300, 300 });
  
  struct gkyl_comm_conn *block_conn[] = { conn_0, conn_1, conn_2 };
  
  for (int bid=0; bid<num_blocks; ++bid) {

    for (int brank=0; brank<num_cuts[bid]; ++brank) {
      struct gkyl_multib_comm_conn *mbcc = gkyl_multib_comm_conn_new_recv_from_connections(bid, brank, nconnected[bid], block_list[bid], dir, decomp);

      TEST_CHECK( num_send_neigh[bid] == mbcc->num_comm_conn );
      if (mbcc->num_comm_conn > 0) {
        for (int ns=0; ns<mbcc->num_comm_conn; ++ns) {
          TEST_CHECK( block_conn[bid][ns].block_id == mbcc->comm_conn[ns].block_id);
          TEST_CHECK( block_conn[bid][ns].rank == mbcc->comm_conn[ns].rank);
          TEST_CHECK( gkyl_range_compare(&block_conn[bid][ns].range, &mbcc->comm_conn[ns].range) );
        }
      }

      gkyl_multib_comm_conn_release(mbcc);
    }
  }

  for (int i=0; i<num_blocks; ++i)
    gkyl_rect_decomp_release(decomp[i]);
  gkyl_free(decomp);
  
  gkyl_block_topo_release(topo);
  gkyl_block_geom_release(geom);
}

TEST_LIST = {
  { "test_L_domain_send_connections_dir0", test_L_domain_send_connections_dir0},
  { "test_L_domain_recv_connections_dir0", test_L_domain_recv_connections_dir0},
  { NULL, NULL },
};
