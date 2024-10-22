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


struct gkyl_comm*
comm_new(bool use_mpi, bool use_gpu, FILE *iostream)
{
  // Construct communicator for use in app.
  struct gkyl_comm *comm = 0;

#ifdef GKYL_HAVE_MPI
  if (use_gpu && use_mpi) {
#ifdef GKYL_HAVE_NCCL
    comm = gkyl_nccl_comm_new( &(struct gkyl_nccl_comm_inp) {
        .mpi_comm = MPI_COMM_WORLD,
      }
    );
#else
    fprintf(iostream, " Using -g and -M together requires NCCL.\n");
    assert(0 == 1);
#endif
  }
  else if (use_mpi) {
    comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
        .mpi_comm = MPI_COMM_WORLD,
      }
    );
  }
  else {
    comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
        .use_gpu = use_gpu
      }
    );
  }
#else
  comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
      .use_gpu = use_gpu
    }
  );
#endif

  return comm;
}

static struct gkyl_block_geom *
create_L_domain_block_geom(int **cuts)
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
  int *cuts0 = cuts[0];
  gkyl_block_geom_set_block(bgeom, 0, &(struct gkyl_block_geom_info) {
      .lower = { 0, 0 },
      .upper = { 1, 1 },
      .cells = { 300, 300 },
      .cuts = { cuts0[0], cuts0[1] },
      
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
  int *cuts1 = cuts[1];
  gkyl_block_geom_set_block(bgeom, 1, &(struct gkyl_block_geom_info) {
      .lower = { 0, 0 },
      .upper = { 1, 1 },
      .cells = { 300, 300 },
      .cuts = { cuts1[0], cuts1[1] },
      
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
  int *cuts2 = cuts[2];
  gkyl_block_geom_set_block(bgeom, 2, &(struct gkyl_block_geom_info) {
      .lower = { 0, 0 },
      .upper = { 1, 1 },
      .cells = { 300, 300 },
      .cuts = { cuts2[0], cuts2[1] },
      
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

static inline int
prod_of_elements_int(int ndim, int *arr)
{
  int pr = 1;
  for (int d=0; d<ndim; ++d) pr *= arr[d];
  return pr;
}


int **
cuts_array_new(int num_blocks, int ndim, int *cuts_all)
{
  // Create an array of cuts from an array with all the cuts listed flat.
  int **cuts_arr = gkyl_malloc(num_blocks * sizeof(int *));
  for (int i=0; i<num_blocks; i++) {
    cuts_arr[i] = gkyl_malloc(ndim * sizeof(int));

    int *cuts = cuts_arr[i];
    for (int d=0; d<ndim; d++)
      cuts[d] = cuts_all[i*ndim + d];
  }
  return cuts_arr;
}

void
cuts_array_release(int num_blocks, int **cuts_arr)
{
  // Release the array of cuts arrays.
  for (int i=0; i<num_blocks; i++)
    gkyl_free(cuts_arr[i]);
  gkyl_free(cuts_arr);
}

static void
test_L_domain_send_connections_dir0_cuts1()
{
  int num_blocks = 3; // L-shaped example.
  int ndim = 2;
  int cuts_flat1[] = {
    1, 1, // Block 0.
    1, 1, // Block 1.
    1, 1, // Block 2.
  };
  int **cuts1 = cuts_array_new(num_blocks, ndim, cuts_flat1);
  struct gkyl_block_geom *geom  = create_L_domain_block_geom(cuts1);
  struct gkyl_block_topo *topo = gkyl_block_geom_topo(geom);
  cuts_array_release(num_blocks, cuts1);

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
  int num_send_neigh[] = { 1, 2, 2 };

  // for testing (these hard-coded values depend on how the algorithm
  // is implemented)  
  struct gkyl_comm_conn conn_0[] = {
    { .block_id = 0, .rank = 0 },
  };
  gkyl_range_init(&conn_0[0].range, 2, (int[]) { 1, 1 }, (int[]) { 300, 300 });

  struct gkyl_comm_conn conn_1[] = {
    { .block_id = 1, .rank = 0 },
    { .block_id = 2, .rank = 0 },    
  };
  gkyl_range_init(&conn_1[0].range, 2, (int[]) { 1, 1 }, (int[]) { 300, 300 });
  gkyl_range_init(&conn_1[1].range, 2, (int[]) { 1, 1 }, (int[]) { 300, 300 });
  
  struct gkyl_comm_conn conn_2[] = {
    { .block_id = 1, .rank = 0 },
    { .block_id = 2, .rank = 0 },
  };
  gkyl_range_init(&conn_2[0].range, 2, (int[]) { 301, 1 }, (int[]) { 600, 300 });
  gkyl_range_init(&conn_2[1].range, 2, (int[]) { 301, 1 }, (int[]) { 600, 300 });
  
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
test_L_domain_recv_connections_dir0_cuts1()
{
  int num_blocks = 3; // L-shaped example.
  int ndim = 2;
  int cuts_flat1[] = {
    1, 1, // Block 0.
    1, 1, // Block 1.
    1, 1, // Block 2.
  };
  int **cuts1 = cuts_array_new(num_blocks, ndim, cuts_flat1);
  struct gkyl_block_geom *geom  = create_L_domain_block_geom(cuts1);
  struct gkyl_block_topo *topo = gkyl_block_geom_topo(geom);
  cuts_array_release(num_blocks, cuts1);

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
  int num_recv_neigh[] = { 1, 2, 2 };

  // for testing (these hard-coded values depend on how the algorithm
  // is implemented)  
  struct gkyl_comm_conn conn_0[] = { 
    { .block_id = 0, .rank = 0 },
  };
  gkyl_range_init(&conn_0[0].range, 2, (int[]) { 1, 1 }, (int[]) { 300, 300 });

  struct gkyl_comm_conn conn_1[] = {
    { .block_id = 1, .rank = 0 },    
    { .block_id = 2, .rank = 0 },    
  };
  gkyl_range_init(&conn_1[0].range, 2, (int[]) { 1, 1 }, (int[]) { 300, 300 });
  gkyl_range_init(&conn_1[1].range, 2, (int[]) { 301, 1 }, (int[]) { 600, 300 });
  
  struct gkyl_comm_conn conn_2[] = {
    { .block_id = 1, .rank = 0 },
    { .block_id = 2, .rank = 0 },
  };
  gkyl_range_init(&conn_2[0].range, 2, (int[]) { 1, 1 }, (int[]) { 300, 300 });
  gkyl_range_init(&conn_2[1].range, 2, (int[]) { 301, 1 }, (int[]) { 600, 300 });
  
  struct gkyl_comm_conn *block_conn[] = { conn_0, conn_1, conn_2 };
  
  for (int bid=0; bid<num_blocks; ++bid) {
    for (int brank=0; brank<num_cuts[bid]; ++brank) {
      struct gkyl_multib_comm_conn *mbcc = gkyl_multib_comm_conn_new_recv_from_connections(bid, brank, nconnected[bid], block_list[bid], dir, decomp);

      TEST_CHECK( num_recv_neigh[bid] == mbcc->num_comm_conn );
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
test_L_domain_send_connections_dir0_cuts2()
{
  int num_blocks = 3; // L-shaped example.
  int ndim = 2;
  int cuts_flat1[] = {
    1, 1, // Block 0.
    1, 1, // Block 1.
    2, 1, // Block 2.
  };
  int **cuts1 = cuts_array_new(num_blocks, ndim, cuts_flat1);
  struct gkyl_block_geom *geom  = create_L_domain_block_geom(cuts1);
  struct gkyl_block_topo *topo = gkyl_block_geom_topo(geom);

  // Construct decomp objects.
  int *branks = gkyl_malloc(sizeof(int[num_blocks]));
  for (int i=0; i<num_blocks; ++i) {
    branks[i] = prod_of_elements_int(ndim, cuts1[i]);
  }
  cuts_array_release(num_blocks, cuts1);
  int num_ranks = 1;
  const struct gkyl_rrobin_decomp* round_robin_decomp = gkyl_rrobin_decomp_new(num_ranks, num_blocks, branks);

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
  int num_send_neigh[] = { 1, 3, 3 };

  // for testing (these hard-coded values depend on how the algorithm
  // is implemented)  
  struct gkyl_comm_conn conn_0[] = { 
    { .block_id = 0, .rank = 0 },
  };
  gkyl_range_init(&conn_0[0].range, 2, (int[]) { 1, 1 }, (int[]) { 300, 300 });

  struct gkyl_comm_conn conn_1[] = {
    { .block_id = 1, .rank = 0 },    // from 0th cut
    { .block_id = 2, .rank = 0 },    // from 0th cut
    { .block_id = 2, .rank = 0 },    // from 0th cut
  };
  gkyl_range_init(&conn_1[0].range, 2, (int[]) { 1, 1 }, (int[]) { 300, 300 });
  gkyl_range_init(&conn_1[1].range, 2, (int[]) { 1, 1 }, (int[]) { 300, 300 });
  gkyl_range_init(&conn_1[2].range, 2, (int[]) { 1, 1 }, (int[]) { 300, 300 });
  
  struct gkyl_comm_conn conn_2[] = {
    { .block_id = 1, .rank = 0 }, // from 0th cut
    { .block_id = 2, .rank = 0 }, // from 0th cut
    { .block_id = 2, .rank = 0 }, // from 0th cut
    { .block_id = 1, .rank = 0 }, // from 1st cut
    { .block_id = 2, .rank = 0 }, // from 1st cut
    { .block_id = 2, .rank = 0 }, // from 1st cut
  };
  gkyl_range_init(&conn_2[0].range, 2, (int[]) { 301, 1 }, (int[]) { 450, 300 });
  gkyl_range_init(&conn_2[1].range, 2, (int[]) { 301, 1 }, (int[]) { 450, 300 });
  gkyl_range_init(&conn_2[2].range, 2, (int[]) { 301, 1 }, (int[]) { 450, 300 });
  gkyl_range_init(&conn_2[3].range, 2, (int[]) { 451, 1 }, (int[]) { 600, 300 });
  gkyl_range_init(&conn_2[4].range, 2, (int[]) { 451, 1 }, (int[]) { 600, 300 });
  gkyl_range_init(&conn_2[5].range, 2, (int[]) { 451, 1 }, (int[]) { 600, 300 });
  
  struct gkyl_comm_conn *block_conn[] = { conn_0, conn_1, conn_2 };

  for (int bid=0; bid<num_blocks; ++bid) {
    int start_ns= 0;
    for (int brank=0; brank<num_cuts[bid]; ++brank) {
      struct gkyl_multib_comm_conn *mbcc = gkyl_multib_comm_conn_new_send_from_connections(bid, brank, nconnected[bid], block_list[bid], dir, decomp);
      TEST_CHECK( num_send_neigh[bid] == mbcc->num_comm_conn );
      if (mbcc->num_comm_conn > 0) {
        for (int ns=0; ns<mbcc->num_comm_conn; ++ns) {
          // need to get the actual rank that owns this cut
          int rank_list[decomp[mbcc->comm_conn[ns].block_id]->ndecomp]; // max number of ranks
          int rank_idx = mbcc->comm_conn[ns].rank;
          gkyl_rrobin_decomp_getranks(round_robin_decomp, mbcc->comm_conn[ns].block_id, rank_list);
          mbcc->comm_conn[ns].rank = rank_list[rank_idx];
          TEST_CHECK( block_conn[bid][start_ns+ns].block_id == mbcc->comm_conn[ns].block_id);
          TEST_CHECK( block_conn[bid][start_ns+ns].rank == mbcc->comm_conn[ns].rank);
          TEST_CHECK( gkyl_range_compare(&block_conn[bid][start_ns+ns].range, &mbcc->comm_conn[ns].range) );
        }
      }

      start_ns+=mbcc->num_comm_conn;
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
test_L_domain_send_connections_dir0_cuts2_par()
{
  printf("\n");
  // Create world comm.
  bool use_mpi = true;
  bool use_gpu = false;
  struct gkyl_comm* comm = comm_new(use_mpi, use_gpu, stderr);

  int my_rank, num_ranks;
  gkyl_comm_get_rank(comm, &my_rank);
  gkyl_comm_get_size(comm, &num_ranks);


  int num_blocks = 3; // L-shaped example.
  int ndim = 2;
  int cuts_flat[] = {
    1, 1, // Block 0.
    1, 1, // Block 1.
    2, 1, // Block 2.
  };
  int **cuts = cuts_array_new(num_blocks, ndim, cuts_flat);
  struct gkyl_block_geom *geom  = create_L_domain_block_geom(cuts);
  struct gkyl_block_topo *topo = gkyl_block_geom_topo(geom);

  // Construct decomp objects.
  int *branks = gkyl_malloc(sizeof(int[num_blocks]));
  for (int i=0; i<num_blocks; ++i) {
    branks[i] = prod_of_elements_int(ndim, cuts[i]);
  }
  cuts_array_release(num_blocks, cuts);
  const struct gkyl_rrobin_decomp* round_robin_decomp = gkyl_rrobin_decomp_new(num_ranks, num_blocks, branks);

  // Check rank info. This is what was used to produce the expected results
  //int num_ranks_pb[3] = {-1};
  //num_ranks_pb[0] = gkyl_rrobin_decomp_nranks(round_robin_decomp, 0);
  //num_ranks_pb[1] = gkyl_rrobin_decomp_nranks(round_robin_decomp, 1);
  //num_ranks_pb[2] = gkyl_rrobin_decomp_nranks(round_robin_decomp, 2);

  //int ranks_0[num_ranks_pb[0]];
  //int ranks_1[num_ranks_pb[1]];
  //int ranks_2[num_ranks_pb[2]];
  //gkyl_rrobin_decomp_getranks(round_robin_decomp, 0, ranks_0);
  //gkyl_rrobin_decomp_getranks(round_robin_decomp, 1, ranks_1);
  //gkyl_rrobin_decomp_getranks(round_robin_decomp, 2, ranks_2);

  //printf("num_ranks_pb = [%d, %d, %d]\n", num_ranks_pb[0], num_ranks_pb[1], num_ranks_pb[2]);
  //printf("Block 0 has ranks: ");
  //for (int i=0; i<num_ranks_pb[0]; i++) printf("%d ", ranks_0[i]);
  //printf("\n");
  //printf("Block 1 has ranks: ");
  //for (int i=0; i<num_ranks_pb[1]; i++) printf("%d ", ranks_1[i]);
  //printf("\n");
  //printf("Block 2 has ranks: ");
  //for (int i=0; i<num_ranks_pb[2]; i++) printf("%d ", ranks_2[i]);
  //printf("\n");


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
  int num_send_neigh[] = { 1, 3, 3 };

  // for testing (these hard-coded values depend on how the algorithm
  // is implemented)  
  struct gkyl_comm_conn conn_0[] = { 
    { .block_id = 0, .rank = 0 },
  };
  gkyl_range_init(&conn_0[0].range, 2, (int[]) { 1, 1 }, (int[]) { 300, 300 });

  struct gkyl_comm_conn conn_1[] = {
    { .block_id = 1, .rank = 1 },    // from 0th cut
    { .block_id = 2, .rank = 0 },    // from 0th cut to 0th cut
    { .block_id = 2, .rank = 1 },    // from 0th cut to 1st cut
  };
  gkyl_range_init(&conn_1[0].range, 2, (int[]) { 1, 1 }, (int[]) { 300, 300 });
  gkyl_range_init(&conn_1[1].range, 2, (int[]) { 1, 1 }, (int[]) { 300, 300 });
  gkyl_range_init(&conn_1[2].range, 2, (int[]) { 1, 1 }, (int[]) { 300, 300 });
  
  struct gkyl_comm_conn conn_2[] = {
    { .block_id = 1, .rank = 1 }, // from 0th cut
    { .block_id = 2, .rank = 0 }, // from 0th cut to 0th cut
    { .block_id = 2, .rank = 1 }, // from 0th cut to 1st cut
    { .block_id = 1, .rank = 1 }, // from 1st cut to 0th cut
    { .block_id = 2, .rank = 0 }, // from 1st cut to 0th cut
    { .block_id = 2, .rank = 1 }, // from 1st cut to 1st cut
  };
  gkyl_range_init(&conn_2[0].range, 2, (int[]) { 301, 1 }, (int[]) { 450, 300 });
  gkyl_range_init(&conn_2[1].range, 2, (int[]) { 301, 1 }, (int[]) { 450, 300 });
  gkyl_range_init(&conn_2[2].range, 2, (int[]) { 301, 1 }, (int[]) { 450, 300 });
  gkyl_range_init(&conn_2[3].range, 2, (int[]) { 451, 1 }, (int[]) { 600, 300 });
  gkyl_range_init(&conn_2[4].range, 2, (int[]) { 451, 1 }, (int[]) { 600, 300 });
  gkyl_range_init(&conn_2[5].range, 2, (int[]) { 451, 1 }, (int[]) { 600, 300 });
  
  struct gkyl_comm_conn *block_conn[] = { conn_0, conn_1, conn_2 };

  for (int bid=0; bid<num_blocks; ++bid) {
    int start_ns= 0;
    for (int brank=0; brank<num_cuts[bid]; ++brank) {
      struct gkyl_multib_comm_conn *mbcc = gkyl_multib_comm_conn_new_send_from_connections(bid, brank, nconnected[bid], block_list[bid], dir, decomp);
      TEST_CHECK( num_send_neigh[bid] == mbcc->num_comm_conn );
      for (int ns=0; ns<mbcc->num_comm_conn; ++ns) {
        // need to get the actual rank that owns this cut
        int rank_list[gkyl_rrobin_decomp_nranks(round_robin_decomp, mbcc->comm_conn[ns].block_id)];
        int rank_idx = mbcc->comm_conn[ns].rank;
        gkyl_rrobin_decomp_getranks(round_robin_decomp, mbcc->comm_conn[ns].block_id, rank_list);
        mbcc->comm_conn[ns].rank = rank_list[rank_idx];
        TEST_CHECK( block_conn[bid][start_ns+ns].block_id == mbcc->comm_conn[ns].block_id);
        TEST_CHECK( block_conn[bid][start_ns+ns].rank == mbcc->comm_conn[ns].rank);
        TEST_CHECK( gkyl_range_compare(&block_conn[bid][start_ns+ns].range, &mbcc->comm_conn[ns].range) );
      }

      start_ns+=mbcc->num_comm_conn;
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
  { "test_L_domain_send_connections_dir0_cuts1", test_L_domain_send_connections_dir0_cuts1},
  { "test_L_domain_recv_connections_dir0_cuts1", test_L_domain_recv_connections_dir0_cuts1},
  { "test_L_domain_send_connections_dir0_cuts2", test_L_domain_send_connections_dir0_cuts2},
  { "test_L_domain_send_connections_dir0_cuts2_par", test_L_domain_send_connections_dir0_cuts2_par},
  { NULL, NULL },
};
