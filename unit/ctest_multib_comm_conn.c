#include <acutest.h>

#include <gkyl_block_geom.h>
#include <gkyl_multib_comm_conn.h>

static struct gkyl_block_geom *
create_L_domain(void)
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
      .cuts = { 1, 1 },
      
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
      .cuts = { 1, 1 },
      
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
      .cuts = { 1, 1 },
      
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

static void
test_0(void)
{
  struct gkyl_comm_conn cclist[] = {
    { .rank = 1 },
    { .rank = 2 },
  };

  struct gkyl_multib_comm_conn *mbcc = gkyl_multib_comm_conn_new(2, cclist);

  TEST_CHECK( mbcc->num_comm_conn == 2 );

  gkyl_multib_comm_conn_release(mbcc);
}

static void
test_L_domain(void)
{
  struct gkyl_block_geom *geom = create_L_domain();
  struct gkyl_block_topo *topo = gkyl_block_geom_topo(geom);

  int num_blocks = topo->num_blocks;
  int num_cuts[num_blocks];
  
  // construct decomp objects
  struct gkyl_rect_decomp **decomp =
    gkyl_malloc(sizeof(struct gkyl_rect_decomp*[num_blocks]));
  for (int i=0; i<num_blocks; ++i) {
    const struct gkyl_block_geom_info *ginfo = gkyl_block_geom_get_block(geom, i);

    num_cuts[i] = 1;
    for (int d=0; d<topo->ndim; ++d)
      num_cuts[i] *= ginfo->cuts[d];
    
    struct gkyl_range range;
    gkyl_range_init_from_shape(&range, 2, ginfo->cells);
    decomp[i] = gkyl_rect_decomp_new_from_cuts(2, ginfo->cuts, &range);
  }

  for (int bid=0; bid<num_blocks; ++bid) {

    for (int brank=0; brank<num_cuts[bid]; ++brank) {
      struct gkyl_multib_comm_conn *mbcc = gkyl_multib_comm_conn_new_send(bid, brank,
        &topo->conn[bid], decomp);

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
  { "test_0", test_0 },
  { "test_L_domain", test_L_domain },
  { NULL, NULL },
};
