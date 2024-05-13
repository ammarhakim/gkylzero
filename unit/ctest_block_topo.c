#include <acutest.h>
#include <gkyl_block_topo.h>

void
test_L_domain()
{
  // 2D with 3 blocks
  struct gkyl_block_topo *btopo = gkyl_block_topo_new(2, 3);

  TEST_CHECK( 3 == btopo->num_blocks );
  // topology is inconsistent at this point!
  TEST_CHECK( 0 == gkyl_block_topo_check_consistency(btopo) );

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
  btopo->conn[0] = (struct gkyl_block_connections) {
    .connections[0] = { // x-direction connections
      { .bid = 0, .dir = 0, .edge = GKYL_BLOCK_EDGE_PHYSICAL }, // physical boundary
      { .bid = 0, .dir = 0, .edge = GKYL_BLOCK_EDGE_PHYSICAL }  // physical boundary
    },
    .connections[1] = { // y-direction connections
      { .bid = 1, .dir = 1, .edge = GKYL_BLOCK_EDGE_UPPER_POSITIVE },
      { .bid = 0, .dir = 1, .edge = GKYL_BLOCK_EDGE_PHYSICAL } // physical boundary
    }
  };
  // topology is inconsistent at this point!
  TEST_CHECK( 0 == gkyl_block_topo_check_consistency(btopo) );
  
  // block 1
  btopo->conn[1] = (struct gkyl_block_connections) {
    .connections[0] = { // x-direction connections
      { .bid = 0, .dir = 0, .edge = GKYL_BLOCK_EDGE_PHYSICAL }, // physical boundary
      { .bid = 2, .dir = 0, .edge = GKYL_BLOCK_EDGE_LOWER_POSITIVE }
    },
    .connections[1] = { // y-direction connections
      { .bid = 0, .dir = 1, .edge = GKYL_BLOCK_EDGE_PHYSICAL }, // physical boundary
      { .bid = 0, .dir = 1, .edge = GKYL_BLOCK_EDGE_LOWER_POSITIVE }
    }
  };
  // topology is inconsistent at this point!
  TEST_CHECK( 0 == gkyl_block_topo_check_consistency(btopo) );
  
  // block 2
  btopo->conn[2] = (struct gkyl_block_connections) {
    .connections[0] = { // x-direction connections
      { .bid = 1, .dir = 0, .edge = GKYL_BLOCK_EDGE_UPPER_POSITIVE },
      { .bid = 0, .dir = 0, .edge = GKYL_BLOCK_EDGE_PHYSICAL } // physical boundary
    },
    .connections[1] = { // y-direction connections
      { .bid = 0, .dir = 1, .edge = GKYL_BLOCK_EDGE_PHYSICAL }, // physical boundary
      { .bid = 0, .dir = 1, .edge = GKYL_BLOCK_EDGE_PHYSICAL } // physical boundary
    }
  };

  // should be fully consistent as all blocks properly specified
  TEST_CHECK( 1 == gkyl_block_topo_check_consistency(btopo) );

  gkyl_block_topo_release(btopo);
}

void
test_mobius_domain()
{
  // 2D with 1 block
  struct gkyl_block_topo *btopo = gkyl_block_topo_new(2, 1);

  TEST_CHECK( 1 == btopo->num_blocks );
  // topology is inconsistent at this point!
  TEST_CHECK( 0 == gkyl_block_topo_check_consistency(btopo) );

  /* Block layout

     periodic
     +------+
     |0     | periodic with twist
     |      |
     +------+
    
  */

  // block 0
  btopo->conn[0] = (struct gkyl_block_connections) {
    .connections[0] = { // x-direction connections
      { .bid = 0, .dir = 0, .edge = GKYL_BLOCK_EDGE_UPPER_NEGATIVE }, // note twist
      { .bid = 0, .dir = 0, .edge = GKYL_BLOCK_EDGE_LOWER_NEGATIVE }
    },
    .connections[1] = { // y-direction connections
      { .bid = 0, .dir = 1, .edge = GKYL_BLOCK_EDGE_UPPER_POSITIVE },
      { .bid = 0, .dir = 1, .edge = GKYL_BLOCK_EDGE_LOWER_POSITIVE }
    }
  };

  TEST_CHECK( 1 == gkyl_block_topo_check_consistency(btopo) );

  gkyl_block_topo_release(btopo);
}



TEST_LIST = {
  { "mobius_domain", test_mobius_domain },
  { "L_domain", test_L_domain },
  { NULL, NULL },
};
